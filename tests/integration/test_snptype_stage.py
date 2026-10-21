"""SNP typing stage test using an in-process fake SNP typer."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CORE_SNP_FASTA
from repgenr.core.errors import UserInputError
from repgenr.core.plugins import ToolCapabilities
from repgenr.maskers.base import Masker
from repgenr.maskers.base import registry as masker_registry
from repgenr.snptypers.base import SnpResult, SnpTyper, registry
from repgenr.stages.snptype import SnptypeParams, run


class _FakeTyper(SnpTyper):
    capabilities = ToolCapabilities(name="faketyper")
    requires_reference = False

    def preflight(self) -> dict[str, str]:
        return {"faketyper": "1.0"}

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core)


@pytest.fixture
def fake_typer():
    registry._load()
    registry.register("faketyper", _FakeTyper, replace=True)
    yield
    registry._classes.pop("faketyper", None)


def test_snptype_writes_core_snp(workdir: Path, genome_files, fake_typer) -> None:
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, SnptypeParams(tool="faketyper", all_genomes=True))

    assert result.core_snp_fasta == ctx.snp_dir / CORE_SNP_FASTA
    assert result.core_snp_fasta.exists()
    # one record per input genome
    assert result.core_snp_fasta.read_text().count(">") == len(genome_files)
    assert ctx.config.stages["snptype"].tool == "faketyper"


class _FullTyper(_FakeTyper):
    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        full = out_dir / "full.fasta"
        core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        full.write_text("".join(f">{g.stem}\nACGTACGTACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core, full_alignment=full)


class _RecordingMasker(Masker):
    capabilities = ToolCapabilities(name="fakemask")
    seen: dict = {}

    def preflight(self) -> dict[str, str]:
        return {"fakemask": "1"}

    def mask(self, full_alignment, out_dir, params, logger):  # noqa: ANN001
        _RecordingMasker.seen = {"input": full_alignment, "threads": params.threads}
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "masked.fasta"
        out.write_text(">a\nAC\n")
        return out


@pytest.fixture
def fake_masker(register_tool):
    register_tool(masker_registry, "fakemask", _RecordingMasker)
    yield


def test_masker_receives_full_alignment(workdir, genome_files, register_tool, fake_masker):
    register_tool(registry, "fulltyper", _FullTyper)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, SnptypeParams(tool="fulltyper", all_genomes=True, mask="fakemask", threads=3))
    assert _RecordingMasker.seen["input"].name == "full.fasta"
    assert _RecordingMasker.seen["threads"] == 3
    assert (ctx.snp_dir / CORE_SNP_FASTA).read_text() == ">a\nAC\n"


def test_mask_refused_without_full_alignment(workdir, genome_files, fake_typer, fake_masker):
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="whole-genome alignment"):
        run(ctx, SnptypeParams(tool="faketyper", all_genomes=True, mask="fakemask"))


class _SplitKmerTyper(_FullTyper):
    """Declares up front that it never produces a whole-genome alignment."""

    produces_full_alignment = False

    def call(self, genomes, reference, out_dir, params, logger):  # noqa: ANN001
        raise AssertionError("the typer must not run when --mask cannot be honoured")


def test_mask_refused_before_the_typer_runs(workdir, genome_files, register_tool, fake_masker):
    register_tool(registry, "splitkmer", _SplitKmerTyper)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="whole-genome alignment"):
        run(ctx, SnptypeParams(tool="splitkmer", all_genomes=True, mask="fakemask"))


class _AbsentMasker(_RecordingMasker):
    def preflight(self) -> dict[str, str]:
        from repgenr.core.errors import MissingBinaryError

        raise MissingBinaryError("absentmask: not found on PATH")


class _MustNotRunTyper(_FullTyper):
    def call(self, genomes, reference, out_dir, params, logger):  # noqa: ANN001
        raise AssertionError("the typer must not run when the masker is missing")


def test_masker_preflight_runs_before_the_typer(workdir, genome_files, register_tool):
    from repgenr.core.errors import MissingBinaryError

    register_tool(masker_registry, "absentmask", _AbsentMasker)
    register_tool(registry, "eagertyper", _MustNotRunTyper)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(MissingBinaryError, match="absentmask"):
        run(ctx, SnptypeParams(tool="eagertyper", all_genomes=True, mask="absentmask"))


class _AllOutputsTyper(_FullTyper):
    """Writes every optional output: VCF, distance matrix and full alignment."""

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        base = super().call(genomes, reference, out_dir, params, logger)
        vcf = out_dir / "calls.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")
        matrix = out_dir / "dist.tsv"
        matrix.write_text("\ta\na\t0\n")
        return SnpResult(
            core_snp_fasta=base.core_snp_fasta,
            full_alignment=base.full_alignment,
            vcf=vcf,
            snp_distance_matrix=matrix,
        )


def test_rerun_with_another_typer_drops_stale_optional_outputs(
    workdir, genome_files, register_tool, fake_typer
):
    """A typer that writes no VCF, matrix or full alignment must not leave the
    previous typer's files in snp/, where they would describe another run."""
    register_tool(registry, "alltyper", _AllOutputsTyper)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, SnptypeParams(tool="alltyper", all_genomes=True))
    for name in ("variants.vcf", "snp_distance_matrix.tsv", "full_alignment.fasta"):
        assert (ctx.snp_dir / name).exists()

    result = run(ctx, SnptypeParams(tool="faketyper", all_genomes=True))
    assert result.vcf is None and result.snp_distance_matrix is None
    for name in ("variants.vcf", "snp_distance_matrix.tsv", "full_alignment.fasta"):
        assert not (ctx.snp_dir / name).exists(), name


class _PlainRefTyper(SnpTyper):
    """A reference typer that cannot read gzip, as parsnp and snippy."""

    capabilities = ToolCapabilities(name="plainreftyper")
    requires_reference = True
    seen: dict = {}

    def preflight(self) -> dict[str, str]:
        return {"plainreftyper": "1.0"}

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        type(self).seen = {"genomes": list(genomes), "reference": reference}
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core)


def test_all_genomes_reference_that_is_also_a_representative(workdir: Path, genome_files) -> None:
    """snptype --all-genomes --reference X, with X also in representatives/,
    types X once: the reference is the genome from genomes/."""
    registry._load()
    registry.register("plainreftyper", _PlainRefTyper, replace=True)
    try:
        reps = workdir / "derep" / "representatives"
        reps.mkdir(parents=True)
        ref = genome_files[0]
        (reps / ref.name).write_bytes(ref.read_bytes())
        ctx = WorkdirContext(workdir, create=True)
        run(
            ctx,
            SnptypeParams(tool="plainreftyper", all_genomes=True, reference=ref.name),
        )
    finally:
        registry._classes.pop("plainreftyper", None)
    seen = _PlainRefTyper.seen
    assert seen["reference"] == ref
    assert seen["genomes"].count(ref) == 1
    assert ctx.snp_dir.joinpath(CORE_SNP_FASTA).read_text().count(">") == len(genome_files)


class _TreeMasker(_RecordingMasker):
    capabilities = ToolCapabilities(name="treemask")

    def mask(self, full_alignment, out_dir, params, logger):  # noqa: ANN001
        out = super().mask(full_alignment, out_dir, params, logger)
        self._tree = out_dir / "its.tree"
        self._tree.write_text("(a,b);\n", encoding="utf-8")
        return out

    def tree_output(self):
        return self._tree


def test_the_maskers_tree_is_kept_beside_the_alignment(
    workdir, genome_files, register_tool
) -> None:
    from repgenr.core.contracts import MASKER_TREE_NWK

    register_tool(registry, "fulltyper", _FullTyper)
    register_tool(masker_registry, "treemask", _TreeMasker)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, SnptypeParams(tool="fulltyper", all_genomes=True, mask="treemask"))
    assert (ctx.snp_dir / MASKER_TREE_NWK).read_text() == "(a,b);\n"
    run(ctx, SnptypeParams(tool="fulltyper", all_genomes=True))
    assert not (ctx.snp_dir / MASKER_TREE_NWK).exists()


def test_a_masker_precheck_runs_before_the_typer(workdir, genome_files, register_tool) -> None:
    from repgenr.core.errors import UserInputError

    class _PrecheckMasker(_RecordingMasker):
        capabilities = ToolCapabilities(name="precheckmask")

        def precheck(self, names, extra, logger):  # noqa: ANN001
            raise UserInputError(f"refused {len(names)} genomes before typing")

    calls: list = []

    class _CountingTyper(_FullTyper):
        def call(self, genomes, reference, out_dir, params, logger):  # noqa: ANN001
            calls.append(len(genomes))
            return super().call(genomes, reference, out_dir, params, logger)

    register_tool(registry, "countingtyper", _CountingTyper)
    register_tool(masker_registry, "precheckmask", _PrecheckMasker)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="before typing"):
        run(ctx, SnptypeParams(tool="countingtyper", all_genomes=True, mask="precheckmask"))
    assert calls == []


def test_full_alignment_typer_writes_constant_site_counts(
    workdir, genome_files, register_tool, fake_typer
) -> None:
    import json

    from repgenr.tree.alignment import ALIGNMENT_STATS

    register_tool(registry, "fulltyper", _FullTyper)
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, SnptypeParams(tool="fulltyper", all_genomes=True))
    # Every record is ACGTACGTACGT: twelve constant columns, three per base.
    assert result.constant_sites == (3, 3, 3, 3)
    stats = json.loads((ctx.snp_dir / ALIGNMENT_STATS).read_text())
    assert stats["constant_sites"] == {"A": 3, "C": 3, "G": 3, "T": 3}
    assert stats["variable_sites"] == 0 and stats["masked"] is False
    # A typer without a whole-genome alignment leaves no stale sidecar.
    run(ctx, SnptypeParams(tool="faketyper", all_genomes=True))
    assert not (ctx.snp_dir / ALIGNMENT_STATS).exists()


class _GappyTyper(_FullTyper):
    """A core where one column is called in every genome and one in a single genome."""

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        full = out_dir / "full.fasta"
        names = [g.stem for g in genomes]
        core.write_text("".join(f">{n}\n{'AC' if i == 0 else 'GN'}\n" for i, n in enumerate(names)))
        full.write_text(
            "".join(f">{n}\n{'AACC' if i == 0 else 'AAGN'}\n" for i, n in enumerate(names))
        )
        return SnpResult(core_snp_fasta=core, full_alignment=full)


def test_core_fraction_filters_the_core_and_is_recorded(
    workdir, genome_files, register_tool, caplog
) -> None:
    import json
    import logging

    from repgenr.tree.alignment import ALIGNMENT_STATS

    register_tool(registry, "gappytyper", _GappyTyper)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-snptype"))
    with caplog.at_level(logging.INFO, logger="test-snptype"):
        result = run(ctx, SnptypeParams(tool="gappytyper", all_genomes=True))
    assert any("--core-fraction 0.95 would keep" in r.getMessage() for r in caplog.records)
    assert all(
        len(ln) == 2
        for ln in result.core_snp_fasta.read_text().splitlines()
        if not ln.startswith(">")
    )
    result = run(ctx, SnptypeParams(tool="gappytyper", all_genomes=True, core_fraction=0.95))
    assert all(
        len(ln) == 1
        for ln in result.core_snp_fasta.read_text().splitlines()
        if not ln.startswith(">")
    )
    stats = json.loads((ctx.snp_dir / ALIGNMENT_STATS).read_text())
    assert stats["core_fraction"] == 0.95 and stats["constant_sites"]["A"] == 2
    assert stats["constant_sites"]["C"] == 0  # the column the filter drops
    assert ctx.config.stages["snptype"].params["core_fraction"] == 0.95
    # With the fraction set, the "would keep" report is not written.
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="test-snptype"):
        run(ctx, SnptypeParams(tool="gappytyper", all_genomes=True, core_fraction=0.95))
    assert not any("would keep" in r.getMessage() for r in caplog.records)


class _AllGappyTyper(_FakeTyper):
    """Every core column lacks a base in one genome; no full alignment (like ska2)."""

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        names = [g.stem for g in genomes]
        core.write_text("".join(f">{n}\n{'AN' if i == 0 else 'NG'}\n" for i, n in enumerate(names)))
        return SnpResult(core_snp_fasta=core)


def test_core_fraction_without_a_full_alignment_filters_and_writes_no_sidecar(
    workdir, genome_files, register_tool
) -> None:
    import logging

    from repgenr.tree.alignment import ALIGNMENT_STATS

    register_tool(registry, "allgappy", _AllGappyTyper)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-snptype"))
    # Three genomes: column 0 is called in one, column 1 in two; 0.5 keeps column 1.
    result = run(ctx, SnptypeParams(tool="allgappy", all_genomes=True, core_fraction=0.5))
    assert all(
        len(ln) == 1 for ln in result.core_snp_fasta.read_text().splitlines() if ln[0] != ">"
    )
    assert not (ctx.snp_dir / ALIGNMENT_STATS).exists()
    with pytest.raises(UserInputError, match="keeps none of the"):
        run(ctx, SnptypeParams(tool="allgappy", all_genomes=True, core_fraction=1.0))
    # The refusal leaves the typer's unfiltered core, not an empty one.
    assert all(
        len(ln) == 2 for ln in result.core_snp_fasta.read_text().splitlines() if ln[0] != ">"
    )


class _AllSitesTyper(_FakeTyper):
    """No whole-genome alignment, but an all-sites one without coordinates (like ska2).

    Column by column over the three genomes: A A A (constant A, in all),
    C C C (constant C), G T N (variable), A N N (constant A in one genome only),
    T T T (constant T).
    """

    seen: dict = {}

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        names = [g.stem for g in genomes]
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{n}\n{'GTN'[i % 3]}\n" for i, n in enumerate(names)))
        rows = ["ACGAT", "ACTNT", "ACNNT"]
        all_sites = out_dir / "all_sites.fasta"
        all_sites.write_text("".join(f">{n}\n{rows[i % 3]}\n" for i, n in enumerate(names)))
        _AllSitesTyper.seen = {"all_sites": all_sites}
        return SnpResult(core_snp_fasta=core, all_sites_alignment=all_sites)


def test_all_sites_alignment_gives_constant_site_counts(
    workdir, genome_files, register_tool, caplog
) -> None:
    import json
    import logging

    from repgenr.tree.alignment import ALIGNMENT_STATS

    assert len(genome_files) == 3
    register_tool(registry, "allsites", _AllSitesTyper)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-snptype"))
    with caplog.at_level(logging.INFO, logger="test-snptype"):
        result = run(ctx, SnptypeParams(tool="allsites", all_genomes=True))
    assert result.constant_sites == (2, 1, 0, 1)
    assert result.full_alignment is None and result.all_sites_alignment is None
    stats = json.loads((ctx.snp_dir / ALIGNMENT_STATS).read_text())
    assert stats["source"] == "all_sites.fasta"
    assert stats["constant_sites"] == {"A": 2, "C": 1, "G": 0, "T": 1}
    assert stats["variable_sites"] == 1 and stats["masked"] is False
    assert stats["alignment"] == CORE_SNP_FASTA
    # Counted, then removed: it is not published, and not left in scratch.
    assert not _AllSitesTyper.seen["all_sites"].exists()
    assert not (ctx.snp_dir / "all_sites.fasta").exists()
    assert not (ctx.snp_dir / "full_alignment.fasta").exists()
    counted = [r.getMessage() for r in caplog.records if "all_sites.fasta" in r.getMessage()]
    # One line with the counts, saying the file is removed.
    assert counted == [
        "Constant sites in all_sites.fasta: A 2, C 1, G 0, T 1 (1 variable columns); "
        "the file is not published and is removed"
    ]
    # --core-fraction applies to the counted columns as it does to the core:
    # at 0.6 the column called in one genome of three is not counted.
    run(ctx, SnptypeParams(tool="allsites", all_genomes=True, core_fraction=0.6))
    stats = json.loads((ctx.snp_dir / ALIGNMENT_STATS).read_text())
    assert stats["constant_sites"] == {"A": 1, "C": 1, "G": 0, "T": 1}
    assert stats["core_fraction"] == 0.6 and stats["variable_sites"] == 1


def test_all_sites_alignment_is_removed_after_a_refusal(
    workdir, genome_files, register_tool
) -> None:
    register_tool(registry, "allsites", _AllSitesTyper)
    ctx = WorkdirContext(workdir, create=True)
    # The core column is called in two genomes of three: 1.0 keeps none.
    with pytest.raises(UserInputError, match="keeps none of the"):
        run(ctx, SnptypeParams(tool="allsites", all_genomes=True, core_fraction=1.0))
    assert not _AllSitesTyper.seen["all_sites"].exists()


def test_neither_alignment_writes_no_stats(workdir, genome_files, register_tool) -> None:
    from repgenr.tree.alignment import ALIGNMENT_STATS

    register_tool(registry, "allsites", _AllSitesTyper)
    register_tool(registry, "faketyper", _FakeTyper)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, SnptypeParams(tool="allsites", all_genomes=True))
    assert (ctx.snp_dir / ALIGNMENT_STATS).exists()
    # A typer with neither alignment leaves no sidecar from the earlier run.
    result = run(ctx, SnptypeParams(tool="faketyper", all_genomes=True))
    assert result.constant_sites is None
    assert not (ctx.snp_dir / ALIGNMENT_STATS).exists()


def test_ska2_warns_on_a_multi_species_set(workdir, register_tool, monkeypatch, caplog) -> None:
    """The stage warns before typing with ska2 when the names span species."""
    import logging

    from repgenr.snptypers import ska2 as ska_mod

    def fake(caps, argv, **kw):  # noqa: ANN001
        argv = [str(a) for a in argv]
        if argv[1] == "build":
            Path(argv[argv.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif argv[1] == "align":
            Path(argv[argv.index("-o") + 1]).write_text(
                ">Fam_Gen_spA_GCA_1\nA\n>Fam_Gen_spB_GCA_2\nC\n", encoding="utf-8"
            )
        return 0

    monkeypatch.setattr(ska_mod, "run_tool", fake)
    monkeypatch.setattr(ska_mod.Ska2Typer, "preflight", lambda self: {"ska": "0.5.1"})
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for n in ("Fam_Gen_spA_GCA_1", "Fam_Gen_spB_GCA_2"):
        (gdir / f"{n}.fasta").write_text(f">{n}\nACGT\n", encoding="utf-8")
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-snptype"))
    with caplog.at_level(logging.WARNING, logger="test-snptype"):
        run(ctx, SnptypeParams(tool="ska2", all_genomes=True))
    assert any("spanning 2 species" in r.getMessage() for r in caplog.records)
