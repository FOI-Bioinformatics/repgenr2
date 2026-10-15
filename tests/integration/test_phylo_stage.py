"""Phylo stage composition tests using in-process fakes (no external tools)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.aligners.base import Aligner, AlignResult
from repgenr.aligners.base import registry as aligner_registry
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import record_name
from repgenr.core.plugins import ToolCapabilities
from repgenr.snptypers.base import SnpResult, SnpTyper
from repgenr.snptypers.base import registry as snp_registry
from repgenr.stages.phylo import PhyloParams, run
from repgenr.treebuilders.base import InputKind, TreeBuilder
from repgenr.treebuilders.base import registry as tb_registry


class _GenomesTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_genomes")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        leaves = [record_name(g) for g in msa_or_genomes]
        tree.write_text("(" + ",".join(leaves) + ");\n")
        return tree


def _tree_from_msa(msa: Path) -> str:
    """A star tree over the MSA's sequence names, as a real builder would label it."""
    lines = msa.read_text(encoding="utf-8").splitlines()
    names = [line[1:].split()[0] for line in lines if line.startswith(">")]
    return "(" + ",".join(names) + ");\n"


class _MsaTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_msa")
    input_kind = InputKind.MSA_FASTA
    seen_extra: dict | None = None

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).seen_extra = dict(params.extra)
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        tree.write_text(_tree_from_msa(Path(msa_or_genomes)))
        return tree


class _FakeAligner(Aligner):
    capabilities = ToolCapabilities(name="fakealigner", accepted_extras=frozenset({"kmer"}))
    seen_extra: dict | None = None

    def preflight(self):
        return {"fakealigner": "1.0"}

    def align(self, genomes, reference, out_dir, params, logger) -> AlignResult:
        type(self).seen_extra = dict(params.extra)
        out_dir.mkdir(parents=True, exist_ok=True)
        msa = out_dir / "msa.fasta"
        msa.write_text("".join(f">{record_name(g)}\nACGT\n" for g in genomes))
        return AlignResult(msa_fasta=msa)


class _FakeSnpTyper(SnpTyper):
    capabilities = ToolCapabilities(name="fakesnptyper", accepted_extras=frozenset({"kmer"}))
    requires_reference = False
    seen_extra: dict | None = None

    def preflight(self):
        return {"fakesnptyper": "1.0"}

    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:
        type(self).seen_extra = dict(params.extra)
        out_dir.mkdir(parents=True, exist_ok=True)
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{record_name(g)}\nACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core)


@pytest.fixture
def fake_phylo_tools():
    tb_registry._load()
    aligner_registry._load()
    tb_registry.register("faketree_genomes", _GenomesTreeBuilder, replace=True)
    tb_registry.register("faketree_msa", _MsaTreeBuilder, replace=True)
    aligner_registry.register("fakealigner", _FakeAligner, replace=True)
    yield
    for n in ("faketree_genomes", "faketree_msa"):
        tb_registry._classes.pop(n, None)
    aligner_registry._classes.pop("fakealigner", None)


@pytest.fixture
def fake_snptyper():
    snp_registry._load()
    snp_registry.register("fakesnptyper", _FakeSnpTyper, replace=True)
    yield
    snp_registry._classes.pop("fakesnptyper", None)


def _make_reps(workdir: Path) -> None:
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        (reps / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGTACGT\n")


def test_alignment_free_path(workdir: Path, fake_phylo_tools) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    tree = run(ctx, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True))
    assert tree.exists()
    assert tree.read_text().startswith("(")
    assert ctx.config.stages["phylo"].tool == "faketree_genomes"


def test_aligner_msa_path(workdir: Path, fake_phylo_tools) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    tree = run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
        ),
    )
    leaves = ",".join(f"Fam_gen_sp_GCA_00000{i}" for i in range(1, 4))
    assert tree.read_text().strip() == f"({leaves});"
    assert (ctx.align_dir / "msa.fasta").exists()


def test_aligner_receives_extra(workdir: Path, fake_phylo_tools) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    _FakeAligner.seen_extra = None
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
            extra={"kmer": "15"},
        ),
    )
    assert _FakeAligner.seen_extra == {"kmer": "15"}


def test_mask_key_stripped_before_aligner_and_treebuilder(workdir: Path, fake_phylo_tools) -> None:
    """'mask' is a phylo-stage-owned key: the aligner and tree builder must
    never see it, even though it rides along in params.extra."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    _FakeAligner.seen_extra = None
    _MsaTreeBuilder.seen_extra = None
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
            extra={"seed_weight": 11, "mask": "gubbins"},
        ),
    )
    assert _FakeAligner.seen_extra == {"seed_weight": 11}
    assert _MsaTreeBuilder.seen_extra == {"seed_weight": 11}


def _extras_warnings(caplog, key: str) -> list[str]:
    return [
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING and key in r.getMessage()
    ]


def test_aligner_key_does_not_warn_from_the_tree_builder(
    workdir: Path, fake_phylo_tools, caplog
) -> None:
    """An extra the aligner reads is consumed by the stage, so the tree builder
    (which declares no extras of its own) must not report it as ignored."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="aligner",
                aligner="fakealigner",
                no_outgroup=True,
                extra={"kmer": "15"},
            ),
        )
    assert _extras_warnings(caplog, "kmer") == []


def test_extra_read_by_no_tool_warns_once(workdir: Path, fake_phylo_tools, caplog) -> None:
    """A key neither the aligner nor the tree builder declares is reported once."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="aligner",
                aligner="fakealigner",
                no_outgroup=True,
                extra={"nosuchkey": "1"},
            ),
        )
    assert len(_extras_warnings(caplog, "nosuchkey")) == 1


def test_snptyper_key_does_not_warn_from_the_tree_builder(
    workdir: Path, fake_phylo_tools, fake_snptyper, caplog
) -> None:
    """Same for the snptype MSA source: the typer's key is consumed, not ignored."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="snptype",
                snptyper="fakesnptyper",
                no_outgroup=True,
                extra={"kmer": "15"},
            ),
        )
    assert _extras_warnings(caplog, "kmer") == []


def test_snptype_receives_extra_without_mask(
    workdir: Path, fake_phylo_tools, fake_snptyper
) -> None:
    """A non-mask extra key (e.g. from --aligner-arg) must reach the SNP
    typer's params.extra when msa_source=snptype; 'mask' must not (it is
    consumed by the phylo stage itself, via SnptypeParams.mask)."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    _FakeSnpTyper.seen_extra = None
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="snptype",
            snptyper="fakesnptyper",
            no_outgroup=True,
            extra={"kmer": "15", "mask": "none"},
        ),
    )
    assert _FakeSnpTyper.seen_extra == {"kmer": "15"}


def test_alignment_free_provenance_records_no_msa_source_tools(
    workdir: Path, fake_phylo_tools
) -> None:
    # A GENOMES builder never runs the aligner or SNP typer, so recording the
    # default aligner would make --aligner part of the resume fingerprint for a
    # run it does not influence.
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True))
    params = ctx.config.stages["phylo"].params
    assert params["msa_source"] is None
    assert params["aligner"] is None
    assert params["snptyper"] is None


class _SideFileTreeBuilder(TreeBuilder):
    """Writes its tree under a tool-specific name, as iqtree/raxmlng do; the stage
    must publish it as tree/tree.nwk without a window where the final path is
    truncated or half-written."""

    capabilities = ToolCapabilities(name="faketree_sidefile")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "builder_output.treefile"
        leaves = ",".join(Path(g).stem for g in msa_or_genomes)
        tree.write_text(f"({leaves});\n", encoding="utf-8")
        return tree


def test_tree_published_through_atomic_path(workdir: Path, monkeypatch, register_tool) -> None:
    from repgenr.stages import phylo as phylo_mod

    register_tool(tb_registry, "faketree_sidefile", _SideFileTreeBuilder)
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    ctx.tree_dir.mkdir(parents=True, exist_ok=True)
    (ctx.tree_dir / "tree.nwk").write_text("(old);\n", encoding="utf-8")

    published_via: list[Path] = []
    real_atomic_path = phylo_mod.atomic_path

    def spy(path: Path):
        published_via.append(path)
        return real_atomic_path(path)

    monkeypatch.setattr(phylo_mod, "atomic_path", spy)
    tree = run(ctx, PhyloParams(treebuilder="faketree_sidefile", no_outgroup=True))
    assert tree.read_text(encoding="utf-8").startswith("(Fam_gen_sp_GCA_000001,")
    assert published_via == [ctx.tree_dir / "tree.nwk"]
    assert not list(ctx.tree_dir.glob("*.part")) and not list(ctx.tree_dir.glob(".*tmp*"))


def test_snptype_source_types_the_outgroup_too(
    workdir: Path, fake_phylo_tools, fake_snptyper
) -> None:
    """The outgroup joins the SNP typing run (live audit: the ska2/simple SNP
    trees had no outgroup leaf to root on)."""
    _make_reps(workdir)
    outgroup = workdir / "outgroup"
    outgroup.mkdir()
    (outgroup / "Fam_gen_og_GCA_000009.fasta").write_text(">og\nACGTACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009\n")
    ctx = WorkdirContext(workdir, create=True)
    run(
        ctx,
        PhyloParams(treebuilder="faketree_msa", msa_source="snptype", snptyper="fakesnptyper"),
    )
    core = (workdir / "tree" / "msa" / "core_snp.fasta").read_text()
    assert core.count(">") == 4 and ">Fam_gen_og_GCA_000009" in core


def _snptype_record_and_tables(workdir: Path) -> Path:
    """A completed snptype stage: its record and its table in snp/."""
    snp = workdir / "snp"
    snp.mkdir()
    (snp / "core_snp.fasta").write_text(">from_snptype\nACGT\n")
    ctx = WorkdirContext(workdir, create=True)
    ctx.config.record_stage("snptype", tool="ska2", completed="2026-10-07T00:00:00")
    ctx.save_config()
    return snp


def test_snptype_source_types_under_tree_msa_and_leaves_snp_alone(
    workdir: Path, fake_phylo_tools, fake_snptyper
) -> None:
    """snp/ belongs to the snptype stage: phylo's typing pass writes its
    alignment and reuse stamp under tree/msa/, so the snptype tables and
    record stay as they were."""
    from repgenr.core.config import Config

    _make_reps(workdir)
    snp = _snptype_record_and_tables(workdir)
    ctx = WorkdirContext(workdir)
    params = PhyloParams(
        treebuilder="faketree_msa",
        msa_source="snptype",
        snptyper="fakesnptyper",
        no_outgroup=True,
    )
    run(ctx, params)
    assert (snp / "core_snp.fasta").read_text() == ">from_snptype\nACGT\n"
    assert sorted(p.name for p in snp.iterdir()) == ["core_snp.fasta"]
    msa_dir = workdir / "tree" / "msa"
    assert (msa_dir / "core_snp.fasta").read_text().count(">") == 3
    assert (msa_dir / "msa_source.json").is_file()
    assert not (workdir / "scratch" / "snptype").exists(), "snptype's scratch is not shared"
    stages = Config.load(workdir).stages
    assert "snptype" in stages and "phylo" in stages


def test_aligner_source_keeps_the_snptype_record(workdir: Path, fake_phylo_tools) -> None:
    from repgenr.core.config import Config

    _make_reps(workdir)
    snp = _snptype_record_and_tables(workdir)
    ctx = WorkdirContext(workdir)
    run(ctx, PhyloParams(treebuilder="faketree_msa", aligner="fakealigner", no_outgroup=True))
    assert "snptype" in Config.load(workdir).stages
    assert (snp / "core_snp.fasta").read_text() == ">from_snptype\nACGT\n"


def _type_calls(monkeypatch) -> list[int]:
    """Count SNP typer invocations across phylo runs."""
    calls: list[int] = []
    original = _FakeSnpTyper.call

    def counting(self, genomes, reference, out_dir, params, logger):
        calls.append(1)
        return original(self, genomes, reference, out_dir, params, logger)

    monkeypatch.setattr(_FakeSnpTyper, "call", counting)
    return calls


def test_snptype_msa_is_reused_from_tree_msa(
    workdir: Path, fake_phylo_tools, fake_snptyper, monkeypatch
) -> None:
    """The typed alignment under tree/msa/ survives the tree builder's cleanup
    of tree/ and is reused when only the tree builder settings change; a
    snptype stage run in between does not invalidate it."""
    from repgenr.stages.snptype import SnptypeParams
    from repgenr.stages.snptype import run as snptype_run

    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _type_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="snptype", snptyper="fakesnptyper")
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1
    snptype_run(ctx, SnptypeParams(tool="fakesnptyper"))
    assert len(calls) == 2
    # The snptype stage's table differs from phylo's alignment (another typer
    # or reference); it must not invalidate the alignment phylo built.
    (workdir / "snp" / "core_snp.fasta").write_text(">other\nTTTT\n")
    run(ctx, PhyloParams(no_outgroup=True, bootstrap=100, **base))
    assert len(calls) == 2, "same inputs and settings: the alignment is reused"


def test_stamp_left_under_snp_by_an_earlier_layout_is_not_reused(
    workdir: Path, fake_phylo_tools, fake_snptyper, monkeypatch
) -> None:
    """A workdir whose phylo typing pass wrote snp/ (before tree/msa/) is typed
    once more, into tree/msa/; the old snp/ files are not read or changed."""
    import shutil

    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _type_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="snptype", snptyper="fakesnptyper")
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1
    msa_dir = workdir / "tree" / "msa"
    shutil.move(str(msa_dir), str(workdir / "snp"))
    legacy = {p.name: p.read_bytes() for p in (workdir / "snp").iterdir()}
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 2
    assert (msa_dir / "msa_source.json").is_file()
    assert {p.name: p.read_bytes() for p in (workdir / "snp").iterdir()} == legacy


def test_alignment_free_builder_warns_that_msa_options_have_no_effect(
    workdir: Path, fake_phylo_tools, caplog
) -> None:
    """--msa-source snptype and --mask are dropped by a builder that reads genomes."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_genomes",
                msa_source="snptype",
                no_outgroup=True,
                extra={"mask": "gubbins"},
            ),
        )
    warned = [r.getMessage() for r in caplog.records if "has no effect" in r.getMessage()]
    assert len(warned) == 1
    assert "--msa-source snptype" in warned[0] and "--mask gubbins" in warned[0]
    assert not (workdir / "snp").exists()
    assert not (workdir / "tree" / "msa").exists()


def _align_calls(monkeypatch) -> list[int]:
    """Count aligner invocations across phylo runs."""
    calls: list[int] = []
    original = _FakeAligner.align

    def counting(self, genomes, reference, out_dir, params, logger):
        calls.append(1)
        return original(self, genomes, reference, out_dir, params, logger)

    monkeypatch.setattr(_FakeAligner, "align", counting)
    return calls


def test_msa_is_reused_when_only_the_tree_builder_changes(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    """Trying another tree builder must not repeat the alignment."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _align_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="aligner", aligner="fakealigner")
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1
    stamp = workdir / "align" / "msa_source.json"
    assert stamp.is_file(), "the alignment is stamped with what produced it"

    run(ctx, PhyloParams(no_outgroup=True, bootstrap=100, **base))
    assert len(calls) == 1, "same inputs and settings: the alignment is reused"
    assert ctx.config.stages["phylo"].tool_versions.get("fakealigner") == "1.0"


def test_msa_is_rebuilt_when_its_own_settings_or_inputs_change(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _align_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="aligner", aligner="fakealigner")
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1

    # An aligner setting the alignment depends on.
    run(ctx, PhyloParams(no_outgroup=True, extra={"kmer": 31}, **base))
    assert len(calls) == 2

    # A genome added to the set.
    reps = workdir / "derep" / "representatives"
    (reps / "Fam_gen_sp_GCA_0000099.fasta").write_text(">s9\nACGTACGT\n")
    run(ctx, PhyloParams(no_outgroup=True, extra={"kmer": 31}, **base))
    assert len(calls) == 3

    # The alignment file itself replaced behind our back.
    (workdir / "align" / "msa.fasta").write_text(">s1\nAAAA\n")
    run(ctx, PhyloParams(no_outgroup=True, extra={"kmer": 31}, **base))
    assert len(calls) == 4, "a stamp that no longer describes the file is not trusted"


def test_force_rebuilds_the_msa(workdir: Path, fake_phylo_tools, monkeypatch) -> None:
    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _align_calls(monkeypatch)
    params = PhyloParams(
        treebuilder="faketree_msa", msa_source="aligner", aligner="fakealigner", no_outgroup=True
    )
    run(ctx, params)
    ctx.force = True
    run(ctx, params)
    assert len(calls) == 2


def test_provenance_records_reference_mask_and_extras(workdir: Path, fake_phylo_tools) -> None:
    """repgenr.yaml must tell a masked tree from an unmasked one, and which
    reference and adapter tuning produced the alignment."""
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    run(
        ctx,
        PhyloParams(
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
            no_outgroup=True,
            reference="Fam_gen_sp_GCA_000002.fasta",
            extra={"kmer": "15"},
        ),
    )
    params = ctx.config.stages["phylo"].params
    assert params["reference"] == "Fam_gen_sp_GCA_000002.fasta"
    assert params["mask"] is None
    assert params["extra"] == {"kmer": "15"}


def test_rebuild_with_another_builder_drops_the_previous_builders_files(
    workdir: Path, fake_phylo_tools, register_tool
) -> None:
    """tree/ holds the current builder's own files: a rebuild with another builder
    must not leave the previous one's side files (matrices, bootstrap trees)."""
    register_tool(tb_registry, "faketree_sidefile", _SideFileTreeBuilder)
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, PhyloParams(treebuilder="faketree_sidefile", no_outgroup=True))
    assert (ctx.tree_dir / "builder_output.treefile").exists()
    (ctx.tree_dir / "signatures").mkdir()
    (ctx.tree_dir / "signatures" / "a.sig").write_text("x")

    tree = run(ctx, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True))
    assert tree.exists()
    assert sorted(p.name for p in ctx.tree_dir.iterdir()) == ["tree.nwk"]


def test_builder_file_cleanup_tolerates_entries_that_vanish(tmp_path: Path, monkeypatch) -> None:
    """On non-HFS volumes macOS removes a file's ``._`` sibling with it, so a
    listed entry can be gone by the time the cleanup reaches it."""
    from repgenr.stages import phylo as phylo_mod

    tree_dir = tmp_path / "tree"
    tree_dir.mkdir()
    (tree_dir / "tree.nwk").write_text("(a);\n")
    (tree_dir / "matrix.tsv").write_text("x")
    real_iterdir = Path.iterdir

    def iterdir_with_ghost(self: Path):
        yield from real_iterdir(self)
        if self == tree_dir:
            yield tree_dir / "._matrix.tsv"  # listed, already gone

    monkeypatch.setattr(Path, "iterdir", iterdir_with_ghost)
    phylo_mod._clear_previous_builder_files(tree_dir)
    assert sorted(p.name for p in real_iterdir(tree_dir)) == ["tree.nwk"]


def _make_n_reps(workdir: Path, n: int) -> None:
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, n + 1):
        (reps / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGTACGT\n")


@pytest.mark.parametrize("n", [1, 2])
def test_phylo_refuses_fewer_than_three_genomes(
    workdir: Path, fake_phylo_tools, monkeypatch, n: int
) -> None:
    from repgenr.core.errors import WorkdirError

    calls: list[int] = []
    monkeypatch.setattr(_GenomesTreeBuilder, "build", lambda *a, **k: calls.append(1))
    _make_n_reps(workdir, n)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(WorkdirError) as exc:
        run(ctx, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True))
    assert str(exc.value) == (
        f"A tree needs at least 3 genomes; {n} found after dereplication "
        "(use --all-genomes, or dereplicate with a higher --secondary-ani, which "
        "keeps more representatives)."
    )
    assert exc.value.exit_code == 3
    assert calls == []


def test_phylo_accepts_exactly_three_genomes(workdir: Path, fake_phylo_tools) -> None:
    _make_n_reps(workdir, 3)
    ctx = WorkdirContext(workdir, create=True)
    assert run(ctx, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True)).exists()


class _DroppingTreeBuilder(TreeBuilder):
    """Writes a tree without the last genome, as mashtree can for a degenerate one."""

    capabilities = ToolCapabilities(name="faketree_drop")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        leaves = [Path(g).stem for g in msa_or_genomes][:-1]
        tree.write_text("(" + ",".join(leaves) + ",extra_leaf);\n")
        return tree


class _RenamingTreeBuilder(TreeBuilder):
    """Writes leaves with '.' replaced by '_', as cactus names its samples."""

    capabilities = ToolCapabilities(name="faketree_rename")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        leaves = [Path(g).stem.replace(".", "_") for g in msa_or_genomes]
        tree.write_text("(" + ",".join(leaves) + ");\n")
        return tree


def test_tree_missing_a_genome_is_refused(workdir: Path, fake_phylo_tools) -> None:
    from repgenr.core.errors import WorkdirError

    tb_registry.register("faketree_drop", _DroppingTreeBuilder, replace=True)
    try:
        _make_reps(workdir)
        ctx = WorkdirContext(workdir, create=True)
        with pytest.raises(WorkdirError) as excinfo:
            run(ctx, PhyloParams(treebuilder="faketree_drop", no_outgroup=True))
    finally:
        tb_registry._classes.pop("faketree_drop", None)
    message = str(excinfo.value)
    assert "faketree_drop" in message
    assert "Fam_gen_sp_GCA_000003" in message
    assert "extra_leaf" in message
    assert excinfo.value.exit_code == 3
    # the tree stays on disk for inspection; the stage is not recorded
    assert (workdir / "tree" / "tree.nwk").exists()
    assert "phylo" not in ctx.config.stages


def test_tree_leaf_check_tolerates_renamed_dots(workdir: Path, fake_phylo_tools) -> None:
    tb_registry.register("faketree_rename", _RenamingTreeBuilder, replace=True)
    try:
        reps = workdir / "derep" / "representatives"
        reps.mkdir(parents=True)
        for i in range(1, 4):
            (reps / f"Fam_gen_sp_GCA_00000{i}.1.fasta").write_text(f">s{i}\nACGTACGT\n")
        ctx = WorkdirContext(workdir, create=True)
        tree = run(ctx, PhyloParams(treebuilder="faketree_rename", no_outgroup=True))
    finally:
        tb_registry._classes.pop("faketree_rename", None)
    assert tree.exists()
    assert ctx.config.stages["phylo"].completed


class _RenamingLengthTreeBuilder(TreeBuilder):
    """A rooted tree with lengths and a support, its leaves renamed as cactus
    (dots to '_') and harvesttools ('.fasta', '.ref') name them."""

    capabilities = ToolCapabilities(name="faketree_rename_len")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        a, b, c = (Path(g).stem for g in msa_or_genomes)
        a, b = a.replace(".", "_"), f"{b}.fasta.ref"
        tree.write_text(f"(({a}:0.1,{b}:0.2)0.95:0.3,{c}.fasta:0.4);\n")
        return tree


def test_renamed_leaves_are_restored_to_the_input_names(workdir: Path, fake_phylo_tools) -> None:
    """tree2tax reads leaves as genome names, so phylo writes them back as the
    inputs were named, keeping branch lengths and supports."""
    import dendropy

    from repgenr.stages.tree2tax import Tree2taxParams
    from repgenr.stages.tree2tax import run as tree2tax_run

    tb_registry.register("faketree_rename_len", _RenamingLengthTreeBuilder, replace=True)
    try:
        reps = workdir / "derep" / "representatives"
        reps.mkdir(parents=True)
        stems = [f"Fam_gen_sp_GCA_00000{i}.1" for i in range(1, 4)]
        for stem in stems:
            (reps / f"{stem}.fasta").write_text(">s\nACGTACGT\n")
        ctx = WorkdirContext(workdir, create=True)
        tree = run(ctx, PhyloParams(treebuilder="faketree_rename_len", no_outgroup=True))
    finally:
        tb_registry._classes.pop("faketree_rename_len", None)
    parsed = dendropy.Tree.get(path=str(tree), schema="newick", preserve_underscores=True)
    assert sorted(n.taxon.label for n in parsed.leaf_node_iter()) == stems
    lengths = sorted(e.length for e in parsed.postorder_edge_iter() if e.length is not None)
    assert lengths == [0.1, 0.2, 0.3, 0.4]
    assert [n.label for n in parsed.internal_nodes() if n.label] == ["0.95"]
    # The taxonomy maps every accession to a leaf of the same name.
    _, gmap = tree2tax_run(WorkdirContext(workdir), Tree2taxParams())
    rows = dict(ln.split("\t") for ln in gmap.read_text().splitlines())
    assert rows == {stem.removeprefix("Fam_gen_sp_"): stem for stem in stems}


def test_restore_leaf_names_leaves_ambiguous_names_alone(tmp_path: Path) -> None:
    """Two inputs that a tool would write as one name are not guessed between."""
    from repgenr.stages.phylo import restore_leaf_names

    tree = tmp_path / "tree.nwk"
    tree.write_text("(x_1,y,z);\n")
    assert restore_leaf_names(tree, ["x.1", "x_1", "y", "z"], logging.getLogger("t")) == 0
    assert tree.read_text() == "(x_1,y,z);\n"


def test_restore_leaf_names_changes_only_the_renamed_labels(tmp_path: Path) -> None:
    """Review of #223: the rewrite kept neither a leading [&R] nor a quoted label
    with a space that was not renamed. Only matched leaf labels change now."""
    from repgenr.stages.phylo import restore_leaf_names

    tree = tmp_path / "tree.nwk"
    original = "[&R] (('GCF 3':0.1,x_GCF_1_1.fasta.ref:0.2)'node, a':0.3,[c;m] x_GCF_2_1:0.4)0.9;\n"
    tree.write_text(original)
    expected = ["GCF 3", "x_GCF_1.1", "x_GCF_2.1"]
    assert restore_leaf_names(tree, expected, logging.getLogger("t")) == 2
    assert tree.read_text() == (
        "[&R] (('GCF 3':0.1,x_GCF_1.1:0.2)'node, a':0.3,[c;m] x_GCF_2.1:0.4)0.9;\n"
    )


def test_rename_quotes_a_new_name_that_needs_it() -> None:
    from repgenr.stages.phylo import _rename_newick_leaves

    text, n = _rename_newick_leaves("(a_b,'c''d',e);", {"a_b": "a b", "c'd": "c'e"}.get)
    assert n == 2
    assert text == "('a b','c''e',e);"


def test_leaf_key_matches_tool_rewritten_names() -> None:
    from repgenr.stages.phylo import _leaf_key

    stem = "Fam_gen_sp_GCA_000001.1"
    for label in (stem, "Fam_gen_sp_GCA_000001_1", f"{stem}.fasta", f"{stem}.fna.ref"):
        assert _leaf_key(label) == _leaf_key(stem)
    assert _leaf_key("Fam_gen_sp_GCA_000002.1") != _leaf_key(stem)


def test_msa_stamped_by_an_earlier_version_is_not_reused(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    """An alignment from before the stamp version changed (for example a snippy
    alignment that still names the reference 'Reference') is rebuilt."""
    from repgenr.stages import phylo as phylo_mod

    # 3: the parsnp and cactus record names changed (#223).
    assert phylo_mod._MSA_STAMP_VERSION >= 3
    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _align_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="aligner", aligner="fakealigner")
    with monkeypatch.context() as m:
        m.setattr(phylo_mod, "_MSA_STAMP_VERSION", 2)
        run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 2, "a stamp from an earlier version is not trusted"


def test_aligner_msa_named_by_path_stem_is_rebuilt(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    """An aligner MSA whose records carry 'x.fasta' for a gzipped genome x.fasta.gz
    (the Path.stem names before the record-name rule) is rebuilt once, without
    a stamp version change that would realign every workdir."""
    import gzip

    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        with gzip.open(reps / f"Fam_gen_sp_GCA_00000{i}.fasta.gz", "wt") as fh:
            fh.write(f">s{i}\nACGTACGT\n")
    ctx = WorkdirContext(workdir)
    calls = _align_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="aligner", aligner="fakealigner")

    def stem_named(self, genomes, reference, out_dir, params, logger) -> AlignResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        msa = out_dir / "msa.fasta"
        # The names Path.stem gave the gzipped inputs before record_name.
        msa.write_text("".join(f">{record_name(g)}.fasta\nACGT\n" for g in genomes))
        return AlignResult(msa_fasta=msa)

    with monkeypatch.context() as m:
        m.setattr(_FakeAligner, "align", stem_named)
        run(ctx, PhyloParams(no_outgroup=True, **base))
    msa = workdir / "align" / "msa.fasta"
    assert ">Fam_gen_sp_GCA_000001.fasta\n" in msa.read_text()
    assert (workdir / "align" / "msa_source.json").is_file()

    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1, "records that are not the genome record names: rebuilt"
    assert ">Fam_gen_sp_GCA_000001\n" in msa.read_text()
    run(ctx, PhyloParams(no_outgroup=True, **base))
    assert len(calls) == 1, "the rebuilt alignment is reused"


def test_fasta_record_names_reads_only_headers(tmp_path: Path) -> None:
    from repgenr.stages.phylo import fasta_record_names

    f = tmp_path / "a.fasta"
    f.write_bytes(b">x desc\r\nACGT\nAC>GT\n>y\n" + b"A" * 200_000 + b"\n>\n>z")
    assert fasta_record_names(f) == ["x", "y", "", "z"]
    f.write_bytes(b"junk\n>x\nA\n")
    assert fasta_record_names(f) == ["x"]
    f.write_bytes(b"")
    assert fasta_record_names(f) == []
    f.write_bytes(b"ACGT\n")
    assert fasta_record_names(f) == []


def test_reused_msa_is_checked_from_the_stamp_without_reading_it(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    """The stamp holds the record names, so a reuse does not scan the MSA; a
    stamp written without them is completed once."""
    import json

    from repgenr.stages import phylo as phylo_mod

    _make_reps(workdir)
    ctx = WorkdirContext(workdir)
    calls = _align_calls(monkeypatch)
    base = dict(treebuilder="faketree_msa", msa_source="aligner", aligner="fakealigner")
    run(ctx, PhyloParams(no_outgroup=True, **base))
    stamp_path = workdir / "align" / "msa_source.json"
    stamp = json.loads(stamp_path.read_text())
    assert stamp["records"] == [f"Fam_gen_sp_GCA_00000{i}" for i in range(1, 4)]

    # A stamp from before the names were recorded is completed once.
    del stamp["records"]
    stamp_path.write_text(json.dumps(stamp))
    run(ctx, PhyloParams(no_outgroup=True, bootstrap=10, **base))
    assert len(calls) == 1
    assert "records" in json.loads(stamp_path.read_text())

    def no_scan(path):
        raise AssertionError("a reuse with recorded names must not scan the MSA")

    monkeypatch.setattr(phylo_mod, "fasta_record_names", no_scan)
    run(ctx, PhyloParams(no_outgroup=True, bootstrap=20, **base))
    assert len(calls) == 1, "reused"


def test_aligner_that_cannot_read_gzip_receives_plain_copies(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    """An aligner without reads_gzip gets decompressed copies from scratch;
    the copies are removed after the alignment, and the tree keeps the names."""
    import gzip

    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        with gzip.open(reps / f"Fam_gen_sp_GCA_00000{i}.fasta.gz", "wt") as fh:
            fh.write(f">s{i}\nACGTACGT\n")
    seen: list[tuple[str, bytes]] = []
    original = _FakeAligner.align

    def recording(self, genomes, reference, out_dir, params, logger):
        seen.extend((str(g), Path(g).read_bytes()[:2]) for g in genomes)
        assert reference in genomes
        return original(self, genomes, reference, out_dir, params, logger)

    monkeypatch.setattr(_FakeAligner, "align", recording)
    ctx = WorkdirContext(workdir)
    tree = run(
        ctx,
        PhyloParams(
            no_outgroup=True,
            treebuilder="faketree_msa",
            msa_source="aligner",
            aligner="fakealigner",
        ),
    )
    assert len(seen) == 3
    assert all(name.endswith(".fasta") and "phylo_inputs" in name for name, _ in seen)
    assert all(head == b">s" for _, head in seen)
    assert not (workdir / "scratch" / "phylo_inputs").exists()
    assert "Fam_gen_sp_GCA_000001" in tree.read_text()
    assert ".fasta" not in (workdir / "align" / "msa.fasta").read_text()


# -- the stage harness: which failures leave an [interrupted] record ---------


def _phylo_record_via_cli(workdir: Path, params: PhyloParams) -> dict | None:
    import typer
    import yaml

    from repgenr.cli import base as cli
    from repgenr.core.config import CONFIG_FILENAME

    with pytest.raises(typer.Exit):
        cli._run("phylo", workdir, lambda: params)
    config = workdir / CONFIG_FILENAME
    if not config.exists():
        return None
    data = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
    return (data.get("stages") or {}).get("phylo")


class _FailingTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_fail")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        from repgenr.core.errors import ToolExecutionError

        raise ToolExecutionError(["faketree", "build"], 1, "boom")


def test_clean_refusal_writing_nothing_leaves_no_record(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    from repgenr.cli import base as cli

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 3):
        (reps / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGTACGT\n")
    record = _phylo_record_via_cli(
        workdir, PhyloParams(treebuilder="faketree_genomes", no_outgroup=True)
    )
    assert record is None


def test_tool_failure_leaves_an_interrupted_record(
    workdir: Path, fake_phylo_tools, monkeypatch, register_tool
) -> None:
    from repgenr.cli import base as cli

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    register_tool(tb_registry, "faketree_fail", _FailingTreeBuilder)
    _make_reps(workdir)
    record = _phylo_record_via_cli(
        workdir, PhyloParams(treebuilder="faketree_fail", no_outgroup=True)
    )
    assert record is not None and not record.get("completed")


def test_leaf_check_failure_leaves_an_interrupted_record(
    workdir: Path, fake_phylo_tools, monkeypatch, register_tool
) -> None:
    from repgenr.cli import base as cli

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    register_tool(tb_registry, "faketree_drop", _DroppingTreeBuilder)
    _make_reps(workdir)
    record = _phylo_record_via_cli(
        workdir, PhyloParams(treebuilder="faketree_drop", no_outgroup=True)
    )
    assert record is not None and not record.get("completed")
    assert (workdir / "tree" / "tree.nwk").exists()


def test_gzipped_genomes_and_outgroup_are_named_without_their_suffix(
    workdir: Path, fake_phylo_tools, fake_snptyper, monkeypatch
) -> None:
    """Leaves, the outgroup leaf and the masker's exclusion use the name without .fasta.gz.

    The typer names records as the simple typer does (strip_fasta_suffix), and
    tree2tax and clusters.tsv use the same form.
    """
    import gzip

    from repgenr.core.contracts import strip_fasta_suffix
    from repgenr.stages import snptype as snptype_stage

    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        with gzip.open(reps / f"Fam_gen_sp_GCA_00000{i}.fasta.gz", "wt") as fh:
            fh.write(f">s{i}\nACGTACGT\n")
    outgroup = workdir / "outgroup"
    outgroup.mkdir()
    with gzip.open(outgroup / "Fam_gen_og_GCA_000009.fasta.gz", "wt") as fh:
        fh.write(">og\nACGTACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009\n")

    def named_call(self, genomes, reference, out_dir, params, logger) -> SnpResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        core = out_dir / "core.fasta"
        core.write_text("".join(f">{strip_fasta_suffix(Path(g).name)}\nACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core)

    monkeypatch.setattr(_FakeSnpTyper, "call", named_call)
    seen: dict = {}
    original_core = snptype_stage.snptype_core

    def recording_core(genomes, reference, snp_dir, scratch, params, logger, **kw):
        seen["mask_exclude"] = params.mask_exclude
        return original_core(genomes, reference, snp_dir, scratch, params, logger, **kw)

    monkeypatch.setattr(snptype_stage, "snptype_core", recording_core)
    original_build = _MsaTreeBuilder.build

    def recording_build(self, msa, out_dir, params, logger):
        seen["outgroup"] = params.outgroup
        return original_build(self, msa, out_dir, params, logger)

    monkeypatch.setattr(_MsaTreeBuilder, "build", recording_build)

    ctx = WorkdirContext(workdir, create=True)
    tree = run(
        ctx,
        PhyloParams(treebuilder="faketree_msa", msa_source="snptype", snptyper="fakesnptyper"),
    )
    text = tree.read_text()
    assert ".fasta" not in text
    assert "Fam_gen_sp_GCA_000001" in text and "Fam_gen_og_GCA_000009" in text
    assert seen["outgroup"] == "Fam_gen_og_GCA_000009"
    # One name form: every typer names records by record_name.
    assert seen["mask_exclude"] == ("Fam_gen_og_GCA_000009",)


class _NoSupportMsaTreeBuilder(_MsaTreeBuilder):
    """An MSA builder whose tree has no support values unless a bootstrap runs."""

    capabilities = ToolCapabilities(name="faketree_nosupport")
    supports_without_bootstrap = False


def _stage_outgroup(workdir: Path) -> None:
    og = workdir / "outgroup"
    og.mkdir()
    (og / "GCF_000009.1.fasta").write_text(">o\nACGTACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCF_000009.1\n")


def _spy_auto_select(monkeypatch, result: str) -> dict:
    """Replace the registry-wide choice with ``result`` and record what the
    stage asked for; the choice itself is covered by the auto_select unit tests."""
    from repgenr.stages import phylo as phylo_mod

    seen: dict = {}

    def fake(registry, n_items, allowed=None):
        seen.update(n=n_items, allowed=None if allowed is None else sorted(allowed))
        return result

    monkeypatch.setattr(phylo_mod, "auto_select", fake)
    return seen


def test_auto_restricts_to_alignment_builders_when_a_snp_source_is_requested(
    workdir: Path, fake_phylo_tools, fake_snptyper, monkeypatch
) -> None:
    """'auto' would otherwise hand the set to an alignment-free builder and
    drop --msa-source snptype with only a warning."""
    seen = _spy_auto_select(monkeypatch, "faketree_msa")
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    run(
        ctx,
        PhyloParams(
            treebuilder="auto", msa_source="snptype", snptyper="fakesnptyper", no_outgroup=True
        ),
    )
    assert "faketree_msa" in seen["allowed"] and "faketree_genomes" not in seen["allowed"]
    assert ctx.config.stages["phylo"].tool == "faketree_msa"


def test_auto_without_an_alignment_request_considers_every_builder(
    workdir: Path, fake_phylo_tools, monkeypatch
) -> None:
    seen = _spy_auto_select(monkeypatch, "faketree_genomes")
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, PhyloParams(treebuilder="auto", no_outgroup=True))
    assert seen["allowed"] is None and seen["n"] == 3
    assert ctx.config.stages["phylo"].tool == "faketree_genomes"


def test_auto_counts_the_outgroup_as_a_leaf(workdir: Path, fake_phylo_tools, monkeypatch) -> None:
    seen = _spy_auto_select(monkeypatch, "faketree_genomes")
    _make_reps(workdir)
    _stage_outgroup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, PhyloParams(treebuilder="auto"))
    assert seen["n"] == 4


def test_bootstrap_zero_warns_when_the_tree_would_carry_no_supports(
    workdir: Path, fake_phylo_tools, register_tool, caplog
) -> None:
    register_tool(tb_registry, "faketree_nosupport", _NoSupportMsaTreeBuilder)
    _make_reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-phylo"))
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_nosupport",
                msa_source="aligner",
                aligner="fakealigner",
                no_outgroup=True,
            ),
        )
    assert any("No bootstrap requested" in r.getMessage() for r in caplog.records)
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="test-phylo"):
        run(
            ctx,
            PhyloParams(
                treebuilder="faketree_msa",
                msa_source="aligner",
                aligner="fakealigner",
                no_outgroup=True,
            ),
        )
    assert not any("No bootstrap requested" in r.getMessage() for r in caplog.records)


def test_split_alignment_refusal_names_the_builders_that_consume_one(
    workdir: Path, fake_phylo_tools
) -> None:
    from repgenr.core.errors import UserInputError
    from repgenr.stages.phylo import PhyloDirs, build_msa

    _make_reps(workdir)
    genomes = sorted((workdir / "derep" / "representatives").glob("*.fasta"))
    dirs = PhyloDirs(
        tree_dir=workdir / "tree",
        align_dir=workdir / "align",
        msa_dir=workdir / "tree" / "msa",
        scratch_dir=workdir / "scratch",
    )
    with pytest.raises(UserInputError, match="faketree_msa"):
        build_msa(
            genomes,
            None,
            dirs,
            PhyloParams(treebuilder="faketree_genomes"),
            logging.getLogger("test-phylo"),
        )
