"""Gubbins is invoked on the whole-genome alignment with a thread budget."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.maskers import gubbins as mod
from repgenr.maskers.base import MaskParams


@pytest.fixture(autouse=True)
def _many_cpus(monkeypatch):
    # The thread cap reads the CPU count where Gubbins runs; pin it so the
    # requested threads pass through on any host (CI runners have 4 CPUs).
    monkeypatch.setattr(mod, "available_cpus", lambda caps: 64)


def test_gubbins_argv(tmp_path: Path, monkeypatch) -> None:
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        out_prefix = str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta"
        Path(out_prefix).write_text(">a\nA\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n", encoding="utf-8")
    masker = mod.GubbinsMasker()
    out = masker.mask(full, tmp_path / "gub", MaskParams(threads=4), logging.getLogger("t"))
    assert out.exists()
    argv = calls[0]
    assert argv[0] == "run_gubbins.py"
    assert argv[argv.index("--threads") + 1] == "4"
    assert "--tree-builder" not in argv, "Gubbins' own default stands when RAxML can run"
    assert argv[-1] == str(full)


def test_extras_select_gubbins_tree_builders(tmp_path: Path, monkeypatch) -> None:
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        Path(str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta").write_text(
            ">a\nA\n", encoding="utf-8"
        )

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: False)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n", encoding="utf-8")
    params = MaskParams(
        threads=8,
        extra={
            "gubbins_tree_builder": "fasttree",
            "gubbins_first_tree_builder": "rapidnj",
            "gubbins_args": "--min-snps 5 --iterations 3",
        },
    )
    mod.GubbinsMasker().mask(full, tmp_path / "gub", params, logging.getLogger("t"))
    argv = calls[0]
    assert argv[argv.index("--tree-builder") + 1] == "fasttree"
    assert argv[argv.index("--first-tree-builder") + 1] == "rapidnj"
    assert argv[argv.index("--min-snps") + 1] == "5"
    assert argv[argv.index("--iterations") + 1] == "3"
    assert argv[argv.index("--threads") + 1] == "8", "a chosen builder keeps the thread budget"
    assert argv[-2] == "--prefix" or argv[-3] == "--prefix"


def test_resolve_tree_builder_falls_back_without_threaded_raxml(monkeypatch, caplog) -> None:
    """Gubbins exits when asked for threads without a PTHREADS RAxML build."""
    logger = logging.getLogger("t")
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: False)

    monkeypatch.setattr(
        mod.shutil, "which", lambda name, path=None: "/bin/x" if name == "iqtree2" else None
    )
    with caplog.at_level(logging.WARNING):
        assert mod.resolve_tree_builder(None, 8, logger, on_host=True) == ("iqtree", 8)
    assert "iqtree" in caplog.text

    monkeypatch.setattr(mod.shutil, "which", lambda name, path=None: None)
    assert mod.resolve_tree_builder(None, 8, logger, on_host=True) == (None, 1)

    # One thread, a container run, or an explicit choice never trigger it.
    assert mod.resolve_tree_builder(None, 1, logger, on_host=True) == (None, 1)
    assert mod.resolve_tree_builder(None, 8, logger, on_host=False) == (None, 8)
    assert mod.resolve_tree_builder("raxmlng", 8, logger, on_host=True) == ("raxmlng", 8)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    assert mod.resolve_tree_builder(None, 8, logger, on_host=True) == (None, 8)


def test_sanitise_replaces_iupac_with_n(tmp_path, caplog) -> None:
    """IUPAC codes from a diploid-style consensus become N before Gubbins."""
    import logging

    from repgenr.maskers.gubbins import sanitise_alignment

    src = tmp_path / "aln.fasta"
    src.write_text(">a\nACGTRYACGT\n>b\nACGTACGT-N\n", encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        out = sanitise_alignment(src, tmp_path / "clean.fasta", logging.getLogger("t"))
    assert out.read_text(encoding="utf-8") == ">a\nACGTNNACGT\n>b\nACGTACGT-N\n"
    assert "2 ambiguous base(s)" in caplog.text


def test_sanitise_keeps_a_clean_alignment(tmp_path) -> None:
    import logging

    from repgenr.maskers.gubbins import sanitise_alignment

    src = tmp_path / "aln.fasta"
    src.write_text(">a\nACGT\n", encoding="utf-8")
    assert sanitise_alignment(src, tmp_path / "clean.fasta", logging.getLogger("t")) == src
    assert not (tmp_path / "clean.fasta").exists()


def test_exclude_runs_gubbins_on_the_ingroup_and_masks_everyone(
    tmp_path: Path, monkeypatch
) -> None:
    """D-10: the outgroup stays out of the scan but keeps its place in the output."""
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        prefix = str(kw["cwd"] / "gubbins")
        # Gubbins predicts a recombinant block in `a` over columns 2-3.
        Path(prefix + ".recombination_predictions.gff").write_text(
            "##gff-version 3\n"
            'SEQUENCE\tGUBBINS\tCDS\t2\t3\t0.0\t.\t.\tnode="a";taxa="  a";snp_count="2"\n',
            encoding="utf-8",
        )

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nATTTA\n>og\nGGGGG\n", encoding="utf-8")
    out = mod.GubbinsMasker().mask(
        full,
        tmp_path / "gub",
        MaskParams(threads=2, exclude=frozenset({"og"})),
        logging.getLogger("t"),
    )
    argv = calls[0]
    ingroup = mod.read_fasta(Path(argv[-1]))
    assert set(ingroup) == {"a", "b", "c"}, "the outgroup is not handed to Gubbins"
    masked = mod.read_fasta(out)
    assert set(masked) == {"a", "b", "c", "og"}, "the outgroup is in the masked alignment"
    # Column 1 (A/G) and column 5 (A/G) vary; columns 2-3 of `a` are N so only
    # b/c/og decide those: T vs G still varies. Column 4: T/T/T/G varies.
    assert masked["og"] == "GGGGG"
    assert masked["a"].startswith("A") and "N" in masked["a"]


def test_gff_regions_and_polymorphic_sites(tmp_path: Path) -> None:
    gff = tmp_path / "p.gff"
    gff.write_text(
        'SEQ\tGUBBINS\tCDS\t1\t2\t0\t.\t.\tnode="x";taxa="  a  b";snp_count="1"\n'
        'SEQ\tGUBBINS\tCDS\t4\t4\t0\t.\t.\tnode="y";taxa="  b";snp_count="1"\n',
        encoding="utf-8",
    )
    regions = mod.read_recombination_gff(gff)
    assert regions == {"a": [(1, 2)], "b": [(1, 2), (4, 4)]}
    masked = mod.apply_masks({"a": "ACGT", "b": "ACGT", "c": "TCGA"}, regions)
    assert masked == {"a": "NNGT", "b": "NNGN", "c": "TCGA"}
    sites = mod.polymorphic_sites(masked)
    assert sites == {"a": "T", "b": "N", "c": "A"}, "only column 4 still varies among ACGT"


def test_variable_fraction_estimates_divergence() -> None:
    """One varying column in eight, whatever the sampling stride."""
    records = {"a": "ACGTACGT", "b": "ACGTACGA", "c": "ACGTACGT"}
    fraction, length = mod.variable_fraction(records)
    assert (round(fraction, 3), length) == (0.125, 8)
    assert mod.variable_fraction({"a": "ACGT"}) == (0.0, 4)
    assert mod.variable_fraction({}) == (0.0, 0)


def test_divergent_alignment_warns_before_gubbins_runs(tmp_path: Path, monkeypatch, caplog) -> None:
    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        Path(str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta").write_text(
            ">a\nA\n>b\nT\n", encoding="utf-8"
        )

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    full = tmp_path / "full.fasta"
    full.write_text(_alignment(variable=110, total=1000), encoding="utf-8")  # 11%: warn only
    with caplog.at_level(logging.WARNING):
        mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    assert "variable" in caplog.text and "--mask none" in caplog.text


def _alignment(*, variable: int, total: int, names: tuple[str, str] = ("a", "b")) -> str:
    """Two records over ``total`` columns of which ``variable`` differ."""
    a = "A" * total
    b = "T" * variable + "A" * (total - variable)
    return f">{names[0]}\n{a}\n>{names[1]}\n{b}\n"


def _fake_gubbins(monkeypatch, calls: list | None = None, with_tree: bool = False):
    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        if calls is not None:
            calls.append([str(a) for a in argv])
        prefix = str(kw["cwd"] / "gubbins")
        Path(prefix + ".filtered_polymorphic_sites.fasta").write_text(
            ">a\nA\n>b\nT\n", encoding="utf-8"
        )
        if with_tree:
            Path(prefix + ".final_tree.tre").write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)


def test_too_divergent_alignment_is_refused_before_gubbins_runs(
    tmp_path: Path, monkeypatch
) -> None:
    from repgenr.core.errors import UserInputError

    calls: list = []
    _fake_gubbins(monkeypatch, calls)
    full = tmp_path / "full.fasta"
    full.write_text(_alignment(variable=200, total=1000), encoding="utf-8")  # 20%
    with pytest.raises(UserInputError, match="20% of the alignment is variable"):
        mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    assert calls == []
    # A higher limit or gubbins_force lets it run.
    mod.GubbinsMasker().mask(
        full,
        tmp_path / "gub2",
        MaskParams(extra={"gubbins_divergence_limit": "0.25"}),
        logging.getLogger("t"),
    )
    mod.GubbinsMasker().mask(
        full, tmp_path / "gub3", MaskParams(extra={"gubbins_force": "true"}), logging.getLogger("t")
    )
    assert len(calls) == 2
    with pytest.raises(UserInputError, match="gubbins_divergence_limit"):
        mod.GubbinsMasker().mask(
            full,
            tmp_path / "gub4",
            MaskParams(extra={"gubbins_divergence_limit": "x"}),
            logging.getLogger("t"),
        )


def test_a_short_alignment_is_neither_refused_nor_warned(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    _fake_gubbins(monkeypatch)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTACGT\n>b\nTGCATGCA\n", encoding="utf-8")  # 8 columns, all variable
    with caplog.at_level(logging.WARNING):
        mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    assert "variable" not in caplog.text


def test_more_than_one_species_is_refused_unless_forced(tmp_path: Path, monkeypatch) -> None:
    from repgenr.core.errors import UserInputError

    calls: list = []
    _fake_gubbins(monkeypatch, calls)
    full = tmp_path / "full.fasta"
    names = ("Fam_Gen_alpha_GCF_1.1", "Fam_Gen_beta_GCF_2.1")
    full.write_text(_alignment(variable=10, total=1000, names=names), encoding="utf-8")
    with pytest.raises(UserInputError, match="spans 2 species"):
        mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    mod.GubbinsMasker().mask(
        full, tmp_path / "gub2", MaskParams(extra={"gubbins_force": "true"}), logging.getLogger("t")
    )
    assert len(calls) == 1


def test_first_tree_defaults_to_rapidnj_when_available(tmp_path: Path, monkeypatch) -> None:
    calls: list = []
    _fake_gubbins(monkeypatch, calls)
    monkeypatch.setattr(mod, "runs_on_host", lambda caps: True)
    monkeypatch.setattr(mod, "host_which", lambda name, tool: "/usr/bin/rapidnj")
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    argv = calls[0]
    assert argv[argv.index("--first-tree-builder") + 1] == "rapidnj"
    assert argv[argv.index("--first-model") + 1] == "JC"
    # An explicit choice, or the flag inside gubbins_args, stands.
    mod.GubbinsMasker().mask(
        full,
        tmp_path / "gub2",
        MaskParams(extra={"gubbins_first_tree_builder": "fasttree"}),
        logging.getLogger("t"),
    )
    assert calls[1][calls[1].index("--first-tree-builder") + 1] == "fasttree"
    assert "--first-model" not in calls[1]
    mod.GubbinsMasker().mask(
        full,
        tmp_path / "gub3",
        MaskParams(extra={"gubbins_args": "--first-tree-builder iqtree --first-model GTR"}),
        logging.getLogger("t"),
    )
    assert calls[2].count("--first-tree-builder") == 1 and calls[2].count("--first-model") == 1
    # Without rapidnj on the host, Gubbins keeps its own default.
    monkeypatch.setattr(mod, "host_which", lambda name, tool: None)
    mod.GubbinsMasker().mask(full, tmp_path / "gub4", MaskParams(), logging.getLogger("t"))
    assert "--first-tree-builder" not in calls[3]


def test_gubbins_tree_is_offered_after_the_run(tmp_path: Path, monkeypatch) -> None:
    _fake_gubbins(monkeypatch, with_tree=True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    masker = mod.GubbinsMasker()
    assert masker.tree_output() is None
    masker.mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    tree = masker.tree_output()
    assert tree is not None and tree.read_text() == "(a,b);\n"


def test_gubbins_failure_reports_the_divergence(tmp_path: Path, monkeypatch) -> None:
    """A crash in the scan is reported with the figure that explains it."""
    from repgenr.core.errors import ToolExecutionError

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        raise ToolExecutionError(["run_gubbins.py"], 1, "Bus error")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTACGT\n>b\nTGCATGCA\n", encoding="utf-8")
    with pytest.raises(ToolExecutionError) as ei:
        mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    assert "within-species" in ei.value.details()
    assert ei.value.returncode == 1  # the tool's status survives for the retry rule


def _fake_gubbins_keeping(kept: set[str], calls: list):
    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        records = mod.read_fasta(Path(argv[-1]))
        prefix = str(kw["cwd"] / "gubbins")
        Path(prefix + ".filtered_polymorphic_sites.fasta").write_text(
            "".join(f">{n}\n{s}\n" for n, s in records.items() if n in kept),
            encoding="utf-8",
        )
        Path(prefix + ".recombination_predictions.gff").write_text(
            "##gff-version 3\n", encoding="utf-8"
        )

    return fake_run_tool


def test_gubbins_keeps_taxa_that_are_mostly_n(tmp_path: Path, monkeypatch) -> None:
    """The simple typer writes N where a genome does not align; Gubbins must keep it."""
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_gubbins_keeping({"a", "b", "c"}, calls))
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nANNNN\n", encoding="utf-8")
    mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(), logging.getLogger("t"))
    argv = calls[0]
    assert argv[argv.index("--filter-percentage") + 1] == "100"

    calls.clear()
    params = MaskParams(extra={"gubbins_args": "--filter-percentage 90"})
    mod.GubbinsMasker().mask(full, tmp_path / "gub2", params, logging.getLogger("t"))
    assert calls[0].count("--filter-percentage") == 1, "the user's choice stands"


@pytest.mark.parametrize("exclude", [frozenset(), frozenset({"og"})])
@pytest.mark.parametrize(
    "gubbins_args", ["--filter-percentage 25", "--filter-percentage=25", "-f 25"]
)
def test_taxa_a_user_filter_would_leave_out_are_refused(
    tmp_path: Path, monkeypatch, exclude, gubbins_args
) -> None:
    """Gubbins would skip such a taxon in its scan yet keep it in its output."""
    from repgenr.core.errors import WorkdirError

    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_gubbins_keeping({"a", "b", "c", "og"}, calls))
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nANNNN\n>og\nGGGGG\n", encoding="utf-8")
    params = MaskParams(exclude=exclude, extra={"gubbins_args": gubbins_args})
    with pytest.raises(WorkdirError) as exc:
        mod.GubbinsMasker().mask(full, tmp_path / "gub", params, logging.getLogger("t"))
    assert ": c." in str(exc.value) and "--filter-percentage" in str(exc.value)
    assert exc.value.exit_code == 3
    assert calls == [], "refused before Gubbins runs"


def test_a_user_filter_that_keeps_every_taxon_runs(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_gubbins_keeping({"a", "b", "c"}, calls))
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nANNNA\n", encoding="utf-8")
    params = MaskParams(extra={"gubbins_args": "--filter-percentage 80"})
    mod.GubbinsMasker().mask(full, tmp_path / "gub", params, logging.getLogger("t"))
    assert len(calls) == 1


def test_a_non_numeric_filter_percentage_is_refused(tmp_path: Path) -> None:
    from repgenr.core.errors import UserInputError

    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n", encoding="utf-8")
    params = MaskParams(extra={"gubbins_args": "--filter-percentage lots"})
    with pytest.raises(UserInputError):
        mod.GubbinsMasker().mask(full, tmp_path / "gub", params, logging.getLogger("t"))


def test_gubbins_threads_are_capped_at_the_available_cpus(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    # The global default is 16 threads; IQ-TREE as Gubbins' tree builder
    # refuses more threads than cores ("more threads than CPU cores available").
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        out_prefix = str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta"
        Path(out_prefix).write_text(">a\nA\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda path=None: True)
    monkeypatch.setattr(mod, "available_cpus", lambda caps: 11)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n", encoding="utf-8")
    with caplog.at_level(logging.INFO):
        mod.GubbinsMasker().mask(
            full, tmp_path / "gub", MaskParams(threads=16), logging.getLogger("t")
        )
    argv = calls[0]
    assert argv[argv.index("--threads") + 1] == "11"
    assert any("16 threads" in r.getMessage() and "11" in r.getMessage() for r in caplog.records)


def test_available_cpus_is_the_engine_count_under_docker(monkeypatch) -> None:
    import subprocess

    from repgenr.core import containers
    from repgenr.core.plugins import ToolCapabilities

    caps = ToolCapabilities(name="gubbins", container="quay.io/x/gubbins:1")
    monkeypatch.setattr(containers, "_CONFIG", containers.ContainerConfig(backend="docker"))
    monkeypatch.setattr(
        containers.subprocess,
        "run",
        lambda argv, **kw: subprocess.CompletedProcess(argv, 0, "6\n", ""),
    )
    containers._ENGINE_CPUS.clear()
    assert containers.available_cpus(caps) == 6
    monkeypatch.setattr(containers, "_CONFIG", containers.ContainerConfig())
    assert containers.available_cpus(caps) >= 1


def test_exclude_name_matching_no_record_is_warned(tmp_path: Path, monkeypatch, caplog) -> None:
    """A typer that names the outgroup 'og.fasta' (Path.stem of og.fasta.gz)
    must not have Gubbins scan the outgroup without a warning."""
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        prefix = str(kw["cwd"] / "gubbins")
        Path(prefix + ".recombination_predictions.gff").write_text(
            "##gff-version 3\n", encoding="utf-8"
        )
        # Nothing was left out, so the masker reads Gubbins' own filtered sites.
        Path(prefix + ".filtered_polymorphic_sites.fasta").write_text(">a\nA\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nATTTA\n>og.fasta\nGGGGG\n", encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        mod.GubbinsMasker().mask(
            full,
            tmp_path / "gub",
            MaskParams(threads=2, exclude=frozenset({"og"})),
            logging.getLogger("t"),
        )
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert any("no alignment record is named og" in w for w in warnings), warnings
    assert "og.fasta" in set(mod.read_fasta(Path(calls[0][-1]))), "nothing was left out"


def test_matched_exclude_name_is_not_warned(tmp_path: Path, monkeypatch, caplog) -> None:
    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        Path(str(kw["cwd"] / "gubbins") + ".recombination_predictions.gff").write_text(
            "##gff-version 3\n", encoding="utf-8"
        )

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nATTTA\n>og\nGGGGG\n", encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        mod.GubbinsMasker().mask(
            full,
            tmp_path / "gub",
            MaskParams(threads=2, exclude=frozenset({"og"})),
            logging.getLogger("t"),
        )
    assert not any("no alignment record" in r.getMessage() for r in caplog.records)


def test_species_are_counted_only_from_accession_like_names() -> None:
    assert (
        mod.species_spread(["Francisella_tularensis_LVS_1", "Francisella_tularensis_SCHU_S4"]) == 0
    )
    assert mod.species_spread(["sample_1_run_A", "sample_2_run_B"]) == 0
    assert mod.species_spread(["Fam_Gen_alpha_GCF_1.1", "Fam_Gen_beta_GCA_2.1"]) == 2
    assert mod.species_spread(["Fam_Gen_alpha_GCF_1.1", "Fam_Gen_alpha_NC_000001.1"]) == 1
    assert mod.species_spread(["Fam_Gen_alpha_SRR1", "Fam_Gen_beta_ERR22"]) == 2


def test_forced_runs_log_what_they_override(tmp_path: Path, monkeypatch, caplog) -> None:
    _fake_gubbins(monkeypatch)
    full = tmp_path / "full.fasta"
    names = ("Fam_Gen_alpha_GCF_1.1", "Fam_Gen_beta_GCF_2.1")
    full.write_text(_alignment(variable=200, total=1000, names=names), encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        mod.GubbinsMasker().mask(
            full,
            tmp_path / "gub",
            MaskParams(extra={"gubbins_force": "true"}),
            logging.getLogger("t"),
        )
    assert "spanning 2 species" in caplog.text and "gubbins_force=true" in caplog.text
    assert "otherwise refused" in caplog.text


def test_precheck_refuses_before_any_typing(caplog) -> None:
    from repgenr.core.errors import UserInputError

    masker = mod.GubbinsMasker()
    with pytest.raises(UserInputError, match="spans 2 species"):
        masker.precheck(
            ["Fam_Gen_alpha_GCF_1.1", "Fam_Gen_beta_GCF_2.1"], {}, logging.getLogger("t")
        )
    with pytest.raises(UserInputError, match="gubbins_divergence_limit"):
        masker.precheck(["a"], {"gubbins_divergence_limit": "2"}, logging.getLogger("t"))
    masker.precheck(["Fam_Gen_alpha_GCF_1.1", "Fam_Gen_alpha_GCF_2.1"], {}, logging.getLogger("t"))


def test_gubbins_force_accepts_only_booleans() -> None:
    from repgenr.core.errors import UserInputError

    for raw, expected in (
        ("true", True),
        ("1", True),
        ("yes", True),
        ("false", False),
        ("", False),
    ):
        assert mod._forced({"gubbins_force": raw}) is expected
    with pytest.raises(UserInputError, match="gubbins_force"):
        mod._forced({"gubbins_force": "ture"})


def test_a_user_first_model_suppresses_the_rapidnj_default(tmp_path: Path, monkeypatch) -> None:
    calls: list = []
    _fake_gubbins(monkeypatch, calls)
    monkeypatch.setattr(mod, "runs_on_host", lambda caps: True)
    monkeypatch.setattr(mod, "host_which", lambda name, tool: "/usr/bin/rapidnj")
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    mod.GubbinsMasker().mask(
        full,
        tmp_path / "gub",
        MaskParams(extra={"gubbins_args": "--first-model GTR"}),
        logging.getLogger("t"),
    )
    assert "--first-tree-builder" not in calls[0] and calls[0].count("--first-model") == 1
