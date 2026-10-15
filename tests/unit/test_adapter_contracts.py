"""Contract tests over every registered dereplicator and treebuilder adapter.

Each adapter runs against a faked ``run_tool`` that records argv and writes
canned tool outputs, so the public contract is pinned without the real
binaries: every input genome ends up in exactly one cluster, normalized
parameters land in argv, and tree builders return an existing Newick file.
"""

from __future__ import annotations

import logging
import shutil
import sys
from pathlib import Path

import pytest

from repgenr.dereplicators.base import DerepParams, DerepResult
from repgenr.dereplicators.base import registry as derep_registry
from repgenr.treebuilders.base import InputKind, TreeParams
from repgenr.treebuilders.base import registry as tree_registry

_LOG = logging.getLogger("test")

_GENOMES = ["g1.fasta", "g2.fasta", "g3.fasta", "g4.fasta"]
# canned clustering: g1 represents {g2}, g3 represents {g4}
_REP_OF = {
    "g1.fasta": "g1.fasta",
    "g2.fasta": "g1.fasta",
    "g3.fasta": "g3.fasta",
    "g4.fasta": "g3.fasta",
}
_NEWICK = "(g1:0.1,g2:0.1,(g3:0.1,g4:0.1):0.1);\n"


@pytest.fixture()
def genomes(tmp_path) -> list[Path]:
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    out = []
    for name in _GENOMES:
        p = gdir / name
        p.write_text(f">{Path(name).stem}\nACGTACGTACGT\n")
        out.append(p)
    return out


def _flag_value(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


def _read_fofn(path: str) -> list[Path]:
    return [Path(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line]


# --- canned tool output writers -----------------------------------------------


def _fake_skder(cmd: list[str]) -> None:
    result_dir = Path(_flag_value(cmd, "-o"))
    # One -g token: a directory staged with one entry per genome.
    genome_dir = Path(_flag_value(cmd, "-g"))
    assert genome_dir.is_dir() and cmd[cmd.index("-g") + 2] == "-o"
    genome_args = sorted(genome_dir.iterdir())
    rep_dir = result_dir / "Dereplicated_Representative_Genomes"
    rep_dir.mkdir(parents=True)
    by_name = {p.name: p for p in genome_args}
    for rep in ("g1.fasta", "g3.fasta"):
        shutil.copy2(by_name[rep], rep_dir / rep)
    edges = result_dir / "Skani_Triangle_Edge_Output.txt"
    lines = ["#query\treference\tANI\tAF_query\tAF_reference"]
    for member, rep in _REP_OF.items():
        if member != rep:
            lines.append(f"{by_name[member]}\t{by_name[rep]}\t99.5\t95.0\t95.0")
    edges.write_text("\n".join(lines) + "\n")


def _fake_drep(cmd: list[str]) -> None:
    drep_wd = Path(cmd[2])
    staged = _read_fofn(_flag_value(cmd, "-g"))
    by_name = {p.name: p for p in staged}
    rep_dir = drep_wd / "dereplicated_genomes"
    rep_dir.mkdir(parents=True)
    for rep in ("g1.fasta", "g3.fasta"):
        shutil.copy2(by_name[rep], rep_dir / rep)
    tables = drep_wd / "data_tables"
    tables.mkdir()
    rows = ["genome,secondary_cluster"]
    cluster_id = {"g1.fasta": "1_1", "g3.fasta": "2_1"}
    for member, rep in _REP_OF.items():
        rows.append(f"{member},{cluster_id[rep]}")
    (tables / "Cdb.csv").write_text("\n".join(rows) + "\n")


def _fake_galah(cmd: list[str]) -> None:
    clusters_file = Path(_flag_value(cmd, "--output-cluster-definition"))
    listed = _read_fofn(_flag_value(cmd, "--genome-fasta-list"))
    by_name = {p.name: p for p in listed}
    rows = [f"{by_name[rep]}\t{by_name[member]}" for member, rep in _REP_OF.items()]
    clusters_file.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _fake_sourmash_sketch(cmd: list[str]) -> None:
    outdir = Path(_flag_value(cmd, "--outdir"))
    for p in _read_fofn(_flag_value(cmd, "--from-file")):
        (outdir / f"{p.name}.sig").write_text("{}")


def _fake_sourmash_compare(cmd: list[str]) -> None:
    csv_path = Path(_flag_value(cmd, "--csv"))
    sigs = _read_fofn(_flag_value(cmd, "--from-file"))
    labels = [p.name.removesuffix(".sig") for p in sigs]

    def sim(a: str, b: str) -> float:
        if a == b:
            return 1.0
        return 0.995 if _REP_OF[a] == _REP_OF[b] else 0.2

    rows = [",".join(labels)]
    for a in labels:
        rows.append(",".join(f"{sim(a, b):.4f}" for b in labels))
    csv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _fake_skani_triangle(cmd: list[str]) -> None:
    paths = _read_fofn(_flag_value(cmd, "-l"))

    def dist(a: Path, b: Path) -> float:
        if a == b:
            return 0.0
        return 0.5 if _REP_OF[a.name] == _REP_OF[b.name] else 12.0

    rows = [str(len(paths))]
    for a in paths:
        rows.append("\t".join([str(a), *(f"{dist(a, b):.2f}" for b in paths)]))
    Path(_flag_value(cmd, "-o")).write_text("\n".join(rows) + "\n", encoding="utf-8")


def _make_fake_run_tool(recorded: list[list[str]]):
    def fake_run_tool(caps, command, *, logger, stdout_path=None, **kwargs):
        cmd = [str(part) for part in command]
        recorded.append(cmd)
        if cmd[:3] == ["sourmash", "scripts", "pairwise"]:
            return 1  # no branchwater plugin -> dense path
        if cmd[0] == "skder":
            _fake_skder(cmd)
        elif cmd[0] == "dRep" and cmd[1] == "dereplicate":
            _fake_drep(cmd)
        elif cmd[0] == "galah":
            _fake_galah(cmd)
        elif cmd[:2] == ["sourmash", "sketch"]:
            _fake_sourmash_sketch(cmd)
        elif cmd[:2] == ["sourmash", "compare"]:
            _fake_sourmash_compare(cmd)
        elif cmd[0] in ("mashtree", "FastTree", "fasttree", "VeryFastTree"):
        elif cmd[:2] == ["skani", "triangle"]:
            _fake_skani_triangle(cmd)
        elif cmd[0] in ("mashtree", "FastTree", "fasttree"):
            Path(stdout_path).write_text(_NEWICK, encoding="utf-8")
        elif cmd[0] == "iqtree":
            msa = Path(_flag_value(cmd, "-s"))
            Path(str(msa) + ".treefile").write_text(_NEWICK, encoding="utf-8")
        elif cmd[0] == "raxml-ng":
            prefix = _flag_value(cmd, "--prefix")
            Path(prefix + ".raxml.bestTree").write_text(_NEWICK, encoding="utf-8")
        else:
            # A third-party adapter's unknown argv must not crash the whole
            # suite; its own (skipped) parametrization covers it.
            return 0
        return 0

    return fake_run_tool


@pytest.fixture()
def recorded(monkeypatch) -> list[list[str]]:
    """Patch run_tool in every adapter module; return the recorded argvs."""
    calls: list[list[str]] = []
    fake = _make_fake_run_tool(calls)
    for reg in (derep_registry, tree_registry):
        for name in reg.names():
            module = sys.modules[reg.get(name).__module__]
            monkeypatch.setattr(module, "run_tool", fake, raising=False)
    import repgenr.dereplicators.sourmash as sm

    sm._BRANCHWATER_CACHE.clear()
    return calls


# --- dereplicator contract ----------------------------------------------------


def _assert_partition(result: DerepResult) -> None:
    covered: list[str] = []
    for rep, members in result.clusters.items():
        covered.append(rep)
        covered.extend(members)
    assert sorted(covered) == sorted(_GENOMES), "every genome in exactly one cluster"
    assert {p.name for p in result.representatives} == set(result.clusters)
    assert set(result.genome_status) == set(_GENOMES)


# tool -> tokens that must land in some recorded argv for these params
_DEREP_PARAM_TOKENS = {
    "skder": ["-i", "99", "-c", "7"],
    "drep": ["-sa", "0.99", "--processors", "7"],
    "galah": [
        "--ani",
        "99",
        "--precluster-ani",
        "90",
        "--min-aligned-fraction",
        "50",
        "--threads",
        "7",
    ],
    "sourmash": ["k=31,scaled=1000", "--processes", "7"],
}


@pytest.mark.parametrize("tool", sorted(_DEREP_PARAM_TOKENS))
def test_dereplicator_contract(tool, genomes, recorded, tmp_path) -> None:
    if tool not in derep_registry.names():
        pytest.skip(f"{tool} not registered")
    adapter = derep_registry.create(tool)
    params = DerepParams(primary_ani=0.90, secondary_ani=0.99, aligned_fraction=0.50, threads=7)
    result = adapter.dereplicate(genomes, tmp_path / "out", params, _LOG)

    _assert_partition(result)
    assert len(result.clusters) == 2
    for rep in result.representatives:
        assert rep.exists(), f"representative {rep} must exist on disk"

    flat = [tok for cmd in recorded for tok in cmd]
    for token in _DEREP_PARAM_TOKENS[tool]:
        assert token in flat, f"{token!r} missing from recorded argv for {tool}"


def test_drep_virus_mode_passes_one_secondary_algorithm(genomes, recorded, tmp_path) -> None:
    """Virus mode selects ANImf once; a user's S_algorithm tool-arg still wins."""
    adapter = derep_registry.create("drep")
    adapter.dereplicate(genomes, tmp_path / "v", DerepParams(extra={"virus": "1"}), _LOG)
    cmd = next(c for c in recorded if c[:2] == ["dRep", "dereplicate"])
    assert cmd.count("--S_algorithm") == 1
    assert cmd[cmd.index("--S_algorithm") + 1] == "ANImf"

    recorded.clear()
    adapter.dereplicate(
        genomes, tmp_path / "u", DerepParams(extra={"virus": "1", "S_algorithm": "ANIn"}), _LOG
    )
    cmd = next(c for c in recorded if c[:2] == ["dRep", "dereplicate"])
    assert cmd.count("--S_algorithm") == 1
    assert cmd[cmd.index("--S_algorithm") + 1] == "ANIn"


def test_galah_returns_the_input_genomes_without_copying(genomes, recorded, tmp_path) -> None:
    """galah's representatives are the input files; no copy is written to out_dir."""
    import repgenr.dereplicators.galah as galah_mod

    out = tmp_path / "out"
    result = galah_mod.GalahDereplicator().dereplicate(genomes, out, DerepParams(), _LOG)
    assert result.representatives
    assert set(result.representatives) <= set(genomes)
    assert not (out / "representatives").exists()


def test_galah_empty_clusters_yields_empty_result(genomes, recorded, tmp_path, monkeypatch) -> None:
    """Pinned as-is: an empty galah clusters.tsv produces an empty DerepResult."""
    import repgenr.dereplicators.galah as galah_mod

    def empty_galah(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        Path(_flag_value(cmd, "--output-cluster-definition")).write_text("", encoding="utf-8")
        return 0

    monkeypatch.setattr(galah_mod, "run_tool", empty_galah)
    result = galah_mod.GalahDereplicator().dereplicate(
        genomes, tmp_path / "out", DerepParams(), _LOG
    )
    assert result.representatives == []
    assert result.clusters == {}
    assert result.genome_status == {}


# --- treebuilder contract -----------------------------------------------------


_TREE_PARAM_TOKENS = {
    "mashtree": ["--numcpus", "5"],
    "fasttree": ["-nt", "-gtr"],
    "veryfasttree": ["-nt", "-gtr", "-threads", "5"],
    "iqtree": ["--threads-max", "5"],
    "raxmlng": ["--threads", "auto{5}"],
    "sourmash": ["k=31,scaled=1000", "--processes", "5"],
    "skani": ["triangle", "-t", "5", "--full-matrix", "--distance"],
}


@pytest.fixture()
def msa(tmp_path) -> Path:
    path = tmp_path / "input_msa.fasta"
    path.write_text("".join(f">{Path(n).stem}\nACGTACGT\n" for n in _GENOMES))
    return path


@pytest.mark.parametrize("tool", sorted(_TREE_PARAM_TOKENS))
def test_treebuilder_contract(tool, genomes, msa, recorded, tmp_path) -> None:
    if tool not in tree_registry.names():
        pytest.skip(f"{tool} not registered")
    builder = tree_registry.create(tool)
    source = genomes if builder.input_kind is InputKind.GENOMES else msa
    tree = builder.build(source, tmp_path / "tree_out", TreeParams(threads=5), _LOG)

    assert tree.exists()
    content = tree.read_text()
    assert content.rstrip().endswith(";"), "output must be a Newick tree"
    for name in _GENOMES:
        assert Path(name).stem in content

    flat = [tok for cmd in recorded for tok in cmd]
    for token in _TREE_PARAM_TOKENS[tool]:
        assert token in flat, f"{token!r} missing from recorded argv for {tool}"


def test_every_registered_dereplicator_has_contract_coverage() -> None:
    """A new in-tree adapter must add itself to the contract parametrization."""
    assert set(derep_registry.names()) >= set(_DEREP_PARAM_TOKENS)
    builtin = {"drep", "galah", "skder", "sourmash"}
    assert builtin & set(derep_registry.names()) <= set(_DEREP_PARAM_TOKENS)


def test_every_registered_treebuilder_has_contract_coverage() -> None:
    builtin = {"fasttree", "iqtree", "mashtree", "raxmlng", "sourmash"}
    assert builtin & set(tree_registry.names()) <= set(_TREE_PARAM_TOKENS)


def test_drep_exit_0_without_results_is_a_tool_failure(genomes, tmp_path, monkeypatch) -> None:
    """dRep without CheckM logs the problem, exits 0 and writes no Cdb.csv."""
    import pytest

    import repgenr.dereplicators.drep as drep_mod
    from repgenr.core.errors import ToolExecutionError

    def drep_without_checkm(caps, command, *, logger, **kwargs):
        wd = Path(str(command[2]))
        (wd / "data_tables").mkdir(parents=True)
        return 0

    monkeypatch.setattr(drep_mod, "run_tool", drep_without_checkm)
    with pytest.raises(ToolExecutionError) as info:
        drep_mod.DrepDereplicator().dereplicate(genomes, tmp_path / "out", DerepParams(), _LOG)
    assert info.value.exit_code == 6
    assert str(info.value) == "dRep exited 0 without writing its results"
    assert "CheckM" in info.value.details()


def test_drep_reports_gzipped_genomes_under_their_input_names(tmp_path, monkeypatch) -> None:
    """dRep sees a decompressed copy (x.fasta); results name the input (x.fasta.gz)."""
    import gzip

    import repgenr.dereplicators.drep as drep_mod

    gdir = tmp_path / "genomes"
    gdir.mkdir()
    genomes = []
    for name in ("g1.fasta", "g2.fasta.gz", "g3.fasta.gz", "g4.fasta"):
        p = gdir / name
        text = f">{name}\nACGTACGTACGT\n".encode()
        p.write_bytes(gzip.compress(text) if name.endswith(".gz") else text)
        genomes.append(p)

    def fake(caps, command, *, logger, **kwargs):
        _fake_drep([str(c) for c in command])
        return 0

    monkeypatch.setattr(drep_mod, "run_tool", fake)
    result = drep_mod.DrepDereplicator().dereplicate(genomes, tmp_path / "out", DerepParams(), _LOG)
    assert set(result.genome_status) == {g.name for g in genomes}
    assert result.clusters == {"g1.fasta": ["g2.fasta.gz"], "g3.fasta.gz": ["g4.fasta"]}
    assert sorted(p.name for p in result.representatives) == ["g1.fasta", "g3.fasta.gz"]


def test_drep_marks_genomes_its_filter_dropped_as_fail_qc(genomes, tmp_path, monkeypatch) -> None:
    """dRep's --length or quality filter leaves a genome out of every table."""
    import repgenr.dereplicators.drep as drep_mod

    extra = genomes[0].parent / "short.fasta"
    extra.write_text(">short\nACGT\n")

    def fake(caps, command, *, logger, **kwargs):
        _fake_drep([str(c) for c in command])  # writes only g1..g4
        return 0

    monkeypatch.setattr(drep_mod, "run_tool", fake)
    result = drep_mod.DrepDereplicator().dereplicate(
        [*genomes, extra], tmp_path / "out", DerepParams(), _LOG
    )
    assert result.genome_status["short.fasta"] == "fail_qc"
    assert all("short.fasta" not in m for m in result.clusters.values())


@pytest.mark.parametrize(
    "label",
    ["/p/x.fasta.gz", "/p/x.fasta", "x.fna.gz", "x"],
)
def test_sourmash_tree_labels_use_the_record_name(label: str) -> None:
    """sourmash labels signatures by the input path; the leaf is the genome's
    record name, so a gzipped genome is 'x' as in every alignment."""
    from repgenr.treebuilders.sourmash import _label_to_genome

    genomes = [Path("/p/x.fasta.gz"), Path("/p/y.fasta")]
    assert _label_to_genome(label, genomes) == "x"
