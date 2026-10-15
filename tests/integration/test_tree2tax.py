"""tree2tax on a small tree (dendropy port): edges, genomes_map, outgroup rooting."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.stages.tree2tax import Tree2taxParams
from repgenr.stages.tree2tax import run as tree2tax_run

_NWK = (
    "((Fam_Gen_sp_GCA_000001.1,Fam_Gen_sp_GCA_000002.1),"
    "(Fam_Gen_sp_GCA_000003.1,Fam_Gen_sp_GCA_000004.1));"
)


def _edges(path: Path) -> set[tuple[str, str]]:
    out = set()
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        c, p = line.split("\t")
        out.add((c, p))
    return out


def _setup(workdir: Path, *, with_outgroup: bool = False, clusters: dict | None = None) -> None:
    (workdir / "tree").mkdir(parents=True)
    (workdir / "tree" / "tree.nwk").write_text(_NWK + "\n")
    if clusters is not None:
        (workdir / "derep").mkdir(parents=True)
        from repgenr.core.contracts import write_clusters

        write_clusters(workdir / "derep" / "clusters.tsv", clusters)
    if with_outgroup:
        og = workdir / "outgroup"
        og.mkdir(parents=True)
        (og / "Fam_Gen_sp_GCA_000004.1.fasta").write_text(">x\nACGT\n")
        (workdir / "outgroup_accession.txt").write_text("GCA_000004.1\n")


def test_tree2tax_unrooted_edges_and_map(workdir: Path) -> None:
    _setup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    t2t, gmap = tree2tax_run(ctx, Tree2taxParams())
    edges = _edges(t2t)
    # every leaf has a path to root; two cherries -> two distinct internal nodes
    leaves = {c for c, _ in edges if c.startswith("Fam_")}
    assert leaves == {
        "Fam_Gen_sp_GCA_000001.1",
        "Fam_Gen_sp_GCA_000002.1",
        "Fam_Gen_sp_GCA_000003.1",
        "Fam_Gen_sp_GCA_000004.1",
    }
    assert any(p == "root" for _, p in edges)
    # genomes_map: each accession -> its leaf
    lines = gmap.read_text().splitlines()
    accs = {ln.split("\t")[0] for ln in lines}
    assert accs == {"GCA_000001.1", "GCA_000002.1", "GCA_000003.1", "GCA_000004.1"}


def test_tree2tax_outgroup_rooting(workdir: Path) -> None:
    _setup(workdir, with_outgroup=True)
    ctx = WorkdirContext(workdir, create=True)
    t2t, _ = tree2tax_run(ctx, Tree2taxParams(remove_outgroup=True))
    edges = _edges(t2t)
    # outgroup removed from the taxonomy after rooting
    assert not any(c == "Fam_Gen_sp_GCA_000004.1" for c, _ in edges)
    assert any(p == "root" for _, p in edges)


def test_tree2tax_include_dereplicated(workdir: Path) -> None:
    # a redundant genome listed under its representative leaf appears in genomes_map
    clusters = {"Fam_Gen_sp_GCA_000001.1": ["Fam_Gen_sp_GCA_000009.1"]}
    _setup(workdir, clusters=clusters)
    ctx = WorkdirContext(workdir, create=True)
    _, gmap = tree2tax_run(ctx, Tree2taxParams(include_dereplicated=True))
    accs = {ln.split("\t")[0] for ln in gmap.read_text().splitlines()}
    assert "GCA_000009.1" in accs  # the redundant member is mapped too


# An unrooted tree as mashtree/fasttree emit it: the seed node is trifurcating and
# the outgroup sits deep inside, next to one ingroup leaf.
_UNROOTED_NWK = (
    "(Fam_Gen_sp_GCA_000001.1:0.01,Fam_Gen_sp_GCA_000002.1:0.01,"
    "((Fam_Gen_sp_GCA_000003.1:0.02,Fam_Out_sp_GCA_000009.1:0.2):0.01,"
    "Fam_Gen_sp_GCA_000004.1:0.01):0.01);"
)


def test_tree2tax_roots_on_outgroup_edge_for_unrooted_tree(workdir: Path) -> None:
    (workdir / "tree").mkdir(parents=True)
    (workdir / "tree" / "tree.nwk").write_text(_UNROOTED_NWK + "\n")
    og = workdir / "outgroup"
    og.mkdir(parents=True)
    (og / "Fam_Out_sp_GCA_000009.1.fasta").write_text(">x\nACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009.1\n")
    ctx = WorkdirContext(workdir, create=True)
    t2t, _ = tree2tax_run(ctx, Tree2taxParams())
    edges = _edges(t2t)
    root_children = {c for c, p in edges if p == "root"}
    # The root splits the outgroup from one ingroup clade; no ingroup leaf may be
    # placed beside the outgroup because the input tree happened to be unrooted.
    assert "Fam_Out_sp_GCA_000009.1" in root_children
    assert len(root_children) == 2
    (ingroup,) = root_children - {"Fam_Out_sp_GCA_000009.1"}
    # The outgroup's former neighbour is inside the ingroup clade, not at the root.
    assert not ingroup.startswith("Fam_")
    assert ("Fam_Gen_sp_GCA_000003.1", "root") not in edges


def test_tree2tax_collapse_length_reparents_leaves_to_grandparent(workdir: Path) -> None:
    # (GCA_000001.1, GCA_000002.1) sit on a near-zero branch inside the ingroup.
    newick = (
        "(Fam_Out_sp_GCA_000009.1:0.3,"
        "((Fam_Gen_sp_GCA_000001.1:0.01,Fam_Gen_sp_GCA_000002.1:0.01):0.00001,"
        "Fam_Gen_sp_GCA_000003.1:0.02):0.3);"
    )
    (workdir / "tree").mkdir(parents=True)
    (workdir / "tree" / "tree.nwk").write_text(newick + "\n")
    og = workdir / "outgroup"
    og.mkdir(parents=True)
    (og / "Fam_Out_sp_GCA_000009.1.fasta").write_text(">x\nACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009.1\n")
    ctx = WorkdirContext(workdir, create=True)
    t2t, _ = tree2tax_run(ctx, Tree2taxParams(collapse_length=0.001))
    edges = _edges(t2t)
    parent = {c: p for c, p in edges}
    # All three ingroup leaves now share one parent: the ingroup clade node.
    assert (
        parent["Fam_Gen_sp_GCA_000001.1"]
        == parent["Fam_Gen_sp_GCA_000002.1"]
        == parent["Fam_Gen_sp_GCA_000003.1"]
    )
    assert parent[parent["Fam_Gen_sp_GCA_000003.1"]] == "root"
    params = ctx.config.stages["tree2tax"].params
    assert params["collapse_length"] == 0.001
    assert params["collapse_support"] is None
    assert params["collapsed_nodes"] == 1


def test_tree2tax_root_name_applies_without_remove_outgroup(workdir: Path) -> None:
    # --root-name labels the top node whether or not the outgroup is removed.
    _setup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    t2t, _ = tree2tax_run(ctx, Tree2taxParams(root_name="TOP"))
    parents = {p for _, p in _edges(t2t)}
    assert "TOP" in parents
    assert "root" not in parents


def test_tree2tax_lists_segment_members_under_their_isolate(workdir: Path) -> None:
    """A segment-grouped isolate genome maps its member accessions, not only its token."""
    from repgenr.core.contracts import SEGMENTS_TSV, write_segments

    _setup(workdir)
    write_segments(workdir / SEGMENTS_TSV, {"GCA_000001.1": ["SEG1.1", "SEG2.1"]})
    ctx = WorkdirContext(workdir, create=True)
    _, gmap = tree2tax_run(ctx, Tree2taxParams())
    rows = {ln.split("\t")[0]: ln.split("\t")[1] for ln in gmap.read_text().splitlines()}
    assert rows["SEG1.1"] == rows["SEG2.1"] == rows["GCA_000001.1"] == "Fam_Gen_sp_GCA_000001.1"
    assert "GCA_000002.1" in rows


def test_tree2tax_rejects_text_after_the_final_semicolon(workdir: Path) -> None:
    """A tree with text after its final ';' is treated as truncated, as doctor does."""
    import pytest

    from repgenr.core.errors import WorkdirError

    _setup(workdir)
    (workdir / "tree" / "tree.nwk").write_text(_NWK + "(extra\n")
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(WorkdirError, match=r"tree\.nwk.*truncated"):
        tree2tax_run(ctx, Tree2taxParams())


def _stray_outgroup(workdir: Path) -> None:
    og = workdir / "outgroup"
    og.mkdir(parents=True)
    (og / "Fam_Gen_sp_GCA_000099.1.fasta").write_text(">x\nACGT\n")
    (workdir / "outgroup_accession.txt").write_text("GCA_000099.1\n")


def test_tree2tax_outgroup_not_a_leaf_is_an_error(workdir: Path) -> None:
    """An outgroup that is not a leaf of the tree exits 3 and names the outgroup."""
    import pytest

    from repgenr.core.errors import WorkdirError

    _setup(workdir)
    _stray_outgroup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(WorkdirError, match=r"Fam_Gen_sp_GCA_000099\.1 is not a leaf"):
        tree2tax_run(ctx, Tree2taxParams())


def test_tree2tax_follows_phylo_built_without_outgroup(workdir: Path) -> None:
    """After `phylo --no-outgroup` the tree has no outgroup leaf; tree2tax leaves it unrooted."""
    _setup(workdir)
    _stray_outgroup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    ctx.config.record_stage(
        "phylo", tool="mashtree", params={"outgroup": None}, completed="2026-10-07T00:00:00"
    )
    t2t, _ = tree2tax_run(ctx, Tree2taxParams())
    assert any(p == "root" for _, p in _edges(t2t))


def test_tree2tax_records_dendropy_as_its_tool(workdir: Path) -> None:
    """The record names the library that does the work and its version."""
    from importlib.metadata import version

    _setup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    tree2tax_run(ctx, Tree2taxParams())
    record = ctx.config.stages["tree2tax"]
    assert record.tool == "dendropy"
    assert record.tool_versions == {"dendropy": version("dendropy")}


def test_tree2tax_rejects_two_concatenated_trees(workdir: Path) -> None:
    """A second complete tree after the first is refused, not silently dropped."""
    import pytest

    from repgenr.core.errors import WorkdirError

    _setup(workdir)
    (workdir / "tree" / "tree.nwk").write_text(_NWK.strip() + "\n" + _NWK)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(WorkdirError, match=r"more than one tree"):
        tree2tax_run(ctx, Tree2taxParams())


def test_tree2tax_warns_when_the_last_phylo_run_did_not_finish(workdir: Path, caplog) -> None:
    """An interrupted phylo keeps the previous tree.nwk; tree2tax says so."""
    import logging

    _setup(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=logging.getLogger("test-t2t"))
    ctx.config.record_stage("phylo", tool="fasttree", completed=None)
    with caplog.at_level(logging.WARNING, logger="test-t2t"):
        tree2tax_run(ctx, Tree2taxParams())
    assert any("did not finish" in r.getMessage() for r in caplog.records)

    caplog.clear()
    ctx.config.record_stage("phylo", tool="fasttree", completed="2026-10-07T00:00:00")
    with caplog.at_level(logging.WARNING, logger="test-t2t"):
        tree2tax_run(ctx, Tree2taxParams())
    assert not any("did not finish" in r.getMessage() for r in caplog.records)


def test_unmatched_outgroup_accession_names_the_outgroup_directory(tmp_path: Path, caplog) -> None:
    """An accession no outgroup file matches is reported as such, not as a
    missing leaf; the tree is left unrooted as phylo-build leaves it."""
    import logging

    from repgenr.stages.tree2tax import Tree2taxStepParams, tree2tax_relations

    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    og = tmp_path / "og"
    og.mkdir()
    (og / "Fam_Gen_sp_GCA_000004.1.fasta").write_text(">x\nACGT\n")
    acc = tmp_path / "acc.txt"
    acc.write_text("GCA_999999.1\n")
    with caplog.at_level(logging.WARNING, logger="test-t2t"):
        tree2tax_relations(
            Tree2taxStepParams(
                tree=tree, out_dir=tmp_path / "out", outgroup_dir=og, outgroup_accession=acc
            ),
            logging.getLogger("test-t2t"),
        )
    messages = [r.getMessage() for r in caplog.records]
    assert any("No file in" in m and "GCA_999999.1" in m for m in messages), messages
    assert not any("among tree leaves" in m for m in messages)


def test_tree2tax_writes_the_node_table_with_lineage(workdir: Path) -> None:
    from repgenr.core.contracts import TREE2TAX_NODES_TSV
    from repgenr.core.manifest import GenomeRecord

    _setup(workdir)
    ctx = WorkdirContext(workdir, create=True)
    ctx.manifest.upsert(
        GenomeRecord(accession="GCA_000001.1", family="Fam", genus="Gen", species="sp")
    )
    tree2tax_run(ctx, Tree2taxParams())
    lines = (workdir / TREE2TAX_NODES_TSV).read_text().splitlines()
    rows = {ln.split("\t")[0]: ln.split("\t") for ln in lines[1:]}
    assert "root" in rows and rows["Fam_Gen_sp_GCA_000001.1"][-1] == "Fam;Gen;sp"
    assert rows["Fam_Gen_sp_GCA_000002.1"][-1] == ""
