"""The node annotation table beside tree2tax.tsv: sizes, depths, supports, RED."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.core.contracts import NODE_TABLE_COLUMNS
from repgenr.stages.tree2tax import _emit_relations

_LOG = logging.getLogger("test-nodes")


def _emit(tmp_path: Path, text: str, **kwargs) -> dict[str, dict[str, str]]:
    kwargs.setdefault("outgroup_leaf", None)
    _emit_relations(
        text,
        kwargs.pop("outgroup_leaf"),
        {},
        node_basename="N",
        root_name="root",
        remove_outgroup=False,
        out_tree2tax=tmp_path / "t.tsv",
        out_map=tmp_path / "m.tsv",
        logger=_LOG,
        out_nodes=tmp_path / "n.tsv",
        **kwargs,
    )
    lines = (tmp_path / "n.tsv").read_text().splitlines()
    assert lines[0].split("\t") == list(NODE_TABLE_COLUMNS)
    return {
        row[0]: dict(zip(NODE_TABLE_COLUMNS, row, strict=True))
        for row in (ln.split("\t") for ln in lines[1:])
    }


def test_red_on_a_symmetric_tree_matches_the_recursion(tmp_path: Path) -> None:
    # Root at the midpoint of a symmetric tree: both inner nodes are 1 from the
    # root and 1 from their leaves, so u = 2 at the root, red = 0 + 1/2 = 0.5.
    rows = _emit(tmp_path, "((A:1,B:1)90:1,(C:1,D:1)80:1);", keep_root=True)
    assert rows["root"]["red"] == "0.0000" and rows["root"]["n_leaves"] == "4"
    inner = [r for r in rows.values() if r["node"].startswith("N")]
    assert {r["red"] for r in inner} == {"0.5000"}
    assert {r["support"] for r in inner} == {"0.9000", "0.8000"}
    assert rows["A"]["red"] == "1.0000" and rows["A"]["depth"] == "2"
    assert rows["A"]["branch_length"] == "1" and rows["A"]["support"] == "NA"
    assert rows["A"]["parent"] in {r["node"] for r in inner}


def test_red_and_depth_are_na_without_lengths_or_a_root(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "((A,B),(C,D));", keep_root=True)
    assert rows["A"]["red"] == "NA" and rows["A"]["depth"] == "NA"
    assert rows["root"]["n_leaves"] == "4"
    rows = _emit(tmp_path, "((A:1,B:1):1,(C:1,D:1):1);")  # unrooted: no outgroup, no kept root
    assert rows["A"]["red"] == "NA" and rows["A"]["branch_length"] == "NA"


def test_outgroup_rooting_makes_red_measurable(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "((A:1,B:1):1,(C:1,OG:1):1);", outgroup_leaf="OG")
    assert rows["OG"]["red"] == "1.0000" and rows["root"]["red"] == "0.0000"
    assert all(r["red"] != "NA" for r in rows.values())
    # Rooted on OG's branch (0.5 each side): the ingroup node (C,(A,B)) has
    # a = 0.5 and b = mean(1, 1 + 1 + 1, 1 + 1 + 1) over C, A, B via the AB
    # node's edge of 1, i.e. (1 + 3 + 3) / 3 = 7/3, red = 0.5/(0.5 + 7/3) =
    # 0.1765; the AB node has a = 1, b = 1, red = 0.1765 + 0.5 * 0.8235 =
    # 0.5882... as PhyloRank's recursion gives 0.7255 with b measured over
    # the node's own leaves (1): red = 0.1765 + (1 / (1 + 1)) * 0.8235.
    inner = sorted(r["red"] for r in rows.values() if r["node"].startswith("N"))
    assert inner == ["0.1765", "0.7255"]


def test_node_table_names_every_node_of_the_relations(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "(((A:1,B:1)95:1,C:1)70:1,(D:1,E:1)50:1);", keep_root=True)
    edges = [ln.split("\t") for ln in (tmp_path / "t.tsv").read_text().splitlines()[1:]]
    nodes = {c for c, _ in edges} | {p for _, p in edges}
    assert set(rows) == nodes
    for child, parent in edges:
        assert rows[child]["parent"] == parent
    assert rows["root"]["parent"] == ""


def test_lineage_is_written_for_leaves_that_have_one(tmp_path: Path) -> None:
    rows = _emit(
        tmp_path,
        "((Fam_Gen_sp_GCA_1.1:1,B:1):1,(C:1,D:1):1);",
        keep_root=True,
        lineage={"GCA_1.1": "Fam;Gen;sp"},
    )
    assert rows["Fam_Gen_sp_GCA_1.1"]["lineage"] == "Fam;Gen;sp" and rows["B"]["lineage"] == ""


def test_unbalanced_sister_clades_get_different_red(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "((A:0.1,B:0.1)90:0.1,(C:5,D:5)80:0.1);", keep_root=True)
    by_support = {r["support"]: r["red"] for r in rows.values() if r["node"].startswith("N")}
    assert by_support["0.9000"] == "0.5000"  # a = 0.1, b = 0.1
    assert by_support["0.8000"] == "0.0196"  # a = 0.1, b = 5


def test_zero_length_clade_sits_at_its_parents_red(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "(((A:0,B:0)70:0,C:1)90:1,(D:1,E:1)80:1);", keep_root=True)
    reds = {r["support"]: r["red"] for r in rows.values() if r["node"].startswith("N")}
    assert reds["0.7000"] == reds["0.9000"]


def test_remove_outgroup_measures_red_from_the_ingroup_ancestor(tmp_path: Path) -> None:
    _emit_relations(
        "(((A:1,B:1)90:1,C:1)95:1,OG:5);",
        "OG",
        {},
        node_basename="N",
        root_name="root",
        remove_outgroup=True,
        out_tree2tax=tmp_path / "t.tsv",
        out_map=tmp_path / "m.tsv",
        logger=_LOG,
        out_nodes=tmp_path / "n.tsv",
    )
    lines = (tmp_path / "n.tsv").read_text().splitlines()
    rows = {
        row[0]: dict(zip(NODE_TABLE_COLUMNS, row, strict=True))
        for row in (ln.split("\t") for ln in lines[1:])
    }
    ingroup = next(r for r in rows.values() if r["parent"] == "root")
    assert ingroup["red"] == "0.0000"
    ab = next(r for r in rows.values() if r["support"] == "0.9000")
    assert ab["red"] == "0.5000"  # a = 1, b = 1 below the ingroup ancestor


def test_percent_supports_are_written_as_fractions(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "((A:1,B:1)95:1,(C:1,D:1)0.5:1);", keep_root=True)
    assert {r["support"] for r in rows.values() if r["node"].startswith("N")} == {
        "0.9500",
        "0.0050",
    }


def test_the_stateless_step_writes_the_node_table(tmp_path: Path) -> None:
    from repgenr.core.contracts import TREE2TAX_NODES_TSV
    from repgenr.stages.tree2tax import Tree2taxStepParams, tree2tax_relations

    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    tree2tax_relations(
        Tree2taxStepParams(tree=tree, out_dir=tmp_path / "out", no_outgroup=True), _LOG
    )
    lines = (tmp_path / "out" / TREE2TAX_NODES_TSV).read_text().splitlines()
    assert lines[0].split("\t") == list(NODE_TABLE_COLUMNS) and len(lines) == 1 + 7
