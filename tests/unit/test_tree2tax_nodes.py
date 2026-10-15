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
    assert {r["support"] for r in inner} == {"90.0000", "80.0000"}
    assert rows["A"]["red"] == "1.0000" and rows["A"]["depth"] == "2.000000"
    assert rows["A"]["branch_length"] == "1.000000" and rows["A"]["support"] == "NA"
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
