"""graft: join the hierarchies of two working directories at the node they share.

A group too divergent for one alignment is built in tiers: a parent working
directory over the whole set, in which one clade holds a subset, and a child
working directory of that subset. graft replaces the clade with the child's
tree. The seam is the parent node whose genomes are exactly the child's.

The parent here: root -> Paeno (its outgroup), N1; N1 -> Ochro, N2;
N2 -> Tele, C; C -> c1, c2, c3. The child is rooted on Tele, its own outgroup,
which tree2tax keeps as a leaf: root -> Tele, X; X -> Y, c3; Y -> c1, c2.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import write_genomes_map, write_tree2tax
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.stages.graft import GraftStepParams, graft

_LOG = logging.getLogger("test")

PARENT_EDGES = [
    ("Paeno", "root"),
    ("N1", "root"),
    ("Ochro", "N1"),
    ("N2", "N1"),
    ("Tele", "N2"),
    ("C", "N2"),
    ("c1", "C"),
    ("c2", "C"),
    ("c3", "C"),
]
PARENT_MAP = [
    ("GCF_pae", "Paeno"),
    ("GCF_och", "Ochro"),
    ("GCF_tel", "Tele"),
    ("GCF_c1", "c1"),
    ("GCF_c2", "c2"),
    ("GCF_c3", "c3"),
]
CHILD_EDGES = [("Tele", "root"), ("X", "root"), ("Y", "X"), ("c3", "X"), ("c1", "Y"), ("c2", "Y")]
CHILD_MAP = [("GCF_tel", "Tele"), ("GCF_c1", "c1"), ("GCF_c2", "c2"), ("GCF_c3", "c3")]


def _workdir(path: Path, edges, mapping, outgroup: str | None) -> Path:
    path.mkdir(parents=True)
    write_tree2tax(path / "tree2tax.tsv", list(edges))
    write_genomes_map(path / "genomes_map.tsv", list(mapping))
    if outgroup is not None:
        (path / "outgroup_accession.txt").write_text(outgroup + "\n")
    return path


def _edges(path: Path) -> set[tuple[str, str]]:
    rows = path.read_text(encoding="utf-8").splitlines()[1:]
    return {tuple(row.split("\t")[:2]) for row in rows}  # type: ignore[misc]


def _map(path: Path) -> set[tuple[str, str]]:
    return {tuple(row.split("\t")) for row in path.read_text(encoding="utf-8").splitlines()}  # type: ignore[misc]


def _graft(
    tmp_path: Path,
    parent_edges=PARENT_EDGES,
    parent_map=PARENT_MAP,
    child_edges=CHILD_EDGES,
    child_map=CHILD_MAP,
) -> tuple[Path, Path]:
    parent = _workdir(tmp_path / "parent", parent_edges, parent_map, "GCF_pae")
    child = _workdir(tmp_path / "child", child_edges, child_map, "GCF_tel")
    return graft(GraftStepParams(parent=parent, children=[child], out_dir=tmp_path / "out"), _LOG)


def test_the_child_tree_replaces_the_clade_at_the_seam(tmp_path: Path) -> None:
    tree, _genomes = _graft(tmp_path)
    assert _edges(tree) == {
        ("Paeno", "root"),
        ("N1", "root"),
        ("Ochro", "N1"),
        ("N2", "N1"),
        ("Tele", "N2"),
        ("C", "N2"),
        # the child below the seam; X held exactly the seam's genomes once Tele
        # was gone, so its children hang from the seam directly
        ("Y", "C"),
        ("c3", "C"),
        ("c1", "Y"),
        ("c2", "Y"),
    }


def test_the_childs_outgroup_is_not_grafted_a_second_time(tmp_path: Path) -> None:
    # Tele roots the child, and tree2tax keeps it as a leaf; in the parent it sits
    # outside the seam, where it stays.
    tree, _genomes = _graft(tmp_path)
    assert [edge for edge in _edges(tree) if edge[0] == "Tele"] == [("Tele", "N2")]


def test_the_genomes_map_joins_both(tmp_path: Path) -> None:
    _tree, genomes = _graft(tmp_path)
    assert _map(genomes) == {
        ("GCF_pae", "Paeno"),
        ("GCF_och", "Ochro"),
        ("GCF_tel", "Tele"),
        ("GCF_c1", "c1"),
        ("GCF_c2", "c2"),
        ("GCF_c3", "c3"),
    }


def test_genomes_spread_over_two_branches_are_refused(tmp_path: Path) -> None:
    # c3 sits outside C in the parent, so no parent node holds exactly the
    # child's genomes: the two trees disagree about the group, which is a
    # finding to report rather than something to resolve by taking an ancestor.
    edges = [e for e in PARENT_EDGES if e != ("c3", "C")] + [("c3", "N1")]
    with pytest.raises(UserInputError, match="no node of the parent") as raised:
        _graft(tmp_path, parent_edges=edges)
    assert "N1" in str(raised.value), "it names the smallest node holding them all"


def test_genomes_the_parent_lacks_are_refused(tmp_path: Path) -> None:
    child_map = [*CHILD_MAP, ("GCF_c4", "c4")]
    child_edges = [*CHILD_EDGES, ("c4", "Y")]
    with pytest.raises(UserInputError, match="GCF_c4"):
        _graft(tmp_path, child_edges=child_edges, child_map=child_map)


def test_node_names_both_trees_use_are_refused(tmp_path: Path) -> None:
    # Two trees named with one --node-basename give BR1, BR2 ... in both. The
    # child's N2 survives the graft (X is folded into the root), and the parent
    # keeps an N2 of its own, above the seam.
    child_edges = [
        ("Tele", "root"),
        ("X", "root"),
        ("N2", "X"),
        ("c3", "X"),
        ("c1", "N2"),
        ("c2", "N2"),
    ]
    with pytest.raises(UserInputError, match="N2") as raised:
        _graft(tmp_path, child_edges=child_edges)
    assert "--node-basename" in str(raised.value), "it says how to avoid it"


def test_a_workdir_without_a_hierarchy_is_named(tmp_path: Path) -> None:
    parent = _workdir(tmp_path / "parent", PARENT_EDGES, PARENT_MAP, "GCF_pae")
    empty = tmp_path / "child"
    empty.mkdir()
    with pytest.raises(WorkdirError, match="tree2tax"):
        graft(GraftStepParams(parent=parent, children=[empty], out_dir=tmp_path / "out"), _LOG)


def test_children_are_grafted_in_turn(tmp_path: Path) -> None:
    # A flat parent, A -> a1 a2 a3 and B -> b1 b2 b3, and a child for each that
    # resolves its clade further.
    parent = _workdir(
        tmp_path / "parent",
        [
            ("A", "root"),
            ("B", "root"),
            *[(f"a{i}", "A") for i in (1, 2, 3)],
            *[(f"b{i}", "B") for i in (1, 2, 3)],
        ],
        [(f"GCF_{x}{i}", f"{x}{i}") for x in "ab" for i in (1, 2, 3)],
        None,
    )
    children = [
        _workdir(
            tmp_path / x,
            [
                (f"{x.upper()}{x.upper()}", "root"),
                (f"{x}3", "root"),
                (f"{x}1", f"{x.upper()}{x.upper()}"),
                (f"{x}2", f"{x.upper()}{x.upper()}"),
            ],
            [(f"GCF_{x}{i}", f"{x}{i}") for i in (1, 2, 3)],
            None,
        )
        for x in "ab"
    ]
    tree, _genomes = graft(
        GraftStepParams(parent=parent, children=children, out_dir=tmp_path / "out"), _LOG
    )
    assert _edges(tree) == {
        ("A", "root"),
        ("B", "root"),
        ("AA", "A"),
        ("a3", "A"),
        ("a1", "AA"),
        ("a2", "AA"),
        ("BB", "B"),
        ("b3", "B"),
        ("b1", "BB"),
        ("b2", "BB"),
    }
