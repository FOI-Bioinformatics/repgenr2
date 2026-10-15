"""Neighbor-joining + Newick helpers (used by the sourmash tree builder)."""

from __future__ import annotations

import re

from repgenr.tree.newick import neighbor_joining


def _leaves(newick: str) -> set[str]:
    # leaf labels are tokens before a ':' that aren't preceded by ')'
    return set(re.findall(r"[(,]([A-Za-z0-9_]+):", newick))


def test_single_leaf() -> None:
    assert neighbor_joining(["A"], [[0.0]]) == "(A:0.0);"


def test_two_leaves_split_distance() -> None:
    assert neighbor_joining(["A", "B"], [[0.0, 4.0], [4.0, 0.0]]) == "(A:2.000000,B:2.000000);"


def test_four_taxa_recovers_expected_grouping() -> None:
    # A,B close; C,D close; the two pairs are far apart.
    labels = ["A", "B", "C", "D"]
    dist = [
        [0.0, 0.1, 1.0, 1.0],
        [0.1, 0.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, 0.1],
        [1.0, 1.0, 0.1, 0.0],
    ]
    nwk = neighbor_joining(labels, dist)
    assert nwk.endswith(");")
    assert _leaves(nwk) == {"A", "B", "C", "D"}
    # A and B should be adjacent in the tree string (a cherry), likewise C and D
    assert "A_" in nwk.replace("A:", "A_") and "B" in nwk
    assert re.search(r"\(A:[0-9.]+,B:[0-9.]+\)", nwk) or re.search(r"\(B:[0-9.]+,A:[0-9.]+\)", nwk)
    assert re.search(r"\(C:[0-9.]+,D:[0-9.]+\)", nwk) or re.search(r"\(D:[0-9.]+,C:[0-9.]+\)", nwk)


def test_labels_with_reserved_characters_are_quoted() -> None:
    import dendropy

    labels = ["a b", "c(d)", "a_b", "it's"]
    dist = [[0.0, 2.0, 3.0, 3.0], [2.0, 0.0, 3.0, 3.0], [3.0, 3.0, 0.0, 2.0], [3.0, 3.0, 2.0, 0.0]]
    nwk = neighbor_joining(labels, dist)
    assert "'a b'" in nwk and "'c(d)'" in nwk and "'it''s'" in nwk and "a_b" in nwk
    tree = dendropy.Tree.get(data=nwk, schema="newick", preserve_underscores=True)
    assert sorted(t.label for t in tree.taxon_namespace) == sorted(labels)
