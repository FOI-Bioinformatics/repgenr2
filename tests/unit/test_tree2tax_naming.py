"""Node naming and edge emission in tree2tax: one pass, stable names."""

from __future__ import annotations

import logging
from pathlib import Path

import dendropy

from repgenr.stages.tree2tax import _edges, _emit_relations, _name_nodes, _prune_leaf

_LOG = logging.getLogger("test-naming")


def _tree(text: str) -> dendropy.Tree:
    return dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True)


def _names(text: str) -> dict[frozenset[str], str]:
    """Internal node name keyed by its leaf set."""
    tree = _tree(text)
    _name_nodes(tree, None)
    return {
        frozenset(lf.taxon.label for lf in node.leaf_iter()): node.label
        for node in tree.internal_nodes()
        if node is not tree.seed_node
    }


def test_name_depends_on_the_leaf_set_only() -> None:
    a = _names("(((A,B),C),(D,E));")
    b = _names("((C,(B,A)),(E,D));")  # same clades, children swapped
    assert a == b
    assert len(set(a.values())) == len(a)  # distinct clades, distinct names


def test_name_is_independent_of_the_resolution_inside_the_clade() -> None:
    resolved = _names("(((A,B),C),(D,E));")
    polytomy = _names("((A,B,C),(D,E));")
    key = frozenset("ABC")
    assert resolved[key] == polytomy[key]


def test_a_unary_node_is_not_named_like_its_child() -> None:
    tree = _tree("(((A,B)),(C,D));")  # ((A,B)) is a node with one child
    _name_nodes(tree, None)
    labels = [n.label for n in tree.internal_nodes() if n is not tree.seed_node]
    assert len(labels) == len(set(labels)) == 3


def test_a_unary_node_does_not_rename_its_ancestors() -> None:
    plain = _names("(((A,B),E),(C,D));")
    with_unary = _names("((((A,B)),E),(C,D));")
    for key in (frozenset("AB"), frozenset("ABE"), frozenset("CD")):
        assert plain[key] == with_unary[key]


def test_an_unlabelled_leaf_does_not_make_a_node_share_its_childs_name() -> None:
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree(taxon_namespace=taxa)
    top = tree.seed_node.new_child()
    inner = top.new_child()
    for label in ("A", "B"):
        inner.new_child().taxon = taxa.require_taxon(label=label)
    top.new_child()  # a leaf without a label contributes nothing to the sum
    other = tree.seed_node.new_child()
    for label in ("C", "D"):
        other.new_child().taxon = taxa.require_taxon(label=label)
    _name_nodes(tree, None)
    assert top.label != inner.label
    tree.seed_node.label = "root"
    assert all(child != parent for child, parent in _edges(tree))


def test_basename_numbers_internal_nodes_in_preorder() -> None:
    tree = _tree("((A,B),(C,D));")
    leaves = _name_nodes(tree, "N")
    assert sorted(leaves) == ["A", "B", "C", "D"]
    assert [n.label for n in tree.preorder_internal_node_iter() if n is not tree.seed_node] == [
        "N0",
        "N1",
    ]


def test_edges_once_per_node_in_a_ladder() -> None:
    # Built node by node: the Newick parser recurses once per nesting level.
    n = 5000
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree(taxon_namespace=taxa)
    node = tree.seed_node
    for i in range(n - 1):
        leaf = node.new_child()
        leaf.taxon = taxa.require_taxon(label=f"A{i}")
        node = node.new_child()
    node.taxon = taxa.require_taxon(label=f"A{n - 1}")
    _name_nodes(tree, None)
    tree.seed_node.label = "root"
    edges = _edges(tree)
    assert len(edges) == len(set(edges)) == 2 * n - 2  # every non-root node once


def test_prune_leaf_removes_the_unifurcation_it_leaves() -> None:
    # The outgroup sits inside a clade, as on a midpoint-rooted tree.
    tree = _tree("((A:1,B:1):1,(C:1,(D:1,OG:1):2):1);")
    _prune_leaf(tree, "OG")
    assert sorted(lf.taxon.label for lf in tree.leaf_node_iter()) == ["A", "B", "C", "D"]
    assert all(n.num_child_nodes() >= 2 for n in tree.internal_nodes())
    d = next(lf for lf in tree.leaf_node_iter() if lf.taxon.label == "D")
    assert d.edge.length == 3  # the collapsed node's length is added


def test_prune_leaf_keeps_the_root_over_one_ingroup_clade() -> None:
    tree = _tree("(OG:1,((A:1,B:1):1,C:1):1);")  # rooted on the outgroup
    _prune_leaf(tree, "OG")
    assert tree.seed_node.num_child_nodes() == 1


def test_remove_outgroup_on_a_nested_outgroup_leaves_no_single_child_node(
    tmp_path: Path,
) -> None:
    text = "((A:1,B:1):1,(C:1,(D:1,OG:1):2):1);"
    t2t, _, _ = _emit_relations(
        text,
        "OG",
        {},
        node_basename=None,
        root_name="root",
        remove_outgroup=True,
        out_tree2tax=tmp_path / "t.tsv",
        out_map=tmp_path / "m.tsv",
        logger=_LOG,
        keep_root=True,
    )
    rows = [line.split("\t") for line in t2t.read_text().splitlines()[1:]]
    children_of: dict[str, set[str]] = {}
    for child, parent in rows:
        children_of.setdefault(parent, set()).add(child)
    assert "OG" not in {c for c, _ in rows}
    assert all(len(kids) >= 2 for kids in children_of.values())


def test_a_clade_keeps_its_name_when_an_unrelated_clade_changes() -> None:
    before = _names("(((A,B),C),(D,E));")
    after = _names("(((A,B),C),((D,E),F));")
    assert before[frozenset("AB")] == after[frozenset("AB")]
    assert before[frozenset("ABC")] == after[frozenset("ABC")]


def test_remove_outgroup_keeps_the_ingroup_clade_under_root(tmp_path: Path) -> None:
    text = "(OG:1,((A:1,B:1):1,(C:1,D:1):1):1);"  # rooted on the outgroup
    t2t, _, _ = _emit_relations(
        text,
        "OG",
        {},
        node_basename=None,
        root_name="root",
        remove_outgroup=True,
        out_tree2tax=tmp_path / "t.tsv",
        out_map=tmp_path / "m.tsv",
        logger=_LOG,
    )
    rows = [tuple(line.split("\t")) for line in t2t.read_text().splitlines()[1:]]
    parent_of = dict(rows)
    under_root = [c for c, p in rows if p == "root"]
    assert len(under_root) == 1  # the ingroup clade, not its two halves
    leaves = {c for c in parent_of if c not in parent_of.values()}

    def ancestors(leaf: str) -> set[str]:
        out, node = set(), leaf
        while node in parent_of:
            node = parent_of[node]
            out.add(node)
        return out

    assert leaves == {"A", "B", "C", "D"}
    assert all(under_root[0] in ancestors(leaf) for leaf in leaves)


def test_stacked_unary_nodes_get_distinct_names_and_leave_ancestors_alone() -> None:
    plain = _names("(((A,B),E),(C,D));")
    for text in ("(((((A,B))),E),(C,D));", "((((((A,B)))),E),(C,D));", "((((A)),B),(C,D));"):
        tree = _tree(text)
        _name_nodes(tree, None)
        labels = [n.label for n in tree.internal_nodes() if n is not tree.seed_node]
        assert len(labels) == len(set(labels)), text
        tree.seed_node.label = "root"
        assert all(child != parent for child, parent in _edges(tree)), text
    with_chain = _names("((((((A,B)))),E),(C,D));")
    assert with_chain[frozenset("ABE")] == plain[frozenset("ABE")]


def test_a_node_over_unlabelled_leaves_only_gets_no_name_and_no_row() -> None:
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree(taxon_namespace=taxa)
    empty = tree.seed_node.new_child()
    empty.new_child()
    empty.new_child()
    other = tree.seed_node.new_child()
    for label in ("C", "D"):
        other.new_child().taxon = taxa.require_taxon(label=label)
    assert _name_nodes(tree, None) == ["C", "D"]
    assert empty.label is None
    tree.seed_node.label = "root"
    assert {c for c, _ in _edges(tree)} == {"C", "D", other.label}
