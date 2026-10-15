"""The optional rank column of tree2tax.tsv (--ranks, --rank-red)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import write_tree2tax
from repgenr.core.errors import UserInputError
from repgenr.stages import tree2tax as mod
from repgenr.stages.tree2tax import _emit_relations, _parse_rank_red

_LOG = logging.getLogger("test-ranks")


def _emit(tmp_path: Path, newick: str, **kwargs) -> dict[str, tuple[str, str]]:
    """child -> (parent, rank) of the written table."""
    kwargs.setdefault("node_basename", "N")
    kwargs.setdefault("root_name", "root")
    kwargs.setdefault("remove_outgroup", False)
    kwargs.setdefault("keep_root", True)
    t2t, _, _ = _emit_relations(
        newick,
        None,
        {},
        out_tree2tax=tmp_path / "t.tsv",
        out_map=tmp_path / "m.tsv",
        logger=_LOG,
        ranks=True,
        **kwargs,
    )
    lines = t2t.read_text().splitlines()
    assert lines[0].split("\t") == ["child", "parent", "rank"]
    return {c: (p, r) for c, p, r in (ln.split("\t") for ln in lines[1:])}


def test_parse_rank_red_sorts_by_threshold() -> None:
    assert _parse_rank_red("species=0.95, genus=0.8") == [("genus", 0.8), ("species", 0.95)]
    assert _parse_rank_red(None) == [] and _parse_rank_red(" ") == []


@pytest.mark.parametrize(
    "text",
    [
        "species",
        "species=",
        "=0.9",
        "species=x",
        "species=0",
        "species=1.5",
        "species=0.9,species=0.8",
        "strain=0.9",
        "species=0.9,genus=0.9",
    ],
)
def test_parse_rank_red_refuses_malformed_bands(text: str) -> None:
    with pytest.raises(UserInputError, match="--rank-red"):
        _parse_rank_red(text)


def test_ranks_without_bands_are_strain_and_no_rank(tmp_path: Path) -> None:
    rows = _emit(tmp_path, "((A,B),(C,D));")
    assert {r for _, r in rows.values()} == {mod.RANK_LEAF, mod.RANK_NONE}
    assert rows["A"][1] == "strain" and rows["N0"][1] == "no rank"
    assert "root" not in rows  # the root has no row; FlexTaxD pre-creates it


def test_rank_red_needs_ranks(tmp_path: Path) -> None:
    with pytest.raises(UserInputError, match="needs --ranks"):
        _emit_relations(
            "((A,B),(C,D));",
            None,
            {},
            node_basename="N",
            root_name="root",
            remove_outgroup=False,
            out_tree2tax=tmp_path / "t.tsv",
            out_map=tmp_path / "m.tsv",
            logger=_LOG,
            rank_red="species=0.9",
        )


def test_rank_red_needs_red(tmp_path: Path) -> None:
    with pytest.raises(UserInputError, match="relative evolutionary divergence"):
        _emit(tmp_path, "((A,B),(C,D));", rank_red="species=0.9")  # no branch lengths
    with pytest.raises(UserInputError, match="relative evolutionary divergence"):
        _emit(tmp_path, "((A:1,B:1):1,(C:1,D:1):1);", rank_red="species=0.9", keep_root=False)


def test_bands_assign_each_rank_once_per_lineage(tmp_path: Path) -> None:
    # RED after Parks et al.: N0 = ((A:1,B:1):1,C:2):2 has a = 2 and a mean
    # of 2 to its leaves, so 0.5; N1 = (A:1,B:1):1 below it has a = b = 1,
    # so 0.5 + 0.5 * 0.5 = 0.75; the sister N2 = (D:1,E:1):2 has 2/(2+1) = 0.667.
    newick = "(((A:1,B:1):1,C:2):2,(D:1,E:1):2);"
    rows = _emit(tmp_path, newick, rank_red="species=0.7,genus=0.5")
    ranks = {c: r for c, (_, r) in rows.items()}
    by_parent = {c: p for c, (p, _) in rows.items()}
    n1 = by_parent["A"]
    n0 = by_parent[n1]
    n2 = by_parent["D"]
    assert ranks[n0] == "genus" and ranks[n1] == "species" and ranks[n2] == "genus"
    assert ranks["A"] == "strain" and ranks["C"] == "strain"
    # A rank is taken once per lineage: N1 at 0.667 under a genus N0 at 0.333
    # stays 'no rank' when only a genus band is given.
    rows = _emit(tmp_path, "(((A:1,B:1):1,C:2):1,(D:1,E:1):1);", rank_red="genus=0.3")
    ranks = {c: r for c, (_, r) in rows.items()}
    inner = [c for c, r in ranks.items() if r != "strain"]
    assert sorted(ranks[c] for c in inner) == ["genus", "genus", "no rank"]


def test_write_tree2tax_refuses_a_missing_rank(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="has no rank"):
        write_tree2tax(tmp_path / "t.tsv", [("A", "root")], {"B": "strain"})
    write_tree2tax(tmp_path / "t.tsv", [("A", "root")], None)
    assert (tmp_path / "t.tsv").read_text() == "child\tparent\nA\troot\n"


def test_a_unary_internal_node_is_not_a_strain(tmp_path: Path) -> None:
    # N1 = (A:1):1 has one leaf below it but is internal: RED 0.5 + 0.5 * 1/2 = 0.75.
    rows = _emit(tmp_path, "(((A:1):1,D:1):1,(B:1,C:1):1);", rank_red="species=0.7")
    ranks = {c: r for c, (_, r) in rows.items()}
    by_parent = {c: p for c, (p, _) in rows.items()}
    unary = by_parent["A"]
    assert unary.startswith("N") and ranks[unary] == "species"
    assert ranks["A"] == "strain"


def test_rank_red_text_is_normalised() -> None:
    assert mod.rank_red_text("species=0.95, genus=0.8") == "genus=0.8,species=0.95"
    assert mod.rank_red_text(None) is None and mod.rank_red_text("") is None
