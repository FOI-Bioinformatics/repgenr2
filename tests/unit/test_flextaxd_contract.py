"""The tree2tax tables satisfy FlexTaxD's reader.

The checks mirror ``flextaxd/modules/ReadTaxonomy.py`` on FlexTaxD's GitHub
master (read 2026-10): the relation file needs a header holding ``parent``
and ``child`` in either order; nodes are keyed by their label, so an empty
label or a label under two parents corrupts the tree; a third column is read
as a rank and an empty one silently inherits the previous row's, so the file
carries two columns only; the top node must be ``root``, which FlexTaxD
pre-creates and validates from; the genome map has no header and two columns
whose second names a node of the relation file. The map check is stricter
than FlexTaxD, which accepts a third column (a reference flag) and only warns
about an unknown node. FlexTaxD itself is not a test dependency; a live test
imports the files where the binary is present.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.stages.tree2tax import FLEXTAXD_ROOT, _emit_relations

_LOG = logging.getLogger("test-flextaxd")

_BINARY = (
    "((Fam_Gen_sp_GCA_000001.1:0.1,Fam_Gen_sp_GCA_000002.1:0.1)0.9:0.2,"
    "(Fam_Gen_sp_GCA_000003.1:0.1,Fam_Gen_sp_GCA_000004.1:0.1)0.4:0.2,"
    "Out_Gen_sp_GCA_000009.1:0.5);"
)


def _mirror_read(t2t: Path, gmap: Path) -> tuple[dict[str, str], list[tuple[str, str]]]:
    """Read the two files the way FlexTaxD's master reader does, raising on
    what it would store wrongly. Returns child -> parent and the map rows."""
    lines = t2t.read_text(encoding="utf-8").splitlines()
    headers = lines[0].strip().split("\t")
    assert "parent" in headers and "child" in headers, "header must name parent and child"
    swap = headers[0] == "parent"
    parent_of: dict[str, str] = {}
    for row in lines[1:]:
        data = row.strip().split("\t")
        if data[0] == "":
            continue  # FlexTaxD skips an empty row
        # A row with an empty field loses that column when FlexTaxD strips
        # it, so two fields after the strip also means no empty label.
        assert len(data) == 2, f"a third column would be read as a rank: {row!r}"
        child, parent = (data[1], data[0]) if swap else (data[0], data[1])
        assert parent_of.setdefault(child, parent) == parent, f"{child!r} under two parents"
    nodes = set(parent_of) | set(parent_of.values())
    tops = nodes - set(parent_of)
    assert tops == {FLEXTAXD_ROOT}, f"FlexTaxD validates from 'root'; top nodes are {tops}"
    rows = []
    for row in gmap.read_text(encoding="utf-8").splitlines():
        if row.strip() == "":
            continue
        genome, node = row.strip().split("\t")  # no header: every row is a genome
        assert node in nodes, f"{node!r} is not a node of the relation file"
        rows.append((genome, node))
    return parent_of, rows


def _emit(tmp_path: Path, **kwargs) -> tuple[Path, Path]:
    t2t, gmap, _ = _emit_relations(
        _BINARY,
        kwargs.pop("outgroup_leaf", "Out_Gen_sp_GCA_000009.1"),
        kwargs.pop("redundant", {"Fam_Gen_sp_GCA_000001.1": ["Fam_Gen_sp_GCA_000011.1"]}),
        node_basename=kwargs.pop("node_basename", None),
        root_name=kwargs.pop("root_name", FLEXTAXD_ROOT),
        remove_outgroup=kwargs.pop("remove_outgroup", False),
        out_tree2tax=tmp_path / "tree2tax.tsv",
        out_map=tmp_path / "genomes_map.tsv",
        logger=_LOG,
        **kwargs,
    )
    return t2t, gmap


def test_default_output_passes_the_mirror_reader(tmp_path: Path) -> None:
    t2t, gmap = _emit(tmp_path)
    parent_of, rows = _mirror_read(t2t, gmap)
    leaves = {c for c in parent_of if c not in parent_of.values()}
    assert leaves == {
        "Fam_Gen_sp_GCA_000001.1",
        "Fam_Gen_sp_GCA_000002.1",
        "Fam_Gen_sp_GCA_000003.1",
        "Fam_Gen_sp_GCA_000004.1",
        "Out_Gen_sp_GCA_000009.1",
    }
    assert ("GCA_000011.1", "Fam_Gen_sp_GCA_000001.1") in rows


@pytest.mark.parametrize(
    "kwargs",
    [
        {"node_basename": "N"},
        {"remove_outgroup": True},
        {"collapse_support": 0.5},
        {"collapse_length": 0.15, "node_basename": "C", "remove_outgroup": True},
        {"outgroup_leaf": None},
    ],
    ids=["basename", "remove-outgroup", "collapse-support", "collapse-both", "no-outgroup"],
)
def test_every_option_keeps_the_contract(tmp_path: Path, kwargs: dict) -> None:
    t2t, gmap = _emit(tmp_path, **kwargs)
    _mirror_read(t2t, gmap)


def test_another_top_name_breaks_the_contract_and_is_announced(tmp_path: Path, caplog) -> None:
    with caplog.at_level(logging.WARNING, logger="test-flextaxd"):
        t2t, gmap = _emit(tmp_path, root_name="top")
    assert any("FlexTaxD" in r.getMessage() for r in caplog.records)
    with pytest.raises(AssertionError, match="validates from"):
        _mirror_read(t2t, gmap)
