"""The partition sidecar beside a concatenated alignment."""

from __future__ import annotations

import json
from pathlib import Path

from repgenr.coregenes.select import Partition
from repgenr.tree.partitions import (
    PARTITIONS,
    PartitionSet,
    copy_partitions,
    read_partitions,
    write_partitions,
)

_PARTS = [Partition("gA", 1, 3), Partition("gB", 4, 9)]


def _aln(path: Path, seq: str = "ACGTACGTA") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">a\n{seq}\n>b\n{seq}\n", encoding="utf-8")
    return path


def test_write_and_read(tmp_path: Path) -> None:
    aln = _aln(tmp_path / "core.fasta")
    path = write_partitions(aln, _PARTS, "dna")
    assert path == tmp_path / PARTITIONS
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["v"] == 1 and data["alignment"] == "core.fasta" and data["alphabet"] == "dna"
    assert data["partitions"] == [["gA", 1, 3], ["gB", 4, 9]]
    assert isinstance(data["alignment_digest"], str) and data["alignment_digest"]
    assert read_partitions(aln) == PartitionSet("dna", _PARTS)


def test_read_refuses_another_file_or_changed_content(tmp_path: Path) -> None:
    aln = _aln(tmp_path / "core.fasta")
    assert read_partitions(aln) is None
    write_partitions(aln, _PARTS, "dna")
    other = _aln(tmp_path / "other.fasta")
    assert read_partitions(other) is None
    _aln(aln, "ACGTACGTT")
    assert read_partitions(aln) is None
    (tmp_path / PARTITIONS).write_text("not json", encoding="utf-8")
    assert read_partitions(other) is None


def test_copy_reissues_the_sidecar_for_the_copy(tmp_path: Path) -> None:
    src = _aln(tmp_path / "a" / "core.fasta")
    assert copy_partitions(src, tmp_path / "b" / "x.fasta") is None
    write_partitions(src, _PARTS, "dna")
    dst = _aln(tmp_path / "b" / "tree_input.fasta")
    assert copy_partitions(src, dst) == tmp_path / "b" / PARTITIONS
    assert read_partitions(dst) == PartitionSet("dna", _PARTS)
