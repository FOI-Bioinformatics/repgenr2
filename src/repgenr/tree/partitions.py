"""The partition sidecar of a concatenated alignment.

When the core-gene route concatenates per-gene alignments itself, it knows
where each gene lies. ``partitions.json`` keeps those column ranges beside
the alignment so a tree builder can fit one model per gene. Like
``alignment_stats.json`` the sidecar names the alignment and holds its
digest, so it is only ever applied to the very file it describes.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import atomic_path
from ..core.inputs import file_digest
from ..coregenes.select import Partition

PARTITIONS = "partitions.json"
_PARTITIONS_VERSION = 1


@dataclass(frozen=True)
class PartitionSet:
    alphabet: str
    partitions: list[Partition]


def write_partitions(alignment: Path, partitions: Sequence[Partition], alphabet: str) -> Path:
    """Write the sidecar beside ``alignment`` and return its path."""
    path = alignment.parent / PARTITIONS
    payload = {
        "v": _PARTITIONS_VERSION,
        "alignment": alignment.name,
        "alignment_digest": file_digest(alignment),
        "alphabet": alphabet,
        "partitions": [[p.gene, p.start, p.end] for p in partitions],
    }
    with atomic_path(path) as tmp:
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_partitions(alignment: Path) -> PartitionSet | None:
    """The partitions beside ``alignment`` when the sidecar describes that very
    file (same name, same content), else None."""
    path = alignment.parent / PARTITIONS
    if not path.is_file() or not alignment.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or data.get("alignment") != alignment.name:
        return None
    if data.get("alignment_digest") != file_digest(alignment):
        return None
    try:
        parts = [Partition(str(g), int(s), int(e)) for g, s, e in data["partitions"]]
        return PartitionSet(str(data["alphabet"]), parts)
    except (KeyError, TypeError, ValueError):
        return None


def copy_partitions(src_alignment: Path, dst_alignment: Path) -> Path | None:
    """Re-issue the sidecar of ``src_alignment`` for a copy of it at
    ``dst_alignment``, or None when there is none."""
    found = read_partitions(src_alignment)
    if found is None:
        return None
    return write_partitions(dst_alignment, found.partitions, found.alphabet)
