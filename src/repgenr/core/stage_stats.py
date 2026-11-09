"""Counts a stage reports that its data files do not hold.

Some numbers of a run exist only while a stage runs: how many candidate
genomes a ``--limit`` cut chose from, how many runs a filter removed. A stage
writes them to a small versioned JSON file beside its outputs, so the report
of the step (``repgenr report``) can show them. The file names the kind of
numbers it holds; a reader that finds another kind or version treats the file
as absent. A stage rewrites its file on every run, so numbers of an earlier
run never remain beside new outputs.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from .contracts import atomic_path

STATS_VERSION = 1


def write_stats(path: Path, kind: str, values: Mapping[str, Any]) -> Path:
    """Write ``values`` (JSON-serialisable) to ``path`` under ``kind``."""
    payload = {"v": STATS_VERSION, "kind": kind, "values": dict(values)}
    with atomic_path(path) as tmp:
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_stats(path: Path, kind: str) -> dict[str, Any] | None:
    """The values written under ``kind``, or None when the file is absent,
    unreadable, of another kind or of another version."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(data, dict) or data.get("v") != STATS_VERSION or data.get("kind") != kind:
        return None
    values = data.get("values")
    return dict(values) if isinstance(values, dict) else None
