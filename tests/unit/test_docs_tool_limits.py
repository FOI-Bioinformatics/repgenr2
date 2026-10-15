"""docs/choosing-tools.md lists each adapter's declared genome limit.

The table is parsed by its header row and compared with the registries, so a
changed limit or a new adapter cannot leave the page out of date.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from repgenr.aligners.base import registry as aligners
from repgenr.annotators.base import registry as annotators
from repgenr.assemblers.base import registry as assemblers
from repgenr.classifiers.base import registry as classifiers
from repgenr.dereplicators.base import registry as dereplicators
from repgenr.maskers.base import registry as maskers
from repgenr.polishers.base import registry as polishers
from repgenr.snptypers.base import registry as snptypers
from repgenr.treebuilders.base import registry as treebuilders

DOC = Path(__file__).resolve().parents[2] / "docs" / "choosing-tools.md"
HEADER = re.compile(r"^\|\s*Family\s*\|\s*Tool\s*\|\s*Declared limit\s*\|\s*$")

# Family names as the table writes them.
REGISTRIES = {
    "dereplicator": dereplicators,
    "aligner": aligners,
    "snptyper": snptypers,
    "masker": maskers,
    "treebuilder": treebuilders,
    "assembler": assemblers,
    "classifier": classifiers,
    "polisher": polishers,
    "annotator": annotators,
}


def _rows() -> list[tuple[str, str, str]]:
    lines = DOC.read_text(encoding="utf-8").splitlines()
    start = next((i for i, line in enumerate(lines) if HEADER.match(line)), None)
    assert start is not None, "limits table header not found in docs/choosing-tools.md"
    rows: list[tuple[str, str, str]] = []
    for line in lines[start + 2 :]:  # skip the header and the separator
        if not line.startswith("|"):
            break
        cells = [c.strip().strip("`") for c in line.strip().strip("|").split("|")]
        assert len(cells) == 3, f"malformed row: {line}"
        rows.append((cells[0], cells[1], cells[2]))
    return rows


def _declared(family: str, tool: str) -> str:
    limit = REGISTRIES[family].get(tool).capabilities.recommended_max_genomes
    return "unbounded" if limit is None else str(limit)


def test_table_matches_declared_limits() -> None:
    rows = _rows()
    assert rows, "limits table is empty"
    for family, tool, limit in rows:
        assert family in REGISTRIES, f"unknown family {family!r}"
        assert limit == _declared(family, tool), f"{family}/{tool}: doc says {limit}"


@pytest.mark.parametrize("family", sorted(REGISTRIES))
def test_every_registered_tool_is_listed(family: str) -> None:
    listed = {tool for fam, tool, _ in _rows() if fam == family}
    assert listed == set(REGISTRIES[family].names())
