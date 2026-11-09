"""Report tables: TSV files and the summary statistics they show."""

from __future__ import annotations

import csv
import math
import statistics
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from ..core.contracts import atomic_replace
from .model import Table


def write_table(
    out_dir: Path,
    id: str,
    caption: str,
    columns: Sequence[str],
    rows: Iterable[Sequence[Any]],
    *,
    numeric: Sequence[str] = (),
    note: str = "",
) -> Table:
    """Write ``<id>.tsv`` and return its :class:`Table`. ``None`` is written
    as an empty cell."""
    file = f"{id}.tsv"
    n = 0
    with atomic_replace(out_dir / file, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(["" if v is None else v for v in row])
            n += 1
    return Table(
        id=id,
        caption=caption,
        file=file,
        columns=tuple(columns),
        n_rows=n,
        numeric=tuple(numeric),
        note=note,
    )


def read_rows(path: Path, limit: int | None = None) -> tuple[list[str], list[list[str]]]:
    """Header and (up to ``limit``) rows of a TSV file."""
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        header = next(reader, [])
        rows: list[list[str]] = []
        for row in reader:
            if limit is not None and len(rows) >= limit:
                break
            rows.append(row)
    return header, rows


def describe(values: Iterable[float | None]) -> dict[str, float | int | None]:
    """n, minimum, first quartile, median, third quartile and maximum of the
    values that are not None or NaN. Quartiles by the inclusive method."""
    data = sorted(v for v in values if v is not None and not math.isnan(v))
    if not data:
        return {"n": 0, "min": None, "q1": None, "median": None, "q3": None, "max": None}
    if len(data) == 1:
        q1 = q3 = data[0]
    else:
        q1, _, q3 = statistics.quantiles(data, n=4, method="inclusive")
    return {
        "n": len(data),
        "min": data[0],
        "q1": q1,
        "median": statistics.median(data),
        "q3": q3,
        "max": data[-1],
    }


def fmt(value: float | int | None, digits: int = 2) -> str:
    """A number for a table cell: integers with thousands separators, floats
    rounded; empty for None."""
    if value is None:
        return ""
    if isinstance(value, int) or (isinstance(value, float) and value.is_integer() and digits == 0):
        return f"{int(value):,}"
    return f"{value:,.{digits}f}"
