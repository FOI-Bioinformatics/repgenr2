"""Column statistics of a nucleotide alignment for ascertainment-bias correction.

A tree built from variable sites only overestimates branch lengths unless the
tree builder is told how many invariant sites the alignment lost. The counts
of constant A, C, G and T columns of the whole-genome alignment are taken in
one pass over the file and kept beside the variable-site alignment in
``alignment_stats.json``; IQ-TREE reads them as ``-fconst`` and RAxML-NG as
``+ASC_STAM``. Both builders also want no invariant column left in the
alignment they are given, so :func:`strip_invariant_columns` removes any that
masking or a typer left and adds them to the counts.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.contracts import atomic_path

ALIGNMENT_STATS = "alignment_stats.json"
_STATS_VERSION = 1
_BASES = "ACGT"

# One bit per base, so a column's bases across all genomes are the OR of its
# bytes; anything else (N, gaps, other ambiguity codes) is 0 and ignored.
_BITS = np.zeros(256, dtype=np.uint8)
for _bit, _base in enumerate(_BASES):
    _BITS[ord(_base)] = _BITS[ord(_base.lower())] = 1 << _bit
_POPCOUNT = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)


def _records(path: Path):
    """Yield (name, sequence bytes) for each FASTA record; lines are joined."""
    name: str | None = None
    parts: list[bytes] = []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(parts)
                name = line[1:].split()[0].decode("utf-8", "replace") if line[1:].split() else ""
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, b"".join(parts)


def column_bases(path: Path, logger: logging.Logger | None = None) -> npt.NDArray[np.uint8]:
    """For each column, the OR of the base bits of every genome (one pass).

    Records of unequal length are compared over the shortest; the rest is
    reported once, since a whole-genome alignment should be rectangular.
    """
    state: npt.NDArray[np.uint8] | None = None
    shortest: int | None = None
    for _name, seq in _records(path):
        bits = _BITS[np.frombuffer(seq, dtype=np.uint8)]
        if state is None:
            state = bits.copy()
            shortest = len(bits)
            continue
        n = min(len(bits), len(state))
        if n != len(state) or n != len(bits):
            shortest = min(shortest or n, n)
        state[:n] |= bits[:n]
    if state is None:
        return np.zeros(0, dtype=np.uint8)
    if shortest is not None and shortest < len(state):
        if logger is not None:
            logger.warning(
                "%s: records differ in length; columns beyond %d are not counted",
                path.name,
                shortest,
            )
        state = state[:shortest]
    return state


def constant_site_counts(
    path: Path, logger: logging.Logger | None = None
) -> tuple[tuple[int, int, int, int], int]:
    """(constant A, C, G, T column counts, variable column count) of ``path``.

    A column is constant when exactly one of A, C, G and T occurs in it, in
    however many genomes; a column with two or more bases is variable; a
    column with no base at all is neither.
    """
    state = column_bases(path, logger)
    counts = tuple(int(np.count_nonzero(state == (1 << bit))) for bit in range(4))
    variable = int(np.count_nonzero(_POPCOUNT[state] >= 2))
    return (counts[0], counts[1], counts[2], counts[3]), variable


def strip_invariant_columns(
    msa: Path, out: Path, logger: logging.Logger | None = None
) -> tuple[int, int, int, int]:
    """Write ``msa`` without its invariant columns to ``out``; return the
    constant A, C, G, T columns removed, to be added to the counts given to
    the tree builder. Columns with no base are dropped and not counted."""
    records = list(_records(msa))
    if not records:
        with atomic_path(out) as tmp:
            tmp.write_text("", encoding="utf-8")
        return (0, 0, 0, 0)
    length = min(len(seq) for _, seq in records)
    state = np.zeros(length, dtype=np.uint8)
    for _, seq in records:
        state |= _BITS[np.frombuffer(seq[:length], dtype=np.uint8)]
    keep = _POPCOUNT[state] >= 2
    removed = tuple(int(np.count_nonzero(state == (1 << bit))) for bit in range(4))
    if logger is not None and not keep.all():
        logger.info(
            "%s: %d of %d columns are invariant and are counted instead of aligned",
            msa.name,
            int(length - keep.sum()),
            length,
        )
    with atomic_path(out) as tmp, open(tmp, "wb") as fo:
        for name, seq in records:
            kept = np.frombuffer(seq[:length], dtype=np.uint8)[keep].tobytes()
            fo.write(b">" + name.encode("utf-8") + b"\n" + kept + b"\n")
    return (removed[0], removed[1], removed[2], removed[3])


def write_alignment_stats(
    beside: Path,
    counts: tuple[int, int, int, int],
    variable: int,
    *,
    source: str,
    masked: bool,
) -> Path:
    """Write the sidecar next to the variable-site alignment ``beside``."""
    path = beside.parent / ALIGNMENT_STATS
    payload = {
        "v": _STATS_VERSION,
        "source": source,
        "constant_sites": dict(zip(_BASES, counts, strict=True)),
        "variable_sites": variable,
        "masked": masked,
    }
    with atomic_path(path) as tmp:
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_constant_sites(beside: Path) -> tuple[int, int, int, int] | None:
    """The constant-site counts recorded beside ``beside``, or None."""
    path = beside.parent / ALIGNMENT_STATS
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        counts = data["constant_sites"]
        return tuple(int(counts[b]) for b in _BASES)  # type: ignore[return-value]
    except (OSError, ValueError, KeyError, TypeError):
        return None


def format_counts(counts: tuple[int, int, int, int], sep: str) -> str:
    return sep.join(str(c) for c in counts)


def add_counts(
    a: tuple[int, int, int, int], b: tuple[int, int, int, int]
) -> tuple[int, int, int, int]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])
