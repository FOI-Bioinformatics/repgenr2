"""Per-gene selection for the core-gene route.

A core-gene builder may hand over one aligned FASTA per gene instead of its
own concatenation. This module reads those alignments, computes a few
per-gene metrics (occupancy, gap fraction, variable and parsimony-informative
sites, mean pairwise identity), ranks and caps the genes, concatenates the
kept ones in a fixed order and records where each gene lies in the result,
so the tree builder can be given one partition per gene.

Pure Python and numpy; no external tool runs here.
"""

from __future__ import annotations

import random
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core.contracts import atomic_path
from ..core.errors import UserInputError, WorkdirError

RANKS = ("occupancy", "informative", "variable", "length", "random")
ALPHABETS = ("dna", "protein")

# The GeneMetrics field each rank sorts on (descending); random has none.
_RANK_FIELDS = {
    "occupancy": "occupancy",
    "informative": "informative_sites",
    "variable": "variable_sites",
    "length": "length",
}
_GAP = ord("-")
# Characters that are no state in a column: gaps and the missing-data symbols.
_NOT_STATE = {"dna": b"-N?.", "protein": b"-X?."}


@dataclass(frozen=True)
class GeneSelection:
    """Which genes of the per-gene alignments go into the concatenation.

    ``max_genes`` 0 keeps every gene that passes the filters; ``min_length``
    and ``max_gaps`` (fraction of gap characters, None for no limit) filter
    before ranking; ``seed`` only matters for ``rank="random"``.
    """

    max_genes: int = 0
    rank: str = "occupancy"
    min_length: int = 0
    max_gaps: float | None = None
    seed: int = 0

    def __post_init__(self) -> None:
        if self.rank not in RANKS:
            raise UserInputError(
                f"Unknown gene rank '{self.rank}'. Choose from: {', '.join(RANKS)}."
            )
        if self.max_genes < 0:
            raise UserInputError(
                f"The maximum number of core genes must be 0 (all) or more, got {self.max_genes}"
            )
        if self.min_length < 0:
            raise UserInputError(
                f"The minimum gene length must be 0 or more, got {self.min_length}"
            )
        if self.max_gaps is not None and not 0 <= self.max_gaps <= 1:
            raise UserInputError(
                f"The maximum gap fraction of a gene must be between 0 and 1, got {self.max_gaps}"
            )


@dataclass(frozen=True)
class GeneMetrics:
    """Metrics of one gene alignment.

    A genome counts as present only when its record holds at least one
    character that is neither a gap nor a missing-data symbol; a record of
    gaps (or N) only is treated as absent, for occupancy and gap fraction.
    """

    length: int  # alignment columns
    present: int  # genomes with a record of at least one called character
    occupancy: float  # present / genomes in the set
    gap_fraction: float  # '-' over all characters of the present records
    variable_sites: int  # columns with at least two states
    informative_sites: int  # columns with at least two states seen twice or more
    mean_identity: float  # mean per-column fraction of identical pairs


@dataclass(frozen=True)
class Partition:
    """One gene's columns in the concatenation, 1-based and inclusive."""

    gene: str
    start: int
    end: int


@dataclass(frozen=True)
class GeneRow:
    """One row of ``core_genes.tsv``. ``rank`` is the 1-based position in the
    rank order and ``start``/``end`` the gene's columns, all None for a gene
    left out."""

    gene: str
    kept: bool
    rank: int | None
    metrics: GeneMetrics
    start: int | None
    end: int | None


def _check_alphabet(alphabet: str) -> None:
    if alphabet not in ALPHABETS:
        raise ValueError(f"Unknown alphabet '{alphabet}'; expected one of {', '.join(ALPHABETS)}")


def _record_name(header: str) -> str:
    """The genome a per-gene record belongs to.

    Panaroo names records ``<isolate>;<gene_id>`` and MAFFT prefixes ``_R_``
    to a sequence it reverse-complemented.
    """
    first = header.split()[0] if header.split() else ""
    return first.removeprefix("_R_").split(";")[0]


def load_gene(path: Path, alphabet: str = "dna") -> dict[str, str]:
    """Read one per-gene alignment as {genome: upper-case sequence}.

    A genome with two records (a paralog the builder kept) is given the one
    with the fewest gaps. Records must all be of one length.
    """
    _check_alphabet(alphabet)
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name, parts = _record_name(line[1:]), []
            else:
                parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(parts)))
    # Checked over every record, a duplicate dropped below included.
    lengths = {len(seq) for _, seq in records}
    if len(lengths) > 1:
        raise WorkdirError(
            f"{path}: the records of a gene alignment differ in length "
            f"({min(lengths)} to {max(lengths)}); the builder's alignment is not usable."
        )
    seqs: dict[str, str] = {}
    for rec, seq in records:
        if not rec:
            raise WorkdirError(f"{path}: a record has no name in its header line")
        seq = seq.upper()
        if rec not in seqs or seq.count("-") < seqs[rec].count("-"):
            seqs[rec] = seq
    return seqs


def gene_metrics(seqs: dict[str, str], n_genomes: int, alphabet: str = "dna") -> GeneMetrics:
    """The metrics of one gene alignment over a set of ``n_genomes`` genomes.

    Gaps and missing-data symbols are no state. A column is variable with two
    or more states and informative when at least two states occur twice or
    more. Mean identity is, over the columns with two or more called
    characters, the fraction of pairs of them that agree; 1.0 when no column
    qualifies.
    """
    _check_alphabet(alphabet)
    if not seqs:
        return GeneMetrics(0, 0, 0.0, 0.0, 0, 0, 1.0)
    flat = np.frombuffer("".join(seqs.values()).encode("ascii", "replace"), dtype=np.uint8)
    length = len(next(iter(seqs.values())))
    arr = flat.reshape(len(seqs), length)

    not_state = np.zeros(256, dtype=bool)
    not_state[list(_NOT_STATE[alphabet])] = True
    called = ~not_state[arr]
    # A record without one called character is no evidence of the gene.
    rows = called.any(axis=1)
    present = int(rows.sum())
    occupancy = present / n_genomes if n_genomes else 0.0
    total = present * length
    gap_fraction = float(np.count_nonzero(arr[rows] == _GAP) / total) if total else 0.0
    # Per-state counts per column, one vectorised pass per symbol in use.
    states = [int(s) for s in np.unique(arr[called])]
    counts = np.zeros((len(states), length), dtype=np.int64)
    for i, s in enumerate(states):
        counts[i] = np.count_nonzero(arr == s, axis=0)
    n_states = np.count_nonzero(counts >= 1, axis=0)
    n_repeated = np.count_nonzero(counts >= 2, axis=0)
    n_col = called.sum(axis=0).astype(np.int64)
    pairs = n_col * (n_col - 1) // 2
    same = (counts * (counts - 1) // 2).sum(axis=0)
    usable = n_col >= 2
    mean_identity = float((same[usable] / pairs[usable]).mean()) if usable.any() else 1.0
    return GeneMetrics(
        length=length,
        present=present,
        occupancy=occupancy,
        gap_fraction=gap_fraction,
        variable_sites=int(np.count_nonzero(n_states >= 2)),
        informative_sites=int(np.count_nonzero(n_repeated >= 2)),
        mean_identity=mean_identity,
    )


def rank_genes(metrics: dict[str, GeneMetrics], selection: GeneSelection) -> list[str]:
    """The genes kept by ``selection``, in rank order.

    Filters first (``min_length``, ``max_gaps``), then sorts by the rank key
    descending with ties broken by length descending and gene id ascending;
    ``random`` shuffles with the selection's seed. Cut to ``max_genes``
    when that is above 0. A gene of no columns is always left out: it would
    add an empty partition.
    """
    passing = sorted(
        g
        for g, m in metrics.items()
        if m.length >= max(selection.min_length, 1)
        and (selection.max_gaps is None or m.gap_fraction <= selection.max_gaps)
    )
    if selection.rank == "random":
        random.Random(selection.seed).shuffle(passing)
        ranked = passing
    else:
        field = _RANK_FIELDS[selection.rank]
        ranked = sorted(
            passing,
            key=lambda g: (-getattr(metrics[g], field), -metrics[g].length, g),
        )
    if selection.max_genes > 0:
        ranked = ranked[: selection.max_genes]
    return ranked


def concatenate(
    kept: Sequence[str],
    genes: Mapping[str, dict[str, str]],
    genomes: Sequence[str],
    out: Path,
) -> list[Partition]:
    """Write the concatenation of the ``kept`` genes to ``out``.

    Genes are laid out by gene id, so the columns do not depend on the rank
    that chose them. Every genome in ``genomes`` gets a record, in that
    order; a genome without a gene is given gaps over that gene's width.
    Returns each gene's columns.
    """
    blocks: dict[str, list[str]] = {g: [] for g in genomes}
    partitions: list[Partition] = []
    start = 1
    for gene in sorted(kept):
        seqs = genes[gene]
        width = len(next(iter(seqs.values()))) if seqs else 0
        gaps = "-" * width
        for genome in genomes:
            blocks[genome].append(seqs.get(genome, gaps))
        partitions.append(Partition(gene, start, start + width - 1))
        start += width
    with atomic_path(out) as tmp, open(tmp, "w", encoding="utf-8") as fo:
        for genome in genomes:
            fo.write(f">{genome}\n{''.join(blocks[genome])}\n")
    return partitions


_TABLE_HEADER = (
    "gene",
    "kept",
    "rank",
    "length",
    "present",
    "occupancy",
    "gap_fraction",
    "variable_sites",
    "informative_sites",
    "mean_identity",
    "start",
    "end",
)


def _cell(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_gene_table(path: Path, rows: Sequence[GeneRow]) -> None:
    """Write ``core_genes.tsv``: one row per gene the builder aligned."""
    lines = ["\t".join(_TABLE_HEADER)]
    for row in rows:
        m = row.metrics
        values: tuple[object, ...] = (
            row.gene,
            row.kept,
            row.rank,
            m.length,
            m.present,
            float(m.occupancy),
            float(m.gap_fraction),
            m.variable_sites,
            m.informative_sites,
            float(m.mean_identity),
            row.start,
            row.end,
        )
        lines.append("\t".join(_cell(v) for v in values))
    with atomic_path(path) as tmp:
        tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
