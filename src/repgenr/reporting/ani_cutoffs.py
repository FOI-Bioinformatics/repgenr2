"""How many clusters a dereplication would give at other ANI thresholds.

The estimate uses the genome sketches (``sketches/<record name>.sig.zip``,
k=31, scaled=1000) and the rule of ``dereplicate --tool sourmash``: the ANI
of a pair is c^(1/k) for the average containment c, a pair sharing too few
hashes has no estimate, and genomes are clustered greedily, the genome with
the most neighbours at or above the threshold first, each claiming its
neighbours that are still free (:func:`greedy_cluster_count`). For the other
dereplicators the numbers are an estimate with this rule, not a rerun of
the tool.

The shared hashes of all pairs are counted by a matrix product of
genome-by-hash indicator matrices, in blocks of hashes, so the cost grows
with the square of the number of genomes; above :data:`MAX_GENOMES` no
estimate is made.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.contracts import record_name
from ..core.sourmash import MIN_SHARED_FRACTION, MIN_SHARED_HASHES
from .member_ani import KSIZE, SKETCH_SUFFIX, load_hashes

# ANI thresholds always tabulated; the thresholds of the run are added.
DEFAULT_THRESHOLDS = (
    0.95,
    0.97,
    0.98,
    0.99,
    0.995,
    0.997,
    0.998,
    0.999,
    0.9993,
    0.9995,
    0.9997,
    0.9999,
)
MAX_GENOMES = 5000
# Hash columns per block of the indicator matrix.
_BLOCK = 20000


@dataclass(frozen=True)
class AniMatrix:
    names: list[str]  # genome file names, in matrix order
    ani: npt.NDArray[np.float64]  # NaN where no estimate
    missing: list[str]  # genomes without a readable sketch


def ani_matrix(genomes: Sequence[str], sketches_dir: Path, ksize: int = KSIZE) -> AniMatrix:
    """Pairwise ANI estimates of ``genomes`` (file names) from their sketches."""
    names: list[str] = []
    sets: list[npt.NDArray[np.uint64]] = []
    missing: list[str] = []
    for name in genomes:
        hashes = load_hashes(sketches_dir / (record_name(name) + SKETCH_SUFFIX), ksize)
        if not hashes:
            missing.append(name)
            continue
        names.append(name)
        sets.append(np.fromiter(hashes, dtype=np.uint64, count=len(hashes)))
    n = len(names)
    if n == 0:
        return AniMatrix(names, np.zeros((0, 0)), missing)
    sizes = np.array([len(s) for s in sets], dtype=np.float64)
    universe, inverse = np.unique(np.concatenate(sets), return_inverse=True)
    rows = np.repeat(np.arange(n), [len(s) for s in sets])
    shared = np.zeros((n, n), dtype=np.float64)
    order = np.argsort(inverse, kind="stable")
    cols_sorted, rows_sorted = inverse[order], rows[order]
    bounds = np.searchsorted(cols_sorted, np.arange(0, len(universe) + _BLOCK, _BLOCK))
    for b in range(len(bounds) - 1):
        lo, hi = bounds[b], bounds[b + 1]
        if lo == hi:
            continue
        block = np.zeros((n, _BLOCK), dtype=np.float32)
        block[rows_sorted[lo:hi], cols_sorted[lo:hi] - b * _BLOCK] = 1.0
        shared += block @ block.T
    containment = (shared / sizes[:, None] + shared / sizes[None, :]) / 2
    with np.errstate(divide="ignore", invalid="ignore"):
        ani = np.power(containment, 1.0 / ksize)
    smaller = np.minimum(sizes[:, None], sizes[None, :])
    unknown = (shared < MIN_SHARED_HASHES) & (shared < MIN_SHARED_FRACTION * smaller)
    unknown |= shared == 0
    ani[unknown] = np.nan
    np.fill_diagonal(ani, 1.0)
    return AniMatrix(names, ani, missing)


def greedy_cluster_count(ani: npt.NDArray[np.float64], threshold: float) -> tuple[int, int, int]:
    """(clusters, singletons, largest cluster) under the greedy rule of the
    sourmash dereplicator (``dereplicators.sourmash._greedy_cluster``)."""
    n = ani.shape[0]
    if n == 0:
        return 0, 0, 0
    adj = np.where(np.isnan(ani), False, ani >= threshold)
    np.fill_diagonal(adj, False)
    order = np.lexsort((np.arange(n), -adj.sum(axis=1)))
    assigned = np.full(n, -1, dtype=np.int64)
    sizes: list[int] = []
    for i in order:
        if assigned[i] != -1:
            continue
        assigned[i] = i
        claim = adj[i] & (assigned == -1)
        assigned[claim] = i
        sizes.append(1 + int(claim.sum()))
    return len(sizes), sum(1 for s in sizes if s == 1), max(sizes)


@dataclass(frozen=True)
class CutoffRow:
    ani: float
    clusters: int
    singletons: int
    largest: int


def cutoff_table(
    matrix: AniMatrix, thresholds: Iterable[float], subset: Sequence[str] | None = None
) -> list[CutoffRow]:
    """Clusters at each threshold, for all genomes of ``matrix`` or ``subset``."""
    ani = matrix.ani
    if subset is not None:
        index = {name: i for i, name in enumerate(matrix.names)}
        keep = [index[g] for g in subset if g in index]
        ani = ani[np.ix_(keep, keep)]
    out = []
    for t in sorted(set(thresholds)):
        clusters, singletons, largest = greedy_cluster_count(ani, t)
        out.append(CutoffRow(t, clusters, singletons, largest))
    return out
