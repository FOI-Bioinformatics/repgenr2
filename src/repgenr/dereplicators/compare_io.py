"""Writers for the glance comparison contract (pairwise table + dendrogram).

Adapters that compute a dense similarity matrix themselves (rather than
calling a tool that already writes these files, as dRep does) use these
helpers to produce the files :class:`~repgenr.dereplicators.base.CompareResult`
points at:

* a pairwise CSV with the columns ``genome1``, ``genome2``, ``similarity``,
  one row per unordered genome pair;
* a dendrogram PDF from average-linkage (UPGMA) clustering on the distance
  ``1 - similarity``, the clustering dRep uses for its primary dendrogram.

The clustering is implemented here with numpy (nearest-neighbour chain, O(N^2)
time and one N x N matrix of memory) so that glance needs no dependency beyond
numpy and matplotlib. The linkage matrix uses the SciPy layout
(``[child_a, child_b, height, size]`` per merge, merges sorted by height,
new clusters numbered from N), so it can be checked against
``scipy.cluster.hierarchy.linkage(method="average")``.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

# A PDF page taller than 200 inches exceeds what common viewers accept.
_MAX_FIG_HEIGHT_IN = 200.0
_LEAF_SPACING_IN = 0.15


def _check_square(sim: npt.NDArray[np.float64], names: Sequence[str]) -> None:
    if sim.ndim != 2 or sim.shape[0] != sim.shape[1]:
        raise ValueError(f"similarity matrix must be square, got shape {sim.shape}")
    if sim.shape[0] != len(names):
        raise ValueError(f"{len(names)} names for a {sim.shape[0]} x {sim.shape[0]} matrix")


def write_pairwise_csv(names: Sequence[str], sim: npt.NDArray[np.float64], path: Path) -> Path:
    """Write one ``genome1,genome2,similarity`` row per unordered pair (i < j).

    An unknown similarity (NaN; for sourmash, a pair whose shared hashes are
    compatible with chance) is written as an empty cell.
    """
    sim = np.asarray(sim, dtype=float)
    _check_square(sim, names)
    rows, cols = np.triu_indices(len(names), k=1)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, lineterminator="\n")
        writer.writerow(["genome1", "genome2", "similarity"])
        for i, j in zip(rows.tolist(), cols.tolist(), strict=True):
            value = sim[i, j]
            writer.writerow([names[i], names[j], "" if np.isnan(value) else f"{value:.6g}"])
    return path


def average_linkage(dist: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Average-linkage (UPGMA) clustering of a square distance matrix.

    Returns an (N-1) x 4 linkage matrix in the SciPy layout. Average linkage
    is reducible, so the nearest-neighbour chain finds the same merges as the
    naive algorithm in O(N^2) time.
    """
    d = np.array(dist, dtype=float, copy=True)
    n = d.shape[0]
    if d.ndim != 2 or d.shape[1] != n:
        raise ValueError(f"distance matrix must be square, got shape {d.shape}")
    if n < 2:
        raise ValueError("average linkage needs at least two items")
    d = (d + d.T) / 2.0  # tolerate tiny asymmetries from the producer
    np.fill_diagonal(d, np.inf)
    size = np.ones(n)
    active = np.ones(n, dtype=bool)
    merges: list[tuple[int, int, float]] = []
    chain: list[int] = []
    while len(merges) < n - 1:
        if not chain:
            chain.append(int(np.flatnonzero(active)[0]))
        while True:
            a = chain[-1]
            row = np.where(active, d[a], np.inf)
            row[a] = np.inf
            b = int(np.argmin(row))
            # Prefer the previous chain element on a tie, so the chain ends
            # at a reciprocal pair instead of cycling between equals.
            if len(chain) > 1 and row[chain[-2]] <= row[b]:
                b = chain[-2]
            if len(chain) > 1 and b == chain[-2]:
                break
            chain.append(b)
        a = chain.pop()
        b = chain.pop()
        i, j = min(a, b), max(a, b)
        height = float(d[i, j])
        merged = (size[i] * d[i] + size[j] * d[j]) / (size[i] + size[j])
        d[i, :] = merged
        d[:, i] = merged
        d[j, :] = np.inf
        d[:, j] = np.inf
        d[i, i] = np.inf
        active[j] = False
        size[i] += size[j]
        merges.append((i, j, height))

    # Slot i always holds the cluster that contains item i, so each merge
    # joins "the cluster of i" with "the cluster of j". Sort by height (stable,
    # so a child created before an equal-height parent stays first) and
    # renumber with a union-find over the items.
    order = sorted(range(len(merges)), key=lambda k: merges[k][2])
    parent = list(range(n))
    cluster_id = list(range(n))  # union-find root -> current cluster number
    cluster_size = [1] * n

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    z = np.zeros((n - 1, 4))
    for step, k in enumerate(order):
        i, j, height = merges[k]
        ri, rj = find(i), find(j)
        ci, cj = cluster_id[ri], cluster_id[rj]
        z[step] = (min(ci, cj), max(ci, cj), height, cluster_size[ri] + cluster_size[rj])
        parent[rj] = ri
        cluster_id[ri] = n + step
        cluster_size[ri] += cluster_size[rj]
    return z


def leaf_order(z: npt.NDArray[np.float64]) -> list[int]:
    """Left-to-right leaf order of a linkage matrix (iterative, no recursion limit)."""
    n = z.shape[0] + 1
    out: list[int] = []
    stack = [2 * n - 2]
    while stack:
        node = stack.pop()
        if node < n:
            out.append(node)
            continue
        left, right = int(z[node - n, 0]), int(z[node - n, 1])
        stack.append(right)
        stack.append(left)
    return out


def write_dendrogram(
    names: Sequence[str],
    sim: npt.NDArray[np.float64],
    path: Path,
    measure: str = "ANI",
) -> list[str]:
    """Draw an average-linkage dendrogram of ``1 - sim`` to ``path`` (PDF).

    Leaves are labelled with ``names``; returns the leaf labels from top to
    bottom of the figure. Needs at least two names. An unknown similarity
    (NaN) is drawn at distance 1.0.
    """
    from matplotlib.collections import LineCollection
    from matplotlib.figure import Figure

    sim = np.asarray(sim, dtype=float)
    _check_square(sim, names)
    dist = np.clip(1.0 - np.nan_to_num(sim, nan=0.0), 0.0, 1.0)
    z = average_linkage(dist)
    n = len(names)
    order = leaf_order(z)

    # Leaf k of the order sits at y = k; an internal node sits midway
    # between its children, at x = its merge height.
    pos: dict[int, tuple[float, float]] = {leaf: (0.0, float(k)) for k, leaf in enumerate(order)}
    segments: list[list[tuple[float, float]]] = []
    for step in range(n - 1):
        a, b, height = int(z[step, 0]), int(z[step, 1]), float(z[step, 2])
        (xa, ya), (xb, yb) = pos[a], pos[b]
        segments.append([(xa, ya), (height, ya)])
        segments.append([(xb, yb), (height, yb)])
        segments.append([(height, ya), (height, yb)])
        pos[n + step] = (height, (ya + yb) / 2.0)

    height_in = min(_MAX_FIG_HEIGHT_IN, max(4.0, _LEAF_SPACING_IN * n + 1.5))
    # Keep labels from overlapping once the page height is capped.
    fontsize = min(8.0, 0.8 * 72.0 * (height_in - 1.5) / n)
    longest = max(len(s) for s in names)
    width_in = 6.0 + min(6.0, 0.06 * longest * fontsize)
    fig = Figure(figsize=(width_in, height_in))
    ax = fig.add_subplot()
    ax.add_collection(LineCollection(segments, colors="black", linewidths=0.6))
    top = max(float(z[-1, 2]), 1e-6)
    ax.set_xlim(top * 1.05, 0.0)  # leaves on the right, root on the left
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    ax.yaxis.tick_right()
    ax.set_yticks(range(n))
    labels = [names[leaf] for leaf in order]
    ax.set_yticklabels(labels, fontsize=fontsize)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(f"1 - {measure} (average linkage)")
    ax.set_title(f"{n} genomes, average-linkage clustering on 1 - {measure}")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    return labels
