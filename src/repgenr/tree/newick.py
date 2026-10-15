"""Newick helpers, including a small neighbor-joining implementation.

Used by the alignment-free sourmash tree builder, which produces a distance
matrix that must be turned into a tree without an external phylogenetics binary.
"""

from __future__ import annotations


def neighbor_joining(labels: list[str], dist: list[list[float]]) -> str:
    """Classic neighbor-joining; returns an unrooted Newick string.

    ``dist`` is a symmetric distance matrix aligned with ``labels``. A
    dict-of-dicts keyed by node id keeps the live distances as nodes merge.
    """
    n = len(labels)
    if n == 1:
        return f"({_sanitize(labels[0])}:0.0);"
    if n == 2:
        half = dist[0][1] / 2
        return f"({_sanitize(labels[0])}:{half:.6f},{_sanitize(labels[1])}:{half:.6f});"

    names: dict[int, str] = {i: _sanitize(label) for i, label in enumerate(labels)}
    d: dict[int, dict[int, float]] = {
        i: {j: dist[i][j] for j in range(n) if j != i} for i in range(n)
    }
    active: set[int] = set(range(n))
    next_id = n

    while len(active) > 2:
        m = len(active)
        net: dict[int, float] = {i: sum(d[i][j] for j in active if j != i) for i in active}

        best_pair: tuple[int, int] | None = None
        best_q = float("inf")
        active_list = sorted(active)
        for ai in range(m):
            for bj in range(ai + 1, m):
                i, j = active_list[ai], active_list[bj]
                q = (m - 2) * d[i][j] - net[i] - net[j]
                if q < best_q:
                    best_q = q
                    best_pair = (i, j)
        if best_pair is None:  # m >= 3 here, so some pair was always scored
            raise RuntimeError("Neighbor-joining found no pair to merge.")
        i, j = best_pair

        dij = d[i][j]
        di = 0.5 * dij + (net[i] - net[j]) / (2 * (m - 2))
        dj = dij - di
        di = max(di, 0.0)
        dj = max(dj, 0.0)

        new = next_id
        next_id += 1
        names[new] = f"({names[i]}:{di:.6f},{names[j]}:{dj:.6f})"

        d[new] = {}
        for k in active:
            if k in (i, j):
                continue
            dist_new_k = 0.5 * (d[i][k] + d[j][k] - dij)
            d[new][k] = dist_new_k
            d[k][new] = dist_new_k
        active.discard(i)
        active.discard(j)
        active.add(new)

    i, j = sorted(active)
    d_final = d[i][j]
    return f"({names[i]}:{d_final / 2:.6f},{names[j]}:{d_final / 2:.6f});"


_NEEDS_QUOTES = frozenset(" ()[]:,;'\t\r\n")


def _sanitize(label: str) -> str:
    """Make a label safe for Newick without changing it.

    A label holding a reserved character (space, parentheses, brackets, colon,
    comma, semicolon, quote) is written in single quotes with an embedded quote
    doubled, as the Newick grammar allows. Replacing the characters instead
    could map two distinct labels onto one leaf name, which no later step can
    undo; quoting is reversible and dendropy reads the label back as written.
    """
    if any(ch in _NEEDS_QUOTES for ch in label):
        return "'" + label.replace("'", "''") + "'"
    return label
