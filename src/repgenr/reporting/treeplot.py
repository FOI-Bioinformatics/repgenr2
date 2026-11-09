"""A rectangular tree figure drawn with matplotlib from a dendropy tree.

The tree is ladderised, so the figure does not depend on the child order
of the Newick file. Branch lengths give the horizontal position (a
phylogram); without them every edge has length one (a cladogram). Lines
are drawn as one ``LineCollection``, so a tree of thousands of leaves
draws in seconds. Support values are written at the nodes for small
trees and shown as filled circles for larger ones.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from ..tree.support import parse_support
from . import style

# Leaf labels are left out above this many leaves (a colour strip remains);
# no figure is drawn above MAX_LEAVES.
MAX_LABELLED = 400
MAX_LEAVES = 20_000
# Support values are written as numbers up to this many leaves.
MAX_SUPPORT_NUMBERS = 100
# Leaf groups (species) coloured in the strip; the others are grey.
MAX_GROUPS = 8


@dataclass(frozen=True)
class TreeLayout:
    """Node coordinates: x is the distance from the root, y the leaf rank."""

    x: dict[int, float]
    y: dict[int, float]
    leaves: list[Any]
    phylogram: bool


def layout(tree: Any) -> TreeLayout:
    """Coordinates of every node of a ladderised ``tree`` (changed in place)."""
    tree.ladderize(ascending=True)
    phylogram = any(
        e.length is not None and e.length > 0 for e in tree.postorder_edge_iter() if e.tail_node
    )
    x: dict[int, float] = {}
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            x[id(node)] = 0.0
            continue
        length = node.edge.length if phylogram else 1.0
        x[id(node)] = x[id(node.parent_node)] + max(float(length or 0.0), 0.0)
    leaves = list(tree.leaf_node_iter())
    y: dict[int, float] = {id(leaf): float(i) for i, leaf in enumerate(leaves)}
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            ys = [y[id(c)] for c in node.child_node_iter()]
            y[id(node)] = (min(ys) + max(ys)) / 2
    return TreeLayout(x=x, y=y, leaves=leaves, phylogram=phylogram)


def segments(tree: Any, lay: TreeLayout) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Line segments of the rectangular tree: one horizontal per edge and one
    vertical per internal node."""
    segs = []
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        if parent is not None:
            segs.append(((lay.x[id(parent)], lay.y[id(node)]), (lay.x[id(node)], lay.y[id(node)])))
        if not node.is_leaf():
            ys = [lay.y[id(c)] for c in node.child_node_iter()]
            segs.append(((lay.x[id(node)], min(ys)), (lay.x[id(node)], max(ys))))
    return segs


def figure_height(n_leaves: int) -> float:
    """Figure height in inches: 0.12 per leaf, between 3 and 20."""
    return min(20.0, max(3.0, 0.12 * n_leaves + 0.8))


def _scale_length(depth: float) -> float:
    """A round length near a tenth of ``depth`` for the scale bar."""
    import math

    if depth <= 0:
        return 0.0
    target = depth / 10
    power = 10 ** math.floor(math.log10(target))
    for step in (1, 2, 5, 10):
        if step * power >= target:
            return step * power
    return 10 * power


def draw_tree(
    fig: Any,
    tree: Any,
    *,
    labels: Mapping[str, str],
    groups: Mapping[str, str],
    outgroup: str | None,
    unit: str,
) -> dict[str, Any]:
    """Draw ``tree`` onto ``fig``. ``labels`` and ``groups`` (species) are
    keyed by leaf taxon label. Returns what was drawn (for captions and
    tests): leaves, labelled, support_mode, groups_shown, scale."""
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    lay = layout(tree)
    n = len(lay.leaves)
    ax = fig.add_subplot()
    ax.add_collection(LineCollection(segments(tree, lay), colors="black", linewidths=0.6))
    depth = max(lay.x.values()) if lay.x else 0.0
    labelled = n <= MAX_LABELLED

    # Species groups, largest first, coloured at the tips.
    counts: dict[str, int] = {}
    for leaf in lay.leaves:
        g = groups.get(_name(leaf), "")
        if g:
            counts[g] = counts.get(g, 0) + 1
    shown = [g for g, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))][:MAX_GROUPS]
    colour = {g: style.color(i) for i, g in enumerate(shown)}
    tip_x = [lay.x[id(leaf)] for leaf in lay.leaves]
    tip_y = [lay.y[id(leaf)] for leaf in lay.leaves]
    tip_c = [colour.get(groups.get(_name(leaf), ""), style.MUTED) for leaf in lay.leaves]
    if len(shown) > 1:
        ax.scatter(tip_x, tip_y, s=8 if n <= 200 else 3, c=tip_c, linewidths=0, zorder=3)

    pad = depth * 0.01
    if labelled:
        size = 6 if n <= 60 else 5 if n <= 150 else 4
        for leaf in lay.leaves:
            name = _name(leaf)
            is_out = outgroup is not None and name == outgroup
            ax.text(
                lay.x[id(leaf)] + pad,
                lay.y[id(leaf)],
                labels.get(name, name),
                va="center",
                ha="left",
                fontsize=size,
                color=style.HIGHLIGHT if is_out else "black",
                fontweight="bold" if is_out else "normal",
            )

    supports: list[tuple[Any, float]] = []
    for node in tree.internal_nodes():
        value = parse_support(node.label) if node.parent_node is not None else None
        if value is not None:
            supports.append((node, value))
    scale100 = any(s > 1.0 for _, s in supports)
    hi, mid = (95.0, 70.0) if scale100 else (0.95, 0.70)
    mode = "none"
    if supports and n <= MAX_SUPPORT_NUMBERS:
        mode = "numbers"
        for node, s in supports:
            ax.text(
                lay.x[id(node)] - pad,
                lay.y[id(node)] + 0.15,
                f"{s:.0f}" if scale100 else f"{s:.2f}",
                ha="right",
                va="bottom",
                fontsize=4.5,
                color="#444444",
            )
    elif supports:
        mode = "dots"
        for threshold, fill in ((hi, "black"), (mid, "#999999")):
            pts = [
                (lay.x[id(node)], lay.y[id(node)])
                for node, s in supports
                if s >= threshold and (fill == "black" or s < hi)
            ]
            if pts:
                ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=4, c=fill, zorder=4)

    scale = _scale_length(depth) if lay.phylogram else 0.0
    if scale:
        y0 = -max(1.0, n * 0.02)
        ax.plot([0, scale], [y0, y0], color="black", lw=0.8)
        ax.text(
            scale / 2,
            y0 - max(0.6, n * 0.012),
            f"{scale:g} {unit}",
            ha="center",
            va="top",
            fontsize=6,
        )
    ax.set_ylim(-max(2.5, n * 0.05), n - 0.5 + max(0.5, n * 0.01))
    right = depth * (1.45 if labelled else 1.05) if depth else 1.0
    ax.set_xlim(-depth * 0.02, right)
    ax.axis("off")
    if len(shown) > 1:
        handles = [
            Line2D([], [], marker="o", ls="", color=colour[g], markersize=4, label=g) for g in shown
        ]
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.0, 0.0),
            ncol=2,
            fontsize=6,
            prop={"style": "italic", "size": 6},
        )
    return {
        "leaves": n,
        "labelled": labelled,
        "support_mode": mode,
        "support_scale": 100 if scale100 else 1,
        "groups_shown": len(shown) if len(shown) > 1 else 0,
        "groups_total": len(counts),
        "scale": scale,
        "phylogram": lay.phylogram,
    }


def _name(leaf: Any) -> str:
    return leaf.taxon.label if leaf.taxon is not None else ""
