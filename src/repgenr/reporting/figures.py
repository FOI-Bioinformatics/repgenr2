"""Drawing and saving report figures.

Figures are built with the object API of matplotlib (``Figure``), never with
``pyplot``, so no display backend is needed. :func:`render` applies the
report style while a figure is drawn and saved, and writes each format
atomically (a PNG at 300 dpi for screens, SVG or PDF for manuscripts).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from ..core.contracts import atomic_path
from . import style

DEFAULT_FORMATS = ("png", "svg")
FORMATS = ("png", "svg", "pdf")
PNG_DPI = 300
# Labels longer than this are cut in figures; tables keep the full text.
MAX_LABEL = 40

# Without these, matplotlib writes the creation date into SVG and PDF files.
_METADATA: dict[str, dict[str, Any]] = {
    "svg": {"Date": None},
    "pdf": {"CreationDate": None},
    "png": {},
}


def render(
    out_dir: Path,
    stem: str,
    draw: Callable[[Any], None],
    *,
    width: float = style.SINGLE_COLUMN,
    height: float = 2.6,
    formats: Sequence[str] = DEFAULT_FORMATS,
) -> dict[str, str]:
    """Draw a figure with ``draw(fig)`` and save it as ``<stem>.<format>``.

    Returns format -> file name, relative to ``out_dir``.
    """
    from matplotlib.figure import Figure

    files: dict[str, str] = {}
    with style.context():
        fig = Figure(figsize=(width, height), layout="constrained")
        draw(fig)
        for fmt in formats:
            name = f"{stem}.{fmt}"
            with atomic_path(out_dir / name) as tmp:
                # The temporary name ends in .tmp, so the format is explicit.
                fig.savefig(tmp, format=fmt, dpi=PNG_DPI, metadata=_METADATA[fmt])
            files[fmt] = name
    return files


def truncate(label: str, limit: int = MAX_LABEL) -> str:
    """``label`` cut to ``limit`` characters, ending in '...' when cut."""
    return label if len(label) <= limit else label[: limit - 3].rstrip() + "..."


def bar_height(n_bars: int, *, per_bar: float = 0.16, base: float = 0.7) -> float:
    """Figure height for a horizontal bar chart of ``n_bars`` bars."""
    return max(1.6, base + per_bar * n_bars)


def horizontal_bars(
    ax: Any,
    labels: Sequence[str],
    values: Sequence[float],
    *,
    colors: Sequence[str] | None = None,
    xlabel: str,
    annotate: bool = True,
) -> None:
    """Largest bar on top, values written beside the bars."""
    y = list(range(len(labels)))[::-1]
    bars = ax.barh(
        y, values, color=list(colors) if colors else style.PRIMARY, height=0.7, linewidth=0
    )
    ax.set_yticks(y)
    ax.set_yticklabels([truncate(label) for label in labels])
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)
    if annotate:
        ax.bar_label(bars, labels=[f"{v:,.0f}" for v in values], padding=2, fontsize=6)
        ax.margins(x=0.12)
