"""One figure style for every report: fonts, sizes and colours.

The palette is that of Okabe and Ito (2008), which stays distinguishable
with the common forms of colour-vision deficiency. Text is set in DejaVu
Sans, which ships with matplotlib, so a figure renders alike on every
machine; SVG keeps text as text (editable in a drawing program) and no
figure file carries a date, so the same data give the same file.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from typing import Any, cast

# Okabe-Ito: blue, orange, bluish green, vermillion, reddish purple, sky
# blue, yellow, grey.
PALETTE = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
    "#F0E442",
    "#999999",
)
PRIMARY = PALETTE[0]
MUTED = "#999999"
HIGHLIGHT = "#D55E00"

# Widths of a journal column and of the full page, in inches.
SINGLE_COLUMN = 3.35
DOUBLE_COLUMN = 6.85

RC = {
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 1.0,
    "axes.prop_cycle": None,  # set in context()
    "svg.fonttype": "none",
    "svg.hashsalt": "repgenr",
    "pdf.fonttype": 42,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
}


@contextmanager
def context() -> Iterator[None]:
    """Apply the report style while figures are created and saved."""
    import matplotlib as mpl
    from cycler import cycler

    rc: dict[str, Any] = dict(RC)
    rc["axes.prop_cycle"] = cycler(color=list(PALETTE))
    # matplotlib types the keys as a Literal of the known rc names.
    with mpl.rc_context(cast(Any, rc)):
        yield


def color(i: int) -> str:
    """The i-th palette colour, cycling."""
    return PALETTE[i % len(PALETTE)]
