"""Report of glance: the distribution of pairwise similarities.

Reads ``glance_stats.json``, which the glance stage writes with a summary
and the histogram counts (the pairwise table itself is removed unless
``--keep-files``), and the dereplicate record for its ANI thresholds.
"""

from __future__ import annotations

from pathlib import Path

from ...core.stage_stats import read_stats
from ...stages.glance import GLANCE_STATS, GLANCE_STATS_KIND
from .. import style
from ..figures import render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import write_table


def sources(workdir: Path) -> list[Path]:
    return [workdir / GLANCE_STATS]


def present(workdir: Path) -> bool:
    return (workdir / GLANCE_STATS).is_file()


def collect(inp: StepInput) -> StepReport:
    stats = read_stats(inp.workdir / GLANCE_STATS, GLANCE_STATS_KIND) or {}
    measure = stats.get("measure") or "similarity"
    pairs = int(stats.get("pairs") or 0)
    report = StepReport(step="glance", title="Pairwise similarity (glance)")
    report.key_numbers = [
        KeyNumber("Genomes compared", stats.get("genomes")),
        KeyNumber("Genome pairs", pairs),
        KeyNumber(f"Median {measure}", _r(stats.get("median"))),
    ]
    if pairs == 0:
        report.summary.append("The comparison gave no pairwise similarity within the range.")
        return report
    report.summary.append(
        f"{stats.get('tool', 'The tool')} compared {stats.get('genomes', 0):,} genomes in "
        f"{pairs:,} pairs; the {measure} ranged from {stats['min']:.4f} to "
        f"{stats['max']:.4f}, with a median of {stats['median']:.4f}."
    )
    report.tables.append(
        write_table(
            inp.out_dir,
            "pairwise_similarity",
            f"Summary of the pairwise {measure} between genomes (each pair counted once).",
            ["measure", "genomes", "pairs", "min", "q1", "median", "q3", "max"],
            [
                (
                    measure,
                    stats.get("genomes"),
                    pairs,
                    *(_r(stats.get(k)) for k in ("min", "q1", "median", "q3", "max")),
                )
            ],
            numeric=["genomes", "pairs", "min", "q1", "median", "q3", "max"],
        )
    )
    hist = stats.get("histogram")
    if hist:
        thresholds = _thresholds(inp)
        report.figures.append(_histogram(inp.out_dir, hist, measure, thresholds, inp.formats))
    return report


def _r(value) -> float | None:
    return None if value is None else round(float(value), 4)


def _thresholds(inp: StepInput) -> dict[str, float]:
    record = inp.records.get("dereplicate")
    if record is None:
        return {}
    out = {}
    for key, label in (("primary_ani", "primary"), ("secondary_ani", "secondary")):
        value = record.params.get(key)
        if isinstance(value, int | float):
            out[label] = float(value)
    return out


def _histogram(
    out: Path, hist: dict, measure: str, thresholds: dict[str, float], formats
) -> FigureRef:
    counts, edges = hist["counts"], hist["edges"]

    def draw(fig) -> None:
        ax = fig.add_subplot()
        widths = [b - a for a, b in zip(edges[:-1], edges[1:], strict=True)]
        ax.bar(edges[:-1], counts, width=widths, align="edge", color=style.PRIMARY, linewidth=0)
        lo, hi = edges[0], edges[-1]
        for i, (label, value) in enumerate(sorted(thresholds.items())):
            if lo <= value <= hi:
                ax.axvline(value, color=style.color(i + 1), lw=0.8, ls="--", label=label)
        if any(lo <= v <= hi for v in thresholds.values()):
            ax.legend(title="dereplication ANI", loc="upper left")
        ax.set_xlabel(measure)
        ax.set_ylabel("Genome pairs")

    files = render(out, "similarity_histogram", draw, formats=formats)
    caption = f"Distribution of the pairwise {measure} ({len(counts)} bins)."
    if thresholds:
        caption += " Dashed lines: the ANI thresholds of the dereplication, where in range."
    return FigureRef(
        id="similarity_histogram",
        caption=caption,
        alt=f"Histogram of pairwise {measure}.",
        files=files,
    )
