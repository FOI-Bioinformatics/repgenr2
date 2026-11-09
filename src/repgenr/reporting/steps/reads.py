"""Report of the run selection: runs, platforms, sketches and the screen.

Reads ``reads.tsv``, ``reads_stats.json`` (runs remaining after each
filter), ``reads_sketches/<variant>/reads_sketches.tsv`` and
``reads_screen/<variant>/reads_screen.tsv``.
"""

from __future__ import annotations

import csv
import statistics
from collections import Counter
from pathlib import Path

from ...core.contracts import (
    READS_SKETCH_STATUSES,
    READS_SKETCHES_DIR,
    READS_SKETCHES_TSV,
    READS_TSV,
    read_reads,
    read_reads_sketches,
)
from ...core.stage_stats import read_stats
from ...stages.reads import READS_STATS, READS_STATS_KIND
from ...stages.reads_screen import READS_SCREEN_DIR, READS_SCREEN_TSV
from .. import style
from ..figures import horizontal_bars, integer_ticks, log_scale, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import write_table


def _variant_files(workdir: Path, folder: str, name: str) -> list[Path]:
    root = workdir / folder
    return sorted(root.glob(f"*/{name}")) if root.is_dir() else []


def sources(workdir: Path) -> list[Path]:
    return [
        workdir / READS_TSV,
        workdir / READS_STATS,
        *_variant_files(workdir, READS_SKETCHES_DIR, READS_SKETCHES_TSV),
        *_variant_files(workdir, READS_SCREEN_DIR, READS_SCREEN_TSV),
    ]


def present(workdir: Path) -> bool:
    return (workdir / READS_TSV).is_file()


def collect(inp: StepInput) -> StepReport:
    wd, out = inp.workdir, inp.out_dir
    rows = read_reads(wd / READS_TSV)
    stats = read_stats(wd / READS_STATS, READS_STATS_KIND) or {}
    gb = [r.bases / 1e9 for r in rows if r.bases]
    platforms = Counter(r.platform or "unknown" for r in rows)

    report = StepReport(step="reads", title="Sequencing runs")
    report.key_numbers = [
        KeyNumber("Runs selected", len(rows)),
        KeyNumber("Biosamples", len({r.biosample for r in rows if r.biosample})),
        KeyNumber("Sequence selected", round(sum(gb), 1), "Gb"),
        KeyNumber("Median per run", round(statistics.median(gb), 2) if gb else None, "Gb"),
    ]
    text = f"{len(rows):,} sequencing runs were selected"
    if stats.get("candidates") is not None:
        text += f" of {stats['candidates']:,} whole-genome runs returned by ENA"
    text += "; " + ", ".join(f"{n:,} {p}" for p, n in platforms.most_common()) + "."
    report.summary.append(text)

    funnel = stats.get("funnel")
    if funnel:
        report.tables.append(
            write_table(
                out,
                "run_selection",
                "Runs remaining after each selection filter, in the order applied.",
                ["filter", "runs"],
                [(label, n) for label, n in funnel],
                numeric=["runs"],
            )
        )
    groups: dict[tuple[str, str, str], list[float]] = {}
    for r in rows:
        key = (r.platform or "unknown", r.instrument_model or "unknown", r.layout or "unknown")
        groups.setdefault(key, []).append(r.bases / 1e9)
    report.tables.append(
        write_table(
            out,
            "runs_by_platform",
            "Selected runs per platform, instrument and library layout, with the median "
            "sequence per run.",
            ["platform", "instrument", "layout", "runs", "median_gb"],
            [
                (*key, len(v), round(statistics.median(v), 3))
                for key, v in sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
            ],
            numeric=["runs", "median_gb"],
        )
    )
    if gb:
        report.figures.append(_bases_figure(out, rows, inp.formats))

    for path in _variant_files(wd, READS_SKETCHES_DIR, READS_SKETCHES_TSV):
        variant = path.parent.name
        counts = Counter(r.status for r in read_reads_sketches(path))
        report.tables.append(
            write_table(
                out,
                f"sketch_status_{_safe(variant)}",
                f"Reads sketches of variant '{variant}': runs per status.",
                ["status", "runs"],
                [(s, counts[s]) for s in READS_SKETCH_STATUSES if counts.get(s)],
                numeric=["runs"],
            )
        )
    for path in _variant_files(wd, READS_SCREEN_DIR, READS_SCREEN_TSV):
        _screen(report, out, path, inp)
    return report


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in name)


def _screen(report: StepReport, out: Path, path: Path, inp: StepInput) -> None:
    variant = path.parent.name
    with open(path, encoding="utf-8", newline="") as fo:
        table = list(csv.DictReader(fo, delimiter="\t"))
    decisions = Counter((r.get("decision", ""), r.get("reason", "")) for r in table)
    report.tables.append(
        write_table(
            out,
            f"screen_{_safe(variant)}",
            f"Reads screen of variant '{variant}': runs per decision and reason.",
            ["decision", "reason", "runs"],
            [(d, reason, n) for (d, reason), n in sorted(decisions.items())],
            numeric=["runs"],
        )
    )
    by_decision = Counter(r.get("decision", "") for r in table)
    passed = by_decision.get("pass", 0)
    report.summary.append(f"The reads screen ({variant}) passed {passed:,} of {len(table):,} runs.")
    fractions = []
    for r in table:
        try:
            fractions.append(float(r.get("fraction") or ""))
        except ValueError:
            continue
    record = inp.records.get("reads_screen")
    settings = (record.params.get("settings") or {}) if record else {}
    threshold = settings.get("min_fraction")
    if fractions:
        report.figures.append(_fraction_figure(out, variant, fractions, threshold, inp.formats))
    items = by_decision.most_common()
    if items:

        def draw(fig) -> None:
            ax = fig.add_subplot()
            horizontal_bars(
                ax,
                [k or "none" for k, _ in items],
                [v for _, v in items],
                xlabel="Runs",
            )

        files = render(
            out, f"screen_decisions_{_safe(variant)}", draw, height=1.8, formats=inp.formats
        )
        report.figures.append(
            FigureRef(
                id=f"screen_decisions_{_safe(variant)}",
                caption=f"Decisions of the reads screen ({variant}).",
                alt="Bar chart of screen decisions.",
                files=files,
            )
        )


def _fraction_figure(
    out: Path, variant: str, fractions: list[float], threshold, formats
) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.hist(fractions, bins=40, range=(0, 1), color=style.PRIMARY, linewidth=0)
        if isinstance(threshold, int | float):
            ax.axvline(float(threshold), color="black", lw=0.8, ls="--")
        ax.set_xlabel("Fraction of the sketch matching the top genome")
        ax.set_ylabel("Runs")
        integer_ticks(ax.yaxis)

    stem = f"screen_fraction_{_safe(variant)}"
    files = render(out, stem, draw, formats=formats)
    caption = (
        f"Fraction of each run's reads sketch matched by its top reference genome ({variant})."
    )
    if isinstance(threshold, int | float):
        caption += f" Dashed line: the minimum fraction for a pass ({threshold:g})."
    return FigureRef(id=stem, caption=caption, alt="Histogram of matched fractions.", files=files)


def _bases_figure(out: Path, rows, formats) -> FigureRef:
    import numpy as np

    by_platform: dict[str, list[float]] = {}
    for r in rows:
        if r.bases:
            by_platform.setdefault(r.platform or "unknown", []).append(r.bases / 1e9)
    allv = [v for vs in by_platform.values() for v in vs]
    lo, hi = min(allv), max(allv)
    bins = np.logspace(np.log10(lo * 0.9), np.log10(hi * 1.1), 30) if hi > lo else 10

    def draw(fig) -> None:
        ax = fig.add_subplot()
        names = sorted(by_platform, key=lambda p: -len(by_platform[p]))
        ax.hist(
            [by_platform[p] for p in names],
            bins=bins,
            stacked=True,
            color=[style.color(i) for i in range(len(names))],
            label=names,
            linewidth=0,
        )
        if hi > lo:
            log_scale(ax, "x")
        ax.set_xlabel("Sequence per run (Gb, log scale)" if hi > lo else "Sequence per run (Gb)")
        ax.set_ylabel("Runs")
        integer_ticks(ax.yaxis)
        if len(names) > 1:
            ax.legend()

    files = render(out, "bases_per_run", draw, formats=formats)
    return FigureRef(
        id="bases_per_run",
        caption="Sequence per selected run, by sequencing platform.",
        alt="Histogram of bases per run.",
        files=files,
    )
