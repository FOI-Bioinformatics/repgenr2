"""Report of the assemblies: contiguity, quality, coverage and excused runs.

Reads ``assembly_stats.tsv`` and ``excused_runs.tsv``.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from ...core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    AssemblyStatsRow,
    read_assembly_stats,
    read_excused_runs,
)
from .. import style
from ..figures import integer_ticks, log_scale, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import describe, write_table

# Excuse reasons are free text (a tool's error tail); grouped by their start.
REASON_CHARS = 80


def sources(workdir: Path) -> list[Path]:
    return [workdir / ASSEMBLY_STATS_TSV, workdir / EXCUSED_RUNS_TSV]


def present(workdir: Path) -> bool:
    return (workdir / ASSEMBLY_STATS_TSV).is_file() or (workdir / EXCUSED_RUNS_TSV).is_file()


def _median_iqr(values, digits: int = 1) -> str:
    d = describe(values)
    if not d["n"]:
        return ""
    return f"{d['median']:.{digits}f} ({d['q1']:.{digits}f}-{d['q3']:.{digits}f})"


def collect(inp: StepInput) -> StepReport:
    wd, out = inp.workdir, inp.out_dir
    rows = (
        read_assembly_stats(wd / ASSEMBLY_STATS_TSV) if (wd / ASSEMBLY_STATS_TSV).is_file() else []
    )
    excused = read_excused_runs(wd / EXCUSED_RUNS_TSV) if (wd / EXCUSED_RUNS_TSV).is_file() else []

    report = StepReport(step="assembly", title="Genome assembly")
    n50 = [r.n50 / 1e3 for r in rows if r.n50]
    length = [r.total_length / 1e6 for r in rows if r.total_length]
    report.key_numbers = [
        KeyNumber("Assemblies", len(rows)),
        KeyNumber("Runs without a genome", len(excused)),
        KeyNumber("Median N50", round(describe(n50)["median"] or 0, 1) if n50 else None, "kb"),
        KeyNumber(
            "Median length", round(describe(length)["median"] or 0, 2) if length else None, "Mb"
        ),
    ]
    text = f"{len(rows):,} runs were assembled"
    if excused:
        text += f" and {len(excused):,} gave no genome"
    text += "."
    if n50:
        text += f" The median N50 was {describe(n50)['median']:.1f} kb"
        text += f" and the median length {describe(length)['median']:.2f} Mb."
    report.summary.append(text)

    if rows:
        by_tool: dict[str, list[AssemblyStatsRow]] = {}
        for r in rows:
            by_tool.setdefault(r.assembler or "unknown", []).append(r)
        report.tables.append(
            write_table(
                out,
                "assembly_metrics",
                "Assembly metrics per assembler: median (first-third quartile).",
                [
                    "assembler",
                    "assemblies",
                    "contigs",
                    "length_mb",
                    "n50_kb",
                    "coverage_x",
                    "completeness_pct",
                    "contamination_pct",
                ],
                [
                    (
                        tool,
                        len(rs),
                        _median_iqr([r.n_contigs for r in rs], 0),
                        _median_iqr([r.total_length / 1e6 for r in rs], 2),
                        _median_iqr([r.n50 / 1e3 for r in rs], 1),
                        _median_iqr([r.est_coverage for r in rs], 1),
                        _median_iqr([r.completeness for r in rs], 1),
                        _median_iqr([r.contamination for r in rs], 2),
                    )
                    for tool, rs in sorted(by_tool.items())
                ],
                numeric=["assemblies"],
            )
        )
        flags = Counter(r.taxonomy_flag or "none" for r in rows)
        if set(flags) != {"none"}:
            report.tables.append(
                write_table(
                    out,
                    "taxonomy_flags",
                    "Assemblies per taxonomy flag (agreement between the run metadata and "
                    "the classification of the assembly).",
                    ["taxonomy_flag", "assemblies"],
                    flags.most_common(),
                    numeric=["assemblies"],
                )
            )
        report.figures.append(_contiguity_figure(out, rows, inp.formats))
        quality = [r for r in rows if r.completeness is not None]
        if len(quality) >= 2:
            report.figures.append(_quality_figure(out, quality, inp.formats))
        coverage = [r.est_coverage for r in rows if r.est_coverage]
        if len(coverage) >= 2:
            report.figures.append(_coverage_figure(out, coverage, inp.formats))
    if excused:
        grouped = Counter((e.step, " ".join(e.reason.split())[:REASON_CHARS]) for e in excused)
        report.tables.append(
            write_table(
                out,
                "excused_runs",
                "Runs that gave no genome, per step and reason (reasons cut to "
                f"{REASON_CHARS} characters).",
                ["step", "reason", "runs"],
                [(s, r, n) for (s, r), n in sorted(grouped.items(), key=lambda kv: -kv[1])],
                numeric=["runs"],
            )
        )
    return report


def _contiguity_figure(out: Path, rows: list[AssemblyStatsRow], formats) -> FigureRef:
    flags = sorted({r.taxonomy_flag or "none" for r in rows})

    def draw(fig) -> None:
        ax = fig.add_subplot()
        for i, flag in enumerate(flags):
            pts = [r for r in rows if (r.taxonomy_flag or "none") == flag and r.n50 > 0]
            ax.scatter(
                [r.total_length / 1e6 for r in pts],
                [r.n50 / 1e3 for r in pts],
                s=12,
                alpha=0.8,
                color=style.color(i),
                linewidths=0,
                label=flag,
            )
        log_scale(ax, "y")
        ax.set_xlabel("Total length (Mb)")
        ax.set_ylabel("N50 (kb, log scale)")
        if len(flags) > 1:
            ax.legend(title="taxonomy flag")

    files = render(out, "contiguity", draw, height=2.6, formats=formats)
    return FigureRef(
        id="contiguity",
        caption="N50 against total length of each assembly"
        + (", coloured by taxonomy flag." if len(flags) > 1 else "."),
        alt="Scatter plot of N50 against assembly length.",
        files=files,
    )


def _quality_figure(out: Path, rows: list[AssemblyStatsRow], formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.scatter(
            [r.completeness for r in rows],
            [r.contamination or 0.0 for r in rows],
            s=12,
            alpha=0.8,
            color=style.PRIMARY,
            linewidths=0,
        )
        ax.set_xlabel("Completeness (%)")
        ax.set_ylabel("Contamination (%)")

    files = render(out, "assembly_quality", draw, formats=formats)
    return FigureRef(
        id="assembly_quality",
        caption="Estimated completeness and contamination of each assembly (CheckM2).",
        alt="Scatter plot of contamination against completeness.",
        files=files,
    )


def _coverage_figure(out: Path, coverage: list[float], formats) -> FigureRef:
    import numpy as np

    lo, hi = min(coverage), max(coverage)

    def draw(fig) -> None:
        ax = fig.add_subplot()
        if hi > lo:
            bins = np.logspace(np.log10(lo * 0.9), np.log10(hi * 1.1), 30)
            ax.hist(coverage, bins=bins, color=style.PRIMARY, linewidth=0)
            log_scale(ax, "x")
            ax.set_xlabel("Estimated coverage (x, log scale)")
        else:
            ax.hist(coverage, bins=10, color=style.PRIMARY, linewidth=0)
            ax.set_xlabel("Estimated coverage (x)")
        ax.set_ylabel("Assemblies")
        integer_ticks(ax.yaxis)

    files = render(out, "coverage", draw, formats=formats)
    return FigureRef(
        id="coverage",
        caption="Estimated read coverage of each assembly.",
        alt="Histogram of estimated coverage.",
        files=files,
    )
