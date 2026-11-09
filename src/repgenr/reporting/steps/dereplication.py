"""Report of the dereplication: clusters, representatives and quality.

Reads ``derep/clusters.tsv``, ``derep/genome_status.tsv`` (with the reason
of each ``fail_qc`` genome), ``derep/cluster_summary.tsv``,
``derep/derep_stats.json`` and, for species names and assembly quality,
``selection.tsv``.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from ...core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    REFINED_TSV,
    SELECTION_TSV,
    ClusterSummaryRow,
    SelectionRow,
    read_cluster_summary,
    read_clusters,
    read_genome_status,
    read_genome_status_reasons,
    read_selection,
)
from ...core.stage_stats import read_stats
from ...dereplicators.base import STATUS_CONTAINED, STATUS_FAIL_QC, STATUS_REPRESENTATIVE
from ...stages.dereplicate import DEREP_STATS, DEREP_STATS_KIND
from .. import style
from ..figures import bar_height, horizontal_bars, integer_ticks, log_scale, render
from ..member_ani import member_anis
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import describe, write_table

TOP_CLUSTERS = 20
TOP_SPECIES = 30
STATUS_ORDER = (STATUS_REPRESENTATIVE, STATUS_CONTAINED, STATUS_FAIL_QC)
STATUS_COLORS = {
    STATUS_REPRESENTATIVE: style.PRIMARY,
    STATUS_CONTAINED: style.MUTED,
    STATUS_FAIL_QC: style.HIGHLIGHT,
}


def _derep(workdir: Path) -> Path:
    return workdir / "derep"


def sources(workdir: Path) -> list[Path]:
    d = _derep(workdir)
    return [
        d / CLUSTERS_TSV,
        d / GENOME_STATUS_TSV,
        d / CLUSTER_SUMMARY_TSV,
        d / REFINED_TSV,
        d / DEREP_STATS,
        workdir / SELECTION_TSV,
        # The genome sketches the member ANI is estimated from (a directory:
        # digested by file names, sizes and times).
        workdir / "sketches",
    ]


def present(workdir: Path) -> bool:
    return (_derep(workdir) / CLUSTERS_TSV).is_file()


def species_label(row: SelectionRow) -> str:
    """'Genus epithet', or the whole token when it already names the genus
    (a viral species token such as 'Hepatovirus-ahepa')."""
    if not row.species:
        return "unclassified"
    if not row.genus or row.species.startswith(row.genus):
        return row.species
    return f"{row.genus} {row.species}"


def collect(inp: StepInput) -> StepReport:
    wd, out = inp.workdir, inp.out_dir
    d = _derep(wd)
    clusters = read_clusters(d / CLUSTERS_TSV)
    status = read_genome_status(d / GENOME_STATUS_TSV)
    reasons = read_genome_status_reasons(d / GENOME_STATUS_TSV)
    summary = (
        read_cluster_summary(d / CLUSTER_SUMMARY_TSV) if (d / CLUSTER_SUMMARY_TSV).is_file() else []
    )
    stats = read_stats(d / DEREP_STATS, DEREP_STATS_KIND) or {}
    selection = (
        {r.filename: r for r in read_selection(wd / SELECTION_TSV)}
        if (wd / SELECTION_TSV).is_file()
        else {}
    )
    record = inp.records.get("dereplicate")
    tool = (record.tool if record else None) or stats.get("tool")

    counts = Counter(status.values())
    n_genomes = len(status) or sum(1 + len(m) for m in clusters.values())
    n_reps = len(clusters)
    sizes = [1 + len(members) for members in clusters.values()]
    singletons = sum(1 for s in sizes if s == 1)
    passed = n_genomes - counts.get(STATUS_FAIL_QC, 0)

    report = StepReport(step="dereplication", title="Dereplication")
    report.key_numbers = [
        KeyNumber("Genomes dereplicated", n_genomes),
        KeyNumber("Representatives", n_reps),
        KeyNumber("Contained genomes", counts.get(STATUS_CONTAINED, 0)),
        KeyNumber("Failed quality filter", counts.get(STATUS_FAIL_QC, 0)),
        KeyNumber("Singleton clusters", singletons),
        KeyNumber("Largest cluster", max(sizes) if sizes else 0, "genomes"),
    ]
    text = f"Of {n_genomes:,} genomes, {n_reps:,} were kept as representatives"
    if passed > 0:
        text += f" ({100 * n_reps / passed:.1f}% of the genomes that passed the quality filter)"
    text += "."
    if counts.get(STATUS_FAIL_QC):
        text += f" {counts[STATUS_FAIL_QC]:,} genomes failed the quality filter of the tool."
    verb = "holds" if singletons == 1 else "hold"
    text += (
        f" {singletons:,} of the {n_reps:,} clusters {verb} a single genome; the largest "
        f"holds {max(sizes) if sizes else 0:,}."
    )
    reduce = stats.get("reduce")
    if reduce:
        text += (
            f" Keeping one representative per {reduce['level']} reduced "
            f"{reduce['before']:,} representatives to {reduce['after']:,}."
        )
    target = stats.get("target_reps")
    if target:
        text += (
            f" A secondary ANI of {target['secondary_ani']:.4f} was chosen in "
            f"{target['iterations']} iterations to approach {target['requested']:,} "
            f"representatives ({target['achieved']:,} reached)."
        )
    if stats.get("keeper_swaps"):
        text += (
            f" The keeper rule replaced the representative chosen by the tool in "
            f"{stats['keeper_swaps']:,} clusters."
        )
    report.summary.append(text)

    # Table: status, with the fail_qc reasons.
    status_rows: list[tuple[str, str, int]] = [
        (s, "", counts[s]) for s in STATUS_ORDER if counts.get(s)
    ]
    reason_counts = Counter(
        part for g, r in reasons.items() if status.get(g) == STATUS_FAIL_QC for part in r.split(",")
    )
    status_rows += [(STATUS_FAIL_QC, reason, n) for reason, n in sorted(reason_counts.items())]
    report.tables.append(
        write_table(
            out,
            "genome_status",
            "Genomes per dereplication status. Rows with a reason break down the genomes "
            "that failed the quality filter; a genome can fail on more than one criterion.",
            ["status", "reason", "genomes"],
            status_rows,
            numeric=["genomes"],
        )
    )
    # Table: largest clusters.
    if summary:
        ranked = sorted(summary, key=lambda r: (-_size(r), r.representative))
        report.tables.append(_clusters_table(out, ranked[:TOP_CLUSTERS], len(summary)))
        multi = [r for r in ranked if r.n_species > 1]
        if multi:
            report.tables.append(
                _clusters_table(
                    out,
                    multi,
                    len(summary),
                    id="multispecies_clusters",
                    caption=(
                        "Clusters whose genomes carry more than one species name in the "
                        "selection, largest first."
                    ),
                )
            )
            report.notes.append(
                f"{len(multi):,} clusters join genomes of more than one species name; "
                "they are listed in a table above."
            )
    # Table + figure: genomes and representatives per species.
    per_species = _per_species(clusters, status, selection)
    if per_species:
        report.tables.append(
            write_table(
                out,
                "species_reduction",
                "Genomes and representatives per species, largest first.",
                ["species", "genomes", "representatives"],
                per_species,
                numeric=["genomes", "representatives"],
            )
        )

    report.figures.append(_size_figure(out, sizes, inp.formats))
    quality = [
        (r.completeness, r.contamination or 0.0, status.get(name, ""))
        for name, r in selection.items()
        if r.completeness is not None and name in status
    ]
    if len(quality) >= 2:
        thresholds = tool == "drep" and not _virus_mode(record)
        report.figures.append(_quality_figure(out, quality, thresholds, inp.formats))
    if reason_counts:
        report.figures.append(_reasons_figure(out, reason_counts, inp.formats))
    if per_species:
        report.figures.append(_species_figure(out, per_species, inp.formats))
    if not selection:
        report.notes.append(
            "selection.tsv was not found, so species names and assembly quality are not shown."
        )
    _member_ani(report, out, clusters, wd / "sketches", record, inp.formats)
    return report


def _member_ani(report: StepReport, out: Path, clusters, sketches: Path, record, formats) -> None:
    """Table and figure of the ANI of each member to its representative,
    estimated from the genome sketches."""
    if not any(clusters.values()):
        return
    if not sketches.is_dir():
        report.notes.append(
            "The ANI of members to their representatives is estimated from the genome "
            "sketches, and this working directory has none (repgenr sketch -wd writes them)."
        )
        return
    rows, missing = member_anis(clusters, sketches)
    known = [r.ani for r in rows if r.ani is not None]
    secondary = record.params.get("secondary_ani") if record else None
    ordered = sorted(rows, key=lambda r: (r.ani is not None, r.ani or 0.0, r.member))
    report.tables.append(
        write_table(
            out,
            "member_ani",
            "ANI of each member to its representative, estimated from the genome sketches "
            "(average containment at k=31, scaled=1000), lowest first; empty where the "
            "sketches share too few hashes or a sketch is missing.",
            ["representative", "member", "ani", "shared_hashes"],
            [
                (
                    r.representative,
                    r.member,
                    None if r.ani is None else round(r.ani, 5),
                    r.shared_hashes,
                )
                for r in ordered
            ],  # fmt: skip
            numeric=["ani", "shared_hashes"],
        )
    )
    if known:
        d = describe(known)
        report.key_numbers.append(KeyNumber("Lowest member ANI", round(d["min"] or 0.0, 4)))
        text = (
            f"Estimated from the genome sketches, members share a median ANI of "
            f"{d['median']:.4f} with their representative (lowest {d['min']:.4f})"
        )
        if isinstance(secondary, int | float):
            below = sum(a < secondary for a in known)
            text += (
                f"; {below:,} of {len(known):,} lie below the secondary threshold of {secondary:g}"
            )
        report.summary.append(text + ".")
    unknown = len(rows) - len(known)
    if unknown:
        report.notes.append(
            f"{unknown:,} of {len(rows):,} members have no sketch ANI to their representative "
            f"({missing:,} genomes without a sketch; the others share too few hashes)."
        )
    if len(known) >= 3:
        report.figures.append(_member_ani_figure(out, known, secondary, formats))


def _member_ani_figure(out: Path, values: list[float], secondary, formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.hist(values, bins=40, color=style.PRIMARY, linewidth=0)
        if isinstance(secondary, int | float):
            ax.axvline(float(secondary), color="black", lw=0.6, ls="--")
        ax.set_xlabel("ANI to the representative (sketch estimate)")
        ax.set_ylabel("Members")
        integer_ticks(ax.yaxis)

    files = render(out, "member_ani", draw, formats=formats)
    caption = (
        "ANI of each member to its representative, estimated from the genome sketches "
        "(average containment at k=31)."
    )
    if isinstance(secondary, int | float):
        caption += " Dashed line: the secondary ANI threshold of the dereplication."
    return FigureRef(
        id="member_ani",
        caption=caption,
        alt="Histogram of member-to-representative ANI.",
        files=files,
    )


def _size(row: ClusterSummaryRow) -> int:
    """Genomes in the cluster, representative included."""
    return row.n_genomes or 1 + row.n_members


def _virus_mode(record) -> bool:
    if record is None:
        return False
    value = (record.params.get("tool_args") or {}).get("virus")
    return str(value).lower() in ("1", "true", "yes")


def _clusters_table(
    out: Path,
    rows: list[ClusterSummaryRow],
    total: int,
    *,
    id: str = "largest_clusters",
    caption: str | None = None,
):
    caption = caption or (
        f"The {min(TOP_CLUSTERS, len(rows))} largest of {total:,} clusters, with the "
        "species names of their genomes and the quality of the representative."
    )
    return write_table(
        out,
        id,
        caption,
        [
            "representative",
            "genomes",
            "species_names",
            "species",
            "rep_completeness",
            "rep_contamination",
        ],
        [
            (
                r.representative,
                _size(r),
                r.n_species,
                r.species,
                r.rep_completeness,
                r.rep_contamination,
            )
            for r in rows
        ],
        numeric=["genomes", "species_names", "rep_completeness", "rep_contamination"],
    )


def _per_species(
    clusters: dict[str, list[str]], status: dict[str, str], selection: dict[str, SelectionRow]
) -> list[tuple[str, int, int]]:
    if not selection:
        return []
    genomes: Counter[str] = Counter()
    reps: Counter[str] = Counter()
    for rep, members in clusters.items():
        for name in (rep, *members):
            row = selection.get(name)
            if row is not None:
                genomes[species_label(row)] += 1
        row = selection.get(rep)
        if row is not None:
            reps[species_label(row)] += 1
    for name, st in status.items():
        if st == STATUS_FAIL_QC and name in selection:
            genomes[species_label(selection[name])] += 1
    return sorted(
        ((sp, n, reps.get(sp, 0)) for sp, n in genomes.items()), key=lambda r: (-r[1], r[0])
    )


def _size_figure(out: Path, sizes: list[int], formats) -> FigureRef:
    counts = Counter(sizes)

    def draw(fig) -> None:
        ax = fig.add_subplot()
        xs = sorted(counts)
        if len(xs) <= 30:
            ax.bar(xs, [counts[x] for x in xs], width=0.8, color=style.PRIMARY, linewidth=0)
            ax.set_xlabel("Genomes per cluster")
            integer_ticks(ax.xaxis)
        else:
            import numpy as np

            bins = np.unique(np.logspace(0, np.log10(max(xs) + 1), 30).astype(int))
            ax.hist(sizes, bins=bins, color=style.PRIMARY, linewidth=0)
            log_scale(ax, "x")
            ax.set_xlabel("Genomes per cluster (log scale)")
        ax.set_ylabel("Clusters")
        integer_ticks(ax.yaxis)
        if max(counts.values()) > 50 * min(counts.values()):
            log_scale(ax, "y")
            ax.set_ylabel("Clusters (log scale)")

    files = render(out, "cluster_sizes", draw, formats=formats)
    return FigureRef(
        id="cluster_sizes",
        caption=(
            "Distribution of cluster sizes: number of clusters per number of genomes in "
            "the cluster, representative included."
        ),
        alt="Bar chart of the number of clusters by cluster size.",
        files=files,
    )


def _quality_figure(
    out: Path, points: list[tuple[float, float, str]], thresholds: bool, formats
) -> FigureRef:
    from ...dereplicators.drep import DREP_MAX_CONTAMINATION, DREP_MIN_COMPLETENESS

    def draw(fig) -> None:
        ax = fig.add_subplot()
        # Distinct markers, so genomes at the same values stay visible: a
        # cross is drawn over a filled circle.
        markers = {STATUS_CONTAINED: "o", STATUS_REPRESENTATIVE: "o", STATUS_FAIL_QC: "x"}
        sizes = {STATUS_CONTAINED: 10, STATUS_REPRESENTATIVE: 16, STATUS_FAIL_QC: 22}
        for st in (STATUS_CONTAINED, STATUS_REPRESENTATIVE, STATUS_FAIL_QC):
            pts = [(c, k) for c, k, s in points if s == st]
            if pts:
                ax.scatter(
                    [p[0] for p in pts],
                    [p[1] for p in pts],
                    s=sizes[st],
                    marker=markers[st],
                    alpha=0.6 if st == STATUS_CONTAINED else 0.9,
                    color=STATUS_COLORS[st],
                    linewidths=0.9 if st == STATUS_FAIL_QC else 0,
                    label=st.replace("_", " "),
                )
        if thresholds:
            ax.axvline(DREP_MIN_COMPLETENESS, color="black", lw=0.6, ls="--")
            ax.axhline(DREP_MAX_CONTAMINATION, color="black", lw=0.6, ls="--")
        ax.set_xlabel("Completeness (%)")
        ax.set_ylabel("Contamination (%)")
        ax.legend(loc="best", markerscale=1.5)

    files = render(out, "quality_by_status", draw, height=2.8, formats=formats)
    caption = (
        "Estimated completeness and contamination of each genome, coloured by its "
        "dereplication status."
    )
    if thresholds:
        caption += (
            f" Dashed lines: the dRep quality filter (completeness at least "
            f"{DREP_MIN_COMPLETENESS:g}%, contamination at most {DREP_MAX_CONTAMINATION:g}%)."
        )
    return FigureRef(
        id="quality_by_status",
        caption=caption,
        alt="Scatter plot of contamination against completeness by status.",
        files=files,
    )


def _reasons_figure(out: Path, reasons: Counter[str], formats) -> FigureRef:
    items = sorted(reasons.items(), key=lambda kv: (-kv[1], kv[0]))

    def draw(fig) -> None:
        ax = fig.add_subplot()
        horizontal_bars(
            ax,
            [k for k, _ in items],
            [v for _, v in items],
            colors=[style.HIGHLIGHT] * len(items),
            xlabel="Genomes",
        )

    files = render(out, "fail_qc_reasons", draw, height=bar_height(len(items)), formats=formats)
    return FigureRef(
        id="fail_qc_reasons",
        caption="Criteria of the quality filter that genomes failed; a genome can fail on more "
        "than one.",
        alt="Bar chart of quality filter criteria.",
        files=files,
    )


def _species_figure(out: Path, rows: list[tuple[str, int, int]], formats) -> FigureRef:
    shown = rows[:TOP_SPECIES]

    def draw(fig) -> None:
        import numpy as np

        ax = fig.add_subplot()
        y = np.arange(len(shown))[::-1]
        ax.barh(y + 0.2, [r[1] for r in shown], height=0.4, color=style.MUTED, label="genomes")
        ax.barh(
            y - 0.2,
            [r[2] for r in shown],
            height=0.4,
            color=style.PRIMARY,
            label="representatives",
        )
        ax.set_yticks(y)
        ax.set_yticklabels([r[0] for r in shown], fontstyle="italic")
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Genomes")
        integer_ticks(ax.xaxis)
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)

    files = render(
        out,
        "species_reduction",
        draw,
        height=bar_height(len(shown), per_bar=0.2),
        formats=formats,
    )
    caption = "Genomes and representatives per species"
    if len(rows) > TOP_SPECIES:
        caption += f", for the {TOP_SPECIES} species with most genomes"
    return FigureRef(
        id="species_reduction",
        caption=caption + ".",
        alt="Paired bar chart of genomes and representatives per species.",
        files=files,
    )
