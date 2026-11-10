"""Report of the genome selection: taxa, counts and assembly quality.

Reads ``selection.tsv`` and, when present, ``metadata_summary.tsv`` (which
names viral species as the entry stage did), ``selection_stats.json``
(candidates before a ``--limit`` cut), ``missing_accessions.txt``,
``outgroup_choice.tsv`` (the candidates of an automatic outgroup choice) and
the genome folder.
"""

from __future__ import annotations

import csv
from pathlib import Path

from ...core.contracts import (
    MISSING_ACCESSIONS_TXT,
    SELECTION_TSV,
    SelectionRow,
    list_fasta,
    read_selection,
)
from ...core.outgroup_choice import OUTGROUP_CHOICE_TSV, read_choice
from ...core.selection_summary import (
    SELECTION_STATS,
    SELECTION_STATS_KIND,
    SUMMARY_TSV,
    count_taxa,
)
from ...core.stage_stats import read_stats
from .. import style
from ..figures import bar_height, horizontal_bars, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import describe, write_table

# Species drawn in the bar chart; the table lists all.
TOP_SPECIES = 30


def sources(workdir: Path) -> list[Path]:
    return [
        workdir / SELECTION_TSV,
        workdir / SUMMARY_TSV,
        workdir / SELECTION_STATS,
        workdir / MISSING_ACCESSIONS_TXT,
        workdir / OUTGROUP_CHOICE_TSV,
    ]


def present(workdir: Path) -> bool:
    return (workdir / SELECTION_TSV).is_file()


def _taxa(workdir: Path, rows: list[SelectionRow]) -> list[tuple[str, str, int]]:
    # The summary written by the entry stage names viral species correctly
    # (the whole species token); counting again would add the genus.
    path = workdir / SUMMARY_TSV
    if path.is_file():
        with open(path, encoding="utf-8", newline="") as fo:
            reader = csv.DictReader(fo, delimiter="\t")
            return [(r["level"], r["name"], int(r["genomes"])) for r in reader]
    return count_taxa(rows)


def _missing(workdir: Path) -> int | None:
    path = workdir / MISSING_ACCESSIONS_TXT
    if not path.is_file():
        return None
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.strip())


def collect(inp: StepInput) -> StepReport:
    wd, out = inp.workdir, inp.out_dir
    rows = read_selection(wd / SELECTION_TSV)
    ingroup = [r for r in rows if not r.is_outgroup]
    outgroup = [r for r in rows if r.is_outgroup]
    taxa = _taxa(wd, rows)
    per_level = {lvl: [(n, c) for lv, n, c in taxa if lv == lvl] for lvl in ("family", "genus")}
    species = [(n, c) for lv, n, c in taxa if lv == "species"]
    stats = read_stats(wd / SELECTION_STATS, SELECTION_STATS_KIND) or {}
    missing = _missing(wd)
    genome_files = {p.name for p in list_fasta(wd / "genomes")}
    downloaded = sum(1 for r in ingroup if r.filename in genome_files) if genome_files else None
    reps = sum(1 for r in ingroup if r.gtdb_representative)
    with_quality = [r for r in ingroup if r.completeness is not None]

    report = StepReport(step="selection", title="Genome selection")
    candidates = stats.get("candidates")
    report.key_numbers = [
        KeyNumber("Genomes selected (ingroup)", len(ingroup)),
        KeyNumber("Outgroup genomes", len(outgroup)),
        KeyNumber("Families", len(per_level["family"])),
        KeyNumber("Genera", len(per_level["genus"])),
        KeyNumber("Species", len(species)),
    ]
    if candidates is not None:
        report.key_numbers.insert(0, KeyNumber("Candidate genomes", candidates))
    text = f"The selection holds {_n(len(ingroup), 'ingroup genome')}"
    text += f" of {_n(len(species), 'species', 'species')}"
    text += f" in {_n(len(per_level['genus']), 'genus', 'genera')}"
    if outgroup:
        text += f" and {_n(len(outgroup), 'outgroup genome')}"
    text += "."
    if candidates is not None and stats.get("limit") and candidates > len(ingroup):
        text += (
            f" They were chosen from {candidates:,} candidate genomes by the limit of "
            f"{stats['limit']:,}, ranked by assembly quality within each species."
        )
    report.summary.append(text)

    funnel: list[tuple[str, int | None]] = []
    if candidates is not None:
        funnel.append(("Candidate genomes", candidates))
    funnel.append(("Selected ingroup genomes", len(ingroup)))
    funnel.append(("Outgroup genomes", len(outgroup)))
    if reps:
        funnel.append(("GTDB species representatives among the selected", reps))
    if downloaded is not None:
        funnel.append(("Ingroup genome files present", downloaded))
    if missing is not None:
        funnel.append(("Accessions not obtained", missing))
    report.tables.append(
        write_table(
            out,
            "selection_counts",
            "Genome counts of the selection.",
            ["quantity", "genomes"],
            funnel,
            numeric=["genomes"],
        )
    )
    report.tables.append(
        write_table(
            out,
            "genomes_per_taxon",
            "Ingroup genomes per family, genus and species, largest first.",
            ["level", "name", "genomes"],
            taxa,
            numeric=["genomes"],
        )
    )
    if with_quality:
        rows_q = []
        quantities = (("Completeness (%)", "completeness"), ("Contamination (%)", "contamination"))
        for label, attr in quantities:
            d = describe(getattr(r, attr) for r in with_quality)
            rows_q.append(
                [label, d["n"], *(_round(d[k]) for k in ("min", "q1", "median", "q3", "max"))]
            )
        report.tables.append(
            write_table(
                out,
                "assembly_quality",
                "Estimated completeness and contamination of the selected ingroup genomes "
                "(CheckM values as given by the source).",
                ["quantity", "genomes", "min", "q1", "median", "q3", "max"],
                rows_q,
                numeric=["genomes", "min", "q1", "median", "q3", "max"],
            )
        )

    choice = read_choice(wd)
    if choice:
        _outgroup_choice(report, out, choice)

    if species:
        report.figures.append(_species_figure(out, species, inp.formats))
    if len(with_quality) >= 2:
        report.figures.append(_quality_figure(out, with_quality, outgroup, inp.formats))
    if not with_quality and ingroup:
        report.notes.append("The selection carries no completeness or contamination values.")
    return report


def _outgroup_choice(report: StepReport, out: Path, rows: list[dict[str, str]]) -> None:
    chosen = rows[0]
    mode = chosen.get("mode", "")
    if mode == "sister":
        how = (
            "the genome of best assembly quality in the group nearest the target in the "
            "GTDB reference tree"
        )
        if chosen.get("group") not in ("", "1"):
            how += f" (group {chosen['group']} out, the nearer groups having no genome of "
            how += "sufficient quality)"
    else:
        how = "the genome of best assembly quality among the representatives of the parent taxon"
    report.summary.append(
        f"The outgroup {chosen['accession']} ({chosen['species']}) was chosen as {how}, "
        f"from {_n(len(rows), 'candidate')}."
    )
    columns = [
        "accession",
        "species",
        "chosen",
        "group",
        "tree_distance",
        "completeness",
        "contamination",
        "n50",
        "contigs",
        "isolate",
        "type_material",
    ]
    report.tables.append(
        write_table(
            out,
            "outgroup_candidates",
            "Outgroup candidates, the chosen genome first, then by assembly quality. "
            "group: 1 for the group sister to the target in the GTDB reference tree, 2 "
            "for the next group out, and so on; tree_distance: path length to the "
            "target's common ancestor (substitutions per site of the GTDB marker "
            "alignment).",
            columns,
            ([r.get(c, "") for c in columns] for r in rows),
            numeric=["group", "tree_distance", "completeness", "contamination", "n50", "contigs"],
        )
    )


def _n(count: int, singular: str, plural: str | None = None) -> str:
    """'1 genus', '3 genera'."""
    word = singular if count == 1 else (plural or singular + "s")
    return f"{count:,} {word}"


def _round(value: float | int | None) -> float | None:
    return None if value is None else round(float(value), 2)


def _species_figure(out: Path, species: list[tuple[str, int]], formats) -> FigureRef:
    shown = species[:TOP_SPECIES]
    rest = species[TOP_SPECIES:]
    labels = [name for name, _ in shown]
    values = [n for _, n in shown]
    colors = [style.PRIMARY] * len(shown)
    if rest:
        labels.append(f"other ({len(rest):,} species)")
        values.append(sum(n for _, n in rest))
        colors.append(style.MUTED)

    def draw(fig) -> None:
        ax = fig.add_subplot()
        horizontal_bars(ax, labels, values, colors=colors, xlabel="Genomes")
        for tick in ax.get_yticklabels():
            if not tick.get_text().startswith("other ("):
                tick.set_fontstyle("italic")

    files = render(
        out,
        "genomes_per_species",
        draw,
        width=style.SINGLE_COLUMN,
        height=bar_height(len(labels)),
        formats=formats,
    )
    caption = "Number of selected ingroup genomes per species"
    if rest:
        caption += f"; the {TOP_SPECIES} largest species are shown and the others are combined"
    return FigureRef(
        id="genomes_per_species",
        caption=caption + ".",
        alt="Horizontal bar chart of genomes per species.",
        files=files,
    )


def _quality_figure(
    out: Path, ingroup: list[SelectionRow], outgroup: list[SelectionRow], formats
) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.scatter(
            [r.completeness for r in ingroup],
            [r.contamination or 0.0 for r in ingroup],
            s=10,
            alpha=0.7,
            color=style.PRIMARY,
            linewidths=0,
            label="ingroup",
        )
        og = [r for r in outgroup if r.completeness is not None]
        if og:
            ax.scatter(
                [r.completeness for r in og],
                [r.contamination or 0.0 for r in og],
                s=18,
                marker="D",
                color=style.HIGHLIGHT,
                linewidths=0,
                label="outgroup",
            )
            ax.legend(loc="upper left")
        ax.set_xlabel("Completeness (%)")
        ax.set_ylabel("Contamination (%)")

    files = render(out, "assembly_quality", draw, height=2.6, formats=formats)
    return FigureRef(
        id="assembly_quality",
        caption=(
            "Estimated completeness and contamination of the selected genomes, one point "
            "per genome, as given by the source of the selection."
        ),
        alt="Scatter plot of contamination against completeness.",
        files=files,
    )
