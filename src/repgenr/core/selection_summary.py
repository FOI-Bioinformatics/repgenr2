"""Summary of a genome selection: genome counts per taxon, as a table and plots.

Every entry stage that writes ``selection.tsv`` (metadata, ingest, vgenome,
assemble) calls :func:`write_selection_summary`, so the summary always
describes the current selection. The file names are those of RepGenR 1, whose
metadata command wrote them:

* ``metadata_summary.tsv``: ``level``, ``name``, ``genomes``; one row per
  family, genus and species of the ingroup, largest first.
* ``metadata_summary_number_in_level.png``: genomes per taxon, one bar per
  family, genus and species.
* ``metadata_summary_number_per_level.png``: distinct taxa per level.

The outgroup is not counted. ``selection_stats.json`` holds the counts of
the selection that ``selection.tsv`` does not (the candidates before a
``--limit`` cut); it is rewritten with every selection, empty when the stage
reports none. ``metadata_selected.tsv`` (written by the
metadata stage with the GTDB and NCBI lineages) is trimmed here to the
accessions still selected, so another entry stage never leaves rows of
genomes that left the selection.
"""

from __future__ import annotations

import csv
import logging
from collections import Counter
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .contracts import SelectionRow, atomic_path, atomic_replace
from .stage_stats import write_stats

SUMMARY_TSV = "metadata_summary.tsv"
PLOT_IN_LEVEL = "metadata_summary_number_in_level.png"
PLOT_PER_LEVEL = "metadata_summary_number_per_level.png"
METADATA_SELECTED_TSV = "metadata_selected.tsv"
SELECTION_STATS = "selection_stats.json"
SELECTION_STATS_KIND = "selection"

LEVELS = ("family", "genus", "species")
UNCLASSIFIED = "unclassified"
# Bars drawn per level in the genomes-per-taxon plot; the table lists all.
MAX_BARS_PER_LEVEL = 60

METADATA_SELECTED_COLUMNS = (
    "accession",
    "accession_ncbi",
    "is_outgroup",
    "source_representative",
    "gtdb_taxonomy",
    "ncbi_taxonomy",
    "ncbi_taxonomy_unfiltered",
    "ncbi_strain_identifiers",
)


@dataclass(frozen=True)
class GenomeDetails:
    """Per-genome metadata from GTDB that ``selection.tsv`` does not carry."""

    accession_ncbi: str = ""
    gtdb_taxonomy: str = ""
    ncbi_taxonomy: str = ""
    ncbi_taxonomy_unfiltered: str = ""
    ncbi_strain_identifiers: str = ""


# True: every genome's species token is the whole species name; a collection:
# the accessions whose token is (an ingest that combines viral and other sets).
FullSpecies = bool | Collection[str]


def _species_name(row: SelectionRow, full_species: FullSpecies) -> str:
    if not row.species:
        return ""
    whole = full_species if isinstance(full_species, bool) else row.accession in full_species
    if whole or not row.genus:
        return row.species
    return f"{row.genus} {row.species}"


def count_taxa(
    rows: Iterable[SelectionRow], *, full_species: FullSpecies = False
) -> list[tuple[str, str, int]]:
    """``(level, name, genomes)`` for the ingroup, per level largest first,
    then by name. An empty token is counted as ``unclassified``.

    A bacterial species token is the epithet, named here with its genus
    ("Yersinia pestis"). A viral token is the whole species name as in the
    file name ("Hepatovirus-ahepa") and is taken unchanged for the genomes
    ``full_species`` names: all of them when True, else those accessions."""
    counters: dict[str, Counter[str]] = {level: Counter() for level in LEVELS}
    for row in rows:
        if row.is_outgroup:
            continue
        names = {
            "family": row.family,
            "genus": row.genus,
            "species": _species_name(row, full_species),
        }
        for level in LEVELS:
            counters[level][names[level] or UNCLASSIFIED] += 1
    out: list[tuple[str, str, int]] = []
    for level in LEVELS:
        for name, n in sorted(counters[level].items(), key=lambda kv: (-kv[1], kv[0])):
            out.append((level, name, n))
    return out


def write_selection_summary(
    workdir: Path,
    rows: Sequence[SelectionRow],
    logger: logging.Logger,
    *,
    full_species: FullSpecies = False,
    stats: Mapping[str, Any] | None = None,
) -> None:
    """Write the summary table and plots of ``rows``; remove them when the
    selection holds no ingroup genome. Trim ``metadata_selected.tsv``.
    ``full_species`` as in :func:`count_taxa` (the viral stages). ``stats``
    (for example ``candidates``) goes to ``selection_stats.json``."""
    _trim_metadata_selected(workdir, {r.accession for r in rows})
    write_stats(workdir / SELECTION_STATS, SELECTION_STATS_KIND, stats or {})
    counts = count_taxa(rows, full_species=full_species)
    paths = [workdir / name for name in (SUMMARY_TSV, PLOT_IN_LEVEL, PLOT_PER_LEVEL)]
    if not counts:
        for path in paths:
            path.unlink(missing_ok=True)
        return
    with atomic_replace(workdir / SUMMARY_TSV, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(["level", "name", "genomes"])
        writer.writerows(counts)
    # The plots are a report beside the selection: a failure to draw them is
    # logged and does not fail the stage that selected the genomes.
    try:
        _plot_in_level(counts, workdir / PLOT_IN_LEVEL)
        _plot_per_level(counts, workdir / PLOT_PER_LEVEL)
    except Exception as exc:  # noqa: BLE001 - see comment above
        for path in paths[1:]:
            path.unlink(missing_ok=True)
        logger.warning("Could not draw the selection summary plots: %s", exc)


def _plot_in_level(counts: list[tuple[str, str, int]], path: Path) -> None:
    from matplotlib.figure import Figure

    names: list[str] = []
    shown: dict[str, int] = {}
    total: dict[str, int] = {}
    fig = Figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot()
    x = 0
    for level in LEVELS:
        level_counts = [(name, n) for lvl, name, n in counts if lvl == level]
        total[level] = len(level_counts)
        kept = level_counts[:MAX_BARS_PER_LEVEL]
        shown[level] = len(kept)
        ax.bar(range(x, x + len(kept)), [n for _, n in kept], label=level)
        names.extend(name for name, _ in kept)
        x += len(kept)
    # Widen the figure with the number of bars, so labels stay legible.
    fig.set_size_inches(max(6.4, 0.22 * len(names) + 2.0), 4.8 + 0.06 * max(map(len, names)))
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90)
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylabel("Genomes")
    ax.legend(loc="upper right")
    title = "Number of genomes in each taxon"
    cut = [f"{shown[lvl]} of {total[lvl]} {lvl}" for lvl in LEVELS if shown[lvl] < total[lvl]]
    if cut:
        title += f" (largest shown: {', '.join(cut)})"
    ax.set_title(title)
    fig.tight_layout()
    _save(fig, path)


def _plot_per_level(counts: list[tuple[str, str, int]], path: Path) -> None:
    from matplotlib.figure import Figure

    per_level = Counter(level for level, _, _ in counts)
    values = [per_level[level] for level in LEVELS]
    fig = Figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot()
    bars = ax.bar(LEVELS, values)
    ax.bar_label(bars)
    ax.set_ylabel("Taxa")
    ax.set_title("Number of taxa per taxonomic level")
    fig.tight_layout()
    _save(fig, path)


def _save(fig, path: Path) -> None:
    # The temporary name ends in .tmp, so the format is given explicitly.
    with atomic_path(path) as tmp:
        fig.savefig(tmp, format="png", dpi=200)


def write_metadata_selected(
    workdir: Path, rows: Sequence[SelectionRow], details: dict[str, GenomeDetails]
) -> None:
    """Write ``metadata_selected.tsv``: one row per selected genome, outgroup
    included, with the GTDB and NCBI lineages and the NCBI strain."""
    with atomic_replace(workdir / METADATA_SELECTED_TSV, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(METADATA_SELECTED_COLUMNS)
        for row in rows:
            d = details.get(row.accession, GenomeDetails())
            writer.writerow(
                [
                    row.accession,
                    d.accession_ncbi,
                    "1" if row.is_outgroup else "0",
                    "1" if row.source_representative else "0",
                    d.gtdb_taxonomy,
                    d.ncbi_taxonomy,
                    d.ncbi_taxonomy_unfiltered,
                    d.ncbi_strain_identifiers,
                ]
            )


def _trim_metadata_selected(workdir: Path, accessions: set[str]) -> None:
    path = workdir / METADATA_SELECTED_TSV
    if not path.exists():
        return
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        header = next(reader, None)
        body = list(reader)
    kept = [r for r in body if r and r[0] in accessions]
    if header is None or not kept:
        path.unlink()
        return
    if len(kept) == len(body):
        return
    with atomic_replace(path, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(kept)
