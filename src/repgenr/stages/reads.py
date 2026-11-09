"""reads stage: select sequencing runs from ENA/SRA and write ``reads.tsv``.

The first stage of the reads chain (``reads -> assemble -> dereplicate ->
phylo -> tree2tax``). Runs are found by a taxon query (``tax_tree`` over the
resolved taxid) or by explicit run, sample or study accessions, filtered by
platform and size, reduced to the best run per sample, and labelled with the
NCBI family, genus and species of their taxid. Network only; nothing is
downloaded here.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core import ena
from ..core.context import WorkdirContext
from ..core.contracts import READS_TSV, ReadRow, sanitise_taxon_tokens, write_reads
from ..core.errors import UserInputError
from ..core.ncbi_taxonomy import get_taxon_data_from_entrez
from ..core.stage_stats import write_stats

PLATFORMS = ("any", "illumina", "ont", "pacbio")
# Runs remaining after each filter of the selection, for the report.
READS_STATS = "reads_stats.json"
READS_STATS_KIND = "reads"

_ENA_PLATFORM = {"illumina": "ILLUMINA", "ont": "OXFORD_NANOPORE", "pacbio": "PACBIO_SMRT"}
_LONG_READ = {"OXFORD_NANOPORE", "PACBIO_SMRT"}


@dataclass
class ReadsParams:
    target_family: str | None = None
    target_genus: str | None = None
    target_species: str | None = None
    accessions: list[str] = field(default_factory=list)
    accession_file: str | None = None
    platform: str = "any"  # any | illumina | ont | pacbio
    max_runs: int | None = None
    min_bases: int = 0
    # Drop runs above this many bases: a whole-host library (tens of Gb for a
    # small endosymbiont) would assemble into a host-dominated genome.
    max_bases: int | None = None
    # ENA library_selection values to drop (case-insensitive). MDA (multiple
    # displacement amplification) gives chimeric, uneven assemblies.
    drop_selection: list[str] = field(default_factory=lambda: ["MDA"])
    # Keep the best run of each sample: a long-read run with enough bases, else
    # the largest run.
    one_per_sample: bool = True


def validate(params: ReadsParams) -> list[str]:
    """Check the selection before any workdir exists; return the accessions.

    Called by the parameter builder, so a rejected invocation creates nothing,
    and again by :func:`run` for callers that bypass the CLI.
    """
    accessions = [*params.accessions, *read_accession_file(params.accession_file)]
    target = params.target_species or params.target_genus or params.target_family
    if not accessions and not target:
        raise UserInputError(
            "Select runs by taxon (-tf/-tg/-ts) or by accession (--accession, --accession-file)."
        )
    return accessions


def run(ctx: WorkdirContext, params: ReadsParams) -> int:
    logger = ctx.logger
    accessions = validate(params)
    target = params.target_species or params.target_genus or params.target_family

    taxid: str | None = None
    records: list[dict] = []
    if target:
        hit = ena.resolve_taxon(target)
        taxid = hit.taxid
        logger.info(
            "Resolved %r to %s (%s, taxid %s)", target, hit.scientific_name, hit.rank, taxid
        )
        query = ena.taxon_query(taxid)
        ena.announce_search(query, hit.scientific_name, logger)
        records += ena.search_runs(query)
    if accessions:
        records += _whole_genome(ena.search_runs(ena.accession_query(accessions)), logger)
    rows = _dedupe(ena.to_read_rows(records))
    candidates = len(rows)
    logger.info("ENA returned %d whole-genome sequencing runs", candidates)

    # Runs remaining after each filter, for reads_stats.json and the report.
    funnel: list[tuple[str, int]] = [("whole-genome runs returned by ENA", candidates)]
    rows = _filter(rows, params, logger, funnel)
    if params.one_per_sample:
        rows = _best_per_sample(rows)
        funnel.append(("one run per biosample", len(rows)))
    # Ties broken by accession, so --max-runs keeps the same runs whatever
    # order the portal returns them in.
    rows.sort(key=lambda r: (-r.bases, r.run_accession))
    if params.max_runs is not None:
        rows = rows[: params.max_runs]
        funnel.append((f"at most {params.max_runs} runs (--max-runs)", len(rows)))
    if not rows:
        raise UserInputError(
            f"No sequencing runs selected ({candidates} candidates before filtering"
            f"{_active_filters(params)}). Loosen these options or check the taxon."
        )

    rows = _label(rows, logger)
    write_reads(ctx.workdir / READS_TSV, rows)
    write_stats(
        ctx.workdir / READS_STATS,
        READS_STATS_KIND,
        {"candidates": candidates, "selected": len(rows), "funnel": funnel},
    )
    ctx.config.record_stage(
        "reads",
        params={
            **asdict(params),
            "taxid": taxid,
            "candidates": candidates,
            "selected": len(rows),
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Selected %d sequencing runs; wrote %s", len(rows), READS_TSV)
    return len(rows)


def _whole_genome(records: list[dict], logger: logging.Logger) -> list[dict]:
    """The WGS runs of genomic source among runs found by accession.

    The taxon query selects these on the server; a study or sample named by
    accession can also hold RNA-Seq, amplicon or metagenomic runs, which
    would assemble into something other than the organism's genome.
    """
    kept = [r for r in records if ena.is_whole_genome(r)]
    other = [r for r in records if not ena.is_whole_genome(r)]
    if other:
        kinds = sorted(
            {f"{r.get('library_strategy') or '?'}/{r.get('library_source') or '?'}" for r in other}
        )
        names = ", ".join(r["run_accession"] for r in other[:5])
        logger.warning(
            "Dropping %d run(s) found by accession that are not whole-genome sequencing of "
            "genomic DNA (%s): %s%s",
            len(other),
            ", ".join(kinds),
            names,
            " ..." if len(other) > 5 else "",
        )
    return kept


def _active_filters(params: ReadsParams) -> str:
    """Name each filter that was in effect, for the no-match message."""
    parts = []
    if params.platform != "any":
        parts.append(f"--platform {params.platform}")
    if params.min_bases:
        parts.append(f"--min-bases {params.min_bases}")
    if params.max_bases is not None:
        parts.append(f"--max-bases {params.max_bases}")
    if params.drop_selection:
        parts.append(f"--drop-selection {','.join(params.drop_selection)}")
    if params.one_per_sample:
        parts.append("--one-per-sample")
    return f" by {', '.join(parts)}" if parts else ""


def read_accession_file(path: str | None) -> list[str]:
    """Accessions from a file, one per line.

    Text from a ``#`` to the end of the line is a comment, whether the ``#``
    opens the line, follows indentation or follows an accession.
    """
    if path is None:
        return []
    file = Path(path)
    if not file.is_file():
        raise UserInputError(f"--accession-file {path} is not a file.")
    out = []
    for line in file.read_text(encoding="utf-8").splitlines():
        accession = line.split("#", 1)[0].strip()
        if accession:
            out.append(accession)
    return out


def _dedupe(rows: list[ReadRow]) -> list[ReadRow]:
    seen: set[str] = set()
    out = []
    for row in rows:
        if row.run_accession not in seen:
            seen.add(row.run_accession)
            out.append(row)
    return out


def _filter(
    rows: list[ReadRow],
    params: ReadsParams,
    logger: logging.Logger,
    funnel: list[tuple[str, int]] | None = None,
) -> list[ReadRow]:
    """The runs that pass the platform, base-count and library filters;
    ``funnel`` receives the runs remaining after each filter that is set."""
    wanted = _ENA_PLATFORM.get(params.platform)
    kept = [
        r for r in rows if (wanted is None or r.platform == wanted) and r.bases >= params.min_bases
    ]
    if funnel is not None and (wanted is not None or params.min_bases):
        label = []
        if wanted is not None:
            label.append(f"platform {params.platform}")
        if params.min_bases:
            label.append(f"at least {params.min_bases:,} bases")
        funnel.append((" and ".join(label), len(kept)))
    if params.max_bases is not None:
        over = [r for r in kept if r.bases > params.max_bases]
        if over:
            logger.info(
                "Dropping %d run(s) above --max-bases %d (largest %d bp, %s): likely "
                "whole-host libraries",
                len(over),
                params.max_bases,
                max(r.bases for r in over),
                max(over, key=lambda r: r.bases).run_accession,
            )
        kept = [r for r in kept if r.bases <= params.max_bases]
        if funnel is not None:
            funnel.append((f"at most {params.max_bases:,} bases", len(kept)))
    drop = {s.upper() for s in params.drop_selection}
    if drop:
        amplified = [r for r in kept if r.library_selection.upper() in drop]
        if amplified:
            logger.info(
                "Dropping %d run(s) by library selection (%s): %s",
                len(amplified),
                ", ".join(sorted(drop)),
                ", ".join(r.run_accession for r in amplified[:5])
                + (" ..." if len(amplified) > 5 else ""),
            )
        kept = [r for r in kept if r.library_selection.upper() not in drop]
        if funnel is not None:
            funnel.append((f"library selection not {', '.join(sorted(drop))}", len(kept)))
    return kept


# A long-read run is preferred over the sample's short-read runs only when it
# carries enough sequence to assemble on its own: at least this many bases,
# and at least this fraction of the largest short-read run. Below that it is
# usually a scaffolding or test run next to the real data (seen on Wolbachia:
# a 45 kb PacBio run and a 338 Mb ONT run beside 9 Gb and 11 Gb Illumina runs).
LONG_READ_MIN_BASES = 100_000_000
LONG_READ_MIN_FRACTION = 0.1


def _best_per_sample(rows: list[ReadRow]) -> list[ReadRow]:
    by_sample: dict[str, list[ReadRow]] = {}
    for row in rows:
        by_sample.setdefault(row.biosample or row.run_accession, []).append(row)
    return [_best_run(runs) for runs in by_sample.values()]


def _best_run(runs: list[ReadRow]) -> ReadRow:
    short_max = max((r.bases for r in runs if r.platform not in _LONG_READ), default=0)

    def rank(r: ReadRow) -> tuple[bool, int]:
        usable_long = (
            r.platform in _LONG_READ
            and r.bases >= LONG_READ_MIN_BASES
            and r.bases >= LONG_READ_MIN_FRACTION * short_max
        )
        return (usable_long, r.bases)

    return max(runs, key=rank)


def ncbi_taxon_names(taxids: list[str], logger: logging.Logger) -> dict[str, tuple[str, str, str]]:
    """The family, genus and species names of each NCBI taxid, as NCBI writes them.

    One Entrez lineage lookup per distinct taxid; a taxid without a lineage
    maps to empty names. Shared by the reads stage and ``repgenr census
    --runs``.
    """
    unique = sorted({t for t in taxids if t})
    data, missing, _alts = get_taxon_data_from_entrez(unique, logger)
    if missing:
        logger.warning("No NCBI lineage for %d taxid(s); their tokens are 'unknown'", len(missing))
    out: dict[str, tuple[str, str, str]] = {}
    for taxid in unique:
        entry = data.get(taxid) or {}
        taxdata = entry.get("taxdata") or {}
        family, genus, species = (
            (taxdata.get(level) or {}).get("name") or "" for level in ("family", "genus", "species")
        )
        out[taxid] = (family, genus, species)
    return out


def ncbi_taxon_tokens(taxids: list[str], logger: logging.Logger) -> dict[str, tuple[str, str, str]]:
    """:func:`ncbi_taxon_names` as filename tokens (``sanitise_taxon_tokens``)."""
    return {
        t: sanitise_taxon_tokens(*names) for t, names in ncbi_taxon_names(taxids, logger).items()
    }


def _label(rows: list[ReadRow], logger: logging.Logger) -> list[ReadRow]:
    """Fill the family/genus/species tokens from each run's NCBI taxid."""
    tokens = ncbi_taxon_tokens([r.taxid for r in rows], logger)
    labelled = []
    for row in rows:
        family, genus, species = tokens.get(row.taxid, ("", "", ""))
        labelled.append(
            ReadRow(
                **{
                    **asdict(row),
                    "family": family or "unknown",
                    "genus": genus or "unknown",
                    "species": species or "unknown",
                }
            )
        )
    return labelled
