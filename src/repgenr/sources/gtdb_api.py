"""GTDB API as a selection source (``--source api``).

Queries the GTDB API (https://gtdb-api.ecogenomic.org) for the target taxon
only: a much smaller transfer than the table. The API serves its current
release, so ``--release`` and ``--gtdb-version`` do not apply; the date of the
query is recorded in their place.
"""

from __future__ import annotations

import csv
import time
import urllib.parse
from collections.abc import Sequence
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

from ..core import gtdb_tree, http
from ..core.contracts import (
    atomic_path,
    sanitise_taxon_tokens,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.manifest import GenomeRecord
from ..core.outgroup_choice import (
    Candidate,
    Choice,
    choose_by_groups,
    choose_by_quality,
    write_choice,
)
from ..core.selection_summary import (
    GenomeDetails,
)
from .base import (
    GenomeSource,
    LimitCandidate,
    MetadataParams,
    Selection,
    log_choice,
    log_limit,
    record_from_tax,
    refuse_outgroup_in_selection,
    refuse_outgroup_in_target,
    species_or_none,
    stratified_limit,
)

MANIFEST_SOURCE = "gtdb"


def _record_from_tax(accession: str, tax: dict, **fields) -> GenomeRecord:
    return record_from_tax(accession, tax, source=MANIFEST_SOURCE, **fields)


GTDB_API_BASE = "https://gtdb-api.ecogenomic.org"


# Single-letter GTDB rank prefixes used by the API and taxon strings.
_RANK_PREFIX = {"family": "f", "genus": "g", "species": "s"}


# The genomes the GTDB API returned for the target, before --limit: what
# `repgenr census -wd` counts as the candidates of an API-source selection.
GTDB_API_GENOMES = "gtdb_api_genomes.tsv"


_API_GENOMES_COLUMNS = (
    "accession",
    "family",
    "genus",
    "species",
    "gtdb_representative",
    "family_name",
    "genus_name",
    "species_name",
)


def _api_get(path: str, params: dict | None = None) -> dict:
    # Shared retry/backoff session; raises WorkdirError naming the URL on failure.
    return http.get_json(f"{GTDB_API_BASE}{path}", params=params)


def target_taxon(params: MetadataParams) -> str:
    """Build the GTDB taxon string for the selection level (e.g. g__Francisella)."""
    prefix = _RANK_PREFIX[params.level]
    if params.level == "family":
        name = params.target_family
    elif params.level == "genus":
        name = params.target_genus
    else:  # species
        name = f"{params.target_genus} {params.target_species}"
    if name is None:
        raise UserInputError(f"Level '{params.level}' needs a --target-{params.level}.")
    return f"{prefix}__{_capitalize_taxon(name)}"


def _capitalize_taxon(name: str) -> str:
    """Spell a target the way GTDB does: the first word starts upper case and
    a species epithet is lower case ('francisella Tularensis' ->
    'Francisella tularensis').

    Only the first letter of the first word is raised and only the epithet
    before a GTDB suffix is lowered; the rest is kept as typed, so suffixes
    and placeholders ('Bacillus_A', 'CAG-74', 'copri_A') survive.
    """

    def epithet(word: str) -> str:
        stem, sep, suffix = word.partition("_")
        return stem.lower() + sep + suffix

    parts = name.strip().split()
    if not parts:
        return name
    parts[0] = parts[0][:1].upper() + parts[0][1:]
    return " ".join([parts[0], *(epithet(p) for p in parts[1:])])


# Parallel card fetches. Measured 2026-09-05 on 1540 Wolbachia cards: eight in
# flight drew 429 responses that outlived the session's five retries for a
# handful of genomes. Four is gentler; the stragglers get a slow second pass.
_API_CARD_WORKERS = 4


_API_CARD_RETRY_PAUSE = 2.0  # seconds between second-pass requests


def _api_quality(card: dict) -> tuple[float | None, float | None]:
    """CheckM quality from a ``/genome/{gid}/card`` response.

    The card nests it under ``metadata_gene``; the ``genomes-detail`` rows carry
    no quality at all. Prefers checkm2_* and falls back to checkm_*, like the TSV
    path. Returns (None, None) when neither is present.
    """
    gene = card.get("metadata_gene")
    source = gene if isinstance(gene, dict) else card

    def val(key: str) -> float | None:
        raw = source.get(key)
        if raw is None or raw == "":
            return None
        try:
            return float(raw)
        except (TypeError, ValueError):
            return None

    for prefix in ("checkm2", "checkm"):
        completeness = val(f"{prefix}_completeness")
        contamination = val(f"{prefix}_contamination")
        if completeness is not None or contamination is not None:
            return completeness, contamination
    return None, None


def _api_card(accession: str) -> dict:
    return _api_get(f"/genome/{urllib.parse.quote(accession, safe='')}/card")


def _api_cards_by_accession(accessions: Sequence[str], logger) -> dict[str, dict]:
    """Fetch each genome's card; an empty dict for a card that failed.

    Cards are fetched a few at a time; any that still fail after the session's
    own retries (the API rate-limits sustained bursts) are retried once more,
    one at a time with a pause. A card that fails both passes leaves that
    genome unscored (a warning names it) rather than failing the whole
    selection; the keeper then treats it like any other genome without quality.
    """
    if not accessions:
        return {}
    logger.info("Fetching assembly quality for %d genomes from the GTDB API", len(accessions))

    def fetch(acc: str) -> tuple[str, dict | None]:
        try:
            return acc, _api_card(acc)
        except WorkdirError:
            return acc, None

    out: dict[str, dict] = {}
    failed: list[str] = []
    for acc, card in parallel_map(fetch, accessions, _API_CARD_WORKERS):
        if card is None:
            failed.append(acc)
        else:
            out[acc] = card
    if failed:
        logger.info("Retrying %d quality card(s) that the API refused, one at a time", len(failed))
    for acc in failed:
        time.sleep(_API_CARD_RETRY_PAUSE)
        try:
            out[acc] = _api_card(acc)
        except WorkdirError as exc:
            logger.warning("No assembly quality for %s: %s", acc, exc)
            out[acc] = {}
    return out


def _api_quality_by_accession(
    accessions: Sequence[str], logger
) -> dict[str, tuple[float | None, float | None]]:
    """CheckM quality of each genome from its card; see :func:`_api_cards_by_accession`."""
    cards = _api_cards_by_accession(accessions, logger)
    return {acc: _api_quality(card) if card else (None, None) for acc, card in cards.items()}


_GTDB_RANK_FIELDS = (
    "gtdbDomain",
    "gtdbPhylum",
    "gtdbClass",
    "gtdbOrder",
    "gtdbFamily",
    "gtdbGenus",
    "gtdbSpecies",
)


def _api_details(row: dict, card: dict) -> GenomeDetails:
    """The metadata_selected.tsv fields of one genome from the API.

    The GTDB lineage comes from the genomes-detail row (or the card's taxonomy
    block); the NCBI fields from the card, empty when the card failed.
    """
    taxonomy = card.get("metadataTaxonomy") or {}
    ncbi = card.get("metadata_ncbi") or {}
    ranks = [row.get(f) or taxonomy.get(f) or "" for f in _GTDB_RANK_FIELDS]

    def text(block: dict, key: str) -> str:
        value = block.get(key)
        return "" if value is None or str(value).lower() == "none" else str(value)

    return GenomeDetails(
        accession_ncbi=text(ncbi, "ncbi_genbank_assembly_accession"),
        gtdb_taxonomy=";".join(ranks) if any(ranks) else "",
        ncbi_taxonomy=text(taxonomy, "ncbi_taxonomy"),
        ncbi_taxonomy_unfiltered=text(taxonomy, "ncbi_taxonomy_unfiltered"),
        ncbi_strain_identifiers=text(ncbi, "ncbi_strain_identifiers"),
    )


def normalize_api_tax(row: dict) -> dict:
    """Normalize an API row's gtdb* fields to the manifest tax dict.

    Strips the rank prefix and applies the same cleanup as the TSV path:
    species has the genus removed, spaces dropped, underscores -> hyphens.
    """

    def strip(value: str) -> str:
        return value.split("__", 1)[1] if "__" in value else value

    family, genus, species = sanitise_taxon_tokens(
        strip(row.get("gtdbFamily", "")),
        strip(row.get("gtdbGenus", "")),
        strip(row.get("gtdbSpecies", "")),
    )
    return {
        "family": family,
        "genus": genus,
        "species": species,
        "gtdbFamily": row.get("gtdbFamily", ""),
        "gtdbGenus": row.get("gtdbGenus", ""),
    }


def write_api_genomes(path: Path, rows: list[dict]) -> None:
    """Keep the API answer (accession, taxonomy tokens, representative flag)."""
    with atomic_path(path) as tmp, open(tmp, "w", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(_API_GENOMES_COLUMNS)
        for row in rows:
            tax = normalize_api_tax(row)
            writer.writerow(
                [
                    row["gid"],
                    tax["family"],
                    tax["genus"],
                    tax["species"],
                    1 if row.get("gtdbIsRep") else 0,
                    tax["gtdbFamily"].split("__", 1)[-1],
                    tax["gtdbGenus"].split("__", 1)[-1],
                    row.get("gtdbSpecies", "").split("__", 1)[-1],
                ]
            )


def read_api_genomes(path: Path) -> list[dict]:
    """Rows of :func:`write_api_genomes` as dicts with a boolean ``is_rep``."""
    with open(path, encoding="utf-8", newline="") as fo:
        return [
            {
                "accession": row["accession"],
                "tax": {
                    "family": row.get("family") or "",
                    "genus": row.get("genus") or "",
                    "species": row.get("species") or "",
                    "family_name": row.get("family_name") or "",
                    "genus_name": row.get("genus_name") or "",
                    "species_name": row.get("species_name") or "",
                },
                "is_rep": (row.get("gtdb_representative") or "0").strip() in ("1", "true"),
            }
            for row in csv.DictReader(fo, delimiter="\t")
        ]


def _select_via_api(
    params: MetadataParams,
    logger,
    workdir: Path | None = None,
    counts: dict[str, int] | None = None,
) -> tuple[list[GenomeRecord], GenomeRecord | None, dict[str, GenomeDetails]]:
    sp_reps = params.dataset == "rep"
    rows: list[dict] = []
    seen: set[str] = set()
    # One query per species given (-ts a,b), else one for the target taxon.
    for species in species_or_none(params):
        taxon = target_taxon(replace(params, target_species=species))
        logger.info("Querying GTDB API for genomes in %s", taxon)
        found = api_genomes_detail(taxon, sp_reps, logger)
        if not found:
            raise UserInputError(
                f"GTDB API returned no genomes for {taxon}. Check the target name."
            )
        for row in found:
            if row["gid"] not in seen:
                seen.add(row["gid"])
                rows.append(row)
    if workdir is not None:
        write_api_genomes(workdir / GTDB_API_GENOMES, rows)

    # The target's representatives before a --limit cut, to place it in the tree.
    target_reps = [row["gid"] for row in rows if row.get("gtdbIsRep")]
    # Quality for every candidate first: the cut below ranks on it.
    cards = _api_cards_by_accession([row["gid"] for row in rows], logger)
    quality = {acc: _api_quality(card) if card else (None, None) for acc, card in cards.items()}
    if counts is not None:
        counts["candidates"] = len(rows)
    if params.limit and len(rows) > params.limit:
        candidates = [
            LimitCandidate(
                accession=row["gid"],
                species=normalize_api_tax(row)["species"],
                is_rep=bool(row.get("gtdbIsRep", False)),
                completeness=quality[row["gid"]][0],
                contamination=quality[row["gid"]][1],
            )
            for row in rows
        ]
        kept = stratified_limit(candidates, params.limit)
        log_limit(params.limit, candidates, kept, logger)
        keep = {c.accession for c in kept}
        rows = [row for row in rows if row["gid"] in keep]

    records: list[GenomeRecord] = []
    selected_acc: set[str] = set()
    details: dict[str, GenomeDetails] = {}
    for row in rows:
        acc = row["gid"]
        completeness, contamination = quality[acc]
        details[acc] = _api_details(row, cards[acc])
        records.append(
            _record_from_tax(
                acc,
                normalize_api_tax(row),
                completeness=completeness,
                contamination=contamination,
                gtdb_representative=bool(row.get("gtdbIsRep", False)),
            )
        )
        selected_acc.add(acc)

    if params.no_outgroup:
        return records, None, details
    outgroup, outgroup_details = _select_outgroup_via_api(
        params, rows, selected_acc, logger, target_reps=target_reps, workdir=workdir
    )
    details[outgroup.accession] = outgroup_details
    return records, outgroup, details


def api_genomes_detail(taxon: str, sp_reps_only: bool, logger) -> list[dict]:
    encoded = urllib.parse.quote(taxon, safe="")
    try:
        data = _api_get(
            f"/taxon/{encoded}/genomes-detail", {"sp_reps_only": str(sp_reps_only).lower()}
        )
    except http.HTTPStatusError as exc:
        if exc.status != 404:
            raise
        # The API answers 404 for a taxon name it does not know.
        raise UserInputError(
            f"GTDB API has no taxon {taxon}. Check the spelling against GTDB "
            "(names are case-sensitive, e.g. Bacillus_A)."
        ) from exc
    return data.get("rows", [])


def _select_outgroup_via_api(
    params: MetadataParams,
    rows: list[dict],
    selected_acc: set[str],
    logger,
    *,
    target_reps: Sequence[str] = (),
    workdir: Path | None = None,
) -> tuple[GenomeRecord, GenomeDetails]:
    # Explicit outgroup: fetch its card for taxonomy.
    if params.outgroup_accession:
        refuse_outgroup_in_selection(params.outgroup_accession, selected_acc)
        try:
            card = _api_card(params.outgroup_accession)
        except http.HTTPStatusError as exc:
            if exc.status != 404:
                raise
            raise UserInputError(
                f"Outgroup accession {params.outgroup_accession} not found in the GTDB API."
            ) from exc
        tax_row = card.get("metadataTaxonomy", {})
        field = _rank_field(params.level)
        for target in sorted({row.get(field, "") for row in rows}):
            refuse_outgroup_in_target(
                params.outgroup_accession, params.level, target, tax_row.get(field, "")
            )
        completeness, contamination = _api_quality(card)
        record = _record_from_tax(
            params.outgroup_accession,
            normalize_api_tax(tax_row),
            is_outgroup=True,
            completeness=completeness,
            contamination=contamination,
        )
        return record, _api_details(tax_row, card)

    level_field = _rank_field(params.level)
    target_vals = {r.get(level_field) for r in rows}

    def outside(gid: str, tax_row: dict) -> bool:
        return gid not in selected_acc and tax_row.get(level_field) not in target_vals

    def record_for(gid: str, tax_row: dict, card: dict) -> tuple[GenomeRecord, GenomeDetails]:
        completeness, contamination = _api_quality(card) if card else (None, None)
        record = _record_from_tax(
            gid,
            normalize_api_tax(tax_row),
            is_outgroup=True,
            completeness=completeness,
            contamination=contamination,
            gtdb_representative=True,
        )
        return record, _api_details(tax_row, card or {})

    mode = params.outgroup_mode
    choice: Choice | None = None
    taxonomy: dict[str, dict] = {}
    cards: dict[str, dict] = {}
    if mode == "sister" and params.gtdb_tree_path:
        tree = gtdb_tree.read_tree(Path(params.gtdb_tree_path).expanduser())
        placement = gtdb_tree.place(tree, target_reps)
        if placement is None:
            logger.warning(
                "None of the target representatives is in %s; choosing the outgroup "
                "by quality within the parent taxon.",
                params.gtdb_tree_path,
            )
        else:

            def groups():
                for number, ring in enumerate(gtdb_tree.sister_rings(tree, placement), 1):
                    near = ring[:_API_OUTGROUP_CARDS]
                    fetched = _api_cards_by_accession([acc for acc, _ in near], logger)
                    group = []
                    for acc, dist in near:
                        card = fetched.get(acc) or {}
                        tax_row = card.get("metadataTaxonomy") or {}
                        if not card or not outside(acc, tax_row):
                            continue
                        cards[acc], taxonomy[acc] = card, tax_row
                        group.append(_card_candidate(acc, card, tax_row, number, dist))
                    yield group

            choice = choose_by_groups(groups())
    elif mode == "sister":
        logger.warning(
            "--source api: the sister-group choice of the outgroup needs the GTDB "
            "reference tree (--gtdb-tree-path); choosing by quality within the parent taxon."
        )

    if choice is None:
        # The parent taxon (one rank up), outside the selected sub-taxa.
        parent_rank = {"species": "genus", "genus": "family", "family": "class"}[params.level]
        parent_field = {"genus": "gtdbGenus", "family": "gtdbFamily", "class": "gtdbClass"}[
            parent_rank
        ]
        parent_taxon = rows[0].get(parent_field)
        if not parent_taxon:
            raise WorkdirError("Could not determine a parent taxon for outgroup selection.")
        logger.info("Selecting outgroup from parent taxon %s", parent_taxon)
        parent_rows = [
            row
            for row in api_genomes_detail(parent_taxon, sp_reps_only=True, logger=logger)
            if row.get("gtdbIsRep", False) and outside(row["gid"], row)
        ]
        if not parent_rows:
            raise WorkdirError(
                "Could not determine an outgroup via the API; specify --outgroup-accession."
            )
        if mode == "first":
            row = parent_rows[0]
            card = _api_cards_by_accession([row["gid"]], logger).get(row["gid"]) or {}
            return record_for(row["gid"], row, card)
        if len(parent_rows) > _API_OUTGROUP_CARDS:
            logger.info(
                "Scoring the first %d of %d representatives of %s for the outgroup",
                _API_OUTGROUP_CARDS,
                len(parent_rows),
                parent_taxon,
            )
        parent_rows = parent_rows[:_API_OUTGROUP_CARDS]
        fetched = _api_cards_by_accession([row["gid"] for row in parent_rows], logger)
        candidates = []
        for row in parent_rows:
            card = fetched.get(row["gid"]) or {}
            cards[row["gid"]], taxonomy[row["gid"]] = card, row
            candidates.append(_card_candidate(row["gid"], card, row))
        choice = choose_by_quality(candidates)
    if choice is None:
        raise WorkdirError(
            "Could not determine an outgroup via the API; specify --outgroup-accession."
        )
    if workdir is not None:
        write_choice(workdir, choice)
    log_choice(choice, logger)
    gid = choice.chosen.accession
    return record_for(gid, taxonomy[gid], cards[gid])


# Cards fetched per group (or from the parent taxon) to score outgroup candidates.
_API_OUTGROUP_CARDS = 100


def _card_candidate(
    acc: str, card: dict, tax_row: dict, group: int | None = None, distance: float | None = None
) -> Candidate:
    completeness, contamination = _api_quality(card) if card else (None, None)
    nucleotide = card.get("metadata_nucleotide") or {}
    ncbi = card.get("metadata_ncbi") or {}

    def as_int(value) -> int | None:
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    species = tax_row.get("gtdbSpecies", "")
    return Candidate(
        accession=acc,
        species=species.split("__", 1)[-1],
        completeness=completeness,
        contamination=contamination,
        n50=as_int(nucleotide.get("n50_contigs")),
        contigs=as_int(nucleotide.get("contig_count")),
        category=str(ncbi.get("ncbi_genome_category") or ""),
        type_material=str(card.get("gtdbTypeDesignation") or ""),
        group=group,
        distance=distance,
    )


def _rank_field(level: str) -> str:
    return {"family": "gtdbFamily", "genus": "gtdbGenus", "species": "gtdbSpecies"}[level]


class GtdbApiSource(GenomeSource):
    name = "api"
    manifest_source = MANIFEST_SOURCE
    record_tool = "gtdb-api"
    stats_source = "gtdb-api"
    workdir_files = (GTDB_API_GENOMES,)

    def prepare(self, ctx, params, logger) -> dict[str, object]:
        # The API exposes no GTDB release number (its /meta/version is the
        # software version), so the date of the query is recorded in its place.
        api_query_date = datetime.now(UTC).isoformat(timespec="seconds")
        ignored = [
            flag
            for flag, value in (
                ("--release", params.release),
                ("--gtdb-version", params.version),
                ("--metadata-path", params.metadata_path),
                ("--nodownload", params.nodownload),
            )
            if value
        ]
        if ignored:
            logger.warning(
                "--source api serves the current GTDB release; ignoring %s. "
                "Use --source tsv to pin a release or to read a local table.",
                ", ".join(ignored),
            )
        return {"api_query_date": api_query_date}

    def select(self, ctx, params, logger, counts) -> Selection:
        return _select_via_api(params, logger, ctx.workdir, counts)
