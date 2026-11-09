"""metadata stage: select a taxon's genomes from GTDB.

Two data sources:

* ``tsv`` (default): download the full GTDB metadata table and parse it. Robust
  and release-pinned, but downloads the whole table.
* ``api``: query the GTDB API (https://gtdb-api.ecogenomic.org), fetching only
  the target taxon's genomes. Much smaller transfer; uses the API's current
  release (``--release``/``--version`` are advisory for this source).

Either way the selection is recorded in the SQLite manifest plus ``repgenr.yaml``
provenance.
"""

from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
import tarfile
import time
import urllib.parse
from collections.abc import Sequence
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

from ..core import http
from ..core.context import WorkdirContext
from ..core.contracts import (
    SELECTION_TSV,
    SelectionRow,
    atomic_path,
    genome_filename,
    sanitise_taxon_tokens,
    write_selection,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.integrity import refuse_foreign_rows
from ..core.manifest import GenomeRecord
from ..core.selection_summary import (
    GenomeDetails,
    write_metadata_selected,
    write_selection_summary,
)

TAXONOMY = ("domain", "phylum", "class", "family", "genus", "species")
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


@dataclass
class MetadataParams:
    dataset: str  # all | rep
    level: str  # family | genus | species
    release: str | None = None  # required for tsv source
    version: str | None = None  # bac120 | ar53; required for tsv source
    source: str = "tsv"  # tsv | api
    target_family: str | None = None
    target_genus: str | None = None
    # One species epithet, or several separated by commas (all of the target
    # genus); see target_species_list.
    target_species: str | None = None
    outgroup_accession: str | None = None
    # Select no outgroup: selection.tsv holds the ingroup only and no
    # outgroup is downloaded (root the tree with phylo --root midpoint).
    no_outgroup: bool = False
    metadata_path: str | None = None
    nodownload: bool = False
    limit: int | None = None
    # Discard genomes another entry path appended (assemble --append) instead
    # of refusing to overwrite the selection that holds them.
    drop_foreign: bool = False


def run(ctx: WorkdirContext, params: MetadataParams) -> int:
    logger = ctx.logger
    refuse_foreign_rows(ctx, "metadata", drop_foreign=params.drop_foreign, logger=logger)
    _validate(params)

    api_query_date: str | None = None
    # Candidate genomes before a --limit cut, for selection_stats.json.
    counts: dict[str, int] = {}
    if params.source == "api":
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
        selected, outgroup, details = _select_via_api(params, logger, ctx.workdir, counts)
    else:
        # A table run leaves no API answer behind for census to read.
        (ctx.workdir / GTDB_API_GENOMES).unlink(missing_ok=True)
        selected, outgroup, details = _select_via_tsv(ctx, params, logger, counts)

    logger.info(
        "Selected %d genomes; outgroup: %s",
        len(selected),
        outgroup.accession if outgroup is not None else "none (--no-outgroup)",
    )
    _populate_manifest(ctx, selected, outgroup)
    _write_outgroup_file(ctx, outgroup.accession if outgroup is not None else None)
    rows = _write_selection(ctx, selected, outgroup)
    write_metadata_selected(ctx.workdir, rows, details)
    write_selection_summary(
        ctx.workdir,
        rows,
        logger,
        stats={
            "source": f"gtdb-{params.source}",
            "candidates": counts.get("candidates", len(selected)),
            "limit": params.limit,
        },
    )

    ctx.config.record_stage(
        "metadata",
        # The GTDB web API or a GTDB metadata table; no external binary.
        tool="gtdb-api" if params.source == "api" else "gtdb-table",
        params={
            "source": params.source,
            "release": params.release,
            "api_query_date": api_query_date,
            "version": params.version,
            "dataset": params.dataset,
            "level": params.level,
            "target_family": params.target_family,
            "target_genus": params.target_genus,
            "target_species": params.target_species,
            "selected_count": len(selected),
            "outgroup": outgroup.accession if outgroup is not None else None,
            "no_outgroup": params.no_outgroup,
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(selected)


def _species_or_none(params: MetadataParams) -> list[str | None]:
    """The species to select one by one; [None] for a genus or family target."""
    return list(target_species_list(params.target_species)) or [None]


def target_species_list(value: str | None) -> list[str]:
    """The species epithets of a -ts value: one, or several separated by
    commas ("tularensis,novicida"), in the order given, without repeats."""
    if not value:
        return []
    out: list[str] = []
    for part in value.split(","):
        name = part.strip()
        if name and name not in out:
            out.append(name)
    return out


def gtdb_provenance(params: dict) -> dict[str, str]:
    """The GTDB reference behind a metadata record, for ``versions`` and ``status``.

    The table path names its release (``gtdb_release``); the API path has no
    release number and names the UTC date of its query (``gtdb_api_query_date``).
    A record from before the query date was recorded yields nothing for the API.
    """
    if params.get("source") == "api":
        queried = params.get("api_query_date")
        return {"gtdb_api_query_date": str(queried)} if queried else {}
    release = params.get("release")
    return {"gtdb_release": str(release)} if release else {}


GTDB_VERSIONS = ("bac120", "ar53")
_RELEASE_RE = re.compile(r"^\d+\.\d+$")


def _validate(params: MetadataParams) -> None:
    if params.source not in ("tsv", "api"):
        raise UserInputError("--source must be 'tsv' or 'api'.")
    if params.source == "tsv":
        if not params.release or not _RELEASE_RE.match(params.release):
            raise UserInputError("tsv source needs --release like '232.0' (major.minor).")
        if not params.version:
            raise UserInputError("tsv source needs --gtdb-version (bac120 or ar53).")
        if params.version not in GTDB_VERSIONS:
            raise UserInputError(
                f"Invalid --gtdb-version '{params.version}'. Choose from: "
                f"{', '.join(GTDB_VERSIONS)}."
            )
        if params.metadata_path and not Path(params.metadata_path).is_file():
            raise UserInputError(f"--metadata-path not found: {params.metadata_path}")
    if not (params.target_genus or params.target_family):
        raise UserInputError("Supply --target-genus or --target-family.")
    if params.level == "species" and not params.target_species:
        raise UserInputError("Level 'species' needs --target-species.")
    species = target_species_list(params.target_species)
    if len(species) > 1 and params.level != "species":
        raise UserInputError(
            f"{len(species)} species were given with -ts; select them with --level species "
            "(at --level genus the whole genus is selected)."
        )
    if len(species) > 1 and not params.target_genus:
        raise UserInputError("Several species (-ts) need their genus (-tg).")
    if params.no_outgroup and params.outgroup_accession:
        raise UserInputError("Give either --outgroup-accession or --no-outgroup, not both.")
    if params.level == "genus" and not params.target_genus:
        raise UserInputError("Level 'genus' needs --target-genus.")
    if params.level == "family" and not params.target_family:
        raise UserInputError("Level 'family' needs --target-family.")


def _select_via_tsv(
    ctx: WorkdirContext, params: MetadataParams, logger, counts: dict[str, int] | None = None
) -> tuple[list[GenomeRecord], GenomeRecord | None, dict[str, GenomeDetails]]:
    metadata_file = _obtain_metadata(ctx, params, logger)
    accessions = _parse_metadata(metadata_file, params, logger)

    targets = []
    for species in _species_or_none(params):
        target_levels = _target_levels(accessions, replace(params, target_species=species))
        if not target_levels:
            raise UserInputError(
                f"Target not found in GTDB: family={params.target_family} "
                f"genus={params.target_genus} species={species} at level {params.level}"
            )
        targets.append(target_levels)
    selected = _select(accessions, targets, params.limit, logger, counts)
    records = [
        _record_from_tax(
            acc,
            data["tax"],
            completeness=data.get("completeness"),
            contamination=data.get("contamination"),
            gtdb_representative=bool(data.get("is_rep")),
        )
        for acc, data in selected.items()
    ]
    details = {acc: data["details"] for acc, data in selected.items()}
    if params.no_outgroup:
        return records, None, details
    outgroup_acc, outgroup_data = _pick_outgroup(accessions, selected, targets, params)

    outgroup = _record_from_tax(
        outgroup_acc,
        outgroup_data["tax"],
        is_outgroup=True,
        completeness=outgroup_data.get("completeness"),
        contamination=outgroup_data.get("contamination"),
        gtdb_representative=bool(outgroup_data.get("is_rep")),
    )
    details[outgroup_acc] = outgroup_data["details"]
    return records, outgroup, details


def _record_from_tax(
    accession: str,
    tax: dict,
    is_outgroup: bool = False,
    completeness: float | None = None,
    contamination: float | None = None,
    gtdb_representative: bool = False,
) -> GenomeRecord:
    return GenomeRecord(
        accession=accession,
        source="gtdb",
        is_outgroup=is_outgroup,
        family=tax["family"],
        genus=tax["genus"],
        species=tax["species"],
        completeness=completeness,
        contamination=contamination,
        gtdb_representative=gtdb_representative,
    )


# Newer GTDB releases (>= r220) ship the metadata as a plain ``.tsv.gz``;
# older ones (<= r214) as a ``.tar.gz``. The modern layout is tried first.
_TABLE_SUFFIXES = (".tsv.gz", ".tar.gz")


def workdir_tables(workdir: Path, release: str, version: str) -> list[Path]:
    """The GTDB table paths a download writes and ``--nodownload`` reuses.

    One per naming scheme, in the order they are tried. The resume fingerprint
    declares them as inputs under ``--nodownload``, so replacing a reused table
    reruns the stage.
    """
    major = int(float(release))
    return [workdir / f"{version}_metadata_r{major}{ext}" for ext in _TABLE_SUFFIXES]


def release_marker(workdir: Path, release: str, version: str) -> Path:
    """The file naming the exact release (e.g. ``232.0``) the workdir table came from.

    Table names carry only the major release, so ``--nodownload`` reads this
    marker to avoid serving a 214.0 table as 214.1.
    """
    major = int(float(release))
    return workdir / f"{version}_metadata_r{major}.release"


def _check_reused_release(workdir: Path, params: MetadataParams, table: Path, logger) -> None:
    marker = release_marker(workdir, str(params.release), str(params.version))
    try:
        recorded = marker.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        logger.warning(
            "Cannot confirm which GTDB release %s came from (no %s); using it as %s.",
            table.name,
            marker.name,
            params.release,
        )
        return
    if recorded != params.release:
        raise UserInputError(
            f"{table.name} in the workdir was downloaded for GTDB release {recorded}, "
            f"not {params.release}. Pass -r {recorded}, or drop --nodownload to download "
            f"release {params.release}."
        )


def _obtain_metadata(ctx: WorkdirContext, params: MetadataParams, logger) -> Path:
    return obtain_metadata_table(ctx.workdir, params, logger)


def obtain_metadata_table(workdir: Path, params: MetadataParams, logger) -> Path:
    """The GTDB metadata table for ``params``: ``--metadata-path``, a table
    reused from ``workdir`` under ``--nodownload``, or a fresh download into
    ``workdir`` (checked against the release's MD5SUM.txt).

    The metadata stage passes its working directory; ``repgenr census``
    passes a cache directory, so both share one download and reuse rule.
    """
    if params.metadata_path:
        logger.info("Using provided metadata: %s", params.metadata_path)
        return Path(params.metadata_path)

    workdir.mkdir(parents=True, exist_ok=True)
    if params.release is None:
        raise UserInputError("tsv source needs --release like '232.0' (major.minor).")
    major = int(float(params.release))
    base = (
        f"https://data.gtdb.ecogenomic.org/releases/release{major}/"
        f"{params.release}/{params.version}_metadata_r{major}"
    )
    tables = workdir_tables(workdir, params.release, str(params.version))
    if params.nodownload:
        for dest in tables:
            if dest.exists():
                _check_reused_release(workdir, params, dest, logger)
                logger.info("Using previously downloaded %s", dest.name)
                return dest
    # Try the modern layout first, then fall back, so current releases
    # resolve on the first request. Only a 404 moves on to the next naming
    # scheme; any other failure (network, checksum) is reported as is.
    for ext, dest in zip(_TABLE_SUFFIXES, tables, strict=True):
        url = base + ext
        logger.info("Downloading %s", url)
        try:
            # Streams through a .part file with a size check (core.http), and
            # retries transient 5xx/429.
            http.download(url, dest, logger=logger)
        except http.HTTPStatusError as exc:
            if exc.status != 404:
                raise
            logger.info("Not found: %s", url)
            continue
        try:
            http.verify_md5_manifest(dest, base.rsplit("/", 1)[0] + "/MD5SUM.txt", logger=logger)
        except WorkdirError:
            # corrupt transfer: remove it so a re-run downloads afresh
            dest.unlink(missing_ok=True)
            raise
        release_marker(workdir, params.release, str(params.version)).write_text(
            params.release + "\n", encoding="utf-8"
        )
        return dest
    raise UserInputError(
        f"GTDB has no {params.version} metadata table for release {params.release} "
        f"({base}{' or '.join(_TABLE_SUFFIXES)} not found); check -r/--release and "
        "--gtdb-version."
    )


def _open_metadata(path: Path, workdir: Path):
    if path.name.endswith(".tar.gz"):
        tsv_gz = workdir / path.name.replace(".tar.gz", ".tsv.gz")
        # Stream the member straight out of the tarball into a gzipped TSV.
        # Avoid extract-to-directory + rmtree, which fails on exFAT/NTFS volumes
        # (shutil.rmtree's dir_fd traversal is unsupported there).
        with tarfile.open(path) as tar:
            member = next((m for m in tar.getmembers() if m.name.endswith(".tsv")), None)
            if member is None:
                raise WorkdirError("No .tsv inside GTDB tarball")
            source = tar.extractfile(member)
            if source is None:
                raise WorkdirError("Could not read .tsv from GTDB tarball")
            with atomic_path(tsv_gz) as tmp, source, gzip.open(tmp, "wb") as fo:
                shutil.copyfileobj(source, fo)
        path = tsv_gz
    return gzip.open(path, "rt")


def _parse_metadata(path: Path, params: MetadataParams, logger) -> dict[str, dict]:
    logger.info("Parsing GTDB metadata")
    accessions: dict[str, dict] = {}
    with _open_metadata(path, path.parent) as fo:
        reader = csv.reader(fo, delimiter="\t")
        header = next(reader, [])
        idx = {name: i for i, name in enumerate(header)}
        for fields in reader:
            if not fields:
                continue
            acc_raw = fields[idx["accession"]]
            accession = acc_raw.replace("GB_", "").replace("RS_", "")
            rep = fields[idx["gtdb_genome_representative"]]
            is_rep = acc_raw == rep

            # Under --dataset rep a non-representative row is kept only when it
            # is the named outgroup, and then only to resolve that outgroup:
            # it never joins the selection.
            outgroup_only = params.dataset == "rep" and not is_rep
            if outgroup_only and accession != params.outgroup_accession:
                continue

            tax = _parse_taxonomy(fields[idx["gtdb_taxonomy"]])
            completeness = _opt_column(fields, idx, "checkm2_completeness", "checkm_completeness")
            contamination = _opt_column(
                fields, idx, "checkm2_contamination", "checkm_contamination"
            )
            accession_ncbi = fields[idx.get("ncbi_genbank_assembly_accession", 0)]
            accessions[accession] = {
                "accession": accession,
                "accession_ncbi": accession_ncbi,
                "details": GenomeDetails(
                    accession_ncbi=accession_ncbi,
                    # Lineage strings repeat across genomes; interning keeps
                    # one copy each when the whole table is parsed.
                    gtdb_taxonomy=sys.intern(fields[idx["gtdb_taxonomy"]]),
                    ncbi_taxonomy=sys.intern(_text_column(fields, idx, "ncbi_taxonomy")),
                    ncbi_taxonomy_unfiltered=sys.intern(
                        _text_column(fields, idx, "ncbi_taxonomy_unfiltered")
                    ),
                    ncbi_strain_identifiers=_text_column(fields, idx, "ncbi_strain_identifiers"),
                ),
                "tax": tax,
                "is_rep": is_rep,
                "outgroup_only": outgroup_only,
                "completeness": completeness,
                "contamination": contamination,
            }
    logger.info("Parsed %d accessions", len(accessions))
    return accessions


def _text_column(fields: list[str], idx: dict[str, int], name: str) -> str:
    """The column ``name`` of a row; empty when absent or ``none``."""
    i = idx.get(name)
    if i is None or i >= len(fields):
        return ""
    raw = fields[i].strip()
    return "" if raw.lower() == "none" else raw


def _opt_column(fields: list[str], idx: dict[str, int], *names: str) -> float | None:
    """First present, non-empty, numeric column among ``names``; else None."""
    for name in names:
        i = idx.get(name)
        if i is None or i >= len(fields):
            continue
        raw = fields[i].strip()
        if not raw or raw.lower() in {"none", "na", "n/a"}:
            continue
        try:
            return float(raw)
        except ValueError:
            continue
    return None


def _parse_taxonomy(raw: str) -> dict[str, str]:
    tax = {level: "" for level in TAXONOMY}
    for chunk in raw.split(";"):
        for level in TAXONOMY:
            key = level[0] + "__"
            if chunk.startswith(key):
                tax[level] = chunk[len(key) :]
    # The names as GTDB writes them ('Bacillus_A', 'Bacillus_A cereus_S'),
    # for display; the levels below become filename tokens.
    tax["family_name"], tax["genus_name"], tax["species_name"] = (
        tax["family"],
        tax["genus"],
        tax["species"],
    )
    tax["family"], tax["genus"], tax["species"] = sanitise_taxon_tokens(
        tax["family"], tax["genus"], tax["species"]
    )
    return tax


def _target_levels(accessions: dict[str, dict], params: MetadataParams) -> dict[str, str]:
    def norm(value: str | None) -> str | None:
        return value.lower().replace("_", "-") if value else None

    tf, tg, ts = norm(params.target_family), norm(params.target_genus), norm(params.target_species)
    for data in accessions.values():
        tax = data["tax"]
        matched = False
        if ts and tg:
            matched = tax["species"].lower() == ts and tax["genus"].lower() == tg
        elif tg:
            matched = tax["genus"].lower() == tg
        elif tf:
            matched = tax["family"].lower() == tf
        if matched:
            levels: dict[str, str] = {}
            for level in TAXONOMY:
                levels[level] = tax[level]
                if level == params.level:
                    break
            return levels
    return {}


@dataclass(frozen=True)
class _Candidate:
    """One genome competing for a place under ``--limit``."""

    accession: str
    species: str
    is_rep: bool
    completeness: float | None
    contamination: float | None

    def score(self) -> float | None:
        if self.completeness is None or self.contamination is None:
            return None
        from .derep_keeper import quality_score

        return quality_score(self.completeness, self.contamination)


def _stratified_limit(candidates: Sequence[_Candidate], limit: int | None) -> list[_Candidate]:
    """Choose up to ``limit`` genomes: the best of every species first, then each
    species' next best, and so on (round-robin), so a heavily sequenced species
    cannot fill the cap and file order plays no part.

    Within a species genomes rank by the CheckM part of the keeper score (no
    N50 exists before download), unscored genomes last, then the GTDB
    species-representative flag, then accession, so the choice is
    deterministic. Species are visited in name order. With no
    limit every candidate is returned in that ranked order.
    """

    def rank(c: _Candidate) -> tuple:
        score = c.score()
        return (score is None, -(score or 0.0), not c.is_rep, c.accession)

    by_species: dict[str, list[_Candidate]] = {}
    for c in candidates:
        by_species.setdefault(c.species, []).append(c)
    queues = [sorted(members, key=rank) for _species, members in sorted(by_species.items())]

    kept: list[_Candidate] = []
    total = len(candidates)
    target = total if limit is None else min(limit, total)
    depth = 0
    while len(kept) < target:
        for queue in queues:
            if depth < len(queue):
                kept.append(queue[depth])
                if len(kept) == target:
                    break
        depth += 1
    return kept


def _log_limit(limit: int | None, candidates: Sequence[_Candidate], kept: Sequence, logger) -> None:
    if limit and len(kept) < len(candidates):
        species = len({c.species for c in candidates})
        logger.info(
            "Limit %d: kept %d of %d candidate genomes across %d species "
            "(best assembly quality per species first)",
            limit,
            len(kept),
            len(candidates),
            species,
        )


def _select(
    accessions: dict[str, dict],
    targets: Sequence[dict[str, str]],
    limit: int | None,
    logger,
    counts: dict[str, int] | None = None,
):
    """Every genome matching one of the targets (one per species given), cut
    to ``limit`` by species-stratified quality ranking (never by file order).
    ``counts["candidates"]`` is set to the number of matching genomes before
    the cut."""
    if isinstance(targets, dict):  # one target, as callers before several species
        targets = [targets]
    matching = {
        acc: data
        for acc, data in accessions.items()
        if not data.get("outgroup_only")
        and any(all(data["tax"][lvl] == val for lvl, val in t.items()) for t in targets)
    }
    if counts is not None:
        counts["candidates"] = len(matching)
    if not limit or len(matching) <= limit:
        return matching
    candidates = [
        _Candidate(
            accession=acc,
            species=data["tax"]["species"],
            is_rep=bool(data.get("is_rep")),
            completeness=data.get("completeness"),
            contamination=data.get("contamination"),
        )
        for acc, data in matching.items()
    ]
    kept = _stratified_limit(candidates, limit)
    _log_limit(limit, candidates, kept, logger)
    return {c.accession: matching[c.accession] for c in kept}


def _pick_outgroup(accessions, selected, targets, params):
    if isinstance(targets, dict):
        targets = [targets]
    target_vals = {t[params.level] for t in targets}
    target_levels = targets[0]
    if params.outgroup_accession:
        if params.outgroup_accession not in accessions:
            raise UserInputError(
                f"Outgroup accession {params.outgroup_accession} not in GTDB metadata."
            )
        _refuse_outgroup_in_selection(params.outgroup_accession, selected)
        data = accessions[params.outgroup_accession]
        for target in sorted(target_vals):
            _refuse_outgroup_in_target(
                params.outgroup_accession, params.level, target, data["tax"][params.level]
            )
        return params.outgroup_accession, data

    # A representative one level above the selection level, outside the target
    # taxon itself: under --limit, genomes of the target that the cap left out
    # are neither selected nor an outgroup.
    levels = list(target_levels)
    upper = levels[-2] if len(levels) >= 2 else levels[-1]
    upper_val = next(iter(selected.values()))["tax"][upper]
    for acc, data in accessions.items():
        if acc in selected or data["tax"][params.level] in target_vals:
            continue
        if data["tax"][upper] == upper_val and data["is_rep"]:
            return acc, data
    raise WorkdirError(
        "Could not determine an outgroup; specify --outgroup-accession, or select none "
        "with --no-outgroup."
    )


def _refuse_outgroup_in_selection(accession: str, selected) -> None:
    if accession in selected:
        raise UserInputError(
            f"--outgroup-accession {accession} is part of the selection itself; "
            "an outgroup must lie outside it."
        )


def _refuse_outgroup_in_target(accession: str, level: str, target: str, value: str) -> None:
    """A named outgroup obeys the automatic rule: it lies outside the target taxon.

    Under --dataset rep or --limit a target genome can be left out of the
    selection; it is still not an outgroup.
    """
    if value and value == target:
        raise UserInputError(
            f"--outgroup-accession {accession} lies inside the target {level} {target}; "
            "an outgroup must lie outside the target taxon."
        )


def _populate_manifest(ctx, selected: list[GenomeRecord], outgroup: GenomeRecord | None) -> None:
    # Replace, not upsert: a re-selection must remove de-selected rows.
    ctx.manifest.replace_genomes([*selected, *([outgroup] if outgroup is not None else [])])


def _write_outgroup_file(ctx, outgroup_acc: str | None) -> None:
    """Name the outgroup; without one (--no-outgroup) remove the file, as
    ingest does, so no earlier outgroup is named."""
    path = ctx.workdir / "outgroup_accession.txt"
    if outgroup_acc is None:
        path.unlink(missing_ok=True)
    else:
        path.write_text(outgroup_acc + "\n")


def _write_selection(
    ctx, selected: list[GenomeRecord], outgroup: GenomeRecord | None
) -> list[SelectionRow]:
    """Publish the portable selection.tsv (the metadata -> genome data-channel hand-off)."""
    rows = []
    for r in (*selected, *([outgroup] if outgroup is not None else [])):
        family, genus, species = r.family or "", r.genus or "", r.species or ""
        rows.append(
            SelectionRow(
                accession=r.accession,
                family=family,
                genus=genus,
                species=species,
                is_outgroup=r.is_outgroup,
                filename=genome_filename(family, genus, species, r.accession),
                completeness=r.completeness,
                contamination=r.contamination,
                gtdb_representative=r.gtdb_representative,
            )
        )
    write_selection(ctx.workdir / SELECTION_TSV, rows)
    return rows


# --- GTDB API source --------------------------------------------------------


def _api_get(path: str, params: dict | None = None) -> dict:
    # Shared retry/backoff session; raises WorkdirError naming the URL on failure.
    return http.get_json(f"{GTDB_API_BASE}{path}", params=params)


def _target_taxon(params: MetadataParams) -> str:
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


def _normalize_api_tax(row: dict) -> dict:
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
            tax = _normalize_api_tax(row)
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
    for species in _species_or_none(params):
        taxon = _target_taxon(replace(params, target_species=species))
        logger.info("Querying GTDB API for genomes in %s", taxon)
        found = _api_genomes_detail(taxon, sp_reps, logger)
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

    # Quality for every candidate first: the cut below ranks on it.
    cards = _api_cards_by_accession([row["gid"] for row in rows], logger)
    quality = {acc: _api_quality(card) if card else (None, None) for acc, card in cards.items()}
    if counts is not None:
        counts["candidates"] = len(rows)
    if params.limit and len(rows) > params.limit:
        candidates = [
            _Candidate(
                accession=row["gid"],
                species=_normalize_api_tax(row)["species"],
                is_rep=bool(row.get("gtdbIsRep", False)),
                completeness=quality[row["gid"]][0],
                contamination=quality[row["gid"]][1],
            )
            for row in rows
        ]
        kept = _stratified_limit(candidates, params.limit)
        _log_limit(params.limit, candidates, kept, logger)
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
                _normalize_api_tax(row),
                completeness=completeness,
                contamination=contamination,
                gtdb_representative=bool(row.get("gtdbIsRep", False)),
            )
        )
        selected_acc.add(acc)

    if params.no_outgroup:
        return records, None, details
    outgroup, outgroup_details = _select_outgroup_via_api(params, rows, selected_acc, logger)
    details[outgroup.accession] = outgroup_details
    return records, outgroup, details


def _api_genomes_detail(taxon: str, sp_reps_only: bool, logger) -> list[dict]:
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
    params: MetadataParams, rows: list[dict], selected_acc: set[str], logger
) -> tuple[GenomeRecord, GenomeDetails]:
    # Explicit outgroup: fetch its card for taxonomy.
    if params.outgroup_accession:
        _refuse_outgroup_in_selection(params.outgroup_accession, selected_acc)
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
            _refuse_outgroup_in_target(
                params.outgroup_accession, params.level, target, tax_row.get(field, "")
            )
        completeness, contamination = _api_quality(card)
        record = _record_from_tax(
            params.outgroup_accession,
            _normalize_api_tax(tax_row),
            is_outgroup=True,
            completeness=completeness,
            contamination=contamination,
        )
        return record, _api_details(tax_row, card)

    # Otherwise pick a representative from the parent taxon (one rank up) that is
    # not part of the selection. Use the same rank order as the TSV path.
    parent_rank = {"species": "genus", "genus": "family", "family": "class"}[params.level]
    parent_field = {"genus": "gtdbGenus", "family": "gtdbFamily", "class": "gtdbClass"}[parent_rank]
    parent_taxon = rows[0].get(parent_field)
    if not parent_taxon:
        raise WorkdirError("Could not determine a parent taxon for outgroup selection.")

    logger.info("Selecting outgroup from parent taxon %s", parent_taxon)
    parent_rows = _api_genomes_detail(parent_taxon, sp_reps_only=True, logger=logger)
    for row in parent_rows:
        if row["gid"] in selected_acc:
            continue
        if not row.get("gtdbIsRep", False):
            continue
        # must be outside the selected sub-taxa (every species given)
        if row.get(_rank_field(params.level)) in {r.get(_rank_field(params.level)) for r in rows}:
            continue
        card = _api_cards_by_accession([row["gid"]], logger)[row["gid"]]
        completeness, contamination = _api_quality(card) if card else (None, None)
        record = _record_from_tax(
            row["gid"],
            _normalize_api_tax(row),
            is_outgroup=True,
            completeness=completeness,
            contamination=contamination,
            gtdb_representative=True,
        )
        return record, _api_details(row, card)
    raise WorkdirError("Could not determine an outgroup via the API; specify --outgroup-accession.")


def _rank_field(level: str) -> str:
    return {"family": "gtdbFamily", "genus": "gtdbGenus", "species": "gtdbSpecies"}[level]
