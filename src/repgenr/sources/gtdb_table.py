"""GTDB metadata table as a selection source (``--source tsv``).

Downloads the full GTDB metadata table of a release and parses it. Robust and
release-pinned, but the whole table is transferred. The outgroup can be chosen
from the GTDB reference tree (the sister group of the target).
"""

from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
import tarfile
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

from ..core import gtdb_tree, http
from ..core.context import WorkdirContext
from ..core.contracts import (
    atomic_path,
    sanitise_taxon_tokens,
)
from ..core.errors import RepGenRError, UserInputError, WorkdirError
from ..core.manifest import GenomeRecord
from ..core.outgroup_choice import (
    DEFAULT_OUTGROUP_MODE,
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


TAXONOMY = ("domain", "phylum", "class", "family", "genus", "species")


GTDB_VERSIONS = ("bac120", "ar53")


RELEASE_RE = re.compile(r"^\d+\.\d+$")


def _select_via_tsv(
    ctx: WorkdirContext, params: MetadataParams, logger, counts: dict[str, int] | None = None
) -> tuple[list[GenomeRecord], GenomeRecord | None, dict[str, GenomeDetails]]:
    metadata_file = _obtain_metadata(ctx, params, logger)
    accessions = parse_metadata(metadata_file, params, logger)

    targets = []
    for species in species_or_none(params):
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
    outgroup_acc, outgroup_data, choice = _pick_outgroup(
        accessions,
        selected,
        targets,
        params,
        tree_loader=lambda: _reference_tree(ctx, params, logger),
        logger=logger,
    )
    write_choice(ctx.workdir, choice)
    if choice is not None:
        log_choice(choice, logger)

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


def _reference_tree(ctx: WorkdirContext, params: MetadataParams, logger) -> gtdb_tree.ReferenceTree:
    path = gtdb_tree.obtain_reference_tree(
        ctx.workdir,
        params.release,
        params.version,
        tree_path=params.gtdb_tree_path,
        nodownload=params.nodownload,
        logger=logger,
    )
    logger.info("Reading the GTDB reference tree %s", path.name)
    return gtdb_tree.read_tree(path)


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


def parse_metadata(path: Path, params: MetadataParams, logger) -> dict[str, dict]:
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

            tax = parse_taxonomy(fields[idx["gtdb_taxonomy"]])
            completeness = _opt_column(fields, idx, "checkm2_completeness", "checkm_completeness")
            contamination = _opt_column(
                fields, idx, "checkm2_contamination", "checkm_contamination"
            )
            # Empty when the table lacks the column (not the GTDB accession).
            accession_ncbi = _text_column(fields, idx, "ncbi_genbank_assembly_accession")
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
                # Kept for the choice of an outgroup among representatives.
                "choice": _choice_fields(fields, idx) if is_rep else None,
            }
    logger.info("Parsed %d accessions", len(accessions))
    return accessions


def _choice_fields(fields: list[str], idx: dict[str, int]) -> dict:
    n50 = _opt_column(fields, idx, "n50_contigs")
    contigs = _opt_column(fields, idx, "contig_count")
    return {
        "n50": int(n50) if n50 is not None else None,
        "contigs": int(contigs) if contigs is not None else None,
        "category": sys.intern(_text_column(fields, idx, "ncbi_genome_category")),
        "type_material": sys.intern(_text_column(fields, idx, "gtdb_type_designation_ncbi_taxa")),
    }


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


def parse_taxonomy(raw: str) -> dict[str, str]:
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
        LimitCandidate(
            accession=acc,
            species=data["tax"]["species"],
            is_rep=bool(data.get("is_rep")),
            completeness=data.get("completeness"),
            contamination=data.get("contamination"),
        )
        for acc, data in matching.items()
    ]
    kept = stratified_limit(candidates, limit)
    log_limit(limit, candidates, kept, logger)
    return {c.accession: matching[c.accession] for c in kept}


def _pick_outgroup(accessions, selected, targets, params, *, tree_loader=None, logger=None):
    """The outgroup accession, its table row and the :class:`Choice` behind an
    automatic pick (None for a named outgroup). See
    :mod:`repgenr.core.outgroup_choice` for the modes."""
    if isinstance(targets, dict):
        targets = [targets]
    target_vals = {t[params.level] for t in targets}
    if params.outgroup_accession:
        if params.outgroup_accession not in accessions:
            raise UserInputError(
                f"Outgroup accession {params.outgroup_accession} not in GTDB metadata."
            )
        refuse_outgroup_in_selection(params.outgroup_accession, selected)
        data = accessions[params.outgroup_accession]
        for target in sorted(target_vals):
            refuse_outgroup_in_target(
                params.outgroup_accession, params.level, target, data["tax"][params.level]
            )
        return params.outgroup_accession, data, None

    def in_target(data) -> bool:
        return any(all(data["tax"].get(k) == v for k, v in t.items()) for t in targets)

    def eligible(acc, data) -> bool:
        # Under --limit, genomes of the target that the cap left out are
        # neither selected nor an outgroup.
        return data["is_rep"] and acc not in selected and not in_target(data)

    levels = list(targets[0])
    upper = levels[-2] if len(levels) >= 2 else levels[-1]
    upper_val = next(iter(selected.values()))["tax"][upper]
    mode = getattr(params, "outgroup_mode", DEFAULT_OUTGROUP_MODE)
    choice = None
    if mode == "sister":
        choice = _sister_choice(accessions, in_target, eligible, tree_loader, logger)
    if choice is None and mode in ("sister", "quality"):
        parent = [
            _choice_candidate(acc, data)
            for acc, data in accessions.items()
            if eligible(acc, data) and data["tax"][upper] == upper_val
        ]
        choice = choose_by_quality(parent)
    if choice is not None:
        return choice.chosen.accession, accessions[choice.chosen.accession], choice
    # --outgroup-mode first: the first representative one rank up, in table order.
    for acc, data in accessions.items():
        if eligible(acc, data) and data["tax"][upper] == upper_val:
            return acc, data, None
    raise WorkdirError(
        "Could not determine an outgroup; specify --outgroup-accession, or select none "
        "with --no-outgroup."
    )


def _choice_candidate(acc: str, data: dict, group: int | None = None, distance=None) -> Candidate:
    q = data.get("choice") or {}
    return Candidate(
        accession=acc,
        species=data["tax"].get("species_name") or data["tax"].get("species", ""),
        completeness=data.get("completeness"),
        contamination=data.get("contamination"),
        n50=q.get("n50"),
        contigs=q.get("contigs"),
        category=q.get("category", ""),
        type_material=q.get("type_material", ""),
        group=group,
        distance=distance,
    )


def _sister_choice(accessions, in_target, eligible, tree_loader, logger) -> Choice | None:
    """The best genome of the nearest group outside the target in the GTDB
    tree; None (with a warning) when the tree cannot place the target."""
    if tree_loader is None:
        return None
    try:
        tree = tree_loader()
    except (RepGenRError, OSError, ValueError) as exc:
        if logger:
            logger.warning(
                "The GTDB reference tree is not available (%s); choosing the outgroup "
                "by quality within the parent taxon. Give the tree with --gtdb-tree-path "
                "for the sister-group choice.",
                exc,
            )
        return None
    target = [acc for acc, data in accessions.items() if data["is_rep"] and in_target(data)]
    placement = gtdb_tree.place(tree, target)
    if placement is None:
        if logger:
            logger.warning(
                "None of the %d target representatives is in the GTDB reference tree; "
                "choosing the outgroup by quality within the parent taxon.",
                len(target),
            )
        return None
    if logger:
        logger.info(
            "Target placed in the GTDB reference tree (%d of %d representatives)",
            placement.found,
            len(target),
        )
        if placement.intruders:
            logger.warning(
                "The target is not a single clade in the GTDB reference tree: %d other "
                "representatives lie within it. The outgroup is taken from outside "
                "that whole clade.",
                placement.intruders,
            )

    def groups():
        for number, ring in enumerate(gtdb_tree.sister_rings(tree, placement), 1):
            yield [
                _choice_candidate(acc, accessions[acc], number, dist)
                for acc, dist in ring
                if acc in accessions and eligible(acc, accessions[acc])
            ]

    return choose_by_groups(groups())


class GtdbTableSource(GenomeSource):
    name = "tsv"
    manifest_source = MANIFEST_SOURCE
    record_tool = "gtdb-table"
    stats_source = "gtdb-tsv"

    def validate(self, params: MetadataParams) -> None:
        if not params.release or not RELEASE_RE.match(params.release):
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

    def select(self, ctx, params, logger, counts) -> Selection:
        return _select_via_tsv(ctx, params, logger, counts)
