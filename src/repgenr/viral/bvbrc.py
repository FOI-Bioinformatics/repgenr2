"""Legacy BV-BRC genome selection (the --source bvbrc back-end).

Selects genomes from the BV-BRC group FASTA plus the Entrez-derived taxonomy
tables (``metadata_base.tsv`` / ``metadata_ncbi.tsv``) written by the vmetadata
stage. Kept for back-compatibility; the default path is NCBI Virus
(:mod:`repgenr.viral.selection`). The BV-BRC header form is
``... | <bvbrc_id>] [<species> | ...``, and the per-taxid length stats come from
the metadata tables rather than from per-record fields.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, median
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..core.context import WorkdirContext
from ..core.contracts import SELECTION_TSV, SelectionRow, write_selection
from ..core.errors import UserInputError, WorkdirError
from ..core.manifest import record_from_selection
from ..core.process import staged_dir
from ..core.sketches import sketch_stage_genomes
from . import _outgroup
from ._common import (
    parse_custom_filter,
    parse_targets,
    select_outgroup_from_matrix,
)
from .entrez import TAXNAMES_ORDERED

if TYPE_CHECKING:
    from ..stages.vgenome import VgenomeParams


def run_select(
    ctx: WorkdirContext,
    params: VgenomeParams,
    fasta: Path,
    base_tsv: Path,
    ncbi_tsv: Path,
    logger: logging.Logger,
) -> int:
    targets = parse_targets(params)
    if not targets:
        raise UserInputError("Select a taxonomy: --target-genus/-species/-serotype/-custom.")

    base = _read_base(base_tsv)
    ncbi = _read_ncbi(ncbi_tsv)

    # One metadata-only pass over the (large) group FASTA; sequences are fetched
    # lazily for just the records actually written. The old path re-parsed the
    # whole FASTA four times.
    records = _scan_fasta(fasta)

    selected, headers = _select_by_taxonomy(records, ncbi, targets, logger)
    if not selected:
        raise UserInputError("No sequences matched the taxonomy selection.")

    length_range = _determine_length_range(selected, base, params, logger)
    kept = _filter_by_length(selected, headers, length_range, params, logger)
    if not kept:
        raise UserInputError("No sequences passed length/discard filtering.")

    if params.print_fasta_headers:
        for taxid in kept:
            for bvbrc_id in kept[taxid]:
                logger.info(headers[taxid][bvbrc_id])
    if params.glance:
        logger.info("--glance specified; stopping before writing genomes")
        return 0

    sequences = _open_sequences(fasta, params.ignore_duplicates)
    n_written, rows = _write_genomes(ctx, records, sequences, kept, ncbi, params, logger)

    tool_versions: dict[str, str] = {}
    outgroup_id: str | None = None
    if params.outgroup_accession:
        outgroup_id = _pinned_outgroup(ctx, records, sequences, kept, params.outgroup_accession)
    elif not params.no_outgroup:
        tool_versions, outgroup_id = _determine_outgroup(
            ctx, records, sequences, base, kept, length_range, params, logger
        )
    _outgroup.clear_stale_outgroup(ctx, f"{outgroup_id}.fasta" if outgroup_id else None)
    if outgroup_id is not None:
        og_taxid = next((r.taxid for r in records if r.name == outgroup_id), "")
        rows.append(_selection_row(outgroup_id, og_taxid, ncbi, is_outgroup=True))
    # The same hand-off the NCBI Virus path and the bacterial stages publish,
    # so dereplicate, tree2tax and doctor see one contract.
    write_selection(ctx.workdir / SELECTION_TSV, rows)
    ctx.manifest.replace_genomes([record_from_selection(r, "bvbrc") for r in rows])
    sketches, sketch_versions = sketch_stage_genomes(ctx, params.sketch, "vgenome", logger)

    ctx.config.record_stage(
        "vgenome",
        # The outgroup search is the one external tool; none when it did not run.
        tool=params.outgroup_treebuilder if tool_versions else None,
        params={
            "selected": n_written,
            "no_outgroup": params.no_outgroup,
            "sketches": sketches,
        },
        tool_versions={**tool_versions, **sketch_versions},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Wrote %d viral genomes", n_written)
    return n_written


def _read_base(path: Path) -> dict[str, dict]:
    base: dict[str, dict] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        for f in reader:
            taxid = f[0]
            base[taxid] = {
                "num": int(f[2]),
                "seq_min": int(f[3]),
                "seq_max": int(f[4]),
                "seq_med": int(f[5]),
                "seq_mean": int(f[6]),
            }
    return base


def _read_ncbi(path: Path) -> dict[str, list[dict]]:
    ncbi: dict[str, list[dict]] = {}
    n = len(TAXNAMES_ORDERED)
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        for f in reader:
            if len(f) < 3 + 2 * n:
                continue
            taxid = f[0]
            names = f[3 : 3 + n]
            taxids = f[3 + n : 3 + 2 * n]
            ncbi[taxid] = [
                {"taxlevelname": TAXNAMES_ORDERED[i], "taxname": names[i], "taxid": taxids[i]}
                for i in range(n)
            ]
    return ncbi


def _matches(ncbi: dict, taxid: str, level: str, values: list[str]) -> bool:
    for entry in ncbi.get(taxid, []):
        if entry["taxlevelname"] == level and (
            entry["taxname"].lower() in values or entry["taxid"] in values
        ):
            return True
    return False


@dataclass
class _Record:
    """Per-record BV-BRC metadata (no sequence held)."""

    name: str  # record id / first header token; also the output filename stem
    bvbrc_id: str
    taxid: str
    description: str
    length: int


def _scan_fasta(fasta: Path) -> list[_Record]:
    """One pass over the group FASTA -> per-record metadata (sequences not held)."""
    records: list[_Record] = []
    for rec in SeqIO.parse(str(fasta), "fasta"):
        bvbrc_id = rec.description.split("| ")[1].split("]")[0]
        records.append(
            _Record(
                name=rec.id,
                bvbrc_id=bvbrc_id,
                taxid=bvbrc_id.split(".")[0],
                description=rec.description,
                length=len(rec.seq),
            )
        )
    return records


def _open_sequences(fasta: Path, ignore_duplicates: bool) -> Mapping[str, SeqRecord]:
    """Lazy name -> SeqRecord index so only the written records load sequences.

    ``SeqIO.index`` rejects duplicate record ids; with --ignore-duplicates fall
    back to an in-memory dict (last record wins), matching the old overwrite
    behavior, otherwise raise a clear error.
    """
    try:
        return SeqIO.index(str(fasta), "fasta")
    except ValueError as exc:
        if not ignore_duplicates:
            raise WorkdirError(
                "Duplicate sequence id in the BV-BRC FASTA. Use --ignore-duplicates."
            ) from exc
        return {rec.id: rec for rec in SeqIO.parse(str(fasta), "fasta")}


def _select_by_taxonomy(records: list[_Record], ncbi: dict, targets: dict, logger):
    selected: dict[str, dict[str, int]] = {}
    headers: dict[str, dict[str, str]] = {}
    custom_pairs = [parse_custom_filter(kv) for kv in targets.get("custom", ())]
    for rec in records:
        ok = True
        for level, values in targets.items():
            if level == "custom":
                hit = any(_matches(ncbi, rec.taxid, key, [val]) for key, val in custom_pairs)
                if not hit:
                    ok = False
            elif not _matches(ncbi, rec.taxid, level, values):
                ok = False
        if not ok:
            continue

        selected.setdefault(rec.taxid, {})[rec.bvbrc_id] = rec.length
        headers.setdefault(rec.taxid, {})[rec.bvbrc_id] = rec.description
    logger.info("Taxonomy selection matched %d taxids", len(selected))
    return selected, headers


def _determine_length_range(selected, base, params: VgenomeParams, logger) -> tuple[int, int]:
    if params.length_range:
        try:
            lo, hi = (int(x) for x in params.length_range.split("-"))
        except ValueError as exc:
            raise UserInputError("--length-range must be 'start-end', e.g. 25000-35000") from exc
        return lo, hi

    stat_values: list[int] = []
    extremes: list[int] = []
    for taxid in selected:
        if taxid not in base:
            continue
        meta = base[taxid]
        stat_values.append(
            meta["seq_med"] if params.length_method == "median_of_medians" else meta["seq_mean"]
        )
        if params.length_all:
            extremes += [meta["seq_min"], meta["seq_max"]]
    if not stat_values:
        raise UserInputError("Could not determine target length; try --length-range.")

    if params.length_all:
        return min(extremes), max(extremes)
    midpoint = (
        median(stat_values) if params.length_method == "median_of_medians" else mean(stat_values)
    )
    dev = params.length_deviation / 100
    rng = (int(midpoint * (1 - dev)), int(midpoint * (1 + dev)))
    logger.info("Target length range: %d-%d", rng[0], rng[1])
    return rng


def _filter_by_length(selected, headers, length_range, params: VgenomeParams, logger):
    lo, hi = length_range
    discard = [x.strip() for x in params.discard.split(",")] if params.discard else None
    kept: dict[str, dict[str, int]] = {}
    for taxid, by_id in selected.items():
        for bvbrc_id, seq_len in by_id.items():
            if not (lo <= seq_len <= hi):
                continue
            if discard and any(tag in headers[taxid][bvbrc_id] for tag in discard):
                continue
            kept.setdefault(taxid, {})[bvbrc_id] = seq_len
    logger.info("After length/discard filtering: %d taxids", len(kept))
    return kept


def _write_record(seq: SeqRecord, dest: Path) -> None:
    dest.write_text(f">{seq.description}\n{seq.seq}\n", encoding="utf-8")


def _tax_name(ncbi: dict, taxid: str, level: str) -> str:
    for entry in ncbi.get(taxid, []):
        if entry["taxlevelname"] == level:
            return entry["taxname"]
    return ""


def _selection_row(name: str, taxid: str, ncbi: dict, *, is_outgroup: bool) -> SelectionRow:
    return SelectionRow(
        name,
        _tax_name(ncbi, taxid, "family"),
        _tax_name(ncbi, taxid, "genus"),
        _tax_name(ncbi, taxid, "species"),
        is_outgroup,
        f"{name}.fasta",
    )


def _write_genomes(
    ctx, records, sequences, kept, ncbi, params: VgenomeParams, logger
) -> tuple[int, list[SelectionRow]]:
    rows: dict[str, SelectionRow] = {}
    # Built beside genomes/ and swapped in only when complete (see staged_dir).
    with staged_dir(ctx.genomes_dir) as genomes_dir:
        for rec in records:
            if rec.taxid not in kept or rec.bvbrc_id not in kept[rec.taxid]:
                continue
            target = genomes_dir / f"{rec.name}.fasta"
            if target.exists():
                if not params.ignore_duplicates:
                    raise WorkdirError(
                        f"Duplicate sequence id {rec.name}. Use --ignore-duplicates to proceed."
                    )
                logger.warning("Duplicate sequence id %s; overwriting", rec.name)
            _write_record(sequences[rec.name], target)
            rows[rec.name] = _selection_row(rec.name, rec.taxid, ncbi, is_outgroup=False)
    return len(rows), list(rows.values())


def _pinned_outgroup(ctx, records, sequences, kept, record_id: str) -> str:
    """The outgroup the user named, by record id or BV-BRC id, outside the selection."""
    rec = next((r for r in records if record_id in (r.name, r.bvbrc_id)), None)
    if rec is None or rec.name not in sequences:
        raise UserInputError(
            f"--outgroup-accession {record_id} is not among the downloaded records; "
            "pick a record id the vmetadata stage fetched, or drop the flag."
        )
    if rec.taxid in kept:
        raise UserInputError(
            f"--outgroup-accession {record_id} belongs to a selected taxon; "
            "an outgroup must lie outside the selection."
        )
    (ctx.workdir / "outgroup_accession.txt").write_text(rec.name + "\n")
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    _write_record(sequences[rec.name], ctx.outgroup_dir / f"{rec.name}.fasta")
    return rec.name


def _determine_outgroup(
    ctx, records, sequences, base, kept, length_range, params: VgenomeParams, logger
) -> tuple[dict[str, str], str]:
    """Pick an outgroup via mashtree; return the resolved tool versions and its id."""
    builder = _outgroup.resolve_outgroup_builder(params.outgroup_treebuilder)
    versions = _outgroup.preflight_outgroup_builder(builder)
    outgroup_wd, genomes_dir = _outgroup.prepare_workdir(ctx)

    candidates = {
        taxid
        for taxid, meta in base.items()
        if taxid not in kept and meta["num"] >= params.outgroup_candidates_taxid_min_genomes
    }
    if not candidates:
        raise WorkdirError(
            "No outgroup candidates found. Lower the minimum genomes per candidate "
            "taxid, or pass --no-outgroup."
        )

    max_per_taxid = 3
    written: dict[str, int] = {}
    for rec in records:
        if written.get(rec.taxid, 0) > max_per_taxid:
            continue
        if not (rec.taxid in kept or rec.taxid in candidates) or rec.taxid not in base:
            continue
        metric = base[rec.taxid]["seq_med"]
        if not _outgroup.within_length_tolerance(
            rec.length, metric, _outgroup.BVBRC_LENGTH_TOLERANCE
        ):
            continue
        prefix = "S" if rec.taxid in kept else "O"
        _write_record(sequences[rec.name], genomes_dir / f"{prefix}_{rec.name}.fasta")
        written[rec.taxid] = written.get(rec.taxid, 0) + 1

    genome_files = sorted(genomes_dir.glob("*.fasta"))
    matrix = _outgroup.run_distance_matrix(
        builder, genome_files, int(mean(length_range)), outgroup_wd, logger
    )
    if not matrix.exists():
        raise WorkdirError(
            f"{builder.capabilities.name} produced no distance matrix: either the tool "
            "failed (check the installation) or no candidate was within its reach "
            "(see the warning above); specify the outgroup manually."
        )

    outgroup_id = select_outgroup_from_matrix(matrix, logger)
    if outgroup_id is None:
        raise WorkdirError("Could not assign an outgroup; specify one manually.")
    if outgroup_id.startswith("O_"):
        outgroup_id = outgroup_id[2:]

    (ctx.workdir / "outgroup_accession.txt").write_text(outgroup_id + "\n")
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    if outgroup_id in sequences:
        _write_record(sequences[outgroup_id], ctx.outgroup_dir / f"{outgroup_id}.fasta")
    logger.info("Selected outgroup: %s", outgroup_id)
    _outgroup.cleanup_workdir(outgroup_wd, params.keep_files)
    return versions, outgroup_id
