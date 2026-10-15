"""NCBI Virus genome selection (the default viral back-end).

From NCBI Virus records (``virus_records.json``) plus the downloaded
``download.fa``, select sequences by taxonomy and length, write canonical
genomes and a ``selection.tsv``, optionally group a segmented virus's segments
per isolate, and pick an outgroup with mashtree. This is the records-based path
of the vgenome stage, kept separate from the legacy BV-BRC path.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, median
from typing import TYPE_CHECKING

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ..core.context import WorkdirContext
from ..core.contracts import (
    SEGMENTS_TSV,
    SelectionRow,
    genome_filename,
    write_segments,
    write_selection,
)
from ..core.errors import UserInputError
from ..core.manifest import record_from_selection
from ..core.process import staged_dir
from ..core.sketches import sketch_stage_genomes
from . import _outgroup
from ._common import (
    parse_custom_filter,
    parse_targets,
    select_outgroup_from_matrix,
)
from .ncbi_virus import UNLABELLED_SEGMENTS, VirusRecord, normalise_segment, read_records

if TYPE_CHECKING:
    from ..stages.vgenome import VgenomeParams


def run_records(
    ctx: WorkdirContext,
    params: VgenomeParams,
    download_wd: Path,
    fasta: Path,
    records_json: Path,
    logger: logging.Logger,
) -> int:
    targets = parse_targets(params)
    if not targets:
        raise UserInputError("Select a taxonomy: --target-genus/-species/-serotype/-custom.")
    records = read_records(records_json)
    targets = _resolve_species_targets(records, targets, logger)
    selected = [r for r in records if _record_matches(r, targets)]
    if not selected:
        raise UserInputError("No sequences matched the taxonomy selection.")

    discard = [x.strip() for x in params.discard.split(",")] if params.discard else None
    seqs = _seq_map(fasta)
    if params.group_segments:
        # Segmented viruses: skip the per-segment length filter (segments differ
        # in length); group whole isolates instead. The kept records' length
        # span still serves as the window for outgroup candidates.
        kept = [r for r in selected if r.accession in seqs]
        lo, hi = min(r.length for r in kept), max(r.length for r in kept)
    else:
        lo, hi = _length_range_records(selected, params, logger)
        kept = [r for r in selected if _passes_length(r, lo, hi, discard, seqs)]
    if not kept:
        raise UserInputError("No sequences passed length/discard filtering.")

    if params.print_fasta_headers:
        for r in kept:
            logger.info(seqs[r.accession].description)
    if params.glance:
        logger.info("--glance specified; stopping before writing genomes")
        return 0

    # The genomes are built beside genomes/ and swapped in only when every
    # file is written; selection.tsv follows. A crash mid-write leaves the
    # previous genome set and its selection table as they were.
    segments: dict[str, list[str]] = {}
    labels: dict[str, tuple[str, str]] = {}
    with staged_dir(ctx.genomes_dir) as genomes_dir:
        if params.group_segments:
            selection_rows = _write_isolate_groups(genomes_dir, kept, seqs, logger, segments)
            # segments.tsv keeps each member's label as submitted beside the
            # normalised segment the grouping used.
            by_acc = {r.accession: r for r in kept}
            labels = {
                acc: (normalise_segment(by_acc[acc].segment), by_acc[acc].segment)
                for accs in segments.values()
                for acc in accs
            }
        else:
            selection_rows = []
            for r in kept:
                name = genome_filename(r.family, r.genus, r.species, r.accession)
                rec = seqs[r.accession]
                (genomes_dir / name).write_text(f">{rec.description}\n{rec.seq}\n")
                selection_rows.append(
                    SelectionRow(r.accession, r.family, r.genus, r.species, False, name)
                )

    tool_versions: dict[str, str] = {}
    og = None
    if params.outgroup_accession:
        og = _pinned_outgroup_record(ctx, records, kept, params.outgroup_accession, seqs, logger)
    elif not params.no_outgroup:
        og, tool_versions = _determine_outgroup_records(
            ctx, records, kept, (lo, hi), params, seqs, logger
        )
    if og is None:
        _outgroup.clear_stale_outgroup(ctx, None)
    else:
        selection_rows.append(
            SelectionRow(
                og.accession,
                og.family,
                og.genus,
                og.species,
                True,
                genome_filename(og.family, og.genus, og.species, og.accession),
            )
        )

    write_selection(ctx.workdir / "selection.tsv", selection_rows)
    # The member accessions behind each grouped isolate, for tree2tax's
    # genomes_map; a run without grouping must not leave a stale table.
    segments_path = ctx.workdir / SEGMENTS_TSV
    if segments:
        write_segments(segments_path, segments, labels)
    else:
        segments_path.unlink(missing_ok=True)
    ctx.manifest.replace_genomes([record_from_selection(r, "ncbi_virus") for r in selection_rows])
    n_written = sum(1 for r in selection_rows if not r.is_outgroup)
    sketches, sketch_versions = sketch_stage_genomes(ctx, params.sketch, "vgenome", logger)
    ctx.config.record_stage(
        "vgenome",
        # The outgroup search is the one external tool; none when it did not run.
        tool=params.outgroup_treebuilder if tool_versions else None,
        params={
            "source": "ncbi_virus",
            "selected": n_written,
            "group_segments": params.group_segments,
            "no_outgroup": params.no_outgroup,
            "outgroup_accession": params.outgroup_accession,
            "sketches": sketches,
        },
        tool_versions={**tool_versions, **sketch_versions},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Wrote %d viral genomes (NCBI Virus)", n_written)
    return n_written


def _norm(value: str) -> str:
    return value.lower().replace("-", "").replace("_", "").replace(" ", "")


def _resolve_species_targets(
    records: list[VirusRecord], targets: dict[str, list[str]], logger: logging.Logger
) -> dict[str, list[str]]:
    """Add the species of records whose organism name a --target-species value is.

    The species of a record is the ICTV binomial from the lineage, so a value
    such as 'Junin virus' (an organism name) matches no species. A value
    stands for the species it names and for the species of every record whose
    organism name it is, so the whole species is selected rather than only the
    records that carry that name. A value that is both ('Hepatovirus A' is a
    species token without a binomial and an organism name under 'Hepatovirus
    ahepa') selects the union, and the log names the added species.
    """
    values = targets.get("species")
    if not values:
        return targets
    resolved: list[str] = []
    for value in values:
        resolved.append(value)
        named = sorted(
            {
                r.species
                for r in records
                if _norm(r.organism) == _norm(value) and _norm(r.species) != _norm(value)
            }
        )
        if named:
            logger.info(
                "--target-species '%s' is an organism name; also selecting its species %s",
                value,
                ", ".join(named),
            )
            resolved.extend(named)
    return {**targets, "species": resolved}


def _record_matches(rec, targets: dict[str, list[str]]) -> bool:
    """All requested taxonomy levels must match the record (normalized compare)."""
    for level, values in targets.items():
        norm_values = [_norm(v) for v in values]
        if level == "genus":
            ok = _norm(rec.genus) in norm_values
        elif level == "species":
            ok = _norm(rec.species) in norm_values
        elif level == "serotype":
            ok = any(_norm(v) in _norm(rec.organism) for v in values)
        elif level == "custom":

            def _custom_ok(kv: str) -> bool:
                key, val = parse_custom_filter(kv)
                return _norm(str(getattr(rec, key, ""))) == _norm(val)

            ok = all(_custom_ok(kv) for kv in values)
        else:
            ok = False
        if not ok:
            return False
    return bool(targets)


def _length_range_records(selected, params: VgenomeParams, logger) -> tuple[int, int]:
    if params.length_range:
        try:
            lo, hi = (int(x) for x in params.length_range.split("-"))
        except ValueError as exc:
            raise UserInputError("--length-range must be 'start-end', e.g. 25000-35000") from exc
        return lo, hi
    by_species: dict[str, list[int]] = {}
    for r in selected:
        by_species.setdefault(r.species, []).append(r.length)
    all_lens = [r.length for r in selected]
    if params.length_all:
        return min(all_lens), max(all_lens)
    # median_of_medians is a deliberate over-representation defense: one
    # median per species, then the median of those, so 900 outbreak records
    # of one species carry a single vote. The 'mean' alternative averages
    # every record and lets the most-sequenced species set the window.
    medians = [int(median(v)) for v in by_species.values()]
    midpoint = median(medians) if params.length_method == "median_of_medians" else mean(all_lens)
    dev = params.length_deviation / 100
    rng = (int(midpoint * (1 - dev)), int(midpoint * (1 + dev)))
    logger.info("Target length range: %d-%d", rng[0], rng[1])
    return rng


def _seq_map(fasta: Path) -> Mapping[str, SeqRecord]:
    """accession -> SeqRecord, from an NCBI Virus genomic.fna (headers '>acc ...').

    Uses an on-disk offset index (``SeqIO.index``) rather than loading every
    record into a dict: only the kept/outgroup accessions are materialized on
    access, so memory stays bounded for a large genus download.
    """
    return SeqIO.index(str(fasta), "fasta")


def _passes_length(rec: VirusRecord, lo: int, hi: int, discard, seqs) -> bool:
    if not (lo <= rec.length <= hi) or rec.accession not in seqs:
        return False
    if discard:
        desc = seqs[rec.accession].description
        if any(tag in desc for tag in discard):
            return False
    return True


def _isolate_token(isolate: str) -> str:
    """Make an isolate name a single safe accession-like token."""
    token = re.sub(r"[^A-Za-z0-9.-]", "", isolate.replace(" ", "-").replace("_", "-"))
    return f"iso-{token or 'NA'}"


def _segment_of(r: VirusRecord) -> str:
    """The normalised segment of a record; empty when it names no segment."""
    segment = normalise_segment(r.segment)
    return "" if segment in UNLABELLED_SEGMENTS else segment


def _segment_labels(recs) -> set[str]:
    """Distinct segments of a record set, by normalised label ('M', 'M; medium'
    and 'middle' are one segment); unlabelled records are ignored."""
    return {s for s in map(_segment_of, recs) if s}


def _segment_rank(r: VirusRecord) -> tuple[bool, int, str]:
    """Which of several records of one segment an isolate keeps: complete
    first, then the longest, then the lowest accession (deterministic)."""
    return (r.completeness != "COMPLETE", -r.length, r.accession)


_IsolateKey = tuple[str, str]  # (species, isolate)


def _isolate_segment_sets(
    records: list[VirusRecord], logger: logging.Logger
) -> tuple[dict[_IsolateKey, list[VirusRecord]], list[VirusRecord]]:
    """Split records into isolate segment sets and singletons.

    Records group by (species, isolate): isolate names repeat across species
    (live, Mammarenavirus: 'Acar 3080' carries a Lassa S and a Mobala L), so
    the name alone concatenated two species. A group is a segment set only
    with at least two distinct segments (by normalised label, see
    :func:`~repgenr.viral.ncbi_virus.normalise_segment`); within it, an unlabelled record
    stays a singleton and a segment submitted more than once (Lassa 'Josiah':
    three L and three S records) is kept once.
    """
    by_key: dict[_IsolateKey, list[VirusRecord]] = {}
    singletons: list[VirusRecord] = []
    for r in records:
        if r.isolate:
            by_key.setdefault((r.species, r.isolate), []).append(r)
        else:
            singletons.append(r)
    groups: dict[_IsolateKey, list[VirusRecord]] = {}
    duplicates = 0
    for key, recs in by_key.items():
        labelled = [r for r in recs if _segment_of(r)]
        if len(_segment_labels(recs)) < 2:
            singletons.extend(recs)
            continue
        singletons.extend(r for r in recs if r not in labelled)
        per_segment: dict[str, list[VirusRecord]] = {}
        for r in labelled:
            per_segment.setdefault(_segment_of(r), []).append(r)
        kept = [min(rs, key=_segment_rank) for rs in per_segment.values()]
        duplicates += len(labelled) - len(kept)
        groups[key] = kept
    if duplicates:
        logger.info(
            "Left out %d repeated segment record(s): each isolate keeps one record per segment",
            duplicates,
        )
    return groups, singletons


def _isolate_tokens(groups: dict[_IsolateKey, list[VirusRecord]]) -> dict[_IsolateKey, str]:
    """A unique accession-like token per isolate group.

    Two groups whose names sanitise alike ('Candid #1' and 'Candid-1') or that
    share an isolate name across species get the first member accession
    appended, so filenames and manifest accessions never collide.
    """
    base = {key: _isolate_token(key[1]) for key in groups}
    counts: dict[str, int] = {}
    for token in base.values():
        counts[token] = counts.get(token, 0) + 1
    out: dict[_IsolateKey, str] = {}
    for key, token in base.items():
        if counts[token] > 1:
            first = min(r.accession for r in groups[key])
            token = f"{token}-{re.sub(r'[^A-Za-z0-9.-]', '', first)}"
        out[key] = token
    return out


def _write_isolate_groups(genomes_dir, records, seqs, logger, segments=None):
    """Combine each isolate's segments into one genome; keep singletons as-is.

    An isolate's segments (one record per segment, see
    :func:`_isolate_segment_sets`) are concatenated in descending length order
    (a deterministic, segment-number-free ordering) into a single canonical
    genome; every other record is written individually.
    """
    groups, singletons = _isolate_segment_sets(records, logger)
    tokens = _isolate_tokens(groups)

    rows: list[SelectionRow] = []
    grouped = 0
    for key, recs in groups.items():
        iso = key[1]
        rep = recs[0]
        acc = tokens[key]
        name = genome_filename(rep.family, rep.genus, rep.species, acc)
        ordered = sorted(recs, key=lambda r: (-r.length, r.accession))
        seq = "".join(str(seqs[r.accession].seq) for r in ordered)
        (genomes_dir / name).write_text(f">{acc} {iso} ({len(ordered)} segments)\n{seq}\n")
        rows.append(SelectionRow(acc, rep.family, rep.genus, rep.species, False, name))
        if segments is not None:
            segments[acc] = [r.accession for r in ordered]
        grouped += len(recs)
    for r in singletons:
        name = genome_filename(r.family, r.genus, r.species, r.accession)
        (genomes_dir / name).write_text(
            f">{seqs[r.accession].description}\n{seqs[r.accession].seq}\n"
        )
        rows.append(SelectionRow(r.accession, r.family, r.genus, r.species, False, name))
    logger.info(
        "Grouped %d segment sequences into %d isolate genomes (+ %d single-record genomes)",
        grouped,
        len(groups),
        len(singletons),
    )
    return rows


def _determine_outgroup_records(ctx, records, kept, length_range, params, seqs, logger):
    """Pick an outgroup from species outside the selection, via mashtree.

    Returns ``(outgroup_record_or_None, tool_versions)`` so the caller can record
    the resolved mashtree version in provenance.
    """
    builder = _outgroup.resolve_outgroup_builder(params.outgroup_treebuilder)
    versions = _outgroup.preflight_outgroup_builder(builder)
    kept_species = {r.species for r in kept}
    kept_acc = {r.accession for r in kept}
    cand_by_species: dict[str, list] = {}
    for r in records:
        if r.species in kept_species or r.accession not in seqs:
            continue
        cand_by_species.setdefault(r.species, []).append(r)
    candidates = {
        sp: rs
        for sp, rs in cand_by_species.items()
        if len(rs) >= params.outgroup_candidates_taxid_min_genomes
    }
    # The tool versions are recorded only when the distance matrix ran.
    if not candidates:
        logger.warning("No outgroup candidates found; proceeding without an outgroup.")
        return None, {}

    outgroup_wd, gdir = _outgroup.prepare_workdir(ctx)
    lo, hi = length_range
    mid = (lo + hi) / 2
    rec_by_acc: dict[str, VirusRecord] = {}

    def _fits(length: int) -> bool:
        # Grouped isolates mix segments of very different lengths, so the
        # window is their whole span widened by the tolerance; a
        # single-record selection uses its midpoint plus/minus the tolerance.
        tolerance = _outgroup.RECORDS_LENGTH_TOLERANCE
        if getattr(params, "group_segments", False):
            return lo * (1 - tolerance) <= length <= hi * (1 + tolerance)
        return _outgroup.within_length_tolerance(length, mid, tolerance)

    def _emit(prefix: str, recs: list, cap: int) -> None:
        written = 0
        for r in recs:
            if written >= cap or not _fits(r.length):
                continue
            (gdir / f"{prefix}_{r.accession}.fasta").write_text(
                f">{r.accession}\n{seqs[r.accession].seq}\n"
            )
            rec_by_acc[r.accession] = r
            written += 1

    _emit("S", list(kept), 12)
    for rs in candidates.values():
        _emit("O", rs, 3)

    genome_files = sorted(gdir.glob("*.fasta"))
    if len([f for f in genome_files if f.name.startswith("O_")]) == 0:
        logger.warning("No length-compatible outgroup candidates; proceeding without one.")
        _outgroup.cleanup_workdir(outgroup_wd, params.keep_files)
        return None, {}
    matrix = _outgroup.run_distance_matrix(builder, genome_files, int(mid), outgroup_wd, logger)
    label = select_outgroup_from_matrix(matrix, logger) if matrix.exists() else None
    if label is None:
        logger.warning("Could not assign an outgroup; proceeding without one.")
        _outgroup.cleanup_workdir(outgroup_wd, params.keep_files)
        return None, versions
    acc = label[2:] if label.startswith(("O_", "S_")) else label
    og = rec_by_acc.get(acc)
    if og is None or acc in kept_acc:
        return None, versions

    _write_outgroup(ctx, og, seqs, logger)
    _outgroup.cleanup_workdir(outgroup_wd, params.keep_files)
    return og, versions


def _pinned_outgroup_record(ctx, records, kept, accession: str, seqs, logger):
    """The outgroup the user named: a downloaded record outside the selection."""
    og = next((r for r in records if r.accession == accession), None)
    if og is None or accession not in seqs:
        raise UserInputError(
            f"--outgroup-accession {accession} is not among the downloaded records; "
            "pick an accession the vmetadata stage fetched, or drop the flag."
        )
    if any(r.accession == accession for r in kept):
        raise UserInputError(
            f"--outgroup-accession {accession} is part of the selection itself; "
            "an outgroup must lie outside it."
        )
    _write_outgroup(ctx, og, seqs, logger)
    return og


def _write_outgroup(ctx, og, seqs, logger) -> None:
    """Stage the outgroup record under outgroup/ and record its accession."""
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    name = genome_filename(og.family, og.genus, og.species, og.accession)
    # Remove outgroups from earlier selections before writing the current one.
    _outgroup.clear_stale_outgroup(ctx, name)
    rec = seqs[og.accession]
    (ctx.outgroup_dir / name).write_text(f">{rec.description}\n{rec.seq}\n")
    (ctx.workdir / "outgroup_accession.txt").write_text(og.accession + "\n")
    logger.info("Selected outgroup: %s (%s)", og.accession, og.species)
