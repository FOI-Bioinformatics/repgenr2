"""assemble stage: fetch the selected sequencing runs, assemble them, and write
the genome contract (``genomes/``, ``selection.tsv``, the manifest).

Each run in ``reads.tsv`` is fetched from the FASTQ locations ENA reported
(a local path is copied instead), verified against its checksum, assembled
with the adapter that accepts its platform and layout, filtered to contigs of
a minimum length, and labelled with the family, genus and species the reads
stage resolved. Runs that cannot be fetched or assembled are written to
``excused_runs.tsv`` so the completeness guard of later stages does not count
them as missing. Finished assemblies carry a marker under ``assemblies/<run>/``
and are skipped on a re-run. ``--max-runs N`` keeps at most N runs with a
finished assembly in the workdir; the remaining runs are recorded as
``deferred`` and are assembled by a later call with a larger N or without it.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import threading
import time
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..assemblers.base import AssembleParams as AdapterParams
from ..assemblers.base import ReadSet, accepting_assemblers, registry, select_assembler
from ..assemblers.contigs import ContigStats, filter_contigs
from ..classifiers.base import Classification, ClassifyParams
from ..classifiers.base import registry as classifier_registry
from ..classifiers.sourmash import gather_workers
from ..core import http, process
from ..core.context import WorkdirContext
from ..core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    FASTA_SUFFIXES,
    READS_TSV,
    SELECTION_TSV,
    AssemblyStatsRow,
    ExcusedRun,
    ReadRow,
    SelectionRow,
    accession_from_filename,
    atomic_replace,
    genome_filename,
    read_reads,
    read_selection,
    record_name,
    write_assembly_stats,
    write_excused_runs,
    write_selection,
)
from ..core.errors import (
    MissingBinaryError,
    RepGenRError,
    ToolExecutionError,
    UserInputError,
    WorkdirError,
)
from ..core.executors import parallel_map
from ..core.manifest import MANIFEST_FILENAME, SketchRecord, record_from_selection
from ..core.process import check_free_disk, link_or_copy, remove_tree, staged_dir
from ..core.sketches import (
    READS_SKETCH_NAME,
    READS_SKETCH_PARAMS,
    SketchSource,
    adapter_sketches,
    clear_partials,
    remove_stale,
    require_sourmash_if_requested,
    sketch_beside,
    sketch_reads,
    sketch_stage_genomes,
    sourmash_for_reads_sketch,
)
from ..polishers.base import PolishParams, accepting_polishers, select_polisher
from ..polishers.base import registry as polisher_registry
from .assemble_qc import (
    CHECKM2_CACHE,
    checkm2_db_from_env,
    load_cached_quality,
    preflight_checkm2,
    quality_cache_key,
    run_checkm2,
    store_cached_quality,
)
from .assemble_screen import (
    GTDB_LINEAGES_ENV,
    GTDB_SKETCH_ENV,
    SCREEN_DIR,
    SCREEN_FAILED,
    SCREEN_JSON,
    SCREEN_READS_TSV,
    SCREEN_STEP,
    ScreenGate,
    ScreenRecord,
    ScreenSettings,
    read_screen_json,
    resolve_settings,
    write_screen_json,
    write_screen_table,
)
from .ingest import OUTGROUP_ACCESSION_TXT
from .taxon_match import AGREE, GENUS_RENAMED, genus_agreement
from .taxon_match import gtdb_tokens as _gtdb_tokens

_DONE_MARKER = "assembly.ok"
_CONTIGS_NAME = "contigs.fasta"
# The sketch of an assembly, beside its contigs, with the parameters of the
# workdir sketches and named by run accession. The classifier gathers with it;
# the sketch step copies it into sketches/ once the genome is named.
_ASSEMBLY_SKETCH = "contigs.sig.zip"
# Rough peak disk per run: the FASTQ files, their uncompressed form and the
# assembler's scratch.
_DISK_FACTOR = 4


@dataclass
class AssembleParams:
    assembler: str = "auto"
    threads: int = 16
    # Concurrent runs; memory, not CPU, bounds this, so it stays small. None
    # resolves to 2, or to 1 when a long-read run is pending.
    jobs: int | None = None
    memory_gb: int = 16
    min_contig_length: int = 500
    # Correct long-read assemblies with the run's reads: auto picks by platform
    # (medaka for ONT, racon for PacBio CLR, none for HiFi and Illumina).
    polisher: str = "auto"
    polish_rounds: int = 1
    # A FASTA file to set aside as the outgroup (ingest semantics).
    outgroup: str | None = None
    # Add the assemblies to a working directory that already holds a selection
    # (metadata + genome, or ingest) instead of replacing it: existing rows,
    # files and outgroup stay, rows of the same accession are replaced.
    append: bool = False
    keep_reads: bool = False
    keep_files: bool = False
    # Quality: CheckM2 runs when a database is configured (flag or CHECKM2DB);
    # assemblies outside the gate are excused rather than selected.
    checkm2_db: str | None = None
    min_completeness: float = 50.0
    max_contamination: float = 10.0
    # Classification: verifies the submitted organism against a GTDB sketch
    # (flag or REPGENR_GTDB_SKETCH / REPGENR_GTDB_LINEAGES). "auto" runs the
    # sourmash classifier when a sketch is configured, "none" never.
    classifier: str = "auto"
    gtdb_sketch: str | None = None
    gtdb_lineages: str | None = None
    # Genome sketches (core.sketches): None sketches when sourmash can run,
    # True requires it, False skips them. --append sketches the new genomes only.
    sketch: bool | None = None
    # The reads sketch of each run (assemblies/<run>/reads.sig.zip), made by
    # sourmash beside the assembler: None sketches when sourmash can run, True
    # requires it, False skips it. Not a genome sketch; never in sketches/.
    reads_sketch: bool | None = None
    # Tool tuning from ``--tool-arg``; each adapter declares the keys it reads.
    extra: dict = field(default_factory=dict)
    # At most this many runs with a finished assembly (reused or new), taken in
    # reads.tsv order; the remaining pending runs are excused as deferred.
    # None assembles every run.
    max_runs: int | None = None
    # The reads screen (assemble_screen): after the fetch and the reads
    # sketch, before an assembler runs; needs the GTDB sketch and lineages.
    screen_reads: bool = False
    # Minimum abundance-weighted fraction of the reads in the top species.
    screen_min_fraction: float = 0.3
    # ANI estimate from which a run contained in an earlier run is a duplicate.
    screen_dup_ani: float = 0.999


@dataclass
class _Outcome:
    row: ReadRow
    assembler: str | None = None
    polisher: str | None = None
    polish_rounds: int = 0
    stats: ContigStats | None = None
    tool_stats: dict = field(default_factory=dict)
    excused: ExcusedRun | None = None
    quality: tuple[float, float] | None = None
    gtdb: Classification | None = None
    # Resolved filename tokens and where they came from.
    label: tuple[str, str, str] | None = None
    label_source: str = "metadata"
    taxonomy_flag: str = ""
    # Tool versions a reused marker recorded (its run needs no preflight).
    versions: dict = field(default_factory=dict)
    # Whether assemblies/<run>/reads.sig.zip holds the run's reads sketch.
    reads_sketch: bool = False
    # The reads screen's decision, when the run was screened.
    screen: ScreenRecord | None = None


def run(ctx: WorkdirContext, params: AssembleParams) -> int:
    logger = ctx.logger
    reads_path = ctx.workdir / READS_TSV
    if not reads_path.exists():
        raise WorkdirError(f"{READS_TSV} not found in {ctx.workdir}. Run the reads stage first.")
    rows = read_reads(reads_path)
    if not rows:
        raise WorkdirError(f"{READS_TSV} lists no runs.")
    _check_outgroup_file(params.outgroup)
    require_sourmash_if_requested(params.sketch)
    require_sourmash_if_requested(params.reads_sketch)
    screen = screen_settings(params)

    assemblies = ctx.workdir / "assemblies"
    scratch = ctx.scratch_dir / "assemble"
    assemblies.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)

    plan = _plan(rows, params, assemblies, logger=logger)
    if params.assembler == "auto":
        _excuse_missing_assemblers(plan, logger)
        _require_something_to_assemble(plan)
    reads_versions = sourmash_for_reads_sketch(
        True if screen is not None else params.reads_sketch, "assemble", logger
    )
    reads_sketch_version = reads_versions.get("sourmash", "") if reads_versions else None
    if reads_sketch_version is not None:
        # Before the screen: a finished run sketched here is screened in this call.
        refill_reads_sketches(
            plan, assemblies, scratch, reads_sketch_version, params.threads, logger
        )
    gate = None
    if screen is not None:
        gate = ScreenGate(
            screen,
            logger,
            gather_slots=gather_workers(len(plan), params.threads, params.memory_gb),
        )
        # Before --max-runs: a run excused from its kept sketch takes no place.
        screen_planned(plan, gate, assemblies, params.threads, logger)
    if params.max_runs is not None:
        defer_beyond(plan, params.max_runs, logger)
    if params.polisher == "auto":
        _warn_missing_polishers(plan, logger)
    checkm2_db = params.checkm2_db or checkm2_db_from_env()
    gtdb_sketch = params.gtdb_sketch or os.environ.get(GTDB_SKETCH_ENV)
    gtdb_lineages = params.gtdb_lineages or os.environ.get(GTDB_LINEAGES_ENV)
    classifier_name = classifier_for(params.classifier, gtdb_sketch)
    versions = {k: v for o in plan for k, v in o.versions.items()}
    versions.update(_preflight(plan, logger))
    # Databases and QC tools are checked before the assemblies, which can take hours.
    versions.update(
        check_quality_inputs(
            checkm2_db=checkm2_db,
            classifier=classifier_name,
            gtdb_sketch=gtdb_sketch,
            gtdb_lineages=gtdb_lineages,
        )
    )
    versions.update(reads_versions)
    pending = [o for o in plan if o.excused is None and o.stats is None]
    if gate is not None:
        gate.set_pending([o.row for o in pending])
    if pending:  # a rerun over finished runs downloads nothing
        check_free_disk(
            ctx.workdir,
            sum(sum(o.row.fastq_bytes) for o in pending) * _DISK_FACTOR,
            logger,
            what=f"assemble {len(pending)} sequencing runs",
        )
    requested = params.jobs if params.jobs is not None else _default_jobs(pending)
    jobs = max(1, min(requested, len(pending) or 1))
    threads_each = max(1, params.threads // jobs)
    n_deferred = sum(1 for o in plan if _is_deferred(o.excused))
    logger.info(
        "Assembling %d runs (%d already done, %d excused, %d deferred) with %d concurrent "
        "jobs, %d threads each",
        len(pending),
        sum(1 for o in plan if o.stats is not None),
        sum(1 for o in plan if o.excused is not None) - n_deferred,
        n_deferred,
        jobs,
        threads_each,
    )

    def work(outcome: _Outcome) -> _Outcome:
        try:
            return _fetch_and_assemble(
                outcome,
                params,
                threads_each,
                assemblies / outcome.row.run_accession,
                scratch / outcome.row.run_accession,
                versions,
                logger,
                reads_sketch=reads_sketch_version,
                gate=gate,
            )
        finally:
            # A run that ends before its screen must not hold up later ones.
            if gate is not None:
                gate.release(outcome.row.run_accession)

    done = parallel_map(work, pending, jobs, logger=logger)
    by_run = {o.row.run_accession: o for o in [*plan, *done]}
    outcomes = [by_run[r.run_accession] for r in rows]

    assembled = [o for o in outcomes if o.stats is not None and o.stats.n_contigs > 0]
    for o in assembled:
        o.label = _taxa(o.row)

    # Every run's contigs file is named alike; QC and classification see them
    # through links named by run accession, which are unique.
    named = named_links(
        {o.row.run_accession: assemblies / o.row.run_accession / _CONTIGS_NAME for o in assembled},
        scratch / "named",
    )
    query_sketches, sketch_records = _assembly_sketches(
        named, classifier_name, params, assemblies, logger
    )
    quality, classified = assess(
        named,
        checkm2_db=checkm2_db,
        classifier=classifier_name,
        gtdb_sketch=gtdb_sketch,
        gtdb_lineages=gtdb_lineages,
        threads=params.threads,
        extra=params.extra,
        scratch=scratch,
        versions=versions,
        logger=logger,
        memory_gb=params.memory_gb,
        cache_dirs={run: assemblies / run for run in named},
        sketches=query_sketches,
    )
    if quality is not None:
        apply_quality(assembled, quality, params.min_completeness, params.max_contamination, logger)
        assembled = [o for o in assembled if o.excused is None]
    # An excused run keeps no reads sketch; a deferred run was not judged, and
    # a run the screen excused keeps it so the reason can be inspected.
    for o in outcomes:
        if keeps_reads_sketch(o.excused):
            continue
        if o.excused is not None:
            drop_reads_sketch(o, assemblies / o.row.run_accession)
    if gate is not None:
        write_screen_table(
            ctx.workdir / SCREEN_READS_TSV,
            [o.screen for o in outcomes if o.screen is not None],
            [r.run_accession for r in rows],
        )
    n_disagree = 0
    if classified is not None:
        n_disagree = apply_classification(assembled, classified, versions, logger)

    excused = [o.excused for o in outcomes if o.excused is not None]
    unfetched = [e.run_accession for e in excused if e.reason.startswith("download_failed")]
    if unfetched:
        logger.warning(
            "%d run(s) could not be fetched (%s); a repeat with the same settings skips this "
            "stage, so rerun assemble with --force to try them again (finished runs are kept).",
            len(unfetched),
            ", ".join(unfetched[:5]) + (" ..." if len(unfetched) > 5 else ""),
        )
    excused_path = ctx.workdir / EXCUSED_RUNS_TSV
    if excused:
        write_excused_runs(excused_path, excused)
    else:
        excused_path.unlink(missing_ok=True)
    if not assembled:
        detail = _total_failure(ctx, params, excused, logger)
        raise WorkdirError(
            f"None of the {len(rows)} runs produced an accepted assembly ({detail}); see "
            f"{EXCUSED_RUNS_TSV}."
        )

    new_rows = [_selection_row(o) for o in assembled]
    if params.append:
        selection_rows = _append_rows(ctx, new_rows, params, logger)
        ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
        for o in assembled:
            _unlink_previous(ctx.genomes_dir, o.row.run_accession)
            link_or_copy(
                assemblies / o.row.run_accession / _CONTIGS_NAME, ctx.genomes_dir / _name(o)
            )
        ctx.manifest.upsert_many([record_from_selection(r, "sra") for r in new_rows])
    else:
        selection_rows = list(new_rows)
        outgroup_row = _stage_outgroup(ctx, params.outgroup, logger)
        if outgroup_row is not None:
            selection_rows.append(outgroup_row)
        with staged_dir(ctx.genomes_dir) as genomes_dir:
            for o in assembled:
                link_or_copy(
                    assemblies / o.row.run_accession / _CONTIGS_NAME, genomes_dir / _name(o)
                )
        ctx.manifest.replace_genomes([record_from_selection(r, "sra") for r in selection_rows])
    write_selection(ctx.workdir / SELECTION_TSV, selection_rows)
    write_assembly_stats(ctx.workdir / ASSEMBLY_STATS_TSV, [_stats_row(o) for o in assembled])
    outgroup_row = next((r for r in selection_rows if r.is_outgroup), None)
    sketches, sketch_versions = sketch_stage_genomes(
        ctx,
        params.sketch,
        "assemble",
        logger,
        threads=params.threads,
        only=[ctx.genomes_dir / _name(o) for o in assembled] if params.append else None,
        # The classifier's sketch of each accepted assembly, renamed to the genome.
        reuse={
            record_name(_name(o)): SketchSource(
                assemblies / o.row.run_accession / _ASSEMBLY_SKETCH,
                sketch_records[o.row.run_accession],
                name=o.row.run_accession,
            )
            for o in assembled
            if o.row.run_accession in sketch_records
        },
    )
    versions.update(sketch_versions)

    assemblers_used = sorted({o.assembler for o in assembled if o.assembler})
    ctx.config.record_stage(
        "assemble",
        # The assemblers that built the accepted genomes; params keep the request.
        tool=",".join(assemblers_used) or params.assembler,
        params={
            **asdict(params),
            # The databases as resolved (flag or environment variable).
            "checkm2_db": checkm2_db,
            "gtdb_sketch": gtdb_sketch,
            "gtdb_lineages": gtdb_lineages,
            "classifier_effective": classifier_name,
            "assemblers_used": assemblers_used,
            "polishers_used": sorted({o.polisher for o in assembled if o.polisher}),
            "n_assembled": len(assembled),
            "n_excused": len(excused) - n_deferred,
            "n_deferred": n_deferred,
            "n_disagree": n_disagree,
            "n_genus_renamed": sum(1 for o in assembled if o.taxonomy_flag == GENUS_RENAMED),
            "outgroup_accession": outgroup_row.accession if outgroup_row else None,
            "sketches": sketches,
            "n_reads_sketches": sum(1 for o in assembled if o.reads_sketch),
            "n_screened_out": sum(1 for e in excused if e.step == SCREEN_STEP),
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Assembled %d of %d runs (%d excused, %d deferred, %d classifier disagreements); wrote %s",
        len(assembled),
        len(rows),
        len(excused) - n_deferred,
        n_deferred,
        n_disagree,
        SELECTION_TSV,
    )
    if n_deferred:
        logger.info(
            "%d run(s) deferred by --max-runs %d are listed in %s; a later call with a "
            "larger --max-runs, or without it, assembles them. Finished runs keep their "
            "places, including ones the quality gate excused, so --force with the same "
            "--max-runs does not.",
            n_deferred,
            params.max_runs,
            EXCUSED_RUNS_TSV,
        )
    return len(assembled)


_LONG_READ_PLATFORMS = frozenset({"OXFORD_NANOPORE", "PACBIO_SMRT"})


def _why(exc: RepGenRError) -> str:
    """The failure for excused_runs.tsv: a tool failure keeps its output tail."""
    if isinstance(exc, ToolExecutionError) and exc.output:
        return f"{exc}: {exc.output}"
    return str(exc)


def planned_layout(row: ReadRow) -> str:
    """The layout a run is assembled as.

    ENA labels some short-read runs PAIRED while serving a single FASTQ file
    (mate-1 reads only, or interleaved pairs the assemblers here do not split).
    Such a run is planned as single-end, so an adapter that needs a read pair
    (shovill) excuses it before anything is downloaded. Long-read runs keep
    their label; their assemblers and polishers take any layout.
    """
    if (
        row.layout == "PAIRED"
        and len(row.fastq_urls) == 1
        and row.platform not in _LONG_READ_PLATFORMS
    ):
        return "SINGLE"
    return row.layout


def read_set(row: ReadRow, files: tuple[Path, ...] = ()) -> ReadSet:
    """The run as an adapter sees it, with the planned layout."""
    return ReadSet(
        row.run_accession, row.platform, row.instrument_model, planned_layout(row), files, row.bases
    )


def _default_jobs(pending: list[_Outcome]) -> int:
    """Two short-read assemblies fit side by side; a long-read one wants the machine."""
    return 1 if any(o.row.platform in _LONG_READ_PLATFORMS for o in pending) else 2


# --- planning -------------------------------------------------------------------


def _plan(
    rows: list[ReadRow],
    params: AssembleParams,
    assemblies: Path,
    *,
    check_settings: bool = True,
    logger: logging.Logger | None = None,
) -> list[_Outcome]:
    """Decide, per run, whether it is done, excused up front, or to be assembled.

    A finished run is reused when its marker's settings agree with ``params``
    (``check_settings``; ``reads-gather`` only collects results and passes
    False). A higher contig floor is applied to the finished contigs; any
    other difference assembles the run again.
    """
    if params.assembler != "auto" and params.assembler not in registry.names():
        raise UserInputError(
            f"Unknown assembler {params.assembler!r}; available: {', '.join(registry.names())}."
        )
    if params.polisher not in ("auto", "none") and params.polisher not in polisher_registry.names():
        raise UserInputError(
            f"Unknown polisher {params.polisher!r}; available: "
            f"{', '.join(polisher_registry.names())}, none."
        )
    plan = []
    for row in rows:
        outcome = _Outcome(row=row)
        if _reuse_finished(outcome, params, assemblies / row.run_accession, check_settings, logger):
            pass
        elif not row.fastq_urls:
            outcome.excused = ExcusedRun(row.run_accession, "fetch", "no_fastq_mirror")
        else:
            reads = read_set(row)
            if reads.layout != row.layout and logger is not None:
                logger.info(
                    "%s: listed as %s with one FASTQ file; planned as single-end",
                    row.run_accession,
                    row.layout,
                )
            reason = "unsupported_platform"
            if params.assembler == "auto":
                outcome.assembler = select_assembler(registry, reads)
            else:
                adapter = registry.create(params.assembler)
                if adapter.accepts(reads):
                    outcome.assembler = params.assembler
                elif reads.platform in adapter.read_types:
                    reason = "unsupported_layout"
            if outcome.assembler is None:
                outcome.excused = ExcusedRun(row.run_accession, "assemble", reason)
            elif params.polisher == "auto":
                outcome.polisher = select_polisher(polisher_registry, reads)
            elif params.polisher != "none":
                probe_cls = polisher_registry.get(params.polisher)
                if probe_cls.__new__(probe_cls).accepts(reads):
                    outcome.polisher = params.polisher
        plan.append(outcome)
    return plan


# --- reuse of finished runs ---------------------------------------------------------


def _read_marker(path: Path) -> dict | None:
    """A finished run's marker, or None when absent or unreadable (cut short by a kill)."""
    try:
        done = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(done, dict) or "assembler" not in done or "stats" not in done:
        return None
    return done


def _settings(params: AssembleParams) -> dict:
    """The result-affecting settings a marker records, as requested."""
    return {
        "assembler": params.assembler,
        "polisher": params.polisher,
        "polish_rounds": params.polish_rounds,
        "min_contig_length": params.min_contig_length,
        "extra": dict(sorted(params.extra.items())),
    }


def _accepted_extras(name: str | None, reg) -> frozenset[str]:  # noqa: ANN001
    if not name or name not in reg.names():
        return frozenset()
    return reg.get(name).capabilities.accepted_extras


def _polisher_request(name: str, row: ReadRow) -> str:
    """A polisher request as it applies to one run: 'none' when nothing would polish it."""
    reads = read_set(row)
    if name == "none" or not accepting_polishers(polisher_registry, reads):
        return "none"
    if name != "auto" and name in polisher_registry.names():
        cls = polisher_registry.get(name)
        if not cls.__new__(cls).accepts(reads):
            return "none"
    return name


def _setting_changes(done: dict, params: AssembleParams, row: ReadRow) -> list[str]:
    """How ``params`` differ from the settings a finished run was built with.

    A request that names the tool the marker records agrees with it, so
    ``--assembler skesa`` reuses a run that ``auto`` assembled with skesa.
    """
    old = done["settings"]
    changes = []
    if old["assembler"] != params.assembler and params.assembler != done["assembler"]:
        changes.append(f"assembler {old['assembler']} -> {params.assembler}")
    old_pol = _polisher_request(old["polisher"], row)
    new_pol = _polisher_request(params.polisher, row)
    if old_pol != new_pol and new_pol != (done.get("polisher") or "none"):
        changes.append(f"polisher {old['polisher']} -> {params.polisher}")
    elif done.get("polisher") and old["polish_rounds"] != params.polish_rounds:
        changes.append(f"polish_rounds {old['polish_rounds']} -> {params.polish_rounds}")
    keys = _accepted_extras(done["assembler"], registry) | _accepted_extras(
        done.get("polisher"), polisher_registry
    )
    old_extra = {k: str(v) for k, v in old["extra"].items() if k in keys}
    new_extra = {k: str(v) for k, v in params.extra.items() if k in keys}
    if old_extra != new_extra:
        changes.append(f"--tool-arg {old_extra} -> {new_extra}")
    return changes


def _reuse_finished(
    outcome: _Outcome,
    params: AssembleParams,
    run_dir: Path,
    check_settings: bool,
    logger: logging.Logger | None,
) -> bool:
    """Fill ``outcome`` from a finished run in ``run_dir``; False when it must be assembled."""
    contigs = run_dir / _CONTIGS_NAME
    done = _read_marker(run_dir / _DONE_MARKER)
    if done is None or not contigs.exists():
        return False
    row = outcome.row
    old_floor = None
    # Markers written before the settings were recorded are reused as they are.
    if check_settings and "settings" in done:
        changes = _setting_changes(done, params, row)
        old_floor = done["settings"]["min_contig_length"]
        if old_floor > params.min_contig_length:
            changes.append(f"min_contig_length {old_floor} -> {params.min_contig_length}")
        if changes:
            if logger is not None:
                logger.info(
                    "%s: finished with other settings (%s); assembling again",
                    row.run_accession,
                    "; ".join(changes),
                )
            return False
    outcome.assembler = done["assembler"]
    outcome.polisher = done.get("polisher") or None
    outcome.polish_rounds = done.get("polish_rounds", 0)
    outcome.stats = ContigStats(**done["stats"])
    outcome.tool_stats = done.get("tool_stats", {})
    outcome.versions = dict(done.get("tool_versions") or {})
    if not outcome.versions and done.get("version"):
        outcome.versions = {done["assembler"]: done["version"]}
    # A marker written before reads sketches existed has no such field.
    sketch_name = done.get("reads_sketch")
    outcome.reads_sketch = bool(sketch_name) and (run_dir / str(sketch_name)).is_file()
    # The screen decision the run was assembled under, when it was screened.
    outcome.screen = ScreenRecord.from_json(done.get("screen"))
    if check_settings and (old_floor is None or old_floor < params.min_contig_length):
        # A marker without settings may hold contigs below the floor; filtering
        # again is a no-op when it does not.
        _refilter(outcome, done, run_dir, params.min_contig_length, logger)
    return True


def _refilter(
    outcome: _Outcome,
    done: dict,
    run_dir: Path,
    min_length: int,
    logger: logging.Logger | None,
) -> None:
    """Apply a higher contig floor to finished contigs, which hold every contig
    above the old floor in assembly order, so the result equals filtering the
    raw assembly. A floor that keeps nothing excuses the run and leaves its
    files as they are, for a later, lower floor."""
    run = outcome.row.run_accession
    contigs = run_dir / _CONTIGS_NAME
    trial = run_dir / f"{_CONTIGS_NAME}.refilter"
    stats = filter_contigs(contigs, trial, min_length=min_length, prefix=run)
    if stats.n_contigs == 0:
        trial.unlink(missing_ok=True)
        outcome.stats = None
        outcome.excused = ExcusedRun(
            run, "assemble", f"assembly_failed: no contig of {min_length} bp or more"
        )
        return
    if "settings" in done:
        done["settings"]["min_contig_length"] = min_length
    elif stats == outcome.stats:
        trial.unlink()
        return
    # The marker goes before the contigs change: a kill in between leaves no
    # marker, so the run is assembled again rather than reused with wrong stats.
    (run_dir / _DONE_MARKER).unlink(missing_ok=True)
    trial.replace(contigs)
    done["stats"] = asdict(stats)
    _write_marker(run_dir / _DONE_MARKER, done)
    outcome.stats = stats
    if logger is not None:
        logger.info("%s: contigs filtered again at %d bp", run, min_length)


def _write_marker(path: Path, marker: dict) -> None:
    with atomic_replace(path) as fo:
        fo.write(json.dumps(marker, indent=1))


def _marker_versions(versions: dict[str, str], names: list[tuple[str, object]]) -> dict[str, str]:
    """The entries of ``versions`` that belong to the named adapters (by tool and binary)."""
    keys: set[str] = set()
    for name, reg in names:
        caps = reg.get(name).capabilities  # type: ignore[attr-defined]
        keys.add(caps.name)
        keys.update(b.name for b in caps.required_binaries)
    return {k: v for k, v in versions.items() if k in keys}


ASSEMBLER_NOT_INSTALLED = "assembler_not_installed"


def _excuse_missing_assemblers(plan: list[_Outcome], logger: logging.Logger) -> None:
    """Re-excuse runs an adapter accepts but whose tool is not installed.

    Under ``auto``, ``_plan`` excuses every run without an available adapter
    as ``unsupported_platform``. That is right when no adapter takes the
    platform; when one would, the tool is missing, so the run is excused as
    ``assembler_not_installed`` and one warning per platform names the
    adapters. The planner itself is unchanged because ``reads-gather`` calls
    it on hosts without assemblers.
    """
    needed: dict[str, list[str]] = {}
    for o in plan:
        if o.excused is None or o.excused.reason != "unsupported_platform":
            continue
        row = o.row
        reads = read_set(row)
        names = accepting_assemblers(registry, reads)
        if names:
            o.excused = ExcusedRun(row.run_accession, "assemble", ASSEMBLER_NOT_INSTALLED)
            needed.setdefault(row.platform, names)
    for platform, names in sorted(needed.items()):
        n = sum(
            1
            for o in plan
            if o.excused is not None
            and o.excused.reason == ASSEMBLER_NOT_INSTALLED
            and o.row.platform == platform
        )
        logger.warning(
            "%d %s run(s) excused as %s: none of %s is installed. Put one on PATH, run "
            "with --container, or narrow the runs with reads --platform, then rerun "
            "assemble with --force.",
            n,
            platform,
            ASSEMBLER_NOT_INSTALLED,
            ", ".join(names),
        )


def _warn_missing_polishers(plan: list[_Outcome], logger: logging.Logger) -> None:
    """Warn when ``--polisher auto`` leaves runs unpolished for want of a tool.

    A run to be assembled with no polisher chosen is either one no adapter
    takes (Illumina, HiFi) or one an adapter would take whose tool is not
    installed. The latter is still assembled, unpolished, and one warning per
    platform names the adapters, as for a missing assembler.
    """
    needed: dict[str, list[str]] = {}
    counts: dict[str, int] = {}
    for o in plan:
        if o.excused is not None or o.stats is not None or o.polisher is not None:
            continue
        row = o.row
        reads = read_set(row)
        names = accepting_polishers(polisher_registry, reads)
        if names:
            needed.setdefault(row.platform, names)
            counts[row.platform] = counts.get(row.platform, 0) + 1
    for platform, names in sorted(needed.items()):
        logger.warning(
            "%d %s run(s) will be assembled without polishing: none of %s is installed. "
            "Put one on PATH or run with --container to polish them, or pass "
            "--polisher none to assemble unpolished without this warning.",
            counts[platform],
            platform,
            ", ".join(names),
        )


DEFERRED = "deferred"


def _is_deferred(excused: ExcusedRun | None) -> bool:
    return excused is not None and excused.reason == DEFERRED


def defer_beyond(plan: list[_Outcome], max_runs: int, logger: logging.Logger) -> None:
    """Excuse as ``deferred`` the pending runs beyond ``max_runs``.

    ``max_runs`` bounds the runs with a finished assembly in the workdir, not
    the runs assembled in one call: finished runs (reused markers) count
    first and are always kept, and the remaining places go to pending runs in
    ``reads.tsv`` order (the reads stage writes the largest by bases first).
    Runs excused up front (no FASTQ mirror, no accepting or installed
    assembler) count toward nothing. The result is therefore the same for a
    repeated call with the same N, which the resume fingerprint relies on.
    """
    if max_runs < 1:
        raise UserInputError(f"--max-runs must be at least 1 (got {max_runs}).")
    n_finished = sum(1 for o in plan if o.stats is not None)
    pending = [o for o in plan if o.excused is None and o.stats is None]
    places = max(0, max_runs - n_finished)
    deferred = pending[places:]
    for o in deferred:
        o.assembler = None
        o.polisher = None
        o.excused = ExcusedRun(o.row.run_accession, "assemble", DEFERRED)
    if n_finished > max_runs:
        logger.info(
            "--max-runs %d: %d finished run(s) are kept beyond N%s",
            max_runs,
            n_finished,
            "; nothing is pending" if not pending else "",
        )
    if deferred:
        logger.info(
            "--max-runs %d: %d run(s) already assembled, %d to assemble, %d deferred",
            max_runs,
            n_finished,
            len(pending) - len(deferred),
            len(deferred),
        )


def _require_something_to_assemble(plan: list[_Outcome]) -> None:
    """Exit 4 when nothing can be assembled and a missing assembler is the cause."""
    missing = [
        o for o in plan if o.excused is not None and o.excused.reason == ASSEMBLER_NOT_INSTALLED
    ]
    if not missing or any(o.excused is None for o in plan):
        return
    by_platform: dict[str, list[str]] = {}
    for o in missing:
        reads = read_set(o.row)
        by_platform.setdefault(o.row.platform, accepting_assemblers(registry, reads))
    detail = "; ".join(
        f"{platform} runs need one of {', '.join(names)}"
        for platform, names in sorted(by_platform.items())
    )
    raise MissingBinaryError(
        f"No run can be assembled: no assembler is installed for --assembler auto ({detail}). "
        "Put one on PATH or run with --container."
    )


def _preflight(plan: list[_Outcome], logger: logging.Logger) -> dict[str, str]:
    """Check the assemblers the pending runs need; finished runs need none."""
    versions: dict[str, str] = {}
    pending = [o for o in plan if o.assembler and o.excused is None and o.stats is None]
    for name in sorted({o.assembler for o in pending if o.assembler}):
        versions.update(registry.create(name).preflight())
    for name in sorted({o.polisher for o in pending if o.polisher}):
        versions.update(polisher_registry.create(name).preflight())
    return versions


# --- one run --------------------------------------------------------------------


def _fetch_and_assemble(
    outcome: _Outcome,
    params: AssembleParams,
    threads: int,
    out_dir: Path,
    run_scratch: Path,
    versions: dict[str, str],
    logger: logging.Logger,
    *,
    reads_sketch: str | None = None,
    gate: ScreenGate | None = None,
) -> _Outcome:
    """Fetch and assemble one run into ``out_dir`` (contigs and the done marker).

    ``reads_sketch`` is the sourmash version when the run's reads are to be
    sketched, None otherwise. The sketch runs in a thread beside the
    assembler, on the verified FASTQ files, and ends before the reads are
    removed. The assembler keeps the run's full thread share and the
    single-threaded sketch runs alongside it, one thread over the share;
    giving that thread up slowed the assembler more than the sketch took.
    A failed sketch is a warning, never an excuse.

    With ``gate`` (``--screen-reads``) the sketch is awaited instead, and the
    reads screen judges the run before any assembler starts; a run it
    excuses keeps its reads sketch. A run whose kept sketch passed the
    screen while planning (``outcome.screen``) is not sketched again.
    """
    row = outcome.row
    sketch_path = out_dir / READS_SKETCH_NAME
    prior = outcome.screen if gate is not None else None
    keep_sketch = prior is not None and sketch_path.is_file()
    # The marker names finished contigs; a run assembled again has none until it ends.
    (out_dir / _DONE_MARKER).unlink(missing_ok=True)
    if not keep_sketch:
        sketch_path.unlink(missing_ok=True)
        _clear_screen(out_dir)
    # A hard exit can leave the temporary file of a sketch thread behind.
    clear_partials(out_dir)
    _clear_scratch(row, run_scratch)
    run_scratch.mkdir(parents=True, exist_ok=True)
    try:
        files = _fetch(row, run_scratch, logger)
    except RepGenRError as exc:
        logger.warning("%s: fetch failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "fetch", f"download_failed: {exc}")
        remove_tree(run_scratch)
        return outcome

    out_dir.mkdir(parents=True, exist_ok=True)
    sketch = None
    if reads_sketch is not None and not keep_sketch:
        sketch = ReadsSketchJob(row.run_accession, files, sketch_path, logger)
    screened_sketch = False
    if gate is not None:
        # Under the screen the sketch is awaited: the screen reads it.
        screened_sketch = keep_sketch or (sketch is not None and sketch.wait())
        sketch = None
        outcome.screen = gate.judge(row, sketch_path if screened_sketch else None, out_dir, prior)
        if not outcome.screen.passed:
            outcome.excused = ExcusedRun(row.run_accession, SCREEN_STEP, outcome.screen.reason)
            outcome.assembler = None
            outcome.polisher = None
            if not (params.keep_files or params.keep_reads):
                remove_tree(run_scratch)
            return outcome
    try:
        marker = _assemble_fetched(
            outcome, params, files, threads, out_dir, run_scratch, versions, logger
        )
    finally:
        # A stopped or failed assembler still waits for the sketch, which
        # reads the files the cleanup below removes.
        if sketch is not None:
            sketch.join()
    # Raises when a stop request cut the sketch short: no marker is written
    # and the reads stay, so a resume keeps them and assembles the run again.
    sketched = sketch.wait() if sketch is not None else screened_sketch
    if marker is None:  # excused
        sketch_path.unlink(missing_ok=True)
        if not params.keep_files:
            remove_tree(run_scratch)
        return outcome
    assert outcome.stats is not None
    outcome.reads_sketch = sketched
    marker.update(reads_sketch_fields(reads_sketch if sketched else None))
    if outcome.screen is not None:
        marker["screen"] = outcome.screen.to_json()
    _write_marker(out_dir / _DONE_MARKER, marker)
    if params.keep_files:
        pass
    elif params.keep_reads:
        remove_tree(run_scratch / "asm")
    else:
        remove_tree(run_scratch)
    logger.info(
        "%s: %d contigs, %d bp, N50 %d (%s)",
        row.run_accession,
        outcome.stats.n_contigs,
        outcome.stats.total_length,
        outcome.stats.n50,
        outcome.assembler,
    )
    return outcome


def _assemble_fetched(
    outcome: _Outcome,
    params: AssembleParams,
    files: tuple[Path, ...],
    threads: int,
    out_dir: Path,
    run_scratch: Path,
    versions: dict[str, str],
    logger: logging.Logger,
) -> dict | None:
    """Assemble, polish and filter fetched reads into ``out_dir``.

    Returns the done marker to write, or None when the run is excused (set on
    ``outcome``). Removes nothing: the caller cleans up once the reads sketch
    has ended.
    """
    row = outcome.row
    assert outcome.assembler is not None
    adapter = registry.create(outcome.assembler)
    reads = read_set(row, files)
    try:
        result = adapter.assemble(
            reads,
            run_scratch / "asm",
            AdapterParams(threads=threads, memory_gb=params.memory_gb, extra=dict(params.extra)),
            logger,
        )
    except RepGenRError as exc:
        logger.warning("%s: assembly failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "assemble", f"assembly_failed: {_why(exc)}")
        return None
    contigs_in = result.contigs
    polish_stats: dict = {}
    if outcome.polisher is not None:
        try:
            polished = polisher_registry.create(outcome.polisher).polish(
                reads,
                result.contigs,
                run_scratch / "polish",
                PolishParams(
                    threads=threads, rounds=params.polish_rounds, extra=dict(params.extra)
                ),
                logger,
            )
        except RepGenRError as exc:
            logger.warning("%s: polishing failed (%s)", row.run_accession, exc)
            outcome.excused = ExcusedRun(
                row.run_accession, "assemble", f"polish_failed: {_why(exc)}"
            )
            return None
        contigs_in = polished.contigs
        outcome.polish_rounds = polished.rounds
        polish_stats = dict(polished.tool_stats)
    try:
        stats = filter_contigs(
            contigs_in,
            out_dir / _CONTIGS_NAME,
            min_length=params.min_contig_length,
            prefix=row.run_accession,
        )
    except RepGenRError as exc:
        logger.warning("%s: contig filtering failed (%s)", row.run_accession, exc)
        outcome.excused = ExcusedRun(row.run_accession, "assemble", f"assembly_failed: {_why(exc)}")
        return None
    if stats.n_contigs == 0:
        outcome.excused = ExcusedRun(
            row.run_accession,
            "assemble",
            f"assembly_failed: no contig of {params.min_contig_length} bp or more",
        )
        (out_dir / _CONTIGS_NAME).unlink(missing_ok=True)
        return None

    outcome.stats = stats
    outcome.tool_stats = dict(result.tool_stats)
    tools: list[tuple[str, object]] = [(outcome.assembler, registry)]
    if outcome.polisher is not None:
        tools.append((outcome.polisher, polisher_registry))
    return {
        "assembler": outcome.assembler,
        "version": versions.get(adapter.capabilities.name, ""),
        "polisher": outcome.polisher or "",
        "polish_rounds": outcome.polish_rounds,
        "polish_stats": polish_stats,
        "stats": asdict(stats),
        "tool_stats": outcome.tool_stats,
        "tool_versions": _marker_versions(versions, tools),
        "settings": _settings(params),
    }


# --- reads sketch -----------------------------------------------------------------


def reads_sketch_fields(version: str | None) -> dict[str, str | None]:
    """The done-marker fields of a reads sketch; all None when there is none.

    ``reads_sketch`` is the file name within the run directory.
    """
    if version is None:
        return {"reads_sketch": None, "reads_sketch_params": None, "reads_sketch_version": None}
    return {
        "reads_sketch": READS_SKETCH_NAME,
        "reads_sketch_params": READS_SKETCH_PARAMS,
        "reads_sketch_version": version,
    }


class ReadsSketchJob:
    """One sourmash call on a run's FASTQ files, in a thread beside the assembler.

    :meth:`wait` joins the thread and reports the outcome: an INFO line with
    the time taken, or a WARNING when the sketch failed. A failure never
    propagates; the temporary file is removed by :func:`sketch_reads`.
    """

    def __init__(self, run: str, files: Sequence[Path], out: Path, logger: logging.Logger) -> None:
        self.run = run
        self.out = out
        self.seconds = 0.0
        self._logger = logger
        self._error: Exception | None = None
        self._thread = threading.Thread(
            target=self._work, args=(tuple(files),), name=f"reads-sketch-{run}", daemon=True
        )
        self._thread.start()

    def _work(self, files: tuple[Path, ...]) -> None:
        start = time.monotonic()
        try:
            sketch_reads(files, self.run, self.out, self._logger)
        except Exception as exc:  # a sketch never fails the assembly
            self._error = exc
        self.seconds = time.monotonic() - start

    def done(self) -> bool:
        return not self._thread.is_alive()

    def join(self) -> None:
        self._thread.join()

    def wait(self) -> bool:
        """Join the thread; True when the sketch was written.

        A sketch cut short by a stop request (SIGTERM or SIGINT) is raised
        rather than reported: the caller must not record the run as finished
        without a sketch, nor remove its reads, since a resume can then sketch
        them again.
        """
        self._thread.join()
        if self._error is not None and process.stop_requested.is_set():
            raise self._error
        if self._error is not None:
            self._logger.warning(
                "%s: reads sketch not written (%s); the assembly is not affected",
                self.run,
                _why(self._error) if isinstance(self._error, RepGenRError) else self._error,
            )
            return False
        self._logger.info("%s: reads sketch written in %.1f s", self.run, self.seconds)
        return True


def _kept_reads(row: ReadRow, run_scratch: Path) -> tuple[Path, ...] | None:
    """The run's FASTQ files left in scratch (--keep-reads), checked against
    their md5 when known; None unless every file is there and matches."""
    if not row.fastq_urls:
        return None
    md5s = list(row.fastq_md5) + [""] * (len(row.fastq_urls) - len(row.fastq_md5))
    files = []
    for url, md5 in zip(row.fastq_urls, md5s, strict=True):
        path = run_scratch / Path(url).name
        if not path.is_file() or (md5 and not _verified(path, md5)):
            return None
        files.append(path)
    return tuple(files)


def refill_reads_sketches(
    plan: list[_Outcome],
    assemblies: Path,
    scratch: Path,
    version: str,
    threads: int,
    logger: logging.Logger,
) -> int:
    """Sketch the reads of finished runs whose marker records no reads sketch.

    Applies to runs whose FASTQ files were kept (``--keep-reads``); the run is
    not assembled again. Runs without kept files keep a null sketch, and one
    INFO line counts them. Returns the number of sketches written.
    """
    missing = [o for o in plan if o.stats is not None and o.excused is None and not o.reads_sketch]
    if not missing:
        return 0
    todo = []
    for o in missing:
        run_dir = assemblies / o.row.run_accession
        clear_partials(run_dir)
        done = _read_marker(run_dir / _DONE_MARKER)
        if done is not None and done.get("reads_sketch"):
            # The marker names a sketch file that is gone: make it say so.
            done.update(reads_sketch_fields(None))
            _write_marker(run_dir / _DONE_MARKER, done)
        files = _kept_reads(o.row, scratch / o.row.run_accession)
        if files is not None:
            todo.append((o, files))
    if len(todo) < len(missing):
        logger.info(
            "%d finished run(s) have no reads sketch and their FASTQ files were not kept "
            "(--keep-reads); their reads sketch stays empty.",
            len(missing) - len(todo),
        )
    if not todo:
        return 0
    logger.info("Sketching the kept reads of %d finished run(s)", len(todo))

    def one(item: tuple[_Outcome, tuple[Path, ...]]) -> bool:
        o, files = item
        run_dir = assemblies / o.row.run_accession
        job = ReadsSketchJob(o.row.run_accession, files, run_dir / READS_SKETCH_NAME, logger)
        if not job.wait():
            return False
        done = _read_marker(run_dir / _DONE_MARKER)
        if done is None:
            (run_dir / READS_SKETCH_NAME).unlink(missing_ok=True)
            return False
        done.update(reads_sketch_fields(version))
        _write_marker(run_dir / _DONE_MARKER, done)
        o.reads_sketch = True
        return True

    results = parallel_map(one, todo, max(1, min(threads, len(todo))), logger=logger)
    return sum(1 for r in results if r)


def drop_reads_sketch(outcome: _Outcome, run_dir: Path) -> None:
    """Remove the reads sketch of an excused run and clear it in its marker."""
    outcome.reads_sketch = False
    sketch = run_dir / READS_SKETCH_NAME
    if not sketch.exists():
        return
    done = _read_marker(run_dir / _DONE_MARKER)
    if done is not None and done.get("reads_sketch"):
        done.update(reads_sketch_fields(None))
        _write_marker(run_dir / _DONE_MARKER, done)
    sketch.unlink(missing_ok=True)


def keeps_reads_sketch(excused: ExcusedRun | None) -> bool:
    """Whether a run keeps its reads sketch at the end of the stage.

    Accepted runs do; of the excused runs, a deferred one (not judged) and one
    the reads screen excused (the sketch shows why, and a later call with other
    thresholds screens it again without fetching).
    """
    return excused is None or _is_deferred(excused) or excused.step == SCREEN_STEP


# --- reads screen -------------------------------------------------------------------


def screen_settings(params: AssembleParams) -> ScreenSettings | None:
    """The reads screen's settings, or None without ``--screen-reads``.

    Refuses, before anything is fetched, a screen without the GTDB sketch and
    lineages (flag or environment variable) or with the reads sketch turned
    off (UserInputError), and one without sourmash (MissingBinaryError).
    """
    if not params.screen_reads:
        return None
    if params.reads_sketch is False:
        raise UserInputError(
            "--screen-reads screens the reads sketch; it cannot be combined with --no-reads-sketch."
        )
    return resolve_settings(
        params.screen_min_fraction,
        params.screen_dup_ani,
        params.gtdb_sketch,
        params.gtdb_lineages,
        command="--screen-reads",
        min_fraction_flag="--screen-min-fraction",
        dup_ani_flag="--screen-dup-ani",
    )


def _clear_screen(run_dir: Path) -> None:
    """Remove a run's screen decision and gather output from an earlier attempt."""
    (run_dir / SCREEN_JSON).unlink(missing_ok=True)
    if (run_dir / SCREEN_DIR).is_dir():
        remove_tree(run_dir / SCREEN_DIR)


def _apply_screen(o: _Outcome, record: ScreenRecord) -> None:
    o.screen = record
    if not record.passed:
        # A finished run keeps its contigs and marker for a later call.
        o.stats = None
        o.assembler = None
        o.polisher = None
        o.excused = ExcusedRun(o.row.run_accession, SCREEN_STEP, record.reason)


def screen_planned(
    plan: list[_Outcome],
    gate: ScreenGate,
    assemblies: Path,
    threads: int,
    logger: logging.Logger,
) -> None:
    """Screen, before any fetch, the runs whose reads sketch is already here.

    Finished runs are screened from their reads sketch, or keep a decision
    made with the same settings; their duplicate check runs in ``reads.tsv``
    order among them, and those accepted are the first runs the pending ones
    are compared with. A finished run without a reads sketch is not screened.
    A pending run with a kept sketch (screened out by an earlier call) keeps
    a decision made with the same settings; with other settings it is
    screened again from the sketch, so it is fetched only when it now passes
    the taxon and fraction checks (the duplicate check follows in the worker).
    """
    settings = gate.settings
    finished: list[tuple[_Outcome, ScreenRecord | None]] = []
    kept: list[tuple[_Outcome, ScreenRecord | None]] = []
    unscreened = 0
    for o in plan:
        if o.excused is not None:
            continue
        run_dir = assemblies / o.row.run_accession
        record = read_screen_json(run_dir)
        has_sketch = (run_dir / READS_SKETCH_NAME).is_file()
        reusable = (
            record is not None
            and record.made_with(settings)
            and not record.reason.startswith(SCREEN_FAILED)
        )
        if o.stats is not None:
            # A pass is reused once its duplicate check was made; an excuse as it is.
            if (
                reusable
                and record is not None
                and (not record.passed or (has_sketch and record.duplicates_checked))
            ):
                finished.append((o, record))
            elif o.reads_sketch and has_sketch:
                finished.append((o, None))
            else:
                o.screen = None
                unscreened += 1
        elif has_sketch and record is not None:
            kept.append((o, record if reusable else None))
    todo = [o for o, rec in [*finished, *kept] if rec is None]
    if todo:
        logger.info("Screening the kept reads sketches of %d run(s)", len(todo))

    def classify(o: _Outcome) -> tuple[str, ScreenRecord]:
        run_dir = assemblies / o.row.run_accession
        return o.row.run_accession, gate.classify(
            o.row, run_dir / READS_SKETCH_NAME, run_dir / SCREEN_DIR
        )

    workers = max(1, min(threads, len(todo)))
    fresh = dict(parallel_map(classify, todo, workers, logger=logger))
    for o, record in finished:
        run_dir = assemblies / o.row.run_accession
        if record is not None:
            if record.passed:
                gate.accept(o.row, run_dir / READS_SKETCH_NAME)
        else:
            record = gate.check_duplicate(
                o.row, fresh[o.row.run_accession], run_dir / READS_SKETCH_NAME
            )
            write_screen_json(run_dir, record)
            done = _read_marker(run_dir / _DONE_MARKER)
            if done is not None and record.passed:
                done["screen"] = record.to_json()
                _write_marker(run_dir / _DONE_MARKER, done)
        if not record.passed:
            logger.warning(
                "%s: finished run excused by the screen (%s)", o.row.run_accession, record.reason
            )
        _apply_screen(o, record)
    for o, record in kept:
        run_dir = assemblies / o.row.run_accession
        if record is None:
            record = fresh[o.row.run_accession]
            write_screen_json(run_dir, record)
        if record.passed:
            o.screen = record  # the worker adds the duplicate check
        else:
            _apply_screen(o, record)
    if unscreened:
        logger.info(
            "%d finished run(s) have no reads sketch and are not screened: their FASTQ "
            "files were not kept (--keep-reads), so they cannot be sketched without "
            "assembling them again.",
            unscreened,
        )


def _clear_scratch(row: ReadRow, run_scratch: Path) -> None:
    """Empty a run's scratch from an earlier attempt, except its FASTQ files.

    A run killed during assembly leaves its downloads behind; :func:`_fetch`
    keeps those that still match their checksum, so a resume does not fetch
    gigabytes again. Tool output and partial downloads are removed.
    """
    if not run_scratch.exists():
        return
    names = {Path(url).name for url in row.fastq_urls if url}
    for entry in run_scratch.iterdir():
        if entry.is_dir() and not entry.is_symlink():
            remove_tree(entry)
        elif entry.name not in names:
            entry.unlink(missing_ok=True)  # an AppleDouble twin goes with its file


def _verified(path: Path, md5: str) -> bool:
    try:
        http.verify_md5(path, md5)
    except WorkdirError:
        return False
    return True


def _fetch(row: ReadRow, run_scratch: Path, logger: logging.Logger) -> tuple[Path, ...]:
    """Bring the run's FASTQ files into scratch, verified when a checksum is known.

    A download is checked against its md5 as it is written; a local copy is
    hashed after the copy. A file left by an interrupted attempt is kept when
    it matches its checksum.
    """
    files = []
    md5s = list(row.fastq_md5) + [""] * (len(row.fastq_urls) - len(row.fastq_md5))
    for url, md5 in zip(row.fastq_urls, md5s, strict=True):
        dest = run_scratch / Path(url).name
        source = Path(url)
        if dest.exists() and md5 and _verified(dest, md5):
            logger.info("%s: keeping %s from an earlier attempt", row.run_accession, dest.name)
            files.append(dest)
            continue
        if source.exists():
            shutil.copy2(source, dest)
            if md5:
                http.verify_md5(dest, md5)
        else:
            # The checksum is computed while the bytes are written.
            http.download(url, dest, logger=logger, md5=md5 or None)
        files.append(dest)
    return tuple(files)


# --- contract rows ----------------------------------------------------------------


def _taxa(row: ReadRow) -> tuple[str, str, str]:
    return row.family or "unknown", row.genus or "unknown", row.species or "unknown"


def named_links(contigs: dict[str, Path], link_dir: Path) -> dict[str, Path]:
    """Run accession -> a link ``<run>.fasta`` to its contigs (unique names for
    tools that key their reports by file name)."""
    if link_dir.exists():
        remove_tree(link_dir)
    link_dir.mkdir(parents=True)
    links = {}
    for run, path in contigs.items():
        link = link_dir / f"{run}.fasta"
        os.symlink(path.resolve(), link)
        links[run] = link
    return links


def _assembly_sketches(
    named: dict[str, Path],
    classifier: str | None,
    params: AssembleParams,
    assemblies: Path,
    logger: logging.Logger,
) -> tuple[dict[Path, Path] | None, dict[str, SketchRecord]]:
    """Sketch each assembly once, beside its contigs, for the classifier and sketches/.

    Applies when the classifier reads genome sketches at parameters they hold
    and genome sketches are not turned off (``--no-sketch``). Returns the
    classifier's query sketches (named link -> file) and the record of each
    sketch by run accession, for the copy into ``sketches/`` at the end of
    the stage. A sketch whose assembly is unchanged is reused from an earlier
    run. An assembly whose sketch fails is left to the classifier.
    """
    records: dict[str, SketchRecord] = {}
    if classifier is None or not named or params.sketch is False:
        return None, records
    adapter = classifier_registry.create(classifier)

    def provide(links: Sequence[Path]) -> dict[Path, Path]:
        def one(link: Path) -> tuple[Path, Path | None]:
            run = record_name(link)
            sig = assemblies / run / _ASSEMBLY_SKETCH
            try:
                records[run] = sketch_beside(link, run, sig, logger)
            except (ToolExecutionError, WorkdirError, OSError) as exc:
                logger.warning("%s: assembly sketch not written (%s)", run, exc)
                return link, None
            return link, sig

        workers = max(1, min(params.threads, len(links)))
        done = parallel_map(one, list(links), workers, logger=logger)
        return {link: sig for link, sig in done if sig is not None}

    sketches = adapter_sketches(
        adapter, params.extra, list(named.values()), provide, logger, "assemble"
    )
    return sketches, records


def classifier_for(classifier: str, gtdb_sketch: str | None) -> str | None:
    """The classifier to run: explicit, or sourmash when a sketch is configured."""
    if classifier == "none":
        return None
    sketch = gtdb_sketch or os.environ.get(GTDB_SKETCH_ENV)
    if classifier == "auto":
        return "sourmash" if sketch else None
    if not sketch:
        raise UserInputError(
            f"--classifier {classifier} needs a reference sketch (--gtdb-sketch or "
            f"{GTDB_SKETCH_ENV})."
        )
    return classifier


def _check_outgroup_file(outgroup: str | None) -> None:
    """Refuse an --outgroup that is not a FASTA file, before any assembly.

    phylo, tree2tax and doctor resolve the outgroup among the FASTA files under
    outgroup/ only, so a file staged under another name would leave the tree
    unrooted.
    """
    if outgroup is None:
        return
    source = Path(outgroup).expanduser()
    if not source.is_file():
        raise UserInputError(f"--outgroup {outgroup} is not a file.")
    if not source.name.endswith(FASTA_SUFFIXES):
        raise UserInputError(
            f"--outgroup {outgroup} has no FASTA suffix; name it with one of "
            f"{', '.join(FASTA_SUFFIXES)}."
        )


def precheck(ctx: WorkdirContext, params: AssembleParams) -> None:
    """Refuse a wrong database path or a missing QC tool before the harness marks
    a finished record incomplete (registered in the CLI's stage prechecks)."""
    _check_outgroup_file(params.outgroup)
    require_sourmash_if_requested(params.sketch)
    require_sourmash_if_requested(params.reads_sketch)
    screen_settings(params)
    check_quality_inputs(
        checkm2_db=params.checkm2_db or checkm2_db_from_env(),
        classifier=classifier_for(params.classifier, params.gtdb_sketch),
        gtdb_sketch=params.gtdb_sketch,
        gtdb_lineages=params.gtdb_lineages,
    )


def check_quality_inputs(
    *,
    checkm2_db: str | None,
    classifier: str | None,
    gtdb_sketch: str | None,
    gtdb_lineages: str | None,
) -> dict[str, str]:
    """Refuse missing databases and find the QC tools; return their versions."""
    versions: dict[str, str] = {}
    # The environment variables stand in for the flags, as in classifier_for().
    gtdb_sketch = gtdb_sketch or os.environ.get(GTDB_SKETCH_ENV)
    gtdb_lineages = gtdb_lineages or os.environ.get(GTDB_LINEAGES_ENV)
    if checkm2_db:
        if not Path(checkm2_db).expanduser().exists():
            raise UserInputError(
                f"--checkm2-db {checkm2_db} does not exist (from the flag or CHECKM2DB)."
            )
        versions.update(preflight_checkm2())
    if classifier:
        if not gtdb_sketch:
            raise UserInputError(
                f"--classifier {classifier} needs a reference sketch (--gtdb-sketch or "
                f"{GTDB_SKETCH_ENV})."
            )
        if not Path(gtdb_sketch).expanduser().exists():
            raise UserInputError(
                f"--gtdb-sketch {gtdb_sketch} does not exist (from the flag or {GTDB_SKETCH_ENV})."
            )
        if gtdb_lineages and not Path(gtdb_lineages).expanduser().exists():
            raise UserInputError(
                f"--gtdb-lineages {gtdb_lineages} does not exist (from the flag or "
                f"{GTDB_LINEAGES_ENV})."
            )
        cls = classifier_registry.get(classifier)
        if cls.needs_lineages and not gtdb_lineages:
            raise UserInputError(
                f"The {classifier} classifier needs the lineages CSV published with the GTDB "
                f"sketch (--gtdb-lineages or {GTDB_LINEAGES_ENV})."
            )
        versions.update(classifier_registry.create(classifier).preflight())
    return versions


def assess(
    named: dict[str, Path],
    *,
    checkm2_db: str | None,
    classifier: str | None,
    gtdb_sketch: str | None,
    gtdb_lineages: str | None,
    threads: int,
    extra: dict[str, str],
    scratch: Path,
    versions: dict[str, str],
    logger: logging.Logger,
    memory_gb: float | None = None,
    cache_dirs: dict[str, Path] | None = None,
    sketches: Mapping[Path, Path] | None = None,
) -> tuple[dict[str, tuple[float, float]] | None, dict[str, Classification] | None]:
    """Score (CheckM2) and classify the assemblies in ``named`` (run -> FASTA).

    Either result is None when the corresponding database is not configured.
    Both are keyed by run accession. With ``cache_dirs`` (run -> directory),
    CheckM2 scores are stored per run and reused while the contigs, the
    database and the CheckM2 version are unchanged. ``memory_gb`` bounds the
    classifier's concurrent work. ``sketches`` (named FASTA -> sourmash
    signature file) are the classifier's query sketches, when it reads them.
    """
    quality: dict[str, tuple[float, float]] | None = None
    classified: dict[str, Classification] | None = None
    if checkm2_db and named:
        checkm2_version = preflight_checkm2()
        versions.update(checkm2_version)
        keys: dict[str, dict[str, object]] = {}
        quality = {}
        for run, link in named.items():
            if cache_dirs is None or run not in cache_dirs:
                continue
            keys[run] = quality_cache_key(link, checkm2_db, checkm2_version.get("checkm2", ""))
            stored = load_cached_quality(cache_dirs[run] / CHECKM2_CACHE, keys[run])
            if stored is not None:
                quality[run] = stored
        to_score = {run: link for run, link in named.items() if run not in quality}
        if quality:
            logger.info(
                "CheckM2: reusing the stored scores of %d of %d assemblies (same contigs, "
                "database and CheckM2 version)",
                len(quality),
                len(named),
            )
        by_name = (
            run_checkm2(
                list(to_score.values()),
                scratch / "checkm2",
                db=Path(checkm2_db),
                threads=threads,
                logger=logger,
            )
            if to_score
            else {}
        )
        for run, link in to_score.items():
            if link.name not in by_name:
                logger.warning(
                    "%s: CheckM2 reported no quality; the assembly is kept without quality values",
                    run,
                )
                continue
            quality[run] = by_name[link.name]
            if run in keys and cache_dirs is not None:
                store_cached_quality(cache_dirs[run] / CHECKM2_CACHE, keys[run], quality[run])
    elif not checkm2_db:
        logger.info(
            "No CheckM2 database configured (--checkm2-db or %s); assemblies are not "
            "quality-scored and --keeper quality will fall back to the tool's pick.",
            "CHECKM2DB",
        )
    if classifier and named:
        adapter = classifier_registry.create(classifier)
        versions.update(adapter.preflight())
        sketch = Path(gtdb_sketch or os.environ[GTDB_SKETCH_ENV])
        lineages = gtdb_lineages or os.environ.get(GTDB_LINEAGES_ENV)
        by_name_cls = adapter.classify(
            list(named.values()),
            scratch / "classify",
            ClassifyParams(
                db=sketch,
                lineages=None if lineages is None else Path(lineages),
                threads=threads,
                memory_gb=memory_gb,
                extra=dict(extra),
                sketches=sketches,
            ),
            logger,
        )
        classified = {
            run: by_name_cls[link.name] for run, link in named.items() if link.name in by_name_cls
        }
    return quality, classified


def apply_quality(
    outcomes: list[_Outcome],
    quality: dict[str, tuple[float, float]],
    min_completeness: float,
    max_contamination: float,
    logger: logging.Logger,
) -> None:
    """Attach CheckM2 values and excuse assemblies outside the gate."""
    for o in outcomes:
        o.quality = quality.get(o.row.run_accession)
        if o.quality is None:
            continue
        completeness, contamination = o.quality
        if completeness < min_completeness or contamination > max_contamination:
            o.excused = ExcusedRun(
                o.row.run_accession,
                "qc",
                f"qc_failed: completeness {completeness:.1f} "
                f"(min {min_completeness:g}), contamination {contamination:.1f} "
                f"(max {max_contamination:g})",
            )
            logger.warning("%s: excused, %s", o.row.run_accession, o.excused.reason)


CLASSIFIER_DISAGREES = "classifier_disagrees"


def apply_classification(
    outcomes: list[_Outcome],
    classified: dict[str, Classification],
    versions: dict[str, str],
    logger: logging.Logger,
) -> int:
    """Name by GTDB tokens where the classifier agrees at genus; flag the rest.

    A genome whose GTDB genus differs while the family and the species
    epithet agree (an NCBI name GTDB has moved to another genus of the same
    family, such as Mycoplasmopsis arginini, GTDB Metamycoplasma arginini) is
    flagged ``genus_renamed``; any other difference is flagged
    ``classifier_disagrees``. Both are warned about and keep the submitted
    name, and neither excuses the genome. The family is required because an
    epithet alone is shared across unrelated genera (Klebsiella pneumoniae,
    Streptococcus pneumoniae).

    Returns the number of disagreements (genus renames not counted).
    """
    n_disagree = 0
    for o in outcomes:
        o.gtdb = classified.get(o.row.run_accession)
        if o.gtdb is None:
            continue
        versions["gtdb_sketch"] = o.gtdb.db_version or versions.get("gtdb_sketch", "")
        gtdb_tokens = _gtdb_tokens(o.gtdb.taxonomy)
        assert o.label is not None
        gtdb_genus = next((c for c in o.gtdb.taxonomy.split(";") if c.startswith("g__")), "g__?")
        agreement = genus_agreement(o.label, gtdb_tokens)
        if agreement == AGREE:
            o.label = gtdb_tokens
            o.label_source = "classifier"
        elif agreement == GENUS_RENAMED:
            o.taxonomy_flag = GENUS_RENAMED
            logger.warning(
                "%s: submitted as %s %s; GTDB places the species in %s of the same family "
                "(%s). Keeping the submitted name and flagging %s",
                o.row.run_accession,
                o.label[1],
                o.label[2],
                gtdb_genus,
                o.gtdb.taxonomy.split(";")[-1],
                GENUS_RENAMED,
            )
        else:
            o.taxonomy_flag = CLASSIFIER_DISAGREES
            n_disagree += 1
            logger.warning(
                "%s: submitted as genus %s (%s) but classified as %s, %s; keeping the "
                "submitted name and flagging %s",
                o.row.run_accession,
                o.label[1],
                o.row.organism,
                gtdb_genus,
                o.gtdb.taxonomy.split(";")[-1],
                CLASSIFIER_DISAGREES,
            )
    return n_disagree


def _name(o: _Outcome) -> str:
    assert o.label is not None
    return genome_filename(*o.label, o.row.run_accession)


def _selection_row(o: _Outcome) -> SelectionRow:
    assert o.label is not None
    family, genus, species = o.label
    completeness, contamination = o.quality if o.quality else (None, None)
    return SelectionRow(
        o.row.run_accession,
        family,
        genus,
        species,
        False,
        _name(o),
        completeness,
        contamination,
    )


def _stats_row(o: _Outcome) -> AssemblyStatsRow:
    assert o.stats is not None and o.label is not None
    coverage = o.row.bases / o.stats.total_length if o.stats.total_length else None
    completeness, contamination = o.quality if o.quality else (None, None)
    return AssemblyStatsRow(
        run_accession=o.row.run_accession,
        filename=_name(o),
        assembler=o.assembler or "",
        n_contigs=o.stats.n_contigs,
        total_length=o.stats.total_length,
        n50=o.stats.n50,
        largest_contig=o.stats.largest_contig,
        est_coverage=None if coverage is None else round(coverage, 2),
        completeness=completeness,
        contamination=contamination,
        ncbi_taxonomy=";".join(_taxa(o.row)),
        gtdb_taxonomy=o.gtdb.taxonomy if o.gtdb else "",
        label_source=o.label_source,
        taxonomy_flag=o.taxonomy_flag,
        polisher=o.polisher or "",
        reads_sketch=o.reads_sketch,
        screen="pass" if o.screen is not None and o.screen.passed else "",
    )


def _append_rows(
    ctx: WorkdirContext, new_rows: list[SelectionRow], params: AssembleParams, logger
) -> list[SelectionRow]:
    """The existing selection with the new rows added (same accession replaced)."""
    from ..core.contracts import read_selection

    selection = ctx.workdir / SELECTION_TSV
    if not selection.exists():
        raise UserInputError(
            f"--append needs a working directory that already holds {SELECTION_TSV} "
            "(from metadata and genome, or ingest); run without --append to start one."
        )
    replaced = {r.accession for r in new_rows}
    kept = [r for r in read_selection(selection) if r.accession not in replaced]
    if params.outgroup is not None:
        kept = [r for r in kept if not r.is_outgroup]
        outgroup_row = _stage_outgroup(ctx, params.outgroup, logger)
        if outgroup_row is not None:
            kept.append(outgroup_row)
    logger.info("Appending %d assemblies to a selection of %d genomes", len(new_rows), len(kept))
    return [*kept, *new_rows]


# Excuses that say nothing about the run's data: it was never fetched, no
# assembler could be run on this host, --max-runs deferred it, or a tool of the
# reads screen failed. A stage in which only these occurred has not judged any
# run, so it does not discard an existing genome set. A run the screen excused
# for its reads (taxon_mismatch, low_match_fraction, duplicate_isolate) was judged.
_UNJUDGED_REASONS = (
    "download_failed",
    "no_fastq_mirror",
    ASSEMBLER_NOT_INSTALLED,
    DEFERRED,
    SCREEN_FAILED,
)


def _total_failure(
    ctx: WorkdirContext,
    params: AssembleParams,
    excused: list[ExcusedRun],
    logger: logging.Logger,
) -> str:
    """Decide what happens to the existing genome set when no run was accepted.

    The set is emptied only when ``assemble`` wrote it (every manifest row,
    outgroup included, has source ``sra``) and at least one run was judged
    (assembled, polished or quality-checked and rejected, or refused by the
    assembler). A set written by ``genome``, ``ingest`` or ``vgenome``, a set
    under ``--append``, or a stage in which every run failed to download is
    left as it is. Returns the clause the error message gives.
    """
    if params.append:
        return "the existing selection is unchanged"
    selection = ctx.workdir / SELECTION_TSV
    manifest_path = ctx.workdir / MANIFEST_FILENAME
    records = ctx.manifest.all_genomes(include_outgroup=True) if manifest_path.exists() else []
    selected = read_selection(selection) if selection.exists() else []
    files = (
        [p for p in ctx.genomes_dir.iterdir() if not p.name.startswith(".")]
        if ctx.genomes_dir.is_dir()
        else []
    )
    if not records and not selected and not files and not ctx.outgroup_dir.exists():
        return "no genome set was written"
    foreign = sorted({r.source or "unknown" for r in records if r.source != "sra"})
    if foreign or not records:
        origin = f"sources: {', '.join(foreign)}" if foreign else "not recorded in the manifest"
        return (
            f"the existing genome set of {len(selected) or len(records)} genome(s) "
            f"({origin}) was kept, since assemble did not write it"
        )
    if all(e.reason.startswith(_UNJUDGED_REASONS) for e in excused):
        return (
            f"the genome set of {len(records)} genome(s) from an earlier assemble call was "
            "kept, since no run could be fetched or assembled in this one; rerun with "
            "--force once the cause is resolved"
        )
    logger.warning(
        "Every run was excused: emptying %s, %s and the manifest (written by an earlier "
        "assemble call) so that later stages do not run on them.",
        ctx.genomes_dir.name,
        SELECTION_TSV,
    )
    _clear_genome_set(ctx, logger)
    return (
        "the previous genome set was cleared, so dereplicate refuses to run until assemble succeeds"
    )


def _clear_genome_set(ctx: WorkdirContext, logger: logging.Logger) -> None:
    """Empty the genome set an earlier ``assemble`` call wrote.

    The selection keeps its header, ``genomes/`` is emptied, the outgroup
    (staged by that call) is removed, and the per-run assemblies and their
    markers stay for a later call.
    """
    with staged_dir(ctx.genomes_dir):
        pass
    write_selection(ctx.workdir / SELECTION_TSV, [])
    (ctx.workdir / ASSEMBLY_STATS_TSV).unlink(missing_ok=True)
    _stage_outgroup(ctx, None, logger)
    ctx.manifest.replace_genomes([])
    try:
        remove_stale(ctx, logger)
    except OSError as exc:  # sketches are not deliverables; the set is already cleared
        logger.warning("Could not remove the sketches of the cleared genome set: %s", exc)


def _unlink_previous(genomes_dir: Path, accession: str) -> None:
    """Remove an earlier genome file of the same run (its label may have changed)."""
    for f in genomes_dir.iterdir():
        if f.is_file() and accession_from_filename(f.name) == accession:
            f.unlink()


def _stage_outgroup(
    ctx: WorkdirContext, outgroup: str | None, logger: logging.Logger
) -> SelectionRow | None:
    acc_file = ctx.workdir / OUTGROUP_ACCESSION_TXT
    if outgroup is None:
        if ctx.outgroup_dir.exists():
            remove_tree(ctx.outgroup_dir)
        acc_file.unlink(missing_ok=True)
        return None
    _check_outgroup_file(outgroup)
    source = Path(outgroup).expanduser()
    family, genus, species, accession = _parse(source.name)
    if ctx.outgroup_dir.exists():
        remove_tree(ctx.outgroup_dir)
    ctx.outgroup_dir.mkdir(parents=True)
    shutil.copy2(source, ctx.outgroup_dir / source.name)
    acc_file.write_text(accession + "\n", encoding="utf-8")
    logger.info("Outgroup: %s (%s)", accession, source.name)
    return SelectionRow(accession, family, genus, species, True, source.name)


def _parse(name: str) -> tuple[str, str, str, str]:
    from ..core.contracts import parse_genome_filename

    family, genus, species, accession = parse_genome_filename(name)
    return family, genus, species, accession or accession_from_filename(name) or name
