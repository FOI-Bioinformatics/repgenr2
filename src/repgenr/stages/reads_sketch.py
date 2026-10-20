"""reads-sketch stage: sketch sequencing runs streamed from the archive.

Follows ``reads``: the runs of ``reads.tsv``, all of them or a selection, are
streamed from NCBI with ``fastq-dump -Z`` and piped into one ``sourmash
sketch dna`` call each, so no FASTQ file is written and a run can be cut
after its first spots (``--max-spots``). Each cap is a variant directory,
``reads_sketches/all/`` for the whole run and ``reads_sketches/spots_N/`` for
a cap, so sketches at several depths coexist. Each run leaves
``<variant>/<run>.sig.zip``, a sketch in the form the assemble stage writes
for a run's reads (k=21, 31 and 51 at scaled=1000, with abundance), and
``<variant>/reads_sketches.tsv`` lists every run of ``reads.tsv`` with its
sketch status. A call adds to its variant: runs sketched before with the
same cap are kept.

The same work without a workdir (``-o DIR`` with bare accessions) is
:func:`run_standalone`; nothing is recorded there.
"""

from __future__ import annotations

import fcntl
import logging
import os
import random
import re
import shutil
import tempfile
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core import process
from ..core.containers import engine_failure_line, ensure_image, run_pipeline
from ..core.context import WorkdirContext
from ..core.contracts import (
    READS_SKETCHES_DIR,
    READS_SKETCHES_TSV,
    READS_TSV,
    ReadRow,
    ReadsSketchRow,
    read_reads,
    read_reads_sketches,
    reads_sketch_variant,
    reads_sketch_variant_spots,
    write_reads_sketches,
)
from ..core.errors import ToolExecutionError, UserInputError, WorkdirError
from ..core.executors import parallel_each
from ..core.plugins import preflight
from ..core.resources import usable_cpus
from ..core.sketches import READS_SKETCH_PARAMS, SKETCH_SUFFIX, clear_partials, partial_path
from ..core.sourmash import SKETCH_TOOL_ENV, SOURMASH_TOOL
from ..core.sratools import SRATOOLS_TOOL, fastq_dump_command
from ..core.versions import write_versions_fragment
from .reads import _ENA_PLATFORM, PLATFORMS, read_accession_file

DEFAULT_JOBS_CAP = 4
# Hidden per-run working directories inside a variant directory (see sketch_one).
_HOME_PREFIX = ".home."
# The lock file of a variant directory: one call writes a variant at a time.
LOCK_NAME = ".lock"
# Attempts per run and the pauses before the second and third (seconds, with
# a quarter of jitter so retries of a batch do not land together). Only a
# failure on the archive side is tried again; see retry_wanted.
ATTEMPTS = 3
RETRY_DELAYS = (30.0, 120.0)
# Per-run time limit with a cap: an allowance for the resolution of the run
# at NCBI (about 10 s measured) and a stream rate well below the slowest
# seen (200000 spots in about 100 s); a stalled connection then frees its
# job slot and is tried again. Without a cap the whole run streams and only
# REPGENR_SUBPROCESS_TIMEOUT applies.
TIMEOUT_BASE = 300.0
TIMEOUT_SPOTS_PER_SECOND = 200.0


@dataclass
class ReadsSketchParams:
    accessions: list[str] = field(default_factory=list)
    accession_file: str | None = None
    platform: str = "any"  # any | illumina | ont | pacbio
    max_spots: int | None = None
    jobs: int | None = None


def default_jobs() -> int:
    """Runs streamed at once: each is one stream from the archive and one
    single-threaded sourmash process, so a job costs about one CPU."""
    return max(1, min(DEFAULT_JOBS_CAP, usable_cpus()))


def run_timeout(max_spots: int | None) -> float | None:
    """The time limit of one run's pipeline; None without a cap."""
    if max_spots is None:
        return None
    return TIMEOUT_BASE + max_spots / TIMEOUT_SPOTS_PER_SECOND


def selected_accessions(params: ReadsSketchParams) -> list[str]:
    """The accessions named on the command line and in the file, in order."""
    return [*params.accessions, *read_accession_file(params.accession_file)]


def select_runs(rows: list[ReadRow], params: ReadsSketchParams) -> list[ReadRow]:
    """The runs of ``reads.tsv`` the selection names; all of them without one.

    An accession absent from ``reads.tsv`` is refused by name, so a typo does
    not quietly sketch a run outside the selection ``reads`` made.
    """
    by_run = {r.run_accession: r for r in rows}
    wanted = selected_accessions(params)
    unknown = [a for a in wanted if a not in by_run]
    if unknown:
        raise UserInputError(
            f"Not in {READS_TSV}: {', '.join(unknown)}. reads-sketch takes runs the reads "
            "stage selected; run 'repgenr reads' with these accessions first."
        )
    chosen = [by_run[a] for a in dict.fromkeys(wanted)] if wanted else list(rows)
    platform = _ENA_PLATFORM.get(params.platform)
    if platform is not None:
        chosen = [r for r in chosen if r.platform == platform]
    return chosen


def precheck(ctx: WorkdirContext, params: ReadsSketchParams) -> None:
    """Refuse before the record is touched: no reads stage, a bad selection."""
    if params.platform not in PLATFORMS:
        raise UserInputError(f"--platform must be one of {', '.join(PLATFORMS)}.")
    reads = ctx.config.stages.get("reads")
    path = ctx.workdir / READS_TSV
    if reads is None or not reads.completed or not path.is_file():
        raise WorkdirError(
            f"reads-sketch needs a completed reads stage ({READS_TSV}) in {ctx.workdir}; "
            "run 'repgenr reads' first."
        )
    select_runs(read_reads(path), params)


def run(ctx: WorkdirContext, params: ReadsSketchParams) -> int:
    """Sketch the selected runs into ``reads_sketches/<variant>/``; return the
    number of runs sketched in this call."""
    logger = ctx.logger
    precheck(ctx, params)
    rows = read_reads(ctx.workdir / READS_TSV)
    chosen = select_runs(rows, params)
    if not chosen:
        logger.warning(
            "No run of %s matches the selection (platform %s); nothing to sketch",
            READS_TSV,
            params.platform,
        )
    out_dir = ctx.workdir / READS_SKETCHES_DIR / reads_sketch_variant(params.max_spots)
    index = read_index(out_dir / READS_SKETCHES_TSV)
    sketched, failed, versions = sketch_runs(
        chosen,
        out_dir,
        index,
        max_spots=params.max_spots,
        jobs=params.jobs or default_jobs(),
        force=ctx.force,
        logger=logger,
        all_runs=[r.run_accession for r in rows],
    )
    if chosen and sketched == 0 and failed == len(chosen):
        # Raised before the record is written, as assemble does when no run
        # was accepted: the stage stays incomplete and a repeat runs again.
        raise WorkdirError(
            f"Every one of the {failed} selected run(s) failed to stream; "
            f"see {out_dir.relative_to(ctx.workdir)}/{READS_SKETCHES_TSV} and the log."
        )
    ctx.config.record_stage(
        "reads_sketch",
        tool="sourmash",
        params={
            **asdict(params),
            "variant": out_dir.name,
            "selected": len(chosen),
            "sketched": sketched,
            "failed": failed,
            "present": sketched_of(index, [r.run_accession for r in rows]),
            "variants": {name: ok for name, _spots, ok, _n in variants(out_dir.parent)},
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Variants under %s: %s", READS_SKETCHES_DIR, describe_variants(out_dir.parent))
    if failed:
        logger.warning(
            "%d run(s) failed; a repeat of the same call streams them again, since a run "
            "without a sketch is a missing deliverable",
            failed,
        )
    return sketched


def run_standalone(
    accessions: list[str],
    out_dir: Path,
    *,
    max_spots: int | None,
    jobs: int | None,
    force: bool,
    logger: logging.Logger,
    reads_tsv: Path | None = None,
    versions_out: Path | None = None,
) -> int:
    """``-o DIR`` without a workdir: sketch bare accessions, or the runs of a
    ``reads.tsv`` (every run, or the named subset), and record nothing.
    ``versions_out`` takes the tool versions as a YAML fragment (Nextflow)."""
    if reads_tsv is not None:
        if not reads_tsv.is_file():
            raise UserInputError(f"--reads-tsv {reads_tsv} is not a file.")
        rows = read_reads(reads_tsv)
        if accessions:
            by_run = {r.run_accession: r for r in rows}
            unknown = [a for a in accessions if a not in by_run]
            if unknown:
                raise UserInputError(f"Not in {reads_tsv}: {', '.join(unknown)}.")
            rows = [by_run[a] for a in dict.fromkeys(accessions)]
        source = str(reads_tsv)
    else:
        if not accessions:
            raise UserInputError(
                "Without -wd, name the runs with --accession or --accession-file, or give "
                "--reads-tsv."
            )
        odd = [a for a in accessions if not _RUN_ACCESSION.match(a)]
        if odd:
            raise UserInputError(
                f"Not a run accession (SRR, ERR or DRR followed by digits): {', '.join(odd)}"
            )
        rows = [
            ReadRow(
                run_accession=a,
                biosample="",
                bioproject="",
                organism="",
                taxid="",
                platform="",
                instrument_model="",
                layout="",
                bases=0,
                read_count=0,
            )
            for a in dict.fromkeys(accessions)
        ]
        source = "the request"
    out_dir = out_dir / reads_sketch_variant(max_spots)
    index = read_index(out_dir / READS_SKETCHES_TSV)
    sketched, failed, versions = sketch_runs(
        rows,
        out_dir,
        index,
        max_spots=max_spots,
        jobs=jobs or default_jobs(),
        force=force,
        logger=logger,
        all_runs=[r.run_accession for r in rows],
        source=source,
    )
    if versions_out is not None:
        write_versions_fragment(versions_out, versions)
    if rows and sketched == 0 and failed == len(rows):
        raise WorkdirError(f"Every one of the {failed} run(s) failed to stream; see {out_dir}.")
    logger.info("Variants under %s: %s", out_dir.parent, describe_variants(out_dir.parent))
    return sketched


def variants(root: Path) -> list[tuple[str, int | None, int, int]]:
    """The variant directories under ``root``: (name, cap, runs with a
    sketch, rows in the index), ``all`` first and then by rising cap."""
    if not root.is_dir():
        return []
    found = []
    for entry in root.iterdir():
        known, spots = reads_sketch_variant_spots(entry.name)
        if not known or not (entry / READS_SKETCHES_TSV).is_file():
            continue
        rows = read_reads_sketches(entry / READS_SKETCHES_TSV)
        ok = sum(1 for r in rows if r.status == "ok" and (entry / r.sketch_file).is_file())
        found.append((entry.name, spots, ok, len(rows)))
    return sorted(found, key=lambda v: (v[1] is not None, v[1] or 0))


def describe_variants(root: Path) -> str:
    """``all (3 of 5 runs sketched), spots_100000 (5 of 5 runs sketched)``, or ``none``."""
    parts = [f"{name} ({ok} of {n} runs sketched)" for name, _spots, ok, n in variants(root)]
    return ", ".join(parts) if parts else "none"


# The run accessions of SRA, ENA and DDBJ; the name is also the sketch file name.
_RUN_ACCESSION = re.compile(r"^[SED]RR[0-9]+$")


# --- the work ----------------------------------------------------------------


def read_index(path: Path) -> dict[str, ReadsSketchRow]:
    """The index rows by run, empty when there is none yet."""
    if not path.is_file():
        return {}
    return {r.run_accession: r for r in read_reads_sketches(path)}


def sketch_file(out_dir: Path, run: str) -> Path:
    return out_dir / f"{run}{SKETCH_SUFFIX}"


def is_current(row: ReadsSketchRow | None, out_dir: Path, max_spots: int | None) -> bool:
    """Whether a run's sketch exists, was made with this cap and these
    parameters, and so need not be streamed again. A run that streamed no
    reads at this cap (``empty``) is current too: streaming it again gives
    the same nothing."""
    if row is None or row.max_spots != max_spots:
        return False
    if row.status == "empty":
        return True
    return (
        row.status == "ok"
        and row.sketch_params == READS_SKETCH_PARAMS
        and sketch_file(out_dir, row.run_accession).is_file()
    )


def sketched_of(index: dict[str, ReadsSketchRow], runs: list[str]) -> int:
    """How many of ``runs`` have a sketch in ``index``."""
    return sum(1 for run in runs if (row := index.get(run)) is not None and row.status == "ok")


@contextmanager
def variant_lock(out_dir: Path) -> Iterator[None]:
    """Hold the variant for one call: a second call on the same directory
    would clear this one's partial sketches and working directories and
    overwrite its index, so it is refused instead."""
    path = out_dir / LOCK_NAME
    fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise WorkdirError(
                f"Another reads-sketch is writing {out_dir} (it holds {path.name}); wait "
                "for it to finish, or stop it, and run again."
            ) from None
        yield
    finally:
        os.close(fd)


def sketch_runs(
    runs: list[ReadRow],
    out_dir: Path,
    index: dict[str, ReadsSketchRow],
    *,
    max_spots: int | None,
    jobs: int,
    force: bool,
    logger: logging.Logger,
    all_runs: list[str],
    source: str = READS_TSV,
) -> tuple[int, int, dict[str, str]]:
    """Stream and sketch ``runs`` that are not current, ``jobs`` at a time.

    ``index`` is updated in place and the file rewritten after every finished
    run, in the order of ``all_runs``, so an interrupted call keeps what was
    done. Returns (sketched, failed, tool versions).
    """
    versions = {**preflight(SRATOOLS_TOOL), **preflight(SOURMASH_TOOL)}
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir.parent / READS_SKETCHES_TSV).is_file():
        logger.warning(
            "%s holds sketches of the earlier flat layout (an index beside the variant "
            "directories); they are not read. Remove them or move them under all/.",
            out_dir.parent,
        )
    with variant_lock(out_dir):
        return _sketch_runs(
            runs,
            out_dir,
            index,
            max_spots=max_spots,
            jobs=jobs,
            force=force,
            logger=logger,
            all_runs=all_runs,
            source=source,
            versions=versions,
        )


def _sketch_runs(
    runs: list[ReadRow],
    out_dir: Path,
    index: dict[str, ReadsSketchRow],
    *,
    max_spots: int | None,
    jobs: int,
    force: bool,
    logger: logging.Logger,
    all_runs: list[str],
    source: str,
    versions: dict[str, str],
) -> tuple[int, int, dict[str, str]]:
    clear_partials(out_dir)
    for entry in out_dir.iterdir():
        if entry.name.startswith(_HOME_PREFIX) and entry.is_dir():
            shutil.rmtree(entry, ignore_errors=True)  # left by a hard exit
    todo = [
        r for r in runs if force or not is_current(index.get(r.run_accession), out_dir, max_spots)
    ]
    skipped = len(runs) - len(todo)
    if skipped:
        logger.info(
            "%d of %d selected run(s) already sketched%s; skipping them (use --force to "
            "sketch them again)",
            skipped,
            len(runs),
            "" if max_spots is None else f" at --max-spots {max_spots}",
        )
    if not todo:
        write_index(out_dir / READS_SKETCHES_TSV, index, all_runs)
        return 0, 0, versions
    # Both images once, before the pool: the first runs would otherwise each
    # pull inside their measured time, and a pull failure is one error.
    ensure_image(SRATOOLS_TOOL, logger)
    ensure_image(SOURMASH_TOOL, logger)
    logger.info(
        "Streaming and sketching %d run(s)%s, %d at a time (%s)",
        len(todo),
        "" if max_spots is None else f", the first {max_spots} spots of each",
        min(jobs, len(todo)),
        READS_SKETCH_PARAMS,
    )
    counts = {"ok": 0, "failed": 0, "empty": 0}

    def one(row: ReadRow) -> ReadsSketchRow:
        return sketch_one(row.run_accession, out_dir, max_spots, logger)

    def record(_row: ReadRow, row: ReadsSketchRow) -> None:
        # On the calling thread, one result at a time (parallel_each).
        index[row.run_accession] = row
        counts[row.status] += 1
        write_index(out_dir / READS_SKETCHES_TSV, index, all_runs)
        if row.status == "ok":
            logger.info(
                "%s: sketched %s sequence(s) from %s spot(s) in %.1f s",
                row.run_accession,
                row.sequences if row.sequences is not None else "?",
                row.spots_read if row.spots_read is not None else "?",
                row.seconds or 0.0,
            )
        elif row.status == "empty":
            logger.warning("%s: no sketch, %s", row.run_accession, row.reason)
        else:
            logger.warning("%s: not sketched (%s)", row.run_accession, row.reason)

    parallel_each(
        one,
        todo,
        min(jobs, len(todo)),
        record,
        logger=logger,
        describe=lambda r: r.run_accession,
    )
    logger.info(
        "Sketched %d run(s), %d failed%s; %d of %d run(s) of %s now have a sketch",
        counts["ok"],
        counts["failed"],
        f", {counts['empty']} streamed no reads" if counts["empty"] else "",
        sketched_of(index, all_runs),
        len(all_runs),
        source,
    )
    return counts["ok"], counts["failed"], versions


def write_index(path: Path, index: dict[str, ReadsSketchRow], all_runs: list[str]) -> None:
    """One row per run of ``all_runs`` in that order; ``pending`` when unseen.
    Rows of runs no longer listed (standalone use, an edited table) follow."""
    rows = [index.get(run) or ReadsSketchRow(run, "pending") for run in all_runs]
    listed = set(all_runs)
    rows += [r for run, r in index.items() if run not in listed]
    write_reads_sketches(path, rows)


_SPOTS_RE = re.compile(r"^(?:Read|Written) (\d+) spots")
# sourmash's note when stdin held no sequence at all (it then exits 0).
_NO_SEQUENCES_RE = re.compile(r"no sequences found")
_SEQUENCES_RE = re.compile(r"calculated \d+ signatures? for (\d+) sequences")


def spots_read(lines: list[str]) -> int | None:
    """The spot count ``fastq-dump`` reports (``Read N spots for ACC``)."""
    for line in lines:
        m = _SPOTS_RE.match(line.strip())
        if m:
            return int(m.group(1))
    return None


def streamed_nothing(lines: list[str]) -> bool:
    """Whether sourmash reported an input without any sequence."""
    return any(_NO_SEQUENCES_RE.search(line) for line in lines)


def sequences_sketched(lines: list[str]) -> int | None:
    """The sequence count sourmash reports (``calculated 1 signature for N sequences``)."""
    for line in lines:
        m = _SEQUENCES_RE.search(line)
        if m:
            return int(m.group(1))
    return None


def sourmash_stdin_command(name: str, out: Path) -> list[str]:
    """``sourmash sketch dna`` over stdin, the reads-sketch parameters."""
    return [
        "sourmash",
        "sketch",
        "dna",
        "-p",
        READS_SKETCH_PARAMS,
        "--name",
        name,
        "-o",
        str(out),
        "-",
    ]


def sketch_one(
    run: str, out_dir: Path, max_spots: int | None, logger: logging.Logger
) -> ReadsSketchRow:
    """Stream one run through sourmash, up to ``ATTEMPTS`` times; the row
    records the outcome of the last attempt.

    A failure on the archive side that is not permanent (see
    :func:`retry_wanted`) is tried again after a pause; a failure of sourmash,
    an unknown or protected accession, and a run that streams no reads are
    not.
    """
    for attempt in range(1, ATTEMPTS + 1):
        row = _attempt(run, out_dir, max_spots, logger)
        if row.status != "failed" or attempt == ATTEMPTS or not retry_wanted(row.reason):
            return row
        pause = RETRY_DELAYS[min(attempt, len(RETRY_DELAYS)) - 1] * random.uniform(0.75, 1.25)  # noqa: S311
        logger.warning(
            "%s: %s; attempt %d of %d in %.0f s", run, row.reason, attempt + 1, ATTEMPTS, pause
        )
        if process.stop_requested.wait(pause):
            raise ToolExecutionError(
                ["reads-sketch", run], -1, output=process.NOT_STARTED, tool="reads-sketch"
            )
    raise RuntimeError("unreachable")


# Failures of fastq-dump a retry cannot fix: the run is unknown, withdrawn or
# protected (dbGaP), or has no data.
_PERMANENT = re.compile(
    r"could not be found|not found|no such|access denied|permission|protected|"
    r"not authorized|no data|invalid accession|withdrawn|suppressed",
    re.IGNORECASE,
)


def retry_wanted(reason: str) -> bool:
    """Whether a failed attempt's reason is the archive side's and transient:
    a fastq-dump failure or a timed-out pipeline, not one the message calls
    permanent."""
    if "[killed after" in reason:
        return True
    return reason.startswith("fastq-dump exited") and not _PERMANENT.search(reason)


def _attempt(
    run: str, out_dir: Path, max_spots: int | None, logger: logging.Logger
) -> ReadsSketchRow:
    """One pass of the pipeline over a run.

    The sketch goes to a hidden temporary file beside the target and is
    renamed into place only when both tools exited cleanly, so a failure
    leaves no partial sketch and a sketch from an earlier call is kept.
    """
    target = sketch_file(out_dir, run)
    tmp = partial_path(target)
    started = time.monotonic()
    # The working directory of the two tools is a disposable hidden directory:
    # under a container backend it is also HOME, where fastq-dump writes its
    # .ncbi settings and cache, which must not land beside the sketches. Its
    # cleanup is best effort: the next call removes what is left.
    with tempfile.TemporaryDirectory(
        dir=out_dir, prefix=f"{_HOME_PREFIX}{run}.", ignore_cleanup_errors=True
    ) as home:
        try:
            result = run_pipeline(
                (SRATOOLS_TOOL, fastq_dump_command(run, max_spots)),
                (SOURMASH_TOOL, sourmash_stdin_command(run, tmp)),
                logger=logger,
                cwd=home,
                timeout=run_timeout(max_spots),
                tool_envs=(None, SKETCH_TOOL_ENV),
            )
        except ToolExecutionError as exc:
            tmp.unlink(missing_ok=True)
            if process.stop_requested.is_set():
                raise
            return ReadsSketchRow(run, "failed", max_spots=max_spots, reason=str(exc))
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    seconds = time.monotonic() - started
    dump_tail, sketch_tail = result.tails
    if process.stop_requested.is_set():
        tmp.unlink(missing_ok=True)
        raise ToolExecutionError(
            ["reads-sketch", run], -1, output=process.NOT_STARTED, tool="reads-sketch"
        )
    engine = engine_failure_line(result.returncodes[0], dump_tail) or engine_failure_line(
        result.returncodes[1], sketch_tail
    )
    if engine:
        # The engine could not start a tool (no image, daemon down): not this
        # run's failure, and the batch stops on it instead of recording it
        # against every run.
        tmp.unlink(missing_ok=True)
        raise ToolExecutionError(
            ["reads-sketch", run],
            result.returncodes[0] or result.returncodes[1],
            output=engine,
            tool="reads-sketch",
            message=f"The container engine could not start a tool for {run}: {engine}",
        )
    spots = spots_read(dump_tail)
    if result.ok and not tmp.is_file() and (spots == 0 or streamed_nothing(sketch_tail)):
        # fastq-dump wrote no biological read (a run of technical reads only,
        # or no spots) and sourmash, given nothing, exits 0 without a file.
        return ReadsSketchRow(
            run,
            "empty",
            max_spots=max_spots,
            spots_read=spots if spots is not None else 0,
            seconds=seconds,
            reason="the stream held no biological reads (technical reads only, or no spots)",
        )
    if not result.ok or not tmp.is_file():
        tmp.unlink(missing_ok=True)
        return ReadsSketchRow(
            run,
            "failed",
            max_spots=max_spots,
            spots_read=spots,
            seconds=seconds,
            reason=failure_reason(result.returncodes, dump_tail, sketch_tail),
        )
    tmp.replace(target)
    return ReadsSketchRow(
        run,
        "ok",
        max_spots=max_spots,
        spots_read=spots_read(dump_tail),
        sequences=sequences_sketched(sketch_tail),
        sketch_file=target.name,  # relative to the index beside it
        sketch_params=READS_SKETCH_PARAMS,
        seconds=seconds,
    )


_GENERIC_EXIT = re.compile(r"quit with error code|exited with|^\s*$")
_ERROR_LINE = re.compile(r"err(or)?\b", re.IGNORECASE)


def last_message(lines: list[str]) -> str:
    """The most telling of a tool's last stderr lines: the last one that
    names an error, else the last line, skipping the generic exit notice
    fastq-dump ends with."""
    telling = [line for line in lines if not _GENERIC_EXIT.search(line)]
    errors = [line for line in telling if _ERROR_LINE.search(line)]
    if errors:
        return errors[-1].strip()
    if telling:
        return telling[-1].strip()
    return lines[-1].strip() if lines else "no output"


def failure_reason(
    returncodes: tuple[int, int], dump_tail: list[str], sketch_tail: list[str]
) -> str:
    """One line naming which tool failed and its last message."""
    dump_rc, sketch_rc = returncodes
    # A producer that failed on its own (an unknown accession, a network
    # error) says so on stderr and leaves the consumer with nothing to
    # sketch, so its failure is the cause. A producer that exited 1 without
    # a word, or ended by a signal, lost its reader: fastq-dump exits 1 with
    # an empty stderr when its pipe breaks, or is stopped once the consumer
    # failed, and which happens first is a race. That is the consumer's doing.
    if dump_rc != 0 and any(_DUMP_ERROR.search(line) for line in dump_tail):
        return f"fastq-dump exited {dump_rc}: {last_message(dump_tail)}"
    if sketch_rc != 0:
        return f"sourmash exited {sketch_rc}: {last_message(sketch_tail)}"
    if dump_rc != 0:
        return f"fastq-dump exited {dump_rc}: {last_message(dump_tail)}"
    return "no sketch was written"


# fastq-dump's own error lines: "fastq-dump.3.4.1 err: ..." and its closing
# "fastq-dump quit with error code N".
_DUMP_ERROR = re.compile(r"\berr: |quit with error code")
