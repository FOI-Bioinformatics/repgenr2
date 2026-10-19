"""reads-sketch stage: sketch sequencing runs streamed from the archive.

Follows ``reads``: the runs of ``reads.tsv``, all of them or a selection, are
streamed from NCBI with ``fastq-dump -Z`` and piped into one ``sourmash
sketch dna`` call each, so no FASTQ file is written and a run can be cut
after its first spots (``--max-spots``). Each run leaves
``reads_sketches/<run>.sig.zip``, a sketch in the form the assemble stage
writes for a run's reads (k=21, 31 and 51 at scaled=1000, with abundance),
and ``reads_sketches/reads_sketches.tsv`` lists every run of ``reads.tsv``
with its sketch status. A call adds to the set: runs sketched before with
the same cap are kept.

The same work without a workdir (``-o DIR`` with bare accessions) is
:func:`run_standalone`; nothing is recorded there.
"""

from __future__ import annotations

import logging
import re
import shutil
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core import process
from ..core.containers import run_pipeline
from ..core.context import WorkdirContext
from ..core.contracts import (
    READS_SKETCHES_DIR,
    READS_SKETCHES_TSV,
    READS_TSV,
    ReadRow,
    ReadsSketchRow,
    read_reads,
    read_reads_sketches,
    write_reads_sketches,
)
from ..core.errors import ToolExecutionError, UserInputError, WorkdirError
from ..core.executors import parallel_each
from ..core.plugins import preflight
from ..core.resources import usable_cpus
from ..core.sketches import READS_SKETCH_PARAMS, SKETCH_SUFFIX, clear_partials, partial_path
from ..core.sourmash import SOURMASH_TOOL
from ..core.sratools import SRATOOLS_TOOL, fastq_dump_command
from .reads import _ENA_PLATFORM, PLATFORMS, read_accession_file

DEFAULT_JOBS_CAP = 4
# Hidden per-run working directories under reads_sketches/ (see sketch_one).
_HOME_PREFIX = ".home."


@dataclass
class ReadsSketchParams:
    accessions: list[str] = field(default_factory=list)
    accession_file: str | None = None
    platform: str = "any"  # any | illumina | ont | pacbio
    max_spots: int | None = None
    jobs: int | None = None


def default_jobs() -> int:
    """Runs streamed at once: each is two single-threaded processes bound by
    the network more than by the CPU."""
    return max(1, min(DEFAULT_JOBS_CAP, usable_cpus()))


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
    """Sketch the selected runs into ``reads_sketches/``; return the number
    of runs sketched in this call."""
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
    out_dir = ctx.workdir / READS_SKETCHES_DIR
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
            f"see {READS_SKETCHES_DIR}/{READS_SKETCHES_TSV} and the log."
        )
    ctx.config.record_stage(
        "reads_sketch",
        tool="sourmash",
        params={
            **asdict(params),
            "selected": len(chosen),
            "sketched": sketched,
            "failed": failed,
            "present": sum(1 for r in index.values() if r.status == "ok"),
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
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
) -> int:
    """``-o DIR`` without a workdir: sketch bare accessions, record nothing."""
    if not accessions:
        raise UserInputError("Without -wd, name the runs with --accession or --accession-file.")
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
    index = read_index(out_dir / READS_SKETCHES_TSV)
    sketched, failed, _versions = sketch_runs(
        rows,
        out_dir,
        index,
        max_spots=max_spots,
        jobs=jobs or default_jobs(),
        force=force,
        logger=logger,
        all_runs=[r.run_accession for r in rows],
        source="the request",
    )
    if sketched == 0 and failed == len(rows):
        raise WorkdirError(f"Every one of the {failed} run(s) failed to stream; see {out_dir}.")
    return sketched


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
    parameters, and so need not be streamed again."""
    return (
        row is not None
        and row.status == "ok"
        and row.max_spots == max_spots
        and row.sketch_params == READS_SKETCH_PARAMS
        and sketch_file(out_dir, row.run_accession).is_file()
    )


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
    logger.info(
        "Streaming and sketching %d run(s)%s, %d at a time (%s)",
        len(todo),
        "" if max_spots is None else f", the first {max_spots} spots of each",
        min(jobs, len(todo)),
        READS_SKETCH_PARAMS,
    )
    counts = {"ok": 0, "failed": 0}

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
        "Sketched %d run(s), %d failed; %d of %d run(s) of %s now have a sketch",
        counts["ok"],
        counts["failed"],
        sum(1 for r in index.values() if r.status == "ok"),
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
_SEQUENCES_RE = re.compile(r"calculated \d+ signatures? for (\d+) sequences")


def spots_read(lines: list[str]) -> int | None:
    """The spot count ``fastq-dump`` reports (``Read N spots for ACC``)."""
    for line in lines:
        m = _SPOTS_RE.match(line.strip())
        if m:
            return int(m.group(1))
    return None


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
    """Stream one run through sourmash; the row records the outcome.

    The sketch goes to a hidden temporary file beside the target and is
    renamed into place only when both tools exited cleanly, so a failure
    leaves no partial sketch and a sketch from an earlier call is kept.
    """
    target = sketch_file(out_dir, run)
    tmp = partial_path(target)
    started = time.monotonic()
    # The working directory of the two tools is a disposable hidden directory:
    # under a container backend it is also HOME, where fastq-dump writes its
    # .ncbi settings and cache, which must not land beside the sketches.
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=_HOME_PREFIX) as home:
        try:
            result = run_pipeline(
                (SRATOOLS_TOOL, fastq_dump_command(run, max_spots)),
                (SOURMASH_TOOL, sourmash_stdin_command(run, tmp)),
                logger=logger,
                cwd=home,
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
    if not result.ok or not tmp.is_file():
        tmp.unlink(missing_ok=True)
        return ReadsSketchRow(
            run,
            "failed",
            max_spots=max_spots,
            spots_read=spots_read(dump_tail),
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
