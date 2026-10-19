"""reads-screen stage: judge sketched runs against GTDB before any download.

Follows ``reads-sketch``: each selected run of ``reads.tsv`` that has a
sketch in one variant of ``reads_sketches/`` is screened with the three
checks of ``assemble --screen-reads`` (see :mod:`.assemble_screen`), so runs
of another taxon, runs whose reads match their species poorly and repeated
sequencing of one isolate can be left out before their FASTQ files are
fetched. The decisions go to ``reads_screen/<variant>/``: one row per run
in ``reads_screen.tsv``, each run's decision in ``<run>/screen.json`` with
the gather output beside it, and ``passed.txt``, the accepted accessions in
``reads.tsv`` order, in the form ``assemble --accession-file`` reads.

A decision made with the same settings on the same sketch is kept on a
repeat; a changed threshold screens the runs again from the kept sketches
without streaming anything. A run the tools failed on (``screen_failed``)
keeps no decision, so a repeat screens it again.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Collection
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    READS_SKETCHES_DIR,
    READS_SKETCHES_TSV,
    READS_TSV,
    ReadRow,
    atomic_replace,
    read_reads,
    reads_sketch_variant,
    reads_sketch_variant_spots,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.plugins import preflight
from ..core.resources import usable_cpus
from ..core.sourmash import SOURMASH_TOOL
from .assemble_screen import (
    DUPLICATE_ISOLATE,
    EXCUSED,
    PASS,
    SCREEN_DIR,
    SCREEN_FAILED,
    SCREEN_JSON,
    ScreenGate,
    ScreenRecord,
    ScreenSettings,
    read_screen_json,
    resolve_settings,
    write_screen_json,
    write_screen_table,
)
from .reads import read_accession_file
from .reads_sketch import DEFAULT_JOBS_CAP, is_current, read_index, sketch_file
from .taxon_match import GENUS_RENAMED

READS_SCREEN_DIR = "reads_screen"
READS_SCREEN_TSV = "reads_screen.tsv"
PASSED_TXT = "passed.txt"


@dataclass
class ReadsScreenParams:
    accessions: list[str] = field(default_factory=list)
    accession_file: str | None = None
    # The sketch variant: the cap of reads-sketch, or None for the one
    # present (``all`` when it is, else the only variant there).
    max_spots: int | None = None
    min_fraction: float = 0.3
    dup_ani: float = 0.999
    gtdb_sketch: str | None = None
    gtdb_lineages: str | None = None
    jobs: int | None = None


def default_jobs() -> int:
    """Gathers at once; each takes about 0.6 GB against the GTDB representatives."""
    return max(1, min(DEFAULT_JOBS_CAP, usable_cpus()))


def selected_accessions(params: ReadsScreenParams) -> list[str]:
    return [*params.accessions, *read_accession_file(params.accession_file)]


def select_runs(rows: list[ReadRow], params: ReadsScreenParams) -> list[ReadRow]:
    """The runs of ``reads.tsv`` the selection names, all of them without
    one; an accession absent from the table is refused by name."""
    by_run = {r.run_accession: r for r in rows}
    wanted = selected_accessions(params)
    if params.accession_file and not wanted:
        raise UserInputError(
            f"--accession-file {params.accession_file} lists no accessions; nothing selected."
        )
    unknown = [a for a in wanted if a not in by_run]
    if unknown:
        raise UserInputError(
            f"Not in {READS_TSV}: {', '.join(unknown)}. reads-screen takes runs the reads "
            "stage selected and reads-sketch sketched."
        )
    return [by_run[a] for a in dict.fromkeys(wanted)] if wanted else list(rows)


def resolve_variant(sketches: Path, max_spots: int | None) -> str:
    """The variant directory name to screen: the cap's, or without a cap
    ``all`` when present, else the only variant there (WorkdirError when
    none or several)."""
    if max_spots is not None:
        name = reads_sketch_variant(max_spots)
        if not (sketches / name / READS_SKETCHES_TSV).is_file():
            raise WorkdirError(
                f"No reads-sketch variant {name} under {sketches}; run "
                f"'repgenr reads-sketch --max-spots {max_spots}' first"
                + (f" (present: {', '.join(present)})." if (present := _present(sketches)) else ".")
            )
        return name
    present = _present(sketches)
    if not present:
        raise WorkdirError(f"No reads sketches under {sketches}; run 'repgenr reads-sketch' first.")
    if reads_sketch_variant(None) in present:
        return reads_sketch_variant(None)
    if len(present) == 1:
        return present[0]
    raise WorkdirError(
        f"Several reads-sketch variants are present ({', '.join(present)}); name one "
        "with --max-spots N."
    )


def _present(sketches: Path) -> list[str]:
    """The variant directories with an index, ``all`` first and then by rising cap."""
    if not sketches.is_dir():
        return []
    found = []
    for entry in sketches.iterdir():
        known, spots = reads_sketch_variant_spots(entry.name)
        if known and (entry / READS_SKETCHES_TSV).is_file():
            found.append((spots is not None, spots or 0, entry.name))
    return [name for _capped, _spots, name in sorted(found)]


def settings_of(params: ReadsScreenParams) -> ScreenSettings:
    return resolve_settings(
        params.min_fraction,
        params.dup_ani,
        params.gtdb_sketch,
        params.gtdb_lineages,
        command="reads-screen",
        min_fraction_flag="--min-fraction",
        dup_ani_flag="--dup-ani",
    )


def _require_reads(ctx: WorkdirContext) -> list[ReadRow]:
    reads = ctx.config.stages.get("reads")
    path = ctx.workdir / READS_TSV
    if reads is None or not reads.completed or not path.is_file():
        raise WorkdirError(
            f"reads-screen needs a completed reads stage ({READS_TSV}) in {ctx.workdir}; "
            "run 'repgenr reads' first."
        )
    return read_reads(path)


def precheck(ctx: WorkdirContext, params: ReadsScreenParams) -> None:
    """Refuse before the record is touched: no reads stage, no sketches, a
    bad selection, a missing GTDB database or sourmash."""
    select_runs(_require_reads(ctx), params)
    resolve_variant(ctx.workdir / READS_SKETCHES_DIR, params.max_spots)
    settings_of(params)
    preflight(SOURMASH_TOOL)


def run(ctx: WorkdirContext, params: ReadsScreenParams) -> int:
    """Screen the selected runs of one sketch variant; return the number of
    runs of ``reads.tsv`` that pass in that variant."""
    logger = ctx.logger
    rows = _require_reads(ctx)
    chosen = select_runs(rows, params)
    variant = resolve_variant(ctx.workdir / READS_SKETCHES_DIR, params.max_spots)
    _known, max_spots = reads_sketch_variant_spots(variant)
    sketch_dir = ctx.workdir / READS_SKETCHES_DIR / variant
    out_dir = ctx.workdir / READS_SCREEN_DIR / variant
    settings = settings_of(params)
    versions = preflight(SOURMASH_TOOL)
    if not chosen:
        logger.warning("No run of %s matches the selection; nothing to screen", READS_TSV)
    counts = screen_runs(
        chosen,
        sketch_dir,
        out_dir,
        settings,
        max_spots=max_spots,
        jobs=params.jobs or default_jobs(),
        force=ctx.force,
        logger=logger,
        all_rows=rows,
    )
    ctx.config.record_stage(
        "reads_screen",
        tool="sourmash",
        params={
            **asdict(params),
            "variant": variant,
            "selected": len(chosen),
            **counts,
            "settings": settings.as_dict(),
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Screens under %s: %s", READS_SCREEN_DIR, describe_screens(out_dir.parent))
    if counts["failed"]:
        logger.warning(
            "%d run(s) could not be screened (screen_failed: gather); a repeat of the same "
            "call screens them again",
            counts["failed"],
        )
    if counts["unsketched"]:
        logger.warning(
            "%d selected run(s) have no sketch in %s; sketch them with reads-sketch, then "
            "screen them",
            counts["unsketched"],
            sketch_dir.relative_to(ctx.workdir),
        )
    return counts["passed"]


# --- the work ----------------------------------------------------------------


def screen_runs(
    runs: list[ReadRow],
    sketch_dir: Path,
    out_dir: Path,
    settings: ScreenSettings,
    *,
    max_spots: int | None,
    jobs: int,
    force: bool,
    logger: logging.Logger,
    all_rows: list[ReadRow],
) -> dict[str, int]:
    """Screen ``runs`` from their sketches in ``sketch_dir`` into ``out_dir``.

    The taxon and fraction checks of a selected run are reused from a
    decision made with ``settings`` on its current sketch, unless ``force``;
    the others are gathered ``jobs`` at a time. The duplicate check is then
    made again on every call, in ``reads.tsv`` order over every run with a
    decision made with these settings (selected or not), so the result of
    several calls is that of one call over every run. Returns the counts of
    this call: ``screened`` (gathered now), ``reused``, ``unsketched``
    (selected without a sketch), and ``passed``, ``excused`` and ``failed``
    (screen_failed: gather) over every selected run.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    index = read_index(sketch_dir / READS_SKETCHES_TSV)
    gate = ScreenGate(
        settings, logger, duplicates=True, gather_slots=jobs, min_fraction_flag="--min-fraction"
    )
    selected = {r.run_accession for r in runs}
    # Every run with a decision in this call: a fresh or reused classification
    # for the selected runs, the kept decision for the others.
    records: dict[str, ScreenRecord] = {}
    unsketched: dict[str, ScreenRecord] = {}
    stale: list[str] = []
    todo: list[ReadRow] = []
    for row in all_rows:
        run = row.run_accession
        sketch = sketch_file(sketch_dir, run)
        if run in selected and not is_current(index.get(run), sketch_dir, max_spots):
            reason = index[run].reason if run in index and index[run].reason else "not sketched"
            unsketched[run] = ScreenRecord(
                run,
                EXCUSED,
                f"{SCREEN_FAILED}: no reads sketch ({reason})",
                settings=settings.as_dict(),
            )
            logger.warning("%s: no reads sketch in %s (%s); not screened", run, sketch_dir, reason)
            continue
        prior = None if force and run in selected else read_screen_json(out_dir / run)
        if prior is not None and _reusable(prior, settings, sketch, out_dir / run):
            records[run] = _before_duplicates(prior, row)
        elif run in selected:
            todo.append(row)
        elif prior is not None or (out_dir / run).is_dir():
            # Decided under other settings, or on an older sketch, and not
            # selected now: its row stays, it takes no part in the duplicate
            # check and it is left out of passed.txt.
            stale.append(run)
    if todo:
        logger.info(
            "Screening %d run(s) of %s against %s, %d at a time (min fraction %g, duplicate "
            "ANI %g)",
            len(todo),
            sketch_dir.name,
            Path(settings.gtdb_sketch).name,
            jobs,
            settings.min_fraction,
            settings.dup_ani,
        )

    def classify(row: ReadRow) -> tuple[str, ScreenRecord]:
        run = row.run_accession
        return run, gate.classify(row, sketch_file(sketch_dir, run), out_dir / run / SCREEN_DIR)

    records.update(parallel_map(classify, todo, max(1, min(jobs, len(todo))), logger=logger))
    reused = sum(1 for run in records if run in selected) - len(todo)
    # The duplicate check, in reads.tsv order over every decided run.
    for row in all_rows:
        run = row.run_accession
        if run not in records:
            continue
        record = records[run]
        if record.passed:
            record = gate.check_duplicate(row, record, sketch_file(sketch_dir, run))
            records[run] = record
        if _is_failed(record):
            (out_dir / run / SCREEN_JSON).unlink(missing_ok=True)
        else:
            write_screen_json(out_dir / run, record)
        if run in selected:
            log = logger.info if record.passed else logger.warning
            log(
                "%s: screen %s (top %s, fraction %s%s)",
                run,
                "passed" if record.passed else "excused the run",
                record.top_match or "-",
                "-" if record.fraction is None else f"{record.fraction:.3f}",
                f"; {record.reason}" if record.reason else "",
            )
    all_runs = [r.run_accession for r in all_rows]
    write_screen_table(
        out_dir / READS_SCREEN_TSV,
        [*records.values(), *unsketched.values()],
        all_runs,
        keep=all_runs,
    )
    passed = write_passed(out_dir, all_runs, exclude=stale)
    if stale:
        logger.warning(
            "%d run(s) of the table were decided with other settings or on an older sketch "
            "and are left out of %s: %s; screen them with this call's settings",
            len(stale),
            PASSED_TXT,
            ", ".join(stale[:5]) + (", ..." if len(stale) > 5 else ""),
        )
    mine = [records[run] for run in selected if run in records]
    counts = {
        "screened": len(todo),
        "reused": reused,
        "unsketched": len(unsketched),
        "stale": len(stale),
        "passed": sum(1 for r in mine if r.passed),
        "excused": sum(1 for r in mine if not r.passed and not _is_failed(r)),
        "failed": sum(1 for r in mine if _is_failed(r)),
    }
    logger.info(
        "Screened %d run(s), %d reused: %d passed, %d excused, %d failed, %d without a sketch; "
        "%d run(s) of %s pass in this variant",
        counts["screened"],
        counts["reused"],
        counts["passed"],
        counts["excused"],
        counts["failed"],
        counts["unsketched"],
        len(passed),
        READS_TSV,
    )
    return counts


def _is_failed(record: ScreenRecord) -> bool:
    return record.reason.startswith(SCREEN_FAILED)


def _reusable(prior: ScreenRecord, settings: ScreenSettings, sketch: Path, run_dir: Path) -> bool:
    """A decision made with these settings, after the sketch was written."""
    if not prior.made_with(settings) or _is_failed(prior):
        return False
    try:
        return (run_dir / SCREEN_JSON).stat().st_mtime_ns >= sketch.stat().st_mtime_ns
    except OSError:
        return False


def _before_duplicates(prior: ScreenRecord, row: ReadRow) -> ScreenRecord:
    """The taxon and fraction outcome of a kept decision, its duplicate
    excuse undone, so the duplicate check can be made again. A pass whose
    GTDB genus differs from the submitted one keeps its ``genus_renamed``
    note, which the duplicate excuse had replaced."""
    record = ScreenRecord(**prior.to_json())
    if record.reason.startswith(DUPLICATE_ISOLATE):
        record.decision = PASS
        record.reason = ""
        record.duplicate_of = ""
        record.ani = None
        if row.genus and record.top_genus and record.top_genus != row.genus:
            record.reason = GENUS_RENAMED
    record.duplicates_checked = False
    return record


def read_decisions(path: Path) -> dict[str, dict[str, str]]:
    """Run accession -> its row of a reads_screen.tsv, by column name."""
    if not path.is_file():
        return {}
    with open(path, encoding="utf-8", newline="") as fo:
        return {rec["run_accession"]: rec for rec in csv.DictReader(fo, delimiter="\t")}


def write_passed(out_dir: Path, order: list[str], *, exclude: Collection[str] = ()) -> list[str]:
    """``passed.txt``: the runs of ``order`` (reads.tsv) that passed, except
    ``exclude`` (runs decided with other settings)."""
    rows = read_decisions(out_dir / READS_SCREEN_TSV)
    passed = [
        run for run in order if run not in exclude and rows.get(run, {}).get("decision") == PASS
    ]
    with atomic_replace(out_dir / PASSED_TXT) as fo:
        fo.write("".join(f"{run}\n" for run in passed))
    return passed


# --- status ------------------------------------------------------------------


def screens(root: Path) -> list[tuple[str, int, int, int]]:
    """The screened variants under ``root``: (name, passed, excused, rows)."""
    if not root.is_dir():
        return []
    found = []
    for entry in root.iterdir():
        known, spots = reads_sketch_variant_spots(entry.name)
        if not known or not (entry / READS_SCREEN_TSV).is_file():
            continue
        rows = read_decisions(entry / READS_SCREEN_TSV).values()
        passed = sum(1 for rec in rows if rec.get("decision") == PASS)
        excused = sum(
            1
            for rec in rows
            if rec.get("decision") != PASS and not rec.get("reason", "").startswith(SCREEN_FAILED)
        )
        found.append((entry.name, spots, passed, excused, len(rows)))
    ordered = sorted(found, key=lambda v: (v[1] is not None, v[1] or 0))
    return [(name, passed, excused, n) for name, _s, passed, excused, n in ordered]


def describe_screens(root: Path) -> str:
    """``spots_100000 (3 passed, 1 excused of 5 runs)``, or ``none``."""
    parts = [
        f"{name} ({passed} passed, {excused} excused of {n} runs)"
        for name, passed, excused, n in screens(root)
    ]
    return ", ".join(parts) if parts else "none"
