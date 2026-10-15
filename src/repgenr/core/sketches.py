"""Genome sketches: one sourmash signature collection per genome of a workdir.

The genome-writing stages (genome, vgenome, ingest, assemble) and the
``repgenr sketch`` command leave ``sketches/<name>.sig.zip`` beside the genome
set, where ``<name>`` is the genome's record name (:func:`record_name`, the
name the aligners and tree builders give it). Each file is a sourmash zip
collection of three DNA signatures, k=21, k=31 and k=51 at scaled=1000 without
abundance, named after the genome, so a consumer can map a signature to its
genome without parsing paths.

Whether a sketch is current is decided from the sha256 of the genome FASTA and
the parameter string recorded with it, not from file times: the manifest holds
both for each genome row (``sketch_file``, ``sketch_params``,
``sketch_digest``); a genome under ``outgroup/`` without a manifest row is
recorded in ``sketches/outgroup.json`` instead.

Sketches are not stage deliverables. A genome set without them is complete;
the stages that compare genomes build their own sketches when these are absent.
Every sketch is written to a hidden temporary file and renamed into place, so
an interrupted run leaves no partial sketch, and the records are written in
batches as the sketches finish, so a rerun sketches only what is missing.
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
import tempfile
import threading
from collections.abc import Callable, Collection, Mapping, Sequence
from concurrent.futures import CancelledError, Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from . import process
from .containers import run_tool
from .contracts import (
    SELECTION_TSV,
    accession_from_filename,
    atomic_replace,
    list_fasta,
    read_selection,
    record_name,
)
from .errors import MissingBinaryError, RepGenRError, ToolExecutionError, WorkdirError
from .inputs import file_digest
from .manifest import SketchRecord
from .plugins import preflight, tool_available
from .resources import usable_cpus
from .sourmash import SOURMASH_TOOL

if TYPE_CHECKING:
    from .context import WorkdirContext

KSIZES: tuple[int, ...] = (21, 31, 51)
SCALED = 1000
# The parameter string passed to ``sourmash sketch dna -p`` and recorded with
# each sketch; a sketch recorded with another string is stale.
SKETCH_PARAMS = ",".join(f"k={k}" for k in KSIZES) + f",scaled={SCALED}"
SKETCHES_DIR = "sketches"
SKETCH_SUFFIX = ".sig.zip"
# The reads sketch of a sequencing run (assemble): the k-mer sizes and scale
# of the genome sketches with abundance tracking, since read k-mers carry
# coverage. It is kept beside the run's assembly and never in sketches/, which
# holds one genome sketch per genome.
READS_SKETCH_PARAMS = SKETCH_PARAMS + ",abund"
READS_SKETCH_NAME = "reads" + SKETCH_SUFFIX
# Sketch records of outgroup genomes that have no manifest row.
OUTGROUP_JSON = "outgroup.json"
# FASTA digests by (size, mtime_ns), so a consumer does not hash an unchanged
# genome again (see resolve_sketches). Hidden: not a sketch, and removed with
# the directory when the last sketch goes.
DIGESTS_JSON = ".digests.json"
# sourmash picks the output format from the file name, so the temporary file
# keeps the .sig.zip suffix; the leading dot hides it from every listing.
_PARTIAL = ".partial"
# The threads a writer stage without -t/--threads gives the sketch step.
_DEFAULT_THREADS = 16

ACTION_PRESENT = "present"
ACTION_WRITTEN = "written"
ACTION_REPLACED = "replaced"
ACTION_FORCED = "forced"
ACTION_COPIED = "copied"


@dataclass(frozen=True)
class SketchTarget:
    """A genome to sketch: its FASTA and where its sketch is recorded.

    ``accession`` names its manifest row; None records it in
    ``sketches/outgroup.json`` under the file name.
    """

    path: Path
    accession: str | None

    @property
    def name(self) -> str:
        return record_name(self.path)


@dataclass(frozen=True)
class SketchSource:
    """A sketch made elsewhere, reusable when its record matches the genome.

    ``name`` is the name its signatures carry when that differs from the
    genome's record name (the assemble stage sketches an assembly under its
    run accession before the genome is named); the copy is then renamed.
    """

    sketch: Path
    record: SketchRecord
    name: str | None = None


@dataclass
class SketchSummary:
    present: int = 0  # up to date, left as they were
    written: int = 0  # sketched for the first time
    replaced: int = 0  # a stale sketch replaced
    forced: int = 0  # a current sketch written again under --force
    copied: int = 0  # reused from another working directory
    removed: int = 0  # sketches of genomes no longer in the set

    def as_dict(self) -> dict[str, int]:
        return {
            "present": self.present,
            "written": self.written,
            "replaced": self.replaced,
            "forced": self.forced,
            "copied": self.copied,
            "removed": self.removed,
        }

    def line(self) -> str:
        forced = f", {self.forced} replaced (--force)" if self.forced else ""
        copied = f", {self.copied} copied" if self.copied else ""
        return (
            f"Sketches: {self.present} present, {self.written} written, "
            f"{self.replaced} stale replaced{forced}{copied}, {self.removed} removed"
        )


@dataclass
class SketchStatus:
    """Record names of the genomes by the state of their sketch."""

    present: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    stale: list[str] = field(default_factory=list)

    @property
    def total(self) -> int:
        return len(self.present) + len(self.missing) + len(self.stale)


def sketches_dir(workdir: Path) -> Path:
    return workdir / SKETCHES_DIR


def sketch_path(ctx: WorkdirContext, genome_path: Path) -> Path:
    """Where the sketch of ``genome_path`` lives in the workdir of ``ctx``."""
    return sketches_dir(ctx.workdir) / (record_name(genome_path) + SKETCH_SUFFIX)


def default_threads() -> int:
    """Sketch threads for a writer stage that takes no -t/--threads."""
    return max(1, min(_DEFAULT_THREADS, usable_cpus()))


def _relative(ctx: WorkdirContext, path: Path) -> str:
    return path.relative_to(ctx.workdir).as_posix()


def _row_paths(ctx: WorkdirContext) -> dict[str, Path]:
    """accession -> genome path of each manifest row, whether or not the file exists.

    The bacterial genome stage records no filename on its outgroup row; that
    row is matched to the FASTA under ``outgroup/`` that carries its accession.
    """
    outgroup_files = list_fasta(ctx.outgroup_dir)
    paths: dict[str, Path] = {}
    for row in ctx.manifest.all_genomes(include_outgroup=True):
        if row.filename:
            directory = ctx.outgroup_dir if row.is_outgroup else ctx.genomes_dir
            paths[row.accession] = directory / row.filename
        elif row.is_outgroup:
            match = [p for p in outgroup_files if accession_from_filename(p.name) == row.accession]
            if len(match) == 1:
                paths[row.accession] = match[0]
    return paths


def sketch_targets(ctx: WorkdirContext) -> list[SketchTarget]:
    """Every genome of the workdir whose FASTA is present, outgroup included.

    The genomes are the manifest rows with a file; a FASTA under
    ``outgroup/`` without a row is added with no accession.
    """
    targets = [
        SketchTarget(path, accession)
        for accession, path in _row_paths(ctx).items()
        if path.is_file()
    ]
    seen = {t.path for t in targets}
    targets += [SketchTarget(p, None) for p in list_fasta(ctx.outgroup_dir) if p not in seen]
    return sorted(targets, key=lambda t: (t.accession is None, t.path.name))


def _expected_names(ctx: WorkdirContext) -> set[str]:
    """Record names of every genome of the set, whether or not its file is present."""
    names = {record_name(path) for path in _row_paths(ctx).values()}
    names |= {record_name(p) for p in list_fasta(ctx.outgroup_dir)}
    return names


def _read_outgroup_json(directory: Path) -> dict[str, SketchRecord]:
    path = directory / OUTGROUP_JSON
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}  # an unreadable record: the outgroup is sketched again
    if not isinstance(data, dict):
        return {}
    out: dict[str, SketchRecord] = {}
    for name, rec in data.items():
        if isinstance(rec, dict):
            out[name] = SketchRecord(rec.get("file"), rec.get("params"), rec.get("digest"))
    return out


def _write_outgroup_json(directory: Path, records: Mapping[str, SketchRecord]) -> None:
    path = directory / OUTGROUP_JSON
    if not records:
        path.unlink(missing_ok=True)
        return
    payload = {
        name: {"file": rec.file, "params": rec.params, "digest": rec.digest}
        for name, rec in sorted(records.items())
    }
    with atomic_replace(path) as fo:
        json.dump(payload, fo, indent=2, sort_keys=True)
        fo.write("\n")


def _record_key(target: SketchTarget) -> str:
    return target.accession if target.accession is not None else target.path.name


def _recorded(ctx: WorkdirContext) -> dict[str, SketchRecord]:
    """Recorded sketches by manifest accession, and by file name for outgroup.json."""
    records = dict(ctx.manifest.sketch_records())
    records.update(_read_outgroup_json(sketches_dir(ctx.workdir)))
    return records


def _is_current(rec: SketchRecord | None, rel: str, digest: str, out: Path) -> bool:
    return (
        rec is not None
        and rec.file == rel
        and rec.params == SKETCH_PARAMS
        and rec.digest == digest
        and out.is_file()
    )


def _partial(out: Path) -> Path:
    """A unique hidden temporary name beside ``out`` (``.<name>.partial.<random>.sig.zip``).

    Unique per call, so two processes writing the same sketch never share a
    temporary file; each renames its own into place. The name is reserved by
    mkstemp and the empty file removed again, because sourmash and a hard
    link both need to create the file themselves.
    """
    stem = out.name[: -len(SKETCH_SUFFIX)]
    fd, name = tempfile.mkstemp(dir=out.parent, prefix=f".{stem}{_PARTIAL}.", suffix=SKETCH_SUFFIX)
    os.close(fd)
    tmp = Path(name)
    tmp.unlink()
    return tmp


def _is_partial(name: str) -> bool:
    # The unique form, and the fixed .<name>.partial.sig.zip of earlier releases.
    return name.startswith(".") and _PARTIAL in name and name.endswith(SKETCH_SUFFIX)


def _clear_partials(directory: Path) -> None:
    """Remove temporary sketches an interrupted run left behind."""
    if not directory.is_dir():
        return
    for entry in directory.iterdir():
        if _is_partial(entry.name):
            entry.unlink(missing_ok=True)


def clear_partials(directory: Path) -> None:
    """Remove temporary sketch files left in ``directory`` (public form, for
    run directories outside ``sketches/``)."""
    _clear_partials(directory)


def sketch_command(genome: Path, name: str, out: Path) -> list[str | os.PathLike[str]]:
    """The one sourmash call that sketches a genome (three k-mer sizes, one file)."""
    return ["sourmash", "sketch", "dna", "-p", SKETCH_PARAMS, "--name", name, "-o", out, genome]


def reads_sketch_command(
    files: Sequence[Path], name: str, out: Path
) -> list[str | os.PathLike[str]]:
    """The one sourmash call that sketches a run's reads: all its FASTQ files,
    one signature per k-mer size, with abundances."""
    return [
        "sourmash",
        "sketch",
        "dna",
        "-p",
        READS_SKETCH_PARAMS,
        "--name",
        name,
        "-o",
        out,
        *files,
    ]


def sketch_reads(files: Sequence[Path], name: str, out: Path, logger: logging.Logger) -> None:
    """Sketch the FASTQ ``files`` of one run into ``out``, atomically.

    The files go to one sourmash call, so the reads of a pair give one
    signature per k-mer size. Gzipped files are read by sourmash itself.
    The temporary file is created beside ``out`` and renamed into place;
    a failure leaves no ``out``.
    """
    if not files:
        raise WorkdirError(f"{name}: no FASTQ files to sketch")
    tmp = _partial(out)
    try:
        run_tool(
            SOURMASH_TOOL,
            reads_sketch_command(files, name, tmp),
            logger=logger,
            log_prefix=f"sourmash {name}",
            # The reads may sit in a scratch directory outside the workdir.
            extra_mounts=sorted({str(Path(f).resolve().parent) for f in files}),
        )
        if not tmp.is_file():
            raise WorkdirError(f"sourmash wrote no reads sketch for {name}")
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def sourmash_for_reads_sketch(
    flag: bool | None, stage: str, logger: logging.Logger
) -> dict[str, str]:
    """Resolve --reads-sketch/--no-reads-sketch to sourmash's versions, or {}.

    ``flag`` None (the default) sketches when sourmash can run and otherwise
    logs one INFO line for the stage; True requires sourmash
    (MissingBinaryError otherwise); False never sketches.
    """
    if flag is False:
        return {}
    if flag is None and not tool_available(SOURMASH_TOOL):
        logger.info(
            "%s: sourmash not found; reads sketches not written (--reads-sketch requires them).",
            stage,
        )
        return {}
    try:
        return preflight(SOURMASH_TOOL)
    except MissingBinaryError as exc:
        if flag:
            raise
        logger.info("%s: reads sketches not written: %s", stage, exc)
        return {}


def sketch_file(genome: Path, name: str, out: Path, logger: logging.Logger) -> None:
    """Sketch one FASTA with the contract parameters into ``out``, atomically."""
    _sketch_one(SketchTarget(genome, None), out, logger, name=name)


def _sketch_one(
    target: SketchTarget, out: Path, logger: logging.Logger, *, name: str | None = None
) -> None:
    tmp = _partial(out)
    try:
        run_tool(
            SOURMASH_TOOL,
            sketch_command(target.path, name or target.name, tmp),
            logger=logger,
            log_prefix="sourmash",
            # A linked genome (ingest) lives elsewhere; its directory must be
            # visible inside a container too.
            extra_mounts=[str(target.path.resolve().parent)],
        )
        if not tmp.is_file():
            raise WorkdirError(f"sourmash wrote no sketch for {target.path.name}")
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _copy_one(source: Path, out: Path, logger: logging.Logger, rename: str | None = None) -> None:
    tmp = _partial(out)
    try:
        if rename is None:
            process.link_or_copy(source, tmp)
        else:
            # The signatures carry another name; a rename rewrites only the
            # zip, it does not read the genome again.
            run_tool(
                SOURMASH_TOOL,
                ["sourmash", "sig", "rename", source, rename, "-o", tmp],
                logger=logger,
                log_prefix="sourmash",
            )
            if not tmp.is_file():
                raise WorkdirError(f"sourmash wrote no renamed sketch for {rename}")
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def sketch_genomes(
    ctx: WorkdirContext,
    genomes: Sequence[SketchTarget],
    threads: int,
    logger: logging.Logger,
    *,
    force: bool = False,
    reuse: Mapping[str, SketchSource] | None = None,
    digests: _DigestCache | None = None,
) -> SketchSummary:
    """Sketch the genomes whose sketch is missing or stale, ``threads`` at a time.

    One sourmash process per genome. An up-to-date sketch (same FASTA sha256,
    same parameters) is left as it is unless ``force``. ``reuse`` maps record
    names to sketches of another workdir; one whose recorded digest and
    parameters match the genome is copied (or linked) instead of sketched.
    Each sketch is recorded on this thread as soon as it is in place, so an
    interrupted run keeps the records of every sketch it finished. After a
    failure no further genome is started; the running ones finish and are
    recorded, then the first failure is raised.
    """
    directory = sketches_dir(ctx.workdir)
    directory.mkdir(parents=True, exist_ok=True)
    _clear_partials(directory)
    records = _recorded(ctx)
    outgroup_records = _read_outgroup_json(directory)
    reuse = reuse or {}
    summary = SketchSummary()
    # Unchanged genomes (size, mtime) are not hashed again; --force hashes all.
    # A caller that hashed already passes its cache (``digests``) and saves it.
    cache = digests if digests is not None else _DigestCache(directory / DIGESTS_JSON)

    def work(target: SketchTarget) -> tuple[SketchTarget, str, SketchRecord]:
        out = sketch_path(ctx, target.path)
        rel = _relative(ctx, out)
        digest = cache.digest(target.path, refresh=force)
        fresh = SketchRecord(rel, SKETCH_PARAMS, digest)
        current = _is_current(records.get(_record_key(target)), rel, digest, out)
        if current and not force:
            return target, ACTION_PRESENT, fresh
        existed = out.is_file()
        source = None if force else reuse.get(target.name)
        if (
            source is not None
            and source.record.params == SKETCH_PARAMS
            and source.record.digest == digest
            and source.sketch.is_file()
        ):
            # A sketch named otherwise is renamed to the genome's record name.
            rename = target.name if source.name not in (None, target.name) else None
            _copy_one(source.sketch, out, logger, rename)
            return target, ACTION_COPIED, fresh
        _sketch_one(target, out, logger)
        if current:
            return target, ACTION_FORCED, fresh
        return target, (ACTION_REPLACED if existed else ACTION_WRITTEN), fresh

    def record(target: SketchTarget, action: str, rec: SketchRecord) -> None:
        setattr(summary, action, getattr(summary, action) + 1)
        if action == ACTION_PRESENT and records.get(_record_key(target)) == rec:
            return
        if target.accession is not None:
            ctx.manifest.set_sketches([(target.accession, rec)])
        else:
            # Records of genomes that left the set are pruned by remove_stale.
            outgroup_records[target.path.name] = rec
            _write_outgroup_json(directory, outgroup_records)

    workers = max(1, min(threads, len(genomes)))
    failures: list[Exception] = []
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        futures: list[Future[tuple[SketchTarget, str, SketchRecord]]] = [
            pool.submit(work, t) for t in genomes
        ]
        for future in as_completed(futures):
            try:
                target, action, rec = future.result()
            except CancelledError:
                continue  # not started after an earlier failure
            except Exception as exc:
                if not failures:
                    for other in futures:
                        other.cancel()  # queued genomes do not start
                failures.append(exc)
                continue
            record(target, action, rec)
    except BaseException:
        # A termination signal: queued genomes must not start; the running
        # ones settle (their tools are stopped by the handler).
        pool.shutdown(wait=True, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
    if digests is None:
        _save_cache(
            cache, [*(t.path for t in sketch_targets(ctx)), *(t.path for t in genomes)], logger
        )
    if failures:
        raise failures[0]
    return summary


def _save_cache(cache: _DigestCache, keep: Collection[Path], logger: logging.Logger) -> None:
    try:
        cache.save(keep)
    except OSError as exc:  # a cache, not a deliverable: the next run hashes again
        logger.debug("Digest cache not written: %s", exc)


def remove_stale(ctx: WorkdirContext, logger: logging.Logger | None = None) -> int:
    """Remove sketches of genomes no longer in the set; return how many.

    Also removes temporary files of an interrupted run and outgroup.json
    records of genomes that left the set or gained a manifest row, and the
    directory itself when nothing is left in it.
    """
    directory = sketches_dir(ctx.workdir)
    if not directory.is_dir():
        return 0
    _clear_partials(directory)
    expected = _expected_names(ctx)
    removed = 0
    for entry in directory.iterdir():
        if entry.name.startswith(".") or not entry.name.endswith(SKETCH_SUFFIX):
            continue
        if entry.is_dir() or entry.name[: -len(SKETCH_SUFFIX)] in expected:
            continue
        entry.unlink()
        removed += 1
        if logger is not None:
            logger.debug("Removed sketch %s (genome no longer in the set)", entry.name)
    unrowed = {p.name for p in list_fasta(ctx.outgroup_dir)} - {
        path.name for path in _row_paths(ctx).values()
    }
    kept = {n: r for n, r in _read_outgroup_json(directory).items() if n in unrowed}
    _write_outgroup_json(directory, kept)
    if not any(not e.name.startswith(".") for e in directory.iterdir()):
        for entry in directory.iterdir():  # dotfiles only (exFAT, Finder)
            if not entry.is_dir():
                entry.unlink(missing_ok=True)
        if not any(directory.iterdir()):
            directory.rmdir()
    if removed and logger is not None:
        logger.info("Removed %d sketch(es) of genomes no longer in the set", removed)
    return removed


def sketch_status(ctx: WorkdirContext, *, verify: bool = True) -> SketchStatus:
    """Present, missing and stale sketches of the genomes of the workdir.

    With ``verify`` each genome FASTA is hashed and compared with its record;
    without it a sketch counts as present when its file exists and its
    record carries the current parameters (no genome is read).
    """
    records = _recorded(ctx)
    status = SketchStatus()
    for target in sketch_targets(ctx):
        out = sketch_path(ctx, target.path)
        rel = _relative(ctx, out)
        rec = records.get(_record_key(target))
        if not out.is_file():
            status.missing.append(target.name)
        elif verify:
            current = _is_current(rec, rel, file_digest(target.path), out)
            (status.present if current else status.stale).append(target.name)
        elif rec is not None and rec.file == rel and rec.params == SKETCH_PARAMS:
            status.present.append(target.name)
        else:
            status.stale.append(target.name)
    return status


def sketch_counts(workdir: Path) -> tuple[int, int] | None:
    """(sketches present, genomes present) from selection.tsv; None without sketches/.

    Reads no genome and opens no manifest, for ``repgenr status``.
    """
    directory = sketches_dir(workdir)
    if not directory.is_dir():
        return None
    selection = workdir / SELECTION_TSV
    try:
        rows = read_selection(selection) if selection.is_file() else []
    except (OSError, ValueError, RepGenRError):
        return None
    genomes = [
        workdir / ("outgroup" if row.is_outgroup else "genomes") / row.filename for row in rows
    ]
    present = [g for g in genomes if g.is_file()]
    sketched = sum(1 for g in present if (directory / (record_name(g) + SKETCH_SUFFIX)).is_file())
    return sketched, len(present)


def expected_sketch_files(workdir: Path) -> list[Path]:
    """The sketch of every genome selection.tsv lists and the workdir holds."""
    selection = workdir / SELECTION_TSV
    if not selection.is_file():
        return []
    try:
        rows = read_selection(selection)
    except (OSError, ValueError, RepGenRError):
        return []
    out = []
    for row in rows:
        genome = workdir / ("outgroup" if row.is_outgroup else "genomes") / row.filename
        if genome.is_file():
            out.append(sketches_dir(workdir) / (record_name(genome) + SKETCH_SUFFIX))
    return out


def require_sourmash_if_requested(flag: bool | None) -> None:
    """Refuse an explicit --sketch without sourmash before a stage writes anything."""
    if flag:
        preflight(SOURMASH_TOOL)


def sourmash_for_sketches(flag: bool | None, stage: str, logger: logging.Logger) -> dict[str, str]:
    """Resolve --sketch/--no-sketch of a writer stage to sourmash's versions, or {}.

    ``flag`` None (the default) sketches when sourmash can run and otherwise
    logs one INFO line naming the flag; True requires sourmash
    (MissingBinaryError otherwise); False never sketches.
    """
    if flag is False:
        return {}
    if flag is None and not tool_available(SOURMASH_TOOL):
        logger.info(
            "%s: sourmash not found; genome sketches not written (--sketch requires them; "
            "'repgenr sketch' adds them later).",
            stage,
        )
        return {}
    try:
        return preflight(SOURMASH_TOOL)
    except MissingBinaryError as exc:
        if flag:
            raise
        logger.info("%s: genome sketches not written: %s", stage, exc)
        return {}


def sketch_stage_genomes(
    ctx: WorkdirContext,
    flag: bool | None,
    stage: str,
    logger: logging.Logger,
    *,
    threads: int | None = None,
    only: Collection[Path] | None = None,
    reuse: Mapping[str, SketchSource] | None = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    """The sketch step at the end of a genome-writing stage.

    Prunes sketches of genomes that left the set (always, even with
    ``--no-sketch``), then sketches the missing and stale ones when sourmash
    can run: every genome of the set, or only the genome files in ``only``
    (assemble --append). Returns the summary for the stage record and
    sourmash's versions. A sketch failure under the default fails nothing:
    sketches are not deliverables, so it is logged and the stage finishes;
    with an explicit ``--sketch`` it is raised.
    """
    try:
        removed = remove_stale(ctx, logger)
    except OSError as exc:
        if flag:
            raise
        _warn_stopped(stage, exc, ctx, logger)
        removed = 0
    versions = sourmash_for_sketches(flag, stage, logger)
    if not versions:
        state = "off" if flag is False else "unavailable"
        return {"state": state, "removed": removed}, {}
    targets = sketch_targets(ctx)
    if only is not None:
        wanted = {Path(p).name for p in only}
        targets = [t for t in targets if t.path.name in wanted]
    try:
        summary = sketch_genomes(ctx, targets, threads or default_threads(), logger, reuse=reuse)
    except (ToolExecutionError, WorkdirError, OSError) as exc:
        if flag or process.stop_requested.is_set():
            raise
        _warn_stopped(stage, exc, ctx, logger)
        return {"state": "failed", "removed": removed}, versions
    summary.removed = removed
    logger.info("%s", summary.line())
    return {"state": "done", **summary.as_dict()}, versions


def _warn_stopped(
    stage: str, exc: BaseException, ctx: WorkdirContext, logger: logging.Logger
) -> None:
    logger.warning(
        "%s: genome sketching stopped (%s); the genome set is complete without "
        "sketches, and 'repgenr sketch -wd %s' builds the missing ones.",
        stage,
        exc,
        ctx.workdir,
    )


def sources_from_workdir(workdir: Path) -> dict[str, SketchSource]:
    """Reusable sketches of an earlier workdir, by record name (for ingest --from-workdir).

    Read from its manifest (read-only); a workdir without sketches, without a
    manifest, or with a pre-v4 manifest offers none.
    """
    from .manifest import MANIFEST_FILENAME, Manifest

    directory = sketches_dir(workdir)
    path = workdir / MANIFEST_FILENAME
    if not directory.is_dir() or not path.is_file():
        return {}
    try:
        manifest = Manifest.open_readonly(path)
    except WorkdirError:
        return {}
    try:
        records = manifest.sketch_records()
    except (sqlite3.Error, WorkdirError):  # an unreadable source manifest: sketch instead
        records = {}
    finally:
        manifest.close()
    out: dict[str, SketchSource] = {}
    for rec in records.values():
        if not rec.file:
            continue
        sketch = workdir / rec.file
        name = Path(rec.file).name
        if name.endswith(SKETCH_SUFFIX):
            out[name[: -len(SKETCH_SUFFIX)]] = SketchSource(sketch, rec)
    return out


# --- consumers: dereplicate, glance, phylo and the classifier -------------------

# A source of genome sketches for a consumer: given genome FASTAs, the
# signature file of each genome it can supply (the others are left out).
SketchProvider = Callable[[Sequence[Path]], dict[Path, Path]]


def contract_mismatch(ksize: int, scaled: int) -> str | None:
    """Why sketches of the contract cannot serve (ksize, scaled); None when they can.

    The contract holds k=21, 31 and 51 at scaled=1000, so a consumer asking for
    any of those k-mer sizes at scaled=1000 selects its signature with ``-k``.
    """
    if scaled != SCALED:
        return f"scaled={scaled} differs from the scaled={SCALED} of the genome sketches"
    if ksize not in KSIZES:
        sizes = ", ".join(str(k) for k in KSIZES)
        return f"k={ksize} is not one of the sketched k-mer sizes ({sizes})"
    return None


def adapter_sketches(
    adapter: object,
    extra: Mapping[str, object],
    genomes: Sequence[Path],
    provider: SketchProvider | None,
    logger: logging.Logger,
    consumer: str,
) -> dict[Path, Path] | None:
    """The genome sketches to give ``adapter``, or None when it sketches itself.

    An adapter reads sketches when it defines ``sketch_request(extra)`` and
    that returns the (ksize, scaled) it compares at. When those parameters are
    not in the contract (:func:`contract_mismatch`), or there is no provider,
    the adapter sketches into its own work directory as before; the reason is
    logged. Genomes the provider cannot supply are sketched by the adapter.
    """
    request = getattr(adapter, "sketch_request", None)
    if request is None or provider is None:
        return None
    wanted = request(extra)
    if wanted is None:
        return None
    ksize, scaled = wanted
    reason = contract_mismatch(ksize, scaled)
    if reason is not None:
        logger.info(
            "%s: genome sketches not used (%s); the tool sketches the genomes itself.",
            consumer,
            reason,
        )
        return None
    return provider(genomes) or None


def resolve_sketches(
    ctx: WorkdirContext,
    genomes: Sequence[Path],
    logger: logging.Logger,
    threads: int,
    *,
    consumer: str,
    write: bool = True,
) -> dict[Path, Path]:
    """Genome -> its sketch under ``sketches/``, sketching what is missing or stale.

    With ``write=False`` only the sketches that are current are returned and
    nothing is written: a consumer that merely reads them (the reference
    choice) must not fill the workdir as a side effect.

    ``genomes`` are matched to the genome set of the workdir by file name, so
    the representatives under ``derep/representatives/`` find the records of
    their genomes. A sketch is reused when its record names the current
    parameters and the sha256 of the given FASTA (taken from
    ``sketches/.digests.json`` while the file's size and mtime are those it
    was hashed at, see :class:`_DigestCache`); the others are written
    through :func:`sketch_genomes` and recorded, so a consumer fills the
    sketches of the workdir as a side effect. Without sourmash, or when
    sketching fails, only the reusable sketches are returned. A genome that is
    not part of the set (no manifest row, not under ``outgroup/``) is left out;
    the consumer sketches it.
    """
    all_targets = sketch_targets(ctx)
    known = {t.path.name: t for t in all_targets}
    known_paths = [t.path for t in all_targets]
    targets: list[SketchTarget] = []
    for genome in genomes:
        hit = known.get(genome.name)
        if hit is not None:
            targets.append(SketchTarget(genome, hit.accession))
    if not targets:
        return {}
    records = _recorded(ctx)
    cache = _DigestCache(sketches_dir(ctx.workdir) / DIGESTS_JSON)
    workers = max(1, min(threads, len(targets)))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        digests = dict(zip(targets, pool.map(lambda t: cache.digest(t.path), targets), strict=True))

    def current(target: SketchTarget, recs: Mapping[str, SketchRecord]) -> bool:
        out = sketch_path(ctx, target.path)
        rel = _relative(ctx, out)
        return _is_current(recs.get(_record_key(target)), rel, digests[target], out)

    reused = [t for t in targets if current(t, records)]
    reused_set = set(reused)
    todo = [t for t in targets if t not in reused_set]
    written: list[SketchTarget] = []
    if todo and not write:
        logger.info(
            "%s: sketches: %d reused, %d missing or stale (not written here)",
            consumer,
            len(reused),
            len(todo),
        )
        return {t.path: sketch_path(ctx, t.path) for t in reused}
    if todo:
        if tool_available(SOURMASH_TOOL):
            try:
                sketch_genomes(ctx, todo, threads, logger, digests=cache)
            except (ToolExecutionError, WorkdirError, OSError) as exc:
                if process.stop_requested.is_set():
                    raise
                logger.warning(
                    "%s: writing genome sketches stopped (%s); the tool sketches the rest itself.",
                    consumer,
                    exc,
                )
            fresh = _recorded(ctx)
            written = [t for t in todo if current(t, fresh)]
        else:
            logger.info(
                "%s: sourmash not available here to write the %d missing genome sketch(es).",
                consumer,
                len(todo),
            )
    if (ctx.workdir / SKETCHES_DIR).is_dir():
        _save_cache(cache, [*known_paths, *(t.path for t in targets)], logger)
    left = len(genomes) - len(reused) - len(written)
    logger.info(
        "%s: sketches: %d reused, %d written%s",
        consumer,
        len(reused),
        len(written),
        f", {left} sketched by the tool" if left else "",
    )
    return {t.path: sketch_path(ctx, t.path) for t in [*reused, *written]}


class _DigestCache:
    """sha256 of genome FASTAs, reused while a file's size and mtime_ns are unchanged.

    The resume fingerprint judges genome directories by the same file
    metadata; an edit that keeps both the size and the modification time is
    therefore not seen (``repgenr --force sketch`` writes every sketch again).
    Keys are resolved paths, so a representative linked to its genome shares
    the entry. On save only the files of the current genome set (and those
    just hashed) are kept.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        self._lock = threading.Lock()
        self.entries: dict[str, tuple[int, int, str]] = {}
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            data = {}
        if isinstance(data, dict):
            for key, value in data.items():
                if isinstance(value, list) and len(value) == 3:
                    try:
                        self.entries[str(key)] = (int(value[0]), int(value[1]), str(value[2]))
                    except (TypeError, ValueError):
                        continue

    def digest(self, genome: Path, *, refresh: bool = False) -> str:
        """The FASTA sha256; ``refresh`` hashes the file whatever the cache holds."""
        key = str(genome.resolve())
        try:
            st = genome.stat()
        except OSError:
            return file_digest(genome)
        with self._lock:
            hit = self.entries.get(key)
        if not refresh and hit is not None and hit[:2] == (st.st_size, st.st_mtime_ns):
            return hit[2]
        value = file_digest(genome)
        with self._lock:
            self.entries[key] = (st.st_size, st.st_mtime_ns, value)
        return value

    def save(self, keep: Collection[Path]) -> None:
        """Write the entries of the files in ``keep`` (genomes of the set) that exist.

        The temporary file has a unique name, so two processes saving at
        once each replace the cache with a whole file.
        """
        wanted = {str(p.resolve()) for p in keep}
        kept = {
            k: list(v) for k, v in sorted(self.entries.items()) if k in wanted and Path(k).is_file()
        }
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(
            dir=self.path.parent, prefix=f"{self.path.name}.", suffix=".tmp"
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fo:
                json.dump(kept, fo, sort_keys=True)
                fo.write("\n")
            os.replace(name, self.path)
        except BaseException:
            Path(name).unlink(missing_ok=True)
            raise


def workdir_provider(
    ctx: WorkdirContext,
    logger: logging.Logger,
    threads: int,
    consumer: str,
    *,
    write: bool = True,
) -> SketchProvider:
    """A provider over the sketches of a workdir (:func:`resolve_sketches`);
    ``write=False`` returns the current sketches only and writes none."""

    def provide(genomes: Sequence[Path]) -> dict[Path, Path]:
        return resolve_sketches(ctx, genomes, logger, threads, consumer=consumer, write=write)

    return provide


def directory_provider(directory: Path, logger: logging.Logger, consumer: str) -> SketchProvider:
    """A provider over a ``sketches/`` directory without a manifest (Nextflow steps).

    A sketch is matched to a genome by record name only; no digest is checked,
    so the directory must have been sketched from the same genome files (the
    pipeline's SKETCH process sketches the genomes it then hands on).
    """

    def provide(genomes: Sequence[Path]) -> dict[Path, Path]:
        out: dict[Path, Path] = {}
        for genome in genomes:
            sketch = directory / (record_name(genome) + SKETCH_SUFFIX)
            if sketch.is_file():
                out[genome] = sketch
        left = len(genomes) - len(out)
        logger.info(
            "%s: sketches: %d reused from %s%s",
            consumer,
            len(out),
            directory,
            f", {left} sketched by the tool" if left else "",
        )
        return out

    return provide


def sketch_beside(genome: Path, name: str, out: Path, logger: logging.Logger) -> SketchRecord:
    """Sketch ``genome`` into ``out`` (signatures named ``name``) unless it is current.

    A stamp ``<out>.json`` records the FASTA sha256 and the parameters; while
    both match, the file is reused. The assemble stage keeps such a sketch
    beside each assembly, gathers with it, and later copies it into
    ``sketches/`` under the genome's record name. Returns the record of the
    sketch (``file`` None: it is not a workdir sketch).
    """
    stamp = out.with_name(out.name + ".json")
    digest = file_digest(genome)
    try:
        stored = json.loads(stamp.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        stored = None
    if (
        isinstance(stored, dict)
        and stored.get("params") == SKETCH_PARAMS
        and stored.get("digest") == digest
        and stored.get("name") == name
        and out.is_file()
    ):
        return SketchRecord(None, SKETCH_PARAMS, digest)
    stamp.unlink(missing_ok=True)
    sketch_file(genome, name, out, logger)
    with atomic_replace(stamp) as fo:
        json.dump({"params": SKETCH_PARAMS, "digest": digest, "name": name}, fo, sort_keys=True)
        fo.write("\n")
    return SketchRecord(None, SKETCH_PARAMS, digest)
