"""ingest stage: populate a working directory from local genome FASTAs.

The offline counterpart of ``metadata`` + ``genome`` (or ``vmetadata`` +
``vgenome``): it stages a directory of FASTA files under ``genomes/``, writes
the same ``selection.tsv`` contract the metadata stage publishes, fills the
SQLite manifest, and optionally sets one genome aside as the outgroup. The
downstream chain (``dereplicate`` -> ``phylo`` -> ``tree2tax``) then runs
without any download, which is what the local Nextflow harness already does
for the data-channel path.

Taxonomy comes from ``--selection`` (a ``selection.tsv``: accession and
filename are required; it may also carry CheckM-style quality for ``--keeper
quality`` and the ``gtdb_representative`` flag for ``--keeper gtdb``, 0 when
the column is absent) or, failing that, from
canonical ``Family_genus_species_ACCESSION.fasta`` filenames; any other name
yields an empty taxonomy and the file stem as accession.

``--from-workdir`` (repeatable, alone or beside ``--genomes-dir``) takes the
genome set of an earlier working directory: its ``selection.tsv`` rows are
read unchanged (taxonomy, quality, ``gtdb_representative``), each file is
resolved under its ``genomes/``, and the manifest source of each genome
(gtdb, sra, ncbi_virus, bvbrc, local) is kept when that workdir has a manifest.
The outgroup of a source workdir is not carried over; ``--outgroup`` sets the
outgroup of the new workdir and may name a genome from any source.

Before staging, the record of the stage that wrote each source's genome set
(the latest of genome, vgenome, ingest and assemble) is evaluated as
``status`` evaluates it: a source whose record is missing, interrupted or
stale is named in one warning, and ``--strict-sources`` refuses it instead.
"""

from __future__ import annotations

import os
import shutil
import sqlite3
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    FASTA_SUFFIXES,
    SELECTION_TSV,
    SelectionRow,
    list_fasta,
    parse_genome_filename,
    read_selection,
    strip_fasta_suffix,
    write_selection,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.integrity import FOREIGN_SOURCES, refuse_foreign_rows
from ..core.manifest import MANIFEST_FILENAME, Manifest, record_from_selection
from ..core.selection_summary import write_selection_summary
from ..core.sketches import (
    SketchSource,
    require_sourmash_if_requested,
    sketch_stage_genomes,
    sources_from_workdir,
)

OUTGROUP_ACCESSION_TXT = "outgroup_accession.txt"


@dataclass
class IngestParams:
    genomes_dir: str | None = None
    selection: str | None = None
    outgroup: str | None = None
    copy: bool = False
    # Discard genomes another entry path appended (assemble --append) instead
    # of refusing to overwrite the selection that holds them.
    drop_foreign: bool = False
    # Earlier working directories whose selection.tsv and genomes/ are merged
    # into this one.
    from_workdirs: list[str] = field(default_factory=list)
    # Genome sketches (core.sketches): None sketches when sourmash can run,
    # True requires it, False skips them. An up-to-date sketch of a
    # --from-workdir genome is copied instead of sketched again.
    sketch: bool | None = None
    # Refuse a --from-workdir whose genome-set record is not done (missing,
    # interrupted or stale) instead of warning about it.
    strict_sources: bool = False


# Stages that write the genome set (selection.tsv and genomes/) of a workdir.
GENOME_SET_STAGES = ("genome", "vgenome", "ingest", "assemble")
SOURCE_DONE = "done"
SOURCE_NO_RECORD = "no stage record"


@dataclass(frozen=True)
class SourceState:
    """The state of the record that wrote a source workdir's genome set.

    ``stage`` is None when no such record exists; ``state`` is done, stale,
    interrupted, ``no stage record``, ``unreadable record`` or ``unchecked``
    (the staleness check itself failed); ``reason`` says why when known.
    """

    path: str
    stage: str | None
    state: str
    reason: str | None = None

    @property
    def done(self) -> bool:
        return self.state == SOURCE_DONE

    def describe(self) -> str:
        if self.stage is None:
            what = self.state
        else:
            what = f"its genome set was written by {self.stage}, which is {self.state}"
        return f"--from-workdir {self.path}: {what}" + (f" ({self.reason})" if self.reason else "")

    def to_dict(self) -> dict[str, str | None]:
        return {"path": self.path, "stage": self.stage, "state": self.state}


def source_state(workdir: Path) -> SourceState:
    """Evaluate the genome-set record of an earlier working directory.

    The record is the latest finished one of GENOME_SET_STAGES, or an
    interrupted one when any is (its re-run may have left the genome set
    partial). Its state comes from the checks ``status`` applies: a record
    without a completion stamp is interrupted, and ``core.doctor.stale_stages``
    (changed inputs, missing deliverables) decides stale.
    """
    from ..core.config import CONFIG_FILENAME, Config
    from ..core.doctor import stale_stages

    path = str(workdir)
    if not (workdir / CONFIG_FILENAME).is_file():
        return SourceState(path, None, SOURCE_NO_RECORD)
    try:
        cfg = Config.load(workdir)
    except WorkdirError as exc:
        return SourceState(path, None, "unreadable record", str(exc))
    records = {n: cfg.stages[n] for n in GENOME_SET_STAGES if n in cfg.stages}
    if not records:
        return SourceState(path, None, SOURCE_NO_RECORD)
    interrupted = [n for n, r in records.items() if r.interrupted]
    if interrupted:
        return SourceState(path, interrupted[-1], "interrupted")
    stage = max(records, key=lambda n: records[n].completed or "")
    try:
        stale = stale_stages(workdir, Config(stages={stage: records[stage]}))
    except Exception as exc:  # a damaged artifact is reported, not raised
        return SourceState(path, stage, "unchecked", str(exc))
    if stage in stale:
        return SourceState(path, stage, "stale", stale[stage])
    return SourceState(path, stage, SOURCE_DONE)


def _source_states(params: IngestParams) -> list[SourceState]:
    return [source_state(Path(wd).expanduser()) for wd in params.from_workdirs]


def _refuse_sources(states: list[SourceState]) -> None:
    """--strict-sources: refuse any source whose genome-set record is not done."""
    bad = [s for s in states if not s.done]
    if bad:
        raise UserInputError(
            "--strict-sources: "
            + "; ".join(s.describe() for s in bad)
            + ". Re-run or finish the stage in that working directory, or drop "
            "--strict-sources to take its genomes as they are."
        )


# Manifest source of a genome that no source manifest describes.
LOCAL_SOURCE = "local"


@dataclass
class _Candidate:
    """A genome ingest can stage: its selection row, its file, its manifest
    source, and where it came from (named in messages)."""

    row: SelectionRow
    path: Path
    source: str
    origin: str

    @property
    def key(self) -> tuple[str, str]:
        return self.origin, self.row.filename


@dataclass
class _Plan:
    """What ingest will stage, resolved and checked before anything is written."""

    origins: list[str]
    ingroup: list[_Candidate]
    outgroup: _Candidate | None

    @property
    def chosen(self) -> list[_Candidate]:
        return [*self.ingroup, *([self.outgroup] if self.outgroup is not None else [])]

    def foreign_accessions(self) -> set[str]:
        """Accessions of genomes with an appended source (sra) that this ingest
        stages again, so that re-running it does not drop them."""
        return {c.row.accession for c in self.chosen if c.source in FOREIGN_SOURCES}


def _not_appended_here(ctx: WorkdirContext, plan: _Plan) -> set[str]:
    """Accessions the foreign-row guard need not protect in this workdir.

    Genomes this ingest stages again are not dropped. Other genomes with
    source sra are left unprotected only when they cannot have been appended
    here: the last ingest of this workdir took ``--from-workdir`` (the one
    other way an sra genome gets in) and no ``assemble`` has been recorded in
    it. Otherwise they are protected, as before ``--from-workdir`` existed.
    """
    keep = plan.foreign_accessions()
    stages = ctx.config.stages
    prior = stages.get("ingest")
    if prior is not None and prior.params.get("from_workdirs") and "assemble" not in stages:
        keep |= {
            g.accession
            for g in ctx.manifest.all_genomes(include_outgroup=True)
            if g.source in FOREIGN_SOURCES
        }
    return keep


def precheck(ctx: WorkdirContext, params: IngestParams) -> None:
    """Refuse an ingest that cannot proceed, before anything is changed.

    The CLI harness calls this before it marks the stage record incomplete, so
    a mistyped directory or a bad genome file does not leave a finished ingest
    (and the stages built on it) looking interrupted.
    """
    plan = _plan(params, ctx.workdir, logger=None)
    if params.strict_sources:
        _refuse_sources(_source_states(params))
    require_sourmash_if_requested(params.sketch)
    # Only an existing manifest can hold appended genomes; opening one in a new
    # workdir would create an empty manifest.sqlite for a refused ingest.
    if not params.drop_foreign and (ctx.workdir / MANIFEST_FILENAME).exists():
        refuse_foreign_rows(
            ctx, "ingest", drop_foreign=False, logger=ctx.logger, keep=_not_appended_here(ctx, plan)
        )


def _plan(params: IngestParams, target: Path, logger) -> _Plan:
    """Resolve and check the genome set; warn through ``logger`` when given."""
    if not params.genomes_dir and not params.from_workdirs:
        raise UserInputError("ingest needs --genomes-dir, --from-workdir, or both.")
    if params.selection and not params.genomes_dir:
        raise UserInputError(
            "--selection names genomes under --genomes-dir; give --genomes-dir with it "
            "(a --from-workdir brings its own selection.tsv)."
        )
    given = [Path(wd).expanduser().resolve() for wd in params.from_workdirs]
    repeated = sorted({str(wd) for wd in given if given.count(wd) > 1})
    if repeated:
        raise UserInputError(f"--from-workdir names {', '.join(repeated)} more than once.")
    origins: list[str] = []
    ingroup: list[_Candidate] = []
    marked: list[_Candidate] = []
    # Outgroup genomes of the source workdirs: not carried over, but
    # --outgroup may name one of them.
    others: list[_Candidate] = []
    if params.genomes_dir:
        origin = f"--genomes-dir {Path(params.genomes_dir).expanduser()}"
        origins.append(origin)
        for c in _genomes_dir_candidates(params, origin, logger):
            (marked if c.row.is_outgroup else ingroup).append(c)
    for wd in params.from_workdirs:
        origin = f"--from-workdir {Path(wd).expanduser()}"
        origins.append(origin)
        rows, outgroups = _workdir_candidates(Path(wd).expanduser(), target, origin)
        ingroup.extend(rows)
        others.extend(outgroups)

    outgroup = _resolve_outgroup(ingroup, marked, others, params.outgroup)
    if outgroup is not None:
        ingroup = [c for c in ingroup if c.key != outgroup.key]
    plan = _Plan(origins, ingroup, outgroup)
    _refuse_duplicates(plan.chosen)
    _refuse_unusable([c.path for c in plan.chosen])
    if logger is not None:
        for c in others:
            if outgroup is None or c.key != outgroup.key:
                logger.info(
                    "%s: its outgroup %s (%s) is not carried over; --outgroup sets the "
                    "outgroup of the new working directory.",
                    c.origin,
                    c.row.accession,
                    c.row.filename,
                )
    return plan


def _genomes_dir_candidates(params: IngestParams, origin: str, logger) -> list[_Candidate]:
    """The genomes of ``--genomes-dir``, from ``--selection`` or the filenames."""
    source = Path(params.genomes_dir or "").expanduser()
    if not source.is_dir():
        raise UserInputError(f"--genomes-dir {source} is not a directory.")
    files = list_fasta(source)
    others = _other_entries(source, files)
    if not files:
        subdirs = [p for p in others if p.is_dir()]
        where = (
            f"; it holds {len(subdirs)} subdirectories, which ingest does not search"
            if subdirs
            else ""
        )
        raise UserInputError(
            f"No genome FASTA files under {source} (expected {', '.join(FASTA_SUFFIXES)}){where}."
        )
    by_name = {f.name: f for f in files}

    if params.selection:
        rows = _rows_from_selection(Path(params.selection), by_name, others)
    else:
        rows = [_row_from_filename(f.name) for f in files]
        skipped = [p.name for p in others if p.is_file()]
        if skipped and logger is not None:
            logger.warning(
                "Skipped %d file(s) under %s without a FASTA suffix (%s): %s",
                len(skipped),
                source,
                ", ".join(FASTA_SUFFIXES),
                _examples(skipped),
            )
    return [_Candidate(r, by_name[r.filename], LOCAL_SOURCE, origin) for r in rows]


def _workdir_candidates(
    wd: Path, target: Path, origin: str
) -> tuple[list[_Candidate], list[_Candidate]]:
    """The ingroup and the outgroup genomes of an earlier working directory.

    Rows are read from its ``selection.tsv`` unchanged; each ingroup file must
    exist under its ``genomes/`` (an outgroup row is looked up under
    ``outgroup/`` only when ``--outgroup`` names it). The manifest source of
    each genome is kept when the workdir has a manifest.
    """
    if not wd.is_dir():
        raise UserInputError(f"{origin} is not a directory.")
    if wd.resolve() == target.resolve():
        raise UserInputError(
            f"{origin} is the working directory being written (-wd); name another one."
        )
    selection = wd / SELECTION_TSV
    genomes = wd / "genomes"
    absent = [
        name
        for name, present in ((SELECTION_TSV, selection.is_file()), ("genomes/", genomes.is_dir()))
        if not present
    ]
    if absent:
        raise UserInputError(
            f"{origin} holds no {' or '.join(absent)}; name a working directory in which "
            "a genome set was written (genome, ingest, vgenome or assemble)."
        )
    try:
        rows = read_selection(selection)
    except WorkdirError as exc:
        raise UserInputError(f"{origin}: {exc}") from exc
    sources = _manifest_sources(wd, origin)

    def candidate(row: SelectionRow, directory: Path) -> _Candidate:
        source = sources.get(row.accession) or LOCAL_SOURCE
        return _Candidate(replace(row, is_outgroup=False), directory / row.filename, source, origin)

    ingroup = [candidate(r, genomes) for r in rows if not r.is_outgroup]
    missing = [c.row.filename for c in ingroup if not c.path.exists()]
    if missing:
        raise UserInputError(
            f"{len(missing)} selection row(s) of {origin} name a file not present under "
            f"{genomes} (e.g. {', '.join(missing[:3])})."
        )
    outgroups = [candidate(r, wd / "outgroup") for r in rows if r.is_outgroup]
    return ingroup, outgroups


def _manifest_sources(wd: Path, origin: str) -> dict[str, str]:
    """accession -> manifest source of a source workdir; empty without a manifest."""
    path = wd / MANIFEST_FILENAME
    if not path.is_file():
        return {}
    try:
        manifest = Manifest.open_readonly(path)
        try:
            return {
                g.accession: g.source
                for g in manifest.all_genomes(include_outgroup=True)
                if g.source
            }
        finally:
            manifest.close()
    except (WorkdirError, sqlite3.Error) as exc:
        raise UserInputError(f"{origin}: cannot read {MANIFEST_FILENAME} ({exc}).") from exc


def run(ctx: WorkdirContext, params: IngestParams) -> int:
    logger = ctx.logger
    plan = _plan(params, ctx.workdir, logger)
    states = _source_states(params)
    if params.strict_sources:
        _refuse_sources(states)
    for state in states:
        if not state.done:
            logger.warning("%s; its genomes are taken as they are.", state.describe())
    require_sourmash_if_requested(params.sketch)
    refuse_foreign_rows(
        ctx,
        "ingest",
        drop_foreign=params.drop_foreign,
        logger=logger,
        keep=_not_appended_here(ctx, plan),
    )
    ingroup, outgroup = plan.ingroup, plan.outgroup

    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    if params.copy:
        logger.info("Copying %d genomes into %s", len(ingroup), ctx.genomes_dir)
    _prune(ctx.genomes_dir, {c.row.filename for c in ingroup}, logger)
    for c in ingroup:
        _stage(c.path, ctx.genomes_dir / c.row.filename, params.copy)

    acc_file = ctx.workdir / OUTGROUP_ACCESSION_TXT
    if outgroup is not None:
        ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
        _prune(ctx.outgroup_dir, {outgroup.row.filename}, logger)
        _stage(outgroup.path, ctx.outgroup_dir / outgroup.row.filename, params.copy)
        acc_file.write_text(outgroup.row.accession + "\n", encoding="utf-8")
        logger.info("Outgroup: %s (%s)", outgroup.row.accession, outgroup.row.filename)
    else:
        if ctx.outgroup_dir.exists():
            _prune(ctx.outgroup_dir, set(), logger)
        acc_file.unlink(missing_ok=True)

    ctx.manifest.replace_genomes([record_from_selection(c.row, c.source) for c in plan.chosen])
    write_selection(ctx.workdir / SELECTION_TSV, [c.row for c in plan.chosen])
    write_selection_summary(ctx.workdir, [c.row for c in plan.chosen], logger)
    reuse: dict[str, SketchSource] = {}
    for wd in params.from_workdirs:
        reuse.update(sources_from_workdir(Path(wd).expanduser()))
    sketches, sketch_versions = sketch_stage_genomes(
        ctx, params.sketch, "ingest", logger, reuse=reuse
    )

    mode = "copied" if params.copy else "linked"
    per_origin = [(o, sum(c.origin == o for c in ingroup)) for o in plan.origins]
    logger.info(
        "Ingested %d genomes (%s): %s",
        len(ingroup),
        mode,
        "; ".join(f"{n} from {o}" for o, n in per_origin),
    )
    by_source: dict[str, int] = {}
    for c in ingroup:
        by_source[c.source] = by_source.get(c.source, 0) + 1
    ctx.config.record_stage(
        "ingest",
        params={
            "genomes_dir": (
                str(Path(params.genomes_dir).expanduser()) if params.genomes_dir else None
            ),
            # The source workdirs, as given (doctor re-derives the resume
            # inputs from them).
            "from_workdirs": [str(Path(wd).expanduser()) for wd in params.from_workdirs],
            "selection": params.selection,
            # The raw flag (doctor re-derives the resume inputs from it) and
            # the accession it resolved to.
            "outgroup": params.outgroup,
            "outgroup_accession": outgroup.row.accession if outgroup is not None else None,
            "copy": params.copy,
            "drop_foreign": params.drop_foreign,
            "total": len(ingroup),
            # Ingroup genomes per manifest source (gtdb, sra, local, ...).
            "sources": dict(sorted(by_source.items())),
            # Per --from-workdir: the stage that wrote its genome set and the
            # state of that record when this ingest read it.
            "source_states": [state.to_dict() for state in states],
            "sketches": sketches,
        },
        tool_versions=sketch_versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(ingroup)


def _rows_from_selection(
    path: Path, by_name: dict[str, Path], others: list[Path]
) -> list[SelectionRow]:
    if not path.is_file():
        raise UserInputError(f"--selection {path} is not a file.")
    try:
        rows = read_selection(path)
    except WorkdirError as exc:
        # The file is the user's input here, not workdir state.
        raise UserInputError(f"--selection: {exc}") from exc
    missing = [r.filename for r in rows if r.filename not in by_name]
    if missing:
        unsupported = {p.name for p in others} & set(missing)
        note = (
            f"; {len(unsupported)} of them exist but lack a FASTA suffix "
            f"({', '.join(FASTA_SUFFIXES)})"
            if unsupported
            else ""
        )
        raise UserInputError(
            f"{len(missing)} selection row(s) name a file not present under --genomes-dir "
            f"(e.g. {', '.join(missing[:3])}){note}."
        )
    return rows


def _other_entries(source: Path, files: list[Path]) -> list[Path]:
    """Entries under ``source`` that are not genome FASTA files (dotfiles excluded)."""
    taken = {f.name for f in files}
    return sorted(p for p in source.iterdir() if not p.name.startswith(".") and p.name not in taken)


def _examples(names: list[str], limit: int = 5) -> str:
    more = f" (+{len(names) - limit} more)" if len(names) > limit else ""
    return ", ".join(names[:limit]) + more


def _refuse_unusable(paths: list[Path]) -> None:
    """Fail on empty or unreadable genome files (a dangling link included) before
    anything is staged.

    One ``stat`` per file: opening each file to check its first bytes cost about
    50 s for 1000 genomes on an exFAT disk, and this check also runs on the
    resume path. ``doctor`` checks the content of the staged genomes.
    """
    bad: list[str] = []
    for path in paths:
        try:
            if path.stat().st_size == 0:
                bad.append(path.name)
        except OSError:
            bad.append(path.name)
    if bad:
        raise UserInputError(
            f"{len(bad)} genome file(s) are empty or unreadable: {_examples(bad)}. "
            "Remove them from --genomes-dir (or from --selection, or from the "
            "selection.tsv of a --from-workdir) and re-run."
        )


def _refuse_duplicates(chosen: list[_Candidate]) -> None:
    """Fail when two genomes share an accession or a filename.

    The manifest is keyed by accession, so a shared accession kept one row of
    several while genomes/ and selection.tsv kept them all. Non-canonical names
    with four or more tokens (``sample_1_run_A.fasta``) are parsed as
    Family_genus_species_ACCESSION and can collide this way. Two sources
    (``--genomes-dir`` and the ``--from-workdir`` directories) that hold the
    same genome are refused as well: neither is given precedence.
    """
    for attr, label in (("accession", "an accession"), ("filename", "a filename")):
        seen: dict[str, _Candidate] = {}
        within: list[str] = []
        across: list[str] = []
        for c in chosen:
            key = getattr(c.row, attr)
            prior = seen.get(key)
            if prior is None:
                seen[key] = c
            elif prior.origin == c.origin:
                within.append(f"{key} ({prior.row.filename}, {c.row.filename})")
            else:
                across.append(
                    f"{key} ({prior.row.filename} from {prior.origin}; "
                    f"{c.row.filename} from {c.origin})"
                )
        if across:
            raise UserInputError(
                f"{len(across)} genome(s) from different sources share {label}: "
                f"{_shown(across)}. Remove each such genome from all but one source "
                "(its selection.tsv or genome directory) and re-run."
            )
        if within:
            hint = (
                " Give each genome its own accession with --selection, or rename the files "
                "(Family_genus_species_ACCESSION.fasta)."
                if attr == "accession"
                else " List each file once in --selection."
            )
            raise UserInputError(f"{len(within)} genome(s) share {label}: {_shown(within)}.{hint}")


def _shown(items: list[str], limit: int = 3) -> str:
    return "; ".join(items[:limit]) + (
        f" (+{len(items) - limit} more)" if len(items) > limit else ""
    )


def _row_from_filename(name: str) -> SelectionRow:
    family, genus, species, accession = parse_genome_filename(name)
    return SelectionRow(accession, family, genus, species, False, name)


def _resolve_outgroup(
    ingroup: list[_Candidate],
    marked: list[_Candidate],
    others: list[_Candidate],
    flag: str | None,
) -> _Candidate | None:
    """The outgroup genome, from the ``--selection`` outgroup row or ``--outgroup``.

    ``marked`` holds the rows ``--selection`` marks as outgroup; ``others`` the
    outgroup genomes of the ``--from-workdir`` directories, which only
    ``--outgroup`` selects. ``--outgroup`` names a genome of any source
    (filename, stem or accession) or points at a FASTA file anywhere. When the
    selection also marks an outgroup, both must name the same genome; a
    conflict is an error rather than a silent drop of the selection's
    outgroup row.
    """
    if len(marked) > 1:
        raise UserInputError("The selection marks more than one genome as outgroup.")
    marked_rows = [c.row for c in marked]
    if flag is None:
        return _as_outgroup(marked[0]) if marked else None

    known = [*ingroup, *marked, *others]
    path = Path(flag).expanduser()
    if path.is_file():
        # phylo, tree2tax and doctor resolve the outgroup among the FASTA
        # files under outgroup/ only; a file staged under another name would
        # leave the tree unrooted.
        if not path.name.endswith(FASTA_SUFFIXES):
            raise UserInputError(
                f"--outgroup {flag} has no FASTA suffix; name it with one of "
                f"{', '.join(FASTA_SUFFIXES)}."
            )
        # A path to a genome of a source is that genome's row, so its
        # selection accession is kept rather than one parsed from the filename.
        resolved = path.resolve()
        for c in known:
            if c.path.resolve() == resolved:
                _refuse_conflict(marked_rows, c.row, flag)
                return _as_outgroup(c)
        row = replace(_row_from_filename(path.name), is_outgroup=True)
        _refuse_conflict(marked_rows, row, flag)
        _refuse_external_clash([c.row for c in ingroup], row, flag)
        return _Candidate(row, path, LOCAL_SOURCE, f"--outgroup {path}")
    for c in known:
        if flag in (c.row.filename, strip_fasta_suffix(c.row.filename), c.row.accession):
            _refuse_conflict(marked_rows, c.row, flag)
            return _as_outgroup(c)
    raise UserInputError(
        f"--outgroup {flag!r} is neither a FASTA file nor a genome of --genomes-dir or "
        "--from-workdir."
    )


def _as_outgroup(c: _Candidate) -> _Candidate:
    return replace(c, row=replace(c.row, is_outgroup=True))


def _refuse_conflict(from_selection: list[SelectionRow], chosen: SelectionRow, flag: str) -> None:
    """Fail when ``--outgroup`` names another genome than the selection's outgroup row."""
    if not from_selection:
        return
    marked = from_selection[0]
    if (marked.filename, marked.accession) == (chosen.filename, chosen.accession):
        return
    raise UserInputError(
        f"The selection marks {marked.filename} ({marked.accession}) as the outgroup, "
        f"but --outgroup {flag!r} names {chosen.filename} ({chosen.accession}). "
        "Name the same genome in both, or drop one of them."
    )


def _refuse_external_clash(rows: list[SelectionRow], chosen: SelectionRow, flag: str) -> None:
    """Fail when an outgroup file from outside the genome sources shares an ingroup
    genome's filename or accession.

    Both name a genome in genomes/, the manifest and the tree, so the outgroup
    would silently replace that ingroup genome rather than sit beside it.
    """
    for row in rows:
        if row.filename == chosen.filename or row.accession == chosen.accession:
            raise UserInputError(
                f"--outgroup {flag!r} is a file outside --genomes-dir and --from-workdir, "
                f"but its name gives {chosen.filename} ({chosen.accession}), which is also "
                f"the ingroup genome {row.filename} ({row.accession}). Rename the outgroup "
                "file, or name the ingroup genome to set that genome aside."
            )


def _stage(src: Path, dst: Path, copy: bool) -> None:
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    if copy:
        shutil.copy2(src, dst)
    else:
        # abspath, not resolve(): a resolved macOS firmlink path is not what the
        # container backend bind-mounts (see write_fofn in core.process).
        os.symlink(os.path.abspath(src), dst)


def _prune(directory: Path, keep: set[str], logger) -> None:
    """Remove genome FASTA files (or links) no longer part of the selection."""
    for f in directory.iterdir():
        if f.name in keep or not f.name.endswith(FASTA_SUFFIXES):
            continue
        if f.is_symlink() or f.is_file():
            f.unlink()
            logger.info("Removed %s (no longer selected)", f.name)
