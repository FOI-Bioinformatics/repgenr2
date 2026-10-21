"""Auxiliary commands: status, doctor, markers-db, census, glance, derep-unpack, derep-stock,
sketch, list-tools."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import typer

from ..core.errors import UserInputError
from .base import (
    _RUN_STATE,
    DEFAULT_THREADS,
    HELP_THREADS,
    HELP_WORKDIR,
    PANEL_ENV,
    PANEL_INSPECT,
    PANEL_PIPELINE,
    PIPELINE_BACTERIAL,
    PIPELINE_LOCAL,
    PIPELINE_READS,
    PIPELINE_VIRAL,
    _echo_json,
    _require_choice,
    _run,
    app,
    require_existing_workdir,
    resolve_threads,
)


@app.command(rich_help_panel=PANEL_ENV)
def versions(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    versions_out: Path | None = typer.Option(
        None, "--versions-out", help="Write a versions.yml fragment here instead of stdout."
    ),
) -> None:
    """Print the external-tool versions recorded in a workdir's repgenr.yaml.

    One 'tool: version' line per tool. A tool that stages recorded with
    different versions (an image in one, the host binary in another) is listed
    once per stage as 'tool (stage): version'. A containerized tool's version
    is its image reference. A stage that did not finish is named on stderr.

    Lets the Nextflow bridge modules (which run a full stage in a scratch
    workdir) surface the resolved tool versions into versions.yml.
    """
    from ..core.config import CONFIG_FILENAME, Config
    from ..core.errors import WorkdirError
    from ..core.logging import configure_logging
    from ..core.versions import merge_stage_versions, write_versions_fragment

    # A warning raised while loading records or plugins carries the standard
    # timestamp and level instead of Python's bare last-resort line.
    configure_logging(None, level=min(_RUN_STATE["log_level"], logging.WARNING))
    if not (workdir / CONFIG_FILENAME).exists():
        # A wrong -wd would otherwise print nothing and exit 0.
        err = WorkdirError(f"No RepGenR2 run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo(f"ERROR {err}", err=True)
        raise typer.Exit(code=err.exit_code)
    from ..stages.metadata import gtdb_provenance

    try:
        cfg = Config.load(workdir)
    except WorkdirError as exc:
        typer.echo(f"ERROR {exc}", err=True)
        raise typer.Exit(code=exc.exit_code) from exc
    for name, record in cfg.stages.items():
        if record.tool_versions and not record.completed:
            typer.echo(
                f"WARNING {name} did not finish; its versions are from its last recorded run.",
                err=True,
            )
    merged = merge_stage_versions({n: r.tool_versions for n, r in cfg.stages.items()})
    metadata_record = cfg.stages.get("metadata")
    if metadata_record is not None:
        merged.update(gtdb_provenance(metadata_record.params))
    if versions_out is not None:
        write_versions_fragment(versions_out, merged)
    else:
        for tool, ver in sorted(merged.items()):
            typer.echo(f"{tool}: {ver}")


def _gtdb_note(params: dict) -> str:
    """The GTDB release, or the date of the API query, for the metadata line."""
    from ..stages.metadata import gtdb_provenance

    gtdb = gtdb_provenance(params)
    if "gtdb_release" in gtdb:
        return f"GTDB release {gtdb['gtdb_release']}"
    if "gtdb_api_query_date" in gtdb:
        return f"GTDB API queried {gtdb['gtdb_api_query_date']}"
    return ""


def _ingest_note(params: dict) -> str:
    """The source workdirs of an ingest that merged earlier working directories.

    A source whose genome-set record was not done when ingest read it is
    named with that stage and state (the ``source_states`` of the record).
    """
    workdirs = params.get("from_workdirs") or []
    if not workdirs:
        return ""
    note = f"from workdirs: {', '.join(workdirs)}"
    flagged = [
        f"{s.get('path')} ({s['stage'] + ' ' if s.get('stage') else ''}{s.get('state')})"
        for s in params.get("source_states") or []
        if isinstance(s, dict) and s.get("state") != "done"
    ]
    if flagged:
        note += f"; sources not done when ingested: {', '.join(flagged)}"
    return note


def _assemble_note(params: dict) -> str:
    """'n assembled, n deferred' for an assemble call limited by --max-runs."""
    if params.get("max_runs") is None and not params.get("n_deferred"):
        return ""
    return f"{params.get('n_assembled', 0)} assembled, {params.get('n_deferred', 0)} deferred"


STATUS_SCHEMA = "repgenr.status/1"
DOCTOR_SCHEMA = "repgenr.doctor/1"
ENTRY_HINT = (
    "Start with 'repgenr metadata' (bacteria), 'vmetadata' (viruses), "
    "'ingest' (local genomes) or 'reads' (sequencing reads), with -wd <wd>."
)
# The stages every lineage shares, followed when no entry stage is recorded.
_SHARED_TAIL = ("dereplicate", "phylo", "tree2tax")
# Shown after a completed record that holds no resume fingerprint.
NO_FINGERPRINT_NOTE = (
    "no resume fingerprint: recorded by an older version or restored by derep-stock "
    "unpack; its next invocation recomputes it"
)
HELP_JSON = (
    "Print one versioned JSON object on stdout instead of the text report "
    "(schema in docs/output.md); exit codes are unchanged."
)


# Stages whose status line reports the sketches of the genome set.
_SKETCH_STAGES = frozenset({"genome", "vgenome", "ingest", "assemble", "sketch"})


def _sketch_note(workdir: Path) -> str:
    """'sketches: n/m' when sketches/ exists: genomes with a sketch / genomes present."""
    from ..core.sketches import sketch_counts

    counts = sketch_counts(workdir)
    if counts is None:
        return ""
    return f"sketches: {counts[0]}/{counts[1]}"


def _reads_sketch_note(workdir: Path) -> str:
    """'variants: all (3 of 5 runs), spots_100000 (5 of 5 runs)' from reads_sketches/."""
    from ..core.contracts import READS_SKETCHES_DIR
    from ..stages.reads_sketch import describe_variants

    found = describe_variants(workdir / READS_SKETCHES_DIR)
    return "" if found == "none" else f"variants: {found}"


def _reads_screen_note(workdir: Path) -> str:
    """'variants: spots_100000 (3 passed, 1 excused of 5 runs)' from reads_screen/."""
    from ..stages.reads_screen import READS_SCREEN_DIR, describe_screens

    found = describe_screens(workdir / READS_SCREEN_DIR)
    return "" if found == "none" else f"variants: {found}"


@app.command(rich_help_panel=PANEL_PIPELINE)
def status(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    as_json: bool = typer.Option(False, "--json", help=HELP_JSON),
) -> None:
    """Show which pipeline stages have completed in a working directory.

    A completed stage is listed as stale when one of its inputs changed or
    one of its outputs is missing since it finished (it re-runs on its next
    invocation), and as interrupted when it did not finish. A -wd that does
    not exist exits 3; an existing directory without repgenr.yaml prints
    which entry stage to run first and exits 0. With --json, a malformed
    record still exits 3 and leaves stdout empty.
    """
    from ..core.config import CONFIG_FILENAME, Config
    from ..core.errors import WorkdirError

    require_existing_workdir(workdir)
    if not (workdir / CONFIG_FILENAME).exists():
        notes = [ENTRY_HINT]
        if as_json:
            _echo_json(_status_envelope(workdir, None, [], None, notes, None))
            raise typer.Exit()
        typer.echo(f"No RepGenR2 run found at {workdir} (no {CONFIG_FILENAME}).")
        typer.echo(notes[0])
        raise typer.Exit()

    try:
        cfg = Config.load(workdir)
    except WorkdirError as exc:
        typer.echo(f"ERROR {exc}", err=True)
        raise typer.Exit(code=exc.exit_code) from exc
    report = _status_report(workdir, cfg)
    if as_json:
        _echo_json(report)
    else:
        _render_status_text(report)


def _status_envelope(
    workdir: Path,
    pipeline: str | None,
    stages: list[dict[str, Any]],
    next_stage: str | None,
    notes: list[str],
    unchecked: str | None,
) -> dict[str, Any]:
    from .. import __version__

    return {
        "schema": STATUS_SCHEMA,
        "repgenr": __version__,
        "workdir": str(workdir),
        "pipeline": pipeline,
        "stages": stages,
        "next": next_stage,
        "notes": notes,
        "unchecked": unchecked,
    }


def _status_report(workdir: Path, cfg: Any) -> dict[str, Any]:
    """The state of each stage of a loaded workdir record, as plain data.

    Shared by the text and JSON renderers of `status`. Each stage entry has
    name, in_chain, state (done, stale, interrupted or pending), reason (the
    stale reason), tool, completed, fingerprint and detail (the GTDB note of
    metadata, the source workdirs of ingest).
    """
    from ..core.config import CONFIG_FILENAME
    from ..core.doctor import holds_outputs, stale_stages

    recorded = cfg.stages
    if not recorded:
        # An emptied or replaced record: no lineage to follow, so the entry
        # stages are named instead of guessing one.
        notes: list[str] = [f"{CONFIG_FILENAME} records no stage.", ENTRY_HINT]
        if holds_outputs(workdir):
            notes.append(
                f"The workdir holds outputs that {CONFIG_FILENAME} does not record; "
                "run repgenr doctor."
            )
        return _status_envelope(workdir, None, [], None, notes, None)
    chain: tuple[str, ...]
    lineage: str | None
    if any(name in recorded for name in ("reads", "assemble")) and "metadata" not in recorded:
        lineage, chain = "reads", PIPELINE_READS
    elif "ingest" in recorded:
        lineage, chain = "local", PIPELINE_LOCAL
    elif any(name in recorded for name in ("vmetadata", "vgenome")):
        lineage, chain = "viral", PIPELINE_VIRAL
    elif any(name in recorded for name in ("metadata", "genome")):
        lineage, chain = "bacterial", PIPELINE_BACTERIAL
    else:
        # Stages are recorded, but none of an entry stage (for example a
        # dereplicate record restored into a fresh workdir): follow the
        # stages every lineage shares.
        lineage, chain = None, _SHARED_TAIL

    unchecked: str | None
    try:
        stale = stale_stages(workdir, cfg)
        unchecked = None
    except Exception as exc:  # a damaged artifact must not stop the report
        stale, unchecked = {}, str(exc)

    def entry(stage: str, in_chain: bool) -> dict[str, Any]:
        rec = recorded.get(stage)
        if rec is None:
            return {
                "name": stage,
                "in_chain": in_chain,
                "state": "pending",
                "reason": None,
                "tool": None,
                "completed": None,
                "fingerprint": False,
                "detail": None,
            }
        if rec.interrupted:
            state = "interrupted"
        elif stage in stale:
            state = "stale"
        else:
            state = "done"
        detail = None
        if stage == "metadata" and not rec.interrupted:
            detail = _gtdb_note(rec.params) or None
        elif stage == "ingest" and not rec.interrupted:
            detail = _ingest_note(rec.params) or None
        elif stage == "assemble" and not rec.interrupted:
            detail = _assemble_note(rec.params) or None
        elif stage == "reads_sketch" and not rec.interrupted:
            detail = _reads_sketch_note(workdir) or None
        elif stage == "reads_screen" and not rec.interrupted:
            detail = _reads_screen_note(workdir) or None
        if stage in _SKETCH_STAGES and not rec.interrupted:
            sketches = _sketch_note(workdir)
            if sketches:
                detail = f"{detail}; {sketches}" if detail else sketches
        return {
            "name": stage,
            "in_chain": in_chain,
            "state": state,
            "reason": stale.get(stage) if state == "stale" else None,
            "tool": rec.tool or None,
            "completed": rec.completed,
            "fingerprint": bool(rec.fingerprint),
            "detail": detail,
        }

    stages = [entry(stage, True) for stage in chain]
    stages += [entry(stage, False) for stage in recorded if stage not in chain]
    next_stage = next((s["name"] for s in stages if s["in_chain"] and s["state"] != "done"), None)
    notes = []
    if next_stage is not None:
        hint = _next_stage_note(workdir, next_stage, recorded)
        if hint:
            notes.append(hint)
    return _status_envelope(workdir, lineage, stages, next_stage, notes, unchecked)


def _render_status_text(report: dict[str, Any]) -> None:
    """The text form of `status`: one line per stage, then the next step."""
    workdir = report["workdir"]
    typer.echo(f"RepGenR2 workdir: {workdir}")
    if not report["stages"]:
        for line in report["notes"]:
            typer.echo(line)
        return
    typer.echo(f"Pipeline: {report['pipeline'] or 'unrecorded entry stage'}\n")

    for s in report["stages"]:
        if not s["in_chain"]:
            continue
        stage = s["name"]
        tool = f" [{s['tool']}]" if s["tool"] else ""
        note = f"  ({s['detail']})" if s["detail"] else ""
        if s["state"] == "done":
            typer.echo(f"  [done]    {stage}{tool}  {s['completed']}{note}{_fingerprint_note(s)}")
        elif s["state"] == "stale":
            # Completed, but an input changed or an output is missing since:
            # the stage re-runs on its next invocation.
            typer.echo(f"  [stale]   {stage}{tool}  {s['completed']}{note}  ({s['reason']})")
        elif s["state"] == "interrupted":
            # The stage started a (re-)run and failed or was killed; outputs
            # may be partial.
            typer.echo(
                f"  [interrupted] {stage}  "
                "(did not finish; outputs may be partial; see repgenr.log)"
            )
        else:
            marker = "next" if stage == report["next"] else "    "
            typer.echo(f"  [{marker}] {stage}")

    extras = [s for s in report["stages"] if not s["in_chain"]]
    if extras:
        typer.echo("\n  optional stages run:")
        for s in extras:
            tool = f" [{s['tool']}]" if s["tool"] else ""
            if s["state"] == "interrupted":
                when = "[interrupted] (did not finish; outputs may be partial; see repgenr.log)"
            elif s["state"] == "stale":
                when = f"{s['completed']}  [stale] ({s['reason']})"
            else:
                when = (s["completed"] or "") + _fingerprint_note(s)
            note = f"  ({s['detail']})" if s["detail"] and s["state"] != "interrupted" else ""
            typer.echo(f"    {s['name']}{tool}  {when}{note}")

    if report["unchecked"] is not None:
        typer.echo(f"\nStale stages were not checked ({report['unchecked']}); run repgenr doctor.")
    if report["next"] is None:
        typer.echo("\nAll stages complete. Deliverables: tree2tax.tsv, genomes_map.tsv.")
    else:
        typer.echo(f"\nNext: repgenr {report['next']} -wd {workdir} ...")
        for line in report["notes"]:
            typer.echo(line)


def _fingerprint_note(stage: dict[str, Any]) -> str:
    """The note after a done record that its next invocation recomputes.

    Only done lines carry it: a stale line already says the stage re-runs.
    """
    if stage["state"] != "done" or stage["fingerprint"]:
        return ""
    return f"  ({NO_FINGERPRINT_NOTE})"


def _next_stage_note(workdir: Path, stage: str, recorded: dict[str, Any]) -> str | None:
    """A known refusal of the suggested stage, said before the user runs it."""
    if stage != "phylo" or "dereplicate" not in recorded:
        return None
    from ..core.contracts import list_fasta
    from ..stages.phylo import MIN_TREE_GENOMES

    phylo = recorded.get("phylo")
    if phylo is not None and phylo.params.get("all_genomes"):
        return None
    count = len(list_fasta(workdir / "derep" / "representatives"))
    if count >= MIN_TREE_GENOMES:
        return None
    return (
        f"Note: derep/representatives holds {count} genome(s) and a tree needs at least "
        f"{MIN_TREE_GENOMES}; run phylo with --all-genomes, or dereplicate again with a "
        "higher --secondary-ani."
    )


@app.command(rich_help_panel=PANEL_ENV)
def doctor(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    quick: bool = typer.Option(
        False,
        "--quick",
        help="Skip reading the first bytes of each genome file (the FASTA content check), "
        "the slowest check on large genome sets. Links, missing and untracked genomes are "
        "still checked.",
    ),
    as_json: bool = typer.Option(False, "--json", help=HELP_JSON),
) -> None:
    """Verify a workdir's outputs against its records (read-only health check).

    `status` lists each stage as done, stale or interrupted; `doctor` also
    checks the outputs themselves: missing, corrupt or untracked genomes,
    dangling links, manifest drift, representative/cluster mismatches,
    truncated tree and tree2tax tables, missing deliverables, stages whose
    inputs changed since completion, a recorded marker database that moved
    or changed version, and leftover temp files.
    Exits 0 when only warnings are found (a stale stage re-runs on its next
    invocation), 7 when any failure is found (including a malformed
    repgenr.yaml or a check that could not complete), 3 when the workdir
    does not exist, and 1 only on an unexpected error.
    """
    from ..core.doctor import diagnose
    from ..core.errors import DOCTOR_FAILURES_EXIT

    require_existing_workdir(workdir)

    findings = diagnose(workdir, quick=quick)
    order = {"fail": 0, "warn": 1, "ok": 2}
    findings = sorted(findings, key=lambda f: (order[f.level], f.area))
    counts = {level: sum(1 for f in findings if f.level == level) for level in order}
    exit_code = DOCTOR_FAILURES_EXIT if counts["fail"] else 0
    if as_json:
        from dataclasses import asdict

        from .. import __version__

        _echo_json(
            {
                "schema": DOCTOR_SCHEMA,
                "repgenr": __version__,
                "workdir": str(workdir),
                "quick": quick,
                "findings": [asdict(f) for f in findings],
                "counts": counts,
                "exit_code": exit_code,
            }
        )
    else:
        label = {"ok": "OK  ", "warn": "WARN", "fail": "FAIL"}
        for f in findings:
            typer.echo(f"[{label[f.level]}] {f.area}: {f.message}")
        typer.echo(f"\n{counts['fail']} failure(s), {counts['warn']} warning(s).")
    if exit_code:
        raise typer.Exit(code=exit_code)


_MARKERS_DB_ACTIONS = ("download", "info")
_MARKERS_DB_SETS = ("bac120", "ar53", "all")


def _markers_db_set(value: str) -> str:
    if value not in _MARKERS_DB_SETS:
        raise typer.BadParameter(f"choose from {', '.join(_MARKERS_DB_SETS)}; got '{value}'.")
    return value


@app.command(name="markers-db", rich_help_panel=PANEL_ENV)
def markers_db(
    action: str = typer.Argument(
        ...,
        metavar="ACTION",
        help="download: fetch the marker models into --output; info: summarise the "
        "database in DIRECTORY.",
    ),
    directory: Path | None = typer.Argument(
        None,
        metavar="[DIRECTORY]",
        help="With info: the database directory (default: $REPGENR_MARKERS_DB).",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        help="With download: the database directory, new, empty or an earlier download "
        "(families already present are kept).",
    ),
    marker_set: str = typer.Option(
        "all",
        "--set",
        callback=_markers_db_set,
        help="With download: bac120, ar53 or all (default) GTDB marker families.",
    ),
    threads: int = typer.Option(
        8,
        "-t",
        "--threads",
        min=1,
        callback=resolve_threads,
        help="With download: concurrent requests (default 8, or the CPU limit when lower).",
    ),
) -> None:
    """Download or describe the GTDB marker-gene database (HMMER3 models).

    download fetches GTDB-Tk's single-copy marker families (bac120: 120,
    ar53: 53, 168 distinct families, about 25 MB) one at a time from the InterPro
    API and writes one HMM file per set, an index (markers.tsv) and
    version.json. The model versions are InterPro's current ones and can be
    newer than those GTDB-Tk pins. The database is written beside --output
    and renamed into place, so an interrupted download leaves none. It is used by
    the markers core-gene builder; set REPGENR_MARKERS_DB to its path.
    info prints the sets, model counts and versions of a database.
    Exits 2 on a wrong argument or a missing database, 3 when a model cannot
    be fetched.
    """
    from ..core.logging import configure_logging
    from ..coregenes import markers_db as mdb
    from .base import stage_errors

    if action not in _MARKERS_DB_ACTIONS:
        raise typer.BadParameter(
            f"choose from {', '.join(_MARKERS_DB_ACTIONS)}; got '{action}'.", param_hint="ACTION"
        )
    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    if action == "download":
        if output is None:
            raise typer.BadParameter(
                "download needs the database directory.", param_hint="--output"
            )
        if directory is not None:
            raise typer.BadParameter(
                "download takes the directory as --output, not as an argument.",
                param_hint="DIRECTORY",
            )
        sets = mdb.MARKER_SETS if marker_set == "all" else (marker_set,)
        with stage_errors(logger):
            mdb.download(output, sets, logger, threads=threads)
        return
    with stage_errors(logger):
        db = mdb.check_markers_db(directory or output)
        info = mdb.summary(db)
    typer.echo(f"Marker database: {info['path']}")
    typer.echo(f"  version:       {info['version']}")
    typer.echo(f"  sets:          {', '.join(info['sets'])} ({info['families']} families)")
    typer.echo(f"  GTDB-Tk lists: {info['gtdbtk_lists']}")
    typer.echo(f"  fetched:       {info['fetched']} from {info['interpro_api']}")
    for name, entry in info["per_set"].items():
        typer.echo(
            f"  {name}: {entry['models']} models ({entry['pfam']} Pfam with GA, "
            f"{entry['tigrfam']} TIGRFAM with NC), {entry['hmm_bytes']} bytes"
        )


def _census_source(value: str | None) -> str | None:
    """Reject an unknown --source before any request is sent."""
    from ..stages.census import SOURCE_CHOICES

    if value is not None and value not in SOURCE_CHOICES:
        raise typer.BadParameter(f"choose from {', '.join(SOURCE_CHOICES)}; got '{value}'.")
    return value


def _census_released_after(value: str | None) -> str | None:
    from .cmd_viral import _validate_released_after

    return _validate_released_after(value)


@app.command(name="census", rich_help_panel=PANEL_INSPECT)
def census(
    workdir: Path | None = typer.Option(
        None,
        "-wd",
        "--workdir",
        help="Count the candidates (after metadata or vmetadata) or the selection (after "
        "genome, vgenome, ingest or assemble) of this working directory; nothing is "
        "written to it.",
    ),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help="Count this family: one row per genus."
    ),
    target_genus: str | None = typer.Option(
        None, "-tg", "--target-genus", help="Count this genus: one row per species."
    ),
    viral: bool = typer.Option(
        False, "--viral", help="Count a viral taxon (NCBI Virus) instead of a GTDB taxon."
    ),
    target: str | None = typer.Option(
        None,
        "--target",
        help="With --viral: the virus taxon, as vmetadata takes it (e.g. picornaviridae).",
    ),
    source: str | None = typer.Option(
        None,
        "--source",
        callback=_census_source,
        help="api (GTDB API, default) or table (GTDB metadata table, also 'tsv'); with "
        "--viral, ncbi_virus (default). bvbrc is counted only from a vmetadata workdir.",
    ),
    release: str | None = typer.Option(
        None, "-r", "--release", help="GTDB release of the table source, e.g. 232.0."
    ),
    gtdb_version: str | None = typer.Option(
        None, "--gtdb-version", help="GTDB table of the table source: bac120 (default) or ar53."
    ),
    metadata_path: str | None = typer.Option(
        None,
        "--metadata-path",
        help="Read this GTDB metadata table (.tsv.gz) instead of downloading the release "
        "table; with -wd, the table to count the candidates from.",
    ),
    runs: bool = typer.Option(
        False,
        "--runs",
        help="Also count the ENA whole-genome sequencing runs under the taxon: runs, "
        "biosamples and runs per platform (bacterial taxa only).",
    ),
    host: str | None = typer.Option(
        None, "--host", help="With --viral: only records from this host species."
    ),
    complete_only: bool = typer.Option(
        False, "--complete-only", help="With --viral: only sequences marked complete."
    ),
    released_after: str | None = typer.Option(
        None,
        "--released-after",
        callback=_census_released_after,
        help="With --viral: only records released after this date (MM/DD/YYYY).",
    ),
    tsv: Path | None = typer.Option(
        None, "--tsv", help="Also write the rows to this TSV file (format in docs/output.md)."
    ),
    as_json: bool = typer.Option(
        False,
        "--json",
        help="Print one JSON object (taxon, mode, totals, rows) on stdout instead of the table.",
    ),
) -> None:
    """Count the genera, species and samples under a taxon (read-only).

    Without -wd, queries a GTDB family (-tf) or genus (-tg) through the GTDB
    API or a GTDB metadata table, optionally with the ENA sequencing runs
    (--runs), or a viral taxon through NCBI Virus (--viral --target; metadata
    only, no sequences). With -wd, counts the candidates the entry stage found
    or the genomes it selected, by source; -tf and -tg narrow the count. A
    family is counted per genus and a genus per species. No stage is recorded
    and nothing is written to a working directory. Exits 2 for a missing taxon
    or an unsupported combination, 3 when a request fails or -wd does not exist.
    """
    from ..core.logging import configure_logging
    from ..stages.census import CensusParams, render, run_census, write_tsv
    from .base import stage_errors

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    if workdir is not None:
        require_existing_workdir(workdir)
    params = CensusParams(
        workdir=workdir,
        target_family=target_family,
        target_genus=target_genus,
        viral=viral,
        target=target,
        source=source,
        release=release,
        gtdb_version=gtdb_version,
        metadata_path=metadata_path,
        runs=runs,
        host=host,
        complete_only=complete_only,
        released_after=released_after,
    )
    with stage_errors(logger):
        result = run_census(params, logger)
        if tsv is not None:
            write_tsv(tsv, result)
    if as_json:
        _echo_json(result.to_json())
    else:
        render(result)


def _glance_tool_help() -> str:
    """Dereplicators that implement the comparison capability, from the registry."""
    from ..dereplicators.base import compare_supporters

    names = compare_supporters()
    return (
        f"Dereplicator with comparison support: auto, {', '.join(names) or '(none registered)'}. "
        "auto uses dRep when it can run (on the PATH or via the container backend), "
        "otherwise sourmash."
    )


@app.command(rich_help_panel=PANEL_INSPECT)
def glance(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    tool: str = typer.Option("auto", "--tool", help=_glance_tool_help()),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    plot_max: float = typer.Option(
        1.0,
        "--plot-max",
        help="Upper bound of the ANI values plotted, as a fraction from 0 to 1.",
    ),
    plot_min: float = typer.Option(
        0.0,
        "--plot-min",
        help="Lower bound of the ANI values plotted, as a fraction from 0 to 1.",
    ),
    keep_files: bool = typer.Option(
        False,
        "--keep-files",
        help="Keep the comparison tool's working directory glance_wd/.",
    ),
) -> None:
    """Quick all-vs-all ANI overview (dendrogram + similarity plots).

    Needs dRep or sourmash (on the PATH, or via the container backend) and
    no dereplication; compares every genome in genomes/ and writes
    glance_clustering_dendrogram.pdf and two ANI plots,
    glance_MASH_ANI_similarity_boxplot.png and
    glance_MASH_ANI_similarity_histogram.png. dRep reports Mash ANI;
    sourmash reports the k-mer ANI estimate that dereplicate --tool
    sourmash thresholds.
    """
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.glance import GlanceParams, log_auto_choice, resolve_auto_tool

    def build() -> GlanceParams:
        _require_choice(tool, {"auto", *_derep_registry.names()}, "--tool")
        if not 0.0 <= plot_min <= plot_max <= 1.0:
            raise UserInputError(
                "--plot-min and --plot-max are ANI fractions with "
                f"0 <= --plot-min <= --plot-max <= 1; got {plot_min} and {plot_max}."
            )
        chosen = tool
        if tool == "auto":
            # Resolved here, before the resume fingerprint is taken, so the
            # fingerprint and the stage record name the concrete tool. When
            # nothing can run, 'auto' passes through and the stage reports it
            # after its workdir and genome checks.
            resolved = resolve_auto_tool()
            if resolved is not None:
                log_auto_choice(logging.getLogger("repgenr"), resolved)
                chosen = resolved
        return GlanceParams(
            tool=chosen,
            threads=threads,
            plot_max=plot_max,
            plot_min=plot_min,
            keep_files=keep_files,
        )

    _run("glance", workdir, build)


@app.command(name="derep-unpack", rich_help_panel=PANEL_INSPECT)
def derep_unpack(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    no_representant: bool = typer.Option(
        False, "--no-representant", help="Leave the representative out of its cluster directory."
    ),
) -> None:
    """Explode clusters into one directory per representative."""
    from ..stages.derep_unpack import DerepUnpackParams

    def build() -> DerepUnpackParams:
        return DerepUnpackParams(no_representant=no_representant)

    _run("derep_unpack", workdir, build)


@app.command(name="cluster-summary", rich_help_panel=PANEL_INSPECT)
def cluster_summary(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
) -> None:
    """Regenerate derep/cluster_summary.tsv (size, species, keeper quality per cluster).

    n_members counts the genomes under the representative and excludes it;
    n_genomes includes it. Species come from the manifest taxonomy, or from
    canonical filenames when the manifest has none, and at most five are
    listed.
    """
    from ..stages.cluster_summary import ClusterSummaryParams

    _run("cluster_summary", workdir, ClusterSummaryParams)


@app.command(name="derep-stock", rich_help_panel=PANEL_INSPECT)
def derep_stock(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    action: str = typer.Option(..., "--action", help="list, pack, unpack or delete."),
    name: str | None = typer.Option(
        None,
        "--name",
        help="Run name for pack/unpack/delete: up to 100 letters, digits, '.', '_' or '-', "
        "starting with a letter or digit.",
    ),
) -> None:
    """Store, load, list or delete named dereplication runs."""
    from ..stages.derep_stock import DerepStockParams
    from .param_builders import DEREP_STOCK_ACTIONS

    def build() -> DerepStockParams:
        _require_choice(action, DEREP_STOCK_ACTIONS, "--action")
        return DerepStockParams(action=action, name=name)

    _run("derep_stock", workdir, build)


@app.command(name="sketch", rich_help_panel=PANEL_INSPECT)
def sketch(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
) -> None:
    """Write the sourmash sketch of each genome to sketches/.

    One file per genome, sketches/<name>.sig.zip, with DNA signatures at
    k=21, 31 and 51 (scaled=1000) named after the genome. Only missing and
    stale sketches are written (stale: the FASTA or the parameters changed);
    sketches of genomes no longer in the set are removed. repgenr --force
    sketch writes every sketch again. Needs sourmash (on the PATH, or via the
    container backend).
    """
    from ..stages.sketch import SketchParams

    _run("sketch", workdir, lambda: SketchParams(threads=threads))


def _tool_label(reg: Any, name: str) -> str:
    """The name, plus the adapter's declared genome limit when it has one."""
    from ..core.plugins import _capabilities_of

    if reg.is_broken(name):
        return f"{name} (broken)"
    cap = _capabilities_of(reg, name)
    limit = None if cap is None else cap.recommended_max_genomes
    return name if limit is None else f"{name} (up to {limit} genomes)"


@app.command(name="list-tools", rich_help_panel=PANEL_ENV)
def list_tools(
    check: bool = typer.Option(
        False,
        "--check",
        help="Run each adapter's preflight and report whether its binaries are present.",
    ),
    strict: bool = typer.Option(
        False,
        "--strict",
        help="With --check, exit 4 when an adapter is missing or errored and 5 when "
        "a plugin is broken, after the full listing.",
    ),
    images: bool = typer.Option(
        False,
        "--images",
        help="With --check under a container backend, also report whether each "
        "tool's image is present locally (docker image inspect or the Singularity "
        ".sif cache; nothing is pulled).",
    ),
) -> None:
    """List the available pluggable tools in each family.

    A tool that declares a recommended scale is shown as 'name (up to N
    genomes)'; auto-selection and the scale warnings use the same limit.
    The last line names the dereplicators that glance can run.
    With --check, every adapter's required binaries are looked up (version
    floors included) and reported per tool, so an environment can be
    verified before a run without a working directory. Under a container
    backend each line names where the tool runs: '[image <ref>]' or
    '[host]'; --images adds, for each tool that passed, whether its images
    (secondary ones such as racon's minimap2 included) are present
    locally. A version query that does not answer within 8 s is stopped,
    and the version is shown as unknown unless the tool's conda package
    record names one. --check alone
    always exits 0, since a host that has only some families installed is
    normal; --check --strict exits 4 when any adapter is missing or errored,
    or 5 when any plugin failed to load, so a script can verify an
    environment.
    """
    from ..core.binaries import version_timeout
    from ..core.errors import MissingBinaryError, PluginError
    from ..core.logging import configure_logging

    if strict and not check:
        raise typer.BadParameter("--strict needs --check.", param_hint="--strict")
    if images and not check:
        raise typer.BadParameter("--images needs --check.", param_hint="--images")
    from ..core.containers import get_config

    if images and not get_config().active:
        raise typer.BadParameter(
            "--images needs a container backend (--container docker or singularity).",
            param_hint="--images",
        )
    # A plugin that fails to load warns through the repgenr logger; give the
    # line the standard timestamp and level.
    configure_logging(None, level=min(_RUN_STATE["log_level"], logging.WARNING))
    from ..aligners.base import registry as aligners
    from ..annotators.base import registry as annotators
    from ..assemblers.base import registry as assemblers
    from ..classifiers.base import registry as classifiers
    from ..coregenes.base import registry as coregenes
    from ..dereplicators.base import registry as dereplicators
    from ..maskers.base import registry as maskers
    from ..polishers.base import registry as polishers
    from ..snptypers.base import registry as snptypers
    from ..treebuilders.base import registry as treebuilders

    statuses: set[str] = set()
    for label, reg in (
        ("dereplicators", dereplicators),
        ("aligners", aligners),
        ("snptypers", snptypers),
        ("maskers", maskers),
        ("treebuilders", treebuilders),
        ("assemblers", assemblers),
        ("classifiers", classifiers),
        ("polishers", polishers),
        ("annotators", annotators),
        ("coregenes", coregenes),
    ):
        entries = [_tool_label(reg, name) for name in reg.names()]
        typer.echo(f"{label}: {', '.join(entries) or '(none)'}")
        if not check:
            continue
        for name in reg.names():
            with version_timeout(_CHECK_VERSION_TIMEOUT):
                state, text = _preflight_summary(reg, name, images=images)
            statuses.add(state)
            typer.echo(f"  {name}: {text}")
    from ..dereplicators.base import compare_supporters

    # glance is not a family of its own: it runs any dereplicator with compare().
    typer.echo(
        f"glance (dereplicators with compare): {', '.join(compare_supporters()) or '(none)'}"
    )
    if strict:
        if "broken" in statuses:
            raise typer.Exit(code=PluginError.exit_code)
        if statuses & {"missing", "error"}:
            raise typer.Exit(code=MissingBinaryError.exit_code)


# Seconds each version query may take under `list-tools --check`: a listing
# of every adapter should not wait the full stage timeout on a tool that hangs.
_CHECK_VERSION_TIMEOUT = 8.0


def _one_line(exc: Exception) -> str:
    """Flatten a multi-line preflight message to one line, dropping its preamble."""
    text = str(exc).replace("Required external tools are missing or outdated:", "")
    return "; ".join(part.strip() for part in text.splitlines() if part.strip())


def _preflight_summary(reg, name: str, *, images: bool = False) -> tuple[str, str]:
    """Status and one line per adapter for `list-tools --check`.

    The status is one of ok, missing, error and broken; the line is ok with
    versions, or why not. Under a container backend the line names where the
    tool runs, '[image <ref>]' or '[host]'. With ``images``, an ok line also
    says whether each image the adapter recorded is present locally, its
    secondary images (racon's minimap2) included.
    """
    from ..core.errors import MissingBinaryError, RepGenRError

    if reg.is_broken(name):
        return "broken", f"broken (failed to load: {_one_line(reg.load_error(name))})"
    try:
        versions = reg.create(name).preflight()
    except MissingBinaryError as exc:
        return "missing", f"missing{_where(reg, name)[0]} ({_one_line(exc)})"
    except RepGenRError as exc:
        return "error", f"error{_where(reg, name)[0]} ({_one_line(exc)})"
    except Exception as exc:  # a third-party adapter must not end the listing
        reason = f"{type(exc).__name__}: {_one_line(exc)}"
        return "error", f"error{_where(reg, name)[0]} ({reason})"
    label, image = _where(reg, name)
    shown = dict(versions)
    if image is not None and shown.get(name) == image:
        # The label names the image; the versions keep the engine.
        del shown[name]
        if images:
            parts = [f"image {image}, {_presence(image)}"]
            for tool, ref in sorted(shown.items()):
                if _looks_like_image(ref):
                    parts.append(f"{tool} image {ref}, {_presence(ref)}")
                    del shown[tool]
            label = f" [{'; '.join(parts)}]"
    text = ", ".join(f"{k} {v}" for k, v in sorted(shown.items())) or "no binaries declared"
    return "ok", f"ok{label} ({text})"


def _looks_like_image(value: str) -> bool:
    """A recorded version that is an image reference (a registry path or a .sif)."""
    return "/" in value or value.endswith(".sif")


def _presence(image: str) -> str:
    from ..core.containers import get_config, image_present

    present = image_present(image, get_config())
    return {True: "present", False: "not pulled", None: "presence unknown"}[present]


def _where(reg, name: str) -> tuple[str, str | None]:
    """Where a tool runs under a container backend: the label and the image.

    The label is ' [image <ref>]' or ' [host]'; without a backend it is empty.
    """
    from ..core.containers import get_config, resolve_image
    from ..core.plugins import _capabilities_of

    config = get_config()
    caps = _capabilities_of(reg, name)
    if not config.active or caps is None:
        return "", None
    try:
        image = resolve_image(caps, config)
    except Exception:  # a failed Wave build is already in the line's reason
        return "", None
    if image is None:
        return " [host]", None
    return f" [image {image}]", image
