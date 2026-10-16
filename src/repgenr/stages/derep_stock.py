"""derep_stock stage: store/load named dereplication runs.

Ports ``derep_stocker.py`` to the new ``derep/`` contract. A packed run keeps
``clusters.tsv`` + ``genome_status.tsv``, symlinks to the representative genome
files and ``record.json``, the ``dereplicate`` stage record (tool, parameters,
tool versions) that described the run when it was packed; unpacking restores
them into the working directory.
"""

from __future__ import annotations

import json
import logging
import re
import shutil
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from ..core.config import StageRecord
from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    atomic_replace,
    list_fasta,
    read_clusters,
    read_genome_status,
    write_cluster_summary,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.process import link_or_copy, remove_tree
from ..dereplicators.base import DerepResult
from .dereplicate import _update_manifest

_FLAT_FILES = (CLUSTERS_TSV, GENOME_STATUS_TSV, CLUSTER_SUMMARY_TSV)
# The dereplicate record of the stored run, written by pack.
_RECORD_JSON = "record.json"
_RECORD_FIELDS = ("tool", "params", "tool_versions", "completed")
# Params key on the dereplicate record while an unpack replaces its outputs;
# removed when the record is re-stamped, so a later interrupted dereplicate
# run is told apart from an interrupted unpack.
_UNPACKING = "unpacking"
# One path component, well below the 255-byte file-name limit.
_NAME_MAX = 100
_NAME_RE = re.compile(rf"[A-Za-z0-9][A-Za-z0-9._-]{{0,{_NAME_MAX - 1}}}")


@dataclass
class DerepStockParams:
    action: str  # list | pack | unpack | delete
    name: str | None = None


def precheck(ctx: WorkdirContext, params: DerepStockParams) -> None:
    """Refuse a pack or unpack that cannot proceed, before anything is changed.

    The CLI harness calls this before it marks the stage record incomplete:
    one ``derep_stock`` record serves every stored run, so a mistyped name or
    an unknown run must not leave the record of the last finished pack or
    unpack looking interrupted.
    """
    if params.action not in ("pack", "unpack"):
        return
    run_path = _run_path(ctx, params)
    if params.action == "pack":
        _check_packable(ctx)
    else:
        _check_unpackable(ctx, run_path)


def run(ctx: WorkdirContext, params: DerepStockParams) -> None:
    if not ctx.workdir.is_dir():
        raise WorkdirError(f"Working directory not found: {ctx.workdir}")
    if params.action == "list":
        _list(ctx.derep_dir / "stock", ctx.logger)
        return
    run_path = _run_path(ctx, params)
    match params.action:
        case "pack":
            _pack(ctx, run_path)
        case "unpack":
            _unpack(ctx, run_path)
        case "delete":
            # Not recorded: a delete leaves nothing to resume, and the CLI runs
            # it as a query so a repeat delete is checked instead of skipped.
            _delete(run_path)
            ctx.logger.info("Deleted stored run '%s'", params.name)
            return
        case _:
            raise UserInputError(f"Unknown action '{params.action}'")
    ctx.config.record_stage(
        "derep_stock", params=asdict(params), completed=datetime.now(UTC).isoformat()
    )
    ctx.save_config()


def _run_path(ctx: WorkdirContext, params: DerepStockParams) -> Path:
    """The stored run's directory, after checking that --name is a safe name."""
    if not params.name:
        raise UserInputError("pack/unpack/delete require --name")
    # The name becomes a directory under the store and is passed to rmtree on
    # pack/delete, so it must be a plain single-component name.
    if not _NAME_RE.fullmatch(params.name):
        raise UserInputError(
            f"Invalid run name '{params.name}': use up to {_NAME_MAX} letters, digits, "
            "'.', '_' or '-', starting with a letter or digit (no path separators)"
        )
    return ctx.derep_dir / "stock" / params.name


def _list(store: Path, logger: logging.Logger) -> None:
    # The names are the command's result, so they go to stdout (one per line,
    # also under --quiet); the empty case is only a log message.
    runs = _stored_runs(store)
    if not runs:
        logger.info("No stored runs")
        return
    for name in runs:
        print(name)


def _check_packable(ctx: WorkdirContext) -> list[Path]:
    """The live representatives, or an error when there is no dereplication."""
    # Checked before touching the store: a workdir without a dereplication
    # would otherwise be stored as an empty run.
    clusters = ctx.derep_dir / CLUSTERS_TSV
    if not clusters.is_file():
        raise WorkdirError(f"Missing {clusters}. Run the dereplicate stage first.")
    reps = list_fasta(ctx.representatives_dir)
    if not reps:
        raise WorkdirError(
            f"No representative genomes under {ctx.representatives_dir}. "
            "Run the dereplicate stage first."
        )
    return reps


def _pack(ctx: WorkdirContext, run_path: Path) -> None:
    reps = _check_packable(ctx)
    if run_path.exists():
        ctx.logger.warning(
            "Replacing stored run '%s' with the current dereplication", run_path.name
        )
        remove_tree(run_path)
    run_path.mkdir(parents=True)
    for name in _FLAT_FILES:
        src = ctx.derep_dir / name
        if src.exists():
            shutil.copy2(src, run_path / name)
    reps_dir = run_path / "representatives"
    reps_dir.mkdir()
    for rep in reps:
        (reps_dir / rep.name).symlink_to((ctx.genomes_dir / rep.name).resolve())
    _store_record(ctx, run_path)
    ctx.logger.info("Packed run to %s", run_path)


def _store_record(ctx: WorkdirContext, run_path: Path) -> None:
    """Write the live dereplicate record into the stored run.

    Only a completed record describes the outputs being packed; an incomplete
    one (an interrupted dereplicate or unpack) or none at all is not stored,
    and unpack then falls back to the record that is live at unpack time.
    """
    record = ctx.config.stages.get("dereplicate")
    if record is None or not record.completed:
        ctx.logger.warning(
            "No completed dereplicate record to store with run '%s'; unpack will "
            "attribute it to the dereplicate record current at that time",
            run_path.name,
        )
        return
    data = {key: record.to_dict()[key] for key in _RECORD_FIELDS}
    with atomic_replace(run_path / _RECORD_JSON) as fo:
        fo.write(json.dumps(data, indent=2, sort_keys=True) + "\n")


def _check_record(data: Any) -> None:
    """Raise ValueError unless ``data`` has the shape pack writes."""
    if not isinstance(data, dict):
        raise ValueError("not a JSON object")
    tool = data.get("tool")
    if tool is not None and not isinstance(tool, str):
        raise ValueError("'tool' is not a string")
    if not isinstance(data.get("params", {}), dict):
        raise ValueError("'params' is not an object")
    versions = data.get("tool_versions", {})
    if not isinstance(versions, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in versions.items()
    ):
        raise ValueError("'tool_versions' is not an object of strings")


def _stored_record(ctx: WorkdirContext, run_path: Path) -> StageRecord | None:
    """The dereplicate record kept with a stored run, or None.

    Runs packed before the record was stored have no ``record.json``; an
    unreadable file is reported and treated the same way.
    """
    path = run_path / _RECORD_JSON
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        _check_record(data)
        return StageRecord.from_dict(data)
    except (ValueError, TypeError, OSError) as exc:
        ctx.logger.warning(
            "Ignoring unreadable %s (%s); the dereplicate record is carried over", path, exc
        )
        return None


def _check_unpackable(ctx: WorkdirContext, run_path: Path) -> list[str]:
    """The stored run's representative names, or an error when it is incomplete."""
    if not run_path.exists():
        stored = _stored_runs(run_path.parent)
        raise UserInputError(
            f"No stored run named '{run_path.name}'; stored runs: "
            f"{', '.join(stored) if stored else 'none'}."
        )
    # Validate the stored run in full before the current dereplication is
    # replaced, so an incomplete run leaves the workdir unchanged.
    stored_clusters = run_path / CLUSTERS_TSV
    if not stored_clusters.is_file():
        raise WorkdirError(f"Stored run '{run_path.name}' is incomplete: missing {stored_clusters}")
    stored_reps = run_path / "representatives"
    if not stored_reps.is_dir():
        raise WorkdirError(f"Stored run '{run_path.name}' is incomplete: missing {stored_reps}")
    rep_names = [rep.name for rep in list_fasta(stored_reps)]
    absent = [name for name in rep_names if not (ctx.genomes_dir / name).is_file()]
    if absent:
        raise WorkdirError(
            f"Stored run '{run_path.name}' names {len(absent)} representative(s) "
            f"not found under {ctx.genomes_dir}, e.g. {absent[0]}"
        )
    return rep_names


def _unpack(ctx: WorkdirContext, run_path: Path) -> None:
    rep_names = _check_unpackable(ctx, run_path)
    prior = ctx.config.stages.get("dereplicate")
    if prior is not None and not prior.completed and _UNPACKING not in prior.params:
        # An incomplete record describes a dereplicate run that did not
        # finish, not the stored run being restored: carry nothing over from
        # it. An interrupted unpack leaves the in-progress marker in the
        # params (set below), and a repeat keeps what that unpack carried.
        prior = None
    stored = _stored_record(ctx, run_path)
    if stored is not None:
        # The run keeps its own record: the restored record names the tool,
        # parameters and versions that produced it.
        carried: tuple[str | None, dict[str, Any], dict[str, str]] | None = (
            stored.tool,
            stored.params,
            stored.tool_versions,
        )
    else:
        # Runs packed without a record: carry over the live record, as before.
        carried = (prior.tool, dict(prior.params), dict(prior.tool_versions)) if prior else None
    if carried is not None:
        carried[1].pop(_UNPACKING, None)
    if prior is not None:
        # The dereplicate stage's outputs are replaced below: mark its record
        # incomplete first, so an unpack that stops half-way is not reported
        # as a finished dereplication.
        prior.completed = None
        prior.fingerprint = None
        prior.params = {**prior.params, _UNPACKING: run_path.name}
        ctx.save_config()
    # The summary is not restored: it is rebuilt below from the restored
    # clusters and the live manifest.
    for name in (CLUSTERS_TSV, GENOME_STATUS_TSV):
        src = run_path / name
        if src.exists():
            shutil.copy2(src, ctx.derep_dir / name)
        elif (ctx.derep_dir / name).exists():
            # Do not leave a file of the replaced dereplication beside the
            # restored ones.
            (ctx.derep_dir / name).unlink()
    if ctx.representatives_dir.exists():
        remove_tree(ctx.representatives_dir)
    ctx.representatives_dir.mkdir(parents=True)
    for name in rep_names:
        link_or_copy(ctx.genomes_dir / name, ctx.representatives_dir / name)
    # The derep contract now describes the stored run: bring the manifest's
    # per-genome status in line with it and re-stamp the dereplicate record
    # without a fingerprint, so `status` reports the run on disk and the next
    # `dereplicate` recomputes instead of skipping on a stale fingerprint.
    clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
    status_path = ctx.derep_dir / GENOME_STATUS_TSV
    genome_status = read_genome_status(status_path) if status_path.exists() else {}
    _update_manifest(
        ctx, DerepResult(representatives=[], clusters=clusters, genome_status=genome_status)
    )
    # Rebuild the summary rather than restore the stored copy: the stored one
    # carries the manifest quality of pack time (or is absent for runs packed
    # before the summary existed), while derep/cluster_summary.tsv describes
    # the live manifest, as `cluster-summary` writes it. The stored copy stays
    # in the store as the pack-time view.
    from ..core.resources import usable_cpus
    from .cluster_summary import gtdb_lookup, summarise_clusters, taxonomy_lookup
    from .derep_keeper import N50_MAX_WORKERS, N50Lookup, warn_missing_n50
    from .dereplicate import quality_lookup

    n50 = N50Lookup([ctx.genomes_dir], workers=min(usable_cpus(), N50_MAX_WORKERS))
    write_cluster_summary(
        ctx.derep_dir / CLUSTER_SUMMARY_TSV,
        summarise_clusters(
            clusters,
            quality_lookup(ctx),
            taxonomy_lookup(ctx),
            n50,
            gtdb_representatives=gtdb_lookup(ctx),
        ),
    )
    warn_missing_n50(n50, ctx.logger)
    tool, params, versions = carried if carried else (None, {}, None)
    ctx.config.record_stage(
        "dereplicate",
        tool=tool,
        params={**params, "stock": run_path.name},
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.logger.info("Unpacked run from %s", run_path)


def _stored_runs(store: Path) -> list[str]:
    if not store.is_dir():
        return []
    return sorted(p.name for p in store.iterdir() if p.is_dir())


def _delete(run_path: Path) -> None:
    if not run_path.exists():
        stored = _stored_runs(run_path.parent)
        raise WorkdirError(
            f"No stored run named '{run_path.name}'; stored runs: "
            f"{', '.join(stored) if stored else 'none'}."
        )
    remove_tree(run_path)
