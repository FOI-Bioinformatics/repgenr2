"""Stage: annotate the genome set once, into a cache the core-gene route reads.

``run`` annotates the representatives (or every genome with
``all_genomes``) and the staged outgroup with one annotator (bakta by
default, prokka without a database), keeps the results under
``annotations/<record>/`` with an index, and on a later run annotates only
what is new or stale (a changed genome, annotator, tool version, database
version or extras); records no longer in the set are removed.
``annotate_core`` is the stateless per-genome step the data-channel
``annotate-build`` command and Nextflow use.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path

from ..annotators.base import AnnotateParams, AnnotateResult, Annotator
from ..annotators.base import registry as annotator_registry
from ..core import annotations as cache
from ..core.annotations import AnnotationRecord
from ..core.context import WorkdirContext
from ..core.contracts import record_name
from ..core.errors import UserInputError
from ..core.plugins import warn_unconsumed_extras
from ..core.process import remove_tree, stage_plain_inputs
from .phylo import _resolve_outgroup
from .snptype import _check_inputs, _genome_set


@dataclass
class AnnotateStageParams:
    tool: str = "bakta"
    threads: int = 16
    all_genomes: bool = False
    # bakta's database directory (flag), else BAKTA_DB; prokka needs none.
    bakta_db: str | None = None
    allow_incomplete: bool = False
    extra: dict = field(default_factory=dict)


@dataclass
class AnnotateBuildParams:
    """The data-channel step: one genome into ``out_dir/<record>/``."""

    genome: Path
    out_dir: Path
    record_name: str | None = None
    tool: str = "bakta"
    bakta_db: str | None = None
    threads: int = 16
    extra: dict = field(default_factory=dict)
    versions_out: Path | None = None


@dataclass
class AnnotateSummary:
    annotated: int
    reused: int
    removed: int
    index: Path


def _database(tool: str, annotator: Annotator, given: str | None) -> Path | None:
    if not annotator.requires_database:
        return None
    from ..annotators.bakta import bakta_db_from_env, check_bakta_db

    raw = given or bakta_db_from_env()
    return check_bakta_db(Path(raw).expanduser() if raw else None)


def annotate_core(
    genome: Path,
    name: str,
    out_dir: Path,
    scratch: Path,
    tool: str,
    params: AnnotateParams,
    logger: logging.Logger,
    *,
    annotator: Annotator | None = None,
    versions: dict[str, str] | None = None,
) -> tuple[AnnotateResult, dict[str, str]]:
    """Annotate one genome into ``out_dir`` (stateless; no config).

    Gzipped input is decompressed into ``scratch`` for a tool that does not
    read it; ``out_dir`` is written through a partial directory and renamed
    into place. Returns the result and the tool versions.
    """
    annotator = annotator or annotator_registry.create(tool)
    # The caller may have preflighted once for the whole set: under a
    # container backend every check is a container run.
    if versions is None:
        versions = annotator.preflight()
        if params.db is not None:
            versions[f"{tool}-db"] = annotator.database_version(params.db)
    parent = out_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.partial_dir(parent, out_dir.name)
    remove_tree(tmp)
    inputs = scratch / "inputs"
    remove_tree(inputs)
    try:
        staged = stage_plain_inputs([genome], annotator.capabilities, inputs, logger)
        result = annotator.call(staged[genome], name, tmp, params, logger)
    except BaseException:
        remove_tree(tmp)
        raise
    finally:
        remove_tree(inputs)
    final = cache.publish(parent, out_dir.name)
    return (
        AnnotateResult(
            gff=final / result.gff.name,
            faa=final / result.faa.name if result.faa else None,
            ffn=final / result.ffn.name if result.ffn else None,
        ),
        versions,
    )


def annotate_build(params: AnnotateBuildParams, logger: logging.Logger) -> AnnotateResult:
    """Annotate one genome from explicit paths (stateless; no config)."""
    from ..core.versions import write_versions_fragment

    if params.tool not in annotator_registry.names():
        raise UserInputError(
            f"Unknown annotator '{params.tool}'. Choose from: "
            f"{', '.join(annotator_registry.names())}."
        )
    if not params.genome.is_file():
        raise UserInputError(f"Genome not found: {params.genome}")
    annotator = annotator_registry.create(params.tool)
    warn_unconsumed_extras(annotator.capabilities, params.extra, logger, family="Annotator")
    db = _database(params.tool, annotator, params.bakta_db)
    name = params.record_name or record_name(params.genome)
    result, versions = annotate_core(
        params.genome,
        name,
        params.out_dir / name,
        params.out_dir / ".scratch",
        params.tool,
        AnnotateParams(threads=params.threads, db=db, extra=dict(params.extra)),
        logger,
        annotator=annotator,
    )
    remove_tree(params.out_dir / ".scratch")
    if params.versions_out is not None:
        write_versions_fragment(params.versions_out, versions)
    logger.info("Annotated %s into %s", name, result.gff.parent)
    return result


def run(ctx: WorkdirContext, params: AnnotateStageParams) -> AnnotateSummary:
    logger = ctx.logger
    if params.tool not in annotator_registry.names():
        raise UserInputError(
            f"Unknown annotator '{params.tool}'. Choose from: "
            f"{', '.join(annotator_registry.names())}."
        )
    _check_inputs(ctx, params.all_genomes, params.allow_incomplete, logger)
    annotator = annotator_registry.create(params.tool)
    warn_unconsumed_extras(annotator.capabilities, params.extra, logger, family="Annotator")
    db = _database(params.tool, annotator, params.bakta_db)
    versions = annotator.preflight()
    db_version = annotator.database_version(db) if db is not None else ""
    if db is not None:
        versions[f"{params.tool}-db"] = db_version
    tool_version = versions.get(params.tool, "")
    extra = cache.extra_text({**annotator.capabilities.default_params, **params.extra})

    genomes = _genome_set(ctx, params.all_genomes)
    sets = [(g, cache.SET_GENOMES) for g in genomes]
    outgroup, _leaf = _resolve_outgroup(ctx, False, logger)
    if outgroup is not None:
        sets.append((outgroup, cache.SET_OUTGROUP))
    if not sets:
        raise UserInputError("No genomes to annotate; run dereplicate first (or --all-genomes).")

    directory = ctx.annotations_dir
    directory.mkdir(parents=True, exist_ok=True)
    cache.clear_partials(directory)
    records = cache.read_index(directory)
    names = [record_name(g) for g, _ in sets]
    if len(set(names)) != len(names):
        seen: dict[str, list[Path]] = {}
        for genome, _which in sets:
            seen.setdefault(record_name(genome), []).append(genome)
        clash = [
            f"{name} ({', '.join(str(p) for p in paths)})"
            for name, paths in sorted(seen.items())
            if len(paths) > 1
        ]
        raise UserInputError(
            "Two genomes share a record name, which the annotation cache keys on: "
            + "; ".join(clash)
        )
    removed = cache.remove_stale(directory, records, names, ctx.workdir, logger)
    if removed:
        cache.write_index(directory, records)

    annotated = reused = 0
    scratch = ctx.scratch_dir / "annotate"
    for genome, which in sets:
        name = record_name(genome)
        digest = cache.genome_digest(genome)
        gff = cache.gff_path(directory, name)
        if (
            cache.is_current(
                records.get(name),
                genome_digest=digest,
                annotator=params.tool,
                annotator_version=tool_version,
                db_version=db_version,
                extra=extra,
                gff=gff,
            )
            and not ctx.force
        ):
            # A genome that moved between the set and the outgroup keeps its
            # annotation; only the column that names the set is updated.
            current = records[name]
            if current.set != which:
                records[name] = replace(current, set=which)
            reused += 1
            continue
        logger.info("Annotating %s with %s", name, params.tool)
        result, _ = annotate_core(
            genome,
            name,
            cache.record_dir(directory, name),
            scratch,
            params.tool,
            AnnotateParams(threads=params.threads, db=db, extra=dict(params.extra)),
            logger,
            annotator=annotator,
            versions=dict(versions),
        )
        records[name] = AnnotationRecord(
            record_name=name,
            set=which,
            genome_path=_relative(ctx, genome),
            genome_digest=digest,
            annotator=params.tool,
            annotator_version=tool_version,
            db_version=db_version,
            extra=extra,
            gff_digest=cache.gff_digest(result.gff),
            completed=datetime.now(UTC).isoformat(),
        )
        # Written after every genome, so an interrupted run keeps what it did.
        cache.write_index(directory, records)
        annotated += 1
    remove_tree(scratch)
    index = cache.write_index(directory, records)
    logger.info(
        "Annotations: %d annotated, %d reused, %d removed (%s)",
        annotated,
        reused,
        removed,
        index,
    )
    ctx.config.record_stage(
        "annotate",
        tool=params.tool,
        params={
            "all_genomes": params.all_genomes,
            "db_version": db_version or None,
            "extra": {k: str(v) for k, v in sorted(params.extra.items())},
            "annotated": annotated,
            "reused": reused,
            "removed": removed,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return AnnotateSummary(annotated=annotated, reused=reused, removed=removed, index=index)


def _relative(ctx: WorkdirContext, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ctx.workdir.resolve()))
    except ValueError:
        return str(path)
