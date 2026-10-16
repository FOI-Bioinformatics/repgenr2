"""Dereplication stage.

Selects a dereplicator adapter, runs it (two-stage chunking when ``--process-size``
is set and exceeded, for any tool), then normalizes the result into the canonical
contract: ``derep/representatives/`` + ``clusters.tsv`` + ``genome_status.tsv``.
Provenance (tool, params, versions) is recorded in ``repgenr.yaml`` and the
manifest derep status is updated.
"""

from __future__ import annotations

import logging
import math
import sqlite3
import time
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    REFINED_TSV,
    RefinedRow,
    accession_from_filename,
    list_fasta,
    write_cluster_summary,
    write_clusters,
    write_genome_status,
    write_refined,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.integrity import check_genome_completeness
from ..core.plugins import (
    auto_select,
    scale_warning,
    warn_ignored_params,
    warn_unconsumed_extras,
)
from ..core.process import link_or_copy, remove_tree
from ..core.resources import usable_cpus
from ..core.sketches import adapter_sketches, workdir_provider
from ..dereplicators.base import (
    STATUS_REPRESENTATIVE,
    DerepParams,
    DerepResult,
    check_result_complete,
    registry,
    run_quality,
)
from .cluster_summary import gtdb_lookup, summarise_clusters, taxonomy_lookup
from .derep_keeper import (
    N50Lookup,
    N50Of,
    choose_keeper,
    log_unpromoted,
    prefer_gtdb_representatives,
    prefetch_n50,
    rescore_representatives,
)
from .derep_refine import (
    Lineage,
    check_refine_ani,
    fine_pass,
    lineage_of,
    manifest_accessions,
    record_params,
    refined_rows,
    resolve_seeds,
    splice,
)


@dataclass
class DereplicateParams:
    tool: str = "skder"
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    # Chunk size. When set and exceeded, the two-stage chunked path runs for ANY
    # tool (native-scaling tools are single-pass by default but can be chunked
    # explicitly, e.g. to bound memory/open-file counts at 1000s-10000s genomes).
    process_size: int | None = None
    num_processes: int = 0  # parallel stage-1 chunk workers (0 = auto); threads split across them
    # Stage-1 (intra-chunk) ANI thresholds. Default to the main thresholds; set
    # them looser to avoid over-collapsing within a chunk before the stage-2 pass.
    pre_primary_ani: float | None = None
    pre_secondary_ani: float | None = None
    # Taxonomy-aware reduction applied after ANI dereplication: collapse the
    # representatives to one per taxon ("species" or "genus") using the manifest
    # taxonomy. "none" keeps the ANI result as-is.
    reduce: str = "none"
    # Target number of representatives: when > 0, search --secondary-ani to land
    # near this count (re-runs dereplication per search step).
    target_reps: int = 0
    extra: dict | None = None
    # Proceed on a partial genome set with a warning instead of refusing.
    allow_incomplete: bool = False
    # Per-cluster representative choice: re-pick by manifest CheckM quality
    # (completeness - 5 x contamination), "gtdb" to keep a GTDB species
    # representative first (quality after it), or "tool" to keep the
    # adapter's own pick.
    keeper: str = "quality"  # quality | gtdb | tool
    # Seeds of a lineage to dereplicate again at refine_ani (derep_refine).
    refine: list[str] = field(default_factory=list)
    refine_ani: float | None = None


def precheck(ctx: WorkdirContext, params: DereplicateParams) -> None:
    """Refusals checked before the harness marks the stage record incomplete.

    A missing selected genome or an empty genomes/ changes nothing on disk, so
    the record of the last finished dereplication must stay complete; the
    stage body repeats these checks for callers that bypass the harness.
    """
    if not params.allow_incomplete:
        check_genome_completeness(
            ctx.genomes_dir, ctx.workdir, logger=ctx.logger, allow_incomplete=False
        )
    genomes = _list_genomes(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(
            f"No genome FASTAs found under {ctx.genomes_dir}. Run the genome stage first."
        )
    tool = params.tool
    if tool == "auto":
        tool = auto_select(registry, len(genomes)) or "skder"
    if tool == "skder":
        from ..dereplicators.skder import check_secondary_ani

        # --target-reps searches the threshold itself, from 0.80 upwards.
        if not params.target_reps:
            check_secondary_ani(params.secondary_ani)
        if params.pre_secondary_ani is not None:
            check_secondary_ani(params.pre_secondary_ani, "--pre-secondary-ani")
        screen = (params.extra or {}).get("screen")
        if screen is not None:
            from ..dereplicators.skder import _as_percent, resolve_screen

            # Against the final cutoff; a chunk or search step below it falls
            # back to the default rule in the adapter. --target-reps searches
            # the cutoff itself, so only the value's form is checked there.
            cutoff = 100.0 if params.target_reps else float(_as_percent(params.secondary_ani))
            resolve_screen(screen, cutoff)


def run(ctx: WorkdirContext, params: DereplicateParams) -> DerepResult:
    logger = ctx.logger
    from ..dereplicators.skder import reset_screen_log

    reset_screen_log()
    check_genome_completeness(
        ctx.genomes_dir,
        ctx.workdir,
        logger=logger,
        allow_incomplete=params.allow_incomplete,
    )
    genomes = _list_genomes(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(
            f"No genome FASTAs found under {ctx.genomes_dir}. Run the genome stage first."
        )

    tool = params.tool
    if tool == "auto":
        tool = auto_select(registry, len(genomes)) or "skder"
        logger.info("Auto-selected dereplicator '%s' for %d genomes", tool, len(genomes))
    else:
        warn = scale_warning(registry, tool, len(genomes))
        if warn:
            limit, alts = warn
            logger.warning(
                "Dereplicator '%s' is tuned for <=%d genomes but you have %d; consider: %s",
                tool,
                limit,
                len(genomes),
                ", ".join(alts) or "none",
            )

    adapter = registry.create(tool)
    versions = adapter.preflight()
    logger.info("Dereplicating %d genomes with %s", len(genomes), tool)

    caps = adapter.capabilities
    extra = {**caps.default_params, **(params.extra or {})}
    warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")
    # Manifest quality is an adapter input whichever keeper rule applies
    # (dRep --genomeInfo, galah --genome-info). Whether the adapters get it is
    # decided once for the run, so all chunks and the merge use one source.
    manifest_quality = quality_lookup(ctx)
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra=extra,
        quality=run_quality(manifest_quality, [g.name for g in genomes], logger, source="manifest"),
    )
    warn_ignored_params(caps, derep_params, logger, family="Dereplicator")
    # An adapter that compares sourmash sketches reads those of the workdir;
    # missing and stale ones are written to sketches/ first.
    derep_params.sketches = adapter_sketches(
        adapter,
        extra,
        genomes,
        workdir_provider(ctx, logger, params.threads, "dereplicate"),
        logger,
        "dereplicate",
    )

    scratch = ctx.scratch_dir / "dereplicate"
    if scratch.exists():
        remove_tree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    if params.target_reps and params.target_reps > 0:
        result = _search_target_reps(
            adapter, genomes, scratch, derep_params, params, params.target_reps, logger
        )
    else:
        result = _dereplicate_to_result(adapter, genomes, scratch, derep_params, params, logger)

    # N50 per genome for the keeper score and the summary, each genome read once.
    n50 = N50Lookup([ctx.genomes_dir], workers=params.threads)
    lineage = None
    if params.refine_ani is not None and not params.refine:
        logger.warning(
            "--refine-ani %s is given without --refine seeds; ignored", params.refine_ani
        )
    if params.refine:
        lineage, result = _refine(
            ctx, adapter, result, genomes, scratch, derep_params, params, logger
        )
    result, keeper_swaps, keeper_effective, quality, gtdb = apply_keeper(
        ctx, result, params.keeper, manifest_quality, n50, logger
    )
    # The refined clusters share a species: --reduce must leave them apart.
    refined_reps = (
        {rep for rep in result.clusters if rep in lineage.coarse_rep_of} if lineage else set()
    )

    if params.reduce != "none":
        before = len(result.representatives)
        # Reuse the quality lookup above instead of re-querying the manifest.
        result = _reduce_by_taxonomy(
            result,
            params.reduce,
            quality,
            logger,
            taxon_of=_taxon_lookup(ctx, params.reduce),
            n50=n50,
            prefer=gtdb,
            keep_apart=refined_reps,
        )
        logger.info(
            "Taxonomy reduction (one per %s): %d -> %d representatives",
            params.reduce,
            before,
            len(result.representatives),
        )

    check_result_complete(result, [g.name for g in genomes])
    rows = refined_rows(result, lineage) if lineage else []
    # The summary reports quality whichever keeper rule was used, so a
    # --keeper tool run still shows where a member outscores the keeper.
    _write_contract(ctx, result, manifest_quality, n50, refined=rows)
    _update_manifest(ctx, result)

    ctx.config.record_stage(
        "dereplicate",
        tool=tool,
        params={
            "requested_tool": params.tool,
            "primary_ani": params.primary_ani,
            "secondary_ani": params.secondary_ani,
            "aligned_fraction": params.aligned_fraction,
            "process_size": params.process_size,
            "num_processes": params.num_processes,
            "pre_primary_ani": params.pre_primary_ani,
            "pre_secondary_ani": params.pre_secondary_ani,
            "reduce": params.reduce,
            "target_reps": params.target_reps,
            "keeper": params.keeper,
            "keeper_effective": keeper_effective,
            "keeper_swaps": keeper_swaps,
            # The tool arguments, also merged below as before; derep-refine
            # reads them from here.
            "tool_args": dict(params.extra or {}),
            **(record_params(lineage, params.refine_ani or 0.0, rows) if lineage else {}),
            **(params.extra or {}),
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Dereplication complete: %d representatives of %d genomes",
        len(result.representatives),
        len(genomes),
    )
    return result


def _list_genomes(genomes_dir: Path) -> list[Path]:
    return list_fasta(genomes_dir)


_MAX_REDUCE_DEPTH = 50  # termination backstop; the shrink-guard stops much sooner
# Warn when the final pass of a chunked run compares more than this many chunks' worth.
_FINAL_PASS_WARN_FACTOR = 4
_MAX_TARGET_ITERS = 12  # binary-search steps for --target-reps


def _resolve_pre_thresholds(params: DereplicateParams, secondary: float) -> tuple[float, float]:
    """Stage-1 (intra-chunk) ANI: explicit --pre-* if set, else the effective values."""
    pre_primary = (
        params.pre_primary_ani if params.pre_primary_ani is not None else params.primary_ani
    )
    pre_secondary = params.pre_secondary_ani if params.pre_secondary_ani is not None else secondary
    return pre_primary, pre_secondary


def _dereplicate_to_result(
    adapter,
    genomes: list[Path],
    scratch: Path,
    derep_params: DerepParams,
    params: DereplicateParams,
    logger,
) -> DerepResult:
    """Produce a DerepResult: chunked (when --process-size is set and exceeded) or
    a single pass. The chunk decision is independent of the search/reduce layers."""
    needs_chunking = (
        params.process_size is not None
        and params.process_size > 0
        and len(genomes) > params.process_size
    )
    if not needs_chunking:
        return adapter.dereplicate(genomes, scratch, derep_params, logger)
    if params.process_size is None:  # guaranteed by needs_chunking
        raise RuntimeError("Chunked dereplication reached without --process-size set.")

    workers = (
        params.num_processes
        if params.num_processes > 0
        else _auto_num_processes(derep_params.threads)
    )
    pre_primary, pre_secondary = _resolve_pre_thresholds(params, derep_params.secondary_ani)
    native = adapter.capabilities.supports_native_scaling
    logger.info(
        "Chunking %d genomes at size %d (%s)%s",
        len(genomes),
        params.process_size,
        "native-scaling tool, explicit chunking" if native else "non-native-scaling tool",
        ""
        if (pre_primary, pre_secondary) == (derep_params.primary_ani, derep_params.secondary_ani)
        else f"; stage-1 ANI pri={pre_primary} sec={pre_secondary}",
    )
    return _dereplicate_chunked(
        adapter,
        genomes,
        scratch,
        derep_params,
        params.process_size,
        workers,
        pre_primary,
        pre_secondary,
        logger,
    )


def apply_keeper(
    ctx: WorkdirContext,
    result: DerepResult,
    keeper: str,
    manifest_quality: Mapping[str, tuple[float, float]],
    n50: N50Of,
    logger: logging.Logger,
) -> tuple[DerepResult, int, str, dict[str, tuple[float, float]], set[str]]:
    """Apply the ``--keeper`` rule to every cluster of ``result``.

    Returns the result, the number of clusters whose representative changed,
    the rule that took effect (``gtdb``, ``quality`` or ``tool``), the
    quality map of the clustered genomes and the GTDB representatives among
    them (both empty when the rule did not apply), which ``--reduce`` reuses.
    """
    swaps = 0
    quality: dict[str, tuple[float, float]] = {}
    gtdb: set[str] = set()
    # Only the clustered genomes count: the outgroup (usually a GTDB
    # representative itself), genomes no longer under genomes/ and genomes
    # rejected on QC re-pick nothing and must not make the record claim a
    # rule applied.
    names = {n for rep, members in result.clusters.items() for n in (rep, *members)}
    if keeper == "gtdb":
        quality = {n: q for n, q in manifest_quality.items() if n in names}
        gtdb = gtdb_lookup(ctx) & names
        if not gtdb:
            logger.warning(
                "No clustered genome is flagged as a GTDB species "
                "representative in the manifest, so "
                "--keeper gtdb falls back to %s. Only GTDB selections from the "
                "metadata stage carry the flag; a workdir selected by an older "
                "version needs `repgenr metadata --force` to record it.",
                "the quality rule" if quality else "the tool's own representatives",
            )
        result, swaps = prefer_gtdb_representatives(result, gtdb, quality, logger, n50)
    elif keeper == "quality":
        quality = {n: q for n, q in manifest_quality.items() if n in names}
        if quality:
            result, swaps = rescore_representatives(result, quality, logger, n50)
        else:
            logger.warning(
                "No assembly quality in the manifest for the genomes under genomes/, "
                "so --keeper quality has no effect: keeping the tool's own "
                "representatives. GTDB selections carry CheckM quality, assemble "
                "--checkm2-db scores assemblies, and ingest --selection reads the "
                "completeness and contamination columns."
            )
    effective = "gtdb" if gtdb else "quality" if quality else "tool"
    return result, swaps, effective, quality, gtdb


def _refine(
    ctx: WorkdirContext,
    adapter,
    result: DerepResult,
    genomes: Sequence[Path],
    scratch: Path,
    derep_params: DerepParams,
    params: DereplicateParams,
    logger: logging.Logger,
) -> tuple[Lineage, DerepResult]:
    """``--refine``: dereplicate the seeds' clusters again at ``refine_ani``."""
    if params.refine_ani is None:
        raise UserInputError("--refine needs --refine-ani.")
    check_refine_ani(params.refine_ani, params.secondary_ani)
    accessions, outgroup = manifest_accessions(ctx)
    seeds = resolve_seeds(params.refine, genomes, accessions=accessions, outgroup=outgroup)
    lineage = lineage_of(result, seeds, logger)
    if not lineage.genomes:
        logger.warning("Every --refine seed is alone in its cluster; nothing was refined")
        return lineage, result
    paths = {g.name: g for g in genomes}
    sub = fine_pass(
        adapter, lineage, paths, scratch / "refine", derep_params, params.refine_ani, logger
    )
    return lineage, splice(result, lineage, sub, paths, logger)


def _auto_num_processes(threads: int) -> int:
    """Default parallel chunk-worker count when --num-processes is 0 (auto).

    Benchmarks (1256 genomes, sourmash) show running chunks in parallel -- each at
    threads/workers threads, so the total thread budget is unchanged -- is ~2.2x
    faster at 4 workers and ~3x at 8 vs serial chunks, because a single tool run
    doesn't use many threads efficiently on one chunk. Target ~4 threads per
    worker (a balanced point, leaving per-chunk threading for larger genomes),
    capped by the CPUs this process may use (affinity mask and cgroup CPU
    quota); never exceed the thread budget (so threads/worker stays >= 1 and
    nothing is oversubscribed).
    """
    cores = usable_cpus()
    return max(1, min(threads // 4, cores))


def _search_target_reps(
    adapter,
    genomes: list[Path],
    scratch: Path,
    derep_params: DerepParams,
    params: DereplicateParams,
    target: int,
    logger,
) -> DerepResult:
    """Binary-search --secondary-ani to land near ``target`` representatives.

    Representative count is monotonic non-decreasing in the ANI threshold (a
    stricter threshold keeps more reps), so we bisect: too few reps -> raise the
    threshold, too many -> lower it. Each step re-runs dereplication; capped at
    ``_MAX_TARGET_ITERS``. Returns the result whose count is closest to target.
    """
    lo, hi = 0.80, 0.9999
    best: DerepResult | None = None
    best_diff = None
    best_ani = None
    seen: set[float] = set()
    # Sketches depend on the genome set + ksize/scaled, not the ANI threshold, so
    # share one sketch dir across iterations: adapters that support it (sourmash)
    # sketch once and every iteration reuses it instead of re-sketching.
    # The adapter keys cache entries to the requested genome set (per-genome
    # signature matching, digest-named zips), so the shared dir stays correct
    # even when an iteration dereplicates in parallel chunks.
    sketch_cache = scratch / "target" / "sketches"
    for i in range(_MAX_TARGET_ITERS):
        mid = round((lo + hi) / 2, 5)
        if mid in seen:
            break
        seen.add(mid)
        dp = replace(
            derep_params,
            secondary_ani=mid,
            extra={**derep_params.extra, "sketch_cache": str(sketch_cache)},
        )
        iter_dir = scratch / "target" / f"iter{i}"
        res = _dereplicate_to_result(adapter, genomes, iter_dir, dp, params, logger)
        n = len(res.representatives)
        diff = abs(n - target)
        logger.info(
            "target-reps search [%d]: secondary-ani=%.5f -> %d reps (target %d)",
            i + 1,
            mid,
            n,
            target,
        )
        if best_diff is None or diff < best_diff:
            best, best_diff, best_ani = res, diff, mid
        if n == target:
            break
        if n < target:
            lo = mid  # too few reps -> stricter threshold
        else:
            hi = mid  # too many reps -> looser threshold
    if best is None:  # the search loop always runs at least one iteration
        raise RuntimeError("target-reps search finished without evaluating any threshold.")
    achieved = len(best.representatives)
    if achieved != target:
        logger.warning(
            "target-reps: requested %d representatives but the closest "
            "achievable count is %d (secondary-ani=%.5f). On low-diversity "
            "sets the count is a step function of the threshold and "
            "intermediate targets cannot be reached.",
            target,
            achieved,
            best_ani,
        )
    logger.info(
        "target-reps: chose secondary-ani=%.5f -> %d representatives (target %d)",
        best_ani,
        achieved,
        target,
    )
    return best


def _recurse_min_shrink(n_input: int, process_size: int) -> int:
    """Genomes a reduction level must remove before another level is attempted.

    A level removes only the close pairs that fall into one chunk. When it
    removes fewer than a chunk's worth (or fewer than a tenth of its input),
    another level of chunks finds few more, and the all-against-all comparison
    that the merge needs cannot be avoided; the final pass is cheaper.
    """
    return max(process_size, math.ceil(n_input / 10))


def _level_pays_above(n_union: int, process_size: int) -> float:
    """Fraction of its genomes a level over ``n_union`` must remove to pay for itself.

    A level over a union of U in chunks of S costs the tool's per-genome work
    for all U plus about S / U of the final pass's pairwise work; removing a
    fraction f of U saves about 2 f of that pairwise work (and f of the
    per-genome work of the final pass). With the pairwise work about U / S
    times the per-genome work, the level breaks even near f = S / U.
    """
    return process_size / n_union


def _final_pass_work_factor(n_union: int, n_level_input: int, process_size: int) -> float:
    """Pairwise work of the final pass over ``n_union`` relative to the level just run.

    A level over n genomes in chunks of S compares about n * S / 2 pairs; one
    pass over the union compares U * U / 2. The ratio U * U / (n * S) says how
    many times the level's pairwise work the final pass is.
    """
    return n_union * n_union / (n_level_input * process_size)


def _format_factor(x: float) -> str:
    """Two significant digits below 10, a whole number above (never scientific notation)."""
    return f"{x:.2g}" if x < 10 else f"{x:.0f}"


def _format_duration(seconds: float) -> str:
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 120 * 60:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def _interleave_by_origin(stage1: list[DerepResult]) -> list[Path]:
    """The union of stage-1 representatives, taken round-robin over the chunks.

    Element 0 of every chunk result, then element 1 of every chunk result, and
    so on. When the union is chunked again, each new chunk holds representatives
    from many previous chunks, so pairs that a fixed order kept apart can meet.
    The order depends only on the input order (no random seed).
    """
    per_chunk = [list(r.representatives) for r in stage1]
    longest = max((len(reps) for reps in per_chunk), default=0)
    return [reps[i] for i in range(longest) for reps in per_chunk if i < len(reps)]


def _dereplicate_chunked(
    adapter,
    genomes: list[Path],
    scratch: Path,
    params: DerepParams,
    process_size: int,
    num_processes: int,
    pre_primary_ani: float,
    pre_secondary_ani: float,
    logger,
    depth: int = 0,
    first_chunk: DerepResult | None = None,
) -> DerepResult:
    """Hierarchical chunked dereplication, usable for any tool.

    ``first_chunk`` is the result of this level's first chunk when the caller
    already ran it as a probe; that chunk is not run again.

    Leaf chunks of <= ``process_size`` are dereplicated independently (in
    parallel, with the stage-1 ANI thresholds). The union of their
    representatives is then reduced:

    * when it is at most two chunks' worth, or when this level removed fewer
      genomes than ``_recurse_min_shrink`` requires, one final pass dereplicates the whole
      union with the final thresholds, whatever its size. The chunk size bounds
      the leaf chunks, not this merge: for a diverse set the union can be most
      of the input.
    * otherwise one chunk of the next level is run first as a probe: the
      first chunk of the union interleaved by chunk of origin, one sample of
      that level (the first representatives of every chunk). When the probe
      removes more than ``_level_pays_above`` of its genomes, the rest of that
      level runs (the probe's result is its first chunk) and the level's union
      is reduced the same way; when it removes less, the union goes to the
      final pass and the probe's result is dropped. Only the leaf chunks of
      level 0 use the stage-1 (``--pre-*``) thresholds; the probe, every later
      level and the final pass use the final thresholds.

    Why a probe: a level over a union of U representatives in chunks of S
    costs the per-genome work of the tool (sketching, N50) for all U plus
    about S / U of the final pass's pairwise work, and removes only the close
    pairs that one chunk happens to hold. How much that is depends on the
    set: a reported 6888-genome run removed 45 percent at level 0 and 5, 1.4
    and 0.5 percent at the next levels, each level costing as much as a
    level-0 wave, while clusters spread a few members per leaf chunk collapse
    almost entirely at the next level. The previous level's removal does not
    predict the next one's, so one chunk of it is measured. The probe costs
    one tool call over S genomes with every thread.

    ``_MAX_REDUCE_DEPTH`` is a termination backstop only. Membership is
    composed level by level so every original genome maps to one final
    representative or keeps its QC failure.
    """
    chunks = [genomes[i : i + process_size] for i in range(0, len(genomes), process_size)]
    if len(chunks) > 1 and len(chunks[-1]) == 1:
        chunks[-2].extend(chunks[-1])
        chunks.pop()

    if len(chunks) == 1:
        # Fits in a single pass: dereplicate directly with the final thresholds.
        if depth > 0:
            # Not reached from the probe path (a union of process_size + 1 goes
            # to the final pass without a level), kept as a guard: a union that
            # collapsed into one chunk makes this call the final pass.
            logger.info(
                "Final pass over %d representatives (chunk size %d bounds the leaf "
                "chunks, not the merge)",
                len(chunks[0]),
                process_size,
            )
        return adapter.dereplicate(chunks[0], scratch / f"level{depth}", params, logger)

    workers = max(1, min(num_processes, len(chunks)))
    threads_per_worker = max(1, params.threads // workers)
    chunk_params = replace(
        params,
        threads=threads_per_worker,
        primary_ani=pre_primary_ani,
        secondary_ani=pre_secondary_ani,
    )
    logger.info(
        "Level %d: %d chunks, %d parallel worker(s) at %d threads each",
        depth,
        len(chunks),
        workers,
        threads_per_worker,
    )

    def run_chunk(indexed: tuple[int, list[Path]]) -> DerepResult:
        idx, chunk = indexed
        chunk_dir = scratch / f"level{depth}" / f"chunk{idx}"
        logger.info("Level %d chunk %d/%d (%d genomes)", depth, idx + 1, len(chunks), len(chunk))
        return adapter.dereplicate(chunk, chunk_dir, chunk_params, logger)

    level_started = time.monotonic()
    if first_chunk is not None:
        # The probe ran chunk 0 (the same genomes, same thresholds); the rest
        # of the level runs now.
        if len(chunks[0]) != len(first_chunk.genome_status):
            raise RuntimeError("The probe chunk does not match this level's first chunk.")
        rest = parallel_map(run_chunk, list(enumerate(chunks))[1:], workers, logger=logger)
        stage1_results = [first_chunk, *rest]
    else:
        stage1_results = parallel_map(run_chunk, list(enumerate(chunks)), workers, logger=logger)
    level_seconds = time.monotonic() - level_started
    # Chunk order for the final pass (as the Nextflow merge step), interleaved
    # only when the union is chunked again.
    union = [rep for r in stage1_results for rep in r.representatives]
    merge_dir = scratch / f"level{depth}" / "merge"
    removed = len(genomes) - len(union)
    min_shrink = _recurse_min_shrink(len(genomes), process_size)
    logger.info(
        "Level %d: %d genomes -> union of %d representatives (%d removed) in %s",
        depth,
        len(genomes),
        len(union),
        removed,
        _format_duration(level_seconds),
    )

    probe: DerepResult | None = None
    # A union of at most two chunks' worth goes to the final pass without a
    # probe: the level could pay only above S / U, one half or more, and a
    # failed probe would cost nearly as much as that final pass.
    if len(union) > 2 * process_size and removed >= min_shrink and depth < _MAX_REDUCE_DEPTH:
        interleaved = _interleave_by_origin(stage1_results)
        probe_chunk = interleaved[:process_size]
        pays_above = _level_pays_above(len(union), process_size)
        logger.info(
            "Level %d removed %d (>= %d); probing level %d with one chunk of %d of the "
            "%d representatives, interleaved by chunk of origin",
            depth,
            removed,
            min_shrink,
            depth + 1,
            len(probe_chunk),
            len(union),
        )
        probe_started = time.monotonic()
        probe = adapter.dereplicate(
            probe_chunk, merge_dir / f"level{depth + 1}" / "chunk0", params, logger
        )
        probe_removed = len(probe_chunk) - len(probe.representatives)
        probe_fraction = probe_removed / len(probe_chunk)
        logger.info(
            "Probe: %d of %d removed (%.1f%%) in %s; a level over %d is expected to pay "
            "above %.1f%%",
            probe_removed,
            len(probe_chunk),
            100 * probe_fraction,
            _format_duration(time.monotonic() - probe_started),
            len(union),
            100 * pays_above,
        )
        if probe_fraction <= pays_above:
            logger.info(
                "Level %d: a further level is expected to cost more than it saves; "
                "the probe's result is dropped",
                depth,
            )
            probe = None

    if probe is not None:
        logger.info(
            "Level %d: chunking the union of %d again (level %d)", depth, len(union), depth + 1
        )
        stage2 = _dereplicate_chunked(
            adapter,
            interleaved,
            merge_dir,
            params,
            process_size,
            num_processes,
            # Levels over representatives compare them at the final thresholds;
            # only the leaf chunks of the original genomes use --pre-*.
            params.primary_ani,
            params.secondary_ani,
            logger,
            depth + 1,
            first_chunk=probe,
        )
    else:
        if len(union) > 2 * process_size and depth >= _MAX_REDUCE_DEPTH:
            logger.info("Level %d: depth limit %d reached", depth, _MAX_REDUCE_DEPTH)
        elif len(union) > 2 * process_size and removed < min_shrink:
            logger.info(
                "Level %d removed %d genomes (another level needs >= %d)",
                depth,
                removed,
                min_shrink,
            )
        logger.info(
            "Final pass over %d representatives (chunk size %d bounds the leaf "
            "chunks, not the merge)",
            len(union),
            process_size,
        )
        if len(union) > process_size:
            # The tool prints no progress that reaches the console; say what to
            # expect from the level just measured. Pairwise work grows with the
            # square of the union; the per-genome work (sketches, N50) is paid
            # once more for every genome of the union.
            logger.info(
                "One call over %d genomes with %d threads: about %s times the pairwise "
                "work of level %d, which took %s with %d worker(s); no progress is "
                "printed while it runs",
                len(union),
                params.threads,
                _format_factor(_final_pass_work_factor(len(union), len(genomes), process_size)),
                depth,
                _format_duration(level_seconds),
                workers,
            )
        if len(union) > _FINAL_PASS_WARN_FACTOR * process_size:
            logger.warning(
                "The final dereplication pass compares %d representatives, more than "
                "%d times the chunk size %d; its time and memory follow the union, "
                "not the chunk size. A larger --process-size, or --pre-secondary-ani "
                "below the final --secondary-ani (so the chunks remove more), "
                "reduces this pass.",
                len(union),
                _FINAL_PASS_WARN_FACTOR,
                process_size,
            )
        final_started = time.monotonic()
        stage2 = adapter.dereplicate(union, merge_dir, params, logger)
        logger.info(
            "Final pass: %d representatives -> %d in %s",
            len(union),
            len(stage2.representatives),
            _format_duration(time.monotonic() - final_started),
        )

    return _compose_two_stage(stage1_results, stage2, logger)


def _compose_two_stage(
    stage1: list[DerepResult], stage2: DerepResult, logger: logging.Logger | None = None
) -> DerepResult:
    from ..dereplicators.base import STATUS_CONTAINED, STATUS_FAIL_QC, STATUS_REPRESENTATIVE

    # Each stage-1 representative -> every original genome in its cluster (the rep
    # itself plus its members). Built once; lets composition run in O(N) instead
    # of scanning all members for every final rep (which was O(reps * N)).
    s1rep_to_members: dict[str, list[str]] = {}
    for r in stage1:
        for rep, members in r.clusters.items():
            s1rep_to_members[rep] = [rep, *members]

    final_clusters: dict[str, list[str]] = {}
    # Genomes a stage-1 pass rejected on QC belong to no cluster, so rebuilding
    # the status from cluster membership alone would drop them entirely and the
    # completeness check would then refuse the whole result. Seed the status with
    # them; the cluster-derived entries below take precedence on any overlap.
    status: dict[str, str] = {
        genome: state
        for r in stage1
        for genome, state in r.genome_status.items()
        if state == STATUS_FAIL_QC
    }

    for final_rep, s1reps_contained in stage2.clusters.items():
        status[final_rep] = STATUS_REPRESENTATIVE
        contained: list[str] = []
        # final_rep's own stage-1 cluster, plus every stage-1 rep stage 2 folded
        # into it. Each original genome belongs to exactly one stage-1 cluster, so
        # these groups are disjoint -- no double counting.
        for s1rep in (final_rep, *s1reps_contained):
            for genome in s1rep_to_members.get(s1rep, [s1rep]):
                if genome == final_rep:
                    continue
                contained.append(genome)
                status[genome] = STATUS_CONTAINED
        final_clusters[final_rep] = contained

    # A genome the later pass rejected on QC is in no stage-2 cluster either.
    # With nested levels, stage 2 is itself a composed result whose QC failures
    # come from its own chunks, so carry them here. A stage-1 representative
    # rejected there takes its stage-1 members with it: their cluster has no
    # representative left. QC is a per-genome check, so in practice this only
    # forwards failures from deeper levels; it keeps every genome accounted for.
    rejected_reps: list[str] = []
    for genome, state in stage2.genome_status.items():
        if state != STATUS_FAIL_QC or genome in status:
            continue
        members = s1rep_to_members.get(genome, [genome])
        if genome in s1rep_to_members:
            rejected_reps.append(genome)
        for member in members:
            status.setdefault(member, STATUS_FAIL_QC)
    if rejected_reps:
        n_members = sum(len(s1rep_to_members[r]) - 1 for r in rejected_reps)
        (logger or logging.getLogger(__name__)).warning(
            "%d representative(s) accepted by an earlier dereplication pass failed QC "
            "in a later pass (e.g. %s); they and their %d cluster member(s) are marked "
            "fail_qc. The QC result differs between calls for the same genome.",
            len(rejected_reps),
            ", ".join(sorted(rejected_reps)[:3]),
            n_members,
        )

    return DerepResult(
        representatives=stage2.representatives,
        clusters=final_clusters,
        genome_status=status,
    )


def _reduce_by_taxonomy(
    result: DerepResult,
    level: str,
    quality: Mapping[str, tuple[float, float]],
    logger,
    *,
    taxon_of: Mapping[str, str],
    n50: N50Of | None = None,
    prefer: Collection[str] = (),
    keep_apart: Collection[str] = (),
) -> DerepResult:
    """Collapse the ANI representatives to one per taxon (species|genus).

    ``taxon_of`` maps a genome filename to its taxon at ``level`` (the workdir
    stage reads it from the manifest, the merge step from selection.tsv or the
    filenames). Representatives sharing a taxon are merged into a single keeper.
    The representative of the largest cluster (then the lexically last name)
    is the default keeper. When ``quality`` is non-empty the quality keeper
    rule (:func:`~.derep_keeper.choose_keeper`, with ``n50`` for the N50
    term) re-picks among the group, exactly as within a cluster; pass ``{}``
    (keeper="tool", or no quality data) to use the largest-cluster rule
    alone. The others
    -- plus their cluster members -- become contained under the keeper.
    A representative named in ``prefer`` (the GTDB species representatives
    under ``--keeper gtdb``) wins over the others before any of these rules.
    Representatives whose taxon is unknown/empty are kept as-is (each its own
    group), so reduction never silently drops an un-annotated genome, and so
    are those named in ``keep_apart`` (the clusters of a refined lineage).
    """
    from ..dereplicators.base import STATUS_CONTAINED, STATUS_FAIL_QC, STATUS_REPRESENTATIVE

    # group rep filename -> taxon key; unknown taxon -> unique per-rep key (kept)
    groups: dict[str, list[str]] = {}
    rep_names = [rep.name for rep in result.representatives]
    for name in rep_names:
        taxon = "" if name in keep_apart else taxon_of.get(name, "")
        key = taxon if taxon else f"\x00{name}"  # sentinel: unknown taxa stay singletons
        groups.setdefault(key, []).append(name)

    rep_by_name = {rep.name: rep for rep in result.representatives}
    new_reps: list = []
    new_clusters: dict[str, list[str]] = {}
    # QC-rejected genomes are in no cluster; carry their status across the
    # reduction instead of losing them with the old membership.
    status: dict[str, str] = {
        genome: state for genome, state in result.genome_status.items() if state == STATUS_FAIL_QC
    }

    # keeper: the largest existing cluster, then lexical name; the quality
    # rule may then replace it, as within a cluster. Preferred genomes
    # (GTDB species representatives under --keeper gtdb), when the group
    # holds any, are the only candidates.
    picks: list[tuple[list[str], list[str], str]] = []
    for members in groups.values():
        pool = [n for n in members if n in prefer] or members
        by_size = max(pool, key=lambda n: (len(result.clusters.get(n, [])), n))
        picks.append((members, pool, by_size))
    if quality:
        prefetch_n50(n50, ((by_size, pool) for _, pool, by_size in picks), quality)

    left_to_size: list[str] = []
    for members, pool, keeper in picks:
        if quality:
            by_size = keeper
            keeper = choose_keeper(by_size, pool, quality, n50)
            if keeper == by_size and by_size not in quality and any(m in quality for m in pool):
                left_to_size.append(by_size)
        new_reps.append(rep_by_name[keeper])
        status[keeper] = STATUS_REPRESENTATIVE
        contained: list[str] = []
        for rep_name in members:
            # the keeper's own members, plus each folded rep and its members
            own = result.clusters.get(rep_name, [])
            folded = [] if rep_name == keeper else [rep_name]
            for g in (*folded, *own):
                contained.append(g)
                status[g] = STATUS_CONTAINED
        new_clusters[keeper] = contained

    log_unpromoted(left_to_size, logger, what="taxon group(s)")
    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
    )


def _taxon_lookup(ctx: WorkdirContext, level: str) -> dict[str, str]:
    """Map each genome filename to its taxon string at ``level`` from the manifest."""
    lookup: dict[str, str] = {}
    try:
        records = ctx.manifest.all_genomes(include_outgroup=True)
    except (sqlite3.OperationalError, OSError):
        # Manifest absent / not yet initialized (e.g. the data-channel path or
        # tests) -> empty taxonomy, everything stays a singleton. A corrupt
        # manifest raises sqlite3.DatabaseError, which is left to propagate.
        return lookup
    for rec in records:
        if not rec.filename:
            continue
        taxon = rec.species if level == "species" else rec.genus
        lookup[rec.filename] = taxon or ""
    return lookup


def quality_lookup(ctx: WorkdirContext) -> dict[str, tuple[float, float]]:
    """Map each genome filename to (completeness, contamination) from the manifest."""
    try:
        return ctx.manifest.quality()
    except (sqlite3.OperationalError, OSError):
        # Manifest absent / not yet initialized (e.g. the data-channel path or
        # tests) -> no quality data, the adapter's own representative stands.
        return {}


def _write_contract(
    ctx: WorkdirContext,
    result: DerepResult,
    quality: Mapping[str, tuple[float, float]],
    n50: N50Of | None = None,
    *,
    refined: Sequence[RefinedRow] = (),
) -> None:
    rep_dir = ctx.representatives_dir
    if rep_dir.exists():
        remove_tree(rep_dir)
    rep_dir.mkdir(parents=True, exist_ok=True)

    for rep in result.representatives:
        # Link from genomes/ when the genome is there: adapters may return a
        # copy in their scratch directory (skDER's output, dRep's staged
        # genomes), and linking that copy would hold a second full copy of
        # every representative on disk.
        staged = ctx.genomes_dir / rep.name
        source = staged if staged.exists() else rep
        if not source.exists():
            raise WorkdirError(f"Representative genome file missing: {rep.name}")
        if source.stat().st_size == 0:
            raise WorkdirError(
                f"Representative genome is empty: {rep.name} ({source}). An "
                "upstream download or staging step produced a zero-length genome."
            )
        link_or_copy(source, rep_dir / rep.name)

    write_clusters(ctx.derep_dir / CLUSTERS_TSV, result.clusters)
    write_genome_status(ctx.derep_dir / GENOME_STATUS_TSV, result.genome_status)
    refined_path = ctx.derep_dir / REFINED_TSV
    if refined:
        write_refined(refined_path, refined)
    else:
        # A result without a refinement leaves no file of an earlier one.
        refined_path.unlink(missing_ok=True)
    write_cluster_summary(
        ctx.derep_dir / CLUSTER_SUMMARY_TSV,
        summarise_clusters(
            result.clusters,
            quality,
            taxonomy_lookup(ctx),
            n50,
            gtdb_representatives=gtdb_lookup(ctx),
            refined={r.representative for r in refined if r.status == STATUS_REPRESENTATIVE},
        ),
    )


def _update_manifest(ctx: WorkdirContext, result: DerepResult) -> None:
    manifest = ctx.manifest
    rep_by_member: dict[str, str] = {}
    for rep, members in result.clusters.items():
        for m in members:
            rep_by_member[m] = rep
    updates: list[tuple[str, str, str | None]] = []
    for genome, status in result.genome_status.items():
        accession = _accession_from_filename(genome)
        if accession is None:
            continue
        representative = None
        if status == "contained":
            rep_file = rep_by_member.get(genome)
            representative = _accession_from_filename(rep_file) if rep_file else None
        updates.append((accession, status, representative))
    # One batched transaction instead of a commit per genome. Accessions absent
    # from the manifest are no-op UPDATEs (not errors), so we no longer swallow
    # exceptions -- a real database error should surface, not corrupt the manifest
    # silently.
    manifest.set_derep_status_many(updates)


def _accession_from_filename(filename: str | None) -> str | None:
    """Recover the accession from a canonical genome filename (None-safe)."""
    if not filename:
        return None
    return accession_from_filename(filename) or None
