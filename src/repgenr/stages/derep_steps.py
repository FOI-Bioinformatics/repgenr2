"""Discrete dereplication steps for scatter-gather orchestration.

The shared-workdir ``dereplicate`` stage runs the whole two-stage pipeline in one
process. For horizontal scaling (Nextflow scatter-gather across nodes) the same
work is exposed as two stateless, file-in / file-out steps:

* :func:`dereplicate_chunk` dereplicates one chunk of genomes and writes a chunk
  result directory in the canonical contract (``representatives/`` +
  ``clusters.tsv`` + ``genome_status.tsv``).
* :func:`dereplicate_merge` takes several chunk result directories, dereplicates
  the union of their representatives with the final thresholds, and composes the
  two-stage membership into a final contract directory.

Neither step touches the SQLite manifest or the run config: their entire state is
the directories they read and write, so Nextflow can stage them between tasks and
cache them independently. A chunk result directory is itself a valid contract, so
``dereplicate_merge`` consumes exactly what ``dereplicate_chunk`` produces. This
is the first level and the final pass of the in-process ``_dereplicate_chunked``,
split into steps the orchestrator scatters. The ``DEREPLICATE_SCATTER``
subworkflow runs one chunk round and one merge: the merge dereplicates the whole
union of chunk representatives in one call, whatever its size, and does not
chunk it again.
"""

from __future__ import annotations

import logging
import os
import shutil
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    atomic_replace,
    parse_genome_filename,
    read_clusters,
    read_genome_status,
    read_selection,
    write_cluster_summary,
    write_clusters,
    write_genome_status,
)
from ..core.errors import WorkdirError
from ..core.plugins import warn_ignored_params, warn_unconsumed_extras
from ..core.process import link_or_copy, remove_tree
from ..core.sketches import adapter_sketches, directory_provider
from ..core.versions import write_versions_fragment
from ..dereplicators.base import (
    DerepParams,
    DerepResult,
    check_result_complete,
    registry,
    run_quality,
)
from .cluster_summary import Taxonomy, summarise_clusters
from .derep_keeper import (
    N50Lookup,
    N50Of,
    prefer_gtdb_representatives,
    rescore_representatives,
)
from .dereplicate import (
    DereplicateParams,
    _compose_two_stage,
    _reduce_by_taxonomy,
    _search_target_reps,
)

_REPRESENTATIVES_DIR = "representatives"
# N50 of the chunk's scored genomes, so the merge step scores members whose
# files it does not receive (only chunk representatives are staged there).
GENOME_N50_TSV = "genome_n50.tsv"


def _write_n50(path: Path, values: Mapping[str, int]) -> None:
    with atomic_replace(path, newline="") as fo:
        fo.write("genome\tn50\n")
        for name in sorted(values):
            fo.write(f"{name}\t{values[name]}\n")


def _read_n50(chunk_dirs: list[Path]) -> dict[str, int]:
    """N50 values the chunks recorded; a chunk without the file adds none."""
    values: dict[str, int] = {}
    for d in chunk_dirs:
        path = d / GENOME_N50_TSV
        if not path.is_file():
            continue
        with open(path, encoding="utf-8") as fo:
            next(fo, None)
            for line in fo:
                name, _, value = line.rstrip("\n").partition("\t")
                if name and value.isdigit():
                    values[name] = int(value)
    return values


def _maybe_write_versions(path: Path | None, versions: dict[str, str]) -> None:
    if path is not None:
        write_versions_fragment(path, versions)


def _gtdb_from_selection(selection_tsv: Path | None) -> set[str]:
    """Filenames selection.tsv flags as GTDB species representatives.

    The outgroup row is left out: it is usually a GTDB representative itself
    but is never dereplicated.
    """
    if selection_tsv is None:
        return set()
    return {
        r.filename
        for r in read_selection(selection_tsv)
        if r.gtdb_representative and not r.is_outgroup
    }


def _warn_no_gtdb_flags(
    flagged: set[str], selection_tsv: Path, step: str, logger: logging.Logger
) -> None:
    if not flagged:
        logger.warning(
            "%s: %s flags no GTDB species representative, so --keeper gtdb falls "
            "back to the quality rule. Only GTDB selections from the metadata "
            "stage carry the gtdb_representative column.",
            step,
            selection_tsv,
        )


def _quality_from_selection(selection_tsv: Path) -> dict[str, tuple[float, float]]:
    """filename -> (completeness, contamination) for rows carrying both values."""
    return {
        r.filename: (r.completeness, r.contamination)
        for r in read_selection(selection_tsv)
        if r.completeness is not None and r.contamination is not None
    }


@dataclass
class ChunkParams:
    tool: str
    genomes: list[Path]
    out_dir: Path
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    extra: dict | None = None
    versions_out: Path | None = None
    # A selection.tsv (with quality columns) that enables the quality-aware
    # keeper below. Unlike at the merge step, every genome in the chunk has a
    # real file here (params.genomes), so any promotion is always resolvable.
    selection_tsv: Path | None = None
    keeper: str = "quality"  # quality | gtdb | tool
    # A sketches/ directory (Nextflow SKETCH) for an adapter that compares
    # sourmash sketches; matched to the genomes by record name.
    sketches_dir: Path | None = None


@dataclass
class MergeParams:
    tool: str
    chunk_dirs: list[Path]
    out_dir: Path
    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    extra: dict | None = None
    versions_out: Path | None = None
    # A selection.tsv (with quality columns) that enables the quality-aware
    # keeper below; without it (or with keeper="tool") the adapter's own
    # merge-level pick stands, as before.
    selection_tsv: Path | None = None
    keeper: str = "quality"  # quality | gtdb | tool
    # Taxonomy-aware reduction after the merge (none | species | genus). The
    # taxonomy comes from selection.tsv when given, else from the canonical
    # genome filenames.
    reduce: str = "none"
    # Search the secondary ANI of the merge pass to land near this many
    # representatives (0 = off).
    target_reps: int = 0
    sketches_dir: Path | None = None


def dereplicate_chunk(params: ChunkParams, logger: logging.Logger) -> DerepResult:
    """Dereplicate one chunk of genomes; write its contract to ``out_dir``."""
    if not params.genomes:
        raise WorkdirError("dereplicate-chunk: no genome paths provided.")
    missing = [g for g in params.genomes if not g.exists()]
    if missing:
        raise WorkdirError(
            f"dereplicate-chunk: {len(missing)} genome file(s) not found, e.g. {missing[0]}"
        )

    adapter = registry.create(params.tool)
    caps = adapter.capabilities
    warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")
    versions = adapter.preflight()
    _maybe_write_versions(params.versions_out, versions)
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra={**caps.default_params, **(params.extra or {})},
        # As in the stage: selection.tsv quality is an adapter input (dRep
        # --genomeInfo, galah --genome-info) whichever keeper rule applies,
        # decided over the whole selection so every chunk agrees.
        quality=_adapter_quality(params.selection_tsv, logger),
    )
    warn_ignored_params(caps, derep_params, logger, family="Dereplicator")
    derep_params.sketches = _step_sketches(
        adapter,
        derep_params.extra,
        params.genomes,
        params.sketches_dir,
        logger,
        "dereplicate-chunk",
    )
    scratch = _fresh(params.out_dir / "scratch")
    result = adapter.dereplicate(params.genomes, scratch, derep_params, logger)
    n50 = N50Lookup(sorted({g.parent for g in params.genomes}))

    if params.keeper == "quality" and params.selection_tsv is not None:
        quality = _quality_from_selection(params.selection_tsv)
        result, keeper_swaps = rescore_representatives(result, quality, logger, n50)
        if keeper_swaps:
            logger.info(
                "dereplicate-chunk: quality-aware keeper changed %d representative(s)",
                keeper_swaps,
            )
    elif params.keeper == "gtdb" and params.selection_tsv is not None:
        flagged = _gtdb_from_selection(params.selection_tsv)
        _warn_no_gtdb_flags(flagged, params.selection_tsv, "dereplicate-chunk", logger)
        result, keeper_swaps = prefer_gtdb_representatives(
            result,
            flagged,
            _quality_from_selection(params.selection_tsv),
            logger,
            n50,
        )
        logger.info("dereplicate-chunk: GTDB keeper changed %d representative(s)", keeper_swaps)
    elif params.keeper == "tool" and params.selection_tsv is not None:
        logger.info(
            "dereplicate-chunk: --selection-tsv given but --keeper is 'tool'; "
            "keeping the adapter's own representatives."
        )

    check_result_complete(result, [g.name for g in params.genomes])
    fallbacks = sorted({g.parent for g in params.genomes})
    summary_quality = _summary_quality(params.selection_tsv)
    _write_step_contract(
        params.out_dir,
        result,
        fallbacks,
        summary_quality,
        _summary_taxonomy(params.selection_tsv),
        n50,
        _gtdb_from_selection(params.selection_tsv),
    )
    # The merge step receives only the representatives' files, so record the
    # N50 of the scored members it may still compare. In a cluster with two or
    # more scored genomes the keeper has read them already (cached); only a
    # scored member under an unscored representative is read here.
    for members in result.clusters.values():
        for name in members:
            if name in summary_quality:
                n50(name)
    _write_n50(params.out_dir / GENOME_N50_TSV, n50.computed())
    shutil.rmtree(scratch, ignore_errors=True)  # drop tool intermediates from the output
    logger.info(
        "dereplicate-chunk: %d genomes -> %d representatives (%s)",
        len(params.genomes),
        len(result.representatives),
        params.tool,
    )
    return result


def dereplicate_merge(params: MergeParams, logger: logging.Logger) -> DerepResult:
    """Dereplicate the union of chunk representatives; compose the final contract."""
    if not params.chunk_dirs:
        raise WorkdirError("dereplicate-merge: no chunk directories provided.")
    stage1 = [_load_chunk(d) for d in params.chunk_dirs]
    union = [rep for r in stage1 for rep in r.representatives]
    if not union:
        raise WorkdirError("dereplicate-merge: the chunk directories hold no representatives.")

    adapter = registry.create(params.tool)
    caps = adapter.capabilities
    warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")
    versions = adapter.preflight()
    _maybe_write_versions(params.versions_out, versions)
    union_names = {rep.name for rep in union}
    derep_params = DerepParams(
        primary_ani=params.primary_ani,
        secondary_ani=params.secondary_ani,
        aligned_fraction=params.aligned_fraction,
        threads=params.threads,
        extra={**caps.default_params, **(params.extra or {})},
        quality={
            name: qual
            for name, qual in _adapter_quality(params.selection_tsv, logger).items()
            if name in union_names
        },
    )
    warn_ignored_params(caps, derep_params, logger, family="Dereplicator")
    derep_params.sketches = _step_sketches(
        adapter, derep_params.extra, union, params.sketches_dir, logger, "dereplicate-merge"
    )
    scratch = _fresh(params.out_dir / "scratch")
    if params.target_reps > 0:
        # The search re-runs the merge pass per step; the stateless step has
        # no chunking of its own, so the single-pass settings apply.
        stage2 = _search_target_reps(
            adapter,
            union,
            scratch,
            derep_params,
            DereplicateParams(tool=params.tool),
            params.target_reps,
            logger,
        )
    else:
        stage2 = adapter.dereplicate(union, scratch, derep_params, logger)
    final = _compose_two_stage(stage1, stage2, logger)
    rep_dirs = [d / _REPRESENTATIVES_DIR for d in params.chunk_dirs]
    n50 = N50Lookup(rep_dirs, known=_read_n50(params.chunk_dirs))

    if params.keeper == "quality" and params.selection_tsv is not None:
        # Only stage-1 representatives are staged at the merge step (each
        # chunk's representatives/ directory); a chunk-level CONTAINED member's
        # file never reaches this step, even though _compose_two_stage's
        # expanded membership lists it in final.clusters. Restrict the quality
        # map to resolvable names so such a member can never be promoted to a
        # representative _write_step_contract cannot then find a file for.
        resolvable = {rep.name for r in stage1 for rep in r.representatives}
        quality = {
            name: qual
            for name, qual in _quality_from_selection(params.selection_tsv).items()
            if name in resolvable
        }
        final, keeper_swaps = rescore_representatives(final, quality, logger, n50)
        if keeper_swaps:
            logger.info(
                "dereplicate-merge: quality-aware keeper changed %d representative(s)",
                keeper_swaps,
            )
    elif params.keeper == "gtdb" and params.selection_tsv is not None:
        # As for the quality rule: only stage-1 representatives have a file
        # at this step, so only they may be promoted.
        resolvable = {rep.name for r in stage1 for rep in r.representatives}
        flagged = _gtdb_from_selection(params.selection_tsv)
        _warn_no_gtdb_flags(flagged, params.selection_tsv, "dereplicate-merge", logger)
        final, keeper_swaps = prefer_gtdb_representatives(
            final,
            flagged & resolvable,
            {
                name: qual
                for name, qual in _quality_from_selection(params.selection_tsv).items()
                if name in resolvable
            },
            logger,
            n50,
        )
        logger.info("dereplicate-merge: GTDB keeper changed %d representative(s)", keeper_swaps)
    elif params.keeper == "tool" and params.selection_tsv is not None:
        logger.info(
            "dereplicate-merge: --selection-tsv given but --keeper is 'tool'; "
            "keeping the adapter's own representatives."
        )

    if params.reduce != "none":
        resolvable = {rep.name for r in stage1 for rep in r.representatives}
        # As in the dereplicate stage: --keeper tool ranks the representatives
        # of a taxon by cluster size alone, without the selection.tsv quality.
        quality = (
            {}
            if params.selection_tsv is None or params.keeper == "tool"
            else {
                name: qual
                for name, qual in _quality_from_selection(params.selection_tsv).items()
                if name in resolvable
            }
        )
        before = len(final.representatives)
        final = _reduce_by_taxonomy(
            final,
            params.reduce,
            quality,
            logger,
            taxon_of=_taxon_from_selection_or_names(
                params.selection_tsv, resolvable, params.reduce
            ),
            n50=n50,
            prefer=(
                _gtdb_from_selection(params.selection_tsv) & resolvable
                if params.keeper == "gtdb"
                else set()
            ),
        )
        logger.info(
            "dereplicate-merge: --reduce %s collapsed %d representatives to %d",
            params.reduce,
            before,
            len(final.representatives),
        )

    # Every genome the chunks saw: cluster members plus the genomes that carry a
    # status without a cluster (QC rejects), so the completeness check covers the
    # whole input set rather than only the clustered part of it.
    all_names = {
        name for r in stage1 for rep, members in r.clusters.items() for name in (rep, *members)
    }
    all_names |= {genome for r in stage1 for genome in r.genome_status}
    check_result_complete(final, all_names)

    # The final representatives are stage-2 representative paths, which live in the
    # chunk representatives/ directories; fall back to those when resolving files.
    fallbacks = rep_dirs
    _write_step_contract(
        params.out_dir,
        final,
        fallbacks,
        _summary_quality(params.selection_tsv),
        _summary_taxonomy(params.selection_tsv),
        n50,
        _gtdb_from_selection(params.selection_tsv),
    )
    _write_n50(params.out_dir / GENOME_N50_TSV, n50.computed())
    shutil.rmtree(scratch, ignore_errors=True)  # drop tool intermediates from the output
    logger.info(
        "dereplicate-merge: %d chunks, union of %d reps -> %d representatives (%s)",
        len(params.chunk_dirs),
        len(union),
        len(final.representatives),
        params.tool,
    )
    return final


def _taxon_from_selection_or_names(
    selection_tsv: Path | None, names: set[str], level: str
) -> dict[str, str]:
    """Genome filename -> taxon at ``level``: from selection.tsv rows when given,
    else parsed from the canonical filename (Family_Genus_species_Accession)."""
    taxon_of: dict[str, str] = {}
    if selection_tsv is not None:
        for row in read_selection(selection_tsv):
            taxon_of[row.filename] = row.species if level == "species" else row.genus
    for name in names:
        if name in taxon_of:
            continue
        _family, genus, species, _acc = parse_genome_filename(name)
        taxon_of[name] = species if level == "species" else genus
    return taxon_of


def _step_sketches(
    adapter: object,
    extra: dict,
    genomes: list[Path],
    sketches_dir: Path | None,
    logger: logging.Logger,
    step: str,
) -> dict[Path, Path] | None:
    """Sketches from a staged sketches/ directory, for an adapter that reads them."""
    if sketches_dir is None:
        return None
    if not sketches_dir.is_dir():
        raise WorkdirError(f"{step}: --sketches-dir {sketches_dir} is not a directory.")
    return adapter_sketches(
        adapter, extra, genomes, directory_provider(sketches_dir, logger, step), logger, step
    )


def _load_chunk(chunk_dir: Path) -> DerepResult:
    """Read a chunk result directory back into a DerepResult.

    Reads the per-genome statuses as well as the clusters: genomes the chunk
    rejected on QC appear only in ``genome_status.tsv``, and dropping them here
    would lose them from the merged contract.
    """
    clusters_path = chunk_dir / CLUSTERS_TSV
    if not clusters_path.exists():
        raise WorkdirError(f"dereplicate-merge: {clusters_path} not found (not a chunk result).")
    clusters = read_clusters(clusters_path)
    status_path = chunk_dir / GENOME_STATUS_TSV
    if not status_path.exists():
        raise WorkdirError(f"dereplicate-merge: {status_path} not found (not a chunk result).")
    status = read_genome_status(status_path)
    rep_dir = chunk_dir / _REPRESENTATIVES_DIR
    reps: list[Path] = []
    for rep_name in clusters:
        rep_path = rep_dir / rep_name
        if not rep_path.exists():
            raise WorkdirError(
                f"dereplicate-merge: representative file missing in chunk: {rep_path}"
            )
        reps.append(rep_path)
    return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


def _write_step_contract(
    out_dir: Path,
    result: DerepResult,
    fallback_dirs: list[Path],
    quality: Mapping[str, tuple[float, float]],
    taxonomy: Taxonomy | None = None,
    n50: N50Of | None = None,
    gtdb_representatives: set[str] | None = None,
) -> None:
    """Write representatives/, clusters.tsv, genome_status.tsv and
    cluster_summary.tsv under ``out_dir``."""
    rep_dir = out_dir / _REPRESENTATIVES_DIR
    if rep_dir.exists():
        remove_tree(rep_dir)
    rep_dir.mkdir(parents=True, exist_ok=True)
    for rep in result.representatives:
        source = rep if rep.exists() else _find(fallback_dirs, rep.name)
        if source is None:
            raise WorkdirError(f"Representative genome file missing: {rep.name}")
        if os.path.getsize(os.path.realpath(source)) == 0:
            raise WorkdirError(
                f"Representative genome is empty: {rep.name}. The source file "
                f"({source}) has zero length -- an upstream download or staging "
                "step produced an empty genome."
            )
        link_or_copy(source, rep_dir / rep.name)
    write_clusters(out_dir / CLUSTERS_TSV, result.clusters)
    write_genome_status(out_dir / GENOME_STATUS_TSV, result.genome_status)
    write_cluster_summary(
        out_dir / CLUSTER_SUMMARY_TSV,
        summarise_clusters(
            result.clusters,
            quality,
            taxonomy,
            n50,
            gtdb_representatives=gtdb_representatives or set(),
        ),
    )


def _adapter_quality(
    selection_tsv: Path | None, logger: logging.Logger
) -> dict[str, tuple[float, float]]:
    """selection.tsv quality for the adapters: all of it when every selected
    (non-outgroup) genome has values, else none. Decided over the whole
    selection, not the chunk, so all chunks and the merge of a run agree."""
    if selection_tsv is None:
        return {}
    names = [r.filename for r in read_selection(selection_tsv) if not r.is_outgroup and r.filename]
    return run_quality(
        _quality_from_selection(selection_tsv), names, logger, source="selection.tsv"
    )


def _summary_quality(selection_tsv: Path | None) -> dict[str, tuple[float, float]]:
    """Quality for the cluster summary: from selection.tsv when given, whatever
    the keeper rule, so the summary can show a member outscoring the keeper."""
    return {} if selection_tsv is None else _quality_from_selection(selection_tsv)


def _summary_taxonomy(selection_tsv: Path | None) -> dict[str, tuple[str, str]]:
    """Taxonomy for the cluster summary from selection.tsv; without one the
    summary parses the canonical filenames."""
    if selection_tsv is None:
        return {}
    return {r.filename: (r.genus, r.species) for r in read_selection(selection_tsv) if r.filename}


def _find(dirs: list[Path], name: str) -> Path | None:
    for d in dirs:
        candidate = d / name
        if candidate.exists():
            return candidate
    return None


def _fresh(path: Path) -> Path:
    if path.exists():
        remove_tree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
