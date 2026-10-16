"""Finer dereplication inside a lineage of interest.

A dereplication applies one secondary ANI to the whole set. A study of a
species often needs coarse coverage of the species and more genomes inside
one lineage (a genotype, an outbreak clade). The lineage is named by seeds:
the user gives one or a few known genomes, every coarse cluster holding a
seed is the lineage, and its genomes are dereplicated again with the same
tool at a finer ANI. The fine clusters replace the coarse ones; the rest of
the result is unchanged.

``dereplicate --refine`` does this in one run (:mod:`.dereplicate` calls the
functions here); ``derep-refine`` (:func:`run`) applies it to the result on
disk without repeating the coarse pass.
"""

from __future__ import annotations

import csv
import logging
import sqlite3
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    REFINED_TSV,
    RefinedRow,
    accession_from_filename,
    list_fasta,
    read_clusters,
    read_genome_status,
    read_refined,
    strip_fasta_suffix,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.process import remove_tree
from ..core.sketches import adapter_sketches, workdir_provider
from ..dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
    check_result_complete,
    registry,
    run_quality,
)

# Recorded dereplicate params that are not tool arguments.
_RECORD_KEYS = frozenset(
    {
        "requested_tool",
        "primary_ani",
        "secondary_ani",
        "aligned_fraction",
        "process_size",
        "num_processes",
        "pre_primary_ani",
        "pre_secondary_ani",
        "reduce",
        "target_reps",
        "keeper",
        "keeper_effective",
        "keeper_swaps",
        "stock",
        "tool_args",
        "refine_seeds",
        "refine_ani",
        "refine_genomes",
        "refine_representatives",
        "refine_lineages",
        "refine_request",
        "refine_file",
        "base_completed",
    }
)


@dataclass
class DerepRefineParams:
    """``derep-refine``: refine the dereplication on disk."""

    seeds: list[str] = field(default_factory=list)
    ani: float | None = None
    # A TSV of seed and ani: several lineages, each with its own ANI.
    refine_file: str | None = None
    threads: int = 16


@dataclass(frozen=True)
class RefineGroup:
    """Seeds refined together at one ANI."""

    seeds: tuple[str, ...]
    ani: float


def read_refine_file(path: Path) -> list[RefineGroup]:
    """The groups of a refine file: a TSV with the columns ``seed`` and
    ``ani``, one seed per row; the seeds of one ANI form one group."""
    if not path.is_file():
        raise UserInputError(f"--refine-file {path} is not a file.")
    by_ani: dict[float, list[str]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo, delimiter="\t")
        if reader.fieldnames is None or not {"seed", "ani"} <= set(reader.fieldnames):
            raise UserInputError(
                f"--refine-file {path} needs the tab-separated columns seed and ani."
            )
        for n, rec in enumerate(reader, start=2):
            seed = (rec.get("seed") or "").strip()
            text = (rec.get("ani") or "").strip()
            if not seed or seed.startswith("#"):
                continue
            try:
                ani = float(text)
            except ValueError:
                raise UserInputError(
                    f"--refine-file {path}, line {n}: ani must be a number, not {text!r}."
                ) from None
            if not 0.0 < ani <= 1.0:
                raise UserInputError(
                    f"--refine-file {path}, line {n}: ani must be in (0, 1], got {ani}."
                )
            group = by_ani.setdefault(ani, [])
            if seed not in group:
                group.append(seed)
    if not by_ani:
        raise UserInputError(f"--refine-file {path} names no seed.")
    return [RefineGroup(tuple(seeds), ani) for ani, seeds in sorted(by_ani.items())]


def refine_groups(
    seeds: Sequence[str], ani: float | None, refine_file: str | None
) -> list[RefineGroup]:
    """The groups of a request: the ``--refine`` seeds at ``--refine-ani``,
    plus those of ``--refine-file``. A seed named twice is refused."""
    groups: list[RefineGroup] = []
    if seeds:
        if ani is None:
            raise UserInputError("--refine needs --refine-ani.")
        groups.append(RefineGroup(tuple(dict.fromkeys(s.strip() for s in seeds)), ani))
    if refine_file:
        groups += read_refine_file(Path(refine_file).expanduser())
    if not groups:
        raise UserInputError(
            "Give at least one --refine seed (with --refine-ani) or --refine-file."
        )
    seen: dict[str, float] = {}
    by_ani: dict[float, list[str]] = {}
    for g in groups:
        for seed in g.seeds:
            if seed in seen and seen[seed] != g.ani:
                raise UserInputError(
                    f"Seed {seed} is given at two ANIs ({seen[seed]} and {g.ani}); a lineage "
                    "has one."
                )
            seen[seed] = g.ani
            merged = by_ani.setdefault(g.ani, [])
            if seed not in merged:
                merged.append(seed)
    # The seeds of one ANI form one group, wherever they were given.
    return [RefineGroup(tuple(seeds), ani) for ani, seeds in sorted(by_ani.items())]


@dataclass
class Lineage:
    """The genomes a refinement works on, and where they came from."""

    genomes: list[str]  # filenames, the coarse representatives included
    coarse_representatives: list[str]
    seed_of: dict[str, str]  # genome -> the seed whose cluster holds it
    coarse_rep_of: dict[str, str]  # genome -> its coarse representative
    skipped: list[str]  # seeds alone in their cluster: nothing to refine
    ani: float = 0.0  # the secondary ANI of this lineage's fine pass


def check_refine_ani(ani: float, secondary_ani: float) -> None:
    """A refinement is finer than the run: its ANI must be above the run's."""
    if not 0.0 < ani <= 1.0:
        raise UserInputError(f"--refine-ani must be in (0, 1], got {ani}.")
    if ani <= secondary_ani:
        raise UserInputError(
            f"--refine-ani {ani} is not finer than the secondary ANI {secondary_ani}; "
            "a refinement keeps more genomes inside the lineage, so its ANI must be higher."
        )


def resolve_seeds(
    seeds: Iterable[str],
    genomes: Sequence[Path],
    *,
    accessions: Mapping[str, str] | None = None,
    outgroup: Collection[str] = (),
) -> list[str]:
    """The genome filename of each seed, which may be a filename, a filename
    without its FASTA suffix, or an accession (from the filename or from
    ``accessions``, accession -> filename). Unknown seeds and the outgroup
    are refused with their names.
    """
    by_key: dict[str, str] = {}
    by_accession: dict[str, set[str]] = {}
    for g in genomes:
        by_key[g.name] = g.name
        by_key[strip_fasta_suffix(g.name)] = g.name
        acc = accession_from_filename(g.name)
        if acc:
            by_accession.setdefault(acc, set()).add(g.name)
    for acc, filename in (accessions or {}).items():
        by_accession.setdefault(acc, set()).add(filename)
    out_keys = {n for o in outgroup for n in (o, strip_fasta_suffix(o))}
    present = {g.name for g in genomes}
    resolved: list[str] = []
    unknown: list[str] = []
    for seed in seeds:
        key = seed.strip()
        name = by_key.get(key)
        if name is None:
            # An accession: refused when two genomes under genomes/ carry it.
            files = sorted(by_accession.get(key, set()) & present)
            if len(files) > 1:
                raise UserInputError(
                    f"--refine {seed} names {len(files)} genomes ({', '.join(files)}); "
                    "give the filename."
                )
            name = files[0] if files else None
        if name is None:
            if key in out_keys:
                raise UserInputError(
                    f"--refine {seed}: the outgroup is not dereplicated, so it cannot seed "
                    "a lineage."
                )
            unknown.append(seed)
        elif name not in resolved:
            resolved.append(name)
    if unknown:
        raise UserInputError(
            "--refine names genome(s) not under genomes/ (give a filename, a filename "
            f"without its suffix, or an accession): {', '.join(unknown)}"
        )
    return resolved


def lineage_of(
    result: DerepResult,
    seeds: Sequence[str],
    logger: logging.Logger,
    prior: Sequence[RefinedRow] = (),
    ani: float = 0.0,
) -> Lineage:
    """The clusters of ``result`` to refine for ``seeds``: every cluster that
    holds a seed, whole, at ``ani``.

    ``prior`` is the ``refined.tsv`` of a result refined before. A seed that
    was refined then brings in its original coarse cluster (every genome the
    rows place under the same ``coarse_representative``), and every current
    cluster holding one of those genomes joins, so refining a seed again, at
    another ANI, starts from the same coarse cluster and not from the fine
    cluster the seed sits in now. ``coarse_rep_of`` keeps that original
    representative for such genomes.
    """
    rep_of: dict[str, str] = {}
    for rep, members in result.clusters.items():
        rep_of[rep] = rep
        for m in members:
            rep_of[m] = rep
    earlier_rep = {r.genome: r.coarse_representative for r in prior}
    earlier_cluster: dict[str, list[str]] = {}
    for genome, rep in earlier_rep.items():
        earlier_cluster.setdefault(rep, []).append(genome)
    genomes: list[str] = []
    coarse: list[str] = []
    seed_of: dict[str, str] = {}
    coarse_rep_of: dict[str, str] = {}
    skipped: list[str] = []
    for seed in seeds:
        if result.genome_status.get(seed) == STATUS_FAIL_QC:
            raise UserInputError(
                f"--refine {seed}: the genome failed the dereplicator's quality check and "
                "belongs to no cluster."
            )
        if seed not in rep_of:
            raise WorkdirError(f"--refine {seed}: the genome is in no cluster of the result.")
        # The seed's original coarse cluster when it was refined before, else
        # the seed; then every current cluster holding any of them.
        start = earlier_cluster.get(earlier_rep.get(seed, ""), [seed])
        reps = sorted({rep_of[g] for g in start if g in rep_of}, key=list(result.clusters).index)
        cluster = [g for rep in reps for g in (rep, *result.clusters[rep])]
        if len(cluster) == 1:
            skipped.append(seed)
            logger.info("Refine: %s is alone in its cluster; nothing to refine there", seed)
            continue
        for rep in reps:
            if rep in coarse:
                continue
            coarse.append(rep)
            for g in (rep, *result.clusters[rep]):
                genomes.append(g)
                seed_of[g] = seed
                coarse_rep_of[g] = earlier_rep.get(g, rep)
    return Lineage(genomes, coarse, seed_of, coarse_rep_of, skipped, ani)


def lineages_of(
    result: DerepResult,
    groups: Sequence[tuple[Sequence[str], float]],
    logger: logging.Logger,
    prior: Sequence[RefinedRow] = (),
) -> list[Lineage]:
    """One lineage per group of (resolved seeds, ani), all against the same
    ``result``, so each cluster belongs to at most one lineage; a cluster two
    groups claim is refused, naming the seeds."""
    lineages: list[Lineage] = []
    claimed: dict[str, str] = {}
    for seeds, ani in groups:
        lineage = lineage_of(result, seeds, logger, prior, ani)
        for rep in lineage.coarse_representatives:
            if rep in claimed:
                raise UserInputError(
                    f"Seeds {claimed[rep]} and {lineage.seed_of[rep]} fall in one cluster "
                    f"(of {rep}) but are given at different ANIs; a lineage has one."
                )
            claimed[rep] = lineage.seed_of[rep]
        if lineage.genomes:
            lineages.append(lineage)
    return lineages


def merge_rows(
    prior: Sequence[RefinedRow], rows: Sequence[RefinedRow], lineages: Sequence[Lineage]
) -> list[RefinedRow]:
    """``rows`` with the earlier refinements of genomes outside ``lineages`` kept."""
    inside = {g for lineage in lineages for g in lineage.genomes}
    return [r for r in prior if r.genome not in inside] + list(rows)


def fine_pass(
    adapter: Dereplicator,
    lineage: Lineage,
    paths: Mapping[str, Path],
    scratch: Path,
    derep_params: DerepParams,
    ani: float,
    logger: logging.Logger,
) -> DerepResult:
    """Dereplicate the lineage genomes with the run's settings at ``ani``."""
    genomes = [paths[name] for name in lineage.genomes]
    names = set(lineage.genomes)
    params = replace(
        derep_params,
        secondary_ani=ani,
        quality={n: q for n, q in derep_params.quality.items() if n in names},
        sketches=(
            None
            if derep_params.sketches is None
            else {p: s for p, s in derep_params.sketches.items() if p.name in names}
        ),
    )
    remove_tree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    logger.info(
        "Refine: dereplicating %d genome(s) of %d coarse cluster(s) at ANI %s",
        len(genomes),
        len(lineage.coarse_representatives),
        ani,
    )
    sub = adapter.dereplicate(genomes, scratch, params, logger)
    check_result_complete(sub, lineage.genomes)
    return sub


def splice(
    result: DerepResult,
    lineage: Lineage,
    sub: DerepResult,
    paths: Mapping[str, Path],
    logger: logging.Logger,
) -> DerepResult:
    """``result`` with the lineage's coarse clusters replaced by ``sub``'s."""
    coarse = set(lineage.coarse_representatives)
    clusters = {rep: list(m) for rep, m in result.clusters.items() if rep not in coarse}
    status = {g: s for g, s in result.genome_status.items() if g not in lineage.coarse_rep_of}
    rejected = sorted(g for g, s in sub.genome_status.items() if s == STATUS_FAIL_QC)
    if rejected:
        logger.warning(
            "Refine: %d genome(s) accepted by the coarse pass failed QC in the fine pass "
            "(e.g. %s); they are marked fail_qc. The QC result differs between calls for "
            "the same genome.",
            len(rejected),
            ", ".join(rejected[:3]),
        )
    for rep, members in sub.clusters.items():
        clusters[rep] = list(members)
    status.update(sub.genome_status)
    reps = [p for p in result.representatives if p.name not in coarse]
    reps += [paths[rep] for rep in sub.clusters]
    out = DerepResult(
        representatives=sorted(reps),
        clusters=clusters,
        genome_status=status,
        genome_information=result.genome_information,
    )
    check_result_complete(out, list(result.genome_status))
    logger.info(
        "Refine: %d coarse cluster(s) became %d cluster(s)",
        len(coarse),
        len(sub.clusters),
    )
    return out


def refined_rows(result: DerepResult, lineage: Lineage) -> list[RefinedRow]:
    """The ``refined.tsv`` rows for ``lineage`` as ``result`` finally places it."""
    rep_of: dict[str, str] = {}
    for rep, members in result.clusters.items():
        rep_of[rep] = rep
        for m in members:
            rep_of[m] = rep
    rows = []
    for g in lineage.genomes:
        status = result.genome_status.get(g, "")
        rows.append(
            RefinedRow(
                genome=g,
                seed=lineage.seed_of[g],
                coarse_representative=lineage.coarse_rep_of[g],
                representative="" if status == STATUS_FAIL_QC else rep_of.get(g, ""),
                status=status,
                ani=repr(lineage.ani),
            )
        )
    return rows


def record_params(lineages: Sequence[Lineage], rows: Sequence[RefinedRow]) -> dict:
    """The refinement's entries for a stage record: the totals, and one entry
    per lineage. ``refine_ani`` is the one ANI when every lineage shares it."""
    anis = sorted({lineage.ani for lineage in lineages})
    reps = {r.genome for r in rows if r.status == STATUS_REPRESENTATIVE}
    return {
        "refine_seeds": sorted({s for lineage in lineages for s in lineage.seed_of.values()}),
        "refine_ani": anis[0] if len(anis) == 1 else None,
        "refine_genomes": sum(len(lineage.genomes) for lineage in lineages),
        "refine_representatives": len(reps),
        "refine_lineages": [
            {
                "seeds": sorted(set(lineage.seed_of.values())),
                "ani": lineage.ani,
                "genomes": len(lineage.genomes),
                "representatives": sum(1 for g in lineage.genomes if g in reps),
            }
            for lineage in lineages
        ],
    }


def request_key(groups: Sequence[RefineGroup]) -> list[list[object]]:
    """The request as recorded, for the repeat check."""
    return [[sorted(g.seeds), g.ani] for g in sorted(groups, key=lambda g: g.ani)]


def manifest_accessions(ctx: WorkdirContext) -> tuple[dict[str, str], set[str]]:
    """(accession -> filename of the ingroup, names and accessions of the
    outgroup) from the manifest; empty without one."""
    try:
        records = ctx.manifest.all_genomes(include_outgroup=True)
    except (sqlite3.OperationalError, OSError):
        return {}, set()
    accessions = {
        r.accession: r.filename for r in records if r.accession and r.filename and not r.is_outgroup
    }
    outgroup = {n for r in records if r.is_outgroup for n in (r.filename, r.accession) if n}
    return accessions, outgroup


def run(ctx: WorkdirContext, params: DerepRefineParams) -> DerepResult:
    """Refine the dereplication on disk with the settings of its record."""
    from .derep_keeper import N50Lookup
    from .dereplicate import (
        _update_manifest,
        _write_contract,
        apply_keeper,
        quality_lookup,
    )

    logger = ctx.logger
    groups = refine_groups(params.seeds, params.ani, params.refine_file)
    record = ctx.config.stages.get("dereplicate")
    if record is None or not record.completed or not record.tool:
        raise UserInputError(
            "derep-refine refines a completed dereplication; none is recorded in this "
            "working directory. Run `repgenr dereplicate` first."
        )
    rec = record.params
    secondary = float(rec.get("secondary_ani", 0.99))
    for g in groups:
        check_refine_ani(g.ani, secondary)
    prior = ctx.config.stages.get("derep_refine")
    if prior is not None and _already_applied(ctx, groups, record.completed):
        # The result on disk is this refinement of this dereplication already;
        # refining the fine clusters again would give the same partition.
        logger.info(
            "derep-refine: the same seeds and ANI were applied to this dereplication "
            "already; nothing to do."
        )
        # The CLI keeps no record a stage did not write itself.
        ctx.config.record_stage(
            "derep_refine",
            tool=prior.tool,
            params=dict(prior.params),
            tool_versions=dict(prior.tool_versions),
            completed=datetime.now(UTC).isoformat(),
        )
        ctx.save_config()
        clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
        return DerepResult(
            representatives=[ctx.genomes_dir / rep for rep in clusters],
            clusters=clusters,
            genome_status=read_genome_status(ctx.derep_dir / GENOME_STATUS_TSV),
        )
    clusters_path = ctx.derep_dir / CLUSTERS_TSV
    if not clusters_path.is_file():
        raise WorkdirError(f"Missing {clusters_path}. Run the dereplicate stage first.")
    clusters = read_clusters(clusters_path)
    status = read_genome_status(ctx.derep_dir / GENOME_STATUS_TSV)
    for rep, members in clusters.items():
        status.setdefault(rep, STATUS_REPRESENTATIVE)
        for m in members:
            status.setdefault(m, STATUS_CONTAINED)
    genomes = list_fasta(ctx.genomes_dir)
    paths = {g.name: g for g in genomes}
    missing = sorted(n for n in status if n not in paths)
    if missing:
        raise WorkdirError(
            f"{len(missing)} genome(s) of the dereplication are not under {ctx.genomes_dir} "
            f"(e.g. {missing[0]}); the result cannot be refined."
        )
    result = DerepResult(
        representatives=[paths[rep] for rep in clusters],
        clusters=clusters,
        genome_status=status,
    )
    accessions, outgroup = manifest_accessions(ctx)
    resolved = [
        (resolve_seeds(g.seeds, genomes, accessions=accessions, outgroup=outgroup), g.ani)
        for g in groups
    ]
    prior_rows = read_refined(ctx.derep_dir / REFINED_TSV)
    lineages = lineages_of(result, resolved, logger, prior_rows)
    if not lineages:
        raise UserInputError(
            "Every --refine seed is alone in its cluster; there is nothing to refine."
        )
    if rec.get("reduce", "none") != "none":
        logger.warning(
            "The dereplication was reduced to one representative per %s, so a cluster on "
            "disk is a whole taxon group; the lineage is refined as such.",
            rec.get("reduce"),
        )

    tool = record.tool
    if tool not in registry.names():
        raise UserInputError(
            f"The dereplication was made with '{tool}', which is not an available "
            f"dereplicator here; choose from {', '.join(registry.names())}."
        )
    adapter = registry.create(tool)
    versions = adapter.preflight()
    caps = adapter.capabilities
    extra = {**caps.default_params, **_tool_args(rec)}
    manifest_quality = quality_lookup(ctx)
    lineage_names = [n for lineage in lineages for n in lineage.genomes]
    lineage_paths = [paths[n] for n in lineage_names]
    derep_params = DerepParams(
        primary_ani=float(rec.get("primary_ani", 0.90)),
        secondary_ani=secondary,
        aligned_fraction=float(rec.get("aligned_fraction", 0.50)),
        threads=params.threads,
        extra=extra,
        quality=run_quality(manifest_quality, lineage_names, logger, source="manifest"),
    )
    derep_params.sketches = adapter_sketches(
        adapter,
        extra,
        lineage_paths,
        workdir_provider(ctx, logger, params.threads, "derep-refine"),
        logger,
        "derep-refine",
    )
    scratch = ctx.scratch_dir / "derep_refine"
    n50 = N50Lookup([ctx.genomes_dir], workers=params.threads)
    keeper = str(rec.get("keeper", "quality"))
    refined = result
    for i, lineage in enumerate(lineages):
        sub = fine_pass(
            adapter, lineage, paths, scratch / str(i), derep_params, lineage.ani, logger
        )
        sub, _swaps, effective, _quality, _gtdb = apply_keeper(
            ctx, sub, keeper, manifest_quality, n50, logger, gtdb_warning=False
        )
        if keeper == "gtdb" and effective != "gtdb":
            logger.info(
                "Refine: no GTDB species representative among the lineage genomes; the %s "
                "rule chose the representatives",
                "quality" if effective == "quality" else "tool's own",
            )
        refined = splice(refined, lineage, sub, paths, logger)
    rows = merge_rows(
        prior_rows,
        [r for lineage in lineages for r in refined_rows(refined, lineage)],
        lineages,
    )
    _write_contract(ctx, refined, manifest_quality, n50, refined=rows)
    _update_manifest(ctx, refined)
    remove_tree(scratch)
    ctx.config.record_stage(
        "derep_refine",
        tool=tool,
        params={
            **record_params(lineages, rows),
            # The request as given and the dereplication it refined, so a
            # repeat of the same request on the same result does nothing.
            "refine_request": request_key(groups),
            "base_completed": record.completed,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info(
        "Refinement complete: %d representatives of %d genomes",
        len(refined.representatives),
        len(refined.genome_status),
    )
    return refined


def _already_applied(ctx: WorkdirContext, groups: Sequence[RefineGroup], base: str | None) -> bool:
    """Whether the result on disk is this request applied to the dereplication
    completed at ``base``.

    The record of the last run (its completion stamp is cleared before this
    run starts, its params are kept) and the files on disk settle it: the
    same seeds, ANI and base, and a ``refined.tsv`` whose rows agree with
    ``clusters.tsv``. Reached when the CLI reruns the stage because the
    genome set or the tool environment changed, or under ``--force``.
    """
    prior = ctx.config.stages.get("derep_refine")
    if prior is None or base is None:
        return False
    p = prior.params
    if p.get("base_completed") != base or p.get("refine_request") != request_key(groups):
        return False
    rows = read_refined(ctx.derep_dir / REFINED_TSV)
    if not rows:
        return False
    clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
    for r in rows:
        if r.status == STATUS_REPRESENTATIVE and r.representative != r.genome:
            return False
        if r.status == STATUS_REPRESENTATIVE and r.genome not in clusters:
            return False
        if r.status == STATUS_CONTAINED and r.genome not in clusters.get(r.representative, []):
            return False
    return True


def _tool_args(rec: Mapping[str, object]) -> dict:
    """The tool arguments of a recorded dereplication: ``tool_args`` when the
    record has it, else the keys that are not the stage's own."""
    args = rec.get("tool_args")
    if isinstance(args, dict):
        return dict(args)
    return {k: v for k, v in rec.items() if k not in _RECORD_KEYS}
