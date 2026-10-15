"""SNP typing stage.

Selects a SNP typer, runs it against a reference (or reference-free),
optionally masks recombination with Gubbins, and writes the canonical SNP
outputs: ``snp/core_snp.fasta`` (+ optional VCF and SNP distance matrix, and
``snp/full_alignment.fasta`` when the typer produces a whole-genome alignment).
The core-SNP alignment is both a standalone typing deliverable and an MSA
source for the phylo stage; the whole-genome alignment is the required input
for recombination maskers.

The compute is factored into :func:`snptype_core`, a stateless engine that takes
explicit input/output directories and never touches the run config or manifest.
The workdir-bound :func:`run` resolves paths from the context, calls the core and
records provenance; the data-channel phylo step reuses the core directly.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTERS_TSV,
    CORE_SNP_FASTA,
    MASKER_TREE_NWK,
    atomic_path,
    list_fasta,
    record_name,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.integrity import check_genome_completeness, check_representatives_consistency
from ..core.plugins import scale_warning, warn_ignored_params, warn_unconsumed_extras
from ..core.process import remove_tree, stage_plain_inputs
from ..snptypers.base import SnpParams, SnpResult
from ..snptypers.base import registry as snp_registry
from ..tree.alignment import ALIGNMENT_STATS, constant_site_counts, write_alignment_stats


@dataclass
class SnptypeParams:
    tool: str = "simple"
    threads: int = 16
    reference: str | None = None
    all_genomes: bool = False
    mask: str = "none"  # none | gubbins
    allow_incomplete: bool = False
    # Record ids the masker leaves out of its scan (the outgroup) but keeps.
    mask_exclude: tuple[str, ...] = ()
    extra: dict = field(default_factory=dict)


def _check_inputs(ctx, all_genomes: bool, allow_incomplete: bool, logger) -> None:
    if all_genomes:
        check_genome_completeness(
            ctx.genomes_dir, ctx.workdir, logger=logger, allow_incomplete=allow_incomplete
        )
    else:
        check_representatives_consistency(
            ctx.representatives_dir,
            ctx.derep_dir / CLUSTERS_TSV,
            logger=logger,
            allow_incomplete=allow_incomplete,
        )


def snptype_core(
    genomes: list[Path],
    reference: Path | None,
    snp_dir: Path,
    scratch: Path,
    params: SnptypeParams,
    logger: logging.Logger,
    *,
    warn_extras: bool = True,
) -> tuple[SnpResult, dict[str, str]]:
    """Run a SNP typer over ``genomes`` into ``snp_dir`` (stateless; no config).

    ``reference`` is the already-resolved reference path (or None). For a
    reference-requiring typer a None reference falls back to ``genomes[0]``;
    reference-free typers ignore it. Returns the SNP result and tool versions.

    ``warn_extras`` reports extras this typer does not read. Callers that drive
    several adapters from one extras dict (the phylo stage) warn once themselves
    and pass False, so a key the tree builder reads is not flagged here.
    """
    if not genomes:
        raise WorkdirError("No genomes found for SNP typing. Run the genome (and derep) stages.")

    warn = scale_warning(snp_registry, params.tool, len(genomes))
    if warn:
        limit, alts = warn
        logger.warning(
            "SNP typer '%s' is tuned for <=%d genomes but you have %d; consider: %s",
            params.tool,
            limit,
            len(genomes),
            ", ".join(alts) or "none",
        )
    typer = snp_registry.create(params.tool)
    masker = None
    if params.mask not in ("none", ""):
        from ..maskers.base import registry as masker_registry

        masker = masker_registry.create(params.mask)
        if not getattr(typer, "produces_full_alignment", True):
            raise UserInputError(_mask_needs_full_alignment(params.mask, params.tool))
    if warn_extras:
        # The typer and the masker read one extras dict; a key is unread
        # only when neither declares it.
        _warn_unread_extras(typer.capabilities, masker, params.extra, logger)
    versions = typer.preflight()
    if masker is not None:
        # Before the typer runs: a missing masker binary must not surface
        # after minutes of SNP calling.
        versions.update(masker.preflight())
        # Likewise a set the masker would refuse by name, or an unusable
        # tuning value: both are known before any SNP is called.
        precheck = getattr(masker, "precheck", None)
        if precheck is not None:
            precheck(
                [record_name(g) for g in genomes if record_name(g) not in params.mask_exclude],
                dict(params.extra),
                logger,
            )

    ref = None
    if typer.requires_reference:
        # Both callers resolve the reference first (the medoid rule, a name
        # or the first genome, each logged there).
        ref = reference if reference is not None else genomes[0]

    snp_dir.mkdir(parents=True, exist_ok=True)
    if scratch.exists():
        remove_tree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    # A typer that cannot read gzip types decompressed copies; they keep the
    # record names, so the outputs name the genomes as the originals would.
    inputs_dir = scratch / "inputs"
    staged = stage_plain_inputs(
        [*genomes, *([ref] if ref is not None else [])], typer.capabilities, inputs_dir, logger
    )
    typed = [staged[g] for g in genomes]
    typed_ref = staged[ref] if ref is not None else None

    snp_params = SnpParams(
        threads=params.threads,
        reference=typed_ref,
        extra=dict(params.extra),
    )
    warn_ignored_params(typer.capabilities, snp_params, logger, family="SNP typer")
    logger.info("SNP typing %d genomes with %s", len(genomes), params.tool)
    try:
        result = typer.call(typed, typed_ref, scratch, snp_params, logger)
    finally:
        remove_tree(inputs_dir)

    core = snp_dir / CORE_SNP_FASTA
    # The sidecar describes the previous core; a failure between the two
    # writes must not leave it beside a new one.
    (snp_dir / ALIGNMENT_STATS).unlink(missing_ok=True)
    masked = False
    if masker is not None:
        from ..maskers.base import MaskParams

        if result.full_alignment is None:
            raise UserInputError(_mask_needs_full_alignment(params.mask, params.tool))
        filtered = masker.mask(
            result.full_alignment,
            scratch / params.mask,
            MaskParams(
                threads=params.threads,
                exclude=frozenset(params.mask_exclude),
                extra=dict(params.extra),
            ),
            logger,
        )
        with atomic_path(core) as tmp:
            shutil.copy2(filtered, tmp)
        masked = True
        masker_tree = masker.tree_output()
        if masker_tree is not None:
            with atomic_path(snp_dir / MASKER_TREE_NWK) as tmp:
                shutil.copy2(masker_tree, tmp)
            logger.info("The masker's own tree is kept as %s", snp_dir / MASKER_TREE_NWK)
        else:
            (snp_dir / MASKER_TREE_NWK).unlink(missing_ok=True)
    else:
        (snp_dir / MASKER_TREE_NWK).unlink(missing_ok=True)
        with atomic_path(core) as tmp:
            shutil.copy2(result.core_snp_fasta, tmp)

    # Optional outputs this typer does not produce are removed: left from an
    # earlier run with another typer or genome set, they would describe that run.
    for name, produced in (
        ("variants.vcf", result.vcf),
        ("snp_distance_matrix.tsv", result.snp_distance_matrix),
        ("full_alignment.fasta", result.full_alignment),
    ):
        if produced is None:
            (snp_dir / name).unlink(missing_ok=True)
    if result.vcf is not None:
        with atomic_path(snp_dir / "variants.vcf") as tmp:
            shutil.copy2(result.vcf, tmp)
    if result.snp_distance_matrix is not None:
        with atomic_path(snp_dir / "snp_distance_matrix.tsv") as tmp:
            shutil.copy2(result.snp_distance_matrix, tmp)
    full_alignment_out: Path | None = None
    constant_sites: tuple[int, int, int, int] | None = None
    if result.full_alignment is not None:
        full_alignment_out = snp_dir / "full_alignment.fasta"
        with atomic_path(full_alignment_out) as tmp:
            shutil.copy2(result.full_alignment, tmp)
        # Constant columns of the whole-genome alignment: the tree builder
        # adds them back as counts, which keeps branch lengths on scale.
        # Counted before masking: masking turns bases into N, which cannot
        # make a constant column variable; a column it makes constant is
        # missed, so after --mask the counts are a slight undercount.
        constant_sites, variable = constant_site_counts(full_alignment_out, logger)
        write_alignment_stats(
            core, constant_sites, variable, source=full_alignment_out.name, masked=masked
        )
        logger.info(
            "Constant sites in %s: A %d, C %d, G %d, T %d (%d variable columns)",
            full_alignment_out.name,
            *constant_sites,
            variable,
        )
    else:
        (snp_dir / ALIGNMENT_STATS).unlink(missing_ok=True)

    return (
        SnpResult(
            core_snp_fasta=core,
            vcf=(snp_dir / "variants.vcf") if result.vcf else None,
            snp_distance_matrix=(snp_dir / "snp_distance_matrix.tsv")
            if result.snp_distance_matrix
            else None,
            masked=masked,
            full_alignment=full_alignment_out,
            constant_sites=constant_sites,
        ),
        versions,
    )


def _mask_needs_full_alignment(mask: str, tool: str) -> str:
    """The refusal text for a typer that cannot feed a masker, naming the ones that can."""
    able = sorted(
        name
        for name in snp_registry.names()
        if not snp_registry.is_broken(name) and snp_registry.get(name).produces_full_alignment
    )
    return (
        f"--mask {mask} needs a whole-genome alignment, which SNP typer '{tool}' does not "
        f"produce. Use {', '.join(able) or 'a typer that writes one'}, or drop --mask."
    )


def _warn_unread_extras(typer_caps, masker, extra: dict, logger: logging.Logger) -> None:
    if masker is None:
        warn_unconsumed_extras(typer_caps, extra, logger, family="SNP typer")
        return
    known = set(typer_caps.accepted_extras) | set(typer_caps.default_params)
    known |= set(masker.capabilities.accepted_extras) | set(masker.capabilities.default_params)
    unread = sorted(set(extra) - known)
    if unread:
        logger.warning(
            "Neither SNP typer '%s' nor masker '%s' reads extra parameter(s): %s",
            typer_caps.name,
            masker.capabilities.name,
            ", ".join(unread),
        )


def run(ctx: WorkdirContext, params: SnptypeParams) -> SnpResult:
    logger = ctx.logger
    _check_inputs(ctx, params.all_genomes, params.allow_incomplete, logger)
    genomes = _genome_set(ctx, params.all_genomes)
    if not genomes:
        raise WorkdirError("No genomes found for SNP typing. Run the genome (and derep) stages.")

    reference = None
    if snp_registry.create(params.tool).requires_reference:
        reference = _reference_path(ctx, params.reference, genomes, logger, params.threads)

    result, versions = snptype_core(
        genomes, reference, ctx.snp_dir, ctx.scratch_dir / "snptype", params, logger
    )

    ctx.config.record_stage(
        "snptype",
        tool=params.tool,
        params={
            "all_genomes": params.all_genomes,
            "reference": reference.name if reference else None,
            "mask": params.mask,
        },
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("SNP typing complete: %s", result.core_snp_fasta)
    return result


def _genome_set(ctx: WorkdirContext, all_genomes: bool) -> list[Path]:
    source = ctx.genomes_dir if all_genomes else ctx.representatives_dir
    return list_fasta(source)


def _reference_path(ctx, reference_name, genomes, logger: logging.Logger, threads: int = 1) -> Path:
    from ..tree.medoid import REFERENCE_RULES

    if reference_name and reference_name not in REFERENCE_RULES:
        # Resolved against the workdir genome dirs only: a path component would
        # let the lookup escape them.
        if Path(reference_name).name != reference_name:
            raise UserInputError(
                f"--reference must be a genome file basename, not a path: {reference_name}"
            )
        # The genome set being typed first: under --all-genomes a representative
        # is in genomes/ too, and representatives/X beside genomes/X would be
        # two inputs with one record name. A genome of the set with the same
        # record name (X.fasta.gz for --reference X.fasta) is the reference too.
        for genome in genomes:
            if genome.name == reference_name:
                return genome
        wanted = record_name(reference_name)
        for genome in genomes:
            if record_name(genome) == wanted:
                logger.info("Using %s as the reference %s", genome.name, reference_name)
                return genome
        for base in (ctx.representatives_dir, ctx.genomes_dir):
            cand = base / reference_name
            if cand.exists():
                return cand
        raise UserInputError(f"Reference genome not found: {reference_name}")
    # No name: the medoid by sketch ANI when the workdir sketches and sourmash
    # are at hand, else the first genome with a warning; the record names it.
    from ..core.sketches import workdir_provider
    from ..tree.medoid import resolve_reference

    return resolve_reference(
        reference_name,
        list(genomes),
        None,
        logger,
        sketch_provider=workdir_provider(ctx, logger, threads, "snptype", write=False),
        scratch_dir=ctx.scratch_dir / "snptype",
        threads=threads,
    ).genome
