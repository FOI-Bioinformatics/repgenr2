"""Single construction point for each pipeline stage's parameters.

Both the manual per-stage commands and the ``run`` orchestrator build stage
parameters through these functions, so the two entry points cannot drift:
identical CLI intent produces identical params objects and therefore identical
resume fingerprints. Defaults live on the params dataclasses -- a builder
argument left unset is simply not passed, so a default changed on the
dataclass takes effect everywhere at once. Shared validation
(tool/threshold checks) also lives here so both entry points reject the same
inputs the same way.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING, Any

from .base import _require_choice, _require_unit_interval

if TYPE_CHECKING:
    from ..stages.assemble import AssembleParams
    from ..stages.dereplicate import DereplicateParams
    from ..stages.genome import GenomeParams
    from ..stages.ingest import IngestParams
    from ..stages.metadata import MetadataParams
    from ..stages.phylo import PhyloParams
    from ..stages.reads import ReadsParams
    from ..stages.tree2tax import Tree2taxParams
    from ..stages.vgenome import VgenomeParams
    from ..stages.vmetadata import VmetadataParams

# Sentinel for "caller did not set this": the argument is not passed to the
# dataclass, whose own default then applies.
_UNSET: Any = object()

# Closed choice sets shared with `run`, which validates the same values under
# its own flag names (--metadata-source, --viral-source).
METADATA_DATASETS = frozenset({"all", "rep"})
METADATA_LEVELS = frozenset({"family", "genus", "species"})
METADATA_SOURCES = frozenset({"tsv", "api"})
VIRAL_SOURCES = frozenset({"ncbi_virus", "bvbrc"})
LENGTH_METHODS = frozenset({"median_of_medians", "mean"})
DEREP_STOCK_ACTIONS = frozenset({"list", "pack", "unpack", "delete"})


def _build(cls: type, **kwargs: Any):
    return cls(**{k: v for k, v in kwargs.items() if v is not _UNSET})


def require_mask(mask: str) -> None:
    """``--mask`` is ``none`` or a registered masker."""
    from ..maskers.base import registry as _mask_registry

    _require_choice(mask, {"none", *_mask_registry.names()}, "--mask")


def metadata_params(
    *,
    dataset: str,
    level: str,
    source: Any = _UNSET,
    release: Any = _UNSET,
    version: Any = _UNSET,
    target_family: Any = _UNSET,
    target_genus: Any = _UNSET,
    target_species: Any = _UNSET,
    outgroup_accession: Any = _UNSET,
    metadata_path: Any = _UNSET,
    nodownload: Any = _UNSET,
    limit: Any = _UNSET,
    drop_foreign: Any = _UNSET,
) -> MetadataParams:
    from ..stages.metadata import MetadataParams

    _require_choice(dataset, METADATA_DATASETS, "--dataset")
    _require_choice(level, METADATA_LEVELS, "--level")
    if source is not _UNSET:
        _require_choice(source, METADATA_SOURCES, "--source")
    return _build(
        MetadataParams,
        dataset=dataset,
        level=level,
        source=source,
        release=release,
        version=version,
        target_family=target_family,
        target_genus=target_genus,
        target_species=target_species,
        outgroup_accession=outgroup_accession,
        metadata_path=metadata_path,
        nodownload=nodownload,
        limit=limit,
        drop_foreign=drop_foreign,
    )


def genome_params(
    *, accession_list_only: Any = _UNSET, keep_files: Any = _UNSET, sketch: Any = _UNSET
) -> GenomeParams:
    from ..stages.genome import GenomeParams

    return _build(
        GenomeParams, accession_list_only=accession_list_only, keep_files=keep_files, sketch=sketch
    )


def _absolute(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def ingest_params(
    *,
    genomes_dir: str | None = None,
    from_workdirs: Any = _UNSET,
    selection: Any = _UNSET,
    outgroup: Any = _UNSET,
    copy: Any = _UNSET,
    drop_foreign: Any = _UNSET,
    sketch: Any = _UNSET,
    strict_sources: Any = _UNSET,
) -> IngestParams:
    from ..core.errors import UserInputError
    from ..stages.ingest import IngestParams

    workdirs = [str(wd) for wd in from_workdirs] if from_workdirs is not _UNSET else []
    if genomes_dir is None and not workdirs:
        raise UserInputError("ingest needs --genomes-dir, --from-workdir, or both.")
    if isinstance(selection, str) and genomes_dir is None:
        raise UserInputError(
            "--selection names genomes under --genomes-dir; give --genomes-dir with it "
            "(a --from-workdir brings its own selection.tsv)."
        )

    # The paths are recorded and digested as resume inputs, and doctor
    # re-derives the inputs from the record, so a relative path must not
    # depend on the directory the command was started from. abspath, not
    # resolve(): see the firmlink note in stages.ingest._stage. --outgroup is
    # made absolute only when it names a file; otherwise it is a genome name.
    if genomes_dir is not None:
        genomes_dir = _absolute(genomes_dir)
    if isinstance(selection, str):
        selection = _absolute(selection)
    if isinstance(outgroup, str) and os.path.isfile(os.path.expanduser(outgroup)):
        outgroup = _absolute(outgroup)
    return _build(
        IngestParams,
        genomes_dir=genomes_dir,
        from_workdirs=[_absolute(wd) for wd in workdirs],
        selection=selection,
        outgroup=outgroup,
        copy=copy,
        drop_foreign=drop_foreign,
        sketch=sketch,
        strict_sources=strict_sources,
    )


def reads_params(
    *,
    target_family: Any = _UNSET,
    target_genus: Any = _UNSET,
    target_species: Any = _UNSET,
    accessions: Any = _UNSET,
    accession_file: Any = _UNSET,
    platform: Any = _UNSET,
    max_runs: Any = _UNSET,
    min_bases: Any = _UNSET,
    max_bases: Any = _UNSET,
    drop_selection: Any = _UNSET,
    one_per_sample: Any = _UNSET,
) -> ReadsParams:
    from ..stages.reads import PLATFORMS, ReadsParams, validate

    if platform is not _UNSET:
        _require_choice(platform, set(PLATFORMS), "--platform")
    params = _build(
        ReadsParams,
        target_family=target_family,
        target_genus=target_genus,
        target_species=target_species,
        accessions=accessions,
        accession_file=accession_file,
        platform=platform,
        max_runs=max_runs,
        min_bases=min_bases,
        max_bases=max_bases,
        drop_selection=drop_selection,
        one_per_sample=one_per_sample,
    )
    # Parse the accession file and check the selection now: the builder runs
    # before the workdir is created, so a rejected call leaves nothing behind.
    validate(params)
    return params


def assemble_params(
    *,
    assembler: Any = _UNSET,
    threads: Any = _UNSET,
    jobs: Any = _UNSET,
    memory_gb: Any = _UNSET,
    min_contig_length: Any = _UNSET,
    polisher: Any = _UNSET,
    polish_rounds: Any = _UNSET,
    outgroup: Any = _UNSET,
    append: Any = _UNSET,
    keep_reads: Any = _UNSET,
    keep_files: Any = _UNSET,
    checkm2_db: Any = _UNSET,
    min_completeness: Any = _UNSET,
    max_contamination: Any = _UNSET,
    classifier: Any = _UNSET,
    gtdb_sketch: Any = _UNSET,
    gtdb_lineages: Any = _UNSET,
    sketch: Any = _UNSET,
    reads_sketch: Any = _UNSET,
    extra: Any = _UNSET,
    max_runs: Any = _UNSET,
    accessions: Any = _UNSET,
    accession_file: Any = _UNSET,
    screen_reads: Any = _UNSET,
    screen_min_fraction: Any = _UNSET,
    screen_dup_ani: Any = _UNSET,
) -> AssembleParams:
    from ..assemblers.base import registry as _asm_registry
    from ..classifiers.base import registry as _cls_registry
    from ..polishers.base import registry as _pol_registry
    from ..stages.assemble import AssembleParams

    if assembler is not _UNSET:
        _require_choice(assembler, {"auto", *_asm_registry.names()}, "--assembler")
    if classifier is not _UNSET:
        _require_choice(classifier, {"auto", "none", *_cls_registry.names()}, "--classifier")
    if polisher is not _UNSET:
        _require_choice(polisher, {"auto", "none", *_pol_registry.names()}, "--polisher")
    return _build(
        AssembleParams,
        assembler=assembler,
        threads=threads,
        jobs=jobs,
        memory_gb=memory_gb,
        min_contig_length=min_contig_length,
        polisher=polisher,
        polish_rounds=polish_rounds,
        outgroup=outgroup,
        append=append,
        keep_reads=keep_reads,
        keep_files=keep_files,
        checkm2_db=checkm2_db,
        min_completeness=min_completeness,
        max_contamination=max_contamination,
        classifier=classifier,
        gtdb_sketch=gtdb_sketch,
        gtdb_lineages=gtdb_lineages,
        sketch=sketch,
        reads_sketch=reads_sketch,
        extra=extra,
        max_runs=max_runs,
        accessions=accessions,
        accession_file=accession_file,
        screen_reads=screen_reads,
        screen_min_fraction=screen_min_fraction,
        screen_dup_ani=screen_dup_ani,
    )


def vmetadata_params(
    *,
    target: Any = _UNSET,
    source: Any = _UNSET,
    filter: Any = _UNSET,
    list_targets: Any = _UNSET,
    host: Any = _UNSET,
    complete_only: Any = _UNSET,
    released_after: Any = _UNSET,
) -> VmetadataParams:
    from ..stages.vmetadata import VmetadataParams

    if source is not _UNSET:
        _require_choice(source, VIRAL_SOURCES, "--source")
    effective_source = source if source is not _UNSET else "ncbi_virus"
    if filter not in (_UNSET, None) and effective_source != "bvbrc":
        from ..core.errors import UserInputError

        raise UserInputError(
            "--filter is a BV-BRC header tag and has no effect on the NCBI Virus source; "
            "use --complete-only (and --host, --released-after) there, or --source bvbrc."
        )
    return _build(
        VmetadataParams,
        target=target,
        source=source,
        filter=filter,
        list_targets=list_targets,
        host=host,
        complete_only=complete_only,
        released_after=released_after,
    )


def vgenome_params(
    *,
    target_genus: Any = _UNSET,
    target_species: Any = _UNSET,
    target_serotype: Any = _UNSET,
    target_custom: Any = _UNSET,
    length_all: Any = _UNSET,
    length_deviation: Any = _UNSET,
    length_method: Any = _UNSET,
    length_range: Any = _UNSET,
    discard: Any = _UNSET,
    no_outgroup: Any = _UNSET,
    outgroup_accession: Any = _UNSET,
    group_segments: Any = _UNSET,
    outgroup_candidates_taxid_min_genomes: Any = _UNSET,
    outgroup_treebuilder: Any = _UNSET,
    glance: Any = _UNSET,
    print_fasta_headers: Any = _UNSET,
    ignore_duplicates: Any = _UNSET,
    keep_files: Any = _UNSET,
    sketch: Any = _UNSET,
) -> VgenomeParams:
    from ..stages.vgenome import VgenomeParams
    from ..viral._outgroup import distance_matrix_builders

    if length_method is not _UNSET:
        _require_choice(length_method, LENGTH_METHODS, "--length-method")
    if outgroup_treebuilder is not _UNSET:
        _require_choice(
            outgroup_treebuilder, set(distance_matrix_builders()), "--outgroup-treebuilder"
        )
    return _build(
        VgenomeParams,
        target_genus=target_genus,
        target_species=target_species,
        target_serotype=target_serotype,
        target_custom=target_custom,
        length_all=length_all,
        length_deviation=length_deviation,
        length_method=length_method,
        length_range=length_range,
        discard=discard,
        no_outgroup=no_outgroup,
        outgroup_accession=outgroup_accession,
        group_segments=group_segments,
        outgroup_candidates_taxid_min_genomes=outgroup_candidates_taxid_min_genomes,
        outgroup_treebuilder=outgroup_treebuilder,
        glance=glance,
        print_fasta_headers=print_fasta_headers,
        ignore_duplicates=ignore_duplicates,
        keep_files=keep_files,
        sketch=sketch,
    )


def dereplicate_params(
    *,
    tool: Any = _UNSET,
    primary_ani: Any = _UNSET,
    secondary_ani: Any = _UNSET,
    aligned_fraction: Any = _UNSET,
    threads: Any = _UNSET,
    process_size: Any = _UNSET,
    num_processes: Any = _UNSET,
    pre_primary_ani: Any = _UNSET,
    pre_secondary_ani: Any = _UNSET,
    reduce: Any = _UNSET,
    target_reps: Any = _UNSET,
    extra: Any = _UNSET,
    allow_incomplete: Any = _UNSET,
    keeper: Any = _UNSET,
    refine: Any = _UNSET,
    refine_ani: Any = _UNSET,
    refine_file: Any = _UNSET,
) -> DereplicateParams:
    from ..core.errors import UserInputError
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.dereplicate import DereplicateParams

    if tool is not _UNSET:
        _require_choice(tool, {"auto", *_derep_registry.names()}, "--tool")
    if reduce is not _UNSET:
        _require_choice(reduce, {"none", "species", "genus"}, "--reduce")
    if keeper is not _UNSET:
        _require_choice(keeper, {"quality", "gtdb", "tool"}, "--keeper")
    if target_reps is not _UNSET and target_reps < 0:
        raise UserInputError(f"--target-reps must be >= 0, got {target_reps}.")
    for value, label in (
        (primary_ani, "--primary-ani"),
        (secondary_ani, "--secondary-ani"),
        (aligned_fraction, "--aligned-fraction"),
        (pre_primary_ani, "--pre-primary-ani"),
        (pre_secondary_ani, "--pre-secondary-ani"),
    ):
        if value is not _UNSET:
            _require_unit_interval(value, label)
    _check_refine(refine, refine_ani)
    refining = (refine is not _UNSET and refine) or (refine_file is not _UNSET and refine_file)
    if refining and target_reps is not _UNSET and target_reps > 0:
        raise UserInputError("--refine cannot be combined with --target-reps.")
    return _build(
        DereplicateParams,
        tool=tool,
        primary_ani=primary_ani,
        secondary_ani=secondary_ani,
        aligned_fraction=aligned_fraction,
        threads=threads,
        process_size=process_size,
        num_processes=num_processes,
        pre_primary_ani=pre_primary_ani,
        pre_secondary_ani=pre_secondary_ani,
        reduce=reduce,
        target_reps=target_reps,
        extra=extra,
        allow_incomplete=allow_incomplete,
        keeper=keeper,
        refine=refine,
        refine_ani=refine_ani,
        refine_file=refine_file,
    )


def _check_refine(refine: Any, refine_ani: Any) -> None:
    """A --refine seed needs --refine-ani, a unit fraction; the stage checks
    that it is finer than the secondary ANI and warns on an ANI without
    seeds."""
    from ..core.errors import UserInputError

    seeds = [] if refine is _UNSET or refine is None else list(refine)
    ani = None if refine_ani is _UNSET else refine_ani
    if seeds and ani is None:
        raise UserInputError("--refine needs --refine-ani.")
    if ani is not None:
        _require_unit_interval(ani, "--refine-ani")


def derep_refine_params(
    *,
    seeds: Any = _UNSET,
    ani: Any = _UNSET,
    refine_file: Any = _UNSET,
    threads: Any = _UNSET,
):
    from ..core.errors import UserInputError
    from ..stages.derep_refine import DerepRefineParams

    _check_refine(seeds, ani)
    if not ((seeds is not _UNSET and seeds) or (refine_file is not _UNSET and refine_file)):
        raise UserInputError(
            "derep-refine needs --refine seeds (with --refine-ani) or --refine-file."
        )
    return _build(DerepRefineParams, seeds=seeds, ani=ani, refine_file=refine_file, threads=threads)


# The per-gene selection of the core-gene route (coregenes/select.py): the
# stage extra each flag fills.
GENE_SELECTION_FLAGS = (
    ("max_core_genes", "--max-core-genes"),
    ("gene_rank", "--gene-rank"),
    ("min_gene_length", "--min-gene-length"),
    ("max_gene_gaps", "--max-gene-gaps"),
)


def gene_selection_extra(
    *,
    max_core_genes: int | None,
    gene_rank: str | None,
    min_gene_length: int | None,
    max_gene_gaps: float | None,
) -> dict[str, Any]:
    """The gene-selection options that were given, as phylo stage extras."""
    given = {
        "max_core_genes": max_core_genes,
        "gene_rank": gene_rank,
        "min_gene_length": min_gene_length,
        "max_gene_gaps": max_gene_gaps,
    }
    return {key: value for key, value in given.items() if value is not None}


def require_gene_rank(rank: str) -> None:
    """Refuse a --gene-rank that is not one of the ranking criteria."""
    from ..coregenes.select import RANKS

    _require_choice(rank, set(RANKS), "--gene-rank")


def phylo_params(
    *,
    treebuilder: Any = _UNSET,
    msa_source: Any = _UNSET,
    aligner: Any = _UNSET,
    snptyper: Any = _UNSET,
    coregene: Any = _UNSET,
    all_genomes: Any = _UNSET,
    no_outgroup: Any = _UNSET,
    root: Any = _UNSET,
    bootstrap: Any = _UNSET,
    reference: Any = _UNSET,
    threads: Any = _UNSET,
    extra: Any = _UNSET,
    allow_incomplete: Any = _UNSET,
    on_run: bool = False,
) -> PhyloParams:
    from ..aligners.base import registry as _aln_registry
    from ..core.errors import UserInputError
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloParams, check_root_request
    from ..tree.rooting import ROOT_METHODS
    from ..treebuilders.base import registry as _tb_registry

    if treebuilder is not _UNSET:
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
    if root is not _UNSET and root is not None:
        _require_choice(root, set(ROOT_METHODS), "--root")
        check_root_request(root, no_outgroup is not _UNSET and bool(no_outgroup))
    effective_source = msa_source if msa_source is not _UNSET else "aligner"
    _require_choice(effective_source, {"aligner", "snptype", "coregene"}, "--msa-source")
    if effective_source == "aligner":
        if aligner is not _UNSET:
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
    elif effective_source == "snptype":
        if snptyper is not _UNSET:
            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
    elif coregene is not _UNSET:
        from ..coregenes.base import registry as _cg_registry

        _require_choice(coregene, set(_cg_registry.names()), "--coregene")
    # The core-gene settings are read by that pass only.
    for key, flag in (
        ("core_threshold", "--core-threshold"),
        ("clean_mode", "--clean-mode"),
        *GENE_SELECTION_FLAGS,
    ):
        if extra is not _UNSET and extra.get(key) is not None and effective_source != "coregene":
            raise UserInputError(f"{flag} applies only with --msa-source coregene.")
    if extra is not _UNSET and extra.get("gene_rank") is not None:
        require_gene_rank(str(extra["gene_rank"]))
    if extra is not _UNSET and extra.get("clean_mode") is not None:
        from ..coregenes.base import CLEAN_MODES

        _require_choice(str(extra["clean_mode"]), set(CLEAN_MODES), "--clean-mode")
    # The phylo stage reads extra["mask"] only on the snptype path; a mask
    # requested with the aligner source would be dropped without a trace.
    mask = extra.get("mask") if extra is not _UNSET else None
    if mask is not None:
        require_mask(mask)
        if effective_source != "snptype":
            hint = ", or on run with --with-snptype" if on_run else ""
            raise UserInputError(f"--mask applies only with --msa-source snptype{hint}.")
    # Likewise --core-fraction: the typing pass filters the core, an aligner
    # alignment is left as it is, so the value would be dropped silently.
    core_fraction = extra.get("core_fraction") if extra is not _UNSET else None
    if core_fraction is not None and effective_source != "snptype":
        hint = ", or on run with --with-snptype" if on_run else ""
        raise UserInputError(f"--core-fraction applies only with --msa-source snptype{hint}.")
    return _build(
        PhyloParams,
        treebuilder=treebuilder,
        msa_source=msa_source,
        aligner=aligner,
        snptyper=snptyper,
        coregene=coregene,
        all_genomes=all_genomes,
        no_outgroup=no_outgroup,
        root=root,
        bootstrap=bootstrap,
        reference=reference,
        threads=threads,
        extra=extra,
        allow_incomplete=allow_incomplete,
    )


def tree2tax_params(
    *,
    node_basename: Any = _UNSET,
    naming: Any = _UNSET,
    root_name: Any = _UNSET,
    remove_outgroup: Any = _UNSET,
    include_dereplicated: Any = _UNSET,
    collapse_support: Any = _UNSET,
    collapse_length: Any = _UNSET,
    collapse_diameter: Any = _UNSET,
    ranks: Any = _UNSET,
    rank_red: Any = _UNSET,
) -> Tree2taxParams:
    from ..stages.tree2tax import NAMING_RULES, Tree2taxParams

    if naming is not _UNSET:
        _require_choice(naming, set(NAMING_RULES), "--naming")
    return _build(
        Tree2taxParams,
        node_basename=node_basename,
        naming=naming,
        root_name=root_name,
        remove_outgroup=remove_outgroup,
        include_dereplicated=include_dereplicated,
        collapse_support=collapse_support,
        collapse_length=collapse_length,
        collapse_diameter=collapse_diameter,
        ranks=ranks,
        rank_red=rank_red,
    )
