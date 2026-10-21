"""Stateless data-channel steps: genome-fetch, dereplicate-chunk, dereplicate-merge,
assemble-run, genome-qc, reads-gather.

These run as discrete Nextflow process steps (no shared workdir); they read
explicit inputs (a selection.tsv or a file-of-filenames) and write a result
directory, rather than going through the workdir-bound :func:`_run` harness.
"""

from __future__ import annotations

from pathlib import Path

import typer

from ..core.errors import UserInputError
from ..core.logging import configure_logging
from .base import (
    _RUN_STATE,
    DEFAULT_THREADS,
    HELP_ALIGNED_FRACTION,
    HELP_ALIGNER_ARG,
    HELP_ANNOTATOR,
    HELP_BAKTA_DB,
    HELP_BOOTSTRAP,
    HELP_CLEAN_MODE,
    HELP_COLLAPSE_DIAMETER,
    HELP_COLLAPSE_LENGTH,
    HELP_COLLAPSE_SUPPORT,
    HELP_CORE_FRACTION,
    HELP_CORE_THRESHOLD,
    HELP_COREGENE,
    HELP_COREGENE_ALPHABET,
    HELP_DEREP_TOOL_ARG,
    HELP_GENE_RANK,
    HELP_GENE_TRIM,
    HELP_INCLUDE_DEREPLICATED,
    HELP_KEEP_FILES,
    HELP_KEEP_ROOT,
    HELP_MAX_CORE_GENES,
    HELP_MAX_GENE_GAPS,
    HELP_MIN_GENE_LENGTH,
    HELP_MSA_SOURCE,
    HELP_NAMING,
    HELP_NO_OUTGROUP,
    HELP_NODE_BASENAME,
    HELP_PARTITIONED,
    HELP_PRESET,
    HELP_PRIMARY_ANI,
    HELP_RANK_RED,
    HELP_RANKS,
    HELP_READS_SKETCH,
    HELP_REFERENCE,
    HELP_REFINE,
    HELP_REFINE_ANI,
    HELP_REFINE_FILE,
    HELP_REMOVE_OUTGROUP,
    HELP_ROOT,
    HELP_ROOT_NAME,
    HELP_SCREEN_MIN_FRACTION,
    HELP_SCREEN_READS,
    HELP_SECONDARY_ANI,
    HELP_SKETCHES_DIR,
    HELP_THREADS,
    HELP_VERSIONS_OUT,
    HELP_VIRUS,
    PANEL_STEPS,
    _aligner_help,
    _assembler_help,
    _classifier_help,
    _derep_help,
    _mask_help_msa,
    _parse_key_values,
    _polisher_help,
    _read_path_fofn,
    _require_choice,
    _require_unit_interval,
    _snp_help,
    _tree_help,
    app,
    gated_extra,
    resolve_threads,
    stage_errors,
)


@app.command(name="genome-fetch", rich_help_panel=PANEL_STEPS)
def genome_fetch_cmd(
    selection: Path = typer.Option(
        ..., "--selection", help="selection.tsv from the metadata stage."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for downloaded genomes."),
    keep_files: bool = typer.Option(False, "--keep-files", help=HELP_KEEP_FILES),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Download genomes listed in a selection.tsv (stateless data-channel step)."""
    from ..stages.genome_steps import GenomeFetchParams, genome_fetch

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        genome_fetch(
            GenomeFetchParams(
                selection_tsv=selection,
                out_dir=out_dir,
                keep_files=keep_files,
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="dereplicate-chunk", rich_help_panel=PANEL_STEPS)
def dereplicate_chunk_cmd(
    genomes_fofn: Path = typer.Option(
        ..., "--genomes-fofn", help="File of genome FASTA paths, one per line."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output directory for the chunk result."),
    tool: str = typer.Option("skder", "--tool", help=_derep_help(auto=False)),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(
        0.50, "-af", "--aligned-fraction", help=HELP_ALIGNED_FRACTION
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    virus: bool = typer.Option(False, "--virus", help=HELP_VIRUS),
    tool_arg: list[str] = typer.Option([], "--tool-arg", help=HELP_DEREP_TOOL_ARG),
    selection_tsv: Path | None = typer.Option(
        None,
        "--selection-tsv",
        help="selection.tsv with quality columns; enables quality-aware representatives.",
    ),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help="Representative choice when --selection-tsv is given: "
        "quality (selection.tsv CheckM values and genome N50), gtdb (GTDB species "
        "representative first, then quality) or tool (adapter's own pick).",
    ),
    sketches_dir: Path | None = typer.Option(None, "--sketches-dir", help=HELP_SKETCHES_DIR),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Dereplicate one chunk of genomes (scatter step; writes a chunk result dir)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import ChunkParams, dereplicate_chunk

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        # A concrete tool (not 'auto') so every scattered chunk and the merge agree.
        _require_choice(tool, set(_derep_registry.names()), "--tool")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_choice(keeper, {"quality", "gtdb", "tool"}, "--keeper")
        genomes = _read_path_fofn(genomes_fofn)
        dereplicate_chunk(
            ChunkParams(
                tool=tool,
                genomes=genomes,
                out_dir=out_dir,
                primary_ani=primary_ani,
                secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction,
                threads=threads,
                extra={
                    **_parse_key_values(tool_arg, "--tool-arg"),
                    **(
                        gated_extra(_derep_registry, tool, "virus", True, flag="--virus")
                        if virus
                        else {}
                    ),
                },
                selection_tsv=selection_tsv,
                keeper=keeper,
                versions_out=versions_out,
                sketches_dir=sketches_dir,
            ),
            logger,
        )


@app.command(name="annotate-build", rich_help_panel=PANEL_STEPS)
def annotate_build_cmd(
    genome: Path = typer.Option(..., "--genome", help="Genome FASTA to annotate (gz accepted)."),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output directory; the record's files go under <out>/<name>/."
    ),
    name: str | None = typer.Option(
        None, "--record-name", help="Record name (prefix of the outputs); the genome's by default."
    ),
    tool: str = typer.Option("bakta", "--tool", help=HELP_ANNOTATOR),
    bakta_db: str | None = typer.Option(None, "--bakta-db", help=HELP_BAKTA_DB),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    tool_arg: list[str] = typer.Option([], "--tool-arg", help="Annotator tuning as key=value."),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Annotate one genome (stateless data-channel step; writes <out>/<name>/<name>.gff3)."""
    from ..annotators.base import registry as _annot_registry
    from ..stages.annotate import AnnotateBuildParams, annotate_build

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(tool, set(_annot_registry.names()), "--tool")
        annotate_build(
            AnnotateBuildParams(
                genome=genome,
                out_dir=out_dir,
                record_name=name,
                tool=tool,
                bakta_db=bakta_db,
                threads=threads,
                extra=_parse_key_values(tool_arg, "--tool-arg"),
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="annotate-index", rich_help_panel=PANEL_STEPS)
def annotate_index_cmd(
    annotations_dir: Path = typer.Option(
        ...,
        "--annotations-dir",
        help="Directory holding the record directories annotate-build wrote "
        "(<dir>/<name>/<name>.gff3).",
    ),
    tool: str = typer.Option(
        "bakta",
        "--tool",
        help="Annotator assumed for a record without annotation.json (output of an "
        "earlier version): bakta or prokka. Each other record is indexed with the "
        "annotator it names.",
    ),
) -> None:
    """Index gathered annotate-build outputs (stateless data-channel step; writes
    <dir>/annotations.tsv, which phylo-build --annotations-dir reads)."""
    from ..annotators.base import registry as _annot_registry
    from ..stages.annotate import AnnotateIndexParams, annotate_index

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(tool, set(_annot_registry.names()), "--tool")
        annotate_index(AnnotateIndexParams(annotations_dir=annotations_dir, tool=tool), logger)


@app.command(name="phylo-build", rich_help_panel=PANEL_STEPS)
def phylo_build_cmd(
    ctx: typer.Context,
    genomes_dir: Path = typer.Option(
        ..., "--genomes-dir", help="Directory of genome FASTA files to build the tree from."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir (writes tree/tree.nwk)."),
    outgroup_dir: Path | None = typer.Option(
        None, "--outgroup-dir", help="Directory holding the outgroup genome file(s)."
    ),
    outgroup_accession: Path | None = typer.Option(
        None, "--outgroup-accession", help="File naming the outgroup accession."
    ),
    preset: str | None = typer.Option(None, "--preset", help=HELP_PRESET),
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help=HELP_MSA_SOURCE),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help=_snp_help()),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help=HELP_NO_OUTGROUP),
    root: str | None = typer.Option(None, "--root", help=HELP_ROOT),
    bootstrap: int = typer.Option(0, "-B", "--bootstrap", min=0, help=HELP_BOOTSTRAP),
    core_fraction: float | None = typer.Option(
        None, "--core-fraction", min=0.0, max=1.0, help=HELP_CORE_FRACTION
    ),
    coregene: str = typer.Option("panaroo", "--coregene", help=HELP_COREGENE),
    core_threshold: float | None = typer.Option(
        None, "--core-threshold", min=0.0, max=1.0, help=HELP_CORE_THRESHOLD
    ),
    clean_mode: str | None = typer.Option(None, "--clean-mode", help=HELP_CLEAN_MODE),
    coregene_alphabet: str | None = typer.Option(
        None, "--coregene-alphabet", help=HELP_COREGENE_ALPHABET
    ),
    gene_trim: str | None = typer.Option(None, "--gene-trim", help=HELP_GENE_TRIM),
    max_core_genes: int | None = typer.Option(
        None, "--max-core-genes", min=0, help=HELP_MAX_CORE_GENES
    ),
    gene_rank: str | None = typer.Option(None, "--gene-rank", help=HELP_GENE_RANK),
    min_gene_length: int | None = typer.Option(
        None, "--min-gene-length", min=0, help=HELP_MIN_GENE_LENGTH
    ),
    max_gene_gaps: float | None = typer.Option(
        None, "--max-gene-gaps", min=0.0, max=1.0, help=HELP_MAX_GENE_GAPS
    ),
    partitioned: bool = typer.Option(False, "--partitioned", help=HELP_PARTITIONED),
    annotations_dir: Path | None = typer.Option(
        None,
        "--annotations-dir",
        help="With --msa-source coregene: an annotations/ directory the annotate stage or "
        "annotate-build wrote (its index names the records); without it the step annotates "
        "the genomes itself with --annotator.",
    ),
    annotator: str = typer.Option("bakta", "--annotator", help=HELP_ANNOTATOR),
    bakta_db: str | None = typer.Option(None, "--bakta-db", help=HELP_BAKTA_DB),
    reference: str | None = typer.Option(None, "--reference", help=HELP_REFERENCE),
    aligner_arg: list[str] = typer.Option([], "--aligner-arg", help=HELP_ALIGNER_ARG),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    mask: str = typer.Option("none", "--mask", help=_mask_help_msa()),
    msa_only: bool = typer.Option(
        False,
        "--msa-only",
        help="Build the alignment and stop, writing msa.fasta (for a separate tree step).",
    ),
    msa: Path | None = typer.Option(
        None, "--msa", help="Build the tree from this alignment instead of constructing one."
    ),
    sketches_dir: Path | None = typer.Option(None, "--sketches-dir", help=HELP_SKETCHES_DIR),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Build a phylogeny from a genomes directory (stateless data-channel step).

    Here --outgroup-accession takes a file that names the accession, not the accession itself.
    """
    from ..aligners.base import registry as _aln_registry
    from ..snptypers.base import registry as _snp_registry
    from ..stages.phylo import PhyloBuildParams, PhyloParams, check_root_request, phylo_build
    from ..tree.rooting import ROOT_METHODS
    from ..treebuilders.base import registry as _tb_registry
    from .param_builders import (
        apply_preset,
        gene_selection_extra,
        refuse_outside_coregene,
        require_gene_rank,
        require_mask,
    )

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        p = apply_preset(ctx, preset)
        treebuilder, msa_source, snptyper = p["treebuilder"], p["msa_source"], p["snptyper"]
        bootstrap, aligner_arg, mask = p["bootstrap"], p["aligner_arg"], p["mask"]
        core_fraction, coregene = p["core_fraction"], p["coregene"]
        core_threshold, max_core_genes = p["core_threshold"], p["max_core_genes"]
        gene_rank, partitioned = p["gene_rank"], p["partitioned"]
        coregene_alphabet = p["coregene_alphabet"]
        _require_choice(treebuilder, {"auto", *_tb_registry.names()}, "--treebuilder")
        _require_choice(msa_source, {"aligner", "snptype", "coregene"}, "--msa-source")
        if root is not None:
            _require_choice(root, set(ROOT_METHODS), "--root")
            check_root_request(root, no_outgroup)
        if msa_source == "aligner":
            _require_choice(aligner, set(_aln_registry.names()), "--aligner")
        elif msa_source == "snptype":
            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
        else:
            from ..annotators.base import registry as _annot_registry
            from ..coregenes.base import registry as _cg_registry

            _require_choice(coregene, set(_cg_registry.names()), "--coregene")
            _require_choice(annotator, set(_annot_registry.names()), "--annotator")
        refuse_outside_coregene(
            msa_source,
            {
                "--core-threshold": core_threshold,
                "--clean-mode": clean_mode,
                "--coregene-alphabet": coregene_alphabet,
                "--gene-trim": gene_trim,
                "--max-core-genes": max_core_genes,
                "--gene-rank": gene_rank,
                "--min-gene-length": min_gene_length,
                "--max-gene-gaps": max_gene_gaps,
            },
        )
        # With --msa the sidecar beside the file decides, whatever the source.
        if partitioned and msa_source != "coregene" and msa is None:
            raise UserInputError(
                "--partitioned applies only with --msa-source coregene, or with --msa FILE "
                "whose partitions.json lies beside it."
            )
        if gene_rank is not None:
            require_gene_rank(gene_rank)
        if clean_mode is not None:
            from ..coregenes.base import CLEAN_MODES

            _require_choice(clean_mode, set(CLEAN_MODES), "--clean-mode")
        if gene_trim is not None:
            from ..coregenes.trim import GENE_TRIMS

            _require_choice(gene_trim, set(GENE_TRIMS), "--gene-trim")
        if coregene_alphabet is not None:
            from ..coregenes.select import ALPHABETS

            _require_choice(coregene_alphabet, set(ALPHABETS), "--coregene-alphabet")
            if msa is None:
                from ..stages.coregene import check_alphabet

                check_alphabet(coregene, coregene_alphabet)
        require_mask(mask)
        if mask != "none" and msa_source != "snptype":
            raise UserInputError("--mask applies only with --msa-source snptype.")
        # Not refused with --msa: Nextflow passes the same phylo_args to the
        # tree half, where the typing pass does not run and the value is idle.
        if core_fraction is not None and msa_source != "snptype" and msa is None:
            raise UserInputError("--core-fraction applies only with --msa-source snptype.")
        if msa_only and msa is not None:
            raise UserInputError("--msa-only builds an alignment; --msa consumes one. Pick one.")
        if msa is not None and not msa.is_file():
            raise UserInputError(f"Alignment not found: {msa}")

        phylo_params = PhyloParams(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            no_outgroup=no_outgroup,
            root=root,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            coregene=coregene,
            extra={
                **_parse_key_values(aligner_arg, "--aligner-arg"),
                **({"mask": mask} if mask != "none" else {}),
                **({"core_fraction": core_fraction} if core_fraction is not None else {}),
                **({"core_threshold": core_threshold} if core_threshold is not None else {}),
                **({"clean_mode": clean_mode} if clean_mode is not None else {}),
                **(
                    {"coregene_alphabet": coregene_alphabet}
                    if coregene_alphabet is not None
                    else {}
                ),
                **({"gene_trim": gene_trim} if gene_trim is not None else {}),
                **gene_selection_extra(
                    max_core_genes=max_core_genes,
                    gene_rank=gene_rank,
                    min_gene_length=min_gene_length,
                    max_gene_gaps=max_gene_gaps,
                ),
                **({"partitioned": True} if partitioned else {}),
            },
            preset=preset,
        )
        phylo_build(
            PhyloBuildParams(
                genomes_dir=genomes_dir,
                out_dir=out_dir,
                outgroup_dir=outgroup_dir,
                outgroup_accession=outgroup_accession,
                phylo=phylo_params,
                versions_out=versions_out,
                msa_only=msa_only,
                msa=msa,
                sketches_dir=sketches_dir,
                annotations_dir=annotations_dir,
                annotator=annotator,
                bakta_db=bakta_db,
            ),
            logger,
        )


@app.command(name="tree2tax-relations", rich_help_panel=PANEL_STEPS)
def tree2tax_relations_cmd(
    tree: Path = typer.Option(..., "--tree", help="Rooted/unrooted tree in Newick (tree.nwk)."),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir (writes tree2tax.tsv + genomes_map.tsv)."
    ),
    clusters: Path | None = typer.Option(
        None, "--clusters", help="derep clusters.tsv (for --include-dereplicated)."
    ),
    outgroup_dir: Path | None = typer.Option(
        None, "--outgroup-dir", help="Directory holding the outgroup genome file(s)."
    ),
    outgroup_accession: Path | None = typer.Option(
        None, "--outgroup-accession", help="File naming the outgroup accession."
    ),
    segments: Path | None = typer.Option(
        None,
        "--segments",
        help="segments.tsv from vgenome --group-segments: member segment accessions "
        "are listed under their isolate's leaf.",
    ),
    node_basename: str | None = typer.Option(None, "--node-basename", help=HELP_NODE_BASENAME),
    naming: str = typer.Option("leafset", "--naming", help=HELP_NAMING),
    root_name: str = typer.Option("root", "--root-name", help=HELP_ROOT_NAME),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help=HELP_REMOVE_OUTGROUP),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help=HELP_NO_OUTGROUP),
    keep_root: bool = typer.Option(False, "--keep-root", help=HELP_KEEP_ROOT),
    include_dereplicated: bool = typer.Option(
        True,
        "--include-dereplicated/--no-include-dereplicated",
        help=HELP_INCLUDE_DEREPLICATED,
    ),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
    collapse_support: float | None = typer.Option(
        None,
        "--collapse-support",
        min=0.0,
        max=1.0,
        help=HELP_COLLAPSE_SUPPORT,
    ),
    collapse_length: float | None = typer.Option(
        None,
        "--collapse-length",
        min=0.0,
        help=HELP_COLLAPSE_LENGTH,
    ),
    collapse_diameter: float | None = typer.Option(
        None,
        "--collapse-diameter",
        min=0.0,
        help=HELP_COLLAPSE_DIAMETER,
    ),
    ranks: bool = typer.Option(False, "--ranks", help=HELP_RANKS),
    rank_red: str | None = typer.Option(None, "--rank-red", help=HELP_RANK_RED),
) -> None:
    """Emit FlexTaxD relations from a tree (stateless data-channel step).

    Here --outgroup-accession takes a file that names the accession, not the accession itself.
    """
    from ..stages.tree2tax import NAMING_RULES, Tree2taxStepParams, tree2tax_relations

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(naming, set(NAMING_RULES), "--naming")
        tree2tax_relations(
            Tree2taxStepParams(
                tree=tree,
                out_dir=out_dir,
                clusters=clusters,
                outgroup_dir=outgroup_dir,
                outgroup_accession=outgroup_accession,
                node_basename=node_basename,
                naming=naming,
                root_name=root_name,
                remove_outgroup=remove_outgroup,
                include_dereplicated=include_dereplicated,
                versions_out=versions_out,
                collapse_support=collapse_support,
                collapse_length=collapse_length,
                collapse_diameter=collapse_diameter,
                ranks=ranks,
                rank_red=rank_red,
                no_outgroup=no_outgroup,
                keep_root=keep_root,
                segments=segments,
            ),
            logger,
        )


@app.command(name="dereplicate-merge", rich_help_panel=PANEL_STEPS)
def dereplicate_merge_cmd(
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for the merged result."),
    chunk_dir: list[Path] = typer.Option(
        [], "--chunk-dir", help="A chunk result directory (repeatable)."
    ),
    chunk_fofn: Path | None = typer.Option(
        None, "--chunk-fofn", help="File listing chunk result directories, one per line."
    ),
    tool: str = typer.Option("skder", "--tool", help=_derep_help(auto=False)),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(
        0.50, "-af", "--aligned-fraction", help=HELP_ALIGNED_FRACTION
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    virus: bool = typer.Option(False, "--virus", help=HELP_VIRUS),
    tool_arg: list[str] = typer.Option([], "--tool-arg", help=HELP_DEREP_TOOL_ARG),
    selection_tsv: Path | None = typer.Option(
        None,
        "--selection-tsv",
        help="selection.tsv with quality columns; enables quality-aware representatives.",
    ),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help="Representative choice when --selection-tsv is given: "
        "quality (selection.tsv CheckM values and genome N50), gtdb (GTDB species "
        "representative first, then quality) or tool (adapter's own pick).",
    ),
    reduce: str = typer.Option(
        "none",
        "--reduce",
        help="Taxonomy-aware reduction after the merge: none, species, or genus "
        "(one representative per taxon; taxonomy from --selection-tsv or the filenames).",
    ),
    target_reps: int = typer.Option(
        0,
        "--target-reps",
        help="Target representative count: search --secondary-ani of the merge pass "
        "to land near it (0 = off; re-runs the merge per search step).",
    ),
    sketches_dir: Path | None = typer.Option(None, "--sketches-dir", help=HELP_SKETCHES_DIR),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Dereplicate the union of chunk representatives (gather step)."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import MergeParams, dereplicate_merge

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(tool, set(_derep_registry.names()), "--tool")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_choice(keeper, {"quality", "gtdb", "tool"}, "--keeper")
        _require_choice(reduce, {"none", "species", "genus"}, "--reduce")
        if target_reps < 0:
            raise UserInputError("--target-reps must be 0 (off) or a positive count.")
        chunk_dirs = list(chunk_dir)
        if chunk_fofn is not None:
            chunk_dirs += _read_path_fofn(chunk_fofn)
        if not chunk_dirs:
            raise UserInputError("Provide at least one --chunk-dir or a --chunk-fofn.")
        dereplicate_merge(
            MergeParams(
                tool=tool,
                chunk_dirs=chunk_dirs,
                out_dir=out_dir,
                primary_ani=primary_ani,
                secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction,
                threads=threads,
                extra={
                    **_parse_key_values(tool_arg, "--tool-arg"),
                    **(
                        gated_extra(_derep_registry, tool, "virus", True, flag="--virus")
                        if virus
                        else {}
                    ),
                },
                selection_tsv=selection_tsv,
                keeper=keeper,
                reduce=reduce,
                target_reps=target_reps,
                versions_out=versions_out,
                sketches_dir=sketches_dir,
            ),
            logger,
        )


@app.command(name="dereplicate-refine", rich_help_panel=PANEL_STEPS)
def dereplicate_refine_cmd(
    contract_dir: Path = typer.Option(
        ..., "--contract-dir", help="A dereplication result directory (the merge step's output)."
    ),
    genomes_fofn: Path = typer.Option(
        ..., "--genomes-fofn", help="File of every clustered genome's FASTA path, one per line."
    ),
    out_dir: Path = typer.Option(..., "-o", "--out", help="Output dir for the refined result."),
    tool: str = typer.Option("skder", "--tool", help=_derep_help(auto=False)),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(
        0.50, "-af", "--aligned-fraction", help=HELP_ALIGNED_FRACTION
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    virus: bool = typer.Option(False, "--virus", help=HELP_VIRUS),
    tool_arg: list[str] = typer.Option([], "--tool-arg", help=HELP_DEREP_TOOL_ARG),
    selection_tsv: Path | None = typer.Option(
        None,
        "--selection-tsv",
        help="selection.tsv with quality columns; enables quality-aware representatives.",
    ),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help="Representative choice when --selection-tsv is given: "
        "quality (selection.tsv CheckM values and genome N50), gtdb (GTDB species "
        "representative first, then quality) or tool (adapter's own pick).",
    ),
    refine: list[str] = typer.Option([], "--refine", help=HELP_REFINE),
    refine_ani: float | None = typer.Option(None, "--refine-ani", help=HELP_REFINE_ANI),
    refine_file: Path | None = typer.Option(None, "--refine-file", help=HELP_REFINE_FILE),
    sketches_dir: Path | None = typer.Option(None, "--sketches-dir", help=HELP_SKETCHES_DIR),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Refine the lineages of the seeds in a dereplication result at finer ANIs."""
    from ..dereplicators.base import registry as _derep_registry
    from ..stages.derep_steps import RefineStepParams, dereplicate_refine

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(tool, set(_derep_registry.names()), "--tool")
        _require_unit_interval(primary_ani, "--primary-ani")
        _require_unit_interval(secondary_ani, "--secondary-ani")
        _require_unit_interval(aligned_fraction, "--aligned-fraction")
        _require_choice(keeper, {"quality", "gtdb", "tool"}, "--keeper")
        if refine_ani is not None:
            _require_unit_interval(refine_ani, "--refine-ani")
        if not refine and refine_file is None:
            raise UserInputError(
                "Give at least one --refine seed (with --refine-ani) or --refine-file."
            )
        genomes = _read_path_fofn(genomes_fofn)
        dereplicate_refine(
            RefineStepParams(
                tool=tool,
                contract_dir=contract_dir,
                genomes=genomes,
                out_dir=out_dir,
                primary_ani=primary_ani,
                secondary_ani=secondary_ani,
                aligned_fraction=aligned_fraction,
                threads=threads,
                extra={
                    **_parse_key_values(tool_arg, "--tool-arg"),
                    **(
                        gated_extra(_derep_registry, tool, "virus", True, flag="--virus")
                        if virus
                        else {}
                    ),
                },
                selection_tsv=selection_tsv,
                keeper=keeper,
                seeds=refine,
                ani=refine_ani,
                refine_file=None if refine_file is None else str(refine_file.resolve()),
                versions_out=versions_out,
                sketches_dir=sketches_dir,
            ),
            logger,
        )


# --- reads chain ------------------------------------------------------------------


@app.command(name="assemble-run", rich_help_panel=PANEL_STEPS)
def assemble_run_cmd(
    reads_tsv: Path = typer.Option(..., "--reads-tsv", help="reads.tsv from the reads stage."),
    run: str = typer.Option(
        ..., "--run", help="The run accession (a row of reads.tsv) to assemble."
    ),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir: contigs.fasta and assembly.ok, or excused_runs.tsv."
    ),
    assembler: str = typer.Option("auto", "--assembler", help=_assembler_help()),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    memory_gb: int = typer.Option(
        16,
        "--memory-gb",
        min=1,
        help="Memory hint for the assembly, in GB, for tools that cap RAM.",
    ),
    min_contig_length: int = typer.Option(
        500, "--min-contig-length", min=0, help="Drop contigs shorter than this many bases."
    ),
    polisher: str = typer.Option("auto", "--polisher", help=_polisher_help()),
    polish_rounds: int = typer.Option(
        1, "--polish-rounds", min=1, help="Polishing rounds (racon; medaka runs one)."
    ),
    keep_reads: bool = typer.Option(
        False, "--keep-reads", help="Keep the downloaded FASTQ files after assembling."
    ),
    keep_files: bool = typer.Option(
        False, "--keep-files", help="Keep the assembler scratch directory."
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Assembler tuning as key=value (repeatable), e.g. mode=nano-raw."
    ),
    reads_sketch: bool | None = typer.Option(
        None,
        "--reads-sketch/--no-reads-sketch",
        help=HELP_READS_SKETCH,
    ),
    screen_reads: bool = typer.Option(
        False,
        "--screen-reads",
        help=HELP_SCREEN_READS + " The duplicate check needs the other runs and is not made "
        "by this step.",
    ),
    screen_min_fraction: float = typer.Option(
        0.3, "--screen-min-fraction", min=0.0, max=1.0, help=HELP_SCREEN_MIN_FRACTION
    ),
    gtdb_sketch: Path | None = typer.Option(
        None,
        "--gtdb-sketch",
        help="GTDB sourmash sketch database for --screen-reads (or set REPGENR_GTDB_SKETCH).",
    ),
    gtdb_lineages: Path | None = typer.Option(
        None,
        "--gtdb-lineages",
        help="The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES).",
    ),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Fetch and assemble one run of a reads.tsv (stateless data-channel step)."""
    from ..assemblers.base import registry as asm_registry
    from ..polishers.base import registry as pol_registry
    from ..stages.assemble_steps import AssembleRunParams, assemble_run

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(assembler, {"auto", *asm_registry.names()}, "--assembler")
        _require_choice(polisher, {"auto", "none", *pol_registry.names()}, "--polisher")
        assemble_run(
            AssembleRunParams(
                reads_tsv=reads_tsv,
                run=run,
                out_dir=out_dir,
                assembler=assembler,
                threads=threads,
                memory_gb=memory_gb,
                min_contig_length=min_contig_length,
                polisher=polisher,
                polish_rounds=polish_rounds,
                keep_reads=keep_reads,
                keep_files=keep_files,
                extra=_parse_key_values(tool_arg, "--tool-arg"),
                reads_sketch=reads_sketch,
                screen_reads=screen_reads,
                screen_min_fraction=screen_min_fraction,
                gtdb_sketch=None if gtdb_sketch is None else str(gtdb_sketch),
                gtdb_lineages=None if gtdb_lineages is None else str(gtdb_lineages),
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="genome-qc", rich_help_panel=PANEL_STEPS)
def genome_qc_cmd(
    assemblies: Path = typer.Option(
        ..., "--assemblies", help="Directory of assemble-run output dirs, one per run."
    ),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir for quality.tsv and classification.tsv."
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    memory_gb: int = typer.Option(
        16,
        "--memory-gb",
        min=1,
        help="Memory budget in GB; caps concurrent classifier gathers at about 0.6 GB each.",
    ),
    checkm2_db: Path | None = typer.Option(
        None,
        "--checkm2-db",
        help="CheckM2 DIAMOND database; enables quality scoring (or set CHECKM2DB).",
    ),
    classifier: str = typer.Option("auto", "--classifier", help=_classifier_help()),
    gtdb_sketch: Path | None = typer.Option(
        None,
        "--gtdb-sketch",
        help="GTDB sourmash sketch database (.sig.zip); enables classification "
        "(or set REPGENR_GTDB_SKETCH).",
    ),
    gtdb_lineages: Path | None = typer.Option(
        None,
        "--gtdb-lineages",
        help="The lineages CSV published with the sketch (or set REPGENR_GTDB_LINEAGES).",
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Classifier tuning as key=value (repeatable)."
    ),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Score (CheckM2) and classify a batch of assemblies (stateless data-channel step)."""
    from ..classifiers.base import registry as cls_registry
    from ..stages.assemble_steps import GenomeQcParams, genome_qc

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        _require_choice(classifier, {"auto", "none", *cls_registry.names()}, "--classifier")
        genome_qc(
            GenomeQcParams(
                assemblies_dir=assemblies,
                out_dir=out_dir,
                threads=threads,
                memory_gb=memory_gb,
                checkm2_db=None if checkm2_db is None else str(checkm2_db),
                classifier=classifier,
                gtdb_sketch=None if gtdb_sketch is None else str(gtdb_sketch),
                gtdb_lineages=None if gtdb_lineages is None else str(gtdb_lineages),
                extra=_parse_key_values(tool_arg, "--tool-arg"),
                versions_out=versions_out,
            ),
            logger,
        )


@app.command(name="reads-gather", rich_help_panel=PANEL_STEPS)
def reads_gather_cmd(
    reads_tsv: Path = typer.Option(..., "--reads-tsv", help="reads.tsv from the reads stage."),
    assemblies: Path = typer.Option(
        ..., "--assemblies", help="Directory of assemble-run output dirs, one per run."
    ),
    out_dir: Path = typer.Option(
        ..., "-o", "--out", help="Output dir for genomes/, selection.tsv and the stats tables."
    ),
    qc: Path | None = typer.Option(
        None, "--qc", help="genome-qc output dir (quality.tsv, classification.tsv), if it ran."
    ),
    min_completeness: float = typer.Option(
        50.0, "--min-completeness", min=0.0, max=100.0, help="CheckM2 completeness floor."
    ),
    max_contamination: float = typer.Option(
        10.0, "--max-contamination", min=0.0, max=100.0, help="CheckM2 contamination ceiling."
    ),
) -> None:
    """Write the genome contract from per-run assemblies (stateless data-channel step)."""
    from ..stages.assemble_steps import ReadsGatherParams, reads_gather

    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        reads_gather(
            ReadsGatherParams(
                reads_tsv=reads_tsv,
                assemblies_dir=assemblies,
                out_dir=out_dir,
                qc_dir=qc,
                min_completeness=min_completeness,
                max_contamination=max_contamination,
            ),
            logger,
        )
