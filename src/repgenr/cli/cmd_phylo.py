"""Phylogenetics commands: snptype, phylo, tree2tax."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import (
    _RUN_STATE,
    DEFAULT_THREADS,
    HELP_ADVISE_JSON,
    HELP_ALIGNER_ARG,
    HELP_ALL_GENOMES,
    HELP_ALLOW_INCOMPLETE,
    HELP_ANNOTATE_PROCESSES,
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
    HELP_GENE_RANK,
    HELP_INCLUDE_DEREPLICATED,
    HELP_MAX_CORE_GENES,
    HELP_MAX_GENE_GAPS,
    HELP_MIN_GENE_LENGTH,
    HELP_MSA_SOURCE,
    HELP_NAMING,
    HELP_NO_OUTGROUP,
    HELP_NODE_BASENAME,
    HELP_PARTITIONED,
    HELP_RANK_RED,
    HELP_RANKS,
    HELP_REFERENCE,
    HELP_REMOVE_OUTGROUP,
    HELP_ROOT,
    HELP_ROOT_NAME,
    HELP_THREADS,
    HELP_WORKDIR,
    PANEL_CORE,
    PANEL_PIPELINE,
    _aligner_help,
    _echo_json,
    _mask_help,
    _mask_help_msa,
    _parse_key_values,
    _require_choice,
    _run,
    _snp_help,
    _tree_help,
    app,
    require_existing_workdir,
    resolve_threads,
    stage_errors,
)


@app.command(rich_help_panel=PANEL_CORE)
def annotate(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    tool: str = typer.Option("bakta", "--tool", help=HELP_ANNOTATOR),
    all_genomes: bool = typer.Option(False, "--all-genomes", help=HELP_ALL_GENOMES),
    bakta_db: str | None = typer.Option(None, "--bakta-db", help=HELP_BAKTA_DB),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Annotator tuning as key=value (repeatable)."
    ),
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help=HELP_ALLOW_INCOMPLETE,
    ),
    num_processes: int = typer.Option(
        0, "-p", "--num-processes", min=0, help=HELP_ANNOTATE_PROCESSES
    ),
) -> None:
    """Annotate the genome set once (GFF3 per genome, cached) for the core-gene route.

    On this command --tool selects the annotator.
    """
    from ..annotators.base import registry as _annot_registry
    from ..stages.annotate import AnnotateStageParams

    def build() -> AnnotateStageParams:
        _require_choice(tool, set(_annot_registry.names()), "--tool")
        return AnnotateStageParams(
            tool=tool,
            threads=threads,
            all_genomes=all_genomes,
            bakta_db=bakta_db,
            allow_incomplete=allow_incomplete,
            extra=_parse_key_values(tool_arg, "--tool-arg"),
            num_processes=num_processes,
        )

    _run("annotate", workdir, build)


@app.command(rich_help_panel=PANEL_CORE)
def snptype(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    tool: str = typer.Option("simple", "--tool", help=_snp_help()),
    reference: str | None = typer.Option(None, "--reference", help=HELP_REFERENCE),
    all_genomes: bool = typer.Option(False, "--all-genomes", help=HELP_ALL_GENOMES),
    mask: str = typer.Option("none", "--mask", help=_mask_help()),
    core_fraction: float | None = typer.Option(
        None, "--core-fraction", min=0.0, max=1.0, help=HELP_CORE_FRACTION
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    tool_arg: list[str] = typer.Option(
        [], "--tool-arg", help="Tool tuning as key=value (repeatable)."
    ),
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help=HELP_ALLOW_INCOMPLETE,
    ),
) -> None:
    """Call SNPs and build a core-SNP alignment.

    On this command --tool selects the SNP typer, whereas on the dereplication
    commands it selects the dereplicator.
    """
    from ..snptypers.base import registry as _snp_registry
    from ..stages.snptype import SnptypeParams
    from .param_builders import require_mask

    def build() -> SnptypeParams:
        _require_choice(tool, set(_snp_registry.names()), "--tool")
        require_mask(mask)
        return SnptypeParams(
            tool=tool,
            threads=threads,
            reference=reference,
            all_genomes=all_genomes,
            mask=mask,
            core_fraction=core_fraction,
            allow_incomplete=allow_incomplete,
            extra=_parse_key_values(tool_arg, "--tool-arg"),
        )

    _run("snptype", workdir, build)


@app.command(rich_help_panel=PANEL_CORE)
def phylo(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help=HELP_MSA_SOURCE),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help=_snp_help()),
    all_genomes: bool = typer.Option(False, "--all-genomes", help=HELP_ALL_GENOMES),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help=HELP_NO_OUTGROUP),
    root: str | None = typer.Option(None, "--root", help=HELP_ROOT),
    bootstrap: int = typer.Option(0, "-B", "--bootstrap", min=0, help=HELP_BOOTSTRAP),
    reference: str | None = typer.Option(None, "--reference", help=HELP_REFERENCE),
    aligner_arg: list[str] = typer.Option(
        [],
        "--aligner-arg",
        help=HELP_ALIGNER_ARG,
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    mask: str = typer.Option(
        "none",
        "--mask",
        help=_mask_help_msa(),
    ),
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
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help=HELP_ALLOW_INCOMPLETE,
    ),
) -> None:
    """Build a phylogenetic tree from an alignment, SNP alignment, or directly."""
    from .param_builders import gene_selection_extra, phylo_params

    def build():
        return phylo_params(
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
                **gene_selection_extra(
                    max_core_genes=max_core_genes,
                    gene_rank=gene_rank,
                    min_gene_length=min_gene_length,
                    max_gene_gaps=max_gene_gaps,
                ),
                **({"partitioned": True} if partitioned else {}),
            },
            allow_incomplete=allow_incomplete,
        )

    _run("phylo", workdir, build)


@app.command(name="phylo-advise", rich_help_panel=PANEL_PIPELINE)
def phylo_advise(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    all_genomes: bool = typer.Option(False, "--all-genomes", help=HELP_ALL_GENOMES),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    as_json: bool = typer.Option(False, "--json", help=HELP_ADVISE_JSON),
) -> None:
    """Suggest a phylo route and tree builder for the genome set of a workdir.

    Reads the genome names, the existing genome sketches (sourmash compare
    for the pairwise ANI), the annotations and the outgroup, applies the
    table in docs/choosing-tools.md and prints the reasons and the full
    phylo command. Writes nothing to the workdir, except that opening a
    manifest of an older layout migrates its schema in place, as any
    command does. A -wd that does not exist exits 3.
    """
    import logging

    from ..core.context import WorkdirContext
    from ..core.logging import configure_logging
    from ..stages.advise import advise, collect_evidence, format_text, to_json

    require_existing_workdir(workdir)
    level = _RUN_STATE["log_level"]
    # The advice is the output; INFO lines from the sketch lookup are not.
    logger = configure_logging(None, level=logging.WARNING if level == logging.INFO else level)
    with stage_errors(logger):
        ctx = WorkdirContext(workdir, logger=logger)
        try:
            evidence = collect_evidence(ctx, all_genomes, threads, logger)
        finally:
            ctx.close()
    rec = advise(evidence)
    if as_json:
        _echo_json(to_json(rec, evidence))
    else:
        typer.echo(format_text(rec, evidence))


@app.command(rich_help_panel=PANEL_CORE)
def tree2tax(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    node_basename: str | None = typer.Option(None, "--node-basename", help=HELP_NODE_BASENAME),
    naming: str = typer.Option("leafset", "--naming", help=HELP_NAMING),
    root_name: str = typer.Option("root", "--root-name", help=HELP_ROOT_NAME),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help=HELP_REMOVE_OUTGROUP),
    include_dereplicated: bool = typer.Option(
        True,
        "--include-dereplicated/--no-include-dereplicated",
        help=HELP_INCLUDE_DEREPLICATED,
    ),
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
    """Emit FlexTaxD-compatible taxonomy relations from the tree."""
    from .param_builders import tree2tax_params

    def build():
        return tree2tax_params(
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

    _run("tree2tax", workdir, build)
