"""One-shot pipeline orchestrator: ``repgenr run``.

Chains the canonical stages so a user need not invoke five commands by hand. It
forwards the common options and relies on stage defaults for the rest; for full
per-stage control use the individual commands. Each stage goes through the same
:func:`_run` harness, so the resume guard applies -- re-running ``run`` skips
stages already completed with the same parameters.
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from ..core.errors import UserInputError
from ..core.logging import configure_logging
from .base import (
    _RUN_STATE,
    DEFAULT_THREADS,
    HELP_ALIGNED_FRACTION,
    HELP_ALIGNER_ARG,
    HELP_ALL_GENOMES,
    HELP_ALLOW_INCOMPLETE,
    HELP_BOOTSTRAP,
    HELP_COLLAPSE_DIAMETER,
    HELP_COLLAPSE_LENGTH,
    HELP_COLLAPSE_SUPPORT,
    HELP_CORE_FRACTION,
    HELP_DATASET,
    HELP_DEREP_TOOL_ARG,
    HELP_GTDB_RELEASE,
    HELP_GTDB_VERSION,
    HELP_INCLUDE_DEREPLICATED,
    HELP_KEEPER,
    HELP_LEVEL,
    HELP_LIMIT,
    HELP_METADATA_PATH,
    HELP_MSA_SOURCE,
    HELP_NAMING,
    HELP_NO_OUTGROUP,
    HELP_NODE_BASENAME,
    HELP_NODOWNLOAD,
    HELP_NUM_PROCESSES,
    HELP_OUTGROUP_ACCESSION,
    HELP_PRE_PRIMARY_ANI,
    HELP_PRE_SECONDARY_ANI,
    HELP_PRIMARY_ANI,
    HELP_PROCESS_SIZE,
    HELP_RANK_RED,
    HELP_RANKS,
    HELP_REDUCE,
    HELP_REFERENCE,
    HELP_REMOVE_OUTGROUP,
    HELP_ROOT,
    HELP_ROOT_NAME,
    HELP_SECONDARY_ANI,
    HELP_TARGET_FAMILY,
    HELP_TARGET_GENUS,
    HELP_TARGET_REPS,
    HELP_TARGET_SPECIES,
    HELP_THREADS,
    HELP_WORKDIR_CREATED,
    PANEL_PIPELINE,
    PIPELINE_BACTERIAL,
    PIPELINE_LOCAL,
    PIPELINE_READS,
    PIPELINE_VIRAL,
    _aligner_help,
    _assembler_help,
    _derep_help,
    _mask_help,
    _parse_key_values,
    _require_choice,
    _run,
    _snp_help,
    _tree_help,
    app,
    gated_extra,
    resolve_threads,
    stage_errors,
)
from .cmd_viral import _validate_released_after


def _virus_extra(derep_tool: str, viral: bool) -> dict:
    """Virus-tuned extras, injected only when the tool actually reads them.

    Passing extra["virus"] to a tool that ignores it would silently change the
    resume fingerprint for nothing (and imply tuning that never happened).
    """
    if not viral:
        return {}
    from ..dereplicators.base import registry as _derep_registry

    return gated_extra(_derep_registry, derep_tool, "virus", True)


def _preflight_tools(
    derep_tool: str,
    treebuilder: str,
    msa_source: str,
    aligner: str,
    snptyper: str,
    mask: str = "none",
    with_snptype: bool = False,
    assembler: str | None = None,
) -> None:
    """Check every external tool the chain will need before the first stage runs.

    Each stage preflights its own adapter, but by then the earlier stages have
    already downloaded and dereplicated; a missing tree builder should fail in
    the first second, not after the genome download. ``auto`` choices are
    resolved from the genome count inside their stage and are skipped here.
    """
    from ..dereplicators.base import registry as derep_registry
    from ..treebuilders.base import InputKind
    from ..treebuilders.base import registry as tb_registry

    if assembler is not None and assembler != "auto":
        from ..assemblers.base import registry as asm_registry

        asm_registry.create(assembler).preflight()
    if derep_tool != "auto":
        derep_registry.create(derep_tool).preflight()
    if with_snptype:
        from ..snptypers.base import registry as snp_registry

        snp_registry.create(snptyper).preflight()
        if mask not in ("none", ""):
            from ..maskers.base import registry as masker_registry

            masker_registry.create(mask).preflight()
    if treebuilder == "auto":
        return
    builder = tb_registry.create(treebuilder)
    builder.preflight()
    if builder.input_kind != InputKind.MSA_FASTA:
        return
    if msa_source == "aligner":
        from ..aligners.base import registry as aligner_registry

        aligner_registry.create(aligner).preflight()
    elif msa_source == "snptype":
        from ..snptypers.base import registry as snp_registry

        snp_registry.create(snptyper).preflight()
        if mask not in ("none", ""):
            from ..maskers.base import registry as masker_registry

            masker_registry.create(mask).preflight()


def _check_rank_request(
    ranks: bool, rank_red: str | None, no_outgroup: bool, root: str | None
) -> None:
    """Refuse a --rank-red request that tree2tax would refuse after every
    other stage has run: malformed bands, bands without --ranks, or a tree
    that will have no root to measure RED from."""
    from ..stages.tree2tax import _parse_rank_red
    from ..tree.rooting import ROOT_MIDPOINT

    bands = _parse_rank_red(rank_red)
    if not bands:
        return
    if not ranks:
        raise UserInputError("--rank-red needs --ranks")
    if no_outgroup and root != ROOT_MIDPOINT:
        raise UserInputError(
            "--rank-red needs a rooted tree: with --no-outgroup add --root midpoint."
        )


def _warn_collapse_without_supports(
    treebuilder: str,
    collapse_support: float | None,
    collapse_length: float | None,
    logger: logging.Logger,
) -> None:
    """Say before any stage runs that --collapse-support cannot act on the tree
    of a builder that declares the bootstrap ignored (mashtree, sourmash)."""
    from ..treebuilders.base import registry as tb_registry

    if collapse_support is None or treebuilder not in tb_registry.names():
        return  # an unknown name is reported by the parameter check that follows
    if "bootstrap" in tb_registry.get(treebuilder).capabilities.ignored_params:
        logger.warning(
            "--collapse-support %g is given, but the %s tree carries no support values, "
            "so tree2tax will collapse nothing on support%s",
            collapse_support,
            treebuilder,
            ""
            if collapse_length is not None
            else "; use --collapse-length, or a tree builder that bootstraps",
        )


def _msa_source_summary(treebuilder: str, msa_source: str, aligner: str, snptyper: str) -> str:
    """The MSA source a dry run will use, or nothing for an alignment-free builder."""
    from ..treebuilders.base import InputKind
    from ..treebuilders.base import registry as tb_registry

    if treebuilder != "auto" and tb_registry.create(treebuilder).input_kind == InputKind.GENOMES:
        return " (alignment-free)"
    if msa_source == "snptype":
        return f", snptyper={snptyper}"
    return f", aligner={aligner}"


@app.command(rich_help_panel=PANEL_PIPELINE)
def run(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR_CREATED),
    viral: bool = typer.Option(
        False, "--viral", help="Run the viral chain (vmetadata -> vgenome) instead of bacterial."
    ),
    # --- selection: local genomes (ingest) ---
    genomes_dir: Path | None = typer.Option(
        None,
        "--genomes-dir",
        help="Start from local genome FASTAs (the ingest chain) instead of downloading.",
    ),
    from_workdir: list[Path] = typer.Option(
        [],
        "--from-workdir",
        help="Start the ingest chain from the genome set of an earlier working directory "
        "(selection.tsv and genomes/); repeatable, and may be combined with --genomes-dir.",
    ),
    selection: Path | None = typer.Option(
        None,
        "--selection",
        help="With --genomes-dir: selection.tsv naming the genomes to take "
        "(accession, taxonomy, filename, outgroup flag, quality).",
    ),
    outgroup: str | None = typer.Option(
        None,
        "--outgroup",
        help="With --genomes-dir or --from-workdir: the outgroup genome, a name found in "
        "those sources or a path to a FASTA file.",
    ),
    copy: bool = typer.Option(
        False,
        "--copy",
        help="With --genomes-dir or --from-workdir: copy the files into genomes/ instead "
        "of linking.",
    ),
    # --- selection: sequencing reads (ENA/SRA) ---
    reads: bool = typer.Option(
        False,
        "--reads",
        help="Run the reads chain (reads -> assemble) from ENA/SRA sequencing runs "
        "selected by -tf/-tg/-ts or --accession-file, instead of downloading assemblies.",
    ),
    accession_file: Path | None = typer.Option(
        None, "--accession-file", help="With --reads: file of run/sample/study accessions."
    ),
    platform: str = typer.Option(
        "any", "--platform", help="With --reads: any, illumina, ont or pacbio."
    ),
    max_runs: int | None = typer.Option(
        None, "--max-runs", min=1, help="With --reads: keep at most N runs, the largest by bases."
    ),
    assembler: str = typer.Option("auto", "--assembler", help=_assembler_help()),
    # --- selection: bacterial (GTDB) ---
    dataset: str = typer.Option("rep", "-d", "--dataset", help=HELP_DATASET),
    level: str | None = typer.Option(None, "-l", "--level", help=HELP_LEVEL),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help=HELP_TARGET_FAMILY
    ),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus", help=HELP_TARGET_GENUS),
    target_species: str | None = typer.Option(
        None, "-ts", "--target-species", help=HELP_TARGET_SPECIES
    ),
    release: str | None = typer.Option(None, "-r", "--release", help=HELP_GTDB_RELEASE),
    gtdb_version: str | None = typer.Option(None, "--gtdb-version", help=HELP_GTDB_VERSION),
    metadata_source: str = typer.Option("tsv", "--metadata-source", help="tsv or api."),
    metadata_path: str | None = typer.Option(None, "--metadata-path", help=HELP_METADATA_PATH),
    nodownload: bool = typer.Option(False, "--nodownload", help=HELP_NODOWNLOAD),
    outgroup_accession: str | None = typer.Option(
        None, "--outgroup-accession", help=HELP_OUTGROUP_ACCESSION
    ),
    limit: int | None = typer.Option(
        None,
        "--limit",
        min=1,
        help=HELP_LIMIT,
    ),
    # --- selection: viral (NCBI Virus) ---
    target: str | None = typer.Option(None, "--target", help="Virus taxon (viral)."),
    viral_source: str = typer.Option("ncbi_virus", "--viral-source", help="ncbi_virus or bvbrc."),
    complete_only: bool = typer.Option(
        False, "--complete-only", help="ncbi_virus: only COMPLETE sequences (viral)."
    ),
    host: str | None = typer.Option(
        None, "--host", help="ncbi_virus: restrict to a host species (viral)."
    ),
    released_after: str | None = typer.Option(
        None,
        "--released-after",
        callback=_validate_released_after,
        help="ncbi_virus: MM/DD/YYYY (viral).",
    ),
    group_segments: bool = typer.Option(False, "--group-segments", help="Group viral segments."),
    # --- genome download ---
    keep_files: bool = typer.Option(
        False, "--keep-files", help="Keep the genome-download intermediates (bacterial chain)."
    ),
    # --- dereplication ---
    derep_tool: str = typer.Option("skder", "--tool", help=_derep_help()),
    primary_ani: float = typer.Option(0.90, "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(0.50, "--aligned-fraction", help=HELP_ALIGNED_FRACTION),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help=HELP_KEEPER,
    ),
    process_size: int | None = typer.Option(None, "-s", "--process-size", help=HELP_PROCESS_SIZE),
    num_processes: int = typer.Option(0, "-p", "--num-processes", help=HELP_NUM_PROCESSES),
    pre_primary_ani: float | None = typer.Option(
        None,
        "--pre-primary-ani",
        help=HELP_PRE_PRIMARY_ANI,
    ),
    pre_secondary_ani: float | None = typer.Option(
        None,
        "--pre-secondary-ani",
        help=HELP_PRE_SECONDARY_ANI,
    ),
    reduce: str = typer.Option("none", "--reduce", help=HELP_REDUCE),
    target_reps: int = typer.Option(0, "--target-reps", help=HELP_TARGET_REPS),
    tool_arg: list[str] = typer.Option([], "--tool-arg", help=HELP_DEREP_TOOL_ARG),
    # --- SNP typing ---
    with_snptype: bool = typer.Option(
        False,
        "--with-snptype",
        help="Run the standalone snptype stage (with --snptyper, --mask, --reference) after "
        "dereplication, so the SNP tables are produced whatever builds the tree.",
    ),
    # --- phylogeny ---
    treebuilder: str = typer.Option("iqtree", "--treebuilder", help=_tree_help()),
    msa_source: str = typer.Option("aligner", "--msa-source", help=HELP_MSA_SOURCE),
    aligner: str = typer.Option("progressivemauve", "--aligner", help=_aligner_help()),
    snptyper: str = typer.Option("simple", "--snptyper", help=_snp_help()),
    no_outgroup: bool = typer.Option(False, "--no-outgroup", help=HELP_NO_OUTGROUP),
    root: str | None = typer.Option(None, "--root", help=HELP_ROOT),
    all_genomes: bool = typer.Option(False, "--all-genomes", help=HELP_ALL_GENOMES),
    bootstrap: int = typer.Option(0, "-B", "--bootstrap", min=0, help=HELP_BOOTSTRAP),
    reference: str | None = typer.Option(None, "--reference", help=HELP_REFERENCE),
    aligner_arg: list[str] = typer.Option([], "--aligner-arg", help=HELP_ALIGNER_ARG),
    mask: str = typer.Option(
        "none", "--mask", help=_mask_help() + " Needs --msa-source snptype or --with-snptype."
    ),
    # --- taxonomy output ---
    node_basename: str | None = typer.Option(
        None,
        "--node-basename",
        help=HELP_NODE_BASENAME,
    ),
    naming: str = typer.Option("leafset", "--naming", help=HELP_NAMING),
    root_name: str = typer.Option("root", "--root-name", help=HELP_ROOT_NAME),
    remove_outgroup: bool = typer.Option(False, "--remove-outgroup", help=HELP_REMOVE_OUTGROUP),
    include_dereplicated: bool = typer.Option(
        True,
        "--include-dereplicated/--no-include-dereplicated",
        help=HELP_INCLUDE_DEREPLICATED,
    ),
    core_fraction: float | None = typer.Option(
        None, "--core-fraction", min=0.0, max=1.0, help=HELP_CORE_FRACTION
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
    # --- common ---
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help=HELP_ALLOW_INCOMPLETE,
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run", help="Print the stages and key parameters, then exit."
    ),
) -> None:
    """Run the whole pipeline end to end (bacterial by default, --viral for viruses).

    Some option names differ from the single commands. Here --tool selects the
    dereplicator, while 'snptype --tool' selects the SNP typer (--snptyper on
    phylo, run and phylo-build). --metadata-source is 'metadata --source' and --viral-source
    is 'vmetadata --source'. --platform filters the sequencing platform with
    --reads, while the global --platform sets the container platform.
    """
    from .param_builders import (
        METADATA_DATASETS,
        METADATA_LEVELS,
        METADATA_SOURCES,
        VIRAL_SOURCES,
        assemble_params,
        dereplicate_params,
        genome_params,
        ingest_params,
        metadata_params,
        phylo_params,
        reads_params,
        require_mask,
        tree2tax_params,
        vgenome_params,
        vmetadata_params,
    )

    logger = configure_logging(workdir if workdir.exists() else None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        # Fail fast before any stage runs; the same validation re-runs inside
        # each stage's builder. The selection sources are checked here under
        # the names `run` gives them (the builders know them as --source).
        _require_choice(metadata_source, METADATA_SOURCES, "--metadata-source")
        _require_choice(viral_source, VIRAL_SOURCES, "--viral-source")
        _require_choice(dataset, METADATA_DATASETS, "--dataset")
        tree2tax_params(naming=naming)
        if level is not None:
            _require_choice(level, METADATA_LEVELS, "--level")
        derep_extra = {
            **_parse_key_values(tool_arg, "--tool-arg"),
            **_virus_extra(derep_tool, viral),
        }
        _warn_collapse_without_supports(treebuilder, collapse_support, collapse_length, logger)
        _check_rank_request(ranks, rank_red, no_outgroup, root)
        # With --with-snptype and the aligner MSA source, --mask belongs to the
        # standalone snptype stage only; phylo would reject it there.
        mask_for_phylo = mask != "none" and (msa_source == "snptype" or not with_snptype)
        # The same holds for --core-fraction, which the typing pass applies.
        core_for_phylo = core_fraction is not None and (msa_source == "snptype" or not with_snptype)
        phylo_extra = {
            **_parse_key_values(aligner_arg, "--aligner-arg"),
            **({"mask": mask} if mask_for_phylo else {}),
            **({"core_fraction": core_fraction} if core_for_phylo else {}),
        }
        dereplicate_params(
            tool=derep_tool,
            primary_ani=primary_ani,
            secondary_ani=secondary_ani,
            aligned_fraction=aligned_fraction,
            pre_primary_ani=pre_primary_ani,
            pre_secondary_ani=pre_secondary_ani,
            keeper=keeper,
            reduce=reduce,
            target_reps=target_reps,
            extra=derep_extra,
        )
        phylo_params(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            no_outgroup=no_outgroup,
            root=root,
            extra=phylo_extra,
            on_run=True,
        )
        local = genomes_dir is not None or bool(from_workdir)
        if local and viral:
            raise UserInputError(
                "--genomes-dir and --from-workdir start the local chain (ingest); they cannot "
                "be combined with --viral."
            )
        if reads and (viral or local):
            raise UserInputError(
                "--reads starts the reads chain; it cannot be combined with --viral, "
                "--genomes-dir or --from-workdir."
            )
        if reads:
            # Also parses --accession-file and checks that a selection is given.
            reads_params(
                target_family=target_family,
                target_genus=target_genus,
                target_species=target_species,
                accession_file=None if accession_file is None else str(accession_file),
                platform=platform,
            )
            assemble_params(assembler=assembler)
        if with_snptype:
            # The standalone snptype stage uses --snptyper whatever the MSA
            # source; check the name here so a typo fails before any stage.
            from ..snptypers.base import registry as _snp_registry

            _require_choice(snptyper, set(_snp_registry.names()), "--snptyper")
            require_mask(mask)
        if not viral and not local and not reads and not level:
            raise UserInputError("The bacterial chain needs -l/--level (family/genus/species).")

    if dry_run:
        chain: tuple[str, ...] = (
            PIPELINE_READS
            if reads
            else PIPELINE_LOCAL
            if local
            else PIPELINE_VIRAL
            if viral
            else PIPELINE_BACTERIAL
        )
        if with_snptype:
            i = chain.index("phylo")
            chain = (*chain[:i], "snptype", *chain[i:])
        lineage = "reads" if reads else "local" if local else "viral" if viral else "bacterial"
        typer.echo(f"[dry-run] {lineage} pipeline in {workdir}:")
        for stage in chain:
            typer.echo(f"  - {stage}")
        selection_summary = (
            f"family={target_family}, genus={target_genus}, species={target_species}, "
            f"accession_file={accession_file}, platform={platform}, max_runs={max_runs}, "
            f"assembler={assembler}"
            if reads
            else f"genomes_dir={genomes_dir}, from_workdir={[str(w) for w in from_workdir]}, "
            f"selection={selection}, outgroup={outgroup}"
            if local
            else f"target={target}, genus={target_genus}, species={target_species}"
            if viral
            else f"dataset={dataset}, level={level}, "
            f"family={target_family}, genus={target_genus}, species={target_species}"
        )
        typer.echo(f"selection: {selection_summary}")
        typer.echo(
            f"dereplicate: tool={derep_tool}, primary_ani={primary_ani}, "
            f"secondary_ani={secondary_ani}; phylo: treebuilder={treebuilder}"
            f"{_msa_source_summary(treebuilder, msa_source, aligner, snptyper)}; "
            f"threads={threads}"
        )
        typer.echo("[dry-run] no work done.")
        return

    with stage_errors(logger):
        _preflight_tools(
            derep_tool,
            treebuilder,
            msa_source,
            aligner,
            snptyper,
            mask,
            with_snptype,
            assembler if reads else None,
        )

    if reads:
        _run(
            "reads",
            workdir,
            lambda: reads_params(
                target_family=target_family,
                target_genus=target_genus,
                target_species=target_species,
                accession_file=None if accession_file is None else str(accession_file),
                platform=platform,
                max_runs=max_runs,
            ),
            create=True,
        )
        _run(
            "assemble",
            workdir,
            lambda: assemble_params(assembler=assembler, threads=threads, outgroup=outgroup),
        )
    elif local:
        _run(
            "ingest",
            workdir,
            lambda: ingest_params(
                genomes_dir=None if genomes_dir is None else str(genomes_dir),
                from_workdirs=[str(wd) for wd in from_workdir],
                selection=None if selection is None else str(selection),
                outgroup=outgroup,
                copy=copy,
            ),
            create=True,
        )
    elif viral:
        _run(
            "vmetadata",
            workdir,
            lambda: vmetadata_params(
                target=target,
                source=viral_source,
                complete_only=complete_only,
                host=host,
                released_after=released_after,
            ),
            create=True,
        )
        _run(
            "vgenome",
            workdir,
            lambda: vgenome_params(
                target_genus=target_genus,
                target_species=target_species,
                no_outgroup=no_outgroup,
                outgroup_accession=outgroup_accession,
                group_segments=group_segments,
            ),
        )
    else:
        _run(
            "metadata",
            workdir,
            lambda: metadata_params(
                dataset=dataset,
                level=level or "",
                source=metadata_source,
                release=release,
                version=gtdb_version,
                target_family=target_family,
                target_genus=target_genus,
                target_species=target_species,
                outgroup_accession=outgroup_accession,
                metadata_path=metadata_path,
                nodownload=nodownload,
                limit=limit,
            ),
            create=True,
        )
        _run("genome", workdir, lambda: genome_params(keep_files=keep_files))

    _run(
        "dereplicate",
        workdir,
        lambda: dereplicate_params(
            tool=derep_tool,
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
            extra=derep_extra,
            keeper=keeper,
            allow_incomplete=allow_incomplete,
        ),
    )

    # phylo's own typing pass (--msa-source snptype) writes under tree/msa/,
    # so the standalone stage's snp/ tables are left as its record describes.
    if with_snptype:
        from ..stages.snptype import SnptypeParams

        _run(
            "snptype",
            workdir,
            lambda: SnptypeParams(
                tool=snptyper,
                threads=threads,
                reference=reference,
                all_genomes=all_genomes,
                mask=mask,
                core_fraction=core_fraction,
                allow_incomplete=allow_incomplete,
                extra=_parse_key_values(aligner_arg, "--aligner-arg"),
            ),
        )
    _run(
        "phylo",
        workdir,
        lambda: phylo_params(
            treebuilder=treebuilder,
            msa_source=msa_source,
            aligner=aligner,
            snptyper=snptyper,
            no_outgroup=no_outgroup,
            root=root,
            all_genomes=all_genomes,
            bootstrap=bootstrap,
            reference=reference,
            threads=threads,
            extra=phylo_extra,
            allow_incomplete=allow_incomplete,
        ),
    )
    _run(
        "tree2tax",
        workdir,
        lambda: tree2tax_params(
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
        ),
    )

    typer.echo(
        f"\nPipeline complete. Deliverables in {workdir}: "
        "tree2tax.tsv, genomes_map.tsv, tree/tree.nwk"
    )
