"""Bacterial lineage commands: metadata, genome, dereplicate."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import (
    DEFAULT_THREADS,
    HELP_ALIGNED_FRACTION,
    HELP_ALLOW_INCOMPLETE,
    HELP_DATASET,
    HELP_DEREP_TOOL_ARG,
    HELP_GTDB_RELEASE,
    HELP_GTDB_VERSION,
    HELP_KEEP_FILES,
    HELP_KEEPER,
    HELP_LEVEL,
    HELP_LIMIT,
    HELP_METADATA_PATH,
    HELP_NODOWNLOAD,
    HELP_NUM_PROCESSES,
    HELP_OUTGROUP_ACCESSION,
    HELP_PRE_PRIMARY_ANI,
    HELP_PRE_SECONDARY_ANI,
    HELP_PRIMARY_ANI,
    HELP_PROCESS_SIZE,
    HELP_REDUCE,
    HELP_REFINE,
    HELP_REFINE_ANI,
    HELP_SECONDARY_ANI,
    HELP_SKETCH,
    HELP_TARGET_FAMILY,
    HELP_TARGET_GENUS,
    HELP_TARGET_REPS,
    HELP_TARGET_SPECIES,
    HELP_THREADS,
    HELP_VIRUS,
    HELP_WORKDIR,
    HELP_WORKDIR_CREATED,
    PANEL_CORE,
    PANEL_ENTRY,
    _derep_help,
    _parse_key_values,
    _run,
    app,
    gated_extra,
    resolve_threads,
)


@app.command(rich_help_panel=PANEL_ENTRY)
def metadata(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR_CREATED),
    dataset: str = typer.Option(..., "-d", "--dataset", help=HELP_DATASET),
    level: str = typer.Option(..., "-l", "--level", help=HELP_LEVEL),
    source: str = typer.Option(
        "tsv", "--source", help="tsv (download full table) or api (GTDB API, target only)."
    ),
    release: str | None = typer.Option(None, "-r", "--release", help=HELP_GTDB_RELEASE),
    gtdb_version: str | None = typer.Option(None, "--gtdb-version", help=HELP_GTDB_VERSION),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help=HELP_TARGET_FAMILY
    ),
    target_genus: str | None = typer.Option(None, "-tg", "--target-genus", help=HELP_TARGET_GENUS),
    target_species: str | None = typer.Option(
        None, "-ts", "--target-species", help=HELP_TARGET_SPECIES
    ),
    outgroup_accession: str | None = typer.Option(
        None, "--outgroup-accession", help=HELP_OUTGROUP_ACCESSION
    ),
    metadata_path: str | None = typer.Option(None, "--metadata-path", help=HELP_METADATA_PATH),
    nodownload: bool = typer.Option(False, "--nodownload", help=HELP_NODOWNLOAD),
    limit: int | None = typer.Option(
        None,
        "--limit",
        min=1,
        help=HELP_LIMIT,
    ),
    drop_foreign: bool = typer.Option(
        False,
        "--drop-foreign",
        help="Discard genomes appended from sequencing runs (assemble --append) instead of "
        "refusing to overwrite the selection that holds them.",
    ),
) -> None:
    """Select a taxon's genomes from GTDB (full table or the GTDB API).

    The --source option is called --metadata-source in 'run'.
    """
    from .param_builders import metadata_params

    def build():
        return metadata_params(
            dataset=dataset,
            level=level,
            source=source,
            release=release,
            version=gtdb_version,
            target_family=target_family,
            target_genus=target_genus,
            target_species=target_species,
            outgroup_accession=outgroup_accession,
            metadata_path=metadata_path,
            nodownload=nodownload,
            limit=limit,
            drop_foreign=drop_foreign,
        )

    _run("metadata", workdir, build, create=True)


@app.command(rich_help_panel=PANEL_ENTRY)
def genome(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    accession_list_only: bool = typer.Option(
        False, "--accession-list-only", help="Write the accession list and stop (no download)."
    ),
    keep_files: bool = typer.Option(False, "--keep-files", help=HELP_KEEP_FILES),
    sketch: bool | None = typer.Option(None, "--sketch/--no-sketch", help=HELP_SKETCH),
) -> None:
    """Download and organize genomes selected by the metadata stage."""
    from .param_builders import genome_params

    def build():
        return genome_params(
            accession_list_only=accession_list_only, keep_files=keep_files, sketch=sketch
        )

    _run("genome", workdir, build)


@app.command(rich_help_panel=PANEL_CORE)
def dereplicate(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    tool: str = typer.Option("skder", "--tool", help=_derep_help()),
    primary_ani: float = typer.Option(0.90, "-pani", "--primary-ani", help=HELP_PRIMARY_ANI),
    secondary_ani: float = typer.Option(0.99, "-sani", "--secondary-ani", help=HELP_SECONDARY_ANI),
    aligned_fraction: float = typer.Option(
        0.50, "-af", "--aligned-fraction", help=HELP_ALIGNED_FRACTION
    ),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    process_size: int | None = typer.Option(
        None,
        "-s",
        "--process-size",
        help=HELP_PROCESS_SIZE,
    ),
    num_processes: int = typer.Option(
        0,
        "-p",
        "--num-processes",
        help=HELP_NUM_PROCESSES,
    ),
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
    reduce: str = typer.Option(
        "none",
        "--reduce",
        help=HELP_REDUCE,
    ),
    target_reps: int = typer.Option(
        0,
        "--target-reps",
        help=HELP_TARGET_REPS,
    ),
    virus: bool = typer.Option(False, "--virus", help=HELP_VIRUS),
    tool_arg: list[str] = typer.Option([], "--tool-arg", help=HELP_DEREP_TOOL_ARG),
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help=HELP_ALLOW_INCOMPLETE,
    ),
    keeper: str = typer.Option(
        "quality",
        "--keeper",
        help=HELP_KEEPER,
    ),
    refine: list[str] = typer.Option([], "--refine", help=HELP_REFINE),
    refine_ani: float | None = typer.Option(None, "--refine-ani", help=HELP_REFINE_ANI),
) -> None:
    """Cluster genomes by ANI and select representatives."""
    from ..dereplicators.base import registry as _derep_registry
    from .param_builders import dereplicate_params

    def build():
        return dereplicate_params(
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
            extra={
                **_parse_key_values(tool_arg, "--tool-arg"),
                **(
                    gated_extra(_derep_registry, tool, "virus", True, flag="--virus")
                    if virus
                    else {}
                ),
            },
            allow_incomplete=allow_incomplete,
            keeper=keeper,
            refine=refine,
            refine_ani=refine_ani,
        )

    _run("dereplicate", workdir, build)


@app.command(name="derep-refine", rich_help_panel=PANEL_CORE)
def derep_refine(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    refine: list[str] = typer.Option(..., "--refine", help=HELP_REFINE),
    refine_ani: float = typer.Option(..., "--refine-ani", help=HELP_REFINE_ANI),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
) -> None:
    """Keep more genomes inside a lineage: dereplicate the seeds' clusters again at a finer ANI.

    Refines the dereplication on disk with the tool and settings of its
    record, without repeating it. Store the current result first with
    derep-stock if it should be kept.
    """
    from .param_builders import derep_refine_params

    _run(
        "derep_refine",
        workdir,
        lambda: derep_refine_params(seeds=refine, ani=refine_ani, threads=threads),
    )
