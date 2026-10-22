"""``repgenr ingest``: start a working directory from local genome FASTAs."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import HELP_SKETCH, HELP_WORKDIR_CREATED, PANEL_ENTRY, _run, app


@app.command(rich_help_panel=PANEL_ENTRY)
def ingest(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR_CREATED),
    genomes_dir: Path | None = typer.Option(
        None,
        "--genomes-dir",
        help="Directory of genome FASTA files (.fasta, .fa, .fna, .fas, .fasta.gz, .fna.gz "
        "or .fa.gz; subdirectories are not searched) to stage under genomes/. "
        "At least one of --genomes-dir and --from-workdir is required.",
    ),
    from_workdir: list[Path] = typer.Option(
        [],
        "--from-workdir",
        help="Earlier working directory whose genome set (selection.tsv and genomes/) is "
        "added with each genome's taxonomy, quality, GTDB representative flag and manifest "
        "source; its outgroup is not carried over. Repeatable, and combines with "
        "--genomes-dir.",
    ),
    strict_sources: bool = typer.Option(
        False,
        "--strict-sources",
        help="Refuse (exit 2, nothing written) a --from-workdir whose genome set has no "
        "stage record or one that status lists as interrupted or stale. Without it such a "
        "source is named in a warning and taken.",
    ),
    selection: Path | None = typer.Option(
        None,
        "--selection",
        help="selection.tsv (accession, taxonomy, filename, outgroup flag, quality) "
        "naming the genomes to take from --genomes-dir. Default: every FASTA under "
        "--genomes-dir, taxonomy parsed from canonical Family_genus_species_ACCESSION names.",
    ),
    outgroup: str | None = typer.Option(
        None,
        "--outgroup",
        help="Genome to set aside as the outgroup: a filename, stem or accession under "
        "--genomes-dir or in a --from-workdir (genomes/ or outgroup/), or a path to a "
        "FASTA file elsewhere.",
    ),
    copy: bool = typer.Option(
        False, "--copy", help="Copy the files into genomes/ instead of symlinking them."
    ),
    drop_foreign: bool = typer.Option(
        False,
        "--drop-foreign",
        help="Discard genomes appended from sequencing runs (assemble --append) instead of "
        "refusing to overwrite the selection that holds them.",
    ),
    sketch: bool | None = typer.Option(
        None,
        "--sketch/--no-sketch",
        help=HELP_SKETCH + " A --from-workdir sketch of the same genome is copied.",
    ),
) -> None:
    """Populate a working directory from local genomes (no download)."""
    from .param_builders import ingest_params

    def build():
        return ingest_params(
            genomes_dir=None if genomes_dir is None else str(genomes_dir),
            from_workdirs=[str(wd) for wd in from_workdir],
            selection=None if selection is None else str(selection),
            outgroup=outgroup,
            copy=copy,
            drop_foreign=drop_foreign,
            sketch=sketch,
            strict_sources=strict_sources,
        )

    _run("ingest", workdir, build, create=True)
