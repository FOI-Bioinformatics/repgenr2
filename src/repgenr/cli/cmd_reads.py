"""``repgenr reads``: select sequencing runs from ENA/SRA for the reads chain."""

from __future__ import annotations

from pathlib import Path

import typer

from .base import (
    DEFAULT_THREADS,
    HELP_MAX_SPOTS,
    HELP_READS_SKETCH,
    HELP_READS_TSV_STANDALONE,
    HELP_SCREEN_MIN_FRACTION,
    HELP_SCREEN_READS,
    HELP_SKETCH,
    HELP_SKETCH_JOBS,
    HELP_SKETCH_OUT,
    HELP_THREADS,
    HELP_VERSIONS_OUT,
    HELP_WORKDIR,
    HELP_WORKDIR_CREATED,
    PANEL_ENTRY,
    _assembler_help,
    _classifier_help,
    _parse_key_values,
    _polisher_help,
    _run,
    app,
    resolve_threads,
    stage_errors,
)

_HELP_TAXON = (
    "Restrict the selection to this {}. Only the most specific of -tf/-tg/-ts "
    "given is used; they are not combined."
)


@app.command(rich_help_panel=PANEL_ENTRY)
def reads(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR_CREATED),
    target_family: str | None = typer.Option(
        None, "-tf", "--target-family", help=_HELP_TAXON.format("family")
    ),
    target_genus: str | None = typer.Option(
        None, "-tg", "--target-genus", help=_HELP_TAXON.format("genus")
    ),
    target_species: str | None = typer.Option(
        None, "-ts", "--target-species", help=_HELP_TAXON.format("species")
    ),
    accession: list[str] = typer.Option(
        [],
        "--accession",
        help="A run (SRR/ERR/DRR), sample (SAMN.., SRS..) or study (PRJNA.., SRP..) "
        "accession to include (repeatable).",
    ),
    accession_file: Path | None = typer.Option(
        None, "--accession-file", help="File of accessions, one per line (# comments allowed)."
    ),
    platform: str = typer.Option(
        "any", "--platform", help="Keep runs of one platform: any, illumina, ont or pacbio."
    ),
    max_runs: int | None = typer.Option(
        None, "--max-runs", min=1, help="Keep at most N runs, the largest by bases."
    ),
    min_bases: int = typer.Option(
        0, "--min-bases", min=0, help="Drop runs with fewer sequenced bases than this."
    ),
    max_bases: int | None = typer.Option(
        None,
        "--max-bases",
        min=1,
        help="Drop runs with more sequenced bases than this (a guard against whole-host "
        "libraries, which would assemble into a host-dominated genome).",
    ),
    drop_selection: list[str] = typer.Option(
        ["MDA"],
        "--drop-selection",
        help="Drop runs whose ENA library selection is this value (repeatable; MDA is "
        "whole-genome amplification). Pass 'none' to keep every selection.",
    ),
    one_per_sample: bool = typer.Option(
        True,
        "--one-per-sample/--all-runs",
        help="Keep the best run of each sample (a long-read run with enough bases, else "
        "the largest run), or every run.",
    ),
) -> None:
    """Select sequencing runs from ENA/SRA by taxon or accession (writes reads.tsv)."""
    from .param_builders import reads_params

    def build():
        return reads_params(
            target_family=target_family,
            target_genus=target_genus,
            target_species=target_species,
            accessions=list(accession),
            accession_file=None if accession_file is None else str(accession_file),
            platform=platform,
            max_runs=max_runs,
            min_bases=min_bases,
            max_bases=max_bases,
            drop_selection=[s for s in drop_selection if s.lower() != "none"],
            one_per_sample=one_per_sample,
        )

    _run("reads", workdir, build, create=True)


@app.command(name="reads-sketch", rich_help_panel=PANEL_ENTRY)
def reads_sketch(
    workdir: Path | None = typer.Option(
        None,
        "-wd",
        "--workdir",
        help="Workdir whose reads stage selected the runs (reads.tsv); the sketches go to "
        "reads_sketches/<variant>/ (all/ or spots_N/) and the stage is recorded.",
    ),
    accession: list[str] = typer.Option(
        [],
        "--accession",
        help="Restrict to these runs of reads.tsv (repeatable); a run not in reads.tsv is "
        "refused. Without -wd: the runs to sketch.",
    ),
    accession_file: Path | None = typer.Option(
        None, "--accession-file", help="File of run accessions, one per line (# comments allowed)."
    ),
    platform: str = typer.Option(
        "any", "--platform", help="Restrict to runs of one platform: any, illumina, ont or pacbio."
    ),
    max_spots: int | None = typer.Option(None, "--max-spots", min=1, help=HELP_MAX_SPOTS),
    jobs: int | None = typer.Option(None, "--jobs", min=1, help=HELP_SKETCH_JOBS),
    out: Path | None = typer.Option(None, "-o", "--out", help=HELP_SKETCH_OUT),
    reads_tsv: Path | None = typer.Option(None, "--reads-tsv", help=HELP_READS_TSV_STANDALONE),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
    allow_all_failed: bool = typer.Option(
        False,
        "--allow-all-failed",
        help="Standalone form: a call where every run failed to stream is a warning, not an "
        "error (the call is one chunk of a batch, as in the Nextflow READS_SKETCH task).",
    ),
) -> None:
    """Sketch the reads of selected runs streamed from SRA, without downloading them."""
    from ..core.logging import configure_logging
    from ..stages.reads_sketch import ReadsSketchParams, run_standalone, selected_accessions
    from .base import _RUN_STATE

    if workdir is None and out is None:
        raise typer.BadParameter("give -wd WORKDIR, or -o DIR for the standalone form")
    if workdir is not None and out is not None:
        raise typer.BadParameter("-o is the standalone form; it does not combine with -wd")
    if workdir is not None and (reads_tsv is not None or versions_out is not None):
        raise typer.BadParameter("--reads-tsv and --versions-out belong to the standalone form")
    if workdir is not None and allow_all_failed:
        raise typer.BadParameter("--allow-all-failed belongs to the standalone form")
    if workdir is None and platform != "any":
        raise typer.BadParameter("--platform selects from reads.tsv; it needs -wd")

    def build():
        return ReadsSketchParams(
            accessions=list(accession),
            accession_file=None if accession_file is None else str(accession_file),
            platform=platform,
            max_spots=max_spots,
            jobs=jobs,
        )

    if workdir is not None:
        _run("reads_sketch", workdir, build)
        return
    assert out is not None
    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        params = build()
        run_standalone(
            selected_accessions(params),
            out.resolve(),
            max_spots=max_spots,
            jobs=jobs,
            force=bool(_RUN_STATE["force"]),
            logger=logger,
            reads_tsv=None if reads_tsv is None else reads_tsv.resolve(),
            versions_out=versions_out,
            allow_all_failed=allow_all_failed,
        )


@app.command(name="reads-screen", rich_help_panel=PANEL_ENTRY)
def reads_screen(
    workdir: Path | None = typer.Option(
        None,
        "-wd",
        "--workdir",
        help="Workdir with reads.tsv and reads_sketches/; the decisions go to "
        "reads_screen/<variant>/ and the stage is recorded.",
    ),
    accession: list[str] = typer.Option(
        [],
        "--accession",
        help="Restrict to these runs of reads.tsv (repeatable); a run not in reads.tsv is refused.",
    ),
    accession_file: Path | None = typer.Option(
        None, "--accession-file", help="File of run accessions, one per line (# comments allowed)."
    ),
    max_spots: int | None = typer.Option(
        None,
        "--max-spots",
        min=1,
        help="Screen the sketches of reads-sketch --max-spots N (reads_sketches/spots_N/). "
        "Default: all/ when present, else the only variant there.",
    ),
    min_fraction: float = typer.Option(
        0.3, "--min-fraction", min=0.0, max=1.0, help=HELP_SCREEN_MIN_FRACTION
    ),
    dup_ani: float = typer.Option(
        0.999,
        "--dup-ani",
        min=0.0,
        max=1.0,
        help="Minimum ANI estimate at which a run contained in an earlier run of the same "
        "biosample is a duplicate (duplicate_isolate).",
    ),
    gtdb_sketch: str | None = typer.Option(
        None,
        "--gtdb-sketch",
        help="GTDB sourmash sketch database (or set REPGENR_GTDB_SKETCH).",
    ),
    gtdb_lineages: str | None = typer.Option(
        None,
        "--gtdb-lineages",
        help="Lineages CSV published with the GTDB sketch (or set REPGENR_GTDB_LINEAGES).",
    ),
    jobs: int | None = typer.Option(
        None,
        "--jobs",
        min=1,
        help="Runs gathered at once; each gather takes about 0.6 GB against the GTDB "
        "representatives sketch. Default: 4, or the usable CPUs when fewer.",
    ),
    reads_tsv: Path | None = typer.Option(None, "--reads-tsv", help=HELP_READS_TSV_STANDALONE),
    sketches_dir: Path | None = typer.Option(
        None,
        "--sketches-dir",
        help="Standalone form: the -o DIR of reads-sketch, holding <variant>/ directories.",
    ),
    out: Path | None = typer.Option(
        None,
        "-o",
        "--out",
        help="Standalone form: write the decisions to DIR/<variant>/ instead of a workdir.",
    ),
    versions_out: Path | None = typer.Option(None, "--versions-out", help=HELP_VERSIONS_OUT),
) -> None:
    """Judge the sketched runs against GTDB (taxon, match fraction, duplicates) before
    downloading any of them."""
    from ..core.logging import configure_logging
    from ..stages.reads_screen import ReadsScreenParams, run_standalone
    from .base import _RUN_STATE

    standalone = [reads_tsv, sketches_dir, out]
    if workdir is None and any(v is None for v in standalone):
        raise typer.BadParameter(
            "give -wd WORKDIR, or --reads-tsv, --sketches-dir and -o DIR for the standalone form"
        )
    if workdir is not None and (any(v is not None for v in standalone) or versions_out):
        raise typer.BadParameter(
            "--reads-tsv, --sketches-dir, -o and --versions-out belong to the standalone form"
        )

    def build():
        return ReadsScreenParams(
            accessions=list(accession),
            accession_file=None if accession_file is None else str(accession_file),
            max_spots=max_spots,
            min_fraction=min_fraction,
            dup_ani=dup_ani,
            gtdb_sketch=gtdb_sketch,
            gtdb_lineages=gtdb_lineages,
            jobs=jobs,
        )

    if workdir is not None:
        _run("reads_screen", workdir, build)
        return
    assert reads_tsv is not None and sketches_dir is not None and out is not None
    logger = configure_logging(None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        run_standalone(
            reads_tsv.resolve(),
            sketches_dir.resolve(),
            out.resolve(),
            build(),
            force=bool(_RUN_STATE["force"]),
            logger=logger,
            versions_out=versions_out,
        )


@app.command(rich_help_panel=PANEL_ENTRY)
def assemble(
    workdir: Path = typer.Option(..., "-wd", "--workdir", help=HELP_WORKDIR),
    assembler: str = typer.Option("auto", "--assembler", help=_assembler_help()),
    threads: int = typer.Option(
        DEFAULT_THREADS, "-t", "--threads", min=1, help=HELP_THREADS, callback=resolve_threads
    ),
    jobs: int | None = typer.Option(
        None,
        "--jobs",
        min=1,
        help="Runs assembled concurrently, threads split across them. Default: 2, or 1 "
        "when a long-read run is pending.",
    ),
    memory_gb: int = typer.Option(
        16,
        "--memory-gb",
        min=1,
        help="Memory hint per assembly, in GB, for tools that cap RAM; also the budget that "
        "caps concurrent classifier gathers (about 0.6 GB each).",
    ),
    min_contig_length: int = typer.Option(
        500, "--min-contig-length", min=0, help="Drop contigs shorter than this many bases."
    ),
    polisher: str = typer.Option("auto", "--polisher", help=_polisher_help()),
    polish_rounds: int = typer.Option(
        1, "--polish-rounds", min=1, help="Polishing rounds (racon; medaka runs one)."
    ),
    outgroup: Path | None = typer.Option(
        None, "--outgroup", help="A FASTA file to set aside as the outgroup for rooting."
    ),
    accession: list[str] = typer.Option(
        [],
        "--accession",
        help="Assemble only these runs of reads.tsv (repeatable); the other pending runs are "
        "deferred, as under --max-runs. A run not in reads.tsv is refused.",
    ),
    accession_file: Path | None = typer.Option(
        None,
        "--accession-file",
        help="File of run accessions to assemble, one per line (# comments allowed), such as "
        "the passed.txt of reads-screen.",
    ),
    max_runs: int | None = typer.Option(
        None,
        "--max-runs",
        min=1,
        help="Keep at most N runs with a finished assembly, in reads.tsv order (largest "
        "first); finished runs count toward N. The others are deferred in "
        "excused_runs.tsv, and a later call with a larger N assembles them.",
    ),
    append: bool = typer.Option(
        False,
        "--append",
        help="Add the assemblies to a working directory that already holds a selection "
        "(metadata and genome, or ingest) instead of replacing it.",
    ),
    keep_reads: bool = typer.Option(
        False, "--keep-reads", help="Keep the downloaded FASTQ files after assembling."
    ),
    keep_files: bool = typer.Option(
        False, "--keep-files", help="Keep each run's assembler scratch directory."
    ),
    checkm2_db: Path | None = typer.Option(
        None,
        "--checkm2-db",
        help="CheckM2 DIAMOND database; enables quality scoring (or set CHECKM2DB).",
    ),
    min_completeness: float = typer.Option(
        50.0, "--min-completeness", min=0.0, max=100.0, help="CheckM2 completeness floor."
    ),
    max_contamination: float = typer.Option(
        10.0, "--max-contamination", min=0.0, max=100.0, help="CheckM2 contamination ceiling."
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
        [], "--tool-arg", help="Assembler tuning as key=value (repeatable), e.g. mode=nano-raw."
    ),
    sketch: bool | None = typer.Option(
        None,
        "--sketch/--no-sketch",
        help=HELP_SKETCH + " With --append only the new genomes are sketched.",
    ),
    reads_sketch: bool | None = typer.Option(
        None,
        "--reads-sketch/--no-reads-sketch",
        help=HELP_READS_SKETCH,
    ),
    screen_reads: bool = typer.Option(False, "--screen-reads", help=HELP_SCREEN_READS),
    screen_min_fraction: float = typer.Option(
        0.3, "--screen-min-fraction", min=0.0, max=1.0, help=HELP_SCREEN_MIN_FRACTION
    ),
    screen_dup_ani: float = typer.Option(
        0.999,
        "--screen-dup-ani",
        min=0.0,
        max=1.0,
        help="With --screen-reads: ANI estimate from which a run contained in an earlier "
        "accepted run of the same biosample is excused as a duplicate.",
    ),
) -> None:
    """Fetch and assemble the selected runs; write genomes/ and selection.tsv."""
    from .param_builders import assemble_params

    def build():
        return assemble_params(
            assembler=assembler,
            threads=threads,
            jobs=jobs,
            memory_gb=memory_gb,
            min_contig_length=min_contig_length,
            polisher=polisher,
            polish_rounds=polish_rounds,
            outgroup=None if outgroup is None else str(outgroup),
            append=append,
            keep_reads=keep_reads,
            keep_files=keep_files,
            checkm2_db=None if checkm2_db is None else str(checkm2_db),
            min_completeness=min_completeness,
            max_contamination=max_contamination,
            classifier=classifier,
            gtdb_sketch=None if gtdb_sketch is None else str(gtdb_sketch),
            gtdb_lineages=None if gtdb_lineages is None else str(gtdb_lineages),
            sketch=sketch,
            reads_sketch=reads_sketch,
            extra=_parse_key_values(tool_arg, "--tool-arg"),
            max_runs=max_runs,
            accessions=list(accession),
            accession_file=None if accession_file is None else str(accession_file),
            screen_reads=screen_reads,
            screen_min_fraction=screen_min_fraction,
            screen_dup_ani=screen_dup_ani,
        )

    _run("assemble", workdir, build)
