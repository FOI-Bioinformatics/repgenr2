"""Shared CLI app, callback and helpers.

The Typer ``app``, the top-level callback (container/logging setup) and the
common stage harness (:func:`_run`) live here so the per-domain command modules
(``cmd_*.py``) can register against a single app without circular imports.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
import os
from collections.abc import Iterator
from collections.abc import Set as AbstractSet
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import typer
from typer.core import TyperGroup

from .. import __version__
from ..core.context import WorkdirContext
from ..core.contracts import (
    ANNOTATIONS_TSV,
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    CORE_SNP_FASTA,
    GENOME_STATUS_TSV,
    GENOMES_MAP_TSV,
    OUTGROUP_ACCESSION_TXT,
    READS_TSV,
    REFINED_TSV,
    SEGMENTS_TSV,
    SELECTION_TSV,
    TREE2TAX_NODES_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    list_fasta,
    read_clusters,
)
from ..core.errors import (
    MissingBinaryError,
    RepGenRError,
    ToolExecutionError,
    UserInputError,
    WorkdirError,
)
from ..core.inputs import dir_stat_digest, inputs_digest, manifest_digest_for_stage
from ..core.logging import configure_logging
from ..core.manifest import MANIFEST_FILENAME

# Top-level run options shared by every subcommand (set in the callback).
_RUN_STATE: dict[str, Any] = {"force": False, "log_level": logging.INFO}

# One default thread count for every stage's -t/--threads, so the CLI is
# consistent (stages previously mixed 16 and 24).
DEFAULT_THREADS = 16

# Help texts shared by several commands, so the same flag reads the same
# everywhere and the CLI matrix (tests/audit/cli_matrix.yaml) can hold one
# sentence per flag.
HELP_THREADS = "Threads for the external tool (default 16, or the CPU limit when lower)."
HELP_PRIMARY_ANI = "Primary (pre-clustering) ANI threshold in (0, 1]."
HELP_SECONDARY_ANI = "Secondary (final cluster) ANI threshold in (0, 1]."
HELP_ALIGNED_FRACTION = "Minimum aligned fraction in (0, 1] for a pair to be compared."
HELP_TARGET_FAMILY = "Restrict the selection to this family."
HELP_TARGET_GENUS = "Restrict the selection to this genus."
HELP_TARGET_SPECIES = "Restrict the selection to this species."
HELP_OUTGROUP_ACCESSION = "Accession to fetch and set aside as the outgroup."
HELP_NO_OUTGROUP = "Do not root with an outgroup."
HELP_ROOT = (
    "How phylo roots tree/tree.nwk: outgroup (on the branch to the outgroup), midpoint, "
    "or none (as the tree builder wrote it). Default: outgroup when an outgroup is "
    "staged and --no-outgroup is not given, else none."
)
HELP_KEEP_ROOT = (
    "Keep the tree's own root (a tree from phylo --root midpoint) instead of rerooting "
    "on the outgroup; the outgroup is still read for --remove-outgroup."
)
HELP_KEEP_FILES = "Keep download and scratch intermediates."
HELP_SKETCH = (
    "Write a sourmash sketch of each genome to sketches/ (k=21,31,51, scaled=1000). "
    "Default: when sourmash can run; --sketch requires it, --no-sketch skips it."
)
HELP_READS_SKETCH = (
    "Sketch each run's reads with sourmash (k=21,31,51, scaled=1000, with abundances) "
    "into assemblies/<run>/reads.sig.zip while the run is assembled; never in sketches/. "
    "Default: when sourmash can run; --reads-sketch requires it, --no-reads-sketch skips it."
)
HELP_SCREEN_READS = (
    "Screen each run's reads sketch before assembling it: a run is excused when its GTDB "
    "genus differs from the submitted genus (taxon_mismatch), when the top species holds "
    "less than --screen-min-fraction of the reads (low_match_fraction), or when an accepted "
    "run of the same biosample contains it (duplicate_isolate). Needs sourmash, "
    "--gtdb-sketch and --gtdb-lineages."
)
HELP_SCREEN_MIN_FRACTION = (
    "With --screen-reads: minimum abundance-weighted fraction of the reads in the top GTDB species."
)
HELP_WORKDIR = "Working directory."
HELP_WORKDIR_CREATED = "Working directory (created)."
HELP_GTDB_RELEASE = "GTDB release (tsv source)."
HELP_GTDB_VERSION = "GTDB table: bac120 or ar53 (tsv source)."
HELP_METADATA_PATH = (
    "Use this GTDB metadata table instead of downloading. With --source tsv, "
    "-r/--release and --gtdb-version are still required."
)
HELP_NODOWNLOAD = "Reuse a GTDB table already present in the workdir."
HELP_LIMIT = "Keep at most N GTDB genomes, round-robin over species by CheckM quality."
HELP_DATASET = "GTDB dataset: all or rep."
HELP_LEVEL = "family, genus or species."
HELP_ALL_GENOMES = "Use all genomes, not only the representatives."
HELP_BOOTSTRAP = (
    "Bootstrap replicates. 0 turns bootstrapping off; IQ-TREE needs at least 1000 when it is on."
)
HELP_REFERENCE = (
    "Reference genome: a filename, 'medoid' (the default: the genome nearest to all others by "
    "sketch ANI, when the genomes are sketched and sourmash is at hand) or 'first'."
)
HELP_CORE_FRACTION = (
    "Keep only the variable columns where at least this fraction of the genomes have a called "
    "base (a soft core, 0 to 1). Unset keeps every variable column and the log reports what "
    "0.95 would keep."
)
HELP_MSA_SOURCE = (
    "aligner, snptype or coregene (a concatenated core-gene alignment from the annotations "
    "the annotate stage wrote)."
)
HELP_COREGENE = (
    "Core-gene alignment builder for --msa-source coregene: panaroo, ppanggolin, scarap "
    "or orthofinder. "
    "Tuning goes through --aligner-arg. panaroo: threshold= is the identity at which it "
    "groups genes (0.98; lower it, 0.95 then 0.90, when the core comes out small on a "
    "divergent set), and also filtered=false, merge_paralogs, remove_invalid_genes, "
    "family_threshold, len_dif_percent. ppanggolin: identity= and coverage= of its "
    "protein clustering (0.8 each), partition=softcore|core|persistent (softcore, at "
    "--core-threshold), single_copy=false, dup_margin, mode, nb_of_partitions, "
    "translation_table. scarap: method (FH), seeds (100), max_core_genes (0, no cap), "
    "core_prefilter (the smaller of 0.90 and --core-threshold). orthofinder: "
    "search=diamond|diamond_ultra_sens|blastp|mmseqs (diamond), inflation (MCL, 1.2)."
)
HELP_CORE_THRESHOLD = (
    "With --msa-source coregene: a gene is core when present in at least this fraction "
    "of the genomes (0.95)."
)
HELP_CLEAN_MODE = (
    "With --msa-source coregene: Panaroo's cleaning of annotation errors, strict (default), "
    "moderate or sensitive."
)
HELP_ALIGNER_ARG = (
    "Aligner, SNP typer, masker or tree builder tuning as key=value (repeatable), "
    "e.g. kmer=15 (sibeliaz), seed_weight=11 (progressivemauve) or model=GTR+G "
    "(iqtree, raxmlng)."
)
HELP_DEREP_TOOL_ARG = "Tool tuning as key=value (repeatable), e.g. mode=greedy."
HELP_KEEPER = (
    "Representative choice per cluster: quality (manifest CheckM values and genome "
    "N50), gtdb (a GTDB species representative first, then quality) or tool "
    "(adapter's own)."
)
HELP_PROCESS_SIZE = "Chunk size; when set and exceeded, two-stage chunking runs for any tool."
HELP_NUM_PROCESSES = (
    "Parallel stage-1 chunk workers (threads split across them). "
    "0 = auto (~threads/4, capped by cores)."
)
HELP_ANNOTATE_PROCESSES = (
    "Genomes annotated at once (threads split across them). 0 = auto (~threads/4, capped by cores)."
)
HELP_RUN_NUM_PROCESSES = (
    "Parallel workers, threads split across them: dereplicate's stage-1 chunk "
    "workers and, with --msa-source coregene, the genomes annotated at once. "
    "0 = auto (~threads/4, capped by cores)."
)
HELP_REFINE = (
    "Seed of a lineage to keep at finer resolution (filename, stem or accession; "
    "repeatable): every cluster holding a seed is dereplicated again at --refine-ani."
)
HELP_REFINE_ANI = "Secondary ANI for the clusters of the --refine seeds; above --secondary-ani."
HELP_PRE_PRIMARY_ANI = "Stage-1 (intra-chunk) primary ANI; defaults to --primary-ani."
HELP_PRE_SECONDARY_ANI = "Stage-1 (intra-chunk) secondary ANI; defaults to --secondary-ani."
HELP_REDUCE = (
    "Taxonomy-aware reduction after ANI: none, species, or genus (one representative per taxon)."
)
HELP_VIRUS = "Pass virus-tuned parameters to dRep (--tool drep); the other tools do not read it."
HELP_TARGET_REPS = (
    "Target representative count: search --secondary-ani to land near it "
    "(0 = off; re-runs dereplication per search step)."
)
HELP_ALLOW_INCOMPLETE = "Proceed with a warning when the input genome set is incomplete."
HELP_NAMING = (
    "How internal nodes are named without --node-basename: 'leafset' (a hash of the whole "
    "leaf set) or 'lca' (a hash of the smallest leaf label under each child, so most names "
    "survive genomes added with later-sorting accessions)."
)
HELP_NODE_BASENAME = (
    "Name internal nodes <basename><n>. Without it, internal nodes receive names "
    "derived from a hash of their descendant leaves."
)
HELP_ROOT_NAME = "Label of the top node. FlexTaxD expects 'root'."
HELP_REMOVE_OUTGROUP = "Leave the outgroup out of the taxonomy after rooting."
HELP_INCLUDE_DEREPLICATED = "List redundant genomes under their representative in the taxonomy."
HELP_COLLAPSE_SUPPORT = "Merge nodes whose support is below this fraction into their parent."
HELP_ANNOTATOR = "Annotator: bakta (default; needs its database, --bakta-db or BAKTA_DB) or prokka."
HELP_BAKTA_DB = (
    "bakta database directory (bakta_db download --output DIR --type light|full); "
    "BAKTA_DB when unset."
)
HELP_COLLAPSE_LENGTH = "Merge nodes whose branch is shorter than this length into their parent."
HELP_COLLAPSE_DIAMETER = (
    "Merge every clade whose leaf-to-leaf distance is within this value into one node, "
    "keeping the fewest such clades (in the tree's own units)."
)
HELP_RANKS = (
    "Add a rank column to tree2tax.tsv: 'strain' on every genome, 'no rank' on internal "
    "nodes unless --rank-red names them."
)
HELP_RANK_RED = (
    "With --ranks: NAME=RED bands, e.g. species=0.95,genus=0.8; an internal node takes the "
    "rank with the largest threshold not above its relative evolutionary divergence, once "
    "per lineage. Needs branch lengths and a rooted tree."
)
HELP_VERSIONS_OUT = "Write resolved tool versions (YAML fragment) here."
HELP_SKETCHES_DIR = (
    "A sketches/ directory of <record name>.sig.zip files (k=21, 31, 51; scaled=1000) "
    "made from these genomes. A sourmash tool reads them instead of sketching; "
    "other tools ignore it."
)

# Canonical stage order per lineage. Used to show progress (`status`) and by
# `run --dry-run` to print the chain.
PIPELINE_BACTERIAL = ("metadata", "genome", "dereplicate", "phylo", "tree2tax")
PIPELINE_VIRAL = ("vmetadata", "vgenome", "dereplicate", "phylo", "tree2tax")
# Offline chain: local genomes staged by `ingest` instead of downloaded.
PIPELINE_LOCAL = ("ingest", "dereplicate", "phylo", "tree2tax")
# Reads chain: sequencing runs selected from ENA/SRA and assembled.
PIPELINE_READS = ("reads", "assemble", "dereplicate", "phylo", "tree2tax")


def _chain(label: str, stages: tuple[str, ...]) -> str:
    return f"{label}: " + " -> ".join(stages)


# Paragraphs are separated by blank lines because Rich joins single newlines.
APP_EPILOG = "\n\n".join(
    [
        "Typical order of commands (or let 'run' chain them; "
        "'status -wd WD' says what comes next):",
        _chain("bacterial", PIPELINE_BACTERIAL),
        _chain("viral", PIPELINE_VIRAL),
        _chain("local genomes", PIPELINE_LOCAL),
        _chain("sequencing reads", PIPELINE_READS),
        "Add 'snptype' between dereplicate and phylo when the SNP tables are a deliverable. "
        "Global options go before the command name: repgenr --container docker dereplicate ...",
    ]
)

PANEL_PIPELINE = "Pipeline"
PANEL_ENTRY = "Entry points: select and fetch genomes"
PANEL_CORE = "Core stages"
PANEL_INSPECT = "Inspect a dereplication"
PANEL_ENV = "Environment and diagnostics"
PANEL_STEPS = "Nextflow data-channel steps"

# Panel -> commands in display order. The single source of truth for the
# grouped --help, the rendered command reference, and the panel test.
COMMAND_PANELS: dict[str, tuple[str, ...]] = {
    PANEL_PIPELINE: ("run", "status"),
    PANEL_ENTRY: ("metadata", "genome", "vmetadata", "vgenome", "ingest", "reads", "assemble"),
    PANEL_CORE: ("dereplicate", "derep-refine", "annotate", "snptype", "phylo", "tree2tax"),
    PANEL_INSPECT: ("census", "glance", "cluster-summary", "derep-unpack", "derep-stock", "sketch"),
    PANEL_ENV: ("list-tools", "doctor", "versions"),
    PANEL_STEPS: (
        "genome-fetch",
        "dereplicate-chunk",
        "dereplicate-merge",
        "annotate-build",
        "annotate-index",
        "phylo-build",
        "tree2tax-relations",
        "assemble-run",
        "genome-qc",
        "reads-gather",
    ),
}

COMMAND_ORDER: tuple[str, ...] = tuple(c for cmds in COMMAND_PANELS.values() for c in cmds)


class _PanelOrderedGroup(TyperGroup):
    """--help lists commands in COMMAND_ORDER (pipeline order), not import order."""

    def list_commands(self, ctx):  # type: ignore[no-untyped-def]
        rank = {n: i for i, n in enumerate(COMMAND_ORDER)}
        return sorted(self.commands, key=lambda n: (rank.get(n, len(rank)), n))


app = typer.Typer(
    cls=_PanelOrderedGroup,
    add_completion=False,
    no_args_is_help=True,
    help="RepGenR2: modular genome dereplication, alignment, SNP typing and phylogenetics.",
    epilog=APP_EPILOG,
)


def _phylo_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    # Neither snp/ (the snptype stage's tables) nor tree/msa/ (phylo's own
    # typing pass, --msa-source snptype) is an input: phylo types the genome
    # set itself, and the MSA stamp under tree/msa/ decides whether the typed
    # alignment is reused. Declaring its own output would force a rerun on
    # every second invocation.
    inputs = [
        ctx.genomes_dir if getattr(params, "all_genomes", False) else ctx.representatives_dir,
        ctx.outgroup_dir,
    ]
    if getattr(params, "msa_source", None) == "coregene":
        # The index carries the GFF digests, so it stands for the whole cache.
        inputs.append(ctx.annotations_dir / ANNOTATIONS_TSV)
    return inputs


def _annotate_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    # The genome set and the outgroup; the per-record cache inside the stage
    # keeps unchanged annotations when the fingerprint changes. The database
    # is an input too: bakta's version.json identifies it, so a database
    # updated in place, or a changed BAKTA_DB, reruns the stage instead of
    # skipping with annotations made against the previous one.
    inputs = [
        ctx.genomes_dir if getattr(params, "all_genomes", False) else ctx.representatives_dir,
        ctx.outgroup_dir,
    ]
    db = _annotator_database(params)
    if db is not None:
        inputs.append(db / "version.json")
    return inputs


def _annotator_database(params: Any) -> Path | None:
    """The annotator's database directory, from the flag or its environment
    variable; None for a tool that needs none or a path that does not resolve."""
    from ..annotators.base import registry as annotator_registry

    try:
        annotator = annotator_registry.create(getattr(params, "tool", ""))
    except Exception:  # noqa: BLE001 - an unknown tool is refused by the stage
        return None
    if not annotator.requires_database:
        return None
    from ..annotators.bakta import bakta_db_from_env

    raw = getattr(params, "bakta_db", None) or bakta_db_from_env()
    return Path(raw).expanduser() if raw else None


def _metadata_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    if getattr(params, "metadata_path", None):
        return [Path(params.metadata_path)]
    release = getattr(params, "release", None)
    version = getattr(params, "version", None)
    reuses = getattr(params, "nodownload", False) and getattr(params, "source", "tsv") == "tsv"
    if reuses and release and version:
        from ..stages.metadata import workdir_tables

        try:
            # Both naming schemes; the absent one digests to a stable sentinel.
            return workdir_tables(ctx.workdir, release, version)
        except ValueError:
            return []  # a malformed --release; the stage reports it
    return []


def _ingest_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    paths: list[Path] = []
    if getattr(params, "genomes_dir", None):
        paths.append(Path(params.genomes_dir).expanduser())
    if getattr(params, "selection", None):
        paths.append(Path(params.selection).expanduser())
    outgroup: str = getattr(params, "outgroup", None) or ""
    outgroup_file = bool(outgroup) and Path(outgroup).expanduser().is_file()
    # Each source workdir's selection (by content) and genomes/ (FASTA files
    # only); its outgroup/ only when --outgroup may name a genome there. The
    # source manifest, which gives only the source label, is not digested.
    for wd in getattr(params, "from_workdirs", None) or []:
        root = Path(wd).expanduser()
        paths += [root / SELECTION_TSV, root / "genomes"]
        if outgroup and not outgroup_file:
            paths.append(root / "outgroup")
    if outgroup_file:
        paths.append(Path(outgroup).expanduser())
    return paths


def _derep_stock_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    action = getattr(params, "action", "")
    name = getattr(params, "name", None) or ""
    if action == "pack":
        return [ctx.derep_dir / CLUSTERS_TSV, ctx.representatives_dir]
    if action == "unpack":
        # The stored run's files are listed one by one: directories are
        # digested from their FASTA files only, and the stored run holds
        # tables. The live derep outputs are inputs too: unpack replaces
        # them, so a repeat unpack after a new dereplicate must restore the
        # run again.
        stored = ctx.derep_dir / "stock" / name
        return [
            stored / CLUSTERS_TSV,
            stored / GENOME_STATUS_TSV,
            stored / "record.json",
            stored / "representatives",
            ctx.derep_dir / CLUSTERS_TSV,
            ctx.representatives_dir,
        ]
    return []


def _tree2tax_inputs(ctx: WorkdirContext, params: Any) -> list[Path]:
    paths = [ctx.tree_dir / TREE_NWK]
    if getattr(params, "include_dereplicated", False):
        paths.append(ctx.derep_dir / CLUSTERS_TSV)
    # Read when present: the outgroup files decide the rooting and
    # --remove-outgroup, segments.tsv adds rows to genomes_map.tsv. A missing
    # one digests to the ABSENT sentinel, so appearing later reruns the stage.
    paths += [
        ctx.workdir / SEGMENTS_TSV,
        ctx.workdir / OUTGROUP_ACCESSION_TXT,
        ctx.outgroup_dir,
    ]
    return paths


# What each stage reads, for the input digests in the resume fingerprint:
# stage -> callable(ctx, params) -> paths (directories are digested from file
# metadata, files by content; see core.inputs). Conditional edges (--all-genomes,
# --msa-source snptype, --include-dereplicated) live in the helpers above.
# Stages not listed digest no inputs and fingerprint on params alone.
STAGE_INPUTS: dict[str, Any] = {
    # metadata downloads its table unless --metadata-path names a local one,
    # or --nodownload reuses the one in the workdir; that table is its input.
    "metadata": _metadata_inputs,
    # ingest reads paths outside the workdir; they are keyed absolute.
    "ingest": _ingest_inputs,
    # reads is network-only; an accession list is its one file input.
    "reads": lambda ctx, p: [Path(p.accession_file)] if p.accession_file else [],
    "assemble": lambda ctx, p: [
        ctx.workdir / READS_TSV,
        *([Path(p.outgroup)] if p.outgroup else []),
    ],
    "vmetadata": lambda ctx, p: [],
    "genome": lambda ctx, p: [ctx.workdir / SELECTION_TSV],
    # vgenome WRITES selection.tsv, so its inputs are the vmetadata download
    # artifacts (records path and legacy BV-BRC tables; absent ones digest to
    # the stable sentinel).
    "vgenome": lambda ctx, p: [
        ctx.workdir / "virus_download_wd" / "download.fa",
        ctx.workdir / "virus_download_wd" / "virus_records.json",
        ctx.workdir / "virus_download_wd" / "metadata_base.tsv",
        ctx.workdir / "virus_download_wd" / "metadata_ncbi.tsv",
    ],
    "dereplicate": lambda ctx, p: [ctx.genomes_dir],
    # Its own output is an input, digested again after the run (below), so a
    # repeat is skipped and any later dereplicate or unpack reruns it. The
    # manifest is left out: the stage rewrites it too.
    "derep_refine": lambda ctx, p: [ctx.derep_dir / CLUSTERS_TSV, ctx.genomes_dir],
    "snptype": lambda ctx, p: [
        ctx.genomes_dir if getattr(p, "all_genomes", False) else ctx.representatives_dir
    ],
    "annotate": _annotate_inputs,
    "phylo": _phylo_inputs,
    "tree2tax": _tree2tax_inputs,
    # Auxiliary stages: recorded so status shows them and a repeat skips.
    "glance": lambda ctx, p: [ctx.genomes_dir],
    "derep_unpack": lambda ctx, p: [ctx.derep_dir / CLUSTERS_TSV, ctx.genomes_dir],
    "cluster_summary": lambda ctx, p: [ctx.derep_dir / CLUSTERS_TSV],
    "derep_stock": _derep_stock_inputs,
    # The genome set: the files and the selection that names them.
    "sketch": lambda ctx, p: [ctx.genomes_dir, ctx.outgroup_dir, ctx.workdir / SELECTION_TSV],
}


def _vmetadata_deliverables(ctx: WorkdirContext, params: Any) -> list[Path]:
    download_wd = ctx.workdir / "virus_download_wd"
    paths = [download_wd / "metadata_base.tsv"]
    source = getattr(params, "source", None)
    if source == "ncbi_virus":
        paths.append(download_wd / "virus_records.json")
    elif source == "bvbrc":
        paths += [download_wd / "download.fa", download_wd / "metadata_ncbi.tsv"]
    return paths


def _derep_stock_deliverables(ctx: WorkdirContext, params: Any) -> list[Path]:
    action = getattr(params, "action", "")
    if action == "pack":
        return [ctx.derep_dir / "stock" / (getattr(params, "name", None) or "") / CLUSTERS_TSV]
    if action == "unpack":
        return [ctx.derep_dir / CLUSTERS_TSV, ctx.representatives_dir]
    # delete runs as a query and is never fingerprinted, so the resume check
    # never reaches it; a delete record from an older version (read by
    # doctor) has nothing to check.
    return []


def _selected_genome_files(ctx: WorkdirContext) -> list[Path]:
    """Every genome file ``selection.tsv`` promises, outgroup included.

    Listed as deliverables of the stage that wrote the genome set, so one
    genome or the outgroup deleted by hand reruns that stage (which restores
    only what is absent) instead of being skipped because genomes/ is not
    empty. Accessions recorded as unavailable are excused.
    """
    from ..core.contracts import read_selection
    from ..core.integrity import excused_accessions

    selection = ctx.workdir / SELECTION_TSV
    if not selection.is_file():
        return []
    try:
        rows = read_selection(selection)
    except (OSError, ValueError, RepGenRError):
        return []  # the stage itself reports an unreadable selection
    excused = excused_accessions(ctx.workdir)
    return [
        (ctx.outgroup_dir if row.is_outgroup else ctx.genomes_dir) / row.filename
        for row in rows
        if row.accession not in excused
    ]


def _genome_deliverables(ctx: WorkdirContext, params: Any) -> list[Path]:
    """The genome directory and manifest, plus every file the selection promises."""
    return [ctx.genomes_dir, ctx.workdir / MANIFEST_FILENAME, *_selected_genome_files(ctx)]


def _genome_set_deliverables(ctx: WorkdirContext) -> list[Path]:
    """What every entry path that writes a genome set leaves in the workdir."""
    return [
        ctx.genomes_dir,
        ctx.workdir / SELECTION_TSV,
        ctx.workdir / MANIFEST_FILENAME,
        *_selected_genome_files(ctx),
    ]


def _sketch_deliverables(ctx: WorkdirContext) -> list[Path]:
    from ..core.sketches import expected_sketch_files

    return expected_sketch_files(ctx.workdir)


def _dereplicate_deliverables(ctx: WorkdirContext, params: Any) -> list[Path]:
    """The four derep tables and directories, plus each listed representative.

    A representative deleted by hand reruns the stage instead of leaving
    phylo to refuse a representatives/ directory out of step with
    clusters.tsv. All four outputs are listed: doctor fails on a missing
    genome_status.tsv and asks for a rerun, which must then not be skipped.
    """
    paths = [
        ctx.derep_dir / CLUSTERS_TSV,
        ctx.derep_dir / GENOME_STATUS_TSV,
        ctx.derep_dir / CLUSTER_SUMMARY_TSV,
        ctx.representatives_dir,
    ]
    clusters = ctx.derep_dir / CLUSTERS_TSV
    if clusters.is_file():
        try:
            paths += [ctx.representatives_dir / name for name in read_clusters(clusters)]
        except (OSError, ValueError, RepGenRError):
            pass  # an unreadable table: the stage reruns on its own terms
    return paths


# What each stage writes that downstream stages or the user rely on: stage ->
# callable(ctx, params) -> paths. A completed stage whose fingerprint matches
# is skipped only when all of these exist (a directory must also be
# non-empty); otherwise it reruns, so outputs deleted by hand are rebuilt
# without --force. The callables read only the ctx attributes that
# core.doctor's stand-in context provides, and params may be the recorded
# params dict as a namespace (doctor), so conditional entries use getattr.
STAGE_DELIVERABLES: dict[str, Any] = {
    "metadata": lambda ctx, p: [ctx.workdir / SELECTION_TSV, ctx.workdir / MANIFEST_FILENAME],
    "ingest": lambda ctx, p: _genome_set_deliverables(ctx),
    "reads": lambda ctx, p: [ctx.workdir / READS_TSV],
    "assemble": lambda ctx, p: _genome_set_deliverables(ctx),
    "vmetadata": _vmetadata_deliverables,
    # genome reads selection.tsv and the manifest; it writes the genome files.
    "genome": _genome_deliverables,
    "vgenome": lambda ctx, p: _genome_set_deliverables(ctx),
    "dereplicate": _dereplicate_deliverables,
    "derep_refine": lambda ctx, p: [ctx.derep_dir / REFINED_TSV, ctx.derep_dir / CLUSTERS_TSV],
    "snptype": lambda ctx, p: [ctx.snp_dir / CORE_SNP_FASTA],
    "annotate": lambda ctx, p: [ctx.annotations_dir / ANNOTATIONS_TSV],
    "phylo": lambda ctx, p: [ctx.tree_dir / TREE_NWK],
    "tree2tax": lambda ctx, p: [
        ctx.workdir / TREE2TAX_TSV,
        ctx.workdir / GENOMES_MAP_TSV,
        ctx.workdir / TREE2TAX_NODES_TSV,
    ],
    # glance's plots are absent when no similarity falls within the plot
    # bounds, so only the dendrogram is checked (a tool that returns none is
    # warned about and reruns each time).
    "glance": lambda ctx, p: [ctx.workdir / "glance_clustering_dendrogram.pdf"],
    # With --no-representant and only singleton clusters, unpacked/ is
    # legitimately left empty, so it is checked only when representatives
    # are unpacked too (every cluster then yields a subdirectory).
    "derep_unpack": lambda ctx, p: (
        [] if getattr(p, "no_representant", False) else [ctx.derep_dir / "unpacked"]
    ),
    "cluster_summary": lambda ctx, p: [ctx.derep_dir / CLUSTER_SUMMARY_TSV],
    "derep_stock": _derep_stock_deliverables,
    # One sketch per genome file selection.tsv names; a sketch deleted by hand
    # reruns the command. The genome-writing stages do not list sketches:
    # a genome set without them is complete.
    "sketch": lambda ctx, p: _sketch_deliverables(ctx),
}


_LOGGED_MISSING = 5  # missing deliverables named one per line before a count


def _deliverable_present(path: Path, *, fasta_dir: bool = False) -> bool:
    if path.is_dir():
        # A genome directory counts only with a genome in it: one holding a
        # leftover x.fasta.tmp is as empty as one holding nothing.
        if fasta_dir:
            return bool(list_fasta(path))
        # Dotfiles do not count: Finder leaves .DS_Store in a directory it
        # showed, and exFAT keeps ._ AppleDouble companions.
        return any(not entry.name.startswith(".") for entry in path.iterdir())
    return path.exists()


def missing_deliverables(ctx: Any, stage_name: str, params: Any) -> list[Path]:
    """Declared deliverables of ``stage_name`` that are absent (or empty dirs)."""
    spec = STAGE_DELIVERABLES.get(stage_name)
    if spec is None:
        return []
    fasta_dirs = {getattr(ctx, "genomes_dir", None), getattr(ctx, "representatives_dir", None)}
    return [
        path
        for path in spec(ctx, params)
        if not _deliverable_present(path, fasta_dir=path in fasta_dirs)
    ]


def _deliverables_state(ctx: Any, stage_name: str, params: Any) -> list[tuple[str, object]]:
    """A cheap snapshot of a stage's deliverables: one stat per file.

    Files by (size, mtime_ns), directories by their flat stat digest, absent
    paths as None; compared before and after a refused run to tell a refusal
    that touched nothing from one that left partial outputs.

    Limitation: the manifest is WAL-mode, and its writes stay in the -wal
    file until the connection closes, so ``manifest.sqlite`` shows no change
    in size or mtime here. A stage that wrote only the manifest and then
    refused would therefore have its record restored. No current stage
    does: each refuses before it writes, or also writes a file deliverable.
    """
    spec = STAGE_DELIVERABLES.get(stage_name)
    state: list[tuple[str, object]] = []
    for path in spec(ctx, params) if spec is not None else []:
        if path.is_dir():
            # Every file, not only FASTA: a refusal that left a partial
            # x.fasta.tmp behind has changed the directory.
            state.append((str(path), dir_stat_digest(path, fasta_only=False)))
        elif path.exists():
            st = path.stat()
            state.append((str(path), (st.st_size, st.st_mtime_ns)))
        else:
            state.append((str(path), None))
    return state


def deliverable_label(workdir: Path, path: Path) -> str:
    """A deliverable path as shown to the user: relative to the workdir."""
    try:
        return str(path.relative_to(workdir))
    except ValueError:
        return str(path)


# Stages whose result also depends on the manifest's genome rows (taxonomy,
# derep status, CheckM quality), digested from ordered query results.
# "dereplicate" reads the manifest for the quality-aware keeper and --reduce
# taxonomy grouping, so a manifest-only edit (no genome file touched) must
# still invalidate a prior resume.
_MANIFEST_INPUT_STAGES = frozenset({"tree2tax", "dereplicate", "cluster_summary"})

# Param flags that turn a stage invocation into a pure query (list/preview
# modes that write no pipeline outputs). Such invocations bypass the resume
# machinery entirely: no skip check, no dirty marker, no fingerprint stamp --
# otherwise a query would either be wrongly skipped or would restamp/dirty the
# record of the last real run.
QUERY_ONLY_FLAGS: dict[str, tuple[str, ...]] = {
    "genome": ("accession_list_only",),
    "vmetadata": ("list_targets",),
    "vgenome": ("glance",),
}


# Invocations that rewrite their own declared inputs (derep-stock unpack
# restores derep/; metadata --nodownload downloads the table it then reuses
# when none is present yet): the record is stamped with digests taken after
# the run, so an identical repeat matches the new state and skips.
_REDIGEST_AFTER_RUN: dict[str, Any] = {
    "derep_stock": lambda p: getattr(p, "action", None) == "unpack",
    "derep_refine": lambda p: True,
    "metadata": lambda p: getattr(p, "nodownload", False),
}


# Refusals a stage can check before the harness marks its record incomplete:
# stage -> callable(ctx, params) raising UserInputError/WorkdirError. Called
# only when the stage will run (after the resume skip check), never on a skip.
# Used where one record serves several invocations (derep-stock's named runs) or
# where a refused re-run would otherwise dirty a finished record (ingest, and
# dereplicate with a selected genome missing), so the record of the last
# finished run stays clean.
def _derep_stock_precheck(ctx: WorkdirContext, params: Any) -> None:
    from ..stages.derep_stock import precheck

    precheck(ctx, params)


def _ingest_precheck(ctx: WorkdirContext, params: Any) -> None:
    from ..stages.ingest import precheck

    precheck(ctx, params)


def _dereplicate_precheck(ctx: WorkdirContext, params: Any) -> None:
    from ..stages.dereplicate import precheck

    precheck(ctx, params)


def _assemble_precheck(ctx: WorkdirContext, params: Any) -> None:
    from ..stages.assemble import precheck

    precheck(ctx, params)


def _phylo_precheck(ctx: WorkdirContext, params: Any) -> None:
    from ..stages.phylo import precheck

    precheck(ctx, params)


_STAGE_PRECHECKS: dict[str, Any] = {
    "derep_stock": _derep_stock_precheck,
    "dereplicate": _dereplicate_precheck,
    "ingest": _ingest_precheck,
    # A wrong database path on a rerun must not leave the finished record dirty.
    "assemble": _assemble_precheck,
    # A core-gene tree over unannotated genomes is refused before any tool runs.
    "phylo": _phylo_precheck,
}

# Query modes keyed on a value rather than a flag.
QUERY_ONLY_PREDICATES: dict[str, Any] = {
    # delete is never skipped: a repeat delete must report the unknown run.
    "derep_stock": lambda p: getattr(p, "action", None) in ("list", "delete"),
}


def _is_query_only(stage_name: str, params: Any) -> bool:
    if any(getattr(params, flag, False) for flag in QUERY_ONLY_FLAGS.get(stage_name, ())):
        return True
    predicate = QUERY_ONLY_PREDICATES.get(stage_name)
    return bool(predicate and predicate(params))


def _stage_input_digests(ctx: WorkdirContext, stage_name: str, params: Any) -> dict[str, str]:
    """Digest a stage's declared inputs; empty for stages with no declaration."""
    spec = STAGE_INPUTS.get(stage_name)
    if spec is None:
        return {}
    digests = inputs_digest(ctx.workdir, spec(ctx, params))
    # Opening the manifest creates it; a missing workdir has none to digest,
    # and the stage itself then reports the missing input.
    if stage_name in _MANIFEST_INPUT_STAGES and ctx.workdir.is_dir():
        digests["manifest"] = manifest_digest_for_stage(stage_name, ctx.manifest)
    return digests


def _env_fragment() -> dict[str, Any]:
    """Result-affecting execution environment for the fingerprint.

    The container backend/platform/wave selection changes which tool builds run
    (and so can change results); the engine binary, cache directory, and extra
    mounts are plumbing and deliberately excluded, like _NON_RESULT_PARAMS.
    """
    from ..core.containers import result_env_fragment

    return result_env_fragment()


def _derep_help(*, auto: bool = True) -> str:
    from ..core.plugins import tool_choices_help
    from ..dereplicators.base import registry

    return tool_choices_help(registry, auto=auto)


def _tree_help(*, auto: bool = True) -> str:
    from ..core.plugins import tool_choices_help
    from ..treebuilders.base import registry

    return tool_choices_help(registry, auto=auto)


def _aligner_help() -> str:
    from ..aligners.base import registry
    from ..core.plugins import tool_choices_help

    return tool_choices_help(registry, auto=False)


def _snp_help() -> str:
    from ..core.plugins import tool_choices_help
    from ..snptypers.base import registry

    return tool_choices_help(registry, auto=False, prefix="SNP typer: ")


def _assembler_help() -> str:
    from ..assemblers.base import registry
    from ..core.plugins import tool_choices_help

    return tool_choices_help(registry, auto=True, prefix="Assembler: ")


def _polisher_help() -> str:
    from ..core.plugins import tool_choices_help
    from ..polishers.base import registry

    return tool_choices_help(
        registry, auto=True, prefix="Polisher for long-read assemblies: none, "
    )


def _classifier_help() -> str:
    from ..classifiers.base import registry
    from ..core.plugins import tool_choices_help

    return tool_choices_help(registry, auto=True, prefix="Classifier: none, ")


def _outgroup_builder_help() -> str:
    from ..viral._outgroup import distance_matrix_builders

    names = ", ".join(distance_matrix_builders()) or "(none registered)"
    return f"Tree builder used for the outgroup distance matrix. Accepted: {names}."


def _mask_help() -> str:
    from ..core.plugins import tool_choices_help
    from ..maskers.base import registry

    return tool_choices_help(
        registry, auto=False, prefix="Recombination masking of the SNP alignment: none, "
    )


def _mask_help_msa() -> str:
    """--mask help on the commands where it applies only with --msa-source snptype."""
    return _mask_help() + " Needs --msa-source snptype."


def _require_choice(value: str, choices: AbstractSet[str], label: str) -> None:
    if value not in choices:
        raise UserInputError(
            f"Invalid {label} {value!r}. Choose from: {', '.join(sorted(choices))}."
        )


def _require_unit_interval(value: float | None, label: str) -> None:
    if value is not None and not (0.0 < value <= 1.0):
        raise UserInputError(f"{label} must be in (0, 1], got {value}.")


# Parameters that change how work is scheduled but not the result, so they are
# excluded from the resume fingerprint -- changing --threads / --num-processes
# must not force an otherwise-identical stage to recompute from scratch.
# allow_incomplete only gates the input-completeness refusal; on complete
# inputs it changes nothing, so it must not invalidate the resume cache.
# sketch (--sketch/--no-sketch of the genome-writing stages) adds sketches/,
# which is not a deliverable, so it changes no result either; a workdir
# finished without sketches gets them from `repgenr sketch`.
# strict_sources (ingest --strict-sources) only refuses a run; the genomes
# staged are the same, so a repeat with it skips like one without.
_NON_RESULT_PARAMS = frozenset(
    {"threads", "num_processes", "allow_incomplete", "sketch", "strict_sources"}
)


def _phylo_fingerprint_params(payload: dict) -> dict:
    from ..stages.phylo import fingerprint_params

    return fingerprint_params(payload)


# Stage-specific views of the parameters before they are hashed, for settings
# that some tools of a family ignore: stage -> callable(dict) -> dict.
_FINGERPRINT_NORMALISERS: dict[str, Any] = {
    "phylo": _phylo_fingerprint_params,
}


# Fingerprint format version. Bumping it guarantees fingerprints from older
# releases never false-match, so old workdirs rerun once instead of skipping
# against semantics they were not computed under.
_FINGERPRINT_VERSION = 2


def _stage_fingerprint(
    stage_name: str, params: object, inputs: dict[str, str], env: dict[str, Any]
) -> str:
    """Stable hash of a stage invocation, used to skip already-completed work.

    Built from the stage name, the parameter object (a dataclass), the digests
    of the stage's declared inputs, and the result-affecting environment
    (container identity), so a skip means "same request, same inputs, same
    execution environment". Paths and other non-JSON values are stringified.
    Non-result params (thread/worker counts) are excluded so they do not
    invalidate the resume cache.
    """
    if dataclasses.is_dataclass(params) and not isinstance(params, type):
        payload: dict = dataclasses.asdict(params)
    else:
        payload = dict(vars(params))
    payload = {k: v for k, v in payload.items() if k not in _NON_RESULT_PARAMS}
    normalise = _FINGERPRINT_NORMALISERS.get(stage_name)
    if normalise is not None:
        payload = normalise(payload)
    blob = json.dumps(
        {
            "fpv": _FINGERPRINT_VERSION,
            "stage": stage_name,
            "params": payload,
            "inputs": inputs,
            "env": env,
        },
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"repgenr {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
    container: str = typer.Option(
        "none",
        "--container",
        envvar="REPGENR_CONTAINER",
        help="Run external tools in containers: none, docker, or singularity.",
    ),
    container_engine: str | None = typer.Option(
        None,
        "--container-engine",
        envvar="REPGENR_CONTAINER_ENGINE",
        help="Engine binary override (e.g. apptainer, podman).",
    ),
    container_cache: str | None = typer.Option(
        None,
        "--container-cache",
        envvar="REPGENR_CONTAINER_CACHE",
        help="Directory for Singularity .sif images and their cache (large; can be external).",
    ),
    platform: str | None = typer.Option(
        None,
        "--platform",
        envvar="REPGENR_CONTAINER_PLATFORM",
        help="Container platform, e.g. linux/amd64 for emulated BioContainers on arm64.",
    ),
    wave: bool = typer.Option(
        False,
        "--wave/--no-wave",
        envvar="REPGENR_WAVE",
        help="Resolve images for multi-tool adapters via the Seqera Wave CLI.",
    ),
    bin_dir: list[str] | None = typer.Option(
        None,
        "--bin-dir",
        metavar="TOOL=DIR",
        help=(
            "Put DIR first on PATH for one tool only, e.g. gubbins=ENV/bin for a satellite "
            "conda environment. Repeatable; also REPGENR_BIN_DIRS='tool=dir,tool=dir'."
        ),
    ),
    force: bool = typer.Option(
        False,
        "--force/--no-force",
        "-f",
        envvar="REPGENR_FORCE",
        help="Re-run a stage even if it already completed with the same parameters.",
    ),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose (DEBUG) logging."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Only warnings and errors."),
) -> None:
    """RepGenR2 top-level entry point."""
    from ..core.containers import configure_container
    from ..core.process import install_termination_handler

    # A terminated repgenr stops the tool it is running instead of leaving it
    # behind (SIGTERM from kill or a scheduler; SIGHUP when the terminal closes).
    install_termination_handler()
    _RUN_STATE["force"] = force
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        env = os.environ.get("REPGENR_LOG_LEVEL")
        level = getattr(logging, env.upper(), logging.INFO) if env else logging.INFO
    _RUN_STATE["log_level"] = level
    try:
        configure_container(
            backend=container,
            engine=container_engine,
            platform=platform,
            cache_dir=container_cache,
            wave_enabled=wave,
        )
    except UserInputError as exc:
        raise typer.BadParameter(str(exc), param_hint="'--container'") from exc
    _warn_ineffective_container_options(
        container, container_engine, container_cache, platform, wave
    )
    _configure_bin_dirs(bin_dir or [], container, wave)


def _configure_bin_dirs(entries: list[str], backend: str, wave: bool) -> None:
    """Set the per-tool directories from REPGENR_BIN_DIRS and --bin-dir (which wins)."""
    from ..core import bindirs

    raw = os.environ.get(bindirs.ENV_VAR, "")
    if not entries and not raw.strip():
        bindirs.configure_bin_dirs({})
        return
    resolved = {}
    # Each source is checked on its own, so an error names where the bad
    # entry came from: the variable or the option.
    for source, items in ((bindirs.ENV_VAR, raw.split(",")), ("--bin-dir", entries)):
        try:
            resolved.update(bindirs.validate(bindirs.parse_entries(items, source), source))
        except UserInputError as exc:
            raise typer.BadParameter(str(exc), param_hint=f"'{source}'") from exc
    bindirs.configure_bin_dirs(resolved)
    if backend != "none":
        for name in bindirs.ineffective_under_backend(wave_enabled=wave):
            typer.echo(
                f"WARNING --bin-dir {name} has no effect: under --container {backend} "
                f"'{name}' runs in its image.",
                err=True,
            )


def resolve_threads(ctx: typer.Context, param: typer.CallbackParam, value: int) -> int:
    """Callback of every ``-t/--threads``: follow the CPU limit of the process.

    Without ``-t`` the count is :data:`DEFAULT_THREADS`, lowered to the CPUs
    this process may use (affinity mask and cgroup CPU quota), so a container
    or pod limited to a few CPUs does not run 16 threads on them. An explicit
    value is kept, with a warning when it exceeds that limit. Logging is not
    configured yet when the callback runs, so this writes to stderr directly.
    """
    from ..core.resources import usable_cpus

    if ctx.resilient_parsing or param.name is None:
        return value
    limit = usable_cpus()
    # Compared by name: Typer vendors its own copy of Click's ParameterSource.
    source = ctx.get_parameter_source(param.name)
    if source is None or source.name in {"DEFAULT", "DEFAULT_MAP"}:
        if limit < value:
            if _RUN_STATE["log_level"] <= logging.INFO:
                typer.echo(
                    f"INFO Using {limit} threads, the CPU limit of this process "
                    f"(default {value}; -t/--threads sets the count).",
                    err=True,
                )
            return limit
        return value
    if value > limit:
        typer.echo(
            f"WARNING -t/--threads {value} exceeds the {limit} CPU(s) this process may "
            "use; the tools then share them.",
            err=True,
        )
    return value


def _warn_ineffective_container_options(
    backend: str, engine: str | None, cache: str | None, platform: str | None, wave: bool
) -> None:
    """Name container options (or their REPGENR_* variables) that change nothing.

    Logging is not configured yet when the callback runs, so this writes to
    stderr directly.
    """
    if backend == "none":
        given = [
            flag
            for flag, value in (
                ("--container-engine", engine),
                ("--container-cache", cache),
                ("--platform", platform),
                ("--wave", wave),
            )
            if value
        ]
        for flag in given:
            typer.echo(
                f"WARNING {flag} has no effect without --container docker or singularity.",
                err=True,
            )
    elif backend == "docker" and cache:
        typer.echo(
            "WARNING --container-cache is used only by --container singularity; "
            "Docker keeps images in its own storage.",
            err=True,
        )


def _tool_exit_code(returncode: int) -> int:
    """Exit code for a failed external tool.

    Interactive use exits with :attr:`ToolExecutionError.exit_code`. Under
    ``REPGENR_PROPAGATE_TOOL_EXIT=1`` (set by the Nextflow modules) the tool's
    code is forwarded, with a signal kill mapped to 128+signum (SIGKILL -> 137),
    so Nextflow's retry-on-exitStatus rule can react to e.g. an OOM kill.
    """
    if os.environ.get("REPGENR_PROPAGATE_TOOL_EXIT", "") in ("", "0"):
        return ToolExecutionError.exit_code
    if returncode < 0:
        code = 128 - returncode
    else:
        code = returncode
    return code if 0 < code <= 255 else ToolExecutionError.exit_code


@contextmanager
def stage_errors(logger: logging.Logger) -> Iterator[None]:
    """Turn errors into clean CLI exits instead of raw tracebacks.

    A :class:`RepGenRError` (expected, user-facing) is logged concisely. Any
    other exception is unexpected: a concise message goes to the console and the
    full traceback is captured -- in the run log when one exists (DEBUG), or on
    the console otherwise (e.g. a data-channel step with no workdir). Both exit
    non-zero. ``--verbose`` shows the traceback on the console too.
    """
    try:
        yield
    except typer.Exit:
        raise
    except ToolExecutionError as exc:
        log_files = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
        if log_files:
            log_name = Path(log_files[0].baseFilename).name
            logger.error("%s; see the run log (%s) for the command and output", exc, log_name)
            logger.debug("%s", exc.details())
        else:
            # No persistent log (data-channel step): keep the tail visible.
            logger.error("%s", exc)
            logger.error("%s", exc.details())
        raise typer.Exit(code=_tool_exit_code(exc.returncode)) from exc
    except RepGenRError as exc:
        logger.error("%s", exc)
        raise typer.Exit(code=exc.exit_code) from exc
    except Exception as exc:
        logger.error("Unexpected error: %s", exc)
        if any(isinstance(h, logging.FileHandler) for h in logger.handlers):
            logger.debug("Full traceback:", exc_info=True)
            logger.error("See the run log for the full traceback, or re-run with --verbose.")
        else:
            # No persistent log (data-channel step): surface the traceback now.
            logger.error("Full traceback:", exc_info=True)
        raise typer.Exit(code=1) from exc


def missing_workdir_message(workdir: Path) -> str:
    """One sentence for a -wd that does not exist, shared by stages and queries."""
    return (
        f"Workdir not found: {workdir}. Create it with an entry stage "
        "(metadata, ingest, reads or vmetadata) first."
    )


def require_existing_workdir(workdir: Path) -> None:
    """Exit 3 (WorkdirError) when a read-only query names a missing -wd.

    `status` and `doctor` run outside the stage harness, so they check here;
    an existing directory without repgenr.yaml is reported by the command itself.
    """
    if not workdir.is_dir():
        typer.echo(missing_workdir_message(workdir), err=True)
        raise typer.Exit(code=WorkdirError.exit_code)


def _run(stage_name: str, workdir: Path, build_params, *, create: bool = False) -> None:
    """Common harness: context, dispatch, clean error handling.

    Resume: a stage that already completed with the same parameters, the same
    inputs (per STAGE_INPUTS digests), and the same container identity is
    skipped, unless ``--force`` is set or one of its STAGE_DELIVERABLES is
    missing. Re-running an upstream stage changes a
    downstream stage's input digests, so the downstream stage reruns
    automatically. A stage that crashed before recording completion has no
    ``completed`` stamp and so always re-runs.
    """
    existed = workdir.exists()
    logger = configure_logging(workdir if existed else None, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        # Parameters are built and validated before an entry stage creates
        # its workdir, so a rejected invocation leaves no directory or log.
        params = build_params()
    if create and not existed:
        logger = configure_logging(workdir, level=_RUN_STATE["log_level"])
    with stage_errors(logger):
        ctx = WorkdirContext(workdir, logger=logger, create=create)
        try:
            _run_stage(stage_name, ctx, params, logger)
        finally:
            ctx.close()


def _run_stage(stage_name: str, ctx: WorkdirContext, params, logger) -> None:
    if not ctx.workdir.is_dir():
        # Only entry stages (create=True) start a workdir; any other stage
        # would otherwise create it as a side effect of opening the manifest.
        raise WorkdirError(missing_workdir_message(ctx.workdir))
    # Stages that cache an intermediate of their own (phylo's MSA) must not
    # reuse it under --force, which means "recompute this stage".
    ctx.force = bool(_RUN_STATE["force"])
    if _is_query_only(stage_name, params):
        # Pure query (list/preview): run the body, leave the resume record
        # of the last real run untouched.
        module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
        module.run(ctx, params)
        return
    # Digested once: upstream inputs are stable while this stage executes,
    # so the same digests are stamped onto the record after the run.
    digests = _stage_input_digests(ctx, stage_name, params)
    fingerprint = _stage_fingerprint(stage_name, params, digests, _env_fragment())
    prior = ctx.config.stages.get(stage_name)
    # Why a completed stage re-runs; logged after the precheck, so that a
    # refused run is not announced as re-running first.
    rerun_notes: list[tuple[Any, ...]] = []

    def rerun_note(*args: Any) -> None:
        rerun_notes.append(args)

    if not _RUN_STATE["force"] and prior is not None and prior.completed:
        if prior.fingerprint == fingerprint:
            missing = missing_deliverables(ctx, stage_name, params)
            if not missing:
                logger.info(
                    "Stage '%s' already completed with the same parameters and "
                    "inputs; skipping (use --force to re-run).",
                    stage_name,
                )
                if getattr(params, "sketch", None) is True:
                    # --sketch is not in the fingerprint, so a skip writes none.
                    logger.info(
                        "A skipped stage writes no sketches; 'repgenr sketch -wd %s' "
                        "writes the missing ones.",
                        ctx.workdir,
                    )
                return
            # The fingerprint excludes outputs, so check them separately: a
            # deliverable deleted by hand must be rebuilt, not skipped.
            # Per-file deliverables (genomes, representatives) can number in
            # the thousands; name a few.
            for path in missing[:_LOGGED_MISSING]:
                rerun_note(
                    "Stage '%s': deliverable %s missing; re-running.",
                    stage_name,
                    deliverable_label(ctx.workdir, path),
                )
            if len(missing) > _LOGGED_MISSING:
                rerun_note(
                    "Stage '%s': %d more deliverable(s) missing.",
                    stage_name,
                    len(missing) - _LOGGED_MISSING,
                )
        # A key present on one side only means the stage now reads a different
        # set of inputs (a flag such as --include-dereplicated, or an outgroup
        # added), not that a file's content changed; say which.
        shared = prior.inputs.keys() & digests.keys()
        changed = sorted(key for key in shared if prior.inputs[key] != digests[key])
        added = sorted(digests.keys() - prior.inputs.keys())
        dropped = sorted(prior.inputs.keys() - digests.keys())
        if changed and prior.inputs:
            rerun_note(
                "Stage '%s': input %s changed since last completion; re-running.",
                stage_name,
                ", ".join(f"'{c}'" for c in changed),
            )
        if (added or dropped) and prior.inputs:
            rerun_note(
                "Stage '%s': reads a different input set than at last completion "
                "(added: %s; no longer read: %s); re-running.",
                stage_name,
                ", ".join(f"'{c}'" for c in added) or "none",
                ", ".join(f"'{c}'" for c in dropped) or "none",
            )
    # The stage will run. A precheck guards the run, not the skip: it refuses
    # here, before the record is dirtied or a provisional record is written,
    # so a refusal leaves the last finished record (or no record) as it was,
    # and a stage that is skipped above is never refused or slowed by it.
    precheck = _STAGE_PRECHECKS.get(stage_name)
    if precheck is not None:
        precheck(ctx, params)
    for args in rerun_notes:
        logger.info(*args)
    finished = (
        (prior.completed, prior.fingerprint) if prior is not None and prior.completed else None
    )
    if prior is not None and prior.completed:
        # Dirty the record before the stage body runs: a crash mid-stage
        # must not leave a completed-looking record over partial outputs.
        prior.completed = None
        prior.fingerprint = None
        ctx.save_config()
    provisional = None
    if prior is None:
        # First run: write a provisional record (no completed stamp) so a stage
        # that fails or is killed mid-run shows as interrupted in `status` and
        # `doctor` instead of as not yet started. The stage's own record on
        # success replaces it.
        provisional = ctx.config.record_stage(
            stage_name,
            tool=_provisional_tool(params),
            params=_provisional_params(params),
        )
        ctx.save_config()
    # The incomplete record now on file (provisional, dirtied, or left by an
    # earlier failure). A stage writes a new record object when it finishes.
    pending = ctx.config.stages.get(stage_name)
    before = (
        _deliverables_state(ctx, stage_name, params)
        if provisional is not None or finished is not None
        else None
    )
    module = __import__(f"repgenr.stages.{stage_name}", fromlist=["run"])
    from ..core import bindirs

    bindirs.reset_used()
    try:
        module.run(ctx, params)
    except (UserInputError, WorkdirError, MissingBinaryError):
        # A refusal that left every deliverable as it was (phylo with too few
        # genomes, a tool missing at preflight) did not start: put the record
        # back as it was, so `status` and `doctor` do not report a stage that
        # never ran as interrupted -- no record on a first run, the last
        # finished one (whose outputs are untouched) on a re-run. A tool
        # failure (exit 6) or a crash keeps the incomplete record.
        if (
            before is not None
            and ctx.config.stages.get(stage_name) is pending
            and _deliverables_state(ctx, stage_name, params) == before
        ):
            if provisional is not None:
                del ctx.config.stages[stage_name]
            elif pending is not None and finished is not None:
                pending.completed, pending.fingerprint = finished
            ctx.save_config()
        raise
    if pending is not None and ctx.config.stages.get(stage_name) is pending:
        # The stage finished without a record of its own: leave no incomplete
        # record behind, or `status` would report it as interrupted.
        del ctx.config.stages[stage_name]
        ctx.save_config()
        return
    # Stamp fingerprint + input digests on the record the stage just wrote,
    # so the next invocation can skip.
    record = ctx.config.stages.get(stage_name)
    redigest = _REDIGEST_AFTER_RUN.get(stage_name)
    if record is not None and redigest is not None and redigest(params):
        digests = _stage_input_digests(ctx, stage_name, params)
        fingerprint = _stage_fingerprint(stage_name, params, digests, _env_fragment())
    if record is not None:
        record.fingerprint = fingerprint
        record.inputs = digests
        record.bin_dirs = bindirs.used()
        ctx.save_config()


def _provisional_params(params: object) -> dict[str, Any]:
    """Parameters of a stage invocation as plain YAML-safe values.

    Paths and other non-JSON values are stringified, as in the fingerprint.
    """
    if dataclasses.is_dataclass(params) and not isinstance(params, type):
        payload: dict = dataclasses.asdict(params)
    else:
        payload = dict(vars(params))
    return json.loads(json.dumps(payload, default=str))


def _provisional_tool(params: object) -> str | None:
    tool = getattr(params, "tool", None)
    return tool if isinstance(tool, str) else None


def gated_extra(registry, tool: str, key: str, value: object, *, flag: str | None = None) -> dict:
    """Return ``{key: value}`` only when ``tool`` reads that extra.

    Injecting a key a tool ignores would change the resume fingerprint without
    changing the result. ``auto`` passes the key through; the stage warns after
    it has picked a concrete tool. ``flag`` names the option the user gave
    (for example ``--virus``): when set and the tool does not read the key, a
    warning says the option has no effect and which tools read it. Callers
    that inject the key themselves (``run --viral``) leave it unset.
    """
    if tool != "auto":
        caps = registry.get(tool).capabilities
        if key not in caps.accepted_extras:
            if flag is not None:
                readers = sorted(
                    name
                    for name in registry.names()
                    if not registry.is_broken(name)
                    and key in registry.get(name).capabilities.accepted_extras
                )
                logging.getLogger("repgenr").warning(
                    "%s has no effect with --tool %s; it is read by: %s.",
                    flag,
                    tool,
                    ", ".join(readers) or "no installed tool",
                )
            return {}
    return {key: value}


def _parse_key_values(items: list[str], label: str) -> dict[str, str]:
    """Parse repeated ``key=value`` options into a dict (used for tool extras)."""
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise UserInputError(f"{label} must be key=value, got '{item}'.")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise UserInputError(f"{label} has an empty key in '{item}'.")
        out[key] = value.strip()
    return out


def _read_path_fofn(path: Path) -> list[Path]:
    """Read a file-of-filenames (one path per line; blank lines ignored)."""
    if not path.exists():
        raise UserInputError(f"File not found: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    return [Path(line.strip()) for line in lines if line.strip()]
