# Developing RepGenR2

This page is for contributors who change the code or add a tool adapter. It
describes how the package is organised, the files that pass between stages,
how the ten tool families are discovered, and the steps and tests a new
adapter needs. Users who only run the pipeline need [usage.md](usage.md)
instead.

Contents:

1. [Architecture](#1-architecture)
2. [Data contracts](#2-data-contracts)
3. [Tool families and the registry](#3-tool-families-and-the-registry)
4. [Adding an adapter](#4-adding-an-adapter)
5. [Testing](#5-testing)
6. [Conventions](#6-conventions)

## 1. Architecture

### Package layout

RepGenR2 is an importable Python package (`src/repgenr`, Python 3.12 or
later). The `repgenr` console script targets the Typer app defined in
`repgenr.cli.base`; `repgenr.cli.main` imports the `cmd_*` modules, which
register their commands on that app. Each command builds a parameter
dataclass and calls an importable function, `repgenr.stages.<name>.run(ctx,
params)` for a stage.

```
cli/
  main.py              entry point; imports the cmd_* modules
  base.py              Typer app, stage harness (_run), resume fingerprint,
                       STAGE_INPUTS / STAGE_DELIVERABLES, COMMAND_PANELS,
                       registry-driven help
  param_builders.py    CLI options -> stage parameter dataclasses
  cmd_*.py             one module per command group (bacterial, viral, ingest,
                       reads, phylo, run, steps, misc)
core/
  context.py           WorkdirContext: paths + services
  config.py            repgenr.yaml (per-stage tool, params, versions,
                       resume fingerprint + input digests)
  manifest.py          SQLite genome inventory
  contracts.py         canonical inter-stage file readers/writers
  inputs.py            resume input digests
  integrity.py         completeness guards
  process.py           subprocess wrapper (no shell, no globs)
  containers.py        docker/singularity/wave backend for tool commands
  binaries.py          binary presence + version preflight
  bindirs.py           per-tool host directories (--bin-dir, REPGENR_BIN_DIRS)
  versions.py          versions.yml writer
  plugins.py           ToolCapabilities + entry-point Registry
  doctor.py            working-directory health check
  sketches.py          genome sketches (sourmash signatures per genome)
  sourmash.py          the shared sourmash tool specification
  annotations.py       annotation cache and its index
  executors.py         in-stage thread pool
  resources.py         CPUs available to the process (affinity, cgroup quota)
  http.py, ena.py, ncbi_taxonomy.py   network clients
  errors.py, logging.py
stages/                one module per stage (metadata, genome, ingest, reads,
                       assemble, dereplicate, annotate, snptype, phylo,
                       tree2tax, glance, sketch, ...) and the stateless step
                       implementations (*_steps.py)
viral/                 NCBI Virus and BV-BRC back-ends of vmetadata/vgenome
dereplicators/  aligners/  snptypers/  maskers/  treebuilders/
annotators/  coregenes/  assemblers/  polishers/  classifiers/
                       one package per tool family: base.py (ABC, parameter
                       and result dataclasses, registry) + one module per adapter
converters/            XMFA/MAF/HAL/GFA -> MSA-FASTA
tree/                  Newick, neighbour joining (rapidNJ when installed),
                       rooting, medoid reference choice, alignment statistics
```

### Stages and steps

The commands fall into two kinds:

* A **stage** writes into a working directory (`-wd`), the directory that
  holds one run's genomes, intermediate files and results. Stages share the
  directory and run through the stage harness `_run` in `cli/base.py`. A
  `WorkdirContext` resolves the canonical path layout, logging, the
  `repgenr.yaml` configuration and provenance file, and the SQLite genome
  manifest. Stages record what they did and can be resumed.
* A **step** is a stateless data-channel command (`genome-fetch`,
  `dereplicate-chunk`, `dereplicate-merge`, `annotate-build`,
  `annotate-index`, `phylo-build`, `tree2tax-relations`, `assemble-run`,
  `genome-qc`, `reads-gather`; see `COMMAND_PANELS` in `cli/base.py`). A step
  reads explicit inputs (a `selection.tsv` or a file of file names) and writes
  a result directory. It does not use a working directory or the stage
  harness. The Nextflow processes call the steps.

Defaults for stage parameters live on the parameter dataclasses and are built
through `cli/param_builders.py`, which both the manual commands and
`repgenr run` share.

### The working directory, configuration and manifest

* `repgenr.yaml` (`core/config.py`) records, per stage, the tool chosen, the
  parameters used, the resolved tool versions, and the resume fingerprint with
  its input digests.
* The manifest (`core/manifest.py`, SQLite, schema version 4) records which
  genomes are selected, where their files are, their taxonomy, their
  dereplication status and their sketches. Genome state lives here, not in
  `str(dict)` or `pickle` blobs, and stages do not rescan directories.

How a stage decides to skip on resume is described under
[Resume semantics](#resume-semantics).

### The Nextflow layer

The pipeline under `nextflow/` runs in one of three modes (`--mode
bacterial`, `viral` or `reads`). It provides the parallelism across genomes
and the HPC and cloud execution; RepGenR2 itself parallelises only within a
stage (`core/executors.py`).

* Every channel carries a run-level meta map (`id`, `mode`) built by
  `run_meta()` (`nextflow/subworkflows/local/run_meta.nf`), which `main.nf`
  calls. Processes exchange `tuple val(meta), path(...)`.
* Processes read tool flags only from `task.ext.args`, which
  `nextflow/conf/modules.config` maps from the pipeline parameters, and they
  are published from configuration rather than from the module files.
* Most modules call a stateless step. The acquisition modules run a stage in
  a task-local working directory (`metadata`, `vmetadata`, `vgenome` and `reads`).
* Resource labels in `nextflow/conf/base.config` (`process_low`,
  `process_medium`, `process_high`, `process_assembly`) place heavy processes,
  such as `PHYLO`, which may run Cactus, on large nodes.
* The pipeline schema does not restrict tool names, so an adapter installed
  as a plugin is selectable from Nextflow as well.

### Scaling

* `process.run` always passes argument vectors and uses files of file names
  instead of shell globs, so large genome sets do not reach `ARG_MAX`. A tool
  without a list-file option gets a staged directory where it reads one
  (skDER, ParSNP). The few that take every path on argv (SibeliaZ,
  snippy-core) call `warn_argv_bytes`. `process.run` measures the final argv
  (container preamble included) plus the environment against
  `os.sysconf("SC_ARG_MAX")` and, on Linux, the 128 KiB single-argument
  limit, and raises `UserInputError` before the exec instead of failing with
  E2BIG.
* Two-stage chunked dereplication is opt-in via `--process-size` for any
  tool. `supports_native_scaling` only changes the log wording. Automatic
  selection keys on `recommended_max_genomes`.

### Phylogenetics: independent choices

The phylo stage combines independent choices:

* genome set: dereplicated representatives or all genomes;
* MSA (multiple sequence alignment) source: a whole-genome aligner, a SNP
  typer's core-SNP alignment, or a concatenated core-gene alignment built
  from the `annotate` stage's annotations (`--msa-source aligner`, `snptype`
  or `coregene`); skipped for alignment-free tree builders;
* tree builder: MSA-based (iqtree, fasttree, veryfasttree, raxmlng) or
  alignment-free (mashtree, sourmash, skani, ska2).

The stage adds the outgroup genome (a genome outside the target taxon, used
to root the tree) to the input set and passes it to the builders that can
root (iqtree, raxmlng). mashtree, fasttree, veryfasttree, sourmash, skani and
ska2 emit an unrooted tree. `phylo --root` then roots the published tree (on
the outgroup edge by default, at the midpoint, or not at all) and keeps the
builder's output as `tree/tree.unrooted.nwk`. tree2tax reroots on the
outgroup only when the tree is not already rooted there, and keeps a midpoint
root as it is. An outgroup that is not a leaf of the tree is an error
(exit 3). A tree that `phylo` built without an outgroup is left unrooted.

## 2. Data contracts

The canonical files that pass between stages are read and written by
`core/contracts.py`. The main ones are listed below; [output.md](output.md)
describes every file a stage writes.

```
selection.tsv                  the selected genomes and their metadata
derep/representatives/         representative genome FASTAs
derep/clusters.tsv             representative <TAB> member
derep/genome_status.tsv        genome <TAB> status
derep/cluster_summary.tsv      one row per representative: size, species, quality
annotations/annotations.tsv    index of the annotation cache (annotate)
align/msa.fasta                aligner output (MSA-FASTA)
snp/core_snp.fasta             SNP typer output (core-SNP alignment, snptype stage)
tree/msa/core_snp.fasta        core-SNP alignment from the typing pass of phylo
tree/msa/core_gene_alignment.fasta   core-gene alignment (--msa-source coregene)
tree/tree.nwk                  Newick tree
tree/tree.unrooted.nwk         the builder's tree, when phylo rerooted it
tree2tax.tsv                   child <TAB> parent (FlexTaxD)
genomes_map.tsv                accession <TAB> leaf
```

An adapter never writes these files. It returns a normalised result
dataclass, and the owning stage writes the contract files. This keeps the
tools within a family interchangeable with no downstream change.

### Record names

A genome has one name in alignments, trees and tables:
`core.contracts.record_name(path)`, its file name without the FASTA suffix
and `.gz` (`x.fasta.gz` gives `x`; `Path.stem` would give `x.fasta`). An
aligner, SNP typer or tree builder adapter names its records and leaves with
it, and maps back to it where the tool rewrites names (ParSNP's `.ref`,
cactus's `.` to `_`). The phylo stage checks tree leaves against the same
names and rebuilds a reused alignment whose records differ from them.

### Genome sketches

Genome sketches (`core/sketches.py`) are a contract written by the
genome-writing stages and by `repgenr sketch`: `sketches/<record name>.sig.zip`
per genome, holding three DNA signatures (k=21, 31, 51; scaled=1000). The
sketch file, parameter string and FASTA SHA-256 are recorded in the manifest
or, for an outgroup without a manifest row, in `sketches/outgroup.json`.
Staleness is decided from those records, not from file times.

* The writers call `sketch_stage_genomes` after the manifest is updated and
  before the stage record is written. Sketches are not in the writers' entries
  of `STAGE_DELIVERABLES`, and `sketch` is excluded from the resume
  fingerprint.
* Every sourmash caller (the dereplicator, the tree builder, the classifier,
  the sketch step and the reads sketch of `assemble`,
  `core.sketches.sketch_reads`) derives its `ToolCapabilities` from the one
  specification in `core/sourmash.py` (`SOURMASH_TOOL`,
  `sourmash_capabilities(...)`). The container image, conda spec and binary
  check therefore live in one place. The Nextflow `SKETCH` module
  (`nextflow/modules/local/dataflow/sketch.nf`) names the same image, which a
  unit test checks.
* The adapters that read the sketches declare it through `sketch_request`
  (see [Optional hooks](#optional-hooks)).

## 3. Tool families and the registry

There are ten tool families. Each is an abstract base class (ABC) plus a
`Registry` bound to an entry-point group. The table lists them in the order a
run uses them, grouped by what they work on: read tools turn sequencing runs
into genome assemblies, genome tools reduce and annotate a genome set,
alignment tools produce the multiple sequence alignment, and tree tools build
the tree. Data sources (GTDB, NCBI Datasets, NCBI Virus, BV-BRC, ENA) are
built into the stages and are not pluggable.

| Group | Family | ABC | Method | Returns | Entry-point group |
|---|---|---|---|---|---|
| Read tools | Assembler | `repgenr.assemblers.base.Assembler` | `assemble(reads, out_dir, params, logger)` | `AssemblyResult` | `repgenr.assemblers` |
| Read tools | Polisher | `repgenr.polishers.base.Polisher` | `polish(reads, draft, out_dir, params, logger)` | `PolishResult` | `repgenr.polishers` |
| Read tools | Classifier | `repgenr.classifiers.base.Classifier` | `classify(genomes, out_dir, params, logger)` | `dict[str, Classification]` | `repgenr.classifiers` |
| Genome tools | Dereplicator | `repgenr.dereplicators.base.Dereplicator` | `dereplicate(genomes, out_dir, params, logger)` | `DerepResult` | `repgenr.dereplicators` |
| Genome tools | Annotator | `repgenr.annotators.base.Annotator` | `call(genome, record_name, out_dir, params, logger)` | `AnnotateResult` | `repgenr.annotators` |
| Alignment tools | Aligner | `repgenr.aligners.base.Aligner` | `align(genomes, reference, out_dir, params, logger)` | `AlignResult` | `repgenr.aligners` |
| Alignment tools | SNP typer | `repgenr.snptypers.base.SnpTyper` | `call(genomes, reference, out_dir, params, logger)` | `SnpResult` | `repgenr.snptypers` |
| Alignment tools | Masker | `repgenr.maskers.base.Masker` | `mask(full_alignment, out_dir, params: MaskParams, logger)` | masked variable-site FASTA `Path` | `repgenr.maskers` |
| Alignment tools | Core-gene builder | `repgenr.coregenes.base.CoreGeneBuilder` | `call(gffs, out_dir, params, logger)` | `CoreGeneResult` | `repgenr.coregenes` |
| Tree tools | Tree builder | `repgenr.treebuilders.base.TreeBuilder` | `build(msa_or_genomes, out_dir, params, logger)` | Newick `Path` | `repgenr.treebuilders` |

The parameter dataclasses are, in the same order, `AssembleParams`,
`PolishParams`, `ClassifyParams`, `DerepParams`, `AnnotateParams`,
`AlignParams`, `SnpParams`, `MaskParams`, `CoreGeneParams` and `TreeParams`,
each in the family's `base.py`.

In-tree and third-party adapters are discovered the same way, through
package entry points, so the core never imports a concrete adapter. A third
party can ship adapters in its own package without modifying RepGenR2: after
`pip install`, the tool appears in `repgenr list-tools` and in the
registry-generated `--help` texts, and it is selectable by name, including
from the Nextflow pipeline. Programmatic registration (tests, embedding
applications) uses the public API: `registry.register("mytool",
MyDereplicator)` and `registry.unregister("mytool")`, or the `register_tool`
pytest fixture.

### Class attributes the stages read

Set these per-family attributes, or inherit the default deliberately:

* `TreeBuilder.input_kind`: `InputKind.GENOMES` for alignment-free builders
  (they receive genome files), else the default `InputKind.MSA_FASTA` (they
  receive one MSA path). Narrow the argument with `as_genome_list` /
  `as_msa_path` from `treebuilders.base`.
* `TreeBuilder.supports_without_bootstrap`: False for a builder whose tree
  has no support values unless a bootstrap is requested (IQ-TREE); the stage
  warns only where that matters.
* `SnpTyper.requires_reference`: when True (the default) the stage resolves a
  reference genome before calling the adapter.
* `SnpTyper.produces_full_alignment`: whether `call` returns a whole-genome
  alignment (`SnpResult.full_alignment`). A `Masker` reads the
  **whole-genome** alignment, not the SNP alignment, and returns the masked
  variable-site FASTA that replaces the typer's core-SNP output. A `SnpTyper`
  should therefore set `SnpResult.full_alignment` whenever it can produce
  one. `--mask` with a typer that declares no full alignment is refused
  before any SNP calling runs; a typer that declares one but returns
  `full_alignment=None` is refused after calling.
* `Masker.tree_output()`: the masker's own tree from the last call, when it
  builds one (Gubbins does); the stage keeps it as `masker_tree.nwk`.
* `Annotator.requires_database`: True for a tool that needs a database
  directory (`AnnotateParams.db`, as for bakta).
* `CoreGeneBuilder.accepted_annotators`: annotator names the builder accepts;
  empty means any GFF3 with sequences.
* `Assembler.read_types` and `Assembler.layouts`, `Polisher.read_types`: the
  ENA instrument platforms and library layouts the tool accepts.
* `Classifier.needs_lineages`: whether `classify` needs
  `ClassifyParams.lineages`; checked before any assembly.

### Optional hooks

* `Dereplicator.compare(...)`: all-against-all similarity for `repgenr
  glance`. Return a `CompareResult` whose CSV has
  `genome1`/`genome2`/`similarity` columns, and set its `measure` to the
  label for the plot axes. An adapter with a dense similarity matrix can
  write both files with `dereplicators/compare_io.py`, as sourmash does. A
  compare-capable adapter ranks after dRep and sourmash for `glance --tool
  auto`.
* `TreeBuilder.distance_matrix(...)`: a pairwise matrix for the viral
  outgroup step.
* `sketch_request(extra)` on a `Dereplicator`, `TreeBuilder` or `Classifier`
  that compares sourmash sketches: return the `(ksize, scaled)` the adapter
  compares at (the default returns None: no sketches). When the
  working-directory sketches hold those parameters, the stage passes genome
  path -> `.sig.zip` in `DerepParams.sketches` (and the `sketches` keyword of
  `compare`), `TreeParams.sketches` or `ClassifyParams.sketches`. The adapter
  selects its k-mer size with `-k` and sketches any genome left out of the
  mapping itself. Each signature is named by the genome's `record_name`. The
  stages resolve the mapping through `core.sketches.adapter_sketches` with a
  provider: `workdir_provider` (writes missing and stale sketches into
  `sketches/` and records them) in the working-directory stages, and
  `directory_provider` (match by record name, no digest) in the stateless
  steps given `--sketches-dir`. The adapters stay unaware of the working
  directory.

A tool without an override of `compare` or `distance_matrix` is rejected
with a clear message when a user selects it for those features.

### ToolCapabilities

Every adapter declares a `ToolCapabilities` (`repgenr.core.plugins`):

```python
capabilities = ToolCapabilities(
    name="mytool",
    required_binaries=(BinarySpec("mytool", version_args=("--version",), min_version="1.2.0"),),
    default_params={"mode": "fast"},          # the adapter's own defaults (see below)
    recommended_max_genomes=2000,             # drives auto-selection and scale warnings
    container="quay.io/biocontainers/mytool:1.2.3--0",  # or:
    conda=("bioconda::mytool",),              # Wave-resolved image when no explicit one
    accepted_extras=frozenset({"mode"}),      # extra keys the adapter actually reads
    ignored_params=frozenset({"primary_ani"}),  # standard params the tool has no use for
    reads_gzip=True,                          # the tool reads gzipped FASTA itself
)
```

* `required_binaries` powers the preflight (presence and minimum version)
  when running natively. Each version query runs in its own process group and
  may take `core.binaries.VERSION_TIMEOUT` seconds (30; `list-tools --check`
  uses 8). A query that does not answer is stopped together with any helper
  it started, and a warning names the binary. A binary that does not answer,
  or prints no version number, is looked up in the conda package records of
  its environment (`conda-meta`). Without such a record its version is
  `unknown`, and a `strict_version` floor then fails with "no version
  reported". A version query is tracked like a running tool, so Ctrl-C or a
  termination signal to repgenr stops it and its helpers as well.
* `default_params` holds the adapter's own tunable defaults. The
  dereplication stage and steps, `glance`, `annotate` and the tree builder in
  `phylo` merge them under `params.extra`. The other callers (aligners, SNP
  typers, maskers, core-gene builders, assemblers, polishers, classifiers)
  pass the user's extras only. Read a default with
  `params.extra.get(key, self.capabilities.default_params[key])` (or
  `parse_extra_int`) rather than assuming the key is present.
* **Set `container` and/or `conda`**, or the tool cannot run under
  `--container`/`--wave` (the run falls back to the host with a warning).
* `recommended_max_genomes` feeds `--tool auto` (the tightest fitting limit
  wins) and the over-scale warnings. Unbounded (`None`) tools are picked when
  nothing bounded fits. The selection rule is described in
  [choosing-tools.md](choosing-tools.md#how-auto-chooses).
* `auto_eligible`: False for a tool that suits one kind of input only (the
  ska2 tree builder, clonal sets), so `auto` never picks it; it must be named.
* `accepted_extras` names the `params.extra` keys the adapter reads. Users
  set them with `--tool-arg key=value` (dereplicate, snptype, annotate,
  assemble and the data-channel steps) or `--aligner-arg key=value` (phylo).
  The stages warn, by name, about extra keys that no adapter they run
  declares in `accepted_extras`. phylo drives two adapters (the MSA source
  and the tree builder) from one dictionary, so it warns only about keys
  neither of them reads. The CLI does not inject stage-level flags (such as
  `--virus`) into the extras of a tool that does not read them. Parse
  integer extras with `repgenr.core.plugins.parse_extra_int`, so that a bad
  value raises a clean `UserInputError`.
* `ignored_params` names the fields of the family's parameter dataclass
  (`DerepParams`, `AlignParams`, `SnpParams`, `TreeParams`,
  `CoreGeneParams`) that the adapter never passes to its tool, for instance
  `primary_ani` on a tool with no pre-clustering step or `bootstrap` on a
  distance-based tree. The stage warns, by name, when a user sets one of them
  to a non-default value, so a request the tool cannot honour is not dropped
  without notice. Declare every standard parameter the adapter drops; leave
  the set empty when the adapter passes them all.
* `reads_gzip` declares that the tool reads gzipped FASTA genomes. It is
  False by default. `snptype`, `phylo` and `annotate` then give an aligner,
  SNP typer, genome-input tree builder or annotator decompressed copies of
  gzipped genomes (`core.process.stage_plain_inputs`), named by
  `record_name`. Set it to True only after checking the tool on a
  `.fasta.gz` genome: a tool that misreads gzip may exit 0 with an empty
  result (SibeliaZ does).

### Running the external tool

**Use `run_tool(self.capabilities, cmd, logger=logger, ...)` from
`repgenr.core.containers` for every external call.** Never call `process.run`
directly and never use `shell=True`. `run_tool` routes the command through the
active container backend (docker, singularity or Wave) using the adapter's
`container` / `conda` declaration, and runs it on the host otherwise.

* Use `write_fofn` (`repgenr.core.process`) instead of a shell glob for large
  input lists. Declare the directories that hold the listed files through
  `extra_mounts=`, since the contents of a file of file names are invisible to
  the container backend.
* Pass a result written to stdout through `stdout_path=`; the file is
  published atomically, and only on success.
* Variables the tool itself must see (a thread count such as
  `OMP_NUM_THREADS` or `RAYON_NUM_THREADS`) go in `tool_env=`, which reaches
  the tool in every backend: Docker gets them as `-e`, Singularity passes
  the client's environment on, the host runs the tool in it. `env=` is the
  environment of the process repgenr starts, which under Docker is the
  `docker run` client, not the tool.
* When one unit of work is several commands in a row, such as a per-genome
  mapping and calling chain, `run_chain(self.capabilities, [(prefix, argv),
  ...], logger=logger)` runs them in order. On the host that is a loop over
  `run_tool`. Under a container backend the whole sequence runs in one
  container, which matters when a chain has half a dozen commands and every
  container start costs time. Steps are plain argument vectors, so a tool
  that writes to stdout needs an output flag, and paths must be absolute to
  be mounted.

The return value of `preflight()` (host versions or the container image
reference) is recorded in the `repgenr.yaml` provenance and appears in
`repgenr versions` and the Nextflow `versions.yml` bridge.

### Resume semantics

A stage's skip decision digests its declared inputs (`STAGE_INPUTS` in
`repgenr/cli/base.py`), its parameters and the container identity. An adapter
therefore needs no resume logic of its own, but it must read only the inputs
the stage declares: a tool that reads undeclared files would resume
incorrectly.

* A declared directory is digested from its genome FASTA files only
  (`core.inputs.dir_stat_digest`, the files `list_fasta` returns). A stage
  that reads anything else from a directory (a table, an index) declares that
  file by its own path. derep-stock's unpack, for example, lists the stored
  run's `clusters.tsv`, `genome_status.tsv`, `record.json` and
  `representatives/` one by one.
* Outputs are not part of the fingerprint. A matching stage is skipped only
  when every path in its entry in `STAGE_DELIVERABLES` (same module) exists,
  a directory counting only when non-empty.
* A new stage adds an entry to both tables.
* A stage whose one record serves several invocations (derep-stock's named
  runs) can also register a check in `_STAGE_PRECHECKS`. The check runs
  before the harness marks the record incomplete, so a refused call leaves
  the last finished one intact.

## 4. Adding an adapter

A new tool is a self-contained adapter module plus one entry-point line. The
steps below are for a built-in (in-tree) adapter. A third-party package
follows steps 1 to 3 in its own repository, adds the entry point to its own
package metadata, and depends on the distribution `repgenr2`.

1. **Write the adapter module** in the family's package
   (`src/repgenr/<family>/mytool.py`). Subclass the family ABC and declare
   `capabilities` (see [ToolCapabilities](#toolcapabilities)) with:
   `container` and/or `conda`; a `BinarySpec.min_version` floor for each
   binary (`tests/unit/test_version_floors.py` requires one, except for tools
   that print no parseable version); and `ignored_params` for any standard
   parameter the tool has no use for. Return the family's result dataclass.
   Do not write `derep/...`, `align/...` or other contract files from the
   adapter. A dereplicator skeleton:

   ```python
   from repgenr.core.binaries import BinarySpec
   from repgenr.core.containers import run_tool
   from repgenr.core.plugins import ToolCapabilities
   from repgenr.dereplicators.base import Dereplicator, DerepResult

   class MyDereplicator(Dereplicator):
       capabilities = ToolCapabilities(
           name="mytool",
           required_binaries=(BinarySpec("mytool", version_args=("--version",)),),
           conda=("bioconda::mytool",),
       )

       def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
           run_tool(self.capabilities, ["mytool", ...], logger=logger, log_prefix="mytool")
           ...  # parse the tool's output
           return DerepResult(representatives=[...], clusters={...}, genome_status={...})
   ```

2. **Register the entry point** in `pyproject.toml`, under the family's group
   (`repgenr.dereplicators`, `repgenr.aligners`, `repgenr.snptypers`,
   `repgenr.maskers`, `repgenr.treebuilders`, `repgenr.annotators`,
   `repgenr.coregenes`, `repgenr.assemblers`, `repgenr.polishers` or
   `repgenr.classifiers`):

   ```toml
   [project.entry-points."repgenr.dereplicators"]
   mytool = "repgenr.dereplicators.mytool:MyDereplicator"
   ```

   Reinstall the package (`pip install -e ".[dev]"`) so that the new entry
   point is visible.

3. **Write the unit tests** (see
   [Unit tests for an adapter](#unit-tests-for-an-adapter)). Add the tool to
   the family's contract suite, including the literal set of built-in
   adapters that suite asserts against, and gate any test that needs the real
   binary with `@pytest.mark.requires_binary`. `test_entry_points_discovered`
   in `tests/unit/test_plugins.py` asserts that the built-in names of six
   families (dereplicators, aligners, SNP typers, tree builders, annotators,
   core-gene builders) are discovered; add the new name to its family's set
   there.

4. **Add the conda package** to `envs/core.yml`, with the same `>=` floor as
   its `BinarySpec.min_version`. `tests/test_env_files.py` checks that every
   preflight floor appears in exactly one env file and that every required
   binary and adapter conda spec is covered by one. When the package does not
   solve together with core, add a satellite file `envs/<tool>.yml`, add its
   name to `EXPECTED_FILES` in `tests/test_env_files.py`, and add it to the
   `env` list of the solve matrix in `.github/workflows/envs.yml` (the test
   checks that the matrix lists exactly those files). A tool distributed
   separately, as cactus is, needs a note instead (`CONTAINER_ONLY` in the
   same test).

5. **Update the documentation tables**:
   * a row in the limits table of [choosing-tools.md](choosing-tools.md#declared-limits)
     (header `| Family | Tool | Declared limit |`), which
     `tests/unit/test_docs_tool_limits.py` checks against the registries;
     update the row also when an adapter's limit changes;
   * the adapter table in `README.md` and, if user-facing behaviour differs,
     a note in [usage.md](usage.md);
   * a citation in `CITATIONS.md`;
   * for a dereplicator, the name in the `derep_tool` description of
     `nextflow/nextflow_schema.json`.

6. **Update the CLI matrix** when the tool adds a flag: a record in
   `tests/audit/cli_matrix.yaml` for every new flag, then
   `python scripts/render_cli_matrix.py` to regenerate
   `docs/audit/cli-matrix.md` and `docs/cli-reference.md`
   (`tests/unit/test_cli_matrix.py` checks both against the command tree). A
   new command must also be added to `COMMAND_PANELS` in
   `src/repgenr/cli/base.py`, which drives the `--help` grouping, the
   rendered reference and a test.

7. **Add a live test** under `tests/live/` so that
   [verification.md](verification.md) can report the adapter as verified
   against the real tool.

Everything else (CLI help, validation, `list-tools`, automatic selection,
containers, provenance, resume and Nextflow) picks the tool up from the
registry.

## 5. Testing

### Setting up and running the checks

```bash
pip install -e ".[dev]"
ruff check src/ tests/ scripts/ benchmarks/
ruff format --check src/ tests/ scripts/ benchmarks/
mypy src/repgenr
pytest -q
```

These are the commands CI runs (`.github/workflows/ci.yml`); CI adds
`--cov=repgenr --cov-report=term-missing --cov-fail-under=85` to pytest.

### pytest

`pytest -q` runs the offline suite (`tests/unit/`, `tests/integration/` and
`tests/test_env_files.py`). No external tools or network are needed: the
`live` marker is deselected by default (`addopts = "-m 'not live'"` in
`pyproject.toml`), and a test marked `@pytest.mark.requires_binary("name")`
is skipped when the binary is not on `PATH` (`tests/conftest.py`). Some tests
read the documentation or the changelog: `test_docs_tool_limits.py`,
`test_cli_matrix.py`, `test_version_consistency.py` (CHANGELOG.md),
`test_genome_stage.py` and
`test_vmetadata_stage.py` under `tests/unit/`.

#### Unit tests for an adapter

Unit-test an adapter by patching the adapter module's `run_tool` (or
`run_chain`) with a fake runner that records the exact argument vector and
writes canned tool output, then assert the returned dataclass. For an
in-tree adapter, add the tool to the family's parametrised contract suite:

| Family | Contract suite |
|---|---|
| dereplicators, tree builders | `tests/unit/test_adapter_contracts.py` |
| SNP typers, aligners | `tests/unit/test_snptyper_aligner_contracts.py` |
| assemblers | `tests/unit/test_assembler_contracts.py` |
| polishers | `tests/unit/test_polisher_contracts.py` |

That means one entry in the family's token dictionary, a canned-output branch
in the fake runner, and the name in the literal set of built-in adapters the
suite asserts against (the suites do not discover new entry points by
themselves). The other families have their own test modules (for example
`test_annotators.py`, `test_coregenes.py`, `test_gubbins_masker.py`).
Integration-test the owning stage with the `register_tool` fixture
(`tests/conftest.py`), which registers an adapter class for one test and
removes it afterwards.

### nf-test stub tests

The Nextflow tests live in `nextflow/tests/` and run under the `test`
profile (`nf-test.config`). Tests tagged `stub` use stub blocks, so they need
no external tools or network; they check channel wiring and the file
contracts. Run them before pushing a Nextflow change:

```bash
nf-test test --tag stub
```

### Real-tool tests

* `tests/live/`: pytest tests marked `live` that run the installed `repgenr`
  console script against real tools. They are never run in CI.
  [verification.md](verification.md#re-running-the-checks-the-live-suite)
  describes the configuration file and the `network` and `container`
  markers. The offline subset runs with:

  ```bash
  pytest tests/live -m "live and not network and not container" \
      --live-config tests/live/live.local.toml -ra
  ```

* Tests marked `requires_binary` in the offline suite run when their binary
  is on `PATH`.
* The nf-test tests tagged `e2e` (`nextflow/tests/local_dataflow.nf.test`) run
  the data-channel pipeline on a small bundled genome set with real tools.

### CI

`.github/workflows/ci.yml` runs on every pull request and on pushes to
`main`:

| Job | What it runs |
|---|---|
| `lint-and-test` | ruff check, ruff format check, mypy and pytest with coverage (gate 85%), on Python 3.12 and 3.13 |
| `nf-test` | `nf-test test --tag stub` (Nextflow 26.04.3, nf-test 0.9.5) |
| `nf-test-e2e` | `nf-test test --tag e2e` in a conda environment with sourmash, mashtree and dendropy |

`.github/workflows/envs.yml` dry-run solves every `envs/*.yml` file on
linux-64 (and core also on osx-arm64) when the env files change, and weekly.

## 6. Conventions

* External tools run only through `run_tool` / `run_chain`, with argument
  vectors and files of file names; no `shell=True` and no shell globs.
* An adapter returns a normalised dataclass; only the owning stage writes the
  contract files.
* Genomes are named by `record_name` in every alignment, tree and table.
* An adapter reads only the inputs its stage declares in `STAGE_INPUTS`.
* New stage options go into `cli/param_builders.py`, so that the manual
  commands and `repgenr run` share them; never into one entry point only.
* Every flag has a record in `tests/audit/cli_matrix.yaml`. The generated
  files (`docs/cli-reference.md`, `docs/audit/cli-matrix.md`) are rendered,
  not edited by hand.
* Declared limits, version floors and env files agree with the adapters; the
  tests named in [Adding an adapter](#4-adding-an-adapter) check this.
