# Developing RepGenR2

This page is for contributors who change the code or add a tool adapter. It
describes how the package is organised, the files that pass between stages,
how the ten tool families are discovered, and the steps and tests a new
adapter needs. Users who only run the pipeline need [usage.md](usage.md)
and the pages it links to instead.

Contents:

1. [Architecture](#1-architecture)
2. [Data contracts](#2-data-contracts)
3. [Tool families and the registry](#3-tool-families-and-the-registry)
4. [Testing](#4-testing)
5. [Conventions](#5-conventions)

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
                       tree2tax, glance, sketch, ...) and the step
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
* A **step** is a stateless command (`genome-fetch`, `dereplicate-chunk`,
  `dereplicate-merge`, `dereplicate-refine`, `annotate-build`,
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
* Most modules call a step. The acquisition modules run a stage in a
  task-local working directory (`metadata`, `vmetadata`, `vgenome` and
  `reads`).
* Resource labels in `nextflow/conf/base.config` (`process_low`,
  `process_medium`, `process_high`, `process_assembly`) place heavy processes,
  such as `PHYLO`, which may run Cactus, on large nodes.
* The pipeline schema does not restrict tool names, so an adapter installed
  as a plugin is selectable from Nextflow as well.

### Scaling

Chunked dereplication ([usage.md](usage.md#chunked-dereplication)) adds
levels over the union of chunk representatives only when a probe suggests
they pay.

* **Shrink rule.** A probe runs when the union U is larger than two chunks
  and the level removed at least a tenth of its input or a chunk's worth
  (N), whichever is larger. Otherwise the union goes to the final pass.
* **Probe.** One chunk of the next level, filled with the first N
  representatives in turn from each chunk of the previous level, so that
  genomes from different chunks meet.
* **Decision.** A level costs the per-genome work (sketching, N50) for every
  representative plus about N/U of the final pass's pairwise work, so it is
  expected to pay when it removes more than about N/U of its genomes. When
  the probe does, the rest of the level runs with the probe as its first
  chunk; otherwise the probe's result is dropped and the union goes to the
  final pass.
* **Thresholds.** The probe, the later levels and the final pass compare at
  the final thresholds, not the `--pre-*` ones.

The final pass logs its pairwise factor, the union squared over the level's
input times the chunk size. It can take longer than the factor suggests,
since the per-genome work is paid again and one call holds all threads.

### The phylo stage

The phylo stage combines three independent choices:

* genome set: dereplicated representatives or all genomes;
* MSA (multiple sequence alignment) source: a whole-genome aligner, a SNP
  typer's core-SNP alignment, or a concatenated core-gene alignment built
  from the `annotate` stage's annotations (`--msa-source aligner`, `snptype`
  or `coregene`, in `stages/phylo.py` and `stages/coregene.py`); skipped for
  alignment-free tree builders;
* tree builder: MSA-based (iqtree, fasttree, veryfasttree, raxmlng, aster) or
  alignment-free (mashtree, sourmash, skani, ska2).

Rooting is done by the stage after the builder returns (`tree/`). The rules
users see are in [phylogeny.md](phylogeny.md#rooting-the-tree), and those of
tree2tax in [usage.md](usage.md#rooting-in-tree2tax).

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
dataclass, and the owning stage writes the contract files. The tools within
a family are therefore interchangeable with no downstream change.

### Record names

A genome has one name in alignments, trees and tables:
`core.contracts.record_name(path)`, its file name without the FASTA suffix
and `.gz` (`x.fasta.gz` gives `x`; `Path.stem` would give `x.fasta`). An
aligner, SNP typer or tree builder adapter names its records and leaves with
it, and maps back to it where the tool rewrites names (ParSNP's `.ref`,
cactus's `.` to `_`). The phylo stage checks tree leaves against the same
names and rebuilds a reused alignment whose records differ from them.

### Genome sketches

The layout and parameters of `sketches/` are in
[output.md](output.md#genome-sketches). The code is in `core/sketches.py`.

* **Writers.** The genome-writing stages and `repgenr sketch` call
  `sketch_stage_genomes` after the manifest is updated and before the stage
  record is written. Sketches are not in `STAGE_DELIVERABLES`, and `sketch`
  is excluded from the resume fingerprint.
* **Records.** The sketch file, parameter string and FASTA SHA-256 are
  recorded in the manifest, or in `sketches/outgroup.json` for an outgroup.
  Staleness is decided from those records, not from file times.
* **One specification.** Every sourmash caller derives its
  `ToolCapabilities` from `core/sourmash.py` (`SOURMASH_TOOL`,
  `sourmash_capabilities(...)`), and a unit test checks that the Nextflow
  `SKETCH` module names the same image.
* **Consumers.** A consumer checks each sketch against its record before
  comparing, writes and records the missing and stale ones, and matches
  `derep/representatives/` to their genomes by file name. The SHA-256 is
  cached in `sketches/.digests.json` and recomputed only when the file's size
  or modification time changes.
* **Reuse.** `assemble --append` sketches only the genomes it adds, and
  `ingest --from-workdir` copies a source sketch whose digest matches.
  `assemble` sketches each accepted assembly once
  (`assemblies/<run>/contigs.sig.zip`); the classifier gathers with it, and
  the sketch step copies it under the genome's record name. An assembly
  excused after classification gets no sketch under `sketches/`. With
  `--no-sketch` the classifier sketches for itself.
* **Forced runs.** `repgenr --force sketch` counts the sketches that were
  current as `replaced (--force)`, apart from stale ones.

The adapters that read the sketches declare it through `sketch_request`
(see [Optional hooks](#optional-hooks)).

## 3. Tool families and the registry

There are ten tool families. Each is an abstract base class (ABC) plus a
`Registry` bound to an entry-point group. The table lists them in the order a
run uses them; `README.md` lists the adapters of each family.

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
applications) uses `registry.register("mytool", MyDereplicator)` and
`registry.unregister("mytool")`, or the `register_tool` pytest fixture.

### Adding an adapter

A new tool is a self-contained adapter module plus one entry-point line. The
steps below are for a built-in (in-tree) adapter. A third-party package
follows steps 1 to 3 in its own repository, adds the entry point to its own
package metadata, and depends on the distribution `repgenr2`.

1. **Write the adapter module** in the family's package
   (`src/repgenr/<family>/mytool.py`). Subclass the family ABC and declare
   `capabilities` (see [ToolCapabilities](#toolcapabilities)) with:
   `container`, `conda` or both; a `BinarySpec.min_version` floor for each
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
   floor appears in exactly one env file and that every binary and conda
   spec is covered. A package that does not solve with core gets a satellite
   `envs/<tool>.yml`, listed in `EXPECTED_FILES` in that test and in the
   `env` matrix of `.github/workflows/envs.yml`, which the test also checks.
   A tool distributed separately, as cactus is, is listed in
   `CONTAINER_ONLY` in the same test.

5. **Update the documentation tables**:
   * a row in the limits table of [choosing-tools.md](choosing-tools.md#declared-limits)
     (header `| Family | Tool | Declared limit |`), which
     `tests/unit/test_docs_tool_limits.py` checks against the registries;
     update the row also when an adapter's limit changes;
   * the adapter table in `README.md` and, if user-facing behaviour differs,
     a note in the user page that owns the stage ([usage.md](usage.md),
     [reads.md](reads.md) or [phylogeny.md](phylogeny.md));
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

7. **Add a live test** under `tests/live/` (see
   [Real-tool tests](#real-tool-tests)) so that
   [verification.md](verification.md) can report the adapter as verified
   against the real tool.

Everything else (CLI help, validation, `list-tools`, automatic selection,
containers, provenance, resume and Nextflow) picks the tool up from the
registry. The subsections below are the reference for steps 1 and 3.

### Class attributes the stages read

Set these per-family attributes, or inherit the default deliberately:

* `TreeBuilder.input_kind`: `InputKind.GENOMES` for alignment-free builders
  (they receive genome files), else the default `InputKind.MSA_FASTA` (they
  receive one MSA path). Narrow the argument with `as_genome_list` or
  `as_msa_path` from `treebuilders.base`.
* `TreeBuilder.supports_without_bootstrap`: False for a builder whose tree
  has no support values unless a bootstrap is requested (IQ-TREE); the stage
  warns only where that matters.
* `SnpTyper.requires_reference`: when True (the default) the stage resolves a
  reference genome before calling the adapter.
* `SnpTyper.produces_full_alignment`: whether `call` returns a whole-genome
  alignment (`SnpResult.full_alignment`), which a `Masker` reads and
  replaces with a masked variable-site FASTA. Set it whenever the tool can
  produce one. `--mask` is refused before SNP calling for a typer that
  declares none, and after calling for one that returns
  `full_alignment=None`.
* `Masker.tree_output()`: the masker's own tree from the last call, when it
  builds one (Gubbins does); the stage keeps it as `masker_tree.nwk`.
* `Annotator.requires_database`: True for a tool that needs a database
  directory (`AnnotateParams.db`, as for bakta).
* `CoreGeneParams.markers_db`: the database directory of the `markers`
  core-gene builder, from `--markers-db` or `REPGENR_MARKERS_DB`. The
  builder refuses a missing or incomplete directory (`check_markers_db`);
  `repgenr markers-db download` creates it. The other builders do not read
  the field.
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
  compares at (the default, None, requests no sketches). The stage then
  passes a mapping from genome path to `.sig.zip` in `DerepParams.sketches`
  (and the `sketches` keyword of `compare`), `TreeParams.sketches` or
  `ClassifyParams.sketches`. The adapter selects its k-mer size with `-k`
  and sketches any genome left out of the mapping itself; each signature is
  named by `record_name`. `core.sketches.adapter_sketches` builds the
  mapping with a provider: `workdir_provider` (writes and records missing
  and stale sketches) in the stages, and `directory_provider` (match by
  record name, no digest) in the steps given `--sketches-dir`. The adapters
  stay unaware of the working directory.

A tool without an override of `compare` or `distance_matrix` is refused
with a message naming the feature when a user selects it for that feature.

### ToolCapabilities

Every adapter declares a `ToolCapabilities` (`repgenr.core.plugins`):

```python
capabilities = ToolCapabilities(
    name="mytool",
    required_binaries=(BinarySpec("mytool", version_args=("--version",), min_version="1.2.0"),),
    default_params={"mode": "fast"},
    recommended_max_genomes=2000,
    container="quay.io/biocontainers/mytool:1.2.3--0",
    conda=("bioconda::mytool",),
    accepted_extras=frozenset({"mode"}),
    ignored_params=frozenset({"primary_ani"}),
    reads_gzip=True,
)
```

| Field | Type | Default | Read by |
|---|---|---|---|
| `name` | `str` | required | the registry, logs and `list-tools` |
| `required_binaries` | tuple of `BinarySpec` | empty | the native preflight (presence and minimum version) |
| `default_params` | `dict` | empty | the adapter; merged under `params.extra` by some callers (note 1) |
| `recommended_max_genomes` | `int` or None | None | `--tool auto` and the over-scale warnings (note 2) |
| `supports_native_scaling` | `bool` | False | the log wording of chunked dereplication only |
| `container` | `str` or None | None | the container backend (note 3) |
| `conda` | tuple of `str` | empty | Wave, when no `container` is set, and the env-file test (note 3) |
| `accepted_extras` | frozenset of `str` | empty | the unknown-key warning (note 4) |
| `ignored_params` | frozenset of `str` | empty | the ignored-parameter warning (note 5) |
| `reads_gzip` | `bool` | False | `snptype`, `phylo` and `annotate` (note 6) |
| `auto_eligible` | `bool` | True | `--tool auto` (note 2) |

1. **Defaults.** The dereplication stage and steps, `glance`, `annotate` and
   the tree builder in `phylo` merge `default_params` under `params.extra`;
   the other callers do not. Read a default with
   `params.extra.get(key, self.capabilities.default_params[key])`.
2. **Automatic selection.** The tightest fitting `recommended_max_genomes`
   wins, and unbounded tools are picked when nothing bounded fits
   ([choosing-tools.md](choosing-tools.md#how-auto-chooses)).
   `auto_eligible=False` marks a tool that suits one kind of input only (the
   ska2 tree builder, clonal sets); it must be named.
3. **Containers.** Set `container`, `conda` or both, or the tool cannot run
   under `--container` or `--wave`; the run then falls back to the host with
   a warning.
4. **Extras.** Users set extras with `--tool-arg key=value` or, on phylo,
   `--aligner-arg key=value`. A stage warns about keys that no adapter it
   runs declares; phylo warns only about keys neither of its two adapters
   reads. The CLI does not inject stage-level flags (such as `--virus`) into
   the extras of a tool that does not read them. Parse integer extras with
   `parse_extra_int`, which raises a `UserInputError` on a bad value.
5. **Ignored parameters.** List every standard parameter the adapter never
   passes to its tool, for instance `primary_ani` on a tool with no
   pre-clustering step. The stage warns when a user sets one to a
   non-default value, so the request is not dropped without notice.
6. **Gzipped input.** With `reads_gzip=False`, the stages give an aligner,
   SNP typer, genome-input tree builder or annotator decompressed copies
   (`core.process.stage_plain_inputs`), named by `record_name`. Set it to
   True only after checking the tool on a `.fasta.gz` genome: a tool that
   misreads gzip may exit 0 with an empty result (SibeliaZ does).

**Version queries.** Each query runs in its own process group with a limit
of `core.binaries.VERSION_TIMEOUT` (30 s; 8 s under `list-tools --check`).
A binary that does not answer, or prints no version number, is looked up in
its environment's `conda-meta` records; without one its version is
`unknown`, and a `strict_version` floor fails.

### Running the external tool

**One entry point.** Use `run_tool(self.capabilities, cmd, logger=logger,
...)` from `repgenr.core.containers` for every external call; never call
`process.run` directly and never use `shell=True`. `run_tool` routes the
command through the active container backend (docker, singularity or Wave)
using the adapter's `container` and `conda` declaration, and runs it on the
host otherwise.

* Use `write_fofn` (`repgenr.core.process`) instead of a shell glob for large
  input lists. Declare the directories that hold the listed files through
  `extra_mounts=`, since the contents of a file of file names are invisible to
  the container backend.
* Pass a result written to stdout through `stdout_path=`; the file is
  published atomically, and only on success.
* Variables the tool itself must see (such as `OMP_NUM_THREADS`) go in
  `tool_env=`, which reaches the tool in every backend. `env=` is the
  environment of the process repgenr starts, which under Docker is the
  `docker run` client, not the tool.
* `run_chain(self.capabilities, [(prefix, argv), ...], logger=logger)` runs
  several commands in a row, such as a per-genome mapping and calling chain:
  on the host as a loop over `run_tool`, under a container backend in one
  container. Steps are plain argument vectors, so a tool that writes to
  stdout needs an output flag, and paths must be absolute to be mounted.

The return value of `preflight()` (host versions or the container image
reference) is recorded in the `repgenr.yaml` provenance and appears in
`repgenr versions` and the Nextflow `versions.yml` bridge.

**Console lines.** A tool command line is shortened on the console to 120
columns; the full text, and under `--container` the engine options, are in
`repgenr.log`. A run of input paths shrinks first (first path and a count,
then the count alone), a path that follows an option counts as that
option's value, and long genome lists show the first five and a count.

**Stopping containers.** Under `--container docker` the container runs with
`--init`, is named `repgenr-<repgenr pid>-<hex>`, and carries the labels
`repgenr.pid` and `repgenr.host`. When repgenr is stopped it forwards
SIGTERM and runs `docker stop -t 0` on the containers of its running tools,
parallel workers included. A second signal ends repgenr at once and starts
one `docker stop -t 0` for all of them, which completes after repgenr has
exited.

### Resume semantics

A stage's skip decision digests its declared inputs (`STAGE_INPUTS` in
`repgenr/cli/base.py`), its parameters and the container identity. An adapter
therefore needs no resume logic of its own, but it must read only the inputs
the stage declares: a tool that reads undeclared files would resume
incorrectly. [usage.md](usage.md#resume-and---force) summarises the rules
for users.

**Inputs.** A declared directory (`genomes/`, `derep/representatives/`,
`outgroup/`, the `ingest --genomes-dir` source) is digested from its genome
FASTA files only (`core.inputs.dir_stat_digest`, the files `list_fasta`
returns: names ending in `.fasta`, `.fa`, `.fna`, `.fas`, `.fasta.gz`,
`.fna.gz` or `.fa.gz`). A leftover `x.fasta.tmp`, a `.fai` index, a
`truth.json` or a README beside the genomes therefore does not rerun the
stages that read the directory. The digest uses file metadata (name, size,
mtime), so an in-place edit that preserves size and mtime is not detected.
A stage that reads anything else from a directory (a table, an index)
declares that file by its own path; derep-stock's unpack, for example, lists
the stored run's `clusters.tsv`, `genome_status.tsv`, `record.json` and
`representatives/` one by one.

**Outputs.** Outputs are not part of the fingerprint. A matching stage is
skipped only when every path in its entry in `STAGE_DELIVERABLES` (same
module) exists, for example `derep/clusters.tsv` and
`derep/representatives/` for `dereplicate`.

* A directory must not be empty, and dotfiles (`.DS_Store`, exFAT `._`
  files) do not count. `genomes/` and `derep/representatives/` count only
  with a genome FASTA file in them.
* Each genome in `selection.tsv` is checked for the stage that wrote the
  genome set (`genome`, `ingest`, `vgenome`, `assemble`), and each
  representative in `derep/clusters.tsv` for `dereplicate`.
* A missing output is logged as
  `Stage 'X': deliverable <path> missing; re-running.` and the stage runs
  again without `--force`, so a genome deleted from `genomes/` is restored.
  To drop a genome, remove it from the source or the `--selection` file and
  rerun the entry stage.
* For `glance` only `glance_clustering_dendrogram.pdf` is checked, since its
  plots are absent when no value falls within the plot bounds.

**Records.** A stage writes its record without a completion stamp before it
starts, and a successful run stamps it. A stage that failed or crashed
mid-run is therefore listed as `[interrupted]` by `status`, reported as a
failure by `doctor`, and always reruns. A failure in parameter validation
writes no record. A stage that refuses (exit 2 or 3, or 4 for a tool missing
at its preflight) without changing any of its main outputs leaves the record
as it was: none on a first run, and the last finished one on a rerun.
Examples are `phylo` with fewer than three genomes and `snptype --mask
gubbins` without Gubbins. A failed external tool (exit 6), or a refusal
after a main output changed, leaves the record interrupted.

**A new stage** adds an entry to both `STAGE_INPUTS` and
`STAGE_DELIVERABLES`. A stage whose one record serves several invocations
(derep-stock's named runs) can also register a check in
`_STAGE_PRECHECKS`. The check runs before the harness marks the record
incomplete, so a refused call leaves the last finished one intact.

**reads-sketch.** The stage records the variant of its last call and the
runs each variant holds.

* A transient failure on the archive side (a dropped connection, a stalled
  transfer, a pipeline past its time limit) is tried again within the call:
  up to three attempts, with pauses of about 30 s and 120 s plus a quarter
  of jitter. An accession NCBI will not serve, and a sourmash failure, are
  not retried.
* With `--max-spots` each run's pipeline has a time limit of 300 s plus one
  second per 200 spots (about 22 min at 200000), or
  `REPGENR_SUBPROCESS_TIMEOUT` when that is lower. A connection that stalls
  then frees its job slot and is tried again. Without a cap only the global
  limit applies, and a run killed by it is tried again the same way (up to
  three times the limit for one run).
* One call writes a variant at a time. A second call on the same variant
  directory is refused while the first holds its `.lock` file, since each
  call clears the partial files of the directory and rewrites its index. On
  a file system without locks the call goes on with a warning, and on a
  network file system the lock may hold on one host only.
* Each sourmash process runs single-threaded (`RAYON_NUM_THREADS=1`); on
  macOS the default thread pool took twice the time for the same signature.
* Under a container back-end both images are pulled once before the first
  run, and under Docker an engine that cannot start a tool stops the call
  with one error instead of failing every run.

**Reads screen.** Each decision records both thresholds and the resolved
paths of the GTDB sketch and lineages CSV.

* `assemble --screen-reads`: a repeat with the same settings excuses runs
  from their recorded decisions without a download. A changed setting
  screens again from the kept reads sketches; an excused run is fetched
  again only when it now passes, and a finished run that now fails is
  excused while its contigs and marker stay. A finished run without a reads
  sketch is not screened. With `--max-runs`, a run excused during the call
  frees its place only in the next call.
* `reads-screen`: a new sketch or a changed threshold screens a run again.
  A run the tools failed on keeps no decision. The duplicate check is
  repeated on every call, in `reads.tsv` order, over every run decided with
  the same settings, so several calls give the result of one.

**Refinement.** `derep-refine` needs a completed `dereplicate` record (one
restored by `derep-stock` counts) and keeps its seeds and ANI in its own
record. A seed refined before brings in its original coarse cluster from
`refined.tsv`. A cluster claimed by two ANIs of a `--refine-file` is
refused. The file's content is a resume input, and its path is stored as an
absolute path. On a result reduced with `--reduce`, a lineage is refined as
the whole taxon group, with a warning.

## 4. Testing

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
`test_cli_matrix.py` and `test_version_consistency.py` (CHANGELOG.md) under
`tests/unit/`.

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

That means an entry in the family's token dictionary, a canned-output branch
in the fake runner, and the name in the suite's literal set of built-in
adapters; the suites do not discover entry points. The other families have
their own modules (for example `test_annotators.py`, `test_coregenes.py`).
Integration-test the owning stage with the `register_tool` fixture
(`tests/conftest.py`), which registers an adapter class for one test.

### nf-test stub tests

The Nextflow tests live in `nextflow/tests/` and run under the `test`
profile (`nf-test.config`). Tests tagged `stub` use stub blocks, so they need
no external tools or network; they check channel wiring and the file
contracts. Run them before pushing a Nextflow change:

```bash
nf-test test --tag stub
```

Nextflow resolves `bin/` and the default config relative to the launched
script, so the harness scripts under `nextflow/tests/` reach the versions
helper through the symlink `nextflow/tests/bin/repgenr_versions_fragment`,
and a script added to `nextflow/bin/` needs its own link there. Run directly,
a harness needs `-c nextflow/nextflow.config`.

Two harnesses start from a directory of genome FASTAs instead of the GTDB or
NCBI front end. `nextflow/tests/local_dataflow.nf` feeds it to
`DEREPLICATE_SCATTER`, as `repgenr ingest` does on the command line, and
`nextflow/tests/dereplicate_scatter.nf` runs the scatter-gather
dereplication alone, which also dereplicates a local collection:

```bash
nextflow run nextflow/tests/dereplicate_scatter.nf -c nextflow/nextflow.config \
    --genomes_dir <DIR> --derep_tool sourmash --derep_process_size 2000 \
    --outdir results -profile standard
```

Add `-stub` to exercise the wiring without running the tools.

### Real-tool tests

`tests/live/` holds pytest tests marked `live` that run the installed
`repgenr` console script against real tools on seeded synthetic genome
sets. They are deselected by default and never run in CI.
[verification.md](verification.md#the-live-suite) records their results.

Copy `tests/live/live.example.toml` to `live.local.toml`, point it at the
conda environments that hold each tool, and run:

```bash
conda run -n repgenr_dev --no-capture-output \
    pytest tests/live -m live --live-config tests/live/live.local.toml -ra
```

The markers `network` and `container` select the tests that need GTDB or
NCBI, or Docker; `-m "live and not network and not container"` is the
offline subset. Tests whose tool is not on `PATH` are skipped, not failed.
`tests/live/README.md` describes the configuration keys.

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

## 5. Conventions

* External tools run only through `run_tool` or `run_chain`, with argument
  vectors and files of file names; no `shell=True` and no shell globs.
* Argument lists stay below `ARG_MAX`. A tool without a list-file option
  gets a staged directory where it reads one (skDER, ParSNP). The few that
  take every path on argv (SibeliaZ, snippy-core) call `warn_argv_bytes`.
  `process.run` measures the final argv (container preamble included) plus
  the environment against `os.sysconf("SC_ARG_MAX")` and, on Linux, the
  128 KiB single-argument limit, and raises `UserInputError` before the
  exec instead of failing with E2BIG.
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
  tests named in [Adding an adapter](#adding-an-adapter) check this.
