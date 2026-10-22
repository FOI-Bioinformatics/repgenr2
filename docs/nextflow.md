# Nextflow

This page describes the Nextflow pipeline, which runs the same stages as the
command line on a cluster or in the cloud. Use it when the work should be
spread over many tasks, for example one assembly per sequencing run or one
dereplication chunk per node.

## Contents

1. [Running the pipeline](#running-the-pipeline)
2. [Pipeline structure](#pipeline-structure)
3. [Parameters](#parameters)
4. [Profiles and containers](#profiles-and-containers)
5. [Reads mode](#reads-mode)
6. [Genome sketches](#genome-sketches), [split alignment and tree tasks](#split-alignment-and-tree-tasks)
   and the [core-gene route](#core-gene-route)
7. [Scaling](#scaling) and [configuring processes](#configuring-processes)

## Running the pipeline

The Nextflow layer runs the stages as typed data channels. Each stage emits
its outputs (the metadata selection, genome FASTAs, per-chunk and merged
representatives, the tree, the taxonomy) as staged files that the next stage
consumes. There is no shared working directory, and results are published
under `--outdir` (see [output.md](output.md#nextflow-output-directory)).
Nextflow owns the fan-out.

```bash
nextflow run nextflow/main.nf \
    --mode bacterial \
    --metadata_args '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella' \
    --outdir results \
    -profile standard
```

The pipeline requires Nextflow 26.04 or later (see
[install.md](install.md#nextflow)), and `repgenr` must be available to every
process. Add `-stub` to any run for a quick wiring check without external
tools.

## Pipeline structure

`nextflow/main.nf` dispatches by `--mode` to one of three subworkflows that
share the dereplication, phylo and tree2tax modules:

```
bacterial: ACQUIRE       (metadata -> genome)          -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
viral:     VACQUIRE      (vmetadata -> vgenome)        -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
reads:     ACQUIRE_READS (reads -> assemble per run
                          -> [genome-qc] -> gather)    -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
```

`metadata` emits a portable `selection.tsv`, and `genome-fetch` downloads
the genomes and emits them as a channel feeding the scatter-gather
dereplication. `phylo` and `tree2tax` run in task-local working directories
and emit `tree.nwk`, `tree2tax.tsv` and `genomes_map.tsv` to `--outdir`.
`PHYLO` (and the split `PHYLO_MSA` and `PHYLO_TREE`) publish the alignment
they built (`phylo/align/` for an aligner, `phylo/tree/msa/` for a SNP
typer, with the reuse stamp) and the tree builder's own files under
`phylo/tree/`, beside the tree.

Every channel carries a meta map built once per run. `id` is a slug of the
selection target (`francisella` for `-tg francisella`), and `mode` is
`bacterial` or `viral`. Task tags and the dereplication chunk names use it.

## Parameters

`nextflow run nextflow/main.nf --help` prints a summary, and
`nextflow/nextflow_schema.json` holds every parameter with its default and
description. Parameters are validated against the schema at launch. The
most used are:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Published results and execution reports. |
| `--mode` | `bacterial` | `bacterial` (GTDB), `viral` (NCBI Virus; BV-BRC with `--vmetadata_args "--source bvbrc"`) or `reads` (sequencing runs, assembled per run). |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--vmetadata_args`, `--vgenome_args` | see config | Arguments for the viral metadata and genome selection stages. |
| `--reads_args` | see config | Arguments for the reads stage. |
| `--derep_tool` | `skder` | Dereplicator for the scatter-gather step. |
| `--derep_process_size` | unset | Genomes per dereplication chunk; a single chunk when unset. |
| `--derep_secondary_ani` | `0.99` | Secondary ANI threshold (also `--derep_primary_ani`, 0.90, and `--derep_aligned_fraction`, 0.50). |
| `--derep_keeper` | `quality` | Representative choice in the chunk and merge steps: `quality`, `gtdb` or `tool`, read from `selection.tsv`. |
| `--phylo_args` | `--treebuilder mashtree` | Aligner or tree builder for the phylogeny, its settings and the rooting. |
| `--tree2tax_args` | (empty) | tree2tax arguments; dereplicated genomes are listed by default (`--no-include-dereplicated` omits them). |

Three parameters forward settings between stages:

- With `--group-segments` in `--vgenome_args`, the acquisition task passes
  `segments.tsv` to `tree2tax-relations`, so `genomes_map.tsv` lists the
  member segment accessions of each grouped isolate. `segments.tsv` itself
  is not published.
- `--tree2tax_args` gains `--no-outgroup` when `--phylo_args` contains it,
  and `--keep-root` when it contains `--root midpoint`.
- `--derep_reduce`, `--derep_target_reps` and the refinement parameters act
  on the merge step (see [Scaling](#scaling)).

### Passing parameters

Pass parameters on the command line (`--name value`) or, for a reusable and
typed configuration, in a params file:

```bash
nextflow run nextflow/main.nf -profile standard -params-file params.json
```

A starting point is in `nextflow/assets/params_example.json`. On the command
line, give a string parameter whose value starts with a dash as
`--key=value` (for example `--tree2tax_args=--include-dereplicated`), since
Nextflow reads `--key --flag` as the boolean `key = true`. A params file
keeps numeric values typed (for example `"derep_process_size": 2000`). On
the command line they arrive as strings, which the schema also accepts for
the numeric options.

## Profiles and containers

Combine an executor profile with an optional container profile, for example
`-profile slurm,singularity`.

- **Executors**: `standard` (local) and `slurm`.
- **Containers**: `docker`, `singularity` and `wave`. They set
  `params.repgenr_opts`, so every stage calls `repgenr --container ...` and
  RepGenR2 runs each tool in its own container (see
  [usage.md](usage.md#running-tools-in-containers)). Under `docker` or
  `singularity` every adapter with a pinned image runs in it. The `simple`
  SNP typer, parsnp, skder and SibeliaZ have no pinned image and run on the
  host, and the log says so. `wave` mints an image from each adapter's conda
  specification instead, which also covers those four. The `wave` profile
  selects Singularity as the container runtime (`--container singularity
  --wave`).
- **`test`**: minimal resources and a small target for a quick smoke run.

The container profiles set no container for the Nextflow processes:
`repgenr` runs on the task host and starts each tool's container itself.
Running `repgenr` inside a Docker process container as well would need one
container engine inside another.

The `standard` profile (the default, local executor) caps every process's
CPU and memory request to what the machine has. `base.config` asks for up to
32 CPUs and 128 GB for the heavy processes, which the local executor would
otherwise refuse. `slurm` leaves the requests as they are and sets only the
executor. The queue or account is site-specific, so pass it in a site config
(`-c site.config` with `process.queue = '...'`).

The same holds for a cloud executor such as AWS Batch. The pipeline ships no
cloud profile, because the region, the job queue and an image that provides
`repgenr` and the tools are values only the site knows. A site config that
sets them is all that is needed.

To give one tool a satellite environment, pass `--bin-dir` through
`--repgenr_opts`, for example `--repgenr_opts '--bin-dir gubbins=/p/bin'`.
This value replaces the one a container profile sets, so with `-profile
docker` give `--repgenr_opts '--container docker --bin-dir gubbins=/p/bin'`.

## Reads mode

`--mode reads` builds the genome set from sequencing runs, as
[reads.md](reads.md) describes for the command line, with one task per run.

```bash
nextflow run nextflow/main.nf --mode reads --outdir results \
    --reads_args "-tg mycoplasmopsis --platform illumina --max-runs 20" \
    --assembler auto --checkm2_db /db/checkm2 \
    --gtdb_sketch /db/gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb_lineages /db/lineages.csv \
    -profile slurm,singularity
```

- `READS_SELECT` runs `repgenr reads` and emits `reads.tsv`. The number of
  runs is bounded on the selection with `--reads_args '--max-runs N'`.
- Every row becomes one `READS_ASSEMBLE` task (`assemble-run`, label
  `process_assembly`: 8 CPUs, 16 GB for a short-read run and 32 GB for a
  long-read one, 6 h, scaled by attempt), so the runs assemble in parallel
  on a cluster. `--assembler` and `--assemble_args` reach it.
- When `--checkm2_db` or `--gtdb_sketch` is set, one `GENOME_QC` task
  (`genome-qc`, `process_medium`) scores and classifies the batch.
- `READS_GATHER` (`reads-gather`) applies the gate and the naming policy and
  emits the genome FASTAs, `selection.tsv` and an empty outgroup accession
  file, the same tuple `ACQUIRE` emits.

`reads.tsv`, `selection.tsv`, `assembly_stats.tsv`, `excused_runs.tsv` and
the QC tables are published under `reads/`. The reads sketch of each
assembly task is turned off with `--assemble_args '--no-reads-sketch'`.
`--assemble_args '--screen-reads'` turns on the in-assembly screen, and
`conf/modules.config` then adds `--gtdb-sketch` and `--gtdb-lineages` from
the pipeline parameters of the same names. That screen applies the taxon and
fraction checks only, since each task sees one run.

### Reads sketch and screen

`--reads_sketch` adds a sketch task between the selection and the
assemblies, and `--reads_screen` a screen task after it; the screen implies
the sketch. Both apply to `--mode reads` only.

- `READS_SKETCH` (`reads-sketch --reads-tsv reads.tsv --accession-file
  chunk.txt -o sketches --allow-all-failed`, `process_low`) streams the runs
  from the archive into sourmash sketches, the first `--reads_sketch_spots`
  spots of each when set. `reads.tsv` is cut into chunks of
  `--reads_sketch_chunk` runs (default 25), one task per chunk. A chunk that
  hits its time limit or is retried streams its own runs again, and the
  others keep their work. On a scheduler a chunk is one queued job, so a
  smaller chunk gives more parallel streams at the cost of more jobs. Every
  chunk task takes the whole `reads.tsv`, so a change to it runs every chunk
  again.
- `READS_SKETCH_MERGE` joins the chunks into one `sketches/` with one index
  per variant in `reads.tsv` order (a run no chunk reached is `pending`). It
  publishes the result under `reads/reads_sketches/<variant>/`. A chunk where
  every run failed is a warning, and the merge fails the batch with exit 3
  when no run of `reads.tsv` streamed at all.
- `READS_SCREEN` (`reads-screen --reads-tsv reads.tsv --sketches-dir
  sketches -o screen`, `process_medium`) judges the sketches with the three
  checks against `--gtdb_sketch` and `--gtdb_lineages`, with
  `--reads_screen_args` for `--min-fraction` and `--dup-ani`. It publishes
  `reads/reads_screen/<variant>/` with `reads_screen.tsv`, one decision per
  run and `passed.txt`.

Only the runs of `passed.txt` become `READS_ASSEMBLE` tasks. The others are
neither fetched nor listed in `excused_runs.tsv`, since they never reach the
assembly step; their reasons are in `reads_screen.tsv`. That includes a run
whose stream failed in the sketch task (`screen_failed: no reads sketch`),
which without the screen would still have been assembled. A transient
archive error therefore drops the run from the batch, so check
`reads_screen.tsv` for such rows. A screen that passes no run stops the
pipeline, and `--reads_screen` without `--gtdb_sketch` and `--gtdb_lineages`
is refused before any task. `--reads_sketch` alone publishes the sketches
without a screen.

Both tasks need `fastq-dump` (sra-tools) and sourmash in the task
environment. The screen also needs the GTDB sketch on a path every task can
read, given as parameters, since the pipeline does not read the
`REPGENR_GTDB_*` environment variables.

```bash
nextflow run nextflow/main.nf --mode reads --outdir results \
    --reads_args "-tg mycoplasmopsis --platform illumina --max-runs 50" \
    --reads_screen --reads_sketch_spots 200000 \
    --gtdb_sketch /db/gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb_lineages /db/lineages.csv
```

## Genome sketches

`--sketch` adds a `SKETCH` task on the bacterial path. It writes
`sketches/<record name>.sig.zip` for every genome and the outgroup (k=21, 31
and 51 at scaled=1000), publishes the directory under `--outdir`, and stages
it into every `DEREP_CHUNK`, `DEREP_MERGE` and `PHYLO` (or `PHYLO_TREE`)
task, which pass it to repgenr as `--sketches-dir`. A sourmash dereplicator
or tree builder then reads the sketches instead of sketching, and other
tools ignore them.

`SKETCH` calls sourmash itself, so the parameter is off by default and a run
without sourmash is unaffected. The viral and reads paths do not stage
sketches. In the stateless steps a sketch is matched to a genome by record
name only, without the digest check of a working directory. This holds
because `SKETCH` sketches the same genome files the steps receive.

## Split alignment and tree tasks

`--phylo_split_msa` splits the phylogeny into `PHYLO_MSA` and `PHYLO_TREE`
(see [Reusing an alignment](phylogeny.md#reusing-an-alignment-across-tree-builders)).
The alignment then keeps its own cache entry, so trying another tree builder
or another bootstrap reruns only the tree. The two halves take their own
resource labels: the alignment `process_high`, the tree `process_medium`.
The option applies to tree builders that consume an alignment; leave it off
for mashtree or sourmash, which build from genomes.

## Core-gene route

With `--msa-source coregene` in `--phylo_args`, or `--preset genus` or
`family`, the bacterial and reads paths annotate before the phylogeny.
`--preset` is expanded by `phylo-build`. `ANNOTATE` runs `repgenr
annotate-build` once per representative, and once for the outgroup unless
`--phylo_args` holds `--no-outgroup`. The annotations therefore run in
parallel, and Nextflow caches them genome by genome. A rerun with another
tree builder or core-gene builder annotates nothing again. The record
directories are published under `annotations/` and gathered into the `PHYLO`
task (or `PHYLO_MSA` with `--phylo_split_msa`). That task indexes them with
`repgenr annotate-index` and passes them to `phylo-build` as
`--annotations-dir`.

- `--annotator` (`bakta` or `prokka`) selects the annotator, and
  `--annotate_args` passes extra `annotate-build` arguments.
- `--bakta_db` is linked into each `ANNOTATE` task; without it bakta reads
  `BAKTA_DB` from the task environment. Under `stageInMode 'copy'`, or on a
  cloud executor without a shared file system, the database would be copied
  once per genome, so there give the tasks a path they can read through
  `BAKTA_DB` instead.
- `--markers_db` is staged into the `PHYLO` (or `PHYLO_MSA`) task as
  `--markers-db` for `--coregene markers`. Without it `phylo-build` reads
  `REPGENR_MARKERS_DB` from the task environment.
- The route is detected from `--phylo_args`: quotes around a value are
  dropped and the last `--msa-source` wins, as on the command line. A
  `withName` override of the phylo task's `ext.args` is not seen.
- Each record directory carries `annotation.json` with the annotator, its
  version and the database version. `annotate-index` indexes each record
  with the tool named there (its `--tool` applies only to records without
  the file), so a set annotated by two tools is refused by the core-gene
  pass rather than read under one name.
- The viral path has no annotation step, since the annotators are
  prokaryotic. With `--msa-source coregene` in `--phylo_args` it logs a
  warning, as the phylo task would then annotate every genome itself.
- The annotator, the core-gene builder and the tree builder must each be
  reachable from their tasks: through `--repgenr_opts '--bin-dir
  bakta=<env>/bin --bin-dir ppanggolin=<env>/bin'` for satellite
  environments, or through a container profile.

## Scaling

For large dereplication inputs (10000 genomes or more), set a chunk size so
the dereplication scatters across tasks, one per chunk:

```bash
nextflow run nextflow/main.nf --outdir results \
    --derep_tool sourmash --derep_process_size 2000 -profile slurm -c site.config
```

The `DEREPLICATE_SCATTER` subworkflow groups the genomes into chunks of
`--derep_process_size`, dereplicates each chunk as a separate task, and
dereplicates the union of the chunk representatives once more in one merge
task. It wraps the `repgenr dereplicate-chunk` and `dereplicate-merge`
steps and is driven by `--derep_tool` and the `--derep_*_ani` and
`--derep_aligned_fraction` thresholds. Compared with
[chunked dereplication](usage.md#chunked-dereplication) on the command line,
it runs one chunk round and one merge task over the whole union. The keeper
is applied within each chunk and again after the merge. `--derep_reduce`
and `--derep_target_reps` act on the merge step, where the final set is
decided.

With `--derep_refine SEEDS --derep_refine_ani F` (comma-separated seeds) or
`--derep_refine_file refine.tsv`, a `DEREP_REFINE` task follows the merge.
The merged result and every genome of the run are staged in, and the
clusters of the seeds are dereplicated again at the finer ANI (`repgenr
dereplicate-refine`). The refined result, published as
`dereplicate/refined/` with its `refined.tsv`, is the one the tree is built
from, and `dereplicate/merged/` stays published beside it as the coarse
result (see [Finer resolution inside a lineage](usage.md#finer-resolution-inside-a-lineage)).
The refined directory holds the contract files and no `genome_n50.tsv`.
Seeds without an ANI, an ANI without seeds, or an ANI not above
`--derep_secondary_ani` are refused at launch. The reads route runs the same
subworkflow, so the refinement applies there too.

Resource labels (`process_low`, `process_medium`, `process_high`, and
`process_assembly` for the per-run assembly tasks of reads mode, whose
memory follows the read platform) scale memory and time with the retry
attempt. A task killed for memory or time is therefore resubmitted with more
headroom. Tune the label values per environment in
`nextflow/conf/base.config`.

## Configuring processes

Each process reads its tool flags from `task.ext.args`, which
`nextflow/conf/modules.config` maps from the user-facing parameters:
`--metadata_args`, `--vmetadata_args` and `--vgenome_args`, `--phylo_args`,
`--tree2tax_args`, and five dereplication parameters (`--derep_tool`,
`--derep_primary_ani`, `--derep_secondary_ani`, `--derep_aligned_fraction`,
`--derep_keeper`). The last five are composed into one string for the two
dereplication processes. `--derep_reduce` and `--derep_target_reps` are
passed to the merge process only, as `ext.args2`. `--derep_process_size` is
read by the scatter subworkflow to size the chunks and is not a tool flag.
Publishing directories live in the same file, and resources and the retry
window in `nextflow/conf/base.config`.

A site can retune one process without changing the pipeline by passing its
own config:

```groovy
// site.config
process {
    withName: 'PHYLO' {
        ext.args = '--treebuilder iqtree -B 1000'
    }
    withName: 'DEREP_CHUNK|DEREP_MERGE' {
        ext.args = '--tool galah --secondary-ani 0.98'
    }
}
```

Run it with `nextflow run nextflow/main.nf -c site.config ...`.
`ext.repgenr_opts` is the hook for the top-level repgenr options every
process prepends; the container profiles set it through `--repgenr_opts`.

An `ext.args` override replaces the whole string for that process. For the
dereplication processes, every flag not repeated in the override then falls
back to the repgenr default rather than to the `--derep_*` parameter.

The test harnesses under `nextflow/tests/` are described in
[developing.md](developing.md#nf-test-stub-tests).
