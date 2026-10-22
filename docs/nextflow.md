# Nextflow

This page describes the Nextflow pipeline, which runs the same stages as the
command line on a cluster or in the cloud. Use it when the work should be
spread over many tasks, such as one assembly per sequencing run or one
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

The stages run as typed data channels: each emits its outputs as staged
files that the next consumes. There is no shared working directory, and
results are published under `--outdir` (see
[output.md](output.md#nextflow-output-directory)).

```bash
nextflow run nextflow/main.nf --mode bacterial \
    --metadata_args '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella' \
    --outdir results -profile standard
```

The pipeline needs Nextflow 26.04 or later ([install.md](install.md#nextflow)),
and `repgenr` must be available to every process. `-stub` checks the wiring
without external tools.

## Pipeline structure

`nextflow/main.nf` dispatches by `--mode` to one of three subworkflows that
share the dereplication, phylo and tree2tax modules:

```
bacterial: ACQUIRE       (metadata -> genome)          -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
viral:     VACQUIRE      (vmetadata -> vgenome)        -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
reads:     ACQUIRE_READS (reads -> assemble per run
                          -> [genome-qc] -> gather)    -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
```

`PHYLO` publishes the alignment it built and the tree builder's files under
`phylo/`. Every channel carries a meta map whose `id` is a slug of the
selection target, used in task tags and chunk names.

## Parameters

`nextflow run nextflow/main.nf --help` prints a summary, and
`nextflow/nextflow_schema.json` holds every parameter with its default and
description; parameters are validated against it at launch. The most used:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Published results and execution reports. |
| `--mode` | `bacterial` | `bacterial` (GTDB), `viral` (NCBI Virus; BV-BRC with `--vmetadata_args "--source bvbrc"`) or `reads`. |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--vmetadata_args`, `--vgenome_args` | see config | Arguments for the viral stages. |
| `--reads_args` | see config | Arguments for the reads stage. |
| `--derep_tool` | `skder` | Dereplicator for the scatter-gather step. |
| `--derep_process_size` | unset | Genomes per dereplication chunk; one chunk when unset. |
| `--derep_secondary_ani` | `0.99` | Secondary ANI (also `--derep_primary_ani`, 0.90, and `--derep_aligned_fraction`, 0.50). |
| `--derep_keeper` | `quality` | Representative choice: `quality`, `gtdb` or `tool`, from `selection.tsv`. |
| `--phylo_args` | `--treebuilder mashtree` | Phylogeny route, settings and rooting. |
| `--tree2tax_args` | (empty) | tree2tax arguments; dereplicated genomes are listed by default. |

`--tree2tax_args` gains `--no-outgroup` when `--phylo_args` contains it, and
`--keep-root` with `--root midpoint`. With `--group-segments` in
`--vgenome_args`, `genomes_map.tsv` lists the segment accessions of each
isolate.

### Passing parameters

A params file (`-params-file params.json`; a starting point is
`nextflow/assets/params_example.json`) keeps values typed. On the command
line, give a string value that starts with a dash as `--key=value` (for
example `--tree2tax_args=--include-dereplicated`), since Nextflow reads
`--key --flag` as a boolean.

## Profiles and containers

Combine an executor profile (`standard`, `slurm`) with an optional
container profile (`docker`, `singularity`, `wave`), for example `-profile
slurm,singularity`. `test` runs a small smoke test.

- The container profiles set `--repgenr_opts`, so every stage runs
  `repgenr --container ...` and each tool runs in its own image while
  `repgenr` runs on the task host. The `simple` typer, parsnp, skder and
  SibeliaZ have no pinned image and run on the host, except under `wave`,
  which builds images from the conda specifications (`--container
  singularity --wave`).
- `standard` caps every process's CPU and memory request to the machine.
  `slurm` sets only the executor; give the queue in a site config (`-c
  site.config` with `process.queue = '...'`). No cloud profile ships, since
  the region, queue and image are site-specific.
- `--repgenr_opts '--bin-dir gubbins=/p/bin'` gives one tool a satellite
  environment. It replaces the value a container profile sets, so with
  `-profile docker` give `--repgenr_opts '--container docker --bin-dir
  gubbins=/p/bin'`.

## Reads mode

`--mode reads` builds the genome set from sequencing runs ([reads.md](reads.md)),
with one task per run:

```bash
nextflow run nextflow/main.nf --mode reads --outdir results \
    --reads_args "-tg mycoplasmopsis --platform illumina --max-runs 20" \
    --assembler auto --checkm2_db /db/checkm2 \
    --gtdb_sketch /db/gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb_lineages /db/lineages.csv \
    -profile slurm,singularity
```

`READS_SELECT` runs `repgenr reads`; bound the number of runs with
`--reads_args '--max-runs N'`. Each run becomes one `READS_ASSEMBLE` task
(`assemble-run`, label `process_assembly`: 8 CPUs, 16 GB for a short-read
run and 32 GB for a long-read one, 6 h, scaled by attempt), with
`--assembler` and `--assemble_args`. `--assemble_args '--no-reads-sketch'`
turns off the per-run reads sketch. With `--checkm2_db` or `--gtdb_sketch`, one `GENOME_QC`
task scores and classifies the batch, and `READS_GATHER` writes the genome
set. The tables are published under `reads/`. `--assemble_args
'--screen-reads'` turns on the per-run screen, which then receives
`--gtdb_sketch` and `--gtdb_lineages`.

### Reads sketch and screen

`--reads_sketch` streams the runs into sourmash sketches before any download
(`READS_SKETCH`, chunks of `--reads_sketch_chunk` runs, default 25, joined by
`READS_SKETCH_MERGE`), capped at `--reads_sketch_spots` spots when set.
A chunk that is retried streams only its own runs again. A chunk where
every run failed is a warning, and the merge exits 3 when no run streamed at
all. `--reads_screen` adds `READS_SCREEN`, which judges the sketches with
`--reads_screen_args` and passes only the runs of `passed.txt` to assembly.
Runs it drops, including runs whose stream failed, are listed in
`reads_screen.tsv` rather than `excused_runs.tsv`, so a transient archive
error drops a run from the batch; check that table. A screen that passes no
run stops the pipeline. Both tasks need `fastq-dump` and sourmash, and the
screen needs
`--gtdb_sketch` and `--gtdb_lineages` as parameters, since the pipeline
does not read the `REPGENR_GTDB_*` variables.

## Genome sketches

`--sketch` adds a `SKETCH` task on the bacterial path. It sketches every
genome and the outgroup, publishes `sketches/`, and passes them to the
dereplication and tree tasks as `--sketches-dir`. It is off by default.

## Split alignment and tree tasks

`--phylo_split_msa` runs the alignment (`PHYLO_MSA`, `process_high`) and
the tree (`PHYLO_TREE`, `process_medium`) as separate tasks, so another tree
builder reruns only the tree (see
[phylogeny.md](phylogeny.md#reusing-an-alignment-across-tree-builders)).
Leave it off for builders that work from genomes.

## Core-gene route

With `--msa-source coregene`, or `--preset genus` or `family`, in
`--phylo_args`, the bacterial and reads paths run one `ANNOTATE` task per
representative and the outgroup, then pass the annotations to `PHYLO`.
`--annotator` chooses bakta or prokka, `--annotate_args` adds arguments, and
`--bakta_db` and `--markers_db` stage the databases (else `BAKTA_DB` and
`REPGENR_MARKERS_DB` are read in the task). Under `stageInMode 'copy'` or on
a cloud executor without a shared file system, the bakta database would be
copied once per genome, so give the tasks a path through `BAKTA_DB` there.
Each record carries `annotation.json`, so a set annotated by two tools is
refused rather than mixed. The route is detected from
`--phylo_args`, not from a `withName` override. The viral path has no
annotation step. The annotator and builders must be reachable through
`--repgenr_opts '--bin-dir ...'` or a container profile.

## Scaling

For 10000 genomes or more, set a chunk size so the dereplication scatters
over tasks:

```bash
nextflow run nextflow/main.nf --outdir results \
    --derep_tool sourmash --derep_process_size 2000 -profile slurm -c site.config
```

`DEREPLICATE_SCATTER` dereplicates each chunk as a task, then the union once
more in one merge task, applying the keeper in both. Unlike `dereplicate
-s`, it runs one chunk round only. Seeds without an ANI, or an ANI not above
`--derep_secondary_ani`, are refused at launch. `--derep_reduce` and
`--derep_target_reps` act on the merge. `--derep_refine SEEDS
--derep_refine_ani F`, or `--derep_refine_file`, add a `DEREP_REFINE` task
(see [usage.md](usage.md#finer-resolution-inside-a-lineage)); its result
under `dereplicate/refined/` is the one the tree is built from. Resource
labels scale memory and time with the retry attempt, and their values are in
`nextflow/conf/base.config`.

## Configuring processes

Each process reads its tool flags from `task.ext.args`, which
`nextflow/conf/modules.config` maps from the parameters. A site config can
retune one process:

```groovy
// site.config
process {
    withName: 'PHYLO' {
        ext.args = '--treebuilder iqtree -B 1000'
    }
}
```

`--derep_reduce` and `--derep_target_reps` reach the merge process as
`ext.args2`, and `ext.repgenr_opts` holds the global repgenr options every
process prepends. An `ext.args` override replaces the whole string, so the dereplication
flags it does not repeat fall back to the repgenr defaults rather than to
the `--derep_*` parameters. The test harnesses are described in
[developing.md](developing.md#nf-test-stub-tests).
