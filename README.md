# RepGenR2

RepGenR2 builds representative genome repositories for a taxonomic group,
together with a phylogeny and a taxonomy for FlexTaxD. This page gives an
overview of what the software does, how it is installed and run, and where the
detailed documentation is.

## Purpose

Public genome collections are uneven: some species and genotypes are
represented by many near-identical assemblies, others by one. A reference
database built directly from such a collection is larger than it needs to be
and weighted towards the over-represented groups. A representative genome set
keeps one genome per cluster of closely related genomes and records the other
genomes as members of that cluster.

RepGenR2 builds such a set in five steps:

1. It selects the genomes of a taxon: from GTDB for bacteria and archaea, from
   NCBI Virus for viruses, from local files, or from sequencing reads that it
   assembles.
2. It downloads or stages the genomes.
3. It dereplicates them. Genomes are clustered by average nucleotide identity
   (ANI), and one representative genome is chosen per cluster.
4. It builds a phylogeny of the representatives.
5. It converts the tree into a taxonomy for FlexTaxD. By default the
   dereplicated genomes are listed as well.

RepGenR2 is the second generation of RepGenR (Representative-Genome
Repositories). It is an importable Python package (Python 3.12 or later,
distributed as `repgenr2`, installed as the `repgenr` command) with ten
pluggable tool families, grouped as read, genome, alignment and tree tools,
and an optional Nextflow pipeline for scatter-gather runs on HPC and cloud
executors. The largest end-to-end runs recorded are 1157 bacterial genomes
(one genus, 12 minutes) and 1256 viral genomes (Hepeviridae). Larger sizes are
measured per tool in the [scaling audit](docs/audit/scaling-audit.md).

## Pipeline routes

Each row below is one route through the pipeline. Each name in a row is a
stage: a `repgenr` subcommand and an importable function.

```
bacterial:         metadata -> genome -> dereplicate -> phylo -> tree2tax
viral:             vmetadata -> vgenome -> dereplicate -> phylo -> tree2tax
local genomes:     ingest -> dereplicate -> phylo -> tree2tax
sequencing reads:  reads -> assemble -> dereplicate -> phylo -> tree2tax
```

`repgenr run` chains a whole row. `repgenr status -wd WD` reports which stages
are done and what comes next. Stages communicate through a single working
directory (given with `-wd`). Its state is recorded in `repgenr.yaml`
(provenance) and in a SQLite genome manifest.

`phylo` builds the tree by one of four routes:

- from a whole-genome multiple sequence alignment (MSA) produced by an aligner;
- from a core-SNP alignment produced by a SNP typer (`--msa-source snptype`),
  optionally with recombination masking (`--mask gubbins`);
- from a core-gene alignment: the `annotate` stage annotates the genomes, then
  `phylo --msa-source coregene` aligns the genes shared by the genomes (the
  core genome) with a core-gene builder;
- directly from the genomes, with an alignment-free tree builder.

`repgenr snptype` runs the SNP typer on its own when the SNP tables are the
deliverable. [docs/choosing-tools.md](docs/choosing-tools.md#phylogeny-routes)
compares the routes.

## Tool families

RepGenR2 runs external tools through adapters. Each of the ten tool families
can be swapped without changing the rest of the pipeline, and adapters are
discovered through Python entry points, so adding one needs no change to the
core (see [docs/developing.md](docs/developing.md)). The families, in the
order a run uses them:

- **Assemblers:** `skesa`, `shovill` (Illumina), `flye` (ONT, PacBio).
- **Polishers:** `medaka` (ONT), `racon` (PacBio CLR).
- **Classifiers:** `sourmash`, which confirms the organism of each assembly.
- **Dereplicators:** `skder` (default), `galah`, `sourmash`, `drep`.
- **Annotators:** `bakta` (default), `prokka`, for the core-gene route.
- **Aligners:** `progressivemauve`, `sibeliaz`, `cactus`.
- **SNP typers:** `simple`, `snippy`, `parsnp`, `ska2`.
- **Maskers:** `gubbins` (`--mask`).
- **Core-gene builders:** `panaroo`, `ppanggolin`, `scarap`, `orthofinder`,
  `markers` (GTDB marker proteins).
- **Tree builders:** `iqtree`, `fasttree`, `veryfasttree`, `raxmlng` and
  `aster` from an alignment; `mashtree`, `sourmash`, `skani` and `ska2` from
  the genomes.

The data sources are built in rather than pluggable: GTDB for bacteria and
archaea, NCBI Virus (default) or BV-BRC for viruses, NCBI Datasets for
genome sequences, and ENA for sequencing runs.
`repgenr list-tools` prints what is available in the current environment.
[docs/choosing-tools.md](docs/choosing-tools.md) compares the tools for a
dataset and its size, and [docs/install.md](docs/install.md) lists the
environments and databases they need.

## Installation

The package and the external tools are installed separately:

```bash
pip install .                          # the package (Python 3.12+); tools are separate
repgenr list-tools --check             # which tools are found, with versions
```

The tools come from conda environments or from containers. `envs/core.yml`
holds RepGenR2 and most tools. Twelve satellite files (`envs/gubbins.yml`,
`envs/mashtree.yml` and others) hold tools that do not solve with it. Append a
satellite's `bin` directory after core on `PATH`, never before it:

```bash
mamba env create -f envs/core.yml && mamba activate repgenr
mamba env create -f envs/gubbins.yml
export PATH="$PATH:$(conda info --base)/envs/repgenr-gubbins/bin"
```

[docs/install.md](docs/install.md) describes the environments, containers
(including Apple Silicon) and the databases some tools need.

## Quick start

```bash
WD=./francisella
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` chains `metadata -> genome -> dereplicate -> phylo -> tree2tax`. Each
stage is also its own subcommand. `repgenr ingest` starts from genomes already
on disk, and `--viral` selects the NCBI Virus route. Re-running a stage does
nothing unless its parameters or inputs changed; `--force` overrides this.

## Nextflow

The same pipeline runs as typed Nextflow data channels, without a shared
working directory, for HPC and cloud executors:

```bash
nextflow run nextflow/main.nf -profile standard --outdir results \
    --metadata_args "-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella" \
    --derep_tool sourmash --phylo_args "--treebuilder mashtree"
```

`--mode viral` and `--mode reads` (ENA/SRA sequencing runs, one assembly task
per run) select the other two front ends. The required Nextflow version is
in [docs/install.md](docs/install.md#nextflow), and the parameters and
profiles are described in [docs/nextflow.md](docs/nextflow.md).

## Documentation

| Page | What it covers |
|------|----------------|
| [docs/install.md](docs/install.md) | Installing the package and the external tools: conda environments, containers, per-tool requirements and databases. |
| [docs/choosing-tools.md](docs/choosing-tools.md) | Which dereplicator, phylogeny route and tree builder to use for a dataset and size, with the measured runs and declared limits behind each. |
| [docs/usage.md](docs/usage.md) | Running the pipeline from the command line: stages, resume, GTDB, local genomes, viruses, representative selection, a summary of the phylogeny routes, tree2tax, containers, troubleshooting. |
| [docs/reads.md](docs/reads.md) | Building the genome set from sequencing runs: run selection, assembly, polishing, quality checks and screening. |
| [docs/phylogeny.md](docs/phylogeny.md) | The phylogeny routes in depth: presets, SNP typing and masking, annotation, the core-gene builders, support values and rooting. |
| [docs/nextflow.md](docs/nextflow.md) | The Nextflow pipeline: structure, parameters, profiles, reads mode, scaling and process configuration. |
| [docs/cli-reference.md](docs/cli-reference.md) | Every command and option, generated from the command tree. |
| [docs/output.md](docs/output.md) | The files each stage writes. |
| [docs/developing.md](docs/developing.md) | Architecture, data contracts, and how to add a tool adapter. |
| [docs/verification.md](docs/verification.md) | Which adapters have been run against their real tools, and measured runs. |
| [docs/audit/](docs/audit/README.md) | Records: the CLI matrix, the command audit, the scaling and bias audit, the phylo methods audit and the phylo and tree2tax gap analysis. |

The changelog is [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use RepGenR2, cite this repository
(<https://github.com/FOI-Bioinformatics/repgenr2>) and the external tools your
run used. [CITATIONS.md](CITATIONS.md) lists the references for the tools,
databases and workflow manager.

## License

MIT; see [LICENSE](LICENSE).
