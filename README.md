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
executors. The largest end-to-end runs recorded are 1157
bacterial genomes (one genus, 12 minutes) and 1256 viral genomes
(Hepeviridae). Larger sizes are measured per tool in the
[scaling audit](docs/audit/scaling-audit.md).

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

RepGenR2 runs external tools through adapters. The adapters belong to ten
tool families, and each family can be swapped without changing the rest of
the pipeline. The families fall into four groups, by what they work on and
in the order a run uses them. A fifth group, background data, is built in
rather than pluggable. Adapters are discovered through Python entry points,
so adding one needs no change to the core (see
[docs/developing.md](docs/developing.md)).

**Background data** (built in, not pluggable): where taxa, genomes and runs
come from, and the reference databases some tools read.

| Kind | Sources |
|------|---------|
| Taxonomy and genome metadata | GTDB (full table or API) for bacteria and archaea, including the species-representative flag that `--keeper gtdb` reads; NCBI Virus (default) or BV-BRC for viruses |
| Genome sequences | NCBI Datasets; local FASTA files (`ingest`) |
| Sequencing runs | ENA (runs mirrored from SRA) |
| Reference databases | bakta database (annotation), CheckM2 database (assembly quality), GTDB sourmash sketch and lineages (organism check of assemblies) |

**Read tools**: sequencing runs to checked genome assemblies (`reads`,
`assemble`).

| Family | Built-in adapters |
|--------|-------------------|
| Assemblers | `skesa`, `shovill` (Illumina), `flye` (ONT, PacBio) |
| Polishers | `medaka` (ONT), `racon` with minimap2 (PacBio CLR); correct long-read assemblies |
| Classifiers | `sourmash` (gather against a GTDB sketch); confirm the organism of each assembly |

**Genome tools**: a genome set to representatives and annotations
(`dereplicate`, `annotate`).

| Family | Built-in adapters |
|--------|-------------------|
| Dereplicators | `drep`, `skder`, `galah`, `sourmash` |
| Annotators | `bakta` (default; needs its database), `prokka`; for the core-gene route |

**Alignment tools**: genomes to the multiple sequence alignment (MSA) a tree
is built from (`phylo`, `snptype`).

| Family | Built-in adapters |
|--------|-------------------|
| Aligners | `progressivemauve`, `sibeliaz`, `cactus` (whole-genome alignment) |
| SNP typers | `simple` (samtools/bcftools), `snippy`, `parsnp`, `ska2` (reference-free) |
| Maskers | `gubbins` (recombination masking of the SNP alignment, `--mask`) |
| Core-gene builders | `panaroo` (with mafft), `ppanggolin`, `scarap`, `orthofinder` (`--msa-source coregene`) |

**Tree tools**: an alignment, or the genomes themselves, to a tree (`phylo`).

| Family | Built-in adapters |
|--------|-------------------|
| Tree builders | from an MSA: `iqtree`, `fasttree`, `veryfasttree`, `raxmlng`; alignment-free, from the genomes: `mashtree`, `sourmash`, `skani`, `ska2` |

`repgenr list-tools` prints what is available in the current environment.
[docs/choosing-tools.md](docs/choosing-tools.md) describes which tools suit a
dataset and its size.

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
per run) select the other two front ends. Nextflow 26.04 or later is required.
Parameters and profiles are described in [docs/usage.md](docs/usage.md#nextflow).

## Documentation

| Page | What it covers |
|------|----------------|
| [docs/install.md](docs/install.md) | Installing the package and the external tools: conda environments, containers, per-tool requirements and databases. |
| [docs/choosing-tools.md](docs/choosing-tools.md) | Which dereplicator, phylogeny route and tree builder to use for a dataset and size, with the measured runs and declared limits behind each. |
| [docs/usage.md](docs/usage.md) | Running the pipeline: CLI stages, local genomes, viruses, resume, representative selection, SNP typing and masking, Nextflow parameters and profiles, containers, troubleshooting. |
| [docs/cli-reference.md](docs/cli-reference.md) | Every command and option, generated from the command tree. |
| [docs/output.md](docs/output.md) | The files each stage writes. |
| [docs/developing.md](docs/developing.md) | Architecture, data contracts, and how to add a tool adapter. |
| [docs/verification.md](docs/verification.md) | Which adapters have been run against their real tools, and measured runs. |
| [docs/audit/](docs/audit/README.md) | Records: the CLI matrix, the scaling and bias audit, and the phylo and tree2tax gap analysis. |

The changelog is [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use RepGenR2, cite this repository
(<https://github.com/FOI-Bioinformatics/repgenr2>) and the external tools your
run used. [CITATIONS.md](CITATIONS.md) lists the references for the tools,
databases and workflow manager.

## License

MIT; see [LICENSE](LICENSE).
