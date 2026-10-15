# RepGenR2

RepGenR2 builds dereplicated, representative genome repositories for
large-scale genomic studies. It selects taxa, downloads genomes, clusters them
by average nucleotide identity (ANI), computes phylogenetic trees, and emits
taxonomy files for downstream tools such as FlexTaxD.

It is the second generation of RepGenR (Representative-Genome Repositories):
a modular, importable Python package (Python 3.12+, distributed as `repgenr2`,
installed as the `repgenr` command) with eight pluggable tool families and an
optional Nextflow pipeline for scatter-gather runs on HPC and cloud. The
largest end-to-end runs recorded are 1157 bacterial genomes (one genus, 12
minutes) and 1256 viral genomes (Hepeviridae). Larger sizes are measured per
tool in the
[scaling audit](docs/audit/scaling-audit.md).

## Pipeline

```
bacterial:         metadata -> genome -> dereplicate -> phylo -> tree2tax
viral:             vmetadata -> vgenome -> dereplicate -> phylo -> tree2tax
local genomes:     ingest -> dereplicate -> phylo -> tree2tax
sequencing reads:  reads -> assemble -> dereplicate -> phylo -> tree2tax
```

`repgenr run` chains a whole row; `repgenr status -wd WD` says what comes next.

`phylo` builds its tree from an alignment it produces with an aligner or a
SNP typer, or directly from the genomes with an alignment-free builder;
`repgenr snptype` runs the SNP typer on its own when the SNP tables are the
deliverable. Each stage is a `repgenr` subcommand and an importable function. Stages
communicate through a single working directory whose state is recorded in
`repgenr.yaml` (provenance) and a SQLite genome manifest.

## Pluggable tools

Tools are discovered through Python entry points; adding one needs no change to
the core (see [docs/developing.md](docs/developing.md)).

| Family | Built-in adapters |
|--------|-------------------|
| Dereplicators | `drep`, `skder`, `galah`, `sourmash` |
| Aligners | `progressivemauve`, `sibeliaz`, `cactus` |
| SNP typers | `simple` (samtools/bcftools), `snippy`, `parsnp`, `ska2` (reference-free) |
| Maskers | `gubbins` (recombination masking via `--mask`) |
| Annotators | `bakta` (default; needs its database), `prokka`; used by `annotate` for the core-gene route |
| Tree builders | `iqtree`, `fasttree`, `veryfasttree`, `raxmlng` (MSA), `mashtree`, `sourmash`, `skani`, `ska2` (alignment-free) |
| Assemblers | `skesa`, `shovill` (Illumina), `flye` (ONT, PacBio); used by the reads chain |
| Classifiers | `sourmash` (gather against a GTDB sketch); verifies the organism of an assembled run |
| Polishers | `medaka` (ONT), `racon` with minimap2 (PacBio CLR); correct long-read assemblies in the reads chain |

`repgenr list-tools` prints what is available in your environment.

## Installation

```bash
pip install .                          # the package (Python 3.12+); tools are separate
repgenr list-tools --check             # which tools are found, with versions
```

The tools come from conda environments or from containers. `envs/core.yml`
holds RepGenR2 and most tools; six satellite files (`envs/gubbins.yml`,
`envs/mashtree.yml` and others) hold tools that do not solve with it. Append a
satellite's `bin` directory after core on `PATH`, never before it:

```bash
mamba env create -f envs/core.yml && mamba activate repgenr
mamba env create -f envs/gubbins.yml
export PATH="$PATH:$(conda info --base)/envs/repgenr-gubbins/bin"
```

The environments, containers (including Apple Silicon) and the databases some
tools need are described in [docs/install.md](docs/install.md).

## Quick start

```bash
WD=./francisella
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` chains `metadata -> genome -> dereplicate -> phylo -> tree2tax`; each stage
is also its own subcommand, and `repgenr ingest` starts from genomes already on
disk. `--viral` selects the NCBI Virus path. Re-running a stage is a no-op
unless its parameters or inputs changed (`--force` overrides).

The same pipeline runs as typed Nextflow data channels (no shared working
directory) for HPC and cloud:

```bash
nextflow run nextflow/main.nf -profile standard --outdir results \
    --metadata_args "-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella" \
    --derep_tool sourmash --phylo_args "--treebuilder mashtree"
```

`--mode viral` and `--mode reads` (ENA/SRA sequencing runs, one assembly task
per run) select the other two front ends. Nextflow 26.04 or later is required.

## Documentation

| Page | What it covers |
|------|----------------|
| [docs/install.md](docs/install.md) | Installing the package and the external tools: conda, several environments, containers, per-tool requirements and databases. |
| [docs/choosing-tools.md](docs/choosing-tools.md) | Which dereplicator, phylogeny route and tree builder to use for a dataset and size, with the measured runs and declared limits behind each. |
| [docs/usage.md](docs/usage.md) | Running the pipeline: CLI stages, local genomes, viruses, resume, representative selection, SNP typing and masking, Nextflow parameters and profiles, containers, troubleshooting. |
| [docs/cli-reference.md](docs/cli-reference.md) | Every command and option, generated from the command tree. |
| [docs/output.md](docs/output.md) | The files each stage writes. |
| [docs/developing.md](docs/developing.md) | Architecture, data contracts, and how to add a tool adapter. |
| [docs/verification.md](docs/verification.md) | Which adapters have been run against their real tools, and measured runs. |
| [docs/audit/](docs/audit/README.md) | Records: the CLI matrix and the scaling and bias audit. |

## License

MIT.
