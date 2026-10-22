# Installing RepGenR2

RepGenR2 is a Python package that calls external bioinformatics tools. The
package does not install those tools. This page describes how to install the
package, how to provide the tools (conda environments or containers), which
databases some tools need, and how to verify the installation. Read it before
the first run; [usage.md](usage.md) then describes how to run the pipeline.

Contents:

- [Requirements](#requirements)
- [Which method to use](#which-method-to-use)
- [The package](#the-package)
- [The external tools](#the-external-tools):
  [1. Core and satellite environments](#1-core-and-satellite-environments),
  [2. Satellites on PATH, after core](#2-satellites-on-path-after-core),
  [3. Containers](#3-containers)
- [Per-tool table](#per-tool-table)
- [Cactus](#cactus)
- [Databases](#databases)
- [Nextflow](#nextflow)
- [macOS on Apple Silicon](#macos-on-apple-silicon)
- [Verify the installation](#verify-the-installation)

## Requirements

- Python 3.12 or later for the package.
- For the tools, either conda (mamba) for the environment files under `envs/`,
  or a container engine: Docker, or Singularity or Apptainer. Containers can
  also be minted with the Seqera Wave CLI.
- For the genome download steps, `ncbi-datasets-cli` (it is in the core
  environment).
- For the Nextflow layer, Nextflow (version in [Nextflow](#nextflow)).
- Disk space for the databases of the tools that need one (see
  [Databases](#databases)) and for the working directory, where scratch data
  and downloads are written.

## Which method to use

| Situation | Recommended method | Why |
|---|---|---|
| Linux x86_64 workstation | The core conda environment plus the satellites you need, appended after core on `PATH` (sections 1 and 2), or containers | All thirteen files under `envs/` solve on linux-64. Cactus and the databases are not included. |
| HPC cluster | `--container singularity` with `--container-cache` on shared storage | Tools run from pinned images without site installs. Images are pulled once and reused. |
| macOS on Apple Silicon | The core environment (osx-arm64) with satellites appended on `PATH`, and containers for the tools that do not run natively | See [macOS on Apple Silicon](#macos-on-apple-silicon). |
| Nextflow on a cluster | A site image or conda per profile; `-profile slurm,singularity` | The `slurm` profile sets only the executor. The container profiles set `--container` for every stage. |
| Nextflow on a cloud executor | A site config with the executor, queue and an image that provides `repgenr` and the tools | No cloud profile ships, because the region, queue and image are site-specific (see [nextflow.md](nextflow.md#profiles-and-containers)). |

## The package

`pip install .` installs the distribution `repgenr2`, which provides the
`repgenr` command and the `repgenr` import package. It needs Python 3.12 or
later.

```bash
pip install .            # the package only
pip install ".[sparse]"  # with the sparse sourmash back-end
pip install -e ".[dev]"  # for development
```

Two extras exist:

- `sparse` adds the `sourmash_plugin_branchwater` back-end for sparse sourmash
  dereplication. Without it the adapter falls back to dense `sourmash compare`.
- `dev` adds the test, lint and type-check tools.

An environment installed from the project under its earlier name holds a
distribution called `repgenr` that owns the same files. Run
`pip uninstall repgenr` there before installing `repgenr2`. The core conda
environment (below) installs the package itself.

## The external tools

The tools of every family are separate programs. A stage that needs a
missing tool stops with exit code 4 and names it (the exceptions for
`assemble` are under [Verify the installation](#verify-the-installation)).
The [per-tool table](#per-tool-table) lists what each adapter needs. There
are three ways to provide the tools.

### 1. Core and satellite environments

The tools do not solve as one conda environment. Gubbins needs Python 3.8 to
3.10, mashtree and snippy need zlib older than 1.3, and progressiveMauve needs
boost 1.74. Other tools are kept apart because of their dependency sets
(CheckM2's TensorFlow, for example, is kept apart from medaka's PyTorch). The
files under `envs/` therefore define a core environment with RepGenR2 and every
tool that shares its solve, plus twelve satellite environments. Each file's
first line names its platforms and why it is separate. CI dry-runs each solve
on every change to `envs/` and once a week (`.github/workflows/envs.yml`).

| File | Environment | Contents | Platforms that solve |
|---|---|---|---|
| `envs/core.yml` | `repgenr` | Python 3.12, RepGenR2, sourmash and branchwater, skDER, skani, galah, dRep, datasets, SibeliaZ, minimap2, samtools, bcftools, ska2, skesa, shovill, flye, medaka, racon, FastTree, VeryFastTree, IQ-TREE, RAxML-NG, ASTER, rapidNJ, ClipKIT, HMMER, sra-tools (`fastq-dump` for `reads-sketch`) | linux-64, osx-arm64 |
| `envs/gubbins.yml` | `repgenr-gubbins` | Gubbins, with its own IQ-TREE 2 and RAxML | linux-64, osx-64, osx-arm64 |
| `envs/mashtree.yml` | `repgenr-mashtree` | mashtree (brings samtools 0.1.19) | linux-64, osx-64, osx-arm64 |
| `envs/snippy.yml` | `repgenr-snippy` | snippy 4.6 | linux-64, osx-64 |
| `envs/parsnp.yml` | `repgenr-parsnp` | parsnp, harvesttools | linux-64, osx-64 |
| `envs/checkm2.yml` | `repgenr-checkm2` | CheckM2 | linux-64 |
| `envs/mauve.yml` | `repgenr-mauve` | progressiveMauve (`mauvealigner`) with boost-cpp 1.74.0 | linux-64 |
| `envs/bakta.yml` | `repgenr-bakta` | bakta (needs its database) | linux-64, osx-arm64 |
| `envs/prokka.yml` | `repgenr-prokka` | prokka (no database) | linux-64, osx-64 |
| `envs/panaroo.yml` | `repgenr-panaroo` | Panaroo and mafft (its dependencies bring prokka) | linux-64, osx-64 |
| `envs/ppanggolin.yml` | `repgenr-ppanggolin` | PPanGGOLiN and MMseqs2 (its dependencies bring mafft) | linux-64, osx-64, osx-arm64 |
| `envs/scarap.yml` | `repgenr-scarap` | SCARAP, with OrthoFinder for its `O-B` and `O-D` methods (its dependencies bring MMseqs2 and mafft) | linux-64, osx-64, osx-arm64 |
| `envs/orthofinder.yml` | `repgenr-orthofinder` | OrthoFinder and DIAMOND (its dependencies bring mafft) | linux-64, osx-64, osx-arm64 |

Create the core environment and only the satellites you need. The core file
installs the package in editable mode from the checkout:

```bash
mamba env create -f envs/core.yml
mamba env create -f envs/gubbins.yml      # one line per satellite you need
```

On linux-64, dRep in the core environment brings CheckM (`checkm-genome`).
Cactus has no conda package (see [Cactus](#cactus)). The osx-arm64
differences are in [macOS on Apple Silicon](#macos-on-apple-silicon).

The versions in the files are lower bounds. They match the minimum versions
RepGenR2 checks before a run (a test keeps them in step) and guard against
known-incompatible old tools. They do not make an environment reproducible.
For a reproducible environment, export a pinned list from a solved
environment, one per platform, for example:

```bash
conda env export -n repgenr --no-builds > repgenr.lock.yml
```

(`conda-lock` serves the same purpose.)

### 2. Satellites on PATH, after core

RepGenR2 finds tools on `PATH`, and the first match wins. Activate the core
environment and append each satellite's `bin` directory after it, never
before:

```bash
mamba activate repgenr
P=$(conda info --base)/envs
export PATH="$PATH:$P/repgenr-gubbins/bin:$P/repgenr-mashtree/bin:$P/repgenr-snippy/bin"
repgenr list-tools --check
```

Appended satellites supply only what core lacks. Core's samtools and
bcftools then come first, ahead of the samtools 0.1.19 in the mashtree
environment, which the `simple` SNP typer would reject; a satellite placed
before core would shadow them.

Two known gaps remain, because a satellite tool calls its helpers by name and
core's come first:

- snippy calls core's samtools and bcftools instead of the versions in its own
  environment. This combination has not been tested.
- Gubbins calls `iqtree` by name, so it runs core's IQ-TREE 3 instead of the
  IQ-TREE 2.4 of its own environment. On a small test set the result was the
  same, but the Gubbins log then records no IQ-TREE version ("for" in its place).

The per-tool `--bin-dir` option closes both gaps. Instead of appending,
each satellite is given to its tool alone with `--bin-dir TOOL=DIR` (or
`REPGENR_BIN_DIRS`). The directory then comes first on `PATH` for that tool
and its subprocesses only, so snippy, for example, runs with the samtools and
bcftools of its own environment:

```bash
export REPGENR_BIN_DIRS="gubbins=$P/repgenr-gubbins/bin,mashtree=$P/repgenr-mashtree/bin,snippy=$P/repgenr-snippy/bin"
repgenr list-tools --check
```

See [usage.md](usage.md#per-tool-directories---bin-dir) for the details.

A single command can also be run inside a named environment without
activating it:

```bash
conda run -n repgenr --no-capture-output repgenr phylo -wd $WD --treebuilder iqtree
```

### 3. Containers

RepGenR2 can run each tool inside a pinned container instead of from `PATH`.
The package itself runs on the host and only the tool subprocess is
containerized, so this works in the plain CLI and inside Nextflow. The default
is `none` (native execution).

```bash
repgenr --container docker dereplicate -wd $WD --tool drep
repgenr --container singularity --container-cache /data/repgenr_sif phylo -wd $WD
```

The options are global and go before the subcommand. Each has an environment
variable: `REPGENR_CONTAINER`, `REPGENR_CONTAINER_ENGINE`,
`REPGENR_CONTAINER_CACHE`, `REPGENR_CONTAINER_PLATFORM` and `REPGENR_WAVE`.

- `--container docker|singularity` selects the backend.
  `--container-engine` replaces the binary (apptainer, podman).
- `--container-cache <dir>` is where Singularity images and their cache
  live. Put it on a large disk. Docker ignores it, and Wave image references
  are resolved once per process, not stored. A container option given without
  `--container` (or its variable set while running natively) has no effect and
  is named in a warning.
- `--platform linux/amd64` runs amd64 images on an arm64 host (see
  [macOS on Apple Silicon](#macos-on-apple-silicon)). On a native amd64
  Linux host no emulation is involved.
- `--wave` mints an image from the adapter's conda spec with the Seqera Wave
  CLI. It is the only image route for the four adapters that carry no pinned
  image: `simple`, parsnp, skder and SibeliaZ. Without `--wave` those four run
  on the host and the log says so.

#### Image sources

Each adapter declares a pinned image, a conda spec, or both:

- The pinned image is used when Wave is off. Most adapters pin a
  `linux/amd64` BioContainer; Cactus pins its project image, and sra-tools
  pins NCBI's `ncbi/sra-tools`, built for amd64 and arm64. The pins are
  listed in the [per-tool table](#per-tool-table).
- With `--wave`, RepGenR2 mints an arm64-native image from the conda spec.
  This is the only image route for `simple` and parsnp, which span several
  packages, and for skder and SibeliaZ. Their BioContainers are
  BusyBox-based, and the GNU-only calls in their shell wrappers fail there,
  so the run yields nothing.

#### Storage location

- **Singularity or Apptainer.** `--container-cache` sets `APPTAINER_CACHEDIR`
  and `SINGULARITY_CACHEDIR` (and the matching `*_TMPDIR`). `docker://` images
  are pulled once to `<cache>/<name>.sif` and reused. Put this on a large or
  external disk.
- **Docker.** The daemon manages image storage (Docker Desktop's disk image).
  It is set there, not per run.
- Large run-time data (the Cactus job store, scratch, downloads) lives under
  `--workdir` and `TMPDIR`.

#### Checking the container setup

- Before a tool runs in an image, RepGenR2 checks that the engine binary is
  on `PATH` and, for Docker, that `docker info` reaches the daemon. A stopped
  Docker Desktop therefore exits 4 before the stage starts.
- `repgenr --container docker list-tools --check` names where each tool
  runs: `[image <ref>]` for a pinned (or `--wave`) image and `[host]` for a
  tool that runs from `PATH`. `--images` adds whether each image is present
  locally (`present`, `not pulled` or `presence unknown`), including
  secondary images such as racon's minimap2 image. Nothing is pulled; an
  image is pulled on first use.
- dRep's CheckM step needs its reference data. Set `CHECKM_DATA_PATH` on the
  host: RepGenR2 passes it into the container and binds that directory at
  the same path. Otherwise run dRep with `--ignoreGenomeQuality`.
- progressiveMauve from bioconda (`mauvealigner`) is linked against boost
  1.74 and fails with an `undefined symbol` error against a newer boost.
  `envs/mauve.yml` and the adapter's Wave spec therefore pin
  `boost-cpp=1.74.0`; the pinned image already holds a compatible boost.

At run time, Docker runs as the host UID and GID, so outputs are owned by
you, and the working directory and `TMPDIR` are bind-mounted at identical
paths. Symlinked inputs, such as a `genomes/` directory staged by
`repgenr ingest`, are followed, and the directory each link points to is
bound as well. A stage that ran a tool in an image records the image
reference as that tool's version, with the engine and its version beside
it, in `repgenr.yaml` and in `repgenr versions`.

## Per-tool table

The table lists the adapters in the order a run uses them; CheckM2, which
scores assemblies in the reads chain, sits with the read tools. "Native on
macOS arm64" records what was run natively on the audit machine
([verification.md](verification.md#adapter-status)); "no" means not run
there, not that the tool cannot work. `none` in the pin column means the
adapter needs `--wave` to run in an image. Minimum versions are in the
`envs/` files, the environment of each tool is in the table of
[section 1](#1-core-and-satellite-environments), and the databases are
described under [Databases](#databases).

| Family | Tool | conda package | Native on macOS arm64 | Container pin | Databases or extras |
|---|---|---|---|---|---|
| Assembler | skesa | `skesa` | yes | `skesa:2.5.1` | none |
| Assembler | shovill | `shovill` | no | `shovill:1.4.2` | none |
| Assembler | flye | `flye` | yes | `flye:2.9.6` | none |
| Polisher | medaka | `medaka` | no | `medaka:2.2.2` | medaka models |
| Polisher | racon | `racon`, `minimap2` | yes | `racon:1.5.0`; minimap2 runs from `minimap2:2.28` | none |
| Classifier | sourmash | `sourmash` | yes | `sourmash:4.9.4` | GTDB sourmash sketch and lineages |
| Quality (reads chain) | checkm2 | `checkm2` | no | `checkm2:1.1.0` | CheckM2 database |
| Dereplicator | drep | `drep` (brings `checkm-genome` on linux-64) | no | `drep:3.7.1` | CheckM reference data |
| Dereplicator | skder | `skder`, `skani` | yes | none (Wave) | none |
| Dereplicator | galah | `galah` | yes | `galah:0.4.2` | none |
| Dereplicator | sourmash | `sourmash` (4.4 or later, for `compare --avg-containment`) | yes | `sourmash:4.9.4` | `sourmash_plugin_branchwater` for the sparse back-end (extra `sparse`) |
| Annotator | bakta | `bakta` | yes | `bakta:1.12.2` | bakta database |
| Annotator | prokka | `prokka` | yes, osx-64 under Rosetta | `prokka:1.15.6` | none |
| Aligner | progressivemauve | `mauvealigner`, `boost-cpp=1.74.0` | no (unpackaged on macOS) | `mauve:2.4.0.snapshot_2015_02_13` | none |
| Aligner | sibeliaz | `sibeliaz` | yes | none (Wave) | none |
| Aligner | cactus | not on conda | no | `cactus:v2.9.3` | none |
| SNP typer | simple | `minimap2`, `samtools`, `bcftools` | yes | none (Wave) | samtools and bcftools 1.10 or later |
| SNP typer | snippy | `snippy` | no | `snippy:4.6.0` | none |
| SNP typer | parsnp | `parsnp`, `harvesttools` | yes, osx-64 under Rosetta | none (Wave) | none |
| SNP typer | ska2 | `ska2` | yes | `ska2:0.5.1` | none |
| Masker | gubbins | `gubbins` | yes | `gubbins:3.4.3` | a multi-threaded RAxML build for more than one thread (see [phylogeny.md](phylogeny.md#recombination-masking-with-gubbins)) |
| Core-gene builder | panaroo | `panaroo`, `mafft` | yes, osx-64 under Rosetta | `panaroo:1.8.0` | none |
| Core-gene builder | ppanggolin | `ppanggolin`, `mmseqs2` (mafft comes with its dependencies) | yes | `ppanggolin:2.3.2` | none |
| Core-gene builder | scarap | `scarap`, `orthofinder`, `pandas<2` (MMseqs2 and mafft come with its dependencies) | yes | `scarap:1.0.2` | none |
| Core-gene builder | orthofinder | `orthofinder`, `diamond` (mafft comes with its dependencies) | yes | `orthofinder:3.1.5` | none |
| Core-gene builder | markers | `hmmsearch`, `hmmalign` (hmmer, in `envs/core.yml`) | yes | `hmmer:3.4` | marker database |
| Core-gene trimming (`--gene-trim clipkit`) | clipkit | `clipkit` (in `envs/core.yml`) | yes | `clipkit:2.14.0` | none |
| Tree builder | iqtree | `iqtree` | yes | `iqtree:3.1.3` | none |
| Tree builder | fasttree | `fasttree` | yes | `fasttree:2.2.0` | none |
| Tree builder | veryfasttree | `veryfasttree` | no record (offline contract test only) | `veryfasttree:4.0.5` | none |
| Tree builder | raxmlng | `raxml-ng` | yes | `raxml-ng:2.0.3` | none |
| Tree builder | mashtree | `mashtree` | yes | `mashtree:1.4.6` | none |
| Tree builder | sourmash | `sourmash` | yes | `sourmash:4.9.4` | optional `rapidnj`; without it more than 5000 leaves are refused |
| Tree builder | skani | `skani` | no record (offline contract test only) | `skani:0.3.2` | optional `rapidnj`, for the joining step above a few thousand leaves |
| Tree builder | ska2 | `ska2` | no record (offline contract test only) | `ska2:0.5.1` | optional `rapidnj`, as for skani |
| Tree builder | aster | `aster` | yes | `aster:1.25` | IQ-TREE or FastTree for the gene trees |

Under a container backend rapidNJ runs from `rapidnj:2.3.3`. `clipkit` is
in `envs/core.yml` and is found by its own name (`--bin-dir clipkit=DIR`),
not through a builder's environment.

## Cactus

Cactus has no conda package in `envs/` and no native install path
in these docs. The adapter pins the project image
(`quay.io/comparative-genomics-toolkit/cactus:v2.9.3`), so use
`--container docker|singularity`. On Apple Silicon this needs the Rosetta
setting described under [macOS on Apple Silicon](#macos-on-apple-silicon).

## Databases

Most tools need no database. The following do.

| Tool | Database | How to provide it |
|---|---|---|
| bakta (annotator) | bakta database; the light type is a 1.4 GB download, about 4 GB unpacked, and the full type 30 GB, about 84 GB unpacked, with more specific CDS annotation | `bakta_db download --output DIR --type light`, then `--bakta-db DIR` or `BAKTA_DB` |
| hmmer (markers core-gene builder) | GTDB marker models (bac120, ar53; about 25 MB) | `repgenr markers-db download --output DIR`, then `REPGENR_MARKERS_DB` |
| checkm2 (reads chain) | CheckM2 database | `checkm2 database --download`, then `assemble --checkm2-db` or `CHECKM2DB` |
| drep (dereplicator) | CheckM reference data | `CHECKM_DATA_PATH`, or run dRep with `--ignoreGenomeQuality` |
| sourmash (classifier) | a GTDB sourmash sketch and its lineages CSV | `assemble --gtdb-sketch` and `--gtdb-lineages` (or `REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`) |
| medaka (polisher) | medaka models | for SRA reads without a basecaller in their headers the adapter assumes ONT's bacterial R10.4.1 model; another is given with `--tool-arg model=...` |

prokka annotates without a database.

## Nextflow

The Nextflow layer needs Nextflow 26.04 or later (`nextflowVersion =
'!>=26.04.0'` in `nextflow/nextflow.config`). The nf-test suite uses nf-test
0.9.5, the version CI installs. The `docker`, `singularity` and `wave`
profiles set `params.repgenr_opts`, so every stage calls `repgenr` with
`--container docker`, `--container singularity` or `--container singularity
--wave`. The `repgenr` command itself must be available to each process. See
[nextflow.md](nextflow.md) for parameters and profiles.

## macOS on Apple Silicon

The core environment and six of the twelve satellites solve for osx-arm64.
The other tools need Rosetta, a container or a Linux host:

- **Native runs.** The [per-tool table](#per-tool-table) records which
  adapters ran natively on the audit machine.
- **osx-64 satellites.** snippy, parsnp, prokka and Panaroo have no
  osx-arm64 build. Create their environments with `CONDA_SUBDIR=osx-64` and
  run them under Rosetta.
- **Containers.** progressiveMauve, Cactus, snippy, dRep, shovill, medaka
  and CheckM2 were run in containers. progressiveMauve is not packaged for
  macOS, and the CheckM2 host builds fail there.
- **dRep.** The osx-arm64 solve gives dRep 3.4.5 without CheckM. Its fastANI
  step fails with pandas releases that no longer accept `delim_whitespace`
  in `read_csv` (seen with pandas 3.0.6), while its ANImf method (`--virus`)
  runs. Run it with `--ignoreGenomeQuality`, or from its image, which holds
  dRep 3.7.1.
- **Images.** BioContainers are `linux/amd64`. Pass `--platform linux/amd64`
  with Docker Desktop's Apple Virtualization framework and "Use Rosetta for
  x86/amd64 emulation" enabled, or use `--wave` for arm64-native images. The
  sra-tools image runs natively. Docker's QEMU emulation cannot run some
  SIMD-heavy binaries: the `vg` bundled in Cactus hangs even on
  `vg version`.
- **exFAT and NTFS volumes.** macOS writes `._*` AppleDouble files there.
  The adapters ignore them, but skDER, dRep and ParSNP must run on a local
  filesystem, since dRep writes such files into its own cache and ParSNP
  reads every file in its input directory.

The SibeliaZ adapter applies a workaround for the BSD tools of macOS, which
it skips when SibeliaZ runs in a Linux container.

## Verify the installation

```bash
repgenr --version
repgenr list-tools --check
```

`list-tools` prints each family with its adapters and, where one is declared,
the genome limit. Its last line names the dereplicators that `glance` can run
(`drep`, `sourmash`).

`--check` runs each adapter's preflight and prints `ok` with the versions
found, `missing` or `error` with the reason, or `broken` for a plugin that
failed to load. It only reports, and exits with status 0. Add `--strict` to
verify an environment from a script: after the full listing it exits 4 when
any adapter is missing or errored, and 5 when any plugin is broken.

A stage run with a missing or outdated tool exits with status 4. The
exception is `assemble --assembler auto`, which excuses runs whose assembler
is missing (reason `assembler_not_installed`, with a warning) and exits 4 only
when no run can be assembled. A missing `auto` polisher is skipped (see
[reads.md](reads.md#polishing-long-read-assemblies)). The full exit-code table
is in [usage.md](usage.md#troubleshooting).
