# Installing RepGenR2

RepGenR2 is a Python package that calls external bioinformatics tools. The
package does not install those tools. This page covers both parts and three
ways to provide the tools.

## What gets installed

- **The package.** `pip install .` installs the distribution `repgenr2`, which
  provides the `repgenr` command and the `repgenr` import package. It needs
  Python 3.12 or later. An environment installed from the project under its
  earlier name holds a distribution called `repgenr` that owns the same files;
  run `pip uninstall repgenr` before installing `repgenr2` there. Two extras exist: `sparse` adds the
  `sourmash_plugin_branchwater` back-end for sparse sourmash dereplication
  (without it the adapter falls back to dense `sourmash compare`), and `dev`
  adds the test, lint and type-check tools.
- **The tools.** Dereplicators, aligners, SNP typers, the masker, tree
  builders, assemblers and polishers are separate programs. A stage that needs
  one that is missing stops with exit code 4 and names it, except under `assemble --assembler auto`, which excuses runs whose assembler is missing (reason `assembler_not_installed`, with a warning) and exits 4 only when no run can be assembled. A missing `auto` polisher is skipped (see [usage.md](usage.md#starting-from-sequencing-reads)). See
  [the per-tool table](#per-tool-table) for what each adapter needs.

```bash
pip install .            # the package only
pip install ".[sparse]"  # with the sparse sourmash back-end
pip install -e ".[dev]"  # for development
```

## Which method to use

| Situation | Recommended method | Why |
|---|---|---|
| Linux x86_64 workstation | The core conda environment plus the satellites you need, appended after core on `PATH` (sections 1 and 2), or containers | All seven files under `envs/` solve on linux-64. Cactus and the databases are not included. |
| HPC cluster | `--container singularity` with `--container-cache` on shared storage | Tools run from pinned images without site installs. Images are pulled once and reused. |
| macOS on Apple Silicon | The core environment (osx-arm64) with satellites appended on `PATH`; containers for progressiveMauve, Cactus, snippy, dRep, shovill, medaka and CheckM2 | These were run inside containers on the audit machine (`verification.md`); the CheckM2 host builds fail on macOS. skder, galah, sourmash, SibeliaZ, `simple`, parsnp, ska2, Gubbins, the tree builders, skesa, flye and racon ran natively. |
| Nextflow on a cluster | A site image or conda per profile; `-profile slurm,singularity` | The `slurm` profile sets only the executor. The container profiles set `--container` for every stage. |
| Nextflow on a cloud executor | A site config with the executor, queue and an image that provides `repgenr` and the tools | No cloud profile ships, because the region, queue and image are site-specific (see [usage.md](usage.md#profiles)). |

## Three ways to provide the tools

### 1. Core and satellite environments

The tools do not solve as one conda environment. Gubbins needs Python 3.8 to
3.10, mashtree and snippy need zlib older than 1.3, and progressiveMauve needs
boost 1.74. The files under `envs/` therefore define a core environment with
RepGenR2 and every tool that shares its solve, plus six satellite environments.
Each file's first line names its platforms and why it is separate. CI
dry-runs each solve on every change to `envs/` and once a week
(`.github/workflows/envs.yml`).

| File | Environment | Contents | Platforms that solve |
|---|---|---|---|
| `envs/core.yml` | `repgenr` | Python 3.12, RepGenR2, sourmash and branchwater, skDER, skani, galah, dRep, datasets, SibeliaZ, minimap2, samtools, bcftools, ska2, skesa, shovill, flye, medaka, racon, FastTree, IQ-TREE, RAxML-NG | linux-64, osx-arm64 |
| `envs/gubbins.yml` | `repgenr-gubbins` | Gubbins (Python 3.10, with its own IQ-TREE 2 and RAxML) | linux-64, osx-64, osx-arm64 |
| `envs/mashtree.yml` | `repgenr-mashtree` | mashtree (brings samtools 0.1.19) | linux-64, osx-64, osx-arm64 |
| `envs/snippy.yml` | `repgenr-snippy` | snippy 4.6 | linux-64, osx-64 |
| `envs/parsnp.yml` | `repgenr-parsnp` | parsnp, harvesttools | linux-64, osx-64 |
| `envs/checkm2.yml` | `repgenr-checkm2` | CheckM2 (TensorFlow, kept apart from medaka's PyTorch) | linux-64 |
| `envs/mauve.yml` | `repgenr-mauve` | progressiveMauve (`mauvealigner`) with boost-cpp 1.74.0 | linux-64 |
| `envs/bakta.yml` | `repgenr-bakta` | bakta (annotation for the core-gene route; needs its database) | linux-64, osx-arm64 |
| `envs/prokka.yml` | `repgenr-prokka` | prokka (annotation without a database; Perl, BioPerl, tbl2asn) | linux-64, osx-64 |
| `envs/panaroo.yml` | `repgenr-panaroo` | Panaroo and mafft (core-gene alignment; Panaroo's dependencies bring prokka) | linux-64, osx-64 |
| `envs/ppanggolin.yml` | `repgenr-ppanggolin` | PPanGGOLiN (core-gene alignment; its dependencies bring MMseqs2 and mafft) | linux-64, osx-64, osx-arm64 |
| `envs/scarap.yml` | `repgenr-scarap` | SCARAP (core-gene alignment; its dependencies bring MMseqs2 and mafft) | linux-64, osx-64, osx-arm64 |
| `envs/orthofinder.yml` | `repgenr-orthofinder` | OrthoFinder and DIAMOND (core-gene alignment; its dependencies bring mafft) | linux-64, osx-64, osx-arm64 |

Create the core environment and only the satellites you need. The core file
installs the package in editable mode from the checkout:

```bash
mamba env create -f envs/core.yml
mamba env create -f envs/gubbins.yml      # one line per satellite you need
```

On linux-64, dRep in the core environment brings CheckM (`checkm-genome`). On
osx-arm64 the solve gives dRep without CheckM, so run dRep there with
`--ignoreGenomeQuality` or from its image. Cactus has no conda package (see
[Cactus](#cactus)). On Apple Silicon, create the osx-64 satellites with
`CONDA_SUBDIR=osx-64` and run them under Rosetta.

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

Appended satellites supply only what core lacks (`run_gubbins.py`, `iqtree2`,
`mashtree`, `snippy`, `parsnp`, `checkm2`, `progressiveMauve`). Core's samtools
and bcftools then come first, ahead of the samtools 0.1.19 in the mashtree
environment, which the `simple` SNP typer would reject. A satellite
placed before core shadows core's tools in this way. Two known gaps remain,
because a satellite tool calls its helpers by name and core's come first:

- snippy calls core's samtools and bcftools instead of the versions in its own
  environment. This combination has not been tested.
- Gubbins calls `iqtree` by name, so it runs core's IQ-TREE 3 instead of the
  IQ-TREE 2.4 of its own environment. On a small test set the result was the
  same, but the Gubbins log then records no IQ-TREE version ("for" in its place).

The per-tool `--bin-dir` option (#247) closes both gaps.

Two further points:

- On exFAT or NTFS volumes, macOS writes `._*` AppleDouble files. The adapters
  ignore them, but skDER and dRep must run on a local filesystem, so keep their
  workdir off an exFAT disk.
- progressiveMauve is not packaged for macOS. Use its image (below) or a Linux
  host.

Instead of appending, each satellite can be given to its tool alone with
`--bin-dir TOOL=DIR` (or `REPGENR_BIN_DIRS`). The directory then comes first
on `PATH` for that tool and its subprocesses only, which also closes the gap
above: snippy runs with the samtools and bcftools of its own environment.

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

The live test suite uses the same idea: its config has a `[bin_dirs]` table
whose entries become `REPGENR_BIN_DIRS` (see `tests/live/README.md`).

### 3. Containers

RepGenR2 can run each tool inside a pinned container instead of on `PATH`. The
package itself runs on the host and only the tool subprocess is containerized,
so this works in the plain CLI and inside Nextflow. The default is `none`
(native execution). The quick start and the option table are in
[usage.md](usage.md#running-tools-in-containers); the rest is here.

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
- `--platform linux/amd64` runs amd64 images on an arm64 host. Docker Desktop
  on Apple Silicon needs the Apple Virtualization framework with "Use Rosetta
  for x86/amd64 emulation" enabled. Docker's QEMU emulation cannot run some
  SIMD-heavy binaries; the bundled `vg` in Cactus hangs even on `vg version`.
  On a native amd64 Linux host no emulation is involved.
- `--wave` mints an image from the adapter's conda spec with the Seqera Wave
  CLI. It is the only image route for the four adapters that carry no pinned
  image: `simple`, parsnp, skder and SibeliaZ. Without `--wave` those four run
  on the host and the log says so.

#### Image sources

Each adapter declares its container metadata in `ToolCapabilities`:

- `container` is a pinned image URI, used as is when Wave is off. The
  single-package adapters pin a BioContainer (galah, sourmash, dRep, snippy,
  ska2, Gubbins, IQ-TREE, FastTree, RAxML-NG, mashtree, and the assemblers and
  polishers). progressiveMauve pins a BioContainer with a compatible boost, and
  Cactus pins its project image. Four adapters have no pin and run in an image
  only under `--wave`: the `simple` SNP typer (minimap2, samtools, bcftools)
  and parsnp (parsnp, harvesttools) span several packages, and the skder and
  SibeliaZ BioContainers are BusyBox-based. The GNU-only calls in their shell
  wrappers (`sort --parallel` in skder's greedy mode, `mktemp --suffix` in
  SibeliaZ) fail there and the run silently yields nothing.
- `conda` is a conda spec (for example `bioconda::skder`). With `--wave`,
  RepGenR2 mints an image for it. This is arm64-native and the only route for
  the multi-package adapters. The pin stays the default without Wave.

BioContainers are `linux/amd64`. On Apple Silicon pass
`--platform linux/amd64` with Rosetta, or use `--wave`.

#### Storage location

- **Singularity or Apptainer.** `--container-cache` sets `APPTAINER_CACHEDIR`
  and `SINGULARITY_CACHEDIR` (and the matching `*_TMPDIR`). `docker://` images
  are pulled once to `<cache>/<name>.sif` and reused. Put this on a large or
  external disk.
- **Docker.** The daemon manages image storage (Docker Desktop's disk image).
  It is set there, not per run.
- Large run-time data (the Cactus job store, scratch, downloads) lives under
  `--workdir` and `TMPDIR`.

#### Notes

- Before a tool runs in an image, RepGenR2 checks that the engine binary is
  on `PATH` and, for Docker, that `docker info` reaches the daemon. A stopped
  Docker Desktop therefore exits 4 before the stage starts, and
  `repgenr --container docker list-tools --check` reports it per tool. Under
  a container backend each line of `list-tools --check` names where the tool
  runs: `[image <ref>]` for a tool with a pinned (or `--wave`) image and
  `[host]` for one that runs from `PATH`, so a host version is not mistaken
  for the image's. The image itself is pulled on first use; `--images` adds,
  for each tool that passed the check, whether its images are present locally
  (`present` or `not pulled`, from `docker image inspect` or the Singularity
  `.sif` cache; nothing is pulled). An adapter's secondary images, such as
  racon's minimap2 image, are listed too. `presence unknown` means the engine
  did not say the image is absent, for example because the daemon is down.
- A stage that ran a tool in an image records the image reference as that
  tool's version and the engine with its version (`docker: 29.5.3`) beside it,
  in `repgenr.yaml` and in `repgenr versions`.
- Docker runs as the host UID and GID so outputs are owned by you. The workdir
  and `TMPDIR` are bind-mounted at identical paths.
- Symlinked inputs (a `genomes/` directory staged by `repgenr ingest`) are
  followed. The directory each link points to is bound as well, so the tool
  sees the same paths inside the container.
- On Apple Silicon most BioContainers are `linux/amd64`. They run under Docker
  emulation, or use `--wave` for arm64-native images.
- The macOS SibeliaZ BSD-wrapper workaround is skipped automatically when
  running in a (Linux) container.
- dRep's CheckM step needs its reference data at run time. Set
  `CHECKM_DATA_PATH` on the host: RepGenR2 passes it into the container and
  binds that directory at the same path. Otherwise run dRep with
  `--ignoreGenomeQuality`.
- progressiveMauve from bioconda (`mauvealigner`) is linked against boost
  1.74 and fails with an `undefined symbol` error against a newer boost.
  `envs/mauve.yml` and the adapter's Wave spec therefore pin `boost-cpp=1.74.0`;
  the pinned image already holds a compatible boost.

## Per-tool table

"Native on macOS arm64" records what was run on the audit machine (Apple
Silicon) and is taken from the adapter status table in
[verification.md](verification.md#adapter-status). "No" means not run natively
there, not that the tool cannot work on another platform. The container pin is
the image in the adapter's declared capabilities; `none` means the adapter
carries no pin and needs `--wave` to run in an image. Minimum versions are in
the `envs/` files, and the environment that holds each tool is listed in
[choosing-tools.md](choosing-tools.md).

| Family | Tool | conda package | Native on macOS arm64 | Container pin | Databases or extras |
|---|---|---|---|---|---|
| Dereplicator | drep | `drep` (brings `checkm-genome` on linux-64) | no | `drep:3.7.1` | CheckM reference data through `CHECKM_DATA_PATH`, or `--ignoreGenomeQuality` |
| Dereplicator | skder | `skder`, `skani` | yes | none (Wave) | none |
| Dereplicator | galah | `galah` | yes | `galah:0.4.2` | none |
| Dereplicator | sourmash | `sourmash` | yes | `sourmash:4.9.4` | `sourmash_plugin_branchwater` for the sparse back-end (extra `sparse`) |
| Aligner | progressivemauve | `mauvealigner`, `boost-cpp=1.74.0` | no (unpackaged on macOS) | `mauve:2.4.0.snapshot_2015_02_13` | none |
| Aligner | sibeliaz | `sibeliaz` | yes | none (Wave) | none |
| Aligner | cactus | not on conda | no | `cactus:v2.9.3` | none |
| SNP typer | simple | `minimap2`, `samtools`, `bcftools` | yes | none (Wave) | samtools and bcftools 1.10 or later |
| SNP typer | snippy | `snippy` | no | `snippy:4.6.0` | none |
| SNP typer | parsnp | `parsnp`, `harvesttools` | yes, from an osx-64 environment | none (Wave) | none |
| SNP typer | ska2 | `ska2` | yes | `ska2:0.5.1` | none |
| Masker | gubbins | `gubbins` | yes | `gubbins:3.4.3` | a multi-threaded RAxML build for more than one thread (see [usage.md](usage.md#snp-typing-and-masking)) |
| Tree builder | iqtree | `iqtree` | yes | `iqtree:3.1.3` | none |
| Tree builder | fasttree | `fasttree` | yes | `fasttree:2.2.0` | none |
| Tree builder | veryfasttree | `veryfasttree` | yes | `veryfasttree:4.0.5` | none |
| Tree builder | raxmlng | `raxml-ng` | yes | `raxml-ng:2.0.3` | none |
| Tree builder | mashtree | `mashtree` | yes | `mashtree:1.4.6` | none |
| Tree builder | sourmash | `sourmash` | yes | `sourmash:4.9.4` | `rapidnj` for the joining step; without it more than 5000 leaves are refused (optional; `rapidnj:2.3.3` under a container backend) |
| Tree builder | skani | `skani` | yes | `skani:0.3.2` | `rapidnj` for the joining step above a few thousand leaves (optional; `rapidnj:2.3.3` under a container backend) |
| Tree builder | ska2 | `ska2` | yes | `ska2:0.5.1` | `rapidnj` for the joining step above a few thousand leaves (optional) |
| Assembler | skesa | `skesa` | yes | `skesa:2.5.1` | none |
| Assembler | shovill | `shovill` | no | `shovill:1.4.2` | none |
| Assembler | flye | `flye` | yes | `flye:2.9.6` | none |
| Classifier | sourmash | `sourmash` | yes | `sourmash:4.9.4` | a GTDB sourmash sketch and its lineages CSV |
| Polisher | medaka | `medaka` | no | `medaka:2.2.2` | medaka models. Reads from SRA have no basecaller in their headers; the adapter then assumes ONT's bacterial R10.4.1 model and records it. Give another with `--tool-arg model=...` |
| Polisher | racon | `racon`, `minimap2` | yes | `racon:1.5.0`; minimap2 runs from `minimap2:2.28` | none |
| Annotator | bakta | `bakta` | yes | `bakta:1.12.2` | bakta database: `bakta_db download --output DIR --type light` (1.4 GB download, about 4 GB unpacked; `--type full` is 30 GB and about 84 GB unpacked, with more specific CDS annotation), then `--bakta-db DIR` or `BAKTA_DB` |
| Annotator | prokka | `prokka` | no (tbl2asn is not packaged for arm64) | `prokka:1.15.6` | none |
| Core-gene builder | panaroo | `panaroo`, `mafft` | no (its dependency prokka is not packaged for arm64; `CONDA_SUBDIR=osx-64` works) | `panaroo:1.8.0` | none |
| Core-gene builder | ppanggolin | `ppanggolin`, `mafft` | yes | `ppanggolin:2.3.2` | none |
| Core-gene builder | scarap | `scarap`, `mmseqs`, `mafft` | yes | `scarap:1.0.2` | none |
| Core-gene builder | orthofinder | `orthofinder`, `diamond`, `mafft` | yes | `orthofinder:3.1.5` | none |
| Quality (reads chain) | checkm2 | `checkm2` | no | `checkm2:1.1.0` | CheckM2 database: `checkm2 database --download`, then `--checkm2-db` or `CHECKM2DB` |

The genome download steps also need `ncbi-datasets-cli`.

## Cactus

Cactus has no conda package in `envs/` and no native install path
in these docs. The adapter pins the project image
(`quay.io/comparative-genomics-toolkit/cactus:v2.9.3`), so use
`--container docker|singularity`. On Apple Silicon this needs the Rosetta
setting described above.

## Nextflow

The Nextflow layer needs Nextflow 26.04 or later (`nextflowVersion =
'!>=26.04.0'` in `nextflow/nextflow.config`). The nf-test suite uses nf-test
0.9.5, the version CI installs. The `docker`, `singularity` and `wave`
profiles set `params.repgenr_opts`, so every stage calls `repgenr` with
`--container docker`, `--container singularity` or `--container singularity
--wave`. The `repgenr` command itself must be available to each process. See
[usage.md](usage.md#nextflow) for parameters and profiles.

## Verify the installation

```bash
repgenr --version
repgenr list-tools --check
```

`list-tools` prints each family with its adapters and, where one is declared,
the genome limit. Its last line names the dereplicators that `glance` can run
(`drep`, `sourmash`). `--check` runs each adapter's preflight and prints `ok` with
the versions found, or `missing` or `error` with the reason, or `broken` for a
plugin that failed to load. It reports and exits with status 0. Add `--strict`
to verify an environment from a script: after the full listing it exits 4 when
any adapter is missing or errored, and 5 when any plugin is broken. A stage run with a missing or outdated tool exits with
status 4, except under `assemble --assembler auto`, which excuses runs whose assembler is missing (reason `assembler_not_installed`, with a warning) and exits 4 only when no run can be assembled. A missing `auto` polisher is skipped (see [usage.md](usage.md#starting-from-sequencing-reads)). The full exit-code table is in
[usage.md](usage.md#troubleshooting).
