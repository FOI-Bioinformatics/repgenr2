# Usage

This page is the user guide for RepGenR2 (the command is `repgenr`). It
follows the order of an analysis: building a genome set, dereplicating it,
inferring a phylogeny, and deriving a FlexTaxD taxonomy from the tree, then
running the same work under Nextflow or with the tools in containers. It is
written for readers who know genomics but not this tool.
[install.md](install.md) explains how to install the package and the external
tools, and [choosing-tools.md](choosing-tools.md) explains which dereplicator,
phylogeny route and tree builder suit a dataset. Every command and option is
listed in [cli-reference.md](cli-reference.md); the files each stage writes are
described in [output.md](output.md).

## Contents

1. [Introduction](#introduction): stages, the working directory, threads,
   running the whole chain, resume and `--force`, checking a working directory
2. [Building the genome set](#building-the-genome-set): bacteria from GTDB,
   local genomes, sequencing reads, viruses, genome sketches, counting taxa
3. [Dereplication and representative selection](#dereplication-and-representative-selection):
   representative selection, chunked dereplication, inspecting a dereplication
4. [Phylogeny](#phylogeny): choosing a route, SNP typing and masking,
   annotation, core-gene phylogeny, rooting, reusing an alignment
5. [tree2tax](#tree2tax): collapsing weak splits, node names, a rank column
6. [Nextflow](#nextflow): pipeline structure, parameters, profiles, scaling,
   configuring processes
7. [Running tools in containers](#running-tools-in-containers)
8. [Troubleshooting](#troubleshooting), including the exit codes

## Introduction

### Stages and entry points

RepGenR2 builds a representative-genome repository for a taxonomic group and
the matching FlexTaxD taxonomy. A representative genome is the genome kept for
a cluster of closely related genomes. The pipeline runs five stages in order:

1. **metadata** -- select genome accessions (GTDB for bacteria/archaea; the
   viral path uses `vmetadata`/`vgenome`, sourcing from NCBI Virus by default).
2. **genome** -- download and organise the selected genomes with NCBI Datasets.
3. **dereplicate** -- cluster genomes by ANI (average nucleotide identity) and
   pick representatives.
4. **phylo** -- build a phylogeny from the representatives.
5. **tree2tax** -- emit a FlexTaxD-compatible taxonomy from the tree.

The first two stages can be replaced by another entry path:

| Entry path | Commands | Section |
|------------|----------|---------|
| Bacteria and archaea from GTDB | `metadata -> genome` | [Bacteria from GTDB](#bacteria-from-gtdb) |
| Genomes already on disk | `ingest` | [Starting from local genomes](#starting-from-local-genomes) |
| Sequencing runs from ENA/SRA | `reads -> assemble` | [Starting from sequencing reads](#starting-from-sequencing-reads) |
| Viruses from NCBI Virus or BV-BRC | `vmetadata -> vgenome` | [Viruses](#viruses) |

Each entry path writes the same genome set (`genomes/`, `selection.tsv` and
the manifest), so `dereplicate`, `phylo` and `tree2tax` are the same for all
of them. An entry path can also stage an outgroup, a genome from outside the
target taxon that is used to root the tree.

The stages above write into a working directory. The Nextflow layer runs the
same work as stateless data channels without a shared working directory (see
[Nextflow](#nextflow)).

### The working directory

Every stage command takes the working directory with `-wd` (long form
`--workdir`). Stages read and write inside it, and `repgenr.yaml` there records
what ran. Large run-time data (scratch, downloads, the Cactus job store) also
lives under the working directory and `TMPDIR`, so put it on a disk with room.

### Threads

Commands that run a threaded tool take `-t/--threads`. Without it the count is
16, or the CPU limit of the repgenr process when that is lower. The CPU limit
is the set of CPUs the process may run on (a Slurm cpuset, `taskset`) and a
cgroup CPU quota (`docker run --cpus`, a Kubernetes CPU limit), rounded up. A
lowered default is reported once on the console. An explicit `-t` is used as
given, with a warning when it exceeds that limit. The thread count is not a
result parameter, so changing it does not make a finished stage rerun.

### Running the whole chain

`repgenr run` runs the whole chain in one command (bacterial by default;
`--viral` for the NCBI Virus path). `repgenr status` shows progress at any
time:

```bash
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` forwards the stage options it shares a name with; for the rest, run
the stage command by hand (`repgenr status` says which comes next). Two
options are worth knowing:

- `--with-snptype` adds the standalone `snptype` stage after dereplication,
  so the SNP tables under `snp/` are produced even when the tree is built
  another way.
- `--genomes-dir` starts the chain from local genomes (see
  [Starting from local genomes](#starting-from-local-genomes)).

With `--msa-source snptype`, `phylo` runs its own typing pass, outgroup
included, into `tree/msa/`. `snp/` holds only the `snptype` stage's tables, so
the two do not replace each other, and a later `phylo` run that changes only
the tree builder reuses its alignment.

### Resume and `--force`

Each stage records its parameters, the digests of its inputs, and the container
identity in `repgenr.yaml`. Re-running a stage is a safe no-op only when all
three are unchanged (it logs that it skipped). Re-running an upstream stage
(e.g. `repgenr --force dereplicate ...`) changes a downstream stage's input
digests, so the downstream stage re-runs automatically the next time it is
invoked. Change a parameter, switch `--container`, or pass `--force` to re-run
explicitly.

`--force` is a global option and goes before the command name
(`repgenr --force dereplicate -wd WD`); after the command name it is rejected
as an unknown option.

**Interrupted and refused stages.** A stage writes its record without a
completion stamp before it starts, so one that failed or crashed mid-run is
listed as `[interrupted]` by `status`, reported as a failure by `doctor`, and
always re-runs; a successful run stamps the record. A failure in parameter
validation writes no record. A stage that refuses (exit 2 or 3, or 4 for a
tool missing at its preflight) without changing any of its main outputs
leaves the record as it was: none on a first run, and the last finished one
on a re-run, since its outputs are untouched. Examples are `phylo` with fewer
than three genomes and `snptype --mask gubbins` without Gubbins. A failed
external tool (exit 6), or a refusal after a main output changed, leaves the
record interrupted.

**Missing outputs.** Before skipping, a stage also checks that its main
outputs exist, for example:

- `genomes/` and `manifest.sqlite` for `ingest`;
- `derep/clusters.tsv` and `derep/representatives/` for `dereplicate`;
- `tree/tree.nwk` for `phylo`;
- `tree2tax.tsv` and `genomes_map.tsv` for `tree2tax`.

A directory must not be empty, and dotfiles such as `.DS_Store` or exFAT `._`
files do not count. `genomes/` and `derep/representatives/` count only with a
genome FASTA file in them, so one holding only a leftover `x.fasta.tmp` is
treated as empty. Each genome listed in `selection.tsv` is checked for the
stage that wrote the genome set (`genome`, `ingest`, `vgenome`, `assemble`),
and each representative listed in `derep/clusters.tsv` for `dereplicate`. If
one was deleted, the stage logs
`Stage 'X': deliverable <path> missing; re-running.` and runs again, so
`--force` is not needed to rebuild it; a genome deleted from `genomes/` is
restored by the next run of that stage. To drop a genome on purpose, remove
it from the source directory or from the `--selection` file and re-run the
entry stage. For `glance` only `glance_clustering_dendrogram.pdf` is checked,
since its plots are absent when no value falls within the plot bounds.

**Input digests.** Input directories (`genomes/`, `derep/representatives/`,
`outgroup/`, the `ingest --genomes-dir` source) are digested from the genome
FASTA files in them, the files the stages read (names ending in `.fasta`,
`.fa`, `.fna`, `.fas`, `.fasta.gz`, `.fna.gz` or `.fa.gz`). A leftover
`x.fasta.tmp`, a `.fai` index, a `truth.json` or a README beside the genomes
therefore does not re-run the stages that read the directory. A working
directory recorded by an earlier version re-runs a stage once only when such
a file was present (for example `ingest` from a source directory that holds a
`truth.json`); `derep-stock --action unpack` repeats once, since its stored
run is now digested file by file.

**Limitations.** Two limitations are both covered by `--force`:

- Input directories are digested from file metadata (name, size, mtime), so
  an in-place edit that preserves size and mtime is not detected.
- Upgrading a natively installed tool binary does not invalidate previous
  results (switching the container backend or platform does).

Working directories created by earlier versions re-run each stage once (the
fingerprint format changed).

### Checking a working directory

`repgenr status -wd <WD>` shows what completed and what is next.
`repgenr doctor -wd <WD>` verifies a working directory's outputs against its
records (missing, corrupt or untracked genomes, manifest drift, truncated
deliverables, interrupted stages).

`doctor` exits 7 on failures and 0 when it finds only warnings; a missing
deliverable or a changed input is a warning under the same path, since the
stage will re-run. `doctor --quick` skips reading the first bytes of each
genome file, the slowest check on a large genome set (about 30 s per 1000
files on an exFAT disk). Links, missing and untracked genomes are still
checked, and the genome line says the content was not read, so a non-FASTA
file saved under a FASTA name goes unnoticed.

`status` reads the same tables as `doctor`, so the two commands agree on what
re-runs:

- A stage that failed is listed as `[interrupted]`, and its outputs may be
  partial until it is re-run.
- A finished stage is listed as `[stale]` when one of its inputs changed or
  one of its outputs is missing since it finished; it re-runs on its next
  invocation.
- A finished stage whose record holds no resume fingerprint (written by an
  older version, or restored by `derep-stock --action unpack`) stays `[done]`
  with a note that its next invocation recomputes it. `doctor` warns about
  every completed record without a fingerprint, stale ones included; `status`
  adds the note only to `[done]` lines, since a `[stale]` line already says
  the stage re-runs.
- When `repgenr.yaml` records no stage, `status` names the entry stages (and
  suggests `doctor` when the working directory holds outputs). When it
  records stages but no entry stage, `status` follows the stages all
  pipelines share (dereplicate, phylo, tree2tax) under "Pipeline: unrecorded
  entry stage".

A malformed `repgenr.yaml` is reported by both: `doctor` as a failure,
`status` and every stage with exit 3.

Both commands accept `--json` and then print one versioned JSON object on
stdout instead of the text report, for scripts and workflow wrappers. The
exit codes are unchanged, and on exit 3 stdout stays empty. The schemas are
in [output.md](output.md#machine-readable-status).

## Building the genome set

This part describes the entry paths that produce a genome set in a working
directory, the sourmash sketches written beside it, and `census`, which counts
what a taxon or a working directory holds.

### Bacteria from GTDB

This route selects bacterial or archaeal genomes from GTDB (the Genome
Taxonomy Database) with `metadata` and downloads them from NCBI with
`genome`. Use it when the target is a prokaryotic taxon in GTDB. `metadata`
reads either the full GTDB metadata table (`--source tsv`, the default) or
the GTDB API (`--source api`).

```bash
WD=./francisella

# Full GTDB table (release-pinned):
repgenr metadata -wd $WD -r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella
# Or query only the target taxon via the GTDB API (no full-table download):
repgenr metadata -wd $WD --source api -d rep -l genus -tg francisella

repgenr genome -wd $WD
repgenr dereplicate -wd $WD --tool skder -t 16
repgenr phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
repgenr tree2tax -wd $WD --include-dereplicated
```

#### The GTDB metadata table

With `--source tsv` (the default), `-r/--release` (major.minor, e.g. `232.0`)
and `--gtdb-version` (`bac120` or `ar53`) are required, because they select
the table's release and domain and are recorded in the provenance. They stay
required with `--metadata-path <table>`, which reads a GTDB metadata table you
already have instead of downloading one; a path that does not exist exits 2.

The downloaded table stays in the working directory
(`bac120_metadata_r232.tsv.gz`), and `--nodownload` reuses it instead of
downloading it again. Table names carry only the major release, so a download
also writes `<version>_metadata_r<major>.release` with the exact release, and
`--nodownload` refuses a table fetched for another minor release (exit 2). A
table without that file, for instance one placed by hand, is reused with a
warning.

The download is checked against its size and against GTDB's `MD5SUM.txt`; a
checksum mismatch removes the file and exits 3. A release or version that GTDB
does not publish exits 2 after trying the current (`.tsv.gz`) and the legacy
(`.tar.gz`) layout. A network failure exits 3 with its cause and does not try
the second layout.

#### The GTDB API

`--source api` serves GTDB's current release and ignores `-r`,
`--gtdb-version`, `--metadata-path` and `--nodownload` (a warning names them).
The API does not report a release number, so the stage record holds
`release: null` and, in its place, `api_query_date`, the UTC time of the query
(ISO 8601). `status` shows it after the metadata line and `versions` prints it
as `gtdb_api_query_date`; for the table path they show the release
(`gtdb_release`). To relate a query date to a GTDB release, compare it with
the release dates on the GTDB website.

#### Taxon names and the outgroup

GTDB taxon names are case-sensitive. The target's first letter is raised and
a species epithet is lowered, so `-tg francisella -ts Tularensis` finds
`s__Francisella tularensis`. Suffixes such as `Bacillus_A` or `copri_A` must
be typed as GTDB spells them. A taxon or `--outgroup-accession` the API does
not know exits 2.

A named outgroup (`--outgroup-accession`) must lie outside the target taxon
on both sources, like the automatic one (exit 2 otherwise). A target genome
that `-d rep` or `--limit` left out of the selection is still not an
outgroup.

#### Limiting the selection

`repgenr metadata --limit N` caps the bacterial selection at N genomes. The
cap is not the first N rows of the GTDB table. Candidates are grouped by
species and taken round-robin: the best-quality genome of every species
first, then each species' next best, until N. A heavily sequenced species
therefore cannot fill the cap on its own. With `-d rep` there is one genome
per species and the rule reduces to a quality ranking across species.

Within a species, genomes rank by the CheckM part of the keeper score
(completeness minus five times contamination; see
[Representative selection](#representative-selection)). The ranking happens
before download, so no N50 exists yet. Unscored genomes come last, and ties
go to the GTDB species-representative flag, then the accession, so the result
is deterministic.

On the `--source api` path the per-genome quality cards are fetched for every
candidate before the cut, one request per genome, four at a time (about five
genomes a second; the 1540 Wolbachia genomes take some five minutes). The API
refuses sustained bursts now and then; refused cards are retried once more,
slowly, and only a genome refused twice is left unscored.

The outgroup is chosen afterwards from the parent taxon, outside the target
taxon (a target genome the cap left out is not a candidate), and never counts
against the limit.

### Starting from local genomes

`repgenr ingest` populates a working directory from local FASTA files instead
of downloading. Use it for genomes you already have, or to combine the genome
sets of earlier working directories. It stands in for `metadata` and
`genome`, so the local chain is `ingest -> dereplicate -> phylo -> tree2tax`.

```bash
repgenr ingest -wd $WD --genomes-dir ./my_genomes --outgroup GCF_003574425.1
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
```

`repgenr run --genomes-dir ./my_genomes ...` runs the same local chain in
one command; `--from-workdir`, `--selection`, `--outgroup` and `--copy` pass
through to `ingest`, and the GTDB selection flags are not needed.

`ingest -wd WD --genomes-dir DIR` links the files (or copies them with
`--copy`) into `genomes/` and writes `selection.tsv` and the SQLite manifest in
the same format `metadata`/`genome` produce. The stage is recorded, so
`status` reports the local chain and a second `ingest` on an unchanged
directory skips. Taxonomy and quality columns come from
`--selection selection.tsv` when given, otherwise from the canonical
`Family_genus_species_ACCESSION.fasta` filename. This is the CLI counterpart
of the Nextflow harness `nextflow/tests/local_dataflow.nf`, which takes a
genome directory straight into `DEREPLICATE_SCATTER`.

#### File names, accessions and the outgroup

Without `--selection`, accession and taxonomy come from the filename:

- `Family_genus_species_ACCESSION.fasta` (four or more `_`-separated tokens)
  gives the first three tokens as taxonomy and the rest as the accession.
  Any name with four or more tokens is read this way, so
  `sample_1_run_A.fasta` gives the accession `A`.
- A name that starts with an NCBI assembly accession
  (`GCF_000008985.1_ASM898v1_genomic.fna`, as NCBI Datasets and the FTP site
  deliver) gives that accession and no taxonomy.
- Any other name gives no taxonomy and the file stem as the accession.

Two genomes with one accession (for example `x.fasta` and `x.fna`, or two
four-token names that end alike) stop `ingest` with exit 2; give them
distinct accessions with `--selection`. The `is_outgroup` column of a
selection takes `1`/`0` (also `true`/`false`, `yes`/`no`), and so does the
optional `gtdb_representative` column, which `--keeper gtdb` reads; without
that column no genome is a GTDB representative.

`--outgroup` sets one genome aside under `outgroup/` and writes
`outgroup_accession.txt`, so `phylo` roots on it exactly as after a download.
It names a genome under `--genomes-dir` or in a `--from-workdir` (filename,
stem or accession), or a FASTA file anywhere. The genome is kept out of the
ingroup. When `--selection` also marks an outgroup row, both must name the
same genome; otherwise `ingest` exits 2 and names both. A file from outside
`--genomes-dir` whose name gives the filename or accession of an ingroup
genome is refused (exit 2), since it would replace that genome.

#### Accepted files and compressed genomes

Only the files directly under `--genomes-dir` with a suffix `.fasta`, `.fa`,
`.fna`, `.fas`, `.fasta.gz`, `.fna.gz` or `.fa.gz` are read, so the `.fna.gz`
files of the NCBI FTP site can be ingested as they are. Other files
(`X.FASTA`, `x.fas.gz`) are listed in a warning and skipped. Subdirectories
are not searched (an NCBI Datasets download keeps each genome in its own
directory, so collect the `.fna` files into one directory first). An empty or
unreadable genome file (a dangling link included) stops `ingest` with exit 2
before anything is staged, and the record of an earlier `ingest` is left as
it was; `doctor` checks that the staged files hold FASTA.

Compressed genomes are staged unchanged, and every tool accepts them. The
dereplicators (sourmash, skDER, galah; dRep through a decompressed copy),
cactus, the `mashtree` and `sourmash` tree builders and the `ska2` and
`simple` typers read gzip themselves. progressiveMauve, SibeliaZ, `parsnp` and
`snippy` cannot (they crashed or wrote an empty alignment before this was
handled), so `snptype` and `phylo` give them a decompressed copy of each
gzipped genome in the stage's scratch directory (`scratch/snptype/inputs/`,
`scratch/phylo_snptype/inputs/` or `scratch/phylo_inputs/`) and remove the
copies when the tool has finished. The copies keep the genome's record name,
so alignments and trees are the same as for uncompressed files.

The scratch directory therefore needs room for the uncompressed genomes while
the tool runs. The stage checks the free space first (about four times the
gzip size) and stops when less than 1 GB is free. A truncated or corrupt gzip
genome stops the stage with exit 2 naming the file.

#### Links, copies and changed sources

By default `genomes/` holds symbolic links to the source files, which takes
under a second for 1000 genomes. If the source is later moved or deleted the
links dangle: `doctor` reports them, and `ingest` from the new location
restages the set. `--copy` makes the working directory independent of the
source (2 GB, 1000 genomes, from an exFAT disk took about 2.5 minutes); the
copies keep the source's permissions, so read-only sources give read-only
copies.

A changed source (files added, removed or rewritten) re-runs `ingest` on its
next invocation, prunes genomes no longer present, and marks `dereplicate`
and later stages for a re-run. An unchanged set keeps the dereplication
status recorded in the manifest.

#### Combining working directories

`--from-workdir WD` (repeatable, alone or together with `--genomes-dir`)
takes the genome set of an earlier working directory: the rows of its
`selection.tsv` and the files they name under its `genomes/`. Use it to
dereplicate together genomes selected from GTDB in one working directory and
genomes assembled from sequencing runs in another:

```bash
# GTDB genomes
repgenr metadata -wd WD1 --source api -d all -l species -tg francisella -ts tularensis
repgenr genome -wd WD1
# Assemblies from sequencing runs
repgenr reads -wd WD2 -ts "Francisella tularensis" --platform illumina
repgenr assemble -wd WD2 --checkm2-db /db/checkm2/uniref100.KO.1.dmnd \
    --gtdb-sketch /db/gtdb/gtdb-rs226-reps.k31.sig.zip
# Both sets in a third working directory
repgenr ingest -wd WD3 --from-workdir WD1 --from-workdir WD2 --outgroup GCF_003574425.1
repgenr dereplicate -wd WD3 --tool skder --keeper gtdb
```

`assemble --append` is the alternative within one working directory: it adds
the assemblies to the GTDB selection already there (see
[Adding assemblies to a GTDB selection](#adding-assemblies-to-a-gtdb-selection)).

**What carries over.** From each source row, `ingest` takes the accession,
taxonomy, filename, CheckM completeness and contamination, and
`gtdb_representative`, so `--keeper quality` and `--keeper gtdb` work on the
combined set. The manifest keeps the source of each genome (`gtdb`, `sra`,
`local` and so on) when the source working directory has a manifest, and
records `local` when it does not; genomes from `--genomes-dir` are `local`.
`--selection` applies to `--genomes-dir` only. Per-run files of the reads
working directory (`assembly_stats.tsv`, `excused_runs.tsv`, `reads.tsv`)
stay there and are not copied.

**The outgroup.** The outgroup of a source working directory is not carried
over; `ingest` logs one line for each source that had one. `--outgroup` sets
the outgroup of the new working directory and may name a genome of any source
(filename, stem or accession, under `genomes/` or `outgroup/` of a source
working directory, or under `--genomes-dir`), or a FASTA file elsewhere.

**Conflicts and resume.** A genome found in two sources (the same accession,
or the same filename with different accessions) stops `ingest` with exit 2
and names both sources; neither takes precedence. A selection row whose file
is missing, or a `--from-workdir` that is the target working directory, also
exits 2. The selection and the `genomes/` directory of each source are resume
inputs, so a change in a source working directory reruns `ingest`. The source
manifest is not, so a change of a source label alone needs `--force`.

**State of the sources.** Before staging, `ingest` reads the record of the
stage that wrote each source's genome set (the latest of `genome`, `vgenome`,
`ingest` and `assemble`) and evaluates it as `repgenr status` would in that
working directory. A source whose record is interrupted, or stale against its
own inputs (an input changed, or a deliverable is missing), is named in one
warning with the stage and its state, and its genomes are taken as they are.
A source with no such record (a `selection.tsv` and `genomes/` written by
other means) is accepted with the warning "no stage record".
`--strict-sources` turns each of these warnings into a refusal (exit 2)
before anything is written. The state of each source is kept in the ingest
record (`source_states`), and `status` names the sources that were not done
when they were ingested.

### Starting from sequencing reads

This route builds the genome set from whole-genome sequencing runs. Use it for
a taxon whose public isolates are available as reads but not as assemblies,
or to add such isolates to a GTDB selection. `repgenr reads` selects runs from
ENA (which mirrors SRA) and writes `reads.tsv`; `repgenr assemble` fetches and
assembles them into `genomes/` with the same `selection.tsv` and manifest the
other entry paths write. The chain is
`reads -> assemble -> dereplicate -> phylo -> tree2tax`.

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina --max-runs 20
repgenr reads -wd $WD --accession PRJNA954307 --accession SRR28800588
repgenr assemble -wd $WD --assembler auto -t 16 --jobs 2
# or the whole chain:
repgenr run -wd $WD --reads -ts "Francisella tularensis" --platform illumina --max-runs 20 \
    --assembler skesa --treebuilder mashtree
```

`run --reads` forwards `--accession-file`, `--platform`, `--max-runs`,
`--assembler`, `--threads` and `--outgroup`; the rest is available on the
stage commands. Which assembler and polisher suit each sequencing platform is
in [choosing-tools.md](choosing-tools.md#assemblers-and-polishers).

#### Selecting runs

Runs are chosen by taxon or by accession.

- **By taxon:** `--target-family`/`-tf`, `--target-genus`/`-tg` or
  `--target-species`/`-ts`, resolved through the ENA taxonomy, synonyms
  included. The number of runs under the taxon is logged before the search,
  with a warning above 100000 runs, such as Salmonella with about 880000,
  whose records take about 3 GB of memory.
- **By accession:** `--accession` takes a run (SRR/ERR/DRR), a sample
  (SAMN.., SRS..) or a study (PRJNA.., SRP..) and repeats; `--accession-file`
  lists them one per line (text from a `#` to the end of the line is a
  comment). Runs found by accession pass the same filter as the taxon query,
  whole-genome sequencing of genomic DNA: RNA-Seq, amplicon or metagenomic
  runs of a named study or sample are dropped with a warning that names them.

The selection is then filtered:

| Option | Effect |
|--------|--------|
| `--platform illumina\|ont\|pacbio` | Keeps one platform. |
| `--min-bases` | Drops small runs. |
| `--max-bases` | Drops runs above a size. An unenriched whole-host library, tens of Gb for a 1 Mb endosymbiont, would assemble into a host-dominated genome. |
| `--drop-selection` | Drops runs by ENA library selection. The default is `MDA` (whole-genome amplification), which assembles into chimeric and uneven contigs. Repeat the flag for more values; `--drop-selection none` keeps every run. The value is kept in `reads.tsv` as `library_selection`. |
| `--one-per-sample` (the default) | Keeps the best run of each sample: a long-read run when it carries at least 100 Mb and a tenth of the sample's largest short-read run, else the largest run. `--all-runs` keeps every run. |
| `--max-runs` | Caps the selection to the largest runs. |

Each run is labelled with the family, genus and species of its NCBI taxid, in
the same filename tokens the GTDB path uses.

#### Sketching runs without downloading them

`reads-sketch` sketches the reads of the selected runs without a download.
Each run is streamed from NCBI with `fastq-dump -Z` and piped into one
`sourmash sketch dna` call, so no FASTQ file is written, and `--max-spots N`
stops the stream after the first N spots of a run (a paired-end spot gives
two reads). The sketches have the parameters of the reads sketch `assemble`
writes (k=21, 31 and 51 at scaled=1000, with abundances), so the same
consumers read them: a `sourmash gather` against the GTDB sketch to see what
a run holds before assembling it, or a comparison between runs.

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina
repgenr reads-sketch -wd $WD --max-spots 200000            # every run of reads.tsv
repgenr reads-sketch -wd $WD --accession SRR1 --accession SRR2   # a subset
repgenr reads-sketch -wd $WD --accession-file picked.txt --platform illumina
repgenr reads-sketch -o sketches/ --accession SRR25474756 --max-spots 100000
```

The runs come from the workdir's `reads.tsv`; without a selection every run
is sketched. `--accession` and `--accession-file` restrict to runs of that
table (an accession not in it is refused by name, so a typo cannot sketch an
unrelated run), and `--platform` to one platform. A call adds to the set: a
run sketched before with the same `--max-spots` is skipped and a later call
with another selection keeps it, so a few runs can be sketched first and the
rest later. Each cap is its own variant with its own directory,
`reads_sketches/all/` for the whole run and `reads_sketches/spots_<N>/` for
`--max-spots N`, so sketches at several depths live side by side and a new
cap fills a new directory without touching the others; `--force` streams the
selected runs of the current variant again. The sketches go to
`reads_sketches/<variant>/<run>.sig.zip` and the index
`reads_sketches/<variant>/reads_sketches.tsv` lists every run of `reads.tsv`
with its status (`ok`, `failed` or `pending`), the cap, the spots and
sequences the tools reported and the reason of a failure; see
[output.md](output.md#reads-sketch). The stage is recorded as `reads_sketch`
(the record names the variant of the last call and the runs each variant
holds), and `status` lists the variants present. A failed run (an accession
NCBI does not hold, a network error) is recorded and the others continue;
since a selected run without a sketch is a missing deliverable, a repeat of
the same call streams the failed runs again and skips the rest; a run that
keeps failing, such as an accession NCBI does not hold, is retried on every
repeat until the selection is narrowed with `--accession`. The command
fails (exit 3, as `assemble` with no accepted run) only when every selected
run failed, and then records nothing.

`--jobs` streams that many runs at once (default 4, or the usable CPUs when
fewer); each is a `fastq-dump` and a `sourmash` process, bound by the
network more than by the CPU. The standalone form `-o DIR` with bare
`--accession` runs (SRR, ERR or DRR accessions), or with `--reads-tsv` the
runs of a `reads.tsv`, sketches into `DIR` without a workdir and records
nothing; `--platform` needs `-wd` and is refused there, and `--versions-out`
writes the tool versions as a YAML fragment (the Nextflow `READS_SKETCH`
task uses this form).

The first N spots of a run are not a random sample: an Illumina file starts
with its first tiles, whose quality profile differs, and an ONT file with
the first reads of the run. For a species-level `gather` on an isolate this
is expected to matter little; a low-abundance organism in a mixed sample may
be missed in a short prefix. A few hundred thousand spots are a reasonable
starting point, not a tested threshold. The command needs
`fastq-dump` (sra-tools 3.0 or later, in the core environment) and reaches
NCBI, not the ENA mirror `assemble` downloads from; a run that only ENA
holds fails with its reason in the index. Measured on a laptop, 5000 spots
of a MiSeq run streamed and sketched in about 12 s. Under `--container
docker` each tool runs in its own image with the pipe between them; on
Apple Silicon the emulated linux/amd64 `fastq-dump` crashed (observed
2026-10-19), so run the command on the host there.

#### Screening the sketches before any download

`reads-screen` judges the sketched runs against GTDB with the three checks
of `assemble --screen-reads` (see [Reads screen](#reads-screen)), so runs of
another taxon, runs whose reads match their species poorly and repeated
sequencing of one isolate can be left out before their FASTQ files are
fetched. It reads one variant of `reads_sketches/` and the submitted
taxonomy and biosample of each run from `reads.tsv`, and needs the GTDB
sketch and its lineages CSV (`--gtdb-sketch` and `--gtdb-lineages`, or
`REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`) and sourmash.

```bash
repgenr reads-sketch -wd $WD --max-spots 200000
repgenr reads-screen -wd $WD --max-spots 200000 --gtdb-sketch $DB --gtdb-lineages $LIN
repgenr reads-screen -wd $WD --max-spots 200000 --min-fraction 0.4 --accession SRR1
repgenr assemble -wd $WD --accession-file $WD/reads_screen/spots_200000/passed.txt
```

`--max-spots N` names the variant to screen (`reads_sketches/spots_N/`);
without it the `all/` variant is used when present, else the only variant
there, and with several variants the flag is required. `--accession` and
`--accession-file` restrict the selection to runs of `reads.tsv`; without
them every run is screened, and a run without a sketch in the variant (not
sketched, or failed in `reads-sketch`) is listed as `screen_failed: no reads
sketch`. The checks and their thresholds are those of the assemble screen:
`taxon_mismatch` (the GTDB genus differs from the submitted one, with the
`genus_renamed` tolerance), `low_match_fraction` (the top species fraction
is below `--min-fraction`, default 0.3) and `duplicate_isolate` (the run is
contained in an earlier run of the same biosample at `--dup-ani` or more,
default 0.999). Note that a capped sketch holds the first N spots only, so
the duplicate check compares the starts of two runs rather than the runs;
with a cap that still covers the genome several times over, the containment
of a true duplicate stays close to 1, but the estimate is less certain than
from whole runs. `--jobs` bounds the gathers run at once (default 4, or the
usable CPUs when fewer; about 0.6 GB each against the GTDB representatives
sketch).

The decisions go to `reads_screen/<variant>/`: `reads_screen.tsv` has one
row per run of `reads.tsv` with the top match, its genus, the fraction, the
duplicate it was contained in, the decision and the reason (the columns of
`screen_reads.tsv`); `<run>/screen.json` holds a run's decision with the
settings it was made with, beside the gather and `tax genome` output; and
`passed.txt` lists the accepted runs in `reads.tsv` order, in the form
`assemble --accession-file` takes, so the assembly can be limited to them.
The stage is recorded as `reads_screen` and `status` lists the variants
screened with their counts. See [output.md](output.md#reads-screen).

The standalone form `reads-screen --reads-tsv reads.tsv --sketches-dir DIR -o
OUT` (the `DIR` of `reads-sketch -o`) screens without a workdir and records
nothing, writing `OUT/<variant>/`; `--versions-out` writes the tool versions
as a YAML fragment. The Nextflow `READS_SCREEN` task uses this form (see
[Reads in Nextflow](#reads-in-nextflow)).

A repeat with the same settings gathers nothing again: a decision made with
the same thresholds and database on the same sketch is kept, and a new
sketch of a run (a later `reads-sketch --force`) or a changed threshold
screens the run again from the kept sketch, without streaming anything. A
run the tools failed on (`screen_failed: gather`) keeps no decision, so the
repeat screens it again; a selected run without a sketch is listed in the
table and needs `reads-sketch` first; `--force` screens every selected run
again. The duplicate check is made again on every call, in `reads.tsv`
order over every run with a decision made with the same settings, selected
or not, so several calls give the result of one call over every run: a run
screened later can excuse a run screened earlier as its duplicate, or stop
being its duplicate when its own decision changes. A run decided with other
settings and not selected in the call keeps its row but takes no part in
the duplicate check and is left out of `passed.txt`, with a warning that
counts such runs; screen every run with the same settings before relying on
the file. `assemble --screen-reads` does not read these decisions: its
sketch is of the whole run, and it screens again from that.

#### Assembly

`assemble` downloads each run's FASTQ files from the locations ENA lists, over
HTTPS. The files are verified against ENA's checksums, which are computed
while the bytes are written, so a file is not read a second time. `assemble`
then assembles the reads, keeps the contigs of at least `--min-contig-length`
bases (500) renamed `<run>_contig<n>`, and names the genome
`Family_genus_species_RUN.fasta` from the tokens the reads stage resolved.

**Assembler choice.** `--assembler auto` (the default) picks by platform and
layout: `skesa` for Illumina, `shovill` (SPAdes) as the alternative for paired
Illumina, `flye` for Oxford Nanopore and PacBio. `--tool-arg` passes tuning
such as `mode=nano-raw` to Flye.

**Layouts.** Long-read runs are assembled whatever layout ENA gives them (some
ONT runs are labelled PAIRED), and the polishers join a run listed as several
FASTQ files into one. A paired run that ENA lists with a third, orphan FASTQ
file is given to skesa as the pair plus the orphan file, and to shovill as the
pair only. A short-read run that ENA labels PAIRED but lists with a single
FASTQ file (mate-1 reads only, for example) is planned as single-end: `auto`
assembles it with SKESA, and `--assembler shovill`, which needs a read pair,
excuses it as `unsupported_layout` before anything is downloaded. The log
names each such run.

**Parallel runs and memory.** `--jobs` runs that many assemblies at once with
`--threads` split across them. Memory, not CPU, is the limit, so the default
is 2, or 1 as soon as a long-read run is pending. `--memory-gb` is the RAM cap
passed to SKESA and shovill (shovill accepts no less than 8).

**Resume.** Each finished run leaves a marker under `assemblies/<run>/`, so an
interrupted stage resumes without refetching. Reads are deleted after a
successful assembly unless `--keep-reads`, and the assembler's scratch unless
`--keep-files`. A run interrupted during its assembly keeps the FASTQ files
that still match their checksum, so the resume does not download them again.

The marker records the settings the run was built with (`--assembler`,
`--polisher`, `--polish-rounds`, `--min-contig-length` and the `--tool-arg`
keys its tools read). A later call, `--force` included, reuses a finished run
only while these agree. A higher `--min-contig-length` filters the finished
contigs again, which gives the same result as filtering the raw assembly. Any
other change assembles the run again, and the log names the change. To
assemble a run again under the same settings, delete `assemblies/<run>/`. A
marker written before the settings were recorded is reused as it is, filtered
at the requested floor.

**Excused runs.** The following runs are written to `excused_runs.tsv` with
the reason, and the rest proceed; the completeness guard of later stages
excuses those runs:

- a run without an ENA FASTQ mirror;
- a run whose download fails its checksum;
- a run no assembler accepts (`unsupported_platform`);
- a run whose assembler is not installed under `--assembler auto`
  (`assembler_not_installed`, with a warning naming the adapters that would
  take it);
- a run whose assembly fails.

When no run can be assembled because no assembler is installed, the stage
exits 4 instead. A repeat with the same settings skips the stage, so a run
whose download failed is retried only with `--force`; the closing warning
names such runs. What happens when every run is excused is described in
[output.md](output.md#when-every-sequencing-run-is-excused).

**Outgroup and statistics.** `--outgroup FASTA` sets a genome aside for
rooting, as `ingest --outgroup` does. Per-assembly metrics (contigs, total
length, N50, coverage from the sequenced bases) are in `assembly_stats.tsv`.

#### Polishing long-read assemblies

Long-read assemblies are polished with the run's own reads before the contig
filter. `--polisher auto` (the default) runs medaka for ONT runs (racon when
medaka is not installed) and racon (minimap2 overlaps, `--polish-rounds`
rounds) for PacBio CLR runs, and nothing for PacBio HiFi or Illumina.
`--polisher none` turns it off.

medaka needs the basecaller model. Reads basecalled with Dorado name it in
their FASTQ headers and medaka resolves it from there. Reads mirrored through
SRA have their headers rewritten and never do, so for them give the model with
`--tool-arg model=...`, or the default applies: ONT's bacterial
methylation-aware model for R10.4.1 400 bps chemistry, the usual case for
public bacterial ONT runs since 2023 (`--tool-arg bacteria=false` uses
medaka's general default model instead). The marker records which of the
three applied.

A polishing failure excuses the run with `polish_failed`;
`assembly_stats.tsv` names the polisher per genome. When `--polisher auto`
finds an adapter for a run but its tool is not installed, the run is assembled
unpolished and one warning per platform names the adapters. Unpolished ONT
assemblies carry indel errors that break genes, which CheckM2 reads as lower
completeness and higher contamination.

#### Quality scoring and classification

Two optional checks run on the assemblies: CheckM2 quality scoring and GTDB
classification with sourmash. Without a database the checks are skipped and
the log says so.

**CheckM2.** With a CheckM2 database (`--checkm2-db`, or the `CHECKM2DB`
variable CheckM2 itself reads; obtain it with `checkm2 database --download`),
every assembly is scored. The completeness and contamination reach
`selection.tsv` and the manifest, so `--keeper quality` works as it does for
GTDB genomes. An assembly below `--min-completeness` (50) or above
`--max-contamination` (10) is excused with `qc_failed`. An assembly for which
CheckM2 reports no result is kept with a warning and without quality values.

The scores are stored per run (`assemblies/<run>/checkm2.json`) with the
contigs' SHA-256, the database's resolved path and size and the CheckM2
version (the image under `--container`). A later call that changes only
`--min-completeness` or `--max-contamination` therefore applies the stored
scores without running CheckM2 again (about 5 minutes for two genomes under
emulation). A run whose contigs, database (path, size or modification time)
or CheckM2 version differ is scored again. `--force` reruns the stage but
does not bypass the stored scores; to score a run again, delete
`assemblies/<run>/checkm2.json`.

**GTDB classification.** With a GTDB sourmash sketch (`--gtdb-sketch` and
`--gtdb-lineages`, or `REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`; the
`gtdb-rs226-reps.k31-sc10k.sig.zip` sketch and its `lineages.csv` from
`https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/` serve), each
assembly is classified by `sourmash gather`. `--classifier auto` runs it when
a sketch is configured; `none` never does.

When the GTDB genus agrees with the submitted organism, the GTDB family, genus
and species name the genome file, so a reads-derived genome groups with
GTDB-downloaded ones. Otherwise the submitted name stays and
`assembly_stats.tsv` flags the genome:

- `genus_renamed` when the genus differs but the family and the species
  epithet agree, so GTDB has moved the species to another genus of the same
  family (NCBI Mycoplasmopsis arginini is GTDB Metamycoplasma arginini);
- `classifier_disagrees` for any other difference, including a shared epithet
  in another family (Klebsiella pneumoniae against Streptococcus pneumoniae).

Both are warned about, neither excuses the genome, and only
`classifier_disagrees` counts as a disagreement in the stage record. Both
lineages are kept in that table, and the sketch release is recorded in
provenance next to the metadata release.

**Preflight and memory.** The database paths and the checkm2 and sourmash
binaries are checked before any run is fetched, so a wrong path exits 2 and a
missing tool exits 4 at once rather than after the assemblies. sourmash runs
one gather per assembly. With the GTDB rs226 representatives sketch each
holds about 0.6 GB of memory (eight at once peaked at 4.3 GB), so as many run
at once as both `--threads` and `--memory-gb` / 0.6 allow (16 GB, the default,
allows 26; at least one runs). The log names the number chosen. `genome-qc`
takes the same `--memory-gb`, and the Nextflow module passes the task's
memory.

#### Adding assemblies to a GTDB selection

`assemble --append` adds the assemblies to a working directory that already
holds a selection from `metadata` and `genome` (or `ingest`). The existing
rows, genome files and outgroup stay, a run assembled earlier is replaced, and
the manifest gains the new genomes with source `sra`, so a mixed GTDB and
reads-derived set dereplicates together.

Because `metadata` and `ingest` replace the selection and the genome stage
prunes what the manifest no longer lists, both refuse to re-run while
appended genomes are present. `--drop-foreign` discards them deliberately, and
`assemble --append` can put them back afterwards. An `ingest` that takes the
same genomes again with `--from-workdir` does not drop them and is not
refused. Genomes with source `sra` that an earlier `ingest --from-workdir`
brought in are not protected this way, unless `assemble` has also run in that
working directory. To keep the two sets in separate working directories
instead, combine them with `ingest --from-workdir` (see
[Combining working directories](#combining-working-directories)).

#### Assembling in parts

`assemble --max-runs N` keeps at most N runs with a finished assembly in the
working directory, so a large `reads.tsv` can be assembled in parts.

Runs are taken in `reads.tsv` order, which the reads stage writes largest by
bases first. Runs already assembled in the working directory count toward N
first and are always kept, also when there are more of them than N; the
remaining places go to the runs still to be assembled. Runs excused before
any download (no FASTQ mirror, no accepting or installed assembler) take no
place. The other runs are listed in `excused_runs.tsv` with reason
`deferred`, so `dereplicate` runs on the assembled part; the selection and the
manifest hold only assembled runs.

N counts runs in the working directory rather than runs in one call, because
the result then depends only on `reads.tsv`, N and the finished assemblies. A
repeat with the same N assembles nothing new and is skipped as finished, while
a larger N, or no `--max-runs`, is a new request that assembles the next
deferred runs and removes their rows. A run counts toward N once assembled,
also when the quality gate later excuses it. A deferred run never counts as
judged when every run is excused (see
[output.md](output.md#when-every-sequencing-run-is-excused)), and `--append`
takes `--max-runs` in the same way. `status` shows the stage as, for example,
`20 assembled, 35 deferred`.

`assemble --accession SRR... --accession-file FILE` names the runs to
assemble instead of a count: the other pending runs are deferred in the same
way, finished runs are kept whether named or not, and an accession not in
`reads.tsv` is refused, as is a file that lists no accessions (an empty
`passed.txt` means the screen excused every run). The `passed.txt` of
`reads-screen` is such a file, so a screened selection can be assembled
without editing `reads.tsv`.

The alternative is to bound the selection itself: `reads --max-runs N` keeps
the N largest runs in `reads.tsv`, and a later `reads` call with a larger N
followed by `assemble` assembles only the new runs, since finished runs are
reused. In the Nextflow pipeline, where each run is one `assemble-run` task,
the bound is set on the selection with `--reads_args '--max-runs N'`;
`assemble-run` takes no `--max-runs`.

#### Reads sketch

While a run is assembled, `assemble` also sketches its reads with sourmash:
one `sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund` call on all
FASTQ files of the run (the two files of a pair give one signature per k-mer
size), named by run accession and written to `assemblies/<run>/reads.sig.zip`.
The k-mer sizes and the scale are those of the genome sketches; abundances are
kept because read k-mers also carry coverage.

The sketch starts once the FASTQ files are verified and runs in a thread
beside the assembler, so the assembler is not delayed. The assembler keeps its
full share of `--threads` and the single-threaded sketch runs alongside it,
so each run uses one thread more than its share while the sketch lasts.
Taking that thread from the assembler was measured to cost more than it saved
(SKESA on a 134 MB MiSeq run took 55 s on three threads against 41 s on four).
The sketch ends before the reads are deleted. It is written to a temporary
file in the run directory and renamed into place, so an interrupted sketch
leaves no `reads.sig.zip`.

`--reads-sketch/--no-reads-sketch` controls it. By default the reads are
sketched when sourmash can run; otherwise one INFO line per stage says so.
`--reads-sketch` requires sourmash and exits 4 before anything is downloaded
when it is missing, as `--sketch` does. A failed sketch never fails the
assembly: it is logged as a warning and the marker records no sketch.

The done marker records the file name (`reads_sketch`), the parameter string
and the sourmash version, or null for each; `assembly_stats.tsv` has a
`reads_sketch` column (1 or 0). An excused run keeps no reads sketch, a
quality-excused one included. A finished run whose marker records no sketch,
for example one assembled before this option existed, is not assembled again.
When its FASTQ files were kept (`--keep-reads`) and still match their
checksums, the next `assemble` call sketches them and updates the marker;
otherwise the sketch stays empty and one INFO line counts such runs. When the
stage is otherwise up to date, that call needs `--force`.

The reads sketch is not a genome sketch. It is kept beside the assembly and
is never copied into `sketches/`, which holds one signature per genome for
the consumers described under [Genome sketches](#genome-sketches).

#### Reads screen

`assemble --screen-reads` (off by default) checks each run's reads sketch
before the run is assembled. Use it to exclude, before the assembly cost is
paid, runs of another taxon, runs whose reads match their species poorly, and
repeated sequencing of one isolate.

Within the per-run worker the order is: fetch and verify the FASTQ files,
sketch the reads, screen the sketch, and only then choose and start the
assembler. Under the screen the reads sketch is therefore awaited rather than
run beside the assembler. A run that fails a check is excused in
`excused_runs.tsv` with step `screen` and one of three reasons:

- `taxon_mismatch`: `sourmash gather` of the reads sketch (k=31) against the
  GTDB sketch, resolved to a lineage with `sourmash tax genome` and the
  lineages CSV (the commands of the post-assembly classifier), gives a genus
  other than the submitted one. The classifier's tolerance applies: a genus
  that differs while the family and the species epithet agree
  (`genus_renamed`) passes. A run with no GTDB match above the gather
  threshold (50 kbp) is also excused with this reason. A run whose submitted
  genus is unknown is not tested for it.
- `low_match_fraction`: the abundance-weighted fraction of the reads assigned
  to the top species (`f_weighted_at_rank` of `tax genome` at species rank) is
  below `--screen-min-fraction` (default 0.3). The excuse reads `top species
  fraction X below --screen-min-fraction Y`.
- `duplicate_isolate`: the run's reads are contained in those of a run
  accepted before it, of the same biosample, at an ANI estimate of
  `--screen-dup-ani` (default 0.999) or more.

**The match fraction.** The fraction mostly measures the distance between the
run's strain and the GTDB representative of its species. Host, contaminant
and repeated sequencing-error k-mers lower it as well, so a low value does not
by itself show contamination. `tax genome` is run with a containment threshold
of 0 here, since most distinct k-mers of a reads sketch are sequencing errors
and the unweighted fraction of a pure isolate stays far below sourmash's
default threshold of 0.1. The fraction of a pure isolate depends on how close
the run's strain is to the GTDB representative of its species, since the
reference sketch holds one genome per species, and on how many error k-mers
are seen more than once. Against the GTDB rs226 representatives sketch (k=31,
scaled=10000), four Illumina MiSeq isolates gave 0.58 (Mycoplasmopsis
arginini), 0.42 and 0.43 (two Mycoplasma mycoides runs) and 0.39 (Escherichia
coli); the default of 0.3 lies below that range. A species far from its
representative can fall lower, so check `screen_reads.tsv` before relying on
the default.

**Duplicate detection.** Runs are compared in `reads.tsv` order, which the
reads stage writes largest first, so of two duplicates the larger run is
kept. A run waits for the decisions on the earlier runs of its biosample, so
the result does not depend on which worker finishes first. The comparison is
limited to runs accepted in the same call and finished runs with a reads
sketch.

It uses the k=31 signature. K-mers seen at least a tenth of the
abundance-weighted median (and at least twice) are taken as genomic, and
their containment C in the k-mers the other run saw at least twice gives the
estimate C^(1/31), as sourmash's containment ANI does. When the two sketches
have different scaled values, both are first reduced to the larger one. At
k=31, `--screen-dup-ani` 0.999 corresponds to a containment of about 0.97
(0.999^31), roughly 0.1 percent divergence. At scaled=1000 a genome of 1 Mb
gives about 1000 such k-mers, so one missing k-mer already moves the estimate
from 1.0 to about 0.99997: the estimate cannot separate 0.9999 from 1.0 for a
genome of that size.

Runs of different biosamples are never compared: clonal isolates of
different samples are distinct samples (two Mycoplasma mycoides runs of
different biosamples gave an estimate of 1.0000 and are both kept). A run
that passes the taxon and fraction checks from a kept sketch is fetched
before its duplicate check, which runs in the worker.

**Requirements and failures.** The screen needs sourmash and the GTDB sketch
and lineages CSV (`--gtdb-sketch` and `--gtdb-lineages`, or
`REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`, as for the classifier).
Without a sketch or lineages file, or with `--no-reads-sketch`, the stage
exits 2 before anything is downloaded, and without sourmash it exits 4. The
classifier setting does not matter: the screen runs with `--classifier none`
as well. The post-assembly classifier runs unchanged on the accepted
assemblies.

A tool failure while screening (a failed reads sketch or gather) excuses the
run with `screen_failed`, which, like a failed download, does not count as a
judgement of the run's data when every run is excused (see
[output.md](output.md#when-every-sequencing-run-is-excused)); the three
reasons above do count.

**Cost and outputs.** The cost is one gather and one `tax genome` call per
run, 30 to 40 s against the GTDB rs226 representatives sketch (measured with
two to four at once), plus the time the sketch would otherwise have
overlapped with the assembler. Gathers of concurrent runs are bounded by the
classifier's memory budget (`--memory-gb` / 0.6 GB). An excused run keeps its
reads sketch, its decision in `assemblies/<run>/screen.json` and the gather
and `tax genome` output under `assemblies/<run>/screen/`, so the reason can
be inspected; its reads are removed as usual. `screen_reads.tsv` in the
working directory lists every screened run (see
[output.md](output.md#working-directory-layout)), and `assembly_stats.tsv`
marks the genomes that passed (`screen` column).

**Resume.** Each decision records its settings (both thresholds and the
resolved paths of the GTDB sketch and lineages CSV). A repeat with the same
settings screens nothing again: an excused run is excused from its recorded
decision without a download, and a finished run keeps its decision. When a
setting changed, the runs are screened again from their kept reads sketches,
without fetching their reads. An excused run is fetched again only when it
now passes, and a finished run that now fails is excused while its contigs
and marker stay for a later call. A finished run assembled without the screen
is screened from its reads sketch; one without a reads sketch is not
screened, and one INFO line counts such runs.

With `--max-runs`, a run excused by the screen while it is processed does not
give its place to a deferred run in the same call, since the places are fixed
before any run is fetched. The next call excuses it from its recorded
decision before the places are counted, so a deferred run takes its place
then.

#### Stateless steps

The same work is available as three stateless steps. The Nextflow reads mode
runs them as separate tasks, and they also serve a scheduler of your own:

```bash
repgenr reads -wd $WD -tg mycoplasmopsis --platform illumina --max-runs 20
# one task per run: fetch, assemble, filter contigs
repgenr assemble-run --reads-tsv $WD/reads.tsv --run SRR25474756 -o asm/SRR25474756 \
    --assembler auto -t 8 --memory-gb 16 --min-contig-length 500
# one batch: CheckM2 and/or the classifier over every finished run directory
repgenr genome-qc --assemblies asm -o qc --checkm2-db $CHECKM2DB \
    --gtdb-sketch gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb-lineages lineages.csv \
    --classifier auto -t 16
# the genome contract: genomes/, selection.tsv, assembly_stats.tsv, excused_runs.tsv
repgenr reads-gather --reads-tsv $WD/reads.tsv --assemblies asm --qc qc -o out \
    --min-completeness 50 --max-contamination 10
```

None of the steps touches a working directory or the manifest.

- **`assemble-run`** writes `contigs.fasta` and the `assembly.ok` marker into
  its `--out` directory, or `excused_runs.tsv` when the run has no FASTQ
  mirror, an unsupported platform, or fails to download, assemble or polish.
  `--polisher`, `--polish-rounds`, `--keep-reads`, `--keep-files`,
  `--tool-arg` and `--reads-sketch/--no-reads-sketch` act as on `assemble`;
  the reads sketch goes to `reads.sig.zip` in the same directory. In
  Nextflow, the parameter `--assemble_args '--no-reads-sketch'` turns it off.
- **`assemble-run --screen-reads`** (with `--screen-min-fraction`,
  `--gtdb-sketch` and `--gtdb-lineages`) applies the taxon and fraction checks
  of the reads screen and writes `screen.json` beside the sketch. The
  duplicate check needs the other runs and is made only by `assemble`, since
  each `assemble-run` task sees one run. In Nextflow,
  `--assemble_args '--screen-reads'` turns it on, and `conf/modules.config`
  then adds `--gtdb-sketch` and `--gtdb-lineages` from the pipeline
  parameters of the same names.
- **`genome-qc`** reads a directory of such run directories (`--assemblies`)
  and writes `quality.tsv` and `classification.tsv` keyed by run accession;
  it needs at least one database.
- **`reads-gather`** applies the quality gate and the naming policy above
  from the optional `--qc` directory, collects the `screen.json` files into
  `screen_reads.tsv`, and writes an empty `outgroup_accession.txt`, since the
  reads chain has no outgroup at this point.

### Viruses

The viral route selects virus genomes with `vmetadata` and stages them with
`vgenome`. It sources from NCBI Virus by default (via the `datasets` CLI);
`vmetadata --source bvbrc` uses the legacy BV-BRC FTP path instead.

```bash
WD=./hav
repgenr vmetadata -wd $WD --target hepatovirus            # NCBI Virus (default)
repgenr vgenome   -wd $WD --target-genus Hepatovirus      # add --group-segments for segmented viruses
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
# or: repgenr run -wd $WD --viral --target hepatovirus -tg Hepatovirus --treebuilder mashtree
```

`--virus` passes virus-tuned settings to dRep (`--tool drep --virus`); the
other dereplicators do not read it, and `dereplicate` warns when it is given
with one of them.

#### The viral outgroup

`vgenome` picks an outgroup by itself: a record of a sister species with
enough genomes, chosen by distance (mashtree by default).
`--outgroup-accession` pins it to a downloaded record instead, an accession
on the NCBI Virus path or a record id on BV-BRC, and
`run --viral --outgroup-accession` forwards the same choice. With
`--group-segments` the search still runs, with the kept records' length span
widened by 15 percent as its window, and the outgroup is one record of the
sister species rather than a grouped isolate. `--no-outgroup` leaves the tree
unrooted, and a run that ends without an outgroup removes the outgroup an
earlier run left in `outgroup/` and `outgroup_accession.txt`.

#### Species assignment

On the NCBI Virus path the species of a record is the current species of its
taxid in NCBI Taxonomy. After the download, `vmetadata` looks up every
distinct taxid once with `datasets summary taxonomy taxon` and stores the
species, genus and family in `virus_records.json`, so `vgenome` needs no
network.

The organism name is not used for this. NCBI organism names are often
strain-level or earlier names (`Sabia virus` in `Mammarenavirus brazilense`,
`Hantaan virus CGAa1011` in `Orthohantavirus hantanense`), and one species
can carry several of them (Junin, taxid 2169991, as
`Argentinian mammarenavirus` and `Mammarenavirus juninense`). The lineage in
the NCBI Virus report is not a reliable substitute either. It nests species
that NCBI Taxonomy keeps as siblings (Maguari virus under
`Orthobunyavirus cacheense` > `Orthobunyavirus maguariense`, while its species
is `Orthobunyavirus maguariense`), and it can lack the current species
altogether (Murutucu virus, species `Orthobunyavirus maritubaense`).

When the lookup fails, or gives no species for a taxid, the species comes from
the report lineage, and the log says so. The species is then the shallowest
lineage name made of a genus (or the genus above a subgenus) and one
lower-case epithet, such as `Mammarenavirus brazilense`. Earlier names below
it that start with the genus (`Hepatovirus A` under `Hepatovirus ahepa`) are
not taken, and all records of a taxid take the binomial found on most of
them. A record without either keeps its organism name as the species. Each
record's `species_source` (`taxonomy`, `lineage` or `organism`) says which
rule applied, and the `vmetadata` record counts them. Organism names of the
form `<Genus> sp.` (`Orthohantavirus sp.`) are a species in NCBI Taxonomy, so
the unclassified records filed under one such name share one species token.

The species sets the species token of the canonical filename, the
`--target-species` match, the per-species medians of the length window and
the grouping of outgroup candidates, so strains of the target species are not
taken as sister species. A `--target-species` value also selects the species
of every record whose organism name it is (`-ts "Argentinian mammarenavirus"`
selects all of `Mammarenavirus juninense`), together with any species token
it matches, and the log says so. The organism name stays in the
`description` column of `virus_metadata_base.tsv` and in `--target-serotype`
matching.

#### Segmented viruses

`--group-segments` groups records that share a species and an isolate name
and carry at least two distinct segments into one genome per isolate. Use it
for segmented viruses, whose segments NCBI Virus holds as separate records.

Segment labels are free text in NCBI Virus and are compared after
normalisation:

- the text after a `;` is dropped;
- words that only name the molecule or the word segment (`RNA`, `DNA`,
  `segment`, `genome`, `component`, `circular`) are skipped;
- the first remaining word is kept in upper case;
- small, medium, middle and large become S, M, M and L.

So `M`, `M; medium` and `middle` are one segment, `S RNA` is `S` and `DNA-A`
is `A`, while `RNA 1` and `RNA 2` stay distinct. A label of `Unknown` counts
as no label.

Each isolate keeps one record per segment (complete before partial, then the
longest, then the lowest accession), since NCBI Virus often holds several
submissions of one segment; records without a segment label stay single
genomes. A grouped isolate's genome carries a synthetic `iso-` token as its
accession (with its first member accession appended when two isolate names
would give the same token). `segments.tsv` records the member accessions
behind it with the normalised segment and the label as submitted, and
`tree2tax` lists the members under the isolate's leaf in `genomes_map.tsv`.

#### Viral length filtering and over-represented species

The viral selection step keeps records whose genome length falls inside a
window around a center value. The default center, `--length-method
median_of_medians`, is a deliberate defense against over-sequenced species:
the median length is computed per species first, and the window center is the
median of those per-species medians, so a species with 900 outbreak records
carries exactly one vote. Switching to `--length-method mean` averages every
record and lets the most-sequenced species set the window; use it only when
the selection is known to be balanced. `--length-range` overrides the window
entirely and `--length-all` disables the filter.

#### Older working directories and source switches

A working directory whose `vmetadata` ran before the species came from NCBI
Taxonomy keeps its genomes and filenames while its stages skip. When
`vgenome` reruns on it (new arguments or `--force`), it refuses the old
records and asks for one `repgenr --force vmetadata` with the same arguments;
that download renames the viral genomes, and `vgenome` and the later stages
then rerun.

`vmetadata` records the source and target of `virus_download_wd/download.fa`
in `download.source`. The BV-BRC source reuses the group FASTA only for the
same target and downloads it again otherwise. Switching a working directory
between the two sources removes the other source's tables, so `vgenome` reads
the latest download. A BV-BRC working directory written by an earlier version
has no `download.source`, so its first re-run downloads the group FASTA
again.

### Genome sketches

`genome`, `vgenome`, `ingest` and `assemble` end by writing a sourmash sketch
of each genome, outgroup included. The sourmash dereplicator, `glance`, tree
builder and classifier read these sketches instead of sketching the genomes
again.

Each sketch is written to `sketches/<name>.sig.zip`, where `<name>` is the
genome's file name without its FASTA suffix (the name it carries in
alignments and trees). Each file holds three DNA signatures, at k=21, k=31
and k=51 with scaled=1000 and no abundance, named after the genome.

**When sketches are written.** By default the step runs when sourmash can run
(on the PATH, or through the container backend) and is otherwise left out
with one INFO line. `--sketch` requires sourmash and exits 4 without it,
before any genome or sketch is written (the run log may already exist), and
`--no-sketch` leaves it out. A failed sketch is a warning under the default
and an error with `--sketch`.

A sketch is written again only when the genome FASTA (its SHA-256) or the
sketch parameters changed, so a rerun of a stage sketches only new and
changed genomes. `assemble --append` sketches only the genomes it adds, and
`ingest --from-workdir` copies the sketch of a source genome whose recorded
digest matches. Sketches of genomes that leave the set are removed with them.

The sketches are not stage outputs: a genome set without them is complete,
the flag does not change the resume fingerprint, and a stage skipped as
already complete does not sketch (with `--sketch` the skip message names
`repgenr sketch`).

**Sketching an existing working directory.** `repgenr sketch -wd $WD [-t N]`
writes the missing and stale sketches, removes those of genomes no longer in
the set, and logs one line with the counts (present, written, stale replaced,
removed). `repgenr --force sketch -wd $WD` writes every sketch again; sketches
that were current are counted as `replaced (--force)`, apart from stale ones.
It is recorded as a stage, so `status` lists it, and `status` shows
`sketches: n/m` (genomes with a sketch, genomes present) on the line of the
stage that wrote the genome set once `sketches/` exists.

**Consumers.** The sourmash tools read these sketches:
`dereplicate --tool sourmash` (both the branchwater and the `sourmash compare`
back-end), `glance --tool sourmash`, `phylo --treebuilder sourmash` (outgroup
included) and the sourmash classifier of `assemble`. Each selects the
signature at its k-mer size with `-k`, so `ksize` 21, 31 or 51 at `scaled`
1000 is served from the sketches. Before comparing, the stage checks each
genome's sketch against its record (the FASTA SHA-256 and the parameters),
writes the missing and stale ones to `sketches/` and records them, and logs
one line such as `dereplicate: sketches: 50 reused, 0 written`. The
representatives under `derep/representatives/` are matched to their genomes
by file name.

The SHA-256 of a genome is kept in `sketches/.digests.json` with the file's
size and modification time and is computed again only when either changes,
the same file metadata the resume fingerprint uses for genome directories. An
edit that keeps both is not seen, and `repgenr --force sketch` hashes and
sketches every genome again.

The tool sketches into its own working directory, as before, when
`--tool-arg ksize=` or `scaled=` asks for other parameters (the log names the
reason), for a genome that is not part of the genome set (no manifest row),
and when sourmash cannot write a sketch.

**Assemblies.** `assemble` sketches each accepted assembly once, beside its
contigs as `assemblies/<run>/contigs.sig.zip` (named by run accession, reused
while the contigs are unchanged). The classifier gathers with that sketch,
and the sketch step copies it to `sketches/` under the genome's record name
once the genome is named. An assembly excused after classification (for
example `qc_failed`) gets no sketch under `sketches/`. With `--no-sketch` the
classifier sketches for itself. The reads sketch is separate (see
[Reads sketch](#reads-sketch)).

### Counting genera, species and samples

`repgenr census` counts what a taxon holds before a run, or what a working
directory holds after its entry stage. Use it to choose a target and a route
before downloading, and to check a selection afterwards. It is read-only: it
records no stage and writes nothing into a working directory. A family is
counted per genus and a genus per species, with the totals on the first line.
The command has three modes, chosen by its arguments.

#### A taxon, without a working directory

`-tf FAMILY` or `-tg GENUS` counts the GTDB genomes of the taxon and its GTDB
species representatives.

- `--source api` (the default) sends one request to the GTDB API, as
  `metadata --source api` does.
- `--source table` reads the GTDB metadata table of `--release`
  (`--gtdb-version` bac120 by default, or ar53) with the parser of the
  metadata stage. The table is downloaded once into
  `$REPGENR_CACHE_DIR/gtdb` (default `~/.cache/repgenr/gtdb`) and reused, or
  `--metadata-path` names a local copy. A table census parses the whole GTDB
  table, which takes about 20 s for release 232.

`--runs` adds the ENA whole-genome sequencing runs under the taxon, as the
reads stage finds them: runs, distinct biosamples, and runs per platform
(Illumina, ONT, PacBio, and `other` for BGISEQ, Ion Torrent and the rest).
The count is logged first, with a warning above 100000 runs; the census
requests only the four run fields it reads, so its search is smaller than
that of the reads stage.

The runs are grouped by the NCBI lineage of their taxid (one lookup per
distinct taxid) and join the GTDB row of the same name. A row that holds runs
only is marked `(NCBI)` (or `(no genus, NCBI)` for taxids without a genus),
listed last and left out of the GTDB totals; such a taxon may be a GTDB taxon
under another name. For a GTDB name with a suffix (`Bacillus_A`), the runs
are those of the NCBI taxon without it (`Bacillus`, the prokaryote genus). In
a genus census they join its species rows by epithet, a note says so, and the
header names the taxon of the runs (`19772 runs of NCBI genus Bacillus`).
Names are shown as GTDB and NCBI write them (`Bacillus_A thuringiensis_S`,
`Francisella sp. LA112445`).

```bash
repgenr census -tg Francisella
repgenr census -tf Francisellaceae --runs
repgenr census -tf Francisellaceae --source table -r 232.0 --gtdb-version bac120
```

```text
Francisellaceae: 9 genera, 48 species, 1205 genomes (48 representatives), 1681 runs, 1473 biosamples
name                    species  genomes  representatives  runs  biosamples  illumina  ont  pacbio  other
Francisella                  26     1157               26  1667        1462      1413   41      56    157
Caedibacter                   9       28                9     0           0         0    0       0      0
M0027                         4        8                4     0           0         0    0       0      0
...
(no genus, NCBI)              0        0                0     3           3         3    0       0      0
Allofrancisella (NCBI)        0        0                0    11           8         7    3       1      0
```

`--viral --target TAXON` counts the NCBI Virus sequence records of a viral
taxon with the request vmetadata sends, but for the metadata report only (no
sequences); `--host`, `--complete-only` and `--released-after` filter it as in
vmetadata. The columns are sequences, complete sequences, isolates (the
genomes `vgenome --group-segments` would form, which equals the number of
sequences for an unsegmented virus) and whether the row holds a segmented
virus. Isolates are counted before vgenome's length filter, so they are an
upper bound on what vgenome selects. A genus target is counted per species
and a family target per genus; `-tf` and `-tg` narrow the report. BV-BRC
offers no metadata-only report, so `--viral --source bvbrc` without `-wd`
exits 2; run `vmetadata --source bvbrc` and count its working directory
instead. `--runs` is not offered for viral taxa.

```text
$ repgenr census --viral --target hepatovirus
Hepatovirus: 1 genus, 31 species, 12807 sequences (369 complete, 12807 isolates)
name                                   sequences  complete  isolates  segmented
Hepatovirus ahepa                          12485       346     12485         no
Rodent hepatovirus                            65         1        65         no
...
```

#### A working directory after metadata or vmetadata

With `-wd`, the census counts the candidates the entry stage chose from:

- the GTDB table under the recorded target (the table named in the metadata
  record, the release table in the working directory, or `--metadata-path`);
- the GTDB API answer that `metadata --source api` keeps in
  `gtdb_api_genomes.tsv`;
- the NCBI Virus records of `virus_records.json`;
- for BV-BRC, the record sets of the taxonomy names in
  `metadata_ncbi_taxnames_data.json` (the complete, isolate and segment
  columns are then empty, and `median_length` is the median of the per-taxid
  median lengths).

The recorded target is used; `-tf` and `-tg` narrow it. An older API working
directory without `gtdb_api_genomes.tsv` is counted from its `selection.tsv`,
with a note.

#### A working directory after genome, vgenome, ingest or assemble

The census counts the genomes of `selection.tsv` that are present in
`genomes/`, one column per manifest source (gtdb, sra, ncbi_virus, bvbrc,
local), the GTDB species representatives, and the dereplication clusters
whose representative falls in the row when `derep/clusters.tsv` exists. The
outgroup is not counted and is named on the first line. When the candidates
of the entry stage are also present, a `candidates` column shows how many
each row had before selection; a species with candidates and no selected
genome appears with 0 genomes. When the entry-stage data cannot be read (for
example a `virus_records.json` written by an earlier version), the census
counts the selection without the `candidates` column and says why in a note.

```text
$ repgenr census -wd hisp
Francisella: 1 genus, 1 species, 4 genomes (1 representative), 4 candidates, 1 cluster; outgroup GCF_000815225.1 excluded
name                      candidates  genomes  gtdb  representatives  clusters
Francisella hispaniensis           4        4     4                1         1
Note: Candidates are the genomes of the GTDB table under the metadata target.
```

#### Output formats and errors

`--tsv PATH` also writes the rows (not the totals) to a TSV file, and
`--json` prints one JSON object with the taxon, mode, totals and rows on
stdout and nothing else; both formats are described in
[output.md](output.md#census-table). The flags that select a taxon query
(`--viral`, `--target`, `--source`, `--release`, `--gtdb-version`, `--runs`,
`--host`, `--complete-only`, `--released-after`) exit 2 with `-wd`, whose
record names the source and target. A command without a taxon or a working
directory exits 2, and a `-wd` that does not exist and a failed request exit
3.

## Dereplication and representative selection

`repgenr dereplicate` clusters the genome set by ANI and keeps one
representative genome per cluster. It runs after any entry path and before
`phylo`. The dereplicator is chosen with `--tool` (drep, skder, galah or
sourmash; see [choosing-tools.md](choosing-tools.md) for which suits a
dataset):

```bash
repgenr dereplicate -wd $WD --tool skder -t 16
```

The clusters are written to `derep/clusters.tsv` and the representatives to
`derep/representatives/`. This part describes how the representative of each
cluster is chosen, how large sets are dereplicated in chunks, and how to
inspect the result.

### Representative selection

After the chosen dereplicator clusters the genomes, a keeper step re-picks
each cluster's representative (`src/repgenr/stages/derep_keeper.py`). Use the
default, `--keeper quality`, unless there is a reason to keep the tool's own
pick or the GTDB species representatives. The keeper step corrects the
tendency of connectivity-based tools to keep the most-sequenced (not the
best-quality) genome in a cluster.

`repgenr dereplicate` and `repgenr run` accept `--keeper quality|gtdb|tool`
(default `quality`), and so do the Nextflow steps `dereplicate-chunk` and
`dereplicate-merge` (`--derep_keeper`). `--keeper tool` restores the
adapter's own pick.

**Quality values.** The values come from the manifest:

- GTDB selections carry CheckM values (`--source api` fetches them from each
  genome's card, one request per selected genome);
- `assemble --checkm2-db` scores assemblies;
- `ingest --selection` reads the `completeness` and `contamination` columns.

The Nextflow steps read the same columns from `selection.tsv`. When no genome
under `genomes/` has quality the stage warns, and `repgenr.yaml` records
`keeper_effective: tool` next to the requested `keeper` and the swap count.

#### The keeper score

The score follows dRep's default weights for the terms that need no further
tool: `completeness - 5 x contamination + 0.5 x log10(N50)`, with the N50
read once per genome from its FASTA (gzipped or not). The N50 term separates
genomes whose CheckM values differ by less than CheckM's precision: on 30
GTDB *F. tularensis* genomes, all scored between 99.44 and 100, a 29-contig
draft at 100/0.00 outscored closed genomes at 100/0.03 on CheckM alone, but
not with the N50 term. Genomes are compared in this order:

1. higher score;
2. higher completeness;
3. lower contamination;
4. higher N50;
5. genome filename, in sort order.

Only genomes with completeness and contamination compete; a genome without
them never becomes the representative. When the adapter's representative has
values, the best genome in the cluster wins, that representative included, so
the result does not depend on which genome the tool picked. When it has none,
a scored genome replaces it only if it is high quality by the MIMAG
thresholds (completeness above 90, contamination below 5). A scored fragment
therefore does not displace an unscored complete genome; the log names the
clusters left with the tool's representative for this reason. A cluster
without any scored genome keeps the adapter's choice.

The working-directory stage applies the keeper once, to the final clusters of
a chunked run. The Nextflow steps apply it within each chunk and again after
the merge, so the merge pass compares the chunk keepers rather than the
tools' picks. Both end with the keeper of each final cluster, but the
clusters themselves can differ slightly, since the merge pass compares
different genomes. Each chunk records the N50 of its scored genomes in
`genome_n50.tsv`, so the merge step scores members whose files it does not
receive.

The keeper rule (N50 term, high-quality condition, tie order) is not part of
the resume fingerprint. A working directory dereplicated by an earlier version
with a different rule keeps its representatives until `dereplicate` is rerun
with `--force`.

#### GTDB species representatives

`--keeper gtdb` keeps the GTDB species representatives. The metadata stage
records which selected genomes GTDB marks as the representative of their
species (`gtdb_genome_representative` equal to the accession in the table,
`gtdbIsRep` in the API), in the `gtdb_representative` column of
`selection.tsv` and in the manifest.

Within each cluster the keeper is then a GTDB representative when the cluster
holds one. A cluster without one is treated as under `--keeper quality`, and
keeps the tool's pick when no member is scored. A cluster can hold several
GTDB representatives when the ANI threshold joins GTDB species (GTDB delimits
species at about 95% ANI, so a `--secondary-ani` below that can do so). One
representative per cluster is kept, chosen among them by the quality rule
(score, tie order and the high-quality condition for an unscored
representative as above), and a warning names the others, which stay
contained; raise `--secondary-ani` to keep them apart. With `--reduce`, a
GTDB representative is also preferred when representatives of one taxon are
merged.

`repgenr.yaml` records `keeper_effective: gtdb` when the manifest flags at
least one genome, and otherwise `quality` or `tool` with a warning. Only GTDB
selections carry the flag: `ingest` (unless its `--selection` has the
column), `reads`/`assemble` and the viral path leave it at 0, so
`--keeper gtdb` acts there as `--keeper quality`.

A working directory selected by an earlier version has no flag recorded;
`repgenr metadata --force` with the same options records it. `genome` then
runs again because `selection.tsv` changed, and downloads nothing already
present; `dereplicate` runs again because the manifest now flags genomes. A
manifest from an earlier version is upgraded in place (schema version 3) when
a stage opens it, and an earlier version then refuses to open it.

#### Quality values passed to the dereplicator

The same manifest values also reach the dereplicator, whichever keeper rule is
chosen (the Nextflow chunk and merge steps read them from `selection.tsv`),
when every genome of the run has both values. dRep receives them as
`--genomeInfo` and then runs without CheckM; galah receives them as
`--genome-info` and ranks genomes by quality. The decision is taken once for
the run, so every chunk and the merge pass use the same source.

When some genomes lack values, none are passed and a warning names them. dRep
then runs CheckM, and galah, which otherwise keeps the first listed genome of
a cluster, receives the genomes by descending file size, so a complete genome
rather than a fragment becomes the representative.

#### One representative per taxon

`--reduce species|genus` collapses the ANI representatives to one per taxon
after dereplication. The representative of the largest cluster is the
default keeper of a taxon; with `--keeper quality` and scores known, the
keeper rule above re-picks among the taxon's representatives.
`--target-reps N` searches the secondary ANI to land near N representatives.

Both options exist on `dereplicate` and, since the merge step is where the
final set is decided, on `dereplicate-merge`, which the Nextflow layer drives
through `--derep_reduce` and `--derep_target_reps`. The merge step takes the
taxonomy from `selection.tsv` when one reaches it and from the canonical
genome filenames otherwise.

#### Finer resolution inside a lineage

One secondary ANI applies to the whole set. A study often needs coarse
coverage of a species and more genomes inside one lineage, such as a
clinically relevant genotype or an outbreak clade. `--refine` keeps that
lineage at a finer resolution without changing the rest of the result.

The lineage is named by seeds. Give one or a few genomes known to belong to
it (`--refine SEED`, repeatable: a filename, a filename without its suffix,
or an accession) and a finer ANI (`--refine-ani`, above `--secondary-ani`).
The dereplication runs as usual; then every cluster that holds a seed is
taken as the lineage, and its genomes, the representative included, are
dereplicated again with the same tool and settings at the finer ANI. The
clusters this gives replace the coarse one. Nothing outside the seed clusters
changes, and a seed alone in its cluster has nothing to refine and is
reported.

```bash
# In one run: the species at 0.99, the lineage of GCA_000123 at 0.999
repgenr dereplicate -wd $WD --tool skder --refine GCA_000123 --refine-ani 0.999

# After the fact, on the result on disk, without repeating the coarse pass
repgenr derep-stock -wd $WD --action pack --name coarse
repgenr derep-refine -wd $WD --refine GCA_000123 --refine-ani 0.999
```

`derep-refine` takes the tool and thresholds from the completed `dereplicate`
record, so it needs one (a record restored by `derep-stock` counts). It
rewrites the dereplication outputs in place and records a `derep_refine`
stage; `phylo` and `tree2tax` then rebuild from the refined set. A later
`dereplicate` (forced, or with changed inputs) recomputes the coarse result
and discards the refinement, so store the result you want to keep with
`derep-stock` first. `repgenr run ... --refine SEED --refine-ani F` passes
both to its `dereplicate` stage. A seed given as an accession needs the
version suffix the filename carries (`GCA_000123.1` when the file is named
so); unknown seeds, the outgroup, a genome the dereplicator rejected and an
ANI that is not finer are refused before the coarse pass runs.

The keeper rule chooses the representative of each refined cluster as of any
other. `--reduce species|genus` leaves the refined clusters apart, since they
share a species by construction; the species then keeps one representative
outside the lineage and the refined ones inside it. `--refine` cannot be
combined with `--target-reps`, whose search sets the secondary ANI itself.
The outgroup is never part of a lineage.

Several lineages, each at its own ANI, are given in a file:
`--refine-file refine.tsv`, a TSV with the columns `seed` and `ani`, one
seed per row; the seeds of one ANI form one lineage, and a cluster claimed
by two ANIs is refused (seeds whose coarse cluster was refined before count
as one cluster). The file can be combined with `--refine` and `--refine-ani`,
which apply to the `--refine` seeds only, and its content is part of the
resume inputs, so editing it reruns the refinement; the path is stored as
an absolute path. Blank lines and lines starting with `#` are skipped. The `ani` column of `refined.tsv` and the
`refine_lineages` entry of the stage record hold each lineage's ANI.

```text
seed	ani
GCA_000123.1	0.999
GCA_000456.1	0.999
GCA_000789.1	0.9995
```

`derep-refine` can be applied more than once, for another lineage or for the
same seed at another ANI. A seed refined before brings in its original
coarse cluster, read from `refined.tsv`, so refining it again at a finer ANI
starts from the same genomes and not from the fine cluster it sits in now;
the rows of lineages left untouched are kept. On a result reduced with
`--reduce`, a cluster on disk is a whole taxon group, and the lineage is
refined as such, with a warning. The seeds and ANI of a `derep-refine` live
in its own stage record; `derep-stock` stores the `dereplicate` record and
`refined.tsv`, so after an unpack the file alone describes the refinement.

The lineage genomes are listed in `derep/refined.tsv` (genome, seed, the
coarse representative, the final representative and status), and the
`refined` column of `cluster_summary.tsv` marks the clusters the refinement
produced (see [output.md](output.md)). The stage record holds the seeds, the
ANI and the counts. `derep-stock` stores and restores the file with the run.

#### The sourmash dereplicator and its sketch cache

`--tool sourmash` reads the genome sketches of the working directory (see
[Genome sketches](#genome-sketches)), so `--target-reps` compares the same
sketches at each search step. Without them (other `ksize` or `scaled`, or
genomes outside the set), it keeps its own signatures in a sketch cache.
`--target-reps` shares one cache across its search steps;
`--tool-arg sketch_cache=DIR` sets a directory that persists across runs.

With the branchwater plugin, cache entries are matched by genome file name,
size and modification time; without it, a per-genome signature is reused
when it is not older than the genome file. A genome replaced under the same
name is therefore sketched again, except when the new file carries an older
timestamp (copied with `cp -p` or `rsync -t`) or, with the plugin, the same
size and timestamp. Empty the directory after such a replacement.

### Chunked dereplication

`dereplicate -s N` (`--process-size`) dereplicates in chunks. Use it for sets
too large for one pass, or sets with many near-identical genomes. It splits
the genomes into chunks of at most N, dereplicates each chunk (in parallel,
`--num-processes`) with the stage-1 thresholds (`--pre-primary-ani`,
`--pre-secondary-ani`, which default to the final ones), and then
dereplicates the union of the chunk representatives with the final
thresholds. Every genome is then assigned to the final representative of its
chunk representative.

**The final pass.** The chunk size bounds the leaf chunks, not the final
pass. A chunk can only remove genomes that have a close relative in the same
chunk, so for a diverse set at a strict threshold the union is most of the
input, and the final pass compares nearly all genomes with each other.

**Further levels and the probe.** When the union is larger than two chunks
and the level removed at least a tenth of its input or a chunk's worth (N),
whichever is larger, one chunk of the next level is run first as a probe: the
first N representatives taken in turn from each chunk of the previous level,
so that genomes from different chunks meet. A level over the union costs the
tool's per-genome work (sketching, N50) for every representative plus about
N/U of the final pass's pairwise work, where U is the union, and is expected
to pay when it removes more than about N/U of its genomes.

- When the probe removes more than that fraction, the rest of the level runs
  (the probe is its first chunk) and its union is treated the same way.
- When it removes less, the union goes to the final pass and the probe's
  result is dropped.
- When a level removes less than the shrink rule asks, the union goes to the
  final pass without a probe.

The probe, the levels after the first and the final pass compare
representatives at the final thresholds, not the `--pre-*` ones.

**The log.** The log names the size of the final pass and what to expect of
it, from the most recently timed level:

```text
Level 0: 6888 genomes -> union of 3807 representatives (3081 removed) in 8.0 min
Level 0 removed 3081 (>= 689); probing level 1 with one chunk of 500 of the 3807 representatives, interleaved by chunk of origin
Probe: 24 of 500 removed (4.8%) in 1.2 min; a level over 3807 is expected to pay above 13.1%
Level 0: a further level is expected to cost more than it saves; the probe's result is dropped
Final pass over 3807 representatives (chunk size 500 bounds the leaf chunks, not the merge)
One call over 3807 genomes with 48 threads: about 4.2 times the pairwise work of level 0, which took 8.0 min with 12 worker(s); no progress is printed while it runs
```

A warning follows when the final pass is more than four times the chunk size.
The pairwise factor is the union squared over the level's input times the
chunk size. The per-genome work is paid once more for every genome of the
union, and the tool's threads are shared by one call rather than by parallel
chunks, so the pass can take longer than the factor alone suggests. The
example is 6888 Bacillaceae genomes with `skder --secondary-ani 0.9995 -s 500`
(the probe line is illustrative): 3807 representatives after level 0, and the
earlier levels of that run removed 5, 1.4 and 0.5 percent each, so the probe
is expected to send the union to the final pass.

**Making the final pass smaller or faster.**

- **Keep the default skani screen with skDER.** For a cutoff at or above 99
  percent the adapter passes skani a screen of the cutoff minus one, so pairs
  whose k-mer estimate is below that are not aligned; skDER's own default is
  the cutoff minus ten, which at 99.95 aligns every pair above about 90
  percent. `--tool-arg screen=PERCENT` sets another value and `screen=skder`
  restores skDER's own; a value above the cutoff is refused. The screen is
  an option of `skani triangle`, which the greedy and dynamic modes run;
  `low_mem_greedy` does not read it.
- **Raise `-s`** so that each chunk holds more close relatives. Chunking pays
  when the set has many near-identical genomes; for a set where most genomes
  are their own representative, one pass without `-s` costs about the same.
- **Set `--pre-secondary-ani` below `--secondary-ani`** (for example 0.999
  against 0.9995) so the chunks remove more. The chunk representatives are
  then compared at the final threshold, so clusters can be coarser than in a
  single pass at the final threshold.

The Nextflow scatter (`--derep_process_size`) runs one chunk round and one
merge task over the whole union, with the same semantics as the first level
and the final pass here (see
[Scatter-gather dereplication](#scatter-gather-dereplication)).

### Inspecting a dereplication

Three commands read a dereplicated working directory without rerunning the
dereplicator, and a fourth inspects the genomes before dereplication:

| Command | What it does |
|---------|--------------|
| `repgenr cluster-summary` | Regenerates `derep/cluster_summary.tsv`, one row per representative (see [output.md](output.md)). |
| `repgenr derep-unpack` | Lays the clusters out as one directory per representative with its members inside. |
| `repgenr derep-stock` | Stores, restores, lists and deletes dereplication results under `derep/stock/`. |
| `repgenr glance` | Compares all genomes in `genomes/` and plots the pairwise values; it does not need a dereplication. |

The commands read `derep/clusters.tsv` and the manifest. They do not check
that the last `dereplicate` finished, so run `repgenr status -wd WD` first: a
failed or refused `dereplicate` shows as `[interrupted]` while `derep/` still
holds the previous result.

#### Which command answers which question

| Question | Where to look |
|----------|---------------|
| How many clusters, and how large is each? | `derep/cluster_summary.tsv` (`repgenr cluster-summary` rebuilds it); `n_members` does not count the representative. |
| Which genomes are in the cluster of representative X? | Rows of `derep/clusters.tsv` whose first column is X, or the directory `derep/unpacked/<X without its extension>/` after `repgenr derep-unpack`. |
| What happened to one genome? | `derep/genome_status.tsv`: `representative`, `contained` or `fail_qc`. |
| Do the genomes fall into clear groups before I pick thresholds? | `repgenr glance`, see below. |
| How do two dereplications differ? | Store each with `derep-stock --action pack`, then compare the stored files, see below. |

#### Finding the members of a cluster

`clusters.tsv` has one row per genome, with the representative in the first
column and the member in the second. The representative also lists itself.

```bash
# members of the third-largest cluster (summary rows are sorted by size)
REP=$(sed -n 4p $WD/derep/cluster_summary.tsv | cut -f1)
awk -F'\t' -v r="$REP" '$1 == r {print $2}' $WD/derep/clusters.tsv
```

#### Unpacking clusters into directories

`repgenr derep-unpack` gives the same grouping as directories, which is
convenient for a tool that takes a directory. The directory name is the
representative's filename without its extension. `--no-representant` leaves
the representative out. A member missing from `genomes/` is left out with a
warning that names it.

The files are hard links to the genomes where the file system allows it, so
editing one in place also changes the genome under `genomes/` (and, after a
default symlinking `ingest`, the original file). Elsewhere, for example on
exFAT or across volumes, they are copies, which for 1000 genomes took minutes
rather than under a second.

The directories describe the dereplication they were made from. After a new
`dereplicate` or a `derep-stock --action unpack`, run `derep-unpack` again (it
reruns because `clusters.tsv` changed, and `doctor` warns until then).

#### Storing and restoring a dereplication

`repgenr derep-stock --action pack --name <run>` stores the current clusters,
statuses, representatives and the completed `dereplicate` record (tool,
parameters, tool versions, as `record.json`) under `derep/stock/<run>`. A run
already stored under that name is replaced, with a warning. `--action list`
prints the stored run names on stdout, one per line in name order, and
`--action delete` removes one. Deleting a run that is not stored exits 3 and
lists the stored runs.

`derep-stock --action unpack --name <run>` makes a stored run the current one
again. It refreshes the manifest and re-stamps the `dereplicate` record.

- Unpack rebuilds `cluster_summary.tsv` from the restored clusters and the
  current manifest quality, and removes a live `genome_status.tsv` the stored
  run lacks. The N50 of a scored genome whose file is no longer in `genomes/`
  cannot be read, so its keeper score in the summary has no N50 term and
  `rep_n50` is blank; unpack and `cluster-summary` name such genomes in a
  warning.
- Unpack replaces `derep/representatives/` as a whole, so other files placed
  there are removed, and it takes the representatives by name from
  `genomes/`.
- The re-stamped record takes the tool, parameters and tool versions from the
  stored `record.json` and adds `stock: <run>`, so `status` and `versions`
  name the tool that produced the stored run; its completion time is that of
  the unpack.
- A run stored without `record.json` (packed before this file was written,
  or packed while no completed `dereplicate` record existed, which pack
  reports with a warning) keeps the tool and parameters of the `dereplicate`
  record current at unpack time instead. In that case a record left
  incomplete by an interrupted `dereplicate` run carries nothing over (no
  tool), and one left by an interrupted unpack keeps what that unpack was
  carrying.

The restored `dereplicate` record holds no resume fingerprint, since unpack
cannot reproduce the conditions of the original run. The next `dereplicate`,
also inside `repgenr run`, therefore recomputes the dereplication; run
`phylo` directly to build on the restored set. `status` and `doctor` note
such a record.

#### Reading `glance`

`repgenr glance` compares all genomes in `genomes/` (at least two) and writes
a dendrogram and two plots. Use it before dereplication to see whether the
genomes fall into clear groups; it does not need a dereplication and sets no
thresholds. `--tool` chooses the comparison tool:

- `auto` (the default) uses dRep when it can run, on the `PATH` or through
  the container backend, and sourmash otherwise. The log names the tool it
  picked, and the stage record holds that tool, not `auto`. With neither
  available, glance exits 4 and names the tools that can compare.
- `drep` runs `dRep compare` (Mash, primary clustering only) and plots Mash
  ANI. The dendrogram is dRep's own.
- `sourmash` compares the genomes with the parameters `dereplicate --tool
  sourmash` uses (k=31, scaled=1000), reading the genome sketches of the
  working directory and writing the missing ones first, runs `sourmash
  compare`, and converts the values to the same ANI estimate the
  dereplication threshold applies to. The dendrogram is average-linkage
  clustering on 1 - ANI. The full N x N matrix is held in memory, so sets
  above 5000 genomes are refused.

The plot axes name the measure, Mash ANI or ANI; the file names stay the same
for both tools. The two estimates differ slightly, so compare plots made with
the same tool.

The histogram shows how the pairwise values are spread. When they fall into
separate groups with an empty gap between them, a threshold placed in the gap
separates them. On the synthetic set `clonal_50_clustered` (three groups of
20, 15 and 15 genomes) the between-group values lie at about 0.95 to 0.96,
the within-group values at 0.995 or higher, and nothing lies between 0.96 and
0.995 (Mash ANI from dRep and the sourmash ANI estimate agree on this). If the
values are spread evenly across the range, the threshold decides the cluster
sizes, so check `cluster_summary.tsv` afterwards.

With `--tool sourmash` and `--keep-files`, `glance_wd/pairwise_ani.csv` holds
the value of every pair and `glance_wd/dendrogram_leaves.txt` the
dendrogram's leaf order. dRep writes into its cache with names beginning
`._` on exFAT volumes and fails there; use an APFS or ext4 working directory
(see [verification.md](verification.md)).

#### Comparing two dereplications

Pack each result under its own name, then compare the stored files. Check
with `status` that `dereplicate` finished before each pack; each stored run's
`record.json` then names the tool and parameters that produced it.

```bash
repgenr dereplicate -wd $WD --tool sourmash
repgenr derep-stock -wd $WD --action pack --name sourmash
repgenr dereplicate -wd $WD --tool galah
repgenr derep-stock -wd $WD --action pack --name galah
for n in sourmash galah; do
  echo "$n: $(tail -n +2 $WD/derep/stock/$n/clusters.tsv | cut -f1 | sort | uniq -c | awk '{print $1}' | tr '\n' ' ')"
done
diff <(ls $WD/derep/stock/sourmash/representatives) <(ls $WD/derep/stock/galah/representatives)
```

Compare the cluster sizes first, then the memberships. Two dereplicators can
cluster the genomes identically and still keep different genomes as
representatives (see the dereplicator table in
[choosing-tools.md](choosing-tools.md)), so a different representative list
alone does not mean a different partition. To compare partitions, label every
genome with the smallest member name of its cluster and diff the labels:

```bash
part() { awk -F'\t' 'NR==FNR { if (NR>1 && (!($1 in m) || $2<m[$1])) m[$1]=$2; next }
                     FNR>1 { print $2 "\t" m[$1] }' "$1" "$1" | sort; }
diff <(part $WD/derep/stock/sourmash/clusters.tsv) \
     <(part $WD/derep/stock/galah/clusters.tsv) && echo "same partition"
```

To make one of the stored runs current again, unpack it (see
[Storing and restoring a dereplication](#storing-and-restoring-a-dereplication)).

## Phylogeny

`repgenr phylo` builds a tree from the representatives and the outgroup and
writes it to `tree/tree.nwk`. It runs after `dereplicate` and before
`tree2tax`. The data-channel step `phylo-build` does the same outside a
working directory. This part covers the choice of route, SNP typing and
recombination masking, annotation and the core-gene route, rooting, and the
reuse of an alignment.

### Choosing a route

A tree is built either from an MSA (multiple sequence alignment) or directly
from the genomes. The MSA comes from one of three sources:

| Route | Option | Section |
|-------|--------|---------|
| Whole-genome alignment | `--aligner` (for example `progressivemauve`) | the bacterial example in [Bacteria from GTDB](#bacteria-from-gtdb) |
| Core-SNP alignment | `--msa-source snptype` | [SNP typing and masking](#snp-typing-and-masking) |
| Concatenated core-gene alignment | `--msa-source coregene` | [Core-gene phylogeny](#core-gene-phylogeny) |
| Alignment-free tree | a tree builder that works from the genomes, such as `--treebuilder mashtree` | below |

[choosing-tools.md](choosing-tools.md#phylogeny-routes) explains when to use
which. Examples of the alternatives to the whole-genome alignment:

```bash
# Scalable dereplication then an alignment-free tree
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree

# SNP-based phylogeny (core-SNP alignment as the MSA source); phylo types the
# genomes itself under tree/msa/, so no separate snptype run is needed
repgenr phylo -wd $WD --msa-source snptype --treebuilder iqtree --mask gubbins

# Core-gene phylogeny: annotate first, then build from the core genes
repgenr annotate -wd $WD
repgenr phylo -wd $WD --msa-source coregene --treebuilder iqtree
```

### Genome count, leaf check and record names

`phylo` and `phylo-build` need at least three ingroup genomes (the outgroup is
not counted) and exit 3 with a message before any tool runs when the set is
smaller.

After the tree is built, its leaves are compared with the input genomes
(ingroup and outgroup). A tree builder can drop a genome it considers
degenerate and still exit 0, as mashtree does; a missing or unexpected leaf
then exits 3 with the names and the builder. The tree is kept in `tree/` for
inspection, and `phylo` is not recorded as completed. Leaf names are compared
without a FASTA extension or ParSNP's `.ref` reference marker, and with
characters other than letters, digits, `_` and `-` read as `_`, since some
tools rewrite them. This matching is looser than `tree2tax`, which uses leaf
names as written, so a tree that passes the check can still leave leaves that
`tree2tax` does not map to a genome.

Records in every alignment (the SNP typers' and the aligners'
`align/msa.fasta`), the leaves of every tree including the sourmash tree, and
the outgroup name given to the masker are the genome file name without its
FASTA suffix and `.gz` (`x.fasta.gz` gives `x`), as `clusters.tsv` and
tree2tax name genomes. `parsnp` names its records by file name and marks the
reference with `.ref`; the typer renames them to the genome names, as the
other typers write them.

### SNP typing and masking

A SNP typer produces a core-SNP alignment, the variable sites among the
genomes. It suits closely related isolates; recombination masking with
Gubbins expects isolates of one species. The `repgenr snptype`
command (and `phylo-build --msa-source snptype`, and `phylo --msa-source
snptype`) calls the typer.

Four typers are built in: `simple` (minimap2, samtools and bcftools, the
default), `snippy`, `parsnp` and `ska2`. The first three map every genome to
one reference and also write the whole-genome alignment a masker needs.
Recombination masking (`--mask gubbins`) runs on the typer's whole-genome
alignment and replaces the core-SNP alignment with Gubbins' filtered
polymorphic sites. Typers that only emit variable sites cannot be masked.

#### Invariant sites and branch lengths

A tree from a variable-site-only alignment overestimates branch lengths
unless the builder is told how many invariant sites the alignment lost. When
the typer writes a whole-genome alignment (`simple`, `snippy`, `parsnp`), the
typing pass counts its constant A, C, G and T columns into
`alignment_stats.json` beside `core_snp.fasta`, and `phylo` passes them on:

- IQ-TREE receives them as `-fconst`, RAxML-NG as `+ASC_STAM{...}` on its
  model. Both are first given the alignment without any invariant column that
  masking left, and those columns join the counts.
- `--aligner-arg asc=lewis` uses the Lewis correction (`+ASC`, `+ASC_LEWIS`)
  instead of the counts, and a `model=` that already carries an ASC term is
  passed as it is, without the counts.
- FastTree applies no correction and says so.

The counts apply only to the alignment they were taken for (the sidecar names
it and holds its digest), so a tree built from `full_alignment.fasta` is not
corrected twice. In the Nextflow split mode the alignment task publishes the
sidecar and the tree task reads it. The counts are taken before masking, so
after `--mask` they are a slight undercount.

`ska2` writes no whole-genome alignment, so its tree has no correction and
`phylo` warns. Compare topologies and supports rather than lengths there, or
build the tree from a whole-genome alignment
(`phylo-build --msa snp/full_alignment.fasta`).

#### Soft core (`--core-fraction`)

The variable-site alignment of the mapping typers keeps every column in which
two bases occur, however many genomes lack a base there, so a position
covered in two genomes out of seventy is in it. `--core-fraction F` on
`snptype`, `phylo`, `phylo-build` and `run` keeps only the columns where at
least the fraction F of the genomes have a called base (a soft core; 0.95
is the usual choice). It is applied after masking, and the constant sites are counted
under the same rule.

It is unset by default: every variable column is kept, and the log reports
how many columns a fraction of 0.95 would keep, so the effect can be seen
before the rule is chosen. The default may become 0.95 in a later release,
after the reported counts have been compared on real data sets. The value is
part of the alignment reuse key and of `alignment_stats.json`.

Three cases differ from the description above:

- `snippy`'s core already holds only the sites called in every genome, so the
  filter changes nothing there, while the constant sites are still counted at
  F.
- `ska2` already keeps the split k-mers present in a fraction `min_freq` of
  the samples (0.9 by default), so `--core-fraction` is a second threshold on
  top of it.
- After `--mask`, the bases the masker replaced count as missing, so a column
  inside a recombination block shared by more than 1 - F of the genomes is
  dropped, while the constant sites come from the unmasked alignment.

`--core-fraction` applies only with `--msa-source snptype` (or, on `run`, to
the standalone `snptype` stage under `--with-snptype`). An aligner alignment
is left as it is and the option is refused there, and a tree builder that
works from the genomes warns that it has no effect.

#### Recombination masking with Gubbins

Gubbins expects isolates of one species. The masker estimates how much of the
alignment is variable, warns above 10%, and repeats the figure if Gubbins
fails: on a genus-level set its scan can die with a bus error (see
[verification.md](verification.md)).

**Refusals.** Masking is refused (exit 2) when more than 12% of the columns
are variable, near where Gubbins has failed (13% in
[verification.md](verification.md)), or when the genome names (canonical
`Family_Genus_species_Accession` names only, read where the last token is an
accession) span more than one species. The species check runs before SNP
typing. `--tool-arg gubbins_divergence_limit=FRACTION` moves the first limit
and `--tool-arg gubbins_force=true` runs it anyway, with a warning.
Alignments shorter than 1000 columns are not judged.

**Gubbins' tree builders.** Gubbins builds a tree in every iteration, by
default with RAxML, and with more than one thread it needs a multi-threaded
RAxML build (`raxmlHPC-PTHREADS*`). Some conda builds ship only the
single-threaded binary, and Gubbins then exits before its first iteration.
When repgenr runs Gubbins natively and finds no such build, it switches to
IQ-TREE (or to one thread when IQ-TREE is missing too) and says so in the
log. `--tool-arg gubbins_tree_builder=raxmlng`,
`--tool-arg gubbins_first_tree_builder=rapidnj` and
`--tool-arg gubbins_args="--min-snps 5"` pass the choice, the first-iteration
builder and any further `run_gubbins.py` arguments through.

Unless a first-iteration builder or model is given, repgenr asks for
`rapidnj` with the JC model when rapidnj is available (it ships with
Gubbins). With IQ-TREE as Gubbins' tree builder, tree building took 76 of 87
minutes on a genus-level set ([verification.md](verification.md)), and a
rapidnj first tree makes the first iteration cheap while later iterations
still use the resolved builder. Gubbins' own final tree, its
recombination-aware maximum-likelihood tree, is kept as `masker_tree.nwk`
beside the masked alignment for comparison with the tree repgenr builds from
the filtered sites.

**Gapped taxa.** Gubbins leaves taxa with more than 25% gaps or N
(`--filter-percentage`) out of its recombination analysis but still writes
them to its outputs, unmasked. Since the `simple` typer writes N where a
genome does not align, repgenr passes `--filter-percentage 100`. A value given
in `gubbins_args` takes precedence; a taxon above it is then refused (exit 3)
before Gubbins runs, with its name.

#### The `simple` typer

**Threads and scratch.** The `simple` typer maps each genome independently,
so `--threads` buys concurrent genomes first and threads inside one genome's
chain only when there are more threads than genomes. Its per-genome
intermediates are written compressed and removed as soon as that genome's
consensus has been read, so scratch stays at a few hundred megabytes whatever
the genome count. A genome whose chain fails keeps its intermediates for
inspection. Under a container backend each genome's chain of tools runs in a
single container, so a genome costs one engine start rather than eight.

**Missing sequence.** The typer builds each genome's consensus by applying
its SNP calls to the reference, and then sets to N every reference position
where none of the genome's primary or supplementary alignments places a base:
positions outside the alignments, and positions within a deletion (CIGAR D
or N operations), as read from the minimap2 SAM. Sequence a genome lacks is
therefore missing data, not the reference base.

A column of the core-SNP alignment (`snp/core_snp.fasta`, or
`tree/msa/core_snp.fasta` from the typing pass of `phylo`) is kept only when
at least two of A, C, G and T occur in it. N, other ambiguity codes and gaps
do not make a column variable, and they are written as N.

**Distances.** The distance matrix (`snp_distance_matrix.tsv` beside it)
counts, for each pair, the differing sites among those where both genomes
have a base, so two genomes are not separated by a region one of them lacks.
A pair that shares no such site has the distance `NA`. On the 50-genome test
set a copy of one genome with 500 kb removed now differs from that genome at
0 sites; before masking it differed at 20539, all but one of them within the
removed region. Distances between genomes of different gene content are
computed over different numbers of sites, so they are not directly
comparable as proportions.

**Coverage.** The log has one line per genome with the fraction of the
reference's bases that genome covers, and a summary line (minimum, median,
maximum). A genome below 50% is warned about. A genome with no base at any
core SNP site, usually one that did not align to the reference at all, is
refused (exit 3), since tree builders refuse a sequence of N only; the
message names it.

**minimap2 preset.** The typer maps assemblies with the minimap2 preset
`asm20`, which suits assembly-to-reference alignment up to several percent
divergence. Three genomes of the 50-genome test set have equal length
(2 Mb), so their true substitution counts are position-wise differences.
With `asm20` the typer's distances were within one site of those counts; with
minimap2's default settings, used before, they differed by 8 to 111 sites. A
genome more divergent than `asm20` tolerates, such as an outgroup from
another species, aligns over less of the reference and is correspondingly
more N; check the coverage lines in the log. `--tool-arg preset=` accepts
`asm5`, `asm10`, `asm20`, `map-ont`, `map-pb`, `map-hifi`, `sr` and `none`
(minimap2's default settings, no `-x`), and refuses any other value (exit 2).
`asm5` suits near-identical genomes, and `none` aligns more divergent
sequence at the cost of the accuracy above.

#### Split k-mers: ska2

`--tool ska2` (split k-mer analysis) is reference-free: every genome is an
ordinary sample, so no assembly's private errors bias the SNP distances, and
the alphabetical-first-genome reference default does not apply. It emits a
variable-site alignment only, so it is not compatible with `--mask`. Tune it
with `--tool-arg ksize=31` and `--tool-arg min_freq=0.9` (the fraction of
samples a split k-mer must occur in).

On the 60 Wolbachia species representatives from GTDB, a genus-level set, it
returned 708 variable sites in about two minutes. It is designed for clonal
and outbreak sets, where the shared k-mer space is far larger, and the stage
says so when the genome names span more than one species.

`--tool-arg lo=true` replaces `ska align` with `ska lo`, which recovers
variants in more diverged samples by local assembly of the split k-mer graph
(Horsfield et al. 2025). The indels it finds are published as
`snp/variants.vcf`, `lo_missing` (0.1) is its limit on missing data per site,
and `min_freq` is not read in this mode.

`--treebuilder ska2` builds a tree without an alignment from the same split
k-mers. `ska distance` counts the SNPs among the k-mers each pair shares
(`metric=snps`, the default, branch lengths in SNPs) or the fraction of
mismatching k-mers (`metric=mismatch`), and the matrix is joined with rapidNJ
when installed. Pairs sharing few k-mers have their SNPs undercounted, so the
same divergence warning applies (the outgroup is left out of the count), and
a pair that shares none is refused, since ska writes it as distance 0.
`min_freq` here is `ska distance`'s own (0 by default). `auto` never picks
this builder; name it.

### Annotation

`repgenr annotate -wd WD` annotates the genome set once, for the core-gene
phylogeny route and any later use of the gene calls. It annotates the
representatives (or every genome with `--all-genomes`) and the outgroup, one
directory per record under `annotations/` with the GFF3 (sequences
appended), the proteins and the gene nucleotides, and an index
`annotations.tsv`.

```bash
repgenr annotate -wd $WD --bakta-db $BAKTA_DB
```

**Annotators.** The annotator is `bakta` by default, which needs its
database. Download it once with `bakta_db download --output DIR --type light`
(a 1.4 GB download, about 4 GB unpacked; `--type full` is 30 GB to download
and about 84 GB unpacked, and annotates CDS more specifically) and pass
`--bakta-db DIR` or set `BAKTA_DB`. The database is part of the stage's
resume inputs through its `version.json`, so replacing it re-annotates.
`--tool prokka` needs no database. Tuning goes through `--tool-arg`
(`skip_crispr=true`, `genus=Francisella`, `translation_table=11` for bakta;
`kingdom=Bacteria`, `gcode=11` for prokka).

**Run time and parallel genomes.** bakta takes about 5 to 15 minutes per
genome. An annotator uses many threads in only some of its steps, so the stage
annotates several genomes at once and splits `-t` between them:
`-p/--num-processes` sets how many, and the default (0) takes about one genome
per four threads, capped by the usable CPUs (`-t 16` gives four genomes with
four threads each). Each running genome is a separate annotator process; allow
a few GB of memory for each with bakta. A genus-level set is still a cluster
job, where the Nextflow route runs one task per genome. `phylo-build`, when it
annotates its own inputs, divides its `-t` in the same way. `repgenr run
--msa-source coregene` passes its `-p` to the annotate stage.

**Resume.** A later run annotates only the genomes that are new or
whose genome, annotator, tool version, database version or tuning changed,
and reuses the rest. The index is written as each genome finishes, so an
interrupted or failed run keeps the genomes that completed; after a failure,
the genomes already running finish, no new ones start, and the first error is
reported. A record whose genome left the set is kept (annotation
is expensive and the genome may return, for instance after another
dereplication); only a record whose genome has left the working directory is
removed. `--force` re-annotates everything.

The data-channel step `annotate-build --genome FILE -o DIR` annotates one
genome.

### Core-gene phylogeny

`phylo --msa-source coregene` builds the tree from a concatenated alignment
of the core genes, the genes present in (nearly) all genomes of the set,
instead of a reference-based SNP alignment or a whole-genome alignment.
[choosing-tools.md](choosing-tools.md#phylogeny-routes) compares it with the
other routes.

```bash
repgenr annotate -wd $WD
repgenr phylo -wd $WD --msa-source coregene --treebuilder iqtree
```

On a working directory, run `annotate` first, as above.
`repgenr run --msa-source coregene` inserts the `annotate` stage before
`phylo` (`--annotator`, `--bakta-db`). The data-channel step
`phylo-build --msa-source coregene` reads `--annotations-dir` or annotates its
inputs itself.

The genomes are annotated (the [annotate](#annotation) stage), a core-gene
builder (`--coregene`) clusters the genes, and the core genes are aligned and
concatenated. Four builders are available: `panaroo` (the default),
`ppanggolin`, `scarap` and `orthofinder`, described below. The alignment
keeps its invariant sites, so the tree builder applies no
ascertainment-bias correction and the branch lengths are substitutions per
site of the core.

**The outgroup.** The outgroup is annotated and included. Being another
species, it lowers the core, so on fewer than about twenty genomes lower
`--core-threshold` or leave the outgroup out with
`--no-outgroup --root midpoint`.

**Reuse and outputs.** The alignment is reused across tree-builder changes
like the other sources, and rebuilt when an annotation changed (the reuse
stamp holds the index digest). Outputs:
`tree/msa/core_gene_alignment.fasta`, the Panaroo tables under
`tree/msa/pangenome/` and `tree/msa/core_gene_stats.json` (see
[output.md](output.md)).

#### Panaroo and its clustering identity

Panaroo (`--coregene panaroo`, the default) clusters the genes into a
pangenome, cleans annotation errors from the graph (`--clean-mode strict`,
the default, `moderate` or `sensitive`), aligns every gene present in at
least `--core-threshold` of the genomes (0.95) with mafft, concatenates them
and removes the genes its entropy filter flags (`--aligner-arg
filtered=false` keeps them).

The parameter to reach for first is the clustering identity. Panaroo groups
genes in a first pass at `--aligner-arg threshold=` identity, 0.98 by
default, so on a divergent set most orthologues of a gene end up in separate
clusters and never reach the core. Lowering it is what recovers them:
measured on three Wolbachia genomes of one supergroup, the genes present in
all three went from 52 at 0.98 to 226 at 0.95, 638 at 0.90 and 796 at 0.80
(see [verification.md](verification.md)). `family_threshold`, the later
protein-family step, did not compensate at 0.5 or 0.3. Lowering
`--core-threshold`, the fraction of genomes a gene must be in, recovers less
and gives a softer core: over 24 Wolbachia genomes at the default identity it
gave 0 genes at 0.95 occupancy against 95 at 0.50.

How much a given setting yields still depends on the set. The same 0.98 gave
1054 genes for two of those genomes, 52 for three and none at all for four
spanning different supergroups, so a wider set needs a lower identity rather
than the same one. RIBAP (Genome Biology 2024) reports the same for other
genera and cautions that too low a value starts to group genes that are not
orthologues, so lower it in steps and check the result.

The panaroo builder warns when the core comes out under 100 genes, unless you
have set `threshold` yourself. `core_gene_stats.json` carries the count and
the identity used, and the count is the number to look at before trusting the
tree. Within an outbreak a SNP tree still draws on far more positions and
costs no annotation pass. Panaroo's declared limit of 500 genomes is not
measured here, and neither is its runtime; the clustering and the per-gene
alignments are the slow parts, after the annotation.

#### PPanGGOLiN

`--coregene ppanggolin` is the alternative for divergent or large sets.
PPanGGOLiN clusters the proteins with MMseqs2 at 80% identity over 80% of
their length (`--aligner-arg identity=` and `coverage=`), partitions the
pangenome graph, and concatenates the DNA alignments of the selected
families.

The selection defaults to the soft core at `--core-threshold`, so the gene
set means the same as with Panaroo. `partition=core` takes the strict core
and `partition=persistent` the tool's statistical partition, which does not
shrink as genomes are added but needs a set of some size (PPanGGOLiN warns
below 15 genomes, and on three it placed every family in it).

PPanGGOLiN aligns only genomes with exactly one copy of a family, so a
multi-copy genome is a gap for it. Families multi-copy in at least 5% of the
genomes (`dup_margin=`) are left out by default (`single_copy=false` keeps
them). A genome that shares no selected family is dropped from the alignment
with a warning, and an outgroup in that state is an error. `--clean-mode` and
`filtered=` are Panaroo's and have no effect here; they are left out of the
resume fingerprint and the alignment stamp, so changing them reruns nothing.

On 24 complete Wolbachia genomes, where Panaroo at its defaults found no gene
in 95% of them, PPanGGOLiN selected 398 single-copy families and the whole
route, IQ-TREE included, took under four minutes (see
[verification.md](verification.md)). Its declared limit of 5000 genomes is
not measured here.

#### SCARAP

`--coregene scarap` infers the core directly, without the whole pangenome.
SCARAP builds orthogroups on a set of seed genomes (`--aligner-arg seeds=`,
100), turns them into MMseqs2 profiles and searches every genome against
them; a gene is core when it is present in a single copy in at least
`--core-threshold` of the genomes. It reads the protein and gene files the
annotate stage writes beside each GFF3, not the GFF3 itself, and aligns each
orthogroup with mafft.

The orthogroups come from profile searches rather than a fixed identity cut,
so there is no identity to lower. `method=` chooses the orthogroup inference
(`FH`, the default) and `max_core_genes=` caps the count. The seed genomes are
drawn at random with no fixed random seed, so on a set larger than `seeds`
two runs can select slightly different cores; the final selection is checked
over all the genomes, not only the seeds. SCARAP aligns with as many
processes as the machine has CPUs, whatever `-t` says.

On the same 24 Wolbachia genomes it found 577 single-copy core genes, and its
IQ-TREE tree shared 20 of 21 splits with the PPanGGOLiN one. It took about
six and a half minutes on this machine, against three for PPanGGOLiN. Its
declared limit of 10000 genomes is not measured here; the authors ran it on
over 31000.

`method=O-D` and `method=O-B` build SCARAP's seed pangenome with OrthoFinder
instead (DIAMOND or BLAST search). They run natively only, from the
environment `envs/scarap.yml` builds: the SCARAP image has no OrthoFinder,
and SCARAP 1.0.2 reads OrthoFinder's result with a pandas call that pandas 2
removed, so that environment keeps pandas older than 2. On the 24 genomes
they found 565 (`O-D`) and 568 (`O-B`) core genes, against 576 for the
default method in the same environment, and took several times longer;
through RepGenR2 `O-D` gave the same tree as the default method.

#### OrthoFinder

`--coregene orthofinder` runs OrthoFinder on its own. It compares every
proteome against every other with DIAMOND (`--aligner-arg search=`, also
`diamond_ultra_sens`, `blastp` or `mmseqs`), clusters the hits into
orthogroups with MCL (`inflation=`, 1.2) and stops there (`-og`); the gene
trees it would build next are not needed.

The adapter keeps the orthogroups present in a single copy in at least
`--core-threshold` of the genomes, aligns each one's proteins with mafft and
back-translates the alignment to codons from the `.ffn` files, so the result
is a DNA alignment like the other builders'. A gene whose length does not
match its protein is a gap, and the log counts them; selenocysteine and
pyrrolysine, which mafft does not accept, are aligned as X. With
`search=mmseqs` OrthoFinder uses fixed file names under `/tmp`, so two such
runs on one host can collide.

On the 24 genomes it kept 585 orthogroups, and its tree shared 20 of 21
splits with the PPanGGOLiN tree and 19 with the SCARAP one. Its core-gene
step took about six minutes, most of it the search, close to SCARAP's. The
search is all against all, so its cost grows roughly with the square of the
genome count; the declared limit of 200 genomes is not measured.

### Rooting the tree

`phylo` (and `phylo-build`, and `run`) roots `tree/tree.nwk` after any tree
builder has finished. The rooting is chosen with `--root`:

- `outgroup` roots on the branch leading to the outgroup, which is split in
  two halves, so the root has two children: the outgroup and the whole
  ingroup. This is the default when an outgroup is staged and
  `--no-outgroup` is not given.
- `midpoint` roots at the midpoint of the longest leaf-to-leaf path. It needs
  branch lengths and does not use the outgroup, so it is the choice for a
  tree built with `--no-outgroup`.
- `none` keeps the tree as the builder wrote it. This is the default without
  an outgroup.

`--root outgroup` without a staged outgroup, or together with
`--no-outgroup`, exits 2 before any tool runs.

Builders that are given the outgroup (`iqtree -o`, `raxmlng --outgroup`) place
it at a root with three children; such a tree is rerooted like any other,
with one INFO line. A tree whose root already has the outgroup as one of two
children is left as written. Branch lengths and support values are kept: each
support value stays with the split it was computed for. A midpoint that falls
on an internal branch splits that branch in two, and its support value is
written on both new root branches, since both describe the same split. When
the tree is rerooted, the builder's own output is kept as
`tree/tree.unrooted.nwk`. The stage record holds the rooting applied (`root`)
and the outgroup leaf (`outgroup`).

**Rooting in tree2tax.** `tree2tax` and `tree2tax-relations` root the tree on
the outgroup when it is not rooted there yet. On a tree from
`phylo --root outgroup` this changes nothing, and the relations are the same
as when `tree2tax` rooted the builder's output itself. They exit 3 when the
outgroup is not a leaf of the tree.

- After `phylo --root midpoint`, `tree2tax` keeps the midpoint root (it reads
  the `root` that phylo recorded); `tree2tax-relations` does so when given
  `--keep-root`. The root is named by `--root-name` in either case, and
  `--remove-outgroup` still removes the outgroup leaf. It is pruned before the
  collapse and naming steps, so no clade is named for a leaf set that held it,
  and a node left with a single child where it sat is merged into its parent.
- After `phylo --no-outgroup` with `--root none`, `tree2tax` leaves the tree
  unrooted (its relations start at the builder's root); `tree2tax-relations`
  does so when given `--no-outgroup`, which ignores `--outgroup-dir` and
  `--outgroup-accession`.

The Nextflow pipeline passes `--no-outgroup` to `tree2tax-relations` when
`phylo_args` contains it, and `--keep-root` when `phylo_args` contains
`--root midpoint`.

### Reusing an alignment across tree builders

`phylo` stamps the alignment it builds (`align/msa_source.json` or
`tree/msa/msa_source.json`) with what produced it: the source and its
settings, the genome set, and the alignment's own digest. A later `phylo` run
that changes only the tree builder, the bootstrap or the thread count reuses
that alignment instead of aligning or SNP-calling again, and says so in the
log. Anything the alignment depends on, such as the aligner, the SNP typer,
`--reference`, `--mask` or the genome set, rebuilds it, as does `--force` or a
change to the alignment file itself.

The same split is available to the stateless step: `phylo-build --msa-only`
builds the alignment and writes `msa.fasta` without a tree, and
`phylo-build --msa <file>` builds a tree from an alignment an earlier call
produced. The Nextflow layer uses these to run the alignment and the tree as
separate tasks. `--msa-only` leaves the SNP typing pass's `tree/msa/`
directory (and `scratch/`) in the step's output directory beside
`msa.fasta`; they are left in place.

## tree2tax

`repgenr tree2tax` turns the rooted tree into a FlexTaxD-compatible taxonomy.
It is the last stage and writes `tree2tax.tsv` (the relations) and
`genomes_map.tsv`, with `tree2tax_nodes.tsv` beside them (see
[output.md](output.md)). By default (`--include-dereplicated`) the genomes
that dereplication left out are listed under their representative. The
`tree2tax-relations` step does the same outside a working directory and is
the one the Nextflow module runs.

```bash
repgenr tree2tax -wd $WD --include-dereplicated
```

How `tree2tax` roots the tree, and how it treats a midpoint root or a tree
without an outgroup, is described under [Rooting the tree](#rooting-the-tree).

### Collapsing weak splits in tree2tax

`repgenr tree2tax` (and the `tree2tax-relations` step) can merge weakly
supported or near-zero-length internal nodes into their parents before
naming them, so a split the data do not support does not become its own
FlexTaxD node. All thresholds are off by default.

**Branch length and support.** Either criterion suffices.

- `--collapse-length L` merges a node whose branch is shorter than L, in the
  tree's own units (mash distance for mashtree, one minus the Jaccard
  similarity of the k-mer sketches for the sourmash builder, substitutions
  per site for the ML builders).
- `--collapse-support S` merges a node whose support is below the fraction S.
  IQ-TREE and RAxML-NG write percentages and are normalised automatically,
  FastTree writes fractions, and mashtree and the sourmash builder write no
  supports, in which case the option warns and does nothing.

**Clade diameter.** `--collapse-diameter D` works on clades rather than
single branches: the fewest clades whose leaf-to-leaf distance stays within D
are kept, and every node inside such a clade is merged into it. A group of
near-identical genomes joined by many short branches therefore becomes one
FlexTaxD node with the genomes as its children (the `max_clade` method of
TreeCluster). D is in the tree's own units:

- mash distance for mashtree (about 0.05 at 95 % ANI for closely related
  genomes);
- one minus the Jaccard similarity at the chosen k for the sourmash builder
  (about 0.88 at 95 % ANI with k = 31, so a mash-scale value is far too strict
  there);
- substitutions per variable site on an ML tree from a SNP alignment.

The diameter rule runs before the other two rules and refuses a tree without
branch lengths. A clade it kept can still be merged into its parent by
`--collapse-length` when its own branch is short, which then puts leaves more
than D apart under one node. An ingroup that fits within D on an
outgroup-rooted tree becomes one node under the root.

**Common rules.** The root, the outgroup/ingroup split under it and the
leaves never collapse. A collapsed node's branch length is added to its
children's. Provenance records the thresholds and the number of nodes
collapsed, so changing one re-runs the stage. On the 171-leaf Wolbachia
mashtree tree, `--collapse-length 0.0005` removes 26 of 168 internal nodes.
`tree2tax_nodes.tsv` beside the two tables describes every node (leaf count,
depth, branch length, support, relative evolutionary divergence, the lineage
of each genome) and is the place to look when choosing a threshold; see
[output.md](output.md).

### Node names

Internal nodes are named by a hash of their leaf set unless
`--node-basename` numbers them. A leaf set changes whenever a genome joins or
leaves it, so every ancestor of the change is renamed on the next run.

`--naming lca` names a node by a hash of the smallest leaf label under each
of its children instead: the node is the last common ancestor of those
leaves and the name is unique per node. Names then change less often, not
never:

- a genome whose label sorts after the present ones (string order, which is
  not always chronological across GCA, GCF and viral prefixes) and attaches by
  splitting a branch leaves every existing name unchanged;
- one that joins a multifurcation (which the collapse options produce) as a
  further child renames that node;
- one that sorts below a child's current smallest label renames the
  ancestors where it becomes the minimum;
- removing a genome renames the nodes whose minimum or child it was.

The tree is also inferred anew on each run, so a topology change elsewhere
renames the nodes it touches. Names are opaque identifiers in either case;
the record carries the rule, or null when a basename numbered the nodes.

### A rank column

`tree2tax.tsv` has two columns by default, since FlexTaxD reads a third as
the rank of every row and an empty value inherits the previous row's.
`--ranks` (on `tree2tax`, `tree2tax-relations` and `run`) adds the column
with `strain` on every genome and `no rank` on every internal node, so a
FlexTaxD database built from it carries a rank column with no taxonomic rank
asserted.

`--rank-red species=0.95,genus=0.8` names internal nodes by their relative
evolutionary divergence (RED, the `red` column of `tree2tax_nodes.tsv`). A
node takes the rank with the largest threshold not above its RED, and
`no rank` when an ancestor already holds that rank, so each named rank
appears at most once on a root-to-leaf path and in threshold order. A lineage
may lack a rank when no ancestor reaches that band.

The thresholds are the user's and the example is illustrative. GTDB's median
RED values (Parks et al. 2018; about 0.92 for genus and 0.77 for family, with
species delimited by ANI rather than RED) come from a bacteria-wide tree and
transfer to a genus-level mashtree or SNP tree only loosely, so inspect
`tree2tax_nodes.tsv` before settling on values.

`--rank-red` needs branch lengths and a rooted tree (an outgroup,
`phylo --root midpoint`, or `--keep-root` on the step) and is refused
otherwise. `run` checks the bands, `--ranks` and the rooting before any stage
runs. The record stores the bands sorted by threshold and carries `ranks` and
`rank_red`, so changing them re-runs the stage.

## Nextflow

The Nextflow pipeline runs the same stages as the command line on a cluster
or in the cloud. Use it when the work should be spread over many tasks, for
example one assembly per sequencing run or one dereplication chunk per node.

The Nextflow layer runs the stages as typed data channels: each stage emits
its outputs (the metadata selection, genome FASTAs, per-chunk and merged
representatives, the tree, the taxonomy) as staged files that the next stage
consumes. There is no shared working directory; results are published under
`--outdir`. Nextflow owns the fan-out (scatter-gather dereplication).

```bash
nextflow run nextflow/main.nf \
    --mode bacterial \
    --metadata_args '-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella' \
    --outdir results \
    -profile standard
```

Run `nextflow run nextflow/main.nf --help` for the parameter summary. The
pipeline requires Nextflow 26.04 or later (`nextflowVersion =
'!>=26.04.0'`). Add `-stub` to any run for a quick wiring check without
external tools.

### Pipeline structure

`nextflow/main.nf` dispatches by `--mode` to one of three data-channel
subworkflows that share the dereplication, phylo and tree2tax modules:

```
bacterial: ACQUIRE       (metadata -> genome)          -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
viral:     VACQUIRE      (vmetadata -> vgenome)        -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
reads:     ACQUIRE_READS (reads -> assemble per run
                          -> [genome-qc] -> gather)    -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
```

`metadata` emits a portable `selection.tsv`; `genome-fetch` downloads the
genomes and emits them as a channel feeding the scatter-gather
dereplication; `phylo` and `tree2tax` run in task-local working directories
and emit `tree.nwk`, `tree2tax.tsv` and `genomes_map.tsv` to `--outdir`.

`PHYLO` (and the split `PHYLO_MSA`/`PHYLO_TREE`) publish the alignment they
built (`phylo/align/` for an aligner, `phylo/tree/msa/` for a SNP typer, with
the reuse stamp) and the tree builder's own files under `phylo/tree/`, beside
the tree.

Every channel carries a meta map built once per run: `id` is a slug of the
selection target (`francisella` for `-tg francisella`), `mode` is
`bacterial` or `viral`. Task tags and the dereplication chunk names use it.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Published results and execution reports. |
| `--mode` | `bacterial` | `bacterial` (GTDB), `viral` (NCBI Virus; BV-BRC via `--vmetadata_args "--source bvbrc"`) or `reads` (ENA/SRA sequencing runs, assembled per run). |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--vmetadata_args` / `--vgenome_args` | see config | Viral metadata / genome selection arguments. With `--group-segments` in `--vgenome_args`, the acquisition task passes `segments.tsv` to `tree2tax-relations`, so `genomes_map.tsv` lists the member segment accessions of each grouped isolate; `segments.tsv` itself is not published. |
| `--reads_args` | see config | Arguments for the reads stage (ENA/SRA run selection by taxon or accession). |
| `--assembler` | `auto` | Assembler for every run in reads mode (`auto` picks by platform). |
| `--assemble_args` | (empty) | Extra `assemble-run` flags in reads mode (`--polisher`, `--min-contig-length`, `--tool-arg`). |
| `--checkm2_db` | `null` | CheckM2 database; switches on quality scoring and the completeness/contamination gate in reads mode. |
| `--gtdb_sketch` / `--gtdb_lineages` | `null` | GTDB sourmash sketch and its lineages CSV; switch on classification and GTDB naming in reads mode. |
| `--reads_sketch` / `--reads_sketch_spots` | `false` / `null` | Sketch every selected run from the archive stream before any download (`READS_SKETCH`), optionally capped at N spots per run. |
| `--reads_screen` / `--reads_screen_args` | `false` / (empty) | Judge the sketches against the GTDB sketch (`READS_SCREEN`) and assemble only the runs that pass; extra `reads-screen` flags. |
| `--derep_tool` | `skder` | Dereplicator for the scatter-gather step. |
| `--derep_process_size` | `null` | Genomes per dereplication chunk (single chunk if unset). |
| `--derep_primary_ani` / `--derep_secondary_ani` / `--derep_aligned_fraction` | `0.90` / `0.99` / `0.50` | ANI / aligned-fraction thresholds. |
| `--derep_keeper` | `quality` | Representative choice in the chunk and merge steps: `quality` (CheckM scores from `selection.tsv`), `gtdb` (a GTDB species representative from `selection.tsv` first, then quality) or `tool`. |
| `--derep_reduce` | `none` | Collapse the merged representatives to one per `species` or `genus` (`dereplicate-merge --reduce`). |
| `--derep_target_reps` | `0` | Search the merge pass's secondary ANI to land near this many representatives (`dereplicate-merge --target-reps`). |
| `--derep_refine` | none | Seed genomes of a lineage to keep at finer resolution, comma-separated; a `DEREP_REFINE` task follows the merge (`dereplicate-refine --refine`). |
| `--derep_refine_ani` | none | Secondary ANI for the seeds' clusters, above `--derep_secondary_ani` (`--refine-ani`). |
| `--derep_refine_file` | none | TSV of `seed` and `ani` for several lineages, each at its own ANI (`--refine-file`). |
| `--phylo_args` | `--treebuilder mashtree` | Aligner or tree builder for the phylogeny, and the rooting (`--root outgroup\|midpoint\|none`). |
| `--phylo_split_msa` | `false` | Run the alignment and the tree as separate tasks. |
| `--sketch` | `false` | Bacterial mode: sketch the genomes and the outgroup with sourmash, publish `sketches/`, and stage them into the dereplication and tree tasks. |
| `--annotator` | `bakta` | With `--msa-source coregene` in `--phylo_args`: the annotator `ANNOTATE` runs on each representative and the outgroup (`bakta` or `prokka`). |
| `--bakta_db` | `null` | bakta database directory, staged into each `ANNOTATE` task; without it bakta reads `BAKTA_DB` from the task environment. |
| `--annotate_args` | (empty) | Extra `annotate-build` arguments (`--tool-arg key=value`). |
| `--tree2tax_args` | (empty) | tree-to-taxonomy (FlexTaxD) arguments; redundant genomes are listed by default (`--no-include-dereplicated` to omit them). `--no-outgroup` is added when `phylo_args` contains it, and `--keep-root` when it contains `--root midpoint`. |

Parameters are validated against `nextflow/nextflow_schema.json` at launch.

#### Passing parameters

Either pass them on the command line (`--name value`) or, for a reusable and
typed configuration, supply a params file:

```bash
nextflow run nextflow/main.nf -profile standard -params-file params.json
```

A starting point is in `nextflow/assets/params_example.json`. On the command
line, give a string parameter whose value starts with a dash as
`--key=value` (for example `--tree2tax_args=--include-dereplicated`):
Nextflow reads `--key --flag` as the boolean `key = true`. A params file
keeps numeric values typed (e.g. `"derep_process_size": 2000`); on the
command line they arrive as strings, which the schema also accepts for the
numeric options.

### Profiles

Combine an executor profile with an optional container profile, e.g.
`-profile slurm,singularity`.

- **Executors**: `standard` (local), `slurm`.
- **Containers**: `docker`, `singularity`, `wave`. These set RepGenR2's own
  adapter-level container backend (`--container ...`; see
  [Container profiles in Nextflow](#container-profiles-in-nextflow)). Under
  `docker` or `singularity` every adapter with a pinned image runs in it (all
  but the `simple` SNP typer, parsnp, skder and SibeliaZ, which run on the
  host and say so in the log); `wave` mints an image from each adapter's
  conda specification instead, which also covers those four. RepGenR2 itself
  must be available to the Nextflow process.
- **`test`**: minimal resources and a small target for a quick smoke run.

The `standard` profile (the default, local executor) caps every process's
CPU and memory request to what the machine has; `base.config` asks for up to
32 CPUs and 128 GB for the heavy processes, which the local executor would
otherwise refuse. `slurm` leaves the requests as they are and sets only the
executor: the queue or account is site-specific, so pass it in a site config
(`-c site.config` with `process.queue = '...'`).

The same holds for a cloud executor such as AWS Batch. The pipeline ships no
cloud profile, because the region, the job queue and an image that provides
`repgenr` and the tools are values only the site knows; a site config that
sets them is all that is needed.

### Reads mode

`--mode reads` builds the genome set from sequencing runs, as
[Starting from sequencing reads](#starting-from-sequencing-reads) does on the
command line, with one task per run.

```bash
nextflow run nextflow/main.nf --mode reads --outdir results \
    --reads_args "-tg mycoplasmopsis --platform illumina --max-runs 20" \
    --assembler auto --checkm2_db /db/checkm2 \
    --gtdb_sketch /db/gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb_lineages /db/lineages.csv \
    -profile slurm,singularity
```

`READS_SELECT` runs `repgenr reads` and emits `reads.tsv`. Every row becomes
one `READS_ASSEMBLE` task (`assemble-run`, label `process_assembly`: 8 CPUs,
16 GB for a short-read run and 32 GB for a long-read one, 6 h, scaled by
attempt), so a hundred runs assemble in parallel on a cluster. When
`--checkm2_db` or `--gtdb_sketch` is set, one `GENOME_QC` task (`genome-qc`,
`process_medium`) scores and classifies the batch. `READS_GATHER`
(`reads-gather`) applies the gate and the naming policy and emits the genome
FASTAs, `selection.tsv` and an empty outgroup accession file, the same tuple
`ACQUIRE` emits. `reads.tsv`, `selection.tsv`, `assembly_stats.tsv`,
`excused_runs.tsv` and the QC tables are published under `reads/`.

The number of runs is bounded on the selection with
`--reads_args '--max-runs N'`, and the reads sketch and screen of the
assembly step are set through `--assemble_args` (see
[Stateless steps](#stateless-steps)).

#### Reads in Nextflow

`--reads_sketch` adds a sketch task between the selection and the assemblies,
and `--reads_screen` a screen task after it (the screen implies the sketch).
`READS_SKETCH` (`reads-sketch --reads-tsv reads.tsv -o sketches`,
`process_low`) streams every run of `reads.tsv` from the archive into a
sourmash sketch, the first `--reads_sketch_spots` spots of each when set,
and publishes them under `reads/reads_sketches/<variant>/`. `READS_SCREEN`
(`reads-screen --reads-tsv reads.tsv --sketches-dir sketches -o screen`,
`process_medium`) judges the sketches with the three checks of the assemble
screen against `--gtdb_sketch` and `--gtdb_lineages`, with
`--reads_screen_args` for `--min-fraction` and `--dup-ani`, and publishes
`reads/reads_screen/<variant>/` with `reads_screen.tsv`, one decision per
run and `passed.txt`. Only the runs of `passed.txt` become `READS_ASSEMBLE`
tasks; the others are neither fetched nor listed in `excused_runs.tsv`,
since they never reach the assembly step (their reasons are in
`reads_screen.tsv`). That includes a run whose stream failed in the sketch
task (`screen_failed: no reads sketch`), which without the screen would
still have been assembled; a transient archive error therefore drops the
run from the batch, so check `reads_screen.tsv` for such rows. A screen that
passes no run stops the pipeline, and `--reads_screen` without
`--gtdb_sketch` and `--gtdb_lineages` is refused before any task.
`--reads_sketch` alone runs the sketch task and publishes the sketches
without a screen; both parameters apply to `--mode reads` only. Both tasks need `fastq-dump`
(sra-tools) and sourmash in the task environment; the screen also needs the
GTDB sketch on a path every task can read.

```bash
nextflow run nextflow/main.nf --mode reads --outdir results \
    --reads_args "-tg mycoplasmopsis --platform illumina --max-runs 50" \
    --reads_screen --reads_sketch_spots 200000 \
    --gtdb_sketch /db/gtdb-rs226-reps.k31-sc10k.sig.zip --gtdb_lineages /db/lineages.csv
```

### Genome sketches in Nextflow

`--sketch` adds a `SKETCH` task on the bacterial path. It writes
`sketches/<record name>.sig.zip` for every genome and the outgroup (k=21, 31
and 51 at scaled=1000), publishes the directory under `--outdir`, and stages
it into every `DEREP_CHUNK`, `DEREP_MERGE` and `PHYLO` (or `PHYLO_TREE`) task,
which pass it to repgenr as `--sketches-dir`. A sourmash dereplicator or tree
builder then reads the sketches instead of sketching; other tools ignore
them.

`SKETCH` calls sourmash itself, so the parameter is off by default and a run
without sourmash is unaffected; the viral and reads paths do not stage
sketches. In the stateless steps a sketch is matched to a genome by record
name only, without the digest check of a working directory, which holds
because `SKETCH` sketches the same genome files the steps receive.

### Split alignment and tree tasks

`--phylo_split_msa` splits the phylogeny into `PHYLO_MSA` and `PHYLO_TREE`
(see [Reusing an alignment across tree builders](#reusing-an-alignment-across-tree-builders)).
The alignment then keeps its own cache entry, so trying another tree builder
or another bootstrap re-runs only the tree, and the two halves take their own
resource labels: the alignment is `process_high`, the tree `process_medium`.
It applies to tree builders that consume an alignment; leave it off for
mashtree or sourmash, which build from genomes.

### Core-gene route in Nextflow

With `--msa-source coregene` in `--phylo_args`, the bacterial and reads
paths annotate before the phylogeny. `ANNOTATE` runs `repgenr annotate-build`
once per representative, and once for the outgroup unless `--phylo_args`
holds `--no-outgroup`, so the annotations run in parallel and Nextflow caches
them genome by genome. A re-run with another tree builder or core-gene
builder annotates nothing again. The record directories are published under
`annotations/` and gathered into the `PHYLO` task (or `PHYLO_MSA` with
`--phylo_split_msa`), which indexes them with `repgenr annotate-index` and
passes them to `phylo-build` as `--annotations-dir`.

- The route is detected from `--phylo_args` (quotes around a value are
  dropped and the last `--msa-source` wins, as on the command line); a
  `withName` override of the phylo task's `ext.args` is not seen.
- `--bakta_db` is linked into each task. Under `stageInMode 'copy'`, or on a
  cloud executor without a shared file system, it would be copied once per
  genome, so there give the tasks a path they can read through `BAKTA_DB`
  instead.
- Each record directory carries `annotation.json` with the annotator, its
  version and the database version, and `annotate-index` indexes each record
  with the tool named there (its `--tool` applies only to records without the
  file), so a set annotated by two tools is refused by the core-gene pass
  rather than read under one name.
- The viral path has no annotation step, since the annotators are
  prokaryotic; with `--msa-source coregene` in `--phylo_args` it logs a
  warning, as the phylo task would then annotate every genome itself.
- The annotator, the core-gene builder and the tree builder must each be
  reachable from their tasks: with `--repgenr_opts '--bin-dir
  bakta=<env>/bin --bin-dir ppanggolin=<env>/bin'` for satellite
  environments, or a container profile.

### Scaling

For large dereplication inputs (10k+ genomes), set a chunk size so the
dereplication scatters across tasks (one per chunk):

```bash
nextflow run nextflow/main.nf --outdir results \
    --derep_tool sourmash --derep_process_size 2000 -profile slurm -c site.config
```

Resource labels (`process_low/medium/high`, and `process_assembly` for the
per-run assembly tasks of reads mode, whose memory follows the read platform)
scale memory and time with the retry attempt, so a task killed for memory or
time is resubmitted with more headroom. Tune the label values per environment
in `nextflow/conf/base.config`.

#### Scatter-gather dereplication

For horizontal scaling, the `DEREPLICATE_SCATTER` subworkflow groups genomes
into chunks of `--derep_process_size`, dereplicates each chunk as a separate
task (one per node on HPC), and dereplicates the union of the chunk
representatives once more (the two-stage reduce-tree expressed as typed data
channels). It wraps the `repgenr dereplicate-chunk` / `dereplicate-merge` CLI
steps and is driven by `--derep_tool` and the `--derep_*_ani` /
`--derep_aligned_fraction` thresholds. Compared with
[Chunked dereplication](#chunked-dereplication) on the command line, it runs
one chunk round and one merge task over the whole union.

With `--derep_refine SEEDS --derep_refine_ani F` (comma-separated seeds) or
`--derep_refine_file refine.tsv`, a `DEREP_REFINE` task follows the merge:
the merged result and every genome of the run are staged in, the clusters of
the seeds are dereplicated again at the finer ANI (`repgenr
dereplicate-refine`), and the refined result, published as
`dereplicate/refined/` with its `refined.tsv`, is the one the tree is built
from; `dereplicate/merged/` stays published beside it as the coarse result
(see [Finer resolution inside a lineage](#finer-resolution-inside-a-lineage)).
The refined directory holds the contract files and no `genome_n50.tsv`.
Seeds without an ANI, an ANI without seeds, or an ANI not above
`--derep_secondary_ani` are refused at launch. The reads route runs the
same subworkflow, so the refinement applies there too.

A standalone harness runs it on a directory of genome FASTAs (no GTDB/NCBI
front-end), useful for testing and for dereplicating a local collection:

```bash
nextflow run nextflow/tests/dereplicate_scatter.nf -c nextflow/nextflow.config \
    --genomes_dir <DIR> --derep_tool sourmash --derep_process_size 2000 \
    --outdir results -profile standard
```

Add `-stub` to exercise the wiring without running the tools.

### Configuring processes

Each Nextflow process reads its tool flags from `task.ext.args`, which
`nextflow/conf/modules.config` maps from the user-facing parameters:
`--metadata_args`, `--vmetadata_args` and `--vgenome_args`, `--phylo_args`,
`--tree2tax_args`, and five dereplication parameters (`--derep_tool`,
`--derep_primary_ani`, `--derep_secondary_ani`, `--derep_aligned_fraction`,
`--derep_keeper`) composed into one string for the two dereplication
processes. `--derep_reduce` and `--derep_target_reps` are passed to the merge
process only, as `ext.args2`. `--derep_process_size` is read by the scatter
subworkflow to size the chunks and is not a tool flag. Publishing directories
live in the same file, and resources and the retry window in
`nextflow/conf/base.config`.

A site can retune one process without touching the pipeline by passing its
own config:

```groovy
// site.config
process {
    withName: 'PHYLO' {
        ext.args = '--treebuilder iqtree --bootstrap 1000'
    }
    withName: 'DEREP_CHUNK|DEREP_MERGE' {
        ext.args = '--tool galah --secondary-ani 0.98'
    }
}
```

`nextflow run nextflow/main.nf -c site.config ...`. `ext.repgenr_opts` is
the hook for the top-level repgenr options every process prepends (the
container profiles set it through `--repgenr_opts`).

An `ext.args` override replaces the whole string for that process; for the
dereplication processes that means every flag not repeated in the override
falls back to the repgenr CLI default rather than to the `--derep_*`
parameter.

Nextflow resolves both `bin/` and the default config relative to the
launched script. The harness scripts under `nextflow/tests/` therefore reach
the versions helper through the symlink
`nextflow/tests/bin/repgenr_versions_fragment`, and running one of them
directly needs `-c nextflow/nextflow.config` (nf-test supplies the config
itself). Any future script added to `nextflow/bin/` needs its own link under
`nextflow/tests/bin/`.

## Running tools in containers

RepGenR2 can run each external tool inside a pinned container instead of
relying on tools installed on `PATH`. Use it to run Linux-only tools on any
host (e.g. macOS) and to make tool versions reproducible. RepGenR2 itself runs
on the host; only the tool subprocess is containerized, so it works in the
plain CLI and inside Nextflow.

### Quick start

```bash
# Run every tool in a container via Docker:
repgenr --container docker dereplicate -wd $WD --tool drep

# Singularity/Apptainer (HPC), with .sif images on a large disk:
repgenr --container singularity --container-cache /Volumes/LaCie/repgenr_sif \
        phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
```

`--container` is a top-level option (place it before the subcommand). The
default is `none` (native execution, unchanged).

### Options and environment variables

Each option is top-level, or set through its environment variable.

| Option | Env var | Meaning |
|--------|---------|---------|
| `--container {none,docker,singularity}` | `REPGENR_CONTAINER` | execution backend |
| `--container-engine <bin>` | `REPGENR_CONTAINER_ENGINE` | engine override (apptainer, podman) |
| `--container-cache <dir>` | `REPGENR_CONTAINER_CACHE` | Singularity `.sif` images and cache (can be external) |
| `--platform <plat>` | `REPGENR_CONTAINER_PLATFORM` | e.g. `linux/amd64` to emulate BioContainers on Apple Silicon |
| `--wave / --no-wave` | `REPGENR_WAVE` | resolve multi-tool/arm64 images via the Seqera Wave CLI |
| `--bin-dir TOOL=DIR` (repeatable) | `REPGENR_BIN_DIRS` (`tool=dir,tool=dir`) | put DIR first on `PATH` for one tool on the host (a satellite conda environment); see below |

Image sources, the storage location, Apple Silicon and Rosetta setup, and the
notes on bind mounts and known image problems are in
[install.md](install.md#3-containers).

### Per-tool directories (`--bin-dir`)

`--bin-dir TOOL=DIR` puts DIR first on `PATH` for one tool only. Use it when
a tool lives in its own conda environment (a satellite environment) rather
than in the environment that runs repgenr. DIR comes first for the tool's
preflight lookup and version query, and for every host subprocess it runs, so
a wrapper such as `run_gubbins.py` calls the helpers of its own environment.
Other tools keep the inherited `PATH`.

```bash
P=$(conda info --base)/envs
repgenr --bin-dir gubbins=$P/repgenr-gubbins/bin --bin-dir mashtree=$P/repgenr-mashtree/bin \
        phylo -wd $WD --msa-source snptype --snptyper simple --mask gubbins
```

TOOL is an adapter name as printed by `repgenr list-tools` (for example
`gubbins`, `mashtree`, `snippy`, `parsnp`, `progressivemauve`), or `checkm2`,
`datasets` or `minimap2` (racon's mapper). The option is repeatable, and
`REPGENR_BIN_DIRS="gubbins=/p/bin,mashtree=/q/bin"` sets the same; an option
given on the command line replaces the variable's entry for that tool. An
unknown tool or a directory that does not exist exits 2.

The directory each tool used is recorded under `bin_dirs` in the stage's
record in `repgenr.yaml`. It is not part of the resume fingerprint, so
setting or changing it does not rerun a finished stage; the tool's version is
recorded as before. Under `--container`, a tool that runs in an image is not
affected, and a warning names it.

In Nextflow, pass the option through `params.repgenr_opts`, for example
`--repgenr_opts '--bin-dir gubbins=/p/bin'`. This value replaces the one a
container profile sets, so with `-profile docker` give
`--repgenr_opts '--container docker --bin-dir gubbins=/p/bin'`.

### Container profiles in Nextflow

`nextflow run nextflow/main.nf -profile docker` (or `singularity`, `wave`)
sets `params.repgenr_opts` so every stage calls `repgenr --container ...`
(see [Profiles](#profiles)). On HPC, Singularity is the more direct choice
(RepGenR2 dispatches each tool to Singularity); stacking with Nextflow's own
Docker engine implies docker-in-docker.

## Troubleshooting

This part lists common failures and how to diagnose them, then the exit codes.

### A tool is not found (`MissingBinaryError`)

The Python package does not install the bioinformatics tools (see
[install.md](install.md)). Put the tool on `PATH` from a conda environment,
or use `--container docker` (or `singularity`) to run tools in pinned images
instead (see [Running tools in containers](#running-tools-in-containers)).
Run `repgenr list-tools` to see the adapters and each tool's declared genome
limit; `list-tools --check` also runs every adapter's preflight and reports
which binaries are missing or too old.

### Apple Silicon / arm64

BioContainers are amd64; pass `--platform linux/amd64` (and enable Rosetta)
so emulated images run.

### A stage failed; where are the details?

Errors print a concise message; the full traceback is in
`<workdir>/repgenr.log`. A failed external tool prints one line naming the
tool and its exit status; the command line and the output tail are in the
same log. Re-run with `--verbose` to see them on the console (a data-channel
step has no log and always prints the tail).

A tool command line is shortened on the console, with or without
`--verbose`, so that the whole line fits 120 columns:

- a run of input paths is shown as its first path and a count, then other
  paths by their last component, then the run as its count alone and other
  paths as `...`;
- a path that directly follows an option (`--reference ref.fasta`) is that
  option's value and is not counted with the inputs;
- long lists of genome names show the first five and a count.

Under `--container` the console shows the tool's own command; the engine
options (mounts, labels, user) are in the run log. The full text of such a
line, and the tool's own output, are in `repgenr.log`.

`repgenr status -wd <WD>` shows what completed and what is next; a stage that
failed is listed as `[interrupted]` and its outputs may be partial until it
is re-run. The other states are described under
[Checking a working directory](#checking-a-working-directory).

### GTDB download fails

Check `--release` (e.g. `232.0`) and `--gtdb-version` (`bac120`/`ar53`);
transient HTTP errors are retried automatically. A host that does not accept
a connection within 15 s counts as unreachable, so a blocked network exits 3
after about two minutes of retries. The `--source api` mode fetches only the
target taxon (no full-table download).

### Genomes NCBI no longer serves

A GTDB release can list assemblies that NCBI has since suppressed. `genome`
records them in `missing_accessions.txt` (also when a whole download batch
consists of them) and later stages excuse them; a re-run asks for them again.
Each rehydrated genome is checked against the package's `md5sum.txt`, and a
mismatch is recorded the same way. An outgroup NCBI does not serve exits 3;
choose another with `metadata --outgroup-accession`.

### NCBI datasets on a blocked network

Before the first `datasets` call, `genome` and `vmetadata` (NCBI Virus
source) send one request to `api.ncbi.nlm.nih.gov` with the 15 s connect
timeout. Before each `datasets rehydrate`, `genome` sends one request to
every host named in the package's `fetch.txt` (with datasets 18.x that is
again `api.ncbi.nlm.nih.gov`). When a host cannot be connected to, the stage
exits 3 at once and names it, instead of waiting through three `datasets`
attempts of several minutes each.

The request uses the proxy settings in `HTTPS_PROXY`, `HTTP_PROXY` and
`NO_PROXY`, as the GTDB and Entrez requests and `datasets` itself do. Only a
failed connection or a connect timeout counts as unreachable. An HTTP error
status, a read timeout and a TLS error count as reachable. A TLS error can
mean only that Python's certificate bundle lacks a proxy's certificate
authority that the system store holds, which `datasets` uses;
`REQUESTS_CA_BUNDLE` points Python at another bundle.
`REPGENR_SKIP_NET_PROBE=1` skips the check. A failure later, inside
`datasets`, still exits 6. `vgenome` works from the `vmetadata` download and
makes no network request.

### NCBI Entrez throttling (viral BV-BRC path)

Set `NCBI_API_KEY` (and optionally `NCBI_EMAIL`) to raise the request-rate
limit. An HTTP error is retried per batch of taxids; a connection error (no
network, or the host does not answer) stops the lookup at the first batch
with exit 3.

### A tool hangs

Set `REPGENR_SUBPROCESS_TIMEOUT=<seconds>` to cap every external tool; on
expiry the process group is killed with a clear error. Under
`--container docker` the tool's container is stopped as well.

### Stopping a run

Each external tool runs in its own process group, together with the helpers
it starts (for example IQ-TREE under `run_gubbins.py`). On SIGTERM (`kill`, a
scheduler), SIGHUP (the terminal closes) or Ctrl-C, repgenr sends SIGTERM to
every running tool group, starts no queued tool, sends SIGKILL to the groups
still present after 5 s, and exits: 143 for SIGTERM, 129 for SIGHUP, and the
usual `KeyboardInterrupt` for Ctrl-C. The stage is then listed as
`[interrupted]` and a partial output is removed. A second signal kills the
tools at once.

A signal repgenr starts with as ignored stays ignored, so
`nohup repgenr ... &` survives the terminal closing. Ctrl-Z suspends the
tools with repgenr, and `fg` or `bg` resumes them; without a terminal
(`setsid`, a workflow manager) SIGTSTP has no effect, as for any process in
that situation. SIGKILL to repgenr cannot be handled, so its tools keep
running; each tool's process group ID is its own PID, so
`pkill -g <tool pid>` stops a tool with its helpers.

**Under `--container docker`.** The container runs with `--init`, so the
SIGTERM the `docker run` client forwards ends the tool in the container. Each
container is named `repgenr-<repgenr pid>-<hex>`, and when repgenr is stopped
it also runs `docker stop -t 0` on the container of the tool it was running,
so a tool that ignores SIGTERM does not keep running. This also applies to
containers started from parallel worker threads, which delays the exit by
about a second (the stops run concurrently). A second signal ends repgenr at
once and starts one `docker stop -t 0` for all its running containers, which
completes after repgenr has exited. Each container carries the labels
`repgenr.pid=<repgenr pid>` and `repgenr.host=<host name>`. SIGKILL to
repgenr leaves its containers running; stop them by label:

```bash
docker ps -q --filter label=repgenr.pid=<repgenr pid> | xargs -r docker stop
```

### Exit codes

A script can tell the failure classes apart without reading the log:

| Code | Meaning |
|------|---------|
| 0 | Success. |
| 1 | An unexpected error (traceback in the run log). |
| 2 | Invalid or missing user input (also Typer's own usage errors). |
| 3 | The working directory does not exist (every command, including `status` and `doctor`) or is missing files or in a bad state, or a request to a remote service (GTDB, NCBI Entrez, NCBI Datasets, BV-BRC, ENA) failed, e.g. because the network is unreachable. `genome` and `vmetadata` on NCBI Virus check that the NCBI Datasets host answers before running the `datasets` CLI; a network failure inside `datasets` is reported as 6 instead. `assemble` and `reads-gather` also exit 3 when every run was excused and nothing was produced; the reasons are in `excused_runs.tsv`. When an earlier `assemble` call wrote the genome set and at least one run was judged (rejected by the assembler, polisher or quality gate), `assemble` empties that set, so `dereplicate` also exits 3 instead of running on stale genomes; a set written by another stage, one kept because every run failed to download, or one under `--append` stays in place. See [output.md](output.md#when-every-sequencing-run-is-excused). |
| 4 | A required external tool is absent or below its version floor, or the Docker daemon cannot be reached under `--container docker`. `list-tools --check --strict` exits 4 when any adapter is missing or errored. |
| 5 | A requested tool adapter could not be found or loaded. `list-tools --check --strict` exits 5 when any plugin failed to load (this takes precedence over 4). |
| 6 | An external tool failed (one console line; the command and output tail are in `repgenr.log`). Under `REPGENR_PROPAGATE_TOOL_EXIT=1` (set by the Nextflow modules) the tool's own status is forwarded instead, a signal kill as 128 plus the signal number. |
| 7 | `doctor` found at least one failure: an output that disagrees with its record, a malformed `repgenr.yaml`, or a check that could not complete. Warnings alone exit 0. |
