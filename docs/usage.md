# Usage

This page is the command-line guide for RepGenR2 (the command is `repgenr`).
It follows the order of an analysis: building a genome set, dereplicating
it, inferring a phylogeny, and deriving a FlexTaxD taxonomy from the tree.
It is written for readers who know genomics but not this tool.

Three pages go into more depth. [reads.md](reads.md) covers genome sets
from sequencing runs, [phylogeny.md](phylogeny.md) the phylogeny routes, and
[nextflow.md](nextflow.md) the Nextflow pipeline. [install.md](install.md)
explains how to install the package and the
external tools, and [choosing-tools.md](choosing-tools.md) explains which
dereplicator, phylogeny route and tree builder suit a dataset. Every command
and option is listed in [cli-reference.md](cli-reference.md), and the files
each stage writes are described in [output.md](output.md).

## Contents

1. [Introduction](#introduction): stages, the working directory, threads,
   running the whole chain, resume and `--force`, checking a working
   directory
2. [Building the genome set](#building-the-genome-set): counting taxa,
   bacteria from GTDB, local genomes, sequencing reads, viruses, genome
   sketches
3. [Dereplication and representative selection](#dereplication-and-representative-selection):
   representative selection, chunked dereplication, inspecting a
   dereplication
4. [Phylogeny](#phylogeny): the routes, with one example each
5. [tree2tax](#tree2tax): rooting, collapsing weak splits, node names, a
   rank column
6. [Running tools in containers](#running-tools-in-containers)
7. [Troubleshooting](#troubleshooting), including working directories from
   earlier versions and the exit codes

## Introduction

### Stages and entry points

RepGenR2 builds a representative-genome repository for a taxonomic group and
the matching FlexTaxD taxonomy. A representative genome is the genome kept
for a cluster of closely related genomes. The pipeline runs five stages in
order:

1. **metadata** selects genome accessions from GTDB.
2. **genome** downloads and organises the selected genomes with NCBI
   Datasets.
3. **dereplicate** clusters the genomes by ANI (average nucleotide identity)
   and picks representatives.
4. **phylo** builds a phylogeny from the representatives.
5. **tree2tax** emits a FlexTaxD-compatible taxonomy from the tree.

The first two stages can be replaced by another entry path (the routes are
listed in the [README](../README.md#pipeline-routes)):

- **Bacteria and archaea from GTDB:** `metadata -> genome`, see
  [Bacteria from GTDB](#bacteria-from-gtdb).
- **Genomes already on disk:** `ingest`, see
  [Starting from local genomes](#starting-from-local-genomes).
- **Sequencing runs from ENA/SRA:** `reads -> assemble`, see
  [reads.md](reads.md).
- **Viruses from NCBI Virus or BV-BRC:** `vmetadata -> vgenome`, see
  [Viruses](#viruses).

Each entry path writes the same genome set (`genomes/`, `selection.tsv` and
the manifest), so `dereplicate`, `phylo` and `tree2tax` are the same for all
of them. An entry path can also stage an outgroup, a genome from outside the
target taxon that is used to root the tree. The Nextflow layer runs the same
work as stateless data channels without a shared working directory (see
[nextflow.md](nextflow.md)).

### The working directory

Every stage command takes the working directory with `-wd` (long form
`--workdir`). Stages read and write inside it, and `repgenr.yaml` there
records what ran. Large run-time data (scratch, downloads, the Cactus job
store) also lives under the working directory and `TMPDIR`, so put it on a
disk with room.

### Threads

Commands that run a threaded tool take `-t/--threads`. Without it the count
is 16, or the CPU limit of the repgenr process when that is lower. The CPU
limit is the set of CPUs the process may run on (a Slurm cpuset, `taskset`)
and a cgroup CPU quota (`docker run --cpus`, a Kubernetes CPU limit),
rounded up. A lowered default is reported once on the console. An explicit
`-t` is used as given, with a warning when it exceeds that limit. The thread
count is not a result parameter, so changing it does not make a finished
stage rerun.

### Running the whole chain

`repgenr run` runs the whole chain in one command (bacterial by default;
`--viral` for the NCBI Virus path). `repgenr status` shows progress at any
time:

```bash
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` forwards the stage options it shares a name with. For the rest, run
the stage command by hand; `repgenr status` says which comes next. Two
options of `run` change the chain:

- `--with-snptype` adds the standalone `snptype` stage after dereplication,
  so the SNP tables under `snp/` are produced even when the tree is built
  another way.
- `--genomes-dir` starts the chain from local genomes (see
  [Starting from local genomes](#starting-from-local-genomes)).

With `--msa-source snptype`, `phylo` runs its own typing pass, outgroup
included, into `tree/msa/`. `snp/` holds only the `snptype` stage's tables,
so the two do not replace each other, and a later `phylo` run that changes
only the tree builder reuses its alignment.

### Resume and `--force`

Each stage records its parameters, the digests of its inputs, and the
container identity in `repgenr.yaml`. Rerunning a stage is a no-op only when
all three are unchanged, and it logs that it skipped. A stage reruns when:

- one of its parameters, or `--container`, changed;
- an input changed, including the output of an upstream stage that was
  rerun (for example after `repgenr --force dereplicate ...`, `phylo` reruns
  the next time it is invoked);
- one of its main outputs is missing, such as a genome deleted from
  `genomes/`, which the next run of the entry stage restores;
- its last run failed or was interrupted.

`--force` reruns a stage explicitly. It is a global option and goes before
the command name (`repgenr --force dereplicate -wd WD`); after the command
name it is rejected as an unknown option.

A stage that refuses (exit 2, 3 or 4) before changing any of its main
outputs leaves its last finished record in place. To drop a genome on
purpose, remove it from the source directory or from the `--selection` file
and rerun the entry stage. Two changes are not detected, and both need
`--force`:

- an in-place edit of an input genome that preserves its size and
  modification time;
- an upgrade of a natively installed tool binary (switching the container
  back-end or platform is detected).

The full rules, including which outputs each stage checks and which files
count as inputs, are in [developing.md](developing.md#resume-semantics).

### Checking a working directory

`repgenr status -wd <WD>` shows what completed and what is next. `repgenr
doctor -wd <WD>` verifies a working directory's outputs against its records:
missing, corrupt or untracked genomes, manifest drift, truncated
deliverables and interrupted stages.

`doctor` exits 7 on failures and 0 when it finds only warnings. A missing
deliverable or a changed input is a warning, since the stage will rerun.
`doctor --quick` skips reading the first bytes of each genome file, the
slowest check on a large genome set. Links, missing and untracked genomes
are still checked, and the genome line says the content was not read, so a
non-FASTA file saved under a FASTA name goes unnoticed.

`status` reads the same tables as `doctor`, so the two commands agree on
what reruns:

- A stage that failed is listed as `[interrupted]`, and its outputs may be
  partial until it is rerun.
- A finished stage is listed as `[stale]` when one of its inputs changed or
  one of its outputs is missing since it finished. It reruns on its next
  invocation.
- A finished stage whose record holds no resume fingerprint (a record from
  an earlier version, or one restored by `derep-stock --action unpack`)
  stays `[done]` with a note that its next invocation recomputes it.
  `doctor` warns about every completed record without a fingerprint, stale
  ones included.
- When `repgenr.yaml` records no stage, `status` names the entry stages, and
  suggests `doctor` when the working directory holds outputs. When it
  records stages but no entry stage, `status` follows the stages all
  pipelines share (dereplicate, phylo, tree2tax) under "Pipeline: unrecorded
  entry stage".

A malformed `repgenr.yaml` is reported by both: `doctor` as a failure,
`status` and every stage with exit 3. Both commands accept `--json` and then
print one versioned JSON object on stdout instead of the text report, for
scripts and workflow wrappers. The exit codes are unchanged, and on exit 3
stdout stays empty. The schemas are in
[output.md](output.md#machine-readable-status).

## Building the genome set

This part describes `census`, which counts what a taxon or a working
directory holds, the entry paths that produce a genome set, and the sourmash
sketches written beside it.

### Counting genera, species and samples

`repgenr census` counts what a taxon holds before a run, or what a working
directory holds after its entry stage. Use it to choose a target and a route
before downloading, and to check a selection afterwards. It is read-only: it
records no stage and writes nothing into a working directory. A family is
counted per genus and a genus per species, with the totals on the first
line. The command has three modes, chosen by its arguments.

#### A taxon, without a working directory

`-tf/--target-family FAMILY` or `-tg/--target-genus GENUS` counts the GTDB
genomes of the taxon and its GTDB species representatives.

- `--source api` (the default) sends one request to the GTDB API, as
  `metadata --source api` does.
- `--source table` reads the GTDB metadata table of `--release`
  (`--gtdb-version` bac120 by default, or ar53) with the parser of the
  metadata stage. The table is downloaded once into
  `$REPGENR_CACHE_DIR/gtdb` (default `~/.cache/repgenr/gtdb`) and reused, or
  `--metadata-path` names a local copy. A table census parses the whole
  GTDB table, which takes about 20 s for release 232.

`--runs` adds the ENA whole-genome sequencing runs under the taxon, as the
reads stage finds them: runs, distinct biosamples, and runs per platform
(Illumina, ONT, PacBio, and `other` for BGISEQ, Ion Torrent and the rest).
The count is logged first, with a warning above 100000 runs. The census
requests only the four run fields it reads, so its search is smaller than
that of the reads stage.

The runs are grouped by the NCBI lineage of their taxid (one lookup per
distinct taxid) and join the GTDB row of the same name.

- A row that holds runs only is marked `(NCBI)` (or `(no genus, NCBI)` for
  taxids without a genus), listed last and left out of the GTDB totals. Such
  a taxon may be a GTDB taxon under another name.
- For a GTDB name with a suffix (`Bacillus_A`), the runs are those of the
  NCBI taxon without it (`Bacillus`, the prokaryote genus). In a genus
  census they join its species rows by epithet, a note says so, and the
  header names the taxon of the runs (`19772 runs of NCBI genus Bacillus`).
- Names are shown as GTDB and NCBI write them (`Bacillus_A thuringiensis_S`,
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
taxon with the request `vmetadata` sends, but for the metadata report only
(no sequences). `--host`, `--complete-only` and `--released-after` filter it
as in `vmetadata`. The isolates column counts the genomes `vgenome
--group-segments` would form before its length filter, an upper bound on
what vgenome selects. The columns are described in
[output.md](output.md#census-table). A genus target is counted per species
and a family target per genus, and `-tf` and `-tg` narrow the report.
BV-BRC offers no metadata-only report, so `--viral --source bvbrc` without
`-wd` exits 2; run `vmetadata --source bvbrc` and count its working
directory instead. `--runs` is not offered for viral taxa.

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
  columns are then empty, and `median_length` is the median of the
  per-taxid median lengths).

The recorded target is used, and `-tf` and `-tg` narrow it. An API working
directory without `gtdb_api_genomes.tsv` is counted from its
`selection.tsv`, with a note.

#### A working directory after genome, vgenome, ingest or assemble

The census counts the genomes of `selection.tsv` that are present in
`genomes/`, one column per manifest source (gtdb, sra, ncbi_virus, bvbrc,
local), the GTDB species representatives, and the dereplication clusters
whose representative falls in the row when `derep/clusters.tsv` exists. The
outgroup is not counted and is named on the first line. When the candidates
of the entry stage are also present, a `candidates` column shows how many
each row had before selection; a species with candidates and no selected
genome appears with 0 genomes. When the entry-stage data cannot be read, the
census counts the selection without the `candidates` column and says why in
a note.

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
stdout and nothing else. Both formats are described in
[output.md](output.md#census-table). The flags that select a taxon query
(`--viral`, `--target`, `--source`, `--release`, `--gtdb-version`, `--runs`,
`--host`, `--complete-only`, `--released-after`) exit 2 with `-wd`, whose
record names the source and target. A command without a taxon or a working
directory exits 2, and a `-wd` that does not exist or a failed request exits
3.

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

With `--source tsv` (the default), `-r/--release` (major.minor, for example
`232.0`) and `--gtdb-version` (`bac120` or `ar53`) are required, because they
select the table's release and domain and are recorded in the provenance.
They stay required with `--metadata-path <table>`, which reads a local GTDB
metadata table instead of downloading one; a path that does not exist exits
2.

The downloaded table stays in the working directory
(`bac120_metadata_r232.tsv.gz`), and `--nodownload` reuses it instead of
downloading it again. Table names carry only the major release, so a
download also writes `<version>_metadata_r<major>.release` with the exact
release. `--nodownload` refuses a table fetched for another minor release
(exit 2). A table without that file, for instance one placed by hand, is
reused with a warning.

The download is checked against its size and against GTDB's `MD5SUM.txt`,
and a checksum mismatch removes the file and exits 3. A release or version
that GTDB does not publish exits 2 after trying the current (`.tsv.gz`) and
the legacy (`.tar.gz`) layout. A network failure exits 3 with its cause and
does not try the second layout.

#### The GTDB API

`--source api` serves GTDB's current release and ignores `-r`,
`--gtdb-version`, `--metadata-path` and `--nodownload` (a warning names
them). The API does not report a release number, so the stage record holds
`release: null` and, in its place, `api_query_date`, the UTC time of the
query (ISO 8601). `status` shows it after the metadata line, and `versions`
prints it as `gtdb_api_query_date`; for the table path they show the release
(`gtdb_release`). To relate a query date to a GTDB release, compare it with
the release dates on the GTDB website.

#### Taxon names and the outgroup

GTDB taxon names are case-sensitive. The target's first letter is raised
and a species epithet is lowered, so `-tg francisella -ts Tularensis` finds
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
per species, and the rule reduces to a quality ranking across species.

Within a species, genomes rank by the CheckM part of the keeper score
(completeness minus five times contamination; see
[Representative selection](#representative-selection)). The ranking happens
before download, so no N50 exists yet. Unscored genomes come last, and ties
go to the GTDB species-representative flag, then the accession, so the
result is deterministic.

On the `--source api` path the quality card of every candidate is fetched
before the cut, one request per genome, four at a time (see the
[genus-scale run](verification.md#genus-scale-run-2026-09-10) for a
duration). The API refuses sustained bursts at times. Refused cards are
retried once more, slowly, and only a genome refused twice is left unscored.

The outgroup is chosen afterwards from the parent taxon, outside the target
taxon (a target genome the cap left out is not a candidate), and never
counts against the limit.

### Starting from local genomes

`repgenr ingest` populates a working directory from local FASTA files
instead of downloading. Use it for genomes already on disk, or to combine
the genome sets of earlier working directories. It stands in for `metadata`
and `genome`, so the local chain is `ingest -> dereplicate -> phylo ->
tree2tax`.

```bash
repgenr ingest -wd $WD --genomes-dir ./my_genomes --outgroup GCF_003574425.1
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
```

`repgenr run --genomes-dir ./my_genomes ...` runs the same local chain in
one command. `--from-workdir`, `--selection`, `--outgroup` and `--copy` pass
through to `ingest`, and the GTDB selection flags are not needed.

`ingest -wd WD --genomes-dir DIR` links the files (or copies them with
`--copy`) into `genomes/` and writes `selection.tsv` and the SQLite manifest
in the same format `metadata` and `genome` produce. The stage is recorded,
so `status` reports the local chain and a second `ingest` on an unchanged
directory skips. Taxonomy and quality columns come from `--selection
selection.tsv` when given, otherwise from the canonical
`Family_genus_species_ACCESSION.fasta` filename.

#### File names, accessions and the outgroup

Without `--selection`, accession and taxonomy come from the filename:

| File name | Accession | Taxonomy |
|---|---|---|
| Four or more `_`-separated tokens, as `Family_genus_species_ACCESSION.fasta` | the tokens after the third (`sample_1_run_A.fasta` gives `A`) | the first three tokens |
| Starts with an NCBI assembly accession, as `GCF_000008985.1_ASM898v1_genomic.fna` from NCBI Datasets or the FTP site | that accession | none |
| Any other name | the file stem | none |

Two genomes with one accession (for example `x.fasta` and `x.fna`, or two
four-token names that end alike) stop `ingest` with exit 2; give them
distinct accessions with `--selection`. The `is_outgroup` column of a
selection takes `1` or `0` (also `true`/`false`, `yes`/`no`). So does the
optional `gtdb_representative` column, which `--keeper gtdb` reads; without
that column no genome is a GTDB representative.

`--outgroup` sets one genome aside under `outgroup/` and writes
`outgroup_accession.txt`, so `phylo` roots on it exactly as after a
download. It names a genome under `--genomes-dir` or in a `--from-workdir`
(filename, stem or accession), or a FASTA file anywhere. The genome is kept
out of the ingroup. When `--selection` also marks an outgroup row, both must
name the same genome; otherwise `ingest` exits 2 and names both. A file from
outside `--genomes-dir` whose name gives the filename or accession of an
ingroup genome is refused (exit 2), since it would replace that genome.

#### Accepted files and compressed genomes

Only the files directly under `--genomes-dir` with a suffix `.fasta`, `.fa`,
`.fna`, `.fas`, `.fasta.gz`, `.fna.gz` or `.fa.gz` are read, so the
`.fna.gz` files of the NCBI FTP site can be ingested as they are. Other
files (`X.FASTA`, `x.fas.gz`) are listed in a warning and skipped.
Subdirectories are not searched. An NCBI Datasets download keeps each
genome in its own directory, so collect the `.fna` files into one directory
first. An empty or unreadable genome file (a dangling link included) stops
`ingest` with exit 2 before anything is staged, and the record of an
earlier `ingest` is left as it was. `doctor` checks that the staged files
hold FASTA.

Compressed genomes are staged unchanged, and every tool accepts them.
Several tools read gzip themselves: the dereplicators (sourmash, skDER,
galah; dRep through a decompressed copy), cactus, the `mashtree` and
`sourmash` tree builders, and the `ska2` and `simple` typers. progressiveMauve, SibeliaZ, `parsnp`
and `snippy` do not, so `snptype` and `phylo` give them a decompressed copy
of each gzipped genome in the stage's scratch directory
(`scratch/snptype/inputs/`, `scratch/phylo_snptype/inputs/` or
`scratch/phylo_inputs/`) and remove the copies when the tool has finished.
The copies keep the genome's record name, so alignments and trees are the
same as for uncompressed files.

The scratch directory therefore needs room for the uncompressed genomes
while the tool runs. The stage checks the free space first (about four
times the gzip size) and stops when less than 1 GB is free. A truncated or
corrupt gzip genome stops the stage with exit 2 naming the file.

#### Links, copies and changed sources

By default `genomes/` holds symbolic links to the source files. If the
source is later moved or deleted the links dangle: `doctor` reports them,
and `ingest` from the new location restages the set. `--copy` makes the
working directory independent of the source. The copies keep the source's
permissions, so read-only sources give read-only copies.

A changed source (files added, removed or rewritten) reruns `ingest` on its
next invocation, prunes genomes no longer present, and marks `dereplicate`
and later stages for a rerun. An unchanged set keeps the dereplication
status recorded in the manifest.

#### Combining working directories

`--from-workdir WD` (repeatable, alone or together with `--genomes-dir`)
takes the genome set of an earlier working directory: the rows of its
`selection.tsv` and the files they name under its `genomes/`. Use it to
dereplicate together genomes selected from GTDB in one working directory
and genomes assembled from sequencing runs in another:

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

`assemble --append` is the alternative within one working directory: it
adds the assemblies to the GTDB selection already there (see
[reads.md](reads.md#adding-assemblies-to-a-gtdb-selection)).

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
over, and `ingest` logs one line for each source that had one. `--outgroup`
sets the outgroup of the new working directory and may name a genome of any
source (filename, stem or accession, under `genomes/` or `outgroup/` of a
source working directory, or under `--genomes-dir`), or a FASTA file
elsewhere.

**Conflicts and resume.** A genome found in two sources (the same
accession, or the same filename with different accessions) stops `ingest`
with exit 2 and names both sources; neither takes precedence. A selection
row whose file is missing, or a `--from-workdir` that is the target working
directory, also exits 2. The selection and the `genomes/` directory of each
source are resume inputs, so a change in a source working directory reruns
`ingest`. The source manifest is not, so a change of a source label alone
needs `--force`.

**State of the sources.** Before staging, `ingest` reads the record of the
stage that wrote each source's genome set (the latest of `genome`,
`vgenome`, `ingest` and `assemble`) and evaluates it as `repgenr status`
would in that working directory. A source whose record is interrupted, or
stale against its own inputs, is named in one warning with the stage and its
state, and its genomes are taken as they are. A source with no such record
(a `selection.tsv` and `genomes/` written by other means) is accepted with
the warning "no stage record". `--strict-sources` turns each of these
warnings into a refusal (exit 2) before anything is written. The state of
each source is kept in the ingest record (`source_states`), and `status`
names the sources that were not done when they were ingested.

### Starting from sequencing reads

This route builds the genome set from whole-genome sequencing runs.
`repgenr reads` selects runs from ENA and writes `reads.tsv`, and `repgenr
assemble` fetches and assembles them into `genomes/`. The chain is `reads ->
assemble -> dereplicate -> phylo -> tree2tax`.

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina --max-runs 20
repgenr assemble -wd $WD --assembler auto -t 16 --jobs 2
# or the whole chain:
repgenr run -wd $WD --reads -ts "Francisella tularensis" --platform illumina --max-runs 20 \
    --assembler skesa --treebuilder mashtree
```

[reads.md](reads.md) describes run selection, assembly and polishing,
quality scoring and classification, screening runs before or during
assembly, assembling in parts, and the stateless steps.

### Viruses

The viral route selects virus genomes with `vmetadata` and stages them with
`vgenome`. It sources from NCBI Virus by default (through the `datasets`
CLI), and `vmetadata --source bvbrc` uses the BV-BRC FTP path instead.

```bash
WD=./hav
repgenr vmetadata -wd $WD --target hepatovirus            # NCBI Virus (default)
repgenr vgenome   -wd $WD --target-genus Hepatovirus      # add --group-segments for segmented viruses
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree
repgenr tree2tax -wd $WD --include-dereplicated
# or: repgenr run -wd $WD --viral --target hepatovirus -tg Hepatovirus --treebuilder mashtree
```

`--virus` passes virus-tuned settings to dRep (`--tool drep --virus`). The
other dereplicators do not read it, and `dereplicate` warns when it is given
with one of them.

#### The viral outgroup

`vgenome` picks an outgroup by itself: a record of a sister species with
enough genomes, chosen by distance (mashtree by default).
`--outgroup-accession` pins it to a downloaded record instead, an accession
on the NCBI Virus path or a record id on BV-BRC, and `run --viral
--outgroup-accession` forwards the same choice. With `--group-segments` the
search still runs, with the kept records' length span widened by 15 percent
as its window, and the outgroup is one record of the sister species rather
than a grouped isolate. `--no-outgroup` leaves the tree unrooted. A run that
ends without an outgroup removes the outgroup an earlier run left in
`outgroup/` and `outgroup_accession.txt`.

#### Species assignment

On the NCBI Virus path the species of a record is the current species of
its taxid in NCBI Taxonomy. After the download, `vmetadata` looks up every
distinct taxid once with `datasets summary taxonomy taxon` and stores the
species, genus and family in `virus_records.json`, so `vgenome` needs no
network.

The organism name is not used for this. NCBI organism names are often
strain-level or earlier names (`Sabia virus` in `Mammarenavirus brazilense`),
and one species can carry several of them (Junin, taxid 2169991, as
`Argentinian mammarenavirus` and `Mammarenavirus juninense`). The lineage in
the NCBI Virus report can nest species that NCBI Taxonomy keeps as siblings,
or lack the current species altogether, so it is not used either.

When the lookup fails, or gives no species for a taxid, the species comes
from the report lineage, and the log says so. It is then the shallowest
lineage name made of a genus (or the genus above a subgenus) and one
lower-case epithet, such as `Mammarenavirus brazilense`, as found on most
records of the taxid.
Earlier names below it that start with the genus (`Hepatovirus A` under
`Hepatovirus ahepa`) are not taken. A
record without either keeps its organism name as the species. Each record's `species_source` (`taxonomy`, `lineage` or
`organism`) says which rule applied, and the `vmetadata` record counts them.
Organism names of the form `<Genus> sp.` (`Orthohantavirus sp.`) are a
species in NCBI Taxonomy, so the unclassified records filed under one such
name share one species token.

The species sets the species token of the canonical filename, the
`--target-species` match, the per-species medians of the length window and
the grouping of outgroup candidates, so strains of the target species are
not taken as sister species. A `--target-species` value also selects the
species of every record whose organism name it is (`-ts "Argentinian
mammarenavirus"` selects all of `Mammarenavirus juninense`), together with
any species token it matches, and the log says so. The organism name stays
in the `description` column of `virus_metadata_base.tsv` and in
`--target-serotype` matching.

#### Segmented viruses

`--group-segments` groups records that share a species and an isolate name
and carry at least two distinct segments into one genome per isolate. Use
it for segmented viruses, whose segments NCBI Virus holds as separate
records.

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

Each isolate keeps one record per segment (complete before partial, then
the longest, then the lowest accession), since NCBI Virus often holds
several submissions of one segment. Records without a segment label stay
single genomes. A grouped isolate's genome carries a synthetic `iso-` token
as its accession, with its first member accession appended when two isolate
names would give the same token. `segments.tsv` records the member
accessions behind it with the normalised segment and the label as
submitted, and `tree2tax` lists the members under the isolate's leaf in
`genomes_map.tsv`.

#### Viral length filtering and over-represented species

The viral selection step keeps records whose genome length falls inside a
window around a center value. The default center, `--length-method
median_of_medians`, limits the influence of over-sequenced species: the
median length is computed per species first, and the window center is the
median of those per-species medians. A species with 900 outbreak records
therefore counts once. `--length-method mean` averages every record and lets
the most-sequenced species set the window; use it only when the selection is
known to be balanced. `--length-range` overrides the window, and
`--length-all` disables the filter.

`vmetadata` records the source and target of `virus_download_wd/download.fa`
in `download.source`. The BV-BRC source reuses the group FASTA only for the
same target and downloads it again otherwise. Switching a working directory
between the two sources removes the other source's tables, so `vgenome`
reads the latest download.

### Genome sketches

`genome`, `vgenome`, `ingest` and `assemble` end by writing a sourmash sketch
of each genome, outgroup included, to `sketches/<name>.sig.zip`. `<name>` is
the genome's file name without its FASTA suffix, the name it carries in
alignments and trees. Each file holds three DNA signatures, at k=21, k=31
and k=51 with scaled=1000 and no abundance. The sourmash dereplicator,
`glance`, the sourmash tree builder and the classifier of `assemble` read
these sketches instead of sketching the genomes again.

- By default the step runs when sourmash can run (on the PATH, or through
  the container back-end) and is otherwise left out with one INFO line.
  `--sketch` requires sourmash and exits 4 without it, before any genome or
  sketch is written. `--no-sketch` leaves the step out. A failed sketch is a
  warning under the default and an error with `--sketch`.
- A sketch is written again only when the genome FASTA (its SHA-256) or the
  sketch parameters changed, so a rerun sketches only new and changed
  genomes. Sketches of genomes that leave the set are removed with them.
- The sketches are not stage outputs. A genome set without them is
  complete, the flag does not change the resume fingerprint, and a stage
  skipped as already complete does not sketch.

`repgenr sketch -wd $WD [-t N]` writes the missing and stale sketches of an
existing working directory, removes those of genomes no longer in the set,
and logs one line with the counts. `repgenr --force sketch -wd $WD` writes
every sketch again. It is recorded as a stage, and `status` shows `sketches:
n/m` (genomes with a sketch, genomes present) on the line of the stage that
wrote the genome set.

A consumer selects the signature at its k-mer size, so `ksize` 21, 31 or 51
at `scaled` 1000 is served from the sketches. It first writes the missing
and stale sketches, and logs a line such as `dereplicate: sketches: 50
reused, 0 written`. The tool sketches into its own working directory when
`--tool-arg ksize=` or `scaled=` asks for other parameters, for a genome
that is not part of the genome set, and when sourmash cannot write a sketch.
How sketches are matched and cached is described in
[developing.md](developing.md#genome-sketches) and
[output.md](output.md#genome-sketches). The reads sketch of `assemble` is
separate and is never copied into `sketches/` (see
[reads.md](reads.md#during-assembly-assemble---screen-reads)).

## Dereplication and representative selection

`repgenr dereplicate` clusters the genome set by ANI and keeps one
representative genome per cluster. It runs after any entry path and before
`phylo`. The dereplicator is chosen with `--tool` (drep, skder, galah or
sourmash; see [choosing-tools.md](choosing-tools.md#dereplicators) for which
suits a dataset):

```bash
repgenr dereplicate -wd $WD --tool skder -t 16
```

The clusters are written to `derep/clusters.tsv` and the representatives to
`derep/representatives/`. This part describes how the representative of
each cluster is chosen, how large sets are dereplicated in chunks, and how
to inspect the result.

### Representative selection

After the dereplicator clusters the genomes, a keeper step re-picks each
cluster's representative. Connectivity-based tools tend to keep the
most-sequenced genome of a cluster rather than the best-quality one, and the
keeper step corrects that. `repgenr dereplicate` and `repgenr run` accept
`--keeper quality|gtdb|tool`. Use the default, `quality`, unless there is a
reason to keep the dereplicator's own pick (`tool`) or the GTDB species
representatives (`gtdb`). The Nextflow steps `dereplicate-chunk` and
`dereplicate-merge` take the same choice (`--derep_keeper`).

**Quality values.** The values come from the manifest:

- GTDB selections carry CheckM values (`--source api` fetches them from each
  genome's card, one request per selected genome);
- `assemble --checkm2-db` scores assemblies;
- `ingest --selection` reads the `completeness` and `contamination` columns.

The Nextflow steps read the same columns from `selection.tsv`. When no
genome under `genomes/` has quality values the stage warns, and
`repgenr.yaml` records `keeper_effective: tool` next to the requested
`keeper` and the swap count.

#### The keeper score

The score follows dRep's default weights for the terms that need no further
tool: `completeness - 5 x contamination + 0.5 x log10(N50)`, with the N50
read once per genome from its FASTA (gzipped or not). The N50 term separates
genomes whose CheckM values differ by less than CheckM's precision, so a
fragmented draft does not outscore a closed genome on CheckM values alone.
Genomes are compared in this order:

1. higher score;
2. higher completeness;
3. lower contamination;
4. higher N50;
5. genome filename, in sort order.

Only genomes with completeness and contamination compete; a genome without
them never becomes the representative. When the dereplicator's
representative has values, the best genome in the cluster wins, that
representative included, so the result does not depend on which genome the
tool picked. When it has none, a scored genome replaces it only if it is
high quality by the MIMAG thresholds (completeness above 90, contamination
below 5). A scored fragment therefore does not displace an unscored complete
genome, and the log names the clusters left with the tool's representative
for this reason. A cluster without any scored genome keeps the
dereplicator's choice.

The working-directory stage applies the keeper once, to the final clusters
of a chunked run. The Nextflow steps apply it within each chunk and again
after the merge, so the merge pass compares the chunk keepers rather than
the tools' picks. Both end with the keeper of each final cluster, but the
clusters themselves can differ slightly, since the merge pass compares
different genomes. Each chunk records the N50 of its scored genomes (see
[output.md](output.md#dereplicate-and-the-dereplication-commands)), so the
merge step scores members whose files it does not receive.

#### GTDB species representatives

`--keeper gtdb` keeps the GTDB species representatives. The metadata stage
records which selected genomes GTDB marks as the representative of their
species (`gtdb_genome_representative` equal to the accession in the table,
`gtdbIsRep` in the API), in the `gtdb_representative` column of
`selection.tsv` and in the manifest.

Within each cluster the keeper is then a GTDB representative when the
cluster holds one. A cluster without one is treated as under `--keeper
quality`, and keeps the tool's pick when no member is scored. GTDB delimits
species at about 95 % ANI, so a `--secondary-ani` below that can join GTDB
species into one cluster. One representative per cluster is kept, chosen
among them by the quality rule above, and a warning names the others, which
stay contained; raise `--secondary-ani` to keep them apart. With `--reduce`,
a GTDB representative is also preferred when representatives of one taxon
are merged.

`repgenr.yaml` records `keeper_effective: gtdb` when the manifest flags at
least one genome, and otherwise `quality` or `tool` with a warning. Only
GTDB selections carry the flag: `ingest` (unless its `--selection` has the
column), `reads`/`assemble` and the viral path leave it at 0, so `--keeper
gtdb` acts there as `--keeper quality`.

#### Quality values passed to the dereplicator

The same manifest values also reach the dereplicator, whichever keeper rule
is chosen. This happens when every genome of the run has both values. dRep receives them
as `--genomeInfo` and then runs without CheckM. galah receives them as
`--genome-info` and ranks genomes by quality. The decision is taken once for
the run, so every chunk and the merge pass use the same source.

When some genomes lack values, none are passed and a warning names them.
dRep then runs CheckM. galah, which otherwise keeps the first listed genome
of a cluster, receives the genomes by descending file size, so a complete
genome rather than a fragment becomes the representative.

#### One representative per taxon

`--reduce species|genus` collapses the ANI representatives to one per taxon
after dereplication. The representative of the largest cluster is the
default keeper of a taxon; with `--keeper quality` and scores known, the
keeper rule re-picks among the taxon's representatives. `--target-reps N`
searches the secondary ANI to land near N representatives. Both options
also exist on `dereplicate-merge`, where the Nextflow pipeline decides the
final set (`--derep_reduce`, `--derep_target_reps`). The merge step takes
the taxonomy from `selection.tsv` when one reaches it and from the
canonical genome filenames otherwise.

#### Finer resolution inside a lineage

One secondary ANI applies to the whole set. A study often needs coarse
coverage of a species and more genomes inside one lineage, such as a
clinically relevant genotype or an outbreak clade. `--refine` keeps that
lineage at a finer resolution without changing the rest of the result.

The lineage is named by seeds: one or a few genomes known to belong to it
(`--refine SEED`, repeatable: a filename, a filename without its suffix, or
an accession), with a finer ANI (`--refine-ani`, above `--secondary-ani`).
The dereplication runs as usual. Every cluster that holds a seed is then
taken as the lineage, and its genomes, the representative included, are
dereplicated again with the same tool and settings at the finer ANI. The
clusters this gives replace the coarse one. Nothing outside the seed
clusters changes, and a seed alone in its cluster has nothing to refine and
is reported.

```bash
# In one run: the species at 0.99, the lineage of GCA_000123 at 0.999
repgenr dereplicate -wd $WD --tool skder --refine GCA_000123 --refine-ani 0.999

# After the fact, on the result on disk, without repeating the coarse pass
repgenr derep-stock -wd $WD --action pack --name coarse
repgenr derep-refine -wd $WD --refine GCA_000123 --refine-ani 0.999
```

- `derep-refine` takes the tool and thresholds from the completed
  `dereplicate` record, so it needs one (a record restored by `derep-stock`
  counts). It rewrites the dereplication outputs in place and records a
  `derep_refine` stage, and `phylo` and `tree2tax` then rebuild from the
  refined set.
- A later `dereplicate` (forced, or with changed inputs) recomputes the
  coarse result and discards the refinement, so store the result to keep
  with `derep-stock` first.
- `repgenr run ... --refine SEED --refine-ani F` passes both to its
  `dereplicate` stage.
- A seed given as an accession needs the version suffix the filename
  carries (`GCA_000123.1` when the file is named so). Unknown seeds, the
  outgroup, a genome the dereplicator rejected and an ANI that is not finer
  are refused before the coarse pass runs.
- The keeper rule chooses the representative of each refined cluster.
  `--reduce species|genus` leaves the refined clusters apart, since they
  share a species by construction; the species then keeps one
  representative outside the lineage and the refined ones inside it.
- `--refine` cannot be combined with `--target-reps`, whose search sets the
  secondary ANI itself. The outgroup is never part of a lineage.

Several lineages, each at its own ANI, are given in a file: `--refine-file
refine.tsv`, a TSV with the columns `seed` and `ani`, one seed per row.
Blank lines and lines starting with `#` are skipped.

```text
seed	ani
GCA_000123.1	0.999
GCA_000456.1	0.999
GCA_000789.1	0.9995
```

The seeds of one ANI form one lineage, and a cluster claimed by two ANIs is
refused (seeds whose coarse cluster was refined before count as one
cluster). The file can be combined with `--refine` and `--refine-ani`, which
apply to the `--refine` seeds only. Its content is part of the resume
inputs, so editing it reruns the refinement, and the path is stored as an
absolute path.

`derep-refine` can be applied more than once, for another lineage or for
the same seed at another ANI. A seed refined before brings in its original
coarse cluster, read from `refined.tsv`, so refining it again at a finer ANI
starts from the same genomes and not from the fine cluster it sits in now.
The rows of lineages left untouched are kept. On a result reduced with
`--reduce`, a cluster on disk is a whole taxon group, and the lineage is
refined as such, with a warning. The seeds and ANI of a `derep-refine` live
in its own stage record. `derep-stock` stores the `dereplicate` record and
`refined.tsv`, so after an unpack the file alone describes the refinement.

The lineage genomes are listed in `derep/refined.tsv`, and the `refined`
column of `cluster_summary.tsv` marks the clusters the refinement produced
(see [output.md](output.md#dereplicate-and-the-dereplication-commands)). The
stage record holds the seeds, the ANI of each lineage (`refine_lineages`)
and the counts.

#### The sourmash dereplicator and its sketch cache

`--tool sourmash` reads the genome sketches of the working directory (see
[Genome sketches](#genome-sketches)), so `--target-reps` compares the same
sketches at each search step. Without them (other `ksize` or `scaled`, or
genomes outside the set), it keeps its own signatures in a sketch cache.
`--target-reps` shares one cache across its search steps, and
`--tool-arg sketch_cache=DIR` sets a directory that persists across runs.

With the branchwater plugin, cache entries are matched by genome file name,
size and modification time. Without it, a per-genome signature is reused
when it is not older than the genome file. A genome replaced under the same
name is therefore sketched again, except when the new file carries an older
timestamp (copied with `cp -p` or `rsync -t`) or, with the plugin, the same
size and timestamp. Empty the directory after such a replacement.

Both back-ends, the branchwater plugin's `pairwise` and `sourmash compare
--avg-containment` without it, estimate ANI as c^(1/k) from the average
containment c of a pair (see
[choosing-tools.md](choosing-tools.md#how-the-tools-estimate-ani)). The
stage record holds `tool_settings`: the back-end used, that estimate, `ksize`
and `scaled` (see [output.md](output.md#files-shared-across-stages)).

### Chunked dereplication

`dereplicate -s N` (`--process-size`) dereplicates in chunks. Use it for
sets too large for one pass, or sets with many near-identical genomes. It
splits the genomes into chunks of at most N and dereplicates each chunk (in
parallel, `--num-processes`) with the stage-1 thresholds
(`--pre-primary-ani`, `--pre-secondary-ani`, which default to the final
ones). It then dereplicates the union of the chunk representatives with the
final thresholds. Every genome is assigned to the final representative of
its chunk representative.

The chunk size bounds the leaf chunks, not the final pass. A chunk can only
remove genomes that have a close relative in the same chunk, so for a
diverse set at a strict threshold the union is most of the input, and the
final pass compares nearly all genomes with each other. When a level removes
enough, one chunk of a further level is first run as a probe, and the level
runs only when the probe suggests it will pay; the rule is in
[developing.md](developing.md#scaling). The log names the size of the final
pass and what to expect of it, and a warning follows when the final pass is
more than four times the chunk size. A measured example on 6888 genomes is
in [choosing-tools.md](choosing-tools.md#by-number-of-genomes).

Three settings make the final pass smaller or faster:

- **Keep the default skani screen with skDER.** For a cutoff at or above 99
  percent the adapter passes skani a screen of the cutoff minus one, so
  pairs whose k-mer estimate is below that are not aligned. skDER's own
  default is the cutoff minus ten, which at 99.95 aligns every pair above
  about 90 percent. `--tool-arg screen=PERCENT` sets another value, and
  `screen=skder` restores skDER's own; a value above the cutoff is refused.
  The screen is an option of `skani triangle`, which the greedy and dynamic
  modes run; `low_mem_greedy` does not read it.
- **Raise `-s`** so that each chunk holds more close relatives. Chunking
  pays when the set has many near-identical genomes. For a set where most
  genomes are their own representative, one pass without `-s` costs about
  the same.
- **Set `--pre-secondary-ani` below `--secondary-ani`** (for example 0.999
  against 0.9995) so the chunks remove more. The chunk representatives are
  then compared at the final threshold, so clusters can be coarser than in a
  single pass at the final threshold.

The Nextflow scatter (`--derep_process_size`) runs one chunk round and one
merge task over the whole union, with the same semantics as the first level
and the final pass here (see [nextflow.md](nextflow.md#scaling)).

### Inspecting a dereplication

Three commands read a dereplicated working directory without rerunning the
dereplicator, and a fourth inspects the genomes before dereplication:

| Command | What it does |
|---------|--------------|
| `repgenr cluster-summary` | Regenerates `derep/cluster_summary.tsv`, one row per representative (see [output.md](output.md#cluster-summary)). |
| `repgenr derep-unpack` | Lays the clusters out as one directory per representative with its members inside. |
| `repgenr derep-stock` | Stores, restores, lists and deletes dereplication results under `derep/stock/`. |
| `repgenr glance` | Compares all genomes in `genomes/` and plots the pairwise values; it does not need a dereplication. |

The commands read `derep/clusters.tsv` and the manifest. They do not check
that the last `dereplicate` finished, so run `repgenr status -wd WD` first:
a failed or refused `dereplicate` shows as `[interrupted]` while `derep/`
still holds the previous result.

| Question | Where to look |
|----------|---------------|
| How many clusters, and how large is each? | `derep/cluster_summary.tsv` (`repgenr cluster-summary` rebuilds it); `n_members` does not count the representative. |
| Which genomes are in the cluster of representative X? | Rows of `derep/clusters.tsv` whose first column is X, or the directory `derep/unpacked/<X without its extension>/` after `repgenr derep-unpack`. |
| What happened to one genome? | `derep/genome_status.tsv`: `representative`, `contained` or `fail_qc`. |
| Do the genomes fall into clear groups before thresholds are chosen? | `repgenr glance`, see below. |
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

`repgenr derep-unpack` gives the same grouping as directories, for tools
that take a directory as input. The directory name is the representative's
filename without its extension. `--no-representant` leaves the
representative out. A member missing from `genomes/` is left out with a
warning that names it.

The files are hard links to the genomes where the file system allows it, so
editing one in place also changes the genome under `genomes/` (and, after a
default symlinking `ingest`, the original file). Elsewhere, for example on
exFAT or across volumes, they are copies, which take longer to write.

The directories describe the dereplication they were made from. After a new
`dereplicate` or a `derep-stock --action unpack`, run `derep-unpack` again;
it reruns because `clusters.tsv` changed, and `doctor` warns until then.

#### Storing and restoring a dereplication

`repgenr derep-stock --action pack --name <run>` stores the current
clusters, statuses, representatives and the completed `dereplicate` record
(tool, parameters, tool versions, as `record.json`) under
`derep/stock/<run>`. A run already stored under that name is replaced, with
a warning. `--action list` prints the stored run names on stdout, one per
line in name order, and `--action delete` removes one. Deleting a run that
is not stored exits 3 and lists the stored runs.

`derep-stock --action unpack --name <run>` makes a stored run the current
one again. It refreshes the manifest, rebuilds `cluster_summary.tsv` and
re-stamps the `dereplicate` record, so `status` and `versions` name the tool
that produced the stored run. What unpack replaces, and how it treats a run
stored without `record.json`, is described in
[output.md](output.md#storing-and-restoring-with-derep-stock).

The restored `dereplicate` record holds no resume fingerprint, since unpack
cannot reproduce the conditions of the original run. The next `dereplicate`,
also inside `repgenr run`, therefore recomputes the dereplication. Run
`phylo` directly to build on the restored set. `status` and `doctor` note
such a record.

#### Reading `glance`

`repgenr glance` compares all genomes in `genomes/` (at least two) and
writes a dendrogram and two plots. Use it before dereplication to see
whether the genomes fall into clear groups; it does not need a
dereplication and sets no thresholds. `--tool` chooses the comparison tool:

- `auto` (the default) uses dRep when it can run, on the `PATH` or through
  the container back-end, and sourmash otherwise. The log names the tool it
  picked, and the stage record holds that tool, not `auto`. With neither
  available, glance exits 4 and names the tools that can compare.
- `drep` runs `dRep compare` (Mash, primary clustering only) and plots Mash
  ANI. The dendrogram is dRep's own.
- `sourmash` compares the genomes with the parameters `dereplicate --tool
  sourmash` uses (k=31, scaled=1000), reading the genome sketches of the
  working directory and writing the missing ones first. It runs `sourmash
  compare --avg-containment` and converts the values to the ANI estimate
  the dereplication threshold applies to, c^(1/k). The dendrogram is
  average-linkage clustering on 1 - ANI. The full N x N matrix is held in
  memory, so sets above 5000 genomes are refused.

The plot axes name the measure, Mash ANI or ANI; the file names stay the
same for both tools. The two estimates differ slightly, so compare plots
made with the same tool.

The histogram shows how the pairwise values are spread. When they fall into
separate groups with an empty gap between them, a threshold placed in the
gap separates them. If the values are spread evenly across the range, the
threshold decides the cluster sizes, so check `cluster_summary.tsv`
afterwards.

With `--tool sourmash` and `--keep-files`, `glance_wd/pairwise_ani.csv`
holds the value of every pair, and `glance_wd/dendrogram_leaves.txt` the
dendrogram's leaf order. An empty cell means the pair is below the sketch's
reach: its sketches share fewer than 30 hashes and under a tenth of the
smaller sketch. dRep writes into its cache with names beginning `._` on
exFAT volumes and fails there; use an APFS or ext4 working directory (see
[install.md](install.md#2-satellites-on-path-after-core)).

#### Comparing two dereplications

Pack each result under its own name, then compare the stored files. Check
with `status` that `dereplicate` finished before each pack; each stored
run's `record.json` then names the tool and parameters that produced it.

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

Compare the cluster sizes first, then the memberships. Two dereplicators
can cluster the genomes identically and still keep different genomes as
representatives (see [choosing-tools.md](choosing-tools.md#dereplicators)),
so a different representative list alone does not mean a different
partition. To compare partitions, label every genome with the smallest
member name of its cluster and diff the labels:

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
`tree2tax`, and the stateless step `phylo-build` does the same outside a
working directory. A tree is built from an MSA (multiple sequence alignment)
or directly from the genomes:

| Route | Option | Suits |
|-------|--------|-------|
| Whole-genome alignment | `--aligner` (for example `progressivemauve`) | small sets, under about 30 genomes |
| Core-SNP alignment | `--msa-source snptype` | one species or a clonal set |
| Concatenated core-gene alignment | `--msa-source coregene`, after `repgenr annotate` | several species of a genus, or wider with `--coregene markers` |
| Alignment-free tree | a tree builder that works from the genomes, such as `--treebuilder mashtree` or `skani` | any set, as a quick look |

```bash
# Whole-genome alignment
repgenr phylo -wd $WD --aligner progressivemauve --treebuilder iqtree

# Core-SNP alignment: phylo types the genomes itself under tree/msa/
repgenr phylo -wd $WD --msa-source snptype --snptyper simple \
  --core-fraction 0.95 --mask gubbins --treebuilder iqtree -B 1000

# Core-gene alignment: annotate first
repgenr annotate -wd $WD
repgenr phylo -wd $WD --msa-source coregene --coregene ppanggolin \
  --treebuilder iqtree -B 1000

# Alignment-free tree
repgenr phylo -wd $WD --treebuilder mashtree
```

`repgenr phylo-advise -wd $WD` reads the working directory and prints the
route, the reasons and the `phylo` command with a `--preset` that sets the
route and its settings. A later `phylo` run that changes only the tree
builder reuses the alignment. Without `-B` or an SH-aLRT test an IQ-TREE
tree carries no supports. [phylogeny.md](phylogeny.md) describes the
presets, the SNP typers and masking, annotation, the core-gene builders,
support values, rooting and alignment reuse, and
[choosing-tools.md](choosing-tools.md#phylogeny-routes) explains which route
and tree builder suit a set.

## tree2tax

`repgenr tree2tax` turns the rooted tree into a FlexTaxD-compatible
taxonomy. It is the last stage and writes `tree2tax.tsv` (the relations) and
`genomes_map.tsv`, with `tree2tax_nodes.tsv` beside them (see
[output.md](output.md#tree2tax)). By default (`--include-dereplicated`) the
genomes that dereplication left out are listed under their representative.
The `tree2tax-relations` step does the same outside a working directory and
is the one the Nextflow module runs.

```bash
repgenr tree2tax -wd $WD --include-dereplicated
```

### Rooting in tree2tax

`tree2tax` and `tree2tax-relations` root the tree on the outgroup when it is
not rooted there yet. On a tree from `phylo --root outgroup` (see
[Rooting the tree](phylogeny.md#rooting-the-tree)) this changes nothing, and
the relations are the same as when `tree2tax` rooted the builder's output
itself. They exit 3 when the outgroup is not a leaf of the tree.

- After `phylo --root midpoint`, `tree2tax` keeps the midpoint root, since
  it reads the `root` that phylo recorded. `tree2tax-relations` does so when
  given `--keep-root`. The root is named by `--root-name` in either case.
  `--remove-outgroup` still removes the outgroup leaf. It is pruned before
  the collapse and naming steps, so no clade is named for a leaf set that
  held it, and a node left with a single child is merged into its parent.
- After `phylo --no-outgroup` with `--root none`, `tree2tax` leaves the tree
  unrooted, and its relations start at the builder's root.
  `tree2tax-relations` does so when given `--no-outgroup`, which ignores
  `--outgroup-dir` and `--outgroup-accession`.

The Nextflow pipeline passes `--no-outgroup` to `tree2tax-relations` when
`phylo_args` contains it, and `--keep-root` when `phylo_args` contains
`--root midpoint`.

### Collapsing weak splits in tree2tax

`repgenr tree2tax` (and the `tree2tax-relations` step) can merge weakly
supported or near-zero-length internal nodes into their parents before
naming them, so a split the data do not support does not become its own
FlexTaxD node. All thresholds are off by default.

**Units.** Lengths and diameters are in the tree's own units:

- mash distance for mashtree (about 0.05 at 95 % ANI for closely related
  genomes);
- one minus ANI as a fraction for the sourmash builder (0.05 at 95 % ANI),
  the ANI being the estimate c^(1/k) from the average containment, as in
  dereplication. A pair with fewer than 30 shared hashes and under a tenth of
  the smaller sketch has an unknown ANI and sits at distance 1.0, so a
  diameter does not reach between such pairs. For genomes of 1 to 1.5 Mb
  that reach ends near 88 to 90 % sketch ANI at k=31 and 83 to 85 % at k=21.
  With `metric=jaccard` the unit is one minus the Jaccard similarity at the
  chosen k (about 0.88 at 95 % ANI with k=31, so a mash-scale value is far
  too strict there);
- one minus ANI as a fraction for skani (0.05 at 95 % ANI);
- for the ska2 builder, the SNP count with `metric=snps` (the default) or
  the fraction of mismatching split k-mers with `metric=mismatch`;
- substitutions per site for the ML builders, and per variable site on an
  ML tree from a SNP alignment.

**Branch length and support.** Either criterion suffices.

- `--collapse-length L` merges a node whose branch is shorter than L.
- `--collapse-support S` merges a node whose support is below the fraction
  S. IQ-TREE (with `-B`, or with `--aligner-arg alrt=N` for SH-aLRT
  supports) and RAxML-NG write percentages, which are normalised. With both
  `-B` and `alrt=` IQ-TREE writes each label as `SH-aLRT/UFBoot` (for
  example `87.3/98`). The last field, the ultrafast bootstrap, is the one
  read, so S applies to the same value as with `-B` alone. FastTree and
  aster write fractions. mashtree, the sourmash builder, skani and the ska2
  builder write no supports, and the option then warns and does nothing.

**Clade diameter.** `--collapse-diameter D` works on clades rather than
single branches. The fewest clades whose leaf-to-leaf distance stays within
D are kept, and every node inside such a clade is merged into it (the
`max_clade` method of TreeCluster). A group of near-identical genomes joined
by many short branches therefore becomes one FlexTaxD node with the genomes
as its children. The diameter rule runs before the other two rules and
refuses a tree without branch lengths. A clade it kept can still be merged
into its parent by `--collapse-length` when its own branch is short, which
then puts leaves more than D apart under one node. An ingroup that fits
within D on an outgroup-rooted tree becomes one node under the root.

**Common rules.** The root, the outgroup/ingroup split under it and the
leaves never collapse. A collapsed node's branch length is added to its
children's. Provenance records the thresholds and the number of nodes
collapsed, so changing one reruns the stage. `tree2tax_nodes.tsv` beside the
two tables describes every node (leaf count, depth, branch length, support,
relative evolutionary divergence, the lineage of each genome) and is the
place to look when choosing a threshold.

### Node names

Internal nodes are named by a hash of their leaf set unless
`--node-basename` numbers them. A leaf set changes whenever a genome joins
or leaves it, so every ancestor of the change is renamed on the next run.

`--naming lca` names a node by a hash of the smallest leaf label under each
of its children instead. The node is the last common ancestor of those
leaves, and the name is unique per node. Names then change less often, but
still change:

- a genome whose label sorts after the present ones (string order, which is
  not always chronological across GCA, GCF and viral prefixes) and attaches
  by splitting a branch leaves every existing name unchanged;
- one that joins a multifurcation (which the collapse options produce) as a
  further child renames that node;
- one that sorts below a child's current smallest label renames the
  ancestors where it becomes the minimum;
- removing a genome renames the nodes whose minimum or child it was.

The tree is also inferred anew on each run, so a topology change elsewhere
renames the nodes it touches. Names are opaque identifiers in either case.
The record carries the rule, or null when a basename numbered the nodes.

### A rank column

`tree2tax.tsv` has two columns by default, since FlexTaxD reads a third as
the rank of every row and an empty value inherits the previous row's.
`--ranks` (on `tree2tax`, `tree2tax-relations` and `run`) adds the column
with `strain` on every genome and `no rank` on every internal node, so a
FlexTaxD database built from it carries a rank column with no taxonomic rank
asserted.

`--rank-red species=0.95,genus=0.8` names internal nodes by their relative
evolutionary divergence (RED, the `red` column of `tree2tax_nodes.tsv`). A
node takes the rank with the largest threshold not above its RED, and `no
rank` when an ancestor already holds that rank. Each named rank therefore
appears at most once on a root-to-leaf path and in threshold order. A
lineage may lack a rank when no ancestor reaches that band.

The thresholds are the user's, and the example is illustrative. GTDB's
median RED values (Parks et al. 2018; about 0.92 for genus and 0.77 for
family, with species delimited by ANI rather than RED) come from a
bacteria-wide tree. They transfer to a genus-level mashtree or SNP tree only
loosely, so inspect `tree2tax_nodes.tsv` before settling on values.

`--rank-red` needs branch lengths and a rooted tree (an outgroup, `phylo
--root midpoint`, or `--keep-root` on the step) and is refused otherwise.
`run` checks the bands, `--ranks` and the rooting before any stage runs. The
record stores the bands sorted by threshold and carries `ranks` and
`rank_red`, so changing them reruns the stage.

## Running tools in containers

RepGenR2 can run each external tool inside a pinned container instead of
relying on tools installed on `PATH`. Use it to run Linux-only tools on any
host (macOS, for example) and to make tool versions reproducible. RepGenR2
itself runs on the host and only the tool subprocess is containerised, so
it works in the plain CLI and inside Nextflow.

```bash
# Run every tool in a container via Docker:
repgenr --container docker dereplicate -wd $WD --tool drep

# Singularity/Apptainer (HPC), with .sif images on a large disk:
repgenr --container singularity --container-cache /Volumes/LaCie/repgenr_sif \
        phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
```

`--container` is a top-level option, placed before the subcommand. The
default is `none` (native execution). The other container options and their
environment variables, the image sources, the storage location and the
setup on Apple Silicon are in [install.md](install.md#3-containers). The
Nextflow container profiles are described in
[nextflow.md](nextflow.md#profiles-and-containers).

### Per-tool directories (`--bin-dir`)

`--bin-dir TOOL=DIR` puts DIR first on `PATH` for one tool only. Use it when
a tool lives in its own conda environment (a satellite environment) rather
than in the environment that runs repgenr. DIR comes first for the tool's
preflight lookup and version query, and for every host subprocess it runs,
so a wrapper such as `run_gubbins.py` calls the helpers of its own
environment. Other tools keep the inherited `PATH`.

```bash
P=$(conda info --base)/envs
repgenr --bin-dir gubbins=$P/repgenr-gubbins/bin --bin-dir mashtree=$P/repgenr-mashtree/bin \
        phylo -wd $WD --msa-source snptype --snptyper simple --mask gubbins
```

TOOL is an adapter name as printed by `repgenr list-tools` (for example
`gubbins`, `mashtree`, `snippy`, `parsnp`, `progressivemauve`), or
`checkm2`, `datasets` or `minimap2` (racon's mapper). The option is
repeatable, and `REPGENR_BIN_DIRS="gubbins=/p/bin,mashtree=/q/bin"` sets the
same; an option given on the command line replaces the variable's entry for
that tool. An unknown tool or a directory that does not exist exits 2.

The directory each tool used is recorded under `bin_dirs` in the stage's
record in `repgenr.yaml`. It is not part of the resume fingerprint, so
setting or changing it does not rerun a finished stage; the tool's version
is recorded as usual. Under `--container`, a tool that runs in an image is
not affected, and a warning names it. In Nextflow, the option is passed
through `--repgenr_opts` (see
[nextflow.md](nextflow.md#profiles-and-containers)).

## Troubleshooting

This part lists common failures and how to diagnose them, then the exit
codes.

### A tool is not found (`MissingBinaryError`)

The Python package does not install the bioinformatics tools (see
[install.md](install.md)). Put the tool on `PATH` from a conda environment,
or use `--container docker` (or `singularity`) to run tools in pinned images
instead (see [Running tools in containers](#running-tools-in-containers)).
Run `repgenr list-tools` to see the adapters and each tool's declared genome
limit; `list-tools --check` also runs every adapter's preflight and reports
which binaries are missing or too old.

### Apple Silicon / arm64

BioContainers are amd64. Pass `--platform linux/amd64` and enable Rosetta so
emulated images run (see [install.md](install.md#3-containers)).

### A stage failed; where are the details?

Errors print a concise message, and the full traceback is in
`<workdir>/repgenr.log`. A failed external tool prints one line naming the
tool and its exit status. The command line and the output tail are in the
same log. Rerun with `--verbose` to see them on the console; a stateless
step has no log and always prints the tail.

Tool command lines are shortened on the console to fit 120 columns (the
rule is in [developing.md](developing.md#running-the-external-tool)). The
full text of such a line, the container engine options and the tool's own
output are in `repgenr.log`.

`repgenr status -wd <WD>` shows what completed and what is next. A stage
that failed is listed as `[interrupted]`, and its outputs may be partial
until it is rerun. The other states are described under
[Checking a working directory](#checking-a-working-directory).

### GTDB download fails

Check `--release` (for example `232.0`) and `--gtdb-version`
(`bac120`/`ar53`). Transient HTTP errors are retried automatically. A host
that does not accept a connection within 15 s counts as unreachable, so a
blocked network exits 3 after about two minutes of retries. The `--source
api` mode fetches only the target taxon, without a full-table download.

### Genomes NCBI no longer serves

A GTDB release can list assemblies that NCBI has since suppressed. `genome`
records them in `missing_accessions.txt` (also when a whole download batch
consists of them), and later stages excuse them; a rerun asks for them
again. Each rehydrated genome is checked against the package's
`md5sum.txt`, and a mismatch is recorded the same way. An outgroup NCBI
does not serve exits 3; choose another with `metadata --outgroup-accession`.

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
limit. An HTTP error is retried per batch of taxids. A connection error (no
network, or the host does not answer) stops the lookup at the first batch
with exit 3.

### BV-BRC behind a proxy

`vmetadata --source bvbrc` reaches BV-BRC over FTPS (`ftp.bvbrc.org`), which
does not use the `HTTPS_PROXY` and `HTTP_PROXY` settings. On a network that
allows only proxied traffic, use the default NCBI Virus source.

### ParSNP on three genomes

ParSNP's internal RAxML step refuses fewer than four genomes, so
`--snptyper parsnp` exits 6 on a set of three. The three-genome minimum that
`phylo` checks before any tool runs does not cover this case.

### A tool hangs

Set `REPGENR_SUBPROCESS_TIMEOUT=<seconds>` to cap every external tool. On
expiry the process group is killed and the stage fails with an error that
names the time limit. Under `--container docker` the tool's container is
stopped as well.

### Stopping a run

Each external tool runs in its own process group, together with the helpers
it starts (for example IQ-TREE under `run_gubbins.py`). On SIGTERM (`kill`,
a scheduler), SIGHUP (the terminal closes) or Ctrl-C, repgenr sends SIGTERM
to every running tool group and starts no queued tool. It sends SIGKILL to
the groups still present after 5 s and exits: 143 for SIGTERM, 129 for
SIGHUP, and the usual `KeyboardInterrupt` for Ctrl-C. The stage is then listed as
`[interrupted]` and a partial output is removed. A second signal kills the
tools at once.

A signal repgenr starts with as ignored stays ignored, so `nohup repgenr ...
&` survives the terminal closing. Ctrl-Z suspends the tools with repgenr,
and `fg` or `bg` resumes them. Without a terminal (`setsid`, a workflow
manager) SIGTSTP has no effect, as for any process in that situation.
SIGKILL to repgenr cannot be handled, so its tools keep running. Each tool's
process group ID is its own PID, so `pkill -g <tool pid>` stops a tool with
its helpers.

**Under `--container docker`.** When repgenr is stopped it also stops the
container of the running tool, so a tool that ignores SIGTERM does not keep
running (details in [developing.md](developing.md#running-the-external-tool)).
SIGKILL to repgenr leaves its containers running. Each container carries the
label `repgenr.pid=<repgenr pid>`, so they can be stopped by label:

```bash
docker ps -q --filter label=repgenr.pid=<repgenr pid> | xargs -r docker stop
```

### Working directories from earlier versions

Three changes affect a working directory written by an earlier version:

- **Each stage reruns once**, since the format of the resume fingerprint
  changed. The manifest is upgraded in place to the current schema when a
  stage opens it, and an earlier version then refuses to open it.
- **The GTDB representative flag is missing** from a selection made before
  the flag was recorded, so `--keeper gtdb` acts as `--keeper quality`.
  `repgenr metadata --force` with the same options records it. `genome` then
  runs again because `selection.tsv` changed, and downloads nothing already
  present; `dereplicate` runs again because the manifest now flags genomes.
- **The viral species of an NCBI Virus working directory** came from the
  organism name before it came from NCBI Taxonomy. Such a working directory
  keeps its genomes and filenames while its stages skip. When `vgenome`
  reruns on it (new arguments or `--force`), it refuses the old records and
  asks for one `repgenr --force vmetadata` with the same arguments. That
  download renames the viral genomes, and `vgenome` and the later stages
  then rerun.

### Exit codes

A script can tell the failure classes apart without reading the log:

| Code | Meaning |
|------|---------|
| 0 | Success. |
| 1 | An unexpected error (traceback in the run log). |
| 2 | Invalid or missing user input (also Typer's own usage errors). |
| 3 | The working directory does not exist (every command, including `status` and `doctor`) or is missing files or in a bad state, or a request to a remote service (GTDB, NCBI Entrez, NCBI Datasets, BV-BRC, ENA) failed, for example because the network is unreachable. `assemble` and `reads-gather` also exit 3 when every run was excused; see [output.md](output.md#when-every-sequencing-run-is-excused). |
| 4 | A required external tool is absent or below its version floor, or the Docker daemon cannot be reached under `--container docker`. `list-tools --check --strict` exits 4 when any adapter is missing or errored. |
| 5 | A requested tool adapter could not be found or loaded. `list-tools --check --strict` exits 5 when any plugin failed to load (this takes precedence over 4). |
| 6 | An external tool failed (one console line; the command and output tail are in `repgenr.log`). A network failure inside the `datasets` CLI is reported as 6. Under `REPGENR_PROPAGATE_TOOL_EXIT=1` (set by the Nextflow modules) the tool's own status is forwarded instead, a signal kill as 128 plus the signal number. |
| 7 | `doctor` found at least one failure: an output that disagrees with its record, a malformed `repgenr.yaml`, or a check that could not complete. Warnings alone exit 0. |
