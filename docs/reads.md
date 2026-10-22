# Genome sets from sequencing runs

This page describes how RepGenR2 builds a genome set from whole-genome
sequencing runs in ENA and SRA: selecting runs, assembling and polishing
them, scoring and classifying the assemblies, and screening runs before or
during assembly. Dereplication, phylogeny and tree2tax then run as for any
other entry path ([usage.md](usage.md)).

## Contents

1. [Overview](#overview)
2. [Selecting runs](#selecting-runs)
3. [Assembly](#assembly) and [polishing long-read assemblies](#polishing-long-read-assemblies)
4. [Quality scoring and classification](#quality-scoring-and-classification)
5. [Screening runs](#screening-runs): the three checks, before download,
   during assembly
6. [Adding assemblies to a GTDB selection](#adding-assemblies-to-a-gtdb-selection)
   and [assembling in parts](#assembling-in-parts)
7. [Stateless steps](#stateless-steps)

## Overview

Use this route for a taxon whose public isolates are available as reads but
not as assemblies, or to add such isolates to a GTDB selection. `repgenr
reads` selects runs from ENA (which mirrors SRA) and writes `reads.tsv`.
`repgenr assemble` fetches and assembles them into `genomes/` with the same
`selection.tsv` and manifest the other entry paths write. The chain is
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
`--assembler`, `--threads` and `--outgroup`. The other options are given on
the stage commands. Which assembler and polisher suit each sequencing
platform is in [choosing-tools.md](choosing-tools.md#assemblers-and-polishers).
The files each step writes are described in [output.md](output.md#assemble).

## Selecting runs

Runs are chosen by taxon or by accession.

- **By taxon:** `--target-family`/`-tf`, `--target-genus`/`-tg` or
  `--target-species`/`-ts`, resolved through the ENA taxonomy, synonyms
  included. The number of runs under the taxon is logged before the search.
  A warning follows above 100000 runs: Salmonella, with about 880000, needs
  about 3 GB of memory for its records.
- **By accession:** `--accession` takes a run (SRR/ERR/DRR), a sample
  (SAMN.., SRS..) or a study (PRJNA.., SRP..) and repeats. `--accession-file`
  lists them one per line; text from a `#` to the end of the line is a
  comment. Runs found by accession pass the same filter as the taxon query,
  whole-genome sequencing of genomic DNA. RNA-Seq, amplicon or metagenomic
  runs of a named study or sample are dropped with a warning that names them.

The selection is then filtered:

- `--platform illumina|ont|pacbio` keeps one platform.
- `--min-bases` drops small runs, and `--max-bases` drops runs above a size.
  An unenriched whole-host library, tens of Gb for a 1 Mb endosymbiont,
  would assemble into a host-dominated genome.
- `--drop-selection` drops runs by ENA library selection. The default is
  `MDA` (whole-genome amplification), which assembles into chimeric and
  uneven contigs. Repeat the flag for more values; `--drop-selection none`
  keeps every run. The value is kept in `reads.tsv` as `library_selection`.
- `--one-per-sample` (the default) keeps the best run of each sample. That
  is a long-read run when it carries at least 100 Mb and a tenth of the
  sample's largest short-read run, else the largest run. `--all-runs` keeps
  every run.
- `--max-runs` caps the selection to the largest runs.

Each run is labelled with the family, genus and species of its NCBI taxid,
in the same filename tokens the GTDB path uses.

## Assembly

`assemble` downloads each run's FASTQ files from the locations ENA lists,
over HTTPS. The files are verified against ENA's checksums while the bytes
are written, so a file is not read a second time. `assemble` then assembles
the reads and keeps the contigs of at least `--min-contig-length` bases
(500), renamed `<run>_contig<n>`. The genome is named
`Family_genus_species_RUN.fasta` from the tokens the reads stage resolved.

**Assembler choice.** `--assembler auto` (the default) picks by platform and
layout: `skesa` for Illumina, `flye` for Oxford Nanopore and PacBio.
`shovill` (SPAdes) is the alternative for paired Illumina runs. `--tool-arg`
passes tuning such as `mode=nano-raw` to Flye.

**Layouts.**

- Long-read runs are assembled whatever layout ENA gives them (some ONT runs
  are labelled PAIRED). The polishers join a run listed as several FASTQ
  files into one.
- A paired run that ENA lists with a third, orphan FASTQ file is given to
  skesa as the pair plus the orphan file, and to shovill as the pair only.
- A short-read run labelled PAIRED but listed with a single FASTQ file
  (mate-1 reads only, for example) is planned as single-end. `auto`
  assembles it with SKESA. `--assembler shovill`, which needs a read pair,
  excuses it as `unsupported_layout` before anything is downloaded. The log
  names each such run.

**Parallel runs and memory.** `--jobs` runs that many assemblies at once,
with `--threads` split across them. Memory, not CPU, is the limit, so the
default is 2, or 1 as soon as a long-read run is pending. `--memory-gb` is
the RAM cap passed to SKESA and shovill (16; shovill accepts no less than
8).

**Resume.** Each finished run leaves a marker under `assemblies/<run>/`, so
an interrupted stage resumes without refetching. Reads are deleted after a
successful assembly unless `--keep-reads` is given, and the assembler's
scratch unless `--keep-files` is given. A run interrupted during its
assembly keeps the FASTQ files that still match their checksum, so the
resume does not download them again.

The marker records the settings the run was built with: `--assembler`,
`--polisher`, `--polish-rounds`, `--min-contig-length` and the `--tool-arg`
keys its tools read. A later call, `--force` included, reuses a finished run
only while these agree. A higher `--min-contig-length` filters the finished
contigs again, which gives the same result as filtering the raw assembly.
Any other change assembles the run again, and the log names the change. To
assemble a run again under the same settings, delete `assemblies/<run>/`. A
marker that records no settings is reused as it is, filtered at the
requested floor.

**Excused runs.** The following runs are written to `excused_runs.tsv` with
the reason, and the rest proceed. The completeness guard of later stages
excuses those runs.

- a run without an ENA FASTQ mirror;
- a run whose download fails its checksum;
- a run no assembler accepts (`unsupported_platform`);
- a run whose assembler is not installed under `--assembler auto`
  (`assembler_not_installed`, with a warning naming the adapters that would
  take it);
- a run whose assembly fails.

When no run can be assembled because no assembler is installed, the stage
exits 4 instead. A repeat with the same settings skips the stage, so a run
whose download failed is retried only with `--force`. The closing warning
names such runs. What happens when every run is excused is described in
[output.md](output.md#when-every-sequencing-run-is-excused).

**Outgroup and statistics.** `--outgroup FASTA` sets a genome aside for
rooting, as `ingest --outgroup` does. Per-assembly metrics (contigs, total
length, N50, coverage from the sequenced bases) are in `assembly_stats.tsv`.

## Polishing long-read assemblies

Long-read assemblies are polished with the run's own reads before the
contig filter. `--polisher auto` (the default) runs medaka for ONT runs
(racon when medaka is not installed). It runs racon (minimap2 overlaps,
`--polish-rounds` rounds) for PacBio CLR runs, and nothing for PacBio HiFi or
Illumina. `--polisher none` turns polishing off.

medaka needs the basecaller model. Reads basecalled with Dorado name it in
their FASTQ headers, and medaka resolves it from there. SRA rewrites the
headers of the reads it mirrors, so for those reads give the model with
`--tool-arg model=...`. Otherwise the default applies: ONT's bacterial
methylation-aware model for R10.4.1 400 bps chemistry, the usual case for
public bacterial ONT runs since 2023. `--tool-arg bacteria=false` uses
medaka's general default model instead. The marker records which of the
three applied.

A polishing failure excuses the run with `polish_failed`, and
`assembly_stats.tsv` names the polisher per genome. When `--polisher auto`
finds an adapter for a run but its tool is not installed, the run is
assembled unpolished and one warning per platform names the adapters.
Unpolished ONT assemblies carry indel errors that break genes, which CheckM2
reads as lower completeness and higher contamination.

## Quality scoring and classification

Two optional checks run on the assemblies: CheckM2 quality scoring and GTDB
classification with sourmash. Without a database the checks are skipped and
the log says so.

**CheckM2.** With a CheckM2 database (`--checkm2-db`, or the `CHECKM2DB`
variable CheckM2 itself reads; obtain it with `checkm2 database
--download`), every assembly is scored. The completeness and contamination
reach `selection.tsv` and the manifest, so `--keeper quality` works as it
does for GTDB genomes. An assembly below `--min-completeness` (50) or above
`--max-contamination` (10) is excused with `qc_failed`. An assembly for
which CheckM2 reports no result is kept with a warning and without quality
values.

The scores are stored per run in `assemblies/<run>/checkm2.json`, with the
contigs' SHA-256, the database's resolved path and size, and the CheckM2
version (the image under `--container`). A later call that changes only
`--min-completeness` or `--max-contamination` applies the stored scores
without running CheckM2 again. A run whose contigs, database (path, size or
modification time) or CheckM2 version differ is scored again. `--force`
reruns the stage but does not bypass the stored scores. To score a run
again, delete `assemblies/<run>/checkm2.json`.

**GTDB classification.** With a GTDB sourmash sketch and its lineages CSV
(`--gtdb-sketch` and `--gtdb-lineages`, or `REPGENR_GTDB_SKETCH` and
`REPGENR_GTDB_LINEAGES`), each assembly is classified by `sourmash gather`.
The `gtdb-rs226-reps.k31-sc10k.sig.zip` sketch and its `lineages.csv` from
`https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/` serve.
`--classifier auto` runs the classifier when a sketch is configured, and
`none` never does.

When the GTDB genus agrees with the submitted organism, the GTDB family,
genus and species name the genome file, so a reads-derived genome groups
with GTDB-downloaded ones. Otherwise the submitted name stays and
`assembly_stats.tsv` flags the genome:

- `genus_renamed` when the genus differs but the family and the species
  epithet agree. GTDB has then moved the species to another genus of the
  same family (NCBI Mycoplasmopsis arginini is GTDB Metamycoplasma
  arginini).
- `classifier_disagrees` for any other difference, including a shared
  epithet in another family (Klebsiella pneumoniae against Streptococcus
  pneumoniae).

Both are warned about, neither excuses the genome, and only
`classifier_disagrees` counts as a disagreement in the stage record. Both
lineages are kept in that table, and the sketch release is recorded in
provenance next to the metadata release.

**Preflight and memory.** The database paths and the checkm2 and sourmash
binaries are checked before any run is fetched. A wrong path therefore
exits 2 and a missing tool exits 4 at once, rather than after the
assemblies. sourmash runs one gather per assembly, and with the GTDB rs226
representatives sketch each gather holds about 0.6 GB of memory. As many run
at once as both `--threads` and `--memory-gb` / 0.6 allow, and at least one
runs (16 GB, the default, allows 26). The log names the number chosen.
`genome-qc` takes the same `--memory-gb`, and the Nextflow module passes the
task's memory.

## Screening runs

A screen judges each run's reads against GTDB before the run is assembled.
Use it to leave out runs of another taxon, runs whose reads match their
species poorly, and repeated sequencing of one isolate before the assembly
cost is paid. The screen reads a sourmash sketch of the run's reads, which
can be made in two places:

- **before any download**, by streaming the runs from NCBI into sketches
  (`reads-sketch`) and judging those sketches (`reads-screen`);
- **during assembly**, from the reads `assemble` has fetched
  (`assemble --screen-reads`).

### The three checks

Both forms apply the same checks with the same defaults. A run that fails
one is excused with its reason:

- `taxon_mismatch`: `sourmash gather` of the reads sketch (k=31) against the
  GTDB sketch, resolved to a lineage with `sourmash tax genome` and the
  lineages CSV (the commands of the classifier), gives a genus other than
  the submitted one. The classifier's tolerance applies: a genus that
  differs while the family and the species epithet agree (`genus_renamed`)
  passes. A run with no GTDB match above the gather threshold (50 kbp) is
  also excused with this reason. A run whose submitted genus is unknown is
  not tested.
- `low_match_fraction`: the abundance-weighted fraction of the reads
  assigned to the top species (`f_weighted_at_rank` of `tax genome` at
  species rank) is below the minimum fraction (default 0.3). The excuse
  reads `top species fraction X below --screen-min-fraction Y`.
- `duplicate_isolate`: the run's reads are contained in those of an earlier
  run of the same biosample, at an ANI estimate of the duplicate threshold
  (default 0.999) or more.

**The match fraction.** The fraction mostly measures the distance between
the run's strain and the GTDB representative of its species, since the
reference sketch holds one genome per species. Host, contaminant and
repeated sequencing-error k-mers lower it as well, so a low value does not
by itself show contamination. `tax genome` is run with a containment
threshold of 0, because most distinct k-mers of a reads sketch are
sequencing errors. The unweighted fraction of a pure isolate therefore
stays far below sourmash's default threshold of 0.1. Against the GTDB rs226
representatives sketch (k=31, scaled=10000), four Illumina MiSeq isolates
gave fractions of 0.39 to 0.58, and the default of 0.3 lies below that
range. A species far from its representative can fall lower, so check the
screen table before relying on the default.

**Duplicate detection.** Runs are compared in `reads.tsv` order, which the
reads stage writes largest first, so of two duplicates the larger run is
kept. Runs of different biosamples are never compared, because clonal
isolates of different samples are distinct samples. The estimate is taken
from the k=31 signature. How the containment is turned into an ANI estimate,
and its resolution, is described in
[choosing-tools.md](choosing-tools.md#how-the-tools-estimate-ani).

### Before download: reads-sketch and reads-screen

`reads-sketch` sketches the reads of the selected runs without a download.
Each run is streamed from NCBI with `fastq-dump -Z` and piped into one
`sourmash sketch dna` call, so no FASTQ file is written. `--max-spots N`
stops the stream after the first N spots of a run (a paired-end spot gives
two reads). The sketches have the parameters of the reads sketch `assemble`
writes (k=21, 31 and 51 at scaled=1000, with abundances), so the same
consumers read them.

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina
repgenr reads-sketch -wd $WD --max-spots 200000            # every run of reads.tsv
repgenr reads-sketch -wd $WD --accession SRR1 --accession SRR2   # a subset
repgenr reads-screen -wd $WD --max-spots 200000 --gtdb-sketch $DB --gtdb-lineages $LIN
repgenr reads-screen -wd $WD --max-spots 200000 --min-fraction 0.4 --accession SRR1
repgenr assemble -wd $WD --accession-file $WD/reads_screen/spots_200000/passed.txt
```

**Sketching.** The runs come from the working directory's `reads.tsv`, and
without a selection every run is sketched.

- `--accession` and `--accession-file` restrict the call to runs of that
  table, and `--platform` to one platform. An accession not in the table is
  refused by name, so a typo cannot sketch an unrelated run.
- A call adds to the set. A run sketched before at the same `--max-spots` is
  skipped, so a few runs can be sketched first and the rest later.
- Each cap is its own variant directory: `reads_sketches/all/` for whole
  runs and `reads_sketches/spots_<N>/` for `--max-spots N`. Sketches at
  several depths therefore live side by side. `--force` streams the selected
  runs of the current variant again.
- The index `reads_sketches/<variant>/reads_sketches.tsv` lists every run of
  `reads.tsv` with its status (`ok`, `failed`, `empty` or `pending`); see
  [output.md](output.md#reads-sketch). The stage is recorded as
  `reads_sketch`, and `status` lists the variants present.
- A failed run (an accession NCBI does not hold, a network error) is
  recorded and the others continue. A repeat of the same call streams the
  failed runs again and skips the rest, so a run that keeps failing is
  retried on every repeat until the selection is narrowed. The command exits
  3 and records nothing only when every selected run failed.
- A run whose stream holds no biological read is marked `empty`, has no
  sketch, and is not streamed again at that cap unless `--force` is given.
- `--jobs` streams that many runs at once (default 4, or the usable CPUs
  when fewer). Each job is a `fastq-dump` stream and one single-threaded
  sourmash process, so a job costs about one CPU.

Transient archive failures are retried within the call, and one variant
directory is written by one call at a time; the retry and lock rules are in
[developing.md](developing.md#resume-semantics).

The first N spots of a run are not a random sample. An Illumina file starts
with its first tiles, whose quality profile differs, and an ONT file with
the first reads of the run. For a species-level `gather` on an isolate this
is expected to matter little, but a low-abundance organism in a mixed sample
may be missed in a short prefix. A few hundred thousand spots is the
suggested starting point; it has not been tested as a threshold.

`reads-sketch` needs `fastq-dump` (sra-tools 3.0 or later, in the core
environment) and reaches NCBI, not the ENA mirror `assemble` downloads
from. A run that only ENA holds fails with its reason in the index. Under
`--container docker` each tool runs in its own image with the pipe between
them. sra-tools then runs from NCBI's own image (`ncbi/sra-tools`), which is
built for amd64 and arm64 and so runs natively on Apple Silicon.

**Screening.** `reads-screen` judges one variant of `reads_sketches/` with
the three checks, using the submitted taxonomy and biosample of each run
from `reads.tsv`. It needs sourmash and the GTDB sketch and lineages CSV
(`--gtdb-sketch` and `--gtdb-lineages`, or the environment variables). The
thresholds are `--min-fraction` and `--dup-ani`.

- `--max-spots N` names the variant to screen. Without it the `all/` variant
  is used when present, else the only variant there. With several variants
  the flag is required.
- `--accession` and `--accession-file` restrict the selection to runs of
  `reads.tsv`; without them every run is screened. A run without a sketch in
  the variant is listed as `screen_failed: no reads sketch`.
- A capped sketch holds the first N spots only, so the duplicate check
  compares the starts of two runs. With a cap that still covers the genome
  several times over, the containment of a true duplicate stays close to 1,
  but the estimate is less certain than from whole runs.
- `--jobs` bounds the gathers run at once (default 4, or the usable CPUs
  when fewer), each holding about 0.6 GB against the GTDB representatives
  sketch.

The decisions go to `reads_screen/<variant>/`: `reads_screen.tsv` with one
row per run of `reads.tsv`, one `screen.json` per run, and `passed.txt`,
which lists the accepted runs in `reads.tsv` order in the form `assemble
--accession-file` takes (see [output.md](output.md#reads-screen)). The stage
is recorded as `reads_screen`, and `status` lists the variants screened with
their counts.

A repeat with the same settings gathers nothing again. A new sketch of a run
or a changed threshold screens the run again from the kept sketch, without
streaming anything. A run the tools failed on (`screen_failed: gather`)
keeps no decision, so the repeat screens it again. `--force` screens every
selected run again.

The duplicate check is made again on every call, in `reads.tsv` order, over
every run with a decision made with the same settings, selected or not.
Several calls therefore give the result of one call over every run: a run
screened later can excuse an earlier run as its duplicate. A run decided
with other settings and not selected in the call keeps its row, takes no
part in the duplicate check and is left out of `passed.txt`, with a warning
that counts such runs. Screen every run with the same settings before
relying on `passed.txt`. `assemble --screen-reads` does not read these
decisions, since it screens again from its sketch of the whole run.

**Without a working directory.** `reads-sketch -o DIR` with bare run
accessions (SRR, ERR or DRR), or with `--reads-tsv` and the runs of a
`reads.tsv`, sketches into `DIR` and records nothing. `--platform` needs
`-wd` and is refused there. `reads-screen --reads-tsv reads.tsv
--sketches-dir DIR -o OUT` screens the sketches of such a directory into
`OUT/<variant>/` and records nothing. Both take `--versions-out`, which
writes the tool versions as a YAML fragment. `reads-sketch
--allow-all-failed` turns the error of a call where every run failed into a
warning, for a call that is one chunk of a larger batch. The Nextflow
`READS_SKETCH` and `READS_SCREEN` tasks use these forms (see
[nextflow.md](nextflow.md#reads-sketch-and-screen)).

### During assembly: assemble --screen-reads

While a run is assembled, `assemble` also sketches its reads: one `sourmash
sketch dna -p k=21,k=31,k=51,scaled=1000,abund` call on all FASTQ files of
the run, named by run accession and written to
`assemblies/<run>/reads.sig.zip`. The two files of a pair give one signature
per k-mer size. Abundances are kept because read k-mers also carry
coverage.

- The sketch starts once the FASTQ files are verified and runs in a thread
  beside the assembler. The assembler keeps its full share of `--threads`,
  so each run uses one thread more than its share while the sketch lasts.
- The sketch ends before the reads are deleted. It is written to a
  temporary file and renamed into place, so an interrupted sketch leaves no
  `reads.sig.zip`.
- `--reads-sketch/--no-reads-sketch` controls it. By default the reads are
  sketched when sourmash can run, and otherwise one INFO line per stage says
  so. `--reads-sketch` requires sourmash and exits 4 before anything is
  downloaded when it is missing.
- A failed sketch never fails the assembly. It is logged as a warning and
  the marker records no sketch. An excused run keeps no reads sketch, a
  quality-excused one included.
- A finished run whose marker records no sketch is not assembled again.
  When its FASTQ files were kept (`--keep-reads`) and still match their
  checksums, the next `assemble` call sketches them and updates the marker.
  Otherwise the sketch stays empty and one INFO line counts such runs. When
  the stage is otherwise up to date, that call needs `--force`.
- The reads sketch is not a genome sketch and is never copied into
  `sketches/` (see [usage.md](usage.md#genome-sketches)).

`assemble --screen-reads` (off by default) screens that sketch before the
run is assembled. The per-run order is: fetch and verify the FASTQ files,
sketch the reads, screen the sketch, and only then choose and start the
assembler. Under the screen the sketch is therefore awaited rather than run
beside the assembler. A run that fails a check is excused in
`excused_runs.tsv` with step `screen`. The thresholds are
`--screen-min-fraction` (0.3) and `--screen-dup-ani` (0.999).

- **Duplicates.** A run waits for the decisions on the earlier runs of its
  biosample, so the result does not depend on which worker finishes first.
  The comparison covers the runs accepted in the same call and finished runs
  with a reads sketch. A run that passes the taxon and fraction checks from a
  kept sketch is fetched before its duplicate check, which runs in the
  worker.
- **Requirements.** The screen needs sourmash and the GTDB sketch and
  lineages CSV, as for the classifier. Without a sketch or lineages file, or
  with `--no-reads-sketch`, the stage exits 2 before anything is
  downloaded, and without sourmash it exits 4. The screen runs with
  `--classifier none` as well, and the classifier runs unchanged on the
  accepted assemblies.
- **Tool failures.** A failed reads sketch or gather excuses the run with
  `screen_failed`. Like a failed download, it does not count as a judgement
  of the run's data when every run is excused (see
  [output.md](output.md#when-every-sequencing-run-is-excused)); the three
  check reasons do count.
- **Cost and outputs.** Each run costs one gather and one `tax genome` call,
  plus the time the sketch would otherwise have overlapped with the
  assembler. Gathers of concurrent runs are bounded by the classifier's
  memory budget (`--memory-gb` / 0.6 GB). An excused run keeps its reads
  sketch, its decision in `assemblies/<run>/screen.json` and the gather and
  `tax genome` output under `assemblies/<run>/screen/`, so the reason can be
  inspected; its reads are removed as usual. `screen_reads.tsv` in the
  working directory lists every screened run, and `assembly_stats.tsv` marks
  the genomes that passed (`screen` column).
- **Resume.** Each decision records both thresholds and the resolved paths
  of the GTDB sketch and lineages CSV. A repeat with the same settings
  screens nothing again: an excused run is excused from its recorded
  decision without a download, and a finished run keeps its decision. When a
  setting changed, the runs are screened again from their kept reads
  sketches without fetching their reads. An excused run is fetched again
  only when it now passes. A finished run that now fails is excused, while
  its contigs and marker stay for a later call. A finished run assembled
  without the screen is screened from its reads sketch. One without a reads
  sketch is not screened, and one INFO line counts such runs.
- **With `--max-runs`.** A run excused by the screen while it is processed
  does not give its place to a deferred run in the same call, since the
  places are fixed before any run is fetched. The next call excuses it from
  its recorded decision before the places are counted, so a deferred run
  takes its place then.

## Adding assemblies to a GTDB selection

`assemble --append` adds the assemblies to a working directory that already
holds a selection from `metadata` and `genome` (or `ingest`). The existing
rows, genome files and outgroup stay, a run assembled earlier is replaced,
and the manifest gains the new genomes with source `sra`. A mixed GTDB and
reads-derived set then dereplicates together.

`metadata` and `ingest` replace the selection, and the genome stage prunes
what the manifest no longer lists. Both therefore refuse to rerun while
appended genomes are present. `--drop-foreign` discards them deliberately,
and `assemble --append` can put them back afterwards. An `ingest` that takes
the same genomes again with `--from-workdir` does not drop them and is not
refused. Genomes with source `sra` that an earlier `ingest --from-workdir`
brought in are not protected this way, unless `assemble` has also run in
that working directory. To keep the two sets in separate working
directories instead, combine them with `ingest --from-workdir` (see
[Combining working directories](usage.md#combining-working-directories)).

## Assembling in parts

`assemble --max-runs N` keeps at most N runs with a finished assembly in the
working directory, so a large `reads.tsv` can be assembled in parts.

- Runs are taken in `reads.tsv` order, which the reads stage writes largest
  by bases first.
- Runs already assembled in the working directory count toward N first and
  are always kept, also when there are more of them than N. The remaining
  places go to the runs still to be assembled.
- Runs excused before any download (no FASTQ mirror, no accepting or
  installed assembler) take no place.
- The other runs are listed in `excused_runs.tsv` with reason `deferred`, so
  `dereplicate` runs on the assembled part. The selection and the manifest
  hold only assembled runs.

N counts runs in the working directory rather than runs in one call, so the
result depends only on `reads.tsv`, N and the finished assemblies. A repeat
with the same N assembles nothing new and is skipped as finished. A larger
N, or no `--max-runs`, assembles the next deferred runs and removes their
rows. A run counts toward N once assembled, also when the quality gate later
excuses it. A deferred run never counts as judged when every run is excused
(see [output.md](output.md#when-every-sequencing-run-is-excused)), and
`--append` takes `--max-runs` in the same way. `status` shows the stage as,
for example, `20 assembled, 35 deferred`.

`assemble --accession SRR... --accession-file FILE` names the runs to
assemble instead of a count. The other pending runs are deferred in the same
way, and finished runs are kept whether named or not. An accession not in
`reads.tsv` is refused, as is a file that lists no accessions (an empty
`passed.txt` means the screen excused every run). The `passed.txt` of
`reads-screen` is such a file, so a screened selection can be assembled
without editing `reads.tsv`.

The alternative is to bound the selection itself. `reads --max-runs N` keeps
the N largest runs in `reads.tsv`. A later `reads` call with a larger N
followed by `assemble` assembles only the new runs, since finished runs are
reused. In the Nextflow pipeline, where each run is one `assemble-run` task,
the bound is set on the selection with `--reads_args '--max-runs N'`, since
`assemble-run` takes no `--max-runs`.

## Stateless steps

The same work is available as three stateless steps. The Nextflow reads mode
runs them as separate tasks, and they also serve another scheduler:

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

- **`assemble-run`** writes `contigs.fasta` and the `assembly.ok` marker
  into its `--out` directory. It writes `excused_runs.tsv` instead when the
  run has no FASTQ mirror or an unsupported platform, or fails to download,
  assemble or polish. `--polisher`, `--polish-rounds`, `--keep-reads`,
  `--keep-files`, `--tool-arg` and `--reads-sketch/--no-reads-sketch` act as
  on `assemble`, and the reads sketch goes to `reads.sig.zip` in the same
  directory.
- **`assemble-run --screen-reads`** (with `--screen-min-fraction`,
  `--gtdb-sketch` and `--gtdb-lineages`) applies the taxon and fraction
  checks and writes `screen.json` beside the sketch. The duplicate check
  needs the other runs, so only `assemble` makes it.
- **`genome-qc`** reads a directory of such run directories
  (`--assemblies`) and writes `quality.tsv` and `classification.tsv` keyed
  by run accession. It needs at least one database.
- **`reads-gather`** applies the quality gate and the naming policy above
  from the optional `--qc` directory, collects the `screen.json` files into
  `screen_reads.tsv`, and writes an empty `outgroup_accession.txt`, since
  the reads chain has no outgroup at this point.

How the Nextflow reads mode passes options to these steps is described in
[nextflow.md](nextflow.md#reads-mode).
