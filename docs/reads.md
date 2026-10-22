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
5. [Screening runs](#screening-runs)
6. [Adding assemblies to a GTDB selection](#adding-assemblies-to-a-gtdb-selection)
   and [assembling in parts](#assembling-in-parts)
7. [Stateless steps](#stateless-steps)

## Overview

Use this route for isolates that are public as reads but not as assemblies.
`repgenr reads` selects runs from ENA (which mirrors SRA) and writes
`reads.tsv`. `repgenr assemble` fetches and assembles them into `genomes/`
with the same `selection.tsv` and manifest as the other entry paths.

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina --max-runs 20
repgenr assemble -wd $WD --assembler auto -t 16 --jobs 2
```

`repgenr run --reads` runs the whole chain and forwards `--accession-file`,
`--platform`, `--max-runs`, `--assembler`, `--threads` and `--outgroup`.
Which assembler and polisher suit each platform is in
[choosing-tools.md](choosing-tools.md#assemblers-and-polishers), and the
files are described in [output.md](output.md#assemble).

## Selecting runs

Runs are chosen by taxon (`--target-family`/`-tf`, `--target-genus`/`-tg`
or `--target-species`/`-ts`, resolved through the ENA taxonomy) or by accession (`--accession` with a run, sample or study,
repeatable, or `--accession-file`). Only whole-genome sequencing of genomic
DNA is kept, also for runs named by accession. A warning follows above
100000 runs under a taxon, since their records need gigabytes of memory.

The selection is then filtered:

- `--platform illumina|ont|pacbio` keeps one platform.
- `--min-bases` and `--max-bases` drop runs outside a size range. An
  unenriched host library for a small endosymbiont assembles into a
  host-dominated genome.
- `--drop-selection` drops runs by ENA library selection; the default,
  `MDA`, drops whole-genome amplification, which gives chimeric contigs.
  The flag repeats, and `--drop-selection none` keeps every run.
- `--one-per-sample` (the default) keeps the best run of each sample: a
  long-read run with at least 100 Mb and a tenth of the largest short-read
  run, else the largest run. `--all-runs` keeps every run.
- `--max-runs` keeps the largest runs.

Each run is labelled with the family, genus and species of its NCBI taxid.

## Assembly

`assemble` downloads each run's FASTQ files from ENA over HTTPS, verifies
them against ENA's checksums, assembles them, and keeps the contigs of at
least `--min-contig-length` bases (500). The genome is named
`Family_genus_species_RUN.fasta`.

- **Assembler.** `--assembler auto` (the default) picks `skesa` for
  Illumina and `flye` for ONT and PacBio; `shovill` is the alternative for
  paired Illumina runs. `--tool-arg` passes tuning such as
  `mode=nano-raw` to Flye.
- **Layouts.** Long-read runs are assembled whatever layout ENA gives them,
  and a run listed as several FASTQ files is joined. A paired run with a
  third, orphan file gives skesa all three and shovill the pair. A
  short-read run listed with a single FASTQ file is assembled as single-end
  by SKESA, and shovill excuses it as `unsupported_layout`.
- **Parallel runs.** `--jobs` runs that many assemblies at once (2, or 1
  when a long-read run is pending), with `--threads` split across them.
  `--memory-gb` (16) caps SKESA and shovill.
- **Resume.** A finished run leaves a marker under `assemblies/<run>/` with
  the settings it was built with, and is reused while they agree. A higher
  `--min-contig-length` filters the finished contigs again; any other change
  assembles the run again. To force a rebuild, delete `assemblies/<run>/`.
  Reads are deleted after a successful assembly unless `--keep-reads` is
  given, and the assembler's scratch unless `--keep-files` is given. A run
  interrupted during assembly keeps the FASTQ files that match their
  checksums, so the resume does not download them again.
- **Excused runs.** A run is written to `excused_runs.tsv`, and the rest
  proceed, when it has no ENA FASTQ mirror, a download that fails its
  checksum, a platform no assembler accepts
  (`unsupported_platform`), an assembler not installed under `auto`
  (`assembler_not_installed`), or a failed assembly. With no assembler installed at
  all the stage exits 4. A failed download is retried only with `--force`.
  When every run is excused, see
  [output.md](output.md#when-every-sequencing-run-is-excused).
- **Outgroup and statistics.** `--outgroup FASTA` sets a genome aside for
  rooting. `assembly_stats.tsv` gives contigs, length, N50 and coverage per
  assembly.

## Polishing long-read assemblies

Long-read assemblies are polished with the run's own reads. `--polisher
auto` (the default) runs medaka for ONT runs (racon when medaka is missing)
and racon for PacBio CLR runs, with `--polish-rounds` rounds, and nothing
for HiFi or Illumina. `--polisher none` turns polishing off.

medaka needs the basecaller model, which it reads from the FASTQ headers of
reads basecalled with Dorado. SRA rewrites the FASTQ headers that name
it, so give it with `--tool-arg model=...`. Otherwise ONT's bacterial
R10.4.1 model is used (`--tool-arg bacteria=false` uses medaka's general
default). A polishing failure excuses the run with `polish_failed`. When the
polisher is not installed the run is assembled unpolished with a warning.
Unpolished ONT assemblies carry indels that CheckM2 reads as lower
completeness.

## Quality scoring and classification

**CheckM2.** With a CheckM2 database (`--checkm2-db` or `CHECKM2DB`), every
assembly is scored, and the values reach `selection.tsv` for `--keeper
quality`. An assembly below `--min-completeness` (50) or above
`--max-contamination` (10) is excused with `qc_failed`. An assembly CheckM2
gives no result for is kept with a warning and without values. Obtain the
database with `checkm2 database --download`. The scores are
stored per run, so changing only the thresholds does not rerun CheckM2; to
score a run again, delete `assemblies/<run>/checkm2.json`.

**GTDB classification.** With a GTDB sourmash sketch and its lineages CSV
(`--gtdb-sketch` and `--gtdb-lineages`, or `REPGENR_GTDB_SKETCH` and
`REPGENR_GTDB_LINEAGES`), each assembly is classified by `sourmash gather`.
The `gtdb-rs226-reps.k31-sc10k.sig.zip` sketch from
`https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs226/` serves.
`--classifier auto` (the default) runs it when a sketch is configured, and
`none` turns it off. When the GTDB genus agrees with the
submitted organism, the GTDB names are used for the file name. Otherwise the
submitted name stays and `assembly_stats.tsv` flags the genome as
`genus_renamed` (same family and epithet, as NCBI Mycoplasmopsis arginini
and GTDB Metamycoplasma arginini) or `classifier_disagrees`. Neither
excuses the genome, and both lineages are kept in the table.

Each gather holds about 0.6 GB of memory, so as many run at once as
`--threads` and `--memory-gb` / 0.6 allow. Database paths and binaries are
checked before any run is fetched, so a wrong path exits 2 and a missing
tool exits 4 at once.

## Screening runs

A screen judges a run's reads against GTDB before the run is assembled, to
leave out runs of another taxon, runs that match their species poorly, and
repeated sequencing of one isolate. It reads a sourmash sketch of the reads,
made either before any download (`reads-sketch`, `reads-screen`) or during
assembly (`assemble --screen-reads`).

### The three checks

- `taxon_mismatch`: the GTDB genus of the reads (`sourmash gather` and
  `sourmash tax genome`) differs from the submitted genus, with the
  classifier's `genus_renamed` tolerance. A run with no GTDB match fails too.
- `low_match_fraction`: the abundance-weighted fraction of the reads in the
  top species is below the minimum fraction (default 0.3).
- `duplicate_isolate`: the run's reads are contained in those of a larger
  run of the same biosample at an ANI estimate of 0.999 or more (default).

The match fraction mostly measures the distance between the run's strain and
the GTDB representative of its species, and host or error k-mers lower it as
well. Four Illumina isolates gave 0.39 to 0.58 against the GTDB rs226
representatives, so check the screen table before relying on the default.
How the duplicate estimate is computed is in
[choosing-tools.md](choosing-tools.md#how-the-tools-estimate-ani).

### Before download: reads-sketch and reads-screen

```bash
repgenr reads-sketch -wd $WD --max-spots 200000
repgenr reads-screen -wd $WD --max-spots 200000 --gtdb-sketch $DB --gtdb-lineages $LIN
repgenr assemble -wd $WD --accession-file $WD/reads_screen/spots_200000/passed.txt
```

`reads-sketch` streams each run of `reads.tsv` from NCBI with `fastq-dump
-Z` into sourmash, without writing FASTQ files. The sketches have the
parameters of the reads sketch of `assemble` (k=21, 31 and 51, with
abundances).

- `--max-spots N` stops after the first N spots, and each cap is its own
  variant directory (`reads_sketches/all/`, `reads_sketches/spots_<N>/`).
- `--accession` and `--accession-file` restrict the call to runs of
  `reads.tsv`, and a later call adds the rest.
- The index `reads_sketches.tsv` gives each run's status (`ok`, `failed`,
  `empty` or `pending`; see [output.md](output.md#reads-sketch)). A failed
  run is streamed again on the next call, and a run that keeps failing is
  retried on every call until the selection is narrowed with `--accession`.
  The command exits 3 only when every selected run failed. `--force`
  streams every selected run again, `empty` ones included.
- `--jobs` (4) streams runs in parallel, each costing about one CPU.

The first spots of a run are not a random sample, and a low-abundance
organism may be missed in a short prefix; a few hundred thousand spots is a
suggested starting point. A capped sketch also makes the duplicate check
compare the starts of two runs. `fastq-dump` needs sra-tools 3.0 or later
and reads from NCBI, so a run that only ENA holds fails. Under `--container
docker` it runs from NCBI's `ncbi/sra-tools` image, which also runs natively
on arm64. A run whose stream holds no biological read is marked `empty` and
not streamed again at that cap. A stalled stream is limited by
`REPGENR_SUBPROCESS_TIMEOUT`, and a refusal naming a `.lock` file means
another call holds the variant directory. Retries and locking are described
in [developing.md](developing.md#resume-semantics).

`reads-screen` judges one variant (`--max-spots N`, required when more than
one variant exists) with the three checks
(`--min-fraction`, `--dup-ani`) and writes `reads_screen/<variant>/`.
`passed.txt` there lists the accepted runs in the form `assemble
--accession-file` takes (see [output.md](output.md#reads-screen)). A run
without a sketch in the variant is listed as `screen_failed: no reads
sketch`; run `reads-sketch` for it first. Repeats reuse earlier decisions
made with the same settings, and `--force` screens every selected run
again. A run decided with
other settings and not selected in the call is left out of `passed.txt`
with a warning, so screen every run with the same settings before using the
file.

Without a working directory, `reads-sketch -o DIR` sketches bare run
accessions or the runs of `--reads-tsv`, and `reads-screen --reads-tsv
reads.tsv --sketches-dir DIR -o OUT` screens them. Both record nothing and
take `--versions-out`. `reads-sketch --allow-all-failed` makes a call where
every run failed a warning, for one chunk of a larger batch. The Nextflow
tasks use these forms.

### During assembly: assemble --screen-reads

`assemble` sketches each run's reads beside the assembler into
`assemblies/<run>/reads.sig.zip`, when sourmash is available
(`--reads-sketch/--no-reads-sketch`). The reads sketch is not a genome
sketch and is never copied into `sketches/`. A failed reads sketch never
fails the assembly. A run finished without a reads sketch gets one only if
its reads were kept with `--keep-reads`, and the call needs `--force` when
the stage is otherwise up to date.

`--screen-reads` (off by default) screens that sketch before the assembler
starts, with `--screen-min-fraction` (0.3) and `--screen-dup-ani` (0.999),
and excuses a failing run in `excused_runs.tsv` with step `screen`. The
reads are then fetched before the check, but the assembly cost is saved.
Duplicates are compared among the runs of one call and the finished runs
with a reads sketch. It
needs sourmash and the GTDB sketch and lineages CSV (`--gtdb-sketch`,
`--gtdb-lineages`), and refuses to start without them. It runs with
`--classifier none` as well. A tool failure excuses the run with
`screen_failed`, which does not count as a judgement of the run's data. An
excused run keeps its reads sketch and decision under `assemblies/<run>/`
for inspection. `screen_reads.tsv` lists every screened run.

`--screen-reads` does not read the decisions of `reads-screen`; it screens
again from its own sketch of the whole run. To assemble only the runs
`reads-screen` passed, give its `passed.txt` to `--accession-file`. On a
repeat with changed thresholds, the runs are screened again from their kept
reads sketches without fetching reads. An excused run is fetched again only
when it now passes, and a finished run that now fails is excused while its
contigs stay. The full rules are in
[developing.md](developing.md#resume-semantics).

## Adding assemblies to a GTDB selection

`assemble --append` adds the assemblies to a working directory that already
holds a GTDB selection, so a mixed set dereplicates together. The existing
rows and outgroup stay, and the new genomes get the manifest source `sra`. `metadata` and
`ingest` then refuse to rerun while appended genomes are present;
`--drop-foreign` discards them, and `assemble --append` can put them back
afterwards. Genomes with source `sra` brought in by an earlier `ingest
--from-workdir` are protected only if `assemble` has also run in that
working directory. With `--append`, a run assembled earlier is replaced, and
`--max-runs` applies in the same way. To keep the two sets apart, combine them
with [`ingest --from-workdir`](usage.md#combining-working-directories)
instead.

## Assembling in parts

`assemble --max-runs N` keeps at most N assembled runs in the working
directory, taken in `reads.tsv` order (largest first). Runs already
assembled count first, and the others are listed as `deferred` in
`excused_runs.tsv`. A larger N, or no `--max-runs`, assembles the next
deferred runs. A run counts toward N once assembled, also when the quality
gate later excuses it, and `status` shows, for example, `20 assembled, 35
deferred`. `assemble --accession-file FILE` names the runs to assemble
instead, for example the `passed.txt` of `reads-screen`; an empty
`passed.txt` means the screen excused every run. `reads --max-runs N` bounds
the selection itself, and a later `reads` with a larger N followed by
`assemble` builds only the new runs. That is the form the Nextflow pipeline
uses.

## Stateless steps

The same work is available as three stateless steps, which the Nextflow
reads mode runs as separate tasks:

```bash
repgenr assemble-run --reads-tsv $WD/reads.tsv --run SRR25474756 -o asm/SRR25474756 --assembler auto -t 8
repgenr genome-qc --assemblies asm -o qc --checkm2-db $CHECKM2DB --gtdb-sketch $DB --gtdb-lineages $LIN
repgenr reads-gather --reads-tsv $WD/reads.tsv --assemblies asm --qc qc -o out
```

None of the steps touches a working directory or the manifest.
`assemble-run` assembles one run and takes the per-run options of
`assemble`, including `--polisher`, `--polish-rounds` and
`--reads-sketch`; its `--screen-reads` applies the taxon and fraction checks
only. `genome-qc` scores and classifies every run directory
(`--assemblies`) and needs at least one database (CheckM2 or the GTDB
sketch); `--memory-gb` and `--classifier` act as on
`assemble`. `reads-gather` applies the quality gate and the naming, and
writes the genome set from `--assemblies` and the optional `--qc`
directory, with an empty `outgroup_accession.txt`.
