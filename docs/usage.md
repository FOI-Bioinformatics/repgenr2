# Usage

This page is the command-line guide for RepGenR2 (the command is `repgenr`).
It follows the order of an analysis: building a genome set, dereplicating
it, inferring a phylogeny, and deriving a FlexTaxD taxonomy from the tree.

Three pages go into more depth. [reads.md](reads.md) covers genome sets
from sequencing runs, [phylogeny.md](phylogeny.md) the phylogeny routes, and
[nextflow.md](nextflow.md) the Nextflow pipeline. [install.md](install.md)
explains the installation, and [choosing-tools.md](choosing-tools.md) which
tools suit a dataset. Every option is listed in
[cli-reference.md](cli-reference.md), and every output file in
[output.md](output.md).

## Contents

1. [Introduction](#introduction): stages, working directory, threads,
   `run`, resume, `status` and `doctor`
2. [Building the genome set](#building-the-genome-set): census, GTDB,
   local genomes, sequencing reads, viruses, genome sketches
3. [Dereplication and representative selection](#dereplication-and-representative-selection)
4. [Phylogeny](#phylogeny)
5. [tree2tax](#tree2tax)
6. [Running tools in containers](#running-tools-in-containers)
7. [Troubleshooting](#troubleshooting), including the exit codes

## Introduction

### Stages and entry points

RepGenR2 builds a representative-genome repository for a taxonomic group and
the matching FlexTaxD taxonomy. A representative genome is the genome kept
for a cluster of closely related genomes. The pipeline runs five stages:

1. **metadata** selects genome accessions from GTDB.
2. **genome** downloads the selected genomes with NCBI Datasets.
3. **dereplicate** clusters the genomes by ANI (average nucleotide identity)
   and picks representatives.
4. **phylo** builds a phylogeny from the representatives.
5. **tree2tax** emits a FlexTaxD-compatible taxonomy from the tree.

The first two stages can be replaced by another entry path (see the
[README](../README.md#pipeline-routes)): `ingest` for genomes on disk,
`reads -> assemble` for sequencing runs, or `vmetadata -> vgenome` for
viruses. Each entry path writes the same genome set (`genomes/`,
`selection.tsv` and the manifest), and can stage an outgroup to root the
tree.

### The working directory

Every stage takes the working directory with `-wd` (`--workdir`). Stages
read and write inside it, and `repgenr.yaml` records what ran. Scratch and
downloads also live there and under `TMPDIR`, so put it on a disk with room.

### Threads

`-t/--threads` defaults to 16, or to the CPU limit of the process when that
is lower (a Slurm cpuset, `taskset`, or a cgroup quota such as `docker run
--cpus`). An explicit `-t` is used as given, with a warning above the limit.
Changing the thread count does not make a finished stage rerun.

### Running the whole chain

```bash
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` runs the bacterial chain, the viral one with `--viral`, or the reads
one with `--reads` ([reads.md](reads.md)), and
forwards the stage options it shares a name with; the rest are given on the
stage commands. `--with-snptype` adds the standalone `snptype` stage, so the
SNP tables under `snp/` are written even when the tree is built another way.
`--genomes-dir` starts from local genomes. `phylo --msa-source snptype` types
the genomes itself into `tree/msa/`, so it needs no `snptype` stage.

### Resume and `--force`

Each stage records its parameters, the digests of its inputs and the
container identity in `repgenr.yaml`, and skips when all three are
unchanged. A stage reruns when a parameter or `--container` changed, an
input changed (including an upstream stage's output), a main output is
missing, or its last run failed. A genome deleted from `genomes/` therefore
returns on the next run of the entry stage; to drop it, remove it from the
source directory or the `--selection` file and rerun that stage. A stage
that refuses before changing its outputs (exit 2, 3, or 4 at preflight)
keeps its last finished record, while a failed tool (exit 6) leaves the
record interrupted. Leftover files beside the genomes (`.fai`, `.tmp`,
README) do not trigger a rerun. `--force` reruns a stage explicitly. It is
a global option and goes before the command name (`repgenr --force
dereplicate -wd WD`).

Two changes are not detected and need `--force`: an in-place edit of a
genome that keeps its size and modification time, and an upgrade of a
natively installed tool. The full rules are in
[developing.md](developing.md#resume-semantics).

### Checking a working directory

`repgenr status -wd WD` shows what completed and what is next. A failed
stage is listed as `[interrupted]` and a stage whose input changed or whose
output is missing as `[stale]`; both rerun on their next invocation.
`repgenr doctor -wd WD` checks the outputs against the records (missing,
corrupt or untracked genomes, manifest drift, truncated deliverables) and
exits 7 on a failure and 0 on warnings alone. `doctor --quick` skips
reading the genome files, the slowest check on a large set. A record without
a resume fingerprint (from an earlier version, or restored by `derep-stock`)
is shown as `[done]` with a note that its next run recomputes it. Both
commands take `--json` (see [output.md](output.md#status-and-doctor)).
A malformed `repgenr.yaml` makes `status` and every stage exit 3.

## Building the genome set

### Counting genera, species and samples

`repgenr census` counts what a taxon or a working directory holds, without
writing anything. Use it to choose a target before downloading, and to check
a selection afterwards.

```bash
repgenr census -tf Francisellaceae --runs
```

```text
Francisellaceae: 9 genera, 48 species, 1205 genomes (48 representatives), 1681 runs, 1473 biosamples
name                    species  genomes  representatives  runs  biosamples  illumina  ont  pacbio  other
Francisella                  26     1157               26  1667        1462      1413   41      56    157
Caedibacter                   9       28                9     0           0         0    0       0      0
...
Allofrancisella (NCBI)        0        0                0    11           8         7    3       1      0
```

- `-tf/--target-family` or `-tg/--target-genus` counts the GTDB genomes and
  species representatives of a taxon, from the GTDB API (`--source api`,
  the default) or from a release table (`--source table -r 232.0`).
  `--runs` adds the ENA sequencing runs per platform. A row with runs only
  is marked `(NCBI)` and left out of the GTDB totals; it may be a GTDB taxon
  under another name. Runs for a suffixed GTDB name such as `Bacillus_A`
  are counted under the NCBI taxon without the suffix (`Bacillus`).
- `--viral --target TAXON` counts NCBI Virus records, complete records and
  isolates (the genomes `vgenome --group-segments` would form, an upper
  bound before the length filter). `--host`, `--complete-only` and
  `--released-after` filter it as in `vmetadata`. BV-BRC has no
  metadata-only report, so count a `vmetadata --source bvbrc` working
  directory instead.
- `-wd WD` counts a working directory: the candidates the entry stage chose
  from, then the selected genomes per source and the dereplication clusters.

`--tsv` and `--json` write the rows for scripts; the columns are described
in [output.md](output.md#census-table). A table census downloads the GTDB
table once into `$REPGENR_CACHE_DIR/gtdb` (default `~/.cache/repgenr/gtdb`).
The taxon flags are refused together with `-wd`, whose record names the
source and target.

### Bacteria from GTDB

`metadata` selects bacterial or archaeal genomes from GTDB, and `genome`
downloads them from NCBI.

```bash
WD=./francisella
repgenr metadata -wd $WD -r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella
repgenr genome -wd $WD
repgenr dereplicate -wd $WD --tool skder -t 16
repgenr phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
repgenr tree2tax -wd $WD --include-dereplicated
```

**The table or the API.** `--source tsv` (the default) reads the GTDB
metadata table of `-r/--release` and `--gtdb-version` (`bac120` or `ar53`),
downloaded into the working directory and checked against GTDB's
checksums. `--nodownload` reuses it and refuses a table fetched for another
minor release. `--metadata-path` reads a local copy, and `-r` and
`--gtdb-version` stay required then. A table placed by hand without its
`.release` file is reused with a warning. `--source api` queries only the target
taxon from GTDB's current release and ignores the table options with a
warning. It records the query date (`api_query_date`) in place of a release
number, and `status` and `versions` show it; compare it with the release
dates on the GTDB website to relate it to a release.

**Taxon names.** GTDB names are case-sensitive. The first letter of the
target is raised and a species epithet lowered, so `-tg francisella -ts
Tularensis` finds `s__Francisella tularensis`. Suffixes such as
`Bacillus_A` must be typed as GTDB spells them. An unknown taxon exits 2.

**The outgroup.** An outgroup is chosen from the parent taxon, outside the
target. `--outgroup-accession` names one instead, which must also lie
outside the target. A target genome that `-d rep` or `--limit` left out of
the selection is still not an outgroup candidate.

**Limiting the selection.** `--limit N` caps the selection at N genomes,
taken round-robin across species: the best genome of every species first,
then each species' next best. Within a species genomes rank by
completeness minus five times contamination, and ties go to the GTDB
representative flag, then the accession. With `-d rep` the rule is a quality
ranking across species. A heavily sequenced species
therefore cannot fill the cap on its own, and the outgroup does not count
against it. On the API path the quality card of every candidate is fetched
before the cut, four at a time. A card the API refuses twice leaves its
genome unscored, and unscored genomes come last.

### Starting from local genomes

`repgenr ingest` builds a working directory from local FASTA files, in
place of `metadata` and `genome`.

```bash
repgenr ingest -wd $WD --genomes-dir ./my_genomes --outgroup GCF_003574425.1
repgenr dereplicate -wd $WD --tool skder
```

`ingest` links the files (or copies them with `--copy`) into `genomes/` and
writes `selection.tsv` and the manifest. `repgenr run --genomes-dir DIR`
runs the same chain in one command, without the GTDB selection flags. Only files directly under the
directory with a suffix `.fasta`, `.fa`, `.fna`, `.fas`, `.fasta.gz`,
`.fna.gz` or `.fa.gz` are read. Suffixes are case-sensitive, so `X.FASTA`
and `x.fas.gz` are skipped with a warning. Taxonomy and quality come from `--selection
selection.tsv` when given, otherwise from the file name:

| File name | Accession | Taxonomy |
|---|---|---|
| Four or more `_`-separated tokens, as `Family_genus_species_ACCESSION.fasta` | the tokens after the third | the first three tokens |
| Starts with an NCBI assembly accession, as `GCF_000008985.1_ASM898v1_genomic.fna` | that accession | none |
| Any other name | the file stem | none |

Every name with four or more tokens is parsed this way, so
`sample_1_run_A.fasta` gives the accession `A`; give such files their
accessions with `--selection`.

Subdirectories are not searched, so collect the `.fna` files of an NCBI
Datasets download into one directory first. An empty or unreadable file
stops `ingest` before anything is staged. Two files with one accession stop
`ingest`; give them distinct accessions with `--selection`, whose
`is_outgroup` and optional `gtdb_representative` columns take `1` or `0`
(also `true`/`false`, `yes`/`no`). Without a `gtdb_representative` column no
genome is a GTDB representative, so `--keeper gtdb` acts as `quality`.
`--outgroup` names a genome under `--genomes-dir`, in a
`--from-workdir`, or a FASTA file anywhere, and sets it aside under
`outgroup/`.

Gzipped genomes are staged unchanged. Tools that cannot read gzip receive a
decompressed copy in the stage's scratch directory, which therefore needs
room for the uncompressed genomes; the stage stops when less than about four
times the gzip size, or 1 GB, is free. `--copy` makes the working directory
independent of the source, and the copies keep its permissions. A corrupt
gzip genome stops `snptype` or `phylo` with exit 2. A changed
source directory reruns
`ingest`, which removes the genomes no longer in the source and marks the
later stages for a rerun. Links dangle when the
source moves; `doctor` reports them, and `ingest` from the new location
restages the set.

#### Combining working directories

`--from-workdir WD` (repeatable) takes the genome set of another working
directory, for example GTDB genomes from one and assemblies from another:

```bash
repgenr ingest -wd WD3 --from-workdir WD1 --from-workdir WD2 --outgroup GCF_003574425.1
```

Accession, taxonomy, CheckM values and `gtdb_representative` carry over,
so `--keeper quality` and `--keeper gtdb` work on the combined set. The
manifest keeps each genome's source (`gtdb`, `sra`, `local`). Per-run files
of a reads working directory stay there. The outgroups of the sources do not
carry over; name one with `--outgroup` (filename, stem or accession, also
from a source's `outgroup/`). To add assemblies to a GTDB selection within
one working directory, use
[`assemble --append`](reads.md#adding-assemblies-to-a-gtdb-selection)
instead. A genome found in two sources stops `ingest`. The selection
and `genomes/` of each source are resume inputs, but its manifest is not, so
a change of a source label alone needs `--force`. A source whose stage record is interrupted, stale or absent
is named in a warning, and `--strict-sources` refuses it instead.

### Starting from sequencing reads

`repgenr reads` selects runs from ENA, and `repgenr assemble` fetches and
assembles them into `genomes/`:

```bash
repgenr reads -wd $WD -ts "Francisella tularensis" --platform illumina --max-runs 20
repgenr assemble -wd $WD --assembler auto -t 16 --jobs 2
```

[reads.md](reads.md) describes run selection, assembly, quality checks and
screening.

### Viruses

`vmetadata` selects virus genomes from NCBI Virus (or BV-BRC with `--source
bvbrc`), and `vgenome` stages them.

```bash
WD=./hav
repgenr vmetadata -wd $WD --target hepatovirus
repgenr vgenome   -wd $WD --target-genus Hepatovirus
repgenr dereplicate -wd $WD --tool skder
# or the whole chain:
repgenr run -wd $WD --viral --target hepatovirus -tg Hepatovirus --treebuilder mashtree
```

`--virus` passes virus-tuned settings to dRep (`--tool drep --virus`); the
other dereplicators do not read it and `dereplicate` warns.

**The outgroup.** `vgenome` picks a record of a sister species with enough
genomes, chosen by distance (mashtree by default). `--outgroup-accession`
names one instead (`run --viral --outgroup-accession` passes it on), and
`--no-outgroup` leaves the tree unrooted. With
`--group-segments` the outgroup is one record of the sister species, found
within the kept records' length span widened by 15 percent.

#### Species assignment

On the NCBI Virus path the species of a record is the current NCBI Taxonomy
species of its taxid, looked up once per taxid. Organism names are not used,
since they are often strain-level or earlier names. When the lookup gives
no species, the species is the shallowest binomial in the report lineage,
else the organism name; `species_source` in `virus_records.json` says which.
The species sets the filename token, the `--target-species` match and the
grouping of outgroup candidates. A `--target-species` value also matches
records by organism name (`-ts "Argentinian mammarenavirus"` selects
`Mammarenavirus juninense`), and `--target-serotype` still matches the
organism name.

**Segmented viruses.** `--group-segments` joins records that share a
species and an isolate name and carry at least two distinct segments into
one genome per isolate. Segment labels are normalised before comparison:
the text after a `;` and words such as `RNA`, `segment` or `genome` are
dropped, the first remaining word is upper-cased, and small, medium, middle
and large become S, M, M and L. `M`, `M; medium` and `middle` are therefore
one segment, while `RNA 1` and `RNA 2` stay distinct. Each isolate keeps one
record per segment (complete before partial, then the longest), and its
genome carries a synthetic `iso-` accession. `segments.tsv` lists the
members, and `tree2tax` lists them under the isolate's leaf.

**Length filter.** `vgenome` keeps records within a length window around
the median of the per-species median lengths, so an over-sequenced species
counts once. `--length-method mean` averages every record instead, which
lets the most-sequenced species set the window. `--length-range` sets the
window, and `--length-all` disables the filter. Switching a working
directory between NCBI Virus and BV-BRC removes the other source's tables.

### Genome sketches

The entry stages end by writing a sourmash sketch of each genome to
`sketches/<name>.sig.zip` (k=21, 31 and 51, scaled=1000). The sourmash
dereplicator, `glance`, the sourmash tree builder and the classifier of
`assemble` read them instead of sketching again. The step runs when
sourmash is available. `--sketch` requires it, and `--no-sketch` turns the
step off; a failed sketch is an error with `--sketch` and a warning
otherwise. Sketches are not stage outputs and do not affect resume, so a
stage skipped as complete writes none; run `repgenr sketch` then. A genome
edited in place with its size and modification time unchanged is not
re-sketched; `repgenr --force sketch` rewrites every sketch.
`repgenr sketch -wd WD` writes missing and stale sketches for an existing
working directory, and `status` then shows `sketches: n/m`. How sketches
are matched and cached is in
[developing.md](developing.md#genome-sketches).

## Dereplication and representative selection

`repgenr dereplicate` clusters the genome set by ANI and keeps one
representative per cluster, chosen with `--tool` (drep, skder, galah or
sourmash; see [choosing-tools.md](choosing-tools.md#dereplicators)):

```bash
repgenr dereplicate -wd $WD --tool skder -t 16
```

The clusters go to `derep/clusters.tsv` and the representatives to
`derep/representatives/`.

### Representative selection

After clustering, a keeper step re-picks each cluster's representative.
`--keeper quality` (the default) takes the best-quality genome, `gtdb`
the GTDB species representative, and `tool` the dereplicator's own pick.
Quality values come from the manifest: CheckM values of GTDB selections,
CheckM2 scores of assemblies, or the columns of `ingest --selection`. When no
genome has values the stage warns and records `keeper_effective: tool`.

#### The keeper score

The score is `completeness - 5 x contamination + 0.5 x log10(N50)`, dRep's
default weights. Ties go to higher completeness, lower contamination,
higher N50, then file name. A genome without quality values never becomes
the representative. When the dereplicator's pick has no values, a scored
genome replaces it only if it meets the MIMAG high-quality thresholds
(completeness above 90, contamination below 5), so a scored fragment does
not displace an unscored complete genome.

The working-directory stage applies the keeper once, to the final clusters.
The Nextflow steps apply it within each chunk and again after the merge, so
the final clusters can differ slightly.

With `--keeper gtdb`, a cluster keeps its GTDB species representative and
otherwise falls back to the quality rule. When a cluster holds several,
the best is kept and the others are named in a warning; raise
`--secondary-ani` to keep them apart. Only GTDB selections carry the flag,
so on the other entry paths `gtdb` acts as `quality`.

When every genome has quality values, they are also passed to dRep
(`--genomeInfo`, which then skips CheckM) and galah (`--genome-info`).
Otherwise galah receives the genomes by descending file size, so a complete
genome rather than a fragment becomes its representative.

#### One representative per taxon

`--reduce species|genus` keeps one representative per taxon after
dereplication: the representative of the taxon's largest cluster, re-picked
by score only with `--keeper quality` and known scores. A GTDB
representative is preferred under `--keeper gtdb`. `--target-reps N` searches the secondary ANI to land near N
representatives.

#### Finer resolution inside a lineage

`--refine SEED --refine-ani F` dereplicates the clusters that hold the
seeds again at a finer ANI, leaving the rest of the result unchanged. Use it
for coarse coverage of a species with more genomes inside one lineage.

```bash
repgenr dereplicate -wd $WD --tool skder --refine GCA_000123.1 --refine-ani 0.999
```

`derep-refine` does the same on a finished result without repeating the
coarse pass, and can be applied again for another lineage. It needs a
completed `dereplicate` record (one restored by `derep-stock` counts) and
records a `derep_refine` stage, after which `phylo` and `tree2tax` rebuild
from the refined set. `run --refine --refine-ani` passes both to
`dereplicate`. Refining a seed again starts from its original coarse
cluster, and on a `--reduce` result a lineage is refined as the whole taxon
group. A later
`dereplicate` discards the refinement, so store the result with
`derep-stock` first. A seed given as an accession needs its version suffix.
`--reduce` leaves refined clusters apart, and `--refine` cannot be combined
with `--target-reps`. Several
lineages, each at its own ANI, are given with `--refine-file`, a TSV with
the columns `seed` and `ani` (lines starting with `#` are skipped). It
combines with `--refine`, whose `--refine-ani` then applies to its own seeds
only, and editing the file reruns the refinement:

```text
seed	ani
GCA_000123.1	0.999
GCA_000789.1	0.9995
```

The lineage genomes are listed in `derep/refined.tsv`. Further rules are in
[developing.md](developing.md#resume-semantics).

#### The sourmash dereplicator and its sketch cache

`--tool sourmash` reads the genome sketches. For other `ksize` or `scaled`
values it keeps its own sketch cache, and `--tool-arg sketch_cache=DIR`
makes that cache persist across runs. Empty that directory after replacing a
genome by a file with the same name, size and timestamp, as `cp -p` or
`rsync -t` produce. Both sourmash back-ends
estimate ANI from the average containment (see
[choosing-tools.md](choosing-tools.md#how-the-tools-estimate-ani)).

### Chunked dereplication

`dereplicate -s N` (`--process-size`) splits the genomes into chunks of at
most N, dereplicates each chunk with the `--pre-*` thresholds, then
dereplicates the union of the chunk representatives with the final ones.
Use it for very large sets or sets with many near-identical genomes. The
chunk size does not bound the final pass: for a diverse set the union is
most of the input. Three settings make the final pass smaller:

- keep skDER's default skani screen. For a cutoff of 99 or more the adapter
  skips the alignment of pairs whose k-mer estimate is more than one point
  below it. `--tool-arg screen=` changes it, `screen=skder` restores skDER's
  own, and `low_mem_greedy` ignores it;
- raise `-s`, so that each chunk holds more close relatives;
- set `--pre-secondary-ani` below `--secondary-ani`, which can give coarser
  clusters than a single pass at the final threshold.

The final pass prints no progress while it runs, and a warning appears when
it exceeds four times the chunk size.

An example run is in
[choosing-tools.md](choosing-tools.md#by-number-of-genomes), and the probe
rule for further levels is in [developing.md](developing.md#scaling).

### Inspecting a dereplication

| Question | Where to look |
|----------|---------------|
| How many clusters, and how large? | `derep/cluster_summary.tsv` (`repgenr cluster-summary` rebuilds it). |
| Which genomes are in the cluster of X? | `derep/clusters.tsv`, or `derep/unpacked/<X>/` after `repgenr derep-unpack`. |
| What happened to one genome? | `derep/genome_status.tsv`. |
| Do the genomes fall into clear groups? | `repgenr glance`. |
| How do two dereplications differ? | `repgenr derep-stock`, see below. |

These commands do not check that the last `dereplicate` finished, so run
`repgenr status` first.

#### Finding the members of a cluster

```bash
REP=$(sed -n 4p $WD/derep/cluster_summary.tsv | cut -f1)   # third-largest cluster
awk -F'\t' -v r="$REP" '$1 == r {print $2}' $WD/derep/clusters.tsv
```

`repgenr derep-unpack` writes the same grouping as one directory per
cluster (`--no-representant` leaves the representative out). The files are
hard links where the file system allows them, so editing one in place also
changes the genome. A member missing from `genomes/` is skipped with a
warning. Rerun it after a new dereplication; `doctor` warns until then.

#### Storing and restoring a dereplication

`repgenr derep-stock --action pack --name RUN` stores the current result
under `derep/stock/RUN`, `--action list` and `--action delete` manage the
stored runs, and `--action unpack` makes one current again
(details in [output.md](output.md#storing-and-restoring-with-derep-stock)).
A restored result has no resume fingerprint, so the next `dereplicate`
recomputes it; run `phylo` directly to build on it.

#### Reading `glance`

`repgenr glance` compares all genomes in `genomes/` and plots the pairwise
ANI values with a dendrogram. `--tool auto` (the default) uses dRep (Mash
ANI) when it can run and sourmash otherwise; sourmash refuses sets above 5000
genomes. Compare plots made with the same tool. dRep fails on exFAT volumes,
so keep the working directory on APFS or ext4. With `--tool sourmash
--keep-files`, `glance_wd/pairwise_ani.csv` holds every pair; an empty cell
means the pair is below the sketch's reach, not missing data. When the
histogram shows separate groups with an empty gap between them, a threshold
in the gap separates them. When the values are spread evenly, the threshold
decides the cluster sizes, so check `cluster_summary.tsv` afterwards.

#### Comparing two dereplications

```bash
repgenr dereplicate -wd $WD --tool sourmash && repgenr derep-stock -wd $WD --action pack --name sourmash
repgenr dereplicate -wd $WD --tool galah && repgenr derep-stock -wd $WD --action pack --name galah
part() { awk -F'\t' 'NR==FNR { if (NR>1 && (!($1 in m) || $2<m[$1])) m[$1]=$2; next }
                     FNR>1 { print $2 "\t" m[$1] }' "$1" "$1" | sort; }
diff <(part $WD/derep/stock/sourmash/clusters.tsv) \
     <(part $WD/derep/stock/galah/clusters.tsv) && echo "same partition"
```

Two dereplicators can form the same clusters and keep different
representatives, so the script labels each genome with the smallest member
name of its cluster and compares the labels.

## Phylogeny

`repgenr phylo` builds `tree/tree.nwk` from the representatives and the
outgroup:

| Route | Option | Suits |
|-------|--------|-------|
| Whole-genome alignment | `--aligner progressivemauve` | under about 30 genomes |
| Core-SNP alignment | `--msa-source snptype` | one species or a clonal set |
| Core-gene alignment | `--msa-source coregene`, after `repgenr annotate` | a genus, or wider with `--coregene markers` |
| Alignment-free | `--treebuilder mashtree` or `skani` | any set, as a quick look |

```bash
repgenr phylo-advise -wd $WD                     # prints a route and a --preset
repgenr phylo -wd $WD --preset species-masked    # the advised route
```

[phylogeny.md](phylogeny.md) describes the routes, presets, typers,
builders, support values and rooting, and
[choosing-tools.md](choosing-tools.md#phylogeny-routes) compares them.

## tree2tax

`repgenr tree2tax` turns the rooted tree into a FlexTaxD taxonomy:
`tree2tax.tsv` (the relations), `genomes_map.tsv` and `tree2tax_nodes.tsv`
(see [output.md](output.md#tree2tax)). By default the genomes that
dereplication left out are listed under their representative.

```bash
repgenr tree2tax -wd $WD --include-dereplicated
```

### Rooting in tree2tax

`tree2tax` roots the tree on the outgroup when it is not rooted there yet.
It exits 3 when the outgroup is not a leaf of the tree. After `phylo --root
midpoint` it keeps the midpoint root, and after `phylo
--no-outgroup --root none` it leaves the tree unrooted. The
`tree2tax-relations` step takes `--keep-root` and `--no-outgroup` for the
same choices; its `--no-outgroup` ignores `--outgroup-dir`. `--remove-outgroup` prunes the outgroup leaf before the collapse and naming
steps.

### Collapsing weak splits in tree2tax

Three options merge internal nodes into their parents before naming, so a
split the data do not support does not become a FlexTaxD node. All are off
by default.

- `--collapse-length L` merges a node whose branch is shorter than L.
- `--collapse-support S` merges a node whose support is below the fraction
  S. Percentages from IQ-TREE and RAxML-NG are normalised. With both `-B`
  and `alrt=`, the ultrafast bootstrap (the last field of the label) is
  read. FastTree and aster write fractions. Builders without
  supports (mashtree, sourmash, skani, ska2) make the option a no-op with a
  warning.
- `--collapse-diameter D` keeps the fewest clades whose leaf-to-leaf
  distance stays within D and merges every node inside them (TreeCluster's
  `max_clade`). A group of near-identical genomes then becomes one node.

Lengths are in the tree's units:

- mash distance for mashtree (about 0.05 at 95 % ANI);
- one minus ANI for the sourmash builder and skani (0.05 at 95 % ANI). For
  sourmash, pairs below the sketch's reach sit at 1.0, so D never joins them;
  the reach ends near 88 to 90 % ANI at k=31 and 83 to 85 % at k=21. With
  sourmash `metric=jaccard` it is one minus the Jaccard similarity (about
  0.88 at 95 % ANI with k=31);
- SNPs (`metric=snps`) or mismatching k-mers (`metric=mismatch`) for ska2;
- substitutions per site for the ML builders.

The diameter rule runs first. `--collapse-length` can still merge a clade
it kept, which puts leaves more than D apart under one node, and an ingroup
within D becomes one node under the root. The root, the outgroup split and
the leaves never collapse, and a collapsed node's length is added to its
children's.
`tree2tax_nodes.tsv` lists every node with its length, support and depth,
and is the place to look when choosing a threshold.

### Node names

Internal nodes are named by a hash of their leaf set, so a node is renamed
whenever a genome joins or leaves it. `--naming lca` names a node by the
smallest leaf label under each child, which changes less often. A genome
that sorts after the others and splits a branch then leaves every name
unchanged. A genome that joins a multifurcation, that sorts below a child's
smallest label, or that is removed still renames the nodes it touches, and
so does a topology change in the re-inferred tree. `--node-basename`
numbers the nodes instead.

### A rank column

`--ranks` adds a third column with `strain` on every genome and `no rank`
on every internal node. `--rank-red species=0.95,genus=0.8` names internal
nodes by relative evolutionary divergence (RED, the `red` column of
`tree2tax_nodes.tsv`); each rank appears at most once on a root-to-leaf
path, and a lineage may lack a rank no node reaches. The thresholds are the
user's. GTDB's median RED values (about 0.92
for genus and 0.77 for family; Parks et al. 2018) transfer to a genus-level
tree only loosely, so inspect `tree2tax_nodes.tsv` first. `--rank-red`
needs branch lengths and a rooted tree: an outgroup, `phylo --root
midpoint`, or `--keep-root` on `tree2tax-relations`.

### FlexTaxD version and depth

The `tree2tax` files follow the reader on FlexTaxD's GitHub master
(`flextaxd/modules/ReadTaxonomy.py`, 2026-10), which accepts either column
order, reads the header to tell, pre-creates `root` and validates from it.
The PyPI 0.1.5 release reads the genome map with a header line, stores the
edges of a child-first file inverted and does not pre-create `root`, so it
cannot consume these files. FlexTaxD walks at most 50 levels below `root`
(`get_children`, `maxdepth=50`): a deeper taxonomy imports and dumps in full
but fails `flextaxd --validate`. An unbalanced tree of more than about 50
genomes can reach that depth; `--collapse-diameter` reduces it. A rank
column is read on every row, and an empty value inherits the previous
row's, so `--ranks` writes a rank on every row.

## Running tools in containers

RepGenR2 can run each external tool in a pinned container, while RepGenR2
itself runs on the host:

```bash
repgenr --container docker dereplicate -wd $WD --tool drep
repgenr --container singularity --container-cache /data/repgenr_sif phylo -wd $WD
```

`--container` is a global option placed before the command; the default is
`none`. The other options, the image sources and the Apple Silicon setup
are in [install.md](install.md#3-containers), and the Nextflow profiles in
[nextflow.md](nextflow.md#profiles-and-containers).

### Per-tool directories (`--bin-dir`)

`--bin-dir TOOL=DIR` puts DIR first on `PATH` for one tool, for a tool kept
in its own conda environment:

```bash
P=$(conda info --base)/envs
repgenr --bin-dir gubbins=$P/repgenr-gubbins/bin phylo -wd $WD --msa-source snptype --mask gubbins
```

TOOL is an adapter name from `repgenr list-tools`, or `checkm2`, `datasets`
or `minimap2`. `REPGENR_BIN_DIRS="gubbins=/p/bin,mashtree=/q/bin"` sets the
same, and a command-line `--bin-dir` replaces the variable's entry for that
tool. Under `--container`, a tool that runs in an image ignores it. The
directory is recorded in `repgenr.yaml` but does not trigger a
rerun.

## Troubleshooting

**A tool is not found** (`MissingBinaryError`, exit 4). The package does
not install the tools ([install.md](install.md)). Put the tool on `PATH`,
or use `--container`. `repgenr list-tools --check` reports which tools are
missing or too old.

**Apple Silicon.** BioContainers are amd64. Pass `--platform linux/amd64`
and enable Rosetta in Docker Desktop.

**Where are the details of a failure?** The console shows one line; the
traceback, the tool's command line and its output are in
`<workdir>/repgenr.log`, as are the container engine options (mounts,
user) under `--container`. `--verbose` prints them on the console. A
stateless step has no log and prints the tool's output tail itself.

**GTDB or NCBI is unreachable.** A host that does not accept a connection
within 15 s counts as unreachable, and the stage exits 3 and names it.
`genome` and `vmetadata` check the NCBI Datasets host before calling
`datasets`. The requests honour `HTTPS_PROXY`, `HTTP_PROXY` and `NO_PROXY`,
and `REQUESTS_CA_BUNDLE` points Python at another certificate bundle.
`REPGENR_SKIP_NET_PROBE=1` skips the check. BV-BRC uses FTPS, which does not
go through these proxies, so on a proxy-only network use the default NCBI
Virus source.

**Genomes NCBI no longer serves** are recorded in `missing_accessions.txt`
and excused by later stages, and a rerun asks for them again. An outgroup
NCBI no longer serves exits 3; choose another with `metadata
--outgroup-accession`.

**NCBI Entrez throttling** (BV-BRC path): set `NCBI_API_KEY`, and
optionally `NCBI_EMAIL`.

**ParSNP on three genomes** fails (exit 6), since ParSNP needs four.

**A tool hangs.** `REPGENR_SUBPROCESS_TIMEOUT=<seconds>` caps every tool;
the tool's process group is killed on expiry.

**Stopping a run.** On SIGTERM, SIGHUP or Ctrl-C, repgenr stops every
running tool and exits (143 for SIGTERM, 129 for SIGHUP), and the stage is
listed as `[interrupted]`. Under `--container docker` it also stops the
tool's container. Ctrl-Z suspends the tools with repgenr; without a
terminal SIGTSTP has no effect. A second signal kills the tools at once. `nohup` keeps a
run alive when the terminal closes. SIGKILL cannot be handled: stop a tool
with `pkill -g <tool pid>`, and leftover containers by label:

```bash
docker ps -q --filter label=repgenr.pid=<repgenr pid> | xargs -r docker stop
```

### Working directories from earlier versions

- Each stage reruns once, since the resume fingerprint changed. The
  manifest is upgraded in place, after which earlier versions refuse it.
- A selection made before the GTDB representative flag existed needs
  `repgenr metadata --force` with the same options for `--keeper gtdb`.
- The keeper rule is not part of the resume fingerprint, so a working
  directory dereplicated under an earlier rule keeps its representatives
  until `repgenr --force dereplicate`.
- An NCBI Virus working directory from before the taxonomy lookup needs one
  `repgenr --force vmetadata` with the same arguments, which renames the
  genomes.

### Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success. |
| 1 | An unexpected error (traceback in the run log). |
| 2 | Invalid or missing user input. |
| 3 | The working directory is missing or in a bad state, or a remote service (GTDB, NCBI, BV-BRC, ENA) failed. `assemble` and `reads-gather` also exit 3 when every run was excused (see [reads.md](reads.md#when-every-sequencing-run-is-excused)). |
| 4 | A required tool is absent or too old, or the Docker daemon is unreachable. |
| 5 | A tool adapter could not be found or loaded. |
| 6 | An external tool failed, including `datasets`. Under `REPGENR_PROPAGATE_TOOL_EXIT=1` (set by the Nextflow modules) the tool's own status is forwarded. |
| 7 | `doctor` found at least one failure. |
