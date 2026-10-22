# Choosing tools

RepGenR2 has ten tool families in four groups: read tools (assemblers,
polishers, classifiers), genome tools (dereplicators, annotators), alignment
tools (aligners, SNP typers, maskers, core-gene builders) and tree tools (tree
builders). Most stages work with their defaults. This page is for users who
need to decide what to use when the defaults do not fit. It starts from four
properties of the dataset and one decision table, then describes the
trade-offs of each family, and ends with the genome limit each adapter
declares. A tool you choose must be installed; the environment that provides
each tool is listed in
[install.md](install.md#1-core-and-satellite-environments).

Contents:

- [Measured and declared figures](#measured-and-declared-figures)
- [Start from the question](#start-from-the-question)
- [By purpose and divergence](#by-purpose-and-divergence) (the decision table)
- [By number of genomes](#by-number-of-genomes)
- [Entry points](#entry-points) (where the genomes come from)
- Read tools: [Assemblers and polishers](#assemblers-and-polishers)
- Genome tools: [Dereplicators](#dereplicators)
- Alignment tools: [Phylogeny routes](#phylogeny-routes) (aligners, SNP
  typers and maskers), [Core-gene builders](#core-gene-builders)
- Tree tools: [Tree builders](#tree-builders)
- [CLI `run` or Nextflow](#cli-run-or-nextflow)
- [How `auto` chooses](#how-auto-chooses)
- [Declared limits](#declared-limits)

Terms used below: ANI is average nucleotide identity; an MSA is a multiple
sequence alignment; a representative is the genome kept for a cluster of
near-identical genomes after dereplication; the core genome is the set of
genes present in all or nearly all genomes of a set; the outgroup is a genome
outside the studied group, used to root the tree; the working directory is
the directory a CLI stage reads from and writes into (`-wd`).

## Measured and declared figures

The page distinguishes two kinds of figures:

- **Measured.** A time, a memory figure or a result from a recorded run, in
  [verification.md](verification.md) or
  [audit/scaling-audit.md](audit/scaling-audit.md). The times in
  [By number of genomes](#by-number-of-genomes) are from one machine (11
  cores, 18 GB). Treat them as orders of magnitude.
- **Declared.** The `recommended_max_genomes` value an adapter sets in its
  code. `auto` selection and the scale warnings use it. Most declared values
  have no recorded measurement behind them. Where the text gives a declared
  figure, it says so.

## Start from the question

Four properties of the dataset decide most choices:

1. **Dataset type.** Genomes from GTDB, local genome files, sequencing runs
   without an assembly, or viral sequences. This decides the
   [entry point](#entry-points).
2. **Number of genomes.** This decides the dereplicator and whether to
   dereplicate in chunks, and it bounds the tree builder
   ([By number of genomes](#by-number-of-genomes)).
3. **Taxonomic span.** A clonal (outbreak) set, one species, a genus, or a
   family or wider. The divergence decides whether a reference, a core
   genome, a set of marker genes or an alignment-free distance is
   appropriate.
4. **Supports.** Alignment-free tree builders give no branch supports. A
   maximum-likelihood (ML) tree with supports needs an alignment: a
   core-SNP, core-gene or whole-genome MSA.

## By purpose and divergence

This table combines the four properties. It is written for 100 to 500
genomes, the range where the whole-genome aligners no longer finish and the
choice of route matters most, and it follows the
[phylo methods audit](audit/phylo-methods-2026-10.md).
`repgenr phylo-advise -wd WD` reads a working directory and prints the row
that applies, the reasons and the `phylo` command
([phylogeny.md](phylogeny.md#asking-for-a-route-phylo-advise)). The Preset
column names the `--preset` of `phylo`, `phylo-build` and `run` that sets
the route and tree builder of the row; its flags and the override rule are
in [phylogeny.md](phylogeny.md#presets).

| Purpose and span | Entry path | Dereplicator | Phylogeny route | Masking | Tree builder | Preset |
|---|---|---|---|---|---|---|
| Outbreak or clonal genotype tree; one lineage, ANI above 99.5 % | `ingest`, or `reads` and `assemble` | `skder` or `galah` | `--msa-source snptype --snptyper ska2` (reference-free) | None; `ska2` cannot be masked | `iqtree` with `--aligner-arg model=GTR+G`, `-B 1000` | `outbreak` |
| Genotype tree within one species; ANI 95 to 99.5 % | `metadata`, or `ingest` for local genomes | `skder` or `galah` | `--msa-source snptype --snptyper simple` (medoid reference) or `parsnp`, `--core-fraction 0.95` | `--mask gubbins` when under 10 % of columns vary | `iqtree -B 1000` | `species`, or `species-masked` |
| Species representatives in a genus; ANI 80 to 95 % | `metadata`, or `--source api` for one taxon | `skder`, or `sourmash` for large sets | `annotate`, then `--msa-source coregene --coregene ppanggolin` (or `scarap`) | None | `iqtree -B 1000`, one model per gene | `genus` |
| Family or wider; ANI below about 80 %, or pairs too distant for a sketch ANI | `metadata` | `skder`, or `sourmash` for large sets | `annotate`, then `--msa-source coregene --coregene markers` | None | `iqtree -B 1000` on protein | `family` |
| Viral | `vmetadata` and `vgenome` | `skder` or `sourmash` | Alignment-free | None | `mashtree` | none |
| Quick look at any set | any | any | Alignment-free | None | `skani` (genus or wider), `mashtree` (diverse or viral), `ska2` (clonal) | `quick` |

Notes on the rows:

- **Outbreak.** The typing pass counts the constant sites of ska2's
  all-sites alignment and passes them to IQ-TREE as `-fconst`. Do not add
  `--aligner-arg lo=true` here: `ska lo` writes variants only, so the tree
  has no constant-site counts and `phylo` warns.
- **Species.** `species-masked` adds `--mask gubbins`. `phylo-advise`
  suggests it for up to 200 genomes with a minimum ANI of 97 % or more. The
  masker warns above 10 % variable columns and refuses above 12 %. `ska2`
  with `--aligner-arg lo=true` is the reference-free alternative when no
  masking is wanted.
- **Genus.** The preset keeps at most 300 core genes ranked by informative
  sites (`--max-core-genes`, `--gene-rank`) and fits one model per gene
  (`--partitioned`); `--coregene-alphabet protein` aligns protein instead.
  Panaroo suits this row only with `--aligner-arg threshold=` lowered to
  about 0.90: at its default of 0.98 no gene was in 95 % of 24 Wolbachia
  genomes. A `skani` tree is a quick check before the core-gene run.
- **Family.** The markers are GTDB's bac120 or ar53 single-copy proteins
  and need the marker database (`repgenr markers-db download`). IQ-TREE
  runs with `-st AA`, with `LG+G` on every partition under `--partitioned`.
  At this span DNA alignments saturate, orthology clustering finds few core
  genes, and mash and Jaccard distances are unreliable. `ppanggolin`,
  `scarap` and `orthofinder` with `--coregene-alphabet protein` remain
  available. `skani`, reliable to about 80 % ANI, gives an overview without
  supports.
- **Quick look.** There are no supports, so collapse the tree in tree2tax
  with `--collapse-diameter` rather than `--collapse-support`.

The evidence behind the rows:

- **Genus.** One `run` over 1157 GTDB Francisella genomes with `skder` and
  `mashtree` took 12 minutes and kept 68 representatives
  ([verification.md](verification.md#genus-scale-run-2026-09-10)). On those
  representatives the `ska2` tree shared 15 of 66 splits with the mashtree
  tree, and Gubbins failed on the 39 % variable genus alignment
  ([verification.md](verification.md#snp-phylogeny-on-the-genus-set-2026-09-10)).
  On 24 Wolbachia genomes the PPanGGOLiN, SCARAP and OrthoFinder trees
  shared 19 to 21 of 21 splits
  ([verification.md](verification.md#builder-comparison-on-24-genomes-2026-10-15)).
- **Species.** Gubbins ran at 1.0 and 8.3 % variable columns and failed
  from 13.2 % upward, on one genus only. On 9 F. tularensis genomes the
  masked and unmasked trees shared no internal split
  ([verification.md](verification.md#gubbins-divergence-limits-2026-09-10)).
- **Outbreak.** There is no measured run on a real outbreak. `ska2` is
  preferred for design reasons: it is reference-free, while mash distances
  saturate near identity. The dereplicators recovered the true clustering of
  synthetic clonal sets ([scaling audit](audit/scaling-audit.md#dereplication)).
- **Family.** On 24 Wolbachia genomes the markers tree was four splits from
  the PPanGGOLiN tree, all four with supports below 95
  ([verification.md](verification.md#gtdb-marker-genes-2026-10-21)).
- **Viral.** The live suite runs the viral chain on Hepatovirus records from
  NCBI Virus with `sourmash` and `mashtree` in 21 seconds. An earlier live
  run on 1256 Hepeviridae genomes from BV-BRC kept 799 representatives with
  `skder`
  ([verification.md](verification.md#viral-run-on-1256-hepeviridae-genomes-2026-06-15)).
  Other viral scaling is not measured.

Rules of thumb:

- A SNP route needs a genome every other genome maps to, so it stops being
  useful past about 5 % divergence.
- A core-gene route needs annotation (about 45 s per genome with bakta's
  light database) and works until the core collapses, which the builder's
  clustering setting governs (see [Core-gene builders](#core-gene-builders)).
- A whole-genome MSA is for sets under about 30 genomes.

With an alignment, the tree builder by leaf count:

- Under 500: `iqtree`.
- 500 to 1000: `iqtree` with `--aligner-arg model=GTR+G`, else `raxmlng`.
  The declared limit of 500 for `iqtree` is not measured. `--treebuilder
  auto` picks `raxmlng` here by its declared limits.
- 1000 to 5000: `fasttree`.
- Above 5000: `veryfasttree`.

For genus and family overviews built from hundreds of core genes, a species
tree summarised from gene trees (`--treebuilder aster` with `--msa-source
coregene --partitioned`) is an alternative to a concatenated tree. At that
span the gene trees disagree, and a concatenated ML tree can give high
support to a topology that many genes contradict. Short or nearly invariant
genes give poorly resolved gene trees, so keep the informative ones with
`--min-gene-length`, or `--gene-rank informative` with `--max-core-genes`.
It does not apply to SNP alignments or clonal sets. Its run time and its
agreement with the concatenated trees have not been measured here.

## By number of genomes

| Genomes | What was measured | What `auto` picks |
|---|---|---|
| Under 100 | `mashtree` 9 s at 100 genomes; sourmash tree 25 s at 100. progressiveMauve about 10 minutes at 20 genomes and about 22 minutes at 50, per build. | `iqtree` (limit 500), with or without an alignment request. |
| 100 to 1000 | `skder` about 7 min and 1.1 to 1.5 GB at 1000; `galah` 4 to 6 min; `sourmash` about 1 min and 0.2 GB. `mashtree` 117 s; sourmash tree 142 s. | Without an alignment request: `iqtree` up to 500, then `sourmash` (limit 2000). With one: `iqtree` up to 500, then `raxmlng` (limit 1000). |
| 1000 to 5000 | `skder` about 69 min and 6.4 to 8.2 GB at 5000 in one pass, about 15 min chunked; `galah` about 63 min; `sourmash` about 5.5 min sparse. `mashtree` 53 min and 7.6 GB at 5000. | Without an alignment request: `sourmash` up to 2000, then `mashtree` (limit 10000). With one: `fasttree` (limit 5000). |
| Over 5000 | Not measured. Single-pass `skder` or `galah` at 10000 would be memory-risky on the 18 GB audit machine. Chunked `skder` or sparse `sourmash` are the safer choices. | Without an alignment request: `mashtree` (limit 10000). With one: `veryfasttree` (limit 20000). |

The last column is the tree builder `--treebuilder auto` picks. Its limits
are declared (see [Declared limits](#declared-limits)). The count is of
leaves, the outgroup included, so the first switch happens at 501 leaves. An
alignment request is `--msa-source snptype` or `coregene`, `--mask`, or an
alignment supplied with `--msa`. Without one, the column assumes the default
aligner, progressiveMauve, whose limit is 500: past that `auto` chooses among
the alignment-free builders only. With `--aligner sibeliaz` or `cactus`
(limit 2000) that switch happens past 2000 genomes. `mashtree` and `skani`
share the limit 10000, and the tie goes to `mashtree` alphabetically. A
builder that cannot run here is passed over (rule 1 of
[How `auto` chooses](#how-auto-chooses)).

Notes on the table:

- **Most limits are declared, not measured.** The aligners and SNP typers
  have not been run at scale; SibeliaZ did not finish within an hour at 20
  genomes. The [status table](#declared-limits) gives the measurement
  behind each limit.
- **Chunks.** Past a few thousand genomes, dereplicate in chunks
  (`--process-size` on the CLI, `--derep_process_size` in Nextflow) to bound
  the time and memory of one tool call. Chunking changes which genome
  represents a cluster, not the clustering
  ([scaling audit](audit/scaling-audit.md#dereplication)).
- **The final pass.** The chunk size does not bound the final pass over the
  union of the chunk representatives. 6888 Bacillaceae genomes at
  `--secondary-ani 0.9995` kept 3807 after the chunks. A larger `-s`, or
  `--pre-secondary-ani` below the final threshold, is the better lever (see
  [Chunked dereplication](usage.md#chunked-dereplication); the probe rule for
  further levels is in [developing.md](developing.md#scaling)).

## Entry points

| Situation | Use |
|---|---|
| Bacteria or archaea in GTDB, many taxa or a whole release | `metadata` with the default `--source tsv`, which downloads the release table once. |
| One taxon only | `metadata --source api`, which asks the GTDB API for the target taxon and also fetches each genome's quality values. |
| Genomes already on disk | `ingest --genomes-dir DIR`. Add a `selection.tsv` for taxonomy and quality. |
| Sequencing runs, no assembly yet | `reads`, then `assemble`. Runs come from ENA, which mirrors SRA. |
| Viruses | `vmetadata` and `vgenome`. NCBI Virus is the default source. Use `--source bvbrc` only for the legacy BV-BRC path. |

`run` chains the stages for the GTDB path, and takes `--viral`, `--reads` or
`--genomes-dir` for the others (see [usage.md](usage.md)).

## Assemblers and polishers

These apply to the reads chain only (`reads`, `assemble`).

| Sequencing platform | Assembler (`--assembler auto`) | Polisher (`--polisher auto`) |
|---|---|---|
| Illumina | `skesa`. `shovill` (SPAdes) is the alternative for paired reads. | None |
| Oxford Nanopore | `flye` | `medaka` (`racon` when medaka is not installed) |
| PacBio CLR | `flye` | `racon` with minimap2 overlaps (`--polish-rounds`) |
| PacBio HiFi | `flye` | None |

Points that bear on the choice:

- Medaka needs the basecaller model. Reads mirrored through SRA lack it, so
  the adapter assumes ONT's bacterial R10.4.1 model and records that; give
  another with `--tool-arg model=...`.
- On one ONT run polishing did not change the CheckM2 scores, because the
  Flye draft was already gene-complete at that resolution
  ([verification.md](verification.md#polishing-an-ont-assembly-2026-09-16)).
  Unpolished long-read assemblies still carry indel errors that matter for
  SNP typing and alignment.
- `shovill` needs a read pair; `auto` assembles a single-file run with
  `skesa`.

Run selection (`--max-bases`, `--drop-selection`, `--one-per-sample`), the
organism check with `--classifier`, and the reuse of CheckM2 scores are
described in [reads.md](reads.md).

## Dereplicators

| Tool | Scaling | Representative chosen by the tool | Needs CheckM data | Verified |
|---|---|---|---|---|
| `skder` (default) | Superlinear: about 7 min at 1000 genomes and 69 min at 5000, up to 8.2 GB; about 15 min chunked at 5000. | Its own aggregate score, not quality-aware; which member of a clone block is arbitrary. | No | Native and container (Wave) |
| `galah` | Built for large sets: about 4 to 6 min at 1000 and 63 min at 5000. | galah's quality score when every genome has manifest quality, else the first listed genome, largest file first. | No; used when present | Native and container |
| `sourmash` | Sparse back-end close to linear in close pairs: about 5.5 min at 5000, 0.3 to 0.7 GB. Dense back-end capped at 5000. | The most-connected genome, alphabetical tie-break; biased toward the most-sequenced genotype. | No | Native and container |
| `drep` | Quadratic within primary clusters; declared limit 2000, chunk-wrapped. | Completeness, contamination, N50 and size; quality-aware, so the least sensitive to block size. | Yes, or manifest quality, or `--virus` | Container only |

Every dereplicator receives the genome set as a list file or a staged
directory, so the genome count does not limit the command line. A command
that would still exceed the operating system's argument-size limit (ARG_MAX)
stops before the tool starts, with a message that names `--process-size`.
skDER receives one staged directory of links. At a cutoff of 99 % or more
its adapter sets skani's screen to the cutoff minus one, so pairs below it
are not aligned (see [Chunked dereplication](usage.md#chunked-dereplication)).
galah, without quality values, lists the genomes by descending file size, so
a fragment is not preferred over a complete genome. The sparse and dense
sourmash back-ends use the same ANI estimate, and the one used is recorded.

dRep needs completeness and contamination for every genome of the run: from
the manifest (or `selection.tsv` for the Nextflow steps), passed as
`--genomeInfo`; else from CheckM on `PATH`. `--virus` passes
`--ignoreGenomeQuality` instead. Without any of these, dRep stops and
`dereplicate` exits 6. When some genomes lack manifest values, none are
passed and the log says so. The decision covers the whole run, so a chunked
run never scores one chunk from the manifest and another with CheckM.
Genomes that dRep's filter removes (`--tool-arg length=N`, or its
completeness and contamination thresholds) are `fail_qc` in
`genome_status.tsv` and in no cluster.

Which dereplicator to use:

- Use `skder`, the default, for most sets.
- Use `sourmash` with `--process-size 2000` (`--derep_process_size 2000` in
  Nextflow) for 10000 genomes or more (not measured at this size).
- Use `drep` when every genome has manifest quality (GTDB selections) or
  CheckM is available, and keep it under about 2000 genomes.

To keep one lineage of a species at a finer resolution than the rest, see
[Finer resolution inside a lineage](usage.md#finer-resolution-inside-a-lineage).

### Which genome represents a clone block

The tools choose the representative of a clone block by their own rules,
not by genome quality. The default `--keeper quality` re-picks each
cluster's representative by
`completeness - 5 x contamination + 0.5 x log10(N50)` (dRep's default
weights), whichever tool clustered; the ties and the rule for unscored
genomes are in [The keeper score](usage.md#the-keeper-score). skDER,
sourmash and galah then end with the same keepers when they form the same
clusters. `--keeper tool` keeps the adapter's choice. `--keeper gtdb` keeps
a GTDB species representative when the cluster holds one, applies the
quality rule otherwise, and needs a GTDB selection from `metadata`.

On the synthetic set `clonal_50_clustered` (groups of 20, 15 and 15
genomes, within-group ANI about 0.995 or higher, between-group about 0.96
or lower), `sourmash`, `galah` and `skder` at the defaults formed the same
three clusters, and each kept a different representative. To check this on
your own data, see
[Comparing two dereplications](usage.md#comparing-two-dereplications).

`repgenr glance` gives an all-against-all overview before choosing
thresholds. It runs on `drep` or `sourmash`, since the other dereplicators
have no comparison mode; `--tool auto`, the default, uses dRep when it can
run and sourmash otherwise. The sourmash glance needs no container on most
systems, plots the ANI estimate that `dereplicate --tool sourmash`
thresholds, and reuses the genome sketches of the working directory. It
holds the full N x N matrix and refuses sets above 5000 genomes. On one
workstation (8 threads, genomes on a USB disk), a glance of 1000 genomes
took about 20 s with the sketches present against 104 s when sketching, and
`dereplicate --tool sourmash` about 6 s against 65 s.

### How the tools estimate ANI

The thresholds are given as fractions (`--secondary-ani 0.99`) and each
adapter converts them to its tool's scale: dRep takes fractions as given,
skDER and galah take percentages (99), and sourmash, which has no ANI
cutoff option, has the cutoff mapped onto its sketch measure as described
below. The tools estimate ANI differently, and the estimates agree best at
the high end of the range where dereplication thresholds sit. dRep uses
fastANI for its secondary clustering (ANImf with `--virus`) and Mash for its
primary clustering. The two estimators compared here are skani, used by
skDER and galah, and sourmash.

skani (Shaw and Yu 2023) estimates ANI by a sparse approximate alignment.
It chains the k-mer seeds two genomes share into colinear blocks and
estimates identity from the shared k-mers within them, with a learned
correction. The value is an identity over mapped regions, as with
mapping-based ANI. The aligned fraction, the share
of each genome covered by chained blocks, is reported beside it, and
`--aligned-fraction` acts on it in skDER and galah. skani documents its
estimates as reliable down to about 80 % identity. Below that the adapter
refuses the cutoff for skDER, which would otherwise stop at a question;
galah has no such guard.

sourmash compares FracMinHash sketches (k=31, scaled=1000 by default). The
containment of one sketch in another is the fraction of its sampled k-mers
that the other holds. Under a model in which each base mutates
independently, without indels, a k-mer survives with probability ANI^k, so
an ANI estimate is containment^(1/k). The two directions of containment
differ when the genomes differ in size. Both back-ends use their average,
c = (|A & B|/|A| + |A & B|/|B|) / 2, and estimate ANI as c^(1/k). The
dense back-end reads it from `sourmash compare --avg-containment`, the
branchwater back-end from the `average_containment` column of
`sourmash scripts pairwise`, which it prefilters at containment 0.99^31
(about 0.73). The back-end used is recorded as `tool_settings` in the
stage record. There is no alignment and no aligned fraction, so
`--primary-ani` and `--aligned-fraction` are ignored. The conversion is
needed because raw similarity is on another scale: a pair at 99.57 % ANI by
skani has a k=31 Jaccard of 0.796 and a sourmash estimate of 99.61 %
([scaling-audit.md](audit/scaling-audit.md)), so a cutoff applied to the
raw value would leave nearly every genome a singleton.

What follows from the difference:

- At 99 % ANI and above the estimates agree well; for the pair above,
  99.61 % against skani's 99.57 %. The containment estimate is a sampling
  estimate, though: with `scaled=1000` a 4 Mb genome keeps about 4000
  k-mers, which gives a standard error of a few hundredths of a percent of
  ANI near 99 %, enough to move a borderline pair across a 99.9 % cutoff.
  skani samples more densely and is steadier there.
- Accessory content lowers the containment estimate at every identity,
  since plasmids and acquired regions add k-mers the other genome lacks,
  while skani reports identity within the mapped blocks and the rest as a
  lower aligned fraction. Below about 95 % ANI the sourmash estimate also
  loses precision: few k-mers survive (0.95^31 is about 0.20, 0.90^31 about
  0.04), so the variance grows and the mutation model fits less well.
- Completeness affects them differently. Containment is asymmetric (a
  fragment is contained in a complete genome, not the reverse), and either
  way of combining the two directions underestimates the ANI of a
  half-complete assembly to a complete genome. skani's identity is largely
  unaffected by completeness, which shows up in the aligned fraction
  instead.

Below 95 % ANI the sourmash estimate runs above skani's, by up to about 5
percentage points at k=31 and less at k=21 on 24 Wolbachia genomes; a
sourmash tree of a set wider than a species is therefore better built with
`--aligner-arg ksize=21`. The per-band figures are in the
[phylo methods audit](audit/phylo-methods-2026-10.md#sketch-ani-audit-2026-10-21).

A handful of shared hashes, from genes conserved across the set, still
gives a high estimate for pairs that have little else in common. Wherever
RepGenR2 converts a dense `sourmash compare` matrix (the dereplicator's
dense back-end, `glance`, the medoid reference, `phylo-advise` and the
sourmash tree builder with `metric=ani`), a pair with fewer than 30 shared
hashes and fewer than a tenth of the smaller sketch has an unknown ANI.
Such a pair never merges in dereplication, is an empty cell in the glance
table and distance 1.0 in the glance dendrogram and in a sourmash tree,
counts as 0 for the medoid, and counts as below 80 % in `phylo-advise`. The
relative part keeps small sketches usable: two identical sketches of a
7.5 kb viral genome share all of their about 8 hashes and keep an ANI of
1.0. For bacterial genomes of 1.0 to 1.5 Mb the sketch's reach ends near
88 to 90 % sketch ANI at k=31 and near 83 to 85 % at k=21; larger genomes
reach lower. The sparse branchwater back-end keeps only pairs above
containment 0.99^k, so the rule does not apply there. The measurements are
in the [phylo methods audit](audit/phylo-methods-2026-10.md#sketch-ani-audit-2026-10-21)
and in [verification.md](verification.md#sourmash-back-ends-on-24-wolbachia-genomes-2026-10-21).

The duplicate check of the reads screen (`duplicate_isolate`, see
[reads.md](reads.md#the-three-checks)) uses the same containment model on
reads sketches. From the k=31 signature, k-mers seen at least a tenth of the
abundance-weighted median (and at least twice) are taken as genomic. Their
containment C in the k-mers the other run saw at least twice gives the
estimate C^(1/31), after both sketches are reduced to the larger scaled
value. The default threshold of 0.999 corresponds to a containment of about
0.97. At scaled=1000 a 1 Mb genome gives about 1000 such k-mers, so the
estimate cannot separate 0.9999 from 1.0 for a genome of that size.

On well-separated data the tools usually form the same clusters at the
same cutoff. The differences appear in fragmented assemblies, in genomes of
unequal size and at the coarse end of the range. `glance` plots the
sourmash estimate (`--tool sourmash`) or dRep's Mash ANI (`--tool drep`) on
the scale a cutoff will be applied to.

## Phylogeny routes

`phylo` combines three independent choices: the genome set (representatives or
all genomes), the source of the alignment
(`--msa-source aligner|snptype|coregene`) and the tree builder
(`--treebuilder`). The alignment source is skipped for alignment-free tree
builders.

| Route | Tools | Use it for | Limits and costs |
|---|---|---|---|
| Alignment-free | `--treebuilder mashtree`, `sourmash`, `skani` or `ska2` | Diverse or large sets and quick trees; `skani` for a genus or wider, `ska2` for a clonal set. | Neighbour joining only, no supports. Mash and sourmash distances compress near identity; `ska2` undercounts past a few per cent divergence. |
| Whole-genome MSA | `--aligner progressivemauve`, `sibeliaz` or `cactus` | Small sets of related genomes. | Declared limits 500, 2000 and 2000. progressiveMauve about 10 to 22 minutes at 20 to 50 genomes; SibeliaZ timed out at 20. Cactus needs its container. |
| Core-gene | `--msa-source coregene` after `annotate`, with `--coregene panaroo`, `ppanggolin`, `scarap`, `orthofinder` or `markers` | A supported tree where a single reference would bias a SNP tree. | Genes in `--core-threshold` of the genomes, aligned with mafft and concatenated, invariant sites kept. Not measured at hundreds of genomes. |
| Core-SNP | `--msa-source snptype` with `--snptyper simple`, `snippy`, `parsnp` or `ska2` | One species or a clonal set. | `simple`, `snippy` and `parsnp` map to one reference, which biases the result; `ska2` is reference-free. |

How much of the core survives on the core-gene route depends on the
builder's clustering setting and on how wide the set is (see
[Core-gene builders](#core-gene-builders)); within a species a SNP tree
draws on far more positions. The `simple` typer maps with minimap2
`-x asm20` and records reference positions a genome does not cover as N, so
a divergent outgroup is largely N; its coverage is in the log.

### Reference bias

The aligners and the mapping typers project every genome onto one
reference. The choice of reference shapes mostly the within-cluster
branching order: swapping the reference changed 72 of 94 splits on a
50-genome set while the deep structure held
([scaling audit](audit/scaling-audit.md#phylogeny)).

How the reference is chosen:

- Unless `--reference` names a genome, the stage uses the medoid: the genome
  with the largest summed sketch ANI to the others. It is computed with
  sourmash from the sketches in the working directory (`--sketch` on the
  genome stages, or `repgenr sketch`); the choice never sketches. Ties go
  to the first filename, so the choice is the same on every run.
- The cached alignment is tied to the chosen genome, so it is rebuilt when
  the medoid changes.
- Under Nextflow, sketches reach `phylo-build` on the bacterial path with
  `params.sketch`. The viral and reads paths have no sketches and use the
  first genome.
- The medoid is expected to reduce the bias, not to remove it: in a set
  with an over-represented genotype it sits inside that genotype.
- Without sketches or sourmash, or with `--reference first`, the
  alphabetically first genome is used, with a warning. On the
  F. tularensis set that genome was another species, and 5.1 % of its
  columns varied inside a species usually described as clonal. For a
  species run, name a genome of that species.
- The chosen genome and the rule are recorded as `reference_genome` and
  `reference_rule`, and in `reference.txt` beside the alignment.

### Masking

`--mask gubbins` needs `--msa-source snptype` (in `run`, `--with-snptype` also
accepts `--mask`) and a typer that writes a whole-genome alignment, which
`simple`, `snippy` and `parsnp` do. It replaces the core-SNP alignment with
Gubbins' filtered polymorphic sites. `ska2` writes variable sites only, so it
cannot be masked. Gubbins is for isolates of one species. The masker
estimates the variable fraction of the alignment, warns above 10 %, refuses
above 12 %, and repeats the figure if Gubbins fails. Gubbins needs a
multi-threaded RAxML build when given more than one thread; see
[phylogeny.md](phylogeny.md#recombination-masking-with-gubbins).

## Core-gene builders

The core-gene route first annotates each genome (`annotate`, with `bakta`,
which needs its database, or `prokka`, which needs none). A core-gene
builder then groups the genes into orthologues and aligns each of those
present in at least `--core-threshold` of the genomes. The four clustering
builders differ in how they group genes, and that grouping, together with
how wide the set is, decides how many core genes a set retains. The fifth,
`markers`, takes a fixed set of 120 (or 53) marker genes and does not group
genes at all. Panaroo's graph cleaning (`--clean-mode`) and entropy filter
remove annotation errors and badly aligned genes; both have no effect with
the other builders.

| Builder | How it groups genes | Setting that governs the core | Measured on 24 Wolbachia genomes | Declared limit |
|---|---|---|---|---|
| `panaroo` (default) | Clusters genes at a fixed identity, then cleans the pangenome graph. | `--aligner-arg threshold=`, 0.98 by default and 0.95 when the names span more than one species. | Core of 52 genes at 0.98 and 638 at 0.90 on three genomes; none across supergroups at the default. | 500, not measured |
| `ppanggolin` | Clusters proteins with MMseqs2 at 80 % identity over 80 % of their length, then partitions the graph. | `identity=`, `coverage=` and `partition=`. | 398 single-copy families. | 5000, not measured |
| `scarap` | Builds orthogroups on seed genomes, then searches every genome against their profiles. | No identity cut; `method=` and `seeds=`. | 577 genes; 20 of 21 splits shared with PPanGGOLiN. | 10000, not measured |
| `orthofinder` | All-against-all DIAMOND search, clustered with MCL. | `search=` and `inflation=`. | 585 genes; 20 of 21 splits shared with PPanGGOLiN. | 200, not measured |
| `markers` | No clustering: GTDB's bac120 or ar53 marker families found with HMMER. Protein only. | `--marker-set`, `--core-threshold` and `min_markers=`. | 99 of 120 markers kept; four splits from PPanGGOLiN. | 20000, not measured |

The measurements are in [verification.md](verification.md#core-gene-route).
`panaroo` and `ppanggolin` write their per-gene alignments themselves
(mafft). For `orthofinder` and `scarap` the adapter aligns each
orthogroup's proteins with mafft and back-translates them to codons, and
`markers` aligns each marker's proteins to its HMM with `hmmalign`. The
stage then ranks, caps and concatenates the genes in the same way for all
five builders.

Trimming (`--gene-trim clipkit`) runs ClipKIT on each per-gene alignment
before the genes are measured and ranked. It may help on genus and wider
sets, where partial genes leave ragged ends and poorly aligned columns that
the ranking would read as signal. Within a species it is expected to change
little. It is off by default, and turning it on or off rebuilds the
core-gene alignment. Its effect on tree topology has not been measured.

Which builder to use:

- `panaroo` is the default. For a wider set, lower `--aligner-arg
  threshold=` in steps; lowering `--core-threshold` instead recovers fewer
  genes and gives a softer core.
- `ppanggolin` is a candidate for a divergent set.
- `scarap` infers the core through profile searches, so there is no
  identity to lower.
- `orthofinder` runs an all-against-all search whose cost grows roughly
  with the square of the genome count, so it suits small and medium sets.
- `markers` suits a set that spans a family or more. Its search is linear
  in the genome count. It gives at most 120 genes, so within a genus the
  clustering builders draw on more of the genome.

The settings and caveats of each builder are in
[phylogeny.md](phylogeny.md#core-gene-builders).

## Tree builders

| Builder | Input | Declared scale | Resolution and supports |
|---|---|---|---|
| `iqtree` | MSA | 500; poor beyond that. | Good on diverse and clonal sets. Supports with `-B 1000` or `--aligner-arg alrt=1000`; with neither the tree has none. |
| `raxmlng` | MSA | 1000; slow beyond that. | Supports from its own bootstrap procedure. |
| `fasttree` | MSA | 5000. | Approximate; weaker on clonal sets. SH-like local supports. |
| `veryfasttree` | MSA | 20000, not measured; `fasttree`'s algorithm on several cores. | As `fasttree`. |
| `mashtree` | Genomes | 10000; 53 min and 7.6 GB at 5000. | Poor on clonal sets. No supports. |
| `sourmash` | Genomes | 2000. | Poor on clonal sets. No supports. |
| `skani` | Genomes | 10000, not measured; ANI reliable to about 80 %. | Expected to resolve closer pairs than mash distance and more distant ones than Jaccard, not measured here. No supports. |
| `ska2` | Genomes | 5000. | Good on clonal and outbreak sets, reference-free; undercounts once pairs share few k-mers. No supports. |
| `aster` | MSA with gene partitions | 5000, not measured. | A species tree from gene trees that allows for their discordance; local posterior probabilities as supports. |

`iqtree` and `raxmlng` take the constant-site counts of a SNP alignment and
the core-gene partitions (`--partitioned`). Without `--aligner-arg model=`,
`iqtree` runs ModelFinder first; `model=GTR+G` skips it. `aster` needs
`--msa-source coregene --partitioned`, and its coalescent branch lengths
refuse `--root midpoint`. The options of each builder are in
[phylogeny.md](phylogeny.md#tree-builder-options).

### Automatic choice of tree builder

`--treebuilder auto` picks the builder whose declared scale fits the leaf
count, the outgroup included. A request that only an alignment can honour
(`--msa-source snptype` or `coregene`, or `--mask`) limits the choice to
builders that consume one. `auto` never picks `aster` or `ska2`. The full
rule is in [How `auto` chooses](#how-auto-chooses).

### Rooting

`phylo` roots `tree/tree.nwk` itself, whatever the builder: on the outgroup
by default when one is staged, at the midpoint with `--root midpoint`, or not
at all with `--root none`. Details are in
[Rooting the tree](phylogeny.md#rooting-the-tree).

### Collapsing weak splits in tree2tax

`--collapse-support` needs supports: IQ-TREE writes them only with `-B` or
`--aligner-arg alrt=N`, FastTree and RAxML-NG write them by default, and the
alignment-free builders write none. `--collapse-length` and
`--collapse-diameter` use branch lengths and so work on every builder; the
diameter maps a distance cutoff, for example about 95 % ANI, onto the
taxonomy without supports. The units of each builder and examples are in
[Collapsing weak splits in tree2tax](usage.md#collapsing-weak-splits-in-tree2tax).

## CLI `run` or Nextflow

| | `repgenr run` | `nextflow run nextflow/main.nf` |
|---|---|---|
| Best for | One machine, up to about a thousand genomes (measured). | Scatter-gather dereplication above a few thousand genomes (a design recommendation, not a measurement), HPC, cloud. |
| Parallelism | Threads within a stage. | Tasks per chunk and per assembled run. |
| Resume | Per stage, from parameter and input fingerprints. `--force` overrides. | Nextflow's `-resume` and task cache. |
| Representative choice | `--keeper quality\|gtdb\|tool` | `--derep_keeper quality\|gtdb\|tool` |
| Options | Forwards the stage options it shares a name with. | `--metadata_args`, `--phylo_args` and the other `*_args` strings, plus typed `--derep_*` parameters. |
| Inspection | `status` and `doctor` are available. | Not available. They are CLI commands that read a working directory. |

The Nextflow layer has no shared working directory, so `status` and `doctor`
have nothing to read there. Parameters and profiles are in
[nextflow.md](nextflow.md).

## How `auto` chooses

Each adapter declares `recommended_max_genomes`, an integer or none.
`--tool auto` (dereplicate) and `--treebuilder auto` (phylo) use it as
follows:

1. Prefer tools that can run here. Under a container backend, a pinned image
   counts as available, and so does a conda spec with `--wave`. A tool with
   neither runs on the host, so its binaries must be on `PATH`, as they must
   be natively.
2. Prefer tools whose limit fits the number of genomes. Among those, the
   tightest limit wins, so the more careful but limited tool is chosen for
   small inputs. Tools without a limit come last among the fitting ones.
3. Break ties with the documented defaults (`skder`, `iqtree`,
   `progressivemauve`, `simple`, `skesa`, `shovill`, `flye`, `medaka`,
   `racon`), then alphabetically.
4. If nothing fits, choose the tool with the largest limit.
5. For `--treebuilder auto` only: when no alignment was requested (no
   `--msa-source snptype` or `coregene`, `--mask` or `--msa`) and the leaf
   count exceeds the selected aligner's limit, steps 1 to 4 run over the
   alignment-free builders alone.

The default for `--tool` on `dereplicate` is `skder`, not `auto`. With `auto`
and all tools available, `drep` (limit 2000) would be preferred for up to 2000
genomes, and it needs manifest quality for every genome or CheckM. Under
`--container docker` or `singularity` dRep always counts as available, because
its image is pinned, so `auto` chooses it for up to 2000 genomes. Without
manifest quality, set `CHECKM_DATA_PATH` on the host (see
[install.md](install.md#databases)) or name another tool.

Passing a tool by name past its limit logs a warning and runs anyway. Adapter
authors set the limit in `ToolCapabilities` (see
[developing.md](developing.md)).

## Declared limits

`repgenr list-tools` prints the declared limits. The table below is checked
against the registered adapters by `tests/unit/test_docs_tool_limits.py`. It
lists what each adapter declares, not what was measured.

| Family | Tool | Declared limit |
|---|---|---|
| dereplicator | drep | 2000 |
| dereplicator | galah | unbounded |
| dereplicator | skder | unbounded |
| dereplicator | sourmash | unbounded |
| aligner | cactus | 2000 |
| aligner | progressivemauve | 500 |
| aligner | sibeliaz | 2000 |
| snptyper | parsnp | 2000 |
| snptyper | simple | 2000 |
| snptyper | ska2 | 5000 |
| snptyper | snippy | 1000 |
| masker | gubbins | unbounded |
| treebuilder | fasttree | 5000 |
| treebuilder | veryfasttree | 20000 |
| treebuilder | iqtree | 500 |
| treebuilder | mashtree | 10000 |
| treebuilder | raxmlng | 1000 |
| treebuilder | skani | 10000 |
| treebuilder | ska2 | 5000 |
| treebuilder | sourmash | 2000 |
| treebuilder | aster | 5000 |
| assembler | flye | unbounded |
| assembler | shovill | unbounded |
| assembler | skesa | unbounded |
| classifier | sourmash | unbounded |
| polisher | medaka | unbounded |
| polisher | racon | unbounded |
| annotator | bakta | unbounded |
| annotator | prokka | unbounded |
| coregene | panaroo | 500 |
| coregene | ppanggolin | 5000 |
| coregene | scarap | 10000 |
| coregene | orthofinder | 200 |
| coregene | markers | 20000 |

How the declared limits relate to the measurements on this page:

| Tool | Status of the declared limit |
|---|---|
| tree builder `sourmash` (2000) | Lowered after the scaling-audit measurement; kept until the rapidNJ path is measured. |
| tree builder `mashtree` (10000) | Measured at 5000 genomes (53 min, 7.6 GB), not above. |
| dereplicators `skder`, `galah`, `sourmash` (unbounded) | Measured up to 5000 genomes in the scaling audit; chunked skder also ran on 6888 Bacillaceae genomes (see [Chunked dereplication](usage.md#chunked-dereplication)). galah and sourmash not measured above 5000. |
| aligner `progressivemauve` (500) | Measured at 20 and 50 genomes only. |
| aligner `sibeliaz` (2000) | Not supported by a run: SibeliaZ did not finish within an hour at 20 genomes. |
| other aligners, SNP typers | Not run at scale. |
| `veryfasttree`, `skani`, `aster`, the core-gene builders | Not measured here. |

These limits count genomes only. They do not model genome size, memory or
thread count. Aligners and SNP typers declare limits that the scale warnings
use, but limit-based `auto` selection applies only to dereplicators and tree
builders.
