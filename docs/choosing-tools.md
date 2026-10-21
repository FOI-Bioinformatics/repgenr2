# Choosing tools

RepGenR2 has ten tool families in four groups: read tools (assemblers,
polishers, classifiers), genome tools (dereplicators, annotators), alignment
tools (aligners, SNP typers, maskers, core-gene builders) and tree tools (tree
builders). Most stages work with their defaults. This page is for users
who need to decide what to use when the defaults do not fit. It starts from
four questions about the dataset, gives a decision table, then describes the
trade-offs of each family, and ends with the genome limit each adapter
declares.

Contents:

- [Measured and declared figures](#measured-and-declared-figures)
- [Start from the question](#start-from-the-question)
- [Decision table](#decision-table)
- [By purpose and divergence](#by-purpose-and-divergence)
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
- [Which environment provides each tool](#which-environment-provides-each-tool)

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
   without an assembly, or viral sequences. This decides the entry point
   ([Entry points](#entry-points)).
2. **Number of genomes.** This decides the dereplicator and whether to
   dereplicate in chunks, and it bounds the tree builder
   ([By number of genomes](#by-number-of-genomes)).
3. **Taxonomic span.** A genus or wider, one species, or a clonal (outbreak)
   set. This decides the phylogeny route: how divergent the genomes are
   determines whether a reference, a core genome or an alignment-free
   distance is appropriate.
4. **Supports or a model-based tree.** Alignment-free tree builders give no
   branch supports. A maximum-likelihood (ML) tree with supports needs an
   alignment: a core-SNP, core-gene or whole-genome MSA.

The flow below combines questions 3 and 4 for the phylogeny route. Each
branch is described in [Phylogeny routes](#phylogeny-routes).

```
Need branch supports or a model-based (ML) tree?
  no  -> alignment-free tree builder: mashtree or sourmash for diverse or
         large sets; ska2 for clonal sets; skani
  yes -> Within one species?
           yes -> Clonal or outbreak set?
                    yes -> core-SNP with ska2 (reference-free)
                    no  -> core-SNP with a mapping typer and --reference,
                           optionally --mask gubbins;
                           or a whole-genome MSA for small sets
           no  -> core-gene route (annotate, then a core-gene builder);
                  the wider the set, the lower the clustering identity it
                  needs, or use ppanggolin, scarap or orthofinder;
                  across a family or wider, the GTDB markers (markers)
```

## Decision table

| Dataset | Entry path | Dereplicator | Phylogeny route | Masking | Tree builder |
|---|---|---|---|---|---|
| Genus or wider (diverse) | `metadata` (GTDB), or `--source api` for one taxon | `skder`, or `sourmash` for large sets | Alignment-free | No | `mashtree` |
| One species | `metadata`, or `ingest` for local genomes | `skder` or `galah` | Core-SNP with `--reference` set to a genome of that species, or whole-genome MSA for small sets | `--mask gubbins` when under 10% of alignment columns vary | `iqtree` |
| Outbreak or clonal set | `ingest`, or `reads` and `assemble` | `skder` or `galah`, with the clone-block caveat in [Dereplicators](#dereplicators) | Core-SNP with `ska2` (reference-free) or a mapping typer | `gubbins` with a mapping typer only; `ska2` cannot be masked | `iqtree` or `fasttree` |
| Viral | `vmetadata` and `vgenome` (NCBI Virus) | `skder` or `sourmash` | Alignment-free | No | `mashtree` |

The core-gene route has no row of its own. It applies where a supported,
model-based tree is wanted and a single reference would bias a SNP tree (see
[Phylogeny routes](#phylogeny-routes) and
[Core-gene builders](#core-gene-builders)).

The evidence behind each row:

- **Genus.** One `run` over all 1157 GTDB Francisella genomes with `skder` and
  `mashtree` took 12 minutes end to end and gave 68 representatives (genus
  run, 2026-09-10, in `verification.md`). On the 68 representatives plus the
  outgroup, `ska2` with IQ-TREE took 7.4 minutes and gave a tree that shared 15
  of 66 splits with the mashtree tree. `simple` with Gubbins failed in
  Gubbins' recombination scan on the genus alignment, which is 39% variable.
- **One species.** Gubbins ran on 9 F. tularensis genomes (1.0% variable
  columns, 1.9 minutes) and on 18 F. philomiragia genomes (8.3%, 4.4 minutes).
  It failed from 13.2% variable columns upward. Masking changed the topology:
  the masked and unmasked IQ-TREE trees on 9 representatives share no internal
  split. These are the only species-level measurements, and they cover one
  genus.
- **Outbreak or clonal.** The dereplicators were measured on synthetic clonal
  sets in the scaling audit (all recovered the true clustering). The SNP typers
  were run on a synthetic within-species set and on the Francisella data
  (`verification.md`). There is no measured run on a real outbreak. The
  reasons to prefer `ska2` here are design reasons from the audit. It is
  reference-free, so no reference assembly biases the distances. Mash
  distances saturate near identity, which makes `mashtree` weak on clonal sets.
- **Viral.** The live suite runs the viral chain on Hepatovirus records from
  NCBI Virus with `sourmash` and `mashtree` (`test_run_viral_chain_end_to_end`,
  21 seconds). The largest viral run recorded is 1256 Hepeviridae (hepatitis E)
  genomes from BV-BRC: `skder` gave 799 representatives, followed by `mashtree`
  and `tree2tax` (recorded in `docs/verification.md` before #150, commit
  7309fe7; the scaling audit keeps the 1256 figure). Other viral scaling is
  not measured.

## By purpose and divergence

The decision table above starts from the dataset. This table starts from the
purpose of the tree and the ANI span of the set, for 100 to 500 genomes, the
range where the whole-genome aligners no longer finish and the choice of route
matters most. It follows the
[phylo methods audit](audit/phylo-methods-2026-10.md).
`repgenr phylo-advise -wd WD` reads a working directory and prints the row
that applies, the reasons and the `phylo` command
([usage](usage.md#asking-for-a-route-phylo-advise)). The Preset column names
the `--preset` of `phylo`, `phylo-build` and `run` that sets the route and
tree builder of the row; its flags and the override rule are listed in
[usage](usage.md#presets).

| Purpose | Preset | Span, ANI | Route for 100 to 500 genomes | Tree builder | Notes |
|---|---|---|---|---|---|
| Outbreak or clonal genotype tree | `outbreak` | One lineage, ANI above 99.5 % | `--msa-source snptype --snptyper ska2` (reference-free) | `iqtree` with `--aligner-arg model=GTR+G` and `-B 1000` (constant sites from ska2's all-sites alignment passed as `-fconst` automatically) | Not with `--aligner-arg lo=true`: `ska lo` writes variants only, so that tree has no constant-site counts and `phylo` warns. The `ska2` tree builder gives a quick look. |
| Genotype tree within a species | `species`, or `species-masked` with `--mask gubbins` | One species, ANI 95 to 99.5 % | `--msa-source snptype --snptyper simple` (medoid reference) or `parsnp`; `--mask gubbins` when recombination is expected and under 12 % of columns vary; `--core-fraction 0.95` | `iqtree` with `-B 1000` (constant sites passed as `-fconst` automatically); `fasttree` as Gubbins' internal builder | `ska2` with `--aligner-arg lo=true` is the reference-free alternative when no masking is wanted. |
| Overview of species representatives in a genus | `genus` (PPanGGOLiN, at most 300 genes ranked by informative sites, `--partitioned`) | Several species, ANI 80 to 95 % | `annotate`, then `--msa-source coregene --coregene ppanggolin` (or `scarap`), `--core-threshold 0.95` | `iqtree -B 1000`, one model per gene with `--partitioned`; a `skani` tree first as a quick check | A gene cap with ranking (`--max-core-genes`, `--gene-rank`), a partitioned model (`--partitioned`) and a protein alignment (`--coregene-alphabet protein`) are available. Panaroo only with `--aligner-arg threshold=` lowered to about 0.90; at 0.98 (Panaroo's own default) no gene was in 95 % of 24 Wolbachia genomes. |
| Overview across a family or wider | `family` (`genus` with `--coregene markers`) | ANI 75 to 85 % or below, or pairs too distant for a sketch ANI | `annotate`, then `--msa-source coregene --coregene markers` (GTDB's bac120 or ar53 single-copy marker proteins, aligned as protein; needs the marker database, `repgenr markers-db download`); `ppanggolin`, `scarap` or `orthofinder` with `--coregene-alphabet protein` remain available | `iqtree` (run with `-st AA`; `LG+G` on every partition with `--partitioned`, ModelFinder without it, or `--aligner-arg model=` to choose) and `-B 1000` | DNA alignments saturate, orthology clustering finds few core genes, and mash and Jaccard distances are unreliable at this span. The markers need no clustering and are the genes of GTDB's reference tree. `skani` (reliable to about 80 % ANI) gives a quick overview without supports. |
| Quick look at any set | `quick` (`skani`) | Any | Alignment-free | `skani` (genus or wider), `mashtree` (diverse or viral), `ska2` (clonal) | No supports, so use `--collapse-diameter` in tree2tax rather than `--collapse-support`. |

Rules of thumb:

- A SNP route needs a genome every other genome maps to, so it stops being
  useful past about 5 % divergence.
- A core-gene route needs annotation (about 45 s per genome with bakta's light
  database) and works until the core collapses, which the builder's
  clustering setting governs (see [Core-gene builders](#core-gene-builders)).
- A whole-genome MSA is for sets under about 30 genomes.

With an alignment, the tree builder by leaf count:

- Under 500: `iqtree`.
- 500 to 1000: `iqtree` with `--aligner-arg model=GTR+G`, else `raxmlng`.
  The declared limit of 500 for `iqtree` is not measured. This differs from
  `--treebuilder auto`, which picks `raxmlng` here by its declared limits.
- 1000 to 5000: `fasttree`.
- Above 5000: `veryfasttree`.

A species tree summarised from gene trees (`--treebuilder aster` with
`--msa-source coregene --partitioned`) is an alternative to a concatenated
tree for genus and family overviews built from hundreds of core genes. At
that span the gene trees disagree with each other, and a concatenated
maximum-likelihood tree can give high support to a topology that many genes
contradict; a summary method such as ASTER takes the discordance into
account. Each gene tree is estimated from one gene only, so short or nearly
invariant genes give poorly resolved gene trees; `--min-gene-length` or
`--gene-rank informative` with `--max-core-genes` keeps the more informative
genes. It does not apply to SNP alignments, which have no genes, or to
clonal sets, where most gene trees are unresolved. Its run time and its
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
the alignment-free builders only. With `--aligner
sibeliaz` or `cactus` (limit 2000) that switch happens past 2000 genomes.
`mashtree` and `skani` share the limit 10000, and the tie goes to `mashtree`
alphabetically. A builder that cannot run here is passed over (rule 1 of
[How `auto` chooses](#how-auto-chooses)).

Notes on the table:

- **Most limits are declared, not benchmarked.** The audit found the declared
  `recommended_max_genomes` values to be round numbers with no recorded
  measurement behind them (500, 1000, 2000, 5000, 10000). The scaling audit
  later measured some of them (for example the sourmash tree builder, whose
  limit was lowered after the measurement). The aligners and SNP typers have
  not been run at scale. In the audit, SibeliaZ did not finish within an hour
  at 20 genomes, so its declared limit of 2000 is not supported by a run. The
  largest recorded end-to-end runs are the 1157-genome Francisella set
  (bacterial) and 1256 Hepeviridae genomes (viral).
- The sourmash tree builder joins its distance matrix with rapidNJ when the
  `rapidnj` binary is installed (`envs/core.yml` has it), and with a
  pure-Python step otherwise. A fit to the measured times of the Python step
  gives about 1.6 hours at 5000 genomes, so without rapidNJ more than 5000
  leaves are refused. With rapidNJ, its authors report that the joining takes
  seconds at that size (not measured here), and the dense `sourmash compare`
  matrix becomes the limit. The declared limit stays 2000 until that is
  measured. `--aligner-arg nj=python` or `nj=rapidnj` forces one backend.
- Past a few thousand genomes, dereplicate in chunks (`--process-size` on the
  CLI, `--derep_process_size` in Nextflow) to bound the time and memory of one
  tool call. Chunking changes which genome represents a cluster, not the
  clustering (scaling audit, finding 3).
- The chunk size bounds the leaf chunks only. The final pass compares the
  union of the chunk representatives, and for a diverse set at a strict
  threshold that union can be most of the input: 6888 Bacillaceae genomes at
  `--secondary-ani 0.9995` kept 3807 after the chunks. A level recomputes the
  per-genome work for every representative and may remove few, so one chunk
  of the next level is probed before that level runs. The union goes to the
  final pass when the probe removes less than the level is expected to cost.
  The final pass then still compares about a third of the single-pass pairs,
  and the per-genome work is paid twice. A larger `-s`, or
  `--pre-secondary-ani` below the final threshold, is the better lever (see
  [Chunked dereplication](usage.md#chunked-dereplication)).
- With skDER, the adapter sets skani's screen to one point below a cutoff of
  99 or more (`--tool-arg screen=` overrides), so most same-species pairs are
  not aligned.
- The genome count does not limit a single call through the command line.
  Every dereplicator receives the genome set as a list file or a staged
  directory. A command that would still exceed the operating system's
  argument-size limit (ARG_MAX) stops before the tool starts, with a message
  that names `--process-size`.

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
- On one ONT run (SRR28800588) polishing did not change the CheckM2 scores,
  because the Flye draft was already gene-complete at that resolution.
  Unpolished long-read assemblies still carry indel errors that matter for SNP
  typing and alignment.
- `shovill` needs a read pair; `auto` assembles a single-file run with
  `skesa`.

Run selection (`--max-bases`, `--drop-selection`, `--one-per-sample`), the
organism check with `--classifier`, and the reuse of CheckM2 scores are
described in [Starting from sequencing reads](usage.md#starting-from-sequencing-reads).

## Dereplicators

| Tool | Scaling shape | Representative criterion | Clone-block behaviour | Needs CheckM data | Verified |
|---|---|---|---|---|---|
| `skder` (default) | Superlinear in practice: 5 times the genomes cost about 10 times the wall time and 6 times the memory. Single pass: about 7 min at 1000 and 69 min at 5000, memory up to 8.2 GB. Chunked: about 15 min at 5000 (measured before the skani screen default below). Receives the genome set as one staged directory of links, so the set size is not bounded by the command line. At a cutoff of 99 percent or more, skani's screen is set to the cutoff minus one (`--tool-arg screen=` overrides), so pairs below it are not aligned. | Its own aggregate score. Not quality-aware. | One representative per block. Which member is arbitrary. | No | Native and container (Wave) |
| `galah` | Built for large sets. About 4 to 6 min at 1000 and 63 min at 5000. | With manifest quality for every genome of the run, galah's quality score (passed as `--genome-info`). Otherwise input position, with the genomes listed by descending file size, so a fragment is not preferred over a complete genome. | Without quality, the first listed member in the three orderings tested. | No; used when present | Native and container |
| `sourmash` | Sparse back-end close to linear in close pairs. About 5.5 min at 5000, 0.3 to 0.7 GB. Dense back-end is capped at 5000. | Most-connected genome, alphabetical tie-break. | Biased toward the most-sequenced genotype. The sparse and dense back-ends use the same ANI estimate (average containment); the one used is recorded. | No | Native and container |
| `drep` | Quadratic within primary clusters. Declared limit 2000, chunk-wrapped. | Completeness, contamination, N50 and size score. Quality-aware. | Best-scored member, so least sensitive to block size. | Yes: completeness and contamination in the manifest (or `selection.tsv` for the Nextflow steps) for every genome of the run, passed as `--genomeInfo`; else CheckM on `PATH`; or `--virus`, which passes `--ignoreGenomeQuality`. Without any of these, dRep stops and `dereplicate` exits 6. | Container only |

To keep one lineage of a species at a finer resolution than the rest,
whichever dereplicator is used, see "Finer resolution inside a lineage" in
[usage.md](usage.md#finer-resolution-inside-a-lineage).

Which dereplicator to use:

- Use `skder`, the default, for most sets.
- Use `sourmash` with `--derep_process_size 2000` (Nextflow) or
  `--process-size 2000` for 10000 genomes or more (not measured at this size).
- Use `drep` when every genome has manifest quality (GTDB selections) or
  CheckM is available, and keep it under about 2000 genomes. With manifest
  quality for every genome, dRep scores genomes with those values instead of
  running CheckM, the same values `--keeper quality` uses. When some genomes
  lack them, the values are not passed (dRep needs all of them) and the log
  says so. The decision covers the whole run, so a chunked run never scores
  one chunk from the manifest and another with CheckM. Genomes that dRep's
  filter removes (`--tool-arg length=N`, or its completeness and
  contamination thresholds) are `fail_qc` in `genome_status.tsv` and in no
  cluster.

Which genome represents a clone block:

- The representative inside a clone block depends on the tool, its back-end
  and the accession names, not on genome quality. The default
  `--keeper quality` re-picks each cluster's representative by
  `completeness - 5 x contamination + 0.5 x log10(N50)` (dRep's default
  weights), using the values in the manifest and the N50 of each genome file,
  whichever tool clustered. Ties go to higher completeness, then lower
  contamination, then higher N50, then the filename, so skDER, sourmash and
  galah end with the same keepers when they form the same clusters. An
  unscored representative is replaced only by a high-quality genome
  (completeness above 90, contamination below 5).
- `--keeper tool` keeps the adapter's choice.
- `--keeper gtdb` keeps a GTDB species representative when the cluster holds
  one and applies the quality rule otherwise, so that the representatives
  coincide with GTDB's where the clustering allows. It needs a GTDB selection
  from `metadata`.
- On the synthetic set `clonal_50_clustered` (groups of 20, 15 and 15 genomes,
  within-group ANI about 0.995 or higher, between-group about 0.96 or lower),
  `sourmash`, `galah` and `skder` at the defaults formed the same three
  clusters, and each kept a different representative. To check this on your own
  data, see "Comparing two dereplications" in
  [usage.md](usage.md#inspecting-a-dereplication).

An overview before choosing thresholds:

- `repgenr glance` (an all-against-all overview before choosing thresholds)
  runs on `drep` or `sourmash`; the other dereplicators have no comparison
  mode. `--tool auto`, the default, uses dRep when it can run and sourmash
  otherwise.
- sourmash needs no container on most systems and plots the ANI estimate that
  `dereplicate --tool sourmash` thresholds. It holds the full N x N matrix and
  refuses sets above 5000 genomes. It reads the genome sketches of the
  working directory, so it does not sketch again when they exist.
- On one workstation (8 threads, genomes on a USB disk), a sourmash glance of
  `clonal_50_clustered` took about 4 s with the sketches against 7 s
  sketching, and of `balanced_1000_clustered` about 20 s against 104 s.
  `dereplicate --tool sourmash` on the 1000 genomes took about 6 s against
  65 s.

### How the tools estimate ANI

The thresholds are given as fractions (`--secondary-ani 0.99`) and each
adapter converts them to its tool's scale: dRep takes fractions as given,
skDER and galah take percentages (99), and sourmash, which has no ANI
cutoff option, has the cutoff mapped onto its sketch measure as described
below. The same nominal cutoff therefore reaches every tool, but the tools
estimate ANI differently, and the estimates agree best at the high end of
the range where dereplication thresholds sit. dRep estimates its secondary
ANI with fastANI by default (ANImf with `--virus`) and its primary
clustering with Mash; the two estimators compared here are those of skDER
and galah on one side and sourmash on the other.

skani (used by skDER and galah; Shaw and Yu 2023) estimates ANI by a sparse
approximate alignment. It samples k-mers from each genome, finds the seeds
two genomes share and chains them into colinear blocks, then estimates
identity from the fraction of shared k-mers within the chained regions,
with a learned correction applied by default. The value is an identity over
mapped regions, as with mapping-based ANI, and the aligned fraction, the
share of each genome covered by chained blocks, is reported beside it; this
is what `--aligned-fraction` acts on in skDER and galah. skani documents its
estimates as reliable down to about 80% identity; below that the adapter
refuses the cutoff for skDER (which would otherwise stop at a question),
while galah has no such guard.

sourmash compares FracMinHash sketches (k=31, scaled=1000 by default) and
stops at set overlap: the containment of one sketch in another is the
fraction of its sampled k-mers that the other holds. Under a model in which
each base mutates independently, without indels, a k-mer survives with
probability ANI^k, so an ANI estimate is containment^(1/k). The two
directions of containment differ when the genomes differ in size; both
back-ends use their average, c = (|A & B|/|A| + |A & B|/|B|) / 2, and
estimate ANI as c^(1/k). The dense back-end reads it from `sourmash compare
--avg-containment`, the branchwater back-end from the `average_containment`
column of `sourmash scripts pairwise` (which it prefilters at containment
0.99^31, about 0.73). Until 2026-10-21 the dense back-end converted the
Jaccard index as 2j/(1+j), the harmonic mean of the two containments, which
is lower than the average for genomes of unequal size, so a clustering
could depend on whether the branchwater plugin was installed; it was
switched so that both back-ends use one quantity. The back-end used is
recorded as `tool_settings` in the stage record. There is no alignment and
no aligned fraction, so `--primary-ani` and `--aligned-fraction` are
ignored. The conversion is needed because raw
similarity is on another scale: a pair at 99.57% ANI by skani has a k=31
Jaccard of 0.796 and a sourmash estimate of 99.61%
([scaling-audit.md](audit/scaling-audit.md)), so a cutoff applied to the
raw value would leave nearly every genome a singleton.

What follows from the difference:

- At 99% ANI and above the estimates agree well; for the synthetic pair
  above, 99.61% against skani's 99.57%. The containment estimate is a
  sampling estimate, though: with `scaled=1000` a 4 Mb genome keeps about
  4000 k-mers, which gives a standard error of a few hundredths of a
  percent of ANI near 99%, enough to move a borderline pair across a 99.9%
  cutoff. skani samples more densely and is steadier there.
- Accessory content lowers the containment estimate at every identity,
  since plasmids and acquired regions add k-mers the other genome lacks,
  while skani reports identity within the mapped blocks and the rest as a
  lower aligned fraction. Below about 95% ANI the sourmash estimate also
  loses precision on its own terms: few k-mers survive (0.95^31 is about
  0.20, 0.90^31 about 0.04), so the variance grows and the mutation model
  fits less well.
- Completeness affects them differently. Containment is asymmetric (a
  fragment is contained in a complete genome, not the reverse), and either
  way of combining the two directions underestimates the ANI of a
  half-complete assembly to a complete genome. skani's identity is largely
  unaffected by completeness, which shows up in the aligned fraction
  instead.

On 24 Wolbachia genomes (276 pairs, two supergroups that skani does not
align against each other) the sourmash estimate (scaled=1000) ran above
skani by the following mean amounts, in percentage points, by skani band; no
pair fell between 90 and 95 %. The last row gives the range of the sourmash
estimate, without the check described below, for the pairs skani does not
report.

| skani band | Pairs | k=31 | k=21 |
|---|---|---|---|
| 95 % and above | 70 | +0.77 | +0.47 |
| 85 to 90 % | 70 | +4.46 | +2.86 |
| 80 to 85 % | 65 | +4.93 | +3.90 |
| below 80 % | 14 | +5.44 | +4.68 |
| not aligned by skani | 57 | 79.9 to 87.3 % | 80.8 to 83.7 % |

k=21 is less biased than k=31 at every band, which is why a sourmash tree of
a set wider than a species is better built with `--aligner-arg ksize=21`
(`phylo-advise` says so, and the sourmash tree builder logs it when the
genome names span more than one genus at k=31; the default stays 31).

A handful of shared hashes, from genes conserved across the set, still
gives an estimate of 0.80 to 0.88 at k=31 for pairs that have little else
in common. Wherever RepGenR2 converts a dense `sourmash compare` matrix (the
dereplicator's dense back-end, `glance`, the medoid reference,
`phylo-advise` and the sourmash tree builder with `metric=ani`), it reads
the signature sizes from `sourmash sig describe` and counts each pair's
shared hashes, j(a + b)/(1 + j) for sketches of a and b hashes. A pair with
fewer than 30 shared hashes and fewer than a tenth of the smaller sketch
has an unknown ANI. Below 30 the Poisson sampling error alone is near 20 %
of the count, and in a bacterial sketch (1000 hashes or more, where a tenth
is 100 or more, so the 30 governs) the shared hashes come mostly from
conserved genes and mobile elements rather than from the whole genome. The
relative part keeps small sketches usable: two identical sketches of a
7.5 kb viral genome share all of their about 8 hashes and keep an ANI of
1.0, while 3 shared hashes of 100 are unknown. Such a pair never merges in dereplication, is an empty cell in the
glance table and distance 1.0 in the glance dendrogram and in a sourmash
tree, counts as 0 for the medoid, and counts as below 80 % in
`phylo-advise`. On the Wolbachia set:

| k | Unknown pairs | Of the 57 skani does not align | Aligned pairs made unknown (skani) | Lowest known sketch ANI |
|---|---|---|---|---|
| 31 | 136 of 276 | 57 | 79 (79.4 to 83.2 %) | 89.85 % |
| 21 | 59 of 276 | 48 | 11 (79.4 to 80.6 %) | 83.38 % |

With sketches of about 1000 to 1600 hashes (genomes of 1.0 to 1.5 Mb at
scaled=1000), the sketch's reach therefore ends near 88 to 90 % sketch ANI
at k=31 (83 to 85 % by skani) and near 83 to 85 % at k=21 (about 80 to 81 %
by skani). Larger genomes reach lower: 30 shared hashes of 5000 (a 5 Mb
genome) is a sketch ANI of about 85 % at k=31 and 78 % at k=21. At k=21, 9 of the
unaligned pairs share 30 to 33 hashes and keep estimates of 83.4 to 83.7 %.
The sparse branchwater back-end keeps only pairs above containment 0.99^k,
hundreds of shared hashes for any bacterial genome, so the rule does not
apply there.

For the same cutoff on well-separated data the tools usually form the same
clusters (see the `clonal_50_clustered` result above); the differences
appear in fragmented assemblies, in genomes of unequal size, and at the
coarse end of the range. `glance` plots the sourmash estimate
(`--tool sourmash`) or dRep's Mash ANI (`--tool drep`), so a set's pairwise
values can be seen on the scale a cutoff will be applied to before choosing
one.

## Phylogeny routes

`phylo` combines three independent choices: the genome set (representatives or
all genomes), the source of the alignment
(`--msa-source aligner|snptype|coregene`) and the tree builder
(`--treebuilder`). The alignment source is skipped for alignment-free tree
builders.

| Route | Tools | Use it for | Limits and costs |
|---|---|---|---|
| Alignment-free | `--treebuilder mashtree`, `sourmash`, `skani` or `ska2` | Diverse sets, large sets, a quick tree; `skani` for a genus or wider, `ska2` for a clonal set. | `mashtree` and `sourmash` distances compress near identity, so clonal sets resolve poorly; `ska2` undercounts past a few per cent divergence. Neighbour joining only. No bootstrap supports. |
| Whole-genome MSA | `--aligner progressivemauve`, `sibeliaz` or `cactus`, then an MSA tree builder | Small sets of related genomes when a whole-genome alignment is wanted. | Declared limits 500, 2000 and 2000. progressiveMauve about 10 to 22 minutes per build at 20 to 50 genomes. SibeliaZ timed out at 20 genomes in the audit. Cactus needs its container and is for same-species sets. |
| Core-gene | `--msa-source coregene` after `annotate` (bakta or prokka), with `--coregene panaroo`, `ppanggolin`, `scarap`, `orthofinder` or `markers` (GTDB marker proteins) | A set where a supported, model-based tree is wanted and the reference biases a SNP tree. How much of the core survives depends on the builder's clustering setting and on how wide the set is (see [Core-gene builders](#core-gene-builders)). | Orthologous genes present in `--core-threshold` of the genomes, each aligned with mafft and concatenated. Invariant sites kept, no ascertainment correction. Slow at hundreds of genomes (not measured here); the annotation is cached per genome. Within a species a SNP tree draws on far more positions. |
| Core-SNP | `--msa-source snptype` with `--snptyper simple`, `snippy`, `parsnp` or `ska2` | One species or a clonal set. | `simple`, `snippy` and `parsnp` map every genome to one reference, so the reference biases the result. `simple` maps with minimap2 `-x asm20` and records reference positions a genome does not cover as N, so a divergent outgroup is largely N (its coverage is in the log). `ska2` is reference-free. |

### Reference bias

The aligners and the mapping typers project every genome onto one reference.
The choice of reference shapes mostly the within-cluster branching order:
swapping the reference changed 72 of 94 splits on a 50-genome set while the
deep structure held (`docs/audit/scaling-audit.md`).

How the reference is chosen:

- Unless `--reference` names a genome, the stage uses the medoid: the genome
  with the largest summed sketch ANI to the others. It is computed with
  sourmash from the sketches in the working directory (`--sketch` on the
  genome stages, or `repgenr sketch`). Ties go to the first filename, so the
  choice is the same on every run. Only existing sketches are read; the
  choice never sketches.
- The cached alignment is tied to the chosen genome, so it is rebuilt when
  the medoid changes.
- Under Nextflow, sketches reach `phylo-build` on the bacterial path with
  `params.sketch` (`PHYLO`, or `PHYLO_MSA` under `--phylo_split_msa`). The
  viral and reads paths have no sketches and use the first genome.
- The medoid is a less arbitrary choice than the first filename and is
  expected to reduce the bias, not to remove it: in a set with an
  over-represented genotype it sits inside that genotype.
- Without sketches or sourmash, or with `--reference first`, the
  alphabetically first genome is used, with a warning. On the F. tularensis
  set that first genome of the whole genus was another species, and 5.1% of
  its columns varied inside a species usually described as clonal. For a
  species run, name a genome of that species.
- The chosen genome and the rule are recorded as `reference_genome` and
  `reference_rule`, and in `reference.txt` beside the alignment (genome,
  rule and, for the medoid, its mean ANI).

### Masking

`--mask gubbins` needs `--msa-source snptype` (in `run`, `--with-snptype` also
accepts `--mask`) and a typer that writes a whole-genome alignment, which
`simple`, `snippy` and `parsnp` do. It replaces the core-SNP alignment with
Gubbins' filtered polymorphic sites. `ska2` writes variable sites only, so it
cannot be masked. Gubbins is for isolates of one species. The masker
estimates the variable fraction of the alignment, warns above 10%, and
repeats the figure if Gubbins fails. Gubbins needs a multi-threaded RAxML
build when given more than one thread; see
[usage.md](usage.md#snp-typing-and-masking).

## Core-gene builders

The core-gene route first annotates each genome (`annotate`, with `bakta` or
`prokka`). A core-gene builder then groups the genes into orthologues and
aligns each of those present in at least `--core-threshold` of the genomes.
The four clustering builders differ in how they group genes, and that
grouping, together with how wide the set is, decides how many core genes a
set retains. The fifth, `markers`, takes a fixed set of 120 (or 53) marker
genes and does not group genes at all. Panaroo's
graph cleaning (`--clean-mode`) and entropy filter remove annotation errors
and badly aligned genes. Both options are Panaroo's and have no effect with
the other builders.

| Builder | How it groups genes | Setting that governs the core | Measured on Wolbachia (`verification.md`) | Declared limit |
|---|---|---|---|---|
| `panaroo` (default) | Clusters genes in a first pass at a fixed identity, then cleans the pangenome graph. | `--aligner-arg threshold=`, the clustering identity, 0.98 by default, and 0.95 when it is not set and the canonical genome names span more than one species. A wider set needs a lower value; lower it in steps. | Three genomes: the core went from 52 genes at 0.98 to 638 at 0.90. Four genomes across supergroups shared none at the default. | 500, not measured |
| `ppanggolin` | Clusters proteins with MMseqs2 at 80% identity over 80% of their length, then partitions the pangenome graph. | `--aligner-arg identity=` and `coverage=`; `partition=` selects the soft core at `--core-threshold` (default), the strict core or the persistent partition. | 24 genomes: 398 single-copy families, where Panaroo's defaults gave none. | 5000, not measured |
| `scarap` | Builds orthogroups on seed genomes, then searches every genome against their MMseqs2 profiles. | No fixed identity cut; `method=` and `seeds=` choose the inference. | 24 genomes: 577 genes; its tree shared 20 of 21 splits with the PPanGGOLiN tree. | 10000, not measured |
| `orthofinder` | Compares all proteomes against all with DIAMOND and clusters the hits with MCL. | `search=` and `inflation=`. | 24 genomes: 585 genes; its tree shared 20 of 21 splits with the PPanGGOLiN tree and 19 with the SCARAP tree. | 200, not measured |
| `markers` | Does not cluster: searches each genome for GTDB's single-copy marker families (bac120 or ar53) with HMMER and aligns each marker's proteins to its model. Protein only. | `--marker-set` and `--core-threshold`; `min_markers=` sets the genome guard. | 24 genomes: 105 to 111 of the 120 bac120 markers per genome, 99 kept at a 0.95 threshold; its tree was four splits from the PPanGGOLiN tree, all four with supports below 95. | 20000, not measured |

Where the per-gene alignments come from:

- `panaroo` and `ppanggolin` write them themselves (mafft), and the adapter
  hands them over.
- `orthofinder` and `scarap` only group the genes; the adapter aligns each
  orthogroup's proteins with mafft and back-translates them to codons.
- `markers` aligns each marker's proteins to its HMM with `hmmalign` and
  keeps the match-state columns.

The stage then ranks, caps and concatenates the genes in the same way for
all five builders.

Trimming (`--gene-trim clipkit`) runs ClipKIT on each per-gene alignment
before the genes are measured and ranked. It may help on genus and wider
sets, where partial genes leave ragged ends and some columns are poorly
aligned; those columns add variable and informative sites that carry little
signal, and the ranking reads them as signal. Within a species the
alignments are rarely ragged and trimming is expected to change little. It
is off by default. It changes the alignment, so turning it on or off
rebuilds the core-gene alignment rather than reusing it. Its effect on tree
topology has not been measured in this repository.

Which builder to use:

- `panaroo` is the default. How much of the core survives depends on
  `--aligner-arg threshold=` and on how wide the set is. Lowering
  `--core-threshold` instead recovers fewer genes and gives a softer core.
- `ppanggolin` is a candidate for a divergent set. Its declared limit for
  large sets is not measured.
- `scarap` infers the core directly through profile searches, so there is no
  identity to lower.
- `orthofinder` runs an all-against-all search whose cost grows roughly with
  the square of the genome count, so it suits small and medium sets.
- `markers` suits a set that spans a family or more. Its search is linear
  in the genome count. It gives at most 120 genes, so within a genus the
  clustering builders draw on more of the genome.

The settings, measurements and caveats of each builder are in
[usage.md](usage.md#core-gene-phylogeny). The two annotators differ in what
they need installed: bakta needs its database and prokka needs none (see
[Which environment provides each tool](#which-environment-provides-each-tool)).

## Tree builders

| Builder | Input | Scaling | Resolution | Rooting by the builder |
|---|---|---|---|---|
| `iqtree` | MSA | Declared 500. Poor beyond that. Without `--aligner-arg model=...` it runs ModelFinder first; `model=GTR+G` skips that pass. With a SNP alignment the constant-site counts are passed as `-fconst`; `asc=lewis` selects `+ASC` instead. Accepts the core-gene partitions (`--partitioned`, `-p`); with partitions and no `model=` it applies `GTR+G` (`LG+G` for protein) to every partition, and `model=MFP` runs ModelFinder per partition. | Good on both diverse and clonal sets. Use at least 1000 ultrafast bootstrap replicates for supports, or `--aligner-arg alrt=N` (N at least 1000) for SH-aLRT supports; with neither the tree has none, and the stage says so. Above 200 genomes without `model=` the log notes that `model=GTR+G` skips ModelFinder. | Roots on the outgroup. |
| `raxmlng` | MSA | Declared 1000. Slow beyond that. | Supports from its own bootstrap procedure. With a SNP alignment the constant-site counts go on the model as `+ASC_STAM`; `asc=lewis` selects `+ASC_LEWIS` instead. Accepts the core-gene partitions (`--partitioned`, `--model partitions.txt`). | Roots on the outgroup. |
| `fasttree` | MSA | Declared 5000. Good to about 5000. `--aligner-arg gamma=true` rescales the branch lengths under the Gamma20 likelihood. | Approximate. Weaker resolution on clonal sets. | Unrooted. |
| `veryfasttree` | MSA | Declared 20000, not measured here; the same algorithm as `fasttree` with multi-core and vectorised arithmetic. `auto` picks it above 5000 leaves when an alignment is requested and never below, since `fasttree`'s tighter limit wins there; without an alignment request `auto` stays with the alignment-free builders once the set exceeds the aligner's limit. `-threads` follows `-t`; `--aligner-arg gamma=true` and `--aligner-arg double_precision=true` are extras. | As `fasttree`. | Unrooted. |
| `mashtree` | Genomes | Declared 10000. 53 min and 7.6 GB measured at 5000. | Poor on clonal sets. No supports. | Unrooted. |
| `sourmash` | Genomes | Declared 2000. The distance is one minus ANI (`metric=ani`, the default), or one minus Jaccard with `metric=jaccard`. The ANI is the k-mer model estimate c^(1/k) from the average containment c (`sourmash compare --avg-containment`), as in dereplication; `--ani` itself is not used because it reports 0 for sketches of small genomes, identical ones included. Pairs with fewer than 30 shared hashes and under a tenth of the smaller sketch (see "How the tools estimate ANI") get distance 1 and are counted in the log. For a set wider than a species, pass `--aligner-arg ksize=21`: k=21 is less biased than the default k=31 at 80 to 95 % identity. The joining step runs rapidNJ when installed (`nj=auto`), else a cubic Python step refused above 5000. | Poor on clonal sets. No supports. | Unrooted. |
| `skani` | Genomes | Declared 10000, not measured here. ANI from sparse chaining of k-mer sketches, reliable to about 80% identity; pairs below the screen (`screen=80`) get distance 1 and are counted in the log. `preset=fast|medium|slow|small-genomes`, `min_af`, `nj` as extras. Also provides the distance matrix of the viral outgroup selection, with the `small-genomes` preset by default there; when no candidate is within skani's reach no matrix is written and the run proceeds without an outgroup. | Expected to resolve closer pairs than mash distance and more distant ones than Jaccard (skani, Shaw and Yu 2023), not measured here; no supports, and a clonal set resolves only as well as ANI does. | Unrooted. |
| `ska2` | Genomes | Declared 5000, as the typer. `ska distance` SNP counts (`metric=snps`) or k-mer mismatch fractions (`metric=mismatch`) over the split k-mers each pair shares, joined with rapidNJ when installed. | Good on clonal and outbreak sets, where it is reference-free; undercounts distances once pairs share few k-mers, so it warns on a set spanning several species. No supports. | Unrooted. |
| `aster` | MSA with gene partitions | Declared 5000, not measured here. Needs `--msa-source coregene --partitioned`: one gene tree per partition (`iqtree --fast`, or `--aligner-arg gene_tree_builder=fasttree`), then a species tree with ASTER (`program=wastral`, the default, `astral4` or `astral-pro3`). Never chosen by `auto`. `rounds=`, `min_gene_columns=` and `model=` (the gene-tree model) are extras. | A summary of the gene trees that allows for discordance among them; node labels are local posterior probabilities in [0, 1], so no bootstrap is needed. `wastral` branch lengths are in coalescent units with no terminal lengths, so `--root midpoint` is refused; `astral4` writes substitutions per site. | Unrooted. |

### Automatic choice of tree builder

`--treebuilder auto` picks the builder whose declared scale fits the leaf
count, the outgroup included. A request that only an alignment can honour
(`--msa-source snptype` or `coregene`, or `--mask`) limits the choice to
builders that consume one. The full rule is in
[How `auto` chooses](#how-auto-chooses).

### Rooting

`phylo` roots `tree/tree.nwk` itself, whatever the builder: on the outgroup
by default when one is staged, at the midpoint with `--root midpoint`, or not
at all with `--root none`. The choice therefore does not depend on whether a
builder can root. Details are in [Rooting the tree](usage.md#rooting-the-tree).

### Collapsing weak splits in tree2tax

The tree builder decides which collapsing options apply. `--collapse-support`
needs supports: IQ-TREE writes them only with bootstrap replicates (`-B`) or
SH-aLRT replicates (`--aligner-arg alrt=N`),
FastTree writes SH-like local supports and RAxML-NG bootstraps by default,
while `mashtree`, the sourmash builder, `skani` and `ska2` write none.
`--collapse-length` and `--collapse-diameter` use branch lengths and so work
on every builder; the diameter maps a distance cutoff (for example about 95 %
ANI) onto the taxonomy without supports. The units are those of the tree:
mash distance for `mashtree`, one minus ANI as a fraction for `sourmash`
(one minus the Jaccard similarity with `metric=jaccard`), one minus ANI as a fraction for `skani`, the SNP count
(`metric=snps`, the default) or the fraction of mismatching k-mers
(`metric=mismatch`) for `ska2`, and substitutions per site for the ML
builders. Units and examples are in
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
[usage.md](usage.md#nextflow). Installation for each environment is in
[install.md](install.md).

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
[install.md](install.md#notes)) or name another tool.

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

## Which environment provides each tool

A tool you choose must be installed. The conda files under `envs/` split the
tools over a core environment and twelve satellites, because they do not
solve as one environment (see
[install.md](install.md#1-core-and-satellite-environments)). Append a
satellite's `bin` directory after core on `PATH`, never before it, or give it
to its tool alone with `--bin-dir TOOL=DIR` (see
[usage.md](usage.md#per-tool-directories---bin-dir)).

| Environment | Tools |
|---|---|
| `repgenr` (`envs/core.yml`) | skder, galah, sourmash, drep, sibeliaz, simple, ska2, skani, iqtree, fasttree, veryfasttree, raxmlng, aster, rapidnj, skesa, shovill, flye, medaka, racon, and `datasets` for the downloads |
| `repgenr-gubbins` | gubbins |
| `repgenr-mashtree` | mashtree |
| `repgenr-snippy` | snippy |
| `repgenr-parsnp` | parsnp |
| `repgenr-checkm2` | checkm2 (the reads chain quality step) |
| `repgenr-mauve` | progressivemauve |
| `repgenr-bakta` | bakta (needs its database) |
| `repgenr-prokka` | prokka (no database to download) |
| `repgenr-panaroo` | panaroo, with mafft |
| `repgenr-ppanggolin` | ppanggolin, with MMseqs2 |
| `repgenr-scarap` | scarap, with OrthoFinder for its `O-B` and `O-D` methods |
| `repgenr-orthofinder` | orthofinder, with DIAMOND |
| none (image only) | cactus |
