# Phylo methods audit (2026-10)

Status: audit complete (2026-10-21); this record is #43. The documentation
fixes listed under "Other defects" were made with it. The improvement items A
to E and H under "What this audit changes" are implemented (#42 to #53; the
pull request of each item is named there). F has shipped since (the marker
database and the `markers` core-gene builder), and G is pending cluster
time. The findings below describe the code at the audited commit.

This record reviews the three MSA sources and the tree builders of the `phylo`
stage at commit e5b66e5, with one question in view: which route and which
builder give a defensible tree for 100 to 500 genomes. Line numbers refer to
that commit. Measurements are cited from `docs/verification.md` and
`docs/audit/scaling-audit.md`; none were repeated for this record. A figure
called declared is a `recommended_max_genomes` value from an adapter, not a
measurement. The earlier gap analysis of the same stage is in
[phylo-tree2tax-gap-2026-10.md](phylo-tree2tax-gap-2026-10.md).

## MSA sources

| Source | What it is today | Fit for 100 to 500 genomes | Defects and gaps |
|---|---|---|---|
| `aligner` (progressiveMauve default, SibeliaZ, cactus) | Reference-anchored whole-genome MSA. | Not usable. progressiveMauve took about 22 min per build at 50 genomes (declared 500); SibeliaZ timed out at 20; cactus needs its container and same-species sets. | The default MSA source is the route least able to reach the target scale. Keep it for small sets (under about 30 genomes); it should not be the default that `phylo` falls into. |
| `snptype` (simple, snippy, parsnp, ska2; Gubbins masker) | Core-SNP alignment. The mapping typers also write the whole-genome alignment, so constant sites are counted (`-fconst` for IQ-TREE, `+ASC_STAM` for RAxML-NG). Medoid reference, `--core-fraction`, Gubbins with refusals above 12 % variable columns or more than one species. | Good within a species. The simple typer took 2.8 min on 68 genomes. parsnp (declared 2000) is a core-genome aligner built for this scale and is under-advertised. | ska2 writes no constant-site counts, so its ML branch lengths are inflated and `phylo` only warns. `--core-fraction` is unset by default, so every column where two different bases occur is kept however few genomes cover it. parsnp is listed as a typer only, although it is the route's natural within-species whole-genome aligner for hundreds of genomes. |
| `coregene` (panaroo default, ppanggolin, scarap, orthofinder) | Concatenated DNA alignment of single-copy genes present in `--core-threshold` (0.95) of the genomes. Per-genome annotation cached. | The only route that reaches a genus overview with supports. Measured only on 3 and 24 Wolbachia genomes. | See the five gaps below. |

### Whole-genome aligners (`--msa-source aligner`)

- progressiveMauve: about 10 min at 20 genomes and about 22 min at 50, per
  build (scaling-audit.md, B4 and B5 runs, Docker on Rosetta). Declared 500.
  Not run above 50 genomes.
- SibeliaZ: every cell from 20 genomes up timed out after one hour (two hours
  at 100) without a tree (scaling-audit.md, timing table). Declared 2000, which
  no run supports.
- cactus: needs its container (Rosetta on Apple Silicon) and is meant for
  same-species sets. Declared 2000, not measured.
- The defaults are `msa_source = "aligner"` and `aligner = "progressivemauve"`
  (`stages/phylo.py:176-177`). A user who runs `phylo` on 300 genomes without
  naming a route is therefore sent to the route that has not been shown to
  finish there. With `--treebuilder auto` and no alignment request, the stage
  avoids the aligner past its declared limit (`stages/phylo.py:642-654`); with
  the default `--treebuilder iqtree` it does not.

### Core-SNP typing (`--msa-source snptype`)

- The mapping typers (simple, snippy, parsnp) write the whole-genome alignment
  beside the core-SNP alignment, and the stage passes the constant-site counts
  to IQ-TREE as `-fconst` (`treebuilders/iqtree.py:103`) and to RAxML-NG as
  `+ASC_STAM` (`treebuilders/raxmlng.py:64`).
- ska2 writes no whole-genome alignment. `phylo` logs that the branch lengths
  carry no ascertainment-bias correction and are inflated
  (`stages/phylo.py:798-810`), and builds the tree anyway.
- `--core-fraction` is unset by default (`stages/phylo.py:95-108`). Without it
  every column where at least two different bases occur is kept (the simple
  typer, `snptypers/simple.py:386`), including columns covered by two genomes
  out of hundreds.
- Gubbins warns above 10 % variable columns and refuses above 12 %
  (`maskers/gubbins.py:182-186`); it also refuses a set whose names span more
  than one species (`maskers/gubbins.py:217`). On Francisella subsets it ran at
  8.3 % and failed from 13.2 % upward (verification.md).
- The simple typer took 2.8 min on the 68 Francisella representatives with
  concurrent calling and compressed BCF, 10.8 min before that change
  (verification.md, SNP typing throughput).
- parsnp has passed the live suite and the gzip tests on small sets only.

### Core genes (`--msa-source coregene`)

Measured on 3 and 24 Wolbachia genomes only (verification.md, core-gene route).
On 24 genomes PPanGGOLiN kept 398 single-copy soft-core families (419047
columns), SCARAP 577 genes, OrthoFinder 585 orthogroups; Panaroo at its defaults
found no gene at 95 % occupancy.

- (a) No gene cap or ranking. Every gene that passes `--core-threshold` enters
  the concatenation. SCARAP alone has an occupancy-sorted `-m`
  (`coregenes/scarap.py:124,153`).
- (b) The per-gene alignments are removed with the scratch directory
  (`stages/coregene.py:203`), so there is no partition file and no per-gene
  quality control.
- (c) DNA only. A genus or family span saturates third codon positions, and no
  protein or codon alphabet is offered; IQ-TREE is always called with `-st DNA`
  (`treebuilders/iqtree.py:87-88`). Shipped (A5): `--coregene-alphabet
  protein` gives a protein core-gene alignment with ppanggolin, scarap or
  orthofinder, and the tree builders use `-st AA`, `LG+G` or `-lg`; codon
  models remain deferred.
- (d) Panaroo's default clustering identity (`-c 0.98`) gives an empty core on
  a divergent set: 0 genes at 95 % occupancy on 24 Wolbachia genomes, 52 genes
  in all three genomes of one supergroup (verification.md). The default
  builder is therefore the one least suited to the genus purpose of the route.
- (e) Invariant sites are kept, so no `-fconst` is needed, but
  `_warn_low_diversity` (`stages/phylo.py:1647`) is the only check on the
  alignment.
- Annotation costs about 45 s per genome with bakta's light database; prokka
  "took 135 seconds against bakta's 45 per genome" on the same three genomes
  (verification.md). The annotation is cached per genome.

## Tree builders

| Builder | Input | Strengths | Weaknesses | Use when |
|---|---|---|---|---|
| iqtree | MSA | ML with ModelFinder, UFBoot, `-fconst`, outgroup; the best resolution of the set; measured in the live suite. | Declared 500, not measured above that. ModelFinder on a 500-taxon, 500 kb concatenation is slow; no partition models; always `-st DNA`. | Default for any alignment up to a few hundred taxa; `model=GTR+G` skips ModelFinder. |
| raxmlng | MSA | Thorough bootstrap (`--all`), `+ASC_STAM`. | Always bootstraps (`autoMRE{200}` when `-B 0`, `treebuilders/raxmlng.py:81-89`); expected to be several times slower than IQ-TREE at 500 taxa (not measured here); declared 1000. | When RAxML-style bootstraps are wanted; otherwise IQ-TREE. |
| fasttree | MSA | Minutes at thousands of taxa; SH-like supports. | Approximate; no ascertainment correction; weaker on clonal sets. | Above the IQ-TREE limit, or as Gubbins' internal builder. |
| veryfasttree | MSA | FastTree algorithm, multi-core; declared 20000. | Same approximations; not measured here. | Above 5000 leaves. |
| mashtree | genomes | Seconds to minutes; 53 min and 7.6 GB measured at 5000. | Mash distance saturates near identity; NJ only; no supports; branch lengths are k-mer distances. | A quick overview of a diverse set; viral sets. |
| sourmash | genomes | Reuses workdir sketches; rapidNJ when installed. | The distance is `1 - Jaccard`, not ANI. It is non-linear in ANI: about 0.88 at 95 % ANI with k = 31 (usage.md), and above 0.99 below about 85 % ANI by the Mash relation between Jaccard and ANI, so a genus tree is compressed at its base. No supports. | Only when sketches exist and the set is within a genus; skani is preferable. |
| skani | genomes | ANI from sparse chaining, reliable to about 80 % identity (Shaw and Yu 2023); declared 10000. | NJ only; no supports; not measured here. | A genus or family overview without supports; a first look before a core-gene run. |
| ska2 (builder) | genomes | Reference-free SNP distances; suited to clonal sets. | Undercounts past a few per cent divergence; no supports; `auto` never picks it (`treebuilders/ska2.py:67`). | A quick tree of an outbreak or clonal set. |

Across the builders: every alignment-free builder is a distance matrix joined by
neighbour joining. None applies minimum-evolution refinement (FastME balanced
minimum evolution with NNI or SPR moves), and none writes supports. The ML
builders have no partition support.

## Other defects

Each item was checked against commit e5b66e5. Items 1 to 3 and 6 to 7 were
corrected with this record. Items 4 and 5 are code defects and are listed under
item B of the ladder below; both have shipped.

1. `docs/choosing-tools.md:131-136` said that `auto` picks raxmlng for 501 to
   1000 genomes and fasttree up to 5000 without an alignment request. The code
   (`stages/phylo.py:642-654`) restricts `auto` to the alignment-free builders
   once the leaf count exceeds the aligner's limit (500 for progressiveMauve),
   so it picks sourmash (limit 2000) and then mashtree (limit 10000; the tie
   with skani is broken alphabetically, `core/plugins.py:347-353`). The table
   predated that change (#15). Rule 5 of "How `auto` chooses" on the same
   page also left out `--msa-source coregene`, which `_wants_alignment`
   (`stages/phylo.py:610-618`) counts as an alignment request.
2. The units of `--collapse-length` and `--collapse-diameter` in
   `docs/usage.md:2492-2515` and `docs/choosing-tools.md:493-496` omitted skani
   (`1 - ANI` as a fraction) and ska2 (SNP count or mismatch fraction).
3. `treebuilders/skani.py:13-14` said the command surface was checked against
   skani 0.4.0, while the image pin is `skani:0.3.2` (`treebuilders/skani.py:42`).
   `tree/nj.py:16-17` said rapidNJ 2.3.2, while the pin is 2.3.3
   (`tree/nj.py:43`).
4. `stages/phylo.py:1332-1337` reads the reference note from `align/` for every
   source other than `snptype`. For a core-gene run a `reference.txt` left by an
   earlier aligner run can therefore be recorded as the `reference_genome` of a
   core-gene tree. The core-gene path clears the note beside its own alignment
   (`stages/phylo.py:532`, `tree/msa/reference.txt`), not `align/reference.txt`.
   Shipped: the note is read from the directory the current source writes
   to (`align/` for the aligner, `tree/msa/` for a SNP typer), and a route
   without a reference (coregene, a reference-free typer) records None.
5. `_msa_key` (`stages/phylo.py:339`) falls back to the `--reference` string
   when no reference genome was resolved, which is the case for coregene and
   ska2. Changing a reference those routes ignore rebuilds the alignment.
   Shipped: the key holds the reference only for the aligner and for a SNP
   typer that maps against one.
6. `docs/audit/phylo-tree2tax-gap-2026-10.md:3-4` still said batch E was
   deferred; E2 to E4 shipped as #15, #16 and #18 (with the skani builder, #17).
7. `docs/verification.md` gave the PPanGGOLiN stage time on 24 Wolbachia
   genomes as 221 s in the end-to-end bullet and 186 s in the summary table,
   without saying why. They are two live runs of the same command: 221 s the
   run in the prose, 186 s the re-run after the review fixes of #22 (commit
   3953525), with the same genes and columns. Both figures stand, and the
   table caption now says so.

## Recommendations

### By purpose, span and ANI

The items named in this table (A5, C, F) have shipped.

| Purpose | Span, ANI | Route for 100 to 500 genomes | Tree builder | Notes |
|---|---|---|---|---|
| Outbreak or clonal genotype tree | one lineage, ANI above 99.5 % | `--msa-source snptype --snptyper ska2` (reference-free) | iqtree with `--aligner-arg model=GTR+G` and `-B 1000` | ska2 writes constant-site counts (item C, shipped), except under `lo=true`. The ska2 tree builder gives a quick look. |
| Genotype tree within a species | one species, ANI 95 to 99.5 % | `--msa-source snptype --snptyper simple` (medoid reference) or `parsnp`; `--mask gubbins` when recombination is expected and under 12 % of columns vary; `--core-fraction 0.95` | iqtree with `-fconst` (automatic) and `-B 1000`; fasttree as Gubbins' internal builder | ska2 with `--aligner-arg lo=true` is the reference-free alternative when no masking is wanted. |
| Overview of species representatives in a genus | several species, ANI 80 to 95 % | `--msa-source coregene --coregene ppanggolin` (or `scarap`), `--core-threshold 0.95` | iqtree with `model=GTR+G` and `-B 1000`; skani NJ first for a quick check | Gene cap and ranking, partitions and a protein alphabet (`--coregene-alphabet protein`, A5, shipped) are available. Panaroo only with `--aligner-arg threshold=` lowered to about 0.90; its default core is empty here. |
| Overview across a family or wider | ANI 75 to 85 % or below | GTDB marker proteins (`--coregene markers`, item F, shipped; `--preset family`) or protein core genes (`--coregene-alphabet protein`, A5, shipped) | iqtree `-st AA` (ModelFinder, or `model=LG+G`) | DNA alignments saturate; mash and Jaccard distances are unreliable; skani reaches about 80 % and gives an overview without supports. |
| Quick look at any set | any | alignment-free | skani (genus or wider), mashtree (diverse or viral), ska2 (clonal) | No supports; use `--collapse-diameter` in tree2tax rather than `--collapse-support`. |

Rules of thumb. A SNP route needs a genome every other genome maps to, so it
stops being useful past about 5 % divergence. A core-gene route needs
annotation (about 45 s per genome with bakta light) and works until the core
collapses, which the builder's clustering setting governs. A whole-genome MSA
is for sets under about 30 genomes.

### Builder by leaf count, with an alignment

- Under 500: iqtree.
- 500 to 1000: iqtree with `model=GTR+G` (the declared limit of 500 is not
  measured; item G settles it), else raxmlng.
- 1000 to 5000: fasttree.
- Above 5000: veryfasttree.

## What this audit changes

One pull request per item, or per part of an item. A to F and H are
implemented; G is pending cluster time.

- **A. Core-gene selection and layout.** A gene cap with ranking (occupancy,
  informative sites, variable sites, length, random), per-gene alignments kept
  and reported, a partition file for IQ-TREE and RAxML-NG, a protein alphabet,
  and optional trimming. Shipped: A1 to A4, that is per-gene alignments,
  the gene selection and `partitions.json` (A1, #42), per-gene alignments
  from all four builders (A2, #44), the selection options (A3, #45) and
  `--partitioned` (A4, #47). The protein alphabet (A5, #49) has shipped as
  `--coregene-alphabet protein`. Optional trimming (A6) has shipped as
  `--gene-trim clipkit` (#51): ClipKIT in smart-gap mode trims each per-gene
  alignment before the genes are measured, ranked and concatenated; off by
  default, and other ClipKIT modes are not offered on the command line.
- **B. Default route and builder hygiene.** A lower Panaroo identity on
  multi-species sets, a recommended IQ-TREE recipe for large alignments,
  `1 - ANI` as the sourmash builder's distance, and defects 4 and 5 above.
  Shipped in #52: B1, Panaroo's `-c` is 0.95 instead of 0.98 when `threshold=` is
  not set and the canonical genome names span more than one species; B2,
  IQ-TREE's `alrt=N` extra (`-alrt`, SH-aLRT supports without `-B`) and a
  log line above 200 genomes without `model=` that `model=GTR+G` skips
  ModelFinder; B3, the sourmash builder's distance is `1 - ANI` by default
  (`metric=ani`; `metric=jaccard` keeps `1 - Jaccard`). The ANI is the k-mer
  model estimate from the Jaccard index, `(2j / (1 + j))^(1/k)`, shared with
  the dereplicator, rather than `sourmash compare --ani`, which reports 0 for identical genomes of 10 kb
  and 30 kb at k=31, scaled=1000 (sourmash 4.9.4). Defects 4 and 5 are
  fixed.
- **C. ska2 constant sites.** Counting constant columns from an unfiltered ska2
  alignment so that the reference-free route gets `-fconst`. Shipped (#46): the
  typer writes an all-sites alignment with `ska align --filter no-filter`, the
  snptype stage counts it into `alignment_stats.json` and removes it; not
  under `lo=true`, where `ska lo` writes variants only.
- **D. Advice and presets.** A `phylo-advise` command that reads a workdir and
  prints the matching row of the table above, and a `--preset` option that
  expands to it.
  Shipped: `phylo-advise` (D1, #48), with the table held as data in
  `src/repgenr/stages/advise.py`, and `--preset` on `phylo`, `phylo-build`
  and `run` (D2, #53), which expands the same rows; an option given on the command
  line overrides the preset. The genus row carries `--partitioned` and the
  family row adds `--coregene-alphabet protein`; per-gene trimming stays an
  explicit choice. `phylo-advise` prints the command with `--preset` and its
  expansion.
- **E. Documentation.** This record, the corrections listed under "Other
  defects", and the recommendation tables in
  [choosing-tools.md](../choosing-tools.md#by-purpose-and-divergence).
  Shipped with this record (#43).
- **H. Species tree from gene trees.** A tree builder `aster` (ASTER,
  github.com/chaoszhang/ASTER, the maintained successor of ASTRAL) that reads
  the partition file of item A, builds one gene tree per partition, and
  summarises them into a species tree with weighted ASTRAL (`wastral`) by
  default. Intended for genus and family overviews, where discordance among
  gene trees makes a concatenated ML tree overconfident. Depends on item A.
  Shipped (#50): the `aster` tree builder (`--treebuilder aster --msa-source
  coregene --partitioned`); gene trees by IQ-TREE (`--fast`, aBayes) or
  FastTree, summarised by `wastral` (default), `astral4` or `astral-pro3`.

Shipped after the audit: a marker-gene route for family-level sets (F),
`repgenr markers-db download` (GTDB's bac120 and ar53 marker models from
InterPro) and the `markers` core-gene builder (`--coregene markers`, HMMER
search and alignment with GTDB-Tk's rules; the builder of `--preset
family`). Pending cluster
time: a benchmark of IQ-TREE and capped core-gene trees at 100 to 500 genomes
on real data (G).

## Sketch ANI audit (2026-10-21)

The ANI that RepGenR2 computes from a `sourmash compare` Jaccard matrix,
`(2j / (1 + j))^(1/k)`, is identical to sourmash's own estimate (maximum
difference 0 on the 276 pairs of 24 Wolbachia genomes at k=21 and k=31).
Against skani on the same genomes (scaled=1000; mean sourmash minus skani,
percentage points; the last row is the range of the unchecked sourmash
estimate for the 57 pairs skani does not align):

| skani band | Pairs | k=31 | k=21 |
|---|---|---|---|
| 95 % and above | 70 | +0.77 | +0.47 |
| 85 to 90 % | 70 | +4.46 | +2.86 |
| 80 to 85 % | 65 | +4.93 | +3.90 |
| below 80 % | 14 | +5.44 | +4.68 |
| not aligned by skani | 57 | 79.9 to 87.3 % | 80.8 to 83.7 % |

Two changes follow.

- **Shared-hash minimum.** The dense conversion reads the signature sizes
  from `sourmash sig describe --csv` and sets the ANI of a pair to unknown
  (NaN) when its sketches share fewer than 30 hashes and fewer than a tenth
  of the smaller sketch (`MIN_SHARED_HASHES` and `MIN_SHARED_FRACTION` in
  `core/sourmash.py`; the count is j(a + b)/(1 + j), equal to the exact
  intersection on these sketches). At these sketch sizes (893 to 1630
  hashes) the 30 governs; the relative part keeps identical sketches of
  small genomes, such as 8-hash viral sketches, at ANI 1.0. The dereplicator
  never merges such a pair, glance writes an empty cell, the medoid counts
  it as 0, the sourmash tree places it at distance 1.0, and `phylo-advise`
  takes the minimum and median over the known pairs, reports the unknown
  count (`ani_unknown_pairs`) and treats an unknown pair as below 80 %, so
  the Wolbachia set is advised as family. Before the change the same
  evidence (canonical names giving one species, `Wolbachia sp`, and a
  minimum sketch ANI of 79.88 % at k=31) gave the species preset, with a
  note to check for mislabelled genomes, since a single species by name is
  tested before the ANI; with names that are not canonical the 79.88 %
  minimum already gave family. The 57 pairs
  skani does not align share 1 to 22 hashes at k=31 and 13 to 33 at k=21;
  the lowest pair skani aligns (79.4 %) shares 8 and 26. Below 30 the
  Poisson sampling error alone is near 20 % of the count. Live
  (`repgenr phylo-build --treebuilder sourmash --no-outgroup --root none
  -t 8`, and `--aligner-arg ksize=21`):

  | k | Unknown pairs | Of the 57 skani does not align | Aligned pairs made unknown (skani) | Lowest known sketch ANI |
  |---|---|---|---|---|
  | 31 | 136 of 276 | 57 | 79 (79.4 to 83.2 %) | 89.85 % |
  | 21 | 59 of 276 | 48 | 11 (79.4 to 80.6 %) | 83.38 % |

  With sketches of about 1000 to 1600 hashes (genomes of 1.0 to 1.5 Mb at
  scaled=1000), the sketch's reach therefore ends near 88 to 90 % sketch ANI
  at k=31 (83 to 85 % by skani) and near 83 to 85 % at k=21 (about 80 to 81 %
  by skani). Larger genomes reach lower: 30 shared hashes of 5000 (a 5 Mb
  genome) is a sketch ANI of about 85 % at k=31 and 78 % at k=21.

  A first version applied sourmash's own probability check
  (`get_exp_probability_nothing_common`, threshold 1e-3). It marked 27 pairs
  at k=31 and none at k=21: it asks whether sharing no hash at all would be
  plausible, not whether a few hashes support an estimate, and was replaced
  by the minimum.
- **k=21 for genus-level sourmash trees.** `phylo-advise` adds, when the set
  spans more than one genus by name, has a minimum ANI below 95 % or has
  unknown pairs, that a sourmash tree should be built with `--aligner-arg
  ksize=21`, and the quick row's summary says so; the sourmash tree builder
  logs the same at k=31 when the genome names span more than one genus. The
  default k stays 31.
