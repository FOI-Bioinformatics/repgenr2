# Phylo methods audit (2026-10)

Status: audit complete (2026-10-21). The documentation fixes listed under
"Other defects" were made with this record. The improvement items A to D and
H under "What this audit changes" are planned and not implemented; nothing in
this record describes them as available.

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
  (`treebuilders/iqtree.py:87-88`).
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
item B of the ladder below.

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
5. `_msa_key` (`stages/phylo.py:339`) falls back to the `--reference` string
   when no reference genome was resolved, which is the case for coregene and
   ska2. Changing a reference those routes ignore rebuilds the alignment.
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

Items marked planned depend on ladder work that is not implemented.

| Purpose | Span, ANI | Route for 100 to 500 genomes | Tree builder | Notes |
|---|---|---|---|---|
| Outbreak or clonal genotype tree | one lineage, ANI above 99.5 % | `--msa-source snptype --snptyper ska2` (reference-free) | iqtree with `--aligner-arg model=GTR+G` and `-B 1000` | ska2 writes constant-site counts (item C, shipped), except under `lo=true`. The ska2 tree builder gives a quick look. |
| Genotype tree within a species | one species, ANI 95 to 99.5 % | `--msa-source snptype --snptyper simple` (medoid reference) or `parsnp`; `--mask gubbins` when recombination is expected and under 12 % of columns vary; `--core-fraction 0.95` | iqtree with `-fconst` (automatic) and `-B 1000`; fasttree as Gubbins' internal builder | ska2 with `--aligner-arg lo=true` is the reference-free alternative when no masking is wanted. |
| Overview of species representatives in a genus | several species, ANI 80 to 95 % | `--msa-source coregene --coregene ppanggolin` (or `scarap`), `--core-threshold 0.95` | iqtree with `model=GTR+G` and `-B 1000`; skani NJ first for a quick check | Gene cap and ranking, partitions and a protein alphabet are planned (item A). Panaroo only with `--aligner-arg threshold=` lowered to about 0.90; its default core is empty here. |
| Overview across a family or wider | ANI 75 to 85 % or below | protein core genes (planned, item A) or a marker-gene route (deferred, item F) | iqtree `-st AA`, `LG+G` (planned) | DNA alignments saturate; mash and Jaccard distances are unreliable; skani reaches about 80 %. Until then, skani gives an overview without supports. |
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

Planned work, one pull request per item; none of it is implemented at the time
of this record.

- **A. Core-gene selection and layout.** A gene cap with ranking (occupancy,
  informative sites, variable sites, length, random), per-gene alignments kept
  and reported, a partition file for IQ-TREE and RAxML-NG, a protein alphabet,
  and optional trimming.
- **B. Default route and builder hygiene.** A lower Panaroo identity on
  multi-species sets, a recommended IQ-TREE recipe for large alignments,
  `1 - ANI` as the sourmash builder's distance, and defects 4 and 5 above.
- **C. ska2 constant sites.** Counting constant columns from an unfiltered ska2
  alignment so that the reference-free route gets `-fconst`. Shipped: the
  typer writes an all-sites alignment with `ska align --filter no-filter`, the
  snptype stage counts it into `alignment_stats.json` and removes it; not
  under `lo=true`, where `ska lo` writes variants only.
- **D. Advice and presets.** A `phylo-advise` command that reads a workdir and
  prints the matching row of the table above, and a `--preset` option that
  expands to it.
  Shipped: `phylo-advise` (D1), with the table held as data in
  `src/repgenr/stages/advise.py`. The genus and family rows leave out the
  partition file and the protein alphabet until those options exist. The
  `--preset` option is open.
- **E. Documentation.** This record, the corrections listed under "Other
  defects", and the recommendation tables in
  [choosing-tools.md](../choosing-tools.md#by-purpose-and-divergence).
- **H. Species tree from gene trees.** A tree builder `aster` (ASTER,
  github.com/chaoszhang/ASTER, the maintained successor of ASTRAL) that reads
  the partition file of item A, builds one gene tree per partition, and
  summarises them into a species tree with weighted ASTRAL (`wastral`) by
  default. Intended for genus and family overviews, where discordance among
  gene trees makes a concatenated ML tree overconfident. Depends on item A.

Deferred: a marker-gene route for family-level sets (F) and a benchmark of
IQ-TREE and capped core-gene trees at 100 to 500 genomes on real data (G).
