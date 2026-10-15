# Choosing tools

RepGenR2 has eight tool families, and most stages work with the defaults. This
page helps when the defaults do not fit: it goes from the kind of dataset, to
its size, to the individual families. Measurements come from
[verification.md](verification.md) and
[audit/scaling-audit.md](audit/scaling-audit.md). Where a figure is a declared
limit and not a measurement, the text says so.

## 1. Decide by dataset first

| Dataset | Entry path | Dereplicator | Phylogeny route | Masking | Tree builder |
|---|---|---|---|---|---|
| Genus or wider (diverse) | `metadata` (GTDB), or `--source api` for one taxon | `skder`, or `sourmash` for large sets | Alignment-free | No | `mashtree` |
| One species | `metadata`, or `ingest` for local genomes | `skder` or `galah` | Core-SNP with `--reference` set to a genome of that species, or whole-genome MSA for small sets | `--mask gubbins` when under 10% of alignment columns vary | `iqtree` |
| Outbreak or clonal set | `ingest`, or `reads` and `assemble` | `skder` or `galah`, with the clone-block caveat in section 4 | Core-SNP with `ska2` (reference-free) or a mapping typer | `gubbins` with a mapping typer only; `ska2` cannot be masked | `iqtree` or `fasttree` |
| Viral | `vmetadata` and `vgenome` (NCBI Virus) | `skder` or `sourmash` | Alignment-free | No | `mashtree` |

What each row rests on:

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
  split. These are the only species-level measurements, on one genus.
- **Outbreak or clonal.** The dereplicators were measured on synthetic clonal
  sets in the scaling audit (all recovered the true clustering). The SNP typers
  were run on a synthetic within-species set and on the Francisella data
  (`verification.md`). There is no measured run on a real outbreak. The
  reasons to prefer `ska2` here are design reasons from the audit: it is
  reference-free, so no reference assembly biases the distances, and Mash
  distances saturate near identity, which makes `mashtree` weak on clonal sets.
- **Viral.** The live suite runs the viral chain on Hepatovirus records from
  NCBI Virus with `sourmash` and `mashtree` (`test_run_viral_chain_end_to_end`,
  21 seconds). The largest viral run recorded is 1256 Hepeviridae (hepatitis E)
  genomes from BV-BRC: `skder` gave 799 representatives, followed by `mashtree`
  and `tree2tax` (recorded in `docs/verification.md` before #150, commit 7309fe7; the scaling audit keeps the 1256 figure). Other viral scaling is not measured.

## 2. Then by size

| Genomes | What was measured | What `auto` picks |
|---|---|---|
| Under 100 | `mashtree` 9 s at 100 genomes; sourmash tree 25 s at 100. progressiveMauve about 10 minutes at 20 genomes and about 22 minutes at 50, per build. | Tree builder: `iqtree` (limit 500). |
| 100 to 1000 | `skder` about 7 min and 1.1 to 1.5 GB at 1000; `galah` 4 to 6 min; `sourmash` about 1 min and 0.2 GB. `mashtree` 117 s; sourmash tree 142 s. | Tree builder: `iqtree` up to 500, then `raxmlng` (limit 1000). |
| 1000 to 5000 | `skder` about 69 min and 6.4 to 8.2 GB at 5000 in one pass, about 15 min chunked; `galah` about 63 min; `sourmash` about 5.5 min sparse. `mashtree` 53 min and 7.6 GB at 5000. | Tree builder: `sourmash` (limit 2000), then `fasttree` (limit 5000). |
| Over 5000 | Not measured. Single-pass `skder` or `galah` at 10000 would be memory-risky on the 18 GB audit machine. Chunked `skder` or sparse `sourmash` are the safer choices. | Tree builder: `mashtree` (limit 10000). |

Notes on the table:

- Times are from one machine (11 cores, 18 GB). Treat them as orders of
  magnitude.
- The sourmash tree builder joins its distance matrix with rapidNJ when the
  `rapidnj` binary is installed (`envs/core.yml` has it) and with a
  pure-Python step otherwise. A fit to the measured times of the Python step
  gives about 1.6 hours at 5000 genomes, so without rapidNJ more than 5000
  leaves are refused; with it the joining is reported by rapidNJ's authors to
  take seconds at that size (not measured here) and the dense `sourmash
  compare` matrix becomes the limit. The declared limit stays 2000
  until that is measured. `--aligner-arg nj=python` or `nj=rapidnj` forces
  one backend.
- **Most limits are declared, not benchmarked.** The audit found the declared
  `recommended_max_genomes` values to be round numbers with no recorded
  measurement behind them (500, 1000, 2000, 5000, 10000). The scaling audit
  later measured some of them (for example the sourmash tree builder, whose
  limit was lowered after the measurement). The aligners and SNP typers have not
  been run at scale. In the audit, SibeliaZ did not finish within an hour at 20
  genomes, so its declared limit of 2000 is not supported by a run. The largest
  recorded end-to-end runs are the 1157-genome Francisella set (bacterial) and
  1256 Hepeviridae genomes (viral).
- Past a few thousand genomes, dereplicate in chunks (`--process-size` on the
  CLI, `--derep_process_size` in Nextflow) to bound the time and memory of one
  tool call. Chunking changes which genome represents a cluster, not the
  clustering (scaling audit, finding 3). The chunk size bounds the leaf
  chunks only. The final pass compares the union of the chunk
  representatives, and for a diverse set at a strict threshold that union can
  be most of the input: 6888 Bacillaceae genomes at `--secondary-ani 0.9995`
  kept 3807 after the chunks. One chunk of the next level is probed before
  that level runs, since a level recomputes the per-genome work for every
  representative and may remove few; the union goes to the final pass when
  the probe removes less than the level is expected to cost. The final pass
  then still compares about a third of the single-pass pairs, and the
  per-genome work is paid twice, so a larger `-s`, or `--pre-secondary-ani`
  below the final threshold, is the better lever (see
  [Chunked dereplication](usage.md#chunked-dereplication)).
  With skDER the adapter sets skani's screen to one point below a cutoff of
  99 or more (`--tool-arg screen=` overrides), so most same-species pairs
  are not aligned. The genome count does not limit a
  single call through the command line: every dereplicator receives the
  genome set as a list file or a staged directory, and a command that would
  still exceed the operating system's argument-size limit (ARG_MAX) stops
  before the tool starts, with a message that names `--process-size`.

## 3. Entry points

| Situation | Use |
|---|---|
| Bacteria or archaea in GTDB, many taxa or a whole release | `metadata` with the default `--source tsv`, which downloads the release table once. |
| One taxon only | `metadata --source api`, which asks the GTDB API for the target taxon and also fetches each genome's quality values. |
| Genomes already on disk | `ingest --genomes-dir DIR`. Add a `selection.tsv` for taxonomy and quality. |
| Sequencing runs, no assembly yet | `reads`, then `assemble`. Runs come from ENA, which mirrors SRA. |
| Viruses | `vmetadata` and `vgenome`. NCBI Virus is the default source. Use `--source bvbrc` only for the legacy BV-BRC path. |

`run` chains the stages for the GTDB path, and takes `--viral`, `--reads` or
`--genomes-dir` for the others (see [usage.md](usage.md)).

## 4. Dereplicators

| Tool | Scaling shape | Representative criterion | Clone-block behaviour | Needs CheckM data | Verified |
|---|---|---|---|---|---|
| `skder` (default) | Superlinear in practice: 5 times the genomes cost about 10 times the wall time and 6 times the memory. Single pass: about 7 min at 1000 and 69 min at 5000, memory up to 8.2 GB. Chunked: about 15 min at 5000 (measured before the skani screen default below). Receives the genome set as one staged directory of links, so the set size is not bounded by the command line. At a cutoff of 99 percent or more, skani's screen is set to the cutoff minus one (`--tool-arg screen=` overrides), so pairs below it are not aligned. | Its own aggregate score. Not quality-aware. | One representative per block. Which member is arbitrary. | No | Native and container (Wave) |
| `galah` | Built for large sets. About 4 to 6 min at 1000 and 63 min at 5000. | With manifest quality for every genome of the run, galah's quality score (passed as `--genome-info`). Otherwise input position, with the genomes listed by descending file size, so a fragment is not preferred over a complete genome. | Without quality, the first listed member in the three orderings tested. | No; used when present | Native and container |
| `sourmash` | Sparse back-end close to linear in close pairs. About 5.5 min at 5000, 0.3 to 0.7 GB. Dense back-end is capped at 5000. | Most-connected genome, alphabetical tie-break. | Biased toward the most-sequenced genotype. The sparse and dense back-ends can pick different members. | No | Native and container |
| `drep` | Quadratic within primary clusters. Declared limit 2000, chunk-wrapped. | Completeness, contamination, N50 and size score. Quality-aware. | Best-scored member, so least sensitive to block size. | Yes: completeness and contamination in the manifest (or `selection.tsv` for the Nextflow steps) for every genome of the run, passed as `--genomeInfo`; else CheckM on `PATH`; or `--virus`, which passes `--ignoreGenomeQuality`. Without any of these, dRep stops and `dereplicate` exits 6. | Container only |

Guidance:

- Use `skder`, the default, for most sets.
- Use `sourmash` with `--derep_process_size 2000` (Nextflow) or `--process-size
  2000` for 10000 genomes or more (not measured at this size).
- Use `drep` when every genome has manifest quality (GTDB selections) or
  CheckM is available, and keep it under about 2000 genomes. With manifest
  quality for every genome, dRep scores genomes with those values instead of
  running CheckM, the same values `--keeper quality` uses; when some genomes
  lack them, the values are not passed (dRep needs all of them) and the log
  says so. The decision covers the whole run, so a chunked run never scores
  one chunk from the manifest and another with CheckM. Genomes dRep's filter removes (`--tool-arg length=N`, or its
  completeness and contamination thresholds) are `fail_qc` in
  `genome_status.tsv` and in no cluster.
- The representative inside a clone block depends on the tool, its back-end and
  the accession names, not on genome quality. The default `--keeper quality`
  re-picks each cluster's representative by `completeness - 5 x contamination
  + 0.5 x log10(N50)` (dRep's default weights), using the values in the
  manifest and the N50 of each genome file, whichever tool clustered. Ties go
  to higher completeness, then lower contamination, then higher N50, then the
  filename, so skDER, sourmash and galah end with the same keepers when they
  form the same clusters. An unscored representative is replaced only by a
  high-quality genome (completeness above 90, contamination below 5).
  `--keeper tool` keeps the adapter's choice. `--keeper gtdb` keeps a GTDB
  species representative when the cluster holds one and applies the quality
  rule otherwise, so that the representatives coincide with GTDB's where the
  clustering allows; it needs a GTDB selection from `metadata`.
- On the synthetic set `clonal_50_clustered` (groups of 20, 15 and 15 genomes,
  within-group ANI about 0.995 or higher, between-group about 0.96 or lower),
  `sourmash`, `galah` and `skder` at the defaults formed the same three
  clusters, and each kept a different representative. To check this on your own
  data, see "Comparing two dereplications" in
  [usage.md](usage.md#inspecting-a-dereplication).
- `repgenr glance` (an all-against-all overview before choosing thresholds)
  runs on `drep` or `sourmash`; the other dereplicators have no comparison
  mode. `--tool auto`, the default, uses dRep when it can run and sourmash
  otherwise. sourmash needs no container on most systems and plots the ANI
  estimate that `dereplicate --tool sourmash` thresholds; it holds the full
  N x N matrix and refuses sets above 5000 genomes. It reads the genome
  sketches of the working directory, so it no longer sketches when they
  exist. On one workstation (8 threads, genomes on a USB disk) a sourmash
  glance of `clonal_50_clustered` took about 4 s with the sketches against
  7 s sketching, and of `balanced_1000_clustered` about 20 s against 104 s;
  `dereplicate --tool sourmash` on the 1000 genomes took about 6 s against
  65 s.

## 5. Phylogeny routes

`phylo` combines three independent choices: the genome set (representatives or
all genomes), the source of the alignment (`--msa-source aligner|snptype`) and
the tree builder (`--treebuilder`). The alignment source is skipped for
alignment-free tree builders.

| Route | Tools | Use it for | Limits and costs |
|---|---|---|---|
| Alignment-free | `--treebuilder mashtree` or `sourmash` | Diverse sets, large sets, a quick tree. | Distances compress near identity, so clonal sets resolve poorly. No bootstrap supports. |
| Whole-genome MSA | `--aligner progressivemauve`, `sibeliaz` or `cactus`, then an MSA tree builder | Small sets of related genomes when a whole-genome alignment is wanted. | Declared limits 500, 2000 and 2000. progressiveMauve about 10 to 22 minutes per build at 20 to 50 genomes. SibeliaZ timed out at 20 genomes in the audit. Cactus needs its container and is for same-species sets. |
| Core-gene | `--msa-source coregene` after `annotate` (bakta or prokka), with `--coregene panaroo`, `ppanggolin`, `scarap` or `orthofinder` | A set where a supported, model-based tree is wanted and the reference biases a SNP tree. How much of the core survives depends on `--aligner-arg threshold=`, the identity of Panaroo's clustering pass (0.98 by default), and on how wide the set is: on three Wolbachia genomes the core went from 52 genes at 0.98 to 638 at 0.90, and four across supergroups shared none at the default (verification.md). `ppanggolin` clusters proteins at 80% identity and selected 398 single-copy families on 24 Wolbachia genomes where Panaroo's defaults gave none, so it is a candidate for a divergent set; its declared limit for large sets is not measured. `scarap` infers the core directly through profile searches and found 577 genes on the same set, with a tree sharing 20 of 21 splits with the PPanGGOLiN one. `orthofinder` compares all proteomes against all and found 585 genes, its tree sharing 19 and 20 of 21 splits with the other two; its all-against-all search grows roughly with the square of the genome count, so it suits small and medium sets. | Orthologous genes present in `--core-threshold` of the genomes, each aligned with mafft and concatenated; Panaroo's graph cleaning (`--clean-mode`) and entropy filter remove annotation errors and badly aligned genes. Invariant sites kept, no ascertainment correction. Slow at hundreds of genomes (not measured here); the annotation is cached per genome. Within a species a SNP tree draws on far more positions. |
| Core-SNP | `--msa-source snptype` with `--snptyper simple`, `snippy`, `parsnp` or `ska2` | One species or a clonal set. | `simple`, `snippy` and `parsnp` map every genome to one reference, so the reference biases the result. `simple` maps with minimap2 `-x asm20` and records reference positions a genome does not cover as N, so a divergent outgroup is largely N (its coverage is in the log). `ska2` is reference-free. |

Reference bias. The aligners and the mapping typers project every genome
onto one reference, and the choice shapes mostly the within-cluster branching
order: swapping the reference changed 72 of 94 splits on a 50-genome set while
the deep structure held (`docs/audit/scaling-audit.md`).
Unless `--reference` names a genome, the stage uses the medoid, the genome
with the largest summed sketch ANI to the others, computed from the workdir
sketches (`--sketch` on the genome stages, or `repgenr sketch`) with
sourmash; ties go to the first filename, so the choice is the same on every
run. Only existing sketches are read; the choice never sketches. The cached
alignment is tied to the chosen genome, so it is rebuilt when the medoid
changes. Under Nextflow, sketches reach `phylo-build` on the bacterial path
with `params.sketch` (`PHYLO`, or `PHYLO_MSA` under `--phylo_split_msa`);
the viral and reads paths have no sketches and use the first genome. The
medoid is a less arbitrary choice than the first filename and is expected to
reduce the bias, not to remove it: in a set with an over-represented
genotype it sits inside that genotype. Without sketches or
sourmash, or with `--reference first`, the alphabetically first genome is
used, with a warning. On the F. tularensis set that first genome of the
whole genus was another species, and 5.1% of its columns varied inside a
species usually described as clonal; for a species run, name a genome of
that species. The chosen genome and the rule are recorded as
`reference_genome` and `reference_rule` and in `reference.txt` beside the
alignment (genome, rule and, for the medoid, its mean ANI).

Masking. `--mask gubbins` needs `--msa-source snptype` (in `run`, `--with-snptype` also accepts `--mask`) and a typer that writes a
whole-genome alignment, which `simple`, `snippy` and `parsnp` do. It replaces
the core-SNP alignment with Gubbins' filtered polymorphic sites. `ska2` writes
variable sites only, so it cannot be masked. Gubbins is for isolates of one
species. The masker estimates the variable fraction of the alignment, warns
above 10%, and repeats the figure if Gubbins fails. Gubbins needs a
multi-threaded RAxML build when given more than one thread; see
[usage.md](usage.md#snp-typing-and-masking).

## 6. Tree builders

| Builder | Input | Scaling | Resolution | Rooting by the builder |
|---|---|---|---|---|
| `iqtree` | MSA | Declared 500. Poor beyond that. Without `--aligner-arg model=...` it runs ModelFinder first; `model=GTR+G` skips that pass. With a SNP alignment the constant-site counts are passed as `-fconst`; `asc=lewis` selects `+ASC` instead. | Good on both diverse and clonal sets. Use at least 1000 ultrafast bootstrap replicates for supports; without `-B` the tree has none, and the stage says so. | Roots on the outgroup. |
| `raxmlng` | MSA | Declared 1000. Slow beyond that. | Supports from its own bootstrap procedure. With a SNP alignment the constant-site counts go on the model as `+ASC_STAM`; `asc=lewis` selects `+ASC_LEWIS` instead. | Roots on the outgroup. |
| `fasttree` | MSA | Declared 5000. Good to about 5000. `--aligner-arg gamma=true` rescales the branch lengths under the Gamma20 likelihood. | Approximate. Weaker resolution on clonal sets. | Unrooted. |
| `veryfasttree` | MSA | Declared 20000, not measured here; the same algorithm as `fasttree` with multi-core and vectorised arithmetic. `auto` picks it above 5000 leaves when an alignment is requested and never below, since `fasttree`'s tighter limit wins there; without an alignment request `auto` stays with the alignment-free builders once the set exceeds the aligner's limit. `-threads` follows `-t`; `--aligner-arg gamma=true` and `--aligner-arg double_precision=true` are extras. | As `fasttree`. | Unrooted. |
| `mashtree` | Genomes | Declared 10000. 53 min and 7.6 GB measured at 5000. | Poor on clonal sets. No supports. | Unrooted. |
| `sourmash` | Genomes | Declared 2000. The joining step runs rapidNJ when installed (`nj=auto`), else a cubic Python step refused above 5000. | Poor on clonal sets. No supports. | Unrooted. |
| `skani` | Genomes | Declared 10000, not measured here. ANI from sparse chaining of k-mer sketches, reliable to about 80% identity; pairs below the screen (`screen=80`) get distance 1 and are counted in the log. `preset=fast|medium|slow|small-genomes`, `min_af`, `nj` as extras. Also provides the distance matrix of the viral outgroup selection, with the `small-genomes` preset by default there; when no candidate is within skani's reach no matrix is written and the run proceeds without an outgroup. | Expected to resolve closer pairs than mash distance and more distant ones than Jaccard (skani, Shaw and Yu 2023), not measured here; no supports, and a clonal set resolves only as well as ANI does. | Unrooted. |
| `ska2` | Genomes | Declared 5000, as the typer. `ska distance` SNP counts (`metric=snps`) or k-mer mismatch fractions (`metric=mismatch`) over the split k-mers each pair shares, joined with rapidNJ when installed. | Good on clonal and outbreak sets, where it is reference-free; undercounts distances once pairs share few k-mers, so it warns on a set spanning several species. No supports. | Unrooted. |

`--treebuilder auto` picks the builder whose declared scale fits the leaf
count, the outgroup included. When `--msa-source snptype` or `--mask` is
given, the choice is made among the builders that consume an alignment, so
the request is honoured rather than dropped. When nothing asks for an
alignment and the leaf count exceeds the selected aligner's declared limit
(500 for the default `progressivemauve`), the choice is made among the
alignment-free builders, since a builder that needs an alignment would make
that aligner run on a set past its scale; the log names the aligner and its
limit when this applies.

The stage adds the outgroup to the input and passes it to builders that can
root. Whatever the builder, `phylo` then roots `tree/tree.nwk` itself
(`--root`): on the outgroup branch by default when an outgroup is staged, so
the trees of `fasttree`, `mashtree` and `sourmash` are rooted as well, at the
midpoint with `--root midpoint`, or not at all with `--root none`. The
builder's own output is kept as `tree/tree.unrooted.nwk`. On an
outgroup-rooted tree the taxonomy from `tree2tax` separates the outgroup from
one ingroup clade. See [usage.md](usage.md#rooting-the-tree).

`tree2tax` can collapse weak splits before it names nodes:
`--collapse-length L` merges a node whose branch is shorter than `L` in the
tree's own units, and `--collapse-support S` merges a node whose support is
below the fraction `S`. IQ-TREE and RAxML-NG supports are normalised
automatically. `mashtree` and the sourmash builder write no supports, so
`--collapse-support` does nothing for them and the stage warns.
`--collapse-diameter D` keeps the fewest clades whose leaf-to-leaf distance
stays within `D` and merges everything inside each, which maps a distance
cutoff onto the taxonomy without supports; the unit is the tree's own (a
mash distance of about 0.05 at 95 % ANI for mashtree, one minus the Jaccard
similarity for the sourmash builder, where 95 % ANI is near 0.88 at k = 31,
substitutions per variable site for the ML builders). The root, the outgroup
split and the leaves never collapse. See
[usage.md](usage.md#collapsing-weak-splits-in-tree2tax).

## 7. Assemblers and polishers

These apply to the reads chain only (`reads`, `assemble`).

| Sequencing platform | Assembler (`--assembler auto`) | Polisher (`--polisher auto`) |
|---|---|---|
| Illumina | `skesa`. `shovill` (SPAdes) is the alternative for paired reads. | None |
| Oxford Nanopore | `flye` | `medaka` (`racon` when medaka is not installed) |
| PacBio CLR | `flye` | `racon` with minimap2 overlaps (`--polish-rounds`) |
| PacBio HiFi | `flye` | None |

Notes:

- Medaka needs the basecaller model. Dorado names it in the read headers. Reads
  mirrored through SRA have rewritten headers, so the adapter assumes ONT's
  bacterial R10.4.1 model and records that. Give another with `--tool-arg
  model=...`.
- On one ONT run (SRR28800588) polishing did not change the CheckM2 scores,
  because the Flye draft was already gene-complete at that resolution.
  Unpolished long-read assemblies still carry indel errors that matter for SNP
  typing and alignment.
- The ENA selection filters that matter most are `--max-bases`, which drops
  whole-host libraries above a size, and `--drop-selection`, which drops
  libraries by ENA library selection. Its default is `MDA`: amplified libraries
  assemble into chimeric, uneven contigs. The 12% contamination seen on one
  library came from MDA amplification, which the default now drops.
- `--classifier auto` (the default) verifies each assembly's organism with
  sourmash gather against a GTDB sketch when one is given; without a sketch the
  check is skipped and the log says so. Each gather against the GTDB rs226
  representatives sketch holds about 0.6 GB, so the number run at once is
  bounded by `--memory-gb` as well as `--threads`. When GTDB has moved a
  species to another genus of the same family (family and epithet agree), the
  genome is flagged `genus_renamed` rather than `classifier_disagrees`.
- `shovill` needs a read pair. A run that ENA labels PAIRED but lists with one
  FASTQ file is planned as single-end, so `--assembler shovill` excuses it
  before downloading; `auto` assembles it with `skesa`.
- Changing only the CheckM2 gate (`--min-completeness`,
  `--max-contamination`) reuses the scores stored per run; CheckM2 runs again
  only for contigs, a database or a CheckM2 version it has not scored.
  `--force` does not bypass the stored scores; delete
  `assemblies/<run>/checkm2.json` to score a run again.
- `--one-per-sample` (the default) keeps one run per sample. A long-read run is
  preferred when it has at least 100 Mb and a tenth of the largest short-read
  run.

## 8. CLI `run` or Nextflow

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

## 9. How `auto` chooses

Each adapter declares `recommended_max_genomes`, an integer or none. `--tool
auto` (dereplicate) and `--treebuilder auto` (phylo) use it as follows:

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
   `--msa-source snptype`, `--mask` or `--msa`) and the leaf count exceeds
   the selected aligner's limit, steps 1 to 4 run over the alignment-free
   builders alone.

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

`repgenr list-tools` prints the declared limits. The table below is checked
against the registered adapters by `tests/unit/test_docs_tool_limits.py`.

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

These limits count genomes only. They do not model genome size, memory or
thread count. Aligners and SNP typers declare limits that the scale warnings
use, but limit-based `auto` selection applies only to dereplicators and tree builders.

## 10. Which environment provides each tool

A tool you choose must be installed. The conda files under `envs/` split the
tools over a core environment and six satellites, because they do not solve
as one environment (see [install.md](install.md#1-core-and-satellite-environments)).
Append a satellite's `bin` directory after core on `PATH`, never before it,
or give it to its tool alone with `--bin-dir TOOL=DIR` (see
[usage.md](usage.md#per-tool-directories---bin-dir)).

| Environment | Tools |
|---|---|
| `repgenr` (`envs/core.yml`) | skder, galah, sourmash, drep, sibeliaz, simple, ska2, iqtree, fasttree, veryfasttree, raxmlng, skesa, shovill, flye, medaka, racon, and `datasets` for the downloads |
| `repgenr-gubbins` | gubbins |
| `repgenr-mashtree` | mashtree |
| `repgenr-snippy` | snippy |
| `repgenr-parsnp` | parsnp |
| `repgenr-checkm2` | checkm2 (the reads chain quality step) |
| `repgenr-mauve` | progressivemauve |
| none (image only) | cactus |
