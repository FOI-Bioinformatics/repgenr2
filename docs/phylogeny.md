# Phylogeny

This page describes the phylogeny routes of `repgenr phylo`: choosing a
route, SNP typing and recombination masking, annotation, the core-gene route
and its builders, support values, rooting, and the reuse of an alignment.
The short form is in [usage.md](usage.md#phylogeny), and the comparison of
routes and tools in [choosing-tools.md](choosing-tools.md#phylogeny-routes).

## Contents

1. [Choosing a route](#choosing-a-route): `phylo-advise` and presets
2. [Genome count, leaf check and record names](#genome-count-leaf-check-and-record-names)
3. [SNP typing and masking](#snp-typing-and-masking)
4. [Annotation](#annotation)
5. [Core-gene phylogeny](#core-gene-phylogeny)
6. [Tree-builder options](#tree-builder-options), [support values](#support-values) and [rooting the tree](#rooting-the-tree)
7. [Reusing an alignment across tree builders](#reusing-an-alignment-across-tree-builders)

## Choosing a route

`phylo` builds the tree from a whole-genome alignment (`--aligner`), a
core-SNP alignment (`--msa-source snptype`), a concatenated core-gene
alignment (`--msa-source coregene`), or directly from the genomes with a
builder such as `mashtree`, `skani` or `ska2`.

### Asking for a route: `phylo-advise`

`repgenr phylo-advise -wd $WD` prints the row of the
[By purpose and divergence](choosing-tools.md#by-purpose-and-divergence)
table that applies, the reasons, and a `phylo` command with `--preset`,
followed by the same command with the preset written out as its flags. It
reads the genome names and the minimum pairwise ANI from the existing genome
sketches, and writes nothing. The rules, in order:

1. a minimum ANI of 99.5 % or more gives `outbreak`;
2. one species (by name, or a minimum ANI of 95 % or more) gives
   `species-masked` for up to 200 genomes at 97 % or more, else `species`;
3. a minimum ANI under 80 % or several genera gives `family`;
4. several species or 80 to 95 % gives `genus`;
5. without ANI or canonical names it gives `quick`.

Without an outgroup it adds `--root midpoint`, and above the declared
IQ-TREE limit of 500 genomes `--treebuilder fasttree` with `--aligner-arg
gamma=true`. A core-gene route without annotations adds a reason to run
`repgenr annotate`, and `family` without `REPGENR_MARKERS_DB` a reason to
download the marker database. Without sketches the ANI is unknown and the
advice rests on the names. `--json` prints the advice as one object (see
[output.md](output.md#machine-readable-status)). The advice is a starting point, not a
replacement for inspecting the set.

### Presets

`--preset NAME` on `phylo`, `phylo-build` and `run` stands for these flags:

| Preset | Flags |
|--------|-------|
| `outbreak` | `--msa-source snptype --snptyper ska2 --treebuilder iqtree --aligner-arg model=GTR+G -B 1000` |
| `species` | `--msa-source snptype --snptyper simple --core-fraction 0.95 --treebuilder iqtree -B 1000` |
| `species-masked` | the `species` flags and `--mask gubbins` |
| `genus` | `--msa-source coregene --coregene ppanggolin --core-threshold 0.95 --max-core-genes 300 --gene-rank informative --partitioned --treebuilder iqtree -B 1000` |
| `family` | the `genus` flags with `--coregene markers` |
| `quick` | `--treebuilder skani` |

The preset is read as if its flags came before the user's own:

- An option given on the command line wins, also at its default value
  (`--preset species -B 0` builds without bootstrap).
- `--aligner-arg` values add to the preset's, and a key given twice takes
  the user's value.
- Another `--treebuilder` drops the preset's `--aligner-arg` keys and
  `--partitioned` where that builder does not read them; `-B` is kept.
- Another `--msa-source` is refused, and a preset setting cannot be removed;
  write the flags out without `--preset` instead.

A preset and the same flags written out reuse each other's alignment and
resume each other, and the stage record names the preset. On `run`, `genus`
and `family` add the `annotate` stage, and `species-masked` masks in
phylo's own typing pass. `--gene-trim clipkit` is in no preset. In
Nextflow, `--preset` goes inside `phylo_args`.

## Genome count, leaf check and record names

`phylo` needs at least three ingroup genomes and exits 3 before any tool
runs on a smaller set. After the tree is built, its leaves are compared with
the input genomes. A builder can drop a genome it considers degenerate and
still succeed, as mashtree does, so a missing or unexpected leaf exits 3
and `phylo` is not recorded as completed; the tree is kept in `tree/` for
inspection. The comparison ignores FASTA
suffixes and ParSNP's `.ref` marker and reads other characters than letters,
digits, `_` and `-` as `_`. This matching is looser than tree2tax's, so a
tree that passes can still have leaves tree2tax does not map. Records and leaves are named by the genome file
name without its FASTA suffix and `.gz` (`x.fasta.gz` gives `x`), as in
`clusters.tsv`.

## SNP typing and masking

A SNP typer produces a core-SNP alignment, the variable sites among closely
related genomes. `phylo --msa-source snptype` types the genomes itself into
`tree/msa/`; `repgenr snptype` writes the same tables to `snp/`.

```bash
repgenr phylo -wd $WD --msa-source snptype --snptyper simple \
  --core-fraction 0.95 --mask gubbins --treebuilder iqtree -B 1000
```

Four typers are built in: `simple` (minimap2, samtools and bcftools, the
default), `snippy`, `parsnp` and `ska2`. The first three map every genome to
one reference: `--reference`, else the medoid genome by sketch ANI, else the
first genome with a warning. They also write the whole-genome alignment that
masking needs. Typer settings are `--tool-arg KEY=VALUE` on `snptype` and
`--aligner-arg KEY=VALUE` on `phylo`, `phylo-build` and `run`.

### Invariant sites and branch lengths

A variable-site alignment overestimates branch lengths unless the builder
knows how many constant sites were left out. The typing pass counts them
into `alignment_stats.json`, and `phylo` passes them to IQ-TREE (`-fconst`)
or RAxML-NG (`+ASC_STAM`). `--aligner-arg asc=lewis` uses the Lewis
correction instead, and a `model=` with its own ASC term is passed as it
is. FastTree applies no correction and says so. The counts are taken before
masking, so after `--mask` they slightly undercount. For `ska2` they come
from an all-sites alignment of the shared split k-mers, which is then
removed. `ska2 lo=true` writes variants only, so its tree has no correction
and `phylo` warns. Compare topologies rather than lengths there, or build
the tree from a whole-genome alignment (`phylo-build --msa
snp/full_alignment.fasta`) when branch lengths matter.

### Soft core (`--core-fraction`)

`--core-fraction F` keeps only the alignment columns where at least the
fraction F of the genomes have a called base; 0.95 is the usual choice. It
is unset by default, and the log reports how many columns 0.95 would keep.
It applies only with `--msa-source snptype` (and to `run --with-snptype`),
after masking, and is refused with an aligner alignment. Bases the masker
replaced count as missing. `snippy`'s core already holds only fully called
sites, and `ska2`'s `min_freq` (0.9) is a first threshold of the same kind.

### Recombination masking with Gubbins

`--mask gubbins` replaces the core-SNP alignment with Gubbins' filtered
polymorphic sites; `ska2` cannot be masked. Gubbins is for isolates of one
species. The masker warns above 10 % variable columns and refuses above
12 % (see [verification.md](verification.md#gubbins-divergence-limits-2026-09-10)),
and refuses genome names that span several species.
`--aligner-arg gubbins_divergence_limit=FRACTION` moves the limit and
`gubbins_force=true` overrides both refusals. Alignments under 1000
columns are not judged.

Gubbins needs a multi-threaded RAxML build for more than one thread. Without
one, repgenr switches Gubbins to IQ-TREE, or to one thread when IQ-TREE is
missing too, and says so. `gubbins_tree_builder=`,
`gubbins_first_tree_builder=` and `gubbins_args=` pass Gubbins' own options,
and repgenr asks for a rapidnj first tree when rapidnj is installed.
A rapidnj first tree makes the first iteration cheap. Gubbins' final,
recombination-aware tree is kept as `masker_tree.nwk` beside the masked
alignment. repgenr passes `--filter-percentage 100`, since the `simple`
typer writes N where a genome does not align. A lower value given in
`gubbins_args` refuses a taxon above it (exit 3).

### The `simple` typer

The `simple` typer maps genomes concurrently and keeps scratch small. It
sets to N every reference position where a genome has no aligned base, so
sequence a genome lacks is missing data rather than the reference base.
A column enters the core-SNP alignment only when at least two of A, C, G
and T occur in it. `snp_distance_matrix.tsv` counts differences only where
both genomes have a base, and `NA` for a pair without such a site. The log
gives each genome's coverage of the reference and warns below 50 %. A genome
with no base at any SNP site is refused (exit 3). A genome failing its chain
keeps its intermediates for inspection, and under a container back-end each
genome's chain runs in one container.

It maps with the minimap2 preset `asm20`, which suits assemblies up to
several percent divergence. `--tool-arg preset=` also accepts `asm5`,
`asm10`, `map-ont`, `map-pb`, `map-hifi`, `sr` and `none` (minimap2's
default settings). `asm5` suits near-identical genomes, and `none` aligns
more divergent sequence at some cost in accuracy. A divergent genome, such as an outgroup of another
species, aligns over less of the reference and is correspondingly more N.

### Split k-mers: ska2

`--snptyper ska2` is reference-free and suits clonal and outbreak sets; it
warns when the names span several species. Tune it with `ksize=31` and
`min_freq=0.9`. `lo=true` uses `ska lo`, which recovers variants in more
diverged samples (`lo_missing`, 0.1, limits missing data per site;
`min_freq` is not read in this mode).
`ska2` writes no whole-genome alignment, so it cannot be masked.
`--treebuilder ska2` builds a tree from `ska distance` without an alignment
(`metric=snps`, the default, or `mismatch`), joined with rapidNJ when
installed. Pairs that share few k-mers are undercounted, and a pair that
shares none is refused. `min_freq` there is `ska distance`'s own (0). `auto`
never picks this builder.

## Annotation

`repgenr annotate -wd WD` annotates the representatives (or every genome
with `--all-genomes`) and the outgroup once, for the core-gene route. Each
record gets a directory under `annotations/` (see
[output.md](output.md#annotate)).

```bash
repgenr annotate -wd $WD --bakta-db $BAKTA_DB
```

`bakta` (the default) needs its database: `bakta_db download --output DIR
--type light` (1.4 GB, about 4 GB unpacked; `--type full` is about 84 GB
unpacked; it annotates CDS more specifically), then `--bakta-db DIR` or
`BAKTA_DB`. Replacing the database
re-annotates. `--tool prokka` needs no database. Tuning goes through
`--tool-arg` (for example `genus=Francisella` for bakta, `kingdom=Bacteria`
for prokka). Several genomes are annotated at once, by default one per four
threads (`-p/--num-processes`). `repgenr run --msa-source coregene` passes
its `-p` to `annotate`, and `phylo-build` divides its `-t` the same way when
it annotates its own inputs. Allow a few GB of memory each with bakta. A
later run annotates only new or changed genomes and keeps records of genomes
that left the set. An interrupted run keeps the genomes that finished, and
`--force` annotates everything again. `annotate-build` annotates one genome
outside a working directory. Times per genome
are in [verification.md](verification.md#annotation-2026-10-15); a
genus-level set is a cluster job.

## Core-gene phylogeny

`--msa-source coregene` builds the tree from a concatenated alignment of
the core genes, the genes present in nearly all genomes. Run `annotate`
first; `repgenr run --msa-source coregene` inserts it (`--annotator`,
`--bakta-db`). `phylo-build --msa-source coregene` reads
`--annotations-dir` or annotates its inputs itself.

```bash
repgenr annotate -wd $WD
repgenr phylo -wd $WD --msa-source coregene --coregene ppanggolin --treebuilder iqtree -B 1000
```

A core-gene builder (`--coregene`) groups the genes, and every builder hands
over one alignment per core gene, which the stage ranks, caps and
concatenates. The alignment keeps its invariant sites. The outgroup is
included and lowers the core, so on fewer than about twenty genomes lower
`--core-threshold` or use `--no-outgroup --root midpoint`. An annotation
change rebuilds the alignment.

### Choosing the genes

- `--min-gene-length N` and `--max-gene-gaps F` drop short or gappy genes.
- `--gene-rank` orders the rest: `occupancy` (the default), `informative`,
  `variable`, `length` or `random`.
- `--max-core-genes N` keeps the first N; 0 keeps all (the default).

The kept genes are concatenated in gene-id order, and
`tree/msa/core_genes.tsv` records the decision for every gene. A genome that
holds none of the kept genes stops the stage; lower `--core-threshold`,
relax the selection, or leave the genome out of the set. Every gene alignment is held
in memory during the selection, roughly 1 GB at 500 genomes and 1500 genes.
The options are refused outside `--msa-source coregene`. For a genus
overview a cap of
a few hundred genes ranked by `informative` is a suggested start; check that
the tree is stable against a larger cap.

### Trimming

`--gene-trim clipkit` trims each gene alignment with ClipKIT (`smart-gap`
mode, codon-aware for DNA; ClipKIT 2.1.2 or later) before the genes are
ranked. Across a genus, ragged gene ends otherwise inflate the site counts
the ranking reads. A gene that loses every column is dropped. The default
is `none`.

### Partitioned models

`--partitioned` fits one model per core gene instead of one for the whole
concatenation, from `tree/msa/partitions.json`. The genes share one tree.

- `iqtree` is run with `-p` and a NEXUS partition file. Every partition gets
  the `--aligner-arg model=`, or `GTR+G` (`LG+G` for protein) without it;
  `model=MFP` runs ModelFinder per partition, which is slow.
- `raxmlng` gets the same models in its partition file.
- `aster` builds one tree per partition (see below).
- Other builders ignore the partitions with a warning.

`--partitioned` is refused without `partitions.json` beside the alignment,
with an ascertainment-bias correction (it would shift the partition
coordinates), and outside `--msa-source coregene` except on `phylo-build
--msa FILE`. It changes the tree, not the alignment, so adding it reuses the
alignment. Gene names the tools cannot parse are rewritten.

### Species tree from gene trees

`--treebuilder aster` (with `--msa-source coregene --partitioned`) infers one
tree per gene, by default IQ-TREE `--fast` with aBayes supports
(`gene_tree_builder=fasttree` uses FastTree), and summarises them with
ASTER. `program=wastral` (the default) weights quartets by support and
gives coalescent branch lengths, so `--root midpoint` is refused; root on
the outgroup instead.
`astral4` gives lengths in substitutions per site, and `astral-pro3` the
same kind of tree for one record per genome. A gene with fewer than
four genomes or `min_gene_columns` columns is skipped. `rounds=N` sets
ASTER's rounds of placement and subsampling (4). Node labels are local
posterior probabilities in [0, 1], which tree2tax `--collapse-support`
reads on that scale. `phylo --root` roots the species tree afterwards, and
`auto` never picks `aster`.

### Protein alignments

`--coregene-alphabet protein` aligns the core genes as amino acids, since
third codon positions saturate across a genus or family. `ppanggolin`,
`scarap` and `orthofinder` support it. `panaroo` is nucleotide only, and
`markers` is protein only. The tree builder then uses an amino-acid model:

| Tree builder | Nucleotide | Protein |
|---|---|---|
| `iqtree` | `-st DNA` | `-st AA`, ModelFinder over protein models unless `model=` is given |
| `raxmlng` | `GTR+G` | `LG+G` |
| `fasttree`, `veryfasttree` | `-nt -gtr` | `-lg` |

`panaroo` with protein, and `markers` with `dna`, are refused when the
command starts. An ascertainment-bias correction is refused on a protein
alignment. The alphabet is recorded in `partitions.json`, so `phylo-build
--msa FILE` takes it from there. A protein `--msa FILE` without that file is read as
nucleotide and IQ-TREE fails on `-st DNA`, so keep `partitions.json` beside
it. Changing the alphabet rebuilds the alignment.

### Core-gene builders

How the builders compare is in
[choosing-tools.md](choosing-tools.md#core-gene-builders), and their
measurements in [verification.md](verification.md#core-gene-route).
Settings are `--aligner-arg KEY=VALUE` unless an option is named.

- **Panaroo** (`--coregene panaroo`, the default) clusters genes at the
  identity `threshold=`, 0.98 by default and 0.95 when the canonical names
  span several species. This identity governs the core: a wider set needs a
  lower value (on three Wolbachia genomes the core grew from 52 genes at
  0.98 to 638 at 0.90). Lower it in steps, since too low a value groups
  non-orthologues. `--clean-mode` (`strict`, `moderate`, `sensitive`) cleans
  the pangenome graph, and `filtered=false` keeps the genes the entropy
  filter flags. Read the core-gene count and the identity used in
  `tree/msa/core_gene_stats.json` before relying on the tree. The builder
  warns when the core has fewer than 100 genes, unless `threshold=` was
  set.
- **PPanGGOLiN** (`--coregene ppanggolin`) clusters proteins with MMseqs2 at
  80 % identity and coverage (`identity=`, `coverage=`). It takes the soft
  core at `--core-threshold`; `partition=core` or `persistent` choose
  otherwise; `persistent` needs a set of some size (PPanGGOLiN warns below
  15 genomes). PPanGGOLiN aligns only single-copy genomes, so families
  multi-copy in at least 5 % of the genomes are left out (`dup_margin=`;
  `single_copy=false` keeps them). A genome sharing no selected family is
  dropped with a warning, and an outgroup in that state is an error.
- **SCARAP and OrthoFinder** read the `.faa` and `.ffn` files beside each
  GFF3, which an `--annotations-dir` given to `phylo-build` must provide.
- **SCARAP** (`--coregene scarap`) builds orthogroups on `seeds=` (100)
  random seed genomes and searches every genome against them; `method=`
  (`FH`) chooses the inference. The `O-D` and `O-B` methods run natively
  only, from `envs/scarap.yml`, which keeps pandas older than 2 for SCARAP
  1.0.2. Seeds are drawn without a fixed random seed, so two runs on a set
  larger than `seeds` can select slightly different cores.
- **OrthoFinder** (`--coregene orthofinder`) searches all proteomes against
  all (`search=`, DIAMOND by default) and clusters with MCL (`inflation=`,
  1.2). The adapter aligns each orthogroup's proteins with mafft and
  back-translates them to codons; selenocysteine and pyrrolysine are aligned
  as X. With `search=mmseqs`, two runs on one host can collide in `/tmp`.
  Its cost grows with the square of the genome count.

### GTDB marker genes

`--coregene markers` aligns GTDB's single-copy marker proteins, bac120 or
ar53, found with HMMER. It needs no clustering, suits a family or wider, and
is the builder of `--preset family`. It needs the marker database:

```bash
repgenr markers-db download --output ~/db/repgenr_markers   # --set bac120|ar53|all
repgenr markers-db info ~/db/repgenr_markers
export REPGENR_MARKERS_DB=~/db/repgenr_markers
repgenr phylo -wd $WD --msa-source coregene --coregene markers --partitioned -B 1000
```

The download takes each family's HMMER3 model from InterPro (`--threads`
requests at once, 8 by default; about 25 MB in all) and fetches only the
missing families on a repeat. It is written beside `--output` and renamed
into place, so an interrupted download leaves no partial database. A
download stopped while replacing an earlier database leaves that database
beside `--output`, and the next download puts it back. The
model versions are InterPro's current ones and can be newer than those
GTDB-Tk pins. `--markers-db` names the database on `phylo`, and `repgenr
doctor` reports a recorded database that moved or changed version.

The builder follows GTDB-Tk's search rules: Pfam families with their
gathering cutoff and TIGRFAM families with their noise cutoff, each protein
assigned to its best model, and a marker found in several different copies
treated as absent. Each marker's proteins are aligned to its model with
`hmmalign`, keeping the match-state columns.
`--marker-set auto` (the default) searches with bac120, adds ar53 when fewer
than half of the genomes reach the bac120 minimum, and keeps the set with
the higher median fraction of its markers found. A marker enters the alignment when it is single-copy in at
least `--core-threshold` of the genomes. A genome with fewer single-copy
markers than `min_markers` (`--aligner-arg min_markers=N`; by default 10 %
of the set, 12 for bac120 and 6 for ar53) is refused, which
usually means a wrong domain or an incomplete annotation. A genome with
fewer than half of the markers is named in a warning. The database version
and the marker set are part of the alignment's reuse key; its path is not.
The proteins come from bakta or prokka rather than Prodigal, and the
alignment is not masked with GTDB's column mask (`--gene-trim clipkit`
trims each marker instead), so the tree can differ from a GTDB-Tk tree.

## Tree-builder options

Each tree builder takes its own settings through `--aligner-arg KEY=VALUE`.
The options of `iqtree` and `raxmlng` are described under
[Invariant sites and branch lengths](#invariant-sites-and-branch-lengths),
[Partitioned models](#partitioned-models) and
[Support values](#support-values), those of `ska2` under
[Split k-mers: ska2](#split-k-mers-ska2), and those of `aster` under
[Species tree from gene trees](#species-tree-from-gene-trees). The other
builders accept the following:

| Builder | Options |
|---|---|
| `fasttree` | `gamma=true` rescales the branch lengths under the Gamma20 likelihood. |
| `veryfasttree` | `gamma=true` as for `fasttree`; `double_precision=true`. Its thread count follows `-t`. |
| `mashtree` | `genomesize=`, passed to mashtree. |
| `sourmash` | `metric=ani` (default) or `jaccard`; `ksize=` (31) and `scaled=` (1000) of the sketches; `nj=auto`, `rapidnj` or `python` for the joining step. |
| `skani` | `screen=` (80, skani's `-s`), `min_af=`, `preset=fast`, `medium`, `slow` or `small-genomes`; `nj=` as for `sourmash`. |

The sourmash distance is one minus the ANI estimate c^(1/k) from the average
containment c, the quantity the sourmash dereplicator thresholds. A pair
with too few shared hashes gets distance 1 and is counted in the log (see
[How the tools estimate ANI](choosing-tools.md#how-the-tools-estimate-ani)).
For a set wider than a species, use `ksize=21`: k=21 is less biased than
k=31 at 80 to 95 % identity, and the builder logs this advice at k=31 when
the genome names span more than one genus. With `nj=auto` the joining step
runs rapidNJ when it is installed and a pure-Python step otherwise; the
Python step refuses more than 5000 leaves. A skani pair below the screen
gets distance 1 and is counted in the log. With
`vgenome --outgroup-treebuilder skani`, the viral outgroup selection uses
the skani distance matrix with the `small-genomes` preset. When no candidate is
within skani's reach, no matrix is written and the log suggests mashtree,
the default.

## Support values

On IQ-TREE, `-B N` (`--bootstrap`) gives ultrafast bootstrap supports and
`--aligner-arg alrt=N` SH-aLRT supports, each with at least 1000
replicates; with both, each node label holds the two values. Without either
the tree has no supports. Without
`--aligner-arg model=` IQ-TREE runs ModelFinder first; `model=GTR+G` skips
it.

## Rooting the tree

`--root` chooses how `tree/tree.nwk` is rooted after any builder:

- `outgroup` roots on the branch to the outgroup (the default when one is
  staged);
- `midpoint` roots at the midpoint of the longest path, for a tree built
  with `--no-outgroup`;
- `none` keeps the builder's root (the default without an outgroup).

`--root outgroup` without an outgroup exits 2. Builders given the outgroup
place it at a three-way root, which is rerooted like any other tree.
Supports stay with their splits, and a midpoint root on an internal branch
carries that branch's support on both root branches. When the tree is
rerooted, the builder's output is kept as `tree/tree.unrooted.nwk`, and the
stage record holds the rooting.

## Reusing an alignment across tree builders

`phylo` stamps the alignment it builds with what produced it. A later run
that changes only the tree builder, the bootstrap or the thread count reuses
it, and says so in the log. Anything the alignment depends on (the
aligner, the typer, `--reference`, `--mask`, the genome set) rebuilds it, as
does `--force` or a change to the alignment file. `phylo-build --msa-only` builds only the alignment, with
its `alignment_stats.json` and `partitions.json` beside it, and
`phylo-build --msa FILE` builds a tree from one and reads those files. The
Nextflow split mode uses the pair.
