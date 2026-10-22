# Phylogeny

This page describes the phylogeny routes of `repgenr phylo` in detail:
choosing a route, SNP typing and recombination masking, annotation, the
core-gene route and its builders, support values, rooting, and the reuse of
an alignment. The short form, with one example per route, is in
[usage.md](usage.md#phylogeny), and the comparison of routes and tools is in
[choosing-tools.md](choosing-tools.md#phylogeny-routes).

## Contents

1. [Choosing a route](#choosing-a-route): `phylo-advise` and presets
2. [Genome count, leaf check and record names](#genome-count-leaf-check-and-record-names)
3. [SNP typing and masking](#snp-typing-and-masking): invariant sites, soft
   core, Gubbins, the `simple` typer, ska2
4. [Annotation](#annotation)
5. [Core-gene phylogeny](#core-gene-phylogeny): gene selection, trimming,
   partitions, species trees, protein alignments, the builders, the marker
   database
6. [Support values](#support-values) and [rooting the tree](#rooting-the-tree)
7. [Reusing an alignment across tree builders](#reusing-an-alignment-across-tree-builders)

## Choosing a route

`phylo` builds the tree either from an MSA (multiple sequence alignment) or
directly from the genomes. The MSA is a whole-genome alignment (`--aligner`),
a core-SNP alignment (`--msa-source snptype`) or a concatenated core-gene
alignment (`--msa-source coregene`). A tree builder such as `mashtree`,
`skani` or `ska2` works from the genomes without an alignment.

### Asking for a route: `phylo-advise`

`repgenr phylo-advise -wd $WD` reads the working directory and prints which
row of the [By purpose and divergence](choosing-tools.md#by-purpose-and-divergence)
table applies, the reasons, and the `repgenr phylo` command with `--preset`.
The same command follows with the preset written out as its flags. It writes
nothing to the working directory and sketches nothing; a manifest of an older
layout is migrated in place when it is opened, as by any other command.

The command counts the genomes (the representatives, or every genome with
`--all-genomes`) and the species and genera in canonical record names. It
takes the minimum and median pairwise ANI by `sourmash compare` over the
existing genome sketches (`--threads` sets its threads). It also checks
whether every genome is annotated and whether an outgroup is staged. Without
sketches the ANI is reported as unknown and the advice rests on the names.

The rules, in order:

1. A minimum ANI of 99.5 % or more gives `outbreak`.
2. One species by name, or, with other names, a minimum ANI of 95 % or more,
   gives `species-masked` for up to 200 genomes with a minimum ANI of 97 % or
   more, else `species`.
3. A minimum ANI under 80 % or several genera gives `family`.
4. Several species or a minimum ANI of 80 to 95 % gives `genus`.
5. With neither ANI nor canonical names the advice is `quick`, an
   alignment-free `skani` tree.

The command then adapts the advice:

- without an outgroup it adds `--root midpoint`;
- above the declared IQ-TREE limit of 500 genomes it adds `--treebuilder
  fasttree`, with `--aligner-arg gamma=true` in place of an IQ-TREE model;
  `-B` stays, as the number of resamples behind FastTree's SH-like supports;
- a core-gene route without annotations adds a reason to run `repgenr
  annotate` first;
- `family` without `REPGENR_MARKERS_DB` set adds a reason to download the
  marker database.

`--json` prints the same as one `repgenr.phylo-advise/1` object. The advice
follows a table of measured and declared figures. It is a starting point and
does not replace an inspection of the set.

```text
$ repgenr phylo-advise -wd $WD --all-genomes
Working directory: /data/wd
Genomes:      4 (all genomes)
Taxonomy:     names are not canonical; species not counted
Pairwise ANI: min 98.81 %, median 99.21 %
Annotations:  missing or incomplete
Outgroup:     none

Preset: species-masked (Genotype tree within a species, recombination masked)
Route:  snptype
As the species row, with recombination masked by Gubbins before the tree is
built.

Reasons:
  - The minimum pairwise ANI is 98.81 %, at or above 95 %: one species.
  - With 4 genomes and a minimum ANI of 98.81 %, Gubbins can mask
    recombination: it is for closely related isolates and is refused above
    12 % variable columns.
  - No outgroup is staged, so the tree is rooted at its midpoint (add `--root
    midpoint`); with an outgroup the default roots on it.

Command:
  repgenr phylo -wd /data/wd --all-genomes --preset species-masked --root midpoint
Expands to:
  repgenr phylo -wd /data/wd --all-genomes --msa-source snptype --snptyper simple --core-fraction 0.95 --treebuilder iqtree -B 1000 --mask gubbins --root midpoint
```

### Presets

`--preset NAME` on `phylo`, `phylo-build` and `run` stands for the flags of
one row of the table that `phylo-advise` applies (`RECOMMENDATIONS` in
`src/repgenr/stages/advise.py`). The reasons for each row are in
[By purpose and divergence](choosing-tools.md#by-purpose-and-divergence).

| Preset | Flags |
|--------|-------|
| `outbreak` | `--msa-source snptype --snptyper ska2 --treebuilder iqtree --aligner-arg model=GTR+G -B 1000` |
| `species` | `--msa-source snptype --snptyper simple --core-fraction 0.95 --treebuilder iqtree -B 1000` |
| `species-masked` | the `species` flags and `--mask gubbins` |
| `genus` | `--msa-source coregene --coregene ppanggolin --core-threshold 0.95 --max-core-genes 300 --gene-rank informative --partitioned --treebuilder iqtree -B 1000` |
| `family` | the `genus` flags with `--coregene markers` in place of `--coregene ppanggolin` (protein is implied) |
| `quick` | `--treebuilder skani` |

```bash
repgenr annotate -wd $WD
repgenr phylo -wd $WD --preset genus
repgenr phylo -wd $WD --preset family --gene-trim clipkit
repgenr run -wd $WD -l species -ts "Francisella tularensis" --preset species-masked
```

The preset is read as if its flags were given before the user's own:

- **An option given on the command line wins** over the preset's value, also
  at its default value (`--preset species -B 0` builds without bootstrap).
  Options left at their defaults take the preset's value.
- **`--aligner-arg` values add to the preset's.** A key given twice takes the
  user's value, so `--preset outbreak --aligner-arg model=GTR+G+I` changes the
  model.
- **A `--treebuilder` other than the preset's drops the tuning that builder
  does not read**: an `--aligner-arg` key outside its accepted parameters,
  and `--partitioned` for a builder without partition support. `--preset
  genus --treebuilder raxmlng` keeps `--partitioned`, while `--treebuilder
  fasttree` drops it and the IQ-TREE model. `-B` is kept.
- **A `--msa-source` other than the preset's is refused**, since the preset's
  settings belong to its route.
- **A preset's settings cannot be removed while keeping its route.** There is
  no `--no-partitioned` and no way to unset `--core-threshold`. To leave one
  out, write the other flags out without `--preset`.

Per-gene trimming (`--gene-trim clipkit`) is not part of any preset; add it
to `genus` or `family` when wanted. On `run`, the `genus` and `family`
presets add the `annotate` stage, as `--msa-source coregene` does.
`species-masked` masks in phylo's own typing pass, so `--with-snptype` is not
needed. The phylo stage record names the preset (`preset`). The MSA key and
the resume fingerprint read only the settings it expanded to, so a preset and
the same flags written out reuse each other's alignment and resume each
other. In the Nextflow pipeline `--preset` is passed inside `phylo_args` (see
[nextflow.md](nextflow.md#core-gene-route)).

## Genome count, leaf check and record names

`phylo` and `phylo-build` need at least three ingroup genomes (the outgroup
is not counted). A smaller set exits 3 with a message before any tool runs.

After the tree is built, its leaves are compared with the input genomes,
ingroup and outgroup. A tree builder can drop a genome it considers
degenerate and still exit 0, as mashtree does. A missing or unexpected leaf
then exits 3 with the names and the builder. The tree is kept in `tree/` for
inspection, and `phylo` is not recorded as completed. Leaf names are compared
without a FASTA extension or ParSNP's `.ref` reference marker, and characters
other than letters, digits, `_` and `-` are read as `_`, since some tools
rewrite them. `tree2tax` uses leaf names as written, so a tree that passes
this check can still leave leaves that `tree2tax` does not map to a genome.

Every alignment record (the SNP typers' and the aligners'
`align/msa.fasta`), every tree leaf and the outgroup name given to the
masker use the record name. That is the genome file name without its FASTA
suffix and `.gz` (`x.fasta.gz` gives `x`), as in `clusters.tsv` and
tree2tax.
`parsnp` names its records by file name and marks the reference with `.ref`;
the typer renames them to the genome names.

## SNP typing and masking

A SNP typer produces a core-SNP alignment, the variable sites among the
genomes. It suits closely related isolates. `repgenr snptype`,
`phylo --msa-source snptype` and `phylo-build --msa-source snptype` call the
typer. `phylo` types the genomes itself under `tree/msa/`, outgroup included,
so no separate `snptype` run is needed.

```bash
repgenr phylo -wd $WD --msa-source snptype --snptyper simple \
  --core-fraction 0.95 --mask gubbins --treebuilder iqtree -B 1000
```

Four typers are built in: `simple` (minimap2, samtools and bcftools, the
default), `snippy`, `parsnp` and `ska2`. The first three map every genome to
one reference and also write the whole-genome alignment a masker needs.
Without `--reference` the reference is the medoid genome by sketch ANI
when the working directory's sketches and sourmash are at hand, else the
first genome with a warning; the stage record names it.
Recombination masking (`--mask gubbins`) runs on the typer's whole-genome
alignment and replaces the core-SNP alignment with Gubbins' filtered
polymorphic sites. Typers that only emit variable sites cannot be masked.

Typer tuning keys are given as `--tool-arg KEY=VALUE` on `snptype`. On
`phylo`, `phylo-build` and `run` the same keys are given as `--aligner-arg
KEY=VALUE`, because those commands pass one set of tuning keys to every tool
they drive. The examples below use the form of the command they name.

### Invariant sites and branch lengths

A tree from a variable-site-only alignment overestimates branch lengths
unless the builder is told how many invariant sites the alignment lost. When
the typer writes a whole-genome alignment (`simple`, `snippy`, `parsnp`),
the typing pass counts its constant A, C, G and T columns into
`alignment_stats.json` beside `core_snp.fasta`, and `phylo` passes them on.
`ska2` has no whole-genome alignment, but it also writes an all-sites
alignment (`ska align --filter no-filter`, every shared split k-mer middle
base). The counts are taken from it under the same `--core-fraction` rule,
and the file is then removed, since it is as large as the core genome.

- IQ-TREE receives the counts as `-fconst`, RAxML-NG as `+ASC_STAM{...}` on
  its model. Both are first given the alignment without any invariant column
  that masking left, and those columns join the counts.
- `--aligner-arg asc=lewis` uses the Lewis correction (`+ASC`, `+ASC_LEWIS`)
  instead of the counts. A `model=` that already carries an ASC term is
  passed as it is, without the counts.
- FastTree applies no correction and says so.

The counts apply only to the alignment they were taken for (the sidecar
names it and holds its digest), so a tree built from `full_alignment.fasta`
is not corrected twice. In the Nextflow split mode the alignment task
publishes the sidecar and the tree task reads it. The counts are taken
before masking, so after `--mask` they are a slight undercount.

The `ska2` counts cover only the split k-mers shared by at least `min_freq`
of the genomes, so they describe the shared part of the genomes. ska's
`no-ambig` filter also drops a variable column that holds an ambiguity code
from `core_snp.fasta`, while a constant column with one still counts as
constant, so the corrected branch lengths are a slight underestimate. With
`lo=true` (`--tool-arg lo=true` on `snptype`, `--aligner-arg lo=true` on
`phylo`), `ska lo` writes variants only. That route has no counts, its tree
has no correction, and `phylo` warns. Compare topologies and supports rather
than lengths there, or build the tree from a whole-genome alignment
(`phylo-build --msa snp/full_alignment.fasta`).

### Soft core (`--core-fraction`)

The variable-site alignment of the mapping typers keeps every column in
which two bases occur, however many genomes lack a base there. A position
covered in two genomes out of seventy is therefore in it. `--core-fraction
F` on `snptype`, `phylo`, `phylo-build` and `run` keeps only the columns
where at least the fraction F of the genomes have a called base (a soft
core; 0.95 is the usual choice). It is applied after masking, and the
constant sites are counted under the same rule.

It is unset by default: every variable column is kept, and the log reports
how many columns a fraction of 0.95 would keep, so the effect can be seen
before the rule is chosen. The value is part of the alignment reuse key and
of `alignment_stats.json`.

Three cases differ from the description above:

- `snippy`'s core already holds only the sites called in every genome, so
  the filter changes nothing there, while the constant sites are still
  counted at F.
- `ska2` already keeps the split k-mers present in a fraction `min_freq` of
  the samples (0.9 by default), so `--core-fraction` is a second threshold
  on top of it.
- After `--mask`, the bases the masker replaced count as missing. A column
  inside a recombination block shared by more than 1 - F of the genomes is
  therefore dropped, while the constant sites come from the unmasked
  alignment.

`--core-fraction` applies only with `--msa-source snptype` (or, on `run`, to
the standalone `snptype` stage under `--with-snptype`). It is refused with an
aligner alignment, and a tree builder that works from the genomes warns that
it has no effect.

### Recombination masking with Gubbins

Gubbins expects isolates of one species. The masker estimates how much of
the alignment is variable, warns above 10 %, and repeats the figure if
Gubbins fails. On a genus-level set its scan can die with a bus error (see
[verification.md](verification.md#gubbins-divergence-limits-2026-09-10)).

**Refusals.** Masking is refused (exit 2) when more than 12 % of the columns
are variable, near where Gubbins has failed. It is also refused when the
genome names (canonical `Family_Genus_species_Accession` names only, read
where the last token is an accession) span more than one species; this
check runs before SNP typing. `--aligner-arg gubbins_divergence_limit=FRACTION`
moves the first limit, and `--aligner-arg gubbins_force=true` runs Gubbins
anyway, with a warning. Alignments shorter than 1000 columns are not judged.

**Gubbins' tree builders.** Gubbins builds a tree in every iteration, by
default with RAxML, and with more than one thread it needs a multi-threaded
RAxML build (`raxmlHPC-PTHREADS*`). Some conda builds ship only the
single-threaded binary, and Gubbins then exits before its first iteration.
When repgenr runs Gubbins natively and finds no such build, it switches to
IQ-TREE (or to one thread when IQ-TREE is missing too) and says so in the
log. `--aligner-arg gubbins_tree_builder=raxmlng`,
`--aligner-arg gubbins_first_tree_builder=rapidnj` and
`--aligner-arg gubbins_args="--min-snps 5"` pass the builder, the
first-iteration builder and further `run_gubbins.py` arguments through.

Unless a first-iteration builder or model is given, repgenr asks for
`rapidnj` with the JC model when rapidnj is available (it ships with
Gubbins). A rapidnj first tree makes the first iteration cheap, while later
iterations still use the resolved builder; with IQ-TREE as Gubbins' builder,
tree building took most of the run time on a genus-level set (see
verification.md). Gubbins' own final tree, its recombination-aware
maximum-likelihood tree, is kept as `masker_tree.nwk` beside the masked
alignment for comparison with the tree repgenr builds from the filtered
sites.

**Gapped taxa.** Gubbins leaves taxa with more than 25 % gaps or N
(`--filter-percentage`) out of its recombination analysis but still writes
them to its outputs, unmasked. The `simple` typer writes N where a genome
does not align, so repgenr passes `--filter-percentage 100`. A value given in
`gubbins_args` takes precedence, and a taxon above it is then refused
(exit 3) before Gubbins runs, with its name.

### The `simple` typer

**Threads and scratch.** The `simple` typer maps each genome independently,
so `--threads` buys concurrent genomes first, and threads inside one genome's
chain only when there are more threads than genomes. Its per-genome
intermediates are written compressed and removed as soon as that genome's
consensus has been read, so scratch stays at a few hundred megabytes
whatever the genome count. A genome whose chain fails keeps its
intermediates for inspection. Under a container back-end each genome's chain
of tools runs in a single container, so a genome costs one engine start
rather than eight.

**Missing sequence.** The typer builds each genome's consensus by applying
its SNP calls to the reference. It then sets to N every reference position
where none of the genome's primary or supplementary alignments places a
base: positions outside the alignments, and positions within a deletion
(CIGAR D or N operations), as read from the minimap2 SAM. Sequence a genome
lacks is therefore missing data, not the reference base. A column of the
core-SNP alignment (`snp/core_snp.fasta`, or `tree/msa/core_snp.fasta` from
the typing pass of `phylo`) is kept only when at least two of A, C, G and T
occur in it. N, other ambiguity codes and gaps do not make a column
variable, and they are written as N.

**Distances.** The distance matrix (`snp_distance_matrix.tsv` beside it)
counts, for each pair, the differing sites among those where both genomes
have a base, so two genomes are not separated by a region one of them lacks.
A pair that shares no such site has the distance `NA`. Distances between
genomes of different gene content are computed over different numbers of
sites, so they are not directly comparable as proportions.

**Coverage.** The log has one line per genome with the fraction of the
reference's bases that genome covers, and a summary line (minimum, median,
maximum). A genome below 50 % is warned about. A genome with no base at any
core SNP site, usually one that did not align to the reference at all, is
refused (exit 3) with its name, since tree builders refuse a sequence of N
only.

**minimap2 preset.** The typer maps assemblies with the minimap2 preset
`asm20`, which suits assembly-to-reference alignment up to several percent
divergence. A genome more divergent than `asm20` tolerates, such as an
outgroup from another species, aligns over less of the reference and is
correspondingly more N; check the coverage lines in the log. `--tool-arg
preset=` accepts `asm5`, `asm10`, `asm20`, `map-ont`, `map-pb`, `map-hifi`,
`sr` and `none` (minimap2's default settings, no `-x`), and refuses any other
value (exit 2). `asm5` suits near-identical genomes. `none` aligns more
divergent sequence at the cost of distance accuracy.

### Split k-mers: ska2

`--snptyper ska2` (`--tool ska2` on `snptype`; split k-mer analysis) is
reference-free. Every genome is an ordinary sample, so no assembly's private
errors bias the SNP distances, and the medoid reference rule does not apply.
It emits a variable-site alignment and no whole-genome alignment, so it is
not compatible with `--mask`. A second `ska align` call with `--filter
no-filter` gives the constant-site counts for the tree builder (see
[Invariant sites and branch lengths](#invariant-sites-and-branch-lengths));
that alignment is counted and removed, not published. Tune it with
`ksize=31` and `min_freq=0.9` (the fraction of samples a split k-mer must
occur in). ska2 is designed for clonal and outbreak sets, where the shared
k-mer space is large, and the stage says so when the genome names span more
than one species.

`lo=true` replaces `ska align` with `ska lo`, which recovers variants in more
diverged samples by local assembly of the split k-mer graph (Horsfield et al.
2025). The indels it finds are published as `snp/variants.vcf`, `lo_missing`
(0.1) is its limit on missing data per site, and `min_freq` is not read in
this mode. `ska lo` writes variants only, so this mode has no constant-site
counts and the tree's branch lengths are not corrected; the log says so.

`--treebuilder ska2` builds a tree without an alignment from the same split
k-mers. `ska distance` counts the SNPs among the k-mers each pair shares
(`metric=snps`, the default, branch lengths in SNPs) or the fraction of
mismatching k-mers (`metric=mismatch`). The matrix is joined with rapidNJ
when installed. Pairs sharing few k-mers have their SNPs undercounted, so
the same divergence warning applies (the outgroup is left out of the count).
A pair that shares none is refused, since ska writes it as distance 0.
`min_freq` here is `ska distance`'s own (0 by default). `auto` never picks
this builder; name it.

## Annotation

`repgenr annotate -wd WD` annotates the genome set once, for the core-gene
route and any later use of the gene calls. It annotates the representatives
(or every genome with `--all-genomes`) and the outgroup, one directory per
record under `annotations/`, with the GFF3 (sequences appended), the
proteins and the gene nucleotides, and an index `annotations.tsv` (see
[output.md](output.md#annotate)).

```bash
repgenr annotate -wd $WD --bakta-db $BAKTA_DB
```

**Annotators.** The annotator is `bakta` by default, which needs its
database. Download it once with `bakta_db download --output DIR --type
light` (a 1.4 GB download, about 4 GB unpacked) and pass `--bakta-db DIR` or
set `BAKTA_DB`. `--type full` is 30 GB to download and about 84 GB
unpacked, and annotates CDS more specifically. The database is part of the
stage's resume inputs through its `version.json`, so replacing it
re-annotates. `--tool prokka` needs no database. Tuning goes through
`--tool-arg` (`skip_crispr=true`, `genus=Francisella`,
`translation_table=11` for bakta; `kingdom=Bacteria`, `gcode=11` for
prokka).

**Parallel genomes.** An annotator uses many threads in only some of its
steps, so the stage annotates several genomes at once and splits `-t`
between them. `-p/--num-processes` sets how many. The default (0) takes
about one genome per four threads, capped by the usable CPUs, so `-t 16`
gives four genomes with four threads each. Each running genome is a
separate annotator process; allow a few GB of memory for each with bakta.
Annotation is the slow step of the core-gene route (times per genome are in
[verification.md](verification.md#annotation-2026-10-15)), so a genus-level
set is a cluster job, where the Nextflow route runs one task per genome.
`phylo-build`, when it annotates its own inputs, divides its `-t` in the same
way. `repgenr run --msa-source coregene` passes its `-p` to the annotate
stage.

**Resume.** A later run annotates only the genomes that are new or whose
genome, annotator, tool version, database version or tuning changed, and
reuses the rest. The index is written as each genome finishes, so an
interrupted or failed run keeps the genomes that completed. After a failure,
the genomes already running finish, no new ones start, and the first error
is reported. A record whose genome left the set is kept, since annotation is
expensive and the genome may return after another dereplication. Only a
record whose genome has left the working directory is removed. `--force`
re-annotates everything.

The stateless step `annotate-build --genome FILE -o DIR` annotates one
genome.

## Core-gene phylogeny

`phylo --msa-source coregene` builds the tree from a concatenated alignment
of the core genes, the genes present in nearly all genomes of the set.

```bash
repgenr annotate -wd $WD
repgenr phylo -wd $WD --msa-source coregene --treebuilder iqtree
```

On a working directory, run `annotate` first, as above. `repgenr run
--msa-source coregene` inserts the `annotate` stage before `phylo`
(`--annotator`, `--bakta-db`). The stateless step `phylo-build --msa-source
coregene` reads `--annotations-dir` or annotates its inputs itself.

A core-gene builder (`--coregene`) clusters the genes of the annotated
genomes, and the core genes are aligned and concatenated. Five builders are
available: `panaroo` (the default), `ppanggolin`, `scarap`, `orthofinder`
and `markers` (see [Core-gene builders](#core-gene-builders)). The alignment
keeps its invariant sites, so the tree builder applies no ascertainment-bias
correction and the branch lengths are substitutions per site of the core.

**Per-gene alignments.** Every builder hands over one alignment per core
gene, and the stage ranks, caps and concatenates them in the same way for
all five. Panaroo's are the genes of its filtered (or, with
`filtered=false`, unfiltered) alignment, from `aligned_gene_sequences/`, and
PPanGGOLiN's are its per-family alignments (`msa_<partition>_dna/`). Where
the others come from is listed in
[choosing-tools.md](choosing-tools.md#core-gene-builders).

**The outgroup.** The outgroup is annotated and included. Being another
species, it lowers the core, so on fewer than about twenty genomes lower
`--core-threshold` or leave the outgroup out with `--no-outgroup --root
midpoint`.

**Reuse and outputs.** The alignment is reused across tree-builder changes
like the other sources, and rebuilt when an annotation changed (the reuse
stamp holds the index digest). The outputs are
`tree/msa/core_gene_alignment.fasta`, the builder's tables under
`tree/msa/pangenome/`, `tree/msa/core_genes.tsv`,
`tree/msa/partitions.json` and `tree/msa/core_gene_stats.json` (see
[output.md](output.md#phylo)).

### Choosing the genes

By default every core gene goes into the alignment. On a large set that can
be tens of thousands of genes, and a tree builder gains little from the
slowly changing ones. Four options on `phylo`, `phylo-build` and `run`
select the genes:

- `--min-gene-length N` drops genes with fewer than N alignment columns.
- `--max-gene-gaps F` drops genes where more than the fraction F of the
  alignment is gap characters.
- `--gene-rank` orders the remaining genes: `occupancy` (the fraction of
  genomes holding the gene, the default), `informative` (parsimony-informative
  sites), `variable` (variable sites), `length` or `random` (seeded, so a
  rerun gives the same genes). Ties are broken by gene length, then gene id.
- `--max-core-genes N` keeps at most the first N genes of that order; 0
  keeps all (the default).

The kept genes are concatenated in gene-id order, not in rank order, so a
different `--gene-rank` that keeps the same set of genes gives the same
alignment. `tree/msa/core_genes.tsv` records the decision for every gene the
builder aligned. Every per-gene alignment is held in memory during the
selection, so the memory it needs is about the size of the alignments:
roughly 1 GB at 500 genomes and 1500 genes.

The options apply only with `--msa-source coregene` and are refused with the
other sources. A genome that holds none of the kept genes would be missing
from the alignment and the tree, so the stage stops before the tree is built
and names it. Lower `--core-threshold`, relax the selection, or leave the
genome out of the set.

For an overview of a genus, a cap of a few hundred genes ranked by
`informative` is a suggested starting point; it is not derived from a
measurement. Check that the tree is stable against a larger cap on the set
at hand.

```bash
repgenr phylo -wd $WD --msa-source coregene \
    --max-core-genes 300 --gene-rank informative --min-gene-length 300
```

### Trimming

`--gene-trim clipkit` (on `phylo`, `phylo-build` and `run`; `none` is the
default) trims each per-gene alignment with
[ClipKIT](https://github.com/JLSteenwyk/ClipKIT) before the genes are
measured, ranked and concatenated. Across a genus or wider set the ends of a
gene alignment are often ragged and some columns poorly aligned. Trimming
removes them, so they do not inflate the variable and informative site
counts that the gene filters and `--gene-rank` read.

- ClipKIT runs in its `smart-gap` mode (`clipkit IN -o OUT -m smart-gap`),
  the only mode offered on the command line.
- ClipKIT is given the sequence type. A nucleotide (codon) gene alignment is
  trimmed with `-s nt -co`, which removes whole codons so the reading frame
  is kept, and a protein one with `-s aa`. This needs ClipKIT 2.1.2 or later.
- A gene whose alignment loses every column is dropped, and the log says how
  many. When every gene is dropped the stage stops.
- Trimming is part of the reuse key of the core-gene alignment, so turning
  it on or off rebuilds the alignment. The phylo stage record holds
  `gene_trim`, and `core_gene_stats.json` the trimming counts.

`clipkit` is in `envs/core.yml`. Under a container back-end it runs from
`quay.io/biocontainers/clipkit:2.14.0--pyhdfd78af_0`.

```bash
repgenr phylo -wd $WD --msa-source coregene --coregene ppanggolin \
    --gene-trim clipkit --max-core-genes 300 --gene-rank informative
```

### Partitioned models

By default the tree builder fits one substitution model to the whole
concatenation. With `--partitioned` (on `phylo`, `phylo-build` and `run`) it
fits one model per core gene instead, from the column ranges in
`tree/msa/partitions.json`. The partitions are edge-linked: the genes share
one tree and one set of branch lengths, scaled per gene. Genes evolve at
different rates and under different substitution patterns, so this can
describe a concatenation better than a single model. It also adds
parameters and run time.

- `iqtree` writes `tree/partitions.txt`, a NEXUS `sets` block with one
  `charset <gene> = <start>-<end>;` per gene and no model, and is run with
  `-p` (proportional branch lengths). The model of every partition is the
  one passed as `-m`: the `--aligner-arg model=` given, or `GTR+G` (`LG+G`
  for a protein alignment) without it, in which case ModelFinder does not
  run. `--aligner-arg model=MFP` selects a model for each partition with
  ModelFinder, which is slow with hundreds of partitions.
- `raxmlng` writes `tree/partitions.txt` with the model on every line
  (`GTR+G, <gene> = <start>-<end>`, `LG+G` for a protein alignment, or the
  `--aligner-arg model=` given) and is run with `--model partitions.txt`.
  Branch lengths are scaled per partition, RAxML-NG's default.
- `aster` does not fit a partition model: it builds one tree per partition
  and summarises them (see [Species tree from gene trees](#species-tree-from-gene-trees)).
- Other tree builders read no partitions; the stage warns that the builder
  ignores them and fits one model.

Gene names that the tools would not parse are rewritten: characters other
than letters, digits, `_` and `-` become `_`, and a repeated name gains a
numeric suffix. `--partitioned` is refused in these cases:

- without `partitions.json` beside the alignment, or when it describes
  another file (partitions exist only for a builder that hands over one
  alignment per gene);
- together with an ascertainment-bias correction (`--aligner-arg asc=`, a
  model with `+ASC`), which removes the invariant columns and would shift
  the partition coordinates;
- outside `--msa-source coregene`, except on `phylo-build --msa FILE`, where
  the `partitions.json` beside the file decides.

`--partitioned` changes the tree, not the alignment. A rerun that only adds
or removes it reuses the core-gene alignment and rebuilds the tree. The
phylo stage record holds `partitioned`.

```bash
repgenr phylo -wd $WD --msa-source coregene --max-core-genes 300 \
    --gene-rank informative --treebuilder iqtree --partitioned -B 1000
```

### Species tree from gene trees

`--treebuilder aster` infers one tree per core gene and summarises the gene
trees into a species tree with ASTER (github.com/chaoszhang/ASTER, the
maintained successor of ASTRAL). It reads the gene partitions, so it needs
`--msa-source coregene --partitioned`. Without them it is refused, and
`auto` never picks it.

1. The alignment is cut into one FASTA per gene (`tree/genes/`). A genome
   with no called character in a gene is left out of that gene. A gene with
   fewer than 4 remaining genomes, or shorter than `--aligner-arg
   min_gene_columns=` columns, is skipped and counted in the log.
2. Each gene tree is inferred with one thread, and the genes run in batches
   over `-t` threads. The default is IQ-TREE with `--fast`, `GTR+G` (`LG+G`
   for a protein alignment, or `--aligner-arg model=`) and aBayes supports.
   `--aligner-arg gene_tree_builder=fasttree` uses FastTree (`-nt -gtr`, or
   `-lg`) with its SH-like local supports. The gene trees are collected in
   `tree/gene_trees.nwk`, one per line.
3. ASTER summarises them (`--aligner-arg program=`):
   - `wastral` (default): weighted ASTRAL, which weights each gene-tree
     quartet by its support and branch length. Internal branch lengths are
     in coalescent units and terminal branches have none, so `--root
     midpoint` is refused before any gene tree is built; root on the
     outgroup instead.
   - `astral4`: unweighted ASTRAL-IV, with branch lengths in substitutions
     per site.
   - `astral-pro3`: ASTRAL-Pro3; with one record per genome it gives the
     same kind of tree as `astral4`.

   `--aligner-arg rounds=N` sets ASTER's rounds of placement and
   subsampling (`-r N -s N`; 4 each by default).

Node labels are local posterior probabilities in [0, 1]. They are written
without `-B`, and tree2tax `--collapse-support` reads them on that scale.
The outgroup, when one is staged, is a genome like the others in the gene
trees, and `phylo --root` roots the species tree afterwards.

```bash
repgenr phylo -wd $WD --msa-source coregene --max-core-genes 500 \
    --gene-rank informative --treebuilder aster --partitioned
```

### Protein alignments

By default the core-gene alignment is nucleotide: each gene is aligned as
codons. Across a genus or a family the third codon positions saturate, and
a nucleotide model then describes the deeper branches poorly.
`--coregene-alphabet protein` (on `phylo`, `phylo-build` and `run`; `dna` is
the default) aligns the core genes as amino acids instead. Codon models are
not offered.

- `ppanggolin` runs `ppanggolin msa --source protein` and hands over the
  family alignments under `msa_<partition>_protein/`.
- `scarap` and `orthofinder` align the proteins of each orthogroup with
  mafft, as for nucleotides, and hand over those alignments without
  back-translating them.
- `panaroo` aligns nucleotides only, so a protein alphabet with `--coregene
  panaroo` (the default builder) is refused when the command starts, before
  any stage runs.
- `markers` aligns protein only, since its models are protein profiles.
  Without `--coregene-alphabet` it uses protein (the log says so), and
  `--coregene-alphabet dna` is refused when the command starts.

The gene selection reads `X` as missing data in a protein alignment (`N` in
a nucleotide one) when it counts variable and informative sites.

The tree builder is told the alphabet and uses an amino-acid model:

| Tree builder | Nucleotide | Protein |
|---|---|---|
| `iqtree` | `-st DNA` | `-st AA`; ModelFinder chooses among protein models unless `--aligner-arg model=` is given |
| `raxmlng` | `GTR+G` | `LG+G`, on every partition line with `--partitioned` |
| `fasttree`, `veryfasttree` | `-nt -gtr` | `-lg` |

An ascertainment-bias correction (`--aligner-arg asc=`, a model with `+ASC`)
is refused on a protein alignment: it applies to nucleotide SNP alignments,
and a core-gene alignment keeps its invariant sites.

The alphabet is recorded in `partitions.json` (field `alphabet`) and in
`core_gene_stats.json`. When `partitions.json` describes the alignment, the
tree builder takes the alphabet from it, so `phylo-build --msa FILE` (the
tree half of the Nextflow split) needs no `--coregene-alphabet` of its own.
Without the sidecar it follows `--coregene-alphabet`. A protein alignment
given as `--msa FILE` with no sidecar and an MSA source other than
`coregene` is treated as nucleotide, and IQ-TREE then fails on `-st DNA`.
The alphabet is part of the reuse key of the core-gene alignment, so
changing it rebuilds the alignment. The phylo stage record holds
`coregene_alphabet`.

```bash
repgenr phylo -wd $WD --msa-source coregene --coregene ppanggolin \
    --coregene-alphabet protein --treebuilder iqtree -B 1000
```

### Core-gene builders

The builders differ in how they group genes, and that grouping decides how
many core genes a set retains. How they compare, with measurements on
Wolbachia, is in
[choosing-tools.md](choosing-tools.md#core-gene-builders) and
[verification.md](verification.md#core-gene-route). The settings below are
given as `--aligner-arg KEY=VALUE` unless an option is named.

**Panaroo** (`--coregene panaroo`, the default) clusters the genes into a
pangenome and cleans annotation errors from the graph (`--clean-mode
strict`, the default, `moderate` or `sensitive`). It aligns every gene
present in at least `--core-threshold` of the genomes (0.95) with mafft,
concatenates them, and removes the genes its entropy filter flags
(`filtered=false` keeps them).

- The clustering identity (`threshold=`) is the setting to adjust first.
  Panaroo groups genes in a first pass at this identity, 0.98 by default, so
  on a divergent set most orthologues of a gene end up in separate clusters
  and never reach the core. On three Wolbachia genomes of one supergroup the
  genes in all three went from 52 at 0.98 to 638 at 0.90.
- When `threshold` is not set and the canonical genome names
  (`Family_Genus_species_Accession`) span more than one species, the builder
  uses 0.95 and logs the species count and the identity chosen. A set of one
  species, or names that are not canonical, keep 0.98. The value used is
  passed to Panaroo on every run and recorded in the stage record.
- A wider set needs a lower identity. Too low a value starts to group genes
  that are not orthologues (RIBAP, Genome Biology 2024), so lower it in
  steps and check the result. Lowering `--core-threshold` instead recovers
  fewer genes and gives a softer core.
- The builder warns when the core comes out under 100 genes, unless
  `threshold` was set explicitly. `core_gene_stats.json` carries the count
  and the identity used; the count is the number to check before relying on
  the tree.

**PPanGGOLiN** (`--coregene ppanggolin`) clusters the proteins with MMseqs2
at 80 % identity over 80 % of their length (`identity=`, `coverage=`),
partitions the pangenome graph, and concatenates the DNA alignments of the
selected families.

- The selection defaults to the soft core at `--core-threshold`, so the gene
  set means the same as with Panaroo. `partition=core` takes the strict
  core, and `partition=persistent` the tool's statistical partition, which
  does not shrink as genomes are added but needs a set of some size
  (PPanGGOLiN warns below 15 genomes).
- PPanGGOLiN aligns only genomes with exactly one copy of a family, so a
  multi-copy genome is a gap for it. Families multi-copy in at least 5 % of
  the genomes (`dup_margin=`) are left out by default, and
  `single_copy=false` keeps them.
- A genome that shares no selected family is dropped from the alignment with
  a warning, and an outgroup in that state is an error.
- `--clean-mode` and `filtered=` are Panaroo's and have no effect here. They
  are left out of the resume fingerprint and the alignment stamp, so
  changing them reruns nothing.

**SCARAP** (`--coregene scarap`) infers the core directly, without the whole
pangenome. It builds orthogroups on a set of seed genomes (`seeds=`, 100),
turns them into MMseqs2 profiles and searches every genome against them. A
gene is core when it is present in a single copy in at least
`--core-threshold` of the genomes.

- It reads the protein and gene files the annotate stage writes beside each
  GFF3, not the GFF3 itself. The adapter aligns and back-translates each
  orthogroup as for OrthoFinder; `scarap concat` is not used, because it
  deletes its per-orthogroup alignments.
- The orthogroups come from profile searches, so there is no identity to
  lower. `method=` chooses the orthogroup inference (`FH`, the default).
- The seed genomes are drawn at random with no fixed random seed, so on a
  set larger than `seeds` two runs can select slightly different cores. The
  final selection is checked over all the genomes, not only the seeds.
- `method=O-D` and `method=O-B` build the seed pangenome with OrthoFinder
  (DIAMOND or BLAST search). They run natively only, from the environment
  `envs/scarap.yml` builds. The SCARAP image has no OrthoFinder, and SCARAP
  1.0.2 reads OrthoFinder's result with a pandas call that pandas 2 removed,
  so that environment keeps pandas older than 2. They take several times
  longer than `FH`.

**OrthoFinder** (`--coregene orthofinder`) compares every proteome against
every other with DIAMOND (`search=`, also `diamond_ultra_sens`, `blastp` or
`mmseqs`), clusters the hits into orthogroups with MCL (`inflation=`, 1.2)
and stops there (`-og`).

- The adapter keeps the orthogroups present in a single copy in at least
  `--core-threshold` of the genomes, aligns each one's proteins with mafft
  and back-translates the alignment to codons from the `.ffn` files.
- A gene whose length does not match its protein is left out of its
  orthogroup (a gap in the concatenation), and the log counts them.
  Selenocysteine and pyrrolysine, which mafft does not accept, are aligned
  as X.
- With `search=mmseqs` OrthoFinder uses fixed file names under `/tmp`, so
  two such runs on one host can collide.
- The search is all against all, so its cost grows roughly with the square
  of the genome count.

### GTDB marker genes

`--coregene markers` builds the alignment from the single-copy marker
proteins of GTDB's reference trees: the 120 bacterial families (`bac120`) or
the 53 archaeal ones (`ar53`). It clusters nothing, so it is not limited by
how far orthology clustering reaches. It suits a set that spans a family or
more, where the clustering builders find few core genes and DNA alignments
saturate. It is the builder of `--preset family`, reads the `.faa` files of
the annotation, and needs the marker database.

**The marker database.** `repgenr markers-db` downloads the HMMER3 model of
each of the 168 distinct families, as listed by GTDB-Tk, from the InterPro
API (`--threads` concurrent requests, 8 by default). The full database is
about 25 MB.

```bash
repgenr markers-db download --output ~/db/repgenr_markers   # --set bac120|ar53|all
export REPGENR_MARKERS_DB=~/db/repgenr_markers
repgenr markers-db info ~/db/repgenr_markers
repgenr annotate -wd $WD
repgenr phylo -wd $WD --msa-source coregene --coregene markers \
    --markers-db ~/db/repgenr_markers --partitioned -B 1000
```

- The model versions are InterPro's current ones, so they can be newer than
  the versions GTDB-Tk pins (for example `PF00380.26` against GTDB-Tk's
  `PF00380.20`). `markers.tsv` records each family's version and the cutoff
  it is searched with: the gathering cutoff (GA) for Pfam and the noise
  cutoff (NC) for TIGRFAM, as GTDB-Tk does.
- The database is written beside `--output` and renamed into place, so an
  interrupted download leaves no partial database. A download stopped while
  it swaps a new database for an earlier one leaves the earlier one beside
  `--output`; the next download puts it back and removes such leftovers
  before it starts.
- A second download into the same directory fetches only the families it
  lacks. `repgenr doctor` reports when a database that a stage recorded has
  moved or changed version.

**The builder** follows GTDB-Tk's rules:

- **Search.** Each genome's proteins are searched with `hmmsearch`: the Pfam
  families with their gathering cutoff (`--cut_ga`) and the TIGRFAM families
  with their noise cutoff (`--cut_nc`). These are two calls per genome, run
  in parallel over the genomes, with the set's HMM file split by cutoff in
  scratch.
- **One model per gene.** A protein that several marker models hit is
  assigned to the model with the highest bit score.
- **One gene per marker.** A marker found once in a genome is *unique*.
  Found several times with identical protein sequences it is *multiple
  unique*, and the best hit is kept. Found several times with different
  sequences it is *multiple* and treated as absent from that genome, as
  GTDB-Tk does. Otherwise it is *missing*.
- **Alignment.** The proteins of each marker are aligned to its model with
  `hmmalign`, all genomes in one call. Only the match-state columns (the `x`
  columns of the `#=GC RF` line) are kept, so each marker contributes as
  many columns as its model has match states. GTDB-Tk aligns one sequence at
  a time; the match-state columns are the same.

A marker enters the alignment when it is single-copy in at least
`--core-threshold` of the genomes (0.95 by default). `--max-core-genes`,
`--gene-rank`, `--gene-trim` and `--partitioned` (one `LG+G` partition per
marker) then apply as for the other builders. The alignment is protein (see
[Protein alignments](#protein-alignments)).

- **The genome guard.** A genome with fewer single-copy markers than
  `min_markers` (`--aligner-arg min_markers=N`; 10 % of the set by default,
  12 for bac120 and 6 for ar53, the marker form of GTDB-Tk's `--min_perc_aa
  10`) is refused by name before any alignment. The message suggests that
  the set may be of the wrong domain or the annotation incomplete. A genome
  with fewer than half of the set's markers is named in a warning.
- **The domain.** `--marker-set` is `bac120`, `ar53` or `auto` (the
  default). With `auto` the genomes are searched with bac120 first. When
  fewer than half of them reach its minimum they are searched with ar53 as
  well, and the set with the higher median fraction of its markers found is
  used. The fraction rather than the count is compared, since the sets
  differ in size. The set used is recorded in `core_gene_stats.json`
  (`builder_settings.marker_set`).
- **Tables.** `tree/msa/pangenome/markers_summary.tsv` and
  `markers_per_genome.tsv` hold the calls per genome (see
  [output.md](output.md#phylo)).
- **Reuse and records.** The reuse key of the alignment holds the database
  version (the sets and a digest of `markers.tsv`, which names every model
  version) and the marker set, but not the database path. A moved copy of
  the same database therefore reuses the alignment. The resume fingerprint
  hashes the version in the same way. The phylo record names the database
  path (`markers_db`) and its version (`markers_db_version`), which `repgenr
  doctor` compares with the database found there.

**Differences from GTDB-Tk.** The proteins are those of the annotation
(bakta or prokka), whereas GTDB-Tk calls genes with Prodigal, so a protein
can differ at its ends. The marker models are InterPro's current versions.
The alignment is not masked with GTDB's column mask; `--gene-trim clipkit`
trims each marker instead. The genome guard counts markers, where GTDB-Tk
counts alignment columns.

## Support values

Support values on an IQ-TREE tree come from `-B N` (long form `--bootstrap
N`), N ultrafast bootstrap replicates (at least 1000). `--aligner-arg
alrt=N` adds the SH-aLRT test with N replicates (`-alrt`, also at least
1000), which gives supports without `-B`. With both, IQ-TREE writes the two
values in each node label. Without either the tree carries no supports, and
`phylo` says so.

Without `--aligner-arg model=` IQ-TREE runs ModelFinder before the tree
search. Above 200 genomes the log notes that `model=GTR+G` skips that pass.
With `--partitioned` a fixed model is applied instead (see
[Partitioned models](#partitioned-models)).

## Rooting the tree

`phylo` (and `phylo-build`, and `run`) roots `tree/tree.nwk` after any tree
builder has finished. The rooting is chosen with `--root`:

- `outgroup` roots on the branch leading to the outgroup, which is split in
  two halves, so the root has two children: the outgroup and the whole
  ingroup. This is the default when an outgroup is staged and
  `--no-outgroup` is not given.
- `midpoint` roots at the midpoint of the longest leaf-to-leaf path. It
  needs branch lengths and does not use the outgroup, so it is the choice
  for a tree built with `--no-outgroup`.
- `none` keeps the tree as the builder wrote it. This is the default without
  an outgroup.

`--root outgroup` without a staged outgroup, or together with
`--no-outgroup`, exits 2 before any tool runs.

Builders that are given the outgroup (`iqtree -o`, `raxmlng --outgroup`)
place it at a root with three children, and such a tree is rerooted like any
other, with one INFO line. A tree whose root already has the outgroup as one
of two children is left as written. Branch lengths and support values are
kept: each support value stays with the split it was computed for. A
midpoint that falls on an internal branch splits that branch in two, and its
support value is written on both new root branches, since both describe the
same split. When the tree is rerooted, the builder's own output is kept as
`tree/tree.unrooted.nwk`. The stage record holds the rooting applied
(`root`) and the outgroup leaf (`outgroup`). How `tree2tax` treats each
rooting is described under [tree2tax](usage.md#tree2tax).

## Reusing an alignment across tree builders

`phylo` stamps the alignment it builds (`align/msa_source.json` or
`tree/msa/msa_source.json`) with what produced it: the source and its
settings, the genome set, and the alignment's own digest. A later `phylo`
run that changes only the tree builder, the bootstrap or the thread count
reuses that alignment instead of aligning or SNP-calling again, and says so
in the log. Anything the alignment depends on, such as the aligner, the SNP
typer, `--reference`, `--mask` or the genome set, rebuilds it, as does
`--force` or a change to the alignment file itself.

The stateless step offers the same split. `phylo-build --msa-only` builds the
alignment and writes `msa.fasta` without a tree, and `phylo-build --msa
<file>` builds a tree from an alignment an earlier call produced. The
Nextflow layer uses these to run the alignment and the tree as separate
tasks. `--msa-only` writes the `alignment_stats.json` and `partitions.json`
of the alignment beside `msa.fasta`, under that file's name and digest, and
`--msa <file>` reads them from beside the file it is given. `--msa-only`
leaves the SNP typing pass's `tree/msa/` directory (and `scratch/`) in the
step's output directory beside `msa.fasta`.
