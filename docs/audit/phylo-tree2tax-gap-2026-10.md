# Phylo and tree2tax gap analysis (2026-10)

Status: analysis complete. Batches A to D are merged. Of batch E, E2 to E4
shipped as #15, #16 and #18 (with the skani tree builder, #17); E1, the
benchmark sweep, needs cluster time and is deferred.

This record compares the `phylo` and `tree2tax` stages at commit 43823e8 with
current practice in bacterial phylogenomics and tree-based taxonomy, lists the
gaps that remain after the earlier audits, and orders the remedies by value
and effort. Evidence is cited by file and line of that commit, by earlier audit
records, or by the external source named in the text. Measurements come from
`docs/verification.md`, `docs/audit/scaling-audit.md` and the changelog; none
were repeated for this record.

## Sources

- Code: `src/repgenr/stages/phylo.py`, `stages/snptype.py`, `stages/tree2tax.py`,
  the adapters under `aligners/`, `snptypers/`, `maskers/`, `treebuilders/`,
  `tree/newick.py`, `tree/rooting.py`, `core/contracts.py`.
- Earlier records: `scaling-audit.md` (appendix B backlog, bias tests B4 and
  B5), `verification.md` (francisella_all runs, Gubbins failures, ska2 against
  mashtree), the changelog entries for PRs 139 to 263.
- External, as of October 2026: the IQ-TREE 3 and RAxML-NG documentation
  (ascertainment bias, `-fconst`, `ASC_STAM`), the VeryFastTree 4 release
  notes, the Gubbins 3 documentation (tree builders, `--first-tree-builder`),
  the SKA2 0.5 documentation and the `ska lo` description, the Core-SNP-filter
  documentation and the soft-core comparison it cites, relative evolutionary
  divergence as defined for GTDB (Parks et al. 2018; PhyloRank), TreeCluster
  (Balaban et al. 2019), Life Identification Numbers, PopPUNK (Lees et al.
  2019), and the FlexTaxD reader (`flextaxd/modules/ReadTaxonomy.py`, GitHub
  master and the PyPI 0.1.5 wheel).

## What the earlier audits already settled

These items are shipped and are not revisited here: rooting in `phylo --root`
with support labels that follow their splits (#263); alignment reuse through
the MSA stamp, the parallel simple typer and the vectorised core-SNP step
(#166 to #168); Gubbins masking of the ingroup only, `--filter-percentage 100`,
the thread cap and the builder fallback (#158, #163, #232, #235); FastTree
bootstrap, warnings for ignored parameters, version floors and container pins
(#147, #172, #178, #179); tree2tax collapse by support and length (#139), the
`--root-name` fix (#171) and the doctor checks of the two tables (#234); the
three-genome refusal, leaf checks and trailing-text handling (#207, #212, #213).

## Gaps in phylo

1. **No ascertainment-bias correction.** Trees are inferred from variable
   columns only. `treebuilders/iqtree.py:51-69` passes neither `-fconst` nor a
   model; `treebuilders/raxmlng.py:44-45` fixes `GTR+G` without an `ASC`
   term. Branch lengths are therefore inflated, as `docs/usage.md` already
   states. The counts of invariant A, C, G and T columns needed for `-fconst`
   fall out of the column scan in `snptypers/simple.py:383-412` and are
   discarded. Relative branch lengths from SNP-only data are known to deviate
   from the true values, and the `-fconst` route is preferred over `+ASC`
   because the latter aborts if any invariant column remains, which masking
   can produce.
2. **IQ-TREE model not settable.** `iqtree.py:18-26` accepts no extras and
   passes no `-m`, so IQ-TREE runs ModelFinder on every call, a model
   selection pass that a user who knows the model cannot skip. RAxML-NG
   already accepts `model=`.
3. **The simple typer's core is not a core.** `simple.py:383-412` keeps any
   column in which two bases are observed, however many genomes have no call
   there, so a position covered in two genomes out of seventy enters the
   alignment, and nothing reports how much of the alignment each genome
   lacks. Current practice applies an explicit completeness threshold (the
   `ska2` typer already does, through `ska align --min-freq`, default 0.9);
   the Core-SNP-filter documentation argues that a soft core with a stated
   threshold such as 95 % keeps many more informative sites than a strict
   core on large, diverse collections.
4. **Reference choice is alphabetical.** `stages/phylo.py:1081-1105` takes the
   first genome and warns. Test B4 measured the effect: swapping the reference
   changed 72 of 94 splits (normalised RF 0.77) on 50 genomes with
   progressiveMauve and FastTree (`scaling-audit.md`, bias section;
   `benchmarks/results/bias_b4.json`). The usual mitigation is the medoid by
   sketch distance; the workdir already holds sketches at k = 21, 31 and 51.
5. **Gubbins defaults are slow and fragile.** Its internal IQ-TREE pass took 76
   of 87 minutes where the rapidnj and fasttree pair took 6; the run ended in
   a bus error at 39 % variable sites, and on subsets the failure appears
   between 8 % and 13 %. Gubbins' own recombination-aware tree is not kept:
   only the filtered alignment is copied out, and the final tree is rebuilt
   from the polymorphic sites (`stages/snptype.py:148-166`).
6. **No validation before work starts.** `phylo_split_msa = true` with the
   shipped `--treebuilder mashtree` fails at run time (`phylo.py:417-422`,
   `nextflow/nextflow.config`); `auto` counts the ingroup only and ignores the
   MSA source (`phylo.py:426-443`); a bootstrap of 0 with a likelihood builder
   produces an unsupported tree without a word.
7. **Newick label sanitising is not injective.** `tree/newick.py:74-79` maps
   every reserved character to an underscore, so two labels can merge.
8. **Scale envelope.** The sourmash builder uses neighbour joining in pure
   Python (cubic, refused above 5000); there is no VeryFastTree adapter; the
   declared limits rest on few measurements (mashtree measured to 5000
   against a declared 10 000; SibeliaZ declared at 2000 and timed out at 20;
   sourmash measured).
9. **ska2 is used only as a typer.** `ska distance` and the 2025 `ska lo` mode
   are unused, and the divergence warning runs on the aligner path only.
10. **Stale statements.** `verification.md` (simple typer filling with the
    reference base, fixed in #232), `scaling-audit.md` (simple typer called
    serial), the `usage.md` SNP-phylogeny example (it runs `snptype` before
    `phylo` as if the two were chained, which they are not),
    `nextflow/nextflow.config` (progressiveMauve called all-pairs).

## Gaps in tree2tax

11. **The FlexTaxD contract is untested.** The reader in FlexTaxD's GitHub
    master requires a header containing `parent` and `child` in either order
    (repgenr's `child` then `parent` is accepted), reads an optional third
    column as a rank, reads the genome map without a header, keys nodes by
    label so identical labels merge, and pre-creates a node named `root` from
    which its tree validation starts. A rank column must be non-empty on
    every row: each row is stripped before it is split, so an empty third
    field disappears and the row silently inherits the previous row's rank.
    The PyPI 0.1.5 wheel is not a usable target: it skips the first line of
    the genome map as a header, so repgenr's headerless `genomes_map.tsv`
    loses its first genome; it treats a third column as another node; and its
    reader passes the parent column as the child of every link, whichever
    header order the file has, so the stored tree is inverted. No document
    names the FlexTaxD version repgenr targets, and with `--root-name` other
    than `root` the tree hangs off a node FlexTaxD does not start from, which
    its optional `--validate` reports as a node-count mismatch.
12. **Quadratic naming and path building.** Each internal node hashes all of
    its descendant leaves (`stages/tree2tax.py:436-438`) and each leaf walks its
    full ancestor chain (`:442-463`); ladder-shaped trees from clonal data make
    this quadratic (`scaling-audit.md`, appendix B).
13. **`--remove-outgroup` prunes a path, not the leaf.** On a midpoint-rooted
    tree the outgroup's parent keeps a single child and the ancestor hashes
    still contain the outgroup label (`tree2tax.py:137-138`).
14. **The stateless step lacks inputs.** `tree2tax-relations` has no
    `--segments` (`tree2tax.py:175-189`), so the Nextflow viral path loses
    segment accessions; `segments.tsv`, `outgroup_accession.txt` and
    `outgroup/` are not in the resume digests (`cli/base.py:306-310`);
    `node_basename` is not recorded.
15. **No clade-level threshold.** Every bifurcation becomes a taxon unless the
    support or length collapse applies. mashtree and sourmash trees carry no
    supports, so `--collapse-support 0.5` on the genus run collapsed nothing
    (`verification.md`). TreeCluster's maximum-diameter criterion yields the
    fewest clades whose patristic diameter stays under a threshold, which maps
    to an ANI-like cutoff; relative evolutionary divergence (Parks et al. 2018)
    is the accepted measure for placing ranks on a rooted tree.
16. **Node names are unstable.** A name is the hash of all descendant leaves,
    so any change of the leaf set renames every ancestor. LIN and PopPUNK keep
    identifiers for existing groups when genomes are added.
17. **No node annotation and no format specification.** Nothing reports a
    node's size, depth, support or divergence; there is no rank column;
    `docs/output.md` describes the two files in one line each, and
    `docs/developing.md` still says tree2tax roots every tree.

## Improvement ladder

Ordered by value over effort within each batch. S, M and L are rough effort
classes.

### Batch A: small, independent

| Id | Change | Gap | Effort |
|---|---|---|---|
| A1 | `_validate_phylo_request()` called from `build_tree`, `build_msa` and `phylo-build`: refuse split MSA with a genome-input builder and name the flag; `auto` counts the outgroup and respects the MSA source; warn on bootstrap 0 with an MSA builder; consistent Nextflow default | 6 | S |
| A2 | IQ-TREE `model=` extra passed as `-m`; the implicit ModelFinder default logged; `model=GTR+G` documented as the fast recipe | 2 | S |
| A3 | Injective label sanitising with collision suffixes and a forward map for `restore_leaf_names` | 7 | S |
| A4 | A unit test mirroring the FlexTaxD master reader (header rule, non-empty rank, headerless map, label-keyed nodes, `root`); a warning when `--root-name` is not `root`; a format block in `docs/output.md` naming the minimum FlexTaxD commit and the 0.1.5 first-row issue; doctor checks (one `root` without a parent, one parent per child, no empty field); a live test behind `requires_binary("flextaxd")` | 11 | S |
| A5 | The four stale statements and the rooting text in `docs/developing.md` | 10, 17 | S |

### Batch B: tree2tax core

| Id | Change | Gap | Effort |
|---|---|---|---|
| B1 | One post-order pass with a commutative 128-bit set hash per node; edges emitted once per non-root node; no empty labels; a 2000-leaf ladder timing test | 12 | S-M |
| B2 | Prune the outgroup leaf and suppress unifurcations before collapse under `--remove-outgroup` | 13 | S |
| B3 | `--segments` on `tree2tax-relations`; the Nextflow viral path emits and stages `segments.tsv`; the three missing inputs join the resume digest; `node_basename` recorded | 14 | S |
| B4 | A warning at `run` argument time when `--collapse-support` is combined with mashtree or sourmash | 15 | S |

### Batch C: phylo methodology (each bumps the MSA stamp version)

| Id | Change | Gap | Effort |
|---|---|---|---|
| C1 | Count invariant A, C, G, T columns in the existing scan; carry them in `SnpResult`, a sidecar next to the alignment and the stamp; `TreeParams.constant_sites`; IQ-TREE `-fconst`, RAxML-NG `ASC_STAM`; FastTree logs once; `asc=lewis` opt-in only | 1 | M |
| C2 | `--core-fraction`: a column is kept when the called-base fraction meets it; per-genome missing fraction logged; the same rule feeds the C1 counts. This round keeps the current behaviour as the default and logs what 0.95 would retain; the default changes in a later release | 3 | M |
| C3 | `--reference medoid\|first\|NAME`; the medoid from the workdir sketches with a deterministic tie-break; fallback to `first` with the present warning when sketches are absent; choice and score in the record and the stamp key | 4 | M |
| C4 | Gubbins `--first-tree-builder rapidnj --first-model JC` by default, through the existing `gubbins_first_tree_builder=` extra, which still overrides it; masking refused above a named variable-site ceiling or when more than one species is present, with an override; Gubbins' own tree kept as `masker_tree.nwk` beside the masked alignment | 5 | S-M |

### Batch D: taxonomy quality

| Id | Change | Gap | Effort |
|---|---|---|---|
| D1 | `--collapse-diameter D`: maximal clades whose patristic diameter is at most D are collapsed to one node, applied before the support and length rules; the unit of D depends on the builder (mash distance, Jaccard distance, substitutions per SNP site) and is documented | 15 | M |
| D2 | `tree2tax_nodes.tsv`: node, parent, leaf count, depth, support, branch length, relative evolutionary divergence (`NA` without lengths or a meaningful root); GTDB lineage per leaf; optional in doctor | 17 | S-M |
| D3 | `--ranks`: a third column, non-empty on every row, with user-supplied divergence bands; off by default | 11, 17 | M |
| D4 | `--naming lca`: a node is named by its two child-minimum leaf labels, so existing names survive the addition of later accessions; the leaf-set hash remains the default for one release | 16 | M |

### Batch E: scale envelope (E2 to E4 shipped; E1 deferred, needs cluster time)

| Id | Change | Gap | Effort |
|---|---|---|---|
| E1 | A phylo benchmark sweep (genome count by builder by MSA source) and reconciliation of every declared limit; mashtree's declared 10 000 reconciled with the measured 5000; SibeliaZ marked unsupported at scale | 8 | M |
| E2 | A VeryFastTree adapter modelled on FastTree, with `-gamma` as an extra for both | 8 | S |
| E3 | rapidNJ or decenttree behind the Python neighbour joining, used when the binary is present | 8 | M |
| E4 | `ska distance` as a distance source, `lo=true` as a typer extra, the divergence warning on the ska2 path; documented as within-species only | 9 | M |

## Not pursued, and why

- GTDB-Tk's de novo workflow inside `phylo`: it answers a different question
  (placement on the reference tree) and needs its database. The pangenome
  route, first listed here as not pursued, shipped later as the `annotate`
  stage (bakta, prokka) and `phylo --msa-source coregene` (Panaroo), since a
  genus- or family-level tree needs it.
- `+ASC_LEWIS` as a default: IQ-TREE aborts under it when any invariant column
  remains, which masking can produce. `-fconst` uses counts the pipeline
  already has and cannot abort on a remaining invariant column. An external
  `snp-sites -C` pass is unnecessary for the same reason.
- Further work on SibeliaZ or cactus at scale: the measurements do not support
  their declared limits, and the engineering would be large. Correcting the
  limits is the useful step.
- ska2 or Gubbins as genus-level defaults: both are within-species tools, and
  the verification record shows why.
- Incremental placement (EPA-ng, UShER): the right long-term route for adding
  genomes to an existing tree, but it needs a reference-tree artifact, a
  placement contract and a new command.
- A Nextflow scatter of phylo: tree inference is global over the genome set;
  the MSA and tree split is the right parallel boundary.
- A different header for `tree2tax.tsv` or a header in `genomes_map.tsv`: the
  FlexTaxD master reader accepts the current forms, and a header row in the
  map would be imported as a genome.
- GTDB taxon strings as internal node names: FlexTaxD merges identical labels,
  and a GTDB species is an ANI set that need not be monophyletic in the tree.
  The lineage belongs in the node sidecar.
- Automatic calibration of divergence to ranks: GTDB derives its rank
  intervals from a tree of thousands of taxa; a genus-level tree gives no
  distribution to calibrate on. User-supplied bands are honest about that.
- FlexTaxD in the `dev` extra or on CI: the mirror reader in A4 covers the
  contract, and a live test covers the real import.
- A second rooting policy in tree2tax: since #263 phylo owns rooting.

## Verification

- CI with fake tools: unit tests for the column scan (invariant counting,
  thresholds), label collisions, the collapse rules (diameter, divergence,
  naming), the FlexTaxD mirror reader, the tree2tax step (segments, outgroup
  pruning under a kept root, root-name warning), the phylo argument assertions
  (`-fconst`, `ASC_STAM`, `-m`, validation errors), a 2000-leaf ladder timing
  test and the nf-test stub of the viral path.
- The four-genome Francisella hispaniensis set on the development machine:
  `phylo` with mashtree and with IQ-TREE and `-fconst`, then `tree2tax` with
  `--collapse-diameter`, `--ranks` and `--remove-outgroup`, followed by a
  FlexTaxD import and dump in a fresh environment.
- The cluster: the Gubbins defaults, the choice of the soft-core default, the
  batch E sweep and the mashtree ceiling.
