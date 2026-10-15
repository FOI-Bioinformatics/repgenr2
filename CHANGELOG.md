# Changelog

All notable changes to RepGenR2 are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `tree2tax --naming lca` (also on `tree2tax-relations` and `run`) names an
  internal node by a hash of the smallest leaf label under each of its
  children, the set whose last common ancestor it is, so a genome added with
  a later-sorting accession that splits a branch leaves existing names
  unchanged; one joining a multifurcation, sorting below a child's minimum,
  or being removed renames the nodes concerned. `leafset`, the default,
  keeps the hash of the whole leaf set. An unknown rule is refused before any
  stage runs; the record carries the rule, or null under `--node-basename`.
- `tree2tax --collapse-diameter D` (also on `tree2tax-relations` and `run`):
  the fewest clades whose leaf-to-leaf distance stays within D are kept and
  every node inside each is merged into it, so a distance cutoff such as a
  mash distance of 0.05 becomes one FlexTaxD node per cluster, with no
  supports needed. It runs before `--collapse-support` and
  `--collapse-length` and refuses a tree without branch lengths; the record
  carries the value and `collapsed_nodes` counts these collapses too; an
  existing workdir reruns `tree2tax` once.
- `tree2tax --ranks` (also on `tree2tax-relations` and `run`) adds a third
  column to `tree2tax.tsv`, `strain` on every genome and `no rank` on every
  internal node, filled on every row as FlexTaxD's reader requires;
  `--rank-red species=0.95,genus=0.8` names internal nodes by relative
  evolutionary divergence instead, each rank once per lineage in threshold
  order, and is refused on a tree without branch lengths or a root. Off by
  default; the record carries both values. The FlexTaxD contract test reads
  the column and rejects an empty rank.

### Documentation
- `docs/audit/phylo-tree2tax-gap-2026-10.md`: gap analysis of the phylo and
  tree2tax stages against current practice, with the improvement ladder
  that later entries implement (batches A to D; E deferred).

### Changed
- `--reference` on `phylo`, `phylo-build`, `snptype` and `run` takes a
  genome filename, `medoid` or `first`. The default is now `medoid`: the
  genome with the largest summed sketch ANI to the others, computed from the
  workdir sketches (`--sketch`, `repgenr sketch`) with sourmash, ties by
  filename; only sketches that exist are read, none are written for this.
  Without a sketch for every genome, or without sourmash, the first genome
  is used as before, with a warning. The medoid is a less arbitrary choice
  than the alphabetically first genome and is expected to reduce the
  reference bias (swapping the reference changed 72 of 94 splits on a
  50-genome set); it does not remove it. The chosen genome, the rule and the
  medoid's mean ANI are kept in `reference.txt` beside the alignment and in
  the stage record as `reference_genome` and `reference_rule`. The alignment
  cache key holds the chosen genome, so a cached alignment is reused only
  for the same reference; an existing workdir whose sketches now give a
  medoid realigns once. Under Nextflow, `PHYLO_MSA` now takes the staged
  sketches directory like `PHYLO` and `PHYLO_TREE` (bacterial path with
  `params.sketch`); the viral and reads paths have no sketches and use the
  first genome.
- `--mask gubbins`: unless `gubbins_first_tree_builder` is given, the first
  iteration uses rapidnj with the JC model when rapidnj is available, and the
  resolved builder thereafter; Gubbins' own final tree is kept as
  `masker_tree.nwk` beside the masked alignment. Masking is refused when more
  than 12% of the columns are variable (`gubbins_divergence_limit=FRACTION`
  moves the limit) or when the genome names span more than one species;
  `gubbins_force=true` runs it anyway, with a warning. The species check
  reads only names whose last token is an accession, so a user's own
  four-token names are not taken for taxonomy, and it runs before SNP typing
  together with the check of `gubbins_divergence_limit`. Alignments shorter
  than 1000 columns are not judged.
  `gubbins_force=true` runs it anyway. Alignments shorter than 1000 columns
  are not judged.
- `tree2tax` writes `tree2tax_nodes.tsv` beside its two tables: one row per
  node with its parent, leaf count, depth, branch length, support (as a
  fraction), relative evolutionary divergence (Parks et al. 2018, computed as
  PhyloRank does; measured from the ingroup's common ancestor when the
  outgroup is removed; NA without branch lengths or a meaningful root) and,
  for a genome, its family, genus and species from the manifest. FlexTaxD
  does not read it; it is for reading the taxonomy and choosing collapse
  thresholds. The stateless step and the Nextflow module write and publish
  it too, without lineages. An existing workdir reruns `tree2tax` once, since
  the file joins its deliverables.
- `snptype` and the typing pass of `phylo --msa-source snptype` count the
  constant A, C, G and T columns of the whole-genome alignment (`simple`,
  `snippy`, `parsnp`) into `alignment_stats.json` beside `core_snp.fasta`.
  `phylo` passes the counts to the tree builder as ascertainment-bias
  correction: IQ-TREE `-fconst`, RAxML-NG `+ASC_STAM{...}`; invariant columns
  that masking left are removed from the alignment the builder sees and join
  the counts. `--aligner-arg asc=lewis` selects the Lewis correction instead.
  FastTree applies none and says so; `ska2` yields no counts and `phylo` warns
  that its branch lengths are uncorrected; a `model=` that already carries an
  ASC term is passed unchanged. The sidecar names its alignment and holds its
  digest, so the counts never apply to another file. An alignment typed
  earlier is typed once more to gain the counts; aligner alignments are
  reused as before. The Nextflow split mode carries the sidecar from the
  alignment task to the tree task.
  that its branch lengths are uncorrected. The alignment reuse stamp is at
  version 4, so an alignment typed earlier is typed once more to gain the
  counts.
- `tree2tax-relations --segments`: the stateless step lists the member
  segment accessions of grouped viral isolates in `genomes_map.tsv`, as the
  workdir stage does; the Nextflow viral path stages `segments.tsv` from the
  acquisition task when `vgenome --group-segments` wrote it. The resume
  fingerprint of `tree2tax` digests `segments.tsv`, `outgroup_accession.txt`
  and `outgroup/` (an existing workdir therefore reruns `tree2tax` once), and
  the stage record carries `node_basename`.
  and `outgroup/` when present, and the stage record carries
  `node_basename`.
- `tree2tax`: internal nodes are named in one pass over the tree (the sum of
  the descendant leaves' label hashes), so a ladder-shaped tree from clonal
  data no longer costs time proportional to the square of its size; the name
  still depends on the leaf set only, but its value differs from the earlier
  md5 of the joined labels. Edges are written once per node. With
  `--remove-outgroup` the outgroup leaf is pruned before collapse and naming,
  so no clade is named for a leaf set that held it and no single-child node
  remains where it sat on a midpoint-rooted tree. The rows of `tree2tax.tsv`
  are now in preorder, a parent before its children; nothing reads the order.
- `run`: with `--collapse-support` and a tree builder that ignores the
  bootstrap (`mashtree`, `sourmash`) a warning is logged before any stage
  runs, since those trees carry no support values.
- `phylo`: `--treebuilder auto` counts the outgroup as a leaf and, when
  `--msa-source snptype` or `--mask` is given, chooses among the builders that
  consume an alignment instead of dropping the request with a warning; the
  same applies when the alignment is built on its own (`phylo-build
  --msa-only`) or given (`--msa`). The log notes an IQ-TREE run without
  bootstrap replicates, since that tree carries no support values.
  `phylo-build --msa-only` with an alignment-free builder names the builders
  that consume one. The Nextflow pipeline refuses `phylo_split_msa`
  with an alignment-free builder before any stage runs.
- `phylo --treebuilder iqtree`: `--aligner-arg model=...` is passed as `-m`;
  without it the log says that IQ-TREE runs ModelFinder first.
- The sourmash tree builder writes a leaf label holding a reserved Newick
  character in quotes instead of replacing the character, so two labels can no
  longer merge into one leaf.
- Documentation: the simple typer's N masking (#232) and concurrency (#166)
  are stated where the audit records still said otherwise; the SNP-phylogeny
  example no longer runs `snptype` first, since `phylo` types the genomes
  itself; the rooting description follows `phylo --root`.
- `tree2tax`: the FlexTaxD contract is stated in `docs/output.md` and held
  by a test that mirrors FlexTaxD's reader (header rule, label-keyed nodes,
  a top node named `root`, no empty field, no header on the genome map).
  `--root-name` other than `root` warns, since FlexTaxD validates its tree
  from a node of that name. `doctor` also fails a taxonomy with an empty
  field, a child under two parents or more than one top node, and warns when
  the top node is not `root`.
- The project is renamed RepGenR2 and lives at
  <https://github.com/FOI-Bioinformatics/repgenr2>. The distribution name is
  `repgenr2`; the `repgenr` command, the `repgenr` import package, the plugin
  entry-point groups (`repgenr.dereplicators` and the others), the workdir
  layout and the conda environment names are unchanged, so existing workdirs
  and plugins keep working. An environment that holds the old `repgenr`
  distribution needs `pip uninstall repgenr` before `pip install .`: the two
  distributions own the same files, `__version__` reads the `repgenr2`
  metadata and reports `0.0.0+unknown` without it, and a later
  `pip uninstall repgenr` would remove files the new install needs. The
  Nextflow manifest name is `repgenr2`; the manifest home page, the schema
  `$id`, the citation and the new `[project.urls]` refer to the new
  repository. User-visible messages name RepGenR2 (for example "No RepGenR2
  run found"), which affects scripts that match on that output.
- `dereplicate --process-size`: before a further level over the union of
  chunk representatives, one chunk of that level is run as a probe. A level
  costs the tool's per-genome work for every representative plus about N/U
  of the final pass's pairwise work (N the chunk size, U the union) and is
  expected to pay when it removes more than about N/U of its genomes; on a
  reported 6888-genome run the levels after the first removed 5, 1.4 and 0.5
  percent, each as long as a level-0 wave, while clusters spread over every
  chunk collapse almost entirely at the next level. When the probe removes
  less than that fraction the union goes to the final pass and the probe's
  result is dropped; otherwise the level runs with the probe as its first
  chunk. The shrink rule (a chunk's worth or a tenth removed) still decides
  whether to probe at all. The log now times each level and the probe,
  states the final pass's pairwise work as a multiple of the level just
  timed, that the tool prints no progress while it runs, and its duration
  afterwards.
- skDER: at an ANI cutoff of 99 percent or more the adapter passes skani a
  screen of the cutoff minus one (`--skani-triangle-parameters="-s PERCENT"`
  to skDER), so pairs whose k-mer estimate is below it are not aligned.
  skDER's own default is the cutoff minus ten, which at 99.95 aligned every
  pair above about 90 percent, most same-species pairs in a genus-level set. `--tool-arg screen=PERCENT` sets
  another value; a value above the final cutoff is refused before any work,
  and a call below the final cutoff (a `--pre-secondary-ani` chunk, a
  `--target-reps` step) keeps the default rule for that call with a warning.
  `screen=skder` restores skDER's own default, which is also kept below a 99
  percent cutoff. The log names the screen and its source once per run.

### Added
- `assemble --screen-reads` (off by default) screens each run's reads sketch
  after the fetch and before an assembler starts, and excuses a run with step
  `screen` when its GTDB genus differs from the submitted one
  (`taxon_mismatch`; `sourmash gather` and `tax genome` as the classifier
  runs them, with its `genus_renamed` tolerance), when the top species holds
  less than `--screen-min-fraction` (0.3) of the abundance-weighted reads
  (`low_match_fraction`), or when an accepted run of the same biosample
  earlier in `reads.tsv` contains it at an ANI estimate of `--screen-dup-ani`
  (0.999) or more (`duplicate_isolate`). A tool failure is `screen_failed`,
  which is not counted as judged. Screened-out runs keep their reads sketch
  and gather output; each decision is stored in `assemblies/<run>/screen.json`
  with its settings, so a repeat screens nothing again and a changed
  threshold screens again from the kept sketch without a download. The
  working directory gains `screen_reads.tsv`, `assembly_stats.tsv` a `screen`
  column (older tables read as empty), the stage record `n_screened_out`.
  The flag needs sourmash and the GTDB sketch and lineages, and exits 2 (4
  without sourmash) before any download otherwise. `assemble-run` takes
  `--screen-reads`, `--screen-min-fraction`, `--gtdb-sketch` and
  `--gtdb-lineages` without the duplicate check; in Nextflow
  `--assemble_args '--screen-reads'` adds the pipeline's GTDB paths, and
  `reads-gather` writes `screen_reads.tsv`. The classifier's gather and
  `tax genome` calls are shared functions now.
- `ingest --from-workdir` checks the record of the stage that wrote each
  source's genome set (the latest of `genome`, `vgenome`, `ingest` and
  `assemble`) with the evaluation `status` applies in that working directory,
  and logs one warning per source whose record is interrupted or stale, or
  that has no record ("no stage record"); the genomes are taken as they are.
  The new `--strict-sources` refuses such a source with exit 2 before
  anything is written. The ingest record keeps `source_states` (path, stage,
  state per source), and the `status` line of `ingest` names the sources that
  were not done when they were ingested.
- `assemble --max-runs N` keeps at most N runs with a finished assembly in
  the working directory, taken in `reads.tsv` order (largest by bases first).
  Finished runs count toward N and are always kept; runs excused before any
  download take no place. The remaining runs are written to
  `excused_runs.tsv` with reason `deferred`, so the completeness guard lets
  `dereplicate` run on the assembled part, and a later call with a larger N,
  or without the flag, assembles them and removes their rows. Deferred runs
  are not counted as judged when every run is excused, N is part of the
  resume fingerprint, `--append` takes it, and `status` shows
  `n assembled, n deferred`. The stage record gains `n_deferred`, and
  `n_excused` no longer counts deferred runs. The Nextflow layer bounds the
  selection instead (`--reads_args '--max-runs N'`).
- `phylo --root outgroup|midpoint|none` (also on `phylo-build` and `run`):
  phylo roots `tree/tree.nwk` after any tree builder, with dendropy. The
  default is `outgroup` when an outgroup is staged and `--no-outgroup` is not
  given, else `none`. `outgroup` roots on the branch to the outgroup, split in
  two halves, as `tree2tax` did; `midpoint` roots at the midpoint and needs
  branch lengths. A tree already rooted on the outgroup (bifurcating root) is
  left as written; the trifurcating root that `iqtree -o` and
  `raxmlng --outgroup` write is rerooted with one INFO line. Support values
  are carried by split, so each stays with the branch it was computed for.
  The builder's output is kept as `tree/tree.unrooted.nwk` when the tree is
  rerooted. `--root outgroup` without an outgroup, or with `--no-outgroup`,
  exits 2. The phylo record gains `root`. The rooting functions are in
  `repgenr.tree.rooting` and `tree2tax` uses the same one.
- `tree2tax-relations --keep-root` keeps the tree's own root (a midpoint
  root) instead of rerooting on the outgroup; `tree2tax` does so when phylo
  recorded `root: midpoint`. The Nextflow `TREE2TAX` process receives
  `--keep-root` when `phylo_args` contains `--root midpoint`.
- Reads sketch: `assemble` sketches each run's reads with sourmash
  (`sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund` on all FASTQ
  files of the run, named by run accession) into
  `assemblies/<run>/reads.sig.zip`. The sketch runs in a thread beside the
  assembler, which keeps its full thread share (the single-threaded sketch
  runs one thread over it), and it ends before the reads are deleted. `--reads-sketch`
  requires sourmash (exit 4 before any download), `--no-reads-sketch` skips
  it, and the default sketches when sourmash can run. A failed sketch is a
  warning and never fails the assembly. The done marker records
  `reads_sketch`, `reads_sketch_params` and `reads_sketch_version` (null
  without a sketch), `assembly_stats.tsv` gains a `reads_sketch` column (1 or
  0; older tables still read), and excused runs keep no sketch. A finished
  run without a sketch is not assembled again; its kept FASTQ files
  (`--keep-reads`) are sketched on the next call. The reads sketch is not a
  genome sketch and is never written to `sketches/`. `assemble-run` takes the
  same flags.
- Downloads check their md5 while the bytes are written: `core.http.download`
  takes an optional `md5` and compares the digest at the end, deleting the
  `.part` file on a mismatch with the error `verify_md5` raises. The
  `assemble` fetch passes ENA's checksum and no longer reads each FASTQ file a
  second time after a fresh download; a file kept from an earlier attempt and
  a local copy are still hashed with `verify_md5`.
- `repgenr census` counts the genera, species and samples under a taxon,
  read-only. Without a working directory it counts the GTDB genomes and
  species representatives of a family (`-tf`, one row per genus) or genus
  (`-tg`, one row per species) through the GTDB API or a cached GTDB metadata
  table (`--source table -r 232.0`), with the ENA whole-genome sequencing runs
  under the taxon on request (`--runs`: runs, biosamples, runs per platform),
  or the NCBI Virus records of a viral taxon (`--viral --target`, metadata
  report only): sequences, complete sequences, isolates after segment grouping
  and whether the virus is segmented. With `-wd` after metadata or vmetadata
  it counts the candidates the entry stage found (GTDB table or API answer,
  NCBI Virus records, BV-BRC taxonomy sets); after genome, vgenome, ingest or
  assemble it counts the selected genomes per manifest source, with the
  candidates and the dereplication clusters beside them. `--tsv` writes the
  rows and `--json` prints one `repgenr.census/1` object. `metadata --source
  api` now keeps the API answer for the target in `gtdb_api_genomes.tsv`.
- Genome sketches, step 2: the sourmash tools read `sketches/` instead of
  sketching each genome again. `dereplicate --tool sourmash` (branchwater
  `pairwise` over a path list of the `.sig.zip` files, or `sourmash compare`),
  `glance --tool sourmash`, `phylo --treebuilder sourmash` (outgroup included)
  and the sourmash classifier of `assemble` select their k-mer size with
  `-k`; k=21, 31 or 51 at scaled=1000 come from the sketches, and other
  `ksize` or `scaled` values fall back to the tool's own sketch with a log
  line naming the reason. A consumer first writes the missing and stale
  sketches of the genomes it compares and logs `sketches: n reused, m
  written`. Whether a sketch is current still rests on the FASTA SHA-256,
  now kept in `sketches/.digests.json` with the file's size and modification
  time, so an unchanged genome is not read again (`repgenr --force sketch`
  hashes every genome). On 1000 genomes `dereplicate --tool sourmash` took
  about 6 s instead of 65 s and `glance` about 20 s instead of 104 s.
  `assemble` sketches each accepted assembly once
  (`assemblies/<run>/contigs.sig.zip`, reused while the contigs are
  unchanged); the classifier gathers with it and the sketch step copies it to
  `sketches/` under the genome's record name (`sourmash sig rename`).
  Adapters declare the parameters they compare at with `sketch_request()`
  and receive the files in `DerepParams.sketches`, the `sketches` argument of
  `compare`, `TreeParams.sketches` or `ClassifyParams.sketches`. The stateless
  steps `dereplicate-chunk`, `dereplicate-merge` and `phylo-build` take
  `--sketches-dir`. Nextflow: `--sketch` (default false) runs `SKETCH` on the
  bacterial path, publishes `sketches/` and stages it into the dereplication
  and tree tasks; `SKETCH` now reads its genome list with a `read` loop
  (file names with spaces stay whole, user arguments are not nested in
  quotes) and fails when any sourmash call fails.
- Genome sketches as a working-directory contract (step 1 of 2):
  `genome`, `vgenome`, `ingest` and `assemble` write
  `sketches/<name>.sig.zip` per genome, outgroup included, with DNA
  signatures at k=21, 31 and 51 (scaled=1000) named after the genome's record
  name, when sourmash can run. `--sketch` requires sourmash (exit 4 before any
  genome or sketch is written) and `--no-sketch` leaves the step out; the flag does
  not change the resume fingerprint. Only missing and stale sketches are
  written (stale: the FASTA SHA-256 or the parameters differ from the
  record), `assemble --append` sketches only the genomes it adds,
  `ingest --from-workdir` copies a matching source sketch, and sketches are
  pruned with the genome set. The new command `repgenr sketch -wd WD [-t N]`
  sketches an existing working directory and logs the counts; `status` shows
  `sketches: n/m`. The manifest moves to schema version 4 with the columns
  `sketch_file`, `sketch_params` and `sketch_digest`; a version 3 manifest is
  migrated when a stage opens it, and its digest for the resume fingerprint is
  unchanged. A Nextflow `SKETCH` module (stub-tested, not yet wired into the
  workflows) sketches a genome directory.
- `ingest --from-workdir WD` (repeatable, alone or with `--genomes-dir`;
  also on `run`) takes the genome set of an earlier working directory: its
  `selection.tsv` rows, with taxonomy, CheckM quality and
  `gtdb_representative` unchanged, and the files under its `genomes/`. The
  manifest keeps each genome's source (`gtdb`, `sra`, ...) from the source
  manifest, or `local` without one, so a GTDB working directory and a
  reads/assemble working directory can be combined and dereplicated together
  with `--keeper gtdb`. The outgroup of a source is not carried over (one log
  line per source); `--outgroup` may name a genome of any source. A genome in
  two sources, a missing file, or the target working directory as a source
  exits 2. The selection and `genomes/` of each source are resume inputs; the
  ingest record lists `from_workdirs` and the genome count per source, and
  `status` shows the source working directories. `--genomes-dir` is no longer
  required, and `--selection` without it exits 2. The new parameter changes
  the resume fingerprint of `ingest`, so an existing ingest record reruns
  once.
- `--keeper gtdb` (#252) for `dereplicate`, `run`, `bacterial`,
  `dereplicate-chunk` and `dereplicate-merge` (Nextflow `--derep_keeper gtdb`):
  within each cluster a GTDB species representative is kept; a cluster without
  one falls back to the quality rule, then to the tool's pick. A cluster that
  holds several GTDB representatives (species joined by the ANI threshold)
  keeps the best-scored one and a warning names the others. `--reduce` also
  prefers a GTDB representative under this rule. `repgenr.yaml` records
  `keeper_effective: gtdb` when the manifest flags at least one genome.
- `selection.tsv` gains a last column `gtdb_representative` (#252), filled by
  `metadata` from the GTDB table (`gtdb_genome_representative`) and the API
  (`gtdbIsRep`), 0 on the other entry paths; it is read with the same 1/0,
  true/false, yes/no parser as `is_outgroup`, and a file without it reads as
  0. `derep/cluster_summary.tsv` gains a last column
  `rep_is_gtdb_representative`.
- `list-tools --check` under a container backend (#248): each line names
  where the tool runs, `[image <ref>]` or `[host]`, so a host version is not
  read as the image's. `--images` adds whether each image of a tool that passed
  the check, secondary images included, is present locally (`docker image
  inspect`, or the Singularity `.sif` cache; nothing is pulled), or `presence
  unknown` when the engine does not say the image is absent.
- Tool versions (#248): a stage that runs a tool in an image records the
  engine and its version beside the image reference (`docker: 29.5.3`), in
  `repgenr.yaml`, `repgenr versions` and the Nextflow `versions.yml`. The
  resume fingerprint does not include versions, so no stage reruns for this.
- `list-tools --check --strict` (#237): after the full listing, exits 4
  when any adapter is missing or errored and 5 when any plugin failed to
  load, so a script can verify an environment. `--check` alone still exits 0,
  since a host with only some families installed is normal. `--strict`
  without `--check` is a usage error.
- Tool versions (#237): a binary that reports no version number itself
  (SibeliaZ has no version flag) is recorded with the version of the conda
  package that installed it, read from `conda-meta` next to the binary's
  `bin` directory; it was recorded as `unknown`. A stage that recorded
  `unknown` is not rerun for this, since versions are not part of the resume
  fingerprint.
- `status --json` and `doctor --json` (#242): one versioned JSON object on
  stdout (`repgenr.status/1`, `repgenr.doctor/1`) in place of the text report,
  so scripts and workflow wrappers need not parse wording; exit codes are
  unchanged and stdout stays empty on exit 3. The schemas are in output.md,
  "Machine-readable status". The text output of both commands is unchanged.
- `doctor --quick` (#242): skips reading the first bytes of each genome file
  (about 30 s per 1000 genomes on an exFAT disk); links, missing and untracked
  genomes are still checked, and the genome line says the content was not
  read.
- Global `--bin-dir TOOL=DIR` (#247), repeatable, and `REPGENR_BIN_DIRS`
  (`tool=dir,tool=dir`): DIR comes first on `PATH` for that tool only, at
  preflight (lookup and version query) and in every host subprocess it runs,
  so a satellite environment's tool uses its own helpers and never shadows
  core's. Gubbins' tree-builder choice, the FastTree binary name and the
  SibeliaZ wrapper lookup use the same per-tool path, and a tool without a
  version flag is read from `conda-meta` of the environment it was found in.
  An unknown tool or a missing directory exits 2, and the error names
  `--bin-dir` or `REPGENR_BIN_DIRS`, whichever held the entry; under `--container` a tool that runs in an image
  is named in a warning. The directory each tool used is recorded under
  `bin_dirs` in the stage record and is not part of the resume fingerprint, so
  existing workdirs do not rerun. The live suite passes its `[bin_dirs]`
  entries through `REPGENR_BIN_DIRS`; a key that is not a tool name (such as
  `nextflow`) is appended to `PATH`.
- Install (#238): `envs/` holds seven conda environment files that solve: a
  core environment (`envs/core.yml`, Python 3.12, RepGenR and every tool that
  shares its solve) and six satellites for tools that do not (Gubbins,
  mashtree, snippy, parsnp with harvesttools, CheckM2, progressiveMauve). Each
  file's first line names the platforms that solve and why it is separate.
  `.github/workflows/envs.yml` dry-runs each solve on linux-64, and core on
  osx-arm64, for changes to `envs/` and weekly. `tests/test_env_files.py`
  checks that every preflight version floor appears in one file with the
  same value, that core pins Python 3.12, and that no package is in two files.
- `status` (#234): a finished stage whose recorded input changed, or whose
  declared output is missing, since it finished is listed as `[stale]` with
  the reason; it re-runs on its next invocation. `Next:` names the first
  stage that is not done (stale, interrupted or not run), and when that is
  `phylo` with fewer than three representatives a note says phylo will refuse
  and names `--all-genomes` and `--secondary-ani`. status and doctor use the
  same checks (`core.doctor.stale_stages`), so they agree; derep-stock
  records are exempt, since a later dereplicate changes their inputs by design.
- `doctor` (#234): checks that `tree2tax.tsv` has its header and names the
  same leaves as `genomes_map.tsv`, warns about genome files under `genomes/`
  that `selection.tsv` does not list, and warns when `repgenr.yaml` records
  no stage beside existing outputs. Its help states the exit codes: 0 with
  only warnings, 1 on failures (7 since #242), 3 for a missing workdir.
- `ingest` (#229): genome files ending in `.fna.gz` (the NCBI FTP default) and
  `.fa.gz` are accepted beside `.fasta.gz` (one suffix list in
  `core.contracts.FASTA_SUFFIXES`), and were previously skipped with a
  warning. They are staged compressed, like `.fasta.gz`; `derep-unpack` names
  their cluster directories without the whole suffix. An earlier `ingest`
  that skipped such files is not rerun on its own, since its source directory
  is unchanged; `repgenr --force ingest` stages them.
- `metadata --source api` (#229): the stage record holds `api_query_date`, the
  UTC time of the GTDB API query, since the API reports no release number
  (`release` stays null). `status` shows the GTDB release or the query date
  after the metadata line, and `versions` prints `gtdb_release` or
  `gtdb_api_query_date`, so the Nextflow `versions.yml` of the metadata step
  carries it too.
- `assemble` stores each run's CheckM2 scores in `assemblies/<run>/checkm2.json`,
  keyed by the contigs' SHA-256, the database's resolved path, size and
  modification time, and the CheckM2 version (or image). `--force` does not
  bypass them; deleting the file scores the run again. A call that changes only `--min-completeness` or
  `--max-contamination` applies the stored scores instead of running CheckM2
  again (about 5 minutes for two genomes under emulation). The first CheckM2
  pass after upgrading stores the scores; no stage reruns because of this
  change. `genome-qc` stores nothing, since Nextflow stages its inputs (#230).
- `assembly_stats.tsv` flags a genome `genus_renamed` when the GTDB genus
  differs from the submitted one but the family and the species epithet agree
  (NCBI Mycoplasmopsis arginini, GTDB Metamycoplasma arginini). The submitted
  name stays, the genome is not excused, a warning is logged, and the stage
  record counts these as `n_genus_renamed` apart from `n_disagree`. A shared
  epithet in another family stays `classifier_disagrees` (#230).
- `genome-qc --memory-gb` (default 16; the Nextflow module passes the task
  memory). With `assemble --memory-gb`, it bounds the number of concurrent
  sourmash gathers at about 0.6 GB each, as well as the thread count; the log
  names the number chosen (#230).
- `glance --tool sourmash`: sourmash is a second comparison backend for
  `glance`. It sketches every genome with the parameters `dereplicate --tool
  sourmash` uses (k=31, scaled=1000), runs `sourmash compare`, and converts the
  values to the ANI estimate the dereplication threshold applies to, so the
  plots show the scale of that threshold. The dendrogram is average-linkage
  clustering on 1 - ANI, computed with numpy (no new dependency). With
  `--keep-files`, `glance_wd/` holds `pairwise_ani.csv` and the dendrogram's
  leaf order. Sets above 5000 genomes are refused (dense matrix). The plot
  axes name the measure, Mash ANI for dRep and ANI for sourmash; the output
  file names are unchanged. `list-tools` ends with a line naming the glance
  backends.
- `cluster_summary.tsv` gains a last column, `n_genomes`, which counts the
  genomes in the cluster including the keeper (`n_members` excludes it).
- `repgenr --help` groups the commands into panels (pipeline, entry points,
  core stages, inspection, environment, Nextflow steps) in pipeline order, and
  its epilog lists the four stage chains (bacterial, viral, local genomes,
  sequencing reads).
- `list-tools` prints each adapter's declared genome limit, for example
  `iqtree (up to 500 genomes)`; a tool without a declared limit shows its name
  only.
- An eighth tool family, polishers (`repgenr.polishers`), corrects long-read
  assemblies with the run's reads inside `assemble` and `assemble-run`:
  medaka for ONT (basecaller model from Dorado read headers, `--tool-arg
  model=...`, or ONT's bacterial methylation model assumed for SRA-mirrored
  reads whose headers name none; recorded in the marker) and racon with minimap2 overlaps for PacBio CLR
  (`--polish-rounds`). `--polisher auto|medaka|racon|none`; HiFi and
  Illumina are never polished. `assembly_stats.tsv` gains a `polisher`
  column and a failed polish excuses the run with `polish_failed`.
- `reads.tsv` records ENA's `library_selection`, and `reads --drop-selection`
  (default `MDA`, repeatable, `none` to keep all) drops amplified libraries,
  which assemble into chimeric and uneven contigs (68 of the 551 Wolbachia
  WGS runs are MDA).
- `reads --max-bases N` drops runs above a size, a guard against unenriched
  whole-host libraries (36 of the 551 Wolbachia runs exceed 5 Gb) that would
  assemble into a host-dominated genome.
- Nextflow `--mode reads`: the reads chain as data-channel tasks. `READS_SELECT`
  runs the reads stage, `READS_ASSEMBLE` assembles one run per task (new
  `process_assembly` label: 8 CPUs, 16 GB for short reads and 32 GB for long
  reads via the run's platform), `GENOME_QC` scores and classifies the batch
  when `--checkm2_db` or `--gtdb_sketch` is set, and `READS_GATHER` writes the
  genome contract that feeds the shared dereplication, phylo and tree2tax
  modules. Parameters `reads_args`, `assembler`, `assemble_args`,
  `checkm2_db`, `gtdb_sketch`, `gtdb_lineages`; tables published under
  `reads/`. Behind the modules, three stateless CLI steps: `assemble-run`,
  `genome-qc` and `reads-gather`. The `assemble` stage now keys its QC and
  classification inputs by run accession rather than by the final genome
  name.

- `assemble --append` adds reads-derived genomes to a working directory that
  already holds a GTDB or local selection, keeping its rows, files and
  outgroup; `metadata --drop-foreign` and `ingest --drop-foreign` are the
  explicit way to discard them, since both stages otherwise refuse to
  overwrite a selection holding appended genomes.

- Assembly quality and organism verification for the reads chain. With a
  CheckM2 database (`--checkm2-db` or `CHECKM2DB`) every assembly is scored,
  the values reach `selection.tsv` and the manifest for the quality keeper,
  and assemblies outside `--min-completeness`/`--max-contamination` are
  excused. A seventh tool family, classifiers (`repgenr.classifiers`), with
  `sourmash` gather against a GTDB sketch (`--gtdb-sketch`,
  `--gtdb-lineages` or the matching environment variables): when the GTDB
  genus agrees with the submitted organism the GTDB tokens name the genome,
  otherwise the submitted name stays and the genome is flagged
  `classifier_disagrees` in `assembly_stats.tsv`. The sketch release is
  recorded in provenance.

- `repgenr assemble`, the second stage of the reads chain: fetches each
  selected run's FASTQ files from ENA over HTTPS with checksum verification,
  assembles them with the adapter that accepts the platform (`--assembler
  auto|skesa|shovill|flye`), filters contigs, names the genome from the
  reads stage's taxonomy and writes `genomes/`, `selection.tsv`, the
  manifest, `assembly_stats.tsv` and `excused_runs.tsv`. Finished runs carry
  a marker and are skipped on a re-run; `--jobs` bounds concurrency. `run
  --reads` runs the whole chain from a taxon or an accession file.

- `repgenr reads`, the first stage of a reads chain (`reads -> assemble ->
  dereplicate -> phylo -> tree2tax`): selects whole-genome sequencing runs
  from ENA (which mirrors SRA) by taxon or by run, sample or study accession,
  filters by platform and size, keeps the best run per sample, labels each
  run with its NCBI family, genus and species, and writes `reads.tsv`.
  `status` recognises the chain. The assemble stage follows.
- A sixth tool family, assemblers (`repgenr.assemblers` entry points), with
  `skesa` and `shovill` for Illumina and `flye` for Oxford Nanopore and
  PacBio reads, each pinning a BioContainers image, plus the shared contig
  filter and summary (contig count, total length, N50, largest contig) the
  assemble stage will use. `list-tools` lists the family.

- `vgenome --outgroup-accession` (and `run --viral --outgroup-accession`)
  pins the viral outgroup to a downloaded record on either back-end; the
  BV-BRC path previously asked the user to "specify one manually" with no
  way to do so. Segment-grouped runs (`--group-segments`) now search for an
  outgroup too, using the kept records' length span as the window; they were
  left unrooted before.

- `repgenr run` starts from local genomes with `--genomes-dir` (the `ingest`
  chain, with `--selection`, `--outgroup` and `--copy` passed through), runs
  the standalone SNP typing stage with `--with-snptype`, and forwards the
  stage options it did not before: `--metadata-path`, `--nodownload`,
  `--keep-files`, `--pre-primary-ani`, `--pre-secondary-ani`,
  `--node-basename`, `--root-name` and `--remove-outgroup`.

- `dereplicate-merge --reduce species|genus` and `--target-reps N`, the two
  representative-selection features that only the workdir `dereplicate`
  command had, so the Nextflow layer can use them through `--derep_reduce`
  and `--derep_target_reps`. The merge step takes the taxonomy from
  `selection.tsv` when one reaches it and from the genome filenames
  otherwise.

- `repgenr list-tools --check` runs every adapter's preflight and prints one
  line per tool: `ok` with the resolved versions, or `missing` with the
  binaries that are absent or below their version floor. An environment can
  be verified before a run, without a working directory.

- `derep/cluster_summary.tsv`: one row per representative with the member
  count, the species the cluster spans and the keeper's CheckM quality against
  its members, including which genome scores best. `dereplicate` writes it
  alongside `clusters.tsv` (the Nextflow chunk and merge steps too), and
  `repgenr cluster-summary -wd <workdir>` regenerates it for an existing
  working directory without rerunning the dereplicator.

### Changed
- `tree/tree.nwk` is rooted on the outgroup by `phylo` when an outgroup is
  staged; it was the builder's output, unrooted for `fasttree`, `mashtree`
  and `sourmash`, and rooted only in `tree2tax`. Without collapsing, the
  `tree2tax` relations are unchanged. The new `root` parameter enters the phylo resume fingerprint, so
  a completed phylo stage runs once more on the next `phylo` or `run`.
- ENA run searches (`reads` by taxon and `census --runs`) count the runs
  with the portal count endpoint first and log the number; above 100000 runs
  a warning names the count, the taxon and the expected memory, and the
  search proceeds. `core.ena.search_runs` reads the portal's TSV format line
  by line instead of one JSON document (the same string values without the
  field names repeated per record, and the body is never held whole), and
  takes a `fields` argument: the census requests the four fields it reads
  (run, sample, taxid and platform) instead of all 17. A TSV line with the
  wrong number of columns, as from a body cut short, is an error.
- The dense sourmash path of `dereplicate` and `glance` orders the
  `sourmash compare` matrix by genome file name before the greedy pick, and
  the sourmash tree builder joins in name order. `sourmash compare` does not
  keep its input order (seen with 4.9.4), so a tie between equally connected
  genomes could pick a different representative from run to run; the pick is
  now the same each run and matches the branchwater path, which already
  took the genomes in name order.
- The sourmash dereplicator, tree builder and classifier take their container
  image, conda spec and binary check from one shared specification
  (`core/sourmash.py`), which the genome sketch step also uses; the pinned
  image is unchanged.
- Quality keeper (#254): the score is `completeness - 5 x contamination +
  0.5 x log10(N50)`, dRep's default weights, with the N50 read once per
  genome from its FASTA (gzip-aware). Equal scores go to higher completeness,
  then lower contamination, then higher N50, then the filename, so skDER,
  sourmash and galah end with the same keepers on the same clusters; the
  tool's pick no longer wins ties. A scored genome replaces an unscored
  representative only when it is high quality (completeness above 90,
  contamination below 5), and an INFO line names the clusters left to the
  tool. The same rule applies to `--reduce` (the largest cluster's
  representative is the default keeper of a taxon) and to `best_member` in
  `cluster_summary.tsv`, which gains `rep_n50` and `best_score` after
  `n_genomes`. `dereplicate-chunk` writes `genome_n50.tsv` for the merge step.
  On 30 GTDB r232 *F. tularensis* genomes the keeper of the largest cluster
  was the 29-contig GCF_016603775.1 (100/0.00) for skDER, sourmash and galah;
  it is now the two-contig GCF_001880245.1 (100/0.01, N50 1.89 Mb, score
  103.09), ahead of the closed GCF_000833375.1 and GCF_000014645.1 (100/0.03,
  102.99). The N50 is read only where two or more scored genomes are
  compared, and for scored representatives in the summary. Representatives of an existing
  workdir change only when dereplicate is rerun with `--force`; the resume
  fingerprint is unchanged. `cluster-summary` and `derep-stock unpack` write
  the new columns.
- Manifest schema version 3 (#252) adds the column `gtdb_representative`. An
  older manifest is upgraded in place when a stage opens it, with every flag
  0, and its resume fingerprints stay valid (the digest adds a field only for
  a set flag). An earlier RepGenR refuses a version 3 manifest. Rerunning
  `metadata --force` in an existing workdir records the flags; `genome` then
  runs again because `selection.tsv` changed (it downloads nothing already
  present), and `dereplicate` runs again because the manifest changed.
- `list-tools`, `versions` (#237): a warning logged while they run (a plugin
  that fails to load) carries the standard timestamp and level instead of
  Python's bare fallback line.
- `status` (#245): a `repgenr.yaml` that records no stage is no longer read
  as the bacterial chain ("Next: repgenr metadata"); `status` names the entry
  stages, and suggests `doctor` when the workdir holds outputs. A record with
  stages but no entry stage (for example only `dereplicate`) is shown as
  "Pipeline: unrecorded entry stage" with the shared stages dereplicate, phylo
  and tree2tax, and `pipeline` is null in `--json`. An `assemble` record
  without `reads` selects the reads chain. A completed record without a resume
  fingerprint (an older version, or `derep-stock --action unpack`) stays
  `[done]` with a note that its next invocation recomputes it (`"fingerprint":
  false` in `--json`). `doctor` warns about every completed record without a
  fingerprint, stale ones included; `status` notes only `[done]` lines, since
  a `[stale]` line already says the stage re-runs.
- `doctor` (#242): exits 7 (`core.errors.DOCTOR_FAILURES_EXIT`) when it finds
  a failure, including a malformed `repgenr.yaml` or a check that could not
  complete; it exited 1, the status of an unexpected error. A script that
  tested for 1 after `doctor` must test for 7. Warnings alone still exit 0.
- CLI (#244): the default of `-t/--threads` is 16, or the CPU limit of the
  process when that is lower: the affinity mask capped by a cgroup v1 or v2
  CPU quota (`docker run --cpus`, a Kubernetes CPU limit), previously not
  seen. A lowered default is reported once; an explicit `-t` above the limit
  is kept with one warning. The same limit now applies to the host CPU count
  that caps Gubbins threads and to the automatic `--num-processes` of
  `dereplicate`. The thread count is not a result parameter, so finished
  stages do not rerun.
- Record names (#241): every aligner, SNP typer and tree builder names a
  genome by `core.contracts.record_name`, its file name without the FASTA
  suffix and `.gz`. progressiveMauve, SibeliaZ, cactus and the sourmash tree
  builder used `Path.stem`, so a gzipped genome `x.fasta.gz` was the record
  `x.fasta` in `align/msa.fasta` and in their intermediates, and the masker
  received both name forms of the outgroup. Plain genomes keep their names.
  An aligner alignment that a previous run left with such names is rebuilt
  once instead of reused. The MSA stamp version is unchanged, so other
  alignments are still reused; the stamp now records the alignment's record
  names, which a reuse compares with the inputs without reading the
  alignment (a stamp written before is completed once from the header
  lines). The Gubbins masker warns when an outgroup name to leave out
  matches no alignment record, since the outgroup is then part of the
  recombination scan.
- Install (#238): `environment.yml` is removed; it did not solve on linux-64
  (Gubbins needs Python 3.8 to 3.10; mashtree and snippy need zlib older than
  1.3). Create `envs/core.yml` and the satellites you need instead. Satellite
  `bin` directories go after core on `PATH`, never before, so that mashtree's
  samtools 0.1.19 and Gubbins' IQ-TREE 2 cannot shadow core's tools; install.md,
  README.md and choosing-tools.md describe this. Two gaps remain under this
  rule: snippy calls core's samtools and bcftools, and Gubbins calls core's
  IQ-TREE 3 instead of its own 2.4 (the Gubbins log then records no IQ-TREE version, "for" in its place);
  `--bin-dir` (#247) closes both. The progressiveMauve adapter's Wave spec
  now names `bioconda::mauvealigner` (with `boost-cpp=1.74.0`), the package
  that provides the binary of its pinned image, instead of `bioconda::mauve`,
  which adds the Java GUI; `envs/mauve.yml` installs the same packages. All
  files list the `nodefaults` channel, since a configured `defaults` channel
  with flexible priority made the core solve run for more than 17 minutes. The live test suite's
  `[bin_dirs]` are now appended to `PATH` instead of prepended, so a local
  live config that relied on a satellite shadowing a legacy samtools in the
  test environment has to run pytest in the core environment.
- `phylo --msa-source snptype` (#228) writes its typing pass under
  `tree/msa/` (alignment, optional tables and the reuse stamp) and its scratch
  under `scratch/phylo_snptype/`; `snp/` and `scratch/snptype/` belong to the
  `snptype` stage alone. Before, the typing pass replaced the `snptype` stage's
  tables in `snp/` while its record stayed, so a repeat `snptype` skipped and
  `doctor` reported nothing. The interim measures of #223 (phylo removing the
  `snptype` record, and `run --with-snptype --msa-source snptype` placing
  `snptype` after `phylo`) are withdrawn: `run` again types before `phylo`.
  A workdir whose alignment stamp is under `snp/` is typed once more, into
  `tree/msa/`; `doctor` warns about the old stamp. The Nextflow `PHYLO_MSA`
  process emits `tree_msa` (was `snp`), `PHYLO` publishes the typing pass
  with `tree/`, and the published alignment moves from `phylo/snp/` to
  `phylo/tree/msa/`.
- `genome`, `vmetadata` (#229): before the first `datasets` call, one request
  with the 15 s connect timeout checks that `api.ncbi.nlm.nih.gov` answers,
  using the proxy settings from the environment; `genome` also checks the
  hosts in each package's `fetch.txt` before `datasets rehydrate`. On a
  blocked network the stage now exits 3 at once naming the host; before,
  `datasets` made three attempts of about 8.5 minutes each and the stage
  exited 6. Only a failed connection or connect timeout counts as
  unreachable; an HTTP error status, a read timeout or a TLS error does not.
  `REPGENR_SKIP_NET_PROBE=1` skips the check.
- `versions --versions-out` writes a date-like value (the GTDB API query
  date) double-quoted, so a YAML 1.1 loader keeps it as a string.
- `assemble` plans a short-read run that ENA labels PAIRED but lists with one
  FASTQ file as single-end. `--assembler shovill` now excuses such a run as
  `unsupported_layout` before downloading it, rather than after; `auto` still
  assembles it with SKESA. An explicit assembler that takes a run's platform
  but not its layout gives the reason `unsupported_layout` in place of
  `unsupported_platform` (#230).
- When every run is excused and at least one was judged (rejected by the
  assembler, polisher or quality gate), `assemble` (without `--append`) no
  longer leaves a genome set an earlier `assemble` call wrote in place: after a
  warning, `genomes/` is emptied, `selection.tsv` keeps only its header, the
  manifest lists no genomes, and `assembly_stats.tsv` and that call's outgroup
  are removed. The stage still exits 3, and `dereplicate` then exits 3 instead
  of running on the earlier genomes. A set written by another stage, or kept
  because every run failed to download, stays and the error names it; a first
  call writes no empty set. Finished runs under `assemblies/` stay (#230).
- `dereplicate --tool galah` (#233): without manifest quality for every genome,
  galah receives the genomes by descending file size instead of by name. galah
  keeps the first listed genome of a cluster when it has no quality, so a
  fragment whose name sorted first became the representative. With quality for
  every genome the order is unchanged and `--keeper quality` picks the
- `dereplicate --tool galah` (#233): when the manifest (or `selection.tsv`)
  has completeness and contamination for every genome of the run, they are
  passed as `--genome-info` (galah 0.4 and later), so galah ranks genomes by
  quality for both the representative and the greedy membership, also under
  `--keeper tool`. Otherwise galah receives the genomes by descending file
  size instead of by name: galah keeps the first listed genome of a cluster
  when it has no quality, so a fragment whose name sorted first became the
  representative. A finished dereplication is not re-made by this change;
  `--force` applies it.
- `dereplicate --tool drep` (#233): when the manifest (or `selection.tsv`, for
  `dereplicate-chunk` and `dereplicate-merge`) has completeness and
  contamination for every genome of the run, they are passed as
  `--genomeInfo`, so dRep runs without CheckM and scores genomes with the
  values `--keeper quality` uses. `--virus` still passes
  `--ignoreGenomeQuality` instead.
- Dereplicators receive genome quality through a new `DerepParams.quality`
  field (#233). Whether a run passes it is decided once, over every genome of
  the run, so all chunks and the merge pass use the same source: with values
  for only some genomes none are passed (dRep then runs CheckM, galah orders
  by size) and a warning names genomes without values. `--keeper quality`
  still uses the partial values.
- `dereplicate --tool sourmash --process-size N --target-reps M` (#233): the
  merge-level signature collection, whose genome set changes with the
  threshold at each search step, is assembled with `sourmash sig cat` from the
  chunk zips in the shared sketch cache instead of sketching the genomes again.
  Each sketched zip gets an index of its signature names with the size and
  modification time of each genome file; a set the indexed zips do not cover,
  or a failed `sig cat`, is sketched as before. Cache entries (zip digests,
  index entries, and per-genome signatures of the dense back-end) now follow
  the genome file's size and modification time as well as its name, so a
  genome replaced under the same name in a persistent `--tool-arg
  sketch_cache=DIR` is sketched again. The sketch CSV and picklist are
  csv-quoted, so a genome name containing a comma is one field.
- `snptype --tool simple` and `phylo --msa-source snptype --snptyper simple`
  (#232): reference positions that none of a genome's primary or supplementary
  minimap2 alignments covers are N in that genome's consensus, instead of the
  reference base. A core-SNP column is variable only when at least two of A,
  C, G and T occur in it, non-ACGT characters are written as N, and the SNP
  distance matrix counts each pair over the sites where both genomes have a
  base. Assemblies are mapped with the minimap2 preset `asm20` (overridable
  with `--tool-arg preset=<preset>`, or `preset=none` for minimap2's default
  settings; other values exit 2). Positions within a deletion in the genome
  are N as well. A pair of genomes with no shared site has the distance `NA`.
  The log gives each genome's covered fraction of the reference and a summary,
  warns below 50%, and a genome with no base at any core site is refused
  (exit 3) instead of reaching the tree builder as N only. SNP counts,
  `core_snp.fasta`, `full_alignment.fasta` and the distances change: on the 50-genome test set a copy of a genome with 500 kb
  removed differed from it at 20539 sites and now at 0. An existing workdir is
  not retyped by the upgrade alone; `snptype` and `phylo` rerun only with
  `--force` or when their inputs or settings change.
- `--mask gubbins` (#232) passes `--filter-percentage 100` to Gubbins unless
  `gubbins_args` sets it. Gubbins 3.4.3 leaves taxa with more than 25% gaps or
  N out of its analysis by default and still writes them to its outputs, so
  such a taxon was never scanned for recombination. With a user-set value, a
  taxon above it is refused (exit 3) before Gubbins runs.
- Genome record names in the SNP typers' alignments, tree leaves, the phylo
  outgroup leaf and the tree2tax outgroup lookup drop the whole FASTA suffix
  (#232): a genome `x.fasta.gz` is `x`, as in `clusters.tsv`, instead of
  `x.fasta`. Names of uncompressed genomes are unchanged.
- `vmetadata`, `vgenome` (#227): on the NCBI Virus path the species of a
  record is the current species of its taxid in NCBI Taxonomy. `vmetadata`
  looks up every distinct taxid once after the download (`datasets summary
  taxonomy taxon --inputfile`) and stores species, genus and family in
  `virus_records.json`. The species was the organism name, which is often
  strain-level or an earlier name, so one species was split across several
  tokens (Orthohantavirus complete genomes: 107 organism names, 76 species;
  Embecovirus: 56 and 7). The NCBI Virus report lineage nests sibling
  species (Maguari virus under `Orthobunyavirus cacheense`) and can lack the
  current species (Murutucu virus, `Orthobunyavirus maritubaense`), so it is
  only the fallback when the lookup fails or gives no species: the
  shallowest lineage name of the form genus plus one lower-case epithet,
  with the genus above a subgenus and one binomial per taxid, else the
  organism name. Each record's `species_source` names the rule, and the
  `vmetadata` record counts them. A record filed at genus level no longer
  takes the family as its genus. The species sets the filename token,
  `--target-species`, the median-of-medians length window and the outgroup
  candidate groups; the organism name stays in the `description` column of
  `virus_metadata_base.tsv`. A `--target-species` value also selects the
  species of the records whose organism name it is. Viral filenames change.
  A finished workdir keeps its old names while its stages skip; when
  `vgenome` reruns on it, it refuses the old records and asks for one
  `repgenr --force vmetadata` with the same arguments, which renames the
  genomes, after which `vgenome` and the later stages rerun once.
- `vgenome --group-segments` (#227) compares normalised segment labels: the
  text after `;` is dropped, molecule words (`RNA`, `DNA`, `segment`,
  `circular`) are skipped, and small, medium, middle and large become S,
  M, M and L, so `M`, `M; medium` and `middle` are one segment while `RNA 1`
  and `RNA 2` stay distinct; `Unknown` counts as no label. `segments.tsv` gains the columns
  `segment` (normalised) and `segment_label` (as submitted).
- `metadata` (#224): requests through RepGenR's HTTP client (GTDB, also NCBI
  Entrez and ENA) use a 15 s connect timeout and a 120 s read timeout, so a
  blocked network exits 3 after about two minutes instead of eight.
- `metadata` writes `<version>_metadata_r<major>.release` beside a downloaded
  GTDB table, and `--nodownload` refuses (exit 2) a table fetched for another
  minor release; a table without the file is reused with a warning.
- `metadata`: `-r/--release` must be major.minor, `--gtdb-version` must be
  `bac120` or `ar53`, and a `--metadata-path` that does not exist is refused,
  all with exit 2. On `--source api` an unknown taxon or outgroup exits 2, and
  an `--outgroup-accession` inside the selection or the target taxon exits 2
  on both sources.
- `vmetadata` records the source and target of `download.fa` in
  `virus_download_wd/download.source`; the BV-BRC source reuses the group
  FASTA only for the same target, and each source removes the other's tables.
  A BV-BRC workdir written by an earlier version has no such file, so its
  first re-run downloads the group FASTA again.
- `vgenome --group-segments` groups records per (species, isolate) and keeps
  one record per segment; its outgroup search uses the kept records' length
  span widened by 15 percent.
- (#223) On SIGTERM or SIGHUP, repgenr stops the external tools it is running,
  starts no queued tool, removes the partial output and exits with 128 plus
  the signal number (143 for SIGTERM, 129 for SIGHUP); the stage stays
  marked as interrupted. Before, the interpreter ended at once and left the
  tool running (seen with FastTree). A pool of parallel tasks
  (`parallel_map`) now cancels its queued items when one task fails.
- (#231) Every external tool starts in its own session and process group,
  and SIGTERM, SIGHUP, Ctrl-C (SIGINT) and a failure inside repgenr stop the
  whole group: SIGTERM first, SIGKILL after 5 s for what is left. Helpers a
  tool starts are now stopped with it (RAxML or IQ-TREE under
  `run_gubbins.py`, skani under skDER, minimap2 under a typer); before, only
  the tool itself was signalled. Because a tool no longer shares repgenr's
  terminal, repgenr forwards Ctrl-C (raising `KeyboardInterrupt` as before)
  and Ctrl-Z (the tools are suspended and resumed with repgenr). A second
  signal kills the remaining tools at once. A signal inherited as ignored
  stays ignored (`nohup` keeps working), and Ctrl-Z is handled only when
  repgenr is the foreground job of a terminal. Singularity and Apptainer run
  the tool as an ordinary process in the group. With the Docker backend the
  group holds the `docker run` client, which forwards SIGTERM to the
  container; `docker run` now passes `--init`, so the tool is no longer
  process 1 in the container, where it would ignore that signal and keep the
  container running after the client is killed. The changed container
  command line does not make a stage rerun.
- `run` (#223): `--with-snptype --msa-source snptype` runs the `snptype` stage after
  `phylo` instead of before it, so the tables left in `snp/` are the ones the
  `snptype` record describes. As a consequence, a later `phylo` run that
  changes only the tree builder types the genomes again, since `snp/` no
  longer holds phylo's own alignment.
- `phylo` (#223): the MSA stamp version is 3: an alignment cached before the ParSNP and
  cactus record names changed is rebuilt once.
- `phylo` (#223) warns when an alignment-free tree builder (mashtree, sourmash) is
  given `--msa-source snptype` or `--mask`, which it does not use.
- `glance --tool` defaults to `auto`: dRep when it can run (on the `PATH` or
  through the container backend, the test `dereplicate --tool auto` uses),
  sourmash otherwise. The log names the tool picked, and the stage record and
  resume fingerprint hold the concrete tool, so `auto` and naming the same
  tool resume each other. With neither tool available, glance exits 4 and
  names the tools that can compare, after the workdir and genome checks
  (which still exit 3); under an active container backend the message says
  that none of them declares a container image. Before, the default was `drep`.
- `glance --help` states that `--plot-min` and `--plot-max` are Mash ANI
  fractions from 0 to 1, and that `--keep-files` keeps the dRep working
  directory `glance_wd/` (it used the shared "download and scratch
  intermediates" wording). Its description says that glance needs dRep (on
  the PATH, or via the container backend) but no dereplication, and names
  the three output files.
- `derep-unpack` logs once when the file system cannot hard-link and genomes
  are copied (about 4.5 minutes and 2.1 GB for 1000 genomes on exFAT, against
  under a second with hard links on APFS), and its closing line counts the
  genome files and cluster directories written rather than the clusters read.
- `derep-stock --action unpack` hardlinks the representatives from `genomes/`
  as `dereplicate` does, copying only where the file system cannot link. At
  1000 representatives on one volume the unpack took 0.5 s instead of 5 s
  and no extra disk instead of 1.9 GB.
- `cluster_summary.tsv` takes the species from the manifest taxonomy (from
  `selection.tsv` in the Nextflow steps) and falls back to the canonical
  filename. Before, `ingest --selection` with non-canonical filenames gave
  wrong species. Species are told apart by genus and epithet. The `species`
  column lists the keeper's species first, then by number of genomes, and
  stops after five names with `+N more`. `n_species` keeps the full count.
- A failed external tool prints one console line that names the tool and its
  exit status and points to `repgenr.log`. The command line and the output
  tail are written to the run log and shown on the console under `--verbose`;
  a data-channel step, which has no run log, still prints the tail.
- Resume now checks a stage's main outputs before skipping it. A stage whose
  parameters, inputs and environment are unchanged but whose deliverable was
  deleted (for example `derep/clusters.tsv` or `genomes_map.tsv`) logs
  `Stage 'X': deliverable <path> missing; re-running.` and runs again instead
  of being skipped, so `--force` is no longer needed for this case. `doctor`
  reports the same paths as warnings.
- One help string per shared flag: duplicated option texts are now shared
  constants, so `run` and the single commands read the same. `--mask` lists
  the registered maskers on snptype, phylo, run and phylo-build, and the three
  alignment commands say it needs `--msa-source snptype`. The known option-name
  differences (`--tool`, `--platform`, `--metadata-source`, `--viral-source`,
  and `--outgroup-accession` on the step commands) are described in the
  command help and the reference. No option was renamed.
- The sourmash classifier runs its per-genome gathers side by side within the
  thread budget and resolves all of them with one `tax genome` call (a gather
  against a GTDB-sized sketch is single-threaded and takes tens of seconds).
  The `assemble` record now stores the CheckM2 database, GTDB sketch and
  lineages as resolved from the flags or the environment; a re-run of
  `assemble` for quality or classification alone no longer requires the
  assembler binaries; the classifier disagreement warning names the
  compared genera.
- `assemble --jobs` defaults to 2, or to 1 as soon as a long-read run is
  pending, since memory rather than CPU bounds concurrent assemblies. The
  Flye path is verified through its pinned image on simulated Nanopore-like
  reads.
- CI lints and format-checks `scripts/` and `benchmarks/` as well as the
  package and the tests; `scripts/fasta_simulate_sequences.py`, which failed
  both, is rewritten with the same command line.
- The CLI exits with a distinct status per failure class instead of 1 for
  everything: 2 for invalid input, 3 for a bad working directory, 4 for a
  missing or outdated tool, 5 for an unknown adapter, 6 for a failed tool
  (the tool's own status under `REPGENR_PROPAGATE_TOOL_EXIT=1`, as before);
  1 stays for unexpected errors and for `doctor` finding failures. A
  Gubbins failure keeps its tool error and return code, with the divergence
  diagnosis appended, instead of being re-raised as a workdir error.
- Every version-checkable binary now carries a minimum-version floor at
  preflight, and `environment.yml` pins the same floors: dRep 3.0, galah
  0.4, cactus 2.5, parsnp 2.0, harvesttools 1.3, snippy 4.6, Gubbins 3.0,
  FastTree 2.1 and IQ-TREE 2.0 join the floors already declared. dRep's
  version is read from its help banner, as it answers no `--version` flag
  (the preflight reported it as unknown). progressiveMauve, SibeliaZ and
  hal2maf print no parseable version and stay lenient; a test guards the
  set.
- The single-package adapters pin a BioContainers image (galah 0.4.2,
  sourmash 4.9.4, dRep 3.7.1, snippy 4.6.0, ska2 0.5.1, Gubbins 3.4.3,
  IQ-TREE 3.1.3, FastTree 2.2.0, RAxML-NG 2.0.3, mashtree 1.4.6), each
  verified by running its stage through the image, so `--container
  docker|singularity` and the Nextflow container profiles run them in an
  image without Wave. Before, only progressiveMauve and cactus had an image
  and the other adapters fell back to the host with a warning. The
  multi-package `simple` SNP typer and parsnp still need `--wave`, and so
  do skder and SibeliaZ: their BioContainer images are BusyBox-based and
  the GNU-only calls in their shell wrappers (`sort --parallel`, `mktemp
  --suffix`) fail there silently.
- `--wave` now takes precedence over a pinned image when the adapter has a
  conda specification: the flag asks for an image minted from the spec (a
  native-architecture build), and the pin is the default without it.
- The Nextflow phylogeny processes publish the alignment they built
  (`phylo/align/` or `phylo/snp/`, with the reuse stamp) and the tree
  builder's own files under `phylo/tree/`; before, only `tree.nwk` left the
  task directory, so the SNP tables and the tree builder's logs were lost.
- The `cloud` Nextflow profile is removed and the `slurm` profile no longer
  names a queue: an executor's queue, account or region is site-specific and
  belongs in a site config (`-c`). The deprecated `nextflow.enable.dsl = 2`
  line is gone from the pipeline and its test harnesses.
- Documentation repairs from the gap analysis: the README counts five tool
  families and shows the chain `run` executes; `docs/output.md` lists
  `selection.tsv`, `outgroup_accession.txt`, `missing_accessions.txt`, the
  viral tables, the SNP distance matrix, the scratch directory, the aligner
  output under `align/` and the alignment stamp; `docs/verification.md`
  gains a per-adapter status table; `docs/developing.md` describes the layout
  as it is and the adapter checklist in full; the links broken by the
  consolidation are repaired; CITATIONS.md adds CheckM, Prodigal, Mash,
  FastANI and the HAL toolkit; `docs/usage.md` covers the inspection commands
  and names every SNP typer.
- `--mask` is refused before any SNP calling when the typer cannot feed a
  masker (SNP typers declare `produces_full_alignment`; ska2 does not), and
  the masker's binaries are checked before the typer runs. `repgenr run`
  preflights the masker with the other tools, so a missing Gubbins fails in
  the first second instead of after SNP typing.
- Adapters declare the standard parameters they do not pass to their tool
  (`ToolCapabilities.ignored_params`), and every stage warns by name when a
  user sets one of them to a non-default value. Before, `--primary-ani` was
  dropped without a message by skder and sourmash, `--aligned-fraction` by
  dRep and sourmash, and `--bootstrap` by mashtree and the sourmash tree
  builder.
- galah now receives `--primary-ani` as its pre-clustering ANI and
  `--aligned-fraction` as its minimum aligned fraction; both were ignored.
- The thread budget now reaches `sourmash compare` (dereplicator dense path
  and tree builder, `--processes`), FastTreeMP (`OMP_NUM_THREADS`) and
  cactus-pangenome (`--maxCores`); each previously ran on whatever the tool
  chose by itself.
- The documentation is consolidated from fifteen pages to six plus two audit
  records. `docs/usage.md` is the one how-to (command line, Nextflow,
  containers, troubleshooting; it absorbs `containers.md` and the README's
  longer sections), `docs/developing.md` merges `architecture.md` and
  `adding-tools.md`, and the scaling audit and its three SWOT pages are one
  record under `docs/audit/`. `docs/README.md` indexes the pages. The
  executed plans under `docs/superpowers/` and the unreferenced legacy figures
  under `docs/images/` are removed; git history keeps them.
- A stage's first run writes a provisional record to `repgenr.yaml` before its
  body starts, so a stage that fails or is killed shows as `[interrupted]` in
  `status` and as a failure in `doctor`, where it previously looked as if it
  had never started. A successful run replaces the record. A failure in
  parameter validation, a query-only invocation, and a stage that refuses its
  input (exit 2 or 3) without changing its main outputs write nothing.
- `phylo` and `phylo-build` compare the leaves of the built tree with the
  input genomes and exit 3, naming the missing or unexpected leaves and the
  tree builder, when they differ; mashtree can drop a degenerate genome and
  exit 0. The tree is kept for inspection and the stage is not recorded as
  completed.
- `tree2tax` and `tree2tax-relations` exit 3 with a message naming the
  outgroup when it is not a leaf of the tree; they logged a warning, left the
  tree unrooted and exited 0. `tree2tax` after `phylo --no-outgroup` still
  leaves the tree unrooted. `tree2tax-relations --no-outgroup` ignores the
  staged outgroup, and the Nextflow pipeline passes it when `phylo_args`
  contains `--no-outgroup`.
- `ingest --outgroup` naming another genome than the outgroup row of
  `--selection` exits 2 and names both; the selection's outgroup row was
  dropped without a message.

### Fixed
- Chunked dereplication (`dereplicate -s N`) no longer descends level by
  level while each level removes only a few genomes. It chunked the union of
  representatives again whenever it had shrunk at all, in the same fixed
  order, so pairs never placed in one chunk stayed apart: 6888 Bacillaceae
  genomes with `skder --secondary-ani 0.9995 -s 500` ran up to fifty levels
  of about five minutes each (3807, 3627, 3577, 3558, ... representatives).
  A further level now runs only when the previous one removed at least
  max(N, a tenth of its input) genomes; it takes the representatives in turn
  from each previous chunk, so genomes from different chunks meet, and
  compares them at the final thresholds (only the leaf chunks use `--pre-*`).
  Otherwise the union goes to one final pass at the final thresholds,
  logged as "Final pass over M representatives (chunk size N bounds the leaf
  chunks, not the merge)", with a warning above four times N that suggests a
  larger `-s` or a lower `--pre-secondary-ani`. Membership composition also
  keeps QC failures reported by a deeper level, with a warning when a later
  pass rejects a representative that an earlier pass accepted. The Nextflow scatter already
  merged in one call and is unchanged. The reported run removes 3081 genomes
  at level 0, so it now runs one further interleaved level and then the final
  pass, instead of fifty levels. A set where few genomes are removed at level 0
  goes straight to the final pass.
- skDER receives the genome set as a staged directory; the genome count no
  longer limits a single call. The adapter links each genome by its basename
  into a scratch directory and passes that one path to `-g` (it passed one
  path per genome before, and warned above 5000), binds the genome
  directories for a container backend, and refuses two genomes with the same
  basename, which skDER could not tell apart. The ARG_MAX warning for skDER
  is removed.
- Every tool command is measured before it starts: argv (including a
  container backend's own arguments) plus the environment against the
  platform's ARG_MAX less a margin, and on Linux each string against the
  128 KiB single-argument limit. A command over the limit stops with exit 2
  and a message that names `--process-size`, instead of failing in exec with
  "Argument list too long". `warn_argv_bytes` (SibeliaZ, snippy-core)
  warns from three quarters of the measured limit instead of a fixed
  900 kB.
- `sourmash tax genome` (classify, and the reads screen of assemble) reads the
  per-genome gather CSVs from a list file (`--from-file`, written to a
  temporary directory) instead of one command-line argument per genome.
- `tree2tax --collapse-support` (and `tree2tax-relations`) judged the wrong
  nodes when the outgroup was not a child of the root of the tree as written.
  Rooting on the outgroup reverses the branches between the old and the new
  root, and dendropy left each support label on its node, so those labels
  described other branches: on
  `(A:0.1,(B:0.2,(O:0.5,C:0.1)0.60:0.03)0.90:0.05,(D:0.1,E:0.3)0.70:0.04);`
  with outgroup `O` and `--collapse-support 0.8`, the clade of A, B, D and E
  (support 0.60) was kept and the split of A, D and E (0.90) was collapsed.
  Supports now follow their splits through the reroot, in `tree2tax` and in
  `phylo --root`.
- `reads --max-runs` kept an arbitrary subset of runs of equal size, which
  depended on the order of the ENA answer; ties are now broken by run accession.
- `census -tg Bacillus_A --runs` printed a header counting the runs of the
  whole NCBI genus Bacillus as if they belonged to the GTDB genus. The header
  now names the taxon of the runs (`19772 runs of NCBI genus Bacillus`), and
  `census --json` gains the key `runs_taxon` (empty when the runs were
  counted under the census taxon itself).
- ENA run searches returned at most 10000 runs: the portal rejects the
  `offset` parameter the client paged with (HTTP 400), so `reads` failed on
  a taxon with more runs. The search now asks for every run in one request
  (`limit=0`).
- Console output (#253): a shortened tool command line now fits 120 columns
  with its timestamp, level and run-log pointer; it was up to about 225. A
  path directly after an option (`--reference /r/ref.fasta`) is kept as that
  option's value and is no longer counted with the input paths that follow.
  Under `--container` the console shows the tool's own command instead of the
  engine options, which alone filled the line. A long adapter prefix is
  capped and the program and its first option are always shown, and a line
  that is still too long ends at a whole token, not inside a count.
- Manifest upgrade (#252): processes that opened one older manifest at the
  same time could both add a new column, and the second stopped with
  "duplicate column name" (observed on exFAT). The upgrade now takes the
  write lock before it reads the schema version and the columns, and a
  column found already present counts as added. This covers the version 1
  to 2 step as well as the new version 3 step.
- Opening a manifest that is not in WAL mode (#252), for example one copied
  or written by hand, from several processes at once failed with "database
  is locked": the lock wait time was set after the switch to WAL mode, and
  SQLite does not wait for that switch. The wait time is now set first, and
  the switch is retried for up to 30 s.
- Console output (#251): a tool command line is shortened on the console, so
  `dereplicate --tool skder`, which passes every genome path on argv, no
  longer prints a line of several thousand characters (6006 for 50 genomes).
  A run of input paths is shown as its first path and a count, other paths by
  their last component, and the line ends with a pointer to the run log, which
  keeps the full command. This applies to every adapter and to chained
  container steps. skDER's warnings for genomes placed under their closest
  representative, and the `derep-unpack` list of missing cluster members,
  show the first five names and a count on the console beyond five and ten
  genomes; the run log lists them all. No stage reruns for this.
- `dereplicate-merge --keeper tool --reduce` (#254) ranked the
  representatives of a taxon by `selection.tsv` quality; it now uses cluster
  size, as the `dereplicate` stage does under `--keeper tool`. On a set with a
  scored 40 percent fragment the merge step kept the fragment where the stage
  kept the 19-member cluster. No workdir stage reruns for this.
- `dereplicate --keeper quality` (#254) recorded `keeper_effective:
  quality` and skipped its warning when the manifest held quality only for
  genomes outside `genomes/`, such as the outgroup. The quality rule now
  counts only the genomes being dereplicated. Representatives are unchanged.
- Docs and help (#254): the keeper applies in the Nextflow chunk and merge
  steps, not only at the merge; its values come from GTDB, `assemble
  --checkm2-db` or `ingest --selection`, not only GTDB; ties, partial quality
  and contiguity are described in usage.md.
- Resume (#239): input directories are digested from their genome FASTA
  files only (the files `list_fasta` returns), so a leftover `x.fasta.tmp`, a
  `.fai` index or a README in `genomes/`, `derep/representatives/` or the
  `ingest` source no longer re-runs dereplicate, snptype, phylo, glance,
  derep-unpack or ingest. A directory holding only FASTA files keeps the
  digest it had, so a clean workdir does not re-run; one recorded with such a
  file present re-runs that stage once (an `ingest` source with `truth.json`,
  for example). `derep-stock --action unpack` digests the stored run's
  tables and representatives file by file and repeats once after the upgrade.
- `phylo`, `tree2tax` and `doctor` (#239): the outgroup is resolved among the
  genome FASTA files under `outgroup/` only; a partial `GCF_x.fasta.tmp` could
  win the substring match before.
- `ingest --outgroup` and `assemble --outgroup` (#239): a file whose name has
  no FASTA suffix is refused with exit 2 naming the accepted suffixes, before
  anything is staged or assembled; phylo, tree2tax and doctor would ignore it
  under `outgroup/` and the tree would stay unrooted.
- Resume (#239): `genomes/` or `derep/representatives/` holding only a
  non-FASTA file such as `x.fasta.tmp` counts as a missing deliverable, so the
  stage that writes it re-runs instead of being skipped.
- Containers (#240): with `--container docker`, a container started from a
  parallel worker thread is now stopped when repgenr receives a termination
  signal, also when its tool ignores SIGTERM; previously only the container of
  the main thread was stopped. A second signal starts one detached
  `docker stop -t 0` for every running container before repgenr exits. A
  container whose client is killed on a timeout is stopped as well.
  Containers carry the labels `repgenr.pid` and `repgenr.host`, so those of a
  repgenr ended by SIGKILL can be found with
  `docker ps --filter label=repgenr.pid=<pid>`.
- `phylo`, `snptype` (#249): gzipped genomes, accepted by `ingest` since
  #229, failed in the tools that cannot read gzip. progressiveMauve crashed
  with signal 11, ParSNP and snippy exited with a parse error, and SibeliaZ
  exited 0 with an empty MAF, which the adapter reported as spoa running out of
  memory. Such a tool is now given decompressed copies of the gzipped genomes
  in the stage's scratch directory, removed when it has finished. Adapters
  declare `ToolCapabilities.reads_gzip` (True for cactus, `ska2`, `simple`,
  and the `sourmash` and `mashtree` tree builders; False by default). The
  alignment reuse key and the stage fingerprint use the original paths, so no
  workdir reruns because of this. When SibeliaZ is called directly with gzip
  input and writes an empty MAF, the error now names the gzipped input. Two
  input genomes with one record name (`x.fasta` and `x.fasta.gz`, or
  `x.fna`) are refused with exit 2 before an aligner, SNP typer or
  genome-input tree builder runs, since their records and tree leaves could
  not be told apart. `snptype --reference X` now resolves X in the genome
  set being typed before `representatives/`, so `--all-genomes` with a
  reference that is also a representative types it once (it was typed twice,
  as two records named X). A truncated or corrupt gzip genome stops with exit
  2 naming the file, and free disk space is checked before the copies are
  written.
- Version queries (#246): a binary that does not answer its version flag
  is stopped after `VERSION_TIMEOUT` (30 s; 8 s under `list-tools --check`)
  together with the helpers it started, which run in its own process group,
  and a warning names the binary; its version is `unknown` unless its conda
  package record names one. Ctrl-C or a termination signal to repgenr also
  stops a running query and its helpers. Before, the query waited 30 s
  without a message, and a helper of a shell wrapper kept running.
- Version probing (#235): a version query that crashed with a Python
  traceback, or a number inside a longer token (`python3.12`, `GLIBC_2.17`),
  is no longer recorded as the tool's version. A broken cactus environment
  was recorded as cactus 3.12.0 and passed the 2.5 floor.
- Containers (#235): without `--wave`, an adapter that declares only a conda
  spec runs on the host, so auto-selection now counts it as available only
  when its binaries are on `PATH`. Preflight runs `docker info` once and exits
  4 when the daemon cannot be reached. An image Docker cannot start (exit 125
  with Docker's own error) is named instead of the tool. A set
  `CHECKM_DATA_PATH` is passed into the container, with its directory bound.
  Containers are named and stopped with `docker stop` when repgenr is stopped,
  so a tool that ignores SIGTERM does not keep running. An unknown
  `--container` value exits 2, and container options without effect are named
  in a warning. racon's minimap2 image is recorded.
- `versions` (#235): a tool recorded with different versions by two stages
  is listed once per stage as `tool (stage)` instead of keeping only the last
  value, and a stage that did not finish is named on stderr.
- `list-tools --check` (#235): shows the import error of a broken plugin and
  reports an adapter whose preflight raises unexpectedly on its own line; the
  auto-select warning about a broken plugin is printed once per run.
- Gubbins (#235): the thread count is capped at the CPUs where Gubbins runs,
  since its IQ-TREE tree builder refuses more threads than cores.
- Resume (#234): a stage that refuses (exit 2, 3, or 4 for a tool missing at
  preflight) without changing any declared output leaves its record as it
  was: none on a first run, the last finished one on a re-run. It was left
  `[interrupted]`, so status pointed at it and doctor failed although the
  outputs were those of the finished run.
- Resume (#234): each genome `selection.tsv` lists is a deliverable of the
  stage that wrote the genome set (ingest, assemble, vgenome, as genome
  already did), and each representative in `clusters.tsv` is a deliverable of
  dereplicate. A genome deleted from `genomes/` is restored by the next
  `ingest` run, and a deleted representative by the next `dereplicate`; both
  were skipped while doctor asked for a rerun. A deliverable directory holding
  only dotfiles (`.DS_Store`, exFAT `._` files) counts as empty.
- A malformed `repgenr.yaml` (unparsable, a list, a field of the wrong type)
  ends `status`, `versions` and every stage with exit 3 naming the file, and
  is a `doctor` failure; it gave a traceback with exit 1 (#234).
- `doctor` (#234): an interrupted record without parameters
  (`cluster_summary`) is a failure, as `status` showed it; advice names the
  stage that wrote the genome set; the integrity guards' text no longer
  reaches the console; exFAT `._` files are not counted as leftovers.
- `phylo` (#234): the refusal for fewer than three representatives named a
  lower ANI threshold, which leaves fewer; it now names `--secondary-ani`.
- `snptype --tool simple` and `phylo --msa-source snptype --snptyper simple`
  (#232) exited 1 with a UnicodeDecodeError when the reference was a gzipped
  genome (`.fasta.gz`). The reference copy is decompressed; query genomes are
  passed to minimap2, which reads gzip itself.
- `metadata` (#224): under `--limit` the automatic outgroup could be a target
  genome the cap left out; it now lies outside the target taxon. The API path
  keeps GTDB suffixes (`Bacillus_A`) and lowers a capitalised epithet. A
  network or checksum failure of the table download names its cause instead
  of "check release/version".
- `genome` (#224): a download batch made only of accessions NCBI no longer
  serves is recorded in `missing_accessions.txt` after one attempt instead of
  failing with exit 6; each rehydrated genome is checked against the
  package's `md5sum.txt`; an outgroup NCBI does not serve, or a package
  without its FASTA, exits 3 instead of exit 6 or a silent success; a present
  outgroup is not downloaded again; a genome or outgroup deleted by hand is
  fetched again on resume.
- `vmetadata --source bvbrc` (#224) reused a group FASTA fetched for another
  target or by the NCBI Virus source.
- `vgenome` (#224): `--group-segments` concatenated repeated segment records
  of one isolate and isolates of different species sharing a name; a run
  without an outgroup left an earlier outgroup in `outgroup/` and
  `outgroup_accession.txt`; mashtree was recorded as the tool when it did not
  run.
- `tree2tax` and `doctor` (#223) refuse a `tree.nwk` that holds more than one tree;
  before, two concatenated trees passed and tree2tax used the first.
- `snptype`, `phylo` (#223): ParSNP alignment records are named by genome stem (harvesttools wrote
  `x.fasta` and `x.fasta.ref`), and cactus alignment records are renamed
  back from cactus sample names ('.' replaced by '_'). With either tool,
  tree2tax could not find the outgroup leaf and exited 3, IQ-TREE's `-o` missed
  a versioned outgroup (cactus), and the Gubbins outgroup exclusion missed it
  (ParSNP). phylo also writes the input names back into `tree.nwk` when a
  tool renamed leaves that its leaf check accepts, changing only those labels.
  ParSNP query genomes are hardlinked into scratch instead of copied.
- `tree2tax` (#223) warns when the last phylo run did not finish, since `tree.nwk` is
  then the tree of an earlier run. An outgroup accession that no file in the
  outgroup directory matches is reported as such, not as a missing leaf.
- `derep-stock --action pack` stores the completed `dereplicate` record
  (tool, parameters, tool versions, completion time) as
  `derep/stock/<name>/record.json`, and `--action unpack` re-stamps the
  `dereplicate` record from it. Before, unpack took the tool and parameters
  of the record live at unpack time, so a sourmash run restored after a skDER
  run was reported as skDER by `status` and `versions`. The restored record
  carries the time of the unpack as its completion time, not the stored one,
  which stays in `record.json`. Runs stored without `record.json` fall back
  to the previous behaviour, and pack warns when there is no completed
  `dereplicate` record to store.
- `ingest` refuses (exit 2) an `--outgroup` file from outside `--genomes-dir`
  whose name gives the filename or accession of an ingroup genome; the
  outgroup replaced that genome without a message.
- NCBI assembly filenames (`GCF_000008985.1_ASM898v1_genomic.fna`) give the
  leading GCA_/GCF_ accession and no taxonomy. They parsed as family `GCF`
  with accession `genomic`, so an ingested NCBI Datasets download kept one
  manifest row of many.
- `ingest` refuses genomes that share an accession (for example `x.fasta`
  and `x.fna`) and a selection that lists one file twice; the manifest kept
  one row while `genomes/` and `selection.tsv` kept all.
- `ingest` refuses empty and unreadable genome files (dangling links
  included) before staging anything, lists the files it skips for lack of a
  FASTA suffix (`x.fna.gz`, `X.FASTA`) in a warning, and says when an empty
  source holds subdirectories, which it does not search.
- A malformed `ingest --selection` exits 2 naming the file and line (it
  exited 3, or 1 for a non-numeric quality value or a non-UTF-8 file), and
  the `is_outgroup` column accepts `true`/`false` and `yes`/`no`; other
  values are an error rather than read as 0.
- A refused `ingest` leaves the record of an earlier finished ingest
  complete; `status` and `doctor` reported it as interrupted. A refused first
  ingest leaves no record and no empty manifest. The check runs only when the
  stage will run, so an unchanged ingest still skips.
- `ingest` records `--genomes-dir`, `--selection` and an outgroup file as
  absolute paths, so `doctor` run from another directory no longer reports
  the source as changed. This changes the ingest fingerprint of a workdir
  made with a relative `--genomes-dir` or `--selection`, so such a workdir
  reruns ingest once after the upgrade.
- A re-selection that leaves a genome unchanged (`repgenr --force ingest` on
  the same directory) keeps its dereplication status in the manifest; it was
  cleared while `dereplicate` was skipped as up to date.
- `doctor` judges `.fasta.gz` genomes by their decompressed content (it
  reported every one as not FASTA), and names links left dangling by a moved
  source instead of reporting them as not FASTA.
- `dereplicate --tool skder` keeps a partial genome (for example a 40 percent
  fragment) under its cluster's representative. The aligned-fraction cutoff
  now applies to the member's own aligned fraction, as in skDER; before, one
  such genome left a member without a representative and the stage exited 3.
- `derep/representatives/` links each representative from `genomes/`. With
  skDER and galah it linked the tool's own copy, so every representative was
  stored twice; galah no longer copies its representatives into scratch.
- `dereplicate --virus` with a tool that does not read it (all but dRep) warns
  and names the tools that do. It was dropped without a message, and the
  viral example in the documentation used `skder --virus`.
- External tools run with stdin closed. A tool that asks a question (skDER
  below 80 percent ANI) waited on the terminal with its prompt hidden in the
  log; it now reads end-of-file.
- `dereplicate --tool skder` (and `auto` resolving to skDER) refuses
  `--secondary-ani` or `--pre-secondary-ani` below 0.80 with exit 2 before
  the run, since skDER stops at an interactive question there.
- A `dereplicate` rerun refused before it starts (a selected genome missing
  from `genomes/`, an empty `genomes/`, or skDER's ANI floor) leaves the
  record of the last finished run complete. Before, `status` showed it
  interrupted and `doctor` failed although `derep/` was untouched.
- `dereplicate --tool drep` without CheckM exits 6 with the cause in the log.
  dRep exits 0 without results in that case, and the run ended as an
  unexpected error (exit 1).
- `dereplicate --tool drep --virus` uses ANImf as its secondary algorithm.
  The adapter's declared default (fastANI) reached the extras through the
  stage and hid the virus default. The resume fingerprint holds the
  user's parameters, not adapter defaults, so it does not change: a finished
  `--tool drep --virus` run is skipped on a repeat and keeps its fastANI
  result. Rerun it with `repgenr --force dereplicate ...` to get ANImf.
- `dereplicate --tool drep` reports a gzipped genome under its input name
  (`x.fasta.gz`); dRep names its decompressed copy, and the stage exited 3.
- Genomes dRep's filter removes (`--tool-arg length=N`, CheckM thresholds)
  are `fail_qc` in `genome_status.tsv`. They had no status and the stage
  exited 3.
- A repeat `dereplicate` rebuilds a deleted `genome_status.tsv` or
  `cluster_summary.tsv`. Only `clusters.tsv` and `representatives/` were
  deliverables, so `doctor` asked for a rerun and the rerun skipped.
- `assemble` reused a finished run whatever the settings, so a rerun with
  another `--min-contig-length`, `--assembler`, `--polisher`,
  `--polish-rounds` or `--tool-arg` kept the old contigs. The `assembly.ok`
  marker now records the settings; a higher contig floor filters the finished
  contigs again and any other change assembles the run again.
- `assemble` recorded no tool versions when it reused finished runs, and its
  record named the tool `auto`; the marker keeps the assembler and polisher
  versions and the record names the assemblers used.
- An unreadable `assembly.ok` (cut short by a kill) stopped `assemble` with a
  JSON traceback; the marker is written atomically and an unreadable one
  means the run is assembled again. A run with no contig above the floor no
  longer leaves its reads in scratch.
- A run interrupted during assembly fetched its FASTQ files again; files
  that still match their checksum are kept.
- Long-read runs that ENA labels PAIRED were excused as
  `unsupported_platform`; Flye now takes any layout. racon polished a run
  listed as several FASTQ files with the first file only and medaka refused
  it after the assembly; both now join the files.
- A wrong `--checkm2-db`, `--gtdb-sketch` or `--gtdb-lineages` path, or an
  absent checkm2 or sourmash binary, was found only after every assembly; it
  is now refused before any download (exit 2 or 4), also in `genome-qc`,
  which reads both from `REPGENR_GTDB_SKETCH` and `REPGENR_GTDB_LINEAGES`
  as well. A rerun refused this way leaves the finished record as it was.
- `assemble` refused to rerun over finished runs on a nearly full disk
  although it had nothing to download.
- An assembly excused by the CheckM2 gate was not logged; it is now warned
  about, and so is a run whose download failed, with how to retry it.
- `reads` took runs found by accession whatever their library strategy, so
  an RNA-Seq run could be assembled as a genome; they now pass the WGS
  genomic filter of the taxon query, and dropped runs are named.
- `status` lists an interrupted optional stage (for example a `glance` run
  killed mid-way) as `[interrupted]` with the same hint as a stage of the
  chain. It showed `(incomplete)`, a word the documentation does not use.
- A repeat `glance` reruns when `glance_clustering_dendrogram.pdf` was
  deleted, as other stages do for their deliverables, and `doctor` warns
  about the missing file. Before, glance declared no deliverable, so the
  repeat skipped and the dendrogram stayed missing until `--force`.
- `glance` warns when the comparison tool returns no dendrogram, and that
  the next run repeats the comparison. Before, the dendrogram was left out
  without a message.
- `glance` rejects `--plot-min` above `--plot-max`, or either bound outside
  0 to 1, with exit 2. Before, such bounds selected no values, the run
  exited 0 and removed the plots of the previous run.
- `glance` on a workdir with one genome exits 3 with "glance needs at least
  two genomes" before running dRep. Before, `dRep compare` failed inside
  scipy (empty distance matrix) and glance exited 6.
- `glance` plots count each genome pair once. dRep's `Mdb.csv` lists every
  pair in both orders, so the histogram counts and the number in the plot
  titles were twice the number of pairs. The histogram's x axis is now
  labelled "MASH ANI" and its y axis "Genome pairs" (the y axis was labelled
  "MASH ANI"), and the box plot no longer shows a "1" tick.
- `derep-unpack` removed the previous `derep/unpacked/` before building the
  new one, so a failed or interrupted run left a partial tree. The tree is now
  built beside the old one and swapped in when complete.
- `derep-unpack` failed with exit 1 and left a partial `derep/unpacked/` when
  two representatives shared a file stem (`x.fasta` and `x.fna`); such
  clusters are now named by the full file name.
- `derep-stock --name` help states the accepted run names.
- `derep-stock --action pack` under a name already in the store warns that
  it replaces the stored run; the run was replaced without notice.
- `derep-stock --action unpack` marks the `dereplicate` record incomplete
  before it replaces the derep outputs, so an unpack that is killed half-way
  shows as an interrupted dereplication in `status` instead of a finished
  one with `phylo` suggested next.
- A refused `derep-stock` pack or unpack (invalid or unknown `--name`, a
  workdir without a dereplication, an incomplete stored run) no longer marks
  the record of the last finished pack or unpack as interrupted, which made
  `doctor` report a failure. The checks now run before the resume harness
  touches the record, and an unknown run on unpack lists the stored runs.
- `derep-stock` refuses a `--name` longer than 100 characters with exit 2;
  a name past the file-system limit ended in an unexpected error (exit 1).
- `derep-stock --action list` prints the stored run names on stdout, one per
  line, so the list can be read by a script and also appears under
  `--quiet`; it was a log message on stderr. `--action delete` logs the run
  it removed.
- `cluster-summary` warns when `derep/clusters.tsv` lists no clusters,
  instead of logging `Summarised 0 clusters` at info level.
- `derep-stock --action unpack` rebuilds `derep/cluster_summary.tsv` from
  the restored clusters and the current manifest instead of restoring the
  stored copy. The stored copy carried the CheckM quality of pack time, and
  `cluster-summary` skipped afterwards (its inputs had not changed), so the
  live summary could disagree with the manifest.
- `cluster_summary.tsv` no longer counts a genome without a species as a
  species of its own: a genome with neither manifest taxonomy nor a canonical
  filename adds nothing to `n_species` and `species` (a cluster of such
  genomes reports 0 and a blank).
- The `reads` no-match message names every active filter, including
  `--max-bases`, and the command reference shows list defaults (for example
  `reads --drop-selection`) as comma-separated values.
- Entrez taxonomy enrichment stops at the first batch of taxids that fails
  with a connection error (no HTTP status) and exits 3 with a message that
  names the host that gave no response. Before, every batch was
  retried three times, about 16 minutes for 1050 taxids with the network
  down. HTTP errors are still retried per batch.
- `glance` removes the previous dendrogram and plots once the comparison
  succeeds, so a plot with no similarity in the `--plot-min`/`--plot-max`
  range is absent instead of left from an earlier run. A dRep failure prints
  one console line; its command and output tail go to the run log (a test now
  covers this for glance).
- `derep-stock --action delete` of a run that is not stored exits 3 with a
  message naming the run and listing the stored runs. A repeat delete was
  skipped by the resume check and exited 0; delete is no longer recorded or
  skipped.
- `derep-unpack` warns about each cluster member that is missing from
  `genomes/` (one line listing them when there are more than ten) instead of
  leaving it out silently. `derep-stock --action unpack` of a stored run
  without `cluster_summary.tsv` rebuilds the summary from the restored
  clusters, and drops a live `genome_status.tsv` the stored run lacks,
  instead of keeping the files of the replaced dereplication.
- `assemble` and `assemble-run` with `--polisher auto` print one warning per
  platform when an adapter would polish the runs but its tool is not
  installed (for example medaka for ONT), naming the adapters and the remedy.
  The runs are still assembled, unpolished; before, this happened without a
  message.
- `reads --accession-file` treats text from a `#` to the end of the line as a
  comment, also after indentation or after an accession. `reads` (and `run
  --reads`) parse the accession file and check that a selection is given
  before the workdir is created, so a rejected call leaves no directory or
  log; entry stages in general build their parameters before creating the
  workdir.
- The `snippy` SNP typer names the reference record by its genome instead of
  snippy-core's `Reference`, so the tree has a leaf for the reference genome
  that `tree2tax` can map to an accession.
- `genome --accession-list-only` writes `ncbi_acc_download_list.txt` with a
  newline after the last accession, so `wc -l` and `while read` loops see
  every accession.
- `vgenome --outgroup-treebuilder` help names the accepted values, taken from
  the registered tree builders with distance-matrix support.
- `phylo` and `phylo-build` exit 3 with "A tree needs at least 3 genomes" when
  the ingroup has fewer than three genomes, before any aligner, SNP typer or
  tree builder runs, instead of failing inside the tree builder (for example
  quicktree on a single representative).
- A BV-BRC group download is written to a temporary file and renamed after
  the size check, so an interrupted transfer no longer leaves a partial
  `download.fa` that the next run reuses as complete. GTDB table downloads
  and the NCBI Virus path were already protected.
- `status` and `doctor` exit 3 with the shared "Workdir not found" message
  when `-wd` does not exist, instead of exiting 0; an existing directory
  without `repgenr.yaml` still prints the entry-stage hint (status) or a
  warning (doctor).
- Paired runs that ENA lists with a third, orphan FASTQ file (260 of the 424
  paired Illumina Wolbachia runs) no longer fail in skesa: the adapters find
  the `_1`/`_2` pair; skesa also takes the orphan file as unpaired input,
  shovill uses the pair only. The per-sample run choice prefers a long-read
  run only when it carries at least 100 Mb and a tenth of the sample's
  largest short-read run (a 45 kb PacBio run was chosen over 9 Gb of
  Illumina reads).
- Container mounts bind the real directory behind a symlinked path, so a
  containerised tool can open inputs that Nextflow staged through a linked
  directory (seen with sourmash `dereplicate-merge` under the docker
  backend). The stateless reads steps resolve their directories before any
  path reaches a tool.
- Segment-grouped viral isolates (`vgenome --group-segments`) reached
  `genomes_map.tsv` only as their synthetic `iso-` token; the member segment
  accessions were lost from the deliverable. `vgenome` now writes
  `segments.tsv` (token to member accessions) and `tree2tax` lists the
  members under the isolate's leaf.
- The `phylo` record in `repgenr.yaml` omitted the reference, the masker and
  the adapter tuning, so it could not tell a masked tree from an unmasked
  one; the three are recorded now (`reference`, `mask`, `extra`).
- `vgenome` (both back-ends) deleted `genomes/` before writing the new set
  and wrote `selection.tsv` last, so a crash mid-write left a partial set
  with no selection table and the next stage ran on it without a guard. The
  genomes are now built in a staging directory beside `genomes/` and swapped
  in only when every file is written; a failed run leaves the previous set
  and its table untouched.
- `tree2tax --root-name` was ignored unless `--remove-outgroup` was also
  given; the top node was always labelled `root`.
- The RAxML-NG tree builder published `.raxml.bestTree`, which carries no
  support values, so `tree2tax --collapse-support` never collapsed anything
  on its trees. It now publishes the support-annotated tree `--all` writes.
- dRep in virus mode passed `--S_algorithm` twice, so a `S_algorithm=`
  tool-arg was silently replaced by `ANImf`.
- `derep-stock --action unpack` restored the derep tables and representatives
  but left the `dereplicate` record with its previous fingerprint, the
  manifest with the previous per-genome status and `cluster_summary.tsv`
  absent. It now stores and restores the summary, refreshes the manifest and
  re-stamps the record without a fingerprint so the next `dereplicate` runs.
- The viral path (`vgenome`, both back-ends) never wrote the SQLite manifest;
  `dereplicate` then created an empty one, `--reduce species|genus` was a
  silent no-op and `doctor` reported every selected genome as missing.
- The contract TSVs (`clusters.tsv`, `genome_status.tsv`, `selection.tsv`,
  `tree2tax.tsv`, `genomes_map.tsv`) ended their rows with `\r\n`, the
  `csv` module's default, so the last column of every row carried a stray
  `\r`: `awk -F'\t' '$2=="x"'` on `clusters.tsv` never matched. They now
  end with `\n`. Readers were unaffected either way.
- Command audit (2026-10-07): `run --with-snptype --snptyper <unknown>` is
  rejected with exit 2 before any stage starts, and `run --with-snptype
  --mask gubbins` reaches the snptype stage instead of exiting 2.
- `versions` on a workdir without `repgenr.yaml` exits 3 and writes no
  fragment. `status` names `metadata`, `vmetadata`, `ingest` and `reads` as
  starting points, and the `-B/--bootstrap` help says what IQ-TREE needs.
- `list-tools` and the version records no longer store a tool's error line
  when the tool rejects its version flag (sibeliaz); the version is recorded
  as unknown.
- A missing workdir exits 3 on `dereplicate`, `tree2tax`, `glance`,
  `derep-unpack`, `cluster-summary`, `derep-stock`, `genome`, `vgenome` and
  `assemble`, instead of a traceback or a silently created workdir. `glance`
  reports a workdir without genomes before it checks for dRep, and records
  the dRep version.
- `genome-fetch` runs the `datasets` preflight without `--versions-out`
  (exit 4 instead of a traceback). A `selection.tsv` or `reads.tsv` that lacks
  required columns exits 3 and names them, instead of a KeyError.
- `excused_runs.tsv` keeps each excused run on one line; a multi-line failure
  reason no longer spans several physical lines.
- `snptype` drops optional outputs (`full_alignment.fasta`,
  `snp_distance_matrix.tsv`, `variants.vcf`) that the current typer does not write, and
  `phylo` clears the previous tree builder's files from `tree/` before a
  rebuild.
- The SibeliaZ wrapper on macOS skips AppleDouble `._*` files and leaves no
  empty block temp files in `align/`.
- `tree2tax` maps a dereplicated member that is also a leaf once, and does
  not log `--include-dereplicated` as a changed input file. A malformed Newick
  tree (`tree2tax`, `tree2tax-relations`) and a nonexistent `--clusters` path
  exit 3.
- `derep-stock`: `pack` refuses a workdir without dereplication outputs,
  `unpack` checks the stored run before it replaces anything, a repeat
  `unpack` reruns when the live dereplication changed, and every action exits
  3 on a nonexistent workdir.
- `docs/output.md` and `docs/usage.md` agree with the runs: `scratch/` is
  written by five stages, tool intermediates of `dereplicate` are under
  `scratch/`, and `glance` compares all genomes.
- `metadata` with a local GTDB table names `--gtdb-version` when it is missing
  (it named `--version`) and reruns when the `--metadata-path` table changes.
  `vmetadata --source bvbrc` reports an unreachable FTP server as exit 3.
- `snptype --tool simple` with no variable sites names the genome count and
  the reference, and suggests a closer reference, more divergent genomes, or
  an alignment-free tree.
- `tree2tax` and `tree2tax-relations` exit 3 on a `tree.nwk` with text after
  its final `;`, the rule `doctor` uses to flag a truncated tree; dendropy
  read the first tree and ignored the rest.
- The `tree2tax` record in `repgenr.yaml` names `dendropy` as its tool with
  the library version, and `tree2tax-relations --versions-out` writes the same
  version; the record had no tool and no versions.
- The entry records in `repgenr.yaml` name a tool, which `status` shows:
  `metadata` records `gtdb-api` or `gtdb-table`, `genome` and the NCBI Virus
  `vmetadata` record `datasets`, the BV-BRC `vmetadata` records `bvbrc`, and
  `vgenome` records the outgroup tree builder when the outgroup search ran.
  The NCBI Virus `vmetadata` record keeps `released_after` in its parameters.
- The `ingest` record in `repgenr.yaml` includes `drop_foreign` in its
  parameters.
- `metadata --nodownload` reruns when the GTDB table it reuses from the
  workdir is replaced, and logs the changed input; the table was not a
  declared resume input, so the stage was skipped. An existing workdir run
  with `--nodownload` reruns once, since its record holds no digest of the
  table.

### Changed
- A run whose assembler is not installed is excused as
  `assembler_not_installed` with a warning, and the other runs proceed. With
  `--assembler auto`, `assemble` exits 4 only when no run can be assembled.
- The `simple` SNP typer's core-SNP reduction is vectorised. It compared every
  pair of genomes character by character in Python, which cost about 8 ms per
  pair at 200000 sites and made the adapter's advertised limit of 2000 genomes
  unreachable. The output is byte for byte what it was; 40 genomes went from
  3.7 s to 0.13, and 500 genomes now take 35 s. The work is still quadratic in
  genomes, but each pair is now a vector comparison.
- Under a container backend, the `simple` typer runs each genome's chain of
  tools in one container instead of starting one per command. The live
  containerised test went from 21.2 s to 12.6 for 8 genomes, and its run log
  shows 9 engine invocations where it previously needed 58.
- The `simple` SNP typer runs genomes concurrently and passes a thread count to
  minimap2, samtools and bcftools; it previously mapped one genome at a time on
  one core whatever `--threads` said. On 69 Francisella genomes with 8 threads
  the typing stage went from 10.8 minutes to 2.8.
- The same typer writes its per-genome pileup and calls as compressed BCF and
  deletes each genome's intermediates once its consensus has been read. The
  Francisella run left 10 GB of scratch behind, 8.3 GB of it uncompressed
  pileup VCF; the same run now leaves 190 MB, nearly all of it the whole-genome
  alignment the stage publishes.
- `phylo` stamps the alignment it builds and reuses it when a later run changes
  only the tree builder, the bootstrap or the thread count. Trying a second
  tree builder repeated the whole alignment or SNP-calling step before, which
  on that same set was 11 minutes per attempt. The stamp records the genome
  set, the source settings and the alignment's digest, so any change to those,
  or `--force`, rebuilds it.

### Added
- The Nextflow layer can run the alignment and the tree as separate tasks
  (`--phylo_split_msa`, off by default). The phylogeny was one task, so
  changing the tree builder or the bootstrap repeated the alignment or the SNP
  calling, and both halves shared one resource label. The split gives the
  alignment its own cache entry and `process_high`, the tree `process_medium`.
- `phylo-build --msa-only` builds the alignment and stops, writing
  `msa.fasta`; `phylo-build --msa <file>` builds the tree from an alignment an
  earlier call produced. These are the two halves the Nextflow processes run.
- The Gubbins masker reports how much of the alignment is variable, warns above
  10%, and repeats the figure when Gubbins fails, instead of leaving a bare
  non-zero exit. Measured on Francisella subsets, a set of one species (1-8%
  variable) is masked in minutes while a two-species or genus-level set
  (13-39%) crashes Gubbins' recombination scan.
- `--mask gubbins` reads `--tool-arg gubbins_tree_builder=...`,
  `gubbins_first_tree_builder=...` and `gubbins_args="..."` (Gubbins'
  `--tree-builder`, `--first-tree-builder` and any further arguments). The
  masker's keys count as read in the unread-extras warning of the snptype
  and phylo stages.

### Fixed
- IQ-TREE refused a recombination-masked alignment with "Unknown sequence
  type": masking replaces recombinant bases with N, and a fifth of the
  characters was enough to defeat its guess. The adapter now states `-st DNA`,
  which every alignment this pipeline produces is.
- Gubbins exited before its first iteration on hosts whose RAxML package has
  no multi-threaded build (`raxmlHPC-PTHREADS*`; the osx-arm64 conda build
  among them) whenever `--threads` was above one, and the only trace was one
  stderr line. A native run without such a build now uses IQ-TREE as the
  Gubbins tree builder, or one thread when IQ-TREE is missing too, and logs
  the switch.
- Clearing a scratch or output directory on a volume without native extended
  attributes (an external exFAT disk) could fail with "No such file or
  directory": macOS drops a file's AppleDouble twin the moment the data file
  goes, and the removal walk then tried to unlink the vanished name. Every
  such removal now tolerates entries that disappear mid-walk (found on a
  genus-scale run whose workdir sat on an external volume).

### Documentation
- Help and wording from the command audit: `glance --tool` lists only
  dereplicators that support comparison, the `reads` taxon options say only the
  most specific is used, `tree2tax --node-basename` says internal nodes
  otherwise get hash-derived names, and `--metadata-path` says `-r` and
  `--gtdb-version` are still required. usage.md notes that variable-site-only
  alignments give inflated branch lengths and that `phylo-build --msa-only`
  leaves `snp/` in place; output.md notes the `ska2` k-mer files in
  `scratch/snptype/`.
- `docs/install.md` covers installation: the package and its extras, three
  ways to provide the external tools (one conda environment, several
  environments on `PATH`, containers), a recommendation per situation, a
  per-tool table with conda packages, container pins and databases, Cactus,
  Nextflow and how to verify an installation. The container image, storage and
  notes sections moved there from `docs/usage.md`.
- `docs/choosing-tools.md` helps choose tools by dataset type and size, with the
  measured runs behind each recommendation, the dereplicator and tree builder
  trade-offs, assemblers and polishers per platform, `run` versus Nextflow, and
  how `auto` chooses. It states which limits are declared and not benchmarked.
  A new unit test checks its table of declared genome limits against the
  registered adapters.
- The README installation section points to the install guide, its pipeline
  section shows the four stage chains, and its scalability claim states the
  largest recorded runs (1157 bacterial and 1256 viral genomes).
- `docs/output.md` lists `derep/unpacked/`, `derep/stock/<name>/`, the
  `glance` plots and `glance_wd/`, and `ncbi_acc_download_list.txt`. The
  exit-code table in `docs/usage.md` states that `assemble` and `reads-gather`
  exit 3 when every run was excused and nothing was produced.

## [2.1.0] - 2026-09-10

The first public release of the fork. It carries the 2.0.0 rewrite (kept
below as its own section, never published) plus everything since: the
quality keeper, ska2, the nf-core Nextflow layer, and the 2026-09-10
pipeline audit with its live suite.

### Added
- `run` exposes the stage flags it used to leave at their defaults (D-3):
  `--limit`; `--complete-only`, `--host`, `--released-after` for the viral
  chain; `--process-size`, `--num-processes`, `--reduce`, `--target-reps`,
  `--tool-arg`, `--allow-incomplete` for dereplication; `--all-genomes`,
  `--bootstrap`, `--reference`, `--aligner-arg`, `--mask` for the
  phylogeny; `--collapse-support` and `--collapse-length` for tree2tax.
  Each goes through the same builder as the per-stage command, so
  validation and the resume fingerprint are identical.
- `phylo-build --mask` (D-6): the stateless step accepts the same
  recombination masker as `phylo` for `--msa-source snptype`, so the
  Nextflow layer can mask through `phylo_args`.
- `docs/cli-reference.md`: every command with its options, defaults and
  help, generated from the command tree by `scripts/render_cli_matrix.py`
  and kept in sync by a test; a second test checks that every flag is
  mentioned in at least one document besides the audit matrix.
  `docs/verification.md` is now a how-to plus a results table that
  `scripts/live_report.py` renders from the live suite's junit output.
- The CLI matrix carries a `nextflow` column: for every flag, the
  `params.*` key or `ext.args` string that reaches the Nextflow layer, or
  the reason it does not; a test checks each key against
  `nextflow_schema.json`. A stub test asserts the dereplication `ext.args`
  composed from `params.derep_*` and the viral mode reaches the process.
  Live Nextflow tests (`tests/live/test_nextflow.py`): the local
  data-channel harness with two dereplicator and tree-builder pairs,
  `main.nf` in bacterial and viral mode, and the docker profile with
  progressiveMauve.
- Container live tests (`-m "live and container"`,
  `tests/live/test_container_runs.py`): skder, sourmash, dRep (`--virus`),
  the simple SNP typer and sibeliaz in Wave-minted images with
  `--container-cache`, the `REPGENR_CONTAINER*` variables and the rule that
  a native result is never reused by a container run; progressiveMauve and
  Cactus in their pinned amd64 images; `glance --tool drep`; the
  `--container-engine podman` negative path. The module skips itself without
  Docker, amd64 emulation and the Wave CLI.
- `dereplicate-chunk --selection-tsv --keeper` covered live.
- Live tests on the Francisella tularensis species set
  (`tests/live/test_species_set.py`): the simple, ska2 and parsnp SNP typers
  with `--reference`, `--all-genomes`, `--tool-arg` and `--allow-incomplete`;
  `--mask gubbins` changing the core alignment (and refused on ska2);
  IQ-TREE with `-B 1000` support labels and `--no-outgroup`, FastTree and
  RAxML-NG from the SNP source, `--mask` through `phylo`; every workdir
  `tree2tax` flag.
- Network live tests (`-m "live and network"`, `tests/live/test_network.py`):
  GTDB API selections at genus, species (`--limit`, `--outgroup-accession`)
  and family level; the GTDB TSV table with `--nodownload` and
  `--metadata-path`; `genome --accession-list-only` and `--keep-files`; NCBI
  Virus with `--complete-only`, `--released-after`, `--host`; every `vgenome`
  selection flag; BV-BRC `--list`, `--source bvbrc` and `--filter`; `run
  --dry-run` and the bacterial and viral chains end to end. Downloads are
  built once under the live cache directory and reused.
- Live tests for the offline stages (`tests/live/`, 29 tests, about three
  minutes): skder, sourmash and galah recover the synthetic partition; every
  `dereplicate` flag has an observable effect; `dereplicate-chunk` x3 plus
  `dereplicate-merge` by directory and by file list; `phylo-build` and
  `tree2tax-relations` with an outgroup, the collapse thresholds and
  `--versions-out`; the mashtree and sourmash tree builders; `derep-unpack`,
  `derep-stock`, `status`, `versions`, `doctor`; `--force`, `--quiet`,
  `--verbose`, `REPGENR_FORCE` and `REPGENR_LOG_LEVEL`; `ingest --outgroup`,
  `--copy` and `--selection`.
- A CLI matrix, `tests/audit/cli_matrix.yaml`, with one record per command
  and flag (aliases, the parameter it sets, how it is validated, the docs
  that mention it, its live test). `tests/unit/test_cli_matrix.py` checks the
  records against the real command tree on every run: flag and alias sets,
  help text, that each flag reaches its parameter, that a bad value is
  rejected naming the flag, and that live and docs references resolve.
  `scripts/render_cli_matrix.py` renders `docs/audit/cli-matrix.md`.
- `list-tools` now lists the maskers family.
- `repgenr ingest -wd WD --genomes-dir DIR [--selection TSV] [--outgroup NAME|FILE]
  [--copy]`: start a working directory from genomes already on disk. It
  writes the `selection.tsv` and manifest the download stages produce, links
  (or copies) the files into `genomes/`, stages an optional outgroup, records
  itself as a stage (resume, `status` with the local chain, `doctor`), and
  takes taxonomy and CheckM quality from the selection table or the canonical
  filename.
- A live verification suite under `tests/live/` (markers `live`, `network`,
  `container`; deselected by default via `addopts`; `--live-config` maps
  tools to bin directories). The first test runs ingest, sourmash, mashtree
  and tree2tax on a seeded synthetic set and checks the recovered partition
  against the generator's truth.
- `tree2tax --collapse-length L` and `--collapse-support S` (also on the
  `tree2tax-relations` step): internal nodes whose branch is shorter than L,
  or whose support is below the fraction S (percentage trees are normalised),
  merge into their parent before nodes are named, so weak splits do not
  become FlexTaxD nodes. The root, the outgroup/ingroup split and leaves
  never collapse. Provenance records the thresholds and the collapsed count.
- `ska2` SNP typer: reference-free split k-mer calling (`repgenr snptype
  --tool ska2`, `phylo --msa-source snptype --snptyper ska2`). No genome is
  privileged as the reference, so reference-private errors do not bias the
  SNP distances. It emits a variable-site alignment only, so `--mask` is
  refused for it. Tuning via `--tool-arg ksize=` and `min_freq=`.
- `repgenr doctor`: read-only workdir health check that verifies outputs
  against the records in `repgenr.yaml` -- interrupted stages, missing or
  corrupt genomes, manifest drift, representative/cluster mismatches,
  truncated deliverables, unresolvable outgroups, leftover temp files, and
  stages whose inputs changed since completion. Exits 1 on failures.
- `repgenr run` accepts `--msa-source` and `--snptyper`, so the SNP-typing
  phylogeny path (previously manual-only via `repgenr phylo`) is reachable from
  the one-shot orchestrator with identical resume fingerprints.
- **Nextflow nf-core rewrite (Phase 4)**: the pipeline is now a typed
  data-channel workflow with no shared working directory. Parameter schema
  (`nextflow_schema.json`) with nf-schema validation, execution reports and
  nf-core template files; stateless `repgenr dereplicate-chunk` /
  `dereplicate-merge` / `genome-fetch` CLI steps and a portable `selection.tsv`
  hand-off; a scatter-gather dereplication subworkflow; data-channel
  `BACTERIAL_DATAFLOW` (metadata -> genome -> derep -> phylo -> tree2tax) and
  `VIRAL_DATAFLOW` pipelines selected by `--mode`; results published under
  `--outdir`. Stub-based nf-test and a CI job run them on Nextflow 26.04. The
  legacy shared-workdir orchestrator and done-signal modules are removed.
- **Sparse sourmash dereplication back-end**: when the optional
  `sourmash_plugin_branchwater` plugin is installed, the sourmash dereplicator
  uses `manysketch` + `pairwise` to compute only above-threshold edges instead of
  the dense N x N `compare` matrix, keeping memory roughly linear in the number of
  close pairs (relevant at 10k+ genomes). Selected automatically; falls back to
  the dense `compare` path when the plugin is absent. Both paths yield the same
  cluster partition and representative count for a given threshold. Install via
  the `sparse` extra or the `sourmash_plugin_branchwater` conda package.
- **Quality-aware representative selection.** `--keeper quality` (default)
  re-picks each cluster's representative by CheckM completeness minus five
  times contamination after any dereplicator runs, at the in-process stage,
  the Nextflow chunk level and taxonomy reduction; `--keeper tool` keeps the
  adapter's own pick. GTDB CheckM values are carried into the manifest
  (schema v2, two quality columns) and `selection.tsv` (two optional
  columns). Provenance records `keeper`, `keeper_effective` and the number
  of representatives changed.
- `snp/full_alignment.fasta`: SNP typers emit the whole-genome alignment
  (snippy `core.full.aln`, ParSNP `harvesttools -M`, the simple typer's
  consensuses) as a deliverable and as the input for recombination maskers.
- Nextflow publishes `pipeline_info/software_versions.yml`, collected from
  every process.

### Changed
- Recombination masking with an outgroup (D-10): Gubbins now scans the
  ingroup only and its predicted regions are masked in the whole-genome
  alignment for every genome, outgroup included, before the variable
  sites are extracted. A species-level outgroup used to abort Gubbins'
  scan ("gubbins ... exit 1"), which forced `--no-outgroup` on masked
  trees. The masked alignment for a run without an outgroup is Gubbins'
  own filtered file, as before.
- The Nextflow `standard` profile caps process resource requests to the
  host's CPUs and memory (D-11). A laptop run used to fail with "Process
  requirement exceeds available CPUs -- req: 32" unless `-profile test`
  was given; `slurm` and `cloud` stay uncapped.
- `glance`, `derep-unpack` and `derep-stock` (pack, unpack, delete) record
  a stage like every other command (D-7): `status` lists them under the
  optional stages, `doctor` sees them, and an identical repeat skips on
  the resume records (genomes and clusters are digested as inputs).
  `derep-stock --action list` stays a query and records nothing.
- `vmetadata --filter` is a BV-BRC option (D-4): it has no default on the
  command line, BV-BRC applies "complete genome" when it is unset, and
  passing it with the NCBI Virus source is an error that points at
  `--complete-only`. It used to be accepted and silently ignored there.
- Short options mean one thing everywhere (D-2): `-t` is `--threads` on
  every command that has threads, including `run`; `--target` on
  `vmetadata` and `run`, `--filter` and `--list` on `vmetadata`, and
  `--root-name` on `tree2tax` and `tree2tax-relations` have no short form
  any more (`-t`, `-f`, `-l` and `-r` used to collide with `--threads`,
  `--force`, `--level` and `--release`). Nextflow's `vmetadata_args`
  default and the run id derivation use `--target`.
- `tree2tax-relations --include-dereplicated` is on by default, as on
  `tree2tax` and `run`; `--no-include-dereplicated` turns it off. The
  Nextflow `tree2tax_args` default is therefore empty (D-1).
- With `--container` but no `--wave`, an adapter that only declares a conda
  spec runs on the host; the warning now says so and names the remedy
  (`pass --wave`), instead of claiming the tool declares no image.
- `Tree2taxParams.all_genomes`, a field no code read, is gone. It was part
  of the tree2tax resume fingerprint, so an existing workdir re-runs
  `tree2tax` once (seconds). `GlanceParams.threads` defaults to 16 like the
  CLI, and the three aligners take the MSA filename from the shared
  contract constant. Synthetic sets at very small n no longer carry an
  empty cluster. `benchmarks/run_bench.py` reads its storage root from
  `REPGENR_BENCH_STORAGE` (or `--storage`) instead of a fixed volume path.
- Closed-choice options are validated when the command is parsed rather
  than deep in the stage: `metadata -d/-l/--source` (and the same values
  under `run --dataset/--level/--metadata-source/--viral-source`, reported
  under those names), `vgenome --length-method` and `--outgroup-treebuilder`
  (any tree builder with distance-matrix support), `glance --tool`,
  `derep-stock --action`, and `--mask` on `snptype` and `phylo`, which is
  now checked against the masker registry instead of a fixed list. `phylo
  --mask` with `--msa-source aligner` is an error; before, the mask was
  silently dropped. `vgenome_params` has an explicit signature, so a
  misspelt keyword is a `TypeError` instead of a dataclass error.
- Every command-line option has help text (48 were blank: the ANI thresholds
  and `--threads` on the dereplication commands, the GTDB target flags, most
  `vgenome` selection flags, `glance` plot bounds, `derep-unpack
  --no-representant`). Shared texts live in `cli/base.py` so the same flag
  reads the same on every command. `run --snptyper` lists the registered SNP
  typers like `snptype --tool` does.
- `metadata --limit N` no longer keeps the first N genomes in GTDB file (or
  API) order. It round-robins over species, taking the best CheckM-scored
  genome of every species first, then each species' next best, until N;
  within a species unscored genomes rank last, then the GTDB representative
  flag, then accession. The same `--limit` therefore returns a different,
  better set than before. On the API path the per-genome quality cards are
  fetched for every candidate before the cut.
- **Nextflow layer in nf-core shape.** Every channel carries a run-level meta
  map (`id` from the selection target, `mode`), processes exchange
  `tuple val(meta), path(...)`, tool flags reach processes as
  `task.ext.args` mapped from the unchanged parameters in
  `nextflow/conf/modules.config` (where publishing now lives too), and
  resources and the retry window sit in `nextflow/conf/base.config`. The
  Nextflow floor is 26.04 (`!>=26.04.0`) with nf-schema 2.6.1, and the layer
  lints clean under the strict parser. Command lines and parameters are
  unchanged; per-stage `versions.yml` fragments are no longer copied under
  `--outdir` (the collected `pipeline_info/software_versions.yml` remains);
  code that included a RepGenR module or subworkflow directly must
  adapt to the tuple shapes.
- **Crash/restart hardening (stage audit)**: a stage that crashes while
  re-running is no longer silently skipped on restart (the resume record is
  dirtied before the stage body runs; `repgenr status` shows `[interrupted]`);
  deliverable writes are atomic (a failing tree builder can no longer truncate
  a previous `tree/tree.nwk`); downloads are validated before landing in
  `genomes/` and stages refuse incomplete input sets (`--allow-incomplete` to
  override); the manifest is reconciled on re-selection instead of accumulating
  stale rows; stale outgroups are pruned and resolved by exact accession.
- **Pluggability**: tool lists in `--help` are generated from the registries;
  `--tool-arg key=value` provides tool tuning on dereplicate/snptype and the
  data-channel steps; `auto` selection is container-aware and prefers the
  tightest-fitting tool (its picks change); the Nextflow `derep_tool` enum is
  dropped so pip-installed dereplicators work from the pipeline; recombination
  masking is a plugin family (`repgenr.maskers`, `phylo` gains a real
  `--mask`); `repgenr glance` gains `--tool` via a `Dereplicator.compare`
  hook and the viral outgroup builder is selectable via
  `--outgroup-treebuilder` (a `TreeBuilder.distance_matrix` hook); public
  `Registry.register/unregister`; the dead `threads_param` and
  `Aligner.output_kind` capability fields are removed.
- **Breaking: resume fingerprints are input-aware.** A completed stage is now
  skipped only when its parameters, the digests of its inputs (upstream stage
  outputs, digested from file metadata for genome directories and content for
  small contract files), and the container identity (backend/platform/wave) are
  all unchanged. Re-running an upstream stage automatically re-runs downstream
  stages; the previous timestamp-based "may be stale" warning is removed.
  Workdirs created by older versions re-run each stage once.
- **Breaking: `repgenr run` and the manual commands now build identical stage
  parameters** through shared builders, so the two entry points share resume
  fingerprints. `tree2tax --include-dereplicated` now defaults to on for both
  (previously only `repgenr run` enabled it); pass
  `--no-include-dereplicated` for the old manual behavior. The bacterial
  `repgenr run` now requires `-l/--level` instead of silently passing an empty
  level.
- **Breaking for third-party adapters:** `Masker.mask(full_alignment,
  out_dir, params, logger)` receives the whole-genome alignment and a
  `MaskParams`; `SnpParams.mask` is removed; `xmfa_to_fasta` loses its unused
  `flank` argument.
- Stage-level flags (`--virus`, `--mask`) enter a tool's extras only when that
  tool declares them, so resume fingerprints change only when behaviour
  does; every stage warns by name about extras no selected tool reads.
- Repository hygiene: `ruff format` is enforced in CI (the mechanical
  reformat is listed in `.git-blame-ignore-revs`), CI caches pip and the
  micromamba environment and installs nf-test from a pinned release tarball,
  the wiki binaries are no longer tracked (figures live in `docs/images/`),
  and a test pins `pyproject.toml`, the Nextflow manifest and this changelog
  to one version. The Nextflow floor is `>=23.10.0`, the minimum nf-schema
  2.3.0 supports.

### Fixed
- Cactus wrote Toil's `.toil/` state under the container's HOME, which
  the backend points at the working directory; the aligner now runs in its
  alignment directory, so nothing lands in the launch directory.
- Containers started in the temp mount when the adapter gave no working
  directory, so a stateless step run with relative paths (`phylo-build -o .`
  under the docker profile) could not open its outputs
  (`align/xmfa/...`). The container now starts in the host process's working
  directory, which is also bound.
- Nextflow: a run without an outgroup (viral `--no-outgroup`, or a GTDB
  selection with no outgroup candidate) aborted in the data-channel
  subworkflows with "Invalid method invocation `call`": the optional
  outgroup join emitted a bare meta map. The join now carries a placeholder
  tuple; the VACQUIRE stub honours `--no-outgroup` so a stub test covers it.
- mashtree received every genome path on its command line and failed with
  "Argument list too long" at about 9500 genomes (a viral `run` without a
  completeness filter). It now reads the paths from a file-of-files, and
  the genome directories are declared as container mounts.
- Container runs on a workdir populated by `repgenr ingest` (symlinked
  genomes) failed inside the container with dangling links: only the
  `genomes/` directory was bound, not the directories the links point to.
  Every mounted directory's symlink targets are bound as well.
- The simple SNP typer called variants with bcftools' diploid default, so
  heterozygous calls on haploid bacteria became IUPAC codes in the consensus
  and Gubbins refused the whole-genome alignment ("contains disallowed
  characters"). Calling is haploid now (`--ploidy 1`), and the Gubbins
  masker replaces any remaining non-ACGTN symbol with N before running.
- `snptype` without `--reference` picked the alphabetically first genome
  silently on the workdir path; it now logs the same warning the stateless
  path always did, naming the genome it chose.
- `phylo --msa-source snptype` typed only the ingroup, so the outgroup never
  reached the SNP alignment and the tree could not be rooted on it. The
  outgroup is now typed with the ingroup, as on the aligner path.
- `phylo --treebuilder fasttree --bootstrap N` was silently ignored; N is
  now FastTree's `-boot` resample count for its local support values.
- `vgenome --group-segments` concatenated every set of records that shared
  an isolate name, segmented or not; on hepatovirus (one segment) it turned
  367 complete genomes into 149 (216 records share the isolate name "RNA"
  and every one is labelled segment "ANONYMOUS"). An isolate is grouped
  only when its records carry at least two distinct real segment labels.
- `vgenome` on the BV-BRC source wrote the genomes but no `selection.tsv`,
  unlike the NCBI Virus path and the bacterial stages, so the downstream
  contract was incomplete. It now publishes the same table (taxonomy from
  the Entrez names, outgroup row included).
- The genome stage logged "Discarding non-FASTA download for <accession>
  (error page?)" for every accession after a successful batch on macOS
  volumes without native extended attributes: the extraction loop picked up
  the AppleDouble `._*.fna` twins. They are skipped now.
- `vmetadata --source bvbrc --target hepatitis_e_virus` failed with "Could
  not find virus group": the group name was capitalised on its first letter
  only, but BV-BRC names groups like `Hepatitis_E_virus`. The name is now
  resolved against the server listing, ignoring case.
- `--versions-out` on the stateless steps failed with "No such file or
  directory" when the fragment path sat in the not-yet-created output
  directory; the parent is now created first.
- Six small defects noted in the 2026-09-01 audit's self-review: `phylo`
  publishes `tree/tree.nwk` through the atomic copy used by every other
  deliverable; the Nextflow retry window no longer includes exit 130 and 131
  (a cancelled task is not resubmitted); an unknown `--container` value is a
  user-input error with the valid choices in the message; IQ-TREE refuses
  `--bootstrap` values below its floor of 1000 before running; the Wave image
  cache is keyed by platform as well as conda spec; and a failure inside the
  tool-output read loop kills the child process before the error propagates.
- `tree2tax` roots on the outgroup edge instead of at the outgroup's parent
  node. For the unrooted trees that mashtree, fasttree and the sourmash tree
  builder emit, the old rooting left the root with three children (the
  outgroup, its nearest ingroup neighbour and everything else), so the
  taxonomy had no node for the whole ingroup and one taxon sat beside the
  outgroup. The phylo docstring and architecture note now say where rooting
  happens.
- `metadata --source api` reads CheckM quality from each genome's GTDB card
  (`metadata_gene.checkm2_*`, `checkm_*` fallback); the genomes-detail rows
  carry none, so API selections had empty quality columns and the quality
  keeper silently kept the adapter's picks. `dereplicate` now warns when the
  manifest has no quality and records `keeper_effective` in `repgenr.yaml`.
- `repgenr run` preflights every external tool (dereplicator, tree builder,
  and the aligner or SNP typer when the builder needs an MSA) before the
  first stage, instead of discovering a missing tree builder after download
  and dereplication.
- Alignment-free `phylo` runs no longer record an aligner in provenance, so
  `--aligner` cannot invalidate a mashtree or sourmash tree; mashtree builds
  and distance matrices share one argument set and come from one call;
  streamed tool output goes to DEBUG (the file log keeps it, the console
  shows the pipeline's own messages), with carriage-return redraws logged
  once in their final state.
- Dereplication refuses incomplete results (a genome without a status, an
  unmarked representative, a contained genome with no or several
  representatives) instead of silently dropping genomes; genomes that fail a
  tool's QC filter are carried through chunked composition, taxonomy
  reduction and the Nextflow merge.
- `doctor` is strictly read-only: it opens the manifest in read-only mode and
  never creates, migrates or journals it.
- Nextflow help text names NCBI Virus; `tree2tax` publishes only its TSVs to
  the outdir root.

### Notes
- With the 26.04 floor, moving version reporting to `eval()` outputs on the
  `versions` topic costs nothing in compatibility and is the natural next step
  for the Nextflow layer.

## [2.0.0] - 2026-06-18

First stable release of the v2 rewrite: a modular `repgenr` Python package
(Typer CLI, entry-point plugin registries for dereplicators / aligners / SNP
typers / tree builders), a SQLite genome manifest, `repgenr.yaml` provenance, a
no-shell subprocess layer, and a Nextflow orchestration layer.

### Added
- **Container execution backend**: run any external tool in a pinned container
  (`--container docker|singularity`, `--wave`), with BioContainers and Seqera
  Wave image resolution; per-call `HOME` and `extra_mounts`; container-cache
  control. All tool families verified in containers.
- **Dereplication scaling**: `--process-size` two-stage chunking for any tool;
  stage-1 ANI thresholds (`--pre-primary-ani`/`--pre-secondary-ani`);
  `--num-processes` parallel chunk workers.
- **Resume/idempotency**: stages that already completed with the same parameters
  are skipped; `--force` to re-run.
- **Validation & logging**: enum/range validation of CLI options; `--verbose`/
  `--quiet`/`REPGENR_LOG_LEVEL`.
- **Tool version floors**: `min_version` preflight enforcement for the tools with
  reliable version strings; lower-bounded `environment.yml`.
- Nextflow: CPU-matched threads, dynamic resources + retry, first-class chunking
  params, alignment-free default tree builder.

### Changed
- Dereplication post-processing is O(n) (two-stage compose) and the sourmash
  clustering uses a numpy matrix; the manifest uses batched transactions + WAL;
  genome staging uses hardlinks (copy fallback). These remove the per-genome
  Python/I-O bottlenecks for 1000-10000 genome sets.

### Fixed
- `maf_to_fasta` produced an empty MSA on versioned accessions; cactus picked a
  per-chromosome HAL and leaked the `_MINIGRAPH_` pseudo-genome; the manifest
  swallowed write errors and had no schema versioning/concurrency timeout.

### Notes
- Reproducibility: tool versions are recorded in `repgenr.yaml`; generate a
  pinned per-platform conda lock for exact reproducibility (see `environment.yml`).
