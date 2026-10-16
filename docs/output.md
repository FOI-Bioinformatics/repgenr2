# Output files

This page describes every file RepGenR2 writes: the contents of a working
directory, stage by stage in pipeline order, followed by the tables of the
inspection commands and the output directory of the Nextflow pipeline. It is
a reference for users who read the results or build on them in scripts.

The CLI writes all stage outputs into one shared working directory
(`--workdir`, `-wd`), the directory every stage reads from and writes into.
The Nextflow data-channel pipeline has no shared working directory: it passes
results between processes as staged channel files and publishes them under
`--outdir` (see [Nextflow output directory](#nextflow-output-directory)).

Terms used below: a representative is the genome kept for a cluster of
near-identical genomes after dereplication; the outgroup is a genome outside
the studied group, used to root the tree; ANI is average nucleotide identity;
an MSA is a multiple sequence alignment; the core genome is the set of genes
present in all or nearly all genomes of a set.

Contents:

- [Working directory layout](#working-directory-layout)
- [Files shared across stages](#files-shared-across-stages)
- [Entry stages: metadata, vmetadata, reads](#entry-stages-metadata-vmetadata-reads)
- [Genome set: genome, ingest, vgenome](#genome-set-genome-ingest-vgenome)
- [assemble](#assemble)
- [Genome sketches](#genome-sketches)
- [glance](#glance)
- [dereplicate and the dereplication commands](#dereplicate-and-the-dereplication-commands)
- [annotate](#annotate)
- [snptype](#snptype)
- [phylo](#phylo)
- [tree2tax](#tree2tax)
- [Census table](#census-table)
- [Machine-readable status](#machine-readable-status)
- [Nextflow output directory](#nextflow-output-directory)

## Working directory layout

The tree below shows where each stage writes in a working directory. Not
every file is present in every run: each stage writes only its own outputs,
and several files depend on options. The sections that follow describe each
file in pipeline order.

```
WD/
  repgenr.yaml                    provenance record (all stages)
  repgenr.log                     run log (all stages)
  manifest.sqlite                 genome manifest
  scratch/                        working files, one subdirectory per stage
  <version>_metadata_r<major>.*   GTDB table (metadata)
  gtdb_api_genomes.tsv            GTDB API candidates (metadata --source api)
  virus_download_wd/              viral records (vmetadata)
  virus_metadata_base.tsv         viral metadata table (vmetadata)
  reads.tsv                       selected sequencing runs (reads)
  selection.tsv                   the selected genomes (entry stages)
  outgroup_accession.txt          outgroup accession
  ncbi_acc_download_list.txt      accessions to download (genome)
  missing_accessions.txt          accessions not returned (genome)
  segments.tsv                    segment groups (vgenome --group-segments)
  genomes/                        genome FASTAs
  outgroup/                       outgroup genome
  assemblies/<run>/               per-run assembly (assemble)
  assembly_stats.tsv              per-assembly metrics (assemble)
  screen_reads.tsv                reads screen (assemble --screen-reads)
  excused_runs.tsv                runs that gave no genome (assemble)
  sketches/                       sourmash sketches, one per genome
  glance_*.pdf, glance_*.png      overview plots (glance)
  glance_wd/                      glance working files (--keep-files)
  derep/                          dereplication (dereplicate)
    representatives/
    clusters.tsv
    genome_status.tsv
    cluster_summary.tsv
    refined.tsv            (only after --refine or derep-refine)
    unpacked/<representative>/    (derep-unpack)
    stock/<name>/                 (derep-stock)
  annotations/<record>/           per-genome annotation (annotate)
  annotations/annotations.tsv     annotation index (annotate)
  snp/                            core-SNP alignment (snptype)
  align/msa.fasta                 whole-genome alignment (phylo)
  tree/                           tree and builder files (phylo)
    tree.nwk
    tree.unrooted.nwk
    msa/                          alignment cache of phylo
  tree2tax.tsv                    taxonomy for FlexTaxD (tree2tax)
  tree2tax_nodes.tsv              per-node table (tree2tax)
  genomes_map.tsv                 accession to node (tree2tax)
```

## Files shared across stages

| Path | Produced by | Description |
|------|-------------|-------------|
| `manifest.sqlite` | metadata, genome, ingest, vgenome, assemble, sketch | Genome manifest (accessions, taxonomy, CheckM quality when known, dereplication status, the genome's sketch record (see [Genome sketches](#genome-sketches)), and the source of each genome: `gtdb`, `ncbi_virus`, `bvbrc`, `sra` or `local`). `ingest` writes `local`, except that a genome taken with `--from-workdir` keeps the source the manifest of that working directory gives it. |
| `repgenr.yaml` | all stages | Provenance: tool name, parameters, resolved tool versions, completion timestamps, and `bin_dirs`, the directory each tool was given with `--bin-dir` (present only when one was used). A tool that reports no version number itself (SibeliaZ) is recorded with the version of the conda package that installed it, read from `conda-meta` in the binary's environment, or as `unknown` without one. A tool run in a container image is recorded with its image reference, and the engine that ran it with the engine's version (`docker: 29.5.3`, or `singularity`, `podman`, `apptainer`); `repgenr versions` and the Nextflow `versions.yml` list the engine once. The metadata record holds the GTDB `release` for the table path; for `--source api`, whose API reports no release number, `release` is null and `api_query_date` holds the UTC time of the query. The ingest record holds `genomes_dir`, `from_workdirs` (the source working directories of `--from-workdir`, empty without it), `sources`, the number of ingroup genomes per manifest source, and `source_states`, one entry per `--from-workdir` with its `path`, the `stage` that wrote its genome set (null without a record) and the `state` of that record when ingest read it (`done`, `stale`, `interrupted`, `no stage record`, `unreadable record` or `unchecked`). |
| `repgenr.log` | all stages | Run log. |
| `scratch/` | genome, assemble, dereplicate, snptype, phylo | Working files, one subdirectory per stage; the typing pass of `phylo --msa-source snptype` uses `phylo_snptype/`, which holds the same kinds of files as `snptype/` and is cleared when phylo types again. The files under `genome_download/` and `assemble/` are removed unless `--keep-files`. `dereplicate/` holds the dereplicator's working files (sketches, pairwise tables); it is kept after the run and cleared when `dereplicate` runs again. Under `snptype/`, each genome's intermediates are removed once its consensus has been read; a genome whose chain failed keeps its own. The `ska2` typer keeps its split k-mer file (`split_kmers.skf`) and sample list there, and `parsnp` its output directory and `input_genomes/` (links to the genomes where the file system allows, copies otherwise); they are left in place and cleared when `snptype` runs again. For gzipped genomes, `parsnp`'s `input_genomes/` holds decompressed copies. A tool that cannot read gzip (progressiveMauve, SibeliaZ, `parsnp`, `snippy`) is given decompressed copies of gzipped genomes in `inputs/` under `snptype/` or `phylo_snptype/`, or in `phylo_inputs/`; they are removed when the tool has finished. |

## Entry stages: metadata, vmetadata, reads

| Path | Produced by | Description |
|------|-------------|-------------|
| `<version>_metadata_r<major>.tsv.gz` (or `.tar.gz`), `<version>_metadata_r<major>.release` | metadata | The downloaded GTDB table (`--source tsv` without `--metadata-path`) and the exact release it was fetched for, which `--nodownload` checks before reusing the table. |
| `gtdb_api_genomes.tsv` | metadata | With `--source api`: the genomes the GTDB API returned for the target before `--limit` (accession, family, genus and species tokens, `gtdb_representative`, and the family, genus and species as GTDB writes them). `repgenr census -wd` counts them as the candidates of the selection. Removed by a run with `--source tsv`. |
| `virus_download_wd/` | vmetadata | Downloaded viral sequences and the metadata tables `vgenome` selects from, and `download.source` naming the source and target of `download.fa`. On the NCBI Virus path `virus_records.json` holds one record per sequence with its report lineage, its species (the current NCBI Taxonomy species of its taxid, else the ICTV binomial from the lineage, else the organism name), `species_source` naming which of the three applied, `species_name` (the species as NCBI writes it, for display by `census`; empty in records written before it), and the segment label as submitted. `virus_metadata_base.tsv` (and `virus_metadata_ncbi.tsv` on the BV-BRC path) at the root of the working directory are copies of those tables; on the NCBI Virus path its `name` column is the species and `description` the organism name. |
| `reads.tsv` | reads | The selected sequencing runs: run, sample and study accessions, organism and taxid, the resolved family/genus/species tokens, platform, instrument, layout, bases, the FASTQ locations, checksums and sizes ENA reports (empty when ENA holds no FASTQ mirror), and ENA's library selection (RANDOM, MDA, PCR, ...). |

## Genome set: genome, ingest, vgenome

These files form the genome set that every later stage reads. On the reads
path, `assemble` writes the same files (see [assemble](#assemble)).

| Path | Produced by | Description |
|------|-------------|-------------|
| `selection.tsv` | metadata, ingest, vgenome | The selected genomes: accession, taxonomy, outgroup flag, canonical filename, quality when known, and, as the last column, `gtdb_representative` (1 for a genome GTDB marks as its species representative, 0 otherwise and on every non-GTDB entry path; read by `--keeper gtdb`). On the NCBI Virus path the species is the current NCBI Taxonomy species of the record's taxid (see usage.md, Viruses). The hand-off every later stage reads, and a resume input. |
| `outgroup_accession.txt` | metadata, ingest, vgenome | Accession of the outgroup genome, read by phylo and tree2tax for rooting. Absent when there is no outgroup. |
| `genomes/` | genome, ingest, vgenome | Genome FASTAs, one per selected accession. `genome` and `vgenome` write `.fasta`; `ingest` stages the source files under their own names, with the suffixes `.fasta`, `.fa`, `.fna`, `.fas`, `.fasta.gz`, `.fna.gz` or `.fa.gz`. |
| `outgroup/` | genome, ingest, vgenome | Outgroup genome for rooting. |
| `ncbi_acc_download_list.txt` | genome | Accessions the download step will fetch, one per line (those without a FASTA in `genomes/`). With `--accession-list-only` it is the only output. |
| `missing_accessions.txt` | genome | Accessions the download did not return, or whose FASTA failed the package checksum; the completeness guard of later stages reads it. |
| `segments.tsv` | vgenome | With `--group-segments`: each grouped isolate's token and its member segment accessions, one row per member, with the normalised segment (`segment`) and the label as submitted (`segment_label`). Absent otherwise. |

## assemble

| Path | Produced by | Description |
|------|-------------|-------------|
| `assemblies/<run>/` | assemble | Each assembled run's filtered contigs and its `assembly.ok` marker (assembler and polisher, their versions, metrics, and the settings the run was built with); a re-run reuses a run whose settings agree and assembles it again otherwise. With a CheckM2 database, `checkm2.json` stores the run's completeness and contamination with what they were computed from (the SHA-256 of the contigs, the database's resolved path, size and modification time, the CheckM2 version or image); a later call with other `--min-completeness` or `--max-contamination` values applies the stored scores, and CheckM2 runs again only for a run where one of these differs. `--force` reruns the stage but does not bypass the stored scores; to score a run again, delete `assemblies/<run>/checkm2.json`. When the sourmash classifier runs and sketches are not turned off, `contigs.sig.zip` is the run's sourmash sketch with the parameters of `sketches/` (signatures named by run accession) and `contigs.sig.zip.json` records the contigs SHA-256 it was built from; the classifier gathers with it, a later call reuses it while the contigs are unchanged, and the sketch step copies it to `sketches/` under the genome's record name. `reads.sig.zip` is the sourmash sketch of the run's reads (all its FASTQ files, k=21, 31 and 51, scaled=1000, with abundances, named by run accession), made beside the assembler unless `--no-reads-sketch` is given or sourmash is missing. The marker records it as `reads_sketch` (the file name, or null), `reads_sketch_params` and `reads_sketch_version` (sourmash). It is not a genome sketch and is never copied to `sketches/`; an excused run keeps none, except one the reads screen excused. With `--screen-reads`, `screen.json` holds the run's screen decision (`decision` `pass` or `excused`, `reason`, `top_match`, `top_genus`, `fraction`, `duplicate_of`, `ani`, `duplicates_checked`) and the `settings` it was made with (`min_fraction`, `dup_ani`, the resolved `gtdb_sketch` and `gtdb_lineages` paths, `ksize`), and `screen/` the gather table (`gather.csv`) and `tax genome` output (`tax.classifications.csv`) it was made from; the marker of a run that passed holds the same record as `screen`. A later call with the same settings reuses the decision; with other settings the run is screened again from `reads.sig.zip`. |
| `assembly_stats.tsv` | assemble | Per-assembly metrics: assembler, contigs, total length, N50, largest contig, estimated coverage, the NCBI taxonomy used for the name, quality and classification columns once those steps run, and the polisher that corrected a long-read assembly. `label_source` is `classifier` when the GTDB tokens name the genome (the genera agree) and `metadata` otherwise. `taxonomy_flag` is empty when the genera agree, `genus_renamed` when the genera differ but the family and the species epithet agree (an NCBI name that GTDB places in another genus of the same family, such as Mycoplasmopsis arginini, GTDB Metamycoplasma arginini; a GTDB suffix such as `_A` is ignored in the comparison), and `classifier_disagrees` for any other difference, including a shared epithet in another family (Klebsiella pneumoniae against Streptococcus pneumoniae). Both flags are warned about, keep the submitted name, and neither excuses the genome; only `classifier_disagrees` counts towards `n_disagree` in the stage record, which records the renames as `n_genus_renamed`. `reads_sketch` (the last column) is 1 when `assemblies/<run>/reads.sig.zip` holds the run's reads sketch and 0 otherwise; a table written before the column existed reads as 0. `screen` (after it) is `pass` when the reads screen passed the run (`assemble --screen-reads`) and empty when the run was not screened; a table written before the column existed reads as empty. |
| `screen_reads.tsv` | assemble with `--screen-reads`, reads-gather | One row per screened run: `run_accession`, `top_match` (the first gather match, as named in the GTDB sketch), `top_genus` (its GTDB genus), `fraction` (the abundance-weighted fraction of the reads in the top species), `duplicate_of` (the run that contains it, for a duplicate), `decision` (`pass` or `excused`) and `reason` (the excuse, or `genus_renamed` for a pass under the genus tolerance). Rewritten by each call: runs screened or reused in the call get their current row, rows of runs the call did not reach are kept. In the Nextflow reads mode `reads-gather` writes it from the `screen.json` files of the run directories. |
| `excused_runs.tsv` | assemble | Runs that produced no genome, with the step that gave up (`fetch`, `screen`, `assemble`, ...) and the reason; the completeness guard excuses them like `missing_accessions.txt`. Reasons include `no_fastq_mirror`, `download_failed`, `unsupported_platform` (no assembler takes the platform), `unsupported_layout` (the requested assembler takes the platform but not the layout, for example shovill and a run with one FASTQ file), `assembler_not_installed`, `assembly_failed`, `polish_failed`, `qc_failed`, the reads screen's `taxon_mismatch`, `low_match_fraction`, `duplicate_isolate` and `screen_failed` (step `screen`, with `--screen-reads`; see usage.md, Reads screen), and `deferred` (step `assemble`: not attempted because `--max-runs` was reached; a later call with a larger `--max-runs`, or without it, assembles the run and removes the row). When every run is excused, `assemble` exits 3 and, in the cases described below, leaves an empty genome set. |

`assemble` also writes the genome set (`selection.tsv`, `genomes/`,
`outgroup/`, `manifest.sqlite`) on the reads path.

### When every sequencing run is excused

When `assemble` accepts no run, it exits 3, `status` shows the stage as
interrupted, and `excused_runs.tsv` gives the reasons. What happens to the
genome set already in the working directory depends on which stage wrote it
and on whether any run was judged:

- **Written by an earlier `assemble` call, and at least one run judged.** A
  run counts as judged when it was assembled, polished or quality-checked and
  rejected, refused by the assembler, or excused by the reads screen for its
  reads (`assembly_failed`, `polish_failed`, `qc_failed`,
  `unsupported_platform`, `unsupported_layout`, `taxon_mismatch`,
  `low_match_fraction`, `duplicate_isolate`). The set is emptied: a warning
  is logged first, then `genomes/` is emptied, `selection.tsv` keeps only its
  header, the manifest lists no genomes, `assembly_stats.tsv` is removed, and
  the outgroup that call staged is removed with `outgroup_accession.txt`.
  `dereplicate` then refuses with exit 3 (no genome files) rather than run on
  the earlier set. The outputs of an earlier `dereplicate` and `phylo` stay
  until those stages run again.
- **Written by an earlier `assemble` call, but no run judged.** When every
  excuse is `download_failed` (ENA unreachable, for example),
  `no_fastq_mirror`, `assembler_not_installed`, `deferred` or
  `screen_failed`, the set is kept and the error says so. Rerun with
  `--force` once the cause is resolved.
- **Written by another stage** (`genome`, `ingest`, `vgenome`; any manifest
  row whose source is not `sra`) or not recorded in the manifest. The set and
  its outgroup are kept, and the error names their sources.
- **Under `--append`.** Nothing is removed; the existing selection stays as
  it was.
- **No set present** (a first call). Nothing is written: no `genomes/`, no
  `selection.tsv`, no manifest.

The finished runs under `assemblies/` always stay, so a later call with other
settings (a looser quality gate, for example) reuses them.

## Genome sketches

| Path | Produced by | Description |
|------|-------------|-------------|
| `sketches/<name>.sig.zip` | genome, vgenome, ingest, assemble, sketch, and dereplicate, glance and phylo with a sourmash tool | One sourmash zip collection per genome, outgroup included, when sourmash could run (or `--sketch`): three DNA signatures at k=21, k=31 and k=51, scaled=1000, no abundance, each named `<name>`, the genome's file name without its FASTA suffix and `.gz` (its record name in alignments and trees). Written to a hidden `.<name>.partial.sig.zip` and renamed, so an interrupted run leaves no partial sketch. Not a stage deliverable. The manifest records each sketch in three columns of the genome's row (schema version 4): `sketch_file` (the path relative to the working directory), `sketch_params` (`k=21,k=31,k=51,scaled=1000`) and `sketch_digest` (the SHA-256 of the FASTA it was built from); a sketch whose FASTA digest or parameter string differs from its record is stale and written again. `dereplicate --tool sourmash`, `glance --tool sourmash` and `phylo --treebuilder sourmash` read these files and first write the missing and stale ones of the genomes they compare. A manifest of schema version 3 gains the columns, empty, when a stage next opens it; they are not part of the manifest digest of the resume fingerprint. |
| `sketches/.digests.json` | every sketch writer (genome, vgenome, ingest, assemble, sketch) and every consumer (dereplicate, glance and phylo with a sourmash tool) | The SHA-256 of each genome FASTA with the size and modification time it was computed at, so an unchanged genome is not read again to decide whether its sketch is current. Only genomes of the current set are kept. A cache: removed with the last sketch, and rebuilt when absent. |
| `sketches/outgroup.json` | genome, vgenome, ingest, assemble, sketch | The same three fields for an outgroup FASTA under `outgroup/` that has no manifest row, keyed by its file name. Every entry path writes the outgroup to the manifest (the bacterial `genome` stage leaves the filename of the outgroup row empty; that row is matched to the file under `outgroup/` by its accession), so the file is normally absent. |

## glance

| Path | Produced by | Description |
|------|-------------|-------------|
| `glance_clustering_dendrogram.pdf` | glance | Clustering dendrogram over all genomes: dRep's own with `--tool drep`; average-linkage clustering on 1 - ANI with `--tool sourmash`. `--tool auto`, the default, uses dRep when it can run and sourmash otherwise. |
| `glance_MASH_ANI_similarity_boxplot.png`, `glance_MASH_ANI_similarity_histogram.png` | glance | Box plot and histogram of the all-against-all ANI values (Mash ANI from dRep, the sketch-based ANI estimate from sourmash; the axes name which) within `--plot-min`/`--plot-max`, one value per genome pair (self-comparisons left out); the histogram counts genome pairs. A run removes the previous glance plots and dendrogram once the comparison succeeds, so a plot with no values in range is absent rather than stale. |
| `glance_wd/` | glance | Working files of the comparison tool; kept only with `--keep-files`. With `--tool sourmash` it holds `pairwise_ani.csv` (`genome1`, `genome2`, `similarity`, one row per genome pair), `dendrogram_leaves.txt` (the dendrogram's leaf order) and the signatures. |

## dereplicate and the dereplication commands

| Path | Produced by | Description |
|------|-------------|-------------|
| `derep/` | dereplicate | Representative genomes and per-tool intermediates. |
| `derep/representatives/` | dereplicate | The representative genomes, one genome file per cluster; the distinct values of the first column of `clusters.tsv`. |
| `derep/clusters.tsv` | dereplicate | `representative<TAB>member`, one row per genome; a representative also lists itself. |
| `derep/genome_status.tsv` | dereplicate | Per-genome status: `representative`, `contained` or `fail_qc`. |
| `derep/refined.tsv` | dereplicate `--refine`, derep-refine | The genomes of a lineage dereplicated again at a finer ANI: `genome`, `seed` (the seed whose coarse cluster brought it in), `coarse_representative`, `representative` (after the refinement; blank for `fail_qc`), `status` and `ani` (the secondary ANI of that lineage's refinement). Absent when the result holds no refinement. |
| `derep/cluster_summary.tsv` | dereplicate, cluster-summary | One row per representative: member count, species spanned and keeper quality against the members (see [Cluster summary](#cluster-summary)). |
| `derep/unpacked/<representative>/` | derep-unpack | One directory per cluster with its member genomes (and the representative unless `--no-representant`), hard-linked from `genomes/` where the file system allows it and copied otherwise. The directory is named by the representative's file name without its genome extension, or by its full file name when two representatives would otherwise share a directory (`x.fasta` and `x.fna`; names differing only in case). The gzip suffixes (`.fasta.gz`, `.fna.gz`, `.fa.gz`) are removed as a whole. Replaced on each run. |
| `derep/stock/<name>/` | derep-stock | A named, stored dereplication run written by `pack`: `clusters.tsv`, `genome_status.tsv`, `cluster_summary.tsv`, `refined.tsv` when the run holds a refinement, a `representatives/` directory of links to the representative genomes, and `record.json`, the completed `dereplicate` record at pack time (tool, params, tool_versions, completed; absent when there was none). `unpack` restores it, re-stamps the `dereplicate` record from `record.json`, and rebuilds `derep/cluster_summary.tsv` from the restored clusters and the current manifest; the stored summary keeps the quality of pack time. The re-stamped record holds no resume fingerprint, so a later `dereplicate` or `repgenr run` recomputes the restored dereplication; run `phylo` directly to build on it. |

### Cluster summary

`derep/cluster_summary.tsv` condenses `clusters.tsv` into one row per
representative, largest cluster first and then by representative name.
`dereplicate` writes it after every run, and
`repgenr cluster-summary -wd <workdir>` regenerates it for an existing
working directory from `clusters.tsv` and the manifest, without rerunning the
dereplicator. The outgroup is set aside before dereplication and has no row.

| Column | Meaning |
|--------|---------|
| `representative` | Keeper filename, as in `clusters.tsv`. |
| `n_members` | Genomes contained under the keeper; the keeper itself is not counted, so a cluster of 20 genomes reads 19 (`clusters.tsv` lists the keeper as a member of itself). |
| `n_species` | Distinct species across keeper and members. The species comes from the manifest taxonomy (`selection.tsv` in the Nextflow steps) and, for a genome without one, from its canonical filename. A genome with neither adds no species, so a cluster of such genomes reports 0. Species are told apart by genus and epithet. |
| `species` | Those species, comma-separated: the keeper's first, then by number of genomes and name. At most five are listed, followed by `+N more`; `n_species` gives the full count. An epithet shared by two genera is written with its genus, for example `Escherichia coli`. |
| `rep_completeness`, `rep_contamination` | CheckM values of the keeper from the manifest; blank when unknown. |
| `member_max_completeness`, `member_min_contamination` | Best values among the scored members; blank when no member is scored. |
| `best_member` | The genome the quality keeper chooses in the cluster (score `completeness - 5 x contamination + 0.5 x log10(N50)`, ties and the high-quality condition as in usage.md, "Representative selection"). Usually equals `representative` under `--keeper quality`, except where `--reduce` or the Nextflow chunk and merge steps compare only representatives; under `--keeper gtdb` it also differs where a GTDB species representative was kept over a better-scored member. Blank when nothing in the cluster is scored. |
| `n_genomes` | Genomes in the cluster, the keeper included (`n_members` + 1). It follows `best_member`, so the positions of the earlier columns are unchanged. |
| `rep_n50` | N50 of the keeper in bases, read from its FASTA. Filled only when the keeper has completeness and contamination, so an unscored representative is not read; blank otherwise or when the file is not found. |
| `best_score` | Keeper score of `best_member`, four decimals; blank when `best_member` is blank. |
| `rep_is_gtdb_representative` | 1 when the keeper is a GTDB species representative (manifest, or `selection.tsv` in the Nextflow steps), else 0. It is 0 throughout for non-GTDB selections and for a working directory whose selection predates the flag. |
| `refined` | 1 when the cluster comes from a refinement at a finer ANI (`dereplicate --refine` or `derep-refine`, see usage.md, "Finer resolution inside a lineage"), read from `refined.tsv`; else 0. It is the last column. |

The species columns come from the manifest taxonomy (`selection.tsv` in the
Nextflow steps), and from the canonical filename
(`Family_genus_species_ACCESSION.fasta`) for a genome without one. A genome
with neither adds no species. In a cluster of one species `n_species` is 1.
A larger value means the cluster joins several species, or, for names such
as the synthetic benchmark sets, that every genome has its own species
token. With no quality in the manifest (for example after `ingest` without
`--selection` columns) the four quality columns, `best_member` and
`best_score` are blank, and the command logs that it left them blank.

A row whose `best_member` differs from its `representative` marks a cluster
where the quality keeper would choose another genome. This happens under
`--keeper tool`, where the adapter's pick stands, including after
`--reduce`. The quality columns come from the manifest in the working
directory CLI and from `--selection-tsv` in the Nextflow steps; without
either they stay blank, and the size and species columns still apply.

To list the members of a cluster, see "Finding the members of a cluster" in
[usage.md](usage.md#inspecting-a-dereplication).

## annotate

| Path | Produced by | Description |
|------|-------------|-------------|
| `annotations/<record>/<record>.gff3`, `.faa`, `.ffn` | annotate | One directory per annotated record (the representatives, or every genome with `--all-genomes`, and the outgroup): the GFF3 with the sequences appended (`##FASTA`), the proteins and the gene nucleotides, as bakta or prokka wrote them (prokka's `.gff` is copied as `.gff3`). The Nextflow pipeline publishes the same directories under `--outdir`, one per `ANNOTATE` task, without the index. A record is annotated again only when its genome, the annotator, the tool version, the database version or the tuning changed; a record whose genome left the set is kept (it is reused if the genome returns) and one whose genome left the working directory is removed. |
| `annotations/<record>/annotation.json` | annotate, `annotate-build` | The annotator, its version, the database version and the tuning that made the record. `annotate-index` reads it, so gathered records are indexed with the tool that made them; records written by an earlier version lack it and take the `--tool` of `annotate-index`. |
| `annotations/annotations.tsv` | annotate (and, in the Nextflow tasks, `annotate-index`) | The index, one row per record: `record_name`, `set` (genomes or outgroup), `genome_path`, `genome_digest`, `annotator`, `annotator_version`, `db_version` (bakta's database version, empty for prokka), `extra` (the tuning as key=value), `gff_digest`, `completed`. Rewritten after every annotated genome, so an interrupted run keeps what it did. |

## snptype

Only the `snptype` stage writes `snp/`. The typing pass of
`phylo --msa-source snptype` writes the same kinds of files under
`tree/msa/` (see [phylo](#phylo)); the rows below that name both
locations apply to both.

| Path | Produced by | Description |
|------|-------------|-------------|
| `snp/core_snp.fasta` | snptype | Core-SNP (variable-site) alignment; masked in place when `--mask` is set. Only the `snptype` stage writes `snp/`; the typing pass of `phylo --msa-source snptype` writes under `tree/msa/`. |
| `snp/full_alignment.fasta` | snptype | Whole-genome alignment in reference coordinates, when the SNP typer produces one (snippy, parsnp, simple); required input for `--mask`. In the `simple` typer's alignment, reference positions where a genome's alignments place no base (outside them, or in a deletion) are N. |
| `snp/snp_distance_matrix.tsv` | snptype | Pairwise SNP distances between genomes; the `simple` typer writes it, the others do not. Each pair is compared only at sites where both genomes have a base (A, C, G or T); a pair without such a site has `NA`. It is computed before masking, so with `--mask` it counts the recombinant sites that `core_snp.fasta` no longer holds. |
| `snp/variants.vcf` | snptype | The typer's variant calls when it writes any (snippy); under `--snptyper ska2 --tool-arg lo=true` the indels `ska lo` recovered, since its SNPs go into the core alignment. |
| `snp/alignment_stats.json`, `tree/msa/alignment_stats.json` | snptype, phylo | Beside `core_snp.fasta` when the typer wrote a whole-genome alignment: the counts of its constant A, C, G and T columns (`constant_sites`), the number of variable columns, the file they were counted from, whether masking followed and the `--core-fraction` applied (null when none). The tree builder adds the counts back as ascertainment-bias correction (IQ-TREE `-fconst`, RAxML-NG `+ASC_STAM`). Absent for `ska2`. |
| `snp/masker_tree.nwk`, `tree/msa/masker_tree.nwk` | snptype, phylo | With `--mask gubbins`: Gubbins' own final tree (recombination-aware, built by its tree builder over its iterations), kept for comparison with the tree repgenr builds from the filtered sites. It describes the alignment beside it, which under `tree/msa/` may come from an earlier run (the reuse cache). Removed when a later run does not mask. |

## phylo

| Path | Produced by | Description |
|------|-------------|-------------|
| `align/msa.fasta` | phylo | Whole-genome alignment from the aligner (`--msa-source aligner`). Each record is named by the genome file name without its FASTA suffix and `.gz`. |
| `align/msa_source.json`, `tree/msa/msa_source.json` | phylo | Stamp beside the alignment phylo built: the genome set, the source settings and the alignment's digest, so a later `phylo` that changes only the tree builder reuses it. A `snp/msa_source.json` was written by phylo versions that typed into `snp/`; it is no longer read, such a working directory is typed once more into `tree/msa/`, and `doctor` warns that the tables in `snp/` may be phylo's. |
| `tree/` | phylo | Phylogeny (`tree.nwk`) and the tree builder's own files (logs, bootstrap trees). `tree.nwk` is rooted as `phylo --root` asks: on the branch to the outgroup (the default with an outgroup), at the midpoint, or not at all (`none`, the default without one). When phylo rerooted it, `tree.unrooted.nwk` holds the builder's output as written. A rebuild removes the previous builder's files but keeps `msa/`, the reuse cache of the SNP typing pass. |
| `tree/msa/` | phylo | With `--msa-source snptype`: the alignment phylo's own typing pass produced, outgroup included (`core_snp.fasta`, with `full_alignment.fasta`, `snp_distance_matrix.tsv` or `variants.vcf` when the typer writes them, masked in place with `--mask`), and its reuse stamp. It is a reuse cache: after a switch to `--msa-source aligner` or an alignment-free builder it may remain from an earlier snptype-source run, and it then does not describe the current tree; a later return to the same source and settings reuses it. `alignment_stats.json` and `masker_tree.nwk` are described under [snptype](#snptype). |
| `tree/msa/core_gene_alignment.fasta`, `tree/msa/pangenome/`, `tree/msa/core_gene_stats.json` | phylo | With `--msa-source coregene`: the concatenated core-gene alignment the builder wrote (Panaroo's entropy-filtered one unless `filtered=false`), its tables (Panaroo: `gene_presence_absence.csv`, `.Rtab`, `summary_statistics.txt`; PPanGGOLiN: `gene_presence_absence.Rtab`, `matrix.csv`, `genomes_statistics.tsv` and one family list per partition; SCARAP: `genes.tsv`, gene, genome and orthogroup per row; OrthoFinder: `Orthogroups.GeneCount.tsv` and `core_orthogroups.txt`, the orthogroups kept) and a summary (core genes, columns, genomes, core threshold, clustering identity, and for Panaroo the clean mode and filtered; `builder_settings` with PPanGGOLiN's partition, coverage and single-copy choice, SCARAP's method and seeds, or OrthoFinder's search and inflation). A genome left as a row of gaps by the builder is removed from the alignment and reported like one the builder dropped. No `alignment_stats.json`: the alignment keeps its invariant sites. The reuse stamp `msa_source.json` holds the annotation index digest, so a re-annotated genome rebuilds it. |

## tree2tax

| Path | Produced by | Description |
|------|-------------|-------------|
| `tree2tax.tsv` | tree2tax | The taxonomy derived from the tree, for FlexTaxD's `--taxonomy_file`: a header row `child<TAB>parent`, then one edge per row (FlexTaxD accepts either column order and reads the header to tell). Nodes are identified by their label: genome file stems for leaves, a hash of the descendant leaf set (or `--node-basename` plus a number) for internal nodes, and `root` for the top node, which FlexTaxD pre-creates and validates from (`--root-name` changes it, with a warning). Every node has one parent and no field is empty. There is no rank column unless `--ranks` is given: FlexTaxD reads a third column as a rank that must then be present on every row, since an empty one silently inherits the previous row's, so with `--ranks` the header is `child<TAB>parent<TAB>rank` and every row carries `strain` (a genome), `no rank` (an internal node) or a rank `--rank-red` assigned by relative evolutionary divergence. The contract follows the reader on FlexTaxD's GitHub master (`flextaxd/modules/ReadTaxonomy.py`, 2026-10); the PyPI 0.1.5 release reads the genome map with a header line, stores the edges of a child-first file inverted and does not pre-create `root`, so it cannot consume these files. FlexTaxD walks at most 50 levels below `root` (`get_children`, `maxdepth=50`): a taxonomy deeper than that imports and dumps in full but fails `flextaxd --validate`, which an unbalanced tree of more than about 50 genomes can reach; `--collapse-diameter` reduces the depth. |
| `tree2tax_nodes.tsv` | tree2tax | One row per node of the taxonomy, the root included, for reading the taxonomy and choosing collapse thresholds; FlexTaxD does not read it. Columns: `node`, `parent`, `n_leaves`, `depth` (path length from the root), `branch_length`, `support` (the builder's value as a fraction, NA for leaves and trees without supports), `red` (relative evolutionary divergence after Parks et al. 2018 as PhyloRank computes it: 0 at the root, or at the ingroup's common ancestor when the outgroup was removed; 1 at the leaves; with the outgroup kept, the ingroup values depend on where the root sits on the outgroup branch; NA on a tree without branch lengths or without an outgroup or midpoint root, which includes the stateless step without `--keep-root` or an outgroup) and `lineage` (family;genus;species of a leaf's genome from the manifest, empty otherwise or on the stateless step). |
| `genomes_map.tsv` | tree2tax | Accession to leaf, for FlexTaxD's `--genomeid2taxid`: no header, two tab-separated columns, the genome accession and the node it belongs to. One row per representative (its own leaf), per dereplicated member (its representative's leaf, with `--include-dereplicated`), and for a grouped viral isolate per member segment accession. |

## Census table

`repgenr census --tsv PATH` writes one row per genus (family census) or
species (genus census), with a header line and no totals line. The totals are
printed on the console and are part of `--json`. The columns depend on the
mode and appear in this order:

| Column | Mode | Meaning |
|--------|------|---------|
| `name` | all | The genus, or the species, as the source writes them: GTDB `Bacillus_A thuringiensis_S`, NCBI `Francisella sp. LA112445`, NCBI Virus `Hepatovirus ahepa`, BV-BRC `Hepatovirus A`. Filename tokens are shown where only those are known (`selection.tsv` without candidates, a `virus_records.json` from before this field). With `--runs`, a row that holds runs only ends in `(NCBI)` and is listed last. |
| `species` | family census | Distinct species of the genus (GTDB species with `--runs`; selected species in a working directory). |
| `candidates` | selection | Genomes, sequence records or sequencing runs the entry stage chose from. |
| `genomes` | GTDB, selection | GTDB genomes; in a working directory, selected genomes present in `genomes/`. |
| `gtdb`, `sra`, `ncbi_virus`, `bvbrc`, `local` | selection | Selected genomes per manifest source; only the sources present are listed. |
| `representatives` | GTDB, selection | GTDB species representatives. |
| `runs`, `biosamples`, `illumina`, `ont`, `pacbio`, `other` | `--runs` | ENA whole-genome sequencing runs, distinct biosamples, and runs per platform; `other` holds every other platform (BGISEQ, Ion Torrent, ...), so the four platform columns add up to `runs`. |
| `sequences`, `complete` | viral | NCBI Virus or BV-BRC sequence records, and those marked complete (empty for BV-BRC). |
| `isolates`, `segmented` | viral | Genomes after grouping segments by isolate (before vgenome's length filter), and `yes` when the row holds two or more distinct segment labels (empty for BV-BRC). |
| `median_length` | BV-BRC | Median of the per-taxid median sequence lengths. |
| `clusters` | selection | Dereplication clusters whose representative falls in the row (`derep/clusters.tsv`). |

`repgenr census --json` prints one object with these keys:

- `schema` (`repgenr.census/1`) and `repgenr` (version);
- `taxon`, `mode` (`taxon`, `candidates` or `selection`) and `source`
  (`gtdb-api`, `gtdb-table`, `ncbi_virus`, `bvbrc`, `selection` or
  `selection.tsv`);
- `rank` (family or genus) and `row_rank` (genus or species);
- `totals` (genera, species, and the counts of the header line; in the
  selection mode also `sources` and `outgroup`);
- `rows` (objects with the TSV columns, an empty BV-BRC field as `null`);
- `notes`;
- `runs_taxon` (with `--runs`, the taxon the runs were counted under when it
  differs from `taxon`, such as `NCBI genus Bacillus` for the GTDB genus
  `Bacillus_A`; otherwise empty).

## Machine-readable status

`repgenr status --json` and `repgenr doctor --json` print one JSON object on
stdout and nothing else; log messages go to stderr. The `schema` key names the
format and its version, and a change that removes or renames a key raises the
version. The exit codes are the same as without `--json`. On exit 3 (a
missing working directory, or for `status` a malformed `repgenr.yaml`) stdout
is empty and the error is on stderr, so a caller checks the exit code before
parsing. Scripts should key on `level`, `area`, `state`, `pipeline` and
`next`. The `message`, `reason`, `detail` and `notes` texts are meant for
people and may change wording between releases.

`status`:

```json
{
  "schema": "repgenr.status/1",
  "repgenr": "2.1.0",
  "workdir": "/data/wd",
  "pipeline": "local",
  "stages": [
    {"name": "ingest", "in_chain": true, "state": "done", "reason": null,
     "tool": null, "completed": "2026-10-08T07:22:34+00:00",
     "fingerprint": true, "detail": null},
    {"name": "dereplicate", "in_chain": true, "state": "stale",
     "reason": "input changed: genomes", "tool": "sourmash",
     "completed": "2026-10-08T07:30:00+00:00", "fingerprint": true,
     "detail": null},
    {"name": "phylo", "in_chain": true, "state": "pending", "reason": null,
     "tool": null, "completed": null, "fingerprint": false, "detail": null}
  ],
  "next": "dereplicate",
  "notes": [],
  "unchecked": null
}
```

- `pipeline` is `bacterial`, `viral`, `local` or `reads` (the lineage `status`
  follows). It is null when the working directory has no `repgenr.yaml` or
  the record holds no stage (`stages` is then empty and `notes` names the
  entry stages), and when stages are recorded but no entry stage (`stages`
  then follows dereplicate, phylo and tree2tax).
- `stages` lists the stages of that lineage in order (`in_chain` true), then
  any other recorded stage (`in_chain` false), such as `snptype` or `glance`.
- `state` is `done`, `stale` (completed, but an input changed or an output is
  missing; `reason` says which), `interrupted` (started and did not finish) or
  `pending` (not run).
- `fingerprint` is true when the record holds a resume fingerprint. A `done`
  stage without one (written by an older version, or restored by
  `derep-stock --action unpack`) is recomputed by its next invocation.
- `detail` holds the GTDB release or API query date of `metadata`, or the
  source working directories of an `ingest --from-workdir` (followed by the
  sources whose record was not done when they were ingested), else null.
- `next` is the first stage of the lineage that is not done, or null when all
  are. `notes` holds advice about that stage, for example that `phylo` will
  refuse too few representatives.
- `unchecked` is null, or the error that stopped the staleness check; the
  states are then based on the records alone.

`doctor`:

```json
{
  "schema": "repgenr.doctor/1",
  "repgenr": "2.1.0",
  "workdir": "/data/wd",
  "quick": false,
  "findings": [
    {"level": "fail", "area": "genomes",
     "message": "1 selected genome(s) missing from /data/wd/genomes (e.g. GCF_2.1); re-run genome."},
    {"level": "ok", "area": "metadata", "message": "completed 2026-10-08T07:22:34+00:00"}
  ],
  "counts": {"fail": 1, "warn": 0, "ok": 1},
  "exit_code": 7
}
```

- `level` is `fail`, `warn` or `ok`; findings are ordered failures first, then
  warnings, then the rest, each group by `area`, as in the text report.
- `area` is a stage name (`dereplicate`, `phylo`) or a topic such as
  `config`, `genomes`, `manifest`, `outgroup`, `leftovers` or `summary`. A
  check that raised an error is reported as a `fail` in an area named after
  the check (for example `tree` or `manifest_drift`), with a message starting
  "Check could not complete".
- `exit_code` is the status the command exits with: 7 when `counts.fail` is
  above zero, else 0.
- `quick` is true under `--quick`; the genome files were then not read, and a
  non-FASTA file under a FASTA name is not reported.

## Nextflow output directory

The Nextflow pipeline publishes under `--outdir` (default `results/`):

- `metadata/selection.tsv`, or, in reads mode, `reads/reads.tsv`,
  `reads/selection.tsv`, `reads/assembly_stats.tsv`,
  `reads/excused_runs.tsv`, `reads/screen_reads.tsv` (with
  `--assemble_args '--screen-reads'`), and the `quality.tsv` and
  `classification.tsv` tables of the QC step;
- the dereplication contract under `dereplicate/`;
- the phylogeny under `phylo/`: the tree and the tree builder's files under
  `phylo/tree/`, and the alignment it was built from under `phylo/align/` or
  `phylo/tree/msa/`;
- with `--msa-source coregene` in `--phylo_args`, the annotation directories,
  one per `ANNOTATE` task, without the index (see [annotate](#annotate));
- `tree2tax.tsv` and `genomes_map.tsv` at the top level;
- the execution reports under `pipeline_info/`.

Files with the same names as in a working directory have the formats
described in the stage sections above.

### Pipeline information

Under `<outdir>/pipeline_info/`, each run writes timestamped Nextflow execution
reports:

- `execution_report_*.html`: resource usage and per-task summary.
- `execution_timeline_*.html`: task timeline.
- `execution_trace_*.txt`: machine-readable trace of every task.
- `pipeline_dag_*.html`: the workflow DAG.
- `software_versions.yml`: every process's resolved tool versions, collected
  and de-duplicated from each process's `versions.yml` fragment.

These reports help to diagnose resource limits (the retry strategy scales
memory and time per attempt) and record provenance.
