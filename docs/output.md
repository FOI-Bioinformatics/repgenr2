# Output files

This page lists the files RepGenR2 writes and the columns or keys of each.
The CLI writes every stage into one working directory (`--workdir`, `-wd`).
The Nextflow pipeline passes files between processes and publishes them
under `--outdir` (see [Nextflow output directory](#nextflow-output-directory));
a file with the same name has the same format in both.

Terms: a representative (keeper) is the genome kept for a cluster of
near-identical genomes after dereplication; the outgroup is a genome outside
the studied group, used to root the tree; ANI is average nucleotide
identity; an MSA is a multiple sequence alignment; the core genome is the
set of genes present in all or nearly all genomes of a set; a record name is
a genome's file name without its FASTA suffix and `.gz`.

Contents:

- [Working directory layout](#working-directory-layout)
- [Files shared across stages](#files-shared-across-stages)
- [Entry stages: metadata, vmetadata, reads](#entry-stages-metadata-vmetadata-reads)
- [Genome set: genome, ingest, vgenome](#genome-set-genome-ingest-vgenome)
- [Reads path](#reads-path): [reads-sketch](#reads-sketch),
  [reads-screen](#reads-screen), [assemble](#assemble)
- [Genome sketches](#genome-sketches)
- [glance](#glance)
- [dereplicate and the dereplication commands](#dereplicate-and-the-dereplication-commands)
- [annotate](#annotate)
- [snptype](#snptype)
- [phylo](#phylo)
- [tree2tax](#tree2tax)
- [Machine-readable output](#machine-readable-output)
- [Nextflow output directory](#nextflow-output-directory)

## Working directory layout

Each stage writes only its own files, and several files depend on options.

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
  reads_sketches/<variant>/       reads sketches and index (reads-sketch)
  reads_screen/<variant>/         reads screen decisions (reads-screen)
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
    refined.tsv                   (only after --refine or derep-refine)
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

| Path | Produced by | Content |
|------|-------------|---------|
| `manifest.sqlite` | metadata, genome, ingest, vgenome, assemble, sketch | SQLite database with one row per genome in the table `genomes` (columns below). |
| `repgenr.yaml` | all stages | One provenance record per stage (keys below). |
| `repgenr.log` | all stages | Run log. |
| `scratch/` | genome, assemble, dereplicate, snptype, phylo | Working files, one subdirectory per stage (list below). |

`manifest.sqlite`, table `genomes` (schema version 4):

| Column | Meaning |
|--------|---------|
| `accession` | Genome accession, the primary key. |
| `filename` | Genome file name. Empty on the outgroup row of the bacterial `genome` stage, which is matched to `outgroup/` by accession. |
| `source` | `gtdb`, `ncbi_virus`, `bvbrc`, `sra` or `local`. `ingest` writes `local`, except that a genome taken with `--from-workdir` keeps the source its own manifest gives. |
| `family`, `genus`, `species` | Taxonomy tokens. |
| `is_outgroup` | 1 for the outgroup, else 0. |
| `derep_status` | `representative`, `contained`, `fail_qc`, or empty before dereplication. |
| `representative` | The representative of a contained genome. |
| `completeness`, `contamination` | CheckM percentages, when known. |
| `gtdb_representative` | 1 for a GTDB species representative, else 0. |
| `sketch_file`, `sketch_params`, `sketch_digest` | The genome's sketch record (see [Genome sketches](#genome-sketches)). A version 3 manifest gains these columns, empty, when a stage next opens it. |

`repgenr.yaml` holds one record per stage with these keys:

- `tool`, `params`, `tool_versions`, `completed` (UTC timestamp, empty
  while the stage runs or after it failed), `fingerprint` and `inputs` (the
  resume fingerprint and input digests).
- `bin_dirs`: the directory each tool was given with `--bin-dir`; present
  only when one was used.
- `tool_versions`: a tool that reports no version (SibeliaZ) is recorded
  with the version of its conda package from `conda-meta`, else `unknown`.
  A tool run in a container is recorded with its image reference, and the
  engine with its version (`docker: 29.5.3`, `singularity`, `podman`,
  `apptainer`).
- `metadata` `params`: `release`, the GTDB release of the table; with
  `--source api` it is null and `api_query_date` holds the UTC time of the
  query.
- `ingest` `params`:
  - `genomes_dir`;
  - `from_workdirs`, the `--from-workdir` directories, empty without it;
  - `sources`, ingroup genomes per manifest source;
  - `source_states`, one entry per `--from-workdir` with `path`, `stage`
    (the stage that wrote its genome set, null without a record) and
    `state` (`done`, `stale`, `interrupted`, `no stage record`,
    `unreadable record` or `unchecked`).
- `dereplicate` `params`: the request (`secondary_ani`, `keeper`,
  `tool_args` and the like) and `tool_settings`, how the tool computed the
  clustering. For sourmash it holds:
  - `backend`: `sparse` (branchwater plugin), `dense` (`sourmash compare`),
    or `dense,sparse` when the passes of a chunked run differed;
  - `ani_estimate` (`average_containment^(1/k)`), `ksize`, `scaled`;
  - `unknown_pairs` on the dense back-end: pairs with too few shared hashes
    for an ANI; left out for a chunked run.

  `tool_settings` is empty for the other dereplicators. A record restored
  by `derep-stock unpack` keeps the `tool_settings` of the packed record and
  adds `stock: <name>`.

`repgenr versions` and the Nextflow `versions.yml` list each container
engine once.

`scratch/` subdirectories:

- `genome_download/`, `assemble/`: removed unless `--keep-files`.
- `dereplicate/`: the dereplicator's working files (sketches, pairwise
  tables); kept, and cleared when `dereplicate` runs again.
- `snptype/`: each genome's intermediates are removed once its consensus is
  read; a genome whose chain failed keeps its own. `ska2` keeps
  `split_kmers.skf` and its sample list; `parsnp` keeps its output directory
  and `input_genomes/` (links, or copies; decompressed copies of gzipped
  genomes). Cleared when `snptype` runs again.
- `phylo_snptype/`: the same kinds of files for the typing pass of
  `phylo --msa-source snptype`; cleared when phylo types again.
- `inputs/` (under `snptype/` or `phylo_snptype/`) and `phylo_inputs/`:
  decompressed copies of gzipped genomes for tools that cannot read gzip
  (progressiveMauve, SibeliaZ, `parsnp`, `snippy`); removed when the tool
  has finished.

## Entry stages: metadata, vmetadata, reads

| Path | Content |
|------|---------|
| `<version>_metadata_r<major>.tsv.gz` (or `.tar.gz`) | The downloaded GTDB table (`--source tsv` without `--metadata-path`). |
| `<version>_metadata_r<major>.release` | The exact release the table was fetched for; `--nodownload` checks it before reuse. |
| `gtdb_api_genomes.tsv` | With `--source api`: the genomes the GTDB API returned for the target before `--limit`. `census -wd` counts them as candidates. Removed by a run with `--source tsv`. |
| `virus_download_wd/` | Downloaded viral sequences, the tables `vgenome` selects from, and `download.source`, the source and target of `download.fa`. |
| `virus_download_wd/virus_records.json` | NCBI Virus path: one record per sequence (keys below). |
| `virus_metadata_base.tsv`, `virus_metadata_ncbi.tsv` | Copies of the viral metadata tables; the second on the BV-BRC path only. On the NCBI Virus path `name` is the species and `description` the organism name. |
| `reads.tsv` | The selected sequencing runs (columns below). |

`gtdb_api_genomes.tsv`:

| Column | Meaning |
|--------|---------|
| `accession` | GTDB genome accession. |
| `family`, `genus`, `species` | Taxonomy as filename tokens. |
| `gtdb_representative` | 1 for a GTDB species representative, else 0. |
| `family_name`, `genus_name`, `species_name` | Family, genus and species as GTDB writes them, without the rank prefix. |

`virus_records.json`, per record: the report lineage; the species (the
current NCBI Taxonomy species of the taxid, else the ICTV binomial from the
lineage, else the organism name); `species_source`, which of the three
applied; `species_name`, the species as NCBI writes it (empty in records
written before it); and the segment label as submitted. See
[usage.md](usage.md#species-assignment).

`reads.tsv`:

| Column | Meaning |
|--------|---------|
| `run_accession` | ENA or SRA run accession. |
| `biosample`, `bioproject` | Sample and project accessions. |
| `organism`, `taxid` | Submitted organism name and NCBI taxid. |
| `platform`, `instrument_model`, `layout` | Sequencing platform, instrument and library layout. |
| `bases`, `read_count` | Bases and reads ENA reports. |
| `family`, `genus`, `species` | Tokens resolved from the taxid. |
| `fastq_urls`, `fastq_md5`, `fastq_bytes` | FASTQ locations, checksums and sizes, `;`-separated per file; empty when ENA holds no FASTQ mirror. |
| `library_selection` | ENA library selection (`RANDOM`, `MDA`, `PCR`, ...). |

## Genome set: genome, ingest, vgenome

These files form the genome set every later stage reads. On the reads path
`assemble` writes the same files.

| Path | Content |
|------|---------|
| `selection.tsv` | The selected genomes (columns below); the hand-off every later stage reads, and a resume input. Written by metadata, ingest, vgenome and assemble. |
| `outgroup_accession.txt` | Accession of the outgroup, read by phylo and tree2tax. Absent without an outgroup. |
| `genomes/` | One FASTA per selected accession. `genome` and `vgenome` write `.fasta`; `ingest` keeps the source names (`.fasta`, `.fa`, `.fna`, `.fas`, `.fasta.gz`, `.fna.gz`, `.fa.gz`). |
| `outgroup/` | The outgroup genome. |
| `ncbi_acc_download_list.txt` | Accessions `genome` will fetch (those without a FASTA in `genomes/`), one per line; the only output with `--accession-list-only`. |
| `missing_accessions.txt` | Accessions the download did not return or whose FASTA failed the package checksum; read by the completeness check of later stages. |
| `segments.tsv` | With `vgenome --group-segments`: columns `isolate` (the grouped isolate's token), `accession`, `segment` (normalised) and `segment_label` (as submitted), one row per member. |

`selection.tsv`:

| Column | Meaning |
|--------|---------|
| `accession` | Genome accession. |
| `family`, `genus`, `species` | Taxonomy tokens. On the NCBI Virus path the species is the current NCBI Taxonomy species of the taxid. |
| `is_outgroup` | 1 for the outgroup, else 0. |
| `filename` | Canonical file name, `Family_genus_species_ACCESSION.fasta`. |
| `completeness`, `contamination` | CheckM percentages, two decimals; empty when unknown. |
| `gtdb_representative` | 1 for a GTDB species representative; 0 otherwise and on every non-GTDB path. Read by `--keeper gtdb`. |

## Reads path

The reads path writes, in order, reads sketches, screen decisions and
assemblies; `assemble` then writes the genome set above. The methods are in
[reads.md](reads.md).

### reads-sketch

| Path | Content |
|------|---------|
| `reads_sketches/<variant>/` | One directory per depth cap: `all/` for the whole run, `spots_<N>/` for `--max-spots N`. With `-o DIR` instead of `-wd`, the same files under `DIR/<variant>/`. |
| `reads_sketches/<variant>/<run>.sig.zip` | Reads sketch of one run, streamed from NCBI with `fastq-dump` into sourmash: k=21, 31 and 51, scaled=1000, with abundances, named by run accession. |
| `reads_sketches/<variant>/reads_sketches.tsv` | Index, one row per run of `reads.tsv` in its order (columns below); rewritten after every finished run. |
| `reads_sketches/<variant>/.lock` | Held by the running call; a second call on the variant is refused while it is held. |

`reads_sketches.tsv`:

| Column | Meaning |
|--------|---------|
| `run_accession` | Run accession. |
| `status` | `ok`; `failed`; `empty` (no biological read at this cap, no sketch, not streamed again); or `pending` (never selected). |
| `max_spots` | Cap the sketch was made with; empty for the whole run. |
| `spots_read`, `sequences` | Spots `fastq-dump` reported and sequences sourmash reported. |
| `sketch_file` | Sketch file name beside the index. |
| `sketch_params` | sourmash parameter string. |
| `seconds` | Wall time of the pipeline, one decimal; for a retried run, the last attempt's. |
| `reason` | Last error line of a failed run, after up to three attempts for a transient archive failure. |

### reads-screen

| Path | Content |
|------|---------|
| `reads_screen/<variant>/reads_screen.tsv` | Decisions on the sketches of one variant, one row per run of `reads.tsv`, with the columns of `screen_reads.tsv`. Rows of runs a call did not select are kept. |
| `reads_screen/<variant>/passed.txt` | The runs that passed, one accession per line in `reads.tsv` order; the input of `assemble --accession-file`. |
| `reads_screen/<variant>/<run>/screen.json` | One run's decision and settings (keys under [assemble](#assemble)). Absent for a run the tools failed on. |
| `reads_screen/<variant>/<run>/screen/` | The `sourmash gather` CSV and `tax genome` output the decision was made from. |

### assemble

| Path | Content |
|------|---------|
| `assemblies/<run>/contigs.fasta` | Filtered contigs of one run. |
| `assemblies/<run>/assembly.ok` | Completion marker (keys below). |
| `assemblies/<run>/checkm2.json` | With a CheckM2 database: completeness, contamination, the contigs SHA-256, the database's resolved path, size and modification time, and the CheckM2 version or image. |
| `assemblies/<run>/contigs.sig.zip`, `contigs.sig.zip.json` | When the sourmash classifier runs: the contig sketch (parameters of `sketches/`, named by run accession) and the contigs SHA-256 it was built from. Copied to `sketches/` under the record name. |
| `assemblies/<run>/reads.sig.zip` | Reads sketch of all FASTQ files of the run: k=21, 31 and 51, scaled=1000, with abundances. Absent with `--no-reads-sketch`, without sourmash, and for an excused run unless the reads screen excused it. Never copied to `sketches/`. |
| `assemblies/<run>/screen.json`, `screen/` | With `--screen-reads`: the run's decision (keys below), and `gather.csv` and `tax.classifications.csv` it was made from. |
| `assembly_stats.tsv` | Per-assembly metrics (columns below). |
| `screen_reads.tsv` | With `--screen-reads`, or from `reads-gather` in the Nextflow reads mode: one row per screened run (columns below). Rows of runs a call did not reach are kept. |
| `excused_runs.tsv` | Runs that produced no genome (columns below); the completeness check excuses them like `missing_accessions.txt`. |

`assembly.ok` records the assembler and polisher, their versions, the
metrics and the settings of the run; `reads_sketch` (file name or null),
`reads_sketch_params` and `reads_sketch_version` (sourmash); and, for a run
that passed the screen, the `screen` record.

`screen.json`:

- `run_accession`, `decision` (`pass` or `excused`) and `reason` (the
  excuse, or `genus_renamed` for a pass under the genus tolerance);
- `top_match`, `top_genus`, `fraction` (see `screen_reads.tsv`);
- `duplicate_of` and `ani`, the containment ANI to that run;
- `duplicates_checked`, whether the duplicate check was part of the
  decision;
- `settings`: `min_fraction`, `dup_ani`, the resolved `gtdb_sketch` and
  `gtdb_lineages` paths, and `ksize` (31).

`assembly_stats.tsv`:

| Column | Meaning |
|--------|---------|
| `run_accession`, `filename` | Run and genome file name. |
| `assembler` | Assembler that built the contigs. |
| `n_contigs`, `total_length`, `n50`, `largest_contig` | Contig count, total length, N50 and largest contig in bases. |
| `est_coverage` | Estimated coverage. |
| `completeness`, `contamination` | CheckM2 percentages; empty when not scored. |
| `ncbi_taxonomy`, `gtdb_taxonomy` | Taxonomy used for the name, and the classifier's GTDB taxonomy. |
| `label_source` | `classifier` when the GTDB tokens name the genome (the genera agree), else `metadata`. |
| `taxonomy_flag` | Empty when the genera agree; `genus_renamed` when genus differs but family and epithet agree; `classifier_disagrees` otherwise. GTDB suffixes such as `_A` are ignored. |
| `polisher` | Polisher that corrected a long-read assembly; empty when none ran. |
| `reads_sketch` | 1 when `assemblies/<run>/reads.sig.zip` holds the reads sketch, else 0. |
| `screen` | `pass` when the reads screen passed the run; empty when it was not screened. The last column. |

Neither taxonomy flag excuses a genome. The stage record counts
`classifier_disagrees` as `n_disagree` and `genus_renamed` as
`n_genus_renamed`. In tables written before `reads_sketch` or `screen`
existed, the column reads as 0 or empty.

`screen_reads.tsv` (and `reads_screen.tsv`):

| Column | Meaning |
|--------|---------|
| `run_accession` | Run accession. |
| `top_match` | First gather match, as named in the GTDB sketch. |
| `top_genus` | Its GTDB genus. |
| `fraction` | Abundance-weighted fraction of the reads in the top species. |
| `duplicate_of` | Earlier run that contains this one, for a duplicate. |
| `decision` | `pass` or `excused`. |
| `reason` | `taxon_mismatch`, `low_match_fraction`, `duplicate_isolate` or `screen_failed`; `genus_renamed` for a pass under the genus tolerance. |

`excused_runs.tsv`:

| Column | Meaning |
|--------|---------|
| `run_accession` | Run accession. |
| `step` | Step that gave up: `fetch`, `screen`, `assemble`, `qc` or `classify`. |
| `reason` | The reason code, with any tool output folded onto one line (codes below). |

Reason codes: `no_fastq_mirror`, `download_failed`, `unsupported_platform`
(no assembler takes the platform), `unsupported_layout` (the assembler takes
the platform but not the layout), `assembler_not_installed`,
`assembly_failed`, `polish_failed`, `qc_failed`; the screen codes
`taxon_mismatch`, `low_match_fraction`, `duplicate_isolate` and
`screen_failed` (step `screen`, see
[reads.md](reads.md#the-three-checks)); and `deferred` (step `assemble`,
not attempted because `--max-runs` was reached).

When every run is excused, `assemble` exits 3 and may empty the genome set
of an earlier `assemble` call. The cases are described in
[reads.md](reads.md#when-every-sequencing-run-is-excused).

## Genome sketches

| Path | Content |
|------|---------|
| `sketches/<name>.sig.zip` | One sourmash zip collection per genome, outgroup included: k=21, 31 and 51 DNA signatures, scaled=1000, no abundance, each named by the record name. |
| `sketches/.digests.json` | SHA-256 of each genome FASTA with the size and modification time it was computed at; genomes of the current set only. A cache, rebuilt when absent. |
| `sketches/outgroup.json` | The three sketch fields for an outgroup FASTA under `outgroup/` without a manifest row, keyed by file name; normally absent. |

Sketches are written by genome, vgenome, ingest, assemble and sketch, and by
`dereplicate --tool sourmash`, `glance --tool sourmash` and
`phylo --treebuilder sourmash`, which first write the missing and stale
sketches of the genomes they compare. A sketch is written to a hidden
`.<name>.partial.sig.zip` and renamed. The manifest records each sketch in
`sketch_file` (path relative to the working directory), `sketch_params`
(`k=21,k=31,k=51,scaled=1000`) and `sketch_digest` (SHA-256 of the FASTA). A
sketch whose digest or parameter string differs from its record is stale.
The sketch columns are not part of the manifest digest of the resume
fingerprint.

## glance

| Path | Content |
|------|---------|
| `glance_clustering_dendrogram.pdf` | Clustering dendrogram over all genomes: dRep's own with `--tool drep`, average linkage on 1 - ANI with `--tool sourmash`. `--tool auto` uses dRep when it can run. |
| `glance_MASH_ANI_similarity_boxplot.png`, `glance_MASH_ANI_similarity_histogram.png` | Box plot and histogram of the pairwise ANI values within `--plot-min`/`--plot-max`, one value per genome pair. Mash ANI from dRep, sketch ANI from sourmash; the axes name which. |
| `glance_wd/` | Working files of the comparison tool, kept only with `--keep-files`. |
| `glance_wd/pairwise_ani.csv` | `--tool sourmash`: columns `genome1`, `genome2`, `similarity`, one row per pair. An empty `similarity` means the pair is below the sketch's reach (under 30 shared hashes and a tenth of the smaller sketch). |
| `glance_wd/dendrogram_leaves.txt` | `--tool sourmash`: leaf order of the dendrogram; the signatures lie beside it. |

A run removes the previous plots once the comparison succeeds, so a plot
with no values in range is absent.

## dereplicate and the dereplication commands

| Path | Content |
|------|---------|
| `derep/` | Representative genomes and per-tool intermediates. |
| `derep/representatives/` | One genome file per cluster, the distinct values of `representative` in `clusters.tsv`. |
| `derep/clusters.tsv` | Columns `representative` and `member`, one row per genome; a representative also lists itself. |
| `derep/genome_status.tsv` | Columns `genome` and `status` (`representative`, `contained` or `fail_qc`). |
| `derep/refined.tsv` | With `dereplicate --refine` or `derep-refine`: the genomes of the refined lineages (columns below). Absent without a refinement. |
| `derep/cluster_summary.tsv` | One row per representative (columns below). Written by dereplicate and `cluster-summary`. |
| `genome_n50.tsv` | In the output of each `dereplicate-chunk` and `dereplicate-merge` step (Nextflow `DEREP_CHUNK`, `DEREP_MERGE`): `genome` and `n50` (bases) of every scored genome the step read. Absent from `dereplicate-refine` output. |
| `derep/unpacked/<representative>/` | `derep-unpack`: the cluster's members, and the representative unless `--no-representant`, hard-linked or copied from `genomes/`. Replaced on each run. |
| `derep/stock/<name>/` | `derep-stock --action pack`: a stored dereplication (contents under [Storing and restoring with derep-stock](#storing-and-restoring-with-derep-stock)). |

An `unpacked/` directory is named by the representative's file name without
its genome extension; the gzip suffixes `.fasta.gz`, `.fna.gz` and `.fa.gz`
are removed as a whole. When two representatives would share a directory
(`x.fasta` and `x.fna`, or names differing only in case), the full file name
is used.

`refined.tsv`:

| Column | Meaning |
|--------|---------|
| `genome` | Genome of the lineage. |
| `seed` | Seed whose coarse cluster brought the genome in. |
| `coarse_representative` | That cluster's representative before the refinement. |
| `representative` | Representative after the refinement; blank for `fail_qc`. |
| `status` | Final status. |
| `ani` | Secondary ANI of the lineage's refinement. |

### Cluster summary

`derep/cluster_summary.tsv` has one row per representative, largest cluster
first, then by representative name. The outgroup has no row. In the
Nextflow steps the taxonomy, quality and GTDB flag come from
`selection.tsv` (`--selection-tsv`) instead of the manifest.

| Column | Meaning |
|--------|---------|
| `representative` | Keeper file name, as in `clusters.tsv`. |
| `n_members` | Genomes contained under the keeper, the keeper not counted. |
| `n_species` | Distinct species (genus and epithet) across keeper and members, from the manifest taxonomy or else the canonical file name. A genome with neither adds none. |
| `species` | Those species, comma-separated: the keeper's first, then by genome count and name; at most five, then `+N more`. An epithet shared by two genera is written with its genus. |
| `rep_completeness`, `rep_contamination` | CheckM values of the keeper; blank when unknown. |
| `member_max_completeness`, `member_min_contamination` | Best values among the scored members; blank when none is scored. |
| `best_member` | Genome the quality keeper chooses in the cluster ([keeper score](usage.md#the-keeper-score)); blank when nothing is scored. |
| `n_genomes` | Genomes in the cluster, keeper included (`n_members` + 1). |
| `rep_n50` | Keeper N50 in bases from its FASTA; filled only for a keeper with completeness and contamination whose file is found. |
| `best_score` | Keeper score of `best_member`, four decimals; blank with `best_member`. |
| `rep_is_gtdb_representative` | 1 when the keeper is a GTDB species representative, else 0. |
| `refined` | 1 for a cluster from a [refinement](usage.md#finer-resolution-inside-a-lineage), read from `refined.tsv`, else 0. The last column. |

Without quality values the four quality columns, `best_member` and
`best_score` are blank and the command logs it. `best_member` differs from
`representative` under `--keeper tool`, after `--reduce` or the Nextflow
merge (only representatives are compared), and under `--keeper gtdb` where a
GTDB representative was kept over a better-scored member. Listing a
cluster's members is described in
[usage.md](usage.md#finding-the-members-of-a-cluster).

### Storing and restoring with derep-stock

`derep/stock/<name>/` holds `clusters.tsv`, `genome_status.tsv`,
`cluster_summary.tsv`, `refined.tsv` when the run holds a refinement,
`representatives/` (links to the representative genomes) and `record.json`,
the completed `dereplicate` record at pack time (`tool`, `params`,
`tool_versions`, `completed`; absent when there was none). The workflow is
in [usage.md](usage.md#storing-and-restoring-a-dereplication).

`derep-stock --action unpack` changes the working directory as follows:

- `derep/cluster_summary.tsv` is rebuilt from the restored clusters and the
  current manifest; the stored summary keeps the quality of pack time. A
  scored genome no longer in `genomes/` has no N50 term in its score and a
  blank `rep_n50`; unpack and `cluster-summary` name it in a warning.
- A live `genome_status.tsv` the stored run lacks is removed.
- `derep/representatives/` is replaced as a whole, from `genomes/` by name.
- The `dereplicate` record takes `tool`, `params` and `tool_versions` from
  `record.json`, adds `stock: <name>`, and is completed at unpack time. It
  holds no resume fingerprint.
- Without `record.json`, the record keeps the tool and parameters of the
  `dereplicate` record current at unpack time; a record left by an
  interrupted `dereplicate` then carries no tool, and one left by an
  interrupted unpack keeps what that unpack carried.

## annotate

| Path | Content |
|------|---------|
| `annotations/<record>/<record>.gff3`, `.faa`, `.ffn` | One directory per annotated record (the representatives, or every genome with `--all-genomes`, and the outgroup): GFF3 with `##FASTA`, proteins and gene nucleotides from bakta or prokka (prokka's `.gff` copied as `.gff3`). |
| `annotations/<record>/annotation.json` | `annotator`, `annotator_version`, `db_version` and `extra` of the record; `annotate-index` indexes the record with this annotator, else with its `--tool`. |
| `annotations/annotations.tsv` | Index, one row per record (columns below); rewritten after every annotated genome. Written by annotate and `annotate-index`. |

A record is annotated again only when its genome, annotator, tool version,
database version or tuning changed. A record whose genome left the set is
kept; one whose genome left the working directory is removed.

`annotations.tsv`:

| Column | Meaning |
|--------|---------|
| `record_name` | Record name. |
| `set` | `genomes` or `outgroup`. |
| `genome_path`, `genome_digest` | Genome path (relative to the working directory, or absolute outside it) and its SHA-256. |
| `annotator`, `annotator_version` | bakta or prokka, and its version. |
| `db_version` | bakta database version; empty for prokka. |
| `extra` | Tuning as sorted `key=value` pairs, comma-separated. |
| `gff_digest` | SHA-256 of the GFF3. |
| `completed` | UTC timestamp. |

## snptype

`snptype` writes `snp/`. The typing pass of `phylo --msa-source snptype`
writes the same files under `tree/msa/` (see [phylo](#phylo)).

| Path | Content |
|------|---------|
| `snp/core_snp.fasta` | Core-SNP (variable-site) alignment; masked in place with `--mask`. |
| `snp/full_alignment.fasta` | Whole-genome alignment in reference coordinates (snippy, parsnp, simple); the input of `--mask`. In the `simple` typer's alignment, positions where a genome has no aligned base are N. |
| `snp/snp_distance_matrix.tsv` | `simple` typer only: pairwise SNP distances over sites where both genomes have A, C, G or T; `NA` without such a site. Computed before masking. |
| `snp/variants.vcf` | The typer's variant calls (snippy); for `ska2 --tool-arg lo=true`, the indels `ska lo` recovered. |
| `snp/alignment_stats.json`, `tree/msa/alignment_stats.json` | Constant-site counts of the alignment beside it (keys below). Absent for `ska2 lo=true`. |
| `snp/masker_tree.nwk`, `tree/msa/masker_tree.nwk` | With `--mask gubbins`: Gubbins' final tree of the alignment beside it. Removed when a later run does not mask. |

`alignment_stats.json` is written when the typer produced a whole-genome
alignment, or for `ska2` an all-sites alignment counted in scratch:

| Key | Meaning |
|-----|---------|
| `v` | Format version. |
| `alignment`, `alignment_digest` | The variable-site alignment the counts belong to, and its SHA-256. |
| `source` | File the counts were taken from. |
| `constant_sites` | Constant A, C, G and T columns; added back by IQ-TREE (`-fconst`) and RAxML-NG (`+ASC_STAM`). |
| `variable_sites` | Variable columns; for `ska2` counted on the all-sites file. |
| `masked` | Whether masking followed. |
| `core_fraction` | `--core-fraction` applied; null when none. |

## phylo

### Tree and alignment

| Path | Content |
|------|---------|
| `align/msa.fasta` | Whole-genome alignment from the aligner (`--msa-source aligner`), records named by record name. |
| `align/msa_source.json`, `tree/msa/msa_source.json` | Reuse stamp of the alignment: `v`, `key` (hash of the genome set and source settings), `artifact`, `artifact_digest`, `versions` and `records`. |
| `tree/tree.nwk` | The tree, rooted as `phylo --root` asks: on the outgroup branch, at the midpoint, or not at all (`none`). |
| `tree/tree.unrooted.nwk` | The builder's output as written, when phylo rerooted it. |
| `tree/` | The builder's own files (logs, bootstrap trees); replaced on a rebuild, `msa/` kept. |
| `tree/partitions.txt` | With `--partitioned`: partition definitions given to `iqtree` (`-p`) or `raxmlng` (`--model`), one gene per line. |
| `tree/gene_trees.nwk`, `tree/aster.nwk` | `--treebuilder aster`: the gene trees ASTER summarised, one Newick per line (IQ-TREE aBayes or FastTree SH-like supports), and ASTER's own output. |
| `tree/genes/`, `tree/gene_trees/` | `--treebuilder aster`: per-gene FASTA cut from the alignment, and the gene-tree builder's files. |

### SNP-typing cache

| Path | Content |
|------|---------|
| `tree/msa/` | With `--msa-source snptype`: the alignment of phylo's own typing pass, outgroup included, and its reuse stamp. |
| `tree/msa/core_snp.fasta` | Core-SNP alignment; masked in place with `--mask`. |
| `tree/msa/full_alignment.fasta`, `snp_distance_matrix.tsv`, `variants.vcf` | As under [snptype](#snptype), when the typer writes them. |

`tree/msa/` is a reuse cache. After a switch to `--msa-source aligner` or an
alignment-free builder it may remain from an earlier run and then does not
describe the current tree.

### Core-gene files

With `--msa-source coregene`:

| Path | Content |
|------|---------|
| `tree/msa/core_gene_alignment.fasta` | Concatenated core-gene alignment. With Panaroo, of the genes in its entropy-filtered alignment unless `filtered=false`. Invariant sites are kept, so there is no `alignment_stats.json`. |
| `tree/msa/pangenome/` | The builder's tables (list below). |
| `tree/msa/core_gene_stats.json` | Summary of the alignment (keys below). |
| `tree/msa/core_genes.tsv` | One row per aligned gene (columns below); only for a builder that hands over one alignment per gene. |
| `tree/msa/partitions.json` | Columns of each kept gene (keys below); same condition. |
| `tree/msa/msa_source.json` | Reuse stamp; it includes the annotation index digest. |

A genome left as a row of gaps by the builder is removed and reported like
one the builder dropped. For a per-gene builder RepGenR2 measures, ranks
and caps the genes, and concatenates the kept ones in gene-id order with
gaps for a genome that lacks a gene. A record of only gaps or
missing-data symbols counts as absent.

Builder tables in `pangenome/`:

| Builder | Tables |
|---------|--------|
| Panaroo | `gene_presence_absence.csv`, `gene_presence_absence.Rtab`, `summary_statistics.txt`. |
| PPanGGOLiN | `gene_presence_absence.Rtab`, `matrix.csv`, `genomes_statistics.tsv`, `mean_persistent_duplication.tsv`, one family list per partition. |
| SCARAP | `genes.tsv`: gene, genome and orthogroup per row. |
| OrthoFinder | `Orthogroups.GeneCount.tsv`, and `core_orthogroups.txt`, the orthogroups kept. |
| markers | `markers_summary.tsv` and `markers_per_genome.tsv` (columns below). |

`markers_summary.tsv`, one row per genome:

| Column | Meaning |
|--------|---------|
| `genome`, `set` | Genome and marker set (bac120 or ar53). |
| `unique`, `multiple_unique`, `multiple`, `missing` | Markers per call. |
| `single_copy_percent` | Single-copy markers as a percentage of the set, one decimal. |

`markers_per_genome.tsv`, one row per genome and marker:

| Column | Meaning |
|--------|---------|
| `genome`, `set`, `marker` | Genome, marker set and marker. |
| `status` | `unique`, `multiple_unique`, `multiple` or `missing`. |
| `gene`, `bitscore` | Protein kept and its bit score; empty for `multiple` and `missing`. |

`core_genes.tsv`:

| Column | Meaning |
|--------|---------|
| `gene` | Gene identifier. |
| `kept` | Whether the gene is in the concatenation. |
| `rank` | 1-based position among the kept genes; empty for a gene left out. |
| `length` | Alignment columns. |
| `present`, `occupancy` | Genomes with at least one called character, and that count over the genomes of the set. |
| `gap_fraction` | Gap characters over all characters of the present records. |
| `variable_sites` | Columns with at least two states (gaps and missing data excluded). |
| `informative_sites` | Columns with at least two states each seen at least twice. |
| `mean_identity` | Mean fraction of identical pairs over the columns with two or more called characters. |
| `start`, `end` | 1-based inclusive columns of the gene; empty for a gene left out. |

`partitions.json`: `v` (1), `alignment`, `alignment_digest`, `alphabet`
(`dna` or `protein`; the tree builder takes its model family from it) and
`partitions`, a list of `[gene, start, end]`. The file is ignored once the
alignment's digest changes.

`core_gene_stats.json`:

| Key | Meaning |
|-----|---------|
| `alignment`, `source` | Alignment file and builder. |
| `alphabet` | `dna`, or `protein` with `--coregene-alphabet protein`. |
| `core_genes` | Genes kept; for a builder that returns only its concatenation, the builder's own count. |
| `core_genes_available` | Genes the builder aligned, after trimming. |
| `variable_sites_total`, `informative_sites_total` | Sums over the kept genes. |
| `gene_table`, `partitions` | `core_genes.tsv` and `partitions.json`. |
| `alignment_length`, `genomes` | Columns and genomes of the alignment. |
| `core_threshold`, `clustering_identity` | Core threshold and the identity the builder clustered at. |
| `selection` | `max_genes`, `rank`, `min_length`, `max_gaps`, `seed`. |
| `trim` | `none`, or `clipkit` with `--gene-trim clipkit`. |
| `trim_mode` | ClipKIT mode (`smart-gap`). |
| `trimmed_columns_removed` | Columns before minus after trimming, over all genes; a dropped gene counts all its columns. |
| `genes_dropped_by_trimming` | Genes that lost every column. |
| `clean_mode`, `filtered` | Panaroo settings; left out for a builder that ignores them. |
| `builder_settings` | Builder-specific settings (list below). |

`core_genes_available`, `variable_sites_total`, `informative_sites_total`,
`gene_table` and `partitions` are null for a builder that returns only its
own concatenation; such a builder accepts only the default selection. The
three trimming counts are null without trimming.

`builder_settings` holds `clustering_on` and:

- PPanGGOLiN: `partition`, `coverage`, `single_copy`;
- SCARAP: `method`, `seeds`, `scarap_core_m`, `core_prefilter`;
- OrthoFinder: `search`, `inflation`;
- markers: `marker_set` (used), `marker_set_requested`, `markers_in_set`,
  `min_markers`, `single_copy_markers_median`, `single_copy_markers_min`
  and `markers_db_version`.

## tree2tax

| Path | Content |
|------|---------|
| `tree2tax.tsv` | Taxonomy for FlexTaxD's `--taxonomy_file`, one edge per row (columns below). |
| `tree2tax_nodes.tsv` | One row per node, root included (columns below); not read by FlexTaxD. |
| `genomes_map.tsv` | For FlexTaxD's `--genomeid2taxid`: no header; genome accession and node, tab-separated. |

`tree2tax.tsv` has a header row; every node has one parent and no field is
empty.

| Column | Meaning |
|--------|---------|
| `child` | Node label: the genome file stem for a leaf; a hash of the descendant leaf set, or `--node-basename` plus a number, for an internal node. |
| `parent` | Parent label; the top node is `root` (`--root-name` changes it). |
| `rank` | With `--ranks` only: `strain` for a genome, `no rank` for an internal node, or a rank `--rank-red` assigned. Present on every row. |

The format follows the FlexTaxD reader on its GitHub master branch
(compatibility in [usage.md](usage.md#flextaxd-version-and-depth)).

`tree2tax_nodes.tsv`:

| Column | Meaning |
|--------|---------|
| `node`, `parent` | Node label and parent label (empty for the root). |
| `n_leaves` | Leaves below the node. |
| `depth` | Path length from the root. |
| `branch_length` | Length of the node's branch. |
| `support` | Builder support as a fraction; NA for leaves and trees without supports. |
| `red` | Relative evolutionary divergence (Parks et al. 2018, as PhyloRank computes it); NA without branch lengths or a root (below). |
| `lineage` | Family;genus;species of a leaf's genome from the manifest; empty otherwise and on the stateless step. |

`red` is 0 at the root, or at the ingroup's common ancestor when the
outgroup was removed, and 1 at the leaves. With the outgroup kept the
ingroup values depend on where the root sits on the outgroup branch. It is
NA on a tree without an outgroup or midpoint root, which includes the
stateless step without `--keep-root` or an outgroup.

`genomes_map.tsv` has one row per representative (its own leaf), per
dereplicated member with `--include-dereplicated` (its representative's
leaf), and, for a grouped viral isolate, per member segment accession.

## Machine-readable output

`repgenr status --json`, `doctor --json`, `census --json` and
`phylo-advise --json` print one JSON object on stdout; log messages go to
stderr. The `schema` key names the format and version, and removing or
renaming a key raises the version. Exit codes are as without `--json`; on
exit 3 (a missing working directory, or for `status` a malformed
`repgenr.yaml`) stdout is empty. The `level`, `area`, `state`, `pipeline`
and `next` values are stable; `message`, `reason`, `detail` and `notes` are
prose.

### status and doctor

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

- `pipeline`: `bacterial`, `viral`, `local` or `reads`. Null without
  `repgenr.yaml` or a recorded stage (`stages` empty, `notes` names the
  entry stages), and when no entry stage is recorded (`stages` then follows
  dereplicate, phylo and tree2tax).
- `stages`: the stages of the lineage in order (`in_chain` true), then any
  other recorded stage such as `snptype` or `glance` (`in_chain` false).
- `state`: `done`, `stale` (an input changed or an output is missing;
  `reason` says which), `interrupted` (started, not finished) or `pending`.
- `fingerprint`: true when the record holds a resume fingerprint. A `done`
  stage without one is recomputed by its next invocation.
- `detail`: the GTDB release or API query date of `metadata`, or the source
  working directories of `ingest --from-workdir` (then the sources whose
  record was not done); else null.
- `next`: the first stage of the lineage not done, or null. `notes` holds
  remarks on that stage.
- `unchecked`: null, or the error that stopped the staleness check; the
  states then rest on the records alone.

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

- `level`: `fail`, `warn` or `ok`; findings are ordered failures, warnings,
  the rest, each group by `area`.
- `area`: a stage name or a topic such as `config`, `genomes`, `manifest`,
  `outgroup`, `leftovers` or `summary`. A check that raised an error is a
  `fail` in an area named after the check (`tree`, `manifest_drift`), with a
  message starting "Check could not complete".
- `exit_code`: 7 when `counts.fail` is above zero, else 0.
- `quick`: true under `--quick`; genome files were not read, so a non-FASTA
  file under a FASTA name is not reported.

### Census table

`repgenr census --tsv PATH` writes one row per genus (family census) or
species (genus census), with a header and no totals line. Columns, in this
order, depend on the mode:

| Column | Mode | Meaning |
|--------|------|---------|
| `name` | all | Genus or species as the source writes it (GTDB `Bacillus_A thuringiensis_S`, NCBI `Francisella sp. LA112445`, NCBI Virus `Hepatovirus ahepa`, BV-BRC `Hepatovirus A`); filename tokens where only those are known. With `--runs`, runs-only rows end in `(NCBI)` and come last. |
| `species` | family census | Distinct species of the genus (GTDB species with `--runs`; selected species in a working directory). |
| `candidates` | selection | Genomes, sequence records or runs the entry stage chose from. |
| `genomes` | GTDB, selection | GTDB genomes; in a working directory, selected genomes present in `genomes/`. |
| `gtdb`, `sra`, `ncbi_virus`, `bvbrc`, `local` | selection | Selected genomes per manifest source; only sources present are listed. |
| `representatives` | GTDB, selection | GTDB species representatives. |
| `runs`, `biosamples`, `illumina`, `ont`, `pacbio`, `other` | `--runs` | ENA whole-genome sequencing runs, distinct biosamples, runs per platform; `other` holds every other platform, so the four add up to `runs`. |
| `sequences`, `complete` | viral | Sequence records, and those marked complete (empty for BV-BRC). |
| `isolates`, `segmented` | viral | Genomes after grouping segments by isolate (before vgenome's length filter); `yes` for two or more segment labels (empty for BV-BRC). |
| `median_length` | BV-BRC | Median of the per-taxid median sequence lengths. |
| `clusters` | selection | Dereplication clusters whose representative falls in the row. |

`census --json` (`repgenr.census/1`) holds:

- `schema` and `repgenr` (version);
- `taxon`, `mode` (`taxon`, `candidates` or `selection`) and `source`
  (`gtdb-api`, `gtdb-table`, `ncbi_virus`, `bvbrc`, `selection` or
  `selection.tsv`);
- `rank` (family or genus) and `row_rank` (genus or species);
- `totals`: genera, species and the header counts; in the selection mode
  also `sources` and `outgroup`;
- `rows`: objects with the TSV columns, an empty BV-BRC field as `null`;
- `notes`;
- `runs_taxon`: with `--runs`, the taxon the runs were counted under when it
  differs from `taxon` (`NCBI genus Bacillus` for GTDB `Bacillus_A`); else
  empty.

### phylo-advise

`repgenr phylo-advise --json` (`repgenr.phylo-advise/1`) holds `schema`,
`repgenr`, `evidence` and `recommendation`.

- `evidence`: `workdir`, `all_genomes`, `n_genomes`, `n_species` and
  `n_genera` (from canonical record names), `names_canonical`, `ani_min`
  and `ani_median` (percent, null without sketches), `ani_unknown_pairs`
  (pairs below the sketch's reach), `has_annotations`, `outgroup` and
  `notes`.
- `recommendation`: `preset`, `purpose`, `route`, `summary`, `reasons`,
  `command` (the `phylo` command with `--preset`, as a list of arguments)
  and `expanded_command` (the same with the preset written out).

## Nextflow output directory

The pipeline publishes under `--outdir` (default `results/`):

| Path | Content |
|------|---------|
| `metadata/` | `selection.tsv` and `outgroup_accession.txt`. |
| `reads/` | Reads mode: `reads.tsv`, `selection.tsv`, `assembly_stats.tsv`, `outgroup_accession.txt`, `excused_runs.tsv`, and `screen_reads.tsv` with `--assemble_args '--screen-reads'`. |
| `reads/qc/` | `quality.tsv` (`run_accession`, `completeness`, `contamination`) and `classification.tsv` (`run_accession`, `taxonomy`, `rank`, `score`, `db_version`) of the QC step. |
| `reads/reads_sketches/<variant>/` | With the reads sketch step: the sketches and `reads_sketches.tsv` ([reads-sketch](#reads-sketch)). |
| `reads/reads_screen/<variant>/` | With the reads screen: `reads_screen.tsv` and the per-run decisions ([reads-screen](#reads-screen)); `passed.txt` is not published. |
| `dereplicate/merged/` | The dereplication contract of the merge step, with `genome_n50.tsv`. |
| `dereplicate/refined/` | With a refinement: the contract of the refine step, with `refined.tsv`. |
| `sketches/` | Genome sketches of the `SKETCH` step. |
| `annotations/<record>/` | With `--msa-source coregene` in `--phylo_args`: one directory per `ANNOTATE` task, without the index ([annotate](#annotate)). |
| `phylo/` | `tree/` (tree and builder files) and the alignment under `align/` or `tree/msa/`. With `phylo_split_msa`, also `msa.fasta`, `alignment_stats.json` and `partitions.json`. |
| `tree2tax.tsv`, `tree2tax_nodes.tsv`, `genomes_map.tsv` | tree2tax output ([tree2tax](#tree2tax)). |
| `pipeline_info/` | Execution reports (below). |

`pipeline_info/` holds, per run, timestamped files:

- `execution_report_*.html`: resource usage and per-task summary;
- `execution_timeline_*.html`: task timeline;
- `execution_trace_*.txt`: machine-readable trace of every task;
- `pipeline_dag_*.html`: the workflow graph;
- `software_versions.yml`: the resolved tool versions of every process,
  collected and de-duplicated from each process's `versions.yml`.
