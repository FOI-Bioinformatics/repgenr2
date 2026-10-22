# Verification record

## Purpose and scope

This page records the evidence behind RepGenR2's tool adapters and routes:
which adapters were run against their real tools, on which data, with which
tool versions where one was recorded, and what was measured. It is meant for
users who judge how far a result can be trusted and for contributors who add
or change an adapter. All records come from one or two Apple Silicon machines
running macOS. They are observations on those machines and are not a
benchmark across platforms.

Contents:

* [Evidence levels](#evidence-levels)
* [Adapter status](#adapter-status)
* [Test machines and data sets](#test-machines-and-data-sets)
* [The live suite](#the-live-suite)
* [Measured runs](#measured-runs): [reads](#reads),
  [dereplication](#dereplication), [SNP route](#snp-route),
  [core-gene route](#core-gene-route), [gzipped inputs](#gzipped-inputs-2026-10-08)
* [Platform notes](#platform-notes)

The command audit of 2026-10-07 is kept as an audit record in
[audit/command-audit-2026-10.md](audit/command-audit-2026-10.md).

## Evidence levels

Each adapter is placed at one of three levels:

| Level | Meaning |
|---|---|
| Real tool, live suite | Run end to end against the installed tool, on real or closely related data, by a test in `tests/live/` that can be re-executed. |
| Real tool, by hand | Run end to end against the installed tool and recorded on this page, without a re-executable live test. |
| Fake binary only | Covered only by offline tests: the contract suite records the argument vector and parses canned output from a fake runner. The real tool has not been run in the records on this page. |

An adapter can be at more than one level; the table below gives the
strongest. Every built-in adapter has at least the offline tests.

## Adapter status

"Native" means the tool on the host `PATH`; "container" means through
`--container docker`, with a pinned image or one built by Wave. "Last
checked" is the date of the most recent record on this page. A tool version
is given only where one was recorded.

| Family | Adapter | Evidence level | Native | Container | Last checked | Where | Tool version recorded |
|---|---|---|---|---|---|---|---|
| dereplicator | skder | real tool, live suite | yes | yes | 2026-09-10 | [live suite](#the-live-suite): test_dereplicators, test_container_runs | not recorded |
| dereplicator | galah | real tool, live suite | yes | yes | 2026-10-08 | [live suite](#the-live-suite): test_dereplicators, test_container_runs; `--genome-info` check (#233) | galah 0.4.2 and 0.5.2 |
| dereplicator | sourmash (sparse and dense back-ends) | real tool, live suite | yes | yes | 2026-10-21 | [live suite](#the-live-suite): test_dereplicators, test_aux_commands, test_container_runs; [both back-ends by hand](#sourmash-back-ends-on-24-wolbachia-genomes-2026-10-21) | not recorded |
| dereplicator | drep | real tool, live suite | no | yes | 2026-09-10 | [live suite](#the-live-suite): test_container_runs (pinned and Wave images, `glance`) | container pin dRep 3.7.1; the conda dRep 3.4.5 on osx-arm64 fails (see [install.md](install.md#2-satellites-on-path-after-core)) |
| aligner | sibeliaz | real tool, live suite | yes | yes | 2026-10-08 | [live suite](#the-live-suite): test_steps, test_container_runs; [gzipped inputs](#gzipped-inputs-2026-10-08) | not recorded |
| aligner | progressivemauve | real tool, live suite | no | yes | 2026-10-08 | [live suite](#the-live-suite): test_container_runs, test_nextflow; [gzipped inputs](#gzipped-inputs-2026-10-08) | not recorded |
| aligner | cactus | real tool, live suite | no | yes | 2026-10-08 | [live suite](#the-live-suite): test_container_runs; [gzipped inputs](#gzipped-inputs-2026-10-08) | cactus v2.9.3 |
| SNP typer | simple | real tool, live suite | yes | yes | 2026-10-08 | [live suite](#the-live-suite): test_species_set, test_container_runs; [typing throughput](#snp-typing-throughput-2026-09-11); [gzipped inputs](#gzipped-inputs-2026-10-08) | not recorded |
| SNP typer | parsnp | real tool, live suite | yes | no | 2026-10-08 | [live suite](#the-live-suite): test_species_set (osx-64 environment), test_gzip_inputs | ParSNP 2.1.5 |
| SNP typer | ska2 | real tool, live suite | yes | yes | 2026-10-08 | [live suite](#the-live-suite): test_species_set, test_container_runs; [SNP phylogeny](#snp-phylogeny-on-the-genus-set-2026-09-10) | ska2 0.5.1 |
| SNP typer | snippy | real tool, live suite | no | yes | 2026-10-08 | test_container_runs (pinned image); [gzipped inputs](#gzipped-inputs-2026-10-08) | snippy 4.6.0 |
| masker | gubbins | real tool, live suite | yes | yes | 2026-09-10 | [live suite](#the-live-suite): test_species_set, test_container_runs; [divergence limits](#gubbins-divergence-limits-2026-09-10) | not recorded |
| tree builder | iqtree | real tool, live suite | yes | yes | 2026-10-21 | [live suite](#the-live-suite): test_species_set, test_container_runs; [core-gene runs](#gene-selection-partitions-trimming-protein-and-aster-2026-10-21) | IQ-TREE 3.1.2 and 3.1.4 (osx-arm64) |
| tree builder | fasttree | real tool, live suite | yes | yes | 2026-10-21 | [live suite](#the-live-suite): test_species_set, test_container_runs; ASTER gene trees | FastTree 2.2.0 (osx-arm64) |
| tree builder | raxmlng | real tool, live suite | yes | yes | 2026-09-10 | [live suite](#the-live-suite): test_species_set, test_container_runs | not recorded |
| tree builder | mashtree | real tool, live suite | yes | yes | 2026-10-08 | [live suite](#the-live-suite): test_treebuilders_offline, test_species_set, test_container_runs; [genus-scale run](#genus-scale-run-2026-09-10) | not recorded |
| tree builder | sourmash | real tool, live suite | yes | yes | 2026-10-08 | [live suite](#the-live-suite): test_treebuilders_offline, test_container_runs | not recorded |
| tree builder | veryfasttree | fake binary only | no | no | none | offline contract test | none |
| tree builder | skani | fake binary only | no | no | none | offline contract test | none |
| tree builder | ska2 | fake binary only | no | no | none | offline contract test | none |
| tree builder | aster | real tool, by hand | yes | no | 2026-10-21 | [ASTER on 24 genomes](#gene-selection-partitions-trimming-protein-and-aster-2026-10-21); also a simulated 8-genome, 6-gene alignment with wastral, astral4 and astral-pro3 | ASTER 1.25.3 (wastral), 1.25.4 (astral4) |
| annotator | bakta | real tool, by hand | yes | no | 2026-10-21 | [annotation](#annotation-2026-10-15) | bakta 1.12.2, light database v6.0 (osx-arm64) |
| annotator | prokka | real tool, by hand | yes | no | 2026-10-15 | [annotation](#annotation-2026-10-15) | prokka 1.15.6 (osx-64 under Rosetta) |
| core-gene builder | panaroo | real tool, by hand | yes | no | 2026-10-15 | [Panaroo identity](#panaroo-clustering-identity-2026-10-15) | Panaroo 1.8.0 with mafft (osx-64 under Rosetta) |
| core-gene builder | ppanggolin | real tool, by hand | yes | no | 2026-10-21 | [PPanGGOLiN](#ppanggolin-2026-10-15); [Nextflow route](#core-gene-route-in-nextflow-2026-10-15) | PPanGGOLiN 2.3.2 (osx-arm64) |
| core-gene builder | scarap | real tool, by hand | yes | no | 2026-10-15 | [SCARAP](#scarap-2026-10-15) | SCARAP 1.0.2 (noarch) |
| core-gene builder | orthofinder | real tool, by hand | yes | no | 2026-10-15 | [OrthoFinder](#orthofinder-2026-10-15) | OrthoFinder 3.1.5 (noarch), DIAMOND 2.2.8 |
| core-gene builder | markers | real tool, by hand | yes | no | 2026-10-21 | [GTDB marker genes](#gtdb-marker-genes-2026-10-21) | HMMER 3.4 |
| assembler | skesa | real tool, live suite | yes | yes | 2026-10-07 | test_reads, test_container_runs; [assembly](#assembly-of-a-public-illumina-run-2026-09-14); native by hand in the reads audit (#226) | not recorded |
| assembler | shovill | real tool, live suite | no | yes | 2026-09-14 | test_container_runs (pinned image, simulated reads) | not recorded |
| assembler | flye | real tool, live suite | yes | yes | 2026-10-07 | test_container_runs (simulated 120 kb genome, 40x ONT-like reads, one contig); [polishing](#polishing-an-ont-assembly-2026-09-16); native by hand in the reads audit (#226) | not recorded |
| classifier | sourmash | real tool, by hand | yes | yes | 2026-09-14 | `genome-qc` and `assemble` on Wolbachia assemblies against the GTDB rs226 sketch, native and pinned image | not recorded |
| polisher | medaka | real tool, live suite | no | yes | 2026-09-16 | test_reads; [polishing](#polishing-an-ont-assembly-2026-09-16) | not recorded |
| polisher | racon | real tool, by hand | yes | no | 2026-10-07 | native by hand in the reads audit (#226); offline contract test | not recorded |
| quality (not an adapter family) | checkm2 | real tool, by hand | no | yes | 2026-10-07 | pinned image with the CheckM2 database in the reads audit (#226); [polishing](#polishing-an-ont-assembly-2026-09-16) | not recorded |

Notes on the table:

* drep, progressivemauve, cactus, snippy, shovill, medaka and CheckM2 have
  been run only inside containers. skder and SibeliaZ have no pinned image,
  because their BioContainer images are BusyBox-based and the GNU-only calls
  in their shell wrappers fail there; they run in Wave images.
* The veryfasttree, skani and ska2 tree builders have offline contract tests
  (`tests/unit/test_adapter_contracts.py`) but no record of a run against
  the real tool.
* The annotators and core-gene builders were run natively, each tool
  reached through `--bin-dir`; no container run is recorded.

## Test machines and data sets

| Machine | Description | Records |
|---|---|---|
| Machine A | macOS on Apple Silicon, 11 cores, 18 GB; Docker Desktop with Rosetta for amd64 images; tools from conda; several working directories on an external volume | Live suite, reads and polishing, genus-scale run, SNP route, gzipped inputs |
| Machine B | An Apple Silicon laptop running macOS; cores and memory not recorded | Core-gene route (2026-10-15 and 2026-10-21), sourmash back-ends |

The records do not state whether machines A and B are the same computer.
Amd64-only images ran under emulation on both, which lengthens their wall
times.

| Data set | Description | Used in |
|---|---|---|
| Synthetic sets | Seeded genome sets from `benchmarks/genomegen.py`, staged with `repgenr ingest`; closely related, for tools that need within-species divergence | Live suite |
| pureclone_20 | Synthetic near-identical set of 20 genomes | Gzipped inputs |
| francisella_all | Every GTDB genome of the genus Francisella: 1157 genomes, 937 of them F. tularensis | Genus-scale run, SNP route |
| Wolbachia, three genomes | Supergroup A: wMel, wRi, wHa | Core-gene route |
| Wolbachia, four genomes | Two supergroup-A genomes with wBm (supergroup D) and wOo (supergroup C) | Panaroo identity |
| Wolbachia, 24 genomes | 24 complete Wolbachia genomes from NCBI | Core-gene route, GTDB markers, sourmash back-ends |
| SRR25474756 | Mycoplasmopsis arginini, Illumina MiSeq paired, 134 MB | Assembly, Nextflow reads mode |
| SRR28800588 | Mycoplasmoides genitalium, MinION, 85 MB | Polishing, Nextflow reads mode |

## The live suite

### Running the suite

`tests/live/` holds the re-executable records: pytest tests marked `live`
that run the installed `repgenr` console script against real tools on
seeded synthetic genome sets. They are deselected by default and never run
in CI. Copy `tests/live/live.example.toml` to `live.local.toml`, point it at
the conda environments that hold each tool, and run:

```bash
conda run -n repgenr_dev --no-capture-output \
    pytest tests/live -m live --live-config tests/live/live.local.toml -ra
```

Markers `network` and `container` select the tests that need GTDB/NCBI or
Docker; `-m "live and not network and not container"` is the offline subset.
Tests whose tool is not on `PATH` are skipped, not failed. See
`tests/live/README.md`.

### Results of the last complete run

The table below is the last complete run, on 2026-09-10 on machine A. It is
rendered from the junit output by `scripts/live_report.py`:

```bash
python scripts/live_report.py /path/to/junit.xml --out docs/verification.md
```

<!-- live-results:start -->
Last run 2026-09-10 13:11 UTC: 75 passed, 0 failed, 0 errors, 0 skipped, 25 min in total.

| module | test | result | seconds |
|---|---|---|---|
| test_aux_commands | test_derep_stock_round_trip | passed | 7 |
| test_aux_commands | test_derep_unpack_with_and_without_representant | passed | 4 |
| test_aux_commands | test_doctor_passes_then_fails_on_a_corrupt_genome | passed | 4 |
| test_aux_commands | test_logging_flags_and_env | passed | 12 |
| test_aux_commands | test_second_run_skips_and_force_reruns | passed | 9 |
| test_aux_commands | test_status_and_versions | passed | 4 |
| test_container_runs | test_cactus_pinned_image | passed | 129 |
| test_container_runs | test_container_engine_podman_is_reported_when_missing | passed | 1 |
| test_container_runs | test_drep_in_a_wave_container_with_virus_extra | passed | 9 |
| test_container_runs | test_glance_drep_compare | passed | 7 |
| test_container_runs | test_native_result_is_not_reused_by_a_container_run | passed | 8 |
| test_container_runs | test_progressivemauve_pinned_image | passed | 9 |
| test_container_runs | test_sibeliaz_in_a_wave_container | passed | 19 |
| test_container_runs | test_simple_typer_multi_tool_image | passed | 23 |
| test_container_runs | test_skder_in_a_wave_container_with_cache_and_env | passed | 7 |
| test_dereplicators | test_adapter_recovers_the_synthetic_partition[galah] | passed | 1 |
| test_dereplicators | test_adapter_recovers_the_synthetic_partition[skder] | passed | 5 |
| test_dereplicators | test_adapter_recovers_the_synthetic_partition[sourmash] | passed | 3 |
| test_dereplicators | test_allow_incomplete_gates_a_missing_genome | passed | 4 |
| test_dereplicators | test_keeper_quality_promotes_the_best_scored_member | passed | 7 |
| test_dereplicators | test_process_size_runs_the_chunked_path | passed | 6 |
| test_dereplicators | test_reduce_species_keeps_one_representative_per_species | passed | 3 |
| test_dereplicators | test_representatives_dir_matches_clusters | passed | 3 |
| test_dereplicators | test_secondary_ani_sweep_changes_representative_count | passed | 6 |
| test_dereplicators | test_target_reps_lands_on_the_requested_count | passed | 5 |
| test_dereplicators | test_tool_arg_reaches_the_tool_command_line | passed | 3 |
| test_dereplicators | test_virus_flag_on_auto_tool_is_reported_when_ignored | passed | 5 |
| test_ingest_flags | test_outgroup_and_copy | passed | 0 |
| test_ingest_flags | test_selection_table_drives_taxonomy_and_subset | passed | 0 |
| test_network | test_api_genus_representatives | passed | 0 |
| test_network | test_api_species_limit_and_explicit_outgroup | passed | 0 |
| test_network | test_api_target_family_widens_the_selection | passed | 13 |
| test_network | test_genome_accession_list_only_is_a_pure_query | passed | 1 |
| test_network | test_genome_fetch_step | passed | 9 |
| test_network | test_genome_keep_files_retains_the_download_scratch | passed | 9 |
| test_network | test_run_bacterial_chain_end_to_end | passed | 15 |
| test_network | test_run_dry_run_prints_the_chain_without_network | passed | 0 |
| test_network | test_run_dry_run_reports_family_and_species_targets | passed | 0 |
| test_network | test_run_viral_chain_end_to_end | passed | 21 |
| test_network | test_tsv_nodownload_and_metadata_path_reuse_the_table | passed | 32 |
| test_network | test_tsv_source_downloads_and_parses_the_release_table | passed | 0 |
| test_network | test_vgenome_bvbrc_needs_ignore_duplicates | passed | 2 |
| test_network | test_vgenome_discard_glance_headers_keep_files | passed | 2 |
| test_network | test_vgenome_selection_flags | passed | 6 |
| test_network | test_vmetadata_bvbrc_source_and_filter | passed | 5 |
| test_network | test_vmetadata_list_targets_reaches_bvbrc | passed | 2 |
| test_network | test_vmetadata_ncbi_virus_complete_only | passed | 0 |
| test_network | test_vmetadata_released_after_and_host_narrow_the_set | passed | 11 |
| test_nextflow | test_docker_profile_with_progressivemauve | passed | 17 |
| test_nextflow | test_local_dataflow_variants[skder---treebuilder sourmash] | passed | 17 |
| test_nextflow | test_local_dataflow_variants[sourmash---treebuilder mashtree] | passed | 12 |
| test_nextflow | test_main_bacterial_test_profile | passed | 18 |
| test_nextflow | test_main_viral_mode | passed | 30 |
| test_smoke | test_offline_chain_sourmash_mashtree | passed | 5 |
| test_species_set | test_fasttree_and_raxmlng_from_snptype | passed | 125 |
| test_species_set | test_gubbins_mask_changes_the_core_alignment | passed | 102 |
| test_species_set | test_iqtree_from_snptype_with_bootstrap_and_outgroup | passed | 127 |
| test_species_set | test_parsnp_typer | passed | 26 |
| test_species_set | test_phylo_mask_gubbins | passed | 124 |
| test_species_set | test_simple_typer_all_genomes_with_explicit_reference | passed | 47 |
| test_species_set | test_simple_typer_on_representatives_only | passed | 12 |
| test_species_set | test_ska2_source_with_reference_and_allow_incomplete | passed | 277 |
| test_species_set | test_ska2_typer_and_tool_arg | passed | 4 |
| test_species_set | test_snptype_allow_incomplete | passed | 34 |
| test_species_set | test_tree2tax_workdir_flags | passed | 10 |
| test_steps | test_chunk_keeper_quality_from_selection_tsv | passed | 9 |
| test_steps | test_chunk_results_carry_the_contract_and_versions | passed | 16 |
| test_steps | test_merge_by_chunk_dir_recovers_the_partition | passed | 19 |
| test_steps | test_merge_by_chunk_fofn | passed | 13 |
| test_steps | test_phylo_build_aligner_and_snp_source_variants | passed | 23 |
| test_steps | test_phylo_build_and_tree2tax_relations_with_outgroup | passed | 1 |
| test_steps | test_tree2tax_relations_collapse_flags | passed | 1 |
| test_treebuilders_offline | test_alignment_free_builder_on_representatives[mashtree] | passed | 5 |
| test_treebuilders_offline | test_alignment_free_builder_on_representatives[sourmash] | passed | 6 |
| test_treebuilders_offline | test_all_genomes_puts_every_genome_in_the_tree | passed | 4 |
<!-- live-results:end -->

The modules `test_reads`, `test_gzip_inputs`, `test_container_signals`,
`test_sketches_live` and `test_reads_sketch_live` are not in this table.
Tests added to the other modules after 2026-09-10 are not in it either, for
example `test_nextflow.py::test_main_reads_mode`. Their results are given in
the measured runs below where this page has them. What each module covers,
with the flags it exercises, is in `docs/audit/cli-matrix.md` (the `live`
column) and `tests/live/README.md`.

## Measured runs

The runs are ordered by pipeline stage. Each states the data, the tools and
versions, the conditions and the result. Versions are given where they were
recorded.

### Reads

#### Assembly of a public Illumina run (2026-09-14)

* **Data:** SRR25474756, selected from ENA.
* **Tools and versions:** SKESA in its pinned image; version not recorded.
* **Conditions:** machine A, amd64 emulation, 4 threads, run by
  `tests/live/test_reads.py` (selection, download with md5 check, assembly).

| Step | Result | Wall time |
|---|---|---|
| reads (ENA, Entrez lineage) | 1 run, labelled Metamycoplasmataceae / Mycoplasmopsis / arginini | 2 s |
| fetch (HTTPS, md5) | 134 MB | 80 s |
| skesa | 24 contigs, 670879 bp, N50 86129, 308x | 70 s |

`test_reads_steps_assemble_a_public_run` runs the stateless steps behind the
Nextflow reads mode (`assemble-run`, `reads-gather`) on the same accession.
`test_nextflow.py::test_main_reads_mode` runs `--mode reads` on SRR25474756
(skesa) and SRR28800588 (flye) through the pinned images with the docker
backend, then sourmash dereplication, a mashtree tree and tree2tax. On
2026-09-14 its two assembly tasks took 1 min 39 s and 2 min 47 s, the whole
run under seven minutes, and the tree had two leaves.

#### Polishing an ONT assembly (2026-09-16)

* **Data:** SRR28800588.
* **Tools and versions:** Flye, medaka and CheckM2 in their pinned images;
  versions not recorded.
* **Conditions:** machine A, amd64 emulation. SRA rewrites FASTQ headers, so
  the reads named no basecaller; the adapter used ONT's bacterial R10.4.1
  model and recorded that assumption.

| Step | Result | Wall time |
|---|---|---|
| flye | 2 contigs, 593908 bp, N50 579957 | about 3 min |
| medaka | 2 contigs, 593932 bp | about 1 min |
| CheckM2, unpolished draft | 98.90 complete, 0.09 contamination | |
| CheckM2, polished | 98.90 complete, 0.09 contamination | 22 min for both |

Polishing changed the total length by 24 bp and did not change the CheckM2
scores. The Flye draft was already gene-complete at CheckM2 resolution.

### Dereplication

#### Genus-scale run (2026-09-10)

* **Data:** francisella_all.
* **Tools and versions:** skDER for dereplication, mashtree for the tree;
  versions not recorded.
* **Conditions:** machine A, one `run`, working directory on the external
  volume:

```bash
repgenr run -wd work/francisella_all -d all -l genus -tg Francisella \
    --metadata-source api --tool skder --treebuilder mashtree -t 8 \
    --keeper quality --collapse-support 0.5
```

| Stage | Result | Wall time |
|---|---|---|
| metadata (API, quality cards) | 1157 genomes, 937 of them F. tularensis; outgroup GCF_003574425.1 | 4 min |
| genome | 1157 downloaded in one batch, 2.2 GB, none missing | 2.5 min |
| dereplicate (skder) | 68 representatives; largest cluster 893 genomes, 36 singletons; quality keeper replaced 11 | 5.7 min |
| phylo (mashtree, 8 threads) | 69 leaves, branch lengths 0.0001 to 0.22 (outgroup) | 1 s |
| tree2tax | 136 relations, 1158 genomes mapped | under 1 s |
| whole run | 25 CPU-minutes, 2.9 GB on disk, `doctor` clean | 12 min |

#### sourmash back-ends on 24 Wolbachia genomes (2026-10-21)

* **Data:** the 24 Wolbachia genomes.
* **Tools and versions:** sourmash with the branchwater plugin; versions not
  recorded.
* **Conditions:** machine B, `repgenr dereplicate --tool sourmash
  --secondary-ani 0.95`. The dense back-end (`sourmash compare
  --avg-containment`) was selected by hiding the plugin behind a
  `--bin-dir` wrapper.

| Back-end | Representatives | Pairs of unknown ANI |
|---|---|---|
| sparse (branchwater) | 6 | not applicable |
| dense (`compare --avg-containment`) | 6, identical membership | 136 of 276 |

The records name the back-end in `tool_settings`.

### SNP route

#### SNP phylogeny on the genus set (2026-09-10)

* **Data:** the 68 Francisella representatives plus the outgroup.
* **Tools and versions:** mashtree; ska2 with IQ-TREE; the `simple` typer
  with Gubbins and IQ-TREE; versions not recorded.
* **Conditions:** machine A, 8 threads. `simple` maps every genome to one
  reference with minimap2 and calls haploid consensus SNPs. The reference
  was the alphabetically first genome, a different species from most of the
  set.

| Path | Wall time | Sites | Outcome |
|---|---|---|---|
| mashtree (alignment-free) | 1 s | mash sketches | 69 leaves |
| ska2 + IQ-TREE | 7.4 min | 1737 variable | 32 zero-length branches, median UFBoot 75, 15 of 66 splits shared with mashtree |
| simple + Gubbins + IQ-TREE | 87 min, then failed | 812015 core SNPs, 2054227 alignment columns | Gubbins stopped with a bus error in its recombination scan |

The genus alignment is 39% variable. Gubbins was designed for isolates of
one species, and its limits on this alignment are given in the next
section.

#### Gubbins divergence limits (2026-09-10)

* **Data:** subsets of the `simple` whole-genome alignment of the section
  above.
* **Tools and versions:** Gubbins and IQ-TREE; versions not recorded.
* **Conditions:** machine A, 8 threads.

| Subset | Genomes | Variable columns | Gubbins |
|---|---|---|---|
| F. tularensis | 9 | 19723 (1.0%) | 5 iterations in 1.9 min, 3234 recombinant regions |
| F. philomiragia | 18 | 170342 (8.3%) | 5 iterations in 4.4 min, 8427 recombinant regions |
| tularensis + one philomiragia | 10 | 271022 (13.2%) | bus error |
| tularensis + philomiragia | 27 | 345308 (16.8%) | bus error |
| whole genus | 68 | 809060 (39.4%) | bus error |

The bus error did not depend on the tree builder or the thread count. The
masker warns above 10% variable columns and refuses above 12%.

On the nine F. tularensis representatives, IQ-TREE with 1000 ultrafast
bootstrap replicates took 27 s on the 105597 unmasked variable sites and
21 s on the 19723 sites Gubbins kept. The two topologies share no internal
split (Robinson-Foulds 12 of a possible 12), with median support 98 unmasked
and 96 masked. These genomes were mapped to a reference of another species,
and 5.1% of that reference's columns vary within F. tularensis.

Gubbins needs a multi-threaded RAxML build (`raxmlHPC-PTHREADS*`) when it is
given more than one thread. The osx-arm64 conda package has none, so the
masker used IQ-TREE, which took 76 of the 87 minutes above. With
`--tool-arg gubbins_tree_builder=fasttree` and
`gubbins_first_tree_builder=rapidnj` the same work took 6 minutes.

#### SNP typing throughput (2026-09-11)

* **Data:** the 68 Francisella representatives.
* **Tools and versions:** the `simple` typer; versions not recorded.
* **Conditions:** machine A, 8 threads, working directory on the external
  volume. The two columns are the serial implementation with plain VCF
  scratch and the concurrent implementation with compressed BCF scratch,
  which is the current one.

| | Serial, plain VCF | Concurrent, compressed BCF |
|---|---|---|
| Mapping, calling and core-SNP reduction | 10.8 min | 2.8 min |
| Scratch left behind | 10 GB | 190 MB |
| Peak resident memory | not recorded | 0.4 GB |

The 190 MB that remain are the whole-genome alignment the stage publishes for
maskers.

### Core-gene route

All core-gene runs used machine B, each tool reached through `--bin-dir`,
and IQ-TREE from `envs/core.yml`.

#### Annotation (2026-10-15)

* **Data:** the three supergroup-A Wolbachia genomes.
* **Tools and versions:** bakta 1.12.2 with the light database
  (osx-arm64); prokka 1.15.6 (osx-64 under Rosetta); Panaroo 1.8.0 with
  mafft for the comparison.
* **Conditions:** 8 threads; Panaroo with `--clean-mode sensitive` and
  otherwise default settings.

| Annotator | Time per genome | Core genes (Panaroo) | Alignment columns |
|---|---|---|---|
| bakta | 38 to 47 s (1270 CDS for wMel) | 45 | 52050 |
| prokka | about 135 s | 45 | 52050 |

On these three genomes the two annotators gave the same 45 genes and branch
lengths equal to six digits. Every bakta GFF3 carried the `##FASTA` section
Panaroo reads, and a second `annotate` reused every record. Panaroo took
about 70 s and IQ-TREE 11 s, and `tree/tree.nwk` named all three genomes.

#### Panaroo clustering identity (2026-10-15)

* **Data:** the three, four and 24 Wolbachia genomes.
* **Tools and versions:** Panaroo 1.8.0 on the bakta annotations.
* **Conditions:** `-c` varied through `--aligner-arg threshold=`, other
  settings at their defaults. Counts are genes present in all genomes, from
  `gene_presence_absence_roary.csv`, before Panaroo's entropy filter (which
  is why the 0.98 run above reports 45 genes in `core_gene_stats.json`
  against 52 here). The alignment step of the 0.80 run was not carried
  through.

| `-c` | Genes in all three supergroup-A genomes |
|---|---|
| 0.98 (Panaroo's default) | 52 |
| 0.95 | 226 |
| 0.90 | 638 |
| 0.80 | 796 |

At 0.98, two supergroup-A genomes shared 1054 genes, three shared 52, and
the four genomes across supergroups shared none (548 genes in one genome,
1050 in two, 5 in three, unchanged by the cleaning mode). On the four-genome
set, `--family_threshold` 0.5 and 0.3 left the result unchanged. Over the
24 genomes (9214 clusters) 0.98 gave 0 genes at 95% occupancy, 1 at 90%, 4
at 80%, 59 at 70% and 95 at 50%, with 831 genes shared by exactly 10
genomes.

These measurements predate the adapter's automatic choice of 0.95 when the
canonical genome names span more than one species. The 24-genome result was
obtained at 0.98 and was not repeated at 0.95.

#### PPanGGOLiN (2026-10-15)

* **Data:** the three and the 24 Wolbachia genomes, bakta annotations.
* **Tools and versions:** PPanGGOLiN 2.3.2 (osx-arm64).
* **Conditions:** PPanGGOLiN defaults (MMseqs2 protein clustering at 80%
  identity and 80% coverage); through RepGenR2, `repgenr phylo
  --all-genomes --msa-source coregene --coregene ppanggolin --treebuilder
  iqtree --no-outgroup --root midpoint -t 8`.

| Set | Panaroo at 0.98, genes in all genomes | PPanGGOLiN strict core | PPanGGOLiN persistent |
|---|---|---|---|
| 3 supergroup-A genomes | 52 | 661 | 1697 (every family) |
| 24 genomes | 0 at 95% occupancy | 370 | 539 |

On three genomes PPanGGOLiN warned that the set is too small to partition
and placed every family in the persistent partition. On 24 genomes it chose
K = 5: 539 persistent, 709 shell (three subpartitions) and 3309 cloud
families, in 53 s, with 97 s for the DNA alignments of the persistent
families. The run through RepGenR2 is in the
[builder comparison](#builder-comparison-on-24-genomes-2026-10-15).

#### SCARAP (2026-10-15)

* **Data:** the three and the 24 Wolbachia genomes, bakta protein and gene
  files.
* **Tools and versions:** SCARAP 1.0.2 (noarch); for the `O-B` and `O-D`
  methods, OrthoFinder 3.1.5 with pandas 1.5.3 (Python 3.11).
* **Conditions:** SCARAP defaults (`core -f 0.95`, method `FH`, 100 seed
  genomes); through RepGenR2, the PPanGGOLiN command with `--coregene
  scarap`.

| Set | Method | Core genes | Time |
|---|---|---|---|
| 3 genomes | `FH` | 782 | `core` 114 s, `concat` 47 s |
| 3 genomes | `O-D` / `O-B` | 769 / 767 | not recorded |
| 24 genomes, direct | `FH` | 576 | `core` 159 s, `concat` 58 s |
| 24 genomes, direct | `O-D` / `O-B` | 565 / 568 | 576 s / 677 s, against 99 s for `FH` in the same session |

All 24 genomes were seeds, so the seed draw does not account for the
one-gene difference between the direct and the RepGenR2 `FH` runs; its cause
was not determined. `FH` gave 576 genes under both pandas 1.5 and pandas
2.3. The three `O-D`, `O-B` and `FH` timings of the second session were taken
beside another job.

#### OrthoFinder (2026-10-15)

* **Data:** the three and the 24 Wolbachia genomes, bakta proteomes.
* **Tools and versions:** OrthoFinder 3.1.5 (noarch) with DIAMOND 2.2.8.
* **Conditions:** directly with `-og` (orthogroups only); through RepGenR2,
  the PPanGGOLiN command with `--coregene orthofinder`.

| Run | Orthogroups single-copy in at least 95% of genomes | Time |
|---|---|---|
| direct, 3 genomes | 790 | 25 s |
| direct, 24 genomes | 584 | 343 s |

Through RepGenR2 (see the comparison below) every gene was back-translated
to codons with no length mismatch. The adapter's selection gives 584 on the
table of the direct run, so the one-orthogroup difference between 584 and
585 lies between the two OrthoFinder runs.

#### Builder comparison on 24 genomes (2026-10-15)

* **Data:** the 24 Wolbachia genomes, from the four sections above.
* **Conditions:** IQ-TREE on each builder's alignment; trees compared by
  unrooted Robinson-Foulds distance (21 non-trivial splits per tree).

| Builder | Core genes | Alignment columns | Stage time |
|---|---|---|---|
| Panaroo at 0.98 | 0 at 95% occupancy | none | refused |
| PPanGGOLiN (soft core 0.95, single copy) | 398 | 419047 | 186 s (221 s on a first run of the same command); PPanGGOLiN 54 s |
| SCARAP `FH` (core 0.95) | 577 | 618594 | 388 s |
| SCARAP `O-D` | 564 | 604140 | 525 s; core-gene step 380 s |
| OrthoFinder | 585 | 619944 | 515 s; search 293 s, 585 mafft alignments 62 s |

| Splits shared | PPanGGOLiN | SCARAP `FH` | SCARAP `O-D` | OrthoFinder |
|---|---|---|---|---|
| PPanGGOLiN | | 20 | 20 | 20 |
| SCARAP `FH` | | | 21 (identical) | 19 |
| SCARAP `O-D` | | | | 19 |

PPanGGOLiN and SCARAP `FH` differ by one split (Robinson-Foulds distance 2,
normalised 0.048).

#### Core-gene route in Nextflow (2026-10-15)

* **Data:** the three supergroup-A Wolbachia genomes.
* **Tools and versions:** bakta with the light database, PPanGGOLiN 2.3.2,
  IQ-TREE.
* **Conditions:** a temporary harness wiring `DEREPLICATE_SCATTER`,
  `ANNOTATE_SCATTER` and `PHYLO`, with `--phylo_args '--msa-source coregene
  --coregene ppanggolin --treebuilder iqtree --no-outgroup --root midpoint'`
  and the satellite environments given through `--repgenr_opts '--bin-dir
  ...'`.

Three `ANNOTATE` tasks ran in parallel, and `PHYLO` indexed the staged
records with `annotate-index` and read them through `--annotations-dir`.
`PHYLO` built 611 genes and 615353 columns and a three-leaf tree in about
one minute. A fresh run without `-resume`, with all six tasks executed, gave
the same 611 genes.

#### Gene selection, partitions, trimming, protein and ASTER (2026-10-21)

* **Data:** the 24 Wolbachia genomes, downloaded from NCBI with `datasets`
  by accession, annotated with `annotate-build` (bakta, light database v6.0,
  three genomes at a time with three threads each, 68 to 119 s per genome,
  1270 CDS for wMel) and indexed with `annotate-index`.
* **Tools and versions:** bakta 1.12.2, PPanGGOLiN 2.3.2, IQ-TREE 3.1.2
  (osx-arm64), ClipKIT 2.14.0, ASTER 1.25 (wastral 1.25.3, astral4 1.25.4).
  Trees compared with dendropy 5.0.8 (unrooted Robinson-Foulds distance on
  the 24 leaves, 21 non-trivial splits per tree).
* **Conditions:** one run at a time, each `repgenr phylo-build
  --genomes-dir ... --annotations-dir ... --msa-source coregene --coregene
  ppanggolin --no-outgroup -t 8` with the flags below, and `--root midpoint`
  except for ASTER (`--root none`).

| Run | Flags added | IQ-TREE model | Genes kept of 398 | Columns | Wall time | RF to run 1 |
|---|---|---|---|---|---|---|
| 1 baseline | `--treebuilder iqtree --aligner-arg model=GTR+G -B 1000` | GTR+G | 398 | 419047 | 146 s | 0 |
| 2 capped | run 1 and `--max-core-genes 100 --gene-rank informative --partitioned` | GTR+G per partition | 100 | 192041 | 195 s | 2 |
| 3 trimmed | run 2 and `--gene-trim clipkit` | GTR+G per partition | 100 | 189408 | 198 s | 2 |
| 4 protein, partitioned | `--treebuilder iqtree -B 1000 --coregene-alphabet protein --partitioned`, capped at 100 genes | LG+G per partition | 100 | not recorded | 469 s | 2 |
| 4c protein, one model | `--treebuilder iqtree -B 1000 --coregene-alphabet protein --aligner-arg model=LG+G` | LG+G | 398 | 139275 | 345 s | 0 |
| 5 ASTER, wastral | run 2 with `--treebuilder aster --root none` | gene trees with `--fast` | 100 | 192041 | 49 s | 0 |
| 5b ASTER, astral4 | run 5 and `--aligner-arg program=astral4` | gene trees with `--fast` | 100 | 192041 | 49 s | 2 |
| 6 preset | `--preset genus` | GTR+F per partition | 300 | 379134 | 211 s | 2 |
| 6b written out | the genus flags without `--preset` | GTR+F per partition | 300 | 379134 | 221 s | 2 |

A first run of row 4 with 398 partitions failed in IQ-TREE on the partition
file. Row 4 is its rerun after #56 (repgenr2), which writes NEXUS charsets
and passes the model as `-m`. Rows 6 and 6b predate #56: their partitions
received IQ-TREE's `DNA` default (GTR+F, no rate heterogeneity), whereas the
current software passes GTR+G, and they were not rerun.

* **Baseline:** the same genes and columns as on 2026-10-15. PPanGGOLiN
  took 14 s, its alignments and the selection 34 s, IQ-TREE 93 s.
* **Selection and partitions:** ranked by parsimony-informative sites, the
  100 kept genes hold 53610 of the 115436 informative sites of the full set.
  The kept rows of `core_genes.tsv`, the entries of `partitions.json` and
  `core_genes` in `core_gene_stats.json` agree in every run, and the last
  partition ends at the alignment length. The total branch length of run 2
  is 0.914 against 0.869 for the baseline.
* **Trimming:** ClipKIT (smart-gap) removed 6376 of 419047 columns over all
  398 genes in 15 s and dropped no gene. The ranking chose the same 100
  genes as run 2, and the tree has the topology of run 2.
* **Protein:** run 4c has the baseline topology. Its total branch length is
  0.590 substitutions per site against 0.869 for DNA, and over the 21
  internal branches the median ratio of protein to DNA length is 0.90.
* **ASTER:** run 5 wrote 100 gene trees (about 0.1 s each) to
  `gene_trees.nwk`, one per entry of `partitions.json`. The wastral tree has
  the baseline topology, with branch lengths in coalescent units on internal
  branches only. The astral4 tree, with lengths in substitutions per site,
  differs from both by one split. `-B 1000` was ignored with a warning.
* **Preset:** `--preset genus` and the same flags written out gave the same
  IQ-TREE command, the same 300 genes in the same order and the same
  topology (log-likelihoods -1463428.03 and -1463428.10). The alignments
  differ in one gene, where PPanGGOLiN's alignment placed one base one
  column apart in 19 genomes. The tree length of run 6 is 0.663 against
  0.869 for the baseline.

Not run: `--preset family`, RAxML-NG with partitions, the other core-gene
builders, the Nextflow route, and a comparison of support values.

#### GTDB marker genes (2026-10-21)

* **Data:** the 24 Wolbachia genomes and bakta annotations of the section
  above; the marker database from `repgenr markers-db download` (InterPro
  models, bac120 and ar53).
* **Tools and versions:** HMMER 3.4, IQ-TREE 3.1.2.
* **Conditions:** machine B, `repgenr phylo-build --genomes-dir ...
  --annotations-dir ... --msa-source coregene --coregene markers
  --markers-db ... --treebuilder iqtree --partitioned -B 1000 --no-outgroup
  --root midpoint -t 8`.

| Quantity | Result |
|---|---|
| Marker set (`--marker-set auto`) | bac120; all 24 genomes above its minimum of 12 |
| Single-copy markers per genome | 105 to 111 of 120 (median 110) |
| Markers in no genome / in two genomes | 6 / 1 (TIGR00964) |
| Markers kept (single-copy in at least 95% of genomes) | 99 of 114 found |
| Alignment | 32351 amino-acid columns, 10309 variable, 8442 parsimony-informative |
| Time | search 7 s, alignment 2 s, IQ-TREE (99 LG+G partitions) about 290 s, run 305 s |
| Distance to the PPanGGOLiN baseline (run 1 above) | four splits (RF 8 of 42) |

TIGR00382 had two different copies in 16 genomes and TIGR00456 in one, so
these count as absent there. Where a genome held a marker in two identical
copies (TIGR00020 in three genomes, TIGR01087 in one), the best copy was
kept. The four splits found only in the markers tree have bootstrap
supports of 45 to 93, and three of the four baseline splits it lacks have
support 100 in the baseline. The tree is also four splits from the
partitioned protein tree of the 398 PPanGGOLiN genes.

### Gzipped inputs (2026-10-08)

* **Data:** gzipped genomes from `pureclone_20` (three or four `.fasta.gz`
  files) with the uncompressed files as a control; four synthetic 50 kb
  genomes for the live tests.
* **Tools and versions:** as in the table.
* **Conditions:** machine A.

The stage passes gzipped genomes unchanged to a tool that reads gzip. For
any other tool it writes decompressed copies, named by record name, under
`scratch/phylo_inputs/` (phylo) or `scratch/inputs/` (snptype) and removes
them after the run.

| Tool | Reads gzip | Evidence |
|---|---|---|
| cactus v2.9.3 (Docker) | yes | HAL written in 125 s |
| ska2 0.5.1 | yes | `ska build` keeps both samples |
| `simple` typer | yes | minimap2 reads gzipped queries; the adapter decompresses the reference |
| mashtree | yes | leaves named without `.fasta` |
| sourmash tree builder | yes | leaves named without `.fasta` |
| progressiveMauve (pinned image) | no, receives copies | gzipped input: signal 11; plain control exit 0 |
| SibeliaZ (native) | no, receives copies | `test_gzip_inputs::test_sibeliaz_aligns_gzipped_genomes` passed; MSA records and leaves without `.fasta`, `scratch/phylo_inputs/` removed |
| ParSNP 2.1.5 | no, receives copies | `test_gzip_inputs::test_parsnp_types_gzipped_genomes` passed; `phylo-build --snptyper parsnp` on four gzipped `pureclone_20` genomes exited 0 in 26 s with four leaves |
| snippy 4.6.0 (pinned image) | no, receives copies | gzipped `--ref`: "Could not guess format", exit 2; gzipped `--ctgs`: BioPerl "not FASTA", exit 25 |

SibeliaZ on four `pureclone_20` genomes was stopped after 7 minutes with
both the gzipped and the uncompressed files. The block alignment step of the
macOS wrapper waited in `xargs` with no `spoa` process and an almost empty
MAF. This stall does not depend on gzip and is unresolved.

## Platform notes

The macOS and Apple Silicon notes are in [install.md](install.md): the
osx-64 satellites under Rosetta and the platforms of each environment
([1. Core and satellite environments](install.md#1-core-and-satellite-environments)),
the order of core and satellites on `PATH`, exFAT volumes and dRep on
osx-arm64 ([2. Satellites on PATH, after core](install.md#2-satellites-on-path-after-core)),
and Docker Desktop with Rosetta for amd64 images ([3. Containers](install.md#3-containers)).
