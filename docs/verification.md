# Adapter verification status

Status of each pluggable tool adapter against its real external binary. "Unit"
means covered by offline tests (mock/parse); "Live" means run end-to-end against
the installed tool on real or closely-related data.

Live verification used a small GTDB/NCBI **Francisella** set (real download) and
a closely-related **synthetic** set (for tools that need within-species
divergence). Tools were installed via conda/mamba on macOS (Apple Silicon).

## Re-running the checks: the live suite

The records below were collected by hand. `tests/live/` holds the
re-executable form: pytest tests marked `live` that run the installed
`repgenr` console script against real tools on seeded synthetic genome sets
(`benchmarks/genomegen.py`) staged with `repgenr ingest`. They are deselected
by default and never run in CI. Copy `tests/live/live.example.toml` to
`live.local.toml`, point it at the conda environments that hold each tool,
and run:

```bash
conda run -n repgenr_dev --no-capture-output \
    pytest tests/live -m live --live-config tests/live/live.local.toml -ra
```

Markers `network` and `container` select the tests that need GTDB/NCBI or
Docker; `-m "live and not network and not container"` is the offline subset.
Tests whose tool is not on PATH are skipped, not failed. See
`tests/live/README.md`.

## Adapter status

Where each built-in adapter has been run against its real tool, from the
live suite as it stands. "Native" means the tool on the host PATH of the
audit machine; "container" means through `--container docker` (a pinned
image or one minted by Wave). An adapter with neither has only the offline
contract test, which checks the argument vector against canned output.

| Family | Adapter | Native | Container | Live module |
|---|---|---|---|---|
| dereplicator | skder | yes | yes | test_dereplicators, test_container_runs |
| dereplicator | galah | yes | yes | test_dereplicators, test_container_runs (pinned image) |
| dereplicator | sourmash | yes | yes | test_dereplicators, test_aux_commands, test_container_runs (pinned image) |
| dereplicator | drep | no | yes | test_container_runs (pinned and Wave images; also `glance`) |
| aligner | sibeliaz | yes | yes | test_steps, test_container_runs |
| aligner | progressivemauve | no | yes | test_container_runs, test_nextflow (unpackaged on macOS) |
| aligner | cactus | no | yes | test_container_runs (pinned image) |
| SNP typer | simple | yes | yes | test_species_set, test_container_runs |
| SNP typer | parsnp | yes | no | test_species_set (osx-64 env) |
| SNP typer | ska2 | yes | yes | test_species_set, test_container_runs (pinned image) |
| SNP typer | snippy | no | yes | test_container_runs (pinned image) |
| masker | gubbins | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | iqtree | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | fasttree | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | raxmlng | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | mashtree | yes | yes | test_treebuilders_offline, test_species_set, test_container_runs (pinned image) |
| tree builder | sourmash | yes | yes | test_treebuilders_offline, test_container_runs (pinned image) |
| assembler | skesa | yes | yes | test_container_runs (pinned image, simulated reads; a `requires_binary` test covers the host); native by hand in the reads audit (#226) |
| assembler | shovill | no | yes | test_container_runs (pinned image, simulated reads) |
| assembler | flye | yes | yes | test_container_runs (pinned image, simulated 120 kb genome at 40x ONT-like reads, one contig); native by hand in the reads audit (#226) |
| classifier | sourmash | yes | yes | genome-qc and assemble on Wolbachia assemblies against the GTDB rs226 sketch, native and pinned image (2026-09-14) |
| polisher | medaka | no | yes | test_reads (SRR28800588 through the pinned image, see below); CheckM2 before/after recorded below |
| polisher | racon | yes | no | offline contract test (minimap2 and racon argv, rounds, stdout capture); native by hand in the reads audit (#226) |
| quality | checkm2 | no | yes | offline test on a canned quality report; pinned image with the CheckM2 database by hand in the reads audit (#226), because the host builds fail on macOS |

drep, progressivemauve, cactus and snippy have been verified only inside
containers. skder and SibeliaZ run in Wave-minted images only: their
BioContainer images are BusyBox-based and the GNU-only calls in their shell
wrappers fail there, which is why they carry no pinned image.

## Polishing an ONT assembly (2026-09-16)

SRR28800588 (Mycoplasmoides genitalium, MinION, 85 MB) assembled with Flye
and polished with medaka through the pinned images under amd64 emulation.
SRA rewrites FASTQ headers, so the reads named no basecaller; the adapter
used ONT's bacterial R10.4.1 model and recorded that assumption.

| Step | Result | Wall time |
|---|---|---|
| flye | 2 contigs, 593908 bp, N50 579957 | about 3 min |
| medaka (assumed bacterial model) | 2 contigs, 593932 bp | about 1 min |
| CheckM2, unpolished draft | 98.90 complete, 0.09 contamination | |
| CheckM2, polished | 98.90 complete, 0.09 contamination | 22 min for both, emulated |

On this run the Flye draft was already gene-complete at CheckM2 resolution
and polishing changed 24 bp of total length without moving the scores. The
12% contamination seen earlier on DRR351706 came from an MDA-amplified
library, which `reads --drop-selection` now excludes by default, not from a
lack of polishing. Polishing still matters for downstream SNP typing and
alignment, where indel errors are not absorbed by gene calling.

## Reads chain on a public run (2026-09-14)

`tests/live/test_reads.py` selects SRR25474756 (Mycoplasmopsis arginini,
Illumina MiSeq paired, 134 MB) from ENA, downloads and verifies both files,
and assembles them with SKESA in its pinned image under amd64 emulation on
the audit machine:

| Step | Result | Wall time |
|---|---|---|
| reads (ENA, Entrez lineage) | 1 run, labelled Metamycoplasmataceae / Mycoplasmopsis / arginini | 2 s |
| fetch (HTTPS, md5) | 134 MB | 80 s |
| skesa (4 threads, emulated) | 24 contigs, 670879 bp, N50 86129, 308x | 70 s |

The stateless steps behind the Nextflow reads mode (`assemble-run`,
`reads-gather`) are run on the same accession by
`test_reads_steps_assemble_a_public_run`; `genome-qc` has only the offline
tests with fakes, since the audit machine holds neither a CheckM2 database
nor a GTDB sketch. The Nextflow layer itself is covered by the nf-test stub
suite (`reads_select_process`, `reads_assemble_process`,
`genome_qc_process`, `reads_gather_process`, `acquire_reads`,
`reads_dataflow` and the `main.nf` reads-mode test) and, with real tools, by
`test_nextflow.py::test_main_reads_mode`, which runs `--mode reads` on
SRR25474756 (Illumina, skesa) and SRR28800588 (ONT, flye) through the pinned
images with the docker backend, then sourmash dereplication, a mashtree tree
and tree2tax (2026-09-14: two assembly tasks of 1m39s and 2m47s, the whole
run under seven minutes; two leaves in the tree). That run surfaced two
fixes: the steps resolve their paths before a containerised tool sees them,
and the container backend binds the real directory behind a symlinked
staging path.

## Results

The table below is rendered from the junit output of the last complete run
on the audit machine (Apple Silicon, Docker Desktop with Rosetta) by
`scripts/live_report.py`:

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

What each module covers, with the flags it exercises, is in
`docs/audit/cli-matrix.md` (the `live` column) and `tests/live/README.md`.

## Genus-scale run (real data, 2026-09-10)

One `run` over every GTDB genome of a genus, on the audit machine (11 cores,
18 GB), workdir on an external volume:

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

12 minutes end to end, 25 CPU-minutes, 2.9 GB on disk, `doctor` clean.

## SNP phylogeny on the genus-scale set (2026-09-10)

The 68 Francisella representatives plus the outgroup, three ways, on the same
machine. `simple` maps every genome to one reference with minimap2 and calls
haploid consensus SNPs; the reference here was the alphabetically first genome,
a different species from most of the set.

| Path | Wall time | Sites | Outcome |
|---|---|---|---|
| mashtree (alignment-free) | 1 s | mash sketches | 69 leaves |
| ska2 + IQ-TREE, 8 threads | 7.4 min | 1737 variable | 32 zero-length branches, median UFBoot 75, 15 of 66 splits shared with mashtree |
| simple + Gubbins + IQ-TREE, 8 threads | 87 min, then failed | 812015 core SNP, 2054227 alignment columns | Gubbins crashed in its recombination scan |

The genus alignment is 39% variable. Gubbins is built for isolates of one
species; its scan allocates per-SNP arrays on the thread stack and dies with a
bus error on input this diverse, whatever tree builder it uses and however many
threads. Where the boundary sits, measured on subsets of the same alignment:

| Subset | Genomes | Variable columns | Gubbins |
|---|---|---|---|
| F. tularensis | 9 | 19723 (1.0%) | 5 iterations in 1.9 min, 3234 recombinant regions |
| F. philomiragia | 18 | 170342 (8.3%) | 5 iterations in 4.4 min, 8427 recombinant regions |
| tularensis + one philomiragia | 10 | 271022 (13.2%) | bus error |
| tularensis + philomiragia | 27 | 345308 (16.8%) | bus error |
| whole genus | 68 | 809060 (39.4%) | bus error |

The masker now estimates this fraction, warns above 10%, and reports the figure
if Gubbins fails. Recombination masking belongs to a within-species run; a
genus-level set is served by mashtree or an alignment-based path.

Masking changes the answer where it does run. On the nine F. tularensis
representatives, IQ-TREE with 1000 ultrafast bootstrap replicates took 27 s on
the 105597 unmasked variable sites and 21 s on the 19723 Gubbins kept. The two
topologies share no internal split (Robinson-Foulds 12 of a possible 12), and
both carry high support: median 98 unmasked, 96 masked. Note what the sites
rest on. These genomes were mapped to the alphabetically first genome of the
whole genus, another species, and 5.1% of that reference's columns vary within
a species usually described as clonal. A within-species run wants `--reference`
pointed at a genome of that species; the stage warns when it falls back.

Gubbins also needs a multi-threaded RAxML build (`raxmlHPC-PTHREADS*`) whenever
it is given more than one thread. The osx-arm64 conda package has none, so the
masker falls back to IQ-TREE, which cost 76 of the 87 minutes above on this
alignment. `--tool-arg gubbins_tree_builder=fasttree` with
`gubbins_first_tree_builder=rapidnj` did the same work in 6 minutes.

## SNP typing throughput (2026-09-11)

The `simple` typer on the 68 Francisella representatives, 8 threads, workdir on
the external volume. The stage maps each genome against the reference, then
reduces the stacked consensuses to variable columns.

| | Serial, plain VCF | Concurrent, compressed BCF |
|---|---|---|
| Mapping, calling and core-SNP reduction | 10.8 min | 2.8 min |
| Scratch left behind | 10 GB | 190 MB |
| Peak resident memory | not recorded | 0.4 GB |

The old scratch was 8.3 GB of uncompressed pileup VCF, which twice filled a
disk during this audit. What remains is the whole-genome alignment the stage
publishes for maskers.

## Command audit (2026-10-07)

Every `repgenr` subcommand was checked against its help text, argument
validation, an offline run on small local data, resume, the error surface
and the matching offline tests. The checks used the offline test suite and the
live offline subset (`-m "live and not network and not container"`), with no
network access and no Docker. Audit data lived on an external volume. The 28
commands were audited in six batches, one per help panel. In the first table a
cell is `ok` (agrees with the documentation), `fixed` (a defect was found and
corrected), `open` (a finding was left for the maintainer), `n/a` (the check
does not apply, for example resume for a stateless step or a live test that
does not exist) or `not run` (a live test exists but needs the network or a
container). The Result column takes the worst cell, in the order open, fixed,
ok, and also counts the docs agreement check.

| Command | Help | Validation | Offline run | Resume | Errors | Live subset | Result |
|---|---|---|---|---|---|---|---|
| run | ok | fixed | fixed | ok | ok | ok | fixed |
| status | ok | ok | ok | n/a | n/a | ok | fixed |
| metadata | ok | fixed | ok | fixed | fixed | not run | open |
| genome | ok | ok | ok | not run | open | not run | open |
| vmetadata | ok | ok | ok | ok | open | not run | open |
| vgenome | open | ok | ok | ok | ok | not run | open |
| ingest | fixed | fixed | fixed | fixed | fixed | ok | fixed |
| reads | ok | fixed | ok | ok | ok | not run | fixed |
| assemble | ok | ok | fixed | ok | fixed | not run | fixed |
| dereplicate | ok | fixed | fixed | ok | ok | ok | fixed |
| snptype | ok | ok | fixed | ok | ok | n/a | fixed |
| phylo | fixed | ok | fixed | ok | ok | ok | fixed |
| tree2tax | ok | fixed | fixed | fixed | ok | n/a | fixed |
| glance | ok | fixed | fixed | ok | fixed | not run | fixed |
| cluster-summary | ok | fixed | fixed | fixed | ok | ok | fixed |
| derep-unpack | ok | ok | ok | fixed | fixed | ok | fixed |
| derep-stock | ok | fixed | fixed | fixed | fixed | ok | fixed |
| list-tools | ok | n/a | fixed | n/a | ok | ok | fixed |
| doctor | ok | ok | ok | n/a | n/a | ok | ok |
| versions | ok | fixed | ok | n/a | n/a | ok | fixed |
| genome-fetch | ok | ok | fixed | n/a | fixed | n/a | fixed |
| dereplicate-chunk | ok | ok | ok | n/a | ok | ok | ok |
| dereplicate-merge | ok | ok | ok | n/a | ok | ok | ok |
| phylo-build | fixed | ok | ok | n/a | ok | ok | fixed |
| tree2tax-relations | ok | ok | fixed | n/a | fixed | ok | fixed |
| assemble-run | ok | ok | fixed | n/a | fixed | n/a | fixed |
| genome-qc | ok | ok | ok | n/a | ok | n/a | ok |
| reads-gather | ok | ok | fixed | n/a | ok | n/a | fixed |

Defects found and fixed, with the commit that fixed each. Several commits
change shared code (the stage runner, the contract readers and the Newick
parser), so they also apply to commands other than the one named.

| Command | Defect | Commit |
|---|---|---|
| run | `--with-snptype --snptyper bogus` passed a dry run and failed the real run with exit 5; it now exits 2 at once | 2722386 |
| run | `--with-snptype --mask gubbins` was rejected with exit 2; the mask now reaches the snptype stage | d16afbd |
| run, phylo | The `--mask` help on `run` did not mention `--with-snptype`; the bootstrap help said ">=1000" although smaller values are accepted | 3c087a7 |
| status | The "no run found" hint named only `metadata` and `vmetadata`; it now names all four starting points | 3c087a7 |
| versions | A nonexistent workdir printed nothing and exited 0; it now exits 3 and writes no fragment | 3608ad9 |
| status, doctor | A nonexistent workdir exited 0 (status printed the no-run hint, doctor a warning); both now exit 3 with the shared message, fixed after the audit at the maintainer's request | follow-up |
| list-tools | A rejected version flag (sibeliaz) was recorded as the tool's error line; it is now recorded as unknown | f1905d7 |
| metadata | The tsv source named `--version` (the global flag) instead of `--gtdb-version` when the version was missing | 5bb9f0c |
| metadata | Replacing the `--metadata-path` table was skipped on resume and kept a stale selection | 4337fab |
| vmetadata | An unreachable BV-BRC FTP server gave a traceback with exit 1; it now exits 3 with a named error | cbf839b |
| docs/usage.md | Exit code 3 was described only as workdir state; it also covers a failed remote request | 8e93630 |
| docs/output.md | `scratch/` was attributed to snptype only; it is written by five stages | 2a9f38e |
| genome-fetch | The `datasets` preflight was skipped without `--versions-out`, giving a traceback; it now exits 4 | 2762f39 |
| genome-fetch | A `selection.tsv` without the needed columns raised KeyError; it now exits 3 and names the columns | faef47b |
| assemble-run, reads-gather | A `reads.tsv` without the needed columns raised KeyError; it now exits 3 and names the columns | 6e7efc6 |
| assemble, assemble-run | A failure reason spanned several lines in `excused_runs.tsv`; each excused run now stays on one line | c091e9f |
| assemble, assemble-run | With `--assembler auto` and no assembler installed, runs were excused as `unsupported_platform`; they are now excused as `assembler_not_installed` with a warning, and the stage exits 4 only when nothing can be assembled | e7afca7, 68bef93 |
| dereplicate, tree2tax | A missing workdir gave a traceback or was created silently; the stage runner now exits 3 (it also covers glance, derep-unpack, cluster-summary, derep-stock, genome, vgenome and assemble) | c04fc22 |
| dereplicate | docs/output.md placed tool intermediates under `derep/`; they are under `scratch/` | 02f44ec |
| snptype | Switching typers left the previous typer's optional outputs (`full_alignment.fasta`, `snp_distance_matrix.tsv`, `variants.vcf`) in `snp/`; they are now dropped | 1429100 |
| phylo | Switching tree builders left the previous builder's files in `tree/`; they are now cleared | d2dfbb1 |
| phylo | The SibeliaZ macOS wrapper concatenated AppleDouble `._*` files into the alignment | a09ede2 |
| phylo | The SibeliaZ macOS wrapper left thousands of empty block temp files | 82dc1c0 |
| tree2tax | Toggling `--include-dereplicated` was logged as a changed input file | 2e01791 |
| tree2tax | A dereplicated member that is also a leaf was mapped twice in `genomes_map.tsv` | 4ffd3ee |
| tree2tax, tree2tax-relations | A malformed Newick tree gave a dendropy traceback; it now exits 3 | ddf48c6 |
| tree2tax-relations | A nonexistent `--clusters` path was skipped silently; it now exits 3 | a66b8ef |
| glance | A workdir without genomes exited 4 (dRep absent) instead of 3 | 645feaa |
| glance | The stage record held no dRep version | e003775 |
| glance | docs/usage.md said glance compares the representatives; it compares all genomes | 8cab3b7 |
| glance | The plots counted each genome pair twice (dRep lists both orders), and the histogram's y axis was labelled MASH ANI instead of a pair count | 05f3cea |
| glance | One genome failed inside dRep with exit 6; it now exits 3 before dRep runs | 3c3fd84 |
| glance | `--plot-min` above `--plot-max`, or a bound outside 0-1, was accepted and removed the earlier plots; it now exits 2 | 27ddeea |
| glance | `--help` did not give the plot bound units or say that `--keep-files` keeps `glance_wd/` | 6825439 |
| glance | A missing dendrogram was not reported | 793d9d0 |
| glance | A deleted dendrogram was not rebuilt by a repeat run (no deliverable declared) | 6ba057d |
| status | An interrupted optional stage was shown as `(incomplete)` instead of `[interrupted]` | c7ba749 |
| cluster-summary | A missing workdir created a manifest, or raised an OSError traceback | 50194c6 |
| derep-stock | `pack` of a workdir without dereplication outputs stored an empty run | f84e11f |
| derep-stock | `unpack` of an incomplete stored run emptied the live representatives before failing | 78d14cc |
| cluster-summary | A genome with a non-canonical filename was counted as a species with a blank name | #219 |
| cluster-summary | After `derep-stock unpack` the live summary kept the quality of pack time, and `cluster-summary` skipped; unpack now rebuilds it from the current manifest | #218 |
| cluster-summary | A `clusters.tsv` with no clusters gave a header-only summary without a warning | #219 |
| cluster-summary | Species came from filenames only, so `ingest --selection` with non-canonical names gave wrong species, and equal epithets of two genera counted once; the summary now uses the manifest taxonomy, adds `n_genomes` and caps the species list | #219 |
| cluster-summary | A deleted `cluster_summary.tsv` was not rebuilt on resume; the deliverables check now reruns the stage (verified on the 50-genome set) | #209 |
| derep-stock | A repeat `unpack` after a new dereplication was skipped and restored nothing | fde2eb9 |
| derep-stock | `list` and `pack` on a nonexistent workdir exited 0 or created the workdir | bec185a |
| derep-stock | `list` wrote the run names to the log on stderr, so `--quiet` hid them; they now go to stdout, and delete logs the removed run | #218 |
| derep-stock | A `--name` past the file-name limit ended in an unexpected error (exit 1); names are limited to 100 characters (exit 2) | #218 |
| derep-stock | A refused pack or unpack (bad or unknown name) marked the last finished pack as interrupted and `doctor` failed | #218 |
| derep-stock | An unpack killed half-way left `status` reporting `dereplicate` as done with `phylo` next | #218 |
| derep-stock | Packing under a stored name replaced that run without notice; it now warns | #218 |
| derep-stock | Unpack copied every representative (1.9 GB at 1000 genomes); it now hardlinks like `dereplicate` | #218 |
| derep-stock | Unpack re-stamped the `dereplicate` record from the record live at unpack time, so a sourmash run restored after a skDER run was shown as `[skder]` by `status` and `versions`; a stored run now keeps its own record (`record.json`), and runs stored without one fall back as before | #220 |
| ingest | An `--outgroup` file outside `--genomes-dir` named like an ingroup genome replaced that genome without a message (49 of 50 genomes left); it now exits 2 | #222 |
| ingest | NCBI assembly filenames (`GCF_000008985.1_ASM898v1_genomic.fna`) parsed as family `GCF` with accession `genomic`, so a Datasets download kept one manifest row of many; the accession now comes from the leading GCA_/GCF_ accession | #222 |
| ingest | Two files with one accession were all staged while the manifest kept one; a selection listing one file twice staged it twice; both now exit 2 | #222 |
| ingest | Empty files and dangling links were staged and failed inside the dereplication tool; they now exit 2 before anything is staged | #222 |
| ingest | Files without a supported suffix (`x.fna.gz`, `X.FASTA`) were skipped without a message; they are now listed in a warning, and an empty source holding subdirectories says they are not searched | #222 |
| ingest | A malformed `--selection` exited 3 or crashed with exit 1, and `is_outgroup` `true` was read as 0 (the outgroup went into the ingroup); now exit 2 with file and line, and true/false, yes/no are accepted | #222 |
| ingest | A refused re-ingest marked the finished ingest as interrupted; the checks now run when the stage will run, before the record is touched, so a refused first ingest leaves no record and no empty manifest, and a skipped ingest is not checked | #222 |
| ingest | Relative paths were recorded as given, so `doctor` run from another directory reported the source as changed; they are now recorded absolute | #222 |
| ingest, metadata | `repgenr --force ingest` on an unchanged set cleared the manifest's dereplication status while `dereplicate` was skipped as up to date; an unchanged genome now keeps it | #222 |
| doctor | Every `.fasta.gz` genome was reported as not FASTA (to be deleted); gzip files are now judged by their decompressed content | #222 |
| doctor | Links left dangling by a moved source were reported as not FASTA with advice to re-run the genome stage; they are now named as dangling links with advice to re-run ingest | #222 |
| metadata | A malformed `--release` (`abc.def`) ended in an unexpected error (exit 1), and an unknown `--gtdb-version` made two failed downloads (exit 3); both now exit 2 | #224 |
| metadata | A `--metadata-path` that does not exist silently downloaded the full GTDB table instead; it now exits 2 | #224 |
| metadata | With `--limit`, the automatic outgroup could be a target genome the cap left out (r232, genus Francisella, `--limit 5`: F. guangzhouensis); it now lies outside the target taxon | #224 |
| metadata | `--outgroup-accession` naming a selected genome wrote it twice to `selection.tsv`, and one naming a target genome left out by `-d rep` or `--limit` was accepted; both now exit 2 (tsv and API) | #224 |
| metadata | `--source api` lowered GTDB suffixes (`Bacillus_A`) and kept a capitalised epithet (`-ts Tularensis`), and an unknown taxon or outgroup was an HTTP error with exit 3; names are spelled as GTDB does and unknown names exit 2 | #224 |
| metadata | A failed table download (network or checksum) was retried on the legacy layout and reported as "check release/version"; only a 404 tries the next layout and other failures name their cause | #224 |
| metadata | `--nodownload -r 232.1` reused the 232.0 table and recorded release 232.1 (table names carry the major release only); the exact release is now recorded beside the table and checked | #224 |
| metadata | Through an unreachable proxy a GTDB request waited 481 s before exit 3; a 15 s connect timeout reports it in 120 s | #224 |
| genome | A download batch made only of accessions NCBI no longer serves failed with exit 6 after three attempts; they are recorded in `missing_accessions.txt` | #224 |
| genome | Rehydrated genomes were not checked against the package's `md5sum.txt` (`datasets rehydrate` does not check it); a mismatch is now discarded and recorded missing | #224 |
| genome | An outgroup package without a FASTA completed the stage without an outgroup, an unserved outgroup exited 6, and a present outgroup was downloaded again on every run | #224 |
| genome | A genome or the outgroup deleted by hand was skipped on resume because `genomes/` was not empty; every promised file is now a deliverable | #224 |
| vmetadata | `--source bvbrc` reused any `download.fa`: a second target was recorded with the first target's sequences (picornaviridae after Hepatitis E), and a workdir switched from NCBI Virus failed with exit 1; reuse now requires the same source and target | #224 |
| vgenome | `--group-segments` concatenated every record of an isolate name (Lassa Josiah: three L and three S segments, 21 kb) and joined isolates of different species that share a name; isolates now group per species with one record per segment | #224 |
| vgenome | The grouped outgroup search admitted no segment of a two-segment virus (span midpoint plus/minus 15 percent); it uses the span widened by 15 percent, as documented | #224 |
| vmetadata, vgenome | The NCBI Virus species was the organism name, often strain-level or an earlier name (Orthohantavirus complete genomes: 107 organism names, 76 lineage species; Embecovirus 56 and 7; Junin under two names), which split one species across filename tokens, `--target-species`, the length window and the outgroup candidates; it is now the current NCBI Taxonomy species of the taxid (Maguari and Playas virus: `Orthobunyavirus maguariense`, which the report lineage nests under `cacheense`; Murutucu virus: `Orthobunyavirus maritubaense`, absent from the lineage), with the shallowest lineage binomial as the offline fallback, and an organism name given to `--target-species` also selects its species | #227 |
| vmetadata | A record filed at genus level (organism `Phlebovirus`, 6 records) took the family as its genus; it takes the genus | #227 |
| vgenome | Free-text segment labels (`M`, `M; medium`, `middle`) counted as different segments in `--group-segments`; labels are normalised, numbered labels such as `RNA 1` and `RNA 2` stay distinct, and `segments.tsv` keeps the submitted label | #227 |
| vgenome | A run without an outgroup left an earlier run's `outgroup/` file and `outgroup_accession.txt` behind (both back-ends), and mashtree was recorded as the tool when it never ran | #224 |
| dereplicate | A partial genome (a 40 percent fragment) made `--tool skder` exit 3: the adapter required both aligned fractions to pass `-af`, while skDER tests the member's own | #225 |
| dereplicate | Representatives were hardlinks of the adapter's copies (skDER output, galah's scratch copy), not of `genomes/`, so each was held twice on disk | #225 |
| dereplicate | `--virus` with a tool that does not read it (`skder --virus`, the viral example in usage.md) was dropped without a word; it now warns and names dRep | #225 |
| dereplicate | `--tool skder -sani` below 0.80 waited on a hidden skDER prompt in an interactive shell; external tools now run with stdin closed and skDER refuses the value with exit 2 | #225 |
| dereplicate | A rerun refused because a selected genome was missing (exit 3) marked the finished dereplication as interrupted and `doctor` failed | #225 |
| dereplicate | dRep without CheckM exits 0 without results; the run ended as an unexpected error (exit 1) and now exits 6 with the cause in the log | #225 |
| dereplicate | `--tool drep --virus` ran fastANI: the adapter's default S_algorithm, merged into the extras by the stage, hid the virus default ANImf | #225 |
| dereplicate | A gzipped genome in `genomes/` made `--tool drep` exit 3, since dRep names its decompressed copy | #225 |
| dereplicate | Genomes dRep's filter removed (`--tool-arg length=N`) had no status and the stage exited 3; they are now `fail_qc` | #225 |
| dereplicate | A deleted `genome_status.tsv` or `cluster_summary.tsv` was not rebuilt: `doctor` asked for a rerun and the rerun skipped | #225 |
| dereplicate | Without manifest quality, `--tool galah` kept a 40 percent fragment whose name sorted first as its cluster's representative (galah prefers the first listed genome); the genomes are now listed by descending file size, and with quality for every genome galah receives it as `--genome-info` (galah 0.4.2 and 0.5.2: the best-scored member represents each cluster under `--keeper tool`; rows are matched by name without the FASTA suffix, checked for all seven suffixes including `.fna.gz` and `.fa.gz`) | #233 |
| dereplicate | `--tool drep` needed CheckM even when the manifest had completeness and contamination for every genome; these are now passed as `--genomeInfo` | #233 |
| dereplicate | Chunked `--tool sourmash --target-reps` sketched the merge-level union again at each search step; it is now assembled from the chunk zips with `sourmash sig cat` | #233 |
| dereplicate | Chunks of one run could be scored from different sources (dRep `--genomeInfo` in a fully scored chunk, CheckM in another); the quality decision now covers the whole run | #233 |
| dereplicate-merge | `--keeper tool --reduce species` ranked the representatives of a species by `selection.tsv` quality, while the stage uses cluster size under `--keeper tool`; on clonal_50 with a scored 40 percent fragment the merge step kept the fragment and the stage the 19-member cluster | (#254) |
| dereplicate | `keeper_effective: quality` was recorded, and the no-quality warning skipped, when the manifest held quality only for genomes outside `genomes/` such as the outgroup | (#254) |
| dereplicate | With quality for part of a cluster, any scored genome replaced an unscored representative: a 40 percent fragment scored 40/0 replaced the complete genome galah had chosen by file size (skDER and galah; unchunked, `-s 25 -p 2` and the Nextflow steps). A scored genome now replaces an unscored representative only at completeness above 90 and contamination below 5; all three run modes keep the complete genome | #254 |
| dereplicate | Equal scores kept the tool's pick, so the keeper depended on the tool (of g000022 at 95/0 and g000030 at 100/1, skDER and galah kept g000030, sourmash and dRep g000022) and on chunking. Ties now go to completeness, contamination, N50, name; skDER, sourmash and galah agree unchunked, chunked and in the steps | #254 |
| dereplicate | On 30 GTDB r232 F. tularensis genomes (scores 99.44 to 100) the keeper replaced the closed one-contig GCF_000833375.1 (skDER) and GCF_000014645.1 (sourmash) with the 29-contig GCF_016603775.1 (100/0.00 against 100/0.03). With the N50 term all three tools keep the two-contig GCF_001880245.1 (100/0.01, N50 1.89 Mb, score 103.09 against 102.99 for the closed genomes and 102.58 for the draft) | #254 |
| tree2tax, doctor | A `tree.nwk` holding two concatenated trees passed the completeness check, and tree2tax used the first; both now refuse it (exit 3 in tree2tax) | #223 |
| snptype, phylo | ParSNP records kept harvesttools' names (`x.fasta`, `x.fasta.ref`); phylo accepted the tree, but tree2tax could not find the outgroup leaf and exited 3, and the Gubbins outgroup exclusion missed it. Records are now named by genome stem (verified on the 50-genome set) | #223 |
| snptype | ParSNP copied every query genome into scratch; they are now hardlinked | #223 |
| phylo | cactus MSA records kept its sample names ('.' replaced by '_'), so IQ-TREE's `-o` and tree2tax missed a versioned outgroup; records are renamed to genome stems (unit test with a fake cactus; cactus itself not run) | #223 |
| phylo, phylo-build | A tree whose leaves a tool renamed (extension, '.ref', characters replaced) passed the leaf check but not tree2tax; phylo now writes the input names back into `tree.nwk` | #223 |
| phylo | `phylo --msa-source snptype` replaced the tables the `snptype` stage wrote in `snp/`, while the `snptype` record stayed; a repeat `snptype` skipped and `doctor` reported nothing. phylo now removes that record with a warning | #223 |
| run | `--with-snptype --msa-source snptype` ran `snptype` before `phylo`, whose typing pass then replaced its tables; `snptype` now runs after `phylo` in that case | #223 |
| phylo, run | The two #223 rows above are superseded: the typing pass of `phylo --msa-source snptype` writes under `tree/msa/`, so `snp/` has one writer, the `snptype` record is kept and `run` types before `phylo` again (verified on clonal_50_clustered) | #228 |
| tree2tax | After an interrupted phylo rebuild, tree2tax used the previous tree without notice; it now warns | #223 |
| phylo | `--msa-source snptype` and `--mask` with an alignment-free builder were dropped without notice; the stage now warns | #223 |
| tree2tax, tree2tax-relations | An outgroup accession matching no file in the outgroup directory was reported as "not present among tree leaves" | #223 |
| all commands | SIGTERM to repgenr left the running tool (FastTree, live) behind; the tools are now stopped and repgenr exits 143 | #223 |
| all commands | SIGTERM stopped the tool but not the helpers it started; each tool now runs in its own process group, which SIGTERM, SIGHUP and Ctrl-C stop as a whole. Unit tests with a shell that starts `sleep`; live on the 50-genome set, SIGTERM during `snptype --mask gubbins` left neither `run_gubbins.py` nor its IQ-TREE running (250 ms, `pgrep` empty). An inherited ignored SIGHUP or SIGINT stays ignored, SIGTSTP without a terminal does not stop repgenr, and `docker run --init` lets the forwarded SIGTERM end the container (live: without it the container kept running) | #231 |
| docs | usage.md and output.md: sourmash tree units for `--collapse-length`, the simple typer's treatment of absent sequence, `snp/` written by phylo, the distance matrix computed before masking | #223 |
| assemble | A finished run was reused whatever the settings, so another `--min-contig-length`, `--assembler`, `--polisher`, `--polish-rounds` or `--tool-arg` kept the old contigs; the marker records the settings, a higher floor refilters and any other change assembles again | #226 |
| assemble | A resumed run recorded no tool versions and the record named the tool `auto`; reused markers supply their versions and the record names the assemblers used | #226 |
| assemble | An unreadable `assembly.ok` gave a JSON traceback; the marker is written atomically and an unreadable one is assembled again; a run without contigs above the floor left its reads in scratch | #226 |
| assemble | A run killed during assembly downloaded its reads again on resume; files matching their checksum are kept (on exFAT the clean-up also failed with Errno 2 on AppleDouble files) | #226 |
| assemble | ONT runs labelled PAIRED by ENA (359 bacterial WGS runs) were excused as `unsupported_platform`; racon polished a multi-file run with its first file and medaka refused it after assembly | #226 |
| assemble, genome-qc | A wrong database path or an absent checkm2 or sourmash binary was found only after every assembly; both are checked before any download (exit 2 or 4) | #226 |
| assemble | A rerun over finished runs refused on a nearly full disk although it downloads nothing | #226 |
| assemble | A CheckM2 `qc_failed` excuse was not logged; it is warned about, and failed downloads are named with how to retry them (`--force`) | #226 |
| reads | Runs found by accession were taken whatever their strategy, so an RNA-Seq run (SRR24576250) was selected as a genome; they now pass the WGS genomic filter | #226 |
| ingest | Only `.fasta.gz` was accepted among compressed suffixes, so `.fna.gz` (the NCBI FTP default) and `.fa.gz` were skipped with a warning; both are genome suffixes in every stage, staged compressed like `.fasta.gz` | #229 |
| metadata | `--source api` recorded `release: null` and nothing else to date the taxonomy; the record holds `api_query_date` (UTC), shown by `status` and `versions` next to the table path's release | #229 |
| genome, vmetadata | On a blocked network `datasets` made three attempts of about 8.5 minutes each before exit 6; a request to `api.ncbi.nlm.nih.gov` with a 15 s connect timeout now stops the stage with exit 3 before `datasets` runs | #229 |
| assemble | Changing only the CheckM2 gate (`--min-completeness`, `--max-contamination`) ran CheckM2 again over every assembly (about 5 min for two genomes under emulation); the scores are stored per run, keyed by the contigs' SHA-256, the database's path, size and modification time and the CheckM2 version, and reapplied | #230 |
| assemble, reads-gather | An NCBI genus that GTDB names differently in the same family with the same species epithet (Mycoplasmopsis arginini, GTDB Metamycoplasma arginini) was flagged `classifier_disagrees`; it is now flagged `genus_renamed`, still warned about, and not counted as a disagreement; a shared epithet in another family stays a disagreement | #230 |
| assemble, genome-qc | Concurrent sourmash gathers were bounded by the threads only (about 0.6 GB each); they are also bounded by `--memory-gb`, and the log names the number chosen | #230 |
| assemble | A short-read run labelled PAIRED with one FASTQ file was downloaded and then excused by `--assembler shovill`; it is planned as single-end and excused as `unsupported_layout` before the download | #230 |
| assemble | When every run was excused the stage exited 3, but `genomes/` and `selection.tsv` of the previous call remained and `dereplicate` ran on them without a warning; a set an earlier `assemble` call wrote is now emptied when at least one run was judged, and `dereplicate` exits 3; a set from another stage, or one kept because every download failed, stays | #230 |
| status, doctor, versions | A malformed `repgenr.yaml` (unparsable, a list, params that are not a mapping) gave a raw traceback with exit 1 from every command; `doctor` now reports a config failure and the others exit 3 | #234 |
| doctor | A `cluster_summary` run killed on its first attempt was listed as interrupted by `status` but passed `doctor`; both now treat any record without a completion stamp as interrupted | #234 |
| ingest, dereplicate | A genome deleted from an ingest workdir, or a representative deleted from `derep/representatives/`, made `doctor` ask for a rerun and `dereplicate` or `phylo` refuse, but the rerun skipped; each selected genome and each representative is now a deliverable | #234 |
| dereplicate | `representatives/` emptied but holding Finder's `.DS_Store` counted as present, so `dereplicate` skipped; dotfiles no longer count | #234 |
| all stages | A tool missing at preflight (exit 4) left a first run `[interrupted]`, and any refusal of a re-run left the finished record `[interrupted]` with its outputs intact (`phylo --treebuilder raxmlng` without progressiveMauve); a refusal that changed no deliverable now leaves the record as it was | #234 |
| doctor | An emptied `tree2tax.tsv` or truncated `genomes_map.tsv` passed; the header and the leaf sets of both tables are now checked | #234 |
| doctor | A genome file under `genomes/` not in `selection.tsv` was reported only as a changed input; it is now named | #234 |
| doctor | Advice named the genome and metadata stages in ingest, vgenome and assemble workdirs; the integrity guards' refusal text reached the console unformatted; exFAT `._` companions counted as leftovers; an emptied `repgenr.yaml` beside outputs passed | #234 |
| status | Every recorded stage showed `[done]` after an input changed or an output was deleted (a 3-leaf tree over 32 representatives gave "All stages complete"); such stages are now `[stale]` with the reason, from the checks `doctor` uses, and `Next:` points at the first stage that is not done | #234 |
| status, phylo | `status` said "Next: repgenr phylo" with two representatives, which phylo refuses; it now adds a note. phylo's refusal advised a lower ANI threshold, which leaves fewer representatives; it now names a higher `--secondary-ani` | #234 |
| snptype, phylo | The `simple` typer filled sequence a genome lacks with the reference base, so a copy of a genome with 500 kb removed differed from that genome at 20539 sites on the 50-genome set. Uncovered reference positions are now N, a column is variable only with two or more of A, C, G, T, and distances count sites where both genomes have a base: the pair now differs at 0 sites. Mapping uses minimap2 `asm20`; on three genomes of that set the distances are within one site of the true substitution counts (default settings: 8 to 111 sites off) | #232 |
| snptype, phylo | With a gzipped reference genome the `simple` typer exited 1 (UnicodeDecodeError); three `.fasta.gz` genomes of the 50-genome set now type to the same distances as the uncompressed files | #232 |
| snptype, phylo | Review of #232: an unaligned genome became all N with distance 0 to every genome, and IQ-TREE refuses an all-N sequence; it is now refused by the typer (exit 3), coverage is logged, and pairs without a shared site are `NA`. Deletions in the genome are N. Gubbins' default `--filter-percentage 25` leaves mostly-N taxa out of its analysis while writing them, unmasked, to its outputs (live: a genome with 1 Mb of 2 Mb missing); repgenr passes 100, and a user-set value that would exclude a taxon exits 3 before Gubbins runs. Gzipped genomes are named without `.fasta.gz` in alignments and trees | #232 |
| list-tools, versions | `cactus-pangenome --version` in a broken environment printed a Python traceback, and the interpreter path in it (`python3.12`) was recorded as cactus 3.12.0, which passed the 2.5 floor; a crashed query and numbers inside longer tokens are no longer read as versions | #235 |
| dereplicate, list-tools | Under `--container` without `--wave`, a conda spec counted as an image, so `--tool auto` above 2000 genomes chose skder, which then ran on the host and exited 4 when absent; availability now mirrors where the tool runs | #235 |
| list-tools, all stages | With the Docker daemon down, `list-tools --check` reported every image-backed tool as ok and a stage failed with exit 6 ("sourmash failed (exit 1)"); preflight now runs `docker info` and exits 4 naming the unreachable daemon | #235 |
| all stages | A pinned image Docker could not pull failed as "sourmash failed (exit 125)"; the message now names the engine, the image and the daemon's error | #235 |
| dereplicate | install.md said to provide dRep's CheckM data through `CHECKM_DATA_PATH`, but the variable never reached the container and its directory was not bound; it is now forwarded and bound | #235 |
| versions | A tool recorded by two stages with different versions kept only the last one (dereplicate's sourmash image was hidden by phylo's host sourmash); each differing version is now listed as `tool (stage)` and an unfinished stage is named on stderr | #235 |
| list-tools | `--check` printed "broken (see list-tools)" for a plugin that failed to import, and an adapter whose preflight raised an unexpected error ended the listing with exit 1; the import error is shown and the adapter is reported on its own line | #235 |
| dereplicate | The auto-select warning about a broken plugin was printed twice per run (precheck and run) | #235 |
| global options | `--container bogus` ended in a traceback with exit 1; it now exits 2. `--platform`, `--wave`, `--container-engine` and `--container-cache` were silently ignored without `--container` (and the cache by Docker); each is named in a warning. The help called `--container-cache` the Wave cache, which is not stored there | #235 |
| assemble | Under a container backend racon's minimap2 image was not recorded | #235 |
| all stages | Under `--container docker`, a tool that ignores the SIGTERM forwarded through `--init` kept its container running after repgenr was stopped (seen with busybox); containers are now named and stopped with `docker stop` when the run is interrupted | #235 |
| snptype | `--mask gubbins` passed the default 16 threads to Gubbins, whose IQ-TREE tree builder refuses more threads than cores (11 on the audit machine); the masker now caps them at the CPUs where Gubbins runs | #235 |
| install | install.md presented `environment.yml` as the Linux single-environment route, but a linux-64 dry-run does not solve (Gubbins needs Python 3.8 to 3.10; mashtree's BioPerl chain needs zlib older than 1.3, against harvesttools and medaka); the docs now say so | #235 |
| doctor, status | `doctor` exited 1 for failures found and for an unexpected error alike, and neither command had machine-readable output; `doctor` now exits 7 for failures, and both accept `--json` | #242 |
| status, doctor | An emptied `repgenr.yaml` was read as the bacterial chain ("Next: repgenr metadata"), and a record without a resume fingerprint (an older version, or restored by `derep-stock unpack`) was shown as done with no hint that its next invocation recomputes it; status now names the entry stages or follows the shared stages, and both commands note such a record | #245 |
| snptype, phylo | Gzipped genomes, accepted by `ingest` since #229, failed in progressiveMauve (signal 11), ParSNP (UnicodeDecodeError), snippy (exit 2 and 25) and SibeliaZ (empty MAF, reported as spoa out of memory); these tools now receive decompressed copies in scratch (see [Gzipped genomes](#gzipped-genomes-in-the-aligners-and-typers-2026-10-08)) | #249 |
| install | `environment.yml` replaced by `envs/core.yml` and six satellites. conda 24.11.3 dry-run solves on 2026-10-08, all passing: core on linux-64 (35 s; samtools 1.24, dRep 4.0.2 with checkm-genome 1.2.3) and osx-arm64 (16 s; dRep 3.4.5 without CheckM); gubbins, mashtree, snippy, parsnp, checkm2 and mauve on linux-64 (19 to 28 s). Satellites are appended after core on `PATH` | #238 |
| install | Review of #238: with a configured `defaults` channel and flexible priority (Miniconda, Anaconda) the core solve ran for more than 17 minutes; every file now lists `nodefaults`, and the eight solves rerun with it all pass (core linux-64 31 s, osx-arm64 14 s; satellites 18 to 27 s). `envs/mauve.yml` and the adapter's Wave spec name `mauvealigner` with `boost-cpp=1.74.0` (1.2.0, the 2022 progressiveMauve build that runs against boost 1.74 in a linux-64 container) instead of the static 2015 `progressivemauve` build, whose boost pin had no effect, and the Java GUI of `mauve`. Under the append rule Gubbins runs core's IQ-TREE 3 (same tree on a small set); `--bin-dir` gives it its own 2.4 | #238 |
| install, all stages | `--bin-dir TOOL=DIR` and `REPGENR_BIN_DIRS`: in unit tests, a tool found only in its directory passes preflight and its subprocess sees that directory first on `PATH` while another tool's does not; a wrapper's helper is taken from the same directory; Gubbins chooses IQ-TREE from its directory with none on `PATH`; the snptype record names the directory used | #247 |

Observations left for the maintainer. None changed a documented behaviour, so
they are recorded here and not fixed.

| Area | Observation |
|---|---|
| Phylogeny | cactus renames its samples ('.' to '_'), while `tree2tax` resolves the outgroup leaf by file stem, so with a versioned accession such as `GCF_000001.1` the outgroup never matched a leaf. tree2tax warned and left the tree unrooted, and later exited 3 naming the outgroup; since the second audit pass the cactus adapter renames the alignment records back to genome stems and phylo restores renamed leaves. |
| Phylogeny | ParSNP's internal RAxML step refuses fewer than four genomes, so `--snptyper parsnp` exits 6 on a three-genome set; the three-genome check in `phylo` does not cover this. |
| SNP typing | The `simple` typer fills sequence a genome lacks with the reference base. A copy of a genome with 500 kb removed differed from that genome at 20539 sites on the 50-genome set. Masking uncovered positions with N, and comparing sites only where both genomes have a base, would change a documented behaviour and the SNP counts, so it is left as a proposal; usage.md states the limitation. |
| Exit codes | A missing tool found by the preflight (exit 4, for example `snptype --mask gubbins` without Gubbins) leaves the stage shown as `[interrupted]` in `status`, although nothing ran; the harness keeps the provisional record for exit 4 and 6 alike. |
| Environment | ParSNP reads every file in its input directory, so the AppleDouble `._*.fasta` files macOS writes on exFAT volumes break it; stage the genomes on an APFS disk. |
| Exit codes | When every assembly fails, `assemble` and `reads-gather` exit 3 and not 6; this is documented behaviour in the exit-code table of docs/usage.md, with the reasons in `excused_runs.tsv`. |
| assemble | A missing CheckM2 result is kept with a warning and not excused; this is documented behaviour in docs/usage.md, since a run CheckM2 could not score is not evidence of a poor assembly. |
| Network | The BV-BRC path uses FTPS directly and does not use HTTP proxy settings, and `vmetadata --list` needs the network whatever `--source` says. A BV-BRC group download is written to a temporary file and renamed after the size check, so an interrupted transfer leaves no partial `download.fa`. |
| cluster-summary | `read_clusters` skips the first line without checking that it is the `representative`/`member` header, so a hand-written `clusters.tsv` without a header loses its first row. |
| Tests | `--live-config <path>` into the main checkout from a worktree loads two conftest files and fails; `--live-config=<path>` works. |
| Environment | dRep 3.4.5 fails on exFAT volumes because macOS writes `._*` files into its cache; use an APFS workdir for glance and `dereplicate --tool drep`. |
| glance | Changing `--plot-min` or `--plot-max` reruns the whole dRep comparison, about 80 s for 1000 genomes, although only the plots change. |
| Environment | The host CheckM2 builds in the RepGenR and repgenr_checkm2 environments fail on macOS (multiprocessing spawn cannot pickle `Predictor.__set_up_prodigal_thread`); the pinned image works. |
| glance | For 1000 genomes the dendrogram PDF is one page about 3.9 m tall; it is readable only when zoomed. |
| glance | The genomes in `outgroup/` are not part of the comparison; glance compares `genomes/` only, as documented. |
| Resume | The skip message says "use --force to re-run", but `--force` is a global option and must come before the command (`repgenr --force glance ...`); `repgenr glance --force` exits 2 with "No such option". |
| derep-stock | A stored run keeps links to the representative files, but unpack restores by name from `genomes/`; a genome replaced under the same name since the pack is restored in its current form. |
| Environment | `status` and `doctor` on a long-running workdir (`francisella_all`) were not exercised, because that workdir was not on the audit machine. |
| ingest | Any name with four or more `_`-separated tokens is read as Family_genus_species_ACCESSION (`sample_1_run_A.fasta` gives accession `A`); documented, with `--selection` as the remedy. |
| Resume | `metadata --metadata-path`, `reads --accession-file` and `assemble --outgroup` record relative paths as given, as ingest did before this audit. |
| vgenome | On the NCBI Virus path the species is the record's organism name, which is often a strain or an older name (Orthohantavirus: `Hantaanvirus-CGAa1011`; Mammarenavirus: `Argentinian mammarenavirus` and `Mammarenavirus juninense` under one taxid). One species then splits into several species tokens, which affects `--target-species`, the median-of-medians window and the outgroup candidates. A lineage-derived binomial is proposed in the deep-audit report. |
| vgenome | NCBI Virus segment labels are not normalised (`M`, `M; medium`, `middle`), so one segment can count as two labels when an isolate mixes them. |
| dereplicate | On `mixed_1000_clustered` (20 truth clusters) skder, galah, sourmash sparse and sourmash dense all recover the truth partition at the defaults (adjusted Rand index 1.0); sparse and dense pick different representatives within clusters, as choosing-tools.md states. |
| dereplicate | Without manifest quality, galah picks a 40 percent fragment as the representative of its cluster (input order), and sourmash keeps the fragment as its own cluster (k-mer similarity counts the missing part). `--keeper quality` with manifest quality corrects the first. |
| dereplicate | `--target-reps` with `--process-size` re-sketches the union of chunk representatives at each search step, since the union changes with the threshold; on 50 genomes the search took 87 s against 22 s unchunked. |
| doctor | The first-bytes FASTA check reads every genome: 28 s for 1000 genomes on an exFAT USB volume (about 35 ms per file, not cached between runs); 8 to 16 threads gave 1.3 to 1.7 times. `status` reads no genome content and takes 0.4 s there. `doctor --quick` skips the check (#242). |
| doctor | Opening the WAL-mode manifest lets SQLite create or touch `manifest.sqlite-shm` and `-wal`; no data changes. |
| Environment | dRep 3.4.5 in the local environment fails in fastANI parsing (`read_csv` no longer accepts `delim_whitespace` in the installed pandas); ANImf (`--virus`) runs. The container pin is dRep 3.7.1. Seen again on 2026-10-08 (pandas 3.0.6; this dRep has no skani option). |
| dereplicate | The workdir stage applies the keeper once to the final clusters; the Nextflow steps apply it per chunk and after the merge. Final representatives agreed on clonal_50 with skDER, sourmash and galah at -s 25, apart from ties among equal scores. |

## Gzipped genomes in the aligners and typers (2026-10-08)

Each tool was given gzipped genomes from `pureclone_20` (three or four
`.fasta.gz` files), with the uncompressed files as a control. Before the
change:

| Tool | Reads gzip | Evidence |
|---|---|---|
| progressiveMauve (pinned BioContainer) | no | signal 11, exit 11; plain control exit 0 |
| SibeliaZ (native) | no | wrapper exits 0; TwoPaCo stops with "Can't read the input file"; the MAF has no alignment rows, and the adapter reported spoa running out of memory (exit 3) |
| cactus v2.9.3 (Docker) | yes | HAL written in 125 s; cactus decompresses `.gz` inputs itself |
| ParSNP 2.1.5 | no | UnicodeDecodeError for a gzipped reference and for gzipped queries |
| snippy 4.6.0 (pinned BioContainer) | no | gzipped `--ref`: "Could not guess format", exit 2; gzipped `--ctgs`: BioPerl "not FASTA", exit 25 |
| ska2 0.5.1 | yes | `ska build` keeps both samples |
| sourmash tree builder | yes | leaves were `x.fasta` (named by `Path.stem`, fixed in #241) |
| mashtree | yes | leaves are `x` |
| simple typer | yes | minimap2 reads gzipped queries; the adapter decompresses the reference |

After the change, the tools marked "no" are given decompressed copies by the
stage. Live checks:

| Check | Result |
|---|---|
| `test_gzip_inputs::test_sibeliaz_aligns_gzipped_genomes` (four synthetic 50 kb genomes, gzipped; `phylo-build --aligner sibeliaz`) | passed; MSA records and leaves without `.fasta`, `scratch/phylo_inputs/` removed. Without the change: the empty-MAF error |
| `test_gzip_inputs::test_parsnp_types_gzipped_genomes` (same set; `--msa-source snptype --snptyper parsnp`) | passed. Without the change: UnicodeDecodeError, exit 6 |
| `phylo-build --snptyper parsnp` on four gzipped `pureclone_20` genomes | exit 0 in 26 s; four leaves named without `.fasta`; the decompressed reference is left out of ParSNP's query directory |

SibeliaZ on four `pureclone_20` genomes was stopped after 7 minutes with the
gzipped and with the uncompressed files alike: the macOS wrapper's block
alignment step waited in `xargs` with no `spoa` process and an almost empty
MAF. The stall does not depend on gzip and is left open.

## Platform notes (macOS / Apple Silicon)

- Several tools lack osx-arm64 builds; some run via an osx-64 (Rosetta) conda env
  (e.g. parsnp/harvesttools). `mauve` is unpackaged on macOS entirely.
- `mashtree` pulls `perl-bio-samtools`, which pins an ancient samtools 0.1.x;
  modern samtools/bcftools for the `simple` SNP typer live in a separate env and
  are placed ahead on PATH.
- On exFAT/NTFS volumes, macOS `._*` AppleDouble files must be ignored (handled
  in the adapters/stages) and skDER must run on a local filesystem.
- amd64-only container tools run under emulation. Docker's QEMU emulation cannot
  run some SIMD-heavy binaries (e.g. Cactus's bundled `vg` hangs even on
  `vg version`). Set Docker Desktop to the **Apple Virtualization framework** with
  **"Use Rosetta for x86/amd64 emulation"** enabled; Rosetta runs these binaries
  correctly. On native amd64/Linux hosts no emulation is involved.
