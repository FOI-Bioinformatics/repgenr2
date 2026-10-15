# CLI matrix

Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;
`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.

33 commands, 398 flags (398 with a live test or an n/a reason, 0 pending).

## Global flags

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--version` |  | n/a: prints the version and exits | none | n/a: prints the version and exits | n/a: unit test test_version_flag_prints_version | docs/developing.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--container` |  | container.backend | choice | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | docs/usage.md, docs/developing.md, docs/install.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--container-engine` |  | container.engine | none | params.repgenr_opts | tests/live/test_container_runs.py::test_container_engine_podman_is_reported_when_missing | docs/install.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--container-cache` |  | container.cache_dir | none | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | docs/install.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--platform` |  | container.platform | none | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | docs/usage.md, docs/install.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--wave` |  | container.wave_enabled | none | params.repgenr_opts | tests/live/test_container_runs.py::test_native_result_is_not_reused_by_a_container_run | docs/developing.md, docs/install.md, docs/cli-reference.md, docs/usage.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--bin-dir` |  | n/a: per-tool host directories in core.bindirs, also read from REPGENR_BIN_DIRS (comma-separated, parsed by repgenr); unit test test_the_cli_option_and_variable_configure_the_directories | callback | params.repgenr_opts | n/a: the live suite passes its [bin_dirs] through REPGENR_BIN_DIRS (tests/live/conftest.py) | docs/usage.md, docs/install.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--force` | -f | state.force | none | n/a: resume is Nextflow's -resume in the data-channel layer | tests/live/test_aux_commands.py::test_second_run_skips_and_force_reruns | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--verbose` | -v | state.log_level | none | params.repgenr_opts | tests/live/test_aux_commands.py::test_logging_flags_and_env | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--quiet` | -q | state.log_level | none | params.repgenr_opts | tests/live/test_aux_commands.py::test_logging_flags_and_env | docs/cli-reference.md, docs/audit/cli-matrix.md |

## derep-stock

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--action` |  | DerepStockParams.action | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--name` |  | DerepStockParams.name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## cluster-summary

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the dereplicate-chunk and dereplicate-merge steps write the summary themselves | tests/live/test_aux_commands.py::test_cluster_summary_regenerates | docs/cli-reference.md, docs/output.md, docs/audit/cli-matrix.md |

## census

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: read-only workdir query; tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-family` | -tf | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--viral` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--source` |  | n/a: query command, no stage parameters | callback | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--release` | -r | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-version` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--metadata-path` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--runs` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--host` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--complete-only` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--released-after` |  | n/a: query command, no stage parameters | callback | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: network query, faked in tests/unit/test_census.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tsv` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: read-only workdir query; tests/unit/test_census.py | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--json` |  | n/a: query command, no stage parameters | none | n/a: query command; the Nextflow layer starts from a workdir or a samplesheet | n/a: read-only workdir query; tests/unit/test_census.py | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |

## derep-unpack

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--no-representant` |  | DerepUnpackParams.no_representant | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant | docs/cli-reference.md, docs/audit/cli-matrix.md |

## sketch

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the SKETCH module sketches a genome directory | tests/live/test_sketches_live.py::test_ingest_sketches_and_sketch_command | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--threads` | -t | SketchParams.threads | range | n/a: workdir command; the SKETCH module sketches a genome directory | tests/live/test_sketches_live.py::test_ingest_sketches_and_sketch_command | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## dereplicate

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | DereplicateParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_adapter_recovers_the_synthetic_partition | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--primary-ani` | -pani | DereplicateParams.primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` | -sani | DereplicateParams.secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_secondary_ani_sweep_changes_representative_count | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` | -af | DereplicateParams.aligned_fraction | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | DereplicateParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--process-size` | -s | DereplicateParams.process_size | none | params.derep_process_size | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/developing.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--num-processes` | -p | DereplicateParams.num_processes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--pre-primary-ani` |  | DereplicateParams.pre_primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--pre-secondary-ani` |  | DereplicateParams.pre_secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--reduce` |  | DereplicateParams.reduce | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_reduce_species_keeps_one_representative_per_species | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--target-reps` |  | DereplicateParams.target_reps | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_target_reps_lands_on_the_requested_count | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--virus` |  | DereplicateParams.extra | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_virus_flag_on_auto_tool_is_reported_when_ignored | docs/usage.md, docs/developing.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | DereplicateParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_tool_arg_reaches_the_tool_command_line | docs/developing.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | DereplicateParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_allow_incomplete_gates_a_missing_genome | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | DereplicateParams.keeper | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_keeper_quality_promotes_the_best_scored_member | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |

## dereplicate-chunk

dispatch: `step:repgenr.stages.derep_steps.dereplicate_chunk`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genomes-fofn` |  | ChunkParams.genomes | stage | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | ChunkParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--tool` |  | ChunkParams.tool | registry | params.derep_tool | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--primary-ani` | -pani | ChunkParams.primary_ani | unit_interval | params.derep_primary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` | -sani | ChunkParams.secondary_ani | unit_interval | params.derep_secondary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` | -af | ChunkParams.aligned_fraction | unit_interval | params.derep_aligned_fraction | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | ChunkParams.threads | range | task.cpus | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--virus` |  | ChunkParams.extra | none | params.mode | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/usage.md, docs/developing.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | ChunkParams.extra | callback | n/a: not exposed; add it to ext.args in modules.config | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/developing.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--selection-tsv` |  | ChunkParams.selection_tsv | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | ChunkParams.keeper | choice | params.derep_keeper | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--sketches-dir` |  | ChunkParams.sketches_dir | stage | module: --sketches-dir sketches when params.sketch stages sketches/ | n/a: unit tests test_sketch_consumers.py; Nextflow stub test with params.sketch | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | ChunkParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |

## dereplicate-merge

dispatch: `step:repgenr.stages.derep_steps.dereplicate_merge`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--out` | -o | MergeParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--chunk-dir` |  | MergeParams.chunk_dirs | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--chunk-fofn` |  | MergeParams.chunk_dirs | stage | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_fofn | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool` |  | MergeParams.tool | registry | params.derep_tool | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--primary-ani` | -pani | MergeParams.primary_ani | unit_interval | params.derep_primary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` | -sani | MergeParams.secondary_ani | unit_interval | params.derep_secondary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` | -af | MergeParams.aligned_fraction | unit_interval | params.derep_aligned_fraction | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | MergeParams.threads | range | task.cpus | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--virus` |  | MergeParams.extra | none | params.mode | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/usage.md, docs/developing.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | MergeParams.extra | callback | n/a: not exposed; add it to ext.args in modules.config | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/developing.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--selection-tsv` |  | MergeParams.selection_tsv | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | MergeParams.keeper | choice | params.derep_keeper | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--reduce` |  | MergeParams.reduce | choice | params.derep_reduce | tests/live/test_steps.py::test_merge_reduce_and_target_reps | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-reps` |  | MergeParams.target_reps | range | params.derep_target_reps | tests/live/test_steps.py::test_merge_reduce_and_target_reps | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--sketches-dir` |  | MergeParams.sketches_dir | stage | module: --sketches-dir sketches when params.sketch stages sketches/ | n/a: unit tests test_sketch_consumers.py; Nextflow stub test with params.sketch | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | MergeParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_fofn | docs/cli-reference.md, docs/audit/cli-matrix.md |

## doctor

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_doctor_passes_then_fails_on_a_corrupt_genome | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--quick` |  | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_doctor_and_status_json_and_quick | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--json` |  | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_doctor_and_status_json_and_quick | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |

## genome

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_api_genus_representatives | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--accession-list-only` |  | GenomeParams.accession_list_only | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_genome_accession_list_only_is_a_pure_query | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GenomeParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_genome_keep_files_retains_the_download_scratch | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--sketch` |  | GenomeParams.sketch | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: the shared sketch step is exercised live through ingest (tests/live/test_sketches_live.py) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## genome-fetch

dispatch: `step:repgenr.stages.genome_steps.genome_fetch`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--selection` |  | GenomeFetchParams.selection_tsv | stage | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | GenomeFetchParams.out_dir | none | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GenomeFetchParams.keep_files | none | n/a: not exposed by the GENOME_FETCH module | tests/live/test_network.py::test_genome_fetch_step | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | GenomeFetchParams.versions_out | none | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | docs/cli-reference.md, docs/audit/cli-matrix.md |

## glance

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | GlanceParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_glance_sourmash | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--threads` | -t | GlanceParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--plot-max` |  | GlanceParams.plot_max | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: plot bounds only change the two plots; the PDF is checked for existence (test_glance_drep_compare) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--plot-min` |  | GlanceParams.plot_min | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: plot bounds only change the two plots; the PDF is checked for existence (test_glance_drep_compare) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GlanceParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/cli-reference.md, docs/audit/cli-matrix.md |

## ingest

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--genomes-dir` |  | IngestParams.genomes_dir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--from-workdir` |  | IngestParams.from_workdirs | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline tests in tests/integration/test_ingest_from_workdir.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--strict-sources` |  | IngestParams.strict_sources | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline tests in tests/integration/test_ingest_from_workdir.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--selection` |  | IngestParams.selection | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_selection_table_drives_taxonomy_and_subset | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup` |  | IngestParams.outgroup | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_outgroup_and_copy | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--copy` |  | IngestParams.copy | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_outgroup_and_copy | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--drop-foreign` |  | IngestParams.drop_foreign | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline guard tests in tests/integration/test_ingest_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--sketch` |  | IngestParams.sketch | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_sketches_live.py::test_ingest_sketches_and_sketch_command | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## list-tools

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--check` |  | n/a: query command, no stage parameters | none | n/a: environment query; the Nextflow layer preflights inside each step | tests/live/test_aux_commands.py::test_list_tools_check_reports_every_adapter | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--strict` |  | n/a: query command, no stage parameters | callback | n/a: environment query; the Nextflow layer preflights inside each step | n/a: exit status from the per-adapter states; tests/unit/test_plugin_visibility.py | docs/cli-reference.md, docs/usage.md, docs/install.md, docs/audit/cli-matrix.md |
| `--images` |  | n/a: query command, no stage parameters | callback | n/a: environment query; the Nextflow layer preflights inside each step | n/a: needs a container engine; tests/unit/test_list_tools_labels.py | docs/cli-reference.md, docs/install.md, docs/audit/cli-matrix.md |

## metadata

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_api_genus_representatives | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--dataset` | -d | MetadataParams.dataset | choice | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--level` | -l | MetadataParams.level | choice | params.metadata_args | tests/live/test_network.py::test_api_target_family_widens_the_selection | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--source` |  | MetadataParams.source | choice | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--release` | -r | MetadataParams.release | none | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-version` |  | MetadataParams.version | none | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--target-family` | -tf | MetadataParams.target_family | none | params.metadata_args | tests/live/test_network.py::test_api_target_family_widens_the_selection | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | MetadataParams.target_genus | none | params.metadata_args | tests/live/test_network.py::test_api_genus_representatives | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | MetadataParams.target_species | none | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--metadata-path` |  | MetadataParams.metadata_path | none | params.metadata_args | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--nodownload` |  | MetadataParams.nodownload | none | params.metadata_args | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--limit` |  | MetadataParams.limit | range | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--drop-foreign` |  | MetadataParams.drop_foreign | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline guard tests in tests/unit/test_foreign_rows.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |

## phylo

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_alignment_free_builder_on_representatives | README.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_iqtree_from_snptype_with_bootstrap_and_outgroup | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--aligner` |  | PhyloParams.aligner | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--all-genomes` |  | PhyloParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_all_genomes_puts_every_genome_in_the_tree | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--root` |  | PhyloParams.root | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_alignment_free_builder_roots_at_the_midpoint | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--bootstrap` | -B | PhyloParams.bootstrap | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_iqtree_from_snptype_with_bootstrap_and_outgroup | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--reference` |  | PhyloParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--aligner-arg` |  | PhyloParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/developing.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | PhyloParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--mask` |  | PhyloParams.extra | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_phylo_mask_gubbins | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/output.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--core-fraction` |  | PhyloParams.extra | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: shares the phylo builder (wiring test); the typing pass is covered offline in tests/integration/test_snptype_stage.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | PhyloParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--coregene` |  | PhyloParams.coregene | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--core-threshold` |  | PhyloParams.extra | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--clean-mode` |  | PhyloParams.extra | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |

## phylo-build

dispatch: `step:repgenr.stages.phylo.phylo_build`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genomes-dir` |  | PhyloBuildParams.genomes_dir | stage | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | PhyloBuildParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--outgroup-dir` |  | PhyloBuildParams.outgroup_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | PhyloBuildParams.outgroup_accession | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--treebuilder` |  | PhyloBuildParams.phylo.treebuilder | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa-source` |  | PhyloBuildParams.phylo.msa_source | choice | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--aligner` |  | PhyloBuildParams.phylo.aligner | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--snptyper` |  | PhyloBuildParams.phylo.snptyper | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | PhyloBuildParams.phylo.no_outgroup | none | params.phylo_args | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--root` |  | PhyloBuildParams.phylo.root | choice | params.phylo_args | tests/live/test_treebuilders_offline.py::test_alignment_free_builder_roots_at_the_midpoint | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--bootstrap` | -B | PhyloBuildParams.phylo.bootstrap | range | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--reference` |  | PhyloBuildParams.phylo.reference | none | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--aligner-arg` |  | PhyloBuildParams.phylo.extra | callback | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/developing.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | PhyloBuildParams.phylo.threads | range | task.cpus | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--sketches-dir` |  | PhyloBuildParams.sketches_dir | stage | module: --sketches-dir sketches when params.sketch stages sketches/ | n/a: unit tests test_sketch_consumers.py; Nextflow stub test with params.sketch | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | PhyloBuildParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--msa-only` |  | PhyloBuildParams.msa_only | stage | module: PHYLO_MSA | tests/live/test_steps.py::test_phylo_build_split_into_msa_and_tree | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa` |  | PhyloBuildParams.msa | stage | module: PHYLO_TREE | tests/live/test_steps.py::test_phylo_build_split_into_msa_and_tree | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--mask` |  | PhyloBuildParams.phylo.extra | registry | params.phylo_args | n/a: the mask path is exercised on phylo (test_phylo_mask_gubbins); the step shares the builder | docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--core-fraction` |  | PhyloBuildParams.phylo.extra | range | params.phylo_args | n/a: shares the phylo builder (wiring test); the typing pass is covered offline in tests/integration/test_snptype_stage.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--coregene` |  | PhyloBuildParams.phylo.coregene | registry | params.phylo_args | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--core-threshold` |  | PhyloBuildParams.phylo.extra | range | params.phylo_args | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--clean-mode` |  | PhyloBuildParams.phylo.extra | choice | params.phylo_args | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--annotations-dir` |  | PhyloBuildParams.annotations_dir | stage | module: the ANNOTATE records, staged and indexed by PHYLO and PHYLO_MSA | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--annotator` |  | PhyloBuildParams.annotator | registry | n/a: ANNOTATE uses params.annotator; phylo-build reads this flag, from params.phylo_args, only when no records are staged | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--bakta-db` |  | PhyloBuildParams.bakta_db | stage | n/a: ANNOTATE uses params.bakta_db; phylo-build reads this flag, from params.phylo_args, only when no records are staged | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/install.md, docs/audit/cli-matrix.md |

## annotate

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | AnnotateStageParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | README.md, docs/cli-reference.md, docs/usage.md, docs/install.md, docs/audit/cli-matrix.md |
| `--all-genomes` |  | AnnotateStageParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--bakta-db` |  | AnnotateStageParams.bakta_db | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | docs/cli-reference.md, docs/usage.md, docs/install.md, docs/audit/cli-matrix.md |
| `--threads` | -t | AnnotateStageParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | AnnotateStageParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | AnnotateStageParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_annotate_stage.py (the live run needs bakta and its database) | docs/cli-reference.md, docs/audit/cli-matrix.md |

## annotate-index

dispatch: `step:repgenr.stages.annotate.annotate_index`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--annotations-dir` |  | AnnotateIndexParams.annotations_dir | stage | module: fixed by the PHYLO and PHYLO_MSA scripts | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_index_step | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool` |  | AnnotateIndexParams.tool | registry | params.annotator | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_index_step | docs/cli-reference.md, docs/audit/cli-matrix.md |

## annotate-build

dispatch: `step:repgenr.stages.annotate.annotate_build`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genome` |  | AnnotateBuildParams.genome | stage | module: fixed by the ANNOTATE script | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | AnnotateBuildParams.out_dir | none | module: fixed by the ANNOTATE script | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--record-name` |  | AnnotateBuildParams.record_name | none | n/a: the ANNOTATE task names the record after the genome file | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool` |  | AnnotateBuildParams.tool | registry | params.annotator | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--bakta-db` |  | AnnotateBuildParams.bakta_db | stage | module: params.bakta_db, staged into the task | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/install.md, docs/audit/cli-matrix.md |
| `--threads` | -t | AnnotateBuildParams.threads | range | module: task.cpus | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | AnnotateBuildParams.extra | callback | params.annotate_args | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | AnnotateBuildParams.versions_out | none | module: fixed by the ANNOTATE script | n/a: covered offline in tests/unit/test_annotate_stage.py::test_annotate_build_step_writes_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |

## assemble

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--assembler` |  | AssembleParams.assembler | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | AssembleParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--jobs` |  | AssembleParams.jobs | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--memory-gb` |  | AssembleParams.memory_gb | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--min-contig-length` |  | AssembleParams.min_contig_length | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--polisher` |  | AssembleParams.polisher | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_reads.py::test_ont_assembly_is_polished_with_medaka | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--polish-rounds` |  | AssembleParams.polish_rounds | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_reads.py::test_ont_assembly_is_polished_with_medaka | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup` |  | AssembleParams.outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-reads` |  | AssembleParams.keep_reads | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | AssembleParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | AssembleParams.extra | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--checkm2-db` |  | AssembleParams.checkm2_db | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--min-completeness` |  | AssembleParams.min_completeness | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--max-contamination` |  | AssembleParams.max_contamination | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--classifier` |  | AssembleParams.classifier | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_stage.py | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-sketch` |  | AssembleParams.gtdb_sketch | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-lineages` |  | AssembleParams.gtdb_lineages | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--append` |  | AssembleParams.append | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests in tests/integration/test_assemble_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--max-runs` |  | AssembleParams.max_runs | range | n/a: workdir command; the Nextflow layer bounds the selection with params.reads_args --max-runs | n/a: offline stage tests with a registered fake assembler in tests/integration/test_assemble_max_runs.py | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--sketch` |  | AssembleParams.sketch | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: the shared sketch step is exercised live through ingest (tests/live/test_sketches_live.py) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--reads-sketch` |  | AssembleParams.reads_sketch | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: offline stage tests with a fake sourmash in tests/integration/test_assemble_reads_sketch.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--screen-reads` |  | AssembleParams.screen_reads | none | n/a: workdir command; the Nextflow layer passes the flag to assemble-run through params.assemble_args | n/a: offline stage tests with fake gather tables in tests/integration/test_assemble_screen.py; a live check is reported in the pull request | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--screen-min-fraction` |  | AssembleParams.screen_min_fraction | unit_interval | n/a: workdir command; the Nextflow layer passes the flag to assemble-run through params.assemble_args | n/a: offline stage tests in tests/integration/test_assemble_screen.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--screen-dup-ani` |  | AssembleParams.screen_dup_ani | unit_interval | n/a: the duplicate check needs every run and is made by the workdir command only | n/a: offline stage tests in tests/integration/test_assemble_screen.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |

## reads

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_reads_selects_runs_from_ena | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-family` | -tf | ReadsParams.target_family | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | ReadsParams.target_genus | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | ReadsParams.target_species | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_reads_selects_runs_from_ena | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--accession` |  | ReadsParams.accessions | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--accession-file` |  | ReadsParams.accession_file | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--platform` |  | ReadsParams.platform | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_reads_selects_runs_from_ena | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--max-runs` |  | ReadsParams.max_runs | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_reads_selects_runs_from_ena | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--min-bases` |  | ReadsParams.min_bases | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--max-bases` |  | ReadsParams.max_bases | range | params.reads_args | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--drop-selection` |  | ReadsParams.drop_selection | none | params.reads_args | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--one-per-sample` |  | ReadsParams.one_per_sample | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: network stage; covered offline on frozen ENA fixtures in tests/unit/test_reads_stage.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |

## run

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--viral` |  | VmetadataParams | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--genomes-dir` |  | n/a: selects the local chain (ingest replaces metadata and genome); wiring test | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_run_local_chain | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--from-workdir` |  | IngestParams.from_workdirs | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered offline in tests/integration/test_ingest_from_workdir.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--reads` |  | n/a: selects the reads chain (reads and assemble replace metadata and genome); wiring test | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--accession-file` |  | ReadsParams.accession_file | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--platform` |  | ReadsParams.platform | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--max-runs` |  | ReadsParams.max_runs | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--assembler` |  | AssembleParams.assembler | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--selection` |  | IngestParams.selection | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup` |  | IngestParams.outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--copy` |  | IngestParams.copy | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--with-snptype` |  | SnptypeParams | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--metadata-path` |  | MetadataParams.metadata_path | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--nodownload` |  | MetadataParams.nodownload | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GenomeParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--pre-primary-ani` |  | DereplicateParams.pre_primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--pre-secondary-ani` |  | DereplicateParams.pre_secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--node-basename` |  | Tree2taxParams.node_basename | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--naming` |  | Tree2taxParams.naming | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_tree2tax_naming.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--root-name` |  | Tree2taxParams.root_name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--remove-outgroup` |  | Tree2taxParams.remove_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--dataset` | -d | MetadataParams.dataset | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--level` | -l | MetadataParams.level | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-family` | -tf | MetadataParams.target_family | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | MetadataParams.target_genus | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | MetadataParams.target_species | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--release` | -r | MetadataParams.release | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to metadata unchanged (wiring test); the TSV path is exercised on metadata in test_network.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-version` |  | MetadataParams.version | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to metadata unchanged (wiring test); the TSV path is exercised on metadata in test_network.py | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--metadata-source` |  | MetadataParams.source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target` |  | VmetadataParams.target | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--viral-source` |  | VmetadataParams.source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--group-segments` |  | VgenomeParams.group_segments | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to vgenome unchanged (wiring test); covered on vgenome in test_vgenome_selection_flags | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool` |  | DereplicateParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--primary-ani` |  | DereplicateParams.primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` |  | DereplicateParams.secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` |  | DereplicateParams.aligned_fraction | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | DereplicateParams.keeper | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_species_set.py covers them on phylo) | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--aligner` |  | PhyloParams.aligner | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_species_set.py covers them on phylo) | docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--root` |  | PhyloParams.root | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_treebuilders_offline.py covers --root on phylo) | docs/usage.md, docs/output.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | DereplicateParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--dry-run` |  | n/a: prints the chain and exits before any stage | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_prints_the_chain_without_network | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--limit` |  | MetadataParams.limit | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--complete-only` |  | VmetadataParams.complete_only | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--host` |  | VmetadataParams.host | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--released-after` |  | VmetadataParams.released_after | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--process-size` | -s | DereplicateParams.process_size | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--num-processes` | -p | DereplicateParams.num_processes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--reduce` |  | DereplicateParams.reduce | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-reps` |  | DereplicateParams.target_reps | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | DereplicateParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--all-genomes` |  | PhyloParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--bootstrap` | -B | PhyloParams.bootstrap | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--reference` |  | PhyloParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--aligner-arg` |  | PhyloParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--mask` |  | PhyloParams.extra | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--core-fraction` |  | PhyloParams.extra | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: shares the phylo builder (wiring test); the typing pass is covered offline in tests/integration/test_snptype_stage.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-support` |  | Tree2taxParams.collapse_support | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--collapse-length` |  | Tree2taxParams.collapse_length | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--collapse-diameter` |  | Tree2taxParams.collapse_diameter | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_tree2tax_collapse.py | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--ranks` |  | Tree2taxParams.ranks | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--rank-red` |  | Tree2taxParams.rank_red | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | DereplicateParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to the stage builder unchanged (wiring test); the stage flag is covered live | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--coregene` |  | PhyloParams.coregene | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--core-threshold` |  | PhyloParams.extra | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--clean-mode` |  | PhyloParams.extra | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--annotator` |  | AnnotateStageParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--bakta-db` |  | AnnotateStageParams.bakta_db | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_phylo_coregene.py (the live run needs bakta, panaroo and the bakta database) | docs/cli-reference.md, docs/install.md, docs/audit/cli-matrix.md |

## snptype

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | SnptypeParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_parsnp_typer | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--reference` |  | SnptypeParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/audit/cli-matrix.md |
| `--all-genomes` |  | SnptypeParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_on_representatives_only | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--mask` |  | SnptypeParams.mask | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_gubbins_mask_changes_the_core_alignment | README.md, docs/developing.md, docs/choosing-tools.md, docs/cli-reference.md, docs/output.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--core-fraction` |  | SnptypeParams.core_fraction | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/integration/test_snptype_stage.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--threads` | -t | SnptypeParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | SnptypeParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_typer_and_tool_arg | docs/developing.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | SnptypeParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_snptype_allow_incomplete | docs/cli-reference.md, docs/audit/cli-matrix.md |

## status

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--json` |  | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_doctor_and_status_json_and_quick | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |

## tree2tax

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--node-basename` |  | Tree2taxParams.node_basename | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--naming` |  | Tree2taxParams.naming | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_tree2tax_naming.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--root-name` |  | Tree2taxParams.root_name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--remove-outgroup` |  | Tree2taxParams.remove_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--collapse-support` |  | Tree2taxParams.collapse_support | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-length` |  | Tree2taxParams.collapse_length | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-diameter` |  | Tree2taxParams.collapse_diameter | range | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/unit/test_tree2tax_collapse.py | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--ranks` |  | Tree2taxParams.ranks | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/integration/test_tree2tax.py and tests/unit/test_tree2tax_ranks.py | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--rank-red` |  | Tree2taxParams.rank_red | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: covered offline in tests/integration/test_tree2tax.py and tests/unit/test_tree2tax_ranks.py | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |

## tree2tax-relations

dispatch: `step:repgenr.stages.tree2tax.tree2tax_relations`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--tree` |  | Tree2taxStepParams.tree | stage | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | Tree2taxStepParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--clusters` |  | Tree2taxStepParams.clusters | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-dir` |  | Tree2taxStepParams.outgroup_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | Tree2taxStepParams.outgroup_accession | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--segments` |  | Tree2taxStepParams.segments | stage | module: staged from VACQUIRE on the viral path; an empty list elsewhere | n/a: covered offline in tests/integration/test_dataflow_phylo_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--node-basename` |  | Tree2taxStepParams.node_basename | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--naming` |  | Tree2taxStepParams.naming | choice | params.tree2tax_args | n/a: covered offline in tests/unit/test_tree2tax_naming.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--root-name` |  | Tree2taxStepParams.root_name | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--remove-outgroup` |  | Tree2taxStepParams.remove_outgroup | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | Tree2taxStepParams.no_outgroup | none | params.phylo_args | n/a: covered offline in tests/integration/test_dataflow_phylo_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-root` |  | Tree2taxStepParams.keep_root | none | params.phylo_args | n/a: covered offline in tests/integration/test_phylo_root.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--include-dereplicated` |  | Tree2taxStepParams.include_dereplicated | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | Tree2taxStepParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--collapse-support` |  | Tree2taxStepParams.collapse_support | range | params.tree2tax_args | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-length` |  | Tree2taxStepParams.collapse_length | range | params.tree2tax_args | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-diameter` |  | Tree2taxStepParams.collapse_diameter | range | params.tree2tax_args | n/a: covered offline in tests/unit/test_tree2tax_collapse.py | docs/cli-reference.md, docs/usage.md, docs/choosing-tools.md, docs/audit/cli-matrix.md |
| `--ranks` |  | Tree2taxStepParams.ranks | none | params.tree2tax_args | n/a: covered offline in tests/integration/test_tree2tax.py and tests/unit/test_tree2tax_ranks.py | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--rank-red` |  | Tree2taxStepParams.rank_red | stage | params.tree2tax_args | n/a: covered offline in tests/integration/test_tree2tax.py and tests/unit/test_tree2tax_ranks.py | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |

## versions

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | module: fixed by the process script | tests/live/test_aux_commands.py::test_status_and_versions | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | n/a: query command, no stage parameters | none | module: fixed by the process script | tests/live/test_aux_commands.py::test_status_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |

## vgenome

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | VgenomeParams.target_genus | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | VgenomeParams.target_species | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-serotype` | -tse | VgenomeParams.target_serotype | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-custom` | -tc | VgenomeParams.target_custom | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--length-all` |  | VgenomeParams.length_all | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--length-deviation` |  | VgenomeParams.length_deviation | range | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--length-method` |  | VgenomeParams.length_method | choice | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--length-range` |  | VgenomeParams.length_range | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--discard` |  | VgenomeParams.discard | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | VgenomeParams.no_outgroup | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | VgenomeParams.outgroup_accession | none | params.vgenome_args | n/a: needs a specific downloaded record; covered offline on synthetic records in tests/integration/test_ncbi_virus.py and test_bvbrc_select.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--group-segments` |  | VgenomeParams.group_segments | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-candidates-taxid-min-genomes` |  | VgenomeParams.outgroup_candidates_taxid_min_genomes | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-treebuilder` |  | VgenomeParams.outgroup_treebuilder | registry | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--glance` |  | VgenomeParams.glance | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--print-fasta-headers` |  | VgenomeParams.print_fasta_headers | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--ignore-duplicates` |  | VgenomeParams.ignore_duplicates | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_bvbrc_needs_ignore_duplicates | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | VgenomeParams.keep_files | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--sketch` |  | VgenomeParams.sketch | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: the shared sketch step is exercised live through ingest (tests/live/test_sketches_live.py) | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## vmetadata

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/cli-reference.md, docs/usage.md, docs/output.md, docs/audit/cli-matrix.md |
| `--target` |  | VmetadataParams.target | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--source` |  | VmetadataParams.source | choice | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--filter` |  | VmetadataParams.filter | stage | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--host` |  | VmetadataParams.host | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--complete-only` |  | VmetadataParams.complete_only | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--released-after` |  | VmetadataParams.released_after | callback | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--list` |  | VmetadataParams.list_targets | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_list_targets_reaches_bvbrc | docs/cli-reference.md, docs/audit/cli-matrix.md |

## assemble-run

dispatch: `step:repgenr.stages.assemble_steps.assemble_run`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--reads-tsv` |  | AssembleRunParams.reads_tsv | stage | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--run` |  | AssembleRunParams.run | stage | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | AssembleRunParams.out_dir | none | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--assembler` |  | AssembleRunParams.assembler | registry | params.assembler | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | AssembleRunParams.threads | range | module: task.cpus | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--memory-gb` |  | AssembleRunParams.memory_gb | range | module: task.memory | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--min-contig-length` |  | AssembleRunParams.min_contig_length | range | params.assemble_args | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--polisher` |  | AssembleRunParams.polisher | choice | params.assemble_args | n/a: offline step tests in tests/integration/test_assemble_steps.py | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--polish-rounds` |  | AssembleRunParams.polish_rounds | range | params.assemble_args | n/a: offline step tests in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-reads` |  | AssembleRunParams.keep_reads | none | n/a: not exposed by the READS_ASSEMBLE module | n/a: offline step tests in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | AssembleRunParams.keep_files | none | n/a: not exposed by the READS_ASSEMBLE module | n/a: offline step tests in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | AssembleRunParams.extra | none | params.assemble_args | n/a: offline step tests in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--reads-sketch` |  | AssembleRunParams.reads_sketch | none | params.assemble_args | n/a: offline step tests with a fake sourmash in tests/integration/test_assemble_reads_sketch.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--screen-reads` |  | AssembleRunParams.screen_reads | none | params.assemble_args | n/a: offline step tests with fake gather tables in tests/integration/test_assemble_screen.py | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--screen-min-fraction` |  | AssembleRunParams.screen_min_fraction | unit_interval | params.assemble_args | n/a: offline step tests in tests/integration/test_assemble_screen.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-sketch` |  | AssembleRunParams.gtdb_sketch | none | module: conf/modules.config adds params.gtdb_sketch when params.assemble_args holds --screen-reads | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_screen.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-lineages` |  | AssembleRunParams.gtdb_lineages | none | module: conf/modules.config adds params.gtdb_lineages when params.assemble_args holds --screen-reads | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_screen.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | AssembleRunParams.versions_out | none | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/cli-reference.md, docs/audit/cli-matrix.md |

## genome-qc

dispatch: `step:repgenr.stages.assemble_steps.genome_qc`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--assemblies` |  | GenomeQcParams.assemblies_dir | stage | module: fixed by the process script | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | GenomeQcParams.out_dir | none | module: fixed by the process script | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | GenomeQcParams.threads | range | module: task.cpus | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--memory-gb` |  | GenomeQcParams.memory_gb | range | module: task.memory | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--checkm2-db` |  | GenomeQcParams.checkm2_db | none | params.checkm2_db | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--classifier` |  | GenomeQcParams.classifier | choice | n/a: sourmash when params.gtdb_sketch is set | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/choosing-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-sketch` |  | GenomeQcParams.gtdb_sketch | none | params.gtdb_sketch | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-lineages` |  | GenomeQcParams.gtdb_lineages | none | params.gtdb_lineages | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | GenomeQcParams.extra | none | n/a: not exposed by the GENOME_QC module | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | GenomeQcParams.versions_out | none | module: fixed by the process script | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |

## reads-gather

dispatch: `step:repgenr.stages.assemble_steps.reads_gather`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--reads-tsv` |  | ReadsGatherParams.reads_tsv | stage | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--assemblies` |  | ReadsGatherParams.assemblies_dir | stage | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | ReadsGatherParams.out_dir | none | module: fixed by the process script | tests/live/test_reads.py::test_reads_steps_assemble_a_public_run | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--qc` |  | ReadsGatherParams.qc_dir | none | module: fixed by the process script | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/usage.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--min-completeness` |  | ReadsGatherParams.min_completeness | range | n/a: site config withName READS_GATHER ext.args | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--max-contamination` |  | ReadsGatherParams.max_contamination | range | n/a: site config withName READS_GATHER ext.args | n/a: needs a reference database; covered offline with fakes in tests/integration/test_assemble_steps.py | docs/cli-reference.md, docs/audit/cli-matrix.md |

## Short-alias collisions

