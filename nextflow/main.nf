#!/usr/bin/env nextflow

// RepGenR2 pipeline entry point (data-channel form).
//
// Stages exchange typed channel files -- metadata selection (or a reads
// selection assembled per run), genome FASTAs,
// per-chunk and merged representatives, the tree and the taxonomy -- rather than
// sharing one working directory. Nextflow owns the fan-out (scatter-gather
// dereplication), and results are published under --outdir.

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { BACTERIAL_DATAFLOW } from './subworkflows/local/bacterial_dataflow'
include { VIRAL_DATAFLOW     } from './subworkflows/local/viral_dataflow'
include { READS_DATAFLOW     } from './subworkflows/local/reads_dataflow'
include { PUBLISH_VERSIONS   } from './subworkflows/local/publish_versions'
include { run_meta           } from './subworkflows/local/run_meta'

workflow {
    // Print usage and exit before validation, so `--help` needs no other params.
    if (params.help) {
        log.info """
        RepGenR2 -- ${workflow.manifest.description} (v${workflow.manifest.version})

        Usage:
          nextflow run main.nf --outdir <DIR> [-profile standard|slurm|docker|singularity|wave|test]

        Key parameters (see nextflow_schema.json and docs/usage.md for all):
          --outdir <DIR>         Published results (default: results)
          --mode bacterial|viral|reads  Pipeline to run (default: bacterial)
          --metadata_args '<str>'  GTDB selection (bacterial)
          --vmetadata_args / --vgenome_args '<str>'  NCBI Virus selection (viral)
          --reads_args '<str>'   ENA/SRA run selection (reads); --assembler, --checkm2_db, --gtdb_sketch
          --reads_sketch, --reads_screen  Sketch the runs from the archive stream; screen them against GTDB before any download
          --derep_tool <tool>    Dereplicator for the scatter-gather step
          --derep_process_size N Genomes per dereplication chunk (large inputs)
          --phylo_args '<str>'   Aligner or tree builder for the phylogeny
          --annotator, --bakta_db  Annotation for phylo_args with --msa-source coregene
        """.stripIndent()
        return
    }

    validateParameters()
    // phylo_split_msa runs the alignment and the tree as separate tasks, which
    // needs a tree builder that consumes an alignment; the alignment-free
    // builders would fail inside the first task.
    if (params.phylo_split_msa && (params.phylo_args =~ /--treebuilder(\s+|=)(mashtree|sourmash)\b/).find()) {
        error "phylo_split_msa needs a tree builder that consumes an alignment (iqtree, raxmlng, fasttree); phylo_args is '${params.phylo_args}'"
    }
    if (params.mode != 'reads' && (params.reads_sketch || params.reads_screen)) {
        log.warn "reads_sketch and reads_screen apply to --mode reads only; ignored in ${params.mode} mode"
    }
    log.info paramsSummaryLog(workflow)

    def ch_meta = channel.value(run_meta(params))

    def ch_versions = channel.empty()
    if (params.mode == 'viral') {
        VIRAL_DATAFLOW(ch_meta)
        ch_versions = VIRAL_DATAFLOW.out.versions
    }
    else if (params.mode == 'reads') {
        READS_DATAFLOW(ch_meta)
        ch_versions = READS_DATAFLOW.out.versions
    }
    else {
        BACTERIAL_DATAFLOW(ch_meta)
        ch_versions = BACTERIAL_DATAFLOW.out.versions
    }

    PUBLISH_VERSIONS(ch_versions)
}
