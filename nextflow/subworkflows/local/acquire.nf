// Acquire genomes (data-channel form): metadata -> genome.
//
// METADATA selects accessions and emits selection.tsv; GENOME downloads them and
// emits the genome FASTAs. Every output carries the run meta, and the genome and
// outgroup lists are normalised to lists so DEREPLICATE_SCATTER can chunk them.

include { METADATA } from '../../modules/local/dataflow/metadata'
include { GENOME   } from '../../modules/local/dataflow/genome'
include { optional_files } from './report_inputs'

workflow ACQUIRE {
    take:
    ch_meta   // value channel: the run meta map

    main:
    def ch_versions = channel.empty()

    METADATA(ch_meta)
    ch_versions = ch_versions.mix(METADATA.out.versions)

    GENOME(METADATA.out.selection)
    ch_versions = ch_versions.mix(GENOME.out.versions)

    def ch_genomes = GENOME.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }

    // The outgroup output is optional: join with remainder so a run without one
    // still emits tuple(meta, []) and downstream joins on the meta keep working.
    def ch_outgroup = ch_genomes
        // A placeholder tuple, not a bare meta: with remainder, an unmatched
        // bare value is emitted as-is and the next closure gets one argument
        // (seen live on a viral run with --no-outgroup).
        .map { meta, _files -> tuple(meta, []) }
        .join(GENOME.out.outgroup, by: 0, remainder: true)
        .map { meta, _placeholder, files ->
            def list = files == null ? [] : (files instanceof List ? files : [files])
            tuple(meta, list)
        }

    emit:
    genomes            = ch_genomes
    outgroup           = ch_outgroup
    selection          = METADATA.out.selection
    outgroup_accession = METADATA.out.outgroup_accession
    // The selection summary and lineage table, for the run report ([] when absent).
    summary            = optional_files(METADATA.out.selection, METADATA.out.report)
    versions           = ch_versions
}
