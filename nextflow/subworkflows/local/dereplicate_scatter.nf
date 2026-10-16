// Scatter-gather dereplication (data-channel style).
//
// Genomes arrive as one tuple(meta, [fasta...]) per run. They are grouped into
// chunks of `derep_process_size`, each chunk is dereplicated independently and
// in parallel (scatter, one task per chunk -- on HPC these land on separate
// nodes), and the union of the chunk representatives is dereplicated once more
// (gather) to produce the final representative set. This is the two-stage
// reduce-tree of the in-process dereplicate stage, expressed as typed channels
// so Nextflow owns the fan-out instead of Python threads.
//
// Chunk and merge metas nest the run meta under `run`, so the result can be
// emitted under the run meta again and joined with the other run-level channels.
//
// With params.derep_refine or params.derep_refine_file, a DEREP_REFINE task
// follows the merge: the clusters of the seed genomes are dereplicated again
// at finer ANIs, and its result is the one emitted.

include { DEREP_CHUNK  } from '../../modules/local/derep_chunk'
include { DEREP_MERGE  } from '../../modules/local/derep_merge'
include { DEREP_REFINE } from '../../modules/local/derep_refine'

workflow DEREPLICATE_SCATTER {
    take:
    ch_genomes    // channel: tuple(meta, [genome FASTA paths])
    ch_selection  // channel: tuple(meta, selection.tsv) or tuple(meta, [])
    ch_sketches   // channel: tuple(meta, sketches/ directory) or tuple(meta, [])

    main:
    def ch_versions = channel.empty()

    // A null/zero process size means a single chunk. Coerce because
    // command-line params arrive as strings.
    def requested_size = params.derep_process_size ? (params.derep_process_size as Integer) : 0
    def chunk_size = requested_size > 0 ? requested_size : 1000000

    def ch_chunks = ch_genomes.flatMap { meta, files ->
        files.collate(chunk_size).withIndex().collect { chunk, i ->
            tuple([id: "${meta.id}.chunk_${i}", mode: meta.mode, run: meta], chunk)
        }
    }

    // The selection file is an auxiliary input shared by every chunk task, so
    // it is passed as a bare value channel (nf-core reference-file style).
    def ch_selection_file = ch_selection.map { _meta, sel -> sel }.first()
    // The genome sketches (SKETCH, params.sketch) likewise: one directory
    // staged into every chunk and the merge; a sourmash dereplicator reads
    // it, the other tools ignore it.
    def ch_sketches_dir = ch_sketches.map { _meta, dir -> dir }.first()

    DEREP_CHUNK(ch_chunks, ch_selection_file, ch_sketches_dir)
    ch_versions = ch_versions.mix(DEREP_CHUNK.out.versions.first())

    // Gather every chunk directory under one merge meta per run.
    def ch_merge_in = DEREP_CHUNK.out.chunk
        .map { meta, dir -> tuple([id: "${meta.run.id}.merged", mode: meta.mode, run: meta.run], dir) }
        .groupTuple()

    DEREP_MERGE(ch_merge_in, ch_selection_file, ch_sketches_dir)
    ch_versions = ch_versions.mix(DEREP_MERGE.out.versions)

    def ch_reps = DEREP_MERGE.out.reps.map { meta, dir -> tuple(meta.run, dir) }
    if (params.derep_refine || params.derep_refine_file) {
        // Refused at launch, not after the chunk and merge tasks.
        if (params.derep_refine && !params.derep_refine_ani) {
            error "--derep_refine needs --derep_refine_ani"
        }
        if (!params.derep_refine && params.derep_refine_ani) {
            error "--derep_refine_ani applies to --derep_refine seeds; a refine file carries its own ani column"
        }
        if (params.derep_refine_ani && (params.derep_refine_ani as Double) <= (params.derep_secondary_ani as Double)) {
            error "--derep_refine_ani ${params.derep_refine_ani} is not finer than --derep_secondary_ani ${params.derep_secondary_ani}"
        }
        // The merged result joined with every genome of its run: the lineage
        // members are not among the merge's representatives.
        def ch_refine_in = ch_reps
            .join(ch_genomes, by: 0)
            .map { meta, dir, files -> tuple([id: "${meta.id}.refined", mode: meta.mode, run: meta], dir, files) }
        def ch_refine_file = params.derep_refine_file
            ? channel.value(file(params.derep_refine_file, checkIfExists: true))
            : channel.value([])
        DEREP_REFINE(ch_refine_in, ch_selection_file, ch_sketches_dir, ch_refine_file)
        ch_versions = ch_versions.mix(DEREP_REFINE.out.versions)
        ch_reps = DEREP_REFINE.out.reps.map { meta, dir -> tuple(meta.run, dir) }
    }

    emit:
    reps     = ch_reps
    versions = ch_versions
}
