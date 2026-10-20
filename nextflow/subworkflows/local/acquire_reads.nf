// Acquire genomes from sequencing runs (data-channel form):
//   READS_SELECT -> [READS_SKETCH -> READS_SCREEN] -> READS_ASSEMBLE (one task per run)
//     -> [GENOME_QC] -> READS_GATHER
//
// With params.reads_sketch the runs are sketched from the archive stream
// (no download), one READS_SKETCH task per chunk of params.reads_sketch_chunk
// runs so a retry streams one chunk again, joined by READS_SKETCH_MERGE and
// published; with params.reads_screen the sketches are judged against the
// GTDB sketch (taxon, match fraction, duplicates) and only the runs of
// passed.txt become assembly tasks.
// READS_SELECT emits reads.tsv; each of its rows becomes one READS_ASSEMBLE task
// (the meta nests the run meta under `run` and carries the accession and the
// platform, which sets the task's memory). The per-run directories are grouped
// under the run meta again, scored and classified in one GENOME_QC task when a
// CheckM2 database or a GTDB sketch is configured, and gathered into the genome
// contract. The emitted tuples match ACQUIRE's, so the dereplication scatter,
// phylo and tree2tax steps run unchanged; the reads chain has no outgroup, so
// the outgroup channel carries an empty list and the accession file is empty.

include { READS_SELECT   } from '../../modules/local/dataflow/reads_select'
include { READS_SKETCH   } from '../../modules/local/dataflow/reads_sketch'
include { READS_SKETCH_MERGE } from '../../modules/local/dataflow/reads_sketch_merge'
include { READS_SCREEN   } from '../../modules/local/dataflow/reads_screen'
include { READS_ASSEMBLE } from '../../modules/local/dataflow/reads_assemble'
include { GENOME_QC      } from '../../modules/local/dataflow/genome_qc'
include { READS_GATHER   } from '../../modules/local/dataflow/reads_gather'

workflow ACQUIRE_READS {
    take:
    ch_meta   // value channel: the run meta map

    main:
    def ch_versions = channel.empty()

    READS_SELECT(ch_meta)
    ch_versions = ch_versions.mix(READS_SELECT.out.versions)

    // The runs to assemble: every run of reads.tsv (null), or with the screen
    // the runs of passed.txt. A screen that passed no run stops the pipeline
    // rather than assembling nothing, or everything, in silence.
    def ch_keep = READS_SELECT.out.reads.map { meta, _tsv -> tuple(meta, null) }
    if (params.reads_screen && !(params.gtdb_sketch && params.gtdb_lineages)) {
        error "--reads_screen needs --gtdb_sketch and --gtdb_lineages (the GTDB sketch and its lineages CSV)"
    }
    def ch_sketches = null
    if (params.reads_sketch || params.reads_screen) {
        // One task per chunk of runs, keyed by a nested meta; a batch
        // without a run is one empty chunk, so the merge still emits.
        def chunk_size = Math.max(1, (params.reads_sketch_chunk ?: 25) as int)
        def ch_chunks = READS_SELECT.out.reads.flatMap { meta, tsv ->
            def runs = tsv.splitCsv(header: true, sep: '\t').collect { row -> row.run_accession }.unique()
            // The accessions go into a shell line of the task unquoted.
            def odd = runs.findAll { run -> !(run ==~ /[SED]RR[0-9]+/) }
            if (odd) {
                error "reads.tsv of ${meta.id} holds a value that is not a run accession: ${odd.join(', ')}"
            }
            def chunks = runs ? runs.collate(chunk_size) : [[]]
            (0..<chunks.size()).collect { i ->
                tuple(
                    [id: "${meta.id}.${i + 1}", mode: meta.mode, run: meta, chunk: i + 1, chunks: chunks.size()],
                    tsv,
                    chunks[i]
                )
            }
        }
        READS_SKETCH(ch_chunks)
        ch_versions = ch_versions.mix(READS_SKETCH.out.versions)
        // The merge of a batch starts when its own chunks are done (the group
        // key carries the count), with the directories in chunk order.
        def ch_chunk_dirs = READS_SKETCH.out.sketches
            .map { meta, dir -> tuple(groupKey(meta.run, meta.chunks), [meta.chunk, dir]) }
            .groupTuple()
            .map { key, pairs ->
                tuple(key.getGroupTarget(), pairs.sort { pair -> pair[0] }.collect { pair -> pair[1] })
            }
        READS_SKETCH_MERGE(READS_SELECT.out.reads.join(ch_chunk_dirs, by: 0))
        ch_sketches = READS_SKETCH_MERGE.out.sketches
    }
    if (params.reads_screen) {
        READS_SCREEN(READS_SELECT.out.reads.join(ch_sketches, by: 0))
        ch_versions = ch_versions.mix(READS_SCREEN.out.versions)
        ch_keep = READS_SCREEN.out.passed.map { meta, passed ->
            def runs = passed.readLines().collect { line -> line.trim() }.findAll { line -> line }
            if (runs.isEmpty()) {
                error "READS_SCREEN (${meta.id}) passed no run: every run of reads.tsv was excused " +
                    "or could not be screened; see reads/reads_screen/<variant>/reads_screen.tsv"
            }
            tuple(meta, runs)
        }
    }

    // One item per run of reads.tsv, keyed by a nested meta.
    def ch_runs = READS_SELECT.out.reads.join(ch_keep, by: 0).flatMap { meta, tsv, keep ->
        tsv.splitCsv(header: true, sep: '\t')
            .findAll { row -> keep == null || keep.contains(row.run_accession) }
            .collect { row ->
                tuple(
                    [id: "${meta.id}.${row.run_accession}", mode: meta.mode, run: meta,
                     accession: row.run_accession, platform: row.platform],
                    tsv
                )
            }
    }

    READS_ASSEMBLE(ch_runs)
    // Every task's fragment: a mixed batch records both assemblers (the
    // publish step de-duplicates identical fragments).
    ch_versions = ch_versions.mix(READS_ASSEMBLE.out.versions)

    // Every run directory back under the run meta.
    def ch_assemblies = READS_ASSEMBLE.out.assembly
        .map { meta, dir -> tuple(meta.run, dir) }
        .groupTuple()

    def ch_qc = null
    if (params.checkm2_db || params.gtdb_sketch) {
        GENOME_QC(ch_assemblies)
        ch_versions = ch_versions.mix(GENOME_QC.out.versions)
        ch_qc = GENOME_QC.out.qc
    }
    else {
        ch_qc = ch_assemblies.map { meta, _dirs -> tuple(meta, []) }
    }

    READS_GATHER(READS_SELECT.out.reads.join(ch_assemblies, by: 0).join(ch_qc, by: 0))
    ch_versions = ch_versions.mix(READS_GATHER.out.versions)

    def ch_genomes = READS_GATHER.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }

    emit:
    genomes            = ch_genomes
    outgroup           = ch_genomes.map { meta, _files -> tuple(meta, []) }
    selection          = READS_GATHER.out.selection
    outgroup_accession = READS_GATHER.out.outgroup_accession
    reads              = READS_SELECT.out.reads
    versions           = ch_versions
}
