#!/usr/bin/env nextflow

// Network-free data-channel pipeline over a local directory of genome FASTAs:
//   DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Used by the real (non-stub) end-to-end CI test with lightweight tools
// (sourmash + mashtree), so it skips the GTDB/NCBI ACQUIRE front and runs with
// no outgroup. Publishes tree2tax.tsv, genomes_map.tsv and the tree to --outdir.

include { DEREPLICATE_SCATTER } from '../subworkflows/local/dereplicate_scatter'
include { PHYLO               } from '../modules/local/dataflow/phylo'
include { TREE2TAX            } from '../modules/local/dataflow/tree2tax'
include { PUBLISH_VERSIONS    } from '../subworkflows/local/publish_versions'

params.genomes_dir = null
params.empty_accession = "${projectDir}/tests/data/empty.txt"

workflow {
    if (!params.genomes_dir) {
        error "Provide --genomes_dir <DIR> with genome FASTAs."
    }
    def meta = [id: 'local', mode: 'bacterial']
    def ch_genomes = channel
        .fromPath("${params.genomes_dir}/*.{fasta,fa,fna,fas,fasta.gz,fna.gz,fa.gz}")
        .filter { f -> !f.name.startsWith('._') }
        .collect()
        .map { files -> tuple(meta, files) }

    DEREPLICATE_SCATTER(ch_genomes, channel.value(tuple(meta, [])), channel.value(tuple(meta, [])))

    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .map { m, dir -> tuple(m, dir, [], file(params.empty_accession)) }
    PHYLO(ch_phylo_in, [])
    def ch_tree2tax_in = PHYLO.out.tree
        .join(ch_phylo_in, by: 0)
        .map { m, t, reps_dir, og, acc -> tuple(m, t, reps_dir, og, acc, []) }
    TREE2TAX(ch_tree2tax_in)

    def ch_versions = DEREPLICATE_SCATTER.out.versions
        .mix(PHYLO.out.versions, TREE2TAX.out.versions)
    PUBLISH_VERSIONS(ch_versions)
}
