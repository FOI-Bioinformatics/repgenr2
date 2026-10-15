// Reads pipeline, data-channel form:
//   ACQUIRE_READS (reads -> assemble per run -> qc -> gather)
//     -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Reuses the same data-channel dereplication, phylo and tree2tax as the
// bacterial path; only the acquisition front differs (ENA/SRA sequencing runs
// assembled here instead of downloaded assemblies).

include { ACQUIRE_READS       } from './acquire_reads'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { PHYLO_MSA           } from '../../modules/local/dataflow/phylo_msa'
include { PHYLO_TREE          } from '../../modules/local/dataflow/phylo_tree'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow READS_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    ACQUIRE_READS(ch_meta)
    ch_versions = ch_versions.mix(ACQUIRE_READS.out.versions)

    // Genome sketches are not staged on the reads path (params.sketch is
    // bacterial only); the sourmash tools sketch for themselves.
    DEREPLICATE_SCATTER(
        ACQUIRE_READS.out.genomes,
        ACQUIRE_READS.out.selection,
        ACQUIRE_READS.out.genomes.map { meta, _files -> tuple(meta, []) }
    )
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ACQUIRE_READS.out.outgroup, by: 0)
        .join(ACQUIRE_READS.out.outgroup_accession, by: 0)

    def ch_tree = null
    if (params.phylo_split_msa) {
        PHYLO_MSA(ch_phylo_in)
        ch_versions = ch_versions.mix(PHYLO_MSA.out.versions)
        // alignment_stats.json exists for a SNP-typer alignment only; an empty
        // list stands in for it otherwise (same remainder pattern as the outgroup).
        def ch_stats = PHYLO_MSA.out.msa
            .map { m, _msa -> tuple(m, []) }
            .join(PHYLO_MSA.out.stats, by: 0, remainder: true)
            .map { m, _placeholder, st -> tuple(m, st == null ? [] : st) }
        PHYLO_TREE(ch_phylo_in.join(PHYLO_MSA.out.msa, by: 0).join(ch_stats, by: 0), [])
        ch_versions = ch_versions.mix(PHYLO_TREE.out.versions)
        ch_tree = PHYLO_TREE.out.tree
    }
    else {
        PHYLO(ch_phylo_in, [])
        ch_versions = ch_versions.mix(PHYLO.out.versions)
        ch_tree = PHYLO.out.tree
    }

    // No segments table on the reads path: an empty list stands in for it.
    def ch_tree2tax_in = ch_tree
        .join(ch_phylo_in, by: 0)
        .map { m, t, reps_dir, og, acc -> tuple(m, t, reps_dir, og, acc, []) }
    TREE2TAX(ch_tree2tax_in)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)

    emit:
    tree        = ch_tree
    tree2tax    = TREE2TAX.out.tree2tax
    genomes_map = TREE2TAX.out.genomes_map
    versions    = ch_versions
}
