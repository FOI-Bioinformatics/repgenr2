// Viral pipeline, data-channel form:
//   VACQUIRE (vmetadata -> vgenome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Reuses the same data-channel dereplication, phylo and tree2tax as the bacterial
// path; only the acquisition front differs (BV-BRC instead of GTDB/NCBI).

include { VACQUIRE            } from '../../modules/local/dataflow/vacquire'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { PHYLO_MSA           } from '../../modules/local/dataflow/phylo_msa'
include { PHYLO_TREE          } from '../../modules/local/dataflow/phylo_tree'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow VIRAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    VACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(VACQUIRE.out.versions)

    def ch_genomes = VACQUIRE.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }
    // vmetadata writes no selection.tsv and viral genomes carry no CheckM
    // quality, so the keeper input is an empty list under the run meta.
    def ch_no_selection = ch_meta.map { meta -> tuple(meta, []) }
    // Genome sketches are not staged on the viral path (params.sketch is
    // bacterial only); the sourmash tools sketch for themselves.
    DEREPLICATE_SCATTER(ch_genomes, ch_no_selection, ch_no_selection)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    def ch_outgroup = ch_genomes
        // A placeholder tuple, not a bare meta: with remainder, an unmatched
        // bare value is emitted as-is and the next closure gets one argument
        // (seen live on a viral run with --no-outgroup).
        .map { meta, _files -> tuple(meta, []) }
        .join(VACQUIRE.out.outgroup, by: 0, remainder: true)
        .map { meta, _placeholder, files ->
            def list = files == null ? [] : (files instanceof List ? files : [files])
            tuple(meta, list)
        }
    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ch_outgroup, by: 0)
        .join(VACQUIRE.out.outgroup_accession, by: 0)
    // No annotation step on the viral path (the annotators are prokaryotic);
    // an empty list stands in for the annotations the phylo tasks can stage.
    def ch_build_in = ch_phylo_in.map { m, reps, og, acc -> tuple(m, reps, og, acc, []) }

    // See bacterial_dataflow: params.phylo_split_msa runs the alignment and
    // the tree as separate tasks.
    def ch_tree = null
    if (params.phylo_split_msa) {
        PHYLO_MSA(ch_build_in, [])
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
        PHYLO(ch_build_in, [])
        ch_versions = ch_versions.mix(PHYLO.out.versions)
        ch_tree = PHYLO.out.tree
    }

    // segments.tsv exists only after vgenome --group-segments; an empty list
    // stands in for it otherwise (same remainder pattern as the outgroup).
    def ch_segments = ch_genomes
        .map { meta, _files -> tuple(meta, []) }
        .join(VACQUIRE.out.segments, by: 0, remainder: true)
        .map { meta, _placeholder, seg -> tuple(meta, seg == null ? [] : seg) }
    // [meta, tree] joined with [meta, reps, outgroup, accession] and [meta, segments]
    def ch_tree2tax_in = ch_tree.join(ch_phylo_in, by: 0).join(ch_segments, by: 0)

    TREE2TAX(ch_tree2tax_in)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)

    emit:
    tree        = ch_tree
    tree2tax    = TREE2TAX.out.tree2tax
    genomes_map = TREE2TAX.out.genomes_map
    versions    = ch_versions
}
