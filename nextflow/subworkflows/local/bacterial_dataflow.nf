// Bacterial pipeline, data-channel form:
//   ACQUIRE (metadata -> genome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//   (PHYLO is PHYLO_MSA + PHYLO_TREE when params.phylo_split_msa is set)
//
// With params.sketch, SKETCH writes a sourmash sketch of every genome and of
// the outgroup (published as sketches/), and the directory is staged into the
// dereplication and tree tasks; a sourmash dereplicator or tree builder reads
// it instead of sketching. SKETCH calls sourmash itself, so it is off by
// default and a run without sourmash is unaffected.
//
// With --msa-source coregene in phylo_args, ANNOTATE_SCATTER annotates every
// representative (and the outgroup) in its own task, and the record
// directories are staged into PHYLO or PHYLO_MSA.
//
// Genomes, representatives, the tree and the taxonomy flow between processes as
// staged channel files; there is no shared working directory. This is the
// data-channel replacement for the legacy done-signal BACTERIAL subworkflow.

include { ACQUIRE             } from './acquire'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { ANNOTATE_SCATTER    } from './annotate_scatter'
include { phylo_uses_coregene } from './annotate_scatter'
include { phylo_uses_markers  } from './annotate_scatter'
include { SKETCH              } from '../../modules/local/dataflow/sketch'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { PHYLO_MSA           } from '../../modules/local/dataflow/phylo_msa'
include { PHYLO_TREE          } from '../../modules/local/dataflow/phylo_tree'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow BACTERIAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    ACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(ACQUIRE.out.versions)

    def ch_sketches = ACQUIRE.out.genomes.map { meta, _files -> tuple(meta, []) }
    if (params.sketch) {
        SKETCH(
            ACQUIRE.out.genomes
                .join(ACQUIRE.out.outgroup, by: 0)
                .map { meta, files, outgroup -> tuple(meta, files + outgroup) }
        )
        ch_versions = ch_versions.mix(SKETCH.out.versions)
        ch_sketches = SKETCH.out.sketches
    }
    // One directory (or nothing) for the tree tasks, as for the dereplication.
    def ch_sketches_dir = ch_sketches.map { _meta, dir -> dir }.first()

    DEREPLICATE_SCATTER(ACQUIRE.out.genomes, ACQUIRE.out.selection, ch_sketches)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ACQUIRE.out.outgroup, by: 0)
        .join(ACQUIRE.out.outgroup_accession, by: 0)

    // The core-gene route's annotations, or an empty list in their place.
    def ch_build_in = ch_phylo_in.map { m, reps, og, acc -> tuple(m, reps, og, acc, []) }
    if (phylo_uses_coregene(params.phylo_args)) {
        ANNOTATE_SCATTER(DEREPLICATE_SCATTER.out.reps, ACQUIRE.out.outgroup)
        ch_versions = ch_versions.mix(ANNOTATE_SCATTER.out.versions)
        // failOnMismatch: a run whose annotations are missing stops here.
        ch_build_in = ch_phylo_in.join(ANNOTATE_SCATTER.out.annotations, by: 0, failOnMismatch: true)
    }

    // The marker database, staged into the alignment task only for the
    // markers builder: phylo-build refuses --markers-db on any other route,
    // so a site-wide params.markers_db must not reach those runs.
    def markers_db = []
    if (params.markers_db && phylo_uses_markers(params.phylo_args)) {
        markers_db = file(params.markers_db, checkIfExists: true)
    }

    // With params.phylo_split_msa the alignment and the tree are separate
    // tasks, so changing the tree builder re-runs only the tree.
    def ch_tree = null
    if (params.phylo_split_msa) {
        PHYLO_MSA(ch_build_in, ch_sketches_dir, markers_db)
        ch_versions = ch_versions.mix(PHYLO_MSA.out.versions)
        // alignment_stats.json exists for a SNP-typer alignment only; an empty
        // list stands in for it otherwise (same remainder pattern as the outgroup).
        def ch_stats = PHYLO_MSA.out.msa
            .map { m, _msa -> tuple(m, []) }
            .join(PHYLO_MSA.out.stats, by: 0, remainder: true)
            .map { m, _placeholder, st -> tuple(m, st == null ? [] : st) }
        // partitions.json exists for a core-gene alignment only (same pattern).
        def ch_partitions = PHYLO_MSA.out.msa
            .map { m, _msa -> tuple(m, []) }
            .join(PHYLO_MSA.out.partitions, by: 0, remainder: true)
            .map { m, _placeholder, pa -> tuple(m, pa == null ? [] : pa) }
        PHYLO_TREE(
            ch_phylo_in.join(PHYLO_MSA.out.msa, by: 0).join(ch_stats, by: 0).join(ch_partitions, by: 0),
            ch_sketches_dir
        )
        ch_versions = ch_versions.mix(PHYLO_TREE.out.versions)
        ch_tree = PHYLO_TREE.out.tree
    }
    else {
        PHYLO(ch_build_in, ch_sketches_dir, markers_db)
        ch_versions = ch_versions.mix(PHYLO.out.versions)
        ch_tree = PHYLO.out.tree
    }

    // [meta, tree] joined with [meta, reps, outgroup, accession]; no segments
    // table on the bacterial path
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
