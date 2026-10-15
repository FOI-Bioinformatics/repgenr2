// Annotation for the core-gene route (phylo_args with --msa-source coregene).
//
// ANNOTATE runs once per representative genome, and once for the outgroup
// unless phylo_args builds the tree without it (--no-outgroup). The record
// directories are gathered per run into one list, which the phylo task
// stages as annotations/ and indexes with `repgenr annotate-index` before
// `repgenr phylo-build --annotations-dir`. Annotating in one task per genome
// lets the annotations run in parallel and cache genome by genome; without
// this the phylo task would annotate every genome itself, one after another.

include { ANNOTATE } from '../../modules/local/dataflow/annotate'

// True when the phylo arguments select the core-gene alignment source.
def phylo_uses_coregene(args) {
    def tokens = (args ?: '').toString().tokenize()
    if (tokens.contains('--msa-source=coregene')) {
        return true
    }
    def at = tokens.indexOf('--msa-source')
    return at >= 0 && at + 1 < tokens.size() && tokens[at + 1] == 'coregene'
}

// True when the phylo arguments keep the outgroup in the tree.
def phylo_keeps_outgroup(args) {
    return !(args ?: '').toString().tokenize().contains('--no-outgroup')
}

workflow ANNOTATE_SCATTER {
    take:
    ch_reps      // [meta, reps_dir] with reps_dir/representatives/*
    ch_outgroup  // [meta, outgroup file or list of files, possibly empty]

    main:
    def ch_genomes = ch_reps.flatMap { meta, reps_dir ->
        files("${reps_dir}/representatives/*").collect { genome -> tuple(meta, genome) }
    }
    if (phylo_keeps_outgroup(params.phylo_args)) {
        ch_genomes = ch_genomes.mix(
            ch_outgroup.flatMap { meta, outgroup ->
                def list = outgroup == null ? [] : (outgroup instanceof List ? outgroup : [outgroup])
                list.collect { genome -> tuple(meta, genome) }
            }
        )
    }
    def bakta_db = params.bakta_db ? file(params.bakta_db, checkIfExists: true) : []

    ANNOTATE(ch_genomes, bakta_db)

    emit:
    // [meta, [record dirs]] once every genome of the run is annotated
    annotations = ANNOTATE.out.annotation.groupTuple(by: 0)
    versions    = ANNOTATE.out.versions
}
