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

// True when the phylo arguments select the core-gene alignment source, by
// --msa-source or, without it, by a core-gene --preset (phylo-build refuses
// a --msa-source other than the preset's). Read from params.phylo_args only
// (a withName ext.args override is not seen). Quotes around a value are
// dropped, and the last --msa-source or --preset wins, as on the repgenr
// command line.
def phylo_uses_coregene(args) {
    def tokens = (args ?: '').toString().tokenize().collect { t -> t.replaceAll(/^["']|["']$/, '') }
    def source = null
    def preset = null
    tokens.eachWithIndex { t, i ->
        if (t.startsWith('--msa-source=')) {
            source = t.substring('--msa-source='.length()).replaceAll(/^["']|["']$/, '')
        }
        else if (t == '--msa-source' && i + 1 < tokens.size()) {
            source = tokens[i + 1]
        }
        else if (t.startsWith('--preset=')) {
            preset = t.substring('--preset='.length()).replaceAll(/^["']|["']$/, '')
        }
        else if (t == '--preset' && i + 1 < tokens.size()) {
            preset = tokens[i + 1]
        }
    }
    if (source != null) {
        return source == 'coregene'
    }
    // The core-gene presets: the genus and family rows of RECOMMENDATIONS in
    // src/repgenr/stages/advise.py (tests/unit/test_cli_presets.py checks
    // that the two agree).
    return preset in ['genus', 'family']
}

// True when the phylo arguments select the markers core-gene builder: the
// core-gene route with the last --coregene being markers or, without an
// explicit --coregene, a preset whose builder is markers. Only then is
// params.markers_db staged into the phylo task; phylo-build refuses
// --markers-db on any other route.
def phylo_uses_markers(args) {
    if (!phylo_uses_coregene(args)) {
        return false
    }
    def tokens = (args ?: '').toString().tokenize().collect { t -> t.replaceAll(/^["']|["']$/, '') }
    def builder = null
    def preset = null
    tokens.eachWithIndex { t, i ->
        if (t.startsWith('--coregene=')) {
            builder = t.substring('--coregene='.length()).replaceAll(/^["']|["']$/, '')
        }
        else if (t == '--coregene' && i + 1 < tokens.size()) {
            builder = tokens[i + 1]
        }
        else if (t.startsWith('--preset=')) {
            preset = t.substring('--preset='.length()).replaceAll(/^["']|["']$/, '')
        }
        else if (t == '--preset' && i + 1 < tokens.size()) {
            preset = tokens[i + 1]
        }
    }
    if (builder != null) {
        return builder == 'markers'
    }
    // The presets whose --coregene is markers in RECOMMENDATIONS
    // (src/repgenr/stages/advise.py; tests/unit/test_cli_presets.py checks
    // that the two agree).
    return preset in markers_presets()
}

def markers_presets() {
    return ['family']
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
    // Listed through the Path API, so a work directory on object storage
    // (s3://, gs://, az://) keeps its scheme.
    def ch_genomes = ch_reps.flatMap { meta, reps_dir ->
        reps_dir.resolve('representatives').listFiles()
            .findAll { f -> f.isFile() && !f.name.startsWith('.') }
            .sort { f -> f.name }
            .collect { genome -> tuple(meta, genome) }
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
