// Build the alignment a tree builder consumes (data-channel form).
//
// The first half of `repgenr phylo-build`, run as its own task so that the
// alignment is not repeated when only the tree builder or the bootstrap
// changes: those belong to PHYLO_TREE, which caches separately. Enabled by
// params.phylo_split_msa; otherwise PHYLO does both halves in one task.
// A staged sketches/ directory is passed as --sketches-dir, as in PHYLO:
// the mapping SNP typers and the aligners pick their reference genome (the
// medoid) from the sketches. Staged annotations/ (the ANNOTATE records,
// --msa-source coregene) are indexed and passed as --annotations-dir. A
// staged marker database (params.markers_db, --coregene markers) is passed
// as --markers-db.

process PHYLO_MSA {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession), path(annotations, stageAs: 'annotations/*')
    path sketches, stageAs: 'sketches'
    path markers_db, stageAs: 'markers_db'

    output:
    tuple val(meta), path("msa.fasta")                   , emit: msa
    // Constant-site counts of a SNP-typer alignment, for the tree task's
    // ascertainment-bias correction; absent for an aligner alignment.
    tuple val(meta), path("alignment_stats.json"), optional: true, emit: stats
    // Column ranges of the core genes in a core-gene alignment, for the tree
    // task's --partitioned; absent for the other sources.
    tuple val(meta), path("partitions.json"), optional: true, emit: partitions
    // The source directory the alignment was built in: the aligner's files
    // (align/) or the SNP typing pass's (tree/msa/), each with the reuse
    // stamp (msa_source.json).
    tuple val(meta), path("align/*")   , optional: true , emit: align
    tuple val(meta), path("tree/msa/*"), optional: true , emit: tree_msa
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    def annotator = task.ext.annotator ?: 'bakta'
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    sk=""
    [ -d sketches ] && sk="--sketches-dir sketches"
    # The marker database of --coregene markers (params.markers_db); without
    # one, phylo-build reads REPGENR_MARKERS_DB from the task environment.
    mk=""
    [ -d markers_db ] && mk="--markers-db markers_db"

    # The ANNOTATE record directories (--msa-source coregene), gathered here:
    # index them so phylo-build reads them instead of annotating every genome.
    an=""
    if [ -d annotations ]; then
        repgenr ${opts} annotate-index --annotations-dir annotations --tool ${annotator}
        an="--annotations-dir annotations"
    fi

    repgenr ${opts} phylo-build \\
        --genomes-dir ${reps_dir}/representatives \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        -o . -t ${task.cpus} ${args} \$sk \$an \$mk \\
        --msa-only \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    echo "ext.args: ${args}"
    echo "ext.repgenr_opts: ${opts}"
    if [ -d sketches ]; then echo "sketches: staged"; else echo "sketches: none"; fi
    if [ -d markers_db ]; then echo "markers_db: staged"; else echo "markers_db: none"; fi
    if [ -d annotations ]; then
        echo "annotations: staged (\$(ls annotations | wc -l | tr -d ' ') records)"
    else
        echo "annotations: none"
    fi
    for genome in ${reps_dir}/representatives/*; do
        name=\$(basename "\$genome" | sed 's/\\.[^.]*\$//')
        echo ">\${name}"
        echo "ACGT"
    done > msa.fasta
    # A core-gene alignment (annotations staged) carries its partitions.
    if [ -d annotations ]; then
        echo '{"v": 1, "alignment": "msa.fasta", "alphabet": "dna", "partitions": [["g1", 1, 4]]}' > partitions.json
    fi
    touch versions.yml
    """
}
