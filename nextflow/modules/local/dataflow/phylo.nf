// Build the phylogeny (data-channel form).
//
// Calls the stateless `repgenr phylo-build` step directly on the staged channel
// files: the merged representatives directory provides the genome set, the
// outgroup files are staged into an outgroup/ directory, and tree/tree.nwk is
// emitted as a channel output. Tool flags arrive as task.ext.args from
// conf/modules.config; publishing is configured there too. A staged sketches/
// directory is passed as --sketches-dir; the sourmash tree builder reads it.
// Staged annotations/ (the ANNOTATE records, --msa-source coregene) are
// indexed and passed as --annotations-dir.

process PHYLO {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession), path(annotations, stageAs: 'annotations/*')
    path sketches, stageAs: 'sketches'

    output:
    tuple val(meta), path("tree/tree.nwk")            , emit: tree
    // The alignment the tree was built from and the tree builder's own files
    // (logs, bootstrap trees), published beside the tree. An aligner writes
    // align/; the SNP typing pass (--msa-source snptype) writes tree/msa/,
    // which tree/* includes. An alignment-free builder writes neither.
    tuple val(meta), path("align/*"), optional: true , emit: align
    path "tree/*"                                     , emit: tree_files
    path "versions.yml"                               , emit: versions

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
        -o . -t ${task.cpus} ${args} \$sk \$an \\
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
    if [ -d annotations ]; then
        echo "annotations: staged (\$(ls annotations | wc -l | tr -d ' ') records)"
    else
        echo "annotations: none"
    fi
    mkdir -p tree
    names=\$(ls ${reps_dir}/representatives | sed 's/\\.[^.]*\$//' | paste -sd, -)
    echo "(\${names});" > tree/tree.nwk
    touch versions.yml
    """
}
