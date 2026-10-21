// Build the tree from a ready alignment (data-channel form).
//
// The second half of `repgenr phylo-build`. It takes the alignment PHYLO_MSA
// produced, so changing the tree builder or the bootstrap re-runs this task
// alone. The genome set is still staged: the tree builder roots on the
// outgroup named in the accession file. A staged sketches/ directory is
// passed as --sketches-dir, as in PHYLO. alignment_stats.json and
// partitions.json, when PHYLO_MSA wrote them, are staged beside the
// alignment, where phylo-build --msa reads them.

process PHYLO_TREE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession), path(msa), path(stats), path(partitions)
    path sketches, stageAs: 'sketches'

    output:
    tuple val(meta), path("tree/tree.nwk"), emit: tree
    path "tree/*"                         , emit: tree_files  // the builder's logs and bootstrap trees
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    sk=""
    [ -d sketches ] && sk="--sketches-dir sketches"

    repgenr ${opts} phylo-build \\
        --genomes-dir ${reps_dir}/representatives \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        --msa ${msa} \\
        -o . -t ${task.cpus} ${args} \$sk \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    echo "ext.args: ${args}"
    echo "ext.repgenr_opts: ${opts}"
    echo "msa: ${msa}"
    if [ -e partitions.json ]; then echo "partitions: staged"; else echo "partitions: none"; fi
    if [ -d sketches ]; then echo "sketches: staged"; else echo "sketches: none"; fi
    mkdir -p tree
    names=\$(grep '^>' ${msa} | sed 's/^>//' | paste -sd, -)
    echo "(\${names});" > tree/tree.nwk
    touch versions.yml
    """
}
