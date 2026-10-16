// Refine the lineages of the seed genomes at finer ANIs (after the gather step).
//
// Data-channel module: the merged result directory from DEREP_MERGE is staged
// in together with every genome of the run (the merge holds only the
// representatives, and a lineage's members are needed too), and the refined
// result is emitted in the same contract shape, with refined.tsv added. Wraps
// `repgenr dereplicate-refine`. The shared dereplication flags arrive as
// task.ext.args, the seeds and ANI as task.ext.args2, from conf/modules.config.
// A staged refine.tsv (params.derep_refine_file) is passed as --refine-file.

process DEREP_REFINE {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(merged, stageAs: 'merged'), path(genomes, stageAs: 'genomes/*')
    path selection, stageAs: 'selection.tsv'
    path sketches, stageAs: 'sketches'
    path refine_file, stageAs: 'refine.tsv'

    output:
    tuple val(meta), path("${task.ext.prefix ?: meta.id}"), emit: reps
    path 'versions.yml'                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    # Every genome of the run, as a file-of-filenames. printf is a shell
    # builtin, so the list never passes through an argv (ARG_MAX).
    printf '%s\n' genomes/* > genomes.fofn

    sel=""
    [ -e selection.tsv ] && sel="--selection-tsv selection.tsv"
    sk=""
    [ -d sketches ] && sk="--sketches-dir sketches"
    rf=""
    [ -e refine.tsv ] && rf="--refine-file refine.tsv"

    repgenr ${opts} dereplicate-refine \\
        --contract-dir merged \\
        --genomes-fofn genomes.fofn \\
        --out ${prefix} \\
        ${args} ${args2} \\
        \$sel \$sk \$rf \\
        --threads ${task.cpus} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "ext.args: ${args}"
    echo "ext.args2: ${args2}"
    if [ -e refine.tsv ]; then echo "refine_file: staged"; else echo "refine_file: none"; fi
    mkdir -p ${prefix}
    cp -r merged/representatives ${prefix}/representatives
    cp merged/clusters.tsv ${prefix}/clusters.tsv
    cp merged/genome_status.tsv ${prefix}/genome_status.tsv
    printf 'genome\\tseed\\tcoarse_representative\\trepresentative\\tstatus\\tani\\n' > ${prefix}/refined.tsv
    touch versions.yml
    """
}
