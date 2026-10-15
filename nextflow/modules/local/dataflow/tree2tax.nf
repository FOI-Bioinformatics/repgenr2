// Turn the rooted tree into FlexTaxD relations (data-channel form).
//
// Calls the stateless `repgenr tree2tax-relations` step directly on the staged
// channel files: the tree, the merged representatives' clusters.tsv (for
// --include-dereplicated), the outgroup directory and, on the viral path, the
// segments.tsv of grouped isolates (an empty list elsewhere). Emits
// tree2tax.tsv and genomes_map.tsv as channel outputs; there is no shared
// working directory.
// Tool flags arrive as task.ext.args from conf/modules.config; publishing is
// configured there too.

process TREE2TAX {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tree), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession), path(segments)

    output:
    tuple val(meta), path("tree2tax.tsv")   , emit: tree2tax
    tuple val(meta), path("genomes_map.tsv"), emit: genomes_map
    tuple val(meta), path("tree2tax_nodes.tsv"), emit: nodes
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    def segments_arg = segments ? "--segments ${segments}" : ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} tree2tax-relations \\
        --tree ${tree} \\
        --clusters ${reps_dir}/clusters.tsv \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} ${segments_arg} \\
        -o . ${args} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    if [ -f segments.tsv ]; then echo "segments: staged"; else echo "segments: none"; fi
    printf 'child\\tparent\\n' > tree2tax.tsv
    for f in ${reps_dir}/representatives/*; do
        leaf=\$(basename \$f | sed 's/\\.[^.]*\$//')
        printf '%s\\troot\\n' "\$leaf" >> tree2tax.tsv
    done
    printf 'node\\tparent\\tn_leaves\\tdepth\\tbranch_length\\tsupport\\tred\\tlineage\\n' > tree2tax_nodes.tsv
    : > genomes_map.tsv
    for f in ${reps_dir}/representatives/*; do
        leaf=\$(basename \$f | sed 's/\\.[^.]*\$//')
        acc=\$(echo \$leaf | awk -F_ '{print \$(NF-1)"_"\$NF}')
        printf '%s\\t%s\\n' "\$acc" "\$leaf" >> genomes_map.tsv
    done
    # As the real step: a staged segments.tsv adds its member accessions.
    if [ -f segments.tsv ]; then
        tail -n +2 segments.tsv | awk -F'\\t' '{print \$2"\\t"\$1}' >> genomes_map.tsv
    fi
    touch versions.yml
    """
}
