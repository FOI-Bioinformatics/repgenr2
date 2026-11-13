// Download selected genomes (data-channel form).
//
// Consumes the selection.tsv channel, downloads the genomes with
// `repgenr genome-fetch` (no manifest, no shared workdir), and emits the genome
// FASTAs (and the outgroup, if any) as channel outputs the dereplication
// scatter consumes.

process GENOME {
    label 'process_medium'
    tag "${meta.id}"
    // Raw genomes are large intermediates (they flow on to dereplication); they
    // are not published. The selected representatives are published downstream.

    input:
    tuple val(meta), path(selection)

    output:
    tuple val(meta), path("out/genomes/*") , emit: genomes
    tuple val(meta), path("out/outgroup/*"), emit: outgroup, optional: true
    // Selected accessions NCBI served no genome for (empty when all arrived).
    tuple val(meta), path("missing_accessions.txt"), emit: missing
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} genome-fetch --selection ${selection} --out out \\
        --versions-out tool_versions.yml
    cp out/missing_accessions.txt missing_accessions.txt

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    """
    echo "ext.args: (none; GENOME takes no tool flags)"
    mkdir -p out/genomes out/outgroup
    : > missing_accessions.txt
    tail -n +2 ${selection} | while IFS=\$'\\t' read -r acc fam gen sp og fname completeness contamination; do
        [ -z "\$fname" ] && continue
        if [ "\$og" = "1" ]; then
            printf '>x\\nACGT\\n' > out/outgroup/\$fname
        else
            printf '>x\\nACGT\\n' > out/genomes/\$fname
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: stub
    END_VERSIONS
    """
}
