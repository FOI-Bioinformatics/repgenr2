// Judge the sketched runs against GTDB before any download (data-channel form).
//
// Runs the standalone form of `repgenr reads-screen` on reads.tsv and the
// sketches directory READS_SKETCH emitted: the three checks of the assemble
// screen (taxon_mismatch, low_match_fraction, duplicate_isolate), with the
// GTDB sketch and lineages and the thresholds in task.ext.args. Emits the
// screen directory (reads_screen.tsv, one screen.json per run with the
// gather output) and passed.txt, the accepted runs, which the acquisition
// uses to restrict the per-run assembly tasks. No shared workdir.

process READS_SCREEN {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_tsv), path(sketches, stageAs: 'sketches')

    output:
    tuple val(meta), path("screen")           , emit: screen
    tuple val(meta), path("screen/*/passed.txt"), emit: passed
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} reads-screen \\
        --reads-tsv ${reads_tsv} \\
        --sketches-dir sketches \\
        -o screen \\
        --jobs ${task.cpus} \\
        ${args} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    variant=\$(ls sketches | head -n 1)
    mkdir -p screen/\$variant
    printf 'run_accession\\ttop_match\\ttop_genus\\tfraction\\tduplicate_of\\tdecision\\treason\\n' > screen/\$variant/reads_screen.tsv
    : > screen/\$variant/passed.txt
    # The stub excuses the second run of the table as a duplicate of the first.
    n=0
    tail -n +2 ${reads_tsv} | cut -f1 | while read -r run; do
        n=\$((n+1))
        if [ "\$n" -eq 2 ]; then
            printf '%s\\tGCF_000000001.1\\tFrancisella\\t0.5\\t%s\\texcused\\tduplicate_isolate: stub\\n' "\$run" "\$first" >> screen/\$variant/reads_screen.tsv
        else
            [ "\$n" -eq 1 ] && first=\$run
            printf '%s\\tGCF_000000001.1\\tFrancisella\\t0.5\\t\\tpass\\t\\n' "\$run" >> screen/\$variant/reads_screen.tsv
            echo "\$run" >> screen/\$variant/passed.txt
        fi
        mkdir -p screen/\$variant/\$run
        printf '{"run_accession": "%s"}\\n' "\$run" > screen/\$variant/\$run/screen.json
    done
    touch versions.yml
    """
}
