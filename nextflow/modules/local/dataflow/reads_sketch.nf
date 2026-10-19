// Sketch the selected runs streamed from SRA (data-channel form).
//
// Runs the standalone form of `repgenr reads-sketch` on the reads.tsv that
// READS_SELECT emitted: each run is streamed with fastq-dump and piped into
// sourmash, so no FASTQ file is written, and params.reads_sketch_spots caps
// the spots read per run (task.ext.args carries --max-spots). Emits the
// sketches directory, which holds one <variant>/ (all/ or spots_N/) with a
// .sig.zip per run and the index reads_sketches.tsv. No shared workdir.

process READS_SKETCH {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_tsv)

    output:
    tuple val(meta), path("sketches"), emit: sketches
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} reads-sketch \\
        --reads-tsv ${reads_tsv} \\
        -o sketches \\
        --jobs ${task.cpus} \\
        ${args} \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def variant = (args =~ /--max-spots\s+(\d+)/) ? "spots_${(args =~ /--max-spots\s+(\d+)/)[0][1]}" : 'all'
    def spots = variant == 'all' ? '' : variant.replace('spots_', '')
    """
    echo "ext.args: ${args}"
    mkdir -p sketches/${variant}
    printf 'run_accession\\tstatus\\tmax_spots\\tspots_read\\tsequences\\tsketch_file\\tsketch_params\\tseconds\\treason\\n' > sketches/${variant}/reads_sketches.tsv
    tail -n +2 ${reads_tsv} | cut -f1 | while read -r run; do
        : > sketches/${variant}/\$run.sig.zip
        printf '%s\\tok\\t%s\\t10\\t20\\t%s.sig.zip\\tk=21,k=31,k=51,scaled=1000,abund\\t1.0\\t\\n' "\$run" "${spots}" "\$run" >> sketches/${variant}/reads_sketches.tsv
    done
    touch versions.yml
    """
}
