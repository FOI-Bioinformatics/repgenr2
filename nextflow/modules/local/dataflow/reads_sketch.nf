// Sketch a chunk of the selected runs streamed from SRA (data-channel form).
//
// Runs the standalone form of `repgenr reads-sketch` on the reads.tsv that
// READS_SELECT emitted, restricted to the accessions of one chunk (the
// subworkflow cuts reads.tsv into chunks of params.reads_sketch_chunk runs,
// so a retry streams one chunk again, not the batch). Each run is streamed
// with fastq-dump and piped into sourmash, so no FASTQ file is written, and
// params.reads_sketch_spots caps the spots read per run (task.ext.args
// carries --max-spots). Each stream is a network-bound fastq-dump and a
// single-threaded sourmash, so at least four run at once. Emits the chunk's
// sketches directory, one <variant>/ (all/ or spots_N/) with a .sig.zip per
// run and an index of the chunk's runs; READS_SKETCH_MERGE joins the chunks.
// No shared workdir.

process READS_SKETCH {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_tsv), val(accessions)

    output:
    tuple val(meta), path("sketches"), emit: sketches
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    // An empty list means every run of reads.tsv: the subworkflow sends one
    // only when reads.tsv has no run, so the task then sketches nothing.
    def selection = accessions ? "--accession-file chunk.txt" : ''
    def chunk = accessions ? "printf '%s\\n' ${accessions.join(' ')} > chunk.txt" : ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    ${chunk}

    repgenr ${opts} reads-sketch \\
        --reads-tsv ${reads_tsv} \\
        ${selection} \\
        -o sketches \\
        --jobs ${Math.max(4, task.cpus)} \\
        --allow-all-failed \\
        ${args} \\
        --versions-out tool_versions.yml

    # The lock file and any working directory of the call are not output.
    # (\\( renders as \( in the command file; a doubled form breaks bash.)
    find sketches -mindepth 2 -maxdepth 2 \\( -name .lock -o -name '.home.*' \\) -exec rm -rf {} +

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def variant = (args =~ /--max-spots\s+(\d+)/) ? "spots_${(args =~ /--max-spots\s+(\d+)/)[0][1]}" : 'all'
    def spots = variant == 'all' ? '' : variant.replace('spots_', '')
    def runs = accessions ? "printf '%s\\n' ${accessions.join(' ')}" : "tail -n +2 ${reads_tsv} | cut -f1"
    """
    echo "ext.args: ${args}"
    mkdir -p sketches/${variant}
    printf 'run_accession\\tstatus\\tmax_spots\\tspots_read\\tsequences\\tsketch_file\\tsketch_params\\tseconds\\treason\\n' > sketches/${variant}/reads_sketches.tsv
    ${runs} | while read -r run; do
        : > sketches/${variant}/\$run.sig.zip
        printf '%s\\tok\\t%s\\t10\\t20\\t%s.sig.zip\\tk=21,k=31,k=51,scaled=1000,abund\\t1.0\\t\\n' "\$run" "${spots}" "\$run" >> sketches/${variant}/reads_sketches.tsv
    done
    touch versions.yml
    """
}
