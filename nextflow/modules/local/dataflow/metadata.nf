// Select a taxon's genomes (data-channel form).
//
// Runs the metadata stage in a task-local working directory and emits the
// portable selection.tsv (accession + taxonomy + outgroup flag + filename) as a
// channel output -- the hand-off the genome step consumes. No shared workdir.
// Tool flags arrive as task.ext.args from conf/modules.config; publishing is
// configured there too.

process METADATA {
    label 'process_low'
    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("selection.tsv")         , emit: selection
    tuple val(meta), path("outgroup_accession.txt"), emit: outgroup_accession
    path "versions.yml"                            , emit: versions
    // Summary table and plots of the selection and the per-genome GTDB/NCBI
    // lineage table; published, not consumed downstream.
    tuple val(meta), path("report/*")              , emit: report, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} metadata -wd metadata_wd ${args}
    cp metadata_wd/selection.tsv selection.tsv
    cp metadata_wd/outgroup_accession.txt outgroup_accession.txt
    mkdir -p report
    for f in metadata_summary.tsv metadata_summary_number_in_level.png \\
             metadata_summary_number_per_level.png metadata_selected.tsv; do
        if [ -f metadata_wd/\$f ]; then cp metadata_wd/\$f report/\$f; fi
    done

    repgenr versions -wd metadata_wd --versions-out tool_versions.yml
    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    printf 'accession\\tfamily\\tgenus\\tspecies\\tis_outgroup\\tfilename\\tcompleteness\\tcontamination\\n' > selection.tsv
    printf 'GCF_000001.1\\tFam\\tGen\\tsp1\\t0\\tFam_Gen_sp1_GCF_000001.1.fasta\\t\\t\\n' >> selection.tsv
    printf 'GCF_000002.1\\tFam\\tGen\\tsp2\\t0\\tFam_Gen_sp2_GCF_000002.1.fasta\\t\\t\\n' >> selection.tsv
    printf 'GCF_000009.1\\tFam\\tOut\\tgrp\\t1\\tFam_Out_grp_GCF_000009.1.fasta\\t\\t\\n' >> selection.tsv
    printf 'GCF_000009.1\\n' > outgroup_accession.txt
    mkdir -p report
    printf 'level\\tname\\tgenomes\\nfamily\\tFam\\t2\\n' > report/metadata_summary.tsv
    printf 'accession\\taccession_ncbi\\tis_outgroup\\n' > report/metadata_selected.tsv
    touch versions.yml
    """
}
