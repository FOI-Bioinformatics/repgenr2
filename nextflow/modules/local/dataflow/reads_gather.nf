// Gather the per-run assemblies into the genome contract (data-channel form).
//
// `repgenr reads-gather` reads reads.tsv, the assemble-run directories and the
// optional genome-qc directory, applies the quality gate and the naming
// policy (GTDB tokens where the classifier agrees at genus), and writes
// genomes/, selection.tsv, assembly_stats.tsv and excused_runs.tsv, the
// selection summary (metadata_summary.tsv and its plots), and
// screen_reads.tsv when assemble-run screened the runs (--screen-reads in
// params.assemble_args; one screen.json per run directory). The
// outputs mirror the ACQUIRE subworkflow's: the genome FASTAs feed the
// dereplication scatter, selection.tsv the quality-aware keeper, and the
// (empty) outgroup accession file the phylogeny and tree2tax steps.

process READS_GATHER {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_tsv), path(assemblies, stageAs: 'assemblies/*'), path(qc, stageAs: 'qc')

    output:
    tuple val(meta), path("out/genomes/*")             , emit: genomes
    tuple val(meta), path("out/selection.tsv")         , emit: selection
    tuple val(meta), path("out/assembly_stats.tsv")    , emit: assembly_stats
    tuple val(meta), path("out/excused_runs.tsv")      , emit: excused, optional: true
    tuple val(meta), path("out/screen_reads.tsv")      , emit: screen, optional: true
    tuple val(meta), path("out/outgroup_accession.txt"), emit: outgroup_accession
    tuple val(meta), path("out/metadata_summary*")     , emit: summary, optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    # The qc directory is optional: without a database GENOME_QC does not run
    # and nothing is staged under that name.
    qc=""
    [ -d qc ] && qc="--qc qc"

    repgenr ${opts} reads-gather \\
        --reads-tsv ${reads_tsv} \\
        --assemblies assemblies \\
        --out out \\
        ${args} \\
        \$qc

    repgenr_versions_fragment "${task.process}"
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    mkdir -p out/genomes
    printf 'accession\\tfamily\\tgenus\\tspecies\\tis_outgroup\\tfilename\\tcompleteness\\tcontamination\\n' > out/selection.tsv
    printf 'run_accession\\tfilename\\tassembler\\tn_contigs\\ttotal_length\\tn50\\tlargest_contig\\test_coverage\\tcompleteness\\tcontamination\\tncbi_taxonomy\\tgtdb_taxonomy\\tlabel_source\\ttaxonomy_flag\\tpolisher\\treads_sketch\\tscreen\\n' > out/assembly_stats.tsv
    : > excused.tmp
    tail -n +2 ${reads_tsv} | while IFS=\$'\\t' read -r run rest; do
        fam=\$(echo "\$rest" | cut -f10); gen=\$(echo "\$rest" | cut -f11); sp=\$(echo "\$rest" | cut -f12)
        if [ -e "assemblies/\$run/assembly.ok" ]; then
            name="\${fam}_\${gen}_\${sp}_\${run}.fasta"
            cp assemblies/\$run/contigs.fasta out/genomes/\$name
            printf '%s\\t%s\\t%s\\t%s\\t0\\t%s\\t\\t\\n' "\$run" "\$fam" "\$gen" "\$sp" "\$name" >> out/selection.tsv
            printf '%s\\t%s\\tstub\\t1\\t12\\t12\\t12\\t\\t\\t\\t%s;%s;%s\\t\\tmetadata\\t\\t\\t0\\t\\n' "\$run" "\$name" "\$fam" "\$gen" "\$sp" >> out/assembly_stats.tsv
        elif [ -e "assemblies/\$run/excused_runs.tsv" ]; then
            tail -n +2 assemblies/\$run/excused_runs.tsv >> excused.tmp
        fi
    done
    if [ -s excused.tmp ]; then
        printf 'run_accession\\tstep\\treason\\n' > out/excused_runs.tsv
        cat excused.tmp >> out/excused_runs.tsv
    fi
    for f in assemblies/*/screen.json; do
        [ -e "\$f" ] || continue
        [ -e out/screen_reads.tsv ] || printf 'run_accession\\ttop_match\\ttop_genus\\tfraction\\tduplicate_of\\tdecision\\treason\\n' > out/screen_reads.tsv
        run=\$(basename \$(dirname "\$f"))
        decision=\$(sed -n 's/.*"decision": *"\\([^"]*\\)".*/\\1/p' "\$f")
        printf '%s\\t\\t\\t\\t\\t%s\\t\\n' "\$run" "\$decision" >> out/screen_reads.tsv
    done
    : > out/outgroup_accession.txt
    if [ -n "\$(ls out/genomes)" ]; then
        printf 'level\\tname\\tgenomes\\n' > out/metadata_summary.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: stub
    END_VERSIONS
    """
}
