// Acquire viral genomes (data-channel form): vmetadata + vgenome.
//
// The viral selection/download writes genomes directly (BV-BRC FASTA, length
// filtering, mashtree outgroup), so this module bridges: it runs vmetadata +
// vgenome in a task-local workdir and emits the genomes, the outgroup, and the
// outgroup accession as channel outputs that feed the same downstream data-channel
// dereplication/phylo/tree2tax as the bacterial path.

process VACQUIRE {
    label 'process_medium'
    tag "${meta.id}"
    // Raw genomes are large intermediates (they flow on to dereplication); they
    // are not published. The selected representatives are published downstream.

    input:
    val meta

    output:
    tuple val(meta), path("out/genomes/*")         , emit: genomes
    tuple val(meta), path("out/outgroup/*")        , emit: outgroup, optional: true
    tuple val(meta), path("outgroup_accession.txt"), emit: outgroup_accession
    tuple val(meta), path("segments.tsv")          , emit: segments, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} vmetadata -wd wd ${args}
    repgenr ${opts} vgenome   -wd wd ${args2}

    mkdir -p out
    cp -r wd/genomes out/genomes
    [ -d wd/outgroup ] && cp -r wd/outgroup out/outgroup || true
    if [ -f wd/outgroup_accession.txt ]; then
        cp wd/outgroup_accession.txt outgroup_accession.txt
    else
        : > outgroup_accession.txt
    fi
    # Written by vgenome --group-segments; absent otherwise (optional output).
    if [ -f wd/segments.tsv ]; then
        cp wd/segments.tsv segments.tsv
    fi

    repgenr versions -wd wd --versions-out tool_versions.yml
    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    echo "ext.args: ${args}"
    echo "ext.args2: ${args2}"
    mkdir -p out/genomes
    printf '>x\\nACGT\\n' > out/genomes/Vir_gen_sp1_iso1.fasta
    printf '>x\\nACGT\\n' > out/genomes/Vir_gen_sp2_iso2.fasta
    # Mirror the real stage: --no-outgroup leaves no outgroup dir and an
    # empty accession file, so the optional output path stays unmatched.
    case "${args2}" in
        *--no-outgroup*) : > outgroup_accession.txt ;;
        *)
            mkdir -p out/outgroup
            printf '>x\\nACGT\\n' > out/outgroup/Vir_out_grp_iso9.fasta
            printf 'iso9\\n' > outgroup_accession.txt ;;
    esac
    case "${args2}" in
        *--group-segments*)
            printf 'isolate\\taccession\\tsegment\\tsegment_label\\n' > segments.tsv
            printf 'iso1\\tiso1_S\\tS\\tS\\n' >> segments.tsv ;;
    esac
    touch versions.yml
    """
}
