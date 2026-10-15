// Annotate one genome for the core-gene route (data-channel form).
//
// One task per genome calls the stateless `repgenr annotate-build`, which
// writes annotations/<record name>/ with the GFF3 (sequences appended), the
// protein (.faa) and gene (.ffn) files. The record name is the genome file
// name without its FASTA suffix, as in repgenr.core.contracts.record_name, so
// the core-gene alignment names its records after the genomes. The tasks run
// in parallel and cache one genome each, so a re-run annotates only what
// changed. The annotator and its arguments arrive as task.ext.args from
// conf/modules.config. A staged bakta database directory is passed as
// --bakta-db; without one, bakta reads BAKTA_DB from the task environment.

process ANNOTATE {
    tag "${meta.id}:${genome.name}"
    label 'process_medium'

    input:
    tuple val(meta), path(genome)
    path bakta_db, stageAs: 'bakta_db'

    output:
    tuple val(meta), path("annotations/*", type: 'dir'), emit: annotation
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    db=""
    [ -d bakta_db ] && db="--bakta-db bakta_db"

    repgenr ${opts} annotate-build \\
        --genome ${genome} \\
        -o annotations -t ${task.cpus} ${args} \$db \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    if [ -d bakta_db ]; then echo "bakta_db: staged"; else echo "bakta_db: none"; fi
    b=\$(basename "${genome}")
    name="\$b"
    for s in .fasta.gz .fna.gz .fa.gz .fasta .fa .fna .fas; do
        case "\$b" in
            *"\$s") name="\${b%"\$s"}"; break ;;
        esac
    done
    mkdir -p "annotations/\$name"
    printf '##gff-version 3\\n##FASTA\\n>%s\\nACGT\\n' "\$name" > "annotations/\$name/\$name.gff3"
    printf '>%s_00001\\nMK\\n' "\$name" > "annotations/\$name/\$name.faa"
    printf '>%s_00001\\nATGAAA\\n' "\$name" > "annotations/\$name/\$name.ffn"
    printf '"%s":\\n    repgenr: stub\\n' "${task.process}" > versions.yml
    """
}
