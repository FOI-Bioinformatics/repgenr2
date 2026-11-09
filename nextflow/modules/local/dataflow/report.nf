// Write the report of a run: per-step tables and figures and one HTML page.
//
// The pipeline keeps no shared working directory, so the files the report
// reads are staged into one here (wd/): the selection, the selection summary,
// the dereplication directory, the tree, the tree2tax tables and, in reads
// mode, the assembly tables. `repgenr report` then writes wd/report/, which is
// published as report/. Every input but the selection may be an empty list;
// a step whose files are absent is left out of the report.
// Failures are ignored (errorStrategy in conf/modules.config): the report is
// written beside the results and must not fail the run.

process REPORT {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(selection, stageAs: 'wd/selection.tsv'), path(summary, stageAs: 'wd/*'), path(derep, stageAs: 'wd/derep'), path(tree, stageAs: 'wd/tree/tree.nwk'), path(tree2tax, stageAs: 'wd/*'), path(assembly, stageAs: 'wd/*')
    path versions

    output:
    tuple val(meta), path("report"), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    def versions_arg = versions ? "--versions-yml ${versions}" : ''
    """
    repgenr ${opts} report -wd wd --title "RepGenR2 report: ${meta.id}" ${versions_arg} ${args}
    mv wd/report report
    """

    stub:
    """
    mkdir -p report/selection
    echo '<!doctype html><title>stub</title>' > report/report.html
    echo '{}' > report/selection/step.json
    """
}
