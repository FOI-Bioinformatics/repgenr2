// Collect every process's versions.yml fragment into one published file.
//
// Fragments are de-duplicated by content and concatenated into
// <outdir>/pipeline_info/software_versions.yml so a run records which tool
// versions produced it. Scattered tasks of the same process emit identical
// fragments from distinct work directories, so the de-duplication reads each
// fragment's text: unique() on the Path objects themselves would compare
// distinct paths and keep every copy.

workflow PUBLISH_VERSIONS {
    take:
    ch_versions   // channel: versions.yml fragments from every process

    main:
    def ch_yml = ch_versions
        .map { f -> f.text }
        .unique()
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
        )

    emit:
    yml = ch_yml   // the published file, for the run report
}
