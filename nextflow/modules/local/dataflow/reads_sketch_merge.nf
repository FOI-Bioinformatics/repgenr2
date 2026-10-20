// Join the sketch chunks of READS_SKETCH into one sketches directory.
//
// Each chunk task emitted a sketches/ with one <variant>/ holding its runs'
// .sig.zip files and an index of those runs. The merge copies the sketches
// of every chunk into one sketches/<variant>/ and writes one index with a
// row per run of reads.tsv in its order (`pending` for a run no chunk
// reached), the layout `reads-sketch` writes for a workdir and the one
// READS_SCREEN reads. Shell and awk only; no tool, no shared workdir. The
// same script serves as the stub.

process READS_SKETCH_MERGE {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads_tsv), path(chunks, stageAs: 'chunk?/sketches')

    output:
    tuple val(meta), path("sketches"), emit: sketches

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p sketches
    for index in chunk*/sketches/*/reads_sketches.tsv; do
        [ -e "\$index" ] || continue
        variant=\$(basename "\$(dirname "\$index")")
        mkdir -p "sketches/\$variant"
        for sketch in "\$(dirname "\$index")"/*.sig.zip; do
            [ -e "\$sketch" ] && cp "\$sketch" "sketches/\$variant/"
        done
    done
    for dir in sketches/*/; do
        [ -d "\$dir" ] || continue
        variant=\$(basename "\$dir")
        awk -F'\\t' -v OFS='\\t' '
            NR == FNR { if (FNR > 1 && \$1 != "") order[++n] = \$1; next }
            FNR == 1 { header = \$0; next }
            \$1 != "" { row[\$1] = \$0 }
            END {
                print header
                for (i = 1; i <= n; i++) {
                    if (order[i] in row) print row[order[i]]
                    else print order[i], "pending", "", "", "", "", "", "", ""
                }
            }' ${reads_tsv} chunk*/sketches/\$variant/reads_sketches.tsv > "sketches/\$variant/reads_sketches.tsv"
    done
    """

    stub:
    """
    mkdir -p sketches
    for index in chunk*/sketches/*/reads_sketches.tsv; do
        [ -e "\$index" ] || continue
        variant=\$(basename "\$(dirname "\$index")")
        mkdir -p "sketches/\$variant"
        for sketch in "\$(dirname "\$index")"/*.sig.zip; do
            [ -e "\$sketch" ] && cp "\$sketch" "sketches/\$variant/"
        done
    done
    for dir in sketches/*/; do
        [ -d "\$dir" ] || continue
        variant=\$(basename "\$dir")
        awk -F'\\t' -v OFS='\\t' '
            NR == FNR { if (FNR > 1 && \$1 != "") order[++n] = \$1; next }
            FNR == 1 { header = \$0; next }
            \$1 != "" { row[\$1] = \$0 }
            END {
                print header
                for (i = 1; i <= n; i++) {
                    if (order[i] in row) print row[order[i]]
                    else print order[i], "pending", "", "", "", "", "", "", ""
                }
            }' ${reads_tsv} chunk*/sketches/\$variant/reads_sketches.tsv > "sketches/\$variant/reads_sketches.tsv"
    done
    """
}
