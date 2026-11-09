// Helpers that gather the files of the run report (modules/local/dataflow/report.nf).

// [meta, files] for every meta of ch_keys: the files of the optional output
// ch_optional for that meta, or an empty list when the process wrote none.
def optional_files(ch_keys, ch_optional) {
    return ch_keys
        .map { meta, _x -> tuple(meta, []) }
        .join(ch_optional, by: 0, remainder: true)
        .map { meta, _placeholder, files ->
            tuple(meta, files == null ? [] : (files instanceof List ? files : [files]))
        }
}

// [meta, [tree2tax.tsv, genomes_map.tsv, tree2tax_nodes.tsv]] from the three
// TREE2TAX outputs.
def tree2tax_files(ch_relations, ch_map, ch_nodes) {
    return ch_relations
        .join(ch_map, by: 0)
        .join(ch_nodes, by: 0)
        .map { meta, relations, map, nodes -> tuple(meta, [relations, map, nodes]) }
}
