# Run report

`repgenr report` writes a report of a run: one HTML page with the results of
each step as numbered tables and figures, a methods paragraph and the
software and parameters used, and one folder per step with the same tables
as TSV and the figures as image files. This page describes what the report
contains, how it is kept up to date, and its limits. The file formats are
listed in [output.md](output.md#report) and the options in
[cli-reference.md](cli-reference.md).

## Contents

1. [Writing the report](#writing-the-report)
2. [The steps](#the-steps)
3. [The HTML page](#the-html-page)
4. [Tables and figures for a manuscript](#tables-and-figures-for-a-manuscript)
5. [Limits](#limits)
6. [Nextflow](#nextflow)

## Writing the report

```bash
repgenr report -wd $WD                       # report/report.html and report/<step>/
repgenr report -wd $WD --steps selection,dereplication --out derep.html
repgenr report -wd $WD --figure-formats png,svg,pdf
```

Each stage refreshes the folders of the steps it feeds when it finishes,
and also when it is skipped because nothing changed, so a working directory
from an earlier version gets its folders on its next run. `repgenr run`
writes the page at the end of the chain. `repgenr report` rebuilds a step
folder only when it is missing or stale: when a file it was built from, the
record of a stage behind it, the figure formats or the repgenr version
changed. `--no-refresh` uses the folders as they are.

A step without results (no reads in a GTDB run, no SNP typing on a
core-gene route) is left out. A step that cannot be reported is named in a
warning and in the page; neither the stage nor the command fails because of
the report.

## The steps

The steps appear in the page in pipeline order:

| Step | Stages | Main content |
|---|---|---|
| `reads` | reads, reads-sketch, reads-screen | runs remaining after each filter; runs per platform; sketch status and screen decisions |
| `assembly` | assemble | contiguity, completeness and coverage per assembler; runs without a genome and why |
| `selection` | metadata, ingest, vgenome, genome, assemble | candidates and selected genomes; genomes per taxon; assembly quality; the outgroup candidates and why the outgroup was chosen |
| `glance` | glance | distribution of the pairwise similarity, with the dereplication thresholds |
| `dereplication` | dereplicate, derep-refine, cluster-summary | status of every genome and the quality-filter reasons; cluster sizes; clusters joining several species; genomes and representatives per species; ANI of each member to its representative (from the sketches) |
| `snptype` | snptype | core-SNP alignment; pairwise SNP distances as a histogram and heatmap; reference coverage per genome (`simple` typer) |
| `alignment` | phylo | source and size of the alignment; occupancy and informative sites of the core genes |
| `phylogeny` | phylo | tree summary (model, support, root-to-tip distances); the tree |
| `tree2tax` | tree2tax | collapsed nodes, nodes per rank, node size, support and RED |

The tables and figures of each step are listed in
[output.md](output.md#report).

## The HTML page

`report/report.html` holds:

- a summary of the key numbers of every step;
- one section per step, with tables numbered (Table 1, 2, ...) and
  captioned above, and figures numbered (Figure 1, 2, ...) and captioned
  below;
- a methods paragraph built from the stage records in `repgenr.yaml`, with
  references for GTDB and NCBI Datasets;
- the stages of the run with their state, the software versions and the
  recorded parameters.

The page is one file: figures (PNG) and the full tables (as TSV downloads,
up to 2 MB each) are embedded, and it loads nothing from the network, so it
can be archived or sent as it is. It prints on A4 with one section per page.
A table shows its first 25 rows in the page.

## Tables and figures for a manuscript

Each folder `report/<step>/` holds the tables as TSV with unformatted
numbers, and each figure as PNG at 300 dpi and SVG; `--figure-formats
png,svg,pdf` adds PDF. The figures share one style: a palette that stays
distinguishable with the common forms of colour-vision deficiency (Okabe
and Ito), DejaVu Sans at 7 to 9 pt, units in the axis labels, and widths of
one journal column (3.35 in) or the full page (6.85 in). In SVG and PDF the
text stays text, so labels can be edited in a drawing program; figure files
carry no date, so the same data give the same files.

## Limits

- The methods paragraph uses only values a stage recorded. A setting that
  was not recorded is left out, not guessed, and it should be checked
  against the parameters in the appendix before it is used in a manuscript.
- The ANI of each member to its representative is estimated from the
  genome sketches (average containment at k=31, scaled=1000), not taken
  from the dereplicator, which does not keep it; it needs `sketches/`, and a
  pair sharing too few hashes has no estimate. It is the quantity
  `dereplicate --tool sourmash` uses, and can differ from the ANI of
  another dereplicator (skani, dRep's ANImf).
- The fraction of the reference each genome covers is reported for the
  `simple` SNP typer only; the other typers do not measure it.
- Large sets are drawn within limits: bar charts show the 30 largest
  groups; the SNP heatmap is drawn for up to 300 genomes; tree leaf labels
  are drawn up to 400 leaves (above that a colour per species remains) and
  no tree figure is drawn above 20,000 leaves.
- The alignment step reads the summary files beside the alignment, not the
  alignment itself.

## Nextflow

The pipeline writes the same report as a final `REPORT` task and publishes
it as `report/` (see [nextflow.md](nextflow.md#parameters)). The pipeline
keeps no `repgenr.yaml`, so that report has no methods paragraph; its
software table comes from `pipeline_info/software_versions.yml`.
`--report false` skips the task.
