# RepGenR2 documentation

The pages are listed in the order a new user is likely to need them. The
project overview is in the [repository README](../README.md).

1. [install.md](install.md): the package, the external tools in conda
   environments or containers, databases, and checking the installation.
2. [choosing-tools.md](choosing-tools.md): which entry path, dereplicator,
   phylogeny route and tree builder suit a dataset, and each adapter's
   declared genome limit.
3. [usage.md](usage.md): running the stages from the command line, from the
   entry paths through dereplication and tree2tax to troubleshooting.
4. [reads.md](reads.md): building the genome set from sequencing runs.
5. [phylogeny.md](phylogeny.md): the phylogeny routes, presets, typers,
   core-gene builders, tree-builder options, supports and rooting.
6. [nextflow.md](nextflow.md): the Nextflow pipeline, its parameters,
   profiles and scaling.
7. [cli-reference.md](cli-reference.md): every command and option,
   generated from the command tree.
8. [output.md](output.md): the files each stage writes and their columns.
9. [developing.md](developing.md): architecture, data contracts and adding
   a tool adapter.
10. [verification.md](verification.md): which adapters have been run against
    their real tools, and the measured runs.
11. [audit/](audit/README.md): the records behind the current design.

The changelog ([CHANGELOG.md](../CHANGELOG.md)) and the citations
([CITATIONS.md](../CITATIONS.md)) are at the repository root.
