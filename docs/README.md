# RepGenR2 documentation

This directory holds the user and developer documentation for RepGenR2. The
pages are listed in the order a new user is likely to need them: installation,
the choice of tools, running the pipeline, and then reference material. The
project overview is in the [repository README](../README.md).

| Page | Read it when |
|------|--------------|
| [install.md](install.md) | You are setting up RepGenR2: the package, the external tools (core and satellite conda environments, or containers), databases, and how to verify the installation. |
| [choosing-tools.md](choosing-tools.md) | You need to pick tools: which method suits a dataset type and size, the dereplicator and tree builder trade-offs, and each adapter's declared genome limit. |
| [usage.md](usage.md) | You want to run the pipeline from the command line: the stages, resume and `--force`, counting taxa, GTDB, local genomes and viruses, genome sketches, representative selection, inspecting and comparing dereplications, a summary of the phylogeny routes, tree2tax, running tools in containers (quick start), and troubleshooting. |
| [reads.md](reads.md) | You build the genome set from sequencing runs: selecting runs, assembly and polishing, quality scoring and classification, screening runs before or during assembly, and the stateless steps. |
| [phylogeny.md](phylogeny.md) | You choose and tune a phylogeny route: `phylo-advise` and presets, SNP typing and masking, annotation, the core-gene route and its builders, support values, rooting and alignment reuse. |
| [nextflow.md](nextflow.md) | You run the Nextflow pipeline: structure, parameters, profiles and containers, reads mode, scaling, and per-process configuration. |
| [cli-reference.md](cli-reference.md) | You need the exact options of a command. Generated from the command tree. |
| [output.md](output.md) | You want to know what a stage wrote and what each file means, including the columns of the cluster summary. |
| [developing.md](developing.md) | You are changing the code: architecture, data contracts, adding a tool adapter. |
| [verification.md](verification.md) | You want to know which adapters have been run against their real tools, on what data, and how long the measured runs took. |
| [audit/](audit/README.md) | Records behind the current design: the CLI matrix, the 2026 scaling and bias audit, and the 2026-10 gap analysis of the phylo and tree2tax stages. |

The changelog ([CHANGELOG.md](../CHANGELOG.md)) and the citations
([CITATIONS.md](../CITATIONS.md)) are at the repository root.
