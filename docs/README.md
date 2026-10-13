# RepGenR2 documentation

| Page | Read it when |
|------|--------------|
| [install.md](install.md) | You are setting up RepGenR2: the package, the external tools (one conda environment, several environments, or containers), databases, and how to verify the installation. |
| [choosing-tools.md](choosing-tools.md) | You need to pick tools: which method suits a dataset type and size, the dereplicator and tree builder trade-offs, and each adapter's declared genome limit. |
| [usage.md](usage.md) | You want to run the pipeline: the CLI stages, starting from local genomes, viruses, resume and `--force`, representative selection, inspecting and comparing dereplications, SNP typing and masking, tree2tax collapsing, the Nextflow layer (parameters, profiles, scaling, per-process configuration), running tools in containers (quick start), and troubleshooting. |
| [cli-reference.md](cli-reference.md) | You need the exact options of a command. Generated from the command tree. |
| [output.md](output.md) | You want to know what a stage wrote and what each file means, including the columns of the cluster summary. |
| [developing.md](developing.md) | You are changing the code: architecture, data contracts, adding a tool adapter. |
| [verification.md](verification.md) | You want to know which adapters have been run against their real tools, on what data, and how long the measured runs took. |
| [audit/](audit/README.md) | Records behind the current design: the CLI matrix and the 2026 scaling and bias audit. |

The changelog is at the repository root ([CHANGELOG.md](../CHANGELOG.md)).
