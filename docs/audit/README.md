# Audit records

This directory holds the records that informed the current design. They are
kept for the citations and measurements they contain. The conclusions that
still apply have been incorporated into the code and the user pages, so a user
of RepGenR2 does not need to read these records.

| File | Date | Status | Content |
|------|------|--------|---------|
| [cli-matrix.md](cli-matrix.md) | current | kept in sync by tests | Every command and flag with the parameter it reaches, how it is validated, and the test and document that cover it. Rendered from `tests/audit/cli_matrix.yaml`; the test suite keeps the two in sync. |
| [scaling-audit.md](scaling-audit.md) | 2026-08 to 2026-09 | complete; later changes noted in its status paragraph | The scaling and bias audit (2026-08 to 2026-09): measured wall time and memory per tool and genome count, bias experiments on over-represented genotypes, and the per-stage SWOT and gap analyses as appendices. |
| [phylo-tree2tax-gap-2026-10.md](phylo-tree2tax-gap-2026-10.md) | 2026-10 | batches A to D shipped; E1 deferred; two defaults open | Gap analysis of the phylo and tree2tax stages against current practice (2026-10): ascertainment bias, soft core, reference choice, Gubbins defaults, FlexTaxD contract, clade thresholds, node naming; the improvement ladder in batches A to E. |
| [phylo-methods-2026-10.md](phylo-methods-2026-10.md) | 2026-10 | items A to F and H shipped; G pending cluster time | Audit of the phylo MSA sources and tree builders for 100 to 500 genomes (2026-10): fit and gaps of each route, strengths and weaknesses of each builder, defects found, recommendations by purpose and divergence, and the planned improvement items. |

PR numbers below about 270 in the August and September records refer to the earlier fork, druvus/repgenr. Later numbers refer to FOI-Bioinformatics/repgenr2.
