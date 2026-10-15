# Audit records

Reference material that informed the current design. It is kept for the
citations and measurements it holds; the conclusions that still apply have
been folded into the code and the user pages.

| File | Content |
|------|---------|
| [cli-matrix.md](cli-matrix.md) | Every command and flag with the parameter it reaches, how it is validated, and the test and document that cover it. Rendered from `tests/audit/cli_matrix.yaml`; the test suite keeps the two in sync. |
| [scaling-audit.md](scaling-audit.md) | The scaling and bias audit (2026-08 to 2026-09): measured wall time and memory per tool and genome count, bias experiments on over-represented genotypes, and the per-stage SWOT and gap analyses as appendices. |
| [phylo-tree2tax-gap-2026-10.md](phylo-tree2tax-gap-2026-10.md) | Gap analysis of the phylo and tree2tax stages against current practice (2026-10): ascertainment bias, soft core, reference choice, Gubbins defaults, FlexTaxD contract, clade thresholds, node naming; the improvement ladder in batches A to E. |
