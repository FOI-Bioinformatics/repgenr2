"""Report of the alignment the tree was built from.

Reads the sidecars in ``tree/msa/`` (``msa_source.json``,
``alignment_stats.json``, and for the core-gene route
``core_gene_stats.json`` and ``core_genes.tsv``); the alignment itself is
not read, so the step stays fast on a large whole-genome alignment.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .. import style
from ..figures import integer_ticks, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import write_table

MSA_DIR = Path("tree") / "msa"
FILES = ("msa_source.json", "alignment_stats.json", "core_gene_stats.json", "core_genes.tsv")


def sources(workdir: Path) -> list[Path]:
    return [workdir / MSA_DIR / name for name in FILES]


def present(workdir: Path) -> bool:
    return (workdir / MSA_DIR / "msa_source.json").is_file()


def _json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def collect(inp: StepInput) -> StepReport:
    msa = inp.workdir / MSA_DIR
    stamp = _json(msa / "msa_source.json")
    stats = _json(msa / "alignment_stats.json")
    core = _json(msa / "core_gene_stats.json")
    record = inp.records.get("phylo")
    params = record.params if record else {}
    source = params.get("msa_source") or ("coregene" if core else "")
    tool = params.get("coregene") or params.get("snptyper") or params.get("aligner") or ""
    genomes = len(stamp.get("records") or []) or core.get("genomes")

    report = StepReport(step="alignment", title="Alignment")
    rows: list[tuple[str, Any]] = [
        ("alignment source", source),
        ("built with", tool),
        ("genomes", genomes),
    ]
    if core:
        rows += [
            ("alphabet", core.get("alphabet")),
            ("alignment columns", core.get("alignment_length")),
            ("variable sites", core.get("variable_sites_total")),
            ("parsimony-informative sites", core.get("informative_sites_total")),
            ("core genes available", core.get("core_genes_available")),
            ("core genes kept", core.get("core_genes")),
            ("core threshold", core.get("core_threshold")),
            ("clustering identity", core.get("clustering_identity")),
            ("gene trimming", core.get("trim")),
        ]
    if stats:
        constant = stats.get("constant_sites") or {}
        rows += [
            ("variable columns", stats.get("variable_sites")),
            ("constant sites (A, C, G, T)", ", ".join(str(constant.get(b, "")) for b in "ACGT")),
            ("recombination masked", stats.get("masked")),
            ("core fraction", stats.get("core_fraction")),
        ]
    report.tables.append(
        write_table(
            inp.out_dir,
            "alignment_summary",
            "Summary of the alignment the tree was built from.",
            ["quantity", "value"],
            [(k, v) for k, v in rows if v not in (None, "")],
        )
    )
    report.key_numbers = [KeyNumber("Aligned genomes", genomes)]
    if core:
        report.key_numbers += [
            KeyNumber("Core genes kept", core.get("core_genes")),
            KeyNumber("Informative sites", core.get("informative_sites_total")),
        ]
        report.summary.append(
            f"The tree was built from {core.get('core_genes', 0):,} of "
            f"{core.get('core_genes_available', 0):,} core genes ({tool}), "
            f"{core.get('alignment_length', 0):,} columns with "
            f"{core.get('informative_sites_total', 0):,} parsimony-informative sites."
        )
        genes = _genes(msa / "core_genes.tsv")
        if genes:
            report.figures += _gene_figures(inp.out_dir, genes, core, inp.formats)
    elif stats:
        report.key_numbers.append(KeyNumber("Variable columns", stats.get("variable_sites")))
        report.summary.append(
            f"The alignment ({source}, {tool}) holds {stats.get('variable_sites', 0):,} "
            "variable columns."
        )
    else:
        report.summary.append(f"The alignment came from {source or 'an unrecorded source'}.")
    return report


def _genes(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def _f(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _gene_figures(out: Path, genes, core, formats) -> list[FigureRef]:
    kept = [g for g in genes if g.get("kept") == "true"]
    other = [g for g in genes if g.get("kept") != "true"]
    threshold = core.get("core_threshold")
    figs: list[FigureRef] = []

    occ = [v for v in (_f(g.get("occupancy")) for g in genes) if v is not None]
    if occ:

        def draw_occ(fig) -> None:
            ax = fig.add_subplot()
            ax.hist(occ, bins=20, range=(0, 1), color=style.PRIMARY, linewidth=0)
            if isinstance(threshold, int | float):
                ax.axvline(float(threshold), color="black", lw=0.6, ls="--")
            ax.set_xlabel("Occupancy (fraction of genomes with the gene)")
            ax.set_ylabel("Genes")
            integer_ticks(ax.yaxis)

        files = render(out, "gene_occupancy", draw_occ, formats=formats)
        caption = "Occupancy of the core genes the builder reported"
        if isinstance(threshold, int | float):
            caption += f"; dashed line: the core threshold ({threshold:g})"
        figs.append(
            FigureRef("gene_occupancy", caption + ".", "Histogram of gene occupancy.", files)
        )

    ranked = sorted(
        (v for v in (_f(g.get("informative_sites")) for g in kept) if v is not None), reverse=True
    )
    if ranked:

        def draw_rank(fig) -> None:
            ax = fig.add_subplot()
            ax.plot(range(1, len(ranked) + 1), ranked, color=style.PRIMARY, lw=1.0)
            ax.set_xlabel("Kept gene, by rank")
            ax.set_ylabel("Parsimony-informative sites")
            integer_ticks(ax.xaxis)

        files = render(out, "informative_by_rank", draw_rank, formats=formats)
        figs.append(
            FigureRef(
                "informative_by_rank",
                "Parsimony-informative sites of each kept core gene, in decreasing order.",
                "Line chart of informative sites by gene rank.",
                files,
            )
        )

    def pts(rows):
        out_pts = []
        for g in rows:
            x, y = _f(g.get("length")), _f(g.get("informative_sites"))
            if x is not None and y is not None:
                out_pts.append((x, y))
        return out_pts

    kp, op = pts(kept), pts(other)
    if kp and op:

        def draw_scatter(fig) -> None:
            ax = fig.add_subplot()
            ax.scatter(
                [p[0] for p in op], [p[1] for p in op], s=6, color=style.MUTED, linewidths=0,
                label="not kept",
            )  # fmt: skip
            ax.scatter(
                [p[0] for p in kp], [p[1] for p in kp], s=6, color=style.PRIMARY, linewidths=0,
                label="kept",
            )  # fmt: skip
            ax.set_xlabel("Gene alignment length (columns)")
            ax.set_ylabel("Parsimony-informative sites")
            ax.legend(loc="upper left")

        files = render(out, "gene_length_informative", draw_scatter, formats=formats)
        figs.append(
            FigureRef(
                "gene_length_informative",
                "Alignment length and parsimony-informative sites of each core gene, kept "
                "and not kept by the gene cap.",
                "Scatter plot of informative sites against gene length.",
                files,
            )
        )
    return figs
