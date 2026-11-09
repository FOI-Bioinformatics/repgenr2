"""Report of tree2tax: the taxonomy built from the tree.

Reads ``tree2tax.tsv`` (with its ``rank`` column under ``--ranks``),
``tree2tax_nodes.tsv`` (leaves, depth, support and RED per node),
``genomes_map.tsv`` and the tree2tax record.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any

from ...core.contracts import GENOMES_MAP_TSV, TREE2TAX_NODES_TSV, TREE2TAX_TSV
from .. import style
from ..figures import integer_ticks, log_scale, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import describe, write_table


def sources(workdir: Path) -> list[Path]:
    return [workdir / TREE2TAX_TSV, workdir / TREE2TAX_NODES_TSV, workdir / GENOMES_MAP_TSV]


def present(workdir: Path) -> bool:
    return (workdir / TREE2TAX_TSV).is_file()


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with open(path, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def _f(value: str | None) -> float | None:
    if not value or value == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def collect(inp: StepInput) -> StepReport:
    wd, out = inp.workdir, inp.out_dir
    relations = _rows(wd / TREE2TAX_TSV)
    nodes = _rows(wd / TREE2TAX_NODES_TSV)
    mapped = 0
    if (wd / GENOMES_MAP_TSV).is_file():
        mapped = sum(1 for line in (wd / GENOMES_MAP_TSV).read_text().splitlines() if line.strip())
    record = inp.records.get("tree2tax")
    params: dict[str, Any] = record.params if record else {}

    leaves = [n for n in nodes if n.get("n_leaves") == "1"]
    internal = [n for n in nodes if n.get("n_leaves") not in (None, "", "1")]
    report = StepReport(step="tree2tax", title="Taxonomy from the tree (tree2tax)")
    report.key_numbers = [
        KeyNumber("Taxonomy nodes", len(relations)),
        KeyNumber("Internal nodes", len(internal) if nodes else None),
        KeyNumber("Genomes mapped", mapped),
    ]
    text = (
        f"tree2tax wrote {len(relations):,} child-parent relations"
        + (f" ({len(internal):,} internal nodes above {len(leaves):,} leaves)" if nodes else "")
        + f" and mapped {mapped:,} genome accessions to nodes."
    )
    if params.get("collapsed_nodes"):
        text += f" {params['collapsed_nodes']:,} internal nodes were collapsed."
    report.summary.append(text)

    settings = [
        ("collapse below support", params.get("collapse_support")),
        ("collapse below branch length", params.get("collapse_length")),
        ("collapse within diameter", params.get("collapse_diameter")),
        ("internal nodes collapsed", params.get("collapsed_nodes")),
        ("outgroup removed", params.get("remove_outgroup")),
        ("dereplicated genomes included", params.get("include_dereplicated")),
        ("node naming", params.get("naming") or params.get("node_basename")),
        ("rank bands (RED)", params.get("rank_red")),
    ]
    rows = [(k, v) for k, v in settings if v not in (None, "", False)]
    if rows:
        report.tables.append(
            write_table(
                out,
                "tree2tax_settings",
                "Settings of the tree2tax run and the nodes it collapsed.",
                ["setting", "value"],
                rows,
            )
        )
    ranks = Counter(r.get("rank", "") for r in relations if r.get("rank"))
    if ranks:
        report.tables.append(
            write_table(
                out,
                "nodes_per_rank",
                "Nodes per rank assigned by the RED bands (--ranks).",
                ["rank", "nodes"],
                sorted(ranks.items(), key=lambda kv: -kv[1]),
                numeric=["nodes"],
            )
        )
    if internal:
        sizes = [int(n["n_leaves"]) for n in internal if n.get("n_leaves", "").isdigit()]
        d = describe(sizes)
        report.tables.append(
            write_table(
                out,
                "internal_nodes",
                "Leaves below the internal nodes of the taxonomy.",
                ["nodes", "min", "q1", "median", "q3", "max"],
                [(d["n"], d["min"], d["q1"], d["median"], d["q3"], d["max"])],
                numeric=["nodes", "min", "q1", "median", "q3", "max"],
            )
        )
        if len(sizes) >= 3:
            report.figures.append(_sizes_figure(out, sizes, inp.formats))
        supports = [v for v in (_f(n.get("support")) for n in internal) if v is not None]
        if len(supports) >= 3:
            report.figures.append(
                _support_figure(out, supports, params.get("collapse_support"), inp.formats)
            )
        red = [(v, n.get("rank", "")) for n in internal if (v := _f(n.get("red"))) is not None]
        if len(red) >= 3:
            report.figures.append(_red_figure(out, red, inp.formats))
    return report


def _sizes_figure(out: Path, sizes: list[int], formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        if max(sizes) > 50:
            import numpy as np

            bins = np.unique(np.logspace(0, np.log10(max(sizes) + 1), 30).astype(int))
            ax.hist(sizes, bins=bins, color=style.PRIMARY, linewidth=0)
            log_scale(ax, "x")
            ax.set_xlabel("Leaves below the node (log scale)")
        else:
            ax.hist(sizes, bins=range(1, max(sizes) + 2), color=style.PRIMARY, linewidth=0)
            ax.set_xlabel("Leaves below the node")
            integer_ticks(ax.xaxis)
        ax.set_ylabel("Internal nodes")
        integer_ticks(ax.yaxis)

    files = render(out, "node_sizes", draw, formats=formats)
    return FigureRef(
        "node_sizes",
        "Number of leaves below each internal node of the taxonomy.",
        "Histogram of node sizes.",
        files,
    )


def _support_figure(out: Path, supports: list[float], threshold, formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.hist(supports, bins=20, color=style.PRIMARY, linewidth=0)
        if isinstance(threshold, int | float):
            # tree2tax reads a threshold below 1 against percent supports as a fraction.
            value = threshold * 100 if threshold <= 1 < max(supports) else threshold
            ax.axvline(value, color="black", lw=0.6, ls="--")
        ax.set_xlabel("Support of the node")
        ax.set_ylabel("Internal nodes")
        integer_ticks(ax.yaxis)

    files = render(out, "node_support", draw, formats=formats)
    caption = "Support values of the internal nodes kept in the taxonomy."
    if isinstance(threshold, int | float):
        caption += " Dashed line: the collapse threshold."
    return FigureRef("node_support", caption, "Histogram of node support.", files)


def _red_figure(out: Path, red: list[tuple[float, str]], formats) -> FigureRef:
    ranks = sorted({r for _, r in red if r and r != "no rank"})

    def draw(fig) -> None:
        ax = fig.add_subplot()
        if ranks:
            groups = [[v for v, r in red if r == rank] for rank in ranks]
            other = [v for v, r in red if r not in ranks]
            if other:
                groups.append(other)
            labels = ranks + (["no rank"] if other else [])
            ax.hist(
                groups,
                bins=20,
                range=(0, 1),
                stacked=True,
                color=[style.color(i) for i in range(len(groups))],
                label=labels,
                linewidth=0,
            )
            ax.legend(title="rank")
        else:
            ax.hist([v for v, _ in red], bins=20, range=(0, 1), color=style.PRIMARY, linewidth=0)
        ax.set_xlabel("Relative evolutionary divergence (RED)")
        ax.set_ylabel("Internal nodes")
        integer_ticks(ax.yaxis)

    files = render(out, "node_red", draw, formats=formats)
    return FigureRef(
        "node_red",
        "Relative evolutionary divergence of the internal nodes"
        + (", by the rank assigned." if ranks else "."),
        "Histogram of RED values.",
        files,
    )
