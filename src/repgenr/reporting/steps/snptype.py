"""Report of the SNP typing: the core-SNP alignment and pairwise distances.

Reads ``snp/alignment_stats.json`` and, when the typer wrote it,
``snp/snp_distance_matrix.tsv``. The heatmap follows the leaf order of
``tree/tree.nwk`` when the tree holds the same genomes, else an average
linkage of the distances.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .. import style
from ..figures import integer_ticks, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import describe, write_table

SNP_DIR = "snp"
STATS = "alignment_stats.json"
MATRIX = "snp_distance_matrix.tsv"
# No heatmap above this many genomes; the histogram is always drawn.
MAX_HEATMAP = 300


def sources(workdir: Path) -> list[Path]:
    return [workdir / SNP_DIR / STATS, workdir / SNP_DIR / MATRIX, workdir / "tree" / "tree.nwk"]


def present(workdir: Path) -> bool:
    return (workdir / SNP_DIR / STATS).is_file() or (workdir / SNP_DIR / MATRIX).is_file()


def read_matrix(path: Path) -> tuple[list[str], list[list[float]]]:
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        header = next(reader)
        names, rows = [], []
        for row in reader:
            if not row:
                continue
            names.append(row[0])
            rows.append([float(v) for v in row[1:]])
    if header[1:] and header[1:] != names:
        raise ValueError(f"{path.name}: row and column names differ")
    return names, rows


def collect(inp: StepInput) -> StepReport:
    snp = inp.workdir / SNP_DIR
    try:
        stats = json.loads((snp / STATS).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        stats = {}
    record = inp.records.get("snptype")
    report = StepReport(step="snptype", title="SNP typing")
    constant = stats.get("constant_sites") or {}
    rows = [
        ("SNP typer", record.tool if record else ""),
        ("reference genome", (record.params.get("reference") if record else "") or ""),
        ("variable columns", stats.get("variable_sites")),
        ("constant sites (A, C, G, T)", ", ".join(str(constant.get(b, "")) for b in "ACGT")),
        ("recombination masked", stats.get("masked")),
        ("core fraction", stats.get("core_fraction")),
    ]
    report.tables.append(
        write_table(
            inp.out_dir,
            "snp_alignment",
            "Summary of the core-SNP alignment.",
            ["quantity", "value"],
            [(k, v) for k, v in rows if v not in (None, "", ", , , ")],
        )
    )
    report.key_numbers = [KeyNumber("Variable columns", stats.get("variable_sites"))]
    if (snp / MATRIX).is_file():
        names, matrix = read_matrix(snp / MATRIX)
        pairs = [matrix[i][j] for i in range(len(names)) for j in range(i + 1, len(names))]
        d = describe(pairs)
        report.key_numbers.append(KeyNumber("Median pairwise SNPs", d["median"]))
        report.summary.append(
            f"The core-SNP alignment of {len(names):,} genomes holds "
            f"{stats.get('variable_sites', 0):,} variable columns; pairwise distances range "
            f"from {d['min']:,.0f} to {d['max']:,.0f} SNPs (median {d['median']:,.0f})."
        )
        report.tables.append(
            write_table(
                inp.out_dir,
                "snp_distances",
                "Pairwise SNP distances between genomes (each pair counted once).",
                ["genomes", "pairs", "min", "q1", "median", "q3", "max"],
                [(len(names), d["n"], d["min"], d["q1"], d["median"], d["q3"], d["max"])],
                numeric=["genomes", "pairs", "min", "q1", "median", "q3", "max"],
            )
        )
        if pairs:
            report.figures.append(_histogram(inp.out_dir, pairs, inp.formats))
        if 2 < len(names) <= MAX_HEATMAP:
            order = _order(inp.workdir, names, matrix)
            report.figures.append(_heatmap(inp.out_dir, names, matrix, order, inp.formats))
        elif len(names) > MAX_HEATMAP:
            report.notes.append(
                f"No distance heatmap is drawn for more than {MAX_HEATMAP} genomes."
            )
    else:
        report.summary.append(
            f"The core-SNP alignment holds {stats.get('variable_sites', 0):,} variable columns. "
            "The SNP typer wrote no distance matrix."
        )
    return report


def _order(workdir: Path, names: list[str], matrix: list[list[float]]) -> list[int]:
    tree_path = workdir / "tree" / "tree.nwk"
    if tree_path.is_file():
        import dendropy

        tree = dendropy.Tree.get(path=str(tree_path), schema="newick", preserve_underscores=True)
        tree.ladderize(ascending=True)
        leaves = [leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon]
        index = {n: i for i, n in enumerate(names)}
        if set(leaves) >= set(names):
            return [index[n] for n in leaves if n in index]
    import numpy as np

    from ...dereplicators.compare_io import average_linkage, leaf_order

    return list(leaf_order(average_linkage(np.asarray(matrix, dtype=float))))


def _histogram(out: Path, pairs: list[float], formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.hist(pairs, bins=40, color=style.PRIMARY, linewidth=0)
        ax.set_xlabel("Pairwise distance (SNPs)")
        ax.set_ylabel("Genome pairs")
        integer_ticks(ax.yaxis)

    files = render(out, "snp_distance_histogram", draw, formats=formats)
    return FigureRef(
        "snp_distance_histogram",
        "Distribution of the pairwise SNP distances.",
        "Histogram of pairwise SNP distances.",
        files,
    )


def _heatmap(out: Path, names, matrix, order, formats) -> FigureRef:
    import numpy as np

    m = np.asarray(matrix)[np.ix_(order, order)]
    labels = [names[i] for i in order]
    n = len(labels)

    def draw(fig) -> None:
        ax = fig.add_subplot()
        image = ax.imshow(m, cmap="viridis", interpolation="nearest")
        if n <= 60:
            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax.set_xticklabels(labels, rotation=90, fontsize=4)
            ax.set_yticklabels(labels, fontsize=4)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        fig.colorbar(image, ax=ax, shrink=0.7, label="SNPs")

    size = min(style.DOUBLE_COLUMN, max(3.0, 0.08 * n + 2.0))
    files = render(out, "snp_distance_heatmap", draw, width=size, height=size, formats=formats)
    return FigureRef(
        "snp_distance_heatmap",
        "Pairwise SNP distances, genomes ordered as the leaves of the tree (or by average "
        "linkage of the distances when no tree holds them).",
        "Heatmap of pairwise SNP distances.",
        files,
    )
