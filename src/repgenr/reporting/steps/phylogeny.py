"""Report of the tree: summary, support, root-to-tip distances and a figure.

Reads ``tree/tree.nwk``, the phylo record (builder, rooting, outgroup) and,
for leaf labels, ``selection.tsv``. The substitution model is taken from an
IQ-TREE report (``*.iqtree``) or a RAxML-NG ``*.bestModel`` when present.
"""

from __future__ import annotations

import re
from pathlib import Path

from ...core.contracts import SELECTION_TSV, glob_visible, read_selection, record_name
from ...tree.support import parse_support
from .. import style
from ..figures import integer_ticks, render
from ..model import FigureRef, KeyNumber, StepReport
from ..registry import StepInput
from ..tables import describe, write_table
from ..treeplot import MAX_LABELLED, MAX_LEAVES, draw_tree, figure_height
from .dereplication import species_label

TREE = Path("tree") / "tree.nwk"
# Branch-length unit per tree builder.
UNITS = {
    "iqtree": "substitutions per site",
    "raxmlng": "substitutions per site",
    "fasttree": "substitutions per site",
    "veryfasttree": "substitutions per site",
    "aster": "coalescent units",
    "mashtree": "Mash distance",
}


def sources(workdir: Path) -> list[Path]:
    tree_dir = workdir / "tree"
    found = glob_visible(tree_dir, "*.iqtree") + glob_visible(tree_dir, "*.bestModel")
    return [workdir / TREE, workdir / SELECTION_TSV, *found]


def present(workdir: Path) -> bool:
    return (workdir / TREE).is_file()


def _load(path: Path):
    import dendropy

    return dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True)


def _model(workdir: Path) -> str:
    tree_dir = workdir / "tree"
    for path in glob_visible(tree_dir, "*.iqtree"):
        text = path.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"Best-fit model according to \w+: (\S+)", text) or re.search(
            r"Model of substitution: (\S+)", text
        )
        if m:
            return m.group(1)
        if "partition" in text.lower():
            return "partitioned (see the IQ-TREE report)"
    for path in glob_visible(tree_dir, "*.bestModel"):
        lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
        # "GTR{1.4/6.3/...}+FU{...}+G4m{0.34}+BU{0.92}, gene = 1-1470": the
        # model with its estimated values left out.
        models = {re.sub(r"\{[^}]*\}", "", ln.split(",", 1)[0]).strip() for ln in lines}
        if len(lines) > 1:
            return f"{', '.join(sorted(models))} per partition ({len(lines)} partitions)"
        if lines:
            return next(iter(models))
    return ""


def collect(inp: StepInput) -> StepReport:
    wd, out = inp.workdir, inp.out_dir
    tree = _load(wd / TREE)
    record = inp.records.get("phylo")
    params = record.params if record else {}
    builder = (record.tool if record else None) or ""
    outgroup = params.get("outgroup")
    selection = (
        {record_name(r.filename): r for r in read_selection(wd / SELECTION_TSV)}
        if (wd / SELECTION_TSV).is_file()
        else {}
    )
    leaves = [leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon is not None]
    labels: dict[str, str] = {}
    groups: dict[str, str] = {}
    for leaf in leaves:
        row = selection.get(leaf)
        if row is not None:
            sp = species_label(row)
            labels[leaf] = f"{sp} {row.accession}" if sp != "unclassified" else row.accession
            groups[leaf] = sp if sp != "unclassified" else ""
    supports = [
        s
        for s in (
            parse_support(n.label) for n in tree.internal_nodes() if n.parent_node is not None
        )
        if s is not None
    ]
    scale100 = any(s > 1.0 for s in supports)
    hi, mid = (95.0, 70.0) if scale100 else (0.95, 0.70)
    lengths = [e.length for e in tree.postorder_edge_iter() if e.length is not None]
    total = sum(lengths) if lengths else None
    root_to_tip = [leaf.distance_from_root() for leaf in tree.leaf_node_iter()] if lengths else []
    rtt = describe(root_to_tip)
    model = _model(wd)

    report = StepReport(step="phylogeny", title="Phylogeny")
    report.key_numbers = [
        KeyNumber("Leaves", len(leaves)),
        KeyNumber("Tree length", round(total, 4) if total is not None else None, _unit(builder)),
    ]
    if supports:
        report.key_numbers.append(
            KeyNumber(
                f"Nodes with support of at least {hi:g}",
                f"{100 * sum(s >= hi for s in supports) / len(supports):.0f}%",
            )
        )
    # Without a stage record (a Nextflow run) the builder is not known.
    text = f"The tree holds {len(leaves):,} genomes"
    if builder:
        text = f"The tree of {len(leaves):,} genomes was built with {_display(builder)}"
    root = params.get("root")
    if root == "outgroup" and outgroup:
        text += f" and rooted on the outgroup {outgroup}"
    elif root == "midpoint":
        text += " and rooted at the midpoint"
    text += "."
    if supports:
        text += (
            f" Of {len(supports):,} internal branches, {sum(s >= hi for s in supports):,} "
            f"have a support of at least {hi:g} and {sum(s < mid for s in supports):,} "
            f"below {mid:g}."
        )
    report.summary.append(text)

    rows = [
        ("leaves", len(leaves)),
        ("tree builder", _display(builder)),
        ("substitution model", model),
        ("rooting", root or ""),
        ("outgroup", outgroup or ""),
        ("bootstrap replicates", params.get("bootstrap") or ""),
        ("internal branches with support", len(supports)),
        ("median support", _r(describe(supports)["median"])),
        (f"branches with support >= {hi:g} (%)", _pct(supports, hi)),
        (f"branches with support < {mid:g} (%)", _pct([-s for s in supports], -mid, strict=True)),
        (f"tree length ({_unit(builder)})", _r(total, 5)),
        ("median root-to-tip distance", _r(rtt["median"], 5)),
        ("maximum root-to-tip distance", _r(rtt["max"], 5)),
    ]
    report.tables.append(
        write_table(
            out,
            "tree_summary",
            "Summary of the tree.",
            ["quantity", "value"],
            [(k, v) for k, v in rows if v not in ("", None)],
        )
    )
    if len(leaves) > MAX_LEAVES:
        report.notes.append(
            f"The tree has {len(leaves):,} leaves; no tree figure is drawn above {MAX_LEAVES:,}."
        )
    else:
        report.figures.append(
            _tree_figure(out, tree, labels, groups, outgroup, builder, inp.formats)
        )
    if len(supports) >= 3:
        report.figures.append(_support_figure(out, supports, hi, mid, inp.formats))
    if len(root_to_tip) >= 3:
        report.figures.append(_rtt_figure(out, root_to_tip, builder, inp.formats))
    return report


def _r(value, digits: int = 2):
    return None if value is None else round(float(value), digits)


def _pct(values: list[float], threshold: float, *, strict: bool = False) -> float | None:
    if not values:
        return None
    n = sum((v > threshold) if strict else (v >= threshold) for v in values)
    return round(100 * n / len(values), 1)


def _unit(builder: str) -> str:
    return UNITS.get(builder, "branch length units")


def _display(builder: str) -> str:
    from ..methods import display

    return display(builder) or "an unrecorded tree builder"


def _tree_figure(out, tree, labels, groups, outgroup, builder, formats) -> FigureRef:
    n = sum(1 for _ in tree.leaf_node_iter())
    drawn: dict = {}

    def draw(fig) -> None:
        drawn.update(
            draw_tree(
                fig, tree, labels=labels, groups=groups, outgroup=outgroup, unit=_unit(builder)
            )
        )

    files = render(
        out, "tree", draw, width=style.DOUBLE_COLUMN, height=figure_height(n), formats=formats
    )
    parts = [
        f"{'Phylogram' if drawn.get('phylogram') else 'Cladogram'} of the {n:,} genomes"
        f" ({_display(builder)}), ladderised."
    ]
    if outgroup and drawn.get("labelled"):
        parts.append("The outgroup is shown in bold.")
    if drawn.get("support_mode") == "numbers":
        parts.append("Numbers at the nodes are support values.")
    elif drawn.get("support_mode") == "dots":
        hi, mid = (95, 70) if drawn.get("support_scale") == 100 else (0.95, 0.7)
        parts.append(
            f"Black circles mark nodes with support of at least {hi:g}, grey circles "
            f"{mid:g} to {hi:g}."
        )
    if drawn.get("groups_shown"):
        parts.append(
            f"Tip colours mark the {drawn['groups_shown']} species with most genomes"
            + (" (others grey)." if drawn["groups_total"] > drawn["groups_shown"] else ".")
        )
    if not drawn.get("labelled"):
        parts.append(f"Leaf labels are left out above {MAX_LABELLED} leaves.")
    if drawn.get("scale"):
        parts.append(f"Scale bar: {drawn['scale']:g} {_unit(builder)}.")
    return FigureRef(
        id="tree",
        caption=" ".join(parts),
        alt="Phylogenetic tree of the genomes.",
        files=files,
    )


def _support_figure(out, supports, hi, mid, formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        top = 100.0 if hi > 1 else 1.0
        ax.hist(supports, bins=20, range=(0, top), color=style.PRIMARY, linewidth=0)
        for value in (mid, hi):
            ax.axvline(value, color="black", lw=0.6, ls="--")
        ax.set_xlabel("Support")
        ax.set_ylabel("Internal branches")
        integer_ticks(ax.yaxis)

    files = render(out, "support", draw, formats=formats)
    return FigureRef(
        id="support",
        caption=f"Distribution of the support values of the internal branches; dashed "
        f"lines at {mid:g} and {hi:g}.",
        alt="Histogram of node support values.",
        files=files,
    )


def _rtt_figure(out, values, builder, formats) -> FigureRef:
    def draw(fig) -> None:
        ax = fig.add_subplot()
        ax.hist(values, bins=30, color=style.PRIMARY, linewidth=0)
        ax.set_xlabel(f"Root-to-tip distance ({_unit(builder)})")
        ax.set_ylabel("Genomes")
        integer_ticks(ax.yaxis)

    files = render(out, "root_to_tip", draw, formats=formats)
    return FigureRef(
        id="root_to_tip",
        caption="Distance from the root to each leaf of the rooted tree.",
        alt="Histogram of root-to-tip distances.",
        files=files,
    )
