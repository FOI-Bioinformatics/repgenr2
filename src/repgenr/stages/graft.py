"""graft: join the hierarchies of working directories at the node they share.

A group too divergent for one alignment is built in tiers. A parent working
directory covers the whole set; in its tree one clade holds a subset as a
group. A child working directory holds that subset alone, aligned and treed
on its own, and so resolves it further. ``graft`` replaces the parent's clade
with the child's tree and writes one ``tree2tax.tsv`` and ``genomes_map.tsv``.

The seam, where the two meet, is the parent node whose genomes are exactly the
child's. Two cases are refused rather than resolved, because each is a finding
about the data: the child's genomes spread over more than one branch of the
parent, so the two trees disagree about the group; and genomes the child holds
that the parent does not. A stateless step: it reads the two working
directories' tree2tax output and writes to ``out_dir``.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from ..core.contracts import (
    GENOMES_MAP_TSV,
    TREE2TAX_TSV,
    read_genomes_map,
    read_tree2tax,
    write_genomes_map,
    write_tree2tax,
)
from ..core.errors import UserInputError, WorkdirError

OUTGROUP_ACCESSION_TXT = "outgroup_accession.txt"

Edges = list[tuple[str, str]]
Mapping = list[tuple[str, str]]


@dataclass
class GraftStepParams:
    parent: Path
    children: list[Path] = field(default_factory=list)
    out_dir: Path = Path("graft")


@dataclass
class _Tier:
    where: Path
    edges: Edges
    mapping: Mapping
    outgroup: str | None


def _read_tier(workdir: Path) -> _Tier:
    workdir = Path(workdir)
    paths = [workdir / TREE2TAX_TSV, workdir / GENOMES_MAP_TSV]
    missing = [p.name for p in paths if not p.is_file()]
    if missing:
        raise WorkdirError(
            f"{workdir} has no {' or '.join(missing)}; "
            "run tree2tax in that working directory first."
        )
    accession_file = workdir / OUTGROUP_ACCESSION_TXT
    outgroup = (
        accession_file.read_text(encoding="utf-8").strip() if accession_file.is_file() else ""
    )
    return _Tier(workdir, read_tree2tax(paths[0]), read_genomes_map(paths[1]), outgroup or None)


def _children(edges: Edges) -> dict[str, list[str]]:
    below: dict[str, list[str]] = defaultdict(list)
    for child, parent in edges:
        below[parent].append(child)
    return below


def _root(edges: Edges, where: Path) -> str:
    children = {child for child, _ in edges}
    roots = sorted({parent for _, parent in edges} - children)
    if len(roots) != 1:
        raise WorkdirError(f"{where}: the hierarchy has {len(roots)} roots, expected one.")
    return roots[0]


def _descendants(edges: Edges, node: str) -> set[str]:
    below = _children(edges)
    found: set[str] = set()
    stack = list(below.get(node, []))
    while stack:
        current = stack.pop()
        if current not in found:
            found.add(current)
            stack.extend(below.get(current, []))
    return found


def _genome_sets(edges: Edges, mapping: Mapping, root: str) -> dict[str, frozenset[str]]:
    """The accessions at or below every node."""
    at: dict[str, set[str]] = defaultdict(set)
    for accession, leaf in mapping:
        at[leaf].add(accession)
    below = _children(edges)
    sets: dict[str, frozenset[str]] = {}
    order: list[str] = []
    stack = [root]
    while stack:  # preorder, then fold it up in reverse
        node = stack.pop()
        order.append(node)
        stack.extend(below.get(node, []))
    for node in reversed(order):
        found = set(at.get(node, ()))
        for child in below.get(node, []):
            found |= sets[child]
        sets[node] = frozenset(found)
    return sets


def _depths(edges: Edges, root: str) -> dict[str, int]:
    below = _children(edges)
    depth = {root: 0}
    stack = [root]
    while stack:
        node = stack.pop()
        for child in below.get(node, []):
            depth[child] = depth[node] + 1
            stack.append(child)
    return depth


def _without_outgroup(child: _Tier) -> tuple[Edges, Mapping]:
    """The child's tree without its own outgroup, which tree2tax keeps as a leaf.

    The child is rooted on it, so once it is gone the root may hold a single
    node with exactly the root's genomes; that node is folded into the root.
    """
    edges, mapping = list(child.edges), list(child.mapping)
    if child.outgroup is not None:
        leaves = {leaf for accession, leaf in mapping if accession == child.outgroup}
        mapping = [(a, leaf) for a, leaf in mapping if a != child.outgroup]
        still_held = {leaf for _, leaf in mapping}
        gone = {leaf for leaf in leaves if leaf not in still_held}
        edges = [(c, p) for c, p in edges if c not in gone]
    root = _root(edges, child.where)
    while True:
        below = _children(edges)
        only = below.get(root, [])
        if len(only) != 1 or not below.get(only[0]):
            break
        single = only[0]
        edges = [(c, root if p == single else p) for c, p in edges if c != single]
    return edges, mapping


def _graft_one(
    edges: Edges, mapping: Mapping, parent_name: Path, child: _Tier, logger: logging.Logger
) -> tuple[Edges, Mapping]:
    child_edges, child_mapping = _without_outgroup(child)
    child_root = _root(child_edges, child.where)
    wanted = frozenset(accession for accession, _ in child_mapping)

    held = {accession for accession, _ in mapping}
    absent = sorted(wanted - held)
    if absent:
        raise UserInputError(
            f"{child.where} holds {len(absent)} genomes that the parent {parent_name} does not, "
            f"so no node of the parent can be the seam: {', '.join(absent[:5])}. Both working "
            "directories must be built from one genome set."
        )

    root = _root(edges, parent_name)
    sets = _genome_sets(edges, mapping, root)
    seams = [node for node, found in sets.items() if found == wanted]
    if not seams:
        holding = [node for node, found in sets.items() if wanted <= found]
        depth = _depths(edges, root)
        smallest = max(holding, key=lambda node: (depth.get(node, 0), node))
        extra = sorted(sets[smallest] - wanted)
        branches = sorted(c for c in _children(edges).get(smallest, []) if sets[c] & wanted)
        raise UserInputError(
            f"no node of the parent {parent_name} holds exactly the {len(wanted)} genomes of "
            f"{child.where}: the smallest that holds them all, {smallest}, also holds "
            f"{len(extra)} others ({', '.join(extra[:5])}), and their genomes are spread over its "
            f"branches {', '.join(branches)}. The two trees disagree about the group."
        )
    depth = _depths(edges, root)
    seam = max(seams, key=lambda node: depth.get(node, 0))

    replaced = _descendants(edges, seam)
    kept = [(c, p) for c, p in edges if c not in replaced]
    kept_nodes = {c for c, _ in kept} | {root}
    grafted = [(c, seam if p == child_root else p) for c, p in child_edges]
    clash = sorted({c for c, _ in grafted} & kept_nodes)
    if clash:
        raise UserInputError(
            f"{child.where} names {len(clash)} nodes that the parent {parent_name} also holds "
            f"outside the seam {seam}: {', '.join(clash[:5])}. Give each working directory its "
            "own tree2tax --node-basename, or name the nodes by --naming leafset."
        )
    logger.info(
        "Grafted %s at %s: %d parent nodes below it replaced by %d from the child.",
        child.where,
        seam,
        len(replaced),
        len({c for c, _ in grafted}),
    )
    merged_mapping = [(a, leaf) for a, leaf in mapping if a not in wanted] + child_mapping
    return kept + grafted, merged_mapping


def _preorder(edges: Edges) -> Edges:
    """Edges with every parent before its children, as tree2tax writes them."""
    root = _root(edges, Path("graft"))
    below = _children(edges)
    ordered: Edges = []
    stack = [root]
    while stack:
        node = stack.pop()
        for child in reversed(below.get(node, [])):
            ordered.append((child, node))
            stack.append(child)
    return ordered


def graft(params: GraftStepParams, logger: logging.Logger) -> tuple[Path, Path]:
    """Graft each child onto the parent in turn; return the written paths."""
    if not params.children:
        raise UserInputError("graft needs at least one --child.")
    parent = _read_tier(params.parent)
    edges, mapping = parent.edges, parent.mapping
    for child_dir in params.children:
        edges, mapping = _graft_one(edges, mapping, params.parent, _read_tier(child_dir), logger)
    out_dir = Path(params.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_path, map_path = out_dir / TREE2TAX_TSV, out_dir / GENOMES_MAP_TSV
    write_tree2tax(tree_path, _preorder(edges))
    write_genomes_map(map_path, mapping)
    logger.info("Wrote %s and %s", tree_path, map_path)
    return tree_path, map_path
