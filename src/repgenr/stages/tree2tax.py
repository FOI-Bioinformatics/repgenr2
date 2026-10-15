"""tree2tax stage: turn a rooted tree into FlexTaxD-compatible relations.

Turns a rooted tree into a FlexTaxD child->parent table (via dendropy): read
``tree/tree.nwk``, root it on the outgroup unless phylo already did so (or
rooted it at the midpoint, which is kept), name internal nodes (a user basename
or a leaf-derived hash), then emit ``tree2tax.tsv`` (child -> parent) and
``genomes_map.tsv`` (accession -> leaf). Optionally lists redundant
(dereplicated) genomes under their representative leaf, read from the derep
``clusters.tsv`` contract.
"""

from __future__ import annotations

import hashlib
import logging
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import dendropy
from dendropy.utility.error import DataParseError

from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTERS_TSV,
    GENOMES_MAP_TSV,
    OUTGROUP_ACCESSION_TXT,
    SEGMENTS_TSV,
    TREE2TAX_NODES_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    NodeRow,
    accession_from_filename,
    list_fasta,
    newick_is_complete,
    read_clusters,
    read_segments,
    strip_fasta_suffix,
    write_genomes_map,
    write_node_table,
    write_tree2tax,
)
from ..core.errors import WorkdirError
from ..tree.rooting import ROOT_MIDPOINT, root_on_outgroup


@dataclass
class Tree2taxParams:
    node_basename: str | None = None
    root_name: str = "root"
    remove_outgroup: bool = False
    # Default on: the genomes_map deliverable lists redundant genomes under
    # their representative, matching what `repgenr run` has always produced.
    include_dereplicated: bool = True
    # Weak-split collapse thresholds (both off by default): a node merges into
    # its parent when its support is below collapse_support (a fraction; percent
    # trees are normalised) or its branch is shorter than collapse_length.
    collapse_support: float | None = None
    collapse_length: float | None = None
    # Clade diameter: the largest leaf-to-leaf distance a named clade may span;
    # the smallest set of such clades is kept and everything inside each is
    # merged into it (TreeCluster's max-diameter criterion).
    collapse_diameter: float | None = None


@dataclass
class Tree2taxStepParams:
    """Inputs for the stateless tree2tax step (explicit paths, no workdir)."""

    tree: Path
    out_dir: Path
    clusters: Path | None = None
    outgroup_dir: Path | None = None
    outgroup_accession: Path | None = None
    node_basename: str | None = None
    root_name: str = "root"
    remove_outgroup: bool = False
    # On by default, as on the workdir command and `run`: the genomes_map
    # deliverable lists redundant genomes under their representative.
    include_dereplicated: bool = True
    versions_out: Path | None = None
    collapse_support: float | None = None
    collapse_length: float | None = None
    collapse_diameter: float | None = None
    # The tree was built without an outgroup (phylo --no-outgroup): ignore the
    # staged outgroup inputs and leave the tree unrooted.
    no_outgroup: bool = False
    # The tree is already rooted as wanted (phylo --root midpoint): do not
    # reroot it on the outgroup. The outgroup is still read for --remove-outgroup.
    keep_root: bool = False
    # vgenome --group-segments writes segments.tsv; with it, every member
    # segment accession of a grouped isolate reaches genomes_map.tsv.
    segments: Path | None = None


# FlexTaxD pre-creates a node with this name and validates the tree from it.
FLEXTAXD_ROOT = "root"


def _emit_relations(
    tree_text: str,
    outgroup_leaf: str | None,
    redundant: dict[str, list[str]],
    *,
    node_basename: str | None,
    root_name: str,
    remove_outgroup: bool,
    out_tree2tax: Path,
    out_map: Path,
    logger: logging.Logger,
    collapse_support: float | None = None,
    collapse_length: float | None = None,
    collapse_diameter: float | None = None,
    segments: dict[str, list[str]] | None = None,
    tree_source: Path | None = None,
    keep_root: bool = False,
    out_nodes: Path | None = None,
    lineage: dict[str, str] | None = None,
) -> tuple[Path, Path, int]:
    """Build FlexTaxD relations from a tree and write the two output tables.

    Stateless core shared by :func:`run` (workdir-bound) and
    :func:`tree2tax_relations` (data-channel step): it roots, collapses weak
    splits when asked, names nodes and emits ``tree2tax.tsv`` +
    ``genomes_map.tsv`` to the given paths. Returns the two paths and the
    number of internal nodes collapsed.
    """
    source = tree_source if tree_source is not None else "the tree"
    if not newick_is_complete(tree_text):
        # Same rule as doctor: dendropy would read the first tree and ignore
        # whatever follows its ';', which hides a truncated or concatenated file.
        raise WorkdirError(
            f"{source} is empty or truncated: it has no terminating ';', has text "
            "after its final ';', or holds more than one tree. Re-run phylo."
        )
    # preserve_underscores: genome leaf names contain '_' (Family_Genus_species_Acc)
    # and newick otherwise turns underscores into spaces.
    try:
        tree = dendropy.Tree.get(data=tree_text, schema="newick", preserve_underscores=True)
    except (DataParseError, ValueError) as exc:
        raise WorkdirError(f"{source} is not a valid Newick tree: {exc}") from exc

    if outgroup_leaf is not None and not keep_root:
        if not _set_outgroup(tree, outgroup_leaf, source):
            logger.info("The tree is already rooted on the outgroup %s", outgroup_leaf)

    # The outgroup leaves the tree before collapse and naming, so no node is
    # named for a leaf set that includes it and no unifurcation is left where
    # it sat (a midpoint-rooted tree can hold it inside a clade).
    if remove_outgroup and outgroup_leaf is not None:
        _prune_leaf(tree, outgroup_leaf)

    # Collapse before naming so node names describe the collapsed topology:
    # clades under the diameter first, then the weak splits that remain.
    by_diameter = _collapse_by_diameter(tree, collapse_diameter, logger, source=source)
    stats = _collapse_weak_nodes(
        tree, support=collapse_support, length=collapse_length, logger=logger
    )
    stats.collapsed += by_diameter

    if root_name != FLEXTAXD_ROOT:
        logger.warning(
            "The top node is named %r; FlexTaxD starts its own tree at a node named %r, "
            "so a taxonomy with another top name does not validate there",
            root_name,
            FLEXTAXD_ROOT,
        )
    # Supports sit on the internal labels, which naming replaces; written as
    # fractions, as --collapse-support reads them.
    supports = {
        id(n): _parse_support(n.label) for n in tree.internal_nodes() if n is not tree.seed_node
    }
    if any(v is not None and v > 1.0 for v in supports.values()):
        supports = {k: (None if v is None else v / 100.0) for k, v in supports.items()}
    leaves = _name_nodes(tree, node_basename)
    tree.seed_node.label = root_name

    write_tree2tax(out_tree2tax, _edges(tree))
    write_genomes_map(out_map, _genome_map(leaves, redundant, segments))
    if out_nodes is not None:
        rooted = keep_root or outgroup_leaf is not None
        write_node_table(
            out_nodes, _node_rows(tree, supports, rooted=rooted, lineage=lineage or {})
        )
        logger.info("Wrote %s, %s and %s", out_tree2tax.name, out_map.name, out_nodes.name)
    else:
        logger.info("Wrote %s and %s", out_tree2tax.name, out_map.name)
    return out_tree2tax, out_map, stats.collapsed


def _node_rows(
    tree: dendropy.Tree,
    supports: dict[int, float | None],
    *,
    rooted: bool,
    lineage: dict[str, str],
) -> list[NodeRow]:
    """One row per node of the taxonomy (the root included), in preorder.

    Depth, branch length and relative evolutionary divergence need every
    branch length; without them, or without a root that means something (an
    outgroup or midpoint root), those columns are NA and the rest stands.
    ``lineage`` maps a genome accession to its family;genus;species.
    """
    lengths_known = all(
        n.edge.length is not None for n in tree.preorder_node_iter() if n is not tree.seed_node
    )
    measured = lengths_known and rooted
    n_leaves: dict[int, int] = {}
    mean_tip: dict[int, float] = {}  # mean distance from a node to its leaves
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            n_leaves[id(node)] = 1
            mean_tip[id(node)] = 0.0
            continue
        kids = node.child_nodes()
        n_leaves[id(node)] = sum(n_leaves[id(c)] for c in kids)
        if measured:
            total = sum((mean_tip[id(c)] + float(c.edge.length)) * n_leaves[id(c)] for c in kids)
            mean_tip[id(node)] = total / n_leaves[id(node)]
    depth: dict[int, float] = {id(tree.seed_node): 0.0}
    red: dict[int, float] = {id(tree.seed_node): 0.0}
    # Under --remove-outgroup the root keeps a single child, the ingroup's
    # most recent common ancestor; RED is measured from there, as GTDB does.
    origin = tree.seed_node
    if origin.num_child_nodes() == 1:
        origin = origin.child_nodes()[0]
        red[id(origin)] = 0.0
    rows: list[NodeRow] = []
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        name = _leaf_label(node) if node.is_leaf() else node.label
        if not name:
            continue
        d = b = r = None
        if parent is not None and measured:
            b = float(node.edge.length)
            d = depth[id(parent)] + b
            depth[id(node)] = d
            if node is origin:
                r = 0.0
            elif node.is_leaf():
                r = 1.0
            else:
                # Parks et al. 2018, as PhyloRank computes it: a is the branch
                # to the parent, b the mean distance from this node to its
                # own leaves; a zero-length clade sits at its parent's value.
                a, mean_below = max(b, 0.0), mean_tip[id(node)]
                parent_red = red[id(parent)]
                r = (
                    parent_red + (a / (a + mean_below)) * (1.0 - parent_red)
                    if a + mean_below > 0
                    else parent_red
                )
                r = min(max(r, parent_red), 1.0)
            red[id(node)] = r
        elif parent is None and measured:
            d, r = 0.0, 0.0
        rows.append(
            NodeRow(
                node=name,
                parent=None if parent is None else parent.label,
                n_leaves=n_leaves[id(node)],
                depth=d,
                branch_length=b,
                support=None if node.is_leaf() else supports.get(id(node)),
                red=r,
                lineage=lineage.get(_accession(name)) if node.is_leaf() else None,
            )
        )
    return rows


def tree2tax_relations(params: Tree2taxStepParams, logger: logging.Logger) -> tuple[Path, Path]:
    """Emit FlexTaxD relations from explicit inputs (stateless; no config)."""
    if not params.tree.exists():
        raise WorkdirError(f"Tree not found: {params.tree}. Run the phylo step first.")
    if params.include_dereplicated and params.clusters is not None and not params.clusters.exists():
        # An explicit path that is absent would otherwise drop every dereplicated
        # member from genomes_map.tsv without notice.
        raise WorkdirError(f"Clusters table not found: {params.clusters}.")
    if params.segments is not None and not params.segments.exists():
        raise WorkdirError(f"Segments table not found: {params.segments}.")
    segments = read_segments(params.segments) if params.segments is not None else None
    outgroup_leaf = None
    if params.no_outgroup and params.keep_root:
        logger.info("No outgroup (--no-outgroup); the tree's own root is kept (--keep-root)")
    elif params.no_outgroup:
        logger.warning("No outgroup (--no-outgroup); tree is left unrooted")
    elif params.outgroup_dir is not None and params.outgroup_accession is not None:
        outgroup_leaf = _resolve_outgroup_leaf_from(
            params.outgroup_dir, params.outgroup_accession, logger
        )
    redundant = (
        _load_redundant_from(params.clusters)
        if params.include_dereplicated and params.clusters is not None
        else {}
    )
    params.out_dir.mkdir(parents=True, exist_ok=True)
    if params.versions_out is not None:
        # No external binaries on this step; the tree work is the dendropy
        # library, recorded as the tool, as the workdir stage does.
        from ..core.versions import write_versions_fragment

        write_versions_fragment(params.versions_out, _dendropy_versions())
    out_tree2tax, out_map, _collapsed = _emit_relations(
        params.tree.read_text().strip(),
        outgroup_leaf,
        redundant,
        node_basename=params.node_basename,
        root_name=params.root_name,
        remove_outgroup=params.remove_outgroup,
        out_tree2tax=params.out_dir / TREE2TAX_TSV,
        out_map=params.out_dir / GENOMES_MAP_TSV,
        logger=logger,
        collapse_support=params.collapse_support,
        collapse_length=params.collapse_length,
        collapse_diameter=params.collapse_diameter,
        segments=segments,
        tree_source=params.tree,
        keep_root=params.keep_root,
        out_nodes=params.out_dir / TREE2TAX_NODES_TSV,
    )
    return out_tree2tax, out_map


def run(ctx: WorkdirContext, params: Tree2taxParams) -> tuple[Path, Path]:
    logger = ctx.logger
    tree_file = ctx.tree_dir / TREE_NWK
    if not tree_file.exists():
        raise WorkdirError(f"Tree not found: {tree_file}. Run the phylo stage first.")
    phylo = ctx.config.stages.get("phylo")
    if phylo is not None and not phylo.completed:
        # tree.nwk is replaced only when a build succeeds, so the file is a
        # whole tree, but from an earlier phylo run than the one recorded.
        logger.warning(
            "The last phylo run did not finish; %s is the tree of an earlier run. "
            "Re-run phylo before tree2tax to use the current settings.",
            tree_file,
        )

    keep_root = phylo is not None and phylo.params.get("root") == ROOT_MIDPOINT
    if keep_root:
        logger.info("phylo rooted the tree at the midpoint; that root is kept")
    outgroup_leaf = _resolve_outgroup_leaf(ctx, logger, quiet=keep_root)
    redundant = _load_redundant(ctx) if params.include_dereplicated else {}
    segments_path = ctx.workdir / SEGMENTS_TSV
    segments = read_segments(segments_path) if segments_path.exists() else None

    out_tree2tax, out_map, collapsed = _emit_relations(
        tree_file.read_text().strip(),
        outgroup_leaf,
        redundant,
        node_basename=params.node_basename,
        root_name=params.root_name,
        remove_outgroup=params.remove_outgroup,
        out_tree2tax=ctx.workdir / TREE2TAX_TSV,
        out_map=ctx.workdir / GENOMES_MAP_TSV,
        logger=logger,
        collapse_support=params.collapse_support,
        collapse_length=params.collapse_length,
        collapse_diameter=params.collapse_diameter,
        segments=segments,
        tree_source=tree_file,
        keep_root=keep_root,
        out_nodes=ctx.workdir / TREE2TAX_NODES_TSV,
        lineage=_leaf_lineages(ctx),
    )

    ctx.config.record_stage(
        "tree2tax",
        tool="dendropy",
        tool_versions=_dendropy_versions(),
        params={
            "remove_outgroup": params.remove_outgroup,
            "include_dereplicated": params.include_dereplicated,
            "root_name": params.root_name,
            "node_basename": params.node_basename,
            "collapse_support": params.collapse_support,
            "collapse_length": params.collapse_length,
            "collapse_diameter": params.collapse_diameter,
            "collapsed_nodes": collapsed,
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return out_tree2tax, out_map


def _dendropy_versions() -> dict[str, str]:
    """The dendropy library version, recorded as this stage's tool version."""
    try:
        return {"dendropy": version("dendropy")}
    except PackageNotFoundError:  # pragma: no cover - dendropy is a hard dependency
        return {}


def _resolve_outgroup_leaf(ctx: WorkdirContext, logger, *, quiet: bool = False) -> str | None:
    """The outgroup leaf to root on, or None.

    ``quiet`` (a tree whose root is kept) drops the 'left unrooted' warnings:
    the outgroup is then needed only for --remove-outgroup.
    """
    phylo = ctx.config.stages.get("phylo")
    if phylo is not None and phylo.completed and phylo.params.get("outgroup", "") is None:
        # phylo ran without an outgroup (--no-outgroup, or none was found), so
        # the tree has no outgroup leaf to root on.
        if not quiet:
            logger.warning("phylo built the tree without an outgroup; tree is left unrooted")
        return None
    acc_file = ctx.workdir / OUTGROUP_ACCESSION_TXT
    if not acc_file.exists() or not ctx.outgroup_dir.exists():
        if not quiet:
            logger.warning("No outgroup available; tree is left unrooted")
        return None
    return _resolve_outgroup_leaf_from(ctx.outgroup_dir, acc_file, logger)


def _resolve_outgroup_leaf_from(
    outgroup_dir: Path, accession_file: Path, logger: logging.Logger
) -> str | None:
    if not accession_file.exists() or not outgroup_dir.exists():
        logger.warning("No outgroup available; tree is left unrooted")
        return None
    accession = accession_file.read_text(encoding="utf-8").strip()
    if not accession:
        logger.warning("No outgroup accession recorded; tree is left unrooted")
        return None
    # The same candidates as phylo's resolve_outgroup_files, so the leaf
    # named here is the one phylo placed in the tree.
    candidates = list_fasta(outgroup_dir)
    for f in candidates:
        if accession_from_filename(f.name) == accession:
            return strip_fasta_suffix(f.name)
    for f in candidates:
        if accession in f.name:
            logger.info("Outgroup resolved by substring match: %s", f.name)
            return strip_fasta_suffix(f.name)
    logger.warning(
        "No file in %s matches outgroup accession %s; tree is left unrooted",
        outgroup_dir,
        accession,
    )
    return None


def _leaf_label(node) -> str:
    return node.taxon.label if node.taxon is not None else ""


def _set_outgroup(tree: dendropy.Tree, leaf_label: str, source: object) -> bool:
    """Root on the outgroup's branch, as phylo does; False when already rooted so.

    Leaving the tree unrooted when the outgroup is missing would give a
    taxonomy rooted at an arbitrary node with exit 0, so that is an error.
    """
    return root_on_outgroup(
        tree,
        leaf_label,
        source=source,
        hint=(
            " Rebuild the tree with this outgroup; if the tree was built without one "
            "(phylo --no-outgroup), pass --no-outgroup to tree2tax-relations."
        ),
    )


@dataclass
class _CollapseStats:
    collapsed: int = 0
    by_support: int = 0
    by_length: int = 0


def _parse_support(label: str | None) -> float | None:
    if label is None:
        return None
    try:
        return float(label)
    except ValueError:
        return None


def _collapse_by_diameter(
    tree: dendropy.Tree,
    diameter: float | None,
    logger: logging.Logger,
    *,
    source: object = "the tree",
) -> int:
    """Merge every clade whose leaf-to-leaf distance stays within ``diameter``
    into one node; return the number of internal nodes removed.

    One post-order pass keeps, per node, the two deepest leaf depths through
    different children; their sum is the clade's diameter. The highest nodes
    under the limit are kept (TreeCluster's max-diameter criterion gives the
    fewest clades), and every internal node strictly below them is collapsed,
    so each clade is one node with its leaves as children. The root and its
    children are never collapsed into, but a root child can be such a clade.
    A tree without branch lengths cannot be measured and is refused.
    """
    if diameter is None:
        return 0
    deepest: dict[int, float] = {}
    span: dict[int, float] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            deepest[id(node)] = 0.0
            span[id(node)] = 0.0
            continue
        depths = []
        for child in node.child_nodes():
            if child.edge.length is None:
                raise WorkdirError(
                    f"{source} has a branch without a length, so --collapse-diameter cannot "
                    "measure clades; build the tree with a builder that writes lengths."
                )
            depths.append(deepest[id(child)] + float(child.edge.length))
        depths.sort(reverse=True)
        deepest[id(node)] = depths[0]
        span[id(node)] = depths[0] + (depths[1] if len(depths) > 1 else 0.0)
    protected = {tree.seed_node}
    removed = 0
    stack = list(tree.seed_node.child_nodes())
    kept = 0
    while stack:
        node = stack.pop()
        if node.is_leaf() or node in protected:
            continue
        if span[id(node)] <= diameter * (1 + 1e-9) + 1e-12:  # float sums of branch lengths
            kept += 1
            inner = [n for n in node.postorder_internal_node_iter() if n is not node]
            for n in inner:
                n.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
            removed += len(inner)
            continue
        stack.extend(node.child_nodes())
    if kept:
        logger.info(
            "Collapsed %d internal node(s) inside %d clade(s) of diameter at most %g",
            removed,
            kept,
            diameter,
        )
    return removed


def _collapse_weak_nodes(
    tree: dendropy.Tree,
    *,
    support: float | None,
    length: float | None,
    logger: logging.Logger,
) -> _CollapseStats:
    """Merge weakly supported or near-zero-length internal nodes into their parents.

    A node collapses when its support is below ``support`` (a fraction; when
    any label in the tree exceeds 1 the tree is read as percentages) or its
    branch is shorter than ``length``; either criterion suffices. The root, the
    root's children (the outgroup/ingroup split) and leaves never collapse.
    Decisions are taken on the tree as read, then applied, so one collapse
    cannot change another node's verdict; a collapsed node's branch length is
    added to its children's. Thresholds of None disable a criterion.
    """
    stats = _CollapseStats()
    if support is None and length is None:
        return stats

    protected = {tree.seed_node, *tree.seed_node.child_nodes()}
    candidates = [n for n in tree.internal_nodes() if n not in protected]

    supports: dict[dendropy.Node, float] = {}
    if support is not None:
        for node in tree.internal_nodes():
            if node is tree.seed_node:
                continue
            value = _parse_support(node.label)
            if value is not None:
                supports[node] = value
        if not supports:
            logger.warning(
                "--collapse-support %.3g given but the tree carries no support values; "
                "nothing collapsed on support (use --collapse-length, or build the tree "
                "with --bootstrap).",
                support,
            )
        elif max(supports.values()) > 1.0:
            logger.info("Support values read as percentages; comparing against %.1f", support * 100)
            supports = {n: v / 100.0 for n, v in supports.items()}

    to_collapse: list[dendropy.Node] = []
    for node in candidates:
        weak_support = support is not None and node in supports and supports[node] < support
        short = length is not None and node.edge.length is not None and node.edge.length < length
        if weak_support:
            stats.by_support += 1
        if short:
            stats.by_length += 1
        if weak_support or short:
            to_collapse.append(node)

    for node in to_collapse:
        node.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
    stats.collapsed = len(to_collapse)
    if to_collapse:
        logger.info(
            "Collapsed %d of %d internal nodes (%d below support %s, %d shorter than %s)",
            stats.collapsed,
            len(candidates),
            stats.by_support,
            "n/a" if support is None else f"{support:g}",
            stats.by_length,
            "n/a" if length is None else f"{length:g}",
        )
    return stats


_HASH_BITS = 128
_HASH_MOD = 1 << _HASH_BITS


def _leaf_hash(label: str) -> int:
    return int.from_bytes(hashlib.md5(label.encode("utf-8")).digest(), "big")


def _prune_leaf(tree: dendropy.Tree, label: str) -> None:
    """Remove the leaf called ``label`` and any unifurcation this leaves below
    the root. The root keeps a single child when the leaf was its other one
    (an outgroup-rooted tree), so the ingroup stays one clade under the root."""
    node = next((n for n in tree.leaf_node_iter() if _leaf_label(n) == label), None)
    if node is None:
        return
    parent = node.parent_node
    tree.prune_subtree(node, suppress_unifurcations=False)
    while parent is not None and parent is not tree.seed_node and parent.num_child_nodes() == 1:
        grandparent = parent.parent_node
        parent.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
        parent = grandparent


def _name_nodes(tree: dendropy.Tree, node_basename: str | None) -> list[str]:
    """Give every internal (non-root) node a stable name; return the leaf names.

    An internal node is named by a hash of its descendant leaf set: the sum
    of the leaves' 128-bit label hashes, which is independent of the order
    the leaves are met in and takes one post-order pass over the tree, so a
    long ladder costs the same per node as a balanced tree. The name does not
    depend on how the clade is resolved inside. A node whose leaf set equals
    one child's (a single child below the root, or one labelled child beside
    unlabelled leaves; left by a tool, never by this stage) would share that child's
    name, so its label is mixed once more while the plain sum still passes
    upward, so the ancestors' names stay those of their leaf sets. With a
    basename, nodes are numbered in preorder instead. Keyed on ``is_leaf``,
    so an internal support label is never taken for a taxon.
    """
    leaves: list[str] = []
    sums: dict[int, int] = {}  # plain leaf-set sum, passed upward
    names: dict[int, int] = {}  # the value a node is labelled with
    with_leaf_set: set[int] = set()  # internal nodes that hold a labelled leaf
    counter = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            name = _leaf_label(node)
            if name:
                leaves.append(name)
                sums[id(node)] = names[id(node)] = _leaf_hash(name)
            continue
        contributing = [child for child in node.child_nodes() if id(child) in sums]
        if not contributing:
            continue  # only unlabelled leaves below: no leaf set, no name, no row
        sums[id(node)] = sum(sums.pop(id(child)) for child in contributing) % _HASH_MOD
        with_leaf_set.add(id(node))
        if len(contributing) >= 2:
            names[id(node)] = sums[id(node)]
        else:
            # The same leaf set as the one child: mix that child's own name,
            # so a chain of such nodes gets a distinct name at every level.
            names[id(node)] = _leaf_hash(f"{names[id(contributing[0])]:032x}/same-leaf-set")
        for child in contributing:
            names.pop(id(child), None)
        if node is not tree.seed_node and not node_basename:
            node.label = f"{names[id(node)]:032x}"
    if node_basename:
        for node in tree.preorder_node_iter():
            if node is not tree.seed_node and not node.is_leaf() and id(node) in with_leaf_set:
                node.label = f"{node_basename}{counter}"
                counter += 1
    return leaves


def _edges(tree: dendropy.Tree) -> list[tuple[str, str]]:
    """One (child, parent) row per non-root node, in preorder. A leaf without
    a label has no row."""
    edges: list[tuple[str, str]] = []
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        child = _leaf_label(node) if node.is_leaf() else node.label
        parent = node.parent_node.label
        if child and parent:
            edges.append((child, parent))
    return edges


def _leaf_lineages(ctx: WorkdirContext) -> dict[str, str]:
    """Accession -> family;genus;species from the manifest, for the node table."""
    out: dict[str, str] = {}
    try:
        records = ctx.manifest.all_genomes(include_outgroup=True)
    except Exception as exc:  # noqa: BLE001 - the table is an annotation, never a failure
        ctx.logger.warning("Lineages left out of the node table: %s", exc)
        return out
    for rec in records:
        if rec.family or rec.genus or rec.species:
            out[rec.accession] = ";".join(x or "" for x in (rec.family, rec.genus, rec.species))
    return out


def _load_redundant(ctx: WorkdirContext) -> dict[str, list[str]]:
    return _load_redundant_from(ctx.derep_dir / CLUSTERS_TSV)


def _load_redundant_from(clusters_file: Path) -> dict[str, list[str]]:
    if not clusters_file.exists():
        return {}
    clusters = read_clusters(clusters_file)
    # map representative leaf-stem -> redundant leaf-stems
    out: dict[str, list[str]] = {}
    for rep, members in clusters.items():
        rep_stem = strip_fasta_suffix(rep)
        redundant = [strip_fasta_suffix(m) for m in members]
        if redundant:
            out[rep_stem] = redundant
    return out


def _genome_map(
    leaves: Sequence[str],
    redundant: dict[str, list[str]],
    segments: dict[str, list[str]] | None = None,
) -> list[tuple[str, str]]:
    """accession -> leaf rows: each leaf, its redundant members, and for a
    segment-grouped isolate (whose accession is a synthetic token) the member
    segment accessions as well, so every real accession reaches the map."""
    members = segments or {}
    mapping: list[tuple[str, str]] = []

    def _add(genome: str, leaf: str) -> None:
        accession = _accession(genome)
        mapping.append((accession, leaf))
        for segment in members.get(accession, []):
            mapping.append((segment, leaf))

    leaf_set = set(leaves)
    for leaf in leaves:
        _add(leaf, leaf)
        for red in redundant.get(leaf, []):
            # A member that is a leaf itself (a tree built with --all-genomes)
            # already maps to its own leaf.
            if red not in leaf_set:
                _add(red, leaf)
    return mapping


def _accession(leaf: str) -> str:
    return accession_from_filename(leaf) or leaf
