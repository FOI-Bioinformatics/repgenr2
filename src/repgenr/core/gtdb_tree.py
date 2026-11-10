"""The GTDB reference tree: where a target taxon sits, and which genomes lie
nearest outside it.

GTDB publishes one tree per release and domain (``bac120_r232.tree``,
``ar53_r232.tree``), inferred from concatenated marker proteins, with one leaf
per species representative (``GB_``/``RS_`` accessions) and the taxa as
internal node labels (``'100:g__Brucella'``). Its ranks are drawn on this
tree, so the closest genus to a genus (or the closest family to a family) is
the group that is sister to the target in it.

:func:`sister_rings` walks up from the most recent common ancestor (MRCA) of
the target's leaves and yields, one ancestor at a time, the leaves of the
subtrees that join the target there: the first ring is the sister group, the
second the group sister to target plus sister, and so on. Each leaf carries
its path length to the MRCA (substitutions per site on the marker
alignment).
"""

from __future__ import annotations

import gzip
import logging
from collections.abc import Collection, Iterator
from dataclasses import dataclass
from pathlib import Path

from . import http
from .errors import UserInputError, WorkdirError

GTDB_BASE = "https://data.gtdb.ecogenomic.org/releases"


def strip_prefix(name: str) -> str:
    """A GTDB leaf name without its ``GB_``/``RS_`` prefix."""
    return name[3:] if name[:3] in ("GB_", "RS_") else name


@dataclass
class ReferenceTree:
    """A rooted tree as parent and branch-length arrays; node 0 is the root."""

    parent: list[int]
    length: list[float]
    children: list[list[int]]
    leaves: dict[str, int]  # accession without prefix -> node

    def depths(self) -> list[float]:
        """Path length of every node from the root (parents precede children)."""
        depth = [0.0] * len(self.parent)
        for node in range(1, len(self.parent)):
            depth[node] = depth[self.parent[node]] + self.length[node]
        return depth

    def leaves_under(self, node: int) -> list[int]:
        out: list[int] = []
        stack = [node]
        while stack:
            n = stack.pop()
            if self.children[n]:
                stack.extend(self.children[n])
            else:
                out.append(n)
        return out


def parse_newick(text: str) -> ReferenceTree:
    """Parse a Newick tree without recursion (the GTDB trees are deep).

    Leaf names and branch lengths are kept; internal labels (support and
    taxa, quoted or not) are read and dropped. Nodes are numbered in the
    order their opening is read, so a parent always precedes its children.
    """
    parent: list[int] = [-1]
    length: list[float] = [0.0]
    children: list[list[int]] = [[]]
    names: dict[int, str] = {}
    current = 0
    i, n = 0, len(text)

    def new_node(par: int) -> int:
        parent.append(par)
        length.append(0.0)
        children.append([])
        children[par].append(len(parent) - 1)
        return len(parent) - 1

    def read_label(pos: int) -> tuple[str, int]:
        if pos < n and text[pos] == "'":
            end = text.index("'", pos + 1)
            return text[pos + 1 : end], end + 1
        start = pos
        while pos < n and text[pos] not in ",():;":
            pos += 1
        return text[start:pos].strip(), pos

    def read_length(pos: int, node: int) -> int:
        if pos < n and text[pos] == ":":
            start = pos = pos + 1
            while pos < n and text[pos] not in ",();":
                pos += 1
            try:
                length[node] = float(text[start:pos])
            except ValueError as exc:
                raise WorkdirError(f"Bad branch length {text[start:pos]!r} in tree") from exc
        return pos

    # The outermost parentheses are node 0; "(" descends into a first child,
    # "," starts a sibling, ")" returns to the parent and reads its label.
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c == "(":
            current = new_node(current)
            i += 1
        elif c == ",":
            current = new_node(parent[current])
            i += 1
        elif c == ")":
            current = parent[current]
            _label, i = read_label(i + 1)
            i = read_length(i, current)
        elif c == ";":
            break
        else:
            name, i = read_label(i)
            if name:
                names[current] = name
            i = read_length(i, current)
    leaves = {strip_prefix(nm): node for node, nm in names.items() if not children[node]}
    return ReferenceTree(parent, length, children, leaves)


def read_tree(path: Path) -> ReferenceTree:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as fo:  # type: ignore[operator]
        return parse_newick(fo.read())


@dataclass(frozen=True)
class Placement:
    """The target's MRCA, its leaves found in the tree, and the leaves under
    the MRCA that are not part of the target (when the target is not a single
    clade in the tree)."""

    mrca: int
    found: int
    intruders: int


def place(tree: ReferenceTree, target: Collection[str]) -> Placement | None:
    """The MRCA of the target accessions present in the tree; None when none is."""
    nodes = [tree.leaves[a] for a in target if a in tree.leaves]
    if not nodes:
        return None
    count = [0] * len(tree.parent)
    for node in nodes:
        count[node] += 1
    for node in range(len(tree.parent) - 1, 0, -1):
        count[tree.parent[node]] += count[node]
    total = len(nodes)
    mrca = 0
    while True:
        nxt = next((c for c in tree.children[mrca] if count[c] == total), None)
        if nxt is None:
            break
        mrca = nxt
    under = len(tree.leaves_under(mrca))
    return Placement(mrca=mrca, found=total, intruders=under - total)


def sister_rings(tree: ReferenceTree, placement: Placement) -> Iterator[list[tuple[str, float]]]:
    """Leaves outside the target, nearest group first, as (accession, path
    length to the target's MRCA)."""
    depth = tree.depths()
    name_of = {node: acc for acc, node in tree.leaves.items()}
    node = placement.mrca
    while tree.parent[node] >= 0:
        anc = tree.parent[node]
        ring: list[tuple[str, float]] = []
        for sib in tree.children[anc]:
            if sib == node:
                continue
            for leaf in tree.leaves_under(sib):
                if leaf in name_of:
                    dist = depth[leaf] + depth[placement.mrca] - 2 * depth[anc]
                    ring.append((name_of[leaf], dist))
        if ring:
            ring.sort(key=lambda x: (x[1], x[0]))
            yield ring
        node = anc


def _fetch(url: str, dest: Path, logger: logging.Logger) -> None:
    """Download ``url`` to ``dest`` (a seam the unit tests replace)."""
    http.download(url, dest, logger=logger)


def tree_filename(release: str, version: str) -> str:
    return f"{version}_r{int(float(release))}.tree"


def obtain_reference_tree(
    workdir: Path,
    release: str | None,
    version: str | None,
    *,
    tree_path: str | None,
    nodownload: bool,
    logger: logging.Logger,
) -> Path:
    """The GTDB reference tree: ``tree_path``, one kept in ``workdir``, or a
    download into ``workdir`` (checked against the release's MD5SUM.txt)."""
    if tree_path:
        path = Path(tree_path).expanduser()
        if not path.is_file():
            raise UserInputError(f"--gtdb-tree-path not found: {tree_path}")
        return path
    if not release or not version:
        raise UserInputError(
            "The GTDB reference tree needs -r/--release and --gtdb-version, or --gtdb-tree-path."
        )
    dest = workdir / tree_filename(release, version)
    if dest.is_file() and (nodownload or dest.stat().st_size > 0):
        logger.info("Using previously downloaded %s", dest.name)
        return dest
    major = int(float(release))
    url = f"{GTDB_BASE}/release{major}/{release}/{dest.name}"
    logger.info("Downloading %s", url)
    try:
        _fetch(url, dest, logger)
    except http.HTTPStatusError as exc:
        if exc.status == 404:
            raise UserInputError(
                f"GTDB has no reference tree {dest.name} for release {release}; "
                "give one with --gtdb-tree-path or use --outgroup-mode quality."
            ) from exc
        raise
    try:
        manifest = f"{GTDB_BASE}/release{major}/{release}/MD5SUM.txt"
        http.verify_md5_manifest(dest, manifest, logger=logger)
    except WorkdirError:
        dest.unlink(missing_ok=True)
        raise
    return dest
