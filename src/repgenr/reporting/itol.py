"""Annotation files for iTOL (https://itol.embl.de) for the tree of a run.

``repgenr itol -wd WD`` writes, into ``tree/itol/``:

``labels.txt`` (LABELS)
    Each leaf as species, strain, accession and the number of genomes of its
    cluster, e.g. ``B. abortus 2308 (GCF_000054005.1, n=12)``.
``colorstrip_species.txt`` (DATASET_COLORSTRIP)
    The species of each leaf (the representative), one colour per species.
``bar_cluster_size.txt`` (DATASET_SIMPLEBAR)
    The genomes of each leaf's cluster, representative included.
``pie_cluster_species.txt`` (DATASET_PIECHART)
    The species of all genomes of each leaf's cluster, the pie area growing
    with the cluster size, in the colours of the colour strip.
``clusters_species.tsv``
    The same per leaf as a table.

Species names are NCBI's (``--names ncbi``, the default: the species of the
NCBI taxonomy, else that of the unfiltered NCBI taxonomy, which names
genomes without a species as ``<Genus> sp. <strain>``) or GTDB's (``--names
gtdb``). They come from ``metadata_selected.tsv`` and, for a genome not
listed there, from the manifest or the genome's file name. Clusters come
from ``derep/clusters.tsv``; a leaf that is no representative (the outgroup,
or a tree of all genomes) counts as a cluster of one.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from ..core.contracts import (
    CLUSTERS_TSV,
    TREE_NWK,
    accession_from_filename,
    parse_genome_filename,
    read_clusters,
    record_name,
)
from ..core.errors import WorkdirError
from ..tree.rooting import read_newick

ITOL_DIR = "itol"
NAME_SOURCES = ("ncbi", "gtdb")
METADATA_SELECTED = "metadata_selected.tsv"
# Okabe-Ito first, then further distinguishable colours; grey after that.
_COLOURS = (
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
    "#882255",
    "#44AA99",
    "#117733",
    "#332288",
    "#AA4499",
    "#999933",
)
_OTHER_COLOUR = "#BBBBBB"
_SPECIES_RE = re.compile(r"s__([^;]+)")


@dataclass
class Leaf:
    name: str  # the leaf label in the tree
    accession: str
    species: str
    strain: str
    members: Counter[str] = field(default_factory=Counter)  # species -> genomes

    @property
    def genomes(self) -> int:
        return sum(self.members.values())


def short_species(name: str) -> str:
    """``Brucella abortus`` -> ``B. abortus``; ``Brucella sp. 2716`` ->
    ``Brucella sp.``; other forms unchanged."""
    parts = name.split()
    if len(parts) >= 2 and parts[1] == "sp.":
        return f"{parts[0]} sp."
    if len(parts) >= 2 and parts[0][:1].isupper():
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return name


def _species_from(row: Mapping[str, str], source: str) -> str:
    if source == "gtdb":
        m = _SPECIES_RE.search(row.get("gtdb_taxonomy", ""))
        return m.group(1).strip() if m else ""
    for column in ("ncbi_taxonomy", "ncbi_taxonomy_unfiltered"):
        m = _SPECIES_RE.search(row.get(column, ""))
        if m:
            return m.group(1).strip()
    return ""


def _metadata(workdir: Path) -> dict[str, dict[str, str]]:
    path = workdir / METADATA_SELECTED
    if not path.is_file():
        return {}
    with open(path, encoding="utf-8", newline="") as fo:
        return {row["accession"]: row for row in csv.DictReader(fo, delimiter="\t")}


def _name_from_file(filename: str) -> str:
    _family, genus, species, _acc = parse_genome_filename(filename)
    if not genus:
        return ""
    return f"{genus} {species}".strip() if species else f"{genus} sp."


def build_leaves(workdir: Path, source: str = "ncbi") -> list[Leaf]:
    """One :class:`Leaf` per leaf of ``tree/tree.nwk``, in tree order."""
    tree_path = workdir / "tree" / TREE_NWK
    if not tree_path.is_file():
        raise WorkdirError(f"No tree at {tree_path}; run phylo first.")
    tree = read_newick(tree_path.read_text(encoding="utf-8"))
    leaf_names = [n.taxon.label for n in tree.leaf_node_iter() if n.taxon is not None]
    meta = _metadata(workdir)
    clusters = read_clusters(workdir / "derep" / CLUSTERS_TSV)
    members_by_leaf = {record_name(rep): [rep, *ms] for rep, ms in clusters.items()}

    def species_of(filename: str) -> tuple[str, str, str]:
        acc = accession_from_filename(filename) or ""
        row = meta.get(acc)
        name = _species_from(row, source) if row else ""
        strain = (row or {}).get("ncbi_strain_identifiers", "").split(";")[0].strip()
        return acc, name or _name_from_file(filename) or "unknown", strain

    leaves = []
    for name in leaf_names:
        acc, species, strain = species_of(name)
        files = members_by_leaf.get(name, [name])
        members: Counter[str] = Counter()
        for f in files:
            members[short_species(species_of(f)[1])] += 1
        # A member listed with itself counts once.
        if files and files[0] in files[1:]:
            members[short_species(species)] -= 1
        leaves.append(Leaf(name, acc or name, short_species(species), strain, members))
    return leaves


def _colours(leaves: list[Leaf]) -> tuple[list[str], dict[str, str]]:
    totals: Counter[str] = Counter()
    for leaf in leaves:
        totals.update(leaf.members)
        totals[leaf.species] += 0
    order = [name for name, _ in totals.most_common()]
    colour = {
        name: _COLOURS[i] if i < len(_COLOURS) else _OTHER_COLOUR for i, name in enumerate(order)
    }
    return order, colour


def _legend(title: str, order: list[str], colour: Mapping[str, str]) -> str:
    tab = "\t"
    return (
        f"LEGEND_TITLE{tab}{title}\n"
        f"LEGEND_SHAPES{tab}{tab.join('1' for _ in order)}\n"
        f"LEGEND_COLORS{tab}{tab.join(colour[s] for s in order)}\n"
        f"LEGEND_LABELS{tab}{tab.join(order)}\n"
    )


def write_itol(workdir: Path, out_dir: Path | None = None, source: str = "ncbi") -> Path:
    """Write the annotation files; returns the folder."""
    if source not in NAME_SOURCES:
        raise ValueError(f"names must be one of {NAME_SOURCES}")
    leaves = build_leaves(workdir, source)
    out = out_dir or workdir / "tree" / ITOL_DIR
    out.mkdir(parents=True, exist_ok=True)
    order, colour = _colours(leaves)
    label = "NCBI" if source == "ncbi" else "GTDB"
    tab = "\t"

    with open(out / "labels.txt", "w", encoding="utf-8") as fo:
        fo.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for leaf in leaves:
            strain = f" {leaf.strain}" if leaf.strain else ""
            fo.write(f"{leaf.name}\t{leaf.species}{strain} ({leaf.accession}, n={leaf.genomes})\n")

    with open(out / "colorstrip_species.txt", "w", encoding="utf-8") as fo:
        fo.write("DATASET_COLORSTRIP\nSEPARATOR TAB\n")
        fo.write(f"DATASET_LABEL\t{label} species\nCOLOR\t#000000\nSTRIP_WIDTH\t25\n")
        fo.write(_legend(f"{label} species (representative)", order, colour))
        fo.write("DATA\n")
        for leaf in leaves:
            fo.write(f"{leaf.name}\t{colour[leaf.species]}\t{leaf.species}\n")

    with open(out / "bar_cluster_size.txt", "w", encoding="utf-8") as fo:
        fo.write("DATASET_SIMPLEBAR\nSEPARATOR TAB\n")
        fo.write("DATASET_LABEL\tGenomes in cluster\nCOLOR\t#0072B2\nWIDTH\t200\n")
        fo.write("SHOW_VALUE\t1\nDATA\n")
        for leaf in leaves:
            fo.write(f"{leaf.name}\t{leaf.genomes}\n")

    with open(out / "pie_cluster_species.txt", "w", encoding="utf-8") as fo:
        fo.write("DATASET_PIECHART\nSEPARATOR TAB\n")
        fo.write("DATASET_LABEL\tCluster composition\nCOLOR\t#000000\n")
        fo.write(f"FIELD_COLORS\t{tab.join(colour[s] for s in order)}\n")
        fo.write(f"FIELD_LABELS\t{tab.join(order)}\n")
        fo.write(_legend(f"{label} species of the cluster's genomes", order, colour))
        fo.write("DATA\n")
        for leaf in leaves:
            # Position -1 places the pie beside the leaf; the radius grows
            # with the square root of the cluster size, so the area with it.
            radius = round(5 + 3 * leaf.genomes**0.5, 1)
            counts = tab.join(str(leaf.members.get(s, 0)) for s in order)
            fo.write(f"{leaf.name}\t-1\t{radius}\t{counts}\n")

    with open(out / "clusters_species.tsv", "w", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(["leaf", "accession", "species", "strain", "genomes", "member_species"])
        for leaf in leaves:
            composition = "; ".join(f"{s} {n}" for s, n in leaf.members.most_common())
            writer.writerow(
                [leaf.name, leaf.accession, leaf.species, leaf.strain, leaf.genomes, composition]
            )
    return out
