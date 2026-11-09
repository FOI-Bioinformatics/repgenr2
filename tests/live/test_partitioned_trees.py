"""Partitioned tree builders with the real binaries (phylo-build --partitioned).

A seeded two-gene alignment of six genomes in three pairs, with the genes
evolving at different rates. The partitioned model must be passed to the
tool (one partition per gene, the unparseable gene name rewritten) and the
tree must recover the three pairs.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest
from live_helpers import newick_leaves

from repgenr.coregenes.select import Partition
from repgenr.tree.partitions import write_partitions

pytestmark = [pytest.mark.live]

_PAIRS = (("gA", "gB"), ("gC", "gD"), ("gE", "gF"))
# Gene names as a core-gene builder may give them; "|" and "." are not
# accepted in a RAxML-NG or NEXUS partition name and are rewritten.
_GENES = (("dnaK|1", 900, 0.03), ("rpoB.2", 1200, 0.10))


def _mutate(seq: list[str], rate: float, rng: random.Random) -> list[str]:
    return [rng.choice([b for b in "ACGT" if b != c]) if rng.random() < rate else c for c in seq]


def _simulate(out: Path) -> tuple[Path, Path]:
    """Write the alignment, its partitions.json and one genome file per record."""
    rng = random.Random(7)
    rows: dict[str, list[str]] = {name: [] for pair in _PAIRS for name in pair}
    partitions: list[Partition] = []
    start = 1
    for gene, length, rate in _GENES:
        root = [rng.choice("ACGT") for _ in range(length)]
        for pair in _PAIRS:
            ancestor = _mutate(root, rate, rng)
            for name in pair:
                rows[name] += _mutate(ancestor, rate / 3, rng)
        partitions.append(Partition(gene, start, start + length - 1))
        start += length
    out.mkdir(parents=True)
    msa = out / "core_gene_alignment.fasta"
    msa.write_text("".join(f">{n}\n{''.join(s)}\n" for n, s in rows.items()), encoding="utf-8")
    write_partitions(msa, partitions, "dna")
    genomes = out / "genomes"
    genomes.mkdir()
    for name, seq in rows.items():
        (genomes / f"{name}.fasta").write_text(f">{name}\n{''.join(seq)}\n", encoding="utf-8")
    return msa, genomes


def _pairs_recovered(newick: str) -> bool:
    import dendropy

    tree = dendropy.Tree.get(data=newick, schema="newick", rooting="force-unrooted")
    tree.encode_bipartitions()
    labels = {t.label: t for t in tree.taxon_namespace}
    splits = {b.split_bitmask for b in tree.bipartition_encoding}
    for a, b in _PAIRS:
        mask = tree.taxon_namespace.taxa_bitmask(taxa=[labels[a], labels[b]])
        if mask not in splits and (tree.taxon_namespace.all_taxa_bitmask() ^ mask) not in splits:
            return False
    return True


@pytest.mark.parametrize(
    ("builder", "bootstrap", "log", "marker"),
    [
        pytest.param(
            "raxmlng",
            "10",
            "raxml.raxml.log",
            "2 partitions",
            marks=pytest.mark.requires_binary("raxml-ng"),
        ),
        pytest.param(
            "iqtree",
            "1000",
            "msa.fasta.log",
            "Reading partition model file",
            marks=pytest.mark.requires_binary("iqtree3"),
        ),
    ],
)
def test_partitioned_builder_fits_one_model_per_gene(
    run_repgenr, tmp_path: Path, builder: str, bootstrap: str, log: str, marker: str
) -> None:
    msa, genomes = _simulate(tmp_path / "in")
    out = tmp_path / "out"
    run_repgenr(
        "phylo-build", "--genomes-dir", genomes, "--msa", msa, "-o", out,
        "--treebuilder", builder, "--partitioned", "-B", bootstrap,
        "--no-outgroup", "--root", "none", "-t", "2",
    )  # fmt: skip
    tree_dir = out / "tree"
    part_file = (tree_dir / "partitions.txt").read_text(encoding="utf-8")
    # Both genes, under names the tools accept, with the alignment's ranges.
    for raw in ("dnaK|1", "rpoB.2"):
        assert raw not in part_file
    assert "1-900" in part_file
    assert "901-2100" in part_file
    assert marker in (tree_dir / log).read_text(encoding="utf-8")

    newick = (tree_dir / "tree.nwk").read_text(encoding="utf-8")
    assert newick_leaves(newick) == {n for pair in _PAIRS for n in pair}
    assert _pairs_recovered(newick)
