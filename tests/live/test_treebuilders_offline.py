"""Alignment-free tree builders through the workdir phylo stage."""

from __future__ import annotations

from pathlib import Path

import pytest
from live_helpers import newick_leaves, representatives

from repgenr.core.contracts import TREE_NWK

pytestmark = [pytest.mark.live, pytest.mark.requires_binary("sourmash")]


@pytest.fixture
def derep_wd(synthetic_set, ingested_workdir, run_repgenr) -> tuple[Path, Path]:
    genomes = synthetic_set("clonal", n=8, length=50_000)
    wd = ingested_workdir(genomes)
    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    return genomes, wd


@pytest.mark.parametrize(
    "builder",
    [
        pytest.param("mashtree", marks=pytest.mark.requires_binary("mashtree")),
        pytest.param("sourmash", marks=pytest.mark.requires_binary("sourmash")),
    ],
)
def test_alignment_free_builder_on_representatives(run_repgenr, derep_wd, builder: str) -> None:
    _, wd = derep_wd
    run_repgenr("phylo", "-wd", wd, "--treebuilder", builder, "--no-outgroup", "-t", "2")
    tree = (wd / "tree" / TREE_NWK).read_text(encoding="utf-8")
    assert newick_leaves(tree) == {Path(r).stem for r in representatives(wd)}


@pytest.mark.parametrize(
    "builder",
    [
        pytest.param("mashtree", marks=pytest.mark.requires_binary("mashtree")),
        pytest.param("sourmash", marks=pytest.mark.requires_binary("sourmash")),
    ],
)
def test_alignment_free_builder_roots_at_the_midpoint(run_repgenr, derep_wd, builder: str) -> None:
    import dendropy

    _, wd = derep_wd
    run_repgenr("phylo", "-wd", wd, "--treebuilder", builder, "--no-outgroup", "--root", "midpoint")
    text = (wd / "tree" / TREE_NWK).read_text(encoding="utf-8")
    tree = dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True)
    assert len(tree.seed_node.child_nodes()) == 2
    assert newick_leaves(text) == {Path(r).stem for r in representatives(wd)}
    unrooted = (wd / "tree" / "tree.unrooted.nwk").read_text(encoding="utf-8")
    assert newick_leaves(unrooted) == newick_leaves(text)


@pytest.mark.requires_binary("mashtree")
def test_all_genomes_puts_every_genome_in_the_tree(run_repgenr, derep_wd) -> None:
    genomes, wd = derep_wd
    run_repgenr("phylo", "-wd", wd, "--treebuilder", "mashtree", "--no-outgroup", "--all-genomes")
    tree = (wd / "tree" / TREE_NWK).read_text(encoding="utf-8")
    assert newick_leaves(tree) == {f.stem for f in genomes.glob("*.fasta")}


@pytest.mark.requires_binary("rapidnj")
def test_rapidnj_joins_a_matrix_with_the_labels_intact(tmp_path: Path) -> None:
    import logging

    from repgenr.tree.nj import nj_tree

    labels = [
        "Fam_Gen_sp_GCA_000001.1",
        "Fam_Gen_sp_GCA_000002.1",
        "Out_Gen_sp_GCA_000009.1",
        "a b",
    ]
    dist = [
        [0.0, 0.1, 0.8, 0.9],
        [0.1, 0.0, 0.8, 0.9],
        [0.8, 0.8, 0.0, 0.3],
        [0.9, 0.9, 0.3, 0.0],
    ]
    newick = nj_tree(labels, dist, tmp_path, logging.getLogger("live"), backend="rapidnj")
    import dendropy

    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    assert {t.label for t in tree.taxon_namespace} == set(labels)
    assert "'a b'" in newick  # quoted by restore_labels, not mangled by rapidnj
