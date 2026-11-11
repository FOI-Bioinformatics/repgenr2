"""Smoke test: the offline chain on synthetic genomes with lightweight tools.

Mirrors ``nextflow/tests/local_dataflow.nf`` (sourmash + mashtree, no outgroup)
through the CLI: ingest -> dereplicate -> phylo -> tree2tax, then a second
pass must skip every stage on the resume records.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from repgenr.core.contracts import (
    CLUSTERS_TSV,
    GENOMES_MAP_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    read_clusters,
)

pytestmark = [pytest.mark.live, pytest.mark.requires_binary("sourmash", "mashtree")]


def _leaves(newick: str) -> set[str]:
    out = set()
    token = ""
    for ch in newick:
        if ch in "(),;":
            if token:
                out.add(token.split(":")[0])
            token = ""
        else:
            token += ch
    return out


def test_offline_chain_sourmash_mashtree(
    run_repgenr, synthetic_set, ingested_workdir, tmp_path: Path
) -> None:
    genomes = synthetic_set("clonal", n=8, length=50_000)
    truth = json.loads((genomes / "truth.json").read_text(encoding="utf-8"))
    wd = ingested_workdir(genomes)

    run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    run_repgenr("phylo", "-wd", wd, "--treebuilder", "mashtree", "--no-outgroup", "-t", "2")
    run_repgenr("tree2tax", "-wd", wd, "--include-dereplicated")

    clusters = read_clusters(wd / "derep" / CLUSTERS_TSV)
    partition = {rep: {rep, *members} for rep, members in clusters.items()}
    assert set().union(*partition.values()) == set(truth["clusters"]), (
        "every ingested genome is assigned to exactly one cluster"
    )
    expected = {}
    for name, cluster in truth["clusters"].items():
        expected.setdefault(cluster, set()).add(name)
    assert sorted(map(sorted, partition.values())) == sorted(map(sorted, expected.values())), (
        "sourmash clustering must recover the synthetic partition"
    )

    reps = {Path(rep).stem for rep in clusters}
    tree = (wd / "tree" / TREE_NWK).read_text(encoding="utf-8")
    assert _leaves(tree) == reps, "mashtree leaves are the representatives"

    rows = (wd / TREE2TAX_TSV).read_text(encoding="utf-8").splitlines()[1:]
    assert {row.split("\t")[0] for row in rows} == reps
    mapped = (wd / GENOMES_MAP_TSV).read_text(encoding="utf-8").splitlines()
    assert len(mapped) == len(truth["clusters"]), (
        "--include-dereplicated maps every genome to its representative"
    )

    status = run_repgenr("status", "-wd", wd).stdout
    assert "Pipeline: local" in status
    for stage in ("ingest", "dereplicate", "phylo", "tree2tax"):
        assert f"[done]    {stage}" in status

    second = run_repgenr("dereplicate", "-wd", wd, "--tool", "sourmash", "-t", "2")
    assert "skipping" in second.stderr + second.stdout

    # The stages wrote the selection report; the command writes the page.
    assert (wd / "report" / "selection" / "step.json").is_file()
    run_repgenr("report", "-wd", wd)
    page = (wd / "report" / "report.html").read_text(encoding="utf-8")
    for text in (">Genome selection</h2>", "Figure 1.", ">Methods</h2>"):
        assert text in page

    # iTOL files: one bar per leaf, the bars summing to every genome.
    run_repgenr("itol", "-wd", wd)
    bars = (wd / "tree" / "itol" / "bar_cluster_size.txt").read_text(encoding="utf-8")
    sizes = [int(line.split("\t")[1]) for line in bars.split("DATA\n", 1)[1].splitlines()]
    assert len(sizes) == len(reps) and sum(sizes) == len(truth["clusters"])
    run_repgenr("itol", "-wd", wd, "--names", "gtdb", "--out", tmp_path / "itol")
    assert (tmp_path / "itol" / "pie_cluster_species.txt").is_file()


def test_run_local_chain(run_repgenr, synthetic_set, tmp_path: Path) -> None:
    """`run --genomes-dir` is the ingest chain in one command, and a second run skips."""
    genomes = synthetic_set("clonal", n=8, length=50_000)
    wd = tmp_path / "wd"
    args = (
        *("run", "-wd", wd, "--genomes-dir", genomes, "--tool", "sourmash"),
        *("--treebuilder", "mashtree", "--no-outgroup", "-t", "2"),
    )
    run_repgenr(*args)
    assert (wd / "tree" / TREE_NWK).exists()
    assert (wd / TREE2TAX_TSV).exists()
    status = run_repgenr("status", "-wd", wd).stdout
    assert "Pipeline: local" in status
    for stage in ("ingest", "dereplicate", "phylo", "tree2tax"):
        assert f"[done]    {stage}" in status
    second = run_repgenr(*args)
    assert (second.stdout + second.stderr).count("skipping") == 4
