"""Clusters at other ANI thresholds, estimated from the sketches, in the
dereplication step of the report."""

from __future__ import annotations

import gzip
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from repgenr.core.contracts import RefinedRow, write_clusters, write_genome_status, write_refined
from repgenr.dereplicators.sourmash import _greedy_cluster
from repgenr.reporting.ani_cutoffs import ani_matrix, cutoff_table, greedy_cluster_count
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.member_ani import KSIZE, pair_ani
from repgenr.reporting.registry import BY_NAME


def _sig(path: Path, hashes) -> None:
    payload = [{"signatures": [{"ksize": KSIZE, "molecule": "DNA", "mins": sorted(hashes)}]}]
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(payload).encode()))


def test_the_matrix_agrees_with_the_pairwise_estimate(tmp_path: Path) -> None:
    sets = {"a.fasta": range(0, 3000), "b.fasta": range(30, 3030), "c.fasta": range(9000, 9100)}
    for name, hashes in sets.items():
        _sig(tmp_path / f"{name[0]}.sig.zip", hashes)
    m = ani_matrix([*sets, "missing.fasta"], tmp_path)
    assert m.names == ["a.fasta", "b.fasta", "c.fasta"] and m.missing == ["missing.fasta"]
    expected, _ = pair_ani(frozenset(sets["a.fasta"]), frozenset(sets["b.fasta"]))
    assert m.ani[0, 1] == pytest.approx(expected) and m.ani[1, 0] == pytest.approx(expected)
    assert np.isnan(m.ani[0, 2])  # nothing shared: no estimate
    assert m.ani[2, 2] == 1.0


def test_the_count_follows_the_sourmash_dereplicator() -> None:
    rng = np.random.default_rng(3)
    for _ in range(20):
        n = 40
        x = rng.uniform(0.98, 1.0, (n, n))
        sim = (x + x.T) / 2
        sim[rng.uniform(size=(n, n)) < 0.05] = np.nan
        for t in (0.985, 0.99, 0.995):
            labels = [str(i) for i in range(n)]
            clusters, _ = _greedy_cluster(labels, sim, {lab: lab for lab in labels}, t)
            count, singletons, largest = greedy_cluster_count(sim, t)
            assert count == len(clusters)
            assert singletons == sum(1 for m in clusters.values() if not m)
            assert largest == 1 + max(len(m) for m in clusters.values())


def test_a_subset_is_counted_on_its_own(tmp_path: Path) -> None:
    for name, hashes in {"a": range(1000), "b": range(1000), "c": range(5000, 6000)}.items():
        _sig(tmp_path / f"{name}.sig.zip", hashes)
    m = ani_matrix(["a.fasta", "b.fasta", "c.fasta"], tmp_path)
    (row,) = cutoff_table(m, [0.99])
    assert (row.clusters, row.singletons, row.largest) == (2, 1, 2)
    (row,) = cutoff_table(m, [0.99], subset=["a.fasta", "c.fasta"])
    assert row.clusters == 2


def test_dereplication_step_tabulates_and_draws_the_cutoffs(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    (wd / "derep").mkdir(parents=True)
    (wd / "sketches").mkdir()
    names = ["r.fasta", "a.fasta", "b.fasta", "c.fasta"]
    write_clusters(wd / "derep" / "clusters.tsv", {"r.fasta": names[1:]})
    write_genome_status(
        wd / "derep" / "genome_status.tsv",
        {"r.fasta": "representative", **{n: "contained" for n in names[1:]}},
    )
    write_refined(
        wd / "derep" / "refined.tsv",
        [
            RefinedRow(n, "r.fasta", "r.fasta", n, "representative", "0.9999")
            for n in ("r.fasta", "a.fasta")
        ],
    )
    _sig(wd / "sketches" / "r.sig.zip", range(2000))
    _sig(wd / "sketches" / "a.sig.zip", range(2000))
    _sig(wd / "sketches" / "b.sig.zip", range(10, 2010))
    _sig(wd / "sketches" / "c.sig.zip", range(400, 2400))
    refresh_step(wd, BY_NAME["dereplication"], formats=("png",))
    frag = wd / "report" / "dereplication"
    rep = read_fragment(frag)
    rows = [line.split("\t") for line in (frag / "ani_cutoffs.tsv").read_text().splitlines()]
    assert rows[0][:4] == ["ani", "clusters", "singletons", "largest_cluster"]
    assert rows[0][4] == "clusters_lineage_r"
    by_ani = {float(r[0]): r for r in rows[1:]}
    assert 0.9999 in by_ani  # the refined lineage's ANI is tabulated
    assert by_ani[0.95][1] == "1"  # every genome within 95%
    # c shares 80% of the hashes: ANI 0.8^(1/31), about 0.9928.
    assert by_ani[0.99][1] == "1" and by_ani[0.995][1] == "2"
    assert "ani_cutoffs" in {f.id for f in rep.figures}
    assert (frag / "ani_cutoffs.png").is_file()
