"""Unit tests for the sourmash greedy clustering (numpy-vectorized)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepParams,
)
from repgenr.dereplicators.sourmash import (
    _greedy_cluster,
    _parse_pairwise_csv,
    _sparse_greedy_cluster,
)


def _name_map(labels):
    return {label: label for label in labels}


def test_two_disjoint_pairs() -> None:
    labels = ["a", "b", "c", "d"]
    # a~b and c~d are similar; cross-pairs are not.
    sim = [
        [1.0, 0.999, 0.1, 0.1],
        [0.999, 1.0, 0.1, 0.1],
        [0.1, 0.1, 1.0, 0.999],
        [0.1, 0.1, 0.999, 1.0],
    ]
    clusters, status = _greedy_cluster(labels, sim, _name_map(labels), 0.99)
    # one representative per pair; members are the other element
    assert len(clusters) == 2
    assert all(len(m) == 1 for m in clusters.values())
    assert set(clusters) == {"a", "c"}  # first index of each pair wins on ties
    assert clusters["a"] == ["b"] and clusters["c"] == ["d"]
    assert status["a"] == STATUS_REPRESENTATIVE and status["b"] == STATUS_CONTAINED


def test_all_connected_collapses_to_one() -> None:
    labels = ["a", "b", "c"]
    sim = np.array([[1.0, 0.999, 0.999], [0.999, 1.0, 0.999], [0.999, 0.999, 1.0]])
    clusters, status = _greedy_cluster(labels, sim, _name_map(labels), 0.99)
    assert len(clusters) == 1
    rep = next(iter(clusters))
    assert sorted(clusters[rep]) == ["b", "c"]


def test_all_distinct_keeps_all() -> None:
    labels = ["a", "b", "c"]
    sim = np.eye(3)  # nothing above threshold off-diagonal
    clusters, status = _greedy_cluster(labels, sim, _name_map(labels), 0.99)
    assert len(clusters) == 3
    assert all(m == [] for m in clusters.values())
    assert all(s == STATUS_REPRESENTATIVE for s in status.values())


# --- sparse (branchwater pairwise) path ---------------------------------------

_PAIRWISE_HEADER = (
    "query_name,query_md5,match_name,match_md5,containment,max_containment,"
    "jaccard,average_containment,intersect_hashes,ksize,scaled,moltype,cosine\n"
)


def _pairwise_csv(tmp_path, rows):
    # rows: (query, match, containment) -- jaccard is set to the same value so
    # the file looks like a real branchwater edge list.
    path = tmp_path / "pairwise.csv"
    lines = [_PAIRWISE_HEADER]
    for q, m, c in rows:
        lines.append(f"{q},md5q,{m},md5m,{c},{c},{c},{c},10,31,1000,DNA,\n")
    path.write_text("".join(lines))
    return path


def test_parse_pairwise_csv_thresholds_and_symmetrizes(tmp_path) -> None:
    known = {"a", "b", "c", "d"}
    csv_path = _pairwise_csv(
        tmp_path,
        [
            ("a", "b", 0.999),  # kept (ANI ~ 1.0)
            ("a", "a", 1.0),  # self-edge dropped
            ("c", "d", 0.5),  # containment 0.5 -> ANI ~ 0.978 < 0.99, dropped
            ("a", "z", 0.999),  # unknown node dropped
        ],
    )
    neighbors = _parse_pairwise_csv(csv_path, 0.99, known, ksize=31)
    assert neighbors == {"a": {"b"}, "b": {"a"}}  # symmetric, single kept edge


def test_parse_pairwise_thresholds_ani_not_jaccard(tmp_path) -> None:
    # Measured on synthetic genomes (skani ground truth 99.57 percent ANI):
    # k=31 jaccard 0.796, average_containment 0.8866 -> ANI estimate 0.9961.
    # The pair must be KEPT at --secondary-ani 0.99 although jaccard < 0.99.
    path = tmp_path / "pairwise.csv"
    path.write_text(
        _PAIRWISE_HEADER
        + "a,md5,b,md5,0.8833,0.8898,0.7962,0.8866,1825,31,1000,DNA,\n"
        + "b,md5,c,md5,0.2574,0.2574,0.1476,0.2572,528,31,1000,DNA,\n"
    )
    neighbors = _parse_pairwise_csv(path, 0.99, {"a", "b", "c"}, ksize=31)
    assert neighbors == {"a": {"b"}, "b": {"a"}}  # 94.96 percent ANI pair dropped


def test_parse_pairwise_falls_back_to_jaccard_conversion(tmp_path) -> None:
    # An edge list without the average_containment column: convert jaccard to
    # its containment equivalent 2j/(1+j) before applying the ANI threshold.
    path = tmp_path / "pairwise.csv"
    path.write_text(
        "query_name,query_md5,match_name,match_md5,jaccard,ksize\n"
        + "a,md5,b,md5,0.7962,31\n"  # 2j/(1+j)=0.8865 -> ANI 0.9961: kept
        + "b,md5,c,md5,0.1476,31\n"  # -> ANI ~ 0.957: dropped
    )
    neighbors = _parse_pairwise_csv(path, 0.99, {"a", "b", "c"}, ksize=31)
    assert neighbors == {"a": {"b"}, "b": {"a"}}


def test_sparse_two_disjoint_pairs() -> None:
    labels = ["a", "b", "c", "d"]
    neighbors = {"a": {"b"}, "b": {"a"}, "c": {"d"}, "d": {"c"}}
    clusters, status = _sparse_greedy_cluster(labels, neighbors, _name_map(labels))
    assert set(clusters) == {"a", "c"}  # first of each pair wins on ties
    assert clusters["a"] == ["b"] and clusters["c"] == ["d"]
    assert status["a"] == STATUS_REPRESENTATIVE and status["b"] == STATUS_CONTAINED


def test_sparse_loner_stays_singleton() -> None:
    labels = ["a", "b", "c"]
    neighbors = {"a": {"b"}, "b": {"a"}}  # c has no edges
    clusters, status = _sparse_greedy_cluster(labels, neighbors, _name_map(labels))
    assert set(clusters) == {"a", "c"}
    assert clusters["c"] == []
    assert status["c"] == STATUS_REPRESENTATIVE


def test_sparse_matches_dense_on_same_graph() -> None:
    # Build one similarity matrix; derive both a dense matrix and the equivalent
    # sparse adjacency from it, and confirm both clusterers agree.
    labels = ["a", "b", "c", "d", "e"]
    edges = {("a", "b"), ("b", "c"), ("d", "e")}  # a-b-c chain + d-e pair
    thr = 0.99
    n = len(labels)
    idx = {lab: i for i, lab in enumerate(labels)}
    sim = np.eye(n)
    neighbors: dict[str, set[str]] = {}
    for u, v in edges:
        sim[idx[u], idx[v]] = sim[idx[v], idx[u]] = 0.999
        neighbors.setdefault(u, set()).add(v)
        neighbors.setdefault(v, set()).add(u)

    dense, _ = _greedy_cluster(labels, sim, _name_map(labels), thr)
    sparse, _ = _sparse_greedy_cluster(labels, neighbors, _name_map(labels))
    assert {k: sorted(v) for k, v in dense.items()} == {k: sorted(v) for k, v in sparse.items()}


# --- threshold semantics on the tool invocations ------------------------------


def test_dense_thresholds_ani_derived_from_containment(tmp_path, monkeypatch) -> None:
    """The dense path reads ``compare --avg-containment`` and converts the
    average containment to ANI estimates, c^(1/k), as the sparse path does.

    The measured 99.57 percent ANI pair has k=31 Jaccard 0.796, average
    containment 2j/(1+j) = 0.8864 for sketches of one size -> ANI estimate
    0.9961, so it must merge at threshold 0.99 although its containment is
    far below 0.99. An identical pair must survive conversion too --
    sourmash's own ``compare --ani`` reports 0 for tiny sketches, which is why
    the adapter converts locally instead.
    """
    import logging

    from repgenr.dereplicators import sourmash as sm

    genomes = []
    for name in ("g1.fasta", "g2.fasta", "g3.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        genomes.append(p)

    argvs: list[list[str]] = []

    def fake_run_tool(caps, argv, **kwargs):
        argvs.append([str(a) for a in argv])
        return None

    containment = [
        [1.0, 0.8864, 0.0],  # g1-g2: the measured 99.57 pct ANI pair
        [0.8864, 1.0, 0.0],
        [0.0, 0.0, 1.0],  # g3 unrelated
    ]
    monkeypatch.setattr(sm, "run_tool", fake_run_tool)
    monkeypatch.setattr(
        sm,
        "_find_signatures",
        lambda sig_dir, gs: {g: tmp_path / f"{g.name}.sig" for g in gs},
    )
    monkeypatch.setattr(
        sm,
        "_read_compare_csv",
        lambda path: (["g1.fasta", "g2.fasta", "g3.fasta"], containment),
    )

    adapter = sm.SourmashDereplicator()
    clusters, _status = adapter._dense_dereplicate(
        genomes, tmp_path / "out", 31, 1000, 0.99, logging.getLogger("t")
    )
    assert set(clusters) == {"g1.fasta", "g3.fasta"}  # g2 merged under g1
    (compare,) = [a for a in argvs if a[:2] == ["sourmash", "compare"]]
    assert "--avg-containment" in compare
    assert compare[compare.index("-k") + 1] == "31"


def test_containment_to_ani_matrix_values() -> None:
    from repgenr.dereplicators.sourmash import _containment_to_ani_matrix

    ani, n_unknown = _containment_to_ani_matrix(np.array([[1.0, 0.8864, 0.0]]), 31)
    assert n_unknown == 0  # no sizes: no shared-hash minimum
    assert ani[0, 0] == 1.0  # identical content stays 1.0
    assert abs(ani[0, 1] - 0.9961) < 5e-4
    assert ani[0, 2] == 0.0  # disjoint stays 0


def test_dense_and_sparse_back_ends_give_one_clustering(tmp_path) -> None:
    """For one matrix of average containment the dense greedy pick and the
    sparse one over the pairwise edge list give the same clusters."""
    from repgenr.dereplicators.sourmash import (
        _containment_to_ani_matrix,
        _greedy_cluster,
        _parse_pairwise_csv,
        _sparse_greedy_cluster,
    )

    labels = ["a", "b", "c", "d"]
    # a-b and c-d close; a 0.6498 fragment pair (ANI 0.986) between b and c.
    c = np.array(
        [
            [1.0, 0.90, 0.10, 0.05],
            [0.90, 1.0, 0.6498, 0.05],
            [0.10, 0.6498, 1.0, 0.95],
            [0.05, 0.05, 0.95, 1.0],
        ]
    )
    names = {x: f"{x}.fasta" for x in labels}
    for threshold in (0.95, 0.99):
        ani, _n = _containment_to_ani_matrix(c, 31)
        dense = _greedy_cluster(labels, ani, names, threshold)
        rows = ["query_name,match_name,average_containment"]
        rows += [f"{labels[i]},{labels[j]},{c[i, j]}" for i in range(4) for j in range(4) if i < j]
        edges = tmp_path / f"pairwise-{threshold}.csv"
        edges.write_text("\n".join(rows) + "\n", encoding="utf-8")
        sparse = _sparse_greedy_cluster(
            labels, _parse_pairwise_csv(edges, threshold, set(labels), 31), names
        )
        assert dense == sparse


def test_sparse_pairwise_prefilter_is_containment_scaled(tmp_path, monkeypatch) -> None:
    """`pairwise -t` is a containment threshold: pass threshold**ksize, not the
    ANI threshold (which would drop every edge below ~99 percent containment)."""
    import logging

    from repgenr.dereplicators import sourmash as sm

    genomes = []
    for name in ("g1.fasta", "g2.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        genomes.append(p)

    calls: list[list[str]] = []

    def fake_run_tool(caps, argv, **kwargs):
        argv = [str(a) for a in argv]
        calls.append(argv)
        if "-o" in argv:  # create whatever output the step promises
            out = Path(argv[argv.index("-o") + 1])
            out.parent.mkdir(parents=True, exist_ok=True)
            header = "query_name,query_md5,match_name,match_md5,jaccard,ksize\n"
            out.write_text(header, encoding="utf-8")

    monkeypatch.setattr(sm, "run_tool", fake_run_tool)

    (tmp_path / "out").mkdir()  # dereplicate() creates out_dir before this call
    adapter = sm.SourmashDereplicator()
    adapter._sparse_dereplicate(
        genomes,
        tmp_path / "out",
        31,
        1000,
        0.99,
        DerepParams(threads=2),
        logging.getLogger("t"),
    )
    pairwise_argv = next(argv for argv in calls if "pairwise" in argv)
    t_value = float(pairwise_argv[pairwise_argv.index("-t") + 1])
    assert abs(t_value - 0.99**31) < 1e-6


# --- the back-end in DerepResult.settings ---------------------------------------


def _settings_run(tmp_path, monkeypatch, *, sparse: bool, sparse_fails: bool = False, extra=None):
    import logging

    from repgenr.core.errors import ToolExecutionError
    from repgenr.dereplicators import sourmash as sm
    from repgenr.dereplicators.base import DerepParams

    genomes = []
    for name in ("g1.fasta", "g2.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        genomes.append(p)
    one = ({"g1.fasta": ["g2.fasta"]}, {"g1.fasta": "representative", "g2.fasta": "contained"})

    def fake_sparse(*args, **kwargs):
        if sparse_fails:
            raise ToolExecutionError(["sourmash"], 1, "pairwise failed", tool="sourmash")
        return one

    monkeypatch.setattr(sm, "_branchwater_available", lambda caps, logger: sparse)
    monkeypatch.setattr(sm.SourmashDereplicator, "_sparse_dereplicate", fake_sparse)
    monkeypatch.setattr(sm.SourmashDereplicator, "_dense_cluster", lambda *a, **k: (*one, 3))
    params = DerepParams(secondary_ani=0.95, extra={"ksize": 21, "scaled": 500, **(extra or {})})
    return sm.SourmashDereplicator().dereplicate(
        genomes, tmp_path / "out", params, logging.getLogger("t")
    )


def test_settings_name_the_sparse_back_end(tmp_path, monkeypatch) -> None:
    result = _settings_run(tmp_path, monkeypatch, sparse=True)
    assert result.settings == {
        "backend": "sparse",
        "ani_estimate": "average_containment^(1/k)",
        "ksize": 21,
        "scaled": 500,
    }


def test_settings_name_the_dense_back_end_and_unknown_pairs(tmp_path, monkeypatch) -> None:
    result = _settings_run(tmp_path, monkeypatch, sparse=False)
    assert result.settings == {
        "backend": "dense",
        "ani_estimate": "average_containment^(1/k)",
        "ksize": 21,
        "scaled": 500,
        "unknown_pairs": 3,
    }


def test_settings_follow_the_dense_fallback(tmp_path, monkeypatch) -> None:
    result = _settings_run(
        tmp_path, monkeypatch, sparse=True, sparse_fails=True, extra={"dense_fallback": True}
    )
    assert result.settings["backend"] == "dense"
    assert result.settings["unknown_pairs"] == 3


def test_dense_cluster_counts_unknown_pairs(tmp_path, monkeypatch) -> None:
    import logging

    from repgenr.dereplicators import sourmash as sm

    ani = np.array([[1.0, np.nan, 0.99], [np.nan, 1.0, np.nan], [0.99, np.nan, 1.0]])
    names = {x: x for x in ("a", "b", "c")}
    monkeypatch.setattr(
        sm.SourmashDereplicator,
        "_dense_ani_matrix",
        lambda *a, **k: (["a", "b", "c"], ani, names),
    )
    clusters, _status, unknown = sm.SourmashDereplicator()._dense_cluster(
        [], tmp_path, 31, 1000, 0.95, logging.getLogger("t")
    )
    assert unknown == 2
    assert clusters == {"a": ["c"], "b": []}
