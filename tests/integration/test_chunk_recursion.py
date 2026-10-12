"""Stop rule, interleaving and composition of the chunked dereplication levels.

A fake dereplicator clusters genomes by a hidden label: every genome with the
same label falls in one cluster, whose representative is the first such genome
in the input order. The number of tool calls, and the size of each, stands in
for the run time, so the tests need no timing.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
    check_result_complete,
    registry,
)
from repgenr.stages import dereplicate as derep_stage
from repgenr.stages.dereplicate import DereplicateParams, _dereplicate_chunked, run


class _Labelled(Dereplicator):
    """Clusters by ``labels[name]``; genomes in ``fail_names`` fail QC."""

    capabilities = ToolCapabilities(name="labelled", supports_native_scaling=True)

    def __init__(self, labels: dict[str, str], fail_names: frozenset[str] = frozenset()):
        self.labels = labels
        self.fail_names = fail_names
        self.calls: list[int] = []

    def preflight(self) -> dict[str, str]:
        return {"labelled": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        self.calls.append(len(genomes))
        status: dict[str, str] = {}
        clusters: dict[str, list[str]] = {}
        first: dict[str, Path] = {}
        reps: list[Path] = []
        for g in genomes:
            if g.name in self.fail_names:
                status[g.name] = STATUS_FAIL_QC
                continue
            label = self.labels[g.name]
            if label not in first:
                first[label] = g
                reps.append(g)
                clusters[g.name] = []
                status[g.name] = STATUS_REPRESENTATIVE
            else:
                clusters[first[label].name].append(g.name)
                status[g.name] = STATUS_CONTAINED
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


def _genomes(n: int) -> list[Path]:
    return [Path(f"/g/Fam_g_s_GCA_{i:06d}.1.fasta") for i in range(n)]


def _params(threads: int = 4) -> DerepParams:
    return DerepParams(primary_ani=0.9, secondary_ani=0.99, threads=threads)


def _chunked(adapter, genomes, process_size, tmp_path, logger=None) -> DerepResult:  # noqa: ANN001
    p = _params()
    return _dereplicate_chunked(
        adapter,
        genomes,
        tmp_path,
        p,
        process_size,
        2,
        p.primary_ani,
        p.secondary_ani,
        logger or logging.getLogger("test_chunk_recursion"),
    )


def _membership(result: DerepResult, labels: dict[str, str]) -> set[frozenset[str]]:
    """Clusters as sets of genome names, independent of the representative."""
    groups = {frozenset([rep, *members]) for rep, members in result.clusters.items()}
    # every cluster is one label, and each label is one cluster
    for group in groups:
        assert len({labels[g] for g in group}) == 1
    return groups


def test_singletons_one_level_then_final_pass(tmp_path: Path, caplog) -> None:
    """A diverse set: no genome contains another at the final threshold.

    The union equals the input, so recursion would only re-chunk the same set.
    One level of chunks and one final pass, nothing more. This is the size of
    the reported set (6888 genomes, -s 500); the reported run itself removed
    3081 genomes at level 0 and is expected to run one further level.
    """
    genomes = _genomes(6888)
    labels = {g.name: g.name for g in genomes}
    adapter = _Labelled(labels)
    with caplog.at_level(logging.INFO, logger="test_chunk_recursion"):
        res = _chunked(adapter, genomes, 500, tmp_path)

    assert len(adapter.calls) == math.ceil(6888 / 500) + 1
    assert max(adapter.calls[:-1]) <= 500
    assert adapter.calls[-1] == 6888  # the final pass sees the whole union
    assert len(res.representatives) == 6888
    check_result_complete(res, labels)
    messages = [r.getMessage() for r in caplog.records]
    assert any(
        m == "Final pass over 6888 representatives (chunk size 500 bounds the leaf chunks, "
        "not the merge)"
        for m in messages
    )
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "6888" in warnings[0] and "--process-size" in warnings[0]
    assert "--pre-secondary-ani" in warnings[0]
    assert not any("Level 1" in m for m in messages)


def test_clonal_clusters_across_chunks_collapse_after_interleaving(tmp_path: Path) -> None:
    """Large clusters spread over every chunk: one interleaved level gathers them.

    2000 genomes in 50 labels (label = index mod 50), chunk size 100. Every leaf
    chunk holds all 50 labels, so level 0 leaves 20 x 50 = 1000 representatives.
    In the fixed order a level-1 chunk would again pair only two chunks' worth;
    interleaved, a level-1 chunk holds five labels from twenty chunks and
    collapses to five genomes, so level 1 ends at 50 and the final pass is small.
    """
    genomes = _genomes(2000)
    labels = {g.name: str(i % 50) for i, g in enumerate(genomes)}
    adapter = _Labelled(labels)
    res = _chunked(adapter, genomes, 100, tmp_path)

    single = _Labelled(labels).dereplicate(genomes, tmp_path / "single", _params(), None)
    # level 0: 20 chunks; level 1: 10 chunks of the interleaved 1000; final pass of 50
    assert len(adapter.calls) == 20 + 10 + 1
    assert adapter.calls[-1] == 50
    assert sorted(p.name for p in res.representatives) == sorted(
        p.name for p in single.representatives
    )
    assert _membership(res, labels) == _membership(single, labels)
    check_result_complete(res, labels)


def test_small_shrink_stops_recursion(tmp_path: Path, caplog) -> None:
    """A level that removes less than a chunk's worth does not recurse.

    1000 genomes, chunk size 100, 40 adjacent pairs share a label: level 0
    removes 40 genomes (fewer than max(100, 100)), so the 960 representatives go
    to one final pass instead of a level-by-level descent.
    """
    genomes = _genomes(1000)
    labels = {g.name: g.name for g in genomes}
    for i in range(0, 80, 2):
        labels[genomes[i + 1].name] = labels[genomes[i].name]
    adapter = _Labelled(labels)
    with caplog.at_level(logging.INFO, logger="test_chunk_recursion"):
        res = _chunked(adapter, genomes, 100, tmp_path)

    assert len(adapter.calls) == 10 + 1
    assert adapter.calls[-1] == 960
    assert len(res.representatives) == 960
    assert any("another level needs >= 100" in r.getMessage() for r in caplog.records)
    check_result_complete(res, labels)


def test_shrink_threshold_is_a_tenth_of_large_inputs() -> None:
    assert derep_stage._recurse_min_shrink(6888, 500) == 689
    assert derep_stage._recurse_min_shrink(3807, 500) == 500
    assert derep_stage._recurse_min_shrink(10, 2) == 2


def test_depth_backstop_ends_in_final_pass(tmp_path: Path, monkeypatch, caplog) -> None:
    """At the depth limit the union goes to the final pass whatever the shrink."""
    monkeypatch.setattr(derep_stage, "_MAX_REDUCE_DEPTH", 0)
    genomes = _genomes(2000)
    labels = {g.name: str(i % 50) for i, g in enumerate(genomes)}
    adapter = _Labelled(labels)
    with caplog.at_level(logging.INFO, logger="test_chunk_recursion"):
        res = _chunked(adapter, genomes, 100, tmp_path)

    assert len(adapter.calls) == 20 + 1
    assert adapter.calls[-1] == 1000
    assert len(res.representatives) == 50
    assert any("depth limit 0 reached" in r.getMessage() for r in caplog.records)
    check_result_complete(res, labels)


def test_no_warning_when_final_pass_is_small(tmp_path: Path, caplog) -> None:
    genomes = _genomes(2000)
    labels = {g.name: str(i % 50) for i, g in enumerate(genomes)}
    with caplog.at_level(logging.INFO, logger="test_chunk_recursion"):
        _chunked(_Labelled(labels), genomes, 100, tmp_path)
    assert not [r for r in caplog.records if r.levelno == logging.WARNING]
    assert any(
        r.getMessage().startswith("Final pass over 50 representatives") for r in caplog.records
    )


def test_interleave_by_origin_round_robin() -> None:
    a, b, c = (_genomes(7)[i : i + 3] for i in (0, 3, 6))
    results = [
        DerepResult(representatives=list(x), clusters={}, genome_status={}) for x in (a, b, c)
    ]
    order = derep_stage._interleave_by_origin(results)
    assert order == [a[0], b[0], c[0], a[1], b[1], a[2], b[2]]


def test_composition_with_qc_failures_over_three_levels(tmp_path: Path) -> None:
    """QC failures keep their status and every genome maps to one representative.

    32 genomes in 4 labels of 8 consecutive genomes, chunk size 2: each level
    pairs genomes of one label, so levels 0, 1 and 2 halve the set (32, 16, 8,
    4). Level 2 removed more than a chunk's worth, so level 3 chunks the four
    labels once more, removes none, and the final pass receives them.
    """
    genomes = _genomes(32)
    labels = {g.name: str(i // 8) for i, g in enumerate(genomes)}
    fail = frozenset({genomes[5].name, genomes[20].name})
    adapter = _Labelled(labels, fail_names=fail)
    res = _chunked(adapter, genomes, 2, tmp_path)

    assert adapter.calls == [2] * 16 + [2] * 8 + [2] * 4 + [2] * 2 + [4]
    for name in fail:
        assert res.genome_status[name] == STATUS_FAIL_QC
        assert name not in {p.name for p in res.representatives}
    assert len(res.genome_status) == 32
    assert len(res.representatives) == 4
    check_result_complete(res, labels)
    kept = {k: v for k, v in labels.items() if k not in fail}
    assert _membership(res, kept) == {
        frozenset(n for n in kept if kept[n] == lab) for lab in "0123"
    }


def test_composition_forwards_qc_failure_from_a_deeper_level(tmp_path: Path, caplog) -> None:
    """A stage-1 representative rejected by a later pass is not lost."""
    s1 = [
        DerepResult(
            representatives=[Path("/g/a"), Path("/g/c")],
            clusters={"a": ["b"], "c": []},
            genome_status={
                "a": STATUS_REPRESENTATIVE,
                "b": STATUS_CONTAINED,
                "c": "representative",
            },
        )
    ]
    s2 = DerepResult(
        representatives=[Path("/g/c")],
        clusters={"c": []},
        genome_status={"c": STATUS_REPRESENTATIVE, "a": STATUS_FAIL_QC},
    )
    log = logging.getLogger("test_chunk_recursion")
    with caplog.at_level(logging.WARNING, logger="test_chunk_recursion"):
        out = derep_stage._compose_two_stage(s1, s2, log)
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "1 representative(s)" in warnings[0] and "1 cluster member(s)" in warnings[0]
    assert out.genome_status == {"a": STATUS_FAIL_QC, "b": STATUS_FAIL_QC, "c": "representative"}
    check_result_complete(out, ["a", "b", "c"])


@pytest.fixture
def labelled_tool():
    registry._load()
    genomes = [f"Fam_g_s_GCA_{i:06d}.1.fasta" for i in range(40)]
    labels = {n: str(i % 5) for i, n in enumerate(genomes)}

    class _Registered(_Labelled):
        def __init__(self) -> None:
            super().__init__(labels, fail_names=frozenset({genomes[3]}))

    registry.register("labelled", _Registered, replace=True)
    yield genomes, labels
    registry._classes.pop("labelled", None)


def test_stage_run_logs_final_pass(workdir: Path, labelled_tool, caplog) -> None:
    genomes, _labels = labelled_tool
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for n in genomes:
        (gdir / n).write_text(">x\nACGT\n")
    ctx = WorkdirContext(workdir, create=True)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        res = run(ctx, DereplicateParams(tool="labelled", process_size=8))
    assert len(res.representatives) == 5
    assert res.genome_status[genomes[3]] == STATUS_FAIL_QC
    assert len(res.genome_status) == 40
    assert any(
        r.getMessage().startswith("Final pass over") and "bounds the leaf chunks" in r.getMessage()
        for r in caplog.records
    )


class _Recording(_Labelled):
    """Records (genome count, primary ANI, secondary ANI) for every call."""

    def __init__(self, labels: dict[str, str]):
        super().__init__(labels)
        self.thresholds: list[tuple[int, float, float]] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        self.thresholds.append((len(list(genomes)), params.primary_ani, params.secondary_ani))
        return super().dereplicate(genomes, out_dir, params, logger)


def test_levels_over_representatives_use_final_thresholds(tmp_path: Path) -> None:
    """Only the level-0 leaf chunks use --pre-*; later levels and the final pass use final."""
    genomes = _genomes(2000)
    labels = {g.name: str(i % 50) for i, g in enumerate(genomes)}
    adapter = _Recording(labels)
    final = _params()
    _dereplicate_chunked(
        adapter,
        genomes,
        tmp_path,
        final,
        100,
        2,
        0.85,
        0.95,
        logging.getLogger("test_chunk_recursion"),
    )
    level0, later = adapter.thresholds[:20], adapter.thresholds[20:]
    assert len(later) == 10 + 1  # one further level of 10 chunks, then the final pass
    assert all((p, s) == (0.85, 0.95) for _n, p, s in level0)
    assert all((p, s) == (final.primary_ani, final.secondary_ani) for _n, p, s in later)


def test_union_collapsing_into_one_chunk_logs_final_pass(tmp_path: Path, caplog) -> None:
    """A union of process_size + 1 becomes one chunk at the next level: still logged.

    110 genomes, chunk size 10, one label per leaf chunk: level 0 leaves 11
    representatives (99 removed), level 1 cuts 10 + 1 and merges the single
    genome into the first chunk, so that one call is the final pass.
    """
    genomes = _genomes(110)
    labels = {g.name: str(i // 10) for i, g in enumerate(genomes)}
    adapter = _Labelled(labels)
    with caplog.at_level(logging.INFO, logger="test_chunk_recursion"):
        res = _chunked(adapter, genomes, 10, tmp_path)
    assert adapter.calls == [10] * 11 + [11]
    finals = [r.getMessage() for r in caplog.records if r.getMessage().startswith("Final pass")]
    assert finals == [
        "Final pass over 11 representatives (chunk size 10 bounds the leaf chunks, not the merge)"
    ]
    assert len(res.representatives) == 11
    check_result_complete(res, labels)
