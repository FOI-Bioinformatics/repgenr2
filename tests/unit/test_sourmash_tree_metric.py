"""sourmash tree builder: distance 1 - ANI by default, 1 - Jaccard on request."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.treebuilders import sourmash as sm
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test-sourmash-tree")
# a and b share most hashes, c shares none with either: the Jaccard matrix,
# and the average containment of the same sketches (equal sizes, 2j/(1+j)).
_CSV = "a.fasta,b.fasta,c.fasta\n1,0.796,0\n0.796,1,0\n0,0,1\n"
_AVG_C = 2 * 0.796 / 1.796
_CSV_AVG = f"a.fasta,b.fasta,c.fasta\n1,{_AVG_C},0\n{_AVG_C},1,0\n0,0,1\n"


def _run(tmp_path: Path, monkeypatch, extra: dict) -> tuple[list[list[str]], dict]:
    calls: list[list[str]] = []
    seen: dict = {}

    def fake_run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        if cmd[:2] == ["sourmash", "sketch"]:
            outdir = Path(cmd[cmd.index("--outdir") + 1])
            fofn = Path(cmd[cmd.index("--from-file") + 1])
            for line in fofn.read_text(encoding="utf-8").splitlines():
                (outdir / (Path(line).name + ".sig")).write_text("{}", encoding="utf-8")
        elif cmd[:2] == ["sourmash", "compare"]:
            body = _CSV_AVG if "--avg-containment" in cmd else _CSV
            Path(cmd[cmd.index("--csv") + 1]).write_text(body, encoding="utf-8")
        return 0

    def fake_nj(labels, dist, out_dir, logger, *, backend, threads=1):
        seen["labels"], seen["dist"] = labels, dist
        return "(a,b,c);"

    monkeypatch.setattr(sm, "run_tool", fake_run_tool)
    monkeypatch.setattr(sm, "nj_tree", fake_nj)
    genomes = []
    for n in "abc":
        p = tmp_path / f"{n}.fasta"
        p.write_text(f">{n}\nACGT\n", encoding="utf-8")
        genomes.append(p)
    out = tmp_path / ("out_" + "_".join(f"{k}{v}" for k, v in extra.items()))
    sm.SourmashBuilder().build(genomes, out, TreeParams(extra={"nj": "python", **extra}), _LOG)
    return calls, seen


def test_default_distance_is_one_minus_ani(tmp_path: Path, monkeypatch, caplog) -> None:
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        calls, seen = _run(tmp_path, monkeypatch, {})
    compare = next(c for c in calls if c[:2] == ["sourmash", "compare"])
    # The ANI is estimated from the average containment, not by compare
    # --ani, which reports 0 for small sketches of identical genomes.
    assert "--ani" not in compare and "--avg-containment" in compare
    dist = seen["dist"]
    # Average containment 0.8864 at k=31 is an ANI estimate of 0.9961.
    assert abs(dist[0][1] - (1 - 0.9961)) < 5e-4
    assert dist[0][0] == 0.0 and dist[0][2] == 1.0 and dist[1][2] == 1.0
    (line,) = [r.getMessage() for r in caplog.records if "share no hash" in r.getMessage()]
    assert line.startswith("sourmash: 2 of 3 pairs")
    assert "1 - ANI" in caplog.text


def test_the_conversion_follows_the_builder_ksize(tmp_path: Path, monkeypatch) -> None:
    """At k=21 the same containment is a lower ANI than at k=31: c^(1/k)."""
    _calls, seen = _run(tmp_path, monkeypatch, {"ksize": 21})
    expected = 1 - _AVG_C ** (1 / 21)
    assert seen["dist"][0][1] == pytest.approx(expected, abs=1e-12)
    assert seen["dist"][0][1] > 1 - 0.9961 + 1e-3  # not the k=31 value


def test_metric_jaccard_keeps_one_minus_jaccard(tmp_path: Path, monkeypatch) -> None:
    calls, seen = _run(tmp_path, monkeypatch, {"metric": "jaccard"})
    compare = next(c for c in calls if c[:2] == ["sourmash", "compare"])
    assert "--avg-containment" not in compare
    assert abs(seen["dist"][0][1] - (1 - 0.796)) < 1e-9
    assert seen["dist"][0][2] == 1.0


def test_an_unknown_metric_is_refused_before_any_work(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(sm, "run_tool", lambda *a, **k: pytest.fail("sourmash ran"))
    with pytest.raises(UserInputError, match="metric"):
        sm.SourmashBuilder().build(
            [tmp_path / "a.fasta"], tmp_path / "out", TreeParams(extra={"metric": "mash"}), _LOG
        )


def test_the_tree_builder_shares_the_dereplicator_converter() -> None:
    import numpy as np

    from repgenr.core.sourmash import containment_to_ani_matrix
    from repgenr.dereplicators import sourmash as derep

    assert derep._containment_to_ani_matrix is containment_to_ani_matrix
    assert sm.containment_to_ani_matrix is containment_to_ani_matrix
    for ksize in (21, 31):
        ani, _unknown = containment_to_ani_matrix(np.array([[1.0, 0.8864, 0.0]]), ksize)
        assert ani[0, 0] == 1.0 and ani[0, 2] == 0.0  # identical stays 1.0, disjoint 0
        assert ani[0, 1] == pytest.approx(0.8864 ** (1 / ksize))
    assert "metric" in sm.SourmashBuilder.capabilities.accepted_extras
    assert sm.SourmashBuilder.capabilities.default_params["metric"] == "ani"
