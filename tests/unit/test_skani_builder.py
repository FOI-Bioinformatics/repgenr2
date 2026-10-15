"""skani tree builder: triangle call, matrix parsing, saturation, NJ hand-off."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.treebuilders import skani as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test-skani")
_FAR = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]


def _genomes(
    tmp_path: Path, names=("Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_1", "Out_Gen_sp_GCA_9")
) -> list[Path]:
    out = []
    for n in names:
        p = tmp_path / f"{n}.fasta"
        p.write_text(f">{n}\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


def _matrix_text(paths: list[Path], cells: list[list[float]]) -> str:
    lines = [str(len(paths))]
    for p, row in zip(paths, cells, strict=True):
        lines.append("\t".join([str(p), *(f"{x:.2f}" for x in row)]))
    return "\n".join(lines) + "\n"


def _fake(calls: list[list[str]], cells: list[list[float]]):
    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        listed = Path(cmd[cmd.index("-l") + 1]).read_text(encoding="utf-8").splitlines()
        paths = [Path(ln) for ln in listed]
        Path(cmd[cmd.index("-o") + 1]).write_text(_matrix_text(paths, cells), encoding="utf-8")
        return 0

    return run_tool


def test_build_runs_triangle_and_joins_in_name_order(tmp_path: Path, monkeypatch, caplog) -> None:
    calls: list[list[str]] = []
    # Input order 2, 1, 9; distances in percent; one undetected pair (1, 9),
    # which skani writes as 100.00 (checked with skani 0.4.0 and -s 95).
    cells = [[0, 1.5, 9.0], [1.5, 0, 100.0], [9.0, 100.0, 0]]
    monkeypatch.setattr(mod, "run_tool", _fake(calls, cells))
    seen: dict = {}

    def fake_nj(labels, dist, out_dir, logger, *, backend, threads=1):
        seen.update(labels=labels, dist=dist, backend=backend, threads=threads)
        return "(a,b,c);"

    monkeypatch.setattr(mod, "nj_tree", fake_nj)
    with caplog.at_level(logging.WARNING, logger="test-skani"):
        tree = mod.SkaniBuilder().build(
            _genomes(tmp_path),
            tmp_path / "out",
            TreeParams(threads=3, extra={"nj": "python"}),
            _LOG,
        )
    (cmd,) = calls
    assert cmd[:2] == ["skani", "triangle"] and cmd[cmd.index("-t") + 1] == "3"
    assert "--full-matrix" in cmd and "--distance" in cmd and "-s" not in cmd
    assert tree.read_text(encoding="utf-8") == "(a,b,c);\n"
    assert seen["labels"] == ["Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Out_Gen_sp_GCA_9"]
    assert seen["dist"][0][1] == pytest.approx(0.015)
    assert seen["dist"][1][2] == pytest.approx(0.09)
    assert seen["dist"][0][2] == mod.SATURATED and seen["backend"] == "python"
    assert seen["threads"] == 3
    assert any("no ANI for 1 of 3 pairs" in r.getMessage() for r in caplog.records)


def test_extras_reach_the_command_and_a_bad_preset_is_refused(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls, _FAR))
    monkeypatch.setattr(mod, "nj_tree", lambda *a, **k: "(a,b,c);")
    params = TreeParams(extra={"screen": "92.5", "min_af": 30, "preset": "slow"})
    mod.SkaniBuilder().build(_genomes(tmp_path), tmp_path / "out", params, _LOG)
    (cmd,) = calls
    assert cmd[cmd.index("-s") + 1] == "92.5" and cmd[cmd.index("--min-af") + 1] == "30"
    assert "--slow" in cmd
    with pytest.raises(UserInputError, match="preset"):
        mod.SkaniBuilder().build(
            _genomes(tmp_path), tmp_path / "out2", TreeParams(extra={"preset": "turbo"}), _LOG
        )


def test_distance_matrix_uses_mashtree_layout(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    cells = [[0, 2.0, 4.0], [2.0, 0, 6.0], [4.0, 6.0, 0]]
    monkeypatch.setattr(mod, "run_tool", _fake(calls, cells))
    matrix = mod.SkaniBuilder().distance_matrix(
        _genomes(tmp_path), tmp_path / "og", TreeParams(), _LOG
    )
    assert "--small-genomes" in calls[0]  # viral-sized inputs by default here
    lines = matrix.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    assert header == ["", "Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Out_Gen_sp_GCA_9"]
    row = dict(zip(header[1:], lines[1].split("\t")[1:], strict=True))
    assert lines[1].split("\t")[0] == "Fam_Gen_sp_GCA_1" and row["Out_Gen_sp_GCA_9"] == "0.060000"


def test_read_full_matrix_rejects_a_torn_file(tmp_path: Path) -> None:
    genomes = _genomes(tmp_path)
    bad = tmp_path / "m.tsv"
    bad.write_text("3\n" + str(genomes[0]) + "\t0\t1\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="3 x 3"):
        mod.read_full_matrix(bad, genomes)
    bad.write_text(_matrix_text([Path("/nowhere/x.fasta"), *genomes[1:]], _FAR), encoding="utf-8")
    with pytest.raises(WorkdirError, match="matches none"):
        mod.read_full_matrix(bad, genomes)


def test_rows_are_matched_by_path_or_record_name(tmp_path: Path) -> None:
    genomes = _genomes(tmp_path)
    # skani echoes the path as given in the list; a relative spelling still maps.
    rel = [Path(os.path.relpath(g, tmp_path)) for g in genomes]
    (tmp_path / "m.tsv").write_text(_matrix_text(rel, _FAR), encoding="utf-8")
    names, dist, saturated = mod.read_full_matrix(tmp_path / "m.tsv", genomes)
    assert names == sorted(names) and saturated == 0 and dist[0][1] == pytest.approx(0.01)


def test_skani_is_a_distance_matrix_builder() -> None:
    from repgenr.viral._outgroup import distance_matrix_builders

    assert "skani" in distance_matrix_builders()


def test_outgroup_matrix_is_withheld_when_no_candidate_is_within_reach(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    genomes = _genomes(tmp_path, ("S_Fam_Gen_sp_GCA_1", "S_Fam_Gen_sp_GCA_2", "O_Out_sp_GCA_9"))
    cells = [[0, 1.0, 100.0], [1.0, 0, 100.0], [100.0, 100.0, 0]]
    monkeypatch.setattr(mod, "run_tool", _fake([], cells))
    with caplog.at_level(logging.WARNING, logger="test-skani"):
        matrix = mod.SkaniBuilder().distance_matrix(genomes, tmp_path / "og", TreeParams(), _LOG)
    assert not matrix.exists()
    assert any(
        "no ANI between the selection and any outgroup" in r.getMessage() for r in caplog.records
    )
    # One candidate within reach is enough for a ranking.
    cells = [[0, 1.0, 100.0], [1.0, 0, 15.0], [100.0, 15.0, 0]]
    monkeypatch.setattr(mod, "run_tool", _fake([], cells))
    assert (
        mod.SkaniBuilder().distance_matrix(genomes, tmp_path / "og2", TreeParams(), _LOG).exists()
    )


def test_build_settles_the_joining_step_before_skani_runs(tmp_path: Path, monkeypatch) -> None:
    from repgenr.tree import nj

    monkeypatch.setattr(nj, "tool_available", lambda caps: False)
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: pytest.fail("skani ran"))
    genomes = [tmp_path / f"g{i}.fasta" for i in range(nj.PYTHON_NJ_MAX + 1)]
    with pytest.raises(WorkdirError, match="Neighbour joining"):
        mod.SkaniBuilder().build(genomes, tmp_path / "out", TreeParams(), _LOG)


def test_auto_prefers_mashtree_on_the_shared_limit() -> None:
    from repgenr.core.plugins import auto_select
    from repgenr.treebuilders.base import registry

    # Both declare 10000; the tie goes to the documented default order, then
    # alphabetical, so mashtree keeps winning auto for large alignment-free sets.
    assert auto_select(registry, 8000, allowed={"mashtree", "skani"}) == "mashtree"
