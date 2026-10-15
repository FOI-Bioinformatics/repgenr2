"""Neighbour joining through rapidNJ or the Python code (repgenr.tree.nj)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.tree import nj as mod

_LOG = logging.getLogger("test-nj")
_DIST = [
    [0.0, 0.1, 0.9, 0.9],
    [0.1, 0.0, 0.9, 0.9],
    [0.9, 0.9, 0.0, 0.1],
    [0.9, 0.9, 0.1, 0.0],
]


def test_auto_prefers_rapidnj_when_available(monkeypatch, caplog) -> None:
    monkeypatch.setattr(mod, "tool_available", lambda caps: True)
    assert mod.choose_backend(mod.NJ_AUTO, 10, _LOG) == mod.NJ_RAPIDNJ
    monkeypatch.setattr(mod, "tool_available", lambda caps: False)
    with caplog.at_level(logging.INFO, logger="test-nj"):
        assert mod.choose_backend(mod.NJ_AUTO, 10, _LOG) == mod.NJ_PYTHON
    assert any("rapidnj is not installed" in r.getMessage() for r in caplog.records)


def test_explicit_backends_and_the_python_ceiling(monkeypatch) -> None:
    monkeypatch.setattr(mod, "tool_available", lambda caps: False)
    assert mod.choose_backend(mod.NJ_RAPIDNJ, 10, _LOG) == mod.NJ_RAPIDNJ  # preflight decides
    assert mod.choose_backend(mod.NJ_PYTHON, 10, _LOG) == mod.NJ_PYTHON
    with pytest.raises(WorkdirError, match="install rapidnj"):
        mod.choose_backend(mod.NJ_PYTHON, mod.PYTHON_NJ_MAX + 1, _LOG)
    monkeypatch.setattr(mod, "tool_available", lambda caps: True)
    with pytest.raises(WorkdirError, match="drop nj=python"):
        mod.choose_backend(mod.NJ_PYTHON, mod.PYTHON_NJ_MAX + 1, _LOG)
    with pytest.raises(UserInputError, match="nj must be one of"):
        mod.choose_backend("fast", 10, _LOG)


def test_write_phylip_uses_stand_in_names_and_clamps(tmp_path: Path) -> None:
    mod.write_phylip(tmp_path / "m.phylip", 2, [[0.0, -0.1], [-0.1, 0.0]])
    assert (tmp_path / "m.phylip").read_text(encoding="utf-8").splitlines() == [
        "2",
        "0 0.000000 0.000000",
        "1 0.000000 0.000000",
    ]


def test_restore_labels_handles_bare_and_quoted_names() -> None:
    labels = ["Fam_Gen_sp_GCA_1", "a b", "c:d", "e"]
    out = mod.restore_labels("(('1':0.05,0:0.05):0.8,'2':0.05,3:0.05);", labels)
    assert out == "(('a b':0.05,Fam_Gen_sp_GCA_1:0.05):0.8,'c:d':0.05,e:0.05);"
    # Branch lengths that look like integers are not names: only tokens before ':'.
    assert mod.restore_labels("(0:1,1:2);", ["x", "y"]) == "(x:1,y:2);"
    # A support value after ')' is not a leaf either.
    assert mod.restore_labels("((0:1,1:2)95:0.5,2:1);", ["x", "y", "z"]) == "((x:1,y:2)95:0.5,z:1);"
    with pytest.raises(WorkdirError, match="missing \\[1\\]"):
        mod.restore_labels("(0:1,0:2);", ["x", "y"])
    with pytest.raises(WorkdirError, match="unknown leaf"):
        mod.restore_labels("(0:1,7:2);", ["x", "y"])


def test_nj_tree_runs_rapidnj_and_restores_labels(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(mod, "tool_available", lambda caps: True)
    seen: dict = {}

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        seen["cmd"] = cmd
        seen["matrix"] = Path(cmd[1]).read_text(encoding="utf-8")
        Path(cmd[cmd.index("-x") + 1]).write_text(
            "(('1':0.05,'0':0.05):0.8,'2':0.05,'3':0.05);\n", encoding="utf-8"
        )
        return 0

    monkeypatch.setattr(mod, "run_tool", fake)
    labels = ["g1", "g 2", "g3", "g4"]
    newick = mod.nj_tree(labels, _DIST, tmp_path, _LOG)
    assert seen["cmd"][:1] == ["rapidnj"] and seen["cmd"][2:5] == ["-i", "pd", "-n"]
    assert seen["cmd"][seen["cmd"].index("-c") + 1] == "1"
    assert not (tmp_path / "distances.phylip").exists()
    assert not (tmp_path / "rapidnj.nwk").exists()
    assert seen["matrix"].startswith("4\n0 0.000000 0.100000")
    assert newick == "(('g 2':0.05,g1:0.05):0.8,g3:0.05,g4:0.05);"


def test_nj_tree_python_path_and_small_inputs(monkeypatch, tmp_path: Path) -> None:
    import dendropy

    monkeypatch.setattr(mod, "tool_available", lambda caps: False)
    newick = mod.nj_tree(["A", "B", "C", "D"], _DIST, tmp_path, _LOG)
    tree = dendropy.Tree.get(data=newick, schema="newick")
    assert sorted(t.label for t in tree.taxon_namespace) == ["A", "B", "C", "D"]
    # Fewer than three leaves never call rapidnj, whatever is installed.
    monkeypatch.setattr(mod, "tool_available", lambda caps: True)
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: pytest.fail("rapidnj called"))
    assert mod.nj_tree(["A", "B"], [[0.0, 0.4], [0.4, 0.0]], tmp_path, _LOG).endswith(");")


def test_nj_tree_rejects_a_mismatched_matrix(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="square"):
        mod.nj_tree(["A", "B"], [[0.0]], tmp_path, _LOG)


def test_rapidnj_without_a_tree_is_an_error(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(mod, "tool_available", lambda caps: True)
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: 0)
    with pytest.raises(WorkdirError, match="wrote no tree"):
        mod.nj_tree(["A", "B", "C"], [[0, 1, 1], [1, 0, 1], [1, 1, 0]], tmp_path, _LOG)


def test_sourmash_builder_passes_the_nj_extra(tmp_path: Path, monkeypatch) -> None:
    from repgenr.treebuilders import sourmash as sm
    from repgenr.treebuilders.base import TreeParams

    seen: dict = {}

    def fake_run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        if cmd[:2] == ["sourmash", "sketch"]:
            outdir = Path(cmd[cmd.index("--outdir") + 1])
            for line in (
                Path(cmd[cmd.index("--from-file") + 1]).read_text(encoding="utf-8").splitlines()
            ):
                (outdir / (Path(line).name + ".sig")).write_text("{}", encoding="utf-8")
        elif cmd[:2] == ["sourmash", "compare"]:
            Path(cmd[cmd.index("--csv") + 1]).write_text(
                "a.fasta,b.fasta,c.fasta\n1,0.9,0.2\n0.9,1,0.2\n0.2,0.2,1\n",
                encoding="utf-8",
            )
        return 0

    def fake_nj(labels, dist, out_dir, logger, *, backend, threads=1):
        seen["backend"] = backend
        seen["labels"] = labels
        return "(a,b,c);"

    monkeypatch.setattr(sm, "run_tool", fake_run_tool)
    monkeypatch.setattr(sm, "nj_tree", fake_nj)
    genomes = []
    for n in "abc":
        p = tmp_path / f"{n}.fasta"
        p.write_text(f">{n}\nACGT\n", encoding="utf-8")
        genomes.append(p)
    sm.SourmashBuilder().build(genomes, tmp_path / "out", TreeParams(extra={"nj": "python"}), _LOG)
    assert seen["backend"] == "python" and seen["labels"] == ["a", "b", "c"]


def test_sourmash_builder_settles_the_backend_before_any_work(tmp_path: Path, monkeypatch) -> None:
    from repgenr.treebuilders import sourmash as sm
    from repgenr.treebuilders.base import TreeParams

    monkeypatch.setattr(mod, "tool_available", lambda caps: False)
    monkeypatch.setattr(sm, "run_tool", lambda *a, **k: pytest.fail("sourmash ran"))
    genomes = [tmp_path / f"g{i}.fasta" for i in range(mod.PYTHON_NJ_MAX + 1)]
    with pytest.raises(WorkdirError, match="Neighbour joining"):
        sm.SourmashBuilder().build(genomes, tmp_path / "out", TreeParams(), _LOG)
    with pytest.raises(UserInputError, match="nj must be one of"):
        sm.SourmashBuilder().build(
            genomes[:3], tmp_path / "out", TreeParams(extra={"nj": "fast"}), _LOG
        )


def test_rapidnj_is_a_known_bin_dir_tool() -> None:
    from repgenr.core.bindirs import known_capabilities

    assert "rapidnj" in known_capabilities()
