"""VeryFastTree builder: the FastTree command family plus -threads and two extras."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.treebuilders import fasttree as ft
from repgenr.treebuilders import veryfasttree as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test")


def _msa(tmp_path: Path) -> Path:
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    return msa


def _run(tmp_path: Path, monkeypatch, params: TreeParams) -> tuple[list[str], dict]:
    calls: list[list[str]] = []
    seen: dict = {}

    def fake(caps, cmd, **k):
        calls.append([str(c) for c in cmd])
        seen["env"] = dict(k.get("tool_env") or {})
        seen["stdout_path"] = k["stdout_path"]
        Path(k["stdout_path"]).write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    tree = mod.VeryFastTreeBuilder().build(_msa(tmp_path), tmp_path / "out", params, _LOG)
    (cmd,) = calls
    assert tree == Path(seen["stdout_path"])
    return cmd, seen


def test_default_command_threads_and_omp(tmp_path: Path, monkeypatch) -> None:
    cmd, seen = _run(tmp_path, monkeypatch, TreeParams(threads=6))
    assert cmd[:3] == ["VeryFastTree", "-nt", "-gtr"]
    assert cmd[cmd.index("-threads") + 1] == "6"
    assert seen["env"]["OMP_NUM_THREADS"] == "6"
    assert "-boot" not in cmd and "-gamma" not in cmd and "-double-precision" not in cmd
    assert cmd[-1].endswith("msa.fasta")


def test_bootstrap_and_extras_reach_the_command(tmp_path: Path, monkeypatch) -> None:
    cmd, _ = _run(
        tmp_path,
        monkeypatch,
        TreeParams(bootstrap=200, extra={"gamma": "true", "double_precision": True}),
    )
    assert cmd[cmd.index("-boot") + 1] == "200"
    assert "-gamma" in cmd and "-double-precision" in cmd


def test_a_bad_boolean_extra_is_a_user_error(tmp_path: Path, monkeypatch) -> None:
    with pytest.raises(UserInputError, match="gamma"):
        _run(tmp_path, monkeypatch, TreeParams(extra={"gamma": "maybe"}))


def test_constant_sites_are_noted_not_applied(tmp_path: Path, monkeypatch, caplog) -> None:
    with caplog.at_level(logging.INFO, logger="test"):
        cmd, _ = _run(tmp_path, monkeypatch, TreeParams(constant_sites=(1, 2, 3, 4)))
    assert "-fconst" not in cmd
    assert any("no ascertainment-bias correction" in r.getMessage() for r in caplog.records)


def test_fasttree_takes_gamma_too(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        calls.append([str(c) for c in cmd])
        Path(k["stdout_path"]).write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(ft, "run_tool", fake)
    ft.FasttreeBuilder().build(_msa(tmp_path), tmp_path / "o", TreeParams(extra={"gamma": 1}), _LOG)
    assert "-gamma" in calls[0]


def test_parse_extra_bool_semantics() -> None:
    from repgenr.core.plugins import parse_extra_bool

    assert parse_extra_bool({}, "k", True) is True and parse_extra_bool({}, "k", False) is False
    assert parse_extra_bool({"k": ""}, "k", True) is True
    assert parse_extra_bool({"k": True}, "k", False) is True
    assert parse_extra_bool({"k": "no"}, "k", True) is False
    assert parse_extra_bool({"k": "YES"}, "k", False) is True
    with pytest.raises(UserInputError, match="'k'"):
        parse_extra_bool({"k": "maybe"}, "k", False)


def test_auto_keeps_alignment_free_builders_past_the_aligner_limit(monkeypatch) -> None:
    """Without an alignment request, auto must not pick an MSA builder for a
    set the default aligner cannot align; with one, veryfasttree is picked
    above fasttree's limit."""
    from repgenr.stages.phylo import PhyloParams, _resolve_treebuilder

    monkeypatch.setattr("repgenr.core.plugins._tool_available", lambda caps: True)
    log = logging.getLogger("test")
    auto = PhyloParams(treebuilder="auto")
    snp = PhyloParams(treebuilder="auto", msa_source="snptype")
    assert _resolve_treebuilder(auto, 12000, log) == "mashtree"
    assert _resolve_treebuilder(auto, 4000, log) == "mashtree"  # aligner limit 500
    assert _resolve_treebuilder(auto, 300, log) == "iqtree"
    assert _resolve_treebuilder(snp, 8000, log) == "veryfasttree"
    assert _resolve_treebuilder(snp, 4000, log) == "fasttree"
