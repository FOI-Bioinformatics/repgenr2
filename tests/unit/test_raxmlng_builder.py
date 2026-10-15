"""RAxML-NG builder: publish the support-annotated tree when --all produced one."""

from __future__ import annotations

import logging
from pathlib import Path

import repgenr.treebuilders.raxmlng as mod
from repgenr.treebuilders.base import TreeParams

_BEST = "((a,b),(c,d));"
_SUPPORT = "((a,b)97,(c,d)100);"


def _fake_run(outputs: dict[str, str]):
    def run_tool(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        prefix = cmd[cmd.index("--prefix") + 1]
        for suffix, text in outputs.items():
            Path(prefix + suffix).write_text(text, encoding="utf-8")
        return 0

    return run_tool


def test_publishes_support_tree_when_present(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        mod, "run_tool", _fake_run({".raxml.bestTree": _BEST, ".raxml.support": _SUPPORT})
    )
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n", encoding="utf-8")
    tree = mod.RaxmlNgBuilder().build(msa, tmp_path / "out", TreeParams(), logging.getLogger())
    assert tree.read_text(encoding="utf-8") == _SUPPORT


def test_falls_back_to_best_tree_without_support(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake_run({".raxml.bestTree": _BEST}))
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n", encoding="utf-8")
    tree = mod.RaxmlNgBuilder().build(msa, tmp_path / "out", TreeParams(), logging.getLogger())
    assert tree.read_text(encoding="utf-8") == _BEST


def _model_of(monkeypatch, tmp_path: Path, params: TreeParams) -> tuple[str, str]:
    seen: dict = {}

    def run_tool(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        seen["model"] = cmd[cmd.index("--model") + 1]
        seen["msa"] = cmd[cmd.index("--msa") + 1]
        Path(cmd[cmd.index("--prefix") + 1] + ".raxml.bestTree").write_text(_BEST, encoding="utf-8")
        return 0

    monkeypatch.setattr(mod, "run_tool", run_tool)
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n>b\nAG\n", encoding="utf-8")
    mod.RaxmlNgBuilder().build(msa, tmp_path / "out", params, logging.getLogger())
    return seen["model"], seen["msa"]


def test_constant_sites_become_asc_stam_over_the_variable_columns(monkeypatch, tmp_path) -> None:
    model, msa = _model_of(monkeypatch, tmp_path, TreeParams(constant_sites=(10, 20, 30, 40)))
    assert model == "GTR+G+ASC_STAM{11/20/30/40}"
    assert Path(msa).name == "msa.variable.fasta"
    assert Path(msa).read_text() == ">a\nC\n>b\nG\n"


def test_asc_lewis_appends_the_lewis_term(monkeypatch, tmp_path) -> None:
    model, _ = _model_of(
        monkeypatch, tmp_path, TreeParams(constant_sites=(1, 1, 1, 1), extra={"asc": "lewis"})
    )
    assert model == "GTR+G+ASC_LEWIS"


def test_without_counts_the_model_is_unchanged(monkeypatch, tmp_path) -> None:
    model, msa = _model_of(monkeypatch, tmp_path, TreeParams(extra={"model": "GTR"}))
    assert model == "GTR" and Path(msa).name == "msa.fasta"
