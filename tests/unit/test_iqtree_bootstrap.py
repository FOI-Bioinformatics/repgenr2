"""IQ-TREE's ultrafast bootstrap needs at least 1000 replicates; refuse fewer up front."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.treebuilders import iqtree as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test")


def _msa(tmp_path: Path) -> Path:
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nACGT\n>b\nACGA\n", encoding="utf-8")
    return msa


def test_bootstrap_below_1000_is_refused_before_running(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", lambda caps, cmd, **k: calls.append([str(c) for c in cmd]))
    with pytest.raises(UserInputError, match="1000"):
        mod.IqtreeBuilder().build(_msa(tmp_path), tmp_path / "out", TreeParams(bootstrap=100), _LOG)
    assert calls == []


def test_bootstrap_of_1000_is_passed_through(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        Path(cmd[cmd.index("-s") + 1] + ".treefile").write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    mod.IqtreeBuilder().build(_msa(tmp_path), tmp_path / "out", TreeParams(bootstrap=1000), _LOG)
    (cmd,) = calls
    assert cmd[cmd.index("-B") + 1] == "1000"


def test_sequence_type_is_stated(tmp_path: Path, monkeypatch) -> None:
    """A masked alignment is a fifth N; left to guess, IQ-TREE refuses it."""
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        Path(cmd[cmd.index("-s") + 1] + ".treefile").write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    mod.IqtreeBuilder().build(_msa(tmp_path), tmp_path / "out", TreeParams(), _LOG)
    assert calls[0][calls[0].index("-st") + 1] == "DNA"


def _run_capturing(tmp_path: Path, monkeypatch, params: TreeParams) -> list[str]:
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        Path(cmd[cmd.index("-s") + 1] + ".treefile").write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    mod.IqtreeBuilder().build(_msa(tmp_path), tmp_path / "out", params, _LOG)
    (cmd,) = calls
    return cmd


def test_model_extra_is_passed_as_m(tmp_path: Path, monkeypatch) -> None:
    cmd = _run_capturing(tmp_path, monkeypatch, TreeParams(extra={"model": "GTR+G"}))
    assert cmd[cmd.index("-m") + 1] == "GTR+G"
    assert "model" in mod.IqtreeBuilder.capabilities.accepted_extras


def test_without_model_modelfinder_is_left_to_iqtree_and_said(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    with caplog.at_level(logging.INFO, logger="test"):
        cmd = _run_capturing(tmp_path, monkeypatch, TreeParams())
    assert "-m" not in cmd
    assert any("ModelFinder" in r.getMessage() for r in caplog.records)


def test_constant_sites_become_fconst_and_invariant_columns_are_stripped(
    tmp_path: Path, monkeypatch
) -> None:
    # One invariant column (A in both) is removed from the alignment and
    # added to the A count; the variable column stays.
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n>b\nAG\n", encoding="utf-8")
    calls: list[list[str]] = []

    def fake(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        Path(cmd[cmd.index("-s") + 1] + ".treefile").write_text("(a,b);\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake)
    mod.IqtreeBuilder().build(
        msa, tmp_path / "out", TreeParams(constant_sites=(10, 20, 30, 40)), _LOG
    )
    (cmd,) = calls
    assert cmd[cmd.index("-fconst") + 1] == "11,20,30,40"
    assert "-m" not in cmd
    assert (tmp_path / "out" / "msa.fasta").read_text() == ">a\nC\n>b\nG\n"


def test_asc_lewis_uses_the_model_term_instead_of_counts(tmp_path: Path, monkeypatch) -> None:
    cmd = _run_capturing(
        tmp_path, monkeypatch, TreeParams(constant_sites=(1, 1, 1, 1), extra={"asc": "lewis"})
    )
    assert cmd[cmd.index("-m") + 1] == "MFP+ASC" and "-fconst" not in cmd
    cmd = _run_capturing(
        tmp_path, monkeypatch, TreeParams(extra={"asc": "lewis", "model": "GTR+G"})
    )
    assert cmd[cmd.index("-m") + 1] == "GTR+G+ASC"


def test_unknown_asc_value_is_refused(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: None)
    with pytest.raises(UserInputError, match="asc="):
        mod.IqtreeBuilder().build(
            _msa(tmp_path), tmp_path / "out", TreeParams(extra={"asc": "x"}), _LOG
        )
