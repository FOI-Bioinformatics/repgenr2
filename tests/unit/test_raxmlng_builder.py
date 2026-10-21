"""RAxML-NG builder: publish the support-annotated tree when --all produced one."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

import repgenr.treebuilders.raxmlng as mod
from repgenr.core.errors import UserInputError
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
    assert Path(msa).read_text(encoding="utf-8") == ">a\nC\n>b\nG\n"


def test_asc_lewis_appends_the_lewis_term(monkeypatch, tmp_path) -> None:
    model, _ = _model_of(
        monkeypatch, tmp_path, TreeParams(constant_sites=(1, 1, 1, 1), extra={"asc": "lewis"})
    )
    assert model == "GTR+G+ASC_LEWIS"


def test_without_counts_the_model_is_unchanged(monkeypatch, tmp_path) -> None:
    model, msa = _model_of(monkeypatch, tmp_path, TreeParams(extra={"model": "GTR"}))
    assert model == "GTR" and Path(msa).name == "msa.fasta"


def test_a_model_with_its_own_asc_term_is_passed_unchanged(monkeypatch, tmp_path) -> None:
    model, msa = _model_of(
        monkeypatch,
        tmp_path,
        TreeParams(constant_sites=(1, 1, 1, 1), extra={"model": "GTR+ASC_LEWIS"}),
    )
    assert model == "GTR+ASC_LEWIS" and Path(msa).name == "msa.variable.fasta"


def test_unknown_asc_value_is_refused(monkeypatch, tmp_path) -> None:
    import pytest

    from repgenr.core.errors import UserInputError

    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n>b\nAG\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="asc="):
        mod.RaxmlNgBuilder().build(
            msa, tmp_path / "out", TreeParams(extra={"asc": "x"}), logging.getLogger()
        )


_PARTS = (("gyrA", 1, 1), ("gyr.A", 2, 2))


def test_partitions_replace_the_model_with_a_file(monkeypatch, tmp_path) -> None:
    model, _ = _model_of(monkeypatch, tmp_path, TreeParams(partitions=_PARTS))
    assert Path(model) == tmp_path / "out" / "partitions.txt"
    assert Path(model).read_text(encoding="utf-8") == "GTR+G, gyrA = 1-1\nGTR+G, gyr_A = 2-2\n"


def test_partitions_carry_the_given_model(monkeypatch, tmp_path) -> None:
    model, msa = _model_of(
        monkeypatch, tmp_path, TreeParams(partitions=_PARTS, extra={"model": "HKY+G"})
    )
    assert Path(model).read_text(encoding="utf-8") == "HKY+G, gyrA = 1-1\nHKY+G, gyr_A = 2-2\n"
    # The alignment is used as given: no invariant-column strip.
    assert Path(msa).name == "msa.fasta"


def test_protein_partitions_default_to_lg(monkeypatch, tmp_path) -> None:
    model, _ = _model_of(
        monkeypatch, tmp_path, TreeParams(partitions=(("g", 1, 2),), alphabet="protein")
    )
    assert Path(model).read_text(encoding="utf-8") == "LG+G, g = 1-2\n"


@pytest.mark.parametrize(
    "params",
    [
        TreeParams(partitions=_PARTS, constant_sites=(1, 1, 1, 1)),
        TreeParams(partitions=_PARTS, extra={"asc": "lewis"}),
        TreeParams(partitions=_PARTS, extra={"model": "GTR+G+ASC_LEWIS"}),
    ],
    ids=["constant_sites", "asc_lewis", "asc_model"],
)
def test_partitions_with_asc_are_refused(monkeypatch, tmp_path, params) -> None:
    calls: list = []
    monkeypatch.setattr(mod, "run_tool", lambda caps, cmd, **k: calls.append(cmd))
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nAC\n>b\nAG\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="partition"):
        mod.RaxmlNgBuilder().build(msa, tmp_path / "out", params, logging.getLogger())
    assert calls == []
    assert mod.RaxmlNgBuilder.supports_partitions


def test_a_protein_alignment_defaults_to_lg(monkeypatch, tmp_path) -> None:
    model, _ = _model_of(monkeypatch, tmp_path, TreeParams(alphabet="protein"))
    assert model == "LG+G"
    (tmp_path / "dna").mkdir()
    model, _ = _model_of(monkeypatch, tmp_path / "dna", TreeParams())
    assert model == "GTR+G"


@pytest.mark.parametrize("extra", [{"asc": "lewis"}, {"model": "LG+ASC_LEWIS"}])
def test_asc_on_a_protein_alignment_is_refused(monkeypatch, tmp_path, extra) -> None:
    calls: list = []
    monkeypatch.setattr(mod, "run_tool", lambda caps, cmd, **k: calls.append(cmd))
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nMKLV\n>b\nMKIV\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="nucleotide SNP alignments"):
        mod.RaxmlNgBuilder().build(
            msa, tmp_path / "out", TreeParams(alphabet="protein", extra=extra), logging.getLogger()
        )
    assert calls == []
