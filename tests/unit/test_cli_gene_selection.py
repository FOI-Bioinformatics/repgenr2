"""The gene-selection options and --partitioned of the core-gene route on phylo,
phylo-build and run."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from repgenr.cli import cmd_phylo, cmd_run
from repgenr.cli.main import app
from repgenr.stages import phylo as phylo_stage
from repgenr.stages.phylo import PhyloParams

_runner = CliRunner()
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_FLAGS = (
    ("--max-core-genes", "5", "max_core_genes", 5),
    ("--gene-rank", "informative", "gene_rank", "informative"),
    ("--min-gene-length", "300", "min_gene_length", 300),
    ("--max-gene-gaps", "0.25", "max_gene_gaps", 0.25),
)


def _text(result) -> str:
    return " ".join(_ANSI.sub("", result.output + str(result.exception or "")).split())


@pytest.fixture
def built(monkeypatch) -> list[Any]:
    """The params objects the commands build instead of running a stage."""
    objects: list[Any] = []

    def fake_run(stage, workdir, build, *, create=False):
        objects.append(build())

    monkeypatch.setattr(cmd_phylo, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    return objects


def _phylo(objects: list[Any]) -> PhyloParams:
    return [o for o in objects if isinstance(o, PhyloParams)][-1]


def _phylo_args(tmp_path: Path) -> list[str]:
    return ["phylo", "-wd", str(tmp_path / "wd")]


def _run_args(tmp_path: Path) -> list[str]:
    return ["run", "-wd", str(tmp_path / "wd"), "-l", "genus"]


def _build_args(tmp_path: Path) -> list[str]:
    genomes = tmp_path / "genomes"
    genomes.mkdir(exist_ok=True)
    (genomes / "a.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    return ["phylo-build", "--genomes-dir", str(genomes), "-o", str(tmp_path / "out")]


COMMANDS = {"phylo": _phylo_args, "phylo-build": _build_args, "run": _run_args}


@pytest.mark.parametrize("command", sorted(COMMANDS))
@pytest.mark.parametrize(("flag", "value", "key", "expected"), _FLAGS)
def test_selection_flag_is_refused_outside_the_coregene_source(
    built, tmp_path, command, flag, value, key, expected
) -> None:
    result = _runner.invoke(app, [*COMMANDS[command](tmp_path), flag, value])
    assert result.exit_code != 0
    text = _text(result)
    assert flag in text and "only with --msa-source coregene" in text


@pytest.mark.parametrize("command", ["phylo", "run"])
def test_selection_flags_reach_the_phylo_extra(built, tmp_path, command) -> None:
    argv = [*COMMANDS[command](tmp_path), "--msa-source", "coregene"]
    for flag, value, _key, _expected in _FLAGS:
        argv += [flag, value]
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, _text(result)
    extra = _phylo(built).extra
    for _flag, _value, key, expected in _FLAGS:
        assert extra[key] == expected


@pytest.mark.parametrize("command", ["phylo", "run"])
def test_unset_selection_adds_nothing_to_the_extra(built, tmp_path, command) -> None:
    result = _runner.invoke(app, [*COMMANDS[command](tmp_path), "--msa-source", "coregene"])
    assert result.exit_code == 0, _text(result)
    assert not {k for _f, _v, k, _e in _FLAGS} & _phylo(built).extra.keys()


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_unknown_gene_rank_is_refused(built, tmp_path, command) -> None:
    argv = [*COMMANDS[command](tmp_path), "--msa-source", "coregene", "--gene-rank", "turbo"]
    result = _runner.invoke(app, argv)
    assert result.exit_code != 0
    assert "--gene-rank" in _text(result) and "occupancy" in _text(result)


@pytest.mark.parametrize("command", sorted(COMMANDS))
@pytest.mark.parametrize(
    ("flag", "bad"),
    [("--max-core-genes", "-1"), ("--min-gene-length", "-1"), ("--max-gene-gaps", "1.5")],
)
def test_out_of_range_values_are_refused(built, tmp_path, command, flag, bad) -> None:
    argv = [*COMMANDS[command](tmp_path), "--msa-source", "coregene", flag, bad]
    result = _runner.invoke(app, argv)
    assert result.exit_code != 0 and flag in _text(result)


def test_phylo_build_passes_the_selection_to_the_stage(monkeypatch, tmp_path) -> None:
    seen: list[Any] = []
    monkeypatch.setattr(phylo_stage, "phylo_build", lambda params, *a, **k: seen.append(params))
    argv = [*_build_args(tmp_path), "--msa-source", "coregene", "--max-core-genes", "7"]
    argv += ["--gene-rank", "random", "--min-gene-length", "90", "--max-gene-gaps", "0.1"]
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, _text(result)
    assert seen[0].phylo.extra["max_core_genes"] == 7
    assert seen[0].phylo.extra["gene_rank"] == "random"
    assert seen[0].phylo.extra["min_gene_length"] == 90
    assert seen[0].phylo.extra["max_gene_gaps"] == 0.1


# --partitioned: a tree-side option of the core-gene route.


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_partitioned_is_refused_outside_the_coregene_source(built, tmp_path, command) -> None:
    result = _runner.invoke(app, [*COMMANDS[command](tmp_path), "--partitioned"])
    assert result.exit_code != 0
    text = _text(result)
    assert "--partitioned" in text and "only with --msa-source coregene" in text


@pytest.mark.parametrize("command", ["phylo", "run"])
def test_partitioned_reaches_the_phylo_extra(built, tmp_path, command) -> None:
    argv = [*COMMANDS[command](tmp_path), "--msa-source", "coregene", "--partitioned"]
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, _text(result)
    assert _phylo(built).extra["partitioned"] is True


@pytest.mark.parametrize("command", ["phylo", "run"])
def test_unset_partitioned_adds_nothing(built, tmp_path, command) -> None:
    result = _runner.invoke(app, [*COMMANDS[command](tmp_path), "--msa-source", "coregene"])
    assert result.exit_code == 0, _text(result)
    assert "partitioned" not in _phylo(built).extra


def test_phylo_build_with_msa_leaves_the_decision_to_the_sidecar(monkeypatch, tmp_path) -> None:
    seen: list[Any] = []
    monkeypatch.setattr(phylo_stage, "phylo_build", lambda params, *a, **k: seen.append(params))
    msa = tmp_path / "msa.fasta"
    msa.write_text(">a\nACGT\n", encoding="utf-8")
    # The default source is not coregene; with --msa the flag is accepted.
    result = _runner.invoke(app, [*_build_args(tmp_path), "--msa", str(msa), "--partitioned"])
    assert result.exit_code == 0, _text(result)
    assert seen[0].phylo.extra["partitioned"] is True
    assert seen[0].msa == msa


def test_phylo_build_coregene_passes_partitioned(monkeypatch, tmp_path) -> None:
    seen: list[Any] = []
    monkeypatch.setattr(phylo_stage, "phylo_build", lambda params, *a, **k: seen.append(params))
    argv = [*_build_args(tmp_path), "--msa-source", "coregene", "--partitioned"]
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, _text(result)
    assert seen[0].phylo.extra["partitioned"] is True
