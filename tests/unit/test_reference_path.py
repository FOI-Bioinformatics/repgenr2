"""snptype --reference resolution must not traverse outside the genome dirs."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pytest

from repgenr.core.errors import UserInputError
from repgenr.stages.snptype import _reference_path


def _ctx(tmp_path: Path) -> SimpleNamespace:
    reps = tmp_path / "derep" / "representatives"
    genomes = tmp_path / "genomes"
    reps.mkdir(parents=True)
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">x\nACGT\n")
    return SimpleNamespace(representatives_dir=reps, genomes_dir=genomes)


_LOG = logging.getLogger("test")


def test_reference_resolves_by_basename(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    assert _reference_path(ctx, "g1.fasta", [], _LOG) == ctx.genomes_dir / "g1.fasta"


def test_reference_with_separator_rejected(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    (tmp_path / "secret.txt").write_text("x")
    with pytest.raises(UserInputError, match="basename"):
        _reference_path(ctx, "../secret.txt", [], _LOG)
    with pytest.raises(UserInputError, match="basename"):
        _reference_path(ctx, str(tmp_path / "secret.txt"), [], _LOG)


def test_default_reference_is_reported(tmp_path, caplog) -> None:
    """Without --reference the stage says which genome it picked (live audit)."""
    from pathlib import Path

    from repgenr.core.context import WorkdirContext

    ctx = WorkdirContext(tmp_path, create=True)
    genomes = [Path("b.fasta"), Path("c.fasta")]
    with caplog.at_level(logging.WARNING):
        assert _reference_path(ctx, None, genomes, logging.getLogger("t")) == Path("b.fasta")
    assert "alphabetically first genome" in caplog.text and "b.fasta" in caplog.text


def test_reference_resolves_in_the_active_genome_set_first(tmp_path: Path) -> None:
    """Under --all-genomes a representative is in genomes/ as well; the
    reference must be that genome, not the copy in representatives/, or the
    typer would get two inputs with one record name."""
    ctx = _ctx(tmp_path)
    (ctx.representatives_dir / "g1.fasta").write_text(">x\nACGT\n")
    active = [ctx.genomes_dir / "g1.fasta"]
    assert _reference_path(ctx, "g1.fasta", active, _LOG) == active[0]


def test_reference_reuses_the_genome_with_its_record_name(tmp_path: Path) -> None:
    """--reference g1.fasta when the set holds g1.fasta.gz is that genome."""
    ctx = _ctx(tmp_path)
    gz = ctx.genomes_dir / "g1.fasta.gz"
    gz.write_bytes(b"\x1f\x8b")
    (ctx.genomes_dir / "g1.fasta").unlink()
    (ctx.representatives_dir / "g1.fasta").write_text(">x\nACGT\n")
    assert _reference_path(ctx, "g1.fasta", [gz], _LOG) == gz


def test_reference_outside_the_active_set_is_still_found(tmp_path: Path) -> None:
    ctx = _ctx(tmp_path)
    (ctx.representatives_dir / "r1.fasta").write_text(">x\nACGT\n")
    active = [ctx.representatives_dir / "r1.fasta"]
    assert _reference_path(ctx, "g1.fasta", active, _LOG) == ctx.genomes_dir / "g1.fasta"
