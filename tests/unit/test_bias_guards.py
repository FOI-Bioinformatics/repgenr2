"""Bias guards from the scaling audit: reference logging, low-diversity warning."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.stages.phylo import _warn_low_diversity
from repgenr.tree.medoid import resolve_reference as _resolve_reference

_LOGGER = logging.getLogger("bias-guards")


def _genomes(tmp_path: Path, names: list[str]) -> list[Path]:
    out = []
    for name in names:
        p = tmp_path / name
        p.write_text(">x\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


# --- reference-choice logging -------------------------------------------------


def test_default_reference_fallback_warns(tmp_path, caplog):
    genomes = _genomes(tmp_path, ["a.fasta", "b.fasta"])
    with caplog.at_level(logging.WARNING):
        ref = _resolve_reference(None, genomes, None, _LOGGER).genome
    assert ref == genomes[0]
    assert any("alphabetically first" in r.message for r in caplog.records)


def test_explicit_reference_logs_info_not_warning(tmp_path, caplog):
    genomes = _genomes(tmp_path, ["a.fasta", "b.fasta"])
    with caplog.at_level(logging.INFO):
        ref = _resolve_reference("b.fasta", genomes, None, _LOGGER).genome
    assert ref.name == "b.fasta"
    assert not [r for r in caplog.records if r.levelno == logging.WARNING]


def test_unknown_reference_still_raises(tmp_path):
    genomes = _genomes(tmp_path, ["a.fasta"])
    with pytest.raises(UserInputError):
        _resolve_reference("missing.fasta", genomes, None, _LOGGER)


# --- low-diversity MSA warning ------------------------------------------------


def _write_msa(tmp_path: Path, seqs: list[str]) -> Path:
    msa = tmp_path / "msa.fasta"
    msa.write_text("".join(f">s{i}\n{seq}\n" for i, seq in enumerate(seqs)), encoding="utf-8")
    return msa


def test_invariant_msa_warns(tmp_path, caplog):
    msa = _write_msa(tmp_path, ["ACGTACGTAC"] * 5)
    with caplog.at_level(logging.WARNING):
        _warn_low_diversity(msa, _LOGGER)
    assert any("variable site" in r.message for r in caplog.records)


def test_near_invariant_msa_warns(tmp_path, caplog):
    base = "A" * 200
    variant = "C" + base[1:]  # a single variable column
    msa = _write_msa(tmp_path, [base, variant, base])
    with caplog.at_level(logging.WARNING):
        _warn_low_diversity(msa, _LOGGER)
    assert any("1 variable site" in r.message for r in caplog.records)


def test_diverse_msa_stays_quiet(tmp_path, caplog):
    base = "A" * 100
    variant = "C" * 20 + base[20:]  # 20 variable columns
    msa = _write_msa(tmp_path, [base, variant])
    with caplog.at_level(logging.WARNING):
        _warn_low_diversity(msa, _LOGGER)
    assert not caplog.records


def test_missing_msa_never_raises(tmp_path, caplog):
    _warn_low_diversity(tmp_path / "absent.fasta", _LOGGER)  # must not raise


# --- median_of_medians is documented ------------------------------------------


def test_length_method_help_explains_the_defense():
    import inspect

    from repgenr.cli import cmd_viral

    source = inspect.getsource(cmd_viral)
    assert "one vote per species" in source
