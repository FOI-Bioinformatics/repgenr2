"""Constant-site counting and invariant-column stripping for the tree builders."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.tree.alignment import (
    ALIGNMENT_STATS,
    add_counts,
    constant_site_counts,
    read_constant_sites,
    strip_invariant_columns,
    write_alignment_stats,
)

_LOG = logging.getLogger("test-alignment")


def _fasta(path: Path, records: dict[str, str]) -> Path:
    # Wrapped lines and lower case, as tools write them.
    text = "".join(f">{name} extra\n{seq[:4]}\n{seq[4:]}\n" for name, seq in records.items())
    path.write_text(text, encoding="utf-8")
    return path


def test_constant_sites_count_one_base_columns_and_ignore_missing(tmp_path: Path) -> None:
    #           columns: 0 A/A/A const A; 1 C/N/C const C; 2 G/T/G variable;
    #           3 N/N/N nothing; 4 t/T/T const T; 5 A/-/C variable; 6 G/G/g const G
    path = _fasta(tmp_path / "full.fasta", {"a": "ACGNtAG", "b": "ANTNT-G", "c": "ACGNTCg"})
    counts, variable = constant_site_counts(path, _LOG)
    assert counts == (1, 1, 1, 1)
    assert variable == 2


def test_unequal_records_are_counted_over_the_shortest(tmp_path: Path, caplog) -> None:
    path = _fasta(tmp_path / "full.fasta", {"a": "AAAAAAAA", "b": "AAAA"})
    with caplog.at_level(logging.WARNING, logger="test-alignment"):
        counts, _ = constant_site_counts(path, _LOG)
    assert counts == (4, 0, 0, 0)
    assert any("differ in length" in r.getMessage() for r in caplog.records)


def test_strip_invariant_columns_keeps_variable_ones_and_returns_the_rest(tmp_path: Path) -> None:
    msa = _fasta(tmp_path / "core.fasta", {"a": "ACGNTAG", "b": "ANTNT-G", "c": "ACGNTCG"})
    out = tmp_path / "variable.fasta"
    removed = strip_invariant_columns(msa, out, _LOG)
    assert removed == (1, 1, 1, 1)  # the all-N column is dropped and not counted
    records = dict(
        zip(
            [line[1:] for line in out.read_text().splitlines() if line.startswith(">")],
            [line for line in out.read_text().splitlines() if not line.startswith(">")],
            strict=True,
        )
    )
    assert records == {"a": "GA", "b": "T-", "c": "GC"}


def test_stats_sidecar_round_trips(tmp_path: Path) -> None:
    core = tmp_path / "core_snp.fasta"
    core.write_text(">a\nA\n")
    path = write_alignment_stats(
        core, (10, 20, 30, 40), 5, source="full_alignment.fasta", masked=True
    )
    assert path.name == ALIGNMENT_STATS
    assert read_constant_sites(core) == (10, 20, 30, 40)
    assert read_constant_sites(tmp_path / "elsewhere" / "x.fasta") is None
    assert add_counts((1, 2, 3, 4), (10, 20, 30, 40)) == (11, 22, 33, 44)


def test_empty_alignment(tmp_path: Path) -> None:
    msa = tmp_path / "empty.fasta"
    msa.write_text("")
    assert constant_site_counts(msa) == ((0, 0, 0, 0), 0)
    assert strip_invariant_columns(msa, tmp_path / "out.fasta") == (0, 0, 0, 0)
