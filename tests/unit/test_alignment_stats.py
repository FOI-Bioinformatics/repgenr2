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


def test_sidecar_is_bound_to_its_alignment_by_name_and_content(tmp_path: Path) -> None:
    from repgenr.tree.alignment import copy_alignment_stats

    core = tmp_path / "core_snp.fasta"
    core.write_text(">a\nA\n>b\nC\n")
    full = tmp_path / "full_alignment.fasta"
    full.write_text(">a\nAAAA\n>b\nAAAC\n")
    write_alignment_stats(core, (3, 0, 0, 0), 1, source=full.name, masked=False)
    assert read_constant_sites(core) == (3, 0, 0, 0)
    assert read_constant_sites(full) is None  # the whole-genome alignment holds its own constants
    core.write_text(">a\nA\n>b\nG\n")  # rewritten core: the old counts no longer apply
    assert read_constant_sites(core) is None
    core.write_text(">a\nA\n>b\nC\n")
    copied = tmp_path / "elsewhere" / "msa.fasta"
    copied.parent.mkdir()
    copied.write_text(core.read_text())
    assert copy_alignment_stats(core, copied) is not None
    assert read_constant_sites(copied) == (3, 0, 0, 0)
    assert copy_alignment_stats(full, tmp_path / "elsewhere" / "other.fasta") is None


def test_records_longer_than_the_first_are_reported(tmp_path: Path, caplog) -> None:
    path = _fasta(tmp_path / "full.fasta", {"a": "AAAA", "b": "AAAAAAAA"})
    with caplog.at_level(logging.WARNING, logger="test-alignment"):
        counts, _ = constant_site_counts(path, _LOG)
    assert counts == (4, 0, 0, 0)
    assert any("4 to 8" in r.getMessage() for r in caplog.records)


def test_a_record_without_a_name_is_refused(tmp_path: Path) -> None:
    import pytest

    from repgenr.core.errors import WorkdirError

    path = tmp_path / "x.fasta"
    path.write_text(">\nACGT\n>b\nACGT\n")
    with pytest.raises(WorkdirError, match="no name"):
        constant_site_counts(path)
    with pytest.raises(WorkdirError, match="no name"):
        strip_invariant_columns(path, tmp_path / "out.fasta")


def test_rna_u_counts_as_t(tmp_path: Path) -> None:
    path = _fasta(tmp_path / "rna.fasta", {"a": "UUUU", "b": "TTTU"})
    assert constant_site_counts(path) == ((0, 0, 0, 4), 0)


def test_completeness_filter_keeps_columns_called_in_enough_genomes(tmp_path: Path) -> None:
    from repgenr.tree.alignment import columns_at_completeness, filter_columns_by_completeness

    # column 0 called in 3/3, column 1 in 2/3, column 2 in 1/3, column 3 in 3/3
    msa = _fasta(tmp_path / "core.fasta", {"a": "ACGT", "b": "ANNT", "c": "ACNT"})
    assert columns_at_completeness(msa, 0.95) == (2, 4)  # columns 0 and 3
    assert columns_at_completeness(msa, 0.6) == (3, 4)  # 2 of 3 is 0.67
    out = tmp_path / "soft.fasta"
    assert filter_columns_by_completeness(msa, out, 0.6, _LOG) == (3, 1)
    seqs = [ln for ln in out.read_text().splitlines() if not ln.startswith(">")]
    assert seqs == ["ACT", "ANT", "ACT"]


def test_constant_sites_follow_the_same_completeness_rule(tmp_path: Path) -> None:
    # Constant A in all three; constant C called in 2 of 3; variable G/T in 2 of 3.
    path = _fasta(tmp_path / "full.fasta", {"a": "ACGT", "b": "ANTT", "c": "ACNT"})
    assert constant_site_counts(path)[0] == (1, 1, 0, 1)
    counts, variable = constant_site_counts(path, min_fraction=0.95)
    assert counts == (1, 0, 0, 1) and variable == 0
