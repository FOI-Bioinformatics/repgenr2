"""Every input genome must leave dereplication with a status and a home."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)


def _ok() -> DerepResult:
    return DerepResult(
        representatives=[Path("a.fasta")],
        clusters={"a.fasta": ["b.fasta"]},
        genome_status={"a.fasta": STATUS_REPRESENTATIVE, "b.fasta": STATUS_CONTAINED},
    )


def test_complete_result_passes() -> None:
    check_result_complete(_ok(), ["a.fasta", "b.fasta"])


def test_missing_status_raises() -> None:
    with pytest.raises(WorkdirError, match="c.fasta"):
        check_result_complete(_ok(), ["a.fasta", "b.fasta", "c.fasta"])


def test_contained_without_cluster_raises() -> None:
    r = _ok()
    r.clusters = {"a.fasta": []}
    with pytest.raises(WorkdirError, match="no representative"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_representative_without_status_raises() -> None:
    r = _ok()
    r.genome_status = {"a.fasta": STATUS_CONTAINED, "b.fasta": STATUS_CONTAINED}
    with pytest.raises(WorkdirError, match="not marked as such"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_unknown_status_value_raises() -> None:
    r = _ok()
    r.genome_status["b.fasta"] = "weird"
    with pytest.raises(WorkdirError, match="weird"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_doubled_membership_raises() -> None:
    r = DerepResult(
        representatives=[Path("a.fasta"), Path("c.fasta")],
        clusters={"a.fasta": ["b.fasta"], "c.fasta": ["b.fasta"]},
        genome_status={
            "a.fasta": STATUS_REPRESENTATIVE,
            "c.fasta": STATUS_REPRESENTATIVE,
            "b.fasta": STATUS_CONTAINED,
        },
    )
    with pytest.raises(WorkdirError, match="more than one cluster"):
        check_result_complete(r, ["a.fasta", "b.fasta", "c.fasta"])


def test_fail_qc_needs_no_cluster() -> None:
    r = DerepResult(
        representatives=[Path("a.fasta")],
        clusters={"a.fasta": []},
        genome_status={"a.fasta": STATUS_REPRESENTATIVE, "b.fasta": "fail_qc"},
    )
    check_result_complete(r, ["a.fasta", "b.fasta"])


# --- the settings of several adapter calls ---------------------------------------


def _with(settings: dict) -> DerepResult:
    result = _ok()
    result.settings = settings
    return result


def test_settings_default_to_empty() -> None:
    assert _ok().settings == {}


def test_combine_settings_keeps_shared_values_and_drops_differing_counts() -> None:
    from repgenr.dereplicators.base import combine_settings

    base = {"backend": "dense", "ksize": 31, "scaled": 1000}
    combined = combine_settings(
        [_with({**base, "unknown_pairs": 2}), _ok(), _with({**base, "unknown_pairs": 5})]
    )
    assert combined == base  # the result without settings takes no part
    same = combine_settings([_with({**base, "unknown_pairs": 2})] * 2)
    assert same == {**base, "unknown_pairs": 2}


def test_combine_settings_joins_differing_back_ends() -> None:
    from repgenr.dereplicators.base import combine_settings

    combined = combine_settings(
        [
            _with({"backend": "sparse", "ksize": 31}),
            _with({"backend": "dense", "ksize": 31, "unknown_pairs": 0}),
        ]
    )
    assert combined == {"backend": "dense,sparse", "ksize": 31}
    assert combine_settings([_ok(), _ok()]) == {}
