"""CLI input validation and typed replacements for production asserts (WP1-5)."""

from __future__ import annotations

import logging
import re
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.errors import UserInputError
from repgenr.core.plugins import parse_extra_int

_runner = CliRunner()
_LOG = logging.getLogger("test")

_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def _plain(output: str) -> str:
    """Strip ANSI color codes (CI terminals get Rich-styled usage errors)."""
    return _ANSI.sub("", output)


# --- Typer range constraints (usage errors exit 2) ----------------------------


def test_threads_zero_rejected(tmp_path) -> None:
    result = _runner.invoke(app, ["dereplicate", "-wd", str(tmp_path), "--threads", "0"])
    assert result.exit_code == 2
    assert "--threads" in _plain(result.output)


def test_bootstrap_negative_rejected(tmp_path) -> None:
    result = _runner.invoke(app, ["phylo", "-wd", str(tmp_path), "--bootstrap", "-1"])
    assert result.exit_code == 2
    assert "--bootstrap" in _plain(result.output)


def test_limit_zero_rejected(tmp_path) -> None:
    result = _runner.invoke(
        app,
        ["metadata", "-wd", str(tmp_path), "-tg", "Francisella", "--limit", "0"],
    )
    assert result.exit_code == 2
    assert "--limit" in _plain(result.output)


def test_length_deviation_negative_rejected(tmp_path) -> None:
    result = _runner.invoke(app, ["vgenome", "-wd", str(tmp_path), "--length-deviation", "-1"])
    assert result.exit_code == 2
    assert "--length-deviation" in _plain(result.output)


# --- --released-after parsed at the CLI boundary ------------------------------


@pytest.mark.parametrize("raw", ["2020-01-01", "31/12/2020", "not-a-date", "13/01/2020"])
def test_released_after_malformed_rejected(tmp_path, raw) -> None:
    result = _runner.invoke(
        app,
        ["vmetadata", "-wd", str(tmp_path), "--target", "coronavirus", "--released-after", raw],
    )
    assert result.exit_code == 2
    assert "MM/DD/YYYY" in _plain(result.output)


# --- parse_extra_int helper ---------------------------------------------------


def test_parse_extra_int_valid() -> None:
    assert parse_extra_int({"ksize": "31"}, "ksize", 21) == 31
    assert parse_extra_int({"ksize": 31}, "ksize", 21) == 31


def test_parse_extra_int_default() -> None:
    assert parse_extra_int({}, "ksize", 21) == 21


def test_parse_extra_int_malformed_raises_user_input_error() -> None:
    with pytest.raises(UserInputError, match="ksize"):
        parse_extra_int({"ksize": "abc"}, "ksize", 21)


# --- typed raises replacing production asserts --------------------------------


def test_metadata_target_taxon_missing_name_raises() -> None:
    from repgenr.sources.gtdb_api import target_taxon
    from repgenr.stages.metadata import MetadataParams

    params = MetadataParams(dataset="all", level="family", target_family=None)
    with pytest.raises(UserInputError, match="family"):
        target_taxon(params)


def test_metadata_obtain_requires_release(tmp_path) -> None:
    from repgenr.sources.gtdb_table import _obtain_metadata
    from repgenr.stages.metadata import MetadataParams

    ctx = SimpleNamespace(workdir=tmp_path)
    params = MetadataParams(dataset="all", level="family", release=None, metadata_path=None)
    with pytest.raises(UserInputError, match="release"):
        _obtain_metadata(ctx, params, _LOG)


def test_vmetadata_ncbi_virus_requires_target(tmp_path) -> None:
    from repgenr.stages.vmetadata import _run_ncbi_virus

    ctx = SimpleNamespace(workdir=tmp_path)
    params = SimpleNamespace(target=None)
    with pytest.raises(UserInputError, match="target"):
        _run_ncbi_virus(ctx, params, tmp_path, _LOG)
