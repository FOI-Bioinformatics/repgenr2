"""The marker database: packaged lists, download contract, checks, CLI and doctor."""

from __future__ import annotations

import gzip
import json
import logging
import re
import threading
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.doctor import diagnose
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.http import HTTPStatusError
from repgenr.coregenes import markers_db as mdb

_runner = CliRunner()
_LOG = logging.getLogger("test_markers_db")
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def _text(result) -> str:  # noqa: ANN001
    """The output without ANSI codes and with wrapped lines joined: Rich colours and
    wraps the usage box under CI, which can split a flag name."""
    return " ".join(_ANSI.sub("", result.output + str(result.exception or "")).split())


def _hmm_text(acc: str, *, version: int = 7, cutoffs: tuple[str, ...] = ("GA", "TC", "NC")) -> str:
    lines = [
        "HMMER3/f [3.3 | Nov 2019]",
        f"NAME  name_{acc}",
        f"ACC   {acc}.{version}",
        "DESC  test family",
        "LENG  42",
        "ALPH  amino",
    ]
    lines += [f"{c}    20.0 20.0;" for c in cutoffs]
    lines += ["HMM          A        C", "  COMPO   2.1 2.2", "//"]
    return "\n".join(lines) + "\n"


class FakeFetcher:
    """Gzip HMM text per URL; records the accessions asked for; can fail some."""

    def __init__(self, fail: set[str] | None = None, status: int = 503) -> None:
        self.calls: list[str] = []
        self.fail = fail or set()
        self.status = status
        self._lock = threading.Lock()

    def __call__(self, url: str) -> bytes:
        acc = url.rsplit("/", 1)[1].split("?")[0]
        with self._lock:
            self.calls.append(acc)
        if acc in self.fail:
            raise HTTPStatusError(f"HTTP request failed: {url}", self.status)
        return gzip.compress(_hmm_text(acc).encode())


def _download(out: Path, sets, fetch, **kw) -> Path:
    return mdb.download(out, sets, _LOG, fetch=fetch, backoff=0, **kw)


# --- marker lists -----------------------------------------------------------------


def test_marker_lists_have_gtdbtk_counts() -> None:
    bac, ar = mdb.load_marker_list("bac120"), mdb.load_marker_list("ar53")
    assert len(bac) == len(set(bac)) == 120
    assert len(ar) == len(set(ar)) == 53
    assert len(set(bac) & set(ar)) == 5
    assert sum(a.startswith("PF") for a in bac) == 6
    assert sum(a.startswith("TIGR") for a in bac) == 114
    assert sum(a.startswith("PF") for a in ar) == 12
    assert sum(a.startswith("TIGR") for a in ar) == 41
    with pytest.raises(UserInputError):
        mdb.load_marker_list("bac121")


def test_url_source_and_cutoff() -> None:
    assert mdb.marker_url("PF00380").endswith("/pfam/PF00380?annotation=hmm")
    assert mdb.marker_url("TIGR00001").endswith("/ncbifam/TIGR00001?annotation=hmm")
    assert mdb.marker_cutoff(mdb.marker_source("PF00380")) == "ga"
    assert mdb.marker_cutoff(mdb.marker_source("TIGR00001")) == "nc"


# --- parsing ----------------------------------------------------------------------


def test_parse_model_reads_fields_and_refuses_bad_input() -> None:
    m = mdb.parse_model("PF00380", gzip.compress(_hmm_text("PF00380", version=26).encode()))
    assert (m.name, m.version, m.length, m.cutoff) == ("name_PF00380", "PF00380.26", 42, "ga")
    assert m.text.endswith("//\n")
    with pytest.raises(WorkdirError, match="not a HMMER3"):
        mdb.parse_model("PF00380", b"<html>error</html>")
    with pytest.raises(WorkdirError, match="ACC"):
        mdb.parse_model("PF00380", _hmm_text("PF00001"))
    with pytest.raises(WorkdirError, match="no NC line"):
        mdb.parse_model("TIGR00001", _hmm_text("TIGR00001", cutoffs=("GA", "TC")))
    with pytest.raises(WorkdirError, match="one HMMER3 model"):
        mdb.parse_model("PF00380", _hmm_text("PF00380").replace("//\n", ""))
    with pytest.raises(WorkdirError, match="gzip"):
        mdb.parse_model("PF00380", gzip.compress(_hmm_text("PF00380").encode())[:30])


# --- download ---------------------------------------------------------------------


def test_download_writes_index_hmms_and_version(tmp_path: Path) -> None:
    out = tmp_path / "db"
    fetch = FakeFetcher()
    _download(out, ["bac120"], fetch, threads=4)
    assert sorted(fetch.calls) == sorted(mdb.load_marker_list("bac120"))
    rows = mdb.read_index(out)
    assert len(rows) == 120
    assert {r["set"] for r in rows} == {"bac120"}
    first = rows[0]
    assert first == {
        "set": "bac120",
        "accession": "PF00380",
        "source": "pfam",
        "cutoff": "ga",
        "name": "name_PF00380",
        "version": "PF00380.7",
        "length": "42",
    }
    assert {r["cutoff"] for r in rows if r["source"] == "tigrfam"} == {"nc"}
    text = mdb.hmm_path(out, "bac120").read_text()
    assert text.count("HMMER3/") == 120 and text.count("//\n") == 120
    version = json.loads((out / "version.json").read_text())
    assert version["v"] == 1
    assert version["sets"] == ["bac120"]
    assert version["gtdbtk_lists"] == mdb.GTDBTK_LISTS
    assert version["interpro_api"] == mdb.INTERPRO_API
    assert version["n_models"] == {"bac120": 120}
    assert version["fetched"]
    assert mdb.check_markers_db(out) == out
    assert not list(tmp_path.glob(".db.*"))


def test_failed_fetch_leaves_no_directory(tmp_path: Path) -> None:
    out = tmp_path / "db"
    fetch = FakeFetcher(fail={"TIGR00006", "PF00380"})
    with pytest.raises(WorkdirError, match="2 of 120.*PF00380.*TIGR00006"):
        _download(out, ["bac120"], fetch, retries=2)
    assert not out.exists()
    assert list(tmp_path.iterdir()) == []
    # 503 is retried: each failing family was asked for twice.
    assert fetch.calls.count("TIGR00006") == 2


def test_not_found_is_not_retried(tmp_path: Path) -> None:
    fetch = FakeFetcher(fail={"PF00380"}, status=404)
    with pytest.raises(WorkdirError, match="PF00380"):
        _download(tmp_path / "db", ["bac120"], fetch, retries=3)
    assert fetch.calls.count("PF00380") == 1


def test_failed_rerun_keeps_earlier_database(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["bac120"], FakeFetcher())
    before = mdb.markers_db_version(out)
    with pytest.raises(WorkdirError):
        _download(out, ["ar53"], FakeFetcher(fail={"PF04919"}), retries=1)
    assert mdb.markers_db_version(out) == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["db"]


def test_rerun_skips_present_families_and_keeps_sets(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["bac120"], FakeFetcher())
    again = FakeFetcher()
    _download(out, ["bac120"], again)
    assert again.calls == []
    more = FakeFetcher()
    _download(out, ["ar53"], more)
    shared = set(mdb.load_marker_list("bac120")) & set(mdb.load_marker_list("ar53"))
    assert sorted(more.calls) == sorted(set(mdb.load_marker_list("ar53")) - shared)
    version = json.loads((out / "version.json").read_text())
    assert version["sets"] == ["bac120", "ar53"]
    assert version["n_models"] == {"bac120": 120, "ar53": 53}
    assert len(mdb.read_index(out)) == 173
    assert mdb.hmm_path(out, "ar53").read_text().count("HMMER3/") == 53


def test_download_refuses_a_foreign_directory(tmp_path: Path) -> None:
    out = tmp_path / "db"
    out.mkdir()
    (out / "something.txt").write_text("x")
    with pytest.raises(UserInputError, match="not empty"):
        _download(out, ["ar53"], FakeFetcher())


def test_version_changes_with_model_versions(tmp_path: Path) -> None:
    a, b = tmp_path / "a", tmp_path / "b"
    _download(a, ["ar53"], FakeFetcher())
    _download(b, ["ar53"], FakeFetcher())
    assert mdb.markers_db_version(a) == mdb.markers_db_version(b)
    assert mdb.markers_db_version(a).startswith("ar53:")
    tsv = b / "markers.tsv"
    tsv.write_text(tsv.read_text().replace("PF04919.7", "PF04919.8"))
    assert mdb.markers_db_version(a) != mdb.markers_db_version(b)


# --- check_markers_db -------------------------------------------------------------


def test_check_markers_db_refusals_and_env(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv(mdb.MARKERS_DB_ENV, raising=False)
    with pytest.raises(UserInputError, match="repgenr markers-db download"):
        mdb.check_markers_db(None)
    with pytest.raises(UserInputError, match="not a directory"):
        mdb.check_markers_db(tmp_path / "missing")
    with pytest.raises(UserInputError, match="no version.json"):
        mdb.check_markers_db(tmp_path)
    db = tmp_path / "db"
    _download(db, ["ar53"], FakeFetcher())
    monkeypatch.setenv(mdb.MARKERS_DB_ENV, str(db))
    assert mdb.check_markers_db(None) == db
    mdb.hmm_path(db, "ar53").unlink()
    with pytest.raises(UserInputError, match="incomplete"):
        mdb.check_markers_db(db)
    (db / "version.json").write_text("{}")
    with pytest.raises(UserInputError, match="version 1"):
        mdb.check_markers_db(db)


# --- CLI --------------------------------------------------------------------------


def test_cli_download_and_info(tmp_path: Path, monkeypatch) -> None:
    fetch = FakeFetcher()
    real = mdb.download
    monkeypatch.setattr(
        mdb, "download", lambda out, sets, logger, **kw: real(out, sets, logger, fetch=fetch, **kw)
    )
    out = tmp_path / "db"
    result = _runner.invoke(app, ["markers-db", "download", "--output", str(out), "--set", "ar53"])
    assert result.exit_code == 0, result.output
    assert len(fetch.calls) == 53
    result = _runner.invoke(app, ["markers-db", "info", str(out)])
    assert result.exit_code == 0, result.output
    assert "ar53: 53 models (12 Pfam with GA, 41 TIGRFAM with NC)" in result.output
    assert mdb.markers_db_version(out) in result.output


def test_cli_refusals(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv(mdb.MARKERS_DB_ENV, raising=False)
    result = _runner.invoke(app, ["markers-db", "fetch"])
    assert result.exit_code == 2
    result = _runner.invoke(app, ["markers-db", "download"])
    assert result.exit_code == 2 and "--output" in _text(result)
    result = _runner.invoke(app, ["markers-db", "download", "--output", "x", "--set", "bac"])
    assert result.exit_code == 2
    result = _runner.invoke(app, ["markers-db", "info", str(tmp_path)])
    assert result.exit_code == 2


# --- doctor -----------------------------------------------------------------------


def _workdir_with_db_record(tmp_path: Path, db: Path, version: str) -> Path:
    wd = tmp_path / "wd"
    wd.mkdir()
    cfg = Config()
    cfg.record_stage(
        "phylo",
        tool="iqtree",
        params={"markers_db": str(db), "markers_db_version": version},
        completed="2026-01-01T00:00:00",
    )
    cfg.save(wd)
    return wd


def test_doctor_checks_the_recorded_markers_db(tmp_path: Path) -> None:
    db = tmp_path / "db"
    _download(db, ["ar53"], FakeFetcher())
    version = mdb.markers_db_version(db)
    wd = _workdir_with_db_record(tmp_path, db, version)
    findings = [f for f in diagnose(wd) if "marker database" in f.message]
    assert [f.level for f in findings] == ["ok"]

    tsv = db / "markers.tsv"
    tsv.write_text(tsv.read_text().replace("PF04919.7", "PF04919.8"))
    findings = [f for f in diagnose(wd) if "marker database" in f.message]
    assert [f.level for f in findings] == ["warn"]
    assert "not the recorded" in findings[0].message

    (db / "version.json").unlink()
    findings = [f for f in diagnose(wd) if "marker database" in f.message]
    assert [f.level for f in findings] == ["warn"]
    assert "no longer resolves" in findings[0].message


# --- leftovers and per-set version --------------------------------------------------


def test_download_removes_partial_leftovers(tmp_path: Path) -> None:
    out = tmp_path / "db"
    stale = tmp_path / ".db.partial-99999"
    (stale / "hmms").mkdir(parents=True)
    (stale / "hmms" / "ar53.hmm").write_text("half")
    _download(out, ["ar53"], FakeFetcher())
    assert sorted(p.name for p in tmp_path.iterdir()) == ["db"]
    assert mdb.check_markers_db(out) == out


def test_download_restores_a_database_moved_aside(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["ar53"], FakeFetcher())
    before = mdb.markers_db_version(out)
    # A run stopped between its two renames: the earlier database moved aside,
    # the new one never renamed into place.
    out.rename(tmp_path / ".db.old-12345")
    (tmp_path / ".db.partial-12345").mkdir()
    again = FakeFetcher()
    _download(out, ["ar53"], again)
    assert again.calls == []  # the restored database already held every family
    assert mdb.markers_db_version(out) == before
    assert sorted(p.name for p in tmp_path.iterdir()) == ["db"]


def test_old_leftover_is_removed_when_the_database_is_in_place(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["ar53"], FakeFetcher())
    (tmp_path / ".db.old-777").mkdir()
    _download(out, ["ar53"], FakeFetcher())
    assert sorted(p.name for p in tmp_path.iterdir()) == ["db"]


def test_per_set_version_ignores_other_sets(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["bac120"], FakeFetcher())
    bac_only = mdb.markers_db_version(out, ["bac120"])
    whole_before = mdb.markers_db_version(out)
    assert bac_only.startswith("bac120:")
    _download(out, ["ar53"], FakeFetcher())
    assert mdb.markers_db_version(out, ["bac120"]) == bac_only
    assert mdb.markers_db_version(out) != whole_before
    assert mdb.markers_db_version(out).startswith("bac120+ar53:")
    assert mdb.markers_db_version(out, ["ar53"]).startswith("ar53:")
    tsv = out / "markers.tsv"
    tsv.write_text(tsv.read_text().replace("TIGR00006.7", "TIGR00006.8"))
    assert mdb.markers_db_version(out, ["bac120"]) != bac_only


def test_per_set_version_refuses_a_missing_set(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["ar53"], FakeFetcher())
    with pytest.raises(UserInputError, match="--set bac120"):
        mdb.markers_db_version(out, ["bac120"])


def test_cli_info_refuses_download_only_flags(tmp_path: Path) -> None:
    out = tmp_path / "db"
    _download(out, ["ar53"], FakeFetcher())
    for extra in (["--set", "ar53"], ["-t", "2"], ["--threads", "2"]):
        result = _runner.invoke(app, ["markers-db", "info", str(out), *extra])
        assert result.exit_code == 2, (extra, result.output)
        assert "download only" in _text(result)
    assert _runner.invoke(app, ["markers-db", "info", str(out)]).exit_code == 0
