"""Annotator adapters: bakta and prokka argv, outputs and the database check."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from repgenr.annotators import bakta as bk
from repgenr.annotators import prokka as pk
from repgenr.annotators.base import AnnotateParams, registry, require_fasta_section
from repgenr.core.errors import UserInputError, WorkdirError

_LOG = logging.getLogger("test-annot")


def _db(tmp_path: Path, kind: str = "light") -> Path:
    db = tmp_path / "db"
    db.mkdir()
    (db / "version.json").write_text(
        json.dumps({"major": 6, "minor": 0, "type": kind}), encoding="utf-8"
    )
    return db


def _genome(tmp_path: Path) -> Path:
    g = tmp_path / "Fam_Gen_sp_GCA_1.fasta"
    g.write_text(">c1\nACGT\n", encoding="utf-8")
    return g


def _fake(calls: list[list[str]], *, gff_suffix: str, fasta: bool = True):
    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        out = Path(cmd[cmd.index("--output" if cmd[0] == "bakta" else "--outdir") + 1])
        prefix = cmd[cmd.index("--prefix") + 1]
        out.mkdir(parents=True, exist_ok=True)
        body = "##gff-version 3\nc1\tx\tCDS\t1\t3\t.\t+\t0\tID=g1\n"
        if fasta:
            body += "##FASTA\n>c1\nACGT\n"
        (out / f"{prefix}{gff_suffix}").write_text(body, encoding="utf-8")
        (out / f"{prefix}.faa").write_text(">g1\nM\n", encoding="utf-8")
        (out / f"{prefix}.ffn").write_text(">g1\nATG\n", encoding="utf-8")
        return 0

    return run_tool


def test_bakta_argv_outputs_and_extras(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(bk, "run_tool", _fake(calls, gff_suffix=".gff3"))
    db = _db(tmp_path)
    params = AnnotateParams(
        threads=3, db=db, extra={"skip_crispr": "true", "genus": "Francisella", "compliant": 1}
    )
    result = bk.BaktaAnnotator().call(
        _genome(tmp_path), "Fam_Gen_sp_GCA_1", tmp_path / "o", params, _LOG
    )
    (cmd,) = calls
    assert cmd[:3] == ["bakta", "--db", str(db)]
    assert cmd[cmd.index("--prefix") + 1] == "Fam_Gen_sp_GCA_1"
    assert cmd[cmd.index("--threads") + 1] == "3"
    assert "--force" in cmd and "--keep-contig-headers" in cmd
    assert "--skip-plot" in cmd and "--skip-crispr" in cmd and "--compliant" in cmd
    assert cmd[cmd.index("--genus") + 1] == "Francisella"
    assert cmd[-1].endswith("Fam_Gen_sp_GCA_1.fasta")
    assert result.gff.name == "Fam_Gen_sp_GCA_1.gff3" and result.faa and result.ffn
    assert bk.BaktaAnnotator().database_version(db) == "6.0-light"


def test_bakta_skip_plot_can_be_turned_off(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(bk, "run_tool", _fake(calls, gff_suffix=".gff3"))
    bk.BaktaAnnotator().call(
        _genome(tmp_path),
        "Fam_Gen_sp_GCA_1",
        tmp_path / "o",
        AnnotateParams(db=_db(tmp_path), extra={"skip_plot": "false"}),
        _LOG,
    )
    assert "--skip-plot" not in calls[0]


def test_bakta_database_is_required_and_checked(tmp_path: Path, monkeypatch) -> None:
    with pytest.raises(UserInputError, match="bakta_db download"):
        bk.check_bakta_db(None)
    with pytest.raises(UserInputError, match="not a directory"):
        bk.check_bakta_db(tmp_path / "missing")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(UserInputError, match="version.json"):
        bk.check_bakta_db(empty)
    assert bk.check_bakta_db(_db(tmp_path)) == tmp_path / "db"
    monkeypatch.setenv(bk.BAKTA_DB_ENV, str(tmp_path / "db"))
    assert bk.bakta_db_from_env() == str(tmp_path / "db")
    monkeypatch.delenv(bk.BAKTA_DB_ENV)
    assert bk.bakta_db_from_env() is None


def test_bakta_without_fasta_section_is_an_error(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(bk, "run_tool", _fake([], gff_suffix=".gff3", fasta=False))
    with pytest.raises(WorkdirError, match="##FASTA"):
        bk.BaktaAnnotator().call(
            _genome(tmp_path), "x", tmp_path / "o", AnnotateParams(db=_db(tmp_path)), _LOG
        )


def test_prokka_argv_and_gff3_copy(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(pk, "run_tool", _fake(calls, gff_suffix=".gff"))
    result = pk.ProkkaAnnotator().call(
        _genome(tmp_path),
        "Fam_Gen_sp_GCA_1",
        tmp_path / "o",
        AnnotateParams(threads=2, extra={"kingdom": "Bacteria", "gcode": 11, "compliant": "yes"}),
        _LOG,
    )
    (cmd,) = calls
    assert cmd[:1] == ["prokka"] and cmd[cmd.index("--cpus") + 1] == "2"
    assert cmd[cmd.index("--kingdom") + 1] == "Bacteria" and cmd[cmd.index("--gcode") + 1] == "11"
    assert "--compliant" in cmd and "--force" in cmd
    assert result.gff.name == "Fam_Gen_sp_GCA_1.gff3" and result.gff.is_file()
    assert (tmp_path / "o" / "Fam_Gen_sp_GCA_1.gff").is_file()
    assert "##FASTA" in result.gff.read_text(encoding="utf-8")


def test_require_fasta_section(tmp_path: Path) -> None:
    g = tmp_path / "a.gff3"
    g.write_text("##gff-version 3\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="##FASTA"):
        require_fasta_section(g, "x")
    g.write_text("##gff-version 3\n##FASTA\n>a\nA\n", encoding="utf-8")
    require_fasta_section(g, "x")


def test_annotators_are_registered_and_never_auto() -> None:
    assert {"bakta", "prokka"} <= set(registry.names())
    assert registry.get("bakta").requires_database is True
    assert registry.get("prokka").requires_database is False
    assert registry.get("bakta").capabilities.auto_eligible is False
