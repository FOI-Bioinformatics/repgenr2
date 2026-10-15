"""`repgenr doctor`: read-only workdir health checks (outputs vs records)."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.contracts import SelectionRow, write_clusters, write_selection
from repgenr.core.doctor import diagnose
from repgenr.core.errors import WorkdirError
from repgenr.core.inputs import dir_stat_digest, manifest_digest_for_stage
from repgenr.core.manifest import GenomeRecord, Manifest

_runner = CliRunner()


def _levels(findings, area: str | None = None) -> set[str]:
    return {f.level for f in findings if area is None or f.area == area}


def _messages(findings, level: str | None = None) -> str:
    return "\n".join(f.message for f in findings if level is None or f.level == level)


def _base_workdir(tmp_path: Path) -> Path:
    """A healthy two-genome workdir: selection + manifest + genomes + outgroup."""
    wd = tmp_path / "wd"
    wd.mkdir()
    rows = [
        SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", False, "Fam_Gen_sp2_GCF_2.1.fasta"),
        SelectionRow("GCF_9.1", "Fam", "Out", "grp", True, "Fam_Out_grp_GCF_9.1.fasta"),
    ]
    write_selection(wd / "selection.tsv", rows)
    (wd / "genomes").mkdir()
    for r in rows[:2]:
        (wd / "genomes" / r.filename).write_text(">x\nACGT\n", encoding="utf-8")
    (wd / "outgroup").mkdir()
    (wd / "outgroup" / rows[2].filename).write_text(">og\nACGT\n", encoding="utf-8")
    (wd / "outgroup_accession.txt").write_text("GCF_9.1\n", encoding="utf-8")
    manifest = Manifest(wd / "manifest.sqlite")
    manifest.replace_genomes(
        [
            GenomeRecord(accession="GCF_1.1", filename=rows[0].filename),
            GenomeRecord(accession="GCF_2.1", filename=rows[1].filename),
            GenomeRecord(accession="GCF_9.1", filename=rows[2].filename, is_outgroup=True),
        ]
    )
    manifest.close()
    cfg = Config()
    cfg.record_stage("metadata", completed="2026-01-01T00:00:00")
    cfg.record_stage("genome", completed="2026-01-01T00:01:00")
    cfg.save(wd)
    return wd


def test_healthy_workdir_has_no_failures(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    findings = diagnose(wd)
    assert "fail" not in _levels(findings), _messages(findings, "fail")


def test_no_workdir_reports_cleanly(tmp_path: Path) -> None:
    findings = diagnose(tmp_path / "nope")
    assert findings and findings[0].level == "warn"


def test_interrupted_stage_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    cfg = Config.load(wd)
    cfg.record_stage("dereplicate", tool="skder", params={"x": 1})  # no completed
    cfg.save(wd)
    findings = diagnose(wd)
    assert any(f.level == "fail" and "dereplicate" in f.area for f in findings)
    assert "re-run" in _messages(findings, "fail")


def test_missing_genome_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    findings = diagnose(wd)
    assert any(f.level == "fail" and "Fam_Gen_sp2_GCF_2.1.fasta" in f.message for f in findings)


def test_excused_missing_accession_is_not_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    (wd / "missing_accessions.txt").write_text("GCF_2.1\n", encoding="utf-8")
    findings = diagnose(wd)
    assert "fail" not in _levels(findings), _messages(findings, "fail")


def test_non_fasta_genome_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp1_GCF_1.1.fasta").write_text(
        "<html>error</html>", encoding="utf-8"
    )
    findings = diagnose(wd)
    assert any(f.level == "fail" and "Fam_Gen_sp1_GCF_1.1.fasta" in f.message for f in findings)


def test_manifest_selection_drift_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    manifest = Manifest(wd / "manifest.sqlite")
    manifest.upsert_many([GenomeRecord(accession="GCF_STALE.1")])
    manifest.close()
    findings = diagnose(wd)
    assert any(f.level == "fail" and "GCF_STALE.1" in f.message for f in findings)


def test_representatives_clusters_mismatch_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    reps = wd / "derep" / "representatives"
    reps.mkdir(parents=True)
    (reps / "Fam_Gen_sp1_GCF_1.1.fasta").write_text(">x\nACGT\n", encoding="utf-8")
    write_clusters(
        wd / "derep" / "clusters.tsv",
        {
            "Fam_Gen_sp1_GCF_1.1.fasta": [],
            "Fam_Gen_sp2_GCF_2.1.fasta": [],  # listed but absent on disk
        },
    )
    findings = diagnose(wd)
    assert any(f.level == "fail" and "Fam_Gen_sp2_GCF_2.1.fasta" in f.message for f in findings)


def test_truncated_tree_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree").mkdir()
    (wd / "tree" / "tree.nwk").write_text("(a,b", encoding="utf-8")  # no ';'
    findings = diagnose(wd)
    assert any(f.level == "fail" and "tree.nwk" in f.message for f in findings)


def test_mismatched_tree2tax_pair_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text("child\tparent\n", encoding="utf-8")
    # genomes_map.tsv missing
    findings = diagnose(wd)
    assert any(f.level == "fail" and "genomes_map.tsv" in f.message for f in findings)


def test_unresolvable_outgroup_is_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "outgroup" / "Fam_Out_grp_GCF_9.1.fasta").unlink()
    findings = diagnose(wd)
    assert any(f.level == "warn" and "GCF_9.1" in f.message for f in findings)


def test_outgroup_resolved_only_by_a_temporary_file_is_a_warning(tmp_path: Path) -> None:
    # doctor applies phylo's rule: a partial GCF_9.1 download is no outgroup.
    wd = _base_workdir(tmp_path)
    og = wd / "outgroup" / "Fam_Out_grp_GCF_9.1.fasta"
    og.rename(og.with_name(og.name + ".tmp"))
    findings = diagnose(wd)
    assert any(f.level == "warn" and f.area == "outgroup" for f in findings)
    assert not any(f.level == "ok" and f.area == "outgroup" for f in findings)


def test_leftover_temp_files_are_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree").mkdir()
    (wd / "tree" / "tree.nwk.part").write_text("(a", encoding="utf-8")
    findings = diagnose(wd)
    assert any(f.level == "warn" and "tree.nwk.part" in f.message for f in findings)


def _legacy_snp_stamp(wd: Path, table: str, stamped: str) -> None:
    """snp/core_snp.fasta and a stamp from the earlier phylo layout naming the
    digest of ``stamped``."""
    import hashlib
    import json

    (wd / "snp").mkdir()
    (wd / "snp" / "core_snp.fasta").write_text(table, encoding="utf-8")
    digest = hashlib.sha256(stamped.encode("utf-8")).hexdigest()
    (wd / "snp" / "msa_source.json").write_text(
        json.dumps({"artifact_digest": digest}), encoding="utf-8"
    )


def test_phylo_stamp_left_in_snp_is_a_warning(tmp_path: Path) -> None:
    """Before tree/msa/, phylo's typing pass wrote snp/ and its stamp there;
    while the stamp still describes snp/core_snp.fasta, the tables are phylo's,
    not the snptype stage's."""
    wd = _base_workdir(tmp_path)
    _legacy_snp_stamp(wd, ">og\nACGT\n", ">og\nACGT\n")
    warned = [f for f in diagnose(wd) if f.level == "warn" and f.area == "snptype"]
    assert len(warned) == 1 and "snp/msa_source.json" in warned[0].message


def test_phylo_stamp_in_snp_replaced_by_snptype_is_not_a_warning(tmp_path: Path) -> None:
    """A snptype run after the old phylo pass rewrote snp/core_snp.fasta; the
    stamp no longer matches it and the tables are the snptype stage's (#223
    ran snptype after phylo, so such workdirs are common)."""
    wd = _base_workdir(tmp_path)
    _legacy_snp_stamp(wd, ">from_snptype\nTTTT\n", ">og\nACGT\n")
    assert not [f for f in diagnose(wd) if f.area == "snptype"]
    (wd / "snp" / "msa_source.json").write_text("not json", encoding="utf-8")
    assert not [f for f in diagnose(wd) if f.area == "snptype"]


def test_phylo_stamp_under_tree_msa_is_not_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree" / "msa").mkdir(parents=True)
    (wd / "tree" / "msa" / "msa_source.json").write_text("{}", encoding="utf-8")
    assert not [f for f in diagnose(wd) if f.area == "snptype"]


def test_changed_inputs_are_a_warning(tmp_path: Path) -> None:
    """A completed stage whose recorded input digests no longer match reality
    is stale (it will re-run) -- doctor should say so."""
    wd = _base_workdir(tmp_path)
    cfg = Config.load(wd)
    cfg.record_stage(
        "dereplicate",
        tool="skder",
        params={},
        completed="2026-01-01T00:02:00",
        fingerprint="f",
        inputs={"genomes": "old-digest"},
    )
    cfg.save(wd)
    findings = diagnose(wd)
    assert any(
        f.level == "warn" and f.area == "dereplicate" and "changed" in f.message for f in findings
    )


def test_missing_deliverable_is_a_warning_named_like_the_resume_message(
    tmp_path: Path,
) -> None:
    """A completed stage whose deliverable was deleted will re-run; doctor names
    the same workdir-relative path that the resume log line does."""
    wd = _base_workdir(tmp_path)
    cfg = Config.load(wd)
    cfg.record_stage("phylo", tool="fasttree", params={}, completed="2026-01-01T00:03:00")
    cfg.save(wd)
    findings = diagnose(wd)
    hits = [f for f in findings if f.area == "phylo" and "deliverable" in f.message]
    assert len(hits) == 1
    assert hits[0].level == "warn"
    assert "deliverable tree/tree.nwk missing" in hits[0].message
    # The healthy base stages report no missing deliverable.
    assert not any("deliverable" in f.message for f in findings if f.area != "phylo")


def test_dereplicate_completion_with_derep_status_is_not_stale(tmp_path: Path) -> None:
    """doctor's stale-input check must use the same per-stage manifest digest
    dereplicate's own resume fingerprint stamps (include_derep=False, since
    dereplicate WRITES derep_status/representative itself) -- otherwise every
    dereplicated workdir falsely reports "input(s) changed" for dereplicate on
    every run, forever, since the two digests can never agree."""
    wd = _base_workdir(tmp_path)
    manifest = Manifest(wd / "manifest.sqlite")
    manifest.set_derep_status_many(
        [
            ("GCF_1.1", "representative", None),
            ("GCF_2.1", "contained", "Fam_Gen_sp1_GCF_1.1.fasta"),
        ]
    )
    manifest.close()

    genomes_digest = dir_stat_digest(wd / "genomes")
    ro = Manifest.open_readonly(wd / "manifest.sqlite")
    try:
        manifest_digest_value = manifest_digest_for_stage("dereplicate", ro)
    finally:
        ro.close()

    cfg = Config.load(wd)
    cfg.record_stage(
        "dereplicate",
        tool="skder",
        params={},
        completed="2026-01-01T00:02:00",
        fingerprint="f",
        inputs={"genomes": genomes_digest, "manifest": manifest_digest_value},
    )
    cfg.save(wd)

    # The record has no derep/ outputs on disk, so only the stale-input
    # findings (not the missing-deliverable ones) are under test here.
    def _stale(findings) -> list:
        return [f for f in findings if f.area == "dereplicate" and "changed" in f.message]

    findings = diagnose(wd)
    assert not _stale(findings)

    # A real quality-only manifest edit must still be caught as stale.
    manifest2 = Manifest(wd / "manifest.sqlite")
    manifest2.upsert_many(
        [
            GenomeRecord(
                accession="GCF_1.1",
                filename="Fam_Gen_sp1_GCF_1.1.fasta",
                completeness=98.0,
                contamination=0.5,
            )
        ]
    )
    manifest2.close()

    findings2 = diagnose(wd)
    dereplicate_warnings = _stale(findings2)
    assert len(dereplicate_warnings) == 1
    assert dereplicate_warnings[0].level == "warn"
    assert "manifest" in dereplicate_warnings[0].message


def test_diagnose_is_read_only(tmp_path: Path) -> None:
    """Doctor must not create files in a workdir (notably no manifest.sqlite)."""
    wd = tmp_path / "wd"
    wd.mkdir()
    Config().save(wd)
    before = sorted(p.name for p in wd.iterdir())
    diagnose(wd)
    assert sorted(p.name for p in wd.iterdir()) == before


# --- CLI ----------------------------------------------------------------------


def _non_wal_listing(workdir: Path) -> list[str]:
    # A WAL-mode database opened read-only may cause sqlite to create a
    # "-shm" (and sometimes "-wal") sidecar purely to read the wal-index;
    # a read-only connection cannot check those back in on close. They are
    # not a write to the manifest itself, so they are excluded here -- the
    # mtime assertion below is what actually proves the manifest untouched.
    return sorted(
        p.name
        for p in workdir.iterdir()
        if not (p.name.endswith("-shm") or p.name.endswith("-wal"))
    )


def test_doctor_leaves_manifest_untouched(workdir: Path) -> None:
    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1", filename="x.fasta"))
    (workdir / "selection.tsv").write_text("accession\tfilename\nGCA_1\tx.fasta\n")
    Config().save(workdir)  # so diagnose() reaches the manifest-drift check
    manifest = workdir / "manifest.sqlite"
    before = (manifest.stat().st_mtime_ns, _non_wal_listing(workdir))

    diagnose(workdir)

    after = (manifest.stat().st_mtime_ns, _non_wal_listing(workdir))
    assert before == after


def test_open_readonly_refuses_writes(workdir: Path) -> None:
    import sqlite3

    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    ro = Manifest.open_readonly(workdir / "manifest.sqlite")
    try:
        assert [g.accession for g in ro.all_genomes()] == ["GCA_1"]
        with pytest.raises(sqlite3.OperationalError):
            ro.upsert(GenomeRecord(accession="GCA_2"))
    finally:
        ro.close()


def test_open_readonly_missing_file_raises_workdir_error(tmp_path: Path) -> None:
    with pytest.raises(WorkdirError):
        Manifest.open_readonly(tmp_path / "absent.sqlite")


def test_open_readonly_unwritable_dir_raises_workdir_error(workdir: Path) -> None:
    """A WAL manifest needs a writable directory for its -shm file even to be
    opened read-only; that failure must surface as WorkdirError, not a raw
    sqlite3.OperationalError, so doctor's callers get the documented contract."""
    if os.geteuid() == 0:
        pytest.skip("root bypasses directory permissions")
    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    manifest_path = workdir / "manifest.sqlite"
    original_mode = stat.S_IMODE(workdir.stat().st_mode)
    workdir.chmod(0o555)
    try:
        with pytest.raises(WorkdirError):
            Manifest.open_readonly(manifest_path)
    finally:
        workdir.chmod(original_mode)


def test_open_readonly_unreadable_file_raises_workdir_error(workdir: Path) -> None:
    """Even connect() itself can fail (manifest file unreadable, e.g. chmod
    0o000); that must also surface as WorkdirError, not a raw sqlite error."""
    if os.geteuid() == 0:
        pytest.skip("root bypasses file permissions")
    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    manifest_path = workdir / "manifest.sqlite"
    original_mode = stat.S_IMODE(manifest_path.stat().st_mode)
    manifest_path.chmod(0o000)
    try:
        with pytest.raises(WorkdirError):
            Manifest.open_readonly(manifest_path)
    finally:
        manifest_path.chmod(original_mode)


def test_cli_doctor_exit_codes(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "OK" in result.stdout or "ok" in result.stdout

    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 7  # failures found; 1 is kept for a crash
    assert "FAIL" in result.stdout


@pytest.mark.parametrize(
    "text",
    [
        "stages: [\n  bad",
        "- a\n- b\n",
        "stages:\n  phylo:\n    params: oops\n",
        "stages:\n  dereplicate: [1, 2]\n",
        "stages:\n  dereplicate:\n    tool_versions: notamap\n",
    ],
)
def test_malformed_record_is_a_config_failure_not_a_traceback(tmp_path: Path, text: str) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "repgenr.yaml").write_text(text, encoding="utf-8")
    findings = diagnose(wd)
    assert [(f.level, f.area) for f in findings] == [("fail", "config")]
    assert "not a readable RepGenR2 record" in findings[0].message
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 7
    assert "[FAIL] config" in result.stdout
    for command in ("status", "versions"):
        result = _runner.invoke(app, [command, "-wd", str(wd)])
        assert result.exit_code == 3, result.output
        assert "not a readable RepGenR2 record" in result.output


def test_interrupted_stage_without_params_is_a_failure(tmp_path: Path) -> None:
    # cluster_summary has no parameters: a run killed on its first attempt
    # leaves a record with neither params nor tool, which status already
    # listed as interrupted while doctor reported nothing.
    wd = _base_workdir(tmp_path)
    cfg = Config.load(wd)
    cfg.record_stage("cluster_summary")
    cfg.record_stage("dereplicate")
    cfg.save(wd)
    findings = diagnose(wd)
    failed = {f.area for f in findings if f.level == "fail"}
    assert {"cluster_summary", "dereplicate"} <= failed
    status = _runner.invoke(app, ["status", "-wd", str(wd)])
    assert "[interrupted] dereplicate" in status.stdout
    assert "cluster_summary  [interrupted]" in status.stdout


def _tree2tax_pair(wd: Path, edges: str, mapping: str) -> None:
    (wd / "tree2tax.tsv").write_text("child\tparent\n" + edges, encoding="utf-8")
    (wd / "genomes_map.tsv").write_text(mapping, encoding="utf-8")


_EDGES = "a\tn1\nb\tn1\nn1\troot\nc\troot\n"
_MAP = "A1\ta\nB1\tb\nB2\tb\nC1\tc\n"


def test_consistent_tree2tax_tables_pass(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    _tree2tax_pair(wd, _EDGES, _MAP)
    assert "fail" not in _levels(diagnose(wd), "tree2tax")


@pytest.mark.parametrize(
    ("edges", "mapping", "expected"),
    [
        ("", _MAP, "holds no edges"),  # tree2tax.tsv emptied to its header
        (_EDGES, "", "is empty"),  # genomes_map.tsv emptied
        (_EDGES, "A1\ta\nB1\tb\n", "name different leaves"),  # map truncated
        ("a\tn1\nb\tn1\n", _MAP, "name different leaves"),  # tree2tax truncated
    ],
)
def test_truncated_tree2tax_tables_are_a_failure(
    tmp_path: Path, edges: str, mapping: str, expected: str
) -> None:
    wd = _base_workdir(tmp_path)
    _tree2tax_pair(wd, edges, mapping)
    failures = _messages([f for f in diagnose(wd) if f.area == "tree2tax"], "fail")
    assert expected in failures


def test_headerless_tree2tax_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text("", encoding="utf-8")
    (wd / "genomes_map.tsv").write_text(_MAP, encoding="utf-8")
    assert "lacks its child/parent header" in _messages(diagnose(wd), "fail")


def test_untracked_genome_is_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp3_GCF_3.1.fasta").write_text(">x\nACGT\n", encoding="utf-8")
    findings = diagnose(wd)
    warned = _messages([f for f in findings if f.area == "genomes"], "warn")
    assert "1 file(s)" in warned and "Fam_Gen_sp3_GCF_3.1.fasta" in warned
    assert "re-run metadata" in warned


def test_advice_names_the_stage_that_wrote_the_genome_set(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    cfg = Config()
    cfg.record_stage("ingest", completed="2026-01-01T00:00:00")
    cfg.save(wd)
    (wd / "genomes" / "Fam_Gen_sp1_GCF_1.1.fasta").unlink()
    failures = _messages(diagnose(wd), "fail")
    assert "re-run ingest" in failures
    assert "genome stage" not in failures


def test_doctor_prints_findings_only(tmp_path: Path) -> None:
    # The completeness guard's refusal text used to reach the console as an
    # unformatted warning above the [FAIL] line that reports the same thing.
    # A subprocess: pytest attaches its own handlers to every logger.
    import subprocess
    import sys

    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp1_GCF_1.1.fasta").unlink()
    result = subprocess.run(
        [sys.executable, "-c", "from repgenr.cli.main import app; app()", "doctor", "-wd", str(wd)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 7, result.stderr
    assert result.stderr == ""
    assert "--allow-incomplete" not in result.stdout


def test_appledouble_companions_are_not_leftovers(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree").mkdir()
    (wd / "tree" / "tree.nwk.part").write_text("(a", encoding="utf-8")
    (wd / "tree" / "._tree.nwk.part").write_bytes(b"\0" * 16)
    warned = _messages([f for f in diagnose(wd) if f.area == "leftovers"], "warn")
    assert "1 temp file(s)" in warned and "._tree" not in warned


def test_emptied_record_beside_outputs_is_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "repgenr.yaml").write_text("", encoding="utf-8")
    warned = _messages([f for f in diagnose(wd) if f.area == "config"], "warn")
    assert "records no stage" in warned


def _doctor_json(wd: Path, *extra: str) -> tuple[int, dict]:
    import json

    result = _runner.invoke(app, ["doctor", "-wd", str(wd), "--json", *extra])
    return result.exit_code, json.loads(result.stdout)


def test_doctor_json_on_a_healthy_workdir(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    code, payload = _doctor_json(wd)
    assert code == 0
    assert payload["schema"] == "repgenr.doctor/1"
    assert payload["workdir"] == str(wd)
    assert payload["quick"] is False
    assert payload["exit_code"] == 0
    assert payload["counts"]["fail"] == 0
    assert set(payload["findings"][0]) == {"level", "area", "message"}
    for level in ("fail", "warn", "ok"):
        listed = sum(1 for f in payload["findings"] if f["level"] == level)
        assert payload["counts"][level] == listed


def test_doctor_json_counts_agree_with_the_exit_code(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").unlink()
    code, payload = _doctor_json(wd)
    assert code == 7 == payload["exit_code"]
    assert payload["counts"]["fail"] >= 1
    assert payload["findings"][0]["level"] == "fail"  # failures first, as in the text
    assert any(f["area"] == "genomes" and f["level"] == "fail" for f in payload["findings"])


def test_doctor_json_on_a_missing_workdir_leaves_stdout_empty(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["doctor", "-wd", str(tmp_path / "missing"), "--json"])
    assert result.exit_code == 3
    assert result.stdout == ""


def test_a_check_that_raises_exits_7_not_1(tmp_path: Path, monkeypatch) -> None:
    import repgenr.core.doctor as doctor_mod

    def _check_tree(workdir: Path, config: Config) -> list:
        raise RuntimeError("unreadable table")

    monkeypatch.setattr(doctor_mod, "_check_tree", _check_tree)
    wd = _base_workdir(tmp_path)
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 7, result.output
    assert "[FAIL] tree: Check could not complete: unreadable table" in result.stdout


def test_quick_check_failure_keeps_the_genomes_area(tmp_path: Path, monkeypatch) -> None:
    # The genome check is wrapped in functools.partial under quick; a failure
    # inside it still names its area.
    import repgenr.core.doctor as doctor_mod

    def boom(*args, **kwargs):
        raise RuntimeError("stat failed")

    monkeypatch.setattr(doctor_mod, "looks_like_fasta", boom)
    findings = diagnose(_base_workdir(tmp_path))
    assert ("fail", "genomes") in {(f.level, f.area) for f in findings}


def test_quick_skips_the_fasta_content_check(tmp_path: Path, monkeypatch) -> None:
    import repgenr.core.doctor as doctor_mod

    wd = _base_workdir(tmp_path)
    (wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").write_text("<html>\n", encoding="utf-8")
    full = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert full.exit_code == 7
    assert "not FASTA" in full.stdout

    def never(path: Path) -> bool:
        raise AssertionError(f"looks_like_fasta called under --quick: {path}")

    monkeypatch.setattr(doctor_mod, "looks_like_fasta", never)
    quick = _runner.invoke(app, ["doctor", "-wd", str(wd), "--quick"])
    assert quick.exit_code == 0, quick.output
    assert "not FASTA" not in quick.stdout
    assert "2 genome file(s) present (content not read: --quick)" in quick.stdout
    code, payload = _doctor_json(wd, "--quick")
    assert code == 0
    assert payload["quick"] is True


def test_quick_still_reports_missing_and_dangling_genomes(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    victim = wd / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta"
    victim.unlink()
    victim.symlink_to(tmp_path / "gone.fasta")
    findings = diagnose(wd, quick=True)
    text = _messages(findings, "fail")
    assert "point at files that no longer exist" in text


def test_a_record_without_a_fingerprint_is_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)  # metadata and genome recorded without fingerprints
    cfg = Config.load(wd)
    cfg.stages["genome"].fingerprint = "abc"
    cfg.save(wd)
    findings = diagnose(wd)
    warned = {f.area for f in findings if f.level == "warn" and "resume fingerprint" in f.message}
    assert warned == {"metadata"}
    result = _runner.invoke(app, ["doctor", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "[WARN] metadata: completed without a resume fingerprint" in result.stdout


def test_tree2tax_child_under_two_parents_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text(
        "child\tparent\na\tA\na\tB\nb\tA\nc\tB\nA\troot\nB\troot\n",
        encoding="utf-8",
    )
    (wd / "genomes_map.tsv").write_text(_MAP, encoding="utf-8")
    assert "under two parents" in _messages(diagnose(wd), "fail")


def test_tree2tax_empty_node_name_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text("child\tparent\na\t\nb\troot\nc\troot\n", encoding="utf-8")
    (wd / "genomes_map.tsv").write_text(_MAP, encoding="utf-8")
    assert "empty node name" in _messages(diagnose(wd), "fail")


def test_tree2tax_two_top_nodes_is_a_failure(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text(
        "child\tparent\na\troot\nb\tother\nc\tother\n",
        encoding="utf-8",
    )
    (wd / "genomes_map.tsv").write_text(_MAP, encoding="utf-8")
    assert "2 top nodes" in _messages(diagnose(wd), "fail")


def test_tree2tax_top_node_not_root_is_a_warning(tmp_path: Path) -> None:
    wd = _base_workdir(tmp_path)
    (wd / "tree2tax.tsv").write_text(
        "child\tparent\na\ttop\nb\ttop\nc\ttop\n",
        encoding="utf-8",
    )
    (wd / "genomes_map.tsv").write_text(_MAP, encoding="utf-8")
    findings = diagnose(wd)
    assert "FlexTaxD validates" in _messages(findings, "warn")
    assert "tree2tax" not in _messages(findings, "fail")
