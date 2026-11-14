"""core.sketches with a fake sourmash: what is written, kept, replaced and removed."""

from __future__ import annotations

import json
import logging
import sqlite3
from pathlib import Path

import pytest

from repgenr.core import sketches
from repgenr.core.context import WorkdirContext
from repgenr.core.errors import MissingBinaryError, ToolExecutionError
from repgenr.core.inputs import file_digest, manifest_digest
from repgenr.core.manifest import SCHEMA_VERSION, GenomeRecord, Manifest, SketchRecord
from repgenr.core.sketches import (
    OUTGROUP_JSON,
    SKETCH_PARAMS,
    SketchSource,
    remove_stale,
    sketch_counts,
    sketch_genomes,
    sketch_stage_genomes,
    sketch_status,
    sketch_targets,
    sources_from_workdir,
)
from repgenr.core.sourmash import SOURMASH_IMAGE, SOURMASH_TOOL, sourmash_capabilities

_LOG = logging.getLogger("test")


def _workdir(root: Path, n: int = 3, *, outgroup: bool = True) -> WorkdirContext:
    ctx = WorkdirContext(root, create=True)
    ctx.genomes_dir.mkdir()
    records = []
    for i in range(n):
        name = f"Fam_Gen_sp_GCA_{i:06d}.1.fasta"
        (ctx.genomes_dir / name).write_text(f">g{i}\n{'ACGT' * (i + 5)}\n", encoding="utf-8")
        records.append(GenomeRecord(f"GCA_{i:06d}.1", name, "local", "Fam", "Gen", "sp"))
    if outgroup:
        ctx.outgroup_dir.mkdir()
        name = "Fam_Out_og_GCF_000009.1.fasta"
        (ctx.outgroup_dir / name).write_text(">o\nTTTTGGGG\n", encoding="utf-8")
        records.append(GenomeRecord("GCF_000009.1", name, "local", "Fam", "Out", "og", True))
    ctx.manifest.replace_genomes(records)
    return ctx


def _sketch_all(ctx: WorkdirContext, threads: int = 2, **kw) -> sketches.SketchSummary:
    return sketch_genomes(ctx, sketch_targets(ctx), threads, _LOG, **kw)


def test_shared_spec_is_used_by_every_sourmash_adapter() -> None:
    from repgenr.classifiers.sourmash import SourmashClassifier
    from repgenr.dereplicators.sourmash import SourmashDereplicator
    from repgenr.treebuilders.sourmash import SourmashBuilder

    for cls in (SourmashClassifier, SourmashDereplicator, SourmashBuilder):
        caps = cls.capabilities
        assert caps.container == SOURMASH_IMAGE == SOURMASH_TOOL.container
        assert caps.required_binaries == SOURMASH_TOOL.required_binaries
        assert caps.conda == SOURMASH_TOOL.conda and caps.name == "sourmash"
    with pytest.raises(ValueError, match="container"):
        sourmash_capabilities(container="other:1")


def test_parameters_and_command() -> None:
    assert SKETCH_PARAMS == "k=21,k=31,k=51,scaled=1000"
    cmd = sketches.sketch_command(Path("/g/x.fasta.gz"), "x", Path("/s/x.sig.zip"))
    assert [str(c) for c in cmd] == [
        "sourmash", "sketch", "dna", "-p", SKETCH_PARAMS,
        "--name", "x", "-o", "/s/x.sig.zip", "/g/x.fasta.gz",
    ]  # fmt: skip


def test_writes_every_missing_sketch_and_records_it(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd")
    summary = _sketch_all(ctx)
    assert summary.as_dict() == {
        "present": 0, "written": 4, "replaced": 0, "forced": 0, "copied": 0, "removed": 0,
    }  # fmt: skip
    # The hidden digest cache (.digests.json) is not a sketch.
    names = {p.name for p in (ctx.workdir / "sketches").iterdir() if not p.name.startswith(".")}
    assert names == {
        "Fam_Gen_sp_GCA_000000.1.sig.zip",
        "Fam_Gen_sp_GCA_000001.1.sig.zip",
        "Fam_Gen_sp_GCA_000002.1.sig.zip",
        "Fam_Out_og_GCF_000009.1.sig.zip",
    }
    # The signature is named after the genome's record name.
    assert sorted(fake_sourmash.calls) == sorted(n.removesuffix(".sig.zip") for n in names)
    records = ctx.manifest.sketch_records()
    rec = records["GCA_000001.1"]
    genome = ctx.genomes_dir / "Fam_Gen_sp_GCA_000001.1.fasta"
    assert rec == SketchRecord(
        "sketches/Fam_Gen_sp_GCA_000001.1.sig.zip", SKETCH_PARAMS, file_digest(genome)
    )
    # The outgroup is a manifest row: no outgroup.json.
    assert "GCF_000009.1" in records
    assert not (ctx.workdir / "sketches" / OUTGROUP_JSON).exists()
    assert sketch_counts(ctx.workdir) == (0, 0)  # no selection.tsv in this workdir


def test_up_to_date_sketches_are_left_alone(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd")
    _sketch_all(ctx)
    fake_sourmash.calls.clear()
    summary = _sketch_all(ctx)
    assert (summary.present, summary.written, summary.replaced) == (4, 0, 0)
    assert fake_sourmash.calls == []
    status = sketch_status(ctx)
    assert len(status.present) == 4 and not status.missing and not status.stale


def test_a_changed_fasta_or_parameter_string_makes_a_sketch_stale(tmp_path, fake_sourmash):
    ctx = _workdir(tmp_path / "wd")
    _sketch_all(ctx)
    genome = ctx.genomes_dir / "Fam_Gen_sp_GCA_000000.1.fasta"
    genome.write_text(">g0\nGGGGCCCC\n", encoding="utf-8")
    assert sketch_status(ctx).stale == ["Fam_Gen_sp_GCA_000000.1"]
    # Without hashing the change is not seen; the cheap form trusts the record.
    assert sketch_status(ctx, verify=False).stale == []
    ctx.manifest.set_sketches([("GCA_000001.1", SketchRecord("sketches/x", "k=31", "d"))])
    assert "Fam_Gen_sp_GCA_000001.1" in sketch_status(ctx, verify=False).stale
    fake_sourmash.calls.clear()
    summary = _sketch_all(ctx)
    assert (summary.present, summary.replaced) == (2, 2)
    assert sorted(fake_sourmash.calls) == ["Fam_Gen_sp_GCA_000000.1", "Fam_Gen_sp_GCA_000001.1"]
    assert ctx.manifest.sketch_records()["GCA_000000.1"].digest == file_digest(genome)


def test_force_sketches_every_genome_again(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd", outgroup=False)
    _sketch_all(ctx)
    summary = _sketch_all(ctx, force=True)
    assert (summary.present, summary.replaced, summary.forced) == (0, 0, 3)
    assert "3 replaced (--force)" in summary.line()


def test_orphans_are_removed_with_the_genome_set(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd")
    _sketch_all(ctx)
    directory = ctx.workdir / "sketches"
    (directory / "._Fam_Gen_sp_GCA_000000.1.sig.zip").write_text("x")  # exFAT companion
    (directory / ".Fam_Gen_sp_GCA_000002.1.partial.sig.zip").write_text("x")
    keep = [
        r for r in ctx.manifest.all_genomes(include_outgroup=True) if r.accession != "GCA_000002.1"
    ]
    ctx.manifest.replace_genomes(keep)
    (ctx.genomes_dir / "Fam_Gen_sp_GCA_000002.1.fasta").unlink()
    assert remove_stale(ctx, _LOG) == 1
    assert not (directory / "Fam_Gen_sp_GCA_000002.1.sig.zip").exists()
    assert not (directory / ".Fam_Gen_sp_GCA_000002.1.partial.sig.zip").exists()
    assert len(list(directory.glob("[!.]*.sig.zip"))) == 3
    # An emptied set leaves no sketches/ (dotfiles do not keep it).
    ctx.manifest.replace_genomes([])
    for f in [*ctx.genomes_dir.iterdir(), *ctx.outgroup_dir.iterdir()]:
        f.unlink()
    assert remove_stale(ctx, _LOG) == 3
    assert not directory.exists()


def test_a_failed_sketch_leaves_no_file_and_the_others_are_recorded(tmp_path, fake_sourmash):
    ctx = _workdir(tmp_path / "wd")
    import threading

    fake_sourmash.fail = {"Fam_Gen_sp_GCA_000001.1"}
    # Four parties for four genomes: no call returns, so none can fail, before
    # all four have started. A started genome cannot be cancelled, so the other
    # three are certain to finish and be recorded.
    fake_sourmash.barrier = threading.Barrier(4)
    with pytest.raises(ToolExecutionError):
        _sketch_all(ctx, threads=4)
    directory = ctx.workdir / "sketches"
    assert not (directory / "Fam_Gen_sp_GCA_000001.1.sig.zip").exists()
    assert not [p for p in directory.iterdir() if "partial" in p.name]
    assert set(ctx.manifest.sketch_records()) == {"GCA_000000.1", "GCA_000002.1", "GCF_000009.1"}
    # The rerun sketches only the one that failed.
    fake_sourmash.fail = set()
    fake_sourmash.barrier = None
    fake_sourmash.calls.clear()
    summary = _sketch_all(ctx)
    assert fake_sourmash.calls == ["Fam_Gen_sp_GCA_000001.1"]
    assert (summary.present, summary.written) == (3, 1)


def test_a_partial_file_of_a_killed_run_is_not_a_sketch(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd", n=1, outgroup=False)
    directory = ctx.workdir / "sketches"
    directory.mkdir()
    (directory / ".Fam_Gen_sp_GCA_000000.1.partial.sig.zip").write_text("half")
    assert sketch_status(ctx).missing == ["Fam_Gen_sp_GCA_000000.1"]
    assert _sketch_all(ctx).written == 1
    assert sorted(p.name for p in directory.iterdir()) == [
        sketches.DIGESTS_JSON,
        "Fam_Gen_sp_GCA_000000.1.sig.zip",
    ]


def test_concurrent_sketches_are_bounded_by_threads(tmp_path, fake_sourmash) -> None:
    import threading

    ctx = _workdir(tmp_path / "wd", n=8, outgroup=False)
    fake_sourmash.barrier = threading.Barrier(2)  # calls pair up: an overlap is certain
    _sketch_all(ctx, threads=3)
    assert 2 <= fake_sourmash.peak <= 3
    assert len(fake_sourmash.calls) == 8


def test_an_interrupt_keeps_the_records_of_finished_sketches(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd", outgroup=False)
    fake_sourmash.interrupt = {"Fam_Gen_sp_GCA_000001.1"}
    with pytest.raises(KeyboardInterrupt):
        _sketch_all(ctx, threads=1)
    # The first genome finished before the interrupt and is recorded.
    assert set(ctx.manifest.sketch_records()) == {"GCA_000000.1"}
    fake_sourmash.interrupt = set()
    fake_sourmash.calls.clear()
    summary = _sketch_all(ctx, threads=1)
    # The worker may have started the third genome before the interrupt was
    # seen; its unrecorded sketch counts as stale, never as present.
    assert summary.present == 1 and summary.written + summary.replaced == 2
    assert "Fam_Gen_sp_GCA_000000.1" not in fake_sourmash.calls


def test_remove_stale_leaves_directories_alone(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd", n=1, outgroup=False)
    _sketch_all(ctx)
    odd = ctx.workdir / "sketches" / "unknown.sig.zip"
    odd.mkdir()
    assert remove_stale(ctx, _LOG) == 0
    assert odd.is_dir()


def test_a_pruning_error_is_a_warning_by_default(tmp_path, fake_sourmash, monkeypatch, caplog):
    ctx = _workdir(tmp_path / "wd", n=1, outgroup=False)

    def unreadable(ctx, logger=None):
        raise PermissionError("sketches/ is read-only")

    monkeypatch.setattr(sketches, "remove_stale", unreadable)
    with caplog.at_level(logging.WARNING):
        summary, _ = sketch_stage_genomes(ctx, None, "genome", _LOG)
    assert summary["state"] == "done" and summary["removed"] == 0
    assert "read-only" in caplog.text
    with pytest.raises(PermissionError):
        sketch_stage_genomes(ctx, True, "genome", _LOG)


def test_an_outgroup_without_a_manifest_row_is_recorded_in_outgroup_json(tmp_path, fake_sourmash):
    ctx = _workdir(tmp_path / "wd")
    ctx.manifest.replace_genomes(
        [r for r in ctx.manifest.all_genomes(include_outgroup=True) if not r.is_outgroup]
    )
    targets = sketch_targets(ctx)
    assert [t.accession for t in targets][-1] is None
    _sketch_all(ctx)
    data = json.loads((ctx.workdir / "sketches" / OUTGROUP_JSON).read_text())
    assert set(data) == {"Fam_Out_og_GCF_000009.1.fasta"}
    assert data["Fam_Out_og_GCF_000009.1.fasta"]["params"] == SKETCH_PARAMS
    fake_sourmash.calls.clear()
    assert _sketch_all(ctx).present == 4 and fake_sourmash.calls == []
    # The outgroup leaves the set: its sketch and record go.
    (ctx.outgroup_dir / "Fam_Out_og_GCF_000009.1.fasta").unlink()
    assert remove_stale(ctx) == 1
    assert not (ctx.workdir / "sketches" / OUTGROUP_JSON).exists()


def test_a_matching_sketch_of_another_workdir_is_copied(tmp_path, fake_sourmash) -> None:
    src = _workdir(tmp_path / "src", outgroup=False)
    _sketch_all(src)
    src.close()
    reuse = sources_from_workdir(src.workdir)
    assert set(reuse) == {f"Fam_Gen_sp_GCA_{i:06d}.1" for i in range(3)}
    dst = _workdir(tmp_path / "dst", outgroup=False)
    # One genome differs from its namesake in the source: it is sketched.
    (dst.genomes_dir / "Fam_Gen_sp_GCA_000002.1.fasta").write_text(">x\nCCCC\n")
    fake_sourmash.calls.clear()
    summary = sketch_genomes(dst, sketch_targets(dst), 2, _LOG, reuse=reuse)
    assert (summary.copied, summary.written) == (2, 1)
    assert fake_sourmash.calls == ["Fam_Gen_sp_GCA_000002.1"]
    copied = dst.workdir / "sketches" / "Fam_Gen_sp_GCA_000000.1.sig.zip"
    assert copied.read_text() == (src.workdir / "sketches" / copied.name).read_text()
    assert (
        dst.manifest.sketch_records()["GCA_000000.1"].digest
        == reuse["Fam_Gen_sp_GCA_000000.1"].record.digest
    )
    # A source without sketches offers none.
    assert sources_from_workdir(tmp_path / "nowhere") == {}
    assert isinstance(next(iter(reuse.values())), SketchSource)


def test_stage_step_default_skips_without_sourmash(tmp_path, caplog) -> None:
    ctx = _workdir(tmp_path / "wd")
    with caplog.at_level(logging.INFO):
        summary, versions = sketch_stage_genomes(ctx, None, "genome", _LOG)
    assert summary == {"state": "unavailable", "removed": 0} and versions == {}
    assert "--sketch" in caplog.text and "sourmash not found" in caplog.text
    assert not (ctx.workdir / "sketches").exists()


def _no_sourmash(caps):
    raise MissingBinaryError("Required binary 'sourmash' not found on PATH.")


def test_stage_step_explicit_sketch_needs_sourmash(tmp_path, monkeypatch, caplog) -> None:
    ctx = _workdir(tmp_path / "wd", n=1, outgroup=False)
    monkeypatch.setattr(sketches, "preflight", _no_sourmash)
    with pytest.raises(MissingBinaryError):
        sketch_stage_genomes(ctx, True, "genome", _LOG)
    # The default tolerates a sourmash that fails its check (too old, say).
    monkeypatch.setattr(sketches, "tool_available", lambda caps: True)
    with caplog.at_level(logging.INFO):
        summary, _ = sketch_stage_genomes(ctx, None, "genome", _LOG)
    assert summary["state"] == "unavailable" and "not found on PATH" in caplog.text


def test_stage_step_with_no_sketch_only_prunes(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd")
    _sketch_all(ctx)
    ctx.manifest.replace_genomes(ctx.manifest.all_genomes(include_outgroup=True)[1:])
    summary, versions = sketch_stage_genomes(ctx, False, "ingest", _LOG)
    assert summary == {"state": "off", "removed": 1} and versions == {}


def test_stage_step_failure_is_a_warning_by_default(tmp_path, fake_sourmash, caplog) -> None:
    ctx = _workdir(tmp_path / "wd", n=2, outgroup=False)
    fake_sourmash.fail = {"Fam_Gen_sp_GCA_000000.1"}
    with caplog.at_level(logging.WARNING):
        summary, versions = sketch_stage_genomes(ctx, None, "genome", _LOG, threads=1)
    assert summary["state"] == "failed" and versions == {"sourmash": "4.9.4"}
    assert "repgenr sketch" in caplog.text
    with pytest.raises(ToolExecutionError):
        sketch_stage_genomes(ctx, True, "genome", _LOG, threads=1)


def test_stage_step_only_sketches_the_named_genomes(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd")
    only = [ctx.genomes_dir / "Fam_Gen_sp_GCA_000001.1.fasta"]
    summary, _ = sketch_stage_genomes(ctx, None, "assemble", _LOG, only=only)
    assert summary["written"] == 1 and fake_sourmash.calls == ["Fam_Gen_sp_GCA_000001.1"]


# --- manifest v4 ----------------------------------------------------------------

_V3_SCHEMA = (
    "CREATE TABLE genomes (accession TEXT PRIMARY KEY, filename TEXT, source TEXT, "
    "family TEXT, genus TEXT, species TEXT, is_outgroup INTEGER DEFAULT 0, "
    "derep_status TEXT, representative TEXT, completeness REAL, contamination REAL, "
    "gtdb_representative INTEGER DEFAULT 0); "
    "PRAGMA user_version=3;"
)


def _v3_database(path: Path) -> None:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.executescript(_V3_SCHEMA)
    conn.execute(
        "INSERT INTO genomes (accession, filename, family, genus, species, completeness, "
        "contamination, gtdb_representative, derep_status) "
        "VALUES ('A', 'a.fasta', 'F', 'G', 's', 99.0, 0.5, 1, 'representative')"
    )
    conn.commit()
    conn.close()


def test_v3_manifest_migrates_with_the_digest_unchanged(tmp_path) -> None:
    path = tmp_path / "m.sqlite"
    _v3_database(path)
    old = Manifest.open_readonly(path)
    before = manifest_digest(old)
    assert old.sketch_records() == {}  # a v3 file read without migration
    old.close()
    m = Manifest(path)
    assert SCHEMA_VERSION == 5
    assert int(m._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION
    cols = {r[1] for r in m._conn.execute("PRAGMA table_info(genomes)")}
    assert {"sketch_file", "sketch_params", "sketch_digest"} <= cols
    # v5 renames the flag column and keeps its values.
    assert "source_representative" in cols and "gtdb_representative" not in cols
    assert manifest_digest(m) == before
    m.set_sketches([("A", SketchRecord("sketches/a.sig.zip", SKETCH_PARAMS, "abc"))])
    # Sketch records are not part of the digest the stages fingerprint on.
    assert manifest_digest(m) == before
    (row,) = m.all_genomes()
    assert row.sketch_file == "sketches/a.sig.zip" and row.source_representative
    # A genome upsert keeps the sketch record.
    m.replace_genomes([GenomeRecord("A", "a.fasta", "gtdb", "F", "G", "s")])
    assert m.sketch_records()["A"].digest == "abc"
    m.set_sketches([("A", None)])
    assert m.sketch_records() == {}
    m.close()


def test_v4_migration_treats_a_duplicate_column_as_done(tmp_path) -> None:
    path = tmp_path / "m.sqlite"
    _v3_database(path)
    conn = sqlite3.connect(path)
    conn.execute("ALTER TABLE genomes ADD COLUMN sketch_file TEXT")  # another writer was first
    conn.commit()
    conn.close()
    m = Manifest(path)
    assert int(m._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION
    m.close()


def test_status_counts_from_the_selection(tmp_path, fake_sourmash) -> None:
    from repgenr.core.contracts import SELECTION_TSV, SelectionRow, write_selection

    ctx = _workdir(tmp_path / "wd")
    write_selection(
        ctx.workdir / SELECTION_TSV,
        [
            SelectionRow(r.accession, "Fam", "Gen", "sp", r.is_outgroup, r.filename)
            for r in ctx.manifest.all_genomes(include_outgroup=True)
        ],
    )
    _sketch_all(ctx)
    (ctx.workdir / "sketches" / "Fam_Gen_sp_GCA_000000.1.sig.zip").unlink()
    assert sketch_counts(ctx.workdir) == (3, 4)
    assert len(sketches.expected_sketch_files(ctx.workdir)) == 4


def test_an_outgroup_row_without_a_filename_is_matched_by_accession(tmp_path, fake_sourmash):
    """The bacterial genome stage records no filename on its outgroup row."""
    ctx = _workdir(tmp_path / "wd", n=1)
    rows = ctx.manifest.all_genomes(include_outgroup=True)
    for row in rows:
        if row.is_outgroup:
            row.filename = None
    ctx.manifest.replace_genomes(rows)
    targets = sketch_targets(ctx)
    assert [t.accession for t in targets] == ["GCA_000000.1", "GCF_000009.1"]
    _sketch_all(ctx)
    assert "GCF_000009.1" in ctx.manifest.sketch_records()
    assert not (ctx.workdir / "sketches" / OUTGROUP_JSON).exists()
    assert remove_stale(ctx) == 0
