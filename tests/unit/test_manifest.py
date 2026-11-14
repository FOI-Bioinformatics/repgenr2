"""Unit tests for the SQLite manifest, incl. batched writes."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from repgenr.core import manifest as manifest_module
from repgenr.core.errors import WorkdirError
from repgenr.core.manifest import (
    BUSY_TIMEOUT_MS,
    SCHEMA_VERSION,
    GenomeRecord,
    Manifest,
)


def _rec(i: int) -> GenomeRecord:
    return GenomeRecord(accession=f"GCF_{i:06d}.1", filename=f"f_{i}.fasta", source="gtdb")


def test_schema_version_and_busy_timeout(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    assert int(m._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION
    assert int(m._conn.execute("PRAGMA busy_timeout").fetchone()[0]) == BUSY_TIMEOUT_MS
    m.close()


def test_adopts_pre_versioning_db(tmp_path: Path) -> None:
    # A pre-versioning DB reports user_version=0; opening it adopts v1 and keeps data.
    p = tmp_path / "manifest.sqlite"
    m = Manifest(p)
    m.upsert_many([_rec(1)])
    m._conn.execute("PRAGMA user_version=0")  # simulate an old, unversioned DB
    m._conn.commit()
    m.close()
    m2 = Manifest(p)
    assert int(m2._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION
    assert m2.count() == 1
    m2.close()


def test_rejects_newer_schema(tmp_path: Path) -> None:
    p = tmp_path / "manifest.sqlite"
    m = Manifest(p)
    m._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION + 1}")
    m._conn.commit()
    m.close()
    with pytest.raises(WorkdirError, match="newer than this"):
        Manifest(p)


def test_readonly_rejects_newer_schema_without_leaking_the_connection(
    tmp_path: Path, monkeypatch
) -> None:
    """A version rejection in the read-only branch must still close the handle."""
    p = tmp_path / "manifest.sqlite"
    m = Manifest(p)
    m._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION + 1}")
    m._conn.commit()
    m.close()

    opened: list[sqlite3.Connection] = []
    real_connect = sqlite3.connect

    def tracking_connect(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        conn = real_connect(*args, **kwargs)
        opened.append(conn)
        return conn

    monkeypatch.setattr(manifest_module.sqlite3, "connect", tracking_connect)
    with pytest.raises(WorkdirError, match="newer than this"):
        Manifest.open_readonly(p)

    assert opened, "the read-only branch never opened a connection"
    for conn in opened:
        with pytest.raises(sqlite3.ProgrammingError):
            conn.execute("SELECT 1")  # a closed connection refuses further use


def test_wal_mode_enabled(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    mode = m._conn.execute("PRAGMA journal_mode").fetchone()[0]
    assert mode.lower() == "wal"
    m.close()


def test_upsert_many_batches_and_persists(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    m.upsert_many([_rec(i) for i in range(50)])
    assert m.count() == 50
    # upsert is an upsert: re-inserting updates, does not duplicate
    again = _rec(0)
    again.species = "updated"
    m.upsert_many([again])
    assert m.count() == 50
    rows = {g.accession: g for g in m.all_genomes()}
    assert rows["GCF_000000.1"].species == "updated"
    m.close()


def test_set_derep_status_many(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "manifest.sqlite")
    m.upsert_many([_rec(i) for i in range(5)])
    m.set_derep_status_many(
        [
            ("GCF_000000.1", "representative", None),
            ("GCF_000001.1", "contained", "GCF_000000.1"),
            ("GCF_999999.1", "contained", "GCF_000000.1"),  # absent -> no-op, not an error
        ]
    )
    reps = {g.accession for g in m.representatives()}
    assert reps == {"GCF_000000.1"}
    by_acc = {g.accession: g for g in m.all_genomes()}
    assert by_acc["GCF_000001.1"].derep_status == "contained"
    assert by_acc["GCF_000001.1"].representative == "GCF_000000.1"
    m.close()


def test_replace_genomes_deletes_deselected(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many(
        [
            GenomeRecord(accession="A1", species="a"),
            GenomeRecord(accession="A2", species="b"),
            GenomeRecord(accession="OG", species="og", is_outgroup=True),
        ]
    )
    # re-selection keeps A1, drops A2, new outgroup
    m.replace_genomes(
        [
            GenomeRecord(accession="A1", species="a"),
            GenomeRecord(accession="OG2", species="og2", is_outgroup=True),
        ]
    )
    accs = {g.accession for g in m.all_genomes(include_outgroup=True)}
    assert accs == {"A1", "OG2"}


def test_replace_genomes_keeps_derep_status_of_unchanged_genomes(tmp_path: Path) -> None:
    """A forced re-ingest of the same set used to clear every derep status, while
    dereplicate (whose digest ignores those columns) stayed up to date."""
    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many(
        [
            GenomeRecord(accession="A1", filename="a1.fasta"),
            GenomeRecord(accession="A2", filename="a2.fasta"),
        ]
    )
    m.set_derep_status_many([("A1", "representative", "A1"), ("A2", "contained", "A1")])
    m.replace_genomes(
        [
            GenomeRecord(accession="A1", filename="a1.fasta", species="x"),
            GenomeRecord(accession="A2", filename="a2_renamed.fasta"),
            GenomeRecord(accession="A3", filename="a3.fasta"),
        ]
    )
    rows = {g.accession: g for g in m.all_genomes()}
    assert (rows["A1"].derep_status, rows["A1"].representative) == ("representative", "A1")
    assert rows["A1"].species == "x"
    assert rows["A2"].derep_status is None  # a new file is not the genome dereplicated
    assert rows["A3"].derep_status is None


def test_set_derep_status_many_clears_stale_rows(tmp_path: Path) -> None:
    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many(
        [
            GenomeRecord(accession="A1"),
            GenomeRecord(accession="A2"),
            GenomeRecord(accession="A3"),
        ]
    )
    m.set_derep_status_many(
        [
            ("A1", "representative", None),
            ("A2", "contained", "A1"),
            ("A3", "representative", None),
        ]
    )
    # a narrower re-dereplication: A3 no longer appears in the result at all
    m.set_derep_status_many(
        [
            ("A1", "representative", None),
            ("A2", "contained", "A1"),
        ]
    )
    by_acc = {g.accession: g for g in m.all_genomes(include_outgroup=True)}
    assert by_acc["A3"].derep_status is None, "stale representative status cleared"
    assert by_acc["A3"].representative is None
    assert by_acc["A1"].derep_status == "representative"


def test_manifest_digest_changes_after_reconciliation(tmp_path: Path) -> None:
    from repgenr.core.inputs import manifest_digest

    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many([GenomeRecord(accession="A1"), GenomeRecord(accession="A2")])
    before = manifest_digest(m)
    m.replace_genomes([GenomeRecord(accession="A1")])
    assert manifest_digest(m) != before


def test_quality_round_trip_and_migration(tmp_path):
    from repgenr.core.manifest import SCHEMA_VERSION, GenomeRecord, Manifest

    m = Manifest(tmp_path / "m.sqlite")
    m.upsert(GenomeRecord(accession="A", filename="a.fasta", completeness=99.1, contamination=0.4))
    m.upsert(GenomeRecord(accession="B", filename="b.fasta"))
    assert m.quality() == {"a.fasta": (99.1, 0.4)}
    assert SCHEMA_VERSION == 5
    m.close()


def test_v1_database_is_migrated(tmp_path):
    import sqlite3

    from repgenr.core.manifest import Manifest

    path = tmp_path / "old.sqlite"
    conn = sqlite3.connect(path)
    conn.executescript(
        "CREATE TABLE genomes (accession TEXT PRIMARY KEY, filename TEXT, source TEXT, "
        "family TEXT, genus TEXT, species TEXT, is_outgroup INTEGER DEFAULT 0, "
        "derep_status TEXT, representative TEXT); PRAGMA user_version=1;"
    )
    conn.execute("INSERT INTO genomes (accession, filename) VALUES ('A', 'a.fasta')")
    conn.commit()
    conn.close()
    m = Manifest(path)
    assert m.quality() == {}
    assert m.all_genomes()[0].completeness is None
    m.close()


def test_quality_on_readonly_v1_database_returns_empty(tmp_path):
    """quality() must not raise on a v1 file opened read-only (never migrated)."""
    import sqlite3

    from repgenr.core.manifest import Manifest

    path = tmp_path / "old.sqlite"
    conn = sqlite3.connect(path)
    conn.executescript(
        "CREATE TABLE genomes (accession TEXT PRIMARY KEY, filename TEXT, source TEXT, "
        "family TEXT, genus TEXT, species TEXT, is_outgroup INTEGER DEFAULT 0, "
        "derep_status TEXT, representative TEXT); PRAGMA user_version=1;"
    )
    conn.execute("INSERT INTO genomes (accession, filename) VALUES ('A', 'a.fasta')")
    conn.commit()
    conn.close()

    m = Manifest.open_readonly(path)
    assert m.quality() == {}
    assert [g.accession for g in m.all_genomes()] == ["A"]
    m.close()

    # never migrated: user_version is still 1 on reopen
    reopened = sqlite3.connect(path)
    assert int(reopened.execute("PRAGMA user_version").fetchone()[0]) == 1
    reopened.close()
