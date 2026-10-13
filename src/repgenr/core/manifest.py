"""SQLite genome manifest.

Replaces the ``str(dict)`` / ``pickle`` state blobs and the repeated
``os.listdir`` scans that do not scale to thousands of genomes. The manifest is
the single source of truth for which genomes are selected, where their files
live, their taxonomy, and their dereplication status.
"""

from __future__ import annotations

import os
import sqlite3
import time
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .errors import WorkdirError

MANIFEST_FILENAME = "manifest.sqlite"
SCHEMA_VERSION = 4  # bump + add a migration step when the table layout changes
BUSY_TIMEOUT_MS = 30000  # wait up to 30s for a competing writer before erroring

_UPSERT_SQL = """
    INSERT INTO genomes (accession, filename, source, family, genus,
                         species, is_outgroup, derep_status, representative,
                         completeness, contamination, gtdb_representative)
    VALUES (:accession, :filename, :source, :family, :genus,
            :species, :is_outgroup, :derep_status, :representative,
            :completeness, :contamination, :gtdb_representative)
    ON CONFLICT(accession) DO UPDATE SET
        filename=excluded.filename,
        source=excluded.source,
        family=excluded.family,
        genus=excluded.genus,
        species=excluded.species,
        is_outgroup=excluded.is_outgroup,
        derep_status=excluded.derep_status,
        representative=excluded.representative,
        completeness=excluded.completeness,
        contamination=excluded.contamination,
        gtdb_representative=excluded.gtdb_representative
"""

_SET_DEREP_SQL = "UPDATE genomes SET derep_status=?, representative=? WHERE accession=?"

# The sketch columns are written only by set_sketches: the genome upserts above
# leave them as they are, so re-selecting a genome keeps its sketch record
# (core.sketches decides from the FASTA digest whether it is still current).
_SKETCH_COLUMNS = ("sketch_file", "sketch_params", "sketch_digest")
_SET_SKETCH_SQL = (
    "UPDATE genomes SET sketch_file=?, sketch_params=?, sketch_digest=? WHERE accession=?"
)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS genomes (
    accession   TEXT PRIMARY KEY,
    filename    TEXT,
    source      TEXT,                  -- gtdb | bvbrc | ncbi
    family      TEXT,
    genus       TEXT,
    species     TEXT,
    is_outgroup INTEGER DEFAULT 0,
    derep_status TEXT,                 -- representative | contained | fail_qc | NULL
    representative TEXT,                -- accession of the representative, if contained
    completeness REAL,                 -- CheckM completeness percentage, if known
    contamination REAL,                -- CheckM contamination percentage, if known
    gtdb_representative INTEGER DEFAULT 0,  -- 1 for a GTDB species representative
    sketch_file TEXT,                  -- sketches/<name>.sig.zip, workdir-relative
    sketch_params TEXT,                -- sourmash parameter string of that sketch
    sketch_digest TEXT                 -- sha256 of the FASTA the sketch was built from
);
CREATE INDEX IF NOT EXISTS idx_genomes_species ON genomes(species);
CREATE INDEX IF NOT EXISTS idx_genomes_derep ON genomes(derep_status);
"""


@dataclass
class GenomeRecord:
    accession: str
    filename: str | None = None
    source: str | None = None
    family: str | None = None
    genus: str | None = None
    species: str | None = None
    is_outgroup: bool = False
    derep_status: str | None = None
    representative: str | None = None
    completeness: float | None = None
    contamination: float | None = None
    gtdb_representative: bool = False
    # Read from the manifest; written only through Manifest.set_sketches.
    sketch_file: str | None = None
    sketch_params: str | None = None
    sketch_digest: str | None = None


@dataclass(frozen=True)
class SketchRecord:
    """The recorded sketch of one genome (see core.sketches)."""

    file: str | None
    params: str | None
    digest: str | None


def record_from_selection(row: Any, source: str) -> GenomeRecord:
    """Build a manifest record from a ``SelectionRow`` (the contract hand-off).

    Shared by the stages that select genomes without the GTDB metadata path
    (ingest, the viral back-ends), so every entry path leaves the manifest
    describing the same genomes as ``selection.tsv``.
    """
    return GenomeRecord(
        accession=row.accession,
        filename=row.filename,
        source=source,
        family=row.family or None,
        genus=row.genus or None,
        species=row.species or None,
        is_outgroup=row.is_outgroup,
        completeness=row.completeness,
        contamination=row.contamination,
        gtdb_representative=bool(getattr(row, "gtdb_representative", False)),
    )


class Manifest:
    """Thin SQLite wrapper for the genome inventory.

    Not thread-safe: the single connection (``check_same_thread`` default) must be
    used from the thread that opened it. Stages write the manifest on the main
    thread; do not call it from inside a ``parallel_map`` worker. WAL + the busy
    timeout cover concurrency across *processes* (Nextflow scatter), not threads.
    """

    def __init__(self, path: str | os.PathLike[str], *, readonly: bool = False):
        self.path = Path(path)
        if readonly:
            if not self.path.exists():
                raise WorkdirError(f"Manifest not found: {self.path}")
            conn = None
            try:
                conn = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)
                conn.row_factory = sqlite3.Row
                conn.execute(f"PRAGMA busy_timeout={BUSY_TIMEOUT_MS}")
                self._conn = conn
                self._check_version()
            except WorkdirError:
                # A rejected schema version must not leak the open handle.
                if conn is not None:
                    conn.close()
                raise
            except sqlite3.OperationalError as exc:
                if conn is not None:
                    conn.close()
                # Every manifest is WAL-mode, and even a read-only connection
                # must create a "-shm" index file for it on first query -- so
                # a directory that is not writable (an archived or
                # permission-locked workdir, exactly where doctor is likely to
                # run), or a manifest file that is itself unreadable, surfaces
                # here as a raw sqlite error, not a missing file.
                raise WorkdirError(
                    f"Manifest at {self.path} cannot be opened read-only: {exc}. "
                    "A WAL-mode database needs a readable file and a writable "
                    "directory for its -shm file; copy the workdir or fix its "
                    "permissions."
                ) from exc
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.path)
        self._conn.row_factory = sqlite3.Row
        # Concurrent writers (parallel stages, two invocations on one workdir)
        # wait for the lock instead of failing immediately with "database is
        # locked". Set before the journal mode below, which needs the lock too.
        self._conn.execute(f"PRAGMA busy_timeout={BUSY_TIMEOUT_MS}")
        # WAL + synchronous=NORMAL: commits no longer pay a full fsync each, which
        # is the dominant cost for many small writes. The manifest is a workdir
        # artifact (regenerable from the stages), so the NORMAL durability
        # trade-off -- a power loss can lose only the last transaction -- is fine.
        _enable_wal(self._conn)
        self._conn.execute("PRAGMA synchronous=NORMAL")
        self._conn.executescript(_SCHEMA)
        self._migrate()
        self._conn.commit()

    def _check_version(self) -> int:
        """Reject a manifest written by a newer, incompatible RepGenR2.

        Returns the current ``user_version`` so callers that need it (e.g.
        ``_migrate``) do not have to issue a second ``PRAGMA`` query.
        """
        version = int(self._conn.execute("PRAGMA user_version").fetchone()[0])
        if version > SCHEMA_VERSION:
            raise WorkdirError(
                f"Manifest at {self.path} has schema version {version}, newer than this "
                f"RepGenR2 supports ({SCHEMA_VERSION}). Upgrade RepGenR2 or use a new workdir."
            )
        return version

    def _migrate(self) -> None:
        """Apply schema migrations keyed on ``PRAGMA user_version``.

        Existing pre-versioning databases report user_version=0; their layout
        already matches v1, so they are adopted as v1. Future schema changes add
        a numbered migration step and bump SCHEMA_VERSION.

        Several processes can open one old manifest at once (Nextflow scatter,
        two invocations on one workdir). The steps run under ``BEGIN
        IMMEDIATE``, which takes the write lock before the version and the
        columns are read, so a second process waits and then finds the
        migration done instead of adding a column twice.
        """
        if self._check_version() >= SCHEMA_VERSION:
            return
        self._conn.commit()  # BEGIN fails inside an open implicit transaction
        self._conn.execute("BEGIN IMMEDIATE")
        try:
            version = self._check_version()
            # (no v0->v1 data change: the CREATE IF NOT EXISTS schema is v1)
            if version < 2:
                for col in ("completeness", "contamination"):
                    _add_column(self._conn, col, "REAL")
            if version < 3:
                # v3 adds the GTDB species-representative flag. Rows of an
                # older database read 0 until the metadata stage runs again.
                _add_column(self._conn, "gtdb_representative", "INTEGER DEFAULT 0")
            if version < 4:
                # v4 records the genome sketch of each row (core.sketches).
                # Rows of an older database have none; the next sketch step
                # builds them.
                for col in _SKETCH_COLUMNS:
                    _add_column(self._conn, col, "TEXT")
            self._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION}")
            self._conn.commit()
        except BaseException:
            self._conn.rollback()
            raise

    @classmethod
    def open(cls, workdir: str | os.PathLike[str]) -> Manifest:
        return cls(Path(workdir) / MANIFEST_FILENAME)

    @classmethod
    def open_readonly(cls, path: str | os.PathLike[str]) -> Manifest:
        """Open an existing manifest without creating, migrating or journaling it.

        Raises ``WorkdirError`` if the file does not exist. Intended for
        read-only callers such as ``repgenr doctor``, which must not create,
        upgrade, or otherwise touch the workdir.

        Every manifest is WAL-mode, and SQLite must create (or open) a
        "-shm" index file for it even for a read-only connection's first
        query. If the manifest's directory is not writable (an archived or
        permission-locked workdir), or the manifest file itself is not
        readable, that fails -- and this raises ``WorkdirError`` rather than
        a raw ``sqlite3.OperationalError``. Copy the workdir or fix its
        permissions first.
        """
        return cls(path, readonly=True)

    def close(self) -> None:
        self._conn.close()

    def __enter__(self) -> Manifest:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    @contextmanager
    def transaction(self) -> Iterator[sqlite3.Connection]:
        try:
            yield self._conn
            self._conn.commit()
        except Exception:
            self._conn.rollback()
            raise

    def upsert(self, record: GenomeRecord) -> None:
        with self.transaction() as conn:
            conn.execute(_UPSERT_SQL, _record_params(record))

    def upsert_many(self, records: list[GenomeRecord]) -> None:
        # One transaction (one commit/fsync) for the whole batch -- committing per
        # record is ~0.5 ms each, i.e. seconds of fsync overhead at 1000s genomes.
        with self.transaction() as conn:
            conn.executemany(_UPSERT_SQL, [_record_params(r) for r in records])

    def replace_genomes(self, records: list[GenomeRecord]) -> None:
        """Make the manifest hold exactly ``records``: delete de-selected rows,
        then upsert, in one transaction. A re-selection (crashed or not) can no
        longer leave the manifest holding the union of old and new selections.

        A genome kept under the same accession and filename keeps its
        dereplication status when the new record carries none: a re-selection
        of an unchanged set leaves dereplicate up to date (its input digest
        ignores these columns), so it does not run again to restore them. A
        changed set re-runs dereplicate, which resets every status."""
        keep = {r.accession for r in records}
        with self.transaction() as conn:
            cur = conn.execute(
                "SELECT accession, filename, derep_status, representative FROM genomes"
            )
            existing = {row[0]: row for row in cur.fetchall()}
            stale = [acc for acc in existing if acc not in keep]
            conn.executemany("DELETE FROM genomes WHERE accession = ?", [(acc,) for acc in stale])
            params = []
            for record in records:
                item = _record_params(record)
                prior = existing.get(record.accession)
                if (
                    prior is not None
                    and prior[1] == record.filename
                    and record.derep_status is None
                    and record.representative is None
                ):
                    item["derep_status"], item["representative"] = prior[2], prior[3]
                params.append(item)
            conn.executemany(_UPSERT_SQL, params)

    def set_derep_status(
        self, accession: str, status: str, representative: str | None = None
    ) -> None:
        with self.transaction() as conn:
            conn.execute(_SET_DEREP_SQL, (status, representative, accession))

    def set_derep_status_many(self, updates: Sequence[tuple[str, str, str | None]]) -> None:
        """Batch derep-status updates (accession, status, representative) in one
        transaction; avoids one commit per genome on large sets.

        Statuses are RESET for every row first: a dereplication result
        describes the whole current genome set, so genomes absent from the new
        result must not keep a stale 'representative' status from a prior run.
        """
        rows = [(status, rep, accession) for accession, status, rep in updates]
        with self.transaction() as conn:
            conn.execute("UPDATE genomes SET derep_status = NULL, representative = NULL")
            conn.executemany(_SET_DEREP_SQL, rows)

    def count(self) -> int:
        cur = self._conn.execute("SELECT COUNT(*) AS n FROM genomes WHERE is_outgroup=0")
        return int(cur.fetchone()["n"])

    def representatives(self) -> list[GenomeRecord]:
        cur = self._conn.execute("SELECT * FROM genomes WHERE derep_status='representative'")
        return [_row_to_record(row) for row in cur.fetchall()]

    def all_genomes(self, include_outgroup: bool = False) -> list[GenomeRecord]:
        query = "SELECT * FROM genomes"
        if not include_outgroup:
            query += " WHERE is_outgroup=0"
        cur = self._conn.execute(query)
        return [_row_to_record(row) for row in cur.fetchall()]

    def quality(self) -> dict[str, tuple[float, float]]:
        """filename -> (completeness, contamination) for genomes with both values.

        A manifest opened read-only never runs ``_migrate``, so a pre-v2 file
        genuinely lacks these columns; querying them would raise
        ``sqlite3.OperationalError``. Check for the columns first and return an
        empty mapping rather than fail, mirroring ``_row_to_record``'s tolerance.
        """
        cols = {r[1] for r in self._conn.execute("PRAGMA table_info(genomes)")}
        if "completeness" not in cols or "contamination" not in cols:
            return {}
        rows = self._conn.execute(
            "SELECT filename, completeness, contamination FROM genomes "
            "WHERE filename IS NOT NULL AND completeness IS NOT NULL AND contamination IS NOT NULL"
        )
        return {r["filename"]: (float(r["completeness"]), float(r["contamination"])) for r in rows}

    def gtdb_representatives(self) -> set[str]:
        """Filenames of the genomes flagged as GTDB species representatives.

        A pre-v3 manifest opened read-only lacks the column; it holds no flag.
        """
        cols = {r[1] for r in self._conn.execute("PRAGMA table_info(genomes)")}
        if "gtdb_representative" not in cols:
            return set()
        rows = self._conn.execute(
            "SELECT filename FROM genomes WHERE gtdb_representative=1 AND filename IS NOT NULL"
        )
        return {r["filename"] for r in rows}

    def sketch_records(self) -> dict[str, SketchRecord]:
        """accession -> recorded sketch, for rows that have one.

        A pre-v4 manifest opened read-only lacks the columns; it records none.
        """
        cols = {r[1] for r in self._conn.execute("PRAGMA table_info(genomes)")}
        if not set(_SKETCH_COLUMNS) <= cols:
            return {}
        rows = self._conn.execute(
            "SELECT accession, sketch_file, sketch_params, sketch_digest FROM genomes "
            "WHERE sketch_file IS NOT NULL"
        )
        return {
            r["accession"]: SketchRecord(r["sketch_file"], r["sketch_params"], r["sketch_digest"])
            for r in rows
        }

    def set_sketches(self, updates: Sequence[tuple[str, SketchRecord | None]]) -> None:
        """Record (or, with None, clear) the sketch of each accession, in one transaction."""
        rows = [
            (
                rec.file if rec else None,
                rec.params if rec else None,
                rec.digest if rec else None,
                accession,
            )
            for accession, rec in updates
        ]
        with self.transaction() as conn:
            conn.executemany(_SET_SKETCH_SQL, rows)


def _enable_wal(conn: sqlite3.Connection) -> None:
    """Switch the database to WAL, retrying while another connection holds it.

    A file already in WAL mode (every manifest RepGenR2 wrote) needs no lock
    for this. A rollback-journal file (a copied or hand-made manifest) needs an
    exclusive lock, and SQLite reports "database is locked" at once rather
    than calling the busy handler when other connections hold a shared lock,
    so the switch is retried for up to ``BUSY_TIMEOUT_MS``.
    """
    deadline = time.monotonic() + BUSY_TIMEOUT_MS / 1000
    delay = 0.01
    while True:
        try:
            conn.execute("PRAGMA journal_mode=WAL")
            return
        except sqlite3.OperationalError as exc:
            if "locked" not in str(exc) or time.monotonic() >= deadline:
                raise
            time.sleep(delay)
            delay = min(delay * 2, 0.5)


def _add_column(conn: sqlite3.Connection, name: str, decl: str) -> None:
    """Add a column to ``genomes`` unless it is there already.

    The columns are read inside the caller's write transaction; a duplicate
    reported by SQLite (another writer added it first) also counts as done.
    """
    cols = {r[1] for r in conn.execute("PRAGMA table_info(genomes)")}
    if name in cols:
        return
    try:
        conn.execute(f"ALTER TABLE genomes ADD COLUMN {name} {decl}")
    except sqlite3.OperationalError as exc:
        if "duplicate column name" not in str(exc):
            raise


def _record_params(record: GenomeRecord) -> dict:
    return {
        "accession": record.accession,
        "filename": record.filename,
        "source": record.source,
        "family": record.family,
        "genus": record.genus,
        "species": record.species,
        "is_outgroup": int(record.is_outgroup),
        "derep_status": record.derep_status,
        "representative": record.representative,
        "completeness": record.completeness,
        "contamination": record.contamination,
        "gtdb_representative": int(record.gtdb_representative),
    }


def _row_to_record(row: sqlite3.Row) -> GenomeRecord:
    keys = row.keys()
    return GenomeRecord(
        accession=row["accession"],
        filename=row["filename"],
        source=row["source"],
        family=row["family"],
        genus=row["genus"],
        species=row["species"],
        is_outgroup=bool(row["is_outgroup"]),
        derep_status=row["derep_status"],
        representative=row["representative"],
        completeness=row["completeness"] if "completeness" in keys else None,
        contamination=row["contamination"] if "contamination" in keys else None,
        gtdb_representative=bool(row["gtdb_representative"])
        if "gtdb_representative" in keys
        else False,
        sketch_file=row["sketch_file"] if "sketch_file" in keys else None,
        sketch_params=row["sketch_params"] if "sketch_params" in keys else None,
        sketch_digest=row["sketch_digest"] if "sketch_digest" in keys else None,
    )
