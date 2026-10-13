"""--keeper gtdb: a GTDB species representative is kept as its cluster's keeper.

Covers the rule itself, the gtdb_representative column of selection.tsv, the
manifest column and its migration, the cluster summary column, the taxonomy
reduction preference, and the metadata stage that fills the flag.
"""

from __future__ import annotations

import logging
import sqlite3
from pathlib import Path

import pytest

from repgenr.core.contracts import (
    ClusterSummaryRow,
    SelectionRow,
    read_cluster_summary,
    read_selection,
    write_cluster_summary,
    write_selection,
)
from repgenr.core.errors import WorkdirError
from repgenr.core.inputs import manifest_digest
from repgenr.core.manifest import (
    SCHEMA_VERSION,
    GenomeRecord,
    Manifest,
    record_from_selection,
)
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)
from repgenr.stages.cluster_summary import summarise_clusters
from repgenr.stages.derep_keeper import prefer_gtdb_representatives

_LOG = logging.getLogger("keeper-gtdb")
_ALL = ["rep.fasta", "m1.fasta", "m2.fasta", "solo.fasta"]


def _result() -> DerepResult:
    return DerepResult(
        representatives=[Path("/g/rep.fasta"), Path("/g/solo.fasta")],
        clusters={"rep.fasta": ["m1.fasta", "m2.fasta"], "solo.fasta": []},
        genome_status={
            "rep.fasta": STATUS_REPRESENTATIVE,
            "solo.fasta": STATUS_REPRESENTATIVE,
            "m1.fasta": STATUS_CONTAINED,
            "m2.fasta": STATUS_CONTAINED,
        },
    )


# --- the rule ---------------------------------------------------------------


def test_gtdb_representative_is_kept_over_a_better_scored_member() -> None:
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2), "m2.fasta": (95.0, 1.0)}
    out, swaps = prefer_gtdb_representatives(_result(), {"m2.fasta"}, quality, _LOG)
    assert swaps == 1
    assert sorted(p.name for p in out.representatives) == ["m2.fasta", "solo.fasta"]
    assert out.clusters["m2.fasta"] == ["m1.fasta", "rep.fasta"]
    assert out.genome_status["m2.fasta"] == STATUS_REPRESENTATIVE
    assert out.genome_status["m1.fasta"] == STATUS_CONTAINED
    assert out.genome_status["rep.fasta"] == STATUS_CONTAINED
    assert out.representatives[0].parent == Path("/g")
    check_result_complete(out, _ALL)


def test_gtdb_representative_without_quality_is_still_kept() -> None:
    quality = {"rep.fasta": (99.0, 0.0), "m1.fasta": (99.0, 0.0)}
    out, swaps = prefer_gtdb_representatives(_result(), {"m2.fasta"}, quality, _LOG)
    assert swaps == 1
    assert "m2.fasta" in out.clusters


def test_tool_pick_that_is_a_gtdb_representative_stays() -> None:
    quality = {"m1.fasta": (99.0, 0.0)}
    out, swaps = prefer_gtdb_representatives(_result(), {"rep.fasta"}, quality, _LOG)
    assert swaps == 0
    assert out.clusters["rep.fasta"] == ["m1.fasta", "m2.fasta"]


def test_cluster_without_gtdb_representative_falls_back_to_quality() -> None:
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2)}
    out, swaps = prefer_gtdb_representatives(_result(), {"solo.fasta"}, quality, _LOG)
    assert swaps == 1
    assert sorted(out.clusters) == ["m1.fasta", "solo.fasta"]
    check_result_complete(out, _ALL)


def test_cluster_without_gtdb_representative_or_quality_keeps_tool_pick() -> None:
    out, swaps = prefer_gtdb_representatives(_result(), set(), {}, _LOG)
    assert swaps == 0
    assert sorted(out.clusters) == ["rep.fasta", "solo.fasta"]


def test_several_gtdb_representatives_keep_the_best_scored_and_warn(caplog) -> None:
    quality = {"m1.fasta": (95.0, 1.0), "m2.fasta": (99.0, 0.1)}
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        out, swaps = prefer_gtdb_representatives(_result(), {"m1.fasta", "m2.fasta"}, quality, _LOG)
    assert swaps == 1
    assert sorted(out.clusters) == ["m2.fasta", "solo.fasta"]
    assert out.genome_status["m1.fasta"] == STATUS_CONTAINED
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "m2.fasta" in warnings[0] and "m1.fasta" in warnings[0]
    check_result_complete(out, _ALL)


def test_several_unscored_gtdb_representatives_prefer_the_tool_pick() -> None:
    out, swaps = prefer_gtdb_representatives(
        _result(), {"rep.fasta", "m2.fasta", "m1.fasta"}, {}, _LOG
    )
    assert swaps == 0
    assert sorted(out.clusters) == ["rep.fasta", "solo.fasta"]


def test_several_unscored_gtdb_members_choose_by_name() -> None:
    out, _ = prefer_gtdb_representatives(_result(), {"m2.fasta", "m1.fasta"}, {}, _LOG)
    assert sorted(out.clusters) == ["m1.fasta", "solo.fasta"]


def test_high_quality_gtdb_representative_replaces_an_unscored_one() -> None:
    quality = {"m2.fasta": (99.0, 1.0)}
    out, _ = prefer_gtdb_representatives(_result(), {"rep.fasta", "m2.fasta"}, quality, _LOG)
    assert sorted(out.clusters) == ["m2.fasta", "solo.fasta"]


def test_unscored_gtdb_representative_stays_against_a_low_quality_one() -> None:
    """The quality rule's guard applies among GTDB representatives too: an
    unscored incumbent is replaced only by a high-quality genome."""
    quality = {"m2.fasta": (50.0, 9.0)}
    out, _ = prefer_gtdb_representatives(_result(), {"rep.fasta", "m2.fasta"}, quality, _LOG)
    assert sorted(out.clusters) == ["rep.fasta", "solo.fasta"]


def test_several_gtdb_representatives_tie_break_on_n50() -> None:
    quality = {"m1.fasta": (99.0, 0.0), "m2.fasta": (99.0, 0.0)}
    n50 = {"m1.fasta": 50_000, "m2.fasta": 2_000_000}.get
    out, _ = prefer_gtdb_representatives(_result(), {"m1.fasta", "m2.fasta"}, quality, _LOG, n50)
    assert sorted(out.clusters) == ["m2.fasta", "solo.fasta"]


def test_gtdb_and_quality_agree_when_every_genome_is_flagged() -> None:
    from repgenr.stages.derep_keeper import rescore_representatives

    quality = {"rep.fasta": (95.0, 1.0), "m1.fasta": (99.0, 0.5), "m2.fasta": (99.0, 0.5)}
    n50 = {"m1.fasta": 10_000, "m2.fasta": 900_000, "rep.fasta": 5_000_000}.get
    by_quality, _ = rescore_representatives(_result(), quality, _LOG, n50)
    by_gtdb, _ = prefer_gtdb_representatives(_result(), set(_ALL), quality, _LOG, n50)
    assert by_gtdb.clusters == by_quality.clusters


# --- selection.tsv ----------------------------------------------------------


def _row(name: str, *, rep: bool, outgroup: bool = False) -> SelectionRow:
    return SelectionRow(
        accession=name.split(".")[0],
        family="Fam",
        genus="g",
        species="s",
        is_outgroup=outgroup,
        filename=name,
        completeness=99.0,
        contamination=0.5,
        gtdb_representative=rep,
    )


def test_selection_round_trip_carries_the_flag(tmp_path) -> None:
    path = tmp_path / "selection.tsv"
    rows = [_row("A.fasta", rep=True), _row("B.fasta", rep=False, outgroup=True)]
    write_selection(path, rows)
    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert header[-1] == "gtdb_representative"
    assert path.read_text(encoding="utf-8").splitlines()[1].split("\t")[-1] == "1"
    assert read_selection(path) == rows


def test_selection_without_the_column_reads_as_not_representative(tmp_path) -> None:
    path = tmp_path / "selection.tsv"
    path.write_text("accession\tfilename\tis_outgroup\nA\tA.fasta\t0\n", encoding="utf-8")
    assert read_selection(path)[0].gtdb_representative is False


@pytest.mark.parametrize(("value", "expected"), [("true", True), ("yes", True), ("", False)])
def test_selection_flag_accepts_the_outgroup_spellings(tmp_path, value, expected) -> None:
    path = tmp_path / "selection.tsv"
    path.write_text(
        f"accession\tfilename\tgtdb_representative\nA\tA.fasta\t{value}\n", encoding="utf-8"
    )
    assert read_selection(path)[0].gtdb_representative is expected


def test_selection_flag_rejects_an_unknown_value(tmp_path) -> None:
    path = tmp_path / "selection.tsv"
    path.write_text("accession\tfilename\tgtdb_representative\nA\tA.fasta\tmaybe\n", "utf-8")
    with pytest.raises(WorkdirError, match="gtdb_representative"):
        read_selection(path)


def test_record_from_selection_carries_the_flag() -> None:
    assert record_from_selection(_row("A.fasta", rep=True), "gtdb").gtdb_representative
    assert not record_from_selection(_row("A.fasta", rep=False), "local").gtdb_representative


# --- manifest ---------------------------------------------------------------

_V2_SCHEMA = (
    "CREATE TABLE genomes (accession TEXT PRIMARY KEY, filename TEXT, source TEXT, "
    "family TEXT, genus TEXT, species TEXT, is_outgroup INTEGER DEFAULT 0, "
    "derep_status TEXT, representative TEXT, completeness REAL, contamination REAL); "
    "PRAGMA user_version=2;"
)


def _v2_database(path: Path) -> None:
    conn = sqlite3.connect(path)
    conn.executescript(_V2_SCHEMA)
    conn.execute(
        "INSERT INTO genomes (accession, filename, completeness, contamination) "
        "VALUES ('A', 'a.fasta', 99.0, 0.5)"
    )
    conn.commit()
    conn.close()


def test_manifest_round_trip_and_lookup(tmp_path) -> None:
    m = Manifest(tmp_path / "m.sqlite")
    m.upsert_many(
        [
            GenomeRecord(accession="A", filename="a.fasta", gtdb_representative=True),
            GenomeRecord(accession="B", filename="b.fasta"),
        ]
    )
    by_acc = {g.accession: g for g in m.all_genomes()}
    assert by_acc["A"].gtdb_representative and not by_acc["B"].gtdb_representative
    assert m.gtdb_representatives() == {"a.fasta"}
    m.close()


def test_v2_manifest_is_migrated_with_the_flag_unset(tmp_path) -> None:
    path = tmp_path / "old.sqlite"
    _v2_database(path)
    m = Manifest(path)
    assert int(m._conn.execute("PRAGMA user_version").fetchone()[0]) == SCHEMA_VERSION == 4
    (record,) = m.all_genomes()
    assert record.gtdb_representative is False and record.completeness == 99.0
    assert m.gtdb_representatives() == set()
    m.upsert(GenomeRecord(accession="A", filename="a.fasta", gtdb_representative=True))
    assert m.gtdb_representatives() == {"a.fasta"}
    m.close()


def test_readonly_v2_manifest_reports_no_representatives(tmp_path) -> None:
    path = tmp_path / "old.sqlite"
    _v2_database(path)
    m = Manifest.open_readonly(path)
    assert m.gtdb_representatives() == set()
    assert m.all_genomes()[0].gtdb_representative is False
    m.close()


def test_manifest_digest_is_unchanged_for_rows_without_the_flag(tmp_path) -> None:
    """Migrated workdirs (every flag 0) keep their dereplicate fingerprint;
    a set flag is a new input and changes it."""
    import hashlib

    m = Manifest(tmp_path / "m.sqlite")
    m.upsert(GenomeRecord(accession="A", filename="a.fasta", species="s"))
    expected = hashlib.sha256(
        ("\0".join(["A", "a.fasta", "", "", "s", "False", "", "", "", ""]) + "\n").encode()
    ).hexdigest()
    assert manifest_digest(m) == expected
    m.upsert(GenomeRecord(accession="A", filename="a.fasta", species="s", gtdb_representative=True))
    assert manifest_digest(m) != expected
    m.close()


# --- cluster summary --------------------------------------------------------


def test_cluster_summary_marks_a_gtdb_representative_keeper(tmp_path) -> None:
    clusters = {"rep.fasta": ["m1.fasta"], "solo.fasta": []}
    rows = summarise_clusters(clusters, {}, {}, gtdb_representatives={"rep.fasta", "m1.fasta"})
    by_rep = {r.representative: r for r in rows}
    assert by_rep["rep.fasta"].rep_is_gtdb_representative is True
    assert by_rep["solo.fasta"].rep_is_gtdb_representative is False
    path = tmp_path / "cluster_summary.tsv"
    write_cluster_summary(path, rows)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t")[-1] == "rep_is_gtdb_representative"
    assert [line.split("\t")[-1] for line in lines[1:]] == ["1", "0"]
    assert read_cluster_summary(path) == rows


def test_cluster_summary_without_the_column_reads_as_false(tmp_path) -> None:
    path = tmp_path / "cluster_summary.tsv"
    path.write_text("representative\tn_members\tn_species\tspecies\nr.fasta\t0\t1\ts\n", "utf-8")
    assert read_cluster_summary(path) == [
        ClusterSummaryRow(representative="r.fasta", n_members=0, n_species=1, species="s")
    ]


# --- taxonomy reduction -----------------------------------------------------


def test_reduce_prefers_a_gtdb_representative() -> None:
    from repgenr.stages.dereplicate import _reduce_by_taxonomy

    result = DerepResult(
        representatives=[Path("/g/a.fasta"), Path("/g/b.fasta")],
        clusters={"a.fasta": ["a1.fasta"], "b.fasta": []},
        genome_status={
            "a.fasta": STATUS_REPRESENTATIVE,
            "b.fasta": STATUS_REPRESENTATIVE,
            "a1.fasta": STATUS_CONTAINED,
        },
    )
    quality = {"a.fasta": (99.0, 0.0), "b.fasta": (90.0, 1.0)}
    taxon_of = {"a.fasta": "sp", "b.fasta": "sp"}
    plain = _reduce_by_taxonomy(result, "species", quality, _LOG, taxon_of=taxon_of)
    assert [p.name for p in plain.representatives] == ["a.fasta"]
    out = _reduce_by_taxonomy(
        result, "species", quality, _LOG, taxon_of=taxon_of, prefer={"b.fasta"}
    )
    assert [p.name for p in out.representatives] == ["b.fasta"]
    assert sorted(out.clusters["b.fasta"]) == ["a.fasta", "a1.fasta"]


# --- command line -----------------------------------------------------------


def test_param_builder_accepts_keeper_gtdb() -> None:
    from repgenr.cli.param_builders import dereplicate_params

    assert dereplicate_params(tool="skder", keeper="gtdb").keeper == "gtdb"


def test_keeper_help_names_gtdb() -> None:
    from repgenr.cli.base import HELP_KEEPER

    assert "gtdb" in HELP_KEEPER


def _open_after_barrier(path: str, barrier, errors) -> None:  # noqa: ANN001
    try:
        barrier.wait()
        Manifest(path).close()
    except Exception as exc:  # noqa: BLE001 - reported to the parent
        errors.put(repr(exc))


@pytest.mark.parametrize("wal", [True, False])
@pytest.mark.parametrize("old_schema", ["v1", "v2", "v3"])
def test_concurrent_opens_of_an_old_manifest_migrate_once(tmp_path, old_schema, wal) -> None:
    """Two processes opening an old manifest at once (Nextflow scatter, two
    invocations on one workdir) both succeed; neither sees a duplicate column."""
    import multiprocessing as mp

    ctx = mp.get_context("spawn")
    schema = _V2_SCHEMA
    if old_schema == "v1":
        schema = schema.replace(", completeness REAL, contamination REAL", "").replace(
            "user_version=2", "user_version=1"
        )
    elif old_schema == "v3":
        schema = schema.replace(
            "contamination REAL)", "contamination REAL, gtdb_representative INTEGER DEFAULT 0)"
        ).replace("user_version=2", "user_version=3")
    for trial in range(6):
        path = tmp_path / f"m{trial}.sqlite"
        conn = sqlite3.connect(path)
        if wal:  # as every manifest this program writes; a copied or hand-made one may not be
            conn.execute("PRAGMA journal_mode=WAL")
        conn.executescript(schema)
        conn.close()
        barrier = ctx.Barrier(4)
        errors = ctx.Queue()
        procs = [
            ctx.Process(target=_open_after_barrier, args=(str(path), barrier, errors))
            for _ in range(4)
        ]
        for p in procs:
            p.start()
        for p in procs:
            p.join(60)
        found = []
        while not errors.empty():
            found.append(errors.get())
        assert found == []
        m = Manifest(path)
        cols = {r[1] for r in m._conn.execute("PRAGMA table_info(genomes)")}
        assert {"completeness", "contamination", "gtdb_representative"} <= cols
        assert {"sketch_file", "sketch_params", "sketch_digest"} <= cols
        m.close()


def test_migration_holds_the_write_lock_between_check_and_alter(tmp_path, monkeypatch) -> None:
    """A second writer that tries to add the column between this connection's
    column check and its ALTER must find the database locked, not win the race
    and leave this connection with a duplicate-column error."""
    import repgenr.core.manifest as manifest_mod

    path = tmp_path / "old.sqlite"
    _v2_database(path)
    racer_results: list[str] = []
    real_connect = sqlite3.connect

    class RacingConnection(sqlite3.Connection):
        def execute(self, sql, *args):  # noqa: ANN001, ANN201
            if "ADD COLUMN gtdb_representative" in sql and not racer_results:
                other = real_connect(path, timeout=0)
                try:
                    other.execute(
                        "ALTER TABLE genomes ADD COLUMN gtdb_representative INTEGER DEFAULT 0"
                    )
                    other.commit()
                    racer_results.append("racer added the column")
                except sqlite3.OperationalError as exc:
                    racer_results.append(f"racer blocked: {exc}")
                finally:
                    other.close()
            return super().execute(sql, *args)

    monkeypatch.setattr(
        manifest_mod.sqlite3,
        "connect",
        lambda *a, **k: real_connect(*a, **{**k, "factory": RacingConnection}),
    )
    m = Manifest(path)
    m.close()
    assert racer_results and racer_results[0].startswith("racer blocked")


def test_migration_treats_a_duplicate_column_as_done(tmp_path) -> None:
    """A column added by another writer whose user_version bump is not yet
    visible is not an error."""
    from repgenr.core.manifest import _add_column

    path = tmp_path / "m.sqlite"
    _v2_database(path)
    conn = sqlite3.connect(path)
    _add_column(conn, "gtdb_representative", "INTEGER DEFAULT 0")
    _add_column(conn, "gtdb_representative", "INTEGER DEFAULT 0")
    conn.close()


def test_gtdb_rule_reads_no_n50_when_one_representative_competes() -> None:
    """One flagged genome per cluster and no quality competition: no genome
    file needs to be read for its N50."""
    reads: list[str] = []

    def n50(name: str) -> int | None:
        reads.append(name)
        return 1_000_000

    # m2 is scored but not high quality, so the quality pass keeps the
    # unscored tool pick without ranking; the GTDB pass then keeps m2.
    quality = {"m2.fasta": (50.0, 9.0)}
    out, _ = prefer_gtdb_representatives(_result(), {"m2.fasta"}, quality, _LOG, n50)
    assert "m2.fasta" in out.clusters
    assert reads == []
