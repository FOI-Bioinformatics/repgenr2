"""ingest --from-workdir: merge the genome sets of earlier working directories.

A GTDB working directory (metadata + genome) and a reads working directory
(reads + assemble) are built here from synthetic selection tables, genome files
and manifests, then combined into a third working directory.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli import base as cli
from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    SELECTION_TSV,
    SelectionRow,
    read_selection,
    write_selection,
)
from repgenr.core.errors import UserInputError
from repgenr.core.manifest import Manifest, record_from_selection
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.ingest import IngestParams, run

_runner = CliRunner()

GTDB_ROWS = [
    SelectionRow(
        "GCF_000001.1",
        "Francisellaceae",
        "Francisella",
        "tularensis",
        False,
        "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta",
        99.5,
        0.2,
        False,
    ),
    SelectionRow(
        "GCF_000002.1",
        "Francisellaceae",
        "Francisella",
        "tularensis",
        False,
        "Francisellaceae_Francisella_tularensis_GCF_000002.1.fasta",
        98.0,
        0.4,
        True,
    ),
    SelectionRow(
        "GCF_000099.1",
        "Francisellaceae",
        "Allofrancisella",
        "guangzhouensis",
        True,
        "Francisellaceae_Allofrancisella_guangzhouensis_GCF_000099.1.fasta",
    ),
]
SRA_ROWS = [
    SelectionRow(
        "SRR0000001",
        "Francisellaceae",
        "Francisella",
        "tularensis",
        False,
        "Francisellaceae_Francisella_tularensis_SRR0000001.fasta",
        97.1,
        1.3,
        False,
    ),
    SelectionRow(
        "SRR0000002",
        "Francisellaceae",
        "Francisella",
        "tularensis",
        False,
        "Francisellaceae_Francisella_tularensis_SRR0000002.fasta",
        96.4,
        0.9,
        False,
    ),
]


def _make_workdir(path: Path, rows: list[SelectionRow], source: str | None) -> Path:
    """A finished working directory: genomes/, outgroup/, selection.tsv and,
    when ``source`` is given, a manifest labelling every genome with it."""
    (path / "genomes").mkdir(parents=True)
    for i, row in enumerate(rows):
        where = path / ("outgroup" if row.is_outgroup else "genomes")
        where.mkdir(exist_ok=True)
        (where / row.filename).write_text(f">{row.accession}\n{'ACGT' * (10 + i)}\n")
    write_selection(path / SELECTION_TSV, rows)
    if source is not None:
        manifest = Manifest.open(path)
        try:
            manifest.replace_genomes([record_from_selection(r, source) for r in rows])
        finally:
            manifest.close()
    return path


@pytest.fixture
def sources(tmp_path: Path) -> tuple[Path, Path]:
    wd1 = _make_workdir(tmp_path / "wd1", GTDB_ROWS, "gtdb")
    wd2 = _make_workdir(tmp_path / "wd2", SRA_ROWS, "sra")
    return wd1, wd2


def _ctx(path: Path, caplog=None) -> WorkdirContext:
    ctx = WorkdirContext(path, create=True)
    if caplog is not None:
        ctx.logger.addHandler(caplog.handler)
    return ctx


def test_two_workdirs_merge_with_rows_and_sources(tmp_path, sources, caplog) -> None:
    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3", caplog)
    with caplog.at_level(logging.INFO):
        n = run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)]))

    assert n == 4
    rows = {r.accession: r for r in read_selection(ctx.workdir / SELECTION_TSV)}
    expected = [r for r in GTDB_ROWS + SRA_ROWS if not r.is_outgroup]
    assert rows == {r.accession: r for r in expected}, "rows carry over unchanged"
    staged = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert staged == sorted(r.filename for r in expected)
    assert all((ctx.genomes_dir / name).is_symlink() for name in staged)
    assert (ctx.genomes_dir / SRA_ROWS[0].filename).resolve() == (
        wd2 / "genomes" / SRA_ROWS[0].filename
    ).resolve()

    records = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert {acc: g.source for acc, g in records.items()} == {
        "GCF_000001.1": "gtdb",
        "GCF_000002.1": "gtdb",
        "SRR0000001": "sra",
        "SRR0000002": "sra",
    }
    assert records["SRR0000001"].completeness == pytest.approx(97.1)
    assert ctx.manifest.gtdb_representatives() == {GTDB_ROWS[1].filename}

    # The outgroup of wd1 is not carried over, and the log says so once.
    assert not (ctx.workdir / "outgroup_accession.txt").exists()
    assert not ctx.outgroup_dir.exists() or not any(ctx.outgroup_dir.iterdir())
    notes = [r.getMessage() for r in caplog.records if "is not carried over" in r.getMessage()]
    assert len(notes) == 1
    assert "GCF_000099.1" in notes[0] and str(wd1) in notes[0]

    params = ctx.config.stages["ingest"].params
    assert params["from_workdirs"] == [str(wd1), str(wd2)]
    assert params["genomes_dir"] is None
    assert params["sources"] == {"gtdb": 2, "sra": 2}


def test_outgroup_names_the_outgroup_of_a_source_workdir(tmp_path, sources, caplog) -> None:
    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3", caplog)
    with caplog.at_level(logging.INFO):
        run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)], outgroup="GCF_000099.1"))

    outgroup = GTDB_ROWS[2]
    assert (ctx.outgroup_dir / outgroup.filename).exists()
    assert (ctx.workdir / "outgroup_accession.txt").read_text().strip() == "GCF_000099.1"
    (marked,) = [r for r in read_selection(ctx.workdir / SELECTION_TSV) if r.is_outgroup]
    assert marked.accession == "GCF_000099.1"
    records = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert records["GCF_000099.1"].is_outgroup and records["GCF_000099.1"].source == "gtdb"
    assert not any("is not carried over" in r.getMessage() for r in caplog.records)


def test_outgroup_names_an_ingroup_genome_of_a_source_workdir(tmp_path, sources) -> None:
    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3")
    n = run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)], outgroup="SRR0000002"))

    assert n == 3
    assert not (ctx.genomes_dir / SRA_ROWS[1].filename).exists()
    assert (ctx.outgroup_dir / SRA_ROWS[1].filename).exists()
    records = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert records["SRR0000002"].is_outgroup and records["SRR0000002"].source == "sra"


def test_unknown_outgroup_names_both_kinds_of_source(tmp_path, sources) -> None:
    wd1, _ = sources
    with pytest.raises(UserInputError, match="--genomes-dir or --from-workdir"):
        run(_ctx(tmp_path / "wd3"), IngestParams(from_workdirs=[str(wd1)], outgroup="nope"))


def test_workdir_without_manifest_gives_local_source(tmp_path) -> None:
    wd = _make_workdir(tmp_path / "plain", SRA_ROWS, source=None)
    ctx = _ctx(tmp_path / "wd3")
    run(ctx, IngestParams(from_workdirs=[str(wd)]))
    assert {g.source for g in ctx.manifest.all_genomes()} == {"local"}


def test_genomes_dir_combined_with_from_workdir(tmp_path, sources) -> None:
    _, wd2 = sources
    local = tmp_path / "local"
    local.mkdir()
    (local / "Fam_Gen_sp_GCA_000777.1.fasta").write_text(">x\nACGTACGT\n")
    ctx = _ctx(tmp_path / "wd3")

    n = run(ctx, IngestParams(genomes_dir=str(local), from_workdirs=[str(wd2)]))

    assert n == 3
    sources_by_acc = {g.accession: g.source for g in ctx.manifest.all_genomes()}
    assert sources_by_acc == {"GCA_000777.1": "local", "SRR0000001": "sra", "SRR0000002": "sra"}
    assert ctx.config.stages["ingest"].params["sources"] == {"local": 1, "sra": 2}


def test_same_accession_in_two_sources_is_refused(tmp_path, sources) -> None:
    wd1, wd2 = sources
    clash = _make_workdir(tmp_path / "clash", [GTDB_ROWS[0]], "local")
    with pytest.raises(UserInputError) as err:
        run(_ctx(tmp_path / "wd3"), IngestParams(from_workdirs=[str(wd1), str(clash)]))
    message = str(err.value)
    assert "share an accession" in message and "GCF_000001.1" in message
    assert str(wd1) in message and str(clash) in message


def test_same_filename_with_other_accessions_is_refused(tmp_path, sources) -> None:
    _, wd2 = sources
    renamed = SelectionRow(
        "SRR9999999", "F", "G", "s", False, SRA_ROWS[0].filename, None, None, False
    )
    clash = _make_workdir(tmp_path / "clash", [renamed], None)
    with pytest.raises(UserInputError) as err:
        run(_ctx(tmp_path / "wd3"), IngestParams(from_workdirs=[str(wd2), str(clash)]))
    message = str(err.value)
    assert "share a filename" in message and SRA_ROWS[0].filename in message
    assert str(wd2) in message and str(clash) in message


def test_missing_genome_file_is_refused(tmp_path, sources) -> None:
    _, wd2 = sources
    (wd2 / "genomes" / SRA_ROWS[1].filename).unlink()
    with pytest.raises(UserInputError) as err:
        run(_ctx(tmp_path / "wd3"), IngestParams(from_workdirs=[str(wd2)]))
    assert str(wd2) in str(err.value) and SRA_ROWS[1].filename in str(err.value)


def test_workdir_without_selection_is_refused(tmp_path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(UserInputError, match="holds no selection.tsv or genomes/"):
        run(_ctx(tmp_path / "wd3"), IngestParams(from_workdirs=[str(empty)]))


def test_from_workdir_equal_to_target_is_refused(tmp_path, sources) -> None:
    wd1, _ = sources
    with pytest.raises(UserInputError, match="working directory being written"):
        run(WorkdirContext(wd1), IngestParams(from_workdirs=[str(wd1)]))


def test_repeated_from_workdir_is_refused(tmp_path, sources) -> None:
    wd1, _ = sources
    with pytest.raises(UserInputError, match="more than once"):
        run(_ctx(tmp_path / "wd3"), IngestParams(from_workdirs=[str(wd1), str(wd1) + "/"]))


def test_selection_without_genomes_dir_is_refused(tmp_path, sources) -> None:
    wd1, _ = sources
    with pytest.raises(UserInputError, match="--selection names genomes under --genomes-dir"):
        run(
            _ctx(tmp_path / "wd3"),
            IngestParams(from_workdirs=[str(wd1)], selection=str(wd1 / SELECTION_TSV)),
        )


def test_reingest_may_drop_sra_genomes_an_earlier_ingest_took(tmp_path, sources) -> None:
    """sra genomes that came in through --from-workdir are not appended genomes."""
    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3")
    run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)]))
    assert run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)])) == 4
    assert run(ctx, IngestParams(from_workdirs=[str(wd1)])) == 2
    assert {g.source for g in ctx.manifest.all_genomes()} == {"gtdb"}


def test_reingest_refuses_to_drop_sra_genomes_after_an_assemble(tmp_path, sources) -> None:
    """Once assemble has run in the workdir, an sra genome may have been appended
    there, so a re-ingest that would drop one is refused as before."""
    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3")
    run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)]))
    ctx.config.record_stage("assemble", params={"append": True}, completed="2026-10-08T00:00:00")
    ctx.save_config()
    assert run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)])) == 4
    with pytest.raises(UserInputError, match="appended from sequencing runs"):
        run(ctx, IngestParams(from_workdirs=[str(wd1)]))
    assert run(ctx, IngestParams(from_workdirs=[str(wd1)], drop_foreign=True)) == 2


# --- CLI: usage errors, resume and status ------------------------------------------


def test_cli_needs_a_source(tmp_path) -> None:
    result = _runner.invoke(app, ["ingest", "-wd", str(tmp_path / "wd3")])
    assert result.exit_code == 2
    assert "--genomes-dir, --from-workdir, or both" in result.output


def test_cli_selection_without_genomes_dir_is_a_usage_error(tmp_path, sources) -> None:
    wd1, _ = sources
    result = _runner.invoke(
        app,
        [
            "ingest",
            "-wd",
            str(tmp_path / "wd3"),
            "--from-workdir",
            str(wd1),
            "--selection",
            str(wd1 / SELECTION_TSV),
        ],
    )
    assert result.exit_code == 2
    assert not (tmp_path / "wd3" / SELECTION_TSV).exists()


def test_cli_resume_skips_then_reruns_on_source_selection_change(
    tmp_path, sources, monkeypatch
) -> None:
    wd1, wd2 = sources
    wd3 = tmp_path / "wd3"
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    argv = ["ingest", "-wd", str(wd3), "--from-workdir", str(wd1), "--from-workdir", str(wd2)]

    first = _runner.invoke(app, argv)
    assert first.exit_code == 0, first.output
    record = Config.load(wd3).stages["ingest"]
    stamp = record.completed
    assert any(key.endswith("wd2/selection.tsv") for key in record.inputs)
    assert any(key.endswith("wd1/genomes") for key in record.inputs)

    second = _runner.invoke(app, argv)
    assert second.exit_code == 0, second.output
    assert Config.load(wd3).stages["ingest"].completed == stamp, "unchanged sources must skip"

    status = _runner.invoke(app, ["status", "-wd", str(wd3)])
    assert status.exit_code == 0, status.output
    assert f"from workdirs: {wd1}, {wd2}" in status.output

    write_selection(wd2 / SELECTION_TSV, SRA_ROWS[:1])
    third = _runner.invoke(app, argv)
    assert third.exit_code == 0, third.output
    assert Config.load(wd3).stages["ingest"].completed != stamp
    assert not (wd3 / "genomes" / SRA_ROWS[1].filename).exists()
    assert {r.accession for r in read_selection(wd3 / SELECTION_TSV)} == {
        "GCF_000001.1",
        "GCF_000002.1",
        "SRR0000001",
    }


# --- dereplicate --keeper gtdb on the merged set -----------------------------------


class _FirstIsRep(Dereplicator):
    """One cluster; the tool's pick is the first genome in name order."""

    capabilities = ToolCapabilities(name="firstrep", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"firstrep": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = sorted(genomes, key=lambda g: g.name)
        rep, members = genomes[0], [g.name for g in genomes[1:]]
        return DerepResult(
            representatives=[rep],
            clusters={rep.name: members},
            genome_status={
                rep.name: STATUS_REPRESENTATIVE,
                **{m: STATUS_CONTAINED for m in members},
            },
        )


@pytest.fixture
def first_rep_tool():
    registry._load()
    registry.register("firstrep", _FirstIsRep, replace=True)
    yield
    registry._classes.pop("firstrep", None)


def test_keeper_gtdb_on_the_merged_workdir(tmp_path, sources, first_rep_tool) -> None:
    from repgenr.stages.dereplicate import DereplicateParams
    from repgenr.stages.dereplicate import run as dereplicate

    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3")
    run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)]))

    result = dereplicate(ctx, DereplicateParams(tool="firstrep", keeper="gtdb"))

    assert [r.name for r in result.representatives] == [GTDB_ROWS[1].filename]
    params = ctx.config.stages["dereplicate"].params
    assert params["keeper_effective"] == "gtdb"


# --- state of the source records ---------------------------------------------------


def _recorded_source(tmp_path: Path, name: str, monkeypatch) -> tuple[Path, Path]:
    """A source workdir written by `repgenr ingest --genomes-dir` (a finished
    ingest record with its resume inputs), and the directory it took."""
    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    raw = tmp_path / f"{name}_fasta"
    raw.mkdir()
    for i, row in enumerate(r for r in GTDB_ROWS if not r.is_outgroup):
        (raw / row.filename).write_text(f">{row.accession}\n{'ACGT' * (10 + i)}\n")
    wd = tmp_path / name
    result = _runner.invoke(app, ["ingest", "-wd", str(wd), "--genomes-dir", str(raw)])
    assert result.exit_code == 0, result.output
    assert not Config.load(wd).stages["ingest"].interrupted
    return wd, raw


def _source_warnings(caplog) -> list[str]:
    return [
        r.getMessage()
        for r in caplog.records
        if r.levelno == logging.WARNING and "--from-workdir" in r.getMessage()
    ]


def test_done_source_gives_no_warning(tmp_path, monkeypatch, caplog) -> None:
    src, _ = _recorded_source(tmp_path, "src", monkeypatch)
    ctx = _ctx(tmp_path / "wd3", caplog)
    with caplog.at_level(logging.INFO):
        run(ctx, IngestParams(from_workdirs=[str(src)], strict_sources=True))
    assert _source_warnings(caplog) == []
    params = ctx.config.stages["ingest"].params
    assert params["source_states"] == [{"path": str(src), "stage": "ingest", "state": "done"}]


def test_stale_source_warns_and_is_taken(tmp_path, monkeypatch, caplog) -> None:
    src, raw = _recorded_source(tmp_path, "src", monkeypatch)
    # A genome file the source ingest took changed after it finished.
    target = raw / GTDB_ROWS[0].filename
    target.write_text(target.read_text() + "ACGTACGT\n")
    ctx = _ctx(tmp_path / "wd3", caplog)
    with caplog.at_level(logging.INFO):
        n = run(ctx, IngestParams(from_workdirs=[str(src)]))
    assert n == 2
    warnings = _source_warnings(caplog)
    assert len(warnings) == 1
    assert str(src) in warnings[0] and "ingest, which is stale" in warnings[0]
    assert "input changed" in warnings[0]
    (state,) = ctx.config.stages["ingest"].params["source_states"]
    assert state == {"path": str(src), "stage": "ingest", "state": "stale"}


def test_interrupted_source_warns(tmp_path, monkeypatch, caplog) -> None:
    src, _ = _recorded_source(tmp_path, "src", monkeypatch)
    cfg = Config.load(src)
    cfg.stages["ingest"].completed = None  # the mark of a run that did not finish
    cfg.save(src)
    ctx = _ctx(tmp_path / "wd3", caplog)
    with caplog.at_level(logging.INFO):
        run(ctx, IngestParams(from_workdirs=[str(src)]))
    (warning,) = _source_warnings(caplog)
    assert "written by ingest, which is interrupted" in warning
    (state,) = ctx.config.stages["ingest"].params["source_states"]
    assert state["state"] == "interrupted"


def test_source_without_record_warns_and_is_taken(tmp_path, sources, caplog) -> None:
    wd1, wd2 = sources
    ctx = _ctx(tmp_path / "wd3", caplog)
    with caplog.at_level(logging.INFO):
        n = run(ctx, IngestParams(from_workdirs=[str(wd1), str(wd2)]))
    assert n == 4
    warnings = _source_warnings(caplog)
    assert len(warnings) == 2
    assert all("no stage record" in w for w in warnings)
    assert [s["state"] for s in ctx.config.stages["ingest"].params["source_states"]] == [
        "no stage record",
        "no stage record",
    ]


def test_latest_genome_set_record_is_evaluated(tmp_path, monkeypatch) -> None:
    from repgenr.stages.ingest import source_state

    src, _ = _recorded_source(tmp_path, "src", monkeypatch)
    cfg = Config.load(src)
    # An older genome record without inputs: the later ingest record is the one read.
    cfg.record_stage("genome", params={}, completed="2000-01-01T00:00:00+00:00")
    cfg.save(src)
    state = source_state(src)
    assert (state.stage, state.state) == ("ingest", "done")


def test_cli_strict_sources_refuses_before_writing(tmp_path, monkeypatch) -> None:
    src, raw = _recorded_source(tmp_path, "src", monkeypatch)
    target = raw / GTDB_ROWS[0].filename
    target.write_text(target.read_text() + "ACGTACGT\n")
    wd3 = tmp_path / "wd3"
    result = _runner.invoke(
        app, ["ingest", "-wd", str(wd3), "--from-workdir", str(src), "--strict-sources"]
    )
    assert result.exit_code == 2, result.output
    assert "--strict-sources" in result.output and "stale" in result.output
    assert not (wd3 / SELECTION_TSV).exists()
    assert not (wd3 / "genomes").exists()
    assert "ingest" not in Config.load(wd3).stages


def test_cli_strict_sources_refuses_a_source_without_record(tmp_path, sources) -> None:
    wd1, _ = sources
    wd3 = tmp_path / "wd3"
    result = _runner.invoke(
        app, ["ingest", "-wd", str(wd3), "--from-workdir", str(wd1), "--strict-sources"]
    )
    assert result.exit_code == 2, result.output
    assert "no stage record" in result.output
    assert not (wd3 / SELECTION_TSV).exists()


def test_status_names_a_source_not_done_when_ingested(tmp_path, monkeypatch, sources) -> None:
    src, _ = _recorded_source(tmp_path, "src", monkeypatch)
    _, wd2 = sources
    wd3 = tmp_path / "wd3"
    argv = ["ingest", "-wd", str(wd3), "--from-workdir", str(src), "--from-workdir", str(wd2)]
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, result.output
    status = _runner.invoke(app, ["status", "-wd", str(wd3)])
    assert status.exit_code == 0, status.output
    assert f"sources not done when ingested: {wd2} (no stage record)" in status.output
    assert f"{src} (" not in status.output, "a done source is not named"


def test_cli_repeat_with_strict_sources_skips(tmp_path, monkeypatch) -> None:
    src, _ = _recorded_source(tmp_path, "src", monkeypatch)
    wd3 = tmp_path / "wd3"
    argv = ["ingest", "-wd", str(wd3), "--from-workdir", str(src)]
    first = _runner.invoke(app, argv)
    assert first.exit_code == 0, first.output
    stamp = Config.load(wd3).stages["ingest"].completed
    second = _runner.invoke(app, [*argv, "--strict-sources"])
    assert second.exit_code == 0, second.output
    assert "skipping" in second.output
    assert Config.load(wd3).stages["ingest"].completed == stamp


def test_cli_strict_refusal_is_not_announced_as_a_rerun(tmp_path, monkeypatch) -> None:
    src, raw = _recorded_source(tmp_path, "src", monkeypatch)
    wd3 = tmp_path / "wd3"
    argv = ["ingest", "-wd", str(wd3), "--from-workdir", str(src)]
    assert _runner.invoke(app, argv).exit_code == 0
    # Changing a source FASTA makes both the source and the ingest of wd3 stale.
    target = raw / GTDB_ROWS[0].filename
    target.write_text(target.read_text() + "ACGTACGT\n")
    result = _runner.invoke(app, [*argv, "--strict-sources"])
    assert result.exit_code == 2, result.output
    assert "re-running" not in result.output
    assert not Config.load(wd3).stages["ingest"].interrupted


VIRAL_ROWS = [
    SelectionRow(
        "OZ389729.1",
        "Picornaviridae",
        "Hepatovirus",
        "Hepatovirus-ahepa",
        False,
        "Picornaviridae_Hepatovirus_Hepatovirus-ahepa_OZ389729.1.fasta",
    ),
]


def test_summary_names_viral_species_by_their_token(tmp_path, sources) -> None:
    # A viral token is the whole species name; a bacterial one the epithet.
    # Each genome of a combined set is named by the rule of its source.
    wd1, _ = sources
    viral = _make_workdir(tmp_path / "viral", VIRAL_ROWS, "ncbi_virus")
    ctx = _ctx(tmp_path / "wd3")
    run(ctx, IngestParams(from_workdirs=[str(wd1), str(viral)]))
    lines = (ctx.workdir / "metadata_summary.tsv").read_text(encoding="utf-8").splitlines()
    species = {line.split("\t")[1] for line in lines if line.startswith("species\t")}
    assert species == {"Francisella tularensis", "Hepatovirus-ahepa"}
