"""reads-screen: judge the sketches of reads-sketch against GTDB before any download."""

from __future__ import annotations

import json
import logging
import os
import threading
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli import cmd_reads
from repgenr.cli.main import app
from repgenr.core.config import Config
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    READS_SKETCHES_DIR,
    READS_SKETCHES_TSV,
    READS_TSV,
    ReadRow,
    ReadsSketchRow,
    write_reads,
    write_reads_sketches,
)
from repgenr.core.errors import ToolExecutionError, UserInputError, WorkdirError
from repgenr.core.sketches import READS_SKETCH_PARAMS
from repgenr.stages import assemble_screen as screen
from repgenr.stages import reads_screen as stage
from repgenr.stages.reads_screen import (
    PASSED_TXT,
    READS_SCREEN_DIR,
    READS_SCREEN_TSV,
    ReadsScreenParams,
    describe_screens,
    precheck,
    resolve_variant,
    run,
    select_runs,
)

_runner = CliRunner()

_LINEAGE = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales;"
    "f__Francisellaceae;g__Francisella;s__Francisella tularensis"
)
_FOREIGN = (
    "d__Bacteria;p__Bacillota;c__Bacilli;o__Mycoplasmatales;f__Mycoplasmataceae;"
    "g__Mycoplasma;s__Mycoplasma mycoides"
)
_TAX_HEADER = (
    "query_name,status,rank,fraction,lineage,query_md5,query_filename,"
    "f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank\n"
)


class FakeGather:
    """``sourmash gather`` and ``tax genome`` from canned hits by run; the
    query sketch is a text file whose first line is the run accession."""

    def __init__(self) -> None:
        self.hits: dict[str, tuple[str, float] | None] = {}
        self.gathers: list[str] = []
        self.fail: set[str] = set()
        self._lock = threading.Lock()

    def run_tool(self, caps, command, *, logger, **kwargs) -> int:  # noqa: ANN001
        argv = [str(c) for c in command]
        if argv[:2] == ["sourmash", "gather"]:
            run = Path(argv[2]).read_text(encoding="utf-8").splitlines()[0]
            with self._lock:
                self.gathers.append(run)
            if run in self.fail:
                raise ToolExecutionError(argv, 1, output="gather died", tool="sourmash")
            out = Path(argv[argv.index("-o") + 1])
            hit = self.hits.get(run)
            lines = ["intersect_bp,f_unique_weighted,name,query_name"]
            if hit is not None:
                lines.append(f"500000,{hit[1]},GCF_000000001.1 {hit[0].split('s__')[-1]},{run}")
            out.write_text("\n".join(lines) + "\n", encoding="utf-8")
            return 0
        assert argv[:3] == ["sourmash", "tax", "genome"]
        csv_list = Path(argv[argv.index("--from-file") + 1])
        gather_csv = Path(csv_list.read_text(encoding="utf-8").split()[0])
        run = gather_csv.read_text(encoding="utf-8").splitlines()[1].split(",")[-1]
        lineage, fraction = self.hits[run]  # type: ignore[misc]
        base = Path(argv[argv.index("--output-base") + 1])
        Path(str(base) + ".classifications.csv").write_text(
            _TAX_HEADER + f'{run},match,species,0.04,"{lineage}",x,y,{fraction},500000,0.9\n',
            encoding="utf-8",
        )
        return 0


class AniTable(dict):
    """ANI estimates by (query run, reference run); unlisted pairs are 0.9."""

    def __init__(self) -> None:
        super().__init__()
        self.compared: list[tuple[str, str]] = []


@pytest.fixture
def gather(monkeypatch) -> FakeGather:
    from repgenr.classifiers import sourmash as sm

    fake = FakeGather()
    monkeypatch.setattr(screen, "run_tool", fake.run_tool)
    monkeypatch.setattr(sm, "run_tool", fake.run_tool)
    monkeypatch.setattr(stage, "preflight", lambda caps: {"sourmash": "4.9.4"})
    return fake


@pytest.fixture
def ani(monkeypatch) -> AniTable:
    table = AniTable()

    def load(path: Path, ksize: int = 31):  # noqa: ANN202
        return path.read_text(encoding="utf-8").splitlines()[0]

    def containment(query, reference, ksize=31):  # noqa: ANN001, ANN202
        table.compared.append((query, reference))
        return table.get((query, reference), 0.9)

    monkeypatch.setattr(screen, "load_reads_hashes", load)
    monkeypatch.setattr(screen, "containment_ani", containment)
    return table


@pytest.fixture
def gtdb(tmp_path: Path, monkeypatch) -> tuple[str, str]:
    db = tmp_path / "db" / "gtdb-rs226-reps.k31.sig.zip"
    db.parent.mkdir()
    db.write_bytes(b"zip")
    lineages = tmp_path / "db" / "lineages.csv"
    lineages.write_text("ident,superkingdom\n", encoding="utf-8")
    monkeypatch.delenv(screen.GTDB_SKETCH_ENV, raising=False)
    monkeypatch.delenv(screen.GTDB_LINEAGES_ENV, raising=False)
    return str(db), str(lineages)


def _row(acc: str, biosample: str | None = None, genus: str = "Francisella") -> ReadRow:
    return ReadRow(
        run_accession=acc,
        biosample=biosample or f"SAMN{acc}",
        bioproject="PRJNA1",
        organism="Francisella tularensis",
        taxid="263",
        platform="ILLUMINA",
        instrument_model="x",
        layout="PAIRED",
        bases=1000,
        read_count=10,
        family="Francisellaceae",
        genus=genus,
        species="tularensis",
    )


def _workdir(tmp_path: Path, rows: list[ReadRow]) -> WorkdirContext:
    ctx = WorkdirContext(tmp_path / "wd", create=True, logger=logging.getLogger("test"))
    write_reads(ctx.workdir / READS_TSV, rows)
    ctx.config.record_stage("reads", params={}, completed="2026-10-19T00:00:00+00:00")
    ctx.save_config()
    return ctx


def _sketch(workdir: Path, variant: str, runs: list[str], *, failed: list[str] = ()) -> Path:
    """A reads-sketch variant: a fake sketch per run (its accession as text) and the index."""
    out = workdir / READS_SKETCHES_DIR / variant
    out.mkdir(parents=True, exist_ok=True)
    spots = None if variant == "all" else int(variant.split("_")[1])
    rows = []
    for acc in runs:
        (out / f"{acc}.sig.zip").write_text(f"{acc}\n", encoding="utf-8")
        rows.append(
            ReadsSketchRow(acc, "ok", spots, 10, 20, f"{acc}.sig.zip", READS_SKETCH_PARAMS, 1.0, "")
        )
    for acc in failed:
        rows.append(ReadsSketchRow(acc, "failed", spots, reason="fastq-dump exited 3: not found"))
    write_reads_sketches(out / READS_SKETCHES_TSV, rows)
    return out


def _params(gtdb: tuple[str, str], **over) -> ReadsScreenParams:  # noqa: ANN003
    return ReadsScreenParams(
        **{"gtdb_sketch": gtdb[0], "gtdb_lineages": gtdb[1], "jobs": 2, **over}
    )


def _table(workdir: Path, variant: str) -> dict[str, dict[str, str]]:
    path = workdir / READS_SCREEN_DIR / variant / READS_SCREEN_TSV
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return {
        cells[0]: dict(zip(header, cells, strict=True))
        for cells in (line.split("\t") for line in lines[1:])
    }


def _passed(workdir: Path, variant: str) -> list[str]:
    path = workdir / READS_SCREEN_DIR / variant / PASSED_TXT
    return path.read_text(encoding="utf-8").split()


# --- selection and the variant ------------------------------------------------------


def test_select_runs_refuses_an_accession_outside_reads_tsv() -> None:
    assert [r.run_accession for r in select_runs([_row("SRR1")], ReadsScreenParams())] == ["SRR1"]
    with pytest.raises(UserInputError, match="Not in reads.tsv: SRR9"):
        select_runs([_row("SRR1")], ReadsScreenParams(accessions=["SRR9"]))


def test_resolve_variant_takes_all_then_the_only_one_and_needs_a_name_otherwise(
    tmp_path: Path,
) -> None:
    sketches = tmp_path / READS_SKETCHES_DIR
    with pytest.raises(WorkdirError, match="No reads sketches"):
        resolve_variant(sketches, None)
    _sketch(tmp_path, "spots_5000", ["SRR1"])
    assert resolve_variant(sketches, None) == "spots_5000"
    assert resolve_variant(sketches, 5000) == "spots_5000"
    with pytest.raises(
        WorkdirError, match="No reads-sketch variant spots_100 .*present: spots_5000"
    ):
        resolve_variant(sketches, 100)
    _sketch(tmp_path, "spots_100", ["SRR1"])
    with pytest.raises(WorkdirError, match="Several .*spots_100, spots_5000.*--max-spots N"):
        resolve_variant(sketches, None)
    _sketch(tmp_path, "all", ["SRR1"])
    assert resolve_variant(sketches, None) == "all"


def test_precheck_refuses_before_anything_is_written(tmp_path: Path, gtdb, gather) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True, logger=logging.getLogger("test"))
    with pytest.raises(WorkdirError, match="completed reads stage"):
        precheck(ctx, _params(gtdb))
    ctx.close()
    ctx = _workdir(tmp_path, [_row("SRR1")])
    with pytest.raises(WorkdirError, match="No reads sketches"):
        precheck(ctx, _params(gtdb))
    _sketch(ctx.workdir, "all", ["SRR1"])
    with pytest.raises(UserInputError, match="Not in reads.tsv: SRR2"):
        precheck(ctx, _params(gtdb, accessions=["SRR2"]))
    with pytest.raises(UserInputError, match="reads-screen needs the GTDB sketch"):
        precheck(ctx, ReadsScreenParams())
    with pytest.raises(UserInputError, match="--gtdb-sketch .* does not exist"):
        precheck(ctx, ReadsScreenParams(gtdb_sketch="/nope", gtdb_lineages=gtdb[1]))
    with pytest.raises(UserInputError, match="--min-fraction must be within"):
        precheck(ctx, _params(gtdb, min_fraction=2.0))
    precheck(ctx, _params(gtdb))
    assert not (ctx.workdir / READS_SCREEN_DIR).exists()


# --- the screen ---------------------------------------------------------------------


def test_run_applies_the_three_checks_and_writes_table_decisions_and_passed(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    rows = [
        _row("SRR1", "SAMN1"),
        _row("SRR2", "SAMN1"),  # a duplicate of SRR1
        _row("SRR3"),  # another genus
        _row("SRR4"),  # a poor match
        _row("SRR5"),  # no GTDB match
        _row("SRR6"),  # sketch failed in reads-sketch
        _row("SRR7"),  # never selected by reads-sketch
    ]
    ctx = _workdir(tmp_path, rows)
    _sketch(ctx.workdir, "spots_5000", ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"], failed=["SRR6"])
    gather.hits = {
        "SRR1": (_LINEAGE, 0.5),
        "SRR2": (_LINEAGE, 0.5),
        "SRR3": (_FOREIGN, 0.6),
        "SRR4": (_LINEAGE, 0.1),
        "SRR5": None,
    }
    ani[("SRR2", "SRR1")] = 0.9995
    assert run(ctx, _params(gtdb)) == 1
    table = _table(ctx.workdir, "spots_5000")
    assert list(table) == [f"SRR{i}" for i in range(1, 8)], "every run of reads.tsv, in order"
    assert table["SRR1"]["decision"] == "pass" and table["SRR1"]["top_genus"] == "Francisella"
    assert table["SRR2"]["reason"].startswith("duplicate_isolate: contained in SRR1")
    assert table["SRR3"]["reason"].startswith("taxon_mismatch: submitted genus Francisella")
    assert table["SRR4"]["reason"].startswith("low_match_fraction: top species fraction 0.100")
    assert table["SRR5"]["reason"].startswith("taxon_mismatch: no GTDB match")
    assert (
        table["SRR6"]["reason"] == "screen_failed: no reads sketch (fastq-dump exited 3: not found)"
    )
    assert table["SRR7"]["reason"] == "screen_failed: no reads sketch (not sketched)"
    assert _passed(ctx.workdir, "spots_5000") == ["SRR1"]
    assert sorted(gather.gathers) == ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"]
    out = ctx.workdir / READS_SCREEN_DIR / "spots_5000"
    assert (out / "SRR1" / "screen.json").is_file() and (
        out / "SRR3" / "screen" / "gather.csv"
    ).is_file()
    assert not (out / "SRR6").exists() and not (out / "SRR7").exists()
    record = ctx.config.stages["reads_screen"]
    assert record.tool == "sourmash" and record.params["variant"] == "spots_5000"
    assert (record.params["selected"], record.params["screened"], record.params["reused"]) == (
        7,
        5,
        0,
    )
    counts = {k: record.params[k] for k in ("passed", "excused", "failed", "unsketched")}
    assert counts == {"passed": 1, "excused": 4, "failed": 0, "unsketched": 2}
    assert record.params["settings"]["gtdb_sketch"] == str(Path(gtdb[0]).resolve())
    assert describe_screens(ctx.workdir / READS_SCREEN_DIR) == (
        "spots_5000 (1 passed, 4 excused of 7 runs)"
    )


def test_the_duplicate_check_compares_only_the_same_biosample_in_reads_order(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1", "SAMN1"), _row("SRR2", "SAMN2"), _row("SRR3", "SAMN1")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2", "SRR3"])
    gather.hits = {r: (_LINEAGE, 0.5) for r in ("SRR1", "SRR2", "SRR3")}
    ani[("SRR2", "SRR1")] = 1.0  # another biosample: never compared
    ani[("SRR3", "SRR1")] = 1.0
    assert run(ctx, _params(gtdb, jobs=3)) == 2
    assert ani.compared == [("SRR3", "SRR1")]
    assert _passed(ctx.workdir, "all") == ["SRR1", "SRR2"]


def test_a_failed_gather_keeps_no_decision_and_is_screened_again(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    gather.fail = {"SRR2"}
    assert run(ctx, _params(gtdb)) == 1
    table = _table(ctx.workdir, "all")
    assert table["SRR2"]["reason"].startswith("screen_failed: gather:")
    assert not (ctx.workdir / READS_SCREEN_DIR / "all" / "SRR2" / "screen.json").exists()
    assert ctx.config.stages["reads_screen"].params["failed"] == 1
    gather.fail = set()
    assert run(ctx, _params(gtdb)) == 2
    assert sorted(gather.gathers[:2]) == ["SRR1", "SRR2"] and gather.gathers[2:] == ["SRR2"]
    assert _passed(ctx.workdir, "all") == ["SRR1", "SRR2"]


def test_a_repeat_reuses_decisions_and_a_changed_threshold_or_sketch_screens_again(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1", "SAMN1"), _row("SRR2", "SAMN1")])
    out = _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.35)}
    assert run(ctx, _params(gtdb)) == 2
    assert run(ctx, _params(gtdb)) == 2
    assert sorted(gather.gathers) == ["SRR1", "SRR2"], "nothing gathered again"
    assert ctx.config.stages["reads_screen"].params["reused"] == 2
    # A stricter fraction screens both again from the kept sketches.
    assert run(ctx, _params(gtdb, min_fraction=0.4)) == 1
    assert sorted(gather.gathers[2:]) == ["SRR1", "SRR2"] and len(gather.gathers) == 4
    assert _table(ctx.workdir, "all")["SRR2"]["decision"] == "excused"
    # Only the sketch of SRR1 is newer than its decision: it alone is screened again.
    decision = out.parent.parent / READS_SCREEN_DIR / "all" / "SRR1" / "screen.json"
    os.utime(out / "SRR1.sig.zip", ns=(decision.stat().st_mtime_ns + 10**9,) * 2)
    assert run(ctx, _params(gtdb, min_fraction=0.4)) == 1
    assert gather.gathers[-1] == "SRR1" and len(gather.gathers) == 5
    # --force screens every selected run again.
    ctx.force = True
    run(ctx, _params(gtdb, min_fraction=0.4))
    assert len(gather.gathers) == 7


def test_a_subset_is_screened_and_the_rest_of_the_table_is_kept(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2"), _row("SRR3")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2", "SRR3"])
    gather.hits = {r: (_LINEAGE, 0.5) for r in ("SRR1", "SRR2", "SRR3")}
    listing = tmp_path / "picked.txt"
    listing.write_text("SRR3\n", encoding="utf-8")
    assert run(ctx, _params(gtdb, accessions=["SRR1"], accession_file=str(listing))) == 2
    assert list(_table(ctx.workdir, "all")) == ["SRR1", "SRR3"]
    assert run(ctx, _params(gtdb, accessions=["SRR2"])) == 1
    assert list(_table(ctx.workdir, "all")) == ["SRR1", "SRR2", "SRR3"], "reads.tsv order"
    assert _passed(ctx.workdir, "all") == ["SRR1", "SRR2", "SRR3"]


def test_variants_are_screened_apart(tmp_path: Path, gtdb, gather, ani) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1")])
    _sketch(ctx.workdir, "spots_100", ["SRR1"])
    _sketch(ctx.workdir, "all", ["SRR1"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5)}
    run(ctx, _params(gtdb, max_spots=100))
    run(ctx, _params(gtdb))
    assert ctx.config.stages["reads_screen"].params["variant"] == "all"
    assert describe_screens(ctx.workdir / READS_SCREEN_DIR) == (
        "all (1 passed, 0 excused of 1 runs), spots_100 (1 passed, 0 excused of 1 runs)"
    )


# --- the command and the harness ------------------------------------------------------


def test_reads_screen_command_builds_its_params(monkeypatch, tmp_path: Path) -> None:
    seen: dict = {}

    def fake_run(stage_name, workdir, build, **kw):  # noqa: ANN001, ANN003, ANN202
        seen["stage"], seen["params"] = stage_name, build()

    monkeypatch.setattr(cmd_reads, "_run", fake_run)
    listing = tmp_path / "picked.txt"
    listing.write_text("SRR2\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "reads-screen",
            "-wd",
            str(tmp_path / "wd"),
            "--accession",
            "SRR1",
            "--accession-file",
            str(listing),
            "--max-spots",
            "5000",
            "--min-fraction",
            "0.4",
            "--dup-ani",
            "0.99",
            "--gtdb-sketch",
            "db.sig.zip",
            "--gtdb-lineages",
            "lin.csv",
            "--jobs",
            "3",
        ],
    )
    assert result.exit_code == 0, result.output
    assert seen["stage"] == "reads_screen"
    assert seen["params"] == ReadsScreenParams(
        accessions=["SRR1"],
        accession_file=str(listing),
        max_spots=5000,
        min_fraction=0.4,
        dup_ani=0.99,
        gtdb_sketch="db.sig.zip",
        gtdb_lineages="lin.csv",
        jobs=3,
    )


def test_harness_skips_a_repeat_and_screens_a_failed_run_again(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    _sketch(ctx.workdir, "spots_5000", ["SRR1", "SRR2"])
    ctx.close()
    wd = str(ctx.workdir)
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    gather.fail = {"SRR2"}
    args = ["reads-screen", "-wd", wd, "--gtdb-sketch", gtdb[0], "--gtdb-lineages", gtdb[1]]
    first = _runner.invoke(app, args)
    assert first.exit_code == 0, first.output
    assert "1 run(s) could not be screened (screen_failed: gather)" in first.output
    # The failed run keeps no decision: the repeat runs the stage for it alone.
    gather.fail = set()
    again = _runner.invoke(app, args)
    assert again.exit_code == 0, again.output
    assert "already completed" not in again.output
    assert sorted(gather.gathers[:2]) == ["SRR1", "SRR2"] and gather.gathers[2:] == ["SRR2"]
    third = _runner.invoke(app, args)
    assert third.exit_code == 0, third.output
    assert "already completed" in third.output and len(gather.gathers) == 3
    record = Config.load(Path(wd)).stages["reads_screen"]
    assert record.completed and any(key.endswith(READS_SKETCHES_TSV) for key in record.inputs)
    status = _runner.invoke(app, ["status", "-wd", wd, "--json"])
    assert status.exit_code == 0, status.output
    report = {s["name"]: s for s in json.loads(status.output)["stages"]}
    assert report["reads_screen"]["state"] == "done"
    assert (
        report["reads_screen"]["detail"] == "variants: spots_5000 (2 passed, 0 excused of 2 runs)"
    )
    # A changed threshold is a new fingerprint: the stage runs again.
    stricter = _runner.invoke(app, [*args, "--min-fraction", "0.6"])
    assert stricter.exit_code == 0, stricter.output
    assert len(gather.gathers) == 5 and _passed(Path(wd), "spots_5000") == []


def test_harness_refuses_without_sketches_and_leaves_no_record(tmp_path: Path, gtdb) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1")])
    ctx.close()
    wd = str(ctx.workdir)
    result = _runner.invoke(
        app, ["reads-screen", "-wd", wd, "--gtdb-sketch", gtdb[0], "--gtdb-lineages", gtdb[1]]
    )
    assert result.exit_code == 3, result.output
    assert "No reads sketches" in result.output
    assert "reads_screen" not in Config.load(Path(wd)).stages


# --- the duplicate check across calls ---------------------------------------------


def test_a_run_that_fails_then_passes_excuses_its_later_duplicate(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1", "SAMN1"), _row("SRR2", "SAMN1")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    ani[("SRR2", "SRR1")] = 1.0
    gather.fail = {"SRR1"}
    assert run(ctx, _params(gtdb)) == 1 and _passed(ctx.workdir, "all") == ["SRR2"]
    gather.fail = set()
    assert run(ctx, _params(gtdb)) == 1
    assert _passed(ctx.workdir, "all") == ["SRR1"]
    assert _table(ctx.workdir, "all")["SRR2"]["duplicate_of"] == "SRR1"


def test_subset_calls_give_the_result_of_one_call(tmp_path: Path, gtdb, gather, ani) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1", "SAMN1"), _row("SRR2", "SAMN1")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    ani[("SRR2", "SRR1")] = 1.0
    run(ctx, _params(gtdb, accessions=["SRR2"]))
    assert _passed(ctx.workdir, "all") == ["SRR2"]
    run(ctx, _params(gtdb, accessions=["SRR1"]))
    assert _passed(ctx.workdir, "all") == ["SRR1"], "SRR2 is now the duplicate of SRR1"
    assert gather.gathers == ["SRR2", "SRR1"], "SRR2 was not gathered again"
    # The duplicate excuse is undone when the earlier run no longer passes.
    gather.hits["SRR1"] = (_LINEAGE, 0.1)
    ctx.force = True
    run(ctx, _params(gtdb, accessions=["SRR1"]))
    table = _table(ctx.workdir, "all")
    assert table["SRR1"]["decision"] == "excused" and table["SRR2"]["decision"] == "pass"
    assert _passed(ctx.workdir, "all") == ["SRR2"]


def test_runs_dropped_from_reads_tsv_leave_the_table_and_passed(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    run(ctx, _params(gtdb))
    write_reads(ctx.workdir / READS_TSV, [_row("SRR1")])
    run(ctx, _params(gtdb))
    assert list(_table(ctx.workdir, "all")) == ["SRR1"]
    assert _passed(ctx.workdir, "all") == ["SRR1"]


def test_the_fraction_reason_names_this_commands_flag(tmp_path: Path, gtdb, gather, ani) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1")])
    _sketch(ctx.workdir, "all", ["SRR1"])
    gather.hits = {"SRR1": (_LINEAGE, 0.1)}
    run(ctx, _params(gtdb))
    reason = _table(ctx.workdir, "all")["SRR1"]["reason"]
    assert "below --min-fraction 0.3" in reason and "--screen-min-fraction" not in reason


def test_an_empty_accession_file_is_refused(tmp_path: Path) -> None:
    empty = tmp_path / "passed.txt"
    empty.write_text("# every run was excused\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="lists no accessions"):
        select_runs([_row("SRR1")], ReadsScreenParams(accession_file=str(empty)))


def test_runs_decided_with_other_settings_leave_passed_with_a_warning(
    tmp_path: Path, gtdb, gather, ani, caplog
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.35), "SRR2": (_LINEAGE, 0.5)}
    run(ctx, _params(gtdb))
    assert _passed(ctx.workdir, "all") == ["SRR1", "SRR2"]
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, _params(gtdb, accessions=["SRR2"], min_fraction=0.4))
    assert _passed(ctx.workdir, "all") == ["SRR2"]
    assert _table(ctx.workdir, "all")["SRR1"]["decision"] == "pass", "its row stays"
    assert "1 run(s) of the table were decided with other settings" in caplog.text
    assert "SRR1" in caplog.text
    assert ctx.config.stages["reads_screen"].params["stale"] == 1
    run(ctx, _params(gtdb, min_fraction=0.4))
    assert _passed(ctx.workdir, "all") == ["SRR2"]
    assert _table(ctx.workdir, "all")["SRR1"]["decision"] == "excused"


def test_a_genus_renamed_pass_keeps_its_note_through_a_redone_duplicate_check(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    renamed = _LINEAGE.replace("g__Francisella;", "g__Allofrancisella;").replace(
        "s__Francisella tularensis", "s__Allofrancisella tularensis"
    )
    ctx = _workdir(tmp_path, [_row("SRR1", "SAMN1"), _row("SRR2", "SAMN1")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (renamed, 0.5), "SRR2": (renamed, 0.5)}
    ani[("SRR2", "SRR1")] = 1.0
    run(ctx, _params(gtdb))
    table = _table(ctx.workdir, "all")
    assert table["SRR1"]["reason"] == "genus_renamed"
    assert table["SRR2"]["reason"].startswith("duplicate_isolate")
    # SRR1 drops out; SRR2's duplicate excuse is undone and its note comes back.
    gather.hits["SRR1"] = (renamed, 0.1)
    ctx.force = True
    run(ctx, _params(gtdb, accessions=["SRR1"]))
    assert _table(ctx.workdir, "all")["SRR2"]["reason"] == "genus_renamed"


def test_a_pass_with_an_unknown_submitted_genus_gets_no_note(
    tmp_path: Path, gtdb, gather, ani
) -> None:
    ctx = _workdir(
        tmp_path, [_row("SRR1", "SAMN1", genus="unknown"), _row("SRR2", "SAMN1", genus="unknown")]
    )
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    ani[("SRR2", "SRR1")] = 1.0
    run(ctx, _params(gtdb))
    assert _table(ctx.workdir, "all")["SRR2"]["reason"].startswith("duplicate_isolate")
    gather.hits["SRR1"] = (_LINEAGE, 0.1)
    ctx.force = True
    run(ctx, _params(gtdb, accessions=["SRR1"]))
    table = _table(ctx.workdir, "all")
    assert table["SRR2"]["decision"] == "pass" and table["SRR2"]["reason"] == ""


def test_a_failed_unselected_run_is_not_called_stale(
    tmp_path: Path, gtdb, gather, ani, caplog
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    _sketch(ctx.workdir, "all", ["SRR1", "SRR2"])
    gather.hits = {"SRR1": (_LINEAGE, 0.5), "SRR2": (_LINEAGE, 0.5)}
    gather.fail = {"SRR2"}
    run(ctx, _params(gtdb))
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, _params(gtdb, accessions=["SRR1"]))
    assert "decided with other settings" not in caplog.text
    assert ctx.config.stages["reads_screen"].params["stale"] == 0
    assert _passed(ctx.workdir, "all") == ["SRR1"]
