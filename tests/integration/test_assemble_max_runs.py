"""assemble --max-runs: a bound on the runs with a finished assembly, the rest
deferred in excused_runs.tsv, with the fake assembler (assemble_fakes.py)."""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.base import _stage_fingerprint
from repgenr.cli.main import app
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    EXCUSED_RUNS_TSV,
    READS_TSV,
    ReadRow,
    read_excused_runs,
    read_selection,
    write_reads,
)
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.integrity import check_genome_completeness
from repgenr.stages import dereplicate
from repgenr.stages.assemble import AssembleParams, run

sys.path.insert(0, str(Path(__file__).resolve().parent))
from assemble_fakes import FakeAssembler as _FakeAssembler  # noqa: E402
from assemble_fakes import read_row as _row  # noqa: E402
from assemble_fakes import register_fake_assembler, unregister_fake_assembler  # noqa: E402

_LOG = logging.getLogger("test")
_RUNS = ["SRR1", "SRR2", "SRR3", "SRR4"]


@pytest.fixture
def fake_assembler():
    register_fake_assembler()
    yield
    unregister_fake_assembler()


def _prepare(workdir: Path, tmp_path: Path, runs: list[str] = _RUNS) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    rows: list[ReadRow] = [_row(tmp_path, r) for r in runs]
    write_reads(workdir / READS_TSV, rows)
    return ctx


def _excused(workdir: Path) -> dict[str, tuple[str, str]]:
    path = workdir / EXCUSED_RUNS_TSV
    if not path.exists():
        return {}
    return {e.run_accession: (e.step, e.reason) for e in read_excused_runs(path)}


def _selected(workdir: Path) -> list[str]:
    return [r.accession for r in read_selection(workdir / "selection.tsv")]


def test_runs_beyond_n_are_deferred_in_reads_order(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    n = run(ctx, AssembleParams(assembler="fakeasm", max_runs=2))
    assert n == 2
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2"]
    assert _selected(workdir) == ["SRR1", "SRR2"]
    assert {g.accession for g in ctx.manifest.all_genomes()} == {"SRR1", "SRR2"}
    assert _excused(workdir) == {"SRR3": ("assemble", "deferred"), "SRR4": ("assemble", "deferred")}
    assert not (workdir / "assemblies" / "SRR3").exists()  # nothing fetched for a deferred run
    params = ctx.config.stages["assemble"].params
    assert params["max_runs"] == 2
    assert params["n_assembled"] == 2 and params["n_deferred"] == 2 and params["n_excused"] == 0


def test_n_larger_than_the_runs_defers_nothing(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    assert run(ctx, AssembleParams(assembler="fakeasm", max_runs=10)) == 4
    assert _excused(workdir) == {}
    assert ctx.config.stages["assemble"].params["n_deferred"] == 0


def test_a_larger_n_assembles_deferred_runs_and_keeps_the_rest_deferred(
    workdir, tmp_path, fake_assembler
) -> None:
    ctx = _prepare(workdir, tmp_path)
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=1))
    assert _FakeAssembler.calls == ["SRR1"]
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=3))
    # Finished runs count toward N and are reused; two places remain.
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2", "SRR3"]
    assert _selected(workdir) == ["SRR1", "SRR2", "SRR3"]
    assert _excused(workdir) == {"SRR4": ("assemble", "deferred")}
    run(ctx, AssembleParams(assembler="fakeasm"))
    assert sorted(_FakeAssembler.calls) == _RUNS
    assert _excused(workdir) == {}
    assert not (workdir / EXCUSED_RUNS_TSV).exists()


def test_the_same_n_again_assembles_nothing_new(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=2))
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=2))
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2"]
    assert _selected(workdir) == ["SRR1", "SRR2"]


def test_finished_runs_beyond_n_are_kept_not_deferred(
    workdir, tmp_path, fake_assembler, caplog
) -> None:
    ctx = _prepare(workdir, tmp_path)
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=3))
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, AssembleParams(assembler="fakeasm", max_runs=1))
    assert "3 finished run(s) are kept beyond N" in caplog.text
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2", "SRR3"]
    assert _selected(workdir) == ["SRR1", "SRR2", "SRR3"]
    assert _excused(workdir) == {"SRR4": ("assemble", "deferred")}


def test_a_deferred_run_with_a_marker_is_reused(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    run(ctx, AssembleParams(assembler="fakeasm"))
    # reads.tsv reordered: SRR4 is now last but already finished.
    write_reads(workdir / READS_TSV, [_row(tmp_path, r) for r in ["SRR5", "SRR1", "SRR4"]])
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=2))
    assert _FakeAssembler.calls.count("SRR4") == 1
    assert "SRR5" not in _FakeAssembler.calls
    assert _selected(workdir) == ["SRR1", "SRR4"]
    assert _excused(workdir) == {"SRR5": ("assemble", "deferred")}


def test_runs_excused_up_front_take_no_place(workdir, tmp_path, fake_assembler) -> None:
    ctx = WorkdirContext(workdir, create=True)
    write_reads(
        workdir / READS_TSV,
        [
            _row(tmp_path, "NOMIRROR", fastq_urls=(), fastq_md5=(), fastq_bytes=()),
            _row(tmp_path, "SRR1"),
            _row(tmp_path, "SRR2"),
        ],
    )
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=1))
    assert _FakeAssembler.calls == ["SRR1"]
    assert _excused(workdir) == {
        "NOMIRROR": ("fetch", "no_fastq_mirror"),
        "SRR2": ("assemble", "deferred"),
    }
    assert ctx.config.stages["assemble"].params["n_excused"] == 1


def test_deferred_runs_never_count_as_judged(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path, ["SRR1", "SRR2"])
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=1))
    assert _selected(workdir) == ["SRR1"]
    # SRR1 gone from its marker and failing download: nothing judged, SRR2 deferred.
    (workdir / "assemblies" / "SRR1" / "assembly.ok").unlink()
    write_reads(
        workdir / READS_TSV,
        [_row(tmp_path, "SRR1", fastq_md5=("0" * 32, "0" * 32)), _row(tmp_path, "SRR2")],
    )
    with pytest.raises(WorkdirError, match="from an earlier assemble call was kept"):
        run(ctx, AssembleParams(assembler="fakeasm", max_runs=1))
    assert _selected(workdir) == ["SRR1"]
    assert _excused(workdir)["SRR2"] == ("assemble", "deferred")


def test_append_with_max_runs(workdir, tmp_path, fake_assembler) -> None:
    from repgenr.core.contracts import SelectionRow, write_selection
    from repgenr.core.manifest import record_from_selection

    ctx = _prepare(workdir, tmp_path, ["SRR1", "SRR2", "SRR3"])
    existing = SelectionRow("GCF_1", "Fam", "Gen", "sp", False, "Fam_Gen_sp_GCF_1.fasta")
    ctx.genomes_dir.mkdir(parents=True)
    (ctx.genomes_dir / existing.filename).write_text(">c\nACGT\n", encoding="utf-8")
    write_selection(workdir / "selection.tsv", [existing])
    ctx.manifest.replace_genomes([record_from_selection(existing, "gtdb")])
    run(ctx, AssembleParams(assembler="fakeasm", append=True, max_runs=1))
    assert _selected(workdir) == ["GCF_1", "SRR1"]
    run(ctx, AssembleParams(assembler="fakeasm", append=True, max_runs=2))
    assert _selected(workdir) == ["GCF_1", "SRR1", "SRR2"]
    assert _excused(workdir) == {"SRR3": ("assemble", "deferred")}
    assert {g.accession for g in ctx.manifest.all_genomes()} == {"GCF_1", "SRR1", "SRR2"}
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []


def test_dereplicate_guard_accepts_deferred_runs(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    run(ctx, AssembleParams(assembler="fakeasm", max_runs=2))
    assert check_genome_completeness(ctx.genomes_dir, workdir, logger=_LOG) == []
    dereplicate.precheck(ctx, dereplicate.DereplicateParams())


def test_max_runs_below_one_is_refused(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    with pytest.raises(UserInputError, match="--max-runs"):
        run(ctx, AssembleParams(assembler="fakeasm", max_runs=0))


def test_max_runs_is_part_of_the_resume_fingerprint() -> None:
    one = _stage_fingerprint("assemble", AssembleParams(max_runs=1), {}, {})
    two = _stage_fingerprint("assemble", AssembleParams(max_runs=2), {}, {})
    none = _stage_fingerprint("assemble", AssembleParams(), {}, {})
    assert len({one, two, none}) == 3


def test_cli_status_detail_and_a_larger_n_rerunning_the_stage(
    workdir, tmp_path, fake_assembler
) -> None:
    _prepare(workdir, tmp_path)
    runner = CliRunner()
    base = ["assemble", "-wd", str(workdir), "--assembler", "fakeasm", "--no-reads-sketch"]
    first = runner.invoke(app, [*base, "--max-runs", "1"])
    assert first.exit_code == 0, first.output
    status = runner.invoke(app, ["status", "-wd", str(workdir)])
    line = next(ln for ln in status.output.splitlines() if "assemble" in ln)
    assert "1 assembled, 3 deferred" in line, status.output
    report = json.loads(runner.invoke(app, ["status", "-wd", str(workdir), "--json"]).output)
    entry = next(s for s in report["stages"] if s["name"] == "assemble")
    assert entry["detail"].startswith("1 assembled, 3 deferred")
    # The same N is a finished stage; a larger N is a new request.
    again = runner.invoke(app, [*base, "--max-runs", "1"])
    assert again.exit_code == 0 and _FakeAssembler.calls == ["SRR1"]
    more = runner.invoke(app, [*base, "--max-runs", "3"])
    assert more.exit_code == 0, more.output
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2", "SRR3"]
    status = runner.invoke(app, ["status", "-wd", str(workdir)])
    line = next(ln for ln in status.output.splitlines() if "assemble" in ln)
    assert "3 assembled, 1 deferred" in line, status.output


def test_cli_refuses_max_runs_zero(workdir, tmp_path) -> None:
    _prepare(workdir, tmp_path)
    result = CliRunner().invoke(app, ["assemble", "-wd", str(workdir), "--max-runs", "0"])
    assert result.exit_code != 0


# --- --accession / --accession-file ------------------------------------------------


def test_a_selection_assembles_the_named_runs_and_defers_the_rest(
    workdir, tmp_path, fake_assembler
) -> None:
    ctx = _prepare(workdir, tmp_path)
    listing = tmp_path / "passed.txt"
    listing.write_text("SRR3  # from reads-screen\n", encoding="utf-8")
    params = AssembleParams(assembler="fakeasm", accessions=["SRR1"], accession_file=str(listing))
    n = run(ctx, params)
    assert n == 2
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR3"]
    assert _selected(workdir) == ["SRR1", "SRR3"]
    assert _excused(workdir) == {"SRR2": ("assemble", "deferred"), "SRR4": ("assemble", "deferred")}
    assert ctx.config.stages["assemble"].params["n_deferred"] == 2
    # Finished runs are kept outside a later selection; the selection adds to them.
    run(ctx, AssembleParams(assembler="fakeasm", accessions=["SRR2"]))
    assert sorted(_FakeAssembler.calls) == ["SRR1", "SRR2", "SRR3"]
    assert _selected(workdir) == ["SRR1", "SRR2", "SRR3"]
    assert _excused(workdir) == {"SRR4": ("assemble", "deferred")}


def test_a_selection_outside_reads_tsv_is_refused(workdir, tmp_path, fake_assembler) -> None:
    ctx = _prepare(workdir, tmp_path)
    with pytest.raises(UserInputError, match="Not in reads.tsv: SRR9"):
        run(ctx, AssembleParams(assembler="fakeasm", accessions=["SRR9"]))
    assert not (workdir / "assemblies").exists()


def test_the_selection_is_part_of_the_resume_fingerprint() -> None:
    base = AssembleParams(assembler="fakeasm")
    named = AssembleParams(assembler="fakeasm", accessions=["SRR1"])
    assert _stage_fingerprint("assemble", base, {}, {}) != _stage_fingerprint(
        "assemble", named, {}, {}
    )
