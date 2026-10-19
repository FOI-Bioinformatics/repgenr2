"""reads-sketch: stream runs through fastq-dump into sourmash, index the sketches."""

from __future__ import annotations

import json
import logging
import os
import re
import stat
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli import cmd_reads
from repgenr.cli.main import app
from repgenr.core import process
from repgenr.core.config import Config
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    READS_SKETCHES_DIR,
    READS_SKETCHES_TSV,
    READS_TSV,
    ReadRow,
    read_reads_sketches,
    reads_sketch_variant,
    reads_sketch_variant_spots,
    write_reads,
)
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.sketches import READS_SKETCH_PARAMS
from repgenr.core.sratools import fastq_dump_command
from repgenr.stages.reads_sketch import (
    ReadsSketchParams,
    failure_reason,
    precheck,
    run,
    run_standalone,
    select_runs,
    sequences_sketched,
    spots_read,
)

_runner = CliRunner()

# A fake fastq-dump: four records per spot requested (two mates), the stderr
# lines of the real tool, and a failure for an accession starting with SRR0.
# Every call appends its argv to calls.jsonl beside the script.
_FASTQ_DUMP = r"""#!/bin/sh
here=$(dirname "$0")
if [ "$1" = "--version" ]; then echo "fastq-dump : 3.4.1"; exit 0; fi
printf '%s\n' "$*" >> "$here/calls.txt"
acc=""; n=100000
while [ $# -gt 0 ]; do
  case "$1" in
    -X) n=$2; shift 2;;
    -*) shift;;
    *) acc=$1; shift;;
  esac
done
case "$acc" in SRR0*)
  echo "fastq-dump.3.4.1 err: item not found while constructing within virtual database" \
    "module - the accession '$acc' could not be found." >&2
  echo "fastq-dump quit with error code 3" >&2
  exit 3;;
esac
if [ "$n" -gt 50 ]; then n=50; fi
i=0
while [ $i -lt $n ]; do
  i=$((i+1))
  printf '@%s.%s.1\nACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$acc" "$i"
  printf '@%s.%s.2\nTTGCATTGCATTGCATTGCATTGCATTGCA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n' "$acc" "$i"
done
echo "Read $n spots for $acc" >&2
echo "Written $n spots for $acc" >&2
"""

# A fake sourmash: counts the records on stdin, writes the -o file, and
# prints the real tool's closing line. A --name starting with ERR0 makes
# it fail after reading a few lines, while fastq-dump is still writing.
_SOURMASH = r"""#!/bin/sh
if [ "$1" = "--version" ]; then echo "sourmash 4.9.4"; exit 0; fi
out=""; name=""
while [ $# -gt 0 ]; do
  case "$1" in
    -o) out=$2; shift 2;;
    --name) name=$2; shift 2;;
    *) shift;;
  esac
done
case "$name" in ERR0*) head -n 8 > /dev/null; echo "error: cannot sketch $name" >&2; exit 2;; esac
lines=$(wc -l < /dev/stdin | tr -d ' ')
n=$((lines / 4))
if [ "$n" -eq 0 ]; then echo "error: no sequences" >&2; exit 1; fi
printf 'fake sketch of %s sequences\n' "$n" > "$out"
echo "calculated 1 signature for $n sequences taken from 1 files" >&2
"""


@pytest.fixture
def fake_tools(tmp_path: Path, monkeypatch) -> Path:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, text in (("fastq-dump", _FASTQ_DUMP), ("sourmash", _SOURMASH)):
        script = bin_dir / name
        script.write_text(text, encoding="utf-8")
        script.chmod(script.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    return bin_dir


def _row(acc: str, platform: str = "ILLUMINA") -> ReadRow:
    return ReadRow(
        run_accession=acc,
        biosample=f"SAMN{acc}",
        bioproject="PRJNA1",
        organism="Francisella tularensis",
        taxid="263",
        platform=platform,
        instrument_model="x",
        layout="PAIRED",
        bases=1000,
        read_count=10,
        family="Francisellaceae",
        genus="Francisella",
        species="tularensis",
    )


def _workdir(tmp_path: Path, rows: list[ReadRow]) -> WorkdirContext:
    ctx = WorkdirContext(tmp_path / "wd", create=True, logger=logging.getLogger("test"))
    write_reads(ctx.workdir / READS_TSV, rows)
    ctx.config.record_stage("reads", params={}, completed="2026-10-19T00:00:00+00:00")
    ctx.save_config()
    return ctx


def _calls(bin_dir: Path) -> list[str]:
    path = bin_dir / "calls.txt"
    return path.read_text(encoding="utf-8").splitlines() if path.exists() else []


# --- selection ---------------------------------------------------------------------


def test_select_runs_takes_every_run_without_a_selection() -> None:
    rows = [_row("SRR1"), _row("SRR2", "OXFORD_NANOPORE")]
    assert [r.run_accession for r in select_runs(rows, ReadsSketchParams())] == ["SRR1", "SRR2"]


def test_select_runs_by_accession_file_and_platform(tmp_path: Path) -> None:
    rows = [_row("SRR1"), _row("SRR2", "OXFORD_NANOPORE"), _row("SRR3")]
    listing = tmp_path / "picked.txt"
    listing.write_text("SRR2  # the long-read run\nSRR3\n", encoding="utf-8")
    params = ReadsSketchParams(accessions=["SRR1"], accession_file=str(listing))
    assert [r.run_accession for r in select_runs(rows, params)] == ["SRR1", "SRR2", "SRR3"]
    params = ReadsSketchParams(accession_file=str(listing), platform="ont")
    assert [r.run_accession for r in select_runs(rows, params)] == ["SRR2"]


def test_select_runs_refuses_an_accession_outside_reads_tsv() -> None:
    with pytest.raises(UserInputError, match="Not in reads.tsv: SRR9"):
        select_runs([_row("SRR1")], ReadsSketchParams(accessions=["SRR9"]))


def test_precheck_needs_a_completed_reads_stage(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True, logger=logging.getLogger("test"))
    with pytest.raises(WorkdirError, match="completed reads stage"):
        precheck(ctx, ReadsSketchParams())


# --- the pipeline --------------------------------------------------------------------


def test_fastq_dump_command_streams_biological_reads_with_the_cap() -> None:
    assert fastq_dump_command("SRR1", None) == [
        "fastq-dump",
        "-Z",
        "--split-spot",
        "--skip-technical",
        "SRR1",
    ]
    assert fastq_dump_command("SRR1", 100)[-3:] == ["-X", "100", "SRR1"]


def test_parsers_read_the_tools_closing_lines() -> None:
    assert spots_read(["Rejected 3 READS because of filtering", "Read 20000 spots for X"]) == 20000
    assert spots_read(["nothing"]) is None
    assert sequences_sketched(
        ["calculated 1 signature for 40000 sequences taken from 1 files"]
    ) == (40000)
    assert sequences_sketched(["calculated 3 signatures for 7 sequences"]) == 7
    tail = ["err: item not found while constructing", "fastq-dump quit with error code 3"]
    assert failure_reason((3, 0), tail, []) == (
        "fastq-dump exited 3: err: item not found while constructing"
    )
    assert failure_reason((3, 0), ["Read 0 spots", "quit with error code 3"], []) == (
        "fastq-dump exited 3: Read 0 spots"
    )
    assert failure_reason((141, 1), ["x"], ["error: bad"]) == "sourmash exited 1: error: bad"
    assert failure_reason((-9, 1), ["x"], ["error: bad"]) == "sourmash exited 1: error: bad"
    assert failure_reason((3, 1), ["err: gone"], ["error: no sequences"]) == (
        "fastq-dump exited 3: err: gone"
    )
    # fastq-dump exits 1 with an empty stderr when its reader goes away.
    assert failure_reason((1, 1), [], ["error: bad"]) == "sourmash exited 1: error: bad"
    assert failure_reason((1, 0), [], []) == "fastq-dump exited 1: no output"
    assert failure_reason((0, 0), [], []) == "no sketch was written"


def test_run_pipeline_joins_two_processes_and_reports_both_exits() -> None:
    logger = logging.getLogger("test")
    result = process.run_pipeline(
        ["sh", "-c", "printf 'a\\nb\\n'; echo made >&2"],
        ["sh", "-c", "wc -l | tr -d ' ' >&2"],
        logger=logger,
    )
    assert result.ok and result.returncodes == (0, 0)
    assert result.tails[0] == ["made"] and result.tails[1] == ["2"]
    # A consumer that stops early leaves the producer with a broken pipe (or
    # killed); the producer's exit is non-zero and the consumer's is reported.
    result = process.run_pipeline(
        ["sh", "-c", "yes"], ["sh", "-c", "head -c 10 > /dev/null; exit 4"], logger=logger
    )
    assert result.returncodes[1] == 4 and result.returncodes[0] != 0
    assert not result.timed_out
    # Each half has its own environment.
    result = process.run_pipeline(
        ["sh", "-c", "echo $WHO >&2"],
        ["sh", "-c", "echo $WHO >&2"],
        logger=logger,
        envs=({"WHO": "first"}, {"WHO": "second"}),
    )
    assert result.tails == (["first"], ["second"])


# --- the workdir form ---------------------------------------------------------------


def test_run_sketches_every_run_and_indexes_failures(tmp_path: Path, fake_tools: Path) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR02"), _row("SRR3")])
    sketched = run(ctx, ReadsSketchParams(max_spots=20, jobs=2))
    assert sketched == 2
    out = ctx.workdir / READS_SKETCHES_DIR / "spots_20"
    rows = {r.run_accession: r for r in read_reads_sketches(out / READS_SKETCHES_TSV)}
    assert [r for r in rows] == ["SRR1", "SRR02", "SRR3"]  # the order of reads.tsv
    ok = rows["SRR1"]
    assert ok.status == "ok" and ok.max_spots == 20 and ok.spots_read == 20
    assert ok.sequences == 40 and ok.sketch_params == READS_SKETCH_PARAMS
    assert ok.sketch_file == "SRR1.sig.zip"
    assert (out / "SRR1.sig.zip").read_text().startswith("fake sketch of 40")
    assert ok.seconds is not None and ok.seconds >= 0
    bad = rows["SRR02"]
    assert bad.status == "failed" and bad.reason.startswith("fastq-dump exited 3:")
    assert "could not be found" in bad.reason and not (out / "SRR02.sig.zip").exists()
    assert not [p for p in out.iterdir() if p.name.startswith(".")], "no partial sketch left"
    record = ctx.config.stages["reads_sketch"]
    assert record.completed and record.tool == "sourmash"
    assert record.params["max_spots"] == 20 and record.params["selected"] == 3
    assert record.params["variant"] == "spots_20" and record.params["variants"] == {"spots_20": 2}
    assert (record.params["sketched"], record.params["failed"], record.params["present"]) == (
        2,
        1,
        2,
    )
    assert record.tool_versions["fastq-dump"].startswith("3.4") and "sourmash" in (
        record.tool_versions
    )
    calls = _calls(fake_tools)
    assert len(calls) == 3 and all("-X 20" in c for c in calls)


def test_a_repeat_skips_current_sketches_and_a_new_cap_is_its_own_variant(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    root = ctx.workdir / READS_SKETCHES_DIR
    run(ctx, ReadsSketchParams(max_spots=10))
    assert len(_calls(fake_tools)) == 2
    assert run(ctx, ReadsSketchParams(max_spots=10)) == 0
    assert len(_calls(fake_tools)) == 2, "nothing streamed again"
    assert run(ctx, ReadsSketchParams(max_spots=30)) == 2
    assert run(ctx, ReadsSketchParams()) == 2
    assert len(_calls(fake_tools)) == 6
    # Three variants side by side, each with its own index and sketches.
    assert sorted(p.name for p in root.iterdir()) == ["all", "spots_10", "spots_30"]
    for name, cap in (("all", None), ("spots_10", 10), ("spots_30", 30)):
        rows = read_reads_sketches(root / name / READS_SKETCHES_TSV)
        assert {r.max_spots for r in rows} == {cap} and (root / name / "SRR1.sig.zip").is_file()
    assert reads_sketch_variant(None) == "all" and reads_sketch_variant(30) == "spots_30"
    assert reads_sketch_variant_spots("all") == (True, None)
    assert reads_sketch_variant_spots("spots_30") == (True, 30)
    for odd in ("spots_0", "spots_007", "spots_\u00b2", "spots_", "other", ".home.x"):
        assert reads_sketch_variant_spots(odd) == (False, None), odd
    # A stray directory and an old flat-layout index are ignored, with a warning.
    (root / "spots_\u00b2").mkdir()
    (root / "reads_sketches.tsv").write_text("x\n", encoding="utf-8")
    assert ctx.config.stages["reads_sketch"].params["variants"] == {
        "all": 2,
        "spots_10": 2,
        "spots_30": 2,
    }
    # --force (ctx.force) streams the current variant again, the others untouched.
    ctx.force = True
    assert run(ctx, ReadsSketchParams(max_spots=30)) == 2
    assert len(_calls(fake_tools)) == 8
    assert ctx.config.stages["reads_sketch"].params["variants"] == {
        "all": 2,
        "spots_10": 2,
        "spots_30": 2,
    }


def test_a_subset_adds_to_the_set_and_the_rest_stay_pending(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2"), _row("SRR3")])
    run(ctx, ReadsSketchParams(accessions=["SRR2"]))
    rows = {
        r.run_accession: r.status
        for r in read_reads_sketches(ctx.workdir / READS_SKETCHES_DIR / "all" / READS_SKETCHES_TSV)
    }
    assert rows == {"SRR1": "pending", "SRR2": "ok", "SRR3": "pending"}
    run(ctx, ReadsSketchParams(accessions=["SRR3"]))
    rows = {
        r.run_accession: r.status
        for r in read_reads_sketches(ctx.workdir / READS_SKETCHES_DIR / "all" / READS_SKETCHES_TSV)
    }
    assert rows == {"SRR1": "pending", "SRR2": "ok", "SRR3": "ok"}
    assert (ctx.workdir / READS_SKETCHES_DIR / "all" / "SRR2.sig.zip").exists()
    assert len(_calls(fake_tools)) == 2
    # No cap: the whole run, recorded as an empty max_spots.
    assert all("-X" not in c for c in _calls(fake_tools))
    assert ctx.config.stages["reads_sketch"].params["present"] == 2


def test_every_run_failing_is_an_error_and_records_nothing(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR01"), _row("SRR02")])
    with pytest.raises(WorkdirError, match="Every one of the 2 selected run"):
        run(ctx, ReadsSketchParams())
    rows = read_reads_sketches(ctx.workdir / READS_SKETCHES_DIR / "all" / READS_SKETCHES_TSV)
    assert [r.status for r in rows] == ["failed", "failed"]
    assert "reads_sketch" not in ctx.config.stages, "no completed record of a failed stage"


def test_a_sourmash_failure_mid_stream_is_blamed_on_sourmash(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("ERR01"), _row("SRR2")])
    assert run(ctx, ReadsSketchParams()) == 1
    out = ctx.workdir / READS_SKETCHES_DIR / "all"
    rows = {r.run_accession: r for r in read_reads_sketches(out / READS_SKETCHES_TSV)}
    assert rows["ERR01"].status == "failed"
    assert rows["ERR01"].reason == "sourmash exited 2: error: cannot sketch ERR01"
    assert not (out / "ERR01.sig.zip").exists()
    assert not [p for p in out.iterdir() if p.name.startswith(".")], "no partial sketch left"


def test_run_refuses_a_selection_outside_reads_tsv_before_any_tool(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1")])
    with pytest.raises(UserInputError, match="SRR7"):
        run(ctx, ReadsSketchParams(accessions=["SRR7"]))
    assert _calls(fake_tools) == [] and "reads_sketch" not in ctx.config.stages


# --- the standalone form -------------------------------------------------------------


def test_standalone_sketches_bare_accessions_into_a_directory(
    tmp_path: Path, fake_tools: Path
) -> None:
    out = tmp_path / "sk"
    logger = logging.getLogger("test")
    n = run_standalone(
        ["SRR1", "SRR1", "SRR02"], out, max_spots=5, jobs=1, force=False, logger=logger
    )
    assert n == 1
    rows = {r.run_accession: r for r in read_reads_sketches(out / "spots_5" / READS_SKETCHES_TSV)}
    assert set(rows) == {"SRR1", "SRR02"} and rows["SRR1"].status == "ok"
    assert (out / "spots_5" / "SRR1.sig.zip").exists()
    with pytest.raises(UserInputError, match="name the runs"):
        run_standalone([], out, max_spots=None, jobs=None, force=False, logger=logger)
    with pytest.raises(UserInputError, match="Not a run accession.*\\.\\./x"):
        run_standalone(["../x"], out, max_spots=None, jobs=None, force=False, logger=logger)


def test_standalone_refuses_platform(tmp_path: Path) -> None:
    result = _runner.invoke(
        app, ["reads-sketch", "-o", str(tmp_path / "o"), "--accession", "SRR1", "--platform", "ont"]
    )
    assert result.exit_code != 0 and "needs -wd" in _plain(result.output)


# --- CLI wiring and the harness ------------------------------------------------------


def test_reads_sketch_command_builds_its_params(monkeypatch, tmp_path: Path) -> None:
    seen = {}

    def fake_run(stage, workdir, build, *, create=False):
        seen["stage"], seen["params"] = stage, build()

    monkeypatch.setattr(cmd_reads, "_run", fake_run)
    listing = tmp_path / "acc.txt"
    listing.write_text("SRR1\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            *("reads-sketch", "-wd", str(tmp_path / "wd"), "--accession", "SRR2"),
            *("--accession-file", str(listing), "--platform", "ont"),
            *("--max-spots", "100000", "--jobs", "3"),
        ],
    )
    assert result.exit_code == 0, result.output
    p = seen["params"]
    assert seen["stage"] == "reads_sketch"
    assert p.accessions == ["SRR2"] and p.accession_file == str(listing)
    assert (p.platform, p.max_spots, p.jobs) == ("ont", 100000, 3)


def _plain(output: str) -> str:
    """Rich's error panel without escapes, box characters and line wrapping."""
    text = re.sub(r"\x1b\[[0-9;]*m", "", output)
    text = re.sub(r"[\u2500-\u257f]", " ", text)
    return " ".join(text.split())


def test_reads_sketch_command_needs_a_workdir_or_an_output_directory(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["reads-sketch", "--accession", "SRR1"])
    assert result.exit_code != 0 and "-o DIR" in _plain(result.output)
    result = _runner.invoke(
        app, ["reads-sketch", "-wd", str(tmp_path / "wd"), "-o", str(tmp_path / "o")]
    )
    assert result.exit_code != 0 and "standalone" in _plain(result.output)


def test_harness_skips_a_repeat_and_reruns_on_a_new_selection(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR2")])
    ctx.close()
    wd = str(ctx.workdir)
    first = _runner.invoke(app, ["reads-sketch", "-wd", wd, "--accession", "SRR1"])
    assert first.exit_code == 0, first.output
    assert len(_calls(fake_tools)) == 1
    again = _runner.invoke(app, ["reads-sketch", "-wd", wd, "--accession", "SRR1"])
    assert again.exit_code == 0, again.output
    assert "already completed" in again.output and len(_calls(fake_tools)) == 1
    more = _runner.invoke(app, ["reads-sketch", "-wd", wd])
    assert more.exit_code == 0, more.output
    assert len(_calls(fake_tools)) == 2, "only the run not yet sketched"
    cfg = Config.load(Path(wd))
    record = cfg.stages["reads_sketch"]
    assert record.completed and record.fingerprint
    assert any(key.endswith(READS_SKETCHES_TSV) for key in record.inputs)
    status = _runner.invoke(app, ["status", "-wd", wd, "--json"])
    assert status.exit_code == 0, status.output
    report = {s["name"]: s for s in json.loads(status.output)["stages"]}
    assert report["reads_sketch"]["state"] == "done"
    assert report["reads_sketch"]["detail"] == "variants: all (2 of 2 runs sketched)"


def test_harness_streams_a_failed_run_again_on_a_repeat(tmp_path: Path, fake_tools: Path) -> None:
    """A selected run without a sketch is a missing deliverable: the repeat
    runs the stage, which skips the current sketch and streams the rest."""
    ctx = _workdir(tmp_path, [_row("SRR1"), _row("SRR02")])
    ctx.close()
    wd = str(ctx.workdir)
    first = _runner.invoke(app, ["reads-sketch", "-wd", wd])
    assert first.exit_code == 0, first.output
    assert "1 run(s) failed; a repeat" in first.output
    assert len(_calls(fake_tools)) == 2
    again = _runner.invoke(app, ["reads-sketch", "-wd", wd])
    assert again.exit_code == 0, again.output
    assert "already completed" not in again.output
    calls = _calls(fake_tools)
    assert len(calls) == 3 and calls[-1].endswith("SRR02"), "only the failed run again"
    # A sketch removed by hand is streamed again the same way.
    (Path(wd) / READS_SKETCHES_DIR / "all" / "SRR1.sig.zip").unlink()
    third = _runner.invoke(app, ["reads-sketch", "-wd", wd, "--accession", "SRR1"])
    assert third.exit_code == 0, third.output
    assert len(_calls(fake_tools)) == 4 and _calls(fake_tools)[-1].endswith("SRR1")


def test_harness_leaves_no_completed_record_when_every_run_fails(
    tmp_path: Path, fake_tools: Path
) -> None:
    ctx = _workdir(tmp_path, [_row("SRR01")])
    ctx.close()
    wd = str(ctx.workdir)
    result = _runner.invoke(app, ["reads-sketch", "-wd", wd])
    assert result.exit_code == 3, result.output
    record = Config.load(Path(wd)).stages.get("reads_sketch")
    assert record is None or not record.completed
    status = _runner.invoke(app, ["status", "-wd", wd, "--json"])
    names = {s["name"]: s["state"] for s in json.loads(status.output)["stages"]}
    assert names.get("reads_sketch", "pending") != "done"


def test_standalone_sketches_the_runs_of_a_reads_tsv(tmp_path: Path, fake_tools: Path) -> None:
    reads = tmp_path / "reads.tsv"
    write_reads(reads, [_row("SRR1"), _row("SRR2")])
    out = tmp_path / "sketches"
    versions = tmp_path / "tool_versions.yml"
    result = _runner.invoke(
        app,
        [
            "reads-sketch",
            "-o",
            str(out),
            "--reads-tsv",
            str(reads),
            "--max-spots",
            "7",
            "--versions-out",
            str(versions),
        ],
    )
    assert result.exit_code == 0, result.output
    index = {r.run_accession: r for r in read_reads_sketches(out / "spots_7" / READS_SKETCHES_TSV)}
    assert sorted(index) == ["SRR1", "SRR2"] and index["SRR1"].status == "ok"
    assert "sourmash: 4.9.4" in versions.read_text(encoding="utf-8")
    subset = _runner.invoke(
        app, ["reads-sketch", "-o", str(out), "--reads-tsv", str(reads), "--accession", "SRR9"]
    )
    assert subset.exit_code == 2 and "Not in" in subset.output
    mixed = _runner.invoke(
        app, ["reads-sketch", "-wd", str(tmp_path / "wd"), "--reads-tsv", str(reads)]
    )
    assert mixed.exit_code != 0 and "standalone form" in mixed.output.replace("\n", " ")
