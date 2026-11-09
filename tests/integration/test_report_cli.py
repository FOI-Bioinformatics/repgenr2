"""`repgenr report` and the report refresh of the stage harness."""

from __future__ import annotations

import shutil
from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app

_runner = CliRunner()


def _genomes(tmp_path: Path) -> Path:
    src = tmp_path / "src"
    src.mkdir()
    for i, sp in enumerate(["tularensis", "tularensis", "novicida"]):
        name = f"Francisellaceae_Francisella_{sp}_GCF_00000{i}.1.fasta"
        (src / name).write_text(f">c{i}\n{'ACGT' * 50}\n")
    return src


def _ingest(tmp_path: Path) -> Path:
    wd = tmp_path / "wd"
    argv = ["ingest", "-wd", str(wd), "--genomes-dir", str(_genomes(tmp_path)), "--no-sketch"]
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, result.output
    return wd


def test_stage_writes_its_report_and_a_skip_restores_it(tmp_path: Path) -> None:
    wd = _ingest(tmp_path)
    frag = wd / "report" / "selection"
    assert (frag / "step.json").is_file()
    assert "species\tFrancisella tularensis\t2" in (frag / "genomes_per_taxon.tsv").read_text()
    shutil.rmtree(wd / "report")
    again = _runner.invoke(
        app, ["ingest", "-wd", str(wd), "--genomes-dir", str(tmp_path / "src"), "--no-sketch"]
    )
    assert again.exit_code == 0, again.output
    assert "skipping" in again.output or "already completed" in (wd / "repgenr.log").read_text()
    assert (frag / "step.json").is_file()


def test_report_command_writes_the_page(tmp_path: Path) -> None:
    wd = _ingest(tmp_path)
    result = _runner.invoke(app, ["report", "-wd", str(wd), "--figure-formats", "png,pdf"])
    assert result.exit_code == 0, result.output
    page = wd / "report" / "report.html"
    text = page.read_text()
    assert ">Genome selection</h2>" in text
    assert "Genomes were taken from a local folder of genome assemblies." in text
    assert "<td>ingest</td><td>done</td>" in text
    assert list((wd / "report" / "selection").glob("*.pdf"))


def test_report_fragments_only_and_from_fragments(tmp_path: Path) -> None:
    wd = _ingest(tmp_path)
    shutil.rmtree(wd / "report")
    only = _runner.invoke(app, ["report", "-wd", str(wd), "--fragments-only"])
    assert only.exit_code == 0, only.output
    assert "selection" in only.output
    assert not (wd / "report" / "report.html").exists()
    out = tmp_path / "nf" / "repgenr_report.html"
    stitched = _runner.invoke(
        app,
        ["report", "--from-fragments", str(wd / "report" / "selection"), "--out", str(out),
         "--title", "Pipeline report"],
    )  # fmt: skip
    assert stitched.exit_code == 0, stitched.output
    text = out.read_text()
    assert "<title>Pipeline report</title>" in text
    assert ">Genome selection</h2>" in text


def test_report_on_an_empty_directory_writes_a_page_with_a_note(tmp_path: Path) -> None:
    wd = tmp_path / "empty"
    wd.mkdir()
    result = _runner.invoke(app, ["report", "-wd", str(wd)])
    assert result.exit_code == 0, result.output
    text = (wd / "report" / "report.html").read_text()
    assert "No step of this run has results to report." in text
    assert "has no repgenr.yaml" in text


def test_report_usage_errors(tmp_path: Path) -> None:
    cases = [
        (["report"], "Give -wd"),
        (["report", "-wd", str(tmp_path), "--steps", "nope"], "Invalid report step 'nope'"),
        (["report", "-wd", str(tmp_path), "--figure-formats", "gif"], "Invalid figure format"),
        (["report", "--from-fragments", str(tmp_path)], "--from-fragments needs --out"),
        (
            ["report", "-wd", str(tmp_path), "--from-fragments", str(tmp_path), "--out", "x"],
            "not both",
        ),
    ]
    for argv, message in cases:
        result = _runner.invoke(app, argv)
        assert result.exit_code == 2, (argv, result.output)
        assert message in result.output, (argv, result.output)
    missing = _runner.invoke(app, ["report", "-wd", str(tmp_path / "absent")])
    assert missing.exit_code == 3
