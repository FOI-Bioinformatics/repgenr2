"""The report layer: stats files, step fragments, the HTML page and the hook."""

from __future__ import annotations

import base64
import logging
import re
import xml.etree.ElementTree as ET
from datetime import UTC, datetime
from pathlib import Path

import pytest

from repgenr import reporting
from repgenr.core.config import Config
from repgenr.core.contracts import SELECTION_TSV, SelectionRow, write_selection
from repgenr.core.selection_summary import SELECTION_STATS, write_selection_summary
from repgenr.core.stage_stats import read_stats, write_stats
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.hooks import refresh_for_stage
from repgenr.reporting.model import KeyNumber, MethodsInfo, StepReport, Table
from repgenr.reporting.registry import BY_NAME, STAGE_TO_STEPS

LOG = logging.getLogger("test")


def _row(i: int, species: str, *, outgroup: bool = False, comp: float | None = 98.0):
    acc = f"GCF_{i:09d}.1"
    return SelectionRow(
        acc,
        "Francisellaceae",
        "Francisella",
        species,
        outgroup,
        f"Francisellaceae_Francisella_{species}_{acc}.fasta",
        comp,
        None if comp is None else 0.5,
        i % 5 == 0,
    )


def _workdir(tmp_path: Path, rows=None, *, stats=None, record: bool = True) -> Path:
    wd = tmp_path / "wd"
    wd.mkdir()
    rows = rows or [_row(i, f"sp{i % 3}") for i in range(12)] + [_row(99, "out", outgroup=True)]
    write_selection(wd / SELECTION_TSV, rows)
    write_selection_summary(wd, rows, LOG, stats=stats)
    if record:
        cfg = Config.load(wd)
        cfg.record_stage(
            "metadata",
            tool="gtdb-table",
            params={
                "source": "tsv",
                "release": "226",
                "target_genus": "Francisella",
                "level": "genus",
                "outgroup": "GCF_000000099.1",
            },
            completed=datetime.now(UTC).isoformat(),
        )
        cfg.save(wd)
    return wd


def test_stage_stats_round_trip_and_kind(tmp_path: Path) -> None:
    path = write_stats(tmp_path / "s.json", "selection", {"candidates": 5})
    assert read_stats(path, "selection") == {"candidates": 5}
    assert read_stats(path, "reads") is None
    assert read_stats(tmp_path / "absent.json", "selection") is None
    path.write_text("{not json")
    assert read_stats(path, "selection") is None


def test_selection_summary_rewrites_its_stats(tmp_path: Path) -> None:
    rows = [_row(1, "a")]
    write_selection_summary(tmp_path, rows, LOG, stats={"candidates": 10, "limit": 1})
    assert read_stats(tmp_path / SELECTION_STATS, "selection") == {"candidates": 10, "limit": 1}
    # Another entry stage without counts must not leave the earlier ones.
    write_selection_summary(tmp_path, rows, LOG)
    assert read_stats(tmp_path / SELECTION_STATS, "selection") == {}


def test_step_report_round_trip() -> None:
    rep = StepReport(
        step="x",
        title="X",
        key_numbers=[KeyNumber("n", 3)],
        tables=[Table("t", "cap", "t.tsv", ("a", "b"), 2, numeric=("b",))],
        methods=[MethodsInfo("metadata", tool="gtdb-table", params={"k": 1})],
    )
    assert StepReport.from_dict(rep.to_dict()) == rep
    with pytest.raises(ValueError):
        StepReport.from_dict({**rep.to_dict(), "schema": "other/9"})


def test_selection_fragment_tables_and_figures(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, stats={"candidates": 40, "limit": 12})
    assert refresh_step(wd, BY_NAME["selection"], formats=("png", "svg")) == "written"
    frag = wd / "report" / "selection"
    rep = read_fragment(frag)
    assert rep is not None and rep.step == "selection"
    assert (frag / "selection_counts.tsv").read_text().splitlines() == [
        "quantity\tgenomes",
        "Candidate genomes\t40",
        "Selected ingroup genomes\t12",
        "Outgroup genomes\t1",
        "GTDB species representatives among the selected\t3",
    ]
    taxa = (frag / "genomes_per_taxon.tsv").read_text().splitlines()
    assert taxa[0] == "level\tname\tgenomes"
    assert "species\tFrancisella sp0\t4" in taxa
    assert {f.id for f in rep.figures} == {"genomes_per_species", "assembly_quality"}
    for fig in rep.figures:
        assert (frag / fig.files["png"]).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        root = ET.parse(frag / fig.files["svg"]).getroot()
        assert root.tag.endswith("svg")
        # No date in the figure file, so the same data give the same file.
        assert "<dc:date>" not in (frag / fig.files["svg"]).read_text()
    assert [m.stage for m in rep.methods] == ["metadata"]
    assert "chosen from 40 candidate genomes by the limit of 12" in rep.summary[0]
    assert "1 outgroup genome." in rep.summary[0]


def test_fragment_is_current_then_stale_then_removed(tmp_path: Path) -> None:
    wd = _workdir(tmp_path)
    spec = BY_NAME["selection"]
    assert refresh_step(wd, spec, formats=("png",)) == "written"
    assert refresh_step(wd, spec, formats=("png",)) == "current"
    # Another figure format rebuilds it.
    assert refresh_step(wd, spec, formats=("png", "svg")) == "written"
    write_selection(wd / SELECTION_TSV, [_row(1, "a"), _row(2, "b")])
    assert refresh_step(wd, spec, formats=("png", "svg")) == "written"
    # A changed stage record rebuilds it too.
    cfg = Config.load(wd)
    cfg.stages["metadata"].params["release"] = "220"
    cfg.save(wd)
    assert refresh_step(wd, spec, formats=("png", "svg")) == "written"
    (wd / SELECTION_TSV).unlink()
    assert refresh_step(wd, spec, formats=("png", "svg")) == "removed"
    assert not (wd / "report" / "selection").exists()
    assert refresh_step(wd, spec, formats=("png", "svg")) == "absent"


def test_failed_collect_keeps_the_earlier_fragment(tmp_path: Path, monkeypatch) -> None:
    wd = _workdir(tmp_path)
    spec = BY_NAME["selection"]
    refresh_step(wd, spec, formats=("png",))
    before = (wd / "report" / "selection" / "step.json").read_text()
    write_selection(wd / SELECTION_TSV, [_row(1, "a")])

    def boom(inp):
        raise RuntimeError("no renderer")

    monkeypatch.setattr("repgenr.reporting.steps.selection.collect", boom)
    with pytest.raises(RuntimeError):
        refresh_step(wd, spec, formats=("png",))
    assert (wd / "report" / "selection" / "step.json").read_text() == before
    assert not any(p.name.startswith(".") for p in (wd / "report").iterdir())


def test_hook_never_raises(tmp_path: Path, monkeypatch, caplog) -> None:
    wd = _workdir(tmp_path)

    def boom(inp):
        raise RuntimeError("no renderer")

    monkeypatch.setattr("repgenr.reporting.steps.selection.collect", boom)
    with caplog.at_level(logging.WARNING):
        refresh_for_stage(wd, "metadata", LOG)
    assert "Could not write the report of step 'selection': no renderer" in caplog.text
    assert "selection" in STAGE_TO_STEPS["ingest"]
    refresh_for_stage(wd, "phylo_build_without_steps", LOG)  # no step: no-op


def test_html_numbering_escaping_and_no_external_resources(tmp_path: Path) -> None:
    rows = [_row(i, "<script>alert(1)</script>") for i in range(3)]
    wd = _workdir(tmp_path, rows)
    reporting.refresh(wd, logger=LOG, strict=True)
    fragments = reporting.collect_fragments([wd / "report"])
    prov = reporting.new_provenance(wd)
    out = reporting.write_html(wd / "report" / "report.html", fragments, prov, LOG)
    text = out.read_text()
    assert "<script>" not in text
    assert "&lt;script&gt;" in text
    assert re.findall(r'class="label">Table (\d+)\.', text) == ["1", "2", "3"]
    figures = re.findall(r'class="label">Figure (\d+)\.', text)
    assert figures == [str(i + 1) for i in range(len(figures))] and figures
    assert 'src="http' not in text and "<link" not in text
    for section in ("Summary", "Genome selection", "Methods", "Software and parameters"):
        assert f">{section}</h2>" in text
    assert "GTDB metadata table of release 226" in text
    assert "doi.org/10.1093/nar/gkab776" in text
    uri = re.search(r'href="data:text/tab-separated-values;base64,([^"]+)"', text).group(1)
    assert base64.b64decode(uri).decode().startswith("quantity\tgenomes")


def test_report_without_quality_values_has_a_note(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, [_row(1, "a", comp=None), _row(2, "b", comp=None)], record=False)
    reporting.refresh(wd, logger=LOG, strict=True)
    rep = read_fragment(wd / "report" / "selection")
    assert [f.id for f in rep.figures] == ["genomes_per_species"]
    assert rep.methods == []
    assert any("no completeness" in n for n in rep.notes)


def test_collect_fragments_from_several_folders(tmp_path: Path) -> None:
    wd = _workdir(tmp_path)
    reporting.refresh(wd, logger=LOG, strict=True)
    single = reporting.collect_fragments([wd / "report" / "selection"])
    assert [r.step for _, r in single] == ["selection"]
    assert reporting.collect_fragments([wd / "report"], ["selection"])[0][1].step == "selection"
    assert reporting.collect_fragments([wd / "report"], ["phylogeny"]) == []
