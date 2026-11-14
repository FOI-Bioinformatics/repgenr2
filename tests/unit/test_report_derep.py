"""Report steps of the dereplication and of glance."""

from __future__ import annotations

import logging
from datetime import UTC, datetime
from pathlib import Path

from repgenr.core.config import Config
from repgenr.core.contracts import (
    SELECTION_TSV,
    SelectionRow,
    write_clusters,
    write_genome_status,
    write_selection,
)
from repgenr.core.stage_stats import write_stats
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.methods import sentence
from repgenr.reporting.model import MethodsInfo
from repgenr.reporting.registry import BY_NAME
from repgenr.reporting.steps.dereplication import species_label
from repgenr.stages.dereplicate import DEREP_STATS, DEREP_STATS_KIND
from repgenr.stages.glance import GLANCE_STATS, GLANCE_STATS_KIND, GlanceParams, similarity_stats

LOG = logging.getLogger("test")


def _name(i: int) -> str:
    return f"Fam_Francisella_tularensis_GCF_{i:06d}.1.fasta"


def _derep_workdir(tmp_path: Path, *, tool: str = "drep", fail_qc: bool = True) -> Path:
    wd = tmp_path / "wd"
    (wd / "derep").mkdir(parents=True)
    rows = [
        SelectionRow(
            f"GCF_{i:06d}.1", "Fam", "Francisella", "tularensis" if i < 6 else "novicida",
            False, _name(i), 60.0 if i == 9 else 99.0, 1.0, False,
        )
        for i in range(10)
    ]  # fmt: skip
    write_selection(wd / SELECTION_TSV, rows)
    # Clusters: {0: 1-4}, {5}, {6: 7}; 8 and 9 failed QC (when fail_qc).
    write_clusters(
        wd / "derep" / "clusters.tsv",
        {_name(0): [_name(i) for i in range(1, 5)], _name(5): [], _name(6): [_name(7)]},
    )
    status = {_name(i): "contained" for i in (1, 2, 3, 4, 7)}
    status.update({_name(i): "representative" for i in (0, 5, 6)})
    reasons = {}
    if fail_qc:
        status.update({_name(8): "fail_qc", _name(9): "fail_qc"})
        reasons = {_name(8): "length", _name(9): "completeness,length"}
    write_genome_status(wd / "derep" / "genome_status.tsv", status, reasons)
    write_stats(
        wd / "derep" / DEREP_STATS,
        DEREP_STATS_KIND,
        {"tool": tool, "input_genomes": 10, "keeper_swaps": 1,
         "reduce": {"level": "species", "before": 4, "after": 3}},
    )  # fmt: skip
    cfg = Config.load(wd)
    cfg.record_stage(
        "dereplicate",
        tool=tool,
        params={"primary_ani": 0.9, "secondary_ani": 0.99, "keeper": "quality",
                "keeper_effective": "quality", "reduce": "species", "tool_args": {}},
        tool_versions={tool: "3.7.1"},
        completed=datetime.now(UTC).isoformat(),
    )  # fmt: skip
    cfg.save(wd)
    return wd


def test_dereplication_tables_and_figures(tmp_path: Path) -> None:
    wd = _derep_workdir(tmp_path)
    assert refresh_step(wd, BY_NAME["dereplication"], formats=("png",)) == "written"
    frag = wd / "report" / "dereplication"
    rep = read_fragment(frag)
    assert (frag / "genome_status.tsv").read_text().splitlines() == [
        "status\treason\tgenomes",
        "representative\t\t3",
        "contained\t\t5",
        "fail_qc\t\t2",
        "fail_qc\tcompleteness\t1",
        "fail_qc\tlength\t2",
    ]
    species = (frag / "species_reduction.tsv").read_text().splitlines()
    assert species[1:] == ["Francisella tularensis\t6\t2", "Francisella novicida\t4\t1"]
    assert {f.id for f in rep.figures} == {
        "cluster_sizes",
        "quality_by_status",
        "fail_qc_reasons",
        "species_reduction",
    }
    quality = next(f for f in rep.figures if f.id == "quality_by_status")
    assert "dRep quality filter" in quality.caption
    keys = {k.label: k.value for k in rep.key_numbers}
    assert keys["Genomes dereplicated"] == 10 and keys["Representatives"] == 3
    assert keys["Singleton clusters"] == 1 and keys["Largest cluster"] == 5
    text = rep.summary[0]
    assert "3 were kept as representatives (37.5% of the genomes that passed" in text
    assert "reduced 4 representatives to 3" in text
    assert "replaced the representative chosen by the tool in 1 clusters" in text


def test_no_fail_qc_figure_and_no_thresholds_for_other_tools(tmp_path: Path) -> None:
    wd = _derep_workdir(tmp_path, tool="skder", fail_qc=False)
    refresh_step(wd, BY_NAME["dereplication"], formats=("png",))
    rep = read_fragment(wd / "report" / "dereplication")
    ids = {f.id for f in rep.figures}
    assert "fail_qc_reasons" not in ids
    quality = next(f for f in rep.figures if f.id == "quality_by_status")
    assert "dRep" not in quality.caption


def test_dereplication_methods_sentence() -> None:
    m = MethodsInfo(
        "dereplicate",
        tool="drep",
        tool_versions={"drep": "quay.io/biocontainers/drep:3.7.1--pyhdfd78af_0"},
        params={"primary_ani": 0.9, "secondary_ani": 0.99, "keeper_effective": "gtdb",
                "reduce": "species", "target_reps": 0},
    )  # fmt: skip
    text, refs = sentence(m)
    assert text.startswith(
        "Genomes were dereplicated with dRep (container image quay.io/biocontainers/drep:"
    )
    assert "primary ANI threshold of 90% and a secondary ANI threshold of 99%" in text
    assert "'source' keeper rule" in text and "per species" in text
    assert refs == []


def test_species_label_keeps_a_viral_token() -> None:
    row = SelectionRow("A", "Picornaviridae", "Hepatovirus", "Hepatovirus-ahepa", False, "a.fa")
    assert species_label(row) == "Hepatovirus-ahepa"
    row = SelectionRow("B", "Fam", "Yersinia", "pestis", False, "b.fa")
    assert species_label(row) == "Yersinia pestis"


def test_glance_stats_and_report(tmp_path: Path) -> None:
    values = [0.95, 0.96, 0.97, 0.99, 0.999]
    stats = similarity_stats(values, GlanceParams(), "ANI", 4, "sourmash")
    assert stats["pairs"] == 5 and stats["median"] == 0.97
    assert sum(stats["histogram"]["counts"]) == 5
    assert len(stats["histogram"]["edges"]) == 101
    wd = _derep_workdir(tmp_path)
    write_stats(wd / GLANCE_STATS, GLANCE_STATS_KIND, stats)
    refresh_step(wd, BY_NAME["glance"], formats=("png",))
    rep = read_fragment(wd / "report" / "glance")
    assert "compared 4 genomes in 5 pairs" in rep.summary[0]
    hist = rep.figures[0]
    assert "ANI thresholds of the dereplication" in hist.caption
    row = (wd / "report" / "glance" / "pairwise_similarity.tsv").read_text().splitlines()[1]
    assert row.split("\t")[:3] == ["ANI", "4", "5"]


def test_glance_without_pairs(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    write_stats(
        wd / GLANCE_STATS,
        GLANCE_STATS_KIND,
        similarity_stats([], GlanceParams(), "ANI", 2, "drep"),
    )
    refresh_step(wd, BY_NAME["glance"], formats=("png",))
    rep = read_fragment(wd / "report" / "glance")
    assert rep.figures == [] and "no pairwise similarity" in rep.summary[0]
