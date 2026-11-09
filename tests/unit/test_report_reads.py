"""Report steps of the run selection and of the assemblies."""

from __future__ import annotations

import random
from pathlib import Path

from repgenr.core.contracts import (
    ASSEMBLY_STATS_TSV,
    EXCUSED_RUNS_TSV,
    READS_TSV,
    AssemblyStatsRow,
    ExcusedRun,
    ReadRow,
    ReadsSketchRow,
    write_assembly_stats,
    write_excused_runs,
    write_reads,
    write_reads_sketches,
)
from repgenr.core.stage_stats import write_stats
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.methods import sentence
from repgenr.reporting.model import MethodsInfo
from repgenr.reporting.registry import BY_NAME
from repgenr.stages.reads import READS_STATS, READS_STATS_KIND


def _reads_wd(tmp_path: Path) -> Path:
    wd = tmp_path / "wd"
    wd.mkdir()
    rng = random.Random(1)
    rows = [
        ReadRow(
            f"SRR{i:06d}", f"SAMN{i}", "PRJ1", "Francisella", "263",
            "ILLUMINA" if i % 4 else "OXFORD_NANOPORE",
            "Illumina MiSeq" if i % 4 else "MinION",
            "PAIRED" if i % 4 else "SINGLE",
            int(rng.uniform(2e8, 5e9)), 1000,
        )
        for i in range(20)
    ]  # fmt: skip
    write_reads(wd / READS_TSV, rows)
    write_stats(
        wd / READS_STATS,
        READS_STATS_KIND,
        {"candidates": 50, "selected": 20,
         "funnel": [["whole-genome runs returned by ENA", 50],
                    ["at most 10,000,000,000 bases", 45], ["one run per biosample", 20]]},
    )  # fmt: skip
    sk = wd / "reads_sketches" / "s200k"
    sk.mkdir(parents=True)
    write_reads_sketches(
        sk / "reads_sketches.tsv",
        [ReadsSketchRow(r.run_accession, "ok" if i else "failed") for i, r in enumerate(rows)],
    )
    return wd


def test_reads_report(tmp_path: Path) -> None:
    wd = _reads_wd(tmp_path)
    assert refresh_step(wd, BY_NAME["reads"], formats=("png",)) == "written"
    frag = wd / "report" / "reads"
    rep = read_fragment(frag)
    assert (frag / "run_selection.tsv").read_text().splitlines()[1:] == [
        "whole-genome runs returned by ENA\t50",
        "at most 10,000,000,000 bases\t45",
        "one run per biosample\t20",
    ]
    by_platform = (frag / "runs_by_platform.tsv").read_text().splitlines()
    assert by_platform[1].startswith("ILLUMINA\tIllumina MiSeq\tPAIRED\t15\t")
    assert (frag / "sketch_status_s200k.tsv").read_text().splitlines()[1:] == [
        "ok\t19",
        "failed\t1",
    ]
    assert "20 sequencing runs were selected of 50" in rep.summary[0]
    assert [f.id for f in rep.figures] == ["bases_per_run"]


def test_assembly_report(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    rows = [
        AssemblyStatsRow(
            f"SRR{i}", f"x{i}.fasta", "spades" if i < 6 else "flye",
            40 + i, 1_900_000 + 1000 * i, 100_000 + 5000 * i, 300_000,
            est_coverage=50.0 + i, completeness=98.0 + i / 10, contamination=0.5,
            taxonomy_flag="" if i else "genus_mismatch",
        )
        for i in range(8)
    ]  # fmt: skip
    write_assembly_stats(wd / ASSEMBLY_STATS_TSV, rows)
    write_excused_runs(
        wd / EXCUSED_RUNS_TSV,
        [ExcusedRun("SRR9", "qc", "completeness below 90"), ExcusedRun("SRR10", "qc",
         "completeness below 90"), ExcusedRun("SRR11", "screen", "top genus Escherichia")],
    )  # fmt: skip
    refresh_step(wd, BY_NAME["assembly"], formats=("png",))
    frag = wd / "report" / "assembly"
    rep = read_fragment(frag)
    metrics = (frag / "assembly_metrics.tsv").read_text().splitlines()
    assert metrics[1].split("\t")[:3] == ["flye", "2", "46 (46-47)"]
    assert metrics[2].split("\t")[:2] == ["spades", "6"]
    excused = (frag / "excused_runs.tsv").read_text().splitlines()
    assert excused[1] == "qc\tcompleteness below 90\t2"
    flags = (frag / "taxonomy_flags.tsv").read_text().splitlines()
    assert "genus_mismatch\t1" in flags
    assert {f.id for f in rep.figures} == {"contiguity", "assembly_quality", "coverage"}
    assert "8 runs were assembled and 3 gave no genome" in rep.summary[0]


def test_reads_and_assemble_methods() -> None:
    reads = MethodsInfo(
        "reads",
        params={"target_genus": "Francisella", "platform": "illumina", "max_bases": 10**10,
                "drop_selection": ["MDA"], "one_per_sample": True},
    )  # fmt: skip
    text, _ = sentence(reads)
    assert text == (
        "Whole-genome sequencing runs of Francisella were selected from the European "
        "Nucleotide Archive (illumina runs; at most 10,000,000,000 bases; library "
        "selection not MDA). One run was kept per biosample."
    )
    asm = MethodsInfo(
        "assemble",
        tool="spades",
        tool_versions={"spades": "4.0.0", "checkm2": "1.0.2"},
        params={"assemblers_used": ["spades"], "checkm2_db": "/db", "polishers_used": []},
    )
    text, _ = sentence(asm)
    assert text.startswith("Runs were assembled with SPAdes 4.0.0.")
    assert "CheckM2 1.0.2" in text
