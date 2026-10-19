"""reads-sketch against SRA: stream two public runs through fastq-dump into
sourmash, with a spot cap, and repeat without streaming again."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.contracts import READS_SKETCHES_TSV, read_reads_sketches

pytestmark = [
    pytest.mark.live,
    pytest.mark.network,
    pytest.mark.requires_binary("fastq-dump"),
    pytest.mark.requires_binary("sourmash"),
]


def test_reads_sketch_streams_two_runs(run_repgenr, tmp_path: Path) -> None:
    """SRR25474756 (MiSeq paired) and SRR28800588 (ONT): 5000 spots each,
    about 15 s per run on a home connection (measured 2026-10-19)."""
    out = tmp_path / "sk"
    result = run_repgenr(
        "reads-sketch",
        "-o",
        out,
        "--accession",
        "SRR25474756",
        "--accession",
        "SRR28800588",
        "--accession",
        "SRR999999999",
        "--max-spots",
        "5000",
    )
    variant = out / "spots_5000"
    rows = {r.run_accession: r for r in read_reads_sketches(variant / READS_SKETCHES_TSV)}
    paired, long_read, bad = rows["SRR25474756"], rows["SRR28800588"], rows["SRR999999999"]
    assert paired.status == "ok" and paired.spots_read == 5000 and paired.sequences == 10000
    assert long_read.status == "ok" and long_read.sequences == 5000
    assert bad.status == "failed" and bad.reason.startswith("fastq-dump exited")
    assert (variant / "SRR25474756.sig.zip").is_file()
    assert not (variant / "SRR999999999.sig.zip").exists()
    assert "Sketched 2 run(s), 1 failed" in result.stderr + result.stdout

    again = run_repgenr(
        "reads-sketch", "-o", out, "--accession", "SRR25474756", "--max-spots", "5000"
    )
    assert "already sketched" in again.stderr + again.stdout
