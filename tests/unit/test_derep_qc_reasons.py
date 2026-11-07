"""The reason of a fail_qc genome survives chunked runs and the step contract."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.contracts import read_genome_status_reasons
from repgenr.dereplicators.base import DerepResult
from repgenr.stages.derep_steps import _load_chunk, _write_step_contract
from repgenr.stages.dereplicate import _compose_two_stage


def test_compose_keeps_reasons_and_names_members_of_a_rejected_rep(tmp_path) -> None:
    chunk = DerepResult(
        representatives=[tmp_path / "a.fasta", tmp_path / "c.fasta"],
        clusters={"a.fasta": ["b.fasta"], "c.fasta": ["d.fasta"]},
        genome_status={
            "a.fasta": "representative",
            "b.fasta": "contained",
            "c.fasta": "representative",
            "d.fasta": "contained",
            "e.fasta": "fail_qc",
        },
        qc_reasons={"e.fasta": "length"},
    )
    # The later pass rejects c.fasta, which takes its member d.fasta with it.
    final = DerepResult(
        representatives=[tmp_path / "a.fasta"],
        clusters={"a.fasta": []},
        genome_status={"a.fasta": "representative", "c.fasta": "fail_qc"},
        qc_reasons={"c.fasta": "completeness"},
    )
    out = _compose_two_stage([chunk], final)
    assert out.genome_status["d.fasta"] == "fail_qc"
    assert out.qc_reasons == {
        "e.fasta": "length",
        "c.fasta": "completeness",
        "d.fasta": "representative_fail_qc",
    }


def test_step_contract_round_trips_reasons(tmp_path: Path) -> None:
    rep = tmp_path / "a.fasta"
    rep.write_text(">a\nACGT\n")
    result = DerepResult(
        representatives=[rep],
        clusters={"a.fasta": []},
        genome_status={"a.fasta": "representative", "b.fasta": "fail_qc"},
        qc_reasons={"b.fasta": "contamination"},
    )
    out = tmp_path / "contract"
    _write_step_contract(out, result, [], {})
    assert read_genome_status_reasons(out / "genome_status.tsv") == {"b.fasta": "contamination"}
    assert _load_chunk(out).qc_reasons == {"b.fasta": "contamination"}
