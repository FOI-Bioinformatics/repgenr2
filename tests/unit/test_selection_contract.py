"""Selection contract (metadata -> genome data-channel hand-off) and genome-fetch."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import (
    SelectionRow,
    genome_filename,
    read_selection,
    write_selection,
)
from repgenr.stages import genome as genome_stage
from repgenr.stages.genome_steps import GenomeFetchParams, genome_fetch

_LOG = logging.getLogger("test")


def test_genome_filename_is_canonical() -> None:
    assert (
        genome_filename("Fam", "Gen", "spec", "GCF_000001.1") == "Fam_Gen_spec_GCF_000001.1.fasta"
    )


def test_selection_round_trip(tmp_path: Path) -> None:
    rows = [
        SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", True, "Fam_Gen_sp2_GCF_2.1.fasta"),
    ]
    path = tmp_path / "selection.tsv"
    write_selection(path, rows)
    back = read_selection(path)
    assert back == rows
    # the outgroup flag survives the round trip
    assert [r.is_outgroup for r in back] == [False, True]


def test_genome_fetch_downloads_selected_and_outgroup(monkeypatch, tmp_path: Path) -> None:
    calls: list[tuple[list[str], Path]] = []

    def fake_download(
        accessions, filenames, dest_dir, scratch_dir, logger, keep_files=False, discarded=None
    ):
        dest_dir.mkdir(parents=True, exist_ok=True)
        for acc in accessions:
            (dest_dir / filenames[acc]).write_text(">x\nACGT\n")
        calls.append((list(accessions), dest_dir))
        return []

    # genome-fetch downloads through genome.download_with_fallback.
    monkeypatch.setattr("repgenr.stages.genome.download_accessions", fake_download)
    monkeypatch.setattr("repgenr.stages.genome_steps.preflight", lambda _caps: {"datasets": "0"})

    rows = [
        SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
        SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", False, "Fam_Gen_sp2_GCF_2.1.fasta"),
        SelectionRow("GCF_9.1", "Fam", "Out", "grp", True, "Fam_Out_grp_GCF_9.1.fasta"),
    ]
    selection = tmp_path / "selection.tsv"
    write_selection(selection, rows)

    out = tmp_path / "out"
    n = genome_fetch(GenomeFetchParams(selection_tsv=selection, out_dir=out), _LOG)

    assert n == 2  # selected count (outgroup excluded)
    genomes = sorted(p.name for p in (out / "genomes").iterdir())
    assert genomes == ["Fam_Gen_sp1_GCF_1.1.fasta", "Fam_Gen_sp2_GCF_2.1.fasta"]
    assert (out / "outgroup" / "Fam_Out_grp_GCF_9.1.fasta").exists()
    # selected and outgroup went to separate destination directories
    dests = {dest.name for _, dest in calls}
    assert dests == {"genomes", "outgroup"}
    # every accession arrived: the record exists and is empty
    assert (out / "missing_accessions.txt").read_text() == ""


def _fetch_with_missing(monkeypatch, tmp_path: Path, not_served: set[str]) -> Path:
    def fake_download(
        accessions, filenames, dest_dir, scratch_dir, logger, keep_files=False, discarded=None
    ):
        dest_dir.mkdir(parents=True, exist_ok=True)
        for acc in accessions:
            if acc not in not_served:
                (dest_dir / filenames[acc]).write_text(">x\nACGT\n")
        return [acc for acc in accessions if acc in not_served]

    monkeypatch.setattr("repgenr.stages.genome.download_accessions", fake_download)
    monkeypatch.setattr("repgenr.stages.genome_steps.preflight", lambda _caps: {"datasets": "0"})
    selection = tmp_path / "selection.tsv"
    write_selection(
        selection,
        [
            SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta"),
            SelectionRow("GCF_2.1", "Fam", "Gen", "sp2", False, "Fam_Gen_sp2_GCF_2.1.fasta"),
            SelectionRow("GCF_9.1", "Fam", "Out", "grp", True, "Fam_Out_grp_GCF_9.1.fasta"),
        ],
    )
    return selection


def test_genome_fetch_records_accessions_ncbi_did_not_serve(monkeypatch, tmp_path: Path) -> None:
    selection = _fetch_with_missing(monkeypatch, tmp_path, {"GCF_2.1", "GCA_2.1"})
    out = tmp_path / "out"
    genome_fetch(GenomeFetchParams(selection_tsv=selection, out_dir=out), _LOG)
    assert (out / "missing_accessions.txt").read_text() == "GCF_2.1\n"
    assert (out / "substituted_accessions.tsv").read_text() == "selected\tdownloaded\n"


def test_genome_fetch_takes_the_genbank_pair_of_an_unserved_accession(
    monkeypatch, tmp_path: Path
) -> None:
    # Neither GCF_2.1 nor the outgroup GCF_9.1 is served; their pairs are.
    selection = _fetch_with_missing(monkeypatch, tmp_path, {"GCF_2.1", "GCF_9.1"})
    out = tmp_path / "out"
    genome_fetch(GenomeFetchParams(selection_tsv=selection, out_dir=out), _LOG)
    assert (out / "genomes" / "Fam_Gen_sp2_GCF_2.1.fasta").exists()
    assert (out / "outgroup" / "Fam_Out_grp_GCF_9.1.fasta").exists()
    assert (out / "missing_accessions.txt").read_text() == ""
    assert (out / "substituted_accessions.tsv").read_text() == (
        "selected\tdownloaded\nGCF_2.1\tGCA_2.1\nGCF_9.1\tGCA_9.1\n"
    )


def test_genome_fetch_without_the_fallback(monkeypatch, tmp_path: Path) -> None:
    selection = _fetch_with_missing(monkeypatch, tmp_path, {"GCF_2.1"})
    out = tmp_path / "out"
    genome_fetch(
        GenomeFetchParams(selection_tsv=selection, out_dir=out, genbank_fallback=False), _LOG
    )
    assert (out / "missing_accessions.txt").read_text() == "GCF_2.1\n"


def test_genome_fetch_stops_when_the_outgroup_is_not_served(monkeypatch, tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    selection = _fetch_with_missing(monkeypatch, tmp_path, {"GCF_9.1", "GCA_9.1"})
    with pytest.raises(WorkdirError, match="outgroup GCF_9.1"):
        genome_fetch(GenomeFetchParams(selection_tsv=selection, out_dir=tmp_path / "out"), _LOG)


def test_genome_fetch_rejects_missing_selection(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    with pytest.raises(WorkdirError):
        genome_fetch(
            GenomeFetchParams(selection_tsv=tmp_path / "nope.tsv", out_dir=tmp_path / "o"), _LOG
        )


def test_genome_fetch_checks_datasets_before_downloading(monkeypatch, tmp_path: Path) -> None:
    # Without --versions-out the step must still run the preflight, so a missing
    # datasets binary is reported as MissingBinaryError (exit 4), not as a
    # FileNotFoundError from the first download (exit 1 with a traceback).
    from repgenr.core.errors import MissingBinaryError

    def missing(_caps):
        raise MissingBinaryError("datasets: not found on PATH")

    def no_download(*_args, **_kwargs):
        raise AssertionError("download attempted before the preflight")

    monkeypatch.setattr("repgenr.stages.genome_steps.preflight", missing)
    monkeypatch.setattr("repgenr.stages.genome.download_accessions", no_download)
    selection = tmp_path / "selection.tsv"
    write_selection(
        selection,
        [SelectionRow("GCF_1.1", "Fam", "Gen", "sp1", False, "Fam_Gen_sp1_GCF_1.1.fasta")],
    )
    with pytest.raises(MissingBinaryError):
        genome_fetch(GenomeFetchParams(selection_tsv=selection, out_dir=tmp_path / "o"), _LOG)


def test_output_name_matches_canonical(tmp_path: Path) -> None:
    # the legacy genome stage and the canonical helper agree on filenames
    from repgenr.core.manifest import GenomeRecord

    rec = GenomeRecord(
        accession="GCF_5.1",
        source="gtdb",
        is_outgroup=False,
        family="Fam",
        genus="Gen",
        species="sp",
    )
    assert genome_stage._output_name(rec) == genome_filename("Fam", "Gen", "sp", "GCF_5.1")


def test_selection_quality_columns_round_trip(tmp_path):
    from repgenr.core.contracts import SelectionRow, read_selection, write_selection

    rows = [
        SelectionRow(
            "A",
            "f",
            "g",
            "s",
            False,
            "f_g_s_A.fasta",
            completeness=98.5,
            contamination=1.0,
        ),
        SelectionRow("B", "f", "g", "s", True, "f_g_s_B.fasta"),
    ]
    write_selection(tmp_path / "s.tsv", rows)
    back = read_selection(tmp_path / "s.tsv")
    assert (back[0].completeness, back[0].contamination) == (98.5, 1.0)
    assert (back[1].completeness, back[1].contamination) == (None, None)


def test_selection_without_quality_columns_still_reads(tmp_path):
    from repgenr.core.contracts import read_selection

    p = tmp_path / "s.tsv"
    p.write_text(
        "accession\tfamily\tgenus\tspecies\tis_outgroup\tfilename\nA\tf\tg\ts\t0\tx.fasta\n"
    )
    assert read_selection(p)[0].completeness is None


def test_read_selection_names_missing_columns(tmp_path: Path) -> None:
    # A table that is not a selection.tsv is an input error naming the file and
    # the absent columns (exit 3), not a KeyError traceback (exit 1).
    from repgenr.core.errors import WorkdirError

    path = tmp_path / "selection.tsv"
    path.write_text("a\tb\nx\ty\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match=r"selection\.tsv.*accession.*filename"):
        read_selection(path)
