"""Unit test for the genome-download batching."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.stages import genome

_LOG = logging.getLogger("test")


def test_download_splits_into_fixed_size_batches(monkeypatch, tmp_path: Path) -> None:
    seen: list[list[str]] = []
    monkeypatch.setattr(genome, "_DOWNLOAD_BATCH_SIZE", 5)
    # hermetic: the real free-disk guard must not fail this test on a full host
    monkeypatch.setattr(genome, "_check_disk", lambda scratch, n, logger: None)

    def fake_batch(batch, filenames, dest_dir, scratch_dir, logger, keep_files, bi):
        seen.append(batch)
        return []  # no missing accessions

    monkeypatch.setattr(genome, "_download_one_batch", fake_batch)
    accs = [f"GCF_{i:03d}.1" for i in range(12)]
    genome.download_accessions(
        accs, {}, tmp_path / "dest", tmp_path / "scratch", _LOG, keep_files=False
    )
    # 12 accessions / batch 5 -> [5, 5, 2], covering every accession exactly once
    assert [len(b) for b in seen] == [5, 5, 2]
    assert [a for b in seen for a in b] == accs


# --- GenBank fallback -----------------------------------------------------------


def _fake_ncbi(monkeypatch, served: set[str]) -> list[list[str]]:
    """Serve only ``served``; returns the batches requested, in order."""
    requested: list[list[str]] = []
    monkeypatch.setattr(genome, "_check_disk", lambda scratch, n, logger: None)
    monkeypatch.setattr(genome, "require_reachable", lambda *_a, **_k: None)

    def fake_batch(batch, filenames, dest_dir, scratch_dir, logger, keep_files, bi, discarded=None):
        requested.append(list(batch))
        dest_dir.mkdir(parents=True, exist_ok=True)
        for acc in batch:
            if acc in served:
                (dest_dir / filenames[acc]).write_text(f">{acc}\nACGT\n")
        return [acc for acc in batch if acc not in served]

    monkeypatch.setattr(genome, "_download_one_batch", fake_batch)
    return requested


def _fallback(tmp_path: Path, accessions: list[str], **kwargs):
    filenames = {acc: f"F_g_s_{acc}.fasta" for acc in accessions}
    return genome.download_with_fallback(
        accessions, filenames, tmp_path / "genomes", tmp_path / "scratch", _LOG, **kwargs
    )


def test_unserved_refseq_accession_is_taken_from_its_genbank_pair(monkeypatch, tmp_path) -> None:
    requested = _fake_ncbi(monkeypatch, {"GCF_000001.1", "GCA_000002.1"})
    missing, substituted = _fallback(tmp_path, ["GCF_000001.1", "GCF_000002.1"])
    assert missing == [] and substituted == {"GCF_000002.1": "GCA_000002.1"}
    assert requested == [["GCF_000001.1", "GCF_000002.1"], ["GCA_000002.1"]]
    # Stored under the name of the selected accession.
    stored = tmp_path / "genomes" / "F_g_s_GCF_000002.1.fasta"
    assert stored.read_text().startswith(">GCA_000002.1")


def test_accession_stays_missing_when_the_pair_is_not_served_either(monkeypatch, tmp_path) -> None:
    _fake_ncbi(monkeypatch, set())
    missing, substituted = _fallback(tmp_path, ["GCF_000002.1"])
    assert missing == ["GCF_000002.1"] and substituted == {}


def test_fallback_can_be_turned_off(monkeypatch, tmp_path) -> None:
    requested = _fake_ncbi(monkeypatch, {"GCA_000002.1"})
    missing, substituted = _fallback(tmp_path, ["GCF_000002.1"], genbank_fallback=False)
    assert missing == ["GCF_000002.1"] and substituted == {}
    assert requested == [["GCF_000002.1"]]


def test_a_genbank_selection_is_not_retried(monkeypatch, tmp_path) -> None:
    requested = _fake_ncbi(monkeypatch, set())
    missing, substituted = _fallback(tmp_path, ["GCA_000003.1"])
    assert missing == ["GCA_000003.1"] and substituted == {} and len(requested) == 1


def test_the_recorded_pair_is_used_when_its_version_differs(monkeypatch, tmp_path) -> None:
    requested = _fake_ncbi(monkeypatch, {"GCA_000227705.3"})
    missing, substituted = _fallback(
        tmp_path, ["GCF_000227705.2"], twins={"GCF_000227705.2": "GCA_000227705.3"}
    )
    assert missing == [] and substituted == {"GCF_000227705.2": "GCA_000227705.3"}
    assert requested[1] == ["GCA_000227705.3"]


def test_a_pair_that_is_itself_selected_is_not_taken_twice(monkeypatch, tmp_path) -> None:
    requested = _fake_ncbi(monkeypatch, {"GCA_000002.1"})
    missing, substituted = _fallback(tmp_path, ["GCF_000002.1", "GCA_000002.1"])
    assert missing == ["GCF_000002.1"] and substituted == {} and len(requested) == 1


def test_genbank_twin() -> None:
    assert genome.genbank_twin("GCF_000001.2") == "GCA_000001.2"
    assert genome.genbank_twin("GCA_000001.2") is None
    assert genome.genbank_twin("NC_045512.2") is None
    # A recorded value that is no GenBank accession is not used.
    assert genome.genbank_twin("GCF_000001.2", {"GCF_000001.2": "none"}) == "GCA_000001.2"


def test_substitutions_round_trip_and_recorded_twins(tmp_path) -> None:
    path = tmp_path / "substituted_accessions.tsv"
    assert genome.read_substituted(path) == {}
    genome.write_substituted(path, {"GCF_2.1": "GCA_2.1", "GCF_1.1": "GCA_1.2"})
    assert path.read_text() == "selected\tdownloaded\nGCF_1.1\tGCA_1.2\nGCF_2.1\tGCA_2.1\n"
    assert genome.read_substituted(path) == {"GCF_1.1": "GCA_1.2", "GCF_2.1": "GCA_2.1"}
    assert genome.recorded_twins(tmp_path) == {}
    (tmp_path / "metadata_selected.tsv").write_text(
        "accession\taccession_ncbi\nGCF_1.1\tGCA_1.2\nGCF_3.1\t\n", encoding="utf-8"
    )
    assert genome.recorded_twins(tmp_path) == {"GCF_1.1": "GCA_1.2"}
