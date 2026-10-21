"""Shared orthogroup alignment: protein files, mafft batches, back-translation."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.coregenes import codon
from repgenr.coregenes.orthofinder import OrthoFinderBuilder

_LOG = logging.getLogger("test-codon")
_CAPS = OrthoFinderBuilder.capabilities


def _sidecars(tmp_path: Path) -> tuple[dict[str, Path], dict[str, Path]]:
    """Three genomes; gene k (MKLM) and gene b (AC), ids suffixed by genome."""
    faas: dict[str, Path] = {}
    ffns: dict[str, Path] = {}
    for genome in ("G1", "G2", "G3"):
        faa = tmp_path / f"{genome}.faa"
        ffn = tmp_path / f"{genome}.ffn"
        faa.write_text(f">k_{genome} product\nMKLM\n>b_{genome}\nAC\n", encoding="utf-8")
        ffn.write_text(f">k_{genome}\nATGAAACTGATGTAA\n>b_{genome}\nGCTTGT\n", encoding="utf-8")
        faas[genome], ffns[genome] = faa, ffn
    return faas, ffns


def _fake_chain(chains: list[list[list[str]]], gap: bool = False):
    """mafft copies its input; with ``gap`` it adds one gap column at the end."""

    def run_chain(caps, steps, **k):
        argvs = [[str(c) for c in argv] for _, argv in steps]
        chains.append(argvs)
        for argv in argvs:
            src, dst = Path(argv[-1]), Path(argv[argv.index("--out") + 1])
            lines = src.read_text(encoding="utf-8").splitlines()
            out = [line if line.startswith(">") or not gap else line + "-" for line in lines]
            dst.write_text("\n".join(out) + "\n", encoding="utf-8")

    return run_chain


def _records(path: Path) -> dict[str, str]:
    parts = path.read_text(encoding="utf-8").split(">")[1:]
    return {p.split("\n")[0]: p.split("\n")[1] for p in parts}


def test_aligns_and_back_translates(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(codon, "run_chain", _fake_chain(chains, gap=True))
    faas, ffns = _sidecars(tmp_path)
    selected = [
        ("OG1", {"G1": "k_G1", "G2": "k_G2", "G3": "k_G3"}),
        ("OG2", {"G2": "b_G2", "G3": "b_G3"}),
    ]
    out, mismatched = codon.align_orthogroups(
        selected, faas, ffns, tmp_path / "o", _CAPS, 2, _LOG, "orthofinder"
    )
    assert mismatched == 0
    assert [(ga.gene, ga.path.name) for ga in out] == [("OG1", "OG1.fasta"), ("OG2", "OG2.fasta")]
    assert all(ga.path.parent == tmp_path / "o" / "genes" for ga in out)
    # The stop codon is dropped and a protein gap is three nucleotide gaps.
    assert _records(out[0].path) == {g: "ATGAAACTGATG---" for g in ("G1", "G2", "G3")}
    # Members only: G1 is not in OG2, so it has no record.
    assert _records(out[1].path) == {"G2": "GCTTGT---", "G3": "GCTTGT---"}
    # Proteins named by genome, members only.
    assert _records(tmp_path / "o" / "orthogroups" / "OG2.faa") == {"G2": "AC", "G3": "AC"}
    argvs = [a for chain in chains for a in chain]
    assert len(argvs) == 2 and all(a[0] == "mafft" and "--auto" in a for a in argvs)
    # One thread per mafft; the batches run in parallel instead.
    assert all(a[a.index("--thread") + 1] == "1" for a in argvs)


def test_a_length_mismatch_is_counted_and_left_out(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(codon, "run_chain", _fake_chain([]))
    faas, ffns = _sidecars(tmp_path)
    ffns["G3"].write_text(">k_G3\nATGAA\n>b_G3\nGCTTGT\n", encoding="utf-8")
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        out, mismatched = codon.align_orthogroups(
            [("OG1", {"G1": "k_G1", "G2": "k_G2", "G3": "k_G3"})],
            faas, ffns, tmp_path / "o", _CAPS, 1, _LOG, "scarap",
        )  # fmt: skip
    assert mismatched == 1 and "1 gene(s) were left out" in caplog.text
    assert set(_records(out[0].path)) == {"G1", "G2"}


def test_a_single_member_is_copied_not_aligned(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(codon, "run_chain", _fake_chain(chains))
    faas, ffns = _sidecars(tmp_path)
    out, _ = codon.align_orthogroups(
        [("OG1", {"G2": "k_G2"})], faas, ffns, tmp_path / "o", _CAPS, 1, _LOG, "scarap"
    )
    assert chains == []
    assert _records(out[0].path) == {"G2": "ATGAAACTGATG"}


def test_selenocysteine_is_sent_as_x(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(codon, "run_chain", _fake_chain([]))
    faas, ffns = _sidecars(tmp_path)
    faas["G1"].write_text(">k_G1\nMKUM\n", encoding="utf-8")
    out, mismatched = codon.align_orthogroups(
        [("OG1", {"G1": "k_G1", "G2": "k_G2"})], faas, ffns, tmp_path / "o", _CAPS, 1, _LOG, "x"
    )
    assert _records(tmp_path / "o" / "orthogroups" / "OG1.faa")["G1"] == "MKXM"
    assert mismatched == 0 and _records(out[0].path)["G1"] == "ATGAAACTGATG"


def test_refusals(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(codon, "run_chain", _fake_chain([]))
    faas, ffns = _sidecars(tmp_path)
    with pytest.raises(WorkdirError, match="scarap placed gene 'zz' of G1 in OG1"):
        codon.align_orthogroups(
            [("OG1", {"G1": "zz", "G2": "k_G2"})], faas, ffns, tmp_path / "a", _CAPS, 1, _LOG,
            "scarap",
        )  # fmt: skip
    with pytest.raises(WorkdirError, match="none of the input genomes"):
        codon.align_orthogroups(
            [("OG1", {"G9": "k_G1"})], faas, ffns, tmp_path / "b", _CAPS, 1, _LOG, "scarap"
        )
    monkeypatch.setattr(codon, "run_chain", lambda caps, steps, **k: None)
    with pytest.raises(WorkdirError, match="mafft wrote no alignment for OG1"):
        codon.align_orthogroups(
            [("OG1", {"G1": "k_G1", "G2": "k_G2"})], faas, ffns, tmp_path / "c", _CAPS, 1, _LOG,
            "scarap",
        )  # fmt: skip


def test_protein_alphabet_stops_after_mafft(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(codon, "run_chain", _fake_chain(chains, gap=True))
    faas, ffns = _sidecars(tmp_path)
    # The gene files are not read: a frame-shifted gene would not matter.
    ffns["G3"].write_text(">k_G3\nATGAA\n>b_G3\nGC\n", encoding="utf-8")
    selected = [
        ("OG1", {"G1": "k_G1", "G2": "k_G2", "G3": "k_G3"}),
        ("OG2", {"G2": "b_G2"}),
    ]
    out, mismatched = codon.align_orthogroups(
        selected, faas, ffns, tmp_path / "o", _CAPS, 2, _LOG, "orthofinder", alphabet="protein"
    )
    assert mismatched == 0
    assert [(ga.gene, ga.path) for ga in out] == [
        ("OG1", tmp_path / "o" / "alignments" / "OG1.aln"),
        ("OG2", tmp_path / "o" / "alignments" / "OG2.aln"),
    ]
    # The mafft protein alignment, records named by genome.
    assert _records(out[0].path) == {g: "MKLM-" for g in ("G1", "G2", "G3")}
    # A single member is its own alignment, not sent to mafft.
    assert _records(out[1].path) == {"G2": "AC"}
    assert not (tmp_path / "o" / "genes").exists() or not any((tmp_path / "o" / "genes").iterdir())
    assert sum(len(chain) for chain in chains) == 1


def test_an_unknown_alphabet_is_a_programming_error(tmp_path: Path) -> None:
    faas, ffns = _sidecars(tmp_path)
    with pytest.raises(ValueError, match="alphabet"):
        codon.align_orthogroups(
            [], faas, ffns, tmp_path / "o", _CAPS, 1, _LOG, "orthofinder", alphabet="rna"
        )


def test_an_empty_protein_alignment_is_refused(tmp_path: Path, monkeypatch) -> None:
    def run_chain(caps, steps, **k):
        for _, argv in steps:
            argv = [str(c) for c in argv]
            # A file that is not empty but holds no record.
            Path(argv[argv.index("--out") + 1]).write_text("\n", encoding="utf-8")

    monkeypatch.setattr(codon, "run_chain", run_chain)
    faas, ffns = _sidecars(tmp_path)
    selected = [("OG1", {"G1": "k_G1", "G2": "k_G2"})]
    with pytest.raises(WorkdirError, match="empty alignment for OG1"):
        codon.align_orthogroups(
            selected, faas, ffns, tmp_path / "o", _CAPS, 1, _LOG, "orthofinder", alphabet="protein"
        )
