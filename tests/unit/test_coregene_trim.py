"""Per-gene trimming of the core-gene route with ClipKIT (coregenes/trim.py)."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes import trim
from repgenr.coregenes.base import GeneAlignment

_LOG = logging.getLogger("test-coregene-trim")


def _fake_chain(chains: list[list[list[str]]], keep: dict[str, int] | None = None, rename=False):
    """clipkit keeps the first ``keep[gene file stem]`` columns (default: all
    but the last); headers are written as they are, or renamed with ``rename``."""

    def run_chain(caps, steps, **k):
        assert caps is trim.CLIPKIT
        argvs = [[str(c) for c in argv] for _, argv in steps]
        chains.append(argvs)
        for argv in argvs:
            src, dst = Path(argv[1]), Path(argv[argv.index("-o") + 1])
            cols = (keep or {}).get(dst.stem)
            out = []
            for line in src.read_text(encoding="utf-8").splitlines():
                if line.startswith(">"):
                    out.append(line + "_x" if rename else line)
                else:
                    out.append(line[:cols] if cols is not None else line[:-1])
            dst.write_text("\n".join(out) + "\n", encoding="utf-8")

    return run_chain


def _genes(tmp_path: Path, n: int, width: int = 4) -> list[GeneAlignment]:
    src = tmp_path / "builder"
    src.mkdir(exist_ok=True)
    out = []
    for i in range(n):
        path = src / f"g{i}.aln.fas"
        # Panaroo-style headers: <genome>;<gene id> and a description.
        path.write_text(
            f">G1;g{i}_1 a gene\n{'A' * width}\n>G2;g{i}_2\n{'C' * width}\n", encoding="utf-8"
        )
        out.append(GeneAlignment(f"g{i}", path))
    return out


def test_trims_in_batches_and_keeps_record_names(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(trim, "run_chain", _fake_chain(chains))
    genes = _genes(tmp_path, 120)
    out = trim.trim_gene_alignments(genes, tmp_path / "s", 4, _LOG)
    assert [ga.gene for ga in out] == [ga.gene for ga in genes]
    assert all(ga.path == tmp_path / "s" / "trimmed" / f"{ga.gene}.fasta" for ga in out)
    # One clipkit per gene, 50 genes per container start.
    assert [len(c) for c in chains] == [50, 50, 20]
    argv = chains[0][0]
    assert argv[0] == "clipkit" and argv[1] == str(genes[0].path)
    assert argv[argv.index("-m") + 1] == "smart-gap"
    # dna: codon-preserving trimming of a nucleotide alignment.
    assert argv[-3:] == ["-s", "nt", "-co"]
    text = out[0].path.read_text(encoding="utf-8")
    assert text == ">G1;g0_1 a gene\nAAA\n>G2;g0_2\nCCC\n"
    assert trim.alignment_width(out[0].path) == 3


def test_a_gene_trimmed_to_nothing_is_dropped_and_counted(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    monkeypatch.setattr(trim, "run_chain", _fake_chain([], keep={"g1": 0}))
    genes = _genes(tmp_path, 3)
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        out = trim.trim_gene_alignments(genes, tmp_path / "s", 1, _LOG)
    assert [ga.gene for ga in out] == ["g0", "g2"]
    assert "1 gene(s) were dropped because ClipKIT removed every column: g1" in caplog.text


def test_mode_is_validated_and_passed(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(trim, "run_chain", _fake_chain(chains))
    trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s", 1, _LOG, mode="kpic-smart-gap")
    assert chains[0][0][-5:] == ["-m", "kpic-smart-gap", "-s", "nt", "-co"]
    with pytest.raises(UserInputError, match="Unknown ClipKIT mode 'c3'"):
        trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s", 1, _LOG, mode="c3")
    assert set(trim.TRIM_MODES) == {"smart-gap", "gappy", "kpic-smart-gap", "kpi-smart-gap"}


def test_renamed_records_or_a_missing_output_are_refused(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(trim, "run_chain", _fake_chain([], rename=True))
    with pytest.raises(WorkdirError, match="changed the record names of gene g0"):
        trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s", 1, _LOG)
    monkeypatch.setattr(trim, "run_chain", lambda *a, **k: None)
    with pytest.raises(WorkdirError, match="wrote no trimmed alignment for gene g0"):
        trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s2", 1, _LOG)


def test_gene_ids_become_safe_unique_file_names() -> None:
    assert trim._file_names(["a/b", "..x", "famA", "FAMA", ""]) == [
        "a_b",
        "x",
        "famA",
        "FAMA.3",
        "gene",
    ]
    # A suffixed name never collides with a later gene id.
    assert trim._file_names(["a", "A", "a.1"]) == ["a", "A.1", "a.1.2"]


def test_check_gene_trim_and_the_capabilities() -> None:
    assert trim.check_gene_trim("clipkit") == "clipkit"
    with pytest.raises(UserInputError, match="--gene-trim must be one of none, clipkit"):
        trim.check_gene_trim("trimal")
    assert trim.CLIPKIT.container == "quay.io/biocontainers/clipkit:2.14.0--pyhdfd78af_0"
    (spec,) = trim.CLIPKIT.required_binaries
    assert spec.name == "clipkit" and spec.version_args == ("--version",)
    assert spec.min_version == "2.1.2"


def test_protein_alignments_are_trimmed_as_amino_acids(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(trim, "run_chain", _fake_chain(chains))
    trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s", 1, _LOG, alphabet="protein")
    argv = chains[0][0]
    assert argv[-2:] == ["-s", "aa"] and "-co" not in argv
    with pytest.raises(ValueError, match="alphabet"):
        trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s", 1, _LOG, alphabet="rna")


def test_an_output_without_records_is_refused_not_dropped(tmp_path: Path, monkeypatch) -> None:
    def run_chain(caps, steps, **k):
        for _, argv in steps:
            Path(str(argv[argv.index("-o") + 1])).write_text("", encoding="utf-8")

    monkeypatch.setattr(trim, "run_chain", run_chain)
    with pytest.raises(WorkdirError, match="without records for gene g0"):
        trim.trim_gene_alignments(_genes(tmp_path, 1), tmp_path / "s", 1, _LOG)


def test_a_failing_batch_propagates(tmp_path: Path, monkeypatch) -> None:
    from repgenr.core.errors import ToolExecutionError

    def run_chain(caps, steps, **k):
        raise ToolExecutionError(["clipkit"], 1)

    monkeypatch.setattr(trim, "run_chain", run_chain)
    # 120 genes in three batches on two threads: the error reaches the caller.
    with pytest.raises(ToolExecutionError):
        trim.trim_gene_alignments(_genes(tmp_path, 120), tmp_path / "s", 2, _LOG)
