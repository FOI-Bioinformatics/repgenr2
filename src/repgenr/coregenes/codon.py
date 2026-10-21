"""Protein alignment and back-translation of orthogroups, for builders that
cluster proteins but hand over DNA gene alignments.

OrthoFinder and SCARAP give orthogroups of genes, not alignments the stage
can use: OrthoFinder with ``-og`` aligns nothing, and SCARAP's ``concat``
deletes its per-orthogroup alignments once it has concatenated them. This
module aligns the proteins of each selected orthogroup with mafft
(``--auto``) and back-translates the alignment to codons from the ``.ffn``
file the annotate stage wrote, so each core gene becomes one aligned DNA
FASTA, named by genome, that the stage measures, ranks and concatenates.

The alignments run in batches, one ``mafft`` per orthogroup, each batch in
one container when a container backend is active. A genome whose gene's
nucleotide length does not match its protein (a partial or frame-shifted
gene) is left out of that gene and counted.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path

from ..core.containers import run_chain
from ..core.errors import WorkdirError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities
from .base import GeneAlignment

# Orthogroups aligned per container start.
_BATCH = 50
# mafft refuses selenocysteine (U) and pyrrolysine (O), which bakta writes into
# the protein of a recoded gene; X keeps the length, so the codons still map.
_MAFFT_UNSAFE = str.maketrans({"U": "X", "O": "X", "u": "x", "o": "x"})


def align_orthogroups(
    selected: Sequence[tuple[str, dict[str, str]]],
    faas: Mapping[str, Path],
    ffns: Mapping[str, Path],
    out_dir: Path,
    capabilities: ToolCapabilities,
    threads: int,
    logger: logging.Logger,
    tool: str,
) -> tuple[list[GeneAlignment], int]:
    """Align each orthogroup of ``selected`` and back-translate it to codons.

    ``selected`` holds (orthogroup, {genome: gene id}); ``faas`` and ``ffns``
    map each genome to its protein and gene FASTA. Writes the proteins to
    ``out_dir/orthogroups/<og>.faa``, the protein alignments to
    ``out_dir/alignments/<og>.aln`` and the codon alignments to
    ``out_dir/genes/<og>.fasta``, records named by genome. A single-member
    orthogroup is its own alignment and is not sent to mafft. Returns the
    codon alignments in the order of ``selected`` and the number of genes
    left out because their nucleotide length does not match the protein.
    """
    wanted: dict[str, set[str]] = {genome: set() for genome in faas}
    for og, members in selected:
        for genome, gene in members.items():
            if genome not in wanted or genome not in ffns:
                raise WorkdirError(
                    f"{tool} placed a gene of {genome} in {og}, but {genome} is none of the "
                    "input genomes."
                )
            wanted[genome].add(gene)
    proteins = {g: _read_fasta(faas[g], ids) for g, ids in wanted.items() if ids}
    nucleotides = {g: _read_fasta(ffns[g], ids) for g, ids in wanted.items() if ids}

    prot_dir = out_dir / "orthogroups"
    aln_dir = out_dir / "alignments"
    gene_dir = out_dir / "genes"
    for d in (prot_dir, aln_dir, gene_dir):
        d.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[str, Path, Path]] = []
    to_align: list[tuple[str, Path, Path]] = []
    for og, members in selected:
        src = prot_dir / f"{og}.faa"
        dst = aln_dir / f"{og}.aln"
        with open(src, "w", encoding="utf-8") as fo:
            for genome, gene in members.items():
                protein = proteins[genome].get(gene)
                if protein is None:
                    raise WorkdirError(
                        f"{tool} placed gene {gene!r} of {genome} in {og}, but the protein "
                        "file has no record of that name; the annotation's identifiers are "
                        "probably not unique in their first word."
                    )
                fo.write(f">{genome}\n{protein.translate(_MAFFT_UNSAFE)}\n")
        jobs.append((og, src, dst))
        if len(members) > 1:
            to_align.append((og, src, dst))
        else:
            # One sequence is its own alignment (possible at a low threshold).
            shutil.copyfile(src, dst)
    _align(capabilities, to_align, out_dir, threads, logger)

    out: list[GeneAlignment] = []
    mismatched = 0
    for (og, _, dst), (_, members) in zip(jobs, selected, strict=True):
        if not dst.is_file() or dst.stat().st_size == 0:
            raise WorkdirError(f"mafft wrote no alignment for {og}")
        aligned = _read_alignment(dst)
        if not aligned:
            raise WorkdirError(f"mafft wrote an empty alignment for {og}")
        path = gene_dir / f"{og}.fasta"
        with open(path, "w", encoding="utf-8") as fo:
            for genome, gene in members.items():
                nt = nucleotides[genome].get(gene)
                if genome not in aligned or nt is None:
                    continue
                codons = _codons(aligned[genome], nt)
                if codons is None:
                    mismatched += 1
                    continue
                fo.write(f">{genome}\n{codons}\n")
        out.append(GeneAlignment(og, path))
    return out, mismatched


def _align(
    capabilities: ToolCapabilities,
    jobs: Sequence[tuple[str, Path, Path]],
    out_dir: Path,
    threads: int,
    logger: logging.Logger,
) -> None:
    batches = [jobs[i : i + _BATCH] for i in range(0, len(jobs), _BATCH)]

    def run(batch: Sequence[tuple[str, Path, Path]]) -> None:
        run_chain(
            capabilities,
            [
                ("mafft", ["mafft", "--quiet", "--auto", "--thread", "1", "--out", dst, src])
                for _, src, dst in batch
            ],
            logger=logger,
            cwd=out_dir,
        )

    logger.info("Aligning %d orthogroups with mafft", len(jobs))
    # The first failure cancels the batches not yet started.
    parallel_map(run, batches, max(1, threads), logger=logger)


def _read_fasta(path: Path, wanted: set[str]) -> dict[str, str]:
    """The records of ``path`` whose identifier is in ``wanted``."""
    seqs: dict[str, list[str]] = {}
    current: list[str] | None = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ident = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                current = seqs.setdefault(ident, []) if ident in wanted else None
            elif current is not None:
                current.append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _read_alignment(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                current = seqs.setdefault(line[1:].strip(), [])
            else:
                current.append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _codons(aligned: str, nucleotides: str) -> str | None:
    """``aligned`` (a gapped protein) as codons from ``nucleotides``.

    None when the lengths do not match, allowing one trailing stop codon.
    """
    residues = len(aligned) - aligned.count("-")
    nt = nucleotides.upper()
    if len(nt) == 3 * residues + 3:
        nt = nt[:-3]
    if len(nt) != 3 * residues:
        return None
    out: list[str] = []
    i = 0
    for aa in aligned:
        if aa == "-":
            out.append("---")
        else:
            out.append(nt[i : i + 3])
            i += 3
    return "".join(out)
