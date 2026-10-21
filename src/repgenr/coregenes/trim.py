"""Optional trimming of the per-gene alignments of the core-gene route.

Across a genus or wider set, the ends of a gene alignment are often ragged
(partial genes, differing start codons) and some columns are poorly
aligned. Those columns add variable and informative sites that carry
little phylogenetic signal, and the gene metrics and ranking read them as
signal. ``trim_gene_alignments`` runs ClipKIT on each gene alignment before
the stage measures, ranks and concatenates the genes.

ClipKIT is told the sequence type rather than left to guess it: a dna
(codon) gene alignment is trimmed with ``-s nt -co``, which removes whole
codons so the reading frame is kept (ClipKIT 2.1.2 or later), and a protein
one with ``-s aa``.

The genes run in batches, one ``clipkit`` per gene, each batch in one
container when a container backend is active. ClipKIT runs on one thread
unless given ``-t`` (its default, and the flag is not passed because ClipKIT
before 2.14 lacks it), so the parallel batches are the only parallelism. A
gene whose alignment loses every column is dropped and counted.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_chain
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities, preflight
from .base import GeneAlignment
from .select import ALPHABETS

# The --gene-trim choices: no trimming, or ClipKIT.
GENE_TRIMS = ("none", "clipkit")
# ClipKIT -m modes accepted here; smart-gap is the one the CLI uses.
TRIM_MODES = ("smart-gap", "gappy", "kpic-smart-gap", "kpi-smart-gap")
DEFAULT_TRIM_MODE = "smart-gap"

CLIPKIT = ToolCapabilities(
    name="clipkit",
    container="quay.io/biocontainers/clipkit:2.14.0--pyhdfd78af_0",
    conda=("bioconda::clipkit",),
    # `clipkit --version` prints "clipkit <version>"; 2.1.2 added -co (codon).
    required_binaries=(BinarySpec("clipkit", version_args=("--version",), min_version="2.1.2"),),
)

# Genes trimmed per container start.
_BATCH = 50
_UNSAFE = re.compile(r"[^A-Za-z0-9._~+-]")


def check_gene_trim(value: str) -> str:
    """``value`` when it is one of GENE_TRIMS; else a UserInputError."""
    if value not in GENE_TRIMS:
        raise UserInputError(f"--gene-trim must be one of {', '.join(GENE_TRIMS)}, not {value!r}")
    return value


def preflight_clipkit() -> dict[str, str]:
    """Check the clipkit binary (or its image); return the resolved versions."""
    return preflight(CLIPKIT)


def alignment_width(path: Path) -> int:
    """Columns of the first record of the FASTA alignment ``path`` (0 when empty).

    Only the first record is read: the records of an alignment share one
    length, which load_gene checks when the gene is read for the metrics.
    """
    seen = False
    width = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if seen:
                    break
                seen = True
            elif seen:
                width += len(line.strip())
    return width


def trim_gene_alignments(
    gene_alignments: Sequence[GeneAlignment],
    out_dir: Path,
    threads: int,
    logger: logging.Logger,
    mode: str = DEFAULT_TRIM_MODE,
    alphabet: str = "dna",
) -> list[GeneAlignment]:
    """Trim each gene alignment with ClipKIT (``clipkit IN -o OUT -m MODE``,
    then ``-s nt -co`` for a dna alphabet or ``-s aa`` for protein).

    Writes ``out_dir/trimmed/<gene>.fasta`` (characters outside letters,
    digits and ``._~+-`` in the gene id are replaced by ``_``). Records keep
    their names; a trimmed file whose record names differ from its input is
    an error. Returns the trimmed alignments in input order, without the
    genes that lost every column (counted in one log line). A trimmed file
    with no records at all is an error, not an emptied gene.
    """
    if alphabet not in ALPHABETS:
        raise ValueError(f"Unknown alphabet '{alphabet}'; expected dna or protein")
    seq_args = ["-s", "nt", "-co"] if alphabet == "dna" else ["-s", "aa"]
    if mode not in TRIM_MODES:
        raise UserInputError(
            f"Unknown ClipKIT mode '{mode}'. Choose from: {', '.join(TRIM_MODES)}."
        )
    trimmed_dir = out_dir / "trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        (ga, trimmed_dir / f"{name}.fasta")
        for ga, name in zip(
            gene_alignments, _file_names(ga.gene for ga in gene_alignments), strict=True
        )
    ]
    batches = [jobs[i : i + _BATCH] for i in range(0, len(jobs), _BATCH)]

    def run(batch: Sequence[tuple[GeneAlignment, Path]]) -> None:
        run_chain(
            CLIPKIT,
            [
                ("clipkit", ["clipkit", ga.path, "-o", dst, "-m", mode, *seq_args])
                for ga, dst in batch
            ],
            logger=logger,
            cwd=out_dir,
        )

    logger.info(
        "Trimming %d gene alignments with ClipKIT (-m %s %s)", len(jobs), mode, " ".join(seq_args)
    )
    # The first failure cancels the batches not yet started.
    parallel_map(run, batches, max(1, threads), logger=logger)

    out: list[GeneAlignment] = []
    dropped: list[str] = []
    for ga, dst in jobs:
        if not dst.is_file():
            raise WorkdirError(f"clipkit wrote no trimmed alignment for gene {ga.gene}")
        after = _records(dst)
        if not after:
            raise WorkdirError(
                f"clipkit wrote a trimmed alignment without records for gene {ga.gene} ({dst})"
            )
        if all(not seq for _, seq in after):
            dropped.append(ga.gene)
            continue
        before = [name for name, _ in _records(ga.path)]
        if sorted(name for name, _ in after) != sorted(before):
            raise WorkdirError(
                f"clipkit changed the record names of gene {ga.gene} ({dst}); the trimmed "
                "alignment no longer names the genomes of the set."
            )
        out.append(GeneAlignment(ga.gene, dst))
    if dropped:
        logger.info(
            "%d gene(s) were dropped because ClipKIT removed every column: %s",
            len(dropped),
            ", ".join(dropped[:5]) + (" ..." if len(dropped) > 5 else ""),
        )
    return out


def _file_names(genes: Iterable[str]) -> list[str]:
    """A file-name stem per gene id, unique without regard to case."""
    used: set[str] = set()
    names: list[str] = []
    for i, gene in enumerate(genes):
        base = _UNSAFE.sub("_", gene).lstrip(".") or "gene"
        name, n = base, 0
        while name.casefold() in used:
            n += 1
            name = f"{base}.{i}" if n == 1 else f"{base}.{i}.{n}"
        used.add(name.casefold())
        names.append(name)
    return names


def _records(path: Path) -> list[tuple[str, str]]:
    """(first word of the header, sequence) of each record of ``path``."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                words = line[1:].split()
                name, parts = (words[0] if words else ""), []
            else:
                parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(parts)))
    return records
