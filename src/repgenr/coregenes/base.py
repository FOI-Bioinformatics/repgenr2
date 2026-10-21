"""CoreGeneBuilder interface: annotated genomes in, core-gene alignments out.

The builder takes the GFF3 files the ``annotate`` stage wrote (sequences
appended), clusters the genes into a pangenome and aligns the core genes.
It hands over either one aligned FASTA per gene, which the stage ranks,
caps and concatenates itself (``coregenes/select.py``), or its own
concatenation. ``phylo --msa-source coregene`` hands the alignment to
an MSA tree builder; it keeps the invariant sites, so no ascertainment-bias
correction applies.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..core.errors import WorkdirError
from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[CoreGeneBuilder] = Registry("repgenr.coregenes")

CLEAN_MODES = ("strict", "moderate", "sensitive")


@dataclass
class CoreGeneParams:
    threads: int = 16
    # A gene is core when present in at least this fraction of the genomes.
    core_threshold: float = 0.95
    # Panaroo's graph cleaning of annotation errors: strict | moderate | sensitive.
    clean_mode: str = "strict"
    # Take the alignment with badly aligned genes removed (the builder's
    # recommended tree input), else the unfiltered one.
    filtered: bool = True
    # The character type of the per-gene alignments: "dna" (codons) or
    # "protein" (amino acids). A builder that cannot hand over protein
    # alignments lists "alphabet" in ignored_params and the stage refuses it.
    alphabet: str = "dna"
    extra: dict = field(default_factory=dict)


@dataclass(frozen=True)
class GeneAlignment:
    """One core gene's alignment as a builder hands it over."""

    # Family or orthogroup id, unique within one result.
    gene: str
    # Aligned FASTA, nucleotide or protein as CoreGeneParams.alphabet asks;
    # records are named by genome record name. A genome that lacks the gene
    # is absent (no gap-only record is required).
    path: Path


@dataclass
class CoreGeneResult:
    """What a builder returns.

    A builder supplies ``gene_alignments`` (the stage then selects and
    concatenates the genes and records a partition per gene) or its own
    concatenated ``alignment``. Both are allowed; ``alignment`` is then
    ignored by the stage.
    """

    # The tool's own concatenation, used only when gene_alignments is None.
    alignment: Path | None = None
    # The builder's tables (presence/absence, summary) to publish beside it.
    tables: list[Path] = field(default_factory=list)
    core_genes: int | None = None
    # The identity the builder clustered genes at, as it ran (default or set),
    # for the stats file: it governs how many genes cluster at all.
    clustering_identity: float | None = None
    # Builder-specific settings that shaped the result, for the stats file.
    settings: dict = field(default_factory=dict)
    # One aligned FASTA per core gene, or None for a builder that returns
    # only its own concatenation.
    gene_alignments: list[GeneAlignment] | None = None


class CoreGeneBuilder(ABC):
    capabilities: ToolCapabilities
    # Annotator names this builder accepts; empty means any GFF3 with sequences.
    accepted_annotators: tuple[str, ...] = ()

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def call(
        self,
        gffs: Sequence[Path],
        out_dir: Path,
        params: CoreGeneParams,
        logger: logging.Logger,
    ) -> CoreGeneResult:
        """Build the core-gene alignment(s) of ``gffs`` under ``out_dir``."""
        raise NotImplementedError


def annotation_sidecars(gffs: Sequence[Path], tool: str) -> tuple[list[Path], list[Path]]:
    """The ``.faa`` and ``.ffn`` files the annotate stage writes beside each GFF3.

    For builders that read proteins rather than GFF3. Refuses, naming the
    first missing files, when any is absent.
    """
    faas = [g.with_suffix(".faa") for g in gffs]
    ffns = [g.with_suffix(".ffn") for g in gffs]
    missing = [str(p) for p in (*faas, *ffns) if not p.is_file()]
    if missing:
        raise WorkdirError(
            f"{tool} reads the protein and gene files the annotate stage writes beside "
            f"each GFF3, and {len(missing)} are missing ({', '.join(missing[:3])}); "
            "rerun `repgenr annotate --force`."
        )
    return faas, ffns
