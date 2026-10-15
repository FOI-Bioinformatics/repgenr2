"""SNP typer interface.

A SNP typer takes genomes (and usually a reference) and produces a core-SNP
alignment plus optional VCF and SNP distance matrix. The core-SNP FASTA is both
a standalone typing deliverable and an alternative MSA source for the phylo
stage (interchangeable with an aligner's MSA).
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[SnpTyper] = Registry("repgenr.snptypers")


@dataclass
class SnpParams:
    threads: int = 16
    reference: Path | None = None
    extra: dict = field(default_factory=dict)


@dataclass
class SnpResult:
    core_snp_fasta: Path
    vcf: Path | None = None
    snp_distance_matrix: Path | None = None
    masked: bool = False
    # Whole-genome alignment in reference coordinates (one record per genome,
    # equal lengths). Required input for recombination maskers; None when the
    # typer only produces variable sites.
    full_alignment: Path | None = None
    # Constant A, C, G, T columns of the whole-genome alignment, for the tree
    # builder's ascertainment-bias correction; None when there is no such
    # alignment to count them from.
    constant_sites: tuple[int, int, int, int] | None = None


class SnpTyper(ABC):
    capabilities: ToolCapabilities
    requires_reference: bool = True
    # Whether ``call`` returns a whole-genome alignment (``SnpResult.full_alignment``),
    # which is what a masker needs. Declared up front so ``--mask`` with a
    # typer that cannot feed it is refused before any SNP calling runs.
    produces_full_alignment: bool = True

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        raise NotImplementedError
