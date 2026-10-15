"""TreeBuilder interface (tree inference, decoupled from alignment).

MSA-based builders (iqtree, FastTree, RAxML-NG) consume an MSA-FASTA from an
aligner or a SNP typer's core-SNP alignment. Alignment-free builders (mashtree,
sourmash) consume genome files directly and need no MSA source.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[TreeBuilder] = Registry("repgenr.treebuilders")


class InputKind(Enum):
    MSA_FASTA = "msa_fasta"  # needs an MSA source (aligner or snptyper)
    GENOMES = "genomes"  # alignment-free; consumes genome files directly


def as_msa_path(value: Path | Sequence[Path]) -> Path:
    """Narrow a build() argument to a single MSA-FASTA path."""
    if isinstance(value, Path):
        return value
    if isinstance(value, str):
        return Path(value)
    raise TypeError("MSA tree builders expect a single MSA-FASTA path")


def as_genome_list(value: Path | Sequence[Path]) -> list[Path]:
    """Narrow a build() argument to a list of genome paths."""
    if isinstance(value, (str, Path)):
        return [Path(value)]
    return [Path(p) for p in value]


@dataclass
class TreeParams:
    threads: int = 16
    outgroup: str | None = None
    bootstrap: int = 0
    extra: dict = field(default_factory=dict)
    # Genome path -> sourmash signature file of the workdir sketches, given
    # only to a builder whose sketch_request() asks for parameters they hold.
    sketches: Mapping[Path, Path] | None = None
    # Constant A, C, G, T columns the variable-site alignment lost, for the
    # builder's ascertainment-bias correction; None for a whole-genome
    # alignment or when nothing counted them.
    constant_sites: tuple[int, int, int, int] | None = None


class TreeBuilder(ABC):
    capabilities: ToolCapabilities
    input_kind: InputKind = InputKind.MSA_FASTA
    # False for a builder whose tree has no support values unless a bootstrap
    # is requested (IQ-TREE). FastTree writes SH-like local supports anyway and
    # RAxML-NG bootstraps by default, so the stage warns only where -B matters.
    supports_without_bootstrap: bool = True

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    def sketch_request(self, extra: Mapping[str, object]) -> tuple[int, int] | None:
        """The (ksize, scaled) of the sourmash sketches this builder compares at.

        None (the default) means the builder reads no genome sketches; one
        that returns a pair receives ``TreeParams.sketches``.
        """
        return None

    def distance_matrix(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        """Optional capability: pairwise distance matrix over ``genomes``.

        Used by the viral outgroup selection. Builders that can emit a matrix
        override this; the default signals the capability is absent.
        """
        raise NotImplementedError(
            f"Tree builder '{self.capabilities.name}' does not produce a "
            "distance matrix (no distance_matrix() implementation)."
        )

    @abstractmethod
    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        """Build a tree; return the path to a Newick file."""
        raise NotImplementedError
