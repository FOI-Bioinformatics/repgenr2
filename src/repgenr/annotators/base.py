"""Annotator interface: one genome in, a GFF3 with its sequences out.

The core-gene route needs every genome annotated once (bakta takes minutes
per genome), so the ``annotate`` stage owns an ``annotations/`` cache and
calls an adapter per genome. Adapters write ``<name>.gff3`` with the FASTA
section appended (Panaroo reads the sequences from it), ``<name>.faa`` and
``<name>.ffn``.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

from ..core.errors import WorkdirError
from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Annotator] = Registry("repgenr.annotators")


@dataclass
class AnnotateParams:
    threads: int = 16
    # The annotator's database directory, for tools that need one (bakta).
    db: Path | None = None
    extra: dict = field(default_factory=dict)


@dataclass
class AnnotateResult:
    gff: Path
    faa: Path | None = None
    ffn: Path | None = None


class Annotator(ABC):
    capabilities: ToolCapabilities
    # True for a tool that needs a database directory (``AnnotateParams.db``).
    requires_database: bool = False

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    def database_version(self, db: Path) -> str:
        """The version string of the database at ``db`` (empty when none)."""
        return ""

    @abstractmethod
    def call(
        self,
        genome: Path,
        record_name: str,
        out_dir: Path,
        params: AnnotateParams,
        logger: logging.Logger,
    ) -> AnnotateResult:
        """Annotate ``genome`` into ``out_dir`` under the prefix ``record_name``."""
        raise NotImplementedError


def require_fasta_section(gff: Path, tool: str) -> None:
    """Panaroo reads the gene sequences from the ``##FASTA`` section of the GFF3."""
    with open(gff, encoding="utf-8", errors="replace") as fo:
        for line in fo:
            if line.startswith("##FASTA"):
                return
    raise WorkdirError(
        f"{tool} wrote {gff.name} without a ##FASTA section; the core-gene step needs the "
        "sequences in the GFF3"
    )
