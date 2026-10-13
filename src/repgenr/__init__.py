"""RepGenR2: Representative-Genome Repositories.

Modular pipeline for genome selection, dereplication, whole-genome alignment,
SNP typing and phylogenetics. See the per-stage subcommands under
``repgenr.cli`` and the pluggable tool families under ``repgenr.dereplicators``,
``repgenr.aligners``, ``repgenr.snptypers`` and ``repgenr.treebuilders``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("repgenr2")
except PackageNotFoundError:  # running from a source tree without install
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
