"""Per-run execution context.

A :class:`WorkdirContext` bundles the working directory, the canonical path
layout, the logger, the loaded config and (lazily) the SQLite manifest. Stages
receive a context plus their parsed parameters; they never re-derive paths or
re-open the logger themselves.
"""

from __future__ import annotations

import logging
import os
from functools import cached_property
from pathlib import Path

from .config import Config
from .contracts import ANNOTATIONS_DIR, PHYLO_MSA_DIR
from .logging import configure_logging
from .manifest import Manifest


class WorkdirContext:
    """Resolved paths and shared services for one RepGenR2 invocation."""

    def __init__(
        self,
        workdir: str | os.PathLike[str],
        *,
        logger: logging.Logger | None = None,
        create: bool = False,
    ):
        self.workdir = Path(workdir).resolve()
        if create:
            self.workdir.mkdir(parents=True, exist_ok=True)
        self.logger = logger or configure_logging(self.workdir)
        self.config = Config.load(self.workdir)
        # Set by the CLI harness from --force: recompute, ignoring caches.
        self.force = False

    # -- canonical path layout (see plan "Canonical data contracts") ----------
    @property
    def genomes_dir(self) -> Path:
        return self.workdir / "genomes"

    @property
    def outgroup_dir(self) -> Path:
        return self.workdir / "outgroup"

    @property
    def derep_dir(self) -> Path:
        return self.workdir / "derep"

    @property
    def representatives_dir(self) -> Path:
        return self.derep_dir / "representatives"

    @property
    def align_dir(self) -> Path:
        return self.workdir / "align"

    @property
    def snp_dir(self) -> Path:
        return self.workdir / "snp"

    @property
    def tree_dir(self) -> Path:
        return self.workdir / "tree"

    @property
    def phylo_msa_dir(self) -> Path:
        """Where phylo's own SNP typing pass writes (snp/ is the snptype stage's)."""
        return self.tree_dir / PHYLO_MSA_DIR

    @property
    def annotations_dir(self) -> Path:
        """The annotate stage's cache: one directory per record and an index."""
        return self.workdir / ANNOTATIONS_DIR

    @property
    def scratch_dir(self) -> Path:
        return self.workdir / "scratch"

    # -- services -------------------------------------------------------------
    def close(self) -> None:
        """Release held resources (the manifest's SQLite connection, if opened)."""
        manifest = self.__dict__.get("manifest")
        if manifest is not None:
            manifest.close()

    @cached_property
    def manifest(self) -> Manifest:
        return Manifest.open(self.workdir)

    def save_config(self) -> None:
        self.config.save(self.workdir)

    def require_dir(self, path: Path, hint: str) -> Path:
        from .errors import WorkdirError

        if not path.exists():
            raise WorkdirError(f"Expected directory not found: {path}\n{hint}")
        return path
