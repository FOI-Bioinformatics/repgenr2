"""RepGenR2 command-line interface (entry point).

Replaces the old ``repgenr.py`` string-rewriting dispatcher with a real Typer
app. The app, top-level callback and shared stage harness live in
:mod:`repgenr.cli.base`; the commands are grouped by domain across the
``cmd_*`` modules. Importing those modules here registers their commands on the
single shared ``app``, which the ``repgenr`` console script targets.
"""

from __future__ import annotations

import sys

# Importing each command module registers its @app.command() subcommands on the
# shared app. Imported for side effects only.
from . import (  # noqa: F401  (registration side effects)
    cmd_bacterial,
    cmd_ingest,
    cmd_misc,
    cmd_phylo,
    cmd_reads,
    cmd_run,
    cmd_steps,
    cmd_viral,
)
from .base import app

__all__ = ["app"]


if __name__ == "__main__":  # pragma: no cover
    sys.exit(app())
