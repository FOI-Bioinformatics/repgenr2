"""Neighbour joining over a distance matrix, with rapidNJ when it is installed.

The pure-Python :func:`repgenr.tree.newick.neighbor_joining` is cubic and
measured at about 1.6 hours for 5000 leaves (docs/audit/scaling-audit.md).
rapidNJ (Simonsen et al. 2008) joins the same matrix in seconds at that size
and runs to tens of thousands of leaves, so a distance-based builder calls
:func:`nj_tree` and gets rapidNJ when the binary is on PATH (or an image is
pinned under a container backend) and the Python code otherwise.

rapidNJ reads a square PHYLIP distance matrix (``n`` on the first line, then
one row per taxon: name and ``n`` distances). It quotes long names itself but
mangles names holding spaces or Newick characters, so the matrix is written
with integer stand-in names and the labels are restored afterwards, quoted as
:func:`repgenr.tree.newick._sanitize` does. ``-n`` adjusts negative branch
lengths to zero, as the Python code does. Command surface checked against
rapidNJ 2.3.2 (bioconda).
"""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import get_config, run_tool
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, tool_available
from .newick import _sanitize, neighbor_joining

NJ_PYTHON = "python"
NJ_RAPIDNJ = "rapidnj"
NJ_AUTO = "auto"
NJ_BACKENDS = (NJ_AUTO, NJ_RAPIDNJ, NJ_PYTHON)

# Above this the Python joining takes hours; refused unless rapidNJ runs.
PYTHON_NJ_MAX = 5000
# Declared, not measured here: rapidNJ's authors report trees of this size.
RAPIDNJ_MAX = 50000

RAPIDNJ = ToolCapabilities(
    name="rapidnj",
    container="quay.io/biocontainers/rapidnj:2.3.3--hd63eeec_0",
    conda=("bioconda::rapidnj",),
    # rapidnj has no version flag; its help exits non-zero without an input.
    required_binaries=(BinarySpec("rapidnj", version_args=()),),
)

_LEAF = re.compile(r"(?<![\w'])(?:'?)(\d+)(?:'?)(?=:)")


def choose_backend(requested: str, n: int, logger: logging.Logger) -> str:
    """The backend to join ``n`` leaves with: ``requested`` when it is explicit,
    else rapidNJ when it is available, else Python (refused above its ceiling)."""
    if requested not in NJ_BACKENDS:
        raise UserInputError(f"nj must be one of {', '.join(NJ_BACKENDS)}, not {requested!r}")
    if requested == NJ_RAPIDNJ:
        return NJ_RAPIDNJ
    if requested == NJ_AUTO and tool_available(RAPIDNJ):
        return NJ_RAPIDNJ
    if n > PYTHON_NJ_MAX:
        raise WorkdirError(
            f"Neighbour joining of {n} leaves in Python extrapolates to many hours "
            f"(measured fit: about 1.6 h at 5000). Install rapidnj (conda: bioconda::rapidnj) "
            "so the fast joining is used, or dereplicate further first."
        )
    if requested == NJ_AUTO:
        logger.info(
            "rapidnj is not installed; joining %d leaves in Python (seconds below a few "
            "hundred leaves, minutes near a thousand, hours near %d)",
            n,
            PYTHON_NJ_MAX,
        )
    return NJ_PYTHON


def nj_tree(
    labels: Sequence[str],
    dist: Sequence[Sequence[float]],
    out_dir: Path,
    logger: logging.Logger,
    *,
    backend: str = NJ_AUTO,
) -> str:
    """An unrooted Newick tree joining ``labels`` by the symmetric ``dist`` matrix."""
    n = len(labels)
    if n != len(dist) or any(len(row) != n for row in dist):
        raise ValueError("the distance matrix must be square and match the labels")
    chosen = choose_backend(backend, n, logger)
    if chosen == NJ_PYTHON or n < 3:
        return neighbor_joining(list(labels), [list(row) for row in dist])
    return _rapidnj(labels, dist, out_dir, logger)


def write_phylip(path: Path, n: int, dist: Sequence[Sequence[float]]) -> None:
    """A square PHYLIP matrix with stand-in names ``0`` to ``n - 1``."""
    with open(path, "w", encoding="utf-8") as fo:
        fo.write(f"{n}\n")
        for i in range(n):
            fo.write(str(i) + " " + " ".join(f"{max(float(x), 0.0):.6f}" for x in dist[i]) + "\n")


def restore_labels(newick: str, labels: Sequence[str]) -> str:
    """Replace the stand-in names rapidNJ wrote (bare or quoted, each before a
    branch length) with the real labels, quoted where Newick needs it."""

    def repl(m: re.Match[str]) -> str:
        return _sanitize(labels[int(m.group(1))])

    return _LEAF.sub(repl, newick)


def _rapidnj(
    labels: Sequence[str], dist: Sequence[Sequence[float]], out_dir: Path, logger: logging.Logger
) -> str:
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix = out_dir / "distances.phylip"
    write_phylip(matrix, len(labels), dist)
    tree = out_dir / "rapidnj.nwk"
    run_tool(
        RAPIDNJ,
        ["rapidnj", matrix, "-i", "pd", "-n", "-x", tree],
        logger=logger,
        log_prefix="rapidnj",
    )
    if not tree.is_file() or tree.stat().st_size == 0:
        raise WorkdirError("rapidnj wrote no tree")
    text = tree.read_text(encoding="utf-8").strip()
    if not text.endswith(";"):
        raise WorkdirError(f"rapidnj wrote an incomplete tree: {text[-40:]!r}")
    logger.info(
        "Joined %d leaves with rapidnj%s",
        len(labels),
        " (container)" if get_config().active else "",
    )
    return restore_labels(text, labels)
