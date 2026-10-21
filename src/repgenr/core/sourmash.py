"""The one sourmash tool specification shared by every sourmash caller.

The dereplicator, the tree builder, the classifier and the genome sketches
(:mod:`repgenr.core.sketches`) all run the same ``sourmash`` binary, from the
host or from one pinned BioContainer. The image, the conda spec and the binary
check are declared here once, so the container pin cannot drift between them;
each adapter derives its own :class:`ToolCapabilities` with
:func:`sourmash_capabilities`, adding only its parameters.

This module lives in ``core`` and imports no adapter, so every adapter family
can import it without an import cycle.
"""

from __future__ import annotations

import csv
import dataclasses
import logging
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .binaries import BinarySpec
from .plugins import ToolCapabilities

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt

SOURMASH_IMAGE = "quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0"
SOURMASH_CONDA: tuple[str, ...] = ("bioconda::sourmash",)
SOURMASH_BINARY = BinarySpec("sourmash", version_args=("--version",), min_version="4.0")

# ``sourmash sketch`` over several k-mer sizes starts a Rayon thread pool.
# Measured on macOS (sourmash 4.9.4, 400000 reads, three k): 14 s and 82 s
# of system time with the default pool against 7 s with one thread, the
# signatures identical; the pool also competes with the other sketches of a
# batch and with an assembler running beside it. The sketch is one thread.
SKETCH_TOOL_ENV = {"RAYON_NUM_THREADS": "1"}

# The base specification. ``name`` is also the key of --bin-dir and of the
# per-tool container settings, so every caller shares it.
SOURMASH_TOOL = ToolCapabilities(
    name="sourmash",
    required_binaries=(SOURMASH_BINARY,),
    container=SOURMASH_IMAGE,
    conda=SOURMASH_CONDA,
)


def sourmash_capabilities(**overrides: Any) -> ToolCapabilities:
    """The shared sourmash specification with adapter-specific fields set.

    The image, conda spec and binary check cannot be overridden here; pass
    only parameter fields (``default_params``, ``accepted_extras``,
    ``ignored_params``, ``reads_gzip``, scale hints).
    """
    fixed = {"name", "required_binaries", "container", "conda"}
    clash = sorted(fixed & overrides.keys())
    if clash:
        raise ValueError(f"sourmash_capabilities: {', '.join(clash)} come from the shared spec")
    return dataclasses.replace(SOURMASH_TOOL, **overrides)


# A pair whose sketches share fewer hashes than this has an unknown sketch
# ANI (NaN). Measured on 24 Wolbachia genomes of 1.0 to 1.5 Mb at
# scaled=1000 (docs/audit/phylo-methods-2026-10.md, "Sketch ANI audit
# (2026-10-21)"): the 57 pairs skani cannot align share 1 to 22 hashes at
# k=31 and 13 to 33 at k=21, while the lowest pair skani aligns (79.4 %)
# shares 8 (k=31) and 26 (k=21). Below 30 hashes the Poisson sampling error
# alone is near 20 % of the count, and the shared hashes come mostly from
# conserved genes and mobile elements rather than from the whole genome.
MIN_SHARED_HASHES = 30


def jaccard_to_ani_matrix(
    sim: npt.NDArray[np.float64],
    ksize: int,
    *,
    n_hashes: Sequence[float] | None = None,
    min_shared: int = MIN_SHARED_HASHES,
) -> tuple[npt.NDArray[np.float64], int]:
    """Elementwise ANI estimate from a k-mer Jaccard similarity matrix, and
    the number of pairs whose ANI is unknown.

    Jaccard's containment equivalent is 2j/(1+j); the ANI estimate of the
    k-mer model is its k-th root (the quantity ``sourmash compare --ani``
    reports; identical to it on 276 Wolbachia pairs at k=21 and k=31).
    Computed here rather than with ``--ani`` because sourmash reports an ANI
    of 0 for sketches of small genomes, identical ones included (sourmash
    4.9.4, 10 kb and 30 kb genomes at k=31, scaled=1000). Used by the
    dereplicator and the tree builder.

    With the hash count of each signature (``n_hashes``, in the order of the
    matrix rows, at the scaled value the matrix was computed at), a pair
    whose number of shared hashes, ``j * (a + b) / (1 + j)`` for sketches of
    ``a`` and ``b`` hashes, is below ``min_shared`` gets an unknown ANI
    (NaN): too few hashes to support an estimate. On distant pairs a few
    hashes from conserved genes would otherwise give an ANI of 0.80 to 0.88
    at k=31 (measured on Wolbachia supergroups that skani does not align).
    A pair without a shared hash is unknown as well; the diagonal stays 1.0.
    Without sizes, zero stays zero and no pair is unknown.
    """
    import numpy as np

    sim = np.asarray(sim, dtype=float)
    containment = 2 * sim / (1 + sim)
    ani = np.where(sim > 0, containment ** (1.0 / ksize), 0.0)
    if n_hashes is None:
        return ani, 0
    sizes = np.asarray(n_hashes, dtype=float)
    if ani.ndim != 2 or ani.shape[0] != ani.shape[1] or sizes.shape != (ani.shape[0],):
        raise ValueError(f"{sizes.size} signature sizes for a matrix of shape {ani.shape}")
    # The intersection size from the Jaccard index and the two sketch sizes;
    # rounded, since it is a count recovered through floating point.
    shared = np.round(sim * (sizes[:, None] + sizes[None, :]) / (1 + sim))
    unknown = shared < min_shared
    np.fill_diagonal(unknown, False)
    ani = np.where(unknown, np.nan, ani)
    np.fill_diagonal(ani, 1.0)
    n_unknown = int(np.count_nonzero(np.triu(unknown, k=1)))
    return ani, n_unknown


def signature_sizes(
    caps: ToolCapabilities,
    fofn: Path,
    out_csv: Path,
    ksize: int,
    labels: Sequence[str],
    logger: logging.Logger,
    run: Callable[..., Any],
    *,
    extra_mounts: Sequence[str] = (),
) -> tuple[list[float], int] | None:
    """Hash count of each signature named in ``labels`` at ``ksize``, and the
    scaled value they are expressed at, from ``sourmash sig describe --csv``.

    ``labels`` are the column labels of ``sourmash compare`` (a signature's
    name, else its file name, else the first 8 characters of its md5, as
    sourmash displays it). Signatures at different scaled values are brought
    to the largest one, which ``compare`` downsamples to. Returns None, after
    one WARNING, when the describe call fails or a label has no row or more
    than one; the caller then converts without the shared-hash minimum. ``run`` is
    the caller's ``run_tool``.
    """
    from .errors import MissingBinaryError, ToolExecutionError

    fallback = "the sketch ANI is computed without the shared-hash minimum"
    out_csv.unlink(missing_ok=True)  # never read a table from an earlier run
    try:
        run(
            caps,
            [
                "sourmash",
                "sig",
                "describe",
                "-q",
                "-k",
                str(ksize),
                "--csv",
                out_csv,
                "--from-file",
                fofn,
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=list(extra_mounts),
        )
        with open(out_csv, encoding="utf-8", newline="") as fo:
            rows = list(csv.DictReader(fo))
    except (ToolExecutionError, MissingBinaryError, OSError) as exc:
        logger.warning("sourmash sig describe failed (%s); %s", exc, fallback)
        return None
    by_label: dict[str, set[tuple[float, int]]] = {}
    try:
        for row in rows:
            if int(row["ksize"]) != ksize:
                continue
            label = row["name"] or row["filename"] or row["md5"][:8]
            by_label.setdefault(label, set()).add((float(row["n_hashes"]), int(row["scaled"])))
    except (KeyError, ValueError) as exc:
        logger.warning("sourmash sig describe wrote an unexpected table (%s); %s", exc, fallback)
        return None
    unmatched = [label for label in labels if len(by_label.get(label, ())) != 1]
    if unmatched:
        logger.warning(
            "sourmash sig describe gives no single size for %d of %d signatures (%s); %s",
            len(unmatched),
            len(labels),
            unmatched[0],
            fallback,
        )
        return None
    pairs = [next(iter(by_label[label])) for label in labels]
    scaled = max(s for _n, s in pairs)
    if scaled <= 1:
        logger.warning("sourmash signatures have no scaled value above 1; %s", fallback)
        return None
    return [n * s / scaled for n, s in pairs], scaled
