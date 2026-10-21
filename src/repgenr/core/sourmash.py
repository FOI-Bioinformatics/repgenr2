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


# The probability above which sourmash treats the shared hashes of a pair as
# compatible with chance (``jaccard_ani``'s ``prob_threshold``).
P_NOTHING_IN_COMMON_THRESHOLD = 1e-3


def exp_n_mutated(n_kmers: Any, ksize: int, mutation_rate: Any) -> Any:
    """Expected number of mutated k-mers among ``n_kmers`` at a per-base
    ``mutation_rate``: ``L * (1 - (1 - r)^k)``.

    Ported from sourmash 4.9.4 (``sourmash/distance_utils.py``,
    ``exp_n_mutated`` and ``r1_to_q``); works on numpy arrays.
    """
    return n_kmers * (1.0 - (1.0 - mutation_rate) ** ksize)


def p_nothing_in_common(ani: Any, ksize: int, scaled: int, n_unique_kmers: Any) -> Any:
    """Expected probability that two FracMinHash sketches share no hash by
    chance alone, for sequences of ``n_unique_kmers`` k-mers at ``ani``.

    ``exp((L - exp_n_mutated(L, k, 1 - ani)) * ln(1 - 1/scaled))``; 1.0 at
    ANI 0 and 0.0 at ANI 1. Ported from sourmash 4.9.4
    (``sourmash/distance_utils.py``, ``get_exp_probability_nothing_common``
    and ``get_expected_log_probability``); works on numpy arrays.
    """
    import numpy as np

    ani_a = np.asarray(ani, dtype=float)
    n_a = np.asarray(n_unique_kmers, dtype=float)
    log_p = (n_a - exp_n_mutated(n_a, ksize, 1.0 - ani_a)) * np.log(1.0 - 1.0 / float(scaled))
    out = np.where(ani_a <= 0.0, 1.0, np.where(ani_a >= 1.0, 0.0, np.exp(log_p)))
    return float(out) if out.ndim == 0 else out


def jaccard_to_ani_matrix(
    sim: npt.NDArray[np.float64],
    ksize: int,
    *,
    scaled: int | None = None,
    n_hashes: Sequence[float] | None = None,
    p_threshold: float = P_NOTHING_IN_COMMON_THRESHOLD,
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

    With ``scaled`` and the hash count of each signature (``n_hashes``, in
    the order of the matrix rows), the chance check of sourmash's
    ``jaccard_ani`` is applied: the number of k-mers is the mean of the two
    hash counts times ``scaled``, and where :func:`p_nothing_in_common` at
    the estimate exceeds ``p_threshold`` the shared hashes are compatible
    with chance and the ANI is set to NaN (unknown). On distant pairs a few
    hashes from conserved genes would otherwise give an ANI of 0.80 to 0.87
    (measured on Wolbachia supergroups that skani does not align). A pair
    without a shared hash is unknown as well; the diagonal stays 1.0.
    Without sizes, zero stays zero and no pair is unknown.
    """
    import numpy as np

    sim = np.asarray(sim, dtype=float)
    containment = 2 * sim / (1 + sim)
    ani = np.where(sim > 0, containment ** (1.0 / ksize), 0.0)
    if scaled is None or n_hashes is None:
        return ani, 0
    sizes = np.asarray(n_hashes, dtype=float)
    if ani.ndim != 2 or ani.shape[0] != ani.shape[1] or sizes.shape != (ani.shape[0],):
        raise ValueError(f"{sizes.size} signature sizes for a matrix of shape {ani.shape}")
    # sourmash rounds the mean k-mer count of the pair (jaccard_ani).
    n_kmers = np.round((sizes[:, None] + sizes[None, :]) / 2.0 * scaled)
    unknown = p_nothing_in_common(ani, ksize, scaled, n_kmers) > p_threshold
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
    than one; the caller then converts without the chance check. ``run`` is
    the caller's ``run_tool``.
    """
    from .errors import MissingBinaryError, ToolExecutionError

    fallback = "the sketch ANI is computed without the check for chance-level sharing"
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
