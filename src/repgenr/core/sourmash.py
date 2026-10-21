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

import dataclasses
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


def jaccard_to_ani_matrix(sim: npt.NDArray[np.float64], ksize: int) -> npt.NDArray[np.float64]:
    """Elementwise ANI estimate from a k-mer Jaccard similarity matrix.

    Jaccard's containment equivalent is 2j/(1+j); the ANI estimate of the
    k-mer model is its k-th root (the quantity ``sourmash compare --ani``
    reports, within 1e-4 on real genomes). Zero stays zero. Computed here
    rather than with ``--ani`` because sourmash reports an ANI of 0 for
    sketches of small genomes, identical ones included (sourmash 4.9.4,
    10 kb and 30 kb genomes at k=31, scaled=1000). Used by the dereplicator
    and the tree builder.
    """
    import numpy as np

    containment = 2 * sim / (1 + sim)
    return np.where(sim > 0, containment ** (1.0 / ksize), 0.0)
