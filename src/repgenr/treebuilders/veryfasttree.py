"""VeryFastTree tree builder (FastTree-2 reimplemented for many cores).

Same command family as FastTree: ``-nt -gtr`` on a nucleotide alignment,
the tree on stdout, SH-like local supports from ``-boot N`` resamples. It
adds ``-threads N`` (FastTreeMP reads OMP_NUM_THREADS only) and
``-double-precision``; its authors report trees of a million sequences
(Pineiro et al. 2020, Bioinformatics 36:4658, doi:10.1093/bioinformatics/btaa582;
version 4 in Pineiro and Pichel 2024, GigaScience 13, giae055).
Command surface checked against VeryFastTree 4.0.5 (bioconda).
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path


class VeryFastTreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="veryfasttree",
        container="quay.io/biocontainers/veryfasttree:4.0.5--h9948957_0",
        conda=("bioconda::veryfasttree",),
        # -h prints "VeryFastTree 4.0.5 (OpenMP)" and the usage, exit 0; there
        # is no --version flag.
        required_binaries=(BinarySpec("VeryFastTree", version_args=("-h",), min_version="4.0"),),
        # Declared, not measured here: the authors report far larger inputs;
        # this is the point past which the alignment itself, not the builder,
        # is the practical limit.
        recommended_max_genomes=20000,
        accepted_extras=frozenset({"gamma", "double_precision"}),
    )
    input_kind = InputKind.MSA_FASTA

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        msa = as_msa_path(msa_or_genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        if params.constant_sites is not None:
            logger.info(
                "veryfasttree applies no ascertainment-bias correction: the branch lengths of "
                "this variable-site alignment are inflated; iqtree or raxmlng use the "
                "constant-site counts"
            )
        cmd: list = ["VeryFastTree", "-nt", "-gtr", "-threads", str(params.threads)]
        if parse_extra_bool(params.extra, "gamma", False):
            cmd.append("-gamma")
        if parse_extra_bool(params.extra, "double_precision", False):
            cmd.append("-double-precision")
        if params.bootstrap > 0:
            cmd += ["-boot", str(params.bootstrap)]
        env = {**os.environ, "OMP_NUM_THREADS": str(params.threads)}
        run_tool(
            self.capabilities,
            [*cmd, msa],
            logger=logger,
            log_prefix="veryfasttree",
            stdout_path=tree,
            env=env,
        )
        return tree
