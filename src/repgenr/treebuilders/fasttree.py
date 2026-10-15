"""FastTree tree builder (approximate ML from an MSA; fast at scale)."""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path


class FasttreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="fasttree",
        container="quay.io/biocontainers/fasttree:2.2.0--h7b50bb2_1",
        conda=("bioconda::fasttree",),
        # -help prints the usage with the version line; no arguments would
        # make FastTree read an alignment from stdin.
        required_binaries=(BinarySpec("FastTree", version_args=("-help",), min_version="2.1"),),
        recommended_max_genomes=5000,
        accepted_extras=frozenset({"gamma"}),
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
        # Only consult the host PATH when running natively; inside a container
        # the image's canonical lowercase name applies regardless of host state.
        from ..core.bindirs import host_which
        from ..core.containers import get_config

        if get_config().active:
            binary = "fasttree"
        else:
            binary = "FastTree" if host_which("FastTree", self.capabilities.name) else "fasttree"
        if params.constant_sites is not None:
            logger.info(
                "fasttree applies no ascertainment-bias correction: the branch lengths of "
                "this variable-site alignment are inflated; iqtree or raxmlng use the "
                "constant-site counts"
            )
        cmd: list = [binary, "-nt", "-gtr"]
        if parse_extra_bool(params.extra, "gamma", False):
            # Rescale the branch lengths under the Gamma20 likelihood after
            # the CAT search; FastTree's own recommendation for comparable lengths.
            cmd.append("-gamma")
        if params.bootstrap > 0:
            # FastTree's SH-like local supports are computed from `-boot N`
            # resamples (default 1000); --bootstrap N sets N instead of being
            # dropped on the floor.
            cmd += ["-boot", str(params.bootstrap)]
        # FastTreeMP takes its thread count from OpenMP, not from a flag.
        env = {**os.environ, "OMP_NUM_THREADS": str(params.threads)}
        run_tool(
            self.capabilities,
            [*cmd, msa],
            logger=logger,
            log_prefix="fasttree",
            stdout_path=tree,
            env=env,
        )
        return tree
