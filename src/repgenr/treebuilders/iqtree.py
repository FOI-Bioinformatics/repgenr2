"""IQ-TREE tree builder (maximum likelihood from an MSA)."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import atomic_path
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path


class IqtreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="iqtree",
        container="quay.io/biocontainers/iqtree:3.1.3--h8471819_0",
        conda=("bioconda::iqtree",),
        # model=: the substitution model passed as -m (for example GTR+G).
        # Without it IQ-TREE runs ModelFinder, a model-selection pass before
        # the tree search.
        accepted_extras=frozenset({"model"}),
        # 2.0: -B (ultrafast bootstrap) and --threads-max as passed below.
        required_binaries=(BinarySpec("iqtree", version_args=("--version",), min_version="2.0"),),
        recommended_max_genomes=500,
    )
    input_kind = InputKind.MSA_FASTA
    # Without -B the tree carries no support values at all.
    supports_without_bootstrap = False

    # IQ-TREE's ultrafast bootstrap (-B) rejects fewer than 1000 replicates.
    MIN_BOOTSTRAP = 1000

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        msa = as_msa_path(msa_or_genomes)
        if 0 < params.bootstrap < self.MIN_BOOTSTRAP:
            raise UserInputError(
                f"IQ-TREE needs at least {self.MIN_BOOTSTRAP} ultrafast bootstrap "
                f"replicates (-B); got {params.bootstrap}. Use --bootstrap "
                f"{self.MIN_BOOTSTRAP} or more, or 0 for none."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        work_msa = out_dir / "msa.fasta"
        if work_msa.resolve() != msa.resolve():
            shutil.copy2(msa, work_msa)

        cmd: list[str | Path] = [
            "iqtree",
            "-T",
            "auto",
            "--threads-max",
            str(params.threads),
            "-s",
            work_msa,
            # Every alignment this pipeline produces is nucleotide. Left to
            # guess, IQ-TREE refuses a recombination-masked alignment, where
            # a fifth of the characters can be N ("Unknown sequence type").
            "-st",
            "DNA",
            "-redo",  # overwrite checkpoints from a previous run at this path
        ]
        model = str(params.extra.get("model", "")).strip()
        if model:
            cmd += ["-m", model]
        else:
            logger.info(
                "iqtree: no model given, so IQ-TREE runs ModelFinder before the tree "
                "search; --aligner-arg model=GTR+G skips the selection"
            )
        if params.outgroup:
            cmd += ["-o", params.outgroup]
        if params.bootstrap > 0:
            cmd += ["-B", str(params.bootstrap)]
        run_tool(self.capabilities, cmd, logger=logger, log_prefix="iqtree")

        treefile = work_msa.with_suffix(".fasta.treefile")
        if not treefile.exists():
            treefile = Path(str(work_msa) + ".treefile")
        if not treefile.exists():
            raise WorkdirError("IQ-TREE did not produce a .treefile")
        tree = out_dir / "tree.nwk"
        with atomic_path(tree) as tmp:
            shutil.copy2(treefile, tmp)
        return tree
