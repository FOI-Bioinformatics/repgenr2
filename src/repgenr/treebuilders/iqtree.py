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
from ..tree.alignment import add_counts, format_counts, strip_invariant_columns
from .base import (
    InputKind,
    TreeBuilder,
    TreeParams,
    as_msa_path,
    partition_names,
    refuse_partitions_with_asc,
)


class IqtreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="iqtree",
        container="quay.io/biocontainers/iqtree:3.1.3--h8471819_0",
        conda=("bioconda::iqtree",),
        # model=: the substitution model passed as -m (for example GTR+G).
        # Without it IQ-TREE runs ModelFinder, a model-selection pass before
        # the tree search. asc=lewis: the +ASC correction instead of the
        # constant-site counts (-fconst) the typing pass provides.
        accepted_extras=frozenset({"model", "asc"}),
        # 2.0: -B (ultrafast bootstrap) and --threads-max as passed below.
        required_binaries=(BinarySpec("iqtree", version_args=("--version",), min_version="2.0"),),
        recommended_max_genomes=500,
    )
    input_kind = InputKind.MSA_FASTA
    # Without -B the tree carries no support values at all.
    supports_without_bootstrap = False
    # -p: one partition per gene, edge-linked with proportional branch lengths.
    supports_partitions = True

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
        refuse_partitions_with_asc(params, "iqtree")
        out_dir.mkdir(parents=True, exist_ok=True)
        work_msa = out_dir / "msa.fasta"
        asc = str(params.extra.get("asc", "")).strip().lower()
        if asc not in ("", "lewis"):
            raise UserInputError(f"asc={asc!r} is not supported; use asc=lewis or leave it unset.")
        model = str(params.extra.get("model", "")).strip()
        model_has_asc = "+asc" in model.lower()
        if model_has_asc:
            # The user's model carries its own correction: no counts, no
            # second +ASC, but the alignment must still hold no invariant column.
            asc = "model"
            logger.info("iqtree: the model %s carries its own ascertainment-bias term", model)
        constant = params.constant_sites
        if asc or constant is not None:
            # Every correction wants no invariant column in the alignment; the
            # columns removed here are constant sites and join the counts.
            removed = strip_invariant_columns(msa, work_msa, logger)
            if constant is not None:
                constant = add_counts(constant, removed)
        elif work_msa.resolve() != msa.resolve():
            shutil.copy2(msa, work_msa)

        partition_args: list[str | Path] = []
        if params.partitions is not None:
            partition_file = _write_partition_file(params, out_dir)
            # Given -p, IQ-TREE names its outputs after the partition file;
            # -pre keeps them where an unpartitioned run writes them.
            partition_args = ["-p", partition_file, "-pre", work_msa]
            logger.info(
                "iqtree: edge-linked partition model over %d partitions (-p %s)",
                len(params.partitions),
                partition_file.name,
            )
        cmd: list[str | Path] = [
            "iqtree",
            "-T",
            "auto",
            "--threads-max",
            str(params.threads),
            "-s",
            work_msa,
            *partition_args,
            # Every alignment this pipeline produces is nucleotide. Left to
            # guess, IQ-TREE refuses a recombination-masked alignment, where
            # a fifth of the characters can be N ("Unknown sequence type").
            "-st",
            "DNA",
            "-redo",  # overwrite checkpoints from a previous run at this path
        ]
        if asc == "lewis":
            # +ASC needs a model term; MFP keeps ModelFinder when none is given.
            cmd += ["-m", f"{model or 'MFP'}+ASC"]
            logger.info("iqtree: ascertainment-bias correction +ASC (asc=lewis)")
        elif model:
            cmd += ["-m", model]
        else:
            logger.info(
                "iqtree: no model given, so IQ-TREE runs ModelFinder before the tree "
                "search; model=GTR+G (--aligner-arg) skips the selection"
            )
        if not asc and constant is not None:
            cmd += ["-fconst", format_counts(constant, ",")]
            logger.info("iqtree: constant sites added as -fconst %s", format_counts(constant, ","))
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


def _write_partition_file(params: TreeParams, out_dir: Path) -> Path:
    """Write ``out_dir/partitions.txt`` in the RAxML-style form IQ-TREE reads."""
    assert params.partitions is not None
    kind = "AA" if params.alphabet == "protein" else "DNA"
    names = partition_names([gene for gene, _, _ in params.partitions])
    lines = [
        f"{kind}, {name} = {start}-{end}"
        for name, (_, start, end) in zip(names, params.partitions, strict=True)
    ]
    path = out_dir / "partitions.txt"
    with atomic_path(path) as tmp:
        tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
