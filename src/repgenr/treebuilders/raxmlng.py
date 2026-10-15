"""RAxML-NG tree builder (maximum likelihood from an MSA)."""

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
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path


class RaxmlNgBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="raxmlng",
        container="quay.io/biocontainers/raxml-ng:2.0.3--h870a6a7_0",
        conda=("bioconda::raxml-ng",),
        # model=: the substitution model (GTR+G). asc=lewis: +ASC_LEWIS instead
        # of the per-base +ASC_STAM counts the typing pass provides.
        accepted_extras=frozenset({"model", "asc"}),
        required_binaries=(BinarySpec("raxml-ng", version_args=("--version",), min_version="1.0"),),
        recommended_max_genomes=1000,
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
        prefix = out_dir / "raxml"
        model = str(params.extra.get("model", "GTR+G")).strip() or "GTR+G"
        asc = str(params.extra.get("asc", "")).strip().lower()
        if asc not in ("", "lewis"):
            raise UserInputError(f"asc={asc!r} is not supported; use asc=lewis or leave it unset.")
        constant = params.constant_sites
        if asc == "lewis" or constant is not None:
            # RAxML-NG refuses an +ASC model over an alignment that still holds
            # an invariant column; the columns removed are constant sites.
            variable_msa = out_dir / "msa.variable.fasta"
            removed = strip_invariant_columns(msa, variable_msa, logger)
            msa = variable_msa
            if asc == "lewis":
                model += "+ASC_LEWIS"
                logger.info("raxml-ng: ascertainment-bias correction +ASC_LEWIS (asc=lewis)")
            else:
                assert constant is not None
                constant = add_counts(constant, removed)
                model += "+ASC_STAM{" + format_counts(constant, "/") + "}"
                logger.info("raxml-ng: constant sites given as %s", model.split("+ASC_")[1])
        cmd: list[str | Path] = [
            "raxml-ng",
            "--all",
            "--msa",
            msa,
            "--model",
            model,
            # auto{N}: let RAxML-NG pick an efficient thread count up to the
            # budget, avoiding its core-oversubscription guard on small alignments.
            "--threads",
            f"auto{{{params.threads}}}",
            "--prefix",
            prefix,
            "--redo",  # overwrite any outputs from a previous run at this prefix
        ]
        # Bound the bootstrap. With `--all` and no `--bs-trees`, RAxML-NG defaults
        # to autoMRE{1000}; on large or divergent alignments that may never reach
        # the MRE convergence criterion and runs for hours past the (already
        # computed) ML tree. Cap the adaptive default at 200 replicates; an
        # explicit `bootstrap` still wins.
        if params.bootstrap > 0:
            cmd += ["--bs-trees", str(params.bootstrap)]
        else:
            cmd += ["--bs-trees", "autoMRE{200}"]
        if params.outgroup:
            cmd += ["--outgroup", params.outgroup]
        run_tool(self.capabilities, cmd, logger=logger, log_prefix="raxml-ng")

        best = Path(str(prefix) + ".raxml.bestTree")
        if not best.exists():
            raise WorkdirError("RAxML-NG did not produce a bestTree")
        # `--all` also writes the best tree annotated with bootstrap support;
        # publish that one so downstream support-based collapsing has values.
        # The bare best tree is the fallback when the support step did not run.
        support = Path(str(prefix) + ".raxml.support")
        source = support if support.exists() else best
        tree = out_dir / "tree.nwk"
        with atomic_path(tree) as tmp:
            shutil.copy2(source, tmp)
        return tree
