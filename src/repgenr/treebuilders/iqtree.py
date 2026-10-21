"""IQ-TREE tree builder (maximum likelihood from an MSA)."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Mapping, Sequence
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
    refuse_asc_on_protein,
    refuse_partitions_with_asc,
)


class IqtreeBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="iqtree",
        container="quay.io/biocontainers/iqtree:3.1.3--h8471819_0",
        conda=("bioconda::iqtree",),
        # model=: the substitution model passed as -m (for example GTR+G).
        # Without it IQ-TREE runs ModelFinder, a model-selection pass before
        # the tree search; with partitions GTR+G (LG+G for protein) is applied
        # to every partition instead, and model=MFP selects one per partition.
        # asc=lewis: the +ASC correction instead of the
        # constant-site counts (-fconst) the typing pass provides. alrt=N:
        # N SH-aLRT replicates (-alrt), supports without the ultrafast
        # bootstrap; with -B as well the tree carries both.
        accepted_extras=frozenset({"model", "asc", "alrt"}),
        # 2.0: -B (ultrafast bootstrap) and --threads-max as passed below.
        required_binaries=(BinarySpec("iqtree", version_args=("--version",), min_version="2.0"),),
        recommended_max_genomes=500,
    )
    input_kind = InputKind.MSA_FASTA
    # Without -B or alrt= the tree carries no support values at all.
    supports_without_bootstrap = False
    # -p: one partition per gene, edge-linked with proportional branch lengths.
    supports_partitions = True
    # The model of every partition when --partitioned is given without model=.
    PARTITION_MODEL_DNA = "GTR+G"
    PARTITION_MODEL_PROTEIN = "LG+G"

    # IQ-TREE's ultrafast bootstrap (-B) rejects fewer than 1000 replicates.
    MIN_BOOTSTRAP = 1000
    # IQ-TREE's SH-aLRT test (-alrt) takes at least 1000 replicates as well.
    MIN_ALRT = 1000
    # Above this many sequences ModelFinder takes a noticeable share of the
    # run, and the INFO line says so.
    MODELFINDER_NOTE_LEAVES = 200

    def writes_supports(self, extra: Mapping[str, object]) -> bool:
        return "alrt" in extra

    def _note_modelfinder(self, msa: Path, logger: logging.Logger, *, quiet_below: bool) -> None:
        """Say that ModelFinder runs, with the sequence count above
        MODELFINDER_NOTE_LEAVES; below it only when ``quiet_below`` is False."""
        leaves = _count_records(msa)
        if leaves > self.MODELFINDER_NOTE_LEAVES:
            logger.info(
                "iqtree: no model given, so IQ-TREE runs ModelFinder on %d sequences "
                "before the tree search, which at this size takes a noticeable share of "
                "the run; model=GTR+G (--aligner-arg) skips the selection",
                leaves,
            )
        elif not quiet_below:
            logger.info(
                "iqtree: no model given, so IQ-TREE runs ModelFinder before the tree "
                "search; model=GTR+G (--aligner-arg) skips the selection"
            )

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
        if params.partitions is not None:
            refuse_partitions_with_asc(params.extra, params.constant_sites, "iqtree")
        refuse_asc_on_protein(params, "iqtree")
        alrt = _alrt(params.extra, self.MIN_ALRT)
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
            # The character type is given, not guessed: left to guess,
            # IQ-TREE refuses a recombination-masked alignment, where a fifth
            # of the characters can be N ("Unknown sequence type"). Every
            # alignment is nucleotide except a protein core-gene alignment.
            "-st",
            "AA" if params.alphabet == "protein" else "DNA",
            "-redo",  # overwrite checkpoints from a previous run at this path
        ]
        if asc == "lewis":
            # +ASC needs a model term; MFP keeps ModelFinder when none is given.
            cmd += ["-m", f"{model or 'MFP'}+ASC"]
            logger.info("iqtree: ascertainment-bias correction +ASC (asc=lewis)")
            if not model:
                self._note_modelfinder(work_msa, logger, quiet_below=True)
        elif model:
            cmd += ["-m", model]
            if params.partitions is not None:
                logger.info("iqtree: model %s applied to every partition", model)
        elif params.partitions is not None:
            # The partition file holds charsets only, so -m sets the model of
            # every partition. A fixed default keeps hundreds of partitions
            # from each running ModelFinder; model=MFP requests that instead.
            default = (
                self.PARTITION_MODEL_PROTEIN
                if params.alphabet == "protein"
                else self.PARTITION_MODEL_DNA
            )
            cmd += ["-m", default]
            logger.info(
                "iqtree: no model given, so %s is applied to every partition; "
                "model=MFP (--aligner-arg) selects a model for each partition, "
                "which is slow with many partitions",
                default,
            )
        else:
            self._note_modelfinder(work_msa, logger, quiet_below=False)
        if not asc and constant is not None:
            cmd += ["-fconst", format_counts(constant, ",")]
            logger.info("iqtree: constant sites added as -fconst %s", format_counts(constant, ","))
        if params.outgroup:
            cmd += ["-o", params.outgroup]
        if params.bootstrap > 0:
            cmd += ["-B", str(params.bootstrap)]
        if alrt is not None:
            cmd += ["-alrt", str(alrt)]
            logger.info("iqtree: SH-aLRT supports from %d replicates (-alrt)", alrt)
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


def _alrt(extra: dict, minimum: int) -> int | None:
    """The SH-aLRT replicate count of ``alrt=``, None when unset."""
    if "alrt" not in extra:
        return None
    raw = extra["alrt"]
    try:
        value = int(str(raw).strip())
    except ValueError:
        raise UserInputError(
            f"iqtree alrt= must be a whole number of replicates, got {raw!r}"
        ) from None
    if value < minimum:
        raise UserInputError(
            f"iqtree alrt= needs at least {minimum} SH-aLRT replicates (-alrt); got {value}."
        )
    return value


def _count_records(path: Path) -> int:
    """The number of FASTA records in ``path``."""
    with open(path, "rb") as fh:
        return sum(1 for line in fh if line.startswith(b">"))


def _write_partition_file(params: TreeParams, out_dir: Path) -> Path:
    """Write ``out_dir/partitions.txt`` as a NEXUS sets block of charsets.

    The file names no model: IQ-TREE reads the first field of a RAxML-style
    line as a model (``DNA`` as GTR, ``AA`` as a model file), so the model of
    every partition is given with ``-m`` instead.
    """
    assert params.partitions is not None
    names = partition_names([gene for gene, _, _ in params.partitions])
    lines = ["#nexus", "begin sets;"]
    lines += [
        f"    charset {name} = {start}-{end};"
        for name, (_, start, end) in zip(names, params.partitions, strict=True)
    ]
    lines.append("end;")
    path = out_dir / "partitions.txt"
    with atomic_path(path) as tmp:
        tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
