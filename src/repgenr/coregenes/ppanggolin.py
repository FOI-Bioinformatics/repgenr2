"""PPanGGOLiN core-gene alignment builder.

Command surface checked against PPanGGOLiN 2.3.2 (bioconda). Three calls:

1. ``ppanggolin workflow --anno LIST.tsv -o DIR -c N --identity I --coverage C
   --no_flat_files`` clusters the proteins with MMseqs2 (80 per cent identity
   and coverage by default), builds the graph and partitions it. ``LIST.tsv``
   names each genome (first column, the record name) and its GFF3 (second);
   the sequences come from the ``##FASTA`` section the annotate stage keeps.
2. ``ppanggolin write_pangenome -p pangenome.h5 --Rtab --csv --stats
   --partitions`` writes the presence tables published beside the alignment.
3. ``ppanggolin msa -p pangenome.h5 --partition P --source dna --phylo
   --single_copy`` aligns each selected family with mafft and concatenates
   them into ``<partition>[_<soft core>]_genome_alignment.aln``, one record
   per genome, gaps where a genome lacks the family.

The partition defaults to ``softcore`` at ``--core-threshold``, so the gene
set means the same as with Panaroo. ``partition=persistent`` selects
PPanGGOLiN's statistical partition instead, which does not shrink as genomes
are added; it needs a set of some size (the tool warns below 15 genomes, and
on three genomes it put every family in it). ``partition=core`` is the strict
core.

``--single_copy`` is passed by default: in the concatenation PPanGGOLiN keeps
the first copy of a duplicated gene and replaces the family with gaps for
that genome, so a multi-copy family adds missing data rather than signal.
``single_copy=false`` keeps them.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import CoreGeneBuilder, CoreGeneParams, CoreGeneResult

PARTITIONS = ("softcore", "core", "persistent")
TABLES = (
    "gene_presence_absence.Rtab",
    "matrix.csv",
    "genomes_statistics.tsv",
    "mean_persistent_duplication.tsv",
)
# PPanGGOLiN's own clustering defaults, named so the messages can quote them.
DEFAULT_IDENTITY = 0.8
DEFAULT_COVERAGE = 0.8
# Below this the statistical partition is not robust (PPanGGOLiN's guidance).
_PERSISTENT_MIN_GENOMES = 15
# Below this a core-gene tree rests on very little.
_SMALL_CORE = 100


class PPanGGOLiNBuilder(CoreGeneBuilder):
    capabilities = ToolCapabilities(
        name="ppanggolin",
        container="quay.io/biocontainers/ppanggolin:2.3.2--py312hfabe715_0",
        conda=("bioconda::ppanggolin",),
        required_binaries=(
            BinarySpec("ppanggolin", version_args=("--version",), min_version="2.0"),
            BinarySpec("mafft", version_args=("--version",), min_version="7.0"),
        ),
        # Declared, not measured here. PPanGGOLiN is built for tens of
        # thousands of genomes; the per-family mafft alignments are the part
        # that grows fastest.
        recommended_max_genomes=5000,
        auto_eligible=False,
        accepted_extras=frozenset(
            {
                "partition",
                # MMseqs2 clustering: protein identity and alignment coverage.
                "identity",
                "coverage",
                "mode",
                "nb_of_partitions",
                "translation_table",
                "single_copy",
                "dup_margin",
            }
        ),
        # Panaroo's graph cleaning and entropy filter have no counterpart.
        ignored_params=frozenset({"clean_mode", "filtered"}),
    )

    def call(
        self,
        gffs: Sequence[Path],
        out_dir: Path,
        params: CoreGeneParams,
        logger: logging.Logger,
    ) -> CoreGeneResult:
        partition = str(params.extra.get("partition", "softcore"))
        if partition not in PARTITIONS:
            raise UserInputError(
                f"ppanggolin partition must be one of {', '.join(PARTITIONS)}, not {partition!r}"
            )
        if not 0.0 < params.core_threshold <= 1.0:
            raise UserInputError(
                f"--core-threshold must be within (0, 1], got {params.core_threshold}"
            )
        identity = _unit(params.extra.get("identity", DEFAULT_IDENTITY), "identity")
        coverage = _unit(params.extra.get("coverage", DEFAULT_COVERAGE), "coverage")
        single_copy = parse_extra_bool(params.extra, "single_copy", True)
        if partition == "persistent" and len(gffs) < _PERSISTENT_MIN_GENOMES:
            logger.warning(
                "partition=persistent on %d genomes: PPanGGOLiN's statistical partition is not "
                "robust below %d, and on small sets it places nearly every family in it. "
                "The default, softcore at --core-threshold, is the safer choice here.",
                len(gffs),
                _PERSISTENT_MIN_GENOMES,
            )

        out_dir.mkdir(parents=True, exist_ok=True)
        listing = out_dir / "annotations.tsv"
        listing.write_text(
            "".join(f"{g.stem}\t{os.path.abspath(g)}\n" for g in gffs), encoding="utf-8"
        )
        mounts = sorted({os.path.dirname(os.path.abspath(g)) for g in gffs})
        pan_dir = out_dir / "pangenome"
        h5 = pan_dir / "pangenome.h5"
        common = ["-c", str(params.threads), "-f", "-d", "--verbose", "1"]

        workflow: list[str | Path] = [
            "ppanggolin", "workflow", "--anno", listing, "-o", pan_dir,
            "--identity", f"{identity:g}", "--coverage", f"{coverage:g}",
            "--no_flat_files", *common,
        ]  # fmt: skip
        for key, flag in (
            ("mode", "--mode"),
            ("nb_of_partitions", "-K"),
            ("translation_table", "--translation_table"),
        ):
            if key in params.extra:
                workflow += [flag, str(params.extra[key])]
        run_tool(
            self.capabilities, workflow, logger=logger, cwd=out_dir,
            log_prefix="ppanggolin", extra_mounts=mounts,
        )  # fmt: skip

        soft_core = f"{params.core_threshold:g}"
        tables_dir = out_dir / "tables"
        run_tool(
            self.capabilities,
            [
                "ppanggolin", "write_pangenome", "-p", h5, "-o", tables_dir,
                "--Rtab", "--csv", "--stats", "--partitions", "--soft_core", soft_core,
                *common,
            ],
            logger=logger, cwd=out_dir, log_prefix="ppanggolin",
        )  # fmt: skip

        msa_dir = out_dir / "msa"
        msa: list[str | Path] = [
            "ppanggolin", "msa", "-p", h5, "-o", msa_dir, "--partition", partition,
            "--source", "dna", "--phylo", *common,
        ]  # fmt: skip
        if partition == "softcore":
            msa += ["--soft_core", soft_core]
        if single_copy:
            msa.append("--single_copy")
            if "dup_margin" in params.extra:
                msa += ["--dup_margin", f"{_unit(params.extra['dup_margin'], 'dup_margin'):g}"]
        run_tool(self.capabilities, msa, logger=logger, cwd=out_dir, log_prefix="ppanggolin")

        alignments = sorted(msa_dir.glob("*_genome_alignment.aln"))
        families = [p for d in msa_dir.glob("msa_*_dna") for p in d.glob("*.aln") if p.is_file()]
        core_genes = len(families)
        if not alignments or core_genes == 0 or not _has_columns(alignments[0]):
            raise WorkdirError(
                f"ppanggolin selected no gene family for the {partition} partition"
                + (f" at --core-threshold {soft_core}" if partition == "softcore" else "")
                + f". Its clustering groups proteins at least {identity:g} identical over "
                f"{coverage:g} of their length (--aligner-arg identity= and coverage=); on a "
                "divergent set lower them, or lower --core-threshold. Check too that the GFF3 "
                "files carry their sequences."
            )
        if core_genes < _SMALL_CORE:
            logger.warning(
                "Only %d gene families in the %s alignment. On a divergent set lower "
                "--aligner-arg identity= (now %g) or --core-threshold.",
                core_genes,
                partition,
                identity,
            )
        tables = [tables_dir / t for t in TABLES if (tables_dir / t).is_file()]
        partitions_dir = tables_dir / "partitions"
        tables += sorted(p for p in partitions_dir.glob("*.txt") if p.is_file())
        return CoreGeneResult(
            alignment=alignments[0],
            tables=tables,
            core_genes=core_genes,
            clustering_identity=identity,
        )


def _unit(raw: object, flag: str) -> float:
    try:
        value = float(str(raw))
    except ValueError:
        raise UserInputError(
            f"ppanggolin {flag} must be a number within (0, 1], got {raw!r}"
        ) from None
    if not 0.0 < value <= 1.0:
        raise UserInputError(f"ppanggolin {flag} must be within (0, 1], got {raw!r}")
    return value


def _has_columns(alignment: Path) -> bool:
    """True when the concatenation holds at least one aligned column."""
    with open(alignment, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            if line.strip():
                return True
    return False
