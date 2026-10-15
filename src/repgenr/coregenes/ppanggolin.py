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

PPanGGOLiN aligns only the genomes that carry exactly one copy of a family,
so a genome with several copies is a row of gaps for that family. With
``--single_copy``, passed by default, families multi-copy in at least
``dup_margin`` (0.05) of the genomes are left out; a genome with several
copies of a family that is kept is still a gap for it. ``single_copy=false``
keeps every family. A genome that shares no selected family at all is a row
of gaps throughout, and the stage drops it from the alignment.
"""

from __future__ import annotations

import contextlib
import logging
import os
import shutil
import tempfile
from collections.abc import Iterator, Sequence
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
        conda=("bioconda::ppanggolin", "bioconda::mmseqs2"),
        required_binaries=(
            BinarySpec("ppanggolin", version_args=("--version",), min_version="2.0"),
            BinarySpec("mmseqs", version_args=("version",), min_version="13.0"),
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
        dup_margin = (
            _unit(params.extra["dup_margin"], "dup_margin")
            if "dup_margin" in params.extra
            else None
        )
        if dup_margin is not None and not single_copy:
            logger.warning("ppanggolin dup_margin applies only with single_copy; it is ignored.")
        if partition != "softcore" and params.core_threshold != CoreGeneParams.core_threshold:
            logger.warning(
                "partition=%s does not read --core-threshold; the value %g is ignored.",
                partition,
                params.core_threshold,
            )
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
        with _partition_tmpdir(out_dir / "tmp", logger) as tmp_dir:
            self._run(
                params, partition, identity, coverage, single_copy, dup_margin,
                listing, pan_dir, h5, out_dir, tmp_dir, [*mounts, str(tmp_dir)], common, logger,
            )  # fmt: skip

        soft_core = f"{params.core_threshold:g}"
        tables_dir = out_dir / "tables"
        msa_dir = out_dir / "msa"
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
            settings={
                "partition": partition,
                "coverage": coverage,
                "single_copy": single_copy,
                # MMseqs2 compares proteins; Panaroo's identity is nucleotide.
                "clustering_on": "protein",
            },
        )

    def _run(
        self,
        params: CoreGeneParams,
        partition: str,
        identity: float,
        coverage: float,
        single_copy: bool,
        dup_margin: float | None,
        listing: Path,
        pan_dir: Path,
        h5: Path,
        out_dir: Path,
        tmp_dir: Path,
        mounts: list[str],
        common: list[str],
        logger: logging.Logger,
    ) -> None:
        """The three PPanGGOLiN calls: workflow, write_pangenome and msa."""
        workflow: list[str | Path] = [
            "ppanggolin", "workflow", "--anno", listing, "-o", pan_dir,
            "--identity", f"{identity:g}", "--coverage", f"{coverage:g}",
            "--no_flat_files", "--tmpdir", tmp_dir, *common,
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
            "--source", "dna", "--phylo", "--tmpdir", tmp_dir, *common,
        ]  # fmt: skip
        if partition == "softcore":
            msa += ["--soft_core", soft_core]
        if single_copy:
            msa.append("--single_copy")
            if dup_margin is not None:
                msa += ["--dup_margin", f"{dup_margin:g}"]
        run_tool(
            self.capabilities, msa, logger=logger, cwd=out_dir,
            log_prefix="ppanggolin", extra_mounts=mounts,
        )  # fmt: skip


# NEM, the partitioning code PPanGGOLiN bundles, keeps file names in buffers of
# 200 characters, and PPanGGOLiN builds them under --tmpdir with about 40 more.
# A longer name makes a partitioning worker fail and the whole run wait for it
# without end (seen under a Nextflow work directory). Beyond this length the
# adapter uses a short temporary directory of its own.
_TMPDIR_MAX = 120


@contextlib.contextmanager
def _partition_tmpdir(preferred: Path, logger: logging.Logger) -> Iterator[Path]:
    """``preferred`` when its path is short enough for NEM, else a short one
    under the system temporary directory, removed afterwards."""
    if len(os.path.abspath(preferred)) <= _TMPDIR_MAX:
        preferred.mkdir(parents=True, exist_ok=True)
        yield preferred
        return
    own = Path(tempfile.mkdtemp(prefix="repgenr-ppg-"))
    if len(str(own)) > _TMPDIR_MAX:
        shutil.rmtree(own, ignore_errors=True)
        raise WorkdirError(
            f"ppanggolin needs a temporary directory path of at most {_TMPDIR_MAX} "
            f"characters (its partitioning code keeps file names in fixed buffers), and "
            f"both {preferred} and the system one are longer; set TMPDIR to a short path."
        )
    logger.info(
        "The work path is too long for ppanggolin's partitioning code; its temporary "
        "files go to %s instead",
        own,
    )
    try:
        yield own
    finally:
        shutil.rmtree(own, ignore_errors=True)


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
