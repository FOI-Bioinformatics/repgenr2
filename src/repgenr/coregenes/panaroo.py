"""Panaroo core-gene alignment builder.

Command surface checked against Panaroo 1.8 (bioconda): ``panaroo -i
LIST.txt -o DIR --clean-mode MODE -a core --aligner mafft --core_threshold
T -t N`` (``-i`` takes the GFF3 paths, or a text file listing them one per
line, which keeps the argument vector short at thousands of genomes). It
writes ``core_gene_alignment.aln``, ``core_gene_alignment_filtered.aln``
(genes failing its entropy filter removed; its recommended tree input),
``gene_presence_absence.csv``, ``gene_presence_absence.Rtab``,
``summary_statistics.txt`` and a header file naming the gene blocks of the
alignment (``core_alignment_header.embl``, and the filtered alignment's own
when it writes one). The genes of the alignment are counted from that
header: ``summary_statistics.txt`` reports Panaroo's fixed 99 to 100 per
cent bin, which is not the set kept at another ``--core_threshold``. mafft
does the per-gene alignments.
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import CLEAN_MODES, CoreGeneBuilder, CoreGeneParams, CoreGeneResult

TABLES = (
    "gene_presence_absence.csv",
    "gene_presence_absence.Rtab",
    "summary_statistics.txt",
    "core_alignment_header.embl",
    "core_alignment_filtered_header.embl",
)
_FEATURE = re.compile(r"^FT\s+feature\s")


class PanarooBuilder(CoreGeneBuilder):
    capabilities = ToolCapabilities(
        name="panaroo",
        container="quay.io/biocontainers/panaroo:1.8.0--pyhdfd78af_0",
        conda=("bioconda::panaroo",),
        required_binaries=(
            BinarySpec("panaroo", version_args=("--version",), min_version="1.3"),
            BinarySpec("mafft", version_args=("--version",), min_version="7.0"),
        ),
        # Declared, not measured here: the all-against-all clustering and the
        # per-gene mafft alignments grow with the genome count.
        recommended_max_genomes=500,
        auto_eligible=False,
        accepted_extras=frozenset(
            {"remove_invalid_genes", "merge_paralogs", "family_threshold", "len_dif_percent"}
        ),
    )

    def call(
        self,
        gffs: Sequence[Path],
        out_dir: Path,
        params: CoreGeneParams,
        logger: logging.Logger,
    ) -> CoreGeneResult:
        if params.clean_mode not in CLEAN_MODES:
            raise UserInputError(
                f"--clean-mode must be one of {', '.join(CLEAN_MODES)}, not {params.clean_mode!r}"
            )
        if not 0.0 < params.core_threshold <= 1.0:
            raise UserInputError(
                f"--core-threshold must be within (0, 1], got {params.core_threshold}"
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        listing = out_dir / "gff_list.txt"
        listing.write_text("".join(f"{os.path.abspath(g)}\n" for g in gffs), encoding="utf-8")
        cmd: list[str | Path] = [
            "panaroo",
            "-i",
            listing,
            "-o",
            out_dir,
            "--clean-mode",
            params.clean_mode,
            "-a",
            "core",
            "--aligner",
            "mafft",
            "--core_threshold",
            f"{params.core_threshold:g}",
            "-t",
            str(params.threads),
        ]
        if parse_extra_bool(params.extra, "remove_invalid_genes", False):
            cmd.append("--remove_invalid_genes")
        if parse_extra_bool(params.extra, "merge_paralogs", False):
            cmd.append("--merge_paralogs")
        for key in ("family_threshold", "len_dif_percent"):
            if key in params.extra:
                cmd += ["--" + key, str(params.extra[key])]
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            cwd=out_dir,
            log_prefix="panaroo",
            extra_mounts=sorted({os.path.dirname(os.path.abspath(g)) for g in gffs}),
        )
        name = "core_gene_alignment_filtered.aln" if params.filtered else "core_gene_alignment.aln"
        alignment = out_dir / name
        if not alignment.is_file() or alignment.stat().st_size == 0:
            if params.filtered and (out_dir / "core_gene_alignment.aln").is_file():
                raise WorkdirError(
                    "panaroo wrote core_gene_alignment.aln but no filtered alignment; pass "
                    "--aligner-arg filtered=false to use the unfiltered one"
                )
            raise WorkdirError(
                f"panaroo wrote no {name}: too few genomes share core genes at "
                f"--core-threshold {params.core_threshold:g}, or the GFF3 files carry no sequences"
            )
        tables = [out_dir / t for t in TABLES if (out_dir / t).is_file()]
        header = out_dir / (
            "core_alignment_filtered_header.embl"
            if params.filtered
            else "core_alignment_header.embl"
        )
        return CoreGeneResult(
            alignment=alignment, tables=tables, core_genes=_core_genes(header, out_dir, logger)
        )


def _core_genes(header: Path, out_dir: Path, logger: logging.Logger) -> int | None:
    """The gene blocks of the alignment, counted in Panaroo's header file.

    Falls back to the unfiltered header, then to nothing: the count in
    ``summary_statistics.txt`` is Panaroo's 99 to 100 per cent bin and would
    misreport the genes kept at another threshold.
    """
    for candidate in (header, out_dir / "core_alignment_header.embl"):
        if candidate.is_file():
            genes = sum(
                1
                for line in candidate.read_text(encoding="utf-8").splitlines()
                if _FEATURE.match(line)
            )
            if candidate != header:
                logger.info(
                    "panaroo wrote no %s, so the gene count (%d) comes from %s and counts the "
                    "genes before the entropy filter",
                    header.name,
                    genes,
                    candidate.name,
                )
            return genes
    logger.info("panaroo wrote no alignment header; the core-gene count is unknown")
    return None
