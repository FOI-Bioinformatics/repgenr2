"""prokka annotator (no database to download).

Command surface checked against prokka 1.15.6 (bioconda): ``prokka --outdir
DIR --prefix NAME --cpus N --force [--kingdom K --genus G --species S
--compliant --centre C --gcode N] GENOME`` writes ``NAME.gff`` (FASTA
appended), ``NAME.faa`` and ``NAME.ffn``; the GFF is copied as
``NAME.gff3`` so every annotator leaves the same file names (prokka's own
``NAME.gff`` is left beside it, as it writes it).

prokka refuses a contig name longer than 37 characters unless it renames
the contigs (``compliant=true`` with ``centre=X``), which assemblies from
the reads route hit (``NODE_1_length_..._cov_...``); the failure is caught
and says so. ``compliant`` also drops contigs under 200 bases.
"""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import AnnotateParams, AnnotateResult, Annotator, require_fasta_section

_VALUE_FLAGS = ("kingdom", "genus", "species", "strain", "centre", "gcode")


_PROKKA_MAX_CONTIG_NAME = 37


def _long_contig_name(genome: Path) -> bool:
    """Whether a FASTA header's first token is too long for prokka."""
    try:
        with open(genome, "rb") as fo:
            for raw in fo:
                if raw.startswith(b">"):
                    name = raw[1:].split()[0] if raw[1:].split() else b""
                    if len(name) > _PROKKA_MAX_CONTIG_NAME:
                        return True
    except OSError:
        return False
    return False


class ProkkaAnnotator(Annotator):
    capabilities = ToolCapabilities(
        name="prokka",
        container="quay.io/biocontainers/prokka:1.15.6--pl5321hdfd78af_0",
        conda=("bioconda::prokka",),
        required_binaries=(BinarySpec("prokka", version_args=("--version",), min_version="1.14"),),
        recommended_max_genomes=None,
        auto_eligible=False,
        accepted_extras=frozenset({*_VALUE_FLAGS, "compliant"}),
    )

    def call(
        self,
        genome: Path,
        record_name: str,
        out_dir: Path,
        params: AnnotateParams,
        logger: logging.Logger,
    ) -> AnnotateResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd: list[str | Path] = [
            "prokka",
            "--outdir",
            out_dir,
            "--prefix",
            record_name,
            "--cpus",
            str(params.threads),
            "--force",
        ]
        if parse_extra_bool(params.extra, "compliant", False):
            cmd.append("--compliant")
        for key in _VALUE_FLAGS:
            if key in params.extra and str(params.extra[key]).strip():
                cmd += ["--" + key, str(params.extra[key])]
        cmd.append(genome)
        try:
            run_tool(
                self.capabilities,
                cmd,
                logger=logger,
                log_prefix=f"prokka {record_name}",
                extra_mounts=[os.path.dirname(os.path.abspath(genome))],
            )
        except ToolExecutionError as exc:
            if _long_contig_name(genome):
                raise WorkdirError(
                    f"prokka failed on {genome.name}, whose contig names are longer than the "
                    "37 characters it accepts (assembler contig names usually are). Pass "
                    "--tool-arg compliant=true --tool-arg centre=X to let prokka rename them "
                    "(it then also drops contigs under 200 bases), or annotate with bakta."
                ) from exc
            raise
        gff = out_dir / f"{record_name}.gff"
        if not gff.is_file():
            raise WorkdirError(f"prokka wrote no {gff.name} in {out_dir}")
        gff3 = out_dir / f"{record_name}.gff3"
        shutil.copy2(gff, gff3)
        require_fasta_section(gff3, "prokka")
        faa = out_dir / f"{record_name}.faa"
        ffn = out_dir / f"{record_name}.ffn"
        return AnnotateResult(
            gff=gff3, faa=faa if faa.is_file() else None, ffn=ffn if ffn.is_file() else None
        )
