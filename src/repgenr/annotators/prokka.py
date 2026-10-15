"""prokka annotator (no database to download).

Command surface checked against prokka 1.14 (bioconda): ``prokka --outdir
DIR --prefix NAME --cpus N --force [--kingdom K --genus G --species S
--compliant --centre C --gcode N] GENOME`` writes ``NAME.gff`` (FASTA
appended), ``NAME.faa`` and ``NAME.ffn``; the GFF is copied as
``NAME.gff3`` so every annotator leaves the same file names.
"""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import AnnotateParams, AnnotateResult, Annotator, require_fasta_section

_VALUE_FLAGS = ("kingdom", "genus", "species", "strain", "centre", "gcode")


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
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="prokka",
            extra_mounts=[os.path.dirname(os.path.abspath(genome))],
        )
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
