"""bakta annotator (bacterial genomes; needs a bakta database).

Command surface checked against bakta 1.12 (bioconda): ``bakta --db DB
--prefix NAME --output DIR --threads N --force --keep-contig-headers
[--skip-plot ...] GENOME`` writes ``NAME.gff3`` (FASTA appended),
``NAME.faa`` and ``NAME.ffn``. The database comes from ``bakta_db download
--output DIR --type light|full`` (light about 4 GB, full about 84 GB);
``DIR/version.json`` names its version.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_bool
from .base import AnnotateParams, AnnotateResult, Annotator, require_fasta_section

BAKTA_DB_ENV = "BAKTA_DB"
# bakta's own --skip-* switches, as extras (skip_plot=true by default: the
# circular plot is slow and not read by anything here).
_SKIP_FLAGS = (
    "skip_trna",
    "skip_tmrna",
    "skip_rrna",
    "skip_ncrna",
    "skip_ncrna_region",
    "skip_crispr",
    "skip_pseudo",
    "skip_sorf",
    "skip_gap",
    "skip_ori",
    "skip_plot",
)
_VALUE_FLAGS = ("genus", "species", "strain", "translation_table")


def bakta_db_from_env() -> str | None:
    return os.environ.get(BAKTA_DB_ENV) or None


def check_bakta_db(db: Path | None) -> Path:
    """The bakta database directory, from the flag or ``BAKTA_DB``; refused
    with the download command when absent or incomplete."""
    hint = (
        "bakta needs its database: run `bakta_db download --output DIR --type light` "
        "(about 4 GB; `--type full` for the complete one) and pass --bakta-db DIR or set "
        f"{BAKTA_DB_ENV}; prokka (--tool prokka) needs no database."
    )
    if db is None:
        raise UserInputError(hint)
    if not db.is_dir():
        raise UserInputError(f"--bakta-db {db} is not a directory. {hint}")
    if not (db / "version.json").is_file():
        raise UserInputError(f"{db} holds no version.json; not a bakta database. {hint}")
    return db


def read_db_version(db: Path) -> str:
    """``<major>.<minor>[-<type>]`` from the database's version.json.

    The ``type`` field (light or full) exists from database version 5.1; an
    older file gives the version alone.
    """
    path = db / "version.json"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise WorkdirError(f"{path} is unreadable: {exc}") from exc
    if not isinstance(data, dict):
        raise WorkdirError(f"{path} does not hold a JSON object; not a bakta database.")
    major, minor = data.get("major", ""), data.get("minor", "")
    kind = data.get("type", "")
    return f"{major}.{minor}{'-' + kind if kind else ''}"


class BaktaAnnotator(Annotator):
    capabilities = ToolCapabilities(
        name="bakta",
        container="quay.io/biocontainers/bakta:1.12.2--pyhdfd78af_0",
        conda=("bioconda::bakta",),
        required_binaries=(BinarySpec("bakta", version_args=("--version",), min_version="1.9"),),
        # Per genome, not per set; the annotate stage runs it once per record.
        recommended_max_genomes=None,
        auto_eligible=False,
        default_params={"skip_plot": True},
        accepted_extras=frozenset({*_SKIP_FLAGS, *_VALUE_FLAGS, "compliant"}),
    )
    requires_database = True

    def database_version(self, db: Path) -> str:
        return read_db_version(db)

    def call(
        self,
        genome: Path,
        record_name: str,
        out_dir: Path,
        params: AnnotateParams,
        logger: logging.Logger,
    ) -> AnnotateResult:
        db = check_bakta_db(params.db)
        out_dir.mkdir(parents=True, exist_ok=True)
        extra = {**self.capabilities.default_params, **params.extra}
        if parse_extra_bool(extra, "skip_cds", False):
            raise UserInputError(
                "bakta skip_cds leaves no protein-coding genes, which the core-gene route "
                "needs; drop it, or annotate with another tool."
            )
        cmd: list[str | Path] = [
            "bakta",
            "--db",
            db,
            "--prefix",
            record_name,
            "--output",
            out_dir,
            "--threads",
            str(params.threads),
            "--force",
            "--keep-contig-headers",
        ]
        for key in _SKIP_FLAGS:
            if parse_extra_bool(extra, key, False):
                cmd.append("--" + key.replace("_", "-"))
        if parse_extra_bool(extra, "compliant", False):
            cmd.append("--compliant")
        for key in _VALUE_FLAGS:
            if key in extra and str(extra[key]).strip():
                cmd += ["--" + key.replace("_", "-"), str(extra[key])]
        cmd.append(genome)
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix=f"bakta {record_name}",
            extra_mounts=[db, os.path.dirname(os.path.abspath(genome))],
        )
        gff = out_dir / f"{record_name}.gff3"
        if not gff.is_file():
            raise WorkdirError(f"bakta wrote no {gff.name} in {out_dir}")
        require_fasta_section(gff, "bakta")
        faa = out_dir / f"{record_name}.faa"
        ffn = out_dir / f"{record_name}.ffn"
        return AnnotateResult(
            gff=gff, faa=faa if faa.is_file() else None, ffn=ffn if ffn.is_file() else None
        )
