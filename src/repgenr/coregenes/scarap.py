"""SCARAP core-gene alignment builder.

Command surface checked against SCARAP 1.0.2 (bioconda). Two calls:

1. ``scarap core FAA_LIST OUT -f F -t N`` infers the core genome directly,
   without the full pangenome: it builds orthogroups on a set of seed
   genomes, turns them into MMseqs2 profiles and searches every genome
   against them. A gene is core when it is present in a single copy in at
   least ``F`` of the genomes (``--core-threshold``). It writes
   ``genes.tsv`` (gene, genome, orthogroup).
2. ``scarap concat FAA_LIST genes.tsv OUT -n FFN_LIST -t N`` aligns each
   orthogroup with mafft and concatenates them; with ``-n`` it writes the
   nucleotide supermatrix ``supermatrix_nucs.fasta`` beside the amino-acid
   one, one record per genome.

SCARAP reads proteins, not GFF3: the adapter takes the ``.faa`` and ``.ffn``
files the annotate stage writes beside each GFF3, and SCARAP names each
genome after its file stem, which is the record name. The orthogroups come
from profile searches rather than a fixed identity cut, so the stats file
records no clustering identity.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
from .base import CoreGeneBuilder, CoreGeneParams, CoreGeneResult

METHODS = ("H", "FH", "FT", "H-nl", "T-nl", "P", "O-B", "O-D", "S")
SUPERMATRIX = "supermatrix_nucs.fasta"
TABLES = ("genes.tsv",)
# Below this a core-gene tree rests on very little.
_SMALL_CORE = 100


class ScarapBuilder(CoreGeneBuilder):
    capabilities = ToolCapabilities(
        name="scarap",
        container="quay.io/biocontainers/scarap:1.0.2--pyhdfd78af_0",
        conda=("bioconda::scarap", "bioconda::mmseqs2"),
        required_binaries=(
            # scarap has no --version; its help opens with the version.
            BinarySpec("scarap", version_args=("-h",), min_version="1.0"),
            BinarySpec("mmseqs", version_args=("version",), min_version="13.0"),
            BinarySpec("mafft", version_args=("--version",), min_version="7.0"),
        ),
        # Declared, not measured here. The authors ran it on over 31000
        # genomes; the per-orthogroup alignments are the part that grows.
        recommended_max_genomes=10000,
        auto_eligible=False,
        accepted_extras=frozenset({"method", "core_prefilter", "max_core_genes", "seeds"}),
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
        if not 0.0 < params.core_threshold <= 1.0:
            raise UserInputError(
                f"--core-threshold must be within (0, 1], got {params.core_threshold}"
            )
        method = str(params.extra.get("method", "FH"))
        if method not in METHODS:
            raise UserInputError(
                f"scarap method must be one of {', '.join(METHODS)}, not {method!r}"
            )
        prefilter = (
            _unit(params.extra["core_prefilter"], "core_prefilter")
            if "core_prefilter" in params.extra
            else None
        )
        max_core = parse_extra_int(params.extra, "max_core_genes", 0)
        seeds = parse_extra_int(params.extra, "seeds", 100)
        if max_core < 0 or seeds < 1:
            raise UserInputError(
                "scarap max_core_genes must be 0 or more and seeds 1 or more, "
                f"got {max_core} and {seeds}"
            )
        faas = [g.with_suffix(".faa") for g in gffs]
        ffns = [g.with_suffix(".ffn") for g in gffs]
        missing = [str(p) for p in (*faas, *ffns) if not p.is_file()]
        if missing:
            raise WorkdirError(
                f"scarap reads the protein and gene files the annotate stage writes beside "
                f"each GFF3, and {len(missing)} are missing ({', '.join(missing[:3])}); "
                "rerun `repgenr annotate --force`."
            )

        out_dir.mkdir(parents=True, exist_ok=True)
        faa_list = out_dir / "faa_files.txt"
        ffn_list = out_dir / "ffn_files.txt"
        faa_list.write_text("".join(f"{os.path.abspath(p)}\n" for p in faas), encoding="utf-8")
        ffn_list.write_text("".join(f"{os.path.abspath(p)}\n" for p in ffns), encoding="utf-8")
        mounts = sorted({os.path.dirname(os.path.abspath(g)) for g in gffs})
        core_dir = out_dir / "core"
        concat_dir = out_dir / "concat"

        core: list[str | Path] = [
            "scarap", "core", faa_list, core_dir,
            "-f", f"{params.core_threshold:g}", "-d", method,
            "-m", str(max_core), "-e", str(seeds), "-t", str(params.threads),
        ]  # fmt: skip
        if prefilter is not None:
            core += ["-p", f"{prefilter:g}"]
        run_tool(
            self.capabilities, core, logger=logger, cwd=out_dir,
            log_prefix="scarap", extra_mounts=mounts,
        )  # fmt: skip

        genes = core_dir / "genes.tsv"
        orthogroups = _orthogroups(genes)
        if not orthogroups:
            raise WorkdirError(
                f"scarap found no core gene present in a single copy in "
                f"{params.core_threshold:g} of the genomes (--core-threshold). Lower it, or "
                "leave the most divergent genomes out."
            )
        run_tool(
            self.capabilities,
            [
                "scarap", "concat", faa_list, genes, concat_dir,
                "-n", ffn_list, "-t", str(params.threads),
            ],
            logger=logger, cwd=out_dir, log_prefix="scarap", extra_mounts=mounts,
        )  # fmt: skip
        alignment = concat_dir / SUPERMATRIX
        if not alignment.is_file() or alignment.stat().st_size == 0:
            raise WorkdirError(
                f"scarap concat wrote no {SUPERMATRIX}; its log under {concat_dir} says why."
            )
        if len(orthogroups) < _SMALL_CORE:
            logger.warning(
                "Only %d core genes in the scarap alignment. On a divergent set lower "
                "--core-threshold, or leave the most divergent genomes out.",
                len(orthogroups),
            )
        tables = [core_dir / t for t in TABLES if (core_dir / t).is_file()]
        return CoreGeneResult(
            alignment=alignment,
            tables=tables,
            core_genes=len(orthogroups),
            settings={
                "method": method,
                "seeds": seeds,
                "max_core_genes": max_core,
                "core_prefilter": prefilter,
                # Orthogroups come from MMseqs2 profile searches on proteins,
                # with no fixed identity cut.
                "clustering_on": "protein profiles",
            },
        )


def _orthogroups(genes: Path) -> set[str]:
    """The orthogroups of ``genes.tsv`` (gene, genome, orthogroup per row)."""
    if not genes.is_file():
        return set()
    out: set[str] = set()
    with open(genes, encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[2]:
                out.add(fields[2])
    return out


def _unit(raw: object, flag: str) -> float:
    try:
        value = float(str(raw))
    except ValueError:
        raise UserInputError(f"scarap {flag} must be a number within (0, 1], got {raw!r}") from None
    if not 0.0 < value <= 1.0:
        raise UserInputError(f"scarap {flag} must be within (0, 1], got {raw!r}")
    return value
