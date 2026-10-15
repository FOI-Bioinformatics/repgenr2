"""SCARAP core-gene alignment builder.

Command surface checked against SCARAP 1.0.2 (bioconda). Two calls:

1. ``scarap core FAA_LIST OUT -p P -f F -e SEEDS -t N`` infers the core
   genome directly, without the full pangenome. It draws ``SEEDS`` genomes
   at random (100 by default, with no fixed random seed, so two runs on a
   larger set can differ), builds their pangenome, keeps the families that
   are single-copy in at least ``P`` of the seeds, turns those with
   single-copy presence ``F`` into MMseqs2 profiles and searches every other
   genome against them. It writes ``genes.tsv`` (gene, genome, orthogroup).
   SCARAP's own ``P`` default is 0.90, which would cap how far a lower
   ``--core-threshold`` reaches, so the adapter passes the smaller of 0.90
   and ``F`` unless ``core_prefilter=`` is set.
2. ``scarap concat FAA_LIST genes.tsv OUT -f F -n FFN_LIST -t N`` keeps the
   orthogroups single-copy in at least ``F`` of all the genomes, not only
   the seeds, aligns each with mafft and concatenates them; with ``-n`` it
   writes the nucleotide supermatrix ``supermatrix_nucs.fasta``. A genome
   missing a gene, or holding several copies, is a gap for it.

SCARAP aligns with as many processes as the machine has CPUs, whatever
``-t`` says, and mafft writes its temporary files to ``$TMPDIR``.

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
from ..core.bindirs import host_which
from ..core.containers import run_tool, runs_on_host
from ..core.errors import ToolExecutionError, UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
from .base import CoreGeneBuilder, CoreGeneParams, CoreGeneResult, annotation_sidecars

METHODS = ("H", "FH", "FT", "H-nl", "T-nl", "P", "S", "O-B", "O-D")
# These build the seed pangenome with OrthoFinder (BLAST or DIAMOND search).
# The SCARAP package and image do not ship it, and SCARAP 1.0.2 reads its
# result with a pandas call removed in pandas 2; envs/scarap.yml provides
# both OrthoFinder and pandas older than 2, so they run natively from it.
ORTHOFINDER_METHODS = ("O-B", "O-D")
SUPERMATRIX = "supermatrix_nucs.fasta"
TABLES = ("genes.tsv",)
# Below this a core-gene tree rests on very little.
_SMALL_CORE = 100
# scarap core's own --core-prefilter default.
_SCARAP_PREFILTER = 0.90


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
            else min(_SCARAP_PREFILTER, params.core_threshold)
        )
        max_core = parse_extra_int(params.extra, "max_core_genes", 0)
        seeds = parse_extra_int(params.extra, "seeds", 100)
        if max_core < 0 or seeds < 1:
            raise UserInputError(
                "scarap max_core_genes must be 0 or more and seeds 1 or more, "
                f"got {max_core} and {seeds}"
            )
        if method in ORTHOFINDER_METHODS:
            _check_orthofinder_method(self.capabilities, method)
        faas, ffns = annotation_sidecars(gffs, "scarap")
        _refuse_shared_gene_ids(faas)

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
            "-f", f"{params.core_threshold:g}", "-p", f"{prefilter:g}", "-d", method,
            "-m", str(max_core), "-e", str(seeds), "-t", str(params.threads),
        ]  # fmt: skip
        try:
            run_tool(
                self.capabilities, core, logger=logger, cwd=out_dir,
                log_prefix="scarap", extra_mounts=mounts,
            )  # fmt: skip
        except ToolExecutionError as exc:
            if method in ORTHOFINDER_METHODS and "has no attribute 'append'" in (exc.output or ""):
                raise WorkdirError(
                    f"scarap method {method} ran OrthoFinder and then failed reading its "
                    "result: SCARAP 1.0.2 uses a pandas call removed in pandas 2. Run it from "
                    "the environment envs/scarap.yml builds, which keeps pandas older than 2."
                ) from exc
            raise

        genes = core_dir / "genes.tsv"
        orthogroups = _core_orthogroups(genes, params.core_threshold, max_core)
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
                # Over all the genomes; the core step judged only the seeds.
                "-f", f"{params.core_threshold:g}", "-m", str(max_core),
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


def _check_orthofinder_method(caps: ToolCapabilities, method: str) -> None:
    """Methods O-B and O-D need OrthoFinder (and BLAST for O-B) beside scarap."""
    if not runs_on_host(caps):
        raise UserInputError(
            f"scarap method {method} needs OrthoFinder, which the SCARAP image does not "
            "include. Run natively from the environment envs/scarap.yml builds (with "
            "--bin-dir scarap=<env>/bin), or choose another method."
        )
    needed = ["orthofinder"] + (["blastp"] if method == "O-B" else ["diamond"])
    missing = [b for b in needed if host_which(b, "scarap") is None]
    if missing:
        raise UserInputError(
            f"scarap method {method} needs {', '.join(missing)} beside scarap, and it is not "
            "on scarap's PATH. The environment envs/scarap.yml builds has it."
        )


def _core_orthogroups(genes: Path, threshold: float, max_core: int) -> list[str]:
    """The orthogroups ``scarap concat -f threshold -m max_core`` keeps.

    Mirrors SCARAP's ``determine_corefams``: an orthogroup is kept when the
    genomes holding it in exactly one copy are at least ``threshold`` of the
    genomes in the table; with ``max_core`` the most frequent are kept.
    """
    if not genes.is_file():
        return []
    copies: dict[tuple[str, str], int] = {}
    genomes: set[str] = set()
    with open(genes, encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[2]:
                continue
            genomes.add(fields[1])
            key = (fields[2], fields[1])
            copies[key] = copies.get(key, 0) + 1
    if not genomes:
        return []
    single: dict[str, int] = {}
    for (group, _genome), n in copies.items():
        single.setdefault(group, 0)
        if n == 1:
            single[group] += 1
    kept = [(g, c / len(genomes)) for g, c in single.items() if c / len(genomes) >= threshold]
    kept.sort(key=lambda item: (-item[1], item[0]))
    if max_core > 0:
        kept = kept[:max_core]
    return [g for g, _ in kept]


def _refuse_shared_gene_ids(faas: Sequence[Path]) -> None:
    """SCARAP joins hits to genomes on the gene identifier; a shared one would
    assign genes to the wrong genome. Annotators that derive locus tags from
    the sequence give two identical assemblies the same identifiers."""
    owner: dict[str, str] = {}
    for faa in faas:
        with open(faa, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                gene = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                other = owner.setdefault(gene, faa.stem)
                if other != faa.stem:
                    raise WorkdirError(
                        f"scarap needs gene identifiers unique across genomes, but {gene} is in "
                        f"both {other} and {faa.stem}; the two are probably the same assembly. "
                        "Leave one out, or use another core-gene builder."
                    )


def _unit(raw: object, flag: str) -> float:
    try:
        value = float(str(raw))
    except ValueError:
        raise UserInputError(f"scarap {flag} must be a number within (0, 1], got {raw!r}") from None
    if not 0.0 < value <= 1.0:
        raise UserInputError(f"scarap {flag} must be within (0, 1], got {raw!r}")
    return value
