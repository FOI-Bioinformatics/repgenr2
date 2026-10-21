"""OrthoFinder core-gene alignment builder.

Command surface checked against OrthoFinder 3.1.5 (bioconda). OrthoFinder
compares every proteome against every other (DIAMOND by default), clusters
the hits into orthogroups with MCL and, with ``-og``, stops there; the gene
trees and species tree it would build next are not needed for a
concatenated alignment, and they are the slow part.

1. ``orthofinder -f FAA_DIR -o RES -og -X -S SEARCH -t N -a N`` on one
   ``<record name>.faa`` per genome. OrthoFinder names each genome after the
   file stem, which is the record name, and each gene after the first word
   of its FASTA header, which is how the ``.faa`` and ``.ffn`` files name it
   (``-X`` affects only the gene-tree outputs ``-og`` skips).
2. The adapter selects the orthogroups present in a single copy in at least
   ``--core-threshold`` of the genomes from ``Orthogroups.tsv``, aligns the
   proteins of each with mafft (``--auto``), and back-translates the
   alignment to codons from the ``.ffn`` file (``codon.py``). It hands over
   one DNA alignment per orthogroup, which the stage ranks, caps and
   concatenates. A genome missing a gene, or holding several copies, is
   absent from it; so is a gene whose nucleotide length does not match its
   protein, which is counted in the log.

With ``search=mmseqs`` OrthoFinder writes to fixed names under ``/tmp``, so
two such runs on one host can collide. A path with a space is refused
before the run: OrthoFinder refuses it and builds its search commands as
shell strings.

The search is all against all, so its cost grows roughly with the square
of the genome count: about 6 minutes for 24 Wolbachia genomes on 8 threads
here, against 25 seconds for three.
"""

from __future__ import annotations

import csv
import logging
import shutil
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities
from .base import CoreGeneBuilder, CoreGeneParams, CoreGeneResult, annotation_sidecars
from .codon import align_orthogroups

SEARCHES = ("diamond", "diamond_ultra_sens", "blastp", "mmseqs")
SELECTED = "core_orthogroups.txt"
# Below this a core-gene tree rests on very little.
_SMALL_CORE = 100
# OrthoFinder's own MCL inflation default, recorded when none is given.
DEFAULT_INFLATION = 1.2


class OrthoFinderBuilder(CoreGeneBuilder):
    capabilities = ToolCapabilities(
        name="orthofinder",
        container="quay.io/biocontainers/orthofinder:3.1.5--hdfd78af_0",
        conda=("bioconda::orthofinder",),
        required_binaries=(
            # OrthoFinder colours its version; the parser strips the codes.
            BinarySpec("orthofinder", version_args=("-v",), min_version="3.0"),
            BinarySpec("diamond", version_args=("version",), min_version="2.0"),
            BinarySpec("mafft", version_args=("--version",), min_version="7.0"),
        ),
        # Declared, not measured past 24 genomes: the all-against-all search
        # grows roughly with the square of the genome count.
        recommended_max_genomes=200,
        auto_eligible=False,
        accepted_extras=frozenset({"search", "inflation"}),
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
        search = str(params.extra.get("search", "diamond"))
        if search not in SEARCHES:
            raise UserInputError(
                f"orthofinder search must be one of {', '.join(SEARCHES)}, not {search!r}"
            )
        inflation = _inflation(params.extra.get("inflation"))
        if " " in str(out_dir.resolve()):
            raise UserInputError(
                f"orthofinder cannot run under a path with a space ({out_dir.resolve()}); "
                "move the workdir or output directory."
            )
        faas, ffns = annotation_sidecars(gffs, "orthofinder")
        names = [g.stem for g in gffs]

        out_dir.mkdir(parents=True, exist_ok=True)
        # OrthoFinder reads every FASTA file in one directory and names the
        # genome after the file; copies keep the run inside out_dir.
        faa_dir = out_dir / "proteomes"
        faa_dir.mkdir(exist_ok=True)
        for name, faa in zip(names, faas, strict=True):
            shutil.copyfile(faa, faa_dir / f"{name}.faa")
        results = out_dir / "orthofinder"
        cmd: list[str | Path] = [
            "orthofinder", "-f", faa_dir, "-o", results, "-og", "-X",
            "-S", search, "-t", str(params.threads), "-a", str(params.threads),
        ]  # fmt: skip
        if inflation is not None:
            cmd += ["-I", f"{inflation:g}"]
        run_tool(self.capabilities, cmd, logger=logger, cwd=out_dir, log_prefix="orthofinder")

        table = _orthogroups_table(results)
        selected = _single_copy_core(table, names, params.core_threshold)
        if not selected:
            raise WorkdirError(
                f"orthofinder found no orthogroup present in a single copy in "
                f"{params.core_threshold:g} of the genomes (--core-threshold). Lower it, or "
                "leave the most divergent genomes out."
            )
        (out_dir / SELECTED).write_text("".join(f"{og}\n" for og, _ in selected), "utf-8")

        gene_alignments, mismatched = align_orthogroups(
            selected,
            dict(zip(names, faas, strict=True)),
            dict(zip(names, ffns, strict=True)),
            out_dir,
            self.capabilities,
            params.threads,
            logger,
            "orthofinder",
        )
        if mismatched:
            logger.info(
                "%d gene(s) were left out because their nucleotide length does not match "
                "the protein (partial or frame-shifted genes)",
                mismatched,
            )
        if len(selected) < _SMALL_CORE:
            logger.warning(
                "Only %d core genes in the orthofinder alignment. On a divergent set lower "
                "--core-threshold, or leave the most divergent genomes out.",
                len(selected),
            )
        tables = [table.parent / "Orthogroups.GeneCount.tsv", out_dir / SELECTED]
        return CoreGeneResult(
            tables=[t for t in tables if t.is_file()],
            core_genes=len(selected),
            settings={
                "search": search,
                "inflation": inflation if inflation is not None else DEFAULT_INFLATION,
                "clustering_on": "protein, all against all",
            },
            gene_alignments=gene_alignments,
        )


def _inflation(raw: object) -> float | None:
    if raw is None:
        return None
    try:
        value = float(str(raw))
    except ValueError:
        raise UserInputError(f"orthofinder inflation must be a number, got {raw!r}") from None
    if value <= 1.0:
        raise UserInputError(f"orthofinder inflation must be above 1, got {raw!r}")
    return value


def _orthogroups_table(results: Path) -> Path:
    found = sorted(results.glob("Results_*/Orthogroups/Orthogroups.tsv"))
    if not found:
        raise WorkdirError(
            f"orthofinder wrote no Orthogroups.tsv under {results}; its log there says why."
        )
    return found[-1]


def _single_copy_core(
    table: Path, names: Sequence[str], threshold: float
) -> list[tuple[str, dict[str, str]]]:
    """Orthogroups single-copy in at least ``threshold`` of ``names``.

    Returns each with its single-copy member per genome, in table order.
    """
    out: list[tuple[str, dict[str, str]]] = []
    with open(table, encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        columns = header[1:]
        unknown = sorted(set(columns) - set(names))
        if unknown:
            raise WorkdirError(
                f"orthofinder named genomes that are none of the inputs ({', '.join(unknown[:3])})"
            )
        for row in reader:
            members: dict[str, str] = {}
            for genome, cell in zip(columns, row[1:], strict=False):
                genes = [g for g in cell.split(", ") if g]
                if len(genes) == 1:
                    members[genome] = genes[0]
            if members and len(members) / len(names) >= threshold:
                out.append((row[0], members))
    return out
