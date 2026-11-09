"""Minigraph-Cactus (cactus-pangenome) aligner.

Builds a pangenome from same-species genomes. We run ``cactus-pangenome`` with a
generated seqFile, then normalize its HAL output to a reference-anchored
MSA-FASTA via HAL -> MAF -> FASTA. Cactus is resource heavy and drives its own
Toil job store; the Nextflow process that wraps this stage requests a high
resource label.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..converters.hal_to_maf import hal_to_maf
from ..converters.maf_to_fasta import maf_to_fasta
from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import MSA_FASTA, glob_visible, record_name, rename_fasta_records
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from .base import Aligner, AlignParams, AlignResult


class CactusAligner(Aligner):
    capabilities = ToolCapabilities(
        name="cactus",
        # cactus decompresses .gz inputs itself (verified on v2.9.3).
        reads_gzip=True,
        container="quay.io/comparative-genomics-toolkit/cactus:v2.9.3",
        required_binaries=(
            # The pinned image is 2.9; 2.5 is the oldest Minigraph-Cactus the
            # cactus-pangenome entry point and --maxCores were verified on.
            BinarySpec("cactus-pangenome", version_args=("--version",), min_version="2.5"),
            BinarySpec("hal2maf", version_args=()),
        ),
        recommended_max_genomes=2000,  # Toil manages its own parallelism
    )

    def align(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: AlignParams,
        logger: logging.Logger,
    ) -> AlignResult:
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        out_dir.mkdir(parents=True, exist_ok=True)

        # Use absolute (not symlink-resolved) paths so they match how the
        # container backend bind-mounts inputs (macOS firmlinks otherwise
        # diverge between /Users and /System/Volumes/Data).
        genome_paths = [os.path.abspath(g) for g in genomes]
        seqfile = out_dir / "seqfile.txt"
        with open(seqfile, "w", encoding="utf-8") as fo:
            for g, p in zip(genomes, genome_paths, strict=True):
                fo.write(f"{_sample_name(g)}\t{p}\n")

        job_store = out_dir / "jobstore"
        results = out_dir / "cactus_out"
        ref_name = _sample_name(reference)
        # Genome paths live inside seqfile.txt, not in argv, so the backend
        # cannot infer their mounts; declare their directories explicitly.
        genome_dirs = sorted({os.path.dirname(p) for p in genome_paths})
        run_tool(
            self.capabilities,
            [
                "cactus-pangenome",
                job_store,
                seqfile,
                "--outDir",
                results,
                "--outName",
                "pangenome",
                "--reference",
                ref_name,
                # Toil otherwise schedules on every host core.
                "--maxCores",
                str(params.threads),
            ],
            logger=logger,
            log_prefix="cactus",
            # Toil writes .toil/ under HOME, which the container backend points
            # at the working directory; keep that inside the alignment dir.
            cwd=out_dir,
            extra_mounts=genome_dirs,
        )

        hal = _find_hal(results)
        if hal is None:
            raise WorkdirError(f"cactus-pangenome produced no HAL under {results}")

        maf = out_dir / "pangenome.maf"
        hal_to_maf(hal, ref_name, maf, logger, caps=self.capabilities)
        msa = out_dir / MSA_FASTA
        raw_msa = out_dir / "cactus_samples.fasta"
        # Drop the Minigraph-Cactus backbone pseudo-genome so it is not a taxon.
        maf_to_fasta(maf, ref_name, raw_msa, exclude={"_MINIGRAPH_"})
        # Back from cactus sample names to record names: a versioned accession
        # (GCF_000001.1) would otherwise reach the tree builder as GCF_000001_1,
        # and neither IQ-TREE's -o nor tree2tax would find the outgroup.
        rename_fasta_records(raw_msa, msa, {_sample_name(g): record_name(g) for g in genomes})
        return AlignResult(msa_fasta=msa, native_format=hal)


def _sample_name(genome: Path) -> str:
    """Cactus sample names must avoid '.' (reserved for haplotype suffixes)."""
    return record_name(genome).replace(".", "_")


def _find_hal(results: Path) -> Path | None:
    candidates = glob_visible(results, "*.hal", recursive=True)
    if not candidates:
        return None
    # Prefer the combined pangenome HAL (``*.full.hal``) over per-chromosome HALs
    # under ``chrom-alignments/`` -- the latter contain only the genomes mapped
    # to that chromosome, so projecting one drops genomes from the MSA.
    full = [c for c in candidates if c.name.endswith(".full.hal")]
    if full:
        return full[0]
    toplevel = [c for c in candidates if c.parent == results]
    return (toplevel or candidates)[0]
