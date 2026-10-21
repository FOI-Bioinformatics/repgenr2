"""ska2 (split k-mer analysis) reference-free SNP typer.

Every genome is compared through split k-mers, so no genome is privileged as
the reference and reference-private errors do not bias distances. The output
is a variable-site alignment; there is no positional whole-genome alignment,
so recombination masking is not available with this typer (the snptype stage
refuses ``--mask`` for it with a message naming the alternatives).

Command surface checked against ska 0.5.1 (bioconda ``ska2``): ``ska build -f
<name TAB path list> -k K --threads N -o <prefix>`` writes ``<prefix>.skf``;
``ska align --min-freq F --filter no-filter --threads N -o all_sites.fasta
<prefix>.skf`` writes the middle base of every shared split k-mer, constant
columns included; the snptype stage counts its constant A, C, G and T
columns for the tree builder's ascertainment-bias correction (IQ-TREE
``-fconst``, RAxML-NG ``+ASC_STAM``) and then removes it. The second call,
``ska align --min-freq F --filter no-ambig-or-const --threads N -o <out>
<prefix>.skf``, writes the variable-site alignment. With ``lo=true`` the
alignment comes from ``ska lo <prefix>.skf <out>`` instead
(``<out>_snps.fas``), which recovers variants in more diverged samples
through local assembly of the split k-mer graph (Horsfield et al. 2025);
``lo_missing`` is its maximum fraction of missing data per site (0.1).
ska lo writes variants only, so that route has no constant-site counts.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import record_name
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_bool, parse_extra_int
from .base import SnpParams, SnpResult, SnpTyper


def _has_columns(fasta: Path) -> bool:
    """True when some record of ``fasta`` has at least one column."""
    with open(fasta, encoding="utf-8") as fo:
        for line in fo:
            if line.startswith(">") or not line.strip():
                continue
            return True
    return False


class Ska2Typer(SnpTyper):
    capabilities = ToolCapabilities(
        name="ska2",
        reads_gzip=True,
        container="quay.io/biocontainers/ska2:0.5.1--h4349ce8_0",
        conda=("bioconda::ska2",),
        required_binaries=(BinarySpec("ska", version_args=("--version",), min_version="0.3"),),
        # Split k-mer building is linear in the input; the alignment step reads
        # the whole .skf into memory, which is the practical ceiling.
        recommended_max_genomes=5000,
        default_params={"ksize": 31, "min_freq": 0.9, "lo": False, "lo_missing": 0.1},
        accepted_extras=frozenset({"ksize", "min_freq", "lo", "lo_missing"}),
    )
    requires_reference = False
    # Split k-mers have no reference coordinates: no whole-genome alignment to
    # mask, only an all-sites alignment for constant-site counts.
    produces_full_alignment = False
    # Shared split k-mers shrink with divergence: the stage warns on a set
    # spanning more than one species.
    divergence_sensitive = True

    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        # ``reference`` is accepted for interface parity and deliberately unused:
        # every genome, including any would-be reference, is an ordinary sample.
        genomes = list(genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        min_freq = str(params.extra.get("min_freq", self.capabilities.default_params["min_freq"]))

        # ska's input list is two tab-separated columns: sample name, FASTA path.
        # Paths are absolute but not symlink-resolved, matching write_fofn, so
        # the container backend can bind the same directories it sees here.
        listing = out_dir / "genomes.tsv"
        listing.write_text(
            "".join(f"{record_name(g)}\t{os.path.abspath(g)}\n" for g in genomes),
            encoding="utf-8",
        )
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})

        prefix = out_dir / "split_kmers"
        run_tool(
            self.capabilities,
            [
                "ska",
                "build",
                "-f",
                listing,
                "-k",
                str(ksize),
                "--threads",
                str(params.threads),
                "-o",
                prefix,
            ],
            logger=logger,
            cwd=out_dir,
            log_prefix="ska-build",
            extra_mounts=genome_dirs,
        )
        skf = Path(str(prefix) + ".skf")
        if not skf.exists():
            raise WorkdirError("ska build did not produce a split k-mer (.skf) file")

        core = out_dir / "core_snp.fasta"
        if parse_extra_bool(params.extra, "lo", False):
            lo_prefix = out_dir / "lo"
            defaults = self.capabilities.default_params
            missing = str(params.extra.get("lo_missing", defaults["lo_missing"]))
            if "min_freq" in params.extra:
                logger.info(
                    "ska lo does not read min_freq; its site filter is lo_missing (%s)", missing
                )
            run_tool(
                self.capabilities,
                ["ska", "lo", "--threads", str(params.threads), "-m", missing, skf, lo_prefix],
                logger=logger,
                cwd=out_dir,
                log_prefix="ska-lo",
            )
            snps = Path(str(lo_prefix) + "_snps.fas")
            if not snps.exists() or not _has_columns(snps):
                raise WorkdirError(
                    "ska lo produced no variable sites (lo_snps.fas is empty); the genomes may "
                    "be identical, or too divergent for split k-mers at this k."
                )
            snps.replace(core)
            indels = Path(str(lo_prefix) + "_indels.vcf")
            logger.info(
                "ska lo: variants recovered by local assembly; the indels are published as "
                "variants.vcf"
            )
            logger.info(
                "ska lo writes variants only; no constant-site counts, so the tree's branch "
                "lengths are not corrected"
            )
            return SnpResult(
                core_snp_fasta=core,
                vcf=indels if indels.exists() else None,
                masked=False,
                full_alignment=None,
            )
        # Every shared split k-mer middle base, constant columns included: the
        # source of the constant-site counts. Written first so the variable-site
        # alignment, and its refusal below, are as without it.
        all_sites = out_dir / "all_sites.fasta"
        self._align(skf, all_sites, "no-filter", min_freq, params, logger, out_dir)
        self._align(skf, core, "no-ambig-or-const", min_freq, params, logger, out_dir)
        if not core.exists() or core.stat().st_size == 0:
            raise WorkdirError(
                "ska align produced no variable-site alignment: the genomes may be "
                "identical, or too divergent for split k-mers at this k (try a smaller "
                "--tool-arg ksize, or lower min_freq)."
            )
        all_sites_alignment: Path | None = all_sites
        if not all_sites.exists():
            logger.warning(
                "ska align --filter no-filter wrote no all_sites.fasta; no constant-site "
                "counts, so the tree's branch lengths are not corrected"
            )
            all_sites_alignment = None
        elif not _has_columns(all_sites):
            logger.warning(
                "ska align --filter no-filter wrote an all_sites.fasta without columns; no "
                "constant-site counts, so the tree's branch lengths are not corrected"
            )
            all_sites_alignment = None
        return SnpResult(
            core_snp_fasta=core,
            masked=False,
            full_alignment=None,
            all_sites_alignment=all_sites_alignment,
        )

    def _align(
        self,
        skf: Path,
        out: Path,
        filter_: str,
        min_freq: str,
        params: SnpParams,
        logger: logging.Logger,
        cwd: Path,
    ) -> None:
        """``ska align`` of ``skf`` into ``out`` with the given site filter."""
        run_tool(
            self.capabilities,
            [
                "ska",
                "align",
                "--min-freq",
                min_freq,
                "--filter",
                filter_,
                "--threads",
                str(params.threads),
                "-o",
                out,
                skf,
            ],
            logger=logger,
            cwd=cwd,
            log_prefix="ska-align",
        )
