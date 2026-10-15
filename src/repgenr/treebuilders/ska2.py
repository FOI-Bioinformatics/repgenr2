"""ska2 tree builder (reference-free; split k-mer SNP distances + neighbour joining).

``ska build`` indexes every genome's split k-mers and ``ska distance``
counts, for each pair, the SNPs among the split k-mers both samples hold
(``Distance``) and the fraction of mismatching k-mers (``Mismatches
(proportion)``). No genome is a reference, so no assembly's private errors
bias the distances, which is the case for an outbreak or clonal set. Pairs
that share few split k-mers have their SNPs undercounted, so the builder
warns on a set spanning several species, as the typer does. The joining step
is :func:`repgenr.tree.nj.nj_tree`. Command surface checked against ska
0.5.1 (bioconda ``ska2``).
"""

from __future__ import annotations

import csv
import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import atomic_replace, record_name, taxonomic_spread
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int, preflight
from ..tree.nj import NJ_AUTO, NJ_RAPIDNJ, RAPIDNJ, choose_backend, nj_tree
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list

METRIC_SNPS = "snps"
METRIC_MISMATCH = "mismatch"
METRICS = (METRIC_SNPS, METRIC_MISMATCH)


def warn_ska2_divergence(
    genomes: Sequence[Path], logger: logging.Logger, *, exclude: Sequence[str] = ()
) -> None:
    """Say when split k-mers are being asked to span more than one species.
    ``exclude`` names records left out of the count (the outgroup, which is
    expected to be another species)."""
    n_genera, n_species = taxonomic_spread(
        [g for g in genomes if record_name(g) not in set(exclude)]
    )
    if n_species > 1:
        logger.warning(
            "ska2 on a set spanning %d species (%d genera): split k-mers shared across the "
            "samples shrink with divergence, so the alignment or the distances cover little "
            "of the genome (708 variable sites on 60 Wolbachia species representatives). "
            "ska2 suits clonal and outbreak sets; for a genus-level set use a mapping typer "
            "within one species, or an alignment-free builder (skani, mashtree, sourmash).",
            n_species,
            n_genera,
        )


class Ska2Builder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="ska2",
        reads_gzip=True,
        container="quay.io/biocontainers/ska2:0.5.1--h4349ce8_0",
        conda=("bioconda::ska2",),
        required_binaries=(BinarySpec("ska", version_args=("--version",), min_version="0.3"),),
        # As the typer: the .skf is read into memory for the distance pass.
        recommended_max_genomes=5000,
        # Clonal and outbreak sets only; `auto` must not reach for it on a
        # genus-level set when fasttree is absent.
        auto_eligible=False,
        # min_freq: ska distance's -m (0 by default there, unlike align's 0.9);
        # metric: snps (SNP count, the default) or mismatch (k-mer mismatch fraction).
        default_params={"ksize": 31, "min_freq": 0, "metric": METRIC_SNPS, "nj": NJ_AUTO},
        accepted_extras=frozenset({"ksize", "min_freq", "metric", "nj"}),
        # A distance-based NJ tree has no bootstrap step.
        ignored_params=frozenset({"bootstrap"}),
    )
    input_kind = InputKind.GENOMES

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        genomes = as_genome_list(msa_or_genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        defaults = self.capabilities.default_params
        ksize = parse_extra_int(params.extra, "ksize", defaults["ksize"])
        min_freq = str(params.extra.get("min_freq", defaults["min_freq"]))
        metric = str(params.extra.get("metric", defaults["metric"]))
        if metric not in METRICS:
            raise UserInputError(f"ska2 metric must be one of {', '.join(METRICS)}, not {metric!r}")
        # Settle the joining step before the split k-mer build (repgenr.tree.nj).
        backend = choose_backend(str(params.extra.get("nj", defaults["nj"])), len(genomes), logger)
        if backend == NJ_RAPIDNJ:
            preflight(RAPIDNJ)
        warn_ska2_divergence(genomes, logger, exclude=(params.outgroup,) if params.outgroup else ())

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
        distances = out_dir / "ska_distances.tsv"
        run_tool(
            self.capabilities,
            [
                "ska",
                "distance",
                "--threads",
                str(params.threads),
                "-m",
                min_freq,
                "-o",
                distances,
                skf,
            ],
            logger=logger,
            cwd=out_dir,
            log_prefix="ska-distance",
        )
        names = sorted(record_name(g) for g in genomes)
        dist = read_distances(distances, names, metric)
        logger.info(
            "ska2 distances: %s over the split k-mers each pair shares (k=%d, min_freq %s)",
            "SNP counts" if metric == METRIC_SNPS else "k-mer mismatch fractions",
            ksize,
            min_freq,
        )
        newick = nj_tree(names, dist, out_dir, logger, backend=backend, threads=params.threads)
        tree = out_dir / "tree.nwk"
        with atomic_replace(tree) as fo:
            fo.write(newick + "\n")
        return tree


def read_distances(path: Path, names: Sequence[str], metric: str) -> list[list[float]]:
    """The symmetric matrix over ``names`` from ``ska distance`` output (tab
    separated: Sample1, Sample2, Distance, Mismatches (proportion), Match
    count, Mismatch count). A pair with no matching split k-mer is written
    by ska as distance 0.00 and mismatch 1.0 (checked with ska 0.5.1): zero
    SNPs because nothing was compared, so such a pair is refused rather than
    joined as identical."""
    column = 2 if metric == METRIC_SNPS else 3
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dist = [[0.0] * n for _ in range(n)]
    seen: set[tuple[int, int]] = set()
    unmatched: list[str] = []
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) <= 4 or header[0] != "Sample1":
            raise WorkdirError(f"{path.name}: not a ska distance table")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            a, b = index.get(row[0]), index.get(row[1])
            if a is None or b is None:
                raise WorkdirError(
                    f"{path.name}: sample {row[0]!r} or {row[1]!r} is not an input genome"
                )
            try:
                value = float(row[column])
                matches = int(float(row[4]))
            except (ValueError, IndexError) as exc:
                raise WorkdirError(f"{path.name}: unreadable row for {row[0]} {row[1]}") from exc
            if matches == 0:
                unmatched.append(f"{row[0]} and {row[1]}")
            dist[a][b] = dist[b][a] = max(value, 0.0)
            seen.add((min(a, b), max(a, b)))
    if unmatched:
        shown = "; ".join(unmatched[:5]) + (" ..." if len(unmatched) > 5 else "")
        raise WorkdirError(
            f"{path.name}: {len(unmatched)} pair(s) share no split k-mer at this k ({shown}); "
            "ska writes them as distance 0, which would join them as identical. Leave such "
            "genomes out (--no-outgroup for a distant outgroup), lower --aligner-arg ksize, "
            "or use an alignment-free builder (skani, mashtree)."
        )
    missing = n * (n - 1) // 2 - len(seen)
    if missing:
        raise WorkdirError(
            f"{path.name}: {missing} of {n * (n - 1) // 2} pairs have no row; the table "
            "is truncated"
        )
    return dist
