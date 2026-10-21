"""sourmash tree builder (alignment-free; k-mer distance + neighbor-joining).

The distance is ``1 - ANI`` by default (``metric=ani``), on the scale of
skani and mashtree, or ``1 - Jaccard`` with ``metric=jaccard``. The ANI is
estimated from the k-mer Jaccard of ``sourmash compare`` as
``(2j / (1 + j)) ** (1 / k)``, the quantity ``sourmash compare --ani``
reports (within 1e-4 on genomes of bacterial size). It is computed here
rather than with ``--ani`` because sourmash 4.9.4 reports an ANI of 0 for
small sketches, identical genomes included (measured on 10 kb and 30 kb
genomes at k=31, scaled=1000), which would place identical genomes at the
largest distance. A pair without a shared hash has ANI 0 and distance 1.0.

With ``TreeParams.sketches`` the genomes that have a workdir sketch are not
sketched again: ``sourmash compare -k`` selects the signature at the
requested k-mer size from each ``.sig.zip``. The joining step runs rapidNJ
when it is installed and the pure-Python code otherwise (``nj=`` extra;
see :mod:`repgenr.tree.nj`).
"""

from __future__ import annotations

import csv
import logging
import os
from collections.abc import Mapping, Sequence
from pathlib import Path

from ..core.containers import run_tool
from ..core.contracts import atomic_replace, record_name
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import parse_extra_int, preflight, tool_available
from ..core.process import write_fofn
from ..core.sourmash import sourmash_capabilities
from ..tree.nj import NJ_AUTO, NJ_RAPIDNJ, RAPIDNJ, choose_backend, nj_tree
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list

METRIC_ANI = "ani"
METRIC_JACCARD = "jaccard"
METRICS = (METRIC_ANI, METRIC_JACCARD)


class SourmashBuilder(TreeBuilder):
    capabilities = sourmash_capabilities(
        reads_gzip=True,
        # metric=ani|jaccard: the distance is 1 - ANI (default) or 1 - Jaccard.
        accepted_extras=frozenset({"ksize", "scaled", "nj", "metric"}),
        # A distance-based NJ tree has no bootstrap step.
        ignored_params=frozenset({"bootstrap"}),
        default_params={"ksize": 31, "scaled": 1000, "nj": NJ_AUTO, "metric": METRIC_ANI},
        recommended_max_genomes=2000,
    )
    input_kind = InputKind.GENOMES

    def preflight(self) -> dict[str, str]:
        """sourmash, plus rapidnj when it is installed: with ``nj=auto`` the
        joining then runs there, so the record names it."""
        versions = super().preflight()
        if tool_available(RAPIDNJ):
            versions.update(preflight(RAPIDNJ))
        return versions

    def sketch_request(self, extra: Mapping[str, object]) -> tuple[int, int]:
        defaults = self.capabilities.default_params
        return (
            parse_extra_int(extra, "ksize", defaults["ksize"]),
            parse_extra_int(extra, "scaled", defaults["scaled"]),
        )

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        genomes = as_genome_list(msa_or_genomes)
        metric = _metric(params.extra, self.capabilities.default_params["metric"])
        # Settle the joining step first: a bad nj value, or more leaves than
        # the Python code can join without rapidnj, must not cost a sketch
        # and a dense compare before failing.
        backend = choose_backend(
            str(params.extra.get("nj", self.capabilities.default_params["nj"])),
            len(genomes),
            logger,
        )
        if backend == NJ_RAPIDNJ:
            preflight(RAPIDNJ)  # the binary or the image, before any work
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        scaled = parse_extra_int(params.extra, "scaled", self.capabilities.default_params["scaled"])

        given = {
            g: Path(params.sketches[g]) for g in genomes if params.sketches and g in params.sketches
        }
        rest = [g for g in genomes if g not in given]
        sigs = [given[g] for g in genomes if g in given]
        if rest:
            sigs += self._sketch(rest, out_dir, ksize, scaled, logger)

        matrix_csv = out_dir / "compare.csv"
        compare_fofn = write_fofn(sigs, out_dir / "signatures.fofn")
        run_tool(
            self.capabilities,
            [
                "sourmash",
                "compare",
                "-k",
                str(ksize),
                "--csv",
                matrix_csv,
                "--from-file",
                compare_fofn,
                "--processes",
                str(params.threads),
            ],
            logger=logger,
            log_prefix="sourmash",
            # The given sketches live outside out_dir; their paths are inside the fofn.
            extra_mounts=sorted({os.path.dirname(os.path.abspath(p)) for p in given.values()}),
        )

        labels, similarity = _read_csv(matrix_csv)
        named = [_label_to_genome(label, genomes) for label in labels]
        # sourmash compare does not keep its input order; join in name order so
        # the same genomes give the same tree whatever order the matrix came in.
        order = sorted(range(len(labels)), key=lambda i: (named[i], i))
        clean_labels = [named[i] for i in order]
        if metric == METRIC_ANI:
            similarity = [[jaccard_to_ani(j, ksize) for j in row] for row in similarity]
        dist = [[1.0 - similarity[i][j] for j in order] for i in order]
        disjoint = sum(
            1 for a in range(len(dist)) for b in range(a + 1, len(dist)) if dist[a][b] >= 1.0
        )
        if disjoint:
            logger.warning(
                "sourmash: %d of %d pairs share no hash at k=%d, scaled=%d; their distance "
                "is set to 1.0, so their placement in the tree is not informed by the data. "
                "A lower scaled= or ksize=, or a tree builder for divergent sets, is the "
                "remedy.",
                disjoint,
                len(dist) * (len(dist) - 1) // 2,
                ksize,
                scaled,
            )
        logger.info(
            "sourmash: distance is 1 - %s (metric=%s)",
            "ANI" if metric == METRIC_ANI else "Jaccard",
            metric,
        )
        newick = nj_tree(
            clean_labels, dist, out_dir, logger, backend=backend, threads=params.threads
        )

        tree = out_dir / "tree.nwk"
        with atomic_replace(tree) as fo:
            fo.write(newick + "\n")
        return tree

    def _sketch(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        logger: logging.Logger,
    ) -> list[Path]:
        """Sketch ``genomes`` into ``out_dir/signatures``; one signature file each."""
        sig_dir = out_dir / "signatures"
        sig_dir.mkdir(exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        # Genome paths live inside the fofn (not argv); declare their dirs so the
        # container backend binds them (un-resolved abspaths, matching write_fofn).
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
        run_tool(
            self.capabilities,
            [
                "sourmash",
                "sketch",
                "dna",
                "-p",
                f"k={ksize},scaled={scaled}",
                "--from-file",
                fofn,
                "--outdir",
                sig_dir,
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=genome_dirs,
        )
        # Skip macOS AppleDouble companions ("._*") that appear on exFAT/NTFS volumes.
        sigs = [
            p
            for p in (sorted(sig_dir.glob("*.sig")) + sorted(sig_dir.glob("*.sig.gz")))
            if not p.name.startswith("._")
        ]
        if not sigs:
            raise WorkdirError("sourmash produced no signatures")
        return sigs


def _metric(extra: Mapping[str, object], default: str) -> str:
    value = str(extra.get("metric", default)).strip().lower()
    if value not in METRICS:
        raise UserInputError(
            f"sourmash metric= must be one of {', '.join(METRICS)}, got {extra.get('metric')!r}"
        )
    return value


def jaccard_to_ani(jaccard: float, ksize: int) -> float:
    """The ANI estimate of a k-mer Jaccard similarity: ``(2j / (1 + j)) **
    (1 / k)``, the containment equivalent of Jaccard to the k-th root, as
    ``sourmash compare --ani`` reports it; 0 stays 0."""
    if jaccard <= 0.0:
        return 0.0
    return float((2.0 * jaccard / (1.0 + jaccard)) ** (1.0 / ksize))


def _read_csv(path: Path) -> tuple[list[str], list[list[float]]]:
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo)
        labels = next(reader)
        matrix = [[float(x) for x in row] for row in reader]
    return labels, matrix


def _label_to_genome(label: str, genomes: Sequence[Path]) -> str:
    """The record name of the genome sourmash labelled ``label`` (its input path)."""
    name = record_name(label)
    if name in {record_name(g) for g in genomes}:
        return name
    return Path(label).name
