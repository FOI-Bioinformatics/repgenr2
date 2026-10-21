"""sourmash dereplication adapter (k-mer / MinHash ANI clustering).

Sketches each genome and greedily picks representatives: walk genomes ordered by
connectivity, and absorb any genome whose *estimated ANI* to a chosen
representative is above the threshold. Fast and low-memory, useful as a scalable
alternative to alignment-based ANI.

The threshold is on sourmash's ANI estimate, the same scale skder and galah
use — not on raw sketch similarity. At k=31 a genome pair at 99.6 percent ANI
has Jaccard only ~0.8, so thresholding the raw similarity at an ANI-style 0.99
would leave nearly everything a singleton (verified against skani on synthetic
sets; see docs/audit/scaling-audit.md).

Two back-ends compute the pairwise similarities, picked automatically:

* **Sparse (preferred at scale)**: when the ``sourmash_plugin_branchwater``
  plugin is installed, ``sourmash scripts manysketch`` + ``sourmash scripts
  pairwise`` produce only the above-threshold edges (an edge list), never the
  dense N x N matrix. This keeps memory roughly linear in the number of close
  pairs rather than quadratic in the number of genomes, which matters at 10k+.
* **Dense (fallback)**: plain ``sourmash sketch`` + ``sourmash compare`` build
  the full N x N similarity matrix. Used when the plugin is absent (e.g. inside a
  stock BioContainer). For the same threshold both back-ends pick the same
  representatives on well-separated inputs.

``compare`` (used by ``repgenr glance``) always takes the dense path, since it
reports every genome pair; it converts to ANI the same way, so the glance plots
use the scale of the dereplication threshold.

Both back-ends and ``compare`` read the workdir genome sketches when the stage
passes them (``DerepParams.sketches``, the ``sketches`` argument of
``compare``): the signature at the requested k-mer size is selected from each
``.sig.zip`` with ``-k``, and only the genomes without such a file are sketched
here. Signatures of the workdir sketches are named by the genome's record
name; those sketched here by the file name (dense) or ``Path.stem`` (sparse).
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import uuid
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..core.containers import run_tool
from ..core.contracts import record_name
from ..core.errors import MissingBinaryError, ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
from ..core.process import write_fofn
from ..core.sourmash import jaccard_to_ani_matrix as _jaccard_to_ani_matrix
from ..core.sourmash import sourmash_capabilities
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    CompareResult,
    Dereplicator,
    DerepParams,
    DerepResult,
)
from .compare_io import write_dendrogram, write_pairwise_csv

# Above this, the dense N x N float64 matrix is too large to hold in memory
# (~0.2 GB at 5k, ~20 GB at 50k); require the sparse branchwater path instead.
_DENSE_MAX_GENOMES = 5000

# Axis label for the glance plots: sourmash's ANI estimate from k-mer sketches.
_MEASURE = "ANI"


class SourmashDereplicator(Dereplicator):
    capabilities = sourmash_capabilities(
        accepted_extras=frozenset({"ksize", "scaled", "sketch_cache", "dense_fallback"}),
        # k-mer containment has no primary/secondary split and no alignment.
        ignored_params=frozenset({"primary_ani", "aligned_fraction"}),
        default_params={"ksize": 31, "scaled": 1000},
        recommended_max_genomes=None,
        supports_native_scaling=True,
    )

    def sketch_request(self, extra: Mapping[str, object]) -> tuple[int, int]:
        defaults = self.capabilities.default_params
        return (
            parse_extra_int(extra, "ksize", defaults["ksize"]),
            parse_extra_int(extra, "scaled", defaults["scaled"]),
        )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        scaled = parse_extra_int(params.extra, "scaled", self.capabilities.default_params["scaled"])
        sani = params.secondary_ani
        threshold = sani if sani <= 1.0 else sani / 100

        # Signatures depend only on the genome set + ksize/scaled, not the ANI
        # threshold. The --target-reps search passes a shared sketch_cache dir so
        # the (expensive) sketching runs once and every threshold iteration
        # reuses it instead of re-sketching the same genomes.
        cache = params.extra.get("sketch_cache")
        sketch_cache = Path(cache) if cache else None
        sketches = params.sketches

        if _branchwater_available(self.capabilities, logger):
            try:
                clusters, status = self._sparse_dereplicate(
                    genomes,
                    out_dir,
                    ksize,
                    scaled,
                    threshold,
                    params,
                    logger,
                    sketch_cache,
                    sketches,
                )
            except (ToolExecutionError, WorkdirError) as exc:
                # Availability was pre-probed, so this is a genuine tool failure
                # (bad parameters, OOM, broken image) -- surface it. Silently
                # switching to the dense back-end can change which
                # representatives are picked, so the fallback is opt-in.
                if not _dense_fallback_requested(params):
                    raise
                logger.warning(
                    "sourmash branchwater sparse path failed (%s); "
                    "dense fallback requested, retrying with dense compare",
                    exc,
                )
                clusters, status = self._dense_dereplicate(
                    genomes,
                    out_dir,
                    ksize,
                    scaled,
                    threshold,
                    logger,
                    sketch_cache,
                    params.threads,
                    sketches,
                )
        else:
            clusters, status = self._dense_dereplicate(
                genomes,
                out_dir,
                ksize,
                scaled,
                threshold,
                logger,
                sketch_cache,
                params.threads,
                sketches,
            )

        rep_paths = [p for p in genomes if p.name in clusters]
        return DerepResult(
            representatives=sorted(rep_paths),
            clusters=clusters,
            genome_status=status,
        )

    def _dense_dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        threshold: float,
        logger: logging.Logger,
        sketch_cache: Path | None = None,
        threads: int = 1,
        sketches: Mapping[Path, Path] | None = None,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        """Stock sourmash sketch + N x N compare (no plugin needed)."""
        labels, ani, name_by_label = self._dense_ani_matrix(
            genomes,
            out_dir,
            ksize,
            scaled,
            logger,
            sketch_cache,
            threads,
            sketches=sketches,
            too_large_hint="Install the branchwater plugin "
            "(pip install sourmash_plugin_branchwater) for the sparse path, or use a "
            "tool that scales better at this size (e.g. --tool skder).",
        )
        return _greedy_cluster(labels, ani, name_by_label, threshold)

    def compare(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        threads: int,
        logger: logging.Logger,
        *,
        sketches: Mapping[Path, Path] | None = None,
    ) -> CompareResult:
        """All-vs-all ANI estimates for ``repgenr glance`` (dense compare).

        Uses the same sketches (default ksize/scaled) and the same
        Jaccard-to-ANI conversion as the dense dereplication path, so the
        glance plots are on the scale the ``--secondary-ani`` threshold of
        ``dereplicate --tool sourmash`` applies to. Writes a pairwise CSV (one
        row per unordered pair) and an average-linkage dendrogram on 1 - ANI.
        """
        out_dir.mkdir(parents=True, exist_ok=True)
        defaults = self.capabilities.default_params
        labels, ani, name_by_label = self._dense_ani_matrix(
            genomes,
            out_dir,
            int(defaults["ksize"]),
            int(defaults["scaled"]),
            logger,
            None,
            threads,
            sketches=sketches,
            too_large_hint="glance plots every genome pair and needs the full matrix; "
            "compare a subset of the genomes instead.",
        )
        names = [name_by_label[label] for label in labels]
        pairwise = write_pairwise_csv(names, ani, out_dir / "pairwise_ani.csv")
        pdf = out_dir / "dendrogram.pdf"
        leaves = write_dendrogram(names, ani, pdf, measure=_MEASURE)
        # The leaf order in plain text, for checking the figure (kept with --keep-files).
        (out_dir / "dendrogram_leaves.txt").write_text("\n".join(leaves) + "\n", encoding="utf-8")
        return CompareResult(similarity_csv=pairwise, dendrogram=pdf, measure=_MEASURE)

    def _dense_ani_matrix(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        logger: logging.Logger,
        sketch_cache: Path | None,
        threads: int,
        *,
        too_large_hint: str,
        sketches: Mapping[Path, Path] | None = None,
    ) -> tuple[list[str], npt.NDArray[np.float64], dict[str, str]]:
        """Sketch (reusing given or cached signatures) and compare all genomes.

        Returns the column labels and the N x N ANI-estimate matrix, both in
        genome file name order (``sourmash compare`` does not keep its input
        order, and the greedy pick breaks ties by position), and the map from
        label to genome basename.
        """
        if len(genomes) > _DENSE_MAX_GENOMES:
            raise WorkdirError(
                f"sourmash dense compare needs an N x N matrix for {len(genomes)} genomes "
                f"(~{len(genomes) ** 2 * 8 / 1e9:.1f} GB). {too_large_hint}"
            )
        given = {g: Path(sketches[g]) for g in genomes if sketches and g in sketches}
        rest = [g for g in genomes if g not in given]
        sig_by_genome: dict[Path, Path] = dict(given)
        if rest:
            sig_dir = sketch_cache if sketch_cache is not None else (out_dir / "signatures")
            sig_by_genome.update(
                self._sketch_dense(rest, sig_dir, out_dir, ksize, scaled, logger, len(given))
            )

        matrix_csv = out_dir / "compare.csv"
        sig_files = [sig_by_genome[g] for g in sorted(genomes, key=lambda g: g.name)]
        # Pass signatures via --from-file, never on argv (ARG_MAX at scale).
        compare_fofn = write_fofn(sig_files, out_dir / "signatures.fofn")
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
                str(threads),
            ],
            logger=logger,
            log_prefix="sourmash",
            # The given sketches live outside out_dir; their paths are inside the fofn.
            extra_mounts=sorted({os.path.dirname(os.path.abspath(p)) for p in given.values()}),
        )

        # Convert the Jaccard matrix to ANI estimates so the threshold applies
        # on the same scale as skder/galah (raw Jaccard at k=31 is ~0.8 for a
        # 99.6 pct ANI pair). Converting here, rather than via ``compare
        # --ani``, keeps identical-content pairs at 1.0 even for tiny sketches,
        # where sourmash's estimator refuses and reports 0.
        labels, sim = _read_compare_csv(matrix_csv)
        name_by_label = _match_labels_to_genomes(labels, genomes)
        order = sorted(range(len(labels)), key=lambda i: (name_by_label[labels[i]], i))
        labels = [labels[i] for i in order]
        sim = np.asarray(sim, dtype=float)[np.ix_(order, order)]
        ani = _jaccard_to_ani_matrix(sim, ksize)
        return labels, ani, name_by_label

    def _sketch_dense(
        self,
        genomes: Sequence[Path],
        sig_dir: Path,
        out_dir: Path,
        ksize: int,
        scaled: int,
        logger: logging.Logger,
        n_given: int,
    ) -> dict[Path, Path]:
        """Signature file per genome in ``sig_dir``, sketching those not cached there."""
        sig_dir.mkdir(parents=True, exist_ok=True)
        # Reuse is decided per genome, never by counting files: the cache dir may
        # be shared with other chunks (disjoint genome sets), so only this call's
        # genomes may be sketched, matched, and compared.
        sig_by_genome = _find_signatures(sig_dir, genomes)
        missing = [g for g in genomes if g not in sig_by_genome]
        if missing:
            if len(missing) < len(genomes) + n_given:
                logger.info(
                    "Reusing %d sourmash signatures, sketching %d",
                    len(genomes) + n_given - len(missing),
                    len(missing),
                )
            fofn = write_fofn(missing, out_dir / "genomes.fofn")
            # The genome paths live inside the fofn, not in argv, so the container
            # backend cannot infer their mounts; declare their directories (using
            # un-resolved abspaths to match write_fofn and the backend's bind logic).
            genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in missing})
            # one signature file per genome, named by genome basename
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
            sig_by_genome = _find_signatures(sig_dir, genomes)
        else:
            logger.info("Reusing %d cached sourmash signatures", len(genomes))

        unsketched = [g.name for g in genomes if g not in sig_by_genome]
        if unsketched:
            raise WorkdirError(
                f"sourmash produced no signatures under {sig_dir} for: "
                + ", ".join(unsketched[:5])
                + ("..." if len(unsketched) > 5 else "")
            )
        return sig_by_genome

    def _sparse_dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        threshold: float,
        params: DerepParams,
        logger: logging.Logger,
        sketch_cache: Path | None = None,
        sketches: Mapping[Path, Path] | None = None,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        """Branchwater manysketch + pairwise: emit only above-threshold edges.

        ``pairwise -t`` is a *containment* threshold for which pairs to report.
        An ANI of ``threshold`` corresponds to containment ``threshold ** ksize``
        (ANI estimate = containment^(1/k)), and average containment never exceeds
        max containment, so prefiltering at ``threshold ** ksize`` yields a
        superset of the wanted edges; the parser then keeps pairs whose
        containment-derived ANI is >= ``threshold``, matching the dense
        ``compare --ani`` graph.

        Genomes with a given sketch are not sketched: ``pairwise`` reads a
        path list of their ``.sig.zip`` files (plus the zip of the others) and
        selects the signatures at ``ksize`` with ``-k``.
        """
        threads = str(params.threads)
        given = {g: Path(sketches[g]) for g in genomes if sketches and g in sketches}
        rest = [g for g in genomes if g not in given]
        collections: list[Path] = []
        if rest:
            collections.append(
                self._sparse_sketch(rest, out_dir, ksize, scaled, threads, logger, sketch_cache)
            )
        mounts: list[str] = []
        if given:
            if rest:
                logger.info(
                    "Reusing %d genome sketches, sketching %d with sourmash",
                    len(given),
                    len(rest),
                )
            source = write_fofn(
                [*(given[g] for g in sorted(given, key=lambda g: g.name)), *collections],
                out_dir / "signatures.txt",
            )
            mounts = sorted({os.path.dirname(os.path.abspath(p)) for p in given.values()})
        else:
            source = collections[0]

        pairwise_csv = out_dir / "pairwise.csv"
        run_tool(
            self.capabilities,
            [
                "sourmash",
                "scripts",
                "pairwise",
                source,
                "-o",
                pairwise_csv,
                "-t",
                f"{threshold**ksize:g}",
                "-k",
                str(ksize),
                "-c",
                threads,
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=mounts,
        )

        # The workdir sketches name their signature by record name; the
        # manysketch CSV below names it by Path.stem.
        label_of = {g: record_name(g) if g in given else g.stem for g in genomes}
        name_by_label = {label_of[g]: g.name for g in genomes}
        # Iterate in genome-basename order so the greedy tie-break (which member of
        # a mutually-similar group becomes the representative) matches the dense
        # ``compare`` path, which orders its matrix by genome file name.
        labels = [label_of[g] for g in sorted(genomes, key=lambda g: g.name)]
        neighbors = _parse_pairwise_csv(pairwise_csv, threshold, set(labels), ksize)
        return _sparse_greedy_cluster(labels, neighbors, name_by_label)

    def _sparse_sketch(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        ksize: int,
        scaled: int,
        threads: str,
        logger: logging.Logger,
        sketch_cache: Path | None,
    ) -> Path:
        """A zip of the signatures of ``genomes`` (manysketch, or the shared cache)."""
        if sketch_cache is not None:
            # The cache dir may be shared with other chunks (disjoint genome
            # sets), so the zip is keyed by the genome set + sketch params: a
            # different set can never be mistaken for this one, and concurrent
            # writers target distinct temp files promoted atomically.
            sketch_cache.mkdir(parents=True, exist_ok=True)
            digest = _genome_set_digest(genomes, ksize, scaled)
            sigs_zip = sketch_cache / f"signatures-{digest}.zip"
        else:
            sigs_zip = out_dir / "signatures.zip"

        if sigs_zip.exists():
            logger.info("Reusing cached sourmash %s", sigs_zip.name)
        elif sketch_cache is not None and self._collect_from_cache(
            genomes, sketch_cache, sigs_zip, out_dir, ksize, scaled, logger
        ):
            pass
        else:
            # manysketch reads a CSV of (name, genome_filename, protein_filename).
            # The name becomes the signature name, which is what pairwise reports --
            # plain sketch leaves it empty, so the edge list would be unlabelled.
            # Fields are csv-quoted, so a name or path with a comma stays one
            # field (manysketch, sig cat and pairwise all accept quoted fields).
            sketch_csv = out_dir / "manysketch.csv"
            _write_csv(
                sketch_csv,
                [
                    ["name", "genome_filename", "protein_filename"],
                    *([g.stem, os.path.abspath(g), ""] for g in genomes),
                ],
            )

            # Unique temp name beside the final zip (same dir, so the replace is
            # atomic and concurrent writers of the same set cannot collide).
            tmp_zip = sigs_zip.parent / f".{sigs_zip.stem}.{uuid.uuid4().hex}.partial.zip"
            genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
            run_tool(
                self.capabilities,
                [
                    "sourmash",
                    "scripts",
                    "manysketch",
                    sketch_csv,
                    "-o",
                    tmp_zip,
                    "-p",
                    f"dna,k={ksize},scaled={scaled}",
                    "-c",
                    threads,
                ],
                logger=logger,
                log_prefix="sourmash",
                extra_mounts=[*genome_dirs, str(sketch_csv), str(tmp_zip.parent)],
            )
            if not tmp_zip.exists():
                raise WorkdirError(f"sourmash manysketch produced no signatures at {tmp_zip}")
            os.replace(tmp_zip, sigs_zip)
            if sketch_cache is not None:
                _write_zip_index(sigs_zip, genomes, ksize, scaled)
        return sigs_zip

    def _collect_from_cache(
        self,
        genomes: Sequence[Path],
        sketch_cache: Path,
        sigs_zip: Path,
        out_dir: Path,
        ksize: int,
        scaled: int,
        logger: logging.Logger,
    ) -> bool:
        """Assemble ``sigs_zip`` from sketched zips in the cache, if they cover the set.

        Under chunked ``--target-reps`` the merge level dereplicates the union
        of the chunk representatives, which changes with the threshold at each
        search step; its signatures are already in the chunk zips. ``sourmash
        sig cat`` with a picklist of the wanted names copies them into a new
        zip, which takes far less time than sketching the genomes again.
        Returns False (the caller sketches) when the cached zips do not cover
        every genome or ``sig cat`` fails.
        """
        sources = _select_cached_zips(
            {g.stem: _file_key(g) for g in genomes}, _read_zip_indexes(sketch_cache, ksize, scaled)
        )
        if sources is None:
            return False
        picklist = out_dir / "picklist.csv"
        # csv quoting keeps a name with a comma one field (sourmash reads the
        # picklist with Python's csv module).
        _write_csv(picklist, [["name"], *([g.stem] for g in sorted(genomes, key=lambda g: g.stem))])
        tmp_zip = sigs_zip.parent / f".{sigs_zip.stem}.{uuid.uuid4().hex}.partial.zip"
        cmd: list[str | Path] = ["sourmash", "sig", "cat", *sources]
        cmd += [
            "--picklist",
            f"{picklist}:name:name",
            "--picklist-require-all",
            "-k",
            str(ksize),
            "-o",
            tmp_zip,
        ]
        try:
            run_tool(
                self.capabilities,
                cmd,
                logger=logger,
                log_prefix="sourmash",
                extra_mounts=[str(sketch_cache), str(picklist)],
            )
        except ToolExecutionError as exc:
            logger.warning(
                "Collecting %d cached sourmash signatures failed (%s); sketching them instead",
                len(genomes),
                exc,
            )
            tmp_zip.unlink(missing_ok=True)
            return False
        if not tmp_zip.exists():
            logger.warning("sourmash sig cat wrote no %s; sketching instead", tmp_zip.name)
            return False
        os.replace(tmp_zip, sigs_zip)
        logger.info(
            "Collected %d cached sourmash signatures from %d sketched zip(s) into %s",
            len(genomes),
            len(sources),
            sigs_zip.name,
        )
        return True


# Signature names are ``Path.stem`` of the genome file (``x.fna.gz`` -> ``x.fna``),
# used alike in the manysketch CSV, the picklist, this index and the pairwise
# labels. They are internal to this adapter and need not follow
# ``strip_fasta_suffix``; genomes are reported under their file names.
#
# Index written beside each zip that manysketch produced in a shared cache:
# which signature names it holds and with which sketch parameters. Zips
# assembled by sig cat get no index, so they are never offered as a source and
# the sources of one collection cannot overlap (sig cat would emit a genome
# twice). A zip without its index (e.g. after an interruption) is not offered.
_INDEX_SUFFIX = ".names.json"


def _write_csv(path: Path, rows: Iterable[Sequence[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fo:
        csv.writer(fo, lineterminator="\n").writerows(rows)


def _file_key(path: Path) -> tuple[int, int]:
    """(size, mtime in ns) of a genome file, following symlinks; (-1, -1) if unreadable.

    Part of every sketch-cache key, so a genome replaced under the same name
    (in a persistent ``--tool-arg sketch_cache=`` directory) is sketched again
    rather than matched to the old signature.
    """
    try:
        st = os.stat(path)
    except OSError:
        return (-1, -1)
    return (st.st_size, st.st_mtime_ns)


def _write_zip_index(sigs_zip: Path, genomes: Sequence[Path], ksize: int, scaled: int) -> None:
    index = sigs_zip.with_name(sigs_zip.stem + _INDEX_SUFFIX)
    tmp = index.with_name(f".{index.name}.{uuid.uuid4().hex}.partial")
    names = {g.stem: list(_file_key(g)) for g in genomes}
    tmp.write_text(
        json.dumps({"ksize": ksize, "scaled": scaled, "names": names}, sort_keys=True),
        encoding="utf-8",
    )
    os.replace(tmp, index)


def _read_zip_indexes(
    sketch_cache: Path, ksize: int, scaled: int
) -> dict[Path, dict[str, tuple[int, int]]]:
    """Sketched zips in ``sketch_cache`` made with these parameters -> name: file key."""
    out: dict[Path, dict[str, tuple[int, int]]] = {}
    for index in sorted(sketch_cache.glob(f"signatures-*{_INDEX_SUFFIX}")):
        zip_path = index.with_name(index.name[: -len(_INDEX_SUFFIX)] + ".zip")
        try:
            data = json.loads(index.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not zip_path.exists() or not isinstance(data, dict):
            continue
        if data.get("ksize") != ksize or data.get("scaled") != scaled:
            continue
        names = data.get("names")
        if not isinstance(names, dict):
            continue
        try:
            out[zip_path] = {str(n): (int(k[0]), int(k[1])) for n, k in names.items()}
        except (TypeError, ValueError, IndexError):
            continue
    return out


def _select_cached_zips(
    wanted: Mapping[str, tuple[int, int]],
    indexes: Mapping[Path, Mapping[str, tuple[int, int]]],
) -> list[Path] | None:
    """Zips that together hold every wanted name, none sharing a wanted name.

    A zip holds a name only when it was sketched from a file of the same size
    and modification time. Greedy, largest overlap first. Returns None when the
    zips cannot cover the set without overlap; the caller then sketches.
    """
    hits = {
        zip_path: {n for n, key in names.items() if wanted.get(n) == key}
        for zip_path, names in indexes.items()
    }
    remaining = set(wanted)
    chosen: list[Path] = []
    for zip_path, hit in sorted(hits.items(), key=lambda kv: (-len(kv[1]), kv[0].name)):
        if not hit or not hit <= remaining:
            continue
        chosen.append(zip_path)
        remaining -= hit
        if not remaining:
            return chosen
    return None


def _find_signatures(sig_dir: Path, genomes: Sequence[Path]) -> dict[Path, Path]:
    """Map each genome to its signature file in ``sig_dir``, by exact name.

    ``sourmash sketch --outdir`` names signatures after the input file
    (``<name>.sig``); some producers use the stem instead, so both are accepted.
    Only exact per-genome matches count -- a shared cache dir may hold
    signatures for other genome sets, which must never satisfy a lookup --
    and only signatures at least as new as the genome file.
    """
    out: dict[Path, Path] = {}
    for g in genomes:
        genome_mtime = _file_key(g)[1]
        for cand in (f"{g.name}.sig", f"{g.name}.sig.gz", f"{g.stem}.sig", f"{g.stem}.sig.gz"):
            p = sig_dir / cand
            # A signature older than its genome file was made from an earlier
            # file of that name (a persistent sketch_cache); sketch again.
            if p.exists() and p.stat().st_mtime_ns >= genome_mtime:
                out[g] = p
                break
    return out


def _genome_set_digest(genomes: Sequence[Path], ksize: int, scaled: int) -> str:
    """Short digest identifying a genome set + sketch params for cache keying.

    Each genome counts by name, size and modification time, so a file
    replaced under the same name gives a different digest.
    """
    entries = sorted(f"{g.name}:{_file_key(g)[0]}:{_file_key(g)[1]}" for g in genomes)
    blob = "\n".join(entries) + f"|k={ksize}|scaled={scaled}"
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:12]


def _dense_fallback_requested(params: DerepParams) -> bool:
    """Opt-in to retrying a failed sparse run with the dense back-end."""
    if params.extra.get("dense_fallback"):
        return True
    return os.environ.get("REPGENR_SOURMASH_DENSE_FALLBACK", "") not in ("", "0")


_BRANCHWATER_CACHE: dict[tuple, bool] = {}


def _branchwater_available(caps: ToolCapabilities, logger: logging.Logger) -> bool:
    """True when the branchwater plugin (``sourmash scripts pairwise``) is usable.

    Probed through ``run_tool`` so the answer reflects where the tool actually
    runs: native (plugin pip-installed) vs a container image that may not bundle
    it. Any failure (missing subcommand, missing image, missing binary) means the
    dense fallback is used instead.

    The probe is a real subprocess/container start, so memoize it per (tool,
    container config) for the process lifetime -- chunked dereplication calls this
    once per chunk otherwise, paying a container cold-start each time.
    """
    from ..core.containers import get_config

    key = (caps.name, get_config().cache_key())
    cached = _BRANCHWATER_CACHE.get(key)
    if cached is not None:
        return cached
    try:
        rc = run_tool(
            caps,
            ["sourmash", "scripts", "pairwise", "--help"],
            logger=logger,
            check=False,
            stdout_path=os.devnull,
            log_prefix="sourmash",
        )
        result = rc == 0
    except (ToolExecutionError, MissingBinaryError, FileNotFoundError, OSError):
        result = False
    _BRANCHWATER_CACHE[key] = result
    return result


def _parse_pairwise_csv(
    path: Path, threshold: float, known: set[str], ksize: int
) -> dict[str, set[str]]:
    """Parse a branchwater ``pairwise`` edge list into a symmetric neighbour map.

    Keeps pairs whose containment-derived ANI estimate (containment^(1/k), the
    quantity ``sourmash compare --ani`` reports) is >= ``threshold``, mirroring
    the dense graph. Rows without an ``average_containment`` column fall back to
    the containment equivalent of Jaccard, 2j/(1+j). Self-edges are dropped, and
    only labels in ``known`` are kept, so a stray name cannot introduce a
    phantom node.
    """
    min_containment = threshold**ksize
    neighbors: dict[str, set[str]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            q = row.get("query_name", "")
            m = row.get("match_name", "")
            if q == m or q not in known or m not in known:
                continue
            try:
                raw = row.get("average_containment")
                if raw:
                    containment = float(raw)
                else:
                    jaccard = float(row.get("jaccard", "") or "nan")
                    containment = 2 * jaccard / (1 + jaccard)
            except (ValueError, ZeroDivisionError):
                continue
            if not containment >= min_containment:  # NaN also fails here
                continue
            neighbors.setdefault(q, set()).add(m)
            neighbors.setdefault(m, set()).add(q)
    return neighbors


def _sparse_greedy_cluster(
    labels: Sequence[str],
    neighbors: Mapping[str, Iterable[str]],
    name_by_label: Mapping[str, str],
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Greedy representative pick over a sparse adjacency map.

    Mirrors :func:`_greedy_cluster`: prefer the most-connected genome (ties broken
    by input order) as the representative, then claim its still-free neighbours.
    Genomes with no above-threshold neighbour become their own representatives.
    """
    order = sorted(
        range(len(labels)),
        key=lambda i: (-len(set(neighbors.get(labels[i], ()))), i),
    )
    assigned: dict[str, str] = {}  # label -> representative label
    reps: list[str] = []
    for i in order:
        lab = labels[i]
        if lab in assigned:
            continue
        reps.append(lab)
        assigned[lab] = lab
        for nb in neighbors.get(lab, ()):
            if nb not in assigned:
                assigned[nb] = lab

    clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    for rep in reps:
        clusters[name_by_label[rep]] = []
        status[name_by_label[rep]] = STATUS_REPRESENTATIVE
    for lab in labels:
        rep = assigned[lab]
        if lab == rep:
            continue
        clusters[name_by_label[rep]].append(name_by_label[lab])
        status[name_by_label[lab]] = STATUS_CONTAINED
    return clusters, status


def _read_compare_csv(path: Path) -> tuple[list[str], npt.NDArray[np.float64]]:
    with open(path, encoding="utf-8", newline="") as fo:
        labels = next(csv.reader(fo))
        # Parse the N x N body straight into a contiguous float array. At 1000s of
        # genomes this is far smaller and faster than a Python list-of-lists
        # (a 10k x 10k matrix is ~0.8 GB as float64 vs several GB of Python floats).
        matrix = np.loadtxt(fo, delimiter=",", ndmin=2)
    return labels, matrix


def _match_labels_to_genomes(labels: Sequence[str], genomes: Sequence[Path]) -> dict[str, str]:
    """Map a sourmash column label to a genome basename.

    sourmash labels are signature names: the file path for a signature
    sketched here, the record name for a workdir sketch. Match by genome
    basename, then record name, then stem.
    """
    out: dict[str, str] = {}
    by_name = {g.name: g.name for g in genomes}
    by_record = {record_name(g): g.name for g in genomes}
    stems = {g.stem: g.name for g in genomes}
    for label in labels:
        base = Path(label).name
        if base in by_name:
            out[label] = base
        elif label in by_record:  # a workdir sketch, named by record name
            out[label] = by_record[label]
        elif Path(label).stem in stems:
            out[label] = stems[Path(label).stem]
        else:
            out[label] = base
    return out


def _greedy_cluster(
    labels: list[str],
    sim: npt.NDArray[np.float64] | list[list[float]],
    name_by_label: dict[str, str],
    threshold: float,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    n = len(labels)
    # Boolean adjacency (>= threshold); the per-row vectorized ops below replace
    # the previous pure-Python O(n^2) double loops (much faster at 1000s genomes).
    adj = np.asarray(sim, dtype=float) >= threshold
    np.fill_diagonal(adj, False)
    connectivity = adj.sum(axis=1)
    # prefer well-connected genomes as representatives (stable: -connectivity, idx)
    order = np.lexsort((np.arange(n), -connectivity))

    assigned = np.full(n, -1, dtype=np.int64)  # member idx -> representative idx (-1 = free)
    reps: list[int] = []
    for i in order:
        if assigned[i] != -1:
            continue
        reps.append(int(i))
        assigned[i] = i
        # claim every still-free neighbour of i in one vectorized step
        claim = adj[i] & (assigned == -1)
        assigned[claim] = i

    clusters: dict[str, list[str]] = {}
    status: dict[str, str] = {}
    for rep_idx in reps:
        rep_name = name_by_label[labels[rep_idx]]
        clusters[rep_name] = []
        status[rep_name] = STATUS_REPRESENTATIVE
    for member_idx in range(n):
        rep_idx = int(assigned[member_idx])
        if member_idx == rep_idx:
            continue
        rep_name = name_by_label[labels[rep_idx]]
        member_name = name_by_label[labels[member_idx]]
        clusters[rep_name].append(member_name)
        status[member_name] = STATUS_CONTAINED
    return clusters, status
