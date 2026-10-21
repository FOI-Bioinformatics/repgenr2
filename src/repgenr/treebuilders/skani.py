"""skani tree builder (alignment-free; ANI distances + neighbour joining).

``skani triangle`` estimates the average nucleotide identity of every pair
from sparse chaining of k-mer sketches, reliable down to about 80 percent ANI
(mash distances saturate below that; Jaccard distances compress near
identity), so this builder suits species through genus sets. The distance is
``1 - ANI`` from the ``--distance --full-matrix`` output, divided by 100. A
pair skani does not report (below its ``-s`` screen, or with too small an
aligned fraction) has ANI 0 in that matrix; its distance is set to 1 and the
count is logged, since such pairs are beyond the method, not identical.
The joining step is :func:`repgenr.tree.nj.nj_tree`. skani also writes an
aligned-fraction matrix beside the distances (``skani_distances.tsv.af``),
left in place for inspection. Command surface checked against skani 0.4.0
(bioconda); the container image below pins the older 0.3.2, whose surface
was not checked separately.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import atomic_replace, record_name
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities, preflight
from ..core.process import write_fofn
from ..tree.nj import NJ_AUTO, NJ_RAPIDNJ, RAPIDNJ, choose_backend, nj_tree
from .base import InputKind, TreeBuilder, TreeParams, as_genome_list

PRESETS = ("fast", "medium", "slow", "small-genomes")
SATURATED = 1.0  # the distance given to a pair skani did not report
DEFAULT_SCREEN = 80.0  # skani's own -s default, in percent identity


class SkaniBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="skani",
        reads_gzip=True,
        container="quay.io/biocontainers/skani:0.3.2--h79ce301_0",
        conda=("bioconda::skani",),
        required_binaries=(BinarySpec("skani", version_args=("-V",), min_version="0.2"),),
        # Declared, not measured here: skani's authors report all-pairs runs
        # of this size; the dense matrix and the joining step then dominate.
        recommended_max_genomes=10000,
        # screen: skani's -s, the approximate identity below which a pair is
        # not aligned (its default 80); min_af: --min-af; preset: one of
        # PRESETS; nj: the joining backend (repgenr.tree.nj).
        accepted_extras=frozenset({"screen", "min_af", "preset", "nj"}),
        default_params={"nj": NJ_AUTO},
        # A distance-based NJ tree has no bootstrap step.
        ignored_params=frozenset({"bootstrap"}),
    )
    input_kind = InputKind.GENOMES

    def _triangle(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
        *,
        preset_default: str | None = None,
        alternatives: str = "mashtree or sourmash",
    ) -> tuple[list[str], list[list[float]], int]:
        """Run ``skani triangle``; return (record names, distance matrix in
        [0, 1], number of pairs set to SATURATED), rows in name order.
        ``preset_default`` applies when no ``preset`` extra is given."""
        out_dir.mkdir(parents=True, exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        matrix = out_dir / "skani_distances.tsv"
        cmd: list[str | Path] = [
            "skani",
            "triangle",
            "-l",
            fofn,
            "-t",
            str(params.threads),
            "--full-matrix",
            "--distance",
            "-o",
            matrix,
        ]
        screen = _extra_float(params.extra, "screen", DEFAULT_SCREEN)
        if "screen" in params.extra:
            cmd += ["-s", f"{screen:g}"]
        if "min_af" in params.extra:
            cmd += ["--min-af", f"{_extra_float(params.extra, 'min_af', 15.0):g}"]
        preset = params.extra.get("preset", preset_default)
        if preset is not None:
            if str(preset) not in PRESETS:
                raise UserInputError(
                    f"skani preset must be one of {', '.join(PRESETS)}, not {preset!r}"
                )
            cmd.append(f"--{preset}")
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="skani",
            extra_mounts=genome_dirs,
        )
        if not matrix.is_file():
            raise WorkdirError("skani triangle wrote no matrix")
        names, dist, saturated = read_full_matrix(matrix, genomes)
        if saturated:
            logger.warning(
                "skani reported no ANI for %d of %d pairs (below its screen, about %g%% "
                "identity, or too little aligned); their distance is set to %g. A set this "
                "divergent is better served by %s.",
                saturated,
                len(names) * (len(names) - 1) // 2,
                screen,
                SATURATED,
                alternatives,
            )
        return names, dist, saturated

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        genomes = as_genome_list(msa_or_genomes)
        # Settle the joining step before the all-pairs run (repgenr.tree.nj).
        backend = choose_backend(
            str(params.extra.get("nj", self.capabilities.default_params["nj"])),
            len(genomes),
            logger,
        )
        if backend == NJ_RAPIDNJ:
            preflight(RAPIDNJ)
        names, dist, _ = self._triangle(genomes, out_dir, params, logger)
        newick = nj_tree(names, dist, out_dir, logger, backend=backend, threads=params.threads)
        tree = out_dir / "tree.nwk"
        with atomic_replace(tree) as fo:
            fo.write(newick + "\n")
        return tree

    def distance_matrix(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        """A square tab-separated matrix in mashtree's layout (an empty first
        header cell, then the record names), as the viral outgroup selection
        reads. Viral genomes are small, so skani's ``small-genomes`` preset
        applies unless a preset is given. When every pair between the selection
        (``S_*``) and the candidates (``O_*``), or every pair at all, is below
        skani's reach, no matrix is written: the candidates cannot be ranked by
        it, and the caller's handling of a missing matrix applies. When some
        candidates are within reach and others are not, the ones within reach
        win (the selection takes the lowest maximum distance), which is the
        intended reading: a saturated candidate cannot be ranked at all."""
        names, dist, _ = self._triangle(
            as_genome_list(genomes),
            out_dir,
            params,
            logger,
            preset_default="small-genomes",
            alternatives="mashtree",
        )
        matrix = out_dir / "distance_matrix.tsv"
        if _all_saturated(names, dist):
            logger.warning(
                "skani reported no ANI between the selection and any outgroup candidate "
                "(all below about %g%% identity); no distance matrix is written. Use "
                "--outgroup-treebuilder mashtree for candidates this distant.",
                _extra_float(params.extra, "screen", DEFAULT_SCREEN),
            )
            matrix.unlink(missing_ok=True)
            return matrix
        with atomic_replace(matrix) as fo:
            fo.write("\t" + "\t".join(names) + "\n")
            for name, row in zip(names, dist, strict=True):
                fo.write(name + "\t" + "\t".join(f"{d:.6f}" for d in row) + "\n")
        return matrix


def _extra_float(extra: Mapping[str, object], key: str, default: float) -> float:
    raw = extra.get(key, default)
    try:
        return float(raw)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise UserInputError(f"Extra parameter '{key}' must be a number, got {raw!r}.") from exc


def _all_saturated(names: Sequence[str], dist: Sequence[Sequence[float]]) -> bool:
    """True when no pair that the outgroup selection ranks on has a distance:
    every selection-to-candidate pair (``S_*`` against ``O_*`` names), or,
    without such names, every pair at all."""
    sel = [i for i, n in enumerate(names) if n.startswith("S_")]
    cand = [i for i, n in enumerate(names) if n.startswith("O_")]
    if sel and cand:
        pairs = [(i, j) for i in sel for j in cand]
    else:
        pairs = [(i, j) for i in range(len(names)) for j in range(i + 1, len(names))]
    return bool(pairs) and all(dist[i][j] >= SATURATED for i, j in pairs)


def read_full_matrix(
    path: Path, genomes: Sequence[Path]
) -> tuple[list[str], list[list[float]], int]:
    """Parse ``skani triangle --full-matrix --distance`` output: the count on
    the first line, then one row per input (its path, then the distances in
    percent). Rows are read one at a time and turned into floats as they come
    (the text of a 10000-genome matrix would otherwise sit in memory three
    times over), joined to the genomes by path and returned in record name
    order; a distance of 100 or more (no ANI reported) becomes SATURATED."""
    by_path = {os.path.abspath(g): record_name(g) for g in genomes}
    by_name = {record_name(g): record_name(g) for g in genomes}
    labels: list[str] = []
    raw_rows: list[list[float]] = []
    n: int | None = None
    with open(path, encoding="utf-8") as fo:
        for line in fo:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if n is None:
                try:
                    n = int(line.strip())
                except ValueError as exc:
                    raise WorkdirError(f"{path.name}: expected a count on the first line") from exc
                continue
            label, _, rest = line.partition("\t")
            name = by_path.get(os.path.abspath(label)) or by_name.get(record_name(label))
            if name is None:
                raise WorkdirError(f"{path.name}: row {label!r} matches none of the input genomes")
            try:
                values = [float(x) for x in rest.split("\t")]
            except ValueError as exc:
                raise WorkdirError(f"{path.name}: unreadable distance in row {label!r}") from exc
            if len(values) != n:
                raise WorkdirError(f"{path.name}: expected a {n} x {n} matrix with a label column")
            labels.append(name)
            raw_rows.append(values)
    if n is None:
        raise WorkdirError(f"{path.name} is empty")
    if len(raw_rows) != n:
        raise WorkdirError(f"{path.name}: expected a {n} x {n} matrix with a label column")
    if len(set(labels)) != len(labels):
        raise WorkdirError(f"{path.name}: two rows map to one genome")
    order = sorted(range(n), key=lambda i: labels[i])
    saturated = 0
    dist: list[list[float]] = []
    for i in order:
        row: list[float] = []
        for j in order:
            value = raw_rows[i][j]
            if i == j:
                row.append(0.0)
            elif value >= 100.0:
                row.append(SATURATED)
                saturated += 1
            else:
                row.append(max(value, 0.0) / 100.0)
        dist.append(row)
    return [labels[i] for i in order], dist, saturated // 2
