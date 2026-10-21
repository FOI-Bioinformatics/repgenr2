"""Choose a reference genome as the medoid of the set by sketch ANI.

Every aligner and every mapping SNP typer projects the set onto one
reference, and the choice shapes mostly the within-cluster branching order:
swapping the reference changed 72 of 94 splits on a 50-genome set while the
deep structure held (docs/audit/scaling-audit.md). The alphabetically first
file is an arbitrary choice; the medoid, the genome with the largest summed
ANI to the others, is the one that is central to the set. It is expected to
reduce reference bias, not to remove it: in a set with an over-represented
genotype the medoid sits inside that genotype.
"""

from __future__ import annotations

import logging
import shutil
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..core.sketches import SketchProvider

REFERENCE_MEDOID = "medoid"
REFERENCE_FIRST = "first"
REFERENCE_RULES = (REFERENCE_MEDOID, REFERENCE_FIRST)


@dataclass(frozen=True)
class MedoidChoice:
    genome: Path
    # Mean ANI estimate of the medoid to the other genomes, in percent.
    mean_ani: float
    n_genomes: int


def ani_matrix(
    genomes: Sequence[Path],
    sketches: dict[Path, Path],
    work_dir: Path,
    threads: int,
    logger: logging.Logger,
) -> tuple[list[str], np.ndarray]:
    """Genome file names (sorted) and their N x N ANI-estimate matrix.

    Runs ``sourmash compare`` through the sourmash dereplicator over the
    given sketches (the workdir's k=31 sketches); nothing is sketched here.
    The dereplicator reads the signature sizes, so a pair with fewer than
    ``MIN_SHARED_HASHES`` shared hashes is NaN.
    """
    from ..dereplicators.sourmash import SourmashDereplicator

    tool = SourmashDereplicator()
    defaults = tool.capabilities.default_params
    labels, ani, name_by_label = tool._dense_ani_matrix(
        genomes,
        work_dir,
        int(defaults["ksize"]),
        int(defaults["scaled"]),
        logger,
        None,
        threads,
        too_large_hint="pass --reference with a genome name instead of the medoid rule.",
        sketches=sketches,
    )
    return [name_by_label[label] for label in labels], ani


def pick_medoid(names: Sequence[str], ani: np.ndarray) -> tuple[str, float]:
    """The name with the largest summed ANI to the others; ties go to the
    name that sorts first, so the choice is the same on every run. An
    unknown ANI (NaN: too few shared hashes to estimate one) counts
    as 0."""
    n = len(names)
    if n == 0:
        raise ValueError("no genomes to choose a medoid from")
    if n == 1:
        return names[0], 100.0
    # Fractions, 1.0 on the diagonal.
    matrix = np.nan_to_num(np.asarray(ani, dtype=float), nan=0.0)
    totals = matrix.sum(axis=1) - np.diag(matrix)
    best = min(range(n), key=lambda i: (-totals[i], names[i]))
    return names[best], float(totals[best] / (n - 1)) * 100.0


def choose_medoid(
    genomes: Sequence[Path],
    provider: SketchProvider | None,
    work_dir: Path,
    threads: int,
    logger: logging.Logger,
) -> MedoidChoice | None:
    """The medoid of ``genomes``, or None when it cannot be computed here
    (no sketch provider, sourmash missing, or the compare failed)."""
    from ..core.errors import MissingBinaryError, ToolExecutionError, WorkdirError
    from ..dereplicators.sourmash import _DENSE_MAX_GENOMES, SourmashDereplicator

    if provider is None or len(genomes) < 2:
        return None
    if len(genomes) > _DENSE_MAX_GENOMES:
        logger.warning(
            "Medoid reference not chosen: the dense ANI matrix is limited to %d genomes and "
            "there are %d; pass --reference with a genome name. The first genome is used.",
            _DENSE_MAX_GENOMES,
            len(genomes),
        )
        return None
    sketches = provider(genomes) or {}
    missing = [g for g in genomes if g not in sketches]
    if missing:
        # Only existing sketches are read (the provider writes none for this
        # consumer); the workdir sketches are made by `repgenr sketch` or the
        # sketch flag of the genome stages.
        logger.info(
            "Medoid reference not chosen: %d of %d genomes have no sketch under sketches/ "
            "(run with --sketch, or repgenr sketch, to enable it)",
            len(missing),
            len(genomes),
        )
        return None
    try:
        SourmashDereplicator().preflight()
    except MissingBinaryError as exc:
        logger.info("Medoid reference not chosen: %s", exc)
        return None
    try:
        names, ani = ani_matrix(genomes, sketches, work_dir, threads, logger)
    except (ToolExecutionError, WorkdirError) as exc:
        logger.warning("Medoid reference not chosen (%s); the first genome is used", exc)
        return None
    by_name = {g.name: g for g in genomes}
    name, mean_ani = pick_medoid(names, ani)
    genome = by_name.get(name)
    if genome is None:
        logger.warning("Medoid reference not chosen: %r is not one of the genomes", name)
        return None
    return MedoidChoice(genome=genome, mean_ani=mean_ani, n_genomes=len(genomes))


REFERENCE_NOTE = "reference.txt"


@dataclass(frozen=True)
class ReferenceChoice:
    """The reference genome and how it was chosen: ``rule`` is ``medoid``,
    ``first``, ``first-fallback`` (the medoid rule could not be applied) or
    ``name``; ``mean_ani`` is the medoid's mean ANI to the others, else None."""

    genome: Path
    rule: str
    mean_ani: float | None = None


def note_reference(directory: Path, choice: ReferenceChoice) -> None:
    """Keep the reference beside the alignment it was built against."""
    directory.mkdir(parents=True, exist_ok=True)
    from ..core.contracts import atomic_path

    with atomic_path(directory / REFERENCE_NOTE) as tmp:
        score = "" if choice.mean_ani is None else f"{choice.mean_ani:.2f}"
        tmp.write_text(f"{choice.genome.name}\t{choice.rule}\t{score}\n", encoding="utf-8")


def clear_reference_note(directory: Path) -> None:
    (directory / REFERENCE_NOTE).unlink(missing_ok=True)


def read_reference_choice(directory: Path) -> ReferenceChoice | None:
    """The reference recorded beside the alignment in ``directory``: the
    genome name (as a bare path), the rule and the medoid's mean ANI."""
    path = directory / REFERENCE_NOTE
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0]:
        return None
    fields = lines[0].split("\t")
    name = fields[0]
    rule = fields[1] if len(fields) > 1 and fields[1] else "name"
    score = fields[2] if len(fields) > 2 else ""
    try:
        mean_ani = float(score) if score else None
    except ValueError:
        mean_ani = None
    return ReferenceChoice(genome=Path(name), rule=rule, mean_ani=mean_ani)


def read_reference_note(directory: Path) -> str | None:
    """The reference genome recorded beside the alignment in ``directory``."""
    choice = read_reference_choice(directory)
    return choice.genome.name if choice is not None else None


def resolve_reference(
    reference_name: str | None,
    genomes: Sequence[Path],
    outgroup_file: Path | None,
    logger: logging.Logger,
    *,
    sketch_provider: SketchProvider | None = None,
    scratch_dir: Path | None = None,
    threads: int = 1,
) -> ReferenceChoice:
    """The reference genome: a name resolved by basename against the genome
    set, ``first`` (the alphabetically first genome), or the medoid by sketch
    ANI (``medoid``, the default) when the workdir sketches and sourmash are
    at hand. The default falls back to the first genome with a warning; an
    explicit ``medoid`` that cannot be computed is an error."""
    from ..core.errors import UserInputError

    requested = (reference_name or "").strip()
    rule = requested or REFERENCE_MEDOID
    if rule not in REFERENCE_RULES:
        pool = list(genomes)
        if outgroup_file is not None:
            pool.append(outgroup_file)
        for p in pool:
            if p.name == rule:
                logger.info("Using reference genome %s", p.name)
                return ReferenceChoice(p, "name")
        raise UserInputError(f"Reference genome not found: {rule}")
    if rule == REFERENCE_MEDOID:
        work = scratch_dir / "medoid" if scratch_dir is not None else Path(tempfile.mkdtemp())
        try:
            choice = choose_medoid(genomes, sketch_provider, work, threads, logger)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        if choice is not None:
            logger.info(
                "Reference genome %s: the medoid of %d genomes by sketch ANI (mean %.2f%% "
                "over its pairs to the others, a pair of unknown ANI counted as 0), a less "
                "arbitrary choice than the first filename and expected "
                "to reduce reference bias; pass --reference NAME to choose deliberately.",
                choice.genome.name,
                choice.n_genomes,
                choice.mean_ani,
            )
            return ReferenceChoice(choice.genome, REFERENCE_MEDOID, choice.mean_ani)
        if requested:
            raise UserInputError(
                "--reference medoid needs a sketch for every genome (--sketch, repgenr sketch) "
                "and sourmash; none could be chosen here. Use --reference first or a genome name."
            )
    logger.warning(
        "Projecting onto the alphabetically first genome '%s'%s. Every aligner and "
        "mapping typer builds the alignment relative to this choice, which can bias "
        "the tree: on sets with an over-represented genotype the first filename is "
        "usually one of its members. Pass --reference NAME to choose deliberately, or "
        "sketch the genomes (sourmash) so the medoid can be chosen.",
        genomes[0].name,
        "" if rule == REFERENCE_FIRST else " (no medoid: sketches or sourmash unavailable)",
    )
    return ReferenceChoice(
        genomes[0], REFERENCE_FIRST if rule == REFERENCE_FIRST else "first-fallback"
    )
