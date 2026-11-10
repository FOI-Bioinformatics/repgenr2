"""Choice of the outgroup genome of a GTDB selection.

Three modes (``metadata --outgroup-mode``):

``sister`` (default)
    The candidates are the species representatives of the group sister to
    the target in the GTDB reference tree (:mod:`repgenr.core.gtdb_tree`):
    the closest genus to a genus, the closest family to a family, the
    closest species to a species. When no genome of the sister group meets
    the minimum quality, the next group out is taken, and so on. When the
    group holds metagenome-assembled or single-cell genomes only, an isolate
    of the next group out is preferred (:func:`choose_by_groups`).
``quality``
    The candidates are the species representatives of the parent taxon (one
    rank up) outside the target, without regard to their distance.
``first``
    The first species representative of the parent taxon in table order (the
    rule before 2.3.0), kept to reproduce earlier selections.

Within the candidates, :func:`quality_key` orders the genomes: assemblies
from an isolate before metagenome-assembled or single-cell genomes, then
higher CheckM2 completeness - 5 x contamination, then type material, then
larger contig N50 and fewer contigs. Scores within :data:`SCORE_TOLERANCE`
of the best count as equal, so type material and N50 decide between
genomes of practically the same completeness. A candidate below
:data:`MIN_COMPLETENESS` or above :data:`MAX_CONTAMINATION` is not chosen
while one that meets both exists.

Every candidate considered is written to ``outgroup_choice.tsv`` with its
group, its path length to the target in the tree and its quality, so the
choice can be checked.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from .contracts import atomic_replace

OUTGROUP_CHOICE_TSV = "outgroup_choice.tsv"
OUTGROUP_MODES = ("sister", "quality", "first")
DEFAULT_OUTGROUP_MODE = "sister"
# The MIMAG medium-quality bounds: a candidate outside them is taken only
# when no candidate of any group meets them.
MIN_COMPLETENESS = 50.0
MAX_CONTAMINATION = 10.0
_NOT_ISOLATE = ("metagenome", "environmental", "single cell")


@dataclass(frozen=True)
class Candidate:
    accession: str
    species: str
    completeness: float | None = None
    contamination: float | None = None
    n50: int | None = None
    contigs: int | None = None
    # NCBI genome category: empty for an isolate, otherwise e.g.
    # "derived from metagenome".
    category: str = ""
    # GTDB type designation, e.g. "type strain of species"; empty when not.
    type_material: str = ""
    group: int | None = None  # 1 = sister group, 2 = next one out, ...
    distance: float | None = None  # path length to the target's MRCA

    @property
    def isolate(self) -> bool:
        cat = self.category.lower()
        return not any(word in cat for word in _NOT_ISOLATE)

    @property
    def is_type(self) -> bool:
        t = self.type_material.lower()
        return bool(t) and not t.startswith("not ")

    @property
    def score(self) -> float | None:
        if self.completeness is None:
            return None
        return self.completeness - 5 * (self.contamination or 0.0)

    @property
    def acceptable(self) -> bool:
        return (
            self.completeness is not None
            and self.completeness >= MIN_COMPLETENESS
            and (self.contamination or 0.0) <= MAX_CONTAMINATION
        )


def quality_key(c: Candidate) -> tuple:
    """Sort key, best first."""
    score = c.score
    return (
        not c.acceptable,
        not c.isolate,
        -(score if score is not None else -math.inf),
        not c.is_type,
        -(c.n50 or 0),
        c.contigs if c.contigs is not None else math.inf,
        c.accession,
    )


# Scores (completeness - 5 x contamination) within this many points of the
# best are treated as equal; type material and contig N50 then decide.
SCORE_TOLERANCE = 1.0


def _pick(ranked: Sequence[Candidate]) -> Candidate:
    """The best of ``ranked`` (sorted by :func:`quality_key`), with scores
    within :data:`SCORE_TOLERANCE` of the best treated as equal."""
    best = ranked[0]
    top = best.score
    if top is None:
        return best
    tied = [
        c
        for c in ranked
        if c.acceptable == best.acceptable
        and c.isolate == best.isolate
        and c.score is not None
        and c.score >= top - SCORE_TOLERANCE
    ]
    return min(tied, key=lambda c: (not c.is_type, -(c.n50 or 0), quality_key(c)))


@dataclass(frozen=True)
class Choice:
    chosen: Candidate
    considered: list[Candidate]
    mode: str
    reason: str


def _where(number: int) -> str:
    return "the sister group" if number == 1 else f"group {number} out from the target"


def choose_by_groups(
    groups: Iterable[Sequence[Candidate]], mode: str = "sister", isolate_lookahead: int = 1
) -> Choice | None:
    """The best candidate of the nearest group that has an acceptable one.

    When that group holds metagenome-assembled or single-cell genomes only,
    up to ``isolate_lookahead`` further groups are searched for an acceptable
    isolate, which is then taken: any group outside the target roots the
    tree, and an isolate assembly is usually more complete and less
    fragmented. When no group has an acceptable genome, the best of the
    nearest non-empty group. ``groups`` is consumed only as far as needed.
    """
    considered: list[Candidate] = []
    first: Candidate | None = None
    fallback: tuple[int, Candidate] | None = None
    for number, group in enumerate(groups, 1):
        if not group:
            continue
        ranked = sorted(group, key=quality_key)
        considered.extend(ranked)
        best = _pick(ranked)
        if first is None:
            first = best
        if best.acceptable and best.isolate:
            reason = f"best quality in {_where(number)}"
            if fallback is not None:
                reason = (
                    f"best isolate in {_where(number)}; {_where(fallback[0])} holds "
                    "metagenome-assembled or single-cell genomes only"
                )
            return Choice(best, considered, mode, reason)
        if best.acceptable and fallback is None:
            fallback = (number, best)
        if fallback is not None and number >= fallback[0] + isolate_lookahead:
            break
    if fallback is not None:
        return Choice(
            fallback[1],
            considered,
            mode,
            f"best quality in {_where(fallback[0])} (no isolate within {isolate_lookahead} "
            "further group)",
        )
    if first is None:
        return None
    return Choice(
        first,
        considered,
        mode,
        "no candidate meets the minimum quality; best of the sister group",
    )


def choose_by_quality(candidates: Sequence[Candidate]) -> Choice | None:
    """The best candidate of the parent taxon (mode ``quality``)."""
    if not candidates:
        return None
    ranked = sorted(candidates, key=quality_key)
    return Choice(_pick(ranked), ranked, "quality", "best quality in the parent taxon")


_COLUMNS = (
    "accession",
    "species",
    "chosen",
    "group",
    "tree_distance",
    "completeness",
    "contamination",
    "n50",
    "contigs",
    "isolate",
    "type_material",
    "mode",
    "reason",
)


def write_choice(workdir: Path, choice: Choice | None) -> None:
    """``outgroup_choice.tsv``: every candidate considered, the chosen one
    first. Removed when there is no automatic choice (a named outgroup or
    --no-outgroup)."""
    path = workdir / OUTGROUP_CHOICE_TSV
    if choice is None:
        path.unlink(missing_ok=True)
        return
    rows = [choice.chosen] + [c for c in choice.considered if c is not choice.chosen]
    with atomic_replace(path, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(_COLUMNS)
        for c in rows:
            chosen = c is choice.chosen
            writer.writerow(
                [
                    c.accession,
                    c.species,
                    "yes" if chosen else "no",
                    "" if c.group is None else c.group,
                    "" if c.distance is None else f"{c.distance:.5f}",
                    "" if c.completeness is None else c.completeness,
                    "" if c.contamination is None else c.contamination,
                    "" if c.n50 is None else c.n50,
                    "" if c.contigs is None else c.contigs,
                    "yes" if c.isolate else "no",
                    c.type_material,
                    choice.mode,
                    choice.reason if chosen else "",
                ]
            )


def read_choice(workdir: Path) -> list[dict[str, str]]:
    path = workdir / OUTGROUP_CHOICE_TSV
    if not path.is_file():
        return []
    with open(path, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


# An outgroup this similar to a genome of the selection probably lies within
# the target species rather than outside the target.
CLOSE_OUTGROUP_ANI = 0.95


def check_outgroup_distance(
    genomes: Sequence[Path], outgroups: Sequence[Path], sketches_dir: Path, logger
) -> float | None:
    """The highest sketch ANI estimate (k=31) of the outgroup to a genome of
    the selection, with a warning when it reaches :data:`CLOSE_OUTGROUP_ANI`.

    None when a sketch is missing or no pair shares enough hashes for an
    estimate, which is expected between genera (below about 80% ANI).
    """
    from ..reporting.member_ani import SKETCH_SUFFIX, load_hashes, pair_ani
    from .contracts import record_name

    best: tuple[float, str, str] | None = None
    for og in outgroups:
        og_hashes = load_hashes(sketches_dir / (record_name(og.name) + SKETCH_SUFFIX))
        if not og_hashes:
            continue
        for genome in genomes:
            hashes = load_hashes(sketches_dir / (record_name(genome.name) + SKETCH_SUFFIX))
            if not hashes:
                continue
            ani, _shared = pair_ani(og_hashes, hashes)
            if ani is not None and (best is None or ani > best[0]):
                best = (ani, og.name, genome.name)
    if best is None:
        logger.info(
            "Outgroup check: no sketch ANI estimate between the outgroup and the selection "
            "(expected when they are of different genera)"
        )
        return None
    ani, og_name, genome_name = best
    if ani >= CLOSE_OUTGROUP_ANI:
        logger.warning(
            "The outgroup %s has an estimated ANI of %.3f to %s (sketches, k=31): it may "
            "lie within the target rather than outside it. Check outgroup_choice.tsv, or "
            "name an outgroup with --outgroup-accession.",
            og_name,
            ani,
            genome_name,
        )
    else:
        logger.info(
            "Outgroup check: highest sketch ANI estimate to the selection %.3f (%s)",
            ani,
            genome_name,
        )
    return ani
