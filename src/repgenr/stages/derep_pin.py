"""Pinned representatives: genomes kept as representatives whatever the clustering.

A user may need given genomes, such as the type or reference strains of
biovars, to be leaves of the tree even where the ANI threshold puts them in
one cluster with other genomes. ``dereplicate --pin`` names those genomes.
After the keeper rule and any refinement, every cluster holding a pinned
genome that is not its representative is split: its representative and each
pinned genome in it become representatives, and every other member joins
the one of them it is most similar to.

Similarity is the estimate ``dereplicate --tool sourmash`` clusters on, read
from the workdir sketches (:mod:`repgenr.reporting.member_ani`); missing
sketches are written first. A member
whose ANI to every candidate is unknown, or equal for several, stays with the
cluster's original representative, so the split depends only on the sketches
and the order of the candidates, not on the dereplicator.

A pinned genome failed by quality control is in no cluster; it is left out
with a warning, as pinning does not override quality control.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Collection, Mapping, Sequence
from pathlib import Path

from ..core.errors import UserInputError
from ..dereplicators.base import STATUS_FAIL_QC, STATUS_REPRESENTATIVE, DerepResult
from ..reporting.member_ani import load_hashes, pair_ani

# Names listed in a log line before "and N more".
_LOG_NAMES_MAX = 5

# (genome, genome) -> ANI estimate, or None when unknown.
AniOf = Callable[[str, str], "float | None"]


def read_pin_file(path: Path) -> list[str]:
    """The genomes named in a pin file: one per line, the first
    tab-separated column; blank lines and lines starting with '#' are
    skipped."""
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise UserInputError(f"--pin-file {path} cannot be read: {error}") from error
    names: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        names.append(line.split("\t")[0].strip())
    if not names:
        raise UserInputError(f"--pin-file {path} names no genome.")
    return names


def sketch_ani(sketches: Mapping[str, Path]) -> AniOf:
    """ANI of two genomes from their sketches (genome filename -> sketch
    file), each sketch read once; unknown when a sketch is missing."""
    cache: dict[str, frozenset[int] | None] = {}

    def hashes(name: str) -> frozenset[int] | None:
        if name not in cache:
            path = sketches.get(name)
            cache[name] = None if path is None else load_hashes(path)
        return cache[name]

    def ani(a: str, b: str) -> float | None:
        ha, hb = hashes(a), hashes(b)
        if ha is None or hb is None:
            return None
        return pair_ani(ha, hb)[0]

    return ani


def split_genomes(result: DerepResult, pinned: Collection[str]) -> list[str]:
    """The genomes of the clusters :func:`pin_representatives` splits, whose
    sketches it reads."""
    pinned_set = set(pinned)
    out: list[str] = []
    for rep_name, members in result.clusters.items():
        if any(m in pinned_set for m in members):
            out.extend([rep_name, *members])
    return out


def _names(names: Sequence[str]) -> str:
    shown = ", ".join(names[:_LOG_NAMES_MAX])
    if len(names) > _LOG_NAMES_MAX:
        shown += f" and {len(names) - _LOG_NAMES_MAX} more"
    return shown


def pin_representatives(
    result: DerepResult,
    pinned: Collection[str],
    ani: AniOf,
    logger: logging.Logger,
) -> tuple[DerepResult, list[str]]:
    """Make every genome in ``pinned`` a representative.

    Returns the new result and the pinned genomes that were not
    representatives before (promoted), in cluster order.
    """
    pinned_set = set(pinned)
    failed = sorted(n for n in pinned_set if result.genome_status.get(n) == STATUS_FAIL_QC)
    if failed:
        logger.warning(
            "Pinned genome(s) failed quality control and are left out: %s", _names(failed)
        )

    rep_paths = {p.name: p for p in result.representatives}
    new_reps: list[Path] = []
    new_clusters: dict[str, list[str]] = {}
    status = dict(result.genome_status)
    promoted: list[str] = []
    unknown_members = 0

    for rep_name, members in result.clusters.items():
        pins = [m for m in members if m in pinned_set]
        if not pins:
            new_reps.append(rep_paths[rep_name])
            new_clusters[rep_name] = list(members)
            continue
        candidates = [rep_name, *pins]
        assigned: dict[str, list[str]] = {c: [] for c in candidates}
        for member in members:
            if member in pinned_set:
                continue
            best, best_ani = rep_name, None
            for candidate in candidates:
                value = ani(member, candidate)
                if value is not None and (best_ani is None or value > best_ani):
                    best, best_ani = candidate, value
            if best_ani is None:
                unknown_members += 1
            assigned[best].append(member)
        for candidate in candidates:
            new_reps.append(rep_paths[rep_name].with_name(candidate))
            new_clusters[candidate] = assigned[candidate]
            status[candidate] = STATUS_REPRESENTATIVE
        promoted.extend(pins)
        logger.info(
            "Pinned: the cluster of %s (%d genomes) split into %d clusters for %s",
            rep_name,
            len(members) + 1,
            len(candidates),
            _names(pins),
        )

    if unknown_members:
        logger.warning(
            "%d member(s) of split clusters have no ANI estimate to any candidate "
            "and stay with the original representative",
            unknown_members,
        )
    return (
        DerepResult(
            representatives=new_reps,
            clusters=new_clusters,
            genome_status=status,
            genome_information=result.genome_information,
            settings=result.settings,
            qc_reasons=result.qc_reasons,
        ),
        promoted,
    )


def pinned_names(
    pins: Sequence[str],
    genomes: Sequence[Path],
    *,
    accessions: Mapping[str, str] | None = None,
    outgroup: Collection[str] = (),
) -> list[str]:
    """The genome filenames of ``pins`` (filename, stem or accession)."""
    from .derep_refine import resolve_seeds

    return resolve_seeds(pins, genomes, accessions=accessions, outgroup=outgroup, flag="--pin")
