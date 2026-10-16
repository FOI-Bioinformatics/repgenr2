"""Quality-aware representative selection applied after any dereplicator.

Dereplicators pick representatives by connectivity, first-seen order or tool
defaults, which favours the most-sequenced genotype in a cluster. This step
re-picks each cluster's representative by assembly quality when the manifest
carries CheckM-style completeness and contamination, and leaves the adapter's
choice in place otherwise.

The score follows dRep's default weights for the terms that can be computed
without running tools: completeness - 5 x contamination + 0.5 x log10(N50),
with the N50 read from the genome FASTA. Equal scores are broken by higher
completeness, then lower contamination, then higher N50, then the genome
name, so the choice does not depend on the dereplicator. A genome without
quality never becomes the representative; a scored genome replaces an
unscored representative only when it is high quality (MIMAG: completeness
above 90, contamination below 5).

``--keeper gtdb`` adds one preference before that: a GTDB species
representative in the cluster is kept (:func:`prefer_gtdb_representatives`).

Reading the N50s is the only part of the step that touches the genome files.
Before the clusters are compared, the N50s they need are read in parallel
(:meth:`N50Lookup.prefetch`), one genome per worker thread; the comparison
itself then runs on cached values, in cluster order, so the choices and the
log lines are the same for any number of threads.
"""

from __future__ import annotations

import gzip
import logging
import math
import zlib
from collections.abc import Callable, Collection, Iterable, Mapping
from pathlib import Path

from ..core.errors import UserInputError
from ..core.executors import parallel_map
from ..core.process import is_gzip
from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE, DerepResult

CONTAMINATION_WEIGHT = 5.0
N50_WEIGHT = 0.5
# MIMAG high-quality draft thresholds (Bowers et al. 2017).
HQ_MIN_COMPLETENESS = 90.0
HQ_MAX_CONTAMINATION = 5.0
# Threads for the N50 reads when no -t is given: each holds about twice a
# genome in memory, and more threads gained little in a measurement.
N50_MAX_WORKERS = 8
# Names listed in a log line before "and N more".
_LOG_NAMES_MAX = 5

Quality = Mapping[str, tuple[float, float]]
# Genome filename -> N50 in bases, or None when the file cannot be found.
N50Of = Callable[[str], "int | None"]


def quality_score(completeness: float, contamination: float) -> float:
    """CheckM part of the score: completeness minus five times contamination."""
    return completeness - CONTAMINATION_WEIGHT * contamination


def keeper_score(completeness: float, contamination: float, n50: int | None) -> float:
    """dRep-style score: ``quality_score`` plus 0.5 x log10(N50).

    An unknown N50 adds nothing.
    """
    score = quality_score(completeness, contamination)
    if n50:
        score += N50_WEIGHT * math.log10(n50)
    return score


def is_high_quality(completeness: float, contamination: float) -> bool:
    """MIMAG high quality: completeness above 90 and contamination below 5."""
    return completeness > HQ_MIN_COMPLETENESS and contamination < HQ_MAX_CONTAMINATION


def genome_n50(path: Path) -> int:
    """N50 of the sequences in a (possibly gzipped) FASTA file.

    The file is read once as bytes; sequence lengths are the bytes between
    header lines less the line breaks, so nothing is decoded.
    """
    try:
        if is_gzip(path):
            with gzip.open(path, "rb") as gz:
                data = gz.read()
        else:
            data = path.read_bytes()
    except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise UserInputError(
            f"Genome {path} is a truncated or corrupt gzip file ({exc}); "
            "replace it with a complete copy."
        ) from exc
    lengths: list[int] = []
    # Records start at ">" at the start of the file or of a line.
    records = data.split(b"\n>")
    if records and records[0].startswith(b">"):
        records[0] = records[0][1:]
    elif records:
        records = records[1:]  # text before the first header is not a sequence
    for record in records:
        newline = record.find(b"\n")
        if newline < 0:
            lengths.append(0)
            continue
        # Counted in place: a slice of the body would copy it.
        start = newline + 1
        breaks = record.count(b"\n", start) + record.count(b"\r", start)
        lengths.append(len(record) - start - breaks - record.count(b" ", start))
    lengths.sort(reverse=True)
    half = sum(lengths) / 2
    running = 0
    for length in lengths:
        running += length
        if running >= half:
            return length
    return 0


class N50Lookup:
    """Genome filename -> N50, read from the first directory holding the file.

    Each genome is read at most once per lookup, and only when asked for, so
    a run reads only the genomes whose scores are compared. ``known`` seeds
    values computed elsewhere (a chunk's ``genome_n50.tsv``). A genome found
    in no directory has no N50 (``None``), which adds nothing to its score.
    :meth:`prefetch` reads a batch of genomes with ``workers`` threads.
    """

    def __init__(
        self,
        dirs: Iterable[Path],
        *,
        known: Mapping[str, int] | None = None,
        workers: int = 1,
    ) -> None:
        self._dirs = list(dirs)
        self._cache: dict[str, int | None] = dict(known or {})
        self._workers = max(1, workers)

    def prefetch(self, names: Iterable[str]) -> None:
        """Read the N50 of each genome not yet known, ``workers`` at a time.

        The reads are mostly file I/O and zlib decompression, which release
        the GIL, so threads run them largely in parallel. The cache is filled from this thread
        only, after the reads.
        """
        todo = sorted({n for n in names if n not in self._cache})
        values = parallel_map(self._read, todo, self._workers)
        self._cache.update(zip(todo, values, strict=True))

    def __call__(self, name: str) -> int | None:
        if name not in self._cache:
            self._cache[name] = self._read(name)
        return self._cache[name]

    def _read(self, name: str) -> int | None:
        for d in self._dirs:
            path = d / name
            if path.is_file():
                return genome_n50(path)
        return None

    def missing(self) -> list[str]:
        """Genomes asked for whose file was found in no directory."""
        return sorted(name for name, n50 in self._cache.items() if n50 is None)

    def computed(self) -> dict[str, int]:
        """The N50 values known so far (genomes that were found)."""
        return {name: n50 for name, n50 in self._cache.items() if n50 is not None}


def warn_missing_n50(lookup: N50Lookup, logger: logging.Logger) -> None:
    """WARNING naming the scored genomes whose file was not found: their keeper
    score and ``rep_n50`` carry no N50 term."""
    missing = lookup.missing()
    if missing:
        logger.warning(
            "No genome file for %d scored genome(s), so their keeper score has no "
            "N50 term and rep_n50 is blank: %s",
            len(missing),
            log_names(missing),
        )


def _no_n50(_name: str) -> int | None:
    return None


def rank_key(name: str, quality: Quality, n50: N50Of) -> tuple[float, float, float, int, str]:
    """Sort key, best first, for a scored genome: score, completeness,
    contamination, N50, name."""
    completeness, contamination = quality[name]
    value = n50(name)
    return (
        -keeper_score(completeness, contamination, value),
        -completeness,
        contamination,
        -(value or 0),
        name,
    )


def _compared(incumbent: str, others: Iterable[str], quality: Quality) -> list[str]:
    """The genomes :func:`choose_keeper` ranks: the scored ones, limited to
    high-quality genomes when the incumbent has no quality."""
    pool = list(dict.fromkeys([incumbent, *others]))
    scored = [n for n in pool if n in quality]
    if incumbent not in quality:
        scored = [n for n in scored if is_high_quality(*quality[n])]
    return scored


def prefetch_n50(
    n50: N50Of | None,
    groups: Iterable[tuple[str, Iterable[str]]],
    quality: Quality,
    *,
    scored_incumbents: bool = False,
) -> None:
    """Read in parallel the N50s the keeper rule needs for ``groups``, each an
    (incumbent, others) pair, when ``n50`` is an :class:`N50Lookup`.

    These are the genomes :func:`choose_keeper` ranks when it compares two or
    more, and the one it promotes alone over an unscored incumbent (its score
    is logged). ``scored_incumbents`` adds each scored incumbent, whose N50
    the cluster summary reports. Any other ``n50`` is called lazily.
    """
    if not isinstance(n50, N50Lookup):
        return
    names: set[str] = set()
    for incumbent, others in groups:
        compared = _compared(incumbent, others, quality)
        if len(compared) > 1 or compared != [incumbent]:
            names.update(compared)
        if scored_incumbents and incumbent in quality:
            names.add(incumbent)
    n50.prefetch(names)


def choose_keeper(
    incumbent: str,
    others: Iterable[str],
    quality: Quality,
    n50: N50Of | None = None,
) -> str:
    """The genome that should represent ``incumbent`` and ``others``.

    Only scored genomes compete. When the incumbent (the tool's pick) is
    scored, the best scored genome wins, the incumbent included, ranked by
    :func:`rank_key`. When it is not, a scored genome replaces it only if it
    is high quality (:func:`is_high_quality`); otherwise the incumbent stays.
    """
    n50_of = n50 or _no_n50
    scored = _compared(incumbent, others, quality)
    if not scored:
        return incumbent
    if len(scored) == 1:
        # Nothing to compare: no N50 is read.
        return scored[0]
    return min(scored, key=lambda n: rank_key(n, quality, n50_of))


def best_score(name: str, quality: Quality, n50: N50Of | None = None) -> float | None:
    """:func:`keeper_score` of a genome, or None when it has no quality."""
    q = quality.get(name)
    if q is None:
        return None
    return keeper_score(q[0], q[1], (n50 or _no_n50)(name))


def log_names(names: list[str]) -> str:
    """At most five names, then "and N more"."""
    shown = ", ".join(names[:_LOG_NAMES_MAX])
    extra = len(names) - _LOG_NAMES_MAX
    return shown if extra <= 0 else f"{shown} and {extra} more"


def log_unpromoted(left: list[str], logger: logging.Logger, what: str = "cluster(s)") -> None:
    """Report groups whose unscored representative was kept although scored
    genomes were available, none of them high quality."""
    if left:
        logger.info(
            "Keeper: %d %s keep the tool's representative, which has no quality, "
            "because no scored genome in them is high quality (completeness > %g, "
            "contamination < %g): %s",
            len(left),
            what,
            HQ_MIN_COMPLETENESS,
            HQ_MAX_CONTAMINATION,
            log_names(sorted(left)),
        )


def rescore_representatives(
    result: DerepResult,
    quality: Quality,
    logger: logging.Logger,
    n50: N50Of | None = None,
) -> tuple[DerepResult, int]:
    """Return a copy of ``result`` whose representatives are chosen by
    :func:`choose_keeper`, and the number of clusters whose representative
    changed. ``n50`` gives each genome's N50; without it the score has no
    N50 term.
    """
    if not quality:
        return result, 0
    prefetch_n50(n50, result.clusters.items(), quality)

    rep_paths = {p.name: p for p in result.representatives}
    new_reps: list[Path] = []
    new_clusters: dict[str, list[str]] = {}
    status = dict(result.genome_status)
    swaps = 0
    left_to_tool: list[str] = []

    for rep_name, members in result.clusters.items():
        best_name = choose_keeper(rep_name, members, quality, n50)
        if best_name == rep_name:
            if rep_name not in quality and any(m in quality for m in members):
                left_to_tool.append(rep_name)
            new_reps.append(rep_paths[rep_name])
            new_clusters[rep_name] = list(members)
            continue
        swaps += 1
        new_reps.append(rep_paths[rep_name].with_name(best_name))
        new_clusters[best_name] = sorted(n for n in [rep_name, *members] if n != best_name)
        status[best_name] = STATUS_REPRESENTATIVE
        status[rep_name] = STATUS_CONTAINED
        current = best_score(rep_name, quality, n50)
        logger.info(
            "Keeper: %s replaces %s (score %.2f vs %s)",
            best_name,
            rep_name,
            best_score(best_name, quality, n50),
            "n/a" if current is None else f"{current:.2f}",
        )

    log_unpromoted(left_to_tool, logger)
    if swaps:
        logger.info(
            "Quality-aware keeper changed %d of %d representatives", swaps, len(new_clusters)
        )
    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
        genome_information=result.genome_information,
    ), swaps


def prefer_gtdb_representatives(
    result: DerepResult,
    gtdb_representatives: Collection[str],
    quality: Quality,
    logger: logging.Logger,
    n50: N50Of | None = None,
) -> tuple[DerepResult, int]:
    """Keep a GTDB species representative as each cluster's representative.

    ``gtdb_representatives`` holds the filenames GTDB marks as species
    representatives. The quality rule (:func:`rescore_representatives`) is
    applied first, so a cluster without a GTDB representative is treated as
    under ``--keeper quality``. A cluster holding GTDB representatives keeps
    one of them, chosen among them alone by :func:`choose_keeper`: the
    incumbent is the current representative when it is one of them, else the
    best-ranked scored one, else the first by name. Several in one cluster
    means the ANI threshold joined GTDB species; one keeper per cluster is
    the contract, so the others stay contained and a warning names them.

    Returns the new result and the number of clusters whose representative
    differs from the tool's pick.
    """
    n50_of = n50 or _no_n50
    tool_picks = set(result.clusters)
    result, _ = rescore_representatives(result, quality, logger, n50)
    flagged = set(gtdb_representatives)
    if isinstance(n50, N50Lookup):
        # The genomes ranked below, read in parallel first: among two or more
        # GTDB representatives, those choose_keeper compares when the
        # representative is one of them, else the scored ones.
        ranked: set[str] = set()
        for rep_name, members in result.clusters.items():
            in_cluster = [n for n in (rep_name, *members) if n in flagged]
            if len(in_cluster) < 2:
                continue
            if rep_name in flagged:
                ranked.update(_compared(rep_name, in_cluster, quality))
            else:
                scored = [n for n in in_cluster if n in quality]
                if len(scored) > 1:
                    ranked.update(scored)
        n50.prefetch(ranked)

    rep_paths = {p.name: p for p in result.representatives}
    new_reps: list[Path] = []
    new_clusters: dict[str, list[str]] = {}
    status = dict(result.genome_status)
    promoted = 0

    for rep_name, members in result.clusters.items():
        candidates = [rep_name, *members]
        in_cluster = [n for n in candidates if n in flagged]
        keeper = rep_name
        if len(in_cluster) == 1:
            # One GTDB representative: kept without ranking, so no N50 is read.
            keeper = in_cluster[0]
        elif in_cluster:
            scored = [n for n in in_cluster if n in quality]
            if rep_name in flagged:
                incumbent = rep_name
            elif len(scored) > 1:
                incumbent = min(scored, key=lambda n: rank_key(n, quality, n50_of))
            elif scored:
                incumbent = scored[0]
            else:
                incumbent = min(in_cluster)
            keeper = choose_keeper(
                incumbent, [n for n in in_cluster if n != incumbent], quality, n50
            )
        if len(in_cluster) > 1:
            logger.warning(
                "Keeper: the cluster of %s holds %d GTDB species representatives; "
                "keeping %s, the others stay contained: %s. The ANI threshold "
                "joins these GTDB species.",
                keeper,
                len(in_cluster),
                keeper,
                log_names(sorted(n for n in in_cluster if n != keeper)),
            )
        if keeper == rep_name:
            new_reps.append(rep_paths[rep_name])
            new_clusters[rep_name] = list(members)
            continue
        promoted += 1
        new_reps.append(rep_paths[rep_name].with_name(keeper))
        new_clusters[keeper] = sorted(n for n in candidates if n != keeper)
        status[keeper] = STATUS_REPRESENTATIVE
        status[rep_name] = STATUS_CONTAINED
        logger.info("Keeper: GTDB species representative %s replaces %s", keeper, rep_name)

    swaps = sum(1 for rep in new_clusters if rep not in tool_picks)
    logger.info(
        "GTDB keeper: %d of %d representatives are GTDB species representatives; "
        "%d differ from the tool's pick",
        sum(1 for rep in new_clusters if rep in flagged),
        len(new_clusters),
        swaps,
    )
    if promoted:
        logger.info("GTDB keeper promoted %d GTDB species representative(s)", promoted)
    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
        genome_information=result.genome_information,
    ), swaps
