"""Dereplicator interface.

An adapter takes a set of genome FASTAs and returns a :class:`DerepResult`
(representatives + cluster membership + per-genome status). The adapter never
writes contract files itself -- the dereplicate stage normalizes the result into
``derep/representatives/`` + ``clusters.tsv`` + ``genome_status.tsv``. This keeps
every dereplicator interchangeable with zero downstream change.
"""

from __future__ import annotations

import csv
import logging
from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..core.errors import WorkdirError
from ..core.plugins import Registry, ToolCapabilities, preflight

registry: Registry[Dereplicator] = Registry("repgenr.dereplicators")

# status values used in DerepResult.genome_status
STATUS_REPRESENTATIVE = "representative"
STATUS_CONTAINED = "contained"
STATUS_FAIL_QC = "fail_qc"


@dataclass
class DerepParams:
    """Normalized dereplication parameters shared across tools.

    ``extra`` carries tool-specific overrides keyed by adapter name.
    ``quality`` maps genome filename to (completeness, contamination) from the
    manifest or selection.tsv, for genomes that carry both values; adapters
    that can use genome quality read it (dRep's ``--genomeInfo``, galah's
    input order) and the others ignore it.
    ``sketches`` maps genome path to its sourmash signature file
    (``sketches/<record name>.sig.zip``: k=21, 31 and 51 at scaled=1000, each
    signature named by the genome's record name). The stage fills it only for
    an adapter whose :meth:`Dereplicator.sketch_request` asks for parameters
    those files hold; the adapter sketches any genome left out.
    """

    primary_ani: float = 0.90
    secondary_ani: float = 0.99
    aligned_fraction: float = 0.50
    threads: int = 16
    extra: dict = field(default_factory=dict)
    quality: dict[str, tuple[float, float]] = field(default_factory=dict)
    sketches: Mapping[Path, Path] | None = None


@dataclass
class CompareResult:
    """Output of an all-vs-all comparison (the glance stage's contract).

    ``similarity_csv`` holds pairwise rows with at least the columns
    ``genome1``, ``genome2``, ``similarity``; ``dendrogram`` is an optional
    pre-rendered clustering figure. ``measure`` names the similarity for the
    plot axes (for example "MASH ANI" for dRep, "ANI" for sourmash's
    sketch-based estimate). The helpers in :mod:`.compare_io` write both files
    from a dense similarity matrix.
    """

    similarity_csv: Path | None = None
    dendrogram: Path | None = None
    measure: str = "ANI"


@dataclass
class DerepResult:
    """Normalized output every dereplicator must return."""

    representatives: list[Path]
    clusters: dict[str, list[str]]  # representative filename -> contained filenames
    genome_status: dict[str, str]  # genome filename -> status
    genome_information: list[dict] | None = None  # optional checkM-like QC rows
    # How the tool computed the result, for the stage record (``tool_settings``):
    # the sourmash back-end and ANI estimate, for example. Empty when the
    # adapter has nothing to add to its parameters.
    settings: dict = field(default_factory=dict)


_VALID_STATUS = frozenset({STATUS_REPRESENTATIVE, STATUS_CONTAINED, STATUS_FAIL_QC})


def combine_settings(results: Iterable[DerepResult]) -> dict:
    """The ``settings`` of several adapter calls (the chunks and merge passes
    of one dereplication) as one dict.

    A key whose value is the same in every call that reports settings is
    kept. Differing back-ends (a chunk that fell back to the dense path) are
    joined, sorted, with a comma; any other differing key, such as the
    per-call count of unknown pairs, is left out, since no single value
    describes the run. Calls without settings (results read back from disk)
    do not take part.
    """
    reported = [r.settings for r in results if r.settings]
    if not reported:
        return {}
    combined: dict = {}
    for key in dict.fromkeys(k for s in reported for k in s):
        values = [s[key] for s in reported if key in s]
        distinct = list(dict.fromkeys(values))
        if len(values) == len(reported) and len(distinct) == 1:
            combined[key] = distinct[0]
        elif key == "backend":
            combined[key] = ",".join(sorted(str(v) for v in distinct))
    return combined


def check_result_complete(result: DerepResult, genome_names: Collection[str]) -> None:
    """Refuse a result that silently drops genomes.

    Every input genome must carry a known status; every representative must be
    marked as such; every contained genome must belong to exactly one cluster.
    Adapters that parse tool tables can otherwise return an empty membership
    (e.g. after an upstream column-layout change) and the stage would complete
    with genomes missing from every deliverable.
    """
    names = set(genome_names)
    status = result.genome_status
    missing = sorted(names - status.keys())
    if missing:
        raise WorkdirError(
            f"Dereplication left {len(missing)} of {len(names)} genome(s) without a "
            f"status (e.g. {', '.join(missing[:3])}). The dereplication result is incomplete."
        )
    bad = sorted(f"{g}={s}" for g, s in status.items() if s not in _VALID_STATUS)
    if bad:
        raise WorkdirError(f"Unknown dereplication status value(s): {', '.join(bad[:3])}")

    rep_names = {p.name for p in result.representatives}
    unmarked = sorted(r for r in rep_names if status.get(r) != STATUS_REPRESENTATIVE)
    if unmarked:
        raise WorkdirError(
            f"{len(unmarked)} representative(s) are not marked as such in genome_status "
            f"(e.g. {', '.join(unmarked[:3])})."
        )

    home: dict[str, int] = {}
    for rep, members in result.clusters.items():
        for m in members:
            if m != rep:
                home[m] = home.get(m, 0) + 1
    orphans = sorted(g for g, s in status.items() if s == STATUS_CONTAINED and home.get(g, 0) == 0)
    if orphans:
        raise WorkdirError(
            f"{len(orphans)} contained genome(s) have no representative "
            f"(e.g. {', '.join(orphans[:3])}). The adapter returned an empty cluster table."
        )
    doubled = sorted(g for g, n in home.items() if n > 1)
    if doubled:
        raise WorkdirError(
            f"{len(doubled)} genome(s) appear in more than one cluster "
            f"(e.g. {', '.join(doubled[:3])})."
        )


def run_quality(
    quality: Mapping[str, tuple[float, float]],
    names: Collection[str],
    logger: logging.Logger,
    *,
    source: str,
) -> dict[str, tuple[float, float]]:
    """The genome quality given to the dereplicator for a whole run.

    Adapters that read ``DerepParams.quality`` (dRep's ``--genomeInfo``,
    galah's ``--genome-info``) need values for every genome they receive. The
    decision is taken once per run, over all genomes in ``names``, so that
    every chunk and the merge pass use the same source: the values when they
    cover every genome, otherwise none (dRep then runs CheckM, galah orders
    the genomes by file size). ``--keeper quality`` uses the partial values
    either way.
    """
    lacking = sorted(n for n in names if n not in quality)
    if not lacking:
        return {n: quality[n] for n in names}
    if len(lacking) < len(names):
        logger.warning(
            "The %s has completeness and contamination for %d of %d genomes (missing e.g. %s); "
            "the dereplicator is given no genome quality for this run, so dRep scores "
            "genomes with CheckM and galah takes them by descending file size.",
            source,
            len(names) - len(lacking),
            len(names),
            ", ".join(lacking[:3]),
        )
    return {}


def write_genome_info(path: Path, rows: Iterable[tuple[str, float, float]]) -> Path:
    """Write a dRep-style genome info table (genome,completeness,contamination).

    dRep and galah both read this layout; each matches the genome column in its
    own way (dRep the basename of the file it reads, galah the name without its
    FASTA suffix, as ``strip_fasta_suffix`` gives).
    """
    with open(path, "w", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, lineterminator="\n")
        writer.writerow(["genome", "completeness", "contamination"])
        for genome, completeness, contamination in rows:
            writer.writerow([genome, f"{completeness:g}", f"{contamination:g}"])
    return path


class Dereplicator(ABC):
    """Base class for dereplication adapters."""

    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        """Confirm required binaries are present; return resolved versions."""
        return preflight(self.capabilities)

    def sketch_request(self, extra: Mapping[str, object]) -> tuple[int, int] | None:
        """The (ksize, scaled) of the sourmash sketches this adapter compares at.

        None (the default) means the adapter reads no genome sketches. An
        adapter that returns a pair receives ``DerepParams.sketches`` (and the
        ``sketches`` argument of :meth:`compare`) when the workdir sketches
        hold that k-mer size at that scaled value.
        """
        return None

    def compare(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        threads: int,
        logger: logging.Logger,
        *,
        sketches: Mapping[Path, Path] | None = None,
    ) -> CompareResult:
        """Optional capability: all-vs-all comparison for ``repgenr glance``.

        Adapters that can produce a pairwise similarity table override this;
        the default signals the capability is absent. ``sketches`` is passed
        only to an adapter whose :meth:`sketch_request` returns a pair.
        """
        raise NotImplementedError(
            f"Dereplicator '{self.capabilities.name}' does not support glance "
            "comparisons (no compare() implementation)."
        )

    @abstractmethod
    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        """Cluster ``genomes`` and return representatives + membership."""
        raise NotImplementedError


def compare_supporters(reg: Registry[Dereplicator] | None = None) -> list[str]:
    """Registered, loadable dereplicators that implement ``compare`` (glance)."""
    reg = registry if reg is None else reg
    return sorted(
        name
        for name in reg.names()
        if not reg.is_broken(name)
        # A third-party adapter need not subclass Dereplicator; without a
        # compare attribute it has no comparison support.
        and getattr(reg.get(name), "compare", Dereplicator.compare) is not Dereplicator.compare
    )
