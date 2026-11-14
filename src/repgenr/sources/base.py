"""Selection-source interface: where the metadata stage finds a taxon's genomes.

A source answers one question: which genomes belong to the target, and which
genome roots them. The metadata stage owns everything after that (manifest,
``selection.tsv``, summary, provenance), so every source yields the same
hand-off. Sources are discovered through the ``repgenr.sources`` entry-point
group, like the tool families; ``metadata --source <name>`` selects one.
"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass

from ..core.context import WorkdirContext
from ..core.errors import UserInputError
from ..core.manifest import GenomeRecord
from ..core.outgroup_choice import (
    DEFAULT_OUTGROUP_MODE,
    Choice,
)
from ..core.plugins import Registry
from ..core.selection_summary import (
    GenomeDetails,
)

Selection = tuple[list[GenomeRecord], GenomeRecord | None, dict[str, GenomeDetails]]


@dataclass
class MetadataParams:
    dataset: str  # all | rep
    level: str  # family | genus | species
    release: str | None = None  # required for tsv source
    version: str | None = None  # bac120 | ar53; required for tsv source
    source: str = "tsv"  # tsv | api
    target_family: str | None = None
    target_genus: str | None = None
    # One species epithet, or several separated by commas (all of the target
    # genus); see target_species_list.
    target_species: str | None = None
    outgroup_accession: str | None = None
    # Select no outgroup: selection.tsv holds the ingroup only and no
    # outgroup is downloaded (root the tree with phylo --root midpoint).
    no_outgroup: bool = False
    # How an automatic outgroup is chosen: sister (the nearest group in the
    # GTDB reference tree), quality (the parent taxon) or first (table order).
    outgroup_mode: str = DEFAULT_OUTGROUP_MODE
    # A local GTDB reference tree (bac120_r232.tree) instead of a download.
    gtdb_tree_path: str | None = None
    metadata_path: str | None = None
    nodownload: bool = False
    limit: int | None = None
    # Discard genomes another entry path appended (assemble --append) instead
    # of refusing to overwrite the selection that holds them.
    drop_foreign: bool = False
    # Options of the NCBI source (see SOURCE_OPTIONS): which assemblies NCBI
    # is asked for. refseq (GCF_), genbank (GCA_) or all (one per assembly,
    # GCF_ where both exist).
    assembly_source: str = "refseq"
    # complete, chromosome, scaffold, contig; several separated by commas.
    assembly_level: str | None = None
    # Keep assemblies NCBI flags as atypical (excluded by default).
    include_atypical: bool = False
    # Leave out metagenome-assembled genomes.
    exclude_mag: bool = False
    # Only assemblies released on or before this date (YYYY-MM-DD).
    released_before: str | None = None
    # Option of the GTDB table source: the taxonomy the target is named in
    # and the genomes are named by. gtdb, or ncbi (the table's ncbi_taxonomy
    # column), which keeps GTDB's release, quality values and representatives.
    taxonomy: str = "gtdb"


# Options only some sources read: field -> (flag, the value that means "not given").
SOURCE_OPTIONS: dict[str, tuple[str, object]] = {
    "assembly_source": ("--assembly-source", "refseq"),
    "assembly_level": ("--assembly-level", None),
    "include_atypical": ("--include-atypical", False),
    "exclude_mag": ("--exclude-mag", False),
    "released_before": ("--released-before", None),
    "taxonomy": ("--taxonomy", "gtdb"),
}

TAXONOMIES = ("gtdb", "ncbi")
# The token of a rank NCBI gives no name for.
UNCLASSIFIED = "unclassified"
_UNSAFE = re.compile(r"[^A-Za-z0-9-]+")


def taxon_token(value: str) -> str:
    """A filename token from an NCBI name: NCBI names hold characters GTDB
    names do not ('sp. FSC1006', '[Clostridium]', 'str. SCHU S4 substr. FTS-634/635')."""
    return re.sub(r"-{2,}", "-", _UNSAFE.sub("-", value)).strip("-") or UNCLASSIFIED


def species_token(genus_name: str, species_name: str) -> str:
    """The species as one filename token: the name without its genus, the
    remaining words joined by hyphens ('Francisella tularensis' ->
    'tularensis', as the GTDB sources write it; 'uncultured Francisella sp.'
    -> 'uncultured-sp'; 'Candidatus Francisella endociliophora' ->
    'Candidatus-endociliophora')."""
    words = [w for w in species_name.split() if w != genus_name]
    return taxon_token("-".join(words))


def species_or_none(params: MetadataParams) -> list[str | None]:
    """The species to select one by one; [None] for a genus or family target."""
    return list(target_species_list(params.target_species)) or [None]


def target_species_list(value: str | None) -> list[str]:
    """The species epithets of a -ts value: one, or several separated by
    commas ("tularensis,novicida"), in the order given, without repeats."""
    if not value:
        return []
    out: list[str] = []
    for part in value.split(","):
        name = part.strip()
        if name and name not in out:
            out.append(name)
    return out


def log_choice(choice: Choice, logger) -> None:
    c = choice.chosen
    where = f", path length {c.distance:.3f} to the target in the GTDB tree" if c.distance else ""
    logger.info(
        "Outgroup %s (%s): %s; completeness %s, contamination %s, N50 %s%s "
        "(%d candidates in outgroup_choice.tsv)",
        c.accession,
        c.species,
        choice.reason,
        c.completeness,
        c.contamination,
        c.n50,
        where,
        len(choice.considered),
    )


@dataclass(frozen=True)
class LimitCandidate:
    """One genome competing for a place under ``--limit``."""

    accession: str
    species: str
    is_rep: bool
    completeness: float | None
    contamination: float | None

    def score(self) -> float | None:
        if self.completeness is None or self.contamination is None:
            return None
        from ..stages.derep_keeper import quality_score

        return quality_score(self.completeness, self.contamination)


def stratified_limit(
    candidates: Sequence[LimitCandidate], limit: int | None
) -> list[LimitCandidate]:
    """Choose up to ``limit`` genomes: the best of every species first, then each
    species' next best, and so on (round-robin), so a heavily sequenced species
    cannot fill the cap and file order plays no part.

    Within a species genomes rank by the CheckM part of the keeper score (no
    N50 exists before download), unscored genomes last, then the GTDB
    species-representative flag, then accession, so the choice is
    deterministic. Species are visited in name order. With no
    limit every candidate is returned in that ranked order.
    """

    def rank(c: LimitCandidate) -> tuple:
        score = c.score()
        return (score is None, -(score or 0.0), not c.is_rep, c.accession)

    by_species: dict[str, list[LimitCandidate]] = {}
    for c in candidates:
        by_species.setdefault(c.species, []).append(c)
    queues = [sorted(members, key=rank) for _species, members in sorted(by_species.items())]

    kept: list[LimitCandidate] = []
    total = len(candidates)
    target = total if limit is None else min(limit, total)
    depth = 0
    while len(kept) < target:
        for queue in queues:
            if depth < len(queue):
                kept.append(queue[depth])
                if len(kept) == target:
                    break
        depth += 1
    return kept


def log_limit(
    limit: int | None, candidates: Sequence[LimitCandidate], kept: Sequence, logger
) -> None:
    if limit and len(kept) < len(candidates):
        species = len({c.species for c in candidates})
        logger.info(
            "Limit %d: kept %d of %d candidate genomes across %d species "
            "(best assembly quality per species first)",
            limit,
            len(kept),
            len(candidates),
            species,
        )


def refuse_outgroup_in_selection(accession: str, selected) -> None:
    if accession in selected:
        raise UserInputError(
            f"--outgroup-accession {accession} is part of the selection itself; "
            "an outgroup must lie outside it."
        )


def refuse_outgroup_in_target(accession: str, level: str, target: str, value: str) -> None:
    """A named outgroup obeys the automatic rule: it lies outside the target taxon.

    Under --dataset rep or --limit a target genome can be left out of the
    selection; it is still not an outgroup.
    """
    if value and value == target:
        raise UserInputError(
            f"--outgroup-accession {accession} lies inside the target {level} {target}; "
            "an outgroup must lie outside the target taxon."
        )


class GenomeSource(ABC):
    """One place the metadata stage can select genomes from."""

    # The --source value.
    name: str
    # ``source`` of the genomes in the manifest.
    manifest_source: str
    # ``tool`` of the stage record (a web service or a table; no binary).
    record_tool: str
    # ``source`` in selection_stats.json.
    stats_source: str
    # Files the source leaves in the workdir for ``census``; a run of another
    # source removes them.
    workdir_files: tuple[str, ...] = ()
    # The fields of SOURCE_OPTIONS this source reads; the others are refused.
    options: frozenset[str] = frozenset()
    # Versions of the binaries a selection ran, set by prepare(); empty for a
    # source that needs none.
    tool_versions: dict[str, str] = {}

    def validate(self, params: MetadataParams) -> None:
        """Refuse parameter values this source cannot use."""
        return None

    def prepare(self, ctx: WorkdirContext, params: MetadataParams, logger) -> dict[str, object]:
        """Run before the selection; the returned keys join the stage record."""
        return {}

    @abstractmethod
    def select(
        self,
        ctx: WorkdirContext,
        params: MetadataParams,
        logger,
        counts: dict[str, int],
    ) -> Selection:
        """The selected genomes, the outgroup (None under --no-outgroup) and
        the per-genome details of ``metadata_selected.tsv``.
        ``counts["candidates"]`` is set to the matching genomes before --limit."""


registry: Registry[GenomeSource] = Registry("repgenr.sources")


def record_from_tax(
    accession: str,
    tax: dict,
    *,
    source: str,
    is_outgroup: bool = False,
    completeness: float | None = None,
    contamination: float | None = None,
    source_representative: bool = False,
) -> GenomeRecord:
    return GenomeRecord(
        accession=accession,
        source=source,
        is_outgroup=is_outgroup,
        family=tax["family"],
        genus=tax["genus"],
        species=tax["species"],
        completeness=completeness,
        contamination=contamination,
        source_representative=source_representative,
    )
