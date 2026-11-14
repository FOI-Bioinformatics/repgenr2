"""NCBI genomes as a selection source (``--source ncbi``).

Selects the assemblies NCBI holds for a taxon of NCBI Taxonomy, through the
``datasets`` command line client: RefSeq (``GCF_``) by default, GenBank or
both on request. It reaches genomes a GTDB release does not contain yet, and
names them as NCBI does.

What it does not offer, compared with the GTDB sources:

* no release to pin: the answer is NCBI's at the time of the query, whose
  date is recorded (``--released-before`` narrows it to older assemblies);
* no reference tree, so the outgroup is chosen by assembly quality among the
  reference genomes of the parent taxon, not as the sister group;
* completeness and contamination are NCBI's CheckM values (taxon-specific
  marker sets), which are absent for some assemblies and are not comparable
  with the CheckM2 values of GTDB.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from datetime import UTC, date, datetime
from pathlib import Path

from ..core import ncbi_datasets
from ..core.contracts import atomic_path
from ..core.errors import UserInputError, WorkdirError
from ..core.http import NCBI_DATASETS_URL, require_reachable
from ..core.manifest import GenomeRecord
from ..core.ncbi_datasets import Taxon
from ..core.outgroup_choice import Candidate, Choice, choose_by_quality, write_choice
from ..core.plugins import preflight
from ..core.selection_summary import GenomeDetails
from .base import (
    UNCLASSIFIED,
    GenomeSource,
    LimitCandidate,
    MetadataParams,
    Selection,
    log_choice,
    log_limit,
    record_from_tax,
    refuse_outgroup_in_selection,
    species_or_none,
    species_token,
    stratified_limit,
    taxon_token,
)

MANIFEST_SOURCE = "ncbi"
# The candidates NCBI returned for the target, before --limit: what
# `repgenr census -wd` counts for a selection of this source.
NCBI_GENOMES = "ncbi_genomes.tsv"
_GENOMES_COLUMNS = (
    "accession",
    "family",
    "genus",
    "species",
    "reference_genome",
    "family_name",
    "genus_name",
    "species_name",
)

# --assembly-source value -> the spelling `datasets` takes.
ASSEMBLY_SOURCES = {"refseq": "RefSeq", "genbank": "GenBank", "all": "all"}
_SOURCE_LABELS = {"refseq": "RefSeq", "genbank": "GenBank", "all": "RefSeq and GenBank"}
ASSEMBLY_LEVELS = ("complete", "chromosome", "scaffold", "contig")
# The ranks tried, in order, for the taxon an automatic outgroup comes from.
_PARENT_RANKS = {
    "species": ("genus", "family", "order"),
    "genus": ("family", "order", "class"),
    "family": ("order", "class", "phylum"),
}
_LINEAGE_RANKS = (
    ("d", "domain"),
    ("p", "phylum"),
    ("c", "class"),
    ("o", "order"),
    ("f", "family"),
    ("g", "genus"),
    ("s", "species"),
)


@dataclass(frozen=True)
class NcbiGenome:
    """What the selection keeps of one line of the genome report."""

    accession: str
    paired_accession: str
    tax_id: str
    organism: str
    strain: str
    completeness: float | None
    contamination: float | None
    n50: int | None
    contigs: int | None
    reference: bool
    type_material: str
    category: str


def _number(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _integer(value) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def genome_from_report(row: dict) -> NcbiGenome | None:
    """One line of ``datasets summary genome`` as an :class:`NcbiGenome`."""
    accession = str(row.get("accession") or "")
    if not accession:
        return None
    info = row.get("assembly_info") or {}
    organism = row.get("organism") or {}
    names = organism.get("infraspecific_names") or {}
    checkm = row.get("checkm_info") or {}
    stats = row.get("assembly_stats") or {}
    type_material = row.get("type_material") or info.get("type_material") or {}
    notes = info.get("genome_notes") or []
    return NcbiGenome(
        accession=accession,
        paired_accession=str(row.get("paired_accession") or ""),
        tax_id=str(organism.get("tax_id") or ""),
        organism=str(organism.get("organism_name") or ""),
        strain=str(names.get("strain") or names.get("isolate") or ""),
        completeness=_number(checkm.get("completeness")),
        contamination=_number(checkm.get("contamination")),
        n50=_integer(stats.get("contig_n50")),
        contigs=_integer(stats.get("number_of_contigs")),
        reference=str(info.get("refseq_category") or "").lower() == "reference genome",
        type_material=str(type_material.get("type_display_text") or ""),
        # "derived from metagenome" and the like; empty for an isolate.
        category="; ".join(str(n) for n in notes),
    )


def tax_of(taxon: Taxon | None) -> dict[str, str]:
    """The manifest tax dict (filename tokens and display names) of a genome
    whose taxid resolved to ``taxon``; every level is 'unclassified' without one."""
    ranks = taxon.ranks if taxon is not None else {}
    family_name = ranks.get("family", "")
    genus_name = ranks.get("genus", "")
    species_name = ranks.get("species", "")
    return {
        "family": taxon_token(family_name),
        "genus": taxon_token(genus_name),
        "species": species_token(genus_name, species_name),
        "family_name": family_name,
        "genus_name": genus_name,
        "species_name": species_name,
    }


def lineage_string(taxon: Taxon | None) -> str:
    """The classification in the GTDB table's notation (d__...;s__...)."""
    if taxon is None:
        return ""
    return ";".join(f"{prefix}__{taxon.ranks.get(rank, '')}" for prefix, rank in _LINEAGE_RANKS)


def collapse_pairs(genomes: dict[str, NcbiGenome]) -> dict[str, NcbiGenome]:
    """One entry per assembly: the GenBank copy is dropped where its RefSeq
    twin is present, so an assembly held by both is counted once, as GCF_."""
    return {
        acc: g
        for acc, g in genomes.items()
        if not (acc.startswith("GCA_") and g.paired_accession in genomes)
    }


def target_names(params: MetadataParams) -> list[str]:
    """The NCBI Taxonomy names to query, one per species given."""

    def capital(name: str) -> str:
        name = name.strip()
        return name[:1].upper() + name[1:]

    if params.level == "family":
        return [capital(params.target_family or "")]
    genus = capital(params.target_genus or "")
    if params.level == "genus":
        return [genus]
    return [f"{genus} {species.strip()}" for species in species_or_none(params) if species]


def genome_flags(params: MetadataParams, *, reference: bool) -> list[str]:
    """The ``datasets summary genome`` filters of a selection."""
    flags = [
        "--assembly-source",
        ASSEMBLY_SOURCES[params.assembly_source],
        "--assembly-version",
        "latest",
    ]
    if reference:
        flags.append("--reference")
    if not params.include_atypical:
        flags.append("--exclude-atypical")
    if params.exclude_mag:
        flags += ["--mag", "exclude"]
    if params.assembly_level:
        flags += ["--assembly-level", params.assembly_level]
    if params.released_before:
        flags += ["--released-before", params.released_before]
    return flags


def _genomes(tax_id: str, flags: list[str], scratch: Path, logger) -> dict[str, NcbiGenome]:
    report = scratch / f"genomes_{tax_id}.jsonl"
    try:
        ncbi_datasets.summary(["genome", "taxon", tax_id, *flags], report, logger)
        found = (genome_from_report(row) for row in ncbi_datasets.read_json_lines(report))
        return {g.accession: g for g in found if g is not None}
    finally:
        report.unlink(missing_ok=True)


def _genome_by_accession(accession: str, scratch: Path, logger) -> NcbiGenome | None:
    report = scratch / "genome_accession.jsonl"
    try:
        ncbi_datasets.summary(["genome", "accession", accession], report, logger)
        for row in ncbi_datasets.read_json_lines(report):
            genome = genome_from_report(row)
            if genome is not None and genome.accession == accession:
                return genome
        return None
    finally:
        report.unlink(missing_ok=True)


def write_genomes(path: Path, genomes: dict[str, NcbiGenome], tax: dict[str, dict]) -> None:
    """Keep NCBI's answer (accession, taxonomy tokens, reference flag)."""
    with atomic_path(path) as tmp, open(tmp, "w", encoding="utf-8", newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(_GENOMES_COLUMNS)
        for acc, genome in genomes.items():
            t = tax[acc]
            writer.writerow(
                [
                    acc,
                    t["family"],
                    t["genus"],
                    t["species"],
                    1 if genome.reference else 0,
                    t["family_name"],
                    t["genus_name"],
                    t["species_name"],
                ]
            )


def read_genomes(path: Path) -> list[dict]:
    """Rows of :func:`write_genomes` as dicts with a boolean ``is_rep``."""
    with open(path, encoding="utf-8", newline="") as fo:
        return [
            {
                "accession": row["accession"],
                "tax": {
                    key: row.get(key) or ""
                    for key in (
                        "family",
                        "genus",
                        "species",
                        "family_name",
                        "genus_name",
                        "species_name",
                    )
                },
                "is_rep": (row.get("reference_genome") or "0").strip() in ("1", "true"),
            }
            for row in csv.DictReader(fo, delimiter="\t")
        ]


def _details(genome: NcbiGenome, taxon: Taxon | None) -> GenomeDetails:
    lineage = lineage_string(taxon)
    return GenomeDetails(
        # The GenBank twin of a RefSeq assembly; the accession itself otherwise.
        accession_ncbi=(
            genome.paired_accession if genome.accession.startswith("GCF_") else genome.accession
        ),
        gtdb_taxonomy="",
        ncbi_taxonomy=lineage,
        # The organism name as submitted carries the strain-level name the
        # classification stops short of.
        ncbi_taxonomy_unfiltered=f"{lineage};x__{genome.organism}" if lineage else "",
        ncbi_strain_identifiers=genome.strain,
    )


def _record(genome: NcbiGenome, tax: dict, *, is_outgroup: bool = False) -> GenomeRecord:
    return record_from_tax(
        genome.accession,
        tax,
        source=MANIFEST_SOURCE,
        is_outgroup=is_outgroup,
        completeness=genome.completeness,
        contamination=genome.contamination,
        # NCBI's reference genome of the species, for --keeper source.
        source_representative=genome.reference,
    )


def _candidate(genome: NcbiGenome) -> Candidate:
    return Candidate(
        accession=genome.accession,
        species=genome.organism,
        completeness=genome.completeness,
        contamination=genome.contamination,
        n50=genome.n50,
        contigs=genome.contigs,
        category=genome.category,
        type_material=genome.type_material,
    )


def candidates(
    params: MetadataParams, scratch: Path, logger
) -> tuple[list[Taxon], dict[str, NcbiGenome], dict[str, Taxon], dict[str, dict[str, str]]]:
    """What NCBI holds for the target of ``params``: the resolved target
    taxa, the genomes by accession, the taxonomy node of each genome's taxid
    and the tax dict of each genome. Shared by the selection and by
    ``repgenr census --source ncbi``."""
    targets = [
        ncbi_datasets.resolve_name(name, params.level, scratch, logger)
        for name in target_names(params)
    ]
    flags = genome_flags(params, reference=params.dataset == "rep")
    genomes: dict[str, NcbiGenome] = {}
    for target in targets:
        logger.info(
            "Querying NCBI for the %s genomes of %s (taxid %s)",
            _SOURCE_LABELS[params.assembly_source],
            target.name,
            target.tax_id,
        )
        found = _genomes(target.tax_id, flags, scratch, logger)
        if not found:
            raise UserInputError(
                f"NCBI returned no genomes for {target.name} (taxid {target.tax_id}) "
                f"with {' '.join(flags)}. Relax the filters, or pass "
                "--assembly-source genbank for assemblies RefSeq does not hold."
            )
        for acc, genome in found.items():
            genomes.setdefault(acc, genome)
    if params.assembly_source == "all":
        genomes = collapse_pairs(genomes)

    taxa = ncbi_datasets.taxonomy({g.tax_id for g in genomes.values()}, scratch, logger)
    unresolved = sorted({g.tax_id for g in genomes.values() if g.tax_id not in taxa})
    if unresolved:
        logger.warning(
            "NCBI Taxonomy gave no classification for %d taxid(s) (%s); their genomes "
            "are named '%s'.",
            len(unresolved),
            ", ".join(unresolved[:5]),
            UNCLASSIFIED,
        )
    tax = {acc: tax_of(taxa.get(g.tax_id)) for acc, g in genomes.items()}
    return targets, genomes, taxa, tax


class NcbiSource(GenomeSource):
    name = "ncbi"
    manifest_source = MANIFEST_SOURCE
    record_tool = "ncbi-datasets"
    stats_source = "ncbi"
    workdir_files = (NCBI_GENOMES,)
    options = frozenset(
        {"assembly_source", "assembly_level", "include_atypical", "exclude_mag", "released_before"}
    )

    def validate(self, params: MetadataParams) -> None:
        if params.assembly_source not in ASSEMBLY_SOURCES:
            raise UserInputError(
                f"--assembly-source must be one of {', '.join(ASSEMBLY_SOURCES)}; "
                f"got {params.assembly_source!r}."
            )
        if params.assembly_level:
            unknown = [
                level
                for level in params.assembly_level.split(",")
                if level.strip() not in ASSEMBLY_LEVELS
            ]
            if unknown:
                raise UserInputError(
                    f"--assembly-level takes {', '.join(ASSEMBLY_LEVELS)} (comma-separated); "
                    f"got {', '.join(unknown)}."
                )
        if params.released_before:
            try:
                date.fromisoformat(params.released_before)
            except ValueError as exc:
                raise UserInputError(
                    f"--released-before takes a date as YYYY-MM-DD; got {params.released_before!r}."
                ) from exc

    def prepare(self, ctx, params, logger) -> dict[str, object]:
        self.tool_versions = preflight(ncbi_datasets.DATASETS_CAPS)
        require_reachable(NCBI_DATASETS_URL, what="NCBI datasets")
        ignored = [
            flag
            for flag, value in (
                ("--release", params.release),
                ("--gtdb-version", params.version),
                ("--metadata-path", params.metadata_path),
                ("--nodownload", params.nodownload),
                ("--gtdb-tree-path", params.gtdb_tree_path),
            )
            if value
        ]
        if ignored:
            logger.warning("--source ncbi reads NCBI, not GTDB; ignoring %s.", ", ".join(ignored))
        return {
            # NCBI has no release to name; the date of the query stands in.
            "query_date": datetime.now(UTC).isoformat(timespec="seconds"),
            "assembly_source": params.assembly_source,
            "assembly_level": params.assembly_level,
            "include_atypical": params.include_atypical,
            "exclude_mag": params.exclude_mag,
            "released_before": params.released_before,
        }

    def select(self, ctx, params, logger, counts) -> Selection:
        scratch = ctx.scratch_dir / "ncbi_metadata"
        scratch.mkdir(parents=True, exist_ok=True)

        targets, genomes, taxa, tax = candidates(params, scratch, logger)
        write_genomes(ctx.workdir / NCBI_GENOMES, genomes, tax)

        counts["candidates"] = len(genomes)
        if params.limit and len(genomes) > params.limit:
            ranked = [
                LimitCandidate(
                    accession=acc,
                    species=tax[acc]["species"],
                    is_rep=g.reference,
                    completeness=g.completeness,
                    contamination=g.contamination,
                )
                for acc, g in genomes.items()
            ]
            kept = stratified_limit(ranked, params.limit)
            log_limit(params.limit, ranked, kept, logger)
            genomes = {c.accession: genomes[c.accession] for c in kept}

        records = [_record(g, tax[acc]) for acc, g in genomes.items()]
        details = {acc: _details(g, taxa.get(g.tax_id)) for acc, g in genomes.items()}
        if params.no_outgroup:
            return records, None, details

        outgroup, outgroup_taxon = self._outgroup(ctx, params, targets, genomes, scratch, logger)
        if params.assembly_source == "genbank" and outgroup.accession.startswith("GCF_"):
            # Reference genomes are listed in RefSeq only; a GenBank selection
            # takes the GenBank copy of the chosen one, so it holds GCA_ alone.
            if outgroup.paired_accession:
                outgroup = replace(
                    outgroup,
                    accession=outgroup.paired_accession,
                    paired_accession=outgroup.accession,
                )
        try:
            scratch.rmdir()  # reports are removed as they are read
        except OSError:
            pass
        details[outgroup.accession] = _details(outgroup, outgroup_taxon)
        return records, _record(outgroup, tax_of(outgroup_taxon), is_outgroup=True), details

    def _outgroup(
        self, ctx, params, targets: list[Taxon], selected: dict[str, NcbiGenome], scratch, logger
    ) -> tuple[NcbiGenome, Taxon | None]:
        target_ids = {t.tax_id for t in targets}

        def taxon_of(genome: NcbiGenome) -> Taxon | None:
            return ncbi_datasets.taxonomy([genome.tax_id], scratch, logger).get(genome.tax_id)

        if params.outgroup_accession:
            refuse_outgroup_in_selection(params.outgroup_accession, selected)
            genome = _genome_by_accession(params.outgroup_accession, scratch, logger)
            if genome is None:
                raise UserInputError(
                    f"Outgroup accession {params.outgroup_accession} not found at NCBI."
                )
            taxon = taxon_of(genome)
            if taxon is not None and taxon.lineage & target_ids:
                inside = next(t for t in targets if t.tax_id in taxon.lineage)
                raise UserInputError(
                    f"--outgroup-accession {params.outgroup_accession} lies inside the target "
                    f"{params.level} {inside.name}; an outgroup must lie outside the target taxon."
                )
            return genome, taxon

        if params.outgroup_mode == "sister":
            logger.warning(
                "--source ncbi has no reference tree for the sister-group choice of the "
                "outgroup; choosing by quality among the reference genomes of the parent taxon."
            )
        # Reference genomes: one per species, so the pool stays small and its
        # members are the assemblies NCBI itself prefers. They exist in RefSeq
        # only, whatever --assembly-source selects.
        pool_flags = [
            "--assembly-source",
            "RefSeq",
            "--assembly-version",
            "latest",
            "--reference",
            "--exclude-atypical",
        ]
        for rank in _PARENT_RANKS[params.level]:
            parent_id = targets[0].rank_ids.get(rank)
            if not parent_id:
                continue
            parent_name = targets[0].ranks.get(rank, parent_id)
            logger.info("Selecting the outgroup from the reference genomes of %s", parent_name)
            pool = _genomes(parent_id, pool_flags, scratch, logger)
            taxa = ncbi_datasets.taxonomy({g.tax_id for g in pool.values()}, scratch, logger)
            outside = {
                acc: g
                for acc, g in pool.items()
                if acc not in selected
                and g.paired_accession not in selected
                and g.tax_id in taxa
                and not (taxa[g.tax_id].lineage & target_ids)
            }
            if not outside:
                logger.info("No reference genome of %s lies outside the target", parent_name)
                continue
            if params.outgroup_mode == "first":
                genome = outside[min(outside)]
            else:
                choice: Choice | None = choose_by_quality([_candidate(g) for g in outside.values()])
                if choice is None:
                    continue
                write_choice(ctx.workdir, choice)
                log_choice(choice, logger)
                genome = outside[choice.chosen.accession]
            return genome, taxa[genome.tax_id]
        raise WorkdirError(
            "Could not determine an outgroup from NCBI; specify --outgroup-accession, or "
            "select none with --no-outgroup."
        )
