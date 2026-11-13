"""census: count the genera, species and samples under a taxon.

A read-only query in three modes that share one table shape:

* taxon (no working directory): the GTDB genomes of a family or genus, from
  the GTDB API or a GTDB metadata table, optionally with the ENA whole-genome
  sequencing runs under the same taxon; or the NCBI Virus sequence records of
  a viral taxon (metadata only, no sequences).
* candidates (``-wd`` after metadata or vmetadata): the genomes or sequence
  records the entry stage found, before genome or vgenome selected from them.
* selection (``-wd`` after genome, vgenome, ingest or assemble): the genomes
  in ``selection.tsv``, split by manifest source, with the dereplication
  clusters when ``derep/clusters.tsv`` exists and the candidates beside them.

For a family the rows are its genera; for a genus they are its species. The
command writes nothing into a working directory and records no stage. A
GTDB table that ``--source table`` downloads is kept in a cache directory
(``$REPGENR_CACHE_DIR/gtdb``, by default ``~/.cache/repgenr/gtdb``).
"""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import tempfile
from collections.abc import Container, Iterable
from dataclasses import dataclass, field, replace
from pathlib import Path
from statistics import median
from typing import Any

from ..core.contracts import (
    CLUSTERS_TSV,
    READS_TSV,
    SELECTION_TSV,
    SelectionRow,
    atomic_path,
    parse_genome_filename,
    read_clusters,
    read_reads,
    read_selection,
    record_name,
    sanitise_taxon_tokens,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.manifest import VIRAL_SOURCES as VIRAL_MANIFEST_SOURCES

FAMILY = "family"
GENUS = "genus"
MODE_TAXON = "taxon"
MODE_CANDIDATES = "candidates"
MODE_SELECTION = "selection"

BACTERIAL_SOURCES = ("api", "table")
VIRAL_SOURCES = ("ncbi_virus", "bvbrc")
# The metadata stage calls the table source 'tsv'; census accepts both names.
SOURCE_ALIASES = {"tsv": "table"}
SOURCE_CHOICES = (*BACTERIAL_SOURCES, "tsv", *VIRAL_SOURCES)

# Manifest sources in the order their columns appear.
MANIFEST_SOURCES = ("gtdb", "sra", "ncbi_virus", "bvbrc", "local", "local_viral")
# ENA instrument platforms counted in their own column.
_PLATFORM_COLUMNS = {"ILLUMINA": "illumina", "OXFORD_NANOPORE": "ont", "PACBIO_SMRT": "pacbio"}
# Runs on any other platform (BGISEQ, ION_TORRENT, ...) are counted as 'other'.
RUN_COLUMNS = ("runs", "biosamples", "illumina", "ont", "pacbio", "other")

CACHE_ENV = "REPGENR_CACHE_DIR"
_CANDIDATE_LABELS = {
    "gtdb-table": "the genomes of the GTDB table under the metadata target",
    "gtdb-api": "the genomes the GTDB API returned for the metadata target",
    "ncbi_virus": "the NCBI Virus sequence records of vmetadata",
    "bvbrc": "the BV-BRC records of vmetadata",
    "reads": "the sequencing runs of reads.tsv",
}
_ENTRY_SELECTION_STAGES = ("genome", "vgenome", "ingest", "assemble")
SCHEMA = "repgenr.census/1"


@dataclass
class CensusParams:
    workdir: Path | None = None
    target_family: str | None = None
    target_genus: str | None = None
    viral: bool = False
    target: str | None = None
    source: str | None = None
    release: str | None = None
    gtdb_version: str | None = None
    metadata_path: str | None = None
    runs: bool = False
    host: str | None = None
    complete_only: bool = False
    released_after: str | None = None


@dataclass
class Census:
    """One census table: the totals line and the rows below it.

    ``rank`` is the rank of the counted taxon (family or genus); the rows are
    one rank below it. ``columns`` orders the row keys for the console table
    and the TSV.
    """

    taxon: str
    mode: str
    source: str
    rank: str
    columns: list[str]
    rows: list[dict[str, Any]]
    totals: dict[str, Any]
    notes: list[str] = field(default_factory=list)
    # The taxon the runs were counted under when it differs from ``taxon``
    # ('NCBI genus Bacillus' for the GTDB genus Bacillus_A); empty otherwise.
    runs_taxon: str = ""

    @property
    def row_rank(self) -> str:
        return "genus" if self.rank == FAMILY else "species"

    def header(self) -> str:
        """'<Taxon>: N genera, N species, N genomes (N representatives)' and its extensions."""
        t = self.totals
        head = f"{self.taxon}: {_n(t['genera'], 'genus', 'genera')}, {t['species']} species"
        if "sequences" in t:
            head += f", {_n(t['sequences'], 'sequence')}"
            extra = []
            if t.get("complete") is not None:
                extra.append(f"{t['complete']} complete")
            if t.get("isolates") is not None:
                extra.append(_n(t["isolates"], "isolate"))
        else:
            head += f", {_n(t['genomes'], 'genome')}"
            extra = [_n(t["representatives"], "representative")]
        if extra:
            head += f" ({', '.join(extra)})"
        for key in ("candidates", "runs", "biosamples", "clusters"):
            if t.get(key) is not None:
                head += f", {_n(t[key], key[:-1])}"
                if key == "runs" and self.runs_taxon:
                    head += f" of {self.runs_taxon}"
        sources = t.get("sources") or {}
        if len(sources) > 1:
            head += "; " + ", ".join(f"{n} {s}" for s, n in sources.items())
        if t.get("outgroup"):
            head += f"; outgroup {t['outgroup']} excluded"
        return head

    def to_json(self) -> dict[str, Any]:
        from .. import __version__

        return {
            "schema": SCHEMA,
            "repgenr": __version__,
            "taxon": self.taxon,
            "mode": self.mode,
            "source": self.source,
            "rank": self.rank,
            "row_rank": self.row_rank,
            "totals": self.totals,
            "rows": self.rows,
            "notes": self.notes,
            "runs_taxon": self.runs_taxon,
        }


def _n(count: int, singular: str, plural: str | None = None) -> str:
    return f"{count} {singular if count == 1 else (plural or singular + 's')}"


# --- rows -------------------------------------------------------------------


def _norm(name: str | None) -> str:
    """Compare taxon names and tokens alike: case, spaces and underscores ignored."""
    return (name or "").strip().lower().replace("_", "-").replace(" ", "-")


def _label(rank: str, genus: str, species: str, *, binomial: bool = True) -> str:
    """The row name: the genus, or the species.

    A bacterial species token is the epithet alone and is written after its
    genus (``binomial``); a viral species token is a whole name.
    """
    if rank == FAMILY:
        return genus or "(no genus)"
    if not species:
        return f"{genus} (no species)".strip()
    if not binomial or not genus or _norm(species).startswith(_norm(genus)):
        return species
    return f"{genus} {species}"


def _row_key(rank: str, genus: str, species: str) -> tuple[str, str]:
    return (_norm(genus), "") if rank == FAMILY else (_norm(genus), _norm(species))


class _Rows:
    """Rows keyed by genus (family census) or species (genus census).

    Counters add up per row; ``distinct`` values are collected as sets and
    written as their sizes. A family census also counts the distinct species
    of each genus.
    """

    def __init__(self, rank: str, *, binomial: bool = True) -> None:
        self.rank = rank
        self.binomial = binomial
        self._rows: dict[tuple[str, str], dict[str, Any]] = {}
        self._species: dict[tuple[str, str], set[str]] = {}
        self._sets: dict[tuple[str, str], dict[str, set[str]]] = {}
        # Rows whose name is a display name rather than built from tokens.
        self._displayed: set[tuple[str, str]] = set()

    def add(
        self,
        genus: str,
        species: str,
        counts: dict[str, int] | None = None,
        distinct: dict[str, str] | None = None,
        *,
        count_species: bool = True,
        label: str = "",
    ) -> dict[str, Any]:
        """Count into the row of (genus, species) tokens.

        Rows are keyed by the tokens, so GTDB genomes, NCBI-labelled runs and
        candidates join; ``label`` is the name to show for the row (the GTDB
        spelling, e.g. 'Bacillus_A'), and the first label given is kept.
        """
        key = _row_key(self.rank, genus, species)
        row = self._rows.get(key)
        if row is None:
            row = {"name": _label(self.rank, genus, species, binomial=self.binomial)}
            self._rows[key] = row
            self._species[key] = set()
            self._sets[key] = {}
        if label and key not in self._displayed:
            row["name"] = label
            self._displayed.add(key)
        if species and count_species:
            self._species[key].add(_norm(species))
        for column, n in (counts or {}).items():
            row[column] = row.get(column, 0) + n
        for column, value in (distinct or {}).items():
            self._sets[key].setdefault(column, set()).add(value)
        return row

    def keys(self) -> list[tuple[str, str]]:
        return list(self._rows)

    def row(self, key: tuple[str, str]) -> dict[str, Any]:
        return self._rows[key]

    def set(self, key: tuple[str, str], column: str, value: Any) -> None:
        self._rows[key][column] = value

    def finish(
        self,
        columns: list[str],
        defaults: dict[str, Any],
        sort_by: str,
        last: Container[tuple[str, str]] = (),
    ) -> list[dict]:
        """The rows, largest ``sort_by`` first, then by name; ``last`` rows go
        to the end in the same order."""
        ordered = []
        for key, row in self._rows.items():
            if self.rank == FAMILY:
                row["species"] = len(self._species[key])
            for column, values in self._sets[key].items():
                row[column] = len(values)
            out = {c: row.get(c, defaults.get(c)) for c in columns}
            ordered.append((key in last, -(out.get(sort_by) or 0), out))
        ordered.sort(key=lambda t: (t[0], t[1], str(t[2]["name"]).lower()))
        return [out for _last, _count, out in ordered]


def _ncbi_label(name: str) -> str:
    """Mark a row that holds only runs, named by the NCBI taxonomy."""
    if name.endswith(")") and "(" in name:
        return name[:-1] + ", NCBI)"
    return f"{name} (NCBI)"


def _distinct_taxa(items: Iterable[tuple[str, str]]) -> tuple[int, int]:
    """Distinct genera and species among (genus, species) pairs."""
    pairs = list(items)
    genera = {_norm(g) for g, _ in pairs if g}
    species = {(_norm(g), _norm(s)) for g, s in pairs if s}
    return len(genera), len(species)


@dataclass(frozen=True)
class Item:
    """One genome or sequence record, by its taxonomy tokens."""

    accession: str
    family: str
    genus: str
    species: str
    representative: bool = False
    # Names as the source writes them (GTDB 'Bacillus_A', 'Bacillus_A
    # thuringiensis_S'), for display; empty where only tokens are known.
    genus_name: str = ""
    species_name: str = ""
    family_name: str = ""


@dataclass(frozen=True)
class Run:
    """One ENA sequencing run, labelled by the NCBI lineage of its taxid."""

    accession: str
    biosample: str
    platform: str
    genus: str
    species: str
    genus_name: str = ""
    species_name: str = ""


def _display(rank: str, genus_name: str, species_name: str) -> str:
    """The display name of a row: the genus, or the species binomial."""
    return genus_name if rank == FAMILY else species_name


def genome_census(
    items: list[Item],
    *,
    taxon: str,
    rank: str,
    mode: str,
    source: str,
    runs: list[Run] | None = None,
    notes: list[str] | None = None,
    runs_taxon: str = "",
) -> Census:
    """Genomes per genus or species, with the ENA runs beside them when given."""
    rows = _Rows(rank)
    for g in items:
        rows.add(
            g.genus,
            g.species,
            {"genomes": 1, "representatives": int(g.representative)},
            label=_display(rank, g.genus_name, g.species_name),
        )
    with_genomes = set(rows.keys())
    columns = ["name", *(["species"] if rank == FAMILY else []), "genomes", "representatives"]
    genera, species = _distinct_taxa((g.genus, g.species) for g in items)
    totals: dict[str, Any] = {
        "genera": genera,
        "species": species,
        "genomes": len(items),
        "representatives": sum(1 for g in items if g.representative),
    }
    if runs is not None:
        for r in runs:
            counts = {"runs": 1}
            counts[_PLATFORM_COLUMNS.get(r.platform.upper(), "other")] = 1
            # The species column counts GTDB species; runs carry NCBI names.
            rows.add(
                r.genus,
                r.species,
                counts,
                {"biosamples": r.biosample or r.accession},
                count_species=False,
                label=_display(rank, r.genus_name, r.species_name),
            )
        columns += list(RUN_COLUMNS)
        totals["runs"] = len(runs)
        totals["biosamples"] = len({r.biosample or r.accession for r in runs})
    # Rows of runs alone carry an NCBI name with no GTDB genome under it.
    ncbi_only = set(rows.keys()) - with_genomes
    for key in ncbi_only:
        rows.row(key)["name"] = _ncbi_label(rows.row(key)["name"])
    defaults = dict.fromkeys(columns[1:], 0)
    return Census(
        taxon=taxon,
        mode=mode,
        source=source,
        rank=rank,
        columns=columns,
        rows=rows.finish(columns, defaults, "genomes", last=ncbi_only),
        totals=totals,
        notes=notes or [],
        runs_taxon=runs_taxon,
    )


# --- viral rows --------------------------------------------------------------


def virus_census(records: list[Any], *, taxon: str, rank: str, mode: str, source: str) -> Census:
    """NCBI Virus records per genus or species: sequences, complete sequences,
    isolates (the genomes ``vgenome --group-segments`` would form) and whether
    the row holds a segmented virus (two or more distinct segment labels)."""
    from ..viral.selection import _isolate_segment_sets, _segment_labels

    # A standalone logger: the grouping's info lines about repeated segments
    # concern vgenome, not a count, and no shared logger is reconfigured.
    quiet = logging.Logger("repgenr.census.segments", logging.WARNING)
    rows = _Rows(rank, binomial=False)
    members: dict[tuple[str, str], list[Any]] = {}
    for r in records:
        complete = int(r.completeness == "COMPLETE")
        label = getattr(r, "species_name", "") if rank == GENUS else ""
        rows.add(r.genus, r.species, {"sequences": 1, "complete": complete}, label=label)
        members.setdefault(_row_key(rank, r.genus, r.species), []).append(r)
    isolates_total = 0
    for key, recs in members.items():
        groups, singletons = _isolate_segment_sets(recs, quiet)
        isolates = len(groups) + len(singletons)
        isolates_total += isolates
        rows.set(key, "isolates", isolates)
        rows.set(key, "segmented", "yes" if len(_segment_labels(recs)) >= 2 else "no")
    columns = [
        "name",
        *(["species"] if rank == FAMILY else []),
        "sequences",
        "complete",
        "isolates",
        "segmented",
    ]
    genera, species = _distinct_taxa((r.genus, r.species) for r in records)
    return Census(
        taxon=taxon,
        mode=mode,
        source=source,
        rank=rank,
        columns=columns,
        rows=rows.finish(columns, {}, "sequences"),
        totals={
            "genera": genera,
            "species": species,
            "sequences": len(records),
            "complete": sum(1 for r in records if r.completeness == "COMPLETE"),
            "isolates": isolates_total,
        },
    )


def viral_scope(
    records: list[Any], family: str | None, genus: str | None, target: str | None
) -> tuple[str, str, list[Any]]:
    """(rank, taxon name, records) of a viral census.

    -tg and -tf narrow the records to a genus or a family. Otherwise the
    target decides: a genus name gives a genus census, a family name a
    family census; a target at another rank gives a genus census when its
    records hold one genus and a family census when they hold several.
    """

    def spelling(attr: str, wanted: str, recs: list[Any]) -> str:
        return next((getattr(r, attr) for r in recs if _norm(getattr(r, attr)) == wanted), "")

    if genus:
        recs = [r for r in records if _norm(r.genus) == _norm(genus)]
        rank, name = GENUS, spelling("genus", _norm(genus), recs) or genus
    elif family:
        recs = [r for r in records if _norm(r.family) == _norm(family)]
        rank, name = FAMILY, spelling("family", _norm(family), recs) or family
    else:
        wanted = _norm(target)
        recs = list(records)
        genera = {_norm(r.genus) for r in recs}
        families = {_norm(r.family) for r in recs}
        if wanted in genera:
            recs = [r for r in recs if _norm(r.genus) == wanted]
            rank, name = GENUS, spelling("genus", wanted, recs)
        elif wanted in families:
            recs = [r for r in recs if _norm(r.family) == wanted]
            rank, name = FAMILY, spelling("family", wanted, recs)
        else:
            rank = GENUS if len(genera) == 1 else FAMILY
            name = target or ""
    if not recs:
        named = f"genus {genus}" if genus else f"family {family}"
        raise UserInputError(f"No records of the {named} under {target or 'the target'}.")
    return rank, name, recs


# --- BV-BRC -------------------------------------------------------------------

BVBRC_TAXNAMES = "metadata_ncbi_taxnames_data.json"


def bvbrc_items(download_wd: Path) -> tuple[list[Item], dict[str, int]]:
    """One item per BV-BRC record id from the taxnames JSON of vmetadata.

    Each taxon name carries its rank and the set of record ids under it, so a
    record's family, genus and species are the names whose sets hold it. Also
    returns the median sequence length of each taxid (metadata_base.tsv).
    """
    from ..viral.bvbrc import _read_base

    path = download_wd / BVBRC_TAXNAMES
    try:
        taxnames = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise WorkdirError(f"{path} could not be read ({exc}). Rerun vmetadata.") from exc
    ranks: dict[str, dict[str, str]] = {}
    for name, entry in taxnames.items():
        level = entry.get("level")
        if level not in ("family", "genus", "species"):
            continue
        for record_id in entry.get("datasets") or []:
            ranks.setdefault(record_id, {})[level] = name
    items = [
        Item(rid, r.get("family", ""), r.get("genus", ""), r.get("species", ""))
        for rid, r in sorted(ranks.items())
    ]
    base_tsv = download_wd / "metadata_base.tsv"
    medians: dict[str, int] = {}
    if base_tsv.is_file():
        medians = {t: v["seq_med"] for t, v in _read_base(base_tsv).items()}
    return items, medians


def bvbrc_census(
    items: list[Item], medians: dict[str, int], *, taxon: str, rank: str, mode: str
) -> Census:
    """BV-BRC records per genus or species, with the median of the per-taxid
    median lengths. BV-BRC records carry no completeness, isolate or segment
    fields, so those columns are empty."""
    rows = _Rows(rank, binomial=False)
    taxids: dict[tuple[str, str], set[str]] = {}
    for it in items:
        rows.add(it.genus, it.species, {"sequences": 1})
        taxids.setdefault(_row_key(rank, it.genus, it.species), set()).add(
            it.accession.split(".", 1)[0]
        )
    for key, ids in taxids.items():
        values = [medians[t] for t in ids if t in medians]
        rows.set(key, "median_length", int(median(values)) if values else None)
    columns = [
        "name",
        *(["species"] if rank == FAMILY else []),
        "sequences",
        "complete",
        "isolates",
        "segmented",
        "median_length",
    ]
    genera, species = _distinct_taxa((it.genus, it.species) for it in items)
    return Census(
        taxon=taxon,
        mode=mode,
        source="bvbrc",
        rank=rank,
        columns=columns,
        rows=rows.finish(columns, {}, "sequences"),
        totals={
            "genera": genera,
            "species": species,
            "sequences": len(items),
            "complete": None,
            "isolates": None,
        },
        notes=[
            "BV-BRC records carry no completeness, isolate or segment fields; "
            "median_length is the median of the per-taxid median lengths."
        ],
    )


def item_scope(
    items: list[Item], family: str | None, genus: str | None, target: str | None
) -> tuple[str, str, list[Item]]:
    """:func:`viral_scope` for items (the BV-BRC and selection paths)."""
    return viral_scope(items, family, genus, target)


# --- GTDB ---------------------------------------------------------------------


def cache_dir() -> Path:
    """Where ``--source table`` keeps the GTDB metadata tables it downloads."""
    root = os.environ.get(CACHE_ENV)
    if root:
        return Path(root).expanduser() / "gtdb"
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base / "repgenr" / "gtdb"


def _items_from_parsed(parsed: Iterable[dict]) -> list[Item]:
    return [
        Item(
            d["accession"],
            d["tax"]["family"],
            d["tax"]["genus"],
            d["tax"]["species"],
            bool(d.get("is_rep")),
            genus_name=d["tax"].get("genus_name", ""),
            species_name=d["tax"].get("species_name", ""),
            family_name=d["tax"].get("family_name", ""),
        )
        for d in parsed
    ]


def _strip_rank(value: str) -> str:
    """'g__Bacillus_A' -> 'Bacillus_A'."""
    return value.split("__", 1)[1] if "__" in value else value


def _matches_target(
    tax: dict[str, str], family: str | None, genus: str | None, species: str | None
) -> bool:
    """The metadata stage's target rule: species within genus, else genus, else family."""
    if species and genus:
        return _norm(tax["genus"]) == _norm(genus) and _norm(tax["species"]) == _norm(species)
    if genus:
        return _norm(tax["genus"]) == _norm(genus)
    if family:
        return _norm(tax["family"]) == _norm(family)
    return False


def gtdb_api_items(rank: str, name: str, logger: logging.Logger) -> list[Item]:
    """The GTDB genomes of a family or genus from the GTDB API (one request)."""
    from ..sources.base import MetadataParams
    from ..sources.gtdb_api import api_genomes_detail, normalize_api_tax, target_taxon

    mp = MetadataParams(
        dataset="all",
        level=rank,
        source="api",
        target_family=name if rank == FAMILY else None,
        target_genus=name if rank == GENUS else None,
    )
    taxon = target_taxon(mp)
    logger.info("Querying the GTDB API for the genomes of %s", taxon)
    rows = api_genomes_detail(taxon, False, logger)
    if not rows:
        raise UserInputError(f"The GTDB API returned no genomes for {taxon}.")
    items = []
    for row in rows:
        tax = normalize_api_tax(row)
        items.append(
            Item(
                row["gid"],
                tax["family"],
                tax["genus"],
                tax["species"],
                bool(row.get("gtdbIsRep")),
                genus_name=_strip_rank(row.get("gtdbGenus", "")),
                species_name=_strip_rank(row.get("gtdbSpecies", "")),
                family_name=_strip_rank(row.get("gtdbFamily", "")),
            )
        )
    return items


def gtdb_table_items(
    table: Path,
    logger: logging.Logger,
    *,
    family: str | None = None,
    genus: str | None = None,
    species: str | None = None,
) -> list[Item]:
    """The genomes of a GTDB metadata table under the target, through the
    metadata stage's parser."""
    from ..sources.base import MetadataParams
    from ..sources.gtdb_table import parse_metadata

    mp = MetadataParams(dataset="all", level=GENUS if genus else FAMILY)
    parsed = parse_metadata(table, mp, logger)
    return _items_from_parsed(
        d for d in parsed.values() if _matches_target(d["tax"], family, genus, species)
    )


def obtain_table(params: CensusParams, logger: logging.Logger) -> Path:
    """The GTDB table for ``--source table``: --metadata-path, or the release
    table in the cache directory (downloaded once, then reused)."""
    from ..sources.base import MetadataParams
    from ..sources.gtdb_table import GTDB_VERSIONS, RELEASE_RE, obtain_metadata_table

    if params.metadata_path:
        path = Path(params.metadata_path)
        if not path.is_file():
            raise UserInputError(f"--metadata-path not found: {params.metadata_path}")
        return path
    if not params.release or not RELEASE_RE.match(params.release):
        raise UserInputError(
            "--source table needs --release like '232.0' (major.minor), or --metadata-path."
        )
    version = params.gtdb_version or "bac120"
    if version not in GTDB_VERSIONS:
        raise UserInputError(
            f"Invalid --gtdb-version '{version}'. Choose from: {', '.join(GTDB_VERSIONS)}."
        )
    mp = MetadataParams(
        dataset="all",
        level=FAMILY,
        release=params.release,
        version=version,
        source="tsv",
        nodownload=True,
    )
    return obtain_metadata_table(cache_dir(), mp, logger)


# The portal fields a census reads from each run (of the 17 in ena.RUN_FIELDS);
# the search is already restricted to whole-genome sequencing runs.
CENSUS_RUN_FIELDS = ("run_accession", "sample_accession", "tax_id", "instrument_platform")


def ena_runs(name: str, logger: logging.Logger) -> list[Run]:
    """The ENA whole-genome sequencing runs under a taxon, each labelled with
    the genus and species of its taxid (one lineage lookup per distinct taxid)."""
    from ..core import ena
    from .reads import ncbi_taxon_names

    # GTDB taxa are prokaryotes: of two taxa that share the name, take the
    # prokaryote (division PRO), not, say, the stick-insect genus Bacillus.
    hit = ena.resolve_taxon(name, division="PRO")
    logger.info("Resolved %r to %s (%s, taxid %s)", name, hit.scientific_name, hit.rank, hit.taxid)
    query = ena.taxon_query(hit.taxid)
    ena.announce_search(query, hit.scientific_name, logger, fields=CENSUS_RUN_FIELDS)
    rows = ena.to_read_rows(ena.search_runs(query, fields=CENSUS_RUN_FIELDS))
    seen: set[str] = set()
    unique = []
    for row in rows:
        if row.run_accession not in seen:
            seen.add(row.run_accession)
            unique.append(row)
    logger.info("ENA lists %d whole-genome sequencing runs under %s", len(unique), name)
    names = ncbi_taxon_names([r.taxid for r in unique], logger)
    out = []
    for row in unique:
        family, genus_name, species_name = names.get(row.taxid, ("", "", ""))
        _family, genus, species = sanitise_taxon_tokens(family, genus_name, species_name)
        out.append(
            Run(
                row.run_accession,
                row.biosample,
                row.platform,
                genus,
                species,
                genus_name=genus_name,
                species_name=species_name,
            )
        )
    return out


# --- taxon mode -------------------------------------------------------------


def validate(params: CensusParams) -> str:
    """Check the flag combination before any request; return the source.

    Raises UserInputError (exit 2) for a census that names neither a
    workdir nor a taxon, BV-BRC without a workdir, --runs on a viral taxon,
    and taxon-query flags given with -wd.
    """
    source = SOURCE_ALIASES.get(params.source or "", params.source)
    if params.workdir is not None:
        given = [
            flag
            for flag, value in (
                ("--viral", params.viral),
                ("--target", params.target),
                ("--source", params.source),
                ("--release", params.release),
                ("--gtdb-version", params.gtdb_version),
                ("--runs", params.runs),
                ("--host", params.host),
                ("--complete-only", params.complete_only),
                ("--released-after", params.released_after),
            )
            if value
        ]
        if given:
            raise UserInputError(
                f"{', '.join(given)} select a taxon query and do not apply with -wd, whose "
                "record names the source and target; -tf, -tg and --metadata-path do."
            )
        return ""
    if params.viral:
        source = source or "ncbi_virus"
        if source not in VIRAL_SOURCES:
            raise UserInputError(
                f"--source {params.source} is a GTDB source; a viral census uses "
                f"{' or '.join(VIRAL_SOURCES)}."
            )
        if source == "bvbrc":
            raise UserInputError(
                "BV-BRC offers no metadata-only report. Run 'repgenr vmetadata --source "
                "bvbrc' and then 'repgenr census -wd <wd>', or use the NCBI Virus source."
            )
        if params.runs:
            raise UserInputError(
                "--runs counts ENA runs of bacterial taxa; it is not offered for viral taxa."
            )
        if not (params.target or params.target_genus or params.target_family):
            raise UserInputError("A viral census needs --target (e.g. picornaviridae), or -wd.")
        return source
    source = source or "api"
    if source not in BACTERIAL_SOURCES:
        raise UserInputError(
            f"--source {params.source} is a viral source; add --viral, or use "
            f"{' or '.join(BACTERIAL_SOURCES)}."
        )
    if params.target:
        raise UserInputError("--target names a viral taxon; add --viral, or use -tf/-tg.")
    if not (params.target_family or params.target_genus):
        raise UserInputError(
            "Name a taxon (-tf FAMILY or -tg GENUS, or --viral --target TAXON) or a "
            "working directory (-wd)."
        )
    return source


def taxon_census(params: CensusParams, source: str, logger: logging.Logger) -> Census:
    if params.viral:
        from ..viral import ncbi_virus

        target = params.target or params.target_genus or params.target_family or ""
        records = ncbi_virus.summary(
            target,
            complete_only=params.complete_only,
            host=params.host,
            released_after=params.released_after,
            logger=logger,
        )
        if not records:
            raise UserInputError(
                f"NCBI Virus holds no records for '{target}' under the given filters "
                "(--host, --complete-only, --released-after)."
            )
        rank, name, recs = viral_scope(
            records, params.target_family, params.target_genus, params.target
        )
        return virus_census(recs, taxon=name, rank=rank, mode=MODE_TAXON, source="ncbi_virus")

    rank = GENUS if params.target_genus else FAMILY
    name = (params.target_genus if rank == GENUS else params.target_family) or ""
    if source == "api":
        items = gtdb_api_items(rank, name, logger)
        label = "gtdb-api"
    else:
        table = obtain_table(params, logger)
        items = gtdb_table_items(
            table,
            logger,
            family=params.target_family if rank == FAMILY else None,
            genus=params.target_genus,
        )
        if not items:
            raise UserInputError(f"No genomes of the {rank} {name} in {table}.")
        label = "gtdb-table"
    if rank == GENUS and params.target_family:
        items = [g for g in items if _norm(g.family) == _norm(params.target_family)]
        if not items:
            raise UserInputError(
                f"No genomes of the genus {name} in the family {params.target_family}."
            )
    taxon = _gtdb_spelling(items, rank, name)
    notes = []
    runs = None
    runs_taxon = ""
    if params.runs:
        ncbi_name = ncbi_spelling(taxon)
        if ncbi_name != taxon:
            notes.append(
                f"Runs are counted for the NCBI {rank} {ncbi_name}, the GTDB name {taxon} "
                "without its suffix; in a genus census they join its species rows by "
                f"epithet. The NCBI {rank} may hold more than the GTDB {rank}."
            )
        runs = ena_runs(ncbi_name, logger)
        if ncbi_name != taxon:
            runs_taxon = f"NCBI {rank} {ncbi_name}"
        if rank == GENUS and ncbi_name != taxon and items:
            # The genus census counts one GTDB genus: NCBI runs of the
            # unsuffixed genus join its species rows by epithet
            # (Bacillus cereus runs under Bacillus_A cereus).
            genus_token = items[0].genus
            runs = [
                replace(r, genus=genus_token) if _norm(r.genus) == _norm(ncbi_name) else r
                for r in runs
            ]
    if runs is not None:
        notes.append(
            "Runs are grouped by the NCBI taxonomy of their taxid and genomes by GTDB. "
            "A row marked NCBI holds runs only: no GTDB genome carries that name, "
            "and the taxon may be a GTDB taxon under another name. Such rows are "
            "listed last and are not part of the GTDB totals."
        )
    return genome_census(
        items,
        taxon=taxon,
        rank=rank,
        mode=MODE_TAXON,
        source=label,
        runs=runs,
        notes=notes,
        runs_taxon=runs_taxon,
    )


_GTDB_SUFFIX = re.compile(r"_[A-Z]+$")


def ncbi_spelling(name: str) -> str:
    """A GTDB name without its suffix ('Bacillus_A' -> 'Bacillus'), which is
    how the NCBI and ENA taxonomies know the taxon."""
    return _GTDB_SUFFIX.sub("", name)


def _gtdb_spelling(items: list[Item], rank: str, name: str) -> str:
    """The taxon as GTDB writes it ('Bacillus_A'), else as its token, else ``name``."""
    attr = "genus" if rank == GENUS else "family"
    for g in items:
        if _norm(getattr(g, attr)) == _norm(name):
            return getattr(g, f"{attr}_name") or getattr(g, attr)
    return name


# --- workdir modes ------------------------------------------------------------


def workdir_census(params: CensusParams, logger: logging.Logger) -> Census:
    """Mode 2 or 3, chosen by the stages the workdir records."""
    from ..core.config import Config

    wd = params.workdir
    assert wd is not None
    cfg = Config.load(wd)
    stages = cfg.stages
    if any(s in stages for s in _ENTRY_SELECTION_STAGES) and (wd / SELECTION_TSV).is_file():
        return selection_census(wd, cfg, params, logger)
    if "metadata" in stages:
        found = gtdb_candidates(wd, cfg, params.metadata_path, logger)
        notes: list[str] = []
        if found is None:
            items = [_item(r) for r in _ingroup(wd)[0]]
            label = "selection"
            notes.append(
                "No GTDB table or API answer is in the workdir; the candidates are the "
                "genomes of selection.tsv. Pass --metadata-path to count the table."
            )
        else:
            items, label = found
        rank, name, items = _bacterial_scope(items, cfg.stages["metadata"].params, params)
        return genome_census(
            items, taxon=name, rank=rank, mode=MODE_CANDIDATES, source=label, notes=notes
        )
    if "vmetadata" in stages:
        return viral_candidates_census(wd, cfg, params)
    raise UserInputError(
        f"{wd} records none of metadata, vmetadata, genome, vgenome, ingest or assemble; "
        "census counts the genomes an entry stage found or selected."
    )


def _bacterial_scope(
    items: list[Item], record: dict, params: CensusParams
) -> tuple[str, str, list[Item]]:
    """Rank and name of a bacterial workdir census: -tg/-tf narrow, else the
    recorded target (a genus or species target gives a genus census)."""
    if params.target_genus or params.target_family:
        return item_scope(items, params.target_family, params.target_genus, None)
    genus, family, species = (
        record.get("target_genus"),
        record.get("target_family"),
        record.get("target_species"),
    )
    if genus:
        name = _gtdb_spelling(items, GENUS, genus)
        if species and record.get("level") == "species":
            name = f"{name} {species}"
        return GENUS, name, items
    if family:
        return FAMILY, _gtdb_spelling(items, FAMILY, family), items
    return item_scope(items, None, None, None)


def _find_table(wd: Path, cfg: Any, metadata_path: str | None) -> Path | None:
    """The GTDB table of a table-source metadata record: --metadata-path, the
    table the record names among its inputs, or the release table in the
    workdir. Only gzipped TSVs are read (a tarball would be unpacked beside
    itself, which census must not do)."""
    from ..sources.gtdb_table import workdir_tables

    if metadata_path:
        path = Path(metadata_path)
        if not path.is_file():
            raise UserInputError(f"--metadata-path not found: {metadata_path}")
        return path
    record = cfg.stages["metadata"]
    for name in record.inputs or {}:
        path = Path(name) if Path(name).is_absolute() else wd / name
        if path.name.endswith(".tsv.gz") and path.is_file():
            return path
    rec = record.params
    if rec.get("release") and rec.get("version"):
        for path in workdir_tables(wd, str(rec["release"]), str(rec["version"])):
            if path.name.endswith(".tsv.gz") and path.is_file():
                return path
    return None


def gtdb_candidates(
    wd: Path, cfg: Any, metadata_path: str | None, logger: logging.Logger
) -> tuple[list[Item], str] | None:
    """The genomes the metadata stage chose from: the GTDB table under the
    recorded target, or the GTDB API answer it kept; None when neither is
    in the workdir."""
    from ..sources.gtdb_api import GTDB_API_GENOMES, read_api_genomes

    rec = cfg.stages["metadata"].params
    if rec.get("source") == "api" and not metadata_path:
        path = wd / GTDB_API_GENOMES
        if not path.is_file():
            return None
        return _items_from_parsed(read_api_genomes(path)), "gtdb-api"
    table = _find_table(wd, cfg, metadata_path)
    if table is None:
        return None
    items = gtdb_table_items(
        table,
        logger,
        family=rec.get("target_family"),
        genus=rec.get("target_genus"),
        species=rec.get("target_species") if rec.get("level") == "species" else None,
    )
    return items, "gtdb-table"


def viral_candidates(wd: Path, cfg: Any) -> tuple[str, list[Any], dict[str, int]]:
    """(source, records or items, BV-BRC median lengths) of a vmetadata workdir."""
    from ..viral.ncbi_virus import read_records

    download_wd = wd / "virus_download_wd"
    records_json = download_wd / "virus_records.json"
    if records_json.is_file():
        return "ncbi_virus", read_records(records_json), {}
    if (download_wd / BVBRC_TAXNAMES).is_file():
        items, medians = bvbrc_items(download_wd)
        return "bvbrc", items, medians
    raise UserInputError(
        f"{download_wd} holds neither virus_records.json nor {BVBRC_TAXNAMES}; rerun vmetadata."
    )


def viral_candidates_census(wd: Path, cfg: Any, params: CensusParams) -> Census:
    source, records, medians = viral_candidates(wd, cfg)
    target = cfg.stages["vmetadata"].params.get("target")
    rank, name, recs = viral_scope(records, params.target_family, params.target_genus, target)
    if source == "bvbrc":
        return bvbrc_census(recs, medians, taxon=name, rank=rank, mode=MODE_CANDIDATES)
    return virus_census(recs, taxon=name, rank=rank, mode=MODE_CANDIDATES, source=source)


def _item(row: SelectionRow) -> Item:
    return Item(row.accession, row.family, row.genus, row.species, row.gtdb_representative)


def _ingroup(wd: Path) -> tuple[list[SelectionRow], list[str]]:
    rows = read_selection(wd / SELECTION_TSV)
    return [r for r in rows if not r.is_outgroup], [r.accession for r in rows if r.is_outgroup]


def manifest_sources(wd: Path) -> dict[str, str]:
    """accession -> manifest source, read from a copy of the manifest.

    SQLite creates a -shm index beside a WAL database even for a read-only
    connection; reading a temporary copy (with its -wal) leaves the workdir
    untouched.
    """
    from ..core.manifest import MANIFEST_FILENAME, Manifest

    path = wd / MANIFEST_FILENAME
    if not path.is_file():
        return {}
    with tempfile.TemporaryDirectory(prefix="repgenr-census-") as tmp:
        copy = Path(tmp) / MANIFEST_FILENAME
        shutil.copyfile(path, copy)
        wal = path.with_name(path.name + "-wal")
        if wal.is_file():
            shutil.copyfile(wal, copy.with_name(copy.name + "-wal"))
        manifest = Manifest.open_readonly(copy)
        try:
            genomes = manifest.all_genomes(include_outgroup=True)
            return {g.accession: g.source or "" for g in genomes}
        finally:
            manifest.close()


def selection_census(wd: Path, cfg: Any, params: CensusParams, logger: logging.Logger) -> Census:
    """Mode 3: the selected genomes per genus or species, by manifest source,
    with representatives, clusters and the candidates beside them."""
    rows, outgroups = _ingroup(wd)
    items = [_item(r) for r in rows]
    if params.target_genus or params.target_family:
        rank, name, items = item_scope(items, params.target_family, params.target_genus, None)
    else:
        genera = {_norm(i.genus) for i in items}
        families = {i.family for i in items if i.family}
        rank = GENUS if len(genera) == 1 else FAMILY
        if rank == GENUS:
            name = items[0].genus if items else ""
        elif len(families) == 1:
            name = next(iter(families))
        else:
            name = _recorded_target(cfg) or wd.name
    keep = {i.accession for i in items}
    rows = [r for r in rows if r.accession in keep]
    notes: list[str] = []

    genomes_dir = wd / "genomes"
    if genomes_dir.is_dir():
        present_names = {p.name for p in genomes_dir.iterdir()}
        present = [r for r in rows if r.filename in present_names]
        absent = len(rows) - len(present)
        if absent:
            notes.append(f"{absent} selected genome(s) are absent from genomes/ and not counted.")
    else:
        present = rows
        notes.append("genomes/ does not exist; every selection.tsv row is counted.")

    sources = manifest_sources(wd)
    viral = any(sources.get(r.accession) in VIRAL_MANIFEST_SOURCES for r in rows)
    table = _Rows(rank, binomial=not viral)
    source_totals: dict[str, int] = {}
    for r in present:
        src = sources.get(r.accession) or "unknown"
        source_totals[src] = source_totals.get(src, 0) + 1
        counts = {"genomes": 1, "representatives": int(r.gtdb_representative), src: 1}
        table.add(r.genus, r.species, counts)
    source_columns = [s for s in MANIFEST_SOURCES if s in source_totals]
    source_columns += sorted(s for s in source_totals if s not in MANIFEST_SOURCES)

    columns = ["name", *(["species"] if rank == FAMILY else [])]
    totals: dict[str, Any] = {}
    try:
        candidates = selection_candidates(wd, cfg, params, rank, name, rows, logger)
    except (UserInputError, WorkdirError) as exc:
        # Unreadable or outdated entry-stage data: count the selection alone.
        candidates = None
        notes.append(f"No candidates column: the entry-stage data could not be read ({exc})")
    if candidates is not None:
        cand_items, cand_label = candidates
        for c in cand_items:
            table.add(
                c.genus,
                c.species,
                {"candidates": 1},
                count_species=False,
                label=_display(rank, c.genus_name, c.species_name),
            )
        columns.append("candidates")
        totals["candidates"] = len(cand_items)
        notes.append(f"Candidates are {_CANDIDATE_LABELS.get(cand_label, cand_label)}.")
    columns += ["genomes", *source_columns, "representatives"]

    clusters_tsv = wd / "derep" / CLUSTERS_TSV
    if clusters_tsv.is_file():
        by_name: dict[str, SelectionRow] = {}
        for r in rows:
            by_name[r.filename] = r
            by_name[record_name(r.filename)] = r
        n_clusters = 0
        for rep in read_clusters(clusters_tsv):
            row = by_name.get(rep) or by_name.get(record_name(rep))
            if row is not None:
                genus, species = row.genus, row.species
            elif params.target_genus or params.target_family:
                continue  # a representative outside the narrowed selection
            else:
                _f, genus, species, _a = parse_genome_filename(rep)
            table.add(genus, species, {"clusters": 1}, count_species=False)
            n_clusters += 1
        columns.append("clusters")
        totals["clusters"] = n_clusters

    n_genera, n_species = _distinct_taxa((r.genus, r.species) for r in present)
    totals = {
        "genera": n_genera,
        "species": n_species,
        "genomes": len(present),
        "representatives": sum(1 for r in present if r.gtdb_representative),
        **totals,
        "sources": {s: source_totals[s] for s in source_columns},
        "outgroup": ", ".join(outgroups) or None,
    }
    return Census(
        taxon=name,
        mode=MODE_SELECTION,
        source="selection.tsv",
        rank=rank,
        columns=columns,
        rows=table.finish(columns, dict.fromkeys(columns[1:], 0), "genomes"),
        totals=totals,
        notes=notes,
    )


def _recorded_target(cfg: Any) -> str | None:
    for stage, keys in (
        ("metadata", ("target_genus", "target_family")),
        ("vmetadata", ("target",)),
    ):
        record = cfg.stages.get(stage)
        if record is None:
            continue
        for key in keys:
            if record.params.get(key):
                return str(record.params[key])
    return None


def selection_candidates(
    wd: Path,
    cfg: Any,
    params: CensusParams,
    rank: str,
    name: str,
    rows: list[SelectionRow],
    logger: logging.Logger,
) -> tuple[list[Item], str] | None:
    """The candidates of the selected genera or species, or None without them.

    Candidates are what the entry stage chose from: the GTDB table or API
    answer (metadata), the NCBI Virus or BV-BRC records (vmetadata), or the
    sequencing runs of reads.tsv (reads). They are narrowed to the genus of a
    genus census and to the families of the selection.
    """
    stages = cfg.stages
    found: tuple[list[Any], str] | None = None
    if "metadata" in stages:
        found = gtdb_candidates(wd, cfg, params.metadata_path, logger)
    elif "vmetadata" in stages:
        source, records, _medians = viral_candidates(wd, cfg)
        if source == "ncbi_virus":
            records = [
                Item(r.accession, r.family, r.genus, r.species, species_name=r.species_name)
                for r in records
            ]
        found = (records, source)
    elif (wd / READS_TSV).is_file():
        runs = read_reads(wd / READS_TSV)
        found = (
            [Item(r.run_accession, r.family, r.genus, r.species) for r in runs],
            "reads",
        )
    if found is None:
        return None
    items, label = found
    if rank == GENUS:
        genera = {_norm(r.genus) for r in rows} or {_norm(name)}
        items = [i for i in items if _norm(i.genus) in genera]
    else:
        families = {_norm(r.family) for r in rows if r.family}
        if params.target_family:
            families = {_norm(params.target_family)}
        if families:
            items = [i for i in items if _norm(i.family) in families]
    return items, label


# --- output -------------------------------------------------------------------


def _cell(value: Any) -> str:
    return "" if value is None else str(value)


def write_tsv(path: Path, census: Census) -> None:
    """The rows of a census as a TSV (columns as on the console; no totals line)."""
    try:
        with atomic_path(path) as tmp, open(tmp, "w", encoding="utf-8", newline="") as fo:
            fo.write("\t".join(census.columns) + "\n")
            for row in census.rows:
                fo.write("\t".join(_cell(row.get(c)) for c in census.columns) + "\n")
    except OSError as exc:
        raise UserInputError(f"--tsv {path} could not be written ({exc}).") from exc


def render(census: Census) -> None:
    """Print the totals line, the table and the notes on stdout."""
    from rich.console import Console
    from rich.table import Table

    cells = [[_cell(row.get(c)) for c in census.columns] for row in census.rows]
    widths = [
        max([len(column), *(len(r[i]) for r in cells)]) for i, column in enumerate(census.columns)
    ]
    # Wide enough for every column at its full width: a narrow terminal (or
    # a pipe, which Rich treats as 80 columns) would otherwise cut headers.
    needed = sum(widths) + 2 * len(widths)
    console = Console(highlight=False)
    console = Console(highlight=False, width=max(console.width, needed))
    console.print(census.header(), markup=False, soft_wrap=True)
    table = Table(show_edge=False, pad_edge=False, box=None, header_style="bold")
    for i, column in enumerate(census.columns):
        table.add_column(column, justify="left" if i == 0 else "right", no_wrap=True)
    for row_cells in cells:
        table.add_row(*row_cells)
    console.print(table)
    for note in census.notes:
        console.print(f"Note: {note}", markup=False, soft_wrap=True)


def run_census(params: CensusParams, logger: logging.Logger) -> Census:
    source = validate(params)
    if params.workdir is not None:
        return workdir_census(params, logger)
    return taxon_census(params, source, logger)
