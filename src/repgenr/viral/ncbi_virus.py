"""NCBI Virus data source via the ``datasets`` CLI.

A single ``datasets download virus genome taxon`` call returns both the sequences
(``genomic.fna``) and structured per-sequence metadata (``data_report.jsonl``):
accession, ranked taxonomy, length, completeness, host, segment, isolate. This
replaces the BV-BRC FTP download + the separate NCBI Entrez taxonomy step with
one maintained dependency RepGenR2 already uses, and yields stable accessions so
viral genomes can adopt the canonical ``Family_Genus_species_Accession`` naming.
"""

from __future__ import annotations

import json
import logging
import re
import shutil
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

from ..core import process
from ..core.binaries import BinarySpec
from ..core.containers import run_tool, run_tool_with_retries
from ..core.contracts import atomic_path
from ..core.errors import MissingBinaryError, ToolExecutionError, WorkdirError
from ..core.http import NCBI_DATASETS_URL, require_reachable
from ..core.plugins import ToolCapabilities
from ..core.process import remove_tree

DATASETS_CAPS = ToolCapabilities(
    name="datasets",
    required_binaries=(BinarySpec("datasets", version_args=("--version",)),),
    conda=("conda-forge::ncbi-datasets-cli",),
)

_ANONYMOUS = "ANONYMOUS"  # NCBI's value for non-segmented virus sequences


@dataclass
class VirusRecord:
    accession: str
    taxid: str
    organism: str
    family: str
    genus: str
    species: str
    length: int
    completeness: str
    segment: str  # the label as submitted; see normalise_segment
    isolate: str
    # Virus lineage names, root to leaf. Absent from records written before
    # the species came from the lineage (see read_records).
    lineage: list[str] = field(default_factory=list)
    # Where the species came from: 'taxonomy' (NCBI Taxonomy, current
    # classification of the taxid), 'lineage' (ICTV binomial in the report
    # lineage) or 'organism' (no binomial; the organism name).
    species_source: str = ""
    # The species as NCBI writes it ('Hepatovirus ahepa'), for display; the
    # species field above is its filename token. Empty in older records.
    species_name: str = ""


def _sanitize(name: str) -> str:
    """Make a taxonomy name a single safe token (no spaces/underscores), so the
    canonical filename round-trips through ``parse_genome_filename``."""
    token = name.strip().replace(" ", "-").replace("_", "-")
    token = re.sub(r"[^A-Za-z0-9.-]", "", token)
    return token or "NA"


def _genus_like(lineage_names: list[str], leaf: str) -> list[str]:
    """Single-word lineage names ending in '-virus', root to leaf: the genus
    and, where one exists, the subgenus (Betacoronavirus > Embecovirus)."""
    return [n for n in lineage_names if " " not in n and n.lower().endswith("virus") and n != leaf]


def _classify(lineage_names: list[str], organism: str) -> tuple[str, str, str | None]:
    """Derive (family, genus, binomial) from an NCBI Virus lineage.

    The report's lineage entries carry names and taxids but no rank labels, so
    ICTV name forms are used: the family ends in ``-viridae``, and genera and
    subgenera are single-word names ending in ``-virus``. The binomial (see
    :func:`species_from_lineage`) is ``None`` when the lineage holds none; its
    first word is then the genus, since the shallowest single-word name is
    the genus only when no subgenus follows it. Without a binomial the genus
    is the shallowest single-word '-virus' name, the organism when it is such
    a name itself, else the name above the leaf unless that is the family.
    """
    leaf = organism or (lineage_names[-1] if lineage_names else "")
    family = next((n for n in lineage_names if n.lower().endswith("viridae")), "")
    candidates = _genus_like(lineage_names, leaf)
    binomial = species_from_lineage(lineage_names, candidates)
    if binomial is not None:
        genus = binomial.split(" ", 1)[0]
    elif candidates:
        genus = candidates[0]
    elif " " not in leaf and leaf.lower().endswith("virus"):
        # The organism is the genus itself (a record filed at genus level).
        genus = leaf
    elif (
        len(lineage_names) >= 2
        and " " not in lineage_names[-2]
        and not lineage_names[-2].lower().endswith("viridae")
    ):
        genus = lineage_names[-2]
    else:
        genus = ""  # never the family
    return family, genus, binomial


def species_from_lineage(lineage_names: list[str], genera: list[str]) -> str | None:
    """The ICTV binomial in a lineage: the shallowest name made of one of
    ``genera`` and a single lower-case epithet, which may hold digits and
    hyphens ('Mammarenavirus brazilense', 'Mammarenavirus dhati-welelense',
    'Lentivirus humimdef1'); ``None`` when there is none.

    NCBI keeps strain-level and earlier names below the binomial, some of
    which also start with the genus: 'Hepatovirus ahepa' > 'Hepatovirus A',
    'Betacoronavirus gravedinis' > 'Betacoronavirus 1', and
    'Orthobunyavirus cacheense' > 'Orthobunyavirus maguariense'. The first
    two are not of the binomial form, and the shallowest match is the
    species, so names below it do not split one species into several.
    """
    for name in lineage_names:
        if any(re.fullmatch(re.escape(g) + r" [a-z][a-z0-9-]*", name) for g in genera):
            return name
    return None


# Segment labels are free text in NCBI Virus: 'M', 'M; medium', 'middle',
# 'S RNA', 'RNA 2', 'DNA-A', 'segment 4' and so on.
_SEGMENT_SIZES = {"SMALL": "S", "MEDIUM": "M", "MIDDLE": "M", "LARGE": "L"}
# Words that name the molecule or the word 'segment' rather than which segment.
_SEGMENT_GENERIC = frozenset(
    {"SEGMENT", "SEG", "GENOME", "COMPONENT", "CIRCULAR", "RNA", "DSRNA", "SSRNA", "DNA"}
)
_MOLECULE_PREFIX = re.compile(r"^(?:DS|SS)?(?:RNA|DNA)(?=\d)")
# Normalised labels that do not identify a segment.
UNLABELLED_SEGMENTS = frozenset({"", _ANONYMOUS, "UNKNOWN"})


def normalise_segment(label: str | None) -> str:
    """The segment a free-text label names, so that variants compare equal.

    The text after a ';' is a comment and is dropped. The rest is split on
    whitespace, '-' and '_'; words that only name the molecule or the word
    segment ('RNA', 'DNA', 'segment', 'genome', 'component') are skipped when
    another word remains ('circular' too), and the first remaining word is
    kept in upper case.
    A molecule prefix before a number is removed ('RNA1' is '1'), and small,
    medium, middle and large become S, M, M and L. So 'M', 'M; medium' and
    'middle' are all 'M', 'S RNA' is 'S', 'DNA-A' is 'A', while 'RNA 1' and
    'RNA 2' stay distinct. Empty, 'ANONYMOUS' and 'Unknown' fall in
    :data:`UNLABELLED_SEGMENTS`.
    """
    head = (label or "").split(";", 1)[0]
    words = [w.upper() for w in re.split(r"[\s_-]+", head.strip()) if w]
    if not words:
        return ""
    specific = [w for w in words if w not in _SEGMENT_GENERIC]
    if not specific:
        return words[0]
    word = _MOLECULE_PREFIX.sub("", specific[0])
    return _SEGMENT_SIZES.get(word, word)


@dataclass(frozen=True)
class TaxonClass:
    """The current NCBI Taxonomy classification of one taxid."""

    family: str
    genus: str
    species: str


TAXONOMY_TIMEOUT = 300.0


def lookup_taxonomy(
    taxids: Iterable[str],
    out_dir: Path,
    *,
    logger: logging.Logger,
    runner: Callable[..., int] | None = None,
) -> dict[str, TaxonClass]:
    """Current family, genus and species of each taxid from NCBI Taxonomy.

    One ``datasets summary taxonomy taxon --inputfile`` call covers every
    distinct taxid. The report lineage of NCBI Virus nests species that NCBI
    Taxonomy keeps as siblings (Maguari virus under 'Orthobunyavirus
    cacheense' > 'Orthobunyavirus maguariense', while its current species is
    'Orthobunyavirus maguariense'), so the classification is preferred over
    the lineage. The lookup runs once, without retries, right after the
    package download, whose reachability check (``require_reachable``) covers
    the same host; a failure is logged and returns an empty mapping, and the
    caller then uses the lineage rule. ``runner`` is injectable for tests.
    """
    ids = sorted({t for t in taxids if t})
    if not ids:
        return {}
    if runner is None:
        runner = run_tool
    out_dir.mkdir(parents=True, exist_ok=True)
    id_file = out_dir / "taxids.txt"
    report = out_dir / "taxonomy_report.jsonl"
    id_file.write_text("\n".join(ids) + "\n", encoding="utf-8")
    cmd = [
        "datasets",
        "summary",
        "taxonomy",
        "taxon",
        "--inputfile",
        str(id_file),
        "--as-json-lines",
    ]
    try:
        runner(
            DATASETS_CAPS,
            cmd,
            logger=logger,
            log_prefix="datasets",
            stdout_path=report,
            timeout=TAXONOMY_TIMEOUT,
        )
        lines = report.read_text(encoding="utf-8").splitlines()
    except (ToolExecutionError, MissingBinaryError, OSError) as exc:
        logger.warning(
            "NCBI Taxonomy lookup of %d taxid(s) failed (%s); the species comes "
            "from the report lineage",
            len(ids),
            exc,
        )
        return {}
    finally:
        id_file.unlink(missing_ok=True)
    out: dict[str, TaxonClass] = {}
    for raw in lines:
        raw = raw.strip()
        if not raw:
            continue
        try:
            row = json.loads(raw)
        except json.JSONDecodeError:
            continue
        taxonomy = row.get("taxonomy") or {}
        cls = taxonomy.get("classification") or {}

        def _name(rank: str, cls: dict = cls) -> str:
            return str((cls.get(rank) or {}).get("name") or "")

        keys = [str(q) for q in (row.get("query") or [])]
        if taxonomy.get("tax_id") is not None:
            keys.append(str(taxonomy["tax_id"]))
        entry = TaxonClass(_name("family"), _name("genus"), _name("species"))
        for key in keys:
            out.setdefault(key, entry)
    report.unlink(missing_ok=True)
    unresolved = [t for t in ids if not (out.get(t) and out[t].species)]
    if unresolved:
        logger.info(
            "NCBI Taxonomy gave no species for %d of %d taxid(s) (%s); their species "
            "comes from the report lineage",
            len(unresolved),
            len(ids),
            ", ".join(unresolved[:5]) + (" ..." if len(unresolved) > 5 else ""),
        )
    return out


def report_taxids(lines: Iterable[str]) -> set[str]:
    """The distinct taxids of a report stream (both field casings)."""
    taxids: set[str] = set()
    for raw in lines:
        raw = raw.strip()
        if not raw:
            continue
        virus = json.loads(raw).get("virus", {}) or {}
        taxid = virus.get("taxId") or virus.get("tax_id")
        if taxid:
            taxids.add(str(taxid))
    return taxids


def parse_report(
    lines: Iterable[str],
    logger: logging.Logger | None = None,
    taxonomy: dict[str, TaxonClass] | None = None,
) -> list[VirusRecord]:
    """Parse a ``data_report.jsonl`` (or ``summary`` JSONL) stream into records.

    Handles both field casings: the download report uses ``organismName``/
    ``taxId``; the summary report uses ``organism_name``/``tax_id``.

    The species of a record is, in order of preference: the current species
    of its taxid in NCBI Taxonomy (``taxonomy``, from :func:`lookup_taxonomy`);
    the ICTV binomial from the report lineage; the organism name. For the
    lineage rule one taxid can carry different lineages across records
    (Mudanjiang phlebovirus under 'Phlebovirus baishanense' in one record and
    'unclassified Phlebovirus' in another), so every record of a taxid takes
    the binomial found on most records of that taxid (ties by name). The
    number of records that keep the organism name is logged, as are taxids
    with conflicting lineage binomials.
    """
    taxonomy = taxonomy or {}
    rows: list[tuple[dict, str, str, list[str], str, str, str | None]] = []
    by_taxid: dict[str, dict[str, int]] = {}
    for raw in lines:
        raw = raw.strip()
        if not raw:
            continue
        row = json.loads(raw)
        virus = row.get("virus", {}) or {}
        organism = virus.get("organismName") or virus.get("organism_name") or ""
        taxid = str(virus.get("taxId") or virus.get("tax_id") or "")
        lineage_names = [e.get("name", "") for e in (virus.get("lineage") or []) if e.get("name")]
        family, genus, binomial = _classify(lineage_names, organism)
        if binomial is not None and taxid:
            counts = by_taxid.setdefault(taxid, {})
            counts[binomial] = counts.get(binomial, 0) + 1
        rows.append((row, organism, taxid, lineage_names, family, genus, binomial))

    chosen = {
        taxid: min(counts, key=lambda name: (-counts[name], name))
        for taxid, counts in by_taxid.items()
    }
    conflicts = sorted(
        t
        for t, counts in by_taxid.items()
        if len(counts) > 1 and not (taxonomy.get(t) and taxonomy[t].species)
    )
    records: list[VirusRecord] = []
    fallback: dict[str, int] = {}  # organism name -> records without a binomial
    for row, organism, taxid, lineage_names, family, genus, binomial in rows:
        current = taxonomy.get(taxid)
        if current is not None and current.species:
            species, source = current.species, "taxonomy"
            family = current.family or family
            genus = current.genus or (species.split(" ", 1)[0] if " " in species else genus)
        else:
            if taxid in chosen:
                binomial = chosen[taxid]
                genus = binomial.split(" ", 1)[0]
            if binomial is None:
                species, source = organism, "organism"
                fallback[organism] = fallback.get(organism, 0) + 1
            else:
                species, source = binomial, "lineage"
            if current is not None:
                family = current.family or family
                genus = current.genus or genus
        records.append(
            VirusRecord(
                accession=row.get("accession", "") or "",
                taxid=taxid,
                organism=organism,
                family=_sanitize(family),
                genus=_sanitize(genus),
                species=_sanitize(species),
                length=int(row.get("length", 0) or 0),
                completeness=(row.get("completeness", "") or "").upper(),
                segment=(row.get("segment", "") or _ANONYMOUS),
                isolate=((row.get("isolate") or {}).get("name", "") or ""),
                lineage=lineage_names,
                species_source=source,
                species_name=species,
            )
        )
    if conflicts and logger is not None:
        logger.info(
            "%d taxid(s) carry different binomials across records; each takes the "
            "most frequent: %s",
            len(conflicts),
            "; ".join(f"{t} -> {chosen[t]}" for t in conflicts[:5]),
        )
    if fallback and logger is not None:
        logger.info(
            "%d of %d record(s) (%d organism name(s)) have no species from NCBI "
            "Taxonomy and no ICTV binomial in their lineage; their species is the "
            "organism name",
            sum(fallback.values()),
            len(records),
            len(fallback),
        )
    return records


def fetch(
    target: str,
    out_dir: Path,
    *,
    complete_only: bool = False,
    host: str | None = None,
    released_after: str | None = None,
    logger: logging.Logger,
    runner: Callable[..., int] | None = None,
    resolve_taxonomy: bool = True,
    taxonomy_runner: Callable[..., int] | None = None,
) -> list[VirusRecord]:
    """Download an NCBI Virus package for ``target`` and return its records.

    The species of each record is resolved through NCBI Taxonomy (see
    :func:`lookup_taxonomy`) unless ``resolve_taxonomy`` is false.

    Writes the package sequences to ``out_dir/download.fa`` (headers
    ``>accession ...``). ``runner`` is injectable for tests.
    """
    if runner is None:

        def runner(caps, cmd, **kw):
            # datasets performs its own network transfers with no built-in
            # retry; cap and retry it like the bacterial download path.
            return run_tool_with_retries(caps, cmd, timeout=3600.0, **kw)

    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / "ncbi_virus.zip"
    extract = out_dir / "ncbi_virus_pkg"
    if extract.exists():
        remove_tree(extract)

    cmd = [
        "datasets",
        "download",
        "virus",
        "genome",
        "taxon",
        target,
        "--filename",
        str(zip_path),
        "--no-progressbar",
    ]
    if complete_only:
        cmd.append("--complete-only")
    if host:
        cmd += ["--host", host]
    if released_after:
        cmd += ["--released-after", released_after]
    require_reachable(NCBI_DATASETS_URL, what="NCBI datasets")
    runner(DATASETS_CAPS, cmd, logger=logger, log_prefix="datasets")

    process.unzip(zip_path, extract)
    report = extract / "ncbi_dataset" / "data" / "data_report.jsonl"
    fna = extract / "ncbi_dataset" / "data" / "genomic.fna"
    if not report.exists() or not fna.exists():
        raise WorkdirError(
            f"NCBI Virus returned no genomes for taxon '{target}'. Check the name "
            "or loosen the filters (--complete-only/--host/--released-after)."
        )
    lines = report.read_text().splitlines()
    taxonomy = (
        lookup_taxonomy(report_taxids(lines), out_dir, logger=logger, runner=taxonomy_runner)
        if resolve_taxonomy
        else {}
    )
    records = parse_report(lines, logger, taxonomy)
    with atomic_path(out_dir / "download.fa") as tmp:
        shutil.copyfile(fna, tmp)
    zip_path.unlink(missing_ok=True)
    shutil.rmtree(extract, ignore_errors=True)
    logger.info("NCBI Virus: %d sequences for taxon '%s'", len(records), target)
    return records


SUMMARY_TIMEOUT = 1800.0


def summary(
    target: str,
    *,
    complete_only: bool = False,
    host: str | None = None,
    released_after: str | None = None,
    logger: logging.Logger,
    runner: Callable[..., int] | None = None,
    taxonomy_runner: Callable[..., int] | None = None,
) -> list[VirusRecord]:
    """The NCBI Virus records of ``target`` without their sequences.

    ``datasets summary virus genome taxon`` returns the same per-sequence
    report as the package :func:`fetch` downloads, with the same filters, so
    the records (species resolved through NCBI Taxonomy as in :func:`fetch`)
    are those vmetadata would select. Nothing is written outside a temporary
    directory. ``runner`` and ``taxonomy_runner`` are injectable for tests.
    """
    import tempfile

    if runner is None:
        runner = run_tool
    cmd = ["datasets", "summary", "virus", "genome", "taxon", target, "--as-json-lines"]
    if complete_only:
        cmd.append("--complete-only")
    if host:
        cmd += ["--host", host]
    if released_after:
        cmd += ["--released-after", released_after]
    require_reachable(NCBI_DATASETS_URL, what="NCBI datasets")
    with tempfile.TemporaryDirectory(prefix="repgenr-virus-summary-") as tmp:
        report = Path(tmp) / "summary.jsonl"
        runner(
            DATASETS_CAPS,
            cmd,
            logger=logger,
            log_prefix="datasets",
            stdout_path=report,
            timeout=SUMMARY_TIMEOUT,
        )
        lines = report.read_text(encoding="utf-8").splitlines() if report.exists() else []
        taxonomy = lookup_taxonomy(
            report_taxids(lines), Path(tmp), logger=logger, runner=taxonomy_runner
        )
    records = parse_report(lines, logger, taxonomy)
    logger.info("NCBI Virus: %d sequence records for taxon '%s'", len(records), target)
    return records


def write_records(path: Path, records: list[VirusRecord]) -> None:
    path.write_text(json.dumps([asdict(r) for r in records]), encoding="utf-8")


def read_records(path: Path) -> list[VirusRecord]:
    """Read ``virus_records.json``.

    Records written before the species was resolved through NCBI Taxonomy
    and the lineage hold no ``species_source`` key and may carry the organism
    name as species; selecting from them would give the earlier filenames and
    species, so they are refused with a request to rerun vmetadata.
    """
    rows = json.loads(path.read_text(encoding="utf-8"))
    if rows and any("species_source" not in row for row in rows):
        raise WorkdirError(
            f"{path} was written by an earlier version, which took the viral species "
            "from the organism name. Rerun the download once with "
            "'repgenr --force vmetadata ...' (same arguments) so that the species "
            "comes from NCBI Taxonomy and the lineage; vgenome and the later stages "
            "then rerun."
        )
    known = {f.name for f in fields(VirusRecord)}
    return [VirusRecord(**{k: v for k, v in row.items() if k in known}) for row in rows]
