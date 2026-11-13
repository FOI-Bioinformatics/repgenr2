"""NCBI Datasets reports: genome summaries and NCBI Taxonomy, as JSON lines.

Thin helpers over ``datasets summary``. The command runs through the shared
tool runner (containers, retries, timeout) and writes its report to a file,
which is then read line by line, so a taxon with many thousands of genomes
is never held in memory as one document.
"""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .binaries import BinarySpec
from .containers import run_tool_with_retries
from .errors import ToolExecutionError, UserInputError
from .plugins import ToolCapabilities

DATASETS_CAPS = ToolCapabilities(
    name="datasets",
    required_binaries=(BinarySpec("datasets", version_args=("--version",)),),
    conda=("conda-forge::ncbi-datasets-cli",),
)

# A report of a large taxon is tens of thousands of lines; a hung transfer
# must still end.
SUMMARY_TIMEOUT = 3600.0
PROKARYOTE_DOMAINS = frozenset({"Bacteria", "Archaea"})

_UNKNOWN_NAME = "is not recognized"
_AMBIGUOUS_NAME = "more than one taxid"
_TAXID_IN_MESSAGE = re.compile(r"taxid:\s*(\d+)")


def _permanent(exc: ToolExecutionError) -> bool:
    """A name NCBI Taxonomy does not know, or knows twice: a retry cannot help."""
    text = (exc.output or "").lower()
    return _UNKNOWN_NAME in text or _AMBIGUOUS_NAME in text


def summary(args: Sequence[str], report: Path, logger: logging.Logger) -> None:
    """Run ``datasets summary <args> --as-json-lines`` into ``report``."""
    report.parent.mkdir(parents=True, exist_ok=True)
    run_tool_with_retries(
        DATASETS_CAPS,
        ["datasets", "summary", *args, "--as-json-lines"],
        logger=logger,
        log_prefix="datasets",
        stdout_path=report,
        timeout=SUMMARY_TIMEOUT,
        permanent=_permanent,
    )


def read_json_lines(path: Path) -> Iterator[dict]:
    """The JSON objects of a report, one per line; other lines are skipped."""
    with open(path, encoding="utf-8", errors="replace") as fo:
        for raw in fo:
            raw = raw.strip()
            if not raw.startswith("{"):
                continue
            try:
                row = json.loads(raw)
            except json.JSONDecodeError:
                continue
            if isinstance(row, dict):
                yield row


@dataclass(frozen=True)
class Taxon:
    """One NCBI Taxonomy node with its current classification."""

    tax_id: str
    rank: str  # lower case: family, genus, species, subspecies, strain, ...
    name: str
    # rank -> name and rank -> taxid of the classification (domain ... species).
    ranks: dict[str, str] = field(default_factory=dict)
    rank_ids: dict[str, str] = field(default_factory=dict)
    # The taxids above this node, and its own.
    lineage: frozenset[str] = frozenset()

    @property
    def domain(self) -> str:
        # 'superkingdom' in reports written before NCBI Taxonomy named the rank domain.
        return self.ranks.get("domain") or self.ranks.get("superkingdom") or ""


def _taxon(row: dict) -> Taxon | None:
    taxonomy = row.get("taxonomy") or {}
    if taxonomy.get("tax_id") is None:
        return None
    tax_id = str(taxonomy["tax_id"])
    classification = taxonomy.get("classification") or {}
    ranks = {
        rank: str(node.get("name") or "")
        for rank, node in classification.items()
        if isinstance(node, dict)
    }
    rank_ids = {
        rank: str(node["id"])
        for rank, node in classification.items()
        if isinstance(node, dict) and node.get("id") is not None
    }
    name = (taxonomy.get("current_scientific_name") or {}).get("name") or ""
    return Taxon(
        tax_id=tax_id,
        rank=str(taxonomy.get("rank") or "").lower().replace("_", " "),
        name=str(name),
        ranks=ranks,
        rank_ids=rank_ids,
        lineage=frozenset([tax_id, *(str(p) for p in taxonomy.get("parents") or [])]),
    )


def taxonomy(queries: Iterable[str], scratch: Path, logger: logging.Logger) -> dict[str, Taxon]:
    """The NCBI Taxonomy node of each query (a taxid or an exact name).

    One ``datasets summary taxonomy taxon --inputfile`` call; the result is
    keyed by every query as given and by the taxid. A query NCBI does not
    answer is absent from the result. Names can be ambiguous, which fails the
    whole call: resolve a name with :func:`resolve_name` and pass taxids here.
    """
    wanted = sorted({str(q).strip() for q in queries if str(q).strip()})
    if not wanted:
        return {}
    scratch.mkdir(parents=True, exist_ok=True)
    id_file = scratch / "taxonomy_queries.txt"
    report = scratch / "taxonomy_report.jsonl"
    id_file.write_text("\n".join(wanted) + "\n", encoding="utf-8")
    try:
        summary(["taxonomy", "taxon", "--inputfile", str(id_file)], report, logger)
        out: dict[str, Taxon] = {}
        for row in read_json_lines(report):
            taxon = _taxon(row)
            if taxon is None:
                continue
            for key in [*(str(q) for q in row.get("query") or []), taxon.tax_id]:
                out.setdefault(key, taxon)
        return out
    finally:
        id_file.unlink(missing_ok=True)
        report.unlink(missing_ok=True)


def resolve_name(name: str, rank: str, scratch: Path, logger: logging.Logger) -> Taxon:
    """The bacterial or archaeal taxon called ``name`` at ``rank``.

    A name that NCBI Taxonomy uses for more than one taxon (Bacillus is also
    a genus of stick insects) is settled by the domain and the rank; a name
    that stays ambiguous, is unknown, or names another rank or domain raises
    :class:`UserInputError`.
    """
    try:
        found = taxonomy([name], scratch, logger).values()
    except ToolExecutionError as exc:
        text = exc.output or ""
        if _UNKNOWN_NAME in text.lower():
            raise UserInputError(
                f"NCBI Taxonomy has no taxon named '{name}'. Check the spelling "
                "(NCBI names, not GTDB names such as Bacillus_A)."
            ) from exc
        if _AMBIGUOUS_NAME not in text.lower():
            raise
        ids = _TAXID_IN_MESSAGE.findall(text)
        found = taxonomy(ids, scratch, logger).values()
    taxa = {t.tax_id: t for t in found}
    fitting = [t for t in taxa.values() if t.domain in PROKARYOTE_DOMAINS and t.rank == rank]
    if len(fitting) == 1:
        return fitting[0]
    if not taxa:
        raise UserInputError(f"NCBI Taxonomy has no taxon named '{name}'.")
    described = "; ".join(
        f"{t.name} (taxid {t.tax_id}, {t.rank or 'no rank'}, {t.domain or 'no domain'})"
        for t in sorted(taxa.values(), key=lambda t: t.tax_id)
    )
    if fitting:
        raise UserInputError(
            f"'{name}' names more than one bacterial or archaeal {rank} in NCBI "
            f"Taxonomy: {described}."
        )
    raise UserInputError(
        f"'{name}' is not a bacterial or archaeal {rank} in NCBI Taxonomy: {described}. "
        "Check the name and --level."
    )
