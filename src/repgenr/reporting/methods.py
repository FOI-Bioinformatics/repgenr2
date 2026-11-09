"""Methods text, software versions and parameters from the stage records.

Sentences are filled only from values a stage recorded; a value that was
not recorded is left out rather than guessed. A stage without a template
gets one generic sentence naming its tool and version.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any

from ..core.versions import merge_stage_versions
from .model import MethodsInfo

# Tool -> (citation, DOI). Only references checked against the publisher.
CITATIONS: dict[str, tuple[str, str]] = {
    "gtdb": (
        "Parks DH et al. (2022) GTDB: an ongoing census of bacterial and archaeal "
        "diversity through a phylogenetically consistent, rank normalized and complete "
        "genome-based taxonomy. Nucleic Acids Res 50:D785-D794.",
        "10.1093/nar/gkab776",
    ),
    "datasets": (
        "O'Leary NA et al. (2024) Exploring and retrieving sequence and metadata for "
        "species across the tree of life with NCBI Datasets. Sci Data 11:732.",
        "10.1038/s41597-024-03571-y",
    ),
}

# Parameters left out of the appendix: bookkeeping, not settings.
_HIDDEN_PARAMS = frozenset(
    {"source_states", "selected_count", "downloaded", "total", "keeper_swaps"}
)


# Tool names as their authors write them.
DISPLAY_NAMES: dict[str, str] = {
    "drep": "dRep",
    "skder": "skDER",
    "galah": "galah",
    "sourmash": "sourmash",
    "iqtree": "IQ-TREE",
    "raxmlng": "RAxML-NG",
    "fasttree": "FastTree",
    "veryfasttree": "VeryFastTree",
    "mashtree": "Mashtree",
    "ska2": "SKA2",
    "skani": "skani",
    "gubbins": "Gubbins",
    "panaroo": "Panaroo",
    "ppanggolin": "PPanGGOLiN",
    "prokka": "Prokka",
    "bakta": "Bakta",
}


def display(tool: str | None) -> str:
    return DISPLAY_NAMES.get(tool or "", tool or "")


def _version(m: MethodsInfo, tool: str | None) -> str:
    """' 3.7.1', or ' (container image quay.io/...)' for a tool run from an image."""
    if not tool:
        return ""
    v = m.tool_versions.get(tool)
    if not v:
        return ""
    if "/" in v or "@sha256:" in v:
        return f" (container image {v})"
    return f" {v}"


def _metadata(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    target = p.get("target_species") or p.get("target_genus") or p.get("target_family")
    level = p.get("level")
    if p.get("source") == "api":
        where = "the GTDB web API"
        if p.get("api_query_date"):
            where += f" (queried on {str(p['api_query_date'])[:10]})"
    else:
        where = "the GTDB metadata table"
        if p.get("release"):
            where += f" of release {p['release']}"
    text = "Genomes were selected from " + where
    if target:
        text += f" for {target}" + (f" ({level})" if level else "")
    text += "."
    if p.get("dataset") == "rep":
        text += " Only GTDB species representatives were considered."
    if p.get("outgroup"):
        text += f" Genome {p['outgroup']} was taken as the outgroup."
    return text, ["gtdb"]


def _genome(m: MethodsInfo) -> tuple[str, list[str]]:
    return (
        f"Genome assemblies were downloaded from NCBI with datasets{_version(m, 'datasets')}.",
        ["datasets"],
    )


def _ingest(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    parts = []
    if p.get("genomes_dir"):
        parts.append("a local folder of genome assemblies")
    if p.get("from_workdirs"):
        n = len(p["from_workdirs"])
        parts.append(f"{n} earlier working director{'y' if n == 1 else 'ies'}")
    text = "Genomes were taken from " + (" and ".join(parts) or "local files") + "."
    if p.get("outgroup"):
        text += f" Genome {p['outgroup']} was taken as the outgroup."
    return text, []


def _ani(value: Any) -> str | None:
    return f"{100 * float(value):g}%" if isinstance(value, int | float) else None


def _dereplicate(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    tool = display(m.tool) or "the chosen tool"
    text = f"Genomes were dereplicated with {tool}{_version(m, m.tool)}"
    parts = []
    if _ani(p.get("primary_ani")):
        parts.append(f"a primary ANI threshold of {_ani(p['primary_ani'])}")
    if _ani(p.get("secondary_ani")):
        parts.append(f"a secondary ANI threshold of {_ani(p['secondary_ani'])}")
    if parts:
        text += " with " + " and ".join(parts)
    text += "."
    if p.get("target_reps"):
        text += (
            f" The secondary threshold was searched to approach {p['target_reps']} representatives."
        )
    keeper = p.get("keeper_effective") or p.get("keeper")
    if keeper and keeper != "tool":
        text += f" The representative of each cluster was chosen by the '{keeper}' keeper rule."
    if p.get("reduce") and p["reduce"] != "none":
        text += f" One representative was then kept per {p['reduce']}."
    return text, []


def _derep_refine(m: MethodsInfo) -> tuple[str, list[str]]:
    text = "Selected lineages were dereplicated again at a finer threshold"
    if _ani(m.params.get("secondary_ani")):
        text += f" ({_ani(m.params['secondary_ani'])} ANI)"
    return text + ".", []


def _glance(m: MethodsInfo) -> tuple[str, list[str]]:
    tool = display(m.tool) or "the chosen tool"
    return (
        f"Pairwise similarities between all genomes were computed with {tool}"
        f"{_version(m, m.tool)}.",
        [],
    )


_TEMPLATES: dict[str, Callable[[MethodsInfo], tuple[str, list[str]]]] = {
    "metadata": _metadata,
    "genome": _genome,
    "ingest": _ingest,
    "glance": _glance,
    "dereplicate": _dereplicate,
    "derep_refine": _derep_refine,
}


def sentence(m: MethodsInfo) -> tuple[str, list[str]]:
    """The methods sentence of one stage record and the citation keys it uses."""
    template = _TEMPLATES.get(m.stage)
    if template is not None:
        return template(m)
    if m.tool:
        return f"The {m.stage} stage used {display(m.tool)}{_version(m, m.tool)}.", []
    return f"The {m.stage} stage was run.", []


def unique_records(records: Iterable[MethodsInfo]) -> list[MethodsInfo]:
    """One record per stage, in first-seen order (several steps may carry the
    same stage)."""
    seen: dict[str, MethodsInfo] = {}
    for m in records:
        seen.setdefault(m.stage, m)
    return list(seen.values())


def software_versions(records: Iterable[MethodsInfo]) -> list[tuple[str, str]]:
    merged = merge_stage_versions({m.stage: m.tool_versions for m in records})
    return sorted(merged.items())


def parameter_rows(records: Iterable[MethodsInfo]) -> list[tuple[str, str, str]]:
    """(stage, parameter, value) for the recorded parameters, nested values
    flattened with dots; empty values left out."""
    rows: list[tuple[str, str, str]] = []
    for m in records:
        for key, value in _flatten(m.params):
            if key.split(".")[0] in _HIDDEN_PARAMS or value in ("", "None", "[]", "{}"):
                continue
            rows.append((m.stage, key, value))
    return rows


def _flatten(data: dict[str, Any], prefix: str = "") -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for key in sorted(data):
        value = data[key]
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.extend(_flatten(value, name + "."))
        elif isinstance(value, list | tuple):
            out.append((name, ", ".join(str(v) for v in value)))
        else:
            out.append((name, str(value)))
    return out
