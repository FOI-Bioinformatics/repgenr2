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
    "spades": "SPAdes",
    "flye": "Flye",
    "medaka": "Medaka",
    "racon": "Racon",
    "checkm2": "CheckM2",
    "dendropy": "DendroPy",
}


# Tool name -> the binary its version is recorded under, where they differ.
BINARY_NAMES: dict[str, str] = {
    "raxmlng": "raxml-ng",
    "gubbins": "run_gubbins.py",
    "ska2": "ska",
    "aster": "wastral",
}


def tool_version(versions: dict[str, str], tool: str) -> str | None:
    """The recorded version of ``tool``: under its name, the name of its
    binary, or either in another letter case (dRep, FastTree)."""
    for key in (tool, BINARY_NAMES.get(tool, "")):
        if key and key in versions:
            return versions[key]
    lowered = {k.lower(): v for k, v in versions.items()}
    for key in (tool, BINARY_NAMES.get(tool, "")):
        if key and key.lower() in lowered:
            return lowered[key.lower()]
    return None


def display(tool: str | None) -> str:
    return DISPLAY_NAMES.get(tool or "", tool or "")


def _version(m: MethodsInfo, tool: str | None) -> str:
    """' 3.7.1', or ' (container image quay.io/...)' for a tool run from an image."""
    if not tool:
        return ""
    v = tool_version(m.tool_versions, tool)
    if not v:
        return ""
    if "/" in v or "@sha256:" in v:
        return f" (container image {v})"
    return f" {v}"


def _metadata(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    target = p.get("target_species") or p.get("target_genus") or p.get("target_family")
    level = p.get("level")
    ncbi = p.get("source") == "ncbi"
    if ncbi:
        held = {"refseq": "RefSeq", "genbank": "GenBank", "all": "RefSeq and GenBank"}
        where = f"the {held.get(str(p.get('assembly_source')), 'NCBI')} assemblies at NCBI"
        where += f", by NCBI taxonomy, with datasets{_version(m, 'datasets')}"
        if p.get("query_date"):
            where += f" (queried on {str(p['query_date'])[:10]})"
    elif p.get("source") == "api":
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
        text += (
            " Only NCBI reference genomes were considered."
            if ncbi
            else " Only GTDB species representatives were considered."
        )
    if p.get("outgroup"):
        text += f" Genome {p['outgroup']} was taken as the outgroup"
        rules = _NCBI_OUTGROUP_RULES if ncbi else _OUTGROUP_RULES
        text += rules.get(str(p.get("outgroup_rule")), "") + "."
    return text, ["datasets"] if ncbi else ["gtdb"]


_NCBI_OUTGROUP_RULES = {
    "quality": (
        ": the NCBI reference genome of best assembly quality (CheckM completeness"
        " - 5 x contamination, isolates before metagenome-assembled genomes) in the parent"
        " taxon of the target"
    ),
}

_OUTGROUP_RULES = {
    "sister": (
        ": the GTDB species representative of best assembly quality (CheckM2 completeness"
        " - 5 x contamination, isolates before metagenome-assembled genomes) in the group"
        " sister to the target in the GTDB reference tree"
    ),
    "quality": (
        ": the GTDB species representative of best assembly quality (CheckM2 completeness"
        " - 5 x contamination, isolates before metagenome-assembled genomes) in the parent"
        " taxon of the target"
    ),
}


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


def _reads(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    target = p.get("target_species") or p.get("target_genus") or p.get("target_family")
    text = "Whole-genome sequencing runs"
    if target:
        text += f" of {target}"
    text += " were selected from the European Nucleotide Archive"
    filters = []
    if p.get("platform") and p["platform"] != "any":
        filters.append(f"{p['platform']} runs")
    if p.get("min_bases"):
        filters.append(f"at least {int(p['min_bases']):,} bases")
    if p.get("max_bases"):
        filters.append(f"at most {int(p['max_bases']):,} bases")
    if p.get("drop_selection"):
        filters.append("library selection not " + ", ".join(p["drop_selection"]))
    if filters:
        text += " (" + "; ".join(filters) + ")"
    text += "."
    if p.get("one_per_sample"):
        text += " One run was kept per biosample."
    return text, []


def _assemble(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    tools = p.get("assemblers_used") or ([m.tool] if m.tool else [])
    names = ", ".join(f"{display(t)}{_version(m, t)}" for t in tools) or "the chosen assembler"
    text = f"Runs were assembled with {names}."
    polishers = p.get("polishers_used") or []
    if polishers:
        text += (
            " Assemblies were polished with "
            + ", ".join(f"{display(t)}{_version(m, t)}" for t in polishers)
            + "."
        )
    if p.get("checkm2_db"):
        text += (
            f" Completeness and contamination were estimated with CheckM2{_version(m, 'checkm2')}."
        )
    if p.get("classifier_effective"):
        c = p["classifier_effective"]
        text += f" Assemblies were classified with {display(c)}{_version(m, c)}."
    return text, []


def _snptype(m: MethodsInfo) -> tuple[str, list[str]]:
    text = (
        f"Core SNPs were called with {display(m.tool) or 'the chosen typer'}{_version(m, m.tool)}"
    )
    if m.params.get("reference"):
        text += f" against {m.params['reference']}"
    text += "."
    if m.params.get("mask"):
        mask = m.params["mask"]
        text += f" Recombination was masked with {display(mask)}{_version(m, mask)}."
    return text, []


def _phylo(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    source = p.get("msa_source")
    tool = display(m.tool) or "the chosen tree builder"
    if source == "coregene":
        builder = p.get("coregene")
        text = (
            f"Core genes were identified with {display(builder)}{_version(m, builder)}"
            if builder
            else "Core genes were identified"
        )
        if isinstance(p.get("core_threshold"), int | float):
            text += f" (present in at least {100 * p['core_threshold']:g}% of the genomes)"
        text += ", aligned and concatenated."
        if p.get("partitioned"):
            text += " Each gene was given its own substitution model."
    elif source == "snptype":
        typer = p.get("snptyper")
        text = f"A core-SNP alignment was built with {display(typer)}{_version(m, typer)}."
    elif source == "aligner":
        aligner = p.get("aligner")
        text = f"A whole-genome alignment was built with {display(aligner)}{_version(m, aligner)}."
    else:
        text = ""
    text += f" The tree was inferred with {tool}{_version(m, m.tool)}"
    if p.get("bootstrap"):
        text += f" with {p['bootstrap']} bootstrap replicates"
    text += "."
    if p.get("root") == "outgroup" and p.get("outgroup"):
        text += f" It was rooted on the outgroup {p['outgroup']}."
    elif p.get("root") == "midpoint":
        text += " It was rooted at the midpoint."
    return text.strip(), []


def _tree2tax(m: MethodsInfo) -> tuple[str, list[str]]:
    p = m.params
    text = "The rooted tree was converted into a node taxonomy for FlexTaxD"
    if m.tool_versions.get("dendropy"):
        text += f" with DendroPy {m.tool_versions['dendropy']}"
    text += "."
    rules = []
    if isinstance(p.get("collapse_support"), int | float):
        rules.append(f"support below {p['collapse_support']:g}")
    if isinstance(p.get("collapse_length"), int | float):
        rules.append(f"branch length below {p['collapse_length']:g}")
    if rules:
        text += " Internal nodes with " + " or ".join(rules) + " were collapsed."
    if isinstance(p.get("collapse_diameter"), int | float):
        text += f" Clades within a leaf-to-leaf distance of {p['collapse_diameter']:g} were merged."
    if p.get("ranks"):
        text += " Ranks were assigned from relative evolutionary divergence (RED) bands."
    return text, []


_TEMPLATES: dict[str, Callable[[MethodsInfo], tuple[str, list[str]]]] = {
    "tree2tax": _tree2tax,
    "snptype": _snptype,
    "phylo": _phylo,
    "reads": _reads,
    "assemble": _assemble,
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
