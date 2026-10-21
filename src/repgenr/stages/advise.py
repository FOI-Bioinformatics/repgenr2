"""Advice on the phylogeny route for the genome set of a working directory.

Reads the genome set, the canonical record names, the workdir sketches, the
annotation index and the staged outgroup, and applies the recommendation
table of docs/choosing-tools.md ("By purpose and divergence"; the audit
behind it is docs/audit/phylo-methods-2026-10.md). Nothing is written to the
workdir and no tool runs except ``sourmash compare`` over sketches that
already exist. ``RECOMMENDATIONS`` holds the table as data, so a preset
option can expand the same rows.
"""

from __future__ import annotations

import logging
import shlex
import shutil
import statistics
import tempfile
import textwrap
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..core.contracts import record_name
from ..core.sketches import SKETCHES_DIR, workdir_provider
from ..maskers.gubbins import canonical_taxon, genus_spread, species_spread
from ..tree.medoid import ani_matrix

ADVISE_SCHEMA = "repgenr.phylo-advise/1"

# The dense sourmash ANI matrix is limited to this many genomes (the same
# limit as the medoid reference choice).
ANI_MAX_GENOMES = 5000
# Gubbins is for closely related isolates: it is refused above 12 % variable
# columns, and its run time grows quickly with the number of genomes.
MASK_MAX_GENOMES = 200
MASK_MIN_ANI = 97.0
OUTBREAK_MIN_ANI = 99.5
SPECIES_MIN_ANI = 95.0
GENUS_MIN_ANI = 80.0
# Used when the iqtree adapter cannot be read from the registry.
_IQTREE_LIMIT_FALLBACK = 500

_SPECIES_FLAGS = (
    "--msa-source",
    "snptype",
    "--snptyper",
    "simple",
    "--core-fraction",
    "0.95",
    "--treebuilder",
    "iqtree",
    "-B",
    "1000",
)
# `--partitioned` (a partition file for IQ-TREE) is not part of this release;
# it joins the genus and family rows when it ships.
_GENUS_FLAGS = (
    "--msa-source",
    "coregene",
    "--coregene",
    "ppanggolin",
    "--core-threshold",
    "0.95",
    "--max-core-genes",
    "300",
    "--gene-rank",
    "informative",
    "--treebuilder",
    "iqtree",
    "-B",
    "1000",
)

# The rows of the recommendation table, by preset name.
RECOMMENDATIONS: dict[str, dict[str, Any]] = {
    "outbreak": {
        "purpose": "Outbreak or clonal genotype tree",
        "route": "snptype",
        "flags": (
            "--msa-source",
            "snptype",
            "--snptyper",
            "ska2",
            "--treebuilder",
            "iqtree",
            "--aligner-arg",
            "model=GTR+G",
            "-B",
            "1000",
        ),
        "summary": (
            "Reference-free split k-mer SNPs (ska2) and an IQ-TREE tree. ska2 writes no "
            "constant-site counts, so branch lengths are inflated and phylo warns."
        ),
        "needs_annotation": False,
    },
    "species": {
        "purpose": "Genotype tree within a species",
        "route": "snptype",
        "flags": _SPECIES_FLAGS,
        "summary": (
            "Core SNPs against the medoid reference (simple typer), columns present in "
            "95 % of genomes, and an IQ-TREE tree with constant sites passed as -fconst."
        ),
        "needs_annotation": False,
    },
    "species-masked": {
        "purpose": "Genotype tree within a species, recombination masked",
        "route": "snptype",
        "flags": (*_SPECIES_FLAGS, "--mask", "gubbins"),
        "summary": (
            "As the species row, with recombination masked by Gubbins before the tree is built."
        ),
        "needs_annotation": False,
    },
    "genus": {
        "purpose": "Overview of species representatives in a genus",
        "route": "coregene",
        "flags": _GENUS_FLAGS,
        "summary": (
            "Core genes from PPanGGOLiN (in 95 % of genomes), capped at 300 ranked by "
            "parsimony-informative sites, and an IQ-TREE tree."
        ),
        "needs_annotation": True,
    },
    "family": {
        "purpose": "Overview across a family or wider",
        "route": "coregene",
        "flags": _GENUS_FLAGS,
        "summary": (
            "The genus route for now. DNA alignments saturate at this span; a protein "
            "alphabet for the core genes is planned. A skani tree (reliable to about "
            "80 % ANI) gives an overview without supports."
        ),
        "needs_annotation": True,
    },
    "quick": {
        "purpose": "Quick look at any set",
        "route": "alignment-free",
        "flags": ("--treebuilder", "skani"),
        "summary": (
            "An alignment-free skani tree without supports; in tree2tax use "
            "--collapse-diameter rather than --collapse-support."
        ),
        "needs_annotation": False,
    },
}


@dataclass(frozen=True)
class Recommendation:
    preset: str
    purpose: str
    route: str
    reasons: tuple[str, ...]
    command: tuple[str, ...]


@dataclass(frozen=True)
class Evidence:
    """What the workdir says about its genome set.

    ``n_genomes`` counts the representatives, or every genome with
    ``all_genomes``; ``n_species`` and ``n_genera`` count the canonical record
    names; ``ani_min`` and ``ani_median`` are percent over the pairs of the
    set (None without sketches or sourmash); ``has_annotations`` is true when
    the annotation index holds every genome of the tree, outgroup included.
    """

    workdir: str
    all_genomes: bool
    n_genomes: int
    n_species: int
    n_genera: int
    names_canonical: bool
    ani_min: float | None
    ani_median: float | None
    has_annotations: bool
    outgroup: bool
    notes: tuple[str, ...] = ()


# A logger that drops what the outgroup lookup reports: the advice says
# itself whether an outgroup is staged.
_SILENT = logging.getLogger("repgenr_advise_silent")
_SILENT.addHandler(logging.NullHandler())
_SILENT.propagate = False


def collect_evidence(ctx, all_genomes: bool, threads: int, logger: logging.Logger) -> Evidence:
    """Read the evidence from the workdir of ``ctx`` without writing to it."""
    from .phylo import _genome_set, resolve_outgroup_files

    notes: list[str] = []
    genomes = _genome_set(ctx, all_genomes)
    if not genomes:
        where = "genomes/" if all_genomes else "derep/representatives/"
        notes.append(f"no genomes found under {where}")
    names = [record_name(g) for g in genomes]
    canonical = bool(names) and all(canonical_taxon(n) is not None for n in names)
    outgroup_file, _leaf = resolve_outgroup_files(
        ctx.outgroup_dir, ctx.workdir / "outgroup_accession.txt", _SILENT
    )
    tree_names = names + ([record_name(outgroup_file)] if outgroup_file is not None else [])
    ani_min, ani_median, ani_note = _ani_summary(ctx, genomes, threads, logger)
    if ani_note:
        notes.append(ani_note)
    return Evidence(
        workdir=str(ctx.workdir),
        all_genomes=all_genomes,
        n_genomes=len(genomes),
        n_species=species_spread(names),
        n_genera=genus_spread(names),
        names_canonical=canonical,
        ani_min=ani_min,
        ani_median=ani_median,
        has_annotations=bool(names) and _annotated(ctx.annotations_dir, tree_names),
        outgroup=outgroup_file is not None,
        notes=tuple(notes),
    )


def _annotated(directory: Path, names: Sequence[str]) -> bool:
    """Whether the annotation index at ``directory`` holds every name with
    its GFF3 and one annotator (the check the core-gene route makes)."""
    from ..core.errors import UserInputError, WorkdirError
    from .coregene import annotated_set

    if not directory.is_dir():
        return False
    try:
        annotated_set(directory, list(names))
    except (UserInputError, WorkdirError):
        return False
    return True


def _sourmash_missing() -> str | None:
    """Why sourmash cannot run here, or None when it can."""
    from ..core.errors import MissingBinaryError
    from ..dereplicators.sourmash import SourmashDereplicator

    try:
        SourmashDereplicator().preflight()
    except MissingBinaryError as exc:
        return str(exc)
    return None


def _ani_summary(
    ctx, genomes: Sequence[Path], threads: int, logger: logging.Logger
) -> tuple[float | None, float | None, str | None]:
    """Minimum and median pairwise ANI (percent) over the existing workdir
    sketches, and a note when it could not be computed."""
    from ..core.errors import ToolExecutionError, WorkdirError
    from ..core.manifest import MANIFEST_FILENAME

    if len(genomes) < 2:
        return None, None, None
    if len(genomes) > ANI_MAX_GENOMES:
        return None, None, f"ANI not computed: {len(genomes)} genomes, above {ANI_MAX_GENOMES}"
    hint = "run `repgenr sketch` to enable the ANI evidence"
    # Opening the manifest would create one; without it there are no sketches.
    if not (ctx.workdir / SKETCHES_DIR).is_dir() or not (ctx.workdir / MANIFEST_FILENAME).is_file():
        return None, None, f"ANI not computed: no genome sketches ({hint})"
    sketches = workdir_provider(ctx, logger, threads, "phylo-advise", write=False)(genomes)
    missing = [g for g in genomes if g not in sketches]
    if missing:
        return (
            None,
            None,
            f"ANI not computed: {len(missing)} of {len(genomes)} genomes have no current "
            f"sketch ({hint})",
        )
    reason = _sourmash_missing()
    if reason is not None:
        return None, None, f"ANI not computed: {reason}"
    work = Path(tempfile.mkdtemp(prefix="repgenr-advise-"))
    try:
        _names, ani = ani_matrix(genomes, sketches, work, threads, logger)
    except (ToolExecutionError, WorkdirError) as exc:
        return None, None, f"ANI not computed: {exc}"
    finally:
        shutil.rmtree(work, ignore_errors=True)
    matrix = np.asarray(ani, dtype=float)
    pairs = matrix[np.triu_indices(matrix.shape[0], k=1)] * 100.0
    if pairs.size == 0:
        return None, None, None
    return float(pairs.min()), float(statistics.median(pairs.tolist())), None


def _iqtree_limit() -> int:
    try:
        from ..treebuilders.base import registry

        limit = registry.create("iqtree").capabilities.recommended_max_genomes
    except Exception:  # noqa: BLE001 - advice falls back to the documented figure
        return _IQTREE_LIMIT_FALLBACK
    return int(limit) if limit else _IQTREE_LIMIT_FALLBACK


def _pct(value: float) -> str:
    return f"{value:.2f} %"


def _choose(ev: Evidence) -> tuple[str, list[str]]:
    """The preset and the reasons for it, by the rules in order."""
    reasons: list[str] = []
    ani = ev.ani_min
    canonical = ev.names_canonical
    if ev.n_genomes == 0:
        reasons.append(
            "No genomes were found in the set; run dereplicate first, or pass --all-genomes."
        )
    if ani is not None and ani >= OUTBREAK_MIN_ANI:
        reasons.append(
            f"The minimum pairwise ANI is {_pct(ani)}, at or above {OUTBREAK_MIN_ANI:g} %: "
            "one lineage, where a reference-free SNP route resolves the most."
        )
        return "outbreak", reasons
    if (canonical and ev.n_species == 1) or (
        not canonical and ani is not None and ani >= SPECIES_MIN_ANI
    ):
        if canonical:
            reasons.append("The genome names give one species.")
            if ani is not None and ani < SPECIES_MIN_ANI:
                reasons.append(
                    f"The minimum pairwise ANI is {_pct(ani)}, below {SPECIES_MIN_ANI:g} %, "
                    "although the names give one species; check for mislabelled genomes."
                )
        else:
            reasons.append(
                f"The minimum pairwise ANI is {_pct(ani or 0.0)}, at or above "
                f"{SPECIES_MIN_ANI:g} %: one species."
            )
        if ev.n_genomes <= MASK_MAX_GENOMES and ani is not None and ani >= MASK_MIN_ANI:
            reasons.append(
                f"With {ev.n_genomes} genomes and a minimum ANI of {_pct(ani)}, Gubbins can "
                "mask recombination: it is for closely related isolates and is refused "
                "above 12 % variable columns."
            )
            return "species-masked", reasons
        if ani is None:
            why = "the divergence is unknown"
        elif ani < MASK_MIN_ANI:
            why = f"the minimum ANI is below {MASK_MIN_ANI:g} %"
        else:
            why = f"there are more than {MASK_MAX_GENOMES} genomes"
        reasons.append(f"Recombination masking is not suggested: {why}.")
        return "species", reasons
    if (ani is not None and ani < GENUS_MIN_ANI) or (canonical and ev.n_genera > 1):
        if canonical and ev.n_genera > 1:
            reasons.append(f"The genome names span {ev.n_genera} genera.")
        if ani is not None and ani < GENUS_MIN_ANI:
            reasons.append(f"The minimum pairwise ANI is {_pct(ani)}, below {GENUS_MIN_ANI:g} %.")
        reasons.append(
            "DNA alignments saturate at this span; protein core genes are planned and "
            "the genus route applies until then."
        )
        return "family", reasons
    if (canonical and ev.n_species > 1) or (ani is not None and ani >= GENUS_MIN_ANI):
        if canonical and ev.n_species > 1:
            reasons.append(f"The genome names span {ev.n_species} species of one genus.")
        if ani is not None:
            reasons.append(
                f"The minimum pairwise ANI is {_pct(ani)}, between {GENUS_MIN_ANI:g} and "
                f"{SPECIES_MIN_ANI:g} %: past the reach of a SNP route."
            )
        return "genus", reasons
    reasons.append(
        "Nothing is known about the divergence: the names are not canonical and there is "
        f"no ANI. Run `repgenr sketch -wd {ev.workdir}`, then phylo-advise again."
    )
    return "quick", reasons


def advise(evidence: Evidence) -> Recommendation:
    """The recommendation for ``evidence``: the preset, the reasons and the
    full ``repgenr phylo`` command."""
    preset, reasons = _choose(evidence)
    entry = RECOMMENDATIONS[preset]
    flags = list(entry["flags"])
    if entry["needs_annotation"] and not evidence.has_annotations:
        annotate = f"repgenr annotate -wd {evidence.workdir}"
        if evidence.all_genomes:
            annotate += " --all-genomes"
        reasons.append(
            f"The core-gene route needs an annotation of every genome; run `{annotate}` first."
        )
    limit = _iqtree_limit()
    if evidence.n_genomes > limit and "iqtree" in flags:
        flags = _to_fasttree(flags)
        reasons.append(
            f"{evidence.n_genomes} genomes are above the declared IQ-TREE limit of {limit}; "
            "use `--treebuilder fasttree` instead (SH-like supports, no -B needed)."
        )
    if not evidence.outgroup:
        flags += ["--root", "midpoint"]
        reasons.append(
            "No outgroup is staged, so the tree is rooted at its midpoint (add `--root "
            "midpoint`); with an outgroup the default roots on it."
        )
    command = ["repgenr", "phylo", "-wd", evidence.workdir]
    if evidence.all_genomes:
        command.append("--all-genomes")
    return Recommendation(
        preset=preset,
        purpose=entry["purpose"],
        route=entry["route"],
        reasons=tuple(reasons),
        command=(*command, *flags),
    )


def _to_fasttree(flags: list[str]) -> list[str]:
    """``flags`` with IQ-TREE swapped for FastTree: the IQ-TREE model and the
    bootstrap count are dropped (FastTree gives SH-like supports by default)."""
    out: list[str] = []
    i = 0
    while i < len(flags):
        flag = flags[i]
        value = flags[i + 1] if i + 1 < len(flags) else ""
        if flag == "--treebuilder" and value == "iqtree":
            out += [flag, "fasttree"]
        elif flag == "--aligner-arg" and value.startswith("model="):
            out += [flag, "gamma=true"]
        elif flag == "-B":
            pass
        else:
            out += [flag, value]
        i += 2
    return out


# Width the summary and the reasons are wrapped at; the command is not
# wrapped, so it can be copied as one line.
_WIDTH = 78


def _wrap(text: str, first: str = "", rest: str = "") -> str:
    """``text`` wrapped at _WIDTH, never between a number and its % sign."""
    filled = textwrap.fill(
        text.replace(" %", "\0%"), _WIDTH, initial_indent=first, subsequent_indent=rest
    )
    return filled.replace("\0", " ")


def format_text(rec: Recommendation, evidence: Evidence) -> str:
    """The advice as plain text: evidence, preset, reasons and command."""
    genomes = "all genomes" if evidence.all_genomes else "representatives"
    if evidence.names_canonical:
        genera = "genus" if evidence.n_genera == 1 else "genera"
        names = f"{evidence.n_species} species, {evidence.n_genera} {genera} (by name)"
    else:
        names = "names are not canonical; species not counted"
    if evidence.ani_min is not None and evidence.ani_median is not None:
        ani = f"min {_pct(evidence.ani_min)}, median {_pct(evidence.ani_median)}"
    else:
        ani = "unknown"
    lines = [
        f"Working directory: {evidence.workdir}",
        f"Genomes:      {evidence.n_genomes} ({genomes})",
        f"Taxonomy:     {names}",
        f"Pairwise ANI: {ani}",
        f"Annotations:  {'complete' if evidence.has_annotations else 'missing or incomplete'}",
        f"Outgroup:     {'staged' if evidence.outgroup else 'none'}",
    ]
    lines += [f"Note:         {n}" for n in evidence.notes]
    lines += [
        "",
        f"Preset: {rec.preset} ({rec.purpose})",
        f"Route:  {rec.route}",
        *_wrap(RECOMMENDATIONS[rec.preset]["summary"]).splitlines(),
        "",
        "Reasons:",
        *(_wrap(r, "  - ", "    ") for r in rec.reasons),
        "",
        "Command:",
        f"  {shlex.join(rec.command)}",
    ]
    return "\n".join(lines)


def to_json(rec: Recommendation, evidence: Evidence) -> dict[str, Any]:
    """The advice as one versioned JSON object."""
    from .. import __version__

    ev = asdict(evidence)
    ev["notes"] = list(evidence.notes)
    return {
        "schema": ADVISE_SCHEMA,
        "repgenr": __version__,
        "evidence": ev,
        "recommendation": {
            "preset": rec.preset,
            "purpose": rec.purpose,
            "route": rec.route,
            "summary": RECOMMENDATIONS[rec.preset]["summary"],
            "reasons": list(rec.reasons),
            "command": list(rec.command),
        },
    }
