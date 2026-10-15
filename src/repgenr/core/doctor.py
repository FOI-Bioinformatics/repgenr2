"""Workdir health checks: verify outputs against records (``repgenr doctor``).

``repgenr status`` reports what ``repgenr.yaml`` claims; ``doctor`` verifies the
claims against the filesystem and the manifest -- interrupted stages, missing or
corrupt genomes, manifest drift, representative/cluster mismatches, truncated
deliverables, unresolvable outgroups, leftover temp files, and stages whose
recorded input digests no longer match reality or whose declared deliverables
are missing (they will re-run).

Read-only: no output, log or record is written in the workdir. Opening the
WAL-mode manifest lets SQLite create or update its ``-shm``/``-wal``
companion files, which hold no data of their own.
"""

from __future__ import annotations

import functools
import logging
import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

from .config import CONFIG_FILENAME, Config
from .contracts import (
    CLUSTERS_TSV,
    CORE_SNP_FASTA,
    GENOME_STATUS_TSV,
    GENOMES_MAP_TSV,
    SELECTION_TSV,
    TREE2TAX_TSV,
    TREE_NWK,
    accession_from_filename,
    list_fasta,
    newick_is_complete,
    read_clusters,
    read_selection,
)
from .errors import WorkdirError
from .inputs import file_digest, inputs_digest, manifest_digest_for_stage
from .integrity import (
    check_genome_completeness,
    check_representatives_consistency,
    looks_like_fasta,
)
from .manifest import MANIFEST_FILENAME, Manifest

# The integrity guards log their refusal text as a warning when told to
# continue; doctor reports the same condition as a finding, so their log
# lines are dropped instead of reaching the console unformatted.
_LOG = logging.getLogger(__name__)
_LOG.addHandler(logging.NullHandler())
_LOG.propagate = False
_MAX_LISTED = 3  # examples shown per finding


@dataclass
class Finding:
    level: str  # "ok" | "warn" | "fail"
    area: str
    message: str


def diagnose(workdir: Path, *, quick: bool = False) -> list[Finding]:
    """Run every health check; never raises for a broken workdir.

    ``quick`` skips reading the first bytes of each genome file (the FASTA
    content check), the slowest check on large genome sets; links, the
    selection shortfall and untracked files are still checked.
    """
    workdir = Path(workdir)
    if not (workdir / CONFIG_FILENAME).exists():
        return [Finding("warn", "config", f"No RepGenR2 run found at {workdir}.")]

    findings: list[Finding] = []
    try:
        config = Config.load(workdir)
    except WorkdirError as exc:
        return [Finding("fail", "config", str(exc))]
    checks: tuple[Callable[[Path, Config], list[Finding]], ...] = (
        _check_stage_records,
        functools.partial(_check_genomes, quick=quick),
        _check_manifest_drift,
        _check_outgroup,
        _check_representatives,
        _check_tree,
        _check_phylo_stamp_in_snp,
        _check_tree2tax_pair,
        _check_stale_inputs,
        _check_deliverables,
        _check_leftovers,
    )
    for check in checks:
        try:
            findings.extend(check(workdir, config))
        except Exception as exc:  # a broken artifact must not abort the report
            findings.append(
                Finding(
                    "fail",
                    getattr(check, "func", check).__name__.removeprefix("_check_"),
                    f"Check could not complete: {exc}",
                )
            )
    if not any(f.level == "fail" for f in findings):
        findings.append(Finding("ok", "summary", "No failures detected."))
    return findings


def _examples(names: list[str]) -> str:
    shown = ", ".join(names[:_MAX_LISTED])
    more = f" (+{len(names) - _MAX_LISTED} more)" if len(names) > _MAX_LISTED else ""
    return shown + more


def _genome_set_stages(config: Config) -> tuple[str, str]:
    """(stage that wrote selection.tsv, stage that placed genomes/) for advice."""
    for name in ("ingest", "vgenome", "assemble"):
        if name in config.stages:
            return name, name
    return "metadata", "genome"


def holds_outputs(workdir: Path) -> bool:
    """Whether the workdir holds pipeline outputs (genomes, selection, derep, tree)."""
    return any((workdir / name).exists() for name in ("genomes", SELECTION_TSV, "derep", "tree"))


def _check_stage_records(workdir: Path, config: Config) -> list[Finding]:
    out: list[Finding] = []
    if not config.stages and holds_outputs(workdir):
        out.append(
            Finding(
                "warn",
                "config",
                f"{CONFIG_FILENAME} records no stage, but the workdir holds outputs; the "
                "record was emptied or replaced, and every stage will re-run.",
            )
        )
    for name, record in config.stages.items():
        if not record.interrupted:
            out.append(Finding("ok", name, f"completed {record.completed}"))
            if not record.fingerprint:
                out.append(
                    Finding(
                        "warn",
                        name,
                        "completed without a resume fingerprint (recorded by an older "
                        "version or restored by derep-stock unpack); its next invocation "
                        "recomputes it.",
                    )
                )
        else:
            out.append(
                Finding(
                    "fail",
                    name,
                    "started a run and did not finish; outputs may be partial; see "
                    "repgenr.log and re-run the stage.",
                )
            )
    return out


def _check_genomes(workdir: Path, config: Config, *, quick: bool = False) -> list[Finding]:
    genomes_dir = workdir / "genomes"
    writer = _genome_set_stages(config)[1]
    out: list[Finding] = []
    shortfall = check_genome_completeness(genomes_dir, workdir, logger=_LOG, allow_incomplete=True)
    if shortfall:
        out.append(
            Finding(
                "fail",
                "genomes",
                f"{len(shortfall)} selected genome(s) missing from {genomes_dir} "
                f"(e.g. {_examples(shortfall)}); re-run {writer}.",
            )
        )
    entries = list_fasta(genomes_dir)
    dangling = [p for p in entries if p.is_symlink() and not p.exists()]
    if dangling:
        out.append(
            Finding(
                "fail",
                "genomes",
                f"{len(dangling)} link(s) under {genomes_dir} point at files that no longer "
                f"exist (e.g. {dangling[0].name} -> {os.readlink(dangling[0])}); the linked "
                "source was moved or deleted. Re-run ingest from its new location, or stage "
                "copies with ingest --copy.",
            )
        )
    bad = (
        [] if quick else [p.name for p in entries if p not in dangling and not looks_like_fasta(p)]
    )
    if bad:
        out.append(
            Finding(
                "fail",
                "genomes",
                f"{len(bad)} file(s) under {genomes_dir} are not FASTA "
                f"(e.g. {_examples(bad)}); delete them and re-run {writer}.",
            )
        )
    untracked = _untracked_genomes(workdir, entries)
    if untracked:
        out.append(
            Finding(
                "warn",
                "genomes",
                f"{len(untracked)} file(s) under {genomes_dir} are not in {SELECTION_TSV} "
                f"(e.g. {_examples(untracked)}); dereplicate would include them. Remove "
                f"them, or re-run {_genome_set_stages(config)[0]} to select them.",
            )
        )
    if not shortfall and not bad and not dangling and genomes_dir.exists():
        count = len(list_fasta(genomes_dir))
        message = (
            f"{count} genome file(s) present (content not read: --quick)"
            if quick
            else f"{count} genome file(s) look sound"
        )
        out.append(Finding("ok", "genomes", message))
    return out


def _untracked_genomes(workdir: Path, entries: list[Path]) -> list[str]:
    """Genome files under genomes/ that selection.tsv does not list."""
    selection = workdir / SELECTION_TSV
    if not selection.exists():
        return []
    selected = {row.filename for row in read_selection(selection)}
    return sorted(p.name for p in entries if p.name not in selected)


def _check_manifest_drift(workdir: Path, config: Config) -> list[Finding]:
    selection = workdir / SELECTION_TSV
    manifest_path = workdir / MANIFEST_FILENAME
    if not selection.exists() or not manifest_path.exists():
        return []
    selected = {row.accession for row in read_selection(selection)}
    manifest = Manifest.open_readonly(manifest_path)
    try:
        recorded = {g.accession for g in manifest.all_genomes(include_outgroup=True)}
    finally:
        manifest.close()
    extra = sorted(recorded - selected)
    missing = sorted(selected - recorded)
    writer = _genome_set_stages(config)[0]
    out: list[Finding] = []
    if extra:
        out.append(
            Finding(
                "fail",
                "manifest",
                f"manifest holds {len(extra)} genome(s) not in selection.tsv "
                f"(e.g. {_examples(extra)}); re-run {writer} to reconcile.",
            )
        )
    if missing:
        out.append(
            Finding(
                "fail",
                "manifest",
                f"manifest is missing {len(missing)} selected genome(s) "
                f"(e.g. {_examples(missing)}); re-run {writer}.",
            )
        )
    if not extra and not missing:
        out.append(Finding("ok", "manifest", "manifest matches selection.tsv"))
    return out


def _check_outgroup(workdir: Path, config: Config) -> list[Finding]:
    acc_file = workdir / "outgroup_accession.txt"
    if not acc_file.exists():
        return []
    accession = acc_file.read_text(encoding="utf-8").strip()
    if not accession:
        return []
    outgroup_dir = workdir / "outgroup"
    # The rule phylo's resolve_outgroup_files applies: genome FASTA files only.
    candidates = list_fasta(outgroup_dir)
    resolved = any(
        accession_from_filename(f.name) == accession or accession in f.name for f in candidates
    )
    if resolved:
        return [Finding("ok", "outgroup", f"outgroup {accession} resolves")]
    return [
        Finding(
            "warn",
            "outgroup",
            f"recorded outgroup {accession} has no matching file under {outgroup_dir}; "
            "phylo/tree2tax will proceed unrooted.",
        )
    ]


def _check_representatives(workdir: Path, config: Config) -> list[Finding]:
    derep_dir = workdir / "derep"
    clusters = derep_dir / CLUSTERS_TSV
    if not clusters.exists():
        return []
    out: list[Finding] = []
    shortfall = check_representatives_consistency(
        derep_dir / "representatives", clusters, logger=_LOG, allow_incomplete=True
    )
    if shortfall:
        out.append(
            Finding(
                "fail",
                "dereplicate",
                f"{len(shortfall)} representative(s) listed in {CLUSTERS_TSV} are "
                f"absent on disk (e.g. {_examples(shortfall)}); re-run dereplicate.",
            )
        )
    extras = sorted(
        {p.name for p in list_fasta(derep_dir / "representatives")} - set(read_clusters(clusters))
    )
    if extras:
        out.append(
            Finding(
                "warn",
                "dereplicate",
                f"{len(extras)} file(s) under representatives/ are not in "
                f"{CLUSTERS_TSV} (e.g. {_examples(extras)}).",
            )
        )
    if not (derep_dir / GENOME_STATUS_TSV).exists():
        out.append(
            Finding(
                "fail",
                "dereplicate",
                f"{CLUSTERS_TSV} present but {GENOME_STATUS_TSV} missing; the "
                "dereplicate stage likely crashed mid-write -- re-run it.",
            )
        )
    if not out:
        out.append(Finding("ok", "dereplicate", "representatives match clusters.tsv"))
    return out


def _check_tree(workdir: Path, config: Config) -> list[Finding]:
    tree = workdir / "tree" / TREE_NWK
    if not tree.exists():
        return []
    if not newick_is_complete(tree.read_text(encoding="utf-8")):
        return [
            Finding(
                "fail",
                "phylo",
                f"{tree} is empty, truncated or holds more than one tree (it must hold one "
                "tree ended by ';'); re-run phylo.",
            )
        ]
    return [Finding("ok", "phylo", f"{TREE_NWK} looks like a complete tree")]


def _check_phylo_stamp_in_snp(workdir: Path, config: Config) -> list[Finding]:
    """A reuse stamp in snp/, left by phylo before its typing pass moved to tree/msa/.

    That phylo version typed into snp/, the snptype stage's directory. Only
    when the stamp's digest still matches snp/core_snp.fasta are the tables
    there phylo's (outgroup included) rather than the snptype stage's; a
    snptype run since then replaced them, and the stamp is then harmless.
    phylo no longer reads or writes snp/.
    """
    import json

    stamp = workdir / "snp" / "msa_source.json"
    if not stamp.is_file():
        return []
    try:
        recorded = json.loads(stamp.read_text(encoding="utf-8")).get("artifact_digest")
    except (OSError, ValueError, AttributeError):
        return []
    if not recorded or recorded != file_digest(workdir / "snp" / CORE_SNP_FASTA):
        return []
    return [
        Finding(
            "warn",
            "snptype",
            "snp/msa_source.json was left by an earlier phylo typing pass (now under "
            "tree/msa/); the tables in snp/ may be phylo's. Delete the stamp, and run "
            "'repgenr --force snptype' if the snptype tables are needed.",
        )
    ]


def _check_tree2tax_pair(workdir: Path, config: Config) -> list[Finding]:
    t2t = workdir / TREE2TAX_TSV
    gmap = workdir / GENOMES_MAP_TSV
    if t2t.exists() and gmap.exists():
        problem = _tree2tax_tables_problem(t2t, gmap)
        if problem:
            return [
                Finding(
                    "fail",
                    "tree2tax",
                    f"{problem}; the files were truncated or edited -- re-run "
                    "repgenr --force tree2tax.",
                )
            ]
        findings = [Finding("ok", "tree2tax", "deliverable pair present and consistent")]
        top = _tree2tax_top_node(t2t)
        if top != "root":
            findings.append(
                Finding(
                    "warn",
                    "tree2tax",
                    f"the top node is named {top!r}; FlexTaxD validates its tree from a "
                    "node named 'root' (tree2tax --root-name).",
                )
            )
        return findings
    if t2t.exists() == gmap.exists():
        return []
    missing = GENOMES_MAP_TSV if t2t.exists() else TREE2TAX_TSV
    return [
        Finding(
            "fail",
            "tree2tax",
            f"deliverables are mismatched: {missing} is missing while its partner "
            "exists; the tree2tax stage likely crashed mid-write -- re-run it.",
        )
    ]


def _tree2tax_top_node(t2t: Path) -> str:
    """The one node of a consistent tree2tax.tsv that is never a child."""
    edges = [line.split("\t") for line in t2t.read_text(encoding="utf-8").splitlines()[1:] if line]
    children = {child for child, _ in edges}
    (top,) = {parent for _, parent in edges} - children
    return top


def _tree2tax_tables_problem(t2t: Path, gmap: Path) -> str | None:
    """Why tree2tax.tsv and genomes_map.tsv do not describe the same leaves.

    tree2tax.tsv holds a ``child``/``parent`` header and one edge per row;
    its leaves are the children that are never a parent. genomes_map.tsv
    maps every leaf to itself (and members to their leaf), so the two leaf
    sets must be equal.
    """
    lines = t2t.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split("\t") != ["child", "parent"]:
        return f"{TREE2TAX_TSV} lacks its child/parent header"
    edges = [line.split("\t") for line in lines[1:] if line]
    if not edges or any(len(edge) != 2 for edge in edges):
        return f"{TREE2TAX_TSV} holds no edges, or a row without two columns"
    if any(not child or not parent for child, parent in edges):
        return f"{TREE2TAX_TSV} holds an empty node name"
    parent_of: dict[str, str] = {}
    for child, parent in edges:
        if parent_of.setdefault(child, parent) != parent:
            return f"{TREE2TAX_TSV} lists {child!r} under two parents"
    parents = {parent for _, parent in edges}
    tops = sorted(parents - set(parent_of))
    if len(tops) != 1:
        return f"{TREE2TAX_TSV} has {len(tops)} top nodes instead of one ({_examples(tops)})"
    leaves = {child for child, _ in edges if child not in parents}
    rows = [line.split("\t") for line in gmap.read_text(encoding="utf-8").splitlines() if line]
    if not rows or any(len(row) != 2 for row in rows):
        return f"{GENOMES_MAP_TSV} is empty, or holds a row without two columns"
    mapped = {leaf for _, leaf in rows}
    if leaves != mapped:
        diff = sorted(leaves ^ mapped)
        return (
            f"{TREE2TAX_TSV} and {GENOMES_MAP_TSV} name different leaves "
            f"({len(diff)} in one only, e.g. {_examples(diff)})"
        )
    return None


def _layout(workdir: Path) -> SimpleNamespace:
    """Stand-in for WorkdirContext's path layout (no manifest is opened)."""
    return SimpleNamespace(
        workdir=workdir,
        genomes_dir=workdir / "genomes",
        outgroup_dir=workdir / "outgroup",
        derep_dir=workdir / "derep",
        representatives_dir=workdir / "derep" / "representatives",
        snp_dir=workdir / "snp",
        tree_dir=workdir / "tree",
    )


# Records exempt from the staleness checks: a derep-stock record logs the
# last pack or unpack, whose inputs are the live derep outputs; a later
# dereplicate changes them by design, and the stored run is unaffected.
_STALENESS_EXEMPT = frozenset({"derep_stock"})


def _missing_by_stage(workdir: Path, config: Config) -> dict[str, list[str]]:
    """Completed stage -> its declared deliverables that are missing.

    Uses the same table (STAGE_DELIVERABLES) and the same workdir-relative
    names as the resume check in the stage harness.
    """
    from ..cli.base import deliverable_label, missing_deliverables  # deferred: core<-cli

    ctx = _layout(workdir)
    out: dict[str, list[str]] = {}
    for name, record in config.stages.items():
        if record.interrupted:
            continue
        params = SimpleNamespace(**record.params)
        labels = [deliverable_label(workdir, p) for p in missing_deliverables(ctx, name, params)]
        if labels:
            out[name] = labels
    return out


def _changed_by_stage(workdir: Path, config: Config) -> dict[str, list[str]]:
    """Completed stage -> its recorded inputs whose digest no longer matches."""
    from ..cli.base import _MANIFEST_INPUT_STAGES, STAGE_INPUTS  # deferred: core<-cli

    ctx = _layout(workdir)
    out: dict[str, list[str]] = {}
    for name, record in config.stages.items():
        if record.interrupted or not record.inputs or name in _STALENESS_EXEMPT:
            continue
        spec = STAGE_INPUTS.get(name)
        if spec is None:
            continue
        params = SimpleNamespace(**record.params)
        digests = inputs_digest(workdir, spec(ctx, params))
        if name in _MANIFEST_INPUT_STAGES:
            manifest_path = workdir / MANIFEST_FILENAME
            if manifest_path.exists():
                manifest = Manifest.open_readonly(manifest_path)
                try:
                    digests["manifest"] = manifest_digest_for_stage(name, manifest)
                finally:
                    manifest.close()
        changed = sorted(
            key for key in {*record.inputs, *digests} if record.inputs.get(key) != digests.get(key)
        )
        if changed:
            out[name] = changed
    return out


def stale_stages(workdir: Path, config: Config) -> dict[str, str]:
    """Completed stages that will re-run on their next invocation, with why.

    A recorded input changed since completion, or a declared deliverable is
    missing. ``status`` shows these as stale and ``doctor`` warns about them,
    from the same two checks, so the commands agree.
    """
    reasons: dict[str, list[str]] = {}
    for name, keys in _changed_by_stage(workdir, config).items():
        reasons.setdefault(name, []).append(f"input changed: {_examples(keys)}")
    for name, labels in _missing_by_stage(workdir, config).items():
        reasons.setdefault(name, []).append(f"missing: {_examples(labels)}")
    return {name: "; ".join(parts) for name, parts in reasons.items()}


def _check_deliverables(workdir: Path, config: Config) -> list[Finding]:
    """Completed stages whose declared deliverables are missing (they will re-run)."""
    out: list[Finding] = []
    for name, labels in _missing_by_stage(workdir, config).items():
        what = (
            f"deliverable {labels[0]} missing"
            if len(labels) == 1
            else f"{len(labels)} deliverables missing ({_examples(labels)})"
        )
        out.append(Finding("warn", name, f"{what}; the stage will re-run on its next invocation."))
    return out


def _check_stale_inputs(workdir: Path, config: Config) -> list[Finding]:
    """Completed stages whose recorded input digests no longer match reality."""
    return [
        Finding(
            "warn",
            name,
            f"input(s) changed since completion ({_examples(changed)}); "
            "the stage will re-run on its next invocation.",
        )
        for name, changed in _changed_by_stage(workdir, config).items()
    ]


def _check_leftovers(workdir: Path, config: Config) -> list[Finding]:
    leftovers = sorted(
        str(p.relative_to(workdir))
        for pattern in ("*.tmp", "*.part")
        for p in workdir.rglob(pattern)
        # exFAT keeps a ._ AppleDouble companion beside each file; it is not
        # a second leftover.
        if "scratch" not in p.parts and not p.name.startswith("._")
    )
    if not leftovers:
        return []
    return [
        Finding(
            "warn",
            "leftovers",
            f"{len(leftovers)} temp file(s) from an interrupted write "
            f"(e.g. {_examples(leftovers)}); safe to delete.",
        )
    ]
