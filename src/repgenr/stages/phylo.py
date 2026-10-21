"""Phylogenetics stage: compose an MSA source with a tree builder.

Three orthogonal choices:
  * genome set    -- dereplicated representatives (default) or all genomes
  * MSA source    -- a whole-genome aligner OR a SNP typer's core-SNP alignment
                     (skipped entirely for alignment-free tree builders)
  * tree builder  -- iqtree / fasttree / veryfasttree / raxmlng (MSA) or
                     mashtree / sourmash / skani / ska2 (alignment-free)
                     (alignment-free)

The outgroup genome is added to the input set here and passed to builders that
can root (iqtree, raxmlng); the others emit an unrooted tree. Whatever the
builder, the tree is then rooted here (--root): on the branch leading to the
outgroup, at the midpoint, or not at all. A rerooted tree is written to
tree/tree.nwk and the builder's own output is kept as tree/tree.unrooted.nwk.

The SNP typing pass (--msa-source snptype) writes under tree/msa/, not snp/:
snp/ belongs to the snptype stage, so the two never replace each other's
tables and each stage's record describes what is on disk.
"""

from __future__ import annotations

import hashlib
import json
import logging
import mmap
import os
import re
import shutil
from collections.abc import Collection, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from ..aligners.base import AlignParams
from ..aligners.base import registry as aligner_registry
from ..core.containers import result_env_fragment
from ..core.context import WorkdirContext
from ..core.contracts import (
    CLUSTERS_TSV,
    CORE_GENE_FASTA,
    CORE_SNP_FASTA,
    MSA_FASTA,
    PHYLO_MSA_DIR,
    TREE_NWK,
    TREE_UNROOTED_NWK,
    accession_from_filename,
    atomic_path,
    list_fasta,
    record_name,
    taxonomic_spread,
)
from ..core.errors import PluginError, UserInputError, WorkdirError
from ..core.inputs import file_digest, paths_stat_digest
from ..core.integrity import check_genome_completeness, check_representatives_consistency
from ..core.plugins import ToolCapabilities, auto_select, scale_warning, warn_ignored_params
from ..core.process import remove_tree, stage_plain_inputs
from ..core.sketches import (
    SketchProvider,
    adapter_sketches,
    directory_provider,
    workdir_provider,
)
from ..coregenes.select import GeneSelection
from ..tree.alignment import ALIGNMENT_STATS, copy_alignment_stats, read_constant_sites
from ..tree.medoid import (
    REFERENCE_RULES as REFERENCE_RULES_ALL,
)
from ..tree.medoid import (
    ReferenceChoice,
    clear_reference_note,
    note_reference,
    read_reference_choice,
    resolve_reference,
)
from ..tree.rooting import (
    ROOT_METHODS,
    ROOT_NONE,
    ROOT_OUTGROUP,
    midpoint_root,
    read_newick,
    root_on_outgroup,
    write_newick,
)
from ..treebuilders.base import InputKind, TreeParams
from ..treebuilders.base import registry as treebuilder_registry

# Re-exported so the data-channel step helpers stay importable from this module.
__all__ = ["PhyloParams", "PhyloBuildParams", "build_tree", "phylo_build", "run", "list_fasta"]

# Keys the phylo stage reads itself; never forwarded to an adapter.
_STAGE_EXTRA_KEYS = frozenset(
    {
        "mask",
        "core_fraction",
        "core_threshold",
        "clean_mode",
        "filtered",
        # The core-gene selection (coregenes/select.py); read by _gene_selection.
        "max_core_genes",
        "gene_rank",
        "min_gene_length",
        "max_gene_gaps",
    }
)
# Extras read by the core-gene pass only; the fingerprint drops them elsewhere.
_COREGENE_ONLY_EXTRA = frozenset(
    {"clean_mode", "filtered", "max_core_genes", "gene_rank", "min_gene_length", "max_gene_gaps"}
)
# Bump when the rule that lays the selected genes out in the core-gene
# concatenation changes, so an alignment built by the old rule is rebuilt.
_COREGENE_LAYOUT = 1


def _core_fraction(extra: Mapping[str, object]) -> float | None:
    """The --core-fraction carried in the stage extras, or None."""
    raw = extra.get("core_fraction")
    if raw is None or str(raw).strip() == "":
        return None
    try:
        value = float(str(raw))
    except ValueError:
        raise UserInputError(
            f"--core-fraction must be a number between 0 and 1, got {raw!r}"
        ) from None
    if not 0 <= value <= 1:
        raise UserInputError(f"--core-fraction must be between 0 and 1, got {raw!r}")
    return value


def _core_threshold(extra: Mapping[str, object]) -> float:
    raw = extra.get("core_threshold", 0.95)
    try:
        value = float(str(raw))
    except ValueError:
        raise UserInputError(
            f"--core-threshold must be a number within (0, 1], got {raw!r}"
        ) from None
    if not 0.0 < value <= 1.0:
        raise UserInputError(f"--core-threshold must be within (0, 1], got {raw!r}")
    return value


def _gene_selection(extra: Mapping[str, object]) -> GeneSelection:
    """The core-gene selection carried in the stage extras; the defaults keep
    every gene, ranked by occupancy. Values are validated by GeneSelection."""

    def raw(key: str) -> str | None:
        value = extra.get(key)
        return None if value is None or str(value).strip() == "" else str(value).strip()

    def integer(key: str, default: int) -> int:
        value = raw(key)
        if value is None:
            return default
        try:
            return int(value)
        except ValueError:
            raise UserInputError(f"{key} must be an integer, got {value!r}") from None

    gaps = raw("max_gene_gaps")
    try:
        max_gaps = float(gaps) if gaps is not None else None
    except ValueError:
        raise UserInputError(
            f"max_gene_gaps must be a number between 0 and 1, got {gaps!r}"
        ) from None
    return GeneSelection(
        max_genes=integer("max_core_genes", 0),
        rank=raw("gene_rank") or "occupancy",
        min_length=integer("min_gene_length", 0),
        max_gaps=max_gaps,
    )


def _filtered(extra: Mapping[str, object]) -> bool:
    from ..core.plugins import parse_extra_bool

    return parse_extra_bool(extra, "filtered", True)


def _adapter_extra(extra: dict) -> dict:
    return {k: v for k, v in extra.items() if k not in _STAGE_EXTRA_KEYS}


def _warn_stage_extras(
    caps: Sequence[ToolCapabilities], extra: Mapping[str, object], logger: logging.Logger
) -> None:
    """Warn once about extras that none of the participating adapters reads.

    The phylo stage drives two adapters from one extras dict (an MSA source and
    a tree builder), so warning per adapter reports a key the other one consumes.
    A key is unused only when no participating tool declares it.
    """
    known: set[str] = set()
    for cap in caps:
        known |= set(cap.accepted_extras) | set(cap.default_params)
    unread = sorted(set(extra) - known)
    if unread:
        logger.warning(
            "No selected phylo tool (%s) reads extra parameter(s): %s",
            ", ".join(cap.name for cap in caps),
            ", ".join(unread),
        )


def _msa_source_capabilities(params: PhyloParams) -> ToolCapabilities | None:
    """Capabilities of the adapter that will build the MSA (aligner or SNP typer).

    None when the source is unknown; ``_build_msa`` raises for that case.
    """
    if params.msa_source == "aligner":
        return aligner_registry.get(params.aligner).capabilities
    if params.msa_source == "snptype":
        from ..snptypers.base import registry as snp_registry

        return snp_registry.get(params.snptyper).capabilities
    if params.msa_source == "coregene":
        from ..coregenes.base import registry as coregene_registry

        return coregene_registry.get(params.coregene).capabilities
    return None


@dataclass
class PhyloParams:
    treebuilder: str = "iqtree"
    msa_source: str = "aligner"  # aligner | snptype | coregene
    aligner: str = "progressivemauve"
    snptyper: str = "simple"
    # The core-gene builder of --msa-source coregene; its settings travel in
    # extra as core_threshold, clean_mode and filtered.
    coregene: str = "panaroo"
    all_genomes: bool = False
    no_outgroup: bool = False
    threads: int = 16
    bootstrap: int = 0
    reference: str | None = None
    allow_incomplete: bool = False
    extra: dict = field(default_factory=dict)
    # outgroup | midpoint | none; None picks outgroup when an outgroup is
    # staged (and not --no-outgroup), else none.
    root: str | None = None


@dataclass
class PhyloDirs:
    """Output directories for the stateless tree-building core."""

    tree_dir: Path
    align_dir: Path
    # Where the SNP typing pass writes its alignment and stamp (tree/msa/).
    msa_dir: Path
    scratch_dir: Path
    # The annotate stage's cache, read by --msa-source coregene.
    annotations_dir: Path | None = None


@dataclass
class PhyloOutcome:
    tree: Path
    treebuilder: str
    versions: dict[str, str]
    outgroup_leaf: str | None
    # The rooting applied: outgroup, midpoint or none.
    root: str = ROOT_NONE


def check_root_request(root: str | None, no_outgroup: bool) -> None:
    """Refuse a --root value, or a combination, that cannot be carried out."""
    if root is None:
        return
    if root not in ROOT_METHODS:
        raise UserInputError(f"--root must be one of {', '.join(ROOT_METHODS)}; got '{root}'.")
    if root == ROOT_OUTGROUP and no_outgroup:
        raise UserInputError(
            "--root outgroup needs the outgroup that --no-outgroup leaves out; "
            "use --root midpoint or --root none."
        )


def resolve_root(root: str | None, outgroup_leaf: str | None, no_outgroup: bool) -> str:
    """The rooting to apply, given the request and whether an outgroup is staged."""
    check_root_request(root, no_outgroup)
    has_outgroup = outgroup_leaf is not None and not no_outgroup
    if root is None:
        return ROOT_OUTGROUP if has_outgroup else ROOT_NONE
    if root == ROOT_OUTGROUP and not has_outgroup:
        raise UserInputError(
            "--root outgroup needs an outgroup, but none is staged (no "
            "outgroup_accession.txt, or its genome is not in outgroup/). Stage one "
            "with the genome stage, or use --root midpoint or --root none."
        )
    return root


def apply_root(
    tree: Path,
    root: str,
    outgroup_leaf: str | None,
    treebuilder: str,
    logger: logging.Logger,
) -> None:
    """Root the tree file in place; keep the builder's output when it changes.

    ``none`` leaves the file as the builder wrote it. A tree that is already
    rooted on the outgroup (a bifurcating root with the outgroup as one child)
    is left alone as well. Otherwise the builder's output is copied to
    tree.unrooted.nwk beside it and the rooted tree replaces tree.nwk.
    """
    if root == ROOT_NONE:
        return
    parsed = read_newick(tree.read_text(encoding="utf-8"))
    if root == ROOT_OUTGROUP:
        assert outgroup_leaf is not None  # resolve_root guarantees it
        if not root_on_outgroup(parsed, outgroup_leaf, source=tree):
            logger.info(
                "The %s tree is already rooted on the outgroup %s; left as written",
                treebuilder,
                outgroup_leaf,
            )
            return
        what = f"on the branch to the outgroup {outgroup_leaf}"
    else:
        midpoint_root(parsed, source=tree)
        what = "at the midpoint"
    unrooted = tree.with_name(TREE_UNROOTED_NWK)
    with atomic_path(unrooted) as tmp:
        shutil.copy2(tree, tmp)
    with atomic_path(tree) as tmp:
        tmp.write_text(write_newick(parsed), encoding="utf-8")
    logger.info(
        "Rooted the %s tree %s; the builder's output is kept as %s",
        treebuilder,
        what,
        unrooted.name,
    )


# Stamp written beside a built MSA, so a later run that only changes the tree
# builder (or the bootstrap) can reuse it instead of aligning or SNP-calling
# again. The stage fingerprint cannot do this: it covers the whole stage.
MSA_STAMP = "msa_source.json"
# 2: the snippy typer names its reference record by genome (was "Reference").
# 3: parsnp and cactus records carry genome stems (were file names with
# '.ref', and '.' replaced by '_'); an older alignment would keep them.
# The SNP typing pass moved from snp/ to tree/msa/ without a version change: a
# stamp left in snp/ is never read, so such a workdir is typed once more.
_MSA_STAMP_VERSION = 3
# The typing pass writes alignment_stats.json (constant-site counts) beside
# its alignment; the key of a snptype-source alignment carries this so one
# typed earlier is typed once more, while aligner alignments are untouched.
_SNPTYPE_STATS_VERSION = 1
# Scratch subdirectory of the SNP typing pass; the snptype stage uses
# scratch/snptype, and each stage clears its own at the start.
_TYPING_SCRATCH = "phylo_snptype"
# Scratch subdirectory for decompressed copies of gzipped genomes, given to an
# aligner or genome-input tree builder that cannot read gzip; removed after use.
_STAGED_INPUTS = "phylo_inputs"


def _msa_artifact(dirs: PhyloDirs, params: PhyloParams) -> Path:
    """Path the MSA source writes, by source."""
    if params.msa_source == "aligner":
        return dirs.align_dir / MSA_FASTA
    if params.msa_source == "coregene":
        return dirs.msa_dir / CORE_GENE_FASTA
    return dirs.msa_dir / CORE_SNP_FASTA


def _msa_key(
    genomes: Sequence[Path],
    outgroup_file: Path | None,
    params: PhyloParams,
    reference: Path | None = None,
    annotations: str | None = None,
) -> str:
    """Hash of everything that determines the MSA, and nothing else.

    Deliberately excludes the tree builder, the bootstrap and the thread count:
    those change the tree, not the alignment the tree is built from.
    """
    inputs = list(genomes) + ([outgroup_file] if outgroup_file is not None else [])
    payload = {
        "v": _MSA_STAMP_VERSION,
        "msa_source": params.msa_source,
        "aligner": params.aligner if params.msa_source == "aligner" else None,
        "snptyper": params.snptyper if params.msa_source == "snptype" else None,
        # The genome the alignment is projected onto, not the rule that chose
        # it: a medoid that changes with the genome set rebuilds the alignment.
        "reference": reference.name if reference is not None else params.reference,
        "all_genomes": params.all_genomes,
        "no_outgroup": params.no_outgroup,
        "mask": str(params.extra.get("mask", "none")),
        # Read by the typing pass only; an aligner alignment is unchanged by it.
        "core_fraction": _core_fraction(params.extra) if params.msa_source == "snptype" else None,
        "stats": _SNPTYPE_STATS_VERSION if params.msa_source == "snptype" else None,
        "extra": {k: str(v) for k, v in sorted(_adapter_extra(params.extra).items())},
        "inputs": paths_stat_digest(inputs),
        "env": result_env_fragment(),
    }
    if params.msa_source == "coregene":
        # Added for this source only, so an aligner or SNP alignment made
        # before the core-gene route existed keeps its key and is reused.
        # It holds the builder, its settings and the annotations the pass
        # read (the index rows carry the GFF digests and the annotator).
        payload["coregene"] = {
            "tool": params.coregene,
            "core_threshold": _core_threshold(params.extra),
            "clean_mode": _coregene_setting(params, "clean_mode"),
            "filtered": _coregene_setting(params, "filtered"),
            "annotations": annotations,
            "layout": _COREGENE_LAYOUT,
            "selection": asdict(_gene_selection(params.extra)),
        }
    blob = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _coregene_ignored(tool: str) -> frozenset[str]:
    """The core-gene parameters ``tool`` declares it does not read."""
    from ..coregenes.base import registry as coregene_registry

    try:
        return frozenset(coregene_registry.get(tool).capabilities.ignored_params)
    except Exception:  # noqa: BLE001 - an unknown builder is refused when it runs
        return frozenset()


def _coregene_setting(params: PhyloParams, key: str) -> str | bool | None:
    """``clean_mode`` or ``filtered`` as the core-gene builder reads it, or None
    when the builder ignores it, so changing a setting the builder never reads
    neither rebuilds the alignment nor reruns the stage."""
    if key in _coregene_ignored(params.coregene):
        return None
    if key == "clean_mode":
        return str(params.extra.get("clean_mode", "strict"))
    return _filtered(params.extra)


def fingerprint_params(payload: dict) -> dict:
    """The phylo parameters as the resume fingerprint sees them.

    ``clean_mode`` and ``filtered`` are removed from ``extra`` when the
    core-gene builder ignores them (PPanGGOLiN, SCARAP, OrthoFinder) or when
    no core-gene pass runs, so a changed value reruns nothing it cannot
    change. The gene-selection keys are removed when no core-gene pass runs.
    """
    extra = payload.get("extra")
    if not isinstance(extra, dict) or not (_COREGENE_ONLY_EXTRA & extra.keys()):
        return payload
    drop: set[str] = set(_COREGENE_ONLY_EXTRA)
    if payload.get("msa_source") == "coregene":
        drop = {"clean_mode", "filtered"} & _coregene_ignored(str(payload.get("coregene")))
    if not drop:
        return payload
    return {**payload, "extra": {k: v for k, v in extra.items() if k not in drop}}


def _read_msa_stamp(artifact: Path, key: str) -> dict | None:
    """The stamp beside ``artifact`` when it still describes that exact file.

    The content digest is checked as well as the key: another stage (or a hand
    edit) may have rewritten the artifact since, and a stale reuse would build
    a tree from an alignment nobody asked for.
    """
    stamp_path = artifact.parent / MSA_STAMP
    if not stamp_path.is_file() or not artifact.is_file():
        return None
    try:
        stamp = json.loads(stamp_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if stamp.get("v") != _MSA_STAMP_VERSION or stamp.get("key") != key:
        return None
    if stamp.get("artifact_digest") != file_digest(artifact):
        return None
    return stamp


def _write_msa_stamp(artifact: Path, key: str, versions: dict[str, str]) -> None:
    stamp = {
        "v": _MSA_STAMP_VERSION,
        "key": key,
        "artifact": artifact.name,
        "artifact_digest": file_digest(artifact),
        "versions": versions,
        # The record names, so a reuse checks them without reading the MSA.
        "records": sorted(set(fasta_record_names(artifact))),
    }
    _save_msa_stamp(artifact, stamp)


def _save_msa_stamp(artifact: Path, stamp: dict) -> None:
    with atomic_path(artifact.parent / MSA_STAMP) as tmp:
        tmp.write_text(json.dumps(stamp, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fasta_record_names(path: Path) -> list[str]:
    """The record names (first header token) of a FASTA file, in file order.

    Searches the memory-mapped file for line starts with '>', so sequence
    lines are skipped without being split into lines or decoded.
    """
    with open(path, "rb") as fh:
        if os.fstat(fh.fileno()).st_size == 0:
            return []
        with mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            names: list[str] = []
            if mm[:1] == b">":
                pos = 0
            else:
                found = mm.find(b"\n>")
                pos = found + 1 if found >= 0 else -1
            while pos >= 0:
                end = mm.find(b"\n", pos)
                if end < 0:
                    end = len(mm)
                fields = mm[pos + 1 : end].decode("utf-8", "replace").split()
                names.append(fields[0] if fields else "")
                found = mm.find(b"\n>", end)
                pos = found + 1 if found >= 0 else -1
            return names


def _annotations_digest(
    genomes: Sequence[Path], outgroup_file: Path | None, dirs: PhyloDirs, params: PhyloParams
) -> str | None:
    """The index digest of the annotations the core-gene pass will read, so a
    re-annotated genome rebuilds the alignment; None for the other sources."""
    if params.msa_source != "coregene" or dirs.annotations_dir is None:
        return None
    from .coregene import annotated_set

    inputs = list(genomes) + ([outgroup_file] if outgroup_file is not None else [])
    return annotated_set(dirs.annotations_dir, [record_name(g) for g in inputs]).digest


def _msa_with_reuse(
    genomes: list[Path],
    outgroup_file: Path | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
    *,
    reuse: bool,
    sketch_provider: SketchProvider | None = None,
) -> tuple[Path, dict[str, str]]:
    """Build the MSA, or reuse the one a previous run left with the same inputs."""
    choice = _reference_for_source(genomes, outgroup_file, dirs, params, logger, sketch_provider)
    reference = choice.genome if choice is not None else None
    annotations = _annotations_digest(genomes, outgroup_file, dirs, params)
    key = _msa_key(genomes, outgroup_file, params, reference, annotations)
    artifact = _msa_artifact(dirs, params)
    if reuse:
        stamp = _read_msa_stamp(artifact, key)
        if stamp is not None and not _msa_records_named(artifact, stamp, genomes, outgroup_file):
            # An aligner MSA from before every adapter used record_name names a
            # gzipped genome x.fasta.gz 'x.fasta'. Rebuilding only such an MSA
            # spares the workdirs that a stamp version change would realign.
            logger.info(
                "The alignment at %s names records other than the genomes' record "
                "names; it is rebuilt.",
                artifact,
            )
            stamp = None
        if stamp is not None:
            logger.info(
                "Reusing the alignment at %s: same source, inputs and settings "
                "as the run that built it (--force rebuilds it).",
                artifact,
            )
            # The same genome under a different rule (the medoid turned out to
            # be the first genome, say) reuses the alignment; the note and the
            # record then follow this run's rule.
            if choice is not None:
                note_reference(artifact.parent, choice)
            return artifact, dict(stamp.get("versions") or {})
    msa, versions = _build_msa(genomes, outgroup_file, dirs, params, logger, reference)
    # The note names the reference the alignment was built against; written
    # after the build, so a failed build cannot leave a note for an older
    # alignment, and cleared when no reference applies.
    if choice is not None:
        note_reference(msa.parent, choice)
    else:
        clear_reference_note(msa.parent)
    if msa.resolve() == artifact.resolve():
        _write_msa_stamp(artifact, key, versions)
    return msa, versions


def _msa_records_named(
    artifact: Path, stamp: dict, genomes: Sequence[Path], outgroup_file: Path | None
) -> bool:
    """True when every record of ``artifact`` is the record name of an input.

    The names come from the stamp, which _read_msa_stamp has matched to the
    artifact's content digest. A stamp written before the stamp recorded them
    is completed once from the header lines, so the next reuse reads nothing.
    A subset, not equality: a typer may leave out a genome it could not type,
    and the tree leaf check reports that.
    """
    records = stamp.get("records")
    if not isinstance(records, list):
        records = sorted(set(fasta_record_names(artifact)))
        _save_msa_stamp(artifact, {**stamp, "records": records})
    inputs = [*genomes, outgroup_file] if outgroup_file is not None else list(genomes)
    return set(records) <= {record_name(g) for g in inputs}


def build_msa(
    genomes: list[Path],
    outgroup_file: Path | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
    *,
    reuse: bool = True,
    sketch_provider: SketchProvider | None = None,
) -> tuple[Path, dict[str, str]]:
    """Build only the alignment a tree builder would consume.

    The alignment half of :func:`build_tree`, for callers that build the tree
    separately: the Nextflow layer runs the two as their own tasks so that
    changing the tree builder does not repeat the alignment.
    """
    n_leaves = len(genomes) + (0 if outgroup_file is None else 1)
    builder = treebuilder_registry.create(
        _resolve_treebuilder(params, n_leaves, logger, alignment_given="--msa-only")
    )
    if builder.input_kind == InputKind.GENOMES:
        raise UserInputError(
            f"Tree builder '{builder.capabilities.name}' builds from genomes, not from an "
            "alignment, so there is no alignment to build separately. Choose a builder "
            f"that consumes one ({', '.join(_msa_builder_names())}) or build the tree in "
            "one step."
        )
    return _msa_with_reuse(
        genomes, outgroup_file, dirs, params, logger, reuse=reuse, sketch_provider=sketch_provider
    )


def _aligner_limit(name: str) -> int | None:
    """The declared genome limit of aligner ``name`` (None when unbounded or unknown)."""
    try:
        return aligner_registry.get(name).capabilities.recommended_max_genomes
    except (PluginError, KeyError):
        return None


def _msa_builder_names() -> list[str]:
    """Registered tree builders that consume an alignment."""
    names = []
    for name in treebuilder_registry.names():
        try:
            cls = treebuilder_registry.get(name)
        except (PluginError, KeyError):
            continue
        if cls.input_kind != InputKind.GENOMES:
            names.append(name)
    return names


def _wants_alignment(params: PhyloParams) -> str | None:
    """The request that only an alignment-consuming builder can honour, or None."""
    if str(params.extra.get("mask", "none")) not in ("none", ""):
        return f"--mask {params.extra['mask']}"
    if params.msa_source == "snptype":
        return "--msa-source snptype"
    if params.msa_source == "coregene":
        return "--msa-source coregene"
    return None


def _resolve_treebuilder(
    params: PhyloParams,
    n_genomes: int,
    logger: logging.Logger,
    *,
    alignment_given: str | None = None,
) -> str:
    """The tree builder to use, resolving 'auto' and warning past its scale.

    ``n_genomes`` counts every leaf the tree will have, the outgroup included.
    With 'auto', a request that only an alignment can honour (a SNP source, a
    masker, or ``alignment_given``, the flag under which the caller builds or
    supplies the alignment itself) restricts the choice to builders that
    consume one; otherwise a large set would be handed to an alignment-free
    builder and the request silently dropped.
    """
    treebuilder = params.treebuilder
    if treebuilder == "auto":
        wants = alignment_given or _wants_alignment(params)
        allowed: Collection[str] | None = _msa_builder_names() if wants else None
        note = f" among builders that consume an alignment ({wants})" if wants else ""
        if not wants:
            # Nothing asked for an alignment, so a builder that needs one would
            # make the aligner run first; past the aligner's own scale that is
            # not a choice `auto` should make on the user's behalf.
            aligner_limit = _aligner_limit(params.aligner)
            if aligner_limit is not None and n_genomes > aligner_limit:
                msa_names = set(_msa_builder_names())
                allowed = [n for n in treebuilder_registry.names() if n not in msa_names]
                note = (
                    f" among alignment-free builders (aligner '{params.aligner}' is tuned "
                    f"for <={aligner_limit} genomes)"
                )
        treebuilder = auto_select(treebuilder_registry, n_genomes, allowed) or "iqtree"
        logger.info(
            "Auto-selected tree builder '%s' for %d genomes%s", treebuilder, n_genomes, note
        )
        return treebuilder
    warn = scale_warning(treebuilder_registry, treebuilder, n_genomes)
    if warn:
        limit, alts = warn
        logger.warning(
            "Tree builder '%s' is tuned for <=%d genomes but you have %d; consider: %s",
            treebuilder,
            limit,
            n_genomes,
            ", ".join(alts) or "none",
        )
    return treebuilder


def build_tree(
    genomes: list[Path],
    outgroup_file: Path | None,
    outgroup_leaf: str | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
    *,
    reuse_msa: bool = True,
    msa: Path | None = None,
    sketch_provider: SketchProvider | None = None,
    reference_provider: SketchProvider | None = None,
) -> PhyloOutcome:
    """Build a phylogeny from explicit inputs into ``dirs`` (stateless; no config).

    Selects (or auto-picks) the tree builder, derives the MSA from an aligner or
    a SNP typer when needed, writes ``tree/tree.nwk`` and roots it as
    ``params.root`` asks (by default on the outgroup when one is staged).
    The reference (for aligner/SNP sources) is resolved by basename against the
    genome set, so the core needs no working directory. ``sketch_provider``
    supplies genome sketches to a builder that compares them.
    """
    if not genomes:
        raise WorkdirError("No genomes found for phylo. Run the genome (and derep) stages first.")
    # Decided before any tool runs, so an impossible request costs nothing.
    root = resolve_root(params.root, outgroup_leaf, params.no_outgroup)

    n_leaves = len(genomes) + (0 if outgroup_file is None else 1)
    treebuilder = _resolve_treebuilder(
        params, n_leaves, logger, alignment_given=None if msa is None else "--msa"
    )
    builder = treebuilder_registry.create(treebuilder)
    if (
        builder.input_kind != InputKind.GENOMES
        and params.bootstrap == 0
        and not builder.supports_without_bootstrap
    ):
        logger.info(
            "No bootstrap replicates requested (-B not set): the %s tree will carry no "
            "support values, so tree2tax --collapse-support has nothing to act on",
            treebuilder,
        )
    versions = builder.preflight()

    # One warning for the whole stage: the tree builder and the MSA source share
    # a single extras dict, so a key is unread only when neither declares it.
    participating = [builder.capabilities]
    if builder.input_kind != InputKind.GENOMES:
        source_caps = _msa_source_capabilities(params)
        if source_caps is not None:
            participating.append(source_caps)
        mask = str(params.extra.get("mask", "none"))
        if params.msa_source == "snptype" and mask not in ("none", ""):
            from ..maskers.base import registry as masker_registry

            participating.append(masker_registry.get(mask).capabilities)
    _warn_stage_extras(participating, _adapter_extra(params.extra), logger)

    tree_params = TreeParams(
        threads=params.threads,
        outgroup=None if params.no_outgroup else outgroup_leaf,
        bootstrap=params.bootstrap,
        extra={**builder.capabilities.default_params, **_adapter_extra(params.extra)},
    )
    warn_ignored_params(builder.capabilities, tree_params, logger, family="Tree builder")
    dirs.tree_dir.mkdir(parents=True, exist_ok=True)

    if builder.input_kind == InputKind.GENOMES:
        unused = []
        if params.msa_source == "snptype":
            unused.append(f"--msa-source snptype (--snptyper {params.snptyper})")
        if params.msa_source == "coregene":
            unused.append(f"--msa-source coregene (--coregene {params.coregene})")
        if str(params.extra.get("mask", "none")) not in ("none", ""):
            unused.append(f"--mask {params.extra['mask']}")
        if _core_fraction(params.extra) is not None:
            unused.append(f"--core-fraction {params.extra['core_fraction']}")
        if unused:
            logger.warning(
                "Tree builder '%s' builds from the genomes without an alignment; %s has no effect.",
                treebuilder,
                " and ".join(unused),
            )
        inputs = list(genomes)
        if outgroup_file is not None:
            inputs.append(outgroup_file)
        logger.info(
            "Building tree with %s (alignment-free) over %d genomes",
            treebuilder,
            len(inputs),
        )
        _clear_previous_builder_files(dirs.tree_dir)
        # Both genome-input builders read gzip today; the rule is the same as
        # for the aligners, so a new builder that does not is given copies.
        inputs_dir = dirs.scratch_dir / _STAGED_INPUTS
        remove_tree(inputs_dir)  # copies a killed run left behind
        sketches = adapter_sketches(
            builder, tree_params.extra, inputs, sketch_provider, logger, "phylo"
        )
        staged = stage_plain_inputs(inputs, builder.capabilities, inputs_dir, logger)
        if sketches is not None:
            tree_params.sketches = {staged[p]: sk for p, sk in sketches.items()}
        try:
            tree = builder.build([staged[p] for p in inputs], dirs.tree_dir, tree_params, logger)
        finally:
            remove_tree(inputs_dir)
    else:
        given_msa = msa
        if msa is not None:
            logger.info("Building tree from the alignment given on the command line: %s", msa)
            source_versions: dict[str, str] = {}
        else:
            msa, source_versions = _msa_with_reuse(
                genomes,
                outgroup_file,
                dirs,
                params,
                logger,
                reuse=reuse_msa,
                sketch_provider=reference_provider
                if reference_provider is not None
                else sketch_provider,
            )
        versions = {**versions, **source_versions}
        _warn_low_diversity(msa, logger)
        tree_params.constant_sites = read_constant_sites(msa)
        if (
            tree_params.constant_sites is None
            and msa is not given_msa
            and params.msa_source == "snptype"
        ):
            logger.warning(
                "No constant-site counts for %s (SNP typer '%s' writes no whole-genome "
                "alignment, or alignment_stats.json is missing or describes another "
                "file): the tree's branch lengths carry no ascertainment-bias correction "
                "and are inflated",
                msa,
                params.snptyper,
            )
        elif tree_params.constant_sites is None and params.msa_source == "coregene":
            logger.info(
                "Core-gene alignment: invariant sites kept, no ascertainment-bias correction"
            )
        elif tree_params.constant_sites is not None:
            logger.info(
                "Constant sites given to the tree builder: A %d, C %d, G %d, T %d",
                *tree_params.constant_sites,
            )
        logger.info("Building tree with %s from MSA %s", treebuilder, msa)
        _clear_previous_builder_files(dirs.tree_dir)
        tree = builder.build(msa, dirs.tree_dir, tree_params, logger)

    final = dirs.tree_dir / TREE_NWK
    if tree.resolve() != final.resolve():
        # Copy through a temp sibling so a previous tree.nwk is never left
        # truncated or half-written if the copy fails.
        with atomic_path(final) as tmp:
            shutil.copy2(tree, tmp)
    logger.info("Phylogenetic tree written to %s", final)
    expected = [*genomes, outgroup_file] if outgroup_file is not None else list(genomes)
    leaf_names = [record_name(g) for g in expected]
    check_tree_leaves(final, leaf_names, treebuilder)
    restore_leaf_names(final, leaf_names, logger)
    apply_root(final, root, outgroup_leaf, treebuilder, logger)
    return PhyloOutcome(
        tree=final,
        treebuilder=treebuilder,
        versions=versions,
        outgroup_leaf=outgroup_leaf,
        root=root,
    )


_LEAF_UNSAFE = re.compile(r"[^A-Za-z0-9_-]")
_LEAF_FILE_SUFFIX = re.compile(r"(\.(fasta|fas|fna|fa|ffn))?(\.gz)?(\.ref)?$", re.IGNORECASE)
_MAX_LEAVES_LISTED = 10


def _leaf_key(label: str) -> str:
    """Compare leaf names the way tools may rewrite them.

    ParSNP names sequences by file name and marks the reference with '.ref',
    so a FASTA extension and that marker are dropped. Tools also replace
    characters they do not accept in a name (cactus turns '.' into '_',
    IQ-TREE rewrites others), so every character outside letters, digits,
    '_' and '-' compares as '_'.
    """
    return _LEAF_UNSAFE.sub("_", _LEAF_FILE_SUFFIX.sub("", label, count=1))


def check_tree_leaves(tree: Path, expected: Sequence[str], treebuilder: str) -> None:
    """Raise WorkdirError when the leaves of ``tree`` differ from ``expected``.

    A tree builder can drop a genome it considers degenerate and still exit 0
    (mashtree does). The tree stays on disk for inspection; the caller records
    nothing, so the stage is not marked completed.
    """
    import dendropy
    from dendropy.utility.error import DataParseError

    try:
        parsed = dendropy.Tree.get(path=str(tree), schema="newick", preserve_underscores=True)
    except DataParseError as exc:
        raise WorkdirError(
            f"Tree builder '{treebuilder}' wrote an unreadable tree {tree}: {exc}"
        ) from exc
    leaves = {
        _leaf_key(node.taxon.label): node.taxon.label
        for node in parsed.leaf_node_iter()
        if node.taxon is not None and node.taxon.label
    }
    wanted = {_leaf_key(name): name for name in expected}
    missing = sorted(wanted[k] for k in wanted.keys() - leaves.keys())
    extra = sorted(leaves[k] for k in leaves.keys() - wanted.keys())
    if not missing and not extra:
        return
    parts = []
    if missing:
        parts.append(f"missing {len(missing)} genome(s): {_listed(missing)}")
    if extra:
        parts.append(f"{len(extra)} unexpected leaf/leaves: {_listed(extra)}")
    raise WorkdirError(
        f"The tree from '{treebuilder}' does not match its input genomes ("
        + "; ".join(parts)
        + f"). The tree is kept at {tree} for inspection; the stage is not recorded "
        "as completed."
    )


def restore_leaf_names(tree: Path, expected: Sequence[str], logger: logging.Logger) -> int:
    """Rename leaves a tool rewrote back to the input genome names.

    check_tree_leaves accepts a leaf that differs from its genome only in the
    way tools rewrite names (a FASTA extension, '.ref', characters replaced by
    '_'). tree2tax and genomes_map.tsv read the leaves as genome names, so the
    tree is rewritten with the input names. Only the renamed leaf labels change
    in the text; every other byte (rooting tags, comments, quoting, lengths,
    supports) stays as the tool wrote it. A name that two inputs would share is
    left as written. Returns the number of leaves renamed.
    """
    keys: dict[str, list[str]] = {}
    for name in expected:
        keys.setdefault(_leaf_key(name), []).append(name)

    def target(label: str) -> str | None:
        names = keys.get(_leaf_key(label), [])
        if len(names) == 1 and names[0] != label:
            return names[0]
        return None

    text = tree.read_text(encoding="utf-8")
    new_text, renamed = _rename_newick_leaves(text, target)
    if renamed:
        with atomic_path(tree) as tmp:
            tmp.write_text(new_text, encoding="utf-8")
        logger.info("Renamed %d tree leaf/leaves back to the input genome names", renamed)
    return renamed


_NEWICK_LABEL_END = set("():,;[") | set(" \t\r\n")
_NEWICK_NEEDS_QUOTES = re.compile(r"[\s():,;\[\]']")


def _rename_newick_leaves(text: str, target) -> tuple[str, int]:  # noqa: ANN001
    """Replace leaf labels in Newick ``text`` for which ``target(label)`` gives a name.

    A leaf label is the label that follows '(' or ','; labels after ')' are
    internal (supports or names) and are never touched. Quoted labels are read
    with '' as an escaped quote; [comments] are skipped. A new name is written
    bare when it holds no character that Newick reserves, and quoted otherwise.
    """
    out: list[str] = []
    renamed = 0
    i, n = 0, len(text)
    expect_leaf = False
    while i < n:
        char = text[i]
        if char == "[":
            close = text.find("]", i)
            close = n - 1 if close < 0 else close
            out.append(text[i : close + 1])
            i = close + 1
            continue
        if char in "(,":
            out.append(char)
            expect_leaf = True
            i += 1
            continue
        if char.isspace():
            out.append(char)
            i += 1
            continue
        if expect_leaf and char not in "():;":
            # Read one label, quoted or bare.
            if char == "'":
                j = i + 1
                parts: list[str] = []
                while j < n:
                    if text[j] == "'":
                        if j + 1 < n and text[j + 1] == "'":
                            parts.append("'")
                            j += 2
                            continue
                        break
                    parts.append(text[j])
                    j += 1
                label, raw_end = "".join(parts), j + 1
            else:
                j = i
                while j < n and text[j] not in _NEWICK_LABEL_END:
                    j += 1
                label, raw_end = text[i:j], j
            new = target(label)
            if new is None:
                out.append(text[i:raw_end])
            else:
                if _NEWICK_NEEDS_QUOTES.search(new):
                    new = "'" + new.replace("'", "''") + "'"
                out.append(new)
                renamed += 1
            i = raw_end
            expect_leaf = False
            continue
        if char == "'":
            # A quoted internal label: copy it whole, so a ',' or '(' inside
            # it is not read as structure.
            j = i + 1
            while j < n:
                if text[j] == "'":
                    if j + 1 < n and text[j + 1] == "'":
                        j += 2
                        continue
                    break
                j += 1
            out.append(text[i : j + 1])
            i = j + 1
            expect_leaf = False
            continue
        expect_leaf = False
        out.append(char)
        i += 1
    return "".join(out), renamed


def _listed(names: list[str]) -> str:
    shown = ", ".join(names[:_MAX_LEAVES_LISTED])
    if len(names) > _MAX_LEAVES_LISTED:
        shown += f" (+{len(names) - _MAX_LEAVES_LISTED} more)"
    return shown


def _clear_previous_builder_files(tree_dir: Path) -> None:
    """Remove an earlier build's side files from ``tree_dir``.

    ``tree/`` holds the current tree builder's own files; a matrix or a set of
    bootstrap trees left by another builder would describe a different run.
    ``tree.nwk`` stays until the new tree replaces it atomically, and
    ``msa/`` stays: it is the reuse cache of the SNP typing pass, which may
    date from an earlier run with another MSA source. An entry may
    vanish after listing (macOS drops a file's ``._`` AppleDouble sibling on
    non-HFS volumes when the file is removed), so missing entries are skipped.
    """
    for entry in list(tree_dir.iterdir()):
        if entry.name in (TREE_NWK, PHYLO_MSA_DIR):
            continue
        if entry.is_dir() and not entry.is_symlink():
            remove_tree(entry)
        else:
            entry.unlink(missing_ok=True)


@dataclass
class PhyloBuildParams:
    """Inputs for the stateless phylo step (explicit paths, no workdir)."""

    genomes_dir: Path
    out_dir: Path
    outgroup_dir: Path | None = None
    outgroup_accession: Path | None = None
    phylo: PhyloParams = field(default_factory=PhyloParams)
    versions_out: Path | None = None
    # Split the step in two: build the alignment and stop (``msa_only``), or
    # build the tree from an alignment a previous call produced (``msa``).
    # The Nextflow layer uses this to keep the alignment when only the tree
    # builder changes; nothing else changes between the two halves.
    msa_only: bool = False
    msa: Path | None = None
    # A sketches/ directory (Nextflow SKETCH) for a builder that compares
    # sourmash sketches; matched to the genomes by record name.
    sketches_dir: Path | None = None
    # --msa-source coregene: an annotations/ directory with an index, or None
    # to annotate the inputs here with ``annotator`` (bakta needs ``bakta_db``).
    annotations_dir: Path | None = None
    annotator: str = "bakta"
    bakta_db: str | None = None


MIN_TREE_GENOMES = 3


def _require_tree_size(genomes: Sequence[Path], hint: str) -> None:
    """Refuse an ingroup too small for a tree before any tool runs.

    The outgroup is not counted. Tree builders fail on fewer than three
    leaves with errors that do not name the cause.
    """
    if len(genomes) < MIN_TREE_GENOMES:
        raise WorkdirError(
            f"A tree needs at least {MIN_TREE_GENOMES} genomes; {len(genomes)} found {hint}."
        )


def phylo_build(params: PhyloBuildParams, logger: logging.Logger) -> Path:
    """Build a phylogeny from an explicit genomes directory (data-channel step).

    With ``msa_only`` it stops once the alignment is built and returns it,
    copied to ``<out_dir>/msa.fasta``. With ``msa`` it skips alignment and
    builds the tree from that file.
    """
    genomes = list_fasta(params.genomes_dir)
    if not genomes:
        raise WorkdirError(f"No genome FASTA files found in {params.genomes_dir}.")
    if params.sketches_dir is not None and not params.sketches_dir.is_dir():
        raise WorkdirError(f"phylo-build: --sketches-dir {params.sketches_dir} is not a directory.")
    _require_tree_size(genomes, f"in {params.genomes_dir}")

    outgroup_file: Path | None = None
    outgroup_leaf: str | None = None
    if (
        not params.phylo.no_outgroup
        and params.outgroup_dir is not None
        and params.outgroup_accession is not None
    ):
        outgroup_file, outgroup_leaf = resolve_outgroup_files(
            params.outgroup_dir, params.outgroup_accession, logger
        )

    dirs = PhyloDirs(
        tree_dir=params.out_dir / "tree",
        align_dir=params.out_dir / "align",
        msa_dir=params.out_dir / "tree" / PHYLO_MSA_DIR,
        scratch_dir=params.out_dir / "scratch",
        annotations_dir=_step_annotations(params, genomes, outgroup_file, logger),
    )
    if params.msa_only:
        msa, versions = build_msa(
            genomes,
            outgroup_file,
            dirs,
            params.phylo,
            logger,
            sketch_provider=(
                None
                if params.sketches_dir is None
                else directory_provider(params.sketches_dir, logger, "phylo-build")
            ),
        )
        published = params.out_dir / MSA_FASTA
        if msa.resolve() != published.resolve():
            with atomic_path(published) as tmp:
                shutil.copy2(msa, tmp)
            # The counts travel with the copy, under the copy's name and digest.
            (published.parent / ALIGNMENT_STATS).unlink(missing_ok=True)
            copy_alignment_stats(msa, published)
        _write_versions(params.versions_out, versions)
        return published

    outcome = build_tree(
        genomes,
        outgroup_file,
        outgroup_leaf,
        dirs,
        params.phylo,
        logger,
        msa=params.msa,
        sketch_provider=(
            None
            if params.sketches_dir is None
            else directory_provider(params.sketches_dir, logger, "phylo-build")
        ),
    )
    _write_versions(params.versions_out, outcome.versions)
    return outcome.tree


def _annotations_record(
    ctx: WorkdirContext, params: PhyloParams, names: Sequence[str]
) -> dict | None:
    """The annotator and the index digest of the tree's own genomes, so the
    record describes what the alignment was built from (the same rows the
    reuse stamp hashes), not whatever else the cache holds."""
    if params.msa_source != "coregene":
        return None
    from .coregene import annotated_set

    try:
        annotated = annotated_set(ctx.annotations_dir, list(names))
    except Exception:  # noqa: BLE001 - the record is an annotation, never a failure
        return None
    return {"annotator": annotated.annotator, "index_digest": annotated.digest}


def _step_annotations(
    params: PhyloBuildParams,
    genomes: Sequence[Path],
    outgroup_file: Path | None,
    logger: logging.Logger,
) -> Path | None:
    """The annotations the step's core-gene pass reads: the given directory,
    or one it builds under out_dir/annotations by annotating every input."""
    if params.phylo.msa_source != "coregene" or params.msa is not None:
        return None
    if params.annotations_dir is not None:
        if not params.annotations_dir.is_dir():
            raise WorkdirError(
                f"phylo-build: --annotations-dir {params.annotations_dir} is not a directory."
            )
        return params.annotations_dir
    from ..annotators.base import registry as annotator_registry
    from ..core.annotations import gff_path, read_index
    from ..core.executors import parallel_each
    from ..core.process import remove_tree
    from .annotate import (
        AnnotateBuildParams,
        _database,
        annotate_build,
        annotate_jobs,
        index_from_directories,
    )

    out = params.out_dir / "annotations"
    inputs = list(genomes) + ([outgroup_file] if outgroup_file is not None else [])
    # Reuse only what this annotator wrote: a rerun with another --annotator
    # must not leave prokka output indexed as bakta's.
    existing = read_index(out)
    todo = [
        g
        for g in inputs
        if not gff_path(out, record_name(g)).is_file()
        or existing.get(record_name(g), None) is None
        or existing[record_name(g)].annotator != params.annotator
    ]
    jobs, per_genome = annotate_jobs(params.phylo.threads, 0, len(todo))
    logger.info(
        "Annotating %d of %d genomes with %s under %s (the rest are already there), "
        "%d at a time with %d thread(s) each",
        len(todo),
        len(inputs),
        params.annotator,
        out,
        jobs,
        per_genome,
    )
    versions: dict[str, str] | None = None
    if todo and params.annotator in annotator_registry.names():
        # One preflight for the whole set: under a container backend each
        # check is a container run.
        annotator = annotator_registry.create(params.annotator)
        db = _database(params.annotator, annotator, params.bakta_db)
        versions = annotator.preflight()
        if db is not None:
            versions[f"{params.annotator}-db"] = annotator.database_version(db)
    try:
        parallel_each(
            lambda genome: annotate_build(
                AnnotateBuildParams(
                    genome=genome,
                    out_dir=out,
                    tool=params.annotator,
                    bakta_db=params.bakta_db,
                    threads=per_genome,
                    extra={},
                ),
                logger,
                versions=versions,
                shared_out_dir=True,
            ),
            todo,
            jobs,
            lambda _genome, _result: None,
            logger=logger,
            describe=record_name,
        )
    except BaseException:
        remove_tree(out / ".scratch")
        # Index what finished, so a rerun reuses it; the build error is the
        # one to report.
        if out.is_dir():
            try:
                index_from_directories(out, params.annotator, logger)
            except Exception as exc:  # noqa: BLE001 - the build error is the one to report
                logger.warning("Could not index the finished annotations: %s", exc)
        raise
    remove_tree(out / ".scratch")
    index_from_directories(out, params.annotator, logger)
    return out


def _write_versions(versions_out: Path | None, versions: dict[str, str]) -> None:
    if versions_out is None:
        return
    from ..core.versions import write_versions_fragment

    write_versions_fragment(versions_out, versions)


def run(ctx: WorkdirContext, params: PhyloParams) -> Path:
    logger = ctx.logger
    if params.all_genomes:
        check_genome_completeness(
            ctx.genomes_dir,
            ctx.workdir,
            logger=logger,
            allow_incomplete=params.allow_incomplete,
        )
    else:
        check_representatives_consistency(
            ctx.representatives_dir,
            ctx.derep_dir / CLUSTERS_TSV,
            logger=logger,
            allow_incomplete=params.allow_incomplete,
        )
    genomes = _genome_set(ctx, params.all_genomes)
    if not genomes:
        raise WorkdirError("No genomes found for phylo. Run the genome (and derep) stages first.")
    _require_tree_size(
        genomes,
        "in the genome set"
        if params.all_genomes
        else "after dereplication (use --all-genomes, or dereplicate with a higher "
        "--secondary-ani, which keeps more representatives)",
    )

    outgroup_file, outgroup_leaf = _resolve_outgroup(ctx, params.no_outgroup, logger)

    dirs = PhyloDirs(
        tree_dir=ctx.tree_dir,
        align_dir=ctx.align_dir,
        msa_dir=ctx.phylo_msa_dir,
        scratch_dir=ctx.scratch_dir,
        annotations_dir=ctx.annotations_dir,
    )
    outcome = build_tree(
        genomes,
        outgroup_file,
        outgroup_leaf,
        dirs,
        params,
        logger,
        # --force means recompute this stage, cached alignment included.
        reuse_msa=not ctx.force,
        sketch_provider=workdir_provider(ctx, logger, params.threads, "phylo"),
        reference_provider=workdir_provider(ctx, logger, params.threads, "reference", write=False),
    )

    is_msa = treebuilder_registry.create(outcome.treebuilder).input_kind == InputKind.MSA_FASTA
    tree_inputs = [record_name(g) for g in genomes]
    if outgroup_file is not None:
        tree_inputs.append(record_name(outgroup_file))
    annotations_record = _annotations_record(ctx, params, tree_inputs) if is_msa else None
    reference_choice = (
        read_reference_choice(
            ctx.phylo_msa_dir if params.msa_source == "snptype" else ctx.align_dir
        )
        if is_msa
        else None
    )
    ctx.config.record_stage(
        "phylo",
        tool=outcome.treebuilder,
        params={
            "requested_treebuilder": params.treebuilder,
            # Alignment-free builders never run an MSA source; recording one
            # would put --aligner into a fingerprint it does not influence.
            "msa_source": params.msa_source if is_msa else None,
            "aligner": params.aligner if is_msa and params.msa_source == "aligner" else None,
            "snptyper": params.snptyper if is_msa and params.msa_source == "snptype" else None,
            "coregene": params.coregene if is_msa and params.msa_source == "coregene" else None,
            # The core-gene pass: its settings and the annotations it read.
            "core_threshold": (
                _core_threshold(params.extra)
                if is_msa and params.msa_source == "coregene"
                else None
            ),
            "clean_mode": (
                _coregene_setting(params, "clean_mode")
                if is_msa and params.msa_source == "coregene"
                else None
            ),
            "gene_selection": (
                asdict(_gene_selection(params.extra))
                if is_msa and params.msa_source == "coregene"
                else None
            ),
            "annotations": annotations_record,
            "all_genomes": params.all_genomes,
            "bootstrap": params.bootstrap,
            "outgroup": None if params.no_outgroup else outcome.outgroup_leaf,
            # The rooting applied to tree.nwk; tree2tax keeps a midpoint root.
            "root": outcome.root,
            # The genome the alignment was projected onto and the rule that
            # chose it (medoid, first, first-fallback, name), from the note the
            # MSA source leaves beside it (None for ska2 or an alignment-free builder).
            "reference_genome": reference_choice.genome.name if reference_choice else None,
            "reference_rule": reference_choice.rule if reference_choice else None,
            # What produced the alignment: the reference a mapping typer or
            # aligner used, the masker if any, and the adapter tuning. Without
            # them the record could not tell a masked tree from an unmasked one.
            "reference": params.reference if is_msa else None,
            "mask": (params.extra.get("mask") or None) if is_msa else None,
            "extra": {k: str(v) for k, v in sorted(_adapter_extra(params.extra).items())},
        },
        tool_versions=outcome.versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return outcome.tree


_SILENT = logging.getLogger(__name__ + ".precheck")
_SILENT.addHandler(logging.NullHandler())
_SILENT.propagate = False


def precheck(ctx: WorkdirContext, params: PhyloParams) -> None:
    """Refusals checked before the harness marks the stage record incomplete.

    With ``--msa-source coregene`` and a builder that reads an alignment, a
    genome of the set (or the outgroup) without an annotation, or a set that
    mixes annotators, is refused here, before any tool runs, so the record of
    the last finished tree stays complete. The core-gene pass repeats the
    check for callers that bypass the harness.
    """
    if params.msa_source != "coregene":
        return
    if params.treebuilder != "auto":
        try:
            kind = treebuilder_registry.get(params.treebuilder).input_kind
        except Exception:  # noqa: BLE001 - an unknown builder is refused by the stage
            return
        if kind != InputKind.MSA_FASTA:
            return
    from .coregene import annotated_set

    names = [record_name(g) for g in _genome_set(ctx, params.all_genomes)]
    if len(names) < MIN_TREE_GENOMES:
        return  # the stage body names the missing genomes or the size limit first
    # The stage body logs the outgroup choice; the precheck stays silent.
    outgroup_file, _ = _resolve_outgroup(ctx, params.no_outgroup, _SILENT)
    if outgroup_file is not None:
        names.append(record_name(outgroup_file))
    annotated_set(ctx.annotations_dir, names)


def _genome_set(ctx: WorkdirContext, all_genomes: bool) -> list[Path]:
    source = ctx.genomes_dir if all_genomes else ctx.representatives_dir
    return list_fasta(source)


def _resolve_outgroup(
    ctx: WorkdirContext, no_outgroup: bool, logger
) -> tuple[Path | None, str | None]:
    if no_outgroup:
        return None, None
    acc_file = ctx.workdir / "outgroup_accession.txt"
    if not acc_file.exists() or not ctx.outgroup_dir.exists():
        logger.warning("No outgroup found; proceeding without one")
        return None, None
    return resolve_outgroup_files(ctx.outgroup_dir, acc_file, logger)


def resolve_outgroup_files(
    outgroup_dir: Path, accession_file: Path, logger: logging.Logger
) -> tuple[Path | None, str | None]:
    """Resolve the outgroup genome file and leaf name from explicit paths."""
    if not accession_file.exists() or not outgroup_dir.exists():
        logger.warning("No outgroup found; proceeding without one")
        return None, None
    accession = accession_file.read_text(encoding="utf-8").strip()
    if not accession:
        logger.warning("No outgroup accession recorded; proceeding without one")
        return None, None
    # Genome FASTA files only: a leftover GCF_x.fasta.tmp or a .fai index
    # must not win the substring match below.
    candidates = list_fasta(outgroup_dir)
    # Exact match on the parsed accession first; a stale file whose name merely
    # contains the accession must not shadow the intended outgroup.
    for f in candidates:
        if accession_from_filename(f.name) == accession:
            logger.info("Using %s as outgroup", f.name)
            return f, record_name(f)
    for f in candidates:
        if accession in f.name:
            logger.info("Using %s as outgroup (substring match)", f.name)
            return f, record_name(f)
    logger.warning("Outgroup accession %s not found in %s", accession, outgroup_dir)
    return None, None


def snp_registry_requires_reference(name: str) -> bool:
    """Whether the SNP typer ``name`` maps against a reference."""
    from ..snptypers.base import registry as snp_registry

    return bool(getattr(snp_registry.get(name), "requires_reference", True))


def _reference_for_source(
    genomes: Sequence[Path],
    outgroup_file: Path | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
    sketch_provider: SketchProvider | None,
) -> ReferenceChoice | None:
    """The reference the MSA source projects onto, or None when it needs none
    (a reference-free typer; a name given there is still checked). Resolved
    before the reuse check, so the cached
    alignment is reused only for the same genome."""
    if params.msa_source == "coregene":
        # Gene alignments project onto no genome; a name given is still checked.
        if params.reference and params.reference not in REFERENCE_RULES_ALL:
            resolve_reference(params.reference, genomes, outgroup_file, logger)
            logger.info(
                "The core-gene source uses no reference; --reference %s is not used",
                params.reference,
            )
        return None
    if params.msa_source == "snptype" and not snp_registry_requires_reference(params.snptyper):
        if params.reference and params.reference not in REFERENCE_RULES_ALL:
            # The name is still checked, so a typo is not swallowed silently.
            resolve_reference(params.reference, genomes, outgroup_file, logger)
            logger.info(
                "SNP typer '%s' is reference-free; --reference %s is not used",
                params.snptyper,
                params.reference,
            )
        return None
    return resolve_reference(
        params.reference,
        genomes,
        outgroup_file,
        logger,
        sketch_provider=sketch_provider,
        scratch_dir=dirs.scratch_dir,
        threads=params.threads,
    )


def _build_msa(
    genomes: list[Path],
    outgroup_file: Path | None,
    dirs: PhyloDirs,
    params: PhyloParams,
    logger: logging.Logger,
    reference: Path | None = None,
) -> tuple[Path, dict[str, str]]:
    inputs = list(genomes)
    if outgroup_file is not None:
        inputs.append(outgroup_file)

    if params.msa_source == "aligner":
        warn = scale_warning(aligner_registry, params.aligner, len(inputs))
        if warn:
            limit, alts = warn
            logger.warning(
                "Aligner '%s' is tuned for <=%d genomes but you have %d; consider: %s",
                params.aligner,
                limit,
                len(inputs),
                ", ".join(alts) or "none",
            )
        aligner = aligner_registry.create(params.aligner)
        versions = aligner.preflight()
        if reference is None:
            choice = _reference_for_source(genomes, outgroup_file, dirs, params, logger, None)
            assert choice is not None
            reference = choice.genome
        _warn_divergence(params.aligner, inputs, logger)
        # An aligner that cannot read gzip aligns decompressed copies, named by
        # record name; the MSA key and stamp keep using the original paths.
        inputs_dir = dirs.scratch_dir / _STAGED_INPUTS
        remove_tree(inputs_dir)  # copies a killed run left behind
        staged = stage_plain_inputs(inputs, aligner.capabilities, inputs_dir, logger)
        align_params = AlignParams(
            threads=params.threads,
            reference=staged[reference],
            extra=_adapter_extra(params.extra),
        )
        warn_ignored_params(aligner.capabilities, align_params, logger, family="Aligner")
        try:
            result = aligner.align(
                [staged[p] for p in inputs],
                staged[reference],
                dirs.align_dir,
                align_params,
                logger,
            )
        finally:
            remove_tree(inputs_dir)
        return result.msa_fasta, versions

    if params.msa_source == "snptype":
        # Reuse the SNP typer's core-SNP alignment as the MSA source. It is
        # written under tree/msa/; snp/ is the snptype stage's directory.
        from .snptype import SnptypeParams, snptype_core

        snp_params = SnptypeParams(
            tool=params.snptyper,
            threads=params.threads,
            reference=params.reference,
            all_genomes=params.all_genomes,
            mask=params.extra.get("mask", "none"),
            core_fraction=_core_fraction(params.extra),
            extra=_adapter_extra(params.extra),
            # A species-level outgroup breaks the recombination scan; the
            # masker runs on the ingroup and applies its regions to all.
            # Every typer names records by record_name.
            mask_exclude=(record_name(outgroup_file),) if outgroup_file is not None else (),
        )
        # None lets a reference-free typer (ska2) ignore it; a mapping typer
        # gets the medoid, a name, or the first genome, as the aligner path does.
        snp_reference: Path | None = reference

        # The outgroup must be typed with the ingroup (`inputs`, as on the
        # aligner path); typing `genomes` alone left it out of the alignment,
        # so the tree had no outgroup leaf to root on.
        snp_result, versions = snptype_core(
            inputs,
            snp_reference,
            dirs.msa_dir,
            dirs.scratch_dir / _TYPING_SCRATCH,
            snp_params,
            logger,
            # build_tree already warned once for the typer and the tree builder
            # together; a second per-typer warning would name keys the builder reads.
            warn_extras=False,
        )
        return snp_result.core_snp_fasta, versions

    if params.msa_source == "coregene":
        from ..coregenes.base import CoreGeneParams
        from .coregene import annotated_set, coregene_core

        if dirs.annotations_dir is None:
            raise UserInputError(
                "--msa-source coregene needs the annotations the annotate stage writes; none "
                "are available on this path."
            )
        names = [record_name(g) for g in inputs]
        annotated = annotated_set(dirs.annotations_dir, names)
        cg_result, versions = coregene_core(
            annotated.gffs,
            dirs.msa_dir,
            dirs.scratch_dir / "phylo_coregene",
            CoreGeneParams(
                threads=params.threads,
                core_threshold=_core_threshold(params.extra),
                clean_mode=str(params.extra.get("clean_mode", "strict")),
                filtered=_filtered(params.extra),
                extra=_adapter_extra(params.extra),
            ),
            params.coregene,
            logger,
            warn_extras=False,
            required=(record_name(outgroup_file),) if outgroup_file is not None else (),
            selection=_gene_selection(params.extra),
        )
        assert cg_result.alignment is not None  # coregene_core always publishes one
        return cg_result.alignment, versions

    raise UserInputError(f"Unknown msa-source '{params.msa_source}' (aligner|snptype|coregene)")


def _taxonomic_spread(genomes: Sequence[Path]) -> tuple[int, int]:
    """Distinct (genera, species) among the inputs; see
    :func:`repgenr.core.contracts.taxonomic_spread`."""
    return taxonomic_spread(genomes)


_LOW_DIVERSITY_SITES = 10  # fewer variable columns than this makes a meaningless tree


def _warn_low_diversity(msa_path: Path, logger: logging.Logger) -> None:
    """Warn when an MSA is (near-)invariant, mirroring _warn_divergence.

    A clone-only input produces an alignment with essentially no variable
    columns; every builder then emits a star-like tree with near-zero branch
    lengths, silently. Streams sequences one at a time and compares each to
    the first, stopping as soon as enough variable sites are seen, so the
    diverse (common) case exits early.
    """
    first: list[str] | None = None
    variable: set[int] = set()
    current: list[str] = []

    def _consume(seq: str) -> bool:
        nonlocal first
        if first is None:
            first = list(seq)
            return False
        for i, (a, b) in enumerate(zip(first, seq, strict=False)):
            if a != b:
                variable.add(i)
                if len(variable) >= _LOW_DIVERSITY_SITES:
                    return True
        return False

    try:
        with open(msa_path, encoding="utf-8") as fo:
            for line in fo:
                if line.startswith(">"):
                    if current and _consume("".join(current)):
                        return
                    current = []
                else:
                    current.append(line.strip())
            if current and _consume("".join(current)):
                return
    except OSError:
        return  # diagnostics must never fail the stage
    logger.warning(
        "The MSA at %s has only %d variable site(s): the input genomes are "
        "nearly identical, and the resulting tree will be star-like with "
        "near-zero branch lengths. Check whether the set is a single clone "
        "before interpreting the topology.",
        msa_path.name,
        len(variable),
    )


def _warn_divergence(aligner_name: str, genomes: Sequence[Path], logger) -> None:
    """Warn when a whole-genome aligner is run on a divergent (genus/family-level)
    set, where the shared collinear core shrinks and the alignment degrades.
    """
    n_genera, n_species = _taxonomic_spread(genomes)
    if aligner_name == "cactus" and n_species > 1:
        logger.warning(
            "Aligner 'cactus' (Minigraph-Cactus) targets same-species genomes, but the input "
            "spans %d species; it will likely drop divergent genomes from the graph. For "
            "genus/family-level data use an alignment-free tree builder (mashtree/sourmash).",
            n_species,
        )
    elif n_genera > 1:
        logger.warning(
            "Whole-genome aligner '%s' on a family-level set (%d genera): the shared collinear "
            "core shrinks sharply with divergence, so the alignment may be small or fragmentary. "
            "Consider an alignment-free tree builder (mashtree/sourmash), or loosen the aligner "
            "seeds (e.g. --aligner-arg kmer=15 for sibeliaz).",
            aligner_name,
            n_genera,
        )
    elif n_species > 1:
        logger.info(
            "Whole-genome aligner '%s' on a genus-level set (%d species): expect a reduced core "
            "alignment as divergence increases; alignment-free builders scale better.",
            aligner_name,
            n_species,
        )
