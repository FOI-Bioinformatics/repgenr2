"""The core-gene alignment pass of ``phylo --msa-source coregene``.

``coregene_core`` is stateless: it takes the GFF3 files of the annotated
genome set (the ``annotate`` stage's cache, or the ``annotate-build``
step's directories), runs a core-gene builder in scratch, and publishes the
concatenated alignment as ``tree/msa/core_gene_alignment.fasta`` with the
builder's tables under ``tree/msa/pangenome/`` and a small
``core_gene_stats.json``. The alignment keeps its invariant sites, so no
``alignment_stats.json`` is written and the tree builder applies no
ascertainment-bias correction.
"""

from __future__ import annotations

import json
import logging
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.annotations import AnnotationRecord, gff_path, index_digest, read_index
from ..core.contracts import CORE_GENE_FASTA, CORE_GENE_STATS, PANGENOME_DIR, atomic_path
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import scale_warning, warn_unconsumed_extras
from ..core.process import remove_tree
from ..coregenes.base import CoreGeneParams, CoreGeneResult
from ..coregenes.base import registry as coregene_registry


@dataclass
class AnnotatedSet:
    """The annotation index rows of the genome set a tree is built from."""

    gffs: list[Path]
    annotator: str
    digest: str  # index_digest over the records


def annotated_set(annotations_dir: Path, names: Sequence[str]) -> AnnotatedSet:
    """The GFF3 files of ``names`` from the cache at ``annotations_dir``.

    Refuses when a record is missing (naming it and ``repgenr annotate``),
    when its GFF3 is gone, or when the records come from different
    annotators (one core-gene alignment, one annotation vocabulary).
    """
    records = read_index(annotations_dir)
    missing = [n for n in names if n not in records]
    if missing:
        shown = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        raise UserInputError(
            f"{len(missing)} genome(s) of the set have no annotation in {annotations_dir} "
            f"({shown}); run `repgenr annotate` (with --all-genomes for the whole set) first."
        )
    gffs: list[Path] = []
    annotators: set[str] = set()
    for name in names:
        rec: AnnotationRecord = records[name]
        gff = gff_path(annotations_dir, name)
        if not gff.is_file():
            raise WorkdirError(
                f"{gff} is listed in the annotation index but missing; rerun annotate."
            )
        gffs.append(gff)
        annotators.add(rec.annotator)
    if len(annotators) > 1:
        raise UserInputError(
            f"The annotations mix annotators ({', '.join(sorted(annotators))}); rerun "
            "`repgenr annotate --force` with one tool."
        )
    return AnnotatedSet(gffs=gffs, annotator=annotators.pop(), digest=index_digest(records, names))


def coregene_core(
    gffs: Sequence[Path],
    msa_dir: Path,
    scratch: Path,
    params: CoreGeneParams,
    tool: str,
    logger: logging.Logger,
    *,
    warn_extras: bool = True,
    required: Sequence[str] = (),
) -> tuple[CoreGeneResult, dict[str, str]]:
    """Run the core-gene builder ``tool`` over ``gffs`` into ``msa_dir``.

    ``required`` names records the tree cannot do without (the outgroup):
    one missing from the alignment is an error here rather than a rooting
    failure two steps later.
    """
    if tool not in coregene_registry.names():
        raise UserInputError(
            f"Unknown core-gene builder '{tool}'. Choose from: "
            f"{', '.join(coregene_registry.names())}."
        )
    builder = coregene_registry.create(tool)
    warn = scale_warning(coregene_registry, tool, len(gffs))
    if warn:
        limit, alts = warn
        logger.warning(
            "Core-gene builder '%s' is declared for <=%d genomes but you have %d (the limit "
            "is declared, not measured); consider: %s",
            tool,
            limit,
            len(gffs),
            ", ".join(alts) or "a smaller set",
        )
    if warn_extras:
        warn_unconsumed_extras(
            builder.capabilities, params.extra, logger, family="Core-gene builder"
        )
    versions = builder.preflight()
    remove_tree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    logger.info("Building the core-gene alignment of %d genomes with %s", len(gffs), tool)
    result = builder.call(gffs, scratch / tool, params, logger)

    msa_dir.mkdir(parents=True, exist_ok=True)
    alignment = msa_dir / CORE_GENE_FASTA
    with atomic_path(alignment) as tmp:
        shutil.copy2(result.alignment, tmp)
    tables_dir = msa_dir / PANGENOME_DIR
    remove_tree(tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    for table in result.tables:
        shutil.copy2(table, tables_dir / table.name)
    names, length = _shape(alignment)
    genomes = len(names)
    expected = {g.stem for g in gffs}
    unknown = sorted(set(names) - expected)
    if unknown:
        raise WorkdirError(
            f"The core-gene alignment names {len(unknown)} record(s) that are none of the "
            f"genomes ({', '.join(unknown[:5])}): the builder did not name its samples after "
            "the annotation files, so the tree leaves would not match the genome set."
        )
    dropped = sorted(expected - set(names))
    missing_required = sorted(set(required) & set(dropped))
    if missing_required:
        raise WorkdirError(
            f"The core-gene alignment does not hold {', '.join(missing_required)}, which the "
            "tree needs (it is the outgroup): it shares no core gene with the rest at "
            f"--core-threshold {params.core_threshold:g}. Lower the threshold, or build the "
            "tree without it (--no-outgroup --root midpoint)."
        )
    if dropped:
        logger.warning(
            "%d genome(s) are not in the core-gene alignment and will be missing from the "
            "tree (%s): they share no core gene at this threshold; lower --core-threshold "
            "or leave them out",
            len(dropped),
            ", ".join(dropped[:5]) + (" ..." if len(dropped) > 5 else ""),
        )
    stats = {
        "alignment": alignment.name,
        "source": tool,
        "core_genes": result.core_genes,
        "alignment_length": length,
        "genomes": genomes,
        "core_threshold": params.core_threshold,
        "clean_mode": params.clean_mode,
        "filtered": params.filtered,
    }
    with atomic_path(msa_dir / CORE_GENE_STATS) as tmp:
        tmp.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    logger.info(
        "Core-gene alignment: %s genes, %d columns, %d genomes (invariant sites kept; no "
        "ascertainment-bias correction applies)",
        result.core_genes if result.core_genes is not None else "?",
        length,
        genomes,
    )
    remove_tree(scratch)
    return CoreGeneResult(
        alignment=alignment, tables=sorted(tables_dir.iterdir()), core_genes=result.core_genes
    ), versions


def _shape(alignment: Path) -> tuple[list[str], int]:
    """(record names, columns of the first record) of a FASTA alignment."""
    names: list[str] = []
    length = 0
    with open(alignment, encoding="utf-8", errors="replace") as fo:
        for line in fo:
            if line.startswith(">"):
                names.append(line[1:].split()[0] if line[1:].split() else "")
            elif len(names) == 1:
                length += len(line.strip())
    return names, length
