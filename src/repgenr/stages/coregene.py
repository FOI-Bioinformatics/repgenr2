"""The core-gene alignment pass of ``phylo --msa-source coregene``.

``coregene_core`` is stateless: it takes the GFF3 files of the annotated
genome set (the ``annotate`` stage's cache, or the ``annotate-build``
step's directories), runs a core-gene builder in scratch, and publishes the
concatenated alignment as ``tree/msa/core_gene_alignment.fasta`` with the
builder's tables under ``tree/msa/pangenome/`` and a small
``core_gene_stats.json``. The alignment keeps its invariant sites, so no
``alignment_stats.json`` is written and the tree builder applies no
ascertainment-bias correction.

A builder that returns one alignment per gene has its genes measured,
ranked and capped by a :class:`~repgenr.coregenes.select.GeneSelection`
and concatenated here; ``core_genes.tsv`` lists every gene with its
metrics and ``partitions.json`` the columns of each kept gene. A builder
that returns only its own concatenation has it published as it is.

``CoreGeneParams.alphabet`` chooses codon (``dna``) or amino-acid
(``protein``) gene alignments. It is recorded in ``partitions.json`` and
``core_gene_stats.json`` so the tree builder is given the matching model;
a builder that declares it does not read the alphabet is refused a protein
alignment rather than warned, because the tree would be fitted with a
model for the wrong character type.

``trim="clipkit"`` trims each per-gene alignment with ClipKIT
(:mod:`repgenr.coregenes.trim`) before the genes are measured, ranked and
concatenated; a builder that returns only its own concatenation is refused
it. Off by default.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from ..core.annotations import AnnotationRecord, gff_path, index_digest, read_index
from ..core.contracts import (
    CORE_GENE_FASTA,
    CORE_GENE_STATS,
    CORE_GENES_TSV,
    PANGENOME_DIR,
    atomic_path,
)
from ..core.errors import UserInputError, WorkdirError
from ..core.plugins import scale_warning, warn_ignored_params, warn_unconsumed_extras
from ..core.process import remove_tree
from ..coregenes.base import CoreGeneParams, CoreGeneResult, GeneAlignment
from ..coregenes.base import registry as coregene_registry
from ..coregenes.select import (
    ALPHABETS,
    GeneMetrics,
    GeneRow,
    GeneSelection,
    Partition,
    concatenate,
    gene_metrics,
    load_gene,
    rank_genes,
    write_gene_table,
)
from ..coregenes.trim import (
    DEFAULT_TRIM_MODE,
    alignment_width,
    check_gene_trim,
    preflight_clipkit,
    trim_gene_alignments,
)
from ..tree.partitions import PARTITIONS, write_partitions

# Every gene, ranked by occupancy: what a builder without per-gene output gets.
_DEFAULT_SELECTION = GeneSelection()


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


def check_alphabet(tool: str, alphabet: str) -> None:
    """Refuse an unknown alphabet, one other than the builder's
    ``required_alphabet`` (the markers builder: protein), or a non-dna one for
    a builder that lists ``alphabet`` in its ignored_params.

    A refusal rather than a warning: the tree model follows the alphabet, so
    an alignment of the other character type would be fitted with the wrong
    model. Called by the stage and, earlier, by the CLI.
    """
    if alphabet not in ALPHABETS:
        raise UserInputError(
            f"Unknown core-gene alphabet '{alphabet}'. Choose from: {', '.join(ALPHABETS)}."
        )
    if tool not in coregene_registry.names():
        return
    required = getattr(coregene_registry.get(tool), "required_alphabet", None)
    if required is not None and alphabet != required:
        raise UserInputError(
            f"The core-gene builder '{tool}' aligns genes as {required} only, so "
            f"--coregene-alphabet {alphabet} cannot apply; leave the alphabet out (it is "
            f"{required} for this builder) or set it to {required}."
        )
    if alphabet == "dna":
        return
    if "alphabet" not in coregene_registry.get(tool).capabilities.ignored_params:
        return
    able = [
        n
        for n in coregene_registry.names()
        if "alphabet" not in coregene_registry.get(n).capabilities.ignored_params
    ]
    raise UserInputError(
        f"The core-gene builder '{tool}' aligns genes as nucleotides only, so "
        f"--coregene-alphabet {alphabet} cannot apply; choose "
        f"{', '.join(able) or 'another builder'} for a protein alignment, or leave the "
        "alphabet at dna."
    )


def default_alphabet(tool: str) -> str:
    """The alphabet of ``tool`` when none is given: the one it requires, else dna."""
    if tool in coregene_registry.names():
        required = getattr(coregene_registry.get(tool), "required_alphabet", None)
        if required is not None:
            return required
    return "dna"


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
    selection: GeneSelection = _DEFAULT_SELECTION,
    trim: str = "none",
) -> tuple[CoreGeneResult, dict[str, str]]:
    """Run the core-gene builder ``tool`` over ``gffs`` into ``msa_dir``.

    ``required`` names records the tree cannot do without (the outgroup):
    one missing from the alignment is an error here rather than a rooting
    failure two steps later. Any other genome missing from the alignment is
    refused as well, since the tree must hold every input genome.
    ``selection`` chooses the genes when the builder returns per-gene
    alignments; on a builder that returns only its own concatenation
    anything but the default is refused. ``trim`` is
    ``none`` or ``clipkit``: with clipkit the per-gene alignments are trimmed
    before the selection (refused for a builder without per-gene output).
    """
    if tool not in coregene_registry.names():
        raise UserInputError(
            f"Unknown core-gene builder '{tool}'. Choose from: "
            f"{', '.join(coregene_registry.names())}."
        )
    builder = coregene_registry.create(tool)
    check_alphabet(tool, params.alphabet)
    check_gene_trim(trim)
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
    warn_ignored_params(builder.capabilities, params, logger, family="Core-gene builder")
    versions = builder.preflight()
    if trim != "none":
        versions.update(preflight_clipkit())
    # Builders pass paths under scratch to tools that run with cwd inside it,
    # so a relative scratch (phylo-build -o .) would not resolve.
    scratch = Path(os.path.abspath(scratch))
    remove_tree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    logger.info("Building the core-gene alignment of %d genomes with %s", len(gffs), tool)
    result = builder.call(gffs, scratch / tool, params, logger)

    msa_dir.mkdir(parents=True, exist_ok=True)
    alignment = msa_dir / CORE_GENE_FASTA
    expected_names = [g.stem for g in gffs]
    # A table or sidecar from an earlier run describes another alignment; the
    # sidecar's digest guard would ignore it, but a failed run should not
    # leave an old table beside a new alignment.
    (msa_dir / CORE_GENES_TSV).unlink(missing_ok=True)
    (msa_dir / PARTITIONS).unlink(missing_ok=True)
    chosen: _Chosen | None = None
    trimming: dict[str, int] | None = None
    if result.gene_alignments is not None:
        gene_alignments = result.gene_alignments
        if trim != "none":
            gene_alignments, trimming = _trim(
                gene_alignments, scratch / "trim", params.threads, tool, logger, params.alphabet
            )
        chosen = _choose_genes(
            gene_alignments,
            expected_names,
            selection,
            scratch,
            tool,
            logger,
            params.alphabet,
        )
        source = chosen.concatenation
    else:
        if trim != "none":
            raise UserInputError(
                f"The core-gene builder '{tool}' returns no per-gene alignments, so "
                f"--gene-trim {trim} cannot apply; leave it at none."
            )
        if selection != _DEFAULT_SELECTION:
            raise UserInputError(
                f"The core-gene builder '{tool}' returns no per-gene alignments, so a gene "
                "selection (maximum gene count, rank, minimum gene length or maximum gap "
                "fraction) cannot apply; leave those settings at their defaults."
            )
        if result.alignment is None:
            raise WorkdirError(
                f"The core-gene builder '{tool}' returned neither per-gene alignments nor an "
                "alignment."
            )
        source = result.alignment
    with atomic_path(alignment) as tmp:
        # A builder that writes every genome (PPanGGOLiN fills a missing family
        # with gaps), or the concatenation here, can leave a genome that shares
        # no selected gene as a row of gaps; dropped here, it meets the same
        # checks as an absent one.
        empty = _copy_without_gap_only_records(source, tmp)
    if empty:
        logger.info(
            "Removed %d record(s) of gaps only from the core-gene alignment: %s",
            len(empty),
            ", ".join(empty[:5]) + (" ..." if len(empty) > 5 else ""),
        )
    tables_dir = msa_dir / PANGENOME_DIR
    remove_tree(tables_dir)
    tables_dir.mkdir(parents=True, exist_ok=True)
    for table in result.tables:
        shutil.copy2(table, tables_dir / table.name)
    names, length = _shape(alignment)
    genomes = len(names)
    expected = set(expected_names)
    _refuse_unknown(names, expected)
    dropped = sorted(expected - set(names))
    missing_required = sorted(set(required) & set(dropped))
    if missing_required:
        raise WorkdirError(
            f"The core-gene alignment does not hold {', '.join(missing_required)}, which the "
            "tree needs (it is the outgroup): it shares no core gene with the rest at "
            f"--core-threshold {params.core_threshold:g}, or it holds none of the genes kept after "
            "the cap and filters (--max-core-genes, --min-gene-length, --max-gene-gaps). Lower the "
            "threshold, relax the selection, or build the tree without it "
            "(--no-outgroup --root midpoint)."
        )
    if dropped:
        # The tree is checked for every input genome after it is built
        # (check_tree_leaves), so a genome missing here would fail only after
        # the tree builder ran; refuse it before.
        raise WorkdirError(
            f"{len(dropped)} genome(s) are not in the core-gene alignment, so the tree would "
            f"lack them: {', '.join(dropped[:10])}{' ...' if len(dropped) > 10 else ''}. They "
            f"share no core gene at --core-threshold {params.core_threshold:g}, or hold none of "
            "the genes kept after the cap and filters. Lower --core-threshold, relax "
            "--max-core-genes, --min-gene-length or --max-gene-gaps, change --gene-rank, or "
            "leave the genome(s) out of the set."
        )
    core_genes = len(chosen.kept) if chosen is not None else result.core_genes
    if chosen is not None:
        write_gene_table(msa_dir / CORE_GENES_TSV, chosen.rows())
        write_partitions(alignment, chosen.partitions, params.alphabet)
    stats = {
        "alignment": alignment.name,
        "source": tool,
        # dna (codons) or protein (amino acids); sets the tree model.
        "alphabet": params.alphabet,
        "core_genes": core_genes,
        # Genes the builder aligned before the selection; None when it hands
        # over only its own concatenation.
        "core_genes_available": len(chosen.metrics) if chosen is not None else None,
        "selection": asdict(selection),
        # Per-gene trimming before the selection: none or clipkit, the ClipKIT
        # mode, the columns it removed over all genes (a dropped gene counts
        # all of its columns) and the genes it emptied. None when not trimmed.
        "trim": trim,
        "trim_mode": DEFAULT_TRIM_MODE if trimming is not None else None,
        "trimmed_columns_removed": trimming["columns_removed"] if trimming else None,
        "genes_dropped_by_trimming": trimming["genes_dropped"] if trimming else None,
        "variable_sites_total": chosen.variable_sites if chosen is not None else None,
        "informative_sites_total": chosen.informative_sites if chosen is not None else None,
        "gene_table": CORE_GENES_TSV if chosen is not None else None,
        "partitions": PARTITIONS if chosen is not None else None,
        "alignment_length": length,
        "genomes": genomes,
        "core_threshold": params.core_threshold,
        # The identity the builder clustered genes at; it governs how many
        # genes cluster at all.
        "clustering_identity": (
            result.clustering_identity
            if result.clustering_identity is not None
            # Only a builder that reads threshold= clusters at it.
            else params.extra.get("threshold")
            if "threshold" in builder.capabilities.accepted_extras
            else None
        ),
    }
    if result.settings:
        stats["builder_settings"] = dict(result.settings)
    # Settings the builder does not read are left out rather than recorded as
    # if they had shaped the result (Panaroo's cleaning under PPanGGOLiN).
    ignored = builder.capabilities.ignored_params
    for key in ("clean_mode", "filtered"):
        if key not in ignored:
            stats[key] = getattr(params, key)
    with atomic_path(msa_dir / CORE_GENE_STATS) as tmp:
        tmp.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    logger.info(
        "Core-gene alignment: %s genes, %d columns, %d genomes (invariant sites kept; no "
        "ascertainment-bias correction applies)",
        core_genes if core_genes is not None else "?",
        length,
        genomes,
    )
    remove_tree(scratch)
    return CoreGeneResult(
        alignment=alignment, tables=sorted(tables_dir.iterdir()), core_genes=core_genes
    ), versions


def _refuse_unknown(names: Sequence[str] | set[str], expected: set[str]) -> None:
    """Refuse record names that are none of the genomes of the set."""
    unknown = sorted(set(names) - expected)
    if unknown:
        raise WorkdirError(
            f"The core-gene alignment names {len(unknown)} record(s) that are none of the "
            f"genomes ({', '.join(unknown[:5])}): the builder did not name its samples after "
            "the annotation files, so the tree leaves would not match the genome set."
        )


def _trim(
    gene_alignments: Sequence[GeneAlignment],
    out_dir: Path,
    threads: int,
    tool: str,
    logger: logging.Logger,
    alphabet: str = "dna",
) -> tuple[list[GeneAlignment], dict[str, int]]:
    """Trim the per-gene alignments with ClipKIT; return them with the counts
    for the stats file (columns removed, genes dropped)."""
    before = sum(alignment_width(ga.path) for ga in gene_alignments)
    trimmed = trim_gene_alignments(gene_alignments, out_dir, threads, logger, alphabet=alphabet)
    after = sum(alignment_width(ga.path) for ga in trimmed)
    counts = {
        "columns_removed": before - after,
        "genes_dropped": len(gene_alignments) - len(trimmed),
    }
    if gene_alignments and not trimmed:
        raise WorkdirError(
            f"ClipKIT removed every column of all {len(gene_alignments)} core gene alignments "
            f"of '{tool}'; build the alignment without --gene-trim."
        )
    logger.info(
        "Gene trimming: %d of %d columns removed over %d genes; %d gene(s) dropped",
        counts["columns_removed"],
        before,
        len(gene_alignments),
        counts["genes_dropped"],
    )
    return trimmed, counts


@dataclass
class _Chosen:
    """The outcome of the gene selection: metrics of every gene, the kept
    ones in rank order, their columns in the concatenation written to
    ``concatenation`` (in scratch, before gap-only rows are dropped)."""

    metrics: dict[str, GeneMetrics]
    kept: list[str]
    partitions: list[Partition]
    concatenation: Path

    @property
    def variable_sites(self) -> int:
        return sum(self.metrics[g].variable_sites for g in self.kept)

    @property
    def informative_sites(self) -> int:
        return sum(self.metrics[g].informative_sites for g in self.kept)

    def rows(self) -> list[GeneRow]:
        """Table rows: kept genes in rank order, then the others by gene id."""
        where = {p.gene: p for p in self.partitions}
        rows = [
            GeneRow(g, True, i, self.metrics[g], where[g].start, where[g].end)
            for i, g in enumerate(self.kept, start=1)
        ]
        kept = set(self.kept)
        rows += [
            GeneRow(g, False, None, self.metrics[g], None, None)
            for g in sorted(self.metrics)
            if g not in kept
        ]
        return rows


def _choose_genes(
    gene_alignments: Sequence[GeneAlignment],
    genomes: Sequence[str],
    selection: GeneSelection,
    scratch: Path,
    tool: str,
    logger: logging.Logger,
    alphabet: str = "dna",
) -> _Chosen:
    """Measure, rank and concatenate the per-gene alignments of a builder.

    ``alphabet`` sets the missing-data symbol the metrics skip (N for dna,
    X for protein)."""
    seen: set[str] = set()
    for ga in gene_alignments:
        if ga.gene in seen:
            raise WorkdirError(
                f"The core-gene builder '{tool}' returned gene '{ga.gene}' more than once."
            )
        seen.add(ga.gene)
    genes = {ga.gene: load_gene(ga.path, alphabet) for ga in gene_alignments}
    _refuse_unknown({name for seqs in genes.values() for name in seqs}, set(genomes))
    metrics = {g: gene_metrics(seqs, len(genomes), alphabet) for g, seqs in genes.items()}
    kept = rank_genes(metrics, selection)
    if not kept:
        if not metrics:
            raise WorkdirError(f"The core-gene builder '{tool}' returned no core gene alignment.")
        raise WorkdirError(
            f"None of the {len(metrics)} core genes passes the gene filters (minimum length "
            f"{selection.min_length}, maximum gap fraction "
            f"{selection.max_gaps if selection.max_gaps is not None else 'none'}); relax them."
        )
    concatenation = scratch / "core_genes_concatenated.fasta"
    partitions = concatenate(kept, genes, genomes, concatenation)
    chosen = _Chosen(metrics, kept, partitions, concatenation)
    logger.info(
        "Core-gene selection: %d genes available, %d kept (rank %s); %d columns, %d variable "
        "and %d informative sites",
        len(metrics),
        len(kept),
        selection.rank,
        partitions[-1].end if partitions else 0,
        chosen.variable_sites,
        chosen.informative_sites,
    )
    return chosen


_GAPS = frozenset("-.?")


def _copy_without_gap_only_records(src: Path, dst: Path) -> list[str]:
    """Copy the FASTA ``src`` to ``dst`` leaving out records of gaps only.

    Streams one record at a time. Returns the names left out.
    """
    removed: list[str] = []
    with (
        open(src, encoding="utf-8", errors="replace") as fi,
        open(dst, "w", encoding="utf-8") as fo,
    ):
        header: str | None = None
        lines: list[str] = []

        def flush() -> None:
            if header is None:
                return
            if any(set(line.strip()) - _GAPS for line in lines):
                fo.write(header)
                fo.writelines(lines)
            else:
                removed.append(header[1:].split()[0] if header[1:].split() else "")

        for line in fi:
            if line.startswith(">"):
                flush()
                header, lines = line, []
            else:
                lines.append(line if line.endswith("\n") else line + "\n")
        flush()
    return removed


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
