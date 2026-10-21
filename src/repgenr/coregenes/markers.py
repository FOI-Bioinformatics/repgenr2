"""GTDB marker-gene builder: single-copy marker proteins found by HMM search.

Command surface checked against HMMER 3.4. The builder needs no orthology
clustering: it searches each genome's proteins with the marker families of
GTDB's reference tree (bac120 or ar53, from the database ``repgenr
markers-db download`` writes) and aligns each marker's proteins to its
model, so it reaches across families and orders. It follows GTDB-Tk's
rules (``gtdbtk/markers.py`` and ``gtdbtk/pipeline/align.py``):

1. Search. ``hmmsearch --tblout HITS --noali --notextw --cut_ga`` over the
   Pfam models and ``--cut_nc`` over the TIGRFAM models, one pair of calls
   per genome. The database keeps one HMM file per set, so the builder
   splits it by source in scratch (Pfam and TIGRFAM files, and one file
   per model for the alignment step).
2. One model per gene. A protein that several models hit is assigned to
   the one with the highest bit score (then the smaller E-value, then the
   accession).
3. One gene per marker. A marker found once in a genome is *unique*. Found
   several times with the same protein sequence each time it is *multiple
   unique* and the best hit is kept. Found several times with different
   sequences it is *multiple* and treated as absent, as GTDB-Tk does. Not
   found it is *missing*.
4. Alignment. ``hmmalign --amino --outformat Pfam`` aligns all the genomes'
   proteins of one marker at once, and only the columns that the ``#=GC
   RF`` line marks as match states are kept, so each marker gives as many
   columns as its model has match states. GTDB-Tk aligns one sequence at a
   time; the match-state columns are the same either way.

A genome with fewer than ``min_markers`` single-copy markers (10 % of the
set by default, GTDB-Tk's ``--min_perc_aa 10`` in markers rather than
columns) is refused by name before any alignment; one below half of the
set is named in a warning. ``marker_set=auto`` (the default) searches
bac120 first, and when fewer than half of the genomes reach its minimum it
searches ar53 as well and keeps the set with the higher median fraction of
markers found. A marker enters the alignment when it is single-copy in at
least ``--core-threshold`` of the genomes; the selection, partitions and
trimming of the stage then apply as for the other builders. The proteins
are the annotate stage's (bakta or prokka); GTDB-Tk calls genes with
Prodigal, so the proteins can differ slightly at gene boundaries.
"""

from __future__ import annotations

import logging
import math
import os
import statistics
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_chain
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities, parse_extra_int
from .base import (
    CoreGeneBuilder,
    CoreGeneParams,
    CoreGeneResult,
    GeneAlignment,
    annotation_sidecars,
)
from .markers_db import (
    DOWNLOAD_HINT,
    MARKER_SETS,
    _split_models,
    check_markers_db,
    hmm_path,
    load_marker_list,
    markers_db_version,
    read_index,
    read_version,
)

MARKER_SET_CHOICES = (*MARKER_SETS, "auto")
SUMMARY_TSV = "markers_summary.tsv"
PER_GENOME_TSV = "markers_per_genome.tsv"
STATUSES = ("unique", "multiple_unique", "multiple", "missing")
# GTDB-Tk's --min_perc_aa default, applied here to the markers found.
DEFAULT_MIN_FRACTION = 0.10
# A genome below this fraction of the set's markers is named in a warning.
WARN_FRACTION = 0.50
# Genomes searched, or markers aligned, per container start.
_SEARCH_BATCH = 25
_ALIGN_BATCH = 50


@dataclass(frozen=True)
class Hit:
    """One row of a hmmsearch --tblout table."""

    gene: str
    marker: str  # the model accession without its version
    evalue: float
    score: float


@dataclass
class SetSearch:
    """The outcome of searching every genome with one marker set."""

    marker_set: str
    markers: list[str]
    min_markers: int
    # genome -> marker -> (status, kept gene or "", bit score of the kept hit)
    calls: dict[str, dict[str, tuple[str, str, float]]]

    def single_copy(self, genome: str) -> int:
        return sum(1 for status, _, _ in self.calls[genome].values() if status in STATUSES[:2])

    def counts(self, genome: str) -> dict[str, int]:
        out = dict.fromkeys(STATUSES, 0)
        for status, _, _ in self.calls[genome].values():
            out[status] += 1
        return out

    def median_fraction(self) -> float:
        return statistics.median(self.single_copy(g) / len(self.markers) for g in self.calls)


class MarkersBuilder(CoreGeneBuilder):
    capabilities = ToolCapabilities(
        name="markers",
        container="quay.io/biocontainers/hmmer:3.4--h7d74f8d_5",
        conda=("bioconda::hmmer",),
        required_binaries=(
            # The HMMER programs print their version in the help banner.
            BinarySpec("hmmsearch", version_args=("-h",), min_version="3.3"),
            BinarySpec("hmmalign", version_args=("-h",), min_version="3.3"),
        ),
        # Declared, not measured past 24 genomes: the search is linear in the
        # genome count and each marker is aligned in one call.
        recommended_max_genomes=20000,
        auto_eligible=False,
        accepted_extras=frozenset({"marker_set", "min_markers"}),
        # Panaroo's graph cleaning and entropy filter have no counterpart.
        ignored_params=frozenset({"clean_mode", "filtered"}),
    )
    # The marker models are protein profiles.
    required_alphabet = "protein"

    def call(
        self,
        gffs: Sequence[Path],
        out_dir: Path,
        params: CoreGeneParams,
        logger: logging.Logger,
    ) -> CoreGeneResult:
        if params.alphabet != "protein":
            raise UserInputError(
                "The markers core-gene builder aligns proteins to the marker models, so it "
                f"cannot hand over {params.alphabet} alignments; use --coregene-alphabet "
                "protein (the default for this builder)."
            )
        if not 0.0 < params.core_threshold <= 1.0:
            raise UserInputError(
                f"--core-threshold must be within (0, 1], got {params.core_threshold}"
            )
        requested = str(params.extra.get("marker_set", "auto")).strip().lower() or "auto"
        if requested not in MARKER_SET_CHOICES:
            raise UserInputError(
                f"marker_set must be one of {', '.join(MARKER_SET_CHOICES)}, not {requested!r}"
            )
        min_given = (
            parse_extra_int(params.extra, "min_markers", 0)
            if "min_markers" in params.extra
            else None
        )
        if min_given is not None and min_given < 1:
            raise UserInputError(f"min_markers must be 1 or more, got {min_given}")
        db = check_markers_db(params.markers_db)
        available = [s for s in MARKER_SETS if s in read_version(db).get("sets", [])]
        if not available or (requested != "auto" and requested not in available):
            raise UserInputError(
                f"The marker database {db} holds {', '.join(available) or 'no set'}, not "
                f"{requested}; add it with `repgenr markers-db download --output {db} --set "
                f"{'all' if requested == 'auto' else requested}` (the families present are kept)."
            )
        if min_given is not None:
            # Checked before any search: a set smaller than min_markers can
            # never pass the guard. auto leaves such a set out, with a warning.
            sizes = {s: len(load_marker_list(s)) for s in available}
            searchable = [s for s in available if sizes[s] >= min_given]
            if requested != "auto" and requested not in searchable:
                raise UserInputError(
                    f"min_markers is {min_given}, but the {requested} set has "
                    f"{sizes[requested]} markers."
                )
            if not searchable:
                raise UserInputError(
                    f"min_markers is {min_given}, more than any marker set holds "
                    f"({', '.join(f'{s}: {sizes[s]}' for s in available)})."
                )
            for s in available:
                if s not in searchable and requested == "auto":
                    logger.warning(
                        "min_markers is %d, more than the %d markers of %s; auto does not "
                        "consider %s.",
                        min_given,
                        sizes[s],
                        s,
                        s,
                    )
            available = searchable
        faas, _ = annotation_sidecars(gffs, "markers")
        names = [g.stem for g in gffs]
        genome_faas = dict(zip(names, faas, strict=True))
        out_dir.mkdir(parents=True, exist_ok=True)
        index = read_index(db)

        def search(marker_set: str) -> SetSearch:
            return _search_set(
                marker_set, db, index, genome_faas, out_dir, self.capabilities,
                params.threads, min_given, logger,
            )  # fmt: skip

        chosen = _choose_set(requested, available, search, len(names), logger)
        _write_tables(chosen, names, out_dir)
        _guard(chosen, names, out_dir, logger)

        threshold_n = math.ceil(params.core_threshold * len(names) - 1e-9)
        kept: list[str] = []
        for marker in chosen.markers:
            n = sum(1 for g in names if chosen.calls[g][marker][0] in STATUSES[:2])
            if n and n >= threshold_n:
                kept.append(marker)
        if not kept:
            raise WorkdirError(
                f"No {chosen.marker_set} marker is single-copy in {params.core_threshold:g} of "
                f"the genomes (--core-threshold); see {out_dir / PER_GENOME_TSV}. Lower the "
                "threshold, or leave the genomes with the fewest markers out."
            )
        dropped = sum(
            1
            for m in chosen.markers
            if m not in kept and any(chosen.calls[g][m][0] in STATUSES[:2] for g in names)
        )
        if dropped:
            logger.info(
                "%d %s marker(s) are single-copy in fewer than %g of the genomes "
                "(--core-threshold) and are left out of the alignment",
                dropped,
                chosen.marker_set,
                params.core_threshold,
            )
        lengths = {
            row["accession"]: int(row["length"]) for row in index if row["set"] == chosen.marker_set
        }
        gene_alignments = _align_markers(
            kept, chosen, genome_faas, lengths, out_dir, self.capabilities, params.threads,
            logger,
        )  # fmt: skip
        singles = [chosen.single_copy(g) for g in names]
        return CoreGeneResult(
            tables=[out_dir / SUMMARY_TSV, out_dir / PER_GENOME_TSV],
            core_genes=len(kept),
            settings={
                "marker_set": chosen.marker_set,
                "marker_set_requested": requested,
                "markers_in_set": len(chosen.markers),
                "min_markers": chosen.min_markers,
                "single_copy_markers_median": statistics.median(singles),
                "single_copy_markers_min": min(singles),
                "markers_db_version": markers_db_version(db, sets=(chosen.marker_set,)),
                "clustering_on": "GTDB marker HMMs",
            },
            gene_alignments=gene_alignments,
        )


def default_min_markers(n_markers: int) -> int:
    """10 % of the set's markers, rounded up."""
    return max(1, math.ceil(DEFAULT_MIN_FRACTION * n_markers - 1e-9))


def _choose_set(
    requested: str,
    available: Sequence[str],
    search,  # noqa: ANN001 - Callable[[str], SetSearch]
    n_genomes: int,
    logger: logging.Logger,
) -> SetSearch:
    """The searched set: the requested one, or the auto choice between them."""
    if requested != "auto":
        return search(requested)
    first = "bac120" if "bac120" in available else available[0]
    result = search(first)
    reach = sum(1 for g in result.calls if result.single_copy(g) >= result.min_markers)
    other = next((s for s in available if s != first), None)
    if 2 * reach >= n_genomes or other is None:
        if 2 * reach < n_genomes:
            logger.warning(
                "Only %d of %d genomes reach the %s minimum of %d markers, and the marker "
                "database holds no other set; download it with `repgenr markers-db download "
                "--output DIR --set all` for archaea.",
                reach,
                n_genomes,
                first,
                result.min_markers,
            )
        logger.info(
            "Marker set %s chosen automatically: %d of %d genomes reach its minimum of %d markers",
            first,
            reach,
            n_genomes,
            result.min_markers,
        )
        return result
    alt = search(other)
    a, b = result.median_fraction(), alt.median_fraction()
    chosen = alt if b > a else result
    logger.info(
        "Only %d of %d genomes reach the %s minimum of %d markers; the median fraction of "
        "markers found is %.2f for %s and %.2f for %s, so %s is used",
        reach,
        n_genomes,
        first,
        result.min_markers,
        a,
        first,
        b,
        other,
        chosen.marker_set,
    )
    return chosen


def _split_set(
    marker_set: str, db: Path, index: Sequence[Mapping[str, str]], out_dir: Path
) -> tuple[list[str], dict[str, Path], dict[str, Path]]:
    """Split the set's HMM file into one file per cutoff and one per model.

    Returns the set's accessions in list order, {cutoff: file} for the
    non-empty cutoff groups, and {accession: model file}.
    """
    markers = load_marker_list(marker_set)
    cutoff = {r["accession"]: r["cutoff"] for r in index if r["set"] == marker_set}
    source = hmm_path(db, marker_set)
    try:
        text = source.read_text(encoding="utf-8")
    except OSError as exc:
        raise WorkdirError(f"{source} is unreadable: {exc}") from exc
    models: dict[str, str] = {}
    for block in _split_models(text):
        acc = next(
            (ln.split()[1].split(".")[0] for ln in block.splitlines() if ln.startswith("ACC ")),
            "",
        )
        if acc:
            models[acc] = block
    missing = [m for m in markers if m not in models or m not in cutoff]
    if missing:
        raise WorkdirError(
            f"The marker database {db} lacks {len(missing)} {marker_set} model(s) "
            f"({', '.join(missing[:5])}); {DOWNLOAD_HINT} again."
        )
    model_dir = out_dir / "hmms" / marker_set
    model_dir.mkdir(parents=True, exist_ok=True)
    per_model: dict[str, Path] = {}
    groups: dict[str, list[str]] = {}
    for acc in markers:
        path = model_dir / f"{acc}.hmm"
        path.write_text(models[acc], encoding="utf-8")
        per_model[acc] = path
        groups.setdefault(cutoff[acc], []).append(models[acc])
    by_cutoff: dict[str, Path] = {}
    for cut, blocks in sorted(groups.items()):
        path = out_dir / "hmms" / f"{marker_set}_{cut}.hmm"
        path.write_text("".join(blocks), encoding="utf-8")
        by_cutoff[cut] = path
    return markers, by_cutoff, per_model


def _search_set(
    marker_set: str,
    db: Path,
    index: Sequence[Mapping[str, str]],
    faas: Mapping[str, Path],
    out_dir: Path,
    caps: ToolCapabilities,
    threads: int,
    min_given: int | None,
    logger: logging.Logger,
) -> SetSearch:
    """Search every genome with ``marker_set`` and classify its markers."""
    markers, by_cutoff, _ = _split_set(marker_set, db, index, out_dir)
    min_markers = min_given if min_given is not None else default_min_markers(len(markers))
    if min_markers > len(markers):
        raise UserInputError(
            f"min_markers is {min_markers}, but the {marker_set} set has {len(markers)} markers."
        )
    search_dir = out_dir / "search" / marker_set
    search_dir.mkdir(parents=True, exist_ok=True)
    genomes = list(faas)
    tables = {g: {cut: search_dir / f"{g}.{cut}.tbl" for cut in by_cutoff} for g in genomes}
    size = max(1, min(_SEARCH_BATCH, math.ceil(len(genomes) / max(1, threads))))
    batches = [genomes[i : i + size] for i in range(0, len(genomes), size)]

    def run(batch: Sequence[str]) -> None:
        steps: list[tuple[str, list[str | os.PathLike[str]]]] = []
        for g in batch:
            for cut, hmm in by_cutoff.items():
                steps.append(
                    (
                        "hmmsearch",
                        [
                            "hmmsearch",
                            "--tblout",
                            tables[g][cut],
                            "--noali",
                            "--notextw",
                            f"--cut_{cut}",
                            "--cpu",
                            "1",
                            hmm,
                            os.path.abspath(faas[g]),
                        ],  # fmt: skip
                    )
                )
        run_chain(caps, steps, logger=logger, cwd=out_dir)

    logger.info(
        "Searching %d genomes for the %d %s markers with hmmsearch",
        len(genomes),
        len(markers),
        marker_set,
    )
    parallel_map(run, batches, max(1, threads), logger=logger)
    calls: dict[str, dict[str, tuple[str, str, float]]] = {}
    for g in genomes:
        hits: list[Hit] = []
        for path in tables[g].values():
            if not path.is_file():
                raise WorkdirError(f"hmmsearch wrote no table for {g} ({path})")
            hits += read_tblout(path)
        calls[g] = classify(markers, hits, lambda ids, g=g: _read_proteins(faas[g], ids))
    return SetSearch(marker_set, markers, min_markers, calls)


def read_tblout(path: Path) -> list[Hit]:
    """The hits of a hmmsearch --tblout table.

    Columns (HMMER 3.4): target name, target accession, query name, query
    accession, full-sequence E-value, full-sequence score, ... The marker is
    the query accession without its version, or the query name when the
    model has no accession.
    """
    hits: list[Hit] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                raise WorkdirError(f"{path}: a table row has {len(fields)} columns: {line!r}")
            query = fields[3] if fields[3] != "-" else fields[2]
            try:
                hits.append(Hit(fields[0], query.split(".")[0], float(fields[4]), float(fields[5])))
            except ValueError as exc:
                raise WorkdirError(f"{path}: unreadable E-value or score in {line!r}") from exc
    return hits


def classify(
    markers: Sequence[str],
    hits: Sequence[Hit],
    proteins,  # noqa: ANN001 - Callable[[set[str]], dict[str, str]]
) -> dict[str, tuple[str, str, float]]:
    """GTDB-Tk's calls for one genome: marker -> (status, kept gene, score).

    Each gene keeps its best model; a marker with several genes is
    ``multiple_unique`` (best kept) when their protein sequences are equal,
    else ``multiple`` (none kept). ``proteins`` returns the sequences of a
    set of gene ids, read only when a marker has several genes.
    """
    wanted = set(markers)
    best: dict[str, Hit] = {}
    for hit in hits:
        if hit.marker not in wanted:
            continue
        cur = best.get(hit.gene)
        if cur is None or (-hit.score, hit.evalue, hit.marker) < (
            -cur.score,
            cur.evalue,
            cur.marker,
        ):
            best[hit.gene] = hit
    by_marker: dict[str, list[Hit]] = {}
    for hit in best.values():
        by_marker.setdefault(hit.marker, []).append(hit)
    several = {h.gene for hs in by_marker.values() if len(hs) > 1 for h in hs}
    seqs = proteins(several) if several else {}
    out: dict[str, tuple[str, str, float]] = {}
    for marker in markers:
        found = sorted(by_marker.get(marker, []), key=lambda h: (-h.score, h.evalue, h.gene))
        if not found:
            out[marker] = ("missing", "", 0.0)
        elif len(found) == 1:
            out[marker] = ("unique", found[0].gene, found[0].score)
        elif len({seqs.get(h.gene) for h in found}) == 1 and seqs.get(found[0].gene):
            out[marker] = ("multiple_unique", found[0].gene, found[0].score)
        else:
            out[marker] = ("multiple", "", 0.0)
    return out


def _read_proteins(path: Path, wanted: set[str]) -> dict[str, str]:
    """The sequences of ``wanted`` gene ids in ``path``, stop symbols removed."""
    seqs: dict[str, list[str]] = {}
    current: list[str] | None = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ident = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                current = seqs.setdefault(ident, []) if ident in wanted else None
            elif current is not None:
                current.append(line.strip())
    return {k: "".join(v).rstrip("*").replace("*", "X") for k, v in seqs.items()}


def _write_tables(result: SetSearch, names: Sequence[str], out_dir: Path) -> None:
    """markers_summary.tsv (counts per genome) and markers_per_genome.tsv."""
    n = len(result.markers)
    rows = ["\t".join(("genome", "set", *STATUSES, "single_copy_percent"))]
    for g in names:
        c = result.counts(g)
        pct = 100.0 * result.single_copy(g) / n
        rows.append("\t".join((g, result.marker_set, *(str(c[s]) for s in STATUSES), f"{pct:.1f}")))
    (out_dir / SUMMARY_TSV).write_text("\n".join(rows) + "\n", encoding="utf-8")
    rows = ["\t".join(("genome", "set", "marker", "status", "gene", "bitscore"))]
    for g in names:
        for marker in result.markers:
            status, gene, score = result.calls[g][marker]
            rows.append(
                "\t".join(
                    (g, result.marker_set, marker, status, gene, f"{score:.1f}" if gene else "")
                )
            )
    (out_dir / PER_GENOME_TSV).write_text("\n".join(rows) + "\n", encoding="utf-8")


def _guard(result: SetSearch, names: Sequence[str], out_dir: Path, logger: logging.Logger) -> None:
    """Refuse genomes below ``min_markers``; warn about those below half the set."""
    n = len(result.markers)
    low = [(g, result.single_copy(g)) for g in names if result.single_copy(g) < result.min_markers]
    if low:
        shown = ", ".join(f"{g} ({k})" for g, k in low[:10]) + (" ..." if len(low) > 10 else "")
        raise WorkdirError(
            f"{len(low)} genome(s) have fewer than {result.min_markers} of the {n} "
            f"{result.marker_set} markers in a single copy: {shown}. The marker set may be "
            "the wrong domain (try --marker-set ar53 or bac120), or the annotation is "
            "incomplete; leave the genome(s) out or lower min_markers. See "
            f"{out_dir / SUMMARY_TSV}."
        )
    half = [(g, result.single_copy(g)) for g in names if result.single_copy(g) < WARN_FRACTION * n]
    if half:
        logger.warning(
            "%d genome(s) have fewer than half of the %d %s markers in a single copy, so they "
            "rest on little data in the tree: %s",
            len(half),
            n,
            result.marker_set,
            ", ".join(f"{g} ({k})" for g, k in half[:10]) + (" ..." if len(half) > 10 else ""),
        )


def _align_markers(
    kept: Sequence[str],
    result: SetSearch,
    faas: Mapping[str, Path],
    lengths: Mapping[str, int],
    out_dir: Path,
    caps: ToolCapabilities,
    threads: int,
    logger: logging.Logger,
) -> list[GeneAlignment]:
    """Align each kept marker's proteins to its model; keep the match columns."""
    names = list(faas)
    wanted = {
        g: {result.calls[g][m][1] for m in kept if result.calls[g][m][0] in STATUSES[:2]}
        for g in names
    }
    proteins = {g: _read_proteins(faas[g], ids) for g, ids in wanted.items() if ids}
    seq_dir = out_dir / "markers"
    aln_dir = out_dir / "alignments"
    gene_dir = out_dir / "genes"
    for d in (seq_dir, aln_dir, gene_dir):
        d.mkdir(parents=True, exist_ok=True)
    model_dir = out_dir / "hmms" / result.marker_set
    jobs: list[tuple[str, Path, Path, list[str]]] = []
    for marker in kept:
        src = seq_dir / f"{marker}.faa"
        members: list[str] = []
        with open(src, "w", encoding="utf-8") as fo:
            for g in names:
                status, gene, _ = result.calls[g][marker]
                if status not in STATUSES[:2]:
                    continue
                seq = proteins[g].get(gene)
                if not seq:
                    raise WorkdirError(
                        f"hmmsearch named gene {gene!r} of {g} for {marker}, but the protein "
                        f"file {faas[g]} has no record of that name."
                    )
                fo.write(f">{g}\n{seq}\n")
                members.append(g)
        jobs.append((marker, src, aln_dir / f"{marker}.sto", members))
    batches = [jobs[i : i + _ALIGN_BATCH] for i in range(0, len(jobs), _ALIGN_BATCH)]

    def run(batch: Sequence[tuple[str, Path, Path, list[str]]]) -> None:
        run_chain(
            caps,
            [
                (
                    "hmmalign",
                    [
                        "hmmalign",
                        "--amino",
                        "--outformat",
                        "Pfam",
                        "-o",
                        dst,
                        model_dir / f"{marker}.hmm",
                        src,
                    ],  # fmt: skip
                )
                for marker, src, dst, _ in batch
            ],
            logger=logger,
            cwd=out_dir,
        )

    logger.info("Aligning %d markers to their models with hmmalign", len(jobs))
    parallel_map(run, batches, max(1, threads), logger=logger)
    out: list[GeneAlignment] = []
    for marker, _, dst, members in jobs:
        if not dst.is_file():
            raise WorkdirError(f"hmmalign wrote no alignment for {marker}")
        aligned = match_columns(dst)
        absent = sorted(set(members) - set(aligned))
        if absent:
            raise WorkdirError(f"hmmalign left {', '.join(absent[:3])} out of {marker} ({dst})")
        width = len(next(iter(aligned.values())))
        if marker in lengths and width != lengths[marker]:
            raise WorkdirError(
                f"The {marker} alignment has {width} match columns, but its model has "
                f"{lengths[marker]} match states ({dst})."
            )
        path = gene_dir / f"{marker}.fasta"
        with open(path, "w", encoding="utf-8") as fo:
            for g in members:
                fo.write(f">{g}\n{aligned[g]}\n")
        out.append(GeneAlignment(marker, path))
    return out


def match_columns(path: Path) -> dict[str, str]:
    """The match-state columns of a hmmalign Pfam-format alignment, by record.

    The ``#=GC RF`` line marks match states with ``x`` and insert columns
    with ``.``; only the match columns are kept, upper-cased, with ``-`` as
    the gap symbol.
    """
    rows: dict[str, str] = {}
    rf: str | None = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#=GC RF"):
                # Pfam format is one block (not interleaved), so there is a
                # single RF line and it replaces rather than extends.
                rf = line.split()[-1]
                continue
            if not line.strip() or line.startswith("#") or line.startswith("//"):
                continue
            name, _, seq = line.rstrip("\n").partition(" ")
            rows[name] = rows.get(name, "") + seq.strip()
    if rf is None:
        raise WorkdirError(f"{path} has no '#=GC RF' line; expected hmmalign --outformat Pfam.")
    keep = [i for i, c in enumerate(rf) if c in "xX"]
    out: dict[str, str] = {}
    for name, seq in rows.items():
        if len(seq) != len(rf):
            raise WorkdirError(
                f"{path}: the row of {name} has {len(seq)} columns, the RF line {len(rf)}."
            )
        out[name] = "".join(seq[i] for i in keep).upper().replace(".", "-")
    return out
