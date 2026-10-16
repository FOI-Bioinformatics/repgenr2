"""Canonical inter-stage file contracts.

Each stage publishes a documented, validated artifact set that the next stage
consumes. These are internal engineering contracts (no legacy backward-compat
obligation), so the names and layout are chosen fresh:

    derep/representatives/      representative genome FASTAs
    derep/clusters.tsv          representative<TAB>member (member==representative for self)
    derep/genome_status.tsv     genome<TAB>status(representative|contained|fail_qc)
    derep/cluster_summary.tsv   one row per representative: size, species, quality
    align/msa.fasta             multiple sequence alignment (aligner output)
    snp/core_snp.fasta          variant-site alignment (snptype stage output)
    tree/msa/core_snp.fasta     variant-site alignment phylo typed for its tree
    tree/tree.nwk               Newick tree
    tree2tax.tsv                child<TAB>parent (FlexTaxD)
    genomes_map.tsv             accession<TAB>leaf

This module owns the writers/readers so producers (adapters) and consumers
(downstream stages) agree on one place.
"""

from __future__ import annotations

import csv
import os
import re
from collections import defaultdict
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

from .errors import WorkdirError

CLUSTERS_TSV = "clusters.tsv"
GENOME_STATUS_TSV = "genome_status.tsv"
CLUSTER_SUMMARY_TSV = "cluster_summary.tsv"
# Written by dereplicate --refine and derep-refine: the lineage genomes.
REFINED_TSV = "refined.tsv"
SELECTION_TSV = "selection.tsv"
MSA_FASTA = "msa.fasta"
CORE_SNP_FASTA = "core_snp.fasta"
MASKER_TREE_NWK = "masker_tree.nwk"
TREE_NWK = "tree.nwk"
# The tree builder's output as written, kept beside tree.nwk when phylo rerooted it.
TREE_UNROOTED_NWK = "tree.unrooted.nwk"
# Under tree/: the alignment phylo's own SNP typing pass (--msa-source snptype)
# writes, with its reuse stamp. snp/ belongs to the snptype stage.
PHYLO_MSA_DIR = "msa"
ANNOTATIONS_DIR = "annotations"
ANNOTATIONS_TSV = "annotations.tsv"
CORE_GENE_FASTA = "core_gene_alignment.fasta"
CORE_GENE_STATS = "core_gene_stats.json"
PANGENOME_DIR = "pangenome"
TREE2TAX_TSV = "tree2tax.tsv"
GENOMES_MAP_TSV = "genomes_map.tsv"
TREE2TAX_NODES_TSV = "tree2tax_nodes.tsv"
# Segment-grouped viral isolates (vgenome --group-segments): the synthetic
# isolate token used as the genome's accession -> its member accessions, so
# tree2tax can list the real accessions under the isolate's leaf.
SEGMENTS_TSV = "segments.tsv"
OUTGROUP_ACCESSION_TXT = "outgroup_accession.txt"
# Accessions the genome stage requested but NCBI returned nothing for;
# the completeness guard excuses them (core.integrity).
MISSING_ACCESSIONS_TXT = "missing_accessions.txt"
# The reads-to-assembly entry path: the sequencing runs the reads stage
# selected, the per-assembly metrics, and the runs the assemble stage could not
# turn into a genome (excused by the completeness guard, like missing
# accessions).
READS_TSV = "reads.tsv"
ASSEMBLY_STATS_TSV = "assembly_stats.tsv"
EXCUSED_RUNS_TSV = "excused_runs.tsv"


# Recognised genome FASTA extensions, the gzip forms first so suffix stripping
# is unambiguous (``.fasta.gz`` before ``.fasta``). ``.fna.gz`` is what the NCBI
# FTP site delivers. One definition shared by every stage and adapter that
# lists or names genome files.
FASTA_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fa", ".fna", ".fas")


def newick_is_complete(text: str) -> bool:
    """True when a Newick text holds exactly one tree, ended by its ';'.

    No ';' at all, text after the final ';', or a second tree before it marks
    a truncated or concatenated file: a Newick reader takes the first tree and
    ignores the rest. A ';' inside a quoted label or a [comment] does not end
    a tree. ``doctor`` and ``tree2tax`` share this rule.
    """
    content = text.strip()
    if not content.endswith(";"):
        return False
    terminators = 0
    quoted = False
    comment_depth = 0
    for char in content:
        if quoted:
            # A doubled quote inside a quoted label toggles twice: no net change.
            if char == "'":
                quoted = False
        elif comment_depth:
            if char == "[":
                comment_depth += 1
            elif char == "]":
                comment_depth -= 1
        elif char == "'":
            quoted = True
        elif char == "[":
            comment_depth = 1
        elif char == ";":
            terminators += 1
    return terminators == 1 and not quoted and not comment_depth


def rename_fasta_records(src: Path, dst: Path, names: dict[str, str]) -> None:
    """Copy a FASTA from ``src`` to ``dst``, renaming the records in ``names``.

    Tools rename their samples (harvesttools adds the file extension, cactus
    turns '.' into '_'); records are renamed back to genome stems so that
    tree leaves, the outgroup lookup and the masker's outgroup exclusion see
    the input names. Streams line by line, since whole-genome alignments are
    large, and removes ``src``. A record not in ``names`` keeps its name.
    """
    with open(src, encoding="utf-8") as fi, open(dst, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith(">"):
                name = line[1:].strip()
                line = f">{names.get(name, name)}\n"
            fo.write(line)
    src.unlink()


def list_fasta(source: Path) -> list[Path]:
    """Sorted genome FASTA files directly under ``source``.

    Skips dotfiles (e.g. macOS ``._`` AppleDouble) and returns an empty list when
    ``source`` does not exist, so callers need no separate existence check.
    """
    if not source.exists():
        return []
    return sorted(
        p
        for p in source.iterdir()
        if not p.name.startswith(".") and p.name.endswith(FASTA_SUFFIXES)
    )


def genome_filename(family: str, genus: str, species: str, accession: str) -> str:
    """Canonical genome FASTA filename. One definition so the metadata selection,
    the genome download and every downstream stage agree on the same names.

    ``family``/``genus``/``species`` must be single tokens (no ``_``) so the name
    round-trips through :func:`parse_genome_filename`; the accession may contain
    underscores (e.g. ``GCF_000001.1``, ``NC_001802.1``).
    """
    return f"{family}_{genus}_{species}_{accession}.fasta"


def strip_fasta_suffix(name: str) -> str:
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def record_name(path: str | Path) -> str:
    """The name of a genome in alignments, trees and tables.

    The file name without its FASTA suffix, gzip included: ``x.fasta.gz`` and
    ``x.fasta`` are both ``x`` (``Path.stem`` would give ``x.fasta`` for the
    first). Every aligner, SNP typer and tree builder adapter names records and
    leaves with this, so one genome carries one name across the stages.
    """
    return strip_fasta_suffix(Path(path).name)


def sanitise_taxon_tokens(family: str, genus: str, species: str) -> tuple[str, str, str]:
    """Turn taxonomy names into the single tokens a canonical filename holds.

    The species drops its genus prefix (``Francisella tularensis`` ->
    ``tularensis``); every token loses its spaces and turns underscores into
    hyphens, so :func:`parse_genome_filename` splits the name back on ``_``.
    One rule for every producer (GTDB TSV and API, the reads path), so genomes
    of one species carry one token whatever selected them.
    """

    def clean(value: str) -> str:
        return value.replace(" ", "").replace("_", "-")

    species = species.replace(genus, "") if genus else species
    return clean(family), clean(genus), clean(species)


def taxonomic_spread(genomes: Sequence[Path]) -> tuple[int, int]:
    """Distinct (genera, species) among ``genomes``, read from the canonical
    ``Family_Genus_species_Accession.fasta`` filenames, so it works without a
    manifest (shared-workdir and data-channel paths alike)."""
    genera: set[str] = set()
    species: set[tuple[str, str]] = set()
    for g in genomes:
        _family, genus, sp, _acc = parse_genome_filename(Path(g).name)
        if genus:
            genera.add(genus.lower())
            species.add((genus.lower(), sp.lower()))
    return len(genera), len(species)


# An NCBI assembly accession leading a filename stem: GCA_/GCF_, digits, version.
_NCBI_ASSEMBLY_PREFIX = re.compile(r"^(GC[AF]_\d+\.\d+)(?:_|$)")


def parse_genome_filename(name: str) -> tuple[str, str, str, str]:
    """Inverse of :func:`genome_filename`. Returns (family, genus, species,
    accession). The first three ``_``-separated tokens are the taxonomy and
    **everything after** is the accession, so accessions with underscores
    (bacterial ``GCF_x.y``, viral ``NC_x.y``) and without (viral ``MN908947.3``)
    all round-trip. A name that starts with an NCBI assembly accession (the
    ``GCF_000008985.1_ASM898v1_genomic.fna`` files NCBI Datasets and the FTP site
    deliver) yields empty taxonomy and that accession. Any other non-canonical
    name (< 4 tokens) yields empty taxonomy and the whole stem as the accession.
    """
    stem = strip_fasta_suffix(Path(name).name)
    ncbi = _NCBI_ASSEMBLY_PREFIX.match(stem)
    if ncbi:
        return "", "", "", ncbi.group(1)
    parts = stem.split("_")
    if len(parts) < 4:
        return "", "", "", stem
    return parts[0], parts[1], parts[2], "_".join(parts[3:])


def accession_from_filename(name: str) -> str:
    """Recover the accession from a canonical genome filename or tree leaf."""
    return parse_genome_filename(name)[3]


@dataclass
class SelectionRow:
    """One selected genome: the portable hand-off from metadata to genome.

    Carries the taxonomy (so downstream taxonomy-aware stages need not re-read the
    SQLite manifest) and the canonical output filename used by the downloader.
    """

    accession: str
    family: str
    genus: str
    species: str
    is_outgroup: bool
    filename: str
    completeness: float | None = None
    contamination: float | None = None
    # GTDB marks one genome per species as its representative; the metadata
    # stage copies that flag here so --keeper gtdb can keep those genomes.
    # Every other entry path leaves it False.
    gtdb_representative: bool = False


def _require_columns(reader: csv.DictReader, path: Path, required: list[str]) -> None:
    """Raise WorkdirError naming the file when its header lacks a required column.

    An empty file (no header) passes; the caller reads it as holding no rows.
    """
    if reader.fieldnames is None:
        return
    missing = [c for c in required if c not in reader.fieldnames]
    if missing:
        raise WorkdirError(f"{path} lacks the column(s) {', '.join(missing)}.")


def _tsv_writer(fo: IO) -> Any:
    """A TSV writer with Unix line endings. ``csv.writer`` defaults to
    ``\\r\\n``, which leaves a stray ``\\r`` on the last column of every row
    for awk/cut consumers of the contract files."""
    return csv.writer(fo, delimiter="\t", lineterminator="\n")


@contextmanager
def atomic_replace(
    path: Path, *, mode: str = "w", encoding: str | None = "utf-8", newline: str | None = None
) -> Iterator[IO]:
    """Open a temp sibling of ``path`` for writing; publish it atomically.

    On clean exit the temp file replaces ``path`` (``os.replace``); on any
    exception the temp is removed and a previous ``path`` is left untouched, so
    a crash mid-write can never truncate a deliverable. Same pattern as
    ``Config.save``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        with open(tmp, mode, encoding=encoding, newline=newline) as fo:
            yield fo
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


@contextmanager
def atomic_path(path: Path) -> Iterator[Path]:
    """Yield a temp sibling path for tools/copies that need a filename.

    On clean exit the temp (if written) replaces ``path``; on exception it is
    removed and the previous ``path`` survives.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        yield tmp
        if tmp.exists():
            os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def write_selection(path: Path, rows: list[SelectionRow]) -> None:
    """Write the metadata selection (accession + taxonomy + filename + outgroup flag).

    ``gtdb_representative`` is the last column, so the positions of the
    earlier columns are those of files written before it was added.
    """
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(
            [
                "accession",
                "family",
                "genus",
                "species",
                "is_outgroup",
                "filename",
                "completeness",
                "contamination",
                "gtdb_representative",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.accession,
                    r.family,
                    r.genus,
                    r.species,
                    "1" if r.is_outgroup else "0",
                    r.filename,
                    "" if r.completeness is None else f"{r.completeness:.2f}",
                    "" if r.contamination is None else f"{r.contamination:.2f}",
                    "1" if r.gtdb_representative else "0",
                ]
            )


def _opt_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


_FLAGS = {"1": True, "true": True, "yes": True, "0": False, "false": False, "no": False}


def _outgroup_flag(value: str | None, column: str = "is_outgroup") -> bool:
    """A 1/0 flag column (is_outgroup, gtdb_representative): 1/0 as written,
    true/false and yes/no accepted, empty or absent read as 0.

    Any other value is an error: reading it as 0 put an intended outgroup into
    the ingroup without a message.
    """
    key = (value or "0").strip().lower() or "0"
    if key not in _FLAGS:
        raise ValueError(f"{column} must be 0 or 1, not {value!r}")
    return _FLAGS[key]


def read_selection(path: Path) -> list[SelectionRow]:
    """Read a selection.tsv back into SelectionRow records.

    A malformed file (undecodable, an unknown outgroup flag, a quality value
    that is not a number) raises WorkdirError naming the file and line.
    """
    rows: list[SelectionRow] = []
    try:
        with open(path, encoding="utf-8", newline="") as fo:
            reader = csv.DictReader(fo, delimiter="\t")
            _require_columns(reader, path, ["accession", "filename"])
            for row in reader:
                try:
                    for column in ("accession", "filename"):
                        if not (row.get(column) or "").strip():
                            raise ValueError(f"the {column} column is empty")
                    rows.append(
                        SelectionRow(
                            accession=row["accession"],
                            family=row.get("family") or "",
                            genus=row.get("genus") or "",
                            species=row.get("species") or "",
                            is_outgroup=_outgroup_flag(row.get("is_outgroup")),
                            filename=row["filename"],
                            completeness=_opt_float(row.get("completeness")),
                            contamination=_opt_float(row.get("contamination")),
                            gtdb_representative=_outgroup_flag(
                                row.get("gtdb_representative"), "gtdb_representative"
                            ),
                        )
                    )
                except ValueError as exc:
                    raise WorkdirError(f"{path} line {reader.line_num}: {exc}.") from exc
    except UnicodeDecodeError as exc:
        raise WorkdirError(f"{path} is not UTF-8 text ({exc.reason}).") from exc
    return rows


def write_clusters(path: Path, clusters: dict[str, list[str]]) -> None:
    """Write representative -> members. Each representative also lists itself."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["representative", "member"])
        for rep, members in clusters.items():
            writer.writerow([rep, rep])
            for member in members:
                if member != rep:
                    writer.writerow([rep, member])


def read_clusters(path: Path) -> dict[str, list[str]]:
    """Read representative -> members (members exclude the representative itself)."""
    clusters: dict[str, list[str]] = defaultdict(list)
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            rep, member = row[0], row[1]
            clusters.setdefault(rep, [])
            if member != rep:
                clusters[rep].append(member)
    return dict(clusters)


@dataclass(frozen=True)
class RefinedRow:
    """One genome of the lineage a dereplication refined at a finer ANI.

    ``seed`` is the seed whose coarse cluster brought the genome in,
    ``coarse_representative`` that cluster's representative before the
    refinement, ``representative`` the genome's representative after it (the
    genome itself when it is one; blank for ``fail_qc``), ``status`` its
    final status, ``ani`` the secondary ANI of the refinement (blank in files
    written before the column existed).
    """

    genome: str
    seed: str
    coarse_representative: str
    representative: str
    status: str
    ani: str = ""


_REFINED_COLUMNS = ("genome", "seed", "coarse_representative", "representative", "status", "ani")


def write_refined(path: Path, rows: Sequence[RefinedRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_REFINED_COLUMNS)
        for r in sorted(rows, key=lambda r: (r.coarse_representative, r.representative, r.genome)):
            writer.writerow(
                [r.genome, r.seed, r.coarse_representative, r.representative, r.status, r.ani]
            )


def read_refined(path: Path) -> list[RefinedRow]:
    """The refined rows; a missing file yields none."""
    if not path.is_file():
        return []
    rows: list[RefinedRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo, delimiter="\t"):
            rows.append(RefinedRow(*(rec.get(c, "") or "" for c in _REFINED_COLUMNS)))
    return rows


def write_genome_status(path: Path, status: dict[str, str]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["genome", "status"])
        for genome, value in sorted(status.items()):
            writer.writerow([genome, value])


def read_genome_status(path: Path) -> dict[str, str]:
    """Read genome -> status. A missing file yields an empty mapping.

    The counterpart of :func:`write_genome_status`, so a step that consumes a
    contract directory (Nextflow scatter-gather) recovers the statuses of
    genomes that belong to no cluster -- ``fail_qc`` genomes would otherwise be
    invisible to the consumer.
    """
    status: dict[str, str] = {}
    if not path.exists():
        return status
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            status[row[0]] = row[1]
    return status


@dataclass
class ClusterSummaryRow:
    """One dereplication cluster seen from its representative.

    ``n_members`` excludes the representative; ``n_genomes`` includes it.
    ``species`` lists the distinct species across representative and members,
    the representative's first, capped at five names plus "+N more".
    Quality columns are ``None`` when the manifest carried no CheckM values;
    the ``member_*`` extremes span scored members only. ``best_member`` is the
    genome the quality keeper chooses (``derep_keeper.choose_keeper``) and
    equals ``representative`` when the keeper is already that genome; it is
    empty when no genome in the cluster is scored. ``best_score`` is its
    keeper score and ``rep_n50`` the representative's N50.
    """

    representative: str
    n_members: int
    n_species: int
    species: str
    rep_completeness: float | None = None
    rep_contamination: float | None = None
    member_max_completeness: float | None = None
    member_min_contamination: float | None = None
    best_member: str = ""
    rep_n50: int | None = None
    best_score: float | None = None
    # Whether the representative is a GTDB species representative (from the
    # manifest or selection.tsv); False when the selection carries no flag.
    rep_is_gtdb_representative: bool = False
    # Whether the cluster comes from a refinement at a finer ANI
    # (``dereplicate --refine`` or ``derep-refine``).
    refined: bool = False

    @property
    def n_genomes(self) -> int:
        """Genomes in the cluster, the representative included."""
        return self.n_members + 1


_CLUSTER_SUMMARY_COLUMNS = (
    "representative",
    "n_members",
    "n_species",
    "species",
    "rep_completeness",
    "rep_contamination",
    "member_max_completeness",
    "member_min_contamination",
    "best_member",
    # Appended last so the positions of the earlier columns stay unchanged.
    "n_genomes",
    "rep_n50",
    "best_score",
    "rep_is_gtdb_representative",
    "refined",
)


def write_cluster_summary(path: Path, rows: list[ClusterSummaryRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_CLUSTER_SUMMARY_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.representative,
                    r.n_members,
                    r.n_species,
                    r.species,
                    _fmt_opt(r.rep_completeness),
                    _fmt_opt(r.rep_contamination),
                    _fmt_opt(r.member_max_completeness),
                    _fmt_opt(r.member_min_contamination),
                    r.best_member,
                    r.n_genomes,
                    "" if r.rep_n50 is None else r.rep_n50,
                    _fmt_opt(r.best_score),
                    "1" if r.rep_is_gtdb_representative else "0",
                    "1" if r.refined else "0",
                ]
            )


def _fmt_opt(value: float | None) -> str:
    return "" if value is None else repr(value)


def read_cluster_summary(path: Path) -> list[ClusterSummaryRow]:
    rows: list[ClusterSummaryRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo, delimiter="\t")
        for rec in reader:
            rows.append(
                ClusterSummaryRow(
                    representative=rec["representative"],
                    n_members=int(rec["n_members"]),
                    n_species=int(rec["n_species"]),
                    species=rec["species"],
                    rep_completeness=_opt_float(rec.get("rep_completeness")),
                    rep_contamination=_opt_float(rec.get("rep_contamination")),
                    member_max_completeness=_opt_float(rec.get("member_max_completeness")),
                    member_min_contamination=_opt_float(rec.get("member_min_contamination")),
                    best_member=rec.get("best_member", ""),
                    rep_n50=int(rec["rep_n50"]) if rec.get("rep_n50") else None,
                    best_score=_opt_float(rec.get("best_score")),
                    rep_is_gtdb_representative=_outgroup_flag(
                        rec.get("rep_is_gtdb_representative"), "rep_is_gtdb_representative"
                    ),
                    refined=_outgroup_flag(rec.get("refined"), "refined"),
                )
            )
    return rows


@dataclass(frozen=True)
class NodeRow:
    """One node of the taxonomy, with what the tree says about it.

    ``depth`` is the path length from the root, ``branch_length`` the edge to
    the parent, ``support`` the builder's value on the node (None for leaves
    and trees without supports), ``red`` the relative evolutionary divergence
    of Parks et al. (2018): 0 at the root (the ingroup's common ancestor
    when the outgroup was removed), 1 at every leaf, and for an inner node
    red(parent) + a/(a + b) * (1 - red(parent)) with a its branch length and
    b the mean distance from the node to its own leaves, as PhyloRank
    computes it. None where the tree has no branch lengths or no meaningful
    root. ``support`` is a fraction. ``lineage`` is the family;genus;species
    of a leaf's genome when the manifest holds it.
    """

    node: str
    parent: str | None
    n_leaves: int
    depth: float | None
    branch_length: float | None
    support: float | None
    red: float | None
    lineage: str | None = None


NODE_TABLE_COLUMNS = (
    "node",
    "parent",
    "n_leaves",
    "depth",
    "branch_length",
    "support",
    "red",
    "lineage",
)


def _na(value: float | None, digits: int) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def _g(value: float | None) -> str:
    # Six significant digits: branch lengths of clonal core-SNP trees sit
    # near 1e-6 and must not print as 0.000000.
    return "NA" if value is None else f"{value:.6g}"


def write_node_table(path: Path, rows: list[NodeRow]) -> None:
    """Write the node annotation table beside tree2tax.tsv (not read by FlexTaxD)."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(NODE_TABLE_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.node,
                    r.parent or "",
                    r.n_leaves,
                    _g(r.depth),
                    _g(r.branch_length),
                    _na(r.support, 4),
                    _na(r.red, 4),
                    r.lineage or "",
                ]
            )


def write_tree2tax(
    path: Path, edges: list[tuple[str, str]], ranks: dict[str, str] | None = None
) -> None:
    """Write child -> parent edges (FlexTaxD), de-duplicated, order preserved.

    With ``ranks`` (child -> rank) a third column ``rank`` is written. FlexTaxD
    reads it on every row and an empty value silently inherits the previous
    row's, so every child must have a non-empty rank; a missing one is an
    error here rather than a wrong taxonomy there.
    """
    seen: set[tuple[str, str]] = set()
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["child", "parent"] + (["rank"] if ranks is not None else []))
        for child, parent in edges:
            if (child, parent) in seen:
                continue
            seen.add((child, parent))
            if ranks is None:
                writer.writerow([child, parent])
                continue
            rank = ranks.get(child)
            if not rank:
                raise ValueError(f"tree2tax: node {child!r} has no rank for the rank column")
            writer.writerow([child, parent, rank])


@dataclass(frozen=True)
class ReadRow:
    """One sequencing run selected by the reads stage (the reads.tsv contract).

    ``fastq_urls``/``fastq_md5``/``fastq_bytes`` are parallel per-file tuples,
    empty when the archive holds no FASTQ mirror for the run. A URL may also
    be a local path, which the assemble stage copies instead of downloading.
    """

    run_accession: str
    biosample: str
    bioproject: str
    organism: str
    taxid: str
    platform: str
    instrument_model: str
    layout: str
    bases: int
    read_count: int
    # Sanitised filename tokens resolved from the taxid by the reads stage.
    family: str = ""
    genus: str = ""
    species: str = ""
    fastq_urls: tuple[str, ...] = ()
    fastq_md5: tuple[str, ...] = ()
    fastq_bytes: tuple[int, ...] = ()
    # ENA's library_selection (RANDOM, MDA, PCR, ...); MDA marks whole-genome
    # amplification, which assembles poorly.
    library_selection: str = ""


_READS_COLUMNS = [
    "run_accession",
    "biosample",
    "bioproject",
    "organism",
    "taxid",
    "platform",
    "instrument_model",
    "layout",
    "bases",
    "read_count",
    "family",
    "genus",
    "species",
    "fastq_urls",
    "fastq_md5",
    "fastq_bytes",
    "library_selection",
]


def write_reads(path: Path, rows: list[ReadRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_READS_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.run_accession,
                    r.biosample,
                    r.bioproject,
                    r.organism,
                    r.taxid,
                    r.platform,
                    r.instrument_model,
                    r.layout,
                    r.bases,
                    r.read_count,
                    r.family,
                    r.genus,
                    r.species,
                    ";".join(r.fastq_urls),
                    ";".join(r.fastq_md5),
                    ";".join(str(b) for b in r.fastq_bytes),
                    r.library_selection,
                ]
            )


def read_reads(path: Path) -> list[ReadRow]:
    rows: list[ReadRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo, delimiter="\t")
        # library_selection is optional: tables written before it existed lack it.
        _require_columns(reader, path, [c for c in _READS_COLUMNS if c != "library_selection"])
        for rec in reader:
            split = lambda s: tuple(x for x in s.split(";") if x)  # noqa: E731
            rows.append(
                ReadRow(
                    run_accession=rec["run_accession"],
                    biosample=rec["biosample"],
                    bioproject=rec["bioproject"],
                    organism=rec["organism"],
                    taxid=rec["taxid"],
                    platform=rec["platform"],
                    instrument_model=rec["instrument_model"],
                    layout=rec["layout"],
                    bases=int(rec["bases"] or 0),
                    read_count=int(rec["read_count"] or 0),
                    family=rec["family"],
                    genus=rec["genus"],
                    species=rec["species"],
                    fastq_urls=split(rec["fastq_urls"]),
                    fastq_md5=split(rec["fastq_md5"]),
                    fastq_bytes=tuple(int(b) for b in split(rec["fastq_bytes"])),
                    # Absent from tables written before the column existed.
                    library_selection=rec.get("library_selection", "") or "",
                )
            )
    return rows


@dataclass(frozen=True)
class AssemblyStatsRow:
    """Per-assembly metrics and labels (assembly_stats.tsv), one row per genome."""

    run_accession: str
    filename: str
    assembler: str
    n_contigs: int
    total_length: int
    n50: int
    largest_contig: int
    est_coverage: float | None = None
    completeness: float | None = None
    contamination: float | None = None
    ncbi_taxonomy: str = ""
    gtdb_taxonomy: str = ""
    label_source: str = "metadata"
    taxonomy_flag: str = ""
    # The polisher that corrected the assembly, empty when none ran.
    polisher: str = ""
    # Whether assemblies/<run>/reads.sig.zip holds a sketch of the run's reads.
    reads_sketch: bool = False
    # "pass" when the reads screen (assemble --screen-reads) passed the run;
    # empty when it was not screened.
    screen: str = ""


_ASSEMBLY_STATS_COLUMNS = [
    "run_accession",
    "filename",
    "assembler",
    "n_contigs",
    "total_length",
    "n50",
    "largest_contig",
    "est_coverage",
    "completeness",
    "contamination",
    "ncbi_taxonomy",
    "gtdb_taxonomy",
    "label_source",
    "taxonomy_flag",
    "polisher",
    "reads_sketch",
    "screen",
]


def write_assembly_stats(path: Path, rows: list[AssemblyStatsRow]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(_ASSEMBLY_STATS_COLUMNS)
        for r in rows:
            writer.writerow(
                [
                    r.run_accession,
                    r.filename,
                    r.assembler,
                    r.n_contigs,
                    r.total_length,
                    r.n50,
                    r.largest_contig,
                    "" if r.est_coverage is None else f"{r.est_coverage:.2f}",
                    "" if r.completeness is None else f"{r.completeness:.2f}",
                    "" if r.contamination is None else f"{r.contamination:.2f}",
                    r.ncbi_taxonomy,
                    r.gtdb_taxonomy,
                    r.label_source,
                    r.taxonomy_flag,
                    r.polisher,
                    1 if r.reads_sketch else 0,
                    r.screen,
                ]
            )


def read_assembly_stats(path: Path) -> list[AssemblyStatsRow]:
    rows: list[AssemblyStatsRow] = []
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo, delimiter="\t"):
            rows.append(
                AssemblyStatsRow(
                    run_accession=rec["run_accession"],
                    filename=rec["filename"],
                    assembler=rec["assembler"],
                    n_contigs=int(rec["n_contigs"]),
                    total_length=int(rec["total_length"]),
                    n50=int(rec["n50"]),
                    largest_contig=int(rec["largest_contig"]),
                    est_coverage=_opt_float(rec["est_coverage"]),
                    completeness=_opt_float(rec["completeness"]),
                    contamination=_opt_float(rec["contamination"]),
                    ncbi_taxonomy=rec["ncbi_taxonomy"],
                    gtdb_taxonomy=rec["gtdb_taxonomy"],
                    label_source=rec["label_source"],
                    taxonomy_flag=rec["taxonomy_flag"],
                    polisher=rec.get("polisher", "") or "",
                    # Absent from tables written before the column existed.
                    reads_sketch=(rec.get("reads_sketch") or "0").strip() == "1",
                    screen=rec.get("screen") or "",
                )
            )
    return rows


@dataclass(frozen=True)
class ExcusedRun:
    """A selected run that produced no genome: which step gave up and why."""

    run_accession: str
    step: str  # fetch | screen | assemble | qc | classify
    reason: str


def write_excused_runs(path: Path, rows: list[ExcusedRun]) -> None:
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["run_accession", "step", "reason"])
        for r in rows:
            # A tool failure's reason carries the output tail on further lines;
            # fold it onto one line so each run stays one TSV row.
            writer.writerow([r.run_accession, r.step, " ".join(r.reason.split())])


def read_excused_runs(path: Path) -> list[ExcusedRun]:
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        return [ExcusedRun(row[0], row[1], row[2]) for row in reader if len(row) >= 3]


def write_segments(
    path: Path,
    members: dict[str, list[str]],
    labels: dict[str, tuple[str, str]] | None = None,
) -> None:
    """Write isolate token -> member accession rows (one row per member).

    ``labels`` maps a member accession to its (normalised segment, segment
    label as submitted); the two columns are empty for an accession without.
    """
    labels = labels or {}
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        writer.writerow(["isolate", "accession", "segment", "segment_label"])
        for token, accessions in members.items():
            for accession in accessions:
                segment, label = labels.get(accession, ("", ""))
                writer.writerow([token, accession, segment, label])


def read_segments(path: Path) -> dict[str, list[str]]:
    """Read ``segments.tsv`` back into isolate token -> member accessions."""
    members: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 2:
                members.setdefault(row[0], []).append(row[1])
    return members


def write_genomes_map(path: Path, mapping: list[tuple[str, str]]) -> None:
    """Write accession -> leaf rows."""
    with atomic_replace(path, newline="") as fo:
        writer = _tsv_writer(fo)
        for accession, leaf in mapping:
            writer.writerow([accession, leaf])
