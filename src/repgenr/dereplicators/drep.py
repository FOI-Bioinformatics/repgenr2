"""dRep dereplication adapter.

Ports the logic of the old ``derep_worker.py``: optionally decompress gzipped
genomes (prodigal inside dRep needs plain FASTA), run ``dRep dereplicate``, then
parse ``Cdb.csv`` (cluster membership) and ``genomeInformation.csv`` (QC) to
build a normalized :class:`DerepResult`.

dRep does not scale natively to very large sets; the dereplicate stage wraps it
with chunking when needed (``supports_native_scaling=False``).
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import list_fasta
from ..core.errors import ToolExecutionError, WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import link_or_copy, write_fofn
from .base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    CompareResult,
    Dereplicator,
    DerepParams,
    DerepResult,
    write_genome_info,
)


class DrepDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="drep",
        container="quay.io/biocontainers/drep:3.7.1--pyhdfd78af_0",
        conda=("bioconda::drep",),
        accepted_extras=frozenset({"virus", "S_algorithm", "length"}),
        # dRep's overlap gate is --cov_thresh with its own default; the shared
        # aligned-fraction setting is not translated onto it.
        ignored_params=frozenset({"aligned_fraction"}),
        required_binaries=(BinarySpec("dRep", version_args=("-h",), min_version="3.0"),),
        # No default_params: the stage merges them into the extras, which
        # would hide whether the user chose S_algorithm, and virus mode's own
        # default (ANImf) would never apply. dereplicate() holds the defaults.
        recommended_max_genomes=2000,
        supports_native_scaling=False,
    )

    def compare(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        threads: int,
        logger: logging.Logger,
        *,
        sketches: Mapping[Path, Path] | None = None,
    ) -> CompareResult:
        """All-vs-all Mash comparison via ``dRep compare`` (primary ANI only).

        dRep sketches with Mash, so the sourmash ``sketches`` are not read.
        """
        # Genome paths go through a fofn (dRep accepts a paths file for -g), so
        # a large set cannot overflow argv; declare the genome dirs for the
        # container backend since the paths are no longer argv tokens.
        out_dir.mkdir(parents=True, exist_ok=True)
        fofn = write_fofn(genomes, out_dir / "genomes.fofn")
        genome_dirs = sorted({str(g.parent) for g in genomes})
        compare_wd = out_dir / "drep_compare_wd"
        run_tool(
            self.capabilities,
            [
                "dRep",
                "compare",
                "--SkipSecondary",
                "-g",
                fofn,
                "--processors",
                str(threads),
                compare_wd,
            ],
            logger=logger,
            log_prefix="drep",
            extra_mounts=genome_dirs,
        )
        dendrogram = compare_wd / "figures" / "Primary_clustering_dendrogram.pdf"
        mdb = compare_wd / "data_tables" / "Mdb.csv"
        return CompareResult(
            similarity_csv=mdb if mdb.exists() else None,
            dendrogram=dendrogram if dendrogram.exists() else None,
            measure="MASH ANI",
        )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        genomes_dir = out_dir / "genomes"
        genomes_dir.mkdir(exist_ok=True)

        # Stage plain (decompressed) copies so dRep/prodigal can read them.
        staged: list[Path] = []
        for src in genomes:
            staged.append(_stage_genome(src, genomes_dir))

        drep_wd = out_dir / "drep_workdir"
        # dRep accepts a text file of genome paths for -g; a large set on argv
        # would hit the OS ARG_MAX limit.
        fofn = write_fofn(staged, out_dir / "genomes.fofn")
        virus = bool(params.extra.get("virus"))
        # Virus mode changes the default secondary algorithm; an explicit
        # S_algorithm tool-arg wins either way, and the flag is passed once.
        s_algorithm = params.extra.get("S_algorithm", "ANImf" if virus else "fastANI")
        cmd: list[str | Path] = [
            "dRep",
            "dereplicate",
            drep_wd,
            "-g",
            fofn,
            "--processors",
            str(params.threads),
            "-sa",
            str(params.secondary_ani),
            "-pa",
            str(params.primary_ani),
            "--S_algorithm",
            s_algorithm,
            "--length",
            str(params.extra.get("length", 0)),
        ]
        genome_info = (
            None if virus else _write_genome_info(genomes, staged, params.quality, out_dir, logger)
        )
        if genome_info is not None:
            cmd += ["--genomeInfo", genome_info]
        if virus:
            cmd += [
                "--cov_thresh",
                "0.5",
                "--N50_weight",
                "0",
                "--size_weight",
                "1",
                "--ignoreGenomeQuality",
                "--clusterAlg",
                "single",
            ]
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="drep",
            extra_mounts=[str(genomes_dir), *([str(genome_info)] if genome_info else [])],
        )

        if not drep_wd.exists():
            raise WorkdirError(
                "dRep working directory was not created; confirm dRep is installed and runs."
            )
        if not (drep_wd / "data_tables" / "Cdb.csv").exists():
            # dRep logs some fatal problems and still exits 0; the most common
            # is CheckM missing when no genome quality is given.
            raise ToolExecutionError(
                [str(c) for c in cmd],
                0,
                "dRep wrote no data_tables/Cdb.csv. dRep stops this way when CheckM "
                "is not on PATH: it scores genome quality with CheckM unless every "
                "genome has completeness and contamination in the manifest (passed "
                "as --genomeInfo) or --ignoreGenomeQuality is set, which --virus "
                "does. The [drep] lines in repgenr.log show the step where it stopped.",
                tool="dRep",
            )
        result = _parse_drep_output(drep_wd, logger)
        result = _restore_input_names(
            result, {st.name: src for st, src in zip(staged, genomes, strict=True)}
        )
        # dRep's filter step (--length, CheckM completeness/contamination)
        # drops genomes from every table it writes; those are QC failures.
        dropped = sorted(g.name for g in genomes if g.name not in result.genome_status)
        if dropped:
            logger.warning(
                "dRep filtered out %d genome(s) before clustering (length or genome "
                "quality); marked fail_qc in genome_status.tsv (e.g. %s).",
                len(dropped),
                ", ".join(dropped[:3]),
            )
            for name in dropped:
                result.genome_status[name] = STATUS_FAIL_QC
        result.qc_reasons = _qc_reasons(
            result.genome_status,
            list(zip(genomes, staged, strict=True)),
            _genome_info_rows(drep_wd),
            params.quality,
            min_length=int(params.extra.get("length", 0) or 0),
            check_quality=not virus,
        )
        return result


def _write_genome_info(
    genomes: Sequence[Path],
    staged: Sequence[Path],
    quality: Mapping[str, tuple[float, float]],
    out_dir: Path,
    logger: logging.Logger,
) -> Path | None:
    """Write dRep's ``--genomeInfo`` table from the manifest quality, if complete.

    With this table dRep uses the given completeness and contamination instead
    of running CheckM, so it runs without CheckM installed and scores genomes
    with the same values as ``--keeper quality``. dRep requires a row for every
    genome, so the table is written only when every genome has both values.
    Rows name the staged file dRep reads (a gzipped input is staged
    decompressed, without ``.gz``).
    """
    lacking = [g.name for g in genomes if g.name not in quality]
    if lacking:
        if len(lacking) < len(genomes):
            logger.warning(
                "dRep is not given the manifest quality (--genomeInfo): %d of %d genome(s) "
                "have no completeness and contamination (e.g. %s), and dRep needs values for "
                "every genome. dRep scores genome quality with CheckM instead.",
                len(lacking),
                len(genomes),
                ", ".join(lacking[:3]),
            )
        else:
            logger.info(
                "No genome quality in the manifest; dRep scores genome quality with CheckM."
            )
        return None
    path = write_genome_info(
        out_dir / "genome_info.csv",
        ((st.name, *quality[src.name]) for src, st in zip(genomes, staged, strict=True)),
    )
    logger.info(
        "dRep uses the manifest completeness and contamination of %d genomes (--genomeInfo)",
        len(genomes),
    )
    return path


def _restore_input_names(result: DerepResult, source_by_staged: dict[str, Path]) -> DerepResult:
    """Report genomes under their input names, not dRep's staged copies.

    A gzipped input is staged decompressed (``x.fasta.gz`` -> ``x.fasta``), so
    dRep names it without ``.gz``; the stage checks every input name for a
    status and links representatives from genomes/, so map the names back.
    """

    def name(staged: str) -> str:
        src = source_by_staged.get(staged)
        return src.name if src is not None else staged

    return DerepResult(
        representatives=[source_by_staged.get(p.name, p) for p in result.representatives],
        clusters={
            name(rep): [name(m) for m in members] for rep, members in result.clusters.items()
        },
        genome_status={name(g): state for g, state in result.genome_status.items()},
        genome_information=[
            {**row, "genome": name(str(row.get("genome", "")))}
            for row in (result.genome_information or [])
        ]
        or result.genome_information,
    )


def _stage_genome(src: Path, dest_dir: Path) -> Path:
    """Stage ``src`` into ``dest_dir`` for dRep, decompressing a .gz on the way.

    A plain FASTA is hardlinked (instant, no extra disk -- dRep only reads it);
    only a .gz genuinely needs a new decompressed file. Staging the whole input
    set previously copied every genome.
    """
    if src.suffix == ".gz":
        import gzip

        target = dest_dir / src.with_suffix("").name
        with gzip.open(src, "rb") as fi, open(target, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return target
    target = dest_dir / src.name
    if not target.exists():
        link_or_copy(src, target)
    return target


def _parse_drep_output(drep_wd: Path, logger: logging.Logger) -> DerepResult:
    derep_genomes = drep_wd / "dereplicated_genomes"
    representatives = list_fasta(derep_genomes)
    rep_names = {p.name for p in representatives}

    # Cluster membership from Cdb.csv (genome, secondary_cluster).
    clusters_by_id: dict[str, set[str]] = {}
    genome_cluster: dict[str, str] = {}
    cdb = drep_wd / "data_tables" / "Cdb.csv"
    with open(cdb, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            genome = row["genome"]
            clust = row.get("secondary_cluster") or row.get("cluster") or ""
            genome_cluster[genome] = clust
            clusters_by_id.setdefault(clust, set()).add(genome)

    status: dict[str, str] = {}
    clusters: dict[str, list[str]] = {}
    for rep in representatives:
        status[rep.name] = STATUS_REPRESENTATIVE
        clust = genome_cluster.get(rep.name)
        members = clusters_by_id.get(clust, set()) if clust is not None else set()
        contained = sorted(m for m in members if m != rep.name)
        clusters[rep.name] = contained
        for m in contained:
            status[m] = STATUS_CONTAINED

    # Genomes filtered out by dRep QC have centrality 0 in genomeInformation.csv.
    info = _parse_genome_information(drep_wd, logger)
    for genome, centrality in info.items():
        if centrality == 0.0 and genome not in rep_names:
            status.setdefault(genome, STATUS_FAIL_QC)

    return DerepResult(
        representatives=representatives,
        clusters=clusters,
        genome_status=status,
        genome_information=[{"genome": g, "centrality": c} for g, c in info.items()],
    )


# dRep's default filter (dRep dereplicate --completeness / --contamination);
# RepGenR2 passes --length (0 unless the length tool-arg is given) and leaves
# these two at dRep's defaults.
DREP_MIN_COMPLETENESS = 75.0
DREP_MAX_CONTAMINATION = 25.0


def _qc_reasons(
    status: Mapping[str, str],
    inputs: Sequence[tuple[Path, Path]],
    info: Mapping[str, dict[str, float | None]],
    quality: Mapping[str, tuple[float, float]],
    *,
    min_length: int,
    check_quality: bool,
) -> dict[str, str]:
    """Name the dRep filter criteria each ``fail_qc`` genome fails.

    ``inputs`` pairs each input genome with the staged copy dRep read. The
    values come from dRep's genomeInformation.csv (staged names), else from
    the manifest quality and the staged FASTA. A genome no criterion explains
    is given ``drep_filter``.
    """
    reasons: dict[str, str] = {}
    for src, st in inputs:
        if status.get(src.name) != STATUS_FAIL_QC:
            continue
        row = info.get(st.name, {})
        completeness = row.get("completeness")
        contamination = row.get("contamination")
        if completeness is None and src.name in quality:
            completeness, contamination = quality[src.name]
        failed: list[str] = []
        if min_length > 0:
            length = row.get("length")
            if length is None:
                length = _sequence_length(st)
            if length < min_length:
                failed.append("length")
        if check_quality:
            if completeness is not None and completeness < DREP_MIN_COMPLETENESS:
                failed.append("completeness")
            if contamination is not None and contamination > DREP_MAX_CONTAMINATION:
                failed.append("contamination")
        reasons[src.name] = ",".join(failed) or "drep_filter"
    return reasons


def _sequence_length(path: Path) -> int:
    """Residues in a plain FASTA (the staged copies are never gzipped)."""
    total = 0
    with open(path, encoding="utf-8", errors="replace") as fo:
        for line in fo:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def _genome_info_rows(drep_wd: Path) -> dict[str, dict[str, float | None]]:
    """genomeInformation.csv as genome -> completeness, contamination, length.

    The 7-column layout holds all three; the 5-column layout of
    ``--ignoreGenomeQuality`` holds the length only. Unreadable values are None.
    """
    path = drep_wd / "data_tables" / "genomeInformation.csv"
    if not path.exists():
        return {}

    def num(value: str) -> float | None:
        try:
            return float(value)
        except ValueError:
            return None

    out: dict[str, dict[str, float | None]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo)
        for row in reader:
            genome = row.get("genome")
            if not genome:
                continue
            out[os.path.basename(genome)] = {
                key: num(row[key]) if row.get(key) not in (None, "") else None
                for key in ("completeness", "contamination", "length")
            }
    return out


def _parse_genome_information(drep_wd: Path, logger: logging.Logger) -> dict[str, float]:
    """Return genome -> centrality, normalizing the viral 5-column variant.

    dRep run with ``--ignoreGenomeQuality`` writes 5 columns
    (genome_path,length,N50,genome,centrality) instead of the usual 7; we read
    the centrality field in either layout.
    """
    path = drep_wd / "data_tables" / "genomeInformation.csv"
    if not path.exists():
        logger.warning("genomeInformation.csv not found; skipping QC status")
        return {}
    out: dict[str, float] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo)
        header = next(reader, None)
        if header is None:
            return out
        ncols = len(header)
        for row in reader:
            if ncols == 7:
                genome, centrality = row[0], row[6]
            elif ncols == 5:
                genome, centrality = row[3], row[4]
            else:
                continue
            try:
                out[genome] = float(centrality)
            except ValueError:
                out[genome] = 0.0
    return out
