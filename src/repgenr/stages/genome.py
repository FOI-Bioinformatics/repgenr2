"""genome stage: download and organize genomes selected by the metadata stage.

Ports ``genome.py``: read selected accessions from the manifest, download them
with NCBI ``datasets`` (dehydrated zip -> rehydrate), and store each as
``genomes/{family}_{genus}_{species}_{accession}.fasta``. The outgroup is
downloaded into ``outgroup/``. Filenames are written back to the manifest.
"""

from __future__ import annotations

import shutil
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlsplit

from ..core import process
from ..core.binaries import BinarySpec
from ..core.containers import run_tool_with_retries
from ..core.context import WorkdirContext
from ..core.contracts import (
    FASTA_SUFFIXES,
    MISSING_ACCESSIONS_TXT,
    atomic_path,
    genome_filename,
    list_fasta,
)
from ..core.errors import ToolExecutionError, WorkdirError
from ..core.http import NCBI_DATASETS_URL, require_reachable, verify_md5
from ..core.integrity import looks_like_fasta
from ..core.outgroup_choice import check_outgroup_distance
from ..core.plugins import ToolCapabilities, preflight
from ..core.process import check_free_disk, remove_tree
from ..core.sketches import require_sourmash_if_requested, sketch_stage_genomes, sketches_dir

_DATASETS = BinarySpec("datasets", version_args=("--version",))
DATASETS_CAPS = ToolCapabilities(
    name="datasets",
    required_binaries=(_DATASETS,),
    conda=("conda-forge::ncbi-datasets-cli",),
)
_DOWNLOAD_BATCH_SIZE = 5000  # accessions per datasets download/rehydrate call
_EST_BYTES_PER_GENOME = 5_000_000  # rough peak-disk estimate (zip + extract + final)


_DATASETS_TIMEOUT = 3600.0  # floor; a hung transfer must not block forever
_PER_ACCESSION_TIMEOUT = 3.0  # a full 5000-accession batch gets ~4 h on a slow link


def _timeout_for(n_items: int) -> float:
    """Per-invocation timeout scaled to the batch size (floor for small calls)."""
    return max(_DATASETS_TIMEOUT, n_items * _PER_ACCESSION_TIMEOUT)


def _run_cmd(cmd, n_items: int = 0, **kwargs):
    """Run the datasets CLI with retries and a timeout (it has neither built in)."""
    return run_tool_with_retries(DATASETS_CAPS, cmd, timeout=_timeout_for(n_items), **kwargs)


def no_assembly_matched(exc: ToolExecutionError) -> bool:
    """datasets found none of the requested accessions (exit 1, no package).

    NCBI suppresses or withdraws assemblies that a GTDB release still lists;
    a batch made only of those fails this way. Retrying cannot help.
    """
    return "no genome assemblies that match" in (exc.output or "").lower()


@dataclass
class GenomeParams:
    accession_list_only: bool = False
    keep_files: bool = False
    # Genome sketches (core.sketches): None sketches when sourmash can run,
    # True requires it, False skips them.
    sketch: bool | None = None


def run(ctx: WorkdirContext, params: GenomeParams) -> int:
    logger = ctx.logger
    versions = preflight(DATASETS_CAPS)
    if not params.accession_list_only:
        require_sourmash_if_requested(params.sketch)

    manifest = ctx.manifest
    selected = [g for g in manifest.all_genomes(include_outgroup=False)]
    outgroup = [g for g in manifest.all_genomes(include_outgroup=True) if g.is_outgroup]
    if not selected:
        raise WorkdirError("No selected genomes in manifest. Run the metadata stage first.")

    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    filenames = {g.accession: _output_name(g) for g in selected}

    # Drop genome files that are no longer selected.
    _prune(ctx.genomes_dir, set(filenames.values()), logger)

    # A present file only counts when it plausibly holds FASTA data: an HTML
    # error page left behind by a crashed run must be re-downloaded, not shipped.
    to_download = [
        acc for acc, name in filenames.items() if not looks_like_fasta(ctx.genomes_dir / name)
    ]
    logger.info(
        "%d to download, %d already present", len(to_download), len(selected) - len(to_download)
    )

    acc_list = ctx.workdir / "ncbi_acc_download_list.txt"
    acc_list.write_text("".join(f"{acc}\n" for acc in to_download))
    if params.accession_list_only:
        logger.info("Accession list written; stopping (--accession-list-only)")
        return 0

    missing: list[str] = []
    if to_download:
        missing = download_accessions(
            to_download,
            filenames,
            ctx.genomes_dir,
            ctx.scratch_dir / "genome_download",
            logger,
            params.keep_files,
        )
    # Record accessions NCBI returned nothing for, so downstream completeness
    # checks can excuse them (and a re-run retries them).
    (ctx.workdir / MISSING_ACCESSIONS_TXT).write_text(
        "".join(f"{acc}\n" for acc in sorted(missing)), encoding="utf-8"
    )

    if outgroup:
        _download_outgroup(ctx, outgroup[0], logger)
    elif ctx.outgroup_dir.is_dir():
        # A selection without an outgroup (metadata --no-outgroup): an
        # outgroup left by an earlier selection would still reach the tree.
        for path in list_fasta(ctx.outgroup_dir):
            path.unlink()
            logger.info("Removed the outgroup %s of an earlier selection", path.name)

    # record filenames back into the manifest (one batched transaction)
    for g in selected:
        g.filename = filenames[g.accession]
    manifest.upsert_many(list(selected))

    sketches, sketch_versions = sketch_stage_genomes(ctx, params.sketch, "genome", logger)
    outgroups = list_fasta(ctx.outgroup_dir) if ctx.outgroup_dir.is_dir() else []
    if outgroups and sketches_dir(ctx.workdir).is_dir():
        check_outgroup_distance(
            list_fasta(ctx.genomes_dir), outgroups, sketches_dir(ctx.workdir), logger
        )
    ctx.config.record_stage(
        "genome",
        tool="datasets",
        params={"downloaded": len(to_download), "total": len(selected), "sketches": sketches},
        tool_versions={**versions, **sketch_versions},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(selected)


def _output_name(g) -> str:
    return genome_filename(g.family, g.genus, g.species, g.accession)


def _prune(genomes_dir: Path, keep: set[str], logger) -> None:
    # Only ever remove genome FASTA files (never directories or user-placed
    # files): a non-selected FASTA is stale, anything else is left untouched.
    for f in genomes_dir.iterdir():
        if f.name in keep or not f.is_file() or not f.name.endswith(FASTA_SUFFIXES):
            continue
        f.unlink()
        logger.info("Removed %s (no longer selected)", f.name)


def _check_disk(scratch_dir: Path, n_accessions: int, logger) -> None:
    """Refuse to start a download with almost no free disk; warn when tight.

    Genome sizes are unknown ahead of time, so the estimate is coarse: a guard
    against the "filled the volume mid-run" failure, not a reservation.
    """
    check_free_disk(
        scratch_dir,
        n_accessions * _EST_BYTES_PER_GENOME * 2,
        logger,
        what=f"download {n_accessions} genomes",
    )


def _assert_fasta(path: Path) -> None:
    """Cheap sanity check that a downloaded file is FASTA, not an error page.

    NCBI/GTDB occasionally serve an HTML error body with a 200 status; catching
    it here names the offending file instead of failing cryptically inside an
    aligner far downstream.
    """
    with open(path, "rb") as fo:
        head = fo.read(4096)
    if not head.lstrip().startswith(b">"):
        raise WorkdirError(
            f"Downloaded genome is not FASTA (does not start with '>'): {path.name}. "
            "The download may be an error page; re-run to retry."
        )


def download_accessions(
    accessions: list[str],
    filenames: dict[str, str],
    dest_dir: Path,
    scratch_dir: Path,
    logger,
    keep_files: bool = False,
) -> list[str]:
    """Download genomes in fixed-size sub-batches into ``dest_dir``.

    Returns the accessions no genome came back for.

    Each sub-batch is downloaded, rehydrated, and moved into ``dest_dir``
    independently, so a failure loses only one batch and a re-run resumes (the
    caller recomputes the still-missing accessions). This also bounds the size of
    any single dehydrated zip / rehydrate at 1000s-100000s of accessions.
    ``filenames`` maps accession -> output filename; scratch (zips, extracts)
    lives under ``scratch_dir``. ctx-free so both the genome stage and the
    stateless ``genome-fetch`` step reuse it.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    total = len(accessions)
    _check_disk(scratch_dir, total, logger)
    if accessions:
        require_reachable(NCBI_DATASETS_URL, what="NCBI datasets")
    n_batches = (total + _DOWNLOAD_BATCH_SIZE - 1) // _DOWNLOAD_BATCH_SIZE
    missing: list[str] = []
    for bi, start in enumerate(range(0, total, _DOWNLOAD_BATCH_SIZE)):
        batch = accessions[start : start + _DOWNLOAD_BATCH_SIZE]
        logger.info("Download batch %d/%d (%d accessions)", bi + 1, n_batches, len(batch))
        missing += _download_one_batch(
            batch, filenames, dest_dir, scratch_dir, logger, keep_files, bi
        )
    return missing


def _download_one_batch(
    batch: list[str],
    filenames: dict[str, str],
    dest_dir: Path,
    scratch_dir: Path,
    logger,
    keep_files: bool,
    bi: int,
) -> list[str]:
    acc_file = scratch_dir / f"ncbi_acc_batch{bi}.txt"
    acc_file.write_text("\n".join(batch))
    zip_path = scratch_dir / f"ncbi_download_{bi}.zip"
    extract = scratch_dir / f"ncbi_extract_{bi}"
    if extract.exists():
        remove_tree(extract)

    try:
        _run_cmd(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                "--dehydrated",
                "--inputfile",
                acc_file,
                "--filename",
                zip_path,
            ],
            n_items=len(batch),
            logger=logger,
            log_prefix="datasets",
            permanent=no_assembly_matched,
        )
    except ToolExecutionError as exc:
        if not no_assembly_matched(exc):
            raise
        logger.warning(
            "NCBI serves none of the %d accessions in batch %d (e.g. %s); recorded as missing",
            len(batch),
            bi,
            ", ".join(batch[:3]),
        )
        acc_file.unlink(missing_ok=True)
        return list(batch)
    process.unzip(zip_path, extract)
    # Rehydrate fetches the files from the URLs in the package's fetch.txt;
    # their hosts depend on the datasets version (api.ncbi.nlm.nih.gov for
    # 18.x), so the hosts named there are probed, not a fixed one.
    for url in _fetch_hosts(extract):
        require_reachable(url, what="datasets rehydrate")
    _run_cmd(
        ["datasets", "rehydrate", "--directory", extract],
        n_items=len(batch),
        logger=logger,
        log_prefix="datasets",
    )

    checksums = _package_checksums(extract)
    produced: set[str] = set()
    for fna in extract.rglob("*.fna"):
        # Skip AppleDouble "._x.fna" side files that macOS writes next to
        # every file on non-native volumes; they are not FASTA and would
        # be reported as a discarded download for the same accession.
        if fna.name.startswith("."):
            continue
        name = filenames.get(fna.parent.name)
        if name:
            # Validate BEFORE moving into genomes/: a non-FASTA body (HTML error
            # page served with HTTP 200) must never land at a genome path, where
            # a re-run would accept it as already present.
            if not looks_like_fasta(fna):
                logger.warning(
                    "Discarding non-FASTA download for %s (error page?)", fna.parent.name
                )
                fna.unlink(missing_ok=True)
                continue
            # rehydrate does not check the package's md5sum.txt; a corrupt
            # transfer that still starts with '>' would otherwise be kept.
            expected = checksums.get(fna.relative_to(extract).as_posix())
            if expected is not None:
                try:
                    verify_md5(fna, expected)
                except WorkdirError as exc:
                    logger.warning("Discarding download for %s: %s", fna.parent.name, exc)
                    fna.unlink(missing_ok=True)
                    continue
            target = dest_dir / name
            shutil.move(str(fna), str(target))
            produced.add(fna.parent.name)

    missing = [acc for acc in batch if acc not in produced]
    if missing:
        logger.warning(
            "NCBI returned no genome for %d of %d accessions in batch %d (e.g. %s)",
            len(missing),
            len(batch),
            bi,
            ", ".join(missing[:3]),
        )

    if not keep_files:
        zip_path.unlink(missing_ok=True)
        acc_file.unlink(missing_ok=True)
        shutil.rmtree(extract, ignore_errors=True)
    return missing


def _fetch_hosts(extract: Path) -> list[str]:
    """One URL per distinct host in a dehydrated package's ``fetch.txt``."""
    fetch = extract / "ncbi_dataset" / "fetch.txt"
    if not fetch.is_file():
        return []
    by_host: dict[str, str] = {}
    for line in fetch.read_text(encoding="utf-8", errors="replace").splitlines():
        url = line.split("\t", 1)[0].strip()
        parts = urlsplit(url)
        if parts.scheme in ("http", "https") and parts.hostname:
            by_host.setdefault(parts.hostname, f"{parts.scheme}://{parts.hostname}/")
    return [by_host[h] for h in sorted(by_host)]


def _package_checksums(extract: Path) -> dict[str, str]:
    """``md5sum.txt`` of a datasets package: path relative to the package -> md5."""
    md5_file = extract / "md5sum.txt"
    if not md5_file.is_file():
        return {}
    out: dict[str, str] = {}
    for line in md5_file.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2 and len(parts[0]) == 32:
            out[parts[1].strip().removeprefix("./")] = parts[0].lower()
    return out


def _download_outgroup(ctx, outgroup, logger) -> None:
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    name = _output_name(outgroup)
    out_path = ctx.outgroup_dir / name
    if looks_like_fasta(out_path):
        logger.info("Outgroup %s already present", outgroup.accession)
    else:
        _fetch_outgroup(ctx, outgroup, out_path, logger)
    # Only the current outgroup may remain: stale outgroups from earlier
    # selections would otherwise accumulate and could shadow this one.
    _prune(ctx.outgroup_dir, {name}, ctx.logger)


def _fetch_outgroup(ctx, outgroup, out_path: Path, logger) -> None:
    zip_path = ctx.workdir / "ncbi_download_outgroup.zip"
    hint = "choose another with 'repgenr metadata --outgroup-accession'."
    require_reachable(NCBI_DATASETS_URL, what="NCBI datasets")
    try:
        _run_cmd(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                outgroup.accession,
                "--filename",
                zip_path,
            ],
            logger=logger,
            log_prefix="datasets",
            permanent=no_assembly_matched,
        )
    except ToolExecutionError as exc:
        if not no_assembly_matched(exc):
            raise
        raise WorkdirError(
            f"NCBI serves no assembly for the outgroup {outgroup.accession}; {hint}"
        ) from exc
    found = False
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for item in zf.namelist():
                if item.endswith(".fna") and outgroup.accession in item:
                    # Written beside the target and renamed, so an interrupted
                    # run never leaves a partial outgroup that looks complete.
                    with atomic_path(out_path) as tmp:
                        tmp.write_bytes(zf.read(item))
                    found = True
                    break
    except zipfile.BadZipFile as exc:
        raise WorkdirError(
            f"Corrupt outgroup download: {zip_path} ({exc}). Re-run to retry."
        ) from exc
    if not found:
        raise WorkdirError(
            f"The outgroup download for {outgroup.accession} holds no genome FASTA; {hint}"
        )
    _assert_fasta(out_path)
    zip_path.unlink(missing_ok=True)
