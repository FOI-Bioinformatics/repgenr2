"""viral metadata stage (BV-BRC + NCBI Entrez).

Ports ``vmetadata.py``: download a virus group's FASTA from the BV-BRC FTP,
derive per-taxid sequence-length statistics from the FASTA headers, enrich the
taxonomy via NCBI Entrez, and write the metadata tables plus a taxname->data
JSON (replacing the old pickle) under ``virus_download_wd/``. The viral genome
stage consumes these files.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from ftplib import FTP, FTP_TLS, all_errors
from pathlib import Path
from statistics import mean, median

from Bio import SeqIO

from ..core.context import WorkdirContext
from ..core.contracts import atomic_path
from ..core.errors import UserInputError, WorkdirError
from ..core.integrity import looks_like_fasta
from ..viral.entrez import TAXNAMES_ORDERED, get_taxon_data_from_entrez

BVBRC_FTP = "ftp.bvbrc.org"
BVBRC_FTP_DIR = "viruses"
# Names the source and target behind virus_download_wd/download.fa, so a
# reuse never serves one target's (or source's) sequences as another's.
DOWNLOAD_SOURCE = "download.source"
# Files only one source writes; the other source removes them so vgenome
# dispatches on the metadata of the latest vmetadata run.
_NCBI_VIRUS_ONLY = ("virus_records.json",)
_BVBRC_ONLY = (
    "metadata_ncbi.tsv",
    "metadata_ncbi_summary.tsv",
    "metadata_ncbi_taxnames_data.json",
)
# The BV-BRC tables copied to the working directory with a virus_ prefix.
_BVBRC_COPIES = ("metadata_base.tsv", "metadata_ncbi.tsv", "metadata_ncbi_summary.tsv")


def _source_key(source: str, target: str) -> str:
    return f"{source}:{target.lower()}"


def _drop(download_wd: Path, names: tuple[str, ...]) -> None:
    for name in names:
        (download_wd / name).unlink(missing_ok=True)


class _ReuseFTP_TLS(FTP_TLS):
    """FTPS client that reuses the control channel's TLS session for data.

    BV-BRC now requires SSL/TLS on the control channel and TLS session reuse on
    the data channel (a common vsftpd ``require_ssl_reuse`` setup); plain FTP and
    a vanilla FTP_TLS both fail ("550 SSL/TLS required" / "425 Unable to build
    data connection").
    """

    def ntransfercmd(self, cmd, rest=None):
        conn, size = FTP.ntransfercmd(self, cmd, rest)
        if self._prot_p:
            conn = self.context.wrap_socket(
                conn, server_hostname=self.host, session=self.sock.session
            )
        return conn, size


def _bvbrc_connect(timeout: int = 120) -> _ReuseFTP_TLS:
    ftp = _ReuseFTP_TLS(BVBRC_FTP, timeout=timeout)
    ftp.login(user="anonymous", passwd="anonymous")
    ftp.prot_p()
    return ftp


@contextmanager
def _bvbrc_session() -> Iterator[_ReuseFTP_TLS]:
    """Open a BV-BRC FTPS session; a connection or protocol failure raises
    :class:`WorkdirError` naming the server, like an HTTP failure in core.http,
    instead of a bare OSError that the CLI would report as unexpected."""
    try:
        with _bvbrc_connect() as ftp:
            yield ftp
    except all_errors as exc:
        raise WorkdirError(f"BV-BRC FTP request failed: {BVBRC_FTP} ({exc})") from exc


@dataclass
class VmetadataParams:
    target: str | None = None
    # BV-BRC header tag; None means the source default ("complete genome").
    # Meaningless for ncbi_virus, which filters with complete_only.
    filter: str | None = None
    list_targets: bool = False
    source: str = "ncbi_virus"  # ncbi_virus | bvbrc
    host: str | None = None  # ncbi_virus: restrict to a host species
    complete_only: bool = False  # ncbi_virus: only COMPLETE sequences
    released_after: str | None = None  # ncbi_virus: MM/DD/YYYY


def run(ctx: WorkdirContext, params: VmetadataParams) -> int:
    logger = ctx.logger
    if params.list_targets:
        _list_targets(logger)
        return 0
    if not params.target:
        raise UserInputError("Supply --target (e.g. adenoviridae) or --list.")

    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True, exist_ok=True)
    if params.source == "ncbi_virus":
        return _run_ncbi_virus(ctx, params, download_wd, logger)
    return _run_bvbrc(ctx, params, download_wd, logger)


def _run_ncbi_virus(ctx, params, download_wd, logger) -> int:
    from ..core.errors import MissingBinaryError
    from ..core.plugins import preflight
    from ..viral import ncbi_virus

    if params.target is None:
        raise UserInputError("Supply --target (e.g. adenoviridae) or --list.")
    records = ncbi_virus.fetch(
        params.target,
        download_wd,
        complete_only=params.complete_only,
        host=params.host,
        released_after=params.released_after,
        logger=logger,
    )
    if not records:
        raise WorkdirError(f"NCBI Virus returned no genomes for '{params.target}'.")
    (download_wd / DOWNLOAD_SOURCE).write_text(
        _source_key("ncbi_virus", params.target) + "\n", encoding="utf-8"
    )
    _drop(download_wd, _BVBRC_ONLY)
    for name in _BVBRC_COPIES[1:]:
        (ctx.workdir / f"virus_{name}").unlink(missing_ok=True)
    ncbi_virus.write_records(download_wd / "virus_records.json", records)
    _write_base_from_records(download_wd / "metadata_base.tsv", records)
    (ctx.workdir / "virus_metadata_base.tsv").write_text(
        (download_wd / "metadata_base.tsv").read_text()
    )

    # Best-effort provenance: fetch already used datasets, so this normally
    # resolves; never fail the stage just to record a version.
    try:
        tool_versions = preflight(ncbi_virus.DATASETS_CAPS)
    except MissingBinaryError:
        tool_versions = {}
    ctx.config.record_stage(
        "vmetadata",
        # The NCBI Virus records come through the NCBI datasets CLI.
        tool="datasets",
        params={
            "source": "ncbi_virus",
            "target": params.target,
            "complete_only": params.complete_only,
            "host": params.host,
            "released_after": params.released_after,
            "sequences": len(records),
            # How many records took their species from NCBI Taxonomy, the
            # report lineage or the organism name.
            "species_source": {
                src: sum(1 for r in records if r.species_source == src)
                for src in sorted({r.species_source for r in records})
            },
        },
        tool_versions=tool_versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Viral metadata (NCBI Virus): %d sequences under %s", len(records), download_wd)
    return len(records)


def _write_base_from_records(path: Path, records) -> None:
    """taxid-keyed length stats (human-readable), grouped from NCBI Virus records."""
    from statistics import mean as _mean
    from statistics import median as _median

    by_taxid: dict[str, dict] = {}
    for r in records:
        g = by_taxid.setdefault(r.taxid, {"lens": [], "species": r.species, "organism": r.organism})
        g["lens"].append(r.length)
    header = ["taxid", "name", "num", "seq_min", "seq_max", "seq_med", "seq_mean", "description"]
    with open(path, "w", encoding="utf-8") as fo:
        fo.write("\t".join(header) + "\n")
        for taxid, g in sorted(by_taxid.items(), key=lambda x: len(x[1]["lens"]), reverse=True):
            lens = g["lens"]
            row = [
                taxid,
                g["species"],
                len(lens),
                min(lens),
                max(lens),
                int(_median(lens)),
                int(_mean(lens)),
                g["organism"],
            ]
            fo.write("\t".join(map(str, row)) + "\n")


def _run_bvbrc(ctx, params, download_wd, logger) -> int:
    target = params.target.lower()
    download_fa = download_wd / "download.fa"
    marker = download_wd / DOWNLOAD_SOURCE
    key = _source_key("bvbrc", target)
    try:
        recorded = marker.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        recorded = None

    # Reuse only the group FASTA this target downloaded: one fetched for
    # another target, or by the NCBI Virus source, is downloaded afresh.
    if looks_like_fasta(download_fa) and recorded == key:
        logger.info("Group FASTA already present; reusing %s", download_fa)
    else:
        if download_fa.exists():
            logger.info(
                "Existing %s was fetched for %s, not %s; downloading again",
                download_fa.name,
                recorded or "an unrecorded target",
                key,
            )
        marker.unlink(missing_ok=True)
        _download_group(target, download_fa, logger)
        marker.write_text(key + "\n", encoding="utf-8")
    _drop(download_wd, _NCBI_VIRUS_ONLY)

    tag = params.filter or "complete genome"
    base, all_taxids, taxid_bvbrc, records = _parse_fasta(download_fa, tag, logger)

    logger.info("Enriching %d taxids via NCBI Entrez", len(all_taxids))
    ncbi_data, missing, _alts = get_taxon_data_from_entrez(all_taxids, logger)
    taxnames_data = _taxnames_data(ncbi_data, missing, taxid_bvbrc)

    _write_base(download_wd / "metadata_base.tsv", base)
    _write_ncbi(download_wd / "metadata_ncbi.tsv", ncbi_data, base)
    _write_ncbi_summary(
        download_wd / "metadata_ncbi_summary.tsv", ncbi_data, missing, taxid_bvbrc, records, base
    )
    serializable = {
        name: {**data, "datasets": sorted(data["datasets"])} for name, data in taxnames_data.items()
    }
    with open(download_wd / "metadata_ncbi_taxnames_data.json", "w", encoding="utf-8") as fo:
        json.dump(serializable, fo)

    # copy human-readable tables to the workdir root with a virus_ prefix
    for name in _BVBRC_COPIES:
        (ctx.workdir / f"virus_{name}").write_text((download_wd / name).read_text())

    ctx.config.record_stage(
        "vmetadata",
        # Group FASTA from the BV-BRC FTP server; no external binary.
        tool="bvbrc",
        params={
            "source": "bvbrc",
            "target": target,
            "filter": tag,
            "taxids": len(all_taxids),
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Viral metadata written under %s", download_wd)
    return len(all_taxids)


def _list_targets(logger) -> None:
    with _bvbrc_session() as ftp:
        ftp.cwd(BVBRC_FTP_DIR)
        targets = sorted(f.replace(".fna", "") for f in ftp.nlst() if f.endswith(".fna"))
    logger.info("Available targets:\n%s", "\n".join(targets))


def resolve_group_name(target: str, listing: list[str]) -> str | None:
    """The BV-BRC group file stem matching ``target``, ignoring case.

    Group files are named like ``Hepatitis_E_virus.fna``; users type the
    target in any case (the CLI lower-cases it), so the listing decides.
    """
    wanted = target.lower()
    for entry in listing:
        stem = entry[:-4] if entry.endswith(".fna") else entry
        if stem.lower() == wanted:
            return stem
    return None


def _download_group(target: str, dest: Path, logger) -> None:
    with _bvbrc_session() as ftp:
        ftp.cwd(BVBRC_FTP_DIR)
        group = resolve_group_name(target, ftp.nlst())
        if group is None:
            raise WorkdirError(
                f"Could not find virus group '{target}' at BV-BRC. "
                "Check the name (try --list), or download all with --target viruses."
            )
        logger.info("Downloading %s from BV-BRC (FTPS)", group)
        remote = f"{group}.fna"  # relative to BVBRC_FTP_DIR after the cwd above
        ftp.sendcmd("TYPE I")
        try:
            remote_size = ftp.size(remote)
        except Exception as exc:
            raise WorkdirError(f"Could not size '{remote}' at BV-BRC: {exc}") from exc
        # Binary transfer with an exact size check: ASCII mode + a 1000-byte
        # tolerance previously let a silently-truncated FASTA pass as complete.
        # The transfer goes to a temporary sibling that replaces ``dest`` only
        # after the size check, so an interrupted run never leaves a partial
        # download.fa for the next run to reuse.
        with atomic_path(dest) as tmp:
            with open(tmp, "wb") as fo:
                ftp.retrbinary("RETR " + remote, fo.write)
            local_size = tmp.stat().st_size
            if remote_size and local_size != remote_size:
                raise WorkdirError(
                    f"Incomplete BV-BRC download: got {local_size} of {remote_size} bytes."
                )
    logger.info("Download finished (%d bytes)", local_size)


def _parse_fasta(path: Path, tag: str, logger):
    """Per-taxid statistics of a BV-BRC group FASTA.

    Returns ``base`` (count and length statistics of the records whose header
    holds ``tag``), every taxid, the BV-BRC genome ids per taxid, and the
    number of records per taxid. A BV-BRC genome id can label several records
    (segments, or contigs of a partial genome), so the last two differ.
    """
    base: dict[str, dict] = {}
    records: dict[str, int] = {}
    all_taxids: set[str] = set()
    taxid_bvbrc: dict[str, set[str]] = {}
    counts: dict[str, int] = {}
    lens: dict[str, list[int]] = {}
    descriptions: dict[str, dict] = {}

    for entry in SeqIO.parse(str(path), "fasta"):
        try:
            bvbrc_id = entry.description.split("| ")[1].split("]")[0]
        except IndexError:
            raise WorkdirError(
                f"{path} is not a BV-BRC group FASTA (header without '[... | id]': "
                f"{entry.description[:80]!r}). Delete it and re-run vmetadata."
            ) from None
        taxid = bvbrc_id.split(".")[0]
        all_taxids.add(taxid)
        taxid_bvbrc.setdefault(taxid, set()).add(bvbrc_id)
        records[taxid] = records.get(taxid, 0) + 1
        if tag not in entry.description:
            continue
        seq_len = len(entry.seq)
        counts[taxid] = counts.get(taxid, 0) + 1
        lens.setdefault(taxid, []).append(seq_len)
        if taxid not in descriptions:
            species = entry.description.split("[")[1].split(" | ")[0]
            full = entry.description.replace(entry.name, "").lstrip().split("[")[0].strip()
            descriptions[taxid] = {"species": species, "full": full}

    for taxid, count in counts.items():
        taxid_lens = lens[taxid]
        base[taxid] = {
            "num": count,
            "min": min(taxid_lens),
            "max": max(taxid_lens),
            "mean": int(mean(taxid_lens)),
            "median": int(median(taxid_lens)),
            "desc": descriptions[taxid],
        }
    logger.info("Found %d taxids; %d passed the '%s' filter", len(all_taxids), len(base), tag)
    return base, all_taxids, taxid_bvbrc, records


def _taxnames_data(ncbi_data, missing, taxid_bvbrc) -> dict[str, dict]:
    taxnames: dict[str, dict] = {}
    for taxid, data in ncbi_data.items():
        if taxid in missing:
            continue
        for level, level_data in data["taxdata"].items():
            name = level_data["name"]
            if name is None:
                continue
            entry = taxnames.setdefault(
                name, {"taxid": level_data["taxid"], "level": level, "datasets": set()}
            )
            for bvbrc_id in taxid_bvbrc.get(taxid, set()):
                entry["datasets"].add(bvbrc_id)
    return taxnames


def _write_base(path: Path, base: dict) -> None:
    header = ["taxid", "name", "num", "seq_min", "seq_max", "seq_med", "seq_mean", "description"]
    with open(path, "w", encoding="utf-8") as fo:
        fo.write("\t".join(header) + "\n")
        for taxid, data in sorted(base.items(), key=lambda x: x[1]["num"], reverse=True):
            row = [
                taxid,
                data["desc"]["species"],
                data["num"],
                data["min"],
                data["max"],
                data["median"],
                data["mean"],
                data["desc"]["full"],
            ]
            fo.write("\t".join(map(str, row)) + "\n")


def _write_ncbi_summary(
    path: Path,
    ncbi_data: dict,
    missing: set[str],
    taxid_bvbrc: dict[str, set[str]],
    records: dict[str, int],
    base: dict,
) -> None:
    """Genomes and sequences per NCBI taxon at each rank (RepGenR 1's
    metadata_ncbi_summary).

    One row per rank and taxon name: the taxon's NCBI taxid, the taxids of the
    group FASTA below it (``taxids``), their BV-BRC genome ids (``genomes``),
    their FASTA records (``sequences``) and the records whose header holds the
    filter tag (``sequences_with_tag``). Ordered by rank as in
    metadata_ncbi.tsv, then by sequences with the tag, then by name.
    """
    totals: dict[tuple[str, str], dict] = {}
    for taxid, data in ncbi_data.items():
        if taxid in missing:
            continue
        for level in TAXNAMES_ORDERED:
            level_data = data["taxdata"].get(level) or {}
            name = level_data.get("name")
            if name is None:
                continue
            entry = totals.setdefault(
                (level, name),
                {
                    "taxid": level_data.get("taxid"),
                    "taxids": 0,
                    "genomes": 0,
                    "seqs": 0,
                    "tagged": 0,
                },
            )
            entry["taxids"] += 1
            entry["genomes"] += len(taxid_bvbrc.get(taxid, ()))
            entry["seqs"] += records.get(taxid, 0)
            entry["tagged"] += base[taxid]["num"] if taxid in base else 0
    rank = {level: i for i, level in enumerate(TAXNAMES_ORDERED)}
    ordered = sorted(totals.items(), key=lambda kv: (rank[kv[0][0]], -kv[1]["tagged"], kv[0][1]))
    with open(path, "w", encoding="utf-8") as fo:
        fo.write("level\tname\ttaxid\ttaxids\tgenomes\tsequences\tsequences_with_tag\n")
        for (level, name), e in ordered:
            row = [
                level,
                name,
                e["taxid"] or "",
                e["taxids"],
                e["genomes"],
                e["seqs"],
                e["tagged"],
            ]
            fo.write("\t".join(map(str, row)) + "\n")


def _write_ncbi(path: Path, ncbi_data: dict, base: dict) -> None:
    header = (
        ["taxid", "name", "num_with_tag"]
        + TAXNAMES_ORDERED
        + [f"{x}_taxid" for x in TAXNAMES_ORDERED]
    )
    with open(path, "w", encoding="utf-8") as fo:
        fo.write("\t".join(header) + "\n")
        ordered = sorted(
            ncbi_data.items(),
            key=lambda x: base[x[0]]["num"] if x[0] in base else 0,
            reverse=True,
        )
        for taxid, data in ordered:
            num = base[taxid]["num"] if taxid in base else None
            row = [taxid, data["name"], num]
            row += [data["taxdata"][name]["name"] for name in TAXNAMES_ORDERED]
            row += [data["taxdata"][name]["taxid"] for name in TAXNAMES_ORDERED]
            fo.write("\t".join("" if v is None else str(v) for v in row) + "\n")
