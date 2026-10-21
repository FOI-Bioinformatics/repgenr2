"""The marker database: GTDB's single-copy marker families as HMMER3 models.

The marker identifiers are GTDB-Tk's (``gtdbtk/config/common.py``,
``BAC120_MARKERS`` and ``AR53_MARKERS``) and ship with the package in
``markers_data/``. The models themselves are fetched one family at a time
from the InterPro API, which serves each as gzip-compressed HMMER3 text with
the family's own GA, TC and NC cutoff lines. The versions are InterPro's
current ones, so they can be newer than the versions GTDB-Tk pins (for
example ``PF00380.26`` against GTDB-Tk's ``PF00380.20``); each one is recorded.

GTDB-Tk searches the Pfam families with the gathering cutoff (GA) and the
TIGRFAM families with the noise cutoff (NC, ``hmmsearch --cut_nc``); the index
records which cutoff applies to each family, so a reader does not need to know.

A database directory holds::

    hmms/<set>.hmm   one concatenated HMMER3 file per marker set
    markers.tsv      set, accession, source, cutoff, name, version, length
    version.json     {"v": 1, "sets", "gtdbtk_lists", "fetched", "interpro_api", "n_models"}

The directory is written in full beside its destination and renamed into
place, so an interrupted download leaves no partial database.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import shutil
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path

from ..core.errors import UserInputError, WorkdirError

MARKER_SETS = ("bac120", "ar53")
MARKERS_DB_ENV = "REPGENR_MARKERS_DB"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/wwwapi/entry"
# The commit date of gtdbtk/config/common.py the packaged lists were taken from
# (Ecogenomics/GTDBTk 6f025f0665).
GTDBTK_LISTS = "2026-10-07"
SCHEMA_VERSION = 1
HMM_DIR = "hmms"
MARKERS_TSV = "markers.tsv"
VERSION_JSON = "version.json"
INDEX_COLUMNS = ("set", "accession", "source", "cutoff", "name", "version", "length")
DOWNLOAD_HINT = "run `repgenr markers-db download --output DIR` (168 families, about 25 MB)"

# url -> response body. The default goes through the shared HTTP session.
Fetcher = Callable[[str], bytes]


@dataclass(frozen=True)
class MarkerModel:
    """One family's HMMER3 model and the fields the index records."""

    accession: str  # without the version suffix, as in the marker lists
    source: str  # pfam | tigrfam
    cutoff: str  # ga | nc
    name: str
    version: str  # the full ACC value, e.g. PF00380.26
    length: int
    text: str  # the model, ending with the "//" line


def load_marker_list(name: str) -> list[str]:
    """The accessions of marker set ``name`` (bac120 or ar53), in GTDB-Tk's order."""
    if name not in MARKER_SETS:
        raise UserInputError(f"Unknown marker set '{name}'; choose from {', '.join(MARKER_SETS)}.")
    text = (resources.files(__package__) / "markers_data" / f"{name}.txt").read_text("utf-8")
    return [s for line in text.splitlines() if (s := line.strip()) and not s.startswith("#")]


def marker_source(accession: str) -> str:
    """``pfam`` or ``tigrfam`` from the accession prefix."""
    if accession.startswith("PF"):
        return "pfam"
    if accession.startswith("TIGR"):
        return "tigrfam"
    raise UserInputError(f"Marker accession '{accession}' is neither Pfam (PF) nor TIGRFAM (TIGR).")


def marker_cutoff(source: str) -> str:
    """GTDB-Tk's threshold: Pfam with the gathering cutoff, TIGRFAM with the noise cutoff."""
    return "ga" if source == "pfam" else "nc"


def marker_url(accession: str) -> str:
    """The InterPro API address of one family's HMM (TIGRFAMs are served under NCBIfam)."""
    member = "pfam" if marker_source(accession) == "pfam" else "ncbifam"
    return f"{INTERPRO_API}/{member}/{accession}?annotation=hmm"


def parse_model(accession: str, data: bytes | str) -> MarkerModel:
    """Check and parse one family's HMM (gzip or plain text); WorkdirError when unusable."""
    if isinstance(data, bytes):
        if data[:2] == b"\x1f\x8b":
            try:
                data = gzip.decompress(data)
            except (OSError, EOFError) as exc:
                raise WorkdirError(f"{accession}: the HMM is not valid gzip ({exc}).") from exc
        text = data.decode("utf-8", errors="replace")
    else:
        text = data
    if not text.startswith("HMMER3/"):
        raise WorkdirError(f"{accession}: the response is not a HMMER3 model.")
    lines = text.splitlines()
    ends = [i for i, line in enumerate(lines) if line.strip() == "//"]
    if len(ends) != 1:
        raise WorkdirError(
            f"{accession}: expected one HMMER3 model ending with '//', found {len(ends)}."
        )
    fields: dict[str, str] = {}
    for line in lines[: ends[0]]:
        if line.startswith("HMM "):
            break
        key, _, value = line.partition(" ")
        fields.setdefault(key, value.strip())
    source = marker_source(accession)
    cutoff = marker_cutoff(source)
    acc = fields.get("ACC", "")
    if acc.split(".")[0] != accession:
        raise WorkdirError(f"{accession}: the model's ACC line reads '{acc or '(none)'}'.")
    if cutoff.upper() not in fields:
        raise WorkdirError(
            f"{accession}: the model has no {cutoff.upper()} line, which the "
            f"{source} families are searched with."
        )
    try:
        length = int(fields.get("LENG", ""))
    except ValueError as exc:
        raise WorkdirError(f"{accession}: the model has no valid LENG line.") from exc
    return MarkerModel(
        accession=accession,
        source=source,
        cutoff=cutoff,
        name=fields.get("NAME", ""),
        version=acc,
        length=length,
        text="\n".join(lines[: ends[0] + 1]) + "\n",
    )


def _split_models(text: str) -> list[str]:
    """The models of a concatenated HMMER3 file, each ending with its '//' line."""
    models, current = [], []
    for line in text.splitlines():
        current.append(line)
        if line.strip() == "//":
            models.append("\n".join(current) + "\n")
            current = []
    return models


# --- reading a database -------------------------------------------------------


def read_version(db: Path) -> dict:
    """The database's version.json as a dict; WorkdirError when unreadable."""
    path = db / VERSION_JSON
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise WorkdirError(f"{path} is unreadable: {exc}") from exc
    if not isinstance(data, dict) or data.get("v") != SCHEMA_VERSION:
        raise WorkdirError(f"{path} is not a version {SCHEMA_VERSION} marker database record.")
    return data


def read_index(db: Path) -> list[dict[str, str]]:
    """The rows of markers.tsv, one per (set, family)."""
    path = db / MARKERS_TSV
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise WorkdirError(f"{path} is unreadable: {exc}") from exc
    if not lines or tuple(lines[0].split("\t")) != INDEX_COLUMNS:
        raise WorkdirError(f"{path} does not start with the markers.tsv header.")
    return [dict(zip(INDEX_COLUMNS, line.split("\t"), strict=True)) for line in lines[1:] if line]


def hmm_path(db: Path, marker_set: str) -> Path:
    """The concatenated HMM file of one set."""
    return db / HMM_DIR / f"{marker_set}.hmm"


def check_markers_db(db: Path | None) -> Path:
    """The marker database directory, from the flag or ``REPGENR_MARKERS_DB``.

    Refused (UserInputError) with the download command when it is not given,
    is not a directory, or holds no version.json, index or HMM files.
    """
    hint = (
        f"The marker database is missing: {DOWNLOAD_HINT} and pass its path or set "
        f"{MARKERS_DB_ENV}."
    )
    if db is None:
        raw = os.environ.get(MARKERS_DB_ENV) or None
        db = Path(raw).expanduser() if raw else None
    if db is None:
        raise UserInputError(hint)
    if not db.is_dir():
        raise UserInputError(f"{db} is not a directory. {hint}")
    if not (db / VERSION_JSON).is_file():
        raise UserInputError(f"{db} holds no {VERSION_JSON}; not a marker database. {hint}")
    try:
        data = read_version(db)
    except WorkdirError as exc:
        raise UserInputError(f"{exc} {hint}") from exc
    missing = [s for s in data.get("sets", []) if not hmm_path(db, s).is_file()]
    if missing or not (db / MARKERS_TSV).is_file():
        what = ", ".join(str(hmm_path(db, s)) for s in missing) or str(db / MARKERS_TSV)
        raise UserInputError(f"{db} is incomplete: {what} is missing. {hint}")
    return db


def markers_db_version(db: Path) -> str:
    """``<sets>:<digest>`` for the fingerprint: the sets and a digest of markers.tsv.

    The index names every family's model version, so a database with any
    newer model gives a different value; a fresh download of the same
    models gives the same one.
    """
    data = read_version(db)
    try:
        digest = hashlib.sha256((db / MARKERS_TSV).read_bytes()).hexdigest()[:12]
    except OSError as exc:
        raise WorkdirError(f"{db / MARKERS_TSV} is unreadable: {exc}") from exc
    return f"{'+'.join(data.get('sets', []))}:{digest}"


def _existing_models(out: Path) -> tuple[list[str], dict[str, MarkerModel]]:
    """The sets and models of a database already at ``out`` (none for an empty dir)."""
    if not out.exists():
        return [], {}
    if not out.is_dir():
        raise UserInputError(f"--output {out} exists and is not a directory.")
    if not (out / VERSION_JSON).is_file():
        if any(out.iterdir()):
            raise UserInputError(
                f"--output {out} is not empty and holds no marker database; choose a new or "
                "empty directory."
            )
        return [], {}
    sets = list(read_version(out).get("sets", []))
    models: dict[str, MarkerModel] = {}
    for name in sets:
        path = hmm_path(out, name)
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            raise WorkdirError(f"{path} is unreadable: {exc}") from exc
        for block in _split_models(text):
            acc = next(
                (ln.split()[1].split(".")[0] for ln in block.splitlines() if ln.startswith("ACC ")),
                "",
            )
            if acc and acc not in models:
                models[acc] = parse_model(acc, block)
    return sets, models


# --- download -------------------------------------------------------------------


def _default_fetcher(timeout: float) -> Fetcher:
    from ..core.http import get_bytes

    return lambda url: get_bytes(url, timeout=(15, timeout))


def _fetch_one(accession: str, fetch: Fetcher, retries: int, backoff: float) -> MarkerModel | str:
    """The family's model, or the reason it could not be fetched after the retries."""
    from ..core.http import HTTPStatusError

    reason = ""
    for attempt in range(1, max(1, retries) + 1):
        try:
            return parse_model(accession, fetch(marker_url(accession)))
        except HTTPStatusError as exc:
            reason = str(exc)
            if exc.status == 404:
                break  # no such entry; a retry gives the same answer
        except WorkdirError as exc:
            reason = str(exc)
        if attempt < retries and backoff > 0:
            time.sleep(backoff * attempt)
    return reason


def _write_db(staging: Path, sets: Sequence[str], models: dict[str, MarkerModel]) -> None:
    (staging / HMM_DIR).mkdir(parents=True)
    rows = ["\t".join(INDEX_COLUMNS)]
    n_models: dict[str, int] = {}
    for name in sets:
        accessions = load_marker_list(name)
        with open(hmm_path(staging, name), "w", encoding="utf-8") as fo:
            for acc in accessions:
                m = models[acc]
                fo.write(m.text)
                rows.append(
                    "\t".join((name, acc, m.source, m.cutoff, m.name, m.version, str(m.length)))
                )
        n_models[name] = len(accessions)
    (staging / MARKERS_TSV).write_text("\n".join(rows) + "\n", encoding="utf-8")
    record = {
        "v": SCHEMA_VERSION,
        "sets": list(sets),
        "gtdbtk_lists": GTDBTK_LISTS,
        "fetched": datetime.now(UTC).date().isoformat(),
        "interpro_api": INTERPRO_API,
        "n_models": n_models,
    }
    (staging / VERSION_JSON).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def _swap_into_place(staging: Path, out: Path) -> None:
    """Rename ``staging`` to ``out``, replacing an earlier database there."""
    if not out.exists():
        staging.rename(out)
        return
    old = out.with_name(f".{out.name}.old-{os.getpid()}")
    out.rename(old)
    try:
        staging.rename(out)
    except BaseException:
        old.rename(out)
        raise
    shutil.rmtree(old, ignore_errors=True)


def download(
    out: Path,
    sets: Sequence[str],
    logger: logging.Logger,
    *,
    retries: int = 3,
    timeout: float = 60,
    threads: int = 8,
    fetch: Fetcher | None = None,
    backoff: float = 2.0,
) -> Path:
    """Fetch the marker models of ``sets`` from InterPro into the database at ``out``.

    Families already present in a database at ``out`` are kept and not
    fetched again; its other sets stay in the rewritten database. The new
    database is written beside ``out`` and renamed into place, so a failure
    leaves the earlier state (or no directory). Any family still missing
    after ``retries`` attempts is named in a WorkdirError.
    """
    from ..core.executors import parallel_map

    requested = list(dict.fromkeys(sets))
    for name in requested:
        load_marker_list(name)  # refuses an unknown set
    out = Path(out).expanduser().absolute()
    have_sets, models = _existing_models(out)
    all_sets = [s for s in MARKER_SETS if s in have_sets or s in requested]
    needed = list(dict.fromkeys(acc for s in all_sets for acc in load_marker_list(s)))
    todo = [acc for acc in needed if acc not in models]
    logger.info(
        "Marker database %s: %d families for %s; %d present, %d to fetch from InterPro",
        out,
        len(needed),
        "+".join(all_sets),
        len(needed) - len(todo),
        len(todo),
    )
    if not todo and set(all_sets) == set(have_sets):
        logger.info("Marker database %s is complete; nothing to fetch.", out)
        return out

    fetcher = fetch if fetch is not None else _default_fetcher(timeout)
    started = time.monotonic()
    results = parallel_map(
        lambda acc: _fetch_one(acc, fetcher, retries, backoff), todo, threads, logger=logger
    )
    failed = {acc: r for acc, r in zip(todo, results, strict=True) if isinstance(r, str)}
    if failed:
        listed = "; ".join(f"{acc} ({why})" for acc, why in sorted(failed.items()))
        raise WorkdirError(
            f"{len(failed)} of {len(todo)} marker models could not be fetched after "
            f"{retries} attempts: {listed}. Nothing was written; re-run the download."
        )
    for acc, r in zip(todo, results, strict=True):
        assert isinstance(r, MarkerModel)
        models[acc] = r

    out.parent.mkdir(parents=True, exist_ok=True)
    # mkdir rather than mkdtemp: the database keeps the umask's permissions
    # (mkdtemp creates 0700). The process id keeps two runs apart.
    staging = out.with_name(f".{out.name}.partial-{os.getpid()}")
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir()
    try:
        _write_db(staging, all_sets, models)
        _swap_into_place(staging, out)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    n_pfam = sum(1 for acc in needed if models[acc].source == "pfam")
    logger.info(
        "Marker database written to %s: %d fetched, %d reused, %d models (%d Pfam with GA, "
        "%d TIGRFAM with NC) in %.1f s",
        out,
        len(todo),
        len(needed) - len(todo),
        len(needed),
        n_pfam,
        len(needed) - n_pfam,
        time.monotonic() - started,
    )
    return out


def summary(db: Path) -> dict:
    """The database summary that ``repgenr markers-db info`` prints."""
    data = read_version(db)
    rows = read_index(db)
    per_set: dict[str, dict[str, int]] = {}
    for row in rows:
        entry = per_set.setdefault(row["set"], {"models": 0, "pfam": 0, "tigrfam": 0})
        entry["models"] += 1
        entry[row["source"]] = entry.get(row["source"], 0) + 1
    for name, entry in per_set.items():
        path = hmm_path(db, name)
        entry["hmm_bytes"] = path.stat().st_size if path.is_file() else 0
    return {
        "path": str(db),
        "sets": data.get("sets", []),
        "gtdbtk_lists": data.get("gtdbtk_lists", ""),
        "fetched": data.get("fetched", ""),
        "interpro_api": data.get("interpro_api", ""),
        "version": markers_db_version(db),
        "per_set": per_set,
        "families": len({row["accession"] for row in rows}),
    }
