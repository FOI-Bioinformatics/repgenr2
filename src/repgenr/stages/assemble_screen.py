"""The reads screen of ``assemble --screen-reads``: three checks on a run's
reads sketch, after the fetch and before an assembler is chosen.

``taxon_mismatch``
    ``sourmash gather`` of the reads sketch (k=31) against the GTDB sketch,
    resolved to a lineage with ``sourmash tax genome`` and the lineages CSV
    (the commands of the post-assembly classifier). The run is excused when
    the GTDB genus differs from the submitted genus, with the classifier's
    ``genus_renamed`` tolerance (family and species epithet agree).
``low_match_fraction``
    The abundance-weighted fraction of the reads assigned to the top species
    (``f_weighted_at_rank`` of ``tax genome`` at species rank, or the top
    match's ``f_unique_weighted`` when the column is absent) is below the
    minimum fraction.
``duplicate_isolate``
    The reads sketch is contained in the reads sketch of a run accepted
    earlier, of the same biosample, at an ANI estimate at or above the
    duplicate threshold. Runs of different biosamples are never compared:
    clonal isolates of different samples are distinct samples. Runs are
    compared in ``reads.tsv`` order (the largest first), so the larger run is
    kept.

A tool failure while screening excuses the run as ``screen_failed``, which
says nothing about the data. Each decision is written to
``assemblies/<run>/screen.json`` with the settings it was made with; the
reads sketch and the gather output under ``assemblies/<run>/screen/`` are
kept for an excused run, so the reason can be inspected and a later call with
other thresholds can screen it again without fetching its reads.
"""

from __future__ import annotations

import csv
import gzip
import json
import logging
import os
import threading
import zipfile
from collections.abc import Collection, Iterable, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from ..classifiers.sourmash import (
    GATHER_KSIZE,
    TaxRow,
    gather_command,
    read_gather_top,
    read_tax_genome,
    tax_genome,
)
from ..core import process
from ..core.containers import run_tool
from ..core.contracts import ReadRow, atomic_replace
from ..core.errors import RepGenRError, ToolExecutionError, UserInputError
from ..core.sourmash import SOURMASH_TOOL
from .taxon_match import DISAGREE, GENUS_RENAMED, genus_agreement, gtdb_tokens

SCREEN_STEP = "screen"
TAXON_MISMATCH = "taxon_mismatch"
LOW_MATCH_FRACTION = "low_match_fraction"
DUPLICATE_ISOLATE = "duplicate_isolate"
SCREEN_FAILED = "screen_failed"
PASS = "pass"
EXCUSED = "excused"

SCREEN_JSON = "screen.json"
SCREEN_DIR = "screen"
SCREEN_READS_TSV = "screen_reads.tsv"
_SCREEN_COLUMNS = [
    "run_accession",
    "top_match",
    "top_genus",
    "fraction",
    "duplicate_of",
    "decision",
    "reason",
]


@dataclass(frozen=True)
class ScreenSettings:
    """The settings a screen decision depends on; a change screens again."""

    min_fraction: float
    dup_ani: float
    gtdb_sketch: str
    gtdb_lineages: str
    ksize: int = GATHER_KSIZE

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class ScreenRecord:
    """One run's screen decision (``screen.json``)."""

    run_accession: str
    decision: str = PASS
    # The excuse for excused_runs.tsv; for a pass, "genus_renamed" or empty.
    reason: str = ""
    top_match: str = ""
    top_genus: str = ""
    fraction: float | None = None
    duplicate_of: str = ""
    ani: float | None = None
    settings: dict = field(default_factory=dict)
    # Whether the duplicate check was part of the decision. A pass from the
    # taxon and fraction checks alone (a kept sketch screened while planning)
    # still needs it.
    duplicates_checked: bool = False

    @property
    def passed(self) -> bool:
        return self.decision == PASS

    def made_with(self, settings: ScreenSettings) -> bool:
        return self.settings == settings.as_dict()

    def to_json(self) -> dict:
        return asdict(self)

    @classmethod
    def from_json(cls, data: object) -> ScreenRecord | None:
        if not isinstance(data, dict) or "run_accession" not in data:
            return None
        known = {f for f in cls.__dataclass_fields__}
        try:
            return cls(**{k: v for k, v in data.items() if k in known})
        except TypeError:
            return None


def write_screen_json(run_dir: Path, record: ScreenRecord) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    with atomic_replace(run_dir / SCREEN_JSON) as fo:
        fo.write(json.dumps(record.to_json(), indent=1))


def read_screen_json(run_dir: Path) -> ScreenRecord | None:
    try:
        data = json.loads((run_dir / SCREEN_JSON).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return ScreenRecord.from_json(data)


# --- screen_reads.tsv -------------------------------------------------------------------


def _row_cells(r: ScreenRecord) -> list[str]:
    return [
        r.run_accession,
        r.top_match,
        r.top_genus,
        "" if r.fraction is None else f"{r.fraction:.4f}",
        r.duplicate_of,
        r.decision,
        " ".join(r.reason.split()),
    ]


def read_screen_table(path: Path) -> dict[str, list[str]]:
    """Run accession -> its row of an existing screen_reads.tsv (cells as written)."""
    if not path.exists():
        return {}
    rows: dict[str, list[str]] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.reader(fo, delimiter="\t")
        next(reader, None)
        for row in reader:
            if row and len(row) >= len(_SCREEN_COLUMNS):
                rows[row[0]] = row[: len(_SCREEN_COLUMNS)]
    return rows


def write_screen_table(
    path: Path,
    records: Sequence[ScreenRecord],
    order: Sequence[str] = (),
    *,
    keep: Collection[str] | None = None,
) -> None:
    """Write screen_reads.tsv: the rows of ``records`` and, for runs not among
    them, the rows of the existing table; in ``order`` (run accessions), then
    the remaining runs in the order they were found. With ``keep``, rows of
    other runs (runs no longer in reads.tsv) are dropped."""
    rows = read_screen_table(path)
    for r in records:
        rows[r.run_accession] = _row_cells(r)
    if keep is not None:
        rows = {run: cells for run, cells in rows.items() if run in keep}
    ranked = {run: i for i, run in enumerate(order)}
    ordered = sorted(rows, key=lambda run: ranked.get(run, len(ranked)))
    with atomic_replace(path, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(_SCREEN_COLUMNS)
        for run in ordered:
            writer.writerow(rows[run])


# --- settings -------------------------------------------------------------------------------


GTDB_SKETCH_ENV = "REPGENR_GTDB_SKETCH"
GTDB_LINEAGES_ENV = "REPGENR_GTDB_LINEAGES"


def resolve_settings(
    min_fraction: float,
    dup_ani: float,
    gtdb_sketch: str | None,
    gtdb_lineages: str | None,
    *,
    command: str,
    min_fraction_flag: str,
    dup_ani_flag: str,
) -> ScreenSettings:
    """The settings of a screen from its flags, with the GTDB paths resolved.

    Refuses (UserInputError) a threshold outside its range, a missing GTDB
    sketch or lineages CSV (flag or environment variable) and a path that
    does not exist, before anything is written. ``command`` names the caller
    in the message. The caller checks sourmash (``preflight``).
    """
    if not 0.0 <= min_fraction <= 1.0:
        raise UserInputError(f"{min_fraction_flag} must be within [0, 1] (got {min_fraction}).")
    if not 0.0 < dup_ani <= 1.0:
        raise UserInputError(f"{dup_ani_flag} must be within (0, 1] (got {dup_ani}).")
    sketch = gtdb_sketch or os.environ.get(GTDB_SKETCH_ENV)
    lineages = gtdb_lineages or os.environ.get(GTDB_LINEAGES_ENV)
    if not sketch or not lineages:
        raise UserInputError(
            f"{command} needs the GTDB sketch and the lineages CSV published with it "
            f"(--gtdb-sketch and --gtdb-lineages, or {GTDB_SKETCH_ENV} and {GTDB_LINEAGES_ENV})."
        )
    for flag, value in (("--gtdb-sketch", sketch), ("--gtdb-lineages", lineages)):
        if not Path(value).expanduser().exists():
            raise UserInputError(
                f"{flag} {value} does not exist (from the flag or the environment variable)."
            )
    return ScreenSettings(
        min_fraction=min_fraction,
        dup_ani=dup_ani,
        gtdb_sketch=str(Path(sketch).expanduser().resolve()),
        gtdb_lineages=str(Path(lineages).expanduser().resolve()),
    )


# --- reads sketch comparison ------------------------------------------------------------


@dataclass(frozen=True)
class ReadsHashes:
    """The k-mer hashes of one reads sketch at one k-mer size, as sorted arrays.

    ``solid`` holds the hashes seen often enough to be taken as genomic (see
    :func:`solid_cutoff`); ``seen`` those seen at least twice. A hash seen
    once is mostly a sequencing error and is left out of both.
    """

    solid: np.ndarray
    seen: np.ndarray
    # The sketch's largest kept hash (2**64 / scaled); 0 when not recorded.
    max_hash: int = 0


def solid_cutoff(abundances: np.ndarray) -> int:
    """The abundance from which a hash counts as genomic.

    A tenth of the abundance-weighted median (about a tenth of the k-mer
    coverage), and at least 2. Error k-mers of an isolate are seen once or
    twice, genomic ones about as often as the coverage.
    """
    if abundances.size == 0:
        return 2
    order = np.sort(abundances)
    cumulative = np.cumsum(order, dtype=np.float64)
    median = order[int(np.searchsorted(cumulative, cumulative[-1] / 2.0))]
    return max(2, int(median) // 10)


def _signatures(path: Path) -> Iterable[dict]:
    """The sketch objects of a sourmash signature file (.sig.zip, .sig, .sig.gz)."""

    def parse(data: bytes) -> Iterable[dict]:
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        loaded = json.loads(data)
        for sig in loaded if isinstance(loaded, list) else [loaded]:
            yield from sig.get("signatures", [])

    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as zf:
            for name in zf.namelist():
                if name.startswith("signatures/"):
                    yield from parse(zf.read(name))
    else:
        yield from parse(path.read_bytes())


def load_reads_hashes(path: Path, ksize: int = GATHER_KSIZE) -> ReadsHashes | None:
    """The hashes of the reads sketch in ``path`` at ``ksize``; None when absent."""
    for sketch in _signatures(path):
        if int(sketch.get("ksize", 0)) != ksize or sketch.get("molecule", "dna").lower() != "dna":
            continue
        mins = np.asarray(sketch.get("mins") or [], dtype=np.uint64)
        max_hash = int(sketch.get("max_hash") or 0)
        abund = sketch.get("abundances")
        if abund is None:  # without abundances every hash counts
            return ReadsHashes(np.sort(mins), np.sort(mins), max_hash)
        counts = np.asarray(abund, dtype=np.int64)
        cutoff = solid_cutoff(counts)
        return ReadsHashes(np.sort(mins[counts >= cutoff]), np.sort(mins[counts >= 2]), max_hash)
    return None


def _downsampled(hashes: np.ndarray, max_hash: int) -> np.ndarray:
    """The hashes a sketch with this ``max_hash`` (a larger scaled) would keep."""
    return hashes[hashes <= np.uint64(max_hash)]


def containment_ani(
    query: ReadsHashes, reference: ReadsHashes, ksize: int = GATHER_KSIZE
) -> float | None:
    """ANI estimate from the containment of ``query``'s genomic hashes in the
    hashes ``reference`` saw at least twice: C ** (1 / k), as sourmash's
    containment ANI. Sketches of different scaled values are first reduced to
    the larger scaled (the smaller max_hash). None when the query has no
    genomic hashes."""
    solid, seen = query.solid, reference.seen
    bounds = [m for m in (query.max_hash, reference.max_hash) if m > 0]
    if len(bounds) == 2 and bounds[0] != bounds[1]:
        bound = min(bounds)
        solid, seen = _downsampled(solid, bound), _downsampled(seen, bound)
    if solid.size == 0:
        return None
    shared = np.isin(solid, seen, assume_unique=True).sum()
    containment = float(shared) / float(solid.size)
    return containment ** (1.0 / ksize) if containment > 0 else 0.0


# --- the gate -----------------------------------------------------------------------------


class ScreenGate:
    """Screens runs, possibly from concurrent workers, with ordered duplicate checks.

    The taxon and fraction checks of a run need only its own sketch. The
    duplicate check compares a run with the runs accepted before it: runs
    accepted while planning (finished runs, in ``reads.tsv`` order) and the
    pending runs earlier in ``reads.tsv`` order, of the same biosample, which
    a worker waits for (see :meth:`set_pending`). The result therefore does
    not depend on which worker finishes first.
    """

    def __init__(
        self,
        settings: ScreenSettings,
        logger: logging.Logger,
        *,
        duplicates: bool = True,
        gather_slots: int = 1,
        min_fraction_flag: str = "--screen-min-fraction",
    ) -> None:
        self.settings = settings
        self.duplicates = duplicates
        self.min_fraction_flag = min_fraction_flag
        self._logger = logger
        self._slots = threading.BoundedSemaphore(max(1, gather_slots))
        self._cond = threading.Condition()
        self._accepted: list[tuple[ReadRow, ReadsHashes]] = []
        self._order: dict[str, int] = {}
        self._rows: dict[str, ReadRow] = {}
        self._undecided: set[str] = set()
        self.records: dict[str, ScreenRecord] = {}

    # -- ordering ---------------------------------------------------------------------

    def set_pending(self, rows: Sequence[ReadRow]) -> None:
        """The runs the workers screen, in ``reads.tsv`` order."""
        with self._cond:
            for i, row in enumerate(rows):
                self._order[row.run_accession] = i
                self._rows[row.run_accession] = row
                self._undecided.add(row.run_accession)

    def release(self, run: str) -> None:
        """Mark a pending run decided (screened, excused before the screen, or failed)."""
        with self._cond:
            self._undecided.discard(run)
            self._cond.notify_all()

    def _wait_for_earlier(self, row: ReadRow) -> None:
        index = self._order.get(row.run_accession)
        if index is None:
            return
        with self._cond:
            while any(
                self._order[r] < index and _related(self._rows[r], row) for r in self._undecided
            ):
                if process.stop_requested.is_set():
                    raise ToolExecutionError(
                        ["screen"], -15, output=process.NOT_STARTED, tool="screen"
                    )
                self._cond.wait(timeout=1.0)

    # -- checks -------------------------------------------------------------------------

    def classify(self, row: ReadRow, sketch: Path, work: Path) -> ScreenRecord:
        """The taxon and fraction checks of one reads sketch (gather, tax genome)."""
        record = ScreenRecord(row.run_accession, settings=self.settings.as_dict())
        try:
            with self._slots:
                top, tax = _gather(row.run_accession, sketch, work, self.settings, self._logger)
        except (RepGenRError, OSError) as exc:
            if process.stop_requested.is_set():
                raise
            record.decision = EXCUSED
            record.reason = f"{SCREEN_FAILED}: gather: {exc}"
            self._logger.warning("%s: screen failed (%s)", row.run_accession, exc)
            return record
        if top is not None:
            record.top_match = top[0]
        if tax is None:
            record.decision = EXCUSED
            record.reason = f"{TAXON_MISMATCH}: no GTDB match above the gather threshold"
            return record
        tokens = gtdb_tokens(tax.lineage)
        record.top_genus = tokens[1]
        fraction = tax.f_weighted
        if fraction is None and top is not None:
            fraction = top[1]
        record.fraction = fraction
        submitted = (row.family or "unknown", row.genus or "unknown", row.species or "unknown")
        if submitted[1] != "unknown":
            agreement = genus_agreement(submitted, tokens)
            if agreement == DISAGREE:
                record.decision = EXCUSED
                record.reason = (
                    f"{TAXON_MISMATCH}: submitted genus {submitted[1]}, GTDB genus "
                    f"{tokens[1] or '?'} ({top[0] if top else tax.lineage})"
                )
                return record
            if agreement == GENUS_RENAMED:
                record.reason = GENUS_RENAMED
        if fraction is None or fraction < self.settings.min_fraction:
            record.decision = EXCUSED
            shown = "unknown" if fraction is None else f"{fraction:.3f}"
            record.reason = (
                f"{LOW_MATCH_FRACTION}: top species fraction {shown} below "
                f"{self.min_fraction_flag} {self.settings.min_fraction:g}"
            )
        return record

    def check_duplicate(self, row: ReadRow, record: ScreenRecord, sketch: Path) -> ScreenRecord:
        """Excuse ``record`` when an accepted run of the same biosample holds it.

        A passed record is registered as accepted, for the runs after it.
        """
        if not record.passed:
            return record
        if not self.duplicates:
            record.duplicates_checked = True
            return record
        self._wait_for_earlier(row)
        try:
            hashes = load_reads_hashes(sketch, self.settings.ksize)
        except (OSError, ValueError, zipfile.BadZipFile) as exc:
            self._logger.warning(
                "%s: reads sketch unreadable for the duplicate check (%s)", row.run_accession, exc
            )
            hashes = None
        with self._cond:
            candidates = [(r, h) for r, h in self._accepted if _related(r, row)]
        record.duplicates_checked = True
        if hashes is None:
            return record
        for other, other_hashes in candidates:
            ani = containment_ani(hashes, other_hashes, self.settings.ksize)
            if ani is not None and ani >= self.settings.dup_ani:
                record.decision = EXCUSED
                record.duplicate_of = other.run_accession
                record.ani = round(ani, 6)
                record.reason = (
                    f"{DUPLICATE_ISOLATE}: contained in {other.run_accession} at ANI "
                    f"{ani:.4f} (min {self.settings.dup_ani:g})"
                )
                return record
        with self._cond:
            self._accepted.append((row, hashes))
        return record

    def accept(self, row: ReadRow, sketch: Path) -> None:
        """Register a run accepted without :meth:`check_duplicate` (a reused decision)."""
        if not self.duplicates:
            return
        try:
            hashes = load_reads_hashes(sketch, self.settings.ksize)
        except (OSError, ValueError, zipfile.BadZipFile):
            hashes = None
        if hashes is not None:
            with self._cond:
                self._accepted.append((row, hashes))

    def judge(
        self,
        row: ReadRow,
        sketch: Path | None,
        run_dir: Path,
        prior: ScreenRecord | None = None,
    ) -> ScreenRecord:
        """The full screen of one pending run, written to ``screen.json``.

        ``prior``, a pass of the taxon and fraction checks with the current
        settings, skips the gather. Releases the run in any case.
        """
        try:
            if sketch is None or not sketch.is_file():
                record = ScreenRecord(
                    row.run_accession,
                    EXCUSED,
                    f"{SCREEN_FAILED}: no reads sketch",
                    settings=self.settings.as_dict(),
                )
            elif prior is not None and prior.passed and prior.made_with(self.settings):
                record = prior
            else:
                record = self.classify(row, sketch, run_dir / SCREEN_DIR)
            if sketch is not None and record.passed:
                record = self.check_duplicate(row, record, sketch)
        finally:
            self.release(row.run_accession)
        self.records[row.run_accession] = record
        write_screen_json(run_dir, record)
        log = self._logger.info if record.passed else self._logger.warning
        log(
            "%s: screen %s (%s%s)",
            row.run_accession,
            "passed" if record.passed else "excused the run",
            f"top {record.top_match or '-'}, fraction "
            + ("-" if record.fraction is None else f"{record.fraction:.3f}"),
            f"; {record.reason}" if record.reason else "",
        )
        return record


def _related(a: ReadRow, b: ReadRow) -> bool:
    """Runs compared by the duplicate check: those of the same biosample."""
    return bool(a.biosample and a.biosample == b.biosample)


def _gather(
    run: str, sketch: Path, work: Path, settings: ScreenSettings, logger: logging.Logger
) -> tuple[tuple[str, float | None] | None, TaxRow | None]:
    """Gather one reads sketch and resolve it; (top match, tax row or None)."""
    work.mkdir(parents=True, exist_ok=True)
    gather_csv = work / "gather.csv"
    gather_csv.unlink(missing_ok=True)
    mounts = sorted(
        {
            str(Path(settings.gtdb_sketch).resolve().parent),
            str(Path(settings.gtdb_lineages).resolve().parent),
            str(sketch.resolve().parent),
            str(work.resolve()),
        }
    )
    run_tool(
        SOURMASH_TOOL,
        gather_command(sketch, settings.gtdb_sketch, gather_csv, ksize=settings.ksize),
        logger=logger,
        log_prefix=f"sourmash {run}",
        extra_mounts=mounts,
    )
    top = read_gather_top(gather_csv)
    if top is None:
        return None, None
    table = tax_genome(
        SOURMASH_TOOL,
        [gather_csv],
        settings.gtdb_lineages,
        work / "tax",
        logger=logger,
        mounts=mounts,
        containment_threshold=0.0,
    )
    return top, read_tax_genome(table).get(run)
