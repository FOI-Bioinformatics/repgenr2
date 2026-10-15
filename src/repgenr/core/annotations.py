"""The annotations cache: ``annotations/<record>/<record>.gff3`` and an index.

One row per annotated record in ``annotations/annotations.tsv``. A row is
current when the genome digest, the annotator, its version, the database
version and the extras match and the GFF3 exists; the annotate stage then
reuses it, annotates what is stale or missing, and removes records no longer
in the set. Partial output goes into a dot-prefixed directory beside the
final one and is renamed into place, so an interrupted run leaves no
half-written record that a later run would take for a finished one.
"""

from __future__ import annotations

import csv
import hashlib
import logging
import os
import shutil
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, fields
from pathlib import Path

from .contracts import ANNOTATIONS_TSV, atomic_replace
from .errors import WorkdirError
from .inputs import file_digest

GFF_SUFFIX = ".gff3"
_PARTIAL = ".partial"
SET_GENOMES = "genomes"
SET_OUTGROUP = "outgroup"


@dataclass(frozen=True)
class AnnotationRecord:
    record_name: str
    set: str  # genomes | outgroup
    genome_path: str  # workdir-relative, or absolute outside it
    genome_digest: str
    annotator: str
    annotator_version: str
    db_version: str
    extra: str  # sorted key=value, comma separated
    gff_digest: str
    completed: str


COLUMNS = tuple(f.name for f in fields(AnnotationRecord))


def extra_text(extra: Mapping[str, object]) -> str:
    """The tuning as a stable string: a boolean is normalised, so ``True`` and
    ``"true"`` describe one setting and do not force a re-annotation. ``1``
    and ``0`` are left alone; they are numbers to the tools that read them.
    """

    def value(raw: object) -> str:
        if isinstance(raw, bool):
            return "true" if raw else "false"
        text = str(raw).strip()
        low = text.lower()
        if low in ("true", "yes"):
            return "true"
        if low in ("false", "no"):
            return "false"
        return text

    return ",".join(f"{k}={value(extra[k])}" for k in sorted(extra))


def index_path(directory: Path) -> Path:
    return directory / ANNOTATIONS_TSV


def record_dir(directory: Path, record_name: str) -> Path:
    return directory / record_name


def gff_path(directory: Path, record_name: str) -> Path:
    return record_dir(directory, record_name) / f"{record_name}{GFF_SUFFIX}"


def read_index(directory: Path) -> dict[str, AnnotationRecord]:
    """record_name -> record; empty when there is no index yet."""
    path = index_path(directory)
    if not path.is_file():
        return {}
    out: dict[str, AnnotationRecord] = {}
    with open(path, encoding="utf-8", newline="") as fo:
        reader = csv.DictReader(fo, delimiter="\t")
        if reader.fieldnames is None or set(COLUMNS) - set(reader.fieldnames):
            raise WorkdirError(f"{path} does not hold the annotation index columns {COLUMNS}")
        for row in reader:
            rec = AnnotationRecord(**{c: row.get(c, "") or "" for c in COLUMNS})
            out[rec.record_name] = rec
    return out


def write_index(directory: Path, records: Mapping[str, AnnotationRecord]) -> Path:
    path = index_path(directory)
    with atomic_replace(path, newline="") as fo:
        writer = csv.writer(fo, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for name in sorted(records):
            rec = records[name]
            writer.writerow([getattr(rec, c) for c in COLUMNS])
    return path


def is_current(
    rec: AnnotationRecord | None,
    *,
    genome_digest: str,
    annotator: str,
    annotator_version: str,
    db_version: str,
    extra: str,
    gff: Path,
) -> bool:
    return (
        rec is not None
        and rec.genome_digest == genome_digest
        and rec.annotator == annotator
        and rec.annotator_version == annotator_version
        and rec.db_version == db_version
        and rec.extra == extra
        and gff.is_file()
    )


def partial_dir(directory: Path, record_name: str) -> Path:
    """A hidden temporary directory beside the record's final one."""
    return directory / f".{record_name}{_PARTIAL}"


def is_partial(name: str) -> bool:
    return name.startswith(".") and name.endswith(_PARTIAL)


def clear_partials(directory: Path) -> None:
    """Remove the temporary directories an interrupted run left behind."""
    if not directory.is_dir():
        return
    for entry in directory.iterdir():
        if entry.is_dir() and is_partial(entry.name):
            shutil.rmtree(entry, ignore_errors=True)


def publish(directory: Path, record_name: str) -> Path:
    """Rename the record's partial directory into place, replacing an older one.

    The previous directory is moved aside first and removed only once the new
    one is in place, so a crash between the two leaves one of them, never
    neither. A leftover ``.<name>.old.partial`` is cleared by
    :func:`clear_partials`.
    """
    final = record_dir(directory, record_name)
    tmp = partial_dir(directory, record_name)
    old = directory / f".{record_name}.old{_PARTIAL}"
    shutil.rmtree(old, ignore_errors=True)
    if final.exists():
        os.replace(final, old)
    os.replace(tmp, final)
    shutil.rmtree(old, ignore_errors=True)
    return final


def remove_stale(
    directory: Path,
    records: dict[str, AnnotationRecord],
    expected: Iterable[str],
    workdir: Path,
    logger: logging.Logger | None = None,
) -> int:
    """Drop the records whose genome file is gone; keep the rest.

    An annotation costs minutes, so a record outside the current set is kept
    (a run on the representatives after one with ``--all-genomes``, or a new
    dereplication, would otherwise throw away work that the next run needs
    again). Only a record whose genome no longer exists is removed, with its
    directory. Consumers read the records of their own set by name.
    """
    keep = set(expected)
    removed = 0
    kept_outside = 0
    for name in sorted(records):
        if name in keep:
            continue
        genome = records[name].genome_path
        path = (workdir / genome) if genome else None
        if path is not None and path.exists():
            kept_outside += 1
            continue
        shutil.rmtree(record_dir(directory, name), ignore_errors=True)
        del records[name]
        removed += 1
    if logger is not None:
        if removed:
            logger.info(
                "Removed %d annotation(s) whose genome is no longer in the workdir", removed
            )
        if kept_outside:
            logger.info(
                "Kept %d annotation(s) of genomes outside this set (they cost minutes each); "
                "they are reused if the set grows again",
                kept_outside,
            )
    return removed


def genome_digest(path: Path) -> str:
    return file_digest(path)


def gff_digest(path: Path) -> str:
    return file_digest(path)


def index_digest(records: Mapping[str, AnnotationRecord], names: Iterable[str]) -> str:
    """A digest of the index rows of ``names`` (the fields a downstream step
    depends on), for a reuse stamp."""
    h = hashlib.sha256()
    for name in sorted(set(names)):
        rec = records.get(name)
        if rec is None:
            h.update(f"{name}\tABSENT\n".encode())
            continue
        h.update(
            "\t".join(
                [
                    rec.record_name,
                    rec.annotator,
                    rec.annotator_version,
                    rec.db_version,
                    rec.extra,
                    rec.gff_digest,
                ]
            ).encode()
            + b"\n"
        )
    return h.hexdigest()
