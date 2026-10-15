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

from .contracts import atomic_replace
from .errors import WorkdirError
from .inputs import file_digest

ANNOTATIONS_DIR = "annotations"
ANNOTATIONS_TSV = "annotations.tsv"
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
    return ",".join(f"{k}={extra[k]}" for k in sorted(extra))


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
    """Rename the record's partial directory into place, replacing an older one."""
    final = record_dir(directory, record_name)
    tmp = partial_dir(directory, record_name)
    if final.exists():
        shutil.rmtree(final)
    os.replace(tmp, final)
    return final


def remove_stale(
    directory: Path,
    records: dict[str, AnnotationRecord],
    expected: Iterable[str],
    logger: logging.Logger | None = None,
) -> int:
    """Drop records (and their directories) that are not in ``expected``."""
    keep = set(expected)
    removed = 0
    for name in sorted(records):
        if name in keep:
            continue
        shutil.rmtree(record_dir(directory, name), ignore_errors=True)
        del records[name]
        removed += 1
    if removed and logger is not None:
        logger.info("Removed %d annotation(s) of genomes no longer in the set", removed)
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
