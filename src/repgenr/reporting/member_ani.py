"""ANI of each cluster member to its representative, from the genome sketches.

The dereplicators do not keep the ANI of a member to its representative, so
the report estimates it from the workdir's genome sketches
(``sketches/<record name>.sig.zip``, FracMinHash at scaled=1000). The
estimate is the one ``dereplicate --tool sourmash`` clusters on: the average
containment c = (|A & B| / |A| + |A & B| / |B|) / 2 at k=31, and ANI = c^(1/k).
A pair whose sketches share too few hashes (fewer than 30 and fewer than a
tenth of the smaller sketch, :func:`repgenr.core.sourmash.containment_to_ani_matrix`)
has no estimate. The sketches are read without sourmash.
"""

from __future__ import annotations

import gzip
import json
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import record_name
from ..core.sourmash import MIN_SHARED_FRACTION, MIN_SHARED_HASHES

KSIZE = 31
SKETCH_SUFFIX = ".sig.zip"


@dataclass(frozen=True)
class MemberAni:
    representative: str
    member: str
    ani: float | None
    shared_hashes: int | None


def _signatures(path: Path) -> Iterable[dict]:
    with zipfile.ZipFile(path) as zf:
        for name in zf.namelist():
            if not name.startswith("signatures/"):
                continue
            data = zf.read(name)
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            loaded = json.loads(data)
            for sig in loaded if isinstance(loaded, list) else [loaded]:
                yield from sig.get("signatures", [])


def load_hashes(path: Path, ksize: int = KSIZE) -> frozenset[int] | None:
    """The hashes of the DNA sketch at ``ksize`` in a ``.sig.zip``; None when
    the file or the sketch is absent or unreadable."""
    try:
        for sketch in _signatures(path):
            if int(sketch.get("ksize", 0)) == ksize and (
                str(sketch.get("molecule", "dna")).lower() == "dna"
            ):
                return frozenset(int(h) for h in sketch.get("mins") or [])
    except (OSError, ValueError, zipfile.BadZipFile, KeyError):
        return None
    return None


def pair_ani(a: frozenset[int], b: frozenset[int], ksize: int = KSIZE) -> tuple[float | None, int]:
    """(ANI estimate or None when unknown, shared hashes) of two sketches."""
    if not a or not b:
        return None, 0
    shared = len(a & b)
    if shared < MIN_SHARED_HASHES and shared < MIN_SHARED_FRACTION * min(len(a), len(b)):
        return None, shared
    if shared == 0:
        return None, 0
    c = (shared / len(a) + shared / len(b)) / 2
    return c ** (1.0 / ksize), shared


def member_anis(
    clusters: Mapping[str, Sequence[str]], sketches_dir: Path
) -> tuple[list[MemberAni], int]:
    """The ANI of every member to its representative, and the number of
    genomes whose sketch was not found. ``clusters`` maps representative
    file names to member file names, as ``clusters.tsv``."""
    cache: dict[str, frozenset[int] | None] = {}
    missing: set[str] = set()

    def hashes(name: str) -> frozenset[int] | None:
        if name not in cache:
            cache[name] = load_hashes(sketches_dir / (record_name(name) + SKETCH_SUFFIX))
            if cache[name] is None:
                missing.add(name)
        return cache[name]

    rows: list[MemberAni] = []
    for rep, members in clusters.items():
        if not members:
            continue
        rep_hashes = hashes(rep)
        for member in members:
            member_hashes = hashes(member)
            if rep_hashes is None or member_hashes is None:
                rows.append(MemberAni(rep, member, None, None))
                continue
            ani, shared = pair_ani(rep_hashes, member_hashes)
            rows.append(MemberAni(rep, member, ani, shared))
    return rows, len(missing)
