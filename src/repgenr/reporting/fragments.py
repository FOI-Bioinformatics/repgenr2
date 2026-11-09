"""Writing, reading and refreshing step fragments (``report/<step>/``).

A fragment is rebuilt when it is missing or stale: when a file it was built
from changed (by digest), when the record of a stage it describes changed,
when the figure formats differ, or when it was written by another version of
repgenr. It is written into a hidden temporary folder and moved into place,
so a failure leaves the previous fragment as it was.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..core.config import Config, StageRecord
from ..core.contracts import atomic_path
from ..core.inputs import file_digest
from .model import STEP_JSON, MethodsInfo, StepReport
from .registry import StepInput, StepSpec

REPORT_DIR = "report"


def report_dir(workdir: Path) -> Path:
    return workdir / REPORT_DIR


def read_fragment(path: Path) -> StepReport | None:
    """The step report in fragment folder ``path``, or None."""
    try:
        data = json.loads((path / STEP_JSON).read_text(encoding="utf-8"))
        return StepReport.from_dict(data)
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _record_digest(record: StageRecord) -> str:
    data = record.to_dict()
    data.pop("fingerprint", None)
    data.pop("inputs", None)
    return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()


def _sources(
    workdir: Path, spec: StepSpec, records: dict[str, StageRecord], formats: tuple[str, ...]
) -> dict[str, str]:
    out: dict[str, str] = {}
    for path in spec.sources(workdir):
        if path.is_file():
            out[path.relative_to(workdir).as_posix()] = file_digest(path)
    for stage in spec.stages:
        if stage in records:
            out[f"#record:{stage}"] = _record_digest(records[stage])
    out["#formats"] = ",".join(formats)
    return out


def _records(workdir: Path) -> dict[str, StageRecord]:
    try:
        return dict(Config.load(workdir).stages)
    except Exception:  # noqa: BLE001 - a damaged record leaves the methods empty
        return {}


def refresh_step(
    workdir: Path,
    spec: StepSpec,
    *,
    formats: tuple[str, ...],
    force: bool = False,
) -> str:
    """Bring ``report/<step>/`` up to date. Returns "written", "current",
    "absent" (the step has no outputs to report) or "removed" (its outputs
    are gone, and so is the fragment)."""
    dest = report_dir(workdir) / spec.name
    if not spec.present(workdir):
        if dest.exists():
            shutil.rmtree(dest)
            return "removed"
        return "absent"
    records = _records(workdir)
    sources = _sources(workdir, spec, records, formats)
    old = read_fragment(dest)
    if (
        not force
        and old is not None
        and old.sources == sources
        and old.repgenr_version == __version__
    ):
        return "current"
    tmp = report_dir(workdir) / f".{spec.name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    try:
        report = spec.collect(StepInput(workdir, tmp, formats, records))
        report.sources = sources
        report.repgenr_version = __version__
        report.generated = datetime.now(UTC).isoformat(timespec="seconds")
        report.methods = [
            MethodsInfo(
                stage=stage,
                tool=records[stage].tool,
                tool_versions=dict(records[stage].tool_versions),
                params=dict(records[stage].params),
                completed=records[stage].completed,
            )
            for stage in spec.stages
            if stage in records
        ]
        write_fragment(tmp, report)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    swap = report_dir(workdir) / f".{spec.name}.old"
    if swap.exists():
        shutil.rmtree(swap)
    if dest.exists():
        dest.rename(swap)
    tmp.rename(dest)
    shutil.rmtree(swap, ignore_errors=True)
    return "written"


def write_fragment(folder: Path, report: StepReport) -> None:
    with atomic_path(folder / STEP_JSON) as tmp:
        tmp.write_text(
            json.dumps(report.to_dict(), indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
