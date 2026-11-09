"""Scientific reports of a run: per-step tables and figures, and one HTML page.

Each report step (:mod:`.registry`) reads the data files of one part of the
pipeline and writes a fragment folder, ``report/<step>/``, with TSV tables,
figures (PNG and SVG) and ``step.json``. The stage harness refreshes the
fragments of a stage when it finishes or is skipped; ``repgenr report``
refreshes all of them and writes ``report/report.html``.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

from .. import __version__
from ..core.contracts import atomic_path
from .figures import DEFAULT_FORMATS
from .fragments import read_fragment, refresh_step, report_dir
from .html import Provenance, build_html
from .model import StepReport
from .registry import BY_NAME, STEP_NAMES, STEPS

REPORT_HTML = "report.html"
# A report above this size is written, with a warning.
LARGE_REPORT = 25 * 1024 * 1024

__all__ = [
    "REPORT_HTML",
    "STEP_NAMES",
    "Provenance",
    "collect_fragments",
    "refresh",
    "write_html",
]


def refresh(
    workdir: Path,
    steps: Sequence[str] | None = None,
    *,
    formats: Sequence[str] = DEFAULT_FORMATS,
    force: bool = False,
    strict: bool = False,
    logger: logging.Logger,
) -> dict[str, str]:
    """Bring the fragments of ``steps`` (all when None) up to date.

    Returns step -> state ("written", "current", "absent", "removed" or
    "failed"). A failing step is logged as a warning, unless ``strict``.
    """
    states: dict[str, str] = {}
    for spec in STEPS:
        if steps is not None and spec.name not in steps:
            continue
        try:
            states[spec.name] = refresh_step(workdir, spec, formats=tuple(formats), force=force)
        except Exception as exc:
            if strict:
                raise
            logger.warning("Could not write the report of step '%s': %s", spec.name, exc)
            states[spec.name] = "failed"
    return states


def collect_fragments(
    folders: Sequence[Path], steps: Sequence[str] | None = None
) -> list[tuple[Path, StepReport]]:
    """The readable fragments among ``folders``, in step order. A folder
    may be a fragment or hold fragments (``report/``)."""
    found: dict[str, tuple[Path, StepReport]] = {}
    for folder in folders:
        candidates = [folder] if (folder / "step.json").is_file() else sorted(folder.iterdir())
        for path in candidates:
            if path.name.startswith(".") or not path.is_dir():
                continue
            rep = read_fragment(path)
            if rep is None or (steps is not None and rep.step not in steps):
                continue
            found.setdefault(rep.step, (path, rep))
    order = {name: i for i, name in enumerate(STEP_NAMES)}
    return sorted(found.values(), key=lambda fr: (order.get(fr[1].step, len(order)), fr[1].step))


def write_html(
    out: Path,
    fragments: list[tuple[Path, StepReport]],
    prov: Provenance,
    logger: logging.Logger,
) -> Path:
    """Write the report page to ``out``."""
    text = build_html(fragments, prov)
    out.parent.mkdir(parents=True, exist_ok=True)
    with atomic_path(out) as tmp:
        tmp.write_text(text, encoding="utf-8")
    size = out.stat().st_size
    if size > LARGE_REPORT:
        logger.warning(
            "The report is %.0f MB; large tables and figures are embedded in it.",
            size / 1e6,
        )
    return out


def now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")


def default_title(workdir: Path | None) -> str:
    return f"RepGenR2 report: {workdir.name}" if workdir else "RepGenR2 report"


def new_provenance(workdir: Path | None, title: str | None = None) -> Provenance:
    return Provenance(
        title=title or default_title(workdir),
        generated=now(),
        repgenr_version=__version__,
        workdir=str(workdir) if workdir else "",
    )


def fragments_dir(workdir: Path) -> Path:
    return report_dir(workdir)


def read_versions_yml(path: Path) -> list[tuple[str, str]]:
    """(software, version) from a Nextflow ``software_versions.yml``: one map
    per process of tool -> version. A tool reported with different versions
    by different processes is listed once per process."""
    import yaml

    from ..core.versions import merge_stage_versions

    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    per_process = {
        str(proc): {str(k): str(v) for k, v in (tools or {}).items()}
        for proc, tools in data.items()
        if isinstance(tools, dict)
    }
    return sorted(merge_stage_versions(per_process).items())


def known_step(name: str) -> bool:
    return name in BY_NAME
