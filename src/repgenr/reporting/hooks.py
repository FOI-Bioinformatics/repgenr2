"""The report refresh the stage harness runs after each stage."""

from __future__ import annotations

import logging
from pathlib import Path

from .figures import DEFAULT_FORMATS
from .fragments import refresh_step
from .registry import BY_NAME, STAGE_TO_STEPS


def refresh_for_stage(workdir: Path, stage: str, logger: logging.Logger) -> None:
    """Refresh the fragments of the steps ``stage`` feeds.

    Never raises: a report is written beside the results, and a failure to
    write it is logged as a warning without failing the stage.
    """
    for name in STAGE_TO_STEPS.get(stage, ()):
        try:
            state = refresh_step(workdir, BY_NAME[name], formats=DEFAULT_FORMATS)
        except Exception as exc:  # noqa: BLE001 - see docstring
            logger.warning("Could not write the report of step '%s': %s", name, exc)
            continue
        if state == "written":
            logger.debug("Report of step '%s' written to report/%s/", name, name)
