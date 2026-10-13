"""Workdir configuration and provenance (``repgenr.yaml``).

A single YAML file at the workdir root records, per stage, the tool chosen, the
parameters used and the resolved tool versions. This replaces the old scattered
``*_parameters.txt`` files and the ``str(dict)`` / ``pickle`` state blobs.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .errors import WorkdirError

CONFIG_FILENAME = "repgenr.yaml"
SCHEMA_VERSION = 1


@dataclass
class StageRecord:
    """Provenance for one completed (or in-progress) stage."""

    tool: str | None = None
    params: dict[str, Any] = field(default_factory=dict)
    tool_versions: dict[str, str] = field(default_factory=dict)
    completed: str | None = None  # ISO timestamp, set by caller
    fingerprint: str | None = None  # hash of the stage invocation, for resume
    inputs: dict[str, str] = field(default_factory=dict)  # input path -> digest
    # Tool -> directory put first on PATH for it (--bin-dir), for the tools
    # that used one. Provenance only: not part of the fingerprint.
    bin_dirs: dict[str, str] = field(default_factory=dict)

    @property
    def interrupted(self) -> bool:
        """True for a record without a completion stamp.

        The stage harness writes a record before a stage runs (a provisional
        one on a first run, or the last finished one with its stamp cleared
        on a re-run) and the stage stamps it when it finishes, so a record
        without a stamp is a run that failed or was killed. ``status`` and
        ``doctor`` both use this test, so they never disagree; a stage with
        no parameters (``cluster_summary``) has an empty record and counts.
        """
        return not self.completed

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "tool": self.tool,
            "params": self.params,
            "tool_versions": self.tool_versions,
            "completed": self.completed,
            "fingerprint": self.fingerprint,
            "inputs": self.inputs,
        }
        if self.bin_dirs:
            data["bin_dirs"] = self.bin_dirs
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StageRecord:
        return cls(
            tool=data.get("tool"),
            params=dict(data.get("params") or {}),
            tool_versions=dict(data.get("tool_versions") or {}),
            completed=data.get("completed"),
            fingerprint=data.get("fingerprint"),
            inputs=dict(data.get("inputs") or {}),
            bin_dirs=dict(data.get("bin_dirs") or {}),
        )


@dataclass
class Config:
    """In-memory view of ``repgenr.yaml``."""

    schema_version: int = SCHEMA_VERSION
    repgenr_version: str = ""
    stages: dict[str, StageRecord] = field(default_factory=dict)

    @classmethod
    def load(cls, workdir: str | os.PathLike[str]) -> Config:
        path = Path(workdir) / CONFIG_FILENAME
        if not path.exists():
            from .. import __version__

            return cls(repgenr_version=__version__)
        try:
            with open(path, encoding="utf-8") as fo:
                data = yaml.safe_load(fo) or {}
            stages = {
                name: StageRecord.from_dict(rec or {})
                for name, rec in (data.get("stages") or {}).items()
            }
            return cls(
                schema_version=data.get("schema_version", SCHEMA_VERSION),
                repgenr_version=data.get("repgenr_version", ""),
                stages=stages,
            )
        except (yaml.YAMLError, AttributeError, TypeError, ValueError, UnicodeDecodeError) as exc:
            # A hand-edited or damaged record would otherwise surface as a raw
            # traceback from every command that reads it, status and doctor
            # included.
            reason = str(exc).splitlines()[0] if str(exc) else type(exc).__name__
            raise WorkdirError(
                f"{path} is not a readable RepGenR2 record ({reason}). Restore it from a "
                "backup, or move it aside and re-run the stages (each re-runs once)."
            ) from exc

    def save(self, workdir: str | os.PathLike[str]) -> Path:
        path = Path(workdir) / CONFIG_FILENAME
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "schema_version": self.schema_version,
            "repgenr_version": self.repgenr_version,
            "stages": {name: rec.to_dict() for name, rec in self.stages.items()},
        }
        # This file carries every stage's provenance and resume fingerprint, so
        # write beside it and rename: a crash mid-write must not truncate it.
        tmp = path.with_name(path.name + ".tmp")
        try:
            with open(tmp, "w", encoding="utf-8") as fo:
                yaml.safe_dump(data, fo, sort_keys=False, default_flow_style=False)
            os.replace(tmp, path)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
        return path

    def record_stage(
        self,
        name: str,
        *,
        tool: str | None = None,
        params: dict[str, Any] | None = None,
        tool_versions: dict[str, str] | None = None,
        completed: str | None = None,
        fingerprint: str | None = None,
        inputs: dict[str, str] | None = None,
    ) -> StageRecord:
        record = StageRecord(
            tool=tool,
            params=dict(params or {}),
            tool_versions=dict(tool_versions or {}),
            completed=completed,
            fingerprint=fingerprint,
            inputs=dict(inputs or {}),
        )
        self.stages[name] = record
        return record
