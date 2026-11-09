"""The report steps, their order in the report and the stages they describe.

A report step reads the data files of one part of the pipeline and writes a
fragment (``report/<step>/``). Several stages can feed one step (metadata,
ingest and vgenome all write the genome selection), and one stage can feed
several steps (assemble writes the assembly table and the selection).
"""

from __future__ import annotations

import importlib
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..core.config import StageRecord
from .model import StepReport


@dataclass(frozen=True)
class StepInput:
    """What a step collector is given."""

    workdir: Path
    out_dir: Path
    formats: tuple[str, ...]
    records: dict[str, StageRecord] = field(default_factory=dict)


@dataclass(frozen=True)
class StepSpec:
    name: str
    title: str
    # Harness stage names whose records describe the step (methods text).
    stages: tuple[str, ...]
    # Module under repgenr.reporting.steps with sources(wd), present(wd) and
    # collect(StepInput) -> StepReport.
    module: str

    def _mod(self) -> Any:
        return importlib.import_module(f"repgenr.reporting.steps.{self.module}")

    def sources(self, workdir: Path) -> list[Path]:
        """The files the step reads (existing or not)."""
        return list(self._mod().sources(workdir))

    def present(self, workdir: Path) -> bool:
        """True when the workdir holds the outputs the step reports on."""
        return bool(self._mod().present(workdir))

    @property
    def collect(self) -> Callable[[StepInput], StepReport]:
        return self._mod().collect


# In pipeline order: the order of the sections of the HTML report.
STEPS: tuple[StepSpec, ...] = (
    StepSpec("reads", "Sequencing runs", ("reads", "reads_sketch", "reads_screen"), "reads"),
    StepSpec("assembly", "Genome assembly", ("assemble",), "assembly"),
    StepSpec(
        "selection",
        "Genome selection",
        ("metadata", "vmetadata", "vgenome", "ingest", "genome", "assemble"),
        "selection",
    ),
    StepSpec("glance", "Pairwise similarity (glance)", ("glance",), "glance"),
    StepSpec(
        "dereplication",
        "Dereplication",
        ("dereplicate", "derep_refine", "cluster_summary"),
        "dereplication",
    ),
)

STEP_NAMES: tuple[str, ...] = tuple(s.name for s in STEPS)
BY_NAME: dict[str, StepSpec] = {s.name: s for s in STEPS}

# Harness stage -> the steps whose fragments it refreshes when it finishes
# (or is skipped).
STAGE_TO_STEPS: dict[str, tuple[str, ...]] = {}
for _spec in STEPS:
    for _stage in _spec.stages:
        STAGE_TO_STEPS.setdefault(_stage, ())
        STAGE_TO_STEPS[_stage] += (_spec.name,)
