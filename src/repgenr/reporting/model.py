"""The report of one step: key numbers, tables, figures and provenance.

A step report is written as ``step.json`` beside its tables and figures (a
"fragment" folder, ``report/<step>/``). The HTML report reads only fragment
folders, so fragments written in separate places (Nextflow tasks) can be put
together without a working directory. Tables and figures carry ids local to
the step; their numbers (Table 3, Figure 7) are given when the HTML report is
assembled.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

SCHEMA = "repgenr.report-step/1"
STEP_JSON = "step.json"


@dataclass(frozen=True)
class KeyNumber:
    """One headline quantity of a step, shown in the report summary."""

    label: str
    value: str | int | float | None
    unit: str = ""


@dataclass(frozen=True)
class Table:
    """A TSV file of the fragment, with its caption."""

    id: str
    caption: str
    file: str
    columns: tuple[str, ...]
    n_rows: int
    # Columns aligned right and set in tabular figures in the HTML.
    numeric: tuple[str, ...] = ()
    note: str = ""


@dataclass(frozen=True)
class FigureRef:
    """A figure of the fragment, written in one or more formats."""

    id: str
    caption: str
    alt: str
    files: dict[str, str]


@dataclass(frozen=True)
class MethodsInfo:
    """The record of one stage that contributed to the step."""

    stage: str
    tool: str | None = None
    tool_versions: dict[str, str] = field(default_factory=dict)
    params: dict[str, Any] = field(default_factory=dict)
    completed: str | None = None


@dataclass
class StepReport:
    step: str
    title: str
    summary: list[str] = field(default_factory=list)
    key_numbers: list[KeyNumber] = field(default_factory=list)
    tables: list[Table] = field(default_factory=list)
    figures: list[FigureRef] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    methods: list[MethodsInfo] = field(default_factory=list)
    # Relative source path -> digest, for the staleness check.
    sources: dict[str, str] = field(default_factory=dict)
    repgenr_version: str = ""
    generated: str = ""
    schema: str = SCHEMA

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StepReport:
        if data.get("schema") != SCHEMA:
            raise ValueError(f"not a {SCHEMA} step report: {data.get('schema')!r}")
        return cls(
            step=str(data["step"]),
            title=str(data["title"]),
            summary=[str(s) for s in data.get("summary", [])],
            key_numbers=[KeyNumber(**k) for k in data.get("key_numbers", [])],
            tables=[
                Table(**{**t, "columns": tuple(t["columns"]), "numeric": tuple(t["numeric"])})
                for t in data.get("tables", [])
            ],
            figures=[FigureRef(**f) for f in data.get("figures", [])],
            notes=[str(n) for n in data.get("notes", [])],
            methods=[MethodsInfo(**m) for m in data.get("methods", [])],
            sources=dict(data.get("sources", {})),
            repgenr_version=str(data.get("repgenr_version", "")),
            generated=str(data.get("generated", "")),
        )
