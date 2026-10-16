"""Grouped --help: every command sits in one panel, in pipeline order."""

from __future__ import annotations

import re

import pytest
import typer
from typer.testing import CliRunner

from repgenr.cli import base
from repgenr.cli.base import COMMAND_ORDER, COMMAND_PANELS
from repgenr.cli.main import app

CLICK = typer.main.get_command(app)


def _help() -> str:
    out = CliRunner().invoke(app, ["--help"]).output
    return re.sub(r"\x1b\[[0-9;]*m", "", out)


def test_every_command_is_in_exactly_one_panel() -> None:
    assigned = [c for cmds in COMMAND_PANELS.values() for c in cmds]
    assert sorted(assigned) == sorted(CLICK.commands)
    assert len(assigned) == len(set(assigned))


def test_registered_panel_matches_mapping() -> None:
    for panel, cmds in COMMAND_PANELS.items():
        for name in cmds:
            assert CLICK.commands[name].rich_help_panel == panel, name


def test_command_order_follows_mapping() -> None:
    assert CLICK.list_commands(None) == list(COMMAND_ORDER)


def test_help_shows_panels_in_order() -> None:
    out = _help()
    positions = [out.index(panel) for panel in COMMAND_PANELS]
    assert positions == sorted(positions)
    rows = [m.group(1) for m in re.finditer(r"^\u2502 (\S+) ", out, re.MULTILINE)]
    listed = [r for r in rows if r in COMMAND_ORDER]
    assert listed == list(COMMAND_ORDER)
    assert "\u2500 Commands " not in out


def test_help_shows_stage_chains() -> None:
    out = _help()
    assert "ingest -> dereplicate" in out
    for panel in COMMAND_PANELS:
        assert panel in out


@pytest.mark.parametrize(
    "chain",
    [base.PIPELINE_BACTERIAL, base.PIPELINE_VIRAL, base.PIPELINE_LOCAL, base.PIPELINE_READS],
)
def test_epilog_uses_pipeline_tuples(chain: tuple[str, ...]) -> None:
    assert " -> ".join(chain) in base.APP_EPILOG


def _help_of(command: str, flag: str) -> str | None:
    for param in CLICK.commands[command].params:
        if flag in param.opts or flag in getattr(param, "secondary_opts", []):
            return param.help
    raise AssertionError(f"{command} has no {flag}")


_STEP_VERSIONS = [
    "dereplicate-chunk",
    "dereplicate-refine",
    "annotate-build",
    "phylo-build",
    "tree2tax-relations",
    "dereplicate-merge",
    "assemble-run",
    "genome-qc",
]

SHARED_FLAGS: list[tuple[str, str, str]] = (
    [("metadata", "run", f) for f in ("--gtdb-version", "--release", "--metadata-path")]
    + [("metadata", "run", f) for f in ("--nodownload", "--limit", "--dataset", "--level")]
    + [("snptype", "phylo", f) for f in ("--all-genomes", "--reference")]
    + [("snptype", "phylo", "--allow-incomplete")]
    + [("phylo", "run", f) for f in ("--all-genomes", "--no-outgroup", "--bootstrap")]
    + [("phylo", "run", f) for f in ("--reference", "--msa-source", "--aligner-arg")]
    + [("phylo", "run", "--allow-incomplete")]
    + [("phylo", "phylo-build", f) for f in ("--no-outgroup", "--bootstrap", "--reference")]
    + [("phylo", "phylo-build", f) for f in ("--msa-source", "--aligner-arg", "--mask")]
    + [("dereplicate", "run", f) for f in ("--tool-arg", "--keeper", "--pre-primary-ani")]
    + [("dereplicate", "run", f) for f in ("--pre-secondary-ani", "--process-size")]
    # --num-processes is not shared: on run it also sets the annotate workers.
    + [("dereplicate", "run", f) for f in ("--reduce", "--target-reps")]
    + [("dereplicate", "run", "--allow-incomplete")]
    + [("dereplicate", "run", f) for f in ("--refine", "--refine-ani", "--refine-file")]
    + [("dereplicate", "derep-refine", f) for f in ("--refine", "--refine-ani", "--refine-file")]
    + [("dereplicate", "dereplicate-chunk", "--tool-arg")]
    + [("dereplicate", "dereplicate-refine", f) for f in ("--tool-arg", "--refine", "--refine-ani")]
    + [("dereplicate", "dereplicate-refine", "--refine-file")]
    + [("annotate", "annotate-build", f) for f in ("--tool", "--bakta-db")]
    + [("dereplicate", "dereplicate-merge", "--tool-arg")]
    + [("tree2tax", "run", f) for f in ("--remove-outgroup", "--node-basename", "--root-name")]
    + [("tree2tax", "run", "--include-dereplicated")]
    + [("tree2tax", "tree2tax-relations", f) for f in ("--remove-outgroup", "--node-basename")]
    + [("tree2tax", "tree2tax-relations", f) for f in ("--root-name", "--include-dereplicated")]
    + [("tree2tax", "run", f) for f in ("--collapse-support", "--collapse-length")]
    + [("tree2tax", "tree2tax-relations", f) for f in ("--collapse-support", "--collapse-length")]
    + [("genome", "genome-fetch", "--keep-files")]
    + [("metadata", "run", "--workdir"), ("dereplicate", "phylo", "--workdir")]
    + [("genome-fetch", c, "--versions-out") for c in _STEP_VERSIONS]
)


@pytest.mark.parametrize(("command_a", "command_b", "flag"), SHARED_FLAGS)
def test_shared_flag_has_one_help_string(command_a: str, command_b: str, flag: str) -> None:
    help_a = _help_of(command_a, flag)
    assert help_a
    assert help_a == _help_of(command_b, flag)


def test_platform_flags_differ() -> None:
    global_param = next((p for p in CLICK.params if "--platform" in p.opts), None)
    assert global_param is not None, "global --platform is missing"
    global_help = global_param.help
    run_help = _help_of("run", "--platform")
    assert global_help and run_help
    assert global_help != run_help
