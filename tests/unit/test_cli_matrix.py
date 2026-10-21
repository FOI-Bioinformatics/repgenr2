"""The CLI matrix (tests/audit/cli_matrix.yaml) against the real command tree.

Every (command, flag) record is checked for: presence in the Click tree with
the same aliases; help text; the parameter it sets (by invoking the command
with the recorder in place of the stage dispatch); the recorded validation
kind (an invalid value exits non-zero naming the flag); live and docs
references that resolve; and short-alias collisions. The rendered Markdown
must be in sync with the YAML.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path
from typing import Any

import pytest
import typer
import yaml
from typer.testing import CliRunner

from repgenr.cli import base as cli_base
from repgenr.cli import (
    cmd_bacterial,
    cmd_ingest,
    cmd_misc,
    cmd_phylo,
    cmd_reads,
    cmd_run,
    cmd_viral,
)
from repgenr.cli.main import app
from repgenr.core import containers

ROOT = Path(__file__).resolve().parents[2]
MATRIX_PATH = ROOT / "tests" / "audit" / "cli_matrix.yaml"
RENDERED_PATH = ROOT / "docs" / "audit" / "cli-matrix.md"
SKIP_OPTS = {"--help", "--install-completion", "--show-completion"}
VALIDATION_KINDS = {"choice", "registry", "unit_interval", "range", "callback", "stage", "none"}
# Flags whose help text is still missing; must stay empty.
ALLOW_MISSING_HELP: set[tuple[str, str]] = set()

_runner = CliRunner()


def _load_matrix() -> dict:
    with open(MATRIX_PATH, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


MATRIX = _load_matrix()
CLICK = typer.main.get_command(app)


def _click_flags(cmd) -> dict[str, Any]:
    out = {}
    for p in cmd.params:
        opts = getattr(p, "opts", None)
        # The matrix records options; a positional argument (markers-db
        # ACTION and DIRECTORY) is listed in the rendered reference only.
        if not opts or opts[0] in SKIP_OPTS or getattr(p, "param_type_name", "") == "argument":
            continue
        long = next(o for o in opts if o.startswith("--"))
        out[long] = p
    return out


def _records() -> list[tuple[str, str, dict]]:
    rows = [("<global>", flag, rec) for flag, rec in MATRIX["global_flags"].items()]
    for name, cmd in MATRIX["commands"].items():
        rows += [(name, flag, rec) for flag, rec in cmd["flags"].items()]
    return rows


RECORDS = _records()
IDS = [f"{c} {f}" for c, f, _ in RECORDS]


# --- structure ----------------------------------------------------------------


def test_command_set_matches_click_tree() -> None:
    assert set(MATRIX["commands"]) == set(CLICK.commands)


@pytest.mark.parametrize("name", sorted(MATRIX["commands"]))
def test_flag_set_matches_click_tree(name: str) -> None:
    assert set(MATRIX["commands"][name]["flags"]) == set(_click_flags(CLICK.commands[name]))


def test_global_flag_set_matches_click_tree() -> None:
    assert set(MATRIX["global_flags"]) == set(_click_flags(CLICK))


def _param(command: str, flag: str):
    cmd = CLICK if command == "<global>" else CLICK.commands[command]
    return _click_flags(cmd)[flag]


@pytest.mark.parametrize(("command", "flag", "rec"), RECORDS, ids=IDS)
def test_record_shape(command: str, flag: str, rec: dict) -> None:
    p = _param(command, flag)
    assert rec.get("aliases", []) == [o for o in p.opts if o != flag], "short aliases drifted"
    assert rec.get("secondary", []) == list(getattr(p, "secondary_opts", []) or [])
    assert rec.get("envvar") == getattr(p, "envvar", None)
    assert rec["validated"] in VALIDATION_KINDS
    assert "live" in rec and "nextflow" in rec and "docs" in rec
    if p.is_flag:
        assert "value" not in rec, "boolean flags take no sample value"
    elif not rec["param"].startswith("n/a") and "wiring" not in rec:
        assert "value" in rec, "a valued option needs a sample value for the wiring test"
    if "wiring" in rec:
        assert rec["wiring"].startswith("n/a: "), "wiring is either absent or an n/a reason"


@pytest.mark.parametrize(("command", "flag", "rec"), RECORDS, ids=IDS)
def test_help_text_present(command: str, flag: str, rec: dict) -> None:
    if (command, flag) in ALLOW_MISSING_HELP:
        pytest.xfail("help text pending")
    assert _param(command, flag).help, f"{command} {flag} has no help text"


# --- wiring ---------------------------------------------------------------------


@pytest.fixture
def recorder(monkeypatch, tmp_path: Path):
    """Capture the params object each command dispatches instead of running it."""
    calls: list[Any] = []

    def fake_run(stage, workdir, build, *, create=False):
        calls.append(("workdir", workdir))
        calls.append(build())

    for mod in (cmd_bacterial, cmd_viral, cmd_phylo, cmd_misc, cmd_run, cmd_ingest, cmd_reads):
        monkeypatch.setattr(mod, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)

    def fake_step(params, logger=None, *a, **k):
        calls.append(params)

    for cmd in MATRIX["commands"].values():
        if cmd["dispatch"].startswith("step:"):
            module, _, func = cmd["dispatch"][5:].rpartition(".")
            monkeypatch.setattr(importlib.import_module(module), func, fake_step)
    return calls


def _placeholders(tmp_path: Path) -> dict[str, str]:
    d = tmp_path / "d"
    d.mkdir(exist_ok=True)
    (d / "a.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    alt = tmp_path / "alt"
    alt.mkdir(exist_ok=True)
    f = tmp_path / "f.txt"
    f.write_text("x\n", encoding="utf-8")
    alt_f = tmp_path / "g.txt"
    alt_f.write_text("y\n", encoding="utf-8")
    fofn = tmp_path / "fofn.txt"
    fofn.write_text(f"{d / 'a.fasta'}\n", encoding="utf-8")
    alt_fofn = tmp_path / "fofn2.txt"
    alt_fofn.write_text(f"{d / 'a.fasta'}\n{d / 'a.fasta'}\n", encoding="utf-8")
    return {
        "{wd}": str(tmp_path / "wd"),
        "{dir}": str(d),
        "{alt_dir}": str(alt),
        "{file}": str(f),
        "{alt_file}": str(alt_f),
        "{fofn}": str(fofn),
        "{alt_fofn}": str(alt_fofn),
        "{out}": str(tmp_path / "out"),
    }


def _fill(args: list[str], ph: dict[str, str]) -> list[str]:
    return [ph.get(a, a) for a in args]


def _argv(command: str, flag: str, rec: dict, ph: dict[str, str], *, use_flag: bool) -> list[str]:
    cmd = MATRIX["commands"].get(command)
    base = _fill(cmd["base_args"], ph) if cmd else ["status", "-wd", ph["{wd}"]]
    # with_args are part of both invocations so only the flag itself differs.
    extra: list[str] = _fill(rec.get("with_args", []), ph)
    if use_flag:
        p = _param(command, flag)
        # A boolean whose default is on is exercised through its --no-X form.
        if p.is_flag and rec.get("secondary") and p.default is True:
            extra = [*extra, rec["secondary"][0]]
        elif p.is_flag:
            extra = [*extra, flag]
        else:
            extra = [*extra, flag, ph.get(str(rec["value"]), str(rec["value"]))]
    if command == "<global>":
        return extra + base
    # Base args first, so a flag repeated in base (e.g. -wd) is overridden.
    return [command, *base, *extra]


def _lookup(calls: list[Any], param: str):
    """Resolve ``Class.field.sub`` (or ``workdir``) on the recorded objects."""
    if param == "workdir":
        return next(v for k, v in (c for c in calls if isinstance(c, tuple)) if k == "workdir")
    cls, _, path = param.partition(".")
    for obj in calls:
        if isinstance(obj, tuple):
            continue
        if type(obj).__name__ == cls:
            if not path:
                return True  # a bare class name asks whether that stage was dispatched
            for attr in path.split("."):
                obj = getattr(obj, attr)
            return obj
    if not path:
        return False
    raise AssertionError(f"no dispatched params object of type {cls}")


WIRED = [
    (c, f, r)
    for c, f, r in RECORDS
    if c != "<global>" and not r["param"].startswith("n/a") and "wiring" not in r
]


@pytest.mark.parametrize(("command", "flag", "rec"), WIRED, ids=[f"{c} {f}" for c, f, _ in WIRED])
def test_flag_reaches_its_parameter(recorder, tmp_path: Path, command, flag, rec) -> None:
    ph = _placeholders(tmp_path)
    baseline = _runner.invoke(app, _argv(command, flag, rec, ph, use_flag=False))
    assert baseline.exit_code == 0, baseline.output
    before = _lookup(recorder, rec["param"])
    recorder.clear()
    result = _runner.invoke(app, _argv(command, flag, rec, ph, use_flag=True))
    assert result.exit_code == 0, result.output
    after = _lookup(recorder, rec["param"])
    assert after != before, f"{command} {flag} did not change {rec['param']}"


GLOBAL_WIRED = [
    (f, r) for f, r in MATRIX["global_flags"].items() if not r["param"].startswith("n/a")
]


@pytest.mark.parametrize(("flag", "rec"), GLOBAL_WIRED, ids=[f for f, _ in GLOBAL_WIRED])
def test_global_flag_reaches_process_state(monkeypatch, tmp_path: Path, flag, rec) -> None:
    ph = _placeholders(tmp_path)
    # The probe runs `status`, which exits 3 on a missing -wd and 0 on an
    # existing directory without repgenr.yaml.
    Path(ph["{wd}"]).mkdir(exist_ok=True)
    monkeypatch.setattr(containers, "_CONFIG", containers.ContainerConfig())
    monkeypatch.setitem(cli_base._RUN_STATE, "force", False)

    def read():
        kind, _, attr = rec["param"].partition(".")
        if kind == "container":
            return getattr(containers.get_config(), attr)
        return cli_base._RUN_STATE[attr]

    assert _runner.invoke(app, _argv("<global>", flag, rec, ph, use_flag=False)).exit_code == 0
    before = read()
    assert _runner.invoke(app, _argv("<global>", flag, rec, ph, use_flag=True)).exit_code == 0
    assert read() != before, f"{flag} did not change {rec['param']}"


def test_version_flag_prints_version() -> None:
    result = _runner.invoke(app, ["--version"])
    assert result.exit_code == 0 and result.output.startswith("repgenr ")


# --- validation -------------------------------------------------------------------

_BAD = {"choice": "__bogus__", "registry": "__bogus__", "unit_interval": "1.5"}
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
VALIDATED = [
    (c, f, r)
    for c, f, r in RECORDS
    if r["validated"] in ("choice", "registry", "unit_interval", "range", "callback")
]


@pytest.mark.parametrize(
    ("command", "flag", "rec"), VALIDATED, ids=[f"{c} {f}" for c, f, _ in VALIDATED]
)
def test_invalid_value_is_rejected_naming_the_flag(monkeypatch, tmp_path, command, flag, rec):
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    ph = _placeholders(tmp_path)
    Path(ph["{wd}"]).mkdir()
    bad = rec.get("invalid", _BAD.get(rec["validated"]))
    assert bad is not None, "range/callback records need an explicit invalid value"
    rec_bad = {**rec, "value": bad}
    result = _runner.invoke(app, _argv(command, flag, rec_bad, ph, use_flag=True))
    assert result.exit_code != 0, f"{command} {flag}={bad} was accepted"
    # Rich colours the usage panel under CI (FORCE_COLOR) and splits the flag
    # name with escape codes; strip them before matching.
    text = _ANSI.sub("", result.output + str(result.exception or ""))
    assert rec.get("reported_as", flag) in text, f"rejection does not name {flag}: {text}"


# --- references -----------------------------------------------------------------


@pytest.mark.parametrize(("command", "flag", "rec"), RECORDS, ids=IDS)
def test_live_reference_resolves(command, flag, rec) -> None:
    live = rec["live"]
    if live.startswith(("n/a:", "todo: PR-")):
        return
    path, _, node = live.partition("::")
    assert path.startswith("tests/live/"), live
    source = (ROOT / path).read_text(encoding="utf-8")
    assert re.search(rf"^def {re.escape(node)}\(", source, re.M), f"{live} not found"


@pytest.mark.parametrize(("command", "flag", "rec"), RECORDS, ids=IDS)
def test_docs_references_resolve(command, flag, rec) -> None:
    for doc in rec["docs"]:
        assert (ROOT / doc).is_file(), doc
        assert flag in (ROOT / doc).read_text(encoding="utf-8"), f"{doc} lacks {flag}"


SCHEMA_PATH = ROOT / "nextflow" / "nextflow_schema.json"


def _schema_params() -> set[str]:
    import json

    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    keys: set[str] = set(schema.get("properties", {}))
    for group in schema.get("$defs", schema.get("definitions", {})).values():
        keys |= set(group.get("properties", {}))
    return keys


@pytest.mark.parametrize(("command", "flag", "rec"), RECORDS, ids=IDS)
def test_nextflow_carrier_resolves(command, flag, rec) -> None:
    carrier = rec["nextflow"]
    assert isinstance(carrier, str) and carrier, f"{command} {flag}: nextflow carrier missing"
    if carrier.startswith("params."):
        assert carrier[len("params.") :] in _schema_params(), f"{carrier} is not a schema parameter"
    else:
        assert carrier == "task.cpus" or carrier.startswith(("module:", "n/a:")), carrier


def test_short_alias_collisions_are_the_known_ones() -> None:
    seen: dict[str, set[str]] = {}
    for _command, flag, rec in RECORDS:
        for alias in rec.get("aliases", []):
            seen.setdefault(alias, set()).add(flag)
    collisions = {a: sorted(f) for a, f in seen.items() if len(f) > 1}
    assert collisions == MATRIX["known_short_alias_collisions"]


REFERENCE_PATH = ROOT / "docs" / "cli-reference.md"
AUDIT_DOCS = {"docs/audit/cli-matrix.md"}


def test_rendered_reference_in_sync() -> None:
    from scripts.render_cli_matrix import render_reference

    assert REFERENCE_PATH.read_text(encoding="utf-8") == render_reference(CLICK), (
        "run: python scripts/render_cli_matrix.py"
    )


def test_rendered_reference_has_panel_and_command_headings() -> None:
    from repgenr.cli.base import COMMAND_PANELS
    from scripts.render_cli_matrix import render_reference

    headings = set(render_reference(CLICK).splitlines())
    for panel, names in COMMAND_PANELS.items():
        assert f"## {panel}" in headings
        for name in names:
            assert f"### {name}" in headings


@pytest.mark.parametrize(("command", "flag", "rec"), RECORDS, ids=IDS)
def test_every_flag_is_documented(command, flag, rec) -> None:
    """Each flag appears in at least one document besides the audit matrix."""
    docs = set(rec["docs"]) - AUDIT_DOCS
    assert docs, f"{command} {flag} is not mentioned in any document"


def test_rendered_matrix_in_sync() -> None:
    from scripts.render_cli_matrix import render

    assert RENDERED_PATH.read_text(encoding="utf-8") == render(MATRIX), (
        "run: python scripts/render_cli_matrix.py"
    )


def test_list_defaults_render_comma_separated() -> None:
    from scripts.render_cli_matrix import _option_rows

    mini = typer.Typer()

    @mini.command()
    def one(drop: list[str] = typer.Option(["MDA", "WGA"], "--drop", help="h")) -> None:  # noqa: B008
        """x"""

    @mini.command()
    def two() -> None:
        """y"""

    cmd = typer.main.get_command(mini).commands["one"]
    row = next(r for r in _option_rows(cmd) if "`--drop`" in r)
    assert "`MDA,WGA`" in row
    assert "[" not in row
