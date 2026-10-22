"""Render tests/audit/cli_matrix.yaml to docs/audit/cli-matrix.md.

Usage: python scripts/render_cli_matrix.py   (from the repository root)
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / "tests" / "audit" / "cli_matrix.yaml"
OUT = ROOT / "docs" / "audit" / "cli-matrix.md"
REFERENCE = ROOT / "docs" / "cli-reference.md"
SKIP_OPTS = {"--help", "--install-completion", "--show-completion"}


def _cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return str(value).replace("|", "\\|")


def _table(flags: dict) -> list[str]:
    rows = [
        "| flag | aliases | param | validated | nextflow | live | docs |",
        "|---|---|---|---|---|---|---|",
    ]
    for flag, rec in flags.items():
        rows.append(
            f"| `{flag}` | {_cell(rec.get('aliases'))} | {_cell(rec['param'])} | "
            f"{_cell(rec['validated'])} | {_cell(rec.get('nextflow'))} | {_cell(rec['live'])} | "
            f"{_cell(rec['docs'])} |"
        )
    return rows


def render(matrix: dict) -> str:
    n_flags = len(matrix["global_flags"]) + sum(
        len(c["flags"]) for c in matrix["commands"].values()
    )
    todo = sum(
        1
        for c in matrix["commands"].values()
        for r in c["flags"].values()
        if str(r["live"]).startswith("todo:")
    ) + sum(1 for r in matrix["global_flags"].values() if str(r["live"]).startswith("todo:"))
    lines = [
        "# CLI matrix",
        "",
        "Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;",
        "`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.",
        "",
        f"{len(matrix['commands'])} commands, {n_flags} flags "
        f"({n_flags - todo} with a live test or an n/a reason, {todo} pending).",
        "",
        "## Global flags",
        "",
        *_table(matrix["global_flags"]),
        "",
    ]
    for name, cmd in matrix["commands"].items():
        lines += [f"## {name}", "", f"dispatch: `{cmd['dispatch']}`", ""]
        lines += _table(cmd["flags"]) if cmd["flags"] else ["(no flags)"]
        lines.append("")
    lines += [
        "## Short-alias collisions",
        "",
        *(
            f"- `{a}`: {', '.join(f'`{f}`' for f in fl)}"
            for a, fl in matrix["known_short_alias_collisions"].items()
        ),
        "",
    ]
    return "\n".join(lines)


def _click_tree():
    import typer

    from repgenr.cli.main import app

    return typer.main.get_command(app)


def _format_default(value) -> str:
    """A default as table text; list defaults are comma-separated."""
    if value in (None, [], ()):
        return ""
    if isinstance(value, (list, tuple)):
        return f"`{','.join(str(v) for v in value)}`"
    return f"`{value}`"


def _option_rows(cmd) -> list[str]:
    rows = ["| option | default | description |", "|---|---|---|"]
    for p in cmd.params:
        opts = getattr(p, "opts", None)
        if not opts or opts[0] in SKIP_OPTS:
            continue
        names = ", ".join(f"`{o}`" for o in [*opts, *(getattr(p, "secondary_opts", None) or [])])
        if getattr(p, "required", False):
            default = "required"
        elif getattr(p, "is_flag", False):
            # A --x/--no-x pair without a default decides at run time.
            if p.default is None and getattr(p, "secondary_opts", None):
                default = "auto"
            else:
                default = "off" if not p.default else "on"
        else:
            default = _format_default(p.default)
        rows.append(f"| {names} | {default} | {_cell(p.help or '')} |")
    return rows


def render_reference(cli=None) -> str:
    """The command reference: every command with its options, defaults and help."""
    from repgenr.cli.base import COMMAND_PANELS

    cli = cli or _click_tree()
    lines = [
        "# Command reference",
        "",
        "This page is generated from the command tree by `scripts/render_cli_matrix.py`,",
        "and a unit test keeps it in sync. Each command's text is the same as",
        "`repgenr <command> --help` prints, and the commands are grouped as in",
        "`repgenr --help`. [usage.md](usage.md) explains the stages, and",
        "[choosing-tools.md](choosing-tools.md) compares the tools.",
        "",
        "The default column reads as follows:",
        "",
        "- A value is the declared default.",
        "- A blank cell means the option is unset; when a rule decides the value at run",
        '  time, the help ends with "Default: <rule>."',
        "- `auto` marks an on/off pair (`--x/--no-x`) decided at run time.",
        "- `on` and `off` give the state of a flag.",
        "- `required` marks an option that must be given.",
        "",
        "Global options go before the command name",
        "(`repgenr --container docker dereplicate ...`).",
        "",
        'The "Nextflow data-channel steps" panel lists the steps that the Nextflow',
        "pipeline calls, one task each. They are not the commands most users run.",
        "",
        "## Global options",
        "",
        *_option_rows(cli),
        "",
    ]
    for panel, names in COMMAND_PANELS.items():
        lines += [f"## {panel}", ""]
        for name in names:
            cmd = cli.commands[name]
            lines += [f"### {name}", "", (cmd.help or "").strip(), ""]
            lines += (
                _option_rows(cmd)
                if any(getattr(p, "opts", None) for p in cmd.params)
                else ["(no options)"]
            )
            lines.append("")
    return "\n".join(lines)


def main() -> int:
    with open(MATRIX, encoding="utf-8") as fh:
        matrix = yaml.safe_load(fh)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render(matrix), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    REFERENCE.write_text(render_reference(), encoding="utf-8")
    print(f"wrote {REFERENCE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
