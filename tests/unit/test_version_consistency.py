"""pyproject, the Nextflow manifest, the changelog and the installed distribution agree."""

from __future__ import annotations

import re
import tomllib
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _pyproject_version() -> str:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return data["project"]["version"]


def test_nextflow_manifest_matches_pyproject() -> None:
    config = (ROOT / "nextflow" / "nextflow.config").read_text(encoding="utf-8")
    match = re.search(r"^\s*version\s*=\s*'([^']+)'", config, flags=re.MULTILINE)
    assert match is not None, "nextflow.config has no manifest version"
    assert match.group(1) == _pyproject_version()


def test_installed_package_matches_pyproject() -> None:
    assert version("repgenr2") == _pyproject_version()


def test_changelog_has_section_for_version() -> None:
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert f"## [{_pyproject_version()}]" in changelog
