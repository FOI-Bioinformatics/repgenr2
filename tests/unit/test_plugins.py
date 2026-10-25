"""Registry + capability + preflight tests."""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

from repgenr.aligners.base import registry as aligners
from repgenr.annotators.base import registry as annotators
from repgenr.assemblers.base import registry as assemblers
from repgenr.classifiers.base import registry as classifiers
from repgenr.core.binaries import BinarySpec
from repgenr.core.errors import MissingBinaryError, PluginError
from repgenr.core.plugins import ToolCapabilities, preflight
from repgenr.coregenes.base import registry as coregenes
from repgenr.dereplicators.base import registry as dereplicators
from repgenr.maskers.base import registry as maskers
from repgenr.polishers.base import registry as polishers
from repgenr.snptypers.base import registry as snptypers
from repgenr.treebuilders.base import registry as treebuilders

ROOT = Path(__file__).resolve().parents[2]

REGISTRIES = {
    "repgenr.dereplicators": dereplicators,
    "repgenr.aligners": aligners,
    "repgenr.snptypers": snptypers,
    "repgenr.annotators": annotators,
    "repgenr.coregenes": coregenes,
    "repgenr.maskers": maskers,
    "repgenr.assemblers": assemblers,
    "repgenr.classifiers": classifiers,
    "repgenr.polishers": polishers,
    "repgenr.treebuilders": treebuilders,
}


def _declared_entry_points() -> dict[str, set[str]]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    groups = data["project"]["entry-points"]
    return {group: set(names) for group, names in groups.items()}


def test_every_family_has_a_registry() -> None:
    assert set(_declared_entry_points()) == set(REGISTRIES)


@pytest.mark.parametrize("group", sorted(REGISTRIES))
def test_entry_points_discovered(group: str) -> None:
    """Every adapter declared in pyproject.toml is discovered by its registry."""
    declared = _declared_entry_points()[group]
    assert declared <= set(REGISTRIES[group].names()), group


def test_unknown_tool_raises() -> None:
    with pytest.raises(PluginError):
        dereplicators.get("does-not-exist")


def test_preflight_missing_binary() -> None:
    caps = ToolCapabilities(
        name="phantom",
        required_binaries=(BinarySpec("definitely-not-a-real-binary-xyz"),),
    )
    with pytest.raises(MissingBinaryError):
        preflight(caps)
