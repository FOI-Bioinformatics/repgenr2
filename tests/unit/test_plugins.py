"""Registry + capability + preflight tests."""

from __future__ import annotations

import pytest

from repgenr.aligners.base import registry as aligners
from repgenr.annotators.base import registry as annotators
from repgenr.core.binaries import BinarySpec
from repgenr.core.errors import MissingBinaryError, PluginError
from repgenr.core.plugins import ToolCapabilities, preflight
from repgenr.coregenes.base import registry as coregenes
from repgenr.dereplicators.base import registry as dereplicators
from repgenr.snptypers.base import registry as snptypers
from repgenr.treebuilders.base import registry as treebuilders


def test_entry_points_discovered() -> None:
    assert {"drep", "skder", "galah", "sourmash"} <= set(dereplicators.names())
    assert {"progressivemauve", "cactus", "sibeliaz"} <= set(aligners.names())
    assert {"simple", "snippy", "parsnp", "ska2"} <= set(snptypers.names())
    assert {"iqtree", "fasttree", "raxmlng", "mashtree", "sourmash"} <= set(treebuilders.names())
    assert {"bakta", "prokka"} <= set(annotators.names())
    assert {"panaroo"} <= set(coregenes.names())


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
