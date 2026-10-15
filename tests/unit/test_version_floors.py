"""Every built-in adapter's binaries carry a version floor, except the tools that
print no parseable version (their floor would never be checked)."""

from __future__ import annotations

import pytest

from repgenr.aligners.base import registry as aligners
from repgenr.annotators.base import registry as annotators
from repgenr.assemblers.base import registry as assemblers
from repgenr.classifiers.base import registry as classifiers
from repgenr.dereplicators.base import registry as dereplicators
from repgenr.maskers.base import registry as maskers
from repgenr.polishers.base import registry as polishers
from repgenr.snptypers.base import registry as snptypers
from repgenr.treebuilders.base import registry as treebuilders

# Binaries that answer no version flag with a dotted version: a floor on them
# could never be enforced, so they are exempt and stay lenient.
_NO_PARSEABLE_VERSION = {"progressiveMauve", "sibeliaz", "hal2maf"}


def _specs():
    for reg in (
        dereplicators,
        aligners,
        snptypers,
        maskers,
        treebuilders,
        assemblers,
        classifiers,
        polishers,
        annotators,
    ):
        for name in reg.names():
            for spec in reg.get(name).capabilities.required_binaries:
                yield f"{name}:{spec.name}", spec


@pytest.mark.parametrize(("label", "spec"), list(_specs()), ids=[lab for lab, _ in _specs()])
def test_binary_declares_a_version_floor(label, spec) -> None:
    if spec.name in _NO_PARSEABLE_VERSION:
        assert spec.min_version is None, f"{label} cannot be version-checked; drop the floor"
        return
    assert spec.min_version, f"{label} declares no min_version"
