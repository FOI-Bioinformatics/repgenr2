"""The conda environment files under envs/ agree with the adapters.

Offline checks only: the solve itself runs in .github/workflows/envs.yml.
Every version floor an adapter enforces at preflight (BinarySpec.min_version)
must appear as the same ``>=`` floor in exactly one envs/*.yml file, so an
environment built from these files passes preflight.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest
import yaml

from repgenr.core.plugins import ToolCapabilities

ENVS = Path(__file__).resolve().parents[1] / "envs"
EXPECTED_FILES = {
    "core",
    "gubbins",
    "mashtree",
    "snippy",
    "parsnp",
    "checkm2",
    "mauve",
    "bakta",
    "prokka",
}
FAMILIES = (
    "aligners",
    "annotators",
    "assemblers",
    "classifiers",
    "dereplicators",
    "maskers",
    "polishers",
    "snptypers",
    "treebuilders",
)
# Binaries with no conda package: cactus is distributed as a container.
CONTAINER_ONLY = {"cactus-pangenome", "hal2maf"}
WORKFLOW = ENVS.parent / ".github" / "workflows" / "envs.yml"
# Packages allowed in more than one file. Empty: a package in two files means
# two environments may put different builds of one tool on PATH.
SHARED_ALLOWED: set[str] = set()

_SPEC_RE = re.compile(r"^([A-Za-z0-9_.\-]+)\s*(>=|==|=)?\s*([^\s#]+)?")


def _load(name: str) -> dict:
    return yaml.safe_load((ENVS / f"{name}.yml").read_text())


def _specs(name: str) -> dict[str, tuple[str | None, str | None]]:
    """Package name -> (operator, version) for the conda entries of one file."""
    out: dict[str, tuple[str | None, str | None]] = {}
    for dep in _load(name)["dependencies"]:
        if not isinstance(dep, str):  # the pip: sub-list
            continue
        match = _SPEC_RE.match(dep.strip())
        assert match, f"{name}.yml: cannot parse {dep!r}"
        out[match.group(1).lower()] = (match.group(2), match.group(3))
    return out


def _env_files() -> list[Path]:
    # macOS writes ._* AppleDouble files beside the real ones on some volumes.
    return sorted(p for p in ENVS.glob("*.yml") if not p.name.startswith("._"))


def _all_specs() -> dict[str, dict[str, tuple[str | None, str | None]]]:
    return {path.stem: _specs(path.stem) for path in _env_files()}


def _capabilities() -> list[ToolCapabilities]:
    caps: list[ToolCapabilities] = []
    for family in FAMILIES:
        registry = importlib.import_module(f"repgenr.{family}.base").registry
        for name in registry.names():
            if not registry.is_broken(name):
                caps.append(registry.get(name).capabilities)
    from repgenr.stages.assemble_qc import CHECKM2_CAPS
    from repgenr.stages.genome import DATASETS_CAPS

    caps.extend([CHECKM2_CAPS, DATASETS_CAPS])
    return caps


def _package_for(binary: str, caps: ToolCapabilities, available: set[str] | None = None) -> str:
    """The conda package that ships ``binary``, from the adapter's conda spec.

    A package whose name contains the binary name or is contained in it
    (gubbins / run_gubbins.py, snippy / snippy-core, ska2 / ska). Otherwise
    the binary's own name when that is a package in ``available`` (skani
    under skder), else the adapter's first package (progressiveMauve from
    mauvealigner).
    """
    lowered = binary.lower()
    packages = [spec.split("::")[-1].split("=")[0].lower() for spec in caps.conda]
    for package in packages:
        if package == lowered:
            return package
    for package in packages:
        if package in lowered or lowered in package:
            return package
    if available is not None and lowered not in available and packages:
        return packages[0]
    return lowered


def test_the_expected_files_exist() -> None:
    assert {path.stem for path in _env_files()} == EXPECTED_FILES
    assert not (ENVS.parent / "environment.yml").exists()


@pytest.mark.parametrize("name", sorted(EXPECTED_FILES))
def test_each_file_names_its_environment_and_channels(name: str) -> None:
    data = _load(name)
    expected = "repgenr" if name == "core" else f"repgenr-{name}"
    assert data["name"] == expected
    # nodefaults: a configured 'defaults' channel (Miniconda, Anaconda) with
    # flexible priority made the core solve run for more than 17 minutes.
    assert data["channels"] == ["conda-forge", "bioconda", "nodefaults"]
    first_line = (ENVS / f"{name}.yml").read_text().splitlines()[0]
    assert first_line.startswith("#") and "Platforms:" in first_line


def test_core_pins_python_312_and_installs_the_package() -> None:
    assert _specs("core")["python"] == ("=", "3.12")
    pip = [dep for dep in _load("core")["dependencies"] if isinstance(dep, dict)]
    assert pip == [{"pip": ["-e .."]}]


def test_no_package_is_in_two_files() -> None:
    seen: dict[str, str] = {}
    for env, specs in _all_specs().items():
        for package in specs:
            if package in SHARED_ALLOWED:
                continue
            assert package not in seen, f"{package} is in {seen[package]}.yml and {env}.yml"
            seen[package] = env


def test_every_preflight_floor_is_in_an_env_file() -> None:
    by_package = {pkg: spec for specs in _all_specs().values() for pkg, spec in specs.items()}
    checked = 0
    for caps in _capabilities():
        for spec in caps.required_binaries:
            if spec.min_version is None or spec.name in CONTAINER_ONLY:
                continue
            package = _package_for(spec.name, caps)
            assert package in by_package, f"{caps.name}: {package} is in no envs/*.yml file"
            assert by_package[package] == (">=", spec.min_version), (
                f"{caps.name}: {package} floor {by_package[package]} differs from "
                f"BinarySpec {spec.name} min_version {spec.min_version}"
            )
            checked += 1
    assert checked >= 20  # guards against the registries silently coming back empty


def test_package_mapping_follows_the_conda_spec() -> None:
    caps = ToolCapabilities(name="x", conda=("bioconda::gubbins",))
    assert _package_for("run_gubbins.py", caps) == "gubbins"
    assert _package_for("skani", ToolCapabilities(name="skder", conda=("bioconda::skder",))) == (
        "skani"
    )
    assert _package_for("ska", ToolCapabilities(name="ska2", conda=("bioconda::ska2",))) == "ska2"
    mauve = ToolCapabilities(name="m", conda=("bioconda::mauvealigner", "boost-cpp=1.74"))
    assert _package_for("progressiveMauve", mauve, {"mauvealigner"}) == "mauvealigner"


def test_every_required_binary_comes_from_an_env_file() -> None:
    available = {pkg for specs in _all_specs().values() for pkg in specs}
    checked = 0
    for caps in _capabilities():
        for spec in caps.required_binaries:
            if spec.name in CONTAINER_ONLY:
                continue
            package = _package_for(spec.name, caps, available)
            assert package in available, f"{caps.name}: {spec.name} ({package}) is in no envs/*.yml"
            checked += 1
    assert checked >= 25


def test_the_adapter_conda_specs_are_in_the_env_files() -> None:
    # --wave builds an image from caps.conda; the env files install the same packages.
    available = {pkg for specs in _all_specs().values() for pkg in specs}
    for caps in _capabilities():
        for spec in caps.conda:
            package = spec.split("::")[-1].split("=")[0].lower()
            assert package in available, f"{caps.name}: {spec} is in no envs/*.yml"


def test_the_ci_matrix_solves_every_file() -> None:
    matrix = yaml.safe_load(WORKFLOW.read_text())["jobs"]["solve"]["strategy"]["matrix"]
    assert set(matrix["env"]) == EXPECTED_FILES
    assert matrix["subdir"] == ["linux-64"]
    assert {"env": "core", "subdir": "osx-arm64"} in matrix["include"]
