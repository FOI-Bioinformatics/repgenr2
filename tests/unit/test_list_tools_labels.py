"""list-tools --check names where each tool runs under a container backend."""

from __future__ import annotations

import importlib
import subprocess

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core import binaries, containers
from repgenr.core.binaries import BinarySpec
from repgenr.core.plugins import Registry, ToolCapabilities, preflight

_DOCKER = ("--container", "docker")
_FAMILIES = (
    "aligners",
    "snptypers",
    "maskers",
    "treebuilders",
    "assemblers",
    "classifiers",
    "polishers",
    "annotators",
    "coregenes",
)


class _ImageTool:
    capabilities = ToolCapabilities(
        name="imagetool",
        container="quay.io/x/imagetool:1",
        required_binaries=(BinarySpec("imagetool"),),
    )

    def preflight(self):
        return preflight(self.capabilities)


class _HostTool:
    capabilities = ToolCapabilities(name="hosttool", required_binaries=(BinarySpec("hosttool"),))

    def preflight(self):
        return preflight(self.capabilities)


@pytest.fixture()
def only_two_tools(monkeypatch):
    """Every family is empty except the dereplicators, which hold the two tools."""
    for family in _FAMILIES:
        module = importlib.import_module(f"repgenr.{family}.base")
        empty = Registry(f"repgenr.test_empty_{family}")
        empty._loaded = True
        monkeypatch.setattr(module, "registry", empty)
    reg = Registry("repgenr.test_labels")
    reg._loaded = True
    reg.register("imagetool", _ImageTool)
    reg.register("hosttool", _HostTool)
    import repgenr.dereplicators.base as derep_base

    monkeypatch.setattr(derep_base, "registry", reg)
    monkeypatch.setattr(binaries.shutil, "which", lambda name, path=None: f"/usr/bin/{name}")
    monkeypatch.setattr(
        binaries,
        "_query_version",
        lambda name, args, timeout=None: "29.5.3" if name == "docker" else "1.0",
    )
    containers._ENGINE_READY.clear()
    inspected: list[list[str]] = []

    def run(argv, **kwargs):
        inspected.append(list(argv))
        # `docker info` answers; only imagetool's image is present locally.
        code = 0 if argv[1] == "info" or argv[-1] == "quay.io/x/imagetool:1" else 1
        err = "" if code == 0 else f"Error: No such image: {argv[-1]}\n"
        return subprocess.CompletedProcess(argv, code, stdout="", stderr=err)

    monkeypatch.setattr(containers.subprocess, "run", run)
    yield inspected
    containers.configure_container("none")


def _lines(output: str) -> dict[str, str]:
    return {
        ln.strip().split(":", 1)[0]: ln.strip().split(": ", 1)[1]
        for ln in output.splitlines()
        if ln.startswith("  ")
    }


def test_labels_name_the_image_or_the_host_under_docker(only_two_tools) -> None:
    result = CliRunner().invoke(app, [*_DOCKER, "list-tools", "--check"])
    assert result.exit_code == 0, result.output
    lines = _lines(result.output)
    assert lines["imagetool"] == "ok [image quay.io/x/imagetool:1] (docker 29.5.3)"
    assert lines["hosttool"] == "ok [host] (hosttool 1.0)"
    # The presence check is opt-in: no `docker image inspect` without --images.
    assert not any(argv[1:3] == ["image", "inspect"] for argv in only_two_tools)


def test_no_label_without_a_backend(only_two_tools) -> None:
    result = CliRunner().invoke(app, ["list-tools", "--check"])
    assert result.exit_code == 0, result.output
    lines = _lines(result.output)
    assert lines["imagetool"] == "ok (imagetool 1.0)"
    assert lines["hosttool"] == "ok (hosttool 1.0)"


def test_images_reports_whether_each_image_is_present(only_two_tools, monkeypatch) -> None:
    other = ToolCapabilities(name="imagetool", container="quay.io/x/imagetool:2")
    result = CliRunner().invoke(app, [*_DOCKER, "list-tools", "--check", "--images"])
    assert result.exit_code == 0, result.output
    assert _lines(result.output)["imagetool"] == (
        "ok [image quay.io/x/imagetool:1, present] (docker 29.5.3)"
    )
    monkeypatch.setattr(_ImageTool, "capabilities", other)
    result = CliRunner().invoke(app, [*_DOCKER, "list-tools", "--check", "--images"])
    assert _lines(result.output)["imagetool"] == (
        "ok [image quay.io/x/imagetool:2, not pulled] (docker 29.5.3)"
    )
    assert ["docker", "image", "inspect", "quay.io/x/imagetool:2"] in only_two_tools


def test_a_missing_host_tool_is_labelled_host(only_two_tools, monkeypatch) -> None:
    monkeypatch.setattr(
        binaries.shutil,
        "which",
        lambda name, path=None: None if name == "hosttool" else f"/usr/bin/{name}",
    )
    result = CliRunner().invoke(app, [*_DOCKER, "list-tools", "--check", "--strict"])
    assert result.exit_code == 4, result.output
    assert _lines(result.output)["hosttool"].startswith("missing [host] (hosttool: not found")


@pytest.mark.parametrize(
    ("args", "backend"),
    [(["--images"], "docker"), (["--check", "--images"], "none")],
)
def test_images_needs_check_and_a_backend(only_two_tools, args, backend) -> None:
    result = CliRunner().invoke(app, ["--container", backend, "list-tools", *args])
    assert result.exit_code == 2, result.output


def test_image_present_asks_docker_without_pulling(monkeypatch) -> None:
    calls: list[list[str]] = []

    def run(argv, **kwargs):
        calls.append(list(argv))
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr="No such image")

    monkeypatch.setattr(containers.subprocess, "run", run)
    config = containers.ContainerConfig(backend="docker", engine="podman")
    assert containers.image_present("quay.io/x/y:1", config) is False
    assert calls == [["podman", "image", "inspect", "quay.io/x/y:1"]]

    def fails(argv, **kwargs):
        raise subprocess.TimeoutExpired(argv, 30)

    monkeypatch.setattr(containers.subprocess, "run", fails)
    assert containers.image_present("quay.io/x/y:1", config) is None


def test_image_present_reads_the_singularity_cache(tmp_path) -> None:
    config = containers.ContainerConfig(backend="singularity", cache_dir=tmp_path)
    assert containers.image_present("quay.io/x/y:1", config) is False
    (tmp_path / f"{containers._sanitize('quay.io/x/y:1')}.sif").write_bytes(b"")
    assert containers.image_present("quay.io/x/y:1", config) is True
    local = tmp_path / "local.sif"
    assert containers.image_present(str(local), config) is False
    no_cache = containers.ContainerConfig(backend="singularity")
    assert containers.image_present("quay.io/x/y:1", no_cache) is None
    assert containers.image_present("quay.io/x/y:1", containers.ContainerConfig()) is None


def test_images_reports_an_adapters_secondary_images(only_two_tools, monkeypatch) -> None:
    # racon records minimap2's own image beside its own; both are reported.
    def with_helper(self):
        versions = preflight(self.capabilities)
        versions["helper"] = "quay.io/x/helper:3"
        return versions

    monkeypatch.setattr(_ImageTool, "preflight", with_helper)
    result = CliRunner().invoke(app, [*_DOCKER, "list-tools", "--check", "--images"])
    assert result.exit_code == 0, result.output
    assert _lines(result.output)["imagetool"] == (
        "ok [image quay.io/x/imagetool:1, present; helper image quay.io/x/helper:3, "
        "not pulled] (docker 29.5.3)"
    )


def test_images_skips_the_presence_query_for_a_missing_tool(only_two_tools, monkeypatch) -> None:
    from repgenr.core.errors import MissingBinaryError

    def missing(self):
        raise MissingBinaryError("imagetool: not available")

    monkeypatch.setattr(_ImageTool, "preflight", missing)
    result = CliRunner().invoke(app, [*_DOCKER, "list-tools", "--check", "--images"])
    assert _lines(result.output)["imagetool"] == (
        "missing [image quay.io/x/imagetool:1] (imagetool: not available)"
    )
    assert not any(argv[1:3] == ["image", "inspect"] for argv in only_two_tools)


@pytest.mark.parametrize(
    ("stderr", "expected"),
    [
        ("Error: No such image: quay.io/x/y:1\n", False),
        ("Error: quay.io/x/y:1: image not known\n", False),
        ("Cannot connect to the Docker daemon at unix:///var/run/docker.sock.\n", None),
        ("", None),
    ],
)
def test_image_present_is_unknown_unless_the_engine_says_absent(
    monkeypatch, stderr, expected
) -> None:
    # With the daemon down every image read "not pulled".
    def run(argv, **kwargs):
        return subprocess.CompletedProcess(argv, 1, stdout="", stderr=stderr)

    monkeypatch.setattr(containers.subprocess, "run", run)
    config = containers.ContainerConfig(backend="docker")
    assert containers.image_present("quay.io/x/y:1", config) is expected
