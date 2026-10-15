"""Tests for capability-based auto tool-selection and scale warnings."""

from __future__ import annotations

from repgenr.core.binaries import BinarySpec
from repgenr.core.plugins import Registry, ToolCapabilities, auto_select, scale_warning


class _Tool:
    capabilities: ToolCapabilities

    def __init__(self):
        pass


def _make_registry(specs: dict[str, ToolCapabilities]) -> Registry:
    reg: Registry = Registry("test.group")
    reg._loaded = True  # bare registry: skip entry-point discovery
    for name, caps in specs.items():
        cls = type(f"Tool_{name}", (_Tool,), {"capabilities": caps})
        reg.register(name, cls, replace=True)
    return reg


def _caps(name, limit, native=False, binaries=()):
    return ToolCapabilities(
        name=name,
        recommended_max_genomes=limit,
        supports_native_scaling=native,
        required_binaries=binaries,
    )


def test_auto_prefers_unbounded_native_tool() -> None:
    reg = _make_registry(
        {
            "drep": _caps("drep", 2000, native=False),
            "skder": _caps("skder", None, native=True),
            "sourmash": _caps("sourmash", None, native=True),
        }
    )
    # 5000 genomes: bounded drep does not fit; unbounded native tools win;
    # tie broken alphabetically -> skder over sourmash
    assert auto_select(reg, 5000) == "skder"


def test_auto_picks_tightest_fitting_tool() -> None:
    reg = _make_registry(
        {
            "iqtree": _caps("iqtree", 500),
            "fasttree": _caps("fasttree", 5000),
        }
    )
    # both fit at 300: the TIGHTEST limit wins (higher-quality tools declare
    # tighter recommended scales); at 3000 only fasttree fits
    assert auto_select(reg, 300) == "iqtree"
    assert auto_select(reg, 3000) == "fasttree"


def test_auto_prefers_bounded_fit_over_unbounded() -> None:
    reg = _make_registry(
        {
            "mashtree": _caps("mashtree", 10000),
            "iqtree": _caps("iqtree", 500),
            "skder": _caps("skder", None, native=True),
        }
    )
    assert auto_select(reg, 100) == "iqtree"  # tightest fitting
    assert auto_select(reg, 8000) == "mashtree"  # iqtree out; bounded beats unbounded
    assert auto_select(reg, 50000) == "skder"  # only unbounded accommodates


def test_auto_tiebreak_follows_documented_defaults() -> None:
    # galah sorts before skder alphabetically; both unbounded -> the documented
    # default (skder) must win the tie.
    reg = _make_registry(
        {
            "galah": _caps("galah", None, native=True),
            "skder": _caps("skder", None, native=True),
        }
    )
    assert auto_select(reg, 100) == "skder"


def test_auto_availability_is_container_aware(monkeypatch) -> None:
    """Under an active container backend, a tool with a pinned image is
    available even when nothing is installed on the host."""
    from repgenr.core import containers, plugins

    reg = _make_registry(
        {
            "hosttool": ToolCapabilities(
                name="hosttool",
                required_binaries=(BinarySpec("python"),),
            ),
            "imagetool": ToolCapabilities(
                name="imagetool",
                required_binaries=(BinarySpec("definitely-missing-xyz"),),
                container="quay.io/x/imagetool:1",
                recommended_max_genomes=10,  # tighter fit than hosttool
            ),
        }
    )
    try:
        containers.configure_container("docker")
        assert plugins.auto_select(reg, 5) == "imagetool"
    finally:
        containers.configure_container("none")
    # natively, the missing binary disqualifies imagetool
    assert plugins.auto_select(reg, 5) == "hosttool"


def test_conda_only_tool_is_not_available_in_a_container_without_wave() -> None:
    """A conda spec yields an image only through --wave. Without it the tool
    runs on the host (run_tool warns and falls back), so availability must be
    the host PATH: skDER was auto-selected for >2000 genomes under --container
    docker and then exited 4 because skder was not installed."""
    from repgenr.core import containers, plugins

    reg = _make_registry(
        {
            "condaonly": ToolCapabilities(
                name="condaonly",
                required_binaries=(BinarySpec("definitely-missing-xyz"),),
                conda=("bioconda::condaonly",),
            ),
            "imagetool": ToolCapabilities(
                name="imagetool",
                required_binaries=(BinarySpec("definitely-missing-xyz"),),
                container="quay.io/x/imagetool:1",
            ),
        }
    )
    try:
        containers.configure_container("docker")
        assert not plugins.tool_available(reg.get("condaonly").capabilities)
        assert plugins.tool_available(reg.get("imagetool").capabilities)
        assert plugins.auto_select(reg, 5) == "imagetool"
        containers.configure_container("docker", wave_enabled=True)
        assert plugins.tool_available(reg.get("condaonly").capabilities)
    finally:
        containers.configure_container("none")


def test_imageless_tool_on_the_host_is_available_in_a_container() -> None:
    from repgenr.core import containers, plugins

    caps = ToolCapabilities(name="hostonly", required_binaries=(BinarySpec("python"),))
    try:
        containers.configure_container("docker")
        assert plugins.tool_available(caps)
    finally:
        containers.configure_container("none")


def test_auto_skips_tools_with_missing_binaries() -> None:
    # 'aaa' sorts first but needs a missing binary; 'zzz' is installed (python).
    reg = _make_registry(
        {
            "aaa": _caps(
                "aaa", None, native=True, binaries=(BinarySpec("definitely-missing-xyz"),)
            ),
            "zzz": _caps("zzz", None, native=True, binaries=(BinarySpec("python"),)),
        }
    )
    assert auto_select(reg, 100) == "zzz"


def test_scale_warning_triggers_over_limit() -> None:
    reg = _make_registry(
        {
            "drep": _caps("drep", 2000),
            "skder": _caps("skder", None, native=True),
        }
    )
    warn = scale_warning(reg, "drep", 5000)
    assert warn is not None
    limit, alts = warn
    assert limit == 2000
    assert "skder" in alts


def test_scale_warning_none_within_limit() -> None:
    reg = _make_registry({"drep": _caps("drep", 2000)})
    assert scale_warning(reg, "drep", 100) is None


def test_scale_warning_none_for_unbounded() -> None:
    reg = _make_registry({"skder": _caps("skder", None, native=True)})
    assert scale_warning(reg, "skder", 100000) is None


def test_allowed_restricts_the_choice_even_past_a_better_fit() -> None:
    reg = _make_registry(
        {
            "iqtree": _caps("iqtree", 500),
            "fasttree": _caps("fasttree", 5000),
            "mashtree": _caps("mashtree", 10000),
        }
    )
    assert auto_select(reg, 300) == "iqtree"
    assert auto_select(reg, 300, allowed={"fasttree", "mashtree"}) == "fasttree"
    assert auto_select(reg, 8000, allowed={"iqtree", "fasttree"}) == "fasttree"
    assert auto_select(reg, 300, allowed=set()) is None
