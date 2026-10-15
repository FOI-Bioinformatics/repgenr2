"""Broken plugins are visible: load warnings, auto-select logging, listings (WP1-6)."""

from __future__ import annotations

import logging
import re

import pytest

from repgenr.core import plugins
from repgenr.core.errors import MissingBinaryError, PluginError
from repgenr.core.plugins import Registry, ToolCapabilities

_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


class _GoodAdapter:
    capabilities = ToolCapabilities(name="goodtool")


class _GoodEP:
    name = "goodtool"

    def load(self):
        return _GoodAdapter


class _BadEP:
    name = "badtool"

    def load(self):
        raise ImportError("boom: missing dependency")


@pytest.fixture()
def registry(monkeypatch) -> Registry:
    # Scope the fake entry points to this test group only: a test that loads
    # the REAL registries while this patch is active (e.g. list-tools) must not
    # cache fake adapters into them for the rest of the pytest process.
    real_entry_points = plugins.entry_points

    def fake_entry_points(group):
        if group == "repgenr.test_tools":
            return [_GoodEP(), _BadEP()]
        return real_entry_points(group=group)

    monkeypatch.setattr(plugins, "entry_points", fake_entry_points)
    return Registry("repgenr.test_tools")


def test_broken_plugin_import_logged_at_warning(registry, caplog) -> None:
    with caplog.at_level(logging.WARNING, logger="repgenr"):
        names = registry.names()
    assert names == ["badtool", "goodtool"]
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert any("badtool" in r.getMessage() for r in warnings)


def test_registry_reports_broken_names(registry) -> None:
    registry.names()  # trigger load
    assert registry.is_broken("badtool")
    assert not registry.is_broken("goodtool")
    assert not registry.is_broken("nonexistent")


def test_get_broken_plugin_still_raises(registry) -> None:
    with pytest.raises(PluginError, match="badtool"):
        registry.get("badtool")


def test_auto_select_logs_skipped_candidates(registry, caplog) -> None:
    registry.names()
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="repgenr"):
        chosen = plugins.auto_select(registry, 10)
    assert chosen == "goodtool"
    skipped = [r for r in caplog.records if "badtool" in r.getMessage()]
    assert skipped, "auto_select should log the skipped broken candidate"


def test_auto_select_warns_once_per_broken_candidate(registry, caplog) -> None:
    # A stage selects twice (precheck, then run); the skipped broken plugin
    # was warned about twice in one dereplicate run.
    registry.names()
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="repgenr"):
        plugins.auto_select(registry, 10)
        plugins.auto_select(registry, 10)
    skipped = [r for r in caplog.records if "auto-select skipping 'badtool'" in r.getMessage()]
    assert len(skipped) == 1


def test_list_tools_check_shows_why_a_plugin_is_broken(registry, monkeypatch) -> None:
    from typer.testing import CliRunner

    import repgenr.dereplicators.base as derep_base
    from repgenr.cli.main import app

    monkeypatch.setattr(derep_base, "registry", registry)
    result = CliRunner().invoke(app, ["list-tools", "--check"])
    assert result.exit_code == 0
    assert "  badtool: broken (failed to load: boom: missing dependency)" in result.output
    # An adapter whose preflight raises something unexpected (here it has
    # none) is reported on its own line instead of ending the listing.
    assert "  goodtool: error (AttributeError:" in result.output
    assert "polishers:" in result.output


def test_list_tools_marks_broken_entries(registry, monkeypatch) -> None:
    from typer.testing import CliRunner

    import repgenr.dereplicators.base as derep_base
    from repgenr.cli.main import app

    monkeypatch.setattr(derep_base, "registry", registry)
    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    assert "badtool (broken)" in result.output
    assert "goodtool" in result.output


def test_register_and_unregister_roundtrip(register_tool) -> None:
    from repgenr.core.errors import PluginError
    from repgenr.core.plugins import Registry, ToolCapabilities

    reg: Registry = Registry("repgenr.test_register")
    reg._loaded = True  # bare registry: skip entry-point discovery

    class Adapter:
        capabilities = ToolCapabilities(name="mytool")

    register_tool(reg, "mytool", Adapter)
    assert "mytool" in reg.names()
    assert reg.get("mytool") is Adapter
    with pytest.raises(PluginError, match="already registered"):
        reg.register("mytool", Adapter)
    reg.register("mytool", Adapter, replace=True)  # explicit override allowed


def test_third_party_adapter_without_run_tool_survives_contract_fixture(
    register_tool, monkeypatch, tmp_path
) -> None:
    """An external adapter module lacking a module-level run_tool symbol must
    not break the contract suite's patching (raising=False)."""
    import types

    from repgenr.dereplicators.base import Dereplicator, registry

    mod = types.ModuleType("thirdparty_derep")

    class ThirdParty(Dereplicator):
        from repgenr.core.plugins import ToolCapabilities

        capabilities = ToolCapabilities(name="thirdparty")

        def dereplicate(self, genomes, out_dir, params, logger):
            raise NotImplementedError

    ThirdParty.__module__ = "thirdparty_derep"
    import sys

    monkeypatch.setitem(sys.modules, "thirdparty_derep", mod)
    register_tool(registry, "thirdparty", ThirdParty)

    # mimic the contract fixture's patch loop over every registered adapter
    for name in registry.names():
        module = sys.modules[registry.get(name).__module__]
        monkeypatch.setattr(module, "run_tool", lambda *a, **k: 0, raising=False)


def test_glance_rejects_tool_without_compare(tmp_path) -> None:
    from repgenr.core.context import WorkdirContext
    from repgenr.core.errors import UserInputError
    from repgenr.stages.glance import GlanceParams
    from repgenr.stages.glance import run as glance_run

    ctx = WorkdirContext(tmp_path, create=True)
    ctx.genomes_dir.mkdir()
    (ctx.genomes_dir / "g.fasta").write_text(">g\nAC\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="drep"):
        glance_run(ctx, GlanceParams(tool="skder"))  # skder has no compare()


def test_outgroup_builder_must_support_distance_matrix() -> None:
    from repgenr.core.errors import UserInputError
    from repgenr.viral._outgroup import resolve_outgroup_builder

    builder = resolve_outgroup_builder("mashtree")
    assert builder.capabilities.name == "mashtree"
    with pytest.raises(UserInputError, match="mashtree"):
        resolve_outgroup_builder("iqtree")  # no distance_matrix support


def test_masker_family_is_registered() -> None:
    from repgenr.maskers.base import registry as masker_registry

    assert "gubbins" in masker_registry.names()


def test_snptype_dispatches_masker_via_registry(tmp_path, monkeypatch, register_tool) -> None:
    """The snptype stage's --mask goes through the masker registry, so a new
    masker is one adapter + one entry point away."""

    from repgenr.core.plugins import ToolCapabilities
    from repgenr.maskers.base import Masker, registry
    from repgenr.snptypers.base import SnpResult
    from repgenr.stages import snptype as snptype_mod
    from repgenr.stages.snptype import SnptypeParams, snptype_core

    class FakeMasker(Masker):
        capabilities = ToolCapabilities(name="fakemask")

        def mask(self, full_alignment, out_dir, params, logger):
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / "filtered.fasta"
            out.write_text(">m\nAC\n", encoding="utf-8")
            return out

    register_tool(registry, "fakemask", FakeMasker)

    class FakeTyper:
        requires_reference = False
        capabilities = ToolCapabilities(name="faketyper")

        def preflight(self):
            return {}

        def call(self, genomes, ref, out_dir, params, logger):
            core = out_dir / "raw.fasta"
            full = out_dir / "full.fasta"
            core.write_text(">a\nACGT\n", encoding="utf-8")
            full.write_text(">a\nACGTACGT\n", encoding="utf-8")
            return SnpResult(core_snp_fasta=core, full_alignment=full)

    monkeypatch.setattr(snptype_mod.snp_registry, "create", lambda name: FakeTyper())
    genomes = [tmp_path / "g1.fasta"]
    genomes[0].write_text(">g1\nACGT\n", encoding="utf-8")
    result, versions = snptype_core(
        genomes,
        None,
        tmp_path / "snp",
        tmp_path / "scratch",
        SnptypeParams(tool="whatever", mask="fakemask"),
        logging.getLogger("test"),
    )
    assert result.masked is True
    assert (tmp_path / "snp" / "core_snp.fasta").read_text(encoding="utf-8") == ">m\nAC\n"


def test_unknown_masker_lists_available(tmp_path, monkeypatch) -> None:
    from repgenr.core.errors import PluginError
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.snptypers.base import SnpResult
    from repgenr.stages import snptype as snptype_mod
    from repgenr.stages.snptype import SnptypeParams, snptype_core

    class FakeTyper:
        requires_reference = False
        capabilities = ToolCapabilities(name="faketyper")

        def preflight(self):
            return {}

        def call(self, genomes, ref, out_dir, params, logger):
            core = out_dir / "raw.fasta"
            core.write_text(">a\nACGT\n", encoding="utf-8")
            return SnpResult(core_snp_fasta=core)

    monkeypatch.setattr(snptype_mod.snp_registry, "create", lambda name: FakeTyper())
    genomes = [tmp_path / "g1.fasta"]
    genomes[0].write_text(">g1\nACGT\n", encoding="utf-8")
    with pytest.raises(PluginError, match="gubbins"):
        snptype_core(
            genomes,
            None,
            tmp_path / "snp",
            tmp_path / "scratch",
            SnptypeParams(tool="x", mask="nosuchmask"),
            logging.getLogger("test"),
        )


def test_snptype_core_hands_the_mask_exclusion_to_the_masker(tmp_path, register_tool) -> None:
    """D-10: SnptypeParams.mask_exclude reaches MaskParams.exclude."""
    import logging

    from repgenr.maskers.base import Masker, MaskParams
    from repgenr.maskers.base import registry as masker_registry
    from repgenr.snptypers.base import SnpResult, SnpTyper
    from repgenr.snptypers.base import registry as snp_registry
    from repgenr.stages.snptype import SnptypeParams, snptype_core

    seen: dict = {}

    class _Typer(SnpTyper):
        capabilities = ToolCapabilities(name="excltyper")
        requires_reference = False

        def preflight(self):
            return {}

        def call(self, genomes, reference, out_dir, params, logger):
            out_dir.mkdir(parents=True, exist_ok=True)
            core = out_dir / "core.fasta"
            core.write_text(">a\nA\n", encoding="utf-8")
            full = out_dir / "full.fasta"
            full.write_text(">a\nAC\n>og\nAT\n", encoding="utf-8")
            return SnpResult(core_snp_fasta=core, full_alignment=full)

    class _Masker(Masker):
        capabilities = ToolCapabilities(name="exclmasker")

        def preflight(self):
            return {}

        def mask(self, full_alignment, out_dir, params: MaskParams, logger):
            seen["exclude"] = params.exclude
            out_dir.mkdir(parents=True, exist_ok=True)
            out = out_dir / "masked.fasta"
            out.write_text(full_alignment.read_text(encoding="utf-8"), encoding="utf-8")
            return out

    register_tool(snp_registry, "excltyper", _Typer)
    register_tool(masker_registry, "exclmasker", _Masker)
    g = tmp_path / "a.fasta"
    g.write_text(">a\nACGT\n", encoding="utf-8")
    snptype_core(
        [g],
        None,
        tmp_path / "snp",
        tmp_path / "scratch",
        SnptypeParams(tool="excltyper", mask="exclmasker", mask_exclude=("og",)),
        logging.getLogger("t"),
    )
    assert seen["exclude"] == frozenset({"og"})


def test_list_tools_check_reports_each_tool(register_tool) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core.binaries import BinarySpec
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.dereplicators.base import Dereplicator
    from repgenr.dereplicators.base import registry as derep_registry

    class Present(Dereplicator):
        capabilities = ToolCapabilities(name="presenttool")

        def preflight(self):
            return {"presenttool": "9.9"}

        def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
            raise NotImplementedError

    class Absent(Dereplicator):
        capabilities = ToolCapabilities(
            name="absenttool",
            required_binaries=(BinarySpec("repgenr-no-such-binary-xyz"),),
        )

        def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
            raise NotImplementedError

    register_tool(derep_registry, "presenttool", Present)
    register_tool(derep_registry, "absenttool", Absent)
    result = CliRunner().invoke(app, ["list-tools", "--check"])
    assert result.exit_code == 0, result.output
    assert "presenttool: ok (presenttool 9.9)" in result.output
    assert "absenttool: missing" in result.output
    assert "repgenr-no-such-binary-xyz" in result.output


def test_list_tools_shows_the_assembler_family() -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    line = next(ln for ln in result.output.splitlines() if ln.startswith("assemblers:"))
    assert "skesa" in line and "shovill" in line and "flye" in line


def test_list_tools_shows_declared_genome_limits() -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    lines = {ln.split(":")[0]: ln for ln in result.output.splitlines() if ": " in ln}
    assert "iqtree (up to 500 genomes)" in lines["treebuilders"]
    assert "mashtree (up to 10000 genomes)" in lines["treebuilders"]
    entries = lines["dereplicators"].split(": ", 1)[1].split(", ")
    assert "skder" in entries


def test_list_tools_names_the_glance_backends() -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    line = next(ln for ln in result.output.splitlines() if ln.startswith("glance"))
    assert "drep" in line and "sourmash" in line and "skder" not in line


def _isolated_families(monkeypatch, registry) -> None:
    """Every family lists only ``registry``'s adapters, so the status is known."""
    import importlib

    for family in (
        "aligners",
        "snptypers",
        "maskers",
        "treebuilders",
        "assemblers",
        "classifiers",
        "polishers",
        "annotators",
    ):
        module = importlib.import_module(f"repgenr.{family}.base")
        empty = Registry(f"repgenr.test_empty_{family}")
        empty._loaded = True
        monkeypatch.setattr(module, "registry", empty)
    import repgenr.dereplicators.base as derep_base

    monkeypatch.setattr(derep_base, "registry", registry)


class _OkAdapter:
    capabilities = ToolCapabilities(name="oktool")

    def preflight(self):
        return {"oktool": "1.0"}


class _MissingAdapter:
    capabilities = ToolCapabilities(name="missingtool")

    def preflight(self):
        raise MissingBinaryError("missingtool: not found on PATH")


def _bare(*adapters) -> Registry:
    reg = Registry("repgenr.test_strict")
    reg._loaded = True
    for adapter in adapters:
        reg.register(adapter.capabilities.name, adapter)
    return reg


@pytest.mark.parametrize(
    ("args", "adapters", "code"),
    [
        (["--check", "--strict"], (_OkAdapter,), 0),
        (["--check", "--strict"], (_OkAdapter, _MissingAdapter), 4),
        (["--check"], (_OkAdapter, _MissingAdapter), 0),
    ],
)
def test_list_tools_strict_exit_code(monkeypatch, args, adapters, code) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _isolated_families(monkeypatch, _bare(*adapters))
    result = CliRunner().invoke(app, ["list-tools", *args])
    assert result.exit_code == code, result.output
    # The full listing is printed before the exit.
    assert "polishers:" in result.output
    assert "glance" in result.output


def test_list_tools_strict_exits_5_on_a_broken_plugin(registry, monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _isolated_families(monkeypatch, registry)
    result = CliRunner().invoke(app, ["list-tools", "--check", "--strict"])
    assert result.exit_code == 5, result.output
    assert "badtool: broken" in result.output
    assert "polishers:" in result.output
    plain = CliRunner().invoke(app, ["list-tools", "--check"])
    assert plain.exit_code == 0


def test_list_tools_strict_exits_4_on_an_errored_adapter(monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    class _Errored:
        capabilities = ToolCapabilities(name="erroredtool")

        def preflight(self):
            raise RuntimeError("unexpected")

    _isolated_families(monkeypatch, _bare(_OkAdapter, _Errored))
    result = CliRunner().invoke(app, ["list-tools", "--check", "--strict"])
    assert result.exit_code == 4, result.output
    assert "erroredtool: error (RuntimeError: unexpected)" in result.output


def test_list_tools_strict_requires_check() -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    result = CliRunner().invoke(app, ["list-tools", "--strict"])
    assert result.exit_code == 2
    # Not an unknown-option error: the command names what --strict needs.
    assert "--strict needs --check" in _ANSI.sub("", result.output)


def test_list_tools_strict_prefers_5_when_broken_and_missing(registry, monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    registry.register("missingtool", _MissingAdapter)
    _isolated_families(monkeypatch, registry)
    result = CliRunner().invoke(app, ["list-tools", "--check", "--strict"])
    assert result.exit_code == 5, result.output
    assert "badtool: broken" in result.output
    assert "missingtool: missing" in result.output


def test_list_tools_verbose_shows_debug_lines(registry, monkeypatch) -> None:
    import logging

    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _isolated_families(monkeypatch, registry)
    real_names = registry.names

    def names():
        logging.getLogger("repgenr").debug("listing %s", registry.group)
        return real_names()

    monkeypatch.setattr(registry, "names", names)
    verbose = CliRunner().invoke(app, ["--verbose", "list-tools"])
    assert verbose.exit_code == 0
    assert "DEBUG listing repgenr.test_tools" in verbose.stderr
    plain = CliRunner().invoke(app, ["list-tools"])
    assert "DEBUG" not in plain.stderr


def test_list_tools_broken_plugin_warning_is_formatted(registry, monkeypatch) -> None:
    import re

    from typer.testing import CliRunner

    from repgenr.cli.main import app

    _isolated_families(monkeypatch, registry)
    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    assert re.search(
        r"^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d WARNING .*badtool", result.stderr, re.MULTILINE
    ), result.stderr
