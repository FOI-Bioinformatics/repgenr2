"""The `repgenr run` orchestrator chains the canonical stages in order."""

from __future__ import annotations

from typer.testing import CliRunner

from repgenr.cli import cmd_run
from repgenr.cli.main import app

_runner = CliRunner()


def _record(monkeypatch) -> list[str]:
    calls: list[str] = []

    def fake_run(stage, workdir, build, *, create=False):
        build()  # exercise the param builder (catches bad kwargs)
        calls.append(stage)

    monkeypatch.setattr(cmd_run, "_run", fake_run)
    # The wiring tests run without external tools; the preflight has its own test.
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    return calls


def test_run_bacterial_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-d",
            "rep",
            "-l",
            "genus",
            "-tg",
            "francisella",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["metadata", "genome", "dereplicate", "phylo", "tree2tax"]
    assert "Pipeline complete" in result.stdout


def test_run_viral_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "--viral",
            "--target",
            "mastadenovirus",
            "-tg",
            "Mastadenovirus",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["vmetadata", "vgenome", "dereplicate", "phylo", "tree2tax"]


def test_run_validates_tool(monkeypatch, tmp_path) -> None:
    _record(monkeypatch)
    result = _runner.invoke(app, ["run", "-wd", str(tmp_path), "--tool", "bogus"])
    assert result.exit_code != 0


def test_run_dry_run_previews_without_executing(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "--dry-run",
            "-l",
            "genus",
            "-tg",
            "francisella",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == []  # no stage executed
    assert "[dry-run]" in result.stdout
    assert "dereplicate" in result.stdout and "tree2tax" in result.stdout


def test_run_preflights_every_tool_before_the_first_stage(
    monkeypatch, tmp_path, register_tool
) -> None:
    # A missing tree builder must surface before metadata/genome download, not
    # after dereplication has finished.
    from repgenr.core.errors import MissingBinaryError
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.dereplicators.base import Dereplicator
    from repgenr.dereplicators.base import registry as derep_registry
    from repgenr.treebuilders.base import InputKind, TreeBuilder
    from repgenr.treebuilders.base import registry as tb_registry

    class OkDerep(Dereplicator):
        capabilities = ToolCapabilities(name="okderep")

        def preflight(self):
            return {"okderep": "1.0"}

        def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    class AbsentBuilder(TreeBuilder):
        capabilities = ToolCapabilities(name="absenttree")
        input_kind = InputKind.GENOMES

        def preflight(self):
            raise MissingBinaryError("absenttree: not found on PATH")

        def build(self, msa_or_genomes, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    register_tool(derep_registry, "okderep", OkDerep)
    register_tool(tb_registry, "absenttree", AbsentBuilder)
    calls: list[str] = []
    monkeypatch.setattr(cmd_run, "_run", lambda stage, *a, **k: calls.append(stage))
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-l",
            "genus",
            "-tg",
            "francisella",
            "--tool",
            "okderep",
            "--treebuilder",
            "absenttree",
        ],
    )
    assert result.exit_code != 0
    assert calls == []  # nothing downloaded, nothing dereplicated
    assert "absenttree" in result.output


def test_run_preflights_the_masker(monkeypatch, tmp_path, register_tool) -> None:
    from repgenr.core.errors import MissingBinaryError
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.dereplicators.base import Dereplicator
    from repgenr.dereplicators.base import registry as derep_registry
    from repgenr.maskers.base import Masker
    from repgenr.maskers.base import registry as masker_registry
    from repgenr.snptypers.base import SnpTyper
    from repgenr.snptypers.base import registry as snp_registry
    from repgenr.treebuilders.base import InputKind, TreeBuilder
    from repgenr.treebuilders.base import registry as tb_registry

    class OkDerep(Dereplicator):
        capabilities = ToolCapabilities(name="okderep")

        def preflight(self):
            return {"okderep": "1.0"}

        def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    class OkTree(TreeBuilder):
        capabilities = ToolCapabilities(name="oktree")
        input_kind = InputKind.MSA_FASTA

        def preflight(self):
            return {"oktree": "1.0"}

        def build(self, msa_or_genomes, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    class OkTyper(SnpTyper):
        capabilities = ToolCapabilities(name="oktyper")

        def preflight(self):
            return {"oktyper": "1.0"}

        def call(self, genomes, reference, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    class AbsentMasker(Masker):
        capabilities = ToolCapabilities(name="absentmask")

        def preflight(self):
            raise MissingBinaryError("absentmask: not found on PATH")

        def mask(self, full_alignment, out_dir, params, logger):  # noqa: ANN001
            raise AssertionError("must not run")

    register_tool(derep_registry, "okderep", OkDerep)
    register_tool(tb_registry, "oktree", OkTree)
    register_tool(snp_registry, "oktyper", OkTyper)
    register_tool(masker_registry, "absentmask", AbsentMasker)
    calls: list[str] = []
    monkeypatch.setattr(cmd_run, "_run", lambda stage, *a, **k: calls.append(stage))
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path),
            "-l",
            "genus",
            "-tg",
            "francisella",
            "--tool",
            "okderep",
            "--treebuilder",
            "oktree",
            "--msa-source",
            "snptype",
            "--snptyper",
            "oktyper",
            "--mask",
            "absentmask",
        ],
    )
    assert result.exit_code != 0
    assert calls == []
    assert "absentmask" in result.output


def test_run_local_chain_from_genomes_dir(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    result = _runner.invoke(
        app, ["run", "-wd", str(tmp_path / "wd"), "--genomes-dir", str(genomes)]
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["ingest", "dereplicate", "phylo", "tree2tax"]


def test_run_rejects_genomes_dir_with_viral(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            str(tmp_path / "wd"),
            "--genomes-dir",
            str(genomes),
            "--viral",
            "--target",
            "x",
        ],
    )
    assert result.exit_code != 0
    assert "--genomes-dir" in result.output
    assert calls == []


def test_run_with_snptype_inserts_the_stage(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        ["run", "-wd", str(tmp_path), "-l", "genus", "-tg", "francisella", "--with-snptype"],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["metadata", "genome", "dereplicate", "snptype", "phylo", "tree2tax"]


def test_run_with_snptype_and_snp_msa_source_keeps_the_natural_order(monkeypatch, tmp_path) -> None:
    """phylo's typing pass writes under tree/msa/, not snp/, so the standalone
    stage keeps its place before phylo with either MSA source."""
    calls = _record(monkeypatch)
    args = ["run", "-wd", str(tmp_path), "-l", "genus", "-tg", "francisella"]
    args += ["--with-snptype", "--msa-source", "snptype"]
    result = _runner.invoke(app, args)
    assert result.exit_code == 0, result.stdout
    assert calls == ["metadata", "genome", "dereplicate", "snptype", "phylo", "tree2tax"]
    dry = _runner.invoke(app, [*args, "--dry-run"])
    assert dry.exit_code == 0, dry.output
    listed = [ln.strip()[2:] for ln in dry.output.splitlines() if ln.strip().startswith("- ")]
    assert listed == ["metadata", "genome", "dereplicate", "snptype", "phylo", "tree2tax"]


def test_run_reads_chain(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        [
            *("run", "-wd", str(tmp_path / "wd"), "--reads", "-ts", "Francisella tularensis"),
            *("--platform", "illumina", "--max-runs", "5", "--assembler", "skesa"),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls == ["reads", "assemble", "dereplicate", "phylo", "tree2tax"]


def test_run_rejects_reads_with_viral(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app, ["run", "-wd", str(tmp_path / "wd"), "--reads", "--viral", "--target", "x"]
    )
    assert result.exit_code != 0 and "--reads" in result.output
    assert calls == []


def test_run_with_snptype_rejects_an_unknown_snptyper_up_front(monkeypatch, tmp_path) -> None:
    """--with-snptype uses --snptyper even with the aligner MSA source, so a bad
    name must exit 2 naming the flag (also under --dry-run), before any stage."""
    calls = _record(monkeypatch)
    base = ["run", "-wd", str(tmp_path / "wd"), "-l", "genus", "-tg", "francisella"]
    for extra in ([], ["--dry-run"]):
        result = _runner.invoke(app, [*base, "--with-snptype", "--snptyper", "bogus", *extra])
        assert result.exit_code == 2, result.output
        assert "--snptyper" in result.output
    assert calls == []


def test_run_with_snptype_passes_mask_to_the_snptype_stage(monkeypatch, tmp_path) -> None:
    """--with-snptype --mask gubbins masks the standalone snptype stage; with the
    aligner MSA source the phylo stage does not see the mask."""
    built: dict[str, object] = {}

    def fake_run(stage, workdir, build, *, create=False):
        built[stage] = build()

    monkeypatch.setattr(cmd_run, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    result = _runner.invoke(
        app,
        [
            *("run", "-wd", str(tmp_path / "wd"), "-l", "genus", "-tg", "francisella"),
            *("--with-snptype", "--mask", "gubbins"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert built["snptype"].mask == "gubbins"
    assert "mask" not in built["phylo"].extra


def test_run_mask_without_snptype_source_is_still_rejected(monkeypatch, tmp_path) -> None:
    calls = _record(monkeypatch)
    result = _runner.invoke(
        app,
        ["run", "-wd", str(tmp_path / "wd"), "-l", "genus", "-tg", "x", "--mask", "gubbins"],
    )
    assert result.exit_code == 2
    assert (
        "--mask applies only with --msa-source snptype, or on run with --with-snptype"
        in result.output
    )
    assert calls == []


def test_run_dry_run_warns_about_collapse_support_with_mashtree(
    monkeypatch, tmp_path, caplog
) -> None:
    _record(monkeypatch)
    import logging

    with caplog.at_level(logging.WARNING):
        result = _runner.invoke(
            app,
            [
                "run",
                "-wd",
                str(tmp_path),
                "--dry-run",
                "-l",
                "genus",
                "-tg",
                "francisella",
                "--treebuilder",
                "mashtree",
                "--collapse-support",
                "0.5",
            ],
        )
    assert result.exit_code == 0, result.stdout
    logged = any("carries no support values" in r.getMessage() for r in caplog.records)
    assert logged or "carries no support values" in result.output
