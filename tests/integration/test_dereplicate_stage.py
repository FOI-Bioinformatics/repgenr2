"""End-to-end dereplicate stage test using an in-process fake adapter."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    read_cluster_summary,
    read_clusters,
)
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.dereplicate import DereplicateParams, run


class _FakeDereplicator(Dereplicator):
    capabilities = ToolCapabilities(name="fake", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"fake": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = list(genomes)
        rep = genomes[0]
        members = [g.name for g in genomes[1:]]
        return DerepResult(
            representatives=[rep],
            clusters={rep.name: members},
            genome_status={
                rep.name: STATUS_REPRESENTATIVE,
                **{m: STATUS_CONTAINED for m in members},
            },
        )


class _NonScalingDereplicator(_FakeDereplicator):
    """Same behaviour but flagged as not scaling natively, to exercise chunking."""

    capabilities = ToolCapabilities(name="chunky", supports_native_scaling=False)


class _RecordingDereplicator(_FakeDereplicator):
    """Native-scaling adapter that records the (size, secondary_ani) of each call."""

    capabilities = ToolCapabilities(name="recording", supports_native_scaling=True)
    calls: list[tuple[int, float]] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = list(genomes)
        type(self).calls.append((len(genomes), params.secondary_ani))
        return super().dereplicate(genomes, out_dir, params, logger)


class _PrimaryBlindDereplicator(_FakeDereplicator):
    """Declares that it never reads primary_ani, as sourmash and skder do."""

    capabilities = ToolCapabilities(
        name="primaryblind", supports_native_scaling=True, ignored_params=frozenset({"primary_ani"})
    )


class _SettingsDereplicator(_NonScalingDereplicator):
    """Reports tool settings as the sourmash adapter does; the count of
    unknown pairs differs between calls of different sizes."""

    capabilities = ToolCapabilities(name="settings", supports_native_scaling=False)

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = list(genomes)
        result = super().dereplicate(genomes, out_dir, params, logger)
        result.settings = {
            "backend": "dense",
            "ani_estimate": "average_containment^(1/k)",
            "ksize": 31,
            "scaled": 1000,
            "unknown_pairs": len(genomes),
        }
        return result


@pytest.fixture
def fake_tool() -> None:
    registry._load()
    registry.register("fake", _FakeDereplicator, replace=True)
    registry.register("chunky", _NonScalingDereplicator, replace=True)
    registry.register("recording", _RecordingDereplicator, replace=True)
    registry.register("primaryblind", _PrimaryBlindDereplicator, replace=True)
    registry.register("settings", _SettingsDereplicator, replace=True)
    _RecordingDereplicator.calls = []
    yield
    registry._classes.pop("fake", None)
    registry._classes.pop("chunky", None)
    registry._classes.pop("recording", None)
    registry._classes.pop("primaryblind", None)
    registry._classes.pop("settings", None)


def test_stage_warns_when_tool_ignores_a_set_parameter(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    import logging

    ctx = WorkdirContext(workdir, create=True)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, DereplicateParams(tool="primaryblind", primary_ani=0.95))
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("does not use primary_ani" in m for m in warnings), warnings


def test_dereplicate_writes_contract(workdir: Path, genome_files, fake_tool) -> None:
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, DereplicateParams(tool="fake"))

    assert len(result.representatives) == 1

    rep_dir = ctx.representatives_dir
    assert rep_dir.is_dir()
    assert [p.name for p in rep_dir.iterdir()] == [genome_files[0].name]

    clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
    assert clusters[genome_files[0].name] == [genome_files[1].name, genome_files[2].name]

    assert (ctx.derep_dir / GENOME_STATUS_TSV).exists()
    (summary,) = read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV)
    assert summary.representative == genome_files[0].name and summary.n_members == 2
    assert (workdir / "repgenr.yaml").exists()
    assert ctx.config.stages["dereplicate"].tool == "fake"


def test_record_has_empty_tool_settings_for_a_tool_without_them(
    workdir: Path, genome_files, fake_tool
) -> None:
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, DereplicateParams(tool="fake", extra={"opt": "1"}))
    params = ctx.config.stages["dereplicate"].params
    assert params["tool_settings"] == {}
    # The existing keys are unchanged.
    assert params["tool_args"] == {"opt": "1"} and params["opt"] == "1"
    assert params["secondary_ani"] == DereplicateParams().secondary_ani


def test_record_has_the_tool_settings(workdir: Path, genome_files, fake_tool) -> None:
    from repgenr.core.config import Config

    ctx = WorkdirContext(workdir, create=True)
    run(ctx, DereplicateParams(tool="settings", keeper="quality"))
    params = Config.load(workdir).stages["dereplicate"].params
    assert params["tool_settings"] == {
        "backend": "dense",
        "ani_estimate": "average_containment^(1/k)",
        "ksize": 31,
        "scaled": 1000,
        "unknown_pairs": 3,
    }
    assert params["tool_args"] == {}


def test_chunked_record_names_the_back_end_once(workdir: Path, genome_files, fake_tool) -> None:
    gdir = workdir / "genomes"
    for i in range(4, 7):
        (gdir / f"Francisellaceae_f_t_GCA_00000{i}.fasta").write_text(f">s{i}\n{'ACGT' * 10}\n")
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, DereplicateParams(tool="settings", process_size=2))
    settings = ctx.config.stages["dereplicate"].params["tool_settings"]
    # Chunks of 2 and a merge of 3: the per-call count of unknown pairs is
    # left out, the settings every call shares are kept.
    assert settings == {
        "backend": "dense",
        "ani_estimate": "average_containment^(1/k)",
        "ksize": 31,
        "scaled": 1000,
    }


def test_chunking_composes_membership(workdir: Path, genome_files, fake_tool) -> None:
    # process_size=2 with 3 genomes -> the lone last chunk is merged, so a
    # single chunk runs. Use a larger set to force a real two-stage pass.
    gdir = workdir / "genomes"
    extra = []
    for i in range(4, 9):
        name = f"Francisellaceae_francisella_tularensis_GCA_00000{i}.fasta"
        (gdir / name).write_text(f">s{i}\n{'ACGT' * 10}\n")
        extra.append(name)

    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, DereplicateParams(tool="chunky", process_size=2))
    # every original genome must be accounted for as rep or contained
    all_names = {p.name for p in gdir.iterdir()}
    accounted = set(result.genome_status)
    assert all_names <= accounted


def test_native_scaling_single_pass_by_default(workdir: Path, genome_files, fake_tool) -> None:
    # No process_size -> native-scaling tool runs in a single pass (one adapter call).
    ctx = WorkdirContext(workdir, create=True)
    result = run(ctx, DereplicateParams(tool="recording"))
    assert len(result.representatives) == 1
    assert len(_RecordingDereplicator.calls) == 1  # one pass over all genomes


def test_native_scaling_can_be_chunked(workdir: Path, genome_files, fake_tool) -> None:
    # With process_size set and exceeded, even a native-scaling tool is chunked
    # (escape hatch for very large sets). 6 genomes / size 2 -> 3 stage-1 chunks
    # + 1 stage-2 pass = 4 adapter calls.
    gdir = workdir / "genomes"
    for i in range(4, 7):
        (gdir / f"Francisellaceae_f_t_GCA_00000{i}.fasta").write_text(f">s{i}\n{'ACGT' * 10}\n")
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, DereplicateParams(tool="recording", process_size=2))
    calls = _RecordingDereplicator.calls
    assert len(calls) == 4  # 3 stage-1 chunks + 1 stage-2
    assert sorted(n for n, _ in calls[:-1]) == [2, 2, 2]  # each stage-1 chunk has 2 genomes
    assert calls[-1][0] == 3  # stage-2 over the 3 chunk representatives


def test_keeper_quality_promotes_best_scoring_member(
    workdir: Path, genome_files, fake_tool
) -> None:
    # _FakeDereplicator clusters all three genomes under genome_files[0]. Seed the
    # manifest with quality so the tool's pick (rep) is the worst scorer, and
    # genome_files[2] the best -- the keeper should promote it to representative.
    ctx = WorkdirContext(workdir, create=True)
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    quality = {
        genome_files[0].name: (80.0, 5.0),  # score 55.0 -- worst
        genome_files[1].name: (95.0, 1.0),  # score 90.0
        genome_files[2].name: (99.0, 0.2),  # score 98.0 -- best
    }
    for filename, (completeness, contamination) in quality.items():
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(filename),
                filename=filename,
                completeness=completeness,
                contamination=contamination,
            )
        )

    result = run(ctx, DereplicateParams(tool="fake", keeper="quality"))

    assert len(result.representatives) == 1
    assert result.representatives[0].name == genome_files[2].name

    rep_dir = ctx.representatives_dir
    assert [p.name for p in rep_dir.iterdir()] == [genome_files[2].name]

    assert ctx.config.stages["dereplicate"].params["keeper"] == "quality"
    assert ctx.config.stages["dereplicate"].params["keeper_swaps"] == 1


def test_keeper_tool_keeps_adapter_pick(workdir: Path, genome_files, fake_tool) -> None:
    ctx = WorkdirContext(workdir, create=True)
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    quality = {
        genome_files[0].name: (80.0, 5.0),
        genome_files[1].name: (95.0, 1.0),
        genome_files[2].name: (99.0, 0.2),
    }
    for filename, (completeness, contamination) in quality.items():
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(filename),
                filename=filename,
                completeness=completeness,
                contamination=contamination,
            )
        )

    result = run(ctx, DereplicateParams(tool="fake", keeper="tool"))

    assert result.representatives[0].name == genome_files[0].name
    assert ctx.config.stages["dereplicate"].params["keeper"] == "tool"
    assert ctx.config.stages["dereplicate"].params["keeper_swaps"] == 0


def test_keeper_quality_without_manifest_quality_warns_and_records_fallback(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    # No completeness/contamination in the manifest (e.g. an API selection that
    # could not fetch quality): the stage must say so loudly and record that the
    # adapter's own picks were kept, so "quality, 0 swaps" is not ambiguous.
    import logging

    ctx = WorkdirContext(workdir, create=True)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, DereplicateParams(tool="fake", keeper="quality"))

    params = ctx.config.stages["dereplicate"].params
    assert params["keeper"] == "quality"
    assert params["keeper_effective"] == "tool"
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("quality" in m.lower() for m in warnings)


def test_keeper_effective_ignores_quality_of_genomes_not_dereplicated(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    # Quality only for the outgroup (kept under outgroup/, never dereplicated):
    # no representative can be re-picked, so the record must not claim the
    # quality rule was applied.
    import logging

    from repgenr.core.manifest import GenomeRecord

    ctx = WorkdirContext(workdir, create=True)
    ctx.logger.addHandler(caplog.handler)
    ctx.manifest.upsert(
        GenomeRecord(
            accession="GCF_999999999.1",
            filename="Fam_Gen_out_GCF_999999999.1.fasta",
            is_outgroup=True,
            completeness=99.0,
            contamination=0.5,
        )
    )
    with caplog.at_level(logging.WARNING):
        run(ctx, DereplicateParams(tool="fake", keeper="quality"))

    assert ctx.config.stages["dereplicate"].params["keeper_effective"] == "tool"
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("--keeper quality has no effect" in m for m in warnings)


def test_cluster_summary_command_regenerates_the_stage_summary(
    workdir: Path, genome_files, fake_tool
) -> None:
    """`cluster-summary` writes the same file as the stage, N50 columns included."""
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord
    from repgenr.stages import cluster_summary

    ctx = WorkdirContext(workdir, create=True)
    for i, f in enumerate(genome_files):
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(f.name),
                filename=f.name,
                completeness=99.0 - i,
                contamination=0.5,
            )
        )
    run(ctx, DereplicateParams(tool="fake", keeper="tool"))
    path = ctx.derep_dir / CLUSTER_SUMMARY_TSV
    written = path.read_text()
    (summary,) = read_cluster_summary(path)
    assert summary.rep_n50 == 40
    assert summary.best_score is not None
    path.unlink()
    cluster_summary.run(ctx, cluster_summary.ClusterSummaryParams())
    assert path.read_text() == written


def test_keeper_quality_with_manifest_quality_records_quality(
    workdir: Path, genome_files, fake_tool
) -> None:
    ctx = WorkdirContext(workdir, create=True)
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    for f in genome_files:
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(f.name),
                filename=f.name,
                completeness=99.0,
                contamination=0.5,
            )
        )
    run(ctx, DereplicateParams(tool="fake", keeper="quality"))
    assert ctx.config.stages["dereplicate"].params["keeper_effective"] == "quality"
    (summary,) = read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV)
    assert (summary.rep_completeness, summary.rep_contamination) == (99.0, 0.5)
    assert summary.member_max_completeness == 99.0


def test_stage1_uses_pre_thresholds(workdir: Path, genome_files, fake_tool) -> None:
    # 3 genomes, process_size=2 -> trailing singleton merges into one chunk, so
    # only stage 1 runs; bump to >=4 genomes to force a real two-stage pass.
    gdir = workdir / "genomes"
    for i in range(4, 7):
        (gdir / f"Francisellaceae_f_t_GCA_00000{i}.fasta").write_text(f">s{i}\n{'ACGT' * 10}\n")

    ctx = WorkdirContext(workdir, create=True)
    run(
        ctx,
        DereplicateParams(
            tool="recording",
            process_size=2,
            secondary_ani=0.99,
            pre_secondary_ani=0.95,
        ),
    )
    calls = _RecordingDereplicator.calls
    # stage-1 chunk calls use the looser pre threshold; the final stage-2 call
    # (on the union of stage-1 reps) uses the main threshold.
    stage1 = [s for s in calls[:-1]]
    assert stage1 and all(sec == 0.95 for _, sec in stage1)
    assert calls[-1][1] == 0.99


@pytest.mark.parametrize("unreachable", [False, True])
def test_missing_workdir_exits_3_without_creating_it(tmp_path: Path, unreachable: bool) -> None:
    """A nonexistent -wd is a workdir error (exit 3), not a traceback or a new directory."""
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    if unreachable:
        blocker = tmp_path / "a_file"
        blocker.write_text("")
        missing = blocker / "wd"  # mkdir under a regular file fails
    else:
        missing = tmp_path / "missing_wd"
    result = CliRunner().invoke(app, ["dereplicate", "-wd", str(missing)])
    assert result.exit_code == 3, result.output
    assert "Traceback" not in result.output
    assert not missing.exists()


class _ScratchCopyDereplicator(_FakeDereplicator):
    """Returns its representative as a copy in its own scratch dir, as skDER does."""

    capabilities = ToolCapabilities(name="scratchcopy", supports_native_scaling=True)

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        import shutil

        result = super().dereplicate(genomes, out_dir, params, logger)
        out_dir.mkdir(parents=True, exist_ok=True)
        copies = []
        for rep in result.representatives:
            dest = out_dir / rep.name
            shutil.copy2(rep, dest)
            copies.append(dest)
        result.representatives = copies
        return result


def test_representatives_link_the_genome_not_the_tool_copy(
    workdir: Path, genome_files, fake_tool
) -> None:
    registry.register("scratchcopy", _ScratchCopyDereplicator, replace=True)
    try:
        ctx = WorkdirContext(workdir, create=True)
        run(ctx, DereplicateParams(tool="scratchcopy"))
    finally:
        registry._classes.pop("scratchcopy", None)
    (rep,) = ctx.representatives_dir.iterdir()
    genome = ctx.genomes_dir / rep.name
    assert rep.stat().st_ino == genome.stat().st_ino


def test_refused_rerun_keeps_the_finished_record(workdir: Path, genome_files, fake_tool) -> None:
    """A selected genome gone from genomes/ refuses the rerun (exit 3) before the
    record of the finished dereplication is marked incomplete."""
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core.config import Config

    rows = ["accession\tfamily\tgenus\tspecies\tis_outgroup\tfilename\tcompleteness\tcontamination"]
    for i, g in enumerate(genome_files):
        rows.append(f"GCA_00000{i + 1}\tFrancisellaceae\tfrancisella\ttularensis\t0\t{g.name}\t\t")
    (workdir / "selection.tsv").write_text("\n".join(rows) + "\n")

    args = ["dereplicate", "-wd", str(workdir), "--tool", "fake"]
    first = CliRunner().invoke(app, args)
    assert first.exit_code == 0, first.output
    genome_files[1].unlink()
    second = CliRunner().invoke(app, args)
    assert second.exit_code == 3, second.output
    assert Config.load(workdir).stages["dereplicate"].completed


def test_drep_virus_mode_reaches_anim_through_the_stage(
    workdir: Path, genome_files, monkeypatch
) -> None:
    """The stage merges the adapter's default_params into the extras; dRep's
    virus default (ANImf) must survive that merge."""
    import repgenr.dereplicators.drep as drep_mod

    seen: list[list[str]] = []

    def record(caps, command, *, logger, **kwargs):
        seen.append([str(c) for c in command])
        raise RuntimeError("stop after recording the command")

    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {"dRep": "3"})
    monkeypatch.setattr(drep_mod, "run_tool", record)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(RuntimeError, match="stop after"):
        run(ctx, DereplicateParams(tool="drep", extra={"virus": True}))
    (cmd,) = seen
    assert cmd[cmd.index("--S_algorithm") + 1] == "ANImf"


@pytest.mark.parametrize("name", ["genome_status.tsv", "cluster_summary.tsv"])
def test_a_deleted_secondary_output_is_rebuilt_without_force(
    workdir: Path, genome_files, fake_tool, name: str
) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    args = ["dereplicate", "-wd", str(workdir), "--tool", "fake"]
    assert CliRunner().invoke(app, args).exit_code == 0
    (workdir / "derep" / name).unlink()
    again = CliRunner().invoke(app, args)
    assert again.exit_code == 0, again.output
    assert (workdir / "derep" / name).exists()


@pytest.mark.parametrize("tool", ["skder", "auto"])
def test_skder_ani_floor_refusal_keeps_the_finished_record(
    workdir: Path, genome_files, fake_tool, tool: str, monkeypatch
) -> None:
    """skDER's 80 percent floor is refused in the precheck (exit 2), before the
    harness marks the finished record incomplete; auto resolves to skDER here."""
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core import plugins
    from repgenr.core.config import Config

    # Hermetic auto: only skDER counts as runnable, whatever is on PATH and
    # whichever test adapters are registered.
    monkeypatch.setattr(plugins, "_tool_available", lambda caps: caps.name == "skder")

    wd = str(workdir)
    assert CliRunner().invoke(app, ["dereplicate", "-wd", wd, "--tool", "fake"]).exit_code == 0
    refused = CliRunner().invoke(app, ["dereplicate", "-wd", wd, "--tool", tool, "-sani", "0.78"])
    assert refused.exit_code == 2, refused.output
    record = Config.load(workdir).stages["dereplicate"]
    assert record.completed and record.tool == "fake"


class _QualityRecorder(_FakeDereplicator):
    capabilities = ToolCapabilities(name="qualityrecorder", supports_native_scaling=True)
    seen: list[dict] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        type(self).seen.append(dict(params.quality))
        return super().dereplicate(genomes, out_dir, params, logger)


@pytest.mark.parametrize("keeper", ["quality", "tool"])
def test_adapters_receive_the_manifest_quality(
    workdir: Path, genome_files, fake_tool, keeper: str
) -> None:
    """dRep's --genomeInfo and galah's input order read DerepParams.quality."""
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    registry.register("qualityrecorder", _QualityRecorder, replace=True)
    _QualityRecorder.seen = []
    ctx = WorkdirContext(workdir, create=True)
    quality = {g.name: (90.0 + i, 0.5 * i) for i, g in enumerate(genome_files)}
    for filename, (completeness, contamination) in quality.items():
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(filename),
                filename=filename,
                completeness=completeness,
                contamination=contamination,
            )
        )
    try:
        run(ctx, DereplicateParams(tool="qualityrecorder", keeper=keeper))
    finally:
        registry._classes.pop("qualityrecorder", None)
    assert _QualityRecorder.seen == [quality]


class _ChunkQualityRecorder(_QualityRecorder):
    capabilities = ToolCapabilities(name="chunkquality", supports_native_scaling=False)


@pytest.mark.parametrize("scored", [7, 6])
def test_every_chunk_and_the_merge_get_the_same_quality_source(
    workdir: Path, genome_files, fake_tool, scored: int
) -> None:
    """With quality for some genomes only, the run gives the adapter none, so
    a chunk whose genomes all happen to be scored is not treated differently
    from one that is not (dRep --genomeInfo in one chunk, CheckM in another)."""
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    gdir = workdir / "genomes"
    for i in range(4, 8):
        (gdir / f"Francisellaceae_francisella_tularensis_GCA_00000{i}.fasta").write_text(
            f">s{i}\n{'ACGT' * 10}\n"
        )
    names = sorted(p.name for p in gdir.iterdir())
    assert len(names) == 7
    registry.register("chunkquality", _ChunkQualityRecorder, replace=True)
    _QualityRecorder.seen = []
    ctx = WorkdirContext(workdir, create=True)
    for i, filename in enumerate(names[:scored]):
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(filename),
                filename=filename,
                completeness=95.0,
                contamination=float(i),
            )
        )
    try:
        run(ctx, DereplicateParams(tool="chunkquality", process_size=3))
    finally:
        registry._classes.pop("chunkquality", None)
    seen = _QualityRecorder.seen
    assert len(seen) == 3  # chunks of 3 and 4 genomes (a lone last genome joins), and the merge
    if scored == 7:
        assert all(len(q) == 7 for q in seen)
    else:
        # the first chunk is fully scored, yet no call gets quality
        assert all(q == {} for q in seen)


def _seed_gtdb(ctx: WorkdirContext, genome_files, reps: set[int]) -> None:
    from repgenr.core.contracts import accession_from_filename
    from repgenr.core.manifest import GenomeRecord

    for i, g in enumerate(genome_files):
        ctx.manifest.upsert(
            GenomeRecord(
                accession=accession_from_filename(g.name),
                filename=g.name,
                completeness=99.0 if i == 0 else 90.0,
                contamination=0.1 if i == 0 else 2.0,
                gtdb_representative=i in reps,
            )
        )


def test_keeper_gtdb_keeps_the_gtdb_representative(workdir: Path, genome_files, fake_tool) -> None:
    ctx = WorkdirContext(workdir, create=True)
    _seed_gtdb(ctx, genome_files, {2})
    result = run(ctx, DereplicateParams(tool="fake", keeper="gtdb"))
    assert [r.name for r in result.representatives] == [genome_files[2].name]
    params = ctx.config.stages["dereplicate"].params
    assert params["keeper"] == "gtdb"
    assert params["keeper_effective"] == "gtdb"
    assert params["keeper_swaps"] == 1
    (summary,) = read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV)
    assert summary.rep_is_gtdb_representative is True
    by_acc = {g.filename: g for g in ctx.manifest.all_genomes()}
    assert by_acc[genome_files[2].name].derep_status == STATUS_REPRESENTATIVE


def test_keeper_gtdb_without_flags_warns_and_uses_quality(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    import logging

    ctx = WorkdirContext(workdir, create=True)
    _seed_gtdb(ctx, genome_files, set())
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        result = run(ctx, DereplicateParams(tool="fake", keeper="gtdb"))
    assert [r.name for r in result.representatives] == [genome_files[0].name]
    assert ctx.config.stages["dereplicate"].params["keeper_effective"] == "quality"
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("GTDB species representative" in m for m in warnings), warnings


def test_cluster_summary_command_regenerates_the_gtdb_column(
    workdir: Path, genome_files, fake_tool
) -> None:
    from repgenr.stages import cluster_summary

    ctx = WorkdirContext(workdir, create=True)
    _seed_gtdb(ctx, genome_files, {2})
    run(ctx, DereplicateParams(tool="fake", keeper="gtdb"))
    path = ctx.derep_dir / CLUSTER_SUMMARY_TSV
    before = path.read_bytes()
    cluster_summary.run(ctx, cluster_summary.ClusterSummaryParams())
    assert path.read_bytes() == before


def test_keeper_gtdb_ignores_a_flagged_genome_outside_genomes_dir(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    """The outgroup is usually a GTDB representative; a flag on a genome that
    is not being dereplicated must not make the record claim the GTDB rule."""
    import logging

    from repgenr.core.manifest import GenomeRecord

    ctx = WorkdirContext(workdir, create=True)
    _seed_gtdb(ctx, genome_files, set())
    ctx.manifest.upsert(
        GenomeRecord(
            accession="GCA_000099.1",
            filename="Francisellaceae_francisella_other_GCA_000099.1.fasta",
            is_outgroup=True,
            completeness=99.0,
            contamination=0.0,
            gtdb_representative=True,
        )
    )
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, DereplicateParams(tool="fake", keeper="gtdb"))
    assert ctx.config.stages["dereplicate"].params["keeper_effective"] == "quality"
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("GTDB species representative" in m for m in warnings), warnings


class _FailFirstDereplicator(_FakeDereplicator):
    """Rejects the first genome on QC and clusters the rest under the second."""

    capabilities = ToolCapabilities(name="failfirst", supports_native_scaling=True)

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        from repgenr.dereplicators.base import STATUS_FAIL_QC

        genomes = list(genomes)
        result = super().dereplicate(genomes[1:], out_dir, params, logger)
        result.genome_status[genomes[0].name] = STATUS_FAIL_QC
        return result


def test_keeper_gtdb_ignores_a_flagged_genome_that_failed_qc(
    workdir: Path, genome_files, fake_tool, caplog
) -> None:
    import logging

    registry.register("failfirst", _FailFirstDereplicator, replace=True)
    try:
        ctx = WorkdirContext(workdir, create=True)
        _seed_gtdb(ctx, genome_files, {0})
        ctx.logger.addHandler(caplog.handler)
        with caplog.at_level(logging.WARNING):
            run(ctx, DereplicateParams(tool="failfirst", keeper="gtdb"))
    finally:
        registry._classes.pop("failfirst", None)
    assert ctx.config.stages["dereplicate"].params["keeper_effective"] == "quality"
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert any("GTDB species representative" in m for m in warnings), warnings
