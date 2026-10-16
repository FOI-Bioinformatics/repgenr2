"""`dereplicate --refine` and `derep-refine` on a working directory."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.cli import base as cli_base
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    REFINED_TSV,
    read_cluster_summary,
    read_clusters,
    read_refined,
)
from repgenr.core.errors import UserInputError
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.derep_refine import DerepRefineParams
from repgenr.stages.derep_refine import run as refine_run
from repgenr.stages.derep_stock import DerepStockParams
from repgenr.stages.derep_stock import run as stock_run
from repgenr.stages.dereplicate import DereplicateParams, run

# Accessions GCA_000XYZ: X is the coarse group, XY the fine group.
NAMES = [
    f"Francisellaceae_francisella_tularensis_GCA_000{a}.fasta"
    for a in (111, 112, 121, 211, 212, 221, 311)
]
A111, A112, A121, B211, B212, B221, C311 = NAMES


class _Tiered(Dereplicator):
    """Clusters by the coarse group at ANI <= 0.99 and by the fine group above."""

    capabilities = ToolCapabilities(name="tiered", supports_native_scaling=True)
    calls: list[tuple[int, float]] = []

    def preflight(self) -> dict[str, str]:
        return {"tiered": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = sorted(genomes, key=lambda p: p.name)
        type(self).calls.append((len(genomes), params.secondary_ani))
        width = 2 if params.secondary_ani > 0.99 else 1
        groups: dict[str, list[Path]] = {}
        for g in genomes:
            groups.setdefault(g.name.split("GCA_000")[1][:width], []).append(g)
        clusters = {grp[0].name: [g.name for g in grp[1:]] for grp in groups.values()}
        status: dict[str, str] = {}
        for rep, members in clusters.items():
            status[rep] = STATUS_REPRESENTATIVE
            status.update({m: STATUS_CONTAINED for m in members})
        return DerepResult(
            representatives=[grp[0] for grp in groups.values()],
            clusters=clusters,
            genome_status=status,
        )


@pytest.fixture
def tiered():
    registry._load()
    registry.register("tiered", _Tiered, replace=True)
    _Tiered.calls = []
    yield
    registry._classes.pop("tiered", None)


@pytest.fixture
def genomes(workdir: Path) -> list[Path]:
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    out = []
    for i, n in enumerate(NAMES):
        p = gdir / n
        p.write_text(f">s{i}\n{'ACGT' * (10 + i)}\n", encoding="utf-8")
        out.append(p)
    return out


def _ctx(workdir: Path) -> WorkdirContext:
    return WorkdirContext(workdir, create=True)


def test_dereplicate_refines_the_seed_cluster_in_one_run(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    result = run(
        ctx,
        DereplicateParams(tool="tiered", refine=["GCA_000112"], refine_ani=0.999, extra={"k": "v"}),
    )
    # The coarse pass over all seven, then the fine pass over cluster A's three.
    assert _Tiered.calls == [(7, 0.99), (3, 0.999)]
    assert read_clusters(ctx.derep_dir / CLUSTERS_TSV) == {
        A111: [A112],
        A121: [],
        B211: [B212, B221],
        C311: [],
    }
    assert sorted(p.name for p in ctx.representatives_dir.iterdir()) == [A111, A121, B211, C311]
    rows = read_refined(ctx.derep_dir / REFINED_TSV)
    assert [(r.genome, r.representative) for r in rows] == [
        (A111, A111),
        (A112, A111),
        (A121, A121),
    ]
    refined = {
        r.representative: r.refined
        for r in read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV)
    }
    assert refined == {A111: True, A121: True, B211: False, C311: False}
    rec = ctx.config.stages["dereplicate"]
    assert rec.params["refine_seeds"] == [A112]
    assert rec.params["refine_ani"] == 0.999
    assert rec.params["refine_genomes"] == 3 and rec.params["refine_representatives"] == 2
    assert rec.params["tool_args"] == {"k": "v"}
    assert len(result.representatives) == 4

    # A plain rerun leaves no refined.tsv of the earlier result behind.
    run(ctx, DereplicateParams(tool="tiered"))
    assert not (ctx.derep_dir / REFINED_TSV).exists()
    assert all(not r.refined for r in read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV))


def test_refine_refusals_in_the_stage(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    with pytest.raises(UserInputError, match="not finer"):
        run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.99))
    with pytest.raises(UserInputError, match="not under genomes/"):
        run(ctx, DereplicateParams(tool="tiered", refine=["GCA_999"], refine_ani=0.999))


def test_derep_refine_refines_the_result_on_disk(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    with pytest.raises(UserInputError, match="none is recorded"):
        refine_run(ctx, DerepRefineParams(seeds=[A111], ani=0.999))
    run(ctx, DereplicateParams(tool="tiered", extra={"k": "v"}))
    _Tiered.calls = []
    seen: dict[str, object] = {}

    class Spy(_Tiered):
        def dereplicate(self, genomes, out_dir, params, logger):
            seen["extra"] = params.extra
            return super().dereplicate(genomes, out_dir, params, logger)

    registry.register("tiered", Spy, replace=True)
    out = refine_run(ctx, DerepRefineParams(seeds=["GCA_000212", C311], ani=0.999))
    # Only the lineage was dereplicated, with the recorded tool arguments.
    assert _Tiered.calls == [(3, 0.999)]
    assert seen["extra"] == {"k": "v"}
    assert read_clusters(ctx.derep_dir / CLUSTERS_TSV) == {
        A111: [A112, A121],
        B211: [B212],
        B221: [],
        C311: [],
    }
    assert sorted(p.name for p in ctx.representatives_dir.iterdir()) == [A111, B211, B221, C311]
    assert sorted(p.name for p in out.representatives) == [A111, B211, B221, C311]
    rows = read_refined(ctx.derep_dir / REFINED_TSV)
    assert {r.genome for r in rows} == {B211, B212, B221}
    rec = ctx.config.stages["derep_refine"]
    assert rec.completed and rec.tool == "tiered"
    assert rec.params == {
        "refine_seeds": [B212],
        "refine_ani": 0.999,
        "refine_genomes": 3,
        "refine_representatives": 2,
        "refine_seeds_given": [C311, "GCA_000212"],
        "base_completed": ctx.config.stages["dereplicate"].completed,
    }
    # The dereplicate record is left as it was.
    assert "refine_ani" not in ctx.config.stages["dereplicate"].params
    # The manifest follows the refined result.
    status = {r.filename: (r.derep_status, r.representative) for r in ctx.manifest.all_genomes()}
    if status:
        assert status[B221][0] == STATUS_REPRESENTATIVE
    # Finer than the run but not finer than the record is refused.
    with pytest.raises(UserInputError, match="not finer"):
        refine_run(ctx, DerepRefineParams(seeds=[A111], ani=0.98))


def test_derep_refine_resumes_and_reruns_on_a_changed_seed(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered"))
    ctx.close()

    def invoke(seed: str) -> None:
        cli_base._run("derep_refine", workdir, lambda: DerepRefineParams(seeds=[seed], ani=0.999))

    _Tiered.calls = []
    invoke(B212)
    invoke(B212)  # the same request is skipped
    assert _Tiered.calls == [(3, 0.999)]
    invoke(A111)  # another seed reruns
    assert _Tiered.calls == [(3, 0.999), (3, 0.999)]


def test_stock_carries_the_refined_file(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999))
    stock_run(ctx, DerepStockParams(action="pack", name="fine"))
    assert (ctx.derep_dir / "stock" / "fine" / REFINED_TSV).is_file()
    run(ctx, DereplicateParams(tool="tiered"))
    stock_run(ctx, DerepStockParams(action="pack", name="coarse"))
    assert not (ctx.derep_dir / "stock" / "coarse" / REFINED_TSV).exists()

    stock_run(ctx, DerepStockParams(action="unpack", name="fine"))
    assert {r.genome for r in read_refined(ctx.derep_dir / REFINED_TSV)} == {A111, A112, A121}
    summary = {
        r.representative: r.refined
        for r in read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV)
    }
    assert summary[A111] and summary[A121] and not summary[B211]
    stock_run(ctx, DerepStockParams(action="unpack", name="coarse"))
    assert not (ctx.derep_dir / REFINED_TSV).exists()
    assert not any(r.refined for r in read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV))
