"""`dereplicate --refine` and `derep-refine` on a working directory."""

from __future__ import annotations

import logging
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
    return WorkdirContext(workdir, create=True, logger=logging.getLogger("test-refine"))


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
        "refine_lineages": [{"seeds": [B212], "ani": 0.999, "genomes": 3, "representatives": 2}],
        "refine_request": [[[C311, "GCA_000212"], 0.999]],
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


def _invoke(workdir: Path, seeds: list[str], ani: float = 0.999) -> None:
    cli_base._run("derep_refine", workdir, lambda: DerepRefineParams(seeds=seeds, ani=ani))


def test_a_later_refinement_keeps_the_earlier_one(workdir: Path, genomes, tiered) -> None:
    """Refining B after A keeps A's fine clusters and its rows; a plain
    dereplicate afterwards drops the derep_refine record with the result."""
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999))
    ctx.close()
    _Tiered.calls = []
    _invoke(workdir, [B212])
    assert _Tiered.calls == [(3, 0.999)]
    ctx = WorkdirContext(workdir)
    assert read_clusters(ctx.derep_dir / CLUSTERS_TSV) == {
        A111: [A112],
        A121: [],
        B211: [B212],
        B221: [],
        C311: [],
    }
    rows = {r.genome: r for r in read_refined(ctx.derep_dir / REFINED_TSV)}
    assert set(rows) == {A111, A112, A121, B211, B212, B221}
    assert rows[A121].seed == A111 and rows[B221].seed == B212
    summary = {
        r.representative: r.refined
        for r in read_cluster_summary(ctx.derep_dir / CLUSTER_SUMMARY_TSV)
    }
    assert summary == {A111: True, A121: True, B211: True, B221: True, C311: False}
    _Tiered.calls = []
    _invoke(workdir, [B212])  # a repeat is skipped by the harness
    assert _Tiered.calls == []
    ctx = WorkdirContext(workdir)
    run(ctx, DereplicateParams(tool="tiered"))
    assert "derep_refine" not in ctx.config.stages
    ctx.close()
    _Tiered.calls = []
    _invoke(workdir, [B212])  # the new dereplication is refined again
    assert _Tiered.calls == [(3, 0.999)]


def test_refining_a_seed_again_starts_from_its_coarse_cluster(
    workdir: Path, genomes, tiered
) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999))
    ctx.close()
    _Tiered.calls = []
    _invoke(workdir, [A111], ani=0.9995)
    # All three genomes of the original coarse cluster, not the two of A111's fine cluster.
    assert _Tiered.calls == [(3, 0.9995)]
    ctx = WorkdirContext(workdir)
    rows = {r.genome: r for r in read_refined(ctx.derep_dir / REFINED_TSV)}
    assert set(rows) == {A111, A112, A121}
    assert {r.coarse_representative for r in rows.values()} == {A111}
    assert ctx.config.stages["derep_refine"].params["refine_ani"] == 0.9995


def test_refine_is_refused_before_the_coarse_pass(workdir: Path, genomes, tiered) -> None:
    from repgenr.stages.dereplicate import precheck

    ctx = _ctx(workdir)
    for params, text in [
        (DereplicateParams(tool="tiered", refine=["GCA_999"], refine_ani=0.999), "not under"),
        (DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.9), "not finer"),
        (
            DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999, target_reps=2),
            "--target-reps",
        ),
    ]:
        # The harness runs the precheck before the record is touched.
        with pytest.raises(UserInputError, match=text):
            precheck(ctx, params)
        with pytest.raises(UserInputError, match=text):
            run(ctx, params)
    assert _Tiered.calls == []


def test_refine_with_reduce_keeps_the_lineage_apart(workdir: Path, genomes, tiered) -> None:
    """Every genome is F. tularensis by filename: --reduce species folds B and
    C together and leaves the two refined clusters."""
    from repgenr.core.manifest import GenomeRecord

    ctx = _ctx(workdir)
    ctx.manifest.upsert_many(
        [
            GenomeRecord(
                accession=f"GCA_000{a}",
                filename=n,
                genus="francisella",
                species="tularensis",
            )
            for a, n in zip((111, 112, 121, 211, 212, 221, 311), NAMES, strict=True)
        ]
    )
    run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999, reduce="species"))
    clusters = read_clusters(ctx.derep_dir / CLUSTERS_TSV)
    assert {A111, A121} < set(clusters) and len(clusters) == 3
    other = next(rep for rep in clusters if rep not in (A111, A121))
    assert sorted([other, *clusters[other]]) == [B211, B212, B221, C311]


def test_refine_after_a_chunked_coarse_pass(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", process_size=4, refine=[B211], refine_ani=0.999))
    assert _Tiered.calls[-1] == (3, 0.999)
    assert read_clusters(ctx.derep_dir / CLUSTERS_TSV) == {
        A111: [A112, A121],
        B211: [B212],
        B221: [],
        C311: [],
    }


def test_derep_refine_with_a_gtdb_keeper_record_warns_nothing_false(
    workdir: Path, genomes, tiered, caplog
) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", keeper="gtdb"))
    caplog.clear()
    with caplog.at_level(logging.INFO):
        refine_run(ctx, DerepRefineParams(seeds=[A111], ani=0.999))
    assert "metadata --force" not in caplog.text
    assert "no GTDB species representative among the lineage genomes" in caplog.text


def test_two_lineages_each_at_their_own_ani(workdir: Path, genomes, tiered, tmp_path) -> None:
    """A refine file with A at 0.999 and B at 0.9995: each lineage is
    dereplicated at its ANI, and the rows carry it."""
    refine_file = tmp_path / "refine.tsv"
    refine_file.write_text(f"seed\tani\n{A112}\t0.999\nGCA_000212\t0.9995\n", encoding="utf-8")
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", refine_file=str(refine_file)))
    assert _Tiered.calls == [(7, 0.99), (3, 0.999), (3, 0.9995)]
    assert read_clusters(ctx.derep_dir / CLUSTERS_TSV) == {
        A111: [A112],
        A121: [],
        B211: [B212],
        B221: [],
        C311: [],
    }
    rows = {r.genome: r for r in read_refined(ctx.derep_dir / REFINED_TSV)}
    assert rows[A121].ani == "0.999" and rows[B221].ani == "0.9995"
    rec = ctx.config.stages["dereplicate"].params
    assert rec["refine_ani"] is None and rec["refine_file"] == str(refine_file)
    assert [lin["ani"] for lin in rec["refine_lineages"]] == [0.999, 0.9995]

    # derep-refine takes a file too, and may add the CLI pair on top.
    run(ctx, DereplicateParams(tool="tiered"))
    ctx.close()
    _Tiered.calls = []
    cli_base._run(
        "derep_refine",
        workdir,
        lambda: DerepRefineParams(seeds=[C311], ani=0.999, refine_file=str(refine_file)),
    )
    assert sorted(_Tiered.calls) == [(3, 0.999), (3, 0.9995)]  # C311 is a singleton
    ctx = WorkdirContext(workdir)
    assert ctx.config.stages["derep_refine"].params["refine_request"] == [
        [[A112, C311], 0.999],
        [["GCA_000212"], 0.9995],
    ]
    # doctor rebuilds the stage inputs from the record: nothing is stale.
    from repgenr.core.doctor import _changed_by_stage

    assert _changed_by_stage(workdir, ctx.config) == {}
    # Editing the file changes the request: the stage reruns.
    refine_file.write_text(f"seed\tani\n{A112}\t0.9995\n", encoding="utf-8")
    ctx.close()
    _Tiered.calls = []
    cli_base._run("derep_refine", workdir, lambda: DerepRefineParams(refine_file=str(refine_file)))
    assert _Tiered.calls == [(3, 0.9995)]


def test_a_second_refinement_counts_its_own_representatives(workdir: Path, genomes, tiered) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999))
    ctx.close()
    _invoke(workdir, [B212])
    ctx = WorkdirContext(workdir)
    rec = ctx.config.stages["derep_refine"].params
    assert rec["refine_genomes"] == 3 and rec["refine_representatives"] == 2
    assert {r.genome for r in read_refined(ctx.derep_dir / REFINED_TSV)} == {
        A111,
        A112,
        A121,
        B211,
        B212,
        B221,
    }


def test_two_seeds_of_one_refined_coarse_cluster_at_two_anis_are_refused(
    workdir: Path, genomes, tiered, tmp_path
) -> None:
    ctx = _ctx(workdir)
    run(ctx, DereplicateParams(tool="tiered", refine=[A111], refine_ani=0.999))
    refine_file = tmp_path / "refine.tsv"
    # A111 and A121 now sit in different fine clusters of one coarse cluster.
    refine_file.write_text(f"seed\tani\n{A111}\t0.9995\n{A121}\t0.9999\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="share a cluster"):
        refine_run(ctx, DerepRefineParams(refine_file=str(refine_file)))
