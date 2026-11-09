"""--target-reps: search --secondary-ani to land near a target rep count."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.dereplicate import DereplicateParams, run


class _AniDep(Dereplicator):
    """Rep count scales with the secondary-ani threshold: keep = round(ani * N)."""

    capabilities = ToolCapabilities(name="anidep", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"anidep": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        n = len(genomes)
        keep = max(1, min(n, round(params.secondary_ani * n)))
        reps = genomes[:keep]
        leftover = [g.name for g in genomes[keep:]]
        clusters: dict[str, list[str]] = {r.name: [] for r in reps}
        clusters[reps[0].name] = leftover  # park leftovers under the first rep
        status = {r.name: STATUS_REPRESENTATIVE for r in reps}
        for m in leftover:
            status[m] = STATUS_CONTAINED
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


@pytest.fixture
def anidep_workdir(workdir: Path):
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for i in range(10):
        (gdir / f"Fam_g_s_GCA_{i:06d}.1.fasta").write_text(">x\nACGT\n")
    registry._load()
    registry.register("anidep", _AniDep, replace=True)
    yield WorkdirContext(workdir, create=True)
    registry._classes.pop("anidep", None)


def test_target_reps_lands_near_target(anidep_workdir) -> None:
    res = run(anidep_workdir, DereplicateParams(tool="anidep", target_reps=9))
    assert abs(len(res.representatives) - 9) <= 1  # search converges close to 9
    assert len(res.genome_status) == 10  # all genomes accounted


def test_target_reps_search_is_written_to_derep_stats(anidep_workdir) -> None:
    from repgenr.core.stage_stats import read_stats
    from repgenr.stages.dereplicate import DEREP_STATS, DEREP_STATS_KIND

    res = run(anidep_workdir, DereplicateParams(tool="anidep", target_reps=9))
    stats = read_stats(anidep_workdir.derep_dir / DEREP_STATS, DEREP_STATS_KIND)
    assert stats is not None
    assert stats["tool"] == "anidep" and stats["input_genomes"] == 10
    search = stats["target_reps"]
    assert search["requested"] == 9
    assert search["achieved"] == len(res.representatives)
    assert 0.8 <= search["secondary_ani"] <= 1.0 and search["iterations"] >= 1
    # A run without the search rewrites the file without it.
    run(anidep_workdir, DereplicateParams(tool="anidep", secondary_ani=0.9))
    stats = read_stats(anidep_workdir.derep_dir / DEREP_STATS, DEREP_STATS_KIND)
    assert "target_reps" not in stats and stats["keeper_swaps"] == 0


def test_target_reps_off_uses_secondary_ani(anidep_workdir) -> None:
    # target_reps=0 -> no search; use --secondary-ani directly (0.9 -> keep 9 of 10)
    res = run(anidep_workdir, DereplicateParams(tool="anidep", secondary_ani=0.9))
    assert len(res.representatives) == 9
