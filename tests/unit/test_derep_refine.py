"""Finer dereplication inside a lineage: seeds, the fine pass and the splice."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.contracts import RefinedRow, read_refined, write_refined
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    DerepParams,
    DerepResult,
    check_result_complete,
)
from repgenr.stages import derep_refine as refine
from repgenr.stages.cluster_summary import summarise_clusters
from repgenr.stages.dereplicate import _reduce_by_taxonomy

_LOG = logging.getLogger("test-refine")

# Accessions GCA_000XYZ: X is the coarse group, XY the fine group.
NAMES = [
    f"Francisellaceae_francisella_tularensis_GCA_000{a}.fasta"
    for a in (111, 112, 121, 211, 212, 221, 311)
]
A111, A112, A121, B211, B212, B221, C311 = NAMES


def _paths(tmp_path: Path) -> dict[str, Path]:
    paths = {}
    for n in NAMES:
        p = tmp_path / n
        p.write_text(">s\nACGT\n", encoding="utf-8")
        paths[n] = p
    return paths


def _coarse(paths: dict[str, Path]) -> DerepResult:
    """Two coarse clusters (A and B) and a singleton C."""
    clusters = {A111: [A112, A121], B211: [B212, B221], C311: []}
    status = {}
    for rep, members in clusters.items():
        status[rep] = STATUS_REPRESENTATIVE
        status.update({m: STATUS_CONTAINED for m in members})
    return DerepResult(
        representatives=[paths[r] for r in clusters], clusters=clusters, genome_status=status
    )


class _Tiered:
    """Clusters by the coarse group at ANI <= 0.99 and by the fine group above."""

    calls: list[tuple[list[str], float]] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        genomes = sorted(genomes, key=lambda p: p.name)
        type(self).calls.append(([g.name for g in genomes], params.secondary_ani))
        width = 2 if params.secondary_ani > 0.99 else 1
        groups: dict[str, list[Path]] = {}
        for g in genomes:
            key = g.name.split("GCA_000")[1][:width]
            groups.setdefault(key, []).append(g)
        clusters = {grp[0].name: [g.name for g in grp[1:]] for grp in groups.values()}
        status = {}
        for rep, members in clusters.items():
            status[rep] = STATUS_REPRESENTATIVE
            status.update({m: STATUS_CONTAINED for m in members})
        return DerepResult(
            representatives=[grp[0] for grp in groups.values()],
            clusters=clusters,
            genome_status=status,
        )


def test_refine_ani_must_be_finer_than_the_run() -> None:
    refine.check_refine_ani(0.999, 0.99)
    with pytest.raises(UserInputError, match="not finer"):
        refine.check_refine_ani(0.99, 0.99)
    with pytest.raises(UserInputError, match="not finer"):
        refine.check_refine_ani(0.95, 0.99)
    with pytest.raises(UserInputError, match=r"\(0, 1\]"):
        refine.check_refine_ani(1.5, 0.99)


def test_seeds_resolve_by_filename_stem_and_accession(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    genomes = list(paths.values())
    assert refine.resolve_seeds([A111], genomes) == [A111]
    assert refine.resolve_seeds([A111.removesuffix(".fasta")], genomes) == [A111]
    assert refine.resolve_seeds(["GCA_000211", " GCA_000211 "], genomes) == [B211]
    # The manifest maps accessions the filename does not carry.
    assert refine.resolve_seeds(["X9"], genomes, accessions={"X9": C311}) == [C311]
    with pytest.raises(UserInputError, match="not under genomes/.*GCA_999999"):
        refine.resolve_seeds(["GCA_999999"], genomes)
    with pytest.raises(UserInputError, match="outgroup"):
        refine.resolve_seeds(["Out_g_s_GCA_7.fasta"], genomes, outgroup={"Out_g_s_GCA_7.fasta"})


def test_lineage_is_every_cluster_holding_a_seed(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    coarse = _coarse(paths)
    # A member seed brings its whole cluster; two seeds in one cluster count once.
    lineage = refine.lineage_of(coarse, [A121, A112], _LOG)
    assert lineage.coarse_representatives == [A111]
    assert lineage.genomes == [A111, A112, A121]
    assert lineage.seed_of == {A111: A121, A112: A121, A121: A121}
    assert lineage.coarse_rep_of[A121] == A111
    # Two clusters, and a singleton seed is skipped.
    lineage = refine.lineage_of(coarse, [B212, C311, A111], _LOG)
    assert lineage.coarse_representatives == [B211, A111]
    assert lineage.skipped == [C311]
    coarse.genome_status[A121] = STATUS_FAIL_QC
    with pytest.raises(UserInputError, match="quality check"):
        refine.lineage_of(coarse, [A121], _LOG)
    with pytest.raises(WorkdirError, match="no cluster"):
        refine.lineage_of(coarse, ["nowhere.fasta"], _LOG)


def test_fine_pass_and_splice_replace_only_the_seed_clusters(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    coarse = _coarse(paths)
    lineage = refine.lineage_of(coarse, [A112], _LOG)
    _Tiered.calls = []
    params = DerepParams(
        secondary_ani=0.99,
        quality={A111: (99.0, 0.1), B211: (90.0, 1.0)},
        sketches={paths[A111]: tmp_path / "a.sig", paths[B211]: tmp_path / "b.sig"},
    )
    sub = refine.fine_pass(_Tiered(), lineage, paths, tmp_path / "scratch", params, 0.999, _LOG)
    # Only the lineage genomes, at the fine ANI, with the quality and sketches of those.
    ((called, ani),) = _Tiered.calls
    assert called == [A111, A112, A121] and ani == 0.999
    assert sub.clusters == {A111: [A112], A121: []}

    out = refine.splice(coarse, lineage, sub, paths, _LOG)
    assert out.clusters == {B211: [B212, B221], C311: [], A111: [A112], A121: []}
    assert out.genome_status[A121] == STATUS_REPRESENTATIVE
    assert [p.name for p in out.representatives] == sorted([A111, A121, B211, C311])
    check_result_complete(out, NAMES)

    rows = refine.refined_rows(out, lineage)
    assert rows == [
        RefinedRow(A111, A112, A111, A111, STATUS_REPRESENTATIVE),
        RefinedRow(A112, A112, A111, A111, STATUS_CONTAINED),
        RefinedRow(A121, A112, A111, A121, STATUS_REPRESENTATIVE),
    ]
    assert refine.record_params(lineage, 0.999, rows) == {
        "refine_seeds": [A112],
        "refine_ani": 0.999,
        "refine_genomes": 3,
        "refine_representatives": 2,
    }


def test_fine_pass_quality_subset_reaches_the_adapter(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    lineage = refine.lineage_of(_coarse(paths), [A111], _LOG)
    seen = {}

    class Spy(_Tiered):
        def dereplicate(self, genomes, out_dir, params, logger):
            seen["quality"] = params.quality
            seen["sketches"] = params.sketches
            seen["primary"] = params.primary_ani
            return super().dereplicate(genomes, out_dir, params, logger)

    params = DerepParams(
        primary_ani=0.8,
        quality={A111: (99.0, 0.1), B211: (90.0, 1.0)},
        sketches={paths[A111]: tmp_path / "a.sig", paths[B211]: tmp_path / "b.sig"},
    )
    refine.fine_pass(Spy(), lineage, paths, tmp_path / "scratch", params, 0.999, _LOG)
    assert seen == {
        "quality": {A111: (99.0, 0.1)},
        "sketches": {paths[A111]: tmp_path / "a.sig"},
        "primary": 0.8,
    }


def test_a_genome_the_fine_pass_rejects_is_marked_fail_qc(tmp_path: Path, caplog) -> None:
    paths = _paths(tmp_path)
    coarse = _coarse(paths)
    lineage = refine.lineage_of(coarse, [A111], _LOG)
    sub = DerepResult(
        representatives=[paths[A111]],
        clusters={A111: [A121]},
        genome_status={A111: STATUS_REPRESENTATIVE, A121: STATUS_CONTAINED, A112: STATUS_FAIL_QC},
    )
    with caplog.at_level(logging.WARNING):
        out = refine.splice(coarse, lineage, sub, paths, _LOG)
    assert out.genome_status[A112] == STATUS_FAIL_QC
    assert "failed QC in the fine pass" in caplog.text
    check_result_complete(out, NAMES)
    (row,) = [r for r in refine.refined_rows(out, lineage) if r.genome == A112]
    assert row.representative == "" and row.status == STATUS_FAIL_QC


def test_reduce_keeps_the_refined_clusters_apart(tmp_path: Path) -> None:
    paths = _paths(tmp_path)
    clusters = {A111: [A112], A121: [], B211: [B212, B221], C311: []}
    status = {}
    for rep, members in clusters.items():
        status[rep] = STATUS_REPRESENTATIVE
        status.update({m: STATUS_CONTAINED for m in members})
    result = DerepResult(
        representatives=[paths[r] for r in clusters], clusters=clusters, genome_status=status
    )
    taxon = {n: "tularensis" for n in NAMES}
    folded = _reduce_by_taxonomy(result, "species", {}, _LOG, taxon_of=taxon)
    assert len(folded.clusters) == 1
    kept = _reduce_by_taxonomy(result, "species", {}, _LOG, taxon_of=taxon, keep_apart={A111, A121})
    # The refined clusters stay; the others fold into the largest of them.
    assert set(kept.clusters) == {A111, A121, B211}
    assert sorted(kept.clusters[B211]) == [B212, B221, C311]


def test_refined_rows_round_trip_and_summary_column(tmp_path: Path) -> None:
    rows = [
        RefinedRow(A112, A112, A111, A111, STATUS_CONTAINED),
        RefinedRow(A111, A112, A111, A111, STATUS_REPRESENTATIVE),
        RefinedRow(A121, A112, A111, A121, STATUS_REPRESENTATIVE),
    ]
    path = tmp_path / "refined.tsv"
    write_refined(path, rows)
    assert read_refined(path) == sorted(rows, key=lambda r: (r.representative, r.genome))
    assert read_refined(tmp_path / "none.tsv") == []
    summary = summarise_clusters({A111: [A112], A121: [], B211: [B212]}, {}, refined={A111, A121})
    assert {r.representative: r.refined for r in summary} == {
        A111: True,
        A121: True,
        B211: False,
    }
