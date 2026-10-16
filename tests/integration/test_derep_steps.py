"""Discrete dereplicate-chunk / dereplicate-merge steps (scatter-gather substrate)."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import (
    CLUSTER_SUMMARY_TSV,
    CLUSTERS_TSV,
    GENOME_STATUS_TSV,
    SelectionRow,
    read_clusters,
    read_genome_status,
    read_selection,
    write_selection,
)
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.derep_steps import (
    ChunkParams,
    MergeParams,
    dereplicate_chunk,
    dereplicate_merge,
)
from repgenr.stages.dereplicate import DereplicateParams
from repgenr.stages.dereplicate import run as dereplicate_run

_LOG = __import__("logging").getLogger("test")


class _Halver(Dereplicator):
    """Collapses each adjacent pair to one rep, so every pass ~halves the set."""

    capabilities = ToolCapabilities(name="halver", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"halver": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        reps = genomes[::2]
        clusters: dict[str, list[str]] = {}
        status: dict[str, str] = {}
        for i, rep in enumerate(reps):
            members = [genomes[2 * i + 1].name] if 2 * i + 1 < len(genomes) else []
            clusters[rep.name] = members
            status[rep.name] = STATUS_REPRESENTATIVE
            for m in members:
                status[m] = STATUS_CONTAINED
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


class _QcHalver(_Halver):
    """Halver that rejects the genomes in ``fail_names`` on QC before clustering."""

    capabilities = ToolCapabilities(name="qchalver", supports_native_scaling=True)
    fail_names: frozenset[str] = frozenset()

    def preflight(self) -> dict[str, str]:
        return {"qchalver": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        failed = [g for g in genomes if g.name in type(self).fail_names]
        kept = [g for g in genomes if g.name not in type(self).fail_names]
        result = super().dereplicate(kept, out_dir, params, logger)
        for g in failed:
            result.genome_status[g.name] = STATUS_FAIL_QC
        return result


@pytest.fixture
def reg():
    registry._load()
    registry.register("halver", _Halver, replace=True)
    registry.register("qchalver", _QcHalver, replace=True)
    yield
    registry._classes.pop("halver", None)
    registry._classes.pop("qchalver", None)
    _QcHalver.fail_names = frozenset()


def _make_genomes(gdir: Path, n: int) -> list[Path]:
    gdir.mkdir(parents=True, exist_ok=True)
    out = []
    for i in range(n):
        p = gdir / f"Fam_g_s_GCA_{i:06d}.1.fasta"
        p.write_text(">x\nACGT\n")
        out.append(p)
    return out


def test_chunk_writes_a_valid_contract(tmp_path: Path, reg) -> None:
    genomes = _make_genomes(tmp_path / "genomes", 4)
    out = tmp_path / "chunk0"
    res = dereplicate_chunk(ChunkParams(tool="halver", genomes=genomes, out_dir=out), _LOG)
    # halver of 4 -> reps [g0, g2]
    assert {r.name for r in res.representatives} == {genomes[0].name, genomes[2].name}
    # contract files + representative FASTAs are present on disk
    assert (out / CLUSTERS_TSV).exists()
    assert (out / GENOME_STATUS_TSV).exists()
    assert (out / CLUSTER_SUMMARY_TSV).exists()
    rep_files = {p.name for p in (out / "representatives").iterdir()}
    assert rep_files == {genomes[0].name, genomes[2].name}
    clusters = read_clusters(out / CLUSTERS_TSV)
    assert clusters[genomes[0].name] == [genomes[1].name]
    # tool intermediates are trimmed; the staged representatives survive
    assert not (out / "scratch").exists()
    assert (out / "representatives" / genomes[0].name).read_text() == ">x\nACGT\n"


def test_merge_composes_membership_over_chunks(tmp_path: Path, reg) -> None:
    genomes = _make_genomes(tmp_path / "genomes", 8)
    # two chunks of four
    c0 = dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    c1 = dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )
    assert len(c0.representatives) == 2 and len(c1.representatives) == 2

    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
        ),
        _LOG,
    )
    # union [g0,g2,g4,g6] halved -> [g0, g4]; composition pulls in every original.
    assert {r.name for r in final.representatives} == {genomes[0].name, genomes[4].name}
    clusters = read_clusters(tmp_path / "merged" / CLUSTERS_TSV)
    assert sorted(clusters[genomes[0].name]) == sorted(
        [genomes[1].name, genomes[2].name, genomes[3].name]
    )
    assert sorted(clusters[genomes[4].name]) == sorted(
        [genomes[5].name, genomes[6].name, genomes[7].name]
    )
    # every original genome accounted for exactly once (8 = 2 reps + 6 members)
    all_members = [m for ms in clusters.values() for m in ms]
    assert len(all_members) + len(clusters) == 8


def test_discrete_matches_in_process_chunked(workdir: Path, reg) -> None:
    """Scatter-gather steps reproduce the in-process two-stage membership."""
    genomes = _make_genomes(workdir / "genomes", 8)

    # Reference: the shared-workdir stage with the same chunk size.
    ctx = WorkdirContext(workdir, create=True)
    ref = dereplicate_run(ctx, DereplicateParams(tool="halver", process_size=4))

    # Discrete: two chunks then a merge, with the same final thresholds.
    dereplicate_chunk(ChunkParams(tool="halver", genomes=genomes[:4], out_dir=workdir / "c0"), _LOG)
    dereplicate_chunk(ChunkParams(tool="halver", genomes=genomes[4:], out_dir=workdir / "c1"), _LOG)
    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[workdir / "c0", workdir / "c1"],
            out_dir=workdir / "merged",
        ),
        _LOG,
    )

    def membership(result: DerepResult) -> dict[str, list[str]]:
        return {rep: sorted(members) for rep, members in result.clusters.items()}

    assert {r.name for r in final.representatives} == {r.name for r in ref.representatives}
    assert membership(final) == membership(ref)


def test_fail_qc_genome_survives_chunk_merge(tmp_path: Path, reg) -> None:
    """A chunk's QC-rejected genome keeps its status through the merge step.

    The genome is in no cluster, so a merge that reads only clusters.tsv loses it
    from both the composed status and the checked name list.
    """
    genomes = _make_genomes(tmp_path / "genomes", 8)
    bad = genomes[3].name
    _QcHalver.fail_names = frozenset({bad})

    c0 = dereplicate_chunk(
        ChunkParams(tool="qchalver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    assert c0.genome_status[bad] == STATUS_FAIL_QC
    assert read_genome_status(tmp_path / "c0" / GENOME_STATUS_TSV)[bad] == STATUS_FAIL_QC
    dereplicate_chunk(
        ChunkParams(tool="qchalver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )

    final = dereplicate_merge(
        MergeParams(
            tool="qchalver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
        ),
        _LOG,
    )

    assert final.genome_status[bad] == STATUS_FAIL_QC
    assert len(final.genome_status) == 8
    on_disk = read_genome_status(tmp_path / "merged" / GENOME_STATUS_TSV)
    assert on_disk[bad] == STATUS_FAIL_QC
    assert len(on_disk) == 8
    # the rejected genome is in no cluster
    clusters = read_clusters(tmp_path / "merged" / CLUSTERS_TSV)
    assert bad not in {m for members in clusters.values() for m in members}
    assert bad not in clusters


def test_chunk_rejects_missing_genome(tmp_path: Path, reg) -> None:
    from repgenr.core.errors import WorkdirError

    bogus = [tmp_path / "nope.fasta"]
    with pytest.raises(WorkdirError):
        dereplicate_chunk(ChunkParams(tool="halver", genomes=bogus, out_dir=tmp_path / "x"), _LOG)


def _write_two_genome_selection(path: Path, genomes: list[Path]) -> None:
    """A minimal selection.tsv giving quality only to ``genomes[0]`` and
    ``genomes[2]`` -- exactly the two chunk-0 halver representatives -- with
    genomes[2] scoring better, so a quality keeper should prefer it."""
    write_selection(
        path,
        [
            SelectionRow(
                accession="GCF_000000.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=genomes[0].name,
                completeness=80.0,
                contamination=5.0,
            ),
            SelectionRow(
                accession="GCF_000002.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=genomes[2].name,
                completeness=99.0,
                contamination=0.5,
            ),
        ],
    )


def test_merge_applies_quality_keeper_from_selection_tsv(tmp_path: Path, reg) -> None:
    """A selection.tsv with quality columns lets the merge step's keeper
    replace the tool's pick with a better-scoring cluster member."""
    genomes = _make_genomes(tmp_path / "genomes", 8)
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )

    selection = tmp_path / "selection.tsv"
    _write_two_genome_selection(selection, genomes)

    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
        ),
        _LOG,
    )

    # genomes[2] (99.0/0.5, score 96.5) outscores the tool's pick genomes[0]
    # (80.0/5.0, score 55.0); genomes[4]'s cluster has no quality data, so the
    # tool's pick stands there.
    assert {r.name for r in final.representatives} == {genomes[2].name, genomes[4].name}
    clusters = read_clusters(tmp_path / "merged" / CLUSTERS_TSV)
    assert sorted(clusters[genomes[2].name]) == sorted(
        [genomes[0].name, genomes[1].name, genomes[3].name]
    )
    assert (tmp_path / "merged" / "representatives" / genomes[2].name).exists()


def test_merge_keeper_tool_ignores_selection_tsv(tmp_path: Path, reg) -> None:
    """keeper='tool' keeps the adapter's own pick even with a selection.tsv."""
    genomes = _make_genomes(tmp_path / "genomes", 8)
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )

    selection = tmp_path / "selection.tsv"
    _write_two_genome_selection(selection, genomes)

    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
            keeper="tool",
        ),
        _LOG,
    )

    assert {r.name for r in final.representatives} == {genomes[0].name, genomes[4].name}


def test_chunk_promotes_best_quality_contained_member(tmp_path: Path, reg) -> None:
    """A selection.tsv with quality columns lets the chunk step's keeper
    replace the tool's pick with a better-scoring cluster member. Every genome
    in the chunk has a real file (params.genomes), so any promotion here is
    always resolvable -- unlike at the merge step (see the deep-member test
    below)."""
    genomes = _make_genomes(tmp_path / "genomes", 4)  # halver: reps [g0, g2]; g0's member g1

    selection = tmp_path / "selection.tsv"
    write_selection(
        selection,
        [
            SelectionRow(
                accession="GCF_000000.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=genomes[0].name,
                completeness=70.0,
                contamination=5.0,
            ),
            SelectionRow(
                accession="GCF_000001.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=genomes[1].name,
                completeness=99.0,
                contamination=0.2,
            ),
        ],
    )

    result = dereplicate_chunk(
        ChunkParams(
            tool="halver",
            genomes=genomes,
            out_dir=tmp_path / "c0",
            selection_tsv=selection,
        ),
        _LOG,
    )

    # genomes[1] (99.0/0.2, score 98.0) outscores the tool's pick genomes[0]
    # (70.0/5.0, score 45.0) and is promoted.
    assert {r.name for r in result.representatives} == {genomes[1].name, genomes[2].name}
    assert (tmp_path / "c0" / "representatives" / genomes[1].name).exists()
    assert not (tmp_path / "c0" / "representatives" / genomes[0].name).exists()
    clusters = read_clusters(tmp_path / "c0" / CLUSTERS_TSV)
    assert clusters[genomes[1].name] == [genomes[0].name]


def test_chunk_keeper_tool_ignores_selection_tsv(tmp_path: Path, reg) -> None:
    """keeper='tool' (the chunk default's counterpart) keeps the adapter's own
    pick even with a selection.tsv."""
    genomes = _make_genomes(tmp_path / "genomes", 4)
    selection = tmp_path / "selection.tsv"
    write_selection(
        selection,
        [
            SelectionRow(
                accession="GCF_000001.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=genomes[1].name,
                completeness=99.0,
                contamination=0.2,
            ),
        ],
    )

    result = dereplicate_chunk(
        ChunkParams(
            tool="halver",
            genomes=genomes,
            out_dir=tmp_path / "c0",
            selection_tsv=selection,
            keeper="tool",
        ),
        _LOG,
    )

    assert {r.name for r in result.representatives} == {genomes[0].name, genomes[2].name}


def test_merge_ignores_quality_for_unresolvable_chunk_member(tmp_path: Path, reg) -> None:
    """A selection.tsv scoring a chunk-level CONTAINED member (never staged at
    the merge step -- only chunk representatives/ files are) must not be
    promoted there: _compose_two_stage's expanded membership includes such
    genomes in final.clusters, but their files live nowhere the merge step can
    see. Reproduces the crash the reviewer found (rescore_representatives would
    pick this genome, then _write_step_contract raised "Representative genome
    file missing") and proves the merge-level resolvable-name filter fixes it
    by leaving the stage-2 pick in place instead."""
    genomes = _make_genomes(tmp_path / "genomes", 8)
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )

    selection = tmp_path / "selection.tsv"
    write_selection(
        selection,
        [
            SelectionRow(
                # genomes[1]: contained under genomes[0] inside chunk c0 --
                # never a chunk representative, so no file reaches the merge step.
                accession="GCF_000001.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=genomes[1].name,
                completeness=99.9,
                contamination=0.0,  # best possible score
            ),
        ],
    )

    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
        ),
        _LOG,
    )

    # genomes[1] cannot be promoted (unresolvable at merge time); the stage-2
    # pick (genomes[0]) stands, and the merge completes without error.
    assert {r.name for r in final.representatives} == {genomes[0].name, genomes[4].name}
    clusters = read_clusters(tmp_path / "merged" / CLUSTERS_TSV)
    assert genomes[1].name in clusters[genomes[0].name]
    assert (tmp_path / "merged" / "representatives" / genomes[0].name).exists()


def test_merge_rejects_chunk_missing_genome_status(tmp_path: Path, reg) -> None:
    """_load_chunk must refuse a chunk directory missing genome_status.tsv, the
    same way it already refuses one missing clusters.tsv."""
    from repgenr.core.errors import WorkdirError

    genomes = _make_genomes(tmp_path / "genomes", 4)
    dereplicate_chunk(ChunkParams(tool="halver", genomes=genomes, out_dir=tmp_path / "c0"), _LOG)
    (tmp_path / "c0" / GENOME_STATUS_TSV).unlink()

    with pytest.raises(WorkdirError):
        dereplicate_merge(
            MergeParams(tool="halver", chunk_dirs=[tmp_path / "c0"], out_dir=tmp_path / "merged"),
            _LOG,
        )


# --- merge-level --reduce and --target-reps ---------------------------------


class _KeepAll(Dereplicator):
    """Every genome is its own representative, so a reduction is observable."""

    capabilities = ToolCapabilities(name="keepall", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"keepall": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        return DerepResult(
            representatives=genomes,
            clusters={g.name: [] for g in genomes},
            genome_status={g.name: STATUS_REPRESENTATIVE for g in genomes},
        )


class _AniDep(Dereplicator):
    """Representative count scales with the secondary ANI: keep = round(ani * N)."""

    capabilities = ToolCapabilities(name="anidep", supports_native_scaling=True)

    def preflight(self) -> dict[str, str]:
        return {"anidep": "1.0"}

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        keep = max(1, min(len(genomes), round(params.secondary_ani * len(genomes))))
        reps = genomes[:keep]
        leftover = [g.name for g in genomes[keep:]]
        clusters: dict[str, list[str]] = {r.name: [] for r in reps}
        clusters[reps[0].name] = leftover
        status = {r.name: STATUS_REPRESENTATIVE for r in reps}
        status.update({m: STATUS_CONTAINED for m in leftover})
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


_TAXA = [
    ("GCF_000001.1", "Fam_aaa_sp1_GCF_000001.1.fasta", "aaa", "aaa-sp1"),
    ("GCF_000002.1", "Fam_aaa_sp1_GCF_000002.1.fasta", "aaa", "aaa-sp1"),
    ("GCF_000003.1", "Fam_aaa_sp2_GCF_000003.1.fasta", "aaa", "aaa-sp2"),
    ("GCF_000004.1", "Fam_bbb_sp3_GCF_000004.1.fasta", "bbb", "bbb-sp3"),
]


def _taxa_chunk(tmp_path: Path, register_tool) -> Path:
    register_tool(registry, "keepall", _KeepAll)
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    genomes = []
    for _acc, fn, _g, _s in _TAXA:
        (gdir / fn).write_text(">x\nACGT\n")
        genomes.append(gdir / fn)
    dereplicate_chunk(ChunkParams(tool="keepall", genomes=genomes, out_dir=tmp_path / "c0"), _LOG)
    return tmp_path / "c0"


def test_merge_reduce_species_uses_selection_taxonomy(tmp_path: Path, register_tool) -> None:
    chunk = _taxa_chunk(tmp_path, register_tool)
    selection = tmp_path / "selection.tsv"
    write_selection(
        selection,
        [SelectionRow(acc, "Fam", g, s, False, fn) for acc, fn, g, s in _TAXA],
    )
    final = dereplicate_merge(
        MergeParams(
            tool="keepall",
            chunk_dirs=[chunk],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
            reduce="species",
        ),
        _LOG,
    )
    assert len(final.representatives) == 3  # sp1 collapsed to one keeper
    status = read_genome_status(tmp_path / "merged" / GENOME_STATUS_TSV)
    assert sorted(status.values()).count(STATUS_CONTAINED) == 1


def test_merge_reduce_keeper_tool_ignores_selection_quality(tmp_path: Path, register_tool) -> None:
    """--keeper tool: the --reduce keeper is the largest cluster, as in the stage.

    The merge step used to rank the representatives by selection.tsv quality
    whatever --keeper said, so a scored genome displaced an unscored one that
    stood for more of the set.
    """
    register_tool(registry, "keepall", _KeepAll)
    register_tool(registry, "anidep", _AniDep)
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    genomes = []
    for _acc, fn, _g, _s in _TAXA[:3]:
        (gdir / fn).write_text(">x\nACGT\n")
        genomes.append(gdir / fn)
    # Chunk 0 folds genome 2 into genome 1 (a cluster of two); chunk 1 keeps
    # genome 3 alone. Genomes 1 and 3 are given the same species below.
    dereplicate_chunk(
        ChunkParams(tool="anidep", genomes=genomes[:2], out_dir=tmp_path / "c0", secondary_ani=0.5),
        _LOG,
    )
    dereplicate_chunk(
        ChunkParams(tool="keepall", genomes=genomes[2:], out_dir=tmp_path / "c1"), _LOG
    )
    selection = tmp_path / "selection.tsv"
    rows = [SelectionRow(acc, "Fam", g, "aaa-sp1", False, fn) for acc, fn, g, _s in _TAXA[:3]]
    # Only the singleton is scored, and high quality, so the quality rule
    # would promote it over the unscored larger cluster.
    rows[2].completeness, rows[2].contamination = 95.0, 1.0
    write_selection(selection, rows)

    def merged(keeper: str) -> list[str]:
        final = dereplicate_merge(
            MergeParams(
                tool="keepall",
                chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
                out_dir=tmp_path / f"m-{keeper}",
                selection_tsv=selection,
                keeper=keeper,
                reduce="species",
            ),
            _LOG,
        )
        return [r.name for r in final.representatives]

    assert merged("tool") == [_TAXA[0][1]]
    assert merged("quality") == [_TAXA[2][1]]


def test_merge_summary_scores_chunk_members_with_the_chunk_n50(
    tmp_path: Path, register_tool
) -> None:
    """The merge step receives only chunk representatives; the N50 of the other
    scored genomes comes from the chunk's genome_n50.tsv."""
    from repgenr.core.contracts import read_cluster_summary
    from repgenr.stages.derep_steps import GENOME_N50_TSV

    register_tool(registry, "anidep", _AniDep)
    register_tool(registry, "keepall", _KeepAll)
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    first, second = (gdir / _TAXA[0][1], gdir / _TAXA[1][1])
    first.write_text(">a\n" + "A" * 10 + "\n>b\n" + "C" * 10 + "\n")  # N50 10
    second.write_text(">a\n" + "A" * 1000 + "\n")  # N50 1000
    selection = tmp_path / "selection.tsv"
    rows = [SelectionRow(acc, "Fam", g, s, False, fn) for acc, fn, g, s in _TAXA[:2]]
    for row in rows:
        row.completeness, row.contamination = 99.0, 0.0
    write_selection(selection, rows)
    dereplicate_chunk(
        ChunkParams(
            tool="anidep",
            genomes=[first, second],
            out_dir=tmp_path / "c0",
            secondary_ani=0.5,
            selection_tsv=selection,
            keeper="tool",
        ),
        _LOG,
    )
    assert (tmp_path / "c0" / GENOME_N50_TSV).read_text().splitlines()[1:] == [
        f"{first.name}\t10",
        f"{second.name}\t1000",
    ]
    shutil.rmtree(gdir)  # the merge must not need the member's file
    dereplicate_merge(
        MergeParams(
            tool="keepall",
            chunk_dirs=[tmp_path / "c0"],
            out_dir=tmp_path / "m",
            selection_tsv=selection,
            keeper="tool",
        ),
        _LOG,
    )
    (row,) = read_cluster_summary(tmp_path / "m" / CLUSTER_SUMMARY_TSV)
    assert (row.representative, row.rep_n50) == (first.name, 10)
    assert row.best_member == second.name
    assert row.best_score == 100.5  # 99 + 0.5 x log10(1000)


def test_merge_reduce_genus_falls_back_to_filename_taxonomy(tmp_path: Path, register_tool) -> None:
    chunk = _taxa_chunk(tmp_path, register_tool)
    final = dereplicate_merge(
        MergeParams(tool="keepall", chunk_dirs=[chunk], out_dir=tmp_path / "m", reduce="genus"),
        _LOG,
    )
    assert len(final.representatives) == 2  # aaa and bbb


def test_merge_target_reps_searches_secondary_ani(tmp_path: Path, register_tool) -> None:
    register_tool(registry, "anidep", _AniDep)
    register_tool(registry, "keepall", _KeepAll)
    genomes = _make_genomes(tmp_path / "genomes", 10)
    dereplicate_chunk(ChunkParams(tool="keepall", genomes=genomes, out_dir=tmp_path / "c0"), _LOG)
    final = dereplicate_merge(
        MergeParams(
            tool="anidep", chunk_dirs=[tmp_path / "c0"], out_dir=tmp_path / "m", target_reps=9
        ),
        _LOG,
    )
    # anidep keeps round(ani * 10); the search must settle on an ANI that yields 9.
    assert len(final.representatives) == 9


def test_chunk_and_merge_summary_take_species_from_selection_tsv(tmp_path: Path, reg) -> None:
    """Non-canonical filenames carry no species; selection.tsv supplies it."""
    from repgenr.core.contracts import read_cluster_summary

    gdir = tmp_path / "genomes"
    gdir.mkdir()
    genomes = []
    for i in range(4):
        p = gdir / f"iso-{i}.fasta"
        p.write_text(">x\nACGT\n")
        genomes.append(p)
    selection = tmp_path / "selection.tsv"
    species = ["tularensis", "tularensis", "novicida", "holarctica"]
    write_selection(
        selection,
        [
            SelectionRow(f"A{i}", "Fam", "Francisella", f"Francisella {species[i]}", False, g.name)
            for i, g in enumerate(genomes)
        ],
    )
    dereplicate_chunk(
        ChunkParams(
            tool="halver", genomes=genomes, out_dir=tmp_path / "c0", selection_tsv=selection
        ),
        _LOG,
    )
    # halver of 4 -> reps [g0, g2]; g0 holds g1, g2 holds g3.
    chunk_rows = {
        r.representative: r for r in read_cluster_summary(tmp_path / "c0" / CLUSTER_SUMMARY_TSV)
    }
    assert (chunk_rows["iso-0.fasta"].n_species, chunk_rows["iso-0.fasta"].species) == (
        1,
        "tularensis",
    )
    assert chunk_rows["iso-2.fasta"].species == "novicida,holarctica"

    dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0"],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
        ),
        _LOG,
    )
    merged = {
        r.representative: r for r in read_cluster_summary(tmp_path / "merged" / CLUSTER_SUMMARY_TSV)
    }
    assert all(r.n_species >= 1 for r in merged.values()), merged
    # Without the selection table the same names carry no species.
    dereplicate_chunk(ChunkParams(tool="halver", genomes=genomes, out_dir=tmp_path / "bare"), _LOG)
    bare = read_cluster_summary(tmp_path / "bare" / CLUSTER_SUMMARY_TSV)
    assert all(r.n_species == 0 for r in bare)


class _QualityRecorder(_Halver):
    capabilities = ToolCapabilities(name="qualityrecorder", supports_native_scaling=True)
    seen: list[dict] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        type(self).seen.append(dict(params.quality))
        return super().dereplicate(genomes, out_dir, params, logger)


def _selection_with_quality(path: Path, genomes: list[Path]) -> dict[str, tuple[float, float]]:
    quality = {g.name: (90.0 + i, 0.5 * i) for i, g in enumerate(genomes)}
    write_selection(
        path,
        [
            SelectionRow(
                accession=f"GCF_{i:06d}.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=g.name,
                completeness=quality[g.name][0],
                contamination=quality[g.name][1],
            )
            for i, g in enumerate(genomes)
        ],
    )
    return quality


def test_chunk_and_merge_pass_selection_quality_to_the_adapter(tmp_path: Path, reg) -> None:
    """The steps hand selection.tsv quality to the adapter, as the stage does
    with the manifest; the merge step passes the union's genomes only."""
    registry.register("qualityrecorder", _QualityRecorder, replace=True)
    _QualityRecorder.seen = []
    genomes = _make_genomes(tmp_path / "genomes", 4)
    selection = tmp_path / "selection.tsv"
    quality = _selection_with_quality(selection, genomes)
    try:
        dereplicate_chunk(
            ChunkParams(
                tool="qualityrecorder",
                genomes=genomes,
                out_dir=tmp_path / "chunk0",
                selection_tsv=selection,
            ),
            _LOG,
        )
        dereplicate_merge(
            MergeParams(
                tool="qualityrecorder",
                chunk_dirs=[tmp_path / "chunk0"],
                out_dir=tmp_path / "merged",
                selection_tsv=selection,
            ),
            _LOG,
        )
    finally:
        registry._classes.pop("qualityrecorder", None)
    chunk_seen, merge_seen = _QualityRecorder.seen
    assert chunk_seen == quality
    # halver keeps genomes 0 and 2 in the chunk; only they reach the merge
    assert merge_seen == {g.name: quality[g.name] for g in (genomes[0], genomes[2])}


def test_partial_selection_quality_reaches_no_chunk(tmp_path: Path, reg) -> None:
    """The steps decide over the whole selection: one unscored genome in
    another chunk means no chunk is given quality."""
    registry.register("qualityrecorder", _QualityRecorder, replace=True)
    _QualityRecorder.seen = []
    genomes = _make_genomes(tmp_path / "genomes", 4)
    selection = tmp_path / "selection.tsv"
    _selection_with_quality(selection, genomes)
    text = selection.read_text(encoding="utf-8").splitlines()
    last = text[-1].split("\t")
    header = text[0].split("\t")
    last[header.index("completeness")] = ""
    last[header.index("contamination")] = ""
    selection.write_text("\n".join([*text[:-1], "\t".join(last)]) + "\n", encoding="utf-8")
    try:
        dereplicate_chunk(
            ChunkParams(
                tool="qualityrecorder",
                genomes=genomes[:2],  # both scored
                out_dir=tmp_path / "chunk0",
                selection_tsv=selection,
            ),
            _LOG,
        )
    finally:
        registry._classes.pop("qualityrecorder", None)
    assert _QualityRecorder.seen == [{}]


def _gtdb_selection(path: Path, genomes: list[Path], reps: set[int]) -> None:
    """selection.tsv flagging ``genomes[i]`` for i in ``reps`` as GTDB species
    representatives; genomes[0] carries the best quality so the quality rule
    alone would keep it."""
    write_selection(
        path,
        [
            SelectionRow(
                accession=f"GCA_{i:06d}.1",
                family="Fam",
                genus="g",
                species="s",
                is_outgroup=False,
                filename=g.name,
                completeness=99.0 if i == 0 else 90.0,
                contamination=0.1 if i == 0 else 2.0,
                gtdb_representative=i in reps,
            )
            for i, g in enumerate(genomes)
        ],
    )


def test_chunk_keeper_gtdb_keeps_the_flagged_member(tmp_path: Path, reg) -> None:
    genomes = _make_genomes(tmp_path / "genomes", 4)  # halver: g0 holds g1, g2 holds g3
    selection = tmp_path / "selection.tsv"
    _gtdb_selection(selection, genomes, {1})
    result = dereplicate_chunk(
        ChunkParams(
            tool="halver",
            genomes=genomes,
            out_dir=tmp_path / "c0",
            selection_tsv=selection,
            keeper="gtdb",
        ),
        _LOG,
    )
    assert {r.name for r in result.representatives} == {genomes[1].name, genomes[2].name}
    assert (tmp_path / "c0" / "representatives" / genomes[1].name).exists()
    summary = (tmp_path / "c0" / CLUSTER_SUMMARY_TSV).read_text(encoding="utf-8").splitlines()
    column = summary[0].split("\t").index("rep_is_gtdb_representative")
    flags = {line.split("\t")[0]: line.split("\t")[column] for line in summary[1:]}
    assert flags == {genomes[1].name: "1", genomes[2].name: "0"}


def test_merge_keeper_gtdb_keeps_the_flagged_chunk_representative(tmp_path: Path, reg) -> None:
    """At the merge, g0 (best quality) and g2 meet in one cluster; g2 is the
    GTDB representative and is kept."""
    genomes = _make_genomes(tmp_path / "genomes", 8)
    selection = tmp_path / "selection.tsv"
    _gtdb_selection(selection, genomes, {2})
    for name, part in (("c0", genomes[:4]), ("c1", genomes[4:])):
        dereplicate_chunk(
            ChunkParams(
                tool="halver",
                genomes=part,
                out_dir=tmp_path / name,
                selection_tsv=selection,
                keeper="gtdb",
            ),
            _LOG,
        )
    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
            keeper="gtdb",
        ),
        _LOG,
    )
    assert genomes[2].name in {r.name for r in final.representatives}
    assert genomes[0].name not in {r.name for r in final.representatives}
    assert (tmp_path / "merged" / "representatives" / genomes[2].name).exists()


def test_merge_keeper_gtdb_without_flags_falls_back_to_quality(tmp_path: Path, reg) -> None:
    genomes = _make_genomes(tmp_path / "genomes", 8)
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[:4], out_dir=tmp_path / "c0"), _LOG
    )
    dereplicate_chunk(
        ChunkParams(tool="halver", genomes=genomes[4:], out_dir=tmp_path / "c1"), _LOG
    )
    selection = tmp_path / "selection.tsv"
    _write_two_genome_selection(selection, genomes)
    final = dereplicate_merge(
        MergeParams(
            tool="halver",
            chunk_dirs=[tmp_path / "c0", tmp_path / "c1"],
            out_dir=tmp_path / "merged",
            selection_tsv=selection,
            keeper="gtdb",
        ),
        _LOG,
    )
    assert {r.name for r in final.representatives} == {genomes[2].name, genomes[4].name}


def test_chunk_keeper_gtdb_warns_when_only_the_outgroup_is_flagged(
    tmp_path: Path, reg, caplog
) -> None:
    import logging

    genomes = _make_genomes(tmp_path / "genomes", 4)
    selection = tmp_path / "selection.tsv"
    _gtdb_selection(selection, genomes, set())
    rows = read_selection(selection)
    rows.append(
        SelectionRow(
            accession="GCA_000099.1",
            family="Fam",
            genus="g",
            species="other",
            is_outgroup=True,
            filename="Fam_g_other_GCA_000099.1.fasta",
            gtdb_representative=True,
        )
    )
    write_selection(selection, rows)
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        dereplicate_chunk(
            ChunkParams(
                tool="halver",
                genomes=genomes,
                out_dir=tmp_path / "c0",
                selection_tsv=selection,
                keeper="gtdb",
            ),
            _LOG,
        )
    assert any("flags no GTDB" in r.getMessage() for r in caplog.records)


# --- dereplicate-refine ---------------------------------------------------------


_TIERED_NAMES = [
    f"Francisellaceae_francisella_tularensis_GCA_000{a}.fasta"
    for a in (111, 112, 121, 211, 212, 221, 311)
]


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


def test_refine_step_refines_a_merged_contract(tmp_path: Path, tiered) -> None:
    """The refine step takes a contract directory and every genome file,
    refines the seeds' clusters at their ANIs and writes the contract shape
    plus refined.tsv; a second pass keeps the first's rows."""
    from repgenr.core.contracts import REFINED_TSV, read_refined
    from repgenr.core.errors import UserInputError, WorkdirError
    from repgenr.stages.derep_steps import RefineStepParams, dereplicate_refine

    gdir = tmp_path / "genomes"
    gdir.mkdir()
    genomes = []
    for i, n in enumerate(_TIERED_NAMES):
        (gdir / n).write_text(f">s{i}\n{'ACGT' * (10 + i)}\n", encoding="utf-8")
        genomes.append(gdir / n)
    merged = dereplicate_chunk(
        ChunkParams(tool="tiered", genomes=genomes, out_dir=tmp_path / "merged"), _LOG
    )
    assert len(merged.representatives) == 3
    refine_file = tmp_path / "refine.tsv"
    refine_file.write_text("seed\tani\nGCA_000212\t0.9995\n", encoding="utf-8")
    _Tiered.calls = []
    out = dereplicate_refine(
        RefineStepParams(
            tool="tiered",
            contract_dir=tmp_path / "merged",
            genomes=genomes,
            out_dir=tmp_path / "refined",
            seeds=[_TIERED_NAMES[2]],
            ani=0.999,
            refine_file=str(refine_file),
            versions_out=tmp_path / "v.yml",
        ),
        _LOG,
    )
    assert _Tiered.calls == [(3, 0.999), (3, 0.9995)]
    assert sorted(p.name for p in out.representatives) == [
        _TIERED_NAMES[0],
        _TIERED_NAMES[2],
        _TIERED_NAMES[3],
        _TIERED_NAMES[5],
        _TIERED_NAMES[6],
    ]
    assert sorted(p.name for p in (tmp_path / "refined" / "representatives").iterdir()) == sorted(
        p.name for p in out.representatives
    )
    rows = {r.genome: r for r in read_refined(tmp_path / "refined" / REFINED_TSV)}
    assert rows[_TIERED_NAMES[2]].ani == "0.999" and rows[_TIERED_NAMES[5]].ani == "0.9995"
    summary = (tmp_path / "refined" / CLUSTER_SUMMARY_TSV).read_text(encoding="utf-8").splitlines()
    column = summary[0].split("\t").index("refined")
    assert {line.split("\t")[0]: line.split("\t")[column] for line in summary[1:]}[
        _TIERED_NAMES[6]
    ] == "0"
    assert (tmp_path / "v.yml").exists()

    # Refining the refined contract again for C keeps A's and B's rows.
    _Tiered.calls = []
    dereplicate_refine(
        RefineStepParams(
            tool="tiered",
            contract_dir=tmp_path / "refined",
            genomes=genomes,
            out_dir=tmp_path / "refined2",
            seeds=[_TIERED_NAMES[0]],
            ani=0.9995,
        ),
        _LOG,
    )
    assert _Tiered.calls == [(3, 0.9995)]
    assert len(read_refined(tmp_path / "refined2" / REFINED_TSV)) == 6

    # Refusals: a missing member file, nothing to refine, an ANI not finer.
    with pytest.raises(WorkdirError, match="not among --genomes-fofn"):
        dereplicate_refine(
            RefineStepParams(
                tool="tiered",
                contract_dir=tmp_path / "merged",
                genomes=genomes[:3],
                out_dir=tmp_path / "x",
                seeds=[_TIERED_NAMES[0]],
                ani=0.999,
            ),
            _LOG,
        )
    with pytest.raises(UserInputError, match="not finer"):
        dereplicate_refine(
            RefineStepParams(
                tool="tiered",
                contract_dir=tmp_path / "merged",
                genomes=genomes,
                out_dir=tmp_path / "x",
                seeds=[_TIERED_NAMES[0]],
                ani=0.95,
            ),
            _LOG,
        )
