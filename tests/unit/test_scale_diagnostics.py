"""Scale-diagnostics guards found by the scaling audit (docs/audit/scaling-audit.md)."""

from __future__ import annotations

import logging

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.core.process import warn_argv_bytes
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepParams,
    DerepResult,
)
from repgenr.stages.dereplicate import DereplicateParams, _search_target_reps
from repgenr.stages.genome import _timeout_for

_LOGGER = logging.getLogger("scale-diagnostics")


# --- argv byte preflight ------------------------------------------------------


def test_warn_argv_bytes_warns_near_arg_max(caplog, monkeypatch):
    from repgenr.core import process

    monkeypatch.setattr(process, "arg_max", lambda: 1 << 20)
    argv = ["tool", *[f"/data/genomes/genome_{i:06d}.fasta" for i in range(40000)]]
    with caplog.at_level(logging.WARNING):
        warn_argv_bytes("mashtree", argv, _LOGGER)
    assert any("ARG_MAX" in rec.message for rec in caplog.records)


def test_warn_argv_bytes_quiet_for_small_argv(caplog):
    with caplog.at_level(logging.WARNING):
        warn_argv_bytes("mashtree", ["tool", "a.fasta", "b.fasta"], _LOGGER)
    assert not caplog.records


def test_argv_over_arg_max_is_refused_before_the_tool_starts(monkeypatch):
    """A command over ARG_MAX fails with advice, not with E2BIG from exec."""
    import subprocess

    import pytest

    from repgenr.core import process
    from repgenr.core.errors import UserInputError

    monkeypatch.setattr(process, "arg_max", lambda: 64 * 1024)

    def must_not_start(*args, **kwargs):
        raise AssertionError("the tool must not be started")

    monkeypatch.setattr(subprocess, "Popen", must_not_start)
    argv = ["mashtree", *[f"/data/genomes/genome_{i:06d}.fasta" for i in range(5000)]]
    with pytest.raises(UserInputError) as exc:
        process.run(argv, logger=_LOGGER, env={"X": "1"})
    assert "--process-size" in str(exc.value) and "ARG_MAX" in str(exc.value)


def test_argv_size_counts_the_environment(monkeypatch):
    import pytest

    from repgenr.core import process
    from repgenr.core.errors import UserInputError

    monkeypatch.setattr(process, "arg_max", lambda: 64 * 1024)
    process.check_argv_size(["tool", "a.fasta"], {"SMALL": "1"})
    with pytest.raises(UserInputError):
        process.check_argv_size(["tool", "a.fasta"], {"BIG": "x" * 60_000})


def test_linux_single_argument_limit(monkeypatch):
    import sys

    import pytest

    from repgenr.core import process
    from repgenr.core.errors import UserInputError

    monkeypatch.setattr(process, "arg_max", lambda: 64 << 20)
    long_arg = "x" * (200 * 1024)
    monkeypatch.setattr(sys, "platform", "darwin")
    process.check_argv_size(["tool", long_arg], {})
    monkeypatch.setattr(sys, "platform", "linux")
    with pytest.raises(UserInputError, match="--process-size"):
        process.check_argv_size(["tool", long_arg], {})


def test_normal_commands_pass_the_argv_check():
    from repgenr.core import process

    assert process.arg_max() > 64 * 1024
    assert process.run(["true"], logger=_LOGGER) == 0


def test_mashtree_reads_genomes_from_a_file_of_files(tmp_path, monkeypatch, caplog):
    """40000 paths never reach argv (live audit: ARG_MAX at ~9500 viral records)."""
    from repgenr.treebuilders import mashtree as mt

    calls = []
    monkeypatch.setattr(mt, "run_tool", lambda caps, cmd, **k: calls.append((cmd, k)))
    genomes = [tmp_path / f"very/long/prefix/genome_{i:06d}.fasta" for i in range(40000)]
    adapter = mt.MashtreeBuilder()
    from repgenr.treebuilders.base import TreeParams

    with caplog.at_level(logging.WARNING):
        adapter.build(genomes, tmp_path / "out", TreeParams(threads=2), _LOGGER)
    assert not any("ARG_MAX" in rec.message for rec in caplog.records)
    ((cmd, kwargs),) = calls
    argv = [str(c) for c in cmd]
    fofn = argv[argv.index("--file-of-files") + 1]
    assert len(open(fofn, encoding="utf-8").read().splitlines()) == 40000
    assert str(tmp_path / "very/long/prefix") in kwargs["extra_mounts"]


# --- sourmash dense memory message --------------------------------------------


def test_dense_refusal_reports_real_gigabytes(tmp_path):
    from repgenr.dereplicators.sourmash import SourmashDereplicator

    genomes = [tmp_path / f"g{i}.fasta" for i in range(5001)]
    with pytest.raises(WorkdirError, match=r"~0\.2 GB"):
        SourmashDereplicator()._dense_dereplicate(
            genomes, tmp_path / "out", 31, 1000, 0.99, _LOGGER
        )


# --- sourmash tree builder measured limit -------------------------------------


def test_sourmash_treebuilder_limit_is_measured():
    from repgenr.treebuilders.sourmash import SourmashBuilder

    # 10000 was an unbenchmarked assertion; the measured n^3 fit extrapolates
    # to ~12 h of pure-Python NJ there (docs/audit/scaling-audit.md).
    assert SourmashBuilder.capabilities.recommended_max_genomes == 2000


def test_python_nj_hard_guard_above_5000(monkeypatch):
    # The guard moved from the sourmash builder to the joining step: without
    # rapidnj, more than 5000 leaves are refused before any matrix is built.
    from repgenr.tree import nj

    monkeypatch.setattr(nj, "tool_available", lambda caps: False)
    with pytest.raises(WorkdirError, match="Neighbour joining of 5001 leaves"):
        nj.choose_backend(nj.NJ_AUTO, 5001, _LOGGER)


# --- target-reps non-convergence warning --------------------------------------


class _StepDerep:
    """Clonal step function: 1 rep below the clone ANI, n at/above."""

    def __init__(self, n: int) -> None:
        self.n = n

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        keep = len(genomes) if params.secondary_ani >= 0.9995 else 1
        reps = genomes[:keep]
        leftover = [g.name for g in genomes[keep:]]
        clusters = {r.name: [] for r in reps}
        clusters[reps[0].name] = leftover
        status = {r.name: STATUS_REPRESENTATIVE for r in reps}
        status.update({m: STATUS_CONTAINED for m in leftover})
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


def test_target_reps_warns_when_target_unreachable(tmp_path, caplog):
    genomes = []
    for i in range(40):
        p = tmp_path / f"g{i:03d}.fasta"
        p.write_text(">x\nACGT\n", encoding="utf-8")
        genomes.append(p)
    with caplog.at_level(logging.WARNING):
        _search_target_reps(
            _StepDerep(40),
            genomes,
            tmp_path / "scratch",
            DerepParams(),
            DereplicateParams(tool="mock"),
            target=20,
            logger=_LOGGER,
        )
    warnings = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("20" in m and "target" in m.lower() for m in warnings)


def test_target_reps_no_warning_on_exact_hit(tmp_path, caplog):
    genomes = []
    for i in range(40):
        p = tmp_path / f"g{i:03d}.fasta"
        p.write_text(">x\nACGT\n", encoding="utf-8")
        genomes.append(p)
    with caplog.at_level(logging.WARNING):
        _search_target_reps(
            _StepDerep(40),
            genomes,
            tmp_path / "scratch",
            DerepParams(),
            DereplicateParams(tool="mock"),
            target=40,
            logger=_LOGGER,
        )
    assert not [r for r in caplog.records if r.levelno == logging.WARNING]


# --- datasets timeout scales with batch size ----------------------------------


def test_datasets_timeout_scales_with_batch():
    assert _timeout_for(10) == 3600.0  # floor for small batches
    assert _timeout_for(5000) == 15000.0  # 3 s per accession at full batch
    assert _timeout_for(0) == 3600.0


# --- mashtree: one argument set for the tree and the distance matrix ----------


def _mashtree_commands(monkeypatch):
    from repgenr.treebuilders import mashtree as mt

    seen: list[list[str]] = []
    monkeypatch.setattr(mt, "run_tool", lambda caps, cmd, **k: seen.append([str(c) for c in cmd]))
    return mt, seen


def test_mashtree_tree_and_matrix_share_flags(tmp_path, monkeypatch):
    from repgenr.treebuilders.base import TreeParams

    mt, seen = _mashtree_commands(monkeypatch)
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    params = TreeParams(threads=4, extra={"genomesize": 1500000})
    adapter = mt.MashtreeBuilder()
    adapter.build(genomes, tmp_path / "tree", params, _LOGGER)
    adapter.distance_matrix(genomes, tmp_path / "dist", params, _LOGGER)
    build_cmd, matrix_cmd = seen

    def flags(cmd: list[str]) -> dict[str, str]:
        out = {}
        for i, tok in enumerate(cmd):
            if tok.startswith("--") and tok != "--outmatrix":
                out[tok] = cmd[i + 1]
        return out

    # threads, mindepth and genomesize are passed identically to both calls
    per_call = {"--outmatrix", "--file-of-files"}  # output paths differ per call
    assert {k: v for k, v in flags(build_cmd).items() if k not in per_call} == {
        k: v for k, v in flags(matrix_cmd).items() if k not in per_call
    }
    assert flags(build_cmd)["--numcpus"] == "4"
    assert flags(build_cmd)["--mindepth"] == "0"
    assert flags(build_cmd)["--genomesize"] == "1500000"


def test_mashtree_build_also_writes_the_distance_matrix(tmp_path, monkeypatch):
    from repgenr.treebuilders.base import TreeParams

    mt, seen = _mashtree_commands(monkeypatch)
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    mt.MashtreeBuilder().build(genomes, tmp_path / "tree", TreeParams(threads=1), _LOGGER)
    (cmd,) = seen
    assert "--outmatrix" in cmd
    assert cmd[cmd.index("--outmatrix") + 1] == str(tmp_path / "tree" / "distance_matrix.tsv")
