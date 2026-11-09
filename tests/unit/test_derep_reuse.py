"""Dereplication reuses what it already has: N50 values kept between runs,
and representatives that did not change (no copy on a file system without
hard links)."""

from __future__ import annotations

import os
import time
from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE, DerepResult
from repgenr.stages import derep_keeper
from repgenr.stages.derep_keeper import N50_CACHE, save_n50_cache, workdir_n50_lookup
from repgenr.stages.dereplicate import _write_contract


def _genome(path: Path, contigs: list[int]) -> Path:
    text = "".join(f">c{i}\n{'A' * n}\n" for i, n in enumerate(contigs))
    path.write_text(text, encoding="utf-8")
    return path


def test_n50_cache_is_reused_until_a_genome_changes(tmp_path: Path, monkeypatch) -> None:
    genomes, scratch = tmp_path / "genomes", tmp_path / "scratch"
    genomes.mkdir()
    a = _genome(genomes / "a.fasta", [100, 50])
    _genome(genomes / "b.fasta", [30, 30, 30])
    first = workdir_n50_lookup(genomes, scratch, workers=1)
    assert first("a.fasta") == 100 and first("b.fasta") == 30
    save_n50_cache(first, genomes, scratch)
    assert (scratch / N50_CACHE).is_file()

    reads: list[Path] = []
    real = derep_keeper.genome_n50
    monkeypatch.setattr(derep_keeper, "genome_n50", lambda p: reads.append(p) or real(p))
    again = workdir_n50_lookup(genomes, scratch, workers=1)
    assert again("a.fasta") == 100 and again("b.fasta") == 30
    assert reads == []  # both from the cache

    # A changed genome (new size and time) is read again.
    time.sleep(0.01)
    _genome(a, [500])
    changed = workdir_n50_lookup(genomes, scratch, workers=1)
    assert changed("a.fasta") == 500 and [p.name for p in reads] == ["a.fasta"]


def test_a_damaged_cache_is_ignored(tmp_path: Path) -> None:
    genomes, scratch = tmp_path / "genomes", tmp_path / "scratch"
    genomes.mkdir()
    scratch.mkdir()
    _genome(genomes / "a.fasta", [7])
    (scratch / N50_CACHE).write_text("not\ta\tcache\n")
    assert workdir_n50_lookup(genomes, scratch, workers=1)("a.fasta") == 7


def _result(genomes: Path, reps: list[str], members: dict[str, list[str]]) -> DerepResult:
    status = {r: STATUS_REPRESENTATIVE for r in reps}
    for ms in members.values():
        status.update({m: STATUS_CONTAINED for m in ms})
    return DerepResult(
        representatives=[genomes / r for r in reps],
        clusters={r: members.get(r, []) for r in reps},
        genome_status=status,
    )


def test_unchanged_representatives_are_not_copied_again(tmp_path: Path, monkeypatch) -> None:
    # As on exFAT: no hard links, every representative is a copy.
    def no_link(src, dst):  # noqa: ANN001
        raise OSError("hard links not supported")

    monkeypatch.setattr(os, "link", no_link)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.genomes_dir.mkdir(parents=True)
    for name in ("a.fasta", "b.fasta", "c.fasta"):
        _genome(ctx.genomes_dir / name, [100])
    _write_contract(
        ctx, _result(ctx.genomes_dir, ["a.fasta", "b.fasta"], {"a.fasta": ["c.fasta"]}), {}
    )
    rep_a = ctx.representatives_dir / "a.fasta"
    ino_a = rep_a.stat().st_ino
    assert sorted(p.name for p in ctx.representatives_dir.iterdir()) == ["a.fasta", "b.fasta"]

    # A refinement keeps a, drops b and adds c: only c is copied.
    _write_contract(
        ctx, _result(ctx.genomes_dir, ["a.fasta", "c.fasta"], {"a.fasta": ["b.fasta"]}), {}
    )
    assert sorted(p.name for p in ctx.representatives_dir.iterdir()) == ["a.fasta", "c.fasta"]
    assert rep_a.stat().st_ino == ino_a

    # A changed source genome is copied again.
    time.sleep(0.03)
    _genome(ctx.genomes_dir / "a.fasta", [100, 100])
    _write_contract(ctx, _result(ctx.genomes_dir, ["a.fasta", "c.fasta"], {}), {})
    assert rep_a.read_text(encoding="utf-8").count(">") == 2
