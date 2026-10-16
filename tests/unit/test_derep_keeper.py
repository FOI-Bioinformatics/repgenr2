"""Quality-aware keeper: the best-scoring cluster member becomes the representative."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)
from repgenr.stages.derep_keeper import (
    N50Lookup,
    choose_keeper,
    genome_n50,
    keeper_score,
    quality_score,
    rescore_representatives,
)

_LOG = logging.getLogger("keeper")


def _result() -> DerepResult:
    return DerepResult(
        representatives=[Path("/g/rep.fasta"), Path("/g/solo.fasta")],
        clusters={"rep.fasta": ["m1.fasta", "m2.fasta"], "solo.fasta": []},
        genome_status={
            "rep.fasta": STATUS_REPRESENTATIVE,
            "solo.fasta": STATUS_REPRESENTATIVE,
            "m1.fasta": STATUS_CONTAINED,
            "m2.fasta": STATUS_CONTAINED,
        },
    )


def test_score_penalises_contamination() -> None:
    assert quality_score(100.0, 0.0) > quality_score(100.0, 2.0)
    assert quality_score(95.0, 0.0) == 95.0


def test_better_member_replaces_representative() -> None:
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2), "m2.fasta": (95.0, 1.0)}
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 1
    assert sorted(p.name for p in out.representatives) == ["m1.fasta", "solo.fasta"]
    assert out.clusters["m1.fasta"] == ["m2.fasta", "rep.fasta"]
    assert out.genome_status["m1.fasta"] == STATUS_REPRESENTATIVE
    assert out.genome_status["rep.fasta"] == STATUS_CONTAINED
    assert out.representatives[0].parent == Path("/g")
    check_result_complete(out, ["rep.fasta", "solo.fasta", "m1.fasta", "m2.fasta"])


def test_high_quality_member_replaces_unscored_representative() -> None:
    """An unscored representative is replaced by a scored member only when
    that member is high quality (completeness > 90, contamination < 5)."""
    quality = {"m1.fasta": (95.0, 2.0)}  # rep.fasta and m2.fasta carry no quality
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 1
    assert sorted(p.name for p in out.representatives) == ["m1.fasta", "solo.fasta"]
    assert out.clusters["m1.fasta"] == ["m2.fasta", "rep.fasta"]
    assert out.genome_status["m1.fasta"] == STATUS_REPRESENTATIVE
    assert out.genome_status["rep.fasta"] == STATUS_CONTAINED


def test_scored_fragment_does_not_replace_unscored_representative(caplog) -> None:
    """A 40 percent fragment with quality must not displace a genome that
    merely lacks quality; the cluster is named in an INFO line."""
    quality = {"m1.fasta": (40.0, 0.0)}
    with caplog.at_level(logging.INFO, logger="keeper"):
        out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 0
    assert out.clusters == _result().clusters
    assert any("rep.fasta" in r.message and "high quality" in r.message for r in caplog.records)


def test_unscored_member_never_wins() -> None:
    quality = {"rep.fasta": (90.0, 3.0)}
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 0
    assert out.clusters == _result().clusters


def test_no_quality_keeps_adapter_choice() -> None:
    out, swaps = rescore_representatives(_result(), {}, _LOG)
    assert swaps == 0
    assert out == _result()


def test_equal_scores_are_broken_by_completeness_then_name() -> None:
    """Ties do not depend on which genome the tool picked: higher completeness
    wins (95/0 and 100/1 both score 95), then the name."""
    tie = {"rep.fasta": (95.0, 0.0), "m1.fasta": (100.0, 1.0), "m2.fasta": (50.0, 0.0)}
    out, _ = rescore_representatives(_result(), tie, _LOG)
    assert "m1.fasta" in out.clusters
    same = {"rep.fasta": (99.0, 0.0), "m1.fasta": (99.0, 0.0), "m2.fasta": (50.0, 0.0)}
    out, _ = rescore_representatives(_result(), same, _LOG)
    assert "m1.fasta" in out.clusters  # "m1.fasta" sorts before "rep.fasta"


def test_n50_term_prefers_the_contiguous_genome() -> None:
    """A closed genome keeps the cluster against a draft 0.15 points better on
    CheckM (F. tularensis: 99.85 at N50 1.9 Mb against 100.0 at 148 kb)."""
    quality = {"rep.fasta": (100.0, 0.03), "m1.fasta": (100.0, 0.0)}
    n50 = {"rep.fasta": 1_910_592, "m1.fasta": 147_757}.get
    out, swaps = rescore_representatives(_result(), quality, _LOG, n50)
    assert swaps == 0
    out, swaps = rescore_representatives(_result(), quality, _LOG)  # no N50 term
    assert swaps == 1


def test_keeper_score_adds_half_log10_n50() -> None:
    assert keeper_score(100.0, 0.0, 1_000_000) == 103.0
    assert keeper_score(100.0, 0.0, None) == 100.0


def test_genome_n50_reads_plain_and_gzip(tmp_path: Path) -> None:
    import gzip

    text = ">a\n" + "A" * 50 + "\n" + "A" * 50 + "\n>b\nCCCC\n>c\n" + "G" * 60 + "\n"
    plain = tmp_path / "g.fasta"
    plain.write_text(text)
    packed = tmp_path / "g.fasta.gz"
    packed.write_bytes(gzip.compress(text.encode()))
    # lengths 100, 60, 4 (total 164): half 82 is reached by the first contig
    assert genome_n50(plain) == genome_n50(packed) == 100


def test_n50_lookup_reads_each_genome_once(tmp_path: Path, monkeypatch) -> None:
    from repgenr.stages import derep_keeper

    (tmp_path / "x.fasta").write_text(">x\nACGT\n")
    calls: list[Path] = []
    real = derep_keeper.genome_n50
    monkeypatch.setattr(derep_keeper, "genome_n50", lambda p: calls.append(p) or real(p))
    lookup = N50Lookup([tmp_path / "missing", tmp_path], known={"k.fasta": 7})
    assert lookup("x.fasta") == lookup("x.fasta") == 4
    assert lookup("k.fasta") == 7
    assert lookup("absent.fasta") is None
    assert len(calls) == 1
    assert lookup.computed() == {"k.fasta": 7, "x.fasta": 4}


def test_choose_keeper_ignores_unscored_genomes() -> None:
    assert choose_keeper("a", ["b", "c"], {"a": (90.0, 3.0)}) == "a"
    assert choose_keeper("a", ["b"], {}) == "a"


def test_fail_qc_genome_survives_rescore_untouched() -> None:
    """A genome with no cluster (rejected on QC upstream) keeps its status and
    still lets the rescored result pass check_result_complete."""
    result = _result()
    result.genome_status["bad.fasta"] = STATUS_FAIL_QC
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2), "m2.fasta": (95.0, 1.0)}

    out, swaps = rescore_representatives(result, quality, _LOG)

    assert swaps == 1
    assert out.genome_status["bad.fasta"] == STATUS_FAIL_QC
    check_result_complete(out, ["rep.fasta", "solo.fasta", "m1.fasta", "m2.fasta", "bad.fasta"])


def test_a_single_competitor_reads_no_n50() -> None:
    """Scored singletons and clusters with one scored genome compare nothing,
    so no genome file is read."""

    def no_read(name: str) -> int | None:
        raise AssertionError(f"N50 read for {name}")

    assert choose_keeper("a", [], {"a": (99.0, 0.0)}, no_read) == "a"
    assert choose_keeper("a", ["b"], {"a": (99.0, 0.0)}, no_read) == "a"
    assert choose_keeper("a", ["b"], {"b": (99.0, 0.0)}, no_read) == "b"
    out, swaps = rescore_representatives(_result(), {"solo.fasta": (99.0, 0.0)}, _LOG, no_read)
    assert swaps == 0


def test_genome_n50_handles_crlf_missing_final_newline_and_leading_text(tmp_path: Path) -> None:
    crlf = tmp_path / "crlf.fasta"
    crlf.write_bytes(b">a desc\r\nAAAA\r\nAAAA\r\n>b\r\nCC")
    assert genome_n50(crlf) == 8  # lengths 8 and 2
    lead = tmp_path / "lead.fasta"
    lead.write_bytes(b"; comment\n>a\nACGTACGT\n>empty\n>c\nAC\n")
    assert genome_n50(lead) == 8


def test_genome_n50_refuses_a_truncated_gzip(tmp_path: Path) -> None:
    import gzip

    import pytest

    from repgenr.core.errors import UserInputError

    path = tmp_path / "g.fasta.gz"
    path.write_bytes(gzip.compress(b">a\n" + b"ACGT" * 1000)[:-20])
    with pytest.raises(UserInputError, match="truncated or corrupt"):
        genome_n50(path)


def _write_genomes(directory: Path, lengths: dict[str, list[int]]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name, contigs in lengths.items():
        body = "".join(f">c{i}\n{'A' * n}\n" for i, n in enumerate(contigs))
        (directory / name).write_text(body)


def test_prefetch_reads_the_compared_genomes_in_parallel(tmp_path: Path, monkeypatch) -> None:
    """The N50s a rescore compares are read before the comparison, on
    several threads at once, and each genome once; a scored singleton and
    the unscored members are not read."""
    import threading
    import time

    from repgenr.stages import derep_keeper

    _write_genomes(
        tmp_path,
        {
            n: [100 * (i + 1)]
            for i, n in enumerate(["rep.fasta", "m1.fasta", "m2.fasta", "solo.fasta"])
        },
    )
    threads: set[int] = set()
    calls: list[str] = []
    real = derep_keeper.genome_n50

    def slow(path: Path) -> int:
        threads.add(threading.get_ident())
        calls.append(path.name)
        time.sleep(0.05)
        return real(path)

    monkeypatch.setattr(derep_keeper, "genome_n50", slow)
    lookup = N50Lookup([tmp_path], workers=4)
    quality = {
        "rep.fasta": (95.0, 1.0),
        "m1.fasta": (95.0, 1.0),
        "m2.fasta": (95.0, 1.0),
        "solo.fasta": (99.0, 0.0),
    }
    result = _result()
    result.clusters["rep.fasta"] = ["m1.fasta", "m2.fasta"]

    out, swaps = rescore_representatives(result, quality, _LOG, lookup)

    assert sorted(calls) == ["m1.fasta", "m2.fasta", "rep.fasta"]
    assert len(threads) > 1
    assert swaps == 1 and "m2.fasta" in out.clusters  # the largest N50


def test_rescore_choices_do_not_depend_on_the_thread_count(tmp_path: Path) -> None:
    """One worker and many give the same representatives and the same log."""
    names = [f"g{i:03d}.fasta" for i in range(60)]
    _write_genomes(tmp_path, {n: [50 + (i * 37) % 400, 30] for i, n in enumerate(names)})
    clusters = {names[i]: names[i + 1 : i + 6] for i in range(0, 60, 6)}
    status = {n: STATUS_REPRESENTATIVE if n in clusters else STATUS_CONTAINED for n in names}
    quality = {n: (92.0 + (i % 3), 1.0) for i, n in enumerate(names) if i % 4}

    def run(workers: int) -> tuple[dict[str, list[str]], int]:
        result = DerepResult(
            representatives=[tmp_path / n for n in clusters],
            clusters={k: list(v) for k, v in clusters.items()},
            genome_status=dict(status),
        )
        out, swaps = rescore_representatives(
            result, quality, _LOG, N50Lookup([tmp_path], workers=workers)
        )
        return out.clusters, swaps

    assert run(1) == run(8)


def test_summary_prefetch_includes_scored_representatives(tmp_path: Path) -> None:
    """The cluster summary reports a scored representative's N50 even when it
    is compared with nothing, so its prefetch reads it."""
    from repgenr.stages.derep_keeper import prefetch_n50

    _write_genomes(tmp_path, {"solo.fasta": [100], "rep.fasta": [10], "m.fasta": [10]})
    lookup = N50Lookup([tmp_path], workers=2)
    groups = [("solo.fasta", []), ("rep.fasta", ["m.fasta"])]
    quality = {"solo.fasta": (99.0, 0.0), "rep.fasta": (99.0, 0.0)}
    prefetch_n50(lookup, groups, quality)
    assert lookup.computed() == {}
    prefetch_n50(lookup, groups, quality, scored_incumbents=True)
    assert lookup.computed() == {"solo.fasta": 100, "rep.fasta": 10}
