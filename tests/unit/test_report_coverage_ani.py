"""Member-to-representative ANI from sketches, and reference coverage of the
simple SNP typer, in the report."""

from __future__ import annotations

import gzip
import json
import logging
import zipfile
from pathlib import Path

import pytest

from repgenr.core.contracts import write_clusters, write_genome_status
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.member_ani import KSIZE, load_hashes, member_anis, pair_ani
from repgenr.reporting.registry import BY_NAME
from repgenr.snptypers.simple import _log_coverage, _write_coverage


def _sig(path: Path, hashes, ksize: int = KSIZE) -> None:
    payload = [{"signatures": [{"ksize": ksize, "molecule": "DNA", "mins": sorted(hashes)}]}]
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(payload).encode()))


def test_pair_ani_is_the_average_containment_estimate() -> None:
    a = frozenset(range(100))
    b = frozenset(range(20, 120))  # 80 shared of 100 each
    ani, shared = pair_ani(a, b)
    assert shared == 80
    assert ani == pytest.approx(0.8 ** (1 / KSIZE))
    assert pair_ani(a, a) == (1.0, 100)
    # Too few shared hashes, absolutely and relative to the sketches: unknown.
    big = frozenset(range(1000))
    other = frozenset(list(range(10)) + list(range(5000, 5990)))
    assert pair_ani(big, other) == (None, 10)
    # A small viral sketch keeps an estimate through the relative rule.
    assert pair_ani(frozenset(range(8)), frozenset(range(8)))[0] == 1.0


def test_load_hashes_reads_the_k31_sketch(tmp_path: Path) -> None:
    _sig(tmp_path / "x.sig.zip", [3, 1, 2])
    assert load_hashes(tmp_path / "x.sig.zip") == frozenset({1, 2, 3})
    _sig(tmp_path / "k21.sig.zip", [1], ksize=21)
    assert load_hashes(tmp_path / "k21.sig.zip") is None
    assert load_hashes(tmp_path / "absent.sig.zip") is None


def test_member_anis_with_a_missing_sketch(tmp_path: Path) -> None:
    sk = tmp_path / "sketches"
    sk.mkdir()
    _sig(sk / "rep.sig.zip", range(100))
    _sig(sk / "m1.sig.zip", range(10, 110))
    rows, missing = member_anis({"rep.fasta": ["m1.fasta", "m2.fasta"], "solo.fasta": []}, sk)
    by_member = {r.member: r for r in rows}
    assert by_member["m1.fasta"].ani == pytest.approx(0.9 ** (1 / KSIZE))
    assert by_member["m2.fasta"].ani is None and missing == 1


def test_dereplication_step_reports_member_ani(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    (wd / "derep").mkdir(parents=True)
    (wd / "sketches").mkdir()
    names = ["r.fasta", "a.fasta", "b.fasta", "c.fasta"]
    write_clusters(wd / "derep" / "clusters.tsv", {"r.fasta": names[1:]})
    write_genome_status(
        wd / "derep" / "genome_status.tsv",
        {"r.fasta": "representative", **{n: "contained" for n in names[1:]}},
    )
    _sig(wd / "sketches" / "r.sig.zip", range(200))
    _sig(wd / "sketches" / "a.sig.zip", range(200))
    _sig(wd / "sketches" / "b.sig.zip", range(10, 210))
    _sig(wd / "sketches" / "c.sig.zip", range(40, 240))
    refresh_step(wd, BY_NAME["dereplication"], formats=("png",))
    frag = wd / "report" / "dereplication"
    rep = read_fragment(frag)
    lines = (frag / "member_ani.tsv").read_text().splitlines()
    assert lines[0] == "representative\tmember\tani\tshared_hashes"
    assert lines[1].split("\t")[1:] == ["c.fasta", f"{0.8 ** (1 / KSIZE):.5f}", "160"]
    assert "member_ani" in {f.id for f in rep.figures}
    assert any(s.startswith("Estimated from the genome sketches") for s in rep.summary)
    # A new sketch makes the fragment stale.
    _sig(wd / "sketches" / "c.sig.zip", range(200))
    assert refresh_step(wd, BY_NAME["dereplication"], formats=("png",)) == "written"


def test_dereplication_step_without_sketches_has_a_note(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    (wd / "derep").mkdir(parents=True)
    write_clusters(wd / "derep" / "clusters.tsv", {"r.fasta": ["a.fasta"]})
    write_genome_status(
        wd / "derep" / "genome_status.tsv", {"r.fasta": "representative", "a.fasta": "contained"}
    )
    refresh_step(wd, BY_NAME["dereplication"], formats=("png",))
    rep = read_fragment(wd / "report" / "dereplication")
    assert any("has none (repgenr sketch -wd writes them)" in n for n in rep.notes)


def test_simple_typer_writes_reference_coverage(tmp_path: Path) -> None:
    consensuses = {"ref": "ACGTACGTAC", "g1": "ACGTANNNNN", "g2": "ACGTACGTAC"}
    fractions = _log_coverage(consensuses, "ref", logging.getLogger("t"))
    assert fractions == {"g1": 0.5, "g2": 1.0}
    path = _write_coverage(fractions, "ref", tmp_path)
    assert path.read_text().splitlines() == [
        "genome\tcovered_fraction\tis_reference",
        "ref\t1.0000\t1",
        "g1\t0.5000\t0",
        "g2\t1.0000\t0",
    ]
    assert _write_coverage({}, "ref", tmp_path / "none") is None


def test_snptype_step_reports_reference_coverage(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    (wd / "snp").mkdir(parents=True)
    (wd / "snp" / "alignment_stats.json").write_text(json.dumps({"variable_sites": 10}))
    (wd / "snp" / "reference_coverage.tsv").write_text(
        "genome\tcovered_fraction\tis_reference\nref\t1.0\t1\ng1\t0.42\t0\ng2\t0.97\t0\n"
    )
    refresh_step(wd, BY_NAME["snptype"], formats=("png",))
    frag = wd / "report" / "snptype"
    rep = read_fragment(frag)
    assert (frag / "reference_coverage.tsv").read_text().splitlines()[1:] == [
        "g1\t0.42",
        "g2\t0.97",
    ]
    assert "reference_coverage" in {f.id for f in rep.figures}
    assert any("1 genomes cover less than 50%" in n for n in rep.notes)
