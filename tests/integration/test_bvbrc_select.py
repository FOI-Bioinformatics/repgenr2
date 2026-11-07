"""Offline tests for the legacy BV-BRC viral path (no network, no mashtree).

Covers the FASTA header parsing (vmetadata._parse_fasta), the length-range
derivation, and an end-to-end bvbrc.run_select with --no-outgroup so the genome
selection/writing runs without the mashtree outgroup step.
"""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.stages.vgenome import VgenomeParams
from repgenr.stages.vmetadata import _parse_fasta
from repgenr.viral import bvbrc
from repgenr.viral.entrez import TAXNAMES_ORDERED

_LOG = logging.getLogger("test")

# BV-BRC header form: ">{name} {desc} | {bvbrc_id}] [{species} | {tag}".
_FASTA = (
    f">acc1 prot | 10535.1] [Human mastadenovirus C | complete genome\n{'A' * 300}\n"
    f">acc2 prot | 10535.2] [Human mastadenovirus C | complete genome\n{'A' * 310}\n"
    f">acc3 prot | 99999.1] [Other virus | complete genome\n{'A' * 50}\n"
)


def _write_metadata(download_wd: Path) -> tuple[Path, Path]:
    base_tsv = download_wd / "metadata_base.tsv"
    base_tsv.write_text(
        "taxid\tname\tnum\tseq_min\tseq_max\tseq_med\tseq_mean\tdescription\n"
        "10535\tHuman mastadenovirus C\t2\t300\t310\t305\t305\tdesc\n"
        "99999\tOther virus\t1\t50\t50\t50\t50\tdesc\n"
    )
    n = len(TAXNAMES_ORDERED)
    gi, si = TAXNAMES_ORDERED.index("genus"), TAXNAMES_ORDERED.index("species")

    def _row(taxid: str, genus: str, genus_taxid: str, species: str) -> str:
        names = [""] * n
        taxids = [""] * n
        names[gi], taxids[gi] = genus, genus_taxid
        names[si], taxids[si] = species, "0"
        return "\t".join([taxid, species, "2", *names, *taxids])

    ncbi_tsv = download_wd / "metadata_ncbi.tsv"
    header = [
        "taxid",
        "name",
        "num_with_tag",
        *TAXNAMES_ORDERED,
        *[f"{x}_taxid" for x in TAXNAMES_ORDERED],
    ]
    ncbi_tsv.write_text(
        "\t".join(header)
        + "\n"
        + _row("10535", "Mastadenovirus", "10509", "Human mastadenovirus C")
        + "\n"
        + _row("99999", "Othervirus", "88888", "Other virus")
        + "\n"
    )
    return base_tsv, ncbi_tsv


def test_parse_fasta_groups_by_taxid(tmp_path: Path) -> None:
    fasta = tmp_path / "download.fa"
    fasta.write_text(_FASTA)
    base, all_taxids, taxid_bvbrc, records = _parse_fasta(fasta, "complete genome", _LOG)
    assert all_taxids == {"10535", "99999"}
    assert base["10535"]["num"] == 2
    assert base["10535"]["min"] == 300 and base["10535"]["max"] == 310
    assert taxid_bvbrc["10535"] == {"10535.1", "10535.2"}
    assert records == {"10535": 2, "99999": 1}


def test_determine_length_range_from_base() -> None:
    base = {"10535": {"seq_med": 300, "seq_mean": 300}}
    selected = {"10535": {"10535.1": 300}}
    params = VgenomeParams(target_genus="mastadenovirus", length_deviation=10)
    lo, hi = bvbrc._determine_length_range(selected, base, params, _LOG)
    assert (lo, hi) == (270, 330)


def test_run_select_writes_matching_genomes(workdir: Path) -> None:
    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)

    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    n = bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)

    assert n == 2  # the two mastadenovirus sequences; the 'Other virus' is excluded
    written = {p.name for p in ctx.genomes_dir.iterdir()}
    assert written == {"acc1.fasta", "acc2.fasta"}


def test_run_select_parses_fasta_once(workdir: Path, monkeypatch) -> None:
    # The legacy path used to re-parse the whole download four times; it should
    # now do a single metadata pass (sequences come from the lazy index).
    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)

    real_parse = bvbrc.SeqIO.parse
    calls = {"n": 0}

    def counting_parse(*a, **k):
        calls["n"] += 1
        return real_parse(*a, **k)

    monkeypatch.setattr(bvbrc.SeqIO, "parse", counting_parse)
    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    assert calls["n"] == 1  # one metadata scan; sequences via SeqIO.index


def test_run_select_publishes_selection_tsv(workdir: Path) -> None:
    """The BV-BRC path writes the same selection.tsv hand-off as NCBI Virus (audit)."""
    from repgenr.core.contracts import SELECTION_TSV, read_selection

    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)
    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    n = bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    rows = read_selection(ctx.workdir / SELECTION_TSV)
    assert len(rows) == n == 2
    assert {r.filename for r in rows} == {p.name for p in (workdir / "genomes").iterdir()}
    assert all(r.genus for r in rows) and not any(r.is_outgroup for r in rows)


def test_run_select_populates_manifest(workdir: Path) -> None:
    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)
    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    records = ctx.manifest.all_genomes()
    assert {g.filename for g in records} == {p.name for p in (workdir / "genomes").iterdir()}
    assert all(g.genus and g.source == "bvbrc" for g in records)


def test_run_select_pinned_outgroup(workdir: Path) -> None:
    """--outgroup-accession names a record id; the mashtree search is skipped."""
    from repgenr.core.contracts import SELECTION_TSV, read_selection

    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)
    params = VgenomeParams(
        target_genus="mastadenovirus", length_range="250-350", outgroup_accession="99999.1"
    )
    bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    # The BV-BRC id resolves to the record; the outgroup is staged under its record id.
    assert (workdir / "outgroup_accession.txt").read_text().strip() == "acc3"
    assert [p.name for p in ctx.outgroup_dir.iterdir()] == ["acc3.fasta"]
    rows = read_selection(ctx.workdir / SELECTION_TSV)
    assert {r.accession for r in rows if r.is_outgroup} == {"acc3"}


def test_run_select_crash_mid_write_leaves_the_previous_genomes_intact(workdir: Path, monkeypatch):
    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)
    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    before = {p.name: p.read_text() for p in ctx.genomes_dir.iterdir()}
    assert len(before) == 2

    real = bvbrc._write_record
    calls = {"n": 0}

    def flaky(record, path):
        calls["n"] += 1
        if calls["n"] == 2:
            raise RuntimeError("disk full")
        return real(record, path)

    monkeypatch.setattr(bvbrc, "_write_record", flaky)
    import pytest

    with pytest.raises(RuntimeError):
        bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    assert {p.name: p.read_text() for p in ctx.genomes_dir.iterdir()} == before
    assert not [p for p in workdir.iterdir() if "staging" in p.name]


def test_run_select_without_outgroup_removes_an_earlier_outgroup(workdir: Path) -> None:
    """A later --no-outgroup run leaves no outgroup behind for phylo to root on."""
    ctx = WorkdirContext(workdir, create=True)
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True)
    fasta = download_wd / "download.fa"
    fasta.write_text(_FASTA)
    base_tsv, ncbi_tsv = _write_metadata(download_wd)
    pinned = VgenomeParams(
        target_genus="mastadenovirus", length_range="250-350", outgroup_accession="99999.1"
    )
    bvbrc.run_select(ctx, pinned, fasta, base_tsv, ncbi_tsv, _LOG)
    stale = ctx.outgroup_dir / "older_run.fasta"
    stale.write_text(">x\nACGT\n")
    bvbrc.run_select(ctx, pinned, fasta, base_tsv, ncbi_tsv, _LOG)
    assert [p.name for p in ctx.outgroup_dir.iterdir()] == ["acc3.fasta"]
    params = VgenomeParams(target_genus="mastadenovirus", no_outgroup=True, length_range="250-350")
    bvbrc.run_select(ctx, params, fasta, base_tsv, ncbi_tsv, _LOG)
    assert not (workdir / "outgroup_accession.txt").exists()
    assert not list(ctx.outgroup_dir.iterdir())
