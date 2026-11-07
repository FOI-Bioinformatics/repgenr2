"""End-to-end tests of the vmetadata stage and the vgenome dispatch (WP2-5).

The NCBI Virus path runs through public ``run()`` against a faked datasets CLI
that writes a real zip package (data_report.jsonl + genomic.fna); the BV-BRC
path runs against a canned group FASTA and canned Entrez taxon data. Network
and FTP never happen.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.stages import vgenome, vmetadata
from repgenr.stages.vgenome import VgenomeParams
from repgenr.stages.vmetadata import VmetadataParams
from repgenr.viral import ncbi_virus
from repgenr.viral.entrez import TAXNAMES_ORDERED


def _report_line(accession: str, organism: str, family: str, length: int) -> str:
    return json.dumps(
        {
            "accession": accession,
            "length": length,
            "completeness": "complete",
            "segment": "",
            "isolate": {"name": f"iso-{accession}"},
            "virus": {
                "organismName": organism,
                "taxId": 10508,
                "lineage": [{"name": "Viruses"}, {"name": family}, {"name": "Mastadenovirus"}],
            },
        }
    )


@pytest.fixture()
def fake_datasets(monkeypatch):
    """Fake run_tool_with_retries writing a real NCBI Virus zip package."""
    calls: list[list[str]] = []

    def fake(caps, cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        zip_path = Path(cmd[cmd.index("--filename") + 1])
        report = "\n".join(
            [
                _report_line("NC_001.1", "Human adenovirus 1", "Adenoviridae", 34000),
                _report_line("NC_002.1", "Human adenovirus 2", "Adenoviridae", 35000),
            ]
        )
        fna = ">NC_001.1 Human adenovirus 1\nACGT\n>NC_002.1 Human adenovirus 2\nACGT\n"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/data_report.jsonl", report)
            zf.writestr("ncbi_dataset/data/genomic.fna", fna)
        return 0

    def fake_taxonomy(caps, cmd, *, logger, stdout_path, **kw):
        # The NCBI Taxonomy lookup after the download (ncbi_virus.run_tool).
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        taxids = Path(cmd[cmd.index("--inputfile") + 1]).read_text(encoding="utf-8").split()
        rows = [
            {
                "query": [taxid],
                "taxonomy": {
                    "tax_id": int(taxid),
                    "classification": {
                        "family": {"id": 10508, "name": "Adenoviridae"},
                        "genus": {"id": 10509, "name": "Mastadenovirus"},
                        "species": {"id": 3, "name": "Mastadenovirus adami"},
                    },
                },
            }
            for taxid in taxids
        ]
        Path(stdout_path).write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")
        return 0

    monkeypatch.setattr(ncbi_virus, "run_tool_with_retries", fake)
    monkeypatch.setattr(ncbi_virus, "run_tool", fake_taxonomy)
    monkeypatch.setattr(vmetadata, "preflight", lambda caps: {"datasets": "16.0"}, raising=False)
    return calls


def test_ncbi_virus_end_to_end(tmp_path, fake_datasets, monkeypatch) -> None:
    import repgenr.core.plugins as plugins

    monkeypatch.setattr(plugins, "check_binaries", lambda specs, path=None: {})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = vmetadata.run(ctx, VmetadataParams(target="adenoviridae"))
    assert count == 2

    wd = ctx.workdir / "virus_download_wd"
    records = json.loads((wd / "virus_records.json").read_text(encoding="utf-8"))
    assert {r["accession"] for r in records} == {"NC_001.1", "NC_002.1"}
    assert (wd / "download.fa").exists()
    assert (ctx.workdir / "virus_metadata_base.tsv").exists()
    assert "vmetadata" in ctx.config.stages
    # filters forwarded to the datasets CLI
    flat = [tok for cmd in fake_datasets for tok in cmd]
    assert "adenoviridae" in flat
    # The species comes from the (faked) NCBI Taxonomy lookup, one call.
    assert {r["species"] for r in records} == {"Mastadenovirus-adami"}
    assert {r["species_source"] for r in records} == {"taxonomy"}
    assert sum(1 for cmd in fake_datasets if "taxonomy" in cmd) == 1
    assert ctx.config.stages["vmetadata"].params["species_source"] == {"taxonomy": 2}


def test_ncbi_virus_filters_forwarded(tmp_path, fake_datasets, monkeypatch) -> None:
    import repgenr.core.plugins as plugins

    monkeypatch.setattr(plugins, "check_binaries", lambda specs, path=None: {})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = VmetadataParams(
        target="adenoviridae", complete_only=True, host="homo sapiens", released_after="01/31/2024"
    )
    vmetadata.run(ctx, params)
    flat = [tok for cmd in fake_datasets for tok in cmd]
    for token in ("--complete-only", "--host", "homo sapiens", "--released-after", "01/31/2024"):
        assert token in flat
    # the record names the downloader and every filter, released_after included
    record = ctx.config.stages["vmetadata"]
    assert record.tool == "datasets"
    assert record.params["released_after"] == "01/31/2024"


def test_ncbi_virus_blocked_network_stops_before_datasets(
    tmp_path, fake_datasets, monkeypatch
) -> None:
    """vmetadata exited 6 after 1550 s: three datasets attempts of about 8.5 minutes."""
    import requests

    import repgenr.core.plugins as plugins
    from repgenr.core import http

    class _Unreachable:
        def get(self, url, **kw):
            raise requests.ConnectionError("Network is unreachable")

        def close(self) -> None:
            pass

    monkeypatch.setattr(plugins, "check_binaries", lambda specs, path=None: {})
    monkeypatch.setattr(http, "_probe_session", _Unreachable)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="api.ncbi.nlm.nih.gov") as info:
        vmetadata.run(ctx, VmetadataParams(target="adenoviridae"))
    assert info.value.exit_code == 3
    assert fake_datasets == []


def test_ncbi_virus_empty_package_raises(tmp_path, monkeypatch) -> None:
    def empty(caps, cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        zip_path = Path(cmd[cmd.index("--filename") + 1])
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/README.md", "empty")
        return 0

    monkeypatch.setattr(ncbi_virus, "run_tool_with_retries", empty)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="no genomes"):
        vmetadata.run(ctx, VmetadataParams(target="nosuchvirus"))


def test_missing_target_raises(tmp_path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--target"):
        vmetadata.run(ctx, VmetadataParams(target=None))


# --- BV-BRC source ------------------------------------------------------------


_BVBRC_FASTA = (
    ">acc1 Strain one, complete genome [Adenovirus A | 111.1]\nACGTACGT\n"
    ">acc2 Strain two, complete genome [Adenovirus A | 111.2]\nACGTACGTAC\n"
    # A second record of genome 111.2 (a contig of the same BV-BRC genome).
    ">acc4 Strain two contig 2, partial sequence [Adenovirus A | 111.2]\nACGTA\n"
    ">acc3 Strain three, partial sequence [Adenovirus B | 222.1]\nACGT\n"
)


def _taxdata(name: str, taxid: str) -> dict:
    levels = {level: {"name": None, "taxid": None} for level in TAXNAMES_ORDERED}
    levels["family"] = {"name": "Adenoviridae", "taxid": "10508"}
    levels["species"] = {"name": name, "taxid": taxid}
    return {"name": name, "taxdata": levels}


def test_bvbrc_end_to_end(tmp_path, monkeypatch) -> None:
    def fake_download_group(target, dest, logger):
        dest.write_text(_BVBRC_FASTA, encoding="utf-8")

    def fake_entrez(taxids, logger):
        data = {t: _taxdata(f"Species {t}", t) for t in taxids}
        return data, set(), {}

    monkeypatch.setattr(vmetadata, "_download_group", fake_download_group)
    monkeypatch.setattr(vmetadata, "get_taxon_data_from_entrez", fake_entrez)

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = VmetadataParams(target="adenoviridae", source="bvbrc")
    count = vmetadata.run(ctx, params)
    assert count == 2  # two distinct taxids seen

    wd = ctx.workdir / "virus_download_wd"
    base = (wd / "metadata_base.tsv").read_text(encoding="utf-8").splitlines()
    # only taxid 111 passes the "complete genome" filter, with 2 sequences
    assert len(base) == 2
    assert base[1].startswith("111\tAdenovirus A\t2\t8\t10")
    assert (wd / "metadata_ncbi.tsv").exists()
    assert (wd / "metadata_ncbi_taxnames_data.json").exists()
    assert (ctx.workdir / "virus_metadata_base.tsv").exists()
    assert ctx.config.stages["vmetadata"].tool == "bvbrc"

    # Sequences per NCBI taxon and rank, copied to the working directory.
    summary = (ctx.workdir / "virus_metadata_ncbi_summary.tsv").read_text(encoding="utf-8")
    rows = [line.split("\t") for line in summary.splitlines()]
    assert rows[0] == [
        "level",
        "name",
        "taxid",
        "taxids",
        "genomes",
        "sequences",
        "sequences_with_tag",
    ]
    by_name = {(r[0], r[1]): r[2:] for r in rows[1:]}
    # Genome 111.2 has two records: 3 genomes, 4 sequences, 2 with the tag.
    assert by_name[("family", "Adenoviridae")] == ["10508", "2", "3", "4", "2"]
    assert by_name[("species", "Species 111")] == ["111", "1", "2", "3", "2"]
    assert by_name[("species", "Species 222")] == ["222", "1", "1", "1", "0"]
    for row in rows[1:]:
        assert int(row[-1]) <= int(row[-2]), row
    # Ranks in order, the family before the species.
    assert [r[0] for r in rows[1:]].index("family") < [r[0] for r in rows[1:]].index("species")


def _refuse_ftp(*args, **kwargs):
    raise ConnectionRefusedError(61, "Connection refused")


@pytest.mark.parametrize(
    "params",
    [
        VmetadataParams(source="bvbrc", list_targets=True),
        VmetadataParams(target="adenoviridae", source="bvbrc"),
    ],
    ids=["list", "download"],
)
def test_bvbrc_unreachable_is_a_named_workdir_error(tmp_path, monkeypatch, params) -> None:
    """An unreachable BV-BRC FTP server raises the same named, expected error
    class as an unreachable HTTP service (core.http), not a bare OSError that
    the CLI reports as an unexpected failure (exit 1)."""
    monkeypatch.setattr(vmetadata, "_ReuseFTP_TLS", _refuse_ftp)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="BV-BRC"):
        vmetadata.run(ctx, params)


class _FakeFTP:
    """A stand-in FTPS session serving one group file."""

    def __init__(self, payload: bytes, *, fail_after: int | None = None, size: int | None = None):
        self.payload = payload
        self.fail_after = fail_after
        self.size_reported = len(payload) if size is None else size

    def cwd(self, d):
        pass

    def nlst(self):
        return ["Adenoviridae.fna"]

    def sendcmd(self, cmd):
        pass

    def size(self, remote):
        return self.size_reported

    def retrbinary(self, cmd, callback):
        if self.fail_after is None:
            callback(self.payload)
            return
        callback(self.payload[: self.fail_after])
        raise ConnectionResetError("connection dropped")


def _fake_session(ftp):
    from contextlib import contextmanager

    @contextmanager
    def session():
        yield ftp

    return session


def _logger():
    import logging

    return logging.getLogger("test")


def test_bvbrc_download_success_publishes_file(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(vmetadata, "_bvbrc_session", _fake_session(_FakeFTP(b">a\nACGT\n")))
    dest = tmp_path / "download.fa"
    vmetadata._download_group("adenoviridae", dest, _logger())
    assert dest.read_bytes() == b">a\nACGT\n"
    assert not (tmp_path / "download.fa.tmp").exists()


def test_bvbrc_interrupted_download_leaves_nothing(tmp_path, monkeypatch) -> None:
    ftp = _FakeFTP(b">a\nACGTACGT\n", fail_after=5)
    monkeypatch.setattr(vmetadata, "_bvbrc_session", _fake_session(ftp))
    dest = tmp_path / "download.fa"
    with pytest.raises(ConnectionResetError):
        vmetadata._download_group("adenoviridae", dest, _logger())
    assert not dest.exists()
    assert not (tmp_path / "download.fa.tmp").exists()


def test_bvbrc_size_mismatch_leaves_nothing(tmp_path, monkeypatch) -> None:
    ftp = _FakeFTP(b">a\nACGT\n", size=999)
    monkeypatch.setattr(vmetadata, "_bvbrc_session", _fake_session(ftp))
    dest = tmp_path / "download.fa"
    with pytest.raises(WorkdirError, match="Incomplete"):
        vmetadata._download_group("adenoviridae", dest, _logger())
    assert not dest.exists()
    assert not (tmp_path / "download.fa.tmp").exists()


def test_bvbrc_complete_existing_file_is_reused(tmp_path, monkeypatch) -> None:
    def no_download(*args, **kwargs):
        raise AssertionError("must not download")

    def fake_entrez(taxids, logger):
        return {t: _taxdata(f"Species {t}", t) for t in taxids}, set(), {}

    monkeypatch.setattr(vmetadata, "_download_group", no_download)
    monkeypatch.setattr(vmetadata, "get_taxon_data_from_entrez", fake_entrez)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(_BVBRC_FASTA, encoding="utf-8")
    (wd / "download.source").write_text("bvbrc:adenoviridae\n", encoding="utf-8")
    vmetadata.run(ctx, VmetadataParams(target="Adenoviridae", source="bvbrc"))
    assert (wd / "metadata_base.tsv").exists()


# --- vgenome dispatch ---------------------------------------------------------


def test_vgenome_missing_metadata_raises(tmp_path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="vmetadata"):
        vgenome.run(ctx, VgenomeParams())


def test_vgenome_dispatches_to_records_backend(tmp_path, monkeypatch) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(">a\nACGT\n", encoding="utf-8")
    (wd / "virus_records.json").write_text("[]", encoding="utf-8")

    import repgenr.viral.selection as selection

    called = {}
    monkeypatch.setattr(
        selection,
        "run_records",
        lambda ctx, params, dwd, fasta, records, logger: called.setdefault("n", 7) or 7,
    )
    assert vgenome.run(ctx, VgenomeParams()) == 7
    assert called["n"] == 7


def test_vgenome_dispatches_to_bvbrc_backend(tmp_path, monkeypatch) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(">a\nACGT\n", encoding="utf-8")
    (wd / "metadata_base.tsv").write_text("taxid\n", encoding="utf-8")
    (wd / "metadata_ncbi.tsv").write_text("taxid\n", encoding="utf-8")

    import repgenr.viral.bvbrc as bvbrc

    monkeypatch.setattr(
        bvbrc,
        "run_select",
        lambda ctx, params, fasta, base, ncbi, logger: 5,
    )
    assert vgenome.run(ctx, VgenomeParams()) == 5


# --- deep audit 2026-10 (entry-net) -------------------------------------------


def _bvbrc_fakes(monkeypatch) -> list[str]:
    downloads: list[str] = []

    def fake_download_group(target, dest, logger):
        downloads.append(target)
        dest.write_text(_BVBRC_FASTA, encoding="utf-8")

    def fake_entrez(taxids, logger):
        return {t: _taxdata(f"Species {t}", t) for t in taxids}, set(), {}

    monkeypatch.setattr(vmetadata, "_download_group", fake_download_group)
    monkeypatch.setattr(vmetadata, "get_taxon_data_from_entrez", fake_entrez)
    return downloads


def test_bvbrc_does_not_reuse_another_targets_group_fasta(tmp_path, monkeypatch) -> None:
    """Live: '--target picornaviridae' reused the Hepatitis E group FASTA and
    recorded it as picornaviridae."""
    downloads = _bvbrc_fakes(monkeypatch)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    vmetadata.run(ctx, VmetadataParams(target="adenoviridae", source="bvbrc"))
    vmetadata.run(ctx, VmetadataParams(target="adenoviridae", source="bvbrc"))
    vmetadata.run(ctx, VmetadataParams(target="picornaviridae", source="bvbrc"))
    assert downloads == ["adenoviridae", "picornaviridae"]
    marker = ctx.workdir / "virus_download_wd" / "download.source"
    assert marker.read_text(encoding="utf-8").strip() == "bvbrc:picornaviridae"


def test_bvbrc_after_ncbi_virus_downloads_and_drops_the_records(
    tmp_path, fake_datasets, monkeypatch
) -> None:
    """Live: switching a workdir to --source bvbrc reused NCBI Virus's
    download.fa and failed with 'list index out of range' (exit 1); the stale
    virus_records.json also kept vgenome on the NCBI Virus back-end."""
    downloads = _bvbrc_fakes(monkeypatch)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    vmetadata.run(ctx, VmetadataParams(target="adenoviridae"))
    wd = ctx.workdir / "virus_download_wd"
    assert (wd / "virus_records.json").exists()
    vmetadata.run(ctx, VmetadataParams(target="adenoviridae", source="bvbrc"))
    assert downloads == ["adenoviridae"]
    assert not (wd / "virus_records.json").exists()
    # and back: NCBI Virus drops the BV-BRC-only tables
    vmetadata.run(ctx, VmetadataParams(target="adenoviridae"))
    assert not (wd / "metadata_ncbi.tsv").exists()
    assert not (ctx.workdir / "virus_metadata_ncbi.tsv").exists()
    assert (wd / "download.source").read_text(encoding="utf-8").strip() == "ncbi_virus:adenoviridae"


def test_bvbrc_unmarked_existing_fasta_is_downloaded_again(tmp_path, monkeypatch) -> None:
    downloads = _bvbrc_fakes(monkeypatch)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(_BVBRC_FASTA, encoding="utf-8")
    vmetadata.run(ctx, VmetadataParams(target="adenoviridae", source="bvbrc"))
    assert downloads == ["adenoviridae"]


def test_bvbrc_parse_names_a_foreign_fasta(tmp_path) -> None:
    fasta = tmp_path / "download.fa"
    fasta.write_text(">NC_004297.1 Lassa virus segment L\nACGT\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="not a BV-BRC group FASTA"):
        vmetadata._parse_fasta(fasta, "complete genome", _logger())
