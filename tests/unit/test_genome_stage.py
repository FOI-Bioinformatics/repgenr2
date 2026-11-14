"""End-to-end tests of the genome stage with a faked datasets CLI (WP2-4).

``_run_cmd`` is replaced by a fake that writes what NCBI datasets would:
a dehydrated zip, rehydrated .fna files, and the outgroup zip. The stage's
organization logic (naming, pruning, manifest write-back, FASTA sanity check)
runs for real.
"""

from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.errors import WorkdirError
from repgenr.core.manifest import GenomeRecord
from repgenr.stages import genome
from repgenr.stages.genome import GenomeParams

_SELECTED = [
    GenomeRecord(
        accession="GCF_000001.1",
        source="gtdb",
        family="Francisellaceae",
        genus="Francisella",
        species="tularensis",
    ),
    GenomeRecord(
        accession="GCF_000002.1",
        source="gtdb",
        family="Francisellaceae",
        genus="Francisella",
        species="tularensis",
    ),
]
_OUTGROUP = GenomeRecord(
    accession="GCF_000010.1",
    source="gtdb",
    is_outgroup=True,
    family="Francisellaceae",
    genus="Francisella",
    species="philomiragia",
)


@pytest.fixture()
def ctx(tmp_path, monkeypatch) -> WorkdirContext:
    monkeypatch.setattr(genome, "preflight", lambda caps: {"datasets": "16.0"})
    monkeypatch.setattr(genome, "_check_disk", lambda *a, **k: None)
    c = WorkdirContext(tmp_path / "wd", create=True)
    c.manifest.upsert_many([*_SELECTED, _OUTGROUP])
    return c


def _fake_run_cmd(monkeypatch, fasta: bytes = b">seq\nACGT\n") -> list[list[str]]:
    """datasets download -> zip with a marker; rehydrate -> per-accession .fna."""
    calls: list[list[str]] = []

    def fake(cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        if cmd[:3] == ["datasets", "download", "genome"] and "--dehydrated" in cmd:
            acc_file = Path(cmd[cmd.index("--inputfile") + 1])
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("README.md", acc_file.read_text(encoding="utf-8"))
        elif cmd[:2] == ["datasets", "rehydrate"]:
            extract = Path(cmd[cmd.index("--directory") + 1])
            readme = extract / "README.md"
            for acc in readme.read_text(encoding="utf-8").splitlines():
                if acc:
                    fna = extract / "data" / acc / f"{acc}_genomic.fna"
                    fna.parent.mkdir(parents=True, exist_ok=True)
                    fna.write_bytes(fasta)
        elif cmd[:3] == ["datasets", "download", "genome"]:  # outgroup, hydrated
            acc = cmd[4]
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr(f"data/{acc}/{acc}_genomic.fna", fasta)
        else:
            raise AssertionError(f"unexpected datasets command: {cmd}")
        return 0

    monkeypatch.setattr(genome, "_run_cmd", fake)
    return calls


def test_run_downloads_and_organizes(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch)
    count = genome.run(ctx, GenomeParams())
    assert count == 2

    fastas = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert fastas == [
        "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta",
        "Francisellaceae_Francisella_tularensis_GCF_000002.1.fasta",
    ]
    out = list(ctx.outgroup_dir.iterdir())
    assert [p.name for p in out] == ["Francisellaceae_Francisella_philomiragia_GCF_000010.1.fasta"]
    # filenames recorded back into the manifest
    by_acc = {g.accession: g for g in ctx.manifest.all_genomes(include_outgroup=False)}
    assert by_acc["GCF_000001.1"].filename == fastas[0]
    assert ctx.config.stages["genome"].tool == "datasets"


def test_run_skips_present_and_prunes_stale(ctx, monkeypatch) -> None:
    calls = _fake_run_cmd(monkeypatch)
    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    present = ctx.genomes_dir / "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta"
    present.write_text(">seq\nACGT\n", encoding="utf-8")
    stale = ctx.genomes_dir / "Old_Stale_genome_GCF_999999.9.fasta"
    stale.write_text(">seq\nACGT\n", encoding="utf-8")
    keepme = ctx.genomes_dir / "notes.txt"
    keepme.write_text("user file", encoding="utf-8")

    genome.run(ctx, GenomeParams())
    assert not stale.exists(), "stale genome FASTA is pruned"
    assert keepme.exists(), "non-FASTA files are never pruned"
    downloaded = [c for c in calls if "--dehydrated" in c]
    acc_list = (ctx.workdir / "ncbi_acc_download_list.txt").read_text(encoding="utf-8")
    assert "GCF_000001.1" not in acc_list, "already-present genome is not re-downloaded"
    assert downloaded, "the missing genome still downloads"


def test_run_accession_list_only_stops_before_download(ctx, monkeypatch) -> None:
    calls = _fake_run_cmd(monkeypatch)
    assert genome.run(ctx, GenomeParams(accession_list_only=True)) == 0
    assert (ctx.workdir / "ncbi_acc_download_list.txt").exists()
    assert calls == []


def test_run_empty_manifest_raises(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(genome, "preflight", lambda caps: {})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="metadata stage"):
        genome.run(ctx, GenomeParams())


def test_non_fasta_download_raises(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch, fasta=b"<html>Service unavailable</html>")
    with pytest.raises(WorkdirError, match="not FASTA"):
        genome.run(ctx, GenomeParams())


def test_missing_accessions_logged_not_fatal(ctx, monkeypatch, caplog) -> None:
    """NCBI returning no genome for an accession warns and continues."""
    import logging

    def fake(cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        if "--dehydrated" in cmd:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("README.md", "")
        elif cmd[:2] == ["datasets", "rehydrate"]:
            extract = Path(cmd[cmd.index("--directory") + 1])
            # only the first accession comes back
            fna = extract / "data" / "GCF_000001.1" / "x.fna"
            fna.parent.mkdir(parents=True, exist_ok=True)
            fna.write_bytes(b">seq\nACGT\n")
        elif cmd[:3] == ["datasets", "download", "genome"]:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("data/GCF_000010.1/x.fna", b">seq\nACGT\n")
        return 0

    monkeypatch.setattr(genome, "_run_cmd", fake)
    with caplog.at_level(logging.WARNING):
        count = genome.run(ctx, GenomeParams())
    assert count == 2
    assert any("no genome" in r.getMessage() for r in caplog.records)


def test_stale_outgroup_files_pruned(ctx, monkeypatch) -> None:
    """Only the current outgroup remains in outgroup/ after the stage runs."""
    _fake_run_cmd(monkeypatch)
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    stale = ctx.outgroup_dir / "Old_Out_grp_GCF_777777.7.fasta"
    stale.write_text(">old\nACGT\n", encoding="utf-8")

    genome.run(ctx, GenomeParams())
    names = sorted(p.name for p in ctx.outgroup_dir.iterdir())
    assert names == ["Francisellaceae_Francisella_philomiragia_GCF_000010.1.fasta"]


def test_selection_without_outgroup_removes_an_earlier_outgroup(ctx, monkeypatch) -> None:
    """metadata --no-outgroup: the outgroup of an earlier selection is removed."""
    _fake_run_cmd(monkeypatch)
    genome.run(ctx, GenomeParams())
    assert list(ctx.outgroup_dir.iterdir())
    ctx.manifest.replace_genomes(_SELECTED)
    genome.run(ctx, GenomeParams())
    assert not [p for p in ctx.outgroup_dir.iterdir() if p.suffix == ".fasta"]


def test_html_download_not_shipped_and_recorded_missing(ctx, monkeypatch) -> None:
    """A non-FASTA body (HTML error page) is never moved into genomes/ and the
    accession lands in missing_accessions.txt instead of failing the stage."""
    from repgenr.core.contracts import MISSING_ACCESSIONS_TXT

    def fake(cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        if "--dehydrated" in cmd:
            acc_file = Path(cmd[cmd.index("--inputfile") + 1])
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("README.md", acc_file.read_text(encoding="utf-8"))
        elif cmd[:2] == ["datasets", "rehydrate"]:
            extract = Path(cmd[cmd.index("--directory") + 1])
            readme = extract / "README.md"
            for i, acc in enumerate(readme.read_text(encoding="utf-8").splitlines()):
                if not acc:
                    continue
                fna = extract / "data" / acc / f"{acc}_genomic.fna"
                fna.parent.mkdir(parents=True, exist_ok=True)
                # first accession downloads fine; second is an HTML error page
                body = b">seq\nACGT\n" if i == 0 else b"<html>error</html>"
                fna.write_bytes(body)
        elif cmd[:3] == ["datasets", "download", "genome"]:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("data/GCF_000010.1/x.fna", b">seq\nACGT\n")
        return 0

    monkeypatch.setattr(genome, "_run_cmd", fake)
    count = genome.run(ctx, GenomeParams())
    assert count == 2
    names = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert names == ["Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta"]
    missing = (ctx.workdir / MISSING_ACCESSIONS_TXT).read_text(encoding="utf-8")
    assert "GCF_000002.1" in missing


def test_present_html_file_is_redownloaded(ctx, monkeypatch) -> None:
    """An HTML file left at a genome path by a crashed run is re-downloaded,
    not accepted as already present."""
    calls = _fake_run_cmd(monkeypatch)
    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    bad = ctx.genomes_dir / "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta"
    bad.write_text("<html>Service unavailable</html>", encoding="utf-8")

    genome.run(ctx, GenomeParams())
    assert bad.read_text(encoding="utf-8").startswith(">"), "bad file replaced"
    acc_list = (ctx.workdir / "ncbi_acc_download_list.txt").read_text(encoding="utf-8")
    assert "GCF_000001.1" in acc_list
    assert any("--dehydrated" in c for c in calls)


def test_accession_list_ends_with_newline(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch)
    genome.run(ctx, GenomeParams(accession_list_only=True))
    text = (ctx.workdir / "ncbi_acc_download_list.txt").read_text(encoding="utf-8")
    assert text.endswith("\n") and not text.endswith("\n\n")


# --- deep audit 2026-10 (entry-net) -------------------------------------------

_NO_MATCH = "Error: There are no genome assemblies that match your query."


def test_batch_of_only_unserved_accessions_is_recorded_missing(ctx, monkeypatch) -> None:
    """Live: a rerun whose remaining accession NCBI no longer serves failed with
    exit 6 after three attempts; it is recorded in missing_accessions.txt."""
    from repgenr.core.contracts import MISSING_ACCESSIONS_TXT

    calls = _fake_run_cmd(monkeypatch)
    inner = genome._run_cmd
    attempts = []

    def fake(cmd, *, logger, permanent=None, **kw):
        if "--dehydrated" in [str(c) for c in cmd]:
            attempts.append(cmd)
            exc = genome.ToolExecutionError([str(c) for c in cmd], 1, output=_NO_MATCH)
            assert permanent is not None and permanent(exc), "no retries for a no-match"
            raise exc
        return inner(cmd, logger=logger, **kw)

    monkeypatch.setattr(genome, "_run_cmd", fake)
    assert genome.run(ctx, GenomeParams()) == 2
    missing = (ctx.workdir / MISSING_ACCESSIONS_TXT).read_text(encoding="utf-8").split()
    assert missing == ["GCF_000001.1", "GCF_000002.1"]
    # One request for the selection and one for its GenBank pairs.
    assert len(attempts) == 2
    assert ctx.config.stages["genome"].completed
    assert calls  # the outgroup was still fetched


def _fake_suppressed(monkeypatch, suppressed: set[str]) -> list[list[str]]:
    """As ``_fake_run_cmd``, but NCBI serves nothing for ``suppressed``."""
    calls = _fake_run_cmd(monkeypatch)
    inner = genome._run_cmd

    def fake(cmd, *, logger, **kw):
        text = [str(c) for c in cmd]
        if text[:2] == ["datasets", "rehydrate"]:
            readme = Path(text[text.index("--directory") + 1]) / "README.md"
            kept = [
                a for a in readme.read_text(encoding="utf-8").splitlines() if a not in suppressed
            ]
            readme.write_text("\n".join(kept), encoding="utf-8")
        elif "--dehydrated" not in text and text[4] in suppressed:  # the outgroup
            calls.append(text)
            raise genome.ToolExecutionError(text, 1, output=_NO_MATCH)
        return inner(cmd, logger=logger, **kw)

    monkeypatch.setattr(genome, "_run_cmd", fake)
    return calls


def _substituted(ctx) -> str:
    from repgenr.core.contracts import SUBSTITUTED_ACCESSIONS_TSV

    return (ctx.workdir / SUBSTITUTED_ACCESSIONS_TSV).read_text(encoding="utf-8")


def test_unserved_refseq_genome_is_taken_from_genbank_and_recorded(ctx, monkeypatch) -> None:
    from repgenr.core.contracts import MISSING_ACCESSIONS_TXT

    _fake_suppressed(monkeypatch, {"GCF_000002.1"})
    assert genome.run(ctx, GenomeParams()) == 2
    # The file keeps the name of the selected accession.
    assert (ctx.genomes_dir / "Francisellaceae_Francisella_tularensis_GCF_000002.1.fasta").exists()
    assert (ctx.workdir / MISSING_ACCESSIONS_TXT).read_text(encoding="utf-8") == ""
    assert _substituted(ctx) == "selected\tdownloaded\nGCF_000002.1\tGCA_000002.1\n"
    assert ctx.config.stages["genome"].params["substituted"] == 1
    # A second run downloads nothing and keeps the record.
    genome.run(ctx, GenomeParams())
    assert _substituted(ctx) == "selected\tdownloaded\nGCF_000002.1\tGCA_000002.1\n"


def test_no_genbank_fallback_records_the_accession_as_missing(ctx, monkeypatch) -> None:
    from repgenr.core.contracts import MISSING_ACCESSIONS_TXT

    _fake_suppressed(monkeypatch, {"GCF_000002.1"})
    genome.run(ctx, GenomeParams(genbank_fallback=False))
    assert (ctx.workdir / MISSING_ACCESSIONS_TXT).read_text(encoding="utf-8") == "GCF_000002.1\n"
    assert _substituted(ctx) == "selected\tdownloaded\n"


def test_unserved_refseq_outgroup_is_taken_from_genbank(ctx, monkeypatch) -> None:
    _fake_suppressed(monkeypatch, {"GCF_000010.1"})
    genome.run(ctx, GenomeParams())
    names = [p.name for p in ctx.outgroup_dir.iterdir()]
    assert names == ["Francisellaceae_Francisella_philomiragia_GCF_000010.1.fasta"]
    assert _substituted(ctx) == "selected\tdownloaded\nGCF_000010.1\tGCA_000010.1\n"
    # The outgroup is on disk now; its record is kept by a second run.
    genome.run(ctx, GenomeParams())
    assert "GCA_000010.1" in _substituted(ctx)


def test_outgroup_without_any_assembly_stops_the_stage(ctx, monkeypatch) -> None:
    from repgenr.core.errors import WorkdirError

    _fake_suppressed(monkeypatch, {"GCF_000010.1", "GCA_000010.1"})
    with pytest.raises(WorkdirError, match="nor for its GenBank pair GCA_000010.1"):
        genome.run(ctx, GenomeParams())


def test_other_datasets_failures_still_fail_the_stage(ctx, monkeypatch) -> None:
    def fake(cmd, *, logger, **kw):
        raise genome.ToolExecutionError([str(c) for c in cmd], 1, output="connection reset")

    monkeypatch.setattr(genome, "_run_cmd", fake)
    with pytest.raises(genome.ToolExecutionError):
        genome.run(ctx, GenomeParams())


def test_rehydrated_genome_failing_its_md5_is_discarded(ctx, monkeypatch) -> None:
    """datasets rehydrate does not check md5sum.txt; the stage does."""
    import hashlib

    from repgenr.core.contracts import MISSING_ACCESSIONS_TXT

    good = b">seq\nACGT\n"
    inner_calls = _fake_run_cmd(monkeypatch, fasta=good)
    inner = genome._run_cmd

    def fake(cmd, *, logger, **kw):
        rc = inner(cmd, logger=logger, **kw)
        if [str(c) for c in cmd][:2] == ["datasets", "rehydrate"]:
            extract = Path(str(cmd[cmd.index("--directory") + 1]))
            digest = hashlib.md5(good).hexdigest()
            (extract / "md5sum.txt").write_text(
                f"{digest}  data/GCF_000001.1/GCF_000001.1_genomic.fna\n"
                f"{'0' * 32}  data/GCF_000002.1/GCF_000002.1_genomic.fna\n",
                encoding="utf-8",
            )
        return rc

    monkeypatch.setattr(genome, "_run_cmd", fake)
    genome.run(ctx, GenomeParams())
    names = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert names == ["Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta"]
    missing = (ctx.workdir / MISSING_ACCESSIONS_TXT).read_text(encoding="utf-8")
    assert "GCF_000002.1" in missing
    assert inner_calls


def test_unserved_outgroup_is_a_named_workdir_error(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch)
    inner = genome._run_cmd

    def fake(cmd, *, logger, **kw):
        if "--dehydrated" not in [str(c) for c in cmd] and str(cmd[1]) == "download":
            raise genome.ToolExecutionError([str(c) for c in cmd], 1, output=_NO_MATCH)
        return inner(cmd, logger=logger, **kw)

    monkeypatch.setattr(genome, "_run_cmd", fake)
    with pytest.raises(WorkdirError, match="no assembly for the outgroup GCF_000010.1"):
        genome.run(ctx, GenomeParams())


def test_outgroup_package_without_fasta_is_an_error_not_silence(ctx, monkeypatch) -> None:
    _fake_run_cmd(monkeypatch)
    inner = genome._run_cmd

    def fake(cmd, *, logger, **kw):
        cmd_s = [str(c) for c in cmd]
        if "--dehydrated" not in cmd_s and cmd_s[1] == "download":
            with zipfile.ZipFile(cmd_s[cmd_s.index("--filename") + 1], "w") as zf:
                zf.writestr("README.md", "no data")
            return 0
        return inner(cmd, logger=logger, **kw)

    monkeypatch.setattr(genome, "_run_cmd", fake)
    with pytest.raises(WorkdirError, match="holds no genome FASTA"):
        genome.run(ctx, GenomeParams())


def test_present_outgroup_is_not_downloaded_again(ctx, monkeypatch) -> None:
    calls = _fake_run_cmd(monkeypatch)
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    present = ctx.outgroup_dir / "Francisellaceae_Francisella_philomiragia_GCF_000010.1.fasta"
    present.write_text(">og\nACGT\n", encoding="utf-8")
    genome.run(ctx, GenomeParams())
    assert not any("GCF_000010.1" in c for c in calls)
    assert present.read_text(encoding="utf-8") == ">og\nACGT\n"


class _Unreachable:
    def get(self, url, **kw):
        import requests

        raise requests.ConnectTimeout("timed out")

    def close(self) -> None:
        pass


def test_blocked_network_stops_before_datasets(ctx, monkeypatch) -> None:
    """datasets on a blocked network made three attempts of about 8.5 minutes each."""
    from repgenr.core import http

    calls = _fake_run_cmd(monkeypatch)
    monkeypatch.setattr(http, "_probe_session", _Unreachable)
    with pytest.raises(WorkdirError, match="api.ncbi.nlm.nih.gov") as info:
        genome.run(ctx, GenomeParams())
    assert info.value.exit_code == 3
    assert calls == []


def test_blocked_network_stops_before_the_outgroup_download(ctx, monkeypatch) -> None:
    from repgenr.core import http

    calls = _fake_run_cmd(monkeypatch)
    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    for g in _SELECTED:
        (ctx.genomes_dir / genome._output_name(g)).write_text(">s\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(http, "_probe_session", _Unreachable)
    with pytest.raises(WorkdirError, match="api.ncbi.nlm.nih.gov"):
        genome.run(ctx, GenomeParams())
    assert calls == []


def test_nothing_to_download_needs_no_network(ctx, monkeypatch) -> None:
    from repgenr.core import http

    calls = _fake_run_cmd(monkeypatch)
    ctx.genomes_dir.mkdir(parents=True, exist_ok=True)
    for g in _SELECTED:
        (ctx.genomes_dir / genome._output_name(g)).write_text(">s\nACGT\n", encoding="utf-8")
    ctx.outgroup_dir.mkdir(parents=True, exist_ok=True)
    (ctx.outgroup_dir / genome._output_name(_OUTGROUP)).write_text(">o\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(http, "_probe_session", _Unreachable)
    assert genome.run(ctx, GenomeParams()) == 2
    assert calls == []


def test_rehydrate_hosts_from_fetch_txt_are_probed_before_rehydrate(ctx, monkeypatch) -> None:
    """Rehydrate fetches from the URLs in fetch.txt; their host is checked first."""
    import requests

    from repgenr.core import http

    calls = _fake_run_cmd(monkeypatch)
    real_fake = genome._run_cmd

    def with_fetch_txt(cmd, **kw):
        result = real_fake(cmd, **kw)
        cmd = [str(c) for c in cmd]
        if "--dehydrated" in cmd:
            zip_path = Path(cmd[cmd.index("--filename") + 1])
            with zipfile.ZipFile(zip_path, "a") as zf:
                zf.writestr(
                    "ncbi_dataset/fetch.txt",
                    "https://files.example.org/a/x.fna\t10\tdata/x.fna\n"
                    "https://files.example.org/a/y.fna\t10\tdata/y.fna\n",
                )
        return result

    monkeypatch.setattr(genome, "_run_cmd", with_fetch_txt)
    probed: list[str] = []

    class _Session:
        def get(self, url, **kw):
            probed.append(url)
            if "files.example.org" in url:
                raise requests.ConnectionError("refused")

        def close(self) -> None:
            pass

    monkeypatch.setattr(http, "_probe_session", _Session)
    with pytest.raises(WorkdirError, match="files.example.org"):
        genome.run(ctx, GenomeParams())
    assert probed == [http.NCBI_DATASETS_URL, "https://files.example.org/"]
    assert not [c for c in calls if c[:2] == ["datasets", "rehydrate"]]
