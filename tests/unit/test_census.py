"""repgenr census: the three modes, the JSON and TSV output and the exit codes.

Every request is faked: the GTDB API (``gtdb_api._api_get``), ENA
(``ena.resolve_taxon``/``ena.search_runs``), the Entrez lineage lookup, and
the ``datasets`` runs behind the NCBI Virus report. The tables are small and
synthetic.
"""

from __future__ import annotations

import gzip
import json
import logging
import re
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core import ena
from repgenr.core import http as core_http
from repgenr.core.config import Config
from repgenr.core.contracts import (
    SELECTION_TSV,
    SelectionRow,
    genome_filename,
    write_clusters,
    write_selection,
)
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.manifest import Manifest, record_from_selection
from repgenr.sources import gtdb_api
from repgenr.stages import census, reads
from repgenr.viral import ncbi_virus

runner = CliRunner()
LOG = logging.getLogger("test.census")
# Rich renders usage errors in colour and wraps them at the terminal width on
# CI; plain, wide output and flattened whitespace keep flag names matchable.
_PLAIN = {"TERM": "dumb", "NO_COLOR": "1", "COLUMNS": "200"}
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_BOX = re.compile(r"[\s\u2500-\u257f]+")


def _plain_error(args: list[str]) -> tuple[int, str]:
    """Exit code and the flattened output of a census invocation."""
    result = runner.invoke(app, ["census", *args], env=_PLAIN)
    text = _ANSI.sub("", result.output + str(result.exception or ""))
    return result.exit_code, _BOX.sub(" ", text)


# --- GTDB fixtures -------------------------------------------------------------

_PREFIX = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Testales"


def _tax(family: str, genus: str, species: str) -> str:
    return f"{_PREFIX};f__{family};g__{genus};s__{genus} {species}"


# accession, representative, taxonomy. Testaceae: 2 genera, 3 species,
# 5 genomes, 3 representatives; one genome of another family.
GTDB_ROWS = [
    ("GCF_000001.1", "GCF_000001.1", _tax("Testaceae", "Alpha", "one")),
    ("GCF_000002.1", "GCF_000001.1", _tax("Testaceae", "Alpha", "one")),
    ("GCF_000003.1", "GCF_000003.1", _tax("Testaceae", "Alpha", "two")),
    ("GCF_000004.1", "GCF_000004.1", _tax("Testaceae", "Beta", "three")),
    ("GCF_000005.1", "GCF_000004.1", _tax("Testaceae", "Beta", "three")),
    ("GCF_000009.1", "GCF_000009.1", _tax("Otheraceae", "Gamma", "four")),
]


def _write_table(path: Path) -> Path:
    header = "accession\tgtdb_genome_representative\tgtdb_taxonomy\tncbi_genbank_assembly_accession"
    lines = [header] + [f"RS_{a}\tRS_{r}\t{t}\t{a}" for a, r, t in GTDB_ROWS]
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    return path


@pytest.fixture
def gtdb_table(tmp_path) -> Path:
    return _write_table(tmp_path / "bac120_metadata_r232.tsv.gz")


def _api_rows() -> list[dict]:
    rows = []
    for acc, rep, tax in GTDB_ROWS:
        parts = dict(p.split("__", 1) for p in tax.split(";"))
        rows.append(
            {
                "gid": acc,
                "gtdbFamily": f"f__{parts['f']}",
                "gtdbGenus": f"g__{parts['g']}",
                "gtdbSpecies": f"s__{parts['s']}",
                "gtdbIsRep": acc == rep,
            }
        )
    return rows


@pytest.fixture
def fake_api(monkeypatch) -> list[str]:
    """The GTDB API answering genomes-detail for f__Testaceae and g__Alpha/g__Beta."""
    calls: list[str] = []

    def api_get(path: str, params: dict | None = None) -> dict:
        calls.append(path)
        taxon = path.split("/")[2].replace("%20", " ").replace("%5F", "_")
        rank, name = taxon.split("__", 1)
        key = {"f": "gtdbFamily", "g": "gtdbGenus"}[rank]
        return {"rows": [r for r in _api_rows() if r[key] == f"{rank}__{name}"]}

    monkeypatch.setattr(gtdb_api, "_api_get", api_get)
    return calls


# --- ENA fixtures --------------------------------------------------------------

# Runs across two species and three platforms; ERR3 and ERR4 share a biosample.
ENA_RECORDS = [
    {"run_accession": "ERR1", "sample_accession": "SAMEA1", "tax_id": "101",
     "instrument_platform": "ILLUMINA"},
    {"run_accession": "ERR2", "sample_accession": "SAMEA2", "tax_id": "101",
     "instrument_platform": "OXFORD_NANOPORE"},
    {"run_accession": "ERR3", "sample_accession": "SAMEA3", "tax_id": "102",
     "instrument_platform": "ILLUMINA"},
    {"run_accession": "ERR4", "sample_accession": "SAMEA3", "tax_id": "102",
     "instrument_platform": "PACBIO_SMRT"},
    {"run_accession": "ERR4", "sample_accession": "SAMEA3", "tax_id": "102",
     "instrument_platform": "PACBIO_SMRT"},
]  # fmt: skip


@pytest.fixture
def fake_ena(monkeypatch) -> dict[str, int]:
    seen = {"entrez": 0}

    def resolve(name: str, **kw) -> ena.TaxonHit:
        assert kw == {"division": "PRO"}
        return ena.TaxonHit("1000", name, "family")

    def search(query: str, **kw) -> list[dict]:
        assert query == ena.taxon_query("1000")
        return [dict(r) for r in ENA_RECORDS]

    def entrez(taxids, logger, **kw):
        seen["entrez"] += 1
        assert sorted(taxids) == ["101", "102"], "one lookup per distinct taxid"
        species = {"101": "Alpha one", "102": "Beta three"}
        data = {
            t: {
                "taxdata": {
                    "family": {"name": "Testaceae"},
                    "genus": {"name": species[t].split()[0]},
                    "species": {"name": species[t]},
                }
            }
            for t in taxids
        }
        return data, set(), {}

    monkeypatch.setattr(ena, "resolve_taxon", resolve)
    monkeypatch.setattr(ena, "search_runs", search)
    monkeypatch.setattr(ena, "count_runs", lambda query: 2)
    monkeypatch.setattr(reads, "get_taxon_data_from_entrez", entrez)
    return seen


# --- NCBI Virus fixtures -------------------------------------------------------


def _virus(acc, genus, species, *, segment="", isolate="", complete=True, length=1000) -> dict:
    lineage = ["Viruses", "Testviridae", genus, f"{genus} {species}"]
    return {
        "accession": acc,
        "completeness": "COMPLETE" if complete else "PARTIAL",
        "length": length,
        "segment": segment,
        "isolate": {"name": isolate} if isolate else {},
        "virus": {
            "lineage": [{"name": n, "tax_id": i} for i, n in enumerate(lineage, 1)],
            "organism_name": f"{species} virus",
            "tax_id": {"alpha": 501, "beta": 502, "gamma": 503}[species],
        },
    }


# Segvirus is segmented (L and S): alpha holds two complete isolates and one
# lone L segment, beta one isolate. Monovirus is not segmented.
VIRUS_REPORT = [
    _virus("SA1", "Segvirus", "alpha", segment="L", isolate="iso1"),
    _virus("SA2", "Segvirus", "alpha", segment="S", isolate="iso1"),
    _virus("SA3", "Segvirus", "alpha", segment="L", isolate="iso2"),
    _virus("SA4", "Segvirus", "alpha", segment="S RNA", isolate="iso2", complete=False),
    _virus("SA5", "Segvirus", "alpha", segment="L", isolate="iso3"),
    _virus("SB1", "Segvirus", "beta", segment="L", isolate="isoX"),
    _virus("SB2", "Segvirus", "beta", segment="S", isolate="isoX"),
    _virus("MG1", "Monovirus", "gamma", isolate="m1"),
    _virus("MG2", "Monovirus", "gamma", isolate="m2", complete=False),
    _virus("MG3", "Monovirus", "gamma"),
]


@pytest.fixture
def fake_datasets(monkeypatch) -> list[list[str]]:
    """``datasets summary virus`` prints VIRUS_REPORT; the taxonomy lookup
    finds nothing, so the species come from the report lineage."""
    commands: list[list[str]] = []

    def run_tool(caps, cmd, *, logger, log_prefix=None, stdout_path=None, timeout=None, **kw):
        commands.append([str(c) for c in cmd])
        text = ""
        if cmd[:3] == ["datasets", "summary", "virus"]:
            text = "\n".join(json.dumps(r) for r in VIRUS_REPORT) + "\n"
        Path(stdout_path).write_text(text, encoding="utf-8")
        return 0

    monkeypatch.setattr(ncbi_virus, "run_tool", run_tool)
    monkeypatch.setattr(ncbi_virus, "require_reachable", lambda url, what: None)
    return commands


# --- helpers ------------------------------------------------------------------


def _json(args: list[str]) -> dict:
    result = runner.invoke(app, ["census", *args, "--json"])
    assert result.exit_code == 0, result.output
    return json.loads(result.stdout)


def _rows(payload: dict) -> dict[str, dict]:
    return {r["name"]: r for r in payload["rows"]}


def _snapshot(root: Path) -> dict[str, tuple[int, int]]:
    return {
        str(p.relative_to(root)): (p.stat().st_mtime_ns, p.stat().st_size)
        for p in root.rglob("*")
        if p.is_file()
    }


def _record(wd: Path, *stages: tuple[str, dict], inputs: dict | None = None) -> None:
    cfg = Config()
    for name, params in stages:
        cfg.record_stage(
            name, params=params, completed="2026-10-09T00:00:00+00:00", inputs=inputs or {}
        )
    cfg.save(wd)


# --- mode 1: GTDB --------------------------------------------------------------


def test_family_from_the_api_counts_genera_species_and_representatives(fake_api) -> None:
    payload = _json(["-tf", "testaceae"])
    assert payload["mode"] == "taxon" and payload["source"] == "gtdb-api"
    assert payload["taxon"] == "Testaceae" and payload["rank"] == "family"
    assert payload["totals"] == {"genera": 2, "species": 3, "genomes": 5, "representatives": 3}
    rows = _rows(payload)
    assert rows["Alpha"] == {"name": "Alpha", "species": 2, "genomes": 3, "representatives": 2}
    assert rows["Beta"] == {"name": "Beta", "species": 1, "genomes": 2, "representatives": 1}
    assert fake_api == ["/taxon/f__Testaceae/genomes-detail"]


def test_genus_from_the_api_has_one_row_per_species(fake_api) -> None:
    payload = _json(["-tg", "Alpha"])
    assert payload["rank"] == "genus" and payload["row_rank"] == "species"
    assert [r["name"] for r in payload["rows"]] == ["Alpha one", "Alpha two"]
    assert payload["rows"][0] == {"name": "Alpha one", "genomes": 2, "representatives": 1}
    assert "species" not in payload["rows"][0]


def test_console_header_and_table(fake_api) -> None:
    result = runner.invoke(app, ["census", "-tf", "Testaceae"], env=_PLAIN)
    assert result.exit_code == 0, result.output
    lines = _ANSI.sub("", result.stdout).splitlines()
    assert lines[0] == "Testaceae: 2 genera, 3 species, 5 genomes (3 representatives)"
    assert lines[1].split() == ["name", "species", "genomes", "representatives"]
    assert lines[2].split() == ["Alpha", "2", "3", "2"]


def test_table_source_reads_metadata_path(gtdb_table) -> None:
    payload = _json(["-tf", "Testaceae", "--source", "table", "--metadata-path", str(gtdb_table)])
    assert payload["source"] == "gtdb-table"
    assert payload["totals"] == {"genera": 2, "species": 3, "genomes": 5, "representatives": 3}


def test_table_source_reuses_the_cached_release_table(tmp_path, monkeypatch) -> None:
    cache = tmp_path / "cache"
    (cache / "gtdb").mkdir(parents=True)
    _write_table(cache / "gtdb" / "bac120_metadata_r232.tsv.gz")
    (cache / "gtdb" / "bac120_metadata_r232.release").write_text("232.0\n", encoding="utf-8")
    monkeypatch.setenv(census.CACHE_ENV, str(cache))

    def no_download(*a, **k):
        raise AssertionError("the cached table must be reused")

    monkeypatch.setattr(core_http, "download", no_download)
    payload = _json(["-tg", "Beta", "--source", "tsv", "-r", "232.0"])
    assert payload["totals"]["genomes"] == 2
    assert census.cache_dir() == cache / "gtdb"


def test_table_source_needs_a_release() -> None:
    code, text = _plain_error(["-tg", "Beta", "--source", "table"])
    assert code == 2
    assert "--release" in text


def test_runs_add_ena_columns_grouped_by_species(fake_api, fake_ena) -> None:
    payload = _json(["-tf", "Testaceae", "--runs"])
    assert payload["totals"]["runs"] == 4
    assert payload["totals"]["biosamples"] == 3
    rows = _rows(payload)
    alpha, beta = rows["Alpha"], rows["Beta"]
    assert (alpha["runs"], alpha["biosamples"], alpha["illumina"], alpha["ont"]) == (2, 2, 1, 1)
    assert (beta["runs"], beta["biosamples"], beta["illumina"], beta["pacbio"]) == (2, 1, 1, 1)
    # The species column counts GTDB species only.
    assert alpha["species"] == 2
    assert fake_ena["entrez"] == 1


def test_runs_per_species_in_a_genus(fake_api, fake_ena, monkeypatch) -> None:
    payload = _json(["-tg", "Beta", "--runs"])
    row = _rows(payload)["Beta three"]
    assert (row["genomes"], row["runs"], row["biosamples"], row["pacbio"]) == (2, 2, 1, 1)
    # Alpha runs sit under the same ENA taxid fake: a row of runs alone,
    # marked NCBI and listed last.
    assert payload["rows"][-1]["name"] == "Alpha one (NCBI)"
    assert payload["rows"][-1]["genomes"] == 0 and payload["rows"][-1]["runs"] == 2
    assert payload["totals"]["genomes"] == 2 and payload["totals"]["species"] == 1


def test_network_failure_exits_3(monkeypatch) -> None:
    def down(path, params=None):
        raise WorkdirError("GTDB API request failed: connection refused")

    monkeypatch.setattr(gtdb_api, "_api_get", down)
    result = runner.invoke(app, ["census", "-tg", "Alpha"])
    assert result.exit_code == 3


# --- mode 1: NCBI Virus ---------------------------------------------------------


def test_viral_family_counts_sequences_isolates_and_segments(fake_datasets) -> None:
    payload = _json(["--viral", "--target", "testviridae", "--complete-only", "--host", "human"])
    summary_cmd = fake_datasets[0]
    assert summary_cmd[:6] == ["datasets", "summary", "virus", "genome", "taxon", "testviridae"]
    assert "--complete-only" in summary_cmd and summary_cmd[-2:] == ["--host", "human"]
    assert not any(c == "download" for c in summary_cmd)
    assert payload["rank"] == "family" and payload["source"] == "ncbi_virus"
    assert payload["totals"] == {
        "genera": 2,
        "species": 3,
        "sequences": 10,
        "complete": 8,
        "isolates": 7,
    }
    rows = _rows(payload)
    assert rows["Segvirus"] == {
        "name": "Segvirus",
        "species": 2,
        "sequences": 7,
        "complete": 6,
        "isolates": 4,
        "segmented": "yes",
    }
    assert rows["Monovirus"]["isolates"] == 3 and rows["Monovirus"]["segmented"] == "no"


def test_viral_genus_target_counts_species(fake_datasets) -> None:
    payload = _json(["--viral", "--target", "segvirus"])
    assert payload["rank"] == "genus" and payload["taxon"] == "Segvirus"
    rows = _rows(payload)
    # Species are shown as NCBI writes them, not as filename tokens.
    assert set(rows) == {"Segvirus alpha", "Segvirus beta"}
    assert rows["Segvirus alpha"]["isolates"] == 3 and rows["Segvirus beta"]["isolates"] == 1


def test_viral_tg_narrows_a_family_report(fake_datasets) -> None:
    payload = _json(["--viral", "--target", "testviridae", "-tg", "Monovirus"])
    assert payload["rank"] == "genus"
    assert payload["totals"]["sequences"] == 3


@pytest.mark.parametrize(
    "args",
    [
        ["--viral", "--target", "x", "--source", "bvbrc"],
        ["--viral", "--target", "x", "--runs"],
        ["--viral"],
        [],
        ["-tf", "Testaceae", "--source", "ncbi_virus"],
        ["--target", "picornaviridae"],
    ],
)
def test_usage_errors_exit_2(args) -> None:
    result = runner.invoke(app, ["census", *args])
    assert result.exit_code == 2, result.output


def test_bvbrc_without_a_workdir_names_vmetadata() -> None:
    with pytest.raises(UserInputError, match="vmetadata"):
        census.validate(census.CensusParams(viral=True, target="x", source="bvbrc"))


def test_unknown_source_is_rejected_naming_the_flag() -> None:
    code, text = _plain_error(["-tf", "X", "--source", "bogus"])
    assert code == 2 and "--source" in text


# --- mode 2 ---------------------------------------------------------------------


def _bacterial_selection(wd: Path, rows: list[tuple]) -> list[SelectionRow]:
    out = []
    for acc, family, genus, species, outgroup, rep in rows:
        out.append(
            SelectionRow(
                acc,
                family,
                genus,
                species,
                outgroup,
                genome_filename(family, genus, species, acc),
                source_representative=rep,
            )
        )
    write_selection(wd / SELECTION_TSV, out)
    return out


def test_candidates_from_the_table_named_in_the_metadata_record(tmp_path, gtdb_table) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    _record(
        wd,
        ("metadata", {"source": "tsv", "release": "232.0", "version": "bac120",
                      "level": "genus", "target_genus": "Alpha"}),
        inputs={str(gtdb_table): "digest"},
    )  # fmt: skip
    _bacterial_selection(wd, [("GCF_000001.1", "Testaceae", "Alpha", "one", False, True)])
    before = _snapshot(wd)
    payload = _json(["-wd", str(wd)])
    assert payload["mode"] == "candidates" and payload["source"] == "gtdb-table"
    assert payload["taxon"] == "Alpha" and payload["rank"] == "genus"
    assert payload["totals"]["genomes"] == 3
    assert _snapshot(wd) == before


def test_candidates_from_the_kept_api_answer(tmp_path) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    _record(wd, ("metadata", {"source": "api", "level": "family", "target_family": "Testaceae"}))
    rows = [r for r in _api_rows() if r["gtdbFamily"] == "f__Testaceae"]
    gtdb_api.write_api_genomes(wd / gtdb_api.GTDB_API_GENOMES, rows)
    payload = _json(["-wd", str(wd)])
    assert payload["source"] == "gtdb-api" and payload["rank"] == "family"
    assert payload["totals"] == {"genera": 2, "species": 3, "genomes": 5, "representatives": 3}
    narrowed = _json(["-wd", str(wd), "-tg", "beta"])
    assert narrowed["rank"] == "genus" and narrowed["totals"]["genomes"] == 2


def test_candidates_fall_back_to_the_selection(tmp_path) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    _record(wd, ("metadata", {"source": "api", "level": "genus", "target_genus": "Alpha"}))
    _bacterial_selection(
        wd,
        [
            ("GCF_000001.1", "Testaceae", "Alpha", "one", False, True),
            ("GCF_000004.1", "Testaceae", "Beta", "three", True, True),
        ],
    )
    payload = _json(["-wd", str(wd)])
    assert payload["source"] == "selection" and payload["totals"]["genomes"] == 1
    assert any("selection.tsv" in n for n in payload["notes"])


def test_api_answer_round_trips() -> None:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / gtdb_api.GTDB_API_GENOMES
        gtdb_api.write_api_genomes(path, _api_rows())
        back = gtdb_api.read_api_genomes(path)
    assert back[0] == {
        "accession": "GCF_000001.1",
        "tax": {
            "family": "Testaceae",
            "genus": "Alpha",
            "species": "one",
            "family_name": "Testaceae",
            "genus_name": "Alpha",
            "species_name": "Alpha one",
        },
        "is_rep": True,
    }
    assert sum(r["is_rep"] for r in back) == 4


def test_ncbi_virus_candidates(tmp_path) -> None:
    wd = tmp_path / "wd"
    download_wd = wd / "virus_download_wd"
    download_wd.mkdir(parents=True)
    _record(wd, ("vmetadata", {"source": "ncbi_virus", "target": "testviridae"}))
    records = ncbi_virus.parse_report([json.dumps(r) for r in VIRUS_REPORT])
    ncbi_virus.write_records(download_wd / "virus_records.json", records)
    payload = _json(["-wd", str(wd)])
    assert payload["mode"] == "candidates" and payload["rank"] == "family"
    assert payload["totals"]["isolates"] == 7


def _bvbrc_workdir(tmp_path: Path) -> Path:
    """A BV-BRC vmetadata workdir: family Testviridae, genera Segvirus (two
    species) and Monovirus (one species); record ids are taxid.n."""
    wd = tmp_path / "bv"
    download_wd = wd / "virus_download_wd"
    download_wd.mkdir(parents=True)
    _record(wd, ("vmetadata", {"source": "bvbrc", "target": "testviridae"}))
    alpha, beta, gamma = ["11.1", "11.2", "11.3"], ["12.1"], ["13.1", "13.2"]
    taxnames = {
        "Testviridae": {"taxid": "1", "level": "family", "datasets": alpha + beta + gamma},
        "Segvirus": {"taxid": "2", "level": "genus", "datasets": alpha + beta},
        "Monovirus": {"taxid": "3", "level": "genus", "datasets": gamma},
        "Segvirus alpha": {"taxid": "11", "level": "species", "datasets": alpha},
        "Segvirus beta": {"taxid": "12", "level": "species", "datasets": beta},
        "Monovirus gamma": {"taxid": "13", "level": "species", "datasets": gamma},
        "Viruses": {"taxid": "0", "level": "superkingdom", "datasets": alpha + beta + gamma},
    }
    (download_wd / census.BVBRC_TAXNAMES).write_text(json.dumps(taxnames), encoding="utf-8")
    header = "taxid\tname\tnum\tseq_min\tseq_max\tseq_med\tseq_mean\tdescription\n"
    base = "11\ta\t3\t900\t1100\t1000\t1000\td\n12\tb\t1\t3000\t3000\t3000\t3000\td\n"
    (download_wd / "metadata_base.tsv").write_text(header + base, encoding="utf-8")
    return wd


def test_bvbrc_candidates_from_the_taxnames_sets(tmp_path) -> None:
    wd = _bvbrc_workdir(tmp_path)
    payload = _json(["-wd", str(wd)])
    assert payload["source"] == "bvbrc" and payload["rank"] == "family"
    assert payload["totals"] == {
        "genera": 2,
        "species": 3,
        "sequences": 6,
        "complete": None,
        "isolates": None,
    }
    rows = _rows(payload)
    assert rows["Segvirus"]["species"] == 2 and rows["Segvirus"]["sequences"] == 4
    assert rows["Segvirus"]["median_length"] == 2000
    assert rows["Monovirus"]["median_length"] is None
    assert rows["Segvirus"]["isolates"] is None
    genus = _json(["-wd", str(wd), "-tg", "segvirus"])
    assert _rows(genus)["Segvirus alpha"]["sequences"] == 3


# --- mode 3 ---------------------------------------------------------------------


def _selection_workdir(tmp_path: Path) -> Path:
    """genome/assemble/ingest genomes from three sources, an outgroup, one
    selected genome missing from genomes/, a kept API answer and clusters."""
    wd = tmp_path / "sel"
    wd.mkdir()
    _record(
        wd,
        ("metadata", {"source": "api", "level": "family", "target_family": "Testaceae"}),
        ("genome", {}),
    )
    gtdb_api.write_api_genomes(
        wd / gtdb_api.GTDB_API_GENOMES,
        [r for r in _api_rows() if r["gtdbFamily"] == "f__Testaceae"],
    )
    rows = _bacterial_selection(
        wd,
        [
            ("GCF_000001.1", "Testaceae", "Alpha", "one", False, True),
            ("GCF_000002.1", "Testaceae", "Alpha", "one", False, False),
            ("SRR1", "Testaceae", "Alpha", "one", False, False),
            ("local1", "Testaceae", "Beta", "three", False, False),
            ("GCF_000004.1", "Testaceae", "Beta", "three", False, True),
            ("GCF_000009.1", "Otheraceae", "Gamma", "four", True, True),
        ],
    )
    sources = {"SRR1": "sra", "local1": "local"}
    manifest = Manifest.open(wd)
    manifest.replace_genomes(
        [record_from_selection(r, sources.get(r.accession, "gtdb")) for r in rows]
    )
    manifest.close()
    genomes = wd / "genomes"
    genomes.mkdir()
    for r in rows[:-2]:  # GCF_000004.1 was not downloaded
        (genomes / r.filename).write_text(">x\nACGT\n", encoding="utf-8")
    (wd / "outgroup").mkdir()
    (wd / "derep").mkdir()
    clusters = {
        rows[0].filename: [rows[0].filename, rows[1].filename],
        rows[3].filename: [rows[3].filename],
        rows[2].filename: [rows[2].filename],
    }
    write_clusters(wd / "derep" / "clusters.tsv", clusters)
    return wd


def test_selection_counts_by_source_with_candidates_and_clusters(tmp_path) -> None:
    wd = _selection_workdir(tmp_path)
    before = _snapshot(wd)
    payload = _json(["-wd", str(wd)])
    assert _snapshot(wd) == before, "census must not write into the workdir"
    assert payload["mode"] == "selection" and payload["rank"] == "family"
    totals = payload["totals"]
    assert totals["genomes"] == 4 and totals["outgroup"] == "GCF_000009.1"
    assert totals["sources"] == {"gtdb": 2, "sra": 1, "local": 1}
    assert totals["candidates"] == 5 and totals["clusters"] == 3
    rows = _rows(payload)
    assert rows["Alpha"] == {
        "name": "Alpha",
        "species": 1,
        "candidates": 3,
        "genomes": 3,
        "gtdb": 2,
        "sra": 1,
        "local": 0,
        "representatives": 1,
        "clusters": 2,
    }
    assert rows["Beta"]["genomes"] == 1 and rows["Beta"]["candidates"] == 2
    assert any("absent from genomes/" in n for n in payload["notes"])


def test_selection_console_and_tsv(tmp_path) -> None:
    wd = _selection_workdir(tmp_path)
    out = tmp_path / "census.tsv"
    result = runner.invoke(app, ["census", "-wd", str(wd), "--tsv", str(out)], env=_PLAIN)
    assert result.exit_code == 0, result.output
    header = _ANSI.sub("", result.stdout).splitlines()[0]
    assert header.startswith("Testaceae: 2 genera, 2 species, 4 genomes (1 representative)")
    assert "outgroup GCF_000009.1 excluded" in header
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "name",
        "species",
        "candidates",
        "genomes",
        "gtdb",
        "sra",
        "local",
        "representatives",
        "clusters",
    ]
    assert lines[1].split("\t") == ["Alpha", "1", "3", "3", "2", "1", "0", "1", "2"]
    assert len(lines) == 3


def test_selection_narrowed_to_a_genus(tmp_path) -> None:
    wd = _selection_workdir(tmp_path)
    payload = _json(["-wd", str(wd), "-tg", "Alpha"])
    assert payload["rank"] == "genus" and payload["taxon"] == "Alpha"
    # A candidate species with no selected genome shows what the selection left out.
    assert payload["rows"] == [
        {
            "name": "Alpha one",
            "candidates": 2,
            "genomes": 3,
            "gtdb": 2,
            "sra": 1,
            "representatives": 1,
            "clusters": 2,
        },
        {
            "name": "Alpha two",
            "candidates": 1,
            "genomes": 0,
            "gtdb": 0,
            "sra": 0,
            "representatives": 0,
            "clusters": 0,
        },
    ]


def test_selection_without_candidates_or_clusters(tmp_path) -> None:
    wd = tmp_path / "ing"
    wd.mkdir()
    _record(wd, ("ingest", {}))
    rows = _bacterial_selection(wd, [("a1", "Testaceae", "Alpha", "one", False, False)])
    (wd / "genomes").mkdir()
    (wd / "genomes" / rows[0].filename).write_text(">a\nAC\n", encoding="utf-8")
    payload = _json(["-wd", str(wd)])
    assert payload["rows"][0] == {
        "name": "Alpha one",
        "genomes": 1,
        "unknown": 1,
        "representatives": 0,
    }
    assert "candidates" not in payload["totals"] and "clusters" not in payload["totals"]


def test_assemble_workdir_takes_reads_as_candidates(tmp_path) -> None:
    from repgenr.core.contracts import READS_TSV, ReadRow, write_reads

    wd = tmp_path / "asm"
    wd.mkdir()
    _record(wd, ("reads", {}), ("assemble", {}))
    rows = _bacterial_selection(wd, [("SRR1", "Testaceae", "Alpha", "one", False, False)])
    (wd / "genomes").mkdir()
    (wd / "genomes" / rows[0].filename).write_text(">a\nAC\n", encoding="utf-8")
    write_reads(
        wd / READS_TSV,
        [
            ReadRow(run_accession=f"SRR{i}", biosample=f"S{i}", bioproject="P", organism="o",
                    taxid="1", platform="ILLUMINA", instrument_model="m", layout="PAIRED",
                    bases=1, read_count=1, fastq_urls=(), fastq_md5=(), fastq_bytes=(),
                    family="Testaceae", genus="Alpha", species="one")
            for i in (1, 2)
        ],
    )  # fmt: skip
    payload = _json(["-wd", str(wd)])
    assert payload["rows"][0]["candidates"] == 2 and payload["rows"][0]["genomes"] == 1


def test_workdir_flags_that_select_a_taxon_query_exit_2(tmp_path) -> None:
    wd = _selection_workdir(tmp_path)
    code, text = _plain_error(["-wd", str(wd), "--source", "api"])
    assert code == 2 and "--source" in text


def test_missing_workdir_exits_3(tmp_path) -> None:
    result = runner.invoke(app, ["census", "-wd", str(tmp_path / "absent")])
    assert result.exit_code == 3


def test_workdir_without_an_entry_record_exits_2(tmp_path) -> None:
    result = runner.invoke(app, ["census", "-wd", str(tmp_path)])
    assert result.exit_code == 2


def test_tsv_that_cannot_be_written_exits_2(tmp_path, fake_api) -> None:
    out = tmp_path / "a-directory"
    out.mkdir()
    result = runner.invoke(app, ["census", "-tf", "Testaceae", "--tsv", str(out)])
    assert result.exit_code == 2


def test_json_prints_only_the_object(fake_api) -> None:
    result = runner.invoke(app, ["census", "-tg", "Beta", "--json"])
    payload = json.loads(result.stdout)
    assert set(payload) == {
        "schema",
        "repgenr",
        "taxon",
        "mode",
        "source",
        "rank",
        "row_rank",
        "totals",
        "rows",
        "notes",
        "runs_taxon",
    }
    assert payload["schema"] == census.SCHEMA


def test_summary_writes_nothing_outside_a_temporary_directory(fake_datasets, tmp_path) -> None:
    import os

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        records = ncbi_virus.summary("testviridae", logger=LOG, released_after="01/31/2024")
    finally:
        os.chdir(cwd)
    assert len(records) == 10 and list(tmp_path.iterdir()) == []
    assert fake_datasets[0][-2:] == ["--released-after", "01/31/2024"]


# --- review fixes: GTDB spellings, NCBI-only rows, stale candidates ----------

_SUFFIX_ROWS = [
    ("GCF_100001.1", "GCF_100001.1", _tax("Bacillaceae", "Bacillus_A", "thuringiensis_S")),
    ("GCF_100002.1", "GCF_100001.1", _tax("Bacillaceae", "Bacillus_A", "thuringiensis_S")),
    ("GCF_100003.1", "GCF_100003.1", _tax("Bacillaceae", "Bacillus_A", "cereus")),
    ("GCF_100004.1", "GCF_100004.1", _tax("Bacillaceae", "Bacillus", "subtilis")),
]


def _suffix_api(monkeypatch) -> None:
    rows = []
    for acc, rep, tax in _SUFFIX_ROWS:
        parts = dict(p.split("__", 1) for p in tax.split(";"))
        rows.append(
            {
                "gid": acc,
                "gtdbFamily": f"f__{parts['f']}",
                "gtdbGenus": f"g__{parts['g']}",
                "gtdbSpecies": f"s__{parts['s']}",
                "gtdbIsRep": acc == rep,
            }
        )

    def api_get(path: str, params: dict | None = None) -> dict:
        taxon = path.split("/")[2].replace("%20", " ")
        rank, name = taxon.split("__", 1)
        key = {"f": "gtdbFamily", "g": "gtdbGenus"}[rank]
        return {"rows": [r for r in rows if r[key] == f"{rank}__{name}"]}

    monkeypatch.setattr(gtdb_api, "_api_get", api_get)


def test_gtdb_names_keep_their_suffixes_from_the_api(monkeypatch) -> None:
    _suffix_api(monkeypatch)
    genus = _json(["-tg", "Bacillus_A"])
    assert genus["taxon"] == "Bacillus_A"
    assert [r["name"] for r in genus["rows"]] == [
        "Bacillus_A thuringiensis_S",
        "Bacillus_A cereus",
    ]
    family = _json(["-tf", "bacillaceae"])
    assert family["taxon"] == "Bacillaceae"
    assert {r["name"] for r in family["rows"]} == {"Bacillus_A", "Bacillus"}
    # Two genera that differ only by the suffix stay two rows.
    assert family["totals"]["genera"] == 2


def test_gtdb_names_keep_their_suffixes_from_the_table(tmp_path) -> None:
    table = tmp_path / "suffix.tsv.gz"
    header = "accession\tgtdb_genome_representative\tgtdb_taxonomy\tncbi_genbank_assembly_accession"
    lines = [header] + [f"RS_{a}\tRS_{r}\t{t}\t{a}" for a, r, t in _SUFFIX_ROWS]
    with gzip.open(table, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    payload = _json(["-tg", "Bacillus_A", "--source", "table", "--metadata-path", str(table)])
    assert _rows(payload)["Bacillus_A thuringiensis_S"]["genomes"] == 2


def test_runs_join_gtdb_rows_and_ncbi_only_rows_are_marked(monkeypatch, fake_api) -> None:
    """A run of a GTDB genus joins its row; runs of a genus GTDB lacks and of
    a taxid without a genus are marked NCBI, listed last and left out of the
    GTDB totals."""
    runs = [
        {"run_accession": "ERR1", "sample_accession": "S1", "tax_id": "101",
         "instrument_platform": "ILLUMINA"},
        {"run_accession": "ERR2", "sample_accession": "S2", "tax_id": "201",
         "instrument_platform": "ILLUMINA"},
        {"run_accession": "ERR3", "sample_accession": "S3", "tax_id": "301",
         "instrument_platform": "OXFORD_NANOPORE"},
    ]  # fmt: skip
    lineages = {
        "101": ("Alpha", "Alpha one"),
        "201": ("Allotest", "Allotest five"),
        "301": ("", "Testaceae bacterium"),
    }

    def entrez(taxids, logger, **kw):
        data = {
            t: {
                "taxdata": {
                    "family": {"name": "Testaceae"},
                    "genus": {"name": lineages[t][0] or None},
                    "species": {"name": lineages[t][1]},
                }
            }
            for t in taxids
        }
        return data, set(), {}

    monkeypatch.setattr(
        ena, "resolve_taxon", lambda name, **kw: ena.TaxonHit("1000", name, "family")
    )
    monkeypatch.setattr(ena, "search_runs", lambda query, **kw: [dict(r) for r in runs])
    monkeypatch.setattr(ena, "count_runs", lambda query: 2)
    monkeypatch.setattr(reads, "get_taxon_data_from_entrez", entrez)
    payload = _json(["-tf", "Testaceae", "--runs"])
    names = [r["name"] for r in payload["rows"]]
    assert names[:2] == ["Alpha", "Beta"]
    assert set(names[2:]) == {"Allotest (NCBI)", "(no genus, NCBI)"}
    assert _rows(payload)["Alpha"]["runs"] == 1
    assert payload["totals"]["genera"] == 2 and payload["totals"]["genomes"] == 5
    assert payload["totals"]["runs"] == 3
    assert any("may be a GTDB taxon under another name" in n for n in payload["notes"])


def test_selection_with_outdated_virus_records_counts_without_candidates(tmp_path) -> None:
    """A virus_records.json from before the taxonomy lookup is refused by
    read_records (WorkdirError); the census still counts the selection."""
    wd = tmp_path / "vg"
    download_wd = wd / "virus_download_wd"
    download_wd.mkdir(parents=True)
    _record(
        wd,
        ("vmetadata", {"source": "ncbi_virus", "target": "testviridae"}),
        ("vgenome", {}),
    )
    old = [{"accession": "MG1", "taxid": "1", "organism": "o", "family": "Testviridae",
            "genus": "Monovirus", "species": "gamma-virus", "length": 1, "completeness": "",
            "segment": "", "isolate": ""}]  # fmt: skip
    (download_wd / "virus_records.json").write_text(json.dumps(old), encoding="utf-8")
    row = SelectionRow(
        "MG1",
        "Testviridae",
        "Monovirus",
        "Monovirus-gamma",
        False,
        genome_filename("Testviridae", "Monovirus", "Monovirus-gamma", "MG1"),
    )
    write_selection(wd / SELECTION_TSV, [row])
    (wd / "genomes").mkdir()
    (wd / "genomes" / row.filename).write_text(">MG1\nAC\n", encoding="utf-8")
    payload = _json(["-wd", str(wd)])
    assert payload["mode"] == "selection" and payload["totals"]["genomes"] == 1
    assert "candidates" not in payload["totals"]
    assert "candidates" not in payload["rows"][0]
    assert any(n.startswith("No candidates column") for n in payload["notes"])


# --- review nits ---------------------------------------------------------------


def test_a_genus_outside_the_named_family_exits_2(fake_api) -> None:
    code, text = _plain_error(["-tg", "Alpha", "-tf", "Otheraceae"])
    assert code == 2 and "Otheraceae" in text


def test_runs_of_a_suffixed_gtdb_genus_use_the_ncbi_genus(monkeypatch) -> None:
    _suffix_api(monkeypatch)
    asked: list[str] = []

    def resolve(name: str, **kw) -> ena.TaxonHit:
        asked.append(name)
        return ena.TaxonHit("1386", name, "genus")

    runs = [
        {"run_accession": "ERR1", "sample_accession": "S1", "tax_id": "1",
         "instrument_platform": "ILLUMINA"},
        {"run_accession": "ERR2", "sample_accession": "S2", "tax_id": "2",
         "instrument_platform": "BGISEQ"},
        {"run_accession": "ERR3", "sample_accession": "S3", "tax_id": "2",
         "instrument_platform": "ION_TORRENT"},
    ]  # fmt: skip
    species = {"1": "Bacillus cereus", "2": "Bacillus sp. LA11-2445"}

    def entrez(taxids, logger, **kw):
        data = {
            t: {
                "taxdata": {
                    "family": {"name": "Bacillaceae"},
                    "genus": {"name": "Bacillus"},
                    "species": {"name": species[t]},
                }
            }
            for t in taxids
        }
        return data, set(), {}

    searched: list[dict] = []

    def search(query, **kw):
        searched.append(kw)
        return [dict(r) for r in runs]

    monkeypatch.setattr(ena, "resolve_taxon", resolve)
    monkeypatch.setattr(ena, "search_runs", search)
    monkeypatch.setattr(ena, "count_runs", lambda query: 3)
    monkeypatch.setattr(reads, "get_taxon_data_from_entrez", entrez)
    payload = _json(["-tg", "Bacillus_A", "--runs"])
    assert asked == ["Bacillus"]
    assert any("NCBI genus Bacillus" in n for n in payload["notes"])
    # The census asks the portal for the four fields it reads, not RUN_FIELDS.
    assert searched == [{"fields": census.CENSUS_RUN_FIELDS}]
    # The runs are those of the NCBI genus, which the header names.
    assert payload["runs_taxon"] == "NCBI genus Bacillus"
    head = runner.invoke(app, ["census", "-tg", "Bacillus_A", "--runs"])
    assert head.exit_code == 0, head.output
    assert "3 runs of NCBI genus Bacillus, 3 biosamples" in head.stdout
    rows = _rows(payload)
    # A Bacillus cereus run joins the GTDB row Bacillus_A cereus by epithet.
    assert rows["Bacillus_A cereus"]["runs"] == 1 and rows["Bacillus_A cereus"]["genomes"] == 1
    # NCBI species keep their spaces; runs on other platforms count as 'other'.
    odd = rows["Bacillus sp. LA11-2445 (NCBI)"]
    assert (odd["runs"], odd["illumina"], odd["other"]) == (2, 0, 2)
    for row in payload["rows"]:
        platforms = row["illumina"] + row["ont"] + row["pacbio"] + row["other"]
        assert platforms == row["runs"]


def test_isolate_grouping_logs_nothing_on_the_shared_loggers(fake_datasets, caplog) -> None:
    with caplog.at_level(logging.INFO):
        _json(["--viral", "--target", "segvirus"])
    assert logging.getLogger("repgenr.census.segments").level == logging.NOTSET


def test_resolve_taxon_prefers_the_named_division(monkeypatch) -> None:
    hits = [
        {"taxId": "55087", "scientificName": "Bacillus", "rank": "genus", "division": "INV"},
        {"taxId": "1386", "scientificName": "Bacillus", "rank": "genus", "division": "PRO"},
    ]
    monkeypatch.setattr(ena.http, "get_json", lambda url, **kw: hits)
    assert ena.resolve_taxon("Bacillus", division="PRO").taxid == "1386"
    with pytest.raises(UserInputError, match="matches 2 taxa"):
        ena.resolve_taxon("Bacillus")


def test_census_fields_are_those_the_census_reads() -> None:
    """Each census field is one of RUN_FIELDS, and ena_runs reads no other."""
    assert set(census.CENSUS_RUN_FIELDS) <= set(ena.RUN_FIELDS)
    record = dict.fromkeys(census.CENSUS_RUN_FIELDS, "x")
    record["run_accession"] = "ERR1"
    (row,) = ena.to_read_rows([record])
    assert (row.run_accession, row.biosample, row.taxid, row.platform) == ("ERR1", "x", "x", "x")


def test_census_header_without_a_runs_taxon_is_unchanged() -> None:
    c = census.Census(
        taxon="Staphylococcus",
        mode="taxon",
        source="api",
        rank="genus",
        columns=[],
        rows=[],
        totals={"genera": 1, "species": 2, "genomes": 3, "representatives": 2,
                "runs": 5, "biosamples": 4},
    )  # fmt: skip
    assert c.header().endswith("3 genomes (2 representatives), 5 runs, 4 biosamples")
    assert c.to_json()["runs_taxon"] == ""


def test_a_large_census_search_is_announced_and_proceeds(monkeypatch, caplog) -> None:
    monkeypatch.setattr(ena, "count_runs", lambda query: 878821)
    monkeypatch.setattr(ena, "resolve_taxon", lambda name, **kw: ena.TaxonHit("590", name, "genus"))
    monkeypatch.setattr(ena, "search_runs", lambda query, **kw: [])
    monkeypatch.setattr(reads, "ncbi_taxon_names", lambda taxids, logger: {})
    with caplog.at_level(logging.INFO):
        assert census.ena_runs("Salmonella", logging.getLogger("repgenr.test")) == []
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "Salmonella has 878821 ENA runs" in warnings[0] and "GB of memory" in warnings[0]
