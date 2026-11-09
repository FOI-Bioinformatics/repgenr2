"""metadata: several species of a genus (-ts a,b) and --no-outgroup."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

from repgenr.cli import base as cli_base
from repgenr.cli.main import app
from repgenr.cli.param_builders import join_species
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SELECTION_TSV
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams, target_species_list

_PREFIX = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales"


def _tax(family: str, genus: str, species: str) -> str:
    return f"{_PREFIX};f__{family};g__{genus};s__{genus} {species}"


_ROWS = [
    ("GCF_000001.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    ("GCF_000002.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    ("GCF_000010.1", "GCF_000010.1", _tax("Francisellaceae", "Francisella", "philomiragia")),
    ("GCF_000011.1", "GCF_000011.1", _tax("Francisellaceae", "Francisella", "philomiragia")),
    ("GCF_000030.1", "GCF_000030.1", _tax("Francisellaceae", "Francisella", "novicida")),
    ("GCF_000020.1", "GCF_000020.1", _tax("Piscirickettsiaceae", "Piscirickettsia", "salmonis")),
]


@pytest.fixture()
def table(tmp_path) -> Path:
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    lines = [
        "accession\tgtdb_genome_representative\tgtdb_taxonomy\tncbi_genbank_assembly_accession"
    ]
    lines += [f"RS_{a}\tRS_{r}\t{t}\t{a}" for a, r, t in _ROWS]
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    return path


def _params(table: Path, **overrides) -> MetadataParams:
    values = dict(
        dataset="all",
        level="species",
        release="232.0",
        version="bac120",
        target_genus="Francisella",
        target_species="tularensis,philomiragia",
        metadata_path=str(table),
    )
    values.update(overrides)
    return MetadataParams(**values)


def _rows(wd: Path) -> list[dict]:
    with open(wd / SELECTION_TSV, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def test_target_species_list_and_join() -> None:
    assert target_species_list("tularensis, philomiragia,tularensis") == [
        "tularensis",
        "philomiragia",
    ]
    assert target_species_list(None) == []
    assert join_species(["tularensis", "philomiragia,novicida"]) == (
        "tularensis,philomiragia,novicida"
    )
    assert join_species("tularensis") == "tularensis"


def test_two_species_are_selected_with_an_outgroup_outside_both(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(table)) == 4
    rows = _rows(ctx.workdir)
    ingroup = {r["accession"] for r in rows if r["is_outgroup"] == "0"}
    assert ingroup == {"GCF_000001.1", "GCF_000002.1", "GCF_000010.1", "GCF_000011.1"}
    outgroup = [r for r in rows if r["is_outgroup"] == "1"]
    # The other species of the genus, not one of the two selected.
    assert [r["accession"] for r in outgroup] == ["GCF_000030.1"]
    assert ctx.config.stages["metadata"].params["target_species"] == "tularensis,philomiragia"


def test_a_named_outgroup_inside_any_selected_species_is_refused(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="lies inside the target species"):
        metadata.run(ctx, _params(table, limit=2, outgroup_accession="GCF_000011.1"))


def test_no_species_left_for_an_outgroup_names_no_outgroup(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="--no-outgroup"):
        metadata.run(ctx, _params(table, target_species="tularensis,philomiragia,novicida"))


def test_no_outgroup_writes_the_ingroup_only(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(table, target_species="tularensis"))
    assert (ctx.workdir / "outgroup_accession.txt").exists()
    metadata.run(ctx, _params(table, target_species="tularensis", no_outgroup=True))
    rows = _rows(ctx.workdir)
    assert {r["is_outgroup"] for r in rows} == {"0"} and len(rows) == 2
    assert not (ctx.workdir / "outgroup_accession.txt").exists()
    assert [g.accession for g in ctx.manifest.all_genomes(include_outgroup=True)] == [
        "GCF_000001.1",
        "GCF_000002.1",
    ]
    record = ctx.config.stages["metadata"].params
    assert record["outgroup"] is None and record["no_outgroup"] is True


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"level": "genus"}, "select them with --level species"),
        ({"target_genus": None, "target_family": "Francisellaceae"}, "need their genus"),
        ({"outgroup_accession": "GCF_000020.1", "no_outgroup": True}, "not both"),
    ],
)
def test_usage_errors(tmp_path, table, overrides, message) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match=message):
        metadata.run(ctx, _params(table, **overrides))


def test_fingerprint_leaves_out_no_outgroup_when_off(table) -> None:
    off = cli_base._stage_fingerprint("metadata", _params(table), {}, {})
    legacy = dict(vars(_params(table)))
    legacy.pop("no_outgroup")
    import types

    assert off == cli_base._stage_fingerprint("metadata", types.SimpleNamespace(**legacy), {}, {})
    on = cli_base._stage_fingerprint("metadata", _params(table, no_outgroup=True), {}, {})
    assert on != off


def test_cli_repeated_ts_and_run_refuses_several_species_off_the_gtdb_chain(
    tmp_path, table
) -> None:
    runner = CliRunner()
    wd = tmp_path / "cli"
    argv = [
        "metadata", "-wd", str(wd), "-d", "all", "-l", "species", "-r", "232.0",
        "--gtdb-version", "bac120", "--metadata-path", str(table), "--nodownload",
        "-tg", "Francisella", "-ts", "tularensis", "-ts", "philomiragia", "--no-outgroup",
    ]  # fmt: skip
    result = runner.invoke(app, argv)
    assert result.exit_code == 0, result.output
    assert len(_rows(wd)) == 4
    viral = runner.invoke(
        app,
        ["run", "-wd", str(tmp_path / "v"), "--viral", "-tg", "Hepatovirus",
         "-ts", "a", "-ts", "b"],
    )  # fmt: skip
    assert viral.exit_code == 2
    assert "GTDB chain only" in viral.output
