"""metadata --source tsv --taxonomy ncbi: the GTDB table selected by NCBI names."""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SELECTION_TSV
from repgenr.core.errors import UserInputError
from repgenr.sources import gtdb_table
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams

LOG = logging.getLogger("test-ncbi-taxonomy")
_PREFIX = "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales"


def _tax(family: str, genus: str, species: str) -> str:
    return f"{_PREFIX};f__{family};g__{genus};s__{species}"


# accession, GTDB representative, GTDB taxonomy, NCBI taxonomy
_ROWS = [
    (
        "GCF_000001.1",
        "GCF_000001.1",
        _tax("Bacillaceae_G", "Bacillus_A", "Bacillus_A cereus"),
        _tax("Bacillaceae", "Bacillus", "Bacillus cereus"),
    ),
    # One NCBI species, another GTDB species.
    (
        "GCF_000002.1",
        "GCF_000002.1",
        _tax("Bacillaceae_G", "Bacillus_A", "Bacillus_A thuringiensis"),
        _tax("Bacillaceae", "Bacillus", "Bacillus cereus"),
    ),
    (
        "GCF_000003.1",
        "GCF_000003.1",
        _tax("Bacillaceae", "Bacillus", "Bacillus subtilis"),
        _tax("Bacillaceae", "Bacillus", "Bacillus subtilis"),
    ),
    # An NCBI name with characters that are unsafe in a file name.
    (
        "GCF_000004.1",
        "GCF_000003.1",
        _tax("Bacillaceae", "Bacillus", "Bacillus sp000004"),
        _tax("Bacillaceae", "Bacillus", "Bacillus sp. X-1 [strain]"),
    ),
    # Another genus of the NCBI family: the outgroup of a genus selection.
    (
        "GCF_000005.1",
        "GCF_000005.1",
        _tax("Bacillaceae_H", "Priestia", "Priestia megaterium"),
        _tax("Bacillaceae", "Priestia", "Priestia megaterium"),
    ),
    # No NCBI names below the order.
    (
        "GCF_000006.1",
        "GCF_000006.1",
        _tax("Bacillaceae_G", "Bacillus_A", "Bacillus_A sp000006"),
        f"{_PREFIX};f__;g__;s__",
    ),
]


def _write_table(path: Path, *, with_ncbi: bool = True) -> Path:
    header = ["accession", "gtdb_genome_representative", "gtdb_taxonomy"]
    header += ["ncbi_taxonomy"] if with_ncbi else []
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\t".join(header) + "\n")
        for acc, rep, gtdb, ncbi in _ROWS:
            fields = [f"RS_{acc}", f"RS_{rep}", gtdb] + ([ncbi] if with_ncbi else [])
            fo.write("\t".join(fields) + "\n")
    return path


@pytest.fixture()
def table(tmp_path) -> Path:
    return _write_table(tmp_path / "bac120_metadata_r232.tsv.gz")


def _params(table: Path, **overrides) -> MetadataParams:
    defaults = dict(
        dataset="all",
        level="genus",
        release="232.0",
        version="bac120",
        target_genus="Bacillus",
        metadata_path=str(table),
        taxonomy="ncbi",
    )
    defaults.update(overrides)
    return MetadataParams(**defaults)


def _selection(workdir: Path) -> dict[str, dict]:
    with open(workdir / SELECTION_TSV, encoding="utf-8", newline="") as fo:
        return {row["accession"]: row for row in csv.DictReader(fo, delimiter="\t")}


def _ingroup(rows: dict[str, dict]) -> set[str]:
    return {acc for acc, row in rows.items() if row["is_outgroup"] == "0"}


def test_genus_is_selected_by_its_ncbi_name(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(table))
    rows = _selection(ctx.workdir)
    # The NCBI genus Bacillus holds the GTDB genera Bacillus and Bacillus_A.
    assert _ingroup(rows) == {"GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCF_000004.1"}
    assert rows["GCF_000001.1"]["filename"] == "Bacillaceae_Bacillus_cereus_GCF_000001.1.fasta"
    # Two GTDB species, one NCBI species.
    assert rows["GCF_000002.1"]["species"] == "cereus"
    assert rows["GCF_000004.1"]["species"] == "sp-X-1-strain"
    # The GTDB representative flag is kept.
    flag = next(k for k in rows["GCF_000004.1"] if k.endswith("_representative"))
    assert rows["GCF_000003.1"][flag] == "1" and rows["GCF_000004.1"][flag] == "0"
    record = ctx.config.stages["metadata"].params
    assert record["taxonomy"] == "ncbi" and record["release"] == "232.0"


def test_the_same_target_by_gtdb_names_is_another_set(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(table, taxonomy="gtdb", no_outgroup=True))
    assert _ingroup(_selection(ctx.workdir)) == {"GCF_000003.1", "GCF_000004.1"}
    # The default is not recorded, so earlier records stay as they were.
    assert "taxonomy" not in ctx.config.stages["metadata"].params


def test_species_is_selected_by_its_ncbi_name(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(table, level="species", target_species="cereus"))
    rows = _selection(ctx.workdir)
    assert _ingroup(rows) == {"GCF_000001.1", "GCF_000002.1"}
    # The outgroup: a GTDB representative of the NCBI genus outside the species.
    (outgroup,) = [acc for acc, row in rows.items() if row["is_outgroup"] == "1"]
    assert outgroup == "GCF_000003.1"


def test_sister_mode_runs_as_quality_without_the_gtdb_tree(
    tmp_path, table, monkeypatch, caplog
) -> None:
    def no_tree(*_args, **_kwargs):
        raise AssertionError("the GTDB reference tree was requested")

    monkeypatch.setattr(gtdb_table, "_reference_tree", no_tree)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.logger.addHandler(caplog.handler)  # the workdir logger does not propagate
    with caplog.at_level(logging.WARNING):
        metadata.run(ctx, _params(table, outgroup_mode="sister"))
    assert "--taxonomy ncbi: the sister-group outgroup follows the GTDB tree" in caplog.text
    rows = _selection(ctx.workdir)
    (outgroup,) = [acc for acc, row in rows.items() if row["is_outgroup"] == "1"]
    assert outgroup == "GCF_000005.1"  # Priestia, of the NCBI family Bacillaceae
    assert ctx.config.stages["metadata"].params["outgroup_rule"] == "quality"


def test_a_gtdb_name_is_not_found_in_the_ncbi_taxonomy(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="not found in the NCBI taxonomy of the GTDB table"):
        metadata.run(ctx, _params(table, target_genus="Bacillus_A"))


def test_rows_without_ncbi_names_are_unclassified() -> None:
    tax = gtdb_table.parse_taxonomy(f"{_PREFIX};f__;g__;s__", "ncbi")
    assert (tax["family"], tax["genus"], tax["species"]) == ("unclassified",) * 3
    named = gtdb_table.parse_taxonomy(_tax("Bacillaceae", "Bacillus", "Bacillus cereus"), "ncbi")
    assert named["species"] == "cereus" and named["species_name"] == "Bacillus cereus"
    # GTDB names are tokenised as before.
    gtdb = gtdb_table.parse_taxonomy(_tax("Bacillaceae_G", "Bacillus_A", "Bacillus_A cereus"))
    assert (gtdb["family"], gtdb["genus"], gtdb["species"]) == (
        "Bacillaceae-G",
        "Bacillus-A",
        "cereus",
    )


def test_a_table_without_the_ncbi_column_is_refused(tmp_path) -> None:
    table = _write_table(tmp_path / "bac120_metadata_r232.tsv.gz", with_ncbi=False)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="has no ncbi_taxonomy column"):
        metadata.run(ctx, _params(table))


def test_an_unknown_taxonomy_is_refused(tmp_path, table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--taxonomy must be one of gtdb, ncbi"):
        metadata.run(ctx, _params(table, taxonomy="silva"))


@pytest.mark.parametrize("source", ["api", "ncbi"])
def test_the_option_belongs_to_the_table_source(source) -> None:
    params = MetadataParams(
        dataset="all", level="genus", source=source, target_genus="Bacillus", taxonomy="ncbi"
    )
    with pytest.raises(UserInputError, match="--taxonomy cannot be used with --source"):
        metadata._validate(params)


def test_census_of_the_workdir_counts_by_ncbi_names(tmp_path, table) -> None:
    from repgenr.stages import census

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(table, limit=2))
    result = census.workdir_census(
        census.CensusParams(workdir=ctx.workdir, metadata_path=str(table)), LOG
    )
    assert result.mode == "candidates"
    assert result.totals["genomes"] == 4  # the NCBI genus, before --limit 2
    assert {row["name"] for row in result.rows} >= {"Bacillus cereus", "Bacillus subtilis"}


def _taxon_census(table: Path, **options):
    from repgenr.stages import census

    params = census.CensusParams(source="table", metadata_path=str(table), **options)
    return census.taxon_census(params, census.validate(params), LOG)


def test_census_of_a_taxon_counts_the_table_by_ncbi_names(table) -> None:
    result = _taxon_census(table, target_genus="Bacillus", taxonomy="ncbi")
    assert result.totals["genomes"] == 4
    assert {row["name"] for row in result.rows} >= {"Bacillus cereus", "Bacillus subtilis"}
    assert any("NCBI taxonomy column" in note for note in result.notes)
    # The same name by GTDB taxonomy is another, smaller genus.
    by_gtdb = _taxon_census(table, target_genus="Bacillus")
    assert by_gtdb.totals["genomes"] == 2
    assert not any("NCBI taxonomy column" in note for note in by_gtdb.notes)


def test_census_does_not_find_a_gtdb_name_in_the_ncbi_taxonomy(table) -> None:
    with pytest.raises(UserInputError, match="NCBI genus Bacillus_A"):
        _taxon_census(table, target_genus="Bacillus_A", taxonomy="ncbi")


@pytest.mark.parametrize(
    "options",
    [
        {"source": "api"},
        {"source": "ncbi"},
        {},
        {"viral": True, "target": "adenoviridae"},
    ],
)
def test_census_takes_the_taxonomy_with_the_table_source_only(options) -> None:
    from repgenr.stages import census

    options.setdefault("target_genus", "Bacillus")
    with pytest.raises(UserInputError, match="--source table only"):
        census.validate(census.CensusParams(taxonomy="ncbi", **options))


def test_census_refuses_an_unknown_taxonomy_and_the_option_with_a_workdir(tmp_path) -> None:
    from repgenr.stages import census

    with pytest.raises(UserInputError, match="--taxonomy silva is not known"):
        census.validate(
            census.CensusParams(source="table", target_genus="Bacillus", taxonomy="silva")
        )
    with pytest.raises(UserInputError, match="--taxonomy.*do not apply with -wd"):
        census.validate(census.CensusParams(workdir=tmp_path, taxonomy="ncbi"))


def test_methods_text_names_the_taxonomy() -> None:
    from repgenr.reporting.methods import _metadata
    from repgenr.reporting.model import MethodsInfo

    info = MethodsInfo(
        stage="metadata",
        tool="gtdb-table",
        params={
            "source": "tsv",
            "release": "232.0",
            "taxonomy": "ncbi",
            "level": "genus",
            "target_genus": "Bacillus",
        },
        tool_versions={},
    )
    text, _cited = _metadata(info)
    assert text.startswith(
        "Genomes were selected from the GTDB metadata table of release 232.0, by the NCBI "
        "taxonomy given in the table for Bacillus (genus)."
    )
