"""Selection summary (metadata_summary.tsv and plots) and metadata_selected.tsv.

RepGenR 1 wrote these files from its metadata command; RepGenR2 writes the
summary from every entry stage and the per-genome lineage table from metadata.
"""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path

import pytest

from repgenr.core import selection_summary as ss
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SelectionRow, genome_filename
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams

LOG = logging.getLogger("test.selection_summary")


def _row(acc: str, family: str, genus: str, species: str, *, outgroup=False, rep=False):
    return SelectionRow(
        acc,
        family,
        genus,
        species,
        outgroup,
        genome_filename(family, genus, species, acc),
        gtdb_representative=rep,
    )


def _read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


ROWS = [
    _row("G1", "Yersiniaceae", "Yersinia", "pestis"),
    _row("G2", "Yersiniaceae", "Yersinia", "pestis"),
    _row("G3", "Yersiniaceae", "Yersinia", "aldovae"),
    _row("G4", "Yersiniaceae", "Serratia", "marcescens"),
    _row("G5", "", "", ""),
    _row("O1", "Morganellaceae", "Xenorhabdus", "anantnagensis", outgroup=True),
]


def test_count_taxa_counts_the_ingroup_per_level_largest_first() -> None:
    counts = ss.count_taxa(ROWS)
    assert counts == [
        ("family", "Yersiniaceae", 4),
        ("family", "unclassified", 1),
        ("genus", "Yersinia", 3),
        ("genus", "Serratia", 1),
        ("genus", "unclassified", 1),
        ("species", "Yersinia pestis", 2),
        ("species", "Serratia marcescens", 1),
        ("species", "Yersinia aldovae", 1),
        ("species", "unclassified", 1),
    ]
    # Each level adds up to the ingroup size; the outgroup is not counted.
    for level in ss.LEVELS:
        assert sum(n for lvl, _, n in counts if lvl == level) == 5


def test_viral_species_token_is_the_whole_name() -> None:
    # vgenome writes the full species name as the token (spaces as hyphens).
    rows = [
        _row("V1", "Picornaviridae", "Hepatovirus", "Hepatovirus-ahepa"),
        _row("V2", "Picornaviridae", "Hepatovirus", "Rodent-hepatovirus"),
    ]
    species = [n for lvl, n, _ in ss.count_taxa(rows, full_species=True) if lvl == "species"]
    assert species == ["Hepatovirus-ahepa", "Rodent-hepatovirus"]
    # The bacterial rule names the epithet with its genus.
    assert ("species", "Hepatovirus Hepatovirus-ahepa", 1) in ss.count_taxa(rows)


def test_summary_table_and_plots_are_written(tmp_path) -> None:
    ss.write_selection_summary(tmp_path, ROWS, LOG)
    table = _read_tsv(tmp_path / ss.SUMMARY_TSV)
    assert list(table[0]) == ["level", "name", "genomes"]
    assert table[0] == {"level": "family", "name": "Yersiniaceae", "genomes": "4"}
    for name in (ss.PLOT_IN_LEVEL, ss.PLOT_PER_LEVEL):
        data = (tmp_path / name).read_bytes()
        assert data.startswith(b"\x89PNG"), name
    assert not list(tmp_path.glob("*.tmp"))


def test_many_species_are_cut_in_the_plot_not_in_the_table(tmp_path) -> None:
    rows = [_row(f"G{i}", "F", "Genus", f"sp{i:03d}") for i in range(ss.MAX_BARS_PER_LEVEL + 15)]
    ss.write_selection_summary(tmp_path, rows, LOG)
    species = [r for r in _read_tsv(tmp_path / ss.SUMMARY_TSV) if r["level"] == "species"]
    assert len(species) == ss.MAX_BARS_PER_LEVEL + 15
    assert (tmp_path / ss.PLOT_IN_LEVEL).exists()


def test_empty_selection_removes_the_summary(tmp_path) -> None:
    ss.write_selection_summary(tmp_path, ROWS, LOG)
    # Only an outgroup left (an assemble call that excused every run).
    ss.write_selection_summary(tmp_path, ROWS[-1:], LOG)
    for name in (ss.SUMMARY_TSV, ss.PLOT_IN_LEVEL, ss.PLOT_PER_LEVEL):
        assert not (tmp_path / name).exists()


def test_plot_failure_is_a_warning_and_keeps_the_table(tmp_path, monkeypatch, caplog) -> None:
    def broken(*args, **kwargs):
        raise RuntimeError("no renderer")

    monkeypatch.setattr(ss, "_plot_per_level", broken)
    with caplog.at_level(logging.WARNING, logger=LOG.name):
        ss.write_selection_summary(tmp_path, ROWS, LOG)
    assert (tmp_path / ss.SUMMARY_TSV).exists()
    assert not (tmp_path / ss.PLOT_IN_LEVEL).exists()
    assert not (tmp_path / ss.PLOT_PER_LEVEL).exists()
    assert "no renderer" in caplog.text


def test_metadata_selected_is_trimmed_to_the_current_selection(tmp_path) -> None:
    details = {"G1": ss.GenomeDetails(accession_ncbi="GCA_1"), "G2": ss.GenomeDetails()}
    ss.write_metadata_selected(tmp_path, ROWS[:2], details)
    # Another entry stage keeps G1 only (ingest --from-workdir, say).
    ss.write_selection_summary(tmp_path, [ROWS[0], _row("L1", "F", "G", "s")], LOG)
    rows = _read_tsv(tmp_path / ss.METADATA_SELECTED_TSV)
    assert [r["accession"] for r in rows] == ["G1"]
    assert rows[0]["accession_ncbi"] == "GCA_1"
    # A selection without any of its genomes removes it.
    ss.write_selection_summary(tmp_path, [_row("L1", "F", "G", "s")], LOG)
    assert not (tmp_path / ss.METADATA_SELECTED_TSV).exists()


# --- metadata stage ---------------------------------------------------------

_PREFIX = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales"
_HEADER = (
    "accession",
    "gtdb_genome_representative",
    "gtdb_taxonomy",
    "ncbi_genbank_assembly_accession",
    "ncbi_taxonomy",
    "ncbi_taxonomy_unfiltered",
    "ncbi_strain_identifiers",
)


def _gtdb_row(acc: str, family: str, genus: str, species: str, strain: str) -> str:
    gtdb = f"{_PREFIX};f__{family};g__{genus};s__{genus} {species}"
    ncbi = f"{_PREFIX};f__Yersiniaceae;g__{genus};s__{genus} {species}"
    unfiltered = ncbi.replace("p__", "x__Pseudomonadati;p__")
    gca = acc.replace("GCF_", "GCA_")
    return "\t".join((f"RS_{acc}", f"RS_{acc}", gtdb, gca, ncbi, unfiltered, strain))


@pytest.fixture()
def gtdb_table(tmp_path) -> Path:
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    lines = [
        "\t".join(_HEADER),
        _gtdb_row("GCF_000001.1", "Enterobacteriaceae", "Yersinia", "pestis", "CO92"),
        _gtdb_row("GCF_000002.1", "Enterobacteriaceae", "Yersinia", "aldovae", "none"),
        _gtdb_row("GCF_000010.1", "Enterobacteriaceae", "Serratia", "marcescens", "Db11"),
    ]
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    return path


def test_metadata_writes_summary_and_lineage_table(tmp_path, gtdb_table) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(
        ctx,
        MetadataParams(
            dataset="all",
            level="genus",
            release="232.0",
            version="bac120",
            target_genus="Yersinia",
            metadata_path=str(gtdb_table),
        ),
    )
    summary = _read_tsv(ctx.workdir / ss.SUMMARY_TSV)
    assert {(r["level"], r["name"], r["genomes"]) for r in summary} == {
        ("family", "Enterobacteriaceae", "2"),
        ("genus", "Yersinia", "2"),
        ("species", "Yersinia pestis", "1"),
        ("species", "Yersinia aldovae", "1"),
    }
    assert (ctx.workdir / ss.PLOT_IN_LEVEL).exists()
    assert (ctx.workdir / ss.PLOT_PER_LEVEL).exists()

    table = {r["accession"]: r for r in _read_tsv(ctx.workdir / ss.METADATA_SELECTED_TSV)}
    assert list(next(iter(table.values()))) == list(ss.METADATA_SELECTED_COLUMNS)
    assert set(table) == {"GCF_000001.1", "GCF_000002.1", "GCF_000010.1"}
    pestis = table["GCF_000001.1"]
    assert pestis["accession_ncbi"] == "GCA_000001.1"
    assert pestis["is_outgroup"] == "0"
    assert pestis["gtdb_representative"] == "1"
    assert pestis["gtdb_taxonomy"].endswith("f__Enterobacteriaceae;g__Yersinia;s__Yersinia pestis")
    assert "f__Yersiniaceae" in pestis["ncbi_taxonomy"]
    assert "x__Pseudomonadati" in pestis["ncbi_taxonomy_unfiltered"]
    assert pestis["ncbi_strain_identifiers"] == "CO92"
    assert table["GCF_000002.1"]["ncbi_strain_identifiers"] == ""  # GTDB writes "none"
    assert table["GCF_000010.1"]["is_outgroup"] == "1"


def test_metadata_without_ncbi_columns_leaves_them_empty(tmp_path) -> None:
    # A reduced table (as in other tests) has no NCBI taxonomy columns.
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    tax = f"{_PREFIX};f__Enterobacteriaceae;g__Yersinia;s__Yersinia pestis"
    other = f"{_PREFIX};f__Enterobacteriaceae;g__Serratia;s__Serratia marcescens"
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("accession\tgtdb_genome_representative\tgtdb_taxonomy\n")
        fo.write(f"RS_GCF_1.1\tRS_GCF_1.1\t{tax}\nRS_GCF_9.1\tRS_GCF_9.1\t{other}\n")
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(
        ctx,
        MetadataParams(
            dataset="all",
            level="genus",
            release="232.0",
            version="bac120",
            target_genus="Yersinia",
            metadata_path=str(path),
        ),
    )
    row = {r["accession"]: r for r in _read_tsv(ctx.workdir / ss.METADATA_SELECTED_TSV)}["GCF_1.1"]
    assert row["gtdb_taxonomy"] == tax
    assert row["ncbi_taxonomy"] == "" and row["ncbi_strain_identifiers"] == ""


def test_api_details_come_from_the_row_and_card() -> None:
    row = {
        "gid": "GCF_000833475.1",
        "gtdbDomain": "d__Bacteria",
        "gtdbPhylum": "p__Pseudomonadota",
        "gtdbClass": "c__Gammaproteobacteria",
        "gtdbOrder": "o__Francisellales",
        "gtdbFamily": "f__Francisellaceae",
        "gtdbGenus": "g__Francisella",
        "gtdbSpecies": "s__Francisella tularensis",
    }
    card = {
        "metadata_ncbi": {
            "ncbi_genbank_assembly_accession": "GCA_000833475.1",
            "ncbi_strain_identifiers": "NIH B-38",
        },
        "metadataTaxonomy": {
            "ncbi_taxonomy": "d__Bacteria;f__Francisellaceae",
            "ncbi_taxonomy_unfiltered": None,
        },
    }
    d = metadata._api_details(row, card)
    assert d.accession_ncbi == "GCA_000833475.1"
    assert d.gtdb_taxonomy.startswith("d__Bacteria;p__Pseudomonadota;")
    assert d.gtdb_taxonomy.endswith(";s__Francisella tularensis")
    assert d.ncbi_taxonomy == "d__Bacteria;f__Francisellaceae"
    assert d.ncbi_taxonomy_unfiltered == ""
    assert d.ncbi_strain_identifiers == "NIH B-38"
    # A failed card leaves the NCBI fields empty, the GTDB lineage from the row.
    empty = metadata._api_details(row, {})
    assert empty.ncbi_taxonomy == "" and empty.gtdb_taxonomy == d.gtdb_taxonomy
