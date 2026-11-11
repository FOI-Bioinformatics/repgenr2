"""iTOL annotation files: species names, cluster sizes and composition."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.contracts import write_clusters
from repgenr.reporting.itol import build_leaves, short_species, write_itol

_HEADER = (
    "accession\taccession_ncbi\tis_outgroup\tgtdb_representative\tgtdb_taxonomy\t"
    "ncbi_taxonomy\tncbi_taxonomy_unfiltered\tncbi_strain_identifiers\n"
)


def _row(acc: str, gtdb: str, ncbi: str, unfiltered: str = "", strain: str = "") -> str:
    tax = f"d__Bacteria;g__Brucella;s__Brucella {gtdb}"
    nc = f"d__Bacteria;g__Brucella;s__{ncbi}" if ncbi else "d__Bacteria;g__Brucella"
    un = f"d__Bacteria;s__{unfiltered}" if unfiltered else ""
    return f"{acc}\t{acc}\t0\t1\t{tax}\t{nc}\t{un}\t{strain}\n"


def _workdir(tmp_path: Path) -> Path:
    wd = tmp_path / "wd"
    (wd / "tree").mkdir(parents=True)
    (wd / "derep").mkdir()
    a, b, c, og = (
        "Rhizobiaceae_Brucella_melitensis_GCF_000001.1",
        "Rhizobiaceae_Brucella_melitensis_GCF_000002.1",
        "Rhizobiaceae_Brucella_melitensis_GCF_000003.1",
        "Rhizobiaceae_Brucella_haematophila_GCF_000009.1",
    )
    (wd / "tree" / "tree.nwk").write_text(f"({og}:0.1,({a}:0.01,{b}:0.02)100:0.01);\n")
    write_clusters(
        wd / "derep" / "clusters.tsv",
        {f"{a}.fasta": [f"{c}.fasta", f"{a}.fasta"], f"{b}.fasta": []},
    )
    (wd / "metadata_selected.tsv").write_text(
        _HEADER
        + _row("GCF_000001.1", "melitensis", "Brucella abortus", strain="2308;x")
        + _row("GCF_000002.1", "melitensis", "", unfiltered="Brucella sp. 2716")
        + _row("GCF_000003.1", "melitensis", "Brucella suis")
    )
    return wd


def test_short_species_names() -> None:
    assert short_species("Brucella abortus") == "B. abortus"
    assert short_species("Brucella sp. 2716") == "Brucella sp."
    assert short_species("unknown") == "unknown"


def test_leaves_carry_species_strain_and_cluster_composition(tmp_path: Path) -> None:
    leaves = {leaf.accession: leaf for leaf in build_leaves(_workdir(tmp_path))}
    a = leaves["GCF_000001.1"]
    assert (a.species, a.strain, a.genomes) == ("B. abortus", "2308", 2)
    assert a.members == {"B. abortus": 1, "B. suis": 1}  # listed with itself once
    # No NCBI species: the unfiltered taxonomy's "Brucella sp." name.
    assert leaves["GCF_000002.1"].species == "Brucella sp."
    # The outgroup is no representative: a cluster of one, named from its file.
    og = leaves["GCF_000009.1"]
    assert (og.species, og.genomes) == ("B. haematophila", 1)
    gtdb = {leaf.accession: leaf for leaf in build_leaves(_workdir(tmp_path / "g"), "gtdb")}
    assert gtdb["GCF_000001.1"].species == "B. melitensis"


def test_the_files_follow_the_itol_formats(tmp_path: Path) -> None:
    wd = _workdir(tmp_path)
    out = write_itol(wd)
    assert out == wd / "tree" / "itol"
    labels = (out / "labels.txt").read_text().splitlines()
    assert labels[:3] == ["LABELS", "SEPARATOR TAB", "DATA"]
    assert (
        "Rhizobiaceae_Brucella_melitensis_GCF_000001.1\tB. abortus 2308 (GCF_000001.1, n=2)"
        in labels
    )
    strip = (out / "colorstrip_species.txt").read_text().splitlines()
    assert strip[0] == "DATASET_COLORSTRIP" and "DATA" in strip
    bars = (out / "bar_cluster_size.txt").read_text().splitlines()
    assert "Rhizobiaceae_Brucella_melitensis_GCF_000001.1\t2" in bars
    pie = (out / "pie_cluster_species.txt").read_text().splitlines()
    fields = next(line for line in pie if line.startswith("FIELD_LABELS")).split("\t")[1:]
    row = next(
        line for line in pie if line.startswith("Rhizobiaceae_Brucella_melitensis_GCF_000001")
    )
    counts = dict(zip(fields, row.split("\t")[3:], strict=True))
    assert counts["B. abortus"] == "1" and counts["B. suis"] == "1"
    table = (out / "clusters_species.tsv").read_text().splitlines()
    assert table[0].split("\t") == [
        "leaf",
        "accession",
        "species",
        "strain",
        "genomes",
        "member_species",
    ]


def test_the_command_writes_the_files_and_refuses_unknown_names(tmp_path: Path) -> None:
    wd = _workdir(tmp_path)
    runner = CliRunner()
    result = runner.invoke(app, ["itol", "-wd", str(wd), "--out", str(tmp_path / "o")])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "o" / "labels.txt").is_file()
    assert runner.invoke(app, ["itol", "-wd", str(wd), "--names", "silva"]).exit_code == 2
    (wd / "tree" / "tree.nwk").unlink()
    result = runner.invoke(app, ["itol", "-wd", str(wd)])
    assert result.exit_code != 0 and "run phylo first" in result.output
