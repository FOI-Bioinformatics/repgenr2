"""Per-cluster summary: one row per representative over clusters.tsv + quality."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.contracts import read_cluster_summary, write_cluster_summary
from repgenr.stages.cluster_summary import summarise_clusters

REP = "Fam_Gen_tularensis_GCF_1.1.fasta"
M1 = "Fam_Gen_tularensis_GCF_2.1.fasta"
M2 = "Fam_Gen_holarctica_GCF_3.1.fasta"
SOLO = "Fam_Gen_novicida_GCF_4.1.fasta"


def test_summary_counts_members_and_species() -> None:
    rows = summarise_clusters({REP: [M1, M2], SOLO: []}, {})
    by_rep = {r.representative: r for r in rows}
    assert by_rep[REP].n_members == 2
    assert by_rep[REP].n_species == 2
    assert by_rep[REP].species == "tularensis,holarctica"
    assert by_rep[SOLO].n_members == 0
    assert by_rep[SOLO].n_species == 1


def test_summary_sorted_by_size_then_name() -> None:
    rows = summarise_clusters({"Zed_Gen_x_GCF_9.1.fasta": [], SOLO: [], REP: [M1, M2]}, {})
    assert [r.representative for r in rows] == [REP, SOLO, "Zed_Gen_x_GCF_9.1.fasta"]


def test_summary_quality_columns_blank_without_quality() -> None:
    (row,) = summarise_clusters({REP: [M1]}, {})
    assert row.rep_completeness is None
    assert row.rep_contamination is None
    assert row.member_max_completeness is None
    assert row.member_min_contamination is None
    assert row.best_member == ""


def test_summary_best_member_is_representative_when_keeper_is_best() -> None:
    quality = {REP: (99.0, 0.1), M1: (95.0, 1.0), M2: (98.0, 0.5)}
    (row,) = summarise_clusters({REP: [M1, M2]}, quality)
    assert row.rep_completeness == 99.0
    assert row.rep_contamination == 0.1
    assert row.member_max_completeness == 98.0
    assert row.member_min_contamination == 0.5
    assert row.best_member == REP


def test_summary_flags_better_member_than_keeper() -> None:
    quality = {REP: (90.0, 3.0), M1: (99.0, 0.2)}
    (row,) = summarise_clusters({REP: [M1]}, quality)
    assert row.best_member == M1


def test_summary_unscored_members_are_ignored_in_extremes() -> None:
    quality = {REP: (99.0, 0.1), M1: (95.0, 1.0)}
    (row,) = summarise_clusters({REP: [M1, M2]}, quality)
    assert row.member_max_completeness == 95.0
    assert row.member_min_contamination == 1.0
    assert row.n_members == 2


def test_summary_round_trip(tmp_path: Path) -> None:
    quality = {REP: (99.0, 0.1), M1: (95.0, 1.0)}
    rows = summarise_clusters({REP: [M1, M2], SOLO: []}, quality)
    path = tmp_path / "cluster_summary.tsv"
    write_cluster_summary(path, rows)
    header = path.read_text().splitlines()[0].split("\t")
    assert header == [
        "representative",
        "n_members",
        "n_species",
        "species",
        "rep_completeness",
        "rep_contamination",
        "member_max_completeness",
        "member_min_contamination",
        "best_member",
        "n_genomes",
        "rep_n50",
        "best_score",
        "rep_is_gtdb_representative",
        "refined",
    ]
    assert read_cluster_summary(path) == rows


def test_summary_does_not_count_unparsed_names_as_a_species() -> None:
    # Non-canonical filenames carry no species; a blank must not be counted
    # or listed as one (it gave n_species 2 and "tularensis," here).
    (row,) = summarise_clusters({REP: ["sample_1.fa", "iso.v2.fasta"]}, {})
    assert (row.n_species, row.species) == (1, "tularensis")
    (solo,) = summarise_clusters({"sample_1.fa": []}, {})
    assert (solo.n_species, solo.species) == (0, "")


def test_summary_prefers_manifest_taxonomy_over_the_filename() -> None:
    # ingest --selection keeps non-canonical filenames but records the
    # taxonomy in the manifest; the summary must use it.
    taxonomy = {
        "iso-1.fasta": ("Francisella", "Francisella tularensis"),
        "iso-2.fasta": ("Francisella", "Francisella novicida"),
    }
    canonical = "Francisellaceae_Francisella_tularensis_GCF_5.1.fasta"
    (row,) = summarise_clusters({"iso-1.fasta": ["iso-2.fasta", canonical]}, {}, taxonomy)
    # ``canonical`` is not in the taxonomy: its filename supplies the species.
    assert (row.n_species, row.species) == (2, "tularensis,novicida")


def test_summary_tells_equal_epithets_of_different_genera_apart() -> None:
    taxonomy = {
        "a.fasta": ("Escherichia", "coli"),
        "b.fasta": ("Campylobacter", "coli"),
    }
    (row,) = summarise_clusters({"a.fasta": ["b.fasta"]}, {}, taxonomy)
    assert (row.n_species, row.species) == (2, "Escherichia coli,Campylobacter coli")


def test_summary_n_genomes_counts_the_keeper(tmp_path: Path) -> None:
    rows = summarise_clusters({REP: [M1, M2], SOLO: []}, {})
    by_rep = {r.representative: r for r in rows}
    assert (by_rep[REP].n_members, by_rep[REP].n_genomes) == (2, 3)
    assert (by_rep[SOLO].n_members, by_rep[SOLO].n_genomes) == (0, 1)
    path = tmp_path / "cluster_summary.tsv"
    write_cluster_summary(path, rows)
    header, first = (line.split("\t") for line in path.read_text().splitlines()[:2])
    assert first[header.index("n_genomes")] == "3"
    assert header.index("n_genomes") == 9  # earlier column positions unchanged


def test_summary_species_list_is_capped_most_frequent_first() -> None:
    # Keeper's species first, then by member count, then name; at most five
    # names and "+N more"; n_species stays exact.
    members = [f"F_G_sp{i}_GCF_{i}.1.fasta" for i in range(1, 8)]
    members += ["F_G_sp7_GCF_70.1.fasta", "F_G_sp7_GCF_71.1.fasta", "F_G_sp3_GCF_30.1.fasta"]
    (row,) = summarise_clusters({"F_G_sp9_GCF_9.1.fasta": members}, {})
    assert row.n_species == 8
    assert row.species == "sp9,sp7,sp3,sp1,sp2,+3 more"


def test_summary_manifest_species_with_blank_genus_matches_canonical_members() -> None:
    # ingest --selection requires neither genus nor species: a row with a
    # species and a blank genus takes the genus from the filename, and the
    # name is never written with a leading genus and a space.
    canonical = "Francisellaceae_Francisella_tularensis_GCF_5.1.fasta"
    for species in ("tularensis", "Francisella tularensis"):
        taxonomy = {canonical: ("", species), "iso-1.fasta": ("", species)}
        (row,) = summarise_clusters({canonical: ["iso-1.fasta"]}, {}, taxonomy)
        assert (row.n_species, row.species) == (1, "tularensis"), species
    (row,) = summarise_clusters({"iso-1.fasta": []}, {}, {"iso-1.fasta": ("", "coli")})
    assert (row.n_species, row.species) == (1, "coli")


def test_taxonomy_lookup_reads_the_manifest(tmp_path: Path) -> None:
    from repgenr.core.context import WorkdirContext
    from repgenr.core.manifest import GenomeRecord
    from repgenr.stages.cluster_summary import taxonomy_lookup

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.manifest.upsert_many(
        [
            GenomeRecord(
                accession="A1", filename="iso-1.fasta", genus="Escherichia", species="coli"
            ),
            GenomeRecord(accession="A2", filename="iso-2.fasta", genus=None, species=None),
            GenomeRecord(
                accession="A3",
                filename="out.fasta",
                genus="Shigella",
                species="flexneri",
                is_outgroup=True,
            ),
        ]
    )
    assert taxonomy_lookup(ctx) == {
        "iso-1.fasta": ("Escherichia", "coli"),
        "iso-2.fasta": ("", ""),
        "out.fasta": ("Shigella", "flexneri"),
    }
    # No manifest at all: the filenames supply the species.
    assert taxonomy_lookup(WorkdirContext(tmp_path / "none")) == {}


def test_summary_reads_no_n50_for_unscored_or_single_scored_clusters() -> None:
    """rep_n50 is filled only for a scored keeper; an unscored cluster and a
    cluster without a second scored genome to compare read nothing else."""
    asked: list[str] = []

    def n50(name: str) -> int | None:
        asked.append(name)
        return 1000

    rows = summarise_clusters({REP: [M1, M2], SOLO: []}, {SOLO: (99.0, 0.0)}, None, n50)
    by_rep = {r.representative: r for r in rows}
    assert by_rep[REP].rep_n50 is None
    assert by_rep[SOLO].rep_n50 == 1000
    assert set(asked) == {SOLO}
