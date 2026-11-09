"""Contract round-trip tests."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.contracts import (
    REASON_UNSPECIFIED,
    glob_visible,
    list_fasta,
    parse_genome_filename,
    read_clusters,
    read_genome_status,
    read_genome_status_reasons,
    record_name,
    strip_fasta_suffix,
    write_clusters,
    write_genome_status,
    write_genomes_map,
    write_tree2tax,
)


def test_list_fasta_includes_gz_and_skips_dotfiles(tmp_path: Path) -> None:
    (tmp_path / "b.fasta").write_text(">b\nAC\n")
    (tmp_path / "a.fasta.gz").write_bytes(b"\x1f\x8b")  # gzipped genome
    (tmp_path / "c.fna").write_text(">c\nGT\n")
    (tmp_path / "._a.fasta").write_text("appledouble")  # macOS metadata
    (tmp_path / "notes.txt").write_text("ignore")
    names = [p.name for p in list_fasta(tmp_path)]
    # sorted, and the .fasta.gz genome is recognised (the bug the dedup fixed)
    assert names == ["a.fasta.gz", "b.fasta", "c.fna"]


def test_list_fasta_missing_dir(tmp_path: Path) -> None:
    assert list_fasta(tmp_path / "nope") == []


def test_strip_fasta_suffix() -> None:
    assert strip_fasta_suffix("Fam_Gen_sp_GCA_1.fasta.gz") == "Fam_Gen_sp_GCA_1"
    assert strip_fasta_suffix("Fam_Gen_sp_GCA_1.fasta") == "Fam_Gen_sp_GCA_1"
    assert strip_fasta_suffix("noext") == "noext"


def test_record_name_drops_every_fasta_suffix() -> None:
    """One name for a genome in alignments, trees and tables, whatever its suffix."""
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fa", ".fna", ".fas"):
        assert record_name(Path(f"/data/Fam_Gen_sp_GCF_1.1{suffix}")) == "Fam_Gen_sp_GCF_1.1"
    # A versioned accession keeps its version; Path.stem would leave '.fasta'
    # on a gzipped name.
    assert record_name("GCF_000008985.1") == "GCF_000008985.1"
    assert record_name(Path("x.fasta.gz")) == "x"
    assert record_name("dir/noext") == "noext"


def test_gzipped_fna_and_fa_are_genome_files(tmp_path: Path) -> None:
    """NCBI FTP delivers .fna.gz; only .fasta.gz was recognised among gzip suffixes."""
    for name in ("a.fna.gz", "b.fa.gz", "c.fasta.gz", "d.fastq.gz", "e.gz"):
        (tmp_path / name).write_bytes(b"\x1f\x8b")
    assert [p.name for p in list_fasta(tmp_path)] == ["a.fna.gz", "b.fa.gz", "c.fasta.gz"]
    assert strip_fasta_suffix("Fam_Gen_sp_GCA_1.fna.gz") == "Fam_Gen_sp_GCA_1"
    assert strip_fasta_suffix("Fam_Gen_sp_GCA_1.fa.gz") == "Fam_Gen_sp_GCA_1"
    assert parse_genome_filename("Fam_Gen_sp_GCA_1.1.fa.gz") == ("Fam", "Gen", "sp", "GCA_1.1")
    assert parse_genome_filename("GCF_000008985.1_ASM898v1_genomic.fna.gz") == (
        "",
        "",
        "",
        "GCF_000008985.1",
    )
    assert parse_genome_filename("MN908947.3.fna.gz") == ("", "", "", "MN908947.3")


def test_clusters_round_trip(tmp_path: Path) -> None:
    clusters = {
        "rep_a.fasta": ["m1.fasta", "m2.fasta"],
        "rep_b.fasta": [],
    }
    path = tmp_path / "clusters.tsv"
    write_clusters(path, clusters)

    # representative always lists itself in the file
    lines = path.read_text().splitlines()
    assert lines[0] == "representative\tmember"
    assert "rep_a.fasta\trep_a.fasta" in lines
    assert "rep_b.fasta\trep_b.fasta" in lines

    # reading back excludes the self-edge
    assert read_clusters(path) == clusters


def test_genome_status_sorted(tmp_path: Path) -> None:
    path = tmp_path / "genome_status.tsv"
    write_genome_status(path, {"z.fasta": "contained", "a.fasta": "representative"})
    lines = path.read_text().splitlines()
    assert lines[0] == "genome\tstatus\treason"
    rows = [line.split("\t") for line in lines[1:]]
    assert rows == [["a.fasta", "representative", ""], ["z.fasta", "contained", ""]]


def test_genome_status_reason_only_for_fail_qc(tmp_path: Path) -> None:
    path = tmp_path / "genome_status.tsv"
    status = {"a.fasta": "representative", "b.fasta": "fail_qc", "c.fasta": "fail_qc"}
    # A reason given for a genome that did not fail QC is not written; a
    # fail_qc genome without one (c.fasta, a plugin adapter) is "unspecified".
    reasons = {"a.fasta": "length", "b.fasta": "completeness,contamination"}
    write_genome_status(path, status, reasons)
    assert read_genome_status(path) == status
    assert read_genome_status_reasons(path) == {
        "b.fasta": "completeness,contamination",
        "c.fasta": REASON_UNSPECIFIED,
    }


def test_genome_status_without_reason_column_reads(tmp_path: Path) -> None:
    # A file written before the reason column existed.
    path = tmp_path / "genome_status.tsv"
    path.write_text("genome\tstatus\na.fasta\tfail_qc\n")
    assert read_genome_status(path) == {"a.fasta": "fail_qc"}
    assert read_genome_status_reasons(path) == {}
    # Written again (a chunk read back), the genome gets the fallback reason.
    write_genome_status(path, read_genome_status(path), read_genome_status_reasons(path))
    assert read_genome_status_reasons(path) == {"a.fasta": REASON_UNSPECIFIED}


def test_genome_status_round_trip(tmp_path: Path) -> None:
    path = tmp_path / "genome_status.tsv"
    status = {
        "a.fasta": "representative",
        "b.fasta": "contained",
        "c.fasta": "fail_qc",
    }
    write_genome_status(path, status)
    assert read_genome_status(path) == status


def test_read_genome_status_missing_file(tmp_path: Path) -> None:
    assert read_genome_status(tmp_path / "absent.tsv") == {}


def test_tree2tax_dedupes(tmp_path: Path) -> None:
    path = tmp_path / "tree2tax.tsv"
    write_tree2tax(path, [("leaf", "n1"), ("leaf", "n1"), ("n1", "root")])
    body = path.read_text().splitlines()[1:]
    assert body == ["leaf\tn1", "n1\troot"]


def test_genomes_map(tmp_path: Path) -> None:
    path = tmp_path / "genomes_map.tsv"
    write_genomes_map(path, [("GCA_000001", "Fam_gen_sp_GCA_000001")])
    assert path.read_text().strip() == "GCA_000001\tFam_gen_sp_GCA_000001"


def test_contract_tsvs_use_unix_line_endings(tmp_path: Path) -> None:
    """csv.writer defaults to \\r\\n, which leaves a stray \\r on the last
    column for awk/cut consumers; every contract writer must emit \\n."""
    write_clusters(tmp_path / "c.tsv", {"a.fasta": ["b.fasta"]})
    write_genome_status(tmp_path / "s.tsv", {"a.fasta": "representative"})
    write_tree2tax(tmp_path / "t.tsv", [("child", "parent")])
    write_genomes_map(tmp_path / "m.tsv", [("acc", "leaf")])
    for name in ("c.tsv", "s.tsv", "t.tsv", "m.tsv"):
        assert b"\r" not in (tmp_path / name).read_bytes(), name


def test_newick_is_complete_accepts_one_tree_only() -> None:
    """One tree ended by ';' passes; a second tree or trailing text does not.

    A Newick reader takes the first of two concatenated trees and ignores the
    second, so a file holding two trees must not pass as one.
    """
    from repgenr.core.contracts import newick_is_complete

    assert newick_is_complete("((a:1,b:1)0.9:1,c:1);\n")
    # ';' inside a quoted label or a comment does not end the tree.
    assert newick_is_complete("('a;b':1,[x;y]c:1,'d''e':1);")
    assert not newick_is_complete("")
    assert not newick_is_complete("((a,b),c)")
    assert not newick_is_complete("((a,b),c); extra")
    assert not newick_is_complete("((a,b),c);\n((a,c),b);\n")
    assert not newick_is_complete("(a,[open comment b);")


def test_glob_visible_skips_dotfiles(tmp_path: Path) -> None:
    (tmp_path / "sub").mkdir()
    for name in ("b.aln", "a.aln", "._a.aln", "sub/c.aln", "sub/._c.aln"):
        (tmp_path / name).write_bytes(b"x")
    assert [p.name for p in glob_visible(tmp_path, "*.aln")] == ["a.aln", "b.aln"]
    found = glob_visible(tmp_path, "*.aln", recursive=True)
    assert [p.relative_to(tmp_path).as_posix() for p in found] == ["a.aln", "b.aln", "sub/c.aln"]
