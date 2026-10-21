"""The per-gene layer of the core-gene route: reading, metrics, ranking,
concatenation and the gene table. Pure Python; no tools."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes.select import (
    GeneMetrics,
    GeneRow,
    GeneSelection,
    Partition,
    concatenate,
    gene_metrics,
    load_gene,
    rank_genes,
    write_gene_table,
)


def _fasta(path: Path, records: dict[str, str] | list[tuple[str, str]]) -> Path:
    items = records.items() if isinstance(records, dict) else records
    path.write_text("".join(f">{n}\n{s}\n" for n, s in items), encoding="utf-8")
    return path


def _m(**kw) -> GeneMetrics:
    base = dict(
        length=100,
        present=5,
        occupancy=1.0,
        gap_fraction=0.0,
        variable_sites=0,
        informative_sites=0,
        mean_identity=1.0,
    )
    base.update(kw)
    return GeneMetrics(**base)  # type: ignore[arg-type]


# --- GeneSelection ------------------------------------------------------------


def test_selection_defaults_and_validation() -> None:
    sel = GeneSelection()
    assert (sel.max_genes, sel.rank, sel.min_length, sel.max_gaps, sel.seed) == (
        0,
        "occupancy",
        0,
        None,
        0,
    )
    for bad in (
        dict(rank="best"),
        dict(max_genes=-1),
        dict(min_length=-5),
        dict(max_gaps=1.5),
        dict(max_gaps=-0.1),
    ):
        with pytest.raises(UserInputError):
            GeneSelection(**bad)  # type: ignore[arg-type]
    assert GeneSelection(max_gaps=0.0).max_gaps == 0.0
    assert GeneSelection(max_gaps=1.0).max_gaps == 1.0


# --- load_gene ----------------------------------------------------------------


def test_load_gene_parses_names_and_upper_cases(tmp_path: Path) -> None:
    path = _fasta(
        tmp_path / "g.fasta",
        [
            ("_R_genomeA;gene_7 some description", "acgt"),
            ("genomeB;gene_9", "AC-T"),
            ("genomeC", "ACGT"),
        ],
    )
    got = load_gene(path)
    assert got == {"genomeA": "ACGT", "genomeB": "AC-T", "genomeC": "ACGT"}


def test_load_gene_joins_lines_and_keeps_the_record_with_fewest_gaps(tmp_path: Path) -> None:
    path = tmp_path / "g.fasta"
    path.write_text(">a;1\nAC\n--\n>a;2\nAC\nG-\n>a;3\n----\n>b\nACGT\n", encoding="utf-8")
    assert load_gene(path) == {"a": "ACG-", "b": "ACGT"}


def test_load_gene_refuses_unequal_lengths(tmp_path: Path) -> None:
    path = _fasta(tmp_path / "odd.fasta", {"a": "ACGT", "b": "ACG"})
    with pytest.raises(WorkdirError, match="odd.fasta"):
        load_gene(path)


def test_load_gene_refuses_an_unknown_alphabet(tmp_path: Path) -> None:
    path = _fasta(tmp_path / "g.fasta", {"a": "ACGT"})
    with pytest.raises(ValueError, match="alphabet"):
        load_gene(path, alphabet="rna")


# --- gene_metrics -------------------------------------------------------------


def test_gene_metrics_counts_sites_gaps_and_occupancy() -> None:
    seqs = {
        "a": "ACGT-",
        "b": "ACGA-",
        "c": "ATGAN",
        "d": "ATGAA",
    }
    m = gene_metrics(seqs, n_genomes=5)
    assert m.length == 5 and m.present == 4
    assert m.occupancy == pytest.approx(0.8)
    # Two '-' of 20 characters; the N is missing data, not a gap.
    assert m.gap_fraction == pytest.approx(0.1)
    # Column 1 (C,C,T,T) is variable and informative; column 3 (T,A,A,A) is
    # variable only; column 4 holds one state (A) once gaps and N are left out.
    assert m.variable_sites == 2 and m.informative_sites == 1
    # Identity over columns 0-3 (column 4 has one called base): 1, 2/6, 1, 3/6.
    assert m.mean_identity == pytest.approx((1 + 2 / 6 + 1 + 3 / 6) / 4)


def test_gene_metrics_mean_identity_by_hand() -> None:
    # Column 0: A,A,C -> C(2,2)/C(3,2) = 1/3. Column 1: A,C,G -> 0.
    m = gene_metrics({"x": "AA", "y": "AC", "z": "CG"}, n_genomes=3)
    assert m.mean_identity == pytest.approx(1 / 6)
    assert m.variable_sites == 2 and m.informative_sites == 0


def test_gene_metrics_excludes_missing_symbols() -> None:
    # N, ? and . are not states for DNA; X is not a state for protein only.
    m = gene_metrics({"a": "AN?.", "b": "CNNN", "c": "CA.A"}, n_genomes=3)
    # Column 0: A,C,C variable, not informative. Column 1: only A. Column 3: only A.
    assert m.variable_sites == 1 and m.informative_sites == 0
    assert m.gap_fraction == 0.0
    dna = gene_metrics({"a": "X", "b": "A"}, n_genomes=2)
    assert dna.variable_sites == 1
    prot = gene_metrics({"a": "X", "b": "A"}, n_genomes=2, alphabet="protein")
    assert prot.variable_sites == 0


def test_gene_metrics_without_qualifying_columns() -> None:
    m = gene_metrics({"a": "A-", "b": "-N"}, n_genomes=4)
    assert m.mean_identity == 1.0 and m.variable_sites == 0
    empty = gene_metrics({}, n_genomes=4)
    assert empty.present == 0 and empty.length == 0 and empty.occupancy == 0.0
    assert empty.gap_fraction == 0.0 and empty.mean_identity == 1.0


# --- rank_genes ---------------------------------------------------------------


def test_rank_by_each_metric_with_tie_breaks() -> None:
    metrics = {
        "g1": _m(occupancy=0.8, informative_sites=5, variable_sites=9, length=300),
        "g2": _m(occupancy=1.0, informative_sites=1, variable_sites=2, length=100),
        "g3": _m(occupancy=1.0, informative_sites=9, variable_sites=1, length=200),
        "g4": _m(occupancy=1.0, informative_sites=9, variable_sites=1, length=200),
    }
    # Ties on the key break by length descending, then gene id ascending.
    assert rank_genes(metrics, GeneSelection()) == ["g3", "g4", "g2", "g1"]
    assert rank_genes(metrics, GeneSelection(rank="informative")) == ["g3", "g4", "g1", "g2"]
    assert rank_genes(metrics, GeneSelection(rank="variable")) == ["g1", "g2", "g3", "g4"]
    assert rank_genes(metrics, GeneSelection(rank="length")) == ["g1", "g3", "g4", "g2"]


def test_random_rank_is_reproducible_by_seed() -> None:
    metrics = {f"g{i:02d}": _m() for i in range(20)}
    a = rank_genes(metrics, GeneSelection(rank="random", seed=1))
    assert a == rank_genes(metrics, GeneSelection(rank="random", seed=1))
    assert sorted(a) == sorted(metrics)
    assert a != rank_genes(metrics, GeneSelection(rank="random", seed=2))
    # Independent of the order the metrics were given in.
    reversed_metrics = dict(reversed(list(metrics.items())))
    assert a == rank_genes(reversed_metrics, GeneSelection(rank="random", seed=1))


def test_filters_and_cap() -> None:
    metrics = {
        "short": _m(length=50, occupancy=1.0),
        "gappy": _m(length=500, gap_fraction=0.4),
        "ok1": _m(length=400, gap_fraction=0.1),
        "ok2": _m(length=300, gap_fraction=0.2),
    }
    assert rank_genes(metrics, GeneSelection(min_length=100, max_gaps=0.3)) == ["ok1", "ok2"]
    # max_gaps keeps a gene exactly at the limit.
    assert "ok2" in rank_genes(metrics, GeneSelection(max_gaps=0.2))
    assert rank_genes(metrics, GeneSelection(max_genes=2)) == ["gappy", "ok1"]
    # A cap of 0 keeps every gene that passes the filters.
    assert len(rank_genes(metrics, GeneSelection(max_genes=0))) == 4
    # A gene of no columns would be an empty partition and is never kept.
    assert "empty" not in rank_genes({**metrics, "empty": _m(length=0)}, GeneSelection())


# --- concatenate --------------------------------------------------------------


def _read(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        else:
            out[name] += line
    return out


def test_concatenate_writes_every_genome_with_gaps_for_absent(tmp_path: Path) -> None:
    genes = {
        "gB": {"a": "ACG", "c": "TTT"},
        "gA": {"a": "AA", "b": "CC", "c": "GG"},
        "gC": {"b": "TTTT"},
    }
    out = tmp_path / "msa.fasta"
    parts = concatenate(["gC", "gB", "gA"], genes, ["c", "b", "a", "d"], out)
    assert parts == [Partition("gA", 1, 2), Partition("gB", 3, 5), Partition("gC", 6, 9)]
    text = out.read_text(encoding="utf-8")
    # Genomes in the given order, each one record.
    assert [ln[1:] for ln in text.splitlines() if ln.startswith(">")] == ["c", "b", "a", "d"]
    got = _read(out)
    assert got == {
        "a": "AAACG----",
        "b": "CC---TTTT",
        "c": "GGTTT----",
        "d": "---------",
    }
    # Contiguous and summing to the total width.
    assert parts[0].start == 1
    assert all(p.start == q.end + 1 for q, p in zip(parts, parts[1:], strict=False))
    assert sum(p.end - p.start + 1 for p in parts) == len(got["a"])
    # The layout follows the gene ids, not the rank order given.
    out2 = tmp_path / "msa2.fasta"
    assert concatenate(["gA", "gB", "gC"], genes, ["c", "b", "a", "d"], out2) == parts
    assert out2.read_bytes() == out.read_bytes()
    assert not (tmp_path / "msa.fasta.tmp").exists()


# --- write_gene_table ---------------------------------------------------------


def test_write_gene_table_format(tmp_path: Path) -> None:
    rows = [
        GeneRow(
            gene="gA",
            kept=True,
            rank=1,
            metrics=_m(length=2, present=3, occupancy=1.0, gap_fraction=0.0, mean_identity=0.5),
            start=1,
            end=2,
        ),
        GeneRow(
            gene="gB",
            kept=False,
            rank=None,
            metrics=_m(length=3, present=2, occupancy=2 / 3, gap_fraction=1 / 9),
            start=None,
            end=None,
        ),
    ]
    path = tmp_path / "core_genes.tsv"
    write_gene_table(path, rows)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
        "gene",
        "kept",
        "rank",
        "length",
        "present",
        "occupancy",
        "gap_fraction",
        "variable_sites",
        "informative_sites",
        "mean_identity",
        "start",
        "end",
    ]
    assert lines[1].split("\t") == [
        "gA",
        "true",
        "1",
        "2",
        "3",
        "1.0000",
        "0.0000",
        "0",
        "0",
        "0.5000",
        "1",
        "2",
    ]
    assert lines[2].split("\t") == [
        "gB",
        "false",
        "",
        "3",
        "2",
        "0.6667",
        "0.1111",
        "0",
        "0",
        "1.0000",
        "",
        "",
    ]
