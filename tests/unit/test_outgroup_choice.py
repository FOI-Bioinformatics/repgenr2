"""Outgroup choice: the GTDB reference tree, the quality order, and the
metadata stage in its three --outgroup-mode settings."""

from __future__ import annotations

import csv
import gzip
import logging
from pathlib import Path

import pytest

from repgenr.core import gtdb_tree
from repgenr.core.context import WorkdirContext
from repgenr.core.errors import UserInputError
from repgenr.core.outgroup_choice import (
    OUTGROUP_CHOICE_TSV,
    Candidate,
    choose_by_groups,
    choose_by_quality,
    quality_key,
)
from repgenr.sources import gtdb_api
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams

_LOG = logging.getLogger("test-outgroup")

# ((target A, target B), (sister C, sister D)), far E. The target genus is
# g__Tgt; its sister genus g__Sis; g__Far is one step further out.
_NEWICK = (
    "(((RS_GCF_A.1:0.01,RS_GCF_B.1:0.02)'100:g__Tgt':0.03,"
    "(GB_GCA_C.1:0.04,RS_GCF_D.1:0.05)'99:g__Sis':0.06)'98.5:f__Fam':0.1,"
    "RS_GCF_E.1:0.5)'100:d__Bacteria';"
)


def test_parse_newick_reads_leaves_lengths_and_quoted_labels() -> None:
    tree = gtdb_tree.parse_newick(_NEWICK)
    assert set(tree.leaves) == {"GCF_A.1", "GCF_B.1", "GCA_C.1", "GCF_D.1", "GCF_E.1"}
    depth = tree.depths()
    assert depth[tree.leaves["GCF_A.1"]] == pytest.approx(0.1 + 0.03 + 0.01)
    assert depth[tree.leaves["GCF_E.1"]] == pytest.approx(0.5)


def test_rings_walk_out_from_the_target() -> None:
    tree = gtdb_tree.parse_newick(_NEWICK)
    placement = gtdb_tree.place(tree, ["GCF_A.1", "GCF_B.1", "GCF_X.1"])
    assert placement is not None and placement.found == 2 and placement.intruders == 0
    rings = list(gtdb_tree.sister_rings(tree, placement))
    assert [[acc for acc, _ in ring] for ring in rings] == [["GCA_C.1", "GCF_D.1"], ["GCF_E.1"]]
    # C: up 0.03 to the family node, down 0.06 + 0.04.
    assert rings[0][0][1] == pytest.approx(0.13)
    assert rings[1][0][1] == pytest.approx(0.03 + 0.1 + 0.5)


def test_a_target_split_in_the_tree_counts_the_genomes_inside_it() -> None:
    tree = gtdb_tree.parse_newick(_NEWICK)
    placement = gtdb_tree.place(tree, ["GCF_A.1", "GCA_C.1"])
    assert placement is not None and placement.intruders == 2
    assert [acc for acc, _ in next(gtdb_tree.sister_rings(tree, placement))] == ["GCF_E.1"]
    assert gtdb_tree.place(tree, ["GCF_NONE.1"]) is None


def _c(acc: str, **kw) -> Candidate:
    base = dict(species="x", completeness=99.0, contamination=0.5, n50=100_000, contigs=50)
    base.update(kw)
    return Candidate(accession=acc, **base)  # type: ignore[arg-type]


def test_quality_order_prefers_isolates_then_score_then_type_material() -> None:
    mag = _c("mag", completeness=100.0, contamination=0.0, category="derived from metagenome")
    low = _c("low", completeness=95.0)
    best = _c("best", completeness=99.5)
    typed = _c("typed", completeness=99.5, type_material="type strain of species")
    poor = _c("poor", completeness=40.0)
    assert [c.accession for c in sorted([mag, low, poor, best, typed], key=quality_key)] == [
        "typed",
        "best",
        "low",
        "mag",
        "poor",
    ]
    assert not _c("t", type_material="not type material").is_type


def test_near_equal_scores_are_decided_by_type_material_and_n50() -> None:
    a = _c("a", completeness=99.96, contamination=0.03, n50=108_000)
    typed = _c(
        "typed",
        completeness=100.0,
        contamination=0.17,
        n50=175_000,
        type_material="type strain of species",
    )
    lower = _c("lower", completeness=97.0, n50=900_000, type_material="type strain of species")
    choice = choose_by_quality([a, lower, typed])
    assert choice is not None and choice.chosen.accession == "typed"
    # Without type material in the tie, the larger N50.
    b = _c("b", completeness=99.5, contamination=0.0, n50=300_000)
    assert choose_by_quality([a, b]).chosen.accession == "b"  # type: ignore[union-attr]


def test_an_isolate_one_group_out_is_preferred_over_a_mag_sister() -> None:
    mag = _c("mag", category="derived from metagenome")
    iso = _c("iso", completeness=98.0)
    choice = choose_by_groups([[mag], [iso]])
    assert choice is not None and choice.chosen.accession == "iso"
    assert "metagenome-assembled" in choice.reason
    # Two groups out is too far: the sister MAG is kept.
    choice = choose_by_groups([[mag], [_c("mag2", category="derived from single cell")], [iso]])
    assert choice is not None and choice.chosen.accession == "mag"
    assert "no isolate within 1 further group" in choice.reason


def test_a_group_without_an_acceptable_genome_is_passed_over() -> None:
    choice = choose_by_groups(
        [[_c("weak", completeness=30.0)], [_c("ok"), _c("better", completeness=99.9)]]
    )
    assert choice is not None and choice.chosen.accession == "better"
    assert "group 2" in choice.reason
    assert {c.accession for c in choice.considered} == {"weak", "ok", "better"}
    # No acceptable genome anywhere: the best of the nearest group.
    fallback = choose_by_groups([[_c("weak", completeness=30.0)], [_c("w2", completeness=20.0)]])
    assert fallback is not None and fallback.chosen.accession == "weak"
    assert choose_by_groups([[], []]) is None
    assert choose_by_quality([]) is None


# --- the metadata stage -----------------------------------------------------

_PREFIX = "d__Bacteria;p__P;c__C;o__O;f__Fam"
# accession, genus, species, completeness, contamination, n50, contigs, category, type
_GENOMES = [
    ("GCF_A.1", "Tgt", "a", 99, 0.1, 500000, 10, "none", "type strain of species"),
    ("GCF_B.1", "Tgt", "b", 99, 0.1, 500000, 10, "none", "not type material"),
    # First in the table, far from the target and a MAG: the old rule's pick.
    ("GCA_0E.1", "Far", "e", 97, 2.8, 63000, 150, "derived from metagenome", "not type material"),
    ("GCA_C.1", "Sis", "c", 95, 0.2, 160000, 86, "none", "type strain of species"),
    ("GCF_D.1", "Sis", "d", 99.9, 0.1, 176000, 67, "none", "type strain of species"),
]


def _write_table(path: Path) -> Path:
    header = [
        "accession",
        "gtdb_genome_representative",
        "gtdb_taxonomy",
        "ncbi_genbank_assembly_accession",
        "checkm2_completeness",
        "checkm2_contamination",
        "n50_contigs",
        "contig_count",
        "ncbi_genome_category",
        "gtdb_type_designation_ncbi_taxa",
    ]
    order = ["GCA_0E.1", "GCF_A.1", "GCF_B.1", "GCA_C.1", "GCF_D.1"]
    rows = {g[0]: g for g in _GENOMES}
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\t".join(header) + "\n")
        for acc in order:
            _, genus, sp, comp, cont, n50, contigs, cat, typ = rows[acc]
            gid = ("GB_" if acc.startswith("GCA") else "RS_") + acc
            tax = f"{_PREFIX};g__{genus};s__{genus} {sp}"
            fields = [gid, gid, tax, acc, comp, cont, n50, contigs, cat, typ]
            fo.write("\t".join(str(f) for f in fields) + "\n")
    return path


def _tree(path: Path) -> Path:
    newick = _NEWICK.replace("RS_GCF_E.1", "GB_GCA_0E.1")
    path.write_text(newick + "\n", encoding="utf-8")
    return path


def _params(tmp_path: Path, **kw) -> MetadataParams:
    base = dict(
        dataset="rep",
        level="genus",
        release="232.0",
        version="bac120",
        target_genus="Tgt",
        metadata_path=str(_write_table(tmp_path / "table.tsv.gz")),
    )
    base.update(kw)
    return MetadataParams(**base)  # type: ignore[arg-type]


def _choice_rows(wd: Path) -> list[dict]:
    with open(wd / OUTGROUP_CHOICE_TSV, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def test_sister_mode_takes_the_best_genome_of_the_sister_genus(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    tree = _tree(tmp_path / "bac120_r232.tree")
    metadata.run(ctx, _params(tmp_path, gtdb_tree_path=str(tree)))
    assert (ctx.workdir / "outgroup_accession.txt").read_text().strip() == "GCF_D.1"
    rows = _choice_rows(ctx.workdir)
    assert [r["accession"] for r in rows] == ["GCF_D.1", "GCA_C.1"]
    assert rows[0]["chosen"] == "yes" and rows[0]["group"] == "1"
    assert rows[0]["reason"] == "best quality in the sister group"
    assert float(rows[0]["tree_distance"]) == pytest.approx(0.03 + 0.06 + 0.05)
    record = ctx.config.stages["metadata"].params
    assert record["outgroup_mode"] == "sister" and record["outgroup"] == "GCF_D.1"
    assert record["outgroup_rule"] == "sister"

    from repgenr.reporting.registry import StepInput
    from repgenr.reporting.steps import selection

    out = tmp_path / "report"
    out.mkdir()
    report = selection.collect(StepInput(ctx.workdir, out, ("png",)))
    table = next(t for t in report.tables if t.id == "outgroup_candidates")
    assert table.n_rows == 2
    assert any("GCF_D.1" in line and "reference tree" in line for line in report.summary)


def test_sister_mode_downloads_the_tree_into_the_workdir(tmp_path: Path, monkeypatch) -> None:
    source = _tree(tmp_path / "source.tree")
    fetched: list[str] = []

    def fake_fetch(url, dest, logger):
        fetched.append(url)
        Path(dest).write_bytes(source.read_bytes())

    monkeypatch.setattr(gtdb_tree, "_fetch", fake_fetch)
    monkeypatch.setattr(gtdb_tree.http, "verify_md5_manifest", lambda *a, **k: True)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(tmp_path))
    assert fetched == [f"{gtdb_tree.GTDB_BASE}/release232/232.0/bac120_r232.tree"]
    assert (ctx.workdir / "bac120_r232.tree").is_file()
    assert (ctx.workdir / "outgroup_accession.txt").read_text().strip() == "GCF_D.1"
    # A second run reuses the tree.
    metadata.run(ctx, _params(tmp_path))
    assert len(fetched) == 1


def test_without_the_tree_sister_mode_ranks_the_parent_taxon_by_quality(
    tmp_path: Path, capsys
) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(tmp_path))
    assert "reference tree is not available" in capsys.readouterr().err
    assert (ctx.workdir / "outgroup_accession.txt").read_text().strip() == "GCF_D.1"
    rows = _choice_rows(ctx.workdir)
    assert rows[0]["mode"] == "quality" and rows[0]["group"] == ""
    assert ctx.config.stages["metadata"].params["outgroup_rule"] == "quality"
    assert {r["accession"] for r in rows} == {"GCF_D.1", "GCA_C.1", "GCA_0E.1"}


def test_first_mode_keeps_the_table_order_rule(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    (ctx.workdir / OUTGROUP_CHOICE_TSV).write_text("stale\n", encoding="utf-8")
    metadata.run(ctx, _params(tmp_path, outgroup_mode="first"))
    assert (ctx.workdir / "outgroup_accession.txt").read_text().strip() == "GCA_0E.1"
    assert not (ctx.workdir / OUTGROUP_CHOICE_TSV).exists()
    assert ctx.config.stages["metadata"].params["outgroup_rule"] == "first"


def test_a_named_outgroup_writes_no_choice_table(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    tree = _tree(tmp_path / "t.tree")
    metadata.run(ctx, _params(tmp_path, outgroup_accession="GCA_C.1", gtdb_tree_path=str(tree)))
    assert (ctx.workdir / "outgroup_accession.txt").read_text().strip() == "GCA_C.1"
    assert not (ctx.workdir / OUTGROUP_CHOICE_TSV).exists()


def test_outgroup_mode_and_tree_path_are_validated(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--outgroup-mode"):
        metadata.run(ctx, _params(tmp_path, outgroup_mode="nearest"))
    with pytest.raises(UserInputError, match="--gtdb-tree-path"):
        metadata.run(ctx, _params(tmp_path, gtdb_tree_path=str(tmp_path / "missing.tree")))


def test_obtain_tree_reports_a_release_without_a_tree(tmp_path: Path, monkeypatch) -> None:
    def not_found(url, dest, logger):
        raise gtdb_tree.http.HTTPStatusError("404", 404)

    monkeypatch.setattr(gtdb_tree, "_fetch", not_found)
    with pytest.raises(UserInputError, match="no reference tree"):
        gtdb_tree.obtain_reference_tree(
            tmp_path, "89.0", "bac120", tree_path=None, nodownload=False, logger=_LOG
        )


def test_api_route_uses_the_tree_and_the_cards(tmp_path: Path, monkeypatch) -> None:
    tree = _tree(tmp_path / "t.tree")

    def card(acc: str) -> dict:
        g = next(x for x in _GENOMES if x[0] == acc)
        return {
            "metadata_gene": {"checkm2_completeness": g[3], "checkm2_contamination": g[4]},
            "metadata_nucleotide": {"n50_contigs": g[5], "contig_count": g[6]},
            "metadata_ncbi": {"ncbi_genome_category": None if g[7] == "none" else g[7]},
            "gtdbTypeDesignation": g[8],
            "metadataTaxonomy": {
                "gtdbFamily": "f__Fam",
                "gtdbGenus": f"g__{g[1]}",
                "gtdbSpecies": f"s__{g[1]} {g[2]}",
            },
        }

    monkeypatch.setattr(
        gtdb_api, "_api_cards_by_accession", lambda accs, logger: {a: card(a) for a in accs}
    )
    rows = [
        {"gid": "GCF_A.1", "gtdbIsRep": True, "gtdbGenus": "g__Tgt", "gtdbFamily": "f__Fam"},
        {"gid": "GCF_B.1", "gtdbIsRep": True, "gtdbGenus": "g__Tgt", "gtdbFamily": "f__Fam"},
    ]
    params = MetadataParams(
        dataset="rep", level="genus", source="api", target_genus="Tgt", gtdb_tree_path=str(tree)
    )
    record, _ = gtdb_api._select_outgroup_via_api(
        params,
        rows,
        {"GCF_A.1", "GCF_B.1"},
        _LOG,
        target_reps=["GCF_A.1", "GCF_B.1"],
        workdir=tmp_path,
    )
    assert record.accession == "GCF_D.1"
    assert _choice_rows(tmp_path)[0]["group"] == "1"


def _sig(path: Path, hashes) -> None:
    import json
    import zipfile

    payload = [{"signatures": [{"ksize": 31, "molecule": "DNA", "mins": sorted(hashes)}]}]
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("signatures/a.sig.gz", gzip.compress(json.dumps(payload).encode()))


def test_an_outgroup_close_to_the_selection_is_reported(tmp_path: Path, caplog) -> None:
    from repgenr.core.outgroup_choice import check_outgroup_distance

    sk = tmp_path / "sketches"
    sk.mkdir()
    _sig(sk / "G1.sig.zip", range(1000))
    _sig(sk / "G2.sig.zip", range(5000, 6000))
    _sig(sk / "OG.sig.zip", range(10, 1010))  # 99% shared with G1
    _sig(sk / "FAR.sig.zip", range(9000, 10000))  # nothing shared
    genomes = [tmp_path / "G1.fasta", tmp_path / "G2.fasta"]
    ani = check_outgroup_distance(genomes, [tmp_path / "OG.fasta"], sk, _LOG)
    assert ani is not None and ani > 0.99
    assert "may lie within the target" in caplog.text
    assert check_outgroup_distance(genomes, [tmp_path / "FAR.fasta"], sk, _LOG) is None
    assert check_outgroup_distance(genomes, [tmp_path / "MISSING.fasta"], sk, _LOG) is None
