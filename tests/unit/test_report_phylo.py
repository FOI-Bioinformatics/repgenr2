"""Report steps of the SNP typing, the alignment and the tree; the tree figure."""

from __future__ import annotations

import json
import random
from datetime import UTC, datetime
from pathlib import Path

import dendropy

from repgenr.core.config import Config
from repgenr.core.contracts import SELECTION_TSV, SelectionRow, write_selection
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.methods import sentence, tool_version
from repgenr.reporting.model import MethodsInfo
from repgenr.reporting.registry import BY_NAME, STAGE_TO_STEPS
from repgenr.reporting.treeplot import draw_tree, layout, segments

NEWICK = (
    "((Fam_Gen_alpha_A1:0.1,Fam_Gen_alpha_A2:0.2)98:0.05,(Fam_Gen_beta_B1:0.3,OUT_X:0.4)60:0.1);"
)


def _tree(text: str = NEWICK):
    return dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True)


def test_layout_geometry() -> None:
    tree = _tree()
    lay = layout(tree)
    names = [leaf.taxon.label for leaf in lay.leaves]
    assert sorted(names) == ["Fam_Gen_alpha_A1", "Fam_Gen_alpha_A2", "Fam_Gen_beta_B1", "OUT_X"]
    assert [lay.y[id(leaf)] for leaf in lay.leaves] == [0.0, 1.0, 2.0, 3.0]
    a1 = next(leaf for leaf in lay.leaves if leaf.taxon.label == "Fam_Gen_alpha_A1")
    assert abs(lay.x[id(a1)] - 0.15) < 1e-9
    # One horizontal segment per edge below the root, one vertical per internal node.
    n_nodes = sum(1 for _ in tree.preorder_node_iter())
    n_internal = sum(1 for _ in tree.internal_nodes())
    assert len(segments(tree, lay)) == (n_nodes - 1) + n_internal
    assert lay.phylogram


def test_cladogram_without_lengths() -> None:
    lay = layout(_tree("((a,b),(c,d));"))
    assert not lay.phylogram
    assert max(lay.x.values()) == 2.0


def test_large_tree_has_no_labels_and_support_dots() -> None:
    from matplotlib.figure import Figure

    rng = random.Random(0)
    taxa = [f"t{i}" for i in range(500)]
    tree = dendropy.simulate.treesim.birth_death_tree(
        birth_rate=1.0, death_rate=0.0, num_extant_tips=500, rng=rng,
        taxon_namespace=dendropy.TaxonNamespace(taxa),
    )  # fmt: skip
    for node in tree.internal_nodes():
        node.label = str(rng.randint(40, 100))
    drawn = draw_tree(Figure(), tree, labels={}, groups={}, outgroup=None, unit="u")
    assert drawn["leaves"] == 500 and not drawn["labelled"]
    assert drawn["support_mode"] == "dots" and drawn["support_scale"] == 100


def _phylo_wd(tmp_path: Path) -> Path:
    wd = tmp_path / "wd"
    (wd / "tree" / "msa").mkdir(parents=True)
    (wd / "tree" / "tree.nwk").write_text(NEWICK + "\n")
    rows = [
        SelectionRow("A1", "Fam", "Gen", "alpha", False, "Fam_Gen_alpha_A1.fasta"),
        SelectionRow("A2", "Fam", "Gen", "alpha", False, "Fam_Gen_alpha_A2.fasta"),
        SelectionRow("B1", "Fam", "Gen", "beta", False, "Fam_Gen_beta_B1.fasta"),
        SelectionRow("X", "Out", "Other", "sp", True, "OUT_X.fasta"),
    ]
    write_selection(wd / SELECTION_TSV, rows)
    (wd / "tree" / "msa" / "msa_source.json").write_text(
        json.dumps({"records": ["A1", "A2", "B1", "X"], "v": 3})
    )
    (wd / "tree" / "msa" / "core_gene_stats.json").write_text(
        json.dumps({"core_genes": 2, "core_genes_available": 3, "alignment_length": 300,
                    "informative_sites_total": 40, "variable_sites_total": 60,
                    "core_threshold": 0.95, "genomes": 4, "alphabet": "dna"})
    )  # fmt: skip
    (wd / "tree" / "msa" / "core_genes.tsv").write_text(
        "gene\tkept\trank\tlength\tpresent\toccupancy\tgap_fraction\tvariable_sites\t"
        "informative_sites\tmean_identity\tstart\tend\n"
        "g1\ttrue\t1\t200\t4\t1.0\t0\t40\t30\t0.9\t1\t200\n"
        "g2\ttrue\t2\t100\t4\t1.0\t0\t20\t10\t0.9\t201\t300\n"
        "g3\tfalse\t\t90\t3\t0.75\t0\t5\t2\t0.9\t\t\n"
    )
    (wd / "tree" / "raxml.raxml.bestModel").write_text(
        "GTR{1/2/3/4/5/1}+FU{0.3/0.2/0.2/0.3}+G4m{0.3}+BU{0.9}, g1 = 1-200\n"
        "GTR{1/2/3/4/5/1}+FU{0.3/0.2/0.2/0.3}+G4m{0.4}+BU{1.1}, g2 = 201-300\n"
    )
    cfg = Config.load(wd)
    cfg.record_stage(
        "phylo",
        tool="raxmlng",
        params={"msa_source": "coregene", "coregene": "ppanggolin", "core_threshold": 0.95,
                "partitioned": True, "bootstrap": 100, "root": "outgroup", "outgroup": "OUT_X"},
        tool_versions={"raxml-ng": "2.0.2", "ppanggolin": "2.3.2"},
        completed=datetime.now(UTC).isoformat(),
    )  # fmt: skip
    cfg.save(wd)
    return wd


def test_phylogeny_and_alignment_steps(tmp_path: Path) -> None:
    wd = _phylo_wd(tmp_path)
    assert STAGE_TO_STEPS["phylo"] == ("alignment", "phylogeny")
    refresh_step(wd, BY_NAME["phylogeny"], formats=("png", "svg"))
    rep = read_fragment(wd / "report" / "phylogeny")
    summary = dict(
        line.split("\t")
        for line in (wd / "report" / "phylogeny" / "tree_summary.tsv").read_text().splitlines()[1:]
    )
    assert summary["leaves"] == "4"
    assert summary["tree builder"] == "RAxML-NG"
    assert summary["substitution model"] == "GTR+FU+G4m+BU per partition (2 partitions)"
    assert summary["internal branches with support"] == "2"
    assert summary["branches with support >= 95 (%)"] == "50.0"
    assert "rooted on the outgroup OUT_X" in rep.summary[0]
    tree_fig = next(f for f in rep.figures if f.id == "tree")
    assert "Numbers at the nodes are support values." in tree_fig.caption
    assert "outgroup is shown in bold" in tree_fig.caption
    svg = (wd / "report" / "phylogeny" / "tree.svg").read_text()
    # Leaves are labelled with species and accession from selection.tsv.
    assert "Gen alpha A1" in svg and "Other sp X" in svg

    refresh_step(wd, BY_NAME["alignment"], formats=("png",))
    aln = read_fragment(wd / "report" / "alignment")
    assert {f.id for f in aln.figures} == {
        "gene_occupancy",
        "informative_by_rank",
        "gene_length_informative",
    }
    assert "2 of 3 core genes (ppanggolin)" in aln.summary[0]


def test_phylo_methods_and_version_aliases() -> None:
    m = MethodsInfo(
        "phylo",
        tool="raxmlng",
        tool_versions={"raxml-ng": "2.0.2", "ppanggolin": "2.3.2"},
        params={"msa_source": "coregene", "coregene": "ppanggolin", "core_threshold": 0.95,
                "partitioned": True, "bootstrap": 100, "root": "midpoint"},
    )  # fmt: skip
    text, _ = sentence(m)
    assert text == (
        "Core genes were identified with PPanGGOLiN 2.3.2 (present in at least 95% of the "
        "genomes), aligned and concatenated. Each gene was given its own substitution model. "
        "The tree was inferred with RAxML-NG 2.0.2 with 100 bootstrap replicates. It was "
        "rooted at the midpoint."
    )
    assert tool_version({"dRep": "3.7.1"}, "drep") == "3.7.1"
    assert tool_version({}, "drep") is None


def _snp_wd(tmp_path: Path, n: int) -> Path:
    wd = tmp_path / f"snp{n}"
    (wd / "snp").mkdir(parents=True)
    names = [f"g{i}" for i in range(n)]
    rng = random.Random(n)
    m = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            m[i][j] = m[j][i] = rng.randint(1, 50)
    lines = ["\t" + "\t".join(names)] + [
        names[i] + "\t" + "\t".join(str(v) for v in m[i]) for i in range(n)
    ]
    (wd / "snp" / "snp_distance_matrix.tsv").write_text("\n".join(lines) + "\n")
    (wd / "snp" / "alignment_stats.json").write_text(
        json.dumps({"variable_sites": 120, "constant_sites": {"A": 1, "C": 2, "G": 3, "T": 4},
                    "masked": False, "core_fraction": None})
    )  # fmt: skip
    return wd


def test_snptype_step_heatmap_and_cap(tmp_path: Path) -> None:
    small = _snp_wd(tmp_path, 6)
    refresh_step(small, BY_NAME["snptype"], formats=("png",))
    rep = read_fragment(small / "report" / "snptype")
    assert {f.id for f in rep.figures} == {"snp_distance_histogram", "snp_distance_heatmap"}
    row = (small / "report" / "snptype" / "snp_distances.tsv").read_text().splitlines()[1]
    assert row.split("\t")[:2] == ["6", "15"]
    big = _snp_wd(tmp_path, 301)
    refresh_step(big, BY_NAME["snptype"], formats=("png",))
    rep = read_fragment(big / "report" / "snptype")
    assert [f.id for f in rep.figures] == ["snp_distance_histogram"]
    assert any("No distance heatmap" in n for n in rep.notes)
