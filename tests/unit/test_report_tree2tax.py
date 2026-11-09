"""Report step of tree2tax."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from repgenr.core.config import Config
from repgenr.reporting.fragments import read_fragment, refresh_step
from repgenr.reporting.methods import sentence
from repgenr.reporting.model import MethodsInfo
from repgenr.reporting.registry import BY_NAME


def _wd(tmp_path: Path) -> Path:
    wd = tmp_path / "wd"
    wd.mkdir()
    (wd / "tree2tax.tsv").write_text(
        "child\tparent\trank\n"
        "n1\troot\tgenus\nn2\tn1\tspecies\nn3\tn1\tspecies\n"
        "a\tn2\tstrain\nb\tn2\tstrain\nc\tn3\tstrain\nd\tn3\tstrain\n"
    )
    (wd / "tree2tax_nodes.tsv").write_text(
        "node\tparent\tn_leaves\tdepth\tbranch_length\tsupport\tred\tlineage\n"
        "root\t\t4\tNA\tNA\tNA\tNA\t\n"
        "n1\troot\t4\t0.1\t0.1\t100\t0.4\tFam;Gen\n"
        "n2\tn1\t2\t0.2\t0.1\t95\t0.9\tFam;Gen;sp1\n"
        "n3\tn1\t2\t0.2\t0.1\t60\t0.92\tFam;Gen;sp2\n"
        "a\tn2\t1\t0.3\t0.1\tNA\tNA\tFam;Gen;sp1\n"
    )
    (wd / "genomes_map.tsv").write_text("A\ta\nB\tb\nC\tc\nD\td\n")
    cfg = Config.load(wd)
    cfg.record_stage(
        "tree2tax",
        tool="dendropy",
        tool_versions={"dendropy": "5.0.8"},
        params={"collapse_support": 0.7, "collapsed_nodes": 2, "ranks": True,
                "rank_red": "genus:0.5,species:0.85"},
        completed=datetime.now(UTC).isoformat(),
    )  # fmt: skip
    cfg.save(wd)
    return wd


def test_tree2tax_report(tmp_path: Path) -> None:
    wd = _wd(tmp_path)
    assert refresh_step(wd, BY_NAME["tree2tax"], formats=("png",)) == "written"
    frag = wd / "report" / "tree2tax"
    rep = read_fragment(frag)
    ranks = (frag / "nodes_per_rank.tsv").read_text().splitlines()
    assert ranks[1:] == ["strain\t4", "species\t2", "genus\t1"]
    assert "internal nodes collapsed\t2" in (frag / "tree2tax_settings.tsv").read_text()
    assert {f.id for f in rep.figures} == {"node_sizes", "node_support", "node_red"}
    assert "wrote 7 child-parent relations" in rep.summary[0]
    assert "mapped 4 genome accessions" in rep.summary[0]
    assert "2 internal nodes were collapsed" in rep.summary[0]


def test_tree2tax_methods() -> None:
    text, _ = sentence(
        MethodsInfo(
            "tree2tax",
            tool="dendropy",
            tool_versions={"dendropy": "5.0.8"},
            params={"collapse_support": 0.7, "collapse_diameter": 0.01, "ranks": True},
        )
    )
    assert text == (
        "The rooted tree was converted into a node taxonomy for FlexTaxD with DendroPy 5.0.8. "
        "Internal nodes with support below 0.7 were collapsed. Clades within a leaf-to-leaf "
        "distance of 0.01 were merged. Ranks were assigned from relative evolutionary "
        "divergence (RED) bands."
    )
