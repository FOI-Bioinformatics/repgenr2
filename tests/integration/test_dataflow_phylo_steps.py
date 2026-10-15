"""Stateless data-channel steps for phylo/tree2tax (no shared workdir).

Covers the pure, binary-free parts: the tree2tax-relations step end to end
(dendropy only) and the phylo input resolvers (genome listing, reference and
outgroup resolution) that let build_tree run without a WorkdirContext.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.stages.phylo import list_fasta, resolve_outgroup_files
from repgenr.stages.tree2tax import Tree2taxStepParams, tree2tax_relations
from repgenr.tree.medoid import resolve_reference as _resolve_reference

_LOG = logging.getLogger("test")

_NWK = (
    "((Fam_Gen_sp_GCA_000001.1,Fam_Gen_sp_GCA_000002.1),"
    "(Fam_Gen_sp_GCA_000003.1,Fam_Gen_sp_GCA_000004.1));"
)


def _edges(path: Path) -> set[tuple[str, str]]:
    return {tuple(line.split("\t")) for line in path.read_text(encoding="utf-8").splitlines()[1:]}  # type: ignore[misc]


def test_tree2tax_relations_step(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    out = tmp_path / "out"
    t2t, gmap = tree2tax_relations(Tree2taxStepParams(tree=tree, out_dir=out), _LOG)
    assert t2t == out / "tree2tax.tsv"
    edges = _edges(t2t)
    leaves = {c for c, _ in edges if c.startswith("Fam_")}
    assert leaves == {
        "Fam_Gen_sp_GCA_000001.1",
        "Fam_Gen_sp_GCA_000002.1",
        "Fam_Gen_sp_GCA_000003.1",
        "Fam_Gen_sp_GCA_000004.1",
    }
    assert any(p == "root" for _, p in edges)
    accs = {ln.split("\t")[0] for ln in gmap.read_text().splitlines()}
    assert accs == {"GCA_000001.1", "GCA_000002.1", "GCA_000003.1", "GCA_000004.1"}


def test_tree2tax_relations_step_outgroup_and_dereplicated(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    og = tmp_path / "outgroup"
    og.mkdir()
    (og / "Fam_Gen_sp_GCA_000004.1.fasta").write_text(">x\nACGT\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCA_000004.1\n")

    from repgenr.core.contracts import write_clusters

    clusters = tmp_path / "clusters.tsv"
    write_clusters(clusters, {"Fam_Gen_sp_GCA_000001.1": ["Fam_Gen_sp_GCA_000009.1"]})

    out = tmp_path / "out"
    t2t, gmap = tree2tax_relations(
        Tree2taxStepParams(
            tree=tree,
            out_dir=out,
            clusters=clusters,
            outgroup_dir=og,
            outgroup_accession=acc,
            remove_outgroup=True,
            include_dereplicated=True,
        ),
        _LOG,
    )
    edges = _edges(t2t)
    assert not any(c == "Fam_Gen_sp_GCA_000004.1" for c, _ in edges)  # outgroup dropped
    accs = {ln.split("\t")[0] for ln in gmap.read_text().splitlines()}
    assert "GCA_000009.1" in accs  # redundant member mapped under its representative


def test_list_fasta(tmp_path: Path) -> None:
    (tmp_path / "a.fasta").write_text(">a\nAC\n")
    (tmp_path / "b.fna").write_text(">b\nGT\n")
    (tmp_path / ".hidden.fasta").write_text(">h\nNN\n")
    (tmp_path / "notes.txt").write_text("ignore")
    names = [p.name for p in list_fasta(tmp_path)]
    assert names == ["a.fasta", "b.fna"]
    assert list_fasta(tmp_path / "missing") == []


def test_resolve_reference(tmp_path: Path) -> None:
    g1 = tmp_path / "Fam_Gen_sp_GCA_1.fasta"
    g2 = tmp_path / "Fam_Gen_sp_GCA_2.fasta"
    for g in (g1, g2):
        g.write_text(">x\nAC\n")
    logger = logging.getLogger("test-resolve-reference")
    assert _resolve_reference(None, [g1, g2], None, logger).genome == g1  # default: first
    assert _resolve_reference("Fam_Gen_sp_GCA_2.fasta", [g1, g2], None, logger).genome == g2
    with pytest.raises(UserInputError):
        _resolve_reference("nope.fasta", [g1, g2], None, logger)


def test_resolve_outgroup_files(tmp_path: Path) -> None:
    og = tmp_path / "outgroup"
    og.mkdir()
    f = og / "Fam_Gen_sp_GCA_000004.1.fasta"
    f.write_text(">x\nAC\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCA_000004.1\n")
    found, leaf = resolve_outgroup_files(og, acc, _LOG)
    assert found == f
    assert leaf == "Fam_Gen_sp_GCA_000004.1"
    # missing dir / accession -> no outgroup, no crash
    assert resolve_outgroup_files(tmp_path / "nope", acc, _LOG) == (None, None)


def test_tree2tax_relations_rejects_a_malformed_tree(tmp_path: Path) -> None:
    # A tree file that is not Newick is an input error (exit 3) naming the
    # file, not a dendropy parser traceback (exit 1).
    from repgenr.core.errors import WorkdirError

    tree = tmp_path / "tree.nwk"
    tree.write_text("not a newick tree\n")
    with pytest.raises(WorkdirError, match="tree.nwk"):
        tree2tax_relations(Tree2taxStepParams(tree=tree, out_dir=tmp_path / "out"), _LOG)


def test_tree2tax_relations_rejects_a_missing_clusters_file(tmp_path: Path) -> None:
    # An explicit --clusters path that does not exist must not be skipped
    # silently: the genome map would then lack every dereplicated member.
    from repgenr.core.errors import WorkdirError

    tree = tmp_path / "tree.nwk"
    tree.write_text(_NWK + "\n")
    with pytest.raises(WorkdirError, match="clusters.tsv"):
        tree2tax_relations(
            Tree2taxStepParams(
                tree=tree, out_dir=tmp_path / "out", clusters=tmp_path / "clusters.tsv"
            ),
            _LOG,
        )


@pytest.mark.parametrize("n", [1, 2, 3])
def test_phylo_build_needs_three_genomes(tmp_path: Path, register_tool, n: int) -> None:
    from repgenr.core.errors import WorkdirError
    from repgenr.core.plugins import ToolCapabilities
    from repgenr.stages.phylo import PhyloBuildParams, PhyloParams, phylo_build
    from repgenr.treebuilders.base import InputKind, TreeBuilder
    from repgenr.treebuilders.base import registry as tb_registry

    calls: list[int] = []

    class _Tree(TreeBuilder):
        capabilities = ToolCapabilities(name="faketree_steps")
        input_kind = InputKind.GENOMES

        def preflight(self):
            return {}

        def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
            calls.append(1)
            out_dir.mkdir(parents=True, exist_ok=True)
            tree = out_dir / "tree.nwk"
            leaves = ",".join(Path(g).stem for g in msa_or_genomes)
            tree.write_text(f"({leaves});\n")
            return tree

    tb_registry._load()
    register_tool(tb_registry, "faketree_steps", _Tree)
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    for i in range(n):
        (genomes / f"g{i}.fasta").write_text(">s\nACGT\n")
    params = PhyloBuildParams(
        genomes_dir=genomes,
        out_dir=tmp_path / "out",
        phylo=PhyloParams(treebuilder="faketree_steps", no_outgroup=True),
    )
    logger = logging.getLogger("t")
    if n < 3:
        with pytest.raises(WorkdirError) as exc:
            phylo_build(params, logger)
        assert f"A tree needs at least 3 genomes; {n} found in {genomes}." in str(exc.value)
        assert exc.value.exit_code == 3
        assert calls == []
    else:
        assert phylo_build(params, logger).exists()


def test_tree2tax_relations_rejects_text_after_the_final_semicolon(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    tree = tmp_path / "tree.nwk"
    tree.write_text("(a,(b,c));garbage\n")
    with pytest.raises(WorkdirError, match=r"tree\.nwk.*truncated"):
        tree2tax_relations(Tree2taxStepParams(tree=tree, out_dir=tmp_path / "out"), _LOG)


def test_tree2tax_relations_outgroup_not_a_leaf_is_an_error(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    tree = tmp_path / "tree.nwk"
    tree.write_text(
        "(Fam_Gen_sp_GCA_000001.1,(Fam_Gen_sp_GCA_000002.1,Fam_Gen_sp_GCA_000003.1));\n"
    )
    og = tmp_path / "outgroup"
    og.mkdir()
    (og / "Fam_Gen_sp_GCA_000099.1.fasta").write_text(">x\nACGT\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCA_000099.1\n")
    params = Tree2taxStepParams(
        tree=tree, out_dir=tmp_path / "out", outgroup_dir=og, outgroup_accession=acc
    )
    with pytest.raises(WorkdirError, match=r"Fam_Gen_sp_GCA_000099\.1 is not a leaf"):
        tree2tax_relations(params, _LOG)


def test_tree2tax_relations_no_outgroup_ignores_the_outgroup_inputs(tmp_path: Path) -> None:
    """A tree built with phylo --no-outgroup has no outgroup leaf; --no-outgroup
    leaves it unrooted even when the outgroup files are staged."""
    tree = tmp_path / "tree.nwk"
    tree.write_text(
        "(Fam_Gen_sp_GCA_000001.1,(Fam_Gen_sp_GCA_000002.1,Fam_Gen_sp_GCA_000003.1));\n"
    )
    og = tmp_path / "outgroup"
    og.mkdir()
    (og / "Fam_Gen_sp_GCA_000099.1.fasta").write_text(">x\nACGT\n")
    acc = tmp_path / "outgroup_accession.txt"
    acc.write_text("GCA_000099.1\n")
    params = Tree2taxStepParams(
        tree=tree,
        out_dir=tmp_path / "out",
        outgroup_dir=og,
        outgroup_accession=acc,
        no_outgroup=True,
    )
    t2t, gmap = tree2tax_relations(params, _LOG)
    assert {c for c, _ in _edges(t2t) if c.startswith("Fam_")} == {
        "Fam_Gen_sp_GCA_000001.1",
        "Fam_Gen_sp_GCA_000002.1",
        "Fam_Gen_sp_GCA_000003.1",
    }


def test_tree2tax_relations_lists_segment_accessions_when_given(tmp_path: Path) -> None:
    """The step mirrors the workdir stage: with segments.tsv, every member
    segment accession of a grouped isolate maps to the isolate's leaf."""
    from repgenr.core.contracts import write_segments

    tree = tmp_path / "tree.nwk"
    tree.write_text("((Vir_gen_sp_ISO1,Vir_gen_sp_ISO2),Vir_gen_sp_ISO3);\n")
    segments = tmp_path / "segments.tsv"
    write_segments(segments, {"ISO1": ["ISO1_S", "ISO1_L"]})
    t2t, gmap = tree2tax_relations(
        Tree2taxStepParams(
            tree=tree, out_dir=tmp_path / "out", no_outgroup=True, segments=segments
        ),
        _LOG,
    )
    rows = {tuple(line.split("\t")) for line in gmap.read_text().splitlines()}
    assert ("ISO1_S", "Vir_gen_sp_ISO1") in rows and ("ISO1_L", "Vir_gen_sp_ISO1") in rows


def test_tree2tax_relations_refuses_a_missing_segments_table(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    tree = tmp_path / "tree.nwk"
    tree.write_text("(A,B,C);\n")
    with pytest.raises(WorkdirError, match="Segments table not found"):
        tree2tax_relations(
            Tree2taxStepParams(
                tree=tree,
                out_dir=tmp_path / "out",
                no_outgroup=True,
                segments=tmp_path / "none.tsv",
            ),
            _LOG,
        )
