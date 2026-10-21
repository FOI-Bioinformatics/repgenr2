"""ASTER tree builder: gene slicing, gene-tree batches, the ASTER call."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.treebuilders import aster as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test-aster")

# Five genomes over three genes: g1 columns 1-4, g2 5-8, g3 9-10.
# In g2 genome e has no called nucleotide and genome d only N/gap, so g2
# keeps a, b, c (three records) and is skipped.
_ALN = {
    "a": "ACGT" + "ACGT" + "AC",
    "b": "ACGA" + "ACGA" + "AG",
    "c": "ACTA" + "ACTA" + "TG",
    "d": "TCTA" + "NN--" + "TC",
    "e": "TCTT" + "----" + "TT",
}
_PARTS = (("g1", 1, 4), ("g/2", 5, 8), ("g3", 9, 10))


def _msa(tmp_path: Path) -> Path:
    msa = tmp_path / "aln.fasta"
    msa.write_text("".join(f">{k}\n{v}\n" for k, v in _ALN.items()), encoding="utf-8")
    return msa


class _Fake:
    """Fake run_chain/run_tool: writes each gene tree and the species tree."""

    def __init__(self, label: str = ")/0.9") -> None:
        self.chains: list[list[list[str]]] = []
        self.tools: list[list[str]] = []
        self.label = label

    def run_chain(self, caps, steps, **k) -> None:
        self.chains.append([[str(c) for c in cmd] for _, cmd in steps])
        for _, cmd in steps:
            cmd = [str(c) for c in cmd]
            if cmd[0] == "iqtree":
                out = Path(cmd[cmd.index("-pre") + 1] + ".treefile")
            else:
                out = Path(cmd[cmd.index("-out") + 1])
            names = [
                ln[1:].strip()
                for ln in Path(cmd[-1] if cmd[0] != "iqtree" else cmd[cmd.index("-s") + 1])
                .read_text(encoding="utf-8")
                .splitlines()
                if ln.startswith(">")
            ]
            out.write_text(
                f"({names[0]}:0.1,{names[1]}:0.1,({','.join(names[2:])}){self.label}:0.2);\n",
                encoding="utf-8",
            )

    def run_tool(self, caps, cmd, **k) -> int:
        cmd = [str(c) for c in cmd]
        self.tools.append(cmd)
        Path(cmd[cmd.index("-o") + 1]).write_text("((a,b)0.9:1,(c,(d,e)1:2));\n", encoding="utf-8")
        return 0


@pytest.fixture
def fake(monkeypatch) -> _Fake:
    f = _Fake()
    monkeypatch.setattr(mod, "run_chain", f.run_chain)
    monkeypatch.setattr(mod, "run_tool", f.run_tool)
    monkeypatch.setattr(mod, "preflight", lambda caps: {caps.name: "1.0"})
    monkeypatch.setattr(mod, "_program_version", lambda caps, program: "1.25.3")
    return f


def _build(tmp_path: Path, **kw) -> Path:
    params = TreeParams(partitions=_PARTS, **kw)
    return mod.AsterBuilder().build(_msa(tmp_path), tmp_path / "out", params, _LOG)


def test_partitions_are_cut_into_gene_files(tmp_path: Path, fake: _Fake, caplog) -> None:
    with caplog.at_level(logging.INFO):
        _build(tmp_path)
    genes = tmp_path / "out" / "genes"
    assert sorted(p.name for p in genes.iterdir()) == ["g1.fasta", "g3.fasta"]
    assert (genes / "g1.fasta").read_text(encoding="utf-8") == "".join(
        f">{k}\n{v[0:4]}\n" for k, v in _ALN.items()
    )
    assert (genes / "g3.fasta").read_text(encoding="utf-8").startswith(">a\nAC\n>b\nAG\n")
    assert "2 of 3 gene partitions used" in caplog.text
    assert "1 with fewer than 4 genomes" in caplog.text


def test_records_without_a_called_character_are_left_out(tmp_path: Path) -> None:
    genes, short, few, n = mod._write_genes(
        _msa(tmp_path), (("g2", 5, 8), ("g1", 1, 4)), 0, tmp_path / "genes"
    )
    assert [g.name for g in genes] == ["g1"] and few == 1 and short == 0 and n == 5


def test_four_records_with_called_characters_are_kept(tmp_path: Path) -> None:
    msa = tmp_path / "m.fasta"
    msa.write_text(">a\nAC\n>b\nNA\n>c\n-G\n>d\nT?\n>e\nN-\n", encoding="utf-8")
    (gene,), _, few, _ = mod._write_genes(msa, (("x", 1, 2),), 0, tmp_path / "genes")
    assert gene.records == 4 and few == 0
    assert ">e\n" not in gene.fasta.read_text(encoding="utf-8")


def test_min_gene_length_skips_short_partitions(tmp_path: Path, fake: _Fake) -> None:
    _build(tmp_path, extra={"min_gene_length": 3})
    assert [p.name for p in (tmp_path / "out" / "genes").iterdir()] == ["g1.fasta"]


def test_a_partition_outside_the_alignment_is_an_error(tmp_path: Path) -> None:
    from repgenr.core.errors import WorkdirError

    with pytest.raises(WorkdirError, match="outside"):
        mod._write_genes(_msa(tmp_path), (("g", 5, 11),), 0, tmp_path / "genes")


def test_iqtree_gene_trees_run_in_batches(tmp_path: Path, fake: _Fake) -> None:
    _build(tmp_path, threads=8)
    # Two genes, eight threads: two batches of one tree, one thread each.
    assert len(fake.chains) == 2
    cmd = fake.chains[0][0]
    assert cmd[0] == "iqtree"
    assert cmd[cmd.index("-m") + 1] == "GTR+G"
    assert cmd[cmd.index("-T") + 1] == "1"
    assert cmd[cmd.index("-st") + 1] == "DNA"
    assert "--fast" in cmd and "--abayes" in cmd
    assert cmd[cmd.index("-s") + 1].endswith("genes/g1.fasta")


def test_one_thread_makes_one_batch(tmp_path: Path, fake: _Fake) -> None:
    _build(tmp_path, threads=1)
    assert len(fake.chains) == 1 and len(fake.chains[0]) == 2


def test_model_and_protein_reach_iqtree(tmp_path: Path, fake: _Fake) -> None:
    _build(tmp_path, alphabet="protein", extra={"model": "WAG"})
    cmd = fake.chains[0][0]
    assert cmd[cmd.index("-m") + 1] == "WAG"
    assert cmd[cmd.index("-st") + 1] == "AA"


def test_fasttree_gene_trees(tmp_path: Path, fake: _Fake, monkeypatch) -> None:
    monkeypatch.setattr("repgenr.core.bindirs.host_which", lambda *a: None)
    _build(tmp_path, extra={"gene_tree_builder": "fasttree"})
    cmd = fake.chains[0][0]
    assert cmd[:4] == ["fasttree", "-nt", "-gtr", "-quiet"]
    assert cmd[cmd.index("-out") + 1].endswith("gene_trees/g1.nwk")
    # FastTree supports are SH-like, on a [0, 1] scale.
    assert fake.tools[0][:2] == ["wastral", "-L"]


def test_the_aster_call(tmp_path: Path, fake: _Fake) -> None:
    _build(tmp_path, threads=6)
    (cmd,) = fake.tools
    out = tmp_path / "out"
    assert cmd == [
        "wastral",
        "-B",
        "-i",
        str(out / "gene_trees.nwk"),
        "-o",
        str(out / "aster.nwk"),
        "-t",
        "6",
    ]


@pytest.mark.parametrize("program", ["astral4", "astral-pro3"])
def test_other_programs_take_no_support_preset(tmp_path: Path, fake: _Fake, program) -> None:
    _build(tmp_path, extra={"program": program, "rounds": 8})
    (cmd,) = fake.tools
    assert cmd[0] == program and "-B" not in cmd and "-L" not in cmd
    assert cmd[cmd.index("-r") + 1] == "8" and cmd[cmd.index("-s") + 1] == "8"


def test_gene_trees_are_collected_one_per_line(tmp_path: Path, fake: _Fake) -> None:
    _build(tmp_path)
    lines = (tmp_path / "out" / "gene_trees.nwk").read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    # IQ-TREE's "/0.9" aBayes label is reduced to the number ASTER parses.
    assert all(")0.9:" in line and "/" not in line for line in lines)


def test_the_species_tree_is_published(tmp_path: Path, fake: _Fake) -> None:
    tree = _build(tmp_path)
    assert tree == tmp_path / "out" / "tree.nwk"
    assert tree.read_text(encoding="utf-8") == "((a,b)0.9:1,(c,(d,e)1:2));\n"


def test_without_partitions_it_is_refused(tmp_path: Path, fake: _Fake) -> None:
    with pytest.raises(UserInputError, match="--msa-source coregene --partitioned"):
        mod.AsterBuilder().build(_msa(tmp_path), tmp_path / "out", TreeParams(), _LOG)
    assert fake.chains == [] and fake.tools == []


@pytest.mark.parametrize(
    ("extra", "match"),
    [
        ({"program": "astral3"}, "program='astral3'"),
        ({"gene_tree_builder": "raxml"}, "gene_tree_builder='raxml'"),
        ({"rounds": 0}, "rounds must be 1 or more"),
        ({"min_gene_length": "x"}, "min_gene_length"),
    ],
)
def test_bad_extras_are_refused(tmp_path: Path, fake: _Fake, extra, match) -> None:
    with pytest.raises(UserInputError, match=match):
        _build(tmp_path, extra=extra)
    assert fake.chains == []


def test_no_usable_gene_is_refused(tmp_path: Path, fake: _Fake) -> None:
    with pytest.raises(UserInputError, match="no gene partition"):
        _build(tmp_path, extra={"min_gene_length": 100})


def test_capabilities() -> None:
    caps = mod.AsterBuilder.capabilities
    assert not caps.auto_eligible
    assert caps.recommended_max_genomes == 5000
    assert mod.AsterBuilder.supports_partitions
    assert mod.AsterBuilder.supports_without_bootstrap
    assert {"program", "gene_tree_builder", "model", "rounds", "min_gene_length"} <= (
        caps.accepted_extras
    )
