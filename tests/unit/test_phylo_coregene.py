"""phylo --msa-source coregene: the branch, the stamp, the record, the refusals."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.annotators.base import AnnotateResult, Annotator
from repgenr.annotators.base import registry as annot_registry
from repgenr.core import annotations as cache
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CORE_GENE_FASTA, CORE_GENE_STATS
from repgenr.core.errors import UserInputError
from repgenr.core.plugins import ToolCapabilities
from repgenr.coregenes.base import CoreGeneBuilder, CoreGeneResult
from repgenr.coregenes.base import registry as cg_registry
from repgenr.stages import phylo as phylo_mod
from repgenr.stages.annotate import AnnotateStageParams
from repgenr.stages.annotate import run as annotate_run
from repgenr.stages.phylo import PhyloBuildParams, PhyloParams, phylo_build
from repgenr.stages.phylo import run as phylo_run
from repgenr.treebuilders.base import InputKind, TreeBuilder
from repgenr.treebuilders.base import registry as tb_registry

_LOG = logging.getLogger("test-phylo-coregene")


class _FakeAnnotator(Annotator):
    capabilities = ToolCapabilities(name="fakeannot")

    def preflight(self):
        return {"fakeannot": "1.0"}

    def call(self, genome, record_name, out_dir, params, logger) -> AnnotateResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        gff = out_dir / f"{record_name}.gff3"
        gff.write_text(f"##gff-version 3\n##FASTA\n>{record_name}\nACGT\n", encoding="utf-8")
        return AnnotateResult(gff=gff)


class _FakeCoreGene(CoreGeneBuilder):
    capabilities = ToolCapabilities(name="fakecg", accepted_extras=frozenset({"cgk"}))
    calls: list[list[str]] = []

    def preflight(self):
        return {"fakecg": "2.0"}

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        type(self).calls.append(sorted(g.stem for g in gffs))
        out_dir.mkdir(parents=True, exist_ok=True)
        aln = out_dir / "aln.fasta"
        aln.write_text("".join(f">{g.stem}\nACGTAC\n" for g in gffs), encoding="utf-8")
        return CoreGeneResult(alignment=aln, tables=[], core_genes=5)


class _MsaTree(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_cg")
    input_kind = InputKind.MSA_FASTA
    seen: list = []

    def preflight(self):
        return {"faketree_cg": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).seen.append(params.constant_sites)
        out_dir.mkdir(parents=True, exist_ok=True)
        names = [
            ln[1:].strip()
            for ln in Path(msa_or_genomes).read_text(encoding="utf-8").splitlines()
            if ln.startswith(">")
        ]
        tree = out_dir / "tree.nwk"
        tree.write_text("(" + ",".join(names) + ");\n", encoding="utf-8")
        return tree


class _GenomesTree(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_genomes_cg")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree_genomes_cg": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        names = [Path(g).stem for g in msa_or_genomes]
        tree.write_text("(" + ",".join(names) + ");\n", encoding="utf-8")
        return tree


@pytest.fixture
def tools(register_tool):
    _FakeCoreGene.calls = []
    _MsaTree.seen = []
    register_tool(annot_registry, "fakeannot", _FakeAnnotator)
    register_tool(cg_registry, "fakecg", _FakeCoreGene)
    register_tool(tb_registry, "faketree_cg", _MsaTree)
    register_tool(tb_registry, "faketree_genomes_cg", _GenomesTree)


def _reps(workdir: Path) -> list[Path]:
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    out = []
    for i in range(1, 4):
        p = reps / f"Fam_gen_sp_GCA_00000{i}.fasta"
        p.write_text(f">s{i}\nACGTACGT\n", encoding="utf-8")
        out.append(p)
    return out


def _params(**kw) -> PhyloParams:
    base = dict(
        treebuilder="faketree_cg", msa_source="coregene", coregene="fakecg", no_outgroup=True
    )
    base.update(kw)
    return PhyloParams(**base)


def test_coregene_source_builds_from_the_annotations(workdir: Path, tools) -> None:
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    with pytest.raises(UserInputError, match="repgenr annotate"):
        phylo_run(ctx, _params())
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    tree = phylo_run(ctx, _params(extra={"core_threshold": 0.9, "clean_mode": "moderate"}))
    assert tree.read_text().count("Fam_gen_sp_GCA") == 3
    assert (ctx.phylo_msa_dir / CORE_GENE_FASTA).is_file()
    assert (ctx.phylo_msa_dir / CORE_GENE_STATS).is_file()
    assert _MsaTree.seen == [None]  # invariant sites kept: no correction
    rec = ctx.config.stages["phylo"].params
    assert rec["msa_source"] == "coregene" and rec["coregene"] == "fakecg"
    assert rec["core_threshold"] == 0.9 and rec["clean_mode"] == "moderate"
    assert rec["annotations"]["annotator"] == "fakeannot" and rec["reference_genome"] is None
    assert rec["aligner"] is None and rec["snptyper"] is None
    # A second run reuses the alignment; a re-annotation rebuilds it.
    phylo_run(ctx, _params(extra={"core_threshold": 0.9, "clean_mode": "moderate"}))
    assert len(_FakeCoreGene.calls) == 1
    gff = cache.gff_path(ctx.annotations_dir, "Fam_gen_sp_GCA_000001")
    gff.write_text(gff.read_text() + ">extra\nAC\n", encoding="utf-8")
    records = cache.read_index(ctx.annotations_dir)
    records["Fam_gen_sp_GCA_000001"] = cache.AnnotationRecord(
        **{**records["Fam_gen_sp_GCA_000001"].__dict__, "gff_digest": cache.gff_digest(gff)}
    )
    cache.write_index(ctx.annotations_dir, records)
    phylo_run(ctx, _params(extra={"core_threshold": 0.9, "clean_mode": "moderate"}))
    assert len(_FakeCoreGene.calls) == 2
    # A changed threshold rebuilds too.
    phylo_run(ctx, _params(extra={"core_threshold": 0.8, "clean_mode": "moderate"}))
    assert len(_FakeCoreGene.calls) == 3


def test_coregene_settings_are_refused_on_other_sources_and_validated(tools) -> None:
    from repgenr.cli.param_builders import phylo_params

    with pytest.raises(UserInputError, match="--core-threshold applies only"):
        phylo_params(msa_source="aligner", extra={"core_threshold": 0.9})
    with pytest.raises(UserInputError, match="--clean-mode applies only"):
        phylo_params(msa_source="snptype", extra={"clean_mode": "strict"})
    with pytest.raises(UserInputError, match="--clean-mode"):
        phylo_params(msa_source="coregene", coregene="fakecg", extra={"clean_mode": "turbo"})
    with pytest.raises(UserInputError, match="--coregene"):
        phylo_params(msa_source="coregene", coregene="nope")
    with pytest.raises(UserInputError, match="--mask applies only"):
        phylo_params(msa_source="coregene", coregene="fakecg", extra={"mask": "gubbins"})
    assert phylo_params(msa_source="coregene", coregene="fakecg").coregene == "fakecg"
    with pytest.raises(UserInputError, match="core-threshold"):
        phylo_mod._core_threshold({"core_threshold": "abc"})


def test_stamp_key_follows_the_core_gene_settings(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")
    base = PhyloParams(msa_source="coregene")
    k1 = phylo_mod._msa_key(genomes, None, base, None, "digest1")
    assert k1 != phylo_mod._msa_key(genomes, None, base, None, "digest2")
    assert k1 != phylo_mod._msa_key(
        genomes,
        None,
        PhyloParams(msa_source="coregene", extra={"filtered": False}),
        None,
        "digest1",
    )
    assert k1 != phylo_mod._msa_key(
        genomes, None, PhyloParams(msa_source="coregene", coregene="other"), None, "digest1"
    )
    # The tree builder is not part of the key.
    assert k1 == phylo_mod._msa_key(
        genomes, None, PhyloParams(msa_source="coregene", treebuilder="fasttree"), None, "digest1"
    )


def test_phylo_build_reads_or_makes_the_annotations(tmp_path: Path, tools) -> None:
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    for i in range(1, 4):
        (gdir / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    out = tmp_path / "out"
    tree = phylo_build(
        PhyloBuildParams(genomes_dir=gdir, out_dir=out, phylo=_params(), annotator="fakeannot"),
        _LOG,
    )
    assert tree.is_file() and (out / "annotations" / "annotations.tsv").is_file()
    assert len(cache.read_index(out / "annotations")) == 3
    # Given annotations are read as they are.
    out2 = tmp_path / "out2"
    tree2 = phylo_build(
        PhyloBuildParams(
            genomes_dir=gdir, out_dir=out2, phylo=_params(), annotations_dir=out / "annotations"
        ),
        _LOG,
    )
    assert tree2.is_file() and not (out2 / "annotations").exists()
    with pytest.raises(Exception, match="not a directory"):
        phylo_build(
            PhyloBuildParams(
                genomes_dir=gdir,
                out_dir=tmp_path / "o3",
                phylo=_params(),
                annotations_dir=tmp_path / "nope",
            ),
            _LOG,
        )


def test_an_alignment_free_builder_cannot_silently_drop_the_source(
    workdir: Path, tools, caplog
) -> None:
    """auto must stay with an MSA builder, and an explicit genome-input
    builder must say that the core-gene pass did not run."""
    from repgenr.stages.phylo import _resolve_treebuilder, _wants_alignment

    assert _wants_alignment(_params()) == "--msa-source coregene"
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    # auto with a huge set: an alignment is still required, so an MSA builder wins.
    chosen = _resolve_treebuilder(_params(treebuilder="auto"), 20000, _LOG)
    assert tb_registry.create(chosen).input_kind == InputKind.MSA_FASTA
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        phylo_run(ctx, _params(treebuilder="faketree_genomes_cg"))
    assert any("--msa-source coregene" in r.getMessage() for r in caplog.records)
    assert not _FakeCoreGene.calls


def test_coregene_with_an_outgroup_roots_the_tree(workdir: Path, tools) -> None:
    _reps(workdir)
    og = workdir / "outgroup"
    og.mkdir()
    (og / "Out_gen_sp_GCA_000009.fasta").write_text(">o\nACGTACGT\n", encoding="utf-8")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009\n", encoding="utf-8")
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    tree = phylo_run(
        ctx, PhyloParams(treebuilder="faketree_cg", msa_source="coregene", coregene="fakecg")
    )
    assert "Out_gen_sp_GCA_000009" in tree.read_text(encoding="utf-8")
    # The outgroup is annotated and part of the alignment.
    assert _FakeCoreGene.calls[0] == sorted(
        [
            "Fam_gen_sp_GCA_000001",
            "Fam_gen_sp_GCA_000002",
            "Fam_gen_sp_GCA_000003",
            "Out_gen_sp_GCA_000009",
        ]
    )
    rec = ctx.config.stages["phylo"].params
    assert rec["outgroup"] == "Out_gen_sp_GCA_000009" and rec["root"] == "outgroup"


def test_an_alignment_naming_other_records_is_refused(tmp_path: Path, register_tool) -> None:
    from repgenr.core.errors import WorkdirError
    from repgenr.coregenes.base import CoreGeneParams
    from repgenr.stages.coregene import coregene_core

    class _MisnamingBuilder(_FakeCoreGene):
        def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
            out_dir.mkdir(parents=True, exist_ok=True)
            aln = out_dir / "aln.fasta"
            aln.write_text(">sample_1\nACGT\n>sample_2\nACGT\n", encoding="utf-8")
            return CoreGeneResult(alignment=aln, tables=[], core_genes=1)

    register_tool(cg_registry, "misnaming", _MisnamingBuilder)
    gffs = []
    for n in ("Fam_gen_sp_GCA_1", "Fam_gen_sp_GCA_2"):
        p = tmp_path / f"{n}.gff3"
        p.write_text("##gff-version 3\n##FASTA\n>c\nACGT\n", encoding="utf-8")
        gffs.append(p)
    with pytest.raises(WorkdirError, match="none of the genomes"):
        coregene_core(gffs, tmp_path / "msa", tmp_path / "s", CoreGeneParams(), "misnaming", _LOG)


def test_phylo_build_msa_only_and_given_msa(tmp_path: Path, tools) -> None:
    from repgenr.core.contracts import MSA_FASTA
    from repgenr.tree.alignment import ALIGNMENT_STATS

    gdir = tmp_path / "genomes"
    gdir.mkdir()
    for i in range(1, 4):
        (gdir / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    out = tmp_path / "out"
    msa = phylo_build(
        PhyloBuildParams(
            genomes_dir=gdir, out_dir=out, phylo=_params(), annotator="fakeannot", msa_only=True
        ),
        _LOG,
    )
    assert msa == out / MSA_FASTA and msa.is_file()
    # A core-gene alignment keeps its invariant sites: no sidecar travels with it.
    assert not (out / ALIGNMENT_STATS).exists()
    # With --msa the step neither annotates nor runs the builder.
    _FakeCoreGene.calls.clear()
    out2 = tmp_path / "out2"
    tree = phylo_build(
        PhyloBuildParams(genomes_dir=gdir, out_dir=out2, phylo=_params(), msa=msa), _LOG
    )
    assert tree.is_file() and not _FakeCoreGene.calls
    assert not (out2 / "annotations").exists()


def test_a_dropped_outgroup_is_an_error_not_a_warning(tmp_path: Path, register_tool) -> None:
    from repgenr.core.errors import WorkdirError
    from repgenr.coregenes.base import CoreGeneParams
    from repgenr.stages.coregene import coregene_core

    class _DroppingBuilder(_FakeCoreGene):
        def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
            out_dir.mkdir(parents=True, exist_ok=True)
            aln = out_dir / "aln.fasta"
            # The outgroup shares no core gene and is left out.
            kept = [g for g in gffs if not g.stem.startswith("Out_")]
            aln.write_text("".join(f">{g.stem}\nACGT\n" for g in kept), encoding="utf-8")
            return CoreGeneResult(alignment=aln, tables=[], core_genes=1)

    register_tool(cg_registry, "dropping", _DroppingBuilder)
    gffs = []
    for n in ("Fam_gen_sp_GCA_1", "Fam_gen_sp_GCA_2", "Out_gen_sp_GCA_9"):
        p = tmp_path / f"{n}.gff3"
        p.write_text("##gff-version 3\n##FASTA\n>c\nACGT\n", encoding="utf-8")
        gffs.append(p)
    with pytest.raises(WorkdirError, match="Out_gen_sp_GCA_9.*outgroup"):
        coregene_core(
            gffs,
            tmp_path / "msa",
            tmp_path / "s",
            CoreGeneParams(),
            "dropping",
            _LOG,
            required=("Out_gen_sp_GCA_9",),
        )
    # Without that requirement a dropped ingroup genome is only a warning.
    result, _ = coregene_core(
        gffs, tmp_path / "msa2", tmp_path / "s2", CoreGeneParams(), "dropping", _LOG
    )
    assert result.alignment.is_file()


def test_phylo_build_re_annotates_when_the_annotator_changes(tmp_path: Path, tools) -> None:
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    for i in range(1, 4):
        (gdir / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    out = tmp_path / "out"
    phylo_build(
        PhyloBuildParams(genomes_dir=gdir, out_dir=out, phylo=_params(), annotator="fakeannot"),
        _LOG,
    )
    first = cache.read_index(out / "annotations")
    assert {r.annotator for r in first.values()} == {"fakeannot"}
    # Another annotator must not leave its output indexed under the first one.
    register = annot_registry.create("fakeannot")
    annot_registry.register("fakeannot2", type(register), replace=True)
    try:
        phylo_build(
            PhyloBuildParams(
                genomes_dir=gdir, out_dir=out, phylo=_params(), annotator="fakeannot2"
            ),
            _LOG,
        )
        assert {r.annotator for r in cache.read_index(out / "annotations").values()} == {
            "fakeannot2"
        }
    finally:
        annot_registry.unregister("fakeannot2")
