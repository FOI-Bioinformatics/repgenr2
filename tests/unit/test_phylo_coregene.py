"""phylo --msa-source coregene: the branch, the stamp, the record, the refusals."""

from __future__ import annotations

import json
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
from repgenr.coregenes.base import CoreGeneBuilder, CoreGeneResult, GeneAlignment
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
    capabilities = ToolCapabilities(name="fakecg", accepted_extras=frozenset({"cgk", "threshold"}))
    calls: list[list[str]] = []
    extras: list[dict] = []

    def preflight(self):
        return {"fakecg": "2.0"}

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        type(self).calls.append(sorted(g.stem for g in gffs))
        type(self).extras.append(dict(params.extra))
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


class _PerGeneCoreGene(CoreGeneBuilder):
    """Hands over four genes of 8, 6, 4 and 10 columns, present in every genome."""

    capabilities = ToolCapabilities(name="fakecg_pergene")

    def preflight(self):
        return {"fakecg_pergene": "1.0"}

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        genes = []
        for name, length in (("g1", 8), ("g2", 6), ("g3", 4), ("g4", 10)):
            path = out_dir / f"{name}.fasta"
            path.write_text(
                "".join(
                    f">{g.stem}\n{'ACGT' * 3}"[: len(g.stem) + 2 + length] + "\n" for g in gffs
                ),
                encoding="utf-8",
            )
            genes.append(GeneAlignment(gene=name, path=path))
        return CoreGeneResult(gene_alignments=genes, core_genes=4)


@pytest.fixture
def tools(register_tool):
    register_tool(cg_registry, "fakecg_pergene", _PerGeneCoreGene)
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
    # The clustering identity reaches the builder rather than being eaten as a
    # stage-level setting, is recorded, and rebuilds the alignment on its own.
    phylo_run(ctx, _params(extra={"core_threshold": 0.8, "threshold": "0.9"}))
    assert len(_FakeCoreGene.calls) == 4 and _FakeCoreGene.extras[-1]["threshold"] == "0.9"
    assert ctx.config.stages["phylo"].params["extra"]["threshold"] == "0.9"


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
    assert k1 != phylo_mod._msa_key(
        genomes,
        None,
        PhyloParams(msa_source="coregene", extra={"threshold": "0.9"}),
        None,
        "digest1",
    )
    # The tree builder is not part of the key.
    assert k1 == phylo_mod._msa_key(
        genomes, None, PhyloParams(msa_source="coregene", treebuilder="fasttree"), None, "digest1"
    )


def test_snptype_stamp_key_follows_the_stats_version(tmp_path: Path, monkeypatch) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")
    snp = PhyloParams(msa_source="snptype", snptyper="ska2")
    aligner = PhyloParams(msa_source="aligner")
    k_snp = phylo_mod._msa_key(genomes, None, snp)
    k_aln = phylo_mod._msa_key(genomes, None, aligner)
    # 2: ska2 alignments typed before the all-sites counts are typed once more.
    assert phylo_mod._SNPTYPE_STATS_VERSION == 2
    monkeypatch.setattr(phylo_mod, "_SNPTYPE_STATS_VERSION", 1)
    assert k_snp != phylo_mod._msa_key(genomes, None, snp)
    # An aligner alignment does not carry the stats version.
    assert k_aln == phylo_mod._msa_key(genomes, None, aligner)


def test_stamp_key_follows_the_gene_selection_and_layout(tmp_path: Path, monkeypatch) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")

    def key(**extra) -> str:
        return phylo_mod._msa_key(
            genomes, None, PhyloParams(msa_source="coregene", extra=extra), None, "d"
        )

    k1 = key()
    # The defaults spelled out are the defaults.
    assert k1 == key(max_core_genes="0", gene_rank="occupancy")
    for changed in (
        {"max_core_genes": "50"},
        {"gene_rank": "informative"},
        {"min_gene_length": "300"},
        {"max_gene_gaps": "0.2"},
    ):
        assert k1 != key(**changed), changed
    monkeypatch.setattr(phylo_mod, "_COREGENE_LAYOUT", phylo_mod._COREGENE_LAYOUT + 1)
    assert k1 != key()


def test_gene_selection_is_read_from_the_stage_extras() -> None:
    from repgenr.coregenes.select import GeneSelection

    extra = {
        "max_core_genes": "50",
        "gene_rank": "length",
        "min_gene_length": "300",
        "max_gene_gaps": "0.25",
        "cgk": "1",
    }
    assert phylo_mod._gene_selection(extra) == GeneSelection(
        max_genes=50, rank="length", min_length=300, max_gaps=0.25
    )
    assert phylo_mod._gene_selection({}) == GeneSelection()
    # Stage keys are never forwarded to the builder.
    assert phylo_mod._adapter_extra(extra) == {"cgk": "1"}
    for bad in (
        {"max_core_genes": "many"},
        {"min_gene_length": "1.5"},
        {"max_gene_gaps": "x"},
        {"max_gene_gaps": "2"},
        {"gene_rank": "best"},
    ):
        with pytest.raises(UserInputError):
            phylo_mod._gene_selection(bad)


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


def test_a_dropped_genome_is_refused_before_the_tree(tmp_path: Path, register_tool) -> None:
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
    # Without that requirement a dropped ingroup genome is refused too, before
    # a tree is built that would lack it.
    with pytest.raises(WorkdirError, match="Out_gen_sp_GCA_9.*--core-threshold.*--gene-rank"):
        coregene_core(gffs, tmp_path / "msa2", tmp_path / "s2", CoreGeneParams(), "dropping", _LOG)


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


class _IgnoringCoreGene(_FakeCoreGene):
    """A builder that, like PPanGGOLiN, SCARAP and OrthoFinder, does not read
    Panaroo's cleaning mode or its filtered alignment."""

    capabilities = ToolCapabilities(
        name="fakecg_ignoring", ignored_params=frozenset({"clean_mode", "filtered"})
    )


def test_settings_a_builder_ignores_change_neither_key_nor_fingerprint(
    tmp_path: Path, register_tool
) -> None:
    from repgenr.cli.base import _stage_fingerprint

    register_tool(cg_registry, "fakecg", _FakeCoreGene)
    register_tool(cg_registry, "fakecg_ignoring", _IgnoringCoreGene)
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")

    def key(tool: str, **extra) -> str:
        return phylo_mod._msa_key(
            genomes, None, PhyloParams(msa_source="coregene", coregene=tool, extra=extra), None, "d"
        )

    def fingerprint(tool: str, **extra) -> str:
        params = PhyloParams(msa_source="coregene", coregene=tool, extra=extra)
        return _stage_fingerprint("phylo", params, {}, {})

    for probe in (key, fingerprint):
        ignoring = probe("fakecg_ignoring")
        assert ignoring == probe("fakecg_ignoring", clean_mode="moderate")
        assert ignoring == probe("fakecg_ignoring", filtered=False)
        assert ignoring != probe("fakecg_ignoring", core_threshold=0.9)
        reading = probe("fakecg")
        assert reading != probe("fakecg", clean_mode="moderate")
        assert reading != probe("fakecg", filtered=False)
    assert (
        phylo_mod._coregene_setting(
            PhyloParams(coregene="fakecg_ignoring", extra={"clean_mode": "moderate"}), "clean_mode"
        )
        is None
    )


def test_precheck_refuses_unannotated_genomes_before_any_tool_runs(workdir: Path, tools) -> None:
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    with pytest.raises(UserInputError, match="repgenr annotate"):
        phylo_mod.precheck(ctx, _params())
    # Another source, or a builder that reads genomes, does not read annotations.
    phylo_mod.precheck(ctx, _params(msa_source="aligner"))
    phylo_mod.precheck(ctx, _params(treebuilder="faketree_genomes_cg"))
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    phylo_mod.precheck(ctx, _params())
    phylo_mod.precheck(ctx, _params(treebuilder="auto"))
    assert _FakeCoreGene.calls == []


def test_a_refused_core_gene_rerun_keeps_the_finished_record(workdir: Path, tools) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app
    from repgenr.core.config import Config

    reps = _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    ctx.save_config()
    argv = [
        "phylo", "-wd", str(workdir), "--treebuilder", "faketree_cg",
        "--msa-source", "coregene", "--coregene", "fakecg", "--no-outgroup",
    ]  # fmt: skip
    done = CliRunner().invoke(app, argv)
    assert done.exit_code == 0, done.output
    assert Config.load(workdir).stages["phylo"].completed
    # A new representative without an annotation: the rerun is refused
    # before the record is marked incomplete.
    (reps[0].parent / "Fam_gen_sp_GCA_000009.fasta").write_text(">s9\nACGT\n", encoding="utf-8")
    refused = CliRunner().invoke(app, argv)
    assert refused.exit_code == 2, refused.output
    assert "repgenr annotate" in refused.output
    assert Config.load(workdir).stages["phylo"].completed
    assert len(_FakeCoreGene.calls) == 1


def test_precheck_leaves_small_sets_and_genome_builders_to_the_stage(workdir: Path, tools) -> None:
    reps = _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    # An alignment-free builder named explicitly reads genomes, not annotations.
    phylo_mod.precheck(ctx, _params(treebuilder="faketree_genomes_cg"))
    # Too few genomes for a tree: the stage names that first, not annotate.
    reps[0].unlink()
    phylo_mod.precheck(ctx, _params())


def test_settings_outside_coregene_do_not_change_default_fingerprints() -> None:
    from repgenr.cli.base import _stage_fingerprint

    default = _stage_fingerprint("phylo", PhyloParams(), {}, {})
    assert default == _stage_fingerprint("phylo", PhyloParams(extra={}), {}, {})
    assert default == _stage_fingerprint(
        "phylo", PhyloParams(extra={"clean_mode": "moderate"}), {}, {}
    )
    masked = PhyloParams(msa_source="snptype", extra={"mask": "gubbins"})
    assert _stage_fingerprint("phylo", masked, {}, {}) != default


def test_record_carries_the_gene_selection(workdir: Path, tools) -> None:
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    phylo_run(ctx, _params())
    default = ctx.config.stages["phylo"].params["gene_selection"]
    assert default == {
        "max_genes": 0,
        "rank": "occupancy",
        "min_length": 0,
        "max_gaps": None,
        "seed": 0,
    }
    # A selection is refused for a builder that returns no per-gene alignments.
    with pytest.raises(UserInputError, match="no per-gene alignments"):
        phylo_run(ctx, _params(extra={"max_core_genes": 3, "gene_rank": "informative"}))


def test_record_carries_a_non_default_selection(workdir: Path, tools) -> None:
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    extra = {"max_core_genes": 2, "gene_rank": "length", "min_gene_length": 5, "max_gene_gaps": 0.5}
    phylo_run(ctx, _params(coregene="fakecg_pergene", extra=extra))
    rec = ctx.config.stages["phylo"].params["gene_selection"]
    assert rec == {
        "max_genes": 2,
        "rank": "length",
        "min_length": 5,
        "max_gaps": 0.5,
        "seed": 0,
    }
    table = (ctx.phylo_msa_dir / "core_genes.tsv").read_text(encoding="utf-8").splitlines()
    kept = {ln.split("\t")[0] for ln in table[1:] if ln.split("\t")[1] == "true"}
    assert kept == {"g1", "g4"}


def test_fingerprint_drops_the_selection_when_no_coregene_pass_runs() -> None:
    selection = {"max_core_genes": "5", "gene_rank": "length", "min_gene_length": "9"}
    selection["max_gene_gaps"] = "0.5"
    for source in ("aligner", "snptype"):
        payload = {"msa_source": source, "extra": {**selection, "model": "GTR"}}
        assert phylo_mod.fingerprint_params(payload)["extra"] == {"model": "GTR"}
    payload = {"msa_source": "coregene", "coregene": "fakecg", "extra": dict(selection)}
    assert phylo_mod.fingerprint_params(payload)["extra"] == selection


# --partitioned: the tree builder reads the partitions.json of a core-gene alignment.


class _TwoGeneCoreGene(_FakeCoreGene):
    """Returns two per-gene alignments, so the stage writes partitions.json."""

    capabilities = ToolCapabilities(name="fakecg_genes")

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        from repgenr.coregenes.base import GeneAlignment

        type(self).calls.append(sorted(g.stem for g in gffs))
        out_dir.mkdir(parents=True, exist_ok=True)
        genes = []
        for gene, seq in (("gyrA", "ACGTAC"), ("rec/A", "TTGA")):
            path = out_dir / f"{len(genes)}.fasta"
            path.write_text("".join(f">{g.stem}\n{seq}\n" for g in gffs), encoding="utf-8")
            genes.append(GeneAlignment(gene=gene, path=path))
        return CoreGeneResult(tables=[], core_genes=2, gene_alignments=genes)


class _PartTree(_MsaTree):
    capabilities = ToolCapabilities(name="faketree_part", accepted_extras=frozenset({"asc"}))
    supports_partitions = True
    parts: list = []

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).parts.append((params.partitions, params.alphabet))
        return super().build(msa_or_genomes, out_dir, params, logger)


@pytest.fixture
def part_tools(tools, register_tool):
    _PartTree.parts = []
    register_tool(cg_registry, "fakecg_genes", _TwoGeneCoreGene)
    register_tool(tb_registry, "faketree_part", _PartTree)


def _genomes_dir(tmp_path: Path) -> Path:
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    for i in range(1, 4):
        (gdir / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    return gdir


def _step(tmp_path: Path, **phylo) -> Path:
    return phylo_build(
        PhyloBuildParams(
            genomes_dir=_genomes_dir(tmp_path),
            out_dir=tmp_path / "out",
            phylo=_params(**phylo),
            annotator="fakeannot",
        ),
        _LOG,
    )


def test_partitioned_reaches_the_tree_params(tmp_path: Path, part_tools) -> None:
    _step(
        tmp_path,
        coregene="fakecg_genes",
        treebuilder="faketree_part",
        extra={"partitioned": True},
    )
    assert _PartTree.parts == [((("gyrA", 1, 6), ("rec/A", 7, 10)), "dna")]


class _AsterLikeTree(_PartTree):
    capabilities = ToolCapabilities(
        name="faketree_aster", accepted_extras=frozenset({"min_gene_columns"})
    )
    extras: list = []

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).extras.append(dict(params.extra))
        return super().build(msa_or_genomes, out_dir, params, logger)


def test_min_gene_columns_reaches_the_builder_and_min_gene_length_does_not(
    tmp_path: Path, part_tools, register_tool
) -> None:
    # aster's column filter is min_gene_columns; min_gene_length is the
    # stage's core-gene selection filter and is never forwarded.
    _AsterLikeTree.extras = []
    register_tool(tb_registry, "faketree_aster", _AsterLikeTree)
    _step(
        tmp_path,
        coregene="fakecg_genes",
        treebuilder="faketree_aster",
        extra={"partitioned": True, "min_gene_columns": "5", "min_gene_length": "1"},
    )
    assert _AsterLikeTree.extras[0]["min_gene_columns"] == "5"
    assert "min_gene_length" not in _AsterLikeTree.extras[0]


def test_without_partitioned_the_builder_gets_no_partitions(tmp_path: Path, part_tools) -> None:
    _step(tmp_path, coregene="fakecg_genes", treebuilder="faketree_part")
    assert _PartTree.parts == [(None, "dna")]


def test_a_builder_without_partition_support_warns(tmp_path: Path, part_tools, caplog) -> None:
    with caplog.at_level(logging.WARNING):
        _step(tmp_path, coregene="fakecg_genes", extra={"partitioned": True})
    assert "faketree_cg ignores partitions" in caplog.text
    assert _MsaTree.seen == [None]


def test_partitioned_without_a_sidecar_is_refused(tmp_path: Path, part_tools) -> None:
    # The fake builder returns its own concatenation: no partitions are known.
    with pytest.raises(UserInputError, match="no partition sidecar"):
        _step(tmp_path, treebuilder="faketree_part", extra={"partitioned": True})
    assert _PartTree.parts == []


def test_partitioned_with_asc_is_refused_before_the_alignment(tmp_path: Path, part_tools) -> None:
    _FakeCoreGene.calls.clear()
    with pytest.raises(UserInputError, match="partitioned model cannot be combined"):
        _step(
            tmp_path,
            coregene="fakecg_genes",
            treebuilder="faketree_part",
            extra={"partitioned": True, "asc": "lewis"},
        )
    assert _FakeCoreGene.calls == []


def test_partitioned_is_not_part_of_the_msa_key(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")
    plain = PhyloParams(msa_source="coregene")
    part = PhyloParams(msa_source="coregene", extra={"partitioned": True})
    assert phylo_mod._msa_key(genomes, None, plain, None, "d") == phylo_mod._msa_key(
        genomes, None, part, None, "d"
    )


def test_record_carries_partitioned(workdir: Path, part_tools) -> None:
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    phylo_run(ctx, _params(coregene="fakecg_genes", treebuilder="faketree_part"))
    assert ctx.config.stages["phylo"].params["partitioned"] is False
    calls = len(_TwoGeneCoreGene.calls)
    phylo_run(
        ctx,
        _params(coregene="fakecg_genes", treebuilder="faketree_part", extra={"partitioned": True}),
    )
    assert ctx.config.stages["phylo"].params["partitioned"] is True
    # A tree-side setting: the alignment is reused.
    assert len(_TwoGeneCoreGene.calls) == calls
    assert _PartTree.parts[-1][0] == (("gyrA", 1, 6), ("rec/A", 7, 10))


def test_aster_with_partitioned_reaches_the_builder(
    tmp_path: Path, part_tools, monkeypatch
) -> None:
    """--treebuilder aster --partitioned: one gene tree per partition, then ASTER."""
    from repgenr.treebuilders import aster as aster_mod

    seen: list = []
    chains: list[list[str]] = []
    real_build = aster_mod.AsterBuilder.build

    def build(self, msa, out_dir, params, logger):
        seen.append(params.partitions)
        return real_build(self, msa, out_dir, params, logger)

    def run_chain(caps, steps, **k):
        for _, cmd in steps:
            cmd = [str(c) for c in cmd]
            chains.append(cmd)
            fasta = Path(cmd[cmd.index("-s") + 1])
            names = [
                ln[1:]
                for ln in fasta.read_text(encoding="utf-8").splitlines()
                if ln.startswith(">")
            ]
            Path(cmd[cmd.index("-pre") + 1] + ".treefile").write_text(
                f"({names[0]},{names[1]},({','.join(names[2:])})/1:1);\n",
                encoding="utf-8",
            )

    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        first = Path(cmd[cmd.index("-i") + 1]).read_text(encoding="utf-8").splitlines()[0]
        Path(cmd[cmd.index("-o") + 1]).write_text(first + "\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(aster_mod.AsterBuilder, "build", build)
    monkeypatch.setattr(aster_mod, "run_chain", run_chain)
    monkeypatch.setattr(aster_mod, "run_tool", run_tool)
    monkeypatch.setattr(aster_mod, "preflight", lambda caps: {caps.name: "1.0"})
    monkeypatch.setattr(aster_mod, "_program_version", lambda caps, program: "")
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    for i in range(1, 5):
        (gdir / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    out = phylo_build(
        PhyloBuildParams(
            genomes_dir=gdir,
            out_dir=tmp_path / "out",
            phylo=_params(
                coregene="fakecg_genes", treebuilder="aster", extra={"partitioned": True}
            ),
            annotator="fakeannot",
        ),
        _LOG,
    )
    assert seen == [(("gyrA", 1, 6), ("rec/A", 7, 10))]
    assert len(chains) == 2
    tree_dir = tmp_path / "out" / "tree"
    assert len((tree_dir / "gene_trees.nwk").read_text(encoding="utf-8").splitlines()) == 2
    assert out.is_file()


def test_aster_wastral_with_midpoint_root_is_refused_before_any_work(
    tmp_path: Path, part_tools, monkeypatch
) -> None:
    """wastral writes no terminal lengths: --root midpoint fails up front."""
    from repgenr.treebuilders import aster as aster_mod

    chains: list = []
    monkeypatch.setattr(aster_mod, "run_chain", lambda *a, **k: chains.append(a))
    monkeypatch.setattr(aster_mod, "preflight", lambda caps: {caps.name: "1.0"})
    _TwoGeneCoreGene.calls.clear()
    with pytest.raises(UserInputError, match="--root outgroup"):
        _step(
            tmp_path,
            coregene="fakecg_genes",
            treebuilder="aster",
            root="midpoint",
            extra={"partitioned": True},
        )
    assert chains == []
    assert _TwoGeneCoreGene.calls == []


# --coregene-alphabet: protein core-gene alignments and the tree model.


def test_the_alphabet_is_part_of_the_msa_key(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")

    def key(**extra) -> str:
        return phylo_mod._msa_key(
            genomes, None, PhyloParams(msa_source="coregene", extra=extra), None, "d"
        )

    assert key() == key(coregene_alphabet="dna")
    assert key() != key(coregene_alphabet="protein")


def test_a_dna_key_has_no_alphabet_field(tmp_path: Path, monkeypatch) -> None:
    """A nucleotide core-gene alignment keeps the key it had before the
    alphabet existed, so it is reused rather than rebuilt."""
    genomes = [tmp_path / "a.fasta"]
    genomes[0].write_text(">x\nACGT\n", encoding="utf-8")
    payloads: list[dict] = []
    real = phylo_mod.json.dumps

    def capture(obj, **kw):
        payloads.append(obj)
        return real(obj, **kw)

    monkeypatch.setattr(phylo_mod.json, "dumps", capture)
    for alphabet in ("dna", "protein"):
        params = PhyloParams(msa_source="coregene", extra={"coregene_alphabet": alphabet})
        phylo_mod._msa_key(genomes, None, params, None, "d")
    dna, protein = (p["coregene"] for p in payloads)
    assert "alphabet" not in dna and protein["alphabet"] == "protein"


def test_the_sidecar_alphabet_reaches_the_tree_params(tmp_path: Path, part_tools) -> None:
    gdir = _genomes_dir(tmp_path)
    protein = _params(
        coregene="fakecg_genes",
        treebuilder="faketree_part",
        extra={"coregene_alphabet": "protein"},
    )
    msa = phylo_build(
        PhyloBuildParams(
            genomes_dir=gdir,
            out_dir=tmp_path / "msa",
            phylo=protein,
            annotator="fakeannot",
            msa_only=True,
        ),
        _LOG,
    )
    sidecar = json.loads((msa.parent / "partitions.json").read_text(encoding="utf-8"))
    assert sidecar["alphabet"] == "protein"
    # The tree half (as Nextflow runs it) is not told the alphabet: the
    # sidecar beside the alignment decides it.
    phylo_build(
        PhyloBuildParams(
            genomes_dir=gdir,
            out_dir=tmp_path / "tree",
            phylo=_params(coregene="fakecg_genes", treebuilder="faketree_part"),
            msa=msa,
        ),
        _LOG,
    )
    assert _PartTree.parts == [(None, "protein")]


def test_without_a_sidecar_the_alphabet_comes_from_the_params(tmp_path: Path, part_tools) -> None:
    # The fake builder returns its own concatenation: no partitions.json.
    _step(tmp_path, treebuilder="faketree_part", extra={"coregene_alphabet": "protein"})
    assert _PartTree.parts == [(None, "protein")]


def test_record_carries_the_alphabet(workdir: Path, part_tools) -> None:
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    phylo_run(ctx, _params(coregene="fakecg_genes", treebuilder="faketree_part"))
    assert ctx.config.stages["phylo"].params["coregene_alphabet"] == "dna"
    calls = len(_TwoGeneCoreGene.calls)
    phylo_run(
        ctx,
        _params(
            coregene="fakecg_genes",
            treebuilder="faketree_part",
            extra={"coregene_alphabet": "protein"},
        ),
    )
    assert ctx.config.stages["phylo"].params["coregene_alphabet"] == "protein"
    # The alignment is rebuilt for the other alphabet.
    assert len(_TwoGeneCoreGene.calls) == calls + 1
    stats = json.loads((ctx.phylo_msa_dir / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["alphabet"] == "protein"


def test_the_alphabet_is_refused_outside_coregene_and_validated(tools) -> None:
    from repgenr.cli.param_builders import phylo_params

    with pytest.raises(UserInputError, match="--coregene-alphabet applies only"):
        phylo_params(msa_source="aligner", extra={"coregene_alphabet": "protein"})
    with pytest.raises(UserInputError, match="--coregene-alphabet"):
        phylo_params(msa_source="coregene", coregene="fakecg", extra={"coregene_alphabet": "rna"})
    params = phylo_params(
        msa_source="coregene", coregene="fakecg", extra={"coregene_alphabet": "protein"}
    )
    assert params.extra["coregene_alphabet"] == "protein"
    with pytest.raises(UserInputError, match="--coregene-alphabet"):
        phylo_mod._coregene_alphabet({"coregene_alphabet": "rna"})


def test_fingerprint_drops_the_alphabet_when_no_coregene_pass_runs() -> None:
    for source in ("aligner", "snptype"):
        payload = {"msa_source": source, "extra": {"coregene_alphabet": "protein"}}
        assert phylo_mod.fingerprint_params(payload)["extra"] == {}
    payload = {
        "msa_source": "coregene",
        "coregene": "ppanggolin",
        "extra": {"coregene_alphabet": "protein"},
    }
    assert phylo_mod.fingerprint_params(payload)["extra"] == {"coregene_alphabet": "protein"}


def test_a_stale_aligner_reference_note_is_not_recorded_on_the_coregene_route(
    workdir: Path, tools
) -> None:
    """An align/reference.txt left by an earlier aligner run names no
    reference of a core-gene tree."""
    from repgenr.tree.medoid import ReferenceChoice, note_reference

    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    note_reference(ctx.align_dir, ReferenceChoice(Path("Fam_gen_sp_GCA_000001.fasta"), "medoid"))
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    phylo_run(ctx, _params())
    rec = ctx.config.stages["phylo"].params
    assert rec["reference_genome"] is None and rec["reference_rule"] is None


def test_reference_is_in_the_msa_key_only_where_the_source_uses_one(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")

    def key(**kw) -> str:
        return phylo_mod._msa_key(genomes, None, PhyloParams(**kw), None, "digest1")

    # Routes that use no reference: the key ignores --reference.
    for source in (
        dict(msa_source="coregene"),
        dict(msa_source="snptype", snptyper="ska2"),
    ):
        assert key(**source) == key(**source, reference="a")
        assert key(**source, reference="a") == key(**source, reference="b")
    # Routes that project onto one: the key follows it.
    for source in (
        dict(msa_source="snptype", snptyper="simple"),
        dict(msa_source="aligner"),
    ):
        assert key(**source, reference="a") != key(**source, reference="b")
        assert key(**source) != key(**source, reference="a")


# --gene-trim: per-gene trimming with ClipKIT on the core-gene route.


def _fake_clipkit(monkeypatch) -> list[str]:
    """clipkit drops the last column of each gene; returns the trimmed file stems."""
    from repgenr.coregenes import trim
    from repgenr.stages import coregene as stage

    trimmed: list[str] = []

    def run_chain(caps, steps, **k):
        for _, argv in steps:
            src, dst = Path(argv[1]), Path(argv[argv.index("-o") + 1])
            trimmed.append(dst.stem)
            lines = src.read_text(encoding="utf-8").splitlines()
            out = [ln if ln.startswith(">") else ln[:-1] for ln in lines]
            dst.write_text("\n".join(out) + "\n", encoding="utf-8")

    monkeypatch.setattr(trim, "run_chain", run_chain)
    monkeypatch.setattr(stage, "preflight_clipkit", lambda: {"clipkit": "2.14.0"})
    return trimmed


def test_gene_trim_is_part_of_the_msa_key_only_when_on(tmp_path: Path, monkeypatch) -> None:
    genomes = [tmp_path / "a.fasta"]
    genomes[0].write_text(">x\nACGT\n", encoding="utf-8")

    def key(**extra) -> str:
        return phylo_mod._msa_key(
            genomes, None, PhyloParams(msa_source="coregene", extra=extra), None, "d"
        )

    assert key() == key(gene_trim="none")
    assert key() != key(gene_trim="clipkit")
    payloads: list[dict] = []
    real = phylo_mod.json.dumps

    def capture(obj, **kw):
        payloads.append(obj)
        return real(obj, **kw)

    monkeypatch.setattr(phylo_mod.json, "dumps", capture)
    for value, alphabet in (("none", "dna"), ("clipkit", "dna"), ("clipkit", "protein")):
        extra = {"gene_trim": value, "coregene_alphabet": alphabet}
        phylo_mod._msa_key(
            genomes, None, PhyloParams(msa_source="coregene", extra=extra), None, "d"
        )
    off, dna, protein = (p["coregene"] for p in payloads)
    # The codon-preserving trim is in the dna key, so an alignment trimmed
    # column by column (before -co) is rebuilt.
    assert "gene_trim" not in off and dna["gene_trim"] == "clipkit:smart-gap:codon"
    assert protein["gene_trim"] == "clipkit:smart-gap"


def test_record_carries_gene_trim_and_the_alignment_is_trimmed(
    workdir: Path, part_tools, monkeypatch
) -> None:
    trimmed = _fake_clipkit(monkeypatch)
    _reps(workdir)
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    annotate_run(ctx, AnnotateStageParams(tool="fakeannot"))
    phylo_run(ctx, _params(coregene="fakecg_genes", treebuilder="faketree_part"))
    assert ctx.config.stages["phylo"].params["gene_trim"] == "none"
    assert trimmed == []
    calls = len(_TwoGeneCoreGene.calls)
    phylo_run(
        ctx,
        _params(
            coregene="fakecg_genes", treebuilder="faketree_part", extra={"gene_trim": "clipkit"}
        ),
    )
    assert ctx.config.stages["phylo"].params["gene_trim"] == "clipkit"
    # The alignment is rebuilt, from the trimmed genes (6 + 4 columns, less one each).
    assert len(_TwoGeneCoreGene.calls) == calls + 1
    assert sorted(trimmed) == ["gyrA", "rec_A"]
    stats = json.loads((ctx.phylo_msa_dir / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["trim"] == "clipkit" and stats["trimmed_columns_removed"] == 2
    assert stats["alignment_length"] == 8


def test_gene_trim_is_refused_outside_coregene_and_validated(tools) -> None:
    from repgenr.cli.param_builders import phylo_params

    with pytest.raises(UserInputError, match="--gene-trim applies only"):
        phylo_params(msa_source="aligner", extra={"gene_trim": "clipkit"})
    with pytest.raises(UserInputError, match="--gene-trim"):
        phylo_params(msa_source="coregene", coregene="fakecg", extra={"gene_trim": "trimal"})
    params = phylo_params(msa_source="coregene", coregene="fakecg", extra={"gene_trim": "clipkit"})
    assert params.extra["gene_trim"] == "clipkit"
    with pytest.raises(UserInputError, match="--gene-trim"):
        phylo_mod._gene_trim({"gene_trim": "trimal"})


def test_fingerprint_drops_gene_trim_when_no_coregene_pass_runs() -> None:
    payload = {"msa_source": "aligner", "extra": {"gene_trim": "clipkit"}}
    assert phylo_mod.fingerprint_params(payload)["extra"] == {}


def test_tree_builder_extras_are_not_part_of_the_coregene_msa_key(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")

    def key(coregene: str = "panaroo", **extra) -> str:
        params = PhyloParams(msa_source="coregene", coregene=coregene, extra=extra)
        return phylo_mod._msa_key(genomes, None, params, None, "d")

    base = key()
    # Tree builder settings change the tree, not the alignment.
    assert key(alrt="1000") == base
    assert key(model="GTR+G") == base
    assert key(gamma="true") == base
    assert key(alrt="1000", model="GTR+G", gamma="true") == base
    # A setting the core-gene builder reads rebuilds the alignment.
    pp = key("ppanggolin")
    assert key("ppanggolin", identity="0.9") != pp
    assert key("ppanggolin", identity="0.9", alrt="1000") == key("ppanggolin", identity="0.9")


def test_tree_builder_extras_are_not_part_of_the_snptype_msa_key(tmp_path: Path) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    for g in genomes:
        g.write_text(">x\nACGT\n", encoding="utf-8")

    def key(**extra) -> str:
        params = PhyloParams(msa_source="snptype", snptyper="ska2", extra=extra)
        return phylo_mod._msa_key(genomes, None, params)

    base = key()
    assert key(model="GTR+G") == base
    assert key(ksize="21") != base
    # A masker extra counts only when that masker runs.
    assert key(gubbins_force="true") == base
    masked = key(mask="gubbins")
    assert key(mask="gubbins", gubbins_force="true") != masked
    assert key(mask="gubbins", model="GTR+G") == masked
