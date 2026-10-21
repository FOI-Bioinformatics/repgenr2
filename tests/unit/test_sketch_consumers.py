"""The consumers of the workdir genome sketches, with a fake sourmash.

The sketch writer is the fake of tests/conftest.py (FakeSourmash): each
``sketches/<name>.sig.zip`` holds the signature name and the genome path. The
consumers' own sourmash calls go through :class:`FakeConsumer`, which answers
``sketch``, ``compare``, ``scripts manysketch`` and ``scripts pairwise`` from
the files it is given, so a test can see which signature files and which
``-k`` each command was given.
"""

from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path

import pytest

from repgenr.core import sketches
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import record_name
from repgenr.core.manifest import GenomeRecord
from repgenr.core.sketches import (
    adapter_sketches,
    contract_mismatch,
    directory_provider,
    resolve_sketches,
    sketch_genomes,
    sketch_targets,
)

_LOG = logging.getLogger("test")

# The two genomes the fake comparisons call close (Jaccard 0.9, containment 0.95).
_CLOSE = frozenset({"Fam_Gen_sp_GCA_000000.1", "Fam_Gen_sp_GCA_000001.1"})


def _workdir(root: Path, n: int = 3, *, outgroup: bool = True) -> WorkdirContext:
    ctx = WorkdirContext(root, create=True)
    ctx.genomes_dir.mkdir()
    records = []
    for i in range(n):
        name = f"Fam_Gen_sp_GCA_{i:06d}.1.fasta"
        (ctx.genomes_dir / name).write_text(f">g{i}\n{'ACGT' * (i + 5)}\n", encoding="utf-8")
        records.append(GenomeRecord(f"GCA_{i:06d}.1", name, "local", "Fam", "Gen", "sp"))
    if outgroup:
        ctx.outgroup_dir.mkdir()
        name = "Fam_Out_og_GCF_000009.1.fasta"
        (ctx.outgroup_dir / name).write_text(">o\nTTTTGGGG\n", encoding="utf-8")
        records.append(GenomeRecord("GCF_000009.1", name, "local", "Fam", "Out", "og", True))
        (ctx.workdir / "outgroup_accession.txt").write_text("GCF_000009.1\n", encoding="utf-8")
    ctx.manifest.replace_genomes(records)
    return ctx


def _genomes(ctx: WorkdirContext) -> list[Path]:
    return sorted(ctx.genomes_dir.glob("*.fasta"))


def _jaccard(a: str, b: str) -> float:
    if a == b:
        return 1.0
    return 0.9 if {a, b} == _CLOSE else 0.01


class FakeConsumer:
    """sourmash as a dereplicator or tree builder calls it."""

    def __init__(self) -> None:
        self.calls: list[list[str]] = []
        # The files each --from-file list named, read at call time (a step
        # removes its scratch directory when it finishes).
        self.lists: dict[tuple[str, ...], list[str]] = {}

    @staticmethod
    def _names(path: Path) -> list[str]:
        """Signature names in a file: a workdir sketch, a .sig, or a manysketch zip."""
        lines = path.read_text(encoding="utf-8").splitlines()
        if path.name.endswith((".sig.zip", ".sig")):
            return [lines[0]]
        return lines

    def run_tool(self, caps, argv, **kwargs) -> int:  # noqa: ANN001
        argv = [str(a) for a in argv]
        self.calls.append(argv)
        if "--from-file" in argv:
            fofn = Path(argv[argv.index("--from-file") + 1])
            self.lists[tuple(argv)] = fofn.read_text(encoding="utf-8").split()
        if "--help" in argv:
            return 0
        if argv[1] == "sketch":
            fofn = Path(argv[argv.index("--from-file") + 1])
            outdir = Path(argv[argv.index("--outdir") + 1])
            for line in fofn.read_text(encoding="utf-8").split():
                # sourmash sketch leaves the name empty: the label is the path.
                (outdir / f"{Path(line).name}.sig").write_text(f"{line}\n", encoding="utf-8")
        elif argv[1] == "compare":
            fofn = Path(argv[argv.index("--from-file") + 1])
            labels = [
                n for f in fofn.read_text(encoding="utf-8").split() for n in self._names(Path(f))
            ]
            labels.reverse()  # sourmash compare does not keep its input order
            with open(argv[argv.index("--csv") + 1], "w", encoding="utf-8", newline="") as fo:
                w = csv.writer(fo)
                w.writerow(labels)
                for a in labels:
                    w.writerow([_jaccard(record_name(a), record_name(b)) for b in labels])
        elif argv[1:3] == ["scripts", "manysketch"]:
            rows = list(csv.reader(open(argv[3], encoding="utf-8")))[1:]
            Path(argv[argv.index("-o") + 1]).write_text(
                "".join(f"{r[0]}\n" for r in rows), encoding="utf-8"
            )
        elif argv[1:3] == ["scripts", "pairwise"]:
            source = Path(argv[3])
            files = (
                [Path(x) for x in source.read_text(encoding="utf-8").split()]
                if source.suffix == ".txt"
                else [source]
            )
            names = [n for f in files for n in self._names(f)]
            with open(argv[argv.index("-o") + 1], "w", encoding="utf-8", newline="") as fo:
                w = csv.writer(fo)
                w.writerow(["query_name", "match_name", "average_containment"])
                for i, a in enumerate(names):
                    for b in names[i + 1 :]:
                        if _jaccard(record_name(a), record_name(b)) > 0.5:
                            w.writerow([a, b, "0.95"])
        return 0

    def sub(self, *words: str) -> list[list[str]]:
        return [c for c in self.calls if c[1 : 1 + len(words)] == list(words)]

    def listed(self, argv: list[str]) -> list[str]:
        return self.lists[tuple(argv)]


@pytest.fixture
def consumer(monkeypatch) -> FakeConsumer:
    import repgenr.dereplicators.sourmash as derep_sm
    import repgenr.treebuilders.sourmash as tree_sm

    fake = FakeConsumer()
    monkeypatch.setattr(derep_sm, "run_tool", fake.run_tool)
    monkeypatch.setattr(tree_sm, "run_tool", fake.run_tool)
    monkeypatch.setattr(derep_sm.SourmashDereplicator, "preflight", lambda self: {"sourmash": "4"})
    monkeypatch.setattr(tree_sm.SourmashBuilder, "preflight", lambda self: {"sourmash": "4"})
    monkeypatch.setattr(derep_sm, "_branchwater_available", lambda caps, logger: False)
    return fake


@pytest.fixture
def branchwater(monkeypatch) -> None:
    import repgenr.dereplicators.sourmash as derep_sm

    monkeypatch.setattr(derep_sm, "_branchwater_available", lambda caps, logger: True)


# --- the contract parameters and the provider -------------------------------------


def test_contract_serves_its_three_ksizes_at_scaled_1000() -> None:
    for k in (21, 31, 51):
        assert contract_mismatch(k, 1000) is None
    assert "k=41" in (contract_mismatch(41, 1000) or "")
    assert "scaled=2000" in (contract_mismatch(31, 2000) or "")


def test_adapter_sketches_falls_back_outside_the_contract(caplog) -> None:
    from repgenr.dereplicators.sourmash import SourmashDereplicator

    asked: list[list[Path]] = []

    def provider(genomes):  # noqa: ANN001, ANN202
        asked.append(list(genomes))
        return {g: Path(f"{g}.sig.zip") for g in genomes}

    adapter = SourmashDereplicator()
    genomes = [Path("a.fasta")]
    with caplog.at_level(logging.INFO, logger="test"):
        assert adapter_sketches(adapter, {"ksize": "41"}, genomes, provider, _LOG, "x") is None
        assert adapter_sketches(adapter, {"scaled": 2000}, genomes, provider, _LOG, "x") is None
    assert asked == []
    assert "x: genome sketches not used (k=41 is not one of the sketched" in caplog.text
    assert "scaled=2000 differs from the scaled=1000" in caplog.text
    # An adapter without sketch_request, or no provider: the adapter sketches itself.
    assert adapter_sketches(object(), {}, genomes, provider, _LOG, "x") is None
    assert adapter_sketches(adapter, {}, genomes, None, _LOG, "x") is None
    got = adapter_sketches(adapter, {"ksize": "21"}, genomes, provider, _LOG, "x")
    assert got == {Path("a.fasta"): Path("a.fasta.sig.zip")}


# --- resolve_sketches ---------------------------------------------------------------


def test_resolve_reuses_current_writes_missing_and_rewrites_stale(
    tmp_path, fake_sourmash, caplog
) -> None:
    ctx = _workdir(tmp_path / "wd")
    g0, g1, g2 = _genomes(ctx)
    targets = {t.path: t for t in sketch_targets(ctx)}
    sketch_genomes(ctx, [targets[g0], targets[g1]], 2, _LOG)
    g1.write_text(">g1\nGGGGCCCCAAAA\n", encoding="utf-8")  # its sketch is now stale
    stray = ctx.genomes_dir / "stray.fasta"  # no manifest row: not part of the set
    stray.write_text(">s\nACGT\n", encoding="utf-8")
    outgroup = next(ctx.outgroup_dir.glob("*.fasta"))
    fake_sourmash.calls.clear()

    with caplog.at_level(logging.INFO, logger="test"):
        got = resolve_sketches(ctx, [g0, g1, g2, outgroup, stray], _LOG, 2, consumer="test")
    assert set(got) == {g0, g1, g2, outgroup}
    assert got[g0] == ctx.workdir / "sketches" / "Fam_Gen_sp_GCA_000000.1.sig.zip"
    assert sorted(fake_sourmash.calls) == sorted(record_name(p) for p in (g1, g2, outgroup))
    assert "test: sketches: 1 reused, 3 written, 1 sketched by the tool" in caplog.text
    # Recorded: a second consumer reuses all of them.
    fake_sourmash.calls.clear()
    assert resolve_sketches(ctx, [g0, g1, g2, outgroup], _LOG, 2, consumer="test") == got
    assert fake_sourmash.calls == []


def test_resolve_hashes_a_genome_again_only_when_its_size_or_mtime_changes(
    tmp_path, fake_sourmash, monkeypatch
) -> None:
    import os

    ctx = _workdir(tmp_path / "wd", outgroup=False)
    genomes = _genomes(ctx)
    resolve_sketches(ctx, genomes, _LOG, 2, consumer="t")
    assert (ctx.workdir / "sketches" / sketches.DIGESTS_JSON).is_file()
    hashed: list[str] = []
    real = sketches.file_digest

    def counting(path):  # noqa: ANN001, ANN202
        hashed.append(Path(path).name)
        return real(path)

    monkeypatch.setattr(sketches, "file_digest", counting)
    fake_sourmash.calls.clear()
    assert len(resolve_sketches(ctx, genomes, _LOG, 2, consumer="t")) == 3
    assert hashed == [] and fake_sourmash.calls == []
    # Same size, new content and a new mtime: hashed again, found stale, sketched.
    g0 = genomes[0]
    g0.write_text(g0.read_text().replace("ACGT", "TGCA", 1), encoding="utf-8")
    st = g0.stat()
    os.utime(g0, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))
    resolve_sketches(ctx, genomes, _LOG, 2, consumer="t")
    assert g0.name in hashed and fake_sourmash.calls == [record_name(g0)]
    # An edit that keeps the size and the mtime is not seen; --force hashes
    # again and corrects the cache, so the next consumer agrees with the record.
    g1 = genomes[1]
    st = g1.stat()
    g1.write_text(g1.read_text().replace("ACGT", "TGCA", 1), encoding="utf-8")
    os.utime(g1, ns=(st.st_atime_ns, st.st_mtime_ns))
    fake_sourmash.calls.clear()
    resolve_sketches(ctx, genomes, _LOG, 2, consumer="t")
    assert fake_sourmash.calls == []
    sketch_genomes(ctx, sketch_targets(ctx), 2, _LOG, force=True)
    assert ctx.manifest.sketch_records()["GCA_000001.1"].digest == real(g1)
    fake_sourmash.calls.clear()
    resolve_sketches(ctx, genomes, _LOG, 2, consumer="t")
    assert fake_sourmash.calls == []
    # The sketches and the cache go together once no genome is left.
    for g in genomes:
        g.unlink()
    ctx.manifest.replace_genomes([])
    sketches.remove_stale(ctx)
    assert not (ctx.workdir / "sketches").exists()


def test_resolve_without_write_leaves_the_workdir_unchanged(tmp_path, fake_sourmash) -> None:
    """write=False reads the current sketches and writes nothing, the digest
    cache included (phylo-advise and the reference choice only read)."""
    ctx = _workdir(tmp_path / "wd", outgroup=False)
    genomes = _genomes(ctx)
    resolve_sketches(ctx, genomes, _LOG, 1, consumer="t")
    cache = ctx.workdir / "sketches" / sketches.DIGESTS_JSON
    before = {p: p.stat().st_mtime_ns for p in (ctx.workdir / "sketches").iterdir()}
    cache.unlink()
    got = resolve_sketches(ctx, genomes, _LOG, 1, consumer="t", write=False)
    assert set(got) == set(genomes)
    assert not cache.exists()
    after = {p: p.stat().st_mtime_ns for p in (ctx.workdir / "sketches").iterdir()}
    assert after == {p: m for p, m in before.items() if p != cache}


def test_resolve_without_sourmash_returns_only_current_sketches(
    tmp_path, fake_sourmash, monkeypatch
) -> None:
    ctx = _workdir(tmp_path / "wd", outgroup=False)
    g0, g1, _g2 = _genomes(ctx)
    targets = {t.path: t for t in sketch_targets(ctx)}
    sketch_genomes(ctx, [targets[g0]], 1, _LOG)
    monkeypatch.setattr(sketches, "tool_available", lambda caps: False)
    assert resolve_sketches(ctx, [g0, g1], _LOG, 1, consumer="t") == {
        g0: sketches.sketch_path(ctx, g0)
    }


def test_resolve_matches_representatives_by_file_name(tmp_path, fake_sourmash) -> None:
    ctx = _workdir(tmp_path / "wd", outgroup=False)
    sketch_genomes(ctx, sketch_targets(ctx), 2, _LOG)
    ctx.representatives_dir.mkdir(parents=True)
    rep = ctx.representatives_dir / _genomes(ctx)[0].name
    shutil.copy2(_genomes(ctx)[0], rep)
    fake_sourmash.calls.clear()
    got = resolve_sketches(ctx, [rep], _LOG, 1, consumer="t")
    assert got == {rep: sketches.sketch_path(ctx, rep)} and fake_sourmash.calls == []


def test_a_failed_sketch_leaves_that_genome_to_the_tool(tmp_path, fake_sourmash, caplog) -> None:
    ctx = _workdir(tmp_path / "wd", outgroup=False)
    g0, g1, g2 = _genomes(ctx)
    fake_sourmash.fail = {record_name(g1)}
    with caplog.at_level(logging.WARNING, logger="test"):
        got = resolve_sketches(ctx, [g0, g1, g2], _LOG, 1, consumer="t")
    assert g1 not in got and len(got) in (0, 1, 2)  # no genome after the failure starts
    assert "writing genome sketches stopped" in caplog.text


def test_directory_provider_matches_by_record_name(tmp_path) -> None:
    directory = tmp_path / "sketches"
    directory.mkdir()
    (directory / "a.sig.zip").write_text("a\n", encoding="utf-8")
    genomes = [tmp_path / "a.fasta.gz", tmp_path / "b.fasta"]
    got = directory_provider(directory, _LOG, "step")(genomes)
    assert got == {genomes[0]: directory / "a.sig.zip"}


# --- dereplicate ---------------------------------------------------------------------


def _derep(ctx: WorkdirContext, **extra):  # noqa: ANN202
    from repgenr.stages.dereplicate import DereplicateParams, run

    return run(
        ctx,
        DereplicateParams(
            tool="sourmash", keeper="tool", allow_incomplete=True, threads=2, extra=extra or None
        ),
    )


def _close_cluster(result) -> None:  # noqa: ANN001
    members = sorted([rep, *m] for rep, m in result.clusters.items() if len(m) == 1)
    assert members == [["Fam_Gen_sp_GCA_000000.1.fasta", "Fam_Gen_sp_GCA_000001.1.fasta"]]
    assert len(result.representatives) == 2


def test_dereplicate_dense_compares_the_workdir_sketches(tmp_path, fake_sourmash, consumer) -> None:
    ctx = _workdir(tmp_path / "wd")
    result = _derep(ctx)
    _close_cluster(result)
    assert consumer.sub("sketch") == []  # nothing sketched by the adapter
    (compare,) = consumer.sub("compare")
    assert compare[compare.index("-k") + 1] == "31"
    listed = consumer.listed(compare)
    assert listed == [str(sketches.sketch_path(ctx, g)) for g in _genomes(ctx)]
    # The three genomes were sketched into the workdir; the outgroup is not compared.
    assert sorted(fake_sourmash.calls) == sorted(record_name(g) for g in _genomes(ctx))
    # A second run reuses them.
    fake_sourmash.calls.clear()
    _close_cluster(_derep(ctx))
    assert fake_sourmash.calls == []


def test_dereplicate_sparse_reads_a_path_list_of_the_sketches(
    tmp_path, fake_sourmash, consumer, branchwater
) -> None:
    ctx = _workdir(tmp_path / "wd")
    _close_cluster(_derep(ctx, ksize="21"))
    assert consumer.sub("scripts", "manysketch") == []
    (pairwise,) = consumer.sub("scripts", "pairwise")
    assert pairwise[pairwise.index("-k") + 1] == "21"
    source = Path(pairwise[3])
    assert source.name == "signatures.txt"
    assert source.read_text(encoding="utf-8").split() == [
        str(sketches.sketch_path(ctx, g)) for g in _genomes(ctx)
    ]


def test_dereplicate_sparse_sketches_a_genome_outside_the_set(
    tmp_path, fake_sourmash, consumer, branchwater
) -> None:
    ctx = _workdir(tmp_path / "wd")
    (ctx.genomes_dir / "zz.fasta").write_text(">z\nACGT\n", encoding="utf-8")
    result = _derep(ctx)
    assert result.genome_status["zz.fasta"] == "representative"
    (many,) = consumer.sub("scripts", "manysketch")
    assert [r[0] for r in csv.reader(open(many[3], encoding="utf-8"))][1:] == ["zz"]
    (pairwise,) = consumer.sub("scripts", "pairwise")
    listed = Path(pairwise[3]).read_text(encoding="utf-8").split()
    assert listed[:3] == [str(sketches.sketch_path(ctx, g)) for g in _genomes(ctx)[:3]]
    assert listed[3].endswith("signatures.zip")


def test_dereplicate_falls_back_to_its_own_sketch_on_other_scaled(
    tmp_path, fake_sourmash, consumer
) -> None:
    ctx = _workdir(tmp_path / "wd")
    _close_cluster(_derep(ctx, scaled="2000"))
    (sketch,) = consumer.sub("sketch")
    assert sketch[sketch.index("-p") + 1] == "k=31,scaled=2000"
    assert fake_sourmash.calls == []  # the workdir sketches are not touched
    assert not (ctx.workdir / "sketches").exists()


def test_dereplicate_target_reps_reuses_the_sketches_each_step(
    tmp_path, fake_sourmash, consumer
) -> None:
    from repgenr.stages.dereplicate import DereplicateParams, run

    ctx = _workdir(tmp_path / "wd")
    run(
        ctx, DereplicateParams(tool="sourmash", keeper="tool", allow_incomplete=True, target_reps=3)
    )
    assert consumer.sub("sketch") == []
    assert len(consumer.sub("compare")) >= 2
    assert len(fake_sourmash.calls) == 3  # sketched once, before the search


# --- the stateless steps (Nextflow) ----------------------------------------------------


def test_chunk_and_merge_read_a_sketches_dir_and_work_without_one(
    tmp_path, fake_sourmash, consumer
) -> None:
    from repgenr.stages.derep_steps import (
        ChunkParams,
        MergeParams,
        dereplicate_chunk,
        dereplicate_merge,
    )

    ctx = _workdir(tmp_path / "wd", outgroup=False)
    sketch_genomes(ctx, sketch_targets(ctx), 2, _LOG)
    staged = tmp_path / "sketches"
    shutil.copytree(ctx.workdir / "sketches", staged)
    genomes = _genomes(ctx)

    chunk = dereplicate_chunk(
        ChunkParams(tool="sourmash", genomes=genomes, out_dir=tmp_path / "c0", sketches_dir=staged),
        _LOG,
    )
    _close_cluster(chunk)
    assert consumer.sub("sketch") == []
    assert consumer.listed(consumer.sub("compare")[0]) == [
        str(staged / (record_name(g) + ".sig.zip")) for g in genomes
    ]

    consumer.calls.clear()
    merged = dereplicate_merge(
        MergeParams(tool="sourmash", chunk_dirs=[tmp_path / "c0"], out_dir=tmp_path / "m",
                    sketches_dir=staged),
        _LOG,
    )  # fmt: skip
    assert len(merged.representatives) == 2
    assert consumer.sub("sketch") == []
    assert all(p.startswith(str(staged)) for p in consumer.listed(consumer.sub("compare")[0]))

    # Without a sketches directory the adapter sketches, as before.
    consumer.calls.clear()
    _close_cluster(
        dereplicate_chunk(
            ChunkParams(tool="sourmash", genomes=genomes, out_dir=tmp_path / "c1"), _LOG
        )
    )
    assert len(consumer.sub("sketch")) == 1


def test_a_missing_sketches_dir_is_refused(tmp_path, consumer) -> None:
    from repgenr.core.errors import WorkdirError
    from repgenr.stages.derep_steps import ChunkParams, dereplicate_chunk

    ctx = _workdir(tmp_path / "wd", outgroup=False)
    with pytest.raises(WorkdirError, match="--sketches-dir"):
        dereplicate_chunk(
            ChunkParams(
                tool="sourmash", genomes=_genomes(ctx), out_dir=tmp_path / "c",
                sketches_dir=tmp_path / "absent",
            ),
            _LOG,
        )  # fmt: skip


# --- glance and phylo -------------------------------------------------------------------


def test_glance_compares_the_workdir_sketches(tmp_path, fake_sourmash, consumer) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from repgenr.stages.glance import GlanceParams
    from repgenr.stages.glance import run as glance_run

    ctx = _workdir(tmp_path / "wd")
    glance_run(ctx, GlanceParams(tool="sourmash", threads=2, keep_files=True))
    assert consumer.sub("sketch") == []
    (compare,) = consumer.sub("compare")
    assert consumer.listed(compare) == [str(sketches.sketch_path(ctx, g)) for g in _genomes(ctx)]
    rows = list(
        csv.DictReader(open(ctx.workdir / "glance_wd" / "pairwise_ani.csv", encoding="utf-8"))
    )
    assert {r["genome1"] for r in rows} | {r["genome2"] for r in rows} == {
        g.name for g in _genomes(ctx)
    }


def test_phylo_sourmash_tree_reads_the_sketches_with_the_outgroup(
    tmp_path, fake_sourmash, consumer
) -> None:
    from repgenr.stages.phylo import PhyloParams
    from repgenr.stages.phylo import run as phylo_run

    ctx = _workdir(tmp_path / "wd")
    tree = phylo_run(
        ctx,
        PhyloParams(treebuilder="sourmash", all_genomes=True, allow_incomplete=True, threads=2),
    )
    assert consumer.sub("sketch") == []
    (compare,) = consumer.sub("compare")
    outgroup = next(ctx.outgroup_dir.glob("*.fasta"))
    expected = [sketches.sketch_path(ctx, g) for g in [*_genomes(ctx), outgroup]]
    assert consumer.listed(compare) == [str(p) for p in expected]
    text = tree.read_text(encoding="utf-8")
    for g in [*_genomes(ctx), outgroup]:
        assert record_name(g) in text


def test_phylo_build_reads_a_sketches_dir(tmp_path, fake_sourmash, consumer) -> None:
    from repgenr.stages.phylo import PhyloBuildParams, PhyloParams, phylo_build

    ctx = _workdir(tmp_path / "wd")
    sketch_genomes(ctx, sketch_targets(ctx), 2, _LOG)
    tree = phylo_build(
        PhyloBuildParams(
            genomes_dir=ctx.genomes_dir,
            out_dir=tmp_path / "out",
            outgroup_dir=ctx.outgroup_dir,
            outgroup_accession=ctx.workdir / "outgroup_accession.txt",
            phylo=PhyloParams(treebuilder="sourmash", threads=1),
            sketches_dir=ctx.workdir / "sketches",
        ),
        _LOG,
    )
    assert tree.is_file()
    assert consumer.sub("sketch") == []
    assert len(consumer.listed(consumer.sub("compare")[0])) == 4


# --- the classifier --------------------------------------------------------------------


def test_classifier_gathers_with_a_given_sketch(tmp_path, monkeypatch) -> None:
    import repgenr.classifiers.sourmash as cls_sm
    from repgenr.classifiers.base import ClassifyParams, registry

    chains: list[list[list[str]]] = []

    def run_chain(caps, steps, *, logger, **kwargs):  # noqa: ANN001, ANN202
        chains.append([[str(a) for a in argv] for _tool, argv in steps])
        gather = chains[-1][-1]
        Path(gather[gather.index("-o") + 1]).write_text("header\n", encoding="utf-8")

    def run_tool(caps, argv, *, logger, **kwargs):  # noqa: ANN001, ANN202
        base = argv[argv.index("--output-base") + 1]
        Path(str(base) + ".classifications.csv").write_text(
            "query_name,status,rank,fraction,lineage\n"
            "SRR1,match,species,0.9,d__Bacteria;s__X\n"
            "SRR2,match,species,0.8,d__Bacteria;s__Y\n",
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(cls_sm, "run_chain", run_chain)
    monkeypatch.setattr(cls_sm, "run_tool", run_tool)
    genomes = [tmp_path / "SRR1.fasta", tmp_path / "SRR2.fasta"]
    for g in genomes:
        g.write_text(">c\nACGT\n", encoding="utf-8")
    sig = tmp_path / "SRR1.sig.zip"
    sig.write_text("SRR1\n", encoding="utf-8")
    db, lineages = tmp_path / "gtdb.sig.zip", tmp_path / "lineages.csv"
    result = registry.create("sourmash").classify(
        genomes,
        tmp_path / "work",
        ClassifyParams(db=db, lineages=lineages, threads=1, sketches={genomes[0]: sig}),
        _LOG,
    )
    assert set(result) == {"SRR1.fasta", "SRR2.fasta"}
    by_genome = {c[-1][2]: c for c in chains}
    given = by_genome[str(sig)]
    assert len(given) == 1 and given[0][:2] == ["sourmash", "gather"]
    assert given[0][given[0].index("-k") + 1] == "31"
    other = next(c for c in chains if c[0][1] == "sketch")
    assert other[0][other[0].index("--name") + 1] == "SRR2"
    assert other[1][2] == other[0][other[0].index("-o") + 1]


def test_classifier_sketch_request_follows_its_extras() -> None:
    from repgenr.classifiers.base import registry

    adapter = registry.create("sourmash")
    assert adapter.sketch_request({}) == (31, 1000)
    assert adapter.sketch_request({"ksize": "21", "scaled": "500"}) == (21, 500)


def test_a_first_consumer_hashes_each_genome_once(tmp_path, fake_sourmash, monkeypatch) -> None:
    ctx = _workdir(tmp_path / "wd")
    genomes = [*_genomes(ctx), next(ctx.outgroup_dir.glob("*.fasta"))]
    hashed: list[str] = []
    real = sketches.file_digest

    def counting(path):  # noqa: ANN001, ANN202
        hashed.append(Path(path).name)
        return real(path)

    monkeypatch.setattr(sketches, "file_digest", counting)
    got = resolve_sketches(ctx, genomes, _LOG, 2, consumer="t")
    assert len(got) == 4 and len(fake_sourmash.calls) == 4
    assert sorted(hashed) == sorted(g.name for g in genomes)


def test_temporary_sketch_names_are_unique_and_cleared(tmp_path) -> None:
    directory = tmp_path / "sketches"
    directory.mkdir()
    out = directory / "x.sig.zip"
    first, second = sketches._partial(out), sketches._partial(out)
    assert first != second and not first.exists()
    assert first.name.startswith(".x.partial.") and first.name.endswith(".sig.zip")
    first.write_text("half", encoding="utf-8")
    (directory / ".y.partial.sig.zip").write_text("half", encoding="utf-8")  # earlier form
    out.write_text("done", encoding="utf-8")
    sketches._clear_partials(directory)
    assert sorted(p.name for p in directory.iterdir()) == ["x.sig.zip"]


def test_the_digest_cache_keeps_only_genomes_of_the_set(tmp_path, fake_sourmash) -> None:
    import json

    ctx = _workdir(tmp_path / "wd", outgroup=False)
    outside = tmp_path / "linked" / "Fam_Gen_sp_GCA_000003.1.fasta"  # a linked genome
    outside.parent.mkdir()
    outside.write_text(">x\nACGTACGT\n", encoding="utf-8")
    link = ctx.genomes_dir / outside.name
    link.symlink_to(outside)
    records = [
        GenomeRecord(
            f"GCA_{i:06d}.1", f"Fam_Gen_sp_GCA_{i:06d}.1.fasta", "local", "Fam", "Gen", "sp"
        )
        for i in range(4)
    ]
    ctx.manifest.replace_genomes(records)
    resolve_sketches(ctx, _genomes(ctx), _LOG, 2, consumer="t")
    cache = ctx.workdir / "sketches" / sketches.DIGESTS_JSON
    assert str(outside.resolve()) in json.loads(cache.read_text())
    # The linked genome leaves the set; its file still exists elsewhere.
    link.unlink()
    ctx.manifest.replace_genomes(records[:3])
    resolve_sketches(ctx, _genomes(ctx), _LOG, 2, consumer="t")
    entries = json.loads(cache.read_text())
    assert str(outside.resolve()) not in entries and len(entries) == 3
    assert not [p for p in cache.parent.iterdir() if p.name.endswith(".tmp")]
