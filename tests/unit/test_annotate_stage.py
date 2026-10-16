"""The annotate stage: cache, index, reuse, stale removal and the record."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from repgenr.annotators import bakta as bk
from repgenr.annotators.base import AnnotateParams, AnnotateResult, Annotator, registry
from repgenr.core import annotations as cache
from repgenr.core.context import WorkdirContext
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.plugins import ToolCapabilities
from repgenr.stages.annotate import AnnotateStageParams, annotate_core, run


class _FakeAnnotator(Annotator):
    capabilities = ToolCapabilities(name="fakeannot", accepted_extras=frozenset({"k"}))
    calls: list[str] = []
    version = "1.0"

    def preflight(self):
        return {"fakeannot": self.version}

    def call(self, genome, record_name, out_dir, params, logger) -> AnnotateResult:
        type(self).calls.append(record_name)
        out_dir.mkdir(parents=True, exist_ok=True)
        gff = out_dir / f"{record_name}.gff3"
        gff.write_text(f"##gff-version 3\n##FASTA\n>{record_name}\nACGT\n", encoding="utf-8")
        return AnnotateResult(gff=gff)


@pytest.fixture
def fake(register_tool):
    _FakeAnnotator.calls = []
    _FakeAnnotator.version = "1.0"
    register_tool(registry, "fakeannot", _FakeAnnotator)
    return _FakeAnnotator


def _ctx(workdir: Path) -> WorkdirContext:
    return WorkdirContext(workdir, create=True, logger=logging.getLogger("test-annotate"))


def test_index_round_trip_and_currency(tmp_path: Path) -> None:
    rec = cache.AnnotationRecord(
        record_name="a",
        set="genomes",
        genome_path="genomes/a.fasta",
        genome_digest="d1",
        annotator="bakta",
        annotator_version="1.12",
        db_version="6.0-light",
        extra="k=1",
        gff_digest="g1",
        completed="now",
    )
    cache.write_index(tmp_path, {"a": rec})
    assert cache.read_index(tmp_path) == {"a": rec}
    gff = cache.gff_path(tmp_path, "a")
    gff.parent.mkdir()
    gff.write_text("x", encoding="utf-8")
    base = dict(
        genome_digest="d1",
        annotator="bakta",
        annotator_version="1.12",
        db_version="6.0-light",
        extra="k=1",
        gff=gff,
    )
    assert cache.is_current(rec, **base)
    for key, value in [
        ("genome_digest", "d2"),
        ("annotator", "prokka"),
        ("annotator_version", "1.13"),
        ("db_version", "7.0-light"),
        ("extra", "k=2"),
    ]:
        assert not cache.is_current(rec, **{**base, key: value})
    gff.unlink()
    assert not cache.is_current(rec, **base)
    assert cache.read_index(tmp_path / "none") == {}


def test_stage_annotates_reuses_and_removes(workdir: Path, genome_files, fake, caplog) -> None:
    ctx = _ctx(workdir)
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 1}))
    assert (summary.annotated, summary.reused, summary.removed) == (3, 0, 0)
    index = cache.read_index(ctx.annotations_dir)
    assert set(index) == {g.stem for g in genome_files}
    assert all(
        r.set == "genomes" and r.extra == "k=1" and r.annotator == "fakeannot"
        for r in index.values()
    )
    assert all(cache.gff_path(ctx.annotations_dir, n).is_file() for n in index)
    rec = ctx.config.stages["annotate"]
    assert (
        rec.tool == "fakeannot"
        and rec.params["annotated"] == 3
        and rec.params["db_version"] is None
    )
    # Second run: everything current.
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 1}))
    assert (summary.annotated, summary.reused) == (0, 3) and len(fake.calls) == 3
    # A changed genome and a changed extra are redone; the rest reused.
    genome_files[0].write_text(">seq0\nTTTT\n", encoding="utf-8")
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 1}))
    assert (summary.annotated, summary.reused) == (1, 2)
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 2}))
    assert summary.annotated == 3
    # A genome leaving the set drops its record and directory.
    genome_files[2].unlink()
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 2}))
    assert summary.removed == 1 and len(cache.read_index(ctx.annotations_dir)) == 2
    assert not cache.record_dir(ctx.annotations_dir, genome_files[2].stem).exists()
    # A new tool version re-annotates; --force re-annotates all.
    fake.version = "2.0"
    assert (
        run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 2})).annotated
        == 2
    )
    ctx.force = True
    assert (
        run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, extra={"k": 2})).annotated
        == 2
    )


def test_stage_includes_the_outgroup(workdir: Path, genome_files, fake) -> None:
    og = workdir / "outgroup"
    og.mkdir()
    (og / "Out_gen_sp_GCA_000009.fasta").write_text(">o\nACGT\n", encoding="utf-8")
    (workdir / "outgroup_accession.txt").write_text("GCA_000009\n", encoding="utf-8")
    ctx = _ctx(workdir)
    run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))
    index = cache.read_index(ctx.annotations_dir)
    assert index["Out_gen_sp_GCA_000009"].set == "outgroup"
    assert sum(r.set == "genomes" for r in index.values()) == 3


def test_partial_directories_are_cleared_and_never_published(
    workdir: Path, genome_files, fake, monkeypatch
) -> None:
    ctx = _ctx(workdir)
    ctx.annotations_dir.mkdir(parents=True)
    leftover = cache.partial_dir(ctx.annotations_dir, "stale")
    leftover.mkdir()
    run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))
    assert not leftover.exists()

    def failing(self, genome, record_name, out_dir, params, logger):
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "half.gff3").write_text("x", encoding="utf-8")
        raise RuntimeError("tool crashed")

    monkeypatch.setattr(_FakeAnnotator, "call", failing)
    genome_files[1].write_text(">seq1\nGGGG\n", encoding="utf-8")
    with pytest.raises(RuntimeError):
        run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))
    name = genome_files[1].stem
    assert not cache.partial_dir(ctx.annotations_dir, name).exists()
    # The previous annotation of that genome survives the failed attempt.
    assert cache.gff_path(ctx.annotations_dir, name).is_file()


def test_bakta_needs_its_database_before_any_work(workdir: Path, genome_files, monkeypatch) -> None:
    monkeypatch.delenv(bk.BAKTA_DB_ENV, raising=False)
    ctx = _ctx(workdir)
    with pytest.raises(UserInputError, match="bakta_db download"):
        run(ctx, AnnotateStageParams(tool="bakta", all_genomes=True))
    with pytest.raises(UserInputError, match="Unknown annotator"):
        run(ctx, AnnotateStageParams(tool="nope", all_genomes=True))


def test_annotate_core_stages_gz_and_publishes_atomically(tmp_path: Path, fake) -> None:
    import gzip

    g = tmp_path / "Fam_Gen_sp_GCA_1.fasta.gz"
    with gzip.open(g, "wt") as fo:
        fo.write(">c\nACGT\n")
    result, versions = annotate_core(
        g,
        "Fam_Gen_sp_GCA_1",
        tmp_path / "out" / "Fam_Gen_sp_GCA_1",
        tmp_path / "scratch",
        "fakeannot",
        AnnotateParams(),
        logging.getLogger("t"),
    )
    assert result.gff == tmp_path / "out" / "Fam_Gen_sp_GCA_1" / "Fam_Gen_sp_GCA_1.gff3"
    assert result.gff.is_file() and versions == {"fakeannot": "1.0"}
    assert not (tmp_path / "scratch" / "inputs").exists()
    assert not any(p.name.startswith(".") for p in (tmp_path / "out").iterdir())


def test_index_digest_follows_the_gff_digests(tmp_path: Path) -> None:
    rec = cache.AnnotationRecord("a", "genomes", "p", "d", "bakta", "1", "6", "", "g1", "t")
    other = cache.AnnotationRecord("a", "genomes", "p", "d", "bakta", "1", "6", "", "g2", "t")
    assert cache.index_digest({"a": rec}, ["a"]) != cache.index_digest({"a": other}, ["a"])
    assert cache.index_digest({"a": rec}, ["a", "b"]) != cache.index_digest({"a": rec}, ["a"])
    assert json.dumps(cache.COLUMNS)  # the columns are plain strings


def test_annotate_build_step_writes_versions(tmp_path: Path, fake) -> None:
    from repgenr.stages.annotate import AnnotateBuildParams, annotate_build

    g = tmp_path / "Fam_Gen_sp_GCA_7.fasta"
    g.write_text(">c\nACGT\n", encoding="utf-8")
    result = annotate_build(
        AnnotateBuildParams(
            genome=g, out_dir=tmp_path / "out", tool="fakeannot", versions_out=tmp_path / "v.yml"
        ),
        logging.getLogger("t"),
    )
    assert result.gff == tmp_path / "out" / "Fam_Gen_sp_GCA_7" / "Fam_Gen_sp_GCA_7.gff3"
    assert "fakeannot: 1.0" in (tmp_path / "v.yml").read_text(encoding="utf-8")
    assert not (tmp_path / "out" / "scratch").exists()
    with pytest.raises(UserInputError, match="Genome not found"):
        annotate_build(
            AnnotateBuildParams(
                genome=tmp_path / "no.fa", out_dir=tmp_path / "o2", tool="fakeannot"
            ),
            logging.getLogger("t"),
        )


def test_a_genome_outside_the_set_keeps_its_annotation(workdir: Path, genome_files, fake) -> None:
    """Annotation costs minutes: leaving the set must not throw it away."""
    ctx = _ctx(workdir)
    run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for g in genome_files[:2]:
        (reps / g.name).write_text(g.read_text(encoding="utf-8"), encoding="utf-8")
    (workdir / "derep" / "clusters.tsv").write_text(
        "representative\tmember\n" + "".join(f"{g.name}\t{g.name}\n" for g in genome_files[:2]),
        encoding="utf-8",
    )
    summary = run(ctx, AnnotateStageParams(tool="fakeannot"))
    assert summary.removed == 0 and summary.reused == 2
    index = cache.read_index(ctx.annotations_dir)
    assert len(index) == 3  # the third genome's annotation is still there
    assert cache.gff_path(ctx.annotations_dir, genome_files[2].stem).is_file()


def test_the_set_column_follows_a_genome_that_becomes_the_outgroup(
    workdir: Path, genome_files, fake
) -> None:
    ctx = _ctx(workdir)
    run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))
    name = genome_files[2].stem
    og = workdir / "outgroup"
    og.mkdir()
    (og / genome_files[2].name).write_text(
        genome_files[2].read_text(encoding="utf-8"), encoding="utf-8"
    )
    (workdir / "outgroup_accession.txt").write_text("GCA_000003\n", encoding="utf-8")
    genome_files[2].unlink()
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))
    assert summary.annotated == 0 and summary.reused == 3
    assert cache.read_index(ctx.annotations_dir)[name].set == "outgroup"


def test_two_genomes_with_one_record_name_are_named_in_the_error(
    workdir: Path, genome_files, fake
) -> None:
    og = workdir / "outgroup"
    og.mkdir()
    (og / genome_files[0].name).write_text(">x\nACGT\n", encoding="utf-8")
    (workdir / "outgroup_accession.txt").write_text("GCA_000001\n", encoding="utf-8")
    ctx = _ctx(workdir)
    with pytest.raises(UserInputError, match=r"share a record name.*GCA_000001"):
        run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True))


def test_a_changed_database_version_re_annotates(workdir: Path, genome_files, fake) -> None:
    class _DbAnnotator(_FakeAnnotator):
        db_version = "6.0-light"

        def database_version(self, db) -> str:  # noqa: ANN001
            return type(self).db_version

    _DbAnnotator.calls = []
    _DbAnnotator.requires_database = True
    ctx = _ctx(workdir)
    from repgenr.annotators.base import registry as annot_registry

    annot_registry.register("fakedb", _DbAnnotator, replace=True)
    try:
        db = workdir / "db"
        db.mkdir(parents=True)
        (db / "version.json").write_text('{"major": 6, "minor": 0}', encoding="utf-8")
        params = AnnotateStageParams(tool="fakedb", all_genomes=True, bakta_db=str(db))
        assert run(ctx, params).annotated == 3
        assert run(ctx, params).reused == 3
        _DbAnnotator.db_version = "7.0-light"
        summary = run(ctx, params)
        assert summary.annotated == 3 and summary.reused == 0
        rec = ctx.config.stages["annotate"]
        assert rec.params["db_version"] == "7.0-light"
    finally:
        annot_registry.unregister("fakedb")


def test_the_database_is_part_of_the_resume_inputs(workdir: Path, monkeypatch) -> None:
    from repgenr.cli.base import STAGE_INPUTS

    ctx = _ctx(workdir)
    db = workdir / "db"
    db.mkdir(parents=True)
    (db / "version.json").write_text("{}", encoding="utf-8")
    flagged = STAGE_INPUTS["annotate"](
        ctx, AnnotateStageParams(tool="bakta", bakta_db=str(db), all_genomes=True)
    )
    assert db / "version.json" in flagged
    monkeypatch.setenv("BAKTA_DB", str(db))
    from_env = STAGE_INPUTS["annotate"](ctx, AnnotateStageParams(tool="bakta", all_genomes=True))
    assert db / "version.json" in from_env
    # prokka needs none, so nothing is added.
    assert STAGE_INPUTS["annotate"](ctx, AnnotateStageParams(tool="prokka")) == [
        ctx.representatives_dir,
        ctx.outgroup_dir,
    ]


def test_boolean_tuning_is_normalised_in_the_index(workdir: Path, genome_files, fake) -> None:
    ctx = _ctx(workdir)

    def once(**extra) -> AnnotateStageParams:
        return AnnotateStageParams(tool="fakeannot", all_genomes=True, extra=extra)

    assert run(ctx, once(k=True)).annotated == 3
    # "true" describes the same setting as True and must not re-annotate.
    assert run(ctx, once(k="true")).reused == 3
    assert run(ctx, once(k="no")).annotated == 3
    # A number stays a number, so 1 is not read as "on".
    assert run(ctx, once(k=1)).annotated == 3
    assert run(ctx, once(k="1")).reused == 3


def test_a_genome_that_stops_being_a_representative_keeps_its_annotation(
    workdir: Path, genome_files, fake
) -> None:
    """Re-dereplication moves a genome out of representatives/ but not out of
    the workdir; its annotation is expensive and must survive."""
    ctx = _ctx(workdir)
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    clusters = workdir / "derep" / "clusters.tsv"
    for g in genome_files:
        (reps / g.name).write_text(g.read_text(encoding="utf-8"), encoding="utf-8")
    clusters.write_text(
        "representative\tmember\n" + "".join(f"{g.name}\t{g.name}\n" for g in genome_files),
        encoding="utf-8",
    )
    assert run(ctx, AnnotateStageParams(tool="fakeannot")).annotated == 3
    # A new dereplication keeps two representatives; the third genome stays in genomes/.
    (reps / genome_files[2].name).unlink()
    clusters.write_text(
        "representative\tmember\n"
        + "".join(f"{g.name}\t{g.name}\n" for g in genome_files[:2])
        + f"{genome_files[0].name}\t{genome_files[2].name}\n",
        encoding="utf-8",
    )
    summary = run(ctx, AnnotateStageParams(tool="fakeannot"))
    assert summary.removed == 0 and summary.reused == 2
    assert cache.gff_path(ctx.annotations_dir, genome_files[2].stem).is_file()


def test_annotate_index_step(tmp_path: Path, monkeypatch) -> None:
    """Gathered annotate-build records, as Nextflow stages them (links to the
    record directories), are indexed without probing the annotator."""
    from repgenr.stages.annotate import AnnotateIndexParams, annotate_index
    from repgenr.stages.coregene import annotated_set

    probed: list[str] = []
    monkeypatch.setattr(
        registry.get("bakta"), "preflight", lambda self: probed.append("x") or {}, raising=False
    )
    elsewhere = tmp_path / "work"
    gathered = tmp_path / "annotations"
    gathered.mkdir()
    for name in ("Fam_Gen_sp_GCA_1.1", "Fam_Gen_sp_GCA_2.1"):
        record = elsewhere / name
        record.mkdir(parents=True)
        (record / f"{name}.gff3").write_text(f"##gff-version 3\n##FASTA\n>{name}\nACGT\n", "utf-8")
        (gathered / name).symlink_to(record, target_is_directory=True)
    index = annotate_index(
        AnnotateIndexParams(annotations_dir=gathered, tool="bakta"), logging.getLogger("t")
    )
    assert index == gathered / cache.ANNOTATIONS_TSV and probed == []
    rows = cache.read_index(gathered)
    assert sorted(rows) == ["Fam_Gen_sp_GCA_1.1", "Fam_Gen_sp_GCA_2.1"]
    assert {r.annotator for r in rows.values()} == {"bakta"}
    # The core-gene pass reads it as it reads a workdir's.
    got = annotated_set(gathered, ["Fam_Gen_sp_GCA_2.1", "Fam_Gen_sp_GCA_1.1"])
    assert [g.name for g in got.gffs] == ["Fam_Gen_sp_GCA_2.1.gff3", "Fam_Gen_sp_GCA_1.1.gff3"]

    with pytest.raises(UserInputError, match="not a directory"):
        annotate_index(
            AnnotateIndexParams(annotations_dir=tmp_path / "nope"), logging.getLogger("t")
        )
    with pytest.raises(UserInputError, match="Unknown annotator"):
        annotate_index(
            AnnotateIndexParams(annotations_dir=gathered, tool="nope"), logging.getLogger("t")
        )
    unfinished = tmp_path / "unfinished"
    (unfinished / "Fam_Gen_sp_GCA_3.1").mkdir(parents=True)
    with pytest.raises(WorkdirError, match="did not finish"):
        annotate_index(AnnotateIndexParams(annotations_dir=unfinished), logging.getLogger("t"))
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(UserInputError, match="No annotated records"):
        annotate_index(AnnotateIndexParams(annotations_dir=empty), logging.getLogger("t"))


def test_each_gathered_record_is_indexed_with_the_tool_that_made_it(
    tmp_path: Path, fake, caplog
) -> None:
    """annotate-build writes annotation.json; annotate-index reads it per
    record, and uses --tool only for a record without one."""
    from repgenr.stages.annotate import (
        AnnotateBuildParams,
        AnnotateIndexParams,
        annotate_build,
        annotate_index,
    )
    from repgenr.stages.coregene import annotated_set

    gathered = tmp_path / "annotations"
    for name in ("Fam_Gen_sp_GCA_1.1", "Fam_Gen_sp_GCA_2.1"):
        genome = tmp_path / f"{name}.fasta"
        genome.write_text(">c\nACGT\n", encoding="utf-8")
        annotate_build(
            AnnotateBuildParams(genome=genome, out_dir=gathered, tool="fakeannot", extra={"k": 1}),
            logging.getLogger("t"),
        )
    provenance = cache.read_provenance(gathered, "Fam_Gen_sp_GCA_1.1")
    assert provenance == {
        "annotator": "fakeannot",
        "annotator_version": "1.0",
        "db_version": "",
        "extra": "k=1",
    }
    # Told bakta, the index still names the tool that made the records.
    with caplog.at_level(logging.INFO):
        annotate_index(
            AnnotateIndexParams(annotations_dir=gathered, tool="bakta"), logging.getLogger("t")
        )
    rows = cache.read_index(gathered)
    assert {(r.annotator, r.annotator_version, r.extra) for r in rows.values()} == {
        ("fakeannot", "1.0", "k=1")
    }
    assert "annotated with fakeannot, not bakta" in caplog.text

    # A record from an earlier version (no annotation.json) takes --tool, so
    # a set mixing it with the others is refused by the core-gene pass.
    old = gathered / "Fam_Gen_sp_GCA_3.1"
    old.mkdir()
    (old / "Fam_Gen_sp_GCA_3.1.gff3").write_text("##gff-version 3\n##FASTA\n>x\nAC\n", "utf-8")
    annotate_index(
        AnnotateIndexParams(annotations_dir=gathered, tool="bakta"), logging.getLogger("t")
    )
    rows = cache.read_index(gathered)
    assert rows["Fam_Gen_sp_GCA_3.1"].annotator == "bakta"
    assert rows["Fam_Gen_sp_GCA_3.1"].annotator_version == ""
    with pytest.raises(UserInputError, match="mix annotators"):
        annotated_set(gathered, sorted(rows))
    # An unreadable file counts as absent.
    (old / cache.PROVENANCE).write_text("{not json", encoding="utf-8")
    assert cache.read_provenance(gathered, "Fam_Gen_sp_GCA_3.1") is None


def test_annotate_jobs_splits_the_threads(monkeypatch) -> None:
    from repgenr.stages import annotate as stage

    monkeypatch.setattr(stage, "usable_cpus", lambda: 64)
    assert stage.annotate_jobs(16, 0, 100) == (4, 4)  # auto: one genome per four threads
    assert stage.annotate_jobs(16, 0, 2) == (2, 8)  # never more than the genomes
    assert stage.annotate_jobs(3, 0, 100) == (1, 3)
    assert stage.annotate_jobs(8, 3, 100) == (3, 2)
    assert stage.annotate_jobs(2, 5, 100) == (2, 1)  # at least one thread each
    assert stage.annotate_jobs(16, 0, 0) == (1, 16)
    monkeypatch.setattr(stage, "usable_cpus", lambda: 2)
    assert stage.annotate_jobs(16, 0, 100) == (2, 8)  # capped by the usable CPUs


def test_stage_annotates_genomes_in_parallel(
    workdir: Path, genome_files, fake, monkeypatch
) -> None:
    """Several genomes run at once with the threads split between them, each
    in its own scratch, and every one lands in the index."""
    import threading
    import time

    from repgenr.stages import annotate as stage

    monkeypatch.setattr(stage, "usable_cpus", lambda: 64)
    lock = threading.Lock()
    state = {"running": 0, "peak": 0}
    threads_seen: list[int] = []
    inputs_seen: list[Path] = []
    real_call = _FakeAnnotator.call

    def call(self, genome, record_name, out_dir, params, logger):
        with lock:
            state["running"] += 1
            state["peak"] = max(state["peak"], state["running"])
            threads_seen.append(params.threads)
            inputs_seen.append(genome.parent)
        time.sleep(0.1)
        try:
            return real_call(self, genome, record_name, out_dir, params, logger)
        finally:
            with lock:
                state["running"] -= 1

    monkeypatch.setattr(_FakeAnnotator, "call", call)
    ctx = _ctx(workdir)
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, threads=8))
    assert summary.annotated == 3
    assert state["peak"] == 2
    assert threads_seen == [4, 4, 4]
    assert len(cache.read_index(ctx.annotations_dir)) == 3
    assert not (ctx.scratch_dir / "annotate").exists()


def test_a_failed_genome_keeps_the_others_in_the_index(
    workdir: Path, genome_files, fake, monkeypatch
) -> None:
    """One genome fails while another runs: the running one finishes and is
    indexed, the error is raised, and a rerun annotates only what is missing."""
    import threading

    gate = threading.Event()
    real_call = _FakeAnnotator.call

    def call(self, genome, record_name, out_dir, params, logger):
        if record_name.endswith("000001"):
            gate.wait(5)
            return real_call(self, genome, record_name, out_dir, params, logger)
        if record_name.endswith("000002"):
            gate.set()
            raise RuntimeError("annotator failed")
        return real_call(self, genome, record_name, out_dir, params, logger)

    monkeypatch.setattr(_FakeAnnotator, "call", call)
    ctx = _ctx(workdir)
    with pytest.raises(RuntimeError, match="annotator failed"):
        run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, num_processes=2))
    assert list(cache.read_index(ctx.annotations_dir)) == [
        "Francisellaceae_francisella_tularensis_GCA_000001"
    ]

    monkeypatch.setattr(_FakeAnnotator, "call", real_call)
    _FakeAnnotator.calls = []
    summary = run(ctx, AnnotateStageParams(tool="fakeannot", all_genomes=True, num_processes=2))
    assert summary.reused == 1 and summary.annotated == 2
    assert sorted(_FakeAnnotator.calls) == [
        "Francisellaceae_francisella_tularensis_GCA_000002",
        "Francisellaceae_francisella_tularensis_GCA_000003",
    ]


def test_parallel_builds_share_an_output_directory(tmp_path: Path, fake) -> None:
    """annotate-build runs side by side in one directory (phylo-build) without
    removing each other's scratch."""
    from repgenr.core.executors import parallel_each
    from repgenr.stages.annotate import AnnotateBuildParams, annotate_build

    genomes = []
    for i in range(4):
        g = tmp_path / f"g{i}.fasta.gz"
        import gzip

        g.write_bytes(gzip.compress(f">g{i}\nACGT\n".encode()))
        genomes.append(g)
    out = tmp_path / "annotations"
    parallel_each(
        lambda g: annotate_build(
            AnnotateBuildParams(genome=g, out_dir=out, tool="fakeannot", threads=1),
            logging.getLogger("t"),
        ),
        genomes,
        4,
        lambda _g, _r: None,
    )
    assert sorted(p.name for p in out.iterdir()) == ["g0", "g1", "g2", "g3"]
