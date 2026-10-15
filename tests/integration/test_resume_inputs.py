"""The per-stage input specs feed the v2 fingerprint (STAGE_INPUTS table)."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

from repgenr.cli import base as cli
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import CLUSTERS_TSV, CORE_SNP_FASTA, TREE_NWK
from repgenr.core.inputs import ABSENT
from repgenr.core.manifest import GenomeRecord
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import (
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepResult,
    registry,
)
from repgenr.stages.dereplicate import DereplicateParams


def _digests(ctx: WorkdirContext, stage: str, params) -> dict[str, str]:
    return cli._stage_input_digests(ctx, stage, params)


def test_dereplicate_inputs_track_genomes_dir(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.genomes_dir.mkdir()
    (ctx.genomes_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    params = SimpleNamespace()
    before = _digests(ctx, "dereplicate", params)
    assert set(before) == {"genomes", "manifest"}

    (ctx.genomes_dir / "g2.fasta").write_text(">g2\nACGT\n", encoding="utf-8")
    after = _digests(ctx, "dereplicate", params)
    assert after["genomes"] != before["genomes"]


def test_dereplicate_inputs_track_manifest_quality_and_taxonomy(tmp_path: Path) -> None:
    """A manifest-only change (CheckM quality or taxonomy, no genome files
    touched) must still invalidate a prior dereplicate resume, since keeper
    and --reduce both read the manifest."""
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.genomes_dir.mkdir()
    (ctx.genomes_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    ctx.manifest.upsert_many([GenomeRecord(accession="GCF_1", filename="g1.fasta")])
    params = SimpleNamespace()
    before = _digests(ctx, "dereplicate", params)
    assert "manifest" in before

    ctx.manifest.upsert_many(
        [GenomeRecord(accession="GCF_1", filename="g1.fasta", completeness=98.0, contamination=0.5)]
    )
    after_quality = _digests(ctx, "dereplicate", params)
    assert after_quality["manifest"] != before["manifest"]

    ctx.manifest.upsert_many(
        [
            GenomeRecord(
                accession="GCF_1",
                filename="g1.fasta",
                completeness=98.0,
                contamination=0.5,
                genus="g",
                species="s",
            )
        ]
    )
    after_taxonomy = _digests(ctx, "dereplicate", params)
    assert after_taxonomy["manifest"] != after_quality["manifest"]


def test_dereplicate_manifest_digest_excludes_derep_status_and_representative(
    tmp_path: Path,
) -> None:
    """dereplicate itself writes derep_status/representative back into the
    manifest (_update_manifest); if its own manifest digest hashed those
    fields, every workdir would pay one spurious re-dereplication on its first
    resume (the digest computed before the run would never match the one
    computed after, since the run just changed them). dereplicate must only be
    invalidated by what it reads (taxonomy, quality), not what it writes."""
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.genomes_dir.mkdir()
    (ctx.genomes_dir / "g1.fasta").write_text(">g1\nACGT\n", encoding="utf-8")
    ctx.manifest.upsert_many([GenomeRecord(accession="GCF_1", filename="g1.fasta")])
    params = SimpleNamespace()
    before = _digests(ctx, "dereplicate", params)

    # Simulate what dereplicate.run()'s _update_manifest does after a run.
    ctx.manifest.set_derep_status_many([("GCF_1", "representative", None)])
    after = _digests(ctx, "dereplicate", params)

    assert after["manifest"] == before["manifest"]

    # tree2tax, in contrast, reads derep status, so it must still be tracked there.
    ctx.tree_dir.mkdir()
    (ctx.tree_dir / TREE_NWK).write_text("(a,b);\n", encoding="utf-8")
    tt_params = SimpleNamespace(include_dereplicated=False)
    tt_before = _digests(ctx, "tree2tax", tt_params)
    ctx.manifest.set_derep_status_many([("GCF_1", "contained", "other.fasta")])
    tt_after = _digests(ctx, "tree2tax", tt_params)
    assert tt_after["manifest"] != tt_before["manifest"]


def test_dereplicate_resume_skips_then_reruns_on_quality_edit(tmp_path: Path, monkeypatch) -> None:
    """End-to-end resume parity for the real dereplicate stage: an identical
    second invocation is skipped, and a quality-only manifest edit between
    invocations forces a rerun (catches the self-invalidation regression: if
    derep_status leaked into the digest, run 2 would rerun instead of skip)."""

    class _ResumeNoRep(Dereplicator):
        """Keeps every input genome as its own representative (no clustering)."""

        capabilities = ToolCapabilities(name="resumenorep", supports_native_scaling=True)

        def preflight(self) -> dict[str, str]:
            return {"resumenorep": "1.0"}

        def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
            genomes = list(genomes)
            return DerepResult(
                representatives=list(genomes),
                clusters={g.name: [] for g in genomes},
                genome_status={g.name: STATUS_REPRESENTATIVE for g in genomes},
            )

    registry._load()
    registry.register("resumenorep", _ResumeNoRep, replace=True)
    try:
        gdir = tmp_path / "genomes"
        gdir.mkdir(parents=True)
        for i in range(2):
            (gdir / f"Fam_g_s_GCF_10000{i}.1.fasta").write_text(">x\nACGT\n", encoding="utf-8")

        seed = WorkdirContext(tmp_path, create=True)
        seed.manifest.upsert_many(
            [
                GenomeRecord(accession=f"GCF_10000{i}.1", filename=f"Fam_g_s_GCF_10000{i}.1.fasta")
                for i in range(2)
            ]
        )
        seed.close()

        def _flush_and_read_log() -> str:
            for h in logging.getLogger("repgenr").handlers:
                h.flush()
            return (tmp_path / "repgenr.log").read_text(encoding="utf-8")

        build_params = lambda: DereplicateParams(tool="resumenorep")  # noqa: E731
        monkeypatch.setitem(cli._RUN_STATE, "force", False)
        monkeypatch.setitem(cli._RUN_STATE, "log_level", logging.INFO)

        cli._run("dereplicate", tmp_path, build_params, create=True)  # run 1: cold
        log1 = _flush_and_read_log()
        assert log1.count("Dereplicating 2 genomes") == 1
        assert "already completed" not in log1

        cli._run("dereplicate", tmp_path, build_params, create=True)  # run 2: identical -> skip
        log2 = _flush_and_read_log()
        assert log2.count("already completed") == 1
        assert log2.count("Dereplicating 2 genomes") == 1  # did not rerun

        edit = WorkdirContext(tmp_path)
        edit.manifest.upsert_many(
            [
                GenomeRecord(
                    accession="GCF_100000.1",
                    filename="Fam_g_s_GCF_100000.1.fasta",
                    completeness=98.0,
                    contamination=0.5,
                )
            ]
        )
        edit.close()

        # run 3: quality changed -> rerun
        cli._run("dereplicate", tmp_path, build_params, create=True)
        log3 = _flush_and_read_log()
        assert log3.count("already completed") == 1  # unchanged: run 3 did not skip
        assert log3.count("Dereplicating 2 genomes") == 2  # reran
    finally:
        registry._classes.pop("resumenorep", None)


def test_phylo_inputs_never_include_own_outputs(tmp_path: Path) -> None:
    """core_snp.fasta is REGENERATED by phylo under msa_source=snptype, so it
    must not be declared as an input (a self-reference would force a spurious
    rerun on the second identical invocation)."""
    ctx = WorkdirContext(tmp_path, create=True)
    aligner = _digests(ctx, "phylo", SimpleNamespace(all_genomes=False, msa_source="aligner"))
    snptype = _digests(ctx, "phylo", SimpleNamespace(all_genomes=False, msa_source="snptype"))
    core_snp_key = f"snp/{CORE_SNP_FASTA}"
    assert core_snp_key not in aligner
    assert core_snp_key not in snptype
    assert aligner == snptype
    assert "derep/representatives" in aligner and "outgroup" in aligner


def test_vgenome_inputs_are_vmetadata_outputs_not_its_own(tmp_path: Path) -> None:
    """vgenome WRITES selection.tsv, so its declared inputs must be the
    vmetadata download artifacts instead."""
    ctx = WorkdirContext(tmp_path, create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(">a\nACGT\n", encoding="utf-8")
    (wd / "virus_records.json").write_text("[]", encoding="utf-8")

    digests = _digests(ctx, "vgenome", SimpleNamespace())
    assert "selection.tsv" not in digests
    assert "virus_download_wd/download.fa" in digests
    assert "virus_download_wd/virus_records.json" in digests

    # writing selection.tsv (what vgenome does) must not change its own digests
    (ctx.workdir / "selection.tsv").write_text("accession\n", encoding="utf-8")
    assert _digests(ctx, "vgenome", SimpleNamespace()) == digests

    # but a changed vmetadata download does
    (wd / "download.fa").write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    assert _digests(ctx, "vgenome", SimpleNamespace()) != digests


def test_phylo_all_genomes_switches_input_dir(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    reps = _digests(ctx, "phylo", SimpleNamespace(all_genomes=False, msa_source="aligner"))
    allg = _digests(ctx, "phylo", SimpleNamespace(all_genomes=True, msa_source="aligner"))
    assert "derep/representatives" in reps and "genomes" not in reps
    assert "genomes" in allg and "derep/representatives" not in allg


def test_tree2tax_inputs_track_tree_manifest_and_clusters(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.tree_dir.mkdir()
    tree = ctx.tree_dir / TREE_NWK
    tree.write_text("(a,b);\n", encoding="utf-8")

    plain = _digests(ctx, "tree2tax", SimpleNamespace(include_dereplicated=False))
    assert set(plain) == {
        f"tree/{TREE_NWK}",
        "manifest",
        "segments.tsv",
        "outgroup_accession.txt",
        "outgroup",
    }

    withderep = _digests(ctx, "tree2tax", SimpleNamespace(include_dereplicated=True))
    clusters_key = f"derep/{CLUSTERS_TSV}"
    assert clusters_key in withderep
    assert withderep[clusters_key] == ABSENT  # not written yet -> stable sentinel

    tree.write_text("(a,b,c);\n", encoding="utf-8")
    changed = _digests(ctx, "tree2tax", SimpleNamespace(include_dereplicated=False))
    assert changed[f"tree/{TREE_NWK}"] != plain[f"tree/{TREE_NWK}"]


def test_genome_inputs_track_selection_tsv(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    sel = ctx.workdir / "selection.tsv"
    sel.write_text("accession\nGCF_1\n", encoding="utf-8")
    before = _digests(ctx, "genome", SimpleNamespace())
    sel.write_text("accession\nGCF_1\nGCF_2\n", encoding="utf-8")
    after = _digests(ctx, "genome", SimpleNamespace())
    assert after["selection.tsv"] != before["selection.tsv"]


def test_metadata_has_no_local_inputs(tmp_path: Path) -> None:
    ctx = WorkdirContext(tmp_path, create=True)
    assert _digests(ctx, "metadata", SimpleNamespace()) == {}


def test_metadata_inputs_track_a_local_metadata_path(tmp_path: Path) -> None:
    """--metadata-path is a local file the selection is computed from; replacing
    the table must invalidate a prior metadata resume, as an accession list
    does for reads."""
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    table = tmp_path / "bac120_metadata_r232.tsv"
    table.write_text("accession\nRS_GCF_1\n", encoding="utf-8")
    params = SimpleNamespace(metadata_path=str(table))
    before = _digests(ctx, "metadata", params)
    assert before
    table.write_text("accession\nRS_GCF_1\nRS_GCF_2\n", encoding="utf-8")
    after = _digests(ctx, "metadata", params)
    assert after != before


def _write_gtdb_table(path: Path, species: list[str]) -> None:
    import gzip

    prefix = "d__Bacteria;p__P;c__C;o__O;f__Francisellaceae;g__Francisella;s__Francisella"
    lines = [
        "accession\tgtdb_genome_representative\tgtdb_taxonomy\tncbi_genbank_assembly_accession"
    ]
    for i, sp in enumerate(species, start=1):
        acc = f"GCF_00000{i}.1"
        lines.append(f"RS_{acc}\tRS_{acc}\t{prefix} {sp}\t{acc}")
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")


def test_metadata_nodownload_inputs_track_the_reused_table(tmp_path: Path) -> None:
    """--nodownload reuses the table in the workdir, so it is a declared input."""
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = SimpleNamespace(
        metadata_path=None, nodownload=True, release="232.0", version="bac120", source="tsv"
    )
    table = ctx.workdir / "bac120_metadata_r232.tsv.gz"
    _write_gtdb_table(table, ["tularensis", "philomiragia"])
    before = _digests(ctx, "metadata", params)
    assert before
    _write_gtdb_table(table, ["tularensis", "tularensis", "philomiragia"])
    assert _digests(ctx, "metadata", params) != before
    # without --nodownload the table is downloaded afresh: no file input
    assert _digests(ctx, "metadata", SimpleNamespace(**{**vars(params), "nodownload": False})) == {}


def test_metadata_nodownload_reruns_when_the_reused_table_is_replaced(
    tmp_path: Path, monkeypatch
) -> None:
    from repgenr.stages.metadata import MetadataParams

    wd = tmp_path / "wd"
    wd.mkdir()
    table = wd / "bac120_metadata_r232.tsv.gz"
    _write_gtdb_table(table, ["tularensis", "philomiragia"])

    def build_params() -> MetadataParams:
        return MetadataParams(
            dataset="all",
            level="species",
            release="232.0",
            version="bac120",
            target_genus="Francisella",
            target_species="tularensis",
            nodownload=True,
        )

    def read_log() -> str:
        for h in logging.getLogger("repgenr").handlers:
            h.flush()
        return (wd / "repgenr.log").read_text(encoding="utf-8")

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(cli._RUN_STATE, "log_level", logging.INFO)
    cli._run("metadata", wd, build_params, create=True)
    cli._run("metadata", wd, build_params, create=True)
    assert read_log().count("already completed") == 1

    _write_gtdb_table(table, ["tularensis", "tularensis", "philomiragia"])
    cli._run("metadata", wd, build_params, create=True)
    log = read_log()
    assert "input 'bac120_metadata_r232.tsv.gz' changed since last completion" in log
    assert log.count("already completed") == 1
    ctx = WorkdirContext(wd)
    assert ctx.config.stages["metadata"].params["selected_count"] == 2


def test_metadata_nodownload_first_run_that_downloads_skips_next_time(
    tmp_path: Path, monkeypatch
) -> None:
    """With no table yet, --nodownload downloads it; the inputs are digested after
    the run, so the next identical run skips instead of seeing a new input."""
    from repgenr.core import http
    from repgenr.stages.metadata import MetadataParams

    wd = tmp_path / "wd"
    wd.mkdir()
    downloads: list[str] = []

    def fake_download(url, dest, logger=None, **kwargs):
        downloads.append(url)
        _write_gtdb_table(dest, ["tularensis", "philomiragia"])

    monkeypatch.setattr(http, "download", fake_download)
    monkeypatch.setattr(http, "verify_md5_manifest", lambda *a, **k: None)

    def build_params() -> MetadataParams:
        return MetadataParams(
            dataset="all",
            level="species",
            release="232.0",
            version="bac120",
            target_genus="Francisella",
            target_species="tularensis",
            nodownload=True,
        )

    def read_log() -> str:
        for h in logging.getLogger("repgenr").handlers:
            h.flush()
        return (wd / "repgenr.log").read_text(encoding="utf-8")

    monkeypatch.setitem(cli._RUN_STATE, "force", False)
    monkeypatch.setitem(cli._RUN_STATE, "log_level", logging.INFO)
    cli._run("metadata", wd, build_params, create=True)
    cli._run("metadata", wd, build_params, create=True)
    assert len(downloads) == 1
    assert read_log().count("already completed") == 1


def test_tree2tax_inputs_include_the_outgroup_and_segments_when_present(tmp_path: Path) -> None:
    """The rooting, --remove-outgroup and the segment rows of genomes_map.tsv
    depend on these files, so a change to them must rerun the stage."""
    ctx = WorkdirContext(tmp_path, create=True)
    ctx.tree_dir.mkdir(parents=True, exist_ok=True)
    (ctx.tree_dir / TREE_NWK).write_text("(a,b);\n", encoding="utf-8")
    params = SimpleNamespace(include_dereplicated=False)
    plain = _digests(ctx, "tree2tax", params)
    assert {plain["outgroup_accession.txt"], plain["outgroup"], plain["segments.tsv"]} == {ABSENT}
    (ctx.workdir / "outgroup_accession.txt").write_text("GCF_9\n", encoding="utf-8")
    (ctx.workdir / "outgroup").mkdir()
    (ctx.workdir / "outgroup" / "GCF_9.fasta").write_text(">o\nACGT\n", encoding="utf-8")
    (ctx.workdir / "segments.tsv").write_text("isolate\taccession\tsegment\tsegment_label\n")
    withfiles = _digests(ctx, "tree2tax", params)
    assert ABSENT not in {
        withfiles["outgroup_accession.txt"],
        withfiles["outgroup"],
        withfiles["segments.tsv"],
    }
    (ctx.workdir / "outgroup_accession.txt").write_text("GCF_8\n", encoding="utf-8")
    assert (
        _digests(ctx, "tree2tax", params)["outgroup_accession.txt"]
        != withfiles["outgroup_accession.txt"]
    )
