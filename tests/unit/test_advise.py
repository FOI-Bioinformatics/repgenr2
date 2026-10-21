"""phylo-advise: the recommendation rules, the evidence read from a workdir,
and the command's text and JSON output."""

from __future__ import annotations

import json
import logging
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.context import WorkdirContext
from repgenr.maskers.gubbins import genus_spread, species_spread
from repgenr.stages import advise as mod

_LOG = logging.getLogger("test-advise")
_runner = CliRunner()


def _ev(**kw) -> mod.Evidence:
    base = dict(
        workdir="/wd",
        all_genomes=False,
        n_genomes=100,
        n_species=1,
        n_genera=1,
        names_canonical=True,
        ani_min=None,
        ani_median=None,
        has_annotations=True,
        outgroup=True,
    )
    base.update(kw)
    return mod.Evidence(**base)


# --- rules -----------------------------------------------------------------------


def test_every_preset_has_flags_summary_and_annotation_need() -> None:
    assert set(mod.RECOMMENDATIONS) == {
        "outbreak",
        "species",
        "species-masked",
        "genus",
        "family",
        "quick",
    }
    for entry in mod.RECOMMENDATIONS.values():
        assert entry["flags"] and entry["summary"]
        assert isinstance(entry["needs_annotation"], bool)
    assert mod.RECOMMENDATIONS["genus"]["needs_annotation"]
    assert not mod.RECOMMENDATIONS["species"]["needs_annotation"]
    masked = mod.RECOMMENDATIONS["species-masked"]["flags"]
    assert masked == (*mod.RECOMMENDATIONS["species"]["flags"], "--mask", "gubbins")


def test_genus_and_family_rows() -> None:
    genus = mod.preset_flags("genus")
    assert genus == [
        "--msa-source",
        "coregene",
        "--coregene",
        "ppanggolin",
        "--core-threshold",
        "0.95",
        "--max-core-genes",
        "300",
        "--gene-rank",
        "informative",
        "--partitioned",
        "--treebuilder",
        "iqtree",
        "-B",
        "1000",
    ]
    # The family row swaps PPanGGOLiN for the GTDB markers, which imply protein.
    family = mod.preset_flags("family")
    assert family == [("markers" if f == "ppanggolin" else f) for f in genus]
    assert "--coregene-alphabet" not in family
    # Trimming stays an explicit choice.
    assert not any("--gene-trim" in mod.preset_flags(name) for name in mod.RECOMMENDATIONS)
    assert not any("planned" in e["summary"] for e in mod.RECOMMENDATIONS.values())


def test_expand_preset_rules() -> None:
    assert mod.expand_preset("quick") == ["--treebuilder", "skani"]
    # An explicit option replaces the preset's value, whichever spelling.
    out = mod.expand_preset("species", ["--bootstrap", "100"])
    assert out.count("-B") == 1 and out[-2:] == ["-B", "100"]
    # --aligner-arg adds to the preset's.
    out = mod.expand_preset("outbreak", ["--aligner-arg", "x=1"])
    assert "model=GTR+G" in out and out[-2:] == ["--aligner-arg", "x=1"]
    # Another tree builder drops the preset's tree-builder tuning.
    out = mod.expand_preset("genus", ["--treebuilder", "fasttree"])
    assert "--partitioned" not in out and out.count("--treebuilder") == 1
    out = mod.expand_preset("outbreak", ["--treebuilder", "fasttree"])
    assert "model=GTR+G" not in out
    # Unless the new builder reads it: raxmlng takes partitions and model=.
    assert "--partitioned" in mod.expand_preset("genus", ["--treebuilder", "raxmlng"])
    assert "model=GTR+G" in mod.expand_preset("outbreak", ["--treebuilder", "raxmlng"])
    assert not mod.builder_reads("auto", "--partitioned", None)
    assert mod.expand_preset("genus", ["--treebuilder", "iqtree"]).count("--partitioned") == 1


def test_outbreak_when_the_minimum_ani_is_at_least_99_5() -> None:
    rec = mod.advise(_ev(ani_min=99.6, ani_median=99.9, names_canonical=False))
    assert rec.preset == "outbreak" and rec.route == "snptype"
    assert rec.command[-2:] == ("--preset", "outbreak")
    assert "--snptyper" in rec.expanded and "ska2" in rec.expanded
    assert any("99.6" in r for r in rec.reasons)


def test_species_masked_for_a_small_close_species_set() -> None:
    rec = mod.advise(_ev(ani_min=98.0, ani_median=99.0))
    assert rec.preset == "species-masked"
    assert rec.command[-2:] == ("--preset", "species-masked")
    assert rec.expanded[-2:] == ("--mask", "gubbins")
    assert any("Gubbins" in r for r in rec.reasons)


@pytest.mark.parametrize(
    "kw",
    [
        {"ani_min": 96.0},  # too divergent for Gubbins
        {"ani_min": 98.0, "n_genomes": 201},  # too many genomes
        {"ani_min": None},  # divergence unknown
    ],
)
def test_species_without_masking(kw) -> None:
    rec = mod.advise(_ev(**kw))
    assert rec.preset == "species"
    assert "gubbins" not in rec.expanded


def test_species_from_ani_when_names_are_not_canonical() -> None:
    rec = mod.advise(_ev(names_canonical=False, n_species=0, n_genera=0, ani_min=96.0))
    assert rec.preset == "species"


def test_one_species_by_name_with_a_low_ani_says_so() -> None:
    rec = mod.advise(_ev(ani_min=90.0))
    assert rec.preset == "species"
    assert any("mislabelled" in r for r in rec.reasons)


def test_genus_for_several_species_or_mid_ani() -> None:
    assert mod.advise(_ev(n_species=4)).preset == "genus"
    rec = mod.advise(_ev(names_canonical=False, n_species=0, n_genera=0, ani_min=88.0))
    assert rec.preset == "genus" and rec.route == "coregene"
    assert "ppanggolin" in rec.expanded and "informative" in rec.expanded
    assert "--partitioned" in rec.expanded and "--coregene-alphabet" not in rec.expanded


def test_family_for_low_ani_or_several_genera(monkeypatch) -> None:
    monkeypatch.delenv("REPGENR_MARKERS_DB", raising=False)
    assert mod.advise(_ev(n_species=5, n_genera=2)).preset == "family"
    rec = mod.advise(_ev(names_canonical=False, n_species=0, n_genera=0, ani_min=78.0))
    assert rec.preset == "family"
    assert any("single-copy marker proteins" in r for r in rec.reasons)
    assert "markers" in rec.expanded and "--coregene-alphabet" not in rec.expanded
    assert "--partitioned" in rec.expanded
    # Without REPGENR_MARKERS_DB the reasons name the download command.
    assert any("repgenr markers-db download --output DIR" in r for r in rec.reasons)
    monkeypatch.setenv("REPGENR_MARKERS_DB", "/db")
    rec = mod.advise(_ev(n_species=5, n_genera=2))
    assert not any("markers-db download" in r for r in rec.reasons)


def test_unknown_pairs_count_as_below_80() -> None:
    """The Wolbachia set: one species by name, the known pairs above 80 %,
    and pairs that share fewer than 30 hashes: family."""
    ev = _ev(ani_min=81.5, ani_median=97.0, ani_unknown_pairs=27)
    rec = mod.advise(ev)
    assert rec.preset == "family"
    assert any(
        "27 pairs share fewer than 30 hashes (and under a tenth of the smaller sketch), "
        "too few to estimate an ANI; they are counted as below 80 %." in r
        for r in rec.reasons
    )
    assert mod.SOURMASH_K21_REASON in rec.reasons
    # Unknown pairs rule out the close presets whatever the known pairs give.
    rec = mod.advise(_ev(ani_min=99.8, names_canonical=False, ani_unknown_pairs=1))
    assert rec.preset == "family"
    # Every pair unknown: no minimum, still family.
    assert mod.advise(_ev(names_canonical=False, ani_unknown_pairs=3)).preset == "family"
    wolbachia = mod.advise(ev)
    text = mod.format_text(wolbachia, ev)
    assert "min 81.50 %, median 97.00 % over the known pairs; 27 pairs unknown" in text
    assert text.isascii()
    assert mod.to_json(wolbachia, ev)["evidence"]["ani_unknown_pairs"] == 27


@pytest.mark.parametrize(
    ("kw", "k21"),
    [
        ({"ani_min": 98.0}, False),  # one species
        ({"ani_min": 99.8}, False),  # outbreak
        ({"ani_min": 90.0}, True),  # below 95 %
        ({"n_species": 5, "n_genera": 2}, True),  # several genera by name
        ({"names_canonical": False, "ani_min": 88.0}, True),
        ({"ani_unknown_pairs": 2, "ani_min": 96.0}, True),
    ],
)
def test_k21_advice_for_genus_or_wider(kw, k21) -> None:
    rec = mod.advise(_ev(**kw))
    assert (mod.SOURMASH_K21_REASON in rec.reasons) is k21
    assert "--aligner-arg ksize=21" in mod.SOURMASH_K21_REASON
    assert "ksize=21" in mod.RECOMMENDATIONS["quick"]["summary"]


def test_quick_when_nothing_is_known() -> None:
    rec = mod.advise(_ev(names_canonical=False, n_species=0, n_genera=0, ani_min=None))
    assert rec.preset == "quick"
    assert rec.expanded[-2:] == ("--treebuilder", "skani")
    assert any("repgenr sketch" in r for r in rec.reasons)


def test_annotation_hint_for_the_coregene_route() -> None:
    rec = mod.advise(_ev(n_species=4, has_annotations=False, all_genomes=True))
    assert any("repgenr annotate -wd /wd --all-genomes" in r for r in rec.reasons)
    assert "--all-genomes" in rec.command
    rec = mod.advise(_ev(n_species=4, has_annotations=True))
    assert not any("repgenr annotate" in r for r in rec.reasons)
    # A route without annotation needs no hint.
    rec = mod.advise(_ev(ani_min=98.0, has_annotations=False))
    assert not any("repgenr annotate" in r for r in rec.reasons)


def test_midpoint_root_without_an_outgroup() -> None:
    rec = mod.advise(_ev(outgroup=False))
    assert rec.command[-2:] == ("--root", "midpoint")
    assert rec.expanded[-2:] == ("--root", "midpoint")
    assert any("--root midpoint" in r for r in rec.reasons)
    assert "--root" not in mod.advise(_ev(outgroup=True)).command
    assert "--no-outgroup" not in rec.command


def test_fasttree_above_the_iqtree_limit() -> None:
    rec = mod.advise(_ev(n_genomes=501, ani_min=99.8))
    assert rec.command[-6:] == (
        "--preset",
        "outbreak",
        "--treebuilder",
        "fasttree",
        "--aligner-arg",
        "gamma=true",
    )
    cmd = rec.expanded
    assert cmd.count("--treebuilder") == 1
    assert cmd[cmd.index("--treebuilder") + 1] == "fasttree"
    # FastTree takes its SH-like supports from the -B resamples.
    assert "model=GTR+G" not in cmd and cmd[cmd.index("-B") + 1] == "1000"
    assert any("--treebuilder fasttree" in r for r in rec.reasons)
    rec = mod.advise(_ev(n_genomes=500, ani_min=99.8))
    assert rec.expanded[rec.expanded.index("--treebuilder") + 1] == "iqtree"
    assert "--treebuilder" not in rec.command
    # A genus set: the swap drops the partitions as well; no IQ-TREE model to replace.
    rec = mod.advise(_ev(n_genomes=600, n_species=4))
    assert rec.command[-2:] == ("--treebuilder", "fasttree")
    assert "--partitioned" not in rec.expanded
    # The alignment-free tree builder is not swapped.
    rec = mod.advise(_ev(n_genomes=900, names_canonical=False, n_species=0, n_genera=0))
    assert rec.expanded[rec.expanded.index("--treebuilder") + 1] == "skani"
    assert "--treebuilder" not in rec.command


def test_command_starts_with_the_workdir() -> None:
    rec = mod.advise(_ev())
    assert rec.command[:4] == ("repgenr", "phylo", "-wd", "/wd")


def test_text_and_json_carry_the_same_recommendation() -> None:
    ev = _ev(ani_min=97.5, ani_median=99.0, n_genomes=40)
    rec = mod.advise(ev)
    text = mod.format_text(rec, ev)
    assert "species-masked" in text and "97.50" in text
    assert "repgenr phylo -wd /wd --preset species-masked" in text
    assert "Expands to:\n  repgenr phylo -wd /wd --msa-source snptype" in text
    assert text.isascii()
    data = mod.to_json(rec, ev)
    assert data["schema"] == mod.ADVISE_SCHEMA
    assert data["recommendation"]["preset"] == "species-masked"
    assert data["recommendation"]["command"][:2] == ["repgenr", "phylo"]
    assert data["recommendation"]["expanded_command"] == list(rec.expanded)
    assert data["evidence"]["ani_min"] == 97.5
    json.dumps(data)


# --- names -----------------------------------------------------------------------


def test_genus_spread_counts_genera_of_canonical_names() -> None:
    names = [
        "Fam_Gen_alpha_GCA_000000001.1",
        "Fam_Gen_beta_GCA_000000002.1",
        "Fam_Other_gamma_GCA_000000003.1",
        "Francisella_tularensis_LVS_1",
    ]
    assert species_spread(names) == 3
    assert genus_spread(names) == 2
    assert genus_spread(["sample1"]) == 0


# --- evidence from a workdir -------------------------------------------------------


def _workdir(tmp_path: Path, names: list[str], *, outgroup: str | None = None) -> Path:
    wd = tmp_path / "wd"
    reps = wd / "derep" / "representatives"
    reps.mkdir(parents=True)
    for n in names:
        (reps / f"{n}.fasta").write_text(">x\nACGT\n", encoding="utf-8")
    if outgroup is not None:
        (wd / "outgroup").mkdir()
        (wd / "outgroup" / f"{outgroup}.fasta").write_text(">o\nACGT\n", encoding="utf-8")
        acc = "_".join(outgroup.split("_")[3:])
        (wd / "outgroup_accession.txt").write_text(acc + "\n", encoding="utf-8")
    return wd


_ONE = ["Fam_Gen_sp_GCA_000000001.1", "Fam_Gen_sp_GCA_000000002.1", "Fam_Gen_sp_GCA_000000003.1"]


def test_collect_evidence_one_species_without_sketches(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, _ONE, outgroup="Fam_Out_x_GCA_000000009.1")
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 1, _LOG)
    assert (ev.n_genomes, ev.n_species, ev.n_genera) == (3, 1, 1)
    assert ev.names_canonical and ev.outgroup and not ev.has_annotations
    assert ev.ani_min is None and ev.ani_median is None
    assert any("sketch" in n for n in ev.notes)
    # Read-only: no manifest, no sketches, no log.
    assert not (wd / "manifest.sqlite").exists()
    assert sorted(p.name for p in wd.iterdir()) == ["derep", "outgroup", "outgroup_accession.txt"]


def test_collect_evidence_several_species_and_no_outgroup(tmp_path: Path) -> None:
    names = [
        "Fam_Gen_alpha_GCA_000000001.1",
        "Fam_Gen_beta_GCA_000000002.1",
        "Fam_Other_gamma_GCA_000000003.1",
    ]
    ev = mod.collect_evidence(
        WorkdirContext(_workdir(tmp_path, names), logger=_LOG), False, 1, _LOG
    )
    assert (ev.n_species, ev.n_genera, ev.outgroup) == (3, 2, False)
    assert mod.advise(ev).preset == "family"


def test_collect_evidence_non_canonical_names(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, ["sample1", "sample2"])
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 1, _LOG)
    assert not ev.names_canonical and ev.n_species == 0
    assert mod.advise(ev).preset == "quick"


def test_collect_evidence_reads_min_and_median_ani(tmp_path: Path, monkeypatch) -> None:
    wd = _workdir(tmp_path, _ONE)
    (wd / "sketches").mkdir()
    (wd / "manifest.sqlite").write_bytes(b"")
    reps = sorted((wd / "derep" / "representatives").iterdir())

    def fake_provider(ctx, logger, threads, consumer, *, write=True):
        assert write is False
        return lambda gs: {g: tmp_path / f"{g.name}.sig" for g in gs}

    matrix = np.array([[1.0, 0.99, 0.97], [0.99, 1.0, 0.98], [0.97, 0.98, 1.0]])
    seen: dict = {}

    def fake_ani(genomes, sketches, work_dir, threads, logger):
        seen["n"] = len(genomes)
        return [g.name for g in reps], matrix

    monkeypatch.setattr(mod, "workdir_provider", fake_provider)
    monkeypatch.setattr(mod, "_sourmash_missing", lambda: None)
    monkeypatch.setattr(mod, "ani_matrix", fake_ani)
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 2, _LOG)
    assert seen["n"] == 3
    assert ev.ani_min == pytest.approx(97.0)
    assert ev.ani_median == pytest.approx(98.0)
    assert mod.advise(ev).preset == "species-masked"


def test_collect_evidence_takes_min_and_median_over_known_pairs(
    tmp_path: Path, monkeypatch
) -> None:
    wd = _workdir(tmp_path, _ONE)
    (wd / "sketches").mkdir()
    (wd / "manifest.sqlite").write_bytes(b"")
    reps = sorted((wd / "derep" / "representatives").iterdir())
    matrix = np.array([[1.0, 0.99, np.nan], [0.99, 1.0, 0.97], [np.nan, 0.97, 1.0]])
    monkeypatch.setattr(
        mod,
        "workdir_provider",
        lambda ctx, logger, threads, consumer, *, write=True: (
            lambda gs: {g: tmp_path / f"{g.name}.sig" for g in gs}
        ),
    )
    monkeypatch.setattr(mod, "_sourmash_missing", lambda: None)
    monkeypatch.setattr(mod, "ani_matrix", lambda *a: ([g.name for g in reps], matrix))
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 1, _LOG)
    assert ev.ani_min == pytest.approx(97.0) and ev.ani_median == pytest.approx(98.0)
    assert ev.ani_unknown_pairs == 1
    assert mod.advise(ev).preset == "family"


def test_collect_evidence_skips_the_ani_beyond_the_dense_limit(tmp_path: Path, monkeypatch) -> None:
    wd = _workdir(tmp_path, _ONE)
    (wd / "sketches").mkdir()
    (wd / "manifest.sqlite").write_bytes(b"")
    monkeypatch.setattr(mod, "ANI_MAX_GENOMES", 2)
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 1, _LOG)
    assert ev.ani_min is None and any("above 2" in n for n in ev.notes)


def test_collect_evidence_sees_complete_annotations(tmp_path: Path, monkeypatch) -> None:
    wd = _workdir(tmp_path, _ONE)
    (wd / "annotations").mkdir()
    monkeypatch.setattr(mod, "_annotated", lambda directory, names: len(names) == 3)
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 1, _LOG)
    assert ev.has_annotations


def test_collect_evidence_all_genomes(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, _ONE[:1])
    (wd / "genomes").mkdir()
    for n in _ONE:
        (wd / "genomes" / f"{n}.fasta").write_text(">x\nACGT\n", encoding="utf-8")
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), True, 1, _LOG)
    assert ev.n_genomes == 3 and ev.all_genomes


def test_no_genomes_is_reported(tmp_path: Path) -> None:
    wd = tmp_path / "wd"
    wd.mkdir()
    ev = mod.collect_evidence(WorkdirContext(wd, logger=_LOG), False, 1, _LOG)
    assert ev.n_genomes == 0
    rec = mod.advise(ev)
    assert "0 genome(s)" in rec.reasons[0] and "run dereplicate first" in rec.reasons[0]
    # No pair to compare: no hint to sketch.
    assert not any("repgenr sketch" in r for r in rec.reasons)


def test_fewer_than_three_genomes_comes_first_and_keeps_the_command() -> None:
    rec = mod.advise(_ev(n_genomes=2, all_genomes=True, ani_min=99.8))
    assert rec.reasons[0].startswith("The set has 2 genome(s); phylo needs at least 3")
    assert "dereplicate" not in rec.reasons[0]
    assert rec.preset == "outbreak" and rec.command[:2] == ("repgenr", "phylo")
    assert not mod.advise(_ev(n_genomes=3)).reasons[0].startswith("The set has")


def test_iqtree_limit_is_read_from_the_registry() -> None:
    from repgenr.treebuilders.base import registry

    assert mod._iqtree_limit() == registry.get("iqtree").capabilities.recommended_max_genomes


# --- command ------------------------------------------------------------------------


def test_command_text_output(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, _ONE)
    result = _runner.invoke(app, ["phylo-advise", "-wd", str(wd), "-t", "1"])
    assert result.exit_code == 0, result.output
    assert "Preset: species" in result.stdout
    assert f"repgenr phylo -wd {wd} --preset species --root midpoint" in result.stdout
    assert f"repgenr phylo -wd {wd} --msa-source snptype" in result.stdout


def test_command_json_output(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, _ONE)
    result = _runner.invoke(app, ["phylo-advise", "-wd", str(wd), "--json"])
    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert data["schema"] == mod.ADVISE_SCHEMA
    assert data["evidence"]["n_genomes"] == 3
    assert data["recommendation"]["preset"] == "species"


def test_command_missing_workdir_exits_3(tmp_path: Path) -> None:
    result = _runner.invoke(app, ["phylo-advise", "-wd", str(tmp_path / "nope")])
    assert result.exit_code == 3


def test_recommendation_is_frozen() -> None:
    rec = mod.advise(_ev())
    with pytest.raises(AttributeError):
        rec.preset = "x"  # type: ignore[misc]
    assert replace(rec, preset="quick").preset == "quick"
