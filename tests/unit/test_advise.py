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


def test_outbreak_when_the_minimum_ani_is_at_least_99_5() -> None:
    rec = mod.advise(_ev(ani_min=99.6, ani_median=99.9, names_canonical=False))
    assert rec.preset == "outbreak" and rec.route == "snptype"
    assert "--snptyper" in rec.command and "ska2" in rec.command
    assert any("99.6" in r for r in rec.reasons)


def test_species_masked_for_a_small_close_species_set() -> None:
    rec = mod.advise(_ev(ani_min=98.0, ani_median=99.0))
    assert rec.preset == "species-masked"
    assert rec.command[-2:] == ("--mask", "gubbins")
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
    assert "gubbins" not in rec.command


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
    assert "ppanggolin" in rec.command and "informative" in rec.command


def test_family_for_low_ani_or_several_genera() -> None:
    assert mod.advise(_ev(n_species=5, n_genera=2)).preset == "family"
    rec = mod.advise(_ev(names_canonical=False, n_species=0, n_genera=0, ani_min=78.0))
    assert rec.preset == "family"
    assert any("protein" in r for r in rec.reasons)


def test_quick_when_nothing_is_known() -> None:
    rec = mod.advise(_ev(names_canonical=False, n_species=0, n_genera=0, ani_min=None))
    assert rec.preset == "quick"
    assert rec.command[-2:] == ("--treebuilder", "skani")
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
    assert any("--root midpoint" in r for r in rec.reasons)
    assert "--root" not in mod.advise(_ev(outgroup=True)).command
    assert "--no-outgroup" not in rec.command


def test_fasttree_above_the_iqtree_limit() -> None:
    rec = mod.advise(_ev(n_genomes=501, ani_min=99.8))
    cmd = rec.command
    assert cmd[cmd.index("--treebuilder") + 1] == "fasttree"
    assert "model=GTR+G" not in cmd and "-B" not in cmd
    assert any("--treebuilder fasttree" in r for r in rec.reasons)
    rec = mod.advise(_ev(n_genomes=500, ani_min=99.8))
    assert rec.command[rec.command.index("--treebuilder") + 1] == "iqtree"
    # The alignment-free tree builder is not swapped.
    rec = mod.advise(_ev(n_genomes=900, names_canonical=False, n_species=0, n_genera=0))
    assert rec.command[rec.command.index("--treebuilder") + 1] == "skani"


def test_command_starts_with_the_workdir() -> None:
    rec = mod.advise(_ev())
    assert rec.command[:4] == ("repgenr", "phylo", "-wd", "/wd")


def test_text_and_json_carry_the_same_recommendation() -> None:
    ev = _ev(ani_min=97.5, ani_median=99.0, n_genomes=40)
    rec = mod.advise(ev)
    text = mod.format_text(rec, ev)
    assert "species-masked" in text and "97.50" in text and "repgenr phylo -wd /wd" in text
    assert text.isascii()
    data = mod.to_json(rec, ev)
    assert data["schema"] == mod.ADVISE_SCHEMA
    assert data["recommendation"]["preset"] == "species-masked"
    assert data["recommendation"]["command"][:2] == ["repgenr", "phylo"]
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
    assert any("no genomes" in r.lower() for r in rec.reasons)


# --- command ------------------------------------------------------------------------


def test_command_text_output(tmp_path: Path) -> None:
    wd = _workdir(tmp_path, _ONE)
    result = _runner.invoke(app, ["phylo-advise", "-wd", str(wd), "-t", "1"])
    assert result.exit_code == 0, result.output
    assert "Preset: species" in result.stdout
    assert f"repgenr phylo -wd {wd}" in result.stdout
    assert "--root midpoint" in result.stdout


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
