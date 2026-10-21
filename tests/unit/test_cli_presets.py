"""--preset on phylo, phylo-build and run: the expansion of the rows of
stages/advise.py RECOMMENDATIONS, the override rule, the refusals, the
record and the MSA key."""

from __future__ import annotations

import dataclasses
import logging
import re
from pathlib import Path
from typing import Any

import pytest
import typer
from typer.testing import CliRunner

from repgenr.cli import cmd_phylo, cmd_run
from repgenr.cli.base import _stage_fingerprint
from repgenr.cli.main import app
from repgenr.cli.param_builders import PRESET_PARAMS
from repgenr.core.context import WorkdirContext
from repgenr.core.plugins import ToolCapabilities
from repgenr.stages import advise
from repgenr.stages import phylo as phylo_stage
from repgenr.stages.annotate import AnnotateStageParams
from repgenr.stages.phylo import PhyloParams, _msa_key
from repgenr.treebuilders.base import InputKind, TreeBuilder
from repgenr.treebuilders.base import registry as tb_registry

_runner = CliRunner()
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_LOG = logging.getLogger("test-cli-presets")
PRESETS = sorted(advise.RECOMMENDATIONS)


def _text(result) -> str:
    return " ".join(_ANSI.sub("", result.output + str(result.exception or "")).split())


@pytest.fixture
def built(monkeypatch) -> list[Any]:
    """The params objects the commands build instead of running a stage."""
    objects: list[Any] = []

    def fake_run(stage, workdir, build, *, create=False):
        objects.append(build())

    monkeypatch.setattr(cmd_phylo, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    monkeypatch.setattr(phylo_stage, "phylo_build", lambda params, *a, **k: objects.append(params))
    return objects


def _phylo_args(tmp_path: Path) -> list[str]:
    return ["phylo", "-wd", str(tmp_path / "wd")]


def _run_args(tmp_path: Path) -> list[str]:
    return ["run", "-wd", str(tmp_path / "wd"), "-l", "genus"]


def _build_args(tmp_path: Path) -> list[str]:
    genomes = tmp_path / "genomes"
    genomes.mkdir(exist_ok=True)
    (genomes / "a.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    return ["phylo-build", "--genomes-dir", str(genomes), "-o", str(tmp_path / "out")]


COMMANDS = {"phylo": _phylo_args, "phylo-build": _build_args, "run": _run_args}


def _invoke(objects: list[Any], argv: list[str]) -> PhyloParams:
    """The PhyloParams the command built for ``argv``."""
    objects.clear()
    result = _runner.invoke(app, argv)
    assert result.exit_code == 0, _text(result)
    found = [o.phylo if hasattr(o, "phylo") else o for o in objects]
    return [o for o in found if isinstance(o, PhyloParams)][-1]


def _unlabelled(params: PhyloParams) -> PhyloParams:
    return dataclasses.replace(params, preset=None)


# --- expansion ---------------------------------------------------------------------


@pytest.mark.parametrize("command", sorted(COMMANDS))
@pytest.mark.parametrize("preset", PRESETS)
def test_preset_expands_to_its_documented_flags(built, tmp_path, command, preset) -> None:
    base = COMMANDS[command](tmp_path)
    with_preset = _invoke(built, [*base, "--preset", preset])
    written_out = _invoke(built, [*base, *advise.preset_flags(preset)])
    assert with_preset.preset == preset and written_out.preset is None
    assert _unlabelled(with_preset) == written_out


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_genus_and_family_presets_set_the_core_gene_route(built, tmp_path, command) -> None:
    params = _invoke(built, [*COMMANDS[command](tmp_path), "--preset", "genus"])
    assert (params.msa_source, params.coregene, params.treebuilder) == (
        "coregene",
        "ppanggolin",
        "iqtree",
    )
    assert params.bootstrap == 1000
    assert params.extra == {
        "core_threshold": 0.95,
        "max_core_genes": 300,
        "gene_rank": "informative",
        "partitioned": True,
    }
    params = _invoke(built, [*COMMANDS[command](tmp_path), "--preset", "family"])
    assert params.extra["coregene_alphabet"] == "protein"
    assert "gene_trim" not in params.extra


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_snptype_presets(built, tmp_path, command) -> None:
    base = COMMANDS[command](tmp_path)
    params = _invoke(built, [*base, "--preset", "outbreak"])
    assert (params.msa_source, params.snptyper) == ("snptype", "ska2")
    assert params.extra == {"model": "GTR+G"}
    params = _invoke(built, [*base, "--preset", "species-masked"])
    assert (params.msa_source, params.snptyper, params.bootstrap) == ("snptype", "simple", 1000)
    assert params.extra == {"core_fraction": 0.95, "mask": "gubbins"}
    params = _invoke(built, [*base, "--preset", "quick"])
    assert params.treebuilder == "skani" and params.msa_source == "aligner"


def test_run_inserts_annotate_for_the_core_gene_presets(built, tmp_path) -> None:
    for preset, annotated in (("genus", True), ("family", True), ("species", False)):
        _invoke(built, [*_run_args(tmp_path), "--preset", preset])
        assert any(isinstance(o, AnnotateStageParams) for o in built) is annotated, preset


def test_run_dry_run_shows_the_preset(tmp_path) -> None:
    result = _runner.invoke(app, [*_run_args(tmp_path), "--preset", "genus", "--dry-run"])
    assert result.exit_code == 0, _text(result)
    assert "- annotate" in result.output
    assert "preset=genus, treebuilder=iqtree, coregene (annotate first)" in _text(result)


def test_every_preset_flag_is_a_preset_parameter_of_each_command() -> None:
    click_app = typer.main.get_command(app)
    for command in COMMANDS:
        by_opt = {
            opt: p.name
            for p in click_app.commands[command].params  # type: ignore[attr-defined]
            for opt in getattr(p, "opts", ())
        }
        assert "--preset" in by_opt, command
        for preset in PRESETS:
            for opt, _value in advise.flag_items(advise.preset_flags(preset)):
                assert by_opt.get(opt) in PRESET_PARAMS, (command, preset, opt)


# --- the override rule ----------------------------------------------------------------


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_explicit_flags_override_the_preset(built, tmp_path, command) -> None:
    base = [*COMMANDS[command](tmp_path), "--preset", "species"]
    params = _invoke(built, [*base, "-B", "100", "--core-fraction", "0.9"])
    assert params.bootstrap == 100 and params.extra["core_fraction"] == 0.9
    assert params.snptyper == "simple" and params.preset == "species"
    # Given on the command line at its default value, it still overrides.
    assert _invoke(built, [*base, "--bootstrap", "0"]).bootstrap == 0
    assert _invoke(built, [*base, "--snptyper", "simple", "--treebuilder", "iqtree"]).bootstrap == (
        1000
    )


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_aligner_args_add_to_the_preset(built, tmp_path, command) -> None:
    base = [*COMMANDS[command](tmp_path), "--preset", "outbreak"]
    params = _invoke(built, [*base, "--aligner-arg", "seed=7"])
    assert params.extra == {"model": "GTR+G", "seed": "7"}
    params = _invoke(built, [*base, "--aligner-arg", "model=GTR+G+I"])
    assert params.extra == {"model": "GTR+G+I"}


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_another_tree_builder_drops_the_preset_tuning(built, tmp_path, command) -> None:
    base = COMMANDS[command](tmp_path)
    params = _invoke(built, [*base, "--preset", "outbreak", "--treebuilder", "fasttree"])
    assert params.treebuilder == "fasttree" and params.extra == {}
    assert params.bootstrap == 1000
    # raxmlng reads the IQ-TREE model key and partitions, so both stay.
    params = _invoke(built, [*base, "--preset", "outbreak", "--treebuilder", "raxmlng"])
    assert params.extra == {"model": "GTR+G"}
    params = _invoke(built, [*base, "--preset", "genus", "--treebuilder", "raxmlng"])
    assert params.extra["partitioned"] is True and params.extra["max_core_genes"] == 300
    params = _invoke(built, [*base, "--preset", "genus", "--treebuilder", "fasttree"])
    assert "partitioned" not in params.extra and params.extra["max_core_genes"] == 300
    # The preset's own tree builder, given explicitly, keeps the tuning.
    params = _invoke(built, [*base, "--preset", "genus", "--treebuilder", "iqtree"])
    assert params.extra["partitioned"] is True


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_another_msa_source_is_refused(built, tmp_path, command) -> None:
    argv = [*COMMANDS[command](tmp_path), "--preset", "species", "--msa-source", "aligner"]
    result = _runner.invoke(app, argv)
    assert result.exit_code != 0
    assert "--preset species sets --msa-source snptype" in _text(result)
    # The preset's own source is accepted; quick sets no source.
    _invoke(built, [*COMMANDS[command](tmp_path), "--preset", "species", "--msa-source", "snptype"])
    params = _invoke(
        built, [*COMMANDS[command](tmp_path), "--preset", "quick", "--msa-source", "snptype"]
    )
    assert params.msa_source == "snptype" and params.treebuilder == "skani"


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_unknown_preset_is_refused(built, tmp_path, command) -> None:
    result = _runner.invoke(app, [*COMMANDS[command](tmp_path), "--preset", "strain"])
    assert result.exit_code != 0
    text = _text(result)
    assert "--preset" in text and "outbreak" in text
    assert not any(isinstance(o, PhyloParams) for o in built)


@pytest.mark.parametrize("command", sorted(COMMANDS))
def test_unknown_preset_exits_with_the_usage_code(tmp_path, command) -> None:
    """Through the real _run (and stage_errors) rather than the recorder."""
    from repgenr.core.errors import UserInputError

    (tmp_path / "wd").mkdir()
    result = _runner.invoke(app, [*COMMANDS[command](tmp_path), "--preset", "strain"])
    assert result.exit_code == UserInputError.exit_code == 2, _text(result)
    assert "Invalid --preset 'strain'" in _text(result)
    assert not (tmp_path / "wd" / "tree").exists()


# --- phylo-advise ------------------------------------------------------------------------


def _ev(**kw) -> advise.Evidence:
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
    return advise.Evidence(**base)


@pytest.mark.parametrize(
    "kw",
    [
        {"ani_min": 99.8},
        {"ani_min": 99.8, "n_genomes": 501, "outgroup": False},
        {"ani_min": 98.0},
        {"n_species": 4, "n_genomes": 600},
        {"n_species": 5, "n_genera": 2},
        {"names_canonical": False, "n_species": 0, "n_genera": 0},
    ],
)
def test_advised_command_and_its_expansion_build_the_same_params(built, tmp_path, kw) -> None:
    rec = advise.advise(_ev(**kw))
    wd = str(tmp_path / "wd")
    with_preset = _invoke(built, ["phylo", "-wd", wd, *rec.command[4:]])
    written_out = _invoke(built, ["phylo", "-wd", wd, *rec.expanded[4:]])
    assert with_preset.preset == rec.preset
    assert _unlabelled(with_preset) == written_out


# --- record, MSA key and fingerprint ----------------------------------------------------


def test_msa_key_and_fingerprint_ignore_the_preset_name(tmp_path) -> None:
    genome = tmp_path / "a.fasta"
    genome.write_text(">a\nACGT\n", encoding="utf-8")
    plain = PhyloParams(msa_source="snptype", bootstrap=1000, extra={"core_fraction": 0.95})
    labelled = dataclasses.replace(plain, preset="species")
    assert _msa_key([genome], None, plain) == _msa_key([genome], None, labelled)
    assert _stage_fingerprint("phylo", plain, {}, {}) == _stage_fingerprint(
        "phylo", labelled, {}, {}
    )


class _GenomesTree(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_preset")
    input_kind = InputKind.GENOMES

    def preflight(self):
        return {"faketree_preset": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        out_dir.mkdir(parents=True, exist_ok=True)
        tree = out_dir / "tree.nwk"
        names = [Path(g).stem for g in msa_or_genomes]
        tree.write_text("(" + ",".join(names) + ");\n", encoding="utf-8")
        return tree


def test_record_carries_the_preset(workdir: Path, register_tool) -> None:
    register_tool(tb_registry, "faketree_preset", _GenomesTree)
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        (reps / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    ctx = WorkdirContext(workdir, create=True, logger=_LOG)
    phylo_stage.run(ctx, PhyloParams(treebuilder="faketree_preset", no_outgroup=True))
    assert ctx.config.stages["phylo"].params["preset"] is None
    phylo_stage.run(
        ctx, PhyloParams(treebuilder="faketree_preset", no_outgroup=True, preset="quick")
    )
    assert ctx.config.stages["phylo"].params["preset"] == "quick"


def test_nextflow_knows_the_core_gene_presets() -> None:
    """nextflow/subworkflows/local/annotate_scatter.nf runs ANNOTATE for these."""
    source = (
        Path(__file__).resolve().parents[2]
        / "nextflow"
        / "subworkflows"
        / "local"
        / "annotate_scatter.nf"
    ).read_text(encoding="utf-8")
    listed = re.search(r"return preset in \[([^\]]*)\]", source)
    assert listed is not None
    names = {n.strip().strip("'\"") for n in listed.group(1).split(",")}
    assert names == {n for n, e in advise.RECOMMENDATIONS.items() if e["route"] == "coregene"}
