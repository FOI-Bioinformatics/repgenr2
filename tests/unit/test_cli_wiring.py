"""CLI smoke tests: option -> stage-parameter wiring (WP2-6).

The shared ``_run`` harness is replaced by a recorder in every command module,
so each command's ``build()`` (including its validation) runs for real and the
produced params object is asserted without touching a workdir or a stage.
"""

from __future__ import annotations

import pytest
from typer.testing import CliRunner

from repgenr.cli import cmd_bacterial, cmd_misc, cmd_phylo, cmd_run, cmd_viral
from repgenr.cli.main import app

_runner = CliRunner()


@pytest.fixture()
def dispatched(monkeypatch) -> list[tuple]:
    calls: list[tuple] = []

    def fake_run(stage, workdir, build, *, create=False):
        calls.append((stage, build(), create))

    for mod in (cmd_bacterial, cmd_viral, cmd_phylo, cmd_misc, cmd_run):
        monkeypatch.setattr(mod, "_run", fake_run)
    monkeypatch.setattr(cmd_run, "_preflight_tools", lambda *a, **k: None)
    return calls


def _invoke(dispatched, args: list[str]):
    result = _runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert len(dispatched) == 1
    return dispatched[0]


def test_metadata_wiring(dispatched, tmp_path) -> None:
    stage, params, create = _invoke(
        dispatched,
        [
            "metadata",
            "-wd",
            str(tmp_path),
            "--source",
            "api",
            "--level",
            "species",
            "-tg",
            "Francisella",
            "-ts",
            "tularensis",
            "--limit",
            "5",
            "--dataset",
            "rep",
        ],
    )
    assert stage == "metadata" and create is True
    assert (params.source, params.level) == ("api", "species")
    assert (params.target_genus, params.target_species) == ("Francisella", "tularensis")
    assert (params.limit, params.dataset) == (5, "rep")


def test_genome_wiring(dispatched, tmp_path) -> None:
    stage, params, create = _invoke(
        dispatched,
        [
            "genome",
            "-wd",
            str(tmp_path),
            "--accession-list-only",
            "--keep-files",
        ],
    )
    assert stage == "genome" and create is False
    assert params.accession_list_only is True and params.keep_files is True


def test_dereplicate_wiring(dispatched, tmp_path) -> None:
    stage, params, _ = _invoke(
        dispatched,
        [
            "dereplicate",
            "-wd",
            str(tmp_path),
            "--tool",
            "skder",
            "-sani",
            "0.98",
            "-pani",
            "0.85",
            "-af",
            "0.6",
            "-t",
            "4",
            "--process-size",
            "100",
            "--reduce",
            "species",
            "--target-reps",
            "10",
            "--keeper",
            "tool",
        ],
    )
    assert stage == "dereplicate"
    assert (params.tool, params.secondary_ani, params.primary_ani) == ("skder", 0.98, 0.85)
    assert (params.aligned_fraction, params.threads) == (0.6, 4)
    assert (params.process_size, params.reduce, params.target_reps) == (100, "species", 10)
    assert params.keeper == "tool"


def test_dereplicate_keeper_defaults_to_quality(dispatched, tmp_path) -> None:
    _, params, _ = _invoke(dispatched, ["dereplicate", "-wd", str(tmp_path), "--tool", "skder"])
    assert params.keeper == "quality"


def test_dereplicate_unknown_keeper_fails(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app, ["dereplicate", "-wd", str(tmp_path), "--tool", "skder", "--keeper", "nosuch"]
    )
    assert result.exit_code != 0
    assert dispatched == []


def test_dereplicate_unknown_tool_fails(dispatched, tmp_path) -> None:
    result = _runner.invoke(app, ["dereplicate", "-wd", str(tmp_path), "--tool", "nosuch"])
    assert result.exit_code != 0
    assert dispatched == []


def test_snptype_wiring(dispatched, tmp_path) -> None:
    stage, params, _ = _invoke(
        dispatched,
        [
            "snptype",
            "-wd",
            str(tmp_path),
            "--tool",
            "simple",
            "--mask",
            "gubbins",
            "--reference",
            "ref.fasta",
            "--all-genomes",
            "-t",
            "3",
        ],
    )
    assert stage == "snptype"
    assert (params.tool, params.mask, params.reference) == ("simple", "gubbins", "ref.fasta")
    assert params.all_genomes is True and params.threads == 3


def test_phylo_wiring(dispatched, tmp_path) -> None:
    stage, params, _ = _invoke(
        dispatched,
        [
            "phylo",
            "-wd",
            str(tmp_path),
            "--treebuilder",
            "mashtree",
            "--msa-source",
            "snptype",
            "--snptyper",
            "simple",
            "--bootstrap",
            "1000",
            "--no-outgroup",
            "--aligner-arg",
            "kmer=15",
        ],
    )
    assert stage == "phylo"
    assert (params.treebuilder, params.msa_source, params.snptyper) == (
        "mashtree",
        "snptype",
        "simple",
    )
    assert params.bootstrap == 1000 and params.no_outgroup is True
    assert params.extra == {"kmer": "15"}


def test_phylo_bad_aligner_arg_fails(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "phylo",
            "-wd",
            str(tmp_path),
            "--aligner-arg",
            "not-key-value",
        ],
    )
    assert result.exit_code != 0
    assert dispatched == []


def test_tree2tax_wiring(dispatched, tmp_path) -> None:
    stage, params, _ = _invoke(
        dispatched,
        [
            "tree2tax",
            "-wd",
            str(tmp_path),
            "--node-basename",
            "NODE",
            "--root-name",
            "myroot",
            "--remove-outgroup",
            "--include-dereplicated",
        ],
    )
    assert stage == "tree2tax"
    assert (params.node_basename, params.root_name) == ("NODE", "myroot")
    assert params.remove_outgroup is True and params.include_dereplicated is True


def test_vmetadata_wiring(dispatched, tmp_path) -> None:
    stage, params, create = _invoke(
        dispatched,
        [
            "vmetadata",
            "-wd",
            str(tmp_path),
            "--target",
            "adenoviridae",
            "--source",
            "ncbi_virus",
            "--complete-only",
            "--host",
            "homo sapiens",
            "--released-after",
            "01/31/2024",
        ],
    )
    assert stage == "vmetadata" and create is True
    assert (params.target, params.source) == ("adenoviridae", "ncbi_virus")
    assert params.complete_only is True and params.host == "homo sapiens"
    assert params.released_after == "01/31/2024"


def test_vmetadata_bad_source_fails(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "vmetadata",
            "-wd",
            str(tmp_path),
            "--target",
            "x",
            "--source",
            "nosuch",
        ],
    )
    assert result.exit_code != 0


def test_vgenome_wiring(dispatched, tmp_path) -> None:
    stage, params, _ = _invoke(
        dispatched,
        [
            "vgenome",
            "-wd",
            str(tmp_path),
            "-tg",
            "Mastadenovirus",
            "--length-range",
            "25000-35000",
            "--length-deviation",
            "15",
            "--group-segments",
            "--no-outgroup",
            "--discard",
            "partial,plasmid",
        ],
    )
    assert stage == "vgenome"
    assert params.target_genus == "Mastadenovirus"
    assert (params.length_range, params.length_deviation) == ("25000-35000", 15)
    assert params.group_segments is True and params.no_outgroup is True
    assert params.discard == "partial,plasmid"


def test_glance_wiring(dispatched, tmp_path) -> None:
    stage, params, _ = _invoke(
        dispatched,
        [
            "glance",
            "-wd",
            str(tmp_path),
            "-t",
            "2",
            "--plot-max",
            "0.9",
            "--plot-min",
            "0.1",
            "--keep-files",
        ],
    )
    assert stage == "glance"
    assert (params.threads, params.plot_max, params.plot_min) == (2, 0.9, 0.1)
    assert params.keep_files is True


# --- stateless data-channel steps ---------------------------------------------


@pytest.fixture()
def step_calls(monkeypatch) -> list:
    """Record the params object each data-channel step core receives."""
    import repgenr.stages.derep_steps as derep_steps
    import repgenr.stages.genome_steps as genome_steps
    import repgenr.stages.phylo as phylo_mod
    import repgenr.stages.tree2tax as tree2tax_mod

    calls: list = []

    def record(params, logger):
        calls.append(params)

    monkeypatch.setattr(genome_steps, "genome_fetch", record)
    monkeypatch.setattr(derep_steps, "dereplicate_chunk", record)
    monkeypatch.setattr(derep_steps, "dereplicate_merge", record)
    monkeypatch.setattr(phylo_mod, "phylo_build", record)
    monkeypatch.setattr(tree2tax_mod, "tree2tax_relations", record)
    return calls


def test_genome_fetch_step_wiring(step_calls, tmp_path) -> None:
    sel = tmp_path / "selection.tsv"
    sel.write_text("accession\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "genome-fetch",
            "--selection",
            str(sel),
            "-o",
            str(tmp_path / "out"),
            "--keep-files",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.selection_tsv == sel and params.keep_files is True


def test_dereplicate_chunk_step_wiring(step_calls, tmp_path) -> None:
    g = tmp_path / "g1.fasta"
    g.write_text(">g1\nACGT\n", encoding="utf-8")
    fofn = tmp_path / "genomes.fofn"
    fofn.write_text(f"{g}\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "dereplicate-chunk",
            "--genomes-fofn",
            str(fofn),
            "-o",
            str(tmp_path / "out"),
            "--tool",
            "galah",
            "-sani",
            "0.97",
            "-t",
            "2",
            "--virus",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert (params.tool, params.secondary_ani, params.threads) == ("galah", 0.97, 2)
    assert params.genomes == [g] and params.extra == {}  # galah ignores 'virus'


def test_dereplicate_chunk_step_injects_virus_for_reader(step_calls, tmp_path) -> None:
    g = tmp_path / "g1.fasta"
    g.write_text(">g1\nACGT\n", encoding="utf-8")
    fofn = tmp_path / "genomes.fofn"
    fofn.write_text(f"{g}\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "dereplicate-chunk",
            "--genomes-fofn",
            str(fofn),
            "-o",
            str(tmp_path / "out"),
            "--tool",
            "drep",
            "-t",
            "2",
            "--virus",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.tool == "drep"
    assert params.extra == {"virus": True}  # drep reads it


def test_dereplicate_chunk_step_wiring_selection_and_keeper(step_calls, tmp_path) -> None:
    g = tmp_path / "g1.fasta"
    g.write_text(">g1\nACGT\n", encoding="utf-8")
    fofn = tmp_path / "genomes.fofn"
    fofn.write_text(f"{g}\n", encoding="utf-8")
    sel = tmp_path / "selection.tsv"
    sel.write_text("accession\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "dereplicate-chunk",
            "--genomes-fofn",
            str(fofn),
            "-o",
            str(tmp_path / "out"),
            "--selection-tsv",
            str(sel),
            "--keeper",
            "tool",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.selection_tsv == sel and params.keeper == "tool"


def test_dereplicate_merge_step_requires_chunks(step_calls, tmp_path) -> None:
    result = _runner.invoke(app, ["dereplicate-merge", "-o", str(tmp_path / "out")])
    assert result.exit_code != 0
    assert step_calls == []


def test_dereplicate_merge_step_wiring(step_calls, tmp_path) -> None:
    chunk = tmp_path / "chunk0"
    chunk.mkdir()
    result = _runner.invoke(
        app,
        [
            "dereplicate-merge",
            "-o",
            str(tmp_path / "out"),
            "--chunk-dir",
            str(chunk),
            "--tool",
            "skder",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.chunk_dirs == [chunk] and params.tool == "skder"


def test_dereplicate_merge_step_wiring_selection_and_keeper(step_calls, tmp_path) -> None:
    chunk = tmp_path / "chunk0"
    chunk.mkdir()
    sel = tmp_path / "selection.tsv"
    sel.write_text("accession\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "dereplicate-merge",
            "-o",
            str(tmp_path / "out"),
            "--chunk-dir",
            str(chunk),
            "--tool",
            "skder",
            "--selection-tsv",
            str(sel),
            "--keeper",
            "tool",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.selection_tsv == sel and params.keeper == "tool"


def test_dereplicate_merge_step_rejects_bad_keeper(step_calls, tmp_path) -> None:
    chunk = tmp_path / "chunk0"
    chunk.mkdir()
    result = _runner.invoke(
        app,
        [
            "dereplicate-merge",
            "-o",
            str(tmp_path / "out"),
            "--chunk-dir",
            str(chunk),
            "--keeper",
            "bogus",
        ],
    )
    assert result.exit_code != 0
    assert step_calls == []


def test_phylo_build_step_wiring(step_calls, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "phylo-build",
            "--genomes-dir",
            str(tmp_path),
            "-o",
            str(tmp_path / "out"),
            "--treebuilder",
            "mashtree",
            "--no-outgroup",
            "-t",
            "3",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.phylo.treebuilder == "mashtree"
    assert params.phylo.no_outgroup is True and params.phylo.threads == 3
    assert params.genomes_dir == tmp_path


def test_tree2tax_relations_step_wiring(step_calls, tmp_path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text("(a,b);\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        [
            "tree2tax-relations",
            "--tree",
            str(tree),
            "-o",
            str(tmp_path / "out"),
            "--node-basename",
            "NODE",
            "--remove-outgroup",
        ],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.tree == tree
    assert params.node_basename == "NODE" and params.remove_outgroup is True
    assert params.no_outgroup is False


def test_tree2tax_relations_step_no_outgroup_wiring(step_calls, tmp_path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text("(a,b);\n", encoding="utf-8")
    result = _runner.invoke(
        app,
        ["tree2tax-relations", "--tree", str(tree), "-o", str(tmp_path / "out"), "--no-outgroup"],
    )
    assert result.exit_code == 0, result.output
    (params,) = step_calls
    assert params.no_outgroup is True


def test_run_viral_injects_virus_extra_only_when_accepted(dispatched, tmp_path) -> None:
    """--viral must not fingerprint-churn tools that ignore extra['virus']."""
    wd = str(tmp_path)
    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "--viral",
            "--target",
            "adeno",
            "--tool",
            "skder",
        ],
    )
    assert result.exit_code == 0, result.output
    derep = {s: p for s, p, _ in dispatched}["dereplicate"]
    assert derep.extra == {}  # skder ignores 'virus'
    dispatched.clear()

    result = _runner.invoke(
        app,
        [
            "run",
            "-wd",
            wd,
            "--viral",
            "--target",
            "adeno",
            "--tool",
            "drep",
        ],
    )
    assert result.exit_code == 0, result.output
    derep = {s: p for s, p, _ in dispatched}["dereplicate"]
    assert derep.extra == {"virus": True}  # drep reads it


def test_tool_arg_reaches_extra(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "dereplicate",
            "-wd",
            str(tmp_path),
            "--tool",
            "skder",
            "--tool-arg",
            "mode=greedy",
        ],
    )
    assert result.exit_code == 0, result.output
    params = {s: p for s, p, _ in dispatched}["dereplicate"]
    assert params.extra == {"mode": "greedy"}


def test_help_lists_registered_tools() -> None:
    import re

    from repgenr.dereplicators.base import registry as derep_registry

    result = _runner.invoke(app, ["dereplicate", "--help"])
    assert result.exit_code == 0
    # Rich wraps the help panel mid-token and may color it; compare against
    # the text with ANSI codes and all whitespace removed.
    # (box-drawing borders can interrupt a wrapped token, so keep word chars only)
    plain = re.sub(r"[^A-Za-z0-9/]", "", re.sub(r"\x1b\[[0-9;]*m", "", result.output))
    for name in derep_registry.names():
        assert name in plain


def test_phylo_mask_option_feeds_extra(dispatched, tmp_path) -> None:
    result = _runner.invoke(
        app,
        [
            "phylo",
            "-wd",
            str(tmp_path),
            "--msa-source",
            "snptype",
            "--snptyper",
            "simple",
            "--mask",
            "gubbins",
        ],
    )
    assert result.exit_code == 0, result.output
    params = {s: p for s, p, _ in dispatched}["phylo"]
    assert params.extra["mask"] == "gubbins"


def test_tree2tax_collapse_flags_wiring(dispatched, tmp_path) -> None:
    stage, params, _create = _invoke(
        dispatched,
        [
            "tree2tax",
            "-wd",
            str(tmp_path),
            "--collapse-support",
            "0.7",
            "--collapse-length",
            "0.001",
        ],
    )
    assert stage == "tree2tax"
    assert params.collapse_support == 0.7
    assert params.collapse_length == 0.001


def test_tree2tax_relations_collapse_flags_wiring(monkeypatch, tmp_path) -> None:
    from repgenr.stages import tree2tax as t2t_mod

    seen = {}

    def fake(params, logger):
        seen["params"] = params
        return (tmp_path / "a", tmp_path / "b")

    monkeypatch.setattr(t2t_mod, "tree2tax_relations", fake)
    tree = tmp_path / "tree.nwk"
    tree.write_text("(a,b);\n")
    result = _runner.invoke(
        app,
        [
            "tree2tax-relations",
            "--tree",
            str(tree),
            "-o",
            str(tmp_path / "out"),
            "--collapse-support",
            "0.7",
            "--collapse-length",
            "0.001",
        ],
    )
    assert result.exit_code == 0, result.output
    assert seen["params"].collapse_support == 0.7
    assert seen["params"].collapse_length == 0.001


def test_outgroup_treebuilder_help_names_accepted_values() -> None:
    import re

    from repgenr.viral._outgroup import distance_matrix_builders

    result = _runner.invoke(app, ["vgenome", "--help"])
    assert result.exit_code == 0
    plain = re.sub(r"[^A-Za-z0-9/]", "", re.sub(r"\x1b\[[0-9;]*m", "", result.output))
    for name in distance_matrix_builders():
        assert plain.count(name) >= 1
    assert "Accepted" in plain


def _help_of(command: str, option: str) -> str:
    import typer.main

    cmd = typer.main.get_command(app).commands[command]  # type: ignore[attr-defined]
    return next(p.help for p in cmd.params if option in getattr(p, "opts", ()))


def test_glance_tool_help_lists_only_comparison_capable_dereplicators() -> None:
    from repgenr.dereplicators.base import compare_supporters, registry

    text = _help_of("glance", "--tool")
    supporters = compare_supporters()
    assert supporters
    for name in registry.names():
        assert (name in text) == (name in supporters)


@pytest.mark.parametrize("option", ["-tf", "-tg", "-ts"])
def test_reads_taxon_help_says_most_specific_is_used(option) -> None:
    text = _help_of("reads", option)
    assert "most specific" in text
    assert "not combined" in text


def test_node_basename_help_mentions_hash_names() -> None:
    text = _help_of("tree2tax", "--node-basename")
    assert "hash" in text
    assert "Without it" in text


def test_metadata_path_help_says_release_and_version_still_required() -> None:
    for command in ("metadata", "run"):
        text = _help_of(command, "--metadata-path")
        assert "-r" in text
        assert "--gtdb-version" in text


def test_run_warns_early_when_collapse_support_meets_a_builder_without_supports(caplog) -> None:
    import logging

    from repgenr.cli.cmd_run import _warn_collapse_without_supports

    log = logging.getLogger("test-run-warn")
    with caplog.at_level(logging.WARNING, logger="test-run-warn"):
        _warn_collapse_without_supports("mashtree", 0.5, log)
        _warn_collapse_without_supports("iqtree", 0.5, log)
        _warn_collapse_without_supports("mashtree", None, log)
        _warn_collapse_without_supports("auto", 0.5, log)
    messages = [r.getMessage() for r in caplog.records]
    assert len(messages) == 1 and "mashtree" in messages[0]
