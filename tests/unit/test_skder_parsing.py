"""Offline tests for skDER output parsing + edge-based membership."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.dereplicators.skder import _parse_skder_output

_LOG = logging.getLogger("test")


def _make_skder_out(tmp_path: Path, edges: str) -> Path:
    out = tmp_path / "skder_out"
    reps = out / "Dereplicated_Representative_Genomes"
    reps.mkdir(parents=True)
    for name in ("repA.fasta", "repB.fasta"):
        (reps / name).write_text(">x\nACGT\n")
    (out / "Skani_Triangle_Edge_Output.txt").write_text(edges)
    return out


def test_membership_from_edges(tmp_path: Path) -> None:
    # memX is closest (highest ANI) to repA; memY only similar to repB
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t99.5\t90\t88\trA\tmX\n"
        "/g/repB.fasta\t/g/memX.fasta\t97.0\t80\t80\trB\tmX\n"
        "/g/memY.fasta\t/g/repB.fasta\t99.2\t85\t85\tmY\trB\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX", "memY")]

    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)

    assert {p.name for p in result.representatives} == {"repA.fasta", "repB.fasta"}
    # memX assigned to repA (99.5 > 97.0); memY assigned to repB
    assert "memX.fasta" in result.clusters["repA.fasta"]
    assert "memY.fasta" in result.clusters["repB.fasta"]
    assert result.genome_status["memX.fasta"] == "contained"
    assert result.genome_status["repA.fasta"] == "representative"


def test_edges_below_cutoff_not_clustered(tmp_path: Path) -> None:
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t90.0\t90\t90\trA\tmX\n"  # ANI below 99 cutoff
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX")]
    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)
    # memX not assigned to any cluster, but still marked contained
    assert result.clusters["repA.fasta"] == []
    assert result.genome_status["memX.fasta"] == "contained"


_FAKE_SKDER = """\
import os, pathlib, sys
args = sys.argv[1:]
pathlib.Path(os.environ["FAKE_SKDER_ARGV"]).write_text("\\n".join(args))
gdir = args[args.index("-g") + 1]
out = pathlib.Path(args[args.index("-o") + 1])
names = sorted(os.listdir(gdir))
# Every staged entry must resolve to an input genome (no dangling link).
assert all(os.path.isfile(os.path.join(gdir, n)) for n in names), "dangling entry"
pathlib.Path(os.environ["FAKE_SKDER_LISTED"]).write_text("\\n".join(names))
reps = out / "Dereplicated_Representative_Genomes"
reps.mkdir(parents=True)
(reps / names[0]).write_text(">x\\nACGT\\n")
with open(out / "Skani_Triangle_Edge_Output.txt", "w") as fo:
    fo.write("Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\n")
    for n in names[1:]:
        fo.write(f"{gdir}/{names[0]}\\t{gdir}/{n}\\t99.9\\t90\\t90\\n")
"""


def _install_fake_skder(tmp_path: Path, monkeypatch) -> tuple[Path, Path]:
    import os
    import sys

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "skder"
    exe.write_text(f"#!{sys.executable}\n{_FAKE_SKDER}")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ.get('PATH', '')}")
    argv_log, listed = tmp_path / "argv.txt", tmp_path / "listed.txt"
    monkeypatch.setenv("FAKE_SKDER_ARGV", str(argv_log))
    monkeypatch.setenv("FAKE_SKDER_LISTED", str(listed))
    return argv_log, listed


def test_skder_receives_one_staged_directory_for_10000_genomes(tmp_path, monkeypatch) -> None:
    """The genome count does not reach argv: one -g token, the results map back."""
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.skder import SkderDereplicator

    argv_log, listed = _install_fake_skder(tmp_path, monkeypatch)
    src = tmp_path / "genomes"
    src.mkdir()
    genomes = []
    for i in range(10000):
        g = src / f"Fam_Gen_sp_GCA_{i:09d}.1.fasta"
        g.touch()
        genomes.append(g)
    # An AppleDouble companion beside the genomes is not part of the set.
    (src / f"._{genomes[0].name}").touch()
    params = DerepParams(primary_ani=0.9, secondary_ani=0.999, threads=1)

    result = SkderDereplicator().dereplicate(genomes, tmp_path / "out", params, _LOG)

    args = argv_log.read_text().split("\n")
    # -g DIR -o OUT -i -f -c -d with their values, plus the screen the cutoff 99.9 sets
    assert len(args) == 13
    assert args[-1] == "--skani-triangle-parameters=-s 98.9"
    staged = Path(args[args.index("-g") + 1])
    assert args.count("-g") == 1 and not staged.name.endswith(".fasta")
    names = listed.read_text().split("\n")
    assert names == sorted(g.name for g in genomes)
    assert [p.name for p in result.representatives] == [genomes[0].name]
    assert len(result.clusters[genomes[0].name]) == 9999
    assert set(result.genome_status) == {g.name for g in genomes}
    assert not staged.exists(), "the staging directory is scratch and is removed"


def test_staged_genome_dir_links_by_basename(tmp_path: Path) -> None:
    import os

    from repgenr.dereplicators.skder import stage_genome_dir

    a = tmp_path / "a" / "Fam_Gen_sp_X1.fasta.gz"
    a.parent.mkdir()
    a.write_bytes(b"gz")
    staged = stage_genome_dir([a], tmp_path / "staged")
    entry = staged / a.name
    assert entry.is_symlink() and os.readlink(entry) == os.path.abspath(a)
    assert [p.name for p in staged.iterdir()] == [a.name]


def test_skder_refuses_two_genomes_with_one_basename(tmp_path, monkeypatch) -> None:
    import pytest

    from repgenr.core.errors import UserInputError
    from repgenr.dereplicators import skder as skder_mod
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.skder import SkderDereplicator

    def must_not_run(caps, cmd, **kwargs):
        raise AssertionError("skDER must not run")

    monkeypatch.setattr(skder_mod, "run_tool", must_not_run)
    first, second = tmp_path / "a" / "x.fasta", tmp_path / "b" / "x.fasta"
    for g in (first, second):
        g.parent.mkdir()
        g.write_text(">x\nACGT\n")
    params = DerepParams(primary_ani=0.9, secondary_ani=0.99, threads=1)
    with pytest.raises(UserInputError) as exc:
        SkderDereplicator().dereplicate([first, second], tmp_path / "out", params, _LOG)
    assert str(first) in str(exc.value) and str(second) in str(exc.value)


def test_partial_genome_joins_by_its_own_aligned_fraction(tmp_path: Path) -> None:
    # A 40 percent fragment of a repA-like genome: 99 percent of the fragment
    # aligns, 40 percent of repA does. skDER covers the fragment by repA, so the
    # cutoff applies to the member's aligned fraction, not the smaller of the two.
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/frag.fasta\t/g/repA.fasta\t99.6\t99.3\t39.7\tf\trA\n"
        "/g/repB.fasta\t/g/frag.fasta\t95.0\t39.7\t99.2\trB\tf\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "frag")]
    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)
    assert result.clusters["repA.fasta"] == ["frag.fasta"]
    assert result.genome_status["frag.fasta"] == "contained"


def test_representative_aligned_fraction_does_not_admit_a_member(tmp_path: Path) -> None:
    # The member's own aligned fraction (30) is below the cutoff; only the
    # representative's side (95) passes. The edge does not assign the member by
    # the cutoff rule, so it is placed under its closest representative at or
    # above the ANI cutoff, with a warning naming it.
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t99.5\t95\t30\trA\tmX\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX")]
    log = logging.getLogger("test.skder.nearest")
    records: list[logging.LogRecord] = []
    handler = logging.Handler()
    handler.emit = records.append  # type: ignore[method-assign]
    log.addHandler(handler)
    try:
        result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=log)
    finally:
        log.removeHandler(handler)
    assert result.clusters["repA.fasta"] == ["memX.fasta"]
    assert any("memX.fasta" in r.getMessage() for r in records if r.levelno == logging.WARNING)


def test_skder_refuses_ani_below_80_percent_before_running(tmp_path, monkeypatch) -> None:
    import pytest

    from repgenr.core.errors import UserInputError
    from repgenr.dereplicators import skder as skder_mod
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.skder import SkderDereplicator

    def must_not_run(caps, cmd, **kwargs):
        raise AssertionError("skDER was started")

    monkeypatch.setattr(skder_mod, "run_tool", must_not_run)
    params = DerepParams(secondary_ani=0.78, threads=1)
    with pytest.raises(UserInputError, match="80 percent"):
        SkderDereplicator().dereplicate([tmp_path / "g.fasta"], tmp_path / "out", params, _LOG)


def test_many_placements_are_one_line_capped_on_the_console(tmp_path: Path) -> None:
    # Twelve members pass only the ANI cutoff: one warning, whose console form
    # names the first few and counts the rest, while the full text lists all.
    rows = "".join(f"/g/repA.fasta\t/g/mem{i:02d}.fasta\t99.5\t95\t30\trA\tm\n" for i in range(12))
    edges = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tR\tQ\n" + rows
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path("/g/repA.fasta"), Path("/g/repB.fasta")]
    genomes += [Path(f"/g/mem{i:02d}.fasta") for i in range(12)]
    log = logging.getLogger("test.skder.many")
    records: list[logging.LogRecord] = []
    handler = logging.Handler()
    handler.emit = records.append  # type: ignore[method-assign]
    log.addHandler(handler)
    try:
        result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=log)
    finally:
        log.removeHandler(handler)
    assert len(result.clusters["repA.fasta"]) == 12
    warnings = [r for r in records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert all(f"mem{i:02d}.fasta" in warnings[0].getMessage() for i in range(12))
    console = warnings[0].repgenr_console  # type: ignore[attr-defined]
    assert "and 7 more" in console and "mem11.fasta" not in console


def _skder_argv(tmp_path, monkeypatch, secondary_ani: float, extra: dict) -> list[str]:
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.skder import SkderDereplicator

    argv_log, _listed = _install_fake_skder(tmp_path, monkeypatch)
    src = tmp_path / "genomes"
    src.mkdir(exist_ok=True)
    genomes = [src / "Fam_Gen_sp_GCA_000000001.1.fasta", src / "Fam_Gen_sp_GCA_000000002.1.fasta"]
    for g in genomes:
        g.touch()
    params = DerepParams(primary_ani=0.9, secondary_ani=secondary_ani, threads=1, extra=extra)
    SkderDereplicator().dereplicate(genomes, tmp_path / "out", params, _LOG)
    return argv_log.read_text().split("\n")


def test_screen_defaults_to_one_below_a_high_cutoff(tmp_path, monkeypatch) -> None:
    """At 99.95 the screen is 98.95, so pairs below that skip the alignment."""
    args = _skder_argv(tmp_path, monkeypatch, 0.9995, {})
    assert args[-1] == "--skani-triangle-parameters=-s 98.95"


def test_screen_is_left_to_skder_below_99(tmp_path, monkeypatch) -> None:
    """Below a 99 percent cutoff the adapter keeps skDER's own screen (cutoff - 10)."""
    args = _skder_argv(tmp_path, monkeypatch, 0.95, {})
    assert not any(a.startswith("--skani-triangle-parameters") for a in args)


def test_screen_tool_arg_overrides_the_default(tmp_path, monkeypatch) -> None:
    args = _skder_argv(tmp_path, monkeypatch, 0.9995, {"screen": "90"})
    assert args[-1] == "--skani-triangle-parameters=-s 90"


@pytest.mark.parametrize(
    ("screen", "fragment"),
    [("99.96", "above the ANI cutoff"), ("abc", "ANI percentage"), ("100", "between 0 and 100")],
)
def test_screen_tool_arg_is_validated(screen: str, fragment: str) -> None:
    from repgenr.core.errors import UserInputError
    from repgenr.dereplicators.skder import resolve_screen

    with pytest.raises(UserInputError, match=fragment):
        resolve_screen(screen, 99.95)


def test_resolve_screen_values() -> None:
    from repgenr.dereplicators.skder import resolve_screen

    assert resolve_screen(None, 99.95) == 98.95
    assert resolve_screen(None, 99.0) == 98.0
    assert resolve_screen(None, 98.9) is None
    assert resolve_screen("95", 99.95) == 95.0
    assert resolve_screen("skder", 99.95) is None


def test_screen_above_a_lower_call_cutoff_falls_back_and_warns_once(caplog) -> None:
    """A --pre-* chunk or a --target-reps step runs below the final cutoff by design."""
    from repgenr.dereplicators import skder as skder_mod
    from repgenr.dereplicators.skder import resolve_screen

    skder_mod._screen_fallbacks_logged.clear()
    log = logging.getLogger("test_screen_fallback")
    with caplog.at_level(logging.WARNING, logger="test_screen_fallback"):
        assert resolve_screen("99.9", 99.0, log) == 98.0
        assert resolve_screen("99.9", 99.0, log) == 98.0
        assert resolve_screen("99.9", 90.0, log) is None
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 2
    assert "using a screen of 98" in warnings[0]
    assert "using skDER's own screen" in warnings[1]


def test_stage_precheck_refuses_a_screen_above_the_final_cutoff(tmp_path) -> None:
    from repgenr.core.context import WorkdirContext
    from repgenr.core.errors import UserInputError
    from repgenr.stages.dereplicate import DereplicateParams, precheck

    wd = tmp_path / "wd"
    (wd / "genomes").mkdir(parents=True)
    (wd / "genomes" / "Fam_Gen_sp_GCA_000000001.1.fasta").write_text(">x\nACGT\n")
    ctx = WorkdirContext(wd, create=True)
    params = DereplicateParams(
        tool="skder", secondary_ani=0.9995, extra={"screen": "99.96"}, allow_incomplete=True
    )
    with pytest.raises(UserInputError, match="above the ANI cutoff"):
        precheck(ctx, params)


def test_screen_argument_survives_a_parser_like_skders(tmp_path, monkeypatch) -> None:
    """skDER declares -p/--skani-triangle-parameters and a store-true -s; the one-token
    equals form must reach -p as its value (a separate "-s 98.95" token did not)."""
    import argparse

    args = _skder_argv(tmp_path, monkeypatch, 0.9995, {})
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", nargs="+")
    parser.add_argument("-o")
    parser.add_argument("-i")
    parser.add_argument("-f")
    parser.add_argument("-c")
    parser.add_argument("-d")
    parser.add_argument("-p", "--skani-triangle-parameters", default="-s X")
    parser.add_argument("-s", "--sanity-check", action="store_true")
    parsed = parser.parse_args(args)
    assert parsed.skani_triangle_parameters == "-s 98.95"
    assert parsed.sanity_check is False


def test_stage_precheck_checks_the_form_under_target_reps(tmp_path) -> None:
    from repgenr.core.context import WorkdirContext
    from repgenr.core.errors import UserInputError
    from repgenr.stages.dereplicate import DereplicateParams, precheck

    wd = tmp_path / "wd"
    (wd / "genomes").mkdir(parents=True)
    (wd / "genomes" / "Fam_Gen_sp_GCA_000000001.1.fasta").write_text(">x\nACGT\n")
    ctx = WorkdirContext(wd, create=True)
    with pytest.raises(UserInputError, match="ANI percentage"):
        precheck(
            ctx,
            DereplicateParams(
                tool="skder", target_reps=5, extra={"screen": "abc"}, allow_incomplete=True
            ),
        )
    # A value above the final cutoff is allowed there: the search sets the cutoff.
    precheck(
        ctx,
        DereplicateParams(
            tool="skder", target_reps=5, extra={"screen": "99.96"}, allow_incomplete=True
        ),
    )
