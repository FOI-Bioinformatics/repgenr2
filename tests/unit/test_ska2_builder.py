"""ska2 tree builder: build and distance calls, matrix parsing, warnings."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.treebuilders import ska2 as mod
from repgenr.treebuilders.base import TreeParams

_LOG = logging.getLogger("test-ska2-tree")
_HEADER = "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"


def _genomes(tmp_path: Path, names=("Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_3")):
    out = []
    for n in names:
        p = tmp_path / f"{n}.fasta"
        p.write_text(f">{n}\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


def _fake(calls: list[list[str]], rows: str):
    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        if cmd[1] == "build":
            Path(cmd[cmd.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif cmd[1] == "distance":
            Path(cmd[cmd.index("-o") + 1]).write_text(_HEADER + rows, encoding="utf-8")
        return 0

    return run_tool


_ROWS = (
    "Fam_Gen_sp_GCA_2\tFam_Gen_sp_GCA_1\t2.00\t0.27907\t310\t120\n"
    "Fam_Gen_sp_GCA_2\tFam_Gen_sp_GCA_3\t5.00\t0.57692\t220\t300\n"
    "Fam_Gen_sp_GCA_1\tFam_Gen_sp_GCA_3\t3.00\t0.39130\t280\t180\n"
)


def test_build_then_distance_then_join(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls, _ROWS))
    seen: dict = {}

    def fake_nj(labels, dist, out_dir, logger, *, backend, threads=1):
        seen.update(labels=labels, dist=dist, backend=backend, threads=threads)
        return "(a,b,c);"

    monkeypatch.setattr(mod, "nj_tree", fake_nj)
    tree = mod.Ska2Builder().build(
        _genomes(tmp_path), tmp_path / "out", TreeParams(threads=4, extra={"ksize": 21}), _LOG
    )
    build, distance = calls
    assert build[:2] == ["ska", "build"] and build[build.index("-k") + 1] == "21"
    assert distance[:2] == ["ska", "distance"]
    assert distance[distance.index("--threads") + 1] == "4"
    assert distance[distance.index("-m") + 1] == "0" and distance[-1].endswith(".skf")
    assert tree.read_text(encoding="utf-8") == "(a,b,c);\n"
    assert seen["labels"] == ["Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_3"]
    assert seen["dist"][0][1] == 2.0 and seen["dist"][1][2] == 5.0 and seen["dist"][0][2] == 3.0
    assert seen["backend"] in ("python", "rapidnj") and seen["threads"] == 4


def test_mismatch_metric_and_bad_metric(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls, _ROWS))
    seen: dict = {}
    monkeypatch.setattr(
        mod, "nj_tree", lambda labels, dist, *a, **k: seen.update(dist=dist) or "(a);"
    )
    mod.Ska2Builder().build(
        _genomes(tmp_path), tmp_path / "out", TreeParams(extra={"metric": "mismatch"}), _LOG
    )
    assert seen["dist"][0][1] == pytest.approx(0.27907)
    with pytest.raises(UserInputError, match="metric"):
        mod.Ska2Builder().build(
            _genomes(tmp_path), tmp_path / "out2", TreeParams(extra={"metric": "jaccard"}), _LOG
        )


def test_a_missing_pair_is_refused(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls, _ROWS.splitlines(keepends=True)[0]))
    with pytest.raises(WorkdirError, match="2 of 3 pairs have no row"):
        mod.Ska2Builder().build(_genomes(tmp_path), tmp_path / "out", TreeParams(), _LOG)


def test_a_pair_without_shared_kmers_is_refused(tmp_path: Path, monkeypatch) -> None:
    # ska 0.5.1 writes an unrelated pair as "0.00 1.00000 0 N": zero matches.
    rows = _ROWS.replace(
        "Fam_Gen_sp_GCA_2\tFam_Gen_sp_GCA_3\t5.00\t0.57692\t220\t300",
        "Fam_Gen_sp_GCA_2\tFam_Gen_sp_GCA_3\t0.00\t1.00000\t0\t199940",
    )
    monkeypatch.setattr(mod, "run_tool", _fake([], rows))
    with pytest.raises(WorkdirError, match="share no split k-mer.*GCA_2 and Fam_Gen_sp_GCA_3"):
        mod.Ska2Builder().build(_genomes(tmp_path), tmp_path / "out", TreeParams(), _LOG)


def test_the_outgroup_is_left_out_of_the_species_count(tmp_path: Path, caplog) -> None:
    genomes = _genomes(tmp_path, ("Fam_Gen_spA_GCA_1", "Fam_Gen_spA_GCA_2", "Out_Gen_spB_GCA_9"))
    with caplog.at_level(logging.WARNING, logger="test-ska2-tree"):
        mod.warn_ska2_divergence(genomes, _LOG, exclude=("Out_Gen_spB_GCA_9",))
    assert not caplog.records
    with caplog.at_level(logging.WARNING, logger="test-ska2-tree"):
        mod.warn_ska2_divergence(genomes, _LOG)
    assert len(caplog.records) == 1


def test_ska2_builder_is_never_chosen_by_auto() -> None:
    from repgenr.core.plugins import auto_select
    from repgenr.treebuilders.base import registry

    assert mod.Ska2Builder.capabilities.auto_eligible is False
    assert auto_select(registry, 3000, allowed={"ska2", "mashtree"}) == "mashtree"
    assert auto_select(registry, 3000, allowed={"ska2"}) is None


def test_read_distances_rejects_foreign_samples_and_headers(tmp_path: Path) -> None:
    table = tmp_path / "d.tsv"
    table.write_text(_HEADER + "x\ty\t1\t0.1\t1\t1\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="not an input genome"):
        mod.read_distances(table, ["a", "b"], "snps")
    table.write_text("foo\tbar\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="not a ska distance table"):
        mod.read_distances(table, ["a", "b"], "snps")


def test_warns_across_species(tmp_path: Path, monkeypatch, caplog) -> None:
    calls: list[list[str]] = []
    rows = (
        "Fam_Gen_spA_GCA_1\tFam_Gen_spB_GCA_2\t9.00\t0.5\t1\t1\n"
        "Fam_Gen_spA_GCA_1\tFam_Gen_spB_GCA_3\t9.00\t0.5\t1\t1\n"
        "Fam_Gen_spB_GCA_2\tFam_Gen_spB_GCA_3\t1.00\t0.1\t1\t1\n"
    )
    monkeypatch.setattr(mod, "run_tool", _fake(calls, rows))
    monkeypatch.setattr(mod, "nj_tree", lambda *a, **k: "(a,b,c);")
    genomes = _genomes(tmp_path, ("Fam_Gen_spA_GCA_1", "Fam_Gen_spB_GCA_2", "Fam_Gen_spB_GCA_3"))
    with caplog.at_level(logging.WARNING, logger="test-ska2-tree"):
        mod.Ska2Builder().build(genomes, tmp_path / "out", TreeParams(), _LOG)
    assert any("spanning 2 species" in r.getMessage() for r in caplog.records)
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="test-ska2-tree"):
        mod.warn_ska2_divergence(_genomes(tmp_path), _LOG)
    assert not caplog.records


def test_build_settles_the_joining_step_before_ska_runs(tmp_path: Path, monkeypatch) -> None:
    from repgenr.tree import nj

    monkeypatch.setattr(nj, "tool_available", lambda caps: False)
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: pytest.fail("ska ran"))
    genomes = [tmp_path / f"Fam_Gen_sp_GCA_{i}.fasta" for i in range(nj.PYTHON_NJ_MAX + 1)]
    with pytest.raises(WorkdirError, match="Neighbour joining"):
        mod.Ska2Builder().build(genomes, tmp_path / "out", TreeParams(), _LOG)


def test_gzipped_inputs_keep_one_name_in_listing_and_table(tmp_path: Path, monkeypatch) -> None:
    """The listing names samples by record_name (gz suffix stripped) and the
    table is looked up by the same names, so ska's echoed names match."""
    import gzip

    genomes = []
    for n in ("Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_3"):
        p = tmp_path / f"{n}.fasta.gz"
        with gzip.open(p, "wt") as fo:
            fo.write(f">{n}\nACGT\n")
        genomes.append(p)
    calls: list[list[str]] = []
    seen: dict = {}

    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        if cmd[1] == "build":
            listing = Path(cmd[cmd.index("-f") + 1]).read_text(encoding="utf-8").splitlines()
            seen["names"] = [ln.split("\t")[0] for ln in listing]
            Path(cmd[cmd.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif cmd[1] == "distance":
            names = seen["names"]
            rows = [_HEADER.rstrip("\n")]
            for i, a in enumerate(names):
                for b in names[i + 1 :]:
                    rows.append(f"{a}\t{b}\t1.00\t0.1\t100\t1")
            Path(cmd[cmd.index("-o") + 1]).write_text("\n".join(rows) + "\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(mod, "run_tool", run_tool)
    monkeypatch.setattr(
        mod, "nj_tree", lambda labels, *a, **k: seen.update(labels=labels) or "(a);"
    )
    mod.Ska2Builder().build(genomes, tmp_path / "out", TreeParams(), _LOG)
    assert seen["names"] == ["Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_3"]
    assert seen["labels"] == seen["names"]
