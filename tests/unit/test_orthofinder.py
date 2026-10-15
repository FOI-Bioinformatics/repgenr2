"""OrthoFinder core-gene builder: argv, selection, back-translation, refusals."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.binaries import _parse_version
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes import orthofinder as mod
from repgenr.coregenes.base import CoreGeneParams, registry

_LOG = logging.getLogger("test-orthofinder")

# Two genes per genome: K-L-M (ATG AAA CTG ATG) and a second, A-C.
_PROT = {"a": "MKLM", "b": "AC"}
_NUC = {"a": "ATGAAACTGATGTAA", "b": "GCTTGT"}


def _gffs(tmp_path: Path, n: int = 3) -> list[Path]:
    out = []
    for i in range(1, n + 1):
        name = f"Fam_Gen_sp_GCA_{i}.1"
        d = tmp_path / "ann" / name
        d.mkdir(parents=True, exist_ok=True)
        gff = d / f"{name}.gff3"
        gff.write_text("##gff-version 3\n", encoding="utf-8")
        (d / f"{name}.faa").write_text(
            "".join(f">{k}{i}\n{v}\n" for k, v in _PROT.items()), encoding="utf-8"
        )
        (d / f"{name}.ffn").write_text(
            "".join(f">{k}{i}\n{v}\n" for k, v in _NUC.items()), encoding="utf-8"
        )
        out.append(gff)
    return out


def _fake(calls: list[list[str]], rows: list[list[str]] | None = None):
    """OrthoFinder writes Orthogroups.tsv; mafft copies its input (one gap added)."""

    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        faa_dir = Path(cmd[cmd.index("-f") + 1])
        names = sorted(p.stem for p in faa_dir.glob("*.faa"))
        res = Path(cmd[cmd.index("-o") + 1]) / "Results_Oct15" / "Orthogroups"
        res.mkdir(parents=True)
        body = rows
        if body is None:
            body = [
                ["OG0000000", *[f"a{i + 1}" for i in range(len(names))]],
                # b is duplicated in the first genome: single-copy in 2 of 3.
                ["OG0000001", f"b1, a{9}", *[f"b{i + 1}" for i in range(1, len(names))]],
            ]
        text = "\t".join(["Orthogroup", *names]) + "\n"
        text += "".join("\t".join(r) + "\n" for r in body)
        (res / "Orthogroups.tsv").write_text(text, encoding="utf-8")
        (res / "Orthogroups.GeneCount.tsv").write_text("Orthogroup\n", encoding="utf-8")
        return 0

    return run_tool


def _fake_chain(chains: list[list[list[str]]]):
    def run_chain(caps, steps, **k):
        argvs = [[str(c) for c in argv] for _, argv in steps]
        chains.append(argvs)
        for argv in argvs:
            src, dst = Path(argv[-1]), Path(argv[argv.index("--out") + 1])
            dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")

    return run_chain


def test_argv_selection_and_back_translation(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    monkeypatch.setattr(mod, "run_chain", _fake_chain(chains))
    gffs = _gffs(tmp_path)
    result = mod.OrthoFinderBuilder().call(
        gffs, tmp_path / "out", CoreGeneParams(threads=4, core_threshold=0.6), _LOG
    )
    (cmd,) = calls
    assert cmd[:2] == ["orthofinder", "-f"] and {"-og", "-X"} <= set(cmd)
    assert cmd[cmd.index("-S") + 1] == "diamond" and cmd[cmd.index("-t") + 1] == "4"
    # One proteome per genome, named by the record name.
    faa_dir = Path(cmd[cmd.index("-f") + 1])
    assert sorted(p.stem for p in faa_dir.glob("*.faa")) == sorted(g.stem for g in gffs)
    # OG0000000 in all three, OG0000001 single-copy in two of three: both kept at 0.6.
    assert result.core_genes == 2
    aligned = [argv for chain in chains for argv in chain]
    assert len(aligned) == 2 and all(a[0] == "mafft" and "--auto" in a for a in aligned)
    records = result.alignment.read_text(encoding="utf-8").split(">")[1:]
    seqs = {r.split("\n")[0]: r.split("\n")[1] for r in records}
    assert list(seqs) == [g.stem for g in gffs]
    # Codons, the stop codon dropped; the duplicated genome is a gap for OG0000001.
    assert seqs[gffs[1].stem] == "ATGAAACTGATG" + "GCTTGT"
    assert seqs[gffs[0].stem] == "ATGAAACTGATG" + "------"
    assert result.settings["search"] == "diamond"
    assert {t.name for t in result.tables} == {"Orthogroups.GeneCount.tsv", mod.SELECTED}


def test_threshold_and_extras(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    monkeypatch.setattr(mod, "run_chain", _fake_chain([]))
    result = mod.OrthoFinderBuilder().call(
        _gffs(tmp_path),
        tmp_path / "out",
        CoreGeneParams(core_threshold=1.0, extra={"search": "mmseqs", "inflation": "1.5"}),
        _LOG,
    )
    assert result.core_genes == 1  # only OG0000000 is single-copy everywhere
    cmd = calls[0]
    assert cmd[cmd.index("-S") + 1] == "mmseqs" and cmd[cmd.index("-I") + 1] == "1.5"


def test_refusals_before_any_call(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    for extra, match in (
        ({"search": "blastn"}, "search"),
        ({"inflation": "x"}, "inflation"),
        ({"inflation": "1.0"}, "above 1"),
    ):
        with pytest.raises(UserInputError, match=match):
            mod.OrthoFinderBuilder().call(
                _gffs(tmp_path), tmp_path / "bad", CoreGeneParams(extra=extra), _LOG
            )
    gffs = _gffs(tmp_path / "noffn")
    gffs[0].with_suffix(".ffn").unlink()
    with pytest.raises(WorkdirError, match="annotate --force"):
        mod.OrthoFinderBuilder().call(gffs, tmp_path / "bad", CoreGeneParams(), _LOG)
    assert calls == []


def test_no_core_and_mismatched_genes(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([], rows=[["OG0000000", "a1", "", ""]]))
    monkeypatch.setattr(mod, "run_chain", _fake_chain([]))
    with pytest.raises(WorkdirError, match="no orthogroup"):
        mod.OrthoFinderBuilder().call(_gffs(tmp_path), tmp_path / "o1", CoreGeneParams(), _LOG)

    # A gene whose nucleotide length does not fit its protein becomes a gap.
    gffs = _gffs(tmp_path / "mm")
    ffn = gffs[2].with_suffix(".ffn")
    ffn.write_text(ffn.read_text(encoding="utf-8").replace("ATGAAACTGATGTAA", "ATGAA"), "utf-8")
    monkeypatch.setattr(mod, "run_tool", _fake([]))
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        result = mod.OrthoFinderBuilder().call(gffs, tmp_path / "o2", CoreGeneParams(), _LOG)
    assert "1 gene(s) were left as gaps" in caplog.text
    tail = result.alignment.read_text(encoding="utf-8").split(">")[3]
    assert tail.split("\n")[1].startswith("------------")


def test_a_single_sequence_group_is_not_sent_to_mafft(tmp_path: Path, monkeypatch) -> None:
    chains: list[list[list[str]]] = []
    monkeypatch.setattr(mod, "run_tool", _fake([], rows=[["OG0000000", "a1", "", ""]]))
    monkeypatch.setattr(mod, "run_chain", _fake_chain(chains))
    result = mod.OrthoFinderBuilder().call(
        _gffs(tmp_path), tmp_path / "o", CoreGeneParams(core_threshold=0.3), _LOG
    )
    assert result.core_genes == 1 and chains == []


def test_coloured_version_output_parses() -> None:
    coloured = (
        "\x1b[38;5;136;49mOrthoFinder\x1b[0m\x1b[39;49m:v\x1b[0m"
        "\x1b[38;5;38;49m3.\x1b[0m\x1b[1;38;5;38;49m1.5\x1b[0m\n"
    )
    assert _parse_version(coloured) == (3, 1, 5)


def test_registered_and_never_auto() -> None:
    assert "orthofinder" in registry.names()
    caps = registry.get("orthofinder").capabilities
    assert caps.auto_eligible is False and caps.recommended_max_genomes == 200


def test_selenocysteine_is_aligned_as_x(tmp_path: Path, monkeypatch) -> None:
    """mafft refuses U and O; bakta writes U into a recoded protein."""
    sent: list[str] = []

    def chain(caps, steps, **k):
        for _, argv in steps:
            argv = [str(c) for c in argv]
            text = Path(argv[-1]).read_text(encoding="utf-8")
            sent.append(text)
            Path(argv[argv.index("--out") + 1]).write_text(text, encoding="utf-8")

    gffs = _gffs(tmp_path)
    faa = gffs[0].with_suffix(".faa")
    faa.write_text(faa.read_text(encoding="utf-8").replace("MKLM", "MKUM"), "utf-8")
    monkeypatch.setattr(mod, "run_tool", _fake([]))
    monkeypatch.setattr(mod, "run_chain", chain)
    result = mod.OrthoFinderBuilder().call(gffs, tmp_path / "o", CoreGeneParams(), _LOG)
    assert any("MKXM" in t for t in sent) and not any("U" in t.split("\n", 1)[1] for t in sent)
    # The codons still map: X keeps the protein length.
    first = result.alignment.read_text(encoding="utf-8").split(">")[1].split("\n")[1]
    assert first.startswith("ATGAAACTGATG")


def test_failures_inside_the_run(tmp_path: Path, monkeypatch) -> None:
    # A gene OrthoFinder names that the protein file lacks.
    monkeypatch.setattr(mod, "run_tool", _fake([], rows=[["OG0000000", "zz1", "a2", "a3"]]))
    monkeypatch.setattr(mod, "run_chain", _fake_chain([]))
    with pytest.raises(WorkdirError, match="'zz1'"):
        mod.OrthoFinderBuilder().call(_gffs(tmp_path), tmp_path / "o1", CoreGeneParams(), _LOG)

    # No Orthogroups.tsv at all.
    monkeypatch.setattr(mod, "run_tool", lambda caps, cmd, **k: 0)
    with pytest.raises(WorkdirError, match="no Orthogroups.tsv"):
        mod.OrthoFinderBuilder().call(_gffs(tmp_path), tmp_path / "o2", CoreGeneParams(), _LOG)

    # A failing mafft batch fails the call.
    def broken(caps, steps, **k):
        raise WorkdirError("mafft failed")

    monkeypatch.setattr(mod, "run_tool", _fake([]))
    monkeypatch.setattr(mod, "run_chain", broken)
    with pytest.raises(WorkdirError, match="mafft failed"):
        mod.OrthoFinderBuilder().call(
            _gffs(tmp_path), tmp_path / "o3", CoreGeneParams(threads=4), _LOG
        )


def test_a_path_with_a_space_is_refused(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    with pytest.raises(UserInputError, match="space"):
        mod.OrthoFinderBuilder().call(
            _gffs(tmp_path), tmp_path / "with space", CoreGeneParams(), _LOG
        )
    assert calls == []


def test_the_default_inflation_is_recorded(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([]))
    monkeypatch.setattr(mod, "run_chain", _fake_chain([]))
    result = mod.OrthoFinderBuilder().call(_gffs(tmp_path), tmp_path / "o", CoreGeneParams(), _LOG)
    assert result.settings["inflation"] == mod.DEFAULT_INFLATION
