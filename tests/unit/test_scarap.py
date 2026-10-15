"""SCARAP core-gene builder: argv, outputs and refusals."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes import scarap as mod
from repgenr.coregenes.base import CoreGeneParams, registry

_LOG = logging.getLogger("test-scarap")


def _gffs(tmp_path: Path, n: int = 3, *, ffn: bool = True, shared_ids: bool = False) -> list[Path]:
    out = []
    for i in range(1, n + 1):
        name = f"Fam_Gen_sp_GCA_{i}.1"
        d = tmp_path / "ann" / name
        d.mkdir(parents=True, exist_ok=True)
        gff = d / f"{name}.gff3"
        gff.write_text(f"##gff-version 3\n##FASTA\n>{name}\nACGT\n", encoding="utf-8")
        gene = "g1" if shared_ids else f"G{i}_00001"
        (d / f"{name}.faa").write_text(f">{gene}\nMK\n", encoding="utf-8")
        if ffn:
            (d / f"{name}.ffn").write_text(f">{gene}\nATGAAA\n", encoding="utf-8")
        out.append(gff)
    return out


def _fake(calls: list[list[str]], *, orthogroups: int = 150, supermatrix: bool = True):
    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        names = [Path(p).stem for p in Path(cmd[2]).read_text(encoding="utf-8").split()]
        if cmd[1] == "core":
            out = Path(cmd[3])
            out.mkdir(parents=True)
            rows = [f"g{i}_{n}\t{n}\tF{i:04d}\n" for i in range(orthogroups) for n in names]
            (out / "genes.tsv").write_text("".join(rows), encoding="utf-8")
        elif cmd[1] == "concat":
            out = Path(cmd[4])
            out.mkdir(parents=True)
            (out / "supermatrix_aas.fasta").write_text(">x\nMK\n", encoding="utf-8")
            if supermatrix:
                (out / mod.SUPERMATRIX).write_text(
                    "".join(f">{n}\nATGAAA\n" for n in names), encoding="utf-8"
                )
        return 0

    return run_tool


def test_argv_and_outputs(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    gffs = _gffs(tmp_path)
    result = mod.ScarapBuilder().call(
        gffs, tmp_path / "out", CoreGeneParams(threads=4, core_threshold=0.9), _LOG
    )
    core, concat = calls
    assert core[:2] == ["scarap", "core"]
    faas = Path(core[2]).read_text(encoding="utf-8").split()
    # The file stem is the genome name SCARAP uses, which is the record name.
    assert [Path(p).name for p in faas] == [f"{g.stem}.faa" for g in gffs]
    assert core[core.index("-f") + 1] == "0.9" and core[core.index("-t") + 1] == "4"
    assert core[core.index("-d") + 1] == "FH"
    # SCARAP's 0.90 prefilter would cap a lower threshold; the smaller is passed.
    assert core[core.index("-p") + 1] == "0.9"
    assert concat[:2] == ["scarap", "concat"] and concat[3].endswith("genes.tsv")
    # The core is judged again over all genomes, not only the seeds.
    assert concat[concat.index("-f") + 1] == "0.9" and concat[2] == core[2]
    ffns = Path(concat[concat.index("-n") + 1]).read_text(encoding="utf-8").split()
    assert [Path(p).suffix for p in ffns] == [".ffn"] * 3
    assert result.alignment.name == mod.SUPERMATRIX and result.core_genes == 150
    assert [t.name for t in result.tables] == ["genes.tsv"]
    assert result.clustering_identity is None
    assert result.settings["method"] == "FH" and result.settings["seeds"] == 100


def test_extras_reach_the_core_call(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    mod.ScarapBuilder().call(
        _gffs(tmp_path),
        tmp_path / "out",
        CoreGeneParams(
            extra={"method": "H", "core_prefilter": "0.8", "max_core_genes": 50, "seeds": 20}
        ),
        _LOG,
    )
    core = calls[0]
    assert core[core.index("-d") + 1] == "H"
    assert core[core.index("-p") + 1] == "0.8"
    assert core[core.index("-m") + 1] == "50" and core[core.index("-e") + 1] == "20"


def test_bad_values_and_missing_files_are_refused_before_any_call(
    tmp_path: Path, monkeypatch
) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    for extra, match in (
        ({"method": "BLAST"}, "method"),
        ({"method": "O-D"}, "method"),
        ({"core_prefilter": "1.5"}, "core_prefilter"),
        ({"core_prefilter": "x"}, "core_prefilter"),
        ({"seeds": 0}, "seeds"),
        ({"max_core_genes": -1}, "max_core_genes"),
    ):
        with pytest.raises(UserInputError, match=match):
            mod.ScarapBuilder().call(
                _gffs(tmp_path), tmp_path / "bad", CoreGeneParams(extra=extra), _LOG
            )
    with pytest.raises(UserInputError, match="core-threshold"):
        mod.ScarapBuilder().call(
            _gffs(tmp_path), tmp_path / "bad", CoreGeneParams(core_threshold=1.2), _LOG
        )
    with pytest.raises(WorkdirError, match="annotate --force"):
        mod.ScarapBuilder().call(
            _gffs(tmp_path / "noffn", ffn=False), tmp_path / "bad", CoreGeneParams(), _LOG
        )
    assert calls == []


def test_small_empty_and_failed_cores(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([], orthogroups=12))
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        result = mod.ScarapBuilder().call(_gffs(tmp_path), tmp_path / "o1", CoreGeneParams(), _LOG)
    assert result.core_genes == 12 and "Only 12 core genes" in caplog.text

    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls, orthogroups=0))
    with pytest.raises(WorkdirError, match="no core gene"):
        mod.ScarapBuilder().call(_gffs(tmp_path), tmp_path / "o2", CoreGeneParams(), _LOG)
    assert [c[1] for c in calls] == ["core"]  # no alignment pass over nothing

    monkeypatch.setattr(mod, "run_tool", _fake([], supermatrix=False))
    with pytest.raises(WorkdirError, match="wrote no supermatrix_nucs.fasta"):
        mod.ScarapBuilder().call(_gffs(tmp_path), tmp_path / "o3", CoreGeneParams(), _LOG)


def test_registered_never_auto_and_ignores_panaroo_settings() -> None:
    assert "scarap" in registry.names()
    caps = registry.get("scarap").capabilities
    assert caps.auto_eligible is False
    assert {"clean_mode", "filtered"} <= caps.ignored_params


def test_a_lower_threshold_lowers_the_prefilter(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    result = mod.ScarapBuilder().call(
        _gffs(tmp_path), tmp_path / "o", CoreGeneParams(core_threshold=0.8), _LOG
    )
    assert calls[0][calls[0].index("-p") + 1] == "0.8"
    assert result.settings["core_prefilter"] == 0.8


def test_the_count_follows_scarap_single_copy_rule(tmp_path: Path) -> None:
    genes = tmp_path / "genes.tsv"
    rows = [
        # A in all four genomes once: kept.
        *[f"a{i}\tG{i}\tA\n" for i in range(4)],
        # B in all four but twice in G0: single-copy in 3 of 4.
        *[f"b{i}\tG{i}\tB\n" for i in range(4)],
        "b9\tG0\tB\n",
        # C in two genomes: 2 of 4.
        "c0\tG0\tC\n",
        "c1\tG1\tC\n",
    ]
    genes.write_text("".join(rows), encoding="utf-8")
    assert mod._core_orthogroups(genes, 1.0, 0) == ["A"]
    assert mod._core_orthogroups(genes, 0.75, 0) == ["A", "B"]
    assert mod._core_orthogroups(genes, 0.5, 0) == ["A", "B", "C"]
    assert mod._core_orthogroups(genes, 0.5, 2) == ["A", "B"]
    assert mod._core_orthogroups(tmp_path / "missing.tsv", 0.5, 0) == []


def test_gene_ids_shared_across_genomes_are_refused(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    gffs = _gffs(tmp_path, shared_ids=True)  # every .faa holds ">g1"
    with pytest.raises(WorkdirError, match="g1 is in both"):
        mod.ScarapBuilder().call(gffs, tmp_path / "o", CoreGeneParams(), _LOG)
    assert calls == []
