"""PPanGGOLiN core-gene builder: argv, outputs, partitions and refusals."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes import ppanggolin as mod
from repgenr.coregenes.base import CoreGeneParams, registry

_LOG = logging.getLogger("test-ppanggolin")


def _gffs(tmp_path: Path, n: int = 3) -> list[Path]:
    out = []
    for i in range(1, n + 1):
        name = f"Fam_Gen_sp_GCA_{i}.1"
        p = tmp_path / "ann" / name / f"{name}.gff3"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(f"##gff-version 3\n##FASTA\n>{name}\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


def _fake(calls: list[list[str]], *, families: int = 150, columns: bool = True):
    """Write what each of the three ppanggolin calls writes."""

    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        sub = cmd[1]
        out = Path(cmd[cmd.index("-o") + 1])
        out.mkdir(parents=True, exist_ok=True)
        if sub == "workflow":
            listed = Path(cmd[cmd.index("--anno") + 1]).read_text(encoding="utf-8")
            (out / "pangenome.h5").write_text(listed, encoding="utf-8")
        elif sub == "write_pangenome":
            (out / "gene_presence_absence.Rtab").write_text("Gene\n", encoding="utf-8")
            (out / "matrix.csv").write_text("Gene\n", encoding="utf-8")
            (out / "partitions").mkdir(exist_ok=True)
            (out / "partitions" / "persistent.txt").write_text("f1\n", encoding="utf-8")
        elif sub == "msa":
            h5 = Path(cmd[cmd.index("-p") + 1])
            names = [line.split("\t")[0] for line in h5.read_text(encoding="utf-8").splitlines()]
            part = cmd[cmd.index("--partition") + 1]
            fam_dir = out / f"msa_{part}_dna"
            fam_dir.mkdir(exist_ok=True)
            for i in range(families):
                (fam_dir / f"F{i}.aln").write_text(">x\nAC\n", encoding="utf-8")
            prefix = part
            if part == "softcore":
                prefix = f"softcore_{cmd[cmd.index('--soft_core') + 1]}"
            seq = "ACGTACGT" if columns else ""
            (out / f"{prefix}_genome_alignment.aln").write_text(
                "".join(f">{n}\n{seq}\n" for n in names), encoding="utf-8"
            )
        return 0

    return run_tool


def test_argv_outputs_and_default_softcore(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    gffs = _gffs(tmp_path)
    result = mod.PPanGGOLiNBuilder().call(
        gffs, tmp_path / "out", CoreGeneParams(threads=4, core_threshold=0.9), _LOG
    )
    workflow, tables, msa = calls
    assert workflow[:2] == ["ppanggolin", "workflow"] and "--no_flat_files" in workflow
    assert workflow[workflow.index("--identity") + 1] == "0.8"
    assert workflow[workflow.index("--coverage") + 1] == "0.8"
    assert workflow[workflow.index("-c") + 1] == "4"
    # The record name is the genome name, so the leaves match the genome set.
    listing = Path(workflow[workflow.index("--anno") + 1]).read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[0] for row in listing] == [g.stem for g in gffs]
    assert tables[1] == "write_pangenome" and {"--Rtab", "--csv", "--partitions"} <= set(tables)
    assert msa[msa.index("--partition") + 1] == "softcore"
    assert msa[msa.index("--soft_core") + 1] == "0.9"
    assert msa[msa.index("--source") + 1] == "dna" and "--phylo" in msa
    assert "--single_copy" in msa
    # Temporary files stay under the builder's directory.
    tmp = str(tmp_path / "out" / "tmp")
    assert workflow[workflow.index("--tmpdir") + 1] == tmp == msa[msa.index("--tmpdir") + 1]
    assert result.settings == {
        "partition": "softcore",
        "coverage": 0.8,
        "single_copy": True,
        "clustering_on": "protein",
    }
    assert result.alignment.name == "softcore_0.9_genome_alignment.aln"
    assert result.core_genes == 150
    assert {t.name for t in result.tables} >= {"gene_presence_absence.Rtab", "persistent.txt"}


def test_extras_reach_the_right_call(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    mod.PPanGGOLiNBuilder().call(
        _gffs(tmp_path),
        tmp_path / "out",
        CoreGeneParams(
            extra={
                "partition": "core",
                "identity": "0.5",
                "coverage": 0.6,
                "mode": 1,
                "nb_of_partitions": 3,
                "translation_table": 4,
                "dup_margin": 0.1,
            }
        ),
        _LOG,
    )
    workflow, _, msa = calls
    assert workflow[workflow.index("--identity") + 1] == "0.5"
    assert workflow[workflow.index("--coverage") + 1] == "0.6"
    assert workflow[workflow.index("-K") + 1] == "3"
    assert workflow[workflow.index("--translation_table") + 1] == "4"
    assert msa[msa.index("--partition") + 1] == "core" and "--soft_core" not in msa
    assert msa[msa.index("--dup_margin") + 1] == "0.1"

    calls.clear()
    mod.PPanGGOLiNBuilder().call(
        _gffs(tmp_path), tmp_path / "o2", CoreGeneParams(extra={"single_copy": "false"}), _LOG
    )
    assert "--single_copy" not in calls[-1]


def test_bad_values_are_refused_before_any_call(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    for extra, match in (
        ({"partition": "shell"}, "partition"),
        ({"identity": "high"}, "identity"),
        ({"coverage": "1.5"}, "coverage"),
        ({"identity": "0"}, "identity"),
        ({"dup_margin": "2"}, "dup_margin"),
    ):
        with pytest.raises(UserInputError, match=match):
            mod.PPanGGOLiNBuilder().call(
                _gffs(tmp_path), tmp_path / "bad", CoreGeneParams(extra=extra), _LOG
            )
    with pytest.raises(UserInputError, match="core-threshold"):
        mod.PPanGGOLiNBuilder().call(
            _gffs(tmp_path), tmp_path / "bad", CoreGeneParams(core_threshold=0), _LOG
        )
    assert calls == []


def test_persistent_on_a_small_set_warns(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([]))
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        mod.PPanGGOLiNBuilder().call(
            _gffs(tmp_path), tmp_path / "o", CoreGeneParams(extra={"partition": "persistent"}), _LOG
        )
    assert "not robust below 15" in caplog.text
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        mod.PPanGGOLiNBuilder().call(
            _gffs(tmp_path, 16),
            tmp_path / "o2",
            CoreGeneParams(extra={"partition": "persistent"}),
            _LOG,
        )
    assert caplog.text == ""


def test_a_small_core_warns_and_an_empty_one_refuses(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([], families=12))
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        result = mod.PPanGGOLiNBuilder().call(
            _gffs(tmp_path), tmp_path / "o", CoreGeneParams(), _LOG
        )
    assert result.core_genes == 12 and "Only 12 gene families" in caplog.text

    monkeypatch.setattr(mod, "run_tool", _fake([], families=0))
    with pytest.raises(WorkdirError, match=r"no gene family .* at --core-threshold 0.95"):
        mod.PPanGGOLiNBuilder().call(_gffs(tmp_path), tmp_path / "o2", CoreGeneParams(), _LOG)
    monkeypatch.setattr(mod, "run_tool", _fake([], columns=False))
    with pytest.raises(WorkdirError, match="identity= and coverage="):
        mod.PPanGGOLiNBuilder().call(_gffs(tmp_path), tmp_path / "o3", CoreGeneParams(), _LOG)


def test_registered_never_auto_and_ignores_panaroo_settings() -> None:
    assert "ppanggolin" in registry.names()
    caps = registry.get("ppanggolin").capabilities
    assert caps.auto_eligible is False
    assert {"clean_mode", "filtered"} <= caps.ignored_params


def test_settings_a_partition_ignores_are_named(tmp_path: Path, monkeypatch, caplog) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([]))
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        mod.PPanGGOLiNBuilder().call(
            _gffs(tmp_path),
            tmp_path / "o",
            CoreGeneParams(
                core_threshold=0.9,
                extra={"partition": "core", "single_copy": "false", "dup_margin": 0.1},
            ),
            _LOG,
        )
    assert "does not read --core-threshold" in caplog.text
    assert "dup_margin applies only with single_copy" in caplog.text
