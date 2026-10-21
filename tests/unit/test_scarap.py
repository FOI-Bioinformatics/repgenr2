"""SCARAP core-gene builder: argv, outputs and refusals."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import ToolExecutionError, UserInputError, WorkdirError
from repgenr.coregenes import codon
from repgenr.coregenes import scarap as mod
from repgenr.coregenes.base import CoreGeneParams, registry

_LOG = logging.getLogger("test-scarap")


def _gffs(
    tmp_path: Path, n: int = 3, *, ffn: bool = True, shared_ids: bool = False, genes: int = 150
) -> list[Path]:
    """``genes`` one-codon genes (MK, ATGAAA) per genome, ids unique per genome."""
    out = []
    for i in range(1, n + 1):
        name = f"Fam_Gen_sp_GCA_{i}.1"
        d = tmp_path / "ann" / name
        d.mkdir(parents=True, exist_ok=True)
        gff = d / f"{name}.gff3"
        gff.write_text(f"##gff-version 3\n##FASTA\n>{name}\nACGT\n", encoding="utf-8")
        ids = [f"g{j}" if shared_ids else f"G{i}_{j:05d}" for j in range(genes)]
        (d / f"{name}.faa").write_text("".join(f">{g}\nMK\n" for g in ids), encoding="utf-8")
        if ffn:
            (d / f"{name}.ffn").write_text(
                "".join(f">{g}\nATGAAA\n" for g in ids), encoding="utf-8"
            )
        out.append(gff)
    return out


def _fake(calls: list[list[str]], *, orthogroups: int = 150):
    """scarap core puts the j-th gene of every genome in orthogroup F<j>."""

    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        faas = [Path(p) for p in Path(cmd[2]).read_text(encoding="utf-8").split()]
        if cmd[1] == "core":
            out = Path(cmd[3])
            out.mkdir(parents=True)
            rows = []
            for faa in faas:
                ids = [
                    line[1:].strip()
                    for line in faa.read_text(encoding="utf-8").splitlines()
                    if line.startswith(">")
                ]
                rows += [f"{g}\t{faa.stem}\tF{j:04d}\n" for j, g in enumerate(ids[:orthogroups])]
            (out / "genes.tsv").write_text("".join(rows), encoding="utf-8")
        return 0

    return run_tool


@pytest.fixture(autouse=True)
def _mafft_copies(monkeypatch) -> list[list[list[str]]]:
    """mafft copies its input: the proteins are already aligned."""
    chains: list[list[list[str]]] = []

    def run_chain(caps, steps, **k):
        argvs = [[str(c) for c in argv] for _, argv in steps]
        chains.append(argvs)
        for argv in argvs:
            src, dst = Path(argv[-1]), Path(argv[argv.index("--out") + 1])
            dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")

    monkeypatch.setattr(codon, "run_chain", run_chain)
    return chains


def test_argv_and_outputs(tmp_path: Path, monkeypatch, _mafft_copies) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    gffs = _gffs(tmp_path)
    result = mod.ScarapBuilder().call(
        gffs, tmp_path / "out", CoreGeneParams(threads=4, core_threshold=0.9), _LOG
    )
    # One SCARAP call: its concat would delete the per-orthogroup alignments.
    (core,) = calls
    assert core[:2] == ["scarap", "core"]
    faas = Path(core[2]).read_text(encoding="utf-8").split()
    # The file stem is the genome name SCARAP uses, which is the record name.
    assert [Path(p).name for p in faas] == [f"{g.stem}.faa" for g in gffs]
    assert core[core.index("-f") + 1] == "0.9" and core[core.index("-t") + 1] == "4"
    assert core[core.index("-d") + 1] == "FH"
    # SCARAP's 0.90 prefilter would cap a lower threshold; the smaller is passed.
    assert core[core.index("-p") + 1] == "0.9"
    # The adapter aligns each orthogroup with mafft and hands over codon alignments.
    assert result.alignment is None and result.core_genes == 150
    assert len(result.gene_alignments) == 150
    # Batches of at most 50 orthogroups per mafft chain.
    assert [len(chain) for chain in _mafft_copies] == [50, 50, 50]
    first = result.gene_alignments[0]
    assert first.gene == "F0000" and first.path.name == "F0000.fasta"
    records = first.path.read_text(encoding="utf-8").split(">")[1:]
    assert {r.split("\n")[0]: r.split("\n")[1] for r in records} == {g.stem: "ATGAAA" for g in gffs}
    assert [t.name for t in result.tables] == ["genes.tsv"]
    assert result.clustering_identity is None
    assert result.settings["method"] == "FH" and result.settings["seeds"] == 100
    assert result.settings["scarap_core_m"] == 0 and "max_core_genes" not in result.settings


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
    assert [c[1] for c in calls] == ["core"]


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

    def groups(threshold: float, cap: int) -> list[str]:
        return [g for g, _ in mod._core_orthogroups(genes, threshold, cap)]

    assert groups(1.0, 0) == ["A"]
    assert groups(0.75, 0) == ["A", "B"]
    assert groups(0.5, 0) == ["A", "B", "C"]
    assert groups(0.5, 2) == ["A", "B"]
    assert mod._core_orthogroups(tmp_path / "missing.tsv", 0.5, 0) == []
    # The members are the genomes holding the orthogroup in exactly one copy.
    members = dict(mod._core_orthogroups(genes, 0.5, 0))
    assert members["A"] == {f"G{i}": f"a{i}" for i in range(4)}
    assert members["B"] == {f"G{i}": f"b{i}" for i in range(1, 4)}
    assert members["C"] == {"G0": "c0", "G1": "c1"}


def test_gene_ids_shared_across_genomes_are_refused(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    gffs = _gffs(tmp_path, shared_ids=True)  # every .faa holds ">g0", ">g1", ...
    with pytest.raises(WorkdirError, match="g0 is in both"):
        mod.ScarapBuilder().call(gffs, tmp_path / "o", CoreGeneParams(), _LOG)
    assert calls == []


def test_orthofinder_methods_need_orthofinder_beside_scarap(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake(calls))
    monkeypatch.setattr(mod, "runs_on_host", lambda caps: True)
    on_path = {"orthofinder", "diamond"}
    monkeypatch.setattr(mod, "host_which", lambda b, t: f"/x/{b}" if b in on_path else None)
    # O-D searches with DIAMOND, which is there.
    mod.ScarapBuilder().call(
        _gffs(tmp_path), tmp_path / "o1", CoreGeneParams(extra={"method": "O-D"}), _LOG
    )
    assert calls[0][calls[0].index("-d") + 1] == "O-D"
    # O-B searches with BLAST, which is not.
    with pytest.raises(UserInputError, match="needs blastp, makeblastdb"):
        mod.ScarapBuilder().call(
            _gffs(tmp_path), tmp_path / "o2", CoreGeneParams(extra={"method": "O-B"}), _LOG
        )
    on_path.clear()
    with pytest.raises(UserInputError, match="needs orthofinder, diamond"):
        mod.ScarapBuilder().call(
            _gffs(tmp_path), tmp_path / "o3", CoreGeneParams(extra={"method": "O-D"}), _LOG
        )
    # In the SCARAP image there is no OrthoFinder at all.
    monkeypatch.setattr(mod, "runs_on_host", lambda caps: False)
    with pytest.raises(UserInputError, match="SCARAP image does not include"):
        mod.ScarapBuilder().call(
            _gffs(tmp_path), tmp_path / "o4", CoreGeneParams(extra={"method": "O-D"}), _LOG
        )
    assert len(calls) == 1  # the core call of the first run only


def test_the_pandas_2_failure_is_explained(tmp_path: Path, monkeypatch) -> None:
    def failing(caps, cmd, **k):
        raise ToolExecutionError(
            [str(c) for c in cmd], 1, "AttributeError: 'DataFrame' object has no attribute 'append'"
        )

    monkeypatch.setattr(mod, "run_tool", failing)
    monkeypatch.setattr(mod, "runs_on_host", lambda caps: True)
    monkeypatch.setattr(mod, "host_which", lambda b, t: f"/x/{b}")
    with pytest.raises(WorkdirError, match="pandas older than 2"):
        mod.ScarapBuilder().call(
            _gffs(tmp_path), tmp_path / "o", CoreGeneParams(extra={"method": "O-D"}), _LOG
        )
    # Any other failure, and the same failure from another method, pass through.
    with pytest.raises(ToolExecutionError):
        mod.ScarapBuilder().call(_gffs(tmp_path), tmp_path / "o2", CoreGeneParams(), _LOG)


def test_protein_alphabet_hands_over_the_protein_alignments(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake([], orthogroups=3))
    gffs = _gffs(tmp_path)
    result = mod.ScarapBuilder().call(
        gffs, tmp_path / "out", CoreGeneParams(core_threshold=0.9, alphabet="protein"), _LOG
    )
    assert len(result.gene_alignments) == 3
    first = result.gene_alignments[0]
    assert first.path == tmp_path / "out" / "alignments" / "F0000.aln"
    records = first.path.read_text(encoding="utf-8").split(">")[1:]
    assert {r.split("\n")[0]: r.split("\n")[1] for r in records} == {g.stem: "MK" for g in gffs}
