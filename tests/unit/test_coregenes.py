"""Core-gene builders: the Panaroo adapter and the stateless core-gene pass."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from repgenr.core import annotations as cache
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes import panaroo as mod
from repgenr.coregenes.base import CoreGeneParams, CoreGeneResult, registry
from repgenr.stages.coregene import annotated_set, coregene_core

_LOG = logging.getLogger("test-coregene")


def _gffs(tmp_path: Path, names=("Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_3")):
    out = []
    for n in names:
        p = tmp_path / f"{n}.gff3"
        p.write_text(f"##gff-version 3\n##FASTA\n>{n}\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


def _header(genes: int) -> str:
    """Panaroo's alignment header: one `FT   feature` block per gene."""
    out = "ID   Genome standard; DNA; PRO; 1234 BP.\nFH   Key             Location/Qualifiers\n"
    for i in range(genes):
        out += f"FT   feature         {i * 100 + 1}..{(i + 1) * 100}\n"
        out += f'FT                   /label="gene_{i}"\n'
    return out


def _fake_panaroo(
    calls: list[list[str]], *, filtered: bool = True, unfiltered: bool = True, header: bool = True
):
    def run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        out = Path(cmd[cmd.index("-o") + 1])
        listed = Path(cmd[cmd.index("-i") + 1]).read_text(encoding="utf-8").split()
        names = [Path(p).stem for p in listed]
        aln = "".join(f">{n}\nACGTACGT\n" for n in names)
        if unfiltered:
            (out / "core_gene_alignment.aln").write_text(aln, encoding="utf-8")
        if filtered:
            (out / "core_gene_alignment_filtered.aln").write_text(
                aln.replace("ACGTACGT", "ACGT"), encoding="utf-8"
            )
        (out / "gene_presence_absence.csv").write_text("Gene\n", encoding="utf-8")
        (out / "gene_presence_absence.Rtab").write_text("Gene\n", encoding="utf-8")
        summary = "Core genes\t(99% <= strains <= 100%)\t1234\n"
        summary += "Total genes\t(0% <= strains <= 100%)\t2000\n"
        (out / "summary_statistics.txt").write_text(summary, encoding="utf-8")
        if header:
            # Panaroo's own 99-100% bin above; the alignment holds fewer genes.
            (out / "core_alignment_header.embl").write_text(_header(20), encoding="utf-8")
            if filtered:
                (out / "core_alignment_filtered_header.embl").write_text(
                    _header(17), encoding="utf-8"
                )
        return 0

    return run_tool


def test_panaroo_argv_outputs_and_summary(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo(calls))
    gffs = _gffs(tmp_path)
    result = mod.PanarooBuilder().call(
        gffs,
        tmp_path / "out",
        CoreGeneParams(
            threads=4,
            core_threshold=0.9,
            clean_mode="moderate",
            extra={"merge_paralogs": "true", "family_threshold": 0.7},
        ),
        _LOG,
    )
    (cmd,) = calls
    assert (
        cmd[:2] == ["panaroo", "-i"] and Path(cmd[2]).read_text(encoding="utf-8").count("\n") == 3
    )
    assert cmd[cmd.index("--clean-mode") + 1] == "moderate"
    assert cmd[cmd.index("-a") + 1] == "core" and cmd[cmd.index("--aligner") + 1] == "mafft"
    assert cmd[cmd.index("--core_threshold") + 1] == "0.9" and cmd[cmd.index("-t") + 1] == "4"
    assert "--merge_paralogs" in cmd and cmd[cmd.index("--family_threshold") + 1] == "0.7"
    assert result.alignment.name == "core_gene_alignment_filtered.aln"
    # The genes of the filtered alignment, not the 1234 of Panaroo's 99% bin.
    assert result.core_genes == 17 and len(result.tables) == 5


def test_panaroo_unfiltered_choice_and_errors(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo(calls))
    result = mod.PanarooBuilder().call(
        _gffs(tmp_path), tmp_path / "o1", CoreGeneParams(filtered=False), _LOG
    )
    assert result.alignment.name == "core_gene_alignment.aln" and result.core_genes == 20
    with pytest.raises(UserInputError, match="clean-mode"):
        mod.PanarooBuilder().call(
            _gffs(tmp_path), tmp_path / "o2", CoreGeneParams(clean_mode="x"), _LOG
        )
    with pytest.raises(UserInputError, match="core-threshold"):
        mod.PanarooBuilder().call(
            _gffs(tmp_path), tmp_path / "o3", CoreGeneParams(core_threshold=0), _LOG
        )
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo([], filtered=False))
    with pytest.raises(WorkdirError, match="filtered=false"):
        mod.PanarooBuilder().call(_gffs(tmp_path), tmp_path / "o4", CoreGeneParams(), _LOG)
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo([], filtered=False, unfiltered=False))
    with pytest.raises(WorkdirError, match="wrote no core_gene_alignment"):
        mod.PanarooBuilder().call(
            _gffs(tmp_path), tmp_path / "o5", CoreGeneParams(filtered=False), _LOG
        )


def test_panaroo_is_registered_and_never_auto() -> None:
    assert "panaroo" in registry.names()
    assert registry.get("panaroo").capabilities.auto_eligible is False


def _index(tmp_path: Path, names, annotator="bakta") -> Path:
    d = tmp_path / "annotations"
    records = {}
    for n in names:
        gff = cache.gff_path(d, n)
        gff.parent.mkdir(parents=True, exist_ok=True)
        gff.write_text(f"##gff-version 3\n##FASTA\n>{n}\nACGT\n", encoding="utf-8")
        records[n] = cache.AnnotationRecord(
            n, "genomes", "", "d", annotator, "1", "6", "", cache.gff_digest(gff), "t"
        )
    cache.write_index(d, records)
    return d


def test_annotated_set_resolves_refuses_missing_and_mixed(tmp_path: Path) -> None:
    d = _index(tmp_path, ["a", "b"])
    got = annotated_set(d, ["b", "a"])
    assert [g.name for g in got.gffs] == ["b.gff3", "a.gff3"] and got.annotator == "bakta"
    with pytest.raises(UserInputError, match="no annotation.*repgenr annotate"):
        annotated_set(d, ["a", "c"])
    records = cache.read_index(d)
    records["b"] = cache.AnnotationRecord("b", "genomes", "", "d", "prokka", "1", "", "", "g", "t")
    cache.write_index(d, records)
    with pytest.raises(UserInputError, match="mix annotators"):
        annotated_set(d, ["a", "b"])
    cache.gff_path(d, "a").unlink()
    with pytest.raises(WorkdirError, match="missing; rerun annotate"):
        annotated_set(d, ["a"])


class _FakeBuilder:
    capabilities = mod.PanarooBuilder.capabilities
    accepted_annotators = ()

    def preflight(self):
        return {"fakecg": "1.0"}

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        aln = out_dir / "aln.fasta"
        aln.write_text("".join(f">{g.stem}\nACGTAC\n" for g in gffs), encoding="utf-8")
        table = out_dir / "summary_statistics.txt"
        table.write_text("Core genes\t(99% <= strains <= 100%)\t7\n", encoding="utf-8")
        return CoreGeneResult(alignment=aln, tables=[table], core_genes=7)


def test_coregene_core_publishes_alignment_tables_and_stats(
    tmp_path: Path, register_tool, caplog
) -> None:
    from repgenr.core.contracts import CORE_GENE_FASTA, CORE_GENE_STATS, PANGENOME_DIR

    register_tool(registry, "fakecg", _FakeBuilder)
    gffs = _gffs(tmp_path)
    with caplog.at_level(logging.INFO, logger="test-coregene"):
        result, versions = coregene_core(
            gffs, tmp_path / "msa", tmp_path / "scratch", CoreGeneParams(), "fakecg", _LOG
        )
    assert result.alignment == tmp_path / "msa" / CORE_GENE_FASTA and result.alignment.is_file()
    assert (tmp_path / "msa" / PANGENOME_DIR / "summary_statistics.txt").is_file()
    stats = json.loads((tmp_path / "msa" / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["core_genes"] == 7 and stats["alignment_length"] == 6 and stats["genomes"] == 3
    assert versions == {"fakecg": "1.0"} and not (tmp_path / "scratch").exists()
    assert any("no ascertainment-bias correction" in r.getMessage() for r in caplog.records)
    with pytest.raises(UserInputError, match="Unknown core-gene builder"):
        coregene_core(gffs, tmp_path / "m2", tmp_path / "s2", CoreGeneParams(), "nope", _LOG)


def test_the_gene_count_falls_back_and_then_gives_up(tmp_path: Path, monkeypatch, caplog) -> None:
    """Without the filtered header the unfiltered one is read; without either
    the count is unknown, never Panaroo's 99% bin."""
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo([], header=True))
    out = tmp_path / "out"
    result = mod.PanarooBuilder().call(_gffs(tmp_path), out, CoreGeneParams(), _LOG)
    assert result.core_genes == 17
    (out / "core_alignment_filtered_header.embl").unlink()
    with caplog.at_level(logging.INFO, logger="test-coregene"):
        again = mod._core_genes(out / "core_alignment_filtered_header.embl", out, _LOG)
    assert again == 20
    assert any("before the entropy filter" in r.getMessage() for r in caplog.records)
    (out / "core_alignment_header.embl").unlink()
    assert mod._core_genes(out / "core_alignment_header.embl", out, _LOG) is None
