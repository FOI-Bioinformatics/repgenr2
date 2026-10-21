"""Core-gene builders: the Panaroo adapter and the stateless core-gene pass."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from repgenr.core import annotations as cache
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.coregenes import panaroo as mod
from repgenr.coregenes.base import CoreGeneParams, CoreGeneResult, GeneAlignment, registry
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
    calls: list[list[str]],
    *,
    filtered: bool = True,
    unfiltered: bool = True,
    header: bool = True,
    genes: int = 20,
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
            (out / "core_alignment_header.embl").write_text(_header(genes), encoding="utf-8")
            if filtered:
                (out / "core_alignment_filtered_header.embl").write_text(
                    _header(max(genes - 3, 0)), encoding="utf-8"
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


def test_panaroo_threshold_reaches_argv_and_is_checked(tmp_path: Path, monkeypatch) -> None:
    """``threshold=`` is Panaroo's -c, the identity of its first clustering pass."""
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo(calls))
    mod.PanarooBuilder().call(
        _gffs(tmp_path), tmp_path / "o1", CoreGeneParams(extra={"threshold": "0.90"}), _LOG
    )
    (cmd,) = calls
    assert cmd[cmd.index("-c") + 1] == "0.9"
    assert "threshold" in mod.PanarooBuilder.capabilities.accepted_extras

    calls.clear()
    mod.PanarooBuilder().call(
        _gffs(tmp_path), tmp_path / "o2", CoreGeneParams(extra={"threshold": 1.0}), _LOG
    )
    (cmd,) = calls
    assert cmd[cmd.index("-c") + 1] == "1"  # the upper bound is allowed

    calls.clear()
    mod.PanarooBuilder().call(_gffs(tmp_path), tmp_path / "o3", CoreGeneParams(), _LOG)
    (cmd,) = calls
    # Passed even unset, so the command and the record carry the governing value.
    assert cmd[cmd.index("-c") + 1] == "0.98"

    for bad in ("high", "1.5", "0", "-0.5"):
        with pytest.raises(UserInputError, match="threshold"):
            mod.PanarooBuilder().call(
                _gffs(tmp_path), tmp_path / "bad", CoreGeneParams(extra={"threshold": bad}), _LOG
            )


def test_a_small_core_names_the_identity_gate(tmp_path: Path, monkeypatch, caplog) -> None:
    """Under 100 genes the warning points at the parameter that governs it."""
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo([]))  # the fake header holds 17 genes
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        mod.PanarooBuilder().call(_gffs(tmp_path), tmp_path / "o1", CoreGeneParams(), _LOG)
    assert "threshold=0.95" in caplog.text and "0.98" in caplog.text

    # Not repeated once the user has set it: they know about the parameter.
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        mod.PanarooBuilder().call(
            _gffs(tmp_path), tmp_path / "o2", CoreGeneParams(extra={"threshold": "0.9"}), _LOG
        )
    assert "threshold=0.95" not in caplog.text

    # Nor on a core that is large enough to stand on.
    caplog.clear()
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo([], genes=150))
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        mod.PanarooBuilder().call(_gffs(tmp_path), tmp_path / "o3", CoreGeneParams(), _LOG)
    assert caplog.text == ""


def test_no_core_at_all_refuses_and_names_the_identity_gate(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(mod, "run_tool", _fake_panaroo([], filtered=False, unfiltered=False))
    with pytest.raises(WorkdirError, match=r"least 0.98 identical"):
        mod.PanarooBuilder().call(
            _gffs(tmp_path), tmp_path / "o1", CoreGeneParams(filtered=False), _LOG
        )
    with pytest.raises(WorkdirError, match=r"You already set it to 0.9") as caught:
        mod.PanarooBuilder().call(
            _gffs(tmp_path),
            tmp_path / "o2",
            CoreGeneParams(filtered=False, extra={"threshold": "0.90"}),
            _LOG,
        )
    # The value is normalised, so the message and the command agree.
    assert "0.90" not in str(caught.value)


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
            gffs,
            tmp_path / "msa",
            tmp_path / "scratch",
            CoreGeneParams(extra={"threshold": "0.9"}),
            "fakecg",
            _LOG,
        )
    assert result.alignment == tmp_path / "msa" / CORE_GENE_FASTA and result.alignment.is_file()
    assert (tmp_path / "msa" / PANGENOME_DIR / "summary_statistics.txt").is_file()
    stats = json.loads((tmp_path / "msa" / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["core_genes"] == 7 and stats["alignment_length"] == 6 and stats["genomes"] == 3
    # The parameter that governs the count is kept beside it.
    assert stats["clustering_identity"] == "0.9"
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


class _IgnoringBuilder(_FakeBuilder):
    """A builder without Panaroo's cleaning or entropy filter, like ppanggolin."""

    capabilities = registry.get("ppanggolin").capabilities

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        result = super().call(gffs, out_dir, params, logger)
        result.clustering_identity = 0.8
        return result


def test_settings_a_builder_ignores_warn_and_stay_out_of_the_stats(
    tmp_path: Path, register_tool, caplog
) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS

    register_tool(registry, "fakeign", _IgnoringBuilder)
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        coregene_core(
            _gffs(tmp_path),
            tmp_path / "msa",
            tmp_path / "scratch",
            CoreGeneParams(clean_mode="sensitive"),
            "fakeign",
            _LOG,
        )
    assert "does not use clean_mode" in caplog.text
    stats = json.loads((tmp_path / "msa" / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert "clean_mode" not in stats and "filtered" not in stats
    # The identity the builder reports wins over an extra it does not read.
    assert stats["clustering_identity"] == 0.8


class _GapRowBuilder(_FakeBuilder):
    """Writes every genome, the last as a row of gaps (as PPanGGOLiN does)."""

    def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        aln = out_dir / "aln.fasta"
        rows = [f">{g.stem}\nACGTAC\nACG\n" for g in gffs[:-1]]
        rows.append(f">{gffs[-1].stem}\n------\n---\n")
        aln.write_text("".join(rows), encoding="utf-8")
        return CoreGeneResult(alignment=aln, tables=[], core_genes=7)


def test_gap_only_records_are_dropped_and_a_required_one_refused(
    tmp_path: Path, register_tool, caplog
) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS

    register_tool(registry, "fakegap", _GapRowBuilder)
    gffs = _gffs(tmp_path)
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        result, _ = coregene_core(
            gffs, tmp_path / "msa", tmp_path / "scratch", CoreGeneParams(), "fakegap", _LOG
        )
    text = result.alignment.read_text(encoding="utf-8")
    assert gffs[-1].stem not in text and text.count(">") == 2
    # Multi-line records survive intact.
    assert "ACGTAC\nACG\n" in text
    stats = json.loads((tmp_path / "msa" / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["genomes"] == 2 and stats["alignment_length"] == 9
    assert "Removed 1 record(s) of gaps only" in caplog.text
    assert "not in the core-gene alignment" in caplog.text
    with pytest.raises(WorkdirError, match="does not hold"):
        coregene_core(
            gffs,
            tmp_path / "m2",
            tmp_path / "s2",
            CoreGeneParams(),
            "fakegap",
            _LOG,
            required=(gffs[-1].stem,),
        )


def test_a_relative_scratch_is_resolved(tmp_path: Path, register_tool, monkeypatch) -> None:
    """phylo-build -o . hands a relative scratch to builders that set cwd inside it."""
    seen: list[Path] = []

    class _Seeing(_FakeBuilder):
        def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
            seen.append(out_dir)
            return super().call(gffs, out_dir, params, logger)

    register_tool(registry, "fakesee", _Seeing)
    gffs = _gffs(tmp_path)
    monkeypatch.chdir(tmp_path)
    coregene_core(gffs, Path("msa"), Path("scratch"), CoreGeneParams(), "fakesee", _LOG)
    assert seen and seen[0].is_absolute() and seen[0] == tmp_path / "scratch" / "fakesee"


def test_builder_settings_reach_the_stats(tmp_path: Path, register_tool) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS

    class _Settings(_FakeBuilder):
        def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
            result = super().call(gffs, out_dir, params, logger)
            result.settings = {"partition": "core"}
            return result

    register_tool(registry, "fakeset", _Settings)
    coregene_core(
        _gffs(tmp_path), tmp_path / "msa", tmp_path / "s", CoreGeneParams(), "fakeset", _LOG
    )
    stats = json.loads((tmp_path / "msa" / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["builder_settings"] == {"partition": "core"}


# --- per-gene alignments: selection, table and partitions --------------------


def _gene_builder(genes: dict[str, dict[str, str]], *, also_alignment: bool = False):
    """A builder that hands over one alignment per gene (records ``<genome>;<id>``)."""

    class _PerGene(_FakeBuilder):
        def call(self, gffs, out_dir, params, logger) -> CoreGeneResult:
            out_dir.mkdir(parents=True, exist_ok=True)
            found = []
            for gene, seqs in genes.items():
                path = out_dir / f"{gene}.aln.fas"
                path.write_text(
                    "".join(f">{g};{gene}_1\n{s}\n" for g, s in seqs.items()), encoding="utf-8"
                )
                found.append(GeneAlignment(gene, path))
            aln = None
            if also_alignment:
                aln = out_dir / "tool_concat.fasta"
                aln.write_text("".join(f">{g.stem}\nNNNN\n" for g in gffs), encoding="utf-8")
            return CoreGeneResult(alignment=aln, gene_alignments=found, core_genes=99)

    return _PerGene


_G1, _G2, _G3 = "Fam_Gen_sp_GCA_1", "Fam_Gen_sp_GCA_2", "Fam_Gen_sp_GCA_3"
_GENES = {
    # Every genome; column 1 (C, T, C) is variable but not informative.
    "famB": {_G1: "ACGT", _G2: "ATGT", _G3: "ACGT"},
    # Two genomes; no variation.
    "famA": {_G1: "GGG", _G2: "GGG"},
    # Every genome but short; column 1 (T, A, A) is variable only.
    "famC": {_G1: "TT", _G2: "TA", _G3: "TA"},
}


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        else:
            out[name] += line
    return out


def test_per_gene_alignments_are_selected_and_concatenated(
    tmp_path: Path, register_tool, caplog
) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS, CORE_GENES_TSV
    from repgenr.coregenes.select import Partition
    from repgenr.tree.partitions import PARTITIONS, PartitionSet, read_partitions

    register_tool(registry, "fakegenes", _gene_builder(_GENES, also_alignment=True))
    gffs = _gffs(tmp_path)
    msa = tmp_path / "msa"
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        result, _ = coregene_core(
            gffs, msa, tmp_path / "scratch", CoreGeneParams(), "fakegenes", _LOG
        )
    # The builder's own concatenation is ignored when it hands over genes;
    # the layout follows the gene ids.
    assert _read_fasta(result.alignment) == {
        _G1: "GGGACGTTT",
        _G2: "GGGATGTTA",
        _G3: "---ACGTTA",
    }
    assert result.core_genes == 3
    assert read_partitions(result.alignment) == PartitionSet(
        "dna", [Partition("famA", 1, 3), Partition("famB", 4, 7), Partition("famC", 8, 9)]
    )
    table = (msa / CORE_GENES_TSV).read_text(encoding="utf-8").splitlines()
    # Kept genes in rank order (occupancy, then length, then id).
    assert [ln.split("\t")[:3] for ln in table[1:]] == [
        ["famB", "true", "1"],
        ["famC", "true", "2"],
        ["famA", "true", "3"],
    ]
    stats = json.loads((msa / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["core_genes"] == 3 and stats["core_genes_available"] == 3
    assert stats["alignment_length"] == 9 and stats["genomes"] == 3
    assert stats["variable_sites_total"] == 2 and stats["informative_sites_total"] == 0
    assert stats["gene_table"] == CORE_GENES_TSV and stats["partitions"] == PARTITIONS
    assert stats["selection"] == {
        "max_genes": 0,
        "rank": "occupancy",
        "min_length": 0,
        "max_gaps": None,
        "seed": 0,
    }
    assert "3 genes available, 3 kept (rank occupancy); 9 columns" in caplog.text
    # The sidecar describes that very file: refused once the alignment changes.
    result.alignment.write_text(">x\nAAAA\n", encoding="utf-8")
    assert read_partitions(result.alignment) is None


def test_a_selection_caps_and_filters_the_genes(tmp_path: Path, register_tool) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS, CORE_GENES_TSV
    from repgenr.coregenes.select import GeneSelection, Partition
    from repgenr.tree.partitions import read_partitions

    register_tool(registry, "fakegenes", _gene_builder(_GENES))
    msa = tmp_path / "msa"
    result, _ = coregene_core(
        _gffs(tmp_path),
        msa,
        tmp_path / "scratch",
        CoreGeneParams(),
        "fakegenes",
        _LOG,
        selection=GeneSelection(max_genes=1, min_length=3),
    )
    assert _read_fasta(result.alignment) == {_G1: "ACGT", _G2: "ATGT", _G3: "ACGT"}
    found = read_partitions(result.alignment)
    assert found is not None and found.partitions == [Partition("famB", 1, 4)]
    rows = {
        ln.split("\t")[0]: ln.split("\t")
        for ln in (msa / CORE_GENES_TSV).read_text(encoding="utf-8").splitlines()[1:]
    }
    assert rows["famB"][1:3] == ["true", "1"] and rows["famB"][-2:] == ["1", "4"]
    assert rows["famA"][1:3] == ["false", ""] and rows["famA"][-2:] == ["", ""]
    stats = json.loads((msa / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["core_genes"] == 1 and stats["core_genes_available"] == 3
    assert stats["selection"]["max_genes"] == 1 and stats["selection"]["min_length"] == 3
    with pytest.raises(WorkdirError, match="None of the 3 core genes passes"):
        coregene_core(
            _gffs(tmp_path),
            tmp_path / "m2",
            tmp_path / "s2",
            CoreGeneParams(),
            "fakegenes",
            _LOG,
            selection=GeneSelection(min_length=10),
        )


def test_per_gene_gap_only_genome_is_dropped_and_a_required_one_refused(
    tmp_path: Path, register_tool, caplog
) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS
    from repgenr.tree.partitions import read_partitions

    # _G3 has only a gap record in the one gene: its concatenated row is gaps.
    genes = {"famA": {_G1: "GGG", _G2: "GGA", _G3: "---"}}
    register_tool(registry, "fakegenes", _gene_builder(genes))
    gffs = _gffs(tmp_path)
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        result, _ = coregene_core(
            gffs, tmp_path / "msa", tmp_path / "scratch", CoreGeneParams(), "fakegenes", _LOG
        )
    assert _read_fasta(result.alignment) == {_G1: "GGG", _G2: "GGA"}
    # The sidecar matches the published alignment, after the row was dropped.
    assert read_partitions(result.alignment) is not None
    stats = json.loads((tmp_path / "msa" / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["genomes"] == 2
    assert "Removed 1 record(s) of gaps only" in caplog.text
    assert "not in the core-gene alignment" in caplog.text
    with pytest.raises(WorkdirError, match="does not hold"):
        coregene_core(
            gffs,
            tmp_path / "m2",
            tmp_path / "s2",
            CoreGeneParams(),
            "fakegenes",
            _LOG,
            required=(_G3,),
        )


def test_per_gene_records_must_name_genomes_of_the_set(tmp_path: Path, register_tool) -> None:
    register_tool(registry, "fakegenes", _gene_builder({"famA": {_G1: "GG", "stranger": "GG"}}))
    with pytest.raises(WorkdirError, match="none of the genomes.*stranger"):
        coregene_core(
            _gffs(tmp_path), tmp_path / "msa", tmp_path / "s", CoreGeneParams(), "fakegenes", _LOG
        )


def test_the_legacy_path_is_unchanged_and_refuses_a_selection(
    tmp_path: Path, register_tool
) -> None:
    from repgenr.core.contracts import CORE_GENE_STATS, CORE_GENES_TSV
    from repgenr.coregenes.select import GeneSelection
    from repgenr.tree.partitions import PARTITIONS

    register_tool(registry, "fakecg", _FakeBuilder)
    msa = tmp_path / "msa"
    msa.mkdir()
    # Left by an earlier per-gene run; they describe another alignment.
    (msa / CORE_GENES_TSV).write_text("gene\n", encoding="utf-8")
    (msa / PARTITIONS).write_text("{}\n", encoding="utf-8")
    result, _ = coregene_core(
        _gffs(tmp_path), msa, tmp_path / "scratch", CoreGeneParams(), "fakecg", _LOG
    )
    assert result.core_genes == 7 and result.alignment is not None
    assert result.alignment.read_text(encoding="utf-8").count("ACGTAC") == 3
    assert not (msa / CORE_GENES_TSV).exists() and not (msa / PARTITIONS).exists()
    stats = json.loads((msa / CORE_GENE_STATS).read_text(encoding="utf-8"))
    assert stats["core_genes"] == 7 and stats["alignment_length"] == 6
    assert stats["genomes"] == 3 and stats["core_threshold"] == 0.95
    assert stats["core_genes_available"] is None
    assert stats["variable_sites_total"] is None and stats["informative_sites_total"] is None
    assert stats["gene_table"] is None and stats["partitions"] is None
    assert stats["selection"] == {
        "max_genes": 0,
        "rank": "occupancy",
        "min_length": 0,
        "max_gaps": None,
        "seed": 0,
    }
    with pytest.raises(UserInputError, match="returns no per-gene alignments"):
        coregene_core(
            _gffs(tmp_path),
            tmp_path / "m2",
            tmp_path / "s2",
            CoreGeneParams(),
            "fakecg",
            _LOG,
            selection=GeneSelection(max_genes=10),
        )
