"""markers core-gene builder: search, GTDB-Tk's calls, the guard, the domain
choice, the match columns, and the phylo wiring."""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import pytest

from repgenr.annotators.base import AnnotateResult, Annotator
from repgenr.annotators.base import registry as annot_registry
from repgenr.core.context import WorkdirContext
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.core.plugins import ToolCapabilities
from repgenr.coregenes import markers as mod
from repgenr.coregenes.base import CoreGeneParams, registry
from repgenr.coregenes.markers_db import (
    INDEX_COLUMNS,
    load_marker_list,
    marker_cutoff,
    marker_source,
)
from repgenr.stages import phylo as phylo_mod
from repgenr.stages.annotate import AnnotateStageParams
from repgenr.stages.annotate import run as annotate_run
from repgenr.stages.coregene import check_alphabet, default_alphabet
from repgenr.stages.phylo import PhyloParams
from repgenr.treebuilders.base import InputKind, TreeBuilder
from repgenr.treebuilders.base import registry as tb_registry

_LOG = logging.getLogger("test-markers")
_ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
# Rich frames the usage error in a box under CI and may wrap inside a phrase.
_BOX = re.compile("[\u2500-\u257f]")


def _text(result) -> str:  # noqa: ANN001
    """The CLI output without ANSI codes or box lines, wrapped lines joined."""
    raw = result.output + str(result.exception or "")
    return " ".join(_BOX.sub(" ", _ANSI.sub("", raw)).split())


# Match states of every fake model.
_LENG = 6
BAC = load_marker_list("bac120")
AR = load_marker_list("ar53")


def make_db(root: Path, sets: tuple[str, ...] = ("bac120", "ar53")) -> Path:
    """A marker database with one small fake model per family."""
    db = root / "markers_db"
    if db.exists():
        return db
    (db / "hmms").mkdir(parents=True)
    rows = ["\t".join(INDEX_COLUMNS)]
    for name in sets:
        blocks = []
        for acc in load_marker_list(name):
            src = marker_source(acc)
            cut = marker_cutoff(src)
            blocks.append(
                f"HMMER3/f [3.4 | Aug 2023]\nNAME  n{acc}\nACC   {acc}.1\nLENG  {_LENG}\n"
                f"{cut.upper()}    10.00 10.00;\nHMM\n//\n"
            )
            rows.append("\t".join((name, acc, src, cut, f"n{acc}", f"{acc}.1", str(_LENG))))
        (db / "hmms" / f"{name}.hmm").write_text("".join(blocks), encoding="utf-8")
    (db / "markers.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (db / "version.json").write_text(json.dumps({"v": 1, "sets": list(sets)}), encoding="utf-8")
    return db


def _seq(i: int, j: int) -> str:
    """A protein longer than the models, distinct per (genome, marker)."""
    letters = "ACDEFGHIKLMNPQRSTVWY"
    return "M" + "".join(letters[(i * 7 + j * 3 + k) % 20] for k in range(9))


class Genomes:
    """Annotated genomes and the hits the fake hmmsearch reports for each.

    ``hits[genome]`` holds (gene, marker, score, evalue) rows.
    """

    def __init__(self, root: Path) -> None:
        self.root = root
        self.hits: dict[str, list[tuple[str, str, float, float]]] = {}
        self.proteins: dict[str, dict[str, str]] = {}
        self.gffs: list[Path] = []

    def add(self, markers: list[str], *, extra=(), proteins=None) -> Path:
        i = len(self.gffs) + 1
        name = f"Fam_Gen_sp_GCA_{i:06d}.1"
        prot = {f"G{i}_{m}": _seq(i, j) for j, m in enumerate(markers)}
        prot.update(proteins or {})
        self.proteins[name] = prot
        self.hits[name] = [(f"G{i}_{m}", m, 100.0, 1e-30) for m in markers] + list(extra)
        d = self.root / "ann" / name
        d.mkdir(parents=True)
        gff = d / f"{name}.gff3"
        gff.write_text(f"##gff-version 3\n##FASTA\n>{name}\nACGT\n", encoding="utf-8")
        (d / f"{name}.faa").write_text(
            "".join(f">{g} product\n{s}\n" for g, s in prot.items()), encoding="utf-8"
        )
        (d / f"{name}.ffn").write_text("".join(f">{g}\nATG\n" for g in prot), encoding="utf-8")
        self.gffs.append(gff)
        return gff


def _accessions(hmm: Path) -> set[str]:
    return {
        ln.split()[1].split(".")[0]
        for ln in hmm.read_text(encoding="utf-8").splitlines()
        if ln.startswith("ACC ")
    }


def fake_chain(genomes: Genomes, chains: list[list[list[str]]]):
    """run_chain that writes hmmsearch tables and hmmalign Pfam output."""

    def run_chain(caps, steps, **k):
        argvs = [[str(c) for c in argv] for _, argv in steps]
        chains.append(argvs)
        for argv in argvs:
            if argv[0] == "hmmsearch":
                table = Path(argv[argv.index("--tblout") + 1])
                hmm, faa = Path(argv[-2]), Path(argv[-1])
                accs = _accessions(hmm)
                rows = ["# target name  accession  query name  accession  E-value  score\n"]
                for gene, marker, score, evalue in genomes.hits[faa.stem]:
                    if marker in accs:
                        rows.append(
                            f"{gene:<20} -          n{marker:<12} {marker}.1  {evalue:g}  "
                            f"{score:.1f}   0.3   1e-30  {score:.1f}   0.3   1.1   1   0   0   "
                            "1   1   1   1 a protein\n"
                        )
                rows.append("# [ok]\n")
                table.write_text("".join(rows), encoding="utf-8")
            elif argv[0] == "hmmalign":
                dst = Path(argv[argv.index("-o") + 1])
                src = Path(argv[-1])
                lines = src.read_text(encoding="utf-8").splitlines()
                names, seqs = lines[0::2], lines[1::2]
                out = ["# STOCKHOLM 1.0", ""]
                for name, seq in zip(names, seqs, strict=True):
                    # One insert residue (lower case) before the match states.
                    row = seq[0].lower() + seq[1 : 1 + _LENG].upper().ljust(_LENG, "-") + "."
                    out.append(f"{name[1:]:<40} {row}")
                    out.append(f"#=GR {name[1:]} PP {'*' * len(row)}")
                out.append(f"#=GC RF{'':<34} .{'x' * _LENG}.")
                out.append("//")
                dst.write_text("\n".join(out) + "\n", encoding="utf-8")
            else:
                raise AssertionError(argv)

    return run_chain


@pytest.fixture
def chains() -> list[list[list[str]]]:
    return []


def _run(tmp_path, monkeypatch, genomes: Genomes, chains, *, db=None, **params):
    monkeypatch.setattr(mod, "run_chain", fake_chain(genomes, chains))
    kw = dict(threads=2, alphabet="protein", core_threshold=0.5)
    kw.update(params)
    return mod.MarkersBuilder().call(
        genomes.gffs,
        tmp_path / "out",
        CoreGeneParams(markers_db=db or make_db(tmp_path), **kw),
        _LOG,
    )


# --- the builder's declaration ---------------------------------------------


def test_registered_with_its_capabilities() -> None:
    cls = registry.get("markers")
    caps = cls.capabilities
    assert caps.container == "quay.io/biocontainers/hmmer:3.4--h7d74f8d_5"
    assert caps.conda == ("bioconda::hmmer",)
    assert [(b.name, b.min_version) for b in caps.required_binaries] == [
        ("hmmsearch", "3.3"),
        ("hmmalign", "3.3"),
    ]
    assert caps.accepted_extras == {"marker_set", "min_markers"}
    assert caps.ignored_params == {"clean_mode", "filtered"}
    assert not caps.auto_eligible and caps.recommended_max_genomes == 20000
    assert cls.required_alphabet == "protein"


# --- parsers ----------------------------------------------------------------


def test_read_tblout_takes_the_query_accession_without_version(tmp_path: Path) -> None:
    table = tmp_path / "t.tbl"
    table.write_text(
        "#                                                               --- full sequence ----\n"
        "# target name        accession  query name           accession    E-value  score\n"
        # A HMMER 3.4 row from the Wolbachia check.
        "MGHKIO_00066         -          Ribosomal_S9         PF00380.26   1.1e-45  151.2   0.3"
        "   1.4e-45  150.9   0.3   1.1   1   0   0   1   1   1   1 Small ribosomal subunit\n"
        # A model without an accession: the query name stands in.
        "MGHKIO_00243         -          TIGR00019            -           1.7e-146  483.9   7.4"
        "  1.9e-146  483.8   7.4   1.0   1   0   0   1   1   1   1 prfA\n"
        "#\n# [ok]\n",
        encoding="utf-8",
    )
    assert mod.read_tblout(table) == [
        mod.Hit("MGHKIO_00066", "PF00380", 1.1e-45, 151.2),
        mod.Hit("MGHKIO_00243", "TIGR00019", 1.7e-146, 483.9),
    ]
    table.write_text("g1 - q PF1.1 x 3\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="unreadable E-value"):
        mod.read_tblout(table)


def test_match_columns_keep_the_rf_match_states(tmp_path: Path) -> None:
    sto = tmp_path / "a.sto"
    # hmmalign 3.4 Pfam output (two Wolbachia prfA proteins), cut to 60 columns.
    a = "mdme-------NNLQDLKKKFNDIERNLENPTNL-SQKEFITLSKEYSELRPIIKIIDEYNTLKEEISD"[:60]
    b = "...n---------LQDLKKKFNDIERNLENPTNL-SQKEFITLSKEYSELRPIIKIIDEYNTLKEEISD"[:60]
    rf = "....xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"[:60]
    sto.write_text(
        f"# STOCKHOLM 1.0\n\ngenomeA         {a}\n#=GR genomeA PP {'9' * 60}\n"
        f"genomeB         {b}\n#=GC RF         {rf}\n//\n",
        encoding="utf-8",
    )
    out = mod.match_columns(sto)
    assert set(out) == {"genomeA", "genomeB"}
    width = 60 - 4
    assert {len(s) for s in out.values()} == {width}
    # The four insert columns (lower case and dots) are gone; gaps stay gaps.
    assert out["genomeA"].startswith("-------NNLQ")
    assert out["genomeB"].startswith("---------LQ")
    assert out["genomeA"] == out["genomeA"].upper()
    sto.write_text("# STOCKHOLM 1.0\ng1 ACGT\n//\n", encoding="utf-8")
    with pytest.raises(WorkdirError, match="RF"):
        mod.match_columns(sto)


# --- GTDB-Tk's calls ---------------------------------------------------------


def test_each_gene_keeps_its_best_model() -> None:
    hits = [
        mod.Hit("g1", "M1", 1e-10, 50.0),
        mod.Hit("g1", "M2", 1e-20, 80.0),  # higher score: g1 is M2's
        mod.Hit("g2", "M1", 1e-12, 60.0),
        mod.Hit("g3", "M3", 1e-5, 40.0),
        mod.Hit("g3", "M4", 1e-6, 40.0),  # same score, smaller E-value
        mod.Hit("g4", "M6", 1e-6, 40.0),
        mod.Hit("g4", "M5", 1e-6, 40.0),  # full tie: accession order
        mod.Hit("g9", "OTHER", 1e-99, 999.0),  # a model outside the set
    ]
    calls = mod.classify(["M1", "M2", "M3", "M4", "M5", "M6"], hits, lambda ids: {})
    assert calls == {
        "M1": ("unique", "g2", 60.0),
        "M2": ("unique", "g1", 80.0),
        "M3": ("missing", "", 0.0),
        "M4": ("unique", "g3", 40.0),
        "M5": ("unique", "g4", 40.0),
        "M6": ("missing", "", 0.0),
    }


def test_several_genes_of_a_marker_are_multiple_unless_identical() -> None:
    hits = [
        mod.Hit("a1", "SAME", 1e-9, 70.0),
        mod.Hit("a2", "SAME", 1e-9, 90.0),
        mod.Hit("b1", "DIFF", 1e-9, 70.0),
        mod.Hit("b2", "DIFF", 1e-9, 90.0),
        mod.Hit("c1", "ONE", 1e-9, 70.0),
    ]
    read: list[set[str]] = []

    def proteins(ids: set[str]) -> dict[str, str]:
        read.append(set(ids))
        return {"a1": "MKV", "a2": "MKV", "b1": "MKV", "b2": "MKL"}

    calls = mod.classify(["SAME", "DIFF", "ONE", "NONE"], hits, proteins)
    assert calls["SAME"] == ("multiple_unique", "a2", 90.0)  # the best hit is kept
    assert calls["DIFF"] == ("multiple", "", 0.0)  # treated as absent
    assert calls["ONE"] == ("unique", "c1", 70.0)
    assert calls["NONE"] == ("missing", "", 0.0)
    # Only the genes of markers with several hits are read.
    assert read == [{"a1", "a2", "b1", "b2"}]


# --- the builder end to end ---------------------------------------------------


def test_search_align_and_result(tmp_path: Path, monkeypatch, chains) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    g.add(BAC[:-1])  # one marker missing
    same = _seq(3, 0)
    g.add(
        BAC[1:],
        # The first marker twice with one sequence (multiple unique) and the
        # second twice with different sequences (multiple: absent).
        extra=[
            ("dupA", BAC[0], 90.0, 1e-20),
            ("dupB", BAC[0], 95.0, 1e-20),
            ("G3_x", BAC[1], 50.0, 1e-9),
        ],
        proteins={"dupA": same, "dupB": same, "G3_x": "MPPPPPPPPP"},
    )
    result = _run(tmp_path, monkeypatch, g, chains, threads=3)

    searches = [argv for chain in chains for argv in chain if argv[0] == "hmmsearch"]
    # Two searches per genome: the Pfam models with GA, the TIGRFAM ones with NC.
    assert len(searches) == 6
    for argv in searches:
        cut, hmm = argv[argv.index("--cpu") - 1], Path(argv[-2])
        accs = _accessions(hmm)
        if cut == "--cut_ga":
            assert hmm.name == "bac120_ga.hmm" and all(a.startswith("PF") for a in accs)
        else:
            assert cut == "--cut_nc" and hmm.name == "bac120_nc.hmm"
            assert all(a.startswith("TIGR") for a in accs)
        assert {"--noali", "--notextw"} <= set(argv) and argv[argv.index("--cpu") + 1] == "1"
    # One chain per genome at three genomes on three threads.
    assert [len(c) for c in chains if c[0][0] == "hmmsearch"] == [2, 2, 2]
    aligns = [argv for chain in chains for argv in chain if argv[0] == "hmmalign"]
    assert all(a[1:4] == ["--amino", "--outformat", "Pfam"] for a in aligns)

    names = [p.stem for p in g.gffs]
    summary = (tmp_path / "out" / mod.SUMMARY_TSV).read_text(encoding="utf-8").splitlines()
    assert summary[0].split("\t") == [
        "genome", "set", "unique", "multiple_unique", "multiple", "missing",
        "single_copy_percent",
    ]  # fmt: skip
    assert [ln.split("\t")[:6] for ln in summary[1:]] == [
        [names[0], "bac120", "120", "0", "0", "0"],
        [names[1], "bac120", "119", "0", "0", "1"],
        [names[2], "bac120", "118", "1", "1", "0"],
    ]
    per = (tmp_path / "out" / mod.PER_GENOME_TSV).read_text(encoding="utf-8").splitlines()
    assert per[0] == "genome\tset\tmarker\tstatus\tgene\tbitscore"
    assert f"{names[2]}\tbac120\t{BAC[0]}\tmultiple_unique\tdupB\t95.0" in per
    assert f"{names[2]}\tbac120\t{BAC[1]}\tmultiple\t\t" in per

    assert result.alignment is None and result.core_genes == 120
    assert [ga.gene for ga in result.gene_alignments] == BAC
    first = dict(
        r.split("\n")[:2]
        for r in result.gene_alignments[0].path.read_text(encoding="utf-8").split(">")[1:]
    )
    # Records named by genome; the kept duplicate is in; only match columns.
    assert set(first) == set(names)
    assert first[names[2]] == same[1 : 1 + _LENG]
    second = result.gene_alignments[1].path.read_text(encoding="utf-8")
    assert names[2] not in second  # multiple: absent
    assert [t.name for t in result.tables] == [mod.SUMMARY_TSV, mod.PER_GENOME_TSV]
    s = result.settings
    assert s["marker_set"] == "bac120" and s["marker_set_requested"] == "auto"
    assert s["min_markers"] == 12 and s["markers_in_set"] == 120
    # The version of the set used: adding or removing ar53 leaves it unchanged.
    assert s["markers_db_version"].startswith("bac120:")
    assert s["clustering_on"] == "GTDB marker HMMs"
    assert result.clustering_identity is None


def test_core_threshold_selects_by_single_copy_occupancy(
    tmp_path: Path, monkeypatch, chains
) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    g.add(BAC)
    g.add(BAC[:-1])
    r = _run(tmp_path, monkeypatch, g, chains, core_threshold=1.0)
    assert r.core_genes == 119 and BAC[-1] not in [ga.gene for ga in r.gene_alignments]
    r = _run(tmp_path, monkeypatch, g, chains, core_threshold=0.6)
    assert r.core_genes == 120


def test_genome_guard_refuses_and_warns(tmp_path: Path, monkeypatch, chains, caplog) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    g.add(BAC[:50])  # under half of 120: a warning
    g.add(BAC)
    caplog.set_level(logging.WARNING)
    _run(tmp_path, monkeypatch, g, chains, extra={"marker_set": "bac120"})
    assert "fewer than half of the 120 bac120 markers" in caplog.text
    assert f"{g.gffs[1].stem} (50)" in caplog.text

    g2 = Genomes(tmp_path / "b")
    g2.add(BAC)
    g2.add(BAC[:11])  # under 12 (10 % of 120, rounded up): refused
    g2.add(BAC)
    chains.clear()
    with pytest.raises(WorkdirError) as exc:
        _run(tmp_path / "b", monkeypatch, g2, chains, extra={"marker_set": "bac120"})
    msg = str(exc.value)
    assert f"{g2.gffs[1].stem} (11)" in msg and "fewer than 12 of the 120 bac120" in msg
    assert "wrong domain" in msg and "annotation" in msg
    # Refused before any alignment; the tables are there to look at.
    assert not any(c[0][0] == "hmmalign" for c in chains)
    assert (tmp_path / "b" / "out" / mod.SUMMARY_TSV).is_file()
    # min_markers lowers the bar.
    r = _run(
        tmp_path / "b", monkeypatch, g2, chains, extra={"marker_set": "bac120", "min_markers": 5}
    )
    assert r.settings["min_markers"] == 5


def test_auto_takes_bac120_when_half_the_genomes_reach_it(
    tmp_path: Path, monkeypatch, chains
) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    g.add(BAC)
    g.add(AR)  # archaeon-like: 5 bac120 markers (the shared ones)
    with pytest.raises(WorkdirError, match="wrong domain"):
        _run(tmp_path, monkeypatch, g, chains)
    # Only bac120 was searched: two of three genomes reach its minimum.
    hmms = {Path(a[-2]).name for c in chains for a in c if a[0] == "hmmsearch"}
    assert hmms == {"bac120_ga.hmm", "bac120_nc.hmm"}


def test_auto_switches_to_ar53_for_archaea(tmp_path: Path, monkeypatch, chains, caplog) -> None:
    g = Genomes(tmp_path)
    for _ in range(3):
        g.add(AR)
    caplog.set_level(logging.INFO)
    r = _run(tmp_path, monkeypatch, g, chains)
    assert r.settings["marker_set"] == "ar53" and r.settings["markers_in_set"] == 53
    assert r.settings["min_markers"] == 6  # 10 % of 53, rounded up
    assert r.core_genes == 53
    assert "so ar53 is used" in caplog.text
    hmms = {Path(a[-2]).name for c in chains for a in c if a[0] == "hmmsearch"}
    assert {"bac120_ga.hmm", "ar53_ga.hmm", "ar53_nc.hmm"} <= hmms
    summary = (tmp_path / "out" / mod.SUMMARY_TSV).read_text(encoding="utf-8")
    assert "\tar53\t53\t" in summary


def test_auto_keeps_bac120_when_the_database_has_no_ar53(
    tmp_path: Path, monkeypatch, chains, caplog
) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    g.add(BAC[:11])
    g.add(BAC[:11])
    caplog.set_level(logging.WARNING)
    with pytest.raises(WorkdirError, match="wrong domain"):
        _run(tmp_path, monkeypatch, g, chains, db=make_db(tmp_path, ("bac120",)))
    assert "holds no other set" in caplog.text


def test_refusals_before_any_search(tmp_path: Path, monkeypatch, chains) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    db = make_db(tmp_path / "x", ("bac120",))
    for params, err, match in (
        ({"alphabet": "dna"}, UserInputError, "aligns proteins"),
        ({"extra": {"marker_set": "bac100"}}, UserInputError, "marker_set must be one of"),
        ({"extra": {"min_markers": 0}}, UserInputError, "min_markers must be 1"),
        ({"extra": {"min_markers": "x"}}, UserInputError, "min_markers"),
        ({"extra": {"marker_set": "ar53"}, "db": db}, UserInputError, "--set ar53"),
        ({"extra": {"marker_set": "bac120", "min_markers": 121}}, UserInputError, "121"),
        ({"db": tmp_path / "nowhere"}, UserInputError, "markers-db download"),
    ):
        with pytest.raises(err, match=match):
            _run(tmp_path, monkeypatch, g, chains, **params)
    assert not any(a[0] == "hmmalign" for c in chains for a in c)


def test_missing_database_reads_the_environment(tmp_path: Path, monkeypatch, chains) -> None:
    g = Genomes(tmp_path)
    for _ in range(3):
        g.add(BAC)
    monkeypatch.setattr(mod, "run_chain", fake_chain(g, chains))
    monkeypatch.setenv("REPGENR_MARKERS_DB", str(make_db(tmp_path)))
    r = mod.MarkersBuilder().call(
        g.gffs, tmp_path / "out", CoreGeneParams(alphabet="protein"), _LOG
    )
    assert r.core_genes == 120


def test_an_incomplete_set_file_is_refused(tmp_path: Path, monkeypatch, chains) -> None:
    g = Genomes(tmp_path)
    g.add(BAC)
    db = make_db(tmp_path)
    hmm = db / "hmms" / "bac120.hmm"
    hmm.write_text(hmm.read_text(encoding="utf-8").split("//\n", 1)[1], encoding="utf-8")
    with pytest.raises(WorkdirError, match=f"lacks 1 bac120 model.*{BAC[0]}"):
        _run(tmp_path, monkeypatch, g, chains, db=db)


def test_default_min_markers() -> None:
    assert mod.default_min_markers(120) == 12
    assert mod.default_min_markers(53) == 6
    assert mod.default_min_markers(5) == 1


# --- the phylo stage ----------------------------------------------------------


class _MarkerAnnotator(Annotator):
    """Writes one protein per bac120 marker for every genome."""

    capabilities = ToolCapabilities(name="fakeannot_mk")

    def preflight(self):
        return {"fakeannot_mk": "1.0"}

    def call(self, genome, record_name, out_dir, params, logger) -> AnnotateResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        gff = out_dir / f"{record_name}.gff3"
        gff.write_text(f"##gff-version 3\n##FASTA\n>{record_name}\nACGT\n", encoding="utf-8")
        (out_dir / f"{record_name}.faa").write_text(
            "".join(f">{record_name}_{m}\n{_seq(1, j)}\n" for j, m in enumerate(BAC)),
            encoding="utf-8",
        )
        (out_dir / f"{record_name}.ffn").write_text(">x\nATG\n", encoding="utf-8")
        return AnnotateResult(gff=gff)


class _MsaTree(TreeBuilder):
    capabilities = ToolCapabilities(name="faketree_mk")
    input_kind = InputKind.MSA_FASTA
    alphabets: list[str] = []

    def preflight(self):
        return {"faketree_mk": "1.0"}

    def build(self, msa_or_genomes, out_dir, params, logger) -> Path:
        type(self).alphabets.append(params.alphabet)
        out_dir.mkdir(parents=True, exist_ok=True)
        names = [
            ln[1:].strip()
            for ln in Path(msa_or_genomes).read_text(encoding="utf-8").splitlines()
            if ln.startswith(">")
        ]
        tree = out_dir / "tree.nwk"
        tree.write_text("(" + ",".join(names) + ");\n", encoding="utf-8")
        return tree


class _AllHits:
    """Every genome hits every bac120 marker once (gene <record>_<marker>)."""

    def __getitem__(self, name: str) -> list[tuple[str, str, float, float]]:
        return [(f"{name}_{m}", m, 100.0, 1e-30) for m in BAC]


@pytest.fixture
def stage(workdir: Path, register_tool, monkeypatch, chains):
    register_tool(annot_registry, "fakeannot_mk", _MarkerAnnotator)
    register_tool(tb_registry, "faketree_mk", _MsaTree)
    _MsaTree.alphabets = []
    g = Genomes(workdir.parent)
    g.hits = _AllHits()  # type: ignore[assignment]
    monkeypatch.setattr(mod, "run_chain", fake_chain(g, chains))
    monkeypatch.setattr(mod.MarkersBuilder, "preflight", lambda self: {"hmmsearch": "3.4"})
    monkeypatch.delenv("REPGENR_MARKERS_DB", raising=False)
    reps = workdir / "derep" / "representatives"
    reps.mkdir(parents=True)
    for i in range(1, 4):
        (reps / f"Fam_gen_sp_GCA_00000{i}.fasta").write_text(f">s{i}\nACGT\n", encoding="utf-8")
    return WorkdirContext(workdir, create=True, logger=_LOG)


def _phylo(db: Path | None, **kw) -> PhyloParams:
    base = dict(
        treebuilder="faketree_mk",
        msa_source="coregene",
        coregene="markers",
        no_outgroup=True,
        markers_db=None if db is None else str(db),
    )
    base.update(kw)
    return PhyloParams(**base)


def test_the_alphabet_is_protein_and_dna_is_refused() -> None:
    assert default_alphabet("markers") == "protein"
    assert default_alphabet("panaroo") == "dna"
    assert phylo_mod._coregene_alphabet({}, "markers") == "protein"
    assert phylo_mod._coregene_alphabet({}, "ppanggolin") == "dna"
    check_alphabet("markers", "protein")
    with pytest.raises(UserInputError, match="'markers' aligns genes as protein only"):
        check_alphabet("markers", "dna")


def test_markers_tree_from_the_stage(stage, tmp_path: Path, caplog) -> None:
    from repgenr.core.doctor import _check_external_databases

    db = make_db(tmp_path / "dbs")
    annotate_run(stage, AnnotateStageParams(tool="fakeannot_mk"))
    caplog.set_level(logging.INFO)
    phylo_mod.run(stage, _phylo(db))
    assert "--coregene-alphabet protein is implied" in caplog.text
    assert _MsaTree.alphabets == ["protein"]
    rec = stage.config.stages["phylo"].params
    assert rec["coregene"] == "markers" and rec["coregene_alphabet"] == "protein"
    assert rec["markers_db"] == str(db.absolute())
    assert rec["markers_db_version"].startswith("bac120:")
    stats = json.loads((stage.phylo_msa_dir / "core_gene_stats.json").read_text("utf-8"))
    assert stats["alphabet"] == "protein" and stats["core_genes"] == 120
    assert stats["builder_settings"]["marker_set"] == "bac120"
    assert stats["builder_settings"]["markers_db_version"] == rec["markers_db_version"]
    assert "clean_mode" not in stats and "filtered" not in stats
    tables = sorted(p.name for p in (stage.phylo_msa_dir / "pangenome").iterdir())
    assert tables == [mod.PER_GENOME_TSV, mod.SUMMARY_TSV]
    partitions = json.loads((stage.phylo_msa_dir / "partitions.json").read_text("utf-8"))
    assert partitions["alphabet"] == "protein" and len(partitions["partitions"]) == 120
    # doctor finds the database the record names, at the recorded version.
    stage.save_config()
    findings = _check_external_databases(stage.workdir, stage.config)
    assert [f.level for f in findings] == ["ok"]
    # Other routes record neither field.
    assert phylo_mod._markers_record(_phylo(db, coregene="ppanggolin"), True) == {
        "markers_db": None,
        "markers_db_version": None,
    }


def test_the_record_resolves_the_database_from_the_environment(
    stage, tmp_path: Path, monkeypatch
) -> None:
    db = make_db(tmp_path / "dbs")
    monkeypatch.setenv("REPGENR_MARKERS_DB", str(db))
    annotate_run(stage, AnnotateStageParams(tool="fakeannot_mk"))
    phylo_mod.run(stage, _phylo(None))
    rec = stage.config.stages["phylo"].params
    assert rec["markers_db"] == str(db.absolute())
    assert rec["markers_db_version"].startswith("bac120:")


def test_precheck_refuses_a_missing_database_first(stage, tmp_path: Path) -> None:
    with pytest.raises(UserInputError, match="markers-db download"):
        phylo_mod.precheck(stage, _phylo(tmp_path / "none"))
    with pytest.raises(UserInputError, match="markers-db download"):
        phylo_mod.precheck(stage, _phylo(None))
    # With the database the next refusal is the missing annotation.
    with pytest.raises(UserInputError, match="repgenr annotate"):
        phylo_mod.precheck(stage, _phylo(make_db(tmp_path)))


def test_msa_key_carries_the_database_version_and_the_set(tmp_path: Path) -> None:
    db = make_db(tmp_path / "a")
    genomes = [tmp_path / "g1.fasta"]

    def key(db_path: Path, **extra) -> str:
        return phylo_mod._msa_key(genomes, None, _phylo(db_path, extra=extra), None, "ann")

    k = key(db)
    # A copy of the same models elsewhere gives the same key.
    import shutil

    shutil.copytree(db, tmp_path / "b" / "markers_db")
    assert key(tmp_path / "b" / "markers_db") == k
    # A newer model version changes it; so does the marker set.
    index = tmp_path / "b" / "markers_db" / "markers.tsv"
    index.write_text(index.read_text("utf-8").replace("PF00380.1", "PF00380.2"), "utf-8")
    assert key(tmp_path / "b" / "markers_db") != k
    assert key(db, marker_set="bac120") != k
    # With an explicit set, another set added to the database changes nothing.
    k_bac = key(db, marker_set="bac120")
    only = make_db(tmp_path / "c", ("bac120",))
    assert key(only, marker_set="bac120") == k_bac
    # auto may use either set, so its key covers every set the database holds.
    assert key(only) != k
    # min_markers is a builder extra and reaches the key too.
    assert key(db, min_markers=20) != k


def test_fingerprint_hashes_the_database_version_not_the_path(tmp_path: Path, monkeypatch) -> None:
    import shutil

    from repgenr.cli.base import _stage_fingerprint

    monkeypatch.delenv("REPGENR_MARKERS_DB", raising=False)
    db = make_db(tmp_path / "a")
    shutil.copytree(db, tmp_path / "b" / "markers_db")
    moved = tmp_path / "b" / "markers_db"

    def fp(p: PhyloParams) -> str:
        return _stage_fingerprint("phylo", p, {}, {})

    assert fp(_phylo(db, treebuilder="iqtree")) == fp(_phylo(moved, treebuilder="iqtree"))
    monkeypatch.setenv("REPGENR_MARKERS_DB", str(moved))
    assert fp(_phylo(None, treebuilder="iqtree")) == fp(_phylo(db, treebuilder="iqtree"))
    index = moved / "markers.tsv"
    index.write_text(index.read_text("utf-8").replace("PF00380.1", "PF00380.2"), "utf-8")
    assert fp(_phylo(moved, treebuilder="iqtree")) != fp(_phylo(db, treebuilder="iqtree"))
    # Another route does not see the path at all.
    assert fp(PhyloParams()) == fp(PhyloParams(markers_db=str(db)))


def test_phylo_build_refuses_a_missing_database_before_annotating(
    tmp_path: Path, monkeypatch
) -> None:
    from repgenr.stages.phylo import PhyloBuildParams, phylo_build

    monkeypatch.delenv("REPGENR_MARKERS_DB", raising=False)
    gd = tmp_path / "g"
    gd.mkdir()
    for i in range(3):
        (gd / f"G{i}.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    params = PhyloBuildParams(
        genomes_dir=gd, out_dir=tmp_path / "o", phylo=_phylo(None, treebuilder="iqtree")
    )
    with pytest.raises(UserInputError, match="markers-db download"):
        phylo_build(params, _LOG)
    assert not (tmp_path / "o" / "annotations").exists()


# --- the CLI --------------------------------------------------------------------


def test_cli_refusals(tmp_path: Path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    monkeypatch.delenv("REPGENR_MARKERS_DB", raising=False)
    gd = tmp_path / "g"
    gd.mkdir()
    for i in range(3):
        (gd / f"G{i}.fasta").write_text(">a\nACGT\n", encoding="utf-8")
    base = ["phylo-build", "--genomes-dir", str(gd), "-o", str(tmp_path / "o")]
    for args, match in (
        (["--markers-db", "x"], "--markers-db applies only with --msa-source coregene"),
        (
            ["--msa-source", "coregene", "--marker-set", "ar53"],
            "--marker-set applies only with --msa-source coregene --coregene markers",
        ),
        (["--msa-source", "coregene", "--coregene", "markers", "--marker-set", "x"], "x"),
        (
            ["--msa-source", "coregene", "--coregene", "markers", "--coregene-alphabet", "dna"],
            "aligns genes as protein only",
        ),
        (["--msa-source", "coregene", "--coregene", "markers"], "markers-db download"),
    ):
        out = CliRunner().invoke(app, [*base, *args])
        assert out.exit_code == 2, (args, out.output)
        assert match in _text(out), (args, out.output)
    wd = tmp_path / "wd"
    for cmd in (["phylo", "-wd", str(wd)], ["run", "-wd", str(wd), "-l", "genus", "--dry-run"]):
        out = CliRunner().invoke(app, [*cmd, "--marker-set", "bac120"])
        assert out.exit_code == 2, out.output
        assert "--marker-set applies only" in _text(out)


def test_auto_leaves_out_a_set_smaller_than_min_markers(
    tmp_path: Path, monkeypatch, chains, caplog
) -> None:
    g = Genomes(tmp_path)
    for _ in range(3):
        g.add(BAC)
    caplog.set_level(logging.WARNING)
    r = _run(tmp_path, monkeypatch, g, chains, extra={"min_markers": 60})
    assert r.settings["marker_set"] == "bac120" and r.settings["min_markers"] == 60
    assert "more than the 53 markers of ar53; auto does not consider ar53" in caplog.text
    # Refused before any search when no set is large enough.
    chains.clear()
    with pytest.raises(UserInputError, match="more than any marker set holds"):
        _run(tmp_path, monkeypatch, g, chains, extra={"min_markers": 121})
    with pytest.raises(UserInputError, match="the ar53 set has 53 markers"):
        _run(tmp_path, monkeypatch, g, chains, extra={"marker_set": "ar53", "min_markers": 54})
    assert chains == []


def test_fingerprint_ignores_the_database_without_a_core_gene_pass(
    tmp_path: Path, monkeypatch
) -> None:
    from repgenr.cli.base import _stage_fingerprint

    monkeypatch.delenv("REPGENR_MARKERS_DB", raising=False)
    db = make_db(tmp_path / "a")
    genomes_tree = _phylo(db, treebuilder="skani")

    def fp(p: PhyloParams) -> str:
        return _stage_fingerprint("phylo", p, {}, {})

    before = fp(genomes_tree)
    index = db / "markers.tsv"
    index.write_text(index.read_text("utf-8").replace("PF00380.1", "PF00380.2"), "utf-8")
    assert fp(genomes_tree) == before
    assert fp(_phylo(None, treebuilder="skani")) == before
