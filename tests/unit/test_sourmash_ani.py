"""Sketch ANI: the shared-hash minimum below which a pair's ANI is unknown,
and its consumers (dereplicator, glance table, medoid, tree builder)."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np
import pytest

from repgenr.core import sourmash as core
from repgenr.core.errors import ToolExecutionError

_LOG = logging.getLogger("test-sourmash-ani")
_DESCRIBE_HEADER = (
    "signature_file,md5,ksize,moltype,num,scaled,n_hashes,seed,with_abundance,name,"
    "filename,license,sum_hashes\n"
)


def _describe_row(name: str, n_hashes: int, *, ksize: int = 31, scaled: int = 1000) -> str:
    md5 = (name.encode().hex() + "0" * 32)[:32]
    return f"x.sig,{md5},{ksize},DNA,0,{scaled},{n_hashes},42,0,{name},,CC0,{n_hashes}\n"


# --- the matrix guard ----------------------------------------------------------------

# a-b close; a-c about 5 shared hashes of 1150 each (sourmash ANI estimate
# 0.836795); b-c none.
_SIM = np.array([[1.0, 0.5, 0.002], [0.5, 1.0, 0.0], [0.002, 0.0, 1.0]])


def test_matrix_guard_marks_pairs_with_few_shared_hashes_unknown() -> None:
    assert core.MIN_SHARED_HASHES == 30
    ani, n_unknown = core.jaccard_to_ani_matrix(_SIM, 31, n_hashes=[1150] * 3)
    assert n_unknown == 2
    assert np.isnan(ani[0, 2]) and np.isnan(ani[2, 0]) and np.isnan(ani[1, 2])
    assert np.all(np.diag(ani) == 1.0)
    assert ani[0, 1] == pytest.approx((2 * 0.5 / 1.5) ** (1 / 31))


def test_matrix_without_sizes_is_unchanged() -> None:
    ani, n_unknown = core.jaccard_to_ani_matrix(_SIM, 31)
    assert n_unknown == 0 and not np.isnan(ani).any()
    assert ani[0, 2] == pytest.approx(0.836795, rel=1e-5)  # the sourmash estimate
    assert ani[1, 2] == 0.0 and np.all(np.diag(ani) == 1.0)


def test_matrix_guard_follows_the_size() -> None:
    """The same Jaccard on larger sketches is more shared hashes (a-c: 80)."""
    _ani, n_unknown = core.jaccard_to_ani_matrix(_SIM, 31, n_hashes=[20000] * 3)
    assert n_unknown == 1  # only the pair without a shared hash


@pytest.mark.parametrize(("shared", "unknown"), [(29, True), (30, False), (31, False)])
def test_matrix_guard_boundary(shared, unknown) -> None:
    """Sketches of 1000 hashes each sharing ``shared``: Jaccard s / (2000 - s)."""
    j = shared / (2000 - shared)
    sim = np.array([[1.0, j], [j, 1.0]])
    ani, n_unknown = core.jaccard_to_ani_matrix(sim, 21, n_hashes=[1000, 1000])
    assert n_unknown == int(unknown)
    assert bool(np.isnan(ani[0, 1])) is unknown
    # The threshold can be lowered.
    assert core.jaccard_to_ani_matrix(sim, 21, n_hashes=[1000, 1000], min_shared=10)[1] == 0


def test_matrix_guard_rejects_misaligned_sizes() -> None:
    with pytest.raises(ValueError, match="signature sizes"):
        core.jaccard_to_ani_matrix(_SIM, 31, n_hashes=[1150, 1150])


# --- the signature sizes -------------------------------------------------------------


def _fake_describe(rows: str):
    calls: list[list[str]] = []

    def run(caps, argv, **kwargs):
        argv = [str(a) for a in argv]
        calls.append(argv)
        Path(argv[argv.index("--csv") + 1]).write_text(_DESCRIBE_HEADER + rows, encoding="utf-8")

    return run, calls


def test_signature_sizes_by_label(tmp_path: Path) -> None:
    rows = (
        _describe_row("a", 1000)
        + _describe_row("b", 500, scaled=2000)
        + _describe_row("a", 9, ksize=21)  # another k is ignored
    )
    run, calls = _fake_describe(rows)
    fofn = tmp_path / "s.fofn"
    out = core.signature_sizes(
        core.SOURMASH_TOOL, fofn, tmp_path / "sizes.csv", 31, ["b", "a"], _LOG, run
    )
    assert out == ([500.0, 500.0], 2000)  # a brought to scaled=2000
    argv = calls[0]
    assert argv[:3] == ["sourmash", "sig", "describe"]
    assert argv[argv.index("-k") + 1] == "31" and argv[argv.index("--from-file") + 1] == str(fofn)


def test_signature_sizes_falls_back_with_one_warning(tmp_path: Path, caplog) -> None:
    run, _calls = _fake_describe(_describe_row("a", 1000))
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        assert (
            core.signature_sizes(
                core.SOURMASH_TOOL, tmp_path / "f", tmp_path / "s.csv", 31, ["a", "b"], _LOG, run
            )
            is None
        )

        def failing(caps, argv, **kwargs):
            raise ToolExecutionError([str(a) for a in argv], 1, "boom", tool="sourmash")

        assert (
            core.signature_sizes(
                core.SOURMASH_TOOL, tmp_path / "f", tmp_path / "s.csv", 31, ["a"], _LOG, failing
            )
            is None
        )
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 2
    assert all("without the shared-hash minimum" in w for w in warnings)


# --- the dereplicator's dense path ---------------------------------------------------


def _dense(tmp_path: Path, monkeypatch, describe: str | None):
    from repgenr.dereplicators import sourmash as sm

    genomes = []
    for name in ("g1.fasta", "g2.fasta", "g3.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        genomes.append(p)

    def fake_run_tool(caps, argv, **kwargs):
        argv = [str(a) for a in argv]
        if argv[:3] == ["sourmash", "sig", "describe"]:
            if describe is None:
                raise ToolExecutionError(argv, 1, "describe failed", tool="sourmash")
            Path(argv[argv.index("--csv") + 1]).write_text(
                _DESCRIBE_HEADER + describe, encoding="utf-8"
            )

    monkeypatch.setattr(sm, "run_tool", fake_run_tool)
    monkeypatch.setattr(
        sm, "_find_signatures", lambda sig_dir, gs: {g: tmp_path / f"{g.name}.sig" for g in gs}
    )
    # compare's column order differs from the file name order.
    labels = ["g3.fasta", "g1.fasta", "g2.fasta"]
    jaccard = np.array([[1.0, 0.002, 0.0], [0.002, 1.0, 0.5], [0.0, 0.5, 1.0]])
    monkeypatch.setattr(sm, "_read_compare_csv", lambda path: (labels, jaccard))
    out = tmp_path / "out"
    out.mkdir(exist_ok=True)
    return sm.SourmashDereplicator()._dense_ani_matrix(
        genomes, out, 31, 1000, _LOG, None, 1, too_large_hint=""
    )


def test_dense_path_marks_unknown_pairs(tmp_path: Path, monkeypatch, caplog) -> None:
    describe = "".join(_describe_row(n, 1150) for n in ("g1.fasta", "g2.fasta", "g3.fasta"))
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        labels, ani, _names = _dense(tmp_path, monkeypatch, describe)
    assert labels == ["g1.fasta", "g2.fasta", "g3.fasta"]
    assert np.isnan(ani[0, 2]) and np.isnan(ani[1, 2]) and not np.isnan(ani[0, 1])
    assert any(
        "2 of 3 pairs share fewer than 30 hashes at k=31" in r.getMessage() for r in caplog.records
    )


def test_dense_path_without_describe_converts_unguarded(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    with caplog.at_level(logging.WARNING, logger=_LOG.name):
        _labels, ani, _names = _dense(tmp_path, monkeypatch, None)
    assert not np.isnan(ani).any()
    assert ani[0, 2] == pytest.approx(0.836795, rel=1e-5)
    assert sum(r.levelno == logging.WARNING for r in caplog.records) == 1


def test_greedy_cluster_keeps_unknown_pairs_apart() -> None:
    from repgenr.dereplicators.sourmash import _greedy_cluster

    labels = ["a", "b", "c"]
    ani = np.array([[1.0, 0.995, np.nan], [0.995, 1.0, np.nan], [np.nan, np.nan, 1.0]])
    clusters, _status = _greedy_cluster(labels, ani, {x: x for x in labels}, 0.99)
    assert clusters == {"a": ["b"], "c": []}
    # NaN is below any threshold, the lowest included.
    clusters, _status = _greedy_cluster(labels, ani, {x: x for x in labels}, 0.0)
    assert clusters == {"a": ["b"], "c": []}


# --- glance and medoid ----------------------------------------------------------------


def test_glance_pairwise_csv_leaves_unknown_empty(tmp_path: Path) -> None:
    from repgenr.dereplicators.compare_io import write_pairwise_csv

    ani = np.array([[1.0, 0.99, np.nan], [0.99, 1.0, 0.98], [np.nan, 0.98, 1.0]])
    path = write_pairwise_csv(["a", "b", "c"], ani, tmp_path / "p.csv")
    rows = {
        (r["genome1"], r["genome2"]): r["similarity"]
        for r in csv.DictReader(path.open(encoding="utf-8"))
    }
    assert rows[("a", "c")] == "" and rows[("a", "b")] == "0.99"


def test_medoid_counts_unknown_as_zero() -> None:
    from repgenr.tree.medoid import pick_medoid

    ani = np.array([[1.0, 0.95, np.nan], [0.95, 1.0, 0.90], [np.nan, 0.90, 1.0]])
    name, mean = pick_medoid(["a", "b", "c"], ani)
    assert name == "b" and mean == pytest.approx(92.5)
    name, mean = pick_medoid(["a", "c"], np.array([[1.0, np.nan], [np.nan, 1.0]]))
    assert (name, mean) == ("a", 0.0)


# --- the tree builder -------------------------------------------------------------


def _tree(tmp_path: Path, monkeypatch, names: list[str], extra: dict, describe: bool = True):
    from repgenr.treebuilders import sourmash as tb
    from repgenr.treebuilders.base import TreeParams

    files = [f"{n}.fasta" for n in names]
    seen: dict = {}

    def fake_run_tool(caps, cmd, **k):
        cmd = [str(c) for c in cmd]
        if cmd[:2] == ["sourmash", "sketch"]:
            outdir = Path(cmd[cmd.index("--outdir") + 1])
            for line in (
                Path(cmd[cmd.index("--from-file") + 1]).read_text(encoding="utf-8").splitlines()
            ):
                (outdir / (Path(line).name + ".sig")).write_text("{}", encoding="utf-8")
        elif cmd[:2] == ["sourmash", "compare"]:
            body = "1,0.5,0.002\n0.5,1,0\n0.002,0,1\n"
            Path(cmd[cmd.index("--csv") + 1]).write_text(
                ",".join(files) + "\n" + body, encoding="utf-8"
            )
        elif cmd[:3] == ["sourmash", "sig", "describe"] and describe:
            rows = "".join(_describe_row(f, 1150) for f in files)
            Path(cmd[cmd.index("--csv") + 1]).write_text(_DESCRIBE_HEADER + rows, encoding="utf-8")
        return 0

    def fake_nj(labels, dist, out_dir, logger, *, backend, threads=1):
        seen["dist"] = dist
        return "(a,b,c);"

    monkeypatch.setattr(tb, "run_tool", fake_run_tool)
    monkeypatch.setattr(tb, "nj_tree", fake_nj)
    genomes = []
    for f in files:
        p = tmp_path / f
        p.write_text(">x\nACGT\n", encoding="utf-8")
        genomes.append(p)
    tb.SourmashBuilder().build(
        genomes, tmp_path / "out", TreeParams(extra={"nj": "python", **extra}), _LOG
    )
    return seen


def test_tree_places_unknown_pairs_at_distance_one(tmp_path: Path, monkeypatch, caplog) -> None:
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        seen = _tree(tmp_path, monkeypatch, ["a", "b", "c"], {})
    dist = seen["dist"]
    assert dist[0][2] == 1.0 and dist[1][2] == 1.0 and dist[0][1] < 0.1
    (line,) = [r.getMessage() for r in caplog.records if "distance is set to 1.0" in r.getMessage()]
    assert line.startswith("sourmash: 2 of 3 pairs share fewer than 30 hashes at k=31")
    # Not a multi-genus set: no k-mer size note.
    assert "ksize=21" not in caplog.text


def test_tree_without_sizes_keeps_the_estimate(tmp_path: Path, monkeypatch) -> None:
    seen = _tree(tmp_path, monkeypatch, ["a", "b", "c"], {}, describe=False)
    assert seen["dist"][0][2] == pytest.approx(1 - 0.836795, rel=1e-4)


def test_tree_recommends_k21_for_several_genera(tmp_path: Path, monkeypatch, caplog) -> None:
    names = [
        "Fam_Gen_alpha_GCA_000000001.1",
        "Fam_Gen_beta_GCA_000000002.1",
        "Fam_Other_gamma_GCA_000000003.1",
    ]
    (tmp_path / "k31").mkdir()
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        _tree(tmp_path / "k31", monkeypatch, names, {})
    notes = [r.getMessage() for r in caplog.records if "ksize=21" in r.getMessage()]
    assert len(notes) == 1 and "span 2 genera" in notes[0]
    caplog.clear()
    (tmp_path / "k21").mkdir()
    with caplog.at_level(logging.INFO, logger=_LOG.name):
        _tree(tmp_path / "k21", monkeypatch, names, {"ksize": 21})
    assert not any("ksize=21" in r.getMessage() for r in caplog.records)
