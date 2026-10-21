"""glance with the sourmash backend, and --tool auto resolution.

sourmash is replaced by a fake ``run_tool``: ``sketch`` writes one signature
file per genome named in the fofn, and ``compare`` writes a Jaccard matrix
CSV with the genome paths as labels (the layout ``sourmash compare --csv``
produces). Auto resolution is tested against a bare registry of fake
adapters whose binaries are or are not on the PATH.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display in CI

import numpy as np  # noqa: E402
import pytest  # noqa: E402
from typer.testing import CliRunner  # noqa: E402

import repgenr.dereplicators.drep as drep_mod  # noqa: E402
import repgenr.dereplicators.sourmash as sm  # noqa: E402
from repgenr.cli.main import app  # noqa: E402
from repgenr.core.binaries import BinarySpec  # noqa: E402
from repgenr.core.config import Config  # noqa: E402
from repgenr.core.context import WorkdirContext  # noqa: E402
from repgenr.core.errors import MissingBinaryError, WorkdirError  # noqa: E402
from repgenr.core.plugins import Registry, ToolCapabilities  # noqa: E402
from repgenr.dereplicators.base import CompareResult, Dereplicator  # noqa: E402
from repgenr.stages import glance as glance_mod  # noqa: E402
from repgenr.stages.glance import GlanceParams, resolve_auto_tool  # noqa: E402
from repgenr.stages.glance import run as glance_run  # noqa: E402

_NAMES = ("a.fasta", "b.fasta", "c.fasta", "d.fasta", "e.fasta")
# Two groups: {a, b, c} and {d, e}; Jaccard 0.8 at k=31 is about 99.6 pct ANI.
_JACCARD = np.array(
    [
        [1.0, 0.8, 0.7, 0.01, 0.01],
        [0.8, 1.0, 0.75, 0.01, 0.01],
        [0.7, 0.75, 1.0, 0.01, 0.01],
        [0.01, 0.01, 0.01, 1.0, 0.85],
        [0.01, 0.01, 0.01, 0.85, 1.0],
    ]
)


def _fake_sourmash(calls: list[list[str]]):
    def run_tool(caps, argv, **kwargs):
        argv = [str(a) for a in argv]
        calls.append(argv)
        sub = argv[1]
        if sub == "sketch":
            fofn = Path(argv[argv.index("--from-file") + 1])
            outdir = Path(argv[argv.index("--outdir") + 1])
            for line in fofn.read_text(encoding="utf-8").split():
                (outdir / f"{Path(line).name}.sig").write_text("{}")
        elif sub == "compare":
            fofn = Path(argv[argv.index("--from-file") + 1])
            sigs = fofn.read_text(encoding="utf-8").split()
            labels = [s.removesuffix(".sig") for s in sigs]
            order = [_NAMES.index(Path(lab).name) for lab in labels]
            out = Path(argv[argv.index("--csv") + 1])
            with open(out, "w", encoding="utf-8", newline="") as fo:
                w = csv.writer(fo)
                w.writerow(labels)
                for i in order:
                    w.writerow([_JACCARD[i, j] for j in order])
        return 0

    return run_tool


def _setup(workdir: Path, names=_NAMES) -> WorkdirContext:
    ctx = WorkdirContext(workdir, create=True)
    ctx.genomes_dir.mkdir(parents=True)
    for name in names:
        (ctx.genomes_dir / name).write_text(">x\nACGT\n")
    return ctx


@pytest.fixture
def fake_sourmash(monkeypatch) -> list[list[str]]:
    calls: list[list[str]] = []
    monkeypatch.setattr(sm, "run_tool", _fake_sourmash(calls))
    monkeypatch.setattr(sm.SourmashDereplicator, "preflight", lambda self: {"sourmash": "4.9.4"})
    return calls


def test_glance_sourmash_writes_the_three_outputs(workdir: Path, fake_sourmash) -> None:
    ctx = _setup(workdir)
    out_pdf = glance_run(ctx, GlanceParams(tool="sourmash", threads=3, keep_files=True))

    for name in glance_mod.GLANCE_OUTPUTS:
        assert (ctx.workdir / name).is_file(), name
    assert out_pdf.read_bytes().startswith(b"%PDF")
    wd = ctx.workdir / "glance_wd"
    leaves = (wd / "dendrogram_leaves.txt").read_text().split()
    assert sorted(leaves) == sorted(_NAMES)
    # The two groups are contiguous in the dendrogram.
    assert set(leaves[:2]) == {"d.fasta", "e.fasta"} or set(leaves[-2:]) == {"d.fasta", "e.fasta"}
    values = glance_mod._pair_similarities(wd / "pairwise_ani.csv", 0.0, 1.0)
    assert len(values) == len(_NAMES) * (len(_NAMES) - 1) // 2


def test_glance_sourmash_uses_the_dereplication_sketch_and_ani(
    workdir: Path, fake_sourmash
) -> None:
    ctx = _setup(workdir)
    glance_run(ctx, GlanceParams(tool="sourmash", threads=3, keep_files=True))
    sketch = next(c for c in fake_sourmash if c[1] == "sketch")
    compare = next(c for c in fake_sourmash if c[1] == "compare")
    assert sketch[sketch.index("-p") + 1] == "k=31,scaled=1000"
    assert compare[compare.index("-k") + 1] == "31"
    assert compare[compare.index("--processes") + 1] == "3"
    rows = list(
        csv.DictReader((ctx.workdir / "glance_wd" / "pairwise_ani.csv").open(encoding="utf-8"))
    )
    ab = next(r for r in rows if {r["genome1"], r["genome2"]} == {"a.fasta", "b.fasta"})
    assert "--avg-containment" in compare
    # The fake matrix is read as average containment: ANI = c^(1/k).
    expected = sm._containment_to_ani_matrix(np.array([[0.8]]), 31)[0][0, 0]
    assert expected == pytest.approx(0.8 ** (1 / 31))
    assert float(ab["similarity"]) == pytest.approx(expected, abs=1e-6)


def test_glance_sourmash_record_and_plot_labels(workdir: Path, fake_sourmash, monkeypatch) -> None:
    from matplotlib.figure import Figure

    labels: dict[str, tuple[str, str, str]] = {}
    real_savefig = Figure.savefig

    def _capture(self, fname, *args, **kwargs):
        ax = self.axes[0]
        labels[Path(fname).name] = (ax.get_xlabel(), ax.get_ylabel(), ax.get_title())
        return real_savefig(self, fname, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", _capture)
    ctx = _setup(workdir)
    glance_run(ctx, GlanceParams(tool="sourmash", threads=2))

    hist_x, hist_y, title = labels["glance_MASH_ANI_similarity_histogram.png"]
    assert (hist_x, hist_y) == ("ANI", "Genome pairs")
    assert title == "ANI, all-vs-all (10 genome pairs)"
    assert labels["glance_MASH_ANI_similarity_boxplot.png"][1] == "ANI"
    record = ctx.config.stages["glance"]
    assert record.tool == "sourmash"
    assert record.tool_versions == {"sourmash": "4.9.4"}
    assert not (ctx.workdir / "glance_wd").exists()


def test_glance_sourmash_refuses_a_set_above_the_dense_limit(
    workdir: Path, fake_sourmash, monkeypatch
) -> None:
    monkeypatch.setattr(sm, "_DENSE_MAX_GENOMES", 3)
    ctx = _setup(workdir)
    with pytest.raises(WorkdirError, match="glance plots every genome pair"):
        glance_run(ctx, GlanceParams(tool="sourmash", threads=2))
    assert not fake_sourmash


def test_compare_supporters_include_sourmash() -> None:
    from repgenr.dereplicators.base import compare_supporters

    assert {"drep", "sourmash"} <= set(compare_supporters())
    assert "skder" not in compare_supporters()


# --- --tool auto --------------------------------------------------------------


class _Comparer(Dereplicator):
    def compare(self, genomes, out_dir, threads, logger) -> CompareResult:
        return CompareResult()

    def dereplicate(self, genomes, out_dir, params, logger):
        raise NotImplementedError


class _NoCompare(Dereplicator):
    def dereplicate(self, genomes, out_dir, params, logger):
        raise NotImplementedError


def _registry(spec: dict[str, tuple[bool, bool]]) -> Registry:
    """name -> (implements compare, binary on PATH)."""
    reg: Registry = Registry("test.dereplicators")
    reg._loaded = True  # bare registry: skip entry-point discovery
    for name, (compares, present) in spec.items():
        binary = "sh" if present else "repgenr-test-absent-binary"
        caps = ToolCapabilities(name=name, required_binaries=(BinarySpec(binary),))
        base = _Comparer if compares else _NoCompare
        reg.register(name, type(f"Fake_{name}", (base,), {"capabilities": caps}), replace=True)
    return reg


@pytest.mark.parametrize(
    ("spec", "expected"),
    [
        # dRep runs: it wins over sourmash and over others.
        ({"drep": (True, True), "sourmash": (True, True), "aaa": (True, True)}, "drep"),
        # dRep absent: sourmash.
        ({"drep": (True, False), "sourmash": (True, True), "aaa": (True, True)}, "sourmash"),
        # Neither: another compare-capable tool, alphabetically.
        ({"drep": (True, False), "sourmash": (True, False), "zzz": (True, True)}, "zzz"),
        # A tool without compare() is never picked, even when it runs.
        ({"drep": (True, False), "skder": (False, True)}, None),
        ({"drep": (True, False), "sourmash": (True, False)}, None),
    ],
)
def test_resolve_auto_tool(spec, expected) -> None:
    assert resolve_auto_tool(_registry(spec)) == expected


def test_resolve_auto_tool_counts_a_container_backend(monkeypatch) -> None:
    # Under an active container backend a declared image is enough, as for
    # dereplicate --tool auto.
    import repgenr.core.containers as containers

    reg: Registry = Registry("test.dereplicators")
    reg._loaded = True
    caps = ToolCapabilities(
        name="drep",
        container="quay.io/biocontainers/drep:3.4.5",
        required_binaries=(BinarySpec("repgenr-test-absent-binary"),),
    )
    reg.register("drep", type("Fake_drep", (_Comparer,), {"capabilities": caps}))

    class _Active:
        active = True

    assert resolve_auto_tool(reg) is None
    monkeypatch.setattr(containers, "get_config", lambda: _Active())
    assert resolve_auto_tool(reg) == "drep"


def test_auto_records_the_concrete_tool_and_resumes(
    workdir: Path, fake_sourmash, monkeypatch
) -> None:
    import repgenr.core.plugins as plugins

    # Hermetic: the real resolver runs, but availability is stubbed (dRep
    # absent, sourmash present) so the result does not depend on the PATH.
    monkeypatch.setattr(plugins, "tool_available", lambda caps: caps.name == "sourmash")
    _setup(workdir)
    runner = CliRunner()
    result = runner.invoke(app, ["glance", "-wd", str(workdir), "-t", "2"])
    assert result.exit_code == 0, result.output
    record = Config.load(workdir).stages["glance"]
    assert record.tool == "sourmash"
    assert record.params["tool"] == "sourmash"
    n_calls = len(fake_sourmash)
    # Same request, same resolved tool: the second run is a resume skip, and
    # naming the tool explicitly matches the auto fingerprint.
    result = runner.invoke(app, ["glance", "-wd", str(workdir), "-t", "2"])
    assert result.exit_code == 0 and "skipping" in result.output
    result = runner.invoke(app, ["glance", "-wd", str(workdir), "--tool", "sourmash"])
    assert result.exit_code == 0 and "skipping" in result.output
    assert len(fake_sourmash) == n_calls


def test_auto_prefers_drep_when_it_can_run(workdir: Path, monkeypatch) -> None:
    from tests.integration.test_glance import _fake_drep

    _setup(workdir, names=("a.fasta", "b.fasta"))
    monkeypatch.setattr(drep_mod.DrepDereplicator, "preflight", lambda self: {"dRep": "3.4.5"})
    monkeypatch.setattr(drep_mod, "run_tool", _fake_drep)
    monkeypatch.setattr(glance_mod, "resolve_auto_tool", lambda reg=None: "drep")
    result = CliRunner().invoke(app, ["glance", "-wd", str(workdir), "-t", "2"])
    assert result.exit_code == 0, result.output
    assert "glance --tool auto selected 'drep'" in result.output
    assert Config.load(workdir).stages["glance"].tool == "drep"


def test_auto_without_a_runnable_tool_exits_4_naming_the_supporters(
    workdir: Path, monkeypatch
) -> None:
    _setup(workdir, names=("a.fasta", "b.fasta"))
    monkeypatch.setattr(glance_mod, "resolve_auto_tool", lambda reg=None: None)
    result = CliRunner().invoke(app, ["glance", "-wd", str(workdir)])
    assert result.exit_code == 4, result.output
    assert "drep, sourmash" in result.output
    # As with a named tool that is missing (exit 4), the record stays
    # incomplete, so status reports the stage as not completed.
    record = Config.load(workdir).stages.get("glance")
    assert record is None or record.completed is None


def test_auto_without_a_tool_still_reports_workdir_problems_first(
    workdir: Path, monkeypatch
) -> None:
    monkeypatch.setattr(glance_mod, "resolve_auto_tool", lambda reg=None: None)
    runner = CliRunner()
    # Missing workdir: exit 3, not 4.
    assert runner.invoke(app, ["glance", "-wd", str(workdir / "absent")]).exit_code == 3
    # One genome: exit 3, not 4.
    _setup(workdir, names=("a.fasta",))
    result = runner.invoke(app, ["glance", "-wd", str(workdir)])
    assert result.exit_code == 3, result.output
    assert "at least two genomes" in result.output


def test_stage_resolves_auto_for_direct_callers(workdir: Path, fake_sourmash, monkeypatch) -> None:
    monkeypatch.setattr(glance_mod, "resolve_auto_tool", lambda reg=None: "sourmash")
    ctx = _setup(workdir)
    glance_run(ctx, GlanceParams(threads=2))  # default tool: auto
    assert ctx.config.stages["glance"].tool == "sourmash"
    assert ctx.config.stages["glance"].params["tool"] == "sourmash"

    monkeypatch.setattr(glance_mod, "resolve_auto_tool", lambda reg=None: None)
    with pytest.raises(MissingBinaryError, match="drep, sourmash"):
        glance_run(ctx, GlanceParams(threads=2))


def test_glance_cli_default_tool_is_auto() -> None:
    import typer.main

    cmd = typer.main.get_command(app).commands["glance"]  # type: ignore[attr-defined]
    opt = next(p for p in cmd.params if "--tool" in getattr(p, "opts", ()))
    assert opt.default == "auto"
    assert GlanceParams().tool == "auto"
    assert "auto" in (opt.help or "") and "sourmash" in (opt.help or "")


def test_log_auto_choice_names_the_tool(caplog) -> None:
    with caplog.at_level(logging.INFO, logger="t"):
        glance_mod.log_auto_choice(logging.getLogger("t"), "sourmash")
    assert "selected 'sourmash'" in caplog.text


def test_no_tool_message_without_a_container_backend() -> None:
    msg = glance_mod.no_compare_tool_message(_registry({"drep": (True, False)}))
    assert "none of drep is on the PATH" in msg
    assert "--container docker" in msg


def test_no_tool_message_under_an_active_container_backend(monkeypatch) -> None:
    # With a backend already active, suggesting one is no help: availability
    # there means a declared image.
    import repgenr.core.containers as containers

    class _Active:
        active = True

    monkeypatch.setattr(containers, "get_config", lambda: _Active())
    msg = glance_mod.no_compare_tool_message(_registry({"drep": (True, False)}))
    assert "container backend is active" in msg
    assert "none of drep declares a container image" in msg
    assert "is on the PATH" not in msg and "--container docker" not in msg
