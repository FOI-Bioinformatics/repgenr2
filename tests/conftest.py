"""Shared pytest fixtures."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    # Registered here (an initial conftest under `testpaths`) so the option is
    # known from any invocation; the live suite's own conftest consumes it.
    parser.addoption(
        "--live-config",
        default=None,
        help="TOML mapping tools to bin dirs and setting cache_dir for tests/live "
        "(see tests/live/live.example.toml).",
    )


def pytest_runtest_setup(item: pytest.Item) -> None:
    """Skip tests marked ``requires_binary("name")`` when the tool is absent."""
    for marker in item.iter_markers(name="requires_binary"):
        for binary in marker.args:
            if shutil.which(binary) is None:
                pytest.skip(f"requires external binary '{binary}'")


class _ReachableSession:
    """Stands in for the reachability probe's session: every host answers."""

    def get(self, url, **kwargs):
        return None

    def close(self) -> None:
        pass


@pytest.fixture(autouse=True)
def _hermetic_reachability_probe(request: pytest.FixtureRequest, monkeypatch) -> None:
    """Keep the core.http reachability probe off the network outside the live suite.

    Stages probe the NCBI datasets host before calling ``datasets``; unit and
    integration tests fake ``datasets`` itself, so the probe answers as if the
    host were reachable. Tests of the probe replace ``_probe_session`` again.
    """
    if request.node.get_closest_marker("live"):
        return
    from repgenr.core import http

    monkeypatch.setattr(http, "_probe_session", _ReachableSession)


@pytest.fixture(autouse=True)
def _no_real_sourmash_sketches(request: pytest.FixtureRequest, monkeypatch) -> None:
    """The genome-writing stages sketch their genomes when sourmash can run.

    The development environment has sourmash on the PATH; outside the live
    suite the sketch step sees it as unavailable, so no test runs a real
    sourmash. Tests of the sketch step set availability themselves.
    """
    if request.node.get_closest_marker("live"):
        return
    from repgenr.core import sketches

    monkeypatch.setattr(sketches, "tool_available", lambda caps: False)


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    return tmp_path / "wd"


@pytest.fixture
def genome_files(workdir: Path) -> list[Path]:
    """Create a small genomes/ dir with RepGenR2-style filenames."""
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    names = [
        "Francisellaceae_francisella_tularensis_GCA_000001.fasta",
        "Francisellaceae_francisella_tularensis_GCA_000002.fasta",
        "Francisellaceae_francisella_tularensis_GCA_000003.fasta",
    ]
    out = []
    for i, name in enumerate(names):
        p = gdir / name
        p.write_text(f">seq{i}\n{'ACGT' * 10}\n")
        out.append(p)
    return out


@pytest.fixture()
def register_tool():
    """Register an adapter class on a registry for one test, restoring after."""
    registered: list[tuple] = []

    def _register(registry, name, cls):
        registry.register(name, cls, replace=True)
        registered.append((registry, name))
        return cls

    yield _register
    for registry, name in registered:
        registry.unregister(name)


@pytest.fixture()
def write_deliverables():
    """Create placeholder outputs for each recorded stage of a workdir.

    Tests that build repgenr.yaml by hand and check `status` need the
    declared deliverables on disk, or status reports the stages as stale.
    """
    from types import SimpleNamespace

    from repgenr.cli.base import STAGE_DELIVERABLES
    from repgenr.core.config import Config
    from repgenr.core.doctor import _layout

    def _write(wd: Path) -> None:
        ctx = _layout(wd)
        for name, record in Config.load(wd).stages.items():
            spec = STAGE_DELIVERABLES.get(name)
            for path in spec(ctx, SimpleNamespace(**record.params)) if spec else []:
                if path.suffix:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if not path.exists():
                        path.write_text("x\n", encoding="utf-8")
                else:
                    path.mkdir(parents=True, exist_ok=True)
                    (path / "placeholder.fasta").write_text(">x\nA\n", encoding="utf-8")

    return _write


class FakeSourmash:
    """Stands in for ``sourmash sketch``: writes the ``-o`` file, records calls.

    The file holds the signature name and the input path. ``fail`` names
    genomes (record names) whose call writes a partial file and then fails,
    as a tool killed mid-write would; ``interrupt`` names genomes whose call
    raises KeyboardInterrupt. ``barrier`` (a threading.Barrier) makes calls
    wait for each other, so an overlap is certain and ``peak`` measures it.

    A reads sketch (``-p ...,abund``, several FASTQ inputs) is recorded in
    ``reads_calls`` as (name, inputs) instead of ``calls``; ``fail`` applies
    to it by name. ``reads_hook``, when set, is called with the name before
    the file is written (a test synchronises with the assembler through it).
    """

    def __init__(self) -> None:
        import threading

        self.calls: list[str] = []
        self.reads_calls: list[tuple[str, list[str]]] = []
        self.reads_hook = None
        self.fail: set[str] = set()
        self.interrupt: set[str] = set()
        self.barrier: threading.Barrier | None = None
        self.peak = 0
        self._active = 0
        self._lock = threading.Lock()

    def run_tool(self, caps, command, *, logger, **kwargs) -> int:
        from repgenr.core.errors import ToolExecutionError

        argv = [str(c) for c in command]
        if argv[:3] == ["sourmash", "sig", "rename"]:
            # A rename keeps the file's input line and changes the name.
            source, name, out = Path(argv[3]), argv[4], Path(argv[argv.index("-o") + 1])
            with self._lock:
                self.calls.append(f"rename:{name}")
            first = source.read_text(encoding="utf-8").splitlines()[1]
            out.write_text(f"{name}\n{first}\n", encoding="utf-8")
            return 0
        name = argv[argv.index("--name") + 1]
        out = Path(argv[argv.index("-o") + 1])
        if argv[:5] == ["sourmash", "sketch", "dna", "-p", "k=21,k=31,k=51,scaled=1000,abund"]:
            inputs = argv[argv.index("-o") + 2 :]
            with self._lock:
                self.reads_calls.append((name, inputs))
            if self.reads_hook is not None:
                self.reads_hook(name)
            if name in self.fail:
                out.write_text("partial", encoding="utf-8")
                raise ToolExecutionError(argv, 1, output="killed", tool="sourmash")
            out.write_text("\n".join([name, *inputs]) + "\n", encoding="utf-8")
            return 0
        assert argv[:5] == ["sourmash", "sketch", "dna", "-p", "k=21,k=31,k=51,scaled=1000"]
        with self._lock:
            self.calls.append(name)
            self._active += 1
            self.peak = max(self.peak, self._active)
        try:
            if self.barrier is not None:
                self.barrier.wait(timeout=10)
            if name in self.interrupt:
                raise KeyboardInterrupt
            if name in self.fail:
                out.write_text("partial", encoding="utf-8")
                raise ToolExecutionError(argv, 1, output="killed", tool="sourmash")
            out.write_text(f"{name}\n{argv[-1]}\n", encoding="utf-8")
        finally:
            with self._lock:
                self._active -= 1
        return 0


@pytest.fixture
def fake_sourmash(monkeypatch) -> FakeSourmash:
    """sourmash available to the sketch step, run by :class:`FakeSourmash`."""
    from repgenr.core import sketches
    from repgenr.stages import sketch as sketch_stage

    fake = FakeSourmash()
    monkeypatch.setattr(sketches, "tool_available", lambda caps: True)
    monkeypatch.setattr(sketches, "preflight", lambda caps: {"sourmash": "4.9.4"})
    monkeypatch.setattr(sketch_stage, "preflight", lambda caps: {"sourmash": "4.9.4"})
    monkeypatch.setattr(sketches, "run_tool", fake.run_tool)
    return fake
