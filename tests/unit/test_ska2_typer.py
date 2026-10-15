"""ska2 reference-free typer: build a split k-mer file, then align variable sites."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.snptypers import ska2 as mod
from repgenr.snptypers.base import SnpParams

_LOG = logging.getLogger("test")


def _genomes(tmp_path: Path) -> list[Path]:
    out = []
    for name in ("a", "b"):
        p = tmp_path / f"{name}.fasta"
        p.write_text(">x\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


def _fake_run_tool(calls: list[list[str]], *, align_output: str = ">a\nA\n>b\nC\n"):
    def fake(caps, argv, **kw):  # noqa: ANN001
        argv = [str(a) for a in argv]
        calls.append(argv)
        if argv[1] == "build":
            Path(argv[argv.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif argv[1] == "align":
            Path(argv[argv.index("-o") + 1]).write_text(align_output, encoding="utf-8")
        return 0

    return fake


def test_ska2_builds_then_aligns(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_run_tool(calls))
    genomes = _genomes(tmp_path)

    result = mod.Ska2Typer().call(
        genomes,
        None,
        tmp_path / "out",
        SnpParams(threads=2, extra={"ksize": 21, "min_freq": 0.8}),
        _LOG,
    )

    build, align = calls
    assert build[:2] == ["ska", "build"]
    assert build[build.index("-k") + 1] == "21"
    assert build[build.index("--threads") + 1] == "2"
    # ska build takes a two-column list: sample name, then the FASTA path.
    listing = Path(build[build.index("-f") + 1])
    assert listing.read_text(encoding="utf-8").splitlines() == [
        f"a\t{genomes[0].resolve()}",
        f"b\t{genomes[1].resolve()}",
    ]
    assert align[:2] == ["ska", "align"]
    assert align[align.index("--min-freq") + 1] == "0.8"
    assert align[align.index("--filter") + 1] == "no-ambig-or-const"
    assert align[align.index("--threads") + 1] == "2"
    assert align[-1].endswith(".skf")
    assert result.core_snp_fasta.read_text(encoding="utf-8").startswith(">a")
    # Split k-mer variants have no positional whole-genome alignment.
    assert result.full_alignment is None
    assert result.masked is False


def test_ska2_ignores_the_reference_argument(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_run_tool(calls))
    genomes = _genomes(tmp_path)
    mod.Ska2Typer().call(genomes, genomes[0], tmp_path / "out", SnpParams(threads=1), _LOG)
    build, _align = calls
    # Every genome, the would-be reference included, is an ordinary sample.
    listing = Path(build[build.index("-f") + 1]).read_text(encoding="utf-8")
    assert listing.count("\n") == 2


def test_ska2_default_parameters(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_run_tool(calls))
    mod.Ska2Typer().call(_genomes(tmp_path), None, tmp_path / "out", SnpParams(threads=1), _LOG)
    build, align = calls
    assert build[build.index("-k") + 1] == "31"
    assert align[align.index("--min-freq") + 1] == "0.9"


def test_ska2_empty_alignment_is_an_error(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []
    monkeypatch.setattr(mod, "run_tool", _fake_run_tool(calls, align_output=""))
    with pytest.raises(WorkdirError, match="no variable"):
        mod.Ska2Typer().call(_genomes(tmp_path), None, tmp_path / "out", SnpParams(threads=1), _LOG)


def test_ska2_is_registered_as_reference_free() -> None:
    from repgenr.snptypers.base import registry

    assert registry.get("ska2").requires_reference is False


def test_ska2_lo_replaces_align(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake(caps, argv, **kw):  # noqa: ANN001
        argv = [str(a) for a in argv]
        calls.append(argv)
        if argv[1] == "build":
            Path(argv[argv.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif argv[1] == "lo":
            Path(argv[-1] + "_snps.fas").write_text(">a\nA\n>b\nC\n", encoding="utf-8")
            Path(argv[-1] + "_indels.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(mod, "run_tool", fake)
    result = mod.Ska2Typer().call(
        _genomes(tmp_path),
        None,
        tmp_path / "out",
        SnpParams(threads=3, extra={"lo": "true", "lo_missing": 0.2}),
        _LOG,
    )
    build, lo = calls
    assert lo[:2] == ["ska", "lo"] and lo[lo.index("--threads") + 1] == "3"
    assert lo[lo.index("-m") + 1] == "0.2" and lo[-2].endswith(".skf")
    assert result.core_snp_fasta.read_text(encoding="utf-8").startswith(">a")
    assert result.full_alignment is None
    assert result.vcf is not None and result.vcf.name == "lo_indels.vcf"
    assert not any(c[1] == "align" for c in calls)


def test_ska2_lo_with_only_headers_is_an_error(tmp_path: Path, monkeypatch) -> None:
    def fake(caps, argv, **kw):  # noqa: ANN001
        argv = [str(a) for a in argv]
        if argv[1] == "build":
            Path(argv[argv.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif argv[1] == "lo":
            Path(argv[-1] + "_snps.fas").write_text(">a\n\n>b\n\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(mod, "run_tool", fake)
    with pytest.raises(WorkdirError, match="no variable sites"):
        mod.Ska2Typer().call(
            _genomes(tmp_path), None, tmp_path / "out", SnpParams(extra={"lo": True}), _LOG
        )


def test_ska2_lo_without_output_is_an_error(tmp_path: Path, monkeypatch) -> None:
    def fake(caps, argv, **kw):  # noqa: ANN001
        argv = [str(a) for a in argv]
        if argv[1] == "build":
            Path(argv[argv.index("-o") + 1] + ".skf").write_bytes(b"skf")
        return 0

    monkeypatch.setattr(mod, "run_tool", fake)
    with pytest.raises(WorkdirError, match="ska lo"):
        mod.Ska2Typer().call(
            _genomes(tmp_path), None, tmp_path / "out", SnpParams(extra={"lo": True}), _LOG
        )
