"""The Nextflow layer with real tools: the local data-channel harness across
dereplicator and tree-builder variants, `main.nf` for both modes, and one
container profile run.

Nextflow is launched from the `repgenr_nf` environment named in the live
config (`[bin_dirs] nextflow`); the tools come from PATH as for the other
live tests. String parameters whose value starts with a dash are passed as `--key=value`
(Nextflow reads `--key --flag` as a boolean). Each run gets its own
`-work-dir` under the test's tmp path, and the
`test` profile caps process resources to what a laptop has (base.config asks
for 32 CPUs otherwise).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
NF = ROOT / "nextflow"

pytestmark = [pytest.mark.live, pytest.mark.requires_binary("nextflow")]


def _java_home() -> dict[str, str]:
    nf = shutil.which("nextflow")
    env = dict(os.environ)
    if nf and "JAVA_HOME" not in env:
        # conda-forge openjdk lives under <env>/lib/jvm; conda's activation
        # script sets JAVA_HOME there, which a bare PATH entry does not.
        root = Path(nf).resolve().parents[1]
        env["JAVA_HOME"] = str(root / "lib" / "jvm" if (root / "lib" / "jvm").is_dir() else root)
    env["PATH"] = os.pathsep.join([str(NF / "bin"), env.get("PATH", "")])
    env.setdefault("NXF_ANSI_LOG", "false")
    return env


@pytest.fixture
def run_nextflow(tmp_path: Path):
    def _run(script: Path, *args: str, timeout: float = 3600) -> subprocess.CompletedProcess[str]:
        outdir = tmp_path / "results"
        argv = [
            "nextflow",
            "run",
            str(script),
            "-c",
            str(NF / "nextflow.config"),
            "-work-dir",
            str(tmp_path / "work"),
            "--outdir",
            str(outdir),
            *args,
        ]
        proc = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            env=_java_home(),
            cwd=ROOT,
            timeout=timeout,
            check=False,
        )
        if proc.returncode != 0:
            tail = f"{proc.stdout[-4000:]}\n{proc.stderr[-4000:]}"
            pytest.fail(f"nextflow exited {proc.returncode}: {' '.join(argv)}\n{tail}")
        proc.outdir = outdir  # type: ignore[attr-defined]
        return proc

    return _run


def _leaves(newick: str) -> int:
    return newick.count(",") + 1 if newick.strip() else 0


@pytest.mark.parametrize(
    ("derep_tool", "phylo_args"),
    [
        pytest.param(
            "sourmash",
            "--treebuilder mashtree",
            marks=pytest.mark.requires_binary("sourmash", "mashtree"),
        ),
        pytest.param(
            "skder",
            "--treebuilder sourmash",
            marks=pytest.mark.requires_binary("skder", "sourmash"),
        ),
    ],
)
def test_local_dataflow_variants(
    run_nextflow, synthetic_set, derep_tool: str, phylo_args: str
) -> None:
    genomes = synthetic_set("clonal", n=8, length=50_000)
    proc = run_nextflow(
        NF / "tests" / "local_dataflow.nf",
        "-profile",
        "test",
        "--genomes_dir",
        str(genomes),
        "--derep_tool",
        derep_tool,
        "--phylo_args=" + phylo_args,
        "--derep_process_size",
        "4",
        "--tree2tax_args=--include-dereplicated",
    )
    out = proc.outdir
    reps = list((out / "dereplicate" / "merged" / "representatives").glob("*.fasta"))
    assert len(reps) == 3, "the clonal set has three clusters"
    assert _leaves((out / "phylo" / "tree" / "tree.nwk").read_text(encoding="utf-8")) == 3
    assert len((out / "genomes_map.tsv").read_text(encoding="utf-8").splitlines()) == 8
    assert "root" in (out / "tree2tax.tsv").read_text(encoding="utf-8")
    versions = (out / "pipeline_info" / "software_versions.yml").read_text(encoding="utf-8")
    assert derep_tool in versions and "repgenr" in versions


@pytest.mark.network
@pytest.mark.requires_binary("sourmash", "mashtree", "datasets")
def test_main_bacterial_test_profile(run_nextflow) -> None:
    proc = run_nextflow(
        NF / "main.nf",
        "-profile",
        "test",
        "--metadata_args=--source api -d rep -l genus -tg francisella --limit 5",
    )
    out = proc.outdir
    assert (out / "tree2tax.tsv").is_file() and (out / "genomes_map.tsv").is_file()
    assert (out / "phylo" / "tree" / "tree.nwk").is_file()
    assert len(list((out / "dereplicate" / "merged" / "representatives").glob("*.fasta"))) >= 2
    # The selection summary, its plots and the lineage table are published.
    meta = out / "metadata"
    for name in (
        "metadata_summary.tsv",
        "metadata_summary_number_in_level.png",
        "metadata_summary_number_per_level.png",
        "metadata_selected.tsv",
    ):
        assert (meta / name).is_file(), name
    summary = (meta / "metadata_summary.tsv").read_text(encoding="utf-8")
    assert "species\tFrancisella " in summary


@pytest.mark.network
@pytest.mark.requires_binary("sourmash", "mashtree", "datasets")
def test_main_viral_mode(run_nextflow) -> None:
    proc = run_nextflow(
        NF / "main.nf",
        "-profile",
        "test",
        "--mode",
        "viral",
        "--vmetadata_args=--target hepatovirus --complete-only",
        "--vgenome_args=-tg Hepatovirus --no-outgroup",
    )
    out = proc.outdir
    assert (out / "tree2tax.tsv").is_file()
    assert len(list((out / "dereplicate" / "merged" / "representatives").glob("*.fasta"))) >= 2
    # The selection summary names a viral species by its token, without the
    # genus in front ("Hepatovirus-ahepa", not "Hepatovirus Hepatovirus-ahepa").
    rows = [
        line.split("\t")
        for line in (out / "metadata" / "metadata_summary.tsv").read_text().splitlines()[1:]
    ]
    species = [name for level, name, _ in rows if level == "species"]
    assert species and not any(name.startswith("Hepatovirus Hepatovirus") for name in species)


@pytest.mark.container
@pytest.mark.requires_binary("docker", "skder", "FastTree")
def test_docker_profile_with_progressivemauve(run_nextflow, synthetic_set) -> None:
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        pytest.skip("docker daemon not running")
    genomes = synthetic_set("balanced", n=4, length=20_000)
    proc = run_nextflow(
        NF / "tests" / "local_dataflow.nf",
        "-profile",
        "test,docker",
        "--genomes_dir",
        str(genomes),
        "--derep_tool",
        "skder",
        "--derep_secondary_ani",
        "0.999",
        "--repgenr_opts=--container docker --platform linux/amd64",
        "--phylo_args=--aligner progressivemauve --treebuilder fasttree",
    )
    out = proc.outdir
    assert _leaves((out / "phylo" / "tree" / "tree.nwk").read_text(encoding="utf-8")) == 4


@pytest.mark.network
@pytest.mark.container
@pytest.mark.requires_binary("docker", "sourmash", "mashtree")
def test_main_reads_mode(run_nextflow) -> None:
    """--mode reads on two public runs (Illumina MiSeq and ONT) through the
    pinned assembler images: one READS_ASSEMBLE task per run, no QC database,
    then dereplication, a mashtree tree and tree2tax."""
    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        pytest.skip("docker daemon not running")
    # The standard profile: the test profile's 4 GB and 1 h caps are too tight
    # for a long-read assembly under emulation.
    proc = run_nextflow(
        NF / "main.nf",
        "-profile",
        "standard",
        "--mode",
        "reads",
        "--reads_args=--accession SRR25474756 --accession SRR28800588",
        "--repgenr_opts=--container docker --platform linux/amd64",
        "--derep_tool",
        "sourmash",
        timeout=7200,
    )
    out = proc.outdir
    rows = (out / "reads" / "selection.tsv").read_text(encoding="utf-8").splitlines()
    assert len(rows) == 3  # header + two assembled runs
    assert (out / "reads" / "assembly_stats.tsv").is_file()
    assert (out / "tree2tax.tsv").is_file() and (out / "genomes_map.tsv").is_file()
    assert _leaves((out / "phylo" / "tree" / "tree.nwk").read_text(encoding="utf-8")) == 2
    versions = (out / "pipeline_info" / "software_versions.yml").read_text(encoding="utf-8")
    assert "skesa" in versions and "flye" in versions
