"""SNP typers, recombination masking, ML tree builders and the workdir
tree2tax flags on the Francisella tularensis species set (ten genomes plus the
F. philomiragia outgroup, cached by the network fixtures)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from conftest import copy_workdir
from live_helpers import log_text, newick_leaves

from repgenr.core.config import Config
from repgenr.core.contracts import CORE_SNP_FASTA, GENOMES_MAP_TSV, TREE2TAX_TSV, TREE_NWK

pytestmark = [pytest.mark.live, pytest.mark.network, pytest.mark.requires_binary("sourmash")]


def _records(fasta: Path) -> int:
    return sum(1 for line in fasta.open(encoding="utf-8") if line.startswith(">"))


@pytest.fixture(scope="module")
def species_wd(species_cache: Path, tmp_path_factory, repgenr_cmd) -> Path:
    """A copy of the species set dereplicated to a few representatives."""
    import subprocess

    wd = copy_workdir(species_cache, tmp_path_factory.mktemp("species") / "wd")
    subprocess.run(
        [
            *repgenr_cmd,
            "dereplicate",
            "-wd",
            str(wd),
            "--tool",
            "sourmash",
            "-t",
            "4",
            "--target-reps",
            "4",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return wd


@pytest.fixture
def wd(species_wd: Path, tmp_path: Path) -> Path:
    return copy_workdir(species_wd, tmp_path / "wd")


def _genome_names(wd: Path) -> list[str]:
    return sorted(p.name for p in (wd / "genomes").glob("*.fasta"))


# --- snptype -------------------------------------------------------------------


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools")
def test_simple_typer_all_genomes_with_explicit_reference(run_repgenr, wd: Path) -> None:
    ref = _genome_names(wd)[1]
    run_repgenr(
        "snptype", "-wd", wd, "--tool", "simple", "--all-genomes", "--reference", ref, "-t", "4"
    )
    assert _records(wd / "snp" / CORE_SNP_FASTA) == 10, "--all-genomes types every genome"
    text = log_text(wd)
    assert "SNP typing 10 genomes with simple" in text
    assert ref in text and "No --reference given" not in text
    assert Config.load(wd).stages["snptype"].params["reference"] == ref


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools")
def test_simple_typer_on_representatives_only(run_repgenr, wd: Path) -> None:
    reps = {p.name for p in (wd / "derep" / "representatives").iterdir()}
    run_repgenr("snptype", "-wd", wd, "--tool", "simple", "-t", "4")
    assert _records(wd / "snp" / CORE_SNP_FASTA) == len(reps)
    assert "No --reference given" in log_text(wd), "the default reference is reported"


@pytest.mark.requires_binary("ska")
def test_ska2_typer_and_tool_arg(run_repgenr, wd: Path) -> None:
    run_repgenr(
        "--verbose",
        "snptype",
        "-wd",
        wd,
        "--tool",
        "ska2",
        "--all-genomes",
        "-t",
        "4",
        "--tool-arg",
        "ksize=21",
    )
    assert _records(wd / "snp" / CORE_SNP_FASTA) == 10
    assert re.search(r"ska build.*(-k 21|--k 21|-k=21)", log_text(wd)), "ksize reaches ska build"
    refused = run_repgenr(
        "snptype", "-wd", wd, "--tool", "ska2", "--all-genomes", "--mask", "gubbins", check=False
    )
    assert refused.returncode != 0 and "whole-genome alignment" in refused.stderr + refused.stdout


@pytest.mark.requires_binary("parsnp", "harvesttools")
def test_parsnp_typer(run_repgenr, wd: Path) -> None:
    run_repgenr("snptype", "-wd", wd, "--tool", "parsnp", "--all-genomes", "-t", "4")
    assert _records(wd / "snp" / CORE_SNP_FASTA) == 10
    assert Config.load(wd).stages["snptype"].tool_versions.get("parsnp")


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools", "run_gubbins.py")
def test_gubbins_mask_changes_the_core_alignment(run_repgenr, wd: Path) -> None:
    run_repgenr("snptype", "-wd", wd, "--tool", "simple", "--all-genomes", "-t", "4")
    unmasked = (wd / "snp" / CORE_SNP_FASTA).read_text(encoding="utf-8")
    run_repgenr(
        "snptype", "-wd", wd, "--tool", "simple", "--all-genomes", "--mask", "gubbins", "-t", "4"
    )
    masked = (wd / "snp" / CORE_SNP_FASTA).read_text(encoding="utf-8")
    assert _records(wd / "snp" / CORE_SNP_FASTA) == 10
    assert masked != unmasked, "the masked alignment differs from the unmasked one"
    rec = Config.load(wd).stages["snptype"]
    assert rec.params["mask"] == "gubbins"
    assert any("gubbins" in key for key in rec.tool_versions), rec.tool_versions


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools")
def test_snptype_allow_incomplete(run_repgenr, wd: Path) -> None:
    (wd / "genomes" / _genome_names(wd)[0]).unlink()
    refused = run_repgenr("snptype", "-wd", wd, "--tool", "simple", "--all-genomes", check=False)
    assert refused.returncode != 0
    run_repgenr(
        "snptype", "-wd", wd, "--tool", "simple", "--all-genomes", "--allow-incomplete", "-t", "4"
    )
    assert _records(wd / "snp" / CORE_SNP_FASTA) == 9


# --- phylo from SNPs --------------------------------------------------------------


def _tree(wd: Path) -> str:
    return (wd / "tree" / TREE_NWK).read_text(encoding="utf-8")


def _outgroup_leaf(wd: Path) -> str:
    return next((wd / "outgroup").glob("*.fasta")).stem


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools", "iqtree")
def test_iqtree_from_snptype_with_bootstrap_and_outgroup(run_repgenr, wd: Path) -> None:
    run_repgenr(
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "simple",
        "--all-genomes",
        "--treebuilder",
        "iqtree",
        "-B",
        "1000",
        "-t",
        "4",
    )
    tree = _tree(wd)
    leaves = newick_leaves(tree)
    assert len(leaves) == 11 and _outgroup_leaf(wd) in leaves
    assert re.search(r"\)\d+(\.\d+)?:", tree), "ultrafast bootstrap writes support labels"
    assert _records(wd / "tree" / "msa" / CORE_SNP_FASTA) == 11, "typed under tree/msa/"
    assert not (wd / "snp").exists(), "snp/ belongs to the snptype stage"
    assert Config.load(wd).stages["phylo"].params["bootstrap"] == 1000

    run_repgenr(
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "simple",
        "--all-genomes",
        "--treebuilder",
        "iqtree",
        "--no-outgroup",
        "-t",
        "4",
    )
    assert (
        _outgroup_leaf(wd) not in newick_leaves(_tree(wd)) and len(newick_leaves(_tree(wd))) == 10
    )


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools", "FastTree", "raxml-ng")
def test_fasttree_and_raxmlng_from_snptype(run_repgenr, wd: Path) -> None:
    run_repgenr(
        "--verbose",
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "simple",
        "--all-genomes",
        "--treebuilder",
        "fasttree",
        "-B",
        "200",
        "-t",
        "4",
    )
    assert len(newick_leaves(_tree(wd))) == 11
    assert "-boot 200" in log_text(wd), "--bootstrap reaches FastTree as its resample count"
    run_repgenr(
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "simple",
        "--all-genomes",
        "--treebuilder",
        "raxmlng",
        "-B",
        "50",
        "-t",
        "4",
    )
    assert len(newick_leaves(_tree(wd))) == 11
    assert "--bs-trees 50" in log_text(wd)


@pytest.mark.requires_binary("ska", "FastTree")
def test_ska2_source_with_reference_and_allow_incomplete(run_repgenr, wd: Path) -> None:
    ref = _genome_names(wd)[0]
    run_repgenr(
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "ska2",
        "--all-genomes",
        "--treebuilder",
        "fasttree",
        "--reference",
        ref,
        "-t",
        "4",
    )
    assert len(newick_leaves(_tree(wd))) == 11
    (wd / "genomes" / _genome_names(wd)[1]).unlink()
    refused = run_repgenr(
        "--force",
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "ska2",
        "--all-genomes",
        "--treebuilder",
        "fasttree",
        check=False,
    )
    assert refused.returncode != 0
    run_repgenr(
        "--force",
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "ska2",
        "--all-genomes",
        "--treebuilder",
        "fasttree",
        "--allow-incomplete",
        "-t",
        "4",
    )
    assert len(newick_leaves(_tree(wd))) == 10


@pytest.mark.requires_binary("minimap2", "samtools", "bcftools", "run_gubbins.py", "FastTree")
def test_phylo_mask_gubbins(run_repgenr, wd: Path) -> None:
    # Gubbins scans the ingroup only and its regions are masked in the full
    # alignment, so the F. philomiragia outgroup keeps its place (D-10).
    run_repgenr(
        "phylo",
        "-wd",
        wd,
        "--msa-source",
        "snptype",
        "--snptyper",
        "simple",
        "--all-genomes",
        "--mask",
        "gubbins",
        "--treebuilder",
        "fasttree",
        "-t",
        "4",
    )
    assert len(newick_leaves(_tree(wd))) == 11
    assert "run_gubbins.py" in log_text(wd), "the masker ran inside the phylo stage"
    assert Config.load(wd).stages["phylo"].completed


# --- tree2tax on a workdir -------------------------------------------------------


@pytest.mark.requires_binary("mashtree")
def test_tree2tax_workdir_flags(run_repgenr, wd: Path) -> None:
    run_repgenr("phylo", "-wd", wd, "--treebuilder", "mashtree", "-t", "4")
    reps = {p.stem for p in (wd / "derep" / "representatives").iterdir()}
    og = _outgroup_leaf(wd)

    run_repgenr("tree2tax", "-wd", wd, "--include-dereplicated")
    rel = (wd / TREE2TAX_TSV).read_text(encoding="utf-8")
    assert og in rel and "\troot\n" in rel
    mapped = (wd / GENOMES_MAP_TSV).read_text(encoding="utf-8").splitlines()
    assert len(mapped) == 11, "every genome plus the outgroup maps to a leaf"

    run_repgenr(
        "tree2tax",
        "-wd",
        wd,
        "--no-include-dereplicated",
        "--remove-outgroup",
        "--node-basename",
        "N",
        "--root-name",
        "top",
        "--collapse-support",
        "0.5",
        "--collapse-length",
        "0.5",
    )
    rel = (wd / TREE2TAX_TSV).read_text(encoding="utf-8")
    assert og not in rel and "\ttop\n" in rel and "\troot\n" not in rel
    children = {line.split("\t")[0] for line in rel.splitlines()[1:]}
    assert children >= reps, "every representative is a child"


@pytest.mark.requires_binary("mashtree", "flextaxd")
def test_tree2tax_tables_import_into_flextaxd(run_repgenr, wd: Path, tmp_path: Path) -> None:
    """The real FlexTaxD builds a database from the two tables and dumps every
    leaf back out; the unit test in tests/unit/test_flextaxd_contract.py
    mirrors its reader where the binary is absent."""
    import subprocess

    run_repgenr("phylo", "-wd", wd, "--treebuilder", "mashtree", "-t", "4")
    run_repgenr("tree2tax", "-wd", wd)
    db = tmp_path / "taxonomy.ftd"
    subprocess.run(
        [
            "flextaxd",
            "--taxonomy_file",
            str(wd / TREE2TAX_TSV),
            "--genomeid2taxid",
            str(wd / GENOMES_MAP_TSV),
            "--database",
            str(db),
            "--force",
        ],
        check=True,
    )
    out = tmp_path / "dump"
    out.mkdir()
    subprocess.run(
        ["flextaxd", "--database", str(db), "--dump", "--dump_descriptions", "--outdir", str(out)],
        check=True,
    )
    names = (out / "names.dmp").read_text(encoding="utf-8")
    leaves = {line.split("\t")[0] for line in (wd / TREE2TAX_TSV).read_text().splitlines()[1:]}
    reps = {p.stem for p in (wd / "derep" / "representatives").iterdir()}
    assert reps <= leaves and all(rep in names for rep in reps)
    assert len((wd / GENOMES_MAP_TSV).read_text(encoding="utf-8").splitlines()) == len(reps)
    params = Config.load(wd).stages["tree2tax"].params
    assert params.get("collapse_support") == 0.5
