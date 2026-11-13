"""Network stages on real data: GTDB (API and TSV), NCBI Datasets, NCBI
Virus, BV-BRC, and `run` end to end. Downloads are cached once under the
live cache directory; every test works on a copy.

Reference sets (2026-09-10): genus Francisella representatives from the GTDB
API (26 genomes, outgroup GCF_003574425.1); species Francisella tularensis
limited to 10 with the F. philomiragia representative GCF_000156715.1 as the
outgroup; hepatovirus from NCBI Virus (complete genomes, ~370 records).
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest
from live_helpers import log_text

from repgenr.core.config import Config
from repgenr.core.contracts import SELECTION_TSV, TREE2TAX_TSV, list_fasta, read_selection

pytestmark = [pytest.mark.live, pytest.mark.network]

from conftest import GENUS_OUTGROUP, GTDB_RELEASE, PHILOMIRAGIA_REP  # noqa: E402


def _rows(wd: Path):
    rows = read_selection(wd / SELECTION_TSV)
    return [r for r in rows if not r.is_outgroup], [r for r in rows if r.is_outgroup]


# --- NCBI (datasets genome summary) ---------------------------------------------


def test_ncbi_refseq_species_selection_and_download(run_repgenr, tmp_path: Path) -> None:
    """--source ncbi selects RefSeq assemblies by NCBI taxonomy, with an
    outgroup from the reference genomes of the genus, and genome downloads them."""
    wd = tmp_path / "ncbi"
    run_repgenr(
        "metadata", "-wd", wd, "--source", "ncbi", "-d", "all", "-l", "species",
        "-tg", "Francisella", "-ts", "philomiragia", "--assembly-level", "complete",
        "--limit", "4",
    )  # fmt: skip
    ingroup, outgroup = _rows(wd)
    assert len(ingroup) == 4 and {r.species for r in ingroup} == {"philomiragia"}
    assert all(r.accession.startswith("GCF_") for r in [*ingroup, *outgroup])
    assert all(r.completeness is not None for r in ingroup)
    assert len(outgroup) == 1 and outgroup[0].genus == "Francisella"
    assert outgroup[0].species != "philomiragia"

    record = Config.load(wd).stages["metadata"]
    assert record.tool == "ncbi-datasets" and "datasets" in (record.tool_versions or {})
    assert record.params["query_date"] and record.params["assembly_level"] == "complete"
    assert record.params["outgroup_rule"] == "quality"
    assert (wd / "ncbi_genomes.tsv").is_file()

    run_repgenr("genome", "-wd", wd)
    assert len(list_fasta(wd / "genomes")) == 4
    assert len(list_fasta(wd / "outgroup")) == 1
    assert (wd / "missing_accessions.txt").read_text().strip() == ""


def test_ncbi_assembly_sources_and_filters(run_repgenr, tmp_path: Path) -> None:
    """genbank selects GCA_ only; all lists every assembly once; the filters
    reach NCBI; a name NCBI does not know is refused."""
    species = ("-l", "species", "-tg", "Francisella", "-ts", "philomiragia", "--no-outgroup")
    gb = tmp_path / "genbank"
    run_repgenr(
        "metadata", "-wd", gb, "--source", "ncbi", "-d", "all", *species,
        "--assembly-source", "genbank", "--exclude-mag", "--released-before", "2020-01-01",
    )  # fmt: skip
    genbank, _ = _rows(gb)
    assert genbank and all(r.accession.startswith("GCA_") for r in genbank)

    both = tmp_path / "all"
    run_repgenr(
        "metadata", "-wd", both, "--source", "ncbi", "-d", "all", *species,
        "--assembly-source", "all", "--include-atypical",
    )  # fmt: skip
    rows, _ = _rows(both)
    assemblies = [r.accession[4:] for r in rows]
    assert len(assemblies) == len(set(assemblies)), "an assembly is listed as GCF_ and GCA_"
    assert len(rows) >= len(genbank)
    assert any(r.accession.startswith("GCF_") for r in rows)

    proc = run_repgenr(
        "metadata", "-wd", tmp_path / "bad", "--source", "ncbi", "-d", "all", "-l", "genus",
        "-tg", "Nonexistentus", check=False,
    )  # fmt: skip
    assert proc.returncode == 2, proc.stderr
    assert "no taxon named" in proc.stdout + proc.stderr


# --- GTDB API ----------------------------------------------------------------


def test_api_genus_representatives(genus_cache: Path) -> None:
    ingroup, outgroup = _rows(genus_cache)
    assert len(ingroup) == 26 and {r.genus for r in ingroup} == {"Francisella"}
    assert [r.accession for r in outgroup] == [GENUS_OUTGROUP]
    assert all(r.completeness is not None for r in ingroup), "API cards carry CheckM quality"
    assert len(list_fasta(genus_cache / "genomes")) == 26
    assert (genus_cache / "outgroup").is_dir()
    assert (genus_cache / "outgroup_accession.txt").read_text().strip() == GENUS_OUTGROUP


def test_api_species_limit_and_explicit_outgroup(species_cache: Path) -> None:
    ingroup, outgroup = _rows(species_cache)
    assert len(ingroup) == 10 and {r.species for r in ingroup} == {"tularensis"}
    assert [r.accession for r in outgroup] == [PHILOMIRAGIA_REP]
    assert "philomiragia" in outgroup[0].species
    assert len(list_fasta(species_cache / "genomes")) == 10
    assert "--limit" in log_text(species_cache) or "limit" in log_text(species_cache).lower()


def test_api_several_species_without_outgroup(run_repgenr, tmp_path: Path) -> None:
    """-ts given twice (and as a comma list) selects both species; --no-outgroup
    leaves no outgroup row and no outgroup_accession.txt."""
    wd = tmp_path / "two"
    run_repgenr(
        "metadata", "-wd", wd, "-d", "rep", "-l", "species", "-tg", "Francisella",
        "-ts", "tularensis", "-ts", "philomiragia", "--source", "api", "--no-outgroup",
    )  # fmt: skip
    ingroup, outgroup = _rows(wd)
    species = {r.species for r in ingroup}
    assert {"tularensis", "philomiragia"} <= {s.split("_")[0] for s in species}, species
    assert outgroup == []
    assert not (wd / "outgroup_accession.txt").exists()
    record = Config.load(wd).stages["metadata"].params
    assert record["no_outgroup"] is True and record["outgroup"] is None

    # With an outgroup it lies outside both species.
    wd2 = tmp_path / "two_og"
    run_repgenr(
        "metadata", "-wd", wd2, "-d", "rep", "-l", "species", "-tg", "Francisella",
        "-ts", "tularensis,philomiragia", "--source", "api",
    )  # fmt: skip
    ingroup2, outgroup2 = _rows(wd2)
    assert {r.accession for r in ingroup2} == {r.accession for r in ingroup}
    assert len(outgroup2) == 1
    assert outgroup2[0].species.split("_")[0] not in {"tularensis", "philomiragia"}


def test_api_target_family_widens_the_selection(run_repgenr, tmp_path: Path) -> None:
    wd = tmp_path / "fam"
    run_repgenr(
        "metadata",
        "-wd",
        wd,
        "-d",
        "rep",
        "-l",
        "family",
        "-tf",
        "Francisellaceae",
        "--source",
        "api",
    )
    ingroup, outgroup = _rows(wd)
    assert len(ingroup) > 26 and len({r.genus for r in ingroup}) > 1
    assert len(outgroup) == 1 and outgroup[0].family != "Francisellaceae"


# --- GTDB TSV ------------------------------------------------------------------


def test_tsv_source_downloads_and_parses_the_release_table(tsv_cache: Path) -> None:
    table = tsv_cache / f"bac120_metadata_r{GTDB_RELEASE.split('.')[0]}.tsv.gz"
    assert table.is_file() and table.stat().st_size > 100_000_000
    ingroup, outgroup = _rows(tsv_cache)
    assert len(ingroup) >= 20 and len(outgroup) == 1
    assert Config.load(tsv_cache).stages["metadata"].params["release"] == GTDB_RELEASE


def test_tsv_nodownload_and_metadata_path_reuse_the_table(
    run_repgenr, tsv_cache: Path, copy_of, tmp_path: Path
) -> None:
    table = next(tsv_cache.glob("bac120_metadata_r*.tsv.gz"))
    wd = copy_of(tsv_cache)
    run_repgenr(
        "--force",
        "metadata",
        "-wd",
        wd,
        "-d",
        "rep",
        "-l",
        "genus",
        "-tg",
        "Francisella",
        "--source",
        "tsv",
        "-r",
        GTDB_RELEASE,
        "--gtdb-version",
        "bac120",
        "--nodownload",
    )
    assert "Using previously downloaded" in log_text(wd)

    fresh = tmp_path / "fresh"
    run_repgenr(
        "metadata",
        "-wd",
        fresh,
        "-d",
        "rep",
        "-l",
        "genus",
        "-tg",
        "Francisella",
        "--source",
        "tsv",
        "-r",
        GTDB_RELEASE,
        "--gtdb-version",
        "bac120",
        "--metadata-path",
        table,
    )
    assert "Using provided metadata" in log_text(fresh)
    assert {r.accession for r in _rows(fresh)[0]} == {r.accession for r in _rows(tsv_cache)[0]}


def test_tsv_outgroup_modes_choose_from_the_gtdb_tree(
    run_repgenr, tsv_cache: Path, copy_of, tmp_path: Path
) -> None:
    """Francisella (GTDB r232): --outgroup-mode sister takes a genome near the
    genus in the GTDB tree, quality the best genome of the family, first the
    first representative in table order."""
    import csv

    from repgenr.core.outgroup_choice import OUTGROUP_CHOICE_TSV

    def choice(wd: Path) -> list[dict]:
        with open(wd / OUTGROUP_CHOICE_TSV, encoding="utf-8", newline="") as fo:
            return list(csv.DictReader(fo, delimiter="\t"))

    base = ["metadata", "-d", "rep", "-l", "genus", "-tg", "Francisella", "--source", "tsv"]
    base += ["-r", GTDB_RELEASE, "--gtdb-version", "bac120", "--nodownload"]
    wd = copy_of(tsv_cache)
    run_repgenr("--force", *base[:1], "-wd", wd, *base[1:], "--outgroup-mode", "sister")
    tree = wd / f"bac120_r{GTDB_RELEASE.split('.')[0]}.tree"
    assert tree.is_file()
    rows = choice(wd)
    assert rows[0]["chosen"] == "yes" and rows[0]["mode"] == "sister"
    # The sister group of Francisella in r232 is one metagenome-assembled
    # genome; an isolate of the next group out is taken instead.
    assert rows[0]["accession"] == "GCF_000379445.1"  # Fangia hongkongensis
    assert rows[0]["group"] == "2" and rows[0]["isolate"] == "yes"
    assert "metagenome-assembled" in rows[0]["reason"]
    assert float(rows[0]["tree_distance"]) > 0
    sister = _rows(wd)[1][0].accession
    assert sister == rows[0]["accession"]

    # The same choice from a tree given by path, in a fresh workdir.
    fresh = tmp_path / "fresh"
    table = next(tsv_cache.glob("bac120_metadata_r*.tsv.gz"))
    run_repgenr(
        *base[:1], "-wd", fresh, *base[1:-1], "--metadata-path", table, "--gtdb-tree-path", tree
    )
    assert not (fresh / tree.name).exists()
    assert _rows(fresh)[1][0].accession == sister

    run_repgenr("--force", *base[:1], "-wd", wd, *base[1:], "--outgroup-mode", "quality")
    assert choice(wd)[0]["mode"] == "quality"
    run_repgenr("--force", *base[:1], "-wd", wd, *base[1:], "--outgroup-mode", "first")
    assert not (wd / OUTGROUP_CHOICE_TSV).exists()


# --- genome ---------------------------------------------------------------------


def test_genome_accession_list_only_is_a_pure_query(run_repgenr, species_cache, copy_of) -> None:
    wd = copy_of(species_cache)
    before = Config.load(wd).stages["genome"].completed
    shutil.rmtree(wd / "genomes")
    run_repgenr("genome", "-wd", wd, "--accession-list-only")
    listed = (wd / "ncbi_acc_download_list.txt").read_text(encoding="utf-8").split()
    assert len(listed) == 10 and not list((wd / "genomes").glob("*.fasta"))
    assert Config.load(wd).stages["genome"].completed == before, "query mode leaves the record"


def test_genome_keep_files_retains_the_download_scratch(
    run_repgenr, species_cache, copy_of
) -> None:
    wd = copy_of(species_cache)
    shutil.rmtree(wd / "genomes")
    run_repgenr("--force", "genome", "-wd", wd, "--keep-files")
    assert len(list_fasta(wd / "genomes")) == 10
    scratch = wd / "scratch" / "genome_download"
    assert scratch.is_dir() and any(scratch.iterdir()), "--keep-files keeps the datasets archives"


# --- NCBI Virus -----------------------------------------------------------------


def test_vmetadata_ncbi_virus_complete_only(viral_cache: Path) -> None:
    import json

    records = json.loads(
        (viral_cache / "virus_download_wd" / "virus_records.json").read_text(encoding="utf-8")
    )
    assert len(records) > 300 and {r["completeness"] for r in records} == {"COMPLETE"}
    assert (viral_cache / "virus_metadata_base.tsv").is_file()
    assert Config.load(viral_cache).stages["vmetadata"].tool_versions.get("datasets")


def test_mammarenavirus_species_come_from_the_lineage(
    run_repgenr, mammarenavirus_cache: Path, copy_of
) -> None:
    """Junin (taxid 2169991) carries two organism names; strain-level names
    such as 'Sabia virus' sit below the binomial. Each is one species."""
    import json

    from repgenr.core.contracts import SEGMENTS_TSV

    wd = copy_of(mammarenavirus_cache)
    records = json.loads((wd / "virus_download_wd" / "virus_records.json").read_text())
    organisms = {r["organism"] for r in records}
    species = {r["species"] for r in records}
    assert len(species) < len(organisms)
    junin = {r["species"] for r in records if r["taxid"] == "2169991"}
    assert junin == {"Mammarenavirus-juninense"}
    binomial = [
        r for r in records if any(re.fullmatch(r"Mammarenavirus [a-z-]+", n) for n in r["lineage"])
    ]
    assert binomial and all(re.fullmatch(r"Mammarenavirus-[a-z-]+", r["species"]) for r in binomial)

    run_repgenr(
        "vgenome",
        "-wd",
        wd,
        "-ts",
        "Argentinian mammarenavirus",
        "--group-segments",
        "--no-outgroup",
    )
    ingroup, _ = _rows(wd)
    assert ingroup and {r.species for r in ingroup} == {"Mammarenavirus-juninense"}
    assert all(
        n.startswith("Arenaviridae_Mammarenavirus_Mammarenavirus-juninense_")
        for n in (r.filename for r in ingroup)
    )
    rows = [line.split("\t") for line in (wd / SEGMENTS_TSV).read_text().splitlines()[1:]]
    assert rows and {row[2] for row in rows} <= {"S", "L"}


def test_vmetadata_released_after_and_host_narrow_the_set(
    run_repgenr, viral_cache, tmp_path
) -> None:
    import json

    def count(wd: Path) -> int:
        return len(json.loads((wd / "virus_download_wd" / "virus_records.json").read_text()))

    base = count(viral_cache)
    recent = tmp_path / "recent"
    run_repgenr(
        "vmetadata",
        "-wd",
        recent,
        "--target",
        "hepatovirus",
        "--complete-only",
        "--released-after",
        "01/01/2022",
    )
    assert 0 < count(recent) < base
    hosted = tmp_path / "host"
    run_repgenr(
        "vmetadata",
        "-wd",
        hosted,
        "--target",
        "hepatovirus",
        "--complete-only",
        "--host",
        "Homo sapiens",
    )
    assert 0 < count(hosted) <= base
    assert "--host Homo sapiens" in log_text(hosted)


def test_vmetadata_list_targets_reaches_bvbrc(run_repgenr, tmp_path: Path) -> None:
    out = run_repgenr("vmetadata", "-wd", tmp_path / "wd", "--source", "bvbrc", "--list")
    text = out.stdout + out.stderr
    assert "Adenoviridae" in text and "Picornaviridae" in text


@pytest.mark.requires_binary("mashtree")
def test_maguari_species_comes_from_ncbi_taxonomy(cacheense_cache: Path) -> None:
    """The report lineage nests Maguari virus (taxid 11575) under
    'Orthobunyavirus cacheense'; NCBI Taxonomy has it under
    'Orthobunyavirus maguariense', and the record takes that species."""
    import json

    records = json.loads((cacheense_cache / "virus_download_wd" / "virus_records.json").read_text())
    maguari = [r for r in records if r["taxid"] == "11575"]
    assert maguari and {r["species"] for r in maguari} == {"Orthobunyavirus-maguariense"}
    assert {r["species_source"] for r in maguari} == {"taxonomy"}
    assert {r["species"] for r in records if r["taxid"] == "35304"} <= {"Orthobunyavirus-cacheense"}
    params = Config.load(cacheense_cache).stages["vmetadata"].params
    assert params["species_source"].get("taxonomy") == len(records)


def test_hepatovirus_species_is_the_binomial(viral_cache: Path) -> None:
    """NCBI keeps 'Hepatovirus A' below 'Hepatovirus ahepa'; the binomial is
    the species, and all records of taxid 12092 share it."""
    import json

    records = json.loads((viral_cache / "virus_download_wd" / "virus_records.json").read_text())
    assert {r["species"] for r in records if r["taxid"] == "12092"} == {"Hepatovirus-ahepa"}
    assert "Hepatovirus-A" not in {r["species"] for r in records}


def test_vgenome_selection_flags(run_repgenr, viral_cache: Path, copy_of) -> None:
    wd = copy_of(viral_cache)

    def select(*flags: str) -> tuple[int, int]:
        run_repgenr("vgenome", "-wd", wd, *flags)
        ingroup, outgroup = _rows(wd)
        return len(ingroup), len(outgroup)

    base, _ = select("-tg", "Hepatovirus")
    assert base > 300
    narrow, _ = select("-tg", "Hepatovirus", "--length-range", "7400-7500")
    assert 0 < narrow < base
    tight, _ = select("-tg", "Hepatovirus", "--length-deviation", "1")
    assert 0 < tight < base
    loose, _ = select("-tg", "Hepatovirus", "--length-all")
    assert loose >= base
    mean, _ = select("-tg", "Hepatovirus", "--length-method", "mean")
    assert mean > 0
    # Species-level target leaves the other hepatovirus species as outgroup
    # candidates; mashtree picks one.
    n_species, og = select("-ts", "Hepatovirus A")
    assert 0 < n_species < loose and og == 1
    assert Config.load(wd).stages["vgenome"].tool_versions.get("mashtree")
    _, og_min = select("-ts", "Hepatovirus A", "--outgroup-candidates-taxid-min-genomes", "1000")
    assert og_min == 0 and "No outgroup candidates found" in log_text(wd)
    _, og_none = select("-ts", "Hepatovirus A", "--no-outgroup")
    assert og_none == 0
    n_sero, _ = select("-tse", "hepatitis A")
    assert n_sero > 0
    n_custom, _ = select("-tc", "completeness:COMPLETE", "-tg", "Hepatovirus")
    assert n_custom == base
    n_group, _ = select("-tg", "Hepatovirus", "--group-segments")
    # No real segment labels: nothing is grouped. The grouping path skips the
    # per-segment length window, so the count can exceed the filtered base.
    assert n_group >= base and "Grouped 0 segment sequences into 0 isolate genomes" in log_text(wd)


@pytest.mark.requires_binary("mashtree")
def test_vgenome_discard_glance_headers_keep_files(run_repgenr, viral_cache: Path, copy_of) -> None:
    wd = copy_of(viral_cache)
    base = len(_rows(wd)[0])
    headers = run_repgenr(
        "vgenome", "-wd", wd, "-tg", "Hepatovirus", "--print-fasta-headers", "--glance"
    )
    text = headers.stdout + headers.stderr
    assert "stopping before writing genomes" in text and "hepatovirus" in text.lower()
    # --glance is a query: nothing rewritten, the record untouched.
    assert len(_rows(wd)[0]) == base

    first_header = next(
        line[1:].split()[0]
        for line in (wd / "virus_download_wd" / "download.fa").open()
        if line.startswith(">")
    )
    run_repgenr("vgenome", "-wd", wd, "-tg", "Hepatovirus", "--discard", first_header)
    assert len(_rows(wd)[0]) == base - 1, "--discard drops the record whose header carries the tag"

    run_repgenr(
        "vgenome",
        "-wd",
        wd,
        "-ts",
        "Hepatovirus A",
        "--keep-files",
        "--outgroup-treebuilder",
        "mashtree",
    )
    assert (wd / "virus_outgroup_wd").is_dir(), "--keep-files keeps the outgroup scratch"


# --- BV-BRC ---------------------------------------------------------------------


@pytest.fixture(scope="session")
def bvbrc_cache(cached_workdir) -> Path:
    return cached_workdir(
        "hepatitis_e_bvbrc",
        [["vmetadata", "--target", "hepatitis_e_virus", "--source", "bvbrc"]],
        "vmetadata",
    )


def test_vmetadata_bvbrc_source_and_filter(bvbrc_cache: Path, run_repgenr, copy_of) -> None:
    assert (bvbrc_cache / "virus_download_wd" / "download.fa").is_file()
    assert (bvbrc_cache / "virus_metadata_base.tsv").is_file()
    wd = copy_of(bvbrc_cache)
    run_repgenr(
        "--force",
        "vmetadata",
        "-wd",
        wd,
        "--target",
        "hepatitis_e_virus",
        "--source",
        "bvbrc",
        "--filter",
        "segment",
    )
    assert "Group FASTA already present" in log_text(wd)
    assert Config.load(wd).stages["vmetadata"].params.get("filter") == "segment"


@pytest.mark.requires_binary("mashtree")
def test_vgenome_bvbrc_needs_ignore_duplicates(bvbrc_cache: Path, run_repgenr, copy_of) -> None:
    wd = copy_of(bvbrc_cache)
    refused = run_repgenr("vgenome", "-wd", wd, "-tg", "Paslahepevirus", check=False)
    assert refused.returncode != 0
    assert "Duplicate sequence id" in refused.stderr + refused.stdout
    run_repgenr(
        "vgenome", "-wd", wd, "-tg", "Paslahepevirus", "--ignore-duplicates", "--no-outgroup"
    )
    ingroup, _ = _rows(wd)
    assert len(ingroup) > 0 and len(list_fasta(wd / "genomes")) == len(ingroup)


# --- run -----------------------------------------------------------------------


def test_run_dry_run_prints_the_chain_without_network(run_repgenr, tmp_path: Path) -> None:
    wd = tmp_path / "dry"
    out = run_repgenr(
        "run",
        "-wd",
        wd,
        "-d",
        "rep",
        "-l",
        "genus",
        "-tg",
        "Francisella",
        "--dry-run",
        "--tool",
        "sourmash",
        "--treebuilder",
        "mashtree",
    ).stdout
    for stage in ("metadata", "genome", "dereplicate", "phylo", "tree2tax"):
        assert f"- {stage}" in out
    assert "no work done" in out and not wd.exists()
    viral = run_repgenr(
        "run", "-wd", wd, "--viral", "--target", "hepatovirus", "-tg", "Hepatovirus", "--dry-run"
    ).stdout
    assert "- vmetadata" in viral and "- vgenome" in viral


@pytest.mark.requires_binary("sourmash", "mashtree")
def test_run_bacterial_chain_end_to_end(run_repgenr, tmp_path: Path) -> None:
    wd = tmp_path / "run"
    run_repgenr(
        "run",
        "-wd",
        wd,
        "-d",
        "rep",
        "-l",
        "genus",
        "-tg",
        "Francisella",
        "--metadata-source",
        "api",
        "--outgroup-accession",
        GENUS_OUTGROUP,
        "--tool",
        "sourmash",
        "--treebuilder",
        "mashtree",
        "--threads",
        "4",
        "--keeper",
        "tool",
        "--secondary-ani",
        "0.98",
        "--primary-ani",
        "0.85",
        "--aligned-fraction",
        "0.3",
        "--limit",
        "12",
        "--process-size",
        "6",
        "--collapse-support",
        "0.5",
        "--all-genomes",
    )
    assert (wd / TREE2TAX_TSV).is_file()
    status = run_repgenr("status", "-wd", wd).stdout
    assert "All stages complete" in status
    stages = Config.load(wd).stages
    assert stages["dereplicate"].params["secondary_ani"] == 0.98
    assert stages["dereplicate"].params["process_size"] == 6
    assert stages["phylo"].params["all_genomes"] is True
    assert stages["tree2tax"].params["collapse_support"] == 0.5
    assert len(list_fasta(wd / "genomes")) <= 12, "--limit reached the metadata stage"
    assert stages["dereplicate"].params["keeper"] == "tool"
    assert (wd / "outgroup_accession.txt").read_text(encoding="utf-8").strip() == GENUS_OUTGROUP


@pytest.mark.requires_binary("sourmash", "mashtree")
def test_run_viral_chain_end_to_end(run_repgenr, tmp_path: Path) -> None:
    wd = tmp_path / "vrun"
    # run --viral has no --complete-only (design item D-3): the length window
    # is what keeps this to a few hundred records. Without it (--group-segments
    # skips the window) every hepatovirus record is kept, ~9500 mostly partial,
    # and mashtree's tree step fails on that many taxa.
    run_repgenr(
        "run",
        "-wd",
        wd,
        "--viral",
        "--target",
        "hepatovirus",
        "-tg",
        "Hepatovirus",
        "--viral-source",
        "ncbi_virus",
        "--tool",
        "sourmash",
        "--treebuilder",
        "mashtree",
        "--threads",
        "4",
        "--no-outgroup",
        "--no-include-dereplicated",
    )
    assert (wd / TREE2TAX_TSV).is_file()
    assert Config.load(wd).stages["vgenome"].params["no_outgroup"] is True
    assert Config.load(wd).stages["tree2tax"].params["include_dereplicated"] is False
    assert "Pipeline: viral" in run_repgenr("status", "-wd", wd).stdout


def test_genome_fetch_step(run_repgenr, species_cache: Path, tmp_path: Path) -> None:
    """The stateless download step on a cached selection.tsv."""
    out = tmp_path / "fetched"
    run_repgenr(
        "genome-fetch",
        "--selection",
        species_cache / SELECTION_TSV,
        "-o",
        out,
        "--keep-files",
        "--versions-out",
        tmp_path / "versions.yml",
    )
    assert len(list_fasta(out / "genomes")) == 10
    assert len(list_fasta(out / "outgroup")) == 1
    assert "datasets:" in (tmp_path / "versions.yml").read_text(encoding="utf-8")
    assert any((out / "scratch").iterdir()), "--keep-files leaves the download scratch"


def test_run_dry_run_reports_family_and_species_targets(run_repgenr, tmp_path: Path) -> None:
    out = run_repgenr(
        "run",
        "-wd",
        tmp_path / "dry",
        "-d",
        "rep",
        "-l",
        "species",
        "-tf",
        "Francisellaceae",
        "-tg",
        "Francisella",
        "-ts",
        "tularensis",
        "--dry-run",
        "--aligner",
        "sibeliaz",
        "--treebuilder",
        "iqtree",
    ).stdout
    assert "family=Francisellaceae" in out and "species=tularensis" in out
    assert "aligner=sibeliaz" in out, "the dry run names the aligner the chain would use"


def test_reads_selects_runs_from_ena(run_repgenr, tmp_path: Path) -> None:
    """ENA taxon query, synonym resolution and NCBI labelling, as a real process."""
    from repgenr.core.contracts import READS_TSV, read_reads

    wd = tmp_path / "reads"
    out = run_repgenr(
        "reads",
        "-wd",
        wd,
        "-ts",
        "Mycoplasma genitalium",
        "--platform",
        "illumina",
        "--max-runs",
        "3",
    )
    assert "taxid 2097" in out.stderr + out.stdout
    rows = read_reads(wd / READS_TSV)
    assert len(rows) == 3
    assert all(r.platform == "ILLUMINA" and r.fastq_urls for r in rows)
    assert all(
        (r.family, r.genus, r.species) == ("Mycoplasmoidaceae", "Mycoplasmoides", "genitalium")
        for r in rows
    )
    assert [r.bases for r in rows] == sorted((r.bases for r in rows), reverse=True)
