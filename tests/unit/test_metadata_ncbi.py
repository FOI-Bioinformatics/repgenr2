"""The NCBI selection source (metadata --source ncbi), offline.

``datasets summary`` is faked: the fake answers the taxonomy and genome
queries from a small table, in the JSON-lines layout of datasets 18, and
records every argument list it was given.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import pytest

from repgenr.core import ncbi_datasets
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SELECTION_TSV
from repgenr.core.errors import ToolExecutionError, UserInputError, WorkdirError
from repgenr.sources import gtdb_api, ncbi
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams

LOG = logging.getLogger("test.metadata_ncbi")

_HIGHER = {
    "domain": (2, "Bacteria"),
    "phylum": (1224, "Pseudomonadota"),
    "class": (1236, "Gammaproteobacteria"),
    "order": (72273, "Thiotrichales"),
}
_FAMILY = (34064, "Francisellaceae")
_GENUS = (262, "Francisella")
_ALLO = (1869285, "Allofrancisella")

# taxid -> (rank, name, classification above and including the node, parents)
_TAXA: dict[int, tuple[str, str, dict[str, tuple[int, str]], list[int]]] = {
    34064: ("FAMILY", "Francisellaceae", {**_HIGHER, "family": _FAMILY}, [1, 2, 1224, 1236, 72273]),
    262: (
        "GENUS",
        "Francisella",
        {**_HIGHER, "family": _FAMILY, "genus": _GENUS},
        [1, 2, 1224, 1236, 72273, 34064],
    ),
    263: (
        "SPECIES",
        "Francisella tularensis",
        {**_HIGHER, "family": _FAMILY, "genus": _GENUS, "species": (263, "Francisella tularensis")},
        [1, 2, 1224, 1236, 72273, 34064, 262],
    ),
    119856: (
        "SUBSPECIES",
        "Francisella tularensis subsp. tularensis",
        {**_HIGHER, "family": _FAMILY, "genus": _GENUS, "species": (263, "Francisella tularensis")},
        [1, 2, 1224, 1236, 72273, 34064, 262, 263],
    ),
    28110: (
        "SPECIES",
        "Francisella philomiragia",
        {
            **_HIGHER,
            "family": _FAMILY,
            "genus": _GENUS,
            "species": (28110, "Francisella philomiragia"),
        },
        [1, 2, 1224, 1236, 72273, 34064, 262],
    ),
    573569: (
        "SPECIES",
        "Francisella salina",
        {**_HIGHER, "family": _FAMILY, "genus": _GENUS, "species": (573569, "Francisella salina")},
        [1, 2, 1224, 1236, 72273, 34064, 262],
    ),
    1085644: (
        "SPECIES",
        "Allofrancisella frigidaquae",
        {
            **_HIGHER,
            "family": _FAMILY,
            "genus": _ALLO,
            "species": (1085644, "Allofrancisella frigidaquae"),
        },
        [1, 2, 1224, 1236, 72273, 34064, 1869285],
    ),
    # Bacillus: a bacterial genus and a genus of stick insects share the name.
    1386: (
        "GENUS",
        "Bacillus",
        {"domain": (2, "Bacteria"), "family": (186817, "Bacillaceae"), "genus": (1386, "Bacillus")},
        [1, 2, 186817],
    ),
    55087: (
        "GENUS",
        "Bacillus",
        {"domain": (2759, "Eukaryota"), "genus": (55087, "Bacillus")},
        [1, 2759],
    ),
}
_AMBIGUOUS = {"Bacillus": (55087, 1386)}


def _genome(
    accession: str,
    tax_id: int,
    organism: str,
    *,
    paired: str = "",
    checkm: tuple[float, float] | None = None,
    reference: bool = False,
    n50: int = 100_000,
    strain: str = "",
) -> dict:
    row: dict = {
        "accession": accession,
        "assembly_info": {"assembly_level": "Complete Genome"},
        "assembly_stats": {"contig_n50": n50, "number_of_contigs": 3},
        "organism": {
            "organism_name": organism,
            "tax_id": tax_id,
            "infraspecific_names": {"strain": strain} if strain else {},
        },
    }
    if paired:
        row["paired_accession"] = paired
    if checkm:
        row["checkm_info"] = {"completeness": checkm[0], "contamination": checkm[1]}
    if reference:
        row["assembly_info"]["refseq_category"] = "reference genome"
        row["type_material"] = {"type_display_text": "assembly from type material"}
    return row


_GENOMES = [
    _genome(
        "GCF_000001.1",
        119856,
        "Francisella tularensis subsp. tularensis SCHU S4",
        paired="GCA_000001.1",
        checkm=(99.0, 0.1),
        strain="SCHU S4",
    ),
    _genome("GCA_000001.1", 119856, "Francisella tularensis subsp. tularensis SCHU S4",
            paired="GCF_000001.1", checkm=(99.0, 0.1)),
    _genome("GCF_000002.1", 263, "Francisella tularensis", paired="GCA_000002.1",
            checkm=(98.0, 0.5)),
    _genome("GCA_000002.1", 263, "Francisella tularensis", paired="GCF_000002.1",
            checkm=(98.0, 0.5)),
    # GenBank only, and without CheckM values.
    _genome("GCA_000003.1", 263, "Francisella tularensis"),
    _genome("GCF_000004.1", 28110, "Francisella philomiragia", paired="GCA_000004.1",
            checkm=(97.0, 1.0), reference=True),
    _genome("GCF_000005.1", 573569, "Francisella salina", paired="GCA_000005.1",
            checkm=(99.5, 0.1), reference=True, n50=2_000_000),
    _genome("GCF_000009.1", 1085644, "Allofrancisella frigidaquae", paired="GCA_000009.1",
            checkm=(90.0, 0.5), reference=True),
]  # fmt: skip


def _taxonomy_row(query: str, tax_id: int) -> dict:
    rank, name, classification, parents = _TAXA[tax_id]
    return {
        "query": [query],
        "taxonomy": {
            "tax_id": tax_id,
            "rank": rank,
            "current_scientific_name": {"name": name},
            "classification": {r: {"id": i, "name": n} for r, (i, n) in classification.items()},
            "parents": parents,
        },
    }


@pytest.fixture
def datasets(monkeypatch) -> list[list[str]]:
    """Fake ``datasets summary``; returns the argument lists it received."""
    calls: list[list[str]] = []
    by_name = {name: tid for tid, (_rank, name, _c, _p) in _TAXA.items() if tid != 55087}

    def summary(args, report: Path, logger) -> None:
        args = [str(a) for a in args]
        calls.append(args)
        rows: list[dict] = []
        if args[:2] == ["taxonomy", "taxon"]:
            queries = (
                Path(args[args.index("--inputfile") + 1]).read_text(encoding="utf-8").split("\n")
            )
            for query in (q for q in queries if q):
                if query in _AMBIGUOUS:
                    listed = "\n".join(f"{query} (genus, taxid: {t}, x)" for t in _AMBIGUOUS[query])
                    raise ToolExecutionError(
                        ["datasets"],
                        1,
                        output=f"Error: The taxonomy name '{query}' is an exact match for "
                        f"more than one taxid. Please select one of these taxids:\n{listed}",
                    )
                tax_id = int(query) if query.isdigit() else by_name.get(query)
                if tax_id is None or tax_id not in _TAXA:
                    if query.isdigit():
                        continue
                    raise ToolExecutionError(
                        ["datasets"],
                        1,
                        output=f"Error: The taxonomy name '{query}' is not recognized.",
                    )
                rows.append(_taxonomy_row(query, tax_id))
        elif args[:2] == ["genome", "taxon"]:
            tax_id = int(args[2])
            source = args[args.index("--assembly-source") + 1]
            prefix = {"RefSeq": ("GCF_",), "GenBank": ("GCA_",), "all": ("GCF_", "GCA_")}[source]
            for row in _GENOMES:
                own = row["organism"]["tax_id"]
                if tax_id != own and tax_id not in _TAXA[own][3]:
                    continue
                if not row["accession"].startswith(prefix):
                    continue
                if "--reference" in args and "refseq_category" not in row["assembly_info"]:
                    continue
                rows.append(row)
        elif args[:2] == ["genome", "accession"]:
            rows = [row for row in _GENOMES if row["accession"] == args[2]]
        report.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")

    monkeypatch.setattr(ncbi_datasets, "summary", summary)
    monkeypatch.setattr(ncbi, "preflight", lambda _caps: {"datasets": "18.14.0"})
    monkeypatch.setattr(ncbi, "require_reachable", lambda *_a, **_k: None)
    return calls


def _params(**overrides) -> MetadataParams:
    defaults = dict(
        dataset="all",
        level="species",
        source="ncbi",
        target_genus="Francisella",
        target_species="tularensis",
    )
    defaults.update(overrides)
    return MetadataParams(**defaults)


def _selection(workdir: Path) -> dict[str, dict]:
    with open(workdir / SELECTION_TSV, encoding="utf-8", newline="") as fo:
        return {row["accession"]: row for row in csv.DictReader(fo, delimiter="\t")}


def _genome_queries(calls: list[list[str]]) -> list[list[str]]:
    return [c for c in calls if c[:2] == ["genome", "taxon"]]


def test_species_selection_is_refseq_by_default(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    # What a GTDB API run left for census does not describe an NCBI selection.
    (ctx.workdir / gtdb_api.GTDB_API_GENOMES).write_text("accession\n", encoding="utf-8")

    assert metadata.run(ctx, _params()) == 2
    assert not (ctx.workdir / gtdb_api.GTDB_API_GENOMES).exists()

    rows = _selection(ctx.workdir)
    ingroup = {acc for acc, row in rows.items() if row["is_outgroup"] == "0"}
    assert ingroup == {"GCF_000001.1", "GCF_000002.1"}
    # A subspecies taxid is named by its species, as one token.
    assert rows["GCF_000001.1"]["species"] == "tularensis"
    assert rows["GCF_000001.1"]["filename"] == (
        "Francisellaceae_Francisella_tularensis_GCF_000001.1.fasta"
    )
    assert rows["GCF_000001.1"]["completeness"] == "99.00"
    assert all(row["gtdb_representative"] == "0" for row in rows.values())
    # The outgroup: the best reference genome of the genus outside the species.
    outgroup = [acc for acc, row in rows.items() if row["is_outgroup"] == "1"]
    assert outgroup == ["GCF_000005.1"]
    assert (ctx.workdir / "outgroup_accession.txt").read_text() == "GCF_000005.1\n"

    first = _genome_queries(datasets)[0]
    assert first[2] == "263"  # queried by taxid, not by name
    assert first[first.index("--assembly-source") + 1] == "RefSeq"
    assert "--exclude-atypical" in first and "--reference" not in first

    sources = {g.accession: g.source for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert set(sources.values()) == {"ncbi"}


def test_record_names_the_query_date_the_filters_and_the_client(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(release="232.0"))
    record = ctx.config.stages["metadata"]
    assert record.tool == "ncbi-datasets"
    assert record.tool_versions == {"datasets": "18.14.0"}
    assert record.params["assembly_source"] == "refseq"
    assert record.params["query_date"].startswith("20")
    assert record.params["outgroup_rule"] == "quality"
    # --release is ignored by this source; no GTDB release is reported for it.
    assert metadata.gtdb_provenance(record.params) == {
        "ncbi_query_date": record.params["query_date"]
    }
    stats = json.loads((ctx.workdir / "selection_stats.json").read_text())
    assert stats["values"]["source"] == "ncbi"
    assert stats["values"]["candidates"] == 2


def test_candidates_and_lineage_are_kept_for_census_and_the_report(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params())
    kept = ncbi.read_genomes(ctx.workdir / ncbi.NCBI_GENOMES)
    assert [g["accession"] for g in kept] == ["GCF_000001.1", "GCF_000002.1"]
    assert kept[0]["tax"]["species_name"] == "Francisella tularensis"

    with open(ctx.workdir / "metadata_selected.tsv", encoding="utf-8", newline="") as fo:
        details = {row["accession"]: row for row in csv.DictReader(fo, delimiter="\t")}
    row = details["GCF_000001.1"]
    assert row["accession_ncbi"] == "GCA_000001.1"
    assert row["gtdb_taxonomy"] == ""
    assert row["ncbi_taxonomy"].endswith("g__Francisella;s__Francisella tularensis")
    assert row["ncbi_taxonomy_unfiltered"].endswith(
        "x__Francisella tularensis subsp. tularensis SCHU S4"
    )
    assert row["ncbi_strain_identifiers"] == "SCHU S4"


def test_all_lists_an_assembly_once_as_refseq_where_both_exist(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(assembly_source="all", no_outgroup=True))
    assert set(_selection(ctx.workdir)) == {"GCF_000001.1", "GCF_000002.1", "GCA_000003.1"}
    assert not (ctx.workdir / "outgroup_accession.txt").exists()


def test_genbank_selection_holds_genbank_accessions_only(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(assembly_source="genbank"))
    rows = _selection(ctx.workdir)
    assert set(rows) == {"GCA_000001.1", "GCA_000002.1", "GCA_000003.1", "GCA_000005.1"}
    # Reference genomes are listed in RefSeq; the GenBank copy is the outgroup.
    assert rows["GCA_000005.1"]["is_outgroup"] == "1"
    # No CheckM values at NCBI: the columns stay empty.
    assert rows["GCA_000003.1"]["completeness"] == ""


def test_rep_dataset_asks_for_reference_genomes(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(dataset="rep", level="genus", target_species=None))
    rows = _selection(ctx.workdir)
    ingroup = {acc for acc, row in rows.items() if row["is_outgroup"] == "0"}
    assert ingroup == {"GCF_000004.1", "GCF_000005.1"}
    assert "--reference" in _genome_queries(datasets)[0]
    # A genus target takes its outgroup from the family.
    assert rows["GCF_000009.1"]["is_outgroup"] == "1"
    assert rows["GCF_000009.1"]["genus"] == "Allofrancisella"


def test_several_species_are_queried_one_by_one(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(target_species="tularensis,philomiragia"))
    rows = _selection(ctx.workdir)
    ingroup = {acc for acc, row in rows.items() if row["is_outgroup"] == "0"}
    assert ingroup == {"GCF_000001.1", "GCF_000002.1", "GCF_000004.1"}
    assert [q[2] for q in _genome_queries(datasets)[:2]] == ["263", "28110"]
    # philomiragia is a target now, so it cannot be the outgroup.
    assert rows["GCF_000005.1"]["is_outgroup"] == "1"


def test_limit_keeps_the_best_quality(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(limit=1, no_outgroup=True)) == 1
    assert set(_selection(ctx.workdir)) == {"GCF_000001.1"}
    stats = json.loads((ctx.workdir / "selection_stats.json").read_text())
    assert stats["values"]["candidates"] == 2


def test_first_mode_takes_the_lowest_accession(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(outgroup_mode="first"))
    assert (ctx.workdir / "outgroup_accession.txt").read_text() == "GCF_000004.1\n"


def test_named_outgroup(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(outgroup_accession="GCF_000009.1"))
    assert (ctx.workdir / "outgroup_accession.txt").read_text() == "GCF_000009.1\n"
    assert ctx.config.stages["metadata"].params["outgroup_rule"] == "named"


@pytest.mark.parametrize(
    ("accession", "message"),
    [
        ("GCF_000001.1", "part of the selection"),
        ("GCA_000003.1", "inside the target species Francisella tularensis"),
        ("GCF_999999.1", "not found at NCBI"),
    ],
)
def test_named_outgroup_is_refused(tmp_path, datasets, accession, message) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match=message):
        metadata.run(ctx, _params(outgroup_accession=accession))


def test_no_reference_genome_outside_the_target_asks_for_an_accession(
    tmp_path, datasets, monkeypatch
) -> None:
    monkeypatch.setattr(ncbi, "_PARENT_RANKS", {"family": ("order",)})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="--outgroup-accession"):
        metadata.run(
            ctx,
            _params(level="family", target_family="Francisellaceae", target_genus=None,
                    target_species=None),
        )  # fmt: skip


def test_a_name_shared_with_a_eukaryote_resolves_to_the_bacterial_genus(tmp_path, datasets) -> None:
    taxon = ncbi_datasets.resolve_name("Bacillus", "genus", tmp_path, LOG)
    assert taxon.tax_id == "1386"
    assert taxon.domain == "Bacteria"


@pytest.mark.parametrize(
    ("name", "rank", "message"),
    [
        ("Nonexistentus", "genus", "no taxon named 'Nonexistentus'"),
        ("Francisellaceae", "genus", "not a bacterial or archaeal genus"),
        ("Francisella tularensis subsp. tularensis", "species", "not a bacterial or archaeal"),
    ],
)
def test_unusable_target_names_are_refused(tmp_path, datasets, name, rank, message) -> None:
    with pytest.raises(UserInputError, match=message):
        ncbi_datasets.resolve_name(name, rank, tmp_path, LOG)


def test_no_genomes_names_the_filters(tmp_path, datasets) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="returned no genomes for Francisella salina"):
        metadata.run(ctx, _params(target_species="salina", assembly_source="genbank"))


def test_filters_reach_the_client() -> None:
    flags = ncbi.genome_flags(
        _params(
            assembly_level="complete,chromosome",
            include_atypical=True,
            exclude_mag=True,
            released_before="2024-01-01",
        ),
        reference=True,
    )
    assert flags == [
        "--assembly-source", "RefSeq", "--assembly-version", "latest", "--reference",
        "--mag", "exclude", "--assembly-level", "complete,chromosome",
        "--released-before", "2024-01-01",
    ]  # fmt: skip


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"assembly_source": "ensembl"}, "--assembly-source must be one of"),
        ({"assembly_level": "complete,draft"}, "--assembly-level takes"),
        ({"released_before": "2010"}, "--released-before takes a date"),
    ],
)
def test_option_values_are_checked(overrides, message) -> None:
    with pytest.raises(UserInputError, match=message):
        metadata._validate(_params(**overrides))


@pytest.mark.parametrize("source", ["tsv", "api"])
def test_ncbi_options_are_refused_by_the_gtdb_sources(source) -> None:
    params = _params(source=source, release="232.0", version="bac120", exclude_mag=True)
    with pytest.raises(UserInputError, match="--exclude-mag cannot be used with --source"):
        metadata._validate(params)


@pytest.mark.parametrize(
    ("genus", "species", "token"),
    [
        ("Francisella", "Francisella tularensis", "tularensis"),
        ("Francisella", "uncultured Francisella sp.", "uncultured-sp"),
        ("Francisella", "Candidatus Francisella endociliophora", "Candidatus-endociliophora"),
        ("Francisella", "Francisella sp. 19X1-34", "sp-19X1-34"),
        ("[Clostridium]", "[Clostridium] scindens", "scindens"),
        ("Francisella", "", "unclassified"),
    ],
)
def test_species_tokens_are_safe_in_filenames(genus, species, token) -> None:
    assert ncbi.species_token(genus, species) == token


def test_a_genome_without_a_classification_is_named_unclassified() -> None:
    assert ncbi.tax_of(None)["genus"] == "unclassified"
    assert ncbi.lineage_string(None) == ""
    assert ncbi._token("[Clostridium]") == "Clostridium"


def test_report_lines_that_are_not_json_objects_are_skipped(tmp_path) -> None:
    report = tmp_path / "report.jsonl"
    report.write_text('New version of client available\n{"accession": "GCF_1.1"}\n{broken\n[1]\n')
    rows = list(ncbi_datasets.read_json_lines(report))
    assert rows == [{"accession": "GCF_1.1"}]
    assert ncbi.genome_from_report({}) is None
    genome = ncbi.genome_from_report(rows[0])
    assert genome is not None and genome.completeness is None and genome.n50 is None


def test_methods_text_names_ncbi_and_not_gtdb() -> None:
    from repgenr.reporting.methods import _metadata
    from repgenr.reporting.model import MethodsInfo

    info = MethodsInfo(
        stage="metadata",
        tool="ncbi-datasets",
        tool_versions={"datasets": "18.14.0"},
        params={
            "source": "ncbi",
            "assembly_source": "refseq",
            "query_date": "2026-11-13T19:27:43+00:00",
            "dataset": "rep",
            "level": "genus",
            "target_genus": "Francisella",
            "outgroup": "GCF_000009.1",
            "outgroup_rule": "quality",
        },
    )
    text, cited = _metadata(info)
    assert text.startswith(
        "Genomes were selected from the RefSeq assemblies at NCBI, by NCBI taxonomy, "
        "with datasets 18.14.0 (queried on 2026-11-13) for Francisella (genus)."
    )
    assert "Only NCBI reference genomes were considered." in text
    assert "the NCBI reference genome of best assembly quality (CheckM completeness" in text
    assert "GTDB" not in text and cited == ["datasets"]
