"""End-to-end tests of the metadata stage over a fixture GTDB TSV slice (WP2-4).

The full GTDB table download is replaced by a small gzipped TSV; the API source
runs against canned ``_api_get`` responses. Both paths run through the public
``run()`` and are asserted on the published contract: selection.tsv, the
outgroup file, the manifest, and the recorded config stage.
"""

from __future__ import annotations

import csv
import gzip
import logging
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

from repgenr.core import http as core_http
from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SELECTION_TSV
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.sources import base as sources_base
from repgenr.sources import gtdb_api, gtdb_table
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams

_TAX_PREFIX = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales"


def _tax(family: str, genus: str, species: str) -> str:
    return f"{_TAX_PREFIX};f__{family};g__{genus};s__{genus} {species}"


# accession, representative accession, taxonomy, ncbi accession
_ROWS = [
    ("GCF_000001.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    ("GCF_000002.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    ("GCF_000003.1", "GCF_000003.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    # same genus, different species: the natural outgroup for a species run
    ("GCF_000010.1", "GCF_000010.1", _tax("Francisellaceae", "Francisella", "philomiragia")),
    # different family entirely
    ("GCF_000020.1", "GCF_000020.1", _tax("Piscirickettsiaceae", "Piscirickettsia", "salmonis")),
]


@pytest.fixture()
def gtdb_tsv(tmp_path) -> Path:
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    header = "accession\tgtdb_genome_representative\tgtdb_taxonomy\tncbi_genbank_assembly_accession"
    lines = [header]
    for acc, rep, tax in _ROWS:
        lines.append(f"RS_{acc}\tRS_{rep}\t{tax}\t{acc}")
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    return path


def _params(gtdb_tsv: Path, **overrides) -> MetadataParams:
    defaults = dict(
        dataset="all",
        level="species",
        release="232.0",
        version="bac120",
        target_genus="Francisella",
        target_species="tularensis",
        metadata_path=str(gtdb_tsv),
    )
    defaults.update(overrides)
    return MetadataParams(**defaults)


def _read_selection(workdir: Path) -> list[dict]:
    with open(workdir / SELECTION_TSV, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def test_tsv_species_selection_end_to_end(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    # An API answer left by an earlier --source api run does not outlive a table run.
    (ctx.workdir / gtdb_api.GTDB_API_GENOMES).write_text("accession\n", encoding="utf-8")
    count = metadata.run(ctx, _params(gtdb_tsv))
    assert count == 3
    assert not (ctx.workdir / gtdb_api.GTDB_API_GENOMES).exists()

    rows = _read_selection(ctx.workdir)
    by_acc = {r["accession"]: r for r in rows}
    assert set(by_acc) == {
        "GCF_000001.1",
        "GCF_000002.1",
        "GCF_000003.1",
        "GCF_000010.1",
    }
    # the same-genus representative is the outgroup
    outgroups = [r for r in rows if r["is_outgroup"] in ("1", "True", "true")]
    assert [r["accession"] for r in outgroups] == ["GCF_000010.1"]
    assert (ctx.workdir / "outgroup_accession.txt").read_text(
        encoding="utf-8"
    ).strip() == "GCF_000010.1"
    assert ctx.config.stages["metadata"].tool == "gtdb-table"
    assert ctx.config.stages["metadata"].params["api_query_date"] is None


def test_tsv_rep_dataset_selects_representatives_only(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = metadata.run(ctx, _params(gtdb_tsv, dataset="rep"))
    # GCF_000002.1 is not its own representative
    assert count == 2
    accs = {r["accession"] for r in _read_selection(ctx.workdir)}
    assert "GCF_000002.1" not in accs


def test_tsv_limit_caps_selection(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(gtdb_tsv, limit=1)) == 1


def test_tsv_unknown_target_raises(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="not found"):
        metadata.run(ctx, _params(gtdb_tsv, target_genus="Yersinia", target_species="pestis"))


def test_tsv_explicit_outgroup_accession(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv, outgroup_accession="GCF_000020.1"))
    assert (ctx.workdir / "outgroup_accession.txt").read_text(
        encoding="utf-8"
    ).strip() == "GCF_000020.1"


def test_tsv_unknown_outgroup_accession_raises(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="Outgroup accession"):
        metadata.run(ctx, _params(gtdb_tsv, outgroup_accession="GCF_999999.9"))


def test_obtain_metadata_downloads_and_verifies(tmp_path, gtdb_tsv, monkeypatch) -> None:
    """Without --metadata-path the stage downloads the modern .tsv.gz layout."""
    calls: dict[str, list] = {"download": [], "verify": []}

    def fake_download(url, dest, *, logger, **kw):
        calls["download"].append(url)
        Path(dest).write_bytes(gtdb_tsv.read_bytes())

    def fake_verify(dest, manifest_url, *, logger, **kw):
        calls["verify"].append(manifest_url)

    monkeypatch.setattr(core_http, "download", fake_download)
    monkeypatch.setattr(core_http, "verify_md5_manifest", fake_verify)

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = metadata.run(ctx, _params(gtdb_tsv, metadata_path=None))
    assert count == 3
    assert calls["download"] and calls["download"][0].endswith("232.0/bac120_metadata_r232.tsv.gz")
    assert calls["verify"], "MD5 manifest verification must run on the download"


def test_obtain_metadata_falls_back_to_tarball_layout(tmp_path, gtdb_tsv, monkeypatch) -> None:
    def fake_download(url, dest, *, logger, **kw):
        if url.endswith(".tsv.gz"):
            raise core_http.HTTPStatusError("404", 404)
        # legacy layout: a .tar.gz containing the .tsv
        import tarfile

        tsv = tmp_path / "member.tsv"
        with gzip.open(gtdb_tsv, "rb") as fo:
            tsv.write_bytes(fo.read())
        with tarfile.open(dest, "w:gz") as tar:
            tar.add(tsv, arcname="bac120_metadata_r232.tsv")

    monkeypatch.setattr(core_http, "download", fake_download)
    monkeypatch.setattr(core_http, "verify_md5_manifest", lambda *a, **k: None)

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(gtdb_tsv, metadata_path=None)) == 3


# --- API source ---------------------------------------------------------------


def _api_row(gid: str, genus: str, species: str, family="f__Francisellaceae", is_rep=True):
    return {
        "gid": gid,
        "gtdbFamily": family,
        "gtdbGenus": f"g__{genus}",
        "gtdbSpecies": f"s__{genus} {species}",
        "gtdbIsRep": is_rep,
    }


def test_api_species_selection_end_to_end(tmp_path, monkeypatch) -> None:
    selection_rows = [
        _api_row("GCF_000001.1", "Francisella", "tularensis"),
        _api_row("GCF_000002.1", "Francisella", "tularensis", is_rep=False),
    ]
    parent_rows = selection_rows + [
        _api_row("GCF_000010.1", "Francisella", "philomiragia"),
    ]

    def fake_api_get(path, params=None):
        if "/taxon/" in path and path.endswith("/genomes-detail"):
            if "tularensis" in path:
                return {"rows": selection_rows}
            return {"rows": parent_rows}
        if path.startswith("/genome/") and path.endswith("/card"):
            return {}  # a card without quality; selection must still succeed
        raise AssertionError(f"unexpected API path {path}")

    monkeypatch.setattr(gtdb_api, "_api_get", fake_api_get)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(
        dataset="all",
        level="species",
        source="api",
        target_genus="Francisella",
        target_species="tularensis",
    )
    count = metadata.run(ctx, params)
    assert count == 2
    rows = _read_selection(ctx.workdir)
    outgroups = [r["accession"] for r in rows if r["is_outgroup"] in ("1", "True", "true")]
    assert outgroups == ["GCF_000010.1"]
    assert ctx.config.stages["metadata"].tool == "gtdb-api"
    # The API exposes no release number; the query date stands in for it.
    record = ctx.config.stages["metadata"].params
    assert record["release"] is None
    queried = datetime.fromisoformat(record["api_query_date"])
    assert queried.utcoffset() == timedelta(0)
    assert abs(datetime.now(UTC) - queried) < timedelta(minutes=5)
    # The API answer for the target is kept for `repgenr census -wd`.
    kept = gtdb_api.read_api_genomes(ctx.workdir / gtdb_api.GTDB_API_GENOMES)
    assert [g["accession"] for g in kept] == ["GCF_000001.1", "GCF_000002.1"]
    assert [g["is_rep"] for g in kept] == [True, False]


def _card(completeness=None, contamination=None, *, checkm2: bool = True) -> dict:
    """The shape of ``/genome/{gid}/card``: quality is nested under metadata_gene."""
    prefix = "checkm2" if checkm2 else "checkm"
    gene = {}
    if completeness is not None:
        gene[f"{prefix}_completeness"] = completeness
        gene[f"{prefix}_contamination"] = contamination
    return {"metadataTaxonomy": {}, "metadata_gene": gene}


def _fake_api_with_cards(selection_rows, parent_rows, cards: dict[str, dict]):
    def fake_api_get(path, params=None):
        if "/taxon/" in path and path.endswith("/genomes-detail"):
            return {"rows": selection_rows if "tularensis" in path else parent_rows}
        if path.startswith("/genome/") and path.endswith("/card"):
            gid = path.removeprefix("/genome/").removesuffix("/card")
            card = cards[gid]
            if isinstance(card, Exception):
                raise card
            return card
        raise AssertionError(f"unexpected API path {path}")

    return fake_api_get


def _quality_by_accession(rows: list[dict]) -> dict[str, tuple[str, str]]:
    return {r["accession"]: (r["completeness"], r["contamination"]) for r in rows}


def test_api_selection_reads_quality_from_genome_cards(tmp_path, monkeypatch) -> None:
    # genomes-detail rows carry no quality; it lives on each genome's card.
    selection_rows = [
        _api_row("GCF_000001.1", "Francisella", "tularensis"),
        _api_row("GCF_000002.1", "Francisella", "tularensis", is_rep=False),
    ]
    parent_rows = selection_rows + [_api_row("GCF_000010.1", "Francisella", "philomiragia")]
    cards = {
        "GCF_000001.1": _card(99.96, 0.18),
        "GCF_000002.1": _card(95.1, 1.2, checkm2=False),  # older card: checkm only
        "GCF_000010.1": _card(98.0, 0.5),
    }
    fake = _fake_api_with_cards(selection_rows, parent_rows, cards)
    monkeypatch.setattr(gtdb_api, "_api_get", fake)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(
        ctx,
        MetadataParams(
            dataset="all",
            level="species",
            source="api",
            target_genus="Francisella",
            target_species="tularensis",
        ),
    )
    quality = _quality_by_accession(_read_selection(ctx.workdir))
    assert quality["GCF_000001.1"] == ("99.96", "0.18")
    assert quality["GCF_000002.1"] == ("95.10", "1.20")
    assert quality["GCF_000010.1"] == ("98.00", "0.50")  # the outgroup too


def test_api_card_failure_leaves_that_genome_unscored(tmp_path, monkeypatch, caplog) -> None:
    selection_rows = [
        _api_row("GCF_000001.1", "Francisella", "tularensis"),
        _api_row("GCF_000002.1", "Francisella", "tularensis", is_rep=False),
    ]
    parent_rows = selection_rows + [_api_row("GCF_000010.1", "Francisella", "philomiragia")]
    cards = {
        "GCF_000001.1": _card(99.96, 0.18),
        "GCF_000002.1": WorkdirError("HTTP 404 for /genome/GCF_000002.1/card"),
        "GCF_000010.1": _card(98.0, 0.5),
    }
    fake = _fake_api_with_cards(selection_rows, parent_rows, cards)
    monkeypatch.setattr(gtdb_api, "_api_get", fake)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.logger.addHandler(caplog.handler)  # the workdir logger does not propagate
    with caplog.at_level(logging.WARNING):
        count = metadata.run(
            ctx,
            MetadataParams(
                dataset="all",
                level="species",
                source="api",
                target_genus="Francisella",
                target_species="tularensis",
            ),
        )
    assert count == 2  # one unreadable card does not fail the selection
    quality = _quality_by_accession(_read_selection(ctx.workdir))
    assert quality["GCF_000001.1"] == ("99.96", "0.18")
    assert quality["GCF_000002.1"] == ("", "")
    warnings = [rec.message for rec in caplog.records if rec.levelname == "WARNING"]
    assert any("GCF_000002.1" in message for message in warnings)


def test_api_no_rows_raises(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(gtdb_api, "_api_get", lambda path, params=None: {"rows": []})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(dataset="all", level="genus", source="api", target_genus="Francisella")
    with pytest.raises(UserInputError, match="no genomes"):
        metadata.run(ctx, params)


# --- CheckM quality columns (checkm2_* preferred, checkm_* fallback) ----------


def _write_quality_tsv(
    path: Path, quality_cols: tuple[str, str], quality_values: dict[str, tuple[str, str]]
) -> None:
    """A copy of the ``gtdb_tsv`` fixture layout with two extra quality columns."""
    header = (
        "accession\tgtdb_genome_representative\tgtdb_taxonomy\t"
        f"ncbi_genbank_assembly_accession\t{quality_cols[0]}\t{quality_cols[1]}"
    )
    lines = [header]
    for acc, rep, tax in _ROWS:
        comp, cont = quality_values.get(acc, ("", ""))
        lines.append(f"RS_{acc}\tRS_{rep}\t{tax}\t{acc}\t{comp}\t{cont}")
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")


def test_parse_metadata_carries_checkm2_quality(tmp_path) -> None:
    path = tmp_path / "quality_checkm2.tsv.gz"
    _write_quality_tsv(
        path,
        ("checkm2_completeness", "checkm2_contamination"),
        {"GCF_000001.1": ("99.0", "0.5")},
    )
    accessions = gtdb_table.parse_metadata(path, _params(path), logging.getLogger("test"))
    assert accessions["GCF_000001.1"]["completeness"] == 99.0
    assert accessions["GCF_000001.1"]["contamination"] == 0.5

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(path))
    rec = next(g for g in ctx.manifest.all_genomes() if g.accession == "GCF_000001.1")
    assert rec.completeness == 99.0
    assert rec.contamination == 0.5


def test_parse_metadata_falls_back_to_checkm_quality(tmp_path) -> None:
    path = tmp_path / "quality_checkm.tsv.gz"
    _write_quality_tsv(
        path,
        ("checkm_completeness", "checkm_contamination"),
        {"GCF_000001.1": ("98.0", "1.2")},
    )
    accessions = gtdb_table.parse_metadata(path, _params(path), logging.getLogger("test"))
    assert accessions["GCF_000001.1"]["completeness"] == 98.0
    assert accessions["GCF_000001.1"]["contamination"] == 1.2


def test_narrower_reselection_shrinks_manifest(tmp_path, gtdb_tsv) -> None:
    """Re-running metadata with a narrower selection removes de-selected rows."""
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv))  # 3 tularensis + outgroup
    assert len(ctx.manifest.all_genomes(include_outgroup=True)) == 4

    metadata.run(ctx, _params(gtdb_tsv, limit=1))
    accs = {g.accession for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert len(accs) == 2  # one selected + one outgroup, older rows gone


# --- --limit: species-stratified, quality-ranked within species ---------------


def _cand(acc: str, species: str, comp=None, cont=None, is_rep=False):
    return sources_base.LimitCandidate(
        accession=acc,
        species=species,
        is_rep=is_rep,
        completeness=comp,
        contamination=cont,
    )


def test_stratified_limit_round_robins_over_species_best_first() -> None:
    cands = [
        _cand("A1", "alpha", 90.0, 1.0),
        _cand("A2", "alpha", 99.0, 0.5),  # best alpha
        _cand("A3", "alpha", 95.0, 1.0),
        _cand("B1", "beta", 80.0, 0.0),  # only beta
        _cand("C1", "gamma", 97.0, 0.2),  # best gamma
        _cand("C2", "gamma", 96.0, 0.2),
    ]
    kept = [c.accession for c in sources_base.stratified_limit(cands, 4)]
    # First pass: best of each species (alphabetical species order); second pass: next best.
    assert kept == ["A2", "B1", "C1", "A3"]


def test_stratified_limit_unscored_genomes_rank_last_within_species() -> None:
    cands = [_cand("X1", "s"), _cand("X2", "s", 70.0, 5.0), _cand("X3", "s", 60.0, 0.0)]
    kept = [c.accession for c in sources_base.stratified_limit(cands, 2)]
    assert kept == ["X3", "X2"]  # 60 - 0 = 60 beats 70 - 25 = 45; unscored X1 last


def test_stratified_limit_ties_prefer_gtdb_representative_then_accession() -> None:
    cands = [
        _cand("Z2", "s", 95.0, 1.0),
        _cand("Z1", "s", 95.0, 1.0),
        _cand("Z9", "s", 95.0, 1.0, is_rep=True),
    ]
    kept = [c.accession for c in sources_base.stratified_limit(cands, 3)]
    assert kept == ["Z9", "Z1", "Z2"]


def test_stratified_limit_without_limit_returns_all_in_ranked_order() -> None:
    cands = [_cand("B1", "beta", 50.0, 0.0), _cand("A1", "alpha", 99.0, 0.0)]
    assert [c.accession for c in sources_base.stratified_limit(cands, None)] == ["A1", "B1"]
    assert len(sources_base.stratified_limit(cands, 10)) == 2


def test_tsv_limit_keeps_best_quality_not_first_rows(tmp_path) -> None:
    # tularensis has three genomes; file order is 000001, 000002, 000003.
    # Scores: 000001 -> 80, 000002 -> 96.5, 000003 -> 90. A limit of 2 must keep
    # the two best, dropping the first row the old slice would have kept.
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    _write_quality_tsv(
        path,
        ("checkm2_completeness", "checkm2_contamination"),
        {
            "GCF_000001.1": ("90.0", "2.0"),
            "GCF_000002.1": ("99.0", "0.5"),
            "GCF_000003.1": ("95.0", "1.0"),
            "GCF_000010.1": ("98.0", "0.1"),
        },
    )
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(path, limit=2)) == 2
    rows = _read_selection(ctx.workdir)
    kept = {r["accession"] for r in rows if r["is_outgroup"] not in ("1", "True", "true")}
    assert kept == {"GCF_000002.1", "GCF_000003.1"}


def test_api_limit_fetches_quality_for_every_candidate_then_stratifies(
    tmp_path, monkeypatch
) -> None:
    selection_rows = [
        _api_row("GCF_000001.1", "Francisella", "tularensis"),
        _api_row("GCF_000002.1", "Francisella", "tularensis", is_rep=False),
        _api_row("GCF_000003.1", "Francisella", "tularensis", is_rep=False),
    ]
    parent_rows = selection_rows + [_api_row("GCF_000010.1", "Francisella", "philomiragia")]
    cards = {
        "GCF_000001.1": _card(90.0, 2.0),  # 80
        "GCF_000002.1": _card(99.0, 0.5),  # 96.5, best
        "GCF_000003.1": _card(95.0, 1.0),  # 90
        "GCF_000010.1": _card(98.0, 0.5),
    }
    fake = _fake_api_with_cards(selection_rows, parent_rows, cards)
    card_requests: list[str] = []

    def counting(path, params=None):
        if path.endswith("/card"):
            card_requests.append(path)
        return fake(path, params)

    monkeypatch.setattr(gtdb_api, "_api_get", counting)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(
        ctx,
        MetadataParams(
            dataset="all",
            level="species",
            source="api",
            limit=2,
            target_genus="Francisella",
            target_species="tularensis",
        ),
    )
    rows = _read_selection(ctx.workdir)
    kept = {r["accession"] for r in rows if r["is_outgroup"] not in ("1", "True", "true")}
    assert kept == {"GCF_000002.1", "GCF_000003.1"}
    # Every candidate's card was fetched before the cut (3), plus the outgroup's (1).
    assert len(card_requests) == 4


def test_api_card_that_was_rate_limited_is_retried_in_a_second_pass(monkeypatch, caplog) -> None:
    # A 429 that outlives the session's retries must not leave the genome
    # unscored when a later, slower attempt would have succeeded.
    attempts: dict[str, int] = {}

    def flaky_card(accession: str) -> dict:
        attempts[accession] = attempts.get(accession, 0) + 1
        if accession == "GCF_000002.1" and attempts[accession] == 1:
            raise WorkdirError(
                "HTTP request failed: .../card (429 Client Error: Too Many Requests)"
            )
        return _card(95.0, 1.0)

    monkeypatch.setattr(gtdb_api, "_api_card", flaky_card)
    monkeypatch.setattr(gtdb_api.time, "sleep", lambda s: None)
    logger = logging.getLogger("test.retry.pass")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        quality = gtdb_api._api_quality_by_accession(
            ["GCF_000001.1", "GCF_000002.1", "GCF_000003.1"], logger
        )
    assert quality["GCF_000002.1"] == (95.0, 1.0)
    assert attempts["GCF_000002.1"] == 2
    assert not [r for r in caplog.records if r.levelname == "WARNING"]


def test_api_card_failing_twice_is_unscored_with_one_warning(monkeypatch, caplog) -> None:
    def broken_card(accession: str) -> dict:
        if accession == "GCF_000002.1":
            raise WorkdirError("HTTP request failed: .../card (429 Client Error)")
        return _card(95.0, 1.0)

    monkeypatch.setattr(gtdb_api, "_api_card", broken_card)
    monkeypatch.setattr(gtdb_api.time, "sleep", lambda s: None)
    logger = logging.getLogger("test.retry.fail")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        quality = gtdb_api._api_quality_by_accession(["GCF_000001.1", "GCF_000002.1"], logger)
    assert quality["GCF_000002.1"] == (None, None)
    warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
    assert len(warnings) == 1 and "GCF_000002.1" in warnings[0]


def test_tsv_missing_gtdb_version_names_the_cli_flag(tmp_path, gtdb_tsv) -> None:
    # The CLI flag is --gtdb-version; --version is the global "show version"
    # option, so naming it in the error sends the user to the wrong flag.
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--gtdb-version"):
        metadata.run(ctx, _params(gtdb_tsv, version=None))


# --- deep audit 2026-10 (entry-net) -------------------------------------------


@pytest.mark.parametrize("release", ["abc.def", "232.x", "232", "232.0.1"])
def test_malformed_release_is_user_input(tmp_path, gtdb_tsv, release) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--release like"):
        metadata.run(ctx, _params(gtdb_tsv, release=release))


def test_unknown_gtdb_version_is_user_input(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="Invalid --gtdb-version 'bogus'"):
        metadata.run(ctx, _params(gtdb_tsv, version="bogus"))


def test_missing_metadata_path_is_refused_not_downloaded(tmp_path, monkeypatch) -> None:
    def no_download(*a, **k):
        raise AssertionError("a missing --metadata-path must not trigger a download")

    monkeypatch.setattr(core_http, "download", no_download)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--metadata-path not found"):
        metadata.run(ctx, _params(tmp_path / "typo.tsv.gz"))


def test_api_source_warns_about_every_ignored_table_flag(tmp_path, monkeypatch, caplog) -> None:
    rows = [_api_row("GCF_000001.1", "Francisella", "tularensis")]
    parent = rows + [_api_row("GCF_000010.1", "Francisella", "philomiragia")]
    fake = _fake_api_with_cards(rows, parent, _cards_for(parent))
    monkeypatch.setattr(gtdb_api, "_api_get", fake)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    ctx.logger.addHandler(caplog.handler)  # the workdir logger does not propagate
    params = MetadataParams(
        dataset="all",
        level="species",
        source="api",
        target_genus="Francisella",
        target_species="tularensis",
        nodownload=True,
        metadata_path="/some/table.tsv.gz",
    )
    with caplog.at_level(logging.WARNING):
        metadata.run(ctx, params)
    assert "ignoring --metadata-path, --nodownload" in caplog.text


def _cards_for(rows: list[dict]) -> dict[str, dict]:
    return {r["gid"]: _card(99.0, 1.0) for r in rows}


def test_tsv_limit_never_picks_an_outgroup_from_the_target_taxon(tmp_path) -> None:
    """Under --limit the target genomes the cap left out are not outgroup candidates.

    Live (r232, -l genus -tg Francisella --limit 5) picked F. guangzhouensis, a
    species representative of the target genus, as the outgroup.
    """
    rows = [
        ("GCF_000001.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
        ("GCF_000010.1", "GCF_000010.1", _tax("Francisellaceae", "Francisella", "philomiragia")),
        ("GCF_000030.1", "GCF_000030.1", _tax("Francisellaceae", "Fangia", "hongkongensis")),
    ]
    path = tmp_path / "t.tsv.gz"
    header = "\t".join(
        [
            "accession",
            "gtdb_genome_representative",
            "gtdb_taxonomy",
            "checkm2_completeness",
            "checkm2_contamination",
        ]
    )
    lines = [header]
    for acc, rep, tax in rows:
        quality = "99\t0" if acc == "GCF_000001.1" else "90\t5"
        lines.append(f"RS_{acc}\tRS_{rep}\t{tax}\t{quality}")
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(path, level="genus", target_species=None, limit=1))
    rows = _read_selection(ctx.workdir)
    ingroup = [r for r in rows if r["is_outgroup"] == "0"]
    outgroup = [r for r in rows if r["is_outgroup"] == "1"]
    assert len(ingroup) == 1 and ingroup[0]["genus"] == "Francisella"
    # the other Francisella representative is left out, not made the outgroup
    assert [r["accession"] for r in outgroup] == ["GCF_000030.1"]


def test_tsv_outgroup_accession_inside_the_selection_is_refused(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="part of the selection"):
        metadata.run(ctx, _params(gtdb_tsv, outgroup_accession="GCF_000001.1"))


def test_named_outgroup_inside_the_target_taxon_is_refused(tmp_path, gtdb_tsv) -> None:
    """GCF_000002.1 is a non-representative tularensis genome: under -d rep it is
    not selected, but it lies inside the target species, so it is no outgroup."""
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="inside the target species tularensis"):
        metadata.run(ctx, _params(gtdb_tsv, dataset="rep", outgroup_accession="GCF_000002.1"))


def test_named_outgroup_left_out_by_limit_is_refused(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="inside the target"):
        metadata.run(ctx, _params(gtdb_tsv, limit=1, outgroup_accession="GCF_000003.1"))


def test_named_outgroup_of_a_sister_species_is_accepted(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv, dataset="rep", outgroup_accession="GCF_000010.1"))
    rows = _read_selection(ctx.workdir)
    assert [r["accession"] for r in rows if r["is_outgroup"] == "1"] == ["GCF_000010.1"]


def test_api_named_outgroup_inside_the_target_taxon_is_refused(tmp_path, monkeypatch) -> None:
    rows = [_api_row("GCF_000001.1", "Francisella", "tularensis")]
    cards = _cards_for(rows)
    cards["GCF_000002.1"] = {
        "metadataTaxonomy": {
            "gtdbFamily": "f__Francisellaceae",
            "gtdbGenus": "g__Francisella",
            "gtdbSpecies": "s__Francisella tularensis",
        },
        "metadata_gene": {},
    }
    monkeypatch.setattr(gtdb_api, "_api_get", _fake_api_with_cards(rows, rows, cards))
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(
        dataset="rep",
        level="species",
        source="api",
        target_genus="Francisella",
        target_species="tularensis",
        outgroup_accession="GCF_000002.1",
    )
    with pytest.raises(UserInputError, match="inside the target species"):
        metadata.run(ctx, params)


def test_api_outgroup_accession_inside_the_selection_is_refused(tmp_path, monkeypatch) -> None:
    rows = [_api_row("GCF_000001.1", "Francisella", "tularensis")]
    monkeypatch.setattr(gtdb_api, "_api_get", _fake_api_with_cards(rows, rows, _cards_for(rows)))
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(
        dataset="all",
        level="species",
        source="api",
        target_genus="Francisella",
        target_species="tularensis",
        outgroup_accession="GCF_000001.1",
    )
    with pytest.raises(UserInputError, match="part of the selection"):
        metadata.run(ctx, params)


def _api_404(path, params=None):
    raise core_http.HTTPStatusError(f"HTTP request failed: {path} (404)", 404)


def test_api_unknown_taxon_is_user_input(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(gtdb_api, "_api_get", _api_404)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(dataset="rep", level="genus", source="api", target_genus="Nogenus")
    with pytest.raises(UserInputError, match="GTDB API has no taxon g__Nogenus"):
        metadata.run(ctx, params)


def test_api_unknown_outgroup_is_user_input(tmp_path, monkeypatch) -> None:
    rows = [_api_row("GCF_000001.1", "Francisella", "tularensis")]
    cards = _cards_for(rows)
    inner = _fake_api_with_cards(rows, rows, cards)

    def fake(path, params=None):
        if "GCF_999" in path:
            return _api_404(path)
        return inner(path, params)

    monkeypatch.setattr(gtdb_api, "_api_get", fake)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(
        dataset="all",
        level="species",
        source="api",
        target_genus="Francisella",
        target_species="tularensis",
        outgroup_accession="GCF_999.1",
    )
    with pytest.raises(UserInputError, match="not found in the GTDB API"):
        metadata.run(ctx, params)


def test_api_server_error_stays_a_workdir_error(tmp_path, monkeypatch) -> None:
    def fake(path, params=None):
        raise core_http.HTTPStatusError("HTTP request failed (503)", 503)

    monkeypatch.setattr(gtdb_api, "_api_get", fake)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(dataset="rep", level="genus", source="api", target_genus="Francisella")
    with pytest.raises(WorkdirError) as info:
        metadata.run(ctx, params)
    assert not isinstance(info.value, UserInputError)


def test_unknown_release_names_both_layouts_as_user_input(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(core_http, "download", lambda url, dest, **k: _api_404(url))
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = _params(Path("unused"), metadata_path=None, release="999.0")
    with pytest.raises(UserInputError, match="no bac120 metadata table for release 999.0"):
        metadata.run(ctx, params)


def test_network_failure_is_reported_without_trying_the_next_layout(tmp_path, monkeypatch) -> None:
    calls = []

    def fake_download(url, dest, **k):
        calls.append(url)
        raise WorkdirError(f"Download failed: {url} (Connection refused)")

    monkeypatch.setattr(core_http, "download", fake_download)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="Connection refused"):
        metadata.run(ctx, _params(Path("unused"), metadata_path=None))
    assert len(calls) == 1


def test_checksum_mismatch_is_reported_as_such(tmp_path, gtdb_tsv, monkeypatch) -> None:
    monkeypatch.setattr(core_http, "download", lambda url, dest, **k: Path(dest).write_bytes(b"x"))

    def bad_md5(dest, url, **k):
        raise WorkdirError("Checksum mismatch for bac120_metadata_r232.tsv.gz")

    monkeypatch.setattr(core_http, "verify_md5_manifest", bad_md5)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="Checksum mismatch"):
        metadata.run(ctx, _params(gtdb_tsv, metadata_path=None))
    assert not (ctx.workdir / "bac120_metadata_r232.tsv.gz").exists()


def _download_fixture(gtdb_tsv):
    def fake_download(url, dest, **k):
        Path(dest).write_bytes(gtdb_tsv.read_bytes())

    return fake_download


def test_nodownload_refuses_a_table_from_another_minor_release(
    tmp_path, gtdb_tsv, monkeypatch
) -> None:
    """Table names carry only the major release (r214 serves 214.0 and 214.1)."""
    monkeypatch.setattr(core_http, "download", _download_fixture(gtdb_tsv))
    monkeypatch.setattr(core_http, "verify_md5_manifest", lambda *a, **k: True)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv, metadata_path=None, nodownload=True))
    assert (ctx.workdir / "bac120_metadata_r232.release").read_text().strip() == "232.0"
    with pytest.raises(UserInputError, match="downloaded for GTDB release 232.0, not 232.1"):
        metadata.run(ctx, _params(gtdb_tsv, metadata_path=None, nodownload=True, release="232.1"))
    # the matching release is reused without a download
    monkeypatch.setattr(core_http, "download", _api_404)
    assert metadata.run(ctx, _params(gtdb_tsv, metadata_path=None, nodownload=True)) == 3


def test_nodownload_reuses_an_unmarked_table_with_a_warning(tmp_path, gtdb_tsv, caplog) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    (ctx.workdir / "bac120_metadata_r232.tsv.gz").write_bytes(gtdb_tsv.read_bytes())
    ctx.logger.addHandler(caplog.handler)  # the workdir logger does not propagate
    with caplog.at_level(logging.WARNING):
        assert metadata.run(ctx, _params(gtdb_tsv, metadata_path=None, nodownload=True)) == 3
    assert "Cannot confirm which GTDB release" in caplog.text


def test_nodownload_reuses_a_legacy_tarball_without_requesting_the_modern_layout(
    tmp_path, gtdb_tsv, monkeypatch
) -> None:
    import tarfile

    monkeypatch.setattr(core_http, "download", _api_404)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    member = tmp_path / "bac120_metadata_r232.tsv"
    with gzip.open(gtdb_tsv, "rb") as fo:
        member.write_bytes(fo.read())
    with tarfile.open(ctx.workdir / "bac120_metadata_r232.tar.gz", "w:gz") as tar:
        tar.add(member, arcname=member.name)
    (ctx.workdir / "bac120_metadata_r232.release").write_text("232.0\n")
    assert metadata.run(ctx, _params(gtdb_tsv, metadata_path=None, nodownload=True)) == 3


def test_tsv_selection_records_the_gtdb_representative_flag(tmp_path, gtdb_tsv) -> None:
    """selection.tsv and the manifest carry gtdb_genome_representative == accession."""
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv))
    flags = {r["accession"]: r["gtdb_representative"] for r in _read_selection(ctx.workdir)}
    assert flags == {
        "GCF_000001.1": "1",
        "GCF_000002.1": "0",
        "GCF_000003.1": "1",
        "GCF_000010.1": "1",
    }
    reps = {g.accession for g in ctx.manifest.all_genomes() if g.gtdb_representative}
    assert reps == {"GCF_000001.1", "GCF_000003.1"}


def test_api_selection_records_the_gtdb_representative_flag(tmp_path, monkeypatch) -> None:
    selection_rows = [
        _api_row("GCF_000001.1", "Francisella", "tularensis"),
        _api_row("GCF_000002.1", "Francisella", "tularensis", is_rep=False),
    ]
    parent_rows = selection_rows + [_api_row("GCF_000010.1", "Francisella", "philomiragia")]
    monkeypatch.setattr(
        gtdb_api,
        "_api_get",
        _fake_api_with_cards(selection_rows, parent_rows, _cards_for(parent_rows)),
    )
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(
        ctx,
        MetadataParams(
            dataset="all",
            level="species",
            source="api",
            target_genus="Francisella",
            target_species="tularensis",
        ),
    )
    flags = {r["accession"]: r["gtdb_representative"] for r in _read_selection(ctx.workdir)}
    assert flags == {"GCF_000001.1": "1", "GCF_000002.1": "0", "GCF_000010.1": "1"}
    reps = {g.accession for g in ctx.manifest.all_genomes() if g.gtdb_representative}
    assert reps == {"GCF_000001.1"}


def test_parse_metadata_without_the_genbank_column_leaves_the_ncbi_accession_empty(
    tmp_path,
) -> None:
    # A table without ncbi_genbank_assembly_accession used to yield column 0,
    # the prefixed GTDB accession, as the NCBI accession.
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    tax = _tax("Francisellaceae", "Francisella", "tularensis")
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("accession\tgtdb_genome_representative\tgtdb_taxonomy\n")
        fo.write(f"RS_GCF_000001.1\tRS_GCF_000001.1\t{tax}\n")
    accessions = gtdb_table.parse_metadata(path, _params(path), logging.getLogger("test"))
    assert accessions["GCF_000001.1"]["accession_ncbi"] == ""
    assert accessions["GCF_000001.1"]["details"].accession_ncbi == ""
