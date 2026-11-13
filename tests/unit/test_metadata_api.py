"""Offline tests for the GTDB-API metadata helpers (no network)."""

from __future__ import annotations

import pytest

from repgenr.core.errors import UserInputError
from repgenr.sources.gtdb_api import normalize_api_tax, target_taxon
from repgenr.stages.metadata import MetadataParams, _validate


def testtarget_taxon_strings() -> None:
    p = MetadataParams(dataset="rep", level="genus", source="api", target_genus="francisella")
    assert target_taxon(p) == "g__Francisella"

    p = MetadataParams(dataset="rep", level="family", source="api", target_family="francisellaceae")
    assert target_taxon(p) == "f__Francisellaceae"

    p = MetadataParams(
        dataset="rep",
        level="species",
        source="api",
        target_genus="francisella",
        target_species="tularensis",
    )
    assert target_taxon(p) == "s__Francisella tularensis"


def testnormalize_api_tax() -> None:
    row = {
        "gtdbFamily": "f__Francisellaceae",
        "gtdbGenus": "g__Francisella",
        "gtdbSpecies": "s__Francisella tularensis",
    }
    tax = normalize_api_tax(row)
    assert tax["family"] == "Francisellaceae"
    assert tax["genus"] == "Francisella"
    # genus removed from species, spaces dropped
    assert tax["species"] == "tularensis"


def test_validate_api_does_not_require_release() -> None:
    # api source: release/version optional
    _validate(
        MetadataParams(dataset="rep", level="genus", source="api", target_genus="francisella")
    )


def test_validate_tsv_requires_release() -> None:
    with pytest.raises(UserInputError):
        _validate(MetadataParams(dataset="rep", level="genus", source="tsv", target_genus="x"))


@pytest.mark.parametrize(
    ("genus", "species", "expected"),
    [
        ("francisella", "Tularensis", "s__Francisella tularensis"),
        ("Bacillus_A", "anthracis", "s__Bacillus_A anthracis"),
        ("Prevotella", "Copri_A", "s__Prevotella copri_A"),
    ],
)
def testtarget_taxon_keeps_gtdb_suffixes_and_lowers_the_epithet(genus, species, expected) -> None:
    p = MetadataParams(
        dataset="rep", level="species", source="api", target_genus=genus, target_species=species
    )
    assert target_taxon(p) == expected


def testtarget_taxon_keeps_a_genus_suffix() -> None:
    p = MetadataParams(dataset="rep", level="genus", source="api", target_genus="bacillus_A")
    assert target_taxon(p) == "g__Bacillus_A"
