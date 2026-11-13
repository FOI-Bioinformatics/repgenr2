"""One token sanitiser behind every producer of Family_genus_species names."""

from __future__ import annotations

from repgenr.core.contracts import sanitise_taxon_tokens
from repgenr.sources.gtdb_api import normalize_api_tax
from repgenr.sources.gtdb_table import parse_taxonomy


def test_tokens_drop_genus_prefix_spaces_and_underscores() -> None:
    assert sanitise_taxon_tokens("Francisellaceae", "Francisella", "Francisella tularensis") == (
        "Francisellaceae",
        "Francisella",
        "tularensis",
    )
    assert sanitise_taxon_tokens("Fam_ily A", "Gen us", "Gen us sp. nov_1") == (
        "Fam-ilyA",
        "Genus",
        "sp.nov-1",
    )


def test_tsv_and_api_paths_agree_with_the_sanitiser() -> None:
    raw = "d__Bacteria;f__Francisellaceae;g__Francisella;s__Francisella tularensis"
    tsv = parse_taxonomy(raw)
    api = normalize_api_tax(
        {
            "gtdbFamily": "f__Francisellaceae",
            "gtdbGenus": "g__Francisella",
            "gtdbSpecies": "s__Francisella tularensis",
        }
    )
    expected = sanitise_taxon_tokens("Francisellaceae", "Francisella", "Francisella tularensis")
    assert (tsv["family"], tsv["genus"], tsv["species"]) == expected
    assert (api["family"], api["genus"], api["species"]) == expected
