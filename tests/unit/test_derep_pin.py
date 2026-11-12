"""dereplicate --pin: pinned genomes become representatives."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.core.errors import UserInputError
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_FAIL_QC,
    STATUS_REPRESENTATIVE,
    DerepResult,
)
from repgenr.stages.derep_pin import (
    pin_representatives,
    pinned_names,
    read_pin_file,
    split_genomes,
)

LOG = logging.getLogger("test")


def _result(clusters: dict[str, list[str]], failed: tuple[str, ...] = ()) -> DerepResult:
    status = {}
    for rep, members in clusters.items():
        status[rep] = STATUS_REPRESENTATIVE
        status.update({m: STATUS_CONTAINED for m in members})
    status.update({f: STATUS_FAIL_QC for f in failed})
    return DerepResult(
        representatives=[Path("/g") / rep for rep in clusters],
        clusters={rep: list(members) for rep, members in clusters.items()},
        genome_status=status,
    )


def _ani(table: dict[frozenset[str], float]):
    def ani(a: str, b: str) -> float | None:
        return table.get(frozenset((a, b)))

    return ani


def test_pinned_member_becomes_representative_and_takes_its_nearest_members() -> None:
    result = _result({"r.fa": ["p.fa", "m1.fa", "m2.fa"], "s.fa": ["t.fa"]})
    ani = _ani(
        {
            frozenset(("m1.fa", "r.fa")): 0.9990,
            frozenset(("m1.fa", "p.fa")): 0.9999,
            frozenset(("m2.fa", "r.fa")): 0.9999,
            frozenset(("m2.fa", "p.fa")): 0.9990,
        }
    )
    out, promoted = pin_representatives(result, ["p.fa"], ani, LOG)
    assert promoted == ["p.fa"]
    assert out.clusters == {"r.fa": ["m2.fa"], "p.fa": ["m1.fa"], "s.fa": ["t.fa"]}
    assert [p.name for p in out.representatives] == ["r.fa", "p.fa", "s.fa"]
    assert out.representatives[1] == Path("/g/p.fa")
    assert out.genome_status["p.fa"] == STATUS_REPRESENTATIVE
    assert out.genome_status["m1.fa"] == STATUS_CONTAINED


def test_several_pins_in_one_cluster_all_become_representatives() -> None:
    result = _result({"r.fa": ["p1.fa", "p2.fa", "m.fa"]})
    ani = _ani({frozenset(("m.fa", "p2.fa")): 0.9999, frozenset(("m.fa", "p1.fa")): 0.999})
    out, promoted = pin_representatives(result, ["p1.fa", "p2.fa"], ani, LOG)
    assert promoted == ["p1.fa", "p2.fa"]
    assert out.clusters == {"r.fa": [], "p1.fa": [], "p2.fa": ["m.fa"]}


def test_member_without_ani_stays_with_original_representative(caplog) -> None:
    result = _result({"r.fa": ["p.fa", "m.fa"]})
    with caplog.at_level(logging.WARNING):
        out, _ = pin_representatives(result, ["p.fa"], _ani({}), LOG)
    assert out.clusters == {"r.fa": ["m.fa"], "p.fa": []}
    assert "no ANI estimate" in caplog.text


def test_equal_ani_keeps_the_original_representative() -> None:
    result = _result({"r.fa": ["p.fa", "m.fa"]})
    ani = _ani({frozenset(("m.fa", "r.fa")): 0.9995, frozenset(("m.fa", "p.fa")): 0.9995})
    out, _ = pin_representatives(result, ["p.fa"], ani, LOG)
    assert out.clusters["r.fa"] == ["m.fa"]


def test_pinned_representative_is_unchanged() -> None:
    result = _result({"r.fa": ["m.fa"]})
    out, promoted = pin_representatives(result, ["r.fa"], _ani({}), LOG)
    assert promoted == []
    assert out.clusters == result.clusters


def test_fail_qc_pin_is_left_out_with_a_warning(caplog) -> None:
    result = _result({"r.fa": ["m.fa"]}, failed=("q.fa",))
    with caplog.at_level(logging.WARNING):
        out, promoted = pin_representatives(result, ["q.fa"], _ani({}), LOG)
    assert promoted == []
    assert out.genome_status["q.fa"] == STATUS_FAIL_QC
    assert "q.fa" in caplog.text


def test_every_genome_keeps_one_cluster() -> None:
    result = _result({"r.fa": ["p.fa", "a.fa", "b.fa"], "s.fa": ["c.fa"]})
    out, _ = pin_representatives(result, ["p.fa"], _ani({}), LOG)
    placed = [g for rep, members in out.clusters.items() for g in (rep, *members)]
    assert sorted(placed) == sorted(["r.fa", "p.fa", "a.fa", "b.fa", "s.fa", "c.fa"])


def test_split_genomes_lists_only_clusters_holding_a_pin() -> None:
    result = _result({"r.fa": ["p.fa", "m.fa"], "s.fa": ["t.fa"]})
    assert split_genomes(result, ["p.fa"]) == ["r.fa", "p.fa", "m.fa"]


def test_read_pin_file_skips_comments_and_takes_the_first_column(tmp_path: Path) -> None:
    path = tmp_path / "pins.tsv"
    path.write_text("# strains\nGCF_1.1\tbiovar 1\n\n  GCF_2.1  \n", encoding="utf-8")
    assert read_pin_file(path) == ["GCF_1.1", "GCF_2.1"]


def test_read_pin_file_refuses_an_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "pins.txt"
    path.write_text("# nothing\n", encoding="utf-8")
    with pytest.raises(UserInputError, match="names no genome"):
        read_pin_file(path)


def test_pinned_names_resolves_accessions_and_names_the_flag() -> None:
    genomes = [Path("/g/Fam_Gen_sp_GCF_000007125.1.fasta"), Path("/g/other.fasta")]
    assert pinned_names(["GCF_000007125.1", "other"], genomes) == [
        "Fam_Gen_sp_GCF_000007125.1.fasta",
        "other.fasta",
    ]
    with pytest.raises(UserInputError, match="--pin names genome"):
        pinned_names(["GCF_999.1"], genomes)
