"""The medoid reference: chosen from the sketch ANI matrix, deterministic."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest

from repgenr.tree import medoid as mod

_LOG = logging.getLogger("test-medoid")


def test_pick_medoid_takes_the_largest_summed_ani() -> None:
    names = ["c", "a", "b"]
    ani = np.array([[1.0, 0.90, 0.95], [0.90, 1.0, 0.99], [0.95, 0.99, 1.0]])
    name, mean = mod.pick_medoid(names, ani)
    assert name == "b" and mean == pytest.approx(97.0)


def test_pick_medoid_breaks_ties_by_name() -> None:
    names = ["z", "a", "m"]
    ani = np.full((3, 3), 0.9)
    np.fill_diagonal(ani, 1.0)
    assert mod.pick_medoid(names, ani)[0] == "a"
    assert mod.pick_medoid(["only"], np.array([[1.0]])) == ("only", 100.0)


def test_choose_medoid_needs_every_sketch(tmp_path: Path, caplog) -> None:
    genomes = [tmp_path / "a.fasta", tmp_path / "b.fasta"]
    with caplog.at_level(logging.INFO, logger="test-medoid"):
        assert (
            mod.choose_medoid(
                genomes, lambda gs: {genomes[0]: tmp_path / "a.sig"}, tmp_path, 1, _LOG
            )
            is None
        )
    assert any("have no sketch" in r.getMessage() for r in caplog.records)
    assert mod.choose_medoid(genomes, None, tmp_path, 1, _LOG) is None
    assert mod.choose_medoid(genomes[:1], lambda gs: {}, tmp_path, 1, _LOG) is None


def _genomes(tmp_path: Path) -> list[Path]:
    out = []
    for name in ("Fam_Gen_sp_GCA_1.fasta", "Fam_Gen_sp_GCA_2.fasta", "Fam_Gen_sp_GCA_3.fasta"):
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        out.append(p)
    return out


def test_resolve_reference_uses_the_medoid_by_default(tmp_path: Path, monkeypatch, caplog) -> None:
    genomes = _genomes(tmp_path)
    seen: dict = {}

    def fake_choose(gs, provider, work_dir, threads, logger):
        seen.update(n=len(gs), provider=provider, threads=threads)
        return mod.MedoidChoice(genome=gs[2], mean_ani=99.1, n_genomes=len(gs))

    monkeypatch.setattr(mod, "choose_medoid", fake_choose)
    with caplog.at_level(logging.INFO, logger="test-medoid"):
        ref = mod.resolve_reference(
            None,
            genomes,
            None,
            _LOG,
            sketch_provider=lambda gs: {},
            scratch_dir=tmp_path,
            threads=4,
        )
    assert ref.genome == genomes[2] and ref.rule == "medoid" and ref.mean_ani == 99.1
    assert seen == {"n": 3, "provider": seen["provider"], "threads": 4}
    assert any("the medoid of 3 genomes" in r.getMessage() for r in caplog.records)
    # The explicit rules.
    assert (
        mod.resolve_reference("medoid", genomes, None, _LOG, sketch_provider=lambda gs: {}).genome
        == genomes[2]
    )
    assert (
        mod.resolve_reference("first", genomes, None, _LOG, sketch_provider=lambda gs: {}).rule
        == "first"
    )
    assert mod.resolve_reference(genomes[1].name, genomes, None, _LOG).genome == genomes[1]


def test_resolve_reference_falls_back_to_the_first_genome_with_a_warning(
    tmp_path: Path, monkeypatch, caplog
) -> None:
    genomes = _genomes(tmp_path)
    monkeypatch.setattr(mod, "choose_medoid", lambda *a, **k: None)
    with caplog.at_level(logging.WARNING, logger="test-medoid"):
        assert mod.resolve_reference(None, genomes, None, _LOG).genome == genomes[0]
    warned = [r.getMessage() for r in caplog.records]
    assert warned and "no medoid" in warned[0] and genomes[0].name in warned[0]
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="test-medoid"):
        mod.resolve_reference("first", genomes, None, _LOG)
    assert "no medoid" not in caplog.records[0].getMessage()


def test_reference_note_round_trips(tmp_path: Path) -> None:
    mod.note_reference(
        tmp_path / "align", mod.ReferenceChoice(tmp_path / "Fam_Gen_sp_GCA_2.fasta", "medoid", 99.5)
    )
    assert mod.read_reference_note(tmp_path / "align") == "Fam_Gen_sp_GCA_2.fasta"
    assert (tmp_path / "align" / "reference.txt").read_text().split("\t")[1:] == [
        "medoid",
        "99.50\n",
    ]
    mod.clear_reference_note(tmp_path / "align")
    assert mod.read_reference_note(tmp_path / "align") is None


def test_choose_medoid_refuses_beyond_the_dense_limit(tmp_path: Path, caplog) -> None:
    from repgenr.dereplicators.sourmash import _DENSE_MAX_GENOMES

    genomes = [tmp_path / f"g{i}.fasta" for i in range(_DENSE_MAX_GENOMES + 1)]
    calls: list = []
    with caplog.at_level(logging.WARNING, logger="test-medoid"):
        assert (
            mod.choose_medoid(genomes, lambda gs: calls.append(1) or {}, tmp_path, 1, _LOG) is None
        )
    assert calls == [] and any("dense ANI matrix" in r.getMessage() for r in caplog.records)


def test_msa_key_follows_the_chosen_genome(tmp_path: Path) -> None:
    from repgenr.stages.phylo import PhyloParams, _msa_key

    genomes = _genomes(tmp_path)
    params = PhyloParams(msa_source="aligner")
    assert _msa_key(genomes, None, params, genomes[0]) != _msa_key(
        genomes, None, params, genomes[1]
    )
    assert _msa_key(genomes, None, params, genomes[0]) == _msa_key(
        genomes, None, PhyloParams(msa_source="aligner", reference="medoid"), genomes[0]
    )


def test_explicit_medoid_without_sketches_is_refused(tmp_path: Path, monkeypatch) -> None:
    from repgenr.core.errors import UserInputError

    genomes = _genomes(tmp_path)
    monkeypatch.setattr(mod, "choose_medoid", lambda *a, **k: None)
    with pytest.raises(UserInputError, match="sketch for every genome"):
        mod.resolve_reference("medoid", genomes, None, _LOG, sketch_provider=lambda gs: {})


def test_read_reference_choice_returns_rule_and_score(tmp_path: Path) -> None:
    mod.note_reference(tmp_path, mod.ReferenceChoice(Path("x_GCA_1.fasta"), "medoid", 98.25))
    choice = mod.read_reference_choice(tmp_path)
    assert choice == mod.ReferenceChoice(Path("x_GCA_1.fasta"), "medoid", 98.25)
    mod.note_reference(tmp_path, mod.ReferenceChoice(Path("x_GCA_1.fasta"), "first-fallback"))
    assert mod.read_reference_choice(tmp_path).mean_ani is None
    # A note from before the rule column: the genome alone.
    (tmp_path / mod.REFERENCE_NOTE).write_text("y_GCA_2.fasta\n")
    assert mod.read_reference_choice(tmp_path) == mod.ReferenceChoice(Path("y_GCA_2.fasta"), "name")
