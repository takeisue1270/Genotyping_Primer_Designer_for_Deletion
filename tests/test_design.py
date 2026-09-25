"""End-to-end behaviour: what a run reports, and what it refuses to imply."""

from __future__ import annotations

import dataclasses

import pytest

from genotyping_primers.config import DesignConfig
from genotyping_primers.design import (
    DesignError,
    choose_checker,
    design,
    prepare,
    resolve_span,
    target_from_file,
    with_insert,
)
from genotyping_primers.report import text_report, write_maps, write_tsv
from genotyping_primers.seqio import parse_genbank
from genotyping_primers.specificity import FastaChecker, GGGenomeChecker, LocalChecker, NullChecker


def _design(target, spec, config, insert=""):
    working, deletion = prepare(target, *resolve_span(target, spec))
    deletion = with_insert(deletion, insert, "insert" if insert else "")
    return design(working, deletion, config)


# --- which backend runs ------------------------------------------------------


def test_auto_uses_the_supplied_sequence_when_nothing_else_is_offered(target):
    checker = choose_checker(target, DesignConfig())
    assert isinstance(checker, LocalChecker)
    assert "cannot show a primer is unique in a genome" in checker.describe()


def test_auto_prefers_gggenome_when_the_species_has_a_database(target):
    checker = choose_checker(target, DesignConfig(species="mus_musculus"))
    assert isinstance(checker, GGGenomeChecker)
    assert checker.database == "mm39"


def test_auto_prefers_a_local_genome_fasta_over_the_network(target, tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
    checker = choose_checker(
        target, DesignConfig(species="mus_musculus", genome_fasta=str(fasta))
    )
    assert isinstance(checker, FastaChecker)


def test_an_unknown_species_says_how_to_proceed(target):
    with pytest.raises(DesignError, match="gggenome-db"):
        choose_checker(target, DesignConfig(specificity="gggenome", species="felis_catus"))


def test_none_skips_the_check_without_claiming_success(target):
    assert isinstance(choose_checker(target, DesignConfig(specificity="none")), NullChecker)


# --- what the notes say ------------------------------------------------------


def test_a_resolvable_deletion_produces_no_gel_warning(target, config):
    result = _design(target, "4001-4500", config)
    assert result.pairs
    assert not any("agarose" in note for note in result.notes)


def test_a_small_deletion_is_called_out_as_unresolvable(target, config):
    result = _design(target, "4001-4030", config)
    assert result.pairs
    assert any("will not resolve" in note for note in result.notes)
    assert any("sequencing across the junction" in note for note in result.notes)


def test_an_insert_the_same_length_as_the_deletion_cannot_be_genotyped_by_size(target, config):
    result = _design(target, "4001-4034", config, insert="A" * 34)
    assert any("no PCR across it can" in note for note in result.notes)


def test_an_insert_longer_than_the_deletion_makes_the_edited_band_the_larger_one(target, config):
    result = _design(target, "4001-4100", config, insert="AC" * 150)
    assert result.pairs
    assert result.best.ko_product > result.best.wt_product
    assert any("the *larger*" in note for note in result.notes)


def test_zero_arm_warns_about_a_repair_template(target):
    result = _design(target, "4001-4500", DesignConfig(specificity="local", homology_arm=0))
    assert any("--arm" in note for note in result.notes)


def test_the_specificity_note_names_the_backend(target, config):
    result = _design(target, "4001-4500", config)
    assert any(note.startswith("Specificity:") for note in result.notes)
    assert result.specificity_backend == "local"


# --- outputs -----------------------------------------------------------------


def test_the_report_states_both_band_sizes(target, config):
    result = _design(target, "4001-4500", config)
    report = text_report(result, config)
    assert f"unedited {result.best.wt_product:,} bp" in report
    assert f"edited {result.best.ko_product:,} bp" in report
    assert "Recommended pair" in report


def test_the_tsv_has_one_row_per_primer(tmp_path, target, config):
    result = _design(target, "4001-4500", config)
    path = write_tsv(tmp_path / "primers.tsv", result)
    rows = path.read_text(encoding="utf-8").strip().splitlines()
    assert len(rows) == 1 + 2 * len(result.pairs)
    assert rows[0].split("\t")[:4] == ["pair", "primer", "name", "sequence"]
    assert result.best.forward.sequence in rows[1]


def test_the_maps_carry_the_two_alleles_and_the_primers(tmp_path, target, config):
    result = _design(target, "4001-4500", config, insert="GGATCCAAGCTTGA")
    unedited, edited = write_maps(tmp_path, result, config)

    _, before, _ = parse_genbank(unedited.read_text(encoding="utf-8"))
    _, after, _ = parse_genbank(edited.read_text(encoding="utf-8"))
    assert before == result.target.sequence
    assert after == result.edited_sequence
    assert len(before) - len(after) == result.deletion.net_change

    text = unedited.read_text(encoding="utf-8")
    assert "deleted region (500 bp)" in text
    assert "homology arm left (43 bp)" in text
    assert result.best.forward.name in text


def test_a_plasmid_deletion_over_the_origin_designs_and_reports_in_file_coordinates(plasmid):
    config = DesignConfig(specificity="local", search_span=600, alternatives=0)
    result = _design(plasmid, "3900-200", config)
    assert result.pairs, result.notes

    report = text_report(result, config)
    assert "3,900-200" in report
    assert "circular" in report
    # The primers must be real sites in the plasmid as the user supplied it.
    original = plasmid.sequence + plasmid.sequence
    for primer in (result.best.forward, result.best.reverse):
        assert primer.top_strand() in original


def test_a_rotated_plasmid_map_says_where_base_one_went(tmp_path, plasmid):
    config = DesignConfig(specificity="local", search_span=600, alternatives=0)
    result = _design(plasmid, "3900-200", config)
    unedited, edited = write_maps(tmp_path, result, config)
    for path in (unedited, edited):
        assert "rotated: base 1 here is base" in path.read_text(encoding="utf-8")


def test_a_linear_map_is_not_marked_as_rotated(tmp_path, target, config):
    result = _design(target, "4001-4500", config)
    unedited, _ = write_maps(tmp_path, result, config)
    assert "rotated" not in unedited.read_text(encoding="utf-8")


def test_reading_a_target_from_a_file_keeps_its_name(tmp_path, make_sequence):
    path = tmp_path / "locus.fa"
    path.write_text(">my_locus\n" + make_sequence(3_000, seed=2) + "\n", encoding="utf-8")
    target, notes = target_from_file(str(path))
    assert target.name == "my_locus"
    assert notes == []
    assert len(target) == 3_000


def test_an_ambiguous_base_is_reported_and_never_used(tmp_path, make_sequence):
    sequence = make_sequence(3_000, seed=4)
    spiked = sequence[:1_500] + "N" * 20 + sequence[1_520:]
    path = tmp_path / "gappy.fa"
    path.write_text(">gappy\n" + spiked + "\n", encoding="utf-8")

    target, notes = target_from_file(str(path))
    assert any("ambiguity code" in note for note in notes)

    config = DesignConfig(specificity="local", search_span=800, alternatives=1)
    working, deletion = prepare(target, *resolve_span(target, "1701-1900"))
    result = design(working, deletion, config, extra_notes=notes)
    assert result.pairs
    for pair in result.pairs:
        assert "N" not in pair.forward.sequence
        assert "N" not in pair.reverse.sequence


def test_the_deletion_cannot_cover_the_whole_plasmid(plasmid):
    with pytest.raises(DesignError, match="whole sequence"):
        prepare(plasmid, 0, 4_000)


def test_an_insert_must_be_dna(target):
    _, deletion = prepare(target, *resolve_span(target, "4001-4500"))
    with pytest.raises(DesignError, match="A/C/G/T"):
        with_insert(deletion, "ACGTX")
    assert with_insert(deletion, "acg t").insert == "ACGT"


def test_config_rejects_impossible_settings():
    with pytest.raises(ValueError, match="homology_arm"):
        DesignConfig(homology_arm=-1)
    with pytest.raises(ValueError, match="specificity"):
        DesignConfig(specificity="magic")
    with pytest.raises(ValueError, match="tm_range"):
        DesignConfig(tm_range=(62.0, 58.0))


def test_dataclasses_stay_frozen(target, config):
    result = _design(target, "4001-4500", config)
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.best.forward.name = "nope"
