"""Reading sequences in, and writing GenBank back out."""

from __future__ import annotations

import pytest

from genotyping_primers.seqio import (
    Feature,
    SequenceError,
    format_genbank,
    parse_genbank,
    read_sequence,
)


def test_reads_fasta_and_reports_extra_records(tmp_path):
    path = tmp_path / "two.fa"
    path.write_text(">first desc here\nACGTACGTAC\nGTAC\n>second\nTTTTGGGG\n", encoding="utf-8")

    name, sequence, circular, notes = read_sequence(path)
    assert (name, sequence, circular) == ("first", "ACGTACGTACGTAC", False)
    assert any("2 records" in note for note in notes)

    name, sequence, _, _ = read_sequence(path, record="second")
    assert (name, sequence) == ("second", "TTTTGGGG")


def test_unknown_record_names_the_alternatives(tmp_path):
    path = tmp_path / "two.fa"
    path.write_text(">a\nACGT\n>b\nTGCA\n", encoding="utf-8")
    with pytest.raises(SequenceError, match="no record named"):
        read_sequence(path, record="c")


def test_reads_raw_sequence_and_uppercases(tmp_path):
    path = tmp_path / "plain.seq"
    path.write_text("acgt acgt\nACGTacgt\n", encoding="utf-8")
    name, sequence, circular, _ = read_sequence(path)
    assert (name, sequence, circular) == ("plain", "ACGTACGTACGTACGT", False)


def test_rejects_a_file_that_is_not_dna(tmp_path):
    path = tmp_path / "notes.txt"
    path.write_text("Design the KO primers for Nanog, please\n", encoding="utf-8")
    with pytest.raises(SequenceError):
        read_sequence(path)


def test_genbank_topology_is_honoured(tmp_path, make_sequence):
    sequence = make_sequence(300, seed=3)
    text = format_genbank("pUC_test", sequence, [], circular=True)
    path = tmp_path / "plasmid.gb"
    path.write_text(text, encoding="utf-8")

    name, parsed, circular, _ = read_sequence(path)
    assert circular is True
    assert parsed == sequence
    assert name == "pUC_test"


def test_genbank_round_trip_keeps_the_sequence_exactly(make_sequence):
    sequence = make_sequence(1_234, seed=5)
    features = [
        Feature("misc_feature", 100, 200, 1, [("label", "deleted region (100 bp)")]),
        Feature("primer_bind", 10, 34, -1, [("label", "test_R")]),
    ]
    text = format_genbank("target", sequence, features)

    name, parsed, circular = parse_genbank(text)
    assert parsed == sequence
    assert circular is False
    assert name == "target"
    # 1-based inclusive locations, and the reverse primer on the minus strand.
    assert "     misc_feature    101..200" in text
    assert "     primer_bind     complement(11..34)" in text


def test_genbank_locus_line_survives_an_awkward_name():
    text = format_genbank("chr11:69,000,000-69,001,000 window", "ACGT", [])
    assert text.startswith("LOCUS       chr11_69_000_000-69_001")
    assert " bp    DNA     linear" in text
