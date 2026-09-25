"""The command line, end to end, without touching the network."""

from __future__ import annotations

import pytest

from genotyping_primers.cli import main


@pytest.fixture
def locus(tmp_path, make_sequence):
    path = tmp_path / "locus.fa"
    path.write_text(">test_locus\n" + make_sequence(8_000, seed=7) + "\n", encoding="utf-8")
    return path


def test_a_plain_run_prints_both_band_sizes(locus, capsys):
    code = main(["--sequence", str(locus), "--delete", "4001-4500", "--specificity", "local"])
    out = capsys.readouterr().out
    assert code == 0
    assert "Recommended pair" in out
    assert "bands: unedited" in out
    assert "500 bp removed" in out


def test_outputs_are_written_when_asked_for(locus, tmp_path, capsys):
    outdir = tmp_path / "out"
    code = main(
        [
            "--sequence", str(locus),
            "--delete", "4001-4500",
            "--specificity", "local",
            "-o", str(outdir),
        ]
    )
    assert code == 0
    written = sorted(p.name for p in outdir.iterdir())
    assert written == [
        "test_locus_edited.gb",
        "test_locus_genotyping_primers.tsv",
        "test_locus_unedited.gb",
    ]
    assert "Written" in capsys.readouterr().out


def test_no_maps_writes_only_the_table(locus, tmp_path):
    outdir = tmp_path / "out"
    main(
        [
            "--sequence", str(locus), "--delete", "4001-4500",
            "--specificity", "local", "--no-maps", "-o", str(outdir),
        ]
    )
    assert [p.suffix for p in outdir.iterdir()] == [".tsv"]


def test_a_literal_insert_is_shown_and_counted(locus, capsys):
    main(
        [
            "--sequence", str(locus), "--delete", "4001-4500",
            "--insert", "GGATCCAAGCTTGA", "--specificity", "local",
        ]
    )
    out = capsys.readouterr().out
    assert "GGATCCAAGCTTGA" in out
    assert "486 bp shorter" in out  # 500 removed, 14 put back


def test_a_lower_case_insert_is_accepted(locus, capsys):
    main(
        [
            "--sequence", str(locus), "--delete", "4001-4500",
            "--insert", "acgtacgt", "--specificity", "local",
        ]
    )
    out = capsys.readouterr().out
    assert "ACGTACGT" in out
    assert "492 bp shorter" in out


def test_a_bad_deletion_fails_with_a_message_not_a_traceback(locus, capsys):
    code = main(["--sequence", str(locus), "--delete", "the second exon"])
    assert code == 1
    assert "error:" in capsys.readouterr().err


def test_no_sequence_and_no_species_is_explained(capsys):
    code = main(["--delete", "4001-4500"])
    assert code == 1
    assert "--sequence" in capsys.readouterr().err


def test_a_genomic_deletion_without_a_sequence_needs_a_chromosome(capsys):
    code = main(["--species", "mouse", "--delete", "4001-4500"])
    assert code == 1
    assert "must name a chromosome" in capsys.readouterr().err


def test_arm_zero_is_accepted_and_reported(locus, capsys):
    code = main(
        [
            "--sequence", str(locus), "--delete", "4001-4500",
            "--arm", "0", "--specificity", "local",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "0 bp homology arm on each side" in out
    assert "--arm" in out  # the note telling them to set it if a donor is used


def test_an_impossible_tm_window_exits_two_with_advice(locus, capsys):
    code = main(
        [
            "--sequence", str(locus), "--delete", "4001-4500",
            "--tm-range", "89,90", "--specificity", "local",
        ]
    )
    assert code == 2
    assert "No primer pair was found" in capsys.readouterr().out
