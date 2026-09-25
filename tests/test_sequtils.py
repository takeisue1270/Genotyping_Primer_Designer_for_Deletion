"""Sequence arithmetic and thermodynamics."""

from __future__ import annotations

import pytest

from genotyping_primers.sequtils import (
    count_near_matches,
    find_all,
    gc_fraction,
    longest_homopolymer,
    melting_temp_nn,
    reverse_complement,
    three_prime_complementarity,
)


def test_reverse_complement_round_trips():
    sequence = "ACGTTGCAnnRYacgt"
    assert reverse_complement(reverse_complement(sequence)) == sequence
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AACCGGTT") == "AACCGGTT"


def test_nn_tm_matches_published_values():
    """SantaLucia NN Tm, against oligos whose Tm is widely tabulated."""
    # M13 forward (-20): vendors quote ~52 C at 250 nM / 50 mM Na+.
    assert melting_temp_nn("GTAAAACGACGGCCAGT") == pytest.approx(51.7, abs=1.5)
    assert melting_temp_nn("GCGCGCGCGCGCGCGCGCGC") > 75
    assert melting_temp_nn("ATATATATATATATATATAT") < 45


def test_nn_tm_undefined_for_ambiguous_or_tiny_sequences():
    assert melting_temp_nn("ACGTNACGTACGTACGTACG") == 0.0
    assert melting_temp_nn("A") == 0.0


def test_three_prime_complementarity_walks_inward_from_both_ends():
    # Reading inward from both 3' ends: C:G, G:C, T:A, then A:G stops the run.
    assert three_prime_complementarity("AAAAATGC", "TTTTGACG") == 3
    # A mismatch at the very 3' end stops the run before it starts.
    assert three_prime_complementarity("AAAAATGC", "TTTTTGCA") == 0


def test_longest_homopolymer():
    assert longest_homopolymer("") == 0
    assert longest_homopolymer("ACGT") == 1
    assert longest_homopolymer("ACGAAAAT") == 4


def test_gc_fraction_ignores_ambiguity():
    assert gc_fraction("GCAT") == 0.5
    assert gc_fraction("GCNN") == 1.0
    assert gc_fraction("") == 0.0


def test_find_all_counts_overlapping_hits():
    assert find_all("ATATAT", "ATAT") == [0, 2]
    assert find_all("ACGT", "TTT") == []


def test_count_near_matches_counts_both_strands_and_its_own_site():
    haystack = "TTTT" + "ACGTACGTAGGTCAGT" + "TTTT"
    needle = "ACGTACGTAGGTCAGT"
    exact, near = count_near_matches(haystack, needle, 0)
    assert (exact, near) == (1, 1)

    # The reverse complement of the same site is the same site, counted once.
    both = haystack + "GG" + reverse_complement(needle)
    exact, near = count_near_matches(both, needle, 0)
    assert exact == 2


def test_count_near_matches_finds_mismatched_sites():
    site = "ACGTACGTAGGTCAGTACGT"
    mutated = site[:10] + ("A" if site[10] != "A" else "C") + site[11:]
    haystack = "TTTTT" + site + "TTTTT" + mutated + "TTTTT"
    exact, near = count_near_matches(haystack, site, 1)
    assert exact == 1
    assert near == 2
