"""Where primers are allowed to sit, and what the two bands measure.

The point of these tests is that nothing is taken on trust from the tool's own
report: every amplicon is re-measured from the sequences it claims to amplify.
"""

from __future__ import annotations

import dataclasses

import pytest

from genotyping_primers.design import design, edited_sequence, prepare, resolve_span
from genotyping_primers.models import Deletion
from genotyping_primers.primers import (
    PrimerError,
    collect_candidates,
    excluded_zone,
    rank_pairs,
    search_zones,
)
from genotyping_primers.sequtils import find_all, reverse_complement
from genotyping_primers.specificity import LocalChecker


@pytest.fixture
def deletion(target) -> Deletion:
    _, deletion = prepare(target, *resolve_span(target, "4001-4500"))
    return deletion


@pytest.fixture
def result(target, deletion, config):
    return design(target, deletion, config, checker=LocalChecker(target.sequence, 2))


# --- the excluded zone -------------------------------------------------------


def test_the_excluded_zone_is_the_deletion_plus_both_arms(target, deletion, config):
    lo, hi = excluded_zone(target, deletion, config)
    assert lo == deletion.start - config.homology_arm
    assert hi == deletion.end + config.homology_arm


def test_search_zones_sit_outside_the_arms(target, deletion, config):
    (f_lo, f_hi), (r_lo, r_hi) = search_zones(target, deletion, config)
    assert f_hi == deletion.start - config.homology_arm
    assert f_lo == f_hi - config.search_span
    assert r_lo == deletion.end + config.homology_arm
    assert r_hi == r_lo + config.search_span


def test_search_zones_are_clamped_to_the_sequence(target, config):
    deletion = Deletion(start=200, end=700)
    (f_lo, f_hi), _ = search_zones(target, deletion, config)
    assert f_lo == 0
    assert f_hi == 200 - config.homology_arm


def test_no_room_is_an_error_that_says_what_to_do(target, config):
    deletion = Deletion(start=10, end=500)
    with pytest.raises(PrimerError, match="lower --arm"):
        search_zones(target, deletion, config)


# --- the pairs ---------------------------------------------------------------


def test_a_pair_is_found_for_a_plain_deletion(result):
    assert result.pairs, result.notes


def test_no_primer_overlaps_the_excluded_zone(result, config):
    lo, hi = excluded_zone(result.target, result.deletion, config)
    for pair in result.pairs:
        assert pair.forward.end <= lo
        assert pair.reverse.start >= hi


def test_clearance_is_the_distance_to_the_excluded_zone(result, config):
    lo, hi = excluded_zone(result.target, result.deletion, config)
    for pair in result.pairs:
        assert pair.forward.clearance == lo - pair.forward.end
        assert pair.reverse.clearance == pair.reverse.start - hi
        assert pair.min_clearance >= 0


def test_every_primer_anneals_exactly_once_in_both_alleles(result):
    for pair in result.pairs:
        for primer in (pair.forward, pair.reverse):
            site = primer.top_strand()
            assert len(find_all(result.target.sequence, site)) == 1
            assert len(find_all(result.edited_sequence, site)) == 1


def test_band_sizes_are_what_the_sequences_actually_give(result):
    for pair in result.pairs:
        for sequence, expected in (
            (result.target.sequence, pair.wt_product),
            (result.edited_sequence, pair.ko_product),
        ):
            start = sequence.index(pair.forward.top_strand())
            end = sequence.index(pair.reverse.top_strand()) + pair.reverse.length
            assert end - start == expected


def test_the_size_shift_is_the_edit_itself(result):
    for pair in result.pairs:
        assert pair.size_shift == result.deletion.net_change == 500


def test_the_reverse_primer_is_ordered_as_its_reverse_complement(result):
    for pair in result.pairs:
        site = result.target.sequence[pair.reverse.start : pair.reverse.end]
        assert pair.reverse.sequence == reverse_complement(site)
        assert pair.reverse.strand == -1


def test_an_insert_shortens_the_shift_by_its_own_length(target, config):
    working, deletion = prepare(target, *resolve_span(target, "4001-4500"))
    deletion = dataclasses.replace(deletion, insert="A" * 34, insert_name="insert")
    result = design(working, deletion, config, checker=LocalChecker(target.sequence, 2))
    assert result.pairs
    for pair in result.pairs:
        assert pair.size_shift == 500 - 34
        assert pair.wt_product - pair.ko_product == 466


# --- warnings ----------------------------------------------------------------


def test_a_primer_inside_the_asked_for_clearance_is_flagged_not_dropped(target, deletion):
    from genotyping_primers.config import DesignConfig

    # Only 120 bp to search in, but 100 bp of clearance demanded: every
    # candidate is too close, and every pair must say so rather than vanish.
    cramped = DesignConfig(
        specificity="local", search_span=120, min_clearance=100, alternatives=0
    )
    result = design(target, deletion, cramped, checker=LocalChecker(target.sequence, 2))
    assert result.pairs
    warnings = " ".join(result.best.warnings)
    assert "homology arm" in warnings
    assert "junction indel" in warnings


def test_no_arm_declared_says_so_in_the_warning(target, deletion):
    from genotyping_primers.config import DesignConfig

    cramped = DesignConfig(
        specificity="local", homology_arm=0, search_span=120, min_clearance=100,
        alternatives=0,
    )
    result = design(target, deletion, cramped, checker=LocalChecker(target.sequence, 2))
    assert result.pairs
    assert "--arm" in " ".join(result.best.warnings)


def test_an_unchecked_primer_is_reported_as_unknown_not_clean(target, deletion, config):
    from genotyping_primers.specificity import NullChecker

    result = design(target, deletion, config, checker=NullChecker())
    assert result.pairs
    assert "specificity is unknown" in " ".join(result.best.warnings)


def test_alternatives_are_different_designs(target, deletion, config):
    result = design(target, deletion, config, checker=LocalChecker(target.sequence, 2))
    starts = {(p.forward.start, p.reverse.start) for p in result.pairs}
    assert len(starts) == len(result.pairs)


def test_a_shared_primer_keeps_one_name(target, deletion, config):
    result = design(target, deletion, config, checker=LocalChecker(target.sequence, 2))
    by_sequence: dict[str, set[str]] = {}
    for pair in result.pairs:
        for primer in (pair.forward, pair.reverse):
            by_sequence.setdefault(primer.sequence, set()).add(primer.name)
    assert all(len(names) == 1 for names in by_sequence.values())


def test_no_pair_at_all_is_reported_with_a_way_forward(target, deletion):
    from genotyping_primers.config import DesignConfig

    impossible = DesignConfig(specificity="local", tm_range=(89.0, 90.0))
    candidates = collect_candidates(target, deletion, impossible)
    pairs, notes = rank_pairs(
        target, deletion, edited_sequence(target, deletion), impossible, candidates
    )
    assert pairs == []
    assert any("--tm-range" in note or "--gc-range" in note for note in notes)
