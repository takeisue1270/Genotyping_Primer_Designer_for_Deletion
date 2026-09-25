"""Turning what the user typed into offsets, and back again for the report."""

from __future__ import annotations

import pytest

from genotyping_primers.design import DesignError, edited_sequence, prepare, resolve_span
from genotyping_primers.ensembl import parse_locus


@pytest.mark.parametrize(
    "spec,expected",
    [
        ("101-600", ("", 101, 600)),
        ("101..600", ("", 101, 600)),
        ("101+500", ("", 101, 600)),
        ("11:69,000,000-69,001,000", ("11", 69_000_000, 69_001_000)),
        ("chr11:69000000..69001000", ("11", 69_000_000, 69_001_000)),
    ],
)
def test_parse_locus_accepts_the_forms_people_actually_type(spec, expected):
    assert parse_locus(spec) == expected


@pytest.mark.parametrize("spec", ["", "0-100", "the second exon", "11:abc-def", "101"])
def test_parse_locus_rejects_nonsense(spec):
    with pytest.raises(ValueError):
        parse_locus(spec)


def test_local_coordinates_are_one_based_inclusive(target):
    start, end = resolve_span(target, "101-600")
    assert (start, end) == (100, 600)  # 0-based half-open
    assert end - start == 500


def test_genomic_coordinates_are_offset_by_the_window(genomic_target):
    start, end = resolve_span(genomic_target, "11:69,001,001-69,001,500")
    assert (start, end) == (1_000, 1_500)
    # And the report converts them straight back.
    assert genomic_target.describe_span(start, end) == "11:69,001,001-69,001,500"


def test_local_coordinates_can_be_forced_on_a_positioned_target(genomic_target):
    # Small numbers are read as offsets by 'auto' anyway, but saying so is allowed.
    assert resolve_span(genomic_target, "1001-1500", "local") == (1_000, 1_500)
    with pytest.raises(DesignError, match="outside the fetched window"):
        resolve_span(genomic_target, "1001-1500", "genomic")


def test_a_region_name_must_match_the_target(genomic_target):
    with pytest.raises(DesignError, match="is from"):
        resolve_span(genomic_target, "7:69,001,001-69,001,500")


def test_a_region_name_needs_a_positioned_target(target):
    with pytest.raises(DesignError, match="no genomic coordinates"):
        resolve_span(target, "11:69,001,001-69,001,500")


def test_coordinates_outside_the_fetched_window_are_refused(genomic_target):
    with pytest.raises(DesignError, match="outside the fetched window"):
        resolve_span(genomic_target, "11:70,000,000-70,000,100")


def test_a_deletion_past_the_end_is_refused(target):
    with pytest.raises(DesignError, match="past the end"):
        resolve_span(target, "7900-9000")


def test_a_reversed_span_is_only_allowed_on_a_circle(target, plasmid):
    with pytest.raises(DesignError, match="circular"):
        resolve_span(target, "600-101")
    assert resolve_span(plasmid, "3900-200") == (3_899, 200)


def test_linear_targets_are_not_rotated(target):
    working, deletion = prepare(target, 100, 600)
    assert working is target
    assert (deletion.start, deletion.end, deletion.length) == (100, 600, 500)


def test_a_circular_deletion_over_the_origin_is_rotated_into_the_middle(plasmid):
    start, end = resolve_span(plasmid, "3900-200")
    working, deletion = prepare(plasmid, start, end)

    assert deletion.length == 301  # 3900..4000 is 101 bases, plus 1..200
    assert working.sequence[deletion.start : deletion.end] == (
        plasmid.sequence[3_899:] + plasmid.sequence[:200]
    )
    # Rotating must not lose or duplicate a base, and the report has to undo it.
    assert sorted(working.sequence) == sorted(plasmid.sequence)
    assert working.describe_span(deletion.start, deletion.end) == "3,900-200"


def test_rotation_leaves_room_on_both_sides(plasmid):
    _, deletion = prepare(plasmid, *resolve_span(plasmid, "3900-200"))
    assert deletion.start > 1_500
    assert len(plasmid.sequence) - deletion.end > 1_500


def test_editing_removes_exactly_the_named_block(target):
    working, deletion = prepare(target, *resolve_span(target, "1001-1500"))
    edited = edited_sequence(working, deletion)
    assert len(edited) == len(working) - 500
    assert edited == working.sequence[:1_000] + working.sequence[1_500:]
