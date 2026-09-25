"""The four ways of asking whether a primer is unique."""

from __future__ import annotations

import pytest

from genotyping_primers.specificity import (
    FastaChecker,
    GGGenomeChecker,
    GGGenomeError,
    LocalChecker,
    NullChecker,
    normalise_species,
    resolve_database,
)


# --- local -------------------------------------------------------------------


def test_local_checker_counts_a_primers_own_site_once(sequence):
    site = sequence[1_000:1_025]
    hits = LocalChecker(sequence, mismatches=2).counts([site])[site]
    assert hits.perfect == 1
    assert hits.is_unique
    assert "local" in hits.describe()


def test_local_checker_sees_a_duplicated_site(sequence):
    site = sequence[1_000:1_025]
    doubled = sequence + site
    hits = LocalChecker(doubled, mismatches=0).counts([site])[site]
    assert hits.perfect == 2
    assert not hits.is_unique


def test_a_circular_checker_sees_the_site_that_spans_the_origin(sequence):
    # A site made of the last 12 bases and the first 13 exists only on a circle.
    junction = sequence[-12:] + sequence[:13]
    linear = LocalChecker(sequence, mismatches=0, circular=False).counts([junction])
    circular = LocalChecker(sequence, mismatches=0, circular=True).counts([junction])
    assert linear[junction].perfect == 0
    assert circular[junction].perfect == 1


def test_null_checker_reports_unknown_not_clean():
    hits = NullChecker().counts(["ACGTACGTACGTACGTACGT"])["ACGTACGTACGTACGTACGT"]
    assert not hits.checked
    assert not hits.is_unique
    assert "not checked" in hits.describe()


# --- a local genome FASTA ----------------------------------------------------


def _write_fasta(path, records):
    path.write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in records), encoding="utf-8"
    )
    return path


def test_fasta_checker_counts_both_strands_across_records(tmp_path, sequence):
    from genotyping_primers.sequtils import reverse_complement

    site = sequence[2_000:2_025]
    fasta = _write_fasta(
        tmp_path / "genome.fa",
        [("chr1", sequence), ("chr2", "TTTT" + reverse_complement(site) + "TTTT")],
    )
    hits = FastaChecker(fasta, seed=15).counts([site])[site]
    assert hits.perfect == 2
    assert hits.three_prime is not None and hits.three_prime >= 2
    assert "exact matches counted genome-wide" in FastaChecker(fasta).describe()


def test_fasta_checker_counts_a_site_spanning_a_chunk_boundary_once(tmp_path, sequence):
    site = sequence[3_000:3_025]
    fasta = _write_fasta(tmp_path / "genome.fa", [("chr1", sequence)])
    # A chunk size that forces the site to straddle a boundary must not make it
    # vanish, and the overlap must not make it count twice.
    for chunk in (500, 1_000, 3_010):
        hits = FastaChecker(fasta, seed=15, chunk_size=chunk).counts([site])[site]
        assert hits.perfect == 1, chunk


# --- GGGenome ----------------------------------------------------------------


class StubbedGGGenome(GGGenomeChecker):
    """A checker whose HTTP layer is replaced by canned payloads."""

    def __init__(self, payloads, **kwargs):
        super().__init__(delay=0.0, **kwargs)
        self.payloads = payloads
        self.requests: list[tuple[str, int]] = []

    def _get(self, sequence, k):
        self.requests.append((sequence, k))
        payload = self.payloads[k]
        self._validate(payload)
        return payload


def _payload(count, database="Mouse genome, GRCm39/mm39", approx=False):
    return {
        "error": "none",
        "database": database,
        "summary": [{"count": count, "count_is_approx": approx}],
    }


def test_a_unique_primer_scores_one_exact_and_one_near():
    checker = StubbedGGGenome({2: _payload(1), 0: _payload(1)}, database="mm39", expect="Mouse")
    hits = checker.counts(["ACGTACGTACGTACGTACGT"])["ACGTACGTACGTACGTACGT"]
    assert (hits.perfect, hits.near) == (1, 1)
    assert hits.is_unique
    assert checker.requests == [("ACGTACGTACGTACGTACGT", 2), ("ACGTACGTACGTACGTACGT", 0)]


def test_a_primer_with_off_targets_skips_the_second_request():
    checker = StubbedGGGenome({2: _payload(7), 0: _payload(1)}, database="mm39", expect="Mouse")
    hits = checker.counts(["ACGTACGTACGTACGTACGT"])["ACGTACGTACGTACGTACGT"]
    assert hits.near == 7
    assert not hits.is_unique
    assert len(checker.requests) == 1  # the exact count cannot rescue it


def test_the_human_fallback_is_refused_rather_than_believed():
    """GGGenome answers an unknown database with hg38 instead of an error."""
    checker = StubbedGGGenome(
        {2: _payload(1, database="Human genome, GRCh38/hg38")},
        database="mm39_typo",
        expect="Mouse",
    )
    hits = checker.counts(["ACGTACGTACGTACGTACGT"])["ACGTACGTACGTACGTACGT"]
    assert not hits.checked
    assert "hg38" in hits.error
    assert "not checked" in hits.describe()


def test_an_error_payload_is_not_a_count():
    checker = StubbedGGGenome({2: {"error": "too many hits"}}, database="mm39")
    hits = checker.counts(["ACGTACGTACGTACGTACGT"])["ACGTACGTACGTACGTACGT"]
    assert not hits.checked
    assert "too many hits" in hits.error


def test_an_ambiguous_primer_is_never_sent():
    checker = StubbedGGGenome({}, database="mm39")
    hits = checker.counts(["ACGTNACGTACGTACGTACG"])["ACGTNACGTACGTACGTACG"]
    assert not hits.checked
    assert checker.requests == []


def test_validate_rejects_the_wrong_organism_directly():
    checker = GGGenomeChecker(database="mm39", expect="Mouse")
    with pytest.raises(GGGenomeError, match="not the expected Mouse genome"):
        checker._validate({"database": "Human genome, GRCh38/hg38"})


# --- choosing a database -----------------------------------------------------


@pytest.mark.parametrize(
    "species,expected",
    [
        ("mus_musculus", ("mm39", "Mouse")),
        ("mouse", ("mm39", "Mouse")),
        ("Mus musculus", ("mm39", "Mouse")),
        ("human", ("hg38", "Human")),
        ("felis_catus", None),
    ],
)
def test_resolve_database_accepts_the_names_people_use(species, expected):
    assert resolve_database(species, None) == expected


def test_an_explicit_database_is_taken_at_face_value():
    # No organism expectation: the user asked for this one by name.
    assert resolve_database("mus_musculus", "mm10") == ("mm10", "")
    # ...unless it names a species, which is looked up like any other.
    assert resolve_database("", "mouse") == ("mm39", "Mouse")


def test_species_names_are_normalised():
    assert normalise_species("Mus musculus") == "mus_musculus"
    assert normalise_species("MOUSE") == "mus_musculus"
    assert normalise_species("felis_catus") == "felis_catus"
