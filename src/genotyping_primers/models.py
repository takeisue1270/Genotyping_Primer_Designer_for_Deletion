"""The objects a design is built from.

Coordinates are the thing to keep straight. Internally every offset is 0-based
half-open into ``Target.sequence``, which is the *working* sequence: for a
circular target that may be a rotation of what the file held, so that a
deletion spanning the origin still has room on both sides. Everything a user
reads - the report, the TSV, the GenBank files - is converted back by
``Target.describe``/``Target.to_source`` into the numbering they supplied.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

from .sequtils import reverse_complement


@dataclass(frozen=True)
class Target:
    """The sequence a deletion is planned in.

    ``sequence`` is linear and ready to slice. ``rotation`` is the offset in
    the source sequence that working offset 0 corresponds to; it is 0 unless a
    circular target had to be rotated. ``region``/``region_start`` are set when
    the sequence came from a reference genome, and turn offsets into genomic
    coordinates in the report.
    """

    name: str
    sequence: str
    circular: bool = False
    rotation: int = 0
    source: str = ""
    species: str = ""
    assembly: str = ""
    region: str = ""
    region_start: int = 1  # 1-based coordinate of source offset 0

    def __len__(self) -> int:
        return len(self.sequence)

    def to_source(self, offset: int) -> int:
        """Working offset -> 1-based position in the sequence as supplied."""
        return (offset + self.rotation) % len(self.sequence) + 1

    def to_genomic(self, offset: int) -> int:
        """Working offset -> 1-based genomic coordinate (only if ``region`` is set)."""
        return self.region_start + self.to_source(offset) - 1

    def describe(self, offset: int) -> str:
        """The position as a user should read it, in their own coordinate system."""
        if self.region:
            return f"{self.region}:{self.to_genomic(offset):,}"
        return f"{self.to_source(offset):,}"

    def describe_span(self, start: int, end: int) -> str:
        """A 0-based half-open working span, written as a 1-based inclusive one."""
        if self.region:
            return (
                f"{self.region}:{self.to_genomic(start):,}-{self.to_genomic(end - 1):,}"
            )
        return f"{self.to_source(start):,}-{self.to_source(end - 1):,}"

    @property
    def coordinate_system(self) -> str:
        if self.region:
            return f"{self.assembly or self.species or 'genome'} {self.region}"
        return f"{self.name} (1-based)"


@dataclass(frozen=True)
class Deletion:
    """The block to be removed, and whatever replaces it.

    ``start``/``end`` are 0-based half-open working offsets. ``insert`` is the
    sequence left behind - a tag, a restriction site, or nothing at all for a
    plain deletion.
    """

    start: int
    end: int
    insert: str = ""
    insert_name: str = ""

    @property
    def length(self) -> int:
        """Bases removed."""
        return self.end - self.start

    @property
    def insert_length(self) -> int:
        return len(self.insert)

    @property
    def net_change(self) -> int:
        """How much shorter the edited allele is. Negative if the insert is longer."""
        return self.length - self.insert_length


@dataclass(frozen=True)
class SpecificityHits:
    """Where a sequence occurs, counted on both strands including its own site."""

    source: str
    mismatches: int
    perfect: int | None  # exact matches; None when not determined
    near: int | None  # matches within ``mismatches``; None when not determined
    approximate: bool = False
    error: str = ""
    #: Exact occurrences of the primer's 3'-terminal stretch, when the backend
    #: counts it. A primer whose 3' end is common misprimes even if its full
    #: length is unique.
    three_prime: int | None = None
    three_prime_length: int = 0

    @property
    def checked(self) -> bool:
        return self.near is not None and not self.error

    @property
    def is_unique(self) -> bool:
        """One exact site and nothing else within the mismatch budget."""
        return self.checked and self.perfect == 1 and self.near == 1

    def describe(self) -> str:
        if self.error:
            return f"{self.source}: not checked ({self.error})"
        about = "~" if self.approximate else ""
        text = (
            f"{self.source}: {self.perfect} exact"
            if self.mismatches == 0
            else f"{self.source}: {self.perfect} exact, "
            f"{about}{self.near} within {self.mismatches} mismatch(es)"
        )
        if self.three_prime is not None:
            text += f", {self.three_prime} sites match its 3' {self.three_prime_length} nt"
        return text


@dataclass(frozen=True)
class Primer:
    """One primer, positioned in the unedited target.

    ``start``/``end`` are 0-based half-open working offsets of the annealing
    site regardless of strand; ``sequence`` is always written 5'->3' as it
    would be ordered, so a reverse primer is the reverse complement of
    ``target.sequence[start:end]``.

    ``clearance`` is how many bases separate the primer from the zone that is
    off limits (the deletion plus its homology arms). 0 means it abuts the arm.
    """

    name: str
    sequence: str
    start: int
    end: int
    strand: int  # +1 forward, -1 reverse
    tm: float
    gc: float
    clearance: int = 0
    specificity: SpecificityHits | None = None

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def strand_symbol(self) -> str:
        return "+" if self.strand >= 0 else "-"

    def binding_site(self, sequence: str) -> str:
        """The top-strand stretch this primer anneals to."""
        return sequence[self.start : self.end]

    def renamed(self, name: str) -> Primer:
        return dataclasses.replace(self, name=name)

    def with_specificity(self, hits: SpecificityHits | None) -> Primer:
        return dataclasses.replace(self, specificity=hits)

    def with_clearance(self, clearance: int) -> Primer:
        return dataclasses.replace(self, clearance=clearance)

    def is_specific(self, max_hits: int) -> bool:
        """True when the specificity check ran and the primer cleared it.

        An unchecked primer is not specific *or* non-specific - it is unknown,
        and reported as such rather than quietly counted as a pass.
        """
        hits = self.specificity
        if hits is None or not hits.checked:
            return False
        return hits.perfect == 1 and (hits.near or 0) <= max_hits

    def specificity_note(self) -> str:
        return self.specificity.describe() if self.specificity else "not checked"

    def ordered_sequence(self) -> str:
        """The 5'->3' sequence to put on a synthesis order."""
        return self.sequence

    def top_strand(self) -> str:
        """The same site read along the target, whichever strand the primer is on."""
        return self.sequence if self.strand >= 0 else reverse_complement(self.sequence)


@dataclass(frozen=True)
class PrimerPair:
    """A genotyping pair and the two band sizes it is read by."""

    forward: Primer
    reverse: Primer
    #: Amplicon length on the unedited allele, in bp.
    wt_product: int
    #: Amplicon length on the edited (deleted) allele, in bp.
    ko_product: int
    deletion_length: int
    insert_length: int
    #: Warnings that do not invalidate the pair but must be read before it is
    #: ordered: tight clearance, unresolvable bands, unchecked specificity.
    warnings: tuple[str, ...] = ()
    score: float = 0.0

    @property
    def size_shift(self) -> int:
        """How much shorter the edited band is; negative if it is longer."""
        return self.wt_product - self.ko_product

    @property
    def tm_difference(self) -> float:
        return abs(self.forward.tm - self.reverse.tm)

    @property
    def min_clearance(self) -> int:
        return min(self.forward.clearance, self.reverse.clearance)

    def resolvable(self, min_difference: int) -> bool:
        """Whether a normal agarose gel would separate the two products."""
        if min(self.wt_product, self.ko_product) <= 0:
            return False
        return abs(self.size_shift) >= min_difference

    def is_specific(self, max_hits: int) -> bool:
        return self.forward.is_specific(max_hits) and self.reverse.is_specific(max_hits)


@dataclass(frozen=True)
class DesignResult:
    """Everything one run produced, ready to be reported or written out."""

    target: Target
    deletion: Deletion
    edited_sequence: str
    pairs: tuple[PrimerPair, ...]
    notes: tuple[str, ...] = ()
    specificity_backend: str = "none"

    @property
    def best(self) -> PrimerPair | None:
        return self.pairs[0] if self.pairs else None


__all__ = [
    "Deletion",
    "DesignResult",
    "Primer",
    "PrimerPair",
    "SpecificityHits",
    "Target",
]
