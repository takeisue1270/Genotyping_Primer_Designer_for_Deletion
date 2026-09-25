"""Choosing the pair, and sizing the two bands it produces.

One PCR across the edit separates all three genotypes: the unedited allele
gives the long product, the edited allele a shorter one, and a heterozygote
both bands. The two lengths are the deliverable, so they are measured from the
actual sequences rather than asserted.

Where a primer may sit
----------------------

::

    ... [ search span ][ arm ][   deleted region   ][ arm ][ search span ] ...
          ^ forward                                          reverse ^
                      |<- clearance ->|

The ``arm`` blocks are excluded outright. They are the sequence a repair
template carries, so a primer inside one anneals to the donor as readily as to
the chromosome: leftover transfection reagent would amplify and every clone
would read as edited. Beyond the arm, a primer that sits too close is accepted
but flagged, because the junction is where indels accumulate and a primer site
a few bases from a Cas9 cut can be damaged along with it.
"""

from __future__ import annotations

import dataclasses
import logging

from .config import (
    MIN_PRODUCT,
    PRIMER_MAX_DIMER,
    PRIMER_MAX_HOMOPOLYMER,
    DesignConfig,
)
from .models import Deletion, Primer, PrimerPair, SpecificityHits, Target
from .sequtils import (
    find_all,
    gc_fraction,
    is_unambiguous,
    longest_homopolymer,
    melting_temp_nn,
    reverse_complement,
    three_prime_complementarity,
)
from .specificity import SpecificityChecker

logger = logging.getLogger(__name__)

#: Above this many exact matches to its 3' end, a primer's business end is
#: sitting in a repeat even if its full length is unique.
THREE_PRIME_MAX_HITS = 50

#: Two alternatives whose primers both start within this many bases of an
#: already-listed pair are the same design, not a second opinion.
ALTERNATIVE_MIN_GAP = 15

#: Clearance a primer is ranked towards. Above ``min_clearance`` a primer is
#: acceptable, but a junction indel is more likely to reach a primer 30 bp away
#: than one 150 bp away, so distance is preferred where it is free.
COMFORTABLE_CLEARANCE = 100

#: Products longer than this amplify slowly enough to be worth trading a
#: little primer quality against.
LONG_PRODUCT = 1_500


class PrimerError(RuntimeError):
    """Raised when the sequence cannot accommodate a primer search at all."""


def excluded_zone(target: Target, deletion: Deletion, config: DesignConfig) -> tuple[int, int]:
    """The half-open span no primer may overlap: the deletion plus both arms."""
    lo = max(0, deletion.start - config.homology_arm)
    hi = min(len(target), deletion.end + config.homology_arm)
    return lo, hi


def search_zones(
    target: Target, deletion: Deletion, config: DesignConfig
) -> tuple[tuple[int, int], tuple[int, int]]:
    """The two windows outside the excluded zone, as 0-based half-open spans.

    Both are clamped to the sequence, so a deletion near the end of a supplied
    contig yields a shorter window rather than an error - until there is no
    room at all, which is an error worth raising.
    """
    lo, hi = excluded_zone(target, deletion, config)
    span = config.search_span

    forward = (max(0, lo - span), lo)
    reverse = (hi, min(len(target), hi + span))

    min_len = config.length_range[0]
    if forward[1] - forward[0] < min_len:
        raise PrimerError(
            f"{target.name}: only {forward[1] - forward[0]} bp available upstream of the "
            f"excluded zone ({config.homology_arm} bp arm); at least {min_len} bp are "
            "needed. Supply more flanking sequence, or lower --arm."
        )
    if reverse[1] - reverse[0] < min_len:
        raise PrimerError(
            f"{target.name}: only {reverse[1] - reverse[0]} bp available downstream of the "
            f"excluded zone ({config.homology_arm} bp arm); at least {min_len} bp are "
            "needed. Supply more flanking sequence, or lower --arm."
        )
    return forward, reverse


def _gc_clamp_ok(sequence: str) -> bool:
    """1-3 G/C in the last five bases: enough to seat the 3' end, not enough to stick."""
    tail = sequence[-5:].upper()
    return 1 <= tail.count("G") + tail.count("C") <= 3


def _quality(tm: float, gc: float, config: DesignConfig) -> float:
    """Higher is better. Tm proximity dominates; GC near 50 % breaks ties."""
    return -abs(tm - config.tm_target) - 2.0 * abs(gc - 0.5)


def enumerate_candidates(
    sequence: str,
    zone: tuple[int, int],
    strand: int,
    config: DesignConfig,
) -> list[Primer]:
    """Every acceptable primer in ``zone``, best first.

    ``strand`` is +1 for the forward primer (reads along the sequence) and -1
    for the reverse primer (the reverse complement is what gets ordered).
    ``start``/``end`` are offsets into ``sequence`` either way.
    """
    lo, hi = zone
    min_len, max_len = config.length_range
    tm_lo, tm_hi = config.tm_range
    gc_lo, gc_hi = config.gc_range

    candidates: list[tuple[float, Primer]] = []
    for start in range(lo, hi - min_len + 1):
        for length in range(min_len, max_len + 1):
            end = start + length
            if end > hi:
                break
            window = sequence[start:end].upper()
            if not is_unambiguous(window):
                continue
            ordered = window if strand >= 0 else reverse_complement(window)
            if longest_homopolymer(ordered) >= PRIMER_MAX_HOMOPOLYMER:
                continue
            gc = gc_fraction(ordered)
            if not gc_lo <= gc <= gc_hi:
                continue
            if not _gc_clamp_ok(ordered):
                continue
            tm = melting_temp_nn(ordered)
            if not tm_lo <= tm <= tm_hi:
                continue
            candidates.append(
                (
                    _quality(tm, gc, config),
                    Primer(
                        name="",
                        sequence=ordered,
                        start=start,
                        end=end,
                        strand=1 if strand >= 0 else -1,
                        tm=tm,
                        gc=gc,
                    ),
                )
            )

    candidates.sort(key=lambda item: item[1].start)
    candidates.sort(key=lambda item: item[0], reverse=True)
    return [primer for _, primer in candidates]


def _spread(primers: list[Primer], limit: int, min_gap: int = ALTERNATIVE_MIN_GAP) -> list[Primer]:
    """Trim to ``limit`` candidates that do not all start at the same base.

    Consecutive offsets give near-identical primers, so if the best one fails
    the specificity check its neighbours almost certainly will too. Spacing the
    shortlist out means each request buys an independent answer.
    """
    chosen: list[Primer] = []
    for primer in primers:
        if len(chosen) >= limit:
            break
        if all(abs(primer.start - other.start) >= min_gap for other in chosen):
            chosen.append(primer)
    for primer in primers:  # backfill if the zone was too cramped to space them
        if len(chosen) >= limit:
            break
        if primer not in chosen:
            chosen.append(primer)
    return chosen


def _with_clearance(primers: list[Primer], boundary: int, forward: bool) -> list[Primer]:
    """Record each primer's distance from the excluded zone."""
    return [
        primer.with_clearance(
            boundary - primer.end if forward else primer.start - boundary
        )
        for primer in primers
    ]


def collect_candidates(
    target: Target, deletion: Deletion, config: DesignConfig
) -> tuple[list[Primer], list[Primer]]:
    """The shortlist put forward for the specificity check, both sides."""
    forward_zone, reverse_zone = search_zones(target, deletion, config)
    lo, hi = excluded_zone(target, deletion, config)

    forward = _spread(
        enumerate_candidates(target.sequence, forward_zone, 1, config),
        config.candidates_per_side,
    )
    reverse = _spread(
        enumerate_candidates(target.sequence, reverse_zone, -1, config),
        config.candidates_per_side,
    )
    return _with_clearance(forward, lo, True), _with_clearance(reverse, hi, False)


def _locate(sequence: str, primer: Primer) -> int | None:
    """Where ``primer`` anneals in ``sequence``, or ``None`` if not exactly once.

    Two sites are as useless as none here: the amplicon length would be
    ambiguous, so the pair is rejected rather than sized from the first hit.
    """
    target = primer.top_strand()
    hits = find_all(sequence.upper(), target)
    return hits[0] if len(hits) == 1 else None


def _primer_warnings(primer: Primer, config: DesignConfig) -> list[str]:
    side = "forward" if primer.strand >= 0 else "reverse"
    warnings: list[str] = []

    if primer.clearance < config.min_clearance:
        if config.homology_arm:
            warnings.append(
                f"the {side} primer anneals {primer.clearance} bp from the "
                f"{config.homology_arm} bp homology arm (under the {config.min_clearance} bp "
                "asked for): it is clear of the repair template, but a junction indel or "
                "resection past the intended cut could still reach it"
            )
        else:
            warnings.append(
                f"the {side} primer anneals {primer.clearance} bp from the deleted region "
                f"(under the {config.min_clearance} bp asked for), and no homology arm was "
                "declared: if a repair template is used, check that this primer does not "
                "overlap it, or re-run with --arm set"
            )

    hits = primer.specificity
    if hits is None or not hits.checked:
        reason = hits.error if hits and hits.error else "not checked"
        warnings.append(f"the {side} primer's specificity is unknown ({reason})")
    elif not primer.is_specific(config.max_hits):
        warnings.append(f"the {side} primer is not unique - {hits.describe()}")
    elif hits.three_prime is not None and hits.three_prime > THREE_PRIME_MAX_HITS:
        warnings.append(
            f"the {side} primer is unique over its full length, but its 3' "
            f"{hits.three_prime_length} nt occur {hits.three_prime} times in the genome, "
            "which is where mispriming starts"
        )
    return warnings


def rank_pairs(
    target: Target,
    deletion: Deletion,
    edited: str,
    config: DesignConfig,
    candidates: tuple[list[Primer], list[Primer]],
    name: str = "",
) -> tuple[list[PrimerPair], list[str]]:
    """Every workable combination, best first, each carrying its own warnings.

    Nothing is filtered on the strength of a warning. A pair that sits close to
    the junction, or whose specificity could not be established, is ranked
    below a clean one but still offered, because on a cramped locus it may be
    the only thing available - and a warning that is read beats a design that
    silently did not appear.
    """
    forward, reverse = candidates
    notes: list[str] = []

    if not forward or not reverse:
        side = "upstream" if not forward else "downstream"
        notes.append(
            f"no primer passed the Tm/GC filters on the {side} side; "
            "widen --tm-range or --gc-range, or --span for more sequence to search"
        )
        return [], notes

    label = name or target.name
    scored: list[tuple[float, PrimerPair]] = []
    rejected_dimer = 0

    for fwd in forward:
        for rev in reverse:
            delta = abs(fwd.tm - rev.tm)
            if delta > config.max_tm_diff:
                continue
            if three_prime_complementarity(fwd.sequence, rev.sequence) > PRIMER_MAX_DIMER:
                rejected_dimer += 1
                continue

            # Size both products from the sequences themselves. The unedited
            # amplicon is read off the target; the edited one off the edited
            # allele, where both primers must still occur exactly once.
            wt_product = rev.end - fwd.start
            edited_forward = _locate(edited, fwd)
            edited_reverse = _locate(edited, rev)
            if edited_forward is None or edited_reverse is None:
                continue
            ko_product = edited_reverse + rev.length - edited_forward
            if ko_product != wt_product - deletion.net_change:
                # The two routes to the same number disagree, so neither is
                # trustworthy; drop the pair rather than report either.
                continue

            warnings = _primer_warnings(fwd, config) + _primer_warnings(rev, config)
            if ko_product < MIN_PRODUCT:
                warnings.append(
                    f"the edited band is only {ko_product} bp and may run off the gel; "
                    "increase --span to push the primers further out"
                )
            if wt_product > config.max_product_warn:
                warnings.append(
                    f"the unedited band is {wt_product:,} bp, long enough that a failed "
                    "reaction and a deleted allele can look alike; keep a positive control "
                    "on the gel, or lower --span"
                )

            # Ranking, in order of what it costs to get wrong: a primer that is
            # not unique, or that sits inside the clearance the user asked for,
            # is pushed below every clean pair. Among clean pairs, Tm and GC
            # decide, then distance from the junction, then amplicon length.
            score = (
                _quality(fwd.tm, fwd.gc, config)
                + _quality(rev.tm, rev.gc, config)
                - delta
                - 100.0 * (not fwd.is_specific(config.max_hits))
                - 100.0 * (not rev.is_specific(config.max_hits))
                - 50.0 * (fwd.clearance < config.min_clearance)
                - 50.0 * (rev.clearance < config.min_clearance)
                - 0.01 * max(0, COMFORTABLE_CLEARANCE - fwd.clearance)
                - 0.01 * max(0, COMFORTABLE_CLEARANCE - rev.clearance)
                - 0.001 * max(0, wt_product - LONG_PRODUCT)
            )
            scored.append(
                (
                    score,
                    PrimerPair(
                        forward=fwd,
                        reverse=rev,
                        wt_product=wt_product,
                        ko_product=ko_product,
                        deletion_length=deletion.length,
                        insert_length=deletion.insert_length,
                        warnings=tuple(warnings),
                        score=score,
                    ),
                )
            )

    if not scored:
        detail = (
            f" ({rejected_dimer} combination(s) were dropped for 3'-end dimerisation)"
            if rejected_dimer
            else ""
        )
        notes.append(
            f"no primer combination satisfied the {config.max_tm_diff} C Tm match "
            f"without dimerising{detail}; raise --max-tm-diff or widen --tm-range"
        )
        return [], notes

    scored.sort(key=lambda item: item[0], reverse=True)

    # Alternatives should be independent designs, not the same pair shifted a
    # base or two: at least one of the two primers has to move meaningfully.
    chosen: list[PrimerPair] = []
    for _, pair in scored:
        if any(
            abs(pair.forward.start - other.forward.start) < ALTERNATIVE_MIN_GAP
            and abs(pair.reverse.start - other.reverse.start) < ALTERNATIVE_MIN_GAP
            for other in chosen
        ):
            continue
        chosen.append(pair)
        if len(chosen) > config.alternatives:
            break

    # Name by sequence, not by position in the list. Neighbouring pairs often
    # share one primer, and giving the same oligo two names is how it ends up
    # on the order form twice.
    names: dict[str, str] = {}
    used_by: dict[str, list[int]] = {}
    for index, pair in enumerate(chosen, start=1):
        suffix = "" if index == 1 else str(index)
        renamed: dict[str, Primer] = {}
        for role, tag, primer in (
            ("forward", "F", pair.forward),
            ("reverse", "R", pair.reverse),
        ):
            name = names.get(primer.sequence)
            if name is None:
                name = f"{label}_{tag}{suffix}"
                names[primer.sequence] = name
            used_by.setdefault(name, []).append(index)
            renamed[role] = primer.renamed(name)
        chosen[index - 1] = dataclasses.replace(pair, **renamed)

    for name, pairs_using in used_by.items():
        if len(pairs_using) > 1:
            listed = ", ".join(str(i) for i in pairs_using)
            notes.append(f"{name} is the same oligo in pairs {listed}; order it once.")

    return chosen, notes


def annotate_specificity(
    candidates: tuple[list[Primer], list[Primer]], checker: SpecificityChecker
) -> tuple[list[Primer], list[Primer]]:
    """Run the check once for every candidate and attach the counts."""
    forward, reverse = candidates
    sequences = [primer.sequence for primer in forward + reverse]
    if not sequences:
        return forward, reverse
    logger.info(
        "checking %d primer candidate(s) against %s", len(set(sequences)), checker.source
    )
    hits: dict[str, SpecificityHits] = checker.counts(sequences)
    return (
        [p.with_specificity(hits.get(p.sequence.upper())) for p in forward],
        [p.with_specificity(hits.get(p.sequence.upper())) for p in reverse],
    )


__all__ = [
    "PrimerError",
    "THREE_PRIME_MAX_HITS",
    "annotate_specificity",
    "collect_candidates",
    "enumerate_candidates",
    "excluded_zone",
    "rank_pairs",
    "search_zones",
]
