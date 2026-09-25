"""Small sequence helpers.

Everything here works on plain ``str``. The nearest-neighbour Tm is the one
part that has to be right rather than approximately right, since primers are
selected *by* their Tm.
"""

from __future__ import annotations

import math

_COMPLEMENT = str.maketrans(
    "ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb"
)


def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA string, preserving case."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    """GC content as a fraction of unambiguous bases; 0.0 for an empty string."""
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    if acgt == 0:
        return 0.0
    return (seq.count("G") + seq.count("C")) / acgt


def is_unambiguous(seq: str) -> bool:
    """True when every base is A/C/G/T (assembly gaps show up as N)."""
    return bool(seq) and set(seq.upper()) <= set("ACGT")


def find_all(haystack: str, needle: str) -> list[int]:
    """All 0-based start offsets of ``needle`` in ``haystack``, overlaps included."""
    if not needle:
        return []
    hits: list[int] = []
    start = haystack.find(needle)
    while start != -1:
        hits.append(start)
        start = haystack.find(needle, start + 1)
    return hits


def count_occurrences(haystack: str, needle: str, both_strands: bool = True) -> int:
    """Occurrences of ``needle`` in ``haystack``, optionally counting both strands."""
    n = len(find_all(haystack, needle))
    if both_strands:
        rc = reverse_complement(needle)
        if rc != needle:
            n += len(find_all(haystack, rc))
    return n


def longest_homopolymer(seq: str) -> int:
    """Length of the longest run of one base; 0 for an empty string."""
    best = run = 0
    previous = ""
    for base in seq.upper():
        run = run + 1 if base == previous else 1
        previous = base
        best = max(best, run)
    return best


_WATSON_CRICK = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}


def three_prime_complementarity(left: str, right: str, window: int = 6) -> int:
    """How many bases pair when the two 3' ends are laid against each other.

    Two primers whose 3' ends base-pair extend on one another and burn the
    reaction as primer-dimer. In an antiparallel duplex the terminal base of
    one primer faces the terminal base of the other, so this walks inwards from
    both 3' ends and returns the length of the uninterrupted run. 0 means the
    ends do not pair at all.
    """
    a = left.upper()[-window:][::-1]  # read 3'->5'
    b = right.upper()[-window:][::-1]
    length = 0
    for base_a, base_b in zip(a, b):
        if (base_a, base_b) not in _WATSON_CRICK:
            break
        length += 1
    return length


# --- Nearest-neighbour Tm ----------------------------------------------------
#
# SantaLucia (1998) unified parameters, PNAS 95:1460. dH in kcal/mol, dS in
# cal/(mol.K).

_NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

#: Helix initiation, applied once per duplex end according to that end's base.
_INIT_GC = (0.1, -2.8)
_INIT_AT = (2.3, 4.1)

_GAS_CONSTANT = 1.987  # cal/(mol.K)

#: Defaults chosen to match a typical genotyping PCR, and the numbers most
#: oligo vendors quote: 50 mM monovalent salt, 250 nM primer.
DEFAULT_NA_MM = 50.0
DEFAULT_PRIMER_NM = 250.0


def melting_temp_nn(
    seq: str,
    primer_nm: float = DEFAULT_PRIMER_NM,
    na_mm: float = DEFAULT_NA_MM,
) -> float:
    """Nearest-neighbour Tm in °C for a primer annealing to its perfect complement.

    Uses the SantaLucia (1998) unified parameter set with that paper's
    monovalent-salt entropy correction. Returns 0.0 for a sequence shorter than
    two bases or one carrying anything but A/C/G/T, since the model is not
    defined there.
    """
    seq = seq.upper()
    if len(seq) < 2 or not is_unambiguous(seq):
        return 0.0

    dh, ds = 0.0, 0.0
    for i in range(len(seq) - 1):
        step = _NN_PARAMS.get(seq[i : i + 2])
        if step is None:  # pragma: no cover - guarded by is_unambiguous
            return 0.0
        dh += step[0]
        ds += step[1]

    for end in (seq[0], seq[-1]):
        init_h, init_s = _INIT_GC if end in "GC" else _INIT_AT
        dh += init_h
        ds += init_s

    # Monovalent salt correction on the entropy term.
    ds += 0.368 * (len(seq) - 1) * math.log(na_mm / 1000.0)

    # Non-self-complementary duplex: the effective concentration is CT/4.
    ct = (primer_nm * 1e-9) / 4.0
    denominator = ds + _GAS_CONSTANT * math.log(ct)
    if denominator == 0:  # pragma: no cover - defensive
        return 0.0
    return (dh * 1000.0) / denominator - 273.15


# --- Mismatch-tolerant counting ----------------------------------------------


def count_near_matches(haystack: str, needle: str, max_mismatches: int) -> tuple[int, int]:
    """``(exact, within max_mismatches)`` occurrences of ``needle``, both strands.

    Vectorised: the Hamming distance at every offset is accumulated one column
    at a time, so a 20 kb sequence and a 25 nt primer cost 25 passes over a 20k
    array rather than a Python loop over every position.
    """
    import numpy as np

    m, n = len(needle), len(haystack)
    if m == 0 or n < m:
        return 0, 0

    text = np.frombuffer(haystack.upper().encode("ascii", "replace"), dtype=np.uint8)
    positions = n - m + 1
    exact = near = 0

    for query in {needle.upper(), reverse_complement(needle.upper())}:
        pattern = np.frombuffer(query.encode("ascii", "replace"), dtype=np.uint8)
        mismatches = np.zeros(positions, dtype=np.int32)
        for offset in range(m):
            mismatches += text[offset : offset + positions] != pattern[offset]
        exact += int(np.count_nonzero(mismatches == 0))
        near += int(np.count_nonzero(mismatches <= max_mismatches))

    return exact, near


__all__ = [
    "count_near_matches",
    "count_occurrences",
    "find_all",
    "gc_fraction",
    "is_unambiguous",
    "longest_homopolymer",
    "melting_temp_nn",
    "reverse_complement",
    "three_prime_complementarity",
]
