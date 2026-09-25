"""Design constants and the run-wide configuration object.

The defaults describe one PCR read out on a gel: both primers sit outside the
deleted region, so the untargeted allele gives a long product, the deleted
allele a shorter one, and a heterozygote both bands.

Two numbers govern where a primer is allowed to sit, and they are the reason
this tool exists rather than a generic primer picker:

``homology_arm``
    A repair template - an ssODN donor, or the arms of a targeting plasmid -
    carries the sequence immediately flanking the cut. A primer that overlaps
    that stretch anneals to the donor as readily as to the chromosome, so
    leftover transfection reagent amplifies and every clone reads as edited.
    Primers are therefore excluded from ``homology_arm`` bases on each side of
    the deletion. Set it to 0 when there is no donor.

``min_clearance``
    Beyond the arm, the junction itself is unreliable: Cas9 cutting leaves
    indels, and resection can eat further than intended. A primer that anneals
    within ``min_clearance`` of the excluded zone is usable but flagged, since
    a junction indel can shift its position or destroy the site outright.
"""

from __future__ import annotations

from dataclasses import dataclass

# --- Where a primer may sit --------------------------------------------------

#: Bases on each side of the deletion that primers must avoid, i.e. the length
#: of the repair template's homology arm. 43 matches a 120 nt ssODN carrying a
#: 34 bp insert; 0 means no donor is involved.
DEFAULT_HOMOLOGY_ARM = 43

#: Below this distance from the excluded zone a primer is reported with a
#: warning rather than silently accepted.
DEFAULT_MIN_CLEARANCE = 20

#: How far beyond each excluded zone the primer search reaches.
DEFAULT_SEARCH_SPAN = 1_000

# --- Primer properties -------------------------------------------------------

#: Acceptable primer length, in bases.
PRIMER_LENGTH_RANGE = (20, 30)

#: Target Tm and the window around it, in °C (nearest-neighbour, SantaLucia 1998).
PRIMER_TM_TARGET = 60.0
PRIMER_TM_RANGE = (58.0, 62.0)

#: Primers are rejected outside this GC range.
PRIMER_GC_RANGE = (0.40, 0.60)

#: The two primers of a pair must anneal within this many °C of each other.
PRIMER_MAX_TM_DIFF = 3.0

#: Runs of this many identical bases misprime on repeats; reject them.
PRIMER_MAX_HOMOPOLYMER = 5

#: Longest tolerated 3'-end base-pairing between the two primers of a pair.
PRIMER_MAX_DIMER = 4

#: How many ranked candidates per side are put through the specificity check.
PRIMER_CANDIDATES_PER_SIDE = 10

# --- Reading the gel ---------------------------------------------------------

#: Two bands closer than this will not be told apart on a normal agarose gel.
MIN_BAND_DIFFERENCE = 50

#: A product shorter than this runs off the front of a gel and is easily missed.
MIN_PRODUCT = 100

#: Above this, a product amplifies slowly enough that a failed reaction and a
#: deleted allele start to look alike. A warning, not a filter.
MAX_PRODUCT_WARN = 3_000

# --- Specificity -------------------------------------------------------------

#: ``auto``      - genome FASTA if given, else GGGenome if a database resolves,
#:                 else the supplied sequence alone.
#: ``gggenome``  - genome-wide near-match counts from the DBCLS service.
#: ``fasta``     - exact-match counts over a local genome FASTA.
#: ``local``     - near-match counts within the supplied sequence only. This is
#:                 the right check for a plasmid, where genome-wide uniqueness
#:                 is not the question being asked.
#: ``none``      - skip the check entirely.
SPECIFICITY_BACKENDS = ("auto", "gggenome", "fasta", "local", "none")

#: Mismatches allowed when counting off-target primer sites.
PRIMER_MISMATCHES = 2

#: Near-matches a primer may have, counting its own site. 1 means "nothing else
#: within the mismatch budget".
PRIMER_MAX_HITS = 1

#: 3'-terminal stretch counted separately by the FASTA backend. A primer whose
#: 3' end is common in the genome misprimes even when its full length is unique.
PRIMER_3PRIME_SEED = 15

#: Seconds between GGGenome requests; the service is a free public one.
GGGENOME_DELAY = 1.0

# --- Sequence retrieval ------------------------------------------------------

ENSEMBL_REST_URL = "https://rest.ensembl.org"
DEFAULT_SPECIES = "mus_musculus"


@dataclass(frozen=True)
class DesignConfig:
    """Everything a run needs beyond the target sequence and the deletion."""

    # where a primer may sit
    homology_arm: int = DEFAULT_HOMOLOGY_ARM
    min_clearance: int = DEFAULT_MIN_CLEARANCE
    search_span: int = DEFAULT_SEARCH_SPAN

    # primer properties
    length_range: tuple[int, int] = PRIMER_LENGTH_RANGE
    tm_target: float = PRIMER_TM_TARGET
    tm_range: tuple[float, float] = PRIMER_TM_RANGE
    gc_range: tuple[float, float] = PRIMER_GC_RANGE
    max_tm_diff: float = PRIMER_MAX_TM_DIFF
    candidates_per_side: int = PRIMER_CANDIDATES_PER_SIDE

    # reading the gel
    min_band_difference: int = MIN_BAND_DIFFERENCE
    max_product_warn: int = MAX_PRODUCT_WARN

    # specificity
    specificity: str = "auto"
    species: str = ""
    gggenome_db: str | None = None
    genome_fasta: str | None = None
    mismatches: int = PRIMER_MISMATCHES
    max_hits: int = PRIMER_MAX_HITS
    gggenome_delay: float = GGGENOME_DELAY

    # reporting
    alternatives: int = 2

    def __post_init__(self) -> None:
        if self.homology_arm < 0:
            raise ValueError("homology_arm must be >= 0")
        if self.min_clearance < 0:
            raise ValueError("min_clearance must be >= 0")
        if self.search_span < 1:
            raise ValueError("search_span must be >= 1")
        lo_len, hi_len = self.length_range
        if not 0 < lo_len <= hi_len:
            raise ValueError("length_range must satisfy 0 < low <= high")
        lo_tm, hi_tm = self.tm_range
        if not lo_tm < hi_tm:
            raise ValueError("tm_range must satisfy low < high")
        lo_gc, hi_gc = self.gc_range
        if not 0 <= lo_gc < hi_gc <= 1:
            raise ValueError("gc_range must satisfy 0 <= low < high <= 1")
        if self.specificity not in SPECIFICITY_BACKENDS:
            raise ValueError(f"specificity must be one of {list(SPECIFICITY_BACKENDS)}")
        if self.mismatches < 0:
            raise ValueError("mismatches must be >= 0")
        if self.alternatives < 0:
            raise ValueError("alternatives must be >= 0")


__all__ = [
    "DesignConfig",
    "GGGENOME_DELAY",
    "MAX_PRODUCT_WARN",
    "MIN_BAND_DIFFERENCE",
    "MIN_PRODUCT",
    "PRIMER_3PRIME_SEED",
    "PRIMER_MAX_DIMER",
    "PRIMER_MAX_HOMOPOLYMER",
    "SPECIFICITY_BACKENDS",
]
