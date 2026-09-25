"""How often a primer occurs somewhere it should not.

A primer will prime on a site that differs by a base or two, especially when
the mismatches sit away from the 3' end, so "is this 25-mer unique?" is not an
exact-match question. Four backends answer it, all through the same ``counts``
call so the caller does not care which one is in use:

``GGGenomeChecker``
    Counts sites within *k* mismatches anywhere in a reference assembly, using
    the DBCLS GGGenome service. This is the real check for a genome edit.

``FastaChecker``
    Exact matches over a local genome FASTA, plus a separate count for the
    primer's 3'-terminal stretch. Offline and fast, but exact-match only: it
    will not find the 1-2 mismatch sites GGGenome reports.

``LocalChecker``
    Mismatch-tolerant counts within the supplied sequence alone. For a plasmid
    this is not a fallback but the correct question - a primer used on a
    construct only has to be unique within that construct. For a genomic
    target it can show a primer is *not* unique, never that it is.

``NullChecker``
    Skips the check. Primers come back marked "not checked", never as clean.

Every backend counts both strands and counts the primer's own site, so a
perfectly specific primer scores exactly 1.
"""

from __future__ import annotations

import gzip
import json
import logging
import pathlib
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterator

from .models import SpecificityHits
from .sequtils import count_near_matches, find_all, is_unambiguous, reverse_complement

logger = logging.getLogger(__name__)

GGGENOME_URL = "https://gggenome.dbcls.jp"

#: Ensembl species -> (GGGenome database, a word that must appear in the name
#: the service reports back).
#:
#: The second element is not decoration. GGGenome answers an unknown database
#: name with the *human* genome instead of an error, so a typo or a species
#: guessed wrong would silently return hg38 hit counts for a mouse primer.
#: Every response is checked against the expected organism before its counts
#: are believed.
GGGENOME_DATABASES: dict[str, tuple[str, str]] = {
    "homo_sapiens": ("hg38", "Human"),
    "mus_musculus": ("mm39", "Mouse"),
    "rattus_norvegicus": ("rn6", "Rat"),
    "danio_rerio": ("danRer11", "Zebrafish"),
    "drosophila_melanogaster": ("dm6", "melanogaster"),
    "caenorhabditis_elegans": ("ce11", "elegans"),
    "gallus_gallus": ("galGal5", "Chicken"),
    "xenopus_tropicalis": ("xenTro9", "tropicalis"),
    "xenopus_laevis": ("xenLae2", "clawed frog"),
    "canis_lupus_familiaris": ("canFam3", "Dog"),
    "canis_familiaris": ("canFam3", "Dog"),
    "bos_taurus": ("bosTau8", "Cow"),
    "macaca_mulatta": ("rheMac8", "Rhesus"),
    "sus_scrofa": ("susScr11", "Pig"),
    "oryctolagus_cuniculus": ("oryCun2", "Rabbit"),
    "oryzias_latipes": ("oryLat2", "Medaka"),
    "ciona_intestinalis": ("ci3", "Ciona"),
    "saccharomyces_cerevisiae": ("sacCer3", "cerevisiae"),
    "arabidopsis_thaliana": ("TAIR10", "Arabidopsis"),
}

#: Common short names accepted wherever a species is expected.
SPECIES_ALIASES: dict[str, str] = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "zebrafish": "danio_rerio",
    "fly": "drosophila_melanogaster",
    "worm": "caenorhabditis_elegans",
    "yeast": "saccharomyces_cerevisiae",
    "chicken": "gallus_gallus",
    "pig": "sus_scrofa",
    "cow": "bos_taurus",
    "rabbit": "oryctolagus_cuniculus",
    "medaka": "oryzias_latipes",
}


def normalise_species(species: str) -> str:
    """Accept ``mouse`` or ``Mus musculus`` where Ensembl wants ``mus_musculus``."""
    key = species.strip().lower().replace(" ", "_")
    return SPECIES_ALIASES.get(key, key)


class SpecificityChecker:
    """Interface: map sequences to hit counts, caching repeats."""

    source = "none"

    def __init__(self) -> None:
        self._cache: dict[str, SpecificityHits] = {}

    def counts(self, sequences: list[str]) -> dict[str, SpecificityHits]:
        for sequence in dict.fromkeys(sequences):
            key = sequence.upper()
            if key not in self._cache:
                self._cache[key] = self._count(key)
        return {s.upper(): self._cache[s.upper()] for s in sequences}

    def _count(self, sequence: str) -> SpecificityHits:  # pragma: no cover - abstract
        raise NotImplementedError

    def describe(self) -> str:
        """One sentence for the report, saying what was and was not covered."""
        return "no specificity check was run"


class NullChecker(SpecificityChecker):
    """Skips the check; every primer comes back unchecked."""

    source = "none"

    def _count(self, sequence: str) -> SpecificityHits:
        return SpecificityHits(
            source=self.source, mismatches=0, perfect=None, near=None,
            error="check disabled",
        )


# --- within the supplied sequence --------------------------------------------


class LocalChecker(SpecificityChecker):
    """Near-match counts within the supplied sequence only."""

    source = "local"

    def __init__(self, sequence: str, mismatches: int, circular: bool = False,
                 label: str = "") -> None:
        super().__init__()
        self.sequence = sequence.upper()
        self.mismatches = mismatches
        self.circular = circular
        self.label = label

    def _haystack(self, length: int) -> str:
        # On a circle, a site spanning the origin is a real site. Appending the
        # first length-1 bases exposes it without duplicating any other match,
        # since no match can start in the appended tail.
        if self.circular and length > 1:
            return self.sequence + self.sequence[: length - 1]
        return self.sequence

    def _count(self, sequence: str) -> SpecificityHits:
        exact, near = count_near_matches(
            self._haystack(len(sequence)), sequence, self.mismatches
        )
        return SpecificityHits(
            source=self.source, mismatches=self.mismatches, perfect=exact, near=near
        )

    def describe(self) -> str:
        where = self.label or "the supplied sequence"
        return (
            f"counted within {where} only, allowing {self.mismatches} mismatch(es). "
            "This cannot show a primer is unique in a genome - only that it is not."
        )


# --- a local genome FASTA -----------------------------------------------------

CHUNK_SIZE = 4_000_000


def _open_fasta(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def iter_fasta_chunks(
    path: pathlib.Path, overlap: int, chunk_size: int = CHUNK_SIZE
) -> Iterator[tuple[str, str, int]]:
    """Yield ``(record_name, sequence, carry_length)`` over a (possibly gzipped) FASTA.

    ``carry_length`` is how many leading bases were already emitted as the tail
    of the previous chunk of the same record, so a scanner can skip the sites
    it has already counted. It is 0 at the start of every record.
    """
    name = ""
    buffer: list[str] = []
    pending = 0
    carry = ""

    with _open_fasta(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if buffer:
                    yield name, carry + "".join(buffer), len(carry)
                name = line[1:].strip().split()[0] if len(line) > 1 else ""
                buffer, pending, carry = [], 0, ""
                continue
            stripped = line.strip()
            if not stripped:
                continue
            buffer.append(stripped)
            pending += len(stripped)
            if pending >= chunk_size:
                text = carry + "".join(buffer)
                yield name, text, len(carry)
                carry = text[-overlap:] if overlap else ""
                buffer, pending = [], 0
        if buffer:
            yield name, carry + "".join(buffer), len(carry)


class FastaChecker(SpecificityChecker):
    """Exact-match counts over a local genome FASTA.

    Streams the assembly once and counts, for every candidate at the same time,
    the exact occurrences of the primer and of its 3'-terminal ``seed`` bases.
    Exact matching is the compromise that makes a genome-scale scan finish in
    under a minute; the 3' count is what catches a primer whose business end
    sits in a repeat even though its full length is unique.
    """

    def __init__(self, path: str | pathlib.Path, seed: int = 15,
                 chunk_size: int = CHUNK_SIZE) -> None:
        super().__init__()
        self.path = pathlib.Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"genome FASTA not found: {self.path}")
        self.seed = seed
        self.chunk_size = chunk_size

    @property
    def source(self) -> str:  # type: ignore[override]
        return f"fasta:{self.path.name}"

    def counts(self, sequences: list[str]) -> dict[str, SpecificityHits]:
        pending = [s.upper() for s in dict.fromkeys(s.upper() for s in sequences)
                   if s.upper() not in self._cache]
        if pending:
            self._scan(pending)
        return {s.upper(): self._cache[s.upper()] for s in sequences}

    def _scan(self, sequences: list[str]) -> None:
        # Patterns: every primer and every 3' seed, on both strands. Counting
        # with str.find keeps the inner loop in C rather than in Python.
        patterns: dict[str, list[str]] = {}
        for sequence in sequences:
            seed = sequence[-self.seed :] if self.seed and len(sequence) > self.seed else ""
            for query in filter(None, (sequence, seed)):
                for oriented in {query, reverse_complement(query)}:
                    patterns.setdefault(oriented, []).append(query)

        tally: dict[str, int] = dict.fromkeys(patterns, 0)
        overlap = max(len(p) for p in patterns) - 1
        scanned = 0
        next_report = 5e8

        for record, chunk, carry in iter_fasta_chunks(
            self.path, overlap=overlap, chunk_size=self.chunk_size
        ):
            upper = chunk.upper()
            scanned += len(chunk) - carry
            if scanned >= next_report:
                logger.info("genome scan: %.1f Gb read (at %s)", scanned / 1e9, record)
                next_report += 5e8
            for pattern in patterns:
                length = len(pattern)
                # A hit fully inside the carried tail was counted last time.
                tally[pattern] += sum(
                    1 for start in find_all(upper, pattern) if start + length > carry
                )

        logger.info("genome scan finished: %.2f Gb", scanned / 1e9)

        for sequence in sequences:
            seed = sequence[-self.seed :] if self.seed and len(sequence) > self.seed else ""
            full = sum(
                tally[o] for o in {sequence, reverse_complement(sequence)}
            )
            seed_hits = (
                sum(tally[o] for o in {seed, reverse_complement(seed)}) if seed else None
            )
            self._cache[sequence] = SpecificityHits(
                source=self.source,
                mismatches=0,
                perfect=full,
                near=full,
                three_prime=seed_hits,
                three_prime_length=len(seed),
            )

    def describe(self) -> str:
        return (
            f"exact matches counted genome-wide in {self.path.name}; sites differing "
            "by one or two bases were not counted, so a near-identical paralogue "
            "would not show up here"
        )


# --- GGGenome ----------------------------------------------------------------


class GGGenomeError(RuntimeError):
    """The service could not be reached, or answered for the wrong genome."""


class GGGenomeChecker(SpecificityChecker):
    """Genome-wide near-match counts from the DBCLS GGGenome service."""

    def __init__(
        self,
        database: str,
        mismatches: int = 2,
        expect: str = "",
        delay: float = 1.0,
        timeout: float = 60.0,
        url: str = GGGENOME_URL,
        max_hits: int = 1,
    ) -> None:
        super().__init__()
        self.database = database
        self.mismatches = mismatches
        self.expect = expect
        self.delay = delay
        self.timeout = timeout
        self.url = url.rstrip("/")
        #: Above this many near-matches the exact count is not worth a request.
        self.max_hits = max_hits
        self._last_request = 0.0

    @property
    def source(self) -> str:  # type: ignore[override]
        return f"gggenome:{self.database}"

    def _get(self, sequence: str, k: int) -> dict:
        # DBCLS ask for gentle use of the public service.
        wait = self.delay - (time.monotonic() - self._last_request)
        if wait > 0:
            time.sleep(wait)

        url = (
            f"{self.url}/{urllib.parse.quote(self.database)}/{k}/"
            f"{urllib.parse.quote(sequence)}.json"
        )
        request = urllib.request.Request(
            url, headers={"User-Agent": "genotyping-primer-designer"}
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as handle:
                payload = json.load(handle)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
            raise GGGenomeError(str(exc)) from exc
        finally:
            self._last_request = time.monotonic()

        self._validate(payload)
        return payload

    def _validate(self, payload: dict) -> None:
        """Reject an error response, or one answered for the wrong genome."""
        error = payload.get("error")
        if error and error != "none":
            raise GGGenomeError(str(error))

        # GGGenome answers an unrecognised database name with the *human*
        # genome rather than an error, so a wrong database would silently score
        # a mouse primer against hg38. Confirm the organism before believing
        # any count that comes back.
        name = str(payload.get("database", ""))
        if self.expect and self.expect.lower() not in name.lower():
            raise GGGenomeError(
                f"database {self.database!r} returned {name!r}, which is not the "
                f"expected {self.expect} genome; pass --gggenome-db explicitly"
            )

    @staticmethod
    def _summary_total(payload: dict) -> tuple[int, bool]:
        """Total hits over the query and its reverse complement."""
        summary = payload.get("summary") or []
        if not summary:
            return len(payload.get("results") or []), False
        total = sum(int(entry.get("count", 0)) for entry in summary)
        approximate = any(bool(entry.get("count_is_approx")) for entry in summary)
        return total, approximate

    def _count(self, sequence: str) -> SpecificityHits:
        if not is_unambiguous(sequence):
            return SpecificityHits(
                source=self.source, mismatches=self.mismatches, perfect=None, near=None,
                error="sequence contains non-ACGT bases",
            )
        try:
            # Ask the expensive question first: if the primer already has more
            # near-matches than we tolerate, the exact count cannot rescue it.
            near, approximate = self._summary_total(self._get(sequence, self.mismatches))
            if near > self.max_hits:
                return SpecificityHits(
                    source=self.source, mismatches=self.mismatches,
                    perfect=None, near=near, approximate=approximate,
                )
            if self.mismatches == 0:
                perfect = near  # the same query; no second request needed
            else:
                perfect, _ = self._summary_total(self._get(sequence, 0))
        except GGGenomeError as exc:
            logger.warning("GGGenome check failed for %s: %s", sequence, exc)
            return SpecificityHits(
                source=self.source, mismatches=self.mismatches,
                perfect=None, near=None, error=str(exc),
            )
        return SpecificityHits(
            source=self.source, mismatches=self.mismatches,
            perfect=perfect, near=near, approximate=approximate,
        )

    def describe(self) -> str:
        return (
            f"counted genome-wide against {self.database} via GGGenome, "
            f"allowing {self.mismatches} mismatch(es)"
        )


def resolve_database(species: str, override: str | None) -> tuple[str, str] | None:
    """Pick the GGGenome database, or ``None`` if the species has none.

    An override naming a species is looked up like any other species; one
    naming a database is taken at face value - the caller asked for it - so it
    carries no organism expectation and skips the hg38 fallback guard.
    """
    if override:
        known = GGGENOME_DATABASES.get(normalise_species(override))
        return known if known else (override, "")
    return GGGENOME_DATABASES.get(normalise_species(species))


__all__ = [
    "CHUNK_SIZE",
    "GGGENOME_DATABASES",
    "FastaChecker",
    "GGGenomeChecker",
    "GGGenomeError",
    "LocalChecker",
    "NullChecker",
    "SpecificityChecker",
    "iter_fasta_chunks",
    "normalise_species",
    "resolve_database",
]
