"""Putting a run together: sequence in, deletion resolved, primers out.

This is the module to call from a script or a notebook. The CLI is a thin
wrapper over ``design``.
"""

from __future__ import annotations

import logging

from .config import PRIMER_3PRIME_SEED, DesignConfig
from .ensembl import EnsemblClient, parse_locus
from .models import Deletion, DesignResult, Target
from .primers import annotate_specificity, collect_candidates, rank_pairs
from .seqio import read_sequence
from .specificity import (
    FastaChecker,
    GGGenomeChecker,
    LocalChecker,
    NullChecker,
    SpecificityChecker,
    normalise_species,
    resolve_database,
)

logger = logging.getLogger(__name__)


class DesignError(RuntimeError):
    """A run could not be set up: bad coordinates, or nothing to design against."""


# --- getting a target --------------------------------------------------------


def target_from_file(
    path: str, name: str | None = None, record: str | None = None, circular: bool | None = None
) -> tuple[Target, list[str]]:
    """Read a FASTA, GenBank or raw sequence file into a ``Target``."""
    parsed_name, sequence, file_circular, notes = read_sequence(path, record=record)
    return (
        Target(
            name=name or parsed_name,
            sequence=sequence,
            circular=file_circular if circular is None else circular,
            source=str(path),
        ),
        notes,
    )


def target_from_genome(
    species: str,
    region: str,
    start: int,
    end: int,
    name: str | None = None,
    client: EnsemblClient | None = None,
) -> tuple[Target, list[str]]:
    """Fetch ``region:start..end`` (1-based inclusive) from a reference genome."""
    species = normalise_species(species)
    client = client or EnsemblClient()
    sequence = client.sequence_region(species, region, start, end)
    if len(sequence) != end - start + 1:
        raise DesignError(
            f"Ensembl returned {len(sequence):,} bp for {region}:{start}-{end}, "
            f"expected {end - start + 1:,}; the region may run off the end of {region}"
        )
    assembly = client.assembly_name(species)
    notes = []
    if not assembly:
        notes.append(
            "Ensembl did not report an assembly name (the /info/assembly endpoint is "
            "the flakiest part of the service); the sequence itself is unaffected"
        )
    return (
        Target(
            name=name or f"{region}_{start}_{end}",
            sequence=sequence,
            circular=False,
            source=f"Ensembl {species} {region}:{start:,}-{end:,}",
            species=species,
            assembly=assembly,
            region=region,
            region_start=start,
        ),
        notes,
    )


# --- resolving the deletion --------------------------------------------------


def resolve_span(
    target: Target, spec: str, coordinates: str = "auto"
) -> tuple[int, int]:
    """Turn ``--delete`` into 0-based offsets into the sequence as supplied.

    ``coordinates`` is ``local`` (1-based into the supplied sequence),
    ``genomic`` (assembly coordinates, only if the sequence was fetched from
    one), or ``auto``: genomic when the span falls inside the fetched window,
    local otherwise. The returned ``end`` may exceed the sequence length for a
    circular target, meaning the deletion crosses the origin.
    """
    region, start, end = parse_locus(spec)
    length = len(target.sequence)

    if region and not target.region:
        raise DesignError(
            f"--delete names region {region!r}, but the supplied sequence carries no "
            "genomic coordinates; drop the region name or fetch the target with --species"
        )
    if region and target.region and region.lstrip("chr") != target.region.lstrip("chr"):
        raise DesignError(
            f"--delete is on {region!r} but the target sequence is from "
            f"{target.region!r}"
        )

    window_end = target.region_start + length - 1
    in_window = bool(target.region) and target.region_start <= start <= window_end
    use_genomic = coordinates == "genomic" or (
        coordinates == "auto" and (bool(region) or in_window)
    )

    if use_genomic:
        if not target.region:
            raise DesignError(
                "genomic coordinates were requested, but the sequence was read from a "
                "file with no genomic position; use local coordinates"
            )
        if not (target.region_start <= start <= window_end):
            raise DesignError(
                f"{start:,} is outside the fetched window "
                f"{target.region}:{target.region_start:,}-{window_end:,}"
            )
        if not (target.region_start <= end <= window_end):
            raise DesignError(
                f"{end:,} is outside the fetched window "
                f"{target.region}:{target.region_start:,}-{window_end:,}"
            )
        start -= target.region_start - 1
        end -= target.region_start - 1

    if start > length or end > length:
        raise DesignError(
            f"the deletion {start:,}-{end:,} runs past the end of a {length:,} bp sequence"
        )
    if end < start and not target.circular:
        raise DesignError(
            f"the deletion end ({end:,}) is before its start ({start:,}); that only "
            "makes sense on a circular sequence, so pass --circular if it is a plasmid"
        )
    return start - 1, end


def prepare(target: Target, start: int, end: int) -> tuple[Target, Deletion]:
    """Rotate a circular target so the deletion has room on both sides.

    Working coordinates are linear. A plasmid deletion sitting near base 1 has
    the rest of the plasmid available to search - but only if the sequence is
    rotated first, so that is done here and recorded in ``Target.rotation`` for
    the report to undo.
    """
    length = len(target.sequence)
    if target.circular:
        deleted = end - start if end > start else length - start + end
        if deleted >= length:
            raise DesignError("the deletion covers the whole sequence")
        rotation = (start - (length - deleted) // 2) % length
        rotated = target.sequence[rotation:] + target.sequence[:rotation]
        working = Target(
            name=target.name,
            sequence=rotated,
            circular=True,
            rotation=rotation,
            source=target.source,
            species=target.species,
            assembly=target.assembly,
            region=target.region,
            region_start=target.region_start,
        )
        offset = (start - rotation) % length
        return working, Deletion(start=offset, end=offset + deleted)

    if end <= start:
        raise DesignError(f"the deletion is empty ({start + 1}-{end})")
    return target, Deletion(start=start, end=end)


def with_insert(deletion: Deletion, insert: str, insert_name: str = "") -> Deletion:
    """Attach the sequence that replaces the deleted block, if any."""
    cleaned = "".join(insert.split()).upper()
    if cleaned and set(cleaned) - set("ACGT"):
        raise DesignError("the insert must be A/C/G/T only")
    return Deletion(
        start=deletion.start, end=deletion.end, insert=cleaned, insert_name=insert_name
    )


def edited_sequence(target: Target, deletion: Deletion) -> str:
    """The allele after the deletion, with the insert in its place."""
    return (
        target.sequence[: deletion.start]
        + deletion.insert
        + target.sequence[deletion.end :]
    )


# --- picking a specificity backend -------------------------------------------


def choose_checker(target: Target, config: DesignConfig) -> SpecificityChecker:
    """Build the backend named by ``config.specificity``, resolving ``auto``."""
    backend = config.specificity
    species = normalise_species(config.species or target.species)

    if backend == "auto":
        if config.genome_fasta:
            backend = "fasta"
        elif resolve_database(species, config.gggenome_db):
            backend = "gggenome"
        else:
            backend = "local"

    if backend == "none":
        return NullChecker()
    if backend == "local":
        label = f"the supplied {'plasmid' if target.circular else 'sequence'}"
        return LocalChecker(
            target.sequence, config.mismatches, circular=target.circular, label=label
        )
    if backend == "fasta":
        if not config.genome_fasta:
            raise DesignError("--specificity fasta needs --genome-fasta")
        return FastaChecker(config.genome_fasta, seed=PRIMER_3PRIME_SEED)

    resolved = resolve_database(species, config.gggenome_db)
    if not resolved:
        raise DesignError(
            f"no GGGenome database is known for {species or 'this target'}; pass "
            "--gggenome-db with one from https://gggenome.dbcls.jp/, or use "
            "--specificity local (within the supplied sequence) or --specificity none"
        )
    database, expect = resolved
    return GGGenomeChecker(
        database=database,
        mismatches=config.mismatches,
        expect=expect,
        delay=config.gggenome_delay,
        max_hits=config.max_hits,
    )


# --- the run -----------------------------------------------------------------


def design(
    target: Target,
    deletion: Deletion,
    config: DesignConfig,
    checker: SpecificityChecker | None = None,
    extra_notes: list[str] | None = None,
) -> DesignResult:
    """Design the genotyping PCR for one deletion."""
    edited = edited_sequence(target, deletion)
    candidates = collect_candidates(target, deletion, config)

    checker = checker or choose_checker(target, config)
    candidates = annotate_specificity(candidates, checker)

    pairs, notes = rank_pairs(target, deletion, edited, config, candidates)

    notes = list(extra_notes or []) + notes
    notes.append(f"Specificity: {checker.describe().rstrip('.')}.")

    shift = deletion.net_change
    if shift == 0:
        notes.append(
            "The edit does not change the length of the allele, so no PCR across it can "
            "separate the two genotypes by size. Genotype by sequencing across the "
            "junction, or by a primer that sits inside the insert."
        )
    elif abs(shift) < config.min_band_difference:
        notes.append(
            f"The two bands differ by only {abs(shift)} bp, which a normal agarose gel "
            "will not resolve. Run a higher-percentage gel or a capillary sizer, delete "
            "more sequence, or genotype by sequencing across the junction."
        )
    elif shift < 0:
        notes.append(
            f"The insert is longer than the deletion, so the edited band is the *larger* "
            f"one, by {-shift} bp."
        )

    if config.homology_arm == 0:
        notes.append(
            "--arm is 0, so primers were only kept out of the deleted region itself. If a "
            "repair template (ssODN or targeting vector) is used, set --arm to its homology "
            "arm length; a primer inside an arm amplifies the donor and scores every clone "
            "as edited."
        )

    return DesignResult(
        target=target,
        deletion=deletion,
        edited_sequence=edited,
        pairs=tuple(pairs),
        notes=tuple(notes),
        specificity_backend=checker.source,
    )


__all__ = [
    "DesignError",
    "choose_checker",
    "design",
    "edited_sequence",
    "prepare",
    "resolve_span",
    "target_from_file",
    "target_from_genome",
    "with_insert",
]
