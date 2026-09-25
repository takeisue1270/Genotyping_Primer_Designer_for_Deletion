"""Command-line entry point.

    genotyping-primers --sequence locus.fa --delete 4800-5300
    genotyping-primers --species mouse --delete 11:69,000,000-69,001,000
    genotyping-primers --sequence plasmid.gb --circular --delete 1200-1700 \\
                       --specificity local
"""

from __future__ import annotations

import argparse
import logging
import pathlib
import sys

from . import __version__
from .config import DEFAULT_SPECIES, SPECIFICITY_BACKENDS, DesignConfig
from .design import (
    DesignError,
    design,
    prepare,
    resolve_span,
    target_from_file,
    target_from_genome,
    with_insert,
)
from .ensembl import EnsemblError, parse_locus
from .primers import PrimerError
from .report import text_report, write_maps, write_tsv
from .seqio import SequenceError
from .specificity import GGGENOME_DATABASES, normalise_species

logger = logging.getLogger("genotyping_primers")


def _pair(text: str, cast, name: str) -> tuple:
    parts = [p for p in text.replace("..", ",").replace("-", ",").split(",") if p.strip()]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f"{name} needs two values, e.g. '58,62'")
    try:
        return cast(parts[0]), cast(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"{name}: {exc}") from exc


def _int_pair(text: str) -> tuple[int, int]:
    return _pair(text, int, "length range")


def _float_pair(text: str) -> tuple[float, float]:
    return _pair(text, float, "range")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genotyping-primers",
        description=(
            "Design the PCR that tells a deleted allele from an unedited one by size. "
            "Give it a sequence and the region to remove; it returns a primer pair "
            "outside the repair template's homology arms, and both band sizes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    source = parser.add_argument_group("target sequence")
    source.add_argument(
        "--sequence", metavar="PATH",
        help="FASTA, GenBank/SnapGene or raw sequence file holding the target",
    )
    source.add_argument(
        "--record", metavar="NAME",
        help="which record to use from a multi-record FASTA (default: the first)",
    )
    source.add_argument(
        "--species", metavar="NAME", default="",
        help=(
            "Ensembl species to fetch from and to check specificity against, "
            f"e.g. mouse, human, {DEFAULT_SPECIES}"
        ),
    )
    source.add_argument(
        "--region", metavar="CHR:START-END",
        help="fetch this window from the reference genome instead of reading a file",
    )
    source.add_argument(
        "--pad", type=int, metavar="BP",
        help=(
            "when the window is derived from --delete, how much flanking sequence to "
            "fetch on each side (default: arm + span + 200)"
        ),
    )
    source.add_argument(
        "--circular", action="store_true",
        help="treat the sequence as circular (a plasmid); GenBank topology is honoured anyway",
    )
    source.add_argument("--name", metavar="NAME", help="name for the outputs")

    edit = parser.add_argument_group("the edit")
    edit.add_argument(
        "--delete", required=True, metavar="SPEC",
        help=(
            "the region to remove, 1-based inclusive: '4800-5300', '4800..5300', "
            "'4800+500', or genomic '11:69,000,000-69,001,000'"
        ),
    )
    edit.add_argument(
        "--coordinates", choices=("auto", "local", "genomic"), default="auto",
        help=(
            "how to read --delete: local (1-based into the supplied sequence), genomic "
            "(assembly coordinates), or auto - genomic when the span falls inside a "
            "fetched window, local otherwise (default: %(default)s)"
        ),
    )
    edit.add_argument(
        "--insert", metavar="SEQ", default="",
        help=(
            "bases left in place of the deleted block, e.g. a tag or a restriction "
            "site; the band sizes account for it (default: a plain deletion)"
        ),
    )

    placement = parser.add_argument_group("where primers may sit")
    placement.add_argument(
        "--arm", type=int, default=DesignConfig.homology_arm, metavar="BP",
        help=(
            "homology arm of the repair template; primers are excluded from this many "
            "bases on each side of the deletion (default: %(default)s, 0 if no donor)"
        ),
    )
    placement.add_argument(
        "--min-clearance", type=int, default=DesignConfig.min_clearance, metavar="BP",
        help="warn when a primer is closer than this to the excluded zone (default: %(default)s)",
    )
    placement.add_argument(
        "--span", type=int, default=DesignConfig.search_span, metavar="BP",
        help="how far beyond the excluded zone to search (default: %(default)s)",
    )

    properties = parser.add_argument_group("primer properties")
    properties.add_argument(
        "--tm", type=float, default=DesignConfig.tm_target, metavar="C",
        help="target melting temperature (default: %(default)s)",
    )
    properties.add_argument(
        "--tm-range", type=_float_pair, metavar="LO,HI",
        help="acceptable Tm window (default: 58,62)",
    )
    properties.add_argument(
        "--max-tm-diff", type=float, default=DesignConfig.max_tm_diff, metavar="C",
        help="largest Tm difference within a pair (default: %(default)s)",
    )
    properties.add_argument(
        "--length", type=_int_pair, metavar="LO,HI",
        help="acceptable primer length (default: 20,30)",
    )
    properties.add_argument(
        "--gc", type=_float_pair, metavar="LO,HI",
        help="acceptable GC fraction (default: 0.4,0.6)",
    )

    specificity = parser.add_argument_group("specificity")
    specificity.add_argument(
        "--specificity", choices=SPECIFICITY_BACKENDS, default="auto",
        help=(
            "auto: genome FASTA if given, else GGGenome if the species has a database, "
            "else within the supplied sequence. Use 'local' for plasmid work and 'none' "
            "to skip the check (default: %(default)s)"
        ),
    )
    specificity.add_argument(
        "--gggenome-db", metavar="DB",
        help=(
            "GGGenome database to check against, e.g. mm10 for a line still on GRCm38; "
            "see https://gggenome.dbcls.jp/"
        ),
    )
    specificity.add_argument(
        "--genome-fasta", metavar="PATH",
        help="local genome FASTA for an offline, exact-match genome-wide count",
    )
    specificity.add_argument(
        "--mismatches", type=int, default=DesignConfig.mismatches, metavar="N",
        help="mismatches allowed when counting off-target sites (default: %(default)s)",
    )

    output = parser.add_argument_group("output")
    output.add_argument(
        "-o", "--outdir", metavar="DIR",
        help="write the TSV and the GenBank maps here (default: print only)",
    )
    output.add_argument(
        "--alternatives", type=int, default=DesignConfig.alternatives, metavar="N",
        help="how many runner-up pairs to report (default: %(default)s)",
    )
    output.add_argument("--no-maps", action="store_true", help="skip the GenBank maps")
    output.add_argument(
        "--min-band-difference", type=int, default=DesignConfig.min_band_difference,
        metavar="BP",
        help="below this the two bands are called unresolvable (default: %(default)s)",
    )
    output.add_argument("-v", "--verbose", action="store_true", help="log progress")
    output.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def _resolve_insert(text: str) -> tuple[str, str]:
    """``("insert", <bases>)``, or ``("", "")`` for a plain deletion."""
    if not text:
        return "", ""
    return "insert", text


def _config_from(args: argparse.Namespace) -> DesignConfig:
    kwargs: dict = {
        "homology_arm": args.arm,
        "min_clearance": args.min_clearance,
        "search_span": args.span,
        "tm_target": args.tm,
        "max_tm_diff": args.max_tm_diff,
        "specificity": args.specificity,
        "species": normalise_species(args.species) if args.species else "",
        "gggenome_db": args.gggenome_db,
        "genome_fasta": args.genome_fasta,
        "mismatches": args.mismatches,
        "alternatives": args.alternatives,
        "min_band_difference": args.min_band_difference,
    }
    if args.tm_range:
        kwargs["tm_range"] = args.tm_range
    if args.length:
        kwargs["length_range"] = (int(args.length[0]), int(args.length[1]))
    if args.gc:
        kwargs["gc_range"] = args.gc
    return DesignConfig(**kwargs)


def _load_target(args: argparse.Namespace, config: DesignConfig):
    """The target sequence, from a file or from a reference genome."""
    if args.sequence:
        return target_from_file(
            args.sequence, name=args.name, record=args.record,
            circular=True if args.circular else None,
        )

    species = normalise_species(args.species) if args.species else ""
    if not species:
        raise DesignError(
            "give either --sequence (a local file) or --species (to fetch from Ensembl); "
            f"species with a specificity database: {', '.join(sorted(GGGENOME_DATABASES))}"
        )

    if args.region:
        region, start, end = parse_locus(args.region)
        if not region:
            raise DesignError("--region needs a chromosome, e.g. 11:69000000-69010000")
    else:
        region, start, end = parse_locus(args.delete)
        if not region:
            raise DesignError(
                "without --sequence or --region, --delete must name a chromosome, "
                "e.g. --delete 11:69,000,000-69,001,000"
            )
        pad = args.pad if args.pad is not None else config.homology_arm + config.search_span + 200
        start, end = max(1, start - pad), end + pad
    return target_from_genome(species, region, start, end, name=args.name)


def main(argv: list[str] | None = None) -> int:
    # A Japanese or Central-European console is not UTF-8, and a design is not
    # worth losing to an encoding error in the report.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")  # type: ignore[union-attr]
        except (AttributeError, ValueError):  # pragma: no cover - exotic streams
            pass

    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(message)s",
    )

    try:
        config = _config_from(args)
        target, notes = _load_target(args, config)
        start, end = resolve_span(target, args.delete, args.coordinates)
        target, deletion = prepare(target, start, end)
        insert_name, insert = _resolve_insert(args.insert)
        deletion = with_insert(deletion, insert, insert_name)
        result = design(target, deletion, config, extra_notes=notes)
    except (DesignError, PrimerError, SequenceError, EnsemblError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(text_report(result, config))

    if args.outdir:
        directory = pathlib.Path(args.outdir)
        written = [write_tsv(directory / f"{result.target.name}_genotyping_primers.tsv", result)]
        if not args.no_maps:
            written += write_maps(directory, result, config)
        print("\n  Written")
        for path in written:
            print(f"    {path}")

    if not result.pairs:
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
