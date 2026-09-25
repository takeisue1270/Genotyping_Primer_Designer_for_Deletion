"""Turning a design into something a person can read, order from, or open.

Three outputs, each aimed at a different moment:

* the printed report, read once while deciding whether the design is sound;
* a TSV, pasted into an oligo order and kept with the lab record;
* two GenBank maps - before and after the edit - opened in SnapGene when a
  band on the gel does not look like it should.
"""

from __future__ import annotations

import csv
import pathlib

from .config import DesignConfig
from .models import DesignResult, Primer, PrimerPair
from .seqio import Feature, write_genbank

TSV_COLUMNS = [
    "pair",
    "primer",
    "name",
    "sequence",
    "length",
    "tm_c",
    "gc_percent",
    "strand",
    "position",
    "clearance_bp",
    "specificity",
    "unedited_bp",
    "edited_bp",
    "difference_bp",
    "warnings",
]


def _position(result: DesignResult, primer: Primer) -> str:
    return result.target.describe_span(primer.start, primer.end)


# --- the printed report ------------------------------------------------------


def _name_width(result: DesignResult) -> int:
    """Wide enough for the longest primer name, so the columns stay columns."""
    names = [p.name for pair in result.pairs for p in (pair.forward, pair.reverse)]
    return max([18, *(len(name) for name in names)])


def _primer_line(result: DesignResult, primer: Primer, width: int) -> str:
    return (
        f"    {primer.name:<{width}} {primer.sequence}\n"
        f"    {'':<{width}} {primer.length} nt  Tm {primer.tm:.1f} C  "
        f"GC {primer.gc * 100:.0f} %  {primer.strand_symbol} strand  "
        f"{_position(result, primer)}  "
        f"{primer.clearance} bp clear of the edit\n"
        f"    {'':<{width}} {primer.specificity_note()}"
    )


def _pair_block(result: DesignResult, pair: PrimerPair, index: int, width: int) -> str:
    heading = "  Recommended pair" if index == 1 else f"  Alternative {index - 1}"
    lines = [
        heading,
        _primer_line(result, pair.forward, width),
        _primer_line(result, pair.reverse, width),
        (
            f"    bands: unedited {pair.wt_product:,} bp / edited {pair.ko_product:,} bp"
            f"   (difference {abs(pair.size_shift):,} bp)"
        ),
    ]
    for warning in pair.warnings:
        lines.append(f"    ! {warning}")
    return "\n".join(lines)


def text_report(result: DesignResult, config: DesignConfig) -> str:
    """The human-readable summary printed at the end of a run."""
    target, deletion = result.target, result.deletion
    lines = [
        f"{target.name}   {target.source or 'local sequence'}",
        f"  {len(target):,} bp {'circular' if target.circular else 'linear'}"
        f"   coordinates: {target.coordinate_system}",
        "",
        f"  deletion   {target.describe_span(deletion.start, deletion.end)}"
        f"   {deletion.length:,} bp removed",
    ]
    if deletion.insert:
        shown = (
            deletion.insert if deletion.insert_length <= 60 else deletion.insert[:57] + "..."
        )
        # The name is only worth a column when it says more than "insert" does.
        named = (
            f"{deletion.insert_name}   "
            if deletion.insert_name and deletion.insert_name != "insert"
            else ""
        )
        lines.append(f"  insert     {named}{deletion.insert_length} bp   {shown}")
    if deletion.net_change:
        direction = "shorter" if deletion.net_change > 0 else "longer"
        lines.append(
            f"  net        the edited allele is {abs(deletion.net_change):,} bp {direction}"
        )
    else:
        lines.append("  net        the edited allele is the same length")
    lines.append(
        f"  excluded   {config.homology_arm} bp homology arm on each side; "
        f"{config.min_clearance} bp of clearance asked for beyond it"
    )
    lines.append("")

    if not result.pairs:
        lines.append("  No primer pair was found.")
    width = _name_width(result)
    for index, pair in enumerate(result.pairs, start=1):
        lines.append(_pair_block(result, pair, index, width))
        lines.append("")

    if result.notes:
        lines.append("  Notes")
        for note in result.notes:
            lines.append(f"    - {note}")
    return "\n".join(lines)


# --- TSV ---------------------------------------------------------------------


def write_tsv(path: str | pathlib.Path, result: DesignResult) -> pathlib.Path:
    """One row per primer, carrying its pair's band sizes and every warning."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(TSV_COLUMNS)
        for index, pair in enumerate(result.pairs, start=1):
            for role, primer in (("forward", pair.forward), ("reverse", pair.reverse)):
                writer.writerow(
                    [
                        index,
                        role,
                        primer.name,
                        primer.sequence,
                        primer.length,
                        f"{primer.tm:.1f}",
                        f"{primer.gc * 100:.0f}",
                        primer.strand_symbol,
                        _position(result, primer).replace(",", ""),  # a table, not prose
                        primer.clearance,
                        primer.specificity_note(),
                        pair.wt_product,
                        pair.ko_product,
                        pair.size_shift,
                        " | ".join(pair.warnings),
                    ]
                )
    return path


# --- GenBank -----------------------------------------------------------------


def _primer_feature(primer: Primer, start: int) -> Feature:
    return Feature(
        key="primer_bind",
        start=start,
        end=start + primer.length,
        strand=primer.strand,
        qualifiers=[
            ("label", primer.name),
            ("note", f"Tm {primer.tm:.1f} C, GC {primer.gc * 100:.0f}%"),
            ("note", f"sequence: {primer.sequence}"),
        ],
    )


def unedited_features(result: DesignResult, config: DesignConfig) -> list[Feature]:
    deletion = result.deletion
    features = [
        Feature(
            key="misc_feature",
            start=deletion.start,
            end=deletion.end,
            qualifiers=[
                ("label", f"deleted region ({deletion.length} bp)"),
                ("note", "removed in the edited allele"),
            ],
        )
    ]
    arm = config.homology_arm
    if arm:
        for start, end, side in (
            (max(0, deletion.start - arm), deletion.start, "left"),
            (deletion.end, min(len(result.target), deletion.end + arm), "right"),
        ):
            features.append(
                Feature(
                    key="misc_feature",
                    start=start,
                    end=end,
                    qualifiers=[
                        ("label", f"homology arm {side} ({end - start} bp)"),
                        ("note", "no primer may anneal here: it would amplify the donor"),
                    ],
                )
            )
    for pair in result.pairs:
        features.append(_primer_feature(pair.forward, pair.forward.start))
        features.append(_primer_feature(pair.reverse, pair.reverse.start))
    return features


def edited_features(result: DesignResult) -> list[Feature]:
    deletion = result.deletion
    sequence = result.edited_sequence
    features: list[Feature] = []

    if deletion.insert:
        features.append(
            Feature(
                key="misc_feature",
                start=deletion.start,
                end=deletion.start + deletion.insert_length,
                qualifiers=[
                    ("label", deletion.insert_name or "insert"),
                    ("note", "left in place of the deleted region"),
                ],
            )
        )
    else:
        features.append(
            Feature(
                key="misc_feature",
                start=max(0, deletion.start - 1),
                end=min(len(sequence), deletion.start + 1),
                qualifiers=[("label", "deletion junction")],
            )
        )

    for pair in result.pairs:
        for primer in (pair.forward, pair.reverse):
            hit = sequence.upper().find(primer.top_strand())
            if hit >= 0:
                features.append(_primer_feature(primer, hit))
    return features


def write_maps(
    directory: str | pathlib.Path, result: DesignResult, config: DesignConfig
) -> list[pathlib.Path]:
    """Write the before and after GenBank maps; returns the paths written."""
    directory = pathlib.Path(directory)
    name = result.target.name

    # A circular target is rotated internally so the deletion has room on both
    # sides. The molecule is the same one, but base 1 of these maps is not base
    # 1 of the file the user supplied, and a map that does not say so is a trap.
    rotation = result.target.rotation
    origin = (
        f" [rotated: base 1 here is base {rotation + 1:,} of the supplied sequence]"
        if rotation
        else ""
    )

    written = [
        write_genbank(
            directory / f"{name}_unedited.gb",
            f"{name}_unedited",
            result.target.sequence,
            unedited_features(result, config),
            circular=result.target.circular,
            definition=f"{name} before the deletion, with genotyping primers{origin}",
        ),
        write_genbank(
            directory / f"{name}_edited.gb",
            f"{name}_edited",
            result.edited_sequence,
            edited_features(result),
            circular=result.target.circular,
            definition=(
                f"{name} after deleting {result.deletion.length} bp"
                + (
                    f" and inserting {result.deletion.insert_length} bp"
                    if result.deletion.insert
                    else ""
                )
                + origin
            ),
        ),
    ]
    return written


__all__ = [
    "TSV_COLUMNS",
    "edited_features",
    "text_report",
    "unedited_features",
    "write_maps",
    "write_tsv",
]
