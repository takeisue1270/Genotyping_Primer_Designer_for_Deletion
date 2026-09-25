"""Reading a target sequence in, and writing annotated maps back out.

Input is whatever the user already has: a FASTA file, a GenBank/SnapGene map,
or a bare stretch of bases pasted into a file. Output is GenBank, because that
is what opens in SnapGene with the primers already drawn on it.

Both are handled here in plain Python. The parser only needs the sequence, its
name and - for a plasmid - whether it is circular, so pulling in a full
sequence-format library would cost a dependency for three fields.
"""

from __future__ import annotations

import datetime
import gzip
import pathlib
import re
import textwrap
from dataclasses import dataclass, field

_DNA_LINE = re.compile(r"^[ACGTURYKMSWBDHVNacgturykmswbdhvn\s\d]+$")


class SequenceError(ValueError):
    """The file could not be read as a nucleotide sequence."""


def _open(path: pathlib.Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


# --- reading -----------------------------------------------------------------


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """``[(name, sequence), ...]`` for every record in a FASTA string."""
    records: list[tuple[str, str]] = []
    name = ""
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if chunks or name:
                records.append((name, "".join(chunks)))
            header = line[1:].strip()
            name = header.split()[0] if header else ""
            chunks = []
        elif line.strip():
            chunks.append("".join(line.split()))
    if chunks or name:
        records.append((name, "".join(chunks)))
    return records


def parse_genbank(text: str) -> tuple[str, str, bool]:
    """``(name, sequence, circular)`` from the first record of a GenBank string."""
    name = ""
    circular = False
    chunks: list[str] = []
    in_origin = False

    for line in text.splitlines():
        if line.startswith("LOCUS"):
            fields = line.split()
            if len(fields) > 1:
                name = fields[1]
            circular = " circular " in f" {line.lower()} "
        elif line.startswith("ORIGIN"):
            in_origin = True
        elif line.startswith("//"):
            break
        elif in_origin:
            chunks.append("".join(c for c in line if c.isalpha()))

    if not chunks:
        raise SequenceError("GenBank record has no ORIGIN sequence")
    # GenBank stores the sequence in lower case; every caller wants ACGT.
    return name, "".join(chunks).upper(), circular


def read_sequence(
    path: str | pathlib.Path, record: str | None = None
) -> tuple[str, str, bool, list[str]]:
    """``(name, sequence, circular, notes)`` for a FASTA, GenBank or raw file.

    The format is taken from the content, not the extension, because ``.txt``
    and ``.seq`` are used for all three. ``record`` selects one entry of a
    multi-record FASTA by name; without it the first is used and the rest are
    reported in ``notes``.
    """
    path = pathlib.Path(path)
    if not path.exists():
        raise SequenceError(f"sequence file not found: {path}")
    with _open(path) as handle:
        text = handle.read()

    notes: list[str] = []
    stripped = text.lstrip()

    if stripped.startswith(">"):
        records = parse_fasta(text)
        if not records:
            raise SequenceError(f"{path}: no FASTA records")
        if record is not None:
            matches = [r for r in records if r[0] == record]
            if not matches:
                raise SequenceError(
                    f"{path}: no record named {record!r}; found "
                    + ", ".join(r[0] for r in records[:10])
                )
            name, sequence = matches[0]
        else:
            name, sequence = records[0]
            if len(records) > 1:
                notes.append(
                    f"{path.name} holds {len(records)} records; used {name!r}. "
                    "Pass --record to choose another."
                )
        circular = False
    elif stripped.startswith("LOCUS"):
        name, sequence, circular = parse_genbank(text)
    else:
        first = stripped.splitlines()[0] if stripped.splitlines() else ""
        if not first or not _DNA_LINE.match(first):
            raise SequenceError(
                f"{path}: not FASTA, not GenBank, and the first line is not DNA"
            )
        name = path.stem
        sequence = "".join(c for c in text if c.isalpha())
        circular = False

    sequence = sequence.upper()
    if not sequence:
        raise SequenceError(f"{path}: the record is empty")
    ambiguous = sorted(set(sequence) - set("ACGTU"))
    if ambiguous:
        notes.append(
            f"{name}: the sequence carries ambiguity code(s) {''.join(ambiguous)} "
            "(assembly gaps show up as N); no primer is placed over those positions"
        )
    sequence = sequence.replace("U", "T")
    return name or path.stem, sequence, circular, notes


# --- writing -----------------------------------------------------------------


@dataclass
class Feature:
    """One GenBank feature over a 0-based half-open span."""

    key: str
    start: int
    end: int
    strand: int = 1
    qualifiers: list[tuple[str, str]] = field(default_factory=list)

    def location(self) -> str:
        span = f"{self.start + 1}..{self.end}"
        return f"complement({span})" if self.strand < 0 else span


def _format_qualifier(name: str, value: str) -> list[str]:
    quoted = value.replace('"', "'")
    text = f'/{name}="{quoted}"'
    return textwrap.wrap(text, width=79, subsequent_indent=" " * 21) or [text]


def format_genbank(
    name: str,
    sequence: str,
    features: list[Feature],
    *,
    circular: bool = False,
    definition: str = "",
    date: str | None = None,
) -> str:
    """A GenBank record as text, ready for SnapGene or Benchling."""
    stamp = date or datetime.date.today().strftime("%d-%b-%Y").upper()
    topology = "circular" if circular else "linear"
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)[:24] or "sequence"

    lines = [
        f"LOCUS       {safe:<24}{len(sequence):>8} bp    DNA     {topology:<8} SYN {stamp}",
        f"DEFINITION  {definition or safe}",
        f"ACCESSION   {safe}",
        "KEYWORDS    .",
        "SOURCE      synthetic DNA construct",
        "  ORGANISM  synthetic DNA construct",
        "FEATURES             Location/Qualifiers",
    ]

    for feature in sorted(features, key=lambda f: (f.start, f.end)):
        lines.append(f"     {feature.key:<16}{feature.location()}")
        for qualifier, value in feature.qualifiers:
            lines.extend(f"{'':21}{line.lstrip()}" for line in _format_qualifier(qualifier, value))

    lines.append("ORIGIN")
    lower = sequence.lower()
    for offset in range(0, len(lower), 60):
        block = lower[offset : offset + 60]
        groups = " ".join(block[i : i + 10] for i in range(0, len(block), 10))
        lines.append(f"{offset + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def write_genbank(
    path: str | pathlib.Path,
    name: str,
    sequence: str,
    features: list[Feature],
    *,
    circular: bool = False,
    definition: str = "",
) -> pathlib.Path:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        format_genbank(
            name, sequence, features, circular=circular, definition=definition
        ),
        encoding="utf-8",
    )
    return path


__all__ = [
    "Feature",
    "SequenceError",
    "format_genbank",
    "parse_fasta",
    "parse_genbank",
    "read_sequence",
    "write_genbank",
]
