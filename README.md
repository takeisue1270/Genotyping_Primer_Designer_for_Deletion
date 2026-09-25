# Genotyping Primer Designer for Deletion

*[日本語版 README はこちら](README.ja.md)*

Tell it which bases you are deleting. It gives you back the PCR that tells an
edited clone from an unedited one by the size of a band on a gel — both primers
kept clear of the repair template, and both product sizes measured from the
sequences themselves.

```
 ... [ search span ][ arm ][   region you are deleting   ][ arm ][ search span ] ...
       ^ forward primer                                           reverse primer ^
                   |<---- clearance ---->|

        unedited allele  ──────────────────────────────────────────  1,427 bp
        edited allele    ───────────────/\───────────────────────     926 bp
                                        deletion (± an insert)
```

Nothing about this is specific to knocking out an exon. The deleted block is
whatever span you name — an exon, a promoter, a regulatory element, a stuffer
in a plasmid — and the same run works on a genomic locus fetched from Ensembl,
a FASTA on your disk, or a circular construct from SnapGene.

## Why the primers sit where they do

This is the part a general-purpose primer picker gets wrong, and the reason
this tool exists.

**Neither primer may overlap the repair template's homology arms.** An ssODN
donor or a targeting vector carries the sequence immediately flanking the cut.
A primer that anneals inside an arm anneals to the *donor* as readily as to the
chromosome, so leftover transfection reagent amplifies and every clone you
screen — including the untransfected ones — reads as edited. `--arm` is the
length of that arm, and primers are excluded from that many bases on each side
of the deletion. Set `--arm 0` when there is no donor.

**Beyond the arm, distance still buys you something.** Cas9 leaves indels at
the junction, and resection sometimes runs further than intended. A primer 15 bp
past the arm can be destroyed along with the junction; one 200 bp out will not
be. So a primer closer than `--min-clearance` (default 20 bp) is still offered
— on a cramped locus it may be all there is — but it comes with a warning that
says exactly why it is a risk. Nothing is silently dropped, and nothing that
carries a risk is silently accepted.

Ranking follows the same order of priority: a primer that is not unique, or
that sits inside the clearance you asked for, is pushed below every clean pair.
Among clean pairs, Tm and GC decide, then distance from the junction, then
amplicon length.

## Install

Python 3.10 or newer. The only dependency is numpy.

```bash
git clone https://github.com/takeisue1270/Genotyping_Primer_Designer_for_Deletion.git
cd Genotyping_Primer_Designer_for_Deletion
pip install -e .
```

Or run straight from the source tree without installing:

```bash
PYTHONPATH=src python -m genotyping_primers.cli --help
```

## Quick start

### A locus you already have on disk

```bash
genotyping-primers --sequence examples/example_locus.fa --delete 5501-6000
```

```
example_locus   examples/example_locus.fa
  12,000 bp linear   coordinates: example_locus (1-based)

  deletion   5,501-6,000   500 bp removed
  net        the edited allele is 500 bp shorter
  excluded   43 bp homology arm on each side; 20 bp of clearance asked for beyond it

  Recommended pair
    example_locus_F    ACGACGTAAGGTGAACTTTCGGTCT
                       25 nt  Tm 60.0 C  GC 48 %  + strand  5,175-5,199  258 bp clear of the edit
                       local: 1 exact, 1 within 2 mismatch(es)
    example_locus_R    GCCACCCAGATTCATCCACTTAATGC
                       26 nt  Tm 60.0 C  GC 50 %  - strand  6,595-6,620  551 bp clear of the edit
                       local: 1 exact, 1 within 2 mismatch(es)
    bands: unedited 1,446 bp / edited 946 bp   (difference 500 bp)
```

### A region of a reference genome

Name the species and give the deletion in genome coordinates. The flanking
sequence is fetched from Ensembl automatically, and the primers are checked
genome-wide against that assembly through
[GGGenome](https://gggenome.dbcls.jp/).

```bash
genotyping-primers --species mouse --delete 11:69,000,000-69,001,000 -o designs
```

Supported species names include `mouse`, `human`, `rat`, `zebrafish`, `fly`,
`worm`, `yeast` and their Ensembl equivalents (`mus_musculus`, …).

### A plasmid

A construct is circular, and genome-wide uniqueness is not the question being
asked of it — uniqueness *within the plasmid* is. `--specificity local` (the
default when no genome is named) answers that one and never touches the
network. A deletion may cross the origin:

```bash
genotyping-primers --sequence pMyVector.gb --circular --delete 4950-250
```

```
  deletion   4,950-250   501 bp removed
  Recommended pair
    ...  + strand  4,574-4,597  309 bp clear of the edit
    ...  - strand    775-800    481 bp clear of the edit
    bands: unedited 1,427 bp / edited 926 bp   (difference 501 bp)
```

Topology is read from a GenBank `LOCUS` line, so `--circular` is only needed
for a FASTA.

### Leaving something behind

`--insert` puts a sequence in place of the deleted block — a tag, a
restriction site, a barcode — and the band sizes account for it. The size
difference between the two bands is then *deletion − insert*, which is what
the report calls the net change.

```bash
genotyping-primers --sequence locus.fa --delete 5501-6000 --insert GGATCCAAGCTTGA
#   net   the edited allele is 486 bp shorter    (500 removed, 14 put back)
```

## What comes out

Without `-o` the report is printed and nothing is written. With it:

```
designs/
├── example_locus_genotyping_primers.tsv   one row per primer, ready to order
├── example_locus_unedited.gb              the target, primers + arms annotated
└── example_locus_edited.gb                the allele after the edit
```

The TSV carries the sequence, length, Tm, GC, position, clearance, specificity
result, both band sizes and every warning. The two GenBank files open directly
in SnapGene: the unedited map shows the deleted region, both homology arms and
the primer binding sites; the edited map shows the same primers on the allele
they are supposed to detect.

Coordinates in the report and the TSV are the ones you supplied — genome
coordinates for a fetched region, file coordinates for a file, and original
numbering for a plasmid even though a circular target is rotated internally to
give the search room on both sides. The GenBank maps are the one exception: a
rotated plasmid is written in its rotated frame, and the `DEFINITION` line
states which base of your sequence became base 1.

## The warnings, and what to do about them

Each is attached to the pair it applies to, or listed under **Notes** when it
is a property of the edit rather than of a particular pair.

| Warning | What it means | What to do |
| --- | --- | --- |
| *anneals N bp from the homology arm* | Closer to the junction than `--min-clearance`. Clear of the donor, but a junction indel could reach it. | Raise `--span` so there is more sequence to choose from, or accept it knowingly. |
| *no homology arm was declared* | `--arm 0`, so only the deleted block itself was avoided. | If a repair template is used, re-run with `--arm` set to its arm length. |
| *specificity is unknown* | The check did not run, or the service failed. | Re-run, or check by hand. An unchecked primer is never reported as clean. |
| *is not unique* | The primer has other genomic sites within the mismatch budget. | Take an alternative pair, or widen `--span`. |
| *3' N nt occur N times* | The full primer is unique but its 3' end is repetitive, which is where mispriming starts. | Prefer an alternative. |
| *the two bands differ by only N bp* | Under `--min-band-difference` (50 bp): a normal agarose gel will not resolve them. | Delete more, run a higher-percentage gel, or sequence across the junction. |
| *the edit does not change the length* | The insert is the same size as the deletion, so no PCR across it can separate the genotypes. | Sequence across the junction, or use a primer inside the insert. |
| *the unedited band is N bp* | Long enough (> 3 kb) that a failed reaction and a deleted allele look alike. | Keep a positive control on the gel, or lower `--span`. |
| *the edited band is only N bp* | Under 100 bp; easy to run off the end of a gel. | Raise `--span` to push the primers outwards. |

## Specificity

`--specificity` picks how "is this primer unique?" gets answered.

| Backend | Covers | Finds mismatched sites | Network |
| --- | --- | --- | --- |
| `gggenome` | the whole assembly | yes (`--mismatches`, default 2) | yes |
| `fasta` | a local genome FASTA | no — exact matches, plus a separate count for the primer's 3' 15 nt | no |
| `local` | the supplied sequence only | yes | no |
| `none` | nothing | — | no |
| `auto` *(default)* | `fasta` if `--genome-fasta` is given, else `gggenome` if the species has a database, else `local` | | |

Two things worth knowing:

* **`local` on a genomic target can only disprove uniqueness, never establish
  it**, and the report says so rather than implying a clean result. On a
  plasmid it is the appropriate check, not a fallback.
* **GGGenome answers an unknown database name with the human genome instead of
  an error.** Every response is checked against the expected organism before
  its counts are believed, so a typo in `--gggenome-db` fails loudly instead of
  scoring your mouse primer against hg38. Use `--gggenome-db` for an assembly
  that is not the current default — `mm10` for a line still on GRCm38, say.

## Options

| Option | Default | |
| --- | --- | --- |
| `--sequence PATH` | | FASTA, GenBank/SnapGene or raw sequence |
| `--record NAME` | first | which record of a multi-record FASTA |
| `--species NAME` | | fetch from Ensembl, and check specificity against that assembly |
| `--region CHR:START-END` | | explicit fetch window |
| `--pad BP` | arm + span + 200 | flanking sequence fetched around `--delete` |
| `--circular` | GenBank topology | treat the sequence as a circle |
| `--delete SPEC` | *required* | `4800-5300`, `4800..5300`, `4800+500`, `11:69,000,000-69,001,000` |
| `--coordinates` | auto | read `--delete` as `local`, `genomic`, or decide automatically |
| `--insert SEQ` | none | bases left in place of the deleted block |
| `--arm BP` | 43 | homology arm; primers excluded from this much on each side |
| `--min-clearance BP` | 20 | warn below this distance from the excluded zone |
| `--span BP` | 1000 | how far beyond the excluded zone to search |
| `--tm C` | 60 | target melting temperature (nearest-neighbour, SantaLucia 1998) |
| `--tm-range LO,HI` | 58,62 | acceptable Tm window |
| `--max-tm-diff C` | 3 | largest Tm difference within a pair |
| `--length LO,HI` | 20,30 | acceptable primer length |
| `--gc LO,HI` | 0.4,0.6 | acceptable GC fraction |
| `--specificity` | auto | `auto`, `gggenome`, `fasta`, `local`, `none` |
| `--gggenome-db DB` | by species | GGGenome database name |
| `--genome-fasta PATH` | | local assembly for an offline genome-wide count |
| `--mismatches N` | 2 | mismatches allowed when counting off-target sites |
| `-o, --outdir DIR` | print only | where to write the TSV and the maps |
| `--alternatives N` | 2 | runner-up pairs to report |
| `--no-maps` | | skip the GenBank files |
| `--min-band-difference BP` | 50 | below this the bands are called unresolvable |
| `-v, --verbose` | | log progress |

Exit codes: `0` a design was produced, `1` the run could not be set up (bad
coordinates, unreadable file, no database for the species), `2` the run was
valid but no pair was found.

Beyond the options, candidates are also rejected for a run of 5 identical
bases, a 3' end that pairs over more than 4 bases with the other primer, an
absent or excessive GC clamp (1–3 G/C in the last five bases), and any overlap
with an ambiguity code.

## As a library

```python
from genotyping_primers import DesignConfig, design, prepare, resolve_span, target_from_file

target, notes = target_from_file("locus.fa")
target, deletion = prepare(target, *resolve_span(target, "5501-6000"))
result = design(target, deletion, DesignConfig(homology_arm=43, specificity="local"))

best = result.best
print(best.forward.sequence, best.reverse.sequence)
print(best.wt_product, best.ko_product, best.warnings)
```

`design` never raises for a design that merely turned out badly — it returns
the pairs it found with their warnings, and the notes explaining anything it
could not do. It raises only when the run could not be set up at all.

## What it does not do

It picks primers; it does not design guide RNAs, and it has no opinion on where
your cut sites should be. It predicts annealing, not amplification: a pair that
passes every filter here can still fail on a GC-rich template or a secondary
structure it cannot see. And the specificity check counts sites — it does not
know which of them a polymerase will actually extend. For anything expensive
downstream, run the final pair past your usual in-silico PCR check as well.

## Tests

```bash
pip install -e ".[dev]"
pytest
```

The suite is entirely offline. It builds synthetic sequences from a fixed seed
so the expected coordinates are exact, and it re-measures every reported
amplicon from the sequences the design claims to amplify rather than trusting
the tool's own arithmetic. It also covers the coordinate conversions (local,
genomic and rotated-circular), the GenBank round trip, the chunk-boundary
handling of the FASTA scanner, and the GGGenome client against stubbed
responses — including the hg38 fallback.

## License

MIT. See [LICENSE](LICENSE).
