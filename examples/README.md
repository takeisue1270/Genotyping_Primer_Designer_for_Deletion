# Examples

Two synthetic sequences, generated from a fixed random seed. They carry no real
biology — they exist so the commands in the README can be run and checked
without a network connection or a genome on disk.

| File | What it is |
| --- | --- |
| `example_locus.fa` | 12 kb linear target |
| `example_plasmid.fa` | 5.2 kb circular construct |

## A linear deletion

```bash
genotyping-primers --sequence examples/example_locus.fa --delete 5501-6000
```

Expect a pair around 1,446 bp / 946 bp. Add `--insert GGATCCAAGCTTGA` and the
edited band becomes 960 bp, because 14 of the 500 deleted bases are put back.

## A deletion that crosses the origin

```bash
genotyping-primers --sequence examples/example_plasmid.fa --circular --delete 4950-250
```

The deletion runs from 4,950 through the origin to 250. The forward primer
lands near 4,574 and the reverse near 800 — on opposite sides of base 1 — and
the report keeps using the plasmid's own numbering throughout.

Note that the GenBank maps written with `-o` are rotated so the deletion sits
in the middle; their `DEFINITION` line says which base of the original sequence
became base 1.

## Reproducing the files

```python
import random, textwrap

random.seed(20260926)
def seq(n, gc=0.48):
    out = []
    for _ in range(n):
        r = random.random()
        out.append("G" if r < gc / 2 else "C" if r < gc else "A" if r < gc + (1 - gc) / 2 else "T")
    return "".join(out)

locus, plasmid = seq(12000), seq(5200)   # generated in that order
```
