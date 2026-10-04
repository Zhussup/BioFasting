# Tutorial 2 — Parsing: iteration, indexes, grids, tables, shims

The parsing surface is the package's core competency, and this tutorial is
the one that explains its four access patterns, the compression and malformed
file rules, the writers, the `SeqRecord` shims and the Arrow tables.  Facts
about speed are *stated and cited* — the measured numbers and their gates
live in `bench/targets.md` (the phase‑1 ranking passes, the writer pass, and
the Arrow pass), never re-quoted here as tables.

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the
owner commits and pushes.

## Four ways in

| what you want | what to call |
|---|---|
| one record after another | `biofasting.open_fastq(path)` / `biofasting.open_fasta(path)` |
| records by name, cheaply, many times | `open_fastq_index(path)` / `open_fasta(path)` |
| whole file as one zero-copy array | `open_fastq_grid(path)` / `open_fasta_grid(path)` |
| the whole file for a data frame | `read_fastq_table(path)` / `read_fasta_table(path)` |

All four map the file instead of reading it; the differences are what each
one *keeps*.

## Iteration

`open_fastq` yields `(title, sequence, quality)` triples of `bytes`.  Gzip is
sniffed from the magic bytes, not the file name — a `.fastq.gz` behaves like a
plain file, and the decompression is libdeflate, whole-file then scan (the
bench passes say what that buys against Biopython's `gzip.open`+parse).  The
parsing tutorial's flagship pipeline (`examples/fastq_qc.py`) reads exactly
this way:

```python
for title, sequence, quality in biofasting.open_fastq(args.path):
    total += 1
    if sequence.count(b"N"):
        failed_n += 1
        continue
    gc = biofasting.gc_fraction(sequence) * 100.0
    if not args.gc_low <= gc <= args.gc_high:
        failed_gc += 1
        continue
    mean_q = sum(quality) / len(quality) - 33.0
    if mean_q < args.min_q:
        failed_q += 1
        continue
    kept += 1
    bases += len(sequence)
    kept_records.append((title, sequence, quality))
biofasting.write_fastq(kept_records, args.out)
```

(copied from `examples/fastq_qc.py`.  Its report, on the quick corpus:

```
kept 1,644 records (16.4%); dropped: N 0, quality 5,797, GC 2,559
kept reads: 246,600 bp, mean length 150.0
```

— a line `tests/test_examples.py` asserts.)

FASTA iteration is the same shape, but the file *is* the index: `open_fasta`
returns a mapping you can iterate keys over, and `read_fasta(path)` is the
convenience that materialises `(name, title, sequence)` in one call — for
tests and small files, per its own docstring.

## Malformed files: the refusal rules

The scanner refuses, never warns-and-continues.  A truncated record raises
with the message and the position:

```python
bad = path with "@good ... @truncated\nAAG\n+\nII\n"
list(biofasting.open_fastq(bad))
```

```
Lengths of sequence and quality values differs for truncated (3 and 2). (record 1, byte offset 18)
```

That is the reference's own message, plus the record number and byte offset
in parentheses.  The whole grammar is the reference's `FastqGeneralIterator`
— the edge-file corpus pins both libraries accepting the same set of files,
including the blank-lines case Biopython's *index* (not its parser) rejects.

## Indexes: fetch by name

`open_fastq_index(path)` behaves like `Bio.SeqIO.index(path, "fastq")` in
what it answers, and like nothing in Biopython in what it costs: the index
stores where each field *is*, so a fetch is a memcpy and
`quality_length(name)` does not touch the file.

```python
idx = biofasting.open_fastq_index("bench/data/fastq/reads_10k.fastq")
name = idx.keys()[0]          # 'SRR000001.0' — the header's first word
one = idx[name]               # the sequence, as bytes
idx.sequence_slice(name, 0, 12)
idx.quality_length(name)
idx.quality(name)
```

```
N type = bytes | value head = b'CCACTAGGACAGTCTCAGTGTGTAGT'
N sequence_length = 150 | quality_length = 150
N sequence_slice 0:12 = b'CCACTAGGACAG'
```

Keys keep file order; a duplicate name raises `ValueError`; a missing name
raises `KeyError`.  A gzipped file is accepted (BGZF included) and indexed in
its inflated form — the docstring records what that costs in memory and buys
per fetch, versus `SeqIO.index`'s 274 µs BGZF fetch measured in
`bench/targets.md`.

The FASTA index (`open_fasta`) is the same mapping shape: `index[name]`,
`index.title(name)`, `index.sequence_slice(...)`, `sequence_length(...)`.
There is no FASTA analogue of *re-parsing on fetch*, and that difference is
the first-pass ranking's headline: choosing 150 bp out of the middle of the
million-base contig is one gather, not a re-parse of the whole record.

## Grids: zero-copy, one array

`open_fastq_grid(path)` makes no Python objects at all: `(records, length)`
`uint8` arrays that are *views of the file's bytes*, strided.  It requires
every record to be one shape — checked by reading, not assumed — and a
record that wraps or a file whose records differ raises `ValueError` with
the rule that failed:

```python
biofasting.open_fastq_grid("bench/data/fastq/reads_10k.fastq")
```

```
ValueError: record 1 does not have the same shape as the first one
```

— and that refusal is the *right* answer here, because the quick corpus is
not uniform: 15% of its reads are 3'-adapter-trimmed, so shorter reads sit
between the 150-base ones (this is a corpus property the generator
documents, and a fact the example above leans on when it says "mean length
150.0" only for the reads that *kept*).  A uniform file has a grid:

```python
import pathlib, tempfile

path = pathlib.Path(tempfile.mkdtemp()) / "uniform.fastq"
path.write_text("".join(f"@read{i}\nACGTACGT\n+\nIIIIIIII\n" for i in range(4)))
grid = biofasting.open_fastq_grid(path)
```

```
uniform demo: (4, 8) (4, 8) title_length = 6
```

`grid.sequences[:, :4]` is an index like numpy slicing — nothing is copied
until you ask.  `open_fasta_grid(path, name)` does the same for one FASTA
record whose lines are one width.

## The writers

`write_fasta(records, handle, wrap=60)`, `write_fastq(records, handle)` and
`write_qual(records, handle)` take their records in the readers' shapes —
`(title, sequence)` pairs, `(title, sequence, quality)` triples, `(title,
quality)` pairs — and are byte-identical to `SeqIO.write` under a gate.  The
kernel builds the whole file in one buffer and writes once.  The measured
multiples against the reference's writer are in `bench/targets.md`'s writer
pass (the quality writer's per-base `%i`-format cost is a stated finding
there — the worst of the three reference writers).

## The SeqRecord shims

`biofasting.interop` is where Biopython objects meet this package:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
shim = biofasting.from_seqrecord(SeqRecord(Seq("ATTCCGGACCGC"), id="demo"))
```

```
J from_seqrecord = tuple | title demo <unknown description> | seq ATTCCGGACCGC
```

`(title, sequence)` or `(title, sequence, quality)` — a triple or a pair
exactly as the reader/writer shapes want.  The other direction:

```python
sr = biofasting.to_seqrecord("SRR000001 title here", b"ATTCCGG", quality=None)
```

```
J to_seqrecord: SeqRecord SRR000001 | SRR000001 title here | ATTCCGG
```

And for whole files, `fastq_seqrecords(path)` and `fasta_seqrecords(path)`
give `Bio.SeqIO.parse`-shaped streams with our scanner underneath — "the same
objects, one at a time, without Biopython parsing the file again".

## The Arrow tables

`read_fastq_table(path)` lays the whole file into Arrow's own buffers —
offsets and one concatenated bytes buffer — in one pass from C++, and the
result is a real `pyarrow.Table` whose columns are `polars-bio`'s columns:
`name, description, sequence, quality`.  (`description` is the title *after*
the key with the whitespace trimmed — `>rec1 a description` gives
`a description`; that is `polars-bio`'s reading, deliberately, and not
`SeqRecord.description`'s.)

`examples/reads_to_polars.py` is the pipeline: read, hand to polars,
aggregate in polars.  The final frame it prints on the quick corpus:

```
│ records ┆ bp total ┆ mean length ┆ pooled GC ┆ min GC ┆ max GC │
│ 10000   ┆ 1500000  ┆ 150.0       ┆ 0.446777  ┆ 45     ┆ 88     │
```

and its last line names the table the frame stands on:

```
the table still holds the bytes the frame reads from (10000 rows, 4 columns,
schema name, description, sequence, quality)
```

The hand-off's price is measured and documented in `bench/targets.md`'s
Arrow pass: `pl.from_arrow` copies no strings but is O(rows) building
polars' own index — a caller who stays in pyarrow does not pay it, and the
frame pays once, not per column.

## What is deliberately absent

* No `SeqIO.parse`-style format guessing: the format is the function you
  called.  Flat files (GenBank/EMBL/SwissProt) have their own
  [tutorial](tutorial_flatfiles.md).
* No compressed *FASTA* indexing: `open_fasta` refuses a gzip file with a
  `ValueError` naming the alternative, rather than silently indexing
  nothing.
* No automatic `str`/`bytes` re-encoding: sequences are `bytes` end to end,
  and `str` is accepted where the API says so.