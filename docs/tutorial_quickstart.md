# Tutorial 1 — Quickstart: the first ten minutes

Install, then read two files, write two files, and call it a session. This
tutorial does not measure anything; it builds the confidence the other
tutorials spend. Every block that shows a pipeline is copied from an example
script in `examples/`, which the test suite runs — the block and the script
cannot drift without a red test.

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the owner
commits and pushes.

## Install

```console
$ pip install .            # the bare package
$ pip install .[arrow]     # + pyarrow, for the record tables
$ pip install .[alignment] # + parasail, for the aligner
```

Extras are optional and independent of each other. Everything this tutorial
and most of the others use — the parsers, seq operations, restriction layer,
protein measurements, flat-file readers, writers — is in the bare package.

## The build identity

The package knows what it was compiled with, and the answer is a dict:

```python
import biofasting

info = biofasting.build_info()
```

which prints, on this repository's machine:

```
    libdeflate = 1.23 (vendored, static)
    nanobind = 3.1.0
    relax_min_size = yes (-O3, not -Os)
    stack_protector = yes (-fstack-protector-strong)
    version = 0.0.1
    python = 3.13.5
    platform = Linux-6.12.111+deb13-amd64-x86_64-with-glibc2.41
    machine = x86_64
```

If anything about the running build is in question, this is the answer, not a
guess — `machine` is the value the CPU-feature dispatch is keyed on.

## Reading a FASTQ

`biofasting.open_fastq(path)` yields `(title, sequence, quality)` triples of
`bytes` — the same shape `Bio.SeqIO.QualityIO.FastqGeneralIterator` yields,
with the `@` removed from the title. Gzip is *sniffed from the magic bytes*,
not from the file name, so a `.gz` read exactly like a plain file:

```python
for title, sequence, quality in biofasting.open_fastq(path):
    ...
```

The first record of the quick corpus, straight:

```python
records = biofasting.open_fastq("bench/data/fastq/reads_10k.fastq")
first = next(records)
```

```
A title   = b'SRR000001.0 HWI-ST1:1:FCBAR:1:1101:1000:'
A quality len = 150 sequence len = 150
```

(The `A` lines are the print statements a demo would add.) Three things worth
knowing before the second file:

* the quality is the **ASCII string the file carries** (phred+33), not a list
  of integers — `sum(quality) - 33 * len(quality)` is the usual per-record
  arithmetic;
* the file is memory-mapped, so opening costs microseconds and no resident
  memory beyond the pages actually touched; there is nothing to close;
* a malformed record raises `ValueError` with the record number and the byte
  offset where it was found — see the parsing tutorial for the exact refusal
  rules.

## Reading a FASTA

FASTA reads build an index once and answer from it:

```python
index = biofasting.open_fasta("bench/data/fasta/genome_1mb.fasta")
name = next(iter(index))

index.title(name)          # the whole header line, without '>'
index[name]                # the sequence, as bytes
index.sequence_slice(name, 1000, 1050)   # just the window asked for
index.sequence_length(name)              # free
```

```
B name = chrS
B title = chrS
B first record length = 1000000
B slice 1000:1050 = CCAATTACTGTCTCCCGCGAAAAATTAGATAACATGTCTCAAATAGAGTC
B sequence_length = 1000000
```

Unlike `Bio.SeqIO.index`, a fetch is a copy out of the mapping — it does not
re-parse the record. The same file's per-access cost is the number the
parsing tutorial cites to `bench/targets.md`'s first-pass ranking.

## Writing

Writers go the other way, byte-identical to `SeqIO.write`, and take the same
shapes the readers yield:

```python
readable = biofasting.read_fasta("bench/data/fasta/genome_1mb.fasta")
window = readable[0][2][:6000]
biofasting.write_fasta([("chrS window", window)], destination, wrap=60)
```

```
W wrote snippet.fasta | read back: 1 record(s), 6000 bp; head:
CGCGGTGCTTATGTCGACGTGCAACTTGAACGGCGACCGTTATTATAACGGATTGAAATC
```

`write_fasta` takes `(title, sequence)` pairs, wraps at 60 by default, and
writes the whole file in one kernel buffer. `write_fastq(records, handle)`
takes the reader's `(title, sequence, quality)` triples and does the same
thing for FASTQ. The full pipeline that combines reading, filtering and
writing is `examples/fastq_qc.py` — next tutorial.

## Measuring the records

The seq operations are module functions, and they take what the readers
yield (`bytes`):

```python
biofasting.gc_fraction(sequence)          # fraction, by default
biofasting.gc_fraction(sequence, ambiguous="ignore")
biofasting.reverse_complement(sequence)
biofasting.count_kmers(sequence, k)       # overlapping k-mers, no stride
```

The protein side is `ProteinAnalysis` (twelve measurements, one constructor)
and `IsoelectricPoint`; the restriction layer is a real port of
`Bio.Restriction`'s data and engine. These have their own tutorials:
[parsing](tutorial_parsing.md),
[flat files](tutorial_flatfiles.md),
[restriction](tutorial_restriction.md),
[protein](tutorial_protein.md),
[alignment](tutorial_alignment.md).

## Where the examples come in

Six scripts in `examples/` are whole pipelines — QC filtering, a plasmid
digest, a protein profile, codon adaptation, primer scoring, the polars
hand-off — each with argparse knobs and each asserting in
`tests/test_examples.py` what it prints. `examples/fastq_qc.py` runs like:

```console
$ python3 examples/fastq_qc.py
read 10,000 records from .../bench/data/fastq/reads_10k.fastq
kept 1,644 records (16.4%); dropped: N 0, quality 5,797, GC 2,559
kept reads: 246,600 bp, mean length 150.0
wrote .../some/temp/dir/qc_passed.fastq
```

The kept reads land in a fresh temp directory by default, so an example run
never dirties the repository; `--out` points it anywhere else.