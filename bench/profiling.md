# Profiling notes (Phase 0, step 0.4)

The runner answers *how much* slower Biopython is.  This answers *why*, and
therefore which parts of it a C core can actually delete.

Reproduce:

```sh
.venv/bin/python bench/profile_paths.py --all --repeat 1 --top 8
.venv/bin/py-spy record --format raw -o /tmp/prof.raw --rate 200 \
    -- .venv/bin/python bench/profile_paths.py --workload fastq:biopython-seqio
.venv/bin/python bench/profile_paths.py --pyspy /tmp/prof.raw --pyspy-mark consume
```

`profile_paths.py` builds its paths from the **runner's own workload
factories**, so the loop being profiled is the loop that was measured and the
path names are `<group>:<implementation>` — the same names as in
`bench/results/latest.json`.

## How to read these numbers

**Shares, not seconds.**  `cProfile` instruments every Python call.  On
call-heavy code that inflates wall time roughly 3x (`fastq:biopython-seqio`:
10.8 s profiled against 3.19 s measured).  A profiled second is not a real
second; a profiled *percentage* is still a fair statement of where the time
goes, and no number here should be compared against `bench/targets.md`.

**"Biopython frames" is the decisive column.**  It is the share of time in
`site-packages/Bio/**`.  For a pure-Python library that is interpreter work,
which a C implementation removes outright rather than merely matching.  For
pyfaidx, zlib and pysam it is ~0%, and that is not a compliment to them: they
already keep their work in C.  The point of the column is to separate *work to
delete* from *work to match*.

**Driver frames are counted separately.**  The benchmark's own loop
(`run.py`'s `consume`/`digest`) does identical work in every row; leaving it
inside the totals would flatter whichever implementation is slowest, because a
fixed overhead is a smaller share of a bigger number.

**Where cProfile is the wrong instrument.**  It cannot see inside a C extension
and charges that time to whichever Python frame called it.  For pysam the
hotspot table therefore shows `run.py:factory` at 65% and says nothing useful;
those rows are read with `py-spy`, which samples the real process and reports
wall-clock self time with no instrumentation at all.

## The table

`Biopython frames` is the cProfile share; `measured` is the same
implementation's median from the full run stored in `bench/results/latest.json`,
given only to show how far instrumentation inflates the profiled seconds.  The
two columns come from separate runs, so the ratio between them carries a few
percent of run-to-run noise on top of the instrumentation factor.

| Path | measured | profiled | Biopython frames | headline hotspot |
|---|---:|---:|---:|---|
| `fastq:biopython-seqio` | 3.188 s | 10.818 s | **55.7%** | `QualityIO.__next__` 37.4% |
| `fastq:biopython-lowlevel` | 0.639 s | 2.532 s | **52.0%** | `FastqGeneralIterator` 52.0% |
| `fastq:naive-python` | 0.448 s | 2.195 s | 0.0% | its own loop 47.7% |
| `fastq:pysam-htslib` | 0.771 s | 1.528 s | 0.0% | the C call 64.7% |
| `fastq-gz:biopython-seqio` | 5.316 s | 15.058 s | **48.4%** | `__next__` 33.9%, zlib 11.2% |
| `fastq-gz:biopython-lowlevel` | 2.551 s | 5.759 s | **36.1%** | FGCI 36.1%, zlib 27.7% |
| `fastq-gz:naive-python` | 2.631 s | 6.730 s | 0.0% | zlib 23.9%, `gzip.readline` 19.1% |
| `fastq-gz:pysam-htslib` | 2.145 s | 3.027 s | 0.0% | the C call 81.5% |
| `fastq-gz:zlib-raw-only` | 1.401 s | 1.460 s | 0.0% | `zlib.decompress` 97.6% |
| `fasta:biopython-seqio` | 0.303 s | 0.591 s | **49.5%** | `FastaIO.__next__` 49.5% |
| `fasta:biopython-lowlevel` | 0.273 s | 1.143 s | **55.4%** | `SimpleFastaParser` 55.4% |
| `fasta:naive-python` | 0.342 s | 1.251 s | 0.0% | its own loop 56.7% |
| `fasta:pyfaidx-sequential` | 0.153 s | 0.152 s | 0.0% | `str.replace` 67.1% |
| `fasta-random:biopython-index-random` | 15.956 s | 39.920 s | **49.2%** | `_index.get_raw` 27.9% |
| `fasta-random:pyfaidx-random` | 0.0003 s | 0.001 s | 0.0% | nothing measurable |
| `ops:biopython-revcomp` | 0.280 s | 1.349 s | **63.9%** | `SeqRecord.reverse_complement` 32.9% |
| `ops:biopython-gc` | 0.697 s | 3.389 s | **37.5%** | `Seq.count` 26.0%, ABC 39% |
| `ops:naive-strtranslate` | 0.041 s | 0.140 s | 15.8% | `str.translate` 18.6% |
| `ops:naive-strcount` | 0.107 s | 0.242 s | 9.6% | `str.count` 36.4% |

![Deletable share per implementation: Python bytecode in Bio/** against rows
that are already C-bound](../misc/profile-share.png)

(chart rendered by `bench/plot_profile.py`; the table above stays the data source)

## What the hotspots say

**`fastq:biopython-seqio` — 55.7% deletable.**  `QualityIO.__next__` alone is
37.4%: a Python state machine reading lines, validating characters, slicing and
re-wrapping.  `SeqRecord._from_validated` (9.1%) and `array.tolist` (4.9%) are
object construction on top of it.  `str.isprintable` at 3.2% is worth naming:
Biopython checks every quality character, per record, in Python.

The low-level `FastqGeneralIterator` is *also* 52% interpreter time — so
"Biopython without `SeqRecord`" is not an escape from Python, only from object
construction.  That is why it is 5x faster than `SeqIO.parse` (0.639 s vs
3.188 s) while still leaving a pure-Python parser (0.448 s) ahead of it.

**`fasta-random:biopython-index-random` — 49.2% deletable, and mostly re-work.**
`_index.get_raw` (27.9%) reads the record's bytes back off disk, and
`FastaIO.__next__` (16.5%) then re-parses it with a regex (`re.Pattern.match`,
14.9%), accumulating lines into a list (11.5%) and translating/joining bytes
(4.3% + 4.1%).  None of this would exist if the index stored a byte offset and
a line layout instead of only whole-record extents.

**`ops:biopython-gc` — the surprising one.**  `Seq.count` is 26%, but type
dispatch is another **39%**: `isinstance` 18.0%, `abc.__instancecheck__` 11.5%
and `_abc._abc_instancecheck` 9.5% — three self-times that add up, because the
ABC fast path is built out of all three frames.  More than a third of a sequence
operation is spent deciding what type its argument is.

**`ops:biopython-revcomp` — 63.9%, the highest share measured.**  The operation
itself (`Seq.reverse_complement`) is only 4.5%; the cost is
`SeqRecord.reverse_complement` (32.9%) walking and rebuilding
`letter_annotations` (6.6%), with `__setitem__` (5.1%) and `isinstance` (7.3%)
along the way.

**`fasta:pyfaidx-sequential` — 0% Biopython frames and still slow-ish.**  Its
67.1% is in `str.replace`, which is C but is doing newline removal by
substitution.  A `mmap`ed reader that never builds a string would not pay it.

## py-spy cross-check

Sampled with no instrumentation, marked to the measured loop, so these are
wall-clock self times rather than cProfile's attribution.

**`fastq:biopython-seqio`** — the time is spread across the FASTQ state machine
with no single dominant line: 9.6%, 6.2%, 5.2%, 3.8%, 3.6% on five different
lines of `Bio/SeqIO/QualityIO.py:__next__`, plus `Seq.__len__` 3.5% and the
harness factory 4.1%.  This is the signature of interpreter overhead rather than
one slow operation, and it is why a C parser can win big here: there is no
algorithm to fix, only work to remove.

**`ops:biopython-gc`** — confirms the cProfile picture and sharpens it: **26.5%**
of samples on the two `abc.__instancecheck__` lines and **56.5%** across six
lines of `Bio/Seq.py:count`, with `SeqUtils/__init__.py` generator expressions at
8.9%.  (py-spy's type-dispatch share is lower than cProfile's 39% because the
sampler charges the C fast path beneath `__instancecheck__` to the frame that is
executing, not to the three separate frames cProfile times.)  Either way, a
large fraction of `gc_fraction` is ABC type dispatch.

**`fastq:pysam-htslib`** — 38.5% + 24.4% + 12.7% of samples sit on three lines
of the harness loop that call into htslib, and `full_header` takes another
11.5%.  Zero samples land in Python inside pysam itself: this row is C-bound,
and the 11.5% is harness-mandated header reconstruction (needed to compare
headers like with like), so pysam is penalised slightly here rather than
flattered.

## Incidental findings worth keeping

- **`gzip.GzipFile` overhead is visible.**  `gzip.py:closed` costs 3.3% on the
  SeqIO gz row, 7.4% on the low-level row, 4.3% on the naive row, and
  `_compression.py:_check_not_closed` another 8.4% on the naive row — every
  `readline()` re-checks a `closed` property that does an `isinstance`.  For a
  parser that owns its own buffered reader and `libdeflate` calls, this
  disappears.
- **Biopython 1.88 does not auto-decompress `.gz`.**  `SeqIO.parse("reads.fastq.gz",
  "fastq")` raises `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b`,
  so the natural one-liner fails on the format most real FASTQ ships in; the
  caller must wrap `gzip.open` themselves.  Measured, not assumed — the runner
  opens the handle explicitly for this reason.

## What this means for Phase 1

The gaps in `bench/targets.md` are **deletable, not merely matchable**.  On
every Biopython row between a third and two thirds of the time is Python
bytecode inside `Bio/**`, and the py-spy samples agree: the time is spread
across interpreter-level work rather than concentrated in one algorithm.  The
exceptions are the gzip decompression floor and pysam, both already C-bound —
for those the goal is parity plus `libdeflate`, not a rewrite of anything.
