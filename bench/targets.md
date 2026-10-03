# Ranked targets (Phase 0, step 0.5)

The point of Phase 0 is a list of targets ordered by **time actually lost**, so
Phase 1 implements what pays, in the order it pays.  Every number here comes
from one reproducible run; nothing is estimated.

Reproduce:

```sh
python3 bench/gen_data.py                 # byte-reproducible corpus
.venv/bin/python bench/run.py --repeat 5  # -> bench/results/latest.json
```

Host: Intel i5-13420H, single-threaded, Debian, Biopython 1.88, pysam 0.24.1,
pyfaidx 0.9.0.4.  All rows below are correctness-gated (see `bench/README.md`):
an implementation whose output does not match Biopython's is printed `INVALID`
and cannot appear in this table.  Every row here validated.

**Run-to-run variance.**  Two full runs of the identical command, same host,
same corpus: `fastq-plain` 3.176 s then 3.188 s (+0.4%), `fasta-random`
150.9 then 159.6 ms/access (+5.8%).  Sub-second rows move by a few percent
between runs; that is normal for a single-threaded laptop measurement, and no
claim here rests on a difference smaller than that.  The numbers below are from
the run currently stored in `bench/results/latest.json`; a re-run on your machine
will land in the same place but not on the same digit.

## The ranking

`Gap` is Biopython's time divided by the best alternative's.  `Stolen` is the
absolute seconds a perfect implementation would remove at this workload size,
which is what says how much a real dataset would save.

| # | Workload | Biopython | Best alternative | Gap | Stolen |
|---|---|---|---|---:|---:|
| 1 | `fastq-plain` — parse 1M × 150 bp, 369 MB | 3.188 s (313,700 rec/s) | naive pure-Python 0.448 s | **7.1×** | 2.74 s |
| 2 | `fasta-random` — 100 × 150 bp slices | 15.956 s (**159.6 ms/access**) | pyfaidx 0.003 ms/access | **~50,000×** | 15.96 s |
| 3 | `fastq-gzip` — same reads, gzipped | 5.316 s | zlib floor 1.401 s | 3.8× | 3.92 s |
| 4 | `ops-gc` — GC fraction over 100k reads | 0.697 s | `bytes.count` 0.107 s | 6.5× | 0.59 s |
| 5 | `ops-revcomp` — reverse complement, 100k reads | 0.280 s | `bytes.translate` 0.041 s | 6.8× | 0.24 s |
| 6 | `fasta-scan` — read 100 Mbp in 5 records | 0.303 s | pyfaidx sequential 0.153 s | 2.0× | 0.15 s |

![The ranking as a chart: each row is Biopython's time vs. the best
alternative's, log scale; the columns to the right are Gap and
Stolen](../misc/g1-ranking.png)

(chart rendered by `bench/plot_g1.py`; the table above stays the data source)

Two rows pull in opposite directions and both matter:

- **Rank 1 is the scalable one.** 1M reads is a small file; production runs are
  10M–1B reads, so ~2.7 s per million is 27 s–45 min per file, and FASTQ parsing
  is on the critical path of nearly every pipeline.  This is where absolute time
  is recovered.
- **Rank 2 is the absurd one.** ~160 ms to fetch 150 bases is not a tuning
  problem, it is a design gap: `SeqIO.index` stores whole-record offsets, so a
  slice costs a full record parse (see `bench/README.md`).  The fix is a real
  `.fai`-style index, which is a small, well-understood piece of work with a
  five-figure speedup at the end of it.

## What each "best alternative" actually is

The gap column is only meaningful with the caveat attached to each winner:

- **naive pure-Python (0.448 s)** — a floor, not a product.  It handles
  single-line FASTQ only: no wrapped sequences, no multi-line records, no
  malformed-input policy, and it does not build record objects at all.  It is
  quoted because it proves the speed is available without C, which means a C
  core should clear it comfortably.  Note it also beats `pysam` (0.771 s) on
  plain FASTQ, because pysam builds an object per record.
- **pyfaidx (0.003 ms/access)** — an equal-work comparison: both rows return the
  same 150 bp and pass the same digest.  pyfaidx wins because it seeks to a byte
  offset from its `.fai`; Biopython has no equivalent structure to seek with.
- **zlib floor (1.401 s)** — decompression only, no parsing.  This is the part
  of the gzip row that no parser can remove; `pysam` (2.145 s) is the honest
  competing parse, only 2.5× behind Biopython because htslib decompresses in C.
  A parser with `libdeflate` should land between the floor and pysam.
- **`bytes.count` / `bytes.translate` (0.107 s / 0.041 s)** — the same caveat as
  the naive parser: they do not implement Biopython's semantics (ambiguous-base
  handling, `N` policies, wrapped input).  They establish the floor a real
  implementation has to approach, not a number to quote as achievable.

## Where the time goes (step 0.4)

Profiling says the gap is deletable rather than a hard floor — see
`bench/profiling.md`.  In short: `Bio.SeqIO.parse` on FASTQ spends **55.7%** of
its time in library Python bytecode, the GC op spends **26.5%** of its py-spy
samples in `abc.__instancecheck__` alone, and only the gzip floor and pysam are
genuinely C-bound already.

## Conclusion for Phase 1 (Gate G1)

**Gate G1 passes: FASTQ/FASTA is confirmed as the flagship from data, not from
intuition.**  The two largest gaps in the table are both in that scope, and they
are large for different reasons — FASTQ parsing because Biopython does it in
Python, FASTA random access because Biopython has no slice-capable index.

Ordering inside Phase 1:

1. **FASTQ parse, plain** (rank 1) — largest recoverable absolute time, largest
   interpreter share, widest audience.  Step 1.2.
2. **FASTA indexed access** (rank 2) — largest ratio, smallest implementation
   surface.  Step 1.2, after the parser kernel.
3. **FASTQ over gzip** (rank 3) — bounded by the 1.401 s decompression floor;
   the win here is `libdeflate`, not parser cleverness, so it is a packaging
   decision more than a parsing one.
4. **Sequence ops** (ranks 4–5) — real but small per call; they become worth
   doing when the parsed data is already in a zero-copy buffer (step 1.3),
   which is the only way `bytes.count`-class speeds reach object-shaped results.

`fasta-scan` (rank 6) is measured and deliberately **not** prioritised: 0.15 s
per 100 Mbp does not justify a kernel on its own.

## Phase 1 result (step 1.6, 2026-10-03)

The ranking above is what Phase 0 asked for; this is what Phase 1 delivered
against it. Same corpus, same machine, same harness (`bench/run.py`), median of
5 after one untimed warmup, every row digest-gated against Biopython 1.88 before
its time is quoted. Reported against the build that produced them:
`GNU 14.2.0`, C++20, Release, `-O3` (not `-Os`), `seqops_level=avx2` on an
Intel i5-13420H — see `bench/results/m8-full.json` for the full identity.

| Workload | Biopython 1.88 | Best alternative then | BioFasting | vs `SeqIO` | vs alt |
|---|---:|---:|---:|---:|---:|
| FASTQ parse, 1M×150 bp plain | 3.164 s | pysam 0.683 s | **0.224 s** | **14.1×** | 3.05× |
| FASTQ parse, 1M×150 bp gzipped | 5.242 s | pysam 2.122 s | **0.864 s** | **6.07×** | 2.46× |
| FASTA scan, 100 Mbp | 0.284 s | pyfaidx 0.138 s | **0.070 s** | **4.06×** | 1.97× |
| FASTA random 150 bp slice | 146.290 ms | pyfaidx 2.572 µs | **0.167 µs** | **~876,000×** | 15.4× |
| reverse complement, 100k reads | 0.290 s | `str.translate` 0.041 s | **0.034 s** | **8.5×** | 1.21× |
| GC fraction, 100k reads | 0.679 s | `str.count`×2 0.105 s | **0.044 s** | **15.4×** | 2.39× |

**Gate G1's exit condition is met:** ≥10× vs `SeqIO` on plain FASTQ (14.1×) and
≥1× vs pysam (3.05×). The gzipped row is the interesting one: 0.864 s for
inflating *and* parsing the file is 0.60× the 1.446 s that `gzip`-class
decompression alone costs, because `libdeflate` is not `zlib`. The whole
pipeline fits inside what the reference-only decompression rows spend before a
single record is looked at.

The zero-copy views cannot be measured by this harness — its model of an
implementation is "yield each record, then do the least possible work with it",
and a view never touches the bytes, so a row for it would win by doing nothing.
`bench/grid.py` measures them separately, and gates the array's *content*
against `SeqIO.parse` before quoting anything:

| What | Cost |
|---|---:|
| make the 102 MB FASTA addressable | 38.1 ms vs `SeqIO.index` 636.2 ms (**16.7×**) |
| get 40 Mbp of chr1 as an array, index held | **0.61 µs** vs `SeqIO.index` 15.19 ms (**~24,900×**), pyfaidx 68.3 ms |
| RSS paid to hold that array | **0 KiB** vs 39,064 KiB — a view is the page cache's copy, not ours |
| GC count over chr1 | 26.8 ms vs `SeqIO`'s `gc_fraction` 468.8 ms (**17.5×**), `str.count` 164.8 ms (6.1×) |

Two rows of that table were wrong when first written and are recorded here
because the correction *is* the finding: the array benchmark initially built our
index inside the timed region while `SeqIO.index` had been built at setup, which
reported the copying implementation as twice as fast; and the GC row compared a
count over the grid's whole-lines-only array against a count over the full
chromosome, which flagged a correct answer as INVALID. Both are fixed, and the
fix in both cases was to compare the same thing on both sides.

The FASTQ grid does not appear above, and the reason is measured rather than
asserted: the corpus reads have **8 distinct header widths (58–65 bytes) across
1,000,000 records**, so their records do not share a stride and no
`(records, length)` array can describe them. `open_fastq_grid` refuses, naming
record 1. Padding the headers would make it a grid; that is a corpus change, not
a kernel one, and it would be a different file.

## Phase 2 ranking (2026-10-03)

Phase 1 spent the ranking above.  This is the same exercise applied to what
comes next, run *before* any of it is built, so that "what Phase 2 works on" is
a measurement rather than a preference.

Translation was added to the runner as a fourth in-memory workload
(`ops-translate`, `make_translate_impls` in `bench/run.py`): every
implementation translates the same already-loaded reads and yields the protein
strings, and the harness digests them against `Seq.translate` before it times
anything.  Median of 5, plain `--quick` corpus, same build identity as Phase 1
(`GNU 14.2.0`, `-O3`, `seqops_level=avx2`).

The row below was measured, as the rule requires, *before* the kernel was
written -- over 10,000 reads, a tenth of what the runner now uses, because the
question at that moment was whether the gap was worth a kernel at all and a
tenth of the work answers it just as well.  The per-read cost is the same at
either size (13.6 µs then, 13.4 µs at 100,000), so the ranking stands as
recorded.

| # | Workload | Biopython 1.88 | Best alternative | Gap | Stolen |
|---|---|---|---:|---:|---:|
| 1 | `ops-translate` — translate 10k × 150 bp (1.5 Mbp) | 0.136 s (11.0 MB/s) | naive dict-per-codon 0.036 s | **3.8×** | 0.100 s |
| — | `ops-gc` (Phase 1, delivered) | 0.697 s | `bytes.count` 0.107 s | 6.5× | 0.59 s |
| — | `ops-revcomp` (Phase 1, delivered) | 0.280 s | `bytes.translate` 0.041 s | 6.8× | 0.24 s |

The alternatives column is what makes this row worth doing.  A vectorised
`numpy` translator — bases to two bits, three bases to one of 64 entries — is
*slower* than the plain dict loop here (0.058 s against 0.036 s), because at
150 bp per read the per-call numpy overhead is never amortised.  So the
alternative to beat is a Python loop, and it is already 3.8× faster than the
reference: 100 Mbp of CDS costs 9.1 s in Biopython and 2.4 s in a five-line
Python function.  That is the same shape as Phase 1's rank 1 — the reference is
paying for its own object model, and the time is there to be taken.

At 100 Mbp the reference spends 9.1 s; a kernel that does not build a `Seq` per
call and does not look codons up in a dict has no reason to stop at 2.4 s.
`Bio.Seq` is also the single largest `perf` block in the triage (86 rows), and
this workload is the part of it with a measured price.

### Delivered — `biofasting.translate` (2026-10-03)

Same workload at the size the runner now uses, 100,000 reads / 15.0 Mbp, which
is what `ops-translate` prints.  Median of 5, `--quick` corpus, same build
identity.

| Implementation | Time | Throughput | vs Biopython |
|---|---:|---:|---:|
| `Seq.translate`, once per record | 1.340 s | 11.2 MB/s | 1.00× |
| `naive-dict-per-codon` | 0.340 s | 44.1 MB/s | 3.9× |
| `numpy-64-entry` | 0.585 s | 25.7 MB/s | 2.3× |
| **`biofasting.translate`** | **0.088 s** | **171.2 MB/s** | **15.2×** |

The target was 3.8× — beating the *best alternative*, not the reference — and
the delivered row is 3.9× against that alternative and 15.2× against the
reference.  15.0 Mbp in 88 ms is 171 Mbase/s, or about 666 ns per 150 bp read.
A second independent run of the same three rows reproduced it at 1.334 s and
0.086 s, 15.5× against the reference and 4.2× against the loop, which is the
spread to read these numbers with.

That 666 ns is worth writing down, because it is now almost all Python and
almost none of it is the kernel.  Per call: `str.upper().encode("latin-1")`
127 ns, resolving `table="Standard"` to a table 62 ns, fetching the kernel's two
tables off it 83 ns, `_core.translate` 205 ns, and about 190 ns of argument
handling and branches.  The kernel is therefore roughly a third of a read-sized
call — which is the honest reason the row is 15.2× rather than the 49× a single
10 Mbp sequence gives (8.51 ms against 419 ms), where the per-call work is
amortised over ten thousand codons instead of fifty.

Two of those costs were taken during this pass and are worth naming.  Table
resolution was 878 ns of the 1,720 ns call, because `int("Standard")` costs
432 ns *by raising* and "Standard" is what every caller who does not care
passes; the two hashable spellings are memoised now, which is 62 ns.  And the
frame handed to the kernel was `buffer[start:end]`, a copy of the entire input
for every call that is not `cds` — invisible at 150 bp, a 100 MB memcpy on a
chromosome, and now skipped when the frame is the whole buffer.

## Phase 2 ranking, second pass: FASTQ random access (2026-10-03)

The `Bio.Seq` block was measured first and is **not** a kernel target.  Its 86
`perf` rows are thin wrappers over `bytes` that already run at C speed —
`Seq(150).count_overlap` 0.45 µs, `Seq(1 Mbp).complement()` 539.72 µs,
`MutableSeq.append ×200` 10.85 µs — so there is no interpreter time in them to
recover and no row was written for one.  That is a negative result and is
recorded as one.

What the same sweep did surface is rank 2 of the Phase 0 table, still open:
`SeqIO.index` on **FASTQ**.  Phase 0 measured it on FASTA (159.6 ms per 150 bp
slice, ~50,000× behind `pyfaidx`) and Phase 1 delivered a FASTA index; the
FASTQ half of `SeqIO.index` was never touched.  Two new harness workloads
measure it (`fastq-index`, `fastq-random` in `bench/run.py`), both on
`fastq/reads_1m.fastq`, 369.3 MB, 1,000,000 records:

| # | Workload | Biopython 1.88 | Floor | Gap to floor | Stolen |
|---|---|---|---:|---:|---:|
| 1 | `fastq-index` — index 1M records by name | 1.207 s (305.9 MB/s, 828,499 rec/s) | pure-Python build 0.879 s | 1.4× | 0.328 s |
| 2 | `fastq-random` — fetch a record by name | **7.174 µs/access** (139,388/s) | dict lookup + 180-byte slice 0.213 µs | **33.7×** | 6.96 µs/access |

Read the two floors differently, because they say different things.

The **random-access floor is the target.**  0.213 µs is a Python dict lookup
and a slice out of a mapping — no parsing, just the two things a real index
does.  `SeqIO.index` is 34× off that because a fetch re-parses the record from
its offset and builds a `SeqRecord`: the same design gap as the FASTA row, in
the one format where the record is four short lines and the parse is therefore
*most* of the cost.  The kernel to beat is a memcpy; nothing in the way is
fundamental.

The **build floor is not the target** — it is a warning.  0.879 s of pure
Python to walk the file and fill a dict is only 1.4× better than the reference,
which means most of the reference's 1.207 s is not the Python it looks like:
`FastqRandomAccess.__iter__` does one `readline` per *line*, four per record,
and 4,000,000 `readline` calls is where the second goes.  A Python-level
"better index" would spend its 0.879 s and buy almost nothing.  The floor that
matters is our own scanner, which already reads this file in **0.217 s**
(1704.4 MB/s, `fastq-plain`) against `biopython-seqio`'s 3.194 s — so the build
has to be done by the kernel that already exists, and the target is the 5.6×
between 1.207 s and that scan, not the 1.4× between 1.207 s and a Python loop.

Combined, the two rows are worth doing for the same reason the FASTA index was:
the second one is a design gap with a five-figure ratio behind it, and the first
one is only cheap if it is the *same* kernel rather than a second one.

### What the reference's index actually promises, measured

Before writing the build, the reference's index was probed for its contract
rather than read for it.  Keys are `line[1:].rstrip().split(None, 1)[0]` — the
rstripped header's first word; duplicates raise `ValueError: Duplicate key
'x'`; a missing key raises `KeyError`; `keys()` is in file order.  Against every
FASTQ file in the corpus, our scanner's names, order and record bytes agree
with `SeqIO.index` on every file the two both accept.

They do not accept the same set, and the file where they differ is the finding:

```
bench/data/edge/blank_lines.fastq   SeqIO.index: ValueError "Problem with line b'\n'"
                                    SeqIO.parse: reads both records
```

`FastqRandomAccess.__iter__` is a *third* implementation of the FASTQ grammar,
separate from `FastqGeneralIterator` that `SeqIO.parse` uses, and it is
stricter: it rejects a blank line between records that the parser folds away.
The two reference paths disagree with each other on a file the corpus keeps
precisely because it is legal.  The malformed files diverge the other way
round, in wording only — `Problem with quality section` from the index against
`Lengths of sequence and quality values differs for read1 (8 and 4)` from the
parser.

The index built here follows **the parser**, because the parser is this
package's FASTQ contract (see the header of `src/core/fastq.hpp`) and a package
whose index and parser disagree about which files are readable has two grammars
in it, which is the bug the reference is demonstrating.  So
`open_fastq_index` accepts `blank_lines.fastq` and reports a malformed file with
the parser's message plus its offset and record number.  Both halves of that
divergence are asserted in `tests/test_fastq_index.py` rather than left to be
discovered.

### Delivered — `biofasting.open_fastq_index` (2026-10-03)

Same two workloads, same corpus, same harness, median of 5 after one untimed
warmup, both rows digest-gated against Biopython before their times are quoted.
Two independent runs, to show the spread:

| Implementation | Run 1 | Run 2 | vs Biopython |
|---|---:|---:|---:|
| `SeqIO.index`, build | 1.226 s | 1.235 s | 1.00× |
| **`open_fastq_index`, build** | **0.407 s** | **0.437 s** | **2.9×** |
| `SeqIO.index`, fetch | 7.498 µs | 7.451 µs | 1.00× |
| **`open_fastq_index`, fetch** | **0.322 µs** | **0.316 µs** | **23.4×** |

The fetch is the row that was worth doing, and it lands where the floor said it
would: 0.32 µs against the 0.213 µs of a bare dict lookup and a slice, 23×
against the reference.  The build is 2.9× against a scan that costs 0.217 s, and
the gap between 0.217 s and 0.41 s is one hash table — so it is worth recording
what that table costs and what it would cost anywhere else.

Measured, not guessed, by two standalone programs (`std::pmr::unordered_map`,
1,000,000 inserts of 60-byte string keys, same compiler and flags):

| | Cost per insert | For 1M records |
|---|---:|---:|
| keys packed in a vector | 87 ns | 0.087 s |
| keys strided through a 363 MB buffer, as a FASTQ file has them | 137 ns | 0.137 s |
| this index, per-node allocator | 260 ns | 0.260 s |
| this index, arena allocator (**shipped**) | 200 ns | 0.200 s |

Two things follow.  First, the 50 ns between the first two rows is the file's
own memory layout and no container can remove it: a hash table whose probes and
whose key bytes are both random over 369 MB is a cache-miss machine, and 0.137 s
is the floor for any name→record index over this file.  So the shipped build is
0.407 s against a perfect-container 0.337 s, and a hand-rolled open-addressing
table would be buying about 17% — which is why there is not one, and why the
arena allocator (four lines, 60 ms, `std::pmr::monotonic_buffer_resource`) is
where the tuning stopped.

Second, **memory**, which the first version of this table did not report and
should have.  Peak RSS for the 1M-record index is **472 MiB** against
`SeqIO.index`'s **162 MiB**, which looks like a loss until the file is taken out
of it: mapping and touching the 369 MB file is 371 MiB on its own, so this
index's own structures are **~101 MiB against the reference's 162 MiB**.  We hold
*less* bookkeeping and *also* hold the file mapped.  The mapping is page cache
charged to us and reclaimable under pressure, but `ru_maxrss` does not
distinguish, so a caller with a memory budget needs to know it is there.


## Phase 2 ranking, third pass: random access into a *compressed* FASTQ (2026-10-03)

Phase 0's rank 2 was random access, and Phase 1 took the FASTA half while M14
took the plain-FASTQ half.  What neither touched is the form the file is
actually in: production FASTQ is gzipped, and the corpus's own largest artifact
is `fastq/reads_1m.fastq.gz`.  M14's `open_fastq_index` says so itself — it
raises rather than indexing a compressed file — so this is the one hole the
package documents about itself.

Two new workloads, `fastq-gz-index` and `fastq-gz-random`, same access plan and
same digest shape as the plain pair so that all four numbers are comparable.
They are not the same *file* on both sides, and that is the finding rather than
a shortcut:

```
SeqIO.index("reads_1m.fastq.gz", "fastq")
    ValueError: Gzipped files are not suitable for indexing, please use
                BGZF (blocked gzip format) instead.
```

Biopython will not index a plain gzip FASTQ at all.  So the reference row is
handed a BGZF copy of the same 369.3 MB of reads — written once, from the plain
corpus file, into `build/bench-cache/` (never tracked, never corpus), and the
harness's digest gate is what proves the two rows read identical content.  Our
row is handed the `.gz` the corpus ships.

Measured before the kernel exists, with our rows standing SKIPPED so that the
target is a measurement and not a justification:

| # | Workload | Biopython 1.88 | Floor | Gap to floor |
|---|---|---|---:|---:|
| 1 | `fastq-gz-index` — index 1M records by name | 3.236 s (114.1 MB/s, 309,002 rec/s) | inflate 0.658 s + plain index 0.407 s ≈ 1.07 s | 3.0× |
| 2 | `fastq-gz-random` — fetch a record by name | **274.257 µs/access** (3,646/s) | dict lookup + 180-byte slice, 0.213 µs | **1,288×** |

Read the two rows differently, and this time the *fetch* is the story by an
order of magnitude more than it was for the plain file.

The build is 3.2 s against a plain-file build of 1.226 s.  BgzfReader is a
seekable binary handle with a block cache, and `FastqRandomAccess.__iter__` asks
it for one *line* at a time — four million `readline` calls, now through a
decompressor rather than through the page cache.  That is 2.6× slower than
parsing the uncompressed file, for a file 2.1× smaller on disk.  The floor under
it is our own: inflate the whole 172 MB into 369 MB (0.658 s, libdeflate) and
build the index M14 already builds (0.407 s), which is 3.0× — the same modest
ratio M14's build showed, for the same reason (a build is one sequential pass
and there is not much of it to remove).

The fetch is where the file *format* stops being a detail.  **274 µs per
access, against 6.8 µs on the uncompressed file**: every fetch seeks to a byte
offset inside a BGZF stream, which means finding the enclosing block, inflating
it, and throwing away the cache for the next seek.  Random access and
compression are in direct tension, and Biopython resolves it by paying the
decompressor on every access rather than indexing the decompressed bytes.  1,288×
off a dict lookup and a slice is not a tuning gap; it is the same design gap as
Phase 0's rank 2, in the case where the cost of re-reading the record is a
decompressor instead of a parser.

Also worth recording, because it is part of the honest price: **producing the
BGZF file takes 18.7 s** for these 369.3 MB with Biopython's own `BgzfWriter`
(compression level 6).  A pipeline that wants random access to a compressed
FASTQ with the reference today pays that once per file, on top of the index.
That is not an argument against the reference — `bgzip` in C is far faster — but
it is the reason the row's ratio understates the end-to-end difference rather
than overstating it.

The target is therefore two numbers, and the second one is the point: a build
within 3× of a floor that is essentially "inflate it and scan it", and a fetch
that should land where M14's landed, at a memcpy.

### Delivered — `open_fastq_index` on gzip (2026-10-03)

Same two workloads, same access plan, same digest gate.  Two independent runs,
median of 5:

| Implementation | Run 1 | Run 2 | vs Biopython |
|---|---:|---:|---:|
| `SeqIO.index`, build (on BGZF) | 3.456 s | 3.220 s | 1.00× |
| **`open_fastq_index`, build (on the plain `.gz`)** | **1.058 s** | **1.060 s** | **3.2×** |
| `SeqIO.index`, fetch | 272.108 µs | 273.881 µs | 1.00× |
| **`open_fastq_index`, fetch** | **0.332 µs** | **0.362 µs** | **~800×** |

The build lands where the floor said it would — 1.06 s against an inflate of
0.658 s plus an index of 0.407 s, which is to say the two costs simply add and
neither has anything left in it.  `libdeflate` is doing the whole 172 MB → 369 MB
in 0.658 s and the scan-and-index is the same 0.4 s M14 already measured on the
uncompressed file; there is no third cost hiding in the combination.

The fetch is the row the milestone exists for.  **0.33 µs from a gzip file**,
which is the same number M14 gets from an uncompressed one (0.322 / 0.316 µs)
and the same order as the 0.213 µs floor.  That is the whole argument for
inflating up front: the cost of a fetch does not depend on how the file is
stored, because by the time an index exists the file is not stored that way any
more.  The reference's 272 µs is 820× that, and the cause is structural rather
than a matter of tuning — a seek into a BGZF stream invalidates the block cache,
so every access pays a decompressor for a block it will not reuse.

Two things this table does not say, and should:

- **The two rows read different files.**  Ours reads `reads_1m.fastq.gz`;
  the reference reads a BGZF copy of the same bytes because it refuses the
  plain `.gz` outright.  The digest gate is what makes the comparison legal —
  both rows produce identical `(name, sequence, quality length)` for 1,000,000
  records and for the 2,000 fetched — but the reader is entitled to know that
  the reference's side of it is a file it had to be given.
- **The BGZF copy costs 18.7 s to produce** with Biopython's own writer, once
  per file, before the reference's index build can start.  So the end-to-end
  ratio for a user starting from the file they have is not 3.2× but two orders
  of magnitude, and the table's 3.2× is the conservative number.

And one thing it does say that is easy to miss: `biofasting-index-random`
reports **451.4 MB/s**.  That is not a throughput — it is 2,000 memcpys of 150
bytes with a dict lookup each, which is why the row is reported per access.  A
reader comparing MB/s columns across these tables would conclude the compressed
fetch is *faster* than the uncompressed one, and it is the same speed; the
column is meaningless at this size and the µs column is not.

## Package result: the wheel carries its own inflate (2026-10-03)

Not a new target — every kernel in this file was already built and measured.
This is the one number in the project that was **about the wrong program**, and
fixing it is the reason the gz rows were re-measured rather than carried over.

Every `.gz` row above was produced by an extension linked against Debian's
`libdeflate.so.0`, which `build_info()` reported as `1.23 (shared)`.  A shared
link means the wheel imports only on a machine that also has libdeflate
installed, and CI could not see it because CI installed libdeflate too — the
step that looked like it was preventing the problem was hiding it.  The host's
static archive is not an escape either, and that was measured before it was
assumed:

```
ld: libdeflate.a(deflate_decompress.c.o): relocation R_X86_64_PC32 against
    symbol `libdeflate_x86_cpu_features' can not be used when making a shared
    object; recompile with -fPIC
```

So release 1.23 is vendored in `third_party/libdeflate` and compiled by this
build into a static PIC archive; `build_info()` now says
`1.23 (vendored, static)`, and `readelf -d` on the installed extension lists
`libstdc++`, `libgcc_s`, `libc` and `libm` and no libdeflate at all.

What it costs and what it changes, both measured:

| | system `1.23 (shared)` | vendored `1.23 (vendored, static)` |
|---|---:|---:|
| extension size | 344,960 B | 377,792 B (**+32 KB, +9.5%**) |
| `fastq-gzip` — `SeqIO.parse` | 5.242 s | 5.240 s |
| `fastq-gzip` — `biofasting` | 0.864 s | 0.863 s (**6.07×**) |
| `fastq-gz-index` — `SeqIO.index` build | 3.456 / 3.220 s | 3.226 s |
| `fastq-gz-index` — `open_fastq_index` build | 1.058 / 1.060 s | 1.049 s (**3.1×**) |
| `fastq-gz-random` — `SeqIO.index` fetch | 272.108 / 273.881 µs | 271.915 µs |
| `fastq-gz-random` — `open_fastq_index` fetch | 0.332 / 0.362 µs | 0.306 µs (**~890×**) |

The rows do not move, which is the result worth having: the vendored copy is the
same release compiled with upstream's own `-O2 -DNDEBUG` rather than this
project's `-O3`, and the difference between the two is inside run-to-run
variance.  It also means every gz number reported before this date stays
comparable, and `-O3` was deliberately *not* tried — an unmeasured speedup in the
dependency would have made every row above it incomparable instead.

The escape hatch stays a supported configuration rather than a debug switch, and
it was exercised rather than asserted: with
`-DBIOFASTING_VENDORED_LIBDEFLATE=OFF` the suite runs **503 passed, 2 skipped**,
and the two skips say which configuration they are about.  What is skipped is
asking a deliberately host-linked build to prove the wheel is self-contained;
what is not skipped is the check that a build *claiming* `vendored, static`
really has no libdeflate in its dynamic section — which is the check that fails
on the old binary, where `readelf` listed `libdeflate.so.0`.

## Phase 2 ranking, fourth pass: `Bio.SeqUtils` (2026-10-03)

`Bio.SeqUtils` is the largest module in the package that is entirely absent
here and entirely in scope: it is the half of the sequence core that *measures*
a sequence rather than changing it.  Ranking it meant answering two questions
per function -- what would a competent pure-Python implementation cost, and is
there anything to port at all -- and the second one has a real answer for one
of them.

There are two tables below and they answer two different questions.  The first
is the ranking pass, taken on a scratch harness before anything was built: it
asks what a *competent pure-Python rewrite* of each function would cost, which
is the question that decides whether a kernel is needed at all.  `crc64` at
1.0× and `gcg` at 1.2× are the justification for writing those two in C++ —
there is no Python rewrite to be had, so the only way to move them is to leave
the interpreter.  The second table is the delivered measurement from this
project's own harness, against `build_info()`, digest-gated.

One number from the ranking pass did not survive.  The `GC123` whole-chromosome
row was printed as **90.674 s**; measured again by the harness it is **8.10 s**
(median of five, not truncated, reproducible twice by hand at 7.93 s).  The 90 s
figure could not be reproduced and is withdrawn rather than explained — it was
almost certainly taken under load during a long run, which is exactly what
`build_info()` and a median-of-five are for.  Nothing else in the ranking pass
was affected: its other rows reproduce to within a few percent against the
harness's own baseline columns.

**Ranking pass.** Corpus: 1,000 reads of 150 bp from
`bench/data/fastq/reads_10k.fastq`, and `chr1` at 40,000,000 bp from
`bench/data/fasta/genome.fasta`.  Median of five, one warm-up, each alternative
asserted equal to the reference on 200 reads and on chr1 before any row was
timed.  That assertion is not ceremony: two of the five alternatives were wrong
on chr1 and right on the reads, and both times the row had already printed a
plausible number.

| row, per 1,000 reads of 150 bp | reference | best pure Python | ratio |
|---|---:|---:|---:|
| `GC123` | 330.67 ms | 2.88 ms (`s[frame::3]` + `str.count`) | **114.7×** |
| `molecular_weight` | 4.60 ms | 3.38 ms (`set` + `count`) | 1.4× |
| `GC_skew` | 4.03 ms | 1.33 ms (windowed `str.count`) | 3.0× |
| `crc64` | 24.30 ms | 23.24 ms (the same loop) | 1.0× |
| `gcg` | 8.16 ms | 6.54 ms (block sum) | 1.2× |
| `CodonAdaptationIndex.calculate` | 59.08 ms | — | — |
| `six_frame_translations` | 122.91 ms | — | — |
| `seq1` / `seq3` | 7.45 / 6.39 ms | — | — |
| `nt_search` | 1.47 ms | — | — |
| `seguid` / `crc32` | 1.31 / 0.30 ms | — | — |

| row, one call | reference | note |
|---|---:|---|
| `GC123`, chr1 | 7.93 s (measured directly; the earlier 90.674 s was withdrawn) | 13.3 M codons, 12 character comparisons each |
| `molecular_weight`, chr1 | — | raises `ValueError: ''N'' is not a valid unambiguous letter for DNA` |

**Delivered measurement.** Corpus: the first 20,000 reads of
`bench/data/fastq/reads_1m.fastq` (3.06 Mbp of 150 bp reads), median of five,
`--max-pass-seconds 60`, every row's full output compared against the
reference's by SHA-256 before it is quoted — `bench/results/sequtils.json` and
`bench/results/sequtils-genome.json`, produced by `bench/run.py`.  Build:
GNU 14.2.0, `-O3`, C++20, x86_64, dispatched to `avx2` (the CPU tops out there),
nanobind 3.1.0, libdeflate 1.23 vendored, Python 3.13.5, 13th Gen Intel
i5-13420H.  `naive-*` rows are the same pure-Python baselines as above, kept so
that the two tables can be read against each other.

| row, 20,000 reads unless stated | Biopython 1.88 | this project | ratio | pure-Python baseline |
|---|---:|---:|---:|---:|
| `GC123` | 0.5982 s | 0.0180 s | **33.2×** | 0.0577 s (10.4×) |
| `GC_skew(window=30)` | 0.0366 s | 0.0087 s | 4.2× | 0.0368 s (1.0×) |
| `molecular_weight` | 0.0968 s | 0.0217 s | 4.5× | — (1.4×, from the ranking pass) |
| `crc64` | 0.4510 s | 0.0165 s | 27.3× | — (1.0×, from the ranking pass) |
| `gcg` | 0.1676 s | 0.0070 s | 24.0× | — (1.2×, from the ranking pass) |
| `translate` → `seq3` → `seq1` | 0.4879 s | 0.2248 s | 2.2× | — |
| `CodonAdaptationIndex.calculate` | 0.1292 s | 0.1299 s | 1.0× | — |
| `six_frame_translations` | 1.7581 s | 0.3854 s | 4.6× | — |
| `GC123`, one 40 Mbp contig | **8.1026 s** | 0.0179 s | **452×** | — |

Three things in those tables are worth more than the ratios.

**`GC123` is the target, and the chromosome row is the reason.**  Its 330 µs per
read is not 150 bases of work, it is 50 codons × 3 positions × 4 letters of
`==`, and it grows linearly with the sequence: 8.10 **seconds** for one
chromosome, against which the best rewrite available in Python is a three-way
stride and eight `str.count` calls, and a kernel answers it in 18 ms.  The
per-read ratio is 33× and the per-chromosome ratio is 452×, because at that size
the reference's cost is the interpreter and ours is memory bandwidth — 167 MB/s
in, against 0.4 MB/s.  A kernel has no pure-Python alternative to beat, and it
is the only row in the module where that is true at this scale.

**`molecular_weight` cannot be ported with a running total.**  CPython 3.12 and
later sum floats with Neumaier compensation, and the reference is written as
`sum(weight_table[x] for x in seq)`.  Measured on this interpreter: over 20,000
random draws from the DNA weight alphabet, a Neumaier reconstruction of that
expression matched it 20,000 times and a naive left-to-right loop mismatched
**5,337 times (26.7%)**.  The kernel sums with compensation or it is not a port.

**`xGC_skew` has no row because it cannot have one.**  It imports `tkinter`,
creates a `Canvas`, draws the two skew curves, calls `canvas.update()`, and
returns `None`; with `DISPLAY` set it takes **237 ms per call**.  There is no
result to compare and no algorithm to port.  It is the one entry in
`inventory/biopython_triage.csv` this pass found to be mis-triaged, and the row
has been corrected rather than left to be re-derived later.

## Phase 2 ranking, fifth pass: the satellites of `Bio.SeqUtils` (2026-10-03)

The fourth pass took the half of `Bio.SeqUtils` that measures a *sequence*.  This
pass takes what is left after it — `ProtParam`, `MeltingTemp`, `IsoelectricPoint`
and `lcc`, 24 of the module's remaining `perf` rows — and asks the same two
questions: what does the reference cost per call, and what would the best
pure-Python rewrite cost.

The corpus question is new here, and it has to be answered in the open: `bench/data/`
is DNA.  There is no protein in it and no primer, so the inputs below are
**generated, not sampled** — SHA-256 in counter mode over the standard SwissProt
amino-acid frequencies, so the same command yields the same residues on any
machine — at 150, 300, 1024 and 10,000 residues.  Oligos are 25, 60 and 1,000
bases.  Python 3.13.5, Linux 6.12.107+deb13-amd64, 13th Gen Intel i5-13420H, median
of seven calls after one warm-up.

**What the reference costs per call**, in microseconds.  Every row is the *first*
call on a fresh object, because that is what a caller pays: `count_amino_acids`
memoises its result into the instance, so timing it on an object that has already
been used measures a cache hit — **0.19 µs**, against the 20 `str.count` calls
behind it.  A ranking pass that missed that would have ranked the cheapest method
in the module first.

| method, per call | 150 res | 300 res | 1024 res | 10,000 res |
|---|---:|---:|---:|---:|
| `flexibility` | 47.84 | 106.41 | **329.44** | **3,501.60** |
| `instability_index` | 12.13 | 25.31 | 85.79 | 873.87 |
| `molecular_weight` | 4.84 | 7.78 | 24.07 | 220.63 |
| `gravy` | 3.95 | 7.24 | 22.13 | 213.14 |
| `isoelectric_point` | 21.31 | 21.79 | 27.33 | 98.91 |
| `secondary_structure_fraction` | 5.45 | 6.40 | 10.66 | 85.14 |
| `count_amino_acids` | 2.91 | 3.60 | 8.46 | 82.29 |
| `aromaticity` | 4.47 | 5.32 | 9.83 | 80.15 |
| `amino_acids_percent` | 3.82 | 4.85 | 9.23 | 75.44 |
| `molar_extinction_coefficient` | 2.85 | 3.78 | 7.95 | 74.25 |
| `charge_at_pH(7.0)` | 2.39 | 2.02 | 2.16 | 4.74 |
| `ProteinAnalysis(seq)` | 0.57 | 0.39 | 0.70 | 3.03 |
| `IsoelectricPoint(seq).pi()` | 21.43 | 21.70 | 26.14 | 96.09 |

| call, per call | 25 bases | 60 bases | 1,000 bases |
|---|---:|---:|---:|
| `Tm_NN(seq, Na=50, Mg=1.5)` | 15.47 | 23.40 | 293.65 |
| `Tm_GC(seq)` | 5.50 | 6.16 | 30.06 |
| `Tm_Wallace(seq)` | 2.45 | 2.40 | 7.59 |
| `salt_correction(Na=50, seq)` | 0.53 | 0.56 | 0.56 |
| `lcc_mult(seq, wsize=20)` | 2.56 | 5.39 | 84.50 |
| `lcc_simp(seq)` | 0.49 | 0.54 | 1.64 |

**What the best pure-Python rewrite costs.**  Each rewrite is asserted equal to
the reference on 200 random proteins of 12–400 residues *before* its time is
quoted, because this module has already produced one rewrite that was right
where it was checked and wrong elsewhere.  Times from the same run as the ratios.

| per call | reference | rewrite | ratio | reference | rewrite | ratio |
|---|---:|---:|---:|---:|---:|---:|
| | **300 residues** | | | **1024 residues** | | |
| `flexibility` | 106.80 µs | 51.72 µs | 2.07× | 376.59 µs | 186.70 µs | 2.02× |
| `instability_index` | 33.03 µs | 40.58 µs | **0.81×** | 89.49 µs | 76.13 µs | 1.18× |
| `gravy` | 7.00 µs | 3.38 µs | 2.07× | 21.72 µs | 10.64 µs | 2.04× |
| | **10,000 residues** | | | | | |
| `flexibility` | 3,605.02 µs | 1,820.49 µs | 1.98× | | | |
| `instability_index` | 899.76 µs | 475.62 µs | 1.89× | | | |
| `gravy` | 229.93 µs | 100.63 µs | 2.28× | | | |

Four things in those three tables are worth more than the numbers.

**`flexibility` is the target, and nothing in Python gets near it.**  It is the
heaviest call in the module at every size — 329 µs on a 1,024-residue protein,
3.5 ms on a 10,000-residue one — and it is a per-residue Python loop that reads
nine window positions through two dict lookups each.  The best rewrite available
in Python is **2.0×**: unrolling the window into nine list indices and one weight
vector.  A kernel is the only thing that changes the order of magnitude, which is
the same argument that justified `crc64` and `gcg` in the fourth pass, and the
opposite of the one that would have justified a kernel for the `GC_skew` baseline.

**The window weights cannot be redistributed, even though the arithmetic says
they can.**  The reference accumulates four *mirror pairs* —
`(flex[front] + flex[back]) * weights[j]` — and then adds index 5 a second time,
so index 4 is never read at all and index 5 carries `0.8125 + 1.0`.  Folding the
four pairs into a symmetric nine-element weight vector gives the same number and a
different float: the rewrite that did it differed from the reference in the last
bit on **200 of 200** proteins, while the rewrite that keeps the grouping differed
on **0 of 200**.  Two multiplications that are equal on paper are not equal in
IEEE-754, and a kernel that reassociates this sum is a porting bug with a
plausible-looking output — which is the whole reason `bench/run.py` compares
digests before it quotes a time.

**The module needs two summations, not one.**  `molecular_weight` and `gravy` are
written as `sum(...)` and CPython 3.12+ sums floats with Neumaier compensation;
`instability_index` is written as `score += value` in a loop and is *not*
compensated.  A port that used one summation for both would be wrong in one of
them, and the existing `molecular_weight_mass` kernel is therefore exactly the
right kernel for `gravy` — it is the same `sum` over a per-letter weight table —
and exactly the wrong one for `instability_index`, which needs the naive
left-to-right accumulation the reference performs.

**Below about 500 residues the Python rewrite of `instability_index` is a
regression.**  It is 0.81× at 300 residues and only 1.89× at 10,000: replacing the
reference's two-character slice with a two-character concatenation costs about
the same, and the flattening of the dict-of-dicts into one dict does not pay for
itself until the sequence is long.  A protein of 300 residues is a normal protein.
This row is worth a kernel for the same reason `flexibility` is — the reference is
a per-residue Python loop and the ceiling in Python is ~1.2× — but it is not worth
quoting a speedup for without measuring the size that was measured.

**Not claimed here.**  `MeltingTemp`'s heaviest realistic call is `Tm_NN` on a
60-base primer at 23.40 µs, and `IsoelectricPoint.pi()` is 26.14 µs on a
1,024-residue protein; both are dominated by a handful of `str.count` calls and
neither has a rewrite that beats it yet.  They stay `perf` rows until a pass
measures a target for them.  The two `lcc` functions are O(n) over a window and
already incremental; `lcc_simp` is 0.54 µs on a 60-mer.

## Delivered result: `Bio.SeqUtils.ProtParam` (M18, 2026-10-03)

The fifth pass ranked this module and named three targets — `flexibility`,
`instability_index` and `protein_scale` — plus one method it did not rank
because it looked cheap: `count_amino_acids`, twenty `str.count` calls, which is
the first statement of six other methods.  This is what was built against that
ranking, measured by `bench/run.py` on the corpus that pass defined.

Reproduce:

```sh
.venv/bin/python bench/run.py --repeat 9 \
    --only protparam-count-amino-acids --only protparam-amino-acids-percent \
    --only protparam-aromaticity --only protparam-secondary-structure-fraction \
    --only protparam-molar-extinction-coefficient --only protparam-flexibility \
    --only protparam-instability-index --only protparam-gravy \
    --only protparam-protein-scale --only protparam-isoelectric-point \
    --only protparam-charge-at-ph
# -> bench/results/protparam.json
```

Corpus: 400 generated proteins — 200 at 150 residues, 120 at 300, 75 at 1,024 and
5 at 10,000, **192,800 residues in all** — SHA-256 in counter mode over the
SwissProt amino-acid frequencies, so the same command yields the same residues
anywhere.  Every row is a **fresh `ProteinAnalysis` object per call**, which is
what a caller pays; `count_amino_acids` memoises into the instance, so a row
timed on a used object would be measuring a cache hit.  Median of nine, one
warm-up, every row's full output compared against Biopython 1.88 by SHA-256
before it is quoted.  Build: GNU 14.2.0, `-O3`, C++20, x86_64, `seqops_level =
avx2`, nanobind 3.1.0, libdeflate 1.23 vendored, Python 3.13.5, 13th Gen Intel
i5-13420H.  All 26 rows validated.

| row, per call on a fresh object | Biopython 1.88 | this project | ratio | pure-Python ceiling |
|---|---:|---:|---:|---:|
| `flexibility` | 193.75 µs | **5.25 µs** | **36.9×** | 176.00 µs (1.10×) |
| `instability_index` | 43.75 µs | **1.00 µs** | **43.8×** | 23.75 µs (1.84×) |
| `protein_scale(kd, 9)` | 235.50 µs | **6.25 µs** | **37.7×** | 205.50 µs (1.15×) |
| `gravy(KyteDoolitle)` | 14.00 µs | 1.50 µs | 9.3× | 8.25 µs (1.70×) |
| `count_amino_acids` | 10.25 µs | 1.75 µs | 5.9× | — (see below) |
| `molar_extinction_coefficient` | 10.25 µs | 1.75 µs | 5.9× | — |
| `amino_acids_percent` | 11.50 µs | 3.00 µs | 3.8× | — |
| `aromaticity` | 11.75 µs | 3.50 µs | 3.4× | — |
| `charge_at_pH(7.0)` | 12.00 µs | 4.00 µs | 3.0× | — |
| `secondary_structure_fraction` | 12.75 µs | 4.50 µs | 2.8× | — |
| `isoelectric_point` | 28.00 µs | 20.50 µs | 1.4× | — |

**The three ranked targets were built, and the fourth was found by measuring
what was built.**  `flexibility` and `protein_scale` came out at 36.9× and
37.7×, against a pure-Python ceiling of 1.10× and 1.15×: the window is a
per-residue interpreter loop and there is no rewrite that makes it not one.
`instability_index` is 43.8× against 1.84×, and the ratio is larger than the
other two for the reason the fifth pass predicted — the reference's per-residue
cost there is a slice plus two dict lookups, and the ceiling in Python is two
short strings and one lookup.

The surprise is the row beneath them.  After the three kernels landed,
`isoelectric_point` and `charge_at_pH` were still at **0.98× and 1.00×** — level
with the reference and no faster — because both begin with
`count_amino_acids()`, which was still the reference's twenty `str.count` calls.
A `str.count` is a C-level scan and a kernel cannot beat it by much *per call*,
but twenty of them read the sequence twenty times and cross into Python twenty
times, and that is the whole of the cost.  One byte-histogram pass fixed it, and
it fixed six rows at once:

| `count_amino_acids`, per call | 150 res | 300 res | 1024 res | 10,000 res |
|---|---:|---:|---:|---:|
| Biopython 1.88 | 4.15 µs | 6.71 µs | 19.41 µs | 174.89 µs |
| this project | 1.19 µs | 1.33 µs | 2.00 µs | 9.47 µs |
| ratio | 3.5× | 5.0× | 9.7× | **18.5×** |

The ratio grows with the sequence because the reference's cost is twenty passes
over the bytes and ours is one; the residual is the call itself.  This is the
one kernel in `protparam.cpp` that **never declines**: `str.count('A')` looks for
one letter and is uninterested in what it does not find, so a residue outside the
twenty is a zero and not an error, and the "a kernel counts, Python decides"
rule has nothing to decide.  Every other kernel there declines rather than
raising, because every other one is reproducing a reference that fails in a way
only Python can spell.

**Two of the three summations are the reference's, exactly.**  `molecular_weight`
and `gravy` are `sum(...)` and are Neumaier-compensated by CPython 3.12+;
`instability_index` is a plain `score += value`.  The `gravy` row reuses the
`molecular_weight_mass` kernel and the `instability_index` row does not, and
`bench/run.py`'s digest gate is what makes that a checked statement rather than
a remembered one: the shortest sequence separating the two summations is
`"EPCM"`, where the naive total is `47.32000000000001` and the compensated one
`47.32`, and `ProteinAnalysis("EPCM").instability_index()` is
`118.30000000000001`, the naive one.

**A `naive-*` row that fails the digest gate is the gate working.**  The
`flexibility` baseline keeps the reference's mirror-pair grouping —
`(front + back) * weight`, index 5 read twice, index 4 never — and passes.  The
obvious rewrite, folding the four pairs into one symmetric nine-element weight
vector, is the same arithmetic on paper and a different double: 200 of 200
proteins differ in the last bit.  It is not in the table because it does not
validate.

**Not claimed.**  `isoelectric_point` at 1.4× is the weakest row here and it is
reported as it stands.  Its cost is now almost entirely `count_amino_acids` plus
a twenty-step binary search whose per-iteration work is nine additions; the count
is amortised differently because a 10,000-residue protein pays 9.47 µs of counting
to resolve a pH to four decimal places.  `MeltingTemp` stays unclaimed: the fifth
pass measured `Tm_NN` at 23.40 µs on a 60-mer and found no rewrite that beats it,
and no kernel has been built against it.
