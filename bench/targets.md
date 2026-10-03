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

## Phase 2 ranking, sixth pass: the remaining `SeqIO` formats (PLAN 2.4, 2026-10-03)

PLAN 2.4 says "remaining `SeqIO` formats, hot-path-first, based on Phase 0
ranking".  The Phase 0 ranking had already answered it, and the answer was
*no*: every format left in `Bio.SeqIO` is filed `parity` in
`inventory/biopython_triage.csv` — "worth having for drop-in compatibility, no
large win expected" — and the `perf` rows there are all closed.  `FastaIO` (13
rows), `QualityIO` (28) and `_index` (21), 62 of them, were delivered in Phase 1.

A prefill is not a decision.  This pass asks the triage's question again with a
timer, over the three flat-file formats that carry essentially all of the
world's sequence data: **GenBank**, **EMBL** and **SwissProt**.

```sh
python3 bench/rank_seqio.py
```

`bench/data/` is DNA and holds none of them.  The corpus is
`bench/seqio_corpus.py`, generated rather than stored, on the same contract as
`bench/gen_data.py`: SHA-256 in counter mode, every decision in integer
arithmetic, no float threshold.  GenBank and EMBL are emitted by
`Bio.SeqIO.write` — those two have a writer, and hand-writing a parser-correct
INSDC record is a week spent learning where column 22 is.  SwissProt has **no
writer** upstream (`Reading format 'swiss' is supported, but not writing`), so
that emitter is written out here and the parser that accepts it is the check on
it.  `tests/test_seqio_corpus.py` holds the corpus to its own promises —
record counts, lengths, ids, alphabets, byte-identical rebuilds — before any
number below rests on it.

### The question that decides

A parse spends its time in three places: moving bytes, building Python objects,
and its own logic.  C can take the first two and not the third, so each format
is measured against both floors:

* **objects** — build the same `SeqRecord` tree from scalars already in hand,
  locations, qualifiers and annotations included.  This is Python work no parser
  avoids and no kernel can go below.
* **scan** — a whole-file `bytes.translate`/`upper`: what moving the bytes costs
  at C speed with nothing above it.

Everything the reference spends above the sum of the two is its own
line-by-line logic, and that is what a kernel removes.

| corpus | records | MB | ref µs/rec | objects | scan | lines | ref / floor |
|---|---:|---:|---:|---:|---:|---:|---:|
| `genbank-1kb` | 3,000 | 5.19 | **48.6** | 5.2 | 1.9 | 5.0 | **6.9×** |
| `genbank-10kb` | 300 | 4.38 | 205.5 | 28.8 | 16.8 | 34.6 | 4.5× |
| `genbank-1kb-bare` | 300 | 0.47 | **28.4** | 1.0 | 1.5 | 2.0 | **11.1×** |
| `embl-1kb` | 3,000 | 5.38 | **46.8** | 4.5 | 1.9 | 5.3 | **7.3×** |
| `embl-10kb` | 300 | 4.63 | 223.0 | 28.5 | 16.8 | 37.4 | 5.0× |
| `swiss-300aa` | 5,000 | 4.68 | 47.5 | 8.3 | 1.1 | 2.8 | 5.1× |

Python 3.13.5, Linux 6.12.107+deb13-amd64, 13th Gen Intel i5-13420H, median of
seven calls after one warm-up; `objects` is a median of five, because building a
5 MB file's worth of `SeqRecord` objects five times is already the slowest thing
in the script.

> **Re-measured later the same day, and this column does not reproduce.**  On an
> idle machine this script returns SwissProt 36.4 rather than 47.5,
> `genbank-10kb` 183.9 rather than 205.5, `genbank-1kb` 51.0 rather than 48.6,
> and an object floor of 3.3 for `genbank-1kb` rather than 5.2.  Every row moved,
> in both directions, the widest by 1.6×.  **Read the absolute µs
> below as approximate and the `ref/floor` column as indicative.**  What
> survives — and what the pass was for — is the shape: the reference is still
> 5.4–10.3× above its own floor.  The delivered numbers at the end of this file
> are timed in one harness in one run and do not rest on this column.

**The triage's verdict does not survive this.**  The reference spends 47–49 µs
per kilobase GenBank or EMBL record and 28 µs on the same record with no
features, where the bytes cost 1.5 and the object tree costs 1.0.  Between
**two-thirds and nine-tenths** of the reference's time is its own logic — a
line-at-a-time `startswith` dispatch, an `ORIGIN` block assembled by a nested
loop, a location string re-parsed per feature.  None of that is object
construction, and none of it is data movement, which is why the headroom column
is 4.5×–11× rather than the 1.0×–1.2× that `parity` predicted.

The feature cost is measured directly rather than inferred.  `genbank-1kb` and
`genbank-1kb-bare` are the same sequences from the same seed, built with feature
density as the only difference: 48.6 µs against 28.4 µs says one feature per
kilobase costs **20.2 µs**, and the `objects` column says building a `SeqFeature`
with its location and three qualifiers costs 4.2 µs of it.  The other 16 µs is
the reference reading the location string and the qualifier block.

The same measurement found two things about the corpus itself, both of which
were my errors and both of which the gates caught.  The feature layout first
drew from the **same** keystream as the sequences, so `genbank-1kb` and
`genbank-1kb-bare` diverged from their second record onward and the one pair the
whole decomposition rests on was not a pair at all; it draws from its own stream
now.  And EMBL right-aligns the residue count on the last line of every `SQ`
block, which a naive extraction appends to the sequence — the corpus's own
comparison column is what caught that one.

### The smallest rewrite there is

The floor under a *Python* rewrite, measured because it is the number that says
whether the win needs C at all: split the file on `//`, take the id from the
first line, strip the sequence out of the `ORIGIN`/`SQ` block.  It reproduces
**ids and sequences exactly** — asserted on every row before its time is
quoted — and produces nothing else: no features, no description, none of the
nine annotations.

| corpus | ref µs/rec | rewrite µs/rec | ratio | ids and sequences identical |
|---|---:|---:|---:|---|
| `genbank-1kb` | 48.6 | 5.5 | **8.9×** | yes |
| `genbank-10kb` | 205.5 | 34.9 | 5.9× | yes |
| `genbank-1kb-bare` | 28.4 | 4.6 | 6.0× | yes |
| `embl-1kb` | 46.8 | 5.3 | **8.8×** | yes |
| `embl-10kb` | 223.0 | 34.0 | 6.6× | yes |
| `swiss-300aa` | 47.5 | 3.6 | **13.2×** | yes |

This is a **floor and not a replacement**, and the distinction is the whole
point of printing it: the rewrite is not asked to produce what the reference
produces, so its ratio is not a speedup anyone could ship.  What it establishes
is that the reference's cost is not in the data — six to thirteen times the
whole job is available before a single line of C is written.

### Verdict

**`GenBankIterator`, `EmblIterator` and `SwissIterator` are `perf`, not
`parity`.**  The triage filed them as compatibility work on the reasoning that a
flat-file parse is "no large win"; measurement says the reference is 4.5×–11×
above the floor that C cannot cross and 6–13× above a rewrite that drops half
the output.  That is a real gap, and it is in the reference's own logic.

The writers (`GenBankWriter`, `EmblWriter`, `ImgtWriter`) stay `parity`:
nothing here measured them, and writing is a different job from parsing.

**Not yet a plan to build.**  What this pass produces is the measured target the
project's first rule requires, three of them, with the binding constraint named:
the object floor is 1.0–8.3 µs per record and it is what a kernel cannot get
below.  `genbank-1kb-bare` — 11.1× over the floor, 6.0× over a rewrite, and no
feature objects in the way — is the row to build first.

## Delivered result: the flat-file kernel — GenBank, EMBL, SwissProt (M20, 2026-10-03)

The sixth pass ranked the three flat-file iterators as `perf` and named the row
to build first.  This is what was built against that ranking: one pass over the
mapped file, in C++, with no record tree on the Python side until a caller asks
for one.

```sh
.venv/bin/python bench/bench_genbank.py --repeat 9
```

The kernel is `FlatFileIndex` in [`src/core/flatfile.hpp`](../src/core/flatfile.hpp)
and [`flatfile.cpp`](../src/core/flatfile.cpp), reached from Python as
`biofasting.open_genbank()` / `read_genbank()`.  It reads a file once, records
where each record's header and residues are, and reproduces **four fields**:
`id`, `name`, `description` and the sequence.  Corpus: the sixth pass's own, all
six rows, **11,900 records and 24.7 MB**, generated rather than stored on the
same contract as `bench/gen_data.py`.  Every row is digest-gated on all four
fields against `Bio.SeqIO.parse` **before** any time is quoted, and the gate
compares all-or-nothing rather than record by record: a mismatching row is a bug
for `tests/test_genbank.py` to localise, not something to be read out of a table.

| corpus | records | MB | ref µs/rec | rewrite | ours | ref/ours | ours MB/s | gate |
|---|---:|---:|---:|---:|---:|---:|---:|:--:|
| `genbank-1kb` | 3,000 | 5.188 | 43.32 | 5.13 | **3.52** | **12.3×** | 491.6 | OK |
| `genbank-10kb` | 300 | 4.381 | 167.07 | 33.86 | **24.97** | 6.7× | 585.0 | OK |
| `genbank-1kb-bare` | 300 | 0.470 | 25.09 | 4.30 | **3.22** | 7.8× | 487.5 | OK |
| `embl-1kb` | 3,000 | 5.375 | 42.49 | 5.39 | **3.92** | **10.8×** | 457.2 | OK |
| `embl-10kb` | 300 | 4.628 | 183.19 | 34.94 | **27.99** | 6.5× | 551.2 | OK |
| `swiss-300aa` | 5,000 | 4.679 | 26.95 | 3.73 | **1.86** | **14.5×** | 502.5 | OK |

Median of nine passes after one warm-up, the reference and both alternatives
timed in the same harness on the same file; the table above is that run.  Build:
GNU 14.2.0, `-O3`, C++20, x86_64, `seqops_level = avx2`, nanobind 3.1.0,
libdeflate 1.23 vendored, Python 3.13.5, 13th Gen Intel i5-13420H.

**How much of those digits is real.**  The script was run three times — twice at
five passes, once at nine — and the per-record cost moves between them by 2–16%
depending on the row, the two shortest rows the widest: `embl-1kb` 3.92/4.29/4.53
and `swiss-300aa` 1.86/1.86/2.11 µs, against `genbank-1kb-bare` at
3.22/3.22/3.29.  So the ratios carry roughly that much uncertainty, the decimal
places above are a courtesy rather than a precision, and a difference of 10%
between two rows here is not a difference.  What is well outside the spread is
the thing being claimed: the reference is 6–14× away, not 1.1×.  (The
nine-pass run is tabulated because it is the most-sampled, not because it is
the fastest — it is not the fastest on every row.)

**All six rows gate OK**, which is the result that matters most here: the
reference's four fields and this reader's are identical on all 11,900 records.
That covers the shape the sixth pass's own comparison column had already caught a
naive extraction getting wrong — EMBL's right-aligned residue count on the last
line of every `SQ` block, which an extraction that strips spaces and not digits
appends to the sequence — and it covers the rules no corpus row can show, because
the reference's writer puts the same string on both the ID and the `AC` lines and
writes every residue in capitals.  Those are asserted by fixtures spliced out of
the writer's own output instead: an EMBL record whose `AC` disagrees with its ID
line, a lowercase `SQ` block, a `VERSION` carrying a suffix.  The kernel is
deliberately **stricter** than
the reference in one place and refuses rather than reproduces a guess: a GenBank
`ORIGIN` line whose column 10 is not a space makes Biopython warn and shift the
line by a byte; this reader raises and names the byte offset.  That refusal, the
`CONTIG` refusal, and the duplicate-id refusal are asserted in
`tests/test_genbank.py` beside the agreements.

**`ours` beats the pure-Python rewrite on every row while producing more.**  The
rewrite extracts ids and sequences; this reader also produces the name and the
description, and is still 1.2–2.0× faster than it (3.22 against 4.30 on the bare
row, 1.86 against 3.73 on SwissProt).  So C was worth writing here — but the
rewrite column is why that sentence is worth saying at all: a pure-Python parse
of the same bytes was already 5–7× the reference, which is where most of the
distance came from.

**What this ratio is not.**  `ours` produces four fields and nothing else — no
FEATURES table, no annotations dict, no `SeqRecord` — so `ref/ours` is a ratio
against a reference doing strictly more work, exactly as the rewrite's is.  The
script prints that next to the table rather than in a footnote.  It is also why
the delivered numbers cross below the sixth pass's floor: `genbank-1kb-bare` was
ranked at 28.4 µs against a floor of 2.5 (objects 1.0 + scan 1.5), and this
reader does the same record in 3.22 — above the scan alone, below the sum.  That
is not a kernel beating a floor.  The floor's `objects` term built the
reference's *full* `SeqRecord` tree, features and annotations included, and this
kernel deliberately builds none of it; the floor stops bounding this design the
moment the output changes shape.  What still bounds it is the 1.5 µs whole-file
scan, and 3.22 against 1.5 — 2.1× — is the honest remaining headroom on that
row.  Throughput barely moves with record size — 396–585 MB/s over the three
runs, on records that differ by 30× — so the reader is close to
bandwidth-bound and the remaining win is per-line dispatch, not data movement.

One number from the sixth pass did not reproduce, and the honest thing is to say
so rather than to reconcile it.  That pass recorded the reference at
26.4–223.0 µs/rec; re-measured here it is 25.1–183.2, with SwissProt the widest
at 26.95 against 47.5.  The obvious suspect is the harness — the ranking pass
parses an in-memory `StringIO` where this one opens the file — and that was
measured directly: it accounts for **2–7%**, not 76%.  So the sixth pass's own
script was re-run on the now-idle machine, and **its own column does not
reproduce either**: SwissProt 47.5 → 36.4, `genbank-10kb` 205.5 → 183.9,
`genbank-1kb` 48.6 → 51.0, and the object floor for `genbank-1kb` 5.2 → 3.3.
Every row moved, in both directions, by up to 1.6×, and the machine that
produced the recorded column is the machine that just disagreed with it.  A
uniform shift would mean a changed corpus; a shift this shape means the column
was taken under conditions — frequency state, or load from whatever else was
running — that this machine no longer reproduces.  The consequence is narrow and
worth stating plainly: **the sixth pass's absolute µs are approximate and its
`ref/floor` column should be read as indicative, not as a number to hold a
kernel to.**  What the pass established is a shape — that the reference's cost on
a flat file is its own logic rather than the data — and that shape survives: the
re-run still puts the reference 5.4–10.3× above its own floor, and every row's
rewrite still reproduces ids and sequences exactly.  The delivery table above
does not depend on the older column: all three implementations are timed in one
harness, on one file, in one run.

## Delivered result: the FEATURES table — GenBank and EMBL (M21, 2026-10-03)

M19's decomposition ended on a named target rather than a verdict: a flat-file
record with one feature per kilobase costs the reference **20.2 µs** more than
the same record with none, and only **4.2 µs** of that is building the
`SeqFeature` object — the other **16 µs is the reference reading the location
string and the qualifier block**.  M20 then wrote down, in its own scope
statement, that the FEATURES table is what it does not produce.  This is the
reader built against that target.

```sh
.venv/bin/python bench/bench_features.py --repeat 9
```

The kernel is [`src/core/location.{hpp,cpp}`](../src/core/location.hpp) — the
location grammar, positions preserved rather than flattened — and
[`src/core/feature.{hpp,cpp}`](../src/core/feature.hpp), reached from Python as
`biofasting.read_features()` and, for callers who want Biopython's types,
`biofasting.to_seqfeature()`.  Corpus: the sixth pass's own, the four rows that
carry features, **6,600 records, 12,000 features and 19.6 MB**; `swiss-300aa` is
not a row, because SwissProt's `FT` block is `Bio.SwissProt._read_ft`, a second
grammar sharing nothing with the INSDC scanner that GenBank and EMBL share.

**The gate is over the features, not over id and sequence.**  Every feature of
every record — the reference's and ours — is reduced to its type, its location
as the canonical nested tuples `tests/test_location.py` defines, its qualifiers
and its location status, and the whole file is hashed: a row that disagrees is
printed `INVALID` and no time from it may be quoted.  The location is reduced to
kind-and-edges rather than compared as objects, because `==` between two
positions in Biopython compares integer values and would not see a `<` that had
been flattened into an exact coordinate.  `ours+interop` is hashed against the
same reference rows, which is what makes it the same work plus the object.

| corpus | records | features | ref µs/rec | ours | ref/ours | ours+interop | gate |
|---|---:|---:|---:|---:|---:|---:|:--:|
| `genbank-1kb` | 3,000 | 3,000 | 48.83 | **9.82** | **4.97×** | 15.06 | OK |
| `genbank-10kb` | 300 | 3,000 | 205.70 | **84.18** | 2.44× | 134.60 | OK |
| `embl-1kb` | 3,000 | 3,000 | 47.80 | **10.00** | **4.78×** | 14.98 | OK |
| `embl-10kb` | 300 | 3,000 | 219.90 | **86.98** | 2.53× | 139.53 | OK |

Median of nine passes after one warm-up, all three implementations timed in the
same harness on the same file, build as above.  `ref` is `Bio.SeqIO.parse` over
the whole record, which also produces the annotations dict and the references,
so `ref/ours` is a ratio against a reference doing strictly more work and is not
the number to quote.

**The number to quote is the marginal one**, because it is the measurement M19
made, re-derived here rather than trusted.  `genbank-1kb` and
`genbank-1kb-bare` are the same 1 kb records from the same seed with feature
density as the only difference (3,000 records against 300 — the corpus sizes the
counts to make the two files about equally long in bytes, so the subtraction is
of two *per-record* times and never of two file times):

| | per feature | against M19's target |
|---|---:|---|
| reference, measured now | 23.71 µs | recorded 20.2 |
| **ours** — the reading, no Python object | **7.62 µs** | 16 µs of "location + qualifiers" |
| ours + `to_seqfeature` | 12.41 µs | — |

Run three times at nine passes: reference 23.62/23.91/23.71, ours
7.97/7.56/7.62, interop 12.88/12.66/12.41.  So the reference's own 20.2 does not
reproduce either — it lands 23.6–23.9, which is inside the ±20% this file has
already recorded for that pass and is the same caveat M19 wrote above, not a new
one.  **The shape reproduces**, and the shape is what M21 was built against: the
reference spends ~16 µs per feature on the location string and the qualifier
block, and this reader spends **7.62 µs** producing both plus the feature key —
**2.1×** on the component the milestone named, **3.1×** against the whole 23.7 µs
the reference actually spends.  The conversion is a separate and much smaller
story: `to_seqfeature` adds **4.79 µs**, against the **4.2 µs** M19 measured for
the reference's own object construction — that is parity, and it is why the win
is entirely in the reading.  A 2–3× on one component is a real number and a
modest one, and it is reported as it stands.

**Scope, stated rather than discovered.**  This reader reproduces the feature
key, the location and the qualifiers, and it does **not** reproduce the
annotations dict — accessions, taxonomy, references, source, comment, date,
keywords — which stays a `SeqIO.parse` job for now.  Three outcomes are distinguished rather than collapsed into one: `ok`;
`parser_error`, which is the reference catching its own `LocationParserError`
and carrying on with `location = None` plus a warning, so the feature keeps its
type and its qualifiers; and a **refusal** — a FEATURES table this reader cannot
reproduce, which raises and names the byte offset instead of returning a
partially parsed table.  Refusing is the deliberate divergence: the reference's
consumer warns and guesses in places where a wrong feature would then travel
silently, and a feature that is wrong is worse than a feature that is absent.

**Warnings are reported by the kernel and spoken by Python.**  The reader hands
back an ordered list of findings — one entry per repaired origin-wrap part
carrying the text the reference quotes, one per dropped `bond` — and
`to_seqfeature` emits the reference's own words in the reference's own order,
which is observable: `join(bond(1),30..5)` on a circle warns about the bond
first and the wrap second.  Reading features and holding them emits nothing at
all, so a caller who never builds a `SeqFeature` never installs a warning
filter; the corpus emits no warnings on either side, which is why the suite's
`filterwarnings = ["error"]` does not trip on it.  All four rows gate OK on both
paths, and the differential tests behind them are `tests/test_location.py` (764
comparisons plus a seeded 4,000-case fuzz), `tests/test_features.py` and
`tests/test_seqfeature.py`, the three of them comparing against
`Bio.SeqFeature`/`Bio.GenBank` field for field.

## Measured target: the annotations dict — GenBank and EMBL (M22, 2026-10-03)

M21's scope statement named what it does not produce: the annotations dict —
accessions, taxonomy, references, source, comment, date, keywords.  This is the
measurement of that, made before the kernel exists, because a kernel without a
measured target is a guess with a build step.

```sh
.venv/bin/python bench/rank_seqio.py      # the last two blocks
```

Three rows were added to `bench/seqio_corpus.py` for this axis, in their own
`ANNOT_CORPORA` tuple so that no number already recorded above was measured on a
different file: `genbank-1kb-annot` and `embl-1kb-annot` (300 records × 1 kb,
annotated, no features) and `embl-1kb-bare`, the EMBL twin that did not exist.
Each annotated row draws its residues from the same seed and the same domain
string as its bare twin, and the blocks are applied **after** the draw, so the
two files differ in nothing but the header — `tests/test_seqio_corpus.py`
asserts that rather than trusting it, in both directions.

| corpus | bare µs/rec | annotated | **pair delta** | +taxonomy/keywords | +comment | +references |
|---|---:|---:|---:|---:|---:|---:|
| GenBank | 25.06 | 57.01 | **31.95** | 1.21 | 16.44 | 13.28 |
| EMBL | 24.24 | 42.13 | **17.89** | 3.75 | 3.60 | 11.53 |

Median of seven passes after one warm-up.  The three block columns are each
measured **one block at a time** against the bare row, not cumulatively, so a
column is a price and not a difference of prices; they sum to 30.93 and 18.87
against pair deltas of 31.95 and 17.89, which is the run spread this file has
already recorded elsewhere.  Run three times: GenBank 31.95/31.39/31.41 (bare
25.06/25.29/24.88), EMBL 17.89/17.50/18.07 (bare 24.24/24.32/24.13).  The
comment column is the stable one — GenBank 16.44/16.40/16.49, EMBL
3.60/3.15/2.79 — and the taxonomy column is the noisy one, because it is one to
three short lines.

### The comment column is not the price of reading a comment

GenBank pays **16.4 µs** to carry a three-line comment and EMBL pays **3.2 µs**
for the same text over the same three lines, which is a five-fold difference
that the text cannot explain.  It is not the reader: both consumers do
`"\n".join(content)` and nothing else.  It is the *probe* in front of it.  The
GenBank path runs

```python
re.search(rf"([^#]+){self.STRUCTURED_COMMENT_START}$", data)
```

on the first line of every `COMMENT` block, to decide whether the record carries
a structured comment.  `[^#]+` is greedy, `-START##` is not there, and the
engine walks the line backwards one character at a time trying to make the
suffix fit.  Measured against line length, with the pattern compiled once so the
column is the matching and not the `re` cache lookup the reference pays on top:

| first comment line | chars | µs | µs / char² |
|---|---:|---:|---:|
| | 20 | 1.18 | 0.00295 |
| | 40 | 4.38 | 0.00273 |
| | 68 | 11.68 | 0.00253 |
| | 120 | 34.21 | 0.00238 |

Quadratic, and paid **per record**: a GenBank record whose comment's first line
is 68 characters — an ordinary sentence — spends ~12 µs failing to find a
structure it does not have, and the record that has no comment at all pays
nothing, which is exactly the 16.4 µs between `genbank-1kb-bare` and
`genbank-1kb-annot`.  EMBL's `CC` path has no such probe, and pays 3.2 µs.

This is a target with an unusual property: it is a cost the reference pays for a
feature almost no record uses, and a reader that parses the comment without
probing for a structure pays for neither.  It also fixes the shape of what a
kernel must reproduce — the comment **text**, with the line breaks the writer
chose, and not the search that decides it is not a table.

### What the reference actually produces, per format

The two formats disagree about the dict, and the kernel is held to both.  On the
corpus rows, field for field:

| | GenBank | EMBL |
|---|---|---|
| `date` | `'01-JAN-2026'` | **absent** |
| `source` | `''` | **absent** |
| `keywords`, thin header | `['']` | **absent** |
| `keywords`, annotated | `['synthetic', 'benchmark']` | same |
| `accessions` | `[record.id]` | `[record.id]` |
| `taxonomy`, thin header | `[]` | `[]` |
| `organism` | `'synthetic construct'` | same |
| `molecule_type`, `topology`, `data_file_division` | `'DNA'`, `'linear'`, `'PLN'` | same |
| `Reference.comment` | `'primary'` | `''` |

So EMBL's thin header carries six keys and GenBank's nine, the missing three are
`date`, `source` and `keywords`, and the same input text that gives GenBank
`Reference(comment='primary')` gives EMBL `Reference(comment='')`.  The comment
itself comes back wrapped at each format's own width, so the kernel reproduces
the writer's `\n` positions rather than re-wrapping the text.

**What is not in scope here**, stated before the kernel rather than after it:
the `date` field is `annotations["date"]` as the header wrote it and not
`record.annotations["date"]` re-derived from `SeqRecord`; `structured_comment`
is a separate key the reference fills only for records that have one (none of
the corpus does, and it is the probe above that costs the money, not the dict);
and `source` for EMBL is absent rather than empty because the reference never
sets it there.

The target the annotations reader is held to is therefore, per record: **31.9 µs**
of GenBank header and **17.9 µs** of EMBL header, of which the comment's 16.4 µs
is a regex that finds nothing in nearly every file.

## Delivered result: the annotations dict — GenBank and EMBL (M22, 2026-10-03)

M21's scope statement named the annotations dict as what it does not produce,
and the measurement above named the price: **31.9 µs** of GenBank header and
**17.9 µs** of EMBL header per record, most of the GenBank figure being a
structured-comment regex that finds nothing.  This is the reader built against
that target.

```sh
.venv/bin/python bench/bench_annotations.py --repeat 9
```

The kernel is [`src/core/annotations.{hpp,cpp}`](../src/core/annotations.hpp),
reached from Python as `biofasting.read_annotations()`, with
`biofasting.to_reference()` for callers who want Biopython's
`Bio.SeqFeature.Reference`.  Corpus: the two bare/annotated pairs the target was
measured on — `genbank-1kb-bare`/`genbank-1kb-annot` and
`embl-1kb-bare`/`embl-1kb-annot`, 1,200 records and 8.2 MB, each pair drawing
its residues from the same seed and its blocks from the same text, so the header
is the only difference.

**The gate is over the dict, keys and order included.**  Every record's
annotations — the reference's and ours — is reduced to its key list *in
insertion order* and to plain values, references through `to_reference` so that
both sides are the reference's own objects, and the whole file is hashed:
`dict.__eq__` is order-blind and would pass a reader that emitted a fixed key
list, which is exactly what this kernel exists not to do.  `ours+interop` is
hashed against the same rows, which is what makes it the same work plus the
object.  The one place the reduction has to be spelled carefully is a reference:
its `location` holds `SimpleLocation` objects on one side and `(start, end)`
pairs on the other, so both are reduced to `(start, end, strand)` with our own
pairs given the reference's `None` strand — a difference in the answer must not
hide behind a difference in how the answer is written.

| corpus | records | keys | ref µs/rec | ours | ref/ours | ours+interop | gate |
|---|---:|---:|---:|---:|---:|---:|:--:|
| `genbank-1kb-bare` | 300 | 9 | 25.25 | **4.48** | 5.64× | 4.72 | OK |
| `genbank-1kb-annot` | 300 | 11 | 58.19 | **9.09** | 6.40× | 11.40 | OK |
| `embl-1kb-bare` | 300 | 6 | 24.72 | **4.13** | 5.99× | 4.32 | OK |
| `embl-1kb-annot` | 300 | 9 | 41.76 | **9.50** | 4.40× | 11.43 | OK |

Median of nine passes after one warm-up, all three implementations timed in the
same harness on the same file.  `ref` is `SeqIO.parse` over the whole record —
sequences, header, references — so `ref/ours` is a ratio against a reference
doing strictly more work and is not the number to quote.  The timed unit is the
parse alone on both sides: the reduction that the gate compares is run
afterwards, because it costs more on the annotated row than on the bare one
(eleven keys and three `Reference` objects against nine and none) and a timing
that included it would put that difference inside the delta.

**The number to quote is the marginal one**, because it is the measurement the
target was made with, re-derived here rather than trusted:

| | GenBank | EMBL |
|---|---:|---:|
| reference, measured now | 32.94 µs | 17.03 µs |
| recorded target | 31.95 | 17.89 |
| **ours** — the reading | **4.61 µs** | **5.37 µs** |
| ours + `to_reference` | 6.69 µs | 7.12 µs |

Run three times at nine passes: GenBank reference 32.94/31.20/31.14, ours
4.61/4.54/5.07, interop 6.69/6.45/6.56; EMBL reference
17.03/17.60/18.76, ours 5.37/5.36/5.52, interop 7.12/6.87/6.96.  The shape
reproduces and the target lands: GenBank 31.1–32.9 against a recorded
31.39–31.95, EMBL 17.0–18.8 against a recorded 17.50–18.07.  So the header costs
the reference what M22 said it costs, and this reader produces it in **4.6 µs**
on GenBank and **5.4 µs** on EMBL — **6.9×** and **3.3×** on the two formats
respectively.  The GenBank figure is the larger win because the larger part of
the reference's column is the quadratic probe the kernel never runs: it reads
the comment text and does not search it for a structure.

`to_reference` adds **2.1 µs** (GenBank) and **1.8 µs** (EMBL) — the cost of
building Biopython's `Reference` objects for three references, which is the
object half of the story and not this milestone's work.

**What the kernel decides and what Python decides.**  The kernel walks the
header's lines, so it alone knows which lines were there and in what order, and
it hands back the keys it created *in the order it created them*; Python names
them and builds the dict.  That division is not decoration: the reference
inserts a key when the line stating it is consumed, so a `COMMENT` above the
first `REFERENCE` puts `comment` before `references` and one below puts it
after, and both are asserted in `tests/test_annotations.py`.  Membership in the
order is also what says a key *exists* — GenBank's thin header has nine keys and
EMBL's six, and the missing three are `date`, `source` and `keywords`, because
EMBL's scanner reads no `DT` line, has no `SOURCE` consumer, and creates
`keywords` only when a `KW` line was there.

**Scope, stated rather than discovered.**  The reader reproduces the header keys
the two INSDC formats share, for the corpus and for the shapes the tests splice
in.  SwissProt's header is `Bio.SwissProt`, a different reader with different
keys, and is not reproduced — a SwissProt record's annotations refuse.  A header
line that sets a key this reader does not produce — `NID`, `PID`, `DBSOURCE`,
`SEGMENT`, a structured comment, a GenBank `LOCUS` line in the pre-229.0 layout
— makes the *record's* annotations refuse and name the shape, rather than return
a dict that is missing a key the reference would have set; the same contract the
features table has, and for the same reason.  `gi` **is** reproduced, because
the files in circulation carry it and a reader that refused them would refuse
most of what `SeqIO.parse` is pointed at.  Refusals are per record and not per
file: the sequence and the other records stay usable, which
`tests/test_annotations.py` asserts alongside the refusal itself.

## Phase 2 ranking, seventh pass: the writers (PLAN 2.4, 2026-10-03)

The sixth pass re-asked the triage's question for the *readers* of the remaining
`SeqIO` formats and the answer was no: the reference is 4.5×–11× above the floor
C cannot cross.  The one direction it left untouched is the other one, and a
prefill is not a decision for `write` any more than it was for `parse`.  So this
pass asks it: `bench/rank_writers.py`, over the same corpus, three floors deep.

* **reference** — `SeqIO.write`, which is what a caller would otherwise use.
* **rewrite** — the identical bytes from the smallest pure-Python writer that
  can produce them, with the per-character loops replaced by a C-level call
  where one exists: `bytes.translate` for the FASTQ quality string, `str.join`
  for the FASTA and QUAL wraps, one `write` for the whole file instead of one
  per line.  Verified byte-for-byte against the reference **before** any time is
  quoted, so it is a floor under the same output and not a different program.
* **payload** — the raw sequence and quality bytes with no formatting above
  them: what a kernel must at least pay.

| row | format | records | reference µs/rec | rewrite | payload | writes/rec |
|---|---|---:|---:|---:|---:|---:|
| `fasta-1kb` | fasta | 3,000 | **2.18** | 2.00 (1.09×) | 0.42 | 18.0 |
| `fastq-150bp` | fastq | 10,000 | **5.26** | 1.09 (**4.81×**) | 0.60 | 1.0 |
| `qual-150bp` | qual | 10,000 | **15.40** | 8.82 (1.75×) | 0.65 | 9.0 |

Two of the three say something a summary would not.  **FASTQ writing is 4.8×
reachable without C at all**: the reference walks the quality list one base at a
time through a dict (`_phred_to_sanger_quality_str[qp]`), and
`bytes(qualities).translate(table)` is the same mapping in one C call.  **FASTA
writing has no Python-reachable headroom** — a byte-identical pure-Python writer
is *slower* than the reference (2.00 against 2.18), so the reference is already
at the floor Python can reach; what is left is its eighteen `write` calls per
record and the slice objects behind them, which only a buffer-building kernel
removes, and the 0.42 µs payload is what such a kernel is bounded by.  **QUAL is
the worst writer in the family in absolute terms**: 15.4 µs for 150 bases, 103 ns
a base, because every quality is formatted with `"%i" % round(q, 0)` and the
lines are packed by popping from the front of a list — O(n) a pop, O(n²) a
record — and the payload is 0.65.

### The flat-file writers, attributed

A byte-identical GenBank record is a page of column rules, so there is no
hand-written twin here; the reference's own methods are timed instead against a
handle that discards.  Times are µs per record.

| | `genbank-1kb` | `genbank-1kb-annot` | `embl-1kb` | `embl-1kb-annot` |
|---|---:|---:|---:|---:|
| `write_record` | **29.00** | 27.83 | **37.70** | 36.94 |
| `_write_the_first_line(s)` | 2.61 | 2.51 | 1.31 | 1.24 |
| `_write_references` | – | 5.92 | – | 4.03 |
| `_write_comment` | – | 1.31 | – | 1.03 |
| `_write_sequence` | **12.61** | 12.89 | **26.88** | 26.75 |
| rest (DEFINITION→SOURCE, keywords, ID) | 13.78 | 5.20 | 9.51 | 3.53 |
| `write` calls/record | **148.7** | 166.0 | **169.6** | 191.0 |

A method with nothing to write on a row is `–`, which is why the bare and the
annotated row are both here: `_write_references` raises `KeyError: 'references'`
and `_write_comment` raises `IndexError` on a record carrying neither, so a
table over the bare rows alone would have said the header is cheap.  The residue
block is the largest single item on both formats and the reason is the write
count: GenBank emits `ORIGIN` **ten bases at a time** (`f" {data[words:words+10]}"`,
170 writes a record) and EMBL six blocks a line, so the 12.6 and 26.9 µs are
almost entirely interpreter crossings rather than formatting.

**Verdict.**  `FastqPhredWriter`, `QualPhredWriter` and the FASTA writer are
`perf`, and measured rather than assumed for the first time.  The flat-file
writers are `perf` too, on a target that is larger than any of them — but a
complete GenBank record means *writing* the annotations dict and the FEATURES
table, which M21 and M22 only read, so that one is recorded here and deferred
rather than started.  What this pass delivers is the sequence-only half:
`write_fasta`, `write_fastq` and `write_qual`, against the floors above.

**Targets the kernels are held to.**  FASTQ ≤ 1.0 µs/record (a 5× over the
reference and at the measured pure-Python floor, with the quality encoding in
C), FASTA ≤ 1.0 µs/record on a 1 kb record (a 2.2× against a reference no
Python rewrite beats), QUAL ≤ 2.0 µs/record (a 7.7× against a reference whose
per-base cost is a Python format call).  The payload column is the bound no
kernel may claim to have gone below.

## Delivered result: the writers (M23, 2026-10-03)

The targets above, answered on the same corpus by the same harness —
`bench/bench_writers.py`, which gates the whole file against `SeqIO.write`
before it quotes anything.  Every row is the writer alone: the records are
converted once, outside the timer, and the reference column is `SeqIO.write`
over the equivalent `SeqRecord` objects.

| row | records | reference µs/rec | ours | ratio | +interop | payload | target | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `fasta-1kb` | 3,000 | 2.12 | **0.46** | 4.6× | 1.09 | 0.12 | 1.0 | PASS |
| `fastq-150bp` | 10,000 | 5.23 | **0.43** | 12.1× | 2.53 | 0.15 | 1.0 | PASS |
| `qual-150bp` | 10,000 | 15.22 | **0.93** | 16.4× | 2.95 | 0.15 | 2.0 | PASS |

The quality encoder on its own — `phred_to_sanger` against
`_get_sanger_quality_str`, both starting from `letter_annotations` so the row is
the encoding and nothing else — is **0.48 µs** a record against **4.51**, a
**9.3×**.  That is the 4.81× the ranking pass found reachable with
`bytes.translate`, so the kernel beats the pure-Python floor rather than
matching it, which is what it was written for.

Three columns need a word each.  `+interop` is the same output reached through
`from_seqrecord`, and it is **above** the writer by 2.4×–6× because that
conversion is per-record Python: a title rule, `str(record.seq).encode` and,
for FASTQ, the quality encoding.  It is reported rather than folded in because
it is the path a caller holding `SeqRecord` objects actually takes, and a table
that showed only the pair path would be describing an API nobody calls.  It is
still below the reference on every row.  `payload` is re-measured here and is
lower than the ranking pass's column (0.42/0.60) because that one converted from
`SeqRecord` objects while this one copies the bytes we already hold — the
number is a floor under *this* harness, which is why it is measured in it rather
than quoted.  And `ours` for QUAL is above `ours` for FASTQ despite writing less
text, because the scores must be rendered as decimal: the kernel does that one
digit at a time into the buffer, which is the cost the reference pays a Python
format call for.

**What the kernels reproduce rather than tidy.**  A wrapped FASTA record with no
sequence writes nothing after the title (`for i in range(0, 0, 60)` never runs)
while the unwrapped branch writes a blank line — two outputs from one record,
and both are produced.  A QUAL width of one to five takes the reference's
`pop(0)` "safe wrapping" branch, which the kernel does not have, so it is
reproduced in Python rather than approximated; six and up is the kernel's
`rfind` cut.  Where that cut would find no space in a window — the place
Biopython 1.88 writes the same string forever — the kernel raises instead, and
the tests assert the *property* that makes it unreachable (a score is at most
three digits and every token but the first is preceded by a space, while the
kernel only sees widths of six and up) rather than pretending to exercise it.
A score above PHRED 93 is truncated to `~` **and warned about**, so a caller who
would see the reference raise under warnings-as-errors sees this raise too.

**Still open.**  The flat-file writers (`GenBankWriter`, `EmblWriter`) are
`parity`, measured and deferred: a byte-identical record means writing the
annotations dict and the FEATURES table, which M21 and M22 only read.  Their
attribution above — 29.0 and 37.7 µs a record, 148.7 and 169.6 `write` calls —
is the target when that work starts, and nothing in it has been delivered.

## Phase 2 ranking, eighth pass: pairwise alignment (PLAN 2.1, 2026-10-04)

The plan files alignment as a *wrapper* — parasail or WFA2-lib behind the
dispatch interface, reusing an optimized Smith-Waterman rather than writing one.
That premise was checked before anything was built, because it can be false
here: `Bio.Align.PairwiseAligner` in 1.88 is **not** Python.  It is a C
extension (`Bio/Align/_pairwisealigner.cpython-313-x86_64-linux-gnu.so`), a
scalar DP at 3.1–3.5 ns a cell, so the question is where the difference is and
whether it survives a wrapper's own call overhead.

`bench/rank_align.py`, six rows, median of seven after warm-up.  Every row is
**gated before any time is quoted**: both sides align every pair of the row
with the same scheme and the scores must be equal or the run stops with the
mismatch printed.  All six gate equal.

| row | pairs | width | reference µs/align | parasail µs/align | speedup |
|---|---:|---:|---:|---:|---:|
| `dna-global-150x512` | 512 | int16 | 70.39 | 11.74 | 6.00× |
| `dna-global-1k` | 1 | int16 | 3511.10 | 170.93 | 20.54× |
| `dna-global-10k` | 1 | int16 | 344746.09 | 19931.11 | 17.30× |
| `dna-semiglobal-150-vs-10k` | 1 | int16 | 5179.27 | 547.74 | 9.46× |
| `dna-local-1k-vs-20k` | 1 | int16 | 131863.34 | 3813.32 | 34.58× |
| `protein-blosum62-300` | 1 | int16 | 279.47 | 32.64 | 8.56× |

With traceback (the reference's `align()` iterator against parasail's
`*_trace_scan`): 16.23×, 24.35×, 9.17×, 37.58×, 17.28×, 16.15× in the same row
order.  Per cell the DP itself is **7–17×**: reference 3.128 / 3.518 / 3.447 /
3.452 / 6.600 / 3.116 ns against parasail 0.522 / 0.171 / 0.199 / 0.365 /
0.191 / 0.364 ns.  The local row's 34.6× is partly the reference's own local
mode being slower (6.6 ns a cell, not 3.4).

**The finding that made the table real — which parasail function is called.**
The first run of this pass used the *unsuffixed* `nw_scan`/`sg_scan`/`sw_scan`
bindings and reported a weak 1.4×–8.6×.  Those are the generic dispatch at
**2.2 ns a cell**.  The explicitly sized bindings are the SIMD kernels — `_8` is
a full sixteen-lane int8 register, `_16` eight lanes — at **0.17–0.52 ns a
cell**, ten to twenty-eight times the unsuffixed call.  So the wrapper calls no
unsuffixed binding at all: it picks the width from the scheme's own score bound,
`n * (largest |matrix| + extend) + open`, and the gate is what proves the choice
did not overflow — a saturated width returns a wrong score, not an error.
Every row here lands on int16; `dna-global-10k` at 2/−2/10/1 bounds at 30,010
and fits.  The bound is deliberately loose, because one width too wide costs a
factor of two and one too narrow is a wrong answer.

The wrapper's overhead was measured too, so that it is a known quantity before
the wrapper exists: **~2 µs a call** (a 64-cell call is 2.07 µs against the
reference's 2.23).  That is why the 150×150 row is the weakest — its DP is
11.7 µs, of which about two are the crossing.

**Floors, and what was left out.**  The floor is parasail itself: a pure-Python
DP is minutes on the 10 kb pair, so it is stated and not carried as a column.
**WFA2-lib is not installable from PyPI** — no package, not a missing wheel — so
it is a *vendoring* decision and not a measurement, and it is absent from the
pass rather than estimated.  The row it exists for is `dna-global-10k` at 1%
divergence: an edit-distance algorithm against a full DP.  That is the row to
measure if it is ever vendored.

**Verdict: `align`/`score` behind a wrapper is `perf`, 6×–35× score-only and
9×–38× with traceback**, with one packaging caveat measured here: `parasail`
1.3.4 ships wheel for x86_64 linux, x86_64 macOS, win32 and win_amd64 **and no
aarch64 of any kind**, so on the CI's `macos-14` (Apple Silicon) and in the
aarch64 wheel job it can only come from the sdist.  The wrapper is therefore an
optional extra, and the number above is honestly an x86-64 number.

**Still open.**  The delivered wrapper's own numbers (through
`biofasting.alignment`, not through parasail directly, which is what the table
measures) are recorded in the delivered section below; the arm64 packaging
decision — require parasail and pay a source build, or vendor its C into
`third_party/` as libdeflate was — is the owner's.
