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
