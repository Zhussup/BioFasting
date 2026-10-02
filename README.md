# BioFasting

<p align="center">
  <img src="misc/biofasting.gif" alt="BioFasting" width="620">
</p>

<p align="center">
  <img src="https://img.shields.io/badge/C%2B%2B-20-00599C?logo=cplusplus&logoColor=white" alt="C++20">
  <img src="https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white" alt="Python 3.9+">
  <img src="https://img.shields.io/badge/status-pre--alpha-orange" alt="pre-alpha">
  <img src="https://img.shields.io/badge/tested--against-Biopython%201.88-4c9349" alt="tested against Biopython 1.88">
  <img src="https://img.shields.io/badge/license-TBD-lightgrey" alt="License TBD">
</p>

A C/C++-backed core library for sequence bioinformatics, built to close the
7–10× performance gap Biopython leaves open on hot paths — aiming to become the
*numpy of bioinformatics*: a fast substrate that other tools build on, not
another API clone on top of it.

**Status:** pre-alpha, Phase 0 (evidence) **complete**. No library code yet — the
API inventory, a reproducible benchmark corpus, and a profiled, ranked list of
targets came first. Phase 1 (the FASTQ/FASTA core) starts from that ranking:
see [bench/targets.md](bench/targets.md).

## Why this exists

On a draft benchmark (100k FASTQ reads, 30 MB, Intel i5-13420H, single-threaded,
informal single-machine numbers):

| Workload | Approach | Time | Throughput |
|---|---|---:|---:|
| Parse 100k reads | `Bio.SeqIO.parse` (1.88) | 0.31 s | 319k rec/s |
| Parse 100k reads | naive pure-Python parser | 0.03 s | 3.42M rec/s |
| Parse 100k reads | `pysam.FastxFile` (C/htslib) | 0.05 s | 2.21M rec/s |
| Reverse complement (same 100k) | `SeqRecord` loop | 0.35 s | — |
| Reverse complement (same 100k) | `str.translate` | 0.06 s | — |

![Time per 100k FASTQ reads: Biopython vs. alternatives](misc/why-this-exists.png)

(chart rendered by `bench/plot_whythis.py`; the table stays as the data source)

Pure Python already beats `SeqIO` by ~10× on parsing; a C core with runtime SIMD
dispatch and `libdeflate` for `.gz` targets the larger gap. Real-world FASTQ is
compressed, so decompression — not tokenization — is usually the first
bottleneck to defeat.

### Measured on the real corpus (Phase 0)

That draft is what started the project. Phase 0 replaced it with a
correctness-gated corpus and a ranked list of targets. On the 1M-read corpus,
the two largest gaps are FASTQ parsing — `SeqIO.parse` 3.2 s vs 0.45 s for a
naive pure-Python loop, **7.1×**, with **55.7%** of Biopython's time in its own
Python frames — and FASTA random access, where 150 bp costs **~160 ms** through
`SeqIO.index` against **3 µs** through pyfaidx, because a `SeqIO` index stores
whole-record offsets and so re-parses the record for every slice.

![Gate G1 ranking — Biopython vs. the best available alternative, per workload,
log scale](misc/g1-ranking.png)

(chart rendered by `bench/plot_g1.py`; the table in `bench/targets.md` stays the
data source)

The gap is not a floor to be matched, it is interpreter work to be deleted:
between a third and two thirds of every Biopython row is Python bytecode inside
`Bio/**`, while the rows already in C — pysam, pyfaidx, the zlib decompression
floor — have nothing left to remove.

![Share of profiled time spent in Biopython's own Python frames, per
implementation](misc/profile-share.png)

(chart rendered by `bench/plot_profile.py`; the profiles are in
`bench/profiling.md`)

The full ranking, with the caveat attached to each alternative and the Phase 1
ordering, is in [bench/targets.md](bench/targets.md); the profiles behind it are
in [bench/profiling.md](bench/profiling.md).

## Principles

1. **Benchmark first.** No code until a bench harness ranks workloads by
   *time lost* (Biopython time minus best-available-alternative time).
2. **Substrate, not a clone.** Biopython stays; BioFasting sits underneath as
   fast cores + interop shims, not as a 1:1 rewrite of ~2,400 hand-written
   API objects.
3. **Correctness is the entry ticket.** In bioinformatics, fast-but-wrong is
   worse than slow: all outputs are cross-validated against Biopython
   (differential + fuzz testing) before a module ships.
4. **One dispatch interface, backends later.** Runtime CPU feature dispatch
   (baseline / AVX2 / AVX-512) from day one; GPU backends hang off the same
   interface if ever needed.

### Decisions

- **C/C++ core** (thin C ABI, `nanobind`/`pybind11` bindings) — chosen over
  Rust mainly for direct access to the C ecosystem: htslib, libdeflate,
  parasail, WFA2.
- **No CUDA for now** (decided 2026-10-02). Parsing is I/O-bound; CPU SIMD +
  threads takes the whole parsing win. GPU remains a deferred, optional
  backend idea for compute-kernel workloads (pairwise alignment, k-mer
  counting), never a hard runtime requirement.

## Repository layout

```
inventory/
  scan_biopython.py           # rerunnable API inventory scanner
  INDEX.md                    # per-package summary of the scan
  biopython_inventory.csv     # 3,595 public API objects, one per row
  biopython_inventory.json    # full dump incl. private members and __all__
bench/
  gen_data.py                 # seeded, byte-reproducible corpus generator
  run.py                      # benchmark runner (SeqIO vs naive vs pysam/pyfaidx)
  profile_paths.py            # cProfile/py-spy over the runner's own paths
  plot_g1.py                  # renders misc/g1-ranking.png from the ranking
  plot_profile.py             # renders misc/profile-share.png from the profiles
  requirements.txt            # pinned bench environment
  README.md                   # what the corpus contains and why it is honest
  profiling.md                # step 0.4: hotspots with % share
  targets.md                  # step 0.5: the ranked target list (Gate G1)
  data/                       # generated (gitignored) + MANIFEST.json checksums
misc/                         # rendered charts referenced by this README
```

Regenerate the inventory any time (updates automatically with your installed
Biopython version):

```sh
python3 inventory/scan_biopython.py
```

Regenerate the benchmark corpus — FASTQ 1M×150bp plus a gzip variant and a
test-sized prefix, a 100 Mbp genome FASTA, and edge/malformed sets:

```sh
python3 bench/gen_data.py
python3 bench/gen_data.py --verify   # re-check every SHA-256 in MANIFEST.json
```

Measure it — the runner refuses to report a time for an implementation whose
output does not match Biopython's:

```sh
python3 -m venv .venv && .venv/bin/pip install -r bench/requirements.txt
.venv/bin/python bench/run.py                     # the comparison table
.venv/bin/python bench/profile_paths.py --all     # where the time goes
.venv/bin/python bench/plot_g1.py                 # -> misc/g1-ranking.png
.venv/bin/python bench/plot_profile.py            # -> misc/profile-share.png
```

The numbers these produce are charted in
[Measured on the real corpus](#measured-on-the-real-corpus-phase-0) above, and
tabulated in [bench/targets.md](bench/targets.md); the profiles behind them are
in [bench/profiling.md](bench/profiling.md).

## Key inventory facts (Biopython 1.88)

- 284 modules, 156,755 lines of code; 3,595 public API objects.
- 1,112 of those objects are auto-generated `Restriction` enzyme classes —
  replaced by data + a parser here, never ported as classes.
- `pairwise2` is deprecated upstream; `cpairwise2` is already a C extension —
  evidence that upstream accepts C accelerators.

## Roadmap

Phase 0 (evidence) → Phase 1 (FASTQ/FASTA flagship core) → Phase 2 (expansion)
→ Phase 3 (ecosystem). Detailed steps, gates, decisions log and risk table:
see [PLAN.md](PLAN.md).

## License

TBD.