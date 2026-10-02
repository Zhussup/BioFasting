# BioFasting — Action Plan

Working document. Owner: @zhus. Created 2026-10-02.
See [README.md](README.md) for the pitch and principles.

## Current state

- [x] 0.1 API inventory of Biopython 1.88 (`inventory/`, 3,595 public
      objects, 284 modules, 156,755 LOC) + draft single-machine benchmarks
      showing a **7–10× gap** on hot paths.
- [x] 0.2 Seeded corpus generator (`bench/gen_data.py`, `bench/README.md`):
      FASTQ 1M×150bp + gzip variant + 10k prefix, genome FASTA, edge and
      malformed sets; byte-reproducible (SHA-256 counter-mode randomness,
      integer-only sampling) with a verified `MANIFEST.json`.
- [x] 0.3 Benchmark runner (`bench/run.py`): Biopython vs naive Python vs
      pysam/pyfaidx, median-of-5, correctness-gated — 20 rows, all validated.
- [x] 0.4 Profiling (`bench/profile_paths.py`, `bench/profiling.md`): a third
      to two thirds of every Biopython row is Python bytecode in `Bio/**`;
      charted as `misc/profile-share.png`.
- [x] 0.5 Ranked targets (`bench/targets.md`); **Gate G1 passed** — FASTQ/FASTA
      confirmed as the flagship from data.  Charted as `misc/g1-ranking.png`.

**Phase 0 is closed.**  All of the above is committed and reproducible from the
commands in [bench/README.md](bench/README.md); the numbers are the ones stored
in `bench/results/latest.json`, and `bench/targets.md` carries the run-to-run
variance so a re-measurement is not mistaken for a regression.

**Next:** Phase 1, step 1.1 — the C++ core scaffold.  No library code exists
yet, by design: the harness came first.

## Phase 0 — Evidence (weeks 1–2)

Goal: replace intuition with a ranked, reproducible list of targets.

| # | Step | Output | Done when |
|---|---|---|---|
| 0.2 | Dataset generator (seeded): FASTQ 1M×150bp, gzipped variant, multi-line FASTQ edge corpus, quality-encoding edge cases, genome-sized FASTA | `bench/data/` | anyone can regenerate the exact corpus from one command ✅ |
| 0.3 | Benchmark runner: Biopython `SeqIO` vs naive pure-Python vs installed C-backed libs (`pysam`, `pyfaidx`); median-of-N runs, MB/s | `bench/run.py`, `bench/results/*.json` | one command produces the full comparison table ✅ |
| 0.4 | Profiling of the slowest paths (cProfile + py-spy): interpreter overhead share per path | `bench/profile_paths.py`, `bench/profiling.md`, `bench/plot_profile.py` | top hotspots named with % share ✅ |
| 0.5 | Ranked targets: *workload → Biopython time → best alternative → stolen gap* | `bench/targets.md`, `bench/plot_g1.py` | ranking exists with evidence behind every row ✅ |

**Gate G1:** confirm flagship choice from data (expected: FASTQ/FASTA).
**Passed 2026-10-02.** The two largest measured gaps are both in that scope —
FASTQ parsing (7.1×, 2.74 s per million reads, 55.7% of it Python bytecode) and
FASTA random access (160 ms vs 3 µs per slice, because `SeqIO.index` has no
slice-capable index). Ranking and ordering: `bench/targets.md`.

## Phase 1 — Flagship: FASTQ/FASTA core (weeks 3–8)

Chosen 2026-10-02 (largest measured gap, widest audience, clearest
competitors). Scope boundary: sequence core only — **no SAM/CRAM**
(that is pysam's territory, we will not duplicate it).

| # | Step | Output | Done when |
|---|---|---|---|
| 1.1 | C++ core scaffold: `src/`, nanobind bindings, CMake build, `tests/`, GitHub Actions CI, cibuildwheel for linux-x86_64 / linux-aarch64 / macos-arm64 | buildable importable package locally | `pytest` green in CI on ≥1 platform |
| 1.2 | FASTQ parser kernel (plain + gz via vendored libdeflate), mmap-backed FASTA with index | first kernels behind the dispatch interface | parse output byte-identical to Biopython on the corpus |
| 1.3 | Sequence ops: reverse complement, GC content, k-mer counting; zero-copy numpy views of parsed buffers | op layer | verified vs Biopython on corpus |
| 1.4 | Interop shims: our records ↔ `Bio.SeqRecord` (both directions) | shim layer | existing `SeqIO`-based snippets work unchanged on our objects |
| 1.5 | Differential + fuzz tests vs Biopython (golden corpus, malformed files) | `tests/` corpus | corpus grows continuously; known upstream bugs (from CSV) are cases |
| 1.6 | Bench vs baseline, publish numbers | updated `bench/` targets table | ≥10× vs SeqIO on plain FASTQ, ≥1× vs pysam |

**Gate G2:** release 0.0.x — PyPI once wheels pass CI (local only before that).

## Phase 2 — Expansion (months 3–6)

- [ ] 2.1 Alignment kernels as wrappers behind our dispatch interface
      (parasail / WFA2-lib — reuse, don't rewrite SW).
- [ ] 2.2 Codon tables / Restriction as **data + parser** instead of 27k LOC
      of auto-generated classes.
- [ ] 2.3 Arrow interop; check compatibility story with polars-bio rather
      than colliding with it.
- [ ] 2.4 Remaining `SeqIO` formats hot-path-first, based on Phase 0
      ranking.

## Phase 3 — Ecosystem

- [ ] 3.1 Docs, tutorials, real-pipeline examples.
- [ ] 3.2 Contributions of C accelerators *into* Biopython itself (they
      already accepted `cpairwise2`) — win-win, grows trust.
- [ ] 3.3 Community building; target the "gold standard" position.

## Standing decisions (log)

| Date | Decision |
|---|---|
| 2026-10-02 | **No CUDA for now.** C/C++ core; GPU stays a deferred, optional backend idea behind the kernel-dispatch interface. GPU candidates later: pairwise alignment, k-mer counting — never parsing. |
| 2026-10-02 | **C/C++ core** over Rust — access to htslib / libdeflate / parasail / WFA2. |
| 2026-10-02 | **Benchmark-first** discipline: no kernel gets built without a measured target from Phase 0/1 benches. |
| 2026-10-02 | **Flagship = FASTQ/FASTA stack**; harness before core code. |
| 2026-10-02 | **CSV triage split:** agent prefills `verdict` (perf / parity / skip) with heuristics (deprecation flags, external-program wrappers, Restriction autogeneration), owner validates the contested rows only. |
| 2026-10-02 | **Repo artifacts in English or Chinese only** (docs, comments, docstrings); Russian is chat-only. |

## Risks and counters

| Risk | Counter |
|---|---|
| Silent data corruption in bio pipelines | differential/fuzz gating before every release (step 1.5) |
| Low adoption vs "just use Biopython" | shim strategy; drop-in at interop layer |
| Collision with pysam / polars-bio | explicit scope boundary (no SAM/CRAM); interop-first thinking |
| Scope creep (Biopython's 2,400-object surface) | flagship gate, verdict column, "data over classes" principle |
| Single machine / single maintainer | CI from day one; seeded public bench corpus |