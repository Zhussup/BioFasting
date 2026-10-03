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

- [ ] 1.1 C++ core scaffold (`src/`, nanobind, CMake, `tests/`, CI,
      cibuildwheel): **written and verified locally, not yet confirmed in CI.**
      The extension builds through scikit-build-core, imports, and reports its
      own build identity (`version`, compiler, flags, target arch) alongside a
      runtime ISA ladder; 16 tests pass; a real wheel (~61 KB) and sdist
      (~67 KB) both build, and each installs and passes the suite in a fresh
      virtualenv, including a wheel rebuilt *from* the sdist.  Those sizes are
      approximate on purpose: the sdist carries the docs, so it moves with
      every README edit.  What is missing is the second half of the Done-when:
      the workflows are written but have never executed, and that needs a push.
      Nothing here should be read as "CI is green".

**Phase 0 is closed.**  All of the above is committed and reproducible from the
commands in [bench/README.md](bench/README.md); the numbers are the ones stored
in `bench/results/latest.json`, and `bench/targets.md` carries the run-to-run
variance so a re-measurement is not mistaken for a regression.

**Next:** step 1.2 — the FASTQ parser kernel, the first thing to use the
dispatch interface step 1.1 put in place.  No library code exists yet, by
design: the harness came first, and now the scaffold it hangs off.

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
| 1.1 | C++ core scaffold: `src/`, nanobind bindings, CMake build, `tests/`, GitHub Actions CI, cibuildwheel for linux-x86_64 / linux-aarch64 / macos-arm64 | buildable importable package locally | `pytest` green in CI on ≥1 platform ✅ |
| 1.2 | FASTQ parser kernel (plain + gz via vendored libdeflate), mmap-backed FASTA with index | first kernels behind the dispatch interface | parse output byte-identical to Biopython on the corpus ✅ |
| 1.3 | Sequence ops: reverse complement, GC content, k-mer counting; zero-copy numpy views of parsed buffers | op layer | verified vs Biopython on corpus ✅ |
| 1.4 | Interop shims: our records ↔ `Bio.SeqRecord` (both directions) | shim layer | existing `SeqIO`-based snippets work unchanged on our objects ✅ |
| 1.5 | Differential + fuzz tests vs Biopython (golden corpus, malformed files) | `tests/` corpus | corpus grows continuously; known upstream bugs (from CSV) are cases ✅ |
| 1.6 | Bench vs baseline, publish numbers | updated `bench/` targets table | ≥10× vs SeqIO on plain FASTQ, ≥1× vs pysam — **measured 14.1× / 3.05×** ✅ |

**Phase 1 closed 2026-10-03.** All 297 tests pass with Biopython and the corpus;
the delivered numbers and the caveat attached to each one are in
`bench/targets.md` ("Phase 1 result") and in the run report
[long-run-night-time-report.md](long-run-night-time-report.md). One measured
limit is worth carrying forward rather than burying: the corpus FASTQ files have
8 distinct header widths across 1,000,000 records, so they are not a grid and
the zero-copy view does not apply to the file the FASTQ benchmark measures.

**Gate G2:** release 0.0.x — PyPI once wheels pass CI (local only before that).

## Phase 2 — Expansion (months 3–6)

- [ ] 2.1 Alignment kernels as wrappers behind our dispatch interface
      (parasail / WFA2-lib — reuse, don't rewrite SW).
- [x] 2.2a Restriction as **data + parser** instead of 27k LOC of
      auto-generated classes.  Done 2026-10-03: `biofasting.restriction` is
      1,088 enzymes as one class plus a table generated from the reference by
      `tools/gen_restriction_data.py`, which refuses to write anything unless it
      can rebuild every one of Biopython's own search patterns from the
      recognition sites.  `tests/test_restriction.py` compares every field,
      every flag, `search`, `catalyze` and `elucidate` against
      `Bio.Restriction` for all 1,088, on linear and circular sequences,
      including the sites that fall off the ends.  Closed 2026-10-03 by M10:
      `_print_format.py` carries `PrintFormat` and `restriction.py` carries
      `Analysis`, so the list, the count and the map are all produced here —
      the first two byte-identical over all 1,088 enzymes, the third identical
      line for line once the names on a line are sorted, which is as far as
      anything can go given that the reference's own map differs from run to
      run.  18 docstring examples are executed as tests.
- [x] 2.2b Codon tables as the same shape: data plus a lookup, not generated
      code.  Done 2026-10-03 by M11: `Bio.Data.CodonTable`'s 1,308 lines, about
      900 of which are 27 genetic codes written out one codon at a time, are
      `src/biofasting/_codon_tables_data.py` instead — 27 rows of NCBI's own
      64-character amino-acid strings, 5 KiB — read by `codon_table.py`, which
      fills all twelve of the reference's dictionaries at import.  The classes
      are behaviour and stay written out, in `_codon_table.py`.  `tools/
      gen_codon_tables.py` refuses to write unless the derived string rebuilds
      the reference's own table for all 27 codes, the source literal it read
      names from equals the table being run, and all 162 rebuilt tables agree
      with the reference's on every codon over their whole alphabet — 371,115
      codons asked one at a time, 82,760 of them refusals, compared by exception
      name.  Three codes whose stops are context-dependent carry the fact in a
      `dual` field that no string of that shape can hold.
- [x] 2.2c Substitution matrices as the same shape: the 30 NCBI data files plus
      the `Array` class that reads them, not generated code.  Done 2026-10-03 by
      M12: the thirty tables are `src/biofasting/_substitution_matrices_data.py`
      — 8,388 numbers, since all thirty are symmetric and each is stored as its
      lower triangle — and `substitution_matrices.py` loads them.  `Array` is
      behaviour and stays written out, in `_substitution_matrices.py`, including
      the `alphabet` property whose three refusals live in Biopython's
      `_arraycore.c` and were reimplemented from that source.  `tools/
      gen_substitution_matrices.py` refuses to write unless both parsers read all
      thirty of the reference's files identically, the rows rebuild the matrices
      entry for entry, and the module loads them back.  With this the triage's
      `data` verdict is discharged in full: 1,254 rows of 1,254.
- [x] 2.2d The `Bio.Seq` operation family completed: **translation**, the fourth
      kernel behind `ops-*` after reverse complement, GC and k-mers.  Done
      2026-10-03 by M13.  This one is not a data verdict but a contract —
      `Bio.Seq.translate` resolves its table three ways, answers `"X"` for a
      codon of legal nucleotides it does not name, warns about and refuses
      `to_stop` for a table whose stops are also residues, never checks
      `stop_symbol` at all, and has `cds` check the start, the length and the
      final stop before calling the first codon `M`.  So the kernel is
      table-driven (a 256-byte base-code table and a 64-byte amino-acid table
      built once per genetic code, in Python) and does the whole job or
      declines, exactly as the FASTQ fast path does; every rule above is applied
      in `_translate.py`.  The target was measured before the kernel existed
      (`bench/targets.md`, "Phase 2 ranking": 3.8× against the best alternative,
      a plain dict-per-codon loop) and the delivered row is **15.2×** against
      `Seq.translate` and 3.9× against that loop, 171 Mbase/s on the benchmark's
      100,000 reads.  `tests/test_translate.py` compares every codon of all 27
      genetic codes in both `to_stop` settings, randomised sequences over eight
      alphabets, and every way the reference refuses a call — and compares the
      kernel against the fallback directly, which is what caught a cache read
      through `AmbiguousCodonTable.__getattr__`'s delegation to the unambiguous
      table it wraps.
- [x] 2.2e The other half of Phase 0's rank 2: **random access into FASTQ**.
      Done 2026-10-03 by M14.  `SeqIO.index(path, "fastq")` stores a byte offset
      per record and re-parses the record from it on every fetch — at four short
      lines per record the parse *is* the fetch, 7.45 µs to return 150 bases.
      `FastqIndex` (`fastq.{hpp,cpp}`) is built by the scanner that already
      exists, so the build costs one scan, and it records where each field *is*
      rather than copying it: a fetch is a memcpy and `quality_length` does not
      touch the file.  Measured **23.4×** per access (0.32 µs against a 0.21 µs
      floor of a dict lookup and a slice) and 2.9× on the build.  The build's
      remaining cost is one hash table, and that was measured rather than
      assumed: 0.137 s of the 0.41 s is the floor for *any* name→record index
      over this file, because probes and key bytes are both random over 369 MB
      (two standalone micro-benchmarks in `bench/targets.md`), so the tuning
      stopped at an arena allocator instead of a hand-rolled table.  The
      milestone's real finding is not the speedup but that Biopython ships two
      FASTQ grammars that disagree about which files are legal: this index
      follows the parser, and both halves of the divergence are asserted in
      `tests/test_fastq_index.py`.
- [x] 2.2f The form the file is actually in: **random access into a compressed
      FASTQ**.  Done 2026-10-03 by M15.  `open_fastq_index` indexed plain files
      only, and said so; production FASTQ is gzipped, so this was the one hole
      the package documented about itself.  It now maps the `.gz`, inflates it
      once with libdeflate, and indexes the inflated bytes — the same one-pass
      scanner over a different span, so a fetch from a 172 MB gzip costs the
      0.33 µs it costs from the 369 MB it inflates to.  Measured **3.2×** on
      the build and **~800×** per fetch against `SeqIO.index`, both digest-gated
      — and the comparison needed a BGZF copy of the corpus bytes, because
      `SeqIO.index` refuses a plain gzip FASTQ outright ("Gzipped files are not
      suitable for indexing, please use BGZF instead"), which is the format most
      FASTQ files are in.  The reference's 272 µs per access is the same design
      gap as Phase 0's rank 2 in its worst case: a seek into a BGZF stream
      invalidates the block cache, so every access pays a decompressor for a
      block it will not reuse.  BGZF is gzip with extra fields, so both formats
      are one code path; `tests/test_fastq_index.py` builds a BGZF stream by hand
      to prove it rather than to assume it.
- [x] 2.2g **The wheel carries its own inflate.**  Done 2026-10-03 by M16, and it
      is not a kernel but the last thing in the tree that was knowingly wrong.
      Every `.gz` row was produced by an extension linked against Debian's
      `libdeflate.so.0`: the wheel imported only where libdeflate was installed,
      and CI could not see it because CI installed libdeflate too.  The host's
      static archive does not fix it — measured, `relocation R_X86_64_PC32
      against symbol libdeflate_x86_cpu_features can not be used when making a
      shared object; recompile with -fPIC`.  So libdeflate 1.23 is vendored whole
      in `third_party/libdeflate` by `tools/vendor_libdeflate.py`, which refuses
      to write unless the tarball's SHA-256, the version in its header and
      upstream's own `LIB_SOURCES` list all agree, and compiled here with
      upstream's own `-O2 -DNDEBUG` rather than this project's `-O3` so that the
      gz rows stay comparable with the ones already published.  The rows do not
      move (6.07× on `fastq-gzip`, 3.1× on the compressed index build, ~890× per
      compressed fetch), `readelf -d` on the installed extension shows no
      libdeflate, and the price is 32 KB of extension size.  Two artifact tests
      hold the claim to the file that gets installed rather than to the build
      system; `-DBIOFASTING_VENDORED_LIBDEFLATE=OFF` goes back to the host's
      library and is exercised in CI, because a supported configuration that
      nothing builds is a configuration that rots.
- [x] 2.2h **The other half of `Bio.SeqUtils`: measurement instead of operation.**
      Done 2026-10-03 by M17.  2.2a–2.2g changed sequences; this one measures them,
      and the split is not cosmetic, because a measurement is where exactness
      stops being free.  Six kernels in `src/core/sequtils.{hpp,cpp}`
      (`gc123`, `gc_skew`, `gc_counts`, `molecular_weight_mass`, `gcg`,
      `crc64`) and two Python modules over them, `biofasting/sequtils.py` and
      `biofasting/checksum.py` — `gc_fraction`'s two missing modes arrive with
      `gc_counts` rather than as a third kernel.  The rule that made the module
      work is **a kernel counts, Python decides**: `gc123` returns twelve
      *counts* and not four percentages, because the reference divides three
      times and raises on the fourth — `GC123("")` is a `ZeroDivisionError`, an
      empty *frame* is the integer `0`, and `gc_fraction("")` is the integer `0`
      rather than `0.0`, which is why the C++ `gc_fraction` returns `-1.0` as a
      sentinel meaning "nothing counted" and the Python layer turns that into
      the integer.  None of that survives a `double`, and the difference is
      observable from the outside, so it is not a detail.
      Two more places where the exact answer is not the obvious one.
      `molecular_weight` is `sum(weight_table[x] for x in seq)`, and CPython 3.12
      sums floats with Neumaier compensation: a naive running total mismatched
      the reference on **5,337 of 20,000** random draws (26.7%), so the kernel
      implements the compensated update step for step.  And `.upper()`/`ord()`
      are *code-point* operations, not byte operations — `gcg("é")` is 201 and
      `gcg("ß")` raises `TypeError: ord() expected a character` — so `gcg`
      declines on any byte ≥ 0x80 and falls back to the reference's own loop,
      while `crc64` needs no decline at all because `ord(c) & 0xFF` *is* the
      Latin-1 byte; the asymmetry is written down in both files.  The
      differential layer earned its keep: `GC_skew` was wrong in C++ — `g - c`
      on unsigned counters wraps to `2^64 - d`, so a window with more C than G
      returned a huge positive number instead of a negative fraction (the
      smallest case, a window of `CC`, is `2^63 - 1`, printed
      `9.223372036854776e+18`) — and it was the differential sweep, not the
      hand-written cases, that found it: replayed against the shipped formula
      that sweep is wrong on **884 of its 1,200 draws**, starting with the
      second.  Delivered: **33.2×** on `GC123`
      (0.5982 s → 0.0180 s per 20,000 reads) and **452×** on one 40 Mbp contig
      (8.1026 s → 0.0179 s), 27.3× on `crc64`, 24.0× on `gcg`, 4.5× on mass, 4.6×
      on the six-frame text, 4.2× on `GC_skew`, 2.2× on the
      translate→`seq3`→`seq1` round trip.  Two rows are below 2× and both are
      left in: `CodonAdaptationIndex.calculate` is 1.0× because it is a faithful
      Python port and no kernel was written for it (59 µs per read is a measured
      target that has not been spent), and the pure-Python `GC_skew` baseline is
      1.0× against the reference because `str.count` over windows is already the
      right algorithm — the kernel's win there is memory, not arithmetic.
      `CodonAdaptationIndex.optimize` is deliberately absent: it returns a
      `Bio.Seq.Seq`, and this package has no `Seq` type to return.  The ranking
      pass also withdrew one of its own numbers — `GC123` on chr1 was published
      as 90.674 s and is 8.10 s when measured under `build_info()` — which is
      recorded in `bench/targets.md` rather than quietly corrected.  It also
      corrected the triage: `Bio.SeqUtils.xGC_skew` was filed `perf`, and it
      imports `tkinter`, draws two curves and returns `None`, so no measurement
      could ever have produced a kernel from it.  It is a `skip` now, on a sixth
      mechanical reason, and the counts are 610 `perf` / 1,579 `parity` /
      1,254 `data` with 152 skipped — 611 `perf` before.  **Open:** those 610 and
      1,579 rows, plus PLAN 2.1, 2.3 and 2.4, none of which has a measured target
      yet.  Tests: 585 pass, 80 of them new.
- [x] 2.2i **`Bio.SeqUtils.ProtParam`: the per-residue half of the module.**
      Done 2026-10-03 by M18.  This is where `Bio.SeqUtils` stops being a
      collection of scans and becomes a class.  `ProteinAnalysis` is
      `biofasting/protparam.py` over four kernels in `src/core/protparam.{hpp,cpp}`
      — `protein_scale_scores`, `flexibility_scores`, `instability_index_sum` and
      `count_residues` — with `IsoelectricPoint` in its own module beside it and
      the 32 scales generated into `biofasting/_protparam_data.py` by
      `tools/gen_protparam_data.py`, which refuses to write unless every table
      equals the reference's value for value *and type for type*.
      The type is not cosmetic, and that is the finding of this milestone.
      `gravy` is `sum(selected_scale[aa] for aa in sequence)`, and CPython 3.12+
      compensates a float item in `sum()` and adds an `int` one *plainly*; **ten**
      of the twenty-eight gravy scales are written with integer values, and on
      four of them — Parker, GoldSack, Engelman, Roseman — a compensated sum over
      the same twenty numbers answers differently on **4% to 29%** of random
      proteins.  The count had been written down as nine and three; it was
      measured, corrected in every file that carried it, and the generator now
      probes four hand-found sequences (`"IYAR"` under Parker, `"HSV"` under
      Engelman, `"TNRC"` under GoldSack, `"THTQG"` under Roseman) on every run, so
      a table that had coerced `2` to `2.0` cannot pass.
      The same care lands in three other places, and each one is the reference's
      own behaviour rather than a tidier version of it: `flexibility` produces
      `size - 9` windows and reads its middle at index **5**, so index 4 is never
      read and index 5 is read twice, and folding the four mirror pairs into one
      symmetric weight vector — the obvious rewrite — differs in the last bit on
      **200 of 200** proteins; `molecular_weight` and `gravy` are `sum(...)` and
      so Neumaier-compensated, while `instability_index` is `score += value` and
      is not, so the module needs **two** summations and neither kernel may be
      borrowed for the other (shortest separating sequence, `"EPCM"`); and
      `protein_scale` *warns and skips* a window term rather than raising, which
      no kernel can reproduce because the message goes to Python's `stderr`.
      Delivered, digest-gated on 400 generated proteins / 192,800 residues, a
      fresh object per call, median of nine: **43.8×** on `instability_index`,
      **37.7×** on `protein_scale`, **36.9×** on `flexibility`, 9.3× on `gravy`,
      and — after a fourth kernel was built for the method the ranking pass had
      passed over — **5.9×** on `count_amino_acids`, 5.9× on
      `molar_extinction_coefficient`, 3.8× on `amino_acids_percent`, 3.4× on
      `aromaticity`, 3.0× on `charge_at_pH` and 2.8× on
      `secondary_structure_fraction`.  That fourth kernel is the one worth
      reading about: `count_amino_acids` is twenty `str.count` calls, and it is
      the first statement of **six** other methods, which sat at 0.98×–1.00×
      until it was written.  It is also the only kernel here that never declines
      — `str.count` is not interested in what it does not find, so a residue
      outside the twenty is a zero and not an error.
      `isoelectric_point` is 1.4× and is reported as it stands.
      **Open:** `MeltingTemp` is measured but not claimed (`Tm_NN` on a 60-mer is
      23.40 µs, with no rewrite that beats it); `ProteinAnalysis.optimize` is
      absent because it returns a `Bio.Seq.Seq` and this package has no `Seq`
      type to return; those 610 `perf` and 1,579 `parity` rows; PLAN 2.1, 2.3 and
      2.4.  Tests: 629 pass, 44 of them in `tests/test_protparam.py`.
- [ ] 2.3 Arrow interop; check compatibility story with polars-bio rather
      than colliding with it.
- [x] 2.4 Remaining `SeqIO` formats hot-path-first, based on Phase 0
      ranking.  **Delivered 2026-10-03** for the three formats that carry
      essentially all of the world's sequence data — GenBank, EMBL and
      SwissProt — as one pass over the mapped file in C++ (`FlatFileIndex`,
      [`src/core/flatfile.hpp`](src/core/flatfile.hpp)), reached from Python as
      `biofasting.open_genbank()` / `read_genbank()`.  It reproduces four
      fields of a record — `id`, `name`, `description` and the sequence — and
      is **not** a drop-in for `SeqIO.parse`: no FEATURES table, no annotations
      dict, no `SeqRecord`.  Digest-gated on all four fields against the
      reference on every record of all six corpus rows, **11,900 records and
      24.7 MB**, and measured in `bench/bench_genbank.py`: **3.52 µs** per 1 kb
      GenBank record against 43.32 (**12.3×**), **3.22** on the bare row
      against 25.09 (7.8×), **1.86** per SwissProt record against 26.95
      (**14.5×**), at 396–585 MB/s across record sizes.  It beats the
      pure-Python rewrite on every row while producing *more* than it does
      (that rewrite extracts ids and sequences only), so C was worth writing —
      but the rewrite was already 5–7× the reference, which is where most of
      the distance came from.  The reader is deliberately **stricter** than the
      reference where the reference warns and guesses: a wrongly indented
      GenBank `ORIGIN` line, a `CONTIG` line inside the sequence block and a
      duplicate id all raise rather than produce a sequence that looks parsed.
      **Measured earlier the same day** by the sixth ranking pass
      (`bench/targets.md`), and the Phase 0 ranking's answer — that every
      format left in `Bio.SeqIO` is `parity`, "no large win expected" — does
      not survive it.  GenBank, EMBL and SwissProt are generated by
      `bench/seqio_corpus.py`, emitted by the reference's own writer where one
      exists and by hand where it does not (SwissProt has no writer upstream),
      and held to their promises by `tests/test_seqio_corpus.py` before anything
      is measured on them.  The reference costs **48.6 µs** per 1 kb GenBank
      record and **46.8** per EMBL one, against **28.4** for the same record
      with no features — where a whole-file C-speed scan costs 1.5 µs and
      building the same `SeqRecord` tree from scalars costs 1.0.  Two-thirds to
      nine-tenths of the time is the reference's own line-by-line logic, which
      is the part C removes, and the headroom is **4.5×–11×**.  A rewrite that
      extracts ids and sequences and nothing else — verified identical on every
      row, and deliberately not a replacement — already reaches 6–13×.
      **`GenBankIterator`, `EmblIterator` and `SwissIterator` are therefore
      `perf`, not `parity`**; the three writers stay `parity`, unmeasured, and
      so do the formats no pass has ranked.  The delivered numbers sit beside
      that ranking in [`bench/targets.md`](bench/targets.md), including the one
      place the ranking's own reference column did not reproduce.  **Delivered
      later the same day (M21)** for the row that statement left open: the
      **FEATURES table** of GenBank and EMBL — feature key, location and
      qualifiers, with positions preserved rather than flattened — as our own
      `Feature` value type in [`src/biofasting/seqfeature.py`](src/biofasting/seqfeature.py)
      and an `interop` conversion to `Bio.SeqFeature`, over
      [`src/core/location.{hpp,cpp}`](src/core/location.hpp) and
      [`src/core/feature.{hpp,cpp}`](src/core/feature.hpp).  Digest-gated on
      type, location and qualifiers on all 12,000 features of the four corpus
      rows that carry them: **9.82 µs** per 1 kb GenBank record against 48.83
      (**5.0×**), 10.00 against 47.80 on EMBL (4.8×).  The measurement M19
      named is the one quoted, re-derived: **7.62 µs** per feature against the
      reference's re-measured **23.71** (recorded 20.2, and that did not
      reproduce — the run notes say so), so the 16 µs M19 attributed to reading
      the location string and the qualifier block is **2.1×**, and the
      conversion costs 4.79 against the reference's own 4.2 for the object —
      parity, so the win is entirely in the reading.  A location the reference's
      parser rejects is reproduced as the behaviour it is (`parser_error`,
      location `None`, warning), while a table this reader cannot reproduce
      **raises** with the byte offset; warnings are reported by the kernel and
      spoken by the Python layer, in the reference's own order.  SwissProt `FT`
      is a second grammar and is out of scope.  **Delivered later the same day
      (M22)** for the row *that* statement left open: the **annotations dict** —
      the header keys of `SeqRecord.annotations`, as plain values, in
      [`src/core/annotations.{hpp,cpp}`](src/core/annotations.hpp), reached as
      `biofasting.read_annotations()` with `to_reference()` for Biopython's own
      `Bio.SeqFeature.Reference`.  The kernel walks the header's lines and hands
      back the keys it created **in the order it created them**, because the
      reference inserts a key when the line stating it is consumed — a `COMMENT`
      above the first `REFERENCE` puts `comment` before `references` — and
      membership in that order is also what says a key *exists*, so an empty
      value and an absent key stay apart.  The two formats disagree about the key
      set and the reader is held to both: EMBL's thin header has six keys and
      GenBank's nine, the missing three being `date`, `source` and `keywords`.
      Digest-gated on every key, its position and every value, on the four
      corpus rows of the two bare/annotated pairs: **4.48 µs** per 1 kb GenBank
      record against 25.25 (**5.6×**) and 4.13 against 24.72 on EMBL (6.0×); the
      marginal numbers answer the recorded target — **4.61 µs** of GenBank
      header against the reference's re-measured 32.94 (recorded 31.9), and
      **5.37** against 17.03 on EMBL (recorded 17.9), so **6.9×** and **3.3×**.
      `gi` is reproduced rather than refused, because the files in circulation
      carry it; the refusal contract is kept for a key that is rare *and*
      unreproduced (`NID`, `PID`, `DBSOURCE`, `SEGMENT`, a structured comment, a
      pre-229.0 `LOCUS`), and a refusal is per record rather than per file.
      SwissProt's header is `Bio.SwissProt`, a second grammar, and refuses.
      **Delivered later the same day (M23)** for the row *that* statement left
      open: the **writers**, and they were measured before they were written.
      The seventh ranking pass ([`bench/targets.md`](bench/targets.md)) asked the
      triage's `parity` verdict a timer's question — three formats, three floors
      deep, `reference` against a byte-identical pure-Python `rewrite` against
      the raw `payload` — and got three different answers: the FASTA reference is
      already at the floor a Python rewrite can reach (2.18 µs a record against
      its own rewrite's 2.00, so the only headroom is its eighteen `write` calls
      a record), FASTQ is **4.81×** reachable without C at all (the reference
      walks the quality list through a dict, a base at a time), and QUAL is the
      worst writer in the family in absolute terms — **15.40 µs** for 150 bases,
      103 ns a base, because every score goes through `"%i" % round(q, 0)` and
      the lines are packed by popping from the front of a list.  The kernels are
      [`src/core/write.{hpp,cpp}`](src/core/write.hpp), reached as
      `write_fasta`, `write_fastq` and `write_qual` (plus `phred_to_sanger` for
      the one place a `SeqRecord` becomes the string they take), each building
      the whole file in one buffer so the handle is crossed **once** rather than
      once a line.  Held to the pass's targets on the same corpus and gated over
      the whole file against `SeqIO.write`: **0.46 µs** a record on 1 kb FASTA
      against 2.12 (4.6×), **0.43** on 150 bp FASTQ against 5.23 (12.1×),
      **0.93** on 150 bp QUAL against 15.22 (16.4×), with the encoder alone at
      **9.3×**.  The reference's awkward branches are reproduced rather than
      tidied: the unwrapped FASTA branch writes a newline where a wrapped one
      writes nothing, and QUAL's width of one to five takes the `pop(0)` branch
      the kernel does not have, in Python.  The flat-file writers stay `parity`
      and unmeasured — a complete GenBank record means *writing* the annotations
      dict and the FEATURES table, which M21 and M22 only read, so that target is
      recorded and deferred rather than started.
      **Open:** the delivered reader is 2.1× above the whole-file scan on the
      bare row and that is the remaining headroom; the flat-file writers are
      measured-and-deferred, and the conversion in front of the sequence writers
      (`from_seqrecord`, 2.5 µs a record) costs more than the writers themselves.

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
| 2026-10-03 | **Python floor is 3.10, not 3.9.** nanobind 3.x declares `Requires-Python >=3.10`, and 3.9 has been EOL since Oct 2025. Holding nanobind at 2.x to keep a dead interpreter supported was the worse trade. The README badge, `requires-python` and the CI matrix all say 3.10–3.13; the badge was corrected rather than left claiming support nothing tests. |
| 2026-10-03 | **Two nanobind defaults overridden** in `CMakeLists.txt`, both named and reported by `build_info()` so a surprising benchmark number can be traced to a build decision instead of argued from memory: `NOMINSIZE` (nanobind adds `-Os` on top of Release's `-O3`, backwards for a throughput library) and `PROTECT_STACK`, plus an explicit `-fstack-protector-strong` so the security claim is true by construction rather than by distro default. Each is a `-DBIOFASTING_*=OFF` away from nanobind's behaviour. |
| 2026-10-03 | **Wheels build on tags only, never per push**, and musllinux is kept even though it doubles the Linux matrix — Alpine is a common base for bioinformatics containers and the matrix only runs on a tag. The per-push signal is `ci.yml`; `wheels.yml` is release infrastructure. |
| 2026-10-03 | **Build outputs never enter git.** A commit made while the scaffold was being built swept in 43 files under `build/` — object files, the compiled `.so`, CMake's compiler probe, two `a.out` binaries, and a `CMakeCache.txt` holding absolute local paths. `.gitignore` now covers `build/`, `dist/`, `wheelhouse/`, `_skbuild/`, `*.so` and `*.egg-info/`, and a build directory is treated as disposable: `rm -rf build/` must never cost anything. The tracked copies still have to be removed from the index once (`git rm -r --cached build/`); deleting them from the worktree is not enough. |
| 2026-10-03 | **The inventory gets its verdict column, and `skip` must be earned mechanically.** The triage decided on 2026-10-02 was never run: the CSV has no `verdict` column and no triage code ever existed in the history, so Phase 1 was ordered by `bench/targets.md` instead. `inventory/triage.py` now prefills it into `inventory/biopython_triage.csv`, one `reason` per row — **3,444 of 3,595 rows are on the remaster** (611 perf, 1,579 parity, 1,254 data) and 151 are skipped. The vocabulary gains `data` for what the plan already calls "data + parser" (Restriction, codon tables, substitution matrices). Skip requires a mechanical reason — deprecated upstream, a network request, launching an external program, already C upstream, third-party numerics — and anything arguable is written down as work. The first draft broke that rule itself, skipping the local-file parsers filed under the network packages; that is recorded in `inventory/TRIAGE.md` rather than quietly fixed. One row moved later the same day: `Bio.SeqUtils.xGC_skew` left `perf` for a sixth skip reason, "it draws a window instead of returning a value" — the count became 610/1,579/1,254 with 152 skipped after M17 re-derived the row by hand. |

## Risks and counters

| Risk | Counter |
|---|---|
| Silent data corruption in bio pipelines | differential/fuzz gating before every release (step 1.5) |
| Low adoption vs "just use Biopython" | shim strategy; drop-in at interop layer |
| Collision with pysam / polars-bio | explicit scope boundary (no SAM/CRAM); interop-first thinking |
| Scope creep (Biopython's 2,400-object surface) | flagship gate, per-object verdict column (built 2026-10-03, `inventory/TRIAGE.md`), "data over classes" principle |
| Single machine / single maintainer | CI from day one; seeded public bench corpus |