# bench/ — the benchmark corpus and harness

Phase 0 of [PLAN.md](../PLAN.md) exists to replace intuition with a ranked,
reproducible list of targets.  That requires an input corpus that everyone can
regenerate byte-for-byte, so a number quoted in one machine's report means the
same thing as a number quoted in another's.

## Regenerating

```sh
python3 bench/gen_data.py          # full corpus  -> bench/data/   (~30 s, ~650 MB)
python3 bench/gen_data.py --quick  # CI smoke set -> bench/data-quick/
python3 bench/gen_data.py --verify # re-check existing artifacts against MANIFEST.json
```

The generated data is **not committed** (`bench/data*/` is gitignored): it is
large, and it is cheaper to regenerate than to store and clone.  `MANIFEST.json`
records a SHA-256 for every artifact plus the parameters that produced it.

## Why the corpus is byte-reproducible

Reproducibility here is a contract, not an aspiration — it is what lets a
benchmark result be re-derived years later from the same command.

- Randomness is SHA-256 in counter mode keyed by `(seed, domain)`, so the byte
  stream does not depend on the platform, Python version or CPU.
- Every sampling decision is integer arithmetic; no float ever reaches a
  threshold, so nothing can drift with libm differences.
- Each dataset draws from its own domain string, so changing one parameter does
  not shift the bytes of unrelated datasets.
- `.gz` members are written with `mtime=0` at a fixed compression level.
  The compressed *container* can differ by a byte or two across zlib versions,
  so the manifest hashes the uncompressed content and `--verify` decompresses
  before comparing.

Verified: two independent full runs (`bench/data/` and a second output
directory) produced byte-identical `MANIFEST.json` and passed `--verify`.

## Contents

| Artifact | Size (default) | Purpose |
|---|---|---|
| `fastq/reads_1m.fastq` | 369 MB, 1,000,000 × 150 bp | the flagship parse benchmark |
| `fastq/reads_1m.fastq.gz` | 172 MB (2.15×) | the compressed-path benchmark |
| `fastq/reads_10k.fastq` | 3.7 MB, 10,000 × 150 bp | exact prefix of the large file; unit tests |
| `fasta/genome.fasta` | 102 MB, 100 Mbp in 5 records | genome-scale scanning and indexing |
| `fasta/genome_1mb.fasta` | 1 MB, 1 Mbp in 1 record | test-sized indexed access |
| `edge/fastq/*` | tiny | legal but awkward FASTQ |
| `edge/fasta/*` | tiny | legal but awkward FASTA |
| `malformed/*` | tiny | inputs a correct parser must reject cleanly |

Total: 476 MB of content, 647 MB on disk.

Two input families are **generated, not stored**: the FASTA/FASTQ reads come from
`gen_data.py`'s SHA-256 keystream (the files above are its output, kept on disk so
a run does not have to re-derive them), and the protein and flat-file corpora are
derived in memory.  `run.py` builds the protein corpus for `protparam-*` through
`make_protein_corpus()`, and `seqio_corpus.py` builds the GenBank, EMBL and
SwissProt files that `rank_seqio.py` measures.  All of them are integer-only
counter-mode constructions with no float threshold anywhere, so the same command
yields the same bytes and the same residues on any machine.

`seqio_corpus.py` is the one place the corpus is not wholly ours: GenBank and
EMBL are emitted by `Bio.SeqIO.write`, because those two have a writer and a
hand-written INSDC record is a week of finding out where the columns are.
SwissProt has no writer upstream, so that emitter is written out here.  The
module says so at the top rather than leaving it to be discovered.

### What "realistic" FASTQ means here

A corpus of uniform random bytes would compress at ~1.25× and would flatter any
parser.  The default `--profile realistic` instead models what production data
actually looks like:

- **Quality decay.** Read cycles split into 5' / middle / 3' bands with falling
  modal quality (measured: Q34.8 → Q32.2 → Q27.6).  Trimming tools exist because
  of this decay, so a parser that mangles quality order fails here.
- **3' adapter read-through** on 15% of reads (Illumina TruSeq), the reason
  adapter trimming exists at all.
- **PCR duplicates** at 5%: an exact copy of a recently seen molecule.
- **GC content** 45%, with N-runs and a soft-masked region in the genome FASTA.

One honest caveat: this corpus compresses to about **2.15×**, where production
FASTQ reaches 3.5–4× (real quality strings are correlated along the read; these
draws are independent).  The gzip payload is therefore *harder* to decompress
than real data, which biases against, never in favour of, a decompressor's
measured win.  See the note on the `.gz` entry in `MANIFEST.json`.

### Edge cases worth knowing

- `edge/quality_contains_at_and_plus.fastq` — a quality string made of `@`.
  Any parser that scans for `@` to find the next header corrupts this record.
- `edge/empty_read.fastq` — zero-length sequence and quality are legal.
- `edge/phred64.fastq` — legacy Phred+64.  Not an error anywhere: Biopython 1.88
  reads the first base as **Q71 instead of Q40** and reports nothing.  Wrong
  numbers that look right are the failure mode this project must not ship.
- `edge/blank_lines.fastq` — stray blank lines between records.  Biopython
  tolerates them; whether we should is a compatibility decision, not obviously
  "malformed".
- `edge/fasta/comment_lines.fasta` — Biopython's strict `fasta` parser *rejects*
  a leading `;` comment; `fasta-pearson` accepts it.  Another deliberate call.

The `malformed/*` notes record Biopython 1.88's exact reaction (error type and
message) to each file.  Those notes are the expected-error baseline the
differential tests in step 1.5 will assert against.

## Benchmark runner

```sh
python3 -m venv .venv
.venv/bin/pip install -r bench/requirements.txt
.venv/bin/python bench/run.py              # full corpus, median of 5
.venv/bin/python bench/run.py --list       # workloads and implementations
.venv/bin/python bench/run.py --only fastq-gzip --repeat 9
.venv/bin/python bench/run.py --only fasta-random --random-accesses 20
```

Progress is printed a row at a time, so a long run is visibly moving rather
than silent until the end.

`bench/run.py` measures, on the corpus above:

- `fastq-plain` / `fastq-gzip` — `Bio.SeqIO.parse` (what people write),
  `FastqGeneralIterator` (Biopython without `SeqRecord` construction), a naive
  pure-Python parser, `pysam.FastxFile`, plus decompression-only baselines
  (`gzip` and raw `zlib`) for the compressed case;
- `fasta-scan` — the same shape for plain multi-record FASTA;
- `fasta-random` — fixed 150 bp slices, pyfaidx vs `SeqIO.index`, reported as
  **latency per access** rather than MB/s (see below);
- `ops-revcomp`, `ops-gc` — the in-memory operations from the draft table in the
  top-level README;
- `ops-translate` — translation over the same reads, against `Seq.translate` and
  a dict-per-codon baseline;
- `ops-gc123`, `ops-gc-skew`, `ops-molecular-weight`, `ops-crc64`, `ops-gcg`,
  `ops-protein-roundtrip`, `ops-cai-calculate`, `ops-six-frame` — the
  `Bio.SeqUtils` measurement functions, plus `sequtils-gc123-genome`: the same
  call on one 40 Mbp contig, where the reference takes eight seconds and the
  per-read ratio stops describing anything.  These rows run over a *smaller*
  slice of the same reads (`SEQUTILS_READS`, 20,000) because their reference
  rows cost three orders of magnitude more per call, and the count is in each
  row's description so that a row over a subset is never mistaken for one over
  the whole corpus.
- `protparam-*` — `Bio.SeqUtils.ProtParam.ProteinAnalysis`, eleven methods, over
  a corpus of 400 generated proteins (200 × 150, 120 × 300, 75 × 1,024, 5 ×
  10,000 residues; SwissProt residue frequencies, `PROTEIN_SEED`), **a fresh
  object per call** because `count_amino_acids` memoises into the instance and a
  row timed on a used object measures a cache hit.  Four of the eleven also carry
  the best pure-Python rewrite, kept because it is the number that says whether a
  kernel was worth writing at all — 1.10× to 1.84×, against 36.9×, 43.8× and
  37.7× for the three that got one.

The flat-file corpus is measured by its own script rather than by `run.py`,
because a corpus built in memory at run time is not a file a reader can map.
`bench/bench_genbank.py` writes each generated row to disk and measures three
implementations on it — `Bio.SeqIO.parse`, the smallest pure-Python parser there
is (`rank_seqio.minimal_parse`, **imported rather than copied** so the two cannot
drift), and `biofasting.read_genbank` — under the same all-or-nothing digest gate
as everything else, on the four fields the kernel produces rather than on the
two the rewrite does.  `bench/rank_seqio.py` is the ranking pass that came
before it: the same corpus split into the reference's own time, the object floor
and a whole-file C-speed scan, which is where the target for the kernel came
from.  `bench/bench_features.py` is the other half: the same corpus rows, but the
gate is over every **feature** rather than over the record — type, location as
canonical nested tuples (kind and edges, not objects, because Biopython's
positions compare by integer value and an object comparison cannot see a `<` that
was flattened), qualifiers and status — and it measures three paths: the
reference, the feature reader, and the reader plus `to_seqfeature`, so the column
that is this milestone's work is separate from the column that is object
construction.  It also re-derives M19's per-feature target as the difference
between two rows that differ only in feature density, per *record* and not per
file.  `bench/bench_annotations.py` is the third half: the two bare/annotated
pairs the M22 target was measured on, gated on the **annotations dict** — every
key, its position in the insertion order and every value, references through
`to_reference` so both sides are the reference's own objects.  The key order is
in that digest and has to be, because `dict.__eq__` is order-blind and would pass
a reader that emitted a fixed key list.  Its marginal block subtracts two
per-*record* times within a pair, and the timed unit is the parse alone on both
sides: the reduction that the gate compares costs more on the annotated row than
on the bare one, so a timing that included it would put that difference inside
the delta.

The other direction — writing — has its own pair of scripts, because the triage
filed every writer as `parity` on an assumption nobody had timed.
`bench/rank_writers.py` is the seventh ranking pass: three formats and three
floors deep (`reference`, a byte-identical pure-Python `rewrite`, and the raw
`payload`), gated byte-for-byte before any time is quoted.  It found three
different answers — the FASTA reference is already at the floor a Python rewrite
can reach, FASTQ is 4.81× reachable without C at all, and QUAL is the worst
writer in the family in absolute terms, 103 ns a base — and the flat-file
writers get a per-method attribution instead of a twin, because a byte-identical
GenBank record is a page of column rules.  `bench/bench_writers.py` is the
delivered kernels held to the targets that pass recorded: `write_fasta`,
`write_fastq` and `write_qual`, gated over the whole file against
`SeqIO.write`, with the `SeqRecord` path (`from_seqrecord`) and the quality
encoder timed as their own columns so the conversion in front of a writer is
never mistaken for the writer.  The scripts that measure on disk —
`bench_genbank`, `bench_features`, `bench_annotations` and `bench_writers` —
write their corpus into `build/bench-cache/seqio/` —
gitignored and disposable, a cache, never a tracked artifact.  (It was
`bench/data/seqio/` first, which `tests/test_corpus.py` correctly failed: it walks
`bench/data/` and requires the walk to equal the generator's manifest, so a file
written there by anything else is a failure by design.)

Each workload also carries a **`biofasting`** row once the package is
importable, so the thing being built is measured by the same harness, against
the same reference, under the same digest gate as everything else.  The runner
does not need it: without an installed package the rows are absent, the rest of
the table is unchanged, and the run is still the Phase 0 comparison.

A time is only meaningful against the binary that produced it, so the results
JSON records `biofasting.build_info()` next to the host facts — compiler,
C++ standard, optimisation flags, the SIMD rung `seqops_level()` actually
dispatched to, and what `best_cpu_level()` says the CPU could support.  A
machine with AVX-512 can be running the scalar kernels, and a number quoted
without the rung behind it is a number about an unknown program.

### The zero-copy views are measured elsewhere, on purpose

`bench/grid.py` measures the `FastaGrid`/`FastqGrid` arrays, and it is a
separate file because the runner's model of an implementation — yield each
record, then do `len()` on it — cannot measure a view at all.  A grid array is
two integer fields over a mapped file: the parse harness never touches the
bytes, so the file is never read, and the row would win by doing nothing.

`grid.py` therefore measures the three things a caller buys — the array, an
operation over it, and the RSS the process pays to hold it — and gates the
array's *content* against `SeqIO.parse` before quoting any time.  It also
measures why the corpus FASTQ files have no grid: the width of their headers.

### Random access is measured in latency, not throughput

`Bio.SeqIO.index` has no slice-level random access: an index entry is a
whole-record byte offset, so fetching 150 bp out of the middle of a record
parses **the entire record** and builds a `SeqRecord` from it.  On this corpus
`chr1` is 40 Mbp, which makes one 150 bp fetch cost ~160 ms — and the naive
`--random-accesses 20000` that this harness first used meant 20,000 × 0.16 s ≈
50 minutes *per pass*, six passes per run.  That is what the default access
count and the per-access column exist to prevent.

So the workload reports milliseconds per access, at a default of 100 accesses,
and the number it produces is the real one: fetching a 150 bp slice costs
~160 ms through `SeqIO.index` and single-digit microseconds through pyfaidx
(~1.5 µs for the fetch itself, measured directly).  Both rows return the same
bytes and are digest-gated, so this is a design difference, not a harness
artefact — a FASTA index that stores record offsets cannot serve a slice
without rebuilding the record, and that gap is a Phase 1 target.

Every implementation is additionally held to a per-pass budget
(`--max-pass-seconds`, default 60).  If one pass exceeds it, the repeats are
skipped, the single warmup sample is reported, and the row is flagged
`TRUNCATED` — a degraded measurement is labelled as one instead of being
presented as a median of N.

For the compressed pair (`fastq-gz-index`, `fastq-gz-random`) the runner writes
one derived file and caches it: a **BGZF copy of the reads**, under
`build/bench-cache/`, because `SeqIO.index` refuses a plain `.gz` outright and
the reference row needs a compressed file it will open.  It is written from the
plain corpus file, keyed on the source's name and size, made once and reused; it
is never corpus and never tracked.  It costs about twenty seconds and 174 MB the
first time, and the run asks for it only when those two groups are actually
selected (`--only`).

### How to read the numbers

Every implementation yields the same `(name, sequence, quality_len)` and the
harness does identical minimal work on top of each, so differences come from
parsing and not from the harness.  Times are the median of `--repeat` runs after
one untimed warmup (which also warms the page cache).  `MB/s` is input
consumed per second — for `fastq-gzip`, the *uncompressed* size, so its rows are
on the same scale as `fastq-plain` and the decompression-only rows above them.

**Rows are correctness-gated.**  Before anything is timed, each
implementation's full output is hashed (quality compared as raw Phred values)
and checked against the Biopython reference.  A row that disagrees is printed as
`INVALID` and must never be quoted as a speedup.  If a parser is fast and wrong,
the table says so instead of rewarding it.

The two decompression rows are the exception, and are labelled `baseline`
instead: they emit 1 MiB chunks rather than records, so there is no parse output
to compare against the reference.  Digest-gating them would report a permanent
false `INVALID`, which would in turn make the end-of-run warning meaningless.
They place a floor under the parsing rows and are never quotable as an
alternative to Biopython.

Two measurement traps this harness went out of its way to avoid, because both
would have produced a flattering fake result:

- **Header identity.** `SeqRecord.id` stops at the first space, while
  `FastqGeneralIterator`, the naive parser and a rebuilt pysam header carry the
  whole line.  Comparing them directly made the naive parser look 20% faster
  *and* subtly wrong; the harness now compares full headers everywhere.
- **pyfaidx coordinates.** pyfaidx 0.9 removed `__getitem__`/`__iter__` and made
  `fetch(name, start, end)` **1-based and inclusive**, where 0.8 and earlier
  took 0-based slices.  `PyfaidxAccess` probes which API is present and
  normalises, so "pyfaidx" means the same work in either version.

Results are written to `bench/results/latest.json` (gitignored — they are
machine-specific).  The durable, human-readable ranking is `bench/targets.md`,
produced in step 0.5.

## Profiling (step 0.4)

The runner says how much slower Biopython is; `bench/profile_paths.py` says why,
by naming the functions that hold the time:

```sh
.venv/bin/python bench/profile_paths.py --list
.venv/bin/python bench/profile_paths.py --all --repeat 1 --top 12
.venv/bin/py-spy record --format raw -o /tmp/prof.raw --rate 200 \
    -- .venv/bin/python bench/profile_paths.py --workload fastq:biopython-seqio
.venv/bin/python bench/profile_paths.py --pyspy /tmp/prof.raw \
    --pyspy-mark consume
```

Path names are `<group>:<implementation>` and come straight from the runner's
own factories, so the loop being profiled is the loop that was measured.

It reuses the runner's own workload factories, so the code being profiled is the
code that was measured.  Two views are available and they answer different
questions: `cProfile` attributes time to Python frames and reports the share
that is interpreter work — the part a C core deletes outright — but its seconds
are inflated by instrumentation and must never be compared to `run.py`; `py-spy`
samples the uninstrumented process and so sees time inside C extensions, which
`cProfile` charges to whichever Python frame called them.  Findings are written
up in `bench/profiling.md`.

## Charts

The images in the top-level README are rendered from these same numbers, so a
picture cannot drift away from the table it illustrates — the tables stay the
data source:

```sh
.venv/bin/python bench/plot_whythis.py    # -> misc/why-this-exists.png
.venv/bin/python bench/plot_g1.py         # -> misc/g1-ranking.png
.venv/bin/python bench/plot_profile.py    # -> misc/profile-share.png
```

`plot_g1.py` carries the ranking rows and `plot_profile.py` the profiling shares
verbatim from `targets.md` and `profiling.md`.  matplotlib is a chart
dependency only — `run.py` does not need it (see `requirements.txt`).

Charts are written to `misc/`, which is gitignored except for the files the
README actually references: an unreferenced design asset stays local.  Adding a
chart to the README therefore means adding a `!misc/<name>.png` exception to
`.gitignore` as well, or the image renders as a broken link on GitHub while
looking fine on the machine that produced it.

## Environment notes

The corpus is generated with the standard library only (`hashlib`, `gzip`,
`json`) — no numpy, no Biopython, no network — so it can be regenerated on a
bare Python 3 installation.  A full run takes about 30 s single-threaded on an
i5-13420H.
