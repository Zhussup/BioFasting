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
  top-level README.

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
