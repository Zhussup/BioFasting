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

## Environment notes

The corpus is generated with the standard library only (`hashlib`, `gzip`,
`json`) — no numpy, no Biopython, no network — so it can be regenerated on a
bare Python 3 installation.  A full run takes about 30 s single-threaded on an
i5-13420H.
