#!/usr/bin/env python3
"""BioFasting benchmark corpus generator (PLAN.md step 0.2).

Every benchmark, differential test and fuzz case in this project draws from the
corpus this script produces, so the corpus must be reproducible: same command
in, byte-identical files out, on any machine.

Reproducibility contract
------------------------
* Randomness comes from SHA-256 in counter mode (``Keystream``) keyed by
  ``(seed, domain)``.  SHA-256 has a frozen definition, so the byte stream does
  not depend on the platform, Python version or CPU.
* Every sampling decision is integer arithmetic.  No floating point value ever
  reaches a threshold, so results cannot drift with libm differences.
* Each dataset has its own domain string, so changing one parameter (say the
  read length) does not shift the bytes of an unrelated dataset.
* ``.gz`` members are written with ``mtime=0`` and a fixed compression level,
  so the *content* is identical everywhere.  (The compressed container may
  differ by one or two bytes across zlib versions; the manifest therefore hashes
  the uncompressed content and ``--verify`` decompresses before hashing.)

The gzip concern is not theoretical: real-world FASTQ is compressed, and
decompression -- not tokenization -- is usually the first bottleneck in a
pipeline.  The corpus models realistic compression behaviour (see
``--profile realistic``) rather than uniform random bytes, which would not
compress at all and would flatter a naive parser.

Usage
-----
    python3 bench/gen_data.py             # full corpus -> bench/data/
    python3 bench/gen_data.py --quick     # CI-sized smoke corpus
    python3 bench/gen_data.py --verify    # re-check checksums, generate nothing
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

GENERATOR_VERSION = 1
DEFAULT_SEED = 20261002

#: Illumina TruSeq Read 1 adapter.  Present at the 3' end of a fraction of
#: reads in every real dataset, and the reason trimming tools exist.
ILLUMINA_ADAPTER = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

QUAL_OFFSET_PHRED33 = 33  # Sanger / Illumina 1.8+
QUAL_OFFSET_PHRED64 = 64  # Illumina 1.3-1.7, edge corpus only

#: Read cycles are split into three bands (5' / middle / 3') with different
#: quality profiles.  Real instruments degrade toward the 3' end, and quality
#: trimming tools exist precisely because of it, so a corpus without this
#: decay would let a parser look correct while breaking downstream trimming.
QUAL_BANDS = ((37, 9, 1), (35, 10, 4), (30, 14, 12))

BASE_A, BASE_C, BASE_G, BASE_T = 65, 67, 71, 84


# --------------------------------------------------------------------------- #
# Deterministic randomness
# --------------------------------------------------------------------------- #


class Keystream:
    """A reproducible byte stream: SHA-256 in counter mode over (seed, domain).

    ``take(n)`` returns the next ``n`` bytes.  Two streams built from the same
    ``(seed, domain)`` pair yield the same bytes forever; streams with
    different domains are independent for all practical purposes.
    """

    _CHUNK_BLOCKS = 512  # 16 KiB generated per hashing burst

    def __init__(self, seed: int, domain: bytes) -> None:
        seed_bytes = seed.to_bytes(8, "big", signed=False)
        self._prefix = hashlib.sha256(
            b"biofasting-corpus-v1\x00" + seed_bytes + b"\x00" + domain
        ).digest()
        self._counter = 0
        self._buf = b""
        self._pos = 0

    def take(self, n: int) -> bytes:
        if n < 0:
            raise ValueError("keystream length must be non-negative")
        if self._pos + n > len(self._buf):
            tail = self._buf[self._pos:]
            missing = n - len(tail)
            blocks = []
            remaining = (missing + 31) // 32
            while remaining > 0:
                burst = min(remaining, self._CHUNK_BLOCKS)
                for _ in range(burst):
                    blocks.append(
                        hashlib.sha256(
                            self._prefix + self._counter.to_bytes(8, "big")
                        ).digest()
                    )
                    self._counter += 1
                remaining -= burst
            self._buf = tail + b"".join(blocks)
            self._pos = 0
        out = self._buf[self._pos : self._pos + n]
        self._pos += n
        return out


def base_table(gc_percent: int) -> bytes:
    """256-byte table mapping a keystream byte to a base at the given GC%."""
    at_each = (100 - gc_percent) * 256 // 200
    cg_each = gc_percent * 256 // 200
    table = bytearray(256)
    pos = 0
    for count, base in (
        (at_each, BASE_A),
        (cg_each, BASE_C),
        (cg_each, BASE_G),
    ):
        if pos + count > 256:
            raise ValueError(f"GC={gc_percent}% does not fit a 256-byte table")
        table[pos : pos + count] = bytes([base]) * count
        pos += count
    table[pos:] = bytes([BASE_T]) * (256 - pos)
    return bytes(table)


def quality_table(weights: dict[int, int], offset: int) -> bytes:
    """256-byte table mapping a keystream byte to a quality character."""
    pairs = sorted(weights.items())
    starts: list[tuple[int, int]] = []
    cumulative = 0
    for score, weight in pairs:
        starts.append((cumulative, score))
        cumulative += weight
    total = cumulative
    table = bytearray(256)
    bucket = 0
    for i in range(256):
        target = (i * total + total // 2) // 256
        while bucket + 1 < len(starts) and target >= starts[bucket + 1][0]:
            bucket += 1
        table[i] = offset + starts[bucket][1]
    return bytes(table)


def uniform_quality_table(low: int, high: int, offset: int) -> bytes:
    """Flat distribution over Phred ``low..high`` (used for fuzz corpora)."""
    return quality_table({score: 1 for score in range(low, high + 1)}, offset)


def ramp_quality_table(peak: int, spread: int, floor: int, offset: int) -> bytes:
    """Triangular quality profile with an error floor, in integer arithmetic.

    ``peak`` is the modal Phred score, ``spread`` the half-width of the
    triangle, and ``floor`` a uniform weight added to every score -- the
    sequencing-error tail that real instruments always show.
    """
    weights = {}
    for score in range(2, 42):
        height = max(1, 100 - (abs(score - peak) * 100) // spread)
        weights[score] = height + floor
    return quality_table(weights, offset)


# --------------------------------------------------------------------------- #
# Output plumbing
# --------------------------------------------------------------------------- #


class Sink:
    """Writes one artifact while tracking its uncompressed size and SHA-256."""

    def __init__(
        self,
        path: Path,
        *,
        relative: str,
        gzip_out: bool = False,
        level: int = 6,
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.relative = relative
        self.gzip_out = gzip_out
        self._raw = path.open("wb")
        self._fh = (
            gzip.GzipFile(fileobj=self._raw, mode="wb", compresslevel=level, mtime=0)
            if gzip_out
            else self._raw
        )
        self._hash = hashlib.sha256()
        self.size = 0

    def write(self, data: bytes) -> None:
        if not data:
            return
        self._fh.write(data)
        self._hash.update(data)
        self.size += len(data)

    def close(self, *, kind: str, records: int | None = None, notes: str = "") -> dict:
        self._fh.close()
        if self.gzip_out:
            self._raw.close()
        entry = {
            "path": self.relative,
            "kind": kind,
            "bytes": self.size,
            "sha256": self._hash.hexdigest(),
            "gzip": self.gzip_out,
            "stored_bytes": self.path.stat().st_size,
        }
        if records is not None:
            entry["records"] = records
        if notes:
            entry["notes"] = notes
        return entry


def gzip_copy(source: Path, target: Path, level: int = 6) -> None:
    """Re-encode an existing artifact as gzip with deterministic headers."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as src, target.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=level, mtime=0) as out:
            while chunk := src.read(1 << 20):
                out.write(chunk)


# --------------------------------------------------------------------------- #
# FASTQ
# --------------------------------------------------------------------------- #


def _illumina_header(prefix: str, index: int) -> bytes:
    """Build a realistic Casava 1.8 header with arithmetic (no RNG) fields."""
    tile = 1101 + index % 4
    x = 1000 + (index * 7919) % 29000
    y = 1000 + (index * 104729) % 29000
    return (
        f"{prefix}.{index} HWI-ST1:1:FCBAR:1:{tile}:{x}:{y} 1:N:0:ATCACG".encode()
    )


def write_fastq(
    sink: Sink,
    *,
    reads: int,
    length: int,
    seed: int,
    domain: bytes,
    gc_percent: int,
    profile: str,
    adapter_rate: int,
    dup_rate: int,
    id_prefix: str,
    small_sink: Sink | None = None,
    small_reads: int = 0,
    history_size: int = 4096,
    flush_every: int = 16384,
) -> None:
    """Write ``reads`` FASTQ records; optionally tee the first ``small_reads``."""
    ks_seq = Keystream(seed, domain + b"/seq")
    ks_qual = Keystream(seed, domain + b"/qual")
    ks_meta = Keystream(seed, domain + b"/meta")
    bases = base_table(gc_percent)
    if profile == "realistic":
        bands = [
            ramp_quality_table(peak, spread, floor, QUAL_OFFSET_PHRED33)
            for peak, spread, floor in QUAL_BANDS
        ]
    else:
        flat = uniform_quality_table(0, 41, QUAL_OFFSET_PHRED33)
        bands = [flat, flat, flat]
    band_cuts = (length // 3, 2 * length // 3)
    dup_cut = dup_rate * 256 // 100
    adapter_cut = adapter_rate * 256 // 100
    adapter_max = min(len(ILLUMINA_ADAPTER), max(length - 40, 0))

    history: list[tuple[bytes, bytes]] = []
    history_next = 0
    buf: list[bytes] = []
    small_buf: list[bytes] = []

    for index in range(reads):
        seq = ks_seq.take(length).translate(bases)
        raw_qual = ks_qual.take(length)
        c0, c1 = band_cuts
        qual = (
            raw_qual[:c0].translate(bands[0])
            + raw_qual[c0:c1].translate(bands[1])
            + raw_qual[c1:].translate(bands[2])
        )
        meta = ks_meta.take(4)

        if history and meta[0] < dup_cut:
            # PCR duplicate: an exact copy of a recently seen molecule.
            seq, qual = history[((meta[1] << 8) | meta[2]) % len(history)]
        else:
            if adapter_max > 12 and meta[3] < adapter_cut:
                # 3' adapter read-through: the adapter starts mid-read.
                plen = 12 + meta[2] % (adapter_max - 12 + 1)
                seq = seq[: length - plen] + ILLUMINA_ADAPTER[:plen]
            if len(history) < history_size:
                history.append((seq, qual))
            else:
                # Ring buffer: overwriting in place keeps this O(1) per read,
                # where deleting the oldest entry would be O(len(history)).
                history[history_next] = (seq, qual)
                history_next = (history_next + 1) % history_size

        record = (
            b"@"
            + _illumina_header(id_prefix, index)
            + b"\n"
            + seq
            + b"\n+\n"
            + qual
            + b"\n"
        )
        buf.append(record)
        if small_sink is not None and index < small_reads:
            small_buf.append(record)
        if len(buf) >= flush_every:
            sink.write(b"".join(buf))
            buf.clear()
        if len(small_buf) >= flush_every:
            small_sink.write(b"".join(small_buf))
            small_buf.clear()

    sink.write(b"".join(buf))
    if small_sink is not None:
        small_sink.write(b"".join(small_buf))


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #


def _write_wrapped(
    sink: Sink, data: bytes, carry: bytearray, width: int
) -> None:
    """Append ``data`` to the stream, wrapping lines at ``width`` bases."""
    carry += data
    whole = len(carry) - (len(carry) % width)
    if whole == 0:
        return
    block = bytes(carry[:whole])
    del carry[:whole]
    sink.write(
        b"\n".join(block[i : i + width] for i in range(0, whole, width)) + b"\n"
    )


def write_fasta_records(
    sink: Sink,
    *,
    records: list[tuple[str, int]],
    seed: int,
    gc_percent: int,
    width: int,
    n_runs_per_record: int = 2,
    n_run_length: int = 1000,
    soft_mask_record: str | None = None,
    chunk: int = 1 << 20,
) -> None:
    """Stream a multi-record genome.  N-runs and soft-masking are overlaid at
    fixed positions so the result stays reproducible and exercises the paths
    real parsers must special-case."""
    bases = base_table(gc_percent)
    for name, size in records:
        sink.write(b">" + name.encode() + b"\n")
        ks = Keystream(seed, b"fasta/" + name.encode())
        # Deterministic N-runs: telomere-ish blocks at fixed fractional offsets.
        runs = []
        for k in range(n_runs_per_record):
            start = size * (2 * k + 1) // (2 * n_runs_per_record + 1)
            runs.append((start, min(start + n_run_length, size)))
        mask_start = mask_end = -1
        if soft_mask_record == name:
            mask_start = size // 4
            mask_end = mask_start + size // 20
        carry = bytearray()
        offset = 0
        while offset < size:
            take = min(chunk, size - offset)
            block = bytearray(ks.take(take).translate(bases))
            for start, end in runs:
                lo, hi = max(start - offset, 0), min(end - offset, take)
                if lo < hi:
                    block[lo:hi] = b"N" * (hi - lo)
            lo, hi = max(mask_start - offset, 0), min(mask_end - offset, take)
            if lo < hi:
                block[lo:hi] = bytes(block[lo:hi]).lower()
            _write_wrapped(sink, bytes(block), carry, width)
            offset += take
        if carry:
            sink.write(bytes(carry) + b"\n")


# --------------------------------------------------------------------------- #
# Edge and malformed corpora (hand-written literals, deliberately tiny)
# --------------------------------------------------------------------------- #

EDGE_FASTQ: dict[str, tuple[bytes, str]] = {
    "edge/multiline.fastq": (
        b"@read1 wrapped sequence and quality\n"
        b"ACGTACGTACGTACGTACGT\n"
        b"ACGTACGTAC\n"
        b"+\n"
        b"IIIIIIIIIIIIIIIIIIII\n"
        b"IIIIIIIIII\n"
        b"@read2 single line for contrast\n"
        b"TTTTGGGGCCCCAAAA\n"
        b"+\n"
        b"JJJJHHHHFFFFDDDD\n",
        "FASTQ allows sequence and quality to wrap across lines.",
    ),
    "edge/quality_contains_at_and_plus.fastq": (
        b"@read1 quality is all '@' -- a parser that scans for '@' to find the "
        b"next header corrupts this record\n"
        b"ACGTACGTACGT\n"
        b"+\n"
        b"@@@@@@@@@@@@\n"
        b"@read2 quality starts with '+'\n"
        b"ACGTACGTACGT\n"
        b"+\n"
        b"++IIIIIIIIII\n",
        "The classic '@'-in-quality trap; also '+' as a quality character.",
    ),
    "edge/empty_read.fastq": (
        b"@empty\n\n+\n\n@nonempty\nACGT\n+\nIIII\n",
        "Zero-length sequence and quality are legal and must round-trip.",
    ),
    "edge/phred64.fastq": (
        b"@legacy1 Illumina 1.3-1.7 Phred+64 qualities\n"
        b"ACGTACGTACGTACGTACGT\n"
        b"+\n"
        b"hhhhhhhhhhhhhhhhhhaa\n"
        b"@legacy2\n"
        b"TTTTGGGGCCCCAAAA\n"
        b"+\n"
        b"hhhhhhaaaaaaBBBB\n",
        "Legacy Phred+64; naive readers silently mis-score these. Biopython 1.88 "
        "reads the first base as Q71 instead of the true Q40.",
    ),
    "edge/crlf.fastq": (
        b"@read1 CRLF line endings\r\nACGTACGT\r\n+\r\nIIIIIIII\r\n"
        b"@read2\r\nTTTTGGGG\r\n+\r\nJJJJHHHH\r\n",
        "Windows line endings: '\\r' must not leak into sequence or quality.",
    ),
    "edge/no_trailing_newline.fastq": (
        b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTGGGG\n+\nJJJJHHHH",
        "Last record ends at EOF with no terminating newline.",
    ),
    "edge/empty_file.fastq": (
        b"",
        "Empty input; the parser must yield zero records, not crash.",
    ),
    "edge/full_quality_range.fastq": (
        b"@q_lowest\nACGT\n+\n!!!!\n"
        b"@q_highest\nACGT\n+\nJJJJ\n"
        b"@q_diverse\nACGTACGTACGT\n+\n!\"#$%&'()*+,\n"
        b"@q_to_tilde\nACGT\n+\n~~~~\n",
        "Q0 through Q93 present; exercises the full printable quality range.",
    ),
    "edge/blank_lines.fastq": (
        b"@read1\nACGTACGT\n+\nIIIIIIII\n\n\n@read2\nTTTTGGGG\n+\nJJJJHHHH\n",
        "Stray blank lines between records. Biopython's FastqGeneralIterator "
        "tolerates them (2 records); whether we should is a compatibility "
        "decision, not an obvious malformedness call.",
    ),
    "edge/long_read.fastq": (
        b"@long_read length=4096\n"
        + b"ACGT" * 1024
        + b"\n+\n"
        + b"I" * 4096
        + b"\n",
        "One nanopore-scale read; catches fixed-size buffer assumptions.",
    ),
}

MALFORMED_FASTQ: dict[str, tuple[bytes, str]] = {
    # Biopython 1.88 behaviour is recorded in each note: these become the
    # expected-error cases for the differential tests in step 1.5.
    "malformed/truncated_mid_quality.fastq": (
        b"@read1\nACGTACGT\n+\nIIII\n@read2\nTTTT\n+\nJJJ",
        "File ends mid-record. Biopython 1.88 raises a length mismatch for "
        "read1 (8 vs 18) rather than an EOF error.",
    ),
    "malformed/missing_plus_line.fastq": (
        b"@read1\nACGTACGT\nIIIIIIII\n",
        "Third line is not '+'. Biopython 1.88: 'End of file without quality "
        "information.'",
    ),
    "malformed/length_mismatch.fastq": (
        b"@read1\nACGTACGT\n+\nIIII\n",
        "Sequence and quality lengths disagree; silent acceptance corrupts data. "
        "Biopython 1.88 rejects it naming the read.",
    ),
    "malformed/header_only.fastq": (
        b"@read1\n",
        "Header with no sequence, quality or '+' line. Biopython 1.88: "
        "'Unexpected end of file.'",
    ),
}

EDGE_FASTA: dict[str, tuple[bytes, str]] = {
    "edge/fasta/blank_lines.fasta": (
        b">rec1\nACGTACGT\n\nACGTACGT\n\n>rec2\nTTTTGGGG\n",
        "Blank lines inside and between records.",
    ),
    "edge/fasta/empty_record.fasta": (
        b">empty\n>nonempty\nACGT\n",
        "A record with an empty sequence.",
    ),
    "edge/fasta/multiline_widths.fasta": (
        b">wrapped\n"
        + b"ACGTACGTAC\n" * 3
        + b"ACGT\n"
        + b">flat\nACGTACGTACGTACGT\n",
        "Varying wrap widths in the same file.",
    ),
    "edge/fasta/soft_masked.fasta": (
        b">masked\nACGTacgtACGTACGTacgtACGT\n",
        "Lower-case soft-masking must survive parsing.",
    ),
    "edge/fasta/comment_lines.fasta": (
        b"; a comment line before the first record\n"
        b">rec1 description with spaces\n"
        b"; an in-record comment\n"
        b"ACGTACGT\n",
        "';' comment lines, legacy FASTA extension. Biopython's strict 'fasta' "
        "parser rejects a leading comment; 'fasta-pearson' accepts it.",
    ),
    "edge/fasta/no_trailing_newline.fasta": (
        b">rec1\nACGTACGT",
        "Last record ends at EOF with no terminating newline.",
    ),
    "edge/fasta/crlf.fasta": (
        b">rec1\r\nACGTACGT\r\n>rec2\r\nTTTTGGGG\r\n",
        "Windows line endings.",
    ),
    "edge/fasta/single_long_line.fasta": (
        b">one_line\n" + b"ACGT" * 2500 + b"\n",
        "10 kb on one line; catches fixed-width line buffers.",
    ),
}


def write_literal_corpus(sink_factory, entries: dict[str, tuple[bytes, str]], kind: str) -> list[dict]:
    manifest = []
    for name, (payload, note) in sorted(entries.items()):
        sink = sink_factory(name)
        sink.write(payload)
        records = payload.count(b"\n@") + (1 if payload.startswith(b"@") else 0)
        if kind == "fasta":
            records = payload.count(b"\n>") + (1 if payload.startswith(b">") else 0)
        manifest.append(sink.close(kind=kind, records=records, notes=note))
    return manifest


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def build(outdir: Path, args: argparse.Namespace) -> dict:
    started = time.time()
    # The generator owns the corpus tree.  Wipe the managed subdirectories first
    # so a dataset that was renamed or dropped in a later version cannot linger
    # on disk as a stale file the manifest no longer describes.
    for sub in ("fastq", "fasta", "edge", "malformed"):
        shutil.rmtree(outdir / sub, ignore_errors=True)
    entries: list[dict] = []

    def make(rel: str, *, gz: bool = False) -> Sink:
        return Sink(outdir / rel, relative=rel, gzip_out=gz)

    def collect(sink: Sink, **kw) -> None:
        entries.append(sink.close(**kw))

    # --- FASTQ: the flagship corpus ------------------------------------
    big = make("fastq/reads_1m.fastq")
    small = make("fastq/reads_10k.fastq")
    write_fastq(
        big,
        reads=args.reads,
        length=args.read_length,
        seed=args.seed,
        domain=b"fastq/main",
        gc_percent=args.gc,
        profile=args.profile,
        adapter_rate=args.adapter_rate,
        dup_rate=args.dup_rate,
        id_prefix="SRR000001",
        small_sink=small,
        small_reads=args.small_reads,
    )
    collect(
        big,
        kind="fastq",
        records=args.reads,
        notes=(
            f"{args.reads} reads x {args.read_length} bp, {args.profile} profile, "
            f"GC={args.gc}%, adapter={args.adapter_rate}%, duplicates={args.dup_rate}%"
        ),
    )
    collect(small, kind="fastq", records=args.small_reads,
            notes="Exact prefix of the large file (shared keystream); test-sized.")
    gzip_copy(outdir / "fastq/reads_1m.fastq", outdir / "fastq/reads_1m.fastq.gz")
    stored = (outdir / "fastq/reads_1m.fastq.gz").stat().st_size
    ratio = entries[0]["bytes"] / stored if stored else 0.0
    entries.append(
        {
            "path": "fastq/reads_1m.fastq.gz",
            "kind": "fastq",
            "bytes": entries[0]["bytes"],
            "sha256": entries[0]["sha256"],
            "gzip": True,
            "records": args.reads,
            "stored_bytes": stored,
            "notes": (
                f"gzip of reads_1m.fastq, mtime=0, ratio {ratio:.2f}x. Production "
                "FASTQ compresses to about 3.5-4x; independent draws here compress "
                "less, so this payload is a harder-than-real decompression case "
                "(the bias understates, never overstates, a decompressor's win)."
            ),
        }
    )

    # --- FASTA: genome-sized -------------------------------------------
    genome = make("fasta/genome.fasta")
    total = args.genome_size
    weights = [40, 30, 20, 8, 2]
    records = []
    allocated = 0
    for i, weight in enumerate(weights):
        size = total * weight // 100 if i < len(weights) - 1 else total - allocated
        allocated += size
        records.append((f"chr{i + 1}", size))
    write_fasta_records(
        genome,
        records=records,
        seed=args.seed,
        gc_percent=args.gc,
        width=60,
        n_runs_per_record=2,
        n_run_length=1000,
        soft_mask_record="chr3",
    )
    collect(genome, kind="fasta", records=len(records),
            notes=f"{total} bp across {len(records)} records, 60-col wrap, N-runs, soft-masked chr3")

    small_genome = make("fasta/genome_1mb.fasta")
    write_fasta_records(
        small_genome,
        records=[("chrS", args.small_genome_size)],
        seed=args.seed,
        gc_percent=args.gc,
        width=60,
        n_runs_per_record=1,
        n_run_length=1000,
    )
    collect(small_genome, kind="fasta", records=1,
            notes="Single-record genome for indexed-access tests.")

    # --- Edge and malformed corpora ------------------------------------
    entries += write_literal_corpus(make, EDGE_FASTQ, "fastq")
    entries += write_literal_corpus(make, MALFORMED_FASTQ, "fastq")
    entries += write_literal_corpus(make, EDGE_FASTA, "fasta")

    entries.sort(key=lambda item: item["path"])
    manifest = {
        "generator": "bench/gen_data.py",
        "generator_version": GENERATOR_VERSION,
        "seed": args.seed,
        "params": {
            "reads": args.reads,
            "read_length": args.read_length,
            "small_reads": args.small_reads,
            "genome_size": total,
            "small_genome_size": args.small_genome_size,
            "gc_percent": args.gc,
            "profile": args.profile,
            "adapter_rate": args.adapter_rate,
            "dup_rate": args.dup_rate,
        },
        "files": entries,
        "totals": {
            "files": len(entries),
            # Content size counts each artifact once; the .gz entry stores the
            # same bytes as its plain counterpart, so it must not be added twice.
            "content_bytes": sum(e["bytes"] for e in entries if not e.get("gzip")),
            "stored_bytes": sum(e["stored_bytes"] for e in entries),
        },
        "elapsed_seconds": round(time.time() - started, 1),
    }
    return manifest


def verify(outdir: Path, seed: int | None) -> int:
    manifest_path = outdir / "MANIFEST.json"
    if not manifest_path.exists():
        print(f"error: {manifest_path} not found; run the generator first", file=sys.stderr)
        return 2
    manifest = json.loads(manifest_path.read_text())
    if seed is not None and manifest.get("seed") != seed:
        print(
            f"error: manifest seed {manifest.get('seed')} != requested --seed {seed}",
            file=sys.stderr,
        )
        return 2
    bad = 0
    for entry in manifest["files"]:
        path = outdir / entry["path"]
        if not path.exists():
            print(f"MISSING  {entry['path']}")
            bad += 1
            continue
        digest = hashlib.sha256()
        size = 0
        opener = gzip.open if entry.get("gzip") else open
        with opener(path, "rb") as fh:
            while chunk := fh.read(1 << 20):
                digest.update(chunk)
                size += len(chunk)
        if digest.hexdigest() != entry["sha256"] or size != entry["bytes"]:
            print(f"MISMATCH {entry['path']} ({size} bytes)")
            bad += 1
    if bad:
        print(f"\n{bad} of {len(manifest['files'])} artifacts failed verification")
        return 1
    print(f"OK: all {len(manifest['files'])} artifacts match MANIFEST.json")
    return 0


def parse_size(text: str) -> int:
    """Parse a size such as ``100M`` or ``3.1G`` into a base count."""
    suffixes = {"k": 1000, "K": 1000, "m": 1000**2, "M": 1000**2, "g": 1000**3, "G": 1000**3}
    if text[-1] in suffixes:
        return int(float(text[:-1]) * suffixes[text[-1]])
    return int(text)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="gen_data.py",
        description="Generate the seeded BioFasting benchmark corpus.",
    )
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "data",
                        help="output directory (default: bench/data)")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="seed for every keystream (default: %(default)s)")
    parser.add_argument("--reads", type=int, default=1_000_000,
                        help="reads in the large FASTQ (default: %(default)s)")
    parser.add_argument("--small-reads", type=int, default=10_000,
                        help="reads in the test-sized FASTQ prefix (default: %(default)s)")
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--genome-size", type=parse_size, default=100_000_000,
                        help="total genome bases, e.g. 100M or 3.1G (default: 100M)")
    parser.add_argument("--small-genome-size", type=parse_size, default=1_000_000)
    parser.add_argument("--gc", type=int, default=45, help="GC content percent")
    parser.add_argument("--profile", choices=("realistic", "uniform"), default="realistic",
                        help="realistic models adapters/duplicates/illumina qualities")
    parser.add_argument("--adapter-rate", type=int, default=15,
                        help="percent of reads carrying 3' adapter read-through")
    parser.add_argument("--dup-rate", type=int, default=5,
                        help="percent of reads that are PCR duplicates")
    parser.add_argument("--quick", action="store_true",
                        help="CI-sized smoke corpus (implies --out data-quick)")
    parser.add_argument("--verify", action="store_true",
                        help="verify existing artifacts against MANIFEST.json and exit")
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing corpus without asking")
    args = parser.parse_args(argv)

    if args.quick:
        args.reads = 20_000
        args.small_reads = 2_000
        args.genome_size = 2_000_000
        args.small_genome_size = 200_000
        if args.out == Path(__file__).parent / "data":
            args.out = Path(__file__).parent / "data-quick"

    if args.verify:
        return verify(args.out, args.seed)

    outdir: Path = args.out
    if (outdir / "MANIFEST.json").exists() and not args.force:
        print(f"error: {outdir}/MANIFEST.json already exists; pass --force to overwrite",
              file=sys.stderr)
        return 2

    outdir.mkdir(parents=True, exist_ok=True)
    manifest = build(outdir, args)
    # Wall-clock is printed, never recorded: MANIFEST.json must stay byte-stable.
    elapsed = manifest.pop("elapsed_seconds")
    (outdir / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=False) + "\n"
    )
    totals = manifest["totals"]
    print(
        f"wrote {totals['files']} artifacts, "
        f"{totals['content_bytes'] / 1e6:.1f} MB of content, "
        f"{totals['stored_bytes'] / 1e6:.1f} MB on disk "
        f"in {elapsed}s -> {outdir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
