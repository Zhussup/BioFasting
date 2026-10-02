#!/usr/bin/env python3
"""BioFasting benchmark runner (PLAN.md step 0.3).

Compares, on the frozen corpus from ``bench/gen_data.py``:

* what people actually write -- ``Bio.SeqIO.parse`` into ``SeqRecord`` objects;
* Biopython's own low-level iterators, which isolate how much of the cost is
  ``SeqRecord`` construction rather than parsing;
* a naive pure-Python line parser -- the floor any C core has to beat;
* the installed C-backed libraries (``pysam``/htslib, ``pyfaidx``);
* for gzipped input, decompression-only baselines, so the table shows how much
  of a compressed pipeline is decompression rather than tokenization.

Method
------
Every implementation yields the same thing -- ``(name, sequence, quality_len)``
-- and the harness does identical, minimal work on top of each: count records,
sum sequence lengths.  Differences in the table therefore come from parsing, not
from the harness.  Quality is reported by length, never re-serialised: doing so
would tax ``SeqRecord`` for a conversion a real user mostly does not pay.

* Timings are the **median of N runs** after one untimed warmup pass.  The
  warmup also pulls the file into the page cache, so numbers measure CPU-bound
  parsing rather than disk.  They are single-machine numbers, not a claim about
  anyone else's machine.
* Correctness is gated: before timing, every implementation's full output is
  hashed and compared against the Biopython reference.  A row whose digest
  differs is reported **INVALID**, and a speedup must never be quoted from it.
  That is this project's principle #3 -- fast-but-wrong is the failure mode we
  exist to remove, so the harness refuses to hand out a number for one.
* ``MB/s`` is bytes of *input* consumed per second.  For the ``.gz`` group that
  is the uncompressed size, so every row in a group is on the same scale and the
  decompression-only baselines are directly comparable with the parsers above
  them.

Usage
-----
    .venv/bin/python bench/run.py                 # full corpus, median of 5
    .venv/bin/python bench/run.py --quick         # test-sized corpus
    .venv/bin/python bench/run.py --only fastq-gzip --repeat 9
    .venv/bin/python bench/run.py --list          # show workloads and exit
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib
import json
import platform
import statistics
import sys
import time
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator

HERE = Path(__file__).parent
DEFAULT_DATA = HERE / "data"
DEFAULT_RESULTS = HERE / "results"

#: Complement table for the naive reverse-complement baseline.  Built with
#: bytes.maketrans so the op itself is one C-level call, exactly as a user would
#: write it -- the comparison is against the Biopython loop, not against SIMD.
_RC_TABLE = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")

#: ASCII quality character -> raw Phred byte.  Gives every implementation the
#: same canonical quality representation for the verification digest, at
#: C level (``bytes.translate``), instead of a per-base Python loop.
_PHRED_TO_BYTE = bytes(max(0, code - 33) for code in range(256))

#: How many reads the in-memory ops workloads load (untimed) and operate on.
_OPS_READS = 100_000


def try_import(module: str):
    """Import an optional backend, returning ``(module_or_None, error)``."""
    try:
        return importlib.import_module(module), None
    except Exception as exc:  # pragma: no cover - environment dependent
        return None, f"{type(exc).__name__}: {exc}"


def _b(value) -> bytes:
    if isinstance(value, bytes):
        return value
    if isinstance(value, str):
        return value.encode()
    return str(value).encode()


def _open_binary(path: Path):
    """Open a corpus file for binary reading, transparently handling ``.gz``."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _open_text(path: Path):
    """Text-mode twin of :func:`_open_binary`.

    Biopython's low-level iterators insist on text handles
    (``FastqGeneralIterator`` raises ``StreamModeError`` on a binary one), which
    is why the corpus readers cannot share a single opener.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def full_header(name: str, comment: str | None) -> str:
    """Rebuild a complete FASTQ header from pysam's split name/comment pair."""
    return f"{name} {comment}" if comment else name


# --------------------------------------------------------------------------- #
# Workload model
# --------------------------------------------------------------------------- #


@dataclass
class Impl:
    """One way of doing a workload.

    ``fast`` is timed and yields ``(name, sequence, quality_len)``; ``verify``
    runs once and yields ``(name, sequence, quality)`` in full.  File readers
    normalise that to ``bytes`` with quality as raw Phred values, so an
    implementation carrying quality as a string and one carrying the numeric
    list are compared in the same units.
    """

    label: str
    fast: Callable[[], Iterator]
    verify: Callable[[], Iterator]
    kind: str = "records"  # "records" | "stream"
    note: str = ""


@dataclass
class Workload:
    group: str
    description: str
    logical_bytes: int
    impls: list[Impl] = field(default_factory=list)
    expected_items: int | None = None
    source: str = ""


# --------------------------------------------------------------------------- #
# FASTQ implementations
# --------------------------------------------------------------------------- #


def make_fastq_impls(path: Path, backends: dict) -> list[Impl]:
    """Readers for a FASTQ file, plain or gzipped (detected from the suffix)."""
    Bio = backends["Bio"]
    compressed = str(path).endswith(".gz")
    impls: list[Impl] = []

    def seqio(full: bool):
        def factory():
            handle = gzip.open(path, "rt") if compressed else str(path)
            close = hasattr(handle, "close")
            try:
                for rec in Bio.SeqIO.parse(handle, "fastq"):
                    quals = rec.letter_annotations["phred_quality"]
                    # description, not id: the low-level parsers hand back the
                    # whole header line, and the digest compares like with like.
                    if full:
                        yield (rec.description.encode(), str(rec.seq).encode(), bytes(quals))
                    else:
                        yield rec.description, rec.seq, len(quals)
            finally:
                if close:
                    handle.close()

        return factory

    impls.append(Impl(
        "biopython-seqio", seqio(False), seqio(True),
        note="Bio.SeqIO.parse -> SeqRecord; what most users write",
    ))

    def lowlevel(full: bool):
        def factory():
            from Bio.SeqIO.QualityIO import FastqGeneralIterator

            with _open_text(path) as fh:
                for title, seq, qual in FastqGeneralIterator(fh):
                    if full:
                        yield (title.encode(), seq.encode(),
                               qual.encode().translate(_PHRED_TO_BYTE))
                    else:
                        yield title, seq, len(qual)

        return factory

    impls.append(Impl(
        "biopython-lowlevel", lowlevel(False), lowlevel(True),
        note="FastqGeneralIterator: no SeqRecord construction",
    ))

    def naive(full: bool):
        def factory():
            with _open_binary(path) as fh:
                while True:
                    header = fh.readline()
                    if not header:
                        return
                    if header in (b"\n", b"\r\n"):
                        continue
                    seq = fh.readline().rstrip(b"\r\n")
                    fh.readline()  # the '+' line
                    qual = fh.readline().rstrip(b"\r\n")
                    title = header[1:].rstrip(b"\r\n")
                    if full:
                        yield title, seq, qual.translate(_PHRED_TO_BYTE)
                    else:
                        yield title, seq, len(qual)

        return factory

    impls.append(Impl(
        "naive-python", naive(False), naive(True),
        note="single-line reads only; the floor a C core must beat",
    ))

    pysam = backends.get("pysam")
    if pysam is not None:

        def via_pysam(full: bool):
            def factory():
                with pysam.FastxFile(str(path)) as fh:
                    for rec in fh:
                        # pysam splits the header; rebuild it for the digest.
                        title = full_header(rec.name, getattr(rec, "comment", None))
                        if full:
                            yield (title.encode(), rec.sequence.encode(),
                                   rec.quality.encode().translate(_PHRED_TO_BYTE))
                        else:
                            yield title, rec.sequence, len(rec.quality)

            return factory

        impls.append(Impl(
            "pysam-htslib", via_pysam(False), via_pysam(True),
            note="pysam.FastxFile (C/htslib)",
        ))
    return impls


def make_decompression_impls(path: Path) -> list[Impl]:
    """Baselines that decompress without parsing anything."""

    def gzip_only():
        def read():
            with gzip.open(path, "rb") as fh:
                while chunk := fh.read(1 << 20):
                    yield chunk

        return read

    def zlib_only():
        def read():
            dec = zlib.decompressobj(31)  # 31 = gzip container
            with open(path, "rb") as fh:
                while chunk := fh.read(1 << 20):
                    out = dec.decompress(chunk)
                    if out:
                        yield out
                tail = dec.flush()
                if tail:
                    yield tail

        return read

    return [
        Impl("gzip-stream-only", gzip_only(), gzip_only(), kind="stream",
             note="decompression with no parsing above it: this file's floor"),
        Impl("zlib-raw-only", zlib_only(), zlib_only(), kind="stream",
             note="same, minus the gzip.GzipFile Python layer"),
    ]


# --------------------------------------------------------------------------- #
# FASTA implementations
# --------------------------------------------------------------------------- #


def _as_sequence(value):
    """Normalise a pyfaidx record/slice to something with ``len()``."""
    seq = getattr(value, "seq", None)
    return seq if seq is not None and not callable(seq) else value


class PyfaidxAccess:
    """Adapter over pyfaidx's two incompatible access APIs.

    pyfaidx <= 0.8 gives ``fa[name]`` a 0-based sliceable record and iterates
    records directly.  pyfaidx 0.9 removed ``__getitem__`` and ``__iter__`` and
    replaced them with ``fetch(name, start, end)``, which is **1-based and
    inclusive**.  The capability probe runs once, outside every timed region.

    This is not incidental detail: a benchmark of "pyfaidx" that silently used
    the wrong coordinates would compare two different amounts of work.
    """

    def __init__(self, fasta) -> None:
        self.fasta = fasta
        self.slicing_api = hasattr(fasta, "__getitem__")

    def names(self) -> list[str]:
        if hasattr(self.fasta, "keys"):
            return list(self.fasta.keys())
        return list(self.fasta.index.keys())

    def slice(self, name: str, start: int, span: int):
        if self.slicing_api:  # pyfaidx <= 0.8: 0-based
            value = self.fasta[name][start : start + span]
        else:  # pyfaidx >= 0.9: 1-based, inclusive
            value = self.fasta.fetch(name, start + 1, start + span)
        return _as_sequence(value)

    def whole(self, name: str):
        if self.slicing_api:
            value = self.fasta[name]
        else:
            value = self.fasta.fetch(name, 1, self.fasta.index[name].rlen)
        return _as_sequence(value)


def make_fasta_impls(path: Path, backends: dict) -> list[Impl]:
    Bio = backends["Bio"]
    impls: list[Impl] = []

    def seqio(full: bool):
        def factory():
            for rec in Bio.SeqIO.parse(str(path), "fasta"):
                if full:
                    yield rec.description.encode(), str(rec.seq).encode(), b""
                else:
                    yield rec.description, rec.seq, ""

        return factory

    impls.append(Impl("biopython-seqio", seqio(False), seqio(True),
                      note="Bio.SeqIO.parse -> SeqRecord"))

    def lowlevel(full: bool):
        def factory():
            from Bio.SeqIO.FastaIO import SimpleFastaParser

            with open(path) as fh:
                for title, seq in SimpleFastaParser(fh):
                    if full:
                        yield title.encode(), seq.encode(), b""
                    else:
                        yield title, seq, ""

        return factory

    impls.append(Impl("biopython-lowlevel", lowlevel(False), lowlevel(True),
                      note="SimpleFastaParser: no SeqRecord construction"))

    def naive(full: bool):
        def factory():
            name = None
            parts: list[bytes] = []
            with open(path, "rb") as fh:
                for line in fh:
                    if line.startswith(b">"):
                        if name is not None:
                            yield name, b"".join(parts), ""
                        name = line[1:].rstrip(b"\r\n")
                        parts = []
                    else:
                        parts.append(line.rstrip(b"\r\n"))
            if name is not None:
                yield name, b"".join(parts), ""

        return factory

    impls.append(Impl("naive-python", naive(False), naive(True),
                      note="line-accumulating reader"))

    pyfaidx = backends.get("pyfaidx")
    if pyfaidx is not None:

        def via_pyfaidx(full: bool):
            def factory():
                fa = pyfaidx.Faidx(str(path))
                try:
                    access = PyfaidxAccess(fa)
                    for name in access.names():
                        value = access.whole(name)
                        if full:
                            yield name.encode(), _b(value), b""
                        else:
                            yield name, value, ""
                finally:
                    fa.close()

            return factory

        impls.append(Impl("pyfaidx-sequential", via_pyfaidx(False), via_pyfaidx(True),
                          note="indexed reader; .fai exists after the warmup pass"))
    return impls


def fasta_lengths(path: Path) -> dict[str, int]:
    """Sequence lengths per record, for building a valid random-access plan."""
    lengths: dict[str, int] = {}
    name = None
    total = 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0].rstrip(b"\r\n").decode()
                total = 0
            else:
                total += len(line.rstrip(b"\r\n"))
    if name is not None:
        lengths[name] = total
    return lengths


def make_fasta_random_impls(path: Path, backends: dict, accesses: int,
                            span: int) -> list[Impl]:
    """Random slice access -- the reason FASTA indexes exist at all.

    Slices are clamped to real record boundaries so that every implementation
    returns the same bytes; otherwise pyfaidx (which clamps) and Biopython
    (which returns empty) would disagree and the rows would be INVALID for a
    reason that says nothing about speed.
    """
    Bio = backends["Bio"]
    impls: list[Impl] = []
    lengths = fasta_lengths(path)
    names = sorted(lengths)
    if not names:
        return impls
    # Fixed, corpus-independent access pattern (no dependency on the corpus rng).
    plan = []
    for i in range(accesses):
        chrom = names[(i * 7919) % len(names)]
        room = max(1, lengths[chrom] - span)
        plan.append((chrom, (i * 104729) % room))

    def bi_index(full: bool):
        # The index is built once and reused: building it is setup, not the
        # access being measured (pyfaidx's equivalent lives in its .fai).
        cache: dict = {}

        def factory():
            idx = cache.get("index")
            if idx is None:
                idx = cache["index"] = Bio.SeqIO.index(str(path), "fasta")
            for chrom, start in plan:
                yield chrom, idx[chrom].seq[start : start + span], ""

        return factory

    impls.append(Impl("biopython-index-random", bi_index(False), bi_index(True),
                      note=f"{accesses} slices of {span} bp; index built once, outside the timing"))

    pyfaidx = backends.get("pyfaidx")
    if pyfaidx is not None:
        cache: dict = {}

        def via_pyfaidx(full: bool):
            def factory():
                access = cache.get("access")
                if access is None:
                    access = cache["access"] = PyfaidxAccess(pyfaidx.Faidx(str(path)))
                for chrom, start in plan:
                    value = access.slice(chrom, start, span)
                    if full:
                        yield chrom.encode(), _b(value), b""
                    else:
                        yield chrom, value, ""

            return factory

        impls.append(Impl("pyfaidx-random", via_pyfaidx(False), via_pyfaidx(True),
                          note=f"{accesses} slices of {span} bp; handle opened once"))
    return impls


# --------------------------------------------------------------------------- #
# In-memory ops (reproducing the draft benchmarks in the README)
# --------------------------------------------------------------------------- #


def make_op_impls(records: list, backends: dict) -> tuple[list[Impl], list[Impl]]:
    """Reverse-complement and GC workloads over already-loaded records.

    Loading happens once outside every timed region and the closures only
    iterate a list, so the measured time is the operation itself.
    """

    def revcomp_bio():
        def factory():
            for rec in records:
                yield rec.reverse_complement().seq

        return factory

    def revcomp_naive():
        def factory():
            for rec in records:
                seq = str(rec.seq)
                yield seq.translate(_RC_TABLE)[::-1]

        return factory

    revcomp = [
        Impl("biopython-revcomp", revcomp_bio(), revcomp_bio(), kind="stream",
             note="SeqRecord.reverse_complement, once per record"),
        Impl("naive-strtranslate", revcomp_naive(), revcomp_naive(), kind="stream",
             note="str.translate + slice"),
    ]

    def gc_bio():
        def factory():
            from Bio.SeqUtils import gc_fraction

            for rec in records:
                yield f"{gc_fraction(rec.seq):.6f}"

        return factory

    def gc_naive():
        def factory():
            for rec in records:
                seq = str(rec.seq)
                yield f"{(seq.count('G') + seq.count('C')) / len(seq):.6f}"

        return factory

    gc = [
        Impl("biopython-gc", gc_bio(), gc_bio(), kind="stream",
             note="SeqUtils.gc_fraction"),
        Impl("naive-strcount", gc_naive(), gc_naive(), kind="stream",
             note="str.count x2; identical semantics on an N-free corpus"),
    ]
    return revcomp, gc


# --------------------------------------------------------------------------- #
# Harness
# --------------------------------------------------------------------------- #


def consume(impl: Impl) -> tuple[int, int]:
    """Iterate once doing minimal, uniform work; returns (items, bases)."""
    items = bases = 0
    if impl.kind == "records":
        for _name, seq, _qlen in impl.fast():
            items += 1
            bases += len(seq)
    else:
        for chunk in impl.fast():
            items += 1
            bases += len(chunk)
    return items, bases


def digest(impl: Impl) -> tuple[int, int, str]:
    """Full-content digest, run once per implementation for correctness."""
    hasher = hashlib.sha256()
    items = bases = 0
    for item in impl.verify():
        if impl.kind == "records":
            name, seq, qual = item
            hasher.update(_b(name))
            hasher.update(b"\x00")
            hasher.update(_b(seq))
            hasher.update(b"\x00")
            hasher.update(_b(qual))
            bases += len(seq)
        else:
            hasher.update(_b(item))
            bases += len(item)
        hasher.update(b"\x01")
        items += 1
    return items, bases, hasher.hexdigest()


def time_impl(impl: Impl, repeats: int, warmups: int = 1) -> list[float]:
    for _ in range(warmups):
        consume(impl)
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        consume(impl)
        times.append(time.perf_counter() - start)
    return sorted(times)


def run_workload(workload: Workload, repeats: int) -> list[dict]:
    print(f"\n== {workload.group}  ({workload.description})")
    reference: str | None = None
    rows: list[dict] = []
    for impl in workload.impls:
        try:
            items, bases, ref_digest = digest(impl)
        except Exception as exc:
            print(f"   {impl.label:26} SKIPPED ({type(exc).__name__}: {exc})")
            rows.append({"group": workload.group, "impl": impl.label, "valid": False,
                         "error": f"{type(exc).__name__}: {exc}", "note": impl.note})
            continue
        if reference is None:
            reference = ref_digest
        valid = ref_digest == reference
        note = impl.note
        if workload.expected_items is not None and items != workload.expected_items:
            valid = False
            note = f"{note} | {items} items != expected {workload.expected_items}"
        times = time_impl(impl, repeats)
        median = statistics.median(times)
        megabytes = workload.logical_bytes / 1e6
        row = {
            "group": workload.group,
            "impl": impl.label,
            "items": items,
            "bases": bases,
            "median_s": round(median, 4),
            "min_s": round(times[0], 4),
            "max_s": round(times[-1], 4),
            "mb_per_s": round(megabytes / median, 1) if median else None,
            "items_per_s": round(items / median) if median else None,
            "valid": valid,
            "note": note,
        }
        rows.append(row)
        print(f"   {impl.label:26} {median:8.3f} s  {row['mb_per_s']:8.1f} MB/s  "
              f"{row['items_per_s']:>13,} /s  {'ok' if valid else 'INVALID'}")
    if reference is not None and not all(r.get("valid") for r in rows):
        print("   !! digest mismatch: never quote a speedup from an INVALID row")
    return rows


# --------------------------------------------------------------------------- #
# Environment and driver
# --------------------------------------------------------------------------- #


def host_info(versions: dict) -> dict:
    flags = ""
    cpu = platform.processor()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
            elif line.startswith("flags"):
                flags = line
                break
    except OSError:
        pass
    return {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cpu": cpu,
        "avx2": "avx2" in flags,
        "avx512": "avx512f" in flags,
        "backends": versions,
    }


def load_manifest(data: Path) -> dict:
    path = data / "MANIFEST.json"
    if not path.exists():
        print(f"error: {path} not found -- run bench/gen_data.py first", file=sys.stderr)
        raise SystemExit(2)
    return json.loads(path.read_text())


def build_workloads(args, backends: dict, manifest: dict) -> list[Workload]:
    sizes = {e["path"]: e["bytes"] for e in manifest["files"]}
    counts = {e["path"]: e.get("records") for e in manifest["files"]}
    data = args.data

    fastq = "fastq/reads_10k.fastq" if args.quick else "fastq/reads_1m.fastq"
    fasta = "fasta/genome_1mb.fasta" if args.quick else "fasta/genome.fasta"
    gz = "fastq/reads_1m.fastq.gz"

    workloads = [Workload(
        "fastq-plain",
        f"{fastq}, {sizes[fastq] / 1e6:.1f} MB, {counts[fastq]:,} records",
        sizes[fastq], make_fastq_impls(data / fastq, backends), counts[fastq],
        source=fastq)]

    if not args.quick and (data / gz).exists():
        workloads.append(Workload(
            "fastq-gzip",
            f"{gz}, {sizes[gz] / 1e6:.1f} MB on disk, {counts[gz]:,} records",
            sizes[gz],  # MB/s on uncompressed size: same scale as fastq-plain
            make_fastq_impls(data / gz, backends) + make_decompression_impls(data / gz),
            counts[gz], source=gz))

    workloads.append(Workload(
        "fasta-scan", f"{fasta}, {sizes[fasta] / 1e6:.1f} MB", sizes[fasta],
        make_fasta_impls(data / fasta, backends), None, source=fasta))

    if not args.quick:
        accesses, span = 20_000, 150
        workloads.append(Workload(
            "fasta-random", f"{fasta}, {accesses:,} slices of {span} bp at fixed positions",
            accesses * span, make_fasta_random_impls(data / fasta, backends, accesses, span),
            accesses, source=fasta))

    Bio = backends["Bio"]
    reads = []
    for rec in Bio.SeqIO.parse(str(data / fastq), "fastq"):
        reads.append(rec)
        if len(reads) >= _OPS_READS:
            break
    ops_bytes = sum(len(rec) for rec in reads)
    revcomp, gc = make_op_impls(reads, backends)
    workloads.append(Workload("ops-revcomp", f"reverse complement over {len(reads):,} reads",
                              ops_bytes, revcomp, len(reads), source="in memory"))
    workloads.append(Workload("ops-gc", f"GC fraction over {len(reads):,} reads",
                              ops_bytes, gc, len(reads), source="in memory"))
    return workloads


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="run.py", description="BioFasting benchmark runner")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--json", type=Path, default=None,
                        help="where to write results (default: bench/results/latest.json)")
    parser.add_argument("--repeat", type=int, default=5, help="timed runs per implementation")
    parser.add_argument("--quick", action="store_true", help="use the test-sized files")
    parser.add_argument("--only", action="append", default=None,
                        help="run only these groups (repeatable)")
    parser.add_argument("--list", action="store_true", help="list workloads and exit")
    args = parser.parse_args(argv)

    backends: dict = {}
    versions: dict = {}
    for name in ("Bio", "pysam", "pyfaidx"):
        module, error = try_import(name)
        if module is not None:
            backends[name] = module
            versions[name] = getattr(module, "__version__", "?")
        else:
            print(f"note: backend {name} unavailable ({error})", file=sys.stderr)
    if "Bio" not in backends:
        print("error: Biopython is required as the reference implementation", file=sys.stderr)
        return 2
    # ``import Bio`` alone does not bind the submodules used below.
    for submodule in ("Bio.SeqIO", "Bio.SeqIO.QualityIO", "Bio.SeqIO.FastaIO", "Bio.SeqUtils"):
        importlib.import_module(submodule)

    manifest = load_manifest(args.data)
    workloads = build_workloads(args, backends, manifest)
    if args.only:
        wanted = set(args.only)
        workloads = [w for w in workloads if w.group in wanted]
        if not workloads:
            print(f"error: no workload matched {sorted(wanted)}", file=sys.stderr)
            return 2

    if args.list:
        for workload in workloads:
            print(f"{workload.group:16} {workload.description}")
            for impl in workload.impls:
                print(f"    {impl.label:26} {impl.note}")
        return 0

    print(f"corpus: {args.data}   repeat: {args.repeat}   "
          f"backends: {', '.join(f'{k} {v}' for k, v in versions.items())}")

    started = time.time()
    results: list[dict] = []
    for workload in workloads:
        results.extend(run_workload(workload, args.repeat))

    out_path = args.json or (DEFAULT_RESULTS / "latest.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        "host": host_info(versions),
        "corpus": {"seed": manifest.get("seed"), "params": manifest.get("params"),
                   "quick": args.quick},
        "repeats": args.repeat,
        "results": results,
    }, indent=2) + "\n")

    print(f"\n{len(results)} rows in {time.time() - started:.1f}s -> {out_path}")
    invalid = [r for r in results if not r.get("valid")]
    if invalid:
        print(f"WARNING: {len(invalid)} row(s) did not validate:")
        for row in invalid:
            print(f"  {row['group']}/{row['impl']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
