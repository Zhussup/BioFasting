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

#: Where a workload may leave a derived file it needs but the corpus does not
#: ship.  Under `build/`, which is disposable and never tracked: anything here
#: must be reproducible from the corpus and cheap enough to lose.
BENCH_CACHE = HERE.parent / "build" / "bench-cache"


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
    # A baseline does not produce records at all -- decompression with no parser
    # above it, say.  It has no parse output to compare against the reference,
    # so it is never digest-gated and never quotable as a speedup; it exists to
    # place a floor under the rows that do parse.
    baseline: bool = False


@dataclass
class Workload:
    group: str
    description: str
    logical_bytes: int
    impls: list[Impl] = field(default_factory=list)
    expected_items: int | None = None
    source: str = ""
    # Set when one "item" is one random access, where latency per access is the
    # meaningful number and throughput in MB/s hides the real cost.
    per_access: bool = False


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

    biofasting = backends.get("biofasting")
    if biofasting is not None:

        def via_biofasting(full: bool):
            def factory():
                # One pass over a memory mapping, one record at a time; for the
                # .gz workload the inflated buffer takes the mapping's place.
                for title, seq, qual in biofasting.open_fastq(path):
                    if full:
                        yield title, seq, qual.translate(_PHRED_TO_BYTE)
                    else:
                        yield title, seq, len(qual)

            return factory

        impls.append(Impl(
            "biofasting", via_biofasting(False), via_biofasting(True),
            note="this project: mmap + a scanner over the bytes; libdeflate for .gz",
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
        Impl("gzip-stream-only", gzip_only(), gzip_only(), kind="stream", baseline=True,
             note="decompression with no parsing above it: this file's floor"),
        Impl("zlib-raw-only", zlib_only(), zlib_only(), kind="stream", baseline=True,
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

    biofasting = backends.get("biofasting")
    if biofasting is not None:

        def via_biofasting(full: bool):
            def factory():
                index = biofasting.open_fasta(path)
                for key in index:
                    # The title is the whole header line here, as `SeqIO`'s
                    # description is, so the digest compares like with like.
                    if full:
                        yield index.title(key).encode(), index[key], b""
                    else:
                        yield index.title(key), index[key], ""

            return factory

        impls.append(Impl("biofasting", via_biofasting(False), via_biofasting(True),
                          note="this project: one scan builds an offset index; "
                               "records are copied out by offset"))

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

    biofasting = backends.get("biofasting")
    if biofasting is not None:
        cache: dict = {}

        def via_biofasting(full: bool):
            def factory():
                index = cache.get("index")
                if index is None:
                    index = cache["index"] = biofasting.open_fasta(path)
                for chrom, start in plan:
                    # A gather out of the mapping: no SeqRecord, no record-wide
                    # copy, which is the difference this row is here to show.
                    value = index.sequence_slice(chrom, start, start + span)
                    if full:
                        yield chrom.encode(), value, b""
                    else:
                        yield chrom, value, ""

            return factory

        impls.append(Impl("biofasting-random", via_biofasting(False), via_biofasting(True),
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
# FASTQ random access by record name
# --------------------------------------------------------------------------- #


def fastq_record_names(path: Path) -> list[str]:
    """The record names, in file order, read once at workload setup.

    Every fourth line is a header, because that is the FASTQ format and not a
    guess about its contents -- a quality line really can begin with ``@``
    (Phred+33 offset 31 is the character ``?``, and 64 is ``@``), which is the
    trap that made `gen_data.py` overcount a corpus file once already.  The same
    four-line assumption is what `SeqIO.index`'s own `FastqRandomAccess` makes.

    The name is everything up to the first space or tab, which is the key the
    reference indexes a record under.
    """
    names = []
    with path.open("rb") as handle:
        for index, line in enumerate(handle):
            if index % 4 == 0:
                rest = line[1:] if line[:1] == b"@" else line
                names.append(rest.split(b" ", 1)[0].split(b"\t", 1)[0]
                             .rstrip(b"\r\n").decode("ascii"))
    return names


def make_fastq_index_impls(path: Path, backends: dict, names: list[str],
                           accesses: int) -> tuple[list[Impl], list[Impl]]:
    """Random access into a FASTQ by record name, and the build that enables it.

    Two different things are timed here and they are deliberately separate rows,
    because for a 369 MB FASTQ they are different sizes and a reader of the
    table has to be able to tell which one they are looking at:

    * `fastq-index` builds the index at all, which for `SeqIO.index` means
      walking the entire file in Python and holding a million keys;
    * `fastq-random` touches it, with the build outside the timing -- the same
      convention `fasta-random` uses, and there it is the right one because the
      build is setup.  Here it is not setup, which is exactly why it gets its
      own row rather than hiding inside this one.

    Returns ``(build_impls, random_impls)``.  `names` arrives already read,
    because reading it means walking the file once and doing that per workload
    would be two passes for one answer.
    """
    Bio = backends["Bio"]
    if not names:
        return [], []
    # Fixed, corpus-independent access pattern, so two runs touch the same
    # records and a digest difference means a real difference.
    plan = [names[(i * 7919) % len(names)] for i in range(accesses)]

    # Both rows report (name, sequence, quality length), and the third field is
    # the interesting one to get right.  It is not a padded slot to make two
    # digests agree: it is the same quantity on both sides, because the FASTQ
    # grammar refuses a record whose quality is not exactly as long as its
    # sequence, so the length of one is the length of the other.  Biopython's
    # index cannot hand over a single field of a record -- it parses the whole
    # record to get any of it -- so its third field costs it nothing; ours reads
    # the length out of the index without touching the file.  Reporting the
    # quality itself instead would charge the reference for rebuilding bytes it
    # never exposes, which would flatter this table rather than inform it.
    #
    # The two rows are therefore "get record X's sequence and its size, two
    # different ways", which is the comparison this workload exists to make.

    build: list[Impl] = []
    random: list[Impl] = []

    def bi_build(full: bool):
        def factory():
            index = Bio.SeqIO.index(str(path), "fastq")
            for name in index:
                yield name, b"", ""
            index.close()

        return factory

    build.append(Impl("biopython-index-build", bi_build(False), bi_build(True),
                      note="SeqIO.index(..., 'fastq') over the whole file, then its keys"))

    def bi_random(full: bool):
        cache: dict = {}

        def factory():
            index = cache.get("index")
            if index is None:
                index = cache["index"] = Bio.SeqIO.index(str(path), "fastq")
            for name in plan:
                record = index[name]
                yield name, record.seq, len(record.seq)

        return factory

    random.append(Impl("biopython-index-random", bi_random(False), bi_random(True),
                       note=f"{accesses} records by name; index built once, outside the timing"))

    biofasting = backends.get("biofasting")
    if biofasting is not None and hasattr(biofasting, "open_fastq_index"):
        def via_biofasting_build(full: bool):
            def factory():
                index = biofasting.open_fastq_index(path)
                for name in index:
                    yield name, b"", ""

            return factory

        build.append(Impl("biofasting-index-build", via_biofasting_build(False),
                          via_biofasting_build(True),
                          note="one scan; no record objects, no quality decoding"))

        cache: dict = {}

        def via_biofasting_random(full: bool):
            def factory():
                index = cache.get("index")
                if index is None:
                    index = cache["index"] = biofasting.open_fastq_index(path)
                for name in plan:
                    # Name, bases and quality length, straight out of the
                    # mapping -- no SeqRecord, no Seq, no per-field decode.
                    yield name, index.sequence(name), index.quality_length(name)
            return factory

        random.append(Impl("biofasting-index-random", via_biofasting_random(False),
                           via_biofasting_random(True),
                           note=f"{accesses} records by name; index built once, outside the timing"))
    return build, random


# --------------------------------------------------------------------------- #
# FASTQ random access by record name, when the file is compressed
# --------------------------------------------------------------------------- #


def bgzf_companion(path: Path) -> Path:
    """A BGZF copy of `path`'s content, made once and cached under `build/`.

    `SeqIO.index` will not index a plain gzip file.  It raises

        ValueError: Gzipped files are not suitable for indexing, please use
        BGZF (blocked gzip format) instead.

    so the only way to put a reference row beside a row of ours on a
    *compressed* FASTQ is to hand the reference a file it will agree to open.
    This makes that file from the corpus file's own bytes, so both rows read the
    same content, and the harness's digest gate proves it rather than assuming
    it.

    It is deliberately not corpus.  The corpus is what the parsers are diffed
    against, every file in it is carried through the differential tests, and a
    third copy of the same 369 MB of reads -- 174 MB on disk, and about twenty
    seconds to write -- would be paid by every test run for the benefit of two
    benchmark rows.  The cache is keyed on the source's name and size, so
    regenerating the corpus invalidates it instead of quietly measuring the
    wrong bytes.
    """
    path = Path(path)
    cache = BENCH_CACHE / f"{path.name}.{path.stat().st_size}.bgz"
    if cache.exists():
        return cache
    from Bio import bgzf

    cache.parent.mkdir(parents=True, exist_ok=True)
    partial = cache.with_suffix(cache.suffix + ".part")
    with open(path, "rb") as source, bgzf.BgzfWriter(str(partial), "wb") as out:
        while chunk := source.read(1 << 20):
            out.write(chunk)
    partial.replace(cache)
    return cache


def make_fastq_gz_index_impls(gz_path: Path, index_path: Path, backends: dict,
                              names: list[str],
                              accesses: int) -> tuple[list[Impl], list[Impl]]:
    """The `fastq-index`/`fastq-random` pair again, on a compressed file.

    Same two rows, same access plan, same digest shape as
    :func:`make_fastq_index_impls`, so the four numbers are comparable across
    the two tables.  The difference is what each side is handed: ours indexes
    the `.gz` the corpus ships -- what a real pipeline actually has -- and the
    reference indexes `index_path`, a BGZF copy of the same content, because
    plain gzip is the one thing it refuses to index.

    That asymmetry is the measurement, not a defect in it.  Indexing a
    compressed FASTQ with Biopython means re-encoding the file as BGZF first,
    and the price of doing so is reported separately in `bench/targets.md`,
    where a table has room to say it.  Here the reference row's time is the
    index build and nothing else, and its `note` says which file it opened, so
    nobody has to remember.

    Both rows still report ``(name, sequence, quality length)`` -- see the
    argument in :func:`make_fastq_index_impls` for why that third field is the
    same quantity on both sides.
    """
    Bio = backends["Bio"]
    if not names:
        return [], []
    plan = [names[(i * 7919) % len(names)] for i in range(accesses)]

    build: list[Impl] = []
    random: list[Impl] = []

    def bi_build(full: bool):
        def factory():
            index = Bio.SeqIO.index(str(index_path), "fastq")
            for name in index:
                yield name, b"", ""
            index.close()

        return factory

    build.append(Impl(
        "biopython-index-build", bi_build(False), bi_build(True),
        note=(f"SeqIO.index(..., 'fastq') over {Path(index_path).name}, a BGZF "
              "copy of the same content: SeqIO.index refuses the plain .gz")))

    def bi_random(full: bool):
        cache: dict = {}

        def factory():
            index = cache.get("index")
            if index is None:
                index = cache["index"] = Bio.SeqIO.index(str(index_path), "fastq")
            for name in plan:
                record = index[name]
                yield name, record.seq, len(record.seq)

        return factory

    random.append(Impl(
        "biopython-index-random", bi_random(False), bi_random(True),
        note=(f"{accesses} records by name from {Path(index_path).name}; "
              "index built once, outside the timing")))

    biofasting = backends.get("biofasting")
    if biofasting is not None and hasattr(biofasting, "open_fastq_index"):
        def via_biofasting_build(full: bool):
            def factory():
                index = biofasting.open_fastq_index(gz_path)
                for name in index:
                    yield name, b"", ""

            return factory

        build.append(Impl(
            "biofasting-index-build", via_biofasting_build(False),
            via_biofasting_build(True),
            note=f"one inflate of {Path(gz_path).name}, then one scan of the inflated bytes"))

        cache: dict = {}

        def via_biofasting_random(full: bool):
            def factory():
                index = cache.get("index")
                if index is None:
                    index = cache["index"] = biofasting.open_fastq_index(gz_path)
                for name in plan:
                    yield name, index.sequence(name), index.quality_length(name)

            return factory

        random.append(Impl(
            "biofasting-index-random", via_biofasting_random(False),
            via_biofasting_random(True),
            note=(f"{accesses} records by name; index built once from the "
                  "compressed file, outside the timing")))

    return build, random


# --------------------------------------------------------------------------- #
# In-memory ops (reproducing the draft benchmarks in the README)
# --------------------------------------------------------------------------- #


def make_op_impls(records: list, backends: dict) -> tuple[list[Impl], list[Impl], list[Impl]]:
    """Reverse-complement, GC and translate workloads over already-loaded records.

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

    biofasting = backends.get("biofasting")
    if biofasting is not None:
        # The operations are already differential-tested against Biopython on
        # every byte value; what this row adds is the cost of calling them from
        # Python, which is what a user actually pays.  The strings round-trip
        # through `str()` so that both rows are fed the same object the
        # reference rows get.
        def revcomp_ours():
            def factory():
                for rec in records:
                    yield biofasting.reverse_complement(str(rec.seq))

            return factory

        def gc_ours():
            def factory():
                for rec in records:
                    yield f"{biofasting.gc_fraction(str(rec.seq)):.6f}"

            return factory

        level = biofasting.seqops_level()
        revcomp.append(Impl("biofasting-revcomp", revcomp_ours(), revcomp_ours(),
                            kind="stream", note=f"this project: dispatched kernel, level={level}"))
        gc.append(Impl("biofasting-gc", gc_ours(), gc_ours(), kind="stream",
                       note=f"this project: dispatched kernel, level={level}"))

    translate = make_translate_impls(records, backends)
    return revcomp, gc, translate


def make_translate_impls(records: list, backends: dict) -> list[Impl]:
    """Protein translation over already-loaded records.

    The corpus reads are 150 bp -- a whole number of codons -- so no row here
    has to decide what a trailing partial codon means.  That question is settled
    by the differential tests, not by a timing; this workload measures one
    reading frame on sequences where the answer is uncontroversial.
    """
    from Bio.Data import CodonTable

    standard = CodonTable.unambiguous_dna_by_id[1]
    #: codon -> amino acid, with stops folded in and everything else "X", which
    #: is what the reference answers for a codon the standard table does not
    #: name.  Built once, outside every timed region.
    codon_table = dict(standard.forward_table)
    codon_table.update({codon: "*" for codon in standard.stop_codons})

    def translate_bio():
        def factory():
            for rec in records:
                yield str(rec.seq.translate())

        return factory

    def translate_naive():
        def factory():
            get = codon_table.get

            for rec in records:
                seq = str(rec.seq)
                yield "".join([get(seq[i:i + 3], "X")
                               for i in range(0, len(seq) - 2, 3)])

        return factory

    impls = [
        Impl("biopython-translate", translate_bio(), translate_bio(), kind="stream",
             note="Seq.translate, once per record"),
        Impl("naive-dict-per-codon", translate_naive(), translate_naive(), kind="stream",
             note="dict lookup per codon; the floor a C core has to beat"),
    ]
    if backends.get("numpy") is not None:
        vectorised = _numpy_translate(records)
        impls.append(Impl("numpy-64-entry", vectorised, vectorised, kind="stream",
                          note="3-base index into a 64-entry amino-acid array"))
    biofasting = backends.get("biofasting")
    if biofasting is not None and hasattr(biofasting, "translate"):
        def translate_ours():
            def factory():
                for rec in records:
                    yield biofasting.translate(str(rec.seq))

            return factory

        impls.append(Impl("biofasting-translate", translate_ours(), translate_ours(),
                          kind="stream",
                          note=f"this project: dispatched kernel, level={biofasting.seqops_level()}"))
    return impls


def _numpy_translate(records: list):
    """A vectorised translator: bases to 0..3, three bases to one of 64 entries."""
    import numpy as np
    from Bio.Data import CodonTable

    standard = CodonTable.unambiguous_dna_by_id[1]
    index = np.full(256, 255, dtype=np.uint8)
    for value, letter in enumerate("TCAG"):
        index[ord(letter)] = value
        index[ord(letter.lower())] = value
    amino = np.full(64, ord("X"), dtype=np.uint8)
    for codon, letter in standard.forward_table.items():
        amino[_codon_index(index, codon)] = ord(letter)
    for codon in standard.stop_codons:
        amino[_codon_index(index, codon)] = ord("*")

    def factory():
        frombuffer = np.frombuffer
        for rec in records:
            raw = frombuffer(str(rec.seq).encode(), dtype=np.uint8)
            whole = (len(raw) // 3) * 3
            triplets = index[raw[:whole]].reshape(-1, 3)
            keys = (triplets[:, 0].astype(np.uint16) * 16
                    + triplets[:, 1] * 4 + triplets[:, 2])
            yield bytes(amino[keys])

    return factory


def _codon_index(index, codon: str) -> int:
    """``"ATG"`` -> 48: two bits per base, most significant first."""
    a, b, c = (index[ord(letter)] for letter in codon)
    return int(a) * 16 + int(b) * 4 + int(c)


# --------------------------------------------------------------------------- #
# Bio.SeqUtils: the measurement workloads
# --------------------------------------------------------------------------- #

#: Fewer reads than the operation workloads above, because the reference's rows
#: here are three orders of magnitude more expensive per read: `GC123` is 330 us
#: against `reverse_complement`'s fraction of a microsecond, so 100,000 reads
#: would be half a minute per pass before the repeats are even counted.
SEQUTILS_READS = 20_000


def make_sequtils_impls(records: list, backends: dict, genome: bytes | None
                        ) -> list[tuple[str, str, list[Impl], int]]:
    """The `Bio.SeqUtils` measurement workloads: composition, mass, checksums.

    Returns ``(group, description, impls, items)`` quadruples, one per workload,
    in the shape :func:`build_workloads` wants.  The item count is carried out
    rather than recomputed there, because the mass row runs over a subset of
    the reads and a count that did not know that would mark a correct row
    INVALID.  Every row streams over the same
    pre-extracted string list, so what is timed is the measurement and not the
    extraction, and both sides are handed the same object -- a `str`, which is
    what `SeqUtils` is defined on.

    One row here is unlike the others.  `sequtils-gc123-genome` is a single call
    on a whole chromosome, where the reference takes ninety seconds; it is here
    because that number is the reason the kernel exists, and the harness will
    flag it as truncated rather than spend five passes on it.
    """
    from Bio.Seq import translate as _translate
    from Bio.SeqUtils import (
        GC123 as _GC123,
        GC_skew as _GC_skew,
        CodonAdaptationIndex as _CAI,
        molecular_weight as _mw,
        seq1 as _seq1,
        seq3 as _seq3,
        six_frame_translations as _six,
    )
    from Bio.SeqUtils.CheckSum import crc64 as _crc64
    from Bio.SeqUtils.CheckSum import gcg as _gcg

    biofasting = backends.get("biofasting")
    texts = [str(rec.seq) for rec in records]

    def stream(fn, items):
        def factory():
            for text in items:
                yield fn(text)

        return factory

    def row(name, description, ours, reference, naive=None, naive_note="",
            items=None, reference_name=None):
        """One workload from a reference, a baseline, and ours."""
        used = texts if items is None else items
        impls = [Impl(f"biopython-{name}", stream(reference, used),
                      stream(reference, used), kind="stream",
                      note=f"{reference_name or 'Bio.SeqUtils.' + name}, "
                           f"once per read")]
        if naive is not None:
            impls.append(Impl(f"naive-{name}", stream(naive, used),
                              stream(naive, used), kind="stream", note=naive_note))
        if biofasting is not None and ours is not None:
            level = biofasting.seqops_level()
            impls.append(Impl(f"biofasting-{name}", stream(ours, used),
                              stream(ours, used), kind="stream",
                              note=f"this project: dispatched kernel, level={level}"))
        return (f"ops-{name}", description, impls, len(used))

    workloads = []

    def gc123_strided(text):
        gc_all = n_all = 0
        frames = []
        for frame in range(3):
            chunk = text[frame::3]
            gc = chunk.count("G") + chunk.count("C")
            gc += chunk.count("g") + chunk.count("c")
            atgc = gc + chunk.count("A") + chunk.count("T")
            atgc += chunk.count("a") + chunk.count("t")
            frames.append(gc * 100.0 / atgc if atgc else 0)
            gc_all += gc
            n_all += atgc
        return (100.0 * gc_all / n_all, *frames)

    workloads.append(row(
        "gc123", f"GC123 over {len(texts):,} reads",
        biofasting.GC123 if biofasting else None, _GC123, gc123_strided,
        "three strided slices and str.count; the whole-sequence row is the sum "
        "of the three frames, because the trailing bases are padded away",
        reference_name="Bio.SeqUtils.GC123"))

    def gc_skew_naive(text):
        values = []
        for i in range(0, len(text), 30):
            window = text[i:i + 30]
            g = window.count("G") + window.count("g")
            c = window.count("C") + window.count("c")
            values.append((g - c) / (g + c) if g + c else 0.0)
        return values

    workloads.append(row(
        "gc-skew", f"GC_skew(window=30) over {len(texts):,} reads",
        (lambda t: biofasting.GC_skew(t, 30)) if biofasting else None,
        lambda t: _GC_skew(t, 30), gc_skew_naive, "windowed str.count",
        reference_name="Bio.SeqUtils.GC_skew"))

    #   Only the reads `molecular_weight` will accept: it is defined on the
    #   unambiguous alphabet, and a single N anywhere in a read is a ValueError
    #   on both sides.  The count is in the description so that the row is not
    #   silently about a smaller corpus than the one above it.
    unambiguous = [text for text in texts if set(text) <= set("ACGT")]
    workloads.append(row(
        "molecular-weight",
        f"molecular_weight over {len(unambiguous):,} unambiguous reads",
        biofasting.molecular_weight if biofasting else None, _mw, None,
        items=unambiguous, reference_name="Bio.SeqUtils.molecular_weight"))

    workloads.append(row(
        "crc64", f"crc64 over {len(texts):,} reads",
        biofasting.crc64 if biofasting else None, _crc64, None,
        reference_name="Bio.SeqUtils.CheckSum.crc64"))

    workloads.append(row(
        "gcg", f"gcg over {len(texts):,} reads",
        biofasting.gcg if biofasting else None, _gcg, None,
        reference_name="Bio.SeqUtils.CheckSum.gcg"))

    def roundtrip(fn_translate, fn_seq3, fn_seq1):
        return lambda t: fn_seq1(fn_seq3(fn_translate(t)))

    workloads.append(row(
        "protein-roundtrip",
        f"translate, then seq3, then seq1 over {len(texts):,} reads",
        roundtrip(biofasting.translate, biofasting.seq3, biofasting.seq1)
        if biofasting else None,
        roundtrip(_translate, _seq3, _seq1),
        reference_name="Bio.Seq.translate + SeqUtils.seq3 + SeqUtils.seq1"))

    #   The index is built once, outside the timed region, from the same reads
    #   it is then used to score: the row is `calculate`, not the table.
    index = (biofasting.CodonAdaptationIndex(texts).calculate
             if biofasting is not None else None)
    workloads.append(row(
        "cai-calculate",
        f"CodonAdaptationIndex.calculate over {len(texts):,} reads",
        index, _CAI(texts).calculate,
        reference_name="Bio.SeqUtils.CodonAdaptationIndex.calculate"))

    workloads.append(row(
        "six-frame", f"six_frame_translations over {len(texts):,} reads",
        biofasting.six_frame_translations if biofasting else None, _six, None,
        reference_name="Bio.SeqUtils.six_frame_translations"))

    if genome:
        text = genome.decode("ascii")
        impls = [Impl("biopython-gc123", stream(_GC123, [text]),
                      stream(_GC123, [text]), kind="stream",
                      note="Bio.SeqUtils.GC123, once, on a whole contig")]
        if biofasting is not None:
            impls.append(Impl("biofasting-gc123", stream(biofasting.GC123, [text]),
                              stream(biofasting.GC123, [text]), kind="stream",
                              note="this project: dispatched kernel"))
        workloads.append(("sequtils-gc123-genome",
                          f"GC123 on one {len(text):,} bp contig", impls, 1))

    return workloads


# --------------------------------------------------------------------------- #
# Bio.SeqUtils.ProtParam: the protein corpus and its workloads
# --------------------------------------------------------------------------- #

#: Amino-acid frequencies in UniProtKB/Swiss-Prot, per 10,000 residues, from the
#: release statistics the reference's own documentation quotes.  They sum to
#: 9,983 rather than 10,000 because the source rounds each to two decimals, which
#: is harmless: the table below is built by integer accumulation and normalises.
PROTEIN_FREQUENCIES = {
    "A": 825, "R": 553, "N": 406, "D": 545, "C": 138,
    "Q": 393, "E": 675, "G": 707, "H": 227, "I": 591,
    "L": 966, "K": 584, "M": 242, "F": 386, "P": 470,
    "S": 656, "T": 534, "W": 108, "Y": 292, "V": 687,
}

#: Sizes and counts of the generated protein corpus.  `bench/data/` is DNA and has
#: no protein in it, so these are generated rather than sampled -- see the note in
#: `bench/targets.md`.  Weighted so that every size is represented and the 10,000
#: residue proteins do not set the whole corpus's price: 400 proteins, 192,800
#: residues, of which a quarter is the five long ones.
PROTEIN_CORPUS = ((150, 200), (300, 120), (1024, 75), (10000, 5))

PROTEIN_SEED = 20261003


def protein_table() -> bytes:
    """256-byte table mapping a keystream byte to an amino-acid letter.

    Integer arithmetic only, for the same reason :class:`Keystream` is SHA-256:
    a corpus that is sampled with a float threshold can differ in the last
    residue between two machines, and a benchmark row that is not reproducible is
    not a measurement.
    """
    total = sum(PROTEIN_FREQUENCIES.values())
    table = bytearray(256)
    position = 0
    for letter, weight in sorted(PROTEIN_FREQUENCIES.items()):
        count = weight * 256 // total
        table[position:position + count] = bytes([ord(letter)]) * count
        position += count
    #   The remainder, if any, goes to the last letter rather than being spread.
    #   It is at most nineteen bytes and which letter they are is a fixed
    #   consequence of the table above, so the corpus is still reproducible.
    table[position:] = bytes([ord(sorted(PROTEIN_FREQUENCIES)[-1])]) * (256 - position)
    return bytes(table)


def make_protein_corpus() -> list[tuple[str, int]]:
    """``(sequence, residues)`` for the generated protein corpus, in size order."""
    #   `bench/run.py` is normally run as a script, which puts `bench/` on
    #   `sys.path` and makes this a plain import.  It is added explicitly anyway,
    #   because the corpus is part of the reproducibility contract and a runner
    #   invoked any other way must still find the one generator that defines it.
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    from gen_data import Keystream

    stream = Keystream(PROTEIN_SEED, b"protein")
    table = protein_table()
    corpus = []
    for size, count in PROTEIN_CORPUS:
        for _ in range(count):
            raw = stream.take(size)
            corpus.append(("".join(map(chr, (table[b] for b in raw))), size))
    return corpus


#: The group names the protein corpus feeds.  Named here so that the corpus can
#: be left ungenerated on a run that asked for none of them.
PROTPARAM_GROUPS = (
    "protparam-count-amino-acids",
    "protparam-amino-acids-percent",
    "protparam-aromaticity",
    "protparam-secondary-structure-fraction",
    "protparam-molar-extinction-coefficient",
    "protparam-flexibility",
    "protparam-instability-index",
    "protparam-gravy",
    "protparam-protein-scale",
    "protparam-isoelectric-point",
    "protparam-charge-at-ph",
)


def make_protparam_impls(corpus: list[tuple[str, int]], backends: dict
                         ) -> list[tuple[str, str, list[Impl], int]]:
    """The `Bio.SeqUtils.ProtParam` workloads, over the generated protein corpus.

    Every row is a fresh object per call, because that is what a caller pays and
    what the ranking pass in `bench/targets.md` measured the reference at: on a
    used object `count_amino_acids` is a cache hit, which is 0.19 us against the
    twenty `str.count` calls behind it.

    The naive rows are the best rewrite measured in that pass, kept here so that
    the delivered number can be read against the ceiling it was chosen over.  They
    are digest-gated like everything else -- and the gate has to hold for them,
    because `flexibility` differs in the last bit on every protein if its mirror
    pairs are folded into one nine-element weight vector, and that fold is exactly
    the rewrite that looks natural.
    """
    from Bio.SeqUtils import ProtParamData as _PD
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as _PA

    biofasting = backends.get("biofasting")
    texts = [text for text, _ in corpus]
    residues = sum(size for _, size in corpus)
    level = biofasting.seqops_level() if biofasting is not None else None

    def call(analysis, method, *args):
        """One stream that builds a fresh object and calls one method on it.

        `analysis` is the `ProteinAnalysis` class itself and not the module it
        lives in, because there are two of them and the row's own class is what
        has to be timed.
        """

        def factory():
            for text in texts:
                yield getattr(analysis(text), method)(*args)

        return factory

    def attr(analysis, name):
        """The same, for a `cached_property` -- `amino_acids_percent` is one.

        Which is not a detail a benchmark may average over: calling it is a
        `TypeError` in the reference and reading it is the call, so the row has
        to read it.  On a fresh object it is not a cache hit either -- it is the
        first and only read, and it is `count_amino_acids` underneath.
        """

        def factory():
            for text in texts:
                yield getattr(analysis(text), name)

        return factory

    def row(name, method, description, ours_impl, naive=None, args=(),
            attr_row=False):
        """One group: the reference, its best rewrite, and this project's kernel.

        ``name`` is a group name and uses hyphens everywhere else in this file;
        ``method`` is the reference's own spelling of the call, so that a row
        says `ProteinAnalysis.instability_index` and not `instability-index`.
        """
        reference_impl = (
            attr(_PA, method) if attr_row else call(_PA, method, *args)
        )
        impls = [Impl(f"biopython-{name}", reference_impl, reference_impl,
                      kind="stream",
                      note=f"Bio.SeqUtils.ProtParam.ProteinAnalysis.{method}, "
                           f"on a fresh object per call")]
        if naive is not None:
            factory, naive_note = naive
            impls.append(Impl(f"naive-{name}", factory, factory, kind="stream",
                              note=naive_note))
        if ours_impl is not None:
            impls.append(Impl(f"biofasting-{name}", ours_impl, ours_impl,
                              kind="stream",
                              note=f"this project: dispatched kernel, level={level}"))
        return (f"protparam-{name}", description, impls, len(texts))

    def flex_naive():
        """The reference's window, with the mirror pairs kept in their grouping.

        Folding the four pairs into one symmetric nine-element weight vector is
        the obvious rewrite and the wrong one: same arithmetic on paper, a
        different double, and the digest gate says so on every protein.
        """
        flex = _PD.Flex
        weights = (0.25, 0.4375, 0.625, 0.8125)

        def factory():
            for text in texts:
                scores = []
                for i in range(len(text) - 9):
                    score = 0.0
                    for j in range(4):
                        score += (flex[text[i + j]]
                                  + flex[text[i + 8 - j]]) * weights[j]
                    score += flex[text[i + 5]]
                    scores.append(score / 5.25)
                yield scores

        return factory

    def instability_naive():
        """The 400 weights flattened into one dict, and a concatenated key.

        `flat[a + b]` and `index[a][b]` hold the same float objects, so the
        accumulation is the reference's own order over the reference's own
        numbers and the digest gate is expected to pass -- which is what makes
        the ratio a ceiling rather than a wish.  The reference slices
        `self.sequence[i : i + 2]`; the concatenation is what was measured.
        """
        flat = {}
        for first, row in _PD.DIWV.items():
            for second, value in row.items():
                flat[first + second] = value

        def factory():
            for text in texts:
                score = 0.0
                previous = text[0]
                for aa in text[1:]:
                    score += flat[previous + aa]
                    previous = aa
                yield (10.0 / len(text)) * score

        return factory

    def gravy_naive():
        """`sum(map(scale.__getitem__, text))`, which is the 2.07x rewrite.

        A generator expression and `map` over a bound `__getitem__` compute the
        same Neumaier sum over the same objects; only the per-item cost differs,
        and that is the whole of the difference between the reference and its
        best rewrite.
        """
        scale = _PD.kd
        get = scale.__getitem__

        def factory():
            for text in texts:
                yield sum(map(get, text)) / len(text)

        return factory

    def scale_naive():
        """The reference's loop for the default call, `window=9, edge=1.0`.

        At `edge=1.0` `_weight_list` returns `[1.0, 1.0, 1.0, 1.0]`, so each
        mirror pair contributes exactly `front + back` and the divisor is
        `sum(weights) * 2 + 1 == 9.0`.  The multiplication is kept in anyway:
        it is exact, and it is where a rewrite that "improves" the weights would
        go wrong.
        """
        scale = _PD.kd
        weights = (1.0, 1.0, 1.0, 1.0)

        def factory():
            for text in texts:
                scores = []
                for i in range(len(text) - 8):
                    sub = text[i:i + 9]
                    score = 0.0
                    for j in range(4):
                        front = scale[sub[j]]
                        back = scale[sub[8 - j]]
                        score += weights[j] * front + weights[j] * back
                    score += scale[sub[4]]
                    scores.append(score / 9.0)
                yield scores

        return factory

    def protein_scale_ours():
        """Ours, on the scale this package ships -- so the kernel, not the loop.

        `_scale_blob` keys on identity, and `protparam_data.kd` is the very dict
        `SCALES_BY_NAME["kd"]` holds; a caller passing a *copy* of the scale gets
        the reference's own loop instead, which is the fallback and not this row.
        """
        scale = biofasting.protparam_data.kd
        analysis = biofasting.ProteinAnalysis

        def factory():
            for text in texts:
                yield analysis(text).protein_scale(scale, 9)

        return factory

    ours = biofasting.ProteinAnalysis if biofasting is not None else None
    ours_flex = call(ours, "flexibility") if ours is not None else None
    ours_inst = call(ours, "instability_index") if ours is not None else None
    ours_gravy = call(ours, "gravy") if ours is not None else None
    ours_pi = call(ours, "isoelectric_point") if ours is not None else None
    ours_charge = call(ours, "charge_at_pH", 7.0) if ours is not None else None
    ours_scale = protein_scale_ours() if biofasting is not None else None

    #   The six that begin by counting the residues.  They are here as their own
    #   rows rather than folded into `count-amino-acids` because they are separate
    #   `perf` rows in the triage, and because what they add on top of the count
    #   is the thing a reader wants to see next to it: for most of them, very
    #   little.
    ours_count = call(ours, "count_amino_acids") if ours is not None else None
    ours_percent = attr(ours, "amino_acids_percent") if ours is not None else None
    ours_arom = call(ours, "aromaticity") if ours is not None else None
    ours_ssf = call(ours, "secondary_structure_fraction") if ours is not None else None
    ours_mec = (call(ours, "molar_extinction_coefficient")
                if ours is not None else None)

    return [
        row("count-amino-acids", "count_amino_acids",
            f"ProteinAnalysis.count_amino_acids on a fresh object, {residues:,} "
            f"residues",
            ours_count),
        row("amino-acids-percent", "amino_acids_percent",
            f"ProteinAnalysis.amino_acids_percent on a fresh object, {residues:,} "
            f"residues",
            ours_percent, attr_row=True),
        row("aromaticity", "aromaticity",
            f"ProteinAnalysis.aromaticity on a fresh object, {residues:,} residues",
            ours_arom),
        row("secondary-structure-fraction", "secondary_structure_fraction",
            f"ProteinAnalysis.secondary_structure_fraction on a fresh object, "
            f"{residues:,} residues",
            ours_ssf),
        row("molar-extinction-coefficient", "molar_extinction_coefficient",
            f"ProteinAnalysis.molar_extinction_coefficient on a fresh object, "
            f"{residues:,} residues",
            ours_mec),
        row("flexibility", "flexibility",
            f"ProteinAnalysis.flexibility on a fresh object, {len(texts)} "
            f"proteins, {residues:,} residues",
            ours_flex,
            (flex_naive(),
             "the reference's window with its dict and weight vector hoisted to "
             "locals, mirror pairs kept in the reference's grouping")),
        row("instability-index", "instability_index",
            f"ProteinAnalysis.instability_index on a fresh object, {residues:,} "
            f"residues",
            ours_inst,
            (instability_naive(),
             "the 400 weights in one flat dict, keyed by a two-character "
             "concatenation instead of the reference's slice")),
        row("gravy", "gravy",
            f"ProteinAnalysis.gravy(KyteDoolitle) on a fresh object, "
            f"{residues:,} residues",
            ours_gravy,
            (gravy_naive(),
             "sum(map(scale.__getitem__, text)): the same sum over the same "
             "objects, without the per-item generator frame")),
        row("protein-scale", "protein_scale",
            f"ProteinAnalysis.protein_scale(kd, window=9, edge=1.0) on a fresh "
            f"object, {residues:,} residues",
            ours_scale,
            (scale_naive(),
             "the reference's loop with the two dict lookups hoisted out of the "
             "expression they are read in"),
            args=(_PD.kd, 9)),
        row("isoelectric-point", "isoelectric_point",
            f"ProteinAnalysis.isoelectric_point on a fresh object, {len(texts)} "
            f"proteins, {residues:,} residues",
            ours_pi),
        row("charge-at-ph", "charge_at_pH",
            f"ProteinAnalysis.charge_at_pH(7.0) on a fresh object, {len(texts)} "
            f"proteins",
            ours_charge, args=(7.0,)),
    ]


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
            #   A stream row usually yields text, but not always: the
            #   `Bio.SeqUtils` measurement rows yield floats and lists of them,
            #   and `len` on a float is a TypeError.  Nothing is done with the
            #   byte count on those rows -- what is timed is the call -- so a
            #   scalar simply contributes nothing to it.
            bases += len(chunk) if isinstance(chunk, (str, bytes, list, tuple)) else 0
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
            #   Coerced first, because a stream row may yield a number or a
            #   list and not only text -- `GC123` yields four floats and
            #   `GC_skew` a list of them.  ``len(item)`` on a float is a
            #   TypeError, and a row that raises inside the digest is reported
            #   as skipped rather than measured, which is how this was found.
            payload = _b(item)
            hasher.update(payload)
            bases += len(payload)
        hasher.update(b"\x01")
        items += 1
    return items, bases, hasher.hexdigest()


def time_impl(impl: Impl, repeats: int, max_pass: float | None = None
              ) -> tuple[list[float], bool, tuple[int, int]]:
    """Time ``impl``; return the sorted samples, a truncation flag, and counts.

    The warmup pass is timed too, but only so that a pathologically slow
    implementation can be caught before it burns the whole run: if one pass
    exceeds ``max_pass`` seconds the warmup sample is reported alone and the
    repeats are skipped.  That is a degraded measurement -- the row is flagged
    as truncated rather than quietly presented as a median of N.
    """
    start = time.perf_counter()
    items, bases = consume(impl)
    warmup = time.perf_counter() - start
    if max_pass is not None and warmup > max_pass:
        return [warmup], True, (items, bases)
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        consume(impl)
        times.append(time.perf_counter() - start)
    return sorted(times), False, (items, bases)


def run_workload(workload: Workload, repeats: int,
                 max_pass: float | None = None) -> list[dict]:
    print(f"\n== {workload.group}  ({workload.description})", flush=True)
    reference: str | None = None
    rows: list[dict] = []
    for impl in workload.impls:
        note = impl.note
        if impl.baseline:
            # Structurally different output (chunks, not records): there is
            # nothing to compare against the parse reference.  Report the time
            # and say plainly that this row is a floor, not a candidate.
            times, truncated, (items, bases) = time_impl(impl, repeats, max_pass=max_pass)
            valid = None
            note = f"{note} | baseline: no parse output to verify"
        else:
            try:
                items, bases, ref_digest = digest(impl)
            except Exception as exc:
                print(f"   {impl.label:26} SKIPPED ({type(exc).__name__}: {exc})", flush=True)
                rows.append({"group": workload.group, "impl": impl.label, "valid": False,
                             "error": f"{type(exc).__name__}: {exc}", "note": impl.note})
                continue
            if reference is None:
                reference = ref_digest
            valid = ref_digest == reference
            if workload.expected_items is not None and items != workload.expected_items:
                valid = False
                note = f"{note} | {items} items != expected {workload.expected_items}"
            times, truncated, _ = time_impl(impl, repeats, max_pass=max_pass)
        if truncated:
            note = f"{note} | TRUNCATED: one pass exceeded {max_pass:.0f}s budget"
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
            "samples": len(times),
            "truncated": truncated,
            "mb_per_s": round(megabytes / median, 1) if median else None,
            "items_per_s": round(items / median) if median else None,
            "valid": valid,
            "note": note,
        }
        if workload.per_access:
            # Microseconds, not milliseconds: the spread on this workload is
            # four orders of magnitude, and millisecond precision renders the
            # fast rows as `0.000` -- which reads as "this did nothing" rather
            # than "this answered in 170 ns".  A row that looks like a no-op is
            # the flattering fake result this harness exists to avoid, so the
            # unit follows the magnitude.
            row["us_per_access"] = round(median * 1e6 / items, 3)
        rows.append(row)
        if workload.per_access:
            microseconds = median * 1e6 / items
            extra = (f"  {microseconds:9.3f} us/access" if microseconds < 1000
                     else f"  {microseconds / 1000:9.3f} ms/access")
        else:
            extra = ""
        status = "baseline" if valid is None else ("ok" if valid else "INVALID")
        print(f"   {impl.label:26} {median:8.3f} s  {row['mb_per_s']:8.1f} MB/s  "
              f"{row['items_per_s']:>13,} /s{extra}  {status}", flush=True)
    if reference is not None and any(r.get("valid") is False for r in rows):
        print("   !! digest mismatch: never quote a speedup from an INVALID row", flush=True)
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
    # A time is only meaningful against the binary that produced it: the
    # compiler, the standard, the optimisation flags, and the SIMD rung the
    # dispatcher picked are all part of the number, and a machine with AVX-512
    # can be running the scalar kernels.  So the build identifies itself into
    # the results file rather than living only in the console output.
    build = None
    try:
        import biofasting

        build = dict(biofasting.build_info())
        build["seqops_level"] = biofasting.seqops_level()
        build["best_cpu_level"] = biofasting.best_cpu_level()
    except Exception as exc:  # pragma: no cover - environment dependent
        build = {"unavailable": f"{type(exc).__name__}: {exc}"}
    return {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cpu": cpu,
        "avx2": "avx2" in flags,
        "avx512": "avx512f" in flags,
        "backends": versions,
        "build": build,
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

    stored = {e["path"]: e.get("stored_bytes") for e in manifest["files"]}

    if not args.quick and (data / gz).exists():
        workloads.append(Workload(
            "fastq-gzip",
            f"{gz}, {sizes[gz] / 1e6:.1f} MB inflated from "
            f"{stored[gz] / 1e6:.1f} MB on disk, {counts[gz]:,} records",
            sizes[gz],  # MB/s on uncompressed size: same scale as fastq-plain
            make_fastq_impls(data / gz, backends) + make_decompression_impls(data / gz),
            counts[gz], source=gz))

    workloads.append(Workload(
        "fasta-scan", f"{fasta}, {sizes[fasta] / 1e6:.1f} MB", sizes[fasta],
        make_fasta_impls(data / fasta, backends), None, source=fasta))

    if not args.quick:
        # Deliberately few accesses: Bio.SeqIO.index has no slice-level random
        # access, so every access materialises the *whole* record -- on this
        # corpus chr1 alone is 40 Mbp.  A large access count measures the cost
        # of rebuilding 40 Mbp of Seq per 150 bp fetch and takes tens of
        # minutes per pass.  Latency per access is the honest metric here, and
        # the gap it shows is real (see bench/README.md).
        accesses, span = args.random_accesses, 150
        workloads.append(Workload(
            "fasta-random", f"{fasta}, {accesses:,} slices of {span} bp at fixed positions",
            accesses * span, make_fasta_random_impls(data / fasta, backends, accesses, span),
            accesses, source=fasta, per_access=True))

    #   FASTQ random access, once per name.  The access count is raised well
    #   above the FASTA default: a FASTQ record is four short lines, so one
    #   access is microseconds where a FASTA access is milliseconds, and a
    #   hundred of them would be lost in the timer's noise.
    index_accesses = max(args.random_accesses, 2000)
    names = fastq_record_names(data / fastq)
    index_build, index_random = make_fastq_index_impls(
        data / fastq, backends, names, index_accesses)
    workloads.append(Workload(
        "fastq-index", f"{fastq}, {len(names):,} records indexed by name",
        sizes[fastq], index_build, len(names), source=fastq))
    workloads.append(Workload(
        "fastq-random", f"{fastq}, {index_accesses:,} records fetched by name",
        index_accesses * 150, index_random, index_accesses, source=fastq,
        per_access=True))

    #   `--only` is applied after this function returns, so the BGZF companion
    #   -- twenty seconds and 174 MB of writing -- is asked for here rather than
    #   discovered later to have been unnecessary for a run that only wanted
    #   `fastq-plain`.
    def wanted(group):
        return args.only is None or group in set(args.only)

    if (not args.quick and (data / gz).exists()
            and (wanted("fastq-gz-index") or wanted("fastq-gz-random"))):
        #   The same pair, on the file a real pipeline actually has.  Only the
        #   reference's copy of it is created here, and it is created as setup
        #   (it is a fixture, not the thing being measured), so the number the
        #   reference row reports is its index build and not its re-encoding.
        #   The companion is written from the *plain* file: a BGZF copy of the
        #   `.gz` would be a copy of already-compressed bytes, which is not the
        #   content either row is supposed to be reading.
        index_path = bgzf_companion(data / fastq)
        gz_build, gz_random = make_fastq_gz_index_impls(
            data / gz, index_path, backends, names, index_accesses)
        #   `sizes[gz]` is the *inflated* size and is what the row's MB/s is
        #   computed against, so that it is on the same scale as `fastq-index`;
        #   the description says both numbers rather than calling the inflated
        #   one "on disk", which it is not and never was.
        workloads.append(Workload(
            "fastq-gz-index",
            f"{gz}, {sizes[gz] / 1e6:.1f} MB inflated from "
            f"{stored[gz] / 1e6:.1f} MB on disk, {len(names):,} records indexed by name",
            sizes[gz], gz_build, len(names), source=gz))
        workloads.append(Workload(
            "fastq-gz-random",
            f"{gz}, {index_accesses:,} records fetched by name from the compressed file",
            index_accesses * 150, gz_random, index_accesses, source=gz,
            per_access=True))

    Bio = backends["Bio"]
    reads = []
    for rec in Bio.SeqIO.parse(str(data / fastq), "fastq"):
        reads.append(rec)
        if len(reads) >= _OPS_READS:
            break
    ops_bytes = sum(len(rec) for rec in reads)
    revcomp, gc, translate = make_op_impls(reads, backends)
    workloads.append(Workload("ops-revcomp", f"reverse complement over {len(reads):,} reads",
                              ops_bytes, revcomp, len(reads), source="in memory"))
    workloads.append(Workload("ops-gc", f"GC fraction over {len(reads):,} reads",
                              ops_bytes, gc, len(reads), source="in memory"))
    workloads.append(Workload("ops-translate", f"translation over {len(reads):,} reads",
                              ops_bytes, translate, len(reads), source="in memory"))

    #   `Bio.SeqUtils`: the measurement workloads.  Their own, smaller slice of
    #   the same reads -- see SEQUTILS_READS -- and one contig of the FASTA, so
    #   that the ninety-second whole-chromosome row is measurable at all.  The
    #   contig is only read when its row was actually asked for: it is 40 Mbp,
    #   and a run that wanted `ops-gc123` should not pay for it.
    sequtils_reads = reads[:SEQUTILS_READS]
    sequtils_bytes = sum(len(rec) for rec in sequtils_reads)
    genome = None
    if wanted("sequtils-gc123-genome"):
        from Bio.SeqIO.FastaIO import SimpleFastaParser

        with _open_text(data / fasta) as handle:
            genome = max((seq.encode() for _, seq in SimpleFastaParser(handle)),
                         key=len)
    for group, description, impls, items in make_sequtils_impls(
            sequtils_reads, backends, genome):
        workloads.append(Workload(group, description, sequtils_bytes, impls,
                                  items, source="in memory"))

    #   `Bio.SeqUtils.ProtParam`: the protein rows, over a corpus that has to be
    #   generated -- `bench/data/` is DNA and holds no protein.  The generator is
    #   the same 192,800 residues for the same command on any machine, but it is
    #   only run when one of its groups was asked for: a `--only fastq-plain` run
    #   should not pay for it, and neither should one that only wanted `ops-gc`.
    if any(wanted(group) for group in PROTPARAM_GROUPS):
        protein_corpus = make_protein_corpus()
        protein_residues = sum(size for _, size in protein_corpus)
        for group, description, impls, items in make_protparam_impls(
                protein_corpus, backends):
            workloads.append(Workload(group, description, protein_residues, impls,
                                      items, source="generated protein corpus"))
    return workloads


def main(argv: list[str] | None = None) -> int:
    # Progress is printed as each row lands; without this, a redirected run
    # shows nothing until it finishes, which makes a slow run look hung.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):
        pass
    parser = argparse.ArgumentParser(prog="run.py", description="BioFasting benchmark runner")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--json", type=Path, default=None,
                        help="where to write results (default: bench/results/latest.json)")
    parser.add_argument("--repeat", type=int, default=5, help="timed runs per implementation")
    parser.add_argument("--quick", action="store_true", help="use the test-sized files")
    parser.add_argument("--random-accesses", type=int, default=100,
                        help="slices for the fasta-random workload (default: 100; the "
                             "Biopython row costs ~0.15 s per access, so this scales "
                             "the run's runtime almost linearly)")
    parser.add_argument("--max-pass-seconds", type=float, default=60.0,
                        help="skip repeats when a single pass exceeds this (default: 60)")
    parser.add_argument("--only", action="append", default=None,
                        help="run only these groups (repeatable)")
    parser.add_argument("--list", action="store_true", help="list workloads and exit")
    args = parser.parse_args(argv)

    backends: dict = {}
    versions: dict = {}
    # biofasting is optional like the rest: the runner still measures the
    # Phase 0 comparison without it, and a row for this project appears only
    # when the package being benchmarked is actually importable.
    for name in ("Bio", "pysam", "pyfaidx", "numpy", "biofasting"):
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
        results.extend(run_workload(workload, args.repeat,
                                    max_pass=args.max_pass_seconds))

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
    # ``valid`` is None for baselines (nothing to verify against a parse), so
    # only an explicit False is a correctness failure.
    invalid = [r for r in results if r.get("valid") is False]
    if invalid:
        print(f"WARNING: {len(invalid)} row(s) did not validate:")
        for row in invalid:
            print(f"  {row['group']}/{row['impl']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
