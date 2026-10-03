"""FASTQ reading: the file-level pieces around the compiled scanner.

The scanner itself is `_core.FastqScanner`, which works on any object exposing
the buffer protocol and therefore on a whole memory-mapped file at once.  What
cannot be expressed in a binding -- opening the file without copying it, and
keeping the mapping alive for exactly as long as the scanner needs it -- is
here.
"""

from typing import NamedTuple

from . import _core
from ._buffer import map_readonly

__all__ = [
    "FastqGrid",
    "open_fastq",
    "open_fastq_grid",
    "open_fastq_index",
    "scan_fastq",
]


class FastqGrid(NamedTuple):
    """A uniform FASTQ file as two zero-copy two-dimensional arrays.

    ``sequences`` and ``qualities`` are ``(records, length)`` arrays of
    ``uint8`` -- the file's own bytes, strided, nothing copied -- and
    ``title_length`` is the header's width from the record's ``@`` to its
    first base, which is the part of the row stride the two arrays do not
    account for.  Both arrays are read-only and both keep the file mapped, so
    there is nothing to close and nothing to copy.
    """

    sequences: object  # numpy.ndarray, (records, length), uint8
    qualities: object
    title_length: int

# The gzip magic number.  Sniffed rather than taken from the file name: an
# extension is a claim about the contents and these two bytes are the format.
_GZIP_MAGIC = b"\x1f\x8b"


def open_fastq(path):
    """Iterate the records of the FASTQ file at `path`, gzipped or not.

    Yields ``(title, sequence, quality)`` triples of :class:`bytes`, with the
    ``@`` removed from the title: the shape
    :func:`Bio.SeqIO.QualityIO.FastqGeneralIterator` returns, and the shape
    the benchmark compares against.

    A plain file is mapped rather than read, so a 369 MB file costs a few
    microseconds and no resident memory beyond the pages actually touched.  A
    gzipped file is mapped, inflated once with libdeflate, and the inflated
    bytes are what gets scanned -- see ``src/core/inflate.hpp`` for why it is
    whole-file inflation rather than a streaming feed.

    The mapping and any inflated buffer are kept alive by the returned object
    and released when it is dropped; there is nothing to close.
    """
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        return _core.FastqScanner.from_gzip(data)
    return _core.FastqScanner(data)


def open_fastq_grid(path):
    """Return a uniform FASTQ file as zero-copy arrays, or explain why not.

    Where :func:`open_fastq` yields a Python object per record, this makes
    none: ``sequences[i]`` is a view of the bytes of the i-th record, so a
    369 MB file costs one mapping and no resident memory beyond the pages
    touched, and a window like ``grid.sequences[:, :12]`` is an index rather
    than a copy.  Requires numpy, which is what the arrays are.

    A file is only a grid if *every* record is the same shape, and that is
    checked by reading them all -- not assumed from the first.  A ragged file
    raises :class:`ValueError` naming the first record that differs, because
    the alternative is an array that is quietly wrong; use :func:`open_fastq`
    for files that are not uniform.  Gzipped files work too, and the arrays
    then point into the inflated buffer.
    """
    data = map_readonly(path)
    scanner = (
        _core.FastqScanner.from_gzip(data)
        if data[:2] == _GZIP_MAGIC
        else _core.FastqScanner(data)
    )
    sequences, qualities, title_length = scanner.grid()
    return FastqGrid(sequences, qualities, title_length)


def scan_fastq(path):
    """Return the records of `path` as a list of triples.

    For tests and interactive use; `open_fastq` is the streaming form and the
    one to reach for on a large file.
    """
    return list(open_fastq(path))


def open_fastq_index(path):
    """Index the FASTQ file at `path` so that records can be fetched by name.

    The result behaves like :func:`Bio.SeqIO.index` with ``"fastq"``: a key is
    the header's first word, keys keep file order, a duplicate key raises
    :class:`ValueError`, and a missing key raises :class:`KeyError`.  What is
    different is what a fetch costs.  ``SeqIO.index`` stores a byte offset and
    re-parses the record from it every time -- 7.17 µs per access on the
    million-read corpus; this stores where each field *is*, so ``index[name]``
    is a memcpy and ``index.quality_length(name)`` does not touch the file at
    all.  Reading a record never builds a ``SeqRecord``; hand the bytes to
    :mod:`biofasting.interop` if that is what you need.

    Which files this accepts is the *parser*'s rule rather than ``SeqIO.index``'s,
    and the two are not the same: Biopython's index runs a second FASTQ grammar
    that rejects a blank line between records, which its own parser folds away.
    So this accepts ``bench/data/edge/blank_lines.fastq`` where ``SeqIO.index``
    refuses it, and it reports a malformed file with the parser's message and
    the record number and byte offset where it was found.  The header of
    ``src/core/fastq.hpp`` argues the choice; ``tests/test_fastq_index.py``
    asserts both halves of it.

    A plain file is mapped rather than read; a gzipped one is mapped, inflated
    once with libdeflate, and the inflated bytes are what gets indexed.  That
    costs the whole decompressed file in memory, and it buys the same thing it
    buys the parser: a fetch that is a memcpy.  The alternative -- seek into the
    compressed stream and inflate the enclosing block on every access, which is
    what ``SeqIO.index`` does on BGZF -- costs **274 µs per fetch** on the
    million-read corpus against 0.32 µs here, because random access into a
    compressed stream defeats the decompressor's block cache every time.

    Any gzip stream is accepted, BGZF included; BGZF is gzip with extra fields,
    so a ``.bgz`` file is the same code path.  ``SeqIO.index`` accepts only
    BGZF and raises ``ValueError: Gzipped files are not suitable for indexing,
    please use BGZF (blocked gzip format) instead`` for everything else, which
    is to say for the format most FASTQ files are actually in.

    The mapping and any inflated buffer are kept alive by the returned object;
    nothing is read lazily and closing the file afterwards is safe.
    """
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        return _core.FastqIndex.from_gzip(data)
    return _core.FastqIndex(data)
