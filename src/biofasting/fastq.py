"""FASTQ reading: the file-level pieces around the compiled scanner.

The scanner itself is `_core.FastqScanner`, which works on any object exposing
the buffer protocol and therefore on a whole memory-mapped file at once.  What
cannot be expressed in a binding -- opening the file without copying it, and
keeping the mapping alive for exactly as long as the scanner needs it -- is
here.
"""

import os
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
    microseconds and no resident memory beyond the pages actually touched --
    and the pages that have been passed are given back to the kernel as the
    scan goes, so iterating a file of any size holds one window of it (32-64 MB)
    rather than the file.  A gzipped file is mapped, inflated once with
    libdeflate, and the inflated bytes are what gets scanned -- see
    ``src/core/inflate.hpp`` for why it is whole-file inflation rather than a
    streaming feed.

    The mapping and any inflated buffer are kept alive by the returned object
    and released when it is dropped; there is nothing to close.
    """
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        return _core.FastqScanner.from_gzip(data)
    return _core.FastqScanner(data, mapped=True)


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


def open_fastq_index(path, *, cache=True):
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

    With ``cache=True`` (the default) the index is also written to
    ``path + ".bfqi"`` and memory-mapped on the next open, so that opening a
    file that has been indexed once costs a mapping rather than a scan: half a
    second becomes a fraction of a millisecond at a million reads.  The file is
    a cache and nothing else -- it carries the data's size and a fingerprint of
    its first and last pages, an index that does not match is ignored and
    rebuilt, and an unwritable directory is not an error.  Pass ``cache=False``
    to keep the filesystem untouched.

    Which files this accepts is the *parser*'s rule rather than ``SeqIO.index``'s,
    and the two are not the same: Biopython's index runs a second FASTQ grammar
    that rejects a blank line between records, which its own parser folds away.
    So this accepts ``bench/data/edge/blank_lines.fastq`` where ``SeqIO.index``
    refuses it, and it reports a malformed file with the parser's message and
    the record number and byte offset where it was found.  The header of
    ``src/core/fastq.hpp`` argues the choice; ``tests/test_fastq_index.py``
    asserts both halves of it.

    What the build costs in memory is the index, not the file: a plain file is
    scanned through its mapping and the pages already scanned are dropped behind
    the scan, so peak memory is the entry table plus a scan window rather than
    the input's size (369 MB of reads build in about 210 MB).  A gzipped file is
    the exception, and it is worth stating plainly: the build holds the whole
    decompressed stream, so a cold index over a 100 GB archive wants 100 GB of
    RAM.  With ``cache=True`` that happens once -- afterwards the index is over
    the ``.bfqd`` copy and opens by mapping it -- but the first build of a very
    large archive is the case this design does not yet cover.

    Inflating the whole archive rather than seeking into it block by block buys
    the same thing for the index as for the parser: a fetch that is a memcpy.
    The alternative -- seek into the compressed stream and inflate the enclosing
    block on every access, which is what ``SeqIO.index`` does on BGZF -- costs
    **274 µs per fetch** on the million-read corpus against 0.32 µs here,
    because random access into a compressed stream defeats the decompressor's
    block cache every time.

    With ``cache=True`` the inflated bytes are what gets written to
    ``path + ".bfqd"``, so that the next open is a mapping rather than a
    decompression: 10 ms becomes 0.18 ms, and a fetch from the ``.gz`` costs what
    a fetch from the plain file costs.  Be aware of the size: that file is the
    *decompressed* data, two to four times the archive, and it exists only
    because an index over inflated offsets is useless without the inflated
    bytes.  The index records the archive's size and fingerprint, so replacing
    the ``.gz`` invalidates both files.  Pass ``cache=False`` for large archives
    you do not want a second copy of.

    Any gzip stream is accepted, BGZF included; BGZF is gzip with extra fields,
    so a ``.bgz`` file is the same code path.  ``SeqIO.index`` accepts only
    BGZF and raises ``ValueError: Gzipped files are not suitable for indexing,
    please use BGZF (blocked gzip format) instead`` for everything else, which
    is to say for the format most FASTQ files are actually in.

    The mapping and any inflated buffer are kept alive by the returned object;
    nothing is read lazily and closing the file afterwards is safe.
    """
    data = map_readonly(path)
    gzipped = data[:2] == _GZIP_MAGIC
    index_path = os.fspath(path) + _INDEX_SUFFIX
    data_path = os.fspath(path) + _DATA_SUFFIX
    if cache:
        cached = _open_cached_index(data, index_path, data_path if gzipped else None)
        if cached is not None:
            return cached
    if gzipped:
        index = _core.FastqIndex.from_gzip(data)
    else:
        # `mapped=True`: the bytes came from `map_readonly`, so the build may
        # give back the pages it has scanned and the fetches that follow are
        # random-access rather than readahead.  Peak memory then follows the
        # index, not the file.
        index = _core.FastqIndex(data, mapped=True)
    if cache:
        try:
            if gzipped:
                index.save_index(index_path, store_data=data_path, source=data)
            else:
                index.save_index(index_path)
        except (OSError, RuntimeError):
            # A cache that cannot be written is not a reason to refuse a file
            # that was read successfully.  Read-only directories, a full disk
            # and a directory the user cannot write are all ordinary.
            pass
    return index


# The suffixes of the two files `open_fastq_index` writes.  Next to the data and
# not in a hidden cache directory, because an index is only valid for exactly
# one file's bytes and moving or copying the data must invalidate it visibly.
_INDEX_SUFFIX = ".bfqi"
_DATA_SUFFIX = ".bfqd"


def _open_cached_index(data, index_path, data_path):
    """The index at `index_path` over `data`, or None if it cannot be used.

    `None` is the answer to several different situations -- a file is not
    there, it describes something else, it is damaged, or the decompressed copy
    it needs is missing -- because the caller does the same thing in all of
    them: build the index and overwrite the cache.
    """
    try:
        index_data = map_readonly(index_path)
    except (FileNotFoundError, NotADirectoryError):
        return None
    if data_path is None:
        try:
            return _core.FastqIndex.from_index(data, index_data)
        except ValueError:
            return None
    try:
        inflated = map_readonly(data_path)
    except (FileNotFoundError, NotADirectoryError):
        return None
    try:
        return _core.FastqIndex.from_index(inflated, index_data, data)
    except ValueError:
        return None
