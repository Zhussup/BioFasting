"""FASTA reading: the file-level pieces around the compiled index.

`_core.FastaIndex` does the work -- it takes any buffer, scans it once, and
then reads records out of it by offset.  What cannot be expressed in a binding
is opening the file without copying it and keeping the mapping alive for
exactly as long as the index needs it, and that is this module.
"""

from typing import NamedTuple

from . import _core
from ._buffer import map_readonly

__all__ = ["FastaGrid", "open_fasta", "open_fasta_grid", "read_fasta"]


class FastaGrid(NamedTuple):
    """A FASTA record as a zero-copy two-dimensional array.

    ``bases`` is a read-only ``(lines, line_bases)`` ``uint8`` view of the
    file's own bytes -- a chromosome costs a mapping and a shape, not a copy.
    It covers whole lines only, and ``length`` is the record's real length, so
    a trailing short line shows up as ``bases.size != length`` instead of
    disappearing.
    """

    bases: object  # numpy.ndarray, (lines, line_bases), uint8
    length: int


_GZIP_MAGIC = b"\x1f\x8b"


def open_fasta(path):
    """Index the FASTA file at `path` and return it.

    The result behaves like a read-only mapping of record key -> sequence
    bytes, following :mod:`Bio.SeqIO`: a key is the title's first word, keys
    keep file order, and a duplicate key raises :class:`ValueError`.  Unlike
    ``SeqIO.index``, reading a record afterwards does not parse it again --
    ``index[name]`` is a copy out of the mapping, ``index.sequence_length(name)``
    is free, and ``index.sequence_slice(name, start, end)`` is a gather.

    A plain file is mapped rather than read.  Compression is *not* handled
    here: a gzipped file raises :class:`ValueError` instead of being indexed as
    if it were text, because the alternative is an index that silently holds no
    records.
    """
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        raise ValueError(
            f"{path!r} is gzip-compressed; the FASTA index does not read "
            "compressed files yet"
        )
    return _core.FastaIndex(data)


def open_fasta_grid(path, name):
    """Return record `name` of `path` as a zero-copy array, or explain why not.

    The record must be a regular grid -- every line the same width, which is
    how FASTA is written in practice and how the corpus genome is written.  The
    result is a :class:`FastaGrid`: ``bases`` is a read-only
    ``(lines, line_bases)`` array of ``uint8`` pointing straight at the mapped
    file, so ``bases.reshape(-1)`` after a ``bases.copy()``, or any vectorised
    work over ``bases``, costs no parse at all.  Requires numpy.

    The one thing to know: ``bases`` holds whole lines.  A record whose length
    is not a multiple of its line width -- 40,000,000 bases at 60 per line --
    ends with a short line that is *not* in the array, and ``length`` is there
    so the difference is visible.  The tail is available from
    ``open_fasta(path).sequence_slice(name, bases.size, length)``.

    A record whose lines are not regular raises :class:`ValueError`; use
    ``open_fasta(path)[name]`` for those.
    """
    bases, length = open_fasta(path).grid(name)
    return FastaGrid(bases, length)


def read_fasta(path):
    """Return the records of `path` as a list of (name, title, sequence).

    For tests and interactive use; `open_fasta` is the form to reach for on a
    large file, since this one materialises every sequence.
    """
    index = open_fasta(path)
    return [(name, index.title(name), index[name]) for name in index]
