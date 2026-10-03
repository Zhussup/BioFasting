"""Flat-file reading: the INSDC formats, around the compiled index.

`_core.FlatFileIndex` does the work -- it takes any buffer, scans it once, and
then answers records by offset.  What cannot be expressed in a binding is
opening the file without copying it and keeping the mapping alive for exactly as
long as the index needs it, and that is this module.

**What this reproduces, and what it does not.**  A record here carries its `id`,
its `name`, its `description` and its sequence.  It does not carry the
annotations dict or the references, so this is not a drop-in for
``SeqIO.parse`` and must not be presented as one.  The FEATURES table sits one
level down, on ``_core.FlatFileIndex.features(id)``: it is read on demand, so a
caller that wants only sequences pays nothing for it, and a caller that wants
features pays only for the records it asks about.  Turning what comes back into
Biopython value types is the interop shim's job (``biofasting.seqfeature``),
which is why these records carry no features of their own.

The measurement behind it (bench/targets.md, sixth ranking pass) is what drew
the line there: on a 1 kb GenBank record with no features the bytes cost about
1.5 us and building the object tree about 1.0, while the reference spends some
25-28 us -- the rest of it its own line-by-line logic over the header block and
the `ORIGIN` block.  (That pass's absolute times were found not to reproduce on
a re-run and are quoted here as a shape rather than as figures to hold anything
to; the delivered numbers beside them were timed in one harness in one run.)

**Where it declines.**  Biopython meets a malformed sequence line by warning and
guessing -- a wrongly indented GenBank line is shifted by one byte and parsed
anyway -- and it meets a malformed feature the same way.  This reader refuses
such a file with the byte offset that made it refuse, because a guess reproduced
differently is a sequence that looks right and is not.  The feature table is
refused the same way but per record rather than per file, so one bad feature
costs its record and not the file.  A caller that meets either refusal should
fall back to ``SeqIO.parse``, which is a design difference and not a bug.
"""

from typing import NamedTuple

from . import _core
from ._buffer import map_readonly

__all__ = ["FLAT_FILE_FORMATS", "FlatFileRecord", "open_genbank", "read_genbank"]

# The names `SeqIO.parse` uses, which are the names this package takes.  The
# format is a parameter and not a guess: a GenBank record and an EMBL record
# share the `ID`/`LOCUS` vocabulary closely enough that a detector would be
# right most of the time, and most of the time is not a contract.
FLAT_FILE_FORMATS = ("genbank", "embl", "swiss")

_GZIP_MAGIC = b"\x1f\x8b"


class FlatFileRecord(NamedTuple):
    """A flat-file record: the four fields this reader reproduces.

    ``sequence`` is ``bytes``, not a ``Seq``, because building a ``Seq`` is the
    interop shim's job and doing it here would put a Biopython dependency in a
    module that does not need one.
    """

    id: str
    name: str
    description: str
    sequence: bytes


def _open(path, format):
    if format not in FLAT_FILE_FORMATS:
        raise ValueError(
            f"unknown flat-file format {format!r}: expected one of "
            + ", ".join(repr(name) for name in FLAT_FILE_FORMATS)
        )
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        raise ValueError(
            f"{path!r} is gzip-compressed; the flat-file reader does not read "
            "compressed files yet"
        )
    return _core.FlatFileIndex(data, format)


def open_genbank(path, format="genbank"):
    """Index the flat file at `path` and return it.

    The result behaves like a read-only mapping of record id -> sequence bytes.
    Ids keep file order and a duplicate id raises :class:`ValueError`, the same
    refusal `SeqIO.index` makes and for the same reason: an id that names two
    records addresses neither.

    ``index[id]`` is a copy out of the mapping, ``index.name(id)``,
    ``index.description(id)`` and ``index.sequence_length(id)`` are free, and
    ``index.sequence_slice(id, start, end)`` gathers just the window asked for.
    ``index.records()`` returns the whole file in one call as
    :class:`FlatFileRecord`s -- a caller iterating everything should not pay a
    lookup per record for an order it already knows.

    `format` is one of :data:`FLAT_FILE_FORMATS`.  The function is named for the
    format that carries most of the data rather than for all three, because
    ``open_genbank(path, format="embl")`` reads better than a name that pretends
    to be format-neutral.

    A plain file is mapped rather than read.  Compression is *not* handled here:
    a gzipped file raises :class:`ValueError` instead of being indexed as if it
    were text, because the alternative is an index that silently holds no
    records.
    """
    return _open(path, format)


def read_genbank(path, format="genbank"):
    """Return the records of `path` as a list of :class:`FlatFileRecord`.

    For tests and interactive use; `open_genbank` is the form to reach for on a
    large file, since this one materialises every sequence.
    """
    index = _open(path, format)
    return [FlatFileRecord(id, name, description, sequence)
            for id, name, description, sequence in index.records()]
