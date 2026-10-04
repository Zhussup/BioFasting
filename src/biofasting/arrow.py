"""The record tables, as ``pyarrow.Table`` -- and the hand-off to polars.

A table can be built from this package without any of this: ``pa.table`` over
:func:`biofasting.open_fastq` works, and it costs one Python call per record.
The measurement that came first (``bench/rank_arrow.py``, recorded in
``bench/targets.md``) is what these two functions are an answer to, and it says
two things.  Against Biopython it is not close -- ``SeqIO.parse`` into lists is
15.8 µs a read on the million-read corpus, and the best pure-Python rewrite is
2.1 µs -- but the tool to beat is ``polars-bio``, whose Rust reader (needletail)
does the same job in **0.97 µs a read**, and whose columns a caller switching
between the two libraries will expect to find again here.  So the columns are
its columns, in its order, and the table is built in one pass in C++ into the
buffers Arrow itself would have used.

What comes back is a real ``pyarrow.Table``: ``pl.from_arrow(table)`` is the
documented hand-off, ``to_pandas()`` and ``to_polars()`` are one call away, and
every value in it is byte for byte what the readers yield -- the differential
tests compare the table against ``SeqIO.parse`` rather than against a fixture,
and the benchmark gates it against four independent producers.

**The hand-off is cheap but it is not free**, and the measurement that says so is
in ``bench/targets.md``.  ``pl.from_arrow`` does not copy the *strings* -- the
frame keeps the table's bytes alive after the table object is gone, which is what
building Arrow's own layout buys -- but polars 1.44 does not adopt an Arrow string
column as it stands either: it indexes it, which is O(rows) at 32-38 ns a record
on the million-read corpus.  A caller who never leaves ``pyarrow`` does not pay
that, and a caller who wants a frame pays it once per frame and not per column.

**The reading of the columns.**  ``name`` is the record's key, the title's first
word, exactly as ``SeqIO.index`` keys it and as this package's own indexes key
it.  ``description`` is the title *after* the key with the whitespace between
them trimmed -- ``>rec1 a description`` gives ``a description``, ``>rec1`` gives
``""`` -- which is ``polars-bio``'s column and **not** ``SeqRecord.description``,
which is the whole title with the key still in it.  That difference is worth a
sentence here rather than a surprise later, because both spellings are
reasonable and only one of them is this function's.  The sequence is the
reference's own bytes: a FASTA record's lines right-stripped, joined, and
stripped of spaces.

Two divergences from ``polars-bio`` are deliberate and are not bugs:

* a record with no description is ``""`` here and a **null** there.  A
  length-zero value and the absence of a value are different things, and this
  package has no nulls anywhere else;
* a header whose first byte after ``>`` is whitespace (``>  rec5 leading``)
  is read here, keyed ``rec5``; ``polars-bio`` refuses the whole file with
  ``FASTA read error: missing name``.  The reference reads it too, and a
  library that cannot open a file Biopython opens is a regression whatever it
  is compatible with.

What these do not do is grow a plugin: there is no polars IO source, no
DataFusion provider and no streaming or chunked reader.  A caller who wants
those wants ``polars-bio`` itself, and the point of the columns being identical
is that the two can be swapped rather than that this replaces it.

``pyarrow`` is an **optional** dependency (``pip install biofasting[arrow]``):
the tables are for a caller who already has Arrow, and nothing else in this
package needs it.  ``polars`` is not imported here at all -- ``pl.from_arrow``
is one call the caller makes.

Nothing was committed and nothing was pushed to produce this file; the working
tree is the deliverable and its owner commits it.
"""

from . import _core
from ._buffer import map_readonly

__all__ = ["available", "read_fastq_table", "read_fasta_table"]

# The gzip magic number.  Sniffed rather than taken from the file name: an
# extension is a claim about the contents and these two bytes are the format.
_GZIP_MAGIC = b"\x1f\x8b"


def available() -> bool:
    """Whether ``pyarrow`` can be imported.

    A caller choosing between this and ``polars-bio`` can ask this instead of
    catching an ``ImportError``; the readers in this package never need it.
    """
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        return False
    return True


def _pyarrow():
    """Import ``pyarrow`` or explain why it cannot be imported.

    The message names the install command rather than the import, for the same
    reason ``biofasting.alignment`` does: the interesting case is a caller who
    installed this package from a wheel and has never had Arrow, and is looking
    at a ``ModuleNotFoundError`` from inside a library they did not know wanted
    one.
    """
    try:
        import pyarrow
    except ImportError as exc:  # pragma: no cover - exercised by the ImportError test
        raise ImportError(
            "biofasting.arrow needs the optional pyarrow dependency, which is "
            "deliberately not a dependency of this package.  Install it with "
            "`pip install biofasting[arrow]` (or `pip install pyarrow`) to get "
            "the tables; the readers (open_fastq, open_fasta) remain the "
            "dependency-free way to read the same records."
        ) from exc
    return pyarrow


def _as_table(built):
    """A built ``_core.ArrowTable`` as a ``pyarrow.Table``, copying nothing.

    The columns are already in Arrow's layout -- a buffer of offsets and a
    buffer of concatenated bytes -- so ``Array.from_buffers`` is not a
    conversion, it is a wrapper around them.  The ``num_buffers``-shaped
    argument is ``[validity, offsets, data]`` and the validity buffer is
    ``None`` on purpose: nothing here has a null in it (see the module
    docstring for what ``polars-bio`` does differently), and a validity buffer
    of all-ones would be a megabyte of nothing.

    The arrays hold the numpy views alive, the views hold the table alive, and
    the table owns the bytes.  Dropping the table while the arrays are in use
    therefore releases nothing, which is the property ``buffers()`` exists for.
    """
    pa = _pyarrow()
    string = pa.large_string()
    arrays = [
        pa.Array.from_buffers(
            string,
            built.length,
            [None, pa.py_buffer(offsets), pa.py_buffer(data)],
        )
        for offsets, data in built.buffers()
    ]
    return pa.Table.from_arrays(arrays, names=built.names)


def read_fastq_table(path):
    """Read the FASTQ file at `path` into a ``pyarrow.Table``.

    Four columns, one row a record: ``name``, ``description``, ``sequence`` and
    ``quality`` -- ``polars-bio``'s ``read_fastq``, except that its quality
    column is called ``quality_scores``.  The quality is the ASCII string the
    file carries, phred+33, and not a list of integers: it is what the file
    holds and what every writer here takes back.

    A gzipped file works, and is sniffed by its magic bytes rather than by its
    name.  The records are parsed by the same scanner ``open_fastq`` uses, so
    the set of files this accepts is exactly the set ``FastqGeneralIterator``
    accepts and a malformed file raises the parser's :class:`ValueError` with
    the record number and byte offset where it was found.  Nothing is read
    twice: the parse, the split of the headers and the layout of the columns are
    one pass.

    The table owns its bytes, so it does not hold the file mapped and can
    outlive it; ``pl.from_arrow(table)`` hands them to polars without copying
    them (see the module docstring for what that call does cost).
    """
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        return _as_table(_core.ArrowTable.from_fastq_gzip(data))
    return _as_table(_core.ArrowTable.from_fastq(data))


def read_fasta_table(path):
    """Read the FASTA file at `path` into a ``pyarrow.Table``.

    Three columns, one row a record: ``name``, ``description`` and ``sequence``
    -- ``polars-bio``'s ``read_fasta``.  The sequences are the records' bytes as
    ``SeqIO.parse(handle, "fasta")`` builds them, with the newlines between the
    lines being the only thing removed that a reader might not expect.

    A duplicate key raises the reference's :class:`ValueError`, exactly as
    ``SeqIO.index`` does, and the whole table is refused rather than the
    duplicate being quietly kept twice.  Compression is not handled, for the
    same reason ``open_fasta`` does not handle it: a gzipped file raises instead
    of being indexed as if it were text, because the alternative is a table that
    silently holds no records.

    The sequences are measured by the index before the first one is copied, so
    the three buffers are sized once and the build does not reallocate -- which
    is where the gigabytes a second on a genome come from.
    """
    data = map_readonly(path)
    if data[:2] == _GZIP_MAGIC:
        raise ValueError(
            f"{path!r} is gzip-compressed; the FASTA table does not read "
            "compressed files yet"
        )
    return _as_table(_core.ArrowTable.from_fasta(data))
