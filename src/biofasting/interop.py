"""Interop with Biopython: our records and ``Bio.SeqRecord``, both directions.

BioFasting is a substrate, not a replacement, so the boundary between its
records and Biopython's has to be cheap and exact in both directions.  That is
what this module is: ``fastq_seqrecords``/``fasta_seqrecords`` turn our readers
into ``SeqRecord`` iterators without a re-parse, ``to_seqrecord`` converts a
record you are already holding, ``from_seqrecord`` goes the other way, and the
two writers turn our records back into FASTQ and FASTA bytes.

Biopython is **not** a dependency of this package, and this module does not
change that: it imports ``Bio`` lazily inside the functions that need it, and
importing :mod:`biofasting.interop` imports nothing.  A user who never touches
a ``SeqRecord`` never pays for one.

The formats are Biopython's, byte for byte: a title excludes its ``@`` or ``>``,
``id`` is the title's first word, ``description`` is the whole title line, and a
quality string is phred+33.  Everything here is differentially tested against
``Bio.SeqIO`` on the corpus, writers included.
"""

import warnings

from . import _core, fasta, fastq
from ._translate import BiopythonWarning

__all__ = [
    "fasta_seqrecords",
    "fastq_seqrecords",
    "from_seqrecord",
    "to_seqrecord",
    "write_fasta",
    "write_fastq",
    "write_qual",
]

_PHRED_OFFSET = 33

#: The last byte a Sanger quality string may hold -- ``~``, PHRED 93 -- and the
#: highest PHRED score that byte encodes.  A score above the ceiling cannot be
#: written, only truncated, and the reference warns when it has to.
_SANGER_MAX = 126
_SANGER_CEILING = 93

#: The reference's own message for that truncation, verbatim.
_DATA_LOSS = "Data loss - max PHRED quality 93 in Sanger FASTQ"

# Biopython's own line width for a written FASTA record.  Named rather than
# written twice, because `SeqIO.write` is what the tests compare against and a
# mismatch here would be a silent one.
_FASTA_WRAP = 60

# And the same for QUAL, which has a width of its own because its lines hold
# numbers rather than bases -- two-digit scores and a space are four characters
# a base, so 60 columns is a line of fifteen scores and not sixty.
_QUAL_WRAP = 60


def _require_biopython():
    """Import ``Bio`` or explain why it cannot be imported.

    The message names the fix rather than the import, because the interesting
    case is a user who installed biofasting from a wheel, has never had
    Biopython, and is looking at a ``ModuleNotFoundError: No module named
    'Bio'`` from inside a library they did not know needed it.
    """
    try:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError as exc:  # pragma: no cover - exercised by the ImportError test
        raise ImportError(
            "biofasting.interop needs Biopython, which is deliberately not a "
            "dependency of this package.  Install it with `pip install "
            "biopython` to convert between our records and Bio.SeqRecord."
        ) from exc
    return Seq, SeqRecord


def _as_bytes(value, what):
    """Accept `str` or bytes-like and return bytes; reject anything else."""
    if isinstance(value, str):
        return value.encode("latin-1")
    if isinstance(value, (bytes, bytearray, memoryview)):
        return bytes(value)
    raise TypeError(f"{what} must be str or bytes-like, not {type(value).__name__}")


def _title_of(record):
    """The title line Biopython's writers would emit for `record`.

    Copied from ``FastaWriter``/``FastqPhredWriter`` deliberately, placeholder
    defaults and all: a record built without an ``id`` writes as ``<unknown id>
    <unknown description>`` in Biopython, and a shim that quietly tidied that
    away would produce bytes the reference never would.
    """
    identifier = record.id or ""
    description = record.description
    if description and description.split(None, 1)[0] == identifier:
        return description
    if description:
        return f"{identifier} {description}"
    return identifier


def to_seqrecord(title, sequence, quality=None):
    """Return a ``Bio.SeqRecord.SeqRecord`` for one of our records.

    `title` is a header line without its ``@``/``>``, `sequence` its bases, and
    `quality` an optional phred+33 string -- omit it for a FASTA record.  Each
    may be `str` or bytes-like; the result carries bytes internally, which is
    what our readers hand out and what ``Bio.Seq`` accepts.

    The result is the record ``Bio.SeqIO.parse`` would have produced from the
    same file: ``id`` and ``name`` are the title's first word, ``description``
    is the whole title, and a FASTQ record carries ``phred_quality``.
    """
    Seq, SeqRecord = _require_biopython()
    title = _as_bytes(title, "title").decode("latin-1")
    # The id is the title's first word, which is the rule SeqIO applies to
    # both formats.  An empty title has no first word and gets an empty id
    # rather than an IndexError.
    words = title.split(None, 1)
    identifier = words[0] if words else ""
    record = SeqRecord(
        Seq(_as_bytes(sequence, "sequence")),
        id=identifier,
        name=identifier,
        description=title,
    )
    if quality is not None:
        record.letter_annotations["phred_quality"] = [
            byte - _PHRED_OFFSET for byte in _as_bytes(quality, "quality")
        ]
    return record


def from_seqrecord(record):
    """Return one of our records for a ``Bio.SeqRecord``: a tuple of bytes.

    ``(title, sequence, quality)`` for a record with ``phred_quality`` -- the
    shape :func:`biofasting.open_fastq` yields -- and ``(title, sequence)``
    for one without, which is the FASTA shape.  The title excludes the marker
    character, as everywhere else in this package.
    """
    title = _title_of(record).encode("latin-1")
    sequence = str(record.seq).encode("latin-1")
    qualities = record.letter_annotations.get("phred_quality")
    if qualities is None:
        return title, sequence
    return title, sequence, _sanger_quality(qualities)


def _sanger_quality(qualities):
    """`phred_quality` as the Sanger ASCII string, in one C pass.

    ``(title, sequence, quality)`` is the shape every writer here takes, so this
    is where a ``SeqRecord`` enters it, and it is the one place in the
    conversion where the reference spends a *dictionary lookup a base*:
    ``_get_sanger_quality_str`` builds its answer with
    ``"".join(_phred_to_sanger_quality_str[qp] for qp in qualities)``.  The
    kernel is the same mapping, ``min(126, qp + 33)``, a byte at a time.

    Two things the kernel cannot do, and both are the reference's own behaviour
    rather than an approximation of it:

    - **The warning.** A Sanger FASTQ holds PHRED 0..93 and no more, so a score
      of 94 or above is truncated to ``~`` *and warned about*
      (``max(qualities) >= 93.5``).  The test suite runs with warnings as
      errors, so a caller who would see the reference refuse sees this refuse
      too.
    - **Scores that are not bytes.** A float, a ``None`` or a value above 255
      cannot go through the kernel, and those keep the reference's slow
      expression, ``chr(min(126, int(round(qp)) + SANGER_SCORE_OFFSET))`` -- the
      rounding included, which is why a float is not simply cast.
    """
    try:
        raw = bytes(qualities)
    except (TypeError, ValueError):
        return _sanger_quality_slow(qualities)
    if raw and max(raw) > _SANGER_CEILING:
        warnings.warn(_DATA_LOSS, BiopythonWarning, stacklevel=3)
    return _core.phred_to_sanger(raw)


def _sanger_quality_slow(qualities):
    """The reference's fallback path, expression for expression.

    Reached where ``bytes(qualities)`` fails, which is exactly where the
    reference's cached table lookup fails: a float, a ``None``, a score outside
    a byte.  The order of the three checks is the reference's, because a list
    holding a ``None`` has to raise its own message rather than a
    ``TypeError`` about comparing ``None`` to an int.
    """
    if None in qualities:
        raise TypeError("A quality value of None was found")
    if max(qualities) >= 93.5:
        warnings.warn(_DATA_LOSS, BiopythonWarning, stacklevel=3)
    return bytes(
        min(_SANGER_MAX, int(round(qp)) + _PHRED_OFFSET) for qp in qualities
    )


def fastq_seqrecords(path):
    """Iterate a FASTQ file's records as ``SeqRecord`` objects.

    ``Bio.SeqIO.parse(path, "fastq")`` with our scanner underneath: the same
    objects, one at a time, without Biopython parsing the file again.
    """
    for title, sequence, quality in fastq.open_fastq(path):
        yield to_seqrecord(title, sequence, quality)


def fasta_seqrecords(path):
    """Iterate a FASTA file's records as ``SeqRecord`` objects.

    The index is built once and each record is copied out of the mapping as it
    is asked for, so a 100 Mbp genome is not materialised to walk it.
    """
    index = fasta.open_fasta(path)
    for name in index:
        yield to_seqrecord(index.title(name), index[name])


def _write(handle, data):
    """Write one `bytes` to a path or to an already-open binary handle.

    It was a loop over chunks until the writers got kernels: the whole point of
    them is that the file is built in one buffer, and a loop that wrote it back
    out in pieces would put the crossings back one level up.
    """
    if hasattr(handle, "write"):
        handle.write(data)
        return
    with open(handle, "wb") as stream:
        stream.write(data)


def _wrap_width(wrap):
    """The reference's own validation of a wrap width, before the kernel sees it.

    ``FastaWriter.__init__`` raises ``ValueError`` for a negative `wrap` and
    treats ``0`` and ``None`` as "do not wrap"; the kernel takes an unsigned
    width, so a negative one has to be refused here or it would arrive as an
    enormous number and produce a single line of nonsense.
    """
    if wrap is None:
        return 0
    if wrap < 0:
        raise ValueError("wrap must not be negative")
    return wrap


def write_fastq(records, handle):
    """Write our FASTQ records to `handle`, byte-identical to ``SeqIO.write``.

    `records` is any iterable of ``(title, sequence, quality)`` triples in our
    form -- what :func:`biofasting.open_fastq` yields, or what
    :func:`from_seqrecord` returns -- and `handle` is a path or an open binary
    stream.  Every record is written on four lines: Biopython does not wrap
    FASTQ, and neither does this.

    The whole file is built by the kernel in one buffer and written once, so a
    million-read file costs one crossing into the handle rather than a million.
    The length check stays here rather than in the kernel for the message it
    carries -- the kernel keeps its own as a backstop for a direct caller.
    """
    def rows():
        for title, sequence, quality in records:
            title = _as_bytes(title, "title")
            sequence = _as_bytes(sequence, "sequence")
            quality = _as_bytes(quality, "quality")
            if len(quality) != len(sequence):
                raise ValueError(
                    f"record {title!r} has sequence length {len(sequence)} "
                    f"but {len(quality)} quality scores"
                )
            yield title, sequence, quality

    _write(handle, _core.write_fastq(rows()))


def write_fasta(records, handle, wrap=_FASTA_WRAP):
    """Write our FASTA records to `handle`, byte-identical to ``SeqIO.write``.

    `records` is any iterable of ``(title, sequence)`` pairs in our form, where
    `title` excludes the ``>``.  Lines are wrapped at `wrap`, 60 by default,
    which is Biopython's default and the width its own writer produces.  A
    `wrap` of 0 or `None` writes each sequence on one line, which is the
    reference's other branch and not a special case of this one: it writes a
    newline for an empty sequence, where a wrapped empty record writes no base
    line at all.
    """
    width = _wrap_width(wrap)

    def rows():
        for title, sequence in records:
            yield _as_bytes(title, "title"), _as_bytes(sequence, "sequence")

    _write(handle, _core.write_fasta(rows(), width))


def write_qual(records, handle, wrap=_QUAL_WRAP):
    """Write our records' quality scores to `handle` as QUAL, byte-identical.

    `records` is any iterable of ``(title, quality)`` pairs where `quality` is
    the phred+33 string a FASTQ record carries -- the third field of what
    :func:`biofasting.open_fastq` yields, and what :func:`from_seqrecord`
    returns.  The scores are written as decimal, which is what QUAL is: a FASTA
    header line and the numbers, with no bases at all.

    This is the writer the reference is worst at: it formats every score with
    ``"%i" % round(q, 0)`` and cuts the lines with ``data.rfind(" ", 0, wrap)``,
    103 ns a base on the corpus.  The kernel does the same arithmetic into a
    buffer, and the cut is the same cut.

    A `wrap` between 1 and 5 takes a branch of the reference the kernel does not
    have -- its "safe wrapping" `pop(0)` loop, kept for widths too narrow for
    the fast one -- and that branch is reproduced here in Python.  It is a shape
    `SeqIO.write` never produces at its default width, and having it in the
    kernel would be a third rule for a case no caller reaches.
    """
    width = _wrap_width(wrap)
    if 0 < width <= 5:
        return _write_qual_narrow(records, handle, width)

    def rows():
        for title, quality in records:
            yield _as_bytes(title, "title"), _as_bytes(quality, "quality")

    _write(handle, _core.write_qual(rows(), width))


def _write_qual_narrow(records, handle, wrap):
    """The reference's `elif wrap:` branch: scores packed by popping a list.

    Reached only for a width of one to five, where `write_record` does not take
    its `rfind` path.  The rule is its own and not a narrower version of the
    fast one, so it is written out rather than approximated.
    """
    out = []
    for title, quality in records:
        title = _as_bytes(title, "title")
        quality = _as_bytes(quality, "quality")
        out.append(b">" + title + b"\n")
        tokens = [str(byte - _PHRED_OFFSET) for byte in quality]
        while tokens:
            line = tokens.pop(0)
            while tokens and len(line) + 1 + len(tokens[0]) < wrap:
                line += " " + tokens.pop(0)
            out.append((line + "\n").encode("latin-1"))
    _write(handle, b"".join(out))
