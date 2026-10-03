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

from . import fasta, fastq

__all__ = [
    "fasta_seqrecords",
    "fastq_seqrecords",
    "from_seqrecord",
    "to_seqrecord",
    "write_fasta",
    "write_fastq",
]

_PHRED_OFFSET = 33

# Biopython's own line width for a written FASTA record.  Named rather than
# written twice, because `SeqIO.write` is what the tests compare against and a
# mismatch here would be a silent one.
_FASTA_WRAP = 60


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
    return title, sequence, bytes(q + _PHRED_OFFSET for q in qualities)


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


def _write(handle, chunks):
    """Write `chunks` to a path or to an already-open binary handle."""
    if hasattr(handle, "write"):
        for chunk in chunks:
            handle.write(chunk)
        return
    with open(handle, "wb") as stream:
        for chunk in chunks:
            stream.write(chunk)


def write_fastq(records, handle):
    """Write our FASTQ records to `handle`, byte-identical to ``SeqIO.write``.

    `records` is any iterable of ``(title, sequence, quality)`` triples in our
    form -- what :func:`biofasting.open_fastq` yields, or what
    :func:`from_seqrecord` returns -- and `handle` is a path or an open binary
    stream.  Every record is written on four lines: Biopython does not wrap
    FASTQ, and neither does this.
    """
    def chunks():
        for title, sequence, quality in records:
            title = _as_bytes(title, "title")
            sequence = _as_bytes(sequence, "sequence")
            quality = _as_bytes(quality, "quality")
            if len(quality) != len(sequence):
                raise ValueError(
                    f"record {title!r} has sequence length {len(sequence)} "
                    f"but {len(quality)} quality scores"
                )
            yield b"@" + title + b"\n" + sequence + b"\n+\n" + quality + b"\n"

    _write(handle, chunks())


def write_fasta(records, handle):
    """Write our FASTA records to `handle`, byte-identical to ``SeqIO.write``.

    `records` is any iterable of ``(title, sequence)`` pairs in our form, where
    `title` excludes the ``>``.  Lines are wrapped at 60 bases, which is
    Biopython's default and the width its own writer produces.
    """
    def chunks():
        for title, sequence in records:
            title = _as_bytes(title, "title")
            sequence = _as_bytes(sequence, "sequence")
            yield b">" + title + b"\n"
            for start in range(0, len(sequence), _FASTA_WRAP):
                yield sequence[start : start + _FASTA_WRAP] + b"\n"

    _write(handle, chunks())
