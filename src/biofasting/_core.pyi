import mmap
from collections.abc import Iterator

__version__: str

def build_info() -> dict[str, str]: ...
def cpu_levels() -> list[str]: ...
def supported_cpu_levels() -> list[str]: ...
def best_cpu_level() -> str: ...

class FastqScanner:
    # Anything exposing the buffer protocol, which is why `mmap` is in the
    # union: mapping the file instead of reading it is the whole point.
    def __init__(self, source: bytes | bytearray | memoryview | mmap.mmap) -> None: ...
    @staticmethod
    def from_gzip(
        source: bytes | bytearray | memoryview | mmap.mmap,
    ) -> FastqScanner: ...
    def __iter__(self) -> FastqScanner: ...
    def __next__(self) -> tuple[bytes, bytes, bytes]: ...
    @property
    def record_index(self) -> int: ...
    @property
    def position(self) -> int: ...
    # Two read-only (records, length) uint8 arrays -- sequences and qualities
    # -- plus the header width, '@' included.  They are numpy arrays, so the
    # stub names them as `object` rather than importing numpy eagerly; the
    # module needs numpy only when this is called.
    def grid(self) -> tuple[object, object, int]: ...

class FastaIndex:
    def __init__(self, source: bytes | bytearray | memoryview | mmap.mmap) -> None: ...
    def __len__(self) -> int: ...
    def __contains__(self, name: str) -> bool: ...
    def __getitem__(self, name: str) -> bytes: ...
    def __iter__(self) -> Iterator[str]: ...
    def keys(self) -> list[str]: ...
    def title(self, name: str) -> str: ...
    def sequence_length(self, name: str) -> int: ...
    def strided(self, name: str) -> bool: ...
    def sequence_slice(self, name: str, start: int, end: int) -> bytes: ...
    def index_entry(self, name: str) -> tuple[str, int, int, int, int]: ...
    # A read-only (lines, line_bases) uint8 view of the record plus the
    # record's real length, which the array can be shorter than when the last
    # line is short.  numpy, likewise on demand.
    def grid(self, name: str) -> tuple[object, int]: ...

# One qualifier as the feature reader assembled it: (key, value, has_value,
# escape_warning, escape_text).  `has_value` false is the `/pseudo` form, whose
# value the reference stores as `""` and whose second occurrence it drops
# outright.  `escape_warning` and `escape_text` are the reference's NCBI escaping
# warning *reported* rather than emitted: it is about the value as it stands
# before the doubled quotes are undone, so the text has to travel with the flag.
QualifierTuple = tuple[str, str, bool, bool, str]

# (type, location_or_None, status, message, warnings, qualifiers).  `location`
# is a `LocationTuple` under the same status vocabulary `parse_location` uses and
# is `None` for `"parser_error"`, which the reference turns into a missing
# location plus a warning.
#
# `warnings` is the location's own in-parse warnings as `(kind, text)` pairs --
# `"origin_wrap"` carrying the part's text the reference quotes, `"bond"` with
# an empty text -- in the reference's own order, one pair per offending *part*.
FeatureTuple = tuple[
    str,
    "LocationTuple | None",
    str,
    str,
    list[tuple[str, str]],
    list[QualifierTuple],
]

class FlatFileIndex:
    # `format` is one of "genbank", "embl" or "swiss" -- the names
    # `SeqIO.parse` uses, so that a caller passes what it already has.
    def __init__(
        self, source: bytes | bytearray | memoryview | mmap.mmap, format: str
    ) -> None: ...
    def __len__(self) -> int: ...
    def __contains__(self, id: str) -> bool: ...
    def __getitem__(self, id: str) -> bytes: ...
    def __iter__(self) -> Iterator[str]: ...
    def keys(self) -> list[str]: ...
    def format(self) -> str: ...
    def name(self, id: str) -> str: ...
    def description(self, id: str) -> str: ...
    def sequence_length(self, id: str) -> int: ...
    def sequence_slice(self, id: str, start: int, end: int) -> bytes: ...
    # Every record in file order, so a stream parse pays no per-record lookup.
    # The same four fields as the per-id calls above, and nothing else: the
    # annotations dict and the FEATURES table are read on demand instead (see
    # `annotations` and `features` below).
    def records(self) -> list[tuple[str, str, str, bytes]]: ...
    # (offset, sequence_offset, sequence_end, length), as byte offsets into the
    # buffer the index was built over.
    def location(self, id: str) -> tuple[int, int, int, int]: ...
    def offsets(self) -> list[int]: ...
    def features(self, id: str) -> tuple[bool, str, list[FeatureTuple]]:
        """The record's FEATURES table, as `(ok, message, features)`.

        `ok` is False where this reader declined the table, and `message` then
        says which line did it.  It declines where the reference only *warns*
        and carries on -- a location that wraps its parentheses without breaking
        at a comma, a line too short to hold a feature, an over-indented
        location column, white space after a qualifier's `=`, a continuation
        with no qualifier above it -- because a repair reproduced differently is
        a feature that looks right and is not.  It also declines a table whose
        header facts it could not establish, since the locations would then be
        read under inputs the reference never used.

        A record with no feature block is `ok` with an empty list, which is an
        answer and not a refusal.  SwissProt's `FT` block is a refusal: it is a
        different grammar, and "no features" and "features not reproduced" are
        not the same statement.
        """
        ...
    def annotations(self, id: str) -> tuple[bool, str, list[str], dict]:
        """The record's header, as `(ok, message, keys, table)`.

        `keys` is the list of annotation keys the header *created*, in the order
        it created them.  That order is the file's, not a fixed list: the
        reference inserts a key when the line stating it is consumed, so a
        `COMMENT` above the first `REFERENCE` puts `comment` before
        `references`, and a dict compares and prints by its insertion order.
        Membership in `keys` is also what says a key *exists* -- EMBL creates
        `data_file_division` for every ID line even when the field is blank,
        where GenBank skips `molecule_type` on an empty column -- so an empty
        value in `table` means what the reference means by it rather than
        standing in for an absence.

        `table` holds the values, and is the reading rather than the naming:
        which of them become keys of `SeqRecord.annotations` is decided one
        level up, because GenBank and EMBL disagree about it.  `references` is a
        list of `(title, authors, consrtm, journal, pubmed_id, medline_id,
        comment, location)` with the location already in Python coordinates.

        `ok` is False where this reader declined the header, and `message` then
        names the shape it declined -- a structured comment, `NID`, `PID`,
        `DBSOURCE`, `SEGMENT`, a `LOCUS` line in the pre-229.0 layout, or a
        format whose header is not reproduced at all (SwissProt's).  It declines
        rather than return a dict missing a key the reference would have set: a
        dict that is wrong is worse than a dict that is absent.  The refusal is
        per record, so the sequence and the other records stay usable.
        """
        ...

# --- feature locations -------------------------------------------------------
#
# `Location.fromstring` and `Position.fromstring`, as kernels.  A location is
# not a pair of integers and these do not flatten one into a pair: `complement`
# carries a strand, `<`/`>` make an end fuzzy, `(3.9)` is a boundary known only
# to lie between two bases, `one-of(...)` is a choice, and `^` is a zero-length
# junction.
#
# A position is (kind, value, left, right, choices): `value` is what the
# reference's `int(position)` returns and what every comparison uses, `left`/
# `right` are the two edges of a `within` boundary, and `choices` is the
# `one_of` set.  `kind` indexes the enum in src/core/location.hpp, where the
# names are `exact before after within one_of uncertain unknown`.
PositionTuple = tuple[int, int, int, int, tuple[int, ...]]
PartTuple = tuple[PositionTuple, PositionTuple, int, str]
LocationTuple = tuple[str, tuple[PartTuple, ...]]

def parse_location(
    text: str,
    length: int | None = ...,
    circular: bool = ...,
    stranded: bool = ...,
) -> tuple[str, LocationTuple | None, str, list[tuple[str, str]]]:
    """Read a feature location.

    Returns `(status, location, message, warnings)`.  `warnings` is what the
    reference emits from inside the parse, *reported* rather than emitted: the
    kernel says what it found, the Python layer says the words.  It is a list of
    `(kind, text)` and not a pair of flags because the reference emits one
    warning per offending part -- `join(30..5,60..2)` on a circle repairs two
    parts and warns twice -- and because the order is observable: a `bond` part
    warns before a part after it that wraps the origin.  `kind` is
    `"origin_wrap"`, whose `text` is the part the reference quotes, or `"bond"`,
    whose wording is fixed and whose `text` is empty.

    The status is `"ok"`, `"parser_error"` -- where the reference raises
    `LocationParserError`, which its feature consumer catches and turns into a
    missing location plus a warning, so it is a behaviour and not just an
    error -- or `"refused"`, where the reference raises something else, or
    accepts the string while silently discarding text from it.  A refusal is a
    design difference: it is this reader saying it will not guess.

    `length` is the record's *declared* size and may be `None`, which is not the
    same fact as zero -- the origin-wrapping repair is skipped for a record with
    no size and taken for one that declares zero.
    """

def parse_position(
    text: str, offset: int = ...
) -> tuple[str, PositionTuple | None, str]:
    """Read one end of a location: `Position.fromstring(text, offset)`.

    `offset` is 0 for an end position and -1 for a start position, the
    reference's own convention; anything else is refused.
    """

class FastqIndex:
    def __init__(self, source: bytes | bytearray | memoryview | mmap.mmap) -> None: ...
    # `source` is a whole gzip stream, BGZF included; it is inflated with
    # libdeflate and the inflated bytes are what get indexed.  The result owns
    # that buffer, so it must not be built from a temporary.
    @staticmethod
    def from_gzip(source: bytes | bytearray | memoryview | mmap.mmap) -> FastqIndex: ...
    def __len__(self) -> int: ...
    def __contains__(self, name: str) -> bool: ...
    def __getitem__(self, name: str) -> bytes: ...
    def __iter__(self) -> Iterator[str]: ...
    def keys(self) -> list[str]: ...
    def title(self, name: str) -> str: ...
    def sequence_length(self, name: str) -> int: ...
    # Equal to sequence_length() by construction: the parser refuses a record
    # whose two lengths differ, so there is one length here, not two.
    def quality_length(self, name: str) -> int: ...
    def quality(self, name: str) -> bytes: ...
    def sequence_slice(self, name: str, start: int, end: int) -> bytes: ...
    # False for a record the scanner had to assemble (a wrapped sequence, a
    # blank line, a '+'-repeated header): still correct, but fetching it
    # re-parses its own span rather than copying it.
    def plain(self, name: str) -> bool: ...
    def index_entry(self, name: str) -> tuple[str, int, int, int, bool]: ...

class ArrowTable:
    """The record tables, built in Arrow's own layout (PLAN 2.3).

    Not built directly by callers: `biofasting.arrow.read_fastq_table` and
    `read_fasta_table` are the file-level pieces around it, the same split as
    `_core.FastqScanner` against `biofasting.fastq.open_fastq`.  This class does
    own its bytes, unlike every other reader here -- a column has to be
    contiguous and a file's records are not -- so it is the table, not the file,
    that `buffers()` hands out views of.
    """

    @staticmethod
    def from_fastq(
        source: bytes | bytearray | memoryview | mmap.mmap,
    ) -> ArrowTable:
        """(name, description, sequence, quality).

        Raises `ValueError` on a malformed record, with the parser's message
        and the record number and byte offset where it was found.
        """
        ...
    @staticmethod
    def from_fastq_gzip(
        source: bytes | bytearray | memoryview | mmap.mmap,
    ) -> ArrowTable: ...
    @staticmethod
    def from_fasta(
        source: bytes | bytearray | memoryview | mmap.mmap,
    ) -> ArrowTable:
        """(name, description, sequence).

        Raises `ValueError` on a duplicate key, with the reference's wording.
        """
        ...
    @property
    def names(self) -> list[str]: ...
    @property
    def length(self) -> int: ...
    @property
    def column_count(self) -> int: ...
    # One (offsets, data) pair a column, in table order: read-only numpy arrays
    # over the table's own bytes, which is exactly what
    # `pyarrow.Array.from_buffers` takes.  They hold the table alive.
    def buffers(self) -> list[tuple[object, object]]: ...

def seqops_level() -> str: ...
def reverse_complement(sequence: bytes | bytearray | memoryview) -> bytes: ...
def gc_fraction(sequence: bytes | bytearray | memoryview) -> float:
    """The GC fraction, or `-1.0` when nothing was counted at all.

    The sentinel exists because the reference's answer for an empty sequence is
    the *integer* `0`, and a C++ `double` cannot hold that.  `biofasting.seqops`
    turns it back into the integer; nothing else should ever see `-1.0`.
    """
    ...
def count_kmers(
    sequence: bytes | bytearray | memoryview, k: int
) -> list[tuple[bytes, int]]: ...

# --- Bio.SeqUtils: the functions that measure a sequence ------------------ #
#
# Every one of these returns *counts* or a raw pair, and the Python layer above
# decides what the answer is.  That is not a style choice: a measurement has a
# type, a rounding and a way of failing, and none of them survives a `double`.
# `gc123` cannot return four percentages because the reference raises on the
# fourth division; `gc_fraction` cannot say 0.0 for an empty sequence because
# the reference says the integer 0.

def gc123(sequence: bytes | bytearray | memoryview) -> list[int]:
    """Twelve counts, frame-major: frame f's A, T, G and C at 4f .. 4f+3.

    Both cases and nothing else, so a partial trailing codon is counted by
    nobody.  `biofasting.sequtils.GC123` turns these into the four numbers the
    reference returns -- and raises where the reference raises.
    """
    ...

def gc_skew(
    sequence: bytes | bytearray | memoryview, window: int
) -> list[float]:
    """`(G-C)/(G+C)` per non-overlapping window; the last window is short.

    A window with neither G nor C is `0.0`, which is the reference's
    `except ZeroDivisionError` arm.  `window == 0` raises `ValueError`, which is
    the reference's `range() arg 3 must not be zero`.
    """
    ...

def gc_counts(sequence: bytes | bytearray | memoryview) -> list[int]:
    """The counts behind `gc_fraction`'s "ignore" and "weighted" modes.

    Twelve numbers: the `"CGScgs"` count, the `"ATWUatwu"` count, then
    `count(x) + count(x.lower())` for x in `"BDHKMNRVXY"` in that order.
    """
    ...

def molecular_weight_mass(
    sequence: bytes | bytearray | memoryview, weights: bytes, valid: bytes
) -> tuple[float, int]:
    """The sum inside `molecular_weight`, Neumaier-compensated.

    `weights` is 256 little-endian doubles and `valid` 256 flags -- a sparse
    weight table handed over dense so that the per-letter branch is not in the
    loop.  Returns `(mass, first_bad)`, where `first_bad` is `-1` on success and
    otherwise the first byte the table has no weight for.

    The compensation is not decoration: CPython 3.12+ sums floats with Neumaier
    compensation and the reference is written as `sum(...)`, so a plain running
    total reproduces it only about three times in four.
    """
    ...

def gcg(sequence: bytes | bytearray | memoryview) -> int | None:
    """The GCG checksum, or `None` when the sequence is not ASCII.

    `None` is a decline, not a result: the reference upper-cases each character,
    and `"ß".upper()` is two characters, so a non-ASCII letter has to go through
    Python's Unicode tables.  The caller runs the reference's own loop then.
    """
    ...

def crc64(sequence: bytes | bytearray | memoryview) -> tuple[int, int]:
    """The crc64 halves, high first.

    No decline, unlike `gcg`: the reference folds each character in with
    `ord(c) & 0xFF`, so a byte is already the whole of what it reads.
    """
    ...

def crc64_table_h() -> list[int]:
    """The 256-entry table the kernel derives at compile time from the
    reference's own recurrence, exposed so a test can compare every entry
    against `CheckSum._table_h`."""
    ...

def translate(
    sequence: bytes | bytearray | memoryview,
    code: bytes,
    amino: bytes,
    stop_symbol: str,
    possible_stop: str,
    to_stop: bool,
    stop_is_error: bool,
) -> tuple[int, int, str]:
    """``(outcome, codon, protein)``.

    `outcome` is 0 for a completed translation, 1 when the kernel declined --
    `protein` then holds only the codons before `codon`, and the caller redoes
    the job -- and 2 when `stop_is_error` met a stop at `codon`.  `code` is 256
    bytes and `amino` 64; both come from `_translate._kernel_tables`.
    """
    ...

# --- the writers -------------------------------------------------------------
#
# One call for the whole file, not one per record: the reference already writes
# record by record, and the cost these remove is the interpreter crossing, so a
# kernel reached once a record would pay it back.  Each returns the whole file
# as one `bytes` and takes an iterable of tuples -- `(title, sequence)` for
# FASTA, `(title, sequence, quality)` for FASTQ, `(title, quality)` for QUAL --
# in this package's own form: bytes, the marker character excluded, and a
# quality string already phred+33.  Every function here is byte-for-byte
# Biopython's own output; the differential tests are the gate.

def write_fasta(
    records: Iterator[tuple[bytes, bytes]] | list[tuple[bytes, bytes]],
    wrap: int = 60,
) -> bytes:
    """The whole FASTA file, lines wrapped at `wrap`.

    `wrap` of 0 is the reference's unwrapped branch, which writes a newline even
    for an empty sequence; a wrapped empty record writes no base line at all.
    """
    ...

def write_fastq(
    records: Iterator[tuple[bytes, bytes, bytes]] | list[tuple[bytes, bytes, bytes]],
) -> bytes:
    """The whole FASTQ file: four lines a record, no wrapping.

    The quality string is copied rather than decoded, so the caller must already
    hold phred+33; a sequence and quality of different lengths is a
    ``ValueError``.
    """
    ...

def write_qual(
    records: Iterator[tuple[bytes, bytes]] | list[tuple[bytes, bytes]],
    wrap: int = 60,
) -> bytes:
    """The whole QUAL file: decimal scores cut at the last space in a window.

    That cut is the reference's *fast* branch, `data.rfind(" ", 0, wrap)` -- the
    one `SeqIO.write` takes -- and not the `pop(0)` loop of `to_string`.  Where
    the reference would loop forever -- no space anywhere in a window -- this
    raises ``ValueError`` instead.  `wrap` of 0 means one line a record.
    """
    ...

def phred_to_sanger(
    scores: bytes | bytearray | memoryview,
) -> bytes:
    """PHRED scores as the Sanger ASCII string, ``min(126, score + 33)`` a byte.

    Truncation at 93 is the reference's own behaviour on the slow path it takes
    for a score its 0..93 table has no entry for, and it is reproduced rather
    than corrected.
    """
    ...
