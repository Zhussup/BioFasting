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
