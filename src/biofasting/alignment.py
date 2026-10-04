"""`Bio.Align.PairwiseAligner`'s alignment, on parasail's SIMD kernels.

The plan files alignment as a wrapper rather than a kernel -- reuse an optimized
Smith-Waterman instead of writing one -- and the premise was measured before
this module was written (`bench/rank_align.py`, the eighth ranking pass).  It is
true, but not for the reason the plan implies: `Bio.Align.PairwiseAligner` in
1.88 is a C extension, a scalar DP at 3.1-3.5 ns a cell, not a Python loop.
What parasail adds is the SIMD width -- 0.17-0.52 ns a cell -- and the measured
result is 6x-35x score-only and 9x-38x with traceback over six rows, gate first.

Two things had to be established by measuring, and both are recorded in
`bench/targets.md` because a wrapper that gets either wrong is worse than none:

* **Which parasail function.**  The unsuffixed `nw_scan`/`sg_scan`/`sw_scan`
  bindings are the generic dispatch at ~2.2 ns a cell; the explicitly sized
  ones are the SIMD kernels.  Nothing here calls an unsuffixed binding: the
  width comes from the scheme's own score bound and `saturated` is checked on
  the result, because an 8-bit register that overflows returns a wrong score
  rather than an error.
* **What Biopython's end-gap scores mean.**  Its global mode prices a terminal
  gap with `end_gap_score`/`extend_end_gap_score`, which may differ from the
  interior ones -- so "global" is two different alignments depending on them.
  parasail has exactly two of those: `nw` (terminal gaps priced like interior)
  and `sg` (both ends free).  A scheme in between is refused rather than
  approximated, and that includes the reference's *own* default, where
  `end_gap_score = -1` and `open_gap_score = -1` but `extend_end_gap_score` is
  not `extend_gap_score`... it is, in fact, exactly equal, which is why the
  default scheme maps to `nw`; a caller who sets only one of them gets an
  error, not a different alignment.  Local mode is insensitive to both
  (measured, 0 of 63 pairs), which is why `sw` is used there.

The free-end scheme carries one correction, and it is measured rather than
assumed.  With all four ends free, the best alignment can be one with **no
substitution column at all** -- `'AAAA'` against `'TTTT'` aligns as
`'AAAA----'` / `'----TTTT'`, every letter in a free gap, score 0 -- and
parasail's recurrence has no such path: an insertion run adjacent to a deletion
run is not representable, so its `sg` answers -2 there.  What it does compute
is the best alignment *with* a substitution column, so the reference's score is
``max(sg, 0)``: over 3,000 random pairs, including 51 whose `sg` score is
negative, that expression equals `Bio.Align.PairwiseAligner`'s score 3,000
times out of 3,000.  A negative `sg` score therefore returns that no-column
alignment explicitly rather than parasail's path, and it is the same alignment
Biopython reports (its own second optimum, in the tests).  Everywhere else --
including every pair with an overlap worth aligning, which is every pair a real
caller has -- `sg` is the score unchanged.

The contract, and its honest limit:

* **Scores are bit-identical to the reference** for every scheme this module
  accepts.  This is what the tests sweep.
* **Tracebacks are optimal and independently tie-broken.**  The reference
  yields every optimal alignment for a pair, in its own DP's order; bit
  identity of a *path* is therefore not a well-defined contract and is not
  claimed here.  `Alignment.counts().score` recomputes the path's score from
  the aligned letters and the matrix -- an independent implementation that must
  equal the kernel's own score -- and wherever the reference has exactly one
  optimal alignment, this module's path equals it column for column.  Measured
  over 3,000 random pairs in all three schemes: our path always *is* one of the
  reference's optima, and its `identities`/`mismatches`/`gaps` always equal that
  optimum's.  That recomputation is a Python step per aligned letter, so it
  happens the first time `counts()` is asked for and not when the alignment is
  built: on a 150-base pair it is half of what `align` costs, and a caller who
  asked for the path did not ask for it.
* **A local alignment with nothing to align is the empty one, not nothing.**
  When no pair of letters scores above zero the reference's `align` yields
  *zero* alignments while its `score` still says 0.0; `sw` reports a zero-length
  alignment instead, so this module returns an `Alignment` with empty gapped
  strings, no coordinates and score 0 -- the same score, an explicit object
  rather than an empty iterator.

`parasail` is an **optional** dependency (`pip install biofasting[alignment]`)
rather than a required one, for a measured packaging reason: parasail 1.3.4
ships wheels for x86_64 Linux, x86_64 macOS, win32 and win_amd64 and **no
aarch64 wheel of any kind**, while this project's CI runs `macos-14` and builds
an aarch64 wheel.  A required dependency would mean compiling parasail from
source on every one of those jobs, so the accelerator is opt-in and the speedup
is honestly an x86-64 speedup for now.

`Bio` is imported only inside :meth:`Alignment.to_biopython`, and `parasail`
only inside the functions that need it, so importing this module costs nothing
and the readers do not start depending on either.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

__all__ = ["Alignment", "Aligner", "Counts", "align", "available", "score"]

# The ops parasail's decoded CIGAR uses: '=' and 'X' consume a letter from both
# sequences, 'I' only from the query (a gap in the reference) and 'D' only from
# the reference (a gap in the query) -- SAM's convention, which is what
# `cigar.decode` prints and what the reference's own gapped strings show.
_CIGAR_RUN = re.compile(r"([0-9]+)([=XIDM])")
# What each op consumes: (first sequence, second sequence).
_CONSUMES = {"=": (1, 1), "X": (1, 1), "M": (1, 1), "I": (1, 0), "D": (0, 1)}

_MODES = ("global", "local")

_END_GAP_REFUSAL = (
    "parasail cannot express this pair of gap prices: a terminal gap costs "
    "(end_gap_score={end_open}, extend_end_gap_score={end_extend}) while an "
    "interior gap costs (open_gap_score={open}, extend_gap_score={extend}).  "
    "Its affine kernels are `nw`, which prices terminal gaps exactly like "
    "interior ones -- set end_gap_score = open_gap_score and "
    "extend_end_gap_score = extend_gap_score -- and `sg`, which makes both "
    "ends free -- set both to 0, the reference's own way of asking for a "
    "semi-global alignment.  A scheme between the two has no parasail "
    "counterpart and Bio.Align.PairwiseAligner is the tool for it."
)


def available() -> bool:
    """Whether the optional accelerator can be imported.

    A caller with a scheme to align and two libraries to choose from can ask
    this instead of catching an ``ImportError``; the readers in this package
    never need it and never ask.
    """
    try:
        import parasail  # noqa: F401
    except ImportError:
        return False
    return True


def _parasail():
    """Import ``parasail`` or explain why it cannot be imported.

    The message names the install command rather than the import: the
    interesting case is a caller who installed biofasting from a wheel, has
    never had parasail, and is looking at a ``ModuleNotFoundError`` from inside
    a library they did not know wanted one -- and `pip install parasail` is not
    the only option they think of, so the extra is named.
    """
    try:
        import parasail
    except ImportError as exc:  # pragma: no cover - exercised by the ImportError test
        raise ImportError(
            "biofasting.alignment needs the optional parasail accelerator, "
            "which is deliberately not a dependency of this package.  Install "
            "it with `pip install biofasting[alignment]` (or `pip install "
            "parasail`) to align with the SIMD kernels; Bio.Align.PairwiseAligner "
            "remains the reference implementation for every scheme here."
        ) from exc
    return parasail


# --------------------------------------------------------------------------
# Scoring schemes, and the refusals that keep them honest.
# --------------------------------------------------------------------------


def _integer(value, what: str) -> int:
    """An int, or a float that is exactly one; anything else refused.

    parasail scores in integer registers of a chosen width, so a fractional
    score has no representation at all -- refusing it is the only honest
    answer, and the refusal names the reference because the reference is what
    accepts it.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{what} must be a number, not {type(value).__name__}")
    if float(value) != int(value):
        raise NotImplementedError(
            f"parasail aligns in integer registers, so {what}={value!r} cannot "
            "be expressed; Bio.Align.PairwiseAligner accepts fractional scores "
            "and is the implementation for that scheme"
        )
    return int(value)


def _gap_cost(value, what: str) -> int:
    """A gap score as the non-negative magnitude parasail subtracts."""
    score = _integer(value, what)
    if score > 0:
        raise NotImplementedError(
            f"parasail cannot express a positive gap score ({what}={value!r}), "
            "which the reference reads as a reward for extending an alignment, "
            "not a penalty"
        )
    return -score


def _as_text(value, what: str) -> str:
    """The sequence as text, one character per byte, refusing anything else.

    The reference formats its argument with ``str()`` and reads it character by
    character, and this package's convention for a measurement is that a byte
    buffer is its Latin-1 characters.  What cannot be accepted is a character
    above U+00FF: parasail's alphabets are C strings indexed by byte, so such a
    letter would be split rather than refused.
    """
    if isinstance(value, (bytes, bytearray, memoryview)):
        text = bytes(value).decode("latin-1")
    elif isinstance(value, str):
        text = value
    else:
        raise TypeError(
            f"{what} must be a str or a bytes-like object, not "
            f"{type(value).__name__}"
        )
    # Encoded rather than scanned: `encode` is one C pass and raises on the
    # first character it cannot carry, where a Python loop over every character
    # would cost more than the kernel on a pair of real reads (measured: 15 us
    # of the 25 us a 150x150 call took before this).
    try:
        text.encode("latin-1")
    except UnicodeEncodeError:
        raise ValueError(
            f"{what} contains a character above U+00FF, which parasail's "
            "byte-indexed alphabets cannot address"
        ) from None
    return text


def _one_byte(value, what: str) -> str:
    """One ASCII letter of an alphabet, checked before it becomes a C string."""
    if not isinstance(value, str) or len(value) != 1:
        raise ValueError(f"{what} must be single characters, not {value!r}")
    if ord(value) > 0x7F:
        raise ValueError(
            f"{what} contains {value!r}, above ASCII; parasail's alphabet is a "
            "C string and cannot carry it"
        )
    return value


def _alphabet(first: str, second: str) -> str:
    """The letters the two sequences actually use, as a sorted alphabet.

    A match/mismatch scheme's alphabet is exactly this and no more, because a
    letter the table does not have would be scored by parasail's unknown slot --
    silently zero against everything.  A letter that is not the same character
    is a mismatch, including a lower-case letter against its upper-case twin,
    because that is what the reference does (measured: ``align("acgt", "ACGT")``
    with match 2 / mismatch -2 scores -8, four mismatches, not four matches).
    """
    alphabet = "".join(sorted(set(first + second)))
    for letter in alphabet:
        _one_byte(letter, "the sequence")
    return alphabet


def _read_table(substitution_matrix):
    """A substitution matrix as an explicit table, refusing what it cannot give.

    The object must be an ``Array`` -- this package's or Biopython's, which are
    both arrays with an ``alphabet`` -- because a matrix that cannot state its
    own alphabet cannot say which letters are unknown, and an unknown letter is
    not a mismatch in this scheme: the reference raises for it, and parasail
    would quietly score it zero against everything.

    Sequence-independent on purpose, so that it can run once in the constructor:
    turning BLOSUM62 into 576 checked integers is 184 us of Python, which is six
    times the kernel it feeds when it is done on every call.
    """
    alphabet = getattr(substitution_matrix, "alphabet", None)
    if alphabet is None or not hasattr(substitution_matrix, "shape"):
        raise TypeError(
            "substitution_matrix must be an Array with an `alphabet` and a "
            f"square shape, not {type(substitution_matrix).__name__}"
        )
    if not isinstance(alphabet, str):
        alphabet = "".join(alphabet)
    for letter in alphabet:
        _one_byte(letter, "the substitution matrix's alphabet")
    size = len(alphabet)
    if tuple(substitution_matrix.shape) != (size, size):
        raise ValueError(
            f"the substitution matrix is {tuple(substitution_matrix.shape)} "
            f"for an alphabet of {size} letters"
        )
    flat = list(substitution_matrix.flat)
    if len(flat) != size * size:
        raise ValueError(
            f"the substitution matrix has {len(flat)} values for an alphabet of "
            f"{size} letters; expected {size * size}"
        )
    # Row-major over the alphabet, the same order `matrix.set_value(i, j, ...)`
    # writes in `_matrix_for` -- one table, read twice, never two tables that
    # happen to agree.  `_integer` does the conversion so that a fractional
    # value is refused rather than truncated to an integer.
    values = tuple(
        _integer(value, f"the substitution matrix's value for index {index}")
        for index, value in enumerate(flat)
    )
    return alphabet, values, size, max(abs(value) for value in values)

    if alphabet is None or not hasattr(substitution_matrix, "shape"):
        raise TypeError(
            "substitution_matrix must be an Array with an `alphabet` and a "
            f"square shape, not {type(substitution_matrix).__name__}"
        )
    if not isinstance(alphabet, str):
        alphabet = "".join(alphabet)
    for letter in alphabet:
        _one_byte(letter, "the substitution matrix's alphabet")
    size = len(alphabet)
    if tuple(substitution_matrix.shape) != (size, size):
        raise ValueError(
            f"the substitution matrix is {tuple(substitution_matrix.shape)} "
            f"for an alphabet of {size} letters"
        )
    flat = list(substitution_matrix.flat)
    if len(flat) != size * size:
        raise ValueError(
            f"the substitution matrix has {len(flat)} values for an alphabet of "
            f"{size} letters; expected {size * size}"
        )
    # Row-major over the alphabet, the same order `matrix.set_value(i, j, ...)`
    # writes in `_matrix_for` -- one table, read twice, never two tables that
    # happen to agree.  `_integer` does the conversion so that a fractional
    # value is refused rather than truncated to an integer.
    values = tuple(
        _integer(value, f"the substitution matrix's value at index {index}")
        for index, value in enumerate(flat)
    )
    return alphabet, values, size, max(abs(value) for value in values)


@lru_cache(maxsize=64)
def _matrix_for(alphabet: str, values: tuple):
    """A parasail matrix for ``(alphabet, values)``, built once.

    Built by hand rather than with ``matrix_create``'s match/mismatch shortcut,
    so that the table parasail scores with and the table this module recomputes
    a path's score with are the same forty numbers instead of two constructions
    that agree today.  ``case_sensitive=True`` is not decoration: it is what
    keeps ``N`` from being read as ``n``.
    """
    parasail = _parasail()
    size = len(alphabet)
    matrix = parasail.matrix_create(alphabet, 0, 0, case_sensitive=True)
    for i in range(size):
        for j in range(size):
            matrix.set_value(i, j, values[i * size + j])
    return matrix


def _width_for(bound: int) -> int:
    """The narrowest exact parasail width a score bound fits in.

    An over-wide width is a factor of two in speed; an under-wide one is a
    wrong answer, which is why this takes a bound and not a guess.  ``_8`` is a
    full sixteen-lane int8 register and ``_16`` eight lanes, which is the whole
    reason the sized bindings are ten to twenty-eight times the unsuffixed
    dispatch.
    """
    for bits in (8, 16, 32):
        if bound <= 2 ** (bits - 1) - 1:
            return bits
    return 64


@dataclass(frozen=True)
class _Scheme:
    """Everything one call needs: which kernel, which table, which prices."""

    prefix: str  # 'nw', 'sg' or 'sw'
    alphabet: str
    values: tuple
    size: int
    open_gap: int  # magnitudes, as parasail subtracts them
    extend_gap: int
    end_open: int
    end_extend: int
    bits: int
    free_ends: bool = False

    def index(self) -> dict:
        return {letter: i for i, letter in enumerate(self.alphabet)}


# --------------------------------------------------------------------------
# The traceback: gapped strings, coordinates, and a score recomputed by hand.
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Counts:
    """What an alignment contains, mirroring the reference's flat names.

    ``identities``, ``mismatches`` and ``gaps`` are the three the reference
    fills in for a pure-Python ``PairwiseAligner``, and they are compared
    against it.  The other three are this module's own arithmetic over the
    aligned letters and the matrix -- the reference leaves them ``None`` for
    these alignments -- and exist because ``substitution_score + gap_score ==
    score`` is the cheapest independent check that the path the SIMD kernel
    traced is the path its own score describes.
    """

    identities: int
    mismatches: int
    gaps: int
    substitution_score: int
    gap_score: int
    score: int


def _decode(cigar_text: str, lengths: tuple, begins: tuple) -> list:
    """Parse the decoded CIGAR, refusing anything this module cannot read.

    ``decode_op``/``decode_len`` are not used: in py-parasail 1.3.4 they are
    static decoders for an *encoded* element, so calling them with an index
    returns a cycling table of nonsense (measured: a rebuild from them
    disagreed with ``cigar.decode`` on 100 of 100 alignments).  ``decode`` is
    the binding's own decoded string, and it is checked against its own runs
    before anything is built from it.
    """
    runs = [(op, int(count)) for count, op in _CIGAR_RUN.findall(cigar_text)]
    if "".join(f"{count}{op}" for op, count in runs) != cigar_text:
        raise RuntimeError(f"unreadable CIGAR from parasail: {cigar_text!r}")
    consumed = [0, 0]
    for op, count in runs:
        first, second = _CONSUMES[op]
        consumed[0] += first * count
        consumed[1] += second * count
    for sequence in (0, 1):
        if begins[sequence] + consumed[sequence] > lengths[sequence]:
            raise RuntimeError(
                f"the CIGAR {cigar_text!r} walks off sequence {sequence}: "
                f"{begins[sequence]} + {consumed[sequence]} > {lengths[sequence]}"
            )
    return runs


def _walk(runs: list, first: str, second: str, begins: tuple, scheme: _Scheme,
          local: bool):
    """One left-to-right pass: the two gapped strings, the coordinates, the gaps.

    A local alignment's optimal segment cannot begin or end with a gap -- a gap
    there would only lose score -- so the leading and trailing single-sequence
    runs parasail reports for `sw` are the parts outside the selected region,
    and the reference's own local alignments drop them.  Trimming them is what
    makes ``align("TTTACGTACGTGGG", "ACGTACGT")`` come back as the eight
    aligned bases and not as fifteen columns.

    The gap half of the counts is finished here, because a price depends on
    whether a run is terminal and only this loop knows that; the other half is
    over the *columns*, which is one Python step per aligned letter, so it lives
    in `_count` and waits until `counts()` asks for it.
    """
    if local:
        while runs and runs[0][0] in "ID":
            begins = (
                begins[0] + runs[0][1] * _CONSUMES[runs[0][0]][0],
                begins[1] + runs[0][1] * _CONSUMES[runs[0][0]][1],
            )
            runs = runs[1:]
        while runs and runs[-1][0] in "ID":
            runs = runs[:-1]

    positions = list(begins)
    coordinate_first: list = []
    coordinate_second: list = []
    gapped_first: list = []
    gapped_second: list = []
    gaps = gap_score = 0

    for number, (op, count) in enumerate(runs):
        terminal = number == 0 or number == len(runs) - 1
        taken_first, taken_second = _CONSUMES[op]
        coordinate_first += [positions[0], positions[0] + taken_first * count]
        coordinate_second += [positions[1], positions[1] + taken_second * count]

        if op in "=XM":
            gapped_first.append(first[positions[0]:positions[0] + count])
            gapped_second.append(second[positions[1]:positions[1] + count])
        else:
            gaps += count
            gap_score -= (
                (scheme.end_open if terminal else scheme.open_gap)
                + (count - 1) * (scheme.end_extend if terminal else scheme.extend_gap)
            )
            if op == "I":
                gapped_first.append(first[positions[0]:positions[0] + count])
                gapped_second.append("-" * count)
            else:
                gapped_first.append("-" * count)
                gapped_second.append(second[positions[1]:positions[1] + count])

        positions[0] += taken_first * count
        positions[1] += taken_second * count

    return (
        "".join(gapped_first),
        "".join(gapped_second),
        coordinate_first,
        coordinate_second,
        gaps,
        gap_score,
    )


def _count(gapped_first: str, gapped_second: str, scheme: _Scheme, gaps: int,
           gap_score: int) -> Counts:
    """What the aligned columns contain, recomputed from the gapped pair.

    One Python step per column, which is why it is the caller's choice and not
    the alignment's: on a 150-base pair this is most of what `align` costs and
    none of what it was asked for.  The gap half arrives already priced, because
    only `_walk` knows which runs were terminal.
    """
    index = scheme.index()
    size = scheme.size
    values = scheme.values
    identities = mismatches = substitution_score = 0
    for letter_first, letter_second in zip(gapped_first, gapped_second):
        if letter_first == "-" or letter_second == "-":
            continue
        if letter_first == letter_second:
            identities += 1
        else:
            mismatches += 1
        substitution_score += values[index[letter_first] * size + index[letter_second]]
    return Counts(
        identities=identities,
        mismatches=mismatches,
        gaps=gaps,
        substitution_score=substitution_score,
        gap_score=gap_score,
        score=substitution_score + gap_score,
    )


def _no_column_alignment(first: str, second: str, scheme: _Scheme) -> Alignment:
    """The alignment with no substitution column, built by hand.

    With every end free and nothing worth substituting, the optimal alignment
    is the whole of one sequence laid against the whole of the other -- an
    insertion run next to a deletion run, so that *no* column pairs two
    letters.  Measured: ``'AAAA'`` against ``'TTTT'`` scores 0 that way and -2
    with any column, and parasail's own answer there is -2, because its `sg`
    kernel cannot trace a path with no diagonal step at all.  The score is 0
    either way once the ``max(·, 0)`` correction is applied, so this function
    supplies the missing *path*.

    A gap priced at zero is still a gap in the coordinates: the reference
    reports this same alignment with coordinates ``[[0, n, n], [0, 0, m]]``,
    which is what `_walk` builds from these two runs.  Both are terminal, so
    both are priced at the (zero) end gap scores, and `_count` finds no column
    with two letters in it -- the whole alignment through the ordinary path
    instead of a second implementation of it.
    """
    gapped_first, gapped_second, coordinate_first, coordinate_second, gaps, gap_score = _walk(
        [("I", len(first)), ("D", len(second))], first, second, (0, 0), scheme, local=False
    )
    return Alignment(
        0, first, second, gapped_first, gapped_second,
        coordinate_first, coordinate_second,
        lambda: _count(gapped_first, gapped_second, scheme, gaps, gap_score),
    )


# --------------------------------------------------------------------------
# The result object, and the aligner.
# --------------------------------------------------------------------------


class Alignment:
    """One optimal alignment: its score, its path, and what the path contains.

    Deliberately not an iterator and not a search: the reference's ``align``
    yields every optimal alignment, and this returns one.  What replaces that
    richness is :meth:`counts`, whose ``score`` recomputes the path's score
    from the aligned letters, so a caller can tell an optimal path from a
    plausible one without a second library.
    """

    __slots__ = ("_first", "_second", "_coordinate_first", "_coordinate_second",
                 "_counts", "_counting", "_gapped_first", "_gapped_second", "score")

    def __init__(self, score, first, second, gapped_first, gapped_second,
                 coordinate_first, coordinate_second, counting):
        self.score = score
        self._first = first
        self._second = second
        self._gapped_first = gapped_first
        self._gapped_second = gapped_second
        self._coordinate_first = coordinate_first
        self._coordinate_second = coordinate_second
        # Not the counts but how to get them: see `counts`.
        self._counting = counting
        self._counts = None

    @property
    def sequences(self) -> tuple:
        """The two sequences as they were given, without gaps."""
        return (self._first, self._second)

    def __getitem__(self, index: int) -> str:
        """The gapped sequence, as the reference's ``alignment[i]`` returns it."""
        if index == 0:
            return self._gapped_first
        if index == 1:
            return self._gapped_second
        raise IndexError(f"index out of range: {index}")

    def __len__(self) -> int:
        return 2

    @property
    def coordinates(self):
        """This alignment as numpy coordinates, one segment per CIGAR run.

        The same alignment Biopython describes, in the finest segmentation --
        where the reference merges a gap into a neighbouring segment, this
        keeps every run separate.  It is a valid ``Alignment.coordinates``
        array, so :meth:`to_biopython` renders it exactly as the reference
        renders its own; it is not byte-identical to the reference's array.
        """
        import numpy

        return numpy.array(
            [self._coordinate_first, self._coordinate_second], dtype=numpy.intp
        )

    def counts(self) -> Counts:
        """What the aligned columns contain, and what they are worth.

        Computed the first time it is asked for and kept afterwards.  A caller
        who wants the path pays for the path: the walk over the columns is one
        Python step per aligned letter -- measured, 24 us of the 46 us a
        150-base pair cost to align in `bench/bench_alignment.py` -- and it is
        not what `align` was asked for.  When it *is* asked for, its ``score``
        is the independent check that the path the kernel traced is the path its
        own score describes.
        """
        counts = self._counts
        if counts is None:
            counts = self._counts = self._counting()
        return counts

    def to_biopython(self):
        """Convert to ``Bio.Align.Alignment``, keeping the reference's surface.

        The conversion is the coordinates this object already holds, so the
        caller gets ``format("clustal")``, ``substitutions()``, ``map()`` and
        the rest from Biopython itself rather than from a reimplementation
        here.  Biopython is imported inside the call, as everywhere else in
        this package.
        """
        try:
            from Bio.Align import Alignment as _ReferenceAlignment
        except ImportError as exc:  # pragma: no cover - exercised by the ImportError test
            raise ImportError(
                "biofasting.alignment needs Biopython to convert to "
                "Bio.Align.Alignment, which is deliberately not a dependency of "
                "this package.  Install it with `pip install biopython`."
            ) from exc
        return _ReferenceAlignment(
            [self._first, self._second], self.coordinates
        )

    def __repr__(self) -> str:
        return (
            f"<Alignment score={self.score} aligned={self._gapped_first!r}/"
            f"{self._gapped_second!r}>"
        )


class Aligner:
    """`Bio.Align.PairwiseAligner`'s scheme, executed on parasail's kernels.

    The attribute names, the scoring formulas and the mode names are the
    reference's, so a scheme written for it can be moved here by changing the
    constructor and nothing else -- and the two defaults that matter are its
    defaults too.  What differs is what is *accepted*: see the module docstring
    for the end-gap rule, which is the one place where the reference expresses
    more than parasail can.

    The scheme is a property of the aligner rather than of a call, so it is
    validated and derived **once, here**, and a caller who gets a price wrong
    hears about it at construction instead of on their first alignment.  It is
    also what makes the wrapper fast enough to be worth having: done per call,
    reading BLOSUM62 into a checked table costs 184 us of Python, six times the
    kernel it feeds (measured, `bench/bench_alignment.py`, before and after).
    """

    def __init__(self, *, mode: str = "global", match_score=1.0, mismatch_score=0.0,
                 open_gap_score=-1.0, extend_gap_score=-1.0, end_gap_score=-1.0,
                 extend_end_gap_score=-1.0, substitution_matrix=None):
        self.mode = mode
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.open_gap_score = open_gap_score
        self.extend_gap_score = extend_gap_score
        self.end_gap_score = end_gap_score
        self.extend_end_gap_score = extend_end_gap_score
        self.substitution_matrix = substitution_matrix

        (
            self._prefix, self._open_gap, self._extend_gap, self._end_open, self._end_extend,
        ) = self._validate()
        if substitution_matrix is None:
            self._table = None
            self._letters = None
            self._match = _integer(match_score, "match_score")
            self._mismatch = _integer(mismatch_score, "mismatch_score")
        else:
            self._table = _read_table(substitution_matrix)
            self._letters = frozenset(self._table[0])
            self._match = self._mismatch = None
        # Per instance, and keyed by the alphabet alone: the numbers are fixed
        # by the constructor, so an alphabet is all a call needs to name its
        # table, and a string hashes far faster than the tuple of values would.
        self._tables: dict[str, tuple] = {}
        self._matrices: dict[str, object] = {}
        self._functions: dict[tuple, object] = {}

    def _validate(self):
        """The mode and the four gap prices, as a kernel choice and magnitudes."""
        if self.mode not in _MODES:
            if self.mode == "fogsaa":
                raise NotImplementedError(
                    "mode='fogsaa' has no parasail counterpart; it is "
                    "Biopython's own algorithm and stays there"
                )
            raise ValueError(
                f"invalid mode {self.mode!r}: expected 'global' or 'local' "
                "(free ends are a global alignment with both end gap scores 0)"
            )

        open_score = _integer(self.open_gap_score, "open_gap_score")
        extend_score = _integer(self.extend_gap_score, "extend_gap_score")
        open_gap = _gap_cost(self.open_gap_score, "open_gap_score")
        extend_gap = _gap_cost(self.extend_gap_score, "extend_gap_score")

        if self.mode == "local":
            # Local scores are insensitive to the end gap scores -- measured,
            # 0 of 63 pairs, for free, interior-priced and the reference's own
            # default -- so `sw` is the kernel whatever they say.
            return "sw", open_gap, extend_gap, open_gap, extend_gap

        end_open_score = _integer(self.end_gap_score, "end_gap_score")
        end_extend_score = _integer(
            self.extend_end_gap_score, "extend_end_gap_score"
        )
        end_open = _gap_cost(self.end_gap_score, "end_gap_score")
        end_extend = _gap_cost(
            self.extend_end_gap_score, "extend_end_gap_score"
        )
        if end_open_score == 0 and end_extend_score == 0:
            prefix = "sg"
        elif end_open_score == open_score and end_extend_score == extend_score:
            prefix = "nw"
        else:
            raise NotImplementedError(
                _END_GAP_REFUSAL.format(
                    end_open=end_open_score, end_extend=end_extend_score,
                    open=open_score, extend=extend_score,
                )
            )
        return prefix, open_gap, extend_gap, end_open, end_extend

    def _match_table(self, alphabet: str):
        """`(values, size, largest)` for a match/mismatch scheme, once per alphabet.

        The same alphabet comes back for every pair of ordinary DNA reads, so
        this is a dictionary hit on all but the first call.
        """
        table = self._tables.get(alphabet)
        if table is None:
            size = len(alphabet)
            values = tuple(
                self._match if i == j else self._mismatch
                for i in range(size)
                for j in range(size)
            )
            table = self._tables[alphabet] = (
                values, size, max(abs(value) for value in values),
            )
        return table

    # -- the scheme, and the two things a call needs from it ---------------

    def _scheme(self, first: str, second: str) -> _Scheme:
        """This pair's kernel, table and width.

        The scheme's own half was settled in the constructor; what is left here
        is the half that depends on the pair -- which letters the alphabet needs,
        and how wide a register its score bound fits in.
        """
        if self._table is None:
            alphabet = _alphabet(first, second)
            values, size, largest = self._match_table(alphabet)
        else:
            alphabet, values, size, largest = self._table
            letters = self._letters
            for sequence, what in (
                (first, "the first sequence"), (second, "the second sequence"),
            ):
                # One C-level pass first and the Python loop only if it found
                # something: an unknown letter is an error, but checking for one
                # letter by letter costs more than the kernel on a short pair.
                unknown = set(sequence) - letters
                if unknown:
                    for letter in sequence:
                        if letter in unknown:
                            raise ValueError(
                                f"{what} contains letters not in the alphabet: {letter!r}"
                            )

        shortest = min(len(first), len(second))
        # An alignment has at most `shortest` columns, each worth at most the
        # matrix's largest value, and the widest gap it can pay for spans the
        # shorter sequence too.  Loose on purpose: one width too wide costs a
        # factor of two, one too narrow is a wrong answer.
        bound = shortest * (largest + self._extend_gap) + self._open_gap
        return _Scheme(
            prefix=self._prefix, alphabet=alphabet, values=values, size=size,
            open_gap=self._open_gap, extend_gap=self._extend_gap,
            end_open=self._end_open, end_extend=self._end_extend,
            bits=_width_for(bound), free_ends=self._prefix == "sg",
        )

    def _matrix(self, scheme: _Scheme):
        """The parasail matrix for a scheme, built once per alphabet.

        Keyed by the alphabet alone because the numbers are fixed by the
        constructor, and memoized per instance so that a repeated alphabet does
        not even hash the tuple of values.
        """
        matrix = self._matrices.get(scheme.alphabet)
        if matrix is None:
            matrix = self._matrices[scheme.alphabet] = _matrix_for(
                scheme.alphabet, scheme.values
            )
        return matrix

    def _function(self, scheme: _Scheme, trace: bool):
        """The sized binding for a scheme and a trace, looked up once."""
        key = (scheme.prefix, scheme.bits, trace)
        function = self._functions.get(key)
        if function is None:
            kind = "trace_scan" if trace else "scan"
            function = self._functions[key] = getattr(
                _parasail(), f"{scheme.prefix}_{kind}_{scheme.bits}"
            )
        return function

    def _run(self, first, second, trace: bool):
        first = _as_text(first, "the first sequence")
        second = _as_text(second, "the second sequence")
        if not first or not second:
            # parasail prints its own diagnostic to the C stderr and returns a
            # NULL result, which py-parasail turns into a ValueError; refusing
            # here gives the reference's message instead of that noise.
            raise ValueError("sequence has zero length")
        # The scheme first, the accelerator second: a scheme this module refuses
        # is refused for its own reason whether or not parasail happens to be
        # installed, so a caller without the extra still gets the useful message.
        scheme = self._scheme(first, second)
        result = self._function(scheme, trace)(
            first, second, scheme.open_gap, scheme.extend_gap, self._matrix(scheme)
        )
        if result.saturated:
            # Unreachable through this class -- the width comes from the score
            # bound -- and kept as the guard that makes it stay unreachable if
            # the bound is ever loosened: parasail saturates silently and the
            # score would be wrong, not an error.
            raise RuntimeError(
                f"parasail saturated its int{scheme.bits} register on "
                f"{len(first)}x{len(second)}; the score bound is too tight"
            )
        return result, scheme

    # -- the two answers ---------------------------------------------------

    def score(self, first, second) -> int:
        """The optimal score, exactly the reference's for a supported scheme.

        With both ends free the answer is ``max(sg, 0)`` -- see the module
        docstring: parasail has no path with zero substitution columns, and with
        all four ends free such a path always scores at least 0.
        """
        result, scheme = self._run(first, second, trace=False)
        return max(result.score, 0) if scheme.free_ends else result.score

    def align(self, first, second) -> Alignment:
        """One optimal alignment, path included.

        The path is optimal and tie-broken by parasail's DP, which is not the
        reference's tie-break: two optimal alignments of the same pair may
        differ in their paths and never in their scores.  What is guaranteed is
        that ``counts().score == score``, so the path is an optimal one.
        """
        first = _as_text(first, "the first sequence")
        second = _as_text(second, "the second sequence")
        result, scheme = self._run(first, second, trace=True)
        if scheme.free_ends and result.score < 0:
            # The optimum is the alignment with no substitution column, which
            # parasail cannot trace; build it here instead of returning the
            # path its own (lower-scoring) answer describes.
            return _no_column_alignment(first, second, scheme)
        cigar = result.cigar
        runs = _decode(
            cigar.decode.decode(),
            (len(first), len(second)),
            (cigar.beg_query, cigar.beg_ref),
        )
        (
            gapped_first, gapped_second, coordinate_first, coordinate_second, gaps, gap_score,
        ) = _walk(
            runs, first, second, (cigar.beg_query, cigar.beg_ref), scheme,
            local=self.mode == "local",
        )
        return Alignment(
            result.score, first, second, gapped_first, gapped_second,
            coordinate_first, coordinate_second,
            lambda: _count(gapped_first, gapped_second, scheme, gaps, gap_score),
        )


def score(first, second, **attributes) -> int:
    """``Aligner(**attributes).score(first, second)``, spelled out."""
    return Aligner(**attributes).score(first, second)


def align(first, second, **attributes) -> Alignment:
    """``Aligner(**attributes).align(first, second)``, spelled out."""
    return Aligner(**attributes).align(first, second)
