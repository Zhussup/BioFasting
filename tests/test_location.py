"""The location kernel against `Bio.SeqFeature.Location.fromstring`.

A feature location is not a pair of integers: `complement` carries a strand,
`<`/`>` make an end fuzzy, `(3.9)` is a boundary known only to lie between two
bases, `one-of(...)` is a choice, and `^` is a zero-length junction.  The kernel
returns that structure, and this file is the gate that says it returns the
*same* structure as the reference -- not a plausible one.

The comparison is not between objects but between canonical forms.  Both sides
are reduced to the same nested tuples of plain integers (`_canonical_reference`
does the reducing for Biopython, `_core.parse_location` returns it directly), so
a difference in the answer cannot hide behind a difference in the spelling of
the answer.  Reading the reference's `WithinPosition` through its private
`_left`/`_right` is deliberate: they are the only names it has for the two edges
of the interval, and `int()` of the same object gives the third value that
matters.

Four outcomes are checked, not two.  The reference has three -- it returns an
object, it raises `LocationParserError`, or it raises something else -- and the
kernel maps them to `ok`, `parser_error` and `refused`.  The mapping is the
point of the `parser_error`/`refused` split: `LocationParserError` is *caught*
by the reference's own feature consumer and turned into a missing location plus
a warning, so a port that treated it as a plain error would lose a recorded
behaviour; anything else kills the parse, and there this reader declines loudly
instead of reproducing the crash.

One divergence is expected and is asserted as such rather than tolerated:
Python's `int` is unbounded and this kernel's is 64 bits, so a coordinate that
does not fit is refused.  A wrapped coordinate would be a wrong coordinate that
looks like an answer, and that is the one outcome a parser may not have.  Those
values are named in `OUT_OF_RANGE` and checked for a refusal, so the day the
kernel starts accepting them the test says so instead of quietly disagreeing.
"""

from __future__ import annotations

import random
import warnings
from pathlib import Path
import importlib.util
import sys

import pytest

from Bio import BiopythonParserWarning
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    Location,
    LocationParserError,
    OneOfPosition,
    Position,
    UncertainPosition,
    UnknownPosition,
    WithinPosition,
)

from biofasting import _core

BENCH = Path(__file__).resolve().parent.parent / "bench"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"_bench_{name}", BENCH / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


seqio_corpus = _load("seqio_corpus")

# The kernel's `PositionKind`, in the order the enum declares it.
EXACT, BEFORE, AFTER, WITHIN, ONE_OF, UNCERTAIN, UNKNOWN = range(7)


# --- the canonical form ------------------------------------------------------


def _canonical_position(position):
    """Biopython's Position as (kind, value, left, right, choices)."""
    if isinstance(position, UnknownPosition):
        return (UNKNOWN, 0, 0, 0, ())
    if isinstance(position, UncertainPosition):
        return (UNCERTAIN, int(position), 0, 0, ())
    if isinstance(position, BeforePosition):
        return (BEFORE, int(position), 0, 0, ())
    if isinstance(position, AfterPosition):
        return (AFTER, int(position), 0, 0, ())
    if isinstance(position, WithinPosition):
        return (WITHIN, int(position), position._left, position._right, ())
    if isinstance(position, OneOfPosition):
        return (
            ONE_OF,
            int(position),
            0,
            0,
            tuple(int(choice) for choice in position.position_choices),
        )
    if isinstance(position, ExactPosition):
        return (EXACT, int(position), 0, 0, ())
    raise AssertionError(f"unknown position type {type(position)!r}")


def _canonical_part(part):
    return (
        _canonical_position(part.start),
        _canonical_position(part.end),
        # `None` is the reference's "no strand at all", which a location on a
        # protein has; the kernel spells it 0, and so does this.
        part.strand if part.strand is not None else 0,
        part.ref or "",
    )


def _canonical_location(location):
    if isinstance(location, CompoundLocation):
        return (location.operator, tuple(_canonical_part(p) for p in location.parts))
    return ("simple", (_canonical_part(location),))


def reference(text, length=None, circular=False, stranded=True):
    """(status, location, warnings) -- the reference's own three outcomes."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            location = Location.fromstring(text, length, circular, stranded)
        except LocationParserError:
            return "parser_error", None, [], [str(w.message) for w in caught]
        except Exception:  # noqa: BLE001 - the taxonomy is the point
            return "refused", None, [], [str(w.message) for w in caught]
        return "ok", _canonical_location(location), [
            str(w.message) for w in caught
        ], []


def kernel(text, length=None, circular=False, stranded=True):
    status, location, _message, warnings = _core.parse_location(
        text, length, circular, stranded
    )
    return status, location, warnings


def compare(text, length=None, circular=False, stranded=True):
    ref_status, ref_location, ref_warnings, _ = reference(
        text, length, circular, stranded
    )
    our_status, our_location, our_warnings = kernel(
        text, length, circular, stranded
    )
    assert our_status == ref_status, (
        f"{text!r} length={length} circular={circular} stranded={stranded}: "
        f"reference said {ref_status}, kernel said {our_status}"
    )
    if ref_status == "ok":
        assert our_location == ref_location, (
            f"{text!r} length={length} circular={circular} stranded={stranded}:\n"
            f"  reference {ref_location}\n  kernel    {our_location}"
        )
        # A warning is something a caller can count, so the count is part of the
        # contract and not just the fact.  The reference emits one warning per
        # offending part and the kernel reports one entry per part, so the
        # totals have to agree.  Only checked for `ok`: on the failing paths the
        # reference's own warning list is not the parsed location's.
        assert len(ref_warnings) == len(our_warnings), (
            f"{text!r} length={length} circular={circular} stranded={stranded}:\n"
            f"  reference warned {ref_warnings}\n"
            f"  kernel    reported {our_warnings}"
        )
    return ref_warnings, our_warnings


# --- the cases ---------------------------------------------------------------

# Every location the generated corpus writes, plus the shapes a writer never
# emits but a file does.  Lengths are varied over every case because the
# declared size is an input to two decisions -- origin wrapping and the `N^1`
# junction -- and a kernel tested at one length has been tested at one length.
LOCATIONS = [
    # the corpus's own shapes
    "1..1000",
    "join(1..100,200..300)",
    "join(1..100,200..300,400..500)",
    "complement(1..100)",
    "complement(join(1..100,200..300))",
    "1..1000",
    "1..100",
    "join(1000..1900,2000..2900)",
    # fuzzy ends
    "<1..100",
    "1..>100",
    "complement(<1..>100)",
    "<1..<100",
    ">1..100",
    # within and one-of
    "(9.10)..(20.25)",
    "1..(20.25)",
    "(9.10)..100",
    "(9.10)..(20.25)",
    "one-of(1,2)..5",
    "5..one-of(1,2)",
    "one-of(1,2)..one-of(3,4)",
    "one-of(5,8,11)",
    "complement(one-of(1,2)..one-of(3,4))",
    # between, including the circular junction
    "123^124",
    "1000^1",
    "1^2",
    "5^6",
    # references
    "AL391218.9:105173..108462",
    "NC_016402.1:6618..6676",
    "join(NC_016402.1:6618..6676,181647..181905)",
    # operators
    "join(1..2,3..4)",
    "order(1..2,3..4)",
    "join(order(1..2,3..4),5..6)",
    "join(1..2)",
    "order(1..2)",
    "bond(196)",
    "join(bond(196),197..200)",
    # a solo position is a zero-length location
    "196",
    "<196",
    ">196",
    "bond(<196)",
    # complements of every shape, and the double complement
    "complement(1..2)",
    "complement(100..50)",
    "complement(join(1..2,3..4))",
    "join(complement(1..2),3..4)",
    "join(complement(1..2),complement(3..4))",
    "complement(complement(1..2))",
    "complement(join(complement(1..2),3..4))",
    # origin wrapping, which needs a declared length and a circular topology
    "100..50",
    "complement(100..50)",
    "join(100..50,200..150)",
    "1000..1",
    # the reference's `int()` leaks into the grammar
    " 12..34",
    "+12..34",
    "1_0..20",
    "5..5",
    "10..10",
    "0..0",
    "0",
    # what a malformed file looks like
    "1..2junk",
    "1..2,3..4",
    "a:1..2,3..4",
    "123^456",
    "123x",
    "123\n",
    "50..10",
    "join()",
    "join(,)",
    "join(1..2,)",
    "order()",
    "bond()",
    "complement(1..2",
    "complement()",
    "one-of(3)",
    "one-of(1+2)..5",
    "one-of(1,2)..",
    "<x..2",
    "?",
    "?5",
    "-5..-2",
    "-5..10",
    "5..-2",
    "a:b:1..2",
    "123:",
    ":1..2",
    "1..2^3",
    "2^4",
    "(9.10)..(20.25)..30",
    "(3.9)..5",
    "5..(3.9)",
    "",
]

# Coordinates Python accepts and a 64-bit kernel does not.  Refused on purpose;
# see the module docstring.
OUT_OF_RANGE = [
    "9223372036854775808",
    "9223372036854775808..9223372036854775809",
    "-9223372036854775808..5",
    "(9223372036854775808.9223372036854775809)..5",
    "one-of(9223372036854775808,1)..5",
]


@pytest.mark.parametrize("text", LOCATIONS)
@pytest.mark.parametrize("length", [None, 0, 5, 1000])
@pytest.mark.parametrize("circular", [False, True])
def test_matches_the_reference(text, length, circular):
    compare(text, length, circular)


@pytest.mark.parametrize("text", ["1..2", "<1..2", "196"])
def test_stranded_false_is_a_strandless_location(text):
    """`stranded=False` is what a protein record asks for, and it is not a
    missing strand: the reference passes `None`, and 0 here means `None`."""
    _, _ = compare(text, 1000, False, False)
    status, location, _ = kernel(text, 1000, False, False)
    assert status == "ok"
    assert all(part[2] == 0 for part in location[1])


@pytest.mark.parametrize("text", OUT_OF_RANGE)
@pytest.mark.parametrize("length", [None, 1000])
def test_out_of_range_coordinates_are_refused(text, length):
    """The one deliberate divergence.  Refusing is the design; what would be a
    bug is a wrapped coordinate, so this asserts the refusal *and* that nothing
    came back at all."""
    status, location, _message, _warnings = _core.parse_location(text, length)
    assert status == "refused", f"{text!r} came back as {status}: {location!r}"
    assert location is None


@pytest.mark.parametrize(
    "text",
    ["9223372036854775808", "-9223372036854775809", "<9223372036854775808"],
)
@pytest.mark.parametrize("offset", [0, -1])
def test_out_of_range_positions_are_refused(text, offset):
    status, position, message = _core.parse_position(text, offset)
    assert status == "refused", f"{text!r} came back as {status}: {position!r}"
    assert position is None and message


def test_the_bottom_of_the_range_does_not_wrap():
    """`-9223372036854775808` is a legal int64 and one below it is not.  Offset
    by -1 the legal one would land one below zero of the range, which in C is
    the other end of it -- a wrong number that looks like an answer."""
    status, position, _message = _core.parse_position("-9223372036854775808", 0)
    assert status == "ok" and position[1] == -(2**63)
    status, position, _message = _core.parse_position("-9223372036854775808", -1)
    assert status == "refused" and position is None


# --- the corpus --------------------------------------------------------------
#
# The generated corpus is what the benchmark measures, so the kernel is held to
# it as well as to the hand-written cases above.  The location strings are read
# out of the file by the format's own column rule rather than taken from the
# generator's objects, because the generator's objects are `SimpleLocation`s and
# what a parser has to read is the text -- and a corpus test that fed the kernel
# objects would have tested the writer twice.

FEATURE_COLUMN = 21  # where a location starts, GenBank and EMBL alike


def _corpus_locations(text: str, fmt: str) -> list[str]:
    """Every feature's location, as the file spells it.

    GenBank: the FEATURES block, a key in column 6 and the location in column
    22.  EMBL: lines beginning `FT   `, with the key occupying columns 6-21 --
    a line whose columns 6-21 are blank continues the qualifiers of the feature
    above it, which is the rule the reference's own consumer uses.
    """
    out: list[str] = []
    in_features = False
    for line in text.splitlines():
        if fmt == "genbank":
            if line.startswith("FEATURES"):
                in_features = True
            elif line.startswith(("ORIGIN", "//", "BASE COUNT", "CONTIG")):
                in_features = False
            elif in_features and line[:5] == "     " and line[5] != " ":
                out.append(line[FEATURE_COLUMN:].rstrip())
        elif line.startswith("FT   ") and line[5:FEATURE_COLUMN].strip():
            out.append(line[FEATURE_COLUMN:].rstrip())
    return out


@pytest.mark.parametrize(
    "name,fmt,count,length,per_kb,annot",
    [row for row in seqio_corpus.CORPORA if row[1] in ("genbank", "embl")],
)
def test_the_corpus_locations_match_the_reference(name, fmt, count, length, per_kb, annot):
    """Each row: extract the locations from the bytes, then hold the kernel to
    the reference on every one of them.

    The extraction is checked before it is used.  Re-parsing each extracted
    string must give the location the reference's own `SeqIO.parse` produced
    from the same file; if the column rule were off by one the two would
    disagree, so a wrong rule fails here rather than quietly testing a
    different set of strings.
    """
    from io import StringIO

    from Bio import SeqIO

    text = seqio_corpus.build(fmt, count, length, per_kb)
    records = list(SeqIO.parse(StringIO(text), fmt))
    extracted = _corpus_locations(text, fmt)
    assert len(extracted) == sum(len(record.features) for record in records)

    checked = 0
    index = 0
    for record in records:
        size = len(record)
        for feature in record.features:
            written = extracted[index]
            index += 1
            assert _canonical_location(Location.fromstring(written, size)) == (
                _canonical_location(feature.location)
            ), f"{name}: {written!r} is not the location of {feature.type}"
            compare(written, size)
            checked += 1
    # A row with `per_kb == 0` has no features at all -- that is what the
    # `-bare` rows are for -- so there is nothing to check there and the
    # count must be zero rather than merely positive.
    assert checked == len(extracted)
    assert checked > 0 if per_kb else checked == 0


# --- fuzz --------------------------------------------------------------------


def test_fuzz_against_the_reference():
    """The cases above are the ones a person thought of.  This is the part that
    is not: strings drawn from a location-ish alphabet, including ones no file
    should contain, compared outcome for outcome.

    The seed is fixed.  A fuzz test whose failures cannot be reproduced is a
    test that reports a mood rather than a bug.
    """
    random.seed(20261003)
    alphabet = list("0123456789.,()<>^:-_ abjoinordercomplementonf") + [
        "join(",
        "order(",
        "bond(",
        "complement(",
        "one-of(",
        "..",
        "^",
    ]
    pieces = [
        "1..2",
        "complement(3..4)",
        "one-of(5,6)..7",
        "8",
        "bond(9)",
        "100..50",
        "<1..>2",
        "(3.9)..(20.25)",
        "",
    ]
    for trial in range(4000):
        kind = trial % 4
        if kind == 0:
            sample = "".join(
                random.choice(alphabet) for _ in range(random.randrange(0, 14))
            )
        elif kind == 1:
            sample = (
                "join("
                + ",".join(
                    random.choice(pieces) for _ in range(random.randrange(0, 4))
                )
                + ")"
            )
        elif kind == 2:
            sample = (
                random.choice(["", "complement(", "join(", "order(", "bond("])
                + "".join(
                    random.choice(alphabet) for _ in range(random.randrange(0, 10))
                )
                + random.choice(["", ")", "))"])
            )
        else:
            a, b = random.randrange(0, 30), random.randrange(0, 30)
            sample = random.choice(["", "complement(", "join(", "order("]) + random.choice(
                [
                    f"{a}..{b}",
                    f"{a}^{b}",
                    str(a),
                    f"(3.{b})..({a}.9)",
                    f"one-of({a},{b})..{a}",
                    f"<{a}..>{b}",
                ]
            ) + random.choice(["", ")", f",{b}..{a})"])
        # The alphabet reaches no coordinate wider than two digits, so no case
        # here can hit the int64 boundary and every one of them has a reference
        # answer.  If that stops being true, `compare` fails loudly.
        compare(
            sample,
            random.choice([None, 0, 5, 1000]),
            random.choice([False, True]),
            random.choice([True, True, False]),
        )


# --- the two warnings, and the two failure statuses --------------------------


def test_warnings_are_reported_as_a_list_of_findings():
    """The two warnings the reference emits from *inside* the parse.  The kernel
    reports what it found and the Python layer is what says the words, so the
    list -- its length, its order, and for the origin wrap the text the
    reference quotes -- is the whole contract here.  One warning per offending
    part, not one per location: a location with two of them is where a
    flag-shaped interface would have silently said `True` where the reference
    said it twice."""
    _status, _location, _message, warnings = _core.parse_location(
        "100..50", 200, True, True
    )
    assert warnings == [("origin_wrap", "100..50")]
    _status, _location, _message, warnings = _core.parse_location("bond(196)", 1000)
    # A top-level `bond(...)` is read as an *operator*, and the operator path
    # does not warn: the warning belongs to `SimpleLocation.fromstring`'s bond
    # category, which a top-level bond never reaches.
    assert warnings == []
    assert reference("bond(196)", 1000)[2] == []
    assert reference("join(bond(196),197..200)", 1000)[2] != []
    _status, _location, _message, warnings = _core.parse_location(
        "join(bond(196),197..200)", 1000
    )
    assert warnings == [("bond", "")]
    _status, _location, _message, warnings = _core.parse_location("1..2", 1000)
    assert warnings == []


def test_the_reference_warns_once_per_part_and_so_do_we():
    """The reason the report is a list and not a flag: on `join(30..5,60..2)`
    over a circle the reference repairs *both* parts and warns twice, and
    `join(bond(1),bond(2))` warns twice more.  A boolean would agree with the
    reference about the kind of warning and disagree about how many times it was
    emitted."""
    assert len(reference("join(30..5,60..2)", 100, True)[2]) == 2
    _status, _location, _message, warnings = _core.parse_location(
        "join(30..5,60..2)", 100, True
    )
    assert warnings == [("origin_wrap", "30..5"), ("origin_wrap", "60..2")]

    assert len(reference("join(bond(1),bond(2))", 30)[2]) == 2
    _status, _location, _message, warnings = _core.parse_location(
        "join(bond(1),bond(2))", 30
    )
    assert warnings == [("bond", ""), ("bond", "")]


def test_the_warnings_of_a_mixed_location_keep_the_references_order():
    """A part that is a `bond(...)` warns where it stands, so a `join` whose
    first part is a bond and whose second wraps the origin has to report the two
    in that order.  Two separate collections would have made the order an
    accident of how the Python layer read them out."""
    assert [str(m)[:6] for m in reference("join(bond(1),30..5)", 100, True)[2]] == [
        "Droppi",
        "Attemp",
    ]
    _status, _location, _message, warnings = _core.parse_location(
        "join(bond(1),30..5)", 100, True
    )
    assert warnings == [("bond", ""), ("origin_wrap", "30..5")]


def test_a_warning_is_not_a_status():
    """`100..50` on a circle is repaired *and* warned about; it is not an
    error, and the repair produces a join whose parts are in the reference's
    order."""
    status, location, _message, warnings = _core.parse_location(
        "100..50", 200, True, True
    )
    assert status == "ok" and warnings == [("origin_wrap", "100..50")]
    assert location[0] == "join"
    (first, second) = location[1]
    assert (first[0][1], first[1][1]) == (99, 200)
    assert (second[0][1], second[1][1]) == (0, 50)


def test_the_wrap_of_a_complemented_span_comes_back_reversed():
    """The reference reverses the parts after assigning the strand, which for
    the origin-wrapped pair is not a no-op: `complement(100..50)` on a circle
    puts the 0..50 part first."""
    status, location, _message, _warnings = _core.parse_location(
        "complement(100..50)", 200, True, True
    )
    assert status == "ok"
    (first, second) = location[1]
    assert (first[0][1], first[1][1], first[2]) == (0, 50, -1)
    assert (second[0][1], second[1][1], second[2]) == (99, 200, -1)


def test_double_complement_reports_the_references_literal_message():
    """The reference's message there is `double complement in '{text}'?` -- its
    f-string is missing the `f`.  A port that "fixed" that would be announcing a
    text the reference's own users never see."""
    status, _location, message, _warnings = _core.parse_location(
        "complement(join(complement(1..2),3..4))", 1000
    )
    assert status == "parser_error"
    assert message == "double complement in '{text}'?"


@pytest.mark.parametrize(
    "text,length,expected",
    [
        # the reference raises a bare AssertionError: its category regex matched
        # a prefix and it asserted the match covered the string
        ("1..2junk", 1000, "refused"),
        ("1..2,3..4", 1000, "refused"),
        ("a:1..2,3..4", 1000, "refused"),
        # a plain ValueError from SimpleLocation, which the consumer does not
        # catch -- and the same string with a declared length is a parse error
        # instead, because the origin-wrapping repair is tried first
        ("50..10", None, "refused"),
        ("50..10", 1000, "parser_error"),
        ("join()", 1000, "refused"),
        ("join(,)", 1000, "refused"),
        ("complement(1..2", 1000, "refused"),
        # LocationParserError proper
        ("123^456", 1000, "parser_error"),
        ("0..0", 1000, "parser_error"),
        ("-5..10", 1000, "parser_error"),
        ("complement(join(complement(1..2),3..4))", 1000, "parser_error"),
    ],
)
def test_parser_error_and_refused_are_different_outcomes(text, length, expected):
    """The split that matters.  A `LocationParserError` is caught by the
    reference's feature consumer and becomes a missing location plus a warning,
    so the kernel reproduces it as `parser_error`; anything else kills the
    parse, and there this reader says no rather than reproducing a crash."""
    assert kernel(text, length)[0] == expected
    assert reference(text, length)[0] == expected
