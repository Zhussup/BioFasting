"""Feature values, and the conversion to Biopython's own types.

The reader one level down -- ``_core.FlatFileIndex.features(id)`` -- is a
kernel: it hands back tuples of numbers and strings and refuses the records it
cannot reproduce.  This module is what turns those tuples into values a caller
can hold (:class:`Position`, :class:`Part`, :class:`Location`,
:class:`Qualifier`, :class:`Feature`) and, for a caller who wants Biopython's
types, :func:`to_seqfeature` builds a ``Bio.SeqFeature.SeqFeature`` out of one.

**Why a location is not a pair of integers.**  ``complement`` carries a strand,
``<``/``>`` make an end fuzzy, ``(3.9)`` is a boundary known only to lie between
two bases, ``one-of(...)`` is a choice, and ``^`` is a zero-length junction
between two bases rather than a span.  A port that flattened those into
``(start, end)`` would lose them silently, so :class:`Position` keeps the kind
and every field the kind has a use for.

**Where the warnings are emitted, and why it is here.**  A kernel that warned
would be a kernel deciding what Python does with it.  The reader therefore
*reports* what it found -- the origin wrap, one entry per repaired part and
carrying the part's text; the dropped ``bond``, one entry per part; the NCBI
escaping, one flag per qualifier and carrying the value the warning is about --
and this module says the words, in the reference's own order and its own
wording.  The suite runs with ``filterwarnings = ["error"]``, so a warning
escaping the kernel would already be a failure rather than a surprise.

The order is not an accident: the reference parses the location first and feeds
the qualifiers after it, so a feature with both kinds of warning warns about its
location before it warns about a qualifier, and a compound location whose second
part wraps the origin warns about the bond in its first part first.

Biopython is **not** a dependency of this package: ``Bio`` is imported inside
:func:`to_seqfeature` and nowhere else, so importing this module imports
nothing, and a caller who only wants our values never pays for a ``SeqFeature``.
"""

from __future__ import annotations

import warnings
from typing import NamedTuple

# The order of `PositionKind` in src/core/location.hpp.  Named rather than
# numeric because the kernel hands back an index and a reader of this file
# should not have to open a C++ header to find out what 3 means.
POSITION_KINDS = (
    "exact",
    "before",
    "after",
    "within",
    "one_of",
    "uncertain",
    "unknown",
)


class Position(NamedTuple):
    """One end of a location, with the kind of fuzziness it carries.

    ``value`` is what the reference's ``int(position)`` returns, which is also
    what every comparison between positions uses -- ``BeforePosition(4) ==
    ExactPosition(4)`` is true in Biopython, and comparing kinds instead of
    values would disagree about whether a feature wraps the origin.  ``left``
    and ``right`` are the two edges of a ``within`` boundary, ``choices`` the
    alternatives of a ``one_of``, and the other kinds use ``value`` alone.
    """

    kind: str
    value: int
    left: int = 0
    right: int = 0
    choices: tuple[int, ...] = ()


class Part(NamedTuple):
    """One span of a location: a start, an end, a strand and a reference name.

    ``strand`` is -1, +1, or 0 for the reference's ``None`` -- a location on a
    protein record has no strand at all, which is a different statement from
    "stranded, strand unknown", and 0 is how the kernel spells that ``None``.
    """

    start: Position
    end: Position
    strand: int
    ref: str


class Location(NamedTuple):
    """A location: one part, or several under an operator.

    ``operator`` is ``"simple"`` for a location that is not compound --
    including ``join(1..2)``, which the reference returns as a plain
    ``SimpleLocation`` because one part is not a compound of anything.
    """

    operator: str
    parts: tuple[Part, ...]


class Qualifier(NamedTuple):
    """One ``/key=value``, with the NCBI escaping warning attached to it.

    ``has_value`` false is the valueless form (``/pseudo``), which the reference
    stores as ``""`` -- and stores only once, however many times the key is
    written, so a repeat is not in the list at all.  ``escape_text`` is the
    value as it stood *before* the doubled quotes were undone: it is what the
    reference's warning quotes, so it travels with the flag rather than being
    recoverable from ``value``.
    """

    key: str
    value: str
    has_value: bool
    escape_warning: bool
    escape_text: str


class Feature(NamedTuple):
    """One feature of a FEATURES table.

    ``status`` is the location's status, and it is not decoration:
    ``"parser_error"`` is the reference catching its own ``LocationParserError``
    and carrying on with ``location = None`` plus a warning, so it is a
    behaviour to reproduce and not a failure -- ``location`` is ``None`` and the
    feature keeps its type and its qualifiers.
    """

    type: str
    location: Location | None
    status: str
    message: str
    warnings: tuple[tuple[str, str], ...]
    qualifiers: tuple[Qualifier, ...]

    def qualifiers_as_dict(self) -> dict[str, list[str]]:
        """The qualifiers as the reference stores them: key -> list of values.

        The reference's dict semantics are already in the list -- a repeated
        valueless key is present once -- so this groups and does not decide.
        """
        out: dict[str, list[str]] = {}
        for qualifier in self.qualifiers:
            out.setdefault(qualifier.key, []).append(qualifier.value)
        return out


def _position(raw) -> Position:
    kind, value, left, right, choices = raw
    return Position(POSITION_KINDS[kind], value, left, right, tuple(choices))


def _part(raw) -> Part:
    start, end, strand, ref = raw
    return Part(_position(start), _position(end), strand, ref)


def _location(raw) -> Location:
    operator, parts = raw
    return Location(operator, tuple(_part(part) for part in parts))


def _qualifier(raw) -> Qualifier:
    return Qualifier(*raw)


def _feature(raw) -> Feature:
    feature_type, location, status, message, warnings_, qualifiers = raw
    return Feature(
        feature_type,
        _location(location) if location is not None else None,
        status,
        message,
        tuple((kind, text) for kind, text in warnings_),
        tuple(_qualifier(qualifier) for qualifier in qualifiers),
    )


def read_features(index, id: str) -> list[Feature]:
    """The features of ``index[id]``, as :class:`Feature` values.

    Raises :class:`ValueError` when the reader declined the record's table --
    the message names the line that made it decline -- and returns an empty list
    for a record with no feature block, which is an answer and not a refusal.
    A feature whose *location* the reference's parser rejects is not a refusal:
    it comes back with ``status == "parser_error"`` and ``location is None``,
    which is what the reference does with it.

    Named like ``read_fasta`` and ``read_genbank``: the reader is a kernel one
    level down, and this is the Python shape a caller holds.
    """
    ok, message, table = index.features(id)
    if not ok:
        raise ValueError(f"{id}: {message}")
    return [_feature(feature) for feature in table]


# --------------------------------------------------------------------------
# Conversion to `Bio.SeqFeature`
# --------------------------------------------------------------------------

# The reference's own wording, character for character -- including the `%r`,
# which is why the text travels with the warning rather than being reconstructed
# from the location: `repr` of the part as the reference saw it is the argument.
_ORIGIN_WRAP = (
    "Attempting to fix invalid location %r as it looks like incorrect origin "
    "wrapping. Please fix input file, this could have unintended behavior."
)
_BOND = "Dropping bond qualifier in feature location"
_NCBI_ESCAPE = (
    'The NCBI states double-quote characters like " should be escaped as "" '
    "(two double - quotes), but here it was not: %r"
)


def _emit_location_warnings(entries) -> None:
    """Say the words for the warns the reader reported, in the order it found them."""
    from Bio import BiopythonParserWarning

    for kind, text in entries:
        if kind == "origin_wrap":
            warnings.warn(_ORIGIN_WRAP % text, BiopythonParserWarning)
        elif kind == "bond":
            warnings.warn(_BOND, BiopythonParserWarning)
        else:  # pragma: no cover - the kernel emits exactly these two kinds
            raise ValueError(f"unknown location warning kind {kind!r}")


def to_location(location: Location):
    """Build the ``Bio.SeqFeature`` location for one of ours.

    The result is a ``SimpleLocation`` or a ``CompoundLocation`` exactly as the
    reference's own parser would have built it -- including ``join(1..2)``,
    which is simple because one part is not a compound.  This function emits no
    warnings; :func:`to_seqfeature` does, from the entries the reader reported.
    """
    from Bio.SeqFeature import (
        AfterPosition,
        BeforePosition,
        CompoundLocation,
        ExactPosition,
        OneOfPosition,
        SimpleLocation,
        UncertainPosition,
        UnknownPosition,
        WithinPosition,
    )

    builders = {
        "exact": ExactPosition,
        "before": BeforePosition,
        "after": AfterPosition,
        "uncertain": UncertainPosition,
        "unknown": UnknownPosition,
    }

    def build(position: Position):
        if position.kind == "within":
            return WithinPosition(position.value, position.left, position.right)
        if position.kind == "one_of":
            return OneOfPosition(
                position.value, [ExactPosition(choice) for choice in position.choices]
            )
        return builders[position.kind](position.value)

    def build_part(part: Part) -> SimpleLocation:
        # The kernel spells the reference's `None` strand as 0, because `None`
        # is not a number a C++ struct can hold; this is where it goes back.
        strand = part.strand if part.strand != 0 else None
        return SimpleLocation(
            build(part.start), build(part.end), strand, ref=part.ref or None
        )

    parts = [build_part(part) for part in location.parts]
    if location.operator == "simple":
        return parts[0]
    return CompoundLocation(parts, operator=location.operator)


def to_seqfeature(feature: Feature):
    """Build the ``Bio.SeqFeature.SeqFeature`` for one of ours.

    The warnings the reference would have emitted come out here, in its own
    order and its own words: the location's in-parse warnings first, then the
    ``parser_error`` note if the location had to be dropped, then one warning
    per qualifier whose NCBI escaping is not standard.
    """
    from Bio import BiopythonParserWarning
    from Bio.SeqFeature import SeqFeature

    _emit_location_warnings(feature.warnings)

    location = None
    if feature.location is not None:
        location = to_location(feature.location)
    elif feature.status == "parser_error":
        # The reference's own consumer catches LocationParserError, warns, and
        # goes on with the feature -- which is a behaviour, so it is reproduced
        # rather than turned into an error here.
        warnings.warn(
            f"{feature.message}; setting feature location to None.",
            BiopythonParserWarning,
        )

    qualifiers: dict[str, list[str]] = {}
    for qualifier in feature.qualifiers:
        if qualifier.escape_warning:
            warnings.warn(_NCBI_ESCAPE % qualifier.escape_text, BiopythonParserWarning)
        qualifiers.setdefault(qualifier.key, []).append(qualifier.value)

    return SeqFeature(location, type=feature.type, qualifiers=qualifiers)
