// Feature locations: the grammar of `Bio.SeqFeature.Location.fromstring`.
//
// A feature location is a string like `complement(join(490883..490885,1..879))`
// or `(9.10)..(20.25)` or `1000^1`, and reading it is what a flat-file parse
// actually spends its time on.  The measurement in bench/targets.md (M19) puts
// the reference at 20.2 us per feature, of which 4.2 is building the
// `SeqFeature` object and 16 is reading the location string and the qualifier
// block.  This file is the 16, split in two: the location half is here, the
// qualifier half is in feature.cpp.
//
// **What a location is, and therefore what this returns.**  It is not a pair of
// integers: `complement` carries a strand, `<`/`>` make an end fuzzy,
// `(3.9)` is a boundary known only to lie between two bases, `one-of(...)` is a
// choice, and `^` is a zero-length junction between two bases rather than a
// span.  A port that flattened those into `(start, end)` would lose them
// silently, so the structure below keeps every one of them.
//
// **Where the boundary is.**  The kernel reproduces the reference's *parser*
// and stops there: it returns positions and parts, and building
// `SimpleLocation`/`CompoundLocation` objects out of them is the interop
// layer's job.  That is the same line M20 drew between the flat-file index and
// `SeqRecord`, and for the same reason -- the core does not depend on
// Biopython.
//
// **Where it declines.**  Three outcomes, not two:
//
//   * `ok` -- the location, reproduced;
//   * `parser_error` -- the reference raises `LocationParserError`, which the
//     consumer catches and turns into `location = None` plus a warning.  That
//     is a *behaviour*, so a port has to reproduce it, and this is the flag
//     that says so;
//   * `refused` -- the reference raises something else, or accepts the string
//     while silently discarding text from it.  Two examples, both real: it
//     raises a bare `AssertionError` on `1..2,3..4` (the category regex matches
//     `1..2` and an `assert value == text` fails), and `_split` hands
//     `1..2junk` two parts, of which the second is parsed as a location and
//     fails.  A refusal is a design difference and not a bug: it is this
//     reader saying it will not guess.
//
// The one place the two really do differ on *valid* input is Unicode: Python's
// `int()` accepts every Unicode digit and this parses ASCII only.  A location
// string is ASCII by the format's own definition, so a Devanagari digit in one
// is a file this reader declines -- and it declines it loudly rather than
// reproducing a Python accident.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace biofasting {

// The reference's `Position` classes, one enumerator each.  `exact` is a plain
// integer, `before`/`after`/`uncertain` are that integer with a sign of
// fuzziness, `within` carries both edges of the interval a boundary lies in,
// `one_of` carries the choices, and `unknown` is the `?` that has no integer
// value at all.
enum class PositionKind : std::uint8_t {
  exact = 0,
  before = 1,
  after = 2,
  within = 3,
  one_of = 4,
  uncertain = 5,
  unknown = 6,
};

struct Position {
  PositionKind kind = PositionKind::exact;
  // What the reference's `int(position)` returns, which is also what every
  // comparison between positions uses: `BeforePosition(4) == ExactPosition(4)`
  // is true in Biopython, and a port that compared kinds instead of values
  // would disagree about whether a feature wraps the origin.
  std::int64_t value = 0;
  std::int64_t left = 0;   // within
  std::int64_t right = 0;  // within
  std::vector<std::int64_t> choices;  // one_of
};

enum class LocationOp : std::uint8_t {
  simple = 0,  // one part, not a compound location
  join = 1,
  order = 2,
  bond = 3,
};

struct LocationPart {
  Position start;
  Position end;
  // -1, +1, or 0 for the reference's `None`: a location on a protein has no
  // strand, and 0 would be a different statement -- "stranded, strand unknown".
  std::int8_t strand = 0;
  std::string ref;  // the `AL391218.9:` prefix, empty when there is none
};

struct Location {
  LocationOp op = LocationOp::simple;
  std::vector<LocationPart> parts;  // never empty when the status is `ok`
};

enum class LocationStatus : std::uint8_t {
  ok = 0,
  parser_error = 1,  // the reference raises LocationParserError
  refused = 2,       // the reference raises something else, or drops text
};

// One warning the reference would have emitted from *inside* the parse, in the
// order it would have emitted it.  `origin_wrap` quotes the part's text and
// `bond` quotes nothing (its wording is fixed), so `text` is empty for a bond
// and carries the argument for a wrap.
enum class LocationWarningKind : std::uint8_t {
  origin_wrap = 0,
  bond = 1,
};

struct LocationWarning {
  LocationWarningKind kind = LocationWarningKind::origin_wrap;
  std::string text;
};

struct LocationResult {
  LocationStatus status = LocationStatus::ok;
  Location location;
  std::string message;  // non-empty unless the status is `ok`
  // The warnings the reference emits while parsing, in its own order.  They are
  // reported rather than spelled because the wording is the Python layer's: the
  // kernel says what it found, and the layer decides what to say about it.
  //
  // A *list* and not a flag because the reference emits one warning per
  // offending part -- `join(30..5,60..2)` on a circle repairs two parts and
  // warns twice -- and because the order is observable: a part that is a
  // `bond(...)` warns before the part after it wraps the origin.  A consumer
  // that emitted one warning where the reference emitted two, or the two in the
  // other order, would differ in something a caller can see.
  std::vector<LocationWarning> warnings;
};

// `Location.fromstring(text, length, circular, stranded)`.  `has_length` is
// `length is not None` in Python and is not the same thing as `length != 0`:
// the origin-wrapping repair is skipped for a record with no declared size and
// taken for a record that declares zero, and those are different files.
LocationResult parse_location(std::string_view text, std::int64_t length,
                              bool has_length, bool circular, bool stranded);

// `Position.fromstring(text, offset)`, exposed because the reference exposes it
// and because a differential test of a *position* is a much sharper instrument
// than a differential test of a whole location.  `offset` is 0 for an end
// position and -1 for a start position; anything else is refused.
LocationStatus parse_position(std::string_view text, int offset, Position& out,
                              std::string& message);

}  // namespace biofasting
