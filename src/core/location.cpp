#include "location.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>

#include "line.hpp"

namespace biofasting {
namespace {

// ---------------------------------------------------------------------------
// Python's int()
// ---------------------------------------------------------------------------
//
// The reference's fast path calls int() *before* any regex runs, so its
// permissiveness is part of the grammar rather than an accident: `" 12..34"`
// and `"+12..34"` are locations the reference accepts, and a kernel that only
// took bare digits would refuse two strings the reference reads.  Python also
// accepts underscores between digits -- `int("1_0") == 10`.  Two things Python
// allows are *not* reproduced: Unicode digits, and integers wider than 64 bits.
// A location is ASCII by the format's own definition, and a coordinate that
// does not fit in an int64 is not a coordinate; both are refusals rather than
// silent truncations.

enum class IntParse { ok, not_an_int, overflow };

constexpr std::uint64_t kMaxPositive = 9223372036854775807ULL;  // 2**63 - 1
constexpr std::uint64_t kMaxNegative = 9223372036854775808ULL;  // 2**63

IntParse parse_int(std::string_view text, std::int64_t& out) {
  std::size_t i = 0;
  while (i < text.size() && is_python_space(text[i])) ++i;
  bool negative = false;
  if (i < text.size() && (text[i] == '+' || text[i] == '-')) {
    negative = text[i] == '-';
    ++i;
  }
  const std::uint64_t limit = negative ? kMaxNegative : kMaxPositive;
  std::uint64_t value = 0;
  bool saw_digit = false;
  bool after_digit = false;
  while (i < text.size()) {
    const char c = text[i];
    if (c >= '0' && c <= '9') {
      const std::uint64_t digit = static_cast<std::uint64_t>(c - '0');
      if (value > (limit - digit) / 10) return IntParse::overflow;
      value = value * 10 + digit;
      saw_digit = true;
      after_digit = true;
      ++i;
    } else if (c == '_' && after_digit && i + 1 < text.size() &&
               text[i + 1] >= '0' && text[i + 1] <= '9') {
      // An underscore separates digits and never trails or doubles, which is
      // Python's own rule rather than a guess about what a file meant.
      after_digit = false;
      ++i;
    } else {
      break;
    }
  }
  if (!saw_digit || !after_digit) return IntParse::not_an_int;
  while (i < text.size() && is_python_space(text[i])) ++i;
  if (i != text.size()) return IntParse::not_an_int;
  if (negative) {
    out = value == kMaxNegative ? std::numeric_limits<std::int64_t>::min()
                                : -static_cast<std::int64_t>(value);
  } else {
    out = static_cast<std::int64_t>(value);
  }
  return IntParse::ok;
}

std::size_t digit_end(std::string_view text, std::size_t i) noexcept {
  while (i < text.size() && text[i] >= '0' && text[i] <= '9') ++i;
  return i;
}

// `value + offset`, with the overflow made explicit.  `offset` is 0 or -1, so
// this can overflow only at the very bottom of the range -- and a coordinate
// that wrapped to the other end of the range is a coordinate that is wrong
// while looking like an answer, which is the one failure mode a parser may not
// have.  Python's int has no bottom of the range; this does, and it says so.
bool add_offset(std::int64_t value, int offset, std::int64_t& out) {
  if (offset < 0 && value == std::numeric_limits<std::int64_t>::min()) {
    return false;
  }
  out = value + offset;
  return true;
}

// ---------------------------------------------------------------------------
// The reference's regular expressions, by hand
// ---------------------------------------------------------------------------
//
// There are two families of them and they are not the same, which is not an
// accident of the reference but the reason for two sets of helpers here.
//
// The **category** regex decides what shape a location string has, and it is
//
//   ^(?P<pair>P)|(?P<between>B)|(?P<within>W)|(?P<oneof>O)|(?P<bond>D)|(?P<solo>S)$
//
// -- note that `^` binds to the first alternative alone and `$` to the last,
// and that `re.match` anchors at the start anyway, so every alternative except
// `solo` is a *prefix* match.  The caller then asserts the match covered the
// whole string.  An input that matches a prefix and continues -- `1..2junk`,
// `1..2,3..4`, `123^456x` -- therefore ends in an `AssertionError` raised from
// inside Biopython, which its feature consumer does not catch.  That is a
// refusal and not a parse error, and it is why the category matchers below
// report where they stopped instead of insisting on the whole string.
//
// The **position** regexes (`_re_within_position`, `_re_oneof_position`) are
// prefix matches with no assert anywhere near them, so `Position.fromstring`
// really does accept `(3.9)junk`.  Nothing in a location string can produce
// that, but `parse_position` is as public as the reference's, and a
// differential test of a position is a much sharper instrument than one of a
// whole location.  So it is reproduced rather than tightened.

// `one\-of\(\d+[,\d+]+\)`.  Two things about it are the regex's rather than the
// format's, and both are reproduced because a parser that disagreed about where
// a match ends would disagree about which strings are locations at all.
//
// The character class in the middle admits a `+`, which no well-formed location
// has.  And `\d+[,\d+]+` needs *two* characters, not one -- so `one-of(3)` is
// not a position (the reference asserts that) while `one-of(77)` is, because
// `\d+` may give up its last digit to `[,\d+]+`.  Written out, the run is: at
// least two characters, all from `{0-9 , +}`, the first of them a digit.
//
// Returns the end of the match, or 0 for no match -- nothing can match and end
// at 0.
std::size_t match_one_of(std::string_view text, std::size_t i) noexcept {
  constexpr std::string_view kPrefix = "one-of(";
  if (text.compare(i, kPrefix.size(), kPrefix) != 0) return 0;
  std::size_t j = i + kPrefix.size();
  if (j >= text.size() || text[j] < '0' || text[j] > '9') return 0;
  const std::size_t begin = j;
  while (j < text.size() && (text[j] == ',' || text[j] == '+' ||
                             (text[j] >= '0' && text[j] <= '9'))) {
    ++j;
  }
  if (j - begin < 2 || j >= text.size() || text[j] != ')') return 0;
  return j + 1;
}

// `[<>]?\d+`, the shape `_within_location` allows at either end.
std::size_t match_plain_end(std::string_view text, std::size_t i) noexcept {
  if (i < text.size() && (text[i] == '<' || text[i] == '>')) ++i;
  const std::size_t digits = i;
  i = digit_end(text, i);
  return i == digits ? 0 : i;
}

// `[<>]?(?:\d+|one-of(...))`, the shape `_oneof_location` allows at either end
// and the shape `_split` allows inside a pair.
std::size_t match_choice_end(std::string_view text, std::size_t i) noexcept {
  const std::size_t plain = match_plain_end(text, i);
  return plain != 0 ? plain : match_one_of(text, i);
}

// `_oneof_location` itself: `_choice_end.._choice_end`.  `_pair_location` is
// this with `[<>]?-?\d+` at both ends, and is written out separately below
// because it is the only shape that admits a negative number.
std::size_t match_choice_pair(std::string_view text, std::size_t i) noexcept {
  const std::size_t left = match_choice_end(text, i);
  if (left == 0 || text.compare(left, 2, "..") != 0) return 0;
  const std::size_t right = match_choice_end(text, left + 2);
  return right == 0 ? 0 : right;
}

// `[<>]?-?\d+$`-shaped end used by the pair category: signed.
std::size_t match_signed_end(std::string_view text, std::size_t i) noexcept {
  if (i < text.size() && (text[i] == '<' || text[i] == '>')) ++i;
  if (i < text.size() && text[i] == '-') ++i;
  const std::size_t digits = i;
  i = digit_end(text, i);
  return i == digits ? 0 : i;
}

// `_pair_location`: `[<>]?-?\d+\.\.[<>]?-?\d+`
std::size_t match_signed_pair(std::string_view text, std::size_t i) noexcept {
  const std::size_t left = match_signed_end(text, i);
  if (left == 0 || text.compare(left, 2, "..") != 0) return 0;
  const std::size_t right = match_signed_end(text, left + 2);
  return right == 0 ? 0 : right;
}

// `_reference`: a name and a colon.  The regex's trailing `[a-zA-Z0-9]?` is
// redundant -- the class before it already admits those -- so this consumes the
// name and asks for the colon.
bool is_name_char(char c) noexcept {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
         (c >= '0' && c <= '9') || c == '_' || c == '.' || c == '|';
}

std::size_t match_reference(std::string_view text, std::size_t i) noexcept {
  if (i >= text.size()) return 0;
  const char c = text[i];
  if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))) return 0;
  ++i;
  while (i < text.size() && is_name_char(text[i])) ++i;
  if (i >= text.size() || text[i] != ':') return 0;
  return i + 1;
}

// One captured group of `_split`, matched at `pos`; returns its end, or 0.
//
// `_split` is `re.compile(_any_location).split` with one capturing group and
// the caller takes every other element, so what matters is *where each match
// starts and ends*:
//
//   ( _reference? _oneof_location | complement(_oneof_location)
//     | [^,]+ | complement\([^,]+\) )
//
// The fourth alternative is unreachable -- `[^,]+` matches wherever it could
// start -- and is left out rather than written and never run.  A pair that
// opens a non-comma run is captured as itself and the rest of the run is *also*
// scanned, which is how `1..2junk` yields two parts and not one; the second
// part is what makes the reference refuse the string, so where a match ends is
// not a cosmetic question.
std::size_t match_piece(std::string_view text, std::size_t pos) noexcept {
  const std::size_t pair = match_choice_pair(text, pos);
  if (pair != 0) return pair;
  const std::size_t ref = match_reference(text, pos);
  if (ref != 0) {
    const std::size_t after = match_choice_pair(text, ref);
    if (after != 0) return after;
  }
  constexpr std::string_view kComplement = "complement(";
  if (text.compare(pos, kComplement.size(), kComplement) == 0) {
    const std::size_t inner = match_choice_pair(text, pos + kComplement.size());
    if (inner != 0 && inner < text.size() && text[inner] == ')') {
      return inner + 1;
    }
  }
  if (pos < text.size() && text[pos] != ',') {
    std::size_t j = pos;
    while (j < text.size() && text[j] != ',') ++j;
    return j;
  }
  return 0;
}

// `_split(text)[1::2]`: every captured group, left to right, with the text
// between them discarded.  A comma matches nothing and is skipped.
void split_pieces(std::string_view text, std::vector<std::string_view>& out) {
  std::size_t pos = 0;
  while (pos < text.size()) {
    const std::size_t end = match_piece(text, pos);
    if (end == 0) {
      ++pos;
      continue;
    }
    out.push_back(text.substr(pos, end - pos));
    pos = end;
  }
}

// `\(\d+\.\d+\)` as `_within_location` spells it: no number is read here, so a
// coordinate too wide for an int64 classifies as a within position and is
// refused later, in parse_position, where the reference reads the two numbers.
std::size_t match_within(std::string_view text, std::size_t i) noexcept {
  if (i >= text.size() || text[i] != '(') return 0;
  std::size_t j = i + 1;
  const std::size_t digits = j;
  j = digit_end(text, j);
  if (j == digits || j >= text.size() || text[j] != '.') return 0;
  ++j;
  const std::size_t second = j;
  j = digit_end(text, j);
  if (j == second || j >= text.size() || text[j] != ')') return 0;
  return j + 1;
}

// `([<>]?\d+|\(\d+\.\d+\))\.\.([<>]?\d+|\(\d+\.\d+\))`
std::size_t match_within_location(std::string_view text,
                                  std::size_t i) noexcept {
  const std::size_t plain = match_plain_end(text, i);
  const std::size_t start = plain != 0 ? plain : match_within(text, i);
  if (start == 0 || text.compare(start, 2, "..") != 0) return 0;
  const std::size_t right = match_plain_end(text, start + 2);
  const std::size_t end = right != 0 ? right : match_within(text, start + 2);
  return end == 0 ? 0 : end;
}

// `([<>]?\d+|one-of(...))\.\.([<>]?\d+|one-of(...))`
std::size_t match_one_of_location(std::string_view text,
                                  std::size_t i) noexcept {
  const std::size_t left = match_choice_end(text, i);
  if (left == 0 || text.compare(left, 2, "..") != 0) return 0;
  const std::size_t right = match_choice_end(text, left + 2);
  return right == 0 ? 0 : right;
}

std::size_t match_between(std::string_view text, std::size_t i) noexcept {
  const std::size_t digits = i;
  i = digit_end(text, i);
  if (i == digits || i >= text.size() || text[i] != '^') return 0;
  ++i;
  const std::size_t second = i;
  i = digit_end(text, i);
  return i == second ? 0 : i;
}

// `bond\([<>]?\d+\)`
std::size_t match_bond(std::string_view text, std::size_t i) noexcept {
  constexpr std::string_view kPrefix = "bond(";
  if (text.compare(i, kPrefix.size(), kPrefix) != 0) return 0;
  const std::size_t end = match_plain_end(text, i + kPrefix.size());
  if (end == 0 || end >= text.size() || text[end] != ')') return 0;
  return end + 1;
}

// `[<>]?\d+$`.  Python's `$` also matches just before a trailing newline, and
// that is the one place the anchor is not the same as "the whole string": on
// `"123\n"` the regex matches `123` and the caller's assert then fails.  The
// distinction is kept here so that the two inputs do not come back as one
// status.
std::size_t match_solo(std::string_view text, std::size_t i) noexcept {
  const std::size_t end = match_plain_end(text, i);
  if (end == 0) return 0;
  if (end != text.size() && text[end] != '\n') return 0;
  return end;
}

enum class Category {
  none,
  pair,
  between,
  within,
  one_of,
  bond,
  solo,
};

struct CategoryMatch {
  Category category = Category::none;
  std::size_t end = 0;  // where the match stopped; == size iff it covered it
};

// The first alternative that matches, in the regex's own order, and where it
// stopped.
CategoryMatch match_category(std::string_view text) {
  if (const std::size_t end = match_signed_pair(text, 0)) {
    return {Category::pair, end};
  }
  if (const std::size_t end = match_between(text, 0)) {
    return {Category::between, end};
  }
  if (const std::size_t end = match_within_location(text, 0)) {
    return {Category::within, end};
  }
  if (const std::size_t end = match_one_of_location(text, 0)) {
    return {Category::one_of, end};
  }
  if (const std::size_t end = match_bond(text, 0)) {
    return {Category::bond, end};
  }
  if (const std::size_t end = match_solo(text, 0)) {
    return {Category::solo, end};
  }
  return {};
}

LocationStatus refuse(std::string message, std::string& out) {
  out = std::move(message);
  return LocationStatus::refused;
}

LocationStatus unparseable(std::string_view text, std::string& out) {
  out = "Could not parse feature location '" + std::string(text) + "'";
  return LocationStatus::parser_error;
}

// What one call to the reference's `SimpleLocation.fromstring` produced, which
// is one part except when the origin-wrapping repair turns one written span
// into two.
//
// The warnings are a list here for the same reason they are in
// `LocationResult`: this call can raise more than one, and the compound path
// below concatenates them across parts without reordering them.
struct SimpleLocationResult {
  LocationOp op = LocationOp::simple;
  std::vector<LocationPart> parts;
  std::vector<LocationWarning> warnings;
};

LocationStatus parse_simple_location(std::string_view text, std::int64_t length,
                                     bool has_length, bool circular,
                                     SimpleLocationResult& out,
                                     std::string& message) {
  std::int8_t strand = 0;
  constexpr std::string_view kComplement = "complement(";
  if (text.size() >= kComplement.size() &&
      text.compare(0, kComplement.size(), kComplement) == 0) {
    // `text[11:-1]`, which on a string shorter than that is empty rather than
    // an error -- Python's slicing, reproduced rather than tightened.
    text = text.size() > kComplement.size() + 1
               ? text.substr(kComplement.size(),
                             text.size() - kComplement.size() - 1)
               : std::string_view();
    strand = -1;
  }

  // The simple case first, for speed: `int(s) - 1` and `int(e)` on a string
  // split at `..`.  A third piece, an empty piece or a non-number all fall
  // through to the general case below, which is where the categories are.
  const std::size_t dots = text.find("..");
  if (dots != std::string_view::npos &&
      text.find("..", dots + 2) == std::string_view::npos) {
    std::int64_t start = 0;
    std::int64_t end = 0;
    const IntParse first = parse_int(text.substr(0, dots), start);
    const IntParse second = parse_int(text.substr(dots + 2), end);
    if (first == IntParse::overflow || second == IntParse::overflow) {
      return refuse("feature position out of range in '" + std::string(text) +
                        "'",
                    message);
    }
    if (first == IntParse::ok && second == IntParse::ok) {
      // `int(s) - 1`, guarded: at the bottom of the range the reference wraps to
      // a huge negative, fails `0 <= s`, and falls through to the general case,
      // which refuses it a moment later.  So falling through here is the same
      // route, without the wrapped number in between.
      if (add_offset(start, -1, start) && start >= 0 && start < end) {
        LocationPart part;
        part.start.kind = PositionKind::exact;
        part.start.value = start;
        part.end.kind = PositionKind::exact;
        part.end.value = end;
        part.strand = strand;
        out.parts.push_back(std::move(part));
        return LocationStatus::ok;
      }
    }
  }

  // `ref, text = text.split(":")`: every piece at once, so a second colon is a
  // ValueError and leaves both the ref and the text untouched.  The ref is not
  // validated here -- `123:1..2` is a location with `ref == "123"` -- because
  // the reference does not validate it either.
  std::string ref;
  {
    const std::size_t colon = text.find(':');
    if (colon != std::string_view::npos &&
        text.find(':', colon + 1) == std::string_view::npos) {
      ref.assign(text.substr(0, colon));
      text = text.substr(colon + 1);
    }
  }

  const CategoryMatch match = match_category(text);
  if (match.category == Category::none) return unparseable(text, message);
  if (match.end != text.size()) {
    // The reference matched a prefix and its `assert value == text` failed.
    return refuse("feature location '" + std::string(text) +
                      "' has text after the location",
                  message);
  }

  Position start_pos;
  Position end_pos;
  if (match.category == Category::bond || match.category == Category::solo) {
    if (match.category == Category::bond) {
      // One warning per `bond(...)` part, and the wording never quotes it.
      out.warnings.push_back({LocationWarningKind::bond, std::string()});
      text = text.substr(5, text.size() - 6);
    }
    const LocationStatus first = parse_position(text, -1, start_pos, message);
    if (first != LocationStatus::ok) return first;
    const LocationStatus second = parse_position(text, 0, end_pos, message);
    if (second != LocationStatus::ok) return second;
  } else if (match.category == Category::between) {
    // A between location like "67^68" is a special case: it has zero length,
    // which in Python slice notation is 67:67.  On a circular genome of length
    // N, "N^1" is the junction at the origin.  `s == length` is compared
    // against the declared size even when it is None -- `0 == None` is false in
    // Python and `0 == 0` is true, so a record declaring a zero length accepts
    // `0^1` and a record declaring none does not.
    const std::size_t caret = text.find('^');
    std::int64_t s = 0;
    std::int64_t e = 0;
    if (parse_int(text.substr(0, caret), s) != IntParse::ok ||
        parse_int(text.substr(caret + 1), e) != IntParse::ok) {
      return refuse("feature position out of range in '" + std::string(text) +
                        "'",
                    message);
    }
    if (!(s != std::numeric_limits<std::int64_t>::max() && s + 1 == e) &&
        !(has_length && s == length && e == 1)) {
      message = "invalid feature location '" + std::string(text) + "'";
      return LocationStatus::parser_error;
    }
    start_pos.kind = PositionKind::exact;
    start_pos.value = s;
    end_pos = start_pos;
  } else {
    const std::size_t at = text.find("..");
    if (at == std::string_view::npos ||
        text.find("..", at + 2) != std::string_view::npos) {
      return unparseable(text, message);
    }
    const LocationStatus first =
        parse_position(text.substr(0, at), -1, start_pos, message);
    if (first != LocationStatus::ok) return first;
    const LocationStatus second =
        parse_position(text.substr(at + 2), 0, end_pos, message);
    if (second != LocationStatus::ok) return second;

    // A feature that spans the origin and was written in one piece.  The
    // *declared* size is what the reference compares against, and its absence
    // is a different fact from a declared zero: `if e_pos <= s_pos and length`
    // skips a record with no LOCUS size and takes one that says zero.
    if (end_pos.value <= start_pos.value && has_length && length != 0 &&
        start_pos.value < length) {
      if (!circular) {
        message = "it appears that '" + std::string(text) +
                  "' is a feature that spans the origin, but the sequence "
                  "topology is undefined";
        return LocationStatus::parser_error;
      }
      // The second half is built as `SimpleLocation(0, e_pos)`, and that raises
      // a plain ValueError for a negative end -- which the consumer does not
      // catch, so in the reference this file kills the parse.  The first half,
      // `SimpleLocation(s_pos, length)`, cannot raise: `s_pos < length` is
      // already known, and SimpleLocation checks only start > end.
      if (end_pos.value < 0) {
        return refuse("End location (" + std::to_string(end_pos.value) +
                          ") must be greater than or equal to start location "
                          "(0)",
                      message);
      }
      // The warning quotes the part's own text, with any `ref:` prefix already
      // taken off -- which is the string the reference's `SimpleLocation`
      // warns with, since it splits the reference name off before it matches.
      out.warnings.push_back(
          {LocationWarningKind::origin_wrap, std::string(text)});
      LocationPart to_end;
      to_end.start = start_pos;
      to_end.end.kind = PositionKind::exact;
      to_end.end.value = length;
      to_end.strand = strand;
      LocationPart from_start;
      from_start.start.kind = PositionKind::exact;
      from_start.start.value = 0;
      from_start.end = end_pos;
      from_start.strand = strand;
      if (strand == -1) {
        out.parts.push_back(std::move(from_start));
        out.parts.push_back(std::move(to_end));
      } else {
        out.parts.push_back(std::move(to_end));
        out.parts.push_back(std::move(from_start));
      }
      out.op = LocationOp::join;
      return LocationStatus::ok;
    }
  }

  if (start_pos.value < 0) {
    message = "negative starting position in feature location '" +
              std::string(text) + "'";
    return LocationStatus::parser_error;
  }
  // `SimpleLocation.__init__` refuses an end before its start, and it raises a
  // plain ValueError rather than a parser error.  The consumer catches only the
  // latter, so in the reference this is a crash; here it is a refusal.
  if (start_pos.value > end_pos.value) {
    return refuse("End location (" + std::to_string(end_pos.value) +
                      ") must be greater than or equal to start location (" +
                      std::to_string(start_pos.value) + ")",
                  message);
  }
  LocationPart part;
  part.start = start_pos;
  part.end = end_pos;
  part.strand = strand;
  part.ref = std::move(ref);
  out.parts.push_back(std::move(part));
  return LocationStatus::ok;
}

}  // namespace

LocationStatus parse_position(std::string_view text, int offset, Position& out,
                              std::string& message) {
  if (offset != 0 && offset != -1) {
    message =
        "To convert one-based indices to zero-based indices, offset must be "
        "either 0 (for end positions) or -1 (for start positions).";
    return LocationStatus::refused;
  }
  if (text == "?") {
    // The one position with no integer value: the reference's UnknownPosition
    // does not subclass int, so `int()` of it and every comparison against it
    // raise.  No location string can reach this -- the category regex has no
    // `?` in it -- so it is here for the same reason the reference has it: the
    // function answers for a position, not only for a whole location.
    out = Position();
    out.kind = PositionKind::unknown;
    return LocationStatus::ok;
  }
  if (!text.empty() && (text[0] == '?' || text[0] == '<' || text[0] == '>')) {
    std::int64_t value = 0;
    const IntParse parsed = parse_int(text.substr(1), value);
    if (parsed == IntParse::overflow) {
      return refuse("feature position out of range in '" + std::string(text) +
                        "'",
                    message);
    }
    if (parsed != IntParse::ok) {
      return refuse("Could not parse feature position '" + std::string(text) +
                        "'",
                    message);
    }
    out = Position();
    out.kind = text[0] == '?'   ? PositionKind::uncertain
               : text[0] == '<' ? PositionKind::before
                                : PositionKind::after;
    if (!add_offset(value, offset, out.value)) {
      return refuse("feature position out of range in '" + std::string(text) +
                        "'",
                    message);
    }
    return LocationStatus::ok;
  }
  // `_re_within_position.match` and `_re_oneof_position.match` are prefix
  // matches with no assert, so `(3.9)junk` is a within position here as it is
  // in the reference.  Only a location string could not produce it.
  if (const std::size_t end = match_within(text, 0)) {
    const std::string_view inner = text.substr(1, end - 2);
    const std::size_t dot = inner.find('.');
    std::int64_t left = 0;
    std::int64_t right = 0;
    if (parse_int(inner.substr(0, dot), left) != IntParse::ok ||
        parse_int(inner.substr(dot + 1), right) != IntParse::ok) {
      return refuse("feature position out of range in '" + std::string(text) +
                        "'",
                    message);
    }
    out = Position();
    out.kind = PositionKind::within;
    if (!add_offset(left, offset, out.left) ||
        !add_offset(right, offset, out.right)) {
      return refuse("feature position out of range in '" + std::string(text) +
                        "'",
                    message);
    }
    // The default a caller gets from `int(position)`: the left edge for a start
    // and the right edge for an end, which is what makes `(9.10)..(20.25)` act
    // like a location with overall start 8 and end 25 rather than the reverse.
    out.value = offset == -1 ? out.left : out.right;
    return LocationStatus::ok;
  }
  if (const std::size_t end = match_one_of(text, 0)) {
    // `_re_oneof_position` captures the class between `one-of(` and the final
    // `)`, which is what gets split on commas.
    const std::string_view klass = text.substr(7, end - 8);
    std::vector<std::int64_t> choices;
    std::size_t begin = 0;
    while (true) {
      const std::size_t comma = klass.find(',', begin);
      const std::string_view piece =
          comma == std::string_view::npos
              ? klass.substr(begin)
              : klass.substr(begin, comma - begin);
      std::int64_t value = 0;
      const IntParse parsed = parse_int(piece, value);
      // The class admits a `+`, so `one-of(1+2)` reaches here; the reference
      // calls int() on it and the ValueError escapes uncaught.
      if (parsed != IntParse::ok) {
        return refuse("Could not parse feature position '" + std::string(text) +
                          "'",
                      message);
      }
      choices.push_back(value);
      if (comma == std::string_view::npos) break;
      begin = comma + 1;
    }
    out = Position();
    out.kind = PositionKind::one_of;
    out.choices.reserve(choices.size());
    std::int64_t lowest = 0;
    std::int64_t highest = 0;
    for (std::size_t i = 0; i < choices.size(); ++i) {
      std::int64_t value = 0;
      if (!add_offset(choices[i], offset, value)) {
        return refuse("feature position out of range in '" + std::string(text) +
                          "'",
                      message);
      }
      out.choices.push_back(value);
      if (i == 0 || value < lowest) lowest = value;
      if (i == 0 || value > highest) highest = value;
    }
    // A start position takes the lowest of the choices and an end the highest,
    // which is the same rule as the within's edges and has the same reason: a
    // location whose bounds are uncertain is as wide as it could be.
    out.value = offset == -1 ? lowest : highest;
    return LocationStatus::ok;
  }
  std::int64_t value = 0;
  if (parse_int(text, value) != IntParse::ok) {
    return refuse("Could not parse feature position '" + std::string(text) + "'",
                  message);
  }
  out = Position();
  out.kind = PositionKind::exact;
  if (!add_offset(value, offset, out.value)) {
    return refuse("feature position out of range in '" + std::string(text) + "'",
                  message);
  }
  return LocationStatus::ok;
}

LocationResult parse_location(std::string_view text, std::int64_t length,
                              bool has_length, bool circular, bool stranded) {
  LocationResult result;
  std::int8_t strand = stranded ? 1 : 0;
  std::string_view body = text;
  constexpr std::string_view kComplement = "complement(";

  if (body.size() >= kComplement.size() &&
      body.compare(0, kComplement.size(), kComplement) == 0) {
    if (body.back() != ')') {
      result.status = LocationStatus::refused;
      result.message = "closing bracket missing in '" + std::string(text) + "'";
      return result;
    }
    body =
        body.substr(kComplement.size(), body.size() - kComplement.size() - 1);
    strand = -1;
  }

  const bool compound = body.compare(0, 5, "join(") == 0 ||
                        body.compare(0, 6, "order(") == 0 ||
                        body.compare(0, 5, "bond(") == 0;
  if (!compound) {
    SimpleLocationResult simple;
    const LocationStatus status = parse_simple_location(
        body, length, has_length, circular, simple, result.message);
    // Merged before the status is checked, not after: a part that warns and
    // *then* fails has already warned in the reference, and a warning dropped on
    // the failing path is a warning a caller counting them would miss.
    result.warnings = std::move(simple.warnings);
    if (status != LocationStatus::ok) {
      result.status = status;
      return result;
    }
    // The reference assigns the outer strand to every part and then reverses
    // the list -- a no-op for one part and *not* a no-op for the origin-wrapped
    // pair it has just built: `complement(100..50)` on a 200 bp circle comes
    // back with the 0..50 part first.
    for (LocationPart& part : simple.parts) part.strand = strand;
    if (strand == -1 && simple.parts.size() > 1) {
      std::reverse(simple.parts.begin(), simple.parts.end());
    }
    result.location.op = simple.op;
    result.location.parts = std::move(simple.parts);
    return result;
  }

  LocationOp op = LocationOp::simple;
  std::string_view inner;
  if (body.compare(0, 5, "join(") == 0) {
    op = LocationOp::join;
    inner = body.substr(5, body.size() - 6);
  } else if (body.compare(0, 6, "order(") == 0) {
    op = LocationOp::order;
    inner = body.substr(6, body.size() - 7);
  } else {
    op = LocationOp::bond;
    inner = body.substr(5, body.size() - 6);
  }

  std::vector<std::string_view> pieces;
  split_pieces(inner, pieces);

  std::vector<LocationPart> parts;
  for (const std::string_view piece : pieces) {
    SimpleLocationResult parsed;
    const LocationStatus status = parse_simple_location(
        piece, length, has_length, circular, parsed, result.message);
    // Concatenated, not merged after the fact: the reference warns part by part
    // as it walks them, so the order across parts is the order of the list.
    result.warnings.insert(result.warnings.end(),
                           std::make_move_iterator(parsed.warnings.begin()),
                           std::make_move_iterator(parsed.warnings.end()));
    if (status != LocationStatus::ok) {
      result.status = status;
      return result;
    }
    // A part complemented inside a complemented parent is written the same way
    // twice, and the reference refuses it rather than picking a reading.  Its
    // message is the literal `double complement in '{text}'?` -- the f-string
    // there is missing its `f`, so the braces are part of the text a user sees.
    if (parsed.parts.front().strand == -1) {
      if (strand == -1) {
        result.status = LocationStatus::parser_error;
        result.message = "double complement in '{text}'?";
        return result;
      }
    } else {
      for (LocationPart& part : parsed.parts) part.strand = strand;
    }
    for (LocationPart& part : parsed.parts) parts.push_back(std::move(part));
  }

  if (parts.empty()) {
    // `CompoundLocation` refuses fewer than two parts, and a join that came
    // apart into none -- `join()`, `join(,)` -- reaches that refusal in the
    // reference as a plain ValueError.
    result.status = LocationStatus::refused;
    result.message =
        "CompoundLocation should have at least 2 parts, not [] in '" +
        std::string(text) + "'";
    return result;
  }
  if (parts.size() == 1) {
    // One part is not a compound location, whichever operator introduced it:
    // `join(1..2)` is a SimpleLocation, and the reference returns the part
    // object itself rather than building anything.
    result.location.op = LocationOp::simple;
    result.location.parts = std::move(parts);
    return result;
  }
  if (strand == -1) {
    // GenBank writes a complemented join backwards, so the whole thing comes
    // back in the other order.  Every part is -1 here, which is what the
    // reference asserts before it reverses them -- the assert cannot fail,
    // because a part that was -1 already took the double-complement branch.
    std::reverse(parts.begin(), parts.end());
  }
  result.location.op = op;
  result.location.parts = std::move(parts);
  return result;
}

}  // namespace biofasting
