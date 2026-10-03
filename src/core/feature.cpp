// The FEATURES table: what the reference's scanner and feature consumer make of
// a feature block.  The grammar is described in feature.hpp; what belongs here
// is the transcription of it, in the reference's own order.
//
// GenBank and EMBL are read by one scanner class in the reference
// (`Bio/GenBank/Scanner.py`), differing in four constants: the header width, the
// feature-key column, the continuation spacer and the end markers.  They are
// named below and nowhere else, so that the two formats cannot drift apart by
// accident.

#include "feature.hpp"

#include <algorithm>
#include <cstring>
#include <limits>

#include "line.hpp"

namespace biofasting {
namespace {

// The column the feature key ends at and a location or a qualifier begins: the
// reference's `FEATURE_QUALIFIER_INDENT`, the same number for both formats.
constexpr std::size_t kQualifierIndent = 21;

// A continuation line: 21 spaces in GenBank, `FT` plus the 19 spaces that fill
// the same column in EMBL.
constexpr char kGenBankSpacer[kQualifierIndent] = {
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};
constexpr char kEmblSpacer[kQualifierIndent] = {
    'F', 'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
    ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

// The line that ends the block.  GenBank's list is the reference's own
// `SEQUENCE_HEADERS`; EMBL's end marker is `XX` and its sequence headers are
// `SQ`/`CO`.  Declared in the header as well, so that the scan that records the
// span and this reader cannot disagree about where one ends.
bool block_ends(FlatFileFormat format, const char* text, std::size_t len) {
  if (format == FlatFileFormat::genbank) {
    const std::string_view prefix(text,
                                  rstripped_size(text, std::min<std::size_t>(len, 12)));
    return prefix == "CONTIG" || prefix == "ORIGIN" || prefix == "BASE COUNT" ||
           prefix == "WGS" || prefix == "TSA" || prefix == "TLS";
  }
  const std::string_view prefix(text,
                                rstripped_size(text, std::min<std::size_t>(len, 5)));
  if (prefix == "SQ" || prefix == "CO") return true;
  return rstripped_size(text, len) == 2 && text[0] == 'X' && text[1] == 'X';
}

bool is_start_marker(FlatFileFormat format, std::string_view line) {
  if (format == FlatFileFormat::genbank) {
    return line == "FEATURES             Location/Qualifiers" ||
           line == "FEATURES";
  }
  return line == "FH   Key             Location/Qualifiers" || line == "FH";
}

// Python's `s[a:b]`, which is empty rather than an error when `a` is past the
// end or `b` is before `a`.  Every column cut in this file goes through it,
// because the reference's cuts are index expressions on a string and a short
// line is a case it simply gets a shorter string for.
std::string_view slice(std::string_view text, std::size_t from,
                       std::size_t to) noexcept {
  if (from >= text.size() || to <= from) return std::string_view();
  return text.substr(from, std::min(to, text.size()) - from);
}

std::string_view strip(std::string_view text) noexcept {
  std::size_t begin = 0;
  std::size_t end = text.size();
  while (begin < end && is_python_space(text[begin])) ++begin;
  while (end > begin && is_python_space(text[end - 1])) --end;
  return text.substr(begin, end - begin);
}

std::string_view lstrip(std::string_view text) noexcept {
  std::size_t begin = 0;
  while (begin < text.size() && is_python_space(text[begin])) ++begin;
  return text.substr(begin);
}

// `"".join(text.split())`, which is the reference's `_clean_location` and its
// `FeatureValueCleaner._clean_translation` alike: every run of Python whitespace
// is a separator and the join adds nothing back, so the result has none left.
// It reads like a squash and is not one -- the join string is empty.
std::string remove_whitespace(std::string_view text) {
  std::string out;
  out.reserve(text.size());
  for (const char c : text) {
    if (!is_python_space(c)) out.push_back(c);
  }
  return out;
}

void replace_all(std::string& text, std::string_view from, std::string_view to) {
  if (from.empty()) return;
  std::string out;
  out.reserve(text.size());
  std::size_t at = 0;
  while (true) {
    const std::size_t found = text.find(from, at);
    if (found == std::string::npos) break;
    out.append(text, at, found - at);
    out.append(to);
    at = found + from.size();
  }
  out.append(text, at, text.size() - at);
  text.swap(out);
}

// The reference's escaping test, transcribed alternative by alternative:
// `re.search(r'[^"]"[^"]|^"[^"]|[^"]"$', value)`.  Every alternative is a
// double quote that is not part of a doubled pair, so one pass over the quotes
// answers all three.
bool has_unescaped_quote(std::string_view value) noexcept {
  const std::size_t size = value.size();
  for (std::size_t i = 0; i < size; ++i) {
    if (value[i] != '"') continue;
    if (i + 1 < size && value[i + 1] != '"' && (i == 0 || value[i - 1] != '"')) {
      return true;  // `[^"]"[^"]` or `^"[^"]`
    }
    if (i + 1 == size && i > 0 && value[i - 1] != '"') {
      return true;  // `[^"]"$`
    }
  }
  return false;
}

// The reference's `int(content)` for the sizes both formats hand it, which are
// decimal digits and nothing else.  Python's `int` also takes a sign, an
// underscore between digits and every Unicode digit; a size written that way is
// one this reader does not reproduce, and the caller declines rather than
// guess at it.
bool decimal_int(std::string_view text, std::int64_t& out) noexcept {
  if (text.empty()) return false;
  std::int64_t value = 0;
  for (const char c : text) {
    if (c < '0' || c > '9') return false;
    if (value > (std::numeric_limits<std::int64_t>::max() - (c - '0')) / 10) {
      return false;
    }
    value = value * 10 + (c - '0');
  }
  out = value;
  return true;
}

// One feature's qualifiers, before the consumer's value handling: the key, the
// value as the scanner assembled it (quotes still on, `\n` between the lines of
// a wrapped value), and whether there was a value at all.
struct RawQualifier {
  std::string key;
  std::string value;
  bool has_value = true;
};

class TableReader {
 public:
  TableReader(FlatFileFormat format, const char* data, std::size_t size,
              std::size_t begin, std::size_t end, std::int64_t declared_size,
              bool has_declared_size, bool circular, bool stranded)
      : format_(format),
        data_(data),
        size_(size),
        limit_(end),
        at_(begin),
        declared_size_(declared_size),
        has_declared_size_(has_declared_size),
        circular_(circular),
        stranded_(stranded) {}

  FeatureTable run() {
    FeatureTable table;
    // The reference's marker skip: `while self.line.rstrip() in
    // FEATURE_START_MARKERS`.  GenBank's block begins at its `FEATURES` line and
    // EMBL's at its first `FH`, which is the line `parse_header` stopped on.
    while (at_ < limit_) {
      const Line line = read_line(data_, size_, at_);
      if (!line.exists) break;
      if (!is_start_marker(format_, rstrip_view(line.text, line.len))) break;
      at_ = line.next;
    }

    while (at_ < limit_) {
      const Line line = read_line(data_, size_, at_);
      if (!line.exists) break;
      if (block_ends(format_, line.text, line.len)) break;
      const std::string_view text(line.text,
                                  rstripped_size(line.text, line.len));

      // A line with nothing in the key column is nothing -- the reference skips
      // it before it looks at its length, which is why a blank line between two
      // features is a blank line and not a short one.
      if (strip(slice(text, 2, kQualifierIndent)).empty()) {
        at_ = line.next;
        continue;
      }
      if (line.len < kQualifierIndent) {
        table.refused = true;
        table.offset = static_cast<std::size_t>(line.text - data_);
        table.message = "feature table, byte offset " + std::to_string(table.offset) +
                        ": a line too short to hold a feature key and a location "
                        "(the reference warns and skips it)";
        return table;
      }
      // A space anywhere past the key column is what the reference reads as a
      // file that over-indented its location: it warns, then re-splits the line
      // on whitespace and takes the second token as the location.  This reader
      // declines instead, because a repair reproduced differently is a location
      // that looks right and is not.
      if (line.text[kQualifierIndent] != ' ' &&
          std::memchr(line.text + kQualifierIndent, ' ',
                      line.len - kQualifierIndent) != nullptr) {
        table.refused = true;
        table.offset = static_cast<std::size_t>(line.text - data_);
        table.message = "feature table, byte offset " + std::to_string(table.offset) +
                        ": the location column contains a space, which the "
                        "reference repairs by re-splitting the line";
        return table;
      }

      const std::size_t feature_offset = static_cast<std::size_t>(line.text - data_);
      const std::string key(strip(slice(text, 2, kQualifierIndent)));

      // The feature's own lines: the location's line as it stands, then every
      // continuation, each reduced to its column-21 text and stripped -- except
      // the first, which the reference leaves as it found it.
      std::vector<std::string> lines;
      lines.emplace_back(slice(text, kQualifierIndent, text.size()));
      std::size_t cursor = line.next;
      while (cursor < limit_) {
        const Line next = read_line(data_, size_, cursor);
        if (!next.exists) break;
        const bool continuation =
            next.len >= kQualifierIndent &&
            std::memcmp(next.text, spacer(), kQualifierIndent) == 0;
        // A blank line in the midst of a feature is a continuation of nothing,
        // which the reference drops when it filters the empty lines out.
        const bool blank = next.len != 0 && rstripped_size(next.text, next.len) == 0;
        if (!continuation && !blank) break;
        lines.emplace_back(strip(slice(std::string_view(next.text, next.len),
                                       kQualifierIndent, next.len)));
        cursor = next.next;
      }
      at_ = cursor;

      if (!read_feature(key, feature_offset, lines, table)) return table;
    }
    return table;
  }

 private:
  const char* spacer() const noexcept {
    return format_ == FlatFileFormat::genbank ? kGenBankSpacer : kEmblSpacer;
  }

  bool refuse(FeatureTable& table, std::size_t offset, std::string why) {
    table.refused = true;
    table.offset = offset;
    table.message = "feature table, byte offset " + std::to_string(offset) +
                    ": " + std::move(why);
    return false;
  }

  // `parse_feature(key, lines)`, and then the consumer's reading of what it
  // returned.  Returns false when the table has been marked refused.
  bool read_feature(std::string_view key, std::size_t offset,
                    std::vector<std::string>& lines, FeatureTable& table) {
    std::vector<RawQualifier> raw;
    std::size_t index = 0;
    const auto next_line = [&]() -> const std::string* {
      // The reference's `(x for x in lines if x)`: an empty line is not a line.
      while (index < lines.size() && lines[index].empty()) ++index;
      if (index >= lines.size()) return nullptr;
      return &lines[index++];
    };

    const std::string* line = next_line();
    if (line == nullptr) {
      return refuse(table, offset, "the feature has no location line");
    }
    std::string location(*line);
    // A location too long for one line is split after a comma, and it is the
    // comma that says so and not the indentation -- which is why the next line
    // is taken whenever this one ends in one.
    while (!location.empty() && location.back() == ',') {
      line = next_line();
      if (line == nullptr) {
        return refuse(table, offset,
                      "the location ends in a comma that no line follows");
      }
      location += strip(*line);
    }
    if (std::count(location.begin(), location.end(), '(') >
        std::count(location.begin(), location.end(), ')')) {
      return refuse(table, offset,
                    "the location wraps its parentheses without breaking at a "
                    "comma (the reference warns and keeps joining lines)");
    }

    std::size_t line_number = 0;
    std::string last_key;
    while ((line = next_line()) != nullptr) {
      const std::string_view text(*line);
      if (line_number == 0 && text[0] == ')') {
        // The location's closing parenthesis on a line of its own.
        location += strip(text);
      } else if (text[0] == '/') {
        const std::size_t equals = text.find('=');
        RawQualifier qualifier;
        std::string_view value;
        if (equals == std::string_view::npos) {
          qualifier.key.assign(text.substr(1));  // the `/pseudo` form
        } else {
          qualifier.key.assign(text.substr(1, equals - 1));
          value = text.substr(equals + 1);
        }
        if (equals != std::string_view::npos && !value.empty() &&
            value[0] == ' ' && !lstrip(value).empty() && lstrip(value)[0] == '"') {
          return refuse(table, offset,
                        "a qualifier has white space between its '=' and its "
                        "opening quote");
        }
        if (equals == std::string_view::npos) {
          qualifier.has_value = false;
        } else if (value.empty()) {
          // `/note=` is an empty value, which is not the same as no value.
        } else if (value == "\"") {
          qualifier.value = "\"";  // one lone quote is left as it stands
        } else if (value[0] == '"') {
          // A quoted value runs until a line ends with a quote.  The lines are
          // joined with newlines here and with spaces by the consumer, which is
          // why a wrapped `translation` ends up with no whitespace at all.
          qualifier.value.assign(value);
          while (qualifier.value.back() != '"') {
            const std::string* more = next_line();
            if (more == nullptr) {
              return refuse(table, offset,
                            "a quoted qualifier value is never closed");
            }
            qualifier.value += '\n';
            qualifier.value += *more;
          }
        } else {
          qualifier.value.assign(value);
        }
        last_key = qualifier.key;
        raw.push_back(std::move(qualifier));
      } else {
        // An unquoted continuation.  The reference asserts that a qualifier came
        // before it and that its key is the one it remembers; both are
        // AssertionErrors there, which is to say a file it crashes on.
        if (raw.empty() || raw.back().key != last_key) {
          return refuse(table, offset,
                        "a line continues a qualifier that is not the one above "
                        "it");
        }
        if (!raw.back().has_value) {
          return refuse(table, offset,
                        "a continuation line follows a qualifier that has no "
                        "value");
        }
        raw.back().value += '\n';
        raw.back().value += text;
      }
      ++line_number;
    }

    Feature feature;
    feature.type.assign(key);
    feature.location_text = clean_location(location);
    feature.location = parse_location(feature.location_text, declared_size_,
                                      has_declared_size_, circular_, stranded_);
    if (feature.location.status == LocationStatus::refused) {
      // The reference raises here -- a bare AssertionError or a ValueError --
      // and the exception leaves `SeqIO.parse` altogether.
      return refuse(table, offset, feature.location.message);
    }

    // The consumer.  The order is the reference's and each step is load
    // bearing: the escaping warning is about the value as it stands *between*
    // the enclosing quotes coming off and the doubled quotes being undone, and
    // `/key` with no value stores an empty string for a key it has not seen
    // before while dropping the key outright if it has.
    std::vector<std::string> seen;
    for (RawQualifier& qualifier : raw) {
      const bool known =
          std::find(seen.begin(), seen.end(), qualifier.key) != seen.end();
      if (!qualifier.has_value) {
        if (known) continue;
        seen.push_back(qualifier.key);
        FeatureQualifier out;
        out.key = std::move(qualifier.key);
        out.has_value = false;
        feature.qualifiers.push_back(std::move(out));
        continue;
      }
      std::string value = std::move(qualifier.value);
      replace_all(value, "\n", " ");
      if (value.size() > 1 && value.front() == '"' && value.back() == '"') {
        value = value.substr(1, value.size() - 2);
      }
      FeatureQualifier out;
      out.key = std::move(qualifier.key);
      out.escape_warning = has_unescaped_quote(value);
      if (out.escape_warning) out.escape_text = value;
      replace_all(value, "\"\"", "\"");
      if (out.key == "translation") value = remove_whitespace(value);
      out.value = std::move(value);
      seen.push_back(out.key);
      feature.qualifiers.push_back(std::move(out));
    }

    table.features.push_back(std::move(feature));
    return true;
  }

  // `_clean_location`, plus the one repair the reference makes for records that
  // predate a proper location syntax: a `replace(266,"c")` becomes the number it
  // names and the rest of the line goes.  The cut is `location_line[8:comma_pos]`
  // in Python, which with no comma at all is a slice to one before the end --
  // reproduced rather than tidied, because it is what the reference hands to its
  // parser.
  std::string clean_location(std::string_view raw) const {
    std::string cleaned = remove_whitespace(raw);
    if (cleaned.find("replace") == std::string::npos) return cleaned;
    const std::size_t size = cleaned.size();
    const std::size_t comma = cleaned.find(',');
    const std::size_t from = std::min<std::size_t>(8, size);
    std::size_t to = comma == std::string::npos ? (size == 0 ? 0 : size - 1) : comma;
    if (to < from) to = from;
    if (to > size) to = size;
    return cleaned.substr(from, to - from);
  }

  FlatFileFormat format_;
  const char* data_;
  std::size_t size_;
  std::size_t limit_;  // the first byte after the block
  std::size_t at_;
  std::int64_t declared_size_;
  bool has_declared_size_;
  bool circular_;
  bool stranded_;
};

}  // namespace

bool feature_block_ends(FlatFileFormat format, const char* text,
                        std::size_t len) {
  return block_ends(format, text, len);
}

FeatureTable parse_feature_table(FlatFileFormat format, const char* data,
                                 std::size_t size, std::size_t begin,
                                 std::size_t end, std::int64_t declared_size,
                                 bool has_declared_size, bool circular,
                                 bool stranded) {
  if (format == FlatFileFormat::swiss) {
    FeatureTable table;
    table.refused = true;
    table.offset = begin;
    table.message =
        "SwissProt's FT lines are a different grammar from the INSDC feature "
        "table and are not reproduced";
    return table;
  }
  if (end <= begin) return FeatureTable();
  TableReader reader(format, data, size, begin, end, declared_size,
                     has_declared_size, circular, stranded);
  return reader.run();
}

}  // namespace biofasting
