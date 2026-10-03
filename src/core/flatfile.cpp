#include "flatfile.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include "annotations.hpp"
#include "feature.hpp"
#include "line.hpp"

namespace biofasting {
namespace {

// The keyword column: the width of the left block that carries the key, and
// therefore where a value starts.  This is Biopython's own HEADER_WIDTH per
// scanner, and it is not a style choice -- a header block's continuation lines
// are indented to exactly this column, so it is also the rule that says which
// lines belong to the keyword above them.
constexpr std::size_t kEmblHeaderWidth = 5;
constexpr std::size_t kGenBankHeaderWidth = 12;

// Where a GenBank ORIGIN line's residues begin.  The position number is
// right-aligned in nine columns and one space follows it, so the cut is at ten.
constexpr std::size_t kGenBankResidueColumn = 10;

constexpr char kFiveSpaces[5] = {' ', ' ', ' ', ' ', ' '};

bool is_record_end(const Line& line) noexcept {
  return rstripped_size(line.text, line.len) == 2 && line.text[0] == '/' &&
         line.text[1] == '/';
}

bool all_digits(const char* text, std::size_t size) noexcept {
  if (size == 0) return false;
  for (std::size_t i = 0; i < size; ++i) {
    if (text[i] < '0' || text[i] > '9') return false;
  }
  return true;
}

// The keyword a header line introduces and its value.  Both live in line.hpp
// rather than here, because the annotations reader walks the same lines and two
// answers to "which keyword is this" would be two answers to "which block is
// this".  `strip_view`, which those two are built on, lives there for the same
// reason.

// `while "  " in data: data = data.replace("  ", " ")`, which both the GenBank
// scanner and the reference's EMBL consumer apply to a version line before this
// reader sees the same text.
std::string collapse_spaces(std::string text) {
  std::size_t at = text.find("  ");
  while (at != std::string::npos) {
    text.replace(at, 2, " ");
    at = text.find("  ", at);
  }
  return text;
}

std::size_t digits_to_size(std::string_view text) noexcept {
  std::size_t value = 0;
  for (const char c : text) {
    value = value * 10 + static_cast<std::size_t>(c - '0');
  }
  return value;
}

// `version()`: the one place a `VERSION` line (GenBank) and an `SV` line (EMBL,
// obsolete) are read.  `AB000001.2` is an accession plus a version suffix, and
// anything else becomes the record's id outright.  Shared because the reference
// shares it -- both lines are mapped to the same consumer -- and because the
// version suffix is what `record_end` later hangs off the id.
void apply_version(const std::string& text, FlatFileRecord& record,
                   std::size_t& sequence_version, bool& has_sequence_version) {
  const std::size_t dot = text.find('.');
  if (dot != std::string::npos && text.find('.', dot + 1) == std::string::npos &&
      all_digits(text.data() + dot + 1, text.size() - dot - 1)) {
    if (record.id.empty()) record.id.assign(text.substr(0, dot));
    sequence_version = digits_to_size(text.substr(dot + 1));
    has_sequence_version = true;
  } else if (!text.empty()) {
    record.id.assign(text);
  }
}

// --- the two facts a feature location is read under ----------------------- //
//
// `Location.fromstring` takes three things besides the location string: the
// record's declared size, whether it is circular, and whether it is stranded.
// All three come from the header line, and the reference reads them by column
// with a branch per layout era.  A layout this reader does not reproduce leaves
// `has_header_facts` false, and a record whose features are then asked for
// refuses rather than parsing them under inputs the reference never used.

// Python's `s[a:b]`, which is short rather than an error at the end of a line.
std::string_view slice_view(std::string_view text, std::size_t from,
                            std::size_t to) noexcept {
  if (from >= text.size() || to <= from) return std::string_view();
  return text.substr(from, std::min(to, text.size()) - from);
}

std::string uppercase(std::string_view text) {
  std::string out(text);
  for (char& c : out) {
    if (c >= 'a' && c <= 'z') c = static_cast<char>(c - ('a' - 'A'));
  }
  return out;
}

bool looks_like_protein(std::string_view residue_type) {
  return uppercase(residue_type).find("PROTEIN") != std::string::npos;
}

// `while "  " in s: s = s.replace("  ", " ")` followed by `s.split(" ")`, which
// is how the reference splits a LOCUS line's name from its length.  Written out
// rather than approximated because the two are not the same rule: a leading
// single space survives the collapse and becomes an empty token, and the
// reference then refuses the line for having three fields.
std::vector<std::string> collapse_and_split(std::string_view text) {
  std::string collapsed(text);
  while (true) {
    const std::size_t at = collapsed.find("  ");
    if (at == std::string::npos) break;
    collapsed.replace(at, 2, " ");
  }
  std::vector<std::string> parts;
  std::size_t start = 0;
  while (true) {
    const std::size_t at = collapsed.find(' ', start);
    if (at == std::string::npos) {
      parts.push_back(collapsed.substr(start));
      break;
    }
    parts.push_back(collapsed.substr(start, at - start));
    start = at + 1;
  }
  return parts;
}

// Python's `text.split()`, which skips runs and both ends.
std::vector<std::string> split_whitespace(std::string_view text) {
  std::vector<std::string> parts;
  std::size_t at = 0;
  while (at < text.size()) {
    while (at < text.size() && is_python_space(text[at])) ++at;
    const std::size_t begin = at;
    while (at < text.size() && !is_python_space(text[at])) ++at;
    if (at > begin) parts.emplace_back(text.substr(begin, at - begin));
  }
  return parts;
}

// The reference's `int(content)` for the sizes it is handed.  A sign, an
// underscore between digits and every Unicode digit are shapes Python's `int`
// also takes and this does not, which is a refusal and not a different number.
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

bool is_size_unit(std::string_view text) noexcept {
  return text == " bp " || text == " aa " || text == " rc ";
}

// The LOCUS line's own statement of the record's size, topology and residue
// type, in the two column layouts the reference knows.  Anything else -- the
// EnsEMBL space-separated form, a truncated line, a strand type it warns about
// -- leaves the facts unset, which is what makes a feature table on such a
// record a refusal rather than a guess.
bool genbank_locus_facts(const Line& line, FlatFileRecord& record) {
  const std::string_view text(line.text, line.len);
  const std::string_view units_at_29 = slice_view(text, 29, 33);

  if (is_size_unit(units_at_29) && slice_view(text, 55, 62) == "       ") {
    // The pre-229.0 layout.
    if (slice_view(text, 41, 42) != " ") return false;
    const std::string_view entry = strip_view(slice_view(text, 42, 51));
    if (!(entry.empty() || entry == "linear" || entry == "circular")) {
      return false;
    }
    if (slice_view(text, 51, 52) != " ") return false;
    if (!strip_view(slice_view(text, 62, 73)).empty() &&
        (slice_view(text, 64, 65) != "-" || slice_view(text, 68, 69) != "-")) {
      return false;  // a date the reference warns about and drops
    }
    const std::vector<std::string> name_and_length =
        collapse_and_split(slice_view(text, 12, 29));
    if (name_and_length.size() != 2) return false;
    if (!decimal_int(name_and_length[1], record.declared_size)) return false;
    std::string residue;
    const std::string_view middle = strip_view(slice_view(text, 33, 51));
    if (middle.empty() && units_at_29 == " aa ") {
      residue = "PROTEIN";
    } else {
      residue.assign(middle);
    }
    record.has_declared_size = true;
    record.circular = entry == "circular";
    record.stranded = !looks_like_protein(residue);
    return true;
  }

  const std::string_view units_at_40 = slice_view(text, 40, 44);
  const std::string_view entry = strip_view(slice_view(text, 54, 64));
  if (!is_size_unit(units_at_40) ||
      !(entry.empty() || entry == "linear" || entry == "circular")) {
    return false;
  }
  // The post-229.0 layout.  The reference pads a short line to 79 columns and
  // warns that it was truncated; a padded line is one this reader declines,
  // because the padding is the reference's repair and not the file's content.
  if (line.len < 79) return false;
  const std::string_view strand = slice_view(text, 44, 47);
  if (!(strand == "   " || strand == "ss-" || strand == "ds-" ||
        strand == "ms-")) {
    return false;  // the reference warns and falls back to name and size only
  }
  const std::string molecule = uppercase(strip_view(slice_view(text, 47, 54)));
  if (!(molecule.empty() || molecule.find("DNA") != std::string::npos ||
        molecule.find("RNA") != std::string::npos)) {
    return false;
  }
  if (slice_view(text, 54, 55) != " ") return false;
  if (slice_view(text, 63, 64) != " ") return false;
  if (slice_view(text, 67, 68) != " ") return false;
  if (!strip_view(slice_view(text, 68, 79)).empty() &&
      (slice_view(text, 70, 71) != "-" || slice_view(text, 74, 75) != "-")) {
    return false;
  }
  const std::vector<std::string> name_and_length =
      collapse_and_split(slice_view(text, 12, 40));
  if (name_and_length.size() != 2) return false;
  if (!decimal_int(name_and_length[1], record.declared_size)) return false;
  std::string residue;
  if (strip_view(slice_view(text, 44, 54)).empty() && units_at_40 == " aa ") {
    residue = "PROTEIN " + std::string(slice_view(text, 54, 63));
    residue = std::string(strip_view(residue));
  } else {
    residue.assign(strip_view(slice_view(text, 44, 63)));
  }
  record.has_declared_size = true;
  record.circular = entry == "circular";
  record.stranded = !looks_like_protein(residue);

  // The four values the annotations dict takes from this line, in the reference's
  // own columns -- `line[44:54]` for the molecule type, which starts inside the
  // strand field, and not the `[47:54]` the residue check above looks at.
  record.molecule_type = std::string(strip_view(slice_view(text, 44, 54)));
  record.topology = std::string(entry);
  if (!strip_view(slice_view(text, 64, 76)).empty()) {
    record.data_file_division = std::string(slice_view(text, 64, 67));
  }
  if (!strip_view(slice_view(text, 68, 79)).empty()) {
    record.date = std::string(slice_view(text, 68, 79));
  }
  record.has_annotation_header = true;
  return true;
}

// The EMBL `ID` line's 2006-and-later layout: seven semicolon-separated fields,
// of which the topology, the molecule type and the size are the three a feature
// location needs.  The older and patent layouts are real, are read by the
// reference, and are not reproduced here -- which is a statement of scope and
// the reason a feature table on such a record refuses.
bool embl_id_facts(std::string_view trimmed, FlatFileRecord& record) {
  if (std::count(trimmed.begin(), trimmed.end(), ';') != 6) return false;
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t at = trimmed.find(';', start);
    if (at == std::string_view::npos) {
      fields.emplace_back(strip_view(trimmed.substr(start)));
      break;
    }
    fields.emplace_back(strip_view(trimmed.substr(start, at - start)));
    start = at + 1;
  }
  if (fields.size() != 7) return false;
  // The reference's topology consumer refuses anything but linear or circular,
  // and an empty field is stored as an empty topology -- which is not circular.
  if (!(fields[2].empty() || fields[2] == "linear" || fields[2] == "circular")) {
    return false;
  }
  const std::vector<std::string> length = split_whitespace(fields[6]);
  if (length.size() != 2) return false;
  const std::string unit = uppercase(length[1]);
  if (!(unit == "BP" || unit == "BP." || unit == "AA" || unit == "AA.")) {
    return false;
  }
  if (!decimal_int(length[0], record.declared_size)) return false;
  record.has_declared_size = true;
  record.circular = fields[2] == "circular";
  record.stranded = !looks_like_protein(fields[2] + " " + fields[3]);

  // The annotations dict's share of this line.  There is no `date` here because
  // the EMBL scanner reads no date at all -- `DT` lines reach
  // `_feed_header_lines`, whose consumer table has no entry for them, and are
  // dropped.  That is a property of the reference and not of this reader.
  record.topology = fields[2];
  record.molecule_type = fields[3];
  record.data_file_division = fields[5];
  record.has_annotation_header = true;
  return true;
}

// One residue line, de-columned the way its format says.  `out` may be null, in
// which case the line is only measured -- and that is not an optimisation but
// the reason the scan and the read cannot disagree: the same function decides
// the record's length and its bytes, so a length can never describe a sequence
// this reader would refuse to produce.
//
// Returns false when the line is not shaped the way its format says it is.  The
// reference copes with several of those shapes by warning and guessing -- a
// wrongly indented GenBank line is shifted by one byte and parsed anyway -- and
// a guess reproduced differently is a sequence that looks right and is not.
bool read_residue_line(FlatFileFormat format, const Line& line, std::string* out,
                       std::size_t& count) {
  count = 0;
  std::size_t begin = 0;
  std::size_t end = rstripped_size(line.text, line.len);

  // A blank line inside a sequence block is a warning in GenBank and an error
  // in EMBL.  Either way it is not data, and neither is a line that says the
  // residues live elsewhere.
  if (end == 0) return false;
  if (format == FlatFileFormat::genbank && line.len >= 6 &&
      std::memcmp(line.text, "CONTIG", 6) == 0) {
    return false;
  }

  if (format == FlatFileFormat::genbank) {
    if (line.len <= kGenBankResidueColumn - 1) return true;  // nothing to add
    if (line.text[kGenBankResidueColumn - 1] != ' ') return false;
    begin = kGenBankResidueColumn;
  } else {
    if (format == FlatFileFormat::embl) {
      if (line.len < kEmblHeaderWidth ||
          std::memcmp(line.text, kFiveSpaces, kEmblHeaderWidth) != 0) {
        return false;
      }
    } else if (line.len >= 2 && !is_python_space(line.text[0])) {
      // SwissProt dispatches on a two-character key, and a residue line's key is
      // two spaces.  Anything else here is a line this reader does not know.
      return false;
    }
    begin = std::min(kEmblHeaderWidth, line.len);
    while (begin < end && line.text[begin] == ' ') ++begin;

    if (format == FlatFileFormat::embl) {
      // EMBL right-aligns the residue count on the last line of every SQ block,
      // so the reference drops the line's last whitespace-delimited token when
      // it is all digits.  Leaving it in appends the record's own length to its
      // sequence -- digits are not residues and nothing downstream would notice.
      std::size_t token = end;
      while (token > begin && line.text[token - 1] != ' ') --token;
      if (token > begin && token < end &&
          all_digits(line.text + token, end - token)) {
        end = token;
      }
    }
  }

  // Uppercasing is the formats' rule and not SwissProt's: the two INSDC
  // consumers uppercase what they are handed and the SwissProt reader does not,
  // which is a difference a corpus of uppercase residues cannot show and a real
  // lowercase TrEMBL file would.
  const bool uppercase = format != FlatFileFormat::swiss;
  for (std::size_t i = begin; i < end; ++i) {
    const char c = line.text[i];
    if (c == ' ') continue;
    ++count;
    if (out != nullptr) {
      out->push_back(uppercase && c >= 'a' && c <= 'z'
                         ? static_cast<char>(c - ('a' - 'A'))
                         : c);
    }
  }
  return true;
}

}  // namespace

const char* flatfile_format_name(FlatFileFormat format) noexcept {
  switch (format) {
    case FlatFileFormat::genbank:
      return "genbank";
    case FlatFileFormat::embl:
      return "embl";
    case FlatFileFormat::swiss:
      return "swiss";
  }
  return "unknown";
}

FlatFileIndex::FlatFileIndex(const char* data, std::size_t size,
                             FlatFileFormat format)
    : data_(data), size_(size), format_(format) {
  const bool genbank = format_ == FlatFileFormat::genbank;
  const std::size_t header_width = genbank ? kGenBankHeaderWidth : kEmblHeaderWidth;

  std::size_t cursor = 0;
  Line line = read_line(data_, size_, cursor);

  while (line.exists) {
    // Between records.  A blank line and a stray `//` are both nothing: the
    // reference's own record loop treats them the same way, and a file with a
    // trailing newline after its last `//` is not a malformed file.
    if (rstripped_size(line.text, line.len) == 0) {
      cursor = line.next;
      line = read_line(data_, size_, cursor);
      continue;
    }
    if (is_record_end(line)) {
      cursor = line.next;
      line = read_line(data_, size_, cursor);
      continue;
    }

    FlatFileRecord record;
    record.offset = static_cast<std::size_t>(line.text - data_);

    std::size_t sequence_version = 0;
    bool has_sequence_version = false;
    bool in_sequence = false;
    bool sequence_started = false;
    bool in_definition = false;
    std::size_t sequence_cursor = 0;

    // A GenBank definition's full stop belongs to the syntax, not to the data,
    // so the reference drops it -- once, when the DEFINITION block has been read
    // to its end.  Doing it here rather than at the end of the record is what
    // keeps a second DEFINITION line from having its own stop eaten as part of
    // the previous one's value.
    const auto close_definition = [&]() {
      if (!in_definition) return;
      in_definition = false;
      if (!record.description.empty() && record.description.back() == '.') {
        record.description.pop_back();
      }
    };

    // The feature block, skipped in one step rather than walked line by line.
    // Its lines are neither keywords nor residues, and reproducing the
    // reference's decision about where the block ends is the whole of this --
    // `feature_block_ends` is the same function the reader stops at, so the
    // span recorded here and the span walked there cannot disagree.
    //
    // It walks its own cursor and leaves `line` where it found it: the line that
    // *ends* the block is the sequence header, and it is the main loop's to
    // read.  Overwriting `line` here would skip that header, and the record
    // would come back with an empty sequence and no error -- a wrong answer of
    // exactly the kind this reader exists not to give.
    const auto scan_feature_block = [&](const Line& marker) {
      record.features_offset = static_cast<std::size_t>(marker.text - data_);
      Line scanned = read_line(data_, size_, marker.next);
      while (true) {
        if (!scanned.exists) {
          fail("Premature end of line during features table", record.offset);
          return;
        }
        if (is_record_end(scanned)) {
          fail("Premature end of features table, marker '//' found",
               static_cast<std::size_t>(scanned.text - data_));
          return;
        }
        if (feature_block_ends(format_, scanned.text, scanned.len)) break;
        scanned = read_line(data_, size_, scanned.next);
      }
      record.features_end = static_cast<std::size_t>(scanned.text - data_);
    };

    while (line.exists && !is_record_end(line)) {
      if (in_sequence) {
        if (!sequence_started) {
          sequence_cursor = static_cast<std::size_t>(line.text - data_);
          sequence_started = true;
        }
        std::size_t count = 0;
        if (!read_residue_line(format_, line, nullptr, count)) {
          fail("sequence line is not shaped the way " +
                   std::string(flatfile_format_name(format_)) +
                   " says it is, so this reader declines rather than guess",
               static_cast<std::size_t>(line.text - data_));
          return;
        }
        record.length += count;
        cursor = line.next;
        line = read_line(data_, size_, cursor);
        continue;
      }

      const std::string_view keyword = keyword_at(line, header_width);
      const std::string_view value = value_at(line, header_width);
      if (genbank && !keyword.empty() && keyword != "DEFINITION") {
        close_definition();
      }

      if (genbank) {
        if (keyword == "LOCUS") {
          record.name = first_word(value.data(), value.size());
          record.has_header_facts = genbank_locus_facts(line, record);
          in_definition = false;
        } else if (keyword == "FEATURES") {
          if (record.header_end == 0) {
            record.header_end = static_cast<std::size_t>(line.text - data_);
          }
          scan_feature_block(line);
          if (failed_) return;
          in_definition = false;
        } else if (keyword == "DEFINITION") {
          // A definition too long for one line continues on lines indented to
          // the keyword column, which is why they arrive with an empty keyword.
          const std::string_view text =
              rstrip_view(value.data(), value.size());
          if (record.description.empty()) {
            record.description.assign(text);
          } else {
            record.description += ' ';
            record.description.append(text);
          }
          in_definition = true;
        } else if (keyword == "ACCESSION") {
          if (record.id.empty()) {
            const std::pair<std::size_t, std::size_t> span =
                first_word_span(value.data(), value.size());
            record.id.assign(value.data() + span.first, span.second);
          }
          in_definition = false;
        } else if (keyword == "VERSION") {
          // The reference reads this line through `version()`, after two edits
          // the scanner makes first: runs of two spaces are collapsed, and a
          // ` GI:...` tail is cut off and handed to a *different* consumer, so
          // `U49845.1  GI:1293613` states a version of `U49845.1` and not a
          // version of the whole line.
          //
          // Two shapes of what is left, and they do different things.
          // `AB000001.2` splits into an accession and a version suffix -- the
          // suffix is what `record_end` hangs off the id, the accession is what
          // the annotations dict records.  A bare `AB000001` with no dot
          // *becomes* the id outright, overwriting what ACCESSION set.
          std::string text = collapse_spaces(std::string(strip_view(value)));
          const std::size_t gi = text.find(" GI:");
          if (gi != std::string::npos) text.resize(gi);
          apply_version(text, record, sequence_version, has_sequence_version);
          in_definition = false;
        } else if (keyword == "ORIGIN") {
          if (record.header_end == 0) {
            record.header_end = static_cast<std::size_t>(line.text - data_);
          }
          in_sequence = true;
        } else if (keyword.empty() && in_definition) {
          const std::size_t kept = rstripped_size(line.text, line.len);
          if (kept > header_width) {
            record.description += ' ';
            record.description.append(line.text + header_width,
                                      kept - header_width);
          }
        } else {
          in_definition = false;
        }
      } else {
        if (keyword == "ID") {
          // The ID line is where the reference takes both the name and -- for
          // EMBL, which has no other source for it -- the id itself.  Its first
          // field ends at the first semicolon; its second carries `SV n`, which
          // is a version suffix and not a field of its own, and which is what
          // later turns `X56734` into `X56734.2` in the record's id.
          const std::string_view trimmed = strip_view(value);
          const std::size_t semicolon = trimmed.find(';');
          const std::string_view first =
              semicolon == std::string_view::npos ? trimmed
                                                  : trimmed.substr(0, semicolon);
          record.name = first_word(first.data(), first.size());
          if (format_ == FlatFileFormat::embl) {
            record.id = record.name;
            record.has_header_facts = embl_id_facts(trimmed, record);
            // What the annotations reader needs from the same line: the
            // accession it states, and whether it stated a version suffix.  The
            // version is a key of the annotations dict and not just a piece of
            // the id, which is why both are kept.
            record.first_line_accession = record.name;
            if (semicolon != std::string_view::npos) {
              const std::string_view rest = trimmed.substr(semicolon + 1);
              const std::size_t next = rest.find(';');
              const std::string_view sv =
                  strip_view(next == std::string_view::npos ? rest
                                                            : rest.substr(0, next));
              const std::size_t space = sv.find(' ');
              if (space != std::string_view::npos) {
                const std::string_view key = strip_view(sv.substr(0, space));
                const std::string_view number = strip_view(sv.substr(space + 1));
                if (key == "SV" && all_digits(number.data(), number.size())) {
                  sequence_version = digits_to_size(number);
                  has_sequence_version = true;
                  record.first_line_sequence_version =
                      static_cast<std::int64_t>(sequence_version);
                  record.has_first_line_sequence_version = true;
                }
              }
            }
          }
        } else if (keyword == "AC") {
          const std::string_view trimmed = strip_view(value);
          const std::size_t semicolon = trimmed.find(';');
          const std::string_view first =
              semicolon == std::string_view::npos ? trimmed
                                                  : trimmed.substr(0, semicolon);
          if (record.id.empty()) {
            record.id = first_word(first.data(), first.size());
          }
        } else if (keyword == "SV" && format_ == FlatFileFormat::embl) {
          // The obsolete version line, removed from the format in June 2006 and
          // still found in old files.  The reference maps it to the same
          // consumer as `VERSION` -- which is what makes its suffix reach the
          // record's id, a fact that surprises until `record_end` is read.
          apply_version(collapse_spaces(std::string(strip_view(value))), record,
                        sequence_version, has_sequence_version);
        } else if (keyword == "DE") {
          // SwissProt repeats the keyword on every line of a description and so
          // does EMBL, so a description is a join and never a continuation.
          const std::string_view text = rstrip_view(value.data(), value.size());
          if (record.description.empty()) {
            record.description.assign(text);
          } else {
            record.description += ' ';
            record.description.append(text);
          }
        } else if (keyword == "FH" && format_ == FlatFileFormat::embl) {
          if (record.header_end == 0) {
            record.header_end = static_cast<std::size_t>(line.text - data_);
          }
          scan_feature_block(line);
          if (failed_) return;
        } else if (keyword == "SQ") {
          if (record.header_end == 0) {
            record.header_end = static_cast<std::size_t>(line.text - data_);
          }
          in_sequence = true;
        }
      }

      cursor = line.next;
      line = read_line(data_, size_, cursor);
    }

    // Finishing touches, in the reference's order.  An id is a promise that a
    // record can be addressed; a GenBank or EMBL record with no ACCESSION falls
    // back to its locus name, and its id with no dot takes the version suffix.
    // SwissProt takes neither, because its reader never asks -- and a record id
    // is the one field a caller addresses a record by, so getting this wrong is
    // a lookup that fails rather than a field that reads oddly.
    close_definition();
    if (format_ != FlatFileFormat::swiss && record.id.empty()) {
      record.id = record.name;
    }
    if (has_sequence_version && record.id.find('.') == std::string::npos) {
      record.id += '.';
      record.id += std::to_string(sequence_version);
    }
    // `cursor` is the offset of the line that ended the record -- the `//`, or
    // the end of the buffer -- because the inner loop re-reads at it and never
    // advances past it.  A record whose ORIGIN block is empty therefore gets an
    // empty span rather than a span reaching into the next record.
    record.sequence_end = cursor;
    record.sequence_offset = sequence_started ? sequence_cursor : cursor;
    // A record with no feature table and no origin line has no header either, and
    // the span is then empty rather than a span reaching past the record.
    if (record.header_end == 0) record.header_end = record.sequence_offset;

    if (!by_id_.emplace(record.id, records_.size()).second) {
      // The reference's own wording for the same refusal, because a duplicate
      // id is a file a caller cannot address rather than an implementation
      // detail of this index.
      fail("Duplicate key '" + record.id + "'", record.offset);
      return;
    }
    records_.push_back(std::move(record));

    if (!line.exists) break;
    cursor = line.next;
    line = read_line(data_, size_, cursor);
  }
}

void FlatFileIndex::fail(std::string message, std::size_t offset) {
  error_message_ = std::move(message) + " (byte offset " + std::to_string(offset) + ")";
  failed_ = true;
  records_.clear();
  by_id_.clear();
}

const FlatFileRecord* FlatFileIndex::find(std::string_view id) const noexcept {
  const auto it = by_id_.find(std::string(id));
  if (it == by_id_.end()) return nullptr;
  return &records_[it->second];
}

std::string FlatFileIndex::sequence(const FlatFileRecord& record) const {
  std::string out;
  out.reserve(record.length);
  std::size_t cursor = record.sequence_offset;
  while (cursor < record.sequence_end) {
    const Line line = read_line(data_, size_, cursor);
    std::size_t count = 0;
    if (!read_residue_line(format_, line, &out, count)) return std::string();
    cursor = line.next;
  }
  return out;
}

std::string FlatFileIndex::sequence_slice(const FlatFileRecord& record,
                                          std::size_t start,
                                          std::size_t end) const {
  start = std::min(start, record.length);
  end = std::min(end, record.length);
  if (end <= start) return std::string();

  // The residues of a flat-file record are not one span of the file: they are
  // interleaved with a coordinate column and, in EMBL, a trailing count.  So a
  // window is walked line by line rather than found by arithmetic -- but only
  // the window is materialised, which is the difference between a slice and a
  // whole-record parse.
  std::string out;
  out.reserve(end - start);
  std::string residues;
  std::size_t position = 0;
  std::size_t cursor = record.sequence_offset;
  while (cursor < record.sequence_end && position < end) {
    const Line line = read_line(data_, size_, cursor);
    residues.clear();
    std::size_t count = 0;
    if (!read_residue_line(format_, line, &residues, count)) return std::string();
    const std::size_t line_end = position + count;
    if (line_end > start) {
      const std::size_t from = position > start ? 0 : start - position;
      const std::size_t to = std::min(count, end - position);
      out.append(residues, from, to - from);
    }
    position = line_end;
    cursor = line.next;
  }
  return out;
}

// The record's feature table.  Three refusals are possible here and they mean
// different things, so they are kept apart:
//
//   * no feature block at all -- an empty table, which is an answer;
//   * a block whose header facts are missing -- a refusal, because the reference
//     would have read the locations under a declared size and topology this
//     reader could not establish;
//   * a block containing a shape the reference only warns about -- a refusal
//     from the reader below, which says which line it was.
//
// SwissProt's `FT` block is a different grammar and is refused outright rather
// than reported as empty, because "no features" and "features not reproduced"
// are not the same statement.
FeatureTable FlatFileIndex::features(const FlatFileRecord& record) const {
  if (format_ == FlatFileFormat::swiss) {
    return parse_feature_table(format_, data_, size_, record.offset,
                               record.offset, record.declared_size,
                               record.has_declared_size, record.circular,
                               record.stranded);
  }
  if (record.features_end <= record.features_offset) return FeatureTable();
  if (!record.has_header_facts) {
    FeatureTable table;
    table.refused = true;
    table.offset = record.offset;
    table.message =
        "the header line's declared size and topology are not in a layout this "
        "reader reproduces, so the feature locations would be read under inputs "
        "the reference never used";
    return table;
  }
  return parse_feature_table(format_, data_, size_, record.features_offset,
                             record.features_end, record.declared_size,
                             record.has_declared_size, record.circular,
                             record.stranded);
}

// The header, read out of the span the scan recorded rather than re-walked from
// the record's first byte: the line that ended the header is already known, and
// finding it twice is how a reader comes to disagree with itself about where a
// block ended.
AnnotationTable FlatFileIndex::annotations(const FlatFileRecord& record) const {
  return read_annotations(format_, data_, size_, record.offset,
                          record.header_end, record);
}

}  // namespace biofasting