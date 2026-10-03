#include "annotations.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "line.hpp"

namespace biofasting {

const char* annotation_key_name(AnnotationKey key) noexcept {
  switch (key) {
    case AnnotationKey::molecule_type: return "molecule_type";
    case AnnotationKey::topology: return "topology";
    case AnnotationKey::data_file_division: return "data_file_division";
    case AnnotationKey::date: return "date";
    case AnnotationKey::accessions: return "accessions";
    case AnnotationKey::sequence_version: return "sequence_version";
    case AnnotationKey::gi: return "gi";
    case AnnotationKey::keywords: return "keywords";
    case AnnotationKey::source: return "source";
    case AnnotationKey::organism: return "organism";
    case AnnotationKey::taxonomy: return "taxonomy";
    case AnnotationKey::references: return "references";
    case AnnotationKey::comment: return "comment";
  }
  return "";
}

void AnnotationTable::note(AnnotationKey key) {
  for (const AnnotationKey seen : order) {
    if (seen == key) return;
  }
  order.push_back(key);
}

bool AnnotationTable::has(AnnotationKey key) const noexcept {
  for (const AnnotationKey seen : order) {
    if (seen == key) return true;
  }
  return false;
}

namespace {

constexpr std::size_t kGenBankIndent = 12;
constexpr std::size_t kEmblIndent = 5;

std::string stripped(std::string_view text) { return std::string(strip_view(text)); }

bool ends_with(std::string_view text, std::string_view suffix) noexcept {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// `str.rstrip(c)`, which the EMBL consumer uses to drop a trailing semicolon.
std::string rstripped_of(std::string_view text, char c) {
  std::size_t end = text.size();
  while (end > 0 && text[end - 1] == c) --end;
  return std::string(text.substr(0, end));
}

bool contains(std::string_view haystack, std::string_view needle) noexcept {
  return haystack.find(needle) != std::string_view::npos;
}

// Python's `str.split(None)`: runs of whitespace, no empty fields.  The
// whitespace set is `is_python_space`, the same one the line reader uses, so a
// value holding a form feed splits here the way it splits there.
std::vector<std::string> split_whitespace(std::string_view text) {
  std::vector<std::string> out;
  std::size_t i = 0;
  while (i < text.size()) {
    while (i < text.size() && is_python_space(text[i])) ++i;
    const std::size_t start = i;
    while (i < text.size() && !is_python_space(text[i])) ++i;
    if (i > start) out.emplace_back(text.substr(start, i - start));
  }
  return out;
}

// `_split_accessions`: semicolons and line feeds are separators, then
// whitespace.  Used for ACCESSION, whose value the scanner has already joined.
std::vector<std::string> split_accessions(std::string_view text) {
  std::string normalised;
  normalised.reserve(text.size());
  for (const char c : text) normalised.push_back(c == ';' ? ' ' : c);
  return split_whitespace(normalised);
}

// `_split_keywords`: an empty field or a lone `.` means one empty keyword rather
// than none, which is why GenBank's thin header carries `['']` and an EMBL
// record with no `KW` line carries no `keywords` key at all.
std::vector<std::string> split_keywords(std::string_view text) {
  std::string_view body = text;
  if (body.empty() || body == ".") {
    body = std::string_view();
  } else if (ends_with(body, ".")) {
    body = body.substr(0, body.size() - 1);
  }
  std::vector<std::string> out;
  std::size_t start = 0;
  while (true) {
    const std::size_t at = body.find(';', start);
    const std::string_view piece = at == std::string_view::npos
                                       ? body.substr(start)
                                       : body.substr(start, at - start);
    out.push_back(stripped(piece));
    if (at == std::string_view::npos) break;
    start = at + 1;
  }
  return out;
}

// `_split_taxonomy`: nothing at all, or a lone `.`, is an empty list -- the
// difference between `taxonomy: []` and a list holding one empty string.
std::vector<std::string> split_taxonomy(std::string_view text) {
  std::string_view body = text;
  if (body.empty() || body == ".") return {};
  if (ends_with(body, ".")) body = body.substr(0, body.size() - 1);
  std::vector<std::string> out;
  std::size_t start = 0;
  while (true) {
    const std::size_t at = body.find(';', start);
    const std::string_view piece = at == std::string_view::npos
                                       ? body.substr(start)
                                       : body.substr(start, at - start);
    out.push_back(stripped(piece));
    if (at == std::string_view::npos) break;
    start = at + 1;
  }
  return out;
}

void append_unique(std::vector<std::string>& into, const std::string& value) {
  if (std::find(into.begin(), into.end(), value) != into.end()) return;
  into.push_back(value);
}

bool to_int64(std::string_view text, std::int64_t& out) noexcept {
  if (text.empty()) return false;
  std::size_t i = 0;
  bool negative = false;
  if (text[0] == '+' || text[0] == '-') {
    negative = text[0] == '-';
    i = 1;
  }
  if (i == text.size()) return false;
  std::int64_t value = 0;
  for (; i < text.size(); ++i) {
    if (text[i] < '0' || text[i] > '9') return false;
    value = value * 10 + (text[i] - '0');
  }
  out = negative ? -value : value;
  return true;
}

// `_split_reference_locations`: `1 to 20; 20 to 100` becomes the pairs the
// reference turns into `SimpleLocation(start - 1, end)`.  The list is *assigned*
// rather than extended, because `reference_bases` does `location = all_locations`
// and an EMBL block with two `RP` lines ends up with the second one's.
bool reference_locations(std::string_view text, AnnotationReference& reference,
                         std::string& message) {
  reference.location.clear();
  std::size_t start = 0;
  while (true) {
    const std::size_t at = text.find(';', start);
    const std::string_view piece = at == std::string_view::npos
                                       ? text.substr(start)
                                       : text.substr(start, at - start);
    const std::size_t to = piece.find("to");
    if (to == std::string_view::npos) {
      message = "reference location without a 'to' in it";
      return false;
    }
    std::int64_t begin = 0;
    std::int64_t end = 0;
    if (!to_int64(strip_view(piece.substr(0, to)), begin) ||
        !to_int64(strip_view(piece.substr(to + 2)), end)) {
      message = "reference location is not two integers";
      return false;
    }
    reference.location.emplace_back(begin - 1, end);
    if (at == std::string_view::npos) break;
    start = at + 1;
  }
  return true;
}

// `reference_bases`, whose argument the scanner has already wrapped in brackets:
// `(bases 1 to 1000)`, `(sites)`, `(bases 1 to 105; 110 to 111)`, or the
// `(residues 1 to 182)` shape an EMBL patent file uses.
bool reference_bases(std::string_view content, AnnotationReference& reference,
                     std::string& message) {
  if (!ends_with(content, ")")) {
    message = "a reference's bases are not bracketed";
    return false;
  }
  const std::string_view body = content.substr(1, content.size() - 2);
  const bool has_to = contains(body, "to");
  if (contains(body, "bases") && has_to) {
    return reference_locations(body.substr(5), reference, message);
  }
  if (contains(body, "residues") && has_to) {
    const std::size_t at = body.find("residues");
    return reference_locations(body.substr(at + 9), reference, message);
  }
  if (body == "sites" || strip_view(body) == "bases") {
    reference.location.clear();
    return true;
  }
  message = "a reference's bases are in a shape this reader does not reproduce";
  return false;
}

// One record's header, walked line by line.  The walk owns the span and the four
// facts the first line stated; every other line is a block.
//
// Cursor discipline is the whole of the difficulty here and it is not this
// reader's invention: the reference's GenBank consumer looks ahead with
// `next(line_iter)` and leaves the line it stopped on for the loop to see again,
// while its EMBL consumer walks a list and never looks ahead at all.  A handler
// below therefore either stops on the line that ended its block -- and the loop
// processes that line next -- or advances past the line it consumed.  Getting
// this wrong skips a line silently, which is the failure mode this file must not
// have.
class HeaderReader {
 public:
  HeaderReader(FlatFileFormat format, const char* data, std::size_t size,
               std::size_t begin, std::size_t end, const FlatFileRecord& record)
      : format_(format),
        data_(data),
        size_(size),
        end_(end),
        indent_(format == FlatFileFormat::genbank ? kGenBankIndent : kEmblIndent),
        record_(record) {
    cursor_ = begin;
    line_ = read_line(data_, size_, cursor_);
  }

  AnnotationTable run() {
    AnnotationTable table;
    table_ = &table;
    if (!record_.has_annotation_header) {
      fail("this record's first line is not in a layout this reader "
           "reproduces, so its annotations are refused rather than guessed");
      return table;
    }
    // The first line's keys, in the order the reference's `_feed_first_line`
    // creates them.  That order is the handler's own and not the file's -- there
    // is one first line -- and it differs between the two formats: EMBL's
    // handler records the accession before the topology, GenBank's cannot,
    // because for GenBank the accession comes from its own `ACCESSION` line
    // further down the header.
    //
    // Each is created only when the reference would have created it: its
    // `molecule_type` and `topology` handlers skip an empty value, GenBank's
    // division and date are behind a non-blank test, and EMBL's division is not
    // -- an EMBL ID line always states a division, so the key always exists.
    if (format_ == FlatFileFormat::embl) {
      table.accessions = split_accessions(record_.first_line_accession);
      table.note(AnnotationKey::accessions);
      if (record_.has_first_line_sequence_version) {
        table.sequence_version = record_.first_line_sequence_version;
        table.note(AnnotationKey::sequence_version);
      }
      if (!record_.topology.empty()) {
        table.topology = record_.topology;
        table.note(AnnotationKey::topology);
      }
      if (!record_.molecule_type.empty()) {
        table.molecule_type = record_.molecule_type;
        table.note(AnnotationKey::molecule_type);
      }
      table.data_file_division = record_.data_file_division;
      table.note(AnnotationKey::data_file_division);
    } else {
      if (!record_.molecule_type.empty()) {
        table.molecule_type = record_.molecule_type;
        table.note(AnnotationKey::molecule_type);
      }
      if (!record_.topology.empty()) {
        table.topology = record_.topology;
        table.note(AnnotationKey::topology);
      }
      if (!record_.data_file_division.empty()) {
        table.data_file_division = record_.data_file_division;
        table.note(AnnotationKey::data_file_division);
      }
      if (!record_.date.empty()) {
        table.date = record_.date;
        table.note(AnnotationKey::date);
      }
    }
    if (format_ == FlatFileFormat::genbank) {
      genbank_header();
    } else if (format_ == FlatFileFormat::embl) {
      embl_header();
    } else {
      fail("the annotations dict of this format is not reproduced");
    }
    return table;
  }

 private:
  bool done() const noexcept {
    if (!line_.exists) return true;
    return static_cast<std::size_t>(line_.text - data_) >= end_;
  }

  void advance() {
    cursor_ = line_.next;
    line_ = read_line(data_, size_, cursor_);
  }

  std::size_t offset() const noexcept {
    return static_cast<std::size_t>(line_.text - data_);
  }

  // The value the reference's header consumer takes: `line[12:].strip()`.
  std::string_view value() const noexcept {
    return line_.len <= indent_
               ? std::string_view()
               : strip_view(std::string_view(line_.text + indent_,
                                             line_.len - indent_));
  }

  // The value with only its right-hand end trimmed.  Continuation lines are
  // joined with this one, so their leading spaces survive into the joined
  // string, exactly as they do in the reference.
  std::string_view raw_value() const noexcept {
    return line_.len <= indent_ ? std::string_view()
                                : std::string_view(line_.text + indent_,
                                                   line_.len - indent_);
  }

  std::string_view keyword() const noexcept { return keyword_at(line_, indent_); }

  bool continuation() const noexcept {
    if (line_.len < indent_) return false;
    for (std::size_t i = 0; i < indent_; ++i) {
      if (line_.text[i] != ' ') return false;
    }
    return true;
  }

  void fail(std::string message) {
    if (table_->failed) return;
    table_->failed = true;
    table_->error_message = std::move(message);
    table_->error_offset = line_.exists ? offset() : end_;
  }

  // The reference's generic continuation loop: keep taking lines indented to the
  // keyword column while they arrive, then stop -- leaving the line that stopped
  // it where the loop will see it again.
  std::string collect_joined() {
    std::string out(value());
    while (true) {
      advance();
      if (done() || !continuation()) break;
      out += ' ';
      out += std::string(raw_value());
    }
    return out;
  }

  static std::string join(const std::string& existing, const std::string& piece) {
    if (existing.empty()) return piece;
    return existing + " " + piece;
  }

  // `while "  " in data: data = data.replace("  ", " ")`, which is how a
  // REFERENCE line's continuation joins without leaving a double space behind.
  static std::string collapse(std::string_view text) {
    std::string out(text);
    std::size_t at = out.find("  ");
    while (at != std::string::npos) {
      out.replace(at, 2, " ");
      at = out.find("  ", at);
    }
    return out;
  }

  // ---- GenBank -----------------------------------------------------------

  void genbank_header() {
    while (!done() && !table_->failed) {
      const std::string_view key = keyword();
      //   The keys the reference's own consumer table names, in the shapes this
      //   reader reproduces.  Everything else is a keyword the reference ignores
      //   -- unless it is one that sets a key of the annotations dict, which is
      //   refused rather than dropped.
      if (key == "LOCUS") {
        advance();
      } else if (key == "ACCESSION") {
        accession(collect_joined());
      } else if (key == "VERSION") {
        genbank_version(value());
        advance();  // this branch is the one the reference does not continue
      } else if (key == "KEYWORDS") {
        keywords(collect_joined());
      } else if (key == "SOURCE") {
        source(collect_joined());
      } else if (key == "ORGANISM") {
        organism();
      } else if (key == "COMMENT") {
        genbank_comment();
      } else if (key == "REFERENCE") {
        genbank_reference();
      } else if (key == "AUTHORS") {
        AnnotationReference& ref = current_reference(key);
        ref.authors = join(ref.authors, collect_joined());
      } else if (key == "TITLE") {
        AnnotationReference& ref = current_reference(key);
        ref.title = join(ref.title, collect_joined());
      } else if (key == "JOURNAL") {
        AnnotationReference& ref = current_reference(key);
        ref.journal = join(ref.journal, collect_joined());
      } else if (key == "CONSRTM") {
        AnnotationReference& ref = current_reference(key);
        ref.consrtm = join(ref.consrtm, collect_joined());
      } else if (key == "PUBMED") {
        current_reference(key).pubmed_id = collect_joined();
      } else if (key == "MEDLINE") {
        current_reference(key).medline_id = collect_joined();
      } else if (key == "REMARK") {
        AnnotationReference& ref = current_reference(key);
        ref.comment = join(ref.comment, collect_joined());
      } else if (key == "DEFINITION" || key == "PROJECT" || key == "DBLINK") {
        //   The description, which the scan has already read, and the two keys
        //   that go to `SeqRecord.dbxrefs` -- a field this package does not
        //   produce, so ignoring them here cannot make the annotations dict
        //   wrong.
        collect_joined();
      } else if (key == "NID" || key == "PID" || key == "DBSOURCE" ||
                 key == "SEGMENT") {
        fail("GenBank keyword '" + std::string(key) +
             "' sets an annotation this reader does not reproduce");
      } else {
        advance();  // a keyword the reference reads and ignores
      }
    }
  }

  // `comment()` appends `"\n" + "\n".join(content)`, so a second COMMENT block
  // continues the first with one newline between them, and the continuation
  // lines keep the writer's own wrapping because they are appended unstripped.
  void genbank_comment() {
    const std::string_view first = value();
    if (ends_with(first, "-START##")) {
      // The structured-comment probe -- the one that costs 16.4 us per record in
      // bench/targets.md -- has fired, and the reference would now build a
      // second annotation key out of the block.
      fail("a structured COMMENT block, which this reader does not reproduce");
      return;
    }
    std::string block(first);
    while (true) {
      advance();
      if (done() || !continuation()) break;
      block += '\n';
      block += std::string(raw_value());
    }
    if (table_->has(AnnotationKey::comment)) {
      table_->comment += '\n';
      table_->comment += block;
    } else {
      table_->comment = block;
    }
    table_->note(AnnotationKey::comment);
  }

  // `REFERENCE    1  (bases 1 to 86436)`, its continuation lines joined with a
  // space, then split at the first space into a number and a bracketed range.
  void genbank_reference() {
    std::string data = collect_joined();
    data = collapse(data);
    const std::size_t space = data.find(' ');
    table_->references.emplace_back();
    table_->note(AnnotationKey::references);
    AnnotationReference& ref = table_->references.back();
    if (space == std::string_view::npos) return;  // a number and no bases
    if (!reference_bases(std::string_view(data).substr(space + 1), ref, error_)) {
      fail(error_);
    }
  }

  // A `TITLE`/`AUTHORS`/... line above any `REFERENCE` is a record the reference
  // warns about and then raises on, so this reader declines instead of writing
  // the value into a reference object that does not exist.
  //
  // One function for both formats, because both loops run in this class and the
  // "no reference yet" state is the same state.  The line it names is spelled
  // the way the file spells it: a reader told about a `REFERENCE` line in an
  // EMBL file that has no `REFERENCE` keyword would go looking for a line that
  // is not there.
  AnnotationReference& current_reference(std::string_view key) {
    if (table_->references.empty()) {
      fail("'" + std::string(key) + "' line before any " +
           (format_ == FlatFileFormat::embl ? "RN" : "REFERENCE") + " line");
      return placeholder_;
    }
    return table_->references.back();
  }

  void accession(std::string_view text) {
    for (const std::string& acc : split_accessions(text)) {
      append_unique(table_->accessions, acc);
    }
    table_->note(AnnotationKey::accessions);
  }

  // `version()`: `AB000001.2` is an accession plus a sequence version, and
  // anything else becomes the record's id outright -- which the scan has already
  // done.  Shared by the GenBank `VERSION` line and the obsolete EMBL `SV` line,
  // because the reference maps both to this one consumer.
  void version(std::string_view text) {
    const std::string data(text);
    const std::size_t dot = data.find('.');
    if (dot == std::string::npos || data.find('.', dot + 1) != std::string::npos) {
      return;
    }
    std::int64_t number = 0;
    if (!to_int64(std::string_view(data).substr(dot + 1), number)) return;
    append_unique(table_->accessions, data.substr(0, dot));
    table_->note(AnnotationKey::accessions);
    table_->sequence_version = number;
    table_->note(AnnotationKey::sequence_version);
  }

  // The GenBank `VERSION` line, which the reference's *scanner* edits before
  // handing it to the consumer above: runs of two spaces collapse, and a
  // ` GI:...` tail is cut off and given to a different consumer, so
  // `U49845.1  GI:1293613` states a version of `U49845.1` and a GI number of its
  // own.  The split lives here rather than in `version()` because the EMBL `SV`
  // line shares the consumer but not this edit.
  void genbank_version(std::string_view text) {
    const std::string data = collapse(text);
    const std::size_t gi = data.find(" GI:");
    if (gi == std::string::npos) {
      version(data);
      return;
    }
    version(std::string_view(data).substr(0, gi));
    // Python's `data.split(" GI:")[1]`: what lies between the first occurrence
    // and the second, which for the one-occurrence shape every real file has is
    // simply the tail.
    const std::size_t next = data.find(" GI:", gi + 1);
    table_->gi = next == std::string::npos
                     ? data.substr(gi + 4)
                     : data.substr(gi + 4, next - gi - 4);
    table_->note(AnnotationKey::gi);
  }

  void keywords(std::string_view text) {
    const std::vector<std::string> split = split_keywords(text);
    if (table_->has(AnnotationKey::keywords)) {
      table_->keywords.insert(table_->keywords.end(), split.begin(), split.end());
    } else {
      table_->keywords = split;
    }
    table_->note(AnnotationKey::keywords);
  }

  // `source()` drops one trailing full stop and nothing else: an empty SOURCE
  // line is the empty string, and a lone `.` becomes one too.
  void source(std::string_view text) {
    std::string info(text);
    if (!info.empty() && info.back() == '.') info.pop_back();
    table_->source = info;
    table_->note(AnnotationKey::source);
  }

  // The ORGANISM block: the first line is the organism and the lines under it
  // are the lineage, but a long species name wraps, and the reference tells the
  // two apart by a semicolon or by one of eight root names.  Reproduced because
  // guessing differently moves a species into its own taxonomy.
  void organism() {
    table_->organism = std::string(value());
    table_->note(AnnotationKey::organism);
    std::string lineage;
    while (true) {
      advance();
      if (done() || !continuation()) break;
      const std::string_view piece = raw_value();
      const std::string_view trimmed = strip_view(piece);
      if (!lineage.empty() || contains(piece, ";") || trimmed == "Bacteria." ||
          trimmed == "Archaea." || trimmed == "Eukaryota." ||
          trimmed == "Unclassified." || trimmed == "Viruses." ||
          trimmed == "cellular organisms." || trimmed == "other sequences." ||
          trimmed == "unclassified sequences.") {
        lineage += ' ';
        lineage += std::string(piece);
      } else if (trimmed == ".") {
        // No lineage data, just a placeholder.
      } else {
        table_->organism += ' ';
        table_->organism += std::string(piece);
      }
    }
    const std::vector<std::string> split = split_taxonomy(strip_view(lineage));
    if (table_->has(AnnotationKey::taxonomy)) {
      table_->taxonomy.insert(table_->taxonomy.end(), split.begin(), split.end());
    } else {
      table_->taxonomy = split;
    }
    table_->note(AnnotationKey::taxonomy);
  }

  // ---- EMBL --------------------------------------------------------------

  // The EMBL consumer is a `for line in lines` loop: it never looks ahead, so
  // every branch here consumes exactly its own line.  A reference block is
  // several independent lines rather than a block with its own loop, which is
  // why `RT` split over two lines is two appends and not one join.
  void embl_header() {
    while (!done() && !table_->failed) {
      const std::string_view key = keyword();
      if (key == "AC") {
        accession(value());
      } else if (key == "DE") {
        ;  // the description, which the scan has already read
      } else if (key == "KW") {
        keywords(rstripped_of(value(), ';'));
      } else if (key == "OS") {
        table_->organism = std::string(value());
        table_->note(AnnotationKey::organism);
      } else if (key == "OC") {
        //   The lineage line, in the same consumer GenBank's ORGANISM block ends
        //   in.  EMBL has one `OC` line per lineage chunk rather than an
        //   indented block, so each line extends the list.
        const std::vector<std::string> split = split_taxonomy(value());
        if (table_->has(AnnotationKey::taxonomy)) {
          table_->taxonomy.insert(table_->taxonomy.end(), split.begin(),
                                  split.end());
        } else {
          table_->taxonomy = split;
        }
        table_->note(AnnotationKey::taxonomy);
      } else if (key == "CC") {
        embl_comment();
      } else if (key == "RN") {
        std::string number(stripped(value()));
        if (number.size() >= 2 && number.front() == '[' && number.back() == ']') {
          number = number.substr(1, number.size() - 2);
        }
        table_->references.emplace_back();
        table_->note(AnnotationKey::references);
      } else if (key == "RP") {
        embl_reference_bases();
      } else if (key == "RX") {
        embl_cross_reference();
      } else if (key == "RA") {
        AnnotationReference& ref = current_reference(key);
        ref.authors = join(ref.authors, rstripped_of(stripped(value()), ';'));
      } else if (key == "RG") {
        // Unlike `RA` and `KW`, the reference passes this one through whole: the
        // consortium keeps its trailing semicolon, which is what `SeqIO.parse`
        // reports and therefore what this reader has to report too.
        AnnotationReference& ref = current_reference(key);
        ref.consrtm = join(ref.consrtm, stripped(value()));
      } else if (key == "RT") {
        AnnotationReference& ref = current_reference(key);
        ref.title = join(ref.title, embl_title());
      } else if (key == "RL") {
        AnnotationReference& ref = current_reference(key);
        ref.journal = join(ref.journal, stripped(value()));
      } else if (key == "SV") {
        // The obsolete version line, which the reference reads with the same
        // consumer the GenBank VERSION line goes to.
        version(value());
      } else if (key == "ID" || key == "XX" || key == "FH") {
        ;  // the first line, which the scan read; the block separator
      } else if (key.empty()) {
        advance();  // a continuation line of a block that has no continuations
      } else {
        //   `DT`, `DR`, `PR` and anything else: the reference's EMBL consumer
        //   drops a header line it has no entry for, so this drops it too.
        //   Nothing it drops sets a key of the annotations dict.
        advance();
      }
      if (table_->failed) return;
      advance();
    }
  }

  void embl_comment() {
    const std::string block(value());
    if (table_->has(AnnotationKey::comment)) {
      table_->comment += '\n';
      table_->comment += block;
    } else {
      table_->comment = block;
    }
    table_->note(AnnotationKey::comment);
  }

  // `RP   1-4639675` becomes `(bases 1 to 4639675)`, and
  // `RP   160-550, 904-1055` becomes `(bases 160 to 550; 904 to 1055)`, because
  // the EMBL consumer reformats the range into the shape the GenBank consumer
  // parses.  `RP   [-]`, which KIPO patent files use, carries nothing.
  void embl_reference_bases() {
    const std::string_view data = strip_view(value());
    if (data == "[-]") return;
    std::string bases = "(bases ";
    std::size_t start = 0;
    bool first = true;
    while (true) {
      const std::size_t at = data.find(',', start);
      const std::string_view piece =
          strip_view(at == std::string_view::npos ? data.substr(start)
                                                  : data.substr(start, at - start));
      if (!piece.empty()) {
        if (!first) bases += "; ";
        std::string converted(piece);
        for (std::size_t i = 0; i < converted.size();) {
          if (converted[i] == '-') {
            converted.replace(i, 1, " to ");
            i += 4;
          } else {
            ++i;
          }
        }
        bases += stripped(converted);
        first = false;
      }
      if (at == std::string_view::npos) break;
      start = at + 1;
    }
    bases += ")";
    if (!reference_bases(bases, current_reference("RP"), error_)) fail(error_);
  }

  // `RX   PUBMED; 264242.` -- the reference reads PUBMED and, despite the
  // `medline_id` field its `Reference` carries, nothing else: the EMBL consumer
  // has `if key == "PUBMED"` and a TODO for the rest, so an `RX   MEDLINE` line
  // is dropped rather than refused, and a record carrying one has an empty
  // `medline_id` in `SeqIO.parse` that this reader has to reproduce.
  void embl_cross_reference() {
    const std::string_view data = strip_view(value());
    const std::size_t semicolon = data.find(';');
    if (semicolon == std::string_view::npos) {
      fail("an RX line with no ';' in it");
      return;
    }
    const std::string_view name = strip_view(data.substr(0, semicolon));
    std::string identifier = stripped(data.substr(semicolon + 1));
    if (!identifier.empty() && identifier.back() == '.') identifier.pop_back();
    identifier = stripped(identifier);
    AnnotationReference& ref = current_reference("RX");
    if (name == "PUBMED") {
      ref.pubmed_id = identifier;
    }
  }

  // `RT   "Synthetic study 1";` loses its enclosing quotes and its trailing
  // semicolon; a title split over two lines gives the second line no quotes, and
  // the two are joined with a space.
  std::string embl_title() const {
    std::string title(stripped(value()));
    if (!title.empty() && title.front() == '"') title.erase(0, 1);
    if (ends_with(title, "\";")) title.resize(title.size() - 2);
    return title;
  }

  FlatFileFormat format_;
  const char* data_;
  std::size_t size_;
  std::size_t end_;
  std::size_t indent_;
  const FlatFileRecord& record_;
  std::size_t cursor_ = 0;
  Line line_;
  AnnotationTable* table_ = nullptr;
  AnnotationReference placeholder_;
  std::string error_;
};

}  // namespace

AnnotationTable read_annotations(FlatFileFormat format, const char* data,
                                 std::size_t size, std::size_t begin,
                                 std::size_t end, const FlatFileRecord& record) {
  HeaderReader reader(format, data, size, begin, end, record);
  return reader.run();
}

}  // namespace biofasting
