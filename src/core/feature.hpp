// The FEATURES table: the qualifier half of what a flat-file parse costs.
//
// location.hpp reads one feature's location string; this reads the block around
// it -- the feature key, the location line as the scanner assembles it, and the
// qualifiers with all of their continuation and quoting rules.  The two halves
// together are M19's measured 16 us of the reference's 20.2 us per feature.
//
// **One grammar, two spellings.**  GenBank and EMBL are read by one scanner
// class in the reference (`Bio/GenBank/Scanner.py`), differing in four
// constants: the header width, the feature-key column, the continuation spacer
// and the end markers.  SwissProt's FT lines are a *different* grammar, read by
// `Bio.SwissProt._read_ft`, and are deliberately not here -- see the header of
// flatfile.hpp.
//
// **Where it declines.**  The reference warns and carries on in a dozen places
// inside this block: a line too short to hold a feature, a feature over-indented
// because its key was too long, a location that wrapped without a trailing
// comma, a qualifier with a space after its `=`.  Every one of those is a shape
// whose reading is a guess about what the writer meant, and this project's rule
// is that a guess reproduced differently is a silent wrong answer.  So a table
// containing one is *refused* -- the whole index fails, naming the line -- and
// the caller falls back to `SeqIO.parse`.  Two things are not in that class and
// are reproduced instead, because they are properties of a *value* rather than
// of the file's shape and real files contain them:
//
//   * a location the reference's own parser rejects (`LocationParserError`),
//     which its consumer turns into a missing location plus a warning.  That is
//     a behaviour, and `LocationResult::status` carries it;
//   * a qualifier value with an unescaped double quote, which the reference
//     warns about and then keeps.  `FeatureQualifier::escape_warning` is the
//     flag, and the Python layer says the words.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "flatfile.hpp"
#include "location.hpp"

namespace biofasting {

struct FeatureQualifier {
  std::string key;
  // The value the reference's consumer ends up with: the continuation join
  // applied, one enclosing quote pair removed, doubled quotes collapsed, and a
  // `translation` stripped of whitespace.  `has_value` false is the `/pseudo`
  // form, where the reference stores `""` for a key it has not seen before and
  // drops the key entirely if it has.
  std::string value;
  bool has_value = true;
  // The reference warns when a value contains a double quote that is not part of
  // a doubled pair, because the NCBI's rule is that a quote inside a value is
  // written twice.  `escape_text` is the value the warning is *about* -- after
  // the enclosing quotes come off and before the doubling is undone -- and it is
  // kept only when the flag is set, so the Python layer can print the
  // reference's own sentence about the string the reference complained of.
  bool escape_warning = false;
  std::string escape_text;
};

struct Feature {
  std::string type;           // the feature key, e.g. "CDS"
  std::string location_text;  // exactly what goes to `Location.fromstring`
  LocationResult location;
  std::vector<FeatureQualifier> qualifiers;
};

struct FeatureTable {
  // True when this reader declines the table.  `message` then names the line
  // and says what shape it could not reproduce; `features` is empty.
  bool refused = false;
  std::string message;
  std::size_t offset = 0;  // the line the refusal is about
  std::vector<Feature> features;
};

// The line that ends a feature block, per format: GenBank's `SEQUENCE_HEADERS`
// (`ORIGIN`, `CONTIG`, ...) and EMBL's `XX` marker or `SQ`/`CO` lines.  Shared
// between the scan that records where the block is and the reader that walks it,
// because two copies of this rule are two chances for a span to be wrong at one
// end.
bool feature_block_ends(FlatFileFormat format, const char* text,
                        std::size_t len);

// The table between `begin` and `end`, read under the record's own header facts.
// `has_declared_size` is `length is not None` in Python; `circular` and
// `stranded` are the two other arguments `Location.fromstring` takes.
FeatureTable parse_feature_table(FlatFileFormat format, const char* data,
                                 std::size_t size, std::size_t begin,
                                 std::size_t end, std::int64_t declared_size,
                                 bool has_declared_size, bool circular,
                                 bool stranded);

}  // namespace biofasting
