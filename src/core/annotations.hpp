// The annotations dict of the INSDC formats: GenBank's and EMBL's headers.
//
// A flat-file record is a sequence, a features table and a *header* -- and the
// header is the part `SeqIO.parse` puts in `SeqRecord.annotations`.  M21 read
// the features table and left this one open; the measurement that justifies
// building it now is in bench/targets.md (M22), and it is an odd one: a GenBank
// record with an ordinary three-line comment costs the reference 16.4 us more
// than the same record without one, and almost none of that is reading a
// comment.  It is `re.search(r"([^#]+)-START##$", first_line)` -- a probe for a
// structured comment, quadratic in the line it fails on, paid by every record
// that has a comment at all.
//
// So the kernel's job here is the same shape as everywhere else in this package:
// one pass over the header's bytes that decides where each block's value starts
// and ends, and no Python object built until the caller asks.  What the Python
// layer then builds is the dict itself, and the per-format key sets -- which
// differ, and differ in ways a corpus alone cannot show: EMBL has no `date`
// (its scanner never reads a `DT` line), no `source` and no `keywords` key
// unless a `KW` line was there, where GenBank's writer emits placeholders that
// become `''` and `['']`.
//
// Scope, stated rather than discovered later.  This reproduces the header keys
// the corpus and the differential tests exercise, for the two INSDC formats that
// share a scanner.  SwissProt's header is `Bio.SwissProt`, a different reader
// with different keys, and is not here.  A header line this reader does not
// reproduce -- `DBLINK`, `PROJECT`, `NID`, `PID`, `SEGMENT`, `DBSOURCE`, a
// structured comment, a GenBank `LOCUS` line in the pre-229.0 layout -- makes
// the *record's* annotations refuse and name the line, rather than produce a
// dict that is missing a key the reference would have set.  That is the same
// contract the features table has, and for the same reason: a dict that is
// wrong is worse than a dict that is absent.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "flatfile.hpp"

namespace biofasting {

// The keys of `SeqRecord.annotations` this reader produces, in one place so that
// the walk can say "this key was created" without spelling it three times.
//
// The names are the reference's own key names, character for character, and the
// Python layer uses them as the dict's keys directly.  That is not a
// coincidence to be relied on silently: it is the naming decision, made once,
// here and in `annotation_key_name`.
enum class AnnotationKey {
  molecule_type,
  topology,
  data_file_division,
  date,
  accessions,
  sequence_version,
  gi,
  keywords,
  source,
  organism,
  taxonomy,
  references,
  comment,
};

const char* annotation_key_name(AnnotationKey key) noexcept;

// One `REFERENCE` block (GenBank) or `RN` block (EMBL), as the reference's
// `Bio.SeqFeature.Reference` holds it.  `location` is already in Python
// coordinates -- the reference subtracts one from each start and leaves the end
// -- because that conversion is the only arithmetic in the block and doing it
// twice is how two implementations come to disagree.
struct AnnotationReference {
  std::string title;
  std::string authors;
  std::string consrtm;
  std::string journal;
  std::string pubmed_id;
  std::string medline_id;
  std::string comment;
  std::vector<std::pair<std::int64_t, std::int64_t>> location;
};

// The header of one record, as plain values, plus **the order they were
// created in**.
//
// The order is not decoration and not a fixed list.  The reference inserts a key
// when the line that states it is consumed, so a `COMMENT` line above the first
// `REFERENCE` puts `comment` before `references` and one below it puts it after;
// a dict is compared with `==` and printed with `repr`, and both show the
// insertion order.  A fixed order would therefore be right about a corpus and
// wrong about a file.
//
// Membership in `order` is also what says a key *exists*: the reference's
// handlers create some keys unconditionally -- `data_file_division` for an EMBL
// ID line -- and skip others when the line is blank, so "not set" and "set to
// the empty string" are different answers and `order` is what tells them apart.
struct AnnotationTable {
  // From the LOCUS line (GenBank) or the ID line (EMBL).  `date` is empty for a
  // record whose format has no place to put one -- EMBL's scanner reads no date
  // at all -- and for a GenBank LOCUS line whose date field is blank.
  std::string molecule_type;
  std::string topology;
  std::string data_file_division;
  std::string date;

  std::string source;
  std::string organism;
  std::string comment;
  std::vector<std::string> accessions;
  std::vector<std::string> keywords;
  std::vector<std::string> taxonomy;
  std::int64_t sequence_version = 0;
  // The GenBank `VERSION ... GI:n` number, a key of its own.  NCBI stopped
  // assigning GI numbers in 2016, so a modern record has none -- but the files
  // in circulation do, and a reader that refused them would refuse most of what
  // `SeqIO.parse` is pointed at.
  std::string gi;
  std::vector<AnnotationReference> references;

  // Every key the walk created, in the order it created it.  A key appears once
  // however many lines state it: `organism()` and `taxonomy()` extend a list
  // rather than replacing it, and a second `CC` line appends to the first.
  std::vector<AnnotationKey> order;

  // Records that `key` was created.  Idempotent, which is what makes it safe to
  // call from a handler that may run several times for one key.
  void note(AnnotationKey key);

  // Whether `key` has been created yet, which is how a handler that *extends* a
  // key tells its first run from its second.
  bool has(AnnotationKey key) const noexcept;

  bool failed = false;
  std::string error_message;
  std::size_t error_offset = 0;
};

// Reads one record's header from its own bytes.  The span is the record's first
// line up to the line that ends the header -- `FEATURES`/`FH`, or the origin
// line for a record with no feature table -- which the scan has already
// recorded.  `record` supplies the four facts the LOCUS/ID line states.
AnnotationTable read_annotations(FlatFileFormat format, const char* data,
                                 std::size_t size, std::size_t begin,
                                 std::size_t end, const FlatFileRecord& record);

}  // namespace biofasting
