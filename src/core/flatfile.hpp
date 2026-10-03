// Flat-file readers for the INSDC formats: GenBank, EMBL and SwissProt.
//
// These three carry essentially all of the world's sequence data, and none of
// them is a line-oriented format in the FASTA sense: a record is a *block* of
// keyworded lines ended by `//`, and the residues sit inside a numbered block
// that has to be de-columned before it is a sequence.
//
// The measurement that justifies this file is in bench/targets.md (sixth
// ranking pass).  The reference spends 28.4 us on a 1 kb GenBank record with no
// features, where moving the bytes costs 1.5 us and building the object tree
// costs 1.0 us.  So the kernel's job is precisely the part that is neither of
// those: one pass over the buffer that records where each record's fields are,
// and a de-columning copy when a caller asks for the residues.
//
// Scope, stated rather than discovered later: this reproduces `id`, `name`,
// `description`, `seq` and the FEATURES table (`FlatFileIndex::features`).  It
// does **not** reproduce the annotations dict beyond the few header facts a
// feature location is parsed under, nor the references, so it is not a drop-in
// for `SeqIO.parse` and must not be presented as one.
//
// The other half of the contract is what it *refuses*.  Biopython reacts to a
// malformed sequence line by warning and guessing -- it shifts a wrongly
// indented GenBank line by one byte and carries on.  A kernel that guessed
// differently would produce a sequence that looks right and is not, which is
// the one failure mode this project must not ship.  So a line shaped in a way
// this reader does not reproduce makes the whole index fail with a message
// naming it, and the Python layer is free to fall back to `SeqIO.parse`.
//
// The feature table is read under the same rule, at the granularity of one
// record rather than of the file: `features()` refuses the record it cannot
// reproduce and says which line did it, leaving the sequences it did read
// usable.  The malformed feature shapes that make it refuse are listed in
// feature.hpp.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace biofasting {

enum class FlatFileFormat { genbank, embl, swiss };

// Defined in feature.hpp, which is the header that knows what a feature is.
// Forward-declared here so that the index can hand one back without this header
// having to pull in the whole feature grammar.
struct FeatureTable;

// The same arrangement for the annotations reader: `annotations.cpp` owns the
// header grammar, and this header only has to be able to name what it returns.
struct AnnotationTable;

// The name `SeqIO.parse` knows the format by, which is also what a caller of
// this package passes.  Here rather than in the binding so that an error
// message can name the format without the binding having to.
const char* flatfile_format_name(FlatFileFormat format) noexcept;

struct FlatFileRecord {
  std::string id;           // SeqRecord.id
  std::string name;         // SeqRecord.name
  std::string description;  // SeqRecord.description
  std::size_t sequence_offset = 0;  // first byte of the residue block
  std::size_t sequence_end = 0;     // one past its last byte
  std::size_t length = 0;           // residues, counted during the scan
  std::size_t offset = 0;           // the record's first line, for error reports

  // The FEATURES/FT block as a span of the file: the first line of the block and
  // one past its last.  Equal for a record with no feature table, which is why
  // the pair is a span and not an "exists" flag -- an empty span is already the
  // statement.  Like the residues, it is read on demand, so a caller that only
  // wants sequences pays nothing for features.
  std::size_t features_offset = 0;
  std::size_t features_end = 0;

  // The header, as a span: the record's first line up to the line that ends the
  // header -- `FEATURES`/`FH`, or the origin line for a record with no feature
  // table.  The annotations dict is read out of exactly this span, on demand.
  std::size_t header_end = 0;

  // The four facts the LOCUS line (GenBank) or the ID line (EMBL) states and the
  // annotations dict carries: `molecule_type`, `topology`, `data_file_division`
  // and -- for GenBank only, since EMBL's scanner reads no date at all -- `date`.
  // Empty when the line did not say, which for `date` is also what a reader of a
  // record that has no date should report.
  std::string molecule_type;
  std::string topology;
  std::string data_file_division;
  std::string date;
  // The accession the *first* line states.  EMBL's ID line carries one and its
  // consumer records it there, before the topology and the molecule type, which
  // is why an EMBL record's annotations begin `accessions, sequence_version,
  // topology, ...`.  GenBank's LOCUS line states none -- its accession arrives
  // later, on `ACCESSION` -- so this is empty for a GenBank record and the walk
  // creates the key where the line is instead.
  std::string first_line_accession;
  // The version suffix the first line stated (`ID   X56734; SV 2; ...`), and
  // whether there was one.  The reference's first-line handler calls
  // `version_suffix` for it, and that handler is what creates the
  // `sequence_version` annotation key -- which is a different statement from the
  // `X56734.2` the record's id becomes, and a record with no suffix has neither.
  std::int64_t first_line_sequence_version = 0;
  bool has_first_line_sequence_version = false;
  // True when those facts were read from a layout this reader reproduces.  False
  // for a GenBank LOCUS line in the pre-229.0 layout, whose columns are laid out
  // differently: the features table may still be readable from such a record,
  // but its annotations are refused rather than built from columns that were
  // never there.
  bool has_annotation_header = false;

  // What the reference's feature consumer takes from the header, and therefore
  // what its `Location.fromstring` call is given: the size the record *declares*
  // (`_expected_size`), whether the record is circular, and whether it is
  // stranded.  `has_declared_size` is `length is not None` in Python, which is
  // not the same fact as a size of zero -- the origin-wrapping repair is skipped
  // for a record with no declared size and taken for one that declares zero.
  std::int64_t declared_size = 0;
  bool has_declared_size = false;
  bool circular = false;
  bool stranded = true;  // false only for a record whose residue type is protein
  // False when the header line is not in a layout this reader reproduces, which
  // leaves the three facts above unset.  A record with a feature table and no
  // header facts refuses its features instead of reading their locations under
  // inputs the reference never used.
  bool has_header_facts = false;
};

// One pass over a whole flat file, with the residues answered by offset
// afterwards.  Keys are record ids, as everywhere else in this package.
class FlatFileIndex {
 public:
  FlatFileIndex(const char* data, std::size_t size, FlatFileFormat format);

  FlatFileFormat format() const noexcept { return format_; }
  std::size_t size() const noexcept { return records_.size(); }
  const std::vector<FlatFileRecord>& records() const noexcept { return records_; }

  const FlatFileRecord* find(std::string_view id) const noexcept;

  // Bytes, not a record: turning these into a SeqRecord is the interop shim's
  // job, and doing it here would put a Biopython dependency in the core.
  std::string sequence(const FlatFileRecord& record) const;
  std::string sequence_slice(const FlatFileRecord& record, std::size_t start,
                             std::size_t end) const;

  // The record's FEATURES table: feature keys, parsed locations and qualifiers,
  // or a refusal naming the line that made this reader decline.  It is a method
  // on the index rather than a free function because a feature's location is
  // read under the record's own declared size, topology and residue type, and
  // those are facts the scan has already established.
  FeatureTable features(const FlatFileRecord& record) const;

  // The record's header: the keys of `SeqRecord.annotations`, as plain values,
  // or a refusal naming the line this reader does not reproduce.  Like the
  // features table it is read on demand out of the span the scan recorded, so a
  // caller that wants sequences pays nothing for a header it never asks about.
  AnnotationTable annotations(const FlatFileRecord& record) const;

  bool failed() const noexcept { return failed_; }
  const std::string& error_message() const noexcept { return error_message_; }

 private:
  void fail(std::string message, std::size_t offset);

  const char* data_;
  std::size_t size_;
  FlatFileFormat format_;
  std::vector<FlatFileRecord> records_;
  std::unordered_map<std::string, std::size_t> by_id_;
  bool failed_ = false;
  std::string error_message_;
};

}  // namespace biofasting
