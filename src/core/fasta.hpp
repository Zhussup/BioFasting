// FASTA reading over a byte span, with an offset index.
//
// The contract is Bio.SeqIO's, in the two shapes it actually has:
//
//   * the records `SeqIO.parse(handle, "fasta")` yields, which is
//     SimpleFastaParser: everything up to the next '>' is sequence, each line
//     is right-stripped, then spaces are removed from the joined result.  Two
//     consequences that are easy to get wrong and are therefore tested: a ';'
//     comment *inside* a record is sequence text, and a tab inside a line is
//     kept while a space is not;
//   * the keys `SeqIO.index(file, "fasta")` exposes, which are the first word
//     of the title after stripping -- so ">  rec1 desc" is keyed "rec1", and a
//     duplicate key is a ValueError there and must be one here.
//
// What the index adds over the reference is that reading a record does not
// mean parsing a record: the build makes one pass and records, for each entry,
// where its sequence bytes are and what shape they have, so `sequence()` is a
// copy and `sequence_slice()` is a gather.  That is only sound while the
// sequence bytes really are a regular grid, so the pass checks for it: an
// entry is `strided` when every line is the same width, nothing needed
// stripping, and no line held a space.  Entries that are not strided are still
// correct -- `sequence()` falls back to assembling line by line, which is what
// the reference does -- they are simply not O(1) to slice.
//
// Nothing here is a validation pass.  A file that is not FASTA at all has no
// records and is not an error, exactly as the reference treats it.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "line.hpp"

namespace biofasting {

struct FastaEntry {
  std::string name;   // first word of the title: the SeqIO.index key
  std::string title;  // the whole title, right-stripped, as SeqIO reports it
  std::size_t offset = 0;  // first byte of the sequence region
  std::size_t end = 0;     // one past its last byte
  std::size_t length = 0;  // sequence length, after the reference's trimming
  std::uint32_t line_bases = 0;  // .fai column 4
  std::uint32_t line_width = 0;  // .fai column 5, terminator included
  bool strided = false;          // the sequence is a regular grid of lines
};

class FastaIndex {
 public:
  // Builds the index.  Never throws: a duplicate key is reported through
  // failed()/error_message() so that the binding decides which Python
  // exception carries it.
  FastaIndex(const char* data, std::size_t size);

  std::size_t size() const noexcept { return entries_.size(); }
  const std::vector<FastaEntry>& entries() const noexcept { return entries_; }
  const FastaEntry* find(std::string_view name) const noexcept;

  // The record's sequence, byte for byte what SimpleFastaParser would build.
  std::string sequence(const FastaEntry& entry) const;

  // The [start, end) window of it.  Both bounds are clamped to the record, so
  // an out-of-range window is short rather than an error; negative indices are
  // not accepted, because there is no slice object here to make them mean
  // anything.
  std::string sequence_slice(const FastaEntry& entry, std::size_t start,
                             std::size_t end) const;

  bool failed() const noexcept { return failed_; }
  const std::string& error_message() const noexcept { return error_message_; }

 private:
  const char* data_;
  std::size_t size_;
  std::vector<FastaEntry> entries_;
  std::unordered_map<std::string, std::size_t> by_name_;
  bool failed_ = false;
  std::string error_message_;
};

}  // namespace biofasting
