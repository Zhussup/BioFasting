// FASTQ record scanning over a byte span.
//
// The contract is Bio.SeqIO.QualityIO.FastqGeneralIterator, byte for byte.  Not
// because it is pretty, but because it is what the benchmark harness gates
// every row against (bench/run.py compares (title, sequence, quality) and the
// digest is the price of admission), so a kernel that differs from it is a bug
// with a speedup attached.
//
// The parts of that contract that surprise people, all of them load-bearing:
//
//   * a quality line may start with '@' and is still quality -- the sequence
//     length, not the character, decides where a record ends;
//   * sequence and quality may be wrapped over several lines and are
//     concatenated with the line breaks removed;
//   * every field is right-stripped, so a blank line between records vanishes
//     into the concatenation instead of ending anything;
//   * the '+' line may repeat the header, and if it does it must match it;
//   * a space or a tab inside the sequence is an error, not a separator.
//
// Two paths, one meaning.  `next()` tries a fast path first: it recognises the
// pristine case -- four single-line fields, nothing to strip, no '+' title, a
// quality line exactly as long as the sequence -- and returns views straight
// into the caller's buffer, no allocation and no copies.  Anything at all out
// of the ordinary falls back to a faithful transliteration of Biopython's
// algorithm, which is slower, allocates, and gets the error messages and their
// order exactly right.  The fast path never reports an error of its own: if
// anything is off it simply declines, so there is exactly one place where a
// malformed file is diagnosed.
//
// The index below is built by that same scanner, and that is a decision with a
// cost attached, so it is written down here rather than discovered later.
// Biopython has *two* FASTQ grammars: `FastqGeneralIterator`, which is what
// `SeqIO.parse` runs, and `FastqRandomAccess.__iter__`, which is what
// `SeqIO.index` runs.  They disagree.  On `bench/data/edge/blank_lines.fastq` --
// a legal file whose records are separated by a blank line -- the parser reads
// two records and the index raises `ValueError: Problem with line b'\n'`.  This
// package has one grammar, the parser's, because a library whose index and
// parser disagree about which files are readable has two grammars in it; that
// is the bug the reference is demonstrating, not a feature to reproduce.  The
// consequence is deliberate and tested: this index accepts that file, and it
// reports a malformed file with the parser's message and its location rather
// than with the index's four different ones.  `tests/test_fastq_index.py`
// asserts both halves against Biopython so that the divergence stays a
// decision.

#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory_resource>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "line.hpp"

namespace biofasting {

enum class FastqStatus { ok, end, error };

struct FastqRecord {
  std::string_view title;     // header line without '@', right-stripped
  std::string_view sequence;  // concatenated across lines, whitespace removed
  std::string_view quality;   // concatenated likewise; same length as sequence
};

// The geometry of a FASTQ file whose records are all the same shape: same
// header length, same sequence length, at a constant stride.  Such a file is a
// regular grid, and every sequence in it is then a strided view of the file's
// own bytes rather than a copy of them -- which is what makes handing a whole
// 369 MB file to numpy cost nothing.
//
// Deciding this cannot be done by looking at the first record and assuming the
// rest match: that is exactly how a parser silently produces wrong answers on
// the one file in a thousand that is ragged.  So the grid is established by
// *scanning* -- with the same scanner that reads the records, and by requiring
// each record's fields to be the views it already returned, at the addresses
// the geometry predicts.  There is no second notion of "pristine" here to
// disagree with the first.
struct FastqGrid {
  std::size_t count = 0;            // records
  std::size_t title_length = 0;     // header bytes, '@' included
  std::size_t length = 0;           // sequence length; the quality length too
  std::size_t stride = 0;           // bytes from one record's '@' to the next
  std::size_t sequence_offset = 0;  // record 0's sequence, from the start of the span
  std::size_t quality_offset = 0;   // record 0's quality, likewise
  std::string reason;               // empty exactly when the file is a grid

  bool ok() const noexcept { return reason.empty(); }
};

// Establishes the grid, or says why the bytes are not one.  Scans the whole
// span: a grid that is only mostly regular is not a grid.
FastqGrid scan_grid(const char* data, std::size_t size);

class FastqScanner {
 public:
  FastqScanner(const char* data, std::size_t size) noexcept
      : data_(data), size_(size) {}

  // The span being scanned.  Needed by the grid, which builds views over the
  // same bytes: for a gzipped file those are the inflated ones, and only the
  // scanner knows where they are.
  const char* data() const noexcept { return data_; }
  std::size_t size() const noexcept { return size_; }

  // Reads the next record.  `end` means clean EOF; `error` means the input is
  // malformed, and error_message()/error_offset() say where.
  FastqStatus next(FastqRecord& out);

  std::size_t record_index() const noexcept { return record_index_; }
  std::size_t error_offset() const noexcept { return error_offset_; }
  const std::string& error_message() const noexcept { return error_message_; }

  // Bytes consumed so far: the offset just past the last record returned.  The
  // benchmark uses it to report MB/s without re-deriving the input size.
  std::size_t position() const noexcept { return position_; }

 private:
  // Returns false when the record is not the pristine four-line shape; leaves
  // position_ untouched.
  bool try_fast_path(FastqRecord& out) noexcept;
  FastqStatus slow_path(FastqRecord& out);

  void fail(std::string message, std::size_t offset);

  const char* data_;
  std::size_t size_;
  std::size_t position_ = 0;
  std::size_t record_index_ = 0;
  std::size_t error_offset_ = 0;
  bool failed_ = false;  // once malformed, always malformed
  std::string error_message_;

  // Where the slow path assembles fields it cannot point at directly because
  // something had to be stripped or joined.
  std::string scratch_title_;
  std::string scratch_sequence_;
  std::string scratch_quality_;
  std::string scratch_other_title_;
};

// A scanner's failure as one string: the parser's own message plus where it was
// found, which the parser alone cannot say and a malformed file of a million
// records makes worth having.
//
// Free, and not a member, because two callers report the same failure -- the
// scanner's own binding and the Arrow table builder -- and two copies of the
// wording would be two ways for one malformed byte to read.
std::string located_error(const FastqScanner& scanner);

// One record, as a place in the file rather than as bytes.
//
// `plain` is the whole design in one bool.  A record whose four fields are
// single lines with nothing stripped has every byte of itself at a fixed
// address, so a fetch is a memcpy from a stored offset and nothing is copied at
// build time.  A record the scanner had to assemble -- a wrapped sequence, a
// blank line, a '+'-repeated header, anything with trailing whitespace -- has no
// such address, and is fetched by re-scanning its own span, which is what the
// reference does for *every* record.
//
// The distinction never changes an answer.  Both paths return the bytes
// `FastqGeneralIterator` would return; they differ only in what a fetch costs.
struct FastqEntry {
  std::size_t offset = 0;           // the record's '@'
  std::size_t end = 0;              // one past the byte that ends it
  std::size_t sequence_offset = 0;  // valid when `plain`
  std::size_t quality_offset = 0;   // valid when `plain`
  const char* name = nullptr;       // the key's first byte: in the buffer, or owned
  std::uint32_t title_length = 0;   // header bytes after the '@', stripped
  std::uint32_t name_length = 0;    // how long the key is
  std::uint32_t length = 0;         // sequence length; the quality's length too
  bool plain = false;
};

// A name-keyed index over a FASTQ file held in memory.
//
// The build is one pass with the scanner above, so the set of files this
// accepts is exactly the set `FastqGeneralIterator` accepts -- see the note at
// the top of this file for why that is a decision and not an accident.  What it
// adds is that a lookup afterwards reads the record's bytes instead of parsing
// them again, which is the difference between 7 µs and 0.2 µs per fetch.
//
// Keys are the reference's: the title's first word, in file order, with a
// duplicate key refused.  Refused through failed()/error_message() rather than
// by throwing, so that the binding owns the exception type -- the same split
// `FastaIndex` makes.
class FastqIndex {
 public:
  FastqIndex(const char* data, std::size_t size);

  std::size_t size() const noexcept { return entries_.size(); }
  const std::vector<FastqEntry>& entries() const noexcept { return entries_; }
  const FastqEntry* find(std::string_view name) const noexcept;

  // The entry's key, as a view: into the caller's buffer for a record the
  // scanner could view, into this index's own store for one it had to assemble.
  // Either way it outlives the index, and that is the whole reason the lookup
  // table can hold views instead of a million copies of a million names.
  std::string_view key_of(const FastqEntry& entry) const noexcept {
    return std::string_view(entry.name, entry.name_length);
  }

  // The record's title, right-stripped, without its '@'.
  std::string title(const FastqEntry& entry) const;
  // The record's sequence, byte for byte what the parser would build.
  std::string sequence(const FastqEntry& entry) const;
  std::string quality(const FastqEntry& entry) const;
  // The [start, end) window of the sequence, clamped to the record.  Negative
  // indices mean nothing here, exactly as in FastaIndex.
  std::string sequence_slice(const FastqEntry& entry, std::size_t start,
                             std::size_t end) const;

  bool failed() const noexcept { return failed_; }
  const std::string& error_message() const noexcept { return error_message_; }

 private:
  // The record's fields, parsed again from their own span.  Only reachable for
  // a record that is not `plain`, which is why it is allowed to be slow.
  //
  // It returns owned strings rather than the `FastqRecord` the scanner hands
  // back, and that is not a style choice: a record the scanner had to assemble
  // is returned as views into the scanner's own scratch buffer, so the record
  // would be dangling the moment the scanner went out of scope.  Copying here
  // is the fix, and the copy is one the reference makes for every record.
  struct FastqText {
    std::string title;
    std::string sequence;
    std::string quality;
  };
  FastqText rescan(const FastqEntry& entry) const;

  const char* data_;
  std::size_t size_;
  std::vector<FastqEntry> entries_;
  // Keyed by view, not by string: a key either sits in the caller's buffer or
  // in `owned_names_` below, and both outlive the index, so the map needs no
  // copy of a name.  A million 60-byte headers is a million allocations the
  // build does not make.
  //
  // The nodes come out of one arena rather than one `new` each, and that is
  // measured rather than assumed: with a per-node allocator the table costs
  // 0.26 s of a 0.49 s build -- 53% of it -- against 0.01 s of bookkeeping for
  // the entries themselves.  The arena is declared first because the map holds
  // a pointer to it, and it is never released early: it lives exactly as long
  // as the index, which is the only thing that reads it.
  std::pmr::monotonic_buffer_resource arena_;
  std::pmr::unordered_map<std::string_view, std::size_t> by_name_{&arena_};
  // Names that have no address in the buffer because the scanner assembled
  // them.  A deque, because its elements do not move when it grows and the map
  // above is holding views of them.
  std::deque<std::string> owned_names_;
  bool failed_ = false;
  std::string error_message_;
};

}  // namespace biofasting
