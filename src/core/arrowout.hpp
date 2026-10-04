// Record tables in Arrow's layout, built where the records are.
//
// The shape of the output is `polars-bio`'s `read_fastq` / `read_fasta`, and
// that is a decision about what this is *for*.  A table can be made from this
// package today -- `pa.table` over `open_fastq` -- one Python call per record;
// what that costs was measured before any of this was written (`bench/
// rank_arrow.py`), and what it costs is the whole gap to the tool in the niche.
// So the table is built in C++, in one pass, into the buffers Arrow itself
// would have used, and the Python side only hands those buffers over.  The
// columns are therefore the other tool's, in its order:
//
//     FASTQ   name, description, sequence, quality
//     FASTA   name, description, sequence
//
// `name` is the reference's key -- the title's first word, the same rule the
// two indexes use -- and the sequence is the reference's own bytes.  The one
// column that is neither is `description`, and it is worth stating exactly
// because it is easy to assume otherwise: **the title after the key, with the
// whitespace between them trimmed**.  `>rec1 a description` gives `a
// description` and `>rec1` gives the empty string.  That is polars-bio's
// reading, and it is *not* `SeqRecord.description`, which is the whole title,
// key included.  The two libraries also disagree about the empty case -- see
// `read_*_table` in biofasting.arrow, which says which way this goes and why.
//
// The layout is Arrow's own, not something Arrow can be handed afterwards:
//
//   * a column is an offsets buffer of n+1 little-endian uint64 and one data
//     buffer holding every value concatenated, with nothing between them --
//     no terminator, no per-value allocation, so a million records cost one
//     allocation and not a million;
//   * the offsets are 64-bit because the values are `large_string`, not
//     `string`.  A 2 GB column is not a hypothetical here: the corpus alone has
//     a 40 Mbp chromosome, and a caller who concatenates a few of those is past
//     the 32-bit limit without doing anything unusual;
//   * a column of bytes may hold any byte, including a NUL, so the data is a
//     `std::string` and never a C string.
//
// What the build does *not* do is keep the file: a table cannot be a view of a
// FASTQ file, because a column has to be contiguous and a file's records are
// not.  Every byte is copied exactly once, into buffers this table owns.  The
// FASTA sequence column is the one place where the copy is not a memcpy of a
// contiguous span -- a sequence's bytes have newlines between the lines -- and
// it is de-columned by the same function the FASTA reader uses rather than by a
// second implementation that could disagree with it.
//
// Nothing here validates.  A malformed FASTQ file fails the way the parser
// fails, with the parser's message and its location; a FASTA file with a
// duplicate key fails the way the index fails; a file that is not FASTA at all
// has no records and is not an error, exactly as the reference treats it.  Like
// the two indexes, a failure is reported through the struct rather than thrown,
// so that the binding owns the choice of Python exception.

#ifndef BIOFASTING_ARROWOUT_HPP
#define BIOFASTING_ARROWOUT_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace biofasting {

// One string column in Arrow's layout, being built.
//
// `offsets` has one entry per value plus a leading zero, and its last entry is
// always the size of `data`: the two are written together, by `append` for a
// value that is already contiguous and by `open_row`/`seal_row` for one that is
// assembled in place.
class ArrowStringColumn {
 public:
  // Enough room for `records` values totalling `bytes`.  Both are estimates:
  // being wrong costs a reallocation, never an answer.
  void reserve(std::size_t records, std::size_t bytes);

  // One value, copied from wherever it already is.
  void append(std::string_view value);
  void append(const char* text, std::size_t size) {
    append(std::string_view(text, size));
  }

  // The pair for a value whose bytes are not contiguous -- a FASTA sequence,
  // whose newlines are not part of it.  `open_row` hands over the data buffer
  // to be appended to, `seal_row` says how many bytes the row was supposed to
  // be and adds the offset.
  //
  // The split exists because the two halves of a row are known by different
  // code: `append_sequence` knows the bytes and never sees an offset, the
  // column knows the offsets and not the bytes.  A row sealed with a length it
  // did not write would be a column whose offsets disagree with its data --
  // every value after it shifted, and no reader that could tell -- so the
  // length is checked rather than trusted, and the two calls are made in one
  // place (`fasta_table`) and nowhere else.
  std::string& open_row() noexcept { return data_; }
  void seal_row(std::size_t length);

  const std::uint64_t* offsets() const noexcept { return offsets_.data(); }
  const char* data() const noexcept { return data_.data(); }
  std::size_t bytes() const noexcept { return data_.size(); }
  std::size_t size() const noexcept { return offsets_.size() - 1; }  // values

 private:
  std::vector<std::uint64_t> offsets_{0};
  std::string data_;
};

// A whole table: the column names, the columns, and how many rows they hold.
struct ArrowTable {
  std::vector<std::string> column_names;
  std::vector<ArrowStringColumn> columns;
  std::size_t records = 0;

  bool failed = false;
  std::string error_message;

  bool ok() const noexcept { return !failed; }
};

// The four FASTQ columns, in one pass with the fastq.hpp scanner -- so the set
// of files this accepts is exactly the set `FastqGeneralIterator` accepts, and
// a malformed record lands in `error_message` with its location.
ArrowTable fastq_table(const char* data, std::size_t size);

// The three FASTA columns, from the fasta.hpp index: the headers are already
// walked and every sequence already measured by the time the first byte is
// copied, so all three buffers are sized exactly and nothing reallocates.
ArrowTable fasta_table(const char* data, std::size_t size);

}  // namespace biofasting

#endif  // BIOFASTING_ARROWOUT_HPP
