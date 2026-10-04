#include "arrowout.hpp"

#include <stdexcept>
#include <utility>

#include "fasta.hpp"
#include "fastq.hpp"
#include "line.hpp"

namespace biofasting {
namespace {

// The two `title`-derived columns of a record, as views into the title.
//
// The key is `first_word`'s, which is the one rule this package keys a record
// by, and the description is what follows it with the whitespace between them
// (and any more of it) taken off.  Both are views, so building these two columns
// for a million records allocates nothing but the columns themselves.
struct TitleParts {
  std::string_view name;
  std::string_view description;
};

TitleParts split_title(std::string_view title) noexcept {
  const std::pair<std::size_t, std::size_t> key =
      first_word_span(title.data(), title.size());
  std::size_t rest = key.first + key.second;
  while (rest < title.size() && is_python_space(title[rest])) ++rest;
  return {title.substr(key.first, key.second), title.substr(rest)};
}

// The four columns in the order polars-bio's read_fastq produces them.
constexpr const char* kFastqColumns[] = {"name", "description", "sequence",
                                         "quality"};
constexpr const char* kFastaColumns[] = {"name", "description", "sequence"};

}  // namespace

void ArrowStringColumn::reserve(std::size_t records, std::size_t bytes) {
  offsets_.reserve(records + 1);
  data_.reserve(bytes);
}

void ArrowStringColumn::append(std::string_view value) {
  // Guarded because an empty view is allowed to have a null `data()`, and
  // `std::string::append(nullptr, 0)` is undefined rather than a no-op.
  if (!value.empty()) data_.append(value.data(), value.size());
  offsets_.push_back(static_cast<std::uint64_t>(data_.size()));
}

void ArrowStringColumn::seal_row(std::size_t length) {
  const std::uint64_t end = static_cast<std::uint64_t>(data_.size());
  const std::uint64_t start = offsets_.back();
  if (end - start != length) {
    // Unreachable in correct code, and silent corruption if it were not: every
    // offset after this row would be off by the difference and the column would
    // still read back as a column.
    throw std::logic_error(
        "ArrowStringColumn::seal_row: the row wrote " +
        std::to_string(end - start) + " bytes, not " + std::to_string(length));
  }
  offsets_.push_back(end);
}

ArrowTable fastq_table(const char* data, std::size_t size) {
  ArrowTable table;
  for (const char* name : kFastqColumns) table.column_names.emplace_back(name);
  table.columns.resize(table.column_names.size());

  FastqScanner scanner(data, size);
  FastqRecord record;
  bool sized = false;

  for (;;) {
    const FastqStatus status = scanner.next(record);
    if (status == FastqStatus::end) break;
    if (status == FastqStatus::error) {
      table.failed = true;
      table.error_message = located_error(scanner);
      return table;
    }

    const TitleParts parts = split_title(record.title);
    if (!sized) {
      // The one thing one pass cannot know is how big the columns will be.  The
      // first record answers it well enough for the files this exists for --
      // records of one shape -- because the count is the input size over the
      // first record's span and each column's share is its own bytes times that
      // count.  The alternative is a second pass over 369 MB to count bytes
      // that were already counted once, which costs more than every
      // reallocation it would save.
      sized = true;
      const std::size_t span = scanner.position();
      const std::size_t rows = span == 0 ? 0 : size / span;
      table.columns[0].reserve(rows, rows * parts.name.size());
      table.columns[1].reserve(rows, rows * parts.description.size());
      table.columns[2].reserve(rows, rows * record.sequence.size());
      table.columns[3].reserve(rows, rows * record.quality.size());
    }

    table.columns[0].append(parts.name);
    table.columns[1].append(parts.description);
    table.columns[2].append(record.sequence);
    table.columns[3].append(record.quality);
    ++table.records;
  }
  return table;
}

ArrowTable fasta_table(const char* data, std::size_t size) {
  ArrowTable table;
  for (const char* name : kFastaColumns) table.column_names.emplace_back(name);
  table.columns.resize(table.column_names.size());

  const FastaIndex index(data, size);
  if (index.failed()) {
    table.failed = true;
    table.error_message = index.error_message();
    return table;
  }

  const std::vector<FastaEntry>& entries = index.entries();

  // Every byte the three columns will hold is known before one is copied: the
  // index has already walked every header and measured every sequence.  So the
  // buffers are sized exactly and nothing reallocates -- the FASTA rows of the
  // benchmark are the ones whose number is in GB/s rather than in microseconds a
  // record, and this is why.
  std::size_t name_bytes = 0;
  std::size_t description_bytes = 0;
  std::size_t sequence_bytes = 0;
  for (const FastaEntry& entry : entries) {
    const TitleParts parts = split_title(entry.title);
    name_bytes += parts.name.size();
    description_bytes += parts.description.size();
    sequence_bytes += entry.length;
  }
  table.columns[0].reserve(entries.size(), name_bytes);
  table.columns[1].reserve(entries.size(), description_bytes);
  table.columns[2].reserve(entries.size(), sequence_bytes);

  for (const FastaEntry& entry : entries) {
    const TitleParts parts = split_title(entry.title);
    table.columns[0].append(parts.name);
    table.columns[1].append(parts.description);
    append_sequence(table.columns[2].open_row(), data, size, entry);
    table.columns[2].seal_row(entry.length);
  }

  table.records = entries.size();
  return table;
}

}  // namespace biofasting
