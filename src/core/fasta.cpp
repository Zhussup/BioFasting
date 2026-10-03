#include "fasta.hpp"

#include <algorithm>
#include <cstring>
#include <utility>

namespace biofasting {
namespace {

// The title's first word, which is what SeqIO.index keys a record by, lives in
// line.cpp: the FASTQ index keys records by the same rule, and two copies of it
// would be two answers to one question.

std::size_t count_spaces(const char* text, std::size_t size) {
  std::size_t count = 0;
  for (std::size_t i = 0; i < size; ++i) {
    if (text[i] == ' ') ++count;
  }
  return count;
}

// Appends the first `size` bytes of `text` with the spaces taken out, which is
// how the reference joins a record's lines: rstrip each, concatenate, then
// `replace(" ", "")`.  Removing them per line gives the same bytes as removing
// them after joining, and needs no second buffer.
void append_without_spaces(std::string& out, const char* text,
                           std::size_t size) {
  const char* const end = text + size;
  const char* chunk = text;
  while (chunk < end) {
    const char* space =
        static_cast<const char*>(std::memchr(chunk, ' ', end - chunk));
    if (space == nullptr) {
      out.append(chunk, static_cast<std::size_t>(end - chunk));
      return;
    }
    out.append(chunk, static_cast<std::size_t>(space - chunk));
    chunk = space + 1;
  }
}

}  // namespace

FastaIndex::FastaIndex(const char* data, std::size_t size)
    : data_(data), size_(size) {
  std::size_t cursor = 0;

  // Skip anything before the first record.  A blank line, a ';' comment and a
  // file that is not FASTA at all are the same thing here: not a record.
  Line line = read_line(data_, size_, cursor);
  while (line.exists && line.first() != '>') {
    cursor = line.next;
    line = read_line(data_, size_, cursor);
  }

  while (line.exists) {
    FastaEntry entry;
    entry.title.assign(line.text + 1,
                       rstripped_size(line.text + 1, line.len - 1));
    entry.name = first_word(line.text + 1, line.len - 1);
    entry.offset = line.next;

    // One pass over the record's lines, doing two jobs at once: summing the
    // sequence length the reference would produce, and finding out whether the
    // bytes are a regular grid that `sequence()` can copy out in bulk.
    //
    // The last line is the trap.  A wrapped record almost always ends on a
    // short line -- 40,000,000 bases at 60 per line ends with 40 -- and that is
    // not a broken grid, it is the tail the bulk copy expects.  So a line is
    // only compared against the grid once a *later* line proves it was not the
    // last one, and the count of lines is checked against the length at the
    // end, which is what catches a last line that is too long for the tail.
    bool strided = true;
    bool seen_a_line = false;
    std::size_t line_count = 0;
    std::size_t previous_length = 0;
    std::size_t previous_width = 0;
    std::size_t sequence_cursor = line.next;
    Line sequence_line = read_line(data_, size_, sequence_cursor);
    while (sequence_line.exists && sequence_line.first() != '>') {
      const std::size_t kept = rstripped_size(sequence_line.text,
                                              sequence_line.len);
      const std::size_t spaces = count_spaces(sequence_line.text, kept);
      const std::size_t width = sequence_line.next - sequence_cursor;
      if (seen_a_line) {
        // `previous_*` is now known not to be the last line of the record.
        if (previous_length != entry.line_bases ||
            previous_width != entry.line_width) {
          strided = false;
        }
      } else {
        entry.line_bases = static_cast<std::uint32_t>(sequence_line.len);
        entry.line_width = static_cast<std::uint32_t>(width);
        seen_a_line = true;
      }
      previous_length = sequence_line.len;
      previous_width = width;
      ++line_count;

      // Trailing whitespace or an interior space both mean this line is not a
      // plain copy of its bytes, and one such line is enough to give up the
      // stride for the whole record.
      if (kept != sequence_line.len || spaces != 0) strided = false;
      entry.length += kept - spaces;

      sequence_cursor = sequence_line.next;
      sequence_line = read_line(data_, size_, sequence_cursor);
    }

    if (entry.line_bases > 0) {
      const std::size_t full_lines = entry.length / entry.line_bases;
      const std::size_t tail = entry.length % entry.line_bases;
      const std::size_t expected_lines = full_lines + (tail == 0 ? 0 : 1);
      if (line_count != expected_lines) strided = false;
    }
    entry.end = sequence_cursor;
    // line_bases/line_width stay as measured even when the grid is broken:
    // they are still the .fai columns for this record, and `sequence()` keys
    // off `strided`, not off them being zero.
    entry.strided = strided && line_count > 0 && entry.line_bases > 0;

    if (!by_name_.emplace(entry.name, entries_.size()).second) {
      // Biopython's wording, because it is a documented refusal rather than an
      // implementation detail of its index.
      error_message_ = "Duplicate key '" + entry.name + "'";
      failed_ = true;
      return;
    }
    entries_.push_back(std::move(entry));

    line = sequence_line;  // either the next header or a line that is not there
  }
}

const FastaEntry* FastaIndex::find(std::string_view name) const noexcept {
  const auto it = by_name_.find(std::string(name));
  if (it == by_name_.end()) return nullptr;
  return &entries_[it->second];
}

std::string FastaIndex::sequence(const FastaEntry& entry) const {
  std::string out;
  out.reserve(entry.length);

  if (entry.strided) {
    const std::size_t bases = entry.line_bases;
    const std::size_t full_lines = entry.length / bases;
    const std::size_t tail = entry.length % bases;
    const char* line = data_ + entry.offset;
    for (std::size_t i = 0; i < full_lines; ++i) {
      out.append(line, bases);
      line += entry.line_width;
    }
    out.append(line, tail);
    return out;
  }

  std::size_t cursor = entry.offset;
  while (cursor < entry.end) {
    const Line line = read_line(data_, size_, cursor);
    append_without_spaces(out, line.text, rstripped_size(line.text, line.len));
    cursor = line.next;
  }
  return out;
}

std::string FastaIndex::sequence_slice(const FastaEntry& entry,
                                       std::size_t start,
                                       std::size_t end) const {
  start = std::min(start, entry.length);
  end = std::min(end, entry.length);
  if (end <= start) return std::string();

  if (!entry.strided) {
    // The grid is not there to walk, so the honest thing is to build the whole
    // sequence and cut it.  Rare, and correct.
    std::string whole = sequence(entry);
    return whole.substr(start, end - start);
  }

  std::string out;
  out.resize(end - start);
  char* destination = out.data();
  std::size_t position = start;
  std::size_t remaining = end - start;
  while (remaining > 0) {
    const std::size_t line = position / entry.line_bases;
    const std::size_t column = position % entry.line_bases;
    const std::size_t take =
        std::min(remaining, static_cast<std::size_t>(entry.line_bases) - column);
    std::memcpy(destination, data_ + entry.offset + line * entry.line_width + column,
                take);
    destination += take;
    position += take;
    remaining -= take;
  }
  return out;
}

}  // namespace biofasting
