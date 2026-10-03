#include "fastq.hpp"

#include <algorithm>

namespace biofasting {

FastqStatus FastqScanner::next(FastqRecord& out) {
  if (failed_) return FastqStatus::error;
  if (position_ >= size_) return FastqStatus::end;
  if (try_fast_path(out)) return FastqStatus::ok;
  return slow_path(out);
}

bool FastqScanner::try_fast_path(FastqRecord& out) noexcept {
  const Line header = read_line(data_, size_, position_);
  if (!header.exists || header.first() != '@') return false;
  const Line sequence_line = read_line(data_, size_, header.next);
  // A line starting with '+' is the quality marker, so it is not sequence.
  if (!sequence_line.exists || sequence_line.first() == '+') return false;
  const Line plus_line = read_line(data_, size_, sequence_line.next);
  if (!plus_line.exists || plus_line.first() != '+') return false;
  const Line quality_line = read_line(data_, size_, plus_line.next);
  if (!quality_line.exists) return false;

  const std::string_view title = rstrip_view(header.text + 1, header.len - 1);
  const std::string_view sequence =
      rstrip_view(sequence_line.text, sequence_line.len);
  const std::string_view quality =
      rstrip_view(quality_line.text, quality_line.len);

  if (sequence.find(' ') != std::string_view::npos ||
      sequence.find('\t') != std::string_view::npos) {
    return false;
  }
  if (sequence.size() != quality.size()) return false;
  const std::string_view plus_title =
      rstrip_view(plus_line.text + 1, plus_line.len - 1);
  if (!plus_title.empty() && plus_title != title) return false;
  // Biopython reads one line past the quality before it can know the record
  // ended.  If that line is not the next header, the quality was longer than
  // the sequence and this record is malformed -- not a fast-path case.
  if (quality_line.next < size_) {
    const Line following = read_line(data_, size_, quality_line.next);
    if (following.first() != '@') return false;
  }

  out.title = title;
  out.sequence = sequence;
  out.quality = quality;
  position_ = quality_line.next;
  ++record_index_;
  return true;
}

FastqStatus FastqScanner::slow_path(FastqRecord& out) {
  const std::size_t record_start = position_;

  Line line = read_line(data_, size_, position_);
  if (line.first() != '@') {
    fail("Records in Fastq files should start with '@' character", record_start);
    return FastqStatus::error;
  }
  scratch_title_.assign(line.text + 1, rstripped_size(line.text + 1, line.len - 1));
  std::size_t cursor = line.next;

  // One or more sequence lines, up to the '+' that marks the quality.
  scratch_sequence_.clear();
  bool found_plus = false;
  while (cursor < size_) {
    line = read_line(data_, size_, cursor);
    if (line.first() == '+') {
      found_plus = true;
      break;
    }
    append_rstripped(scratch_sequence_, line);
    cursor = line.next;
  }
  if (!found_plus) {
    fail(scratch_sequence_.empty() ? "Unexpected end of file"
                                   : "End of file without quality information.",
         size_);
    return FastqStatus::error;
  }

  // `cursor` still points at the '+' line, and `line` is that line.
  const std::string_view second_title =
      rstrip_view(line.text + 1, line.len - 1);
  if (!second_title.empty() && second_title != scratch_title_) {
    fail("Sequence and quality captions differ.", cursor);
    return FastqStatus::error;
  }
  cursor = line.next;
  if (scratch_sequence_.find(' ') != std::string::npos ||
      scratch_sequence_.find('\t') != std::string::npos) {
    fail("Whitespace is not allowed in the sequence.", record_start);
    return FastqStatus::error;
  }
  const std::size_t seq_len = scratch_sequence_.size();

  // Quality lines, until a line that starts with '@' arrives once the quality
  // is already long enough to be complete.  That rule -- not the character --
  // is how a quality line beginning with '@' survives.
  scratch_quality_.clear();
  bool ran = false;
  while (cursor < size_) {
    line = read_line(data_, size_, cursor);
    ran = true;
    if (line.first() == '@' && scratch_quality_.size() >= seq_len) break;
    append_rstripped(scratch_quality_, line);
    cursor = line.next;
  }
  if (!ran) {
    fail("Unexpected end of file", cursor);
    return FastqStatus::error;
  }
  if (scratch_quality_.size() != seq_len) {
    fail("Lengths of sequence and quality values differs for " + scratch_title_ +
             " (" + std::to_string(seq_len) + " and " +
             std::to_string(scratch_quality_.size()) + ").",
         record_start);
    return FastqStatus::error;
  }

  out.title = scratch_title_;
  out.sequence = scratch_sequence_;
  out.quality = scratch_quality_;
  position_ = cursor;
  ++record_index_;
  return FastqStatus::ok;
}

void FastqScanner::fail(std::string message, std::size_t offset) {
  error_message_ = std::move(message);
  error_offset_ = offset;
  failed_ = true;
}

namespace {

// Whether a field points at bytes of the span rather than into a scanner's
// scratch buffer.  Without this the grid would subtract pointers into two
// different allocations, which is not merely wrong but undefined -- and a file
// that needs any stripping at all lands here, so it is not a corner case.
bool points_into(const char* data, std::size_t size, std::string_view view) noexcept {
  return view.data() >= data && view.data() + view.size() <= data + size;
}

}  // namespace

FastqGrid scan_grid(const char* data, std::size_t size) {
  FastqGrid grid;
  FastqScanner scanner(data, size);
  FastqRecord record;

  // Every field of record 0 is checked against the span before it is used as
  // the origin of anything: a record whose fields were assembled in scratch has
  // no address in this file to be an origin.
  bool established = false;
  std::size_t origin = 0;  // record 0's '@'

  while (true) {
    const FastqStatus status = scanner.next(record);
    if (status == FastqStatus::end) break;
    if (status == FastqStatus::error) {
      grid.reason = "record " + std::to_string(scanner.record_index()) +
                    " is malformed (" + scanner.error_message() + ")";
      return grid;
    }

    if (!established) {
      if (!points_into(data, size, record.title) ||
          !points_into(data, size, record.sequence) ||
          !points_into(data, size, record.quality)) {
        grid.reason =
            "record 0 is not a single-line record with nothing to strip, so "
            "there is no grid to build";
        return grid;
      }
      grid.title_length = record.title.size() + 1;
      grid.length = record.sequence.size();
      grid.sequence_offset = static_cast<std::size_t>(record.sequence.data() - data);
      grid.quality_offset = static_cast<std::size_t>(record.quality.data() - data);
      origin = static_cast<std::size_t>(record.title.data() - data) - 1;
      // The scanner stops just past the record it returned, so where it
      // stopped *is* the stride -- no layout has to be spelled out here, and a
      // CRLF file or an unusual '+' line is handled by being measured rather
      // than by being assumed.
      grid.stride = scanner.position() - origin;
      established = true;
      grid.count = 1;
      continue;
    }

    const std::size_t start =
        origin + grid.count * grid.stride;  // record `count` starts here
    const bool same_shape =
        record.title.size() + 1 == grid.title_length &&
        record.sequence.size() == grid.length &&
        record.quality.size() == grid.length;
    const bool same_place =
        record.title.data() == data + start + 1 &&
        record.sequence.data() == data + start + grid.sequence_offset - origin &&
        record.quality.data() == data + start + grid.quality_offset - origin;
    if (!same_shape || !same_place) {
      grid.reason = "record " + std::to_string(grid.count) +
                    " does not have the same shape as the first one";
      return grid;
    }
    ++grid.count;
  }

  if (grid.count == 0) {
    // An empty file is a grid with no records in it, not an error: the same
    // reading open_fasta() gives an empty file.
    grid.title_length = 0;
    grid.length = 0;
    grid.stride = 0;
    grid.sequence_offset = 0;
    grid.quality_offset = 0;
  }
  return grid;
}

FastqIndex::FastqIndex(const char* data, std::size_t size)
    : data_(data), size_(size) {
  FastqScanner scanner(data, size);
  FastqRecord record;

  while (true) {
    // The scanner stops just past the record it returned, so where it stands
    // *before* a call is where the next record's '@' is.  Nothing here has to
    // know the layout of a FASTQ record, which is why a CRLF file, a '+'-titled
    // record and a wrapped one are all handled by being measured.
    const std::size_t start = scanner.position();
    const FastqStatus status = scanner.next(record);
    if (status == FastqStatus::end) break;
    if (status == FastqStatus::error) {
      // The parser's message plus where it was found: the location is the part
      // Biopython cannot give, and a malformed file of a million records is
      // where it is worth having.
      error_message_ = scanner.error_message() + " (record " +
                       std::to_string(scanner.record_index()) +
                       ", byte offset " +
                       std::to_string(scanner.error_offset()) + ")";
      failed_ = true;
      return;
    }

    FastqEntry entry;
    entry.offset = start;
    entry.end = scanner.position();
    entry.plain = points_into(data, size, record.title) &&
                  points_into(data, size, record.sequence) &&
                  points_into(data, size, record.quality);
    if (entry.plain) {
      // The name is not copied out to be found later: it is recorded as a place
      // inside the title, which is itself a place inside the caller's buffer.
      // So the build allocates nothing per record, and the key the map is
      // hashed on is a view of the file's own bytes.
      const std::string_view title = record.title;
      const std::pair<std::size_t, std::size_t> word =
          first_word_span(title.data(), title.size());
      entry.title_length = static_cast<std::uint32_t>(title.size());
      entry.name = title.data() + word.first;
      entry.name_length = static_cast<std::uint32_t>(word.second);
      entry.sequence_offset =
          static_cast<std::size_t>(record.sequence.data() - data);
      entry.quality_offset =
          static_cast<std::size_t>(record.quality.data() - data);
    } else {
      // Assembled fields have no address in the file, so the key has to be
      // kept; the rest of the record is re-parsed on demand.  A deque, because
      // the entry below is about to hold a pointer into it.
      owned_names_.push_back(first_word(record.title.data(), record.title.size()));
      entry.name = owned_names_.back().data();
      entry.name_length = static_cast<std::uint32_t>(owned_names_.back().size());
    }

    // A read or a header wider than 4 GiB does not fit the entry.  Nothing
    // produces one; if something does, saying so is the only honest answer,
    // because the alternative is an index that answers with a truncated length.
    if (record.title.size() > UINT32_MAX || record.sequence.size() > UINT32_MAX) {
      error_message_ =
          "record " + std::to_string(entries_.size()) +
          " is longer than 4 GiB, which this index cannot address";
      failed_ = true;
      return;
    }
    entry.length = static_cast<std::uint32_t>(record.sequence.size());

    const std::string_view key = key_of(entry);
    if (by_name_.empty()) {
      // One reservation, made from the size of the first record, because a
      // million rehashes of a table that is one short of its next prime is a
      // third of the build: at 1M records the difference is 0.65 s against
      // 0.44 s.  The estimate is only an estimate -- a ragged file can be off
      // by a lot -- so it is a reserve and not a bound, and the table grows if
      // it has to, exactly as it would have without this line.
      const std::size_t span = entry.end - entry.offset;
      if (span > 0) by_name_.reserve(size_ / span + 16);
    }
    if (!by_name_.emplace(key, entries_.size()).second) {
      // Biopython's wording, because it is a documented refusal rather than an
      // implementation detail of its index.
      error_message_ = "Duplicate key '" + std::string(key) + "'";
      failed_ = true;
      return;
    }
    entries_.push_back(entry);
  }
}

const FastqEntry* FastqIndex::find(std::string_view name) const noexcept {
  const auto it = by_name_.find(name);
  if (it == by_name_.end()) return nullptr;
  return &entries_[it->second];
}

FastqIndex::FastqText FastqIndex::rescan(const FastqEntry& entry) const {
  // The span is exactly one record long, so the scanner sees the same bytes it
  // saw during the build and stops at the same place.  It cannot fail: it read
  // this record once already, and the build would have refused the file
  // otherwise.
  FastqScanner scanner(data_ + entry.offset, entry.end - entry.offset);
  FastqRecord record;
  scanner.next(record);
  // Copied out before the scanner dies, because for a record that is not
  // `plain` these views point into its scratch buffer.
  FastqText text;
  text.title.assign(record.title.data(), record.title.size());
  text.sequence.assign(record.sequence.data(), record.sequence.size());
  text.quality.assign(record.quality.data(), record.quality.size());
  return text;
}

std::string FastqIndex::title(const FastqEntry& entry) const {
  if (entry.plain) return std::string(data_ + entry.offset + 1, entry.title_length);
  return rescan(entry).title;
}

std::string FastqIndex::sequence(const FastqEntry& entry) const {
  if (entry.plain) return std::string(data_ + entry.sequence_offset, entry.length);
  return rescan(entry).sequence;
}

std::string FastqIndex::quality(const FastqEntry& entry) const {
  if (entry.plain) return std::string(data_ + entry.quality_offset, entry.length);
  return rescan(entry).quality;
}

std::string FastqIndex::sequence_slice(const FastqEntry& entry,
                                       std::size_t start,
                                       std::size_t end) const {
  start = std::min(start, static_cast<std::size_t>(entry.length));
  end = std::min(end, static_cast<std::size_t>(entry.length));
  if (end <= start) return std::string();
  if (entry.plain) {
    return std::string(data_ + entry.sequence_offset + start, end - start);
  }
  // No stride to walk, so the honest thing is to build the record and cut it.
  // Rare, and correct.
  const std::string whole = sequence(entry);
  return whole.substr(start, end - start);
}

}  // namespace biofasting
