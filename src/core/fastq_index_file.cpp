// The file-backed FASTQ index: writing it, validating it, reading through it.
//
// See the header for what is stored and why.  The three parts of this file are
// deliberately separate and touch each other at one point only -- the entry
// encoding -- because they fail in different ways and a bug in one should not
// be able to pass for a bug in another:
//
//   * `write_fastq_index_file` turns the build's heap entries into disk ones,
//     and refuses rather than truncating anything it cannot encode;
//   * `FastqIndexView::open` checks the file against the data it claims to
//     describe, and says which check failed;
//   * `find` never checks anything, because `open` already did, and a lookup
//     that re-validated would cost more than the lookup it is meant to save.

#include "fastq_index_file.hpp"

#include <bit>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

namespace biofasting {
namespace {

std::size_t next_power_of_two(std::size_t value) noexcept {
  std::size_t n = 8;
  while (n < value) n <<= 1;
  return n;
}

bool little_endian(std::string& error) {
  if constexpr (std::endian::native == std::endian::little) {
    return true;
  } else {
    error = "index files are little-endian and this host is not";
    return false;
  }
}

}  // namespace

std::uint64_t fastq_data_fingerprint(const char* data,
                                     std::size_t size) noexcept {
  const std::size_t page = 4096;
  const std::size_t head = size < page ? size : page;
  const std::size_t tail = size < page ? 0 : page;
  const std::uint64_t head_hash = fastq_key_hash(data, head);
  const std::uint64_t tail_hash =
      tail == 0 ? 0 : fastq_key_hash(data + size - tail, tail);
  // Mixed rather than concatenated: two files that differ only in the middle
  // are both caught by the size, and one that was replaced with a file of the
  // same size is caught by the ends.
  return head_hash ^ (tail_hash * 0x9e3779b97f4a7c15ULL) ^
         (static_cast<std::uint64_t>(size) * 0xc2b2ae3d27d4eb4fULL);
}

bool write_fastq_index_file(const std::string& path, const FastqIndex& index,
                            const char* data, std::size_t size,
                            const char* origin, std::size_t origin_size,
                            std::string& error) {
  if (!little_endian(error)) return false;

  const std::vector<FastqEntry>& entries = index.entries();
  const std::size_t count = entries.size();
  const std::size_t slot_count = next_power_of_two(count * 2 + 1);

  std::vector<FastqEntryDisk> disk;
  disk.reserve(count);
  std::string names;
  std::vector<std::uint32_t> slots(slot_count, 0xFFFFFFFFu);

  for (std::size_t i = 0; i < count; ++i) {
    const FastqEntry& entry = entries[i];
    FastqEntryDisk row{};
    row.offset = entry.offset;
    row.end_delta = static_cast<std::uint32_t>(entry.end - entry.offset);
    row.title_length = static_cast<std::uint16_t>(entry.title_length);
    row.name_length = static_cast<std::uint16_t>(entry.name_length);
    row.length = entry.length;
    row.plain = entry.plain ? 1 : 0;

    if (entry.title_length > std::numeric_limits<std::uint16_t>::max() ||
        entry.name_length > std::numeric_limits<std::uint16_t>::max() ||
        entry.end - entry.offset >
            std::numeric_limits<std::uint32_t>::max()) {
      error = "record " + std::to_string(i) +
              " is too long to encode in an index file";
      return false;
    }

    if (entry.plain) {
      const std::ptrdiff_t name_in_title =
          entry.name - (data + entry.offset + 1);
      if (name_in_title < 0 ||
          static_cast<std::size_t>(name_in_title) > entry.title_length) {
        error = "record " + std::to_string(i) +
                " has a key outside its own header";
        return false;
      }
      row.sequence_delta =
          static_cast<std::uint32_t>(entry.sequence_offset - entry.offset);
      row.quality_delta =
          static_cast<std::uint32_t>(entry.quality_offset - entry.offset);
      row.name_offset = static_cast<std::uint32_t>(name_in_title);
    } else {
      row.name_offset = static_cast<std::uint32_t>(names.size());
      names.append(entry.name, entry.name_length);
    }
    if (names.size() > std::numeric_limits<std::uint32_t>::max()) {
      error = "the side table of assembled names is too large to encode";
      return false;
    }
    disk.push_back(row);

    // The same table the build uses, laid out flat: hash, linear probe.  Built
    // here from the entries rather than copied from the build's map, because
    // the map is keyed on pointers and the file needs indices.
    std::uint64_t h = (entry.key_hash != 0
                           ? entry.key_hash
                           : fastq_key_hash(entry.name, entry.name_length)) &
                      (slot_count - 1);
    while (slots[h] != 0xFFFFFFFFu) h = (h + 1) & (slot_count - 1);
    slots[h] = static_cast<std::uint32_t>(i);
  }

  const std::uint64_t entries_offset = sizeof(FastqIndexHeader);
  const std::uint64_t names_offset = entries_offset + count * sizeof(FastqEntryDisk);
  // Padded to a multiple of eight: the lookup table is read as `uint32_t`
  // through a pointer into the mapping, and a misaligned one is undefined
  // behaviour that happens to work on x86 and faults on some ARM cores.
  const std::uint64_t slots_offset = (names_offset + names.size() + 7) & ~7ULL;
  const std::string padding(
      static_cast<std::size_t>(slots_offset - (names_offset + names.size())),
      '\0');

  FastqIndexHeader header{};
  std::memcpy(header.magic, kFastqIndexMagic, sizeof(header.magic));
  header.version = kFastqIndexVersion;
  header.header_size = sizeof(FastqIndexHeader);
  header.data_size = size;
  header.data_fingerprint = fastq_data_fingerprint(data, size);
  header.record_count = count;
  header.entries_offset = entries_offset;
  header.names_offset = names_offset;
  header.names_size = names.size();
  header.slots_offset = slots_offset;
  header.slots_count = slot_count;
  header.source_size = origin == nullptr ? 0 : origin_size;
  header.source_fingerprint =
      origin == nullptr ? 0 : fastq_data_fingerprint(origin, origin_size);

  const std::string temporary = path + ".tmp";
  std::FILE* file = std::fopen(temporary.c_str(), "wb");
  if (file == nullptr) {
    error = "cannot open " + temporary + " for writing";
    return false;
  }
  bool ok = std::fwrite(&header, sizeof(header), 1, file) == 1;
  ok = ok && (disk.empty() ||
              std::fwrite(disk.data(), sizeof(FastqEntryDisk), disk.size(),
                          file) == disk.size());
  ok = ok && (names.empty() ||
              std::fwrite(names.data(), 1, names.size(), file) == names.size());
  ok = ok && (padding.empty() ||
              std::fwrite(padding.data(), 1, padding.size(), file) ==
                  padding.size());
  ok = ok && std::fwrite(slots.data(), sizeof(std::uint32_t), slots.size(),
                         file) == slots.size();
  ok = std::fflush(file) == 0 && ok;
  ok = std::fclose(file) == 0 && ok;
  if (!ok) {
    std::remove(temporary.c_str());
    error = "cannot write " + temporary;
    return false;
  }
  if (std::rename(temporary.c_str(), path.c_str()) != 0) {
    std::remove(temporary.c_str());
    error = "cannot replace " + path;
    return false;
  }
  return true;
}

bool FastqIndexView::open(const char* index_bytes, std::size_t index_size,
                          const char* data, std::size_t data_size,
                          const char* origin, std::size_t origin_size,
                          FastqIndexView& out, std::string& error) {
  if (!little_endian(error)) return false;
  if (index_size < sizeof(FastqIndexHeader)) {
    error = "index file is shorter than its header";
    return false;
  }
  FastqIndexHeader header{};
  std::memcpy(&header, index_bytes, sizeof(header));
  if (std::memcmp(header.magic, kFastqIndexMagic, sizeof(header.magic)) != 0) {
    error = "not an index file";
    return false;
  }
  if (header.version != kFastqIndexVersion ||
      header.header_size != sizeof(FastqIndexHeader)) {
    error = "index file was written by another version";
    return false;
  }
  if (header.data_size != data_size ||
      header.data_fingerprint != fastq_data_fingerprint(data, data_size)) {
    error = "index file describes different data";
    return false;
  }
  // An index over inflated bytes carries the archive it came from; without it
  // the index could outlive the archive it describes and answer from a
  // decompression of data that no longer exists.
  if (header.source_size != 0 &&
      (origin == nullptr || header.source_size != origin_size ||
       header.source_fingerprint != fastq_data_fingerprint(origin,
                                                            origin_size))) {
    error = "index file describes a different compressed source";
    return false;
  }
  const std::uint64_t record_count = header.record_count;
  const std::uint64_t slot_count = header.slots_count;
  if (slot_count == 0 || (slot_count & (slot_count - 1)) != 0) {
    error = "index file has a corrupt lookup table size";
    return false;
  }
  // Bounds checked before any multiplication, because both counts come from
  // the file and a crafted one would otherwise overflow its way past the range
  // check below.  Entries are addressed by a 32-bit slot, so the count cannot
  // exceed that either.
  if (record_count > std::numeric_limits<std::uint32_t>::max() ||
      record_count > index_size / sizeof(FastqEntryDisk) ||
      slot_count > index_size / sizeof(std::uint32_t)) {
    error = "index file is truncated";
    return false;
  }
  // Every range is checked against the mapping before it is used.  An index
  // file is a cache, and a cache that was truncated by a full disk is a
  // realistic input; reading past the end of it is not an acceptable response.
  const auto fits = [index_size](std::uint64_t offset, std::uint64_t length) {
    return offset <= index_size && length <= index_size - offset;
  };
  if (!fits(header.entries_offset,
            record_count * sizeof(FastqEntryDisk)) ||
      !fits(header.names_offset, header.names_size) ||
      !fits(header.slots_offset, slot_count * sizeof(std::uint32_t))) {
    error = "index file is truncated";
    return false;
  }

  out.data_ = data;
  out.data_size_ = data_size;
  out.entries_ = reinterpret_cast<const FastqEntryDisk*>(index_bytes +
                                                          header.entries_offset);
  out.count_ = static_cast<std::size_t>(record_count);
  out.slots_ = reinterpret_cast<const std::uint32_t*>(index_bytes +
                                                      header.slots_offset);
  out.slot_mask_ = static_cast<std::size_t>(slot_count - 1);
  out.names_ = index_bytes + header.names_offset;
  out.names_size_ = static_cast<std::size_t>(header.names_size);
  return true;
}

std::string_view FastqIndexView::key_of_disk(
    const FastqEntryDisk& disk) const noexcept {
  if (disk.plain != 0) {
    return std::string_view(data_ + disk.offset + 1 + disk.name_offset,
                            disk.name_length);
  }
  return std::string_view(names_ + disk.name_offset, disk.name_length);
}

void FastqIndexView::restore(const FastqEntryDisk& disk,
                             FastqEntry& out) const noexcept {
  out.offset = static_cast<std::size_t>(disk.offset);
  out.end = out.offset + disk.end_delta;
  out.plain = disk.plain != 0;
  out.length = disk.length;
  out.title_length = disk.title_length;
  out.name_length = disk.name_length;
  if (out.plain) {
    out.sequence_offset = out.offset + disk.sequence_delta;
    out.quality_offset = out.offset + disk.quality_delta;
    out.name = data_ + out.offset + 1 + disk.name_offset;
  } else {
    out.sequence_offset = 0;
    out.quality_offset = 0;
    out.name = names_ + disk.name_offset;
  }
}

const FastqEntry* FastqIndexView::find(std::string_view name) const noexcept {
  if (count_ == 0) return nullptr;
  std::uint64_t h = fastq_key_hash(name.data(), name.size()) & slot_mask_;
  for (;;) {
    const std::uint32_t slot = slots_[h];
    if (slot == 0xFFFFFFFFu) return nullptr;
    const FastqEntryDisk& disk = entries_[slot];
    const std::string_view key = key_of_disk(disk);
    if (key.size() == name.size() &&
        std::memcmp(key.data(), name.data(), name.size()) == 0) {
      restore(disk, scratch_);
      return &scratch_;
    }
    h = (h + 1) & slot_mask_;
  }
}

void FastqIndexView::entry_at(std::size_t i, FastqEntry& out) const noexcept {
  restore(entries_[i], out);
}

std::string FastqIndexView::title(const FastqEntry& entry) const {
  return fastq_entry_title(data_, entry);
}

std::string FastqIndexView::sequence(const FastqEntry& entry) const {
  return fastq_entry_sequence(data_, entry);
}

std::string FastqIndexView::quality(const FastqEntry& entry) const {
  return fastq_entry_quality(data_, entry);
}

std::string FastqIndexView::sequence_slice(const FastqEntry& entry,
                                           std::size_t start,
                                           std::size_t end) const {
  return fastq_entry_slice(data_, entry, start, end);
}

}  // namespace biofasting
