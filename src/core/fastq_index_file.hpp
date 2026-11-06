// A FASTQ index that lives in a file instead of in the heap.
//
// The in-memory `FastqIndex` costs one pass over the file and about a hundred
// bytes of heap per record, every time it is opened.  On the million-read
// corpus that is half a second and a hundred megabytes before the first fetch;
// a tool with an index on disk (pyfastx writes a SQLite database next to the
// data) opens in under a millisecond, because opening is a `mmap` and nothing
// else.  That difference is not a constant factor -- it makes an "open, fetch
// one record, exit" program a different class of program -- so it is worth a
// file format.
//
// What is stored is deliberately not the build's in-memory structures.  A heap
// entry is 56 bytes of pointers that mean nothing in another process; the
// on-disk entry is 40 bytes of offsets and lengths that do.  The lookup table
// is a flat open-addressed array of indices, so a mapped lookup is two or
// three probes over pages the kernel can hand over without a copy, and the
// index can be larger than RAM without being slower than one that fits.
//
// What is *not* stored is any name string that the data already contains.
// For a plain record the key is a view into the FASTQ bytes at a known offset,
// so the index stores an offset into the header and not a copy of the name; a
// million copies of a million names is the thing that made the old build
// allocation-bound.  Only records the scanner had to assemble -- wrapped
// sequences, blank lines, a repeated '+' header -- have no such address, and
// those keys go into a side blob.
//
// Identity: the index is valid for exactly one file's bytes.  It carries the
// data's size and a fingerprint of its first and last pages, and both are
// checked on open.  A mismatch is not an error to report but a reason to
// rebuild: an index next to a file that has been rewritten is worse than no
// index at all, because it answers confidently and wrongly.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include "fastq.hpp"

namespace biofasting {

// Little-endian on disk, and refused rather than byte-swapped on a big-endian
// host: nothing here targets one, and an unexercised swap path is a second
// implementation with no test behind it.
inline constexpr char kFastqIndexMagic[8] = {'B', 'F', 'Q', 'I',
                                             'D', 'X', '0', '1'};
inline constexpr std::uint32_t kFastqIndexVersion = 1;

struct FastqIndexHeader {
  char magic[8];
  std::uint32_t version;
  std::uint32_t header_size;
  std::uint64_t data_size;
  std::uint64_t data_fingerprint;
  std::uint64_t record_count;
  std::uint64_t entries_offset;
  std::uint64_t names_offset;
  std::uint64_t names_size;
  std::uint64_t slots_offset;
  std::uint64_t slots_count;
  // The file the indexed bytes came from, when it is not the bytes themselves.
  // A gzipped source is indexed after inflation, so an index over the inflated
  // bytes would still "match" after the .gz was replaced -- the inflated size
  // and fingerprint say nothing about the compressed file.  These two fields
  // are what tie such an index to the archive it was built from; both are zero
  // when the indexed bytes are the source.
  std::uint64_t source_size;
  std::uint64_t source_fingerprint;
  std::uint8_t padding[32];  // header occupies 128 bytes, one cache line pair
};
static_assert(sizeof(FastqIndexHeader) == 128, "header must stay 128 bytes");

// One record, as a place in the data rather than as bytes.
//
// `end_delta`, `sequence_delta` and `quality_delta` are relative to `offset`
// because a record is small and an absolute 64-bit offset three times over is
// twelve bytes of high zeroes.  For a record that is not `plain`,
// `name_offset` indexes the side blob; for one that is, it is the key's
// position inside the header (after the '@'), so that even the key needs no
// copy.
struct FastqEntryDisk {
  std::uint64_t offset;
  std::uint32_t end_delta;
  std::uint32_t sequence_delta;
  std::uint32_t quality_delta;
  std::uint32_t name_offset;
  std::uint16_t title_length;
  std::uint16_t name_length;
  std::uint32_t length;
  std::uint8_t plain;
  std::uint8_t reserved[3];
};
static_assert(sizeof(FastqEntryDisk) == 40, "entry must stay 40 bytes");

// FNV-1a over the size and the first and last pages of the data.  A checksum of
// the whole file would be a third pass over it, which is the cost this file
// exists to avoid; a page at each end catches truncation, appending and
// replacement, which are the ways an index goes stale in practice.
std::uint64_t fastq_data_fingerprint(const char* data,
                                     std::size_t size) noexcept;

// Writes `index`, which must have been built over [data, data + size), to
// `path`.  `origin` is the file those bytes came from when they are not the
// bytes themselves -- a gzip archive, whose inflated contents are `data` -- and
// is checked on open so that replacing the archive invalidates the index; pass
// null for a plain file, whose identity is the data.  The write is to
// `path + ".tmp"` and a rename, so a reader either sees the old index or the
// new one and never half of one.  Returns false with `error` set on any
// failure -- including a record too long to encode -- and leaves no partial
// file behind.
bool write_fastq_index_file(const std::string& path, const FastqIndex& index,
                            const char* data, std::size_t size,
                            const char* origin, std::size_t origin_size,
                            std::string& error);

// A read-only index over a mapped index file, and a mapped data file.
//
// Nothing is allocated per record and nothing is copied out of the mapping:
// `find` fills one scratch entry with offsets and returns a pointer to it, and
// the field accessors copy the bytes the caller asked for.  The view borrows
// both buffers and outlives neither; the owner of the mappings keeps them
// alive for as long as this is used.
class FastqIndexView {
 public:
  FastqIndexView() = default;

  // Validates the mapped index against `data` and points the view at it.
  // `origin` is the compressed file the data was inflated from, when there is
  // one; an index that records a source refuses to open without it.  `error`
  // says which check failed, for the caller that wants to log why it rebuilt;
  // the caller is expected to treat every failure the same way.
  static bool open(const char* index_bytes, std::size_t index_size,
                   const char* data, std::size_t data_size,
                   const char* origin, std::size_t origin_size,
                   FastqIndexView& out, std::string& error);

  std::size_t size() const noexcept { return count_; }

  // The entry for `name`, or nullptr.  The pointer is valid until the next
  // call on this view.
  const FastqEntry* find(std::string_view name) const noexcept;

  // The i-th entry in file order, which is the order `keys()` promises.
  void entry_at(std::size_t i, FastqEntry& out) const noexcept;

  std::string_view key_of(const FastqEntry& entry) const noexcept {
    return std::string_view(entry.name, entry.name_length);
  }

  // The same four accessors as FastqIndex, over the same entry layout, so a
  // fetch cannot mean two different things depending on how the index was
  // opened.
  std::string title(const FastqEntry& entry) const;
  std::string sequence(const FastqEntry& entry) const;
  std::string quality(const FastqEntry& entry) const;
  std::string sequence_slice(const FastqEntry& entry, std::size_t start,
                             std::size_t end) const;

 private:
  void restore(const FastqEntryDisk& disk, FastqEntry& out) const noexcept;
  std::string_view key_of_disk(const FastqEntryDisk& disk) const noexcept;

  const char* data_ = nullptr;
  std::size_t data_size_ = 0;
  const FastqEntryDisk* entries_ = nullptr;
  std::size_t count_ = 0;
  const std::uint32_t* slots_ = nullptr;
  std::size_t slot_mask_ = 0;
  const char* names_ = nullptr;
  std::size_t names_size_ = 0;
  // Returned by `find`, overwritten by the next one: the alternative is an
  // allocation per lookup, which is what this index exists to avoid.
  mutable FastqEntry scratch_{};
};

}  // namespace biofasting
