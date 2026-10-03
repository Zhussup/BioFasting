#include "inflate.hpp"

#include <cstdint>
#include <cstring>
#include <utility>

#include <libdeflate.h>

namespace biofasting {
namespace {

struct DecompressorDeleter {
  void operator()(libdeflate_decompressor* d) const noexcept {
    libdeflate_free_decompressor(d);
  }
};

using Decompressor = std::unique_ptr<libdeflate_decompressor, DecompressorDeleter>;

// The gzip trailer's ISIZE field: the uncompressed size of the *last* member,
// modulo 2^32.  For the single-member files gzip produces -- everything the
// benchmark corpus has -- that is the exact size of the whole thing, so the
// buffer is allocated once and the decompressor never has to grow it.  When it
// is wrong (multi-member, or a stream over 4 GB) the growth path below handles
// it, and does so correctly rather than approximately.
std::size_t trailing_size_field(const char* data, std::size_t size) noexcept {
  const auto* p = reinterpret_cast<const unsigned char*>(data + size - 4);
  return static_cast<std::size_t>(p[0]) | (static_cast<std::size_t>(p[1]) << 8) |
         (static_cast<std::size_t>(p[2]) << 16) |
         (static_cast<std::size_t>(p[3]) << 24);
}

std::string describe(libdeflate_result result) {
  switch (result) {
    case LIBDEFLATE_BAD_DATA:
      return "the data is not a valid gzip stream";
    case LIBDEFLATE_SHORT_OUTPUT:
      return "the stream ended before the declared length";
    case LIBDEFLATE_INSUFFICIENT_SPACE:
      return "the output buffer is too small";
    default:
      return "libdeflate reported result " + std::to_string(result);
  }
}

}  // namespace

bool looks_gzipped(const char* data, std::size_t size) noexcept {
  return size >= 2 && static_cast<unsigned char>(data[0]) == 0x1f &&
         static_cast<unsigned char>(data[1]) == 0x8b;
}

bool gunzip(const char* data, std::size_t size, GzipBuffer& out,
            std::string& error) {
  out = GzipBuffer();
  if (!looks_gzipped(data, size)) {
    error = "not a gzip stream (the first two bytes are not 0x1f 0x8b)";
    return false;
  }
  // A gzip member with no content is still 20 bytes: 10 of header, 8 of
  // trailer, and a 2-byte empty deflate block.
  if (size < 20) {
    error = "not a gzip stream (shorter than an empty gzip member)";
    return false;
  }

  Decompressor decompressor(libdeflate_alloc_decompressor());
  if (!decompressor) {
    error = "libdeflate could not allocate a decompressor";
    return false;
  }

  std::size_t capacity = trailing_size_field(data, size);
  // ISIZE is 0 both for an empty member and for a 4 GB multiple, so it cannot
  // be trusted as a starting point on its own.
  if (capacity == 0) capacity = size * 4;
  std::unique_ptr<char[]> buffer(new char[capacity]);

  std::size_t in_pos = 0;
  std::size_t out_pos = 0;
  while (in_pos < size) {
    std::size_t consumed = 0;
    std::size_t produced = 0;
    const libdeflate_result result = libdeflate_gzip_decompress_ex(
        decompressor.get(), data + in_pos, size - in_pos, buffer.get() + out_pos,
        capacity - out_pos, &consumed, &produced);
    if (result == LIBDEFLATE_INSUFFICIENT_SPACE) {
      const std::size_t grown = capacity ? capacity * 2 : 64;
      std::unique_ptr<char[]> bigger(new char[grown]);
      std::memcpy(bigger.get(), buffer.get(), out_pos);
      buffer = std::move(bigger);
      capacity = grown;
      continue;  // retry this member against the larger buffer
    }
    if (result != LIBDEFLATE_SUCCESS) {
      error = describe(result);
      return false;
    }
    if (consumed == 0) {  // cannot happen for a success, but a loop that never
      error = "the gzip stream did not advance";  // advances never ends
      return false;
    }
    in_pos += consumed;
    out_pos += produced;
  }

  out.data_ = std::move(buffer);
  out.size_ = out_pos;
  return true;
}

}  // namespace biofasting
