// gzip decompression, via libdeflate.
//
// The FASTQ kernel scans a span of bytes, and a gzip file is not one until it
// has been decompressed.  So the .gz path is: map the compressed file, inflate
// it into one buffer, hand that buffer to the same scanner.  One code path for
// the parsing, one for the decompression, and neither knows about the other.
//
// Whole-file inflation rather than a streaming feed, deliberately.  Streaming
// would mean the scanner had to survive a record split across block
// boundaries, which is a real design (and a real cost) bought for peak memory
// alone: the compressed file is mapped, not read, and the decompressed buffer
// is the only large allocation.  The corpus is 369 MB inflated.  A file that
// does not fit in memory is a phase-2 problem and saying so is better than
// pretending the streaming version exists.
//
// libdeflate rather than zlib because zlib is the thing being beaten: the
// README's own claim is that on real-world FASTQ the decompression, not the
// tokenising, is the first bottleneck, so using the same inflate as the
// baseline would leave the largest part of the gap on the table.

#pragma once

#include <cstddef>
#include <memory>
#include <string>

namespace biofasting {

// A decompressed stream, owned.
//
// Plain char storage rather than std::vector<char>: growing a vector
// value-initialises the new bytes, and this buffer exists to be immediately
// overwritten by the decompressor, so zeroing 369 MB first is pure cost.
class GzipBuffer {
 public:
  GzipBuffer() = default;
  GzipBuffer(GzipBuffer&&) noexcept = default;
  GzipBuffer& operator=(GzipBuffer&&) noexcept = default;
  GzipBuffer(const GzipBuffer&) = delete;
  GzipBuffer& operator=(const GzipBuffer&) = delete;

  const char* data() const noexcept { return data_.get(); }
  std::size_t size() const noexcept { return size_; }

 private:
  std::unique_ptr<char[]> data_;
  std::size_t size_ = 0;

  friend bool gunzip(const char*, std::size_t, GzipBuffer&, std::string&);
};

// Two magic bytes, not a file extension.  An extension is a claim about the
// contents; these two bytes are the format.
bool looks_gzipped(const char* data, std::size_t size) noexcept;

// Decompresses every concatenated gzip member in [data, data + size) into
// `out`.  Returns false with `error` set when the bytes are not a valid gzip
// stream; `out` is then empty.
bool gunzip(const char* data, std::size_t size, GzipBuffer& out,
            std::string& error);

}  // namespace biofasting
