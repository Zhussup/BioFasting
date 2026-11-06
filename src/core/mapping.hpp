// Telling the kernel what a mapping is for.
//
// Two read paths in this package map a file and then walk it once, front to
// back: an index build, and a streaming read.  The kernel's default guess for a
// mapping -- keep everything you have read -- leaves the whole file resident.
// On the benchmark corpus that is 369 MB for a 369 MB file, and it reads in
// `VmHWM` as "this tool needs the file's size in RAM" even though every page is
// clean and reclaimable.  On a 100 GB file the same number is an obstacle to
// anyone who has a memory limit, whatever the reclaim semantics say.
//
// `mapping_drop_pages` is the other half of that story: it gives back the part
// of a mapping that has already been consumed, so such a read holds one window
// of the file rather than all of it, and the peak stops growing with the file.
//
// What is deliberately *not* here is advice for random access.  `MADV_RANDOM`
// is the obvious candidate for the index's fetches, and it was measured: on
// Linux 6.12 it changed neither the time nor the resident set of a thousand
// random fetches (8.19 against 8.54 us a fetch, 153 MB either way), because
// fault-around, not readahead, is what maps the surrounding pages.  Advice that
// has been measured not to help does not get shipped with a comment claiming it
// does.
//
// These are hints, not promises.  Each can fail -- an address the platform does
// not accept, a buffer that is not a mapping at all, a system without
// `madvise` -- and the right response to that failure is to carry on, because
// wrong advice costs a little performance and ignoring the hint costs nothing
// except the performance it was meant to save.  So nothing here returns
// anything, and nothing here is required for correctness.

#pragma once

#include <cstddef>

namespace biofasting {

// How much of a mapping stays resident behind a scan that gives pages back.
// A window rather than a number derived from the file: what a scan needs is the
// record it is on, and 32 MB is far enough ahead that the readahead the kernel
// has already done is not thrown away immediately, while small enough that the
// resident set of a build or a stream stays a few tens of megabytes rather than
// the file's size.  Pages are released in whole windows, so the resident set of
// a scan runs between one and two of these -- and, importantly, does not grow
// with the file.
inline constexpr std::size_t kMappingWindow = std::size_t{32} << 20;

// Drops the process's own pages for [base + from, base + to), leaving the bytes
// in the page cache and in the file.
//
// **Only for file-backed mappings.**  On anonymous memory the same call zeroes
// the range and frees it, which is why every caller of this one is reading a
// mapped file and never a buffer somebody allocated.  A range that is not a
// whole number of pages is widened to whole pages, and `to <= from` does
// nothing.
void mapping_drop_pages(const char* base, std::size_t from,
                        std::size_t to) noexcept;

}  // namespace biofasting
