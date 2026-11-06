#include "mapping.hpp"

#include <cstdint>

#if defined(__linux__) || defined(__APPLE__)
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace biofasting {
namespace {

#if defined(__linux__) || defined(__APPLE__)

std::size_t page_size() noexcept {
  // Asked once: `_SC_PAGESIZE` is a constant per system, and this is on the
  // path of every record of a streaming read.
  static const std::size_t size = [] {
    const long value = ::sysconf(_SC_PAGESIZE);
    return value > 0 ? static_cast<std::size_t>(value) : std::size_t{4096};
  }();
  return size;
}

// `madvise` wants a page-aligned start, and neither a record's offset nor the
// end of a scan is one.
void advise(const char* base, std::size_t from, std::size_t to,
            int advice) noexcept {
  if (to <= from) return;
  const std::uintptr_t page = static_cast<std::uintptr_t>(page_size());
  const std::uintptr_t begin =
      (reinterpret_cast<std::uintptr_t>(base) + from) & ~(page - 1);
  const std::uintptr_t end = reinterpret_cast<std::uintptr_t>(base) + to;
  if (end <= begin) return;
  ::madvise(reinterpret_cast<void*>(begin), end - begin, advice);
}

#endif

}  // namespace

void mapping_drop_pages(const char* base, std::size_t from,
                        std::size_t to) noexcept {
#if defined(__linux__) || defined(__APPLE__)
  advise(base, from, to, MADV_DONTNEED);
#else
  (void)base;
  (void)from;
  (void)to;
#endif
}

}  // namespace biofasting
