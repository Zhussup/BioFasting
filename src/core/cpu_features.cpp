#include "cpu_features.hpp"

// Detection is written against the compiler and the target, never against the
// machine that happens to be running the build.  Two targets are supported --
// the ones the wheel matrix covers -- and everything else degrades to
// "baseline" rather than reporting a level it cannot verify.

#if defined(__x86_64__) || defined(_M_X64)
#  define BIOFASTING_TARGET_X86_64 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#  define BIOFASTING_TARGET_AARCH64 1
#endif

// __builtin_cpu_supports is GCC and Clang only.  On MSVC (and on any target not
// listed above) there is no detection, so the ladder is just the floor.  That
// under-reports performance and never over-reports capability, which is the
// direction this project wants to be wrong in.
#if defined(BIOFASTING_TARGET_X86_64) && \
    (defined(__GNUC__) || defined(__clang__))
#  define BIOFASTING_HAVE_X86_DETECTION 1
#endif

// SVE exists on Linux-aarch64 servers and not on Apple silicon at all, so the
// query is both platform- and header-dependent.
#if defined(BIOFASTING_TARGET_AARCH64) && defined(__linux__)
#  include <sys/auxv.h>
#  if __has_include(<asm/hwcap.h>)
#    include <asm/hwcap.h>
#  endif
#endif

#if defined(BIOFASTING_TARGET_AARCH64) && defined(__linux__) && \
    defined(HWCAP_SVE)
#  define BIOFASTING_HAVE_SVE_DETECTION 1
#endif

namespace biofasting {
namespace {

#if defined(BIOFASTING_HAVE_X86_DETECTION)

// This has to be a macro, not a function: __builtin_cpu_supports resolves its
// argument at compile time into a CPUID feature bit and rejects anything that
// is not a string literal ("parameter to builtin must be a string constant or
// literal" -- gcc 14).  A function taking const char* therefore does not
// compile, which is how this was discovered.
//
// __builtin_cpu_init() populates the cached feature table; GCC requires it
// before the CPUID builtins in a static-initialisation context, and afterwards
// it is a guarded no-op.
//
// The AVX-family answers account for XCR0, not only CPUID, so they mean "this
// process can execute it" rather than "the silicon has it".  Skipping the XCR0
// half is the classic way to build a wheel that SIGILLs under a hypervisor
// which does not save the upper vector state.
#define BIOFASTING_X86_HAS(feature) \
  (__builtin_cpu_init(), __builtin_cpu_supports(feature))

#endif

#if defined(BIOFASTING_HAVE_SVE_DETECTION)

bool has_sve() { return (getauxval(AT_HWCAP) & HWCAP_SVE) != 0; }

#endif

std::vector<std::string> detect() {
  std::vector<std::string> supported{"baseline"};

#if defined(BIOFASTING_HAVE_X86_DETECTION)
  // Ascending, and each rung probed on its own: a CPU can have AVX2 without
  // AVX-512, and AVX-512F without the later AVX-512 subsets.
  if (BIOFASTING_X86_HAS("sse4.2")) supported.emplace_back("sse4.2");
  if (BIOFASTING_X86_HAS("avx2")) supported.emplace_back("avx2");
  if (BIOFASTING_X86_HAS("avx512f")) supported.emplace_back("avx512f");
#elif defined(BIOFASTING_HAVE_SVE_DETECTION)
  if (has_sve()) supported.emplace_back("sve");
#endif

  return supported;
}

}  // namespace

std::vector<std::string> cpu_levels() {
#if defined(BIOFASTING_HAVE_X86_DETECTION)
  return {"baseline", "sse4.2", "avx2", "avx512f"};
#elif defined(BIOFASTING_TARGET_AARCH64)
  // NEON is absent from this ladder on purpose: ARMv8-A mandates it, so there
  // is no aarch64 build in which it could be the thing being detected.
  return {"baseline", "sve"};
#else
  return {"baseline"};
#endif
}

std::vector<std::string> supported_cpu_levels() {
  // CPUID is not free and the answer cannot change while the process runs, so
  // it is computed once.  Function-local static: initialisation is thread-safe
  // and happens on first use, which also keeps the cost off the import path
  // until something actually asks.
  static const std::vector<std::string> supported = detect();
  return supported;
}

std::string best_cpu_level() {
  // supported_cpu_levels() is built ascending and detect() never reorders it,
  // so the last element is the highest rung.  A test asserts this rather than
  // trusting the comment.
  return supported_cpu_levels().back();
}

}  // namespace biofasting
