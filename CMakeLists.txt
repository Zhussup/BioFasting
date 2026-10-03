# Build of the BioFasting C++ core.
#
# Normal use goes through the PEP 517 frontend, which is what supplies cmake and
# ninja:
#
#     pip install -e .
#
# A bare cmake build is the faster loop while chasing a compiler error, but it
# cannot install a wheel and it needs both tools on PATH plus nanobind's cmake
# directory named explicitly -- verified, not assumed:
#
#     export PATH="$PWD/.venv/bin:$PATH"
#     cmake -S . -B build/raw -G Ninja -DPython_EXECUTABLE="$PWD/.venv/bin/python" \
#           -Dnanobind_DIR="$(.venv/bin/python -c 'import nanobind; print(nanobind.cmake_dir())')"
#     cmake --build build/raw
#
# Deliberate: no -march=native and no per-file ISA flags.  A wheel is built on
# one machine and run on another, so the baseline binary has to be portable; ISA
# specialisation belongs behind the runtime dispatch in src/core/cpu_features.cpp
# and arrives with the first kernel (step 1.2).  A build that quietly required
# AVX2 would pass every test on the developer's laptop and die on a user's
# machine, which is the one failure mode this project cannot afford.

cmake_minimum_required(VERSION 3.15...3.30)

# C as well as C++: libdeflate is C99 and is compiled by this build (see the
# libdeflate block below).  The language is enabled unconditionally rather than
# only when vendoring is on, so there is one configure path and not two.
project(biofasting LANGUAGES C CXX)

# scikit-build-core injects the PEP 621 version from pyproject.toml.  A plain
# cmake build has no wheel metadata to read from, so it says so explicitly
# rather than reporting a version it would be guessing at.
if(DEFINED SKBUILD_PROJECT_VERSION)
  set(BIOFASTING_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
  set(BIOFASTING_VERSION "0.0.0.dev0")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# CMake 4.x scans every translation unit for C++20 module dependencies by
# default under Ninja.  Nothing here uses modules, and the scan adds a
# `-fmodules-ts` pass to every compile, so it is off.
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Two of nanobind's defaults are wrong for this library, and both are named
# rather than hard-coded so that build_info() can report which way they went and
# so that flipping one is a configure flag instead of an edit.
#
# RELAX_MIN_SIZE: nanobind adds -Os on top of Release's -O3 to keep wheels small.
# For a library whose whole purpose is throughput that is backwards.
option(BIOFASTING_RELAX_MIN_SIZE "Compile for speed, not wheel size (-O3 not -Os)" ON)
#
# PROTECT_STACK: nanobind drops -fstack-protector, worth 1-2% at runtime and
# ~10% of binary size.  A parser that will be pointed at untrusted files is the
# wrong place to save that, so it is kept.  -DBIOFASTING_PROTECT_STACK=OFF
# matches nanobind's default.
option(BIOFASTING_PROTECT_STACK "Keep the stack protector on the C++ core" ON)
#
# WERROR: warnings as errors.  Off by default locally, so that a compiler
# version bump does not become a hard blocker mid-debug; CI turns it on.
option(BIOFASTING_WERROR "Treat warnings in our own sources as errors" OFF)

find_package(Python 3.10 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

set(BIOFASTING_CORE_SOURCES
  src/core/module.cpp
  src/core/build_info.cpp
  src/core/cpu_features.cpp
  src/core/line.cpp
  src/core/fastq.cpp
  src/core/fasta.cpp
  src/core/flatfile.cpp
  src/core/annotations.cpp
  src/core/feature.cpp
  src/core/location.cpp
  src/core/seqops.cpp
  src/core/sequtils.cpp
  src/core/protparam.cpp
  src/core/write.cpp
  src/core/inflate.cpp)

# libdeflate: the gz path's decompressor, and the one third-party dependency the
# core takes.  Chosen over zlib because zlib is part of the baseline being
# beaten, so matching its inflate would leave the larger half of the gap where
# it is.
#
# It is **vendored**: release 1.23 lives in third_party/libdeflate, put there and
# kept honest by tools/vendor_libdeflate.py, and this build compiles it into a
# static PIC archive that travels inside the wheel.  Taking it from the host
# instead made the wheel's portability a property of the machine that built it:
# a wheel linked against a shared libdeflate imports only where libdeflate is
# installed, which is the "works on my laptop" failure this project exists not
# to ship.  Linking the host's *static* archive is not the way out either --
# Debian's libdeflate.a is built without -fPIC and fails outright when linked
# into a shared module, with "relocation R_X86_64_PC32 ... can not be used when
# making a shared object".  So the choice is between vendoring the source and
# shipping a wheel that only runs where it was built.
#
# Two options, the second of which exists for a distro packager who wants to
# build against their own libdeflate rather than for anyone benchmarking:
#
#   -DBIOFASTING_VENDORED_LIBDEFLATE=OFF   use the host's libdeflate
#   -DBIOFASTING_LIBDEFLATE_STATIC=ON     ...and link its static archive
#
# Whichever way it goes is reported by build_info(), because "which inflate
# produced this gz number, and does the wheel carry it" is not a detail.
option(BIOFASTING_VENDORED_LIBDEFLATE
  "Compile the pinned libdeflate in third_party/ rather than using the host's" ON)
option(BIOFASTING_LIBDEFLATE_STATIC
  "With BIOFASTING_VENDORED_LIBDEFLATE=OFF, link the host's static archive" OFF)

if(BIOFASTING_VENDORED_LIBDEFLATE)

  set(BIOFASTING_LIBDEFLATE_INCLUDE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/libdeflate")
  if(NOT EXISTS "${BIOFASTING_LIBDEFLATE_INCLUDE_DIR}/libdeflate_sources.cmake")
    message(FATAL_ERROR
      "third_party/libdeflate is missing or incomplete.  Restore it with\n"
      "    python3 tools/vendor_libdeflate.py\n"
      "or configure with -DBIOFASTING_VENDORED_LIBDEFLATE=OFF to use the "
      "host's libdeflate instead.")
  endif()

  # libdeflate selects its SIMD paths with intrinsics guarded by compiler
  # version checks, which is not enough when the compiler is much newer than the
  # assembler it drives: the assembler can then reject an instruction the
  # compiler was happy to emit.  Upstream probes for exactly that and falls back
  # with -DLIBDEFLATE_ASSEMBLER_DOES_NOT_SUPPORT_*.  The probe is carried over
  # unchanged, because the point of vendoring a release is to build what its
  # author builds -- without it a machine where upstream compiles would fail
  # here, which is the portability regression vendoring was meant to remove.
  # What it prevents is a hard assembler error, never wrong output.
  function(biofasting_check_assembler_support feature assembly_code)
    execute_process(COMMAND echo "${assembly_code}"
      COMMAND "${CMAKE_C_COMPILER}" -c -x assembler -o /dev/null -
      RESULT_VARIABLE result ERROR_QUIET)
    if(NOT result STREQUAL "0")
      add_definitions(-DLIBDEFLATE_ASSEMBLER_DOES_NOT_SUPPORT_${feature})
      message(STATUS
        "libdeflate: the assembler rejects ${feature}, which this compiler "
        "supports; building without it.  Upgrading binutils is recommended.")
    endif()
  endfunction()
  if(UNIX AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND "${CMAKE_C_COMPILER}" -dumpmachine
      OUTPUT_VARIABLE biofasting_machine OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(biofasting_machine MATCHES "^(x86_64|i[3-6]86)")
      if(CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 8.1)
        biofasting_check_assembler_support(AVX512VNNI
          "vpdpbusd %zmm0, %zmm0, %zmm0")
        biofasting_check_assembler_support(VPCLMULQDQ
          "vpclmulqdq $0, %zmm0, %zmm0, %zmm0")
      endif()
      if(CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 12.1)
        biofasting_check_assembler_support(AVX_VNNI
          "{vex} vpdpbusd %ymm0, %ymm0, %ymm0")
      endif()
    elseif(biofasting_machine MATCHES "^aarch64")
      if(CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 8.1)
        biofasting_check_assembler_support(DOTPROD
          ".arch armv8.2-a+dotprod\nudot v0.4s, v0.16b, v0.16b")
      endif()
      if(CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL 9.1)
        biofasting_check_assembler_support(SHA3
          ".arch armv8.2-a+sha3\neor3 v0.16b, v0.16b, v0.16b, v0.16b")
      endif()
    endif()
  endif()

  # The source list is upstream's own, written next to the sources by
  # tools/vendor_libdeflate.py out of upstream's CMakeLists rather than copied
  # by hand.  A hand-copied list is what goes stale, and its failure is silent:
  # a new source file sits in the tree, never gets compiled, and the library
  # builds and behaves like the release before it.
  include("${BIOFASTING_LIBDEFLATE_INCLUDE_DIR}/libdeflate_sources.cmake")
  # The paths in that list are upstream's, so they are relative to upstream's
  # root.  Joining them here keeps the generated file a readable copy of
  # upstream's own list rather than a CMake expression that only means anything
  # once a variable is set.
  list(TRANSFORM BIOFASTING_LIBDEFLATE_SOURCES
    PREPEND "${BIOFASTING_LIBDEFLATE_INCLUDE_DIR}/")

  add_library(biofasting_libdeflate STATIC ${BIOFASTING_LIBDEFLATE_SOURCES})
  # PIC because this archive goes into a shared object.  That is exactly the
  # property the host's archive lacks, and the reason this target exists.
  set_target_properties(biofasting_libdeflate PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_STANDARD 99)
  # Upstream's own Release settings (-O2 -DNDEBUG), not this project's -O3.  The
  # point of vendoring a release is that it is the same code compiled the same
  # way, and a quiet -O3 here would put every gz row of bench/targets.md out of
  # reach of a re-measurement.  Stated rather than inherited, for the same
  # reason the stack protector is: a distro could set these differently and then
  # the number would belong to a build nobody described.
  target_compile_definitions(biofasting_libdeflate PRIVATE
    "$<$<CONFIG:Release>:NDEBUG>")
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(biofasting_libdeflate PRIVATE
      "$<$<CONFIG:Release>:-O2>")
  endif()
  target_include_directories(biofasting_libdeflate PRIVATE
    "${BIOFASTING_LIBDEFLATE_INCLUDE_DIR}")

  set(BIOFASTING_LIBDEFLATE_LIBRARY biofasting_libdeflate)
  set(BIOFASTING_LIBDEFLATE_LINKAGE "vendored, static")

else()

  find_path(BIOFASTING_LIBDEFLATE_INCLUDE_DIR NAMES libdeflate.h)
  if(BIOFASTING_LIBDEFLATE_STATIC)
    find_library(BIOFASTING_LIBDEFLATE_STATIC_LIBRARY NAMES libdeflate.a)
    if(NOT BIOFASTING_LIBDEFLATE_STATIC_LIBRARY)
      message(FATAL_ERROR
        "BIOFASTING_LIBDEFLATE_STATIC is ON but no libdeflate.a was found.")
    endif()
    set(BIOFASTING_LIBDEFLATE_LIBRARY "${BIOFASTING_LIBDEFLATE_STATIC_LIBRARY}")
    set(BIOFASTING_LIBDEFLATE_LINKAGE "host, static")
  else()
    find_library(BIOFASTING_LIBDEFLATE_SHARED_LIBRARY NAMES deflate)
    if(NOT BIOFASTING_LIBDEFLATE_SHARED_LIBRARY)
      message(FATAL_ERROR
        "libdeflate was not found.  It is required for the gzip path.  Debian "
        "and Ubuntu: apt install libdeflate-dev.  macOS: brew install "
        "libdeflate.  Or drop -DBIOFASTING_VENDORED_LIBDEFLATE=OFF and let the "
        "vendored copy in third_party/ be used instead.")
    endif()
    set(BIOFASTING_LIBDEFLATE_LIBRARY "${BIOFASTING_LIBDEFLATE_SHARED_LIBRARY}")
    set(BIOFASTING_LIBDEFLATE_LINKAGE "host, shared")
  endif()

  if(NOT BIOFASTING_LIBDEFLATE_INCLUDE_DIR)
    message(FATAL_ERROR
      "libdeflate's headers were not found next to its library.  Install the "
      "development package (libdeflate-dev / brew install libdeflate).")
  endif()

endif()

# A gz benchmark number is only reproducible if it says which inflate produced
# it, so the version travels into build_info() like the compiler does.
file(STRINGS "${BIOFASTING_LIBDEFLATE_INCLUDE_DIR}/libdeflate.h"
  BIOFASTING_LIBDEFLATE_VERSION_LINE REGEX "^#define LIBDEFLATE_VERSION_STRING")
string(REGEX REPLACE ".*\"([^\"]+)\".*" "\\1" BIOFASTING_LIBDEFLATE_VERSION
  "${BIOFASTING_LIBDEFLATE_VERSION_LINE}")
if(NOT BIOFASTING_LIBDEFLATE_VERSION)
  set(BIOFASTING_LIBDEFLATE_VERSION "unknown")
endif()
set(BIOFASTING_LIBDEFLATE_REPORTED
  "${BIOFASTING_LIBDEFLATE_VERSION} (${BIOFASTING_LIBDEFLATE_LINKAGE})")

set(BIOFASTING_NANOBIND_ARGS NB_STATIC)
if(BIOFASTING_RELAX_MIN_SIZE)
  list(APPEND BIOFASTING_NANOBIND_ARGS NOMINSIZE)
endif()
if(BIOFASTING_PROTECT_STACK)
  list(APPEND BIOFASTING_NANOBIND_ARGS PROTECT_STACK)
endif()

nanobind_add_module(_core ${BIOFASTING_NANOBIND_ARGS} ${BIOFASTING_CORE_SOURCES})

# Warnings are ours only.  NB_STATIC compiles nanobind's own sources into this
# target, and holding third-party code to -Werror would turn a nanobind patch
# release into our build failure.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(BIOFASTING_CORE_FLAGS "-Wall;-Wextra;-ffp-contract=off")
  #   FP contract is off, stated in the build and not left to the compiler:
  #   the contract with the reference is bit-identity, and on arm64 -- Apple
  #   Silicon and the cibuildwheel aarch64 matrix alike -- both compilers
  #   contract `w * front + w * back` into a fused multiply-add by default,
  #   which rounds one step less and moves results by one ulp against
  #   Biopython's interpreted arithmetic.  Reproduced on x86-64 built with
  #   -mfma: the same test, the same index-0 value the macos-14 job printed.
  if(BIOFASTING_WERROR)
    list(APPEND BIOFASTING_CORE_FLAGS "-Werror")
  endif()
  # Stated explicitly rather than left to the distro.  GCC's own default is
  # -fno-stack-protector; Debian and Ubuntu patch it to `strong`, macOS does
  # not.  Relying on that would make build_info()'s "stack_protector: yes" true
  # on the machine that built the wheel and false on the machine that reports
  # from it, which is the wrong way round for a security claim.
  if(BIOFASTING_PROTECT_STACK)
    list(APPEND BIOFASTING_CORE_FLAGS "-fstack-protector-strong")
  endif()
  set_source_files_properties(${BIOFASTING_CORE_SOURCES}
    PROPERTIES COMPILE_OPTIONS "${BIOFASTING_CORE_FLAGS}")
endif()

target_include_directories(_core PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core"
  "${CMAKE_CURRENT_BINARY_DIR}/generated"
  "${BIOFASTING_LIBDEFLATE_INCLUDE_DIR}")

target_link_libraries(_core PRIVATE "${BIOFASTING_LIBDEFLATE_LIBRARY}")

# The target architecture of the build, not of the machine running cmake:
# cibuildwheel cross-compiles linux-aarch64, and the C++ preprocessor macros in
# cpu_features.cpp see the target as well, so the two must agree.  A test
# compares this against platform.machine() for exactly that reason.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
  set(BIOFASTING_ARCH "x86_64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
  set(BIOFASTING_ARCH "aarch64")
else()
  set(BIOFASTING_ARCH "${CMAKE_SYSTEM_PROCESSOR}")
endif()

# Facts a bug report needs.  Recorded at build time rather than guessed at
# runtime, because "which compiler and flags produced this binary" is not
# something the binary can recover from the CPU it is running on.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.__version__)"
  OUTPUT_VARIABLE BIOFASTING_NANOBIND_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT BIOFASTING_NANOBIND_VERSION)
  set(BIOFASTING_NANOBIND_VERSION "unknown")
endif()

# Deliberately NOT reported as a single flag string.  Under scikit-build-core
# and nanobind the interesting flags are per-target COMPILE_OPTIONS and
# generator expressions, which resolve to nothing useful at configure time;
# CMAKE_CXX_FLAGS itself is empty on a stock wheel.  Reporting it would put a
# confidently empty "cxx_flags" field in every bug report.  What is reported
# instead is the caller's global flags (non-empty only when someone passed
# CMAKE_ARGS) plus the two build decisions above as named yes/no.
set(BIOFASTING_EXTRA_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
string(STRIP "${BIOFASTING_EXTRA_CXX_FLAGS}" BIOFASTING_EXTRA_CXX_FLAGS)
if(NOT BIOFASTING_EXTRA_CXX_FLAGS)
  set(BIOFASTING_EXTRA_CXX_FLAGS "none")
endif()

if(BIOFASTING_RELAX_MIN_SIZE)
  set(BIOFASTING_RELAX_MIN_SIZE_REPORTED "yes (-O3, not -Os)")
else()
  set(BIOFASTING_RELAX_MIN_SIZE_REPORTED "no (nanobind's -Os)")
endif()

if(BIOFASTING_PROTECT_STACK)
  set(BIOFASTING_PROTECT_STACK_REPORTED "yes (-fstack-protector-strong)")
else()
  set(BIOFASTING_PROTECT_STACK_REPORTED "no (nanobind default)")
endif()

#   The one decision made against the compilers' defaults rather than nanobind's,
#   so it has no option: every value that keeps the contract true is "off", and
#   turning contraction back on changes the answers on arm64.
set(BIOFASTING_FP_CONTRACT_REPORTED "off (-ffp-contract=off)")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core/build_config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/biofasting/build_config.hpp"
  @ONLY)

# nanobind and scikit-build-core agree on this layout: the extension lands next
# to the pure-Python package so that `from . import _core` resolves in the wheel.
install(TARGETS _core LIBRARY DESTINATION biofasting)
