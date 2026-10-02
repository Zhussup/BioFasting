# Build of the BioFasting C++ core.
#
# Normal use goes through the PEP 517 frontend, which is what supplies cmake and
# ninja:
#
#     pip install -e .
#
# A bare `cmake -S . -B build` also works and is the faster loop while chasing a
# compiler error, but it cannot install a wheel.
#
# Deliberate: no -march=native and no per-file ISA flags.  A wheel is built on
# one machine and run on another, so the baseline binary has to be portable; ISA
# specialisation belongs behind the runtime dispatch in src/core/cpu_features.cpp
# and arrives with the first kernel (step 1.2).  A build that quietly required
# AVX2 would pass every test on the developer's laptop and die on a user's
# machine, which is the one failure mode this project cannot afford.

cmake_minimum_required(VERSION 3.15...3.30)

project(biofasting LANGUAGES CXX)

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

find_package(Python 3.10 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

set(BIOFASTING_CORE_SOURCES
  src/core/module.cpp
  src/core/build_info.cpp
  src/core/cpu_features.cpp)

nanobind_add_module(_core
  NB_STATIC
  # nanobind's default is to add -Os on top of the Release -O3, to keep wheels
  # small.  For a library whose entire purpose is throughput that is exactly
  # backwards, so it is switched off here rather than discovered later as an
  # unexplained few-percent gap against the C baseline.
  NOMINSIZE
  # nanobind also drops -fstack-protector by default, worth 1-2% at runtime and
  # ~10% of binary size.  A parser that will be pointed at untrusted files is
  # the wrong place to save that, so it is kept.  Remove this argument to match
  # nanobind's default; the cost is small and measurable.
  PROTECT_STACK
  ${BIOFASTING_CORE_SOURCES})

# Warnings are ours only.  NB_STATIC compiles nanobind's own sources into this
# target, and holding third-party code to -Werror would turn a nanobind patch
# release into our build failure.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set_source_files_properties(${BIOFASTING_CORE_SOURCES}
    PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra")
  if(BIOFASTING_WERROR)
    set_source_files_properties(${BIOFASTING_CORE_SOURCES}
      PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Werror")
  endif()
endif()

target_include_directories(_core PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core"
  "${CMAKE_CURRENT_BINARY_DIR}/generated")

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

set(BIOFASTING_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
string(STRIP "${BIOFASTING_CXX_FLAGS}" BIOFASTING_CXX_FLAGS)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/core/build_config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/biofasting/build_config.hpp"
  @ONLY)

# nanobind and scikit-build-core agree on this layout: the extension lands next
# to the pure-Python package so that `from . import _core` resolves in the wheel.
install(TARGETS _core LIBRARY DESTINATION biofasting)
