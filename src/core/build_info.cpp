#include "build_info.hpp"

#include "biofasting/build_config.hpp"

namespace biofasting {

std::map<std::string, std::string> build_info() {
  return {
      {"version", BIOFASTING_BUILD_VERSION},
      {"compiler", BIOFASTING_BUILD_COMPILER},
      {"cxx_standard", BIOFASTING_BUILD_CXX_STANDARD},
      {"build_type", BIOFASTING_BUILD_TYPE},
      {"arch", BIOFASTING_BUILD_ARCH},
      {"cmake", BIOFASTING_BUILD_CMAKE_VERSION},
      {"nanobind", BIOFASTING_BUILD_NANOBIND_VERSION},
      // Which inflate produced a .gz row, and whether it travels inside the
      // wheel or has to be on the host.  Without it a gz number is not
      // reproducible and a portability failure is invisible.
      {"libdeflate", BIOFASTING_BUILD_LIBDEFLATE},
      // Reported as named decisions rather than one flag string: see the note
      // at the top of build_config.hpp.in for why a flag string would be empty
      // on exactly the builds anyone would ask about.
      {"extra_cxx_flags", BIOFASTING_BUILD_EXTRA_CXX_FLAGS},
      {"relax_min_size", BIOFASTING_BUILD_RELAX_MIN_SIZE},
      {"stack_protector", BIOFASTING_BUILD_PROTECT_STACK},
      //   The one decision against the compilers' defaults: the reference's
      //   arithmetic is interpreted IEEE double operations, and a contracted
      //   multiply-add rounds one step less than the interpreter's own, which
      //   moves the floating kernels off the contract on arm64.  Without this
      //   field a 3.13 answer that differs by one ulp gets argued from memory.
      {"fp_contract", BIOFASTING_BUILD_FP_CONTRACT},
  };
}

}  // namespace biofasting
