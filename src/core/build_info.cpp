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
      // Reported as named decisions rather than one flag string: see the note
      // at the top of build_config.hpp.in for why a flag string would be empty
      // on exactly the builds anyone would ask about.
      {"extra_cxx_flags", BIOFASTING_BUILD_EXTRA_CXX_FLAGS},
      {"relax_min_size", BIOFASTING_BUILD_RELAX_MIN_SIZE},
      {"stack_protector", BIOFASTING_BUILD_PROTECT_STACK},
  };
}

}  // namespace biofasting
