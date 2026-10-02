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
      {"cxx_flags", BIOFASTING_BUILD_CXX_FLAGS},
      {"cmake", BIOFASTING_BUILD_CMAKE_VERSION},
      {"nanobind", BIOFASTING_BUILD_NANOBIND_VERSION},
  };
}

}  // namespace biofasting
