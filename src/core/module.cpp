// The nanobind module.  Deliberately thin: it declares the names and hands off
// to build_info.cpp and cpu_features.cpp, so that adding a kernel in step 1.2
// does not mean editing the place where the ABI is defined.

#include <nanobind/nanobind.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "biofasting/build_config.hpp"
#include "build_info.hpp"
#include "cpu_features.hpp"

namespace nb = nanobind;

NB_MODULE(_core, m) {
  m.doc() =
      "C++ core of BioFasting.\n\n"
      "Private: import `biofasting` instead.  Nothing in this module is a "
      "stable interface, and in particular the extension module name may "
      "change without a deprecation cycle while the project is pre-alpha.";

  m.attr("__version__") = BIOFASTING_BUILD_VERSION;

  m.def("build_info", &biofasting::build_info,
        "How this extension was compiled, as plain strings.  The runtime half "
        "of the same answer is added by biofasting.build_info().");

  m.def("cpu_levels", &biofasting::cpu_levels,
        "The ISA ladder this build defines, ascending, \"baseline\" first.");

  m.def("supported_cpu_levels", &biofasting::supported_cpu_levels,
        "The subset of cpu_levels() the running CPU supports.  Detected once, "
        "on first call.");

  m.def("best_cpu_level", &biofasting::best_cpu_level,
        "The highest supported level: the rung a kernel dispatches to.");
}
