// Runtime ISA dispatch -- the interface, not a kernel.
//
// Principle 4 of the project: one dispatch interface, backends later.  Nothing
// dispatches yet; step 1.2 is the first caller.  What exists here is the ladder
// a kernel will be compiled against and the detection that decides which rung
// that kernel runs on, so that the choice is never made by -march=native and
// never by a guess.
//
// The names are strings rather than an enum on purpose: an enum would make the
// set of levels part of the ABI, and a wheel built by an older compiler must
// still be loadable by a newer Python that knows about more levels.

#pragma once

#include <string>
#include <vector>

namespace biofasting {

// The ISA ladder this build defines, ascending, "baseline" at index 0.
//
// "baseline" is not a feature flag but a promise: the instructions every target
// of this wheel is allowed to assume.  x86-64's psABI guarantees SSE2, and
// ARMv8-A makes NEON mandatory, so on both the compiler already emits vector
// code at -O2 without being asked.  The rungs above it are the ones that need
// detection.
std::vector<std::string> cpu_levels();

// The subset of cpu_levels() that this CPU actually supports, ascending.
// Always contains "baseline": if the extension loaded, its floor is available
// by definition.
std::vector<std::string> supported_cpu_levels();

// The top of supported_cpu_levels() -- the rung a kernel would dispatch to.
// Equals "baseline" wherever detection is not implemented (see the note on
// MSVC and on non-x86/ARM targets in cpu_features.cpp).
std::string best_cpu_level();

}  // namespace biofasting
