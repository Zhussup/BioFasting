// Compile-time facts about this extension, for a bug report or for reproducing
// a benchmark on someone else's machine.
//
// These are the half of the answer that only the build knows.  The other half
// -- interpreter version, OS, running machine -- is added at runtime by
// src/biofasting/__init__.py, because a binary cannot recover from the CPU it
// is executing on which compiler and flags produced it.

#pragma once

#include <map>
#include <string>

namespace biofasting {

// A flat string map, not a struct: it is meant to be printed, pasted into an
// issue and compared between two machines, so a schema would only get in the
// way of that.
std::map<std::string, std::string> build_info();

}  // namespace biofasting
