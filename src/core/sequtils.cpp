#include "sequtils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace biofasting {
namespace {

// ---------------------------------------------------------------------------
// The byte tables.  Each is the *definition* of what its kernel counts, so
// every one of them is also asserted against the reference in the tests --
// a table that drifts is a kernel that answers a different question.
// ---------------------------------------------------------------------------

// Slot codes for the GC123 scan: 0 means "no frame counts this byte", and
// 1..4 are A, T, G and C.  Both cases, and nothing else -- which is how the
// reference's `codon[pos] == nt or codon[pos] == nt.lower()` reads, and why an
// N or a gap is counted by nobody.
constexpr std::array<std::uint8_t, 256> make_base_code() noexcept {
  std::array<std::uint8_t, 256> table{};
  constexpr char letters[] = "ATGCatgc";
  for (std::size_t i = 0; i + 1 < sizeof(letters); ++i) {
    table[static_cast<unsigned char>(letters[i])] =
        static_cast<std::uint8_t>(i % 4 + 1);
  }
  return table;
}

constexpr std::array<std::uint8_t, 256> kBaseCode = make_base_code();

// The same idea for GC_skew, which asks a smaller question: 1 is G, 2 is C.
constexpr std::array<std::uint8_t, 256> make_gc_code() noexcept {
  std::array<std::uint8_t, 256> table{};
  for (char c : {'G', 'g'}) table[static_cast<unsigned char>(c)] = 1;
  for (char c : {'C', 'c'}) table[static_cast<unsigned char>(c)] = 2;
  return table;
}

constexpr std::array<std::uint8_t, 256> kGcCode = make_gc_code();

// gc_fraction's three groups in one table, because one pass can serve all
// three modes and the reference's own definition is a sum of `count` calls
// over exactly these three sets: 1 is the GC set "CGScgs", 2 is the rest of
// the "remove" denominator "ATWUatwu", and 3..12 are the ten ambiguous letters
// "BDHKMNRVXY" that only the "weighted" mode looks at.  Note that S is in the
// GC set and W is in the AT set: both are ambiguous, and both are counted as
// if they were the one thing they can only half be, which is the reference's
// documented behaviour and not an oversight here.
constexpr std::array<std::uint8_t, 256> make_gc_slot() noexcept {
  std::array<std::uint8_t, 256> table{};
  for (char c : {'C', 'G', 'S', 'c', 'g', 's'}) {
    table[static_cast<unsigned char>(c)] = 1;
  }
  for (char c : {'A', 'T', 'W', 'U', 'a', 't', 'w', 'u'}) {
    table[static_cast<unsigned char>(c)] = 2;
  }
  constexpr char ambiguous[] = "BDHKMNRVXY";
  for (std::size_t i = 0; i + 1 < sizeof(ambiguous); ++i) {
    table[static_cast<unsigned char>(ambiguous[i])] =
        static_cast<std::uint8_t>(i + 3);
    table[static_cast<unsigned char>(ambiguous[i] | 0x20)] =
        static_cast<std::uint8_t>(i + 3);
  }
  return table;
}

constexpr std::array<std::uint8_t, 256> kGcSlot = make_gc_slot();

// ---------------------------------------------------------------------------
// The crc64 high table.  This one is *derived* rather than dumped, because the
// reference derives it too: `CheckSum._init_table_h` is fourteen lines of
// shift-and-xor with no data in it, and re-running that recurrence at compile
// time is the same statement of the same algorithm.  The test compares all 256
// entries against the reference's list, so a transcription slip cannot hide.
constexpr std::array<std::uint32_t, 256> make_table_h() noexcept {
  std::array<std::uint32_t, 256> table{};
  for (std::uint32_t i = 0; i < 256; ++i) {
    std::uint32_t part_l = i;
    std::uint32_t part_h = 0;
    for (int j = 0; j < 8; ++j) {
      const std::uint32_t rflag = part_l & 1U;
      part_l >>= 1;
      if (part_h & 1U) part_l |= 1U << 31;
      part_h >>= 1;
      if (rflag) part_h ^= 0xD8000000U;
    }
    table[i] = part_h;
  }
  return table;
}

constexpr std::array<std::uint32_t, 256> kTableH = make_table_h();

// ---------------------------------------------------------------------------
// The compensated sum.
//
// CPython 3.12 and later sum floats with Neumaier compensation, and
// `molecular_weight` is written as `sum(weight_table[x] for x in seq)`.  This
// is that algorithm, step for step: keep the running total and the running
// correction, and return their sum.  Replacing it with `total += x` is not a
// simplification, it is a different function -- see the header.
inline void neumaier_add(double& total, double& correction, double x) noexcept {
  const double next = total + x;
  if (std::fabs(total) >= std::fabs(x)) {
    correction += (total - next) + x;
  } else {
    correction += (x - next) + total;
  }
  total = next;
}

}  // namespace

// ---------------------------------------------------------------------------

Gc123Counts gc123(const char* data, std::size_t size) {
  Gc123Counts out;
  std::uint64_t* frames[3] = {out.counts[0], out.counts[1], out.counts[2]};

  // The frame of a byte is its position modulo three, and the reference works
  // in codons, so the two agree exactly: a partial trailing codon contributes
  // its first one or two bases to frames 0 and 1 and nothing to frame 2, which
  // is what the padding with spaces does on the reference's side.
  //
  // Unrolled by three so that each of the three counters is a fixed array the
  // compiler can keep in registers instead of computing an index per byte.
  std::size_t i = 0;
  for (; i + 3 <= size; i += 3) {
    for (int f = 0; f < 3; ++f) {
      const std::uint8_t code = kBaseCode[static_cast<unsigned char>(data[i + f])];
      if (code != 0) ++frames[f][code - 1];
    }
  }
  for (int f = 0; i < size; ++i, ++f) {
    const std::uint8_t code = kBaseCode[static_cast<unsigned char>(data[i])];
    if (code != 0) ++frames[f][code - 1];
  }
  return out;
}

std::vector<double> gc_skew(const char* data, std::size_t size,
                            std::size_t window) {
  if (window == 0) {
    // The reference reaches this through `range(0, len(seq), window)`, and the
    // caller has already reproduced that call for its message.  Kept here so
    // that the kernel is not the one place where a zero can be divided by.
    throw std::invalid_argument("range() arg 3 must not be zero");
  }
  std::vector<double> out;
  out.reserve(size / window + 1);
  for (std::size_t start = 0; start < size; start += window) {
    std::uint64_t g = 0;
    std::uint64_t c = 0;
    const std::size_t end = std::min(start + window, size);
    for (std::size_t i = start; i < end; ++i) {
      const std::uint8_t code = kGcCode[static_cast<unsigned char>(data[i])];
      g += static_cast<std::uint64_t>(code == 1);
      c += static_cast<std::uint64_t>(code == 2);
    }
    //   The subtraction happens in `double` and not in the 64-bit counters.
    //   `g - c` on unsigned ints wraps to 2**64 - d, so a window with more C
    //   than G -- 28% of the windows in that sweep, and approaching half as the
    //   window grows -- came out as a huge positive number in the 1e17..1e19
    //   range instead of a negative fraction.  The differential sweep in
    //   `tests/test_sequtils.py` is wrong on 884 of its 1,200 draws under the
    //   wrapped formula; the hand-written cases missed it.
    const double dg = static_cast<double>(g);
    const double dc = static_cast<double>(c);
    out.push_back(dg + dc == 0.0 ? 0.0 : (dg - dc) / (dg + dc));
  }
  return out;
}

GcCounts gc_counts(const char* data, std::size_t size) {
  GcCounts out;
  for (std::size_t i = 0; i < size; ++i) {
    const std::uint8_t slot = kGcSlot[static_cast<unsigned char>(data[i])];
    if (slot == 1) {
      ++out.gc;
    } else if (slot == 2) {
      ++out.at;
    } else if (slot > 2) {
      ++out.ambiguous[slot - 3];
    }
  }
  return out;
}

bool molecular_weight_mass(const char* data, std::size_t size,
                           const double* table, const unsigned char* valid,
                           double* total, unsigned char* bad) {
  double sum = 0.0;
  double correction = 0.0;
  //   `sum()` starts with an integer accumulator of zero and stays in it until
  //   the first float arrives; that item is then added to the accumulator
  //   *plainly*, with no compensation to correct, and only the second one onward
  //   is a Neumaier step.  It is not a formality: for "IYAR" under the Parker
  //   scale -- -8, -1.9, 2.1, 4.2 -- the reference answers -3.6 and compensating
  //   that first pair answers -3.5999999999999996.  See the header for what this
  //   reproduces and where its bound is.
  bool leading = true;
  for (std::size_t i = 0; i < size; ++i) {
    const unsigned char byte = static_cast<unsigned char>(data[i]);
    const unsigned char flag = valid[byte];
    if (flag == kAbsent) {
      *bad = byte;
      *total = sum + correction;
      return false;
    }
    const double x = table[byte];
    if (leading || flag == kPlain) {
      //   An `int` item anywhere, or any item at all while the accumulator is
      //   still the reference's integer one.  The flags are documented where
      //   they are declared; what matters here is that this is not an
      //   optimisation but the reference's own arithmetic, and that dropping it
      //   changes the last bit of ten of the twenty-eight gravy scales.
      sum += x;
      if (flag != kPlain) {
        leading = false;
      }
      continue;
    }
    neumaier_add(sum, correction, x);
  }
  *total = sum + correction;
  return true;
}

bool gcg(const char* data, std::size_t size, std::uint32_t* checksum) {
  std::uint32_t index = 0;
  std::uint64_t total = 0;
  for (std::size_t i = 0; i < size; ++i) {
    const unsigned char byte = static_cast<unsigned char>(data[i]);
    if (byte >= 0x80) return false;  // `char.upper()` is not a byte operation
    // `toupper` is only defined for unsigned char values and EOF, and only
    // for the C locale's letters -- which is exactly ASCII here, since the
    // byte is below 0x80.  So this is `char.upper()` on the reference's side
    // and on a Turkish locale's C library it is still ASCII A-Z, because C
    // locales never change the ASCII range.
    const unsigned char upper =
        (byte >= 'a' && byte <= 'z') ? static_cast<unsigned char>(byte - 32)
                                     : byte;
    ++index;
    total += static_cast<std::uint64_t>(index) * upper;
    if (index == 57) index = 0;
  }
  *checksum = static_cast<std::uint32_t>(total % 10000);
  return true;
}

void crc64(const char* data, std::size_t size, std::uint32_t* high,
           std::uint32_t* low) {
  // No decline here, and that is not an oversight: the reference folds each
  // character in with `ord(c) & 0xFF`, so a byte is already the whole of what
  // it looks at.  This is the one kernel in the file whose arithmetic is
  // defined on the low byte rather than on the character.
  std::uint32_t crcl = 0;
  std::uint32_t crch = 0;
  for (std::size_t i = 0; i < size; ++i) {
    const std::uint32_t byte = static_cast<unsigned char>(data[i]);
    const std::uint32_t shr = (crch & 0xFFU) << 24;
    const std::uint32_t temp1h = crch >> 8;
    const std::uint32_t temp1l = (crcl >> 8) | shr;
    crch = temp1h ^ kTableH[(crcl ^ byte) & 0xFFU];
    crcl = temp1l;
  }
  *high = crch;
  *low = crcl;
}

std::vector<std::uint32_t> crc64_table_h() {
  return std::vector<std::uint32_t>(kTableH.begin(), kTableH.end());
}

}  // namespace biofasting
