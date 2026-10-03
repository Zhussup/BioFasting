// The measurement kernels: composition, mass and checksums.
//
// seqops.cpp holds the operations that *produce* a sequence out of another one
// -- reverse complement, translation -- and the two countings that were already
// there.  This file holds the other half of `Bio.SeqUtils`: the functions that
// look at a sequence and return a number.  They are separated because the
// contract is different in kind.  An operation has to decide what an unknown
// byte means and is judged on the bytes it emits; a measurement counts, and the
// arithmetic that turns those counts into an answer belongs where it can be
// read.
//
// So the split here is deliberate and consistent: **a kernel counts, Python
// decides**.  `gc123` returns twelve counters and not four percentages,
// `gc_skew` returns a pair per window and not the ratios, because the reference
// computes those with `/` on Python floats and has three separate ways of
// failing -- `GC123("")` raises ZeroDivisionError, a frame with no A/C/G/T in it
// answers the *integer* 0, `gc_fraction("")` answers 0.0 -- and none of those
// distinctions survives a round trip through a C++ double.
//
// The one kernel that does more than count is molecular_weight's sum, and it
// does more because it has to.  See the comment on Neumaier below.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace biofasting {

// ---------------------------------------------------------------------------
// GC123.
//
// Biopython counts, for each of the three reading frames, how many of the
// A, T, G and C letters sit in it; every other byte -- including the one or two
// bases a partial trailing codon is padded with -- is counted by nobody.  That
// is the whole kernel.  `counts[frame][base]` is in the order A, T, G, C, and
// covers both cases.
struct Gc123Counts {
  std::uint64_t counts[3][4] = {};
};

Gc123Counts gc123(const char* data, std::size_t size);

// ---------------------------------------------------------------------------
// GC_skew.
//
// One ratio per window of `window` bytes, in order, counting both cases and
// ignoring everything that is not a G or a C.  The reference slices
// `seq[i:i+window]` and calls `.count` four times, so the windows are
// non-overlapping and the last one is short; `window` must be positive and the
// caller has already checked it.
//
// Unlike the other kernels this one returns the answer and not the counts, and
// the exception is worth justifying: the reference's zero-window case is
// `except ZeroDivisionError: skew = 0.0`, and every other path is a true
// division of two ints, so the result is a float either way and no distinction
// is being thrown away.  What is being thrown away is 400,000 Python divisions
// on a chromosome, which is the part that costs.
std::vector<double> gc_skew(const char* data, std::size_t size,
                            std::size_t window);

// ---------------------------------------------------------------------------
// gc_fraction's other two modes.
//
// One pass answers all three, but `remove` -- the default -- is answered by the
// vectorised kernel in seqops.cpp and this is what the other two use.  Keeping
// the fast path where it was measured is the whole reason these are two
// functions and not one.
//
// `gc` is the "CGScgs" count, `at` the "ATWUatwu" one, and `ambiguous[i]` is
// `seq.count(x) + seq.count(x.lower())` for x in "BDHKMNRVXY" in that order --
// the reference's own expression, ready to be multiplied by `_gc_values` on the
// Python side, where the Neumaier sum that consumes them can be written the way
// the reference writes it.
struct GcCounts {
  std::uint64_t gc = 0;
  std::uint64_t at = 0;
  std::uint64_t ambiguous[10] = {};
};

GcCounts gc_counts(const char* data, std::size_t size);

// ---------------------------------------------------------------------------
// molecular_weight's and gravy's inner sum.
//
// `table` is 256 doubles indexed by byte, and `valid` is 256 flags saying what
// each of them is.  The reference's weight dicts are four, four, or twenty-two
// entries long and a caller's hydrophobicity scale is twenty, and handing over
// a sparse table in a dense array is what keeps the per-letter branch out of the
// loop.
//
// The flag is not a boolean, because CPython's `sum()` -- which is what both
// callers spell this line as -- does not treat every number the same way.  Since
// 3.12 it sums a *float* item with Neumaier compensation, and an *int* item with
// a plain `+=` that leaves the compensation untouched.  Both are readable in the
// answer: on a 300-residue protein under the Engelman scale, whose C, G and H are
// written as integers, the reference returns 1.3170000000000002 where a
// compensated sum over the same twenty numbers returns 1.317.
//
// Ten of the twenty-eight gravy scales carry such integers -- Engelman, Fasman,
// Fauchere, GoldSack, Jones, Parker, Ponnuswamy, Roseman, Wilson and Zimmerman --
// and on four of those ten the distinction shows: over 200,000 random proteins of
// one to fifty-nine residues, a compensated sum over the same numbers answers
// differently 8,167 times under Roseman, 11,413 under Engelman, 11,788 under
// GoldSack and 58,541 under Parker, and never once under the other six, whose
// single integer value is one whose double form rounds the same way either way.
// None of the weight tables `molecular_weight` uses has an integer in it, so
// there every entry is `kCompensated`.
//
//   kAbsent       0   not in the table -- the decline below
//   kCompensated  1   a float item: Neumaier-compensated
//   kPlain        2   an int item: `total += table[byte]`, compensation untouched
//
// **A kernel counts, Python decides**, and the third thing `sum()` does is the
// same kind of fact.  It starts with an integer accumulator and stays in it
// until the first float arrives; *that* item is added to the accumulator plainly,
// with no compensation to correct.  For "IYAR" under the Parker scale -- -8,
// -1.9, 2.1, 4.2 -- the reference answers -3.6, and compensating that first pair
// answers -3.5999999999999996.  The kernel below carries the same phase, so the
// only place it can part company with the reference is where an integer
// accumulator and a double one can: a run of int items at the head of the
// sequence whose partial sums pass 2**53.  The largest value in any scale this
// package ships is 12.3, so that is a sequence of about 7e14 residues, and a
// buffer cannot be that long on any machine that has run the tests in this
// repository.
//
// Returns false and sets `bad` to the first byte that is not in the table, which
// is the reference's `KeyError` and whose message only Python can spell:
// molecular_weight's caller re-derives `f"'{e}' is not a valid unambiguous
// letter for {seq_type}"`, quotes and all.
//
// The Neumaier compensation is not a flourish.  On this interpreter a plain
// running total disagrees with `sum()` on roughly a quarter of all inputs --
// measured, not assumed: 5,337 mismatches in 20,000 random draws of four-to-twelve
// terms from the DNA weight table.  A kernel that adds naively is a kernel that
// reproduces `molecular_weight` three times out of four, which is worse than one
// that does not claim to.
constexpr unsigned char kAbsent = 0;
constexpr unsigned char kCompensated = 1;
constexpr unsigned char kPlain = 2;

bool molecular_weight_mass(const char* data, std::size_t size,
                           const double* table, const unsigned char* valid,
                           double* total, unsigned char* bad);

// ---------------------------------------------------------------------------
// The two CheckSum kernels that are loops over the sequence.
//
// Both are defined by the reference on *code points*, not on bytes: `gcg` calls
// `char.upper()` and `crc64` calls `ord(c)`, so a non-ASCII letter has to go
// through Python's Unicode tables, and `gcg("ß")` raises `TypeError: ord()
// expected a character, but string of length 2 found` because `"ß".upper()` is
// two characters.  Reproducing that in C++ would be reproducing the Unicode
// database, so these two **decline** on any byte at or above 0x80 and the caller
// runs the reference's own loop instead.  On ASCII -- which is every sequence,
// both cases, plus '-' and '*' -- `char.upper()` is `toupper` and `ord(c) & 0xFF`
// is the byte, exactly.
bool gcg(const char* data, std::size_t size, std::uint32_t* checksum);
void crc64(const char* data, std::size_t size, std::uint32_t* high,
           std::uint32_t* low);

// The 256-entry high table crc64 is built on, derived at compile time from the
// reference's own recurrence rather than dumped from it.  Exposed so that a
// test can compare every entry against `CheckSum._table_h`, which is what makes
// "derived" an acceptable answer to "where did this table come from".
std::vector<std::uint32_t> crc64_table_h();

}  // namespace biofasting
