// The per-residue kernels of `ProteinAnalysis`: the profile, and the two sums.
//
// sequtils.cpp holds the kernels that look at a sequence and return a number
// about the whole of it.  This file holds the three that are loops over every
// *residue* -- `protein_scale`, `flexibility` and `instability_index` -- and they
// are separated because they are the only place in this package where the
// arithmetic is a running accumulation whose *order* is part of the answer.
//
// That is the whole difficulty of this file, and it is worth stating before any
// code.  `flexibility` accumulates four mirror pairs, `(front + back) * weight`,
// and then adds the middle residue a second time; `protein_scale` accumulates
// `weight * front + weight * back` instead, which is the same number on paper
// and a different double in IEEE-754 -- measured, 200 of 200 proteins differ in
// the last bit.  `instability_index` adds its dipeptide weights with a plain
// `+=` loop while `molecular_weight` and `gravy` are written as `sum(...)` and
// are therefore Neumaier-compensated by CPython 3.12 and later.  Two functions in
// the same class need two different summations, and neither can be borrowed from
// the other.  So each kernel below reproduces the reference's own sequence of
// operations, and the comment says which sequence that is.
//
// **A kernel counts, Python decides** still holds, in the form it takes here:
// every one of these three declines rather than raising, because the reference
// fails in ways only Python can spell.  `flexibility` and `instability_index`
// raise `KeyError` naming the residue, and the residue Python names is the first
// one *read*, not the first one in the sequence; `protein_scale` does not raise
// at all, it writes a warning to stderr and drops that window's term --
// reproducing which from a compiled kernel would mean reproducing Python's
// stderr.  So all three hand the sequence back to the reference's own loop.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace biofasting {

// ---------------------------------------------------------------------------
// ProteinAnalysis.protein_scale -- the sliding profile.
//
// One score per window that fits, `size - window + 1` of them, each divided by
// `sum_of_weights`.  The weights are **passed in**, half a window of them,
// because the reference computes them as `edge + unit * i` with a `unit` built
// from the edge value, and recomputing that here would be a second place for a
// rounding difference to enter.  Python computes them with the reference's own
// expression and hands over the result.
//
// `table` is 256 doubles indexed by byte and `valid` is 256 flags, the same
// dense-array-of-a-sparse-dict shape `molecular_weight_mass` takes and for the
// same reason -- here a caller's scale is an arbitrary mapping, so there is no
// fixed letter order to index it by.  Any byte the caller's scale has no entry
// for is a decline: the reference would warn and skip that window's term, and
// only Python can print that.
//
// The accumulation is `score += weights[j] * front + weights[j] * back`, in that
// order, and then `score += table[middle]` with `middle` at `window / 2`.  For an
// odd window that covers every residue of the window exactly once; for an even
// one the residue at `window / 2` is counted twice, which is what the reference
// does and therefore what this does.
//
// `window` must be at least 2: at 1 the reference divides by `window - 1` while
// building its weights, so it raises before this loop ever runs.
bool protein_scale_scores(const char* data, std::size_t size, std::size_t window,
                          const double* weights, const double* table,
                          const unsigned char* valid, double sum_of_weights,
                          std::vector<double>& out);

// ---------------------------------------------------------------------------
// ProteinAnalysis.flexibility -- protein_scale's window of nine, hard-coded.
//
// It is not a special case of the kernel above, and the difference is the kind
// that produces a plausible wrong number:
//
//   * the weights are fixed at 0.25, 0.4375, 0.625 and 0.8125 against mirror
//     pairs two apart from the ends, and the middle term is read at
//     `window / 2 + 1` -- index 5 of a nine-residue window -- not at index 4;
//   * so index 4 of every window is never read at all, and index 5 is read
//     twice, once as the back of the last pair and once as the middle;
//   * the accumulation groups as `(front + back) * weight`, not
//     `weight * front + weight * back`;
//   * the count is `size - 9`, one window fewer than fit, which is the
//     reference's own off-by-one and is reproduced here rather than fixed;
//   * and the divisor is the literal 5.25.
//
// `flex` is 20 doubles in `AMINO_ACIDS` order and `map` is as above.  Any 0xFF
// is a decline, for `KeyError`'s sake.
bool flexibility_scores(const char* data, std::size_t size, const double* flex,
                        const unsigned char* map, std::vector<double>& out);

// ---------------------------------------------------------------------------
// ProteinAnalysis.instability_index's sum.
//
// `size - 1` dipeptides, `table[a * 20 + b]` for each adjacent pair, added with
// a plain `+=` in sequence order.  **Not** Neumaier-compensated, because the
// reference is not: it writes `score += dipeptide_value` in a loop, where
// `molecular_weight` and `gravy` write `sum(...)`.  The caller multiplies by
// `10.0 / length` afterwards, which is where the empty sequence raises.
//
// `table` is 400 doubles, row-major in `AMINO_ACIDS`.  Returns false on any byte
// the map has no position for.
bool instability_index_sum(const char* data, std::size_t size,
                           const double* table, const unsigned char* map,
                           double* total);

// ---------------------------------------------------------------------------
// ProteinAnalysis.count_amino_acids -- how many of each standard residue.
//
// The one kernel in this file with no arithmetic order to get wrong, and the
// reason it exists is a measurement rather than a subtlety.  The reference is
// twenty `str.count` calls, so it walks the sequence twenty times and pays
// twenty Python-level calls to do it; a residue count is one pass and one
// bucket per byte.  Measured on the generated protein corpus before this was
// written: 3.2 us per 150-residue protein and 178 us per 10,000-residue one,
// against 0.8 us and 16.5 us for the same count done in one pass.
//
// It matters out of proportion to its own row, because it is the floor under
// seven others: `charge_at_pH`, `isoelectric_point`, `amino_acids_percent`,
// `aromaticity`, `secondary_structure_fraction` and
// `molar_extinction_coefficient` all begin by calling it, and on a fresh object
// -- which is what a caller pays for -- that call is most of what they cost.
//
// **This kernel never declines.**  The three above do, because the reference
// fails in a way only Python can spell; here the reference does not fail at
// all.  `str.count('A')` looks for one letter and is uninterested in what it
// does not find, so a residue outside the twenty is a zero and not an error.
//
// `counts` is twenty entries in `AMINO_ACIDS` order.  `map` is 256 bytes holding
// each byte's position in that order, or 0xFF for a byte that has none.
void count_residues(const char* data, std::size_t size,
                    const unsigned char* map, std::uint64_t* counts);

}  // namespace biofasting
