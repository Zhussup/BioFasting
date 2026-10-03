// Sequence operations, and the first kernel that actually dispatches.
//
// cpu_features.cpp has carried the ISA ladder since step 1.1 with nothing on
// it.  This is the caller: `reverse_complement` has a baseline implementation
// and an AVX2 one, and which runs is decided once, at runtime, from the ladder
// -- never by -march=native, which would make the wheel build on one machine
// and die on another.
//
// The contracts are Biopython's, because the whole point of the port is that a
// caller can swap the import and get the same answers:
//
//   * revcomp is `bytes.translate(_dna_complement_table)[::-1]`, which is the
//     full IUPAC set in both cases, and leaves every byte it does not know
//     alone -- so '-' is '-', and 'N' is 'N';
//   * gc_fraction is SeqUtils.gc_fraction with its default ambiguous="remove":
//     G, C and S count as GC, A, T, W and U are the rest of the denominator,
//     and everything else is removed from both.
//
// Which of them is worth a vector kernel was a measurement, not a judgement.
// Reverse complement and GC fraction both are, and both have one; k-mer
// counting is not, because on any input where the answer is interesting the
// time is in the hash map it builds rather than in the scan.  Where a kernel is
// not worth writing, the scalar code is the implementation, not a fallback --
// and it is the same code the vector path hands its tail to, so there is never
// a third version to keep correct.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace biofasting {

// Which rung the sequence kernels dispatch to.  Resolved once, on first use.
std::string seqops_level();

// Biopython's Seq.reverse_complement, byte for byte.
std::string reverse_complement(const char* data, std::size_t size);

// Biopython's SeqUtils.gc_fraction(seq), default arguments.
double gc_fraction(const char* data, std::size_t size);

// Overlapping k-mer counts over ACGT, k in [1, 16], case-insensitive, keys
// returned uppercase and sorted.  A window containing anything else -- N
// included -- is not a k-mer and is skipped, which is what makes the counts
// comparable with `Seq.count_overlap` per key.
//
// Throws std::invalid_argument if k is out of range; the binding turns that
// into a ValueError.
std::vector<std::pair<std::string, std::uint64_t>> count_kmers(
    const char* data, std::size_t size, unsigned k);

// ---------------------------------------------------------------------------
// Translation.
//
// This kernel is table-driven on purpose.  Which codon means which amino acid,
// and which codons are stops, is a fact about a genetic code and not about a
// loop, so it arrives from Python as two flat tables -- `code` (256 bytes: the
// 2-bit base value of every byte, or 0x80 for a byte this genetic code does not
// accept) and `amino` (64 bytes: the amino-acid letter for each of the 64
// codons, or one of the three sentinels below).  Building them is a per-table
// job done once in Python; the kernel never asks what table it is running.
enum : unsigned char {
  //: The table names no amino acid for this codon and it is not a stop --
  //: `_translate_str` would raise TranslationError.  The caller takes over.
  kTranslateUnclassified = 0,
  //: A stop codon: `stop_symbol`, or a stop of the output if `to_stop`, or an
  //: error if `stop_is_error`.
  kTranslateStop = 1,
  //: A codon every reading of which agrees except that some are stops, so the
  //: reference answers `pos_stop` ("X"); TAN is the standard example.
  kTranslatePossibleStop = 2,
};

enum class TranslateOutcome {
  kOk,
  //: The kernel met a codon it cannot classify -- 0x80 in `code`, or
  //: kTranslateUnclassified in `amino`.  `codon` says which.  Nothing is
  //: claimed about `protein`; the caller redoes the whole translation.
  kDeclined,
  //: cds=True and an in-frame stop stands at `codon`.  The reference raises
  //: TranslationError naming that codon, which only Python can spell.
  kExtraStop,
};

struct TranslateResult {
  TranslateOutcome outcome = TranslateOutcome::kOk;
  std::size_t codon = 0;
  std::string protein;
};

// Translates in frame from the start.  `stop_symbol` and `possible_stop` are
// strings rather than characters because the reference concatenates them
// without checking -- `stop_symbol="STOP"` is accepted there and produces
// "STOP" once per stop codon, and a kernel that took a char would silently
// disagree with the library it is checked against.
TranslateResult translate(const char* data, std::size_t size, const char* code,
                          const char* amino, std::string_view stop_symbol,
                          std::string_view possible_stop, bool to_stop,
                          bool stop_is_error);

}  // namespace biofasting
