#include "seqops.hpp"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "cpu_features.hpp"

// The AVX2 kernel is guarded on three things at once: the target has to be
// x86-64, the compiler has to support the function target attribute, and it
// must not be MSVC -- where the attribute does not exist and a separate
// compilation unit would be needed.  Where any of that is missing the kernel is
// not compiled and the ladder cannot select it, so the two halves of that
// decision cannot drift apart.
#if (defined(__x86_64__) || defined(_M_X64)) && \
    (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
#  define BIOFASTING_HAVE_X86_AVX2 1
#  include <immintrin.h>
#endif

namespace biofasting {
namespace {

// ---------------------------------------------------------------------------
// The scalar tables.  Both are the contract, not an optimisation: the AVX2
// kernel below is only allowed to run where it reproduces these exactly.
// ---------------------------------------------------------------------------

// Biopython's _dna_complement_table, dumped from the reference rather than
// written from memory: the IUPAC pairs in both cases, and identity everywhere
// else, which is what makes '-' stay '-' and 'N' stay 'N'.
constexpr std::array<char, 256> make_complement_table() noexcept {
  std::array<char, 256> table{};
  for (std::size_t i = 0; i < table.size(); ++i) {
    table[i] = static_cast<char>(i);
  }
  constexpr char from[] = "ABCDGHKMRTUVYabcdghkmrtuvy";
  constexpr char to[] = "TVGHCDMKYAABRtvghcdmkyaabr";
  static_assert(sizeof(from) == sizeof(to), "the two halves must pair up");
  for (std::size_t i = 0; i + 1 < sizeof(from); ++i) {
    table[static_cast<unsigned char>(from[i])] = to[i];
  }
  return table;
}

constexpr std::array<char, 256> kComplement = make_complement_table();

// gc_fraction's categories: 1 is GC, 2 is the rest of the denominator, 0 is
// removed from both.  Written out of SeqUtils.gc_fraction's definition, which
// is `gc = count("CGScgs")` and, for the default ambiguous="remove",
// `length = gc + count("ATWUatwu")`.
constexpr std::array<std::uint8_t, 256> make_gc_category() noexcept {
  std::array<std::uint8_t, 256> table{};
  for (char c : {'C', 'G', 'S', 'c', 'g', 's'}) {
    table[static_cast<unsigned char>(c)] = 1;
  }
  for (char c : {'A', 'T', 'W', 'U', 'a', 't', 'w', 'u'}) {
    table[static_cast<unsigned char>(c)] = 2;
  }
  return table;
}

constexpr std::array<std::uint8_t, 256> kGcCategory = make_gc_category();

// The AVX2 kernel's hash: (c >> 1) & 0x0F.  Shifting right by one moves the
// ASCII case bit (0x20) into bit 4, where the mask throws it away -- so the
// same index serves 'A' and 'a', and the five letters ACGTN land on five
// distinct indexes even though 'T' is 0x54 and 'A' is 0x41.  Anything the table
// does not name is caught by the equality check and the kernel declines.
constexpr std::array<char, 16> make_acgtn_valid() noexcept {
  std::array<char, 16> table{};
  for (char c : {'A', 'C', 'G', 'T', 'N'}) {
    table[(static_cast<unsigned char>(c) >> 1) & 0x0F] = c;
  }
  return table;
}

constexpr std::array<char, 16> make_acgtn_complement() noexcept {
  std::array<char, 16> table{};
  for (char c : {'A', 'C', 'G', 'T', 'N'}) {
    table[(static_cast<unsigned char>(c) >> 1) & 0x0F] =
        kComplement[static_cast<unsigned char>(c)];
  }
  return table;
}

constexpr std::array<char, 16> kAcgtnValid = make_acgtn_valid();
constexpr std::array<char, 16> kAcgtnComplement = make_acgtn_complement();

// ---------------------------------------------------------------------------
// The two scalar implementations.
// ---------------------------------------------------------------------------

void complement_iupac(const char* data, std::size_t size, char* out) noexcept {
  // Filled back to front, so the reversal is free and there is no second pass
  // over the result.
  char* end = out + size;
  for (std::size_t i = 0; i < size; ++i) {
    *--end = kComplement[static_cast<unsigned char>(data[i])];
  }
}

// Adds the GC and AT counts of `size` bytes to `gc` and `at`.  This is both
// the whole implementation on a machine without AVX2 and the tail of the
// vector one, so the two are the same code wherever they meet.
//
// Branchless on purpose, and not for style.  Written as
// `if (category == 1) ++gc; else if (category == 2) ++at;` this runs at 241
// MB/s: real DNA alternates between the two categories at random, so every
// byte is a coin-flip the branch predictor loses.  Two compares and two adds
// per byte is several times faster than one unpredictable branch.
void gc_counts_scalar(const char* data, std::size_t size, std::size_t& gc,
                      std::size_t& at) noexcept {
  for (std::size_t i = 0; i < size; ++i) {
    const std::uint8_t category =
        kGcCategory[static_cast<unsigned char>(data[i])];
    gc += static_cast<std::size_t>(category == 1);
    at += static_cast<std::size_t>(category == 2);
  }
}

int base_code(char c) noexcept {
  switch (c) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      return -1;
  }
}

#ifdef BIOFASTING_HAVE_X86_AVX2

// Complements an ACGTN sequence into `out`, or returns false without a
// complete answer if any byte is outside that alphabet -- in which case the
// caller redoes the whole thing with the IUPAC table.  Declining rather than
// guessing is the same rule the FASTQ fast path follows.
//
// Each 32-byte block is complemented and byte-reversed in registers and then
// stored at the *mirrored* position, so the reversal is not a second pass over
// the result.  The first version did complement-then-reverse-in-place and was
// 1.3x the scalar path: it read and wrote the whole sequence twice to do one
// sequence's worth of work, and at 40 MB the loop was already waiting on
// memory.
__attribute__((target("avx2"))) bool reverse_complement_acgtn_avx2(
    const char* data, std::size_t size, char* out) {
  const __m256i comp_lut = _mm256_broadcastsi128_si256(
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(kAcgtnComplement.data())));
  const __m256i valid_lut = _mm256_broadcastsi128_si256(
      _mm_loadu_si128(reinterpret_cast<const __m128i*>(kAcgtnValid.data())));
  const __m256i nibble = _mm256_set1_epi8(0x0F);
  const __m256i case_bit = _mm256_set1_epi8(0x20);
  const __m256i reverse_16 = _mm256_setr_epi8(
      15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11,
      10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);

  // Blocks are consumed from the front and stored mirrored, so the block read
  // at [front, front+32) lands at [size-front-32, size-front): the far end of
  // the output.  The scalar loop below then covers input [front, size) -- the
  // near end of the output -- and the two meet in the middle with nothing
  // either side of the join left unwritten.  The earlier bound here (a
  // "middle" expression in terms of front) skipped the input tail whenever size
  // was close to a multiple of 64 and left the first bytes of the output as
  // whatever the string was allocated with; it passed every large-n test and
  // failed first at n = 64.
  std::size_t front = 0;
  while (front + 32 <= size) {
    const __m256i c =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + front));
    // The 16-bit shift drags a bit of the neighbouring byte into bit 7 of each
    // byte; the 0x0F mask is what removes it, so no per-byte fixup is needed.
    const __m256i index = _mm256_and_si256(_mm256_srli_epi16(c, 1), nibble);
    const __m256i upper = _mm256_andnot_si256(case_bit, c);
    const __m256i valid = _mm256_shuffle_epi8(valid_lut, index);
    if (_mm256_movemask_epi8(_mm256_cmpeq_epi8(valid, upper)) != -1) {
      return false;
    }
    __m256i complement = _mm256_shuffle_epi8(comp_lut, index);
    complement = _mm256_or_si256(complement, _mm256_and_si256(c, case_bit));
    // Reverse within each 128-bit lane, then swap the lanes: together that is
    // a byte reversal of all 32.
    complement = _mm256_shuffle_epi8(complement, reverse_16);
    complement = _mm256_permute2x128_si256(complement, complement, 0x01);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(out + size - front - 32), complement);
    front += 32;
  }

  // The tail uses the same hash and the same tables as the scalar path, so the
  // two cannot disagree about what ACGTN means.  It runs to `size`, not to
  // `size - front`: the blocks above consumed the front of the input, and the
  // output they wrote sits at the far end, so every remaining input byte is
  // this loop's.
  for (std::size_t i = front; i < size; ++i) {
    const unsigned char c = static_cast<unsigned char>(data[i]);
    const unsigned char index = (c >> 1) & 0x0F;
    if (static_cast<unsigned char>(kAcgtnValid[index]) != (c & ~0x20u)) {
      return false;
    }
    out[size - 1 - i] = static_cast<char>(kAcgtnComplement[index] | (c & 0x20u));
  }
  return true;
}

// Sums the four 64-bit lanes of a byte-histogram vector.  The accumulator
// holds one small count per byte lane, so `sad` (sum of absolute differences
// against zero) reduces each 8-byte group in one instruction.
// The target attribute is not decoration: AVX2 intrinsics are `always_inline`
// and refuse to inline into a function the compiler has not been told may use
// them, so a helper without it fails to compile rather than silently running
// slow.
__attribute__((target("avx2"))) inline std::uint64_t sum_byte_counts(
    __m256i counts) noexcept {
  const __m256i sums = _mm256_sad_epu8(counts, _mm256_setzero_si256());
  const __m128i halves =
      _mm_add_epi64(_mm256_castsi256_si128(sums), _mm256_extracti128_si256(sums, 1));
  return static_cast<std::uint64_t>(_mm_cvtsi128_si64(halves)) +
         static_cast<std::uint64_t>(_mm_extract_epi64(halves, 1));
}

// Counts GC and AT bytes 32 at a time, and returns how many bytes it consumed
// -- always 32-aligned, with the tail left to the scalar path.
//
// Unlike the reverse-complement kernel this one never declines, and that is
// worth being sure of rather than hoping: every byte is folded with `| 0x20`
// and then compared against the seven lower-case letters, and `x | 0x20` equals
// a target letter only for x in {letter, letter & ~0x20} -- which is exactly the
// case pair the scalar table names.  There is no byte it classifies differently
// from the table, so the two cannot disagree, and the differential test sweeps
// all 256 of them anyway.
//
// Measured before it was written: the scalar loop runs at ~1.8 GB/s and gets
// the same number whether the sequence is 1 MB or 40 MB, so it is not waiting
// on memory -- it is doing 2.4 cycles per byte of table lookups.
__attribute__((target("avx2"))) std::size_t gc_counts_avx2(
    const char* data, std::size_t size, std::size_t& gc, std::size_t& at) noexcept {
  const __m256i case_bit = _mm256_set1_epi8(0x20);
  const __m256i gc_c = _mm256_set1_epi8('c');
  const __m256i gc_g = _mm256_set1_epi8('g');
  const __m256i gc_s = _mm256_set1_epi8('s');
  const __m256i at_a = _mm256_set1_epi8('a');
  const __m256i at_t = _mm256_set1_epi8('t');
  const __m256i at_w = _mm256_set1_epi8('w');
  const __m256i at_u = _mm256_set1_epi8('u');

  // Each compare produces 0xFF where it matches, so subtracting adds one to
  // that lane's count.  A lane can gain at most 3 (GC) or 4 (AT) per block, so
  // 16 blocks cannot overflow the byte: 48 and 64 against a limit of 255.
  constexpr std::size_t kBlocksPerFlush = 16;
  constexpr std::size_t kChunk = 32 * kBlocksPerFlush;

  std::size_t i = 0;
  for (; i + kChunk <= size; i += kChunk) {
    __m256i gc_acc = _mm256_setzero_si256();
    __m256i at_acc = _mm256_setzero_si256();
    for (std::size_t b = 0; b < kBlocksPerFlush; ++b) {
      const __m256i folded = _mm256_or_si256(
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + i + 32 * b)),
          case_bit);
      gc_acc = _mm256_sub_epi8(gc_acc, _mm256_cmpeq_epi8(folded, gc_c));
      gc_acc = _mm256_sub_epi8(gc_acc, _mm256_cmpeq_epi8(folded, gc_g));
      gc_acc = _mm256_sub_epi8(gc_acc, _mm256_cmpeq_epi8(folded, gc_s));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_a));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_t));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_w));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_u));
    }
    gc += static_cast<std::size_t>(sum_byte_counts(gc_acc));
    at += static_cast<std::size_t>(sum_byte_counts(at_acc));
  }

  // Fewer than 16 blocks left: same shape, one flush at the end.
  if (i + 32 <= size) {
    __m256i gc_acc = _mm256_setzero_si256();
    __m256i at_acc = _mm256_setzero_si256();
    for (; i + 32 <= size; i += 32) {
      const __m256i folded = _mm256_or_si256(
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data + i)), case_bit);
      gc_acc = _mm256_sub_epi8(gc_acc, _mm256_cmpeq_epi8(folded, gc_c));
      gc_acc = _mm256_sub_epi8(gc_acc, _mm256_cmpeq_epi8(folded, gc_g));
      gc_acc = _mm256_sub_epi8(gc_acc, _mm256_cmpeq_epi8(folded, gc_s));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_a));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_t));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_w));
      at_acc = _mm256_sub_epi8(at_acc, _mm256_cmpeq_epi8(folded, at_u));
    }
    gc += static_cast<std::size_t>(sum_byte_counts(gc_acc));
    at += static_cast<std::size_t>(sum_byte_counts(at_acc));
  }
  return i;
}

#endif  // BIOFASTING_HAVE_X86_AVX2

#ifdef BIOFASTING_HAVE_X86_AVX2

// Whether the AVX2 kernel may run.  Resolved once: the answer cannot change
// while the process runs, and the caller is on the hot path -- which is also
// why this is a cached bool rather than a comparison of the level name, since
// a per-call std::string would cost more than the operation does on a 150 bp
// read.
//
// BIOFASTING_SEQOPS_LEVEL=baseline pins the dispatch to the scalar path.  It
// exists because a benchmark has to be able to measure the thing it is
// claiming to beat, and because a CPU whose AVX2 this build mis-detects needs
// a way out that is not a recompile.  It only ever *lowers* the rung: forcing
// "avx2" on a CPU without it is not a thing this will do, so the variable
// cannot turn a working install into SIGILL.
bool avx2_usable() {
  static const bool usable = [] {
    if (const char* forced = std::getenv("BIOFASTING_SEQOPS_LEVEL")) {
      if (std::string_view(forced) == "baseline") return false;
    }
    for (const std::string& name : supported_cpu_levels()) {
      if (name == "avx2") return true;
    }
    return false;
  }();
  return usable;
}

#endif  // BIOFASTING_HAVE_X86_AVX2

}  // namespace

std::string seqops_level() {
#ifdef BIOFASTING_HAVE_X86_AVX2
  if (avx2_usable()) return std::string("avx2");
#endif
  // "baseline" is not a fallback in the apologetic sense: it is a promise
  // that every target of this wheel can run it.
  return std::string("baseline");
}

std::string reverse_complement(const char* data, std::size_t size) {
  std::string out(size, '\0');
  if (size == 0) return out;

#ifdef BIOFASTING_HAVE_X86_AVX2
  if (avx2_usable() && reverse_complement_acgtn_avx2(data, size, out.data())) {
    return out;
  }
#endif
  complement_iupac(data, size, out.data());
  return out;
}

double gc_fraction(const char* data, std::size_t size) {
  std::size_t gc = 0;
  std::size_t at = 0;
  std::size_t done = 0;
#ifdef BIOFASTING_HAVE_X86_AVX2
  if (avx2_usable()) done = gc_counts_avx2(data, size, gc, at);
#endif
  gc_counts_scalar(data + done, size - done, gc, at);
  const std::size_t length = gc + at;
  //   -1.0 and not 0.0 when nothing was counted, because the reference's
  //   `return 0` there is the *integer* zero and a `double` cannot say so.  A
  //   real fraction is never negative, so the caller -- which can spell the int
  //   -- is the right place to make that distinction, and this is how it is
  //   told that it needs to.
  if (length == 0) return -1.0;
  return static_cast<double>(gc) / static_cast<double>(length);
}

std::vector<std::pair<std::string, std::uint64_t>> count_kmers(
    const char* data, std::size_t size, unsigned k) {  if (k < 1 || k > 16) {
    throw std::invalid_argument("k must be between 1 and 16, got " +
                                std::to_string(k));
  }

  // Two bits per base, so a window of at most 16 fits a 32-bit code and the
  // sliding update is one shift and one mask.
  const std::uint64_t mask = (std::uint64_t{1} << (2 * k)) - 1;
  std::unordered_map<std::uint64_t, std::uint64_t> counts;
  std::uint64_t code = 0;
  std::size_t run = 0;
  for (std::size_t i = 0; i < size; ++i) {
    const int base = base_code(data[i]);
    if (base < 0) {
      // A window containing anything that is not ACGT is not a k-mer; the run
      // has to restart, or N would silently join the letters on either side.
      code = 0;
      run = 0;
      continue;
    }
    code = ((code << 2) | static_cast<std::uint64_t>(base)) & mask;
    if (++run >= k) ++counts[code];
  }

  std::vector<std::pair<std::string, std::uint64_t>> out;
  out.reserve(counts.size());
  constexpr char bases[] = "ACGT";
  for (const auto& [value, count] : counts) {
    std::string text(k, 'A');
    for (unsigned position = 0; position < k; ++position) {
      // Most significant base first: the code was built left to right.
      text[position] = bases[(value >> (2 * (k - 1 - position))) & 0x3];
    }
    out.emplace_back(std::move(text), count);
  }
  std::sort(out.begin(), out.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  return out;
}

TranslateResult translate(const char* data, std::size_t size, const char* code,
                          const char* amino, std::string_view stop_symbol,
                          std::string_view possible_stop, bool to_stop,
                          bool stop_is_error) {
  TranslateResult result;
  const std::size_t codons = size / 3;
  result.protein.reserve(codons + stop_symbol.size());

  // One pass, no lookahead, and no per-codon allocation: the three table
  // lookups and the shift-or that combines them are the whole operation.  The
  // base codes are two bits each, so the codon index is the twelve-bit
  // `a << 4 | b << 2 | c` and the amino table is exactly 4 KiB of nothing.
  for (std::size_t i = 0; i < codons; ++i) {
    const char* triplet = data + 3 * i;
    const auto a = static_cast<unsigned char>(code[static_cast<unsigned char>(triplet[0])]);
    const auto b = static_cast<unsigned char>(code[static_cast<unsigned char>(triplet[1])]);
    const auto c = static_cast<unsigned char>(code[static_cast<unsigned char>(triplet[2])]);
    // 0x80 is the one bit no base code sets, so this is three ORs and a test
    // rather than three branches -- and it is what makes "the kernel refuses
    // anything but this genetic code's own four letters" a property of the
    // table rather than of a character class written here.
    if ((a | b | c) & 0x80) {
      result.outcome = TranslateOutcome::kDeclined;
      result.codon = i;
      return result;
    }
    const auto entry =
        static_cast<unsigned char>(amino[(a << 4) | (b << 2) | c]);
    if (entry >= 3) {
      result.protein.push_back(static_cast<char>(entry));
    } else if (entry == kTranslateStop) {
      if (stop_is_error) {
        result.outcome = TranslateOutcome::kExtraStop;
        result.codon = i;
        return result;
      }
      // to_stop ends the protein here, which is a break and not a return: the
      // codons after it are never examined, exactly as in the reference, so a
      // stop hides nothing that would have been refused.
      if (to_stop) break;
      result.protein.append(stop_symbol);
    } else if (entry == kTranslatePossibleStop) {
      result.protein.append(possible_stop);
    } else {
      // kTranslateUnclassified: the table has no amino acid here and this is
      // not a stop, so the reference raises TranslationError naming the codon.
      result.outcome = TranslateOutcome::kDeclined;
      result.codon = i;
      return result;
    }
  }
  return result;
}

}  // namespace biofasting
