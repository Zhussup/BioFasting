#include "protparam.hpp"

namespace biofasting {

namespace {

// The position of a residue in the caller's table, or a decline.
//
// Scanning for an unmapped byte before doing any arithmetic is deliberate: the
// caller re-runs the whole call in Python when this returns false, so a kernel
// that computed half a profile and then gave up would have bought nothing but a
// half-written vector.
inline bool all_mapped(const char* data, std::size_t size,
                       const unsigned char* map) {
  for (std::size_t i = 0; i < size; ++i) {
    if (map[static_cast<unsigned char>(data[i])] == 0xFF) {
      return false;
    }
  }
  return true;
}

// The same, for the dense tables that carry their own flag array.
inline bool all_valid(const char* data, std::size_t size,
                      const unsigned char* valid) {
  for (std::size_t i = 0; i < size; ++i) {
    if (!valid[static_cast<unsigned char>(data[i])]) {
      return false;
    }
  }
  return true;
}

}  // namespace

bool protein_scale_scores(const char* data, std::size_t size, std::size_t window,
                          const double* weights, const double* table,
                          const unsigned char* valid, double sum_of_weights,
                          std::vector<double>& out) {
  if (window < 2) {
    //   `window == 1` divides by zero in the reference while it builds its
    //   weights, before the loop that this function is; `window == 0` reaches
    //   `subsequence[0]` on an empty slice.  Both are Python's to reproduce.
    return false;
  }
  if (!all_valid(data, size, valid)) {
    return false;
  }
  const std::size_t half = window / 2;
  //   `range(size - window + 1)`, which is empty rather than negative when the
  //   window does not fit.
  const std::size_t windows = size + 1 > window ? size - window + 1 : 0;
  out.clear();
  out.reserve(windows);
  for (std::size_t start = 0; start < windows; ++start) {
    const char* const window_data = data + start;
    double score = 0.0;
    for (std::size_t j = 0; j < half; ++j) {
      const double front = table[static_cast<unsigned char>(window_data[j])];
      const double back =
          table[static_cast<unsigned char>(window_data[window - j - 1])];
      //   Two separate `weight * value` products added together, which is the
      //   reference's grouping and not a rearrangement of it.
      score += weights[j] * front + weights[j] * back;
    }
    score += table[static_cast<unsigned char>(window_data[half])];
    out.push_back(score / sum_of_weights);
  }
  return true;
}

bool flexibility_scores(const char* data, std::size_t size, const double* flex,
                        const unsigned char* map, std::vector<double>& out) {
  if (!all_mapped(data, size, map)) {
    return false;
  }
  constexpr std::size_t kWindow = 9;
  constexpr std::size_t kHalf = 4;
  constexpr double kWeights[kHalf] = {0.25, 0.4375, 0.625, 0.8125};
  //   Index 5, not 4: the reference reads the middle at `window / 2 + 1`, so the
  //   residue at the true centre of the window is not read at all and the one
  //   after it is read twice.  Reproduced, not repaired -- a profile that
  //   disagreed with Biopython's would be a bug however much better it looked.
  constexpr std::size_t kMiddle = 5;
  const std::size_t windows = size > kWindow ? size - kWindow : 0;
  out.clear();
  out.reserve(windows);
  for (std::size_t start = 0; start < windows; ++start) {
    const char* const window_data = data + start;
    double score = 0.0;
    for (std::size_t j = 0; j < kHalf; ++j) {
      const double front =
          flex[map[static_cast<unsigned char>(window_data[j])]];
      const double back = flex[map[static_cast<unsigned char>(
          window_data[kWindow - j - 1])]];
      //   The pair is added first and multiplied once, which is what makes this
      //   a different function from protein_scale and not a slower spelling of
      //   it: distributing the weight over the two residues gives the same
      //   number and a different float.
      score += (front + back) * kWeights[j];
    }
    score += flex[map[static_cast<unsigned char>(window_data[kMiddle])]];
    out.push_back(score / 5.25);
  }
  return true;
}

bool instability_index_sum(const char* data, std::size_t size,
                           const double* table, const unsigned char* map,
                           double* total) {
  if (!all_mapped(data, size, map)) {
    return false;
  }
  double score = 0.0;
  for (std::size_t i = 0; i + 1 < size; ++i) {
    const std::size_t a = map[static_cast<unsigned char>(data[i])];
    const std::size_t b = map[static_cast<unsigned char>(data[i + 1])];
    //   A plain running total.  The reference adds floats here with `+=` and
    //   not with `sum()`, so this is *not* compensated, and a Neumaier sum here
    //   would agree with the reference on the proteins that do not need it and
    //   differ on the ones that do.
    score += table[a * 20 + b];
  }
  *total = score;
  return true;
}

void count_residues(const char* data, std::size_t size,
                    const unsigned char* map, std::uint64_t* counts) {
  //   One pass, one load of the map and one increment per byte.  Twenty
  //   separate `str.count` calls -- what the reference does -- read the whole
  //   sequence twenty times over, and that is the entire difference.
  //
  //   The histogram is sized 21 and not 20 so that the unmapped bucket has
  //   somewhere to go: a byte outside the twenty still has to be looked up, and
  //   dropping it with a branch would make the loop's cost depend on the input
  //   in a way the branch predictor learns and the reader does not expect.  It
  //   is counted and then not copied out.
  std::uint64_t buckets[21] = {};
  for (std::size_t i = 0; i < size; ++i) {
    const unsigned char position = map[static_cast<unsigned char>(data[i])];
    ++buckets[position < 20 ? position : 20];
  }
  for (std::size_t a = 0; a < 20; ++a) counts[a] = buckets[a];
}

}  // namespace biofasting
