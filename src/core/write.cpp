#include "write.hpp"

#include <algorithm>
#include <cstdint>

namespace biofasting {
namespace {

// The QUAL line cut, in the reference's own terms: `data.rfind(" ", 0, wrap)`
// is the last space at an index below `wrap`.
constexpr char kSpace = ' ';

// `SANGER_SCORE_OFFSET`: a Sanger FASTQ stores a PHRED score as `score + 33`,
// so the inverse is the whole of what a QUAL record has to undo.
constexpr unsigned char kSangerOffset = 33;

// The last character a Sanger FASTQ may use, and so the truncation point for a
// score the format cannot hold.
constexpr unsigned char kSangerMax = 126;

// One PHRED score as decimal, without a format call.  A score is a byte, so it
// has at most three digits and the branches are the whole of it -- the
// reference spends a Python-level `"%i" % round(q, 0)` on every score, which is
// what makes its QUAL writer cost 103 ns a base.
void append_score(std::string& out, unsigned char score) {
  if (score >= 100) {
    out.push_back(static_cast<char>('0' + score / 100));
    out.push_back(static_cast<char>('0' + (score / 10) % 10));
    out.push_back(static_cast<char>('0' + score % 10));
  } else if (score >= 10) {
    out.push_back(static_cast<char>('0' + score / 10));
    out.push_back(static_cast<char>('0' + score % 10));
  } else {
    out.push_back(static_cast<char>('0' + score));
  }
}

}  // namespace

void append_fasta(std::string& out, std::string_view title,
                  std::string_view sequence, std::size_t wrap) {
  out.push_back('>');
  out.append(title.data(), title.size());
  out.push_back('\n');

  if (wrap == 0) {
    // `if self.wrap:` is false, so the reference writes `data + "\n"` -- the
    // whole sequence on one line, and one newline even for no sequence.
    out.append(sequence.data(), sequence.size());
    out.push_back('\n');
    return;
  }

  // A wrapped record with an empty sequence writes **nothing** after the title:
  // `for i in range(0, 0, 60)` never runs.  A reader that wrote a blank line
  // would be one byte wrong on a legal record.
  for (std::size_t i = 0; i < sequence.size(); i += wrap) {
    const std::size_t take = std::min(wrap, sequence.size() - i);
    out.append(sequence.data() + i, take);
    out.push_back('\n');
  }
}

void append_fastq(std::string& out, std::string_view title,
                  std::string_view sequence, std::string_view quality) {
  out.push_back('@');
  out.append(title.data(), title.size());
  out.push_back('\n');
  out.append(sequence.data(), sequence.size());
  out.append("\n+\n", 3);
  out.append(quality.data(), quality.size());
  out.push_back('\n');
}

bool append_qual(std::string& out, std::string_view title,
                 std::string_view quality, std::size_t wrap) {
  out.push_back('>');
  out.append(title.data(), title.size());
  out.push_back('\n');

  // The scores as decimal, joined by single spaces: `" ".join(qualities_strs)`.
  // The reference builds this list with a Python format call per score and then
  // joins it; here it is the same characters built once, and the line cut then
  // works on the string rather than on a list that is popped from the front.
  std::string data;
  data.reserve(quality.size() * 4);
  for (std::size_t i = 0; i < quality.size(); ++i) {
    if (i != 0) data.push_back(kSpace);
    // Our quality is the ASCII string a reader produced, so the score is the
    // byte less the offset.  A byte below the offset is not a quality at all
    // and the caller has already rejected it.
    append_score(data, static_cast<unsigned char>(quality[i] - kSangerOffset));
  }

  if (wrap == 0) {
    out.append(data);
    out.push_back('\n');
    return true;
  }

  while (true) {
    if (data.size() <= wrap) {
      out.append(data);
      out.push_back('\n');
      return true;
    }
    const std::size_t cut = data.rfind(kSpace, wrap - 1);
    if (cut == std::string::npos) {
      // `data.rfind(" ", 0, wrap)` returned -1, and the reference then writes
      // `data[:-1]` and slices `data[0:]` -- the same string back, forever.
      // This is a hang and not a line, so there is no output to reproduce.
      return false;
    }
    out.append(data, 0, cut);
    out.push_back('\n');
    data.erase(0, cut + 1);
  }
}

void phred_to_sanger(const unsigned char* scores, std::size_t count,
                     std::string& out) {
  // `min(126, qp + SANGER_SCORE_OFFSET)`: PHRED 93 is the highest score a
  // Sanger FASTQ can hold and 126 (`~`) the last character it can use, so a
  // score above 93 is truncated rather than wrapped -- which is what the
  // reference does on the slow path it takes for any score its 0..93 table has
  // no entry for.
  out.resize(count);
  for (std::size_t i = 0; i < count; ++i) {
    const unsigned int encoded = static_cast<unsigned int>(scores[i]) +
                                 static_cast<unsigned int>(kSangerOffset);
    out[i] = static_cast<char>(encoded > kSangerMax ? kSangerMax : encoded);
  }
}

}  // namespace biofasting
