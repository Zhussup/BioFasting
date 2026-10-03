// Lines and whitespace, as Python's text-mode file object defines them.
//
// This is a separate module because both kernels need exactly the same answer
// and there is no version of "we will keep them in sync" that survives a year.
// Biopython is handed a *text* stream, so its notion of a line and of trailing
// whitespace is Python's, not the file's, and a kernel that works on bytes has
// to reproduce that deliberately:
//
//   * "\n", "\r\n" and a lone "\r" all end a line, and none of them is part of
//     the content, because universal-newline translation has removed it before
//     anything else sees it;
//   * rstrip() strips the ASCII whitespace set below, which is what matters for
//     a file that is ASCII -- and every file here is.
//
// The consequence for both kernels is the same and is worth stating once: a
// line's content never contains '\r' or '\n', so no downstream code has to
// defend against them.

#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

namespace biofasting {

// The characters Python's str.rstrip() removes, for the byte values that can
// actually appear.  Python also strips U+0085 and U+00A0, but reaching those
// needs a decoded stream; the bytes that could produce them are never accepted
// as content here, so they are not part of the relevant behaviour and are not
// pretended to be.
inline bool is_python_space(char c) noexcept {
  switch (c) {
    case ' ':
    case '\t':
    case '\n':
    case '\r':
    case '\v':
    case '\f':
    case '\x1c':
    case '\x1d':
    case '\x1e':
    case '\x1f':
      return true;
    default:
      return false;
  }
}

// The length of the longest prefix of [text, text + size) that str.rstrip()
// would keep.
std::size_t rstripped_size(const char* text, std::size_t size) noexcept;

// The first word of a record's title, which is the key both `SeqIO.index` and
// this package key a record by: `text.strip().split(None, 1)[0]`, so leading
// whitespace is skipped and a run of whitespace ends the word.  That is not
// the same thing as the title, and ">  rec1 desc" is the case that shows it.
//
// Both kernels call this one function rather than each writing the three lines
// out, because the two would otherwise be free to disagree about where a key
// ends -- and a key is what every lookup is addressed by.  The span form is
// what the FASTQ index uses, since it records a key as a place inside a title
// rather than as a copy of it.
//
// One deliberate imprecision, recorded rather than hidden: Biopython runs the
// FASTA index over a text stream and the FASTQ index over a binary one, so its
// two key rules differ on the four ASCII bytes 0x1c-0x1f, which Python strips
// in `str` and keeps in `bytes`.  This function is the `str` rule, because the
// *parser* is this package's contract for both formats; see the header of
// fastq.hpp for why the index follows the parser rather than the reference's
// own index.
std::pair<std::size_t, std::size_t> first_word_span(const char* text,
                                                    std::size_t size) noexcept;

inline std::string first_word(const char* text, std::size_t size) {
  const std::pair<std::size_t, std::size_t> span = first_word_span(text, size);
  return std::string(text + span.first, span.second);
}

inline std::string_view rstrip_view(const char* text, std::size_t size) noexcept {
  return std::string_view(text, rstripped_size(text, size));
}

// One line, as the text stream would hand it over.
struct Line {
  const char* text = nullptr;  // content; terminator excluded
  std::size_t len = 0;
  std::size_t next = 0;  // offset at which the following line starts
  bool exists = false;

  // What Python's line[0] would be.  An empty line still has a terminator, and
  // that terminator is what Python's indexing sees -- which is why a blank line
  // inside a FASTQ record is content and a blank line before a FASTA record is
  // skipped, without either parser writing the case down.
  char first() const noexcept { return len > 0 ? text[0] : '\n'; }
};

// Reads the line starting at `at`.  With `at >= size` the result does not
// exist, which is how a caller tells "no more lines" from "an empty last line".
//
// Both searches are bounded by the end of this line, not by the end of the
// buffer: a file with no carriage return in it must not be scanned twice just
// to discover that.
Line read_line(const char* data, std::size_t size, std::size_t at) noexcept;

// The keyword a flat-file header line introduces: the first token of the line's
// first `width` columns.  A continuation line is indented to exactly that
// column and so has no keyword of its own, which is what makes this the rule
// that says which lines belong to the keyword above them.  Comparing whole
// tokens rather than byte prefixes is what keeps "ID" from matching "IDE" and
// "SQ" from matching "SQX".
//
// Both header readers need this and neither may have its own copy: the record
// scan walks the same lines the annotations reader does, and two answers to
// "which keyword is this" is two answers to "which block is this".
inline std::string_view keyword_at(const Line& line, std::size_t width) noexcept {
  const std::size_t upto = std::min(width, line.len);
  std::size_t end = upto;
  while (end > 0 && is_python_space(line.text[end - 1])) --end;
  std::size_t begin = 0;
  while (begin < end && is_python_space(line.text[begin])) ++begin;
  return std::string_view(line.text + begin, end - begin);
}

// Python's `str.strip()`.  Both flat-file readers take values apart with it, and
// the rule they take them apart by -- which whitespace counts, and that it is
// stripped from both ends and not one -- is the same rule the line reader
// already owns, so it lives here rather than in either reader.
inline std::string_view strip_view(std::string_view text) noexcept {
  std::size_t begin = 0;
  std::size_t end = text.size();
  while (begin < end && is_python_space(text[begin])) ++begin;
  while (end > begin && is_python_space(text[end - 1])) --end;
  return text.substr(begin, end - begin);
}

// A keyworded line's value: everything past the keyword column, right-stripped.
// Left-stripping is deliberately *not* done here -- a value's leading spaces are
// meaningful to the one caller that joins continuation lines, and the callers
// that want Python's `.strip()` say so.
inline std::string_view value_at(const Line& line, std::size_t width) noexcept {
  if (line.len <= width) return std::string_view();
  return rstrip_view(line.text + width, line.len - width);
}

inline void append_rstripped(std::string& out, const Line& line) {
  out.append(line.text, rstripped_size(line.text, line.len));
}

}  // namespace biofasting
