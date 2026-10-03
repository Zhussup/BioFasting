#include "line.hpp"

#include <cstring>

namespace biofasting {

std::size_t rstripped_size(const char* text, std::size_t size) noexcept {
  while (size > 0 && is_python_space(text[size - 1])) --size;
  return size;
}

std::pair<std::size_t, std::size_t> first_word_span(const char* text,
                                                    std::size_t size) noexcept {
  const char* begin = text;
  const char* end = text + size;
  while (begin < end && is_python_space(*begin)) ++begin;
  while (end > begin && is_python_space(end[-1])) --end;
  const char* word = begin;
  while (word < end && !is_python_space(*word)) ++word;
  return {static_cast<std::size_t>(begin - text),
          static_cast<std::size_t>(word - begin)};
}

Line read_line(const char* data, std::size_t size, std::size_t at) noexcept {
  Line line;
  if (at >= size) return line;
  line.exists = true;
  line.text = data + at;
  const char* const end = data + size;
  const std::size_t remaining = size - at;
  const char* const lf =
      static_cast<const char*>(std::memchr(line.text, '\n', remaining));
  const char* const stop = lf ? lf : end;
  const char* const cr =
      static_cast<const char*>(std::memchr(line.text, '\r', stop - line.text));
  if (cr != nullptr) {
    line.len = static_cast<std::size_t>(cr - line.text);
    line.next = static_cast<std::size_t>(cr - data) + 1;
    if (cr + 1 < end && cr[1] == '\n') line.next += 1;
  } else if (lf != nullptr) {
    line.len = static_cast<std::size_t>(lf - line.text);
    line.next = static_cast<std::size_t>(lf - data) + 1;
  } else {
    line.len = remaining;  // last line, unterminated
    line.next = size;
  }
  return line;
}

}  // namespace biofasting
