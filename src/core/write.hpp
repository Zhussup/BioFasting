// The writers: the text of the sequence-only formats, built in one buffer.
//
// A reader turns bytes into records; a writer turns records into bytes, and the
// reference does it one small piece at a time -- eighteen `handle.write` calls
// for a 1 kb FASTA record, nine for a QUAL record, one per record for FASTQ.
// The kernels here build the whole record's text in a `std::string` and the
// caller writes it once, which is why they take the output by reference rather
// than returning a piece: a kernel that returned a fragment would leave the
// concatenation in Python, and the concatenation is half of what it removes.
//
// Every rule reproduced here is the reference's, quoted where it is surprising.
// The wrap widths, the newline that an empty FASTA sequence does *not* get, and
// the QUAL line cut are all behaviours a reader would guess differently.

#ifndef BIOFASTING_WRITE_HPP
#define BIOFASTING_WRITE_HPP

#include <cstddef>
#include <string>
#include <string_view>

namespace biofasting {

// One FASTA record: `>title`, then the sequence cut into lines of `wrap`.
//
// `wrap` of 0 means no wrapping and is not the same program: the reference
// takes a different branch, writes the sequence and a newline even when it is
// empty, where a wrapped record of no length writes no base line at all.
void append_fasta(std::string& out, std::string_view title,
                  std::string_view sequence, std::size_t wrap);

// One Sanger FASTQ record, on the four lines the format has.  `quality` is the
// ASCII quality string as our readers produce it -- phred+33 -- and is copied
// rather than re-encoded, because a kernel has no business decoding what the
// caller already holds.
void append_fastq(std::string& out, std::string_view title,
                  std::string_view sequence, std::string_view quality);

// One QUAL record: `>title`, then the scores as decimal, packed into lines.
//
// Returns false where the reference's own cut cannot be made -- a `wrap`-wide
// window with no space in it -- which is the shape that makes Biopython 1.88
// loop forever rather than write anything.  Declining is the one honest answer
// to a rule whose other outcome is a hang.
bool append_qual(std::string& out, std::string_view title,
                 std::string_view quality, std::size_t wrap);

// PHRED scores as the ASCII string a FASTQ record carries: the reference's
// `min(126, score + SANGER_SCORE_OFFSET)`, truncation at 93 included.  This is
// the encoding a `SeqRecord`'s `phred_quality` needs before any of the writers
// above will take it, and it is where the reference spends a dict lookup a
// base.
void phred_to_sanger(const unsigned char* scores, std::size_t count,
                     std::string& out);

}  // namespace biofasting

#endif  // BIOFASTING_WRITE_HPP
