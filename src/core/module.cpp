// The nanobind module.  Deliberately thin: it declares the names and hands off
// to build_info.cpp and cpu_features.cpp, so that adding a kernel in step 1.2
// does not mean editing the place where the ABI is defined.

#include <cstddef>
#include <string>

#include <cstdint>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "biofasting/build_config.hpp"
#include "build_info.hpp"
#include "cpu_features.hpp"
#include "fasta.hpp"
#include "fastq.hpp"
#include "inflate.hpp"
#include "protparam.hpp"
#include "seqops.hpp"
#include "sequtils.hpp"

namespace nb = nanobind;

namespace {

// Holds the buffer protocol view for as long as the scanner needs the bytes.
// This is what lets a whole 369 MB file be handed over as an `mmap` without a
// copy: the view keeps the exporter alive, the exporter keeps the pages
// mapped, and the kernel only ever sees a pointer and a length.
class BufferRef {
 public:
  BufferRef() = default;
  explicit BufferRef(nb::handle source) {
    if (PyObject_GetBuffer(source.ptr(), &view_, PyBUF_SIMPLE) != 0) {
      throw nb::python_error();
    }
    held_ = true;
  }
  ~BufferRef() {
    if (held_) PyBuffer_Release(&view_);
  }
  BufferRef(const BufferRef&) = delete;
  BufferRef& operator=(const BufferRef&) = delete;
  // Movable so that a scanner can be returned from a factory: the buffer
  // reference is transferred, not duplicated, or the two copies would each
  // release it.
  BufferRef(BufferRef&& other) noexcept : view_(other.view_), held_(other.held_) {
    other.held_ = false;
  }
  BufferRef& operator=(BufferRef&&) = delete;

  const char* data() const noexcept {
    return static_cast<const char*>(view_.buf);
  }
  std::size_t size() const noexcept {
    return static_cast<std::size_t>(view_.len);
  }

 private:
  Py_buffer view_{};
  bool held_ = false;
};

// The Python face of FastqScanner.  Member order matters: the buffer a scanner
// reads from is declared before the scanner that points into it, and a scanner
// either points at the caller's buffer or owns an inflated one -- never both.
class PyFastqScanner {
 public:
  explicit PyFastqScanner(nb::handle source)
      : buffer_(source), scanner_(buffer_.data(), buffer_.size()) {}

  // `source` is a whole gzip file: inflate it and scan the result.  The
  // compressed buffer is only needed until the inflation finishes.
  static PyFastqScanner from_gzip(nb::handle source) {
    BufferRef compressed(source);
    biofasting::GzipBuffer inflated;
    std::string error;
    if (!biofasting::gunzip(compressed.data(), compressed.size(), inflated,
                            error)) {
      throw nb::value_error(error.c_str());
    }
    return PyFastqScanner(std::move(inflated));
  }

  nb::tuple next() {
    biofasting::FastqRecord record;
    switch (scanner_.next(record)) {
      case biofasting::FastqStatus::ok:
        break;
      case biofasting::FastqStatus::end:
        throw nb::stop_iteration();
      case biofasting::FastqStatus::error:
        throw nb::value_error(located_error().c_str());
    }
    return nb::make_tuple(
        nb::bytes(record.title.data(), record.title.size()),
        nb::bytes(record.sequence.data(), record.sequence.size()),
        nb::bytes(record.quality.data(), record.quality.size()));
  }

  std::size_t record_index() const noexcept { return scanner_.record_index(); }
  std::size_t position() const noexcept { return scanner_.position(); }

  // The zero-copy views: two strided two-dimensional arrays over the bytes the
  // scanner is reading, one for the sequences and one for the qualities, plus
  // the header length that the arrays' own strides cannot express.
  //
  // `owner` is the Python object the arrays keep alive.  It has to be: the
  // bytes belong either to the caller's buffer or to a buffer this scanner
  // inflated, and in both cases the scanner is what holds them -- an array that
  // outlived it would be a read of freed memory, and numpy would have no way to
  // know.
  nb::tuple grid(nb::handle owner) const {
    const biofasting::FastqGrid geometry =
        biofasting::scan_grid(scanner_.data(), scanner_.size());
    if (!geometry.ok()) throw nb::value_error(geometry.reason.c_str());

    // `const` in the scalar type is what makes these read-only, and read-only
    // is not a nicety: the bytes may be a mapping of the user's file, so a
    // stray write through the array would be a write to their data.
    using ReadOnlyBytes =
        nb::ndarray<nb::numpy, const std::uint8_t, nb::shape<-1, -1>>;

    const auto row_stride = static_cast<std::int64_t>(geometry.stride);
    const char* bytes = scanner_.data();
    ReadOnlyBytes sequences(
        bytes + geometry.sequence_offset,
        {geometry.count, geometry.length}, owner, {row_stride, 1});
    ReadOnlyBytes qualities(
        bytes + geometry.quality_offset,
        {geometry.count, geometry.length}, owner, {row_stride, 1});
    return nb::make_tuple(sequences, qualities, geometry.title_length);
  }

 private:
  // Only used by from_gzip, where the bytes being scanned are the ones this
  // object inflated rather than the ones the caller passed in.
  explicit PyFastqScanner(biofasting::GzipBuffer&& inflated)
      : inflated_(std::move(inflated)),
        scanner_(inflated_.data(), inflated_.size()) {}

  // Biopython's own wording, because the differential tests read it, plus the
  // location -- which Biopython cannot give and a malformed file of a million
  // records makes valuable.
  std::string located_error() const {
    return scanner_.error_message() + " (record " +
           std::to_string(scanner_.record_index()) + ", byte offset " +
           std::to_string(scanner_.error_offset()) + ")";
  }

  BufferRef buffer_;
  biofasting::GzipBuffer inflated_;
  biofasting::FastqScanner scanner_;
};

// The Python face of FastaIndex.  Same buffer-lifetime rule as the FASTQ
// scanner: the index points into the caller's buffer, so the reference is
// declared first and held for as long as the index lives.
class PyFastaIndex {
 public:
  explicit PyFastaIndex(nb::handle source)
      : buffer_(source), index_(buffer_.data(), buffer_.size()) {
    if (index_.failed()) throw nb::value_error(index_.error_message().c_str());
  }

  std::size_t size() const noexcept { return index_.size(); }

  bool has(const std::string& name) const { return index_.find(name) != nullptr; }

  // The sequence as bytes, not as a record: converting to a SeqRecord is the
  // interop shim's job, and doing it here would put a Biopython dependency in
  // the core.
  nb::bytes get(const std::string& name) const {
    const biofasting::FastaEntry& entry = require(name);
    const std::string sequence = index_.sequence(entry);
    return nb::bytes(sequence.data(), sequence.size());
  }

  nb::bytes slice(const std::string& name, std::size_t start,
                  std::size_t end) const {
    const biofasting::FastaEntry& entry = require(name);
    const std::string window = index_.sequence_slice(entry, start, end);
    return nb::bytes(window.data(), window.size());
  }

  std::size_t sequence_length(const std::string& name) const {
    return require(name).length;
  }

  // Whether reading this record is a bulk copy (true) or a line-by-line
  // assembly (false).  Not an internal detail: it is the difference between a
  // slice costing microseconds and costing a full materialisation, and a test
  // asserting it is what stops the fast path from silently disappearing.
  bool strided(const std::string& name) const { return require(name).strided; }

  std::string title(const std::string& name) const {
    return require(name).title;
  }

  // The .fai row for one record, in the file's own column order, so the test
  // that checks our index against the pyfaidx-produced .fai is comparing the
  // same five things pyfaidx wrote.
  nb::tuple index_entry(const std::string& name) const {
    const biofasting::FastaEntry& entry = require(name);
    return nb::make_tuple(entry.name, entry.length, entry.offset,
                          entry.line_bases, entry.line_width);
  }

  // The record as a strided two-dimensional view of the file's own bytes:
  // (lines, line_bases) at a row stride of line_width.  This is what makes a
  // whole chromosome available to numpy for nothing -- 40 Mbp is one mapping
  // and one shape, not a 40 MB copy.
  //
  // The array covers whole lines only, and says so by being returned next to
  // the record's real length: a chromosome of 40,000,000 bases at 60 per line
  // ends with a 40-base line, and no constant stride can express that last
  // row.  A caller that assumed otherwise would silently analyse 39,999,960
  // bases, which is why the length travels with the array rather than being
  // left to be looked up.
  nb::tuple grid(nb::handle owner, const std::string& name) const {
    const biofasting::FastaEntry& entry = require(name);
    if (!entry.strided) {
      throw nb::value_error(
          ("record '" + name +
           "' is not a regular grid of lines, so it has no strided view")
              .c_str());
    }
    const std::size_t rows =
        entry.line_bases == 0 ? 0 : entry.length / entry.line_bases;
    if (rows == 0) {
      throw nb::value_error(
          ("record '" + name + "' is shorter than one line, so it has no grid")
              .c_str());
    }
    using ReadOnlyBytes =
        nb::ndarray<nb::numpy, const std::uint8_t, nb::shape<-1, -1>>;
    ReadOnlyBytes bases(
        buffer_.data() + entry.offset, {rows, entry.line_bases}, owner,
        {static_cast<std::int64_t>(entry.line_width), 1});
    return nb::make_tuple(bases, entry.length);
  }

  nb::list keys() const {
    nb::list names;
    for (const biofasting::FastaEntry& entry : index_.entries()) {
      names.append(nb::str(entry.name.c_str(), entry.name.size()));
    }
    return names;
  }

 private:
  const biofasting::FastaEntry& require(const std::string& name) const {
    const biofasting::FastaEntry* entry = index_.find(name);
    if (entry == nullptr) throw nb::key_error(name.c_str());
    return *entry;
  }

  BufferRef buffer_;
  biofasting::FastaIndex index_;
};

// The Python face of FastqIndex.  Same buffer-lifetime rule as FastaIndex --
// and here it is load-bearing twice over, because the index does not merely
// point into the buffer: it *keys* itself on views of the buffer's bytes, so a
// buffer that moved would take every key in the map with it.
class PyFastqIndex {
 public:
  explicit PyFastqIndex(nb::handle source)
      : buffer_(source), index_(buffer_.data(), buffer_.size()) {
    if (index_.failed()) throw nb::value_error(index_.error_message().c_str());
  }

  // `source` is a whole gzip file -- any gzip file, BGZF included, since BGZF
  // is gzip with extra fields -- and the index is built over the inflated
  // bytes.  Returns a pointer rather than a value because it has to: the index
  // holds a `std::pmr::monotonic_buffer_resource` and a hash table keyed on
  // views of the inflated buffer, none of which can be relocated, so there is
  // no move constructor to return by.  The caller owns the result.
  //
  // Inflating the whole file rather than reading it block by block is the same
  // decision `from_gzip` on the scanner makes and for the same reason, and it
  // is what makes a fetch from a compressed file cost the same as a fetch from
  // an uncompressed one.  The alternative -- seek into the compressed stream,
  // inflate the enclosing block, throw the cache away on the next seek -- is
  // what `SeqIO.index` does on BGZF, and it costs 274 µs an access against this
  // index's 0.32 µs.
  static PyFastqIndex* from_gzip(nb::handle source) {
    BufferRef compressed(source);
    biofasting::GzipBuffer inflated;
    std::string error;
    if (!biofasting::gunzip(compressed.data(), compressed.size(), inflated,
                            error)) {
      throw nb::value_error(error.c_str());
    }
    return new PyFastqIndex(std::move(inflated));
  }

  std::size_t size() const noexcept { return index_.size(); }

  bool has(const std::string& name) const { return index_.find(name) != nullptr; }

  // Bytes, not a SeqRecord: building one is the interop shim's job, and doing
  // it here would put a Biopython dependency in the core.  The reference's
  // index returns a record and pays for it on every fetch; that difference is
  // the whole point of this class.
  nb::bytes get(const std::string& name) const {
    const std::string sequence = index_.sequence(require(name));
    return nb::bytes(sequence.data(), sequence.size());
  }

  nb::bytes quality(const std::string& name) const {
    const std::string text = index_.quality(require(name));
    return nb::bytes(text.data(), text.size());
  }

  nb::bytes slice(const std::string& name, std::size_t start,
                  std::size_t end) const {
    const std::string window = index_.sequence_slice(require(name), start, end);
    return nb::bytes(window.data(), window.size());
  }

  std::size_t sequence_length(const std::string& name) const {
    return require(name).length;
  }

  // The same number as sequence_length(), and not a second measurement: the
  // parser refuses a record whose quality is not exactly as long as its
  // sequence, so there is one length here and not two that could drift.
  std::size_t quality_length(const std::string& name) const {
    return require(name).length;
  }

  std::string title(const std::string& name) const {
    return index_.title(require(name));
  }

  // Whether reading this record is a memcpy (true) or a re-parse of its own
  // span (false).  Not an internal detail: it is the difference between a fetch
  // costing microseconds and costing what the reference charges for every
  // fetch, and a test asserting it is what stops the fast path from quietly
  // disappearing.
  bool plain(const std::string& name) const { return require(name).plain; }

  // Where the record is, in the same shape FastaIndex exposes, so that the two
  // indexes can be checked against each other and against the file.
  nb::tuple index_entry(const std::string& name) const {
    const biofasting::FastqEntry& entry = require(name);
    return nb::make_tuple(name, entry.length, entry.offset, entry.end,
                          entry.plain);
  }

  nb::list keys() const {
    nb::list names;
    for (const biofasting::FastqEntry& entry : index_.entries()) {
      const std::string_view key = index_.key_of(entry);
      names.append(nb::str(key.data(), key.size()));
    }
    return names;
  }

 private:
  // A compressed source is inflated into `inflated_`, whose bytes become the
  // ones indexed; an uncompressed one is read through `buffer_`.  Exactly one
  // of the two is ever non-empty, and both are declared before `index_`,
  // because `index_` points into whichever one it was given -- and, for the
  // plain case, keys itself on views of it.
  explicit PyFastqIndex(biofasting::GzipBuffer&& inflated)
      : inflated_(std::move(inflated)),
        index_(inflated_.data(), inflated_.size()) {
    if (index_.failed()) throw nb::value_error(index_.error_message().c_str());
  }

  const biofasting::FastqEntry& require(const std::string& name) const {
    const biofasting::FastqEntry* entry = index_.find(name);
    if (entry == nullptr) throw nb::key_error(name.c_str());
    return *entry;
  }

  BufferRef buffer_;
  biofasting::GzipBuffer inflated_;
  biofasting::FastqIndex index_;
};

}  // namespace

NB_MODULE(_core, m) {
  using namespace nanobind::literals;  // "source"_a


  m.doc() =
      "C++ core of BioFasting.\n\n"
      "Private: import `biofasting` instead.  Nothing in this module is a "
      "stable interface, and in particular the extension module name may "
      "change without a deprecation cycle while the project is pre-alpha.";

  m.attr("__version__") = BIOFASTING_BUILD_VERSION;

  m.def("build_info", &biofasting::build_info,
        "How this extension was compiled, as plain strings.  The runtime half "
        "of the same answer is added by biofasting.build_info().");

  m.def("cpu_levels", &biofasting::cpu_levels,
        "The ISA ladder this build defines, ascending, \"baseline\" first.");

  m.def("supported_cpu_levels", &biofasting::supported_cpu_levels,
        "The subset of cpu_levels() the running CPU supports.  Detected once, "
        "on first call.");

  m.def("best_cpu_level", &biofasting::best_cpu_level,
        "The highest supported level: the rung a kernel dispatches to.");

  nb::class_<PyFastqScanner>(m, "FastqScanner",
                             "Scans FASTQ records out of any object exposing "
                             "the buffer protocol.\n\n"
                             "Iterating yields (title, sequence, quality) as "
                             "bytes, with title stripped of its '@'.\n"
                             "Semantics are Bio.SeqIO.QualityIO."
                             "FastqGeneralIterator's, including its "
                             "interpretation of '@' inside a quality line; a "
                             "malformed record raises ValueError carrying the "
                             "original message plus the byte offset and record "
                             "number where it was found.\n"
                             "Use biofasting.fastq.open_fastq() rather than "
                             "building this directly: it is the piece that "
                             "keeps the buffer alive.")
      .def(nb::init<nb::handle>(), "source"_a,
           "Scan `source`, which must expose the buffer protocol.")
      .def_static("from_gzip", &PyFastqScanner::from_gzip, "source"_a,
                  "Inflate a whole gzip file from a buffer and scan the "
                  "result.\n\n"
                  "Every concatenated member is decompressed.  The compressed "
                  "bytes are not needed afterwards, so only the inflated "
                  "buffer is retained.  Raises ValueError if the bytes are not "
                  "a valid gzip stream.")
      // Returning the Python object rather than a C++ reference is not a
      // style choice.  A lambda that returns `PyFastqScanner&` makes nanobind
      // mint a second wrapper for an instance that already has one, and under
      // Release its assertion for that case is compacted into "encountered an
      // unrecoverable error condition" -- found by `list(scanner)` aborting the
      // interpreter, not by reading documentation.
      .def("__iter__", [](nb::object self) { return self; })
      .def("__next__", &PyFastqScanner::next)
      .def_prop_ro("record_index", &PyFastqScanner::record_index,
                   "Number of records returned so far.")
      .def_prop_ro("position", &PyFastqScanner::position,
                   "Bytes consumed so far: the offset just past the last "
                   "record returned.")
      // A lambda rather than a member pointer because the arrays have to keep
      // *this Python object* alive, and a member function is handed only
      // `this`.  Without the distinction the arrays would outlive the mapping
      // they point into.
      .def(
          "grid",
          [](nb::object self) {
            return nb::cast<PyFastqScanner&>(self).grid(self);
          },
          "The whole file as two zero-copy 2-D arrays, "
          "(sequences, qualities, title_length).\n\n"
          "Rows are records, columns are positions within a record, and dtype "
          "is uint8 -- the bytes of the file itself, strided, with no copy "
          "made.  Both arrays hold this scanner alive, so the mapping (or the "
          "inflated buffer) cannot be released while they are readable, and "
          "both are read-only.\n"
          "Raises ValueError if the records are not all the same shape, "
          "naming the first record that differs: a grid is only a grid if "
          "every record is on it, and the alternative is an array that is "
          "silently wrong on the one ragged file in a thousand.");

  nb::class_<PyFastaIndex>(m, "FastaIndex",
                           "Indexes the records of a FASTA file held in any "
                           "object exposing the buffer protocol.\n\n"
                           "Keys, record ids and sequence bytes follow "
                           "Bio.SeqIO: the key is the title's first word, a "
                           "duplicate key raises ValueError, and a sequence is "
                           "the record's lines right-stripped, joined, and "
                           "stripped of spaces.  Reading a record afterwards "
                           "does not re-parse the file.\n"
                           "Use biofasting.fasta.open_fasta() rather than "
                           "building this directly: it is the piece that keeps "
                           "the buffer alive.")
      .def(nb::init<nb::handle>(), "source"_a,
           "Index `source`, which must expose the buffer protocol.")
      .def("__len__", &PyFastaIndex::size, "Number of records.")
      .def("__contains__", &PyFastaIndex::has, "name"_a,
           "Whether a record with this key exists.")
      .def("__getitem__", &PyFastaIndex::get, "name"_a,
           "The record's sequence, as bytes.  Raises KeyError if there is no "
           "such record.")
      .def("__iter__",
           [](const PyFastaIndex& self) { return nb::iter(self.keys()); },
           "Iterate the record keys, in file order.")
      .def("keys", &PyFastaIndex::keys,
           "The record keys, in file order, as a list.")
      .def("title", &PyFastaIndex::title, "name"_a,
           "The whole title line after '>', right-stripped.")
      .def("sequence_length", &PyFastaIndex::sequence_length, "name"_a,
           "The sequence length, from the index rather than from the bytes.")
      .def("strided", &PyFastaIndex::strided, "name"_a,
           "Whether the record's lines form the regular grid that makes a "
           "bulk copy -- and so a cheap slice -- possible.")
      .def("sequence_slice", &PyFastaIndex::slice, "name"_a, "start"_a, "end"_a,
           "The [start, end) window of the sequence.  Both bounds are clamped "
           "to the record.")
      .def("index_entry", &PyFastaIndex::index_entry, "name"_a,
           "The record's .fai row: (name, length, offset, line_bases, "
           "line_width).")
      .def(
          "grid",
          [](nb::object self, const std::string& name) {
            return nb::cast<PyFastaIndex&>(self).grid(self, name);
          },
          "name"_a,
          "The record as (bases, length), where `bases` is a zero-copy "
          "read-only (lines, line_bases) uint8 view over the file itself.\n\n"
          "`bases` covers whole lines only: a record whose length is not a "
          "multiple of its line width has a tail that no constant stride can "
          "express, and `length` is the record's true length so that the "
          "difference is visible rather than assumed.  The array holds the "
          "index alive, and so the mapping with it.\n"
          "Raises ValueError if the record's lines are not a regular grid.");

  nb::class_<PyFastqIndex>(m, "FastqIndex",
                           "Indexes the records of a FASTQ file held in any "
                           "object exposing the buffer protocol, so that "
                           "fetching one by name does not parse the file "
                           "again.\n\n"
                           "Keys follow Bio.SeqIO.index: the header's first "
                           "word, in file order, with a duplicate refused.  A "
                           "fetch is a copy of the record's bytes, not a "
                           "SeqRecord -- converting is the interop shim's job, "
                           "and doing it here would put a Biopython dependency "
                           "in the core.\n"
                           "Which files this accepts is FastqGeneralIterator's "
                           "rule, not SeqIO.index's: Biopython's index runs a "
                           "second FASTQ grammar that disagrees with its own "
                           "parser about blank lines, and this package has one "
                           "grammar.  See the header of src/core/fastq.hpp.\n"
                           "Use biofasting.fastq.open_fastq_index() rather than "
                           "building this directly: it is the piece that keeps "
                           "the buffer alive -- which here matters twice, "
                           "because the lookup table is keyed on views of the "
                           "buffer's own bytes.")
      .def(nb::init<nb::handle>(), "source"_a,
           "Index `source`, which must expose the buffer protocol.")
      .def_static("from_gzip", &PyFastqIndex::from_gzip, "source"_a,
                  nb::rv_policy::take_ownership,
                  "Index a whole gzip file, `source`, by inflating it and "
                  "indexing the inflated bytes.  Any gzip stream will do, "
                  "BGZF included -- BGZF is gzip with extra fields -- so a "
                  "file `SeqIO.index` refuses for being plain gzip is "
                  "readable here.  The owner keeps the inflated bytes alive; "
                  "that is what the returned object holds.")
      .def("__len__", &PyFastqIndex::size, "Number of records.")
      .def("__contains__", &PyFastqIndex::has, "name"_a,
           "Whether a record with this key exists.")
      .def("__getitem__", &PyFastqIndex::get, "name"_a,
           "The record's sequence, as bytes.  Raises KeyError if there is no "
           "such record.")
      .def("__iter__",
           [](const PyFastqIndex& self) { return nb::iter(self.keys()); },
           "Iterate the record keys, in file order.")
      .def("keys", &PyFastqIndex::keys,
           "The record keys, in file order, as a list.")
      .def("title", &PyFastqIndex::title, "name"_a,
           "The whole header line after '@', right-stripped.")
      .def("sequence", &PyFastqIndex::get, "name"_a,
           "The record's sequence, as bytes.  The same call as index[name].")
      .def("sequence_length", &PyFastqIndex::sequence_length, "name"_a,
           "The sequence length, from the index rather than from the bytes.")
      .def("quality_length", &PyFastqIndex::quality_length, "name"_a,
           "The quality length.  Equal to sequence_length() by construction: "
           "the parser refuses a record whose two lengths differ.")
      .def("quality", &PyFastqIndex::quality, "name"_a,
           "The record's quality string, as bytes.")
      .def("sequence_slice", &PyFastqIndex::slice, "name"_a, "start"_a, "end"_a,
           "The [start, end) window of the sequence.  Both bounds are clamped "
           "to the record.")
      .def("plain", &PyFastqIndex::plain, "name"_a,
           "Whether reading this record is a memcpy (true) or a re-parse of "
           "its own span (false).")
      .def("index_entry", &PyFastqIndex::index_entry, "name"_a,
           "The record's place in the file: (name, length, offset, end, "
           "plain).");

  m.def("seqops_level", &biofasting::seqops_level,
        "The rung the sequence kernels dispatch to: \"baseline\" or \"avx2\".\n\n"
        "Not the same question as best_cpu_level(): a CPU can support AVX-512 "
        "while no kernel here uses it, and this is the answer a benchmark has "
        "to report, because it is the code that ran.");

  m.def(
      "reverse_complement",
      [](nb::handle source) {
        const BufferRef buffer(source);
        const std::string result =
            biofasting::reverse_complement(buffer.data(), buffer.size());
        return nb::bytes(result.data(), result.size());
      },
      "sequence"_a,
      "Bio.Seq.reverse_complement, byte for byte: the full IUPAC set in both "
      "cases, and every other byte left alone.");

  m.def(
      "gc_fraction",
      [](nb::handle source) {
        const BufferRef buffer(source);
        return biofasting::gc_fraction(buffer.data(), buffer.size());
      },
      "sequence"_a,
      "Bio.SeqUtils.gc_fraction with its default ambiguous=\"remove\": G, C "
      "and S count as GC, A, T, W and U make up the rest of the denominator, "
      "and everything else is removed from both.\n\n"
      "Returns -1.0 when nothing was counted at all, which is a sentinel and "
      "not a fraction: the reference answers the *integer* 0 there and only the "
      "caller, which has an int, can say so.  A real fraction is never "
      "negative.");

  m.def(
      "count_kmers",
      [](nb::handle source, unsigned k) {
        const BufferRef buffer(source);
        try {
          const auto counts =
              biofasting::count_kmers(buffer.data(), buffer.size(), k);
          nb::list out;
          for (const auto& [kmer, count] : counts) {
            out.append(nb::make_tuple(nb::bytes(kmer.data(), kmer.size()),
                                      nb::int_(count)));
          }
          return out;
        } catch (const std::invalid_argument& error) {
          throw nb::value_error(error.what());
        }
      },
      "sequence"_a, "k"_a,
      "Overlapping k-mer counts over ACGT, ascending by k-mer, as a list of "
      "(bytes, int).\n\n"
      "Case-insensitive, and a window containing anything that is not ACGT is "
      "skipped rather than split around, so the counts are exactly "
      "Seq.count_overlap's per k-mer.");

  m.def(
      "translate",
      [](nb::handle source, nb::bytes code, nb::bytes amino,
         const std::string& stop_symbol, const std::string& possible_stop,
         bool to_stop, bool stop_is_error) {
        const BufferRef buffer(source);
        if (code.size() != 256) {
          throw nb::value_error("code must be 256 bytes: one base code per byte value");
        }
        if (amino.size() != 64) {
          throw nb::value_error("amino must be 64 bytes: one entry per codon");
        }
        const biofasting::TranslateResult result = biofasting::translate(
            buffer.data(), buffer.size(), code.c_str(), amino.c_str(),
            stop_symbol, possible_stop, to_stop, stop_is_error);
        //   The protein comes back as `str`, not as `bytes`, because it is not
        //   necessarily ASCII: `stop_symbol` is an arbitrary Python string, and
        //   the reference concatenates it without a look at it.  The only
        //   non-ASCII bytes the kernel can emit are whole characters copied out
        //   of `stop_symbol` or `possible_stop`, so UTF-8 is the encoding both
        //   sides already agree on.
        return nb::make_tuple(static_cast<int>(result.outcome), result.codon,
                              nb::str(result.protein.data(), result.protein.size()));
      },
      "sequence"_a, "code"_a, "amino"_a, "stop_symbol"_a, "possible_stop"_a,
      "to_stop"_a, "stop_is_error"_a,
      "One reading frame of a genetic code, described by two tables.\n\n"
      "`code` is 256 bytes: the 2-bit base value of every byte value (A=0, "
      "C=1, G=2, T=3), or 0x80 for a byte this genetic code does not accept.  "
      "`amino` is 64 bytes indexed by ``a << 4 | b << 2 | c``: the amino-acid "
      "letter, or 1 for a stop codon, or 2 for a codon the reference answers "
      "with `possible_stop` (TAN and friends), or 0 for a codon it refuses.\n\n"
      "Returns (outcome, codon, protein) where outcome is 0 for a finished "
      "translation, 1 when `codon` holds something the tables cannot classify "
      "-- in which case `protein` means nothing and the caller redoes the "
      "work -- and 2 when cds=True met an in-frame stop at `codon`.\n\n"
      "`protein` is a str: the reference never checks `stop_symbol`, so it may "
      "be any string at all and the kernel copies it through as it stands.");

  // -------------------------------------------------------------------------
  // Bio.SeqUtils: the functions that measure a sequence rather than change one.

  m.def(
      "gc123",
      [](nb::handle source) {
        const BufferRef buffer(source);
        const biofasting::Gc123Counts counts =
            biofasting::gc123(buffer.data(), buffer.size());
        // Flat and frame-major: frame f's A, T, G and C counts are at 4f .. 4f+3.
        // Twelve ints and not four percentages, because the reference divides
        // three times and raises on the fourth case -- see sequtils.hpp.
        std::vector<std::uint64_t> flat;
        flat.reserve(12);
        for (const auto& frame : counts.counts) {
          flat.insert(flat.end(), frame, frame + 4);
        }
        return flat;
      },
      "sequence"_a,
      "G+C per reading frame, counted and not divided.\n\n"
      "Returns twelve counts, frame-major, each frame's four bases in the order "
      "A, T, G, C, over both cases and over nothing else -- which is how the "
      "reference's per-codon `==` comparisons read, and why the one or two "
      "bases of a partial trailing codon are counted by nobody.");

  m.def(
      "gc_skew",
      [](nb::handle source, std::size_t window) {
        const BufferRef buffer(source);
        try {
          return biofasting::gc_skew(buffer.data(), buffer.size(), window);
        } catch (const std::invalid_argument& error) {
          throw nb::value_error(error.what());
        }
      },
      "sequence"_a, "window"_a,
      "GC skew (G-C)/(G+C) per non-overlapping window, as a list of floats.\n\n"
      "Windows with no G and no C are 0.0, which is the reference's "
      "`except ZeroDivisionError` arm; the last window is short.");

  m.def(
      "gc_counts",
      [](nb::handle source) {
        const BufferRef buffer(source);
        const biofasting::GcCounts counts =
            biofasting::gc_counts(buffer.data(), buffer.size());
        std::vector<std::uint64_t> flat;
        flat.reserve(12);
        flat.push_back(counts.gc);
        flat.push_back(counts.at);
        flat.insert(flat.end(), counts.ambiguous, counts.ambiguous + 10);
        return flat;
      },
      "sequence"_a,
      "gc_fraction's counts for its \"ignore\" and \"weighted\" modes.\n\n"
      "Twelve numbers: the \"CGScgs\" count, the \"ATWUatwu\" count, and then "
      "`count(x) + count(x.lower())` for x in \"BDHKMNRVXY\" in that order.  The "
      "default \"remove\" mode is answered by the vectorised gc_fraction "
      "kernel; this is the one pass that serves the other two.");

  m.def(
      "molecular_weight_mass",
      [](nb::handle source, nb::bytes weights, nb::bytes valid) {
        const BufferRef buffer(source);
        if (weights.size() != 256 * sizeof(double)) {
          throw nb::value_error("weights must be 256 doubles");
        }
        if (valid.size() != 256) {
          throw nb::value_error("valid must be 256 bytes");
        }
        double total = 0.0;
        unsigned char bad = 0;
        const bool ok = biofasting::molecular_weight_mass(
            buffer.data(), buffer.size(),
            reinterpret_cast<const double*>(weights.c_str()),
            reinterpret_cast<const unsigned char*>(valid.c_str()), &total,
            &bad);
        return nb::make_tuple(total, ok ? -1 : static_cast<int>(bad));
      },
      "sequence"_a, "weights"_a, "valid"_a,
      "The sum inside molecular_weight and gravy, as CPython's sum() does it.\n\n"
      "`weights` is 256 little-endian doubles and `valid` 256 flags: a sparse "
      "weight table handed over dense, so that the per-letter branch is not in "
      "the loop.  Returns (mass, first_bad), where first_bad is -1 on success "
      "and otherwise the first byte the table has no weight for -- the "
      "reference's KeyError, whose message only Python can spell.\n\n"
      "The flag is not a boolean.  0 is absent, 1 is a float item and 2 an "
      "integer one, because sum() compensates the first and adds the second "
      "plain: ten of the twenty-eight gravy scales carry integers, and four of "
      "those answer differently for it.\n\n"
      "The compensation is not decoration either.  CPython 3.12 and later sum "
      "floats with Neumaier compensation and the reference is written as "
      "`sum(...)`, so a plain running total reproduces it only about three "
      "times in four.");

  // Exposed so that a test can build a flag table by hand and say what each of
  // its entries is, rather than writing a bare 1 or 2 that nothing checks.
  m.attr("kAbsent") = biofasting::kAbsent;
  m.attr("kCompensated") = biofasting::kCompensated;
  m.attr("kPlain") = biofasting::kPlain;

  m.def(
      "gcg",
      [](nb::handle source) {
        const BufferRef buffer(source);
        std::uint32_t checksum = 0;
        if (!biofasting::gcg(buffer.data(), buffer.size(), &checksum)) {
          return nb::object(nb::none());
        }
        return nb::object(nb::int_(checksum));
      },
      "sequence"_a,
      "The GCG checksum, or None if the sequence is not ASCII.\n\n"
      "`None` is a decline, not a result: the reference upper-cases each "
      "character, and `\"ß\".upper()` is two characters, so a non-ASCII letter "
      "has to go through Python's Unicode tables.  The caller runs the "
      "reference's own loop in that case.");

  m.def(
      "crc64",
      [](nb::handle source) {
        const BufferRef buffer(source);
        std::uint32_t high = 0;
        std::uint32_t low = 0;
        biofasting::crc64(buffer.data(), buffer.size(), &high, &low);
        return nb::make_tuple(high, low);
      },
      "sequence"_a,
      "The crc64 halves, high first, as two unsigned 32-bit ints.\n\n"
      "No decline: the reference folds each character in with "
      "`ord(c) & 0xFF`, so a byte is already the whole of what it reads.");

  m.def(
      "crc64_table_h",
      []() { return biofasting::crc64_table_h(); },
      "The 256-entry high table the crc64 kernel uses, derived at compile "
      "time from the reference's own recurrence and exposed so that a test can "
      "compare all 256 entries against `CheckSum._table_h`.");

  // --- ProteinAnalysis: the three per-residue kernels --------------------- //
  //
  // Each returns `None` instead of a number when the reference would raise or
  // warn, and the Python layer then runs the reference's own loop.  That is the
  // same decline `gcg` makes, for the same reason: the message names a residue,
  // and the residue it names is the first one *read*, which is not the first one
  // in the sequence.

  m.def(
      "protein_scale_scores",
      [](nb::handle source, long long window, nb::bytes weights,
         nb::bytes table, nb::bytes valid, double sum_of_weights) {
        const BufferRef buffer(source);
        if (table.size() != 256 * sizeof(double)) {
          throw nb::value_error("table must be 256 doubles");
        }
        if (valid.size() != 256) {
          throw nb::value_error("valid must be 256 bytes");
        }
        if (window < 2) {
          return nb::object(nb::none());
        }
        const std::size_t half = static_cast<std::size_t>(window) / 2;
        if (weights.size() != half * sizeof(double)) {
          throw nb::value_error("weights must be window // 2 doubles");
        }
        std::vector<double> out;
        const bool ok = biofasting::protein_scale_scores(
            buffer.data(), buffer.size(), static_cast<std::size_t>(window),
            reinterpret_cast<const double*>(weights.c_str()),
            reinterpret_cast<const double*>(table.c_str()),
            reinterpret_cast<const unsigned char*>(valid.c_str()), sum_of_weights,
            out);
        if (!ok) {
          return nb::object(nb::none());
        }
        return nb::cast(out);
      },
      "sequence"_a, "window"_a, "weights"_a, "table"_a, "valid"_a,
      "sum_of_weights"_a,
      "One `ProteinAnalysis.protein_scale` score per window that fits, or None "
      "when the reference would warn or raise.\n\n"
      "`weights` is the half-window of weights the reference computed; `table` "
      "is 256 doubles and `valid` 256 flags, a caller's scale handed over dense, "
      "because `protein_scale` takes an arbitrary mapping and there is no fixed "
      "letter order to index it by.  A byte the scale has no entry for is the "
      "decline: the reference writes a warning to stderr and drops that window's "
      "term, and reproducing Python's stderr from a kernel is not a kernel's "
      "job.");

  m.def(
      "flexibility_scores",
      [](nb::handle source, nb::bytes flex, nb::bytes map) {
        const BufferRef buffer(source);
        if (flex.size() != 20 * sizeof(double)) {
          throw nb::value_error("flex must be 20 doubles");
        }
        if (map.size() != 256) {
          throw nb::value_error("map must be 256 bytes");
        }
        std::vector<double> out;
        const bool ok = biofasting::flexibility_scores(
            buffer.data(), buffer.size(),
            reinterpret_cast<const double*>(flex.c_str()),
            reinterpret_cast<const unsigned char*>(map.c_str()), out);
        if (!ok) {
          return nb::object(nb::none());
        }
        return nb::cast(out);
      },
      "sequence"_a, "flex"_a, "map"_a,
      "`ProteinAnalysis.flexibility`'s scores, or None for a residue `Flex` has "
      "no value for.\n\n"
      "Not a special case of `protein_scale_scores`: the pair is added before it "
      "is multiplied, the middle is read at index 5 of the nine-residue window "
      "rather than at 4, index 4 is never read at all, and there is one window "
      "fewer than fit.  All four are the reference's, and all four change the "
      "last bit of the answer.");

  m.def(
      "instability_index_sum",
      [](nb::handle source, nb::bytes table, nb::bytes map) {
        const BufferRef buffer(source);
        if (table.size() != 400 * sizeof(double)) {
          throw nb::value_error("table must be 400 doubles");
        }
        if (map.size() != 256) {
          throw nb::value_error("map must be 256 bytes");
        }
        double total = 0.0;
        const bool ok = biofasting::instability_index_sum(
            buffer.data(), buffer.size(),
            reinterpret_cast<const double*>(table.c_str()),
            reinterpret_cast<const unsigned char*>(map.c_str()), &total);
        if (!ok) {
          return nb::object(nb::none());
        }
        return nb::cast(total);
      },
      "sequence"_a, "table"_a, "map"_a,
      "The sum inside `instability_index`, or None for a residue `DIWV` has no "
      "row for.\n\n"
      "A plain left-to-right accumulation and not a Neumaier one, because the "
      "reference writes `score += value` here where it writes `sum(...)` in "
      "`molecular_weight` and `gravy`.  The caller multiplies by "
      "`10.0 / length`.");

  m.def(
      "count_residues",
      [](nb::handle source, nb::bytes map) {
        const BufferRef buffer(source);
        if (map.size() != 256) {
          throw nb::value_error("map must be 256 bytes");
        }
        std::uint64_t counts[20] = {};
        biofasting::count_residues(buffer.data(), buffer.size(),
                                   reinterpret_cast<const unsigned char*>(
                                       map.c_str()),
                                   counts);
        std::vector<std::uint64_t> out(counts, counts + 20);
        return nb::cast(out);
      },
      "sequence"_a, "map"_a,
      "How many of each of the twenty standard residues, in `AMINO_ACIDS` "
      "order.\n\n"
      "Twenty passes of `str.count` and one pass over the bytes are the same "
      "twenty numbers; the difference is that the reference reads the sequence "
      "twenty times and calls into Python twenty times to get them.  It never "
      "declines: `str.count` is not interested in what it does not find, so a "
      "residue outside the twenty is a zero rather than an error, and `map` says "
      "which bytes those are.");
}
