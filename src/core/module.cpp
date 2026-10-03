// The nanobind module.  Deliberately thin: it declares the names and hands off
// to build_info.cpp and cpu_features.cpp, so that adding a kernel in step 1.2
// does not mean editing the place where the ABI is defined.

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "biofasting/build_config.hpp"
#include "annotations.hpp"
#include "build_info.hpp"
#include "cpu_features.hpp"
#include "fasta.hpp"
#include "fastq.hpp"
#include "feature.hpp"
#include "flatfile.hpp"
#include "inflate.hpp"
#include "location.hpp"
#include "protparam.hpp"
#include "seqops.hpp"
#include "sequtils.hpp"
#include "write.hpp"

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

// The canonical tuple spellings of the location kernel's structures.  Kept as
// plain tuples rather than as bound classes because the Python layer already
// has value types for them (biofasting.seqfeature): a second set of classes
// here would be a second spelling of the same thing, free to drift from the
// first, and the conversion to `Bio.SeqFeature` objects would then have two
// inputs to disagree about.
//
//   position := (kind: int, value: int, left: int, right: int,
//                choices: tuple[int, ...])
//   part     := (start: position, end: position, strand: int, ref: str)
//   location := (operator: str, parts: tuple[part, ...])
//
// `kind` indexes biofasting.seqfeature.POSITION_KINDS and `strand` is -1, +1 or
// 0 for the reference's `None`: a location on a protein has no strand at all,
// which is a different statement from "stranded, strand unknown".
nb::tuple position_tuple(const biofasting::Position& position) {
  nb::list choices;
  for (const std::int64_t choice : position.choices) choices.append(choice);
  return nb::make_tuple(
      static_cast<int>(position.kind), position.value, position.left,
      position.right, nb::steal<nb::tuple>(PyList_AsTuple(choices.ptr())));
}

nb::tuple part_tuple(const biofasting::LocationPart& part) {
  return nb::make_tuple(position_tuple(part.start), position_tuple(part.end),
                        static_cast<int>(part.strand),
                        nb::str(part.ref.c_str(), part.ref.size()));
}

const char* location_op_name(biofasting::LocationOp op) noexcept {
  switch (op) {
    case biofasting::LocationOp::simple: return "simple";
    case biofasting::LocationOp::join: return "join";
    case biofasting::LocationOp::order: return "order";
    case biofasting::LocationOp::bond: return "bond";
  }
  return "simple";
}

nb::tuple location_tuple(const biofasting::Location& location) {
  nb::list parts;
  for (const biofasting::LocationPart& part : location.parts) {
    parts.append(part_tuple(part));
  }
  return nb::make_tuple(nb::str(location_op_name(location.op)),
                        nb::steal<nb::tuple>(PyList_AsTuple(parts.ptr())));
}

const char* location_status_name(biofasting::LocationStatus status) noexcept {
  switch (status) {
    case biofasting::LocationStatus::ok: return "ok";
    case biofasting::LocationStatus::parser_error: return "parser_error";
    case biofasting::LocationStatus::refused: return "refused";
  }
  return "refused";
}

const char* location_warning_name(biofasting::LocationWarningKind kind) noexcept {
  switch (kind) {
    case biofasting::LocationWarningKind::origin_wrap: return "origin_wrap";
    case biofasting::LocationWarningKind::bond: return "bond";
  }
  return "bond";
}

// The warnings a parse would have raised, in the order the reference would have
// raised them: a list of `(kind, text)` and not a pair of flags, because
// `join(30..5,60..2)` on a circular record warns twice and
// `join(bond(1),30..5)` warns in that order.  `text` is the part's own text for
// a wrap -- which is what the reference quotes -- and empty for a bond, whose
// wording quotes nothing.  The Python layer says the words.
nb::list location_warnings(const std::vector<biofasting::LocationWarning>& warnings) {
  nb::list out;
  for (const biofasting::LocationWarning& warning : warnings) {
    out.append(nb::make_tuple(nb::str(location_warning_name(warning.kind)),
                              nb::str(warning.text.c_str(), warning.text.size())));
  }
  return out;
}

// A feature's qualifiers as `(key, value, has_value, escape_warning,
// escape_text)`.  The last two are the reference's NCBI escaping warning: the
// flag, and the value the warning is about -- which is the value as it stood
// before the doubled quotes were undone, not the value kept.  The Python layer
// decides whether to say the words; the kernel is what noticed.
nb::list feature_qualifiers(const biofasting::Feature& feature) {
  nb::list out;
  for (const biofasting::FeatureQualifier& qualifier : feature.qualifiers) {
    out.append(nb::make_tuple(
        nb::str(qualifier.key.c_str(), qualifier.key.size()),
        nb::str(qualifier.value.c_str(), qualifier.value.size()),
        qualifier.has_value, qualifier.escape_warning,
        nb::str(qualifier.escape_text.c_str(), qualifier.escape_text.size())));
  }
  return out;
}

// The format name as `SeqIO.parse` spells it, resolved once at construction so
// that a typo is an error at the call and not a record read as the wrong shape.
biofasting::FlatFileFormat flatfile_format(const std::string& name) {
  if (name == "genbank") return biofasting::FlatFileFormat::genbank;
  if (name == "embl") return biofasting::FlatFileFormat::embl;
  if (name == "swiss") return biofasting::FlatFileFormat::swiss;
  throw nb::value_error(("unknown flat-file format '" + name +
                         "': expected 'genbank', 'embl' or 'swiss'")
                            .c_str());
}

// The Python face of FlatFileIndex.  Same buffer-lifetime rule as FastaIndex:
// the index points into the caller's buffer, so the reference is declared first
// and held for as long as the index lives.
//
// Note what is *not* here: no SeqRecord and no annotations dict.  What this
// exposes is the three fields the sixth ranking pass measured as the part of a
// flat-file parse that is neither data movement nor object construction, plus
// the FEATURES table read on demand, and the interop shim is what turns them
// into a record.
class PyFlatFileIndex {
 public:
  PyFlatFileIndex(nb::handle source, const std::string& format)
      : buffer_(source),
        index_(buffer_.data(), buffer_.size(), flatfile_format(format)) {
    if (index_.failed()) throw nb::value_error(index_.error_message().c_str());
  }

  std::size_t size() const noexcept { return index_.size(); }

  // The format the index was built as, by the name `SeqIO.parse` uses.  Not an
  // internal detail: the two flat-file formats carry different key sets in
  // their annotations dict -- EMBL has no `date`, no `source`, and no `keywords`
  // unless a `KW` line was there -- and the layer that decides the naming is one
  // level up, so it has to be able to ask which format it is naming for rather
  // than carry a second copy of the answer.
  std::string format() const {
    return biofasting::flatfile_format_name(index_.format());
  }

  bool has(const std::string& id) const { return index_.find(id) != nullptr; }

  std::string name(const std::string& id) const { return require(id).name; }

  std::string description(const std::string& id) const {
    return require(id).description;
  }

  std::size_t sequence_length(const std::string& id) const {
    return require(id).length;
  }

  nb::bytes get(const std::string& id) const {
    const std::string sequence = index_.sequence(require(id));
    return nb::bytes(sequence.data(), sequence.size());
  }

  nb::bytes slice(const std::string& id, std::size_t start,
                  std::size_t end) const {
    const std::string window = index_.sequence_slice(require(id), start, end);
    return nb::bytes(window.data(), window.size());
  }

  nb::list keys() const {
    nb::list ids;
    for (const biofasting::FlatFileRecord& record : index_.records()) {
      ids.append(nb::str(record.id.c_str(), record.id.size()));
    }
    return ids;
  }

  // The whole file in one call, in file order: (id, name, description,
  // sequence).  This is what a stream parse wants -- a caller that iterates
  // every record should not pay a lookup per record for an order it already
  // knows, and the C++ side can build the list without a Python call per field.
  nb::list records() const {
    nb::list out;
    for (const biofasting::FlatFileRecord& record : index_.records()) {
      const std::string sequence = index_.sequence(record);
      out.append(nb::make_tuple(
          nb::str(record.id.c_str(), record.id.size()),
          nb::str(record.name.c_str(), record.name.size()),
          nb::str(record.description.c_str(), record.description.size()),
          nb::bytes(sequence.data(), sequence.size())));
    }
    return out;
  }

  // Where a record's parts are, in the file's own terms.  Not an internal
  // detail: it is what says whether a reader is walking the bytes once or
  // re-parsing them, and a test asserting it is what stops the single-pass
  // property from quietly disappearing.
  nb::tuple location(const std::string& id) const {
    const biofasting::FlatFileRecord& record = require(id);
    return nb::make_tuple(record.offset, record.sequence_offset,
                          record.sequence_end, record.length);
  }

  nb::list offsets() const {
    nb::list out;
    for (const biofasting::FlatFileRecord& record : index_.records()) {
      out.append(record.offset);
    }
    return out;
  }

  // `(ok, message, features)`.  Each feature is `(type, location|None, status,
  // message, warnings, qualifiers)`, with the location as
  // `_core.parse_location` returns one -- `None` when the reference's parser
  // rejected the string, which is a behaviour its consumer turns into a missing
  // location and a warning -- and the warnings in the reference's own order.
  //
  // The first element is False when this reader declined the table: the record's
  // header is not in a layout it reproduces, or a line is shaped the way the
  // reference only warns about.  A record with no feature block at all is
  // `ok` with an empty list, which is an answer and not a refusal.
  nb::tuple features(const std::string& id) const {
    const biofasting::FeatureTable table = index_.features(require(id));
    nb::list features;
    for (const biofasting::Feature& feature : table.features) {
      nb::object location = nb::none();
      if (feature.location.status == biofasting::LocationStatus::ok) {
        location = location_tuple(feature.location.location);
      }
      features.append(nb::make_tuple(
          nb::str(feature.type.c_str(), feature.type.size()), location,
          nb::str(location_status_name(feature.location.status)),
          nb::str(feature.location.message.c_str(),
                  feature.location.message.size()),
          location_warnings(feature.location.warnings),
          feature_qualifiers(feature)));
    }
    return nb::make_tuple(!table.refused,
                          nb::str(table.message.c_str(), table.message.size()),
                          features);
  }

  // `(ok, message, keys, table)` for the record's header.
  //
  // `keys` is the list of annotation keys the header *created*, in the order it
  // created them, and it is not a formality: the reference inserts each key when
  // the line stating it is consumed, so a `COMMENT` above the first `REFERENCE`
  // puts `comment` before `references`, and a dict compares and prints by that
  // order.  It is also what says a key exists at all -- `data_file_division` is
  // created for every EMBL ID line even if the field is blank, where GenBank's
  // `molecule_type` is skipped when its column is empty -- so the table carries
  // values and `keys` carries the set, and an empty string in the table means
  // what the reference means by it.
  //
  // The values are the reading, not the naming; the Python layer builds the dict
  // and is where the per-format differences are stated.  The references come back
  // as tuples of the eight fields the reference's `Reference` holds, with the
  // location already in Python coordinates -- `(start - 1, end)`, which is the
  // one piece of arithmetic in the block and is done once, here.
  nb::tuple annotations(const std::string& id) const {
    const biofasting::AnnotationTable table =
        index_.annotations(require(id));
    nb::list keys;
    for (const biofasting::AnnotationKey key : table.order) {
      keys.append(nb::str(biofasting::annotation_key_name(key)));
    }
    nb::dict out;
    out["molecule_type"] = nb::cast(table.molecule_type);
    out["topology"] = nb::cast(table.topology);
    out["data_file_division"] = nb::cast(table.data_file_division);
    out["date"] = nb::cast(table.date);
    out["accessions"] = strings(table.accessions);
    out["keywords"] = strings(table.keywords);
    out["taxonomy"] = strings(table.taxonomy);
    out["source"] = nb::cast(table.source);
    out["organism"] = nb::cast(table.organism);
    out["comment"] = nb::cast(table.comment);
    out["sequence_version"] = nb::cast(table.sequence_version);
    out["gi"] = nb::cast(table.gi);
    nb::list references;
    for (const biofasting::AnnotationReference& reference : table.references) {
      nb::list location;
      for (const std::pair<std::int64_t, std::int64_t>& part : reference.location) {
        location.append(nb::make_tuple(part.first, part.second));
      }
      references.append(nb::make_tuple(
          nb::str(reference.title.c_str(), reference.title.size()),
          nb::str(reference.authors.c_str(), reference.authors.size()),
          nb::str(reference.consrtm.c_str(), reference.consrtm.size()),
          nb::str(reference.journal.c_str(), reference.journal.size()),
          nb::str(reference.pubmed_id.c_str(), reference.pubmed_id.size()),
          nb::str(reference.medline_id.c_str(), reference.medline_id.size()),
          nb::str(reference.comment.c_str(), reference.comment.size()),
          location));
    }
    out["references"] = references;
    return nb::make_tuple(!table.failed,
                          nb::str(table.error_message.c_str(),
                                  table.error_message.size()),
                          keys, out);
  }

 private:
  static nb::list strings(const std::vector<std::string>& values) {
    nb::list out;
    for (const std::string& value : values) {
      out.append(nb::str(value.c_str(), value.size()));
    }
    return out;
  }

  const biofasting::FlatFileRecord& require(const std::string& id) const {
    const biofasting::FlatFileRecord* record = index_.find(id);
    if (record == nullptr) throw nb::key_error(id.c_str());
    return *record;
  }

  BufferRef buffer_;
  biofasting::FlatFileIndex index_;
};

// The Python face of FastqIndex.  Same buffer-lifetime rule as FastaIndex --// and here it is load-bearing twice over, because the index does not merely
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

// One field of one record row, borrowed rather than copied.  The rows are
// tuples of `bytes` -- the shape every reader in this package produces -- and a
// `str` here is a caller error rather than something to encode, because only
// the caller knows which encoding produced its text.
std::string_view tuple_bytes(const nb::tuple& row, Py_ssize_t index,
                             const char* what) {
  if (index >= static_cast<Py_ssize_t>(row.size())) {
    throw nb::value_error("each record must be a tuple of bytes");
  }
  PyObject* cell = PyTuple_GET_ITEM(row.ptr(), index);
  if (!PyBytes_Check(cell)) throw nb::type_error(what);
  return std::string_view(PyBytes_AS_STRING(cell),
                          static_cast<std::size_t>(PyBytes_GET_SIZE(cell)));
}

// Walk any Python iterable of record tuples and build one format's whole output
// in a single buffer.
//
// The iteration is the Python C API's and not a conversion to `std::vector` on
// purpose: a vector would copy every title, sequence and quality into a C++
// string before the first byte of output, and a caller writing a million reads
// would pay for a second copy of all of them.  Each row is borrowed for the
// length of one record and appended straight into the result, so the memory the
// writer needs is the file it is writing and nothing else.
template <typename Emit>
nb::bytes write_all(nb::handle records, Emit emit) {
  std::string out;
  PyObject* iterator = PyObject_GetIter(records.ptr());
  if (iterator == nullptr) throw nb::python_error();
  PyObject* item = nullptr;
  while ((item = PyIter_Next(iterator)) != nullptr) {
    const nb::object holder = nb::steal(item);
    emit(out, nb::cast<nb::tuple>(holder));
  }
  Py_DECREF(iterator);
  if (PyErr_Occurred()) throw nb::python_error();
  return nb::bytes(out.data(), out.size());
}

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

  nb::class_<PyFlatFileIndex>(
      m, "FlatFileIndex",
      "Indexes the records of a GenBank, EMBL or SwissProt file held in any "
      "object exposing the buffer protocol.\n\n"
      "Keys are record ids in file order, with a duplicate refused, so a fetch "
      "is a copy of the record's bytes rather than a re-parse.\n\n"
      "What it reproduces is `id`, `name`, `description`, the sequence, the "
      "FEATURES table and the header values the annotations dict is built from "
      "-- the parts the flat-file measurement says are the reference's own logic "
      "rather than its data movement or its object building.  What it does not "
      "reproduce is `SeqRecord.dbxrefs` or the record's `Seq` object, so it is "
      "not a drop-in for Bio.SeqIO.parse and must not be presented as one.\n\n"
      "It is also stricter than the reference on purpose.  Biopython meets a "
      "malformed sequence line by warning and guessing -- a wrongly indented "
      "GenBank line is shifted by one byte and parsed anyway -- and a guess "
      "reproduced differently is a sequence that looks right and is not.  This "
      "reader refuses such a file with the byte offset that made it refuse, and "
      "the caller falls back to Bio.SeqIO.parse.\n"
      "Use biofasting.genbank.open_genbank() rather than building this "
      "directly: it is the piece that keeps the buffer alive.")
      .def(nb::init<nb::handle, const std::string&>(), "source"_a, "format"_a,
           "Index `source`, which must expose the buffer protocol, as one of "
           "'genbank', 'embl' or 'swiss'.  Raises ValueError for an unknown "
           "format, for a duplicate record id, or for a file this reader "
           "declines to guess at.")
      .def("__len__", &PyFlatFileIndex::size, "Number of records.")
      .def("format", &PyFlatFileIndex::format,
           "The format this index was built as: 'genbank', 'embl' or 'swiss'.")
      .def("__contains__", &PyFlatFileIndex::has, "id"_a,
           "Whether a record with this id exists.")
      .def("__getitem__", &PyFlatFileIndex::get, "id"_a,
           "The record's sequence, as bytes.  Raises KeyError if there is no "
           "such record.")
      .def("__iter__",
           [](const PyFlatFileIndex& self) { return nb::iter(self.keys()); },
           "Iterate the record ids, in file order.")
      .def("keys", &PyFlatFileIndex::keys,
           "The record ids, in file order, as a list.")
      .def("name", &PyFlatFileIndex::name, "id"_a,
           "The record's name: the LOCUS name for GenBank, the ID line's first "
           "field for EMBL and SwissProt.")
      .def("description", &PyFlatFileIndex::description, "id"_a,
           "The record's description: DEFINITION or DE, joined and, for "
           "GenBank, without its trailing full stop.")
      .def("sequence_length", &PyFlatFileIndex::sequence_length, "id"_a,
           "The sequence length, counted during the scan rather than measured "
           "from the bytes handed back.")
      .def("sequence_slice", &PyFlatFileIndex::slice, "id"_a, "start"_a, "end"_a,
           "The [start, end) window of the sequence.  Both bounds are clamped "
           "to the record.")
      .def("records", &PyFlatFileIndex::records,
           "Every record in file order, as (id, name, description, sequence).\n\n"
           "This is what a stream parse wants: a caller iterating the whole "
           "file should not pay a lookup per record for an order it already "
           "knows.")
      .def("location", &PyFlatFileIndex::location, "id"_a,
           "Where the record is: (offset, sequence_offset, sequence_end, "
           "length), as byte offsets into the buffer.")
      .def("offsets", &PyFlatFileIndex::offsets,
           "Every record's starting byte offset, in file order.")
      .def("features", &PyFlatFileIndex::features, "id"_a,
           "The record's FEATURES table, as (ok, message, features).\n\n"
           "Each feature is (type, location_or_None, status, message,\n"
           "warnings, qualifiers), and each qualifier is\n"
           "(key, value, has_value, escape_warning, escape_text).  `warnings` is\n"
           "what the reference would have warned about while parsing the\n"
           "location, as (kind, text) pairs in its own order -- see\n"
           "parse_location.  `ok` is False\n"
           "where this reader declined the table -- a shape the reference only\n"
           "warns about, or a header whose declared size and topology it could\n"
           "not read -- and `message` then says which line it was.")
      .def("annotations", &PyFlatFileIndex::annotations, "id"_a,
           "The record's header, as (ok, message, keys, table).\n\n"
           "`keys` is the list of annotation keys the header *created*, in the\n"
           "order it created them: the reference inserts a key when the line\n"
           "stating it is consumed, so a COMMENT above the first REFERENCE puts\n"
           "`comment` before `references`, and a dict compares and prints by\n"
           "that order.  Membership in `keys` is also what says a key exists at\n"
           "all -- `data_file_division` is created for every EMBL ID line even\n"
           "when the field is blank, where GenBank's `molecule_type` is skipped\n"
           "when its column is empty -- so an empty value in `table` means what\n"
           "the reference means by it rather than standing in for an absence.\n\n"
           "`table` holds the values the header stated.  `references` is a list\n"
           "of (title, authors, consrtm, journal, pubmed_id, medline_id,\n"
           "comment, location) with the location already in Python coordinates.\n"
           "Which of these become keys of SeqRecord.annotations is decided one\n"
           "level up, because it differs between GenBank and EMBL.  `ok` is\n"
           "False where this reader declined the header, and `message` names the\n"
           "shape it declined.");

  // ---------------------------------------------------------------------
  // Feature locations
  // ---------------------------------------------------------------------
  //
  // `Location.fromstring` and `Position.fromstring`, as kernels.  They return a
  // status rather than raising, because the reference's three outcomes are not
  // three exceptions: it catches `LocationParserError` inside the feature
  // consumer and turns it into `location = None` plus a warning, so a kernel
  // that raised for it would erase the difference between "unparseable, and the
  // reference says so quietly" and "this reader declines to guess".
  //
  // `parse_location(text, length, circular, stranded)` returns
  //
  //   (status, location_or_None, message, warnings)
  //
  // with `status` one of "ok", "parser_error" (the reference raises
  // LocationParserError) or "refused" (the reference raises something else, or
  // accepts the string while dropping text from it).  `length` is the record's
  // declared size and may be None, which is not the same fact as zero.
  m.def(
      "parse_location",
      [](const std::string& text, nb::object length, bool circular,
         bool stranded) {
        std::int64_t declared = 0;
        const bool has_length = !length.is_none();
        if (has_length) declared = nb::cast<std::int64_t>(length);
        const biofasting::LocationResult result = biofasting::parse_location(
            text, declared, has_length, circular, stranded);
        nb::object location = nb::none();
        if (result.status == biofasting::LocationStatus::ok) {
          location = location_tuple(result.location);
        }
        return nb::make_tuple(location_status_name(result.status), location,
                              nb::str(result.message.c_str(),
                                      result.message.size()),
                              location_warnings(result.warnings));
      },
      "text"_a, "length"_a = nb::none(), "circular"_a = false,
      "stranded"_a = true,
      "Read a feature location such as "
      "\"complement(join(490883..490885,1..879))\" into canonical positions "
      "and parts.\n\n"
      "Returns (status, location, message, warnings).  "
      "A location is not a pair of integers, and this does not flatten it into "
      "one: `complement` carries a strand, `<`/`>` make an end fuzzy, `(3.9)` "
      "is a boundary known only to lie between two bases, `one-of(...)` is a "
      "choice, and `^` is a zero-length junction.\n\n"
      "The status is \"ok\", \"parser_error\" -- where the reference raises "
      "LocationParserError, which its feature consumer catches and turns into a "
      "missing location plus a warning -- or \"refused\", where the reference "
      "raises something else or silently discards text.  A refusal is a design "
      "difference and not a bug: it is this reader saying it will not guess.\n\n"
      "`warnings` is what the reference would have warned about *inside* the "
      "parse, reported rather than spelled and in its own order: `(kind, text)` "
      "with kind \"origin_wrap\" (one entry per part repaired as origin "
      "wrapping, carrying the part's own text, which the reference quotes) or "
      "\"bond\" (one entry per dropped `bond` part, whose wording is fixed, so "
      "its text is empty).");

  m.def(
      "parse_position",
      [](const std::string& text, int offset) {
        biofasting::Position position;
        std::string message;
        const biofasting::LocationStatus status =
            biofasting::parse_position(text, offset, position, message);
        nb::object out = nb::none();
        if (status == biofasting::LocationStatus::ok) {
          out = position_tuple(position);
        }
        return nb::make_tuple(location_status_name(status), out,
                              nb::str(message.c_str(), message.size()));
      },
      "text"_a, "offset"_a = 0,
      "Read one end of a location: `Position.fromstring(text, offset)`.\n\n"
      "`offset` is 0 for an end position and -1 for a start position, which is "
      "the reference's own convention and not a convenience -- it is what makes "
      "`(9.10)` mean 9 at the start of a location and 10 at its end.  Anything "
      "else is refused.\n\n"
      "Returns (status, position, message) with the same status vocabulary as "
      "parse_location.");

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

  // --- the writers ---------------------------------------------------------
  //
  // Each takes the whole iterable of records and returns the whole file as one
  // `bytes`.  Not one call per record: the reference already writes one record
  // at a time and the cost being removed is the interpreter crossing, so a
  // kernel reached once per record would pay the crossing back.  The rows are
  // `(title, sequence)` for FASTA, `(title, sequence, quality)` for FASTQ and
  // `(title, quality)` for QUAL, in this package's own form -- bytes, the marker
  // character excluded, and a quality string already phred+33.

  m.def(
      "write_fasta",
      [](nb::handle records, std::size_t wrap) {
        return write_all(records, [wrap](std::string& out, const nb::tuple& row) {
          biofasting::append_fasta(
              out, tuple_bytes(row, 0, "a FASTA title must be bytes"),
              tuple_bytes(row, 1, "a FASTA sequence must be bytes"), wrap);
        });
      },
      "records"_a, "wrap"_a = 60,
      "The whole FASTA file as one `bytes`, from `(title, sequence)` rows.\n\n"
      "`title` excludes the `>` and lines are wrapped at `wrap`, which is "
      "Biopython's default and the width its own writer produces.  A `wrap` of "
      "0 is not the same program: the reference takes its unwrapped branch "
      "there and writes a newline even for an empty sequence, where a wrapped "
      "record of no length writes no base line at all.");

  m.def(
      "write_fastq",
      [](nb::handle records) {
        return write_all(records, [](std::string& out, const nb::tuple& row) {
          const std::string_view sequence =
              tuple_bytes(row, 1, "a FASTQ sequence must be bytes");
          const std::string_view quality =
              tuple_bytes(row, 2, "a FASTQ quality string must be bytes");
          if (sequence.size() != quality.size()) {
            throw nb::value_error(
                "a FASTQ record's sequence and quality must be the same "
                "length");
          }
          biofasting::append_fastq(
              out, tuple_bytes(row, 0, "a FASTQ title must be bytes"), sequence,
              quality);
        });
      },
      "records"_a,
      "The whole FASTQ file as one `bytes`, from `(title, sequence, quality)` "
      "rows.\n\n"
      "Four lines a record and no wrapping, which is what Biopython does and "
      "what every consumer of the format expects.  The quality string is copied "
      "rather than decoded: our records already carry phred+33, and a kernel "
      "that re-encoded it would be undoing work the caller did.");

  m.def(
      "write_qual",
      [](nb::handle records, std::size_t wrap) {
        return write_all(records, [wrap](std::string& out, const nb::tuple& row) {
          if (!biofasting::append_qual(
                  out, tuple_bytes(row, 0, "a QUAL title must be bytes"),
                  tuple_bytes(row, 1, "a QUAL quality string must be bytes"),
                  wrap)) {
            throw nb::value_error(
                "this record's scores cannot be wrapped at this width: the "
                "reference's own line cut needs a space inside every window, "
                "and Biopython 1.88 loops forever where there is none");
          }
        });
      },
      "records"_a, "wrap"_a = 60,
      "The whole QUAL file as one `bytes`, from `(title, quality)` rows.\n\n"
      "The scores are decimal, joined by single spaces and cut at the last one "
      "inside each `wrap`-wide window, which is the reference's "
      "`data.rfind(' ', 0, wrap)` -- its *fast* wrapping branch, the one "
      "`SeqIO.write` takes at the default width, and not the `pop(0)` loop its "
      "`to_string` uses.  A `wrap` of 0 means one line a record.");

  m.def(
      "phred_to_sanger",
      [](nb::handle source) {
        const BufferRef buffer(source);
        std::string out;
        biofasting::phred_to_sanger(
            reinterpret_cast<const unsigned char*>(buffer.data()), buffer.size(),
            out);
        return nb::bytes(out.data(), out.size());
      },
      "scores"_a,
      "PHRED scores as the Sanger ASCII string a FASTQ or QUAL record carries.\n\n"
      "`min(126, score + 33)` a byte, which is the reference's own encoding "
      "including its truncation at 93: a Sanger FASTQ cannot hold a score above "
      "93, and the reference warns and truncates rather than wrapping.  This "
      "exists because the reference does it with a dictionary lookup a base, "
      "and it is the one place a `SeqRecord`'s `phred_quality` becomes the "
      "string the writers above take.");
}
