"""`biofasting.arrow`, against `Bio.SeqIO` and against the readers underneath.

The tables are a compatibility surface rather than a kernel, so what is under
test is not "is it fast" (that is `bench/bench_arrow.py`, and it is gated the same
way) but **is it the same table**, and the two things that word covers here:

* the values, against `SeqIO.parse` on the corpus and on hand-written edge
  cases.  The `description` column is the one that needs saying out loud: it is
  the title *after* the key, which is `polars-bio`'s reading and not
  `SeqRecord.description`'s.  A test that compared it against the whole title
  would be testing the wrong contract and would pass on a table no polars-bio
  user recognises;
* the *absence* of a copy.  The whole reason the build is in C++ is that
  `Array.from_buffers` wraps the bytes rather than converting them, so the test
  that matters is that the address pyarrow reports for a column's buffers is the
  address of the numpy views `_core.ArrowTable.buffers()` handed over.  A
  regression that copied would leave every value correct and every number in the
  benchmark wrong, which is exactly the kind of bug this file exists to catch.

The edge cases are all from the reference's own behaviour rather than from the
tables' documentation, and the two divergences from `polars-bio` are pinned by
name in a section of their own -- a deliberate difference that is not tested is
a bug with a story.  That section also builds `polars-bio`'s own table on the
same file and compares them, because "compatible with `polars-bio`" is a claim
about a second library and a test of it that never imports that library would be
a test of this one's opinion of it.
"""

import gzip
import importlib.util
import io
import sys
from pathlib import Path

import pytest

import biofasting
from biofasting import _core
from biofasting.arrow import read_fasta_table, read_fastq_table

pa = pytest.importorskip("pyarrow")

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_corpus = pytest.mark.skipif(
    not (DATA / "fastq" / "reads_10k.fastq").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)
needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

FASTQ_COLUMNS = ["name", "description", "sequence", "quality"]
FASTA_COLUMNS = ["name", "description", "sequence"]


# --------------------------------------------------------------------------
# Helpers: the reference's columns, built the way a Biopython user would.
# --------------------------------------------------------------------------


def reference_columns(path, kind):
    """`SeqIO.parse` into lists, column by column, as `polars-bio` names them.

    `description` is the title minus the key -- and therefore *not*
    `record.description`, which keeps the key.  The subtraction is what
    `polars-bio`'s column means; see the module docstring.

    The quality is the ASCII string, built with the same per-base dictionary
    lookup `SeqIO.QualityIO.FastqPhredWriter` uses, because Biopython exposes
    `letter_annotations["phred_quality"]` as *integers* and a comparison against
    an integer column would pass or fail for the wrong reason.
    """
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import _phred_to_sanger_quality_str

    handle = io.StringIO(path.read_bytes().decode("latin-1"))
    columns = {name: [] for name in (FASTQ_COLUMNS if kind == "fastq" else FASTA_COLUMNS)}
    for record in SeqIO.parse(handle, kind):
        description = record.description
        if description.startswith(record.id):
            description = description[len(record.id) :].lstrip()
        columns["name"].append(record.id)
        columns["description"].append(description)
        columns["sequence"].append(str(record.seq))
        if kind == "fastq":
            columns["quality"].append(
                "".join(
                    _phred_to_sanger_quality_str[score]
                    for score in record.letter_annotations["phred_quality"]
                )
            )
    return columns


def table_columns(table):
    return {name: table.column(name).to_pylist() for name in table.column_names}


def read(path, kind):
    return read_fastq_table(str(path)) if kind == "fastq" else read_fasta_table(str(path))


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_bytes(text if isinstance(text, bytes) else text.encode("latin-1"))
    return path


# --------------------------------------------------------------------------
# The corpus, against the reference.
# --------------------------------------------------------------------------


@needs_corpus
@needs_biopython
def test_fastq_table_matches_the_reference_on_the_corpus():
    path = DATA / "fastq" / "reads_10k.fastq"
    table = read_fastq_table(str(path))
    assert table.column_names == FASTQ_COLUMNS
    assert table.num_rows == 10_000
    assert table_columns(table) == reference_columns(path, "fastq")


@needs_corpus
@needs_biopython
def test_fasta_table_matches_the_reference_on_the_corpus():
    path = DATA / "fasta" / "genome_1mb.fasta"
    table = read_fasta_table(str(path))
    assert table.column_names == FASTA_COLUMNS
    assert table_columns(table) == reference_columns(path, "fasta")


@needs_corpus
def test_the_table_is_the_same_one_the_readers_build():
    """The tables agree with `open_fastq` / `open_fasta`, the non-Arrow path.

    Two implementations of one parse is what this checks for: the table is built
    in C++ in one pass and the readers hand back record by record, so a
    divergence between them would be a second grammar in the package -- the
    thing `fastq.hpp` says it has exactly one of.
    """
    fastq = DATA / "fastq" / "reads_10k.fastq"
    records = list(biofasting.open_fastq(str(fastq)))
    table = table_columns(read_fastq_table(str(fastq)))
    assert table["sequence"] == [sequence.decode("latin-1") for _, sequence, _ in records]
    assert table["quality"] == [quality.decode("latin-1") for _, _, quality in records]
    assert table["name"] == [title.split(b" ")[0].decode("latin-1") for title, _, _ in records]

    fasta = DATA / "fasta" / "genome_1mb.fasta"
    index = biofasting.open_fasta(str(fasta))
    table = table_columns(read_fasta_table(str(fasta)))
    assert table["name"] == [name for name in index]
    assert table["sequence"] == [index[name].decode("latin-1") for name in index]
    assert table["description"] == [
        index.title(name)[len(name) :].lstrip() for name in index
    ]


@needs_corpus
def test_gzipped_fastq_is_sniffed_and_read():
    """A `.gz` is decompressed by its magic bytes, not by its name."""
    path = DATA / "fastq" / "reads_1m.fastq.gz"
    if not path.exists():
        pytest.skip("gzipped corpus row not generated")
    table = read_fastq_table(str(path))
    assert table.num_rows == 1_000_000
    plain = read_fastq_table(str(DATA / "fastq" / "reads_10k.fastq"))
    # The corpus is one generator's output, so the first reads of the 1M file are
    # the first reads of the 10k file; comparing a prefix keeps this honest
    # without inflating a 1M-row table twice.
    assert table.column("name").to_pylist()[:10_000] == plain.column("name").to_pylist()


@needs_corpus
def test_gzip_is_sniffed_by_bytes_and_not_by_extension(tmp_path):
    """A gzipped file called `.fastq` is still read, and vice versa."""
    source = (DATA / "fastq" / "reads_10k.fastq").read_bytes()
    misnamed = write(tmp_path, "reads.fastq", gzip.compress(source))
    table = read_fastq_table(str(misnamed))
    assert table.num_rows == 10_000


# --------------------------------------------------------------------------
# The layout, and the copy that is not there.
# --------------------------------------------------------------------------


def test_the_offsets_are_arrow_offsets(tmp_path):
    """`length + 1` uint64 offsets, ascending, the last equal to the data size."""
    path = write(tmp_path, "a.fastq", "@r1 one\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n")
    built = _core.ArrowTable.from_fastq(path.read_bytes())
    assert built.length == 2
    assert built.column_count == 4
    assert built.names == FASTQ_COLUMNS
    for offsets, data in built.buffers():
        assert offsets.dtype == "uint64"
        assert len(offsets) == built.length + 1
        assert offsets[0] == 0
        assert list(offsets) == sorted(offsets)
        assert offsets[-1] == len(data)
        assert offsets.flags.writeable is False
        assert data.flags.writeable is False


def test_the_arrays_are_the_table_s_own_buffers(tmp_path):
    """`from_buffers` wraps the bytes; it does not convert them.

    The addresses are the test.  A copy would leave every value in this file
    correct and would make the zero-copy hand-off -- the only reason the build
    is in C++ rather than a `pa.table` call -- a claim instead of a fact.

    It has to be the *same* built table on both sides, which is why this one
    goes through `_as_table` rather than `read_fastq_table`: two calls to the
    readers build two tables, and comparing their addresses would compare two
    allocations.
    """
    from biofasting.arrow import _as_table

    path = write(tmp_path, "a.fastq", "@r1 one\nACGT\n+\nIIII\n")
    built = _core.ArrowTable.from_fastq(path.read_bytes())
    table = _as_table(built)
    for index, (offsets, data) in enumerate(built.buffers()):
        array = table.column(index).chunk(0)
        validity, arrow_offsets, arrow_data = array.buffers()
        assert validity is None, "nothing here has a null in it"
        assert arrow_offsets.address == offsets.ctypes.data
        assert arrow_data.address == data.ctypes.data
        assert arrow_offsets.size == (built.length + 1) * 8

    # And the arrays are what keeps those bytes alive: the table is dropped, and
    # the values still read.  numpy holds the views, the views hold the table.
    import gc

    del built
    gc.collect()
    assert table.column("name").to_pylist() == ["r1"]


def test_the_table_owns_its_bytes(tmp_path):
    """Dropping everything but the table leaves every value readable.

    The table is the one reader here that does not point into the caller's
    buffer -- a column has to be contiguous -- so the file, the mapping and the
    built object can all go and the values must not.
    """
    import gc

    path = write(tmp_path, "a.fastq", "@r1 one\nACGTACGT\n+\nIIIIIIII\n")
    table = biofasting.read_fastq_table(str(path))
    path.unlink()
    gc.collect()
    assert table.column("sequence").to_pylist() == ["ACGTACGT"]
    assert table.column("quality").to_pylist() == ["IIIIIIII"]


def test_polars_takes_the_table_without_copying_the_strings():
    """`pl.from_arrow` is the documented hand-off, and the bytes are not copied.

    Not "and it is free": polars 1.44 indexes an Arrow string column rather than
    adopting it, which is 32-38 ns a record (measured in `bench/bench_arrow.py`).
    What is pinned here is the part that matters to the layout -- the frame reads
    correctly after the table object is gone, so the bytes it is reading are the
    ones this package built and not the table wrapper's.

    Skipped rather than required: polars is not a dependency of anything here.
    """
    polars = pytest.importorskip("polars")
    path = DATA / "fastq" / "reads_10k.fastq"
    if not path.exists():
        pytest.skip("benchmark corpus not generated")
    import gc

    table = read_fastq_table(str(path))
    frame = polars.from_arrow(table)
    del table
    gc.collect()
    assert frame.columns == FASTQ_COLUMNS
    assert frame.height == 10_000
    assert frame["name"][-1] == "SRR000001.9999"
    assert len(frame["sequence"][0]) == 150
    assert len(frame["quality"][0]) == 150


# --------------------------------------------------------------------------
# The reading of a header: the key, and what follows it.
# --------------------------------------------------------------------------


def test_the_description_is_the_title_after_the_key(tmp_path):
    path = write(tmp_path, "a.fastq", "@r1 a description\nACGT\n+\nIIII\n")
    table = read_fastq_table(str(path))
    assert table.column("name").to_pylist() == ["r1"]
    assert table.column("description").to_pylist() == ["a description"]


def test_a_header_with_no_description_is_the_empty_string(tmp_path):
    """Not a null, and not the key repeated: `polars-bio` writes null here.

    Both are defensible readings of a header that says only `>rec4`, and this
    one is pinned because switching between the two libraries must not turn a
    value into an absence.
    """
    path = write(tmp_path, "a.fasta", ">rec4\nACGT\n")
    table = read_fasta_table(str(path))
    assert table.column("description").to_pylist() == [""]
    assert table.column("description").null_count == 0


def test_a_run_of_whitespace_between_key_and_description_disappears(tmp_path):
    """A tab separates, and a run of spaces is taken off the front.

    Matched against `polars-bio`'s own reading, which is the point: a caller
    moving a frame from one library to the other must not have to rewrite the
    header.
    """
    path = write(
        tmp_path,
        "a.fasta",
        ">rec1 simple desc\nACGT\n>rec2\ttab desc\nACGT\n>rec3  two  spaces\nACGT\n",
    )
    table = read_fasta_table(str(path))
    assert table.column("name").to_pylist() == ["rec1", "rec2", "rec3"]
    assert table.column("description").to_pylist() == [
        "simple desc",
        "tab desc",
        "two  spaces",
    ]


def test_a_title_of_only_whitespace_has_no_key(tmp_path):
    path = write(tmp_path, "a.fasta", ">   \nACGT\n")
    table = read_fasta_table(str(path))
    assert table.column("name").to_pylist() == [""]
    assert table.column("description").to_pylist() == [""]


def test_a_leading_space_is_read_where_polars_bio_refuses_the_file(tmp_path):
    """`>  rec5 leading` is a record here; `polars-bio` refuses the whole file.

    The reference reads it and keys it `rec5` -- `SeqIO.index` strips first --
    and a library that cannot open a file Biopython opens is a regression
    whatever it is compatible with.  This divergence is therefore deliberate.
    """
    path = write(tmp_path, "a.fasta", ">  rec5 leading\nACGT\n")
    table = read_fasta_table(str(path))
    assert table.column("name").to_pylist() == ["rec5"]
    assert table.column("description").to_pylist() == ["leading"]


# --------------------------------------------------------------------------
# FASTA: the bytes, and the refusals.
# --------------------------------------------------------------------------


def test_a_wrapped_fasta_sequence_loses_its_newlines_and_its_spaces(tmp_path):
    """What `SimpleFastaParser` builds, including the two easy-to-miss rules.

    A space inside a line is removed and a tab is kept -- the reference rstripes
    each line, joins them, then removes the spaces -- and a `;` inside a record
    is sequence text rather than a comment.
    """
    path = write(
        tmp_path,
        "a.fasta",
        ">rec1 desc\nAC GT\nAC\tGT\nAC;GT\n",
    )
    table = read_fasta_table(str(path))
    assert table.column("sequence").to_pylist() == ["ACGTAC\tGTAC;GT"]


def test_a_duplicate_key_refuses_the_whole_table(tmp_path):
    path = write(tmp_path, "a.fasta", ">rec1 a\nACGT\n>rec1 b\nTTTT\n")
    with pytest.raises(ValueError, match="Duplicate key 'rec1'"):
        read_fasta_table(str(path))


def test_a_gzipped_fasta_is_refused_by_name(tmp_path):
    """Compression is sniffed and refused, not indexed as if it were text."""
    path = write(tmp_path, "a.fasta", gzip.compress(b">rec1\nACGT\n"))
    with pytest.raises(ValueError, match="gzip-compressed"):
        read_fasta_table(str(path))


def test_a_file_that_is_not_fasta_has_no_records(tmp_path):
    """Not an error: the reference has no records there and neither does this."""
    path = write(tmp_path, "a.fasta", "not a fasta file at all\n")
    table = read_fasta_table(str(path))
    assert table.num_rows == 0
    assert table.column_names == FASTA_COLUMNS


# --------------------------------------------------------------------------
# FASTQ: the parse, and the failure that carries its location.
# --------------------------------------------------------------------------


def test_quality_is_the_ascii_string_the_file_carries(tmp_path):
    path = write(tmp_path, "a.fastq", "@r1\nACGT\n+\n!\"#$\n")
    table = read_fastq_table(str(path))
    assert table.column("quality").to_pylist() == ['!"#$']


def test_a_wrapped_record_is_joined_the_way_the_parser_joins_it(tmp_path):
    path = write(tmp_path, "a.fastq", "@r1 one\nAC\nGT\n+\nII\nII\n")
    table = read_fastq_table(str(path))
    assert table.column("sequence").to_pylist() == ["ACGT"]
    assert table.column("quality").to_pylist() == ["IIII"]


def test_a_malformed_record_raises_with_its_location(tmp_path):
    """The parser's message, plus the record number and the byte offset.

    The location is what the reference cannot give and what a malformed file of
    a million records makes worth having, and the table builder reports it
    through the same function the scanner's own binding uses -- two copies of
    that wording would be two ways for one byte to read.
    """
    path = write(tmp_path, "a.fastq", "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+")
    with pytest.raises(ValueError) as caught:
        read_fastq_table(str(path))
    assert "record 1, byte offset" in str(caught.value)


def test_an_empty_file_is_an_empty_table(tmp_path):
    path = write(tmp_path, "a.fastq", "")
    table = read_fastq_table(str(path))
    assert table.num_rows == 0
    assert table.column_names == FASTQ_COLUMNS
    assert table.column("sequence").to_pylist() == []


# --------------------------------------------------------------------------
# The compatibility story, against the library it is a story about.
# --------------------------------------------------------------------------


needs_polars_bio = pytest.mark.skipif(
    importlib.util.find_spec("polars_bio") is None,
    reason="polars-bio is not installed; it is the compatibility target, not a dependency",
)


@needs_corpus
@needs_polars_bio
def test_the_columns_are_the_ones_polars_bio_builds():
    """Both libraries' tables, on the same file, column by column.

    This is the claim PLAN 2.3 is about -- a frame moved from one library to the
    other must not have to be rewritten -- and the only honest way to test it is
    to build both and compare, because the columns are the contract and a
    contract tested against a fixture is a contract tested against its author.

    Two differences are known, named in `biofasting.arrow`'s module docstring,
    and normalised here rather than in the library.  The FASTQ quality column is
    `quality_scores` there, and a record with no description is a **null** there
    and `""` here.  Neither is corrected in either direction: the first would
    break this package's own naming, the second would put the only null in a
    package that has none.  The null is filled here so that the rest of the
    comparison can be strict, and it is asserted as a divergence of its own
    below, so that this fill cannot be mistaken for agreement.
    """
    import polars_bio as pb

    fastq = str(DATA / "fastq" / "reads_10k.fastq")
    theirs = pb.read_fastq(fastq).rename({"quality_scores": "quality"})
    ours = read_fastq_table(fastq)
    assert ours.column_names == theirs.columns
    for name in ours.column_names:
        assert ours.column(name).to_pylist() == theirs[name].fill_null("").to_list()

    fasta = str(DATA / "fasta" / "genome_1mb.fasta")
    theirs = pb.read_fasta(fasta)
    ours = read_fasta_table(fasta)
    assert ours.column_names == theirs.columns
    for name in ours.column_names:
        assert ours.column(name).to_pylist() == theirs[name].fill_null("").to_list()


@needs_polars_bio
def test_a_header_with_no_description_is_empty_here_and_null_there(tmp_path):
    """The divergence, pinned by name, so that it is a contract and not a bug.

    A length-zero value and the absence of a value are different things, and a
    caller who moves a frame across libraries meets the difference here first.
    """
    import polars_bio as pb

    path = write(tmp_path, "a.fasta", ">rec4\nACGT\n")
    assert pb.read_fasta(str(path))["description"].to_list() == [None]
    table = read_fasta_table(str(path))
    assert table.column("description").to_pylist() == [""]
    assert table.column("description").null_count == 0


@needs_polars_bio
def test_a_leading_space_is_read_here_where_polars_bio_refuses_the_file(tmp_path):
    """The second divergence, and the reason it is deliberate.

    `polars-bio` refuses the whole file with `FASTA read error: missing name`;
    the reference reads the record and keys it `rec5`.  A library that cannot
    open a file Biopython opens is a regression whatever it is compatible with,
    so this one is a copy of the reference's behaviour and not a departure from
    it.
    """
    import polars.exceptions
    import polars_bio as pb

    path = write(tmp_path, "a.fasta", ">  rec5 leading\nACGT\n")
    with pytest.raises(polars.exceptions.ComputeError, match="missing name"):
        pb.read_fasta(str(path))
    table = read_fasta_table(str(path))
    assert table.column("name").to_pylist() == ["rec5"]
    assert table.column("description").to_pylist() == ["leading"]


# --------------------------------------------------------------------------
# The optional dependency.
# --------------------------------------------------------------------------


def test_available_reports_whether_pyarrow_is_importable():
    assert biofasting.arrow.available() is True


def test_a_missing_pyarrow_is_explained_rather_than_imported(monkeypatch):
    """The message names the extra, not the module that is missing.

    The interesting case is a caller who installed this package from a wheel
    and has never had Arrow: `ModuleNotFoundError: No module named 'pyarrow'`
    from inside a library they did not know wanted one is a worse answer than
    the command that fixes it.
    """
    monkeypatch.setitem(sys.modules, "pyarrow", None)
    with pytest.raises(ImportError, match=r"pip install biofasting\[arrow\]"):
        biofasting.arrow._pyarrow()


def test_the_tables_do_not_import_pyarrow_until_one_is_asked_for():
    """Importing the package is not importing pyarrow.

    Checked in a fresh interpreter rather than by looking at `sys.modules`: the
    point is what a caller pays when they never ask for a table, and the test
    that ran before this one has already imported pyarrow in this process.
    """
    import subprocess

    program = (
        "import sys, biofasting;"
        "assert 'pyarrow' not in sys.modules, 'pyarrow was imported eagerly';"
        "assert 'polars' not in sys.modules, 'polars was imported eagerly'"
    )
    result = subprocess.run(
        [sys.executable, "-c", program], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
