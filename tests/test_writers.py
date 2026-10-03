"""M23: the sequence-format writers, gated against ``SeqIO.write`` byte for byte.

A writer is a different kind of kernel from a reader.  A reader can be checked
field by field, because a field is a small thing and a disagreement about one
is easy to see; a writer's output is one long byte string and the only honest
check is the whole of it.  So every comparison here is
``bytes == SeqIO.write(...)`` and never a hand-written expectation -- the
reference is the specification, and a writer that emits a slightly different
file is a bug the caller cannot detect, because both files parse.

The three formats are not equally interesting and the tests say so:

- **FASTA** is the reference at its Python floor -- one ``handle.write`` per
  60-base slice -- so the kernel is checked mostly on the branches a corpus
  never reaches: an empty sequence with and without wrapping, ``wrap`` of 0 and
  ``None``, and the newline the unwrapped branch writes for no sequence at all.
- **FASTQ** is four lines and no wrapping, and its edge is the length check:
  the message names the record, and it stays in Python because the kernel
  cannot raise a Python exception about a record it never described.
- **QUAL** is where the reference is worst and where the interesting branch
  lives: ``write_record`` cuts lines with ``data.rfind(" ", 0, wrap)`` -- its
  *fast* branch, not the ``pop(0)`` loop ``to_string`` uses -- and a width of
  one to five takes a third branch that is reproduced in Python.  All three are
  compared at every width, because "the kernel is byte-identical at 60" says
  nothing about 5.

The one thing not tested is the kernel's refusal when no space falls inside a
window: it is unreachable from here, and the test that says so is a property
rather than a ``pytest.raises`` -- see `test_a_window_always_holds_a_space`.
"""

import importlib.util
import io
import random
import warnings
from pathlib import Path

import pytest

import biofasting

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)
needs_corpus = pytest.mark.skipif(
    not (DATA / "fasta" / "genome_1mb.fasta").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)
needs_fastq_corpus = pytest.mark.skipif(
    not (DATA / "fastq" / "reads_10k.fastq").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)

pytestmark = needs_biopython

FASTA_1MB = DATA / "fasta" / "genome_1mb.fasta"
FASTQ_10K = DATA / "fastq" / "reads_10k.fastq"

_SAMPLE_FASTQ = (
    b"@SRR000001.0 HWI-ST1:1:FCBAR:1:1101:1000:1000 1:N:0:ATCACG\n"
    b"ACGTACGTN\n+\nIIIIIIIII\n"
    b"@plain\nTTTT\n+\n!!!!\n"
)


# --------------------------------------------------------------------------
# Helpers: the reference writer, and the two ways to reach a wrap width.
# --------------------------------------------------------------------------


_DEFAULT = object()
"""Sentinel for "the width `SeqIO.write` picks", so that a `wrap` of ``None``
means the reference's *no wrapping* branch rather than the default one.  The two
are different outputs and a test that confused them would pass for the wrong
reason."""


def reference(fmt: str, records, wrap=_DEFAULT) -> bytes:
    """``SeqIO.write``, optionally at a wrap width the format string cannot set.

    ``SeqIO.write`` only takes a format *string*, so a non-default width means
    instantiating the writer class itself.  It is the same call ``SeqIO.write``
    makes -- ``writer_class(handle).write_file(sequences)`` -- so the comparison
    stays against the reference and not against a re-implementation of it.
    """
    from Bio import SeqIO

    handle = io.StringIO()
    if wrap is _DEFAULT:
        SeqIO.write(records, handle, fmt)
    elif fmt == "fasta":
        from Bio.SeqIO.FastaIO import FastaWriter

        FastaWriter(handle, wrap=wrap).write_file(records)
    elif fmt == "qual":
        from Bio.SeqIO.QualityIO import QualPhredWriter

        QualPhredWriter(handle, wrap=wrap).write_file(records)
    else:
        raise AssertionError(fmt)
    return handle.getvalue().encode("latin-1")


def seqrecords(rows, qualities=False):
    """Our rows as ``SeqRecord`` objects, for the reference to write."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    out = []
    for row in rows:
        title = row[0].decode("latin-1")
        record = SeqRecord(
            Seq(row[1].decode("latin-1")),
            id=title.split()[0] if title.split() else "",
            description=title,
        )
        if qualities:
            record.letter_annotations["phred_quality"] = [
                byte - 33 for byte in row[-1]
            ]
        out.append(record)
    return out


class CountingHandle(io.BytesIO):
    """A binary handle that counts its ``write`` calls.

    The whole reason the writers are kernels is that the file is one buffer
    written once, so the count is not a detail of the implementation -- it is
    the property being bought, and it is asserted rather than assumed.
    """

    def __init__(self):
        super().__init__()
        self.calls = 0

    def write(self, data):
        self.calls += 1
        return super().write(data)


def phred33(scores):
    """Scores as the phred+33 byte string our records carry."""
    return bytes(min(126, score + 33) for score in scores)


def random_rows(seed, count, length, low=0, high=80):
    rng = random.Random(seed)
    rows = []
    for index in range(count):
        scores = [rng.randrange(low, high) for _ in range(length)]
        rows.append(
            (
                f"read{index} seed={seed}".encode(),
                b"ACGTN" * (length // 5) + b"ACGTN"[: length % 5],
                phred33(scores),
            )
        )
    return rows


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------


def test_write_fasta_matches_the_reference_on_a_one_megabyte_corpus():
    """The stored corpus, through our index and our writer, against SeqIO.

    Not the read-write identity (that is `test_interop`'s) but a comparison with
    the reference's own bytes over a megabyte of wrapped records, which is where
    a wrap boundary that is off by one would show up a few hundred thousand
    times.
    """
    index = biofasting.open_fasta(FASTA_1MB)
    rows = [(index.title(key).encode(), index[key]) for key in index]

    actual = io.BytesIO()
    biofasting.write_fasta(rows, actual)
    assert actual.getvalue() == reference("fasta", seqrecords(rows))


def test_write_fasta_wraps_at_every_width_the_reference_accepts():
    rows = [(b"chr1 note", b"ACGTACGTAC" * 7), (b"chr2", b"TTTT")]
    records = seqrecords(rows)
    for wrap in (1, 2, 5, 6, 7, 59, 60, 61, 70, 700):
        actual = io.BytesIO()
        biofasting.write_fasta(rows, actual, wrap=wrap)
        assert actual.getvalue() == reference("fasta", records, wrap), wrap


def test_write_fasta_wrap_of_zero_and_none_is_one_line_a_record():
    rows = [(b"chr1", b"ACGT" * 30)]
    actual = io.BytesIO()
    biofasting.write_fasta(rows, actual, wrap=0)
    assert actual.getvalue() == b">chr1\n" + b"ACGT" * 30 + b"\n"
    # `None` is the reference's other spelling of the same branch, and both
    # must land on it rather than on "wrap at 60".
    assert actual.getvalue() == reference("fasta", seqrecords(rows), wrap=0)
    again = io.BytesIO()
    biofasting.write_fasta(rows, again, wrap=None)
    assert again.getvalue() == actual.getvalue()


def test_an_empty_sequence_writes_nothing_when_wrapped_and_a_newline_when_not():
    """The reference's two branches, which disagree about an empty record.

    ``for i in range(0, 0, 60)`` never runs, so a wrapped empty record is a
    title line and nothing else; the unwrapped branch writes ``data + "\\n"``
    and so emits a blank line.  A writer that picked one behaviour for both
    would be one byte wrong on a legal file, and the difference is asserted
    against the reference rather than against a memory of it.
    """
    rows = [(b"chr1", b"")]
    wrapped = io.BytesIO()
    biofasting.write_fasta(rows, wrapped)
    assert wrapped.getvalue() == b">chr1\n"
    assert wrapped.getvalue() == reference("fasta", seqrecords(rows))

    unwrapped = io.BytesIO()
    biofasting.write_fasta(rows, unwrapped, wrap=0)
    assert unwrapped.getvalue() == b">chr1\n\n"
    assert unwrapped.getvalue() == reference("fasta", seqrecords(rows), wrap=0)


def test_write_fasta_refuses_a_negative_wrap():
    # `FastaWriter.__init__` raises ValueError for a negative width, and the
    # kernel takes an unsigned one, so the refusal has to happen before it.
    with pytest.raises(ValueError, match="negative"):
        biofasting.write_fasta([(b"chr1", b"ACGT")], io.BytesIO(), wrap=-1)


def test_write_fasta_takes_any_iterable_not_just_a_list():
    """Generators, so a caller never materialises the file to write it."""
    index = biofasting.open_fasta(FASTA_1MB)
    rows = ((index.title(key).encode(), index[key]) for key in index)
    actual = io.BytesIO()
    biofasting.write_fasta(rows, actual)
    assert actual.getvalue() == FASTA_1MB.read_bytes()


# --------------------------------------------------------------------------
# FASTQ
# --------------------------------------------------------------------------


def test_write_fastq_matches_the_reference_on_the_ten_thousand_read_corpus():
    rows = list(biofasting.open_fastq(FASTQ_10K))
    assert len(rows) == 10_000
    actual = io.BytesIO()
    biofasting.write_fastq(rows, actual)
    assert actual.getvalue() == reference("fastq", seqrecords(rows, qualities=True))


def test_write_fastq_is_a_round_trip_on_the_corpus():
    actual = io.BytesIO()
    biofasting.write_fastq(biofasting.open_fastq(FASTQ_10K), actual)
    assert actual.getvalue() == FASTQ_10K.read_bytes()


def test_write_fastq_is_four_lines_and_no_wrapping():
    rows = [(b"r1 a description", b"A" * 200, b"I" * 200)]
    actual = io.BytesIO()
    biofasting.write_fastq(rows, actual)
    lines = actual.getvalue().split(b"\n")
    assert lines[0] == b"@r1 a description"
    assert lines[1] == b"A" * 200
    assert lines[2] == b"+"
    assert lines[3] == b"I" * 200
    assert actual.getvalue() == reference("fastq", seqrecords(rows, qualities=True))


def test_the_fastq_length_check_names_the_record():
    with pytest.raises(ValueError, match="record b'r1' has sequence length 4 but 2"):
        biofasting.write_fastq([(b"r1", b"ACGT", b"II")], io.BytesIO())


def test_the_fastq_length_check_is_the_same_for_every_offending_record():
    # Not just the first: the check is per record, so a later one is caught too
    # and the message is about the one that was wrong.
    with pytest.raises(ValueError, match="record b'r2'"):
        biofasting.write_fastq(
            [(b"r1", b"ACGT", b"IIII"), (b"r2", b"ACGT", b"IIIIIII")], io.BytesIO()
        )


def test_write_fastq_accepts_a_record_with_no_bases():
    rows = [(b"empty", b"", b"")]
    actual = io.BytesIO()
    biofasting.write_fastq(rows, actual)
    assert actual.getvalue() == b"@empty\n\n+\n\n"
    assert actual.getvalue() == reference("fastq", seqrecords(rows, qualities=True))


def test_write_fastq_takes_str_fields_and_rejects_what_is_neither():
    """`str` is encoded for the caller; anything else is refused by name.

    The Python layer accepts text because every other function here does, and
    the kernel below it does not: a caller reaching `_core.write_fastq` directly
    gets told the field must be bytes, which is the contract the kernel actually
    has.
    """
    from biofasting import _core

    actual = io.BytesIO()
    biofasting.write_fastq([("r1", "ACGT", "IIII")], actual)
    assert actual.getvalue() == b"@r1\nACGT\n+\nIIII\n"

    with pytest.raises(TypeError, match="quality must be str or bytes-like, not int"):
        biofasting.write_fastq([(b"r1", b"ACGT", 42)], io.BytesIO())

    with pytest.raises(TypeError, match="a FASTQ sequence must be bytes"):
        _core.write_fastq([(b"r1", "ACGT", b"IIII")])


# --------------------------------------------------------------------------
# QUAL
# --------------------------------------------------------------------------


def test_write_qual_matches_the_reference_at_every_wrap_width():
    """Every branch of the reference's ``write_record``, at every width.

    Zero and ``None`` are the no-wrap branch, one to five the ``pop(0)`` loop
    the kernel does not have, and six and up the ``rfind`` cut.  Widths around
    the boundaries are included because that is where an off-by-one would hide:
    5/6 is the branch change and 59/60/61 is a line exactly one token over.
    """
    rows = random_rows(seed=1, count=3, length=97)
    records = seqrecords(rows, qualities=True)
    pairs = [(title, quality) for title, _, quality in rows]
    for wrap in list(range(0, 66)) + [None, 70, 128]:
        actual = io.BytesIO()
        biofasting.write_qual(pairs, actual, wrap=wrap)
        assert actual.getvalue() == reference("qual", records, wrap=wrap), wrap


def test_write_qual_matches_the_reference_on_the_ten_thousand_read_corpus():
    rows = list(biofasting.open_fastq(FASTQ_10K))
    actual = io.BytesIO()
    biofasting.write_qual([(title, quality) for title, _, quality in rows], actual)
    assert actual.getvalue() == reference("qual", seqrecords(rows, qualities=True))


def test_write_qual_at_sixty_reads_back_as_the_scores_it_started_with():
    """A QUAL file the reference parses gives back the same numbers.

    Byte-identity with the reference already implies this, but it is stated
    separately because it is the property a user has: read FASTQ, write QUAL,
    parse it, and get the quality scores.
    """
    from Bio import SeqIO

    rows = random_rows(seed=2, count=5, length=150)
    actual = io.BytesIO()
    biofasting.write_qual([(title, quality) for title, _, quality in rows], actual)
    parsed = list(SeqIO.parse(io.StringIO(actual.getvalue().decode()), "qual"))
    assert [record.letter_annotations["phred_quality"] for record in parsed] == [
        [byte - 33 for byte in quality] for _, _, quality in rows
    ]


def test_write_qual_with_no_wrapping_is_one_line_a_record():
    rows = random_rows(seed=3, count=2, length=40)
    pairs = [(title, quality) for title, _, quality in rows]
    actual = io.BytesIO()
    biofasting.write_qual(pairs, actual, wrap=0)
    lines = actual.getvalue().split(b"\n")
    assert lines[0] == b">read0 seed=3"
    # Every score but the first is preceded by a space, so one line of forty
    # scores holds thirty-nine of them.
    assert lines[1].count(b" ") == 39
    assert actual.getvalue() == reference(
        "qual", seqrecords(rows, qualities=True), wrap=0
    )
    # `None` is the reference's other spelling of the same branch, and must land
    # on it rather than on the default width.
    again = io.BytesIO()
    biofasting.write_qual(pairs, again, wrap=None)
    assert again.getvalue() == actual.getvalue()


def test_a_record_with_no_qualities_writes_a_title_and_the_reference_s_newline():
    """An empty record, in all three of the reference's branches.

    The two wrapping branches write a blank line -- ``len(data) <= wrap`` holds
    for the empty string -- while the ``pop(0)`` branch never enters its loop
    and writes the title alone.  Three outputs from one record, and each is
    compared with the reference rather than described.
    """
    rows = [(b"r1", b"")]
    records = seqrecords(rows, qualities=True)
    for wrap, expected in ((0, b">r1\n\n"), (None, b">r1\n\n"), (60, b">r1\n\n"), (5, b">r1\n")):
        actual = io.BytesIO()
        biofasting.write_qual(rows, actual, wrap=wrap)
        assert actual.getvalue() == expected, wrap
        assert actual.getvalue() == reference("qual", records, wrap=wrap), wrap


def test_a_window_always_holds_a_space_at_the_widths_the_kernel_sees():
    """The kernel's refusal cannot be reached, and this is why.

    ``append_qual`` returns false where ``data.rfind(" ", 0, wrap)`` finds no
    space in the window -- the place Biopython 1.88 would write the same string
    forever.  A score is at most three digits (a quality byte is at most 222)
    and every token but the first is preceded by a space, so the first space is
    never past index three, while the kernel only sees widths of six and up.
    The property is asserted directly rather than trusted, because the comment
    in the reference -- "unless we have X digit or higher quality scores!" --
    names exactly the input that would falsify it.
    """
    from biofasting import _core

    # Score 222 is 255 less the offset: the widest token a byte can produce.
    rows = [(b"r1", bytes([255]) * 30)]
    text = " ".join(["222"] * 30)
    assert len(max(text.split(" "), key=len)) == 3
    for wrap in range(6, 40):
        assert text.rfind(" ", 0, wrap - 1) != -1, wrap
        assert _core.write_qual(rows, wrap) == reference(
            "qual", _records_with_scores([222] * 30), wrap=wrap
        ), wrap


def _records_with_scores(scores):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("A" * len(scores)), id="r1", description="r1")
    record.letter_annotations["phred_quality"] = list(scores)
    return [record]


def test_write_qual_refuses_a_negative_wrap():
    with pytest.raises(ValueError, match="negative"):
        biofasting.write_qual([(b"r1", b"IIII")], io.BytesIO(), wrap=-1)


# --------------------------------------------------------------------------
# phred_to_sanger
# --------------------------------------------------------------------------


def test_phred_to_sanger_is_the_reference_table():
    from Bio.SeqIO.QualityIO import _phred_to_sanger_quality_str

    for score in range(94):
        assert biofasting._core.phred_to_sanger(bytes([score])) == (
            _phred_to_sanger_quality_str[score].encode()
        ), score


def test_phred_to_sanger_truncates_above_ninety_three_like_the_reference():
    """The reference's slow path, reproduced rather than corrected.

    ``chr(min(126, int(round(qp)) + 33))`` is what ``_get_sanger_quality_str``
    falls back to for a score its 0..93 table has no entry for, and it caps at
    ``~``.  A kernel that raised instead would be refusing input the reference
    accepts.
    """
    assert biofasting._core.phred_to_sanger(bytes([94, 100, 200, 255])) == b"~~~~"


def test_from_seqrecord_encodes_qualities_through_the_kernel():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ACGT"), id="r1", description="r1")
    record.letter_annotations["phred_quality"] = [0, 20, 40, 93]
    assert biofasting.from_seqrecord(record) == (b"r1", b"ACGT", b"!5I~")


def test_from_seqrecord_declines_scores_a_byte_cannot_hold():
    """A float or a ``None`` leaves the kernel and takes the reference's path.

    ``bytes(qualities)`` raises for both, and the reference's own slow
    expression is the only thing that rounds a float the way it rounds one --
    ``int(round(qp))`` -- and that raises the ``TypeError`` a ``None`` deserves.
    """
    from Bio.Seq import Seq
    from Bio.SeqIO.QualityIO import _get_sanger_quality_str
    from Bio.SeqRecord import SeqRecord

    def make(qualities):
        record = SeqRecord(Seq("A" * len(qualities)), id="r1", description="r1")
        record.letter_annotations["phred_quality"] = list(qualities)
        return record

    for qualities in ([0, 20.0, 40, 93], [0, 20.5, 40, 93], [0, 20.4, 40, 93]):
        record = make(qualities)
        # The reference itself is the expectation, warning behaviour included:
        # the differential is not against a memory of the rounding rule.
        assert biofasting.from_seqrecord(record) == (
            b"r1",
            b"A" * len(qualities),
            _get_sanger_quality_str(record).encode(),
        ), qualities

    with pytest.raises(TypeError, match="A quality value of None was found"):
        biofasting.from_seqrecord(make([0, None, 40, 93]))


def test_a_score_above_the_sanger_ceiling_is_truncated_and_warned_about():
    """The reference's data-loss warning, which our writers must also raise.

    PHRED 94 has no Sanger character, so the reference writes ``~`` and warns.
    With warnings as errors -- which is how this suite runs -- reproducing the
    warning is the difference between a caller seeing the truncation and not.
    """
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ACGT"), id="r1", description="r1")
    record.letter_annotations["phred_quality"] = [0, 20, 40, 100]
    with pytest.warns(biofasting.BiopythonWarning, match="Data loss"):
        assert biofasting.from_seqrecord(record) == (b"r1", b"ACGT", b"!5I~")

    # A score exactly at the ceiling is not data loss and must stay quiet.
    record.letter_annotations["phred_quality"] = [0, 20, 40, 93]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert biofasting.from_seqrecord(record) == (b"r1", b"ACGT", b"!5I~")


# --------------------------------------------------------------------------
# The property the kernels exist for, and the paths they write to.
# --------------------------------------------------------------------------


def test_each_writer_reaches_the_handle_once():
    """One crossing for the whole file, which is the point of a kernel here.

    The reference writes a FASTA record in 60-base slices -- eighteen calls for
    a kilobase -- and a QUAL record a line at a time; a kernel that kept the
    loop would have removed the parsing and left the crossings.
    """
    fasta_rows = [(b"chr1", b"ACGT" * 250)]
    fastq_rows = [(b"r1", b"ACGT" * 40, b"I" * 160)]
    qual_rows = [(b"r1", b"I" * 160)]

    for writer, rows in (
        (biofasting.write_fasta, fasta_rows),
        (biofasting.write_fastq, fastq_rows),
        (biofasting.write_qual, qual_rows),
    ):
        handle = CountingHandle()
        writer(rows, handle)
        assert handle.calls == 1, writer.__name__
        assert handle.getvalue()


def test_the_writers_take_a_path_as_well_as_a_handle(tmp_path):
    source = tmp_path / "reads.fastq"
    source.write_bytes(_SAMPLE_FASTQ)

    fasta = tmp_path / "out.fasta"
    biofasting.write_fasta(
        [(title, sequence) for title, sequence, _ in biofasting.open_fastq(source)],
        fasta,
    )
    assert fasta.read_bytes() == b">SRR000001.0 HWI-ST1:1:FCBAR:1:1101:1000:1000 1:N:0:ATCACG\nACGTACGTN\n>plain\nTTTT\n"

    qual = tmp_path / "out.qual"
    biofasting.write_qual(
        [(title, quality) for title, _, quality in biofasting.open_fastq(source)], qual
    )
    assert qual.read_bytes() == b">SRR000001.0 HWI-ST1:1:FCBAR:1:1101:1000:1000 1:N:0:ATCACG\n40 40 40 40 40 40 40 40 40\n>plain\n0 0 0 0\n"

    again = tmp_path / "roundtrip.fastq"
    biofasting.write_fastq(biofasting.open_fastq(source), again)
    assert again.read_bytes() == _SAMPLE_FASTQ
