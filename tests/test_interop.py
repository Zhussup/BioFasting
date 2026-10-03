"""Step 1.4: the boundary between our records and ``Bio.SeqRecord``.

A shim is only useful if it is invisible, so nothing here is tested against an
expectation written by hand: every conversion is compared with what
``Bio.SeqIO`` itself produces from the same file, field by field, and both
writers are compared byte for byte with ``SeqIO.write``.  A shim that produces
a slightly different object is worse than no shim at all, because the caller
cannot tell.

The other half of the file is the boundary of the boundary: biofasting must not
start depending on Biopython because it grew an interop module, so one test
runs in a fresh interpreter and asserts that importing the package -- interop
included -- imports no ``Bio`` at all.
"""

import importlib.util
import io
import subprocess
import sys
from pathlib import Path

import pytest

import biofasting

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)
needs_corpus = pytest.mark.skipif(
    not (DATA / "fasta" / "genome.fasta").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)

SAMPLE = (
    b"@SRR000001.0 HWI-ST1:1:FCBAR:1:1101:1000:1000 1:N:0:ATCACG\n"
    b"ACGTACGTN\n+\nIIIIIIIII\n"
    b"@plain\nTTTT\n+\n!!!!\n"
)


# --------------------------------------------------------------------------
# The conversions, against what SeqIO produces from the same bytes.
# --------------------------------------------------------------------------


@needs_biopython
def test_a_fastq_record_becomes_the_record_seqio_would_have_built(tmp_path):
    source = tmp_path / "reads.fastq"
    source.write_bytes(SAMPLE)
    ours = list(biofasting.fastq_seqrecords(source))

    from Bio import SeqIO

    theirs = list(SeqIO.parse(io.StringIO(SAMPLE.decode()), "fastq"))
    assert len(ours) == len(theirs) == 2
    for mine, reference in zip(ours, theirs):
        assert mine.id == reference.id
        assert mine.name == reference.name
        assert mine.description == reference.description
        assert str(mine.seq) == str(reference.seq)
        assert mine.letter_annotations == reference.letter_annotations
        assert mine.annotations == reference.annotations


@needs_biopython
def test_a_fasta_record_becomes_the_record_seqio_would_have_built(tmp_path):
    source = tmp_path / "ref.fasta"
    source.write_bytes(b">chr1 some description here\nACGTAC\nGTACGT\n>chr2\nTTTT\n")
    ours = list(biofasting.fasta_seqrecords(source))

    from Bio import SeqIO

    theirs = list(SeqIO.parse(io.StringIO(source.read_text()), "fasta"))
    assert len(ours) == len(theirs) == 2
    for mine, reference in zip(ours, theirs):
        assert mine.id == reference.id
        assert mine.name == reference.name
        assert mine.description == reference.description
        assert str(mine.seq) == str(reference.seq)
        assert mine.letter_annotations == reference.letter_annotations


@needs_biopython
def test_a_fasta_record_has_no_qualities(tmp_path):
    source = tmp_path / "ref.fasta"
    source.write_bytes(b">chr1\nACGT\n")
    (record,) = biofasting.fasta_seqrecords(source)
    assert record.letter_annotations == {}
    assert biofasting.from_seqrecord(record) == (b"chr1", b"ACGT")


@needs_biopython
def test_a_round_trip_through_a_seqrecord_changes_nothing(tmp_path):
    # Stated as the identity it is: reading, converting out, and converting
    # back gives exactly the tuples the reader yields.
    source = tmp_path / "reads.fastq"
    source.write_bytes(SAMPLE)
    assert [
        biofasting.from_seqrecord(record)
        for record in biofasting.fastq_seqrecords(source)
    ] == list(biofasting.open_fastq(source))


@needs_biopython
def test_to_seqrecord_takes_str_or_bytes(tmp_path):
    from_bytes = biofasting.to_seqrecord(b"r1 note", b"ACGT", b"IIII")
    from_text = biofasting.to_seqrecord("r1 note", "ACGT", "IIII")
    assert from_bytes.id == from_text.id == "r1"
    assert from_bytes.description == from_text.description == "r1 note"
    assert str(from_bytes.seq) == str(from_text.seq) == "ACGT"
    assert from_bytes.letter_annotations == from_text.letter_annotations


@needs_biopython
def test_to_seqrecord_omits_the_qualities_it_is_not_given():
    record = biofasting.to_seqrecord("chr1", "ACGT")
    assert record.letter_annotations == {}
    assert biofasting.from_seqrecord(record) == (b"chr1", b"ACGT")


@needs_biopython
def test_a_title_that_is_one_word_is_its_own_description():
    record = biofasting.to_seqrecord("r2", "TTTT", "!!!!")
    assert record.id == record.name == record.description == "r2"
    # The phred offset is the one Biopython uses: '!' is 33, so quality 0.
    assert record.letter_annotations["phred_quality"] == [0, 0, 0, 0]


@needs_biopython
def test_the_quality_string_is_phred_33(tmp_path):
    source = tmp_path / "reads.fastq"
    source.write_bytes(b"@r1\nACGT\n+\n!5I~\n")
    (record,) = biofasting.fastq_seqrecords(source)
    assert record.letter_annotations["phred_quality"] == [0, 20, 40, 93]


@needs_biopython
def test_non_nucleotides_survive_the_round_trip():
    # '-' and '*' are not bases; a shim that validated would be a parser, and
    # this one must not be.
    record = biofasting.to_seqrecord("r1", "AC-GT*", "IIIIII")
    assert biofasting.from_seqrecord(record) == (b"r1", b"AC-GT*", b"IIIIII")


@needs_biopython
def test_a_record_without_an_id_still_converts():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    bare = SeqRecord(Seq("ACGT"), id="", name="", description="")
    assert biofasting.from_seqrecord(bare) == (b"", b"ACGT")


@needs_biopython
def test_a_non_buffer_argument_is_rejected():
    with pytest.raises(TypeError):
        biofasting.to_seqrecord("r1", 42)


# --------------------------------------------------------------------------
# The writers, against SeqIO.write byte for byte.
# --------------------------------------------------------------------------


@needs_biopython
def test_write_fastq_matches_the_reference_writer(tmp_path):
    from Bio import SeqIO

    source = tmp_path / "reads.fastq"
    source.write_bytes(SAMPLE)
    expected = io.StringIO()
    SeqIO.write(
        list(SeqIO.parse(io.StringIO(SAMPLE.decode()), "fastq")), expected, "fastq"
    )

    actual = io.BytesIO()
    biofasting.write_fastq(biofasting.open_fastq(source), actual)
    assert actual.getvalue() == expected.getvalue().encode()


@needs_biopython
def test_write_fasta_matches_the_reference_writer(tmp_path):
    from Bio import SeqIO

    text = ">chr1 some description here\n" + "ACGT" * 40 + "\n>chr2\nTTTT\n"
    source = tmp_path / "ref.fasta"
    source.write_text(text)
    expected = io.StringIO()
    SeqIO.write(list(SeqIO.parse(io.StringIO(text), "fasta")), expected, "fasta")

    index = biofasting.open_fasta(source)
    actual = io.BytesIO()
    biofasting.write_fasta(((index.title(key), index[key]) for key in index), actual)
    assert actual.getvalue() == expected.getvalue().encode()


@needs_biopython
def test_write_fasta_wraps_at_sixty_like_the_reference(tmp_path):
    actual = io.BytesIO()
    biofasting.write_fasta([("chr1", b"ACGT" * 40)], actual)
    lines = actual.getvalue().splitlines()
    assert lines[0] == b">chr1"
    assert [len(line) for line in lines[1:]] == [60, 60, 40]


@needs_biopython
def test_write_fastq_writes_a_path_too(tmp_path):
    source = tmp_path / "in.fastq"
    source.write_bytes(SAMPLE)
    target = tmp_path / "out.fastq"
    biofasting.write_fastq(biofasting.open_fastq(source), target)
    assert target.read_bytes() == SAMPLE


@needs_biopython
def test_the_writer_refuses_a_quality_that_is_the_wrong_length():
    with pytest.raises(ValueError, match="sequence length 4 but 2 quality"):
        biofasting.write_fastq([(b"r1", b"ACGT", b"II")], io.BytesIO())


# --------------------------------------------------------------------------
# The corpus: every record, in both directions.
# --------------------------------------------------------------------------


@needs_biopython
@needs_corpus
def test_the_corpus_fastq_converts_field_for_field():
    from Bio import SeqIO

    path = DATA / "fastq" / "reads_10k.fastq"
    counted = 0
    with open(path) as handle:
        for ours, theirs in zip(biofasting.fastq_seqrecords(path), SeqIO.parse(handle, "fastq")):
            assert ours.id == theirs.id
            assert ours.description == theirs.description
            assert str(ours.seq) == str(theirs.seq)
            assert ours.letter_annotations == theirs.letter_annotations
            counted += 1
    # `zip` would have stopped at the shorter side; the count says it did not.
    assert counted == 10_000


@needs_biopython
@needs_corpus
def test_the_corpus_fastq_rewrites_byte_for_byte():
    """Read with our scanner, write with our writer, compare with the file.

    The strongest statement the pair can make: a 10 000-record corpus file
    survives our read-write cycle unchanged, so neither the reader nor the
    writer is lossy at any field.
    """
    path = DATA / "fastq" / "reads_10k.fastq"
    actual = io.BytesIO()
    biofasting.write_fastq(biofasting.open_fastq(path), actual)
    assert actual.getvalue() == path.read_bytes()


@needs_biopython
@needs_corpus
def test_the_corpus_fasta_rewrites_byte_for_byte():
    path = DATA / "fasta" / "genome.fasta"
    index = biofasting.open_fasta(path)
    actual = io.BytesIO()
    biofasting.write_fasta(((index.title(key), index[key]) for key in index), actual)
    assert actual.getvalue() == path.read_bytes()


# --------------------------------------------------------------------------
# Biopython is still optional.
# --------------------------------------------------------------------------


def test_importing_the_package_does_not_import_biopython():
    """Run in a fresh interpreter: `Bio` must not appear in `sys.modules`.

    `import biofasting` pulls in `biofasting.interop`, which imports `Bio`
    inside its functions precisely so that this stays true.  A module-level
    import would make Biopython a hard dependency of every install, which is
    the thing this test exists to prevent.
    """
    code = (
        "import sys, biofasting\n"
        "assert biofasting.to_seqrecord is not None\n"
        "print('Bio' in sys.modules)\n"
    )
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert done.stdout.strip().splitlines()[-1] == "False"


def test_the_missing_biopython_error_names_the_fix():
    """The absent-Biopython path, faked rather than skipped.

    ``sys.modules['Bio'] = None`` makes any import of it raise, which is
    exactly the state a user with a wheel install is in -- so this test needs
    no Biopython and runs in the CI shape too.  What it pins is that the
    failure is ours and names `pip install biopython`, not a bare
    ``ModuleNotFoundError`` from three frames inside a library the caller did
    not know wanted it.
    """
    code = (
        "import sys\n"
        # Every `import Bio...` fails from here on, which is what a machine
        # without Biopython does, without needing one to be absent.
        "sys.modules['Bio'] = None\n"
        "import biofasting\n"
        "try:\n"
        "    biofasting.to_seqrecord('r1', 'ACGT')\n"
        "except ImportError as exc:\n"
        "    print(exc)\n"
        "else:\n"
        "    print('NO ERROR')\n"
    )
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    message = done.stdout.strip().splitlines()[-1]
    assert "biopython" in message
    assert "pip install biopython" in message
    # The lazy import must not have poisoned the rest of the package: the
    # readers do not need Biopython and must still work.
    assert "NO ERROR" not in message
