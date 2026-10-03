"""Step 1.2: the FASTA index, against Bio.SeqIO.

Same two layers as the FASTQ kernel's tests, and for the same reason.  The
first layer asserts behaviour on hand-written byte strings that need neither
Biopython nor the generated corpus, so CI -- which has neither -- still tests
the index.  The second is the real test: the whole corpus, indexed here and
parsed there, compared record for record, plus the pyfaidx-produced ``.fai``
that was written by a third implementation and is therefore the only piece of
evidence here that does not come from Biopython.

Two things the reference does that are easy to assume away, and are pinned
below rather than described: a ';' comment *inside* a record is sequence text,
and a tab inside a line survives while a space does not.
"""

import gzip
import importlib.util
import io
import random
from pathlib import Path

import pytest

import biofasting

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_corpus = pytest.mark.skipif(
    not (DATA / "fasta" / "genome.fasta").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)
needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)


def index(data):
    """Build an index over `data`, letting a duplicate key raise."""
    return biofasting._core.FastaIndex(data)


def load(path):
    """Run the file-level helper, which is what a user actually calls."""
    return biofasting.open_fasta(path)


def reference(text):
    """The reference's records: (title, sequence) pairs, in file order."""
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    return list(SimpleFastaParser(io.StringIO(text, newline=None)))


def key_of(title):
    """SeqIO.index's rule for turning a title into a key, written out.

    Not a re-implementation of the kernel -- it is the *specification* the
    kernel is held to, taken from the reference's own source:
    ``line[1:].strip().split(None, 1)[0]``.
    """
    stripped = title.strip()
    return stripped.split(None, 1)[0] if stripped else ""


# --------------------------------------------------------------------------
# Behaviour, without the reference.
# --------------------------------------------------------------------------


def test_reads_a_plain_record():
    idx = index(b">rec1\nACGTACGT\n")
    assert len(idx) == 1
    assert list(idx) == ["rec1"]
    assert idx["rec1"] == b"ACGTACGT"
    assert idx.title("rec1") == "rec1"
    assert idx.sequence_length("rec1") == 8


def test_the_key_is_the_first_word_and_the_title_is_the_whole_line():
    idx = index(b">rec1 a description with spaces\nACGT\n")
    assert list(idx) == ["rec1"]
    assert idx.title("rec1") == "rec1 a description with spaces"
    assert idx["rec1"] == b"ACGT"


def test_leading_whitespace_is_skipped_before_the_key():
    # `line[1:].strip().split(None, 1)[0]` strips both ends first, so the key
    # is not simply everything up to the first space.
    idx = index(b">   rec2 desc\nACGT\n")
    assert list(idx) == ["rec2"]
    assert idx.title("rec2") == "   rec2 desc"


def test_records_are_wrapped_and_joined():
    idx = index(b">rec1\nACGT\nACGT\nAC\n")
    assert idx["rec1"] == b"ACGTACGTAC"
    entry = idx.index_entry("rec1")
    assert entry == ("rec1", 10, 6, 4, 5)


def test_a_short_last_line_keeps_the_fast_path():
    """The regression that a correctness-only test cannot see.

    A wrapped record normally ends on a short line, and treating that as a
    broken grid is invisible from the outside -- the fallback assembles the
    same bytes, just slowly, so every other test here stays green while the
    index quietly stops being an index.  This one asserts the flag.
    """
    idx = index(b">rec1\nACGTAC\nACGT\n")
    assert idx["rec1"] == b"ACGTACACGT"
    assert idx.strided("rec1") is True


def test_a_short_line_in_the_middle_is_not_a_grid():
    idx = index(b">rec1\nACGTAC\nAC\nACGTAC\n")
    assert idx["rec1"] == b"ACGTACACACGTAC"
    assert idx.strided("rec1") is False


def test_a_last_line_longer_than_the_grid_is_not_a_grid():
    # Same length as a clean two-line record, so nothing but the line count
    # catches it -- and without that check the bulk copy would read past the
    # record's own bytes.
    idx = index(b">rec1\nACGTAC\nACGTACACGTAC\n")
    assert idx["rec1"] == b"ACGTACACGTACACGTAC"
    assert idx.strided("rec1") is False
    whole = idx["rec1"]
    assert idx.sequence_slice("rec1", 0, len(whole)) == whole


def test_a_blank_line_inside_a_record_is_nothing():
    idx = index(b">rec1\nACGT\n\nACGT\n")
    assert idx["rec1"] == b"ACGTACGT"


def test_comment_lines_before_the_first_record_are_skipped():
    idx = index(b"; a comment\n\n>rec1\nACGT\n")
    assert list(idx) == ["rec1"]
    assert idx["rec1"] == b"ACGT"


def test_a_comment_inside_a_record_is_sequence_text():
    # Not a bug and not a nicety: the reference only treats a line starting
    # with '>' as a boundary, so this is what it builds -- with the space in
    # the comment removed, which is the next surprise.
    idx = index(b">rec1\nACGT\n; a comment\nACGT\n")
    assert idx["rec1"] == b"ACGT;acommentACGT"


def test_a_space_is_removed_and_a_tab_is_kept():
    idx = index(b">rec1\nAC GT\nAC\tGT\n")
    assert idx["rec1"] == b"ACGTAC\tGT"


def test_trailing_whitespace_is_stripped_not_sequence():
    idx = index(b">rec1  \nACGT  \nACGT\t\n")
    assert idx.title("rec1") == "rec1"
    assert idx["rec1"] == b"ACGTACGT"


def test_crlf_ends_a_line_like_the_text_stream_does():
    idx = index(b">rec1\r\nACGTACGT\r\n>rec2\r\nTTTT\r\n")
    assert list(idx) == ["rec1", "rec2"]
    assert idx["rec1"] == b"ACGTACGT"
    assert idx["rec2"] == b"TTTT"


def test_a_lone_carriage_return_also_ends_a_line():
    idx = index(b">rec1\rACGTACGT\r>rec2\rTTTT\r")
    assert list(idx) == ["rec1", "rec2"]
    assert idx["rec2"] == b"TTTT"


def test_the_last_line_may_have_no_terminator():
    idx = index(b">rec1\nACGTACGT")
    assert idx["rec1"] == b"ACGTACGT"


def test_a_record_may_be_empty():
    idx = index(b">empty\n>nonempty\nACGT\n")
    assert list(idx) == ["empty", "nonempty"]
    assert idx["empty"] == b""
    assert idx.sequence_length("empty") == 0


def test_an_empty_file_has_no_records_and_is_not_an_error():
    idx = index(b"")
    assert len(idx) == 0
    assert list(idx) == []


def test_text_that_is_not_fasta_has_no_records():
    # The reference does not validate; it looks for '>' and finds none.
    idx = index(b"ACGT\nACGT\n")
    assert len(idx) == 0


def test_a_duplicate_key_is_rejected_with_the_reference_wording():
    with pytest.raises(ValueError, match=r"Duplicate key 'rec1'"):
        index(b">rec1\nACGT\n>rec1\nTTTT\n")


def test_a_missing_key_raises_key_error():
    idx = index(b">rec1\nACGT\n")
    with pytest.raises(KeyError):
        idx["nope"]
    assert "rec1" in idx
    assert "nope" not in idx


def test_sequence_length_agrees_with_the_sequence():
    idx = index(b">rec1\nACGT\nACGT\nAC\n>rec2\nTT\n")
    for name in idx:
        assert idx.sequence_length(name) == len(idx[name])


def test_slicing_walks_the_grid():
    idx = index(b">rec1\nACGTACGTAC\nACGTACGTAC\n")  # 20 bases, 10 per line
    assert idx["rec1"] == b"ACGTACGTACACGTACGTAC"
    assert idx.sequence_slice("rec1", 8, 12) == b"ACAC"
    assert idx.sequence_slice("rec1", 0, 10) == b"ACGTACGTAC"
    assert idx.sequence_slice("rec1", 10, 20) == b"ACGTACGTAC"
    assert idx.sequence_slice("rec1", 0, 20) == idx["rec1"]


def test_slicing_clamps_rather_than_raising():
    idx = index(b">rec1\nACGTACGT\n")
    assert idx.sequence_slice("rec1", 0, 1000) == b"ACGTACGT"
    assert idx.sequence_slice("rec1", 1000, 2000) == b""
    assert idx.sequence_slice("rec1", 5, 2) == b""
    assert idx.sequence_slice("rec1", 4, 4) == b""


def test_slicing_a_record_that_is_not_a_regular_grid():
    # Two lines of different widths: the sequence is still right, it is just
    # not sliceable by arithmetic, and the fallback has to agree with Python.
    idx = index(b">rec1\nACGTAC\nACGT\n")
    whole = idx["rec1"]
    assert whole == b"ACGTACACGT"
    for start in range(len(whole) + 2):
        for end in range(len(whole) + 2):
            assert idx.sequence_slice("rec1", start, end) == whole[start:end]


def test_a_file_helper_maps_and_reads(tmp_path):
    path = tmp_path / "ref.fasta"
    path.write_bytes(b">rec1 desc\nACGT\nACGT\n>rec2\nTTTT\n")
    idx = load(path)
    assert list(idx) == ["rec1", "rec2"]
    assert idx["rec1"] == b"ACGTACGT"
    assert biofasting.read_fasta(path) == [
        ("rec1", "rec1 desc", b"ACGTACGT"),
        ("rec2", "rec2", b"TTTT"),
    ]


def test_the_file_helper_handles_an_empty_file(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_bytes(b"")
    assert biofasting.read_fasta(path) == []


def test_the_file_helper_refuses_a_compressed_file(tmp_path):
    # The alternative is an index that silently holds no records, which is a
    # worse answer than "not yet".
    path = tmp_path / "ref.fasta.gz"
    path.write_bytes(gzip.compress(b">rec1\nACGT\n"))
    with pytest.raises(ValueError, match="gzip"):
        load(path)


def test_the_index_refuses_a_non_buffer():
    with pytest.raises(TypeError):
        biofasting._core.FastaIndex(42)


# --------------------------------------------------------------------------
# Differential: the corpus.
# --------------------------------------------------------------------------

GOOD_FASTA = [
    "edge/fasta/blank_lines.fasta",
    "edge/fasta/comment_lines.fasta",
    "edge/fasta/crlf.fasta",
    "edge/fasta/empty_record.fasta",
    "edge/fasta/multiline_widths.fasta",
    "edge/fasta/no_trailing_newline.fasta",
    "edge/fasta/single_long_line.fasta",
    "edge/fasta/soft_masked.fasta",
]


@needs_biopython
@needs_corpus
@pytest.mark.parametrize("relative", GOOD_FASTA)
def test_the_corpus_parses_identically(relative):
    text = (DATA / relative).read_bytes().decode("ascii")
    expected = reference(text)
    idx = index(text.encode("ascii"))
    assert [idx.title(name) for name in idx] == [title for title, _ in expected]
    assert [idx[name] for name in idx] == [
        sequence.encode("ascii") for _, sequence in expected
    ]


@needs_biopython
@needs_corpus
def test_the_genome_parses_identically_record_for_record():
    text = (DATA / "fasta" / "genome_1mb.fasta").read_bytes().decode("ascii")
    expected = reference(text)
    idx = index(text.encode("ascii"))
    assert list(idx) == [key_of(title) for title, _ in expected]
    assert [idx[name] for name in idx] == [s.encode("ascii") for _, s in expected]


@needs_biopython
@needs_corpus
def test_the_index_matches_seqio_index():
    """The keys, and every sequence, against the reference's own index."""
    from Bio import SeqIO

    path = DATA / "fasta" / "genome.fasta"
    reference_index = SeqIO.index(str(path), "fasta")
    try:
        idx = load(path)
        assert list(idx) == list(reference_index.keys())
        for name in reference_index:
            assert idx.sequence_length(name) == len(reference_index[name])
            assert idx[name] == str(reference_index[name].seq).encode("ascii")
    finally:
        # SeqIO.index leaves its handle open, and pytest rightly turns the
        # resulting unraisable exception into a failure.
        reference_index.close()


@needs_corpus
def test_our_index_rows_match_the_pyfaidx_written_fai():
    """The only check here whose answer came from a third implementation.

    ``genome.fasta.fai`` was written by pyfaidx before this kernel existed, so
    agreeing with it is evidence that does not pass through Biopython or
    through us.  The five columns are name, length, offset, line_bases and
    line_width -- which is why ``index_entry`` returns exactly those.
    """
    path = DATA / "fasta" / "genome.fasta"
    fai = path.with_suffix(path.suffix + ".fai")
    if not fai.exists():
        pytest.skip("the pyfaidx index is not present")
    idx = load(path)
    rows = [line.split("\t") for line in fai.read_text().splitlines() if line]
    assert len(rows) == len(idx)
    for name, length, offset, line_bases, line_width in rows:
        assert idx.index_entry(name) == (
            name,
            int(length),
            int(offset),
            int(line_bases),
            int(line_width),
        )


# --------------------------------------------------------------------------
# Fuzz: inputs nobody wrote by hand.
# --------------------------------------------------------------------------

# '>' decides where a record ends, ';' is the comment character that is not
# one, and the whitespace is all there because the reference treats each kind
# differently.
ALPHABET = b">;ACGT \t\r\n"


@needs_biopython
@pytest.mark.parametrize("seed", range(8))
def test_random_input_gets_the_same_records(seed):
    rng = random.Random(seed)
    for _ in range(300):
        lines = [
            bytes(rng.choice(ALPHABET) for _ in range(rng.randrange(0, 10)))
            for _ in range(rng.randrange(1, 10))
        ]
        data = b"\n".join(lines) + rng.choice([b"", b"\n"])
        text = data.decode("ascii")
        expected = reference(text)
        names = [key_of(title) for title, _ in expected]

        if any(name == "" for name in names):
            # The reference's own index cannot key a nameless record -- it
            # raises IndexError -- so there is nothing to compare it to.
            continue

        if len(set(names)) != len(names):
            with pytest.raises(ValueError, match="Duplicate key"):
                index(data)
            continue

        idx = index(data)
        assert list(idx) == names, f"seed={seed} data={data!r}"
        assert [idx[name] for name in idx] == [
            sequence.encode("ascii") for _, sequence in expected
        ], f"seed={seed} data={data!r}"


@needs_biopython
@pytest.mark.parametrize("seed", range(4))
def test_random_slices_agree_with_python(seed):
    rng = random.Random(seed)
    for _ in range(40):
        width = rng.randrange(1, 12)
        line_count = rng.randrange(1, 6)
        sequence = bytes(rng.choice(b"ACGT") for _ in range(width * line_count))
        lines = [sequence[i : i + width] for i in range(0, len(sequence), width)]
        data = b">rec\n" + b"\n".join(lines) + b"\n"
        idx = index(data)
        whole = idx["rec"]
        assert whole == sequence
        for _ in range(20):
            start = rng.randrange(0, len(whole) + 5)
            end = rng.randrange(0, len(whole) + 5)
            assert idx.sequence_slice("rec", start, end) == whole[start:end]
