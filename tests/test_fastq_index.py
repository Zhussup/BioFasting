"""Step 2.2d: the FASTQ index, against Bio.SeqIO.index.

Same two layers as every other kernel's tests, for the same reason.  The first
asserts behaviour on hand-written byte strings that need neither Biopython nor
the generated corpus, so CI -- which has neither -- still tests the index.  The
second is the real test: every FASTQ file in the corpus, indexed here and
indexed there, compared key for key and record for record.

The divergence this file exists to pin down is not a bug in either
implementation; it is a disagreement *inside* the reference.  Biopython has two
FASTQ grammars -- ``FastqGeneralIterator``, which ``SeqIO.parse`` runs, and
``FastqRandomAccess.__iter__``, which ``SeqIO.index`` runs -- and they do not
accept the same files.  ``bench/data/edge/blank_lines.fastq`` is a legal file
that the parser reads and the index refuses.  This package has one grammar and
it is the parser's, so the index here accepts that file and refuses malformed
input with the parser's message and its location.  Both halves are asserted
below, against Biopython, so that the choice stays a choice.
"""

import gc
import gzip
import importlib.util
import struct
import warnings
import zlib
from pathlib import Path

import pytest

import biofasting

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_corpus = pytest.mark.skipif(
    not (DATA / "fastq" / "reads_10k.fastq").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)
needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)


def index(data):
    """Build an index over `data`, letting a refused file raise."""
    return biofasting._core.FastqIndex(data)


def mine(path):
    """Index `path` and return (index, exception), the way the reference tests do."""
    try:
        return biofasting.open_fastq_index(path), None
    except ValueError as exc:  # exactly the type the reference raises
        return None, exc


def theirs(path):
    """Index `path` with Biopython and return (index, exception).

    The exception caught is `Exception`, not `ValueError`, because the two
    halves of the reference do not agree even on the type: a header with no word
    in it makes its index raise `IndexError` out of `split(None, 1)[0]`, and that
    the reference's index and its parser do not fail the same way is the subject
    of several tests here.

    A refusal leaks.  `SeqIO.index` opens the file and then iterates it to build
    its offsets, so an exception during that leaves the handle to the collector;
    this project turns a stray `ResourceWarning` into a test failure, so the
    collection is forced here, inside a block that ignores it, rather than left
    to fire at a random later moment.  That the reference leaks is the
    reference's business; failing someone else's test over it is not useful.

    The caller closes what it gets.
    """
    from Bio import SeqIO

    try:
        return SeqIO.index(str(path), "fastq"), None
    except Exception as exc:  # noqa: BLE001 - which exception it is, is the answer
        # The exception is handed back, so its frames have to be let go of
        # first: a traceback holds the half-built index, the index holds the open
        # handle, and the handle then outlives every `finally` in this module and
        # surfaces as an unraisable warning at interpreter shutdown.
        #
        # The suppression has to wrap the clearing, not follow it, and that is
        # not obvious: dropping the traceback is the very statement that frees
        # the buffer, so the `ResourceWarning` is raised *there*, under this
        # project's `filterwarnings = ["error"]` it is raised as an exception,
        # and an exception raised while handling one is unraisable -- which
        # pytest then reports as a failure of whatever test happens to be
        # running.  Found by exactly that failure.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResourceWarning)
            exc.__traceback__ = None
            exc.__context__ = None
            exc.__cause__ = None
            gc.collect()
        return None, exc


def fastq_files():
    """Every FASTQ file in the corpus, generated and hand-written.

    `empty_file.fastq` is included rather than filtered out: an empty file is a
    valid (empty) input for both implementations, and the file-level helper has
    to survive `mmap` refusing a zero-length file.
    """
    return sorted(
        path
        for folder in ("fastq", "edge", "malformed")
        for path in (DATA / folder).glob("*.fastq")
    )


# --------------------------------------------------------------------------
# Behaviour, without the reference.
# --------------------------------------------------------------------------


def test_reads_a_plain_record():
    idx = index(b"@read1 comment\nACGT\n+\nIIII\n")
    assert len(idx) == 1
    assert list(idx) == ["read1"]
    assert idx["read1"] == b"ACGT"
    assert idx.quality("read1") == b"IIII"
    assert idx.title("read1") == "read1 comment"


def test_the_key_is_the_first_word_and_the_title_is_the_whole_header():
    idx = index(b"@read1 a description with spaces\nACGT\n+\nIIII\n")
    assert list(idx) == ["read1"]
    assert idx.title("read1") == "read1 a description with spaces"


def test_leading_whitespace_is_skipped_before_the_key():
    # `line[1:].rstrip().split(None, 1)[0]` splits on a run of whitespace and
    # drops leading whitespace, so the key is not simply everything up to the
    # first space.
    idx = index(b"@   read2 desc\nACGT\n+\nIIII\n")
    assert list(idx) == ["read2"]
    assert idx.title("read2") == "   read2 desc"


def test_a_tab_ends_the_key_like_a_space_does():
    idx = index(b"@read1\tdesc\nACGT\n+\nIIII\n")
    assert list(idx) == ["read1"]


def test_the_two_lengths_are_one_number():
    # Not two measurements that happen to agree: the parser refuses a record
    # whose quality is not exactly as long as its sequence, so a record that
    # exists has one length, and `quality_length` cannot drift from it.
    idx = index(b"@read1\nACGTACGT\n+\nIIIIIIII\n")
    assert idx.sequence_length("read1") == 8
    assert idx.quality_length("read1") == 8
    assert len(idx.quality("read1")) == 8


def test_a_wrapped_record_is_assembled_and_still_plain_to_the_caller():
    idx = index(b"@read1 desc\nACGT\nACGT\n+\nIIII\nIIII\n")
    assert idx["read1"] == b"ACGTACGT"
    assert idx.quality("read1") == b"IIIIIIII"
    assert idx.sequence_length("read1") == 8


def test_a_blank_line_between_records_is_folded_away():
    idx = index(b"@read1\nACGT\n+\nIIII\n\n@read2\nTTTT\n+\nJJJJ\n")
    assert list(idx) == ["read1", "read2"]
    assert idx["read2"] == b"TTTT"


def test_a_plus_line_may_repeat_the_header():
    idx = index(b"@read1 desc\nACGT\n+read1 desc\nIIII\n")
    assert idx["read1"] == b"ACGT"


def test_a_quality_line_may_start_with_at():
    # The trap the whole format is famous for: what ends a record is the
    # sequence length, not the character on the next line.
    idx = index(b"@read1\nACGTACGT\n+\n@@@@@@@@\n@read2\nTTTT\n+\nJJJJ\n")
    assert list(idx) == ["read1", "read2"]
    assert idx.quality("read1") == b"@@@@@@@@"
    assert idx["read2"] == b"TTTT"


def test_crlf_ends_a_line_like_the_text_stream_does():
    idx = index(b"@read1\r\nACGT\r\n+\r\nIIII\r\n@read2\r\nTTTT\r\n+\r\nJJJJ\r\n")
    assert list(idx) == ["read1", "read2"]
    assert idx["read1"] == b"ACGT"
    assert idx.quality("read1") == b"IIII"


def test_the_last_record_may_have_no_terminator():
    idx = index(b"@read1\nACGT\n+\nIIII")
    assert idx.quality("read1") == b"IIII"


def test_an_empty_record_is_a_record():
    idx = index(b"@empty\n\n+\n\n@read2\nACGT\n+\nIIII\n")
    assert list(idx) == ["empty", "read2"]
    assert idx["empty"] == b""
    assert idx.sequence_length("empty") == 0
    assert idx.quality("empty") == b""


def test_an_empty_file_has_no_records_and_is_not_an_error():
    idx = index(b"")
    assert len(idx) == 0
    assert list(idx) == []
    assert idx.keys() == []


def test_a_missing_key_raises_keyerror_and_contains_is_false():
    idx = index(b"@read1\nACGT\n+\nIIII\n")
    assert "read1" in idx
    assert "read2" not in idx
    with pytest.raises(KeyError):
        idx["read2"]
    with pytest.raises(KeyError):
        idx.quality("read2")
    with pytest.raises(KeyError):
        idx.title("read2")


def test_a_duplicate_key_is_refused_with_the_reference_wording():
    with pytest.raises(ValueError, match=r"Duplicate key 'read1'"):
        index(b"@read1 a\nACGT\n+\nIIII\n@read1 b\nTTTT\n+\nJJJJ\n")


def test_a_duplicate_is_refused_even_when_the_records_are_shaped_differently():
    # The two records take different paths through the scanner -- the second is
    # assembled, the first is a view -- so their keys come from different places.
    # A duplicate check that only saw one of those places would miss this.
    with pytest.raises(ValueError, match=r"Duplicate key 'read1'"):
        index(b"@read1\nACGT\n+\nIIII\n@read1\nACGT\nACGT\n+\nIIII\nIIII\n")


def test_records_that_are_not_fastq_are_refused_with_the_parsers_message():
    with pytest.raises(ValueError, match="should start with '@' character"):
        index(b"ACGT\nACGT\n+\nIIII\n")


def test_a_short_quality_is_the_parsers_error_and_names_its_place():
    with pytest.raises(ValueError) as caught:
        index(b"@read1\nACGTACGT\n+\nIIII\n")
    message = str(caught.value)
    assert "Lengths of sequence and quality values differs" in message
    # The location is the part the reference cannot give, and the part that
    # matters on a malformed file of a million records.
    assert "record 0" in message
    assert "byte offset" in message


def test_the_plain_flag_says_what_a_fetch_costs():
    """The regression guard a correctness-only test cannot be.

    A record that is a view of the file is a memcpy to fetch; one the scanner
    had to assemble is a re-parse.  Both return the same bytes, so if the fast
    path disappeared every other test here would stay green while the index
    quietly stopped being an index -- the same trap `strided` guards in the FASTA
    index.  So the flag is asserted, not inferred.
    """
    idx = index(b"@plain\nACGT\n+\nIIII\n@wrapped\nACGT\nACGT\n+\nIIII\nIIII\n")
    assert idx["plain"] == b"ACGT"
    assert idx["wrapped"] == b"ACGTACGT"
    assert idx.plain("plain") is True
    assert idx.plain("wrapped") is False


def test_index_entry_locates_the_record_in_the_file():
    data = b"@read1 desc\nACGT\n+\nIIII\n@read2\nTTTTTT\n+\nJJJJJJ\n"
    idx = index(data)
    name, length, offset, end, plain = idx.index_entry("read1")
    assert name == "read1"
    assert length == 4
    assert data[offset : offset + 1] == b"@"
    assert data[end - 1 : end] == b"\n"
    assert idx.plain("read1") is plain


def test_sequence_slice_clamps_and_never_reads_past_the_record():
    idx = index(b"@read1\nACGTACGT\n+\nIIIIIIII\n")
    assert idx.sequence_slice("read1", 2, 5) == b"GTA"
    assert idx.sequence_slice("read1", 0, 8) == b"ACGTACGT"
    assert idx.sequence_slice("read1", 0, 1000) == b"ACGTACGT"
    assert idx.sequence_slice("read1", 6, 1000) == b"GT"
    assert idx.sequence_slice("read1", 5, 5) == b""
    assert idx.sequence_slice("read1", 7, 3) == b""
    assert idx.sequence_slice("read1", 1000, 2000) == b""


def test_sequence_slice_on_an_assembled_record_is_its_window_too():
    # The fallback assembles the whole record and cuts it, which is the same
    # answer by a slower road -- and the only way to tell them apart is speed.
    idx = index(b"@read1\nACGT\nACGT\n+\nIIII\nIIII\n")
    assert idx.plain("read1") is False
    assert idx.sequence_slice("read1", 2, 5) == b"GTA"


# --------------------------------------------------------------------------
# The file-level helper.
# --------------------------------------------------------------------------


def bgzf_bytes(payload, block=65536):
    """A BGZF stream, written out by hand.

    BGZF is gzip with an extra field (BC, the block size) and a 28-byte empty
    block at the end, so hand-writing it here is the point rather than a
    shortcut: it produces a stream that both `gzip.decompress` and Biopython's
    `BgzfReader` accept, which is what this index has to do too.  Built with
    `zlib` only, so that this test does not need Biopython to build its input.
    """
    out = bytearray()
    for start in range(0, len(payload), block):
        chunk = payload[start:start + block]
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = compressor.compress(chunk) + compressor.flush()
        header = (b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00\xff"
                  + struct.pack("<H", 6) + b"BC" + struct.pack("<H", 2)
                  + struct.pack("<H", len(body) + 25))
        out += header + body + struct.pack("<I", zlib.crc32(chunk) & 0xFFFFFFFF)
        out += struct.pack("<I", len(chunk))
    out += (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
            b"\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00")
    return bytes(out)


def test_a_gzipped_file_is_indexed_by_inflating_it_once(tmp_path):
    # A gzip file is mapped, inflated with libdeflate, and the inflated bytes
    # are indexed.  Everything downstream of that is the plain path, so the
    # records have to be identical to the uncompressed file's -- including the
    # `plain` flag, which says the fetch is a memcpy and not a re-parse.
    raw = b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTT\n+\n!!!!\n"
    plain = tmp_path / "reads.fastq"
    plain.write_bytes(raw)
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(raw))

    idx = biofasting.open_fastq_index(path)
    assert idx.keys() == biofasting.open_fastq_index(plain).keys()
    assert idx.sequence("read1") == b"ACGTACGT"
    assert idx.quality("read2") == b"!!!!"
    assert idx.title("read2") == "read2"
    assert idx.plain("read1") is True


def test_bgzf_is_gzip_and_is_indexed_by_the_same_code(tmp_path):
    # BGZF is what `SeqIO.index` demands and what `bgzip` writes, and it is a
    # gzip stream with extra fields -- so if the hand-written stream here did
    # not read back through `gzip.decompress` first, the rest of the test would
    # be proving the wrong thing.
    raw = b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTT\n+\n!!!!\n"
    compressed = bgzf_bytes(raw)
    assert gzip.decompress(compressed) == raw

    path = tmp_path / "reads.fastq.bgz"
    path.write_bytes(compressed)
    idx = biofasting.open_fastq_index(path)
    assert idx.keys() == ["read1", "read2"]
    assert idx.sequence("read1") == b"ACGTACGT"


def test_a_gzip_file_of_several_members_is_one_stream(tmp_path):
    # `gzip -c a b > ab.gz` and a BGZF file both concatenate members, and the
    # parser reads them as one stream -- so the index has to as well, or the
    # two would disagree about the same file.
    first = gzip.compress(b"@read1\nACGT\n+\nIIII\n")
    second = gzip.compress(b"@read2\nTTTT\n+\n!!!!\n")
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(first + second)
    idx = biofasting.open_fastq_index(path)
    assert idx.keys() == ["read1", "read2"]
    assert idx.sequence("read2") == b"TTTT"


def test_a_corrupt_gzip_stream_is_refused_not_indexed_as_text(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(b"@read1\nACGT\n+\nIIII\n")[:-4])
    with pytest.raises(ValueError):
        biofasting.open_fastq_index(path)


def test_the_index_owns_the_inflated_bytes_not_the_caller(tmp_path):
    # The inflated buffer is the index's own, which is the whole reason the
    # returned object keeps it: an index that borrowed the caller's copy would
    # be reading freed memory by the time this line runs.
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(b"@read1\nACGTACGT\n+\nIIIIIIII\n"))
    idx = biofasting.open_fastq_index(path)
    for _ in range(3):
        gc.collect()
    assert idx.sequence("read1") == b"ACGTACGT"


def test_open_fastq_index_refuses_a_file_that_is_not_gzip_and_not_fastq(tmp_path):
    # Refused by the scanner, with the parser's message and location, exactly as
    # for an uncompressed file: compression changes how the bytes are fetched,
    # not which bytes are legal.
    path = tmp_path / "garbage.fastq.gz"
    path.write_bytes(gzip.compress(b"this is not a FASTQ file\n"))
    with pytest.raises(ValueError, match="record 0"):
        biofasting.open_fastq_index(path)


def test_open_fastq_index_maps_the_file_and_survives_the_handle_closing(tmp_path):
    # The mapping is what keeps the bytes alive, not the file object: an index
    # that needed the handle would fail here, and on a 369 MB file it would
    # also have read the whole thing.
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"@read1\nACGT\n+\nIIII\n")
    idx = biofasting.open_fastq_index(path)
    assert idx["read1"] == b"ACGT"


def test_open_fastq_index_accepts_a_nonexistent_path_as_the_reference_does(tmp_path):
    with pytest.raises(FileNotFoundError):
        biofasting.open_fastq_index(tmp_path / "not-there.fastq")


# --------------------------------------------------------------------------
# The reference.
# --------------------------------------------------------------------------


@needs_biopython
def test_the_reference_refuses_the_blank_line_file_that_its_own_parser_reads():
    """The divergence, asserted from both sides so it cannot drift silently.

    `blank_lines.fastq` separates its records with a blank line.  Biopython's
    FASTQ *parser* reads it -- this package's scanner is digest-gated against
    that parser over the whole corpus -- and Biopython's FASTQ *index* raises.
    The two are different implementations of the same format inside one
    library, and this is where they part company.
    """
    path = DATA / "edge" / "blank_lines.fastq"
    if not path.exists():
        pytest.skip("benchmark corpus not generated (python bench/gen_data.py)")

    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    with path.open() as handle:
        parsed = [t for t, _, _ in FastqGeneralIterator(handle)]
    assert parsed == ["read1", "read2"]

    _, error = theirs(path)
    assert error is not None, "Biopython's index stopped refusing this file"
    assert "Problem with line" in str(error)

    idx = biofasting.open_fastq_index(path)
    assert list(idx) == parsed


@needs_biopython
@needs_corpus
def test_the_corpus_agrees_with_the_reference_key_for_key_and_record_for_record():
    """The real test: every FASTQ file in the corpus, indexed twice.

    Both implementations are asked for the same four things per record -- the
    key, the sequence, the quality and the title -- and the two lists are
    compared whole.  A file only one of them accepts is a failure here unless it
    is the one known divergence, which is asserted separately below; anything
    else failing is a bug.
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    # The only file the two are allowed to disagree about, and they disagree the
    # way `test_the_reference_refuses_the_blank_line_file_that_its_own_parser_reads`
    # spells out.  Named rather than pattern-matched, so that a second divergence
    # cannot hide behind a rule.
    known_divergence = {"blank_lines.fastq"}
    divergence_seen = set()
    compared = 0

    for path in fastq_files():
        reference, reference_error = theirs(path)
        ours, our_error = mine(path)
        try:
            if reference_error is not None or our_error is not None:
                if reference_error is not None and our_error is None:
                    assert path.name in known_divergence, (
                        f"{path.name}: this index accepted a file the reference "
                        f"refused ({reference_error!r}), and it is not the known "
                        "divergence"
                    )
                    divergence_seen.add(path.name)
                else:
                    assert our_error is not None, (
                        f"{path.name}: the reference refused it "
                        f"({reference_error!r}) and this index did not"
                    )
                    assert isinstance(our_error, ValueError)
                continue

            compared += 1
            assert list(ours) == list(reference), f"{path.name}: keys or order differ"
            assert len(ours) == len(reference), f"{path.name}: record count differs"

            # The parser's records, independently of both indexes: the index has
            # to agree with the grammar too, not only with the reference's index,
            # or a file the two indexes both mis-read would pass this test.
            with path.open() as handle:
                parsed = [
                    (title.split(None, 1)[0], sequence, quality)
                    for title, sequence, quality in FastqGeneralIterator(handle)
                ]
            assert list(ours) == [name for name, _, _ in parsed], path.name

            for name, sequence, quality in parsed:
                assert ours[name].decode() == sequence, (path.name, name)
                assert ours.quality(name).decode() == quality, (path.name, name)
                assert ours.sequence_length(name) == len(sequence), (path.name, name)
                assert ours.quality_length(name) == len(quality), (path.name, name)
                assert ours.sequence_slice(name, 0, len(sequence)) == ours[name]

            # And the reference's own view of each record, which is what a caller
            # swapping one index for the other would compare against.
            for name in reference:
                record = reference[name]
                assert ours[name].decode() == str(record.seq), (path.name, name)
                assert ours.title(name) == record.description, (path.name, name)
        finally:
            if reference is not None:
                reference.close()

    assert compared >= 8, f"only {compared} files were compared"
    assert divergence_seen == {"blank_lines.fastq"}, divergence_seen


@needs_biopython
@needs_corpus
def test_the_reference_and_this_index_agree_on_the_malformed_files_verdict():
    """Both refuse the malformed corpus, and for once the messages differ too.

    The reference's index has four error strings of its own -- "Problem with
    quality section", "Premature end of file in seq section" and so on -- which
    are not the ones its parser raises.  This package reports the parser's
    message with the location added, so the two differ in wording while agreeing
    on the verdict.  That is the documented divergence; the verdict is what is
    asserted as equal.
    """
    folder = DATA / "malformed"
    if not folder.exists():
        pytest.skip("benchmark corpus not generated (python bench/gen_data.py)")
    for path in sorted(folder.glob("*.fastq")):
        reference, reference_error = theirs(path)
        if reference is not None:
            reference.close()
        ours, our_error = mine(path)
        assert reference_error is not None, f"{path.name}: the reference accepted it"
        assert our_error is not None, f"{path.name}: this index accepted it"
        assert isinstance(our_error, ValueError)

@needs_biopython
@needs_corpus
def test_the_million_read_corpus_indexes_the_same_way_at_scale():
    path = DATA / "fastq" / "reads_1m.fastq"
    if not path.exists():
        pytest.skip("the 1M corpus is not generated (python bench/gen_data.py)")
    ours = biofasting.open_fastq_index(path)
    reference, error = theirs(path)
    assert error is None
    try:
        assert len(ours) == len(reference) == 1_000_000
        assert list(ours) == list(reference)
        # The plain path has to survive at scale, not only on an eight-line
        # string: a corpus of a million identically-shaped records is exactly
        # the case the fast path exists for.  Sampled from the index's own keys
        # rather than from names guessed here, since the corpus names are the
        # generator's business and not this test's.
        names = list(ours)
        for position in (0, 1, len(names) // 2, len(names) - 1):
            name = names[position]
            assert ours.plain(name) is True, name
            assert ours[name].decode() == str(reference[name].seq), name
            assert ours.quality_length(name) == len(reference[name].seq), name
    finally:
        reference.close()


@needs_biopython
def test_the_reference_refuses_the_plain_gzip_that_this_index_reads(tmp_path):
    """The file format the reference will not index, and why the row is legal.

    `SeqIO.index` accepts BGZF and refuses everything else, with a message that
    says so:

        ValueError: Gzipped files are not suitable for indexing, please use
                    BGZF (blocked gzip format) instead.

    Plain gzip is what a FASTQ file is when nobody has run `bgzip` on it, which
    is to say most of them.  This index reads it, so a benchmark row that puts
    our `.gz` beside the reference's `.bgz` is comparing two files -- and this
    test is the statement that the reference has no row at all on the first one
    rather than a worse one.
    """
    from Bio import SeqIO

    raw = b"@read1\nACGT\n+\nIIII\n"
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(raw))

    with pytest.raises(ValueError, match="BGZF"):
        SeqIO.index(str(path), "fastq")

    idx = biofasting.open_fastq_index(path)
    assert idx.keys() == ["read1"]
    assert idx["read1"] == b"ACGT"


@needs_biopython
@needs_corpus
def test_the_two_indexes_agree_on_a_bgzf_file(tmp_path):
    """The differential test for the compressed path, on a BGZF file both accept.

    A real corpus file, re-encoded as BGZF by hand, indexed here and by
    `SeqIO.index` -- which is the only compressed format it will index, and
    therefore the only place the two can be compared on equal terms.  The
    content is a corpus file's, so this is not a toy: it is the ten-thousand
    read file, five hundred blocks of BGZF, and every key and every record is
    compared.
    """
    from Bio import SeqIO

    source = DATA / "fastq" / "reads_10k.fastq"
    payload = source.read_bytes()
    path = tmp_path / "reads_10k.fastq.bgz"
    path.write_bytes(bgzf_bytes(payload))

    ours = biofasting.open_fastq_index(path)
    reference = SeqIO.index(str(path), "fastq")
    try:
        assert list(ours) == list(reference)
        assert len(ours) == 10_000
        names = list(ours)
        for position in (0, 1, len(names) // 2, len(names) - 1):
            name = names[position]
            assert ours[name].decode() == str(reference[name].seq), name
            assert ours.quality(name).decode() == reference[name].format("fastq").splitlines()[3], name
    finally:
        reference.close()


@needs_corpus
def test_the_compressed_corpus_indexes_the_same_way_as_the_uncompressed_one():
    """The million-read file, gzipped, against itself uncompressed.

    The reference cannot be the other half of this comparison -- it refuses the
    `.gz` -- and building a 174 MB BGZF copy inside a test would cost more than
    the test is worth.  The differential half is the line above; this half is
    the scale check, and it is the stronger one for what it is checking: the
    same million records, one file mapped and one file inflated, indexed to the
    same keys and the same bytes.  A compressed path that dropped or shifted a
    record would fail here and nowhere else.
    """
    plain = DATA / "fastq" / "reads_1m.fastq"
    compressed = DATA / "fastq" / "reads_1m.fastq.gz"
    if not (plain.exists() and compressed.exists()):
        pytest.skip("the 1M corpus is not generated (python bench/gen_data.py)")

    here = biofasting.open_fastq_index(plain)
    there = biofasting.open_fastq_index(compressed)
    assert len(here) == len(there) == 1_000_000
    assert list(here) == list(there)
    names = list(there)
    for position in (0, 1, len(names) // 2, len(names) - 1):
        name = names[position]
        assert there.plain(name) is True, name
        assert there.sequence(name) == here.sequence(name), name
        assert there.quality(name) == here.quality(name), name


@needs_biopython
def test_a_header_with_no_word_is_a_key_the_reference_would_refuse():
    """A corner this package answers and the reference does not.

    `b"   ".split(None, 1)[0]` raises `IndexError` in Biopython: a header with
    no word in it has no key.  Here the key is the empty string, which is what
    the FASTA index already does for `>   ` and is therefore the consistent
    answer rather than a second rule.  Nothing in the corpus has such a header;
    this is written down because it is a difference either way.
    """
    reference, reference_error = theirs_bytes(b"@   \nACGT\n+\nIIII\n")
    assert reference_error is not None
    assert isinstance(reference_error, IndexError)

    idx = index(b"@   \nACGT\n+\nIIII\n")
    assert list(idx) == [""]
    assert idx[""] == b"ACGT"


def theirs_bytes(data):
    """Index a byte string with Biopython, which insists on a real path."""
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".fastq", delete=False) as handle:
        handle.write(data)
        path = Path(handle.name)
    try:
        reference, error = theirs(path)
        if reference is not None:
            reference.close()
        return reference, error
    finally:
        path.unlink()
