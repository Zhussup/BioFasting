"""`biofasting.checksum`, against `Bio.SeqUtils.CheckSum`.

Four functions, and the only hard part of any of them is what a character is.
`crc32` and `seguid` are `binascii` and `hashlib` on both sides -- what is under
test is the coercion in front of them, which is ``bytes(seq)`` if that works and
``seq.encode()`` if it does not, so a `Seq` is read as its bytes and a `str` as
UTF-8.  `crc64` folds each character in with ``ord(c) & 0xFF``, which is a byte
operation on Latin-1 text and therefore never needs to fall back.  `gcg`
upper-cases first, and ``str.upper`` is not a byte operation, so it does.

The doctest values are the reference's own, including the one in `gcg`'s
docstring that is easy to misread: ``gcg("ACGTACGTACGT")`` is 5688 and
``gcg("acgtACGTacgt")`` is 5688 as well, because `gcg` is case-insensitive while
`crc32` and `crc64` next to it are not.
"""

import importlib.util
import random
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


# --------------------------------------------------------------------------
# The reference's own doctests, as hand-written expectations.
# --------------------------------------------------------------------------


def test_crc32_hand_checked_values():
    assert biofasting.crc32("ACGTACGTACGT") == 20049947
    assert biofasting.crc32("acgtACGTacgt") == 1688586483
    assert biofasting.crc32("") == 0


def test_crc64_hand_checked_values():
    assert biofasting.crc64("ACGTACGTACGT") == "CRC-C4FBB762C4A87EBD"
    assert biofasting.crc64("acgtACGTacgt") == "CRC-DA4509DC64A87EBD"
    assert biofasting.crc64("") == "CRC-0000000000000000"


def test_gcg_hand_checked_values():
    assert biofasting.gcg("ACGTACGTACGT") == 5688
    # Case-insensitive, unlike the two crc functions next to it.
    assert biofasting.gcg("acgtACGTacgt") == 5688
    assert biofasting.gcg("") == 0
    # The index counts 1..57 and then starts over, so the 58th character is
    # weighted 1 again -- and the two blocks add, which is why a doubled
    # sequence is not simply twice the checksum.
    assert biofasting.gcg("A" * 57) == sum(range(1, 58)) * ord("A") % 10000
    assert biofasting.gcg("A" * 58) == (sum(range(1, 58)) + 1) * ord("A") % 10000


def test_seguid_hand_checked_values():
    assert biofasting.seguid("ACGTACGTACGT") == "If6HIvcnRSQDVNiAoefAzySc6i4"
    assert biofasting.seguid("acgtACGTacgt") == "If6HIvcnRSQDVNiAoefAzySc6i4"
    assert biofasting.seguid("") == "2jmj7l5rSw0yVb/vlWAYkK/YBwk"


# --------------------------------------------------------------------------
# What a character is.
# --------------------------------------------------------------------------


@needs_biopython
def test_a_seq_is_read_as_its_bytes_and_a_str_as_utf8():
    """The one coercions in the module, and the one place they differ.  A
    `Seq` has ``__bytes__``, so ``bytes(seq)`` succeeds and the string is never
    encoded; a `str` does not, so it is encoded as UTF-8.  For ASCII the two
    agree, which is why the difference is easy to miss."""
    from Bio.Seq import Seq

    assert biofasting.crc32("ACGT") == biofasting.crc32(Seq("ACGT"))
    assert biofasting.crc32(Seq(b"ACGT")) == biofasting.crc32("ACGT")
    # Two bytes, because "é".encode() is UTF-8 and not Latin-1.
    assert biofasting.crc32("ACGTé") == biofasting.crc32(b"ACGT\xc3\xa9")
    assert biofasting.seguid("ACGTé") == biofasting.seguid(b"ACGT\xc3\xa9")


def test_crc64_needs_no_fallback_for_latin1_text():
    """It folds each character in with ``ord(c) & 0xFF``, so for Latin-1 the
    character *is* the byte.  Only a character above U+00FF, which has no
    one-byte encoding, goes back to Python."""
    assert biofasting.crc64("ACGTé") == biofasting.crc64(b"ACGT\xe9")
    assert biofasting.crc64("ACGTé") != biofasting.crc64(b"ACGT\xc3\xa9")


def test_gcg_falls_back_to_python_and_keeps_the_reference_failure():
    """`"ß".upper()` is two characters, so ``ord`` refuses it.  Reproducing
    that in C++ would mean reproducing Unicode's case tables; declining and
    running the reference's own loop reproduces the exception instead."""
    with pytest.raises(TypeError, match="ord\\(\\) expected a character"):
        biofasting.gcg("ACGTß")
    with pytest.raises(TypeError, match="ord\\(\\) expected a character"):
        biofasting.gcg("ß")
    # A single-character upper-case mapping is a Latin-1 byte the kernel will
    # not guess at, and Python answers it: "é".upper() is "É", 233 -> 201.
    assert biofasting.gcg("é") == 201


def test_gcg_is_case_insensitive_and_crc64_is_not():
    assert biofasting.gcg("acgt") == biofasting.gcg("ACGT")
    assert biofasting.crc64("acgt") != biofasting.crc64("ACGT")
    assert biofasting.crc32("acgt") != biofasting.crc32("ACGT")


# --------------------------------------------------------------------------
# The differential layer.
# --------------------------------------------------------------------------


@needs_biopython
def test_every_checksum_matches_biopython_over_the_whole_byte_range():
    """Every byte value, at every offset, over 3,000 random sequences -- and
    every exception compared, not just every answer."""
    from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid

    rng = random.Random(41)
    for _ in range(3000):
        data = bytes(rng.randrange(256) for _ in range(rng.randrange(0, 90)))
        text = data.decode("latin-1")
        for ours, reference in (
            (biofasting.crc32, crc32),
            (biofasting.crc64, crc64),
            (biofasting.gcg, gcg),
            (biofasting.seguid, seguid),
        ):
            try:
                expected = reference(text)
            except Exception as error:  # noqa: BLE001 - the message is the point
                with pytest.raises(type(error), match=None):
                    ours(text)
                continue
            assert ours(text) == expected, (ours.__name__, data)


@needs_biopython
def test_the_table_the_kernel_derives_is_the_reference_table():
    """The crc64 high table is 256 entries the kernel computes at compile time
    from the reference's own recurrence rather than being handed.  Derived is
    only acceptable if it can be checked, so this checks every entry."""
    from Bio.SeqUtils.CheckSum import _table_h

    table = biofasting._core.crc64_table_h()
    assert len(table) == 256
    assert list(table) == _table_h


@needs_biopython
@needs_corpus
def test_the_checksums_match_biopython_on_real_reads():
    from Bio.Seq import Seq
    from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid

    for _, sequence, _ in biofasting.open_fastq(DATA / "fastq" / "reads_10k.fastq"):
        assert biofasting.crc32(sequence) == crc32(Seq(sequence)), sequence
        assert biofasting.crc64(sequence) == crc64(Seq(sequence)), sequence
        assert biofasting.gcg(sequence) == gcg(Seq(sequence)), sequence
        assert biofasting.seguid(sequence) == seguid(Seq(sequence)), sequence
        assert biofasting.crc64(sequence.decode()) == crc64(sequence.decode())
