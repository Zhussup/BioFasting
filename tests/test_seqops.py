"""Step 1.3: the sequence kernels, against Bio.Seq and Bio.SeqUtils.

Same two layers as the readers, for the same reason.  The first layer asserts
the behaviour on strings written by hand, against the specification written out
here in Python -- neither Biopython nor the corpus is needed, so CI, which has
neither, still tests the kernels.  The second is the real evidence: the same
operations through `Bio.Seq` and `Bio.SeqUtils` over the whole byte range, over
randomised sequences and over the real genome.

The exhaustive sweeps are not decoration.  Both kernels here have a vector path
that only runs on some inputs and a scalar tail that catches the rest, and the
seam between them is where the bugs are: the reverse-complement kernel shipped
with a bound that skipped the input tail whenever the length was near a
multiple of 64, which no hand-written case and no test shorter than 64 bases
would have caught.  `test_every_length_agrees_with_the_specification` exists
because of it.
"""

import importlib.util
import itertools
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

import biofasting

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_corpus = pytest.mark.skipif(
    not (DATA / "fasta" / "genome_1mb.fasta").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)
needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)


# --------------------------------------------------------------------------
# The specification, written out.  These are transcriptions of the reference's
# own definitions rather than re-implementations of the kernels: if one of them
# is wrong, the differential tests below will disagree with Biopython and say
# so.
# --------------------------------------------------------------------------

# Bio.Seq._dna_complement_table, as a text translation table.
_COMPLEMENT = str.maketrans(
    "ABCDGHKMRTUVYabcdghkmrtuvy", "TVGHCDMKYAABRtvghcdmkyaabr"
)


def spec_reverse_complement(text):
    """`Seq.reverse_complement`: translate through the table, then reverse."""
    return text.translate(_COMPLEMENT)[::-1]


def spec_gc_fraction(text):
    """`SeqUtils.gc_fraction` with its default ambiguous="remove"."""
    gc = sum(text.count(base) for base in "CGScgs")
    length = gc + sum(text.count(base) for base in "ATWUatwu")
    if length == 0:
        return 0
    return gc / length


# --------------------------------------------------------------------------
# Reverse complement.
# --------------------------------------------------------------------------


def test_reverse_complement_of_an_empty_sequence_is_empty():
    assert biofasting.reverse_complement("") == ""
    assert biofasting.reverse_complement(b"") == b""


def test_a_palindrome_comes_back_unchanged():
    """EcoRI's site is its own reverse complement; a kernel that forgot to
    reverse, or reversed twice, cannot pass this one."""
    assert biofasting.reverse_complement("GAATTC") == "GAATTC"


def test_the_iupac_pairs_are_swapped_in_both_cases():
    for base, complement in [
        ("A", "T"),
        ("C", "G"),
        ("G", "C"),
        ("T", "A"),
        ("R", "Y"),
        ("Y", "R"),
        ("S", "S"),
        ("W", "W"),
        ("K", "M"),
        ("M", "K"),
        ("B", "V"),
        ("V", "B"),
        ("D", "H"),
        ("H", "D"),
        ("N", "N"),
        ("U", "A"),
    ]:
        assert biofasting.reverse_complement(base) == complement
        assert biofasting.reverse_complement(base.lower()) == complement.lower()


def test_a_byte_that_is_not_a_nucleotide_is_left_alone():
    """The table is the identity everywhere it is not names, so a gap or a
    wildcard survives -- and it survives *in place after the reversal*, which
    is the part a hand-written case would get wrong."""
    assert biofasting.reverse_complement("AC-GT*") == "*AC-GT"


@pytest.mark.parametrize(
    "length",
    [0, 1, 31, 32, 33, 63, 64, 65, 95, 96, 97, 127, 128, 129, 511, 512, 513],
)
def test_every_length_agrees_with_the_specification(length):
    """The vector kernel consumes 32 bytes at a time and leaves the tail to the
    scalar one; these are the lengths where the two meet, including 64 -- the
    length at which the original bound left the first 32 bytes of the result
    unwritten."""
    rng = random.Random(length)
    text = "".join(rng.choice("ACGTN") for _ in range(length))
    assert biofasting.reverse_complement(text) == spec_reverse_complement(text)


def test_every_byte_agrees_with_the_specification():
    """All 256 bytes in one sequence, so the vector path sees every byte it is
    allowed to reject and the scalar path sees every byte the table names."""
    text = "".join(chr(byte) for byte in range(256))
    assert biofasting.reverse_complement(text) == spec_reverse_complement(text)


def test_str_in_gives_str_out_and_bytes_in_gives_bytes_out():
    assert isinstance(biofasting.reverse_complement("ACGT"), str)
    assert isinstance(biofasting.reverse_complement(b"ACGT"), bytes)
    assert isinstance(biofasting.reverse_complement(bytearray(b"ACGT")), bytes)
    assert isinstance(biofasting.reverse_complement(memoryview(b"ACGT")), bytes)


def test_a_character_that_has_no_byte_is_rejected():
    """The kernel is defined on bytes, so a string it cannot represent is
    refused rather than quietly encoded into something else."""
    with pytest.raises(UnicodeEncodeError):
        biofasting.reverse_complement("ACGTα")


# --------------------------------------------------------------------------
# GC fraction.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("", 0.0),
        ("G", 1.0),
        ("A", 0.0),
        ("GCGC", 1.0),
        ("ATAT", 0.0),
        ("ACTG", 0.5),
        ("ggcc", 1.0),
        ("ACTGSSSS", 0.75),
        ("GGAUCUUCGGAUCU", 0.5),
        ("ACTGN", 0.5),
        ("N", 0.0),
        ("NNNN", 0.0),
        ("-----", 0.0),
        ("ACTG-", 0.5),
    ],
)
def test_gc_fraction_hand_checked_values(text, expected):
    assert biofasting.gc_fraction(text) == pytest.approx(expected)


def test_an_ambiguous_base_is_removed_from_the_denominator_too():
    """N is not "not GC" -- it is not in the sequence at all, for this
    calculation.  GC of "ACTGN" is 0.5, not 0.4."""
    assert biofasting.gc_fraction("ACTGN") == pytest.approx(0.5)
    assert biofasting.gc_fraction("ACTGBDHKMNRVXY") == pytest.approx(0.5)


def test_every_byte_agrees_with_the_specification():
    text = "".join(chr(byte) for byte in range(256))
    assert biofasting.gc_fraction(text) == spec_gc_fraction(text)


# --------------------------------------------------------------------------
# k-mer counting.
# --------------------------------------------------------------------------


def test_count_kmers_counts_overlapping_windows():
    """Every window, not every `k`-th one: "ACGTACGT" has five 4-mers, and the
    two that are ACGT sit one whole turn of the repeat apart."""
    assert biofasting.count_kmers("AAAA", 2) == {"AA": 3}
    assert biofasting.count_kmers("ACGTACGT", 4) == {
        "ACGT": 2,
        "CGTA": 1,
        "GTAC": 1,
        "TACG": 1,
    }
    assert biofasting.count_kmers("ACGT", 1) == {"A": 1, "C": 1, "G": 1, "T": 1}


def test_a_window_containing_anything_else_is_not_a_kmer():
    """The counts have to be comparable with `Seq.count_overlap` per key, and
    that counts an exact substring -- so no window may span an N.  The run
    restarts at the N, which is why the good windows on either side of it still
    count: "AC-GT" has two 2-mers, not zero."""
    assert biofasting.count_kmers("ACGTNACGT", 4) == {"ACGT": 2}
    assert biofasting.count_kmers("NNNN", 1) == {}
    assert biofasting.count_kmers("AC-GT", 2) == {"AC": 1, "GT": 1}
    assert biofasting.count_kmers("A-N", 2) == {}


def test_a_run_shorter_than_k_produces_nothing():
    assert biofasting.count_kmers("ACGT", 5) == {}
    assert biofasting.count_kmers("", 1) == {}


def test_counting_is_case_insensitive():
    assert biofasting.count_kmers("acgtACGT", 4) == {
        "ACGT": 2,
        "CGTA": 1,
        "GTAC": 1,
        "TACG": 1,
    }


def test_the_result_is_sorted_by_kmer():
    counts = biofasting.count_kmers("TTTTAAAACCCCGGGG", 2)
    assert list(counts) == sorted(counts)


@pytest.mark.parametrize("k", [0, 17, 100])
def test_a_k_outside_the_supported_range_is_rejected(k):
    with pytest.raises(ValueError, match="between 1 and 16"):
        biofasting.count_kmers("ACGT", k)


# --------------------------------------------------------------------------
# Dispatch.
# --------------------------------------------------------------------------


def test_the_selected_rung_is_one_the_build_defines():
    assert biofasting.seqops_level() in biofasting.cpu_levels()


def test_the_override_pins_the_scalar_path_and_changes_no_answer():
    """`BIOFASTING_SEQOPS_LEVEL=baseline` is what makes the AVX2 number in the
    report reproducible, so it is tested rather than trusted -- in a subprocess,
    because the choice is cached once per process by design."""
    program = (
        "import biofasting as b;"
        "print(b.seqops_level(),"
        "      b.reverse_complement('ACGTN' * 40),"
        "      b.gc_fraction('ACGTN' * 40),"
        "      b.count_kmers('ACGTTGCA' * 4, 3))"
    )
    def run(**environment):
        result = subprocess.run(
            [sys.executable, "-c", program],
            capture_output=True,
            text=True,
            check=True,
            env={**os.environ, **environment},
        )
        # The editable install rebuilds the extension on import and says so on
        # stdout, so the program's own line is the last one and not the first.
        return result.stdout.strip().splitlines()[-1]

    auto = run()
    pinned = run(BIOFASTING_SEQOPS_LEVEL="baseline")
    assert auto.split()[0] in ("baseline", "avx2")
    assert pinned.startswith("baseline ")
    assert pinned.split(" ", 1)[1] == auto.split(" ", 1)[1]


# --------------------------------------------------------------------------
# Differential, against Biopython.
# --------------------------------------------------------------------------


@needs_biopython
def test_reverse_complement_matches_biopython_for_every_single_byte():
    from Bio.Seq import Seq

    for byte in range(256):
        data = bytes([byte])
        assert biofasting.reverse_complement(data) == Seq(data).reverse_complement()._data


@needs_biopython
def test_reverse_complement_matches_biopython_on_random_bytes():
    from Bio.Seq import Seq

    rng = random.Random(11)
    for _ in range(2000):
        data = bytes(rng.randrange(256) for _ in range(rng.randrange(0, 200)))
        assert (
            biofasting.reverse_complement(data)
            == Seq(data).reverse_complement()._data
        ), data


@needs_biopython
def test_gc_fraction_matches_biopython_on_random_sequences():
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction as reference

    rng = random.Random(12)
    alphabets = [
        b"ACGT",
        b"ACGTN",
        b"ACGTWSUNBDHKMNRVXY",
        b"acgtwsunbdhkmnrvxy",
        b"ACGTacgt-*.",
        bytes(range(256)),
    ]
    for alphabet in alphabets:
        for _ in range(300):
            data = bytes(rng.choice(alphabet) for _ in range(rng.randrange(0, 120)))
            assert biofasting.gc_fraction(data) == reference(Seq(data)), data


@needs_biopython
@pytest.mark.parametrize("k", [1, 2, 3, 4, 5, 6])
def test_count_kmers_matches_count_overlap(k):
    """Per k-mer, against `Seq.count_overlap` -- which is case-sensitive, so
    the input is upper-cased first and the comparison is per key."""
    from Bio.Seq import Seq

    every_kmer = ["".join(letters) for letters in itertools.product("ACGT", repeat=k)]
    rng = random.Random(13 + k)
    for _ in range(30):
        data = bytes(rng.choice(b"ACGTN") for _ in range(rng.randrange(0, 80)))
        sequence = Seq(data.decode("ascii").upper())
        ours = biofasting.count_kmers(data, k)
        # Compared in both directions: a k-mer we report must match the
        # reference's count, and a k-mer we omit must be one the reference
        # cannot find either.  Only the second half would catch a window
        # silently dropped at a run boundary.
        assert set(ours) <= set(every_kmer)
        for kmer in every_kmer:
            assert ours.get(kmer, 0) == sequence.count_overlap(kmer), (data, kmer)


@needs_biopython
@needs_corpus
def test_the_genome_ops_match_biopython():
    """The same three operations over a real chromosome, not a toy string."""
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction as reference_gc

    index = biofasting.open_fasta(DATA / "fasta" / "genome_1mb.fasta")
    for name in index:
        sequence = index[name]
        reference_sequence = Seq(sequence)
        assert (
            biofasting.reverse_complement(sequence)
            == reference_sequence.reverse_complement()._data
        ), name
        assert biofasting.gc_fraction(sequence) == reference_gc(reference_sequence), name
        ours = biofasting.count_kmers(sequence, 6)
        for kmer, count in ours.items():
            assert reference_sequence.count_overlap(kmer) == count, (name, kmer)
