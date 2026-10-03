"""`Bio.SeqUtils`: the functions that measure a sequence, against Biopython.

The same two layers as everywhere else in this project.  The first asserts
behaviour on sequences written by hand, against a specification written out in
Python here, so that a CI without Biopython still tests the kernels.  The second
is the evidence: the same calls through `Bio.SeqUtils` over the whole byte
range, over random sequences from six alphabets, and over a real chromosome.

Three tests in here are not repeats of each other and are worth naming.  The
`GC_skew` sweep exists because the kernel shipped with `g - c` on two unsigned
64-bit counters, which wraps whenever a window holds more C than G -- 28% of the
windows in this sweep -- and answered a huge positive number instead of a
negative fraction; replayed against the shipped formula, this sweep is wrong on
**884 of its 1,200 draws**, and the hand-written cases did not catch it.
`test_molecular_weight_is_exactly_the_reference_sum` exists because
CPython 3.12 and later sum floats with Neumaier compensation and a naive
running total disagrees with the reference on about a quarter of all inputs.
And the `int`/`float` test on the empty sequence exists because ``0 == 0.0``,
so nothing else in the file would notice.
"""

import importlib.util
import random
from pathlib import Path

import pytest

import biofasting
from biofasting import _iupac_data

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
# The specification, written out.
# --------------------------------------------------------------------------


def spec_gc123(text):
    """`SeqUtils.GC123`'s per-codon loop, with the padding and the two zeros."""
    counts = {nt: [0, 0, 0] for nt in "ATGC"}
    for i in range(0, len(text), 3):
        codon = text[i : i + 3]
        if len(codon) < 3:
            codon += "  "
        for pos in range(3):
            for nt in "ATGC":
                if codon[pos] == nt or codon[pos] == nt.lower():
                    counts[nt][pos] += 1
    gc = {}
    gcall = nall = 0
    for i in range(3):
        n = sum(counts[nt][i] for nt in "ATGC")
        gc[i] = (counts["G"][i] + counts["C"][i]) * 100.0 / n if n else 0
        gcall += counts["G"][i] + counts["C"][i]
        nall += n
    return 100.0 * gcall / nall, gc[0], gc[1], gc[2]


def spec_gc_skew(text, window=100):
    """`SeqUtils.GC_skew`, zero-window arm included."""
    values = []
    for i in range(0, len(text), window):
        s = text[i : i + window]
        g = s.count("G") + s.count("g")
        c = s.count("C") + s.count("c")
        values.append((g - c) / (g + c) if g + c else 0.0)
    return values


# --------------------------------------------------------------------------
# Hand-written behaviour, no Biopython required.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("ACTGTN", (40.0, 50.0, 50.0, 0.0)),
        ("", None),  # nothing counted at all: the reference raises
        ("NNNN", None),
        ("A", (0.0, 0.0, 0, 0)),
        ("ACG", (66.66666666666667, 0.0, 100.0, 100.0)),
        ("G", (100.0, 100.0, 0, 0)),
        # A C G in frame 0, then the trailing T alone in frame 0 again: the
        # reference pads the partial codon with spaces, so the fourth base is
        # counted and the two spaces are not.
        ("ACGT", (50.0, 0.0, 100.0, 100.0)),
    ],
)
def test_gc123_hand_checked_values(text, expected):
    if expected is None:
        with pytest.raises(ZeroDivisionError):
            biofasting.GC123(text)
    else:
        assert biofasting.GC123(text) == expected


def test_gc123_a_frame_with_nothing_counted_is_the_integer_zero():
    """``0`` and not ``0.0``, and only where the reference's division actually
    raises.  "A" puts one base in frame 0 -- so frame 0 divides and answers
    ``0.0`` -- and nothing at all in frames 1 and 2, which answer the integer
    ``0`` through the ``except`` arm.  Nothing but ``type()`` can tell those
    apart, so nothing but ``type()`` can test this."""
    gcall, first, second, third = biofasting.GC123("A")
    assert first == 0.0 and type(first) is float
    assert second == 0 and type(second) is int
    assert third == 0 and type(third) is int
    assert type(gcall) is float


def test_gc123_leaves_the_padding_out_of_the_denominator():
    """A base the reference cannot count is not in the denominator either, and
    that includes the spaces it pads a partial trailing codon with -- so adding
    one changes nothing at all."""
    assert biofasting.GC123("ACG") == biofasting.GC123("ACG ")
    assert biofasting.GC123("ACG") == biofasting.GC123("ACG\n")
    # Five bases, five counted: the A that lands in frame 1 counts, even though
    # it is in the padded codon.
    assert biofasting.GC123("ACGT")[0] == 50.0
    assert biofasting.GC123("ACGTA")[0] == 40.0
    # The N in the last codon is counted by nobody, so the denominator is 8 and
    # not 9, and the four GC bases that survive give 50%.
    assert biofasting.GC123("ACGTACGTN") == (50.0, 33.333333333333336,
                                             33.333333333333336, 100.0)


@pytest.mark.parametrize(
    ("text", "window", "expected"),
    [
        ("GGGCCC", 3, [1.0, -1.0]),
        ("CCCGGG", 3, [-1.0, 1.0]),
        ("AAAA", 3, [0.0, 0.0]),
        ("", 100, []),
        ("ACGT", 0, None),
        ("ACGT", -1, []),
        ("ACGT", 3.0, None),
    ],
)
def test_gc_skew_hand_checked_values(text, window, expected):
    if expected is None:
        with pytest.raises((ValueError, TypeError)):
            biofasting.GC_skew(text, window)
    else:
        assert biofasting.GC_skew(text, window) == expected


def test_gc_skew_is_negative_when_there_is_more_c_than_g():
    """The bug this kernel shipped with: ``g - c`` on two unsigned counters
    wraps to ``2**64 - d``, so every window with more C than G came out as a huge
    positive number instead of a negative fraction -- ``8.8e17`` in the sweep's
    first failing draw, and ``2**63 - 1`` for the smallest case, a window of
    ``CC``."""
    assert biofasting.GC_skew("CCCCGGGG", 4) == [-1.0, 1.0]
    assert biofasting.GC_skew("CCCAAA", 3) == [-1.0, 0.0]


def test_gc_skew_window_zero_raises_the_range_error():
    with pytest.raises(ValueError, match="range\\(\\) arg 3 must not be zero"):
        biofasting.GC_skew("ACGT", 0)


@pytest.mark.parametrize(
    ("text", "mode", "expected"),
    [
        ("", "remove", 0),
        ("", "ignore", 0),
        ("", "weighted", 0),
        ("ACGT", "remove", 0.5),
        ("ACTGN", "remove", 0.5),
        ("ACTGN", "ignore", 0.4),
        ("ACTGN", "weighted", 0.5),
        ("GDVV", "ignore", 0.25),
        ("GDVV", "remove", 1.0),
        ("ACTGSSSS", "remove", 0.75),
        ("ACTGSSSS", "ignore", 0.75),
        ("ACTGSSSS", "weighted", 0.75),
    ],
)
def test_gc_fraction_modes_hand_checked_values(text, mode, expected):
    result = biofasting.gc_fraction(text, mode)
    assert result == expected
    assert type(result) is type(expected)


def test_gc_fraction_returns_the_integer_zero_when_nothing_was_counted():
    """``0`` for an empty sequence and ``0.0`` for a real fraction of zero.
    Nothing but ``type()`` can tell these apart, which is the whole reason the
    kernel returns a sentinel instead of a fraction."""
    assert biofasting.gc_fraction("") == 0
    assert type(biofasting.gc_fraction("")) is int
    assert biofasting.gc_fraction("AAAA") == 0.0
    assert type(biofasting.gc_fraction("AAAA")) is float


def test_gc_fraction_rejects_an_unknown_mode():
    with pytest.raises(ValueError, match="ambiguous value 'both' not recognized"):
        biofasting.gc_fraction("ACGT", "both")


# --------------------------------------------------------------------------
# Molecular weight.
# --------------------------------------------------------------------------


def test_molecular_weight_hand_checked_values():
    assert round(biofasting.molecular_weight("AGC"), 2) == 949.61
    assert round(biofasting.molecular_weight("AGC", "RNA"), 2) == 997.61
    assert round(biofasting.molecular_weight("AGC", "protein"), 2) == 249.29
    # Nothing to join, so the mass is one water molecule added rather than
    # subtracted: (len - 1) is -1.
    assert biofasting.molecular_weight("") == 18.0153
    assert biofasting.molecular_weight("", monoisotopic=True) == 18.010565


def test_molecular_weight_names_the_letter_it_refuses():
    """The message has doubled quotes around it, because the reference
    interpolates a `KeyError` -- and `str(KeyError("N"))` is ``"'N'"``.  It
    looks like a typo.  It is the reference's output."""
    with pytest.raises(ValueError) as caught:
        biofasting.molecular_weight("ACGTN")
    assert str(caught.value) == "''N'' is not a valid unambiguous letter for DNA"
    with pytest.raises(ValueError) as caught:
        biofasting.molecular_weight("ACGTACGT", "RNA")
    assert str(caught.value) == "''T'' is not a valid unambiguous letter for RNA"
    with pytest.raises(ValueError) as caught:
        biofasting.molecular_weight("ACGT", "cDNA")
    assert str(caught.value) == "Allowed seq_types are DNA, RNA or protein, not 'cDNA'"


def test_molecular_weight_reports_the_first_bad_letter_in_sequence_order():
    with pytest.raises(ValueError, match="''Z''"):
        biofasting.molecular_weight("ZZN")
    with pytest.raises(ValueError, match="''N''"):
        biofasting.molecular_weight("ACGTNZ")


def test_molecular_weight_strips_whitespace_before_counting():
    assert biofasting.molecular_weight("A G\nC\tT") == biofasting.molecular_weight("AGCT")


def test_molecular_weight_refuses_double_stranded_protein():
    with pytest.raises(ValueError, match="protein sequences cannot be double-stranded"):
        biofasting.molecular_weight("ACDEF", "protein", double_stranded=True)


@needs_biopython
def test_molecular_weight_is_exactly_the_reference_sum():
    """The reason the kernel compensates.  CPython 3.12 and later sum floats
    with Neumaier compensation, and the reference is written ``sum(...)``; over
    these 20,000 draws a naive running total disagrees with it about a quarter
    of the time.  Anything less than ``==`` here would be hiding that."""
    from Bio.SeqUtils import molecular_weight as reference

    rng = random.Random(21)
    letters = "ACGT"
    mismatches = 0
    naive = 0.0
    for _ in range(20000):
        text = "".join(rng.choice(letters) for _ in range(rng.randrange(1, 13)))
        expected = reference(text)
        assert biofasting.molecular_weight(text) == expected, text
        for letter in text:
            naive += _iupac_data.DNA_WEIGHTS[letter]
        naive -= (len(text) - 1) * 18.0153
        if naive != expected:
            mismatches += 1
        naive = 0.0
    assert mismatches > 1000, (
        "a naive sum agreed with the reference too often for this test to be "
        "evidence; the interpreter's summation has changed"
    )


# --------------------------------------------------------------------------
# The protein alphabets and the degenerate search.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("MAIVMGRWKGAR*", "MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer"),
        ("MAIVMGRWKGA--R*", "MetAlaIleValMetGlyArgTrpLysGlyAlaXaaXaaArgTer"),
        ("", ""),
        ("M", "Met"),
    ],
)
def test_seq3_hand_checked_values(text, expected):
    assert biofasting.seq3(text) == expected


def test_seq3_undef_code_and_custom_map():
    assert (
        biofasting.seq3("MAIVMGRWKGA--R*", undef_code="---")
        == "MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer"
    )
    assert (
        biofasting.seq3("MAIVMGRWKGAR*", custom_map={"*": "***"})
        == "MetAlaIleValMetGlyArgTrpLysGlyAlaArg***"
    )


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer", "MAIVMGRWKGAR*"),
        ("METalaIlEValMetGLYArgtRplysGlyAlaARGTer", "MAIVMGRWKGAR*"),
        ("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer", "MAIVMGRWKGAXXR*"),
        ("", ""),
        ("Ala", "A"),
        # Three characters at a time: the trailing pair is not a codon and is
        # not reported as an unknown one either.
        ("AlaAr", "A"),
    ],
)
def test_seq1_hand_checked_values(text, expected):
    assert biofasting.seq1(text) == expected


def test_seq1_undef_code_and_custom_map():
    assert (
        biofasting.seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer", undef_code="?")
        == "MAIVMGRWKGA??R*"
    )
    assert (
        biofasting.seq1("MetAlaIleValMetGlyArgTrpLysGlyAla***", custom_map={"***": "*"})
        == "MAIVMGRWKGA*"
    )


@pytest.mark.parametrize(
    ("seq", "subseq", "expected"),
    [
        ("ACGTACGT", "ACGT", ["ACGT", 0, 4]),
        ("AAAA", "NN", ["[GATC][GATC]", 0, 1, 2]),
        ("AAAA", "ACGTN", ["ACGT[GATC]"]),
        ("AAAA", "C", ["C"]),
        ("ACGTACGT", "NN", ["[GATC][GATC]", 0, 1, 2, 3, 4, 5, 6]),
        ("", "A", ["A"]),
    ],
)
def test_nt_search_hand_checked_values(seq, subseq, expected):
    assert biofasting.nt_search(seq, subseq) == expected


def test_nt_search_builds_the_same_pattern_as_the_reference():
    """The pattern is the first element of the result, so the ambiguity codes
    are checkable without Biopython: every letter of the IUPAC set expands to
    its own alternatives, in the reference's order."""
    for letter, value in _iupac_data.AMBIGUOUS_DNA_VALUES.items():
        pattern = biofasting.nt_search("", letter)[0]
        assert pattern == (value if len(value) == 1 else f"[{value}]"), letter


# --------------------------------------------------------------------------
# The six-frame picture and the codon adaptation index.
# --------------------------------------------------------------------------


def test_six_frame_translations_reproduces_the_reference_doctest():
    expected = (
        "GC_Frame: a:5 t:0 g:8 c:5\n"
        "Sequence: auggccauug ... gggccgcuga, 24 nt, 54.17 %GC\n"
        "\n"
        "\n"
        "1/1\n"
        "  G  H  C  N  G  P  L\n"
        " W  P  L  *  W  A  A\n"
        "M  A  I  V  M  G  R  *\n"
        "auggccauuguaaugggccgcuga   54 %\n"
        "uaccgguaacauuacccggcgacu\n"
        "A  M  T  I  P  R  Q\n"
        " H  G  N  Y  H  A  A  S\n"
        "  P  W  Q  L  P  G  S\n"
        "\n"
    )
    assert biofasting.six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA") == expected


def test_six_frame_translations_of_nothing_is_just_the_header():
    assert biofasting.six_frame_translations("") == (
        "GC_Frame: a:0 t:0 g:0 c:0\nSequence: , 0 nt, 0.00 %GC\n\n\n"
    )


def test_codon_adaptation_index_of_preferred_codons_only_is_one():
    index = biofasting.CodonAdaptationIndex(["GCTGCTGCTGCT"])
    assert index.calculate("GCTGCTGCTGCT") == pytest.approx(1.0)
    assert str(index).endswith("\n")
    assert len(index) == 64


def test_codon_adaptation_index_rejects_an_illegal_codon():
    with pytest.raises(ValueError, match="illegal codon 'NNN'"):
        biofasting.CodonAdaptationIndex(["ACGTTTNNN"])


@needs_biopython
def test_codon_adaptation_index_reports_the_gene_name_of_a_seqrecord():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ACGTTTNNN"), id="gene7")
    with pytest.raises(ValueError, match="illegal codon 'NNN' in gene gene7"):
        biofasting.CodonAdaptationIndex([record])


def test_codon_adaptation_index_calculate_rejects_an_illegal_codon():
    index = biofasting.CodonAdaptationIndex(["GCTGCTGCTGCT"])
    with pytest.raises(TypeError, match="illegal codon in sequence: NNN"):
        index.calculate("GCTNNN")


# --------------------------------------------------------------------------
# The differential layer.
# --------------------------------------------------------------------------


@needs_biopython
def test_everything_matches_biopython_over_six_alphabets():
    from Bio.SeqUtils import GC123, GC_skew, gc_fraction, molecular_weight

    rng = random.Random(31)
    alphabets = [
        "ACGT",
        "ACGTacgt",
        "ACGTN",
        "ACGTWSUNBDHKMNRVXYacgtwsunbdhkmnrvxy",
        "ACGTacgt-*.",
        "ACDEFGHIKLMNPQRSTVWY*-",
    ]
    for alphabet in alphabets:
        for _ in range(200):
            text = "".join(rng.choice(alphabet) for _ in range(rng.randrange(0, 120)))
            window = rng.randrange(1, 40)
            ours = biofasting.GC_skew(text, window)
            theirs = GC_skew(text, window)
            assert ours == theirs, (text, window)
            for mode in ("remove", "ignore", "weighted"):
                a = biofasting.gc_fraction(text, mode)
                b = gc_fraction(text, mode)
                assert a == b and type(a) is type(b), (text, mode)
            if not alphabet.startswith("ACGT") or "N" in alphabet:
                continue
            try:
                expected = GC123(text)
            except ZeroDivisionError:
                # A sequence with nothing countable at all.  The reference
                # raises from its uncaught whole-sequence division, and so does
                # this -- which is the assertion, not a skip.
                with pytest.raises(ZeroDivisionError):
                    biofasting.GC123(text)
            else:
                assert biofasting.GC123(text) == expected, text
            for seq_type, letters in (
                ("DNA", "ACGT"),
                ("RNA", "ACGU"),
                ("protein", "ACDEFGHIKLMNPQRSTVWY"),
            ):
                if not text or not set(text.upper()) <= set(letters):
                    continue
                assert biofasting.molecular_weight(text, seq_type) == (
                    molecular_weight(text, seq_type)
                ), (text, seq_type)


@needs_biopython
def test_GC123_matches_biopython_over_the_whole_byte_range():
    """Every byte a sequence can hold, at every position of a codon, because
    the kernel decides which of them count with a 256-entry table."""
    from Bio.SeqUtils import GC123 as reference

    rng = random.Random(32)
    alphabet = bytes(range(1, 256))
    for _ in range(400):
        data = bytes(rng.choice(alphabet) for _ in range(rng.randrange(1, 60)))
        text = data.decode("latin-1")
        try:
            expected = reference(text)
        except ZeroDivisionError:
            with pytest.raises(ZeroDivisionError):
                biofasting.GC123(text)
            continue
        assert biofasting.GC123(text) == expected, data


@needs_biopython
def test_gcg_matches_biopython_over_the_whole_byte_range():
    from Bio.SeqUtils.CheckSum import gcg as reference

    rng = random.Random(33)
    for _ in range(2000):
        data = bytes(rng.randrange(256) for _ in range(rng.randrange(0, 80)))
        text = data.decode("latin-1")
        try:
            expected = reference(text)
        except TypeError as error:
            with pytest.raises(TypeError, match="ord\\(\\) expected a character"):
                biofasting.gcg(text)
            assert "length 2" in str(error)
            continue
        assert biofasting.gcg(text) == expected, data


@needs_biopython
def test_crc64_matches_biopython_over_the_whole_byte_range():
    """No decline needed anywhere: the reference folds each character in with
    ``ord(c) & 0xFF``, which for Latin-1 text *is* the byte."""
    from Bio.SeqUtils.CheckSum import crc64 as reference

    rng = random.Random(34)
    for _ in range(2000):
        data = bytes(rng.randrange(256) for _ in range(rng.randrange(0, 80)))
        assert biofasting.crc64(data.decode("latin-1")) == reference(
            data.decode("latin-1")
        ), data


@needs_biopython
def test_the_crc64_table_is_the_reference_table():
    """The kernel derives its 256-entry table at compile time from the
    reference's own recurrence instead of being handed a dump of it.  This is
    what makes that acceptable: every entry, compared."""
    from Bio.SeqUtils.CheckSum import _table_h

    assert list(biofasting._core.crc64_table_h()) == _table_h


@needs_biopython
@needs_corpus
def test_everything_matches_biopython_on_the_corpus():
    """Real reads and a real chromosome, not random strings."""
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC123, GC_skew, gc_fraction
    from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid

    reads = [sequence for _, sequence, _ in biofasting.open_fastq(
        DATA / "fastq" / "reads_10k.fastq")]
    assert len(reads) >= 1000
    for data in reads[:1000]:
        text = data.decode("ascii")
        assert biofasting.GC123(text) == GC123(text), data
        assert biofasting.GC_skew(text, 30) == GC_skew(text, 30), data
        assert biofasting.gcg(data) == gcg(Seq(data)), data
        assert biofasting.crc64(data) == crc64(Seq(data)), data
        assert biofasting.crc32(data) == crc32(Seq(data)), data
        assert biofasting.seguid(data) == seguid(Seq(data)), data
        for mode in ("remove", "ignore", "weighted"):
            assert biofasting.gc_fraction(data, mode) == gc_fraction(Seq(data), mode)

    index = biofasting.open_fasta(DATA / "fasta" / "genome_1mb.fasta")
    for name in index:
        sequence = index[name]
        text = sequence.decode("ascii")
        assert biofasting.GC123(text) == GC123(text), name
        assert biofasting.GC_skew(text, 100) == GC_skew(text, 100), name
        assert biofasting.gcg(sequence) == gcg(Seq(sequence)), name
        assert biofasting.crc64(sequence) == crc64(Seq(sequence)), name
        assert biofasting.gc_fraction(sequence, "weighted") == gc_fraction(
            Seq(sequence), "weighted"
        ), name
