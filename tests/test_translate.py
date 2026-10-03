"""Translation, against `Bio.Seq.translate`.

Same two layers as the rest of the suite.  The first asserts the contract on
sequences written by hand, against the specification spelled out here in Python
-- neither Biopython nor the corpus is needed, so CI, which has neither, still
tests the kernel and the fallback.  The second is the real evidence: the same
calls through `Bio.Seq.translate` over every genetic code and every codon of
every one of them.

Why the sweeps are this thorough, rather than a handful of proteins.  `translate`
is two implementations -- a table-driven kernel that answers in one pass, and
the reference's own codon loop for everything the kernel declines -- and the
kernel is handed a *cached* pair of tables derived from the table object.  That
cache reads an attribute off the table, and an ambiguous table answers for
attributes it does not have by delegating to the unambiguous table it was built
from.  So a fresh ambiguous table found the unambiguous table's tables: same 64
codons, plausible-looking values, wrong alphabet.  `"NNN"` became
``TranslationError: Codon 'NNN' is invalid`` where the reference answers ``"X"``.

It reproduced only when the unambiguous table had been used *first*, so a smoke
test on the standard code passed while a sweep over every code did not.  The
ordering that found it is kept below as
`test_a_fresh_ambiguous_table_does_not_inherit_the_unambiguous_one`, which is the
two-line version, and the sweep it was found by is
`test_every_codon_of_every_table_agrees`.
"""

import importlib.util
import random
import warnings

import pytest

import biofasting
from biofasting import codon_table, translate
from biofasting._codon_table import TranslationError

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

#: Every genetic code the reference registers, by id.  1-6, 9-16, 21-31; the
#: gaps are ids the reference itself does not have.
TABLE_IDS = sorted(codon_table.ambiguous_generic_by_id)


# -------------------------------------------------------------------------- #
# The specification, written out
# -------------------------------------------------------------------------- #

#: The standard code, as the reference's own `CodonTable` holds it: this is the
#: source of truth for the hand-written layer, typed out rather than imported
#: so that a change to `_codon_tables_data.py` cannot move both sides at once.
_STANDARD = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def spec_translate(text, to_stop=False, stop_symbol="*"):
    """The unambiguous standard code, which is all the hand-written layer needs."""
    out = []
    for i in range(0, len(text) - len(text) % 3, 3):
        residue = _STANDARD[text[i:i + 3]]
        if residue == "*":
            if to_stop:
                break
            residue = stop_symbol
        out.append(residue)
    return "".join(out)


# -------------------------------------------------------------------------- #
# Layer 1: the contract, without Biopython
# -------------------------------------------------------------------------- #


def test_the_standard_code_is_what_it_says():
    """All 64 codons, read out of the table the kernel is handed."""
    assert translate("".join(sorted(_STANDARD))) == "".join(
        _STANDARD[c] for c in sorted(_STANDARD)
    )


@pytest.mark.parametrize("sequence,expected", [
    ("ATGAAATAG", "MK*"),
    ("atgaaatag", "MK*"),          # case-insensitive
    ("ATGAAATAGAAATAA", "MK*K*"),
    ("", ""),
    ("ATG", "M"),
])
def test_default_translation(sequence, expected):
    assert translate(sequence) == expected


def test_a_trailing_partial_codon_is_dropped_after_a_warning():
    with pytest.warns(biofasting.BiopythonWarning, match="Partial codon"):
        assert translate("AT") == ""
    with pytest.warns(biofasting.BiopythonWarning, match="Partial codon"):
        assert translate("ATGAAAAT") == "MK"


def test_to_stop_truncates_at_the_first_stop():
    assert translate("ATGAAATAGAAATAA", to_stop=True) == "MK"
    assert translate("ATGAAATAG", to_stop=True) == "MK"
    #   A stop as the very first codon is an empty protein, not an error.
    assert translate("TAAATG", to_stop=True) == ""


def test_stop_symbol_is_taken_verbatim():
    """The reference never checks it, so neither do we -- including for ""."""
    assert translate("ATGTAATAG", stop_symbol="X") == "MXX"
    assert translate("ATGTAATAG", stop_symbol="") == "M"
    assert translate("ATGTAATAG", stop_symbol="STOP") == "MSTOPSTOP"


def test_gap_validates_even_when_no_gap_occurs():
    """The reference checks `gap` before it looks at a single codon."""
    with pytest.raises(TypeError):
        translate("ATGAAA", gap=1)
    with pytest.raises(ValueError):
        translate("ATGAAA", gap="--")
    #   A one-character gap is accepted, and unused codons never reach it.
    assert translate("ATGAAA", gap="-") == "MK"
    #   "---" is a gap codon only once a gap is named.
    assert translate("ATG---AAA", gap="-") == "M-K"
    with pytest.raises(TranslationError):
        translate("ATG---AAA", gap=None)


def test_a_codon_of_nucleotides_that_is_not_in_the_table_is_X():
    """`pos_stop`: every reading agrees except that some are stops."""
    assert translate("TAN") == "X"
    assert translate("TAR") == "*"        # every reading of TAR *is* a stop
    assert translate("NNN") == "X"
    assert translate("ATGNNN", to_stop=True) == "MX"


def test_a_codon_that_is_not_a_nucleotide_at_all_is_an_error():
    with pytest.raises(TranslationError, match="Codon 'AT-' is invalid"):
        translate("AT-")


def test_cds_checks_start_length_and_stop():
    assert translate("ATGAAATAA", cds=True) == "MK"
    #   The first codon becomes M whatever it was.  TTG is a start of the
    #   standard code; GTG is not, which is why it is not the example here --
    #   a bacterial start set is table 11's business, not table 1's.
    assert translate("TTGAAATAA", cds=True) == "MK"
    with pytest.raises(TranslationError, match="First codon 'GTG'"):
        translate("GTGAAATAA", cds=True)
    with pytest.raises(TranslationError, match="First codon 'AAA'"):
        translate("AAAAAATAA", cds=True)
    with pytest.raises(TranslationError, match="not a multiple of three"):
        translate("ATGAAATAAx", cds=True)
    with pytest.raises(TranslationError, match="Final codon 'AAA'"):
        translate("ATGAAAAAA", cds=True)


def test_cds_rejects_a_stop_in_the_middle():
    with pytest.raises(TranslationError, match="Extra in frame stop"):
        translate("ATGTAAAAATAA", cds=True)


def test_cds_rejects_a_stop_in_the_middle_even_past_the_kernel():
    """The same error, reached through the fallback rather than the kernel.

    `"NNN"` is `pos_stop`, so the kernel declines and the codon loop runs; the
    stop after it still has to be an error rather than a symbol.
    """
    for sequence in ("ATGNNNTAAAAATAA", "ATGNNNTAATAA"):
        with pytest.raises(TranslationError, match="Extra in frame stop"):
            translate(sequence, cds=True)


def test_a_partial_codon_warns_and_is_dropped():
    with pytest.warns(biofasting.BiopythonWarning, match="Partial codon"):
        assert translate("ATGAAAT") == "MK"
    #   ...but not for cds, which refuses the length outright, nor for an
    #   exactly-divisible sequence.  The suite runs with `filterwarnings =
    #   ["error"]`, so an unrequested warning would fail the line below anyway;
    #   the context manager is here to say that it is meant to be silent.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert translate("ATGAAATAA") == "MK*"


@pytest.mark.parametrize("sequence,expected", [
    (b"ATGAAATAG", b"MK*"),
    (bytearray(b"ATGAAATAG"), b"MK*"),
    (memoryview(b"ATGAAATAG"), b"MK*"),
])
def test_bytes_in_bytes_out(sequence, expected):
    assert translate(sequence) == expected


def test_the_table_argument_has_three_spellings():
    assert translate("ATGAAATAG", table="Standard") == "MK*"
    assert translate("ATGAAATAG", table=1) == "MK*"
    assert translate("ATGAAATAG", table=codon_table.ambiguous_generic_by_id[1]) == "MK*"
    assert translate("AUGAAAUAG", table=2) == "MK*"
    assert translate("ATGTGA", table=1) == "M*"
    #   The two remaining lines are the same genetic code under its two
    #   spellings, and they do *not* behave the same -- in the reference either.
    #   `ambiguous_generic_by_id[27]` calls TGA a stop and also calls it W, so
    #   translating through it warns; `ambiguous_generic_by_name["Vertebrate
    #   Mitochondrial"]` is a *different object* built without TGA among its
    #   stops at all, and translating through it does not.  `is` is False and
    #   `stop_codons` differs.  Both answer "MW", which is what a caller sees,
    #   and this port reproduces both, warning included -- a tidier library here
    #   would be a library that disagrees with the one it is checked against.
    with pytest.warns(biofasting.BiopythonWarning, match="both STOP and an amino acid"):
        assert translate("ATGTGA", table=27) == "MW"
    assert translate("ATGTGA", table="Vertebrate Mitochondrial") == "MW"


def test_a_bad_table_argument_is_refused_the_way_the_reference_refuses_it():
    with pytest.raises(ValueError, match="DO NOT.*character string mapping table"):
        translate("ATG", table="Nonsense")
    with pytest.raises(ValueError, match="Bad table argument"):
        translate("ATG", table=object())
    with pytest.raises(ValueError, match="Bad table argument"):
        translate("ATG", table={})
    with pytest.raises(KeyError):
        translate("ATG", table=9999)


def test_a_table_with_a_dual_coding_codon_warns_and_refuses_to_stop():
    """Table 27's TGA is both a stop and W; the reference will not pick one."""
    with pytest.warns(biofasting.BiopythonWarning, match="both STOP and an amino acid"):
        translate("ATGTGA", table=27, to_stop=False)
    with pytest.raises(ValueError, match="cannot use 'to_stop=True'"):
        translate("ATGTGA", table=27, to_stop=True)
    #   ...and the warning is about the *table*, so a sequence with no such
    #   codon in it still warns.
    with pytest.warns(biofasting.BiopythonWarning, match="both STOP and an amino acid"):
        translate("ATGAAA", table=27)


# -------------------------------------------------------------------------- #
# Layer 2: the differential sweeps
# -------------------------------------------------------------------------- #

_ALPHABETS = [
    "ACGT",
    "ACGTN",
    "ACGTNRYWSMKHBVDN",
    "acgtn",
    "ACGU",
    "ACGT-",
    "ACGT*X",
    "ACDEFGHIKLMNPQRSTVWY",
]


def _random_sequences(rng, count, lengths):
    for _ in range(count):
        length = rng.choice(lengths)
        alphabet = rng.choice(_ALPHABETS)
        yield "".join(rng.choice(alphabet) for _ in range(length))


def _both(sequence, **options):
    """`(ours, the reference's)` for the same call, or the two exceptions."""
    from Bio.Seq import translate as reference_translate

    def call(function):
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                return ("ok", function(sequence, **options),
                        [str(w.message) for w in caught])
        except Exception as error:  # noqa: BLE001 - the exception *is* the answer
            return ("raised", type(error).__name__, str(error))

    return call(translate), call(reference_translate)


@needs_biopython
def test_a_fresh_ambiguous_table_does_not_inherit_the_unambiguous_one():
    """The cache bug, in the calls that reproduce it, in order.

    The kernel's `code` and `amino` are cached on the table object.  An
    ambiguous table answers for attributes it does not have by delegating to the
    unambiguous table it *was* built from, so a `getattr` for that cache found the
    unambiguous table's -- same 64 codons, different alphabet -- and `"NNN"`
    raised `TranslationError` where the reference answers `"X"`.

    Both halves have to be the *right* pair.  `ambiguous_dna_by_id[1]` delegates
    to `unambiguous_dna_by_id[1]`, and `ambiguous_generic_by_id[1]` to
    `generic_by_id[1]`; those are the two pairs, and passing the table objects
    themselves is what pins the pair rather than whatever `table=1` happens to
    resolve to this month.  The order is the other half: the unambiguous table has
    to be used first, or there is nothing cached for the delegation to find.  That
    is why this is a named test and not a line inside the sweep -- the sweep takes
    three seconds and this takes microseconds, and when it breaks, it says which
    bug came back.
    """
    import biofasting.codon_table as tables

    for plain, ambiguous, sequence, expected in [
        (tables.unambiguous_dna_by_id[1], tables.ambiguous_dna_by_id[1], "NNN", "X"),
        (tables.generic_by_id[1], tables.ambiguous_generic_by_id[1], "NNN", "X"),
    ]:
        #   Unambiguous first, which is what warms the cache the delegation finds.
        assert translate("ATGAAA", table=plain) == "MK"
        #   ...and now the ambiguous table built from it must answer its own way.
        assert translate(sequence, table=ambiguous) == expected
        assert translate("TAR", table=ambiguous) == "*"
        assert translate("ATGAAA", table=ambiguous) == "MK"

    #   And the same through the spellings a caller actually types, so the
    #   resolution path is covered too and not just the table objects.
    assert translate("ATGAAA", table=1) == "MK"
    assert translate("NNN", table=1) == "X"
    assert translate("TAN", table=1) == "X"
    assert translate("AUGNNN", table=2) == "MX"


@needs_biopython
def test_every_codon_of_every_table_agrees():
    """Every genetic code, every codon of it, with and without to_stop.

    This is the sweep that found the ambiguity-table cache bug: 27 codes times
    4,096 codons -- all sixteen ambiguity letters in all three positions --
    times two option shapes, 221,184 comparisons, with the tables primed in the
    orders below.
    """
    #   The order matters, and it is the whole point: the unambiguous table
    #   first, then the ambiguous one built from it.  Reversed, the delegation
    #   in `AmbiguousCodonTable.__getattr__` has nothing to find.
    for name in ("unambiguous_dna_by_id", "ambiguous_dna_by_id",
                 "unambiguous_rna_by_id", "ambiguous_rna_by_id",
                 "generic_by_id", "ambiguous_generic_by_id"):
        getattr(codon_table, name)[1]

    letters = "ACGTNRYWSMKHBVDN"
    codons = [a + b + c for a in letters for b in letters for c in letters]
    mismatches = []
    for table_id in TABLE_IDS:
        for to_stop in (False, True):
            for codon in codons:
                ours, theirs = _both(codon, table=table_id, to_stop=to_stop)
                if ours != theirs:
                    mismatches.append((table_id, codon, to_stop, ours, theirs))
    assert not mismatches, f"{len(mismatches)} of {len(TABLE_IDS) * 2 * len(codons)} disagree"


@needs_biopython
def test_random_sequences_agree_over_every_option_set():
    rng = random.Random(20261003)
    sequences = list(_random_sequences(rng, 400, [0, 1, 2, 3, 5, 9, 12, 30, 150, 999]))
    mismatches = []
    for sequence in sequences:
        for options in (
            {},
            {"to_stop": True},
            {"stop_symbol": ""},
            {"stop_symbol": "STOP"},
            {"table": 2},
            {"table": 11},
            {"table": 27},
            {"gap": "-"},
            {"gap": "-", "stop_symbol": "X"},
        ):
            ours, theirs = _both(sequence, **options)
            if ours != theirs:
                mismatches.append((sequence[:40], options, ours, theirs))
    assert not mismatches, f"{len(mismatches)} mismatches: {mismatches[:3]}"


@needs_biopython
def test_random_sequences_agree_under_cds():
    """`cds` is all checks, so it is swept separately over starts and stops."""
    rng = random.Random(1)
    starts = ["ATG", "GTG", "TTG", "CTG", "AAA", "NNN", "AT"]
    stops = ["TAA", "TAG", "TGA", "TAAx", "AAA"]
    middles = ["", "AAA", "NNN", "AT-", "TAN", "TAR", "TAA", "GGG" * 7]
    mismatches = []
    for start in starts:
        for middle in middles:
            for stop in stops:
                ours, theirs = _both(start + middle + stop, cds=True)
                if ours != theirs:
                    mismatches.append((start + middle + stop, ours, theirs))
    assert not mismatches, f"{len(mismatches)} mismatches: {mismatches[:3]}"


@needs_biopython
@pytest.mark.parametrize("stop_symbol", ["*", "", "X", "STOP", "é", 5, ["x"], None, 3.5])
def test_stop_symbol_is_whatever_the_caller_passed(stop_symbol):
    for sequence in ("ATGAAATAG", "ATGAAA", "ATGNNNTAG"):
        ours, theirs = _both(sequence, stop_symbol=stop_symbol)
        assert ours == theirs, (sequence, stop_symbol, ours, theirs)


@needs_biopython
@pytest.mark.parametrize("gap", [None, "-", ".", "--", "", 1, b"-"])
def test_gap_is_refused_the_way_the_reference_refuses_it(gap):
    for sequence in ("ATGAAATAG", "ATG---AAA", "ATGNNN"):
        ours, theirs = _both(sequence, gap=gap)
        assert ours == theirs, (sequence, gap, ours, theirs)


@needs_biopython
def test_a_non_ascii_stop_symbol_survives():
    """The kernel returns a `str`, not bytes, so a multi-byte symbol is intact."""
    from Bio.Seq import translate as reference_translate

    assert translate("ATGAAATAG", stop_symbol="é") == reference_translate(
        "ATGAAATAG", stop_symbol="é"
    )


@needs_biopython
def test_the_two_implementations_in_this_package_agree_with_each_other():
    """The kernel against the fallback: same answer, or a prefix of it.

    The differential sweeps above compare the *public* function with
    Biopython's, which is what a caller cares about; this one opens the package
    and asks whether the two paths inside it are the same function.  Without
    Biopython it still runs, and it is the only test that says so directly: the
    sweeps would pass with a kernel that answered correctly and a fallback that
    did not, on every input that never reaches the fallback.

    A decline is not a disagreement, so it cannot be asserted as one.  The
    kernel gives up on the first codon whose byte is not one of the table's four
    letters, and an *ambiguous* codon is not that -- ``"BCG"`` is three legal
    nucleotides and a residue to be looked up, and the loop answers ``"X"`` for
    it where the kernel has no entry to consult.  So a decline means the loop
    takes over, not that the loop fails; what has to hold is that everything the
    kernel did emit agrees with the loop's own beginning, codon for codon, up to
    where it stopped.  That is the property the fast path is allowed to rely on.
    """
    from biofasting._translate import _kernel_tables, _resolve_table, _translate_loop

    rng = random.Random(7)
    declined = served = 0
    for sequence in _random_sequences(rng, 400, [3, 9, 30, 150]):
        buffer = sequence.upper().encode("latin-1")
        table = _resolve_table("Standard")
        tables = _kernel_tables(table)
        assert tables.code is not None
        for to_stop in (False, True):
            outcome, codon, kernel = biofasting._core.translate(
                buffer, tables.code, tables.amino, "*", "X", to_stop, False
            )
            try:
                loop = _translate_loop(
                    buffer, 0, len(buffer), table, tables.valid_letters, "*",
                    to_stop, False, None,
                )
            except TranslationError:
                loop = None
            if outcome == 0:
                served += 1
                assert loop is not None, (sequence, to_stop)
                assert kernel == loop, (sequence, to_stop, kernel, loop)
            else:
                declined += 1
                #   `kernel` holds the codons before the one that stopped it --
                #   every one of which the kernel decided, so the loop must agree
                #   with it on all of them.  Slicing the input down to those
                #   codons and running the loop over just that is the comparison
                #   that survives the loop *raising* at the decline codon, which
                #   it does when the reason for the decline was a byte that is
                #   not a nucleotide at all: the loop has no answer to return,
                #   but it still has to agree about everything before.
                stem = buffer[:3 * codon]
                try:
                    head = _translate_loop(
                        stem, 0, len(stem), table, tables.valid_letters, "*",
                        to_stop, False, None,
                    )
                except TranslationError:
                    head = None
                assert head is not None, (sequence, to_stop, codon)
                assert kernel == head, (sequence, to_stop, codon, kernel, head, loop)
    #   Both arms have to be exercised, or the loop above proved only one of
    #   them: these alphabets are chosen so that random codons are ambiguous,
    #   invalid, and clean in roughly equal measure.
    assert declined and served, (declined, served)


@needs_biopython
def test_agrees_over_the_real_corpus():
    """The reads the benchmark translates, against the reference."""
    from pathlib import Path

    import Bio.SeqIO
    from Bio.Seq import translate as reference_translate

    fastq = Path(__file__).resolve().parent.parent / "bench" / "data" / "fastq" / "reads_10k.fastq"
    if not fastq.exists():
        pytest.skip("benchmark corpus not generated (python bench/gen_data.py)")
    reads = []
    #   Opened and closed by hand: the suite raises on unraisable exceptions,
    #   and a `SeqIO.parse` left to the collector closes its file too late.
    with fastq.open() as handle:
        for record in Bio.SeqIO.parse(handle, "fastq"):
            reads.append(str(record.seq))
            if len(reads) >= 2000:
                break
    assert len(reads) == 2000
    for read in reads:
        assert translate(read) == reference_translate(read)
