"""`Bio.SeqUtils.ProtParam` and `IsoelectricPoint`, against Biopython.

The same two layers as everywhere else in this project.  The first asserts
behaviour on sequences chosen by hand, against a specification written out in
Python here, so that a CI without Biopython still tests the kernels.  The second
is the evidence: the same calls through `Bio.SeqUtils.ProtParam`, over random
proteins, over every one of the reference's thirty-two scales, and over the
inputs where the reference raises or writes a warning.

Four things in here are not repetitions of each other and are worth naming.

`test_the_gravy_flag_is_not_a_boolean` and
`test_the_first_float_after_an_integer_run_is_added_plainly_too` exist because
CPython's `sum()` does not treat every number alike.  It compensates a float term,
adds an `int` one plainly, and starts in an *integer* accumulator so that the
first float term is a plain add too.  Ten of the reference's twenty-eight gravy
scales are written with integer values, and the generated tables under this
package had kept the reference's values while flattening all of them to floats;
four of those ten scales then disagreed with the reference on between 4% and 29%
of random proteins.  The two tests pin the kernel's half of that; the differential
layer pins the data's.

`test_the_instability_kernel_is_not_the_compensated_one` is the other half of the
same surprise.  `gravy` and `molecular_weight` are written `sum(...)` and
`instability_index` is written `score += value`, and both summations live in the
same class; a kernel that used one for both would be wrong on the second.

`test_a_scale_is_recognised_by_identity_and_not_by_value` records a deliberate
trade.  The kernel path needs the *same dict object* the reference ships, not an
equal one, so a caller who builds their own mapping silently gets the reference's
own loop -- correct, and about twenty times slower.  Stating it here makes that a
decision on the record instead of something to be discovered from a profile.
"""

import contextlib
import doctest
import importlib.util
import io
import random
import struct
import sys

import pytest

from biofasting import _core, isoelectric_point, protparam, protparam_data
from biofasting.protparam import ProteinAnalysis
from biofasting.protparam_data import DIWV, Flex, SCALES_BY_NAME, cw, gravy_scales

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

#: `sum()` compensates a float term from CPython 3.12 on and adds it plainly
#: before that, and the kernel implements the 3.12 rule -- `sequtils.hpp` is
#: where that is written down.  Two of the tests below compare the kernel with
#: the interpreter's own `sum`, so they carry an interpreter in them and the
#: comparison is only true of one that has the rule.  The mode distinction they
#: are about is asserted on every version; the part that is a fact about 3.12 is
#: asserted where it is a fact.  What that leaves behind is pinned by
#: `test_the_pre_3_12_sum_is_not_the_one_this_package_uses` below.
compensated_sum = sys.version_info >= (3, 12)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

#: The protein in the reference's own class docstring, which this package's
#: `ProteinAnalysis` documents too, so that the two docstrings agree.
DOCTEST_PROTEIN = (
    "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT"
    "RDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEEC"
    "LFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILF"
    "LPLPV"
)


# --------------------------------------------------------------------------
# The specification, written out from the reference.
# --------------------------------------------------------------------------


def spec_instability_index(text, index):
    """`instability_index`'s dipeptide loop, plain `+=` and all."""
    score = 0.0
    for i in range(len(text) - 1):
        this, following = text[i : i + 2]
        score += index[this][following]
    return (10.0 / len(text)) * score


def spec_flexibility(text, flex):
    """`flexibility`'s window of nine, whose middle is read at index 5."""
    weights = (0.25, 0.4375, 0.625, 0.8125)
    out = []
    for i in range(len(text) - 9):
        score = 0.0
        for j in range(4):
            score += (flex[text[i + j]] + flex[text[i + 8 - j]]) * weights[j]
        score += flex[text[i + 5]]
        out.append(score / 5.25)
    return out


def spec_gravy(text, scale):
    """`gravy`'s sum, whose compensation depends on each value's type."""
    return sum(scale[aa] for aa in text) / len(text)


def spec_protein_scale(text, scale, window, edge):
    """`protein_scale`'s two-sided window, warning-and-skipping and all."""
    weights = [
        edge + (2 * (1.0 - edge) / (window - 1)) * i for i in range(window // 2)
    ]
    sum_of_weights = sum(weights) * 2 + 1
    out = []
    for i in range(len(text) - window + 1):
        sub = text[i : i + window]
        score = 0.0
        for j in range(window // 2):
            try:
                front = scale[sub[j]]
                back = scale[sub[window - j - 1]]
                score += weights[j] * front + weights[j] * back
            except KeyError:
                continue
        if sub[window // 2] in scale:
            score += scale[sub[window // 2]]
        out.append(score / sum_of_weights)
    return out


def spec_charge(text, pH, pos_pKs, neg_pKs):
    """`charge_at_pH`, in the reference's two loops and its insertion order."""
    counts = {aa: text.count(aa) for aa in AMINO_ACIDS}
    charged = {aa: float(counts[aa]) for aa in "KRHDECY"}
    charged["Nterm"] = 1.0
    charged["Cterm"] = 1.0
    positive = 0.0
    for aa, pK in pos_pKs.items():
        positive += charged[aa] * (1.0 / (10 ** (pH - pK) + 1.0))
    negative = 0.0
    for aa, pK in neg_pKs.items():
        negative += charged[aa] * (1.0 / (10 ** (pK - pH) + 1.0))
    return positive - negative


def spec_pKs(text):
    """The two pK tables for one sequence, with the terminal overrides in place.

    The reference copies the module's dicts and overrides `Nterm`/`Cterm` in
    place, so a residue with no override leaves the module's value where it was
    -- and where it was is the insertion position, which `charge_at_pH` sums in.
    """
    pos = dict(isoelectric_point.positive_pKs)
    neg = dict(isoelectric_point.negative_pKs)
    if text[0] in isoelectric_point.pKnterminal:
        pos["Nterm"] = isoelectric_point.pKnterminal[text[0]]
    if text[-1] in isoelectric_point.pKcterminal:
        neg["Cterm"] = isoelectric_point.pKcterminal[text[-1]]
    return pos, neg


def spec_pi(text, pH=7.775, min_=4.05, max_=12):
    """`pi`'s bisection, reproduced rather than solved for its root."""
    pos, neg = spec_pKs(text)
    charge = spec_charge(text, pH, pos, neg)
    if max_ - min_ > 0.0001:
        if charge > 0.0:
            min_ = pH
        else:
            max_ = pH
        return spec_pi(text, (min_ + max_) / 2, min_, max_)
    return pH


def raw_table(entries):
    """The `(256 doubles, 256 flags)` pair the sum kernel takes, from a list.

    `entries` is `(byte, value, flag)` triples, so that a test can hand the kernel
    a table no scale in this package has and say what each of its entries is --
    rather than writing a bare 1 or 2 into a flags array and hoping.
    """
    values = bytearray(256 * 8)
    flags = bytearray(256)
    for byte, value, flag in entries:
        struct.pack_into("<d", values, byte * 8, value)
        flags[byte] = flag
    return bytes(values), bytes(flags)


def call(function):
    """A call's value, or its exception as a comparable tuple.

    Nearly every case below is an input where the reference raises, and *which*
    exception, with which arguments, is the answer being compared rather than a
    failure of the test.
    """
    try:
        return ("ok", function())
    except Exception as error:  # noqa: BLE001 - the exception *is* the answer
        return ("raised", type(error).__name__, error.args)


# --------------------------------------------------------------------------
# Hand-chosen sequences, against the specification above.
# --------------------------------------------------------------------------


def test_the_reference_doctest_protein_answers_what_the_reference_prints():
    """One protein, every method, the numbers Biopython prints for it."""
    protein = ProteinAnalysis(DOCTEST_PROTEIN)
    assert len(DOCTEST_PROTEIN) == 152
    assert protein.length == 152
    assert protein.count_amino_acids()["A"] == 6
    assert protein.count_amino_acids()["E"] == 12
    assert protein.count_amino_acids()["W"] == 1
    assert set(protein.count_amino_acids()) == set(AMINO_ACIDS)
    assert round(protein.amino_acids_percent["A"], 2) == 3.95
    assert round(protein.amino_acids_percent["L"], 2) == 11.84
    assert round(protein.molecular_weight(), 2) == 17103.16
    assert round(ProteinAnalysis(DOCTEST_PROTEIN, True).molecular_weight(), 2) == (
        17092.61
    )
    assert round(protein.aromaticity(), 2) == 0.10
    assert round(protein.instability_index(), 2) == 41.98
    assert round(protein.isoelectric_point(), 2) == 7.72
    assert round(protein.gravy(), 2) == -0.60
    helix, turn, sheet = protein.secondary_structure_fraction()
    assert (round(helix, 2), round(turn, 2), round(sheet, 2)) == (0.33, 0.29, 0.37)
    assert protein.molar_extinction_coefficient() == (17420, 17545)
    assert len(protein.flexibility()) == 152 - 9
    assert round(protein.charge_at_pH(7.0), 2) == 0.93


def test_monoisotopic_reaches_the_mass_and_nothing_else():
    average = ProteinAnalysis(DOCTEST_PROTEIN)
    mono = ProteinAnalysis(DOCTEST_PROTEIN, monoisotopic=True)
    assert mono.molecular_weight() != average.molecular_weight()
    assert mono.gravy() == average.gravy()
    assert mono.instability_index() == average.instability_index()
    assert mono.aromaticity() == average.aromaticity()
    assert mono.isoelectric_point() == average.isoelectric_point()


def test_the_sequence_is_upper_cased_once_at_construction():
    lower = ProteinAnalysis("acdefghik")
    assert lower.sequence == "ACDEFGHIK"
    assert lower.count_amino_acids()["A"] == 1
    assert lower.length == 9
    assert lower.instability_index() == (
        ProteinAnalysis("ACDEFGHIK").instability_index()
    )


def test_count_amino_acids_caches_into_the_object():
    """The cache is the reference's, and it is visible from outside."""
    protein = ProteinAnalysis("ACDEFG")
    counts = protein.count_amino_acids()
    assert counts is protein.count_amino_acids()
    counts["A"] = 999
    assert protein.count_amino_acids()["A"] == 999


def test_count_residues_counts_the_twenty_and_ignores_everything_else():
    """The kernel never declines, and this is why it does not have to.

    `str.count('A')` looks for one letter and is uninterested in what it does
    not find, so a residue outside the twenty is a zero and not an error -- the
    opposite of every other kernel in `protparam.cpp`, and the one place where
    "a kernel counts, Python decides" has nothing to decide.
    """
    counts = _core.count_residues("AXBZJUO*".encode("ascii"), protparam._AA_MAP)
    assert list(counts) == [1 if aa == "A" else 0 for aa in AMINO_ACIDS]
    assert _core.count_residues(b"", protparam._AA_MAP) == [0] * 20
    #   Every byte a `str` can hold in Latin-1, in one call.  All 256 of them are
    #   in there and each of the twenty letters is one of them, so an exact count
    #   of one apiece is the whole of what the buckets may contain -- a 0xFF from
    #   the map landing in the bucket that is copied out would show up here as a
    #   2 somewhere, and a byte silently dropped as a 0.
    assert _core.count_residues(bytes(range(256)), protparam._AA_MAP) == [1] * 20
    #   And a byte of every value, none of which is one of the twenty: `bytes`
    #   cannot hold anything else, so this is the whole of the space the kernel
    #   can be handed.
    assert _core.count_residues(bytes(range(256)) * 3, protparam._AA_MAP) == [3] * 20


def test_count_amino_acids_keeps_the_reference_s_letters_and_order():
    """The dict a caller gets: the twenty standard letters, in that order."""
    counts = ProteinAnalysis("AXBZJUO*acdefg").count_amino_acids()
    assert list(counts) == list(AMINO_ACIDS)
    assert counts["A"] == 2 and counts["C"] == 1 and counts["G"] == 1
    assert counts["X"] if "X" in counts else True
    assert set(counts) == set(AMINO_ACIDS)


def test_amino_acids_percent_is_a_property_and_not_a_method():
    protein = ProteinAnalysis("AACC")
    assert protein.amino_acids_percent["A"] == 50.0
    with pytest.raises(TypeError):
        protein.amino_acids_percent()


def test_the_percentages_of_the_empty_sequence_raise_the_integer_division():
    """The reference divides a count by the length, so this is `0 / 0` on ints."""
    with pytest.raises(ZeroDivisionError) as info:
        ProteinAnalysis("").amino_acids_percent
    assert str(info.value) == "division by zero"


def test_gravy_of_the_empty_sequence_raises_the_integer_division():
    """And here the reference's `sum` over nothing is the *integer* zero."""
    with pytest.raises(ZeroDivisionError) as info:
        ProteinAnalysis("").gravy()
    assert str(info.value) == "division by zero"


def test_instability_index_of_the_empty_sequence_raises_the_float_one():
    """Where the divisor is the literal `10.0`, so the message says float."""
    with pytest.raises(ZeroDivisionError) as info:
        ProteinAnalysis("").instability_index()
    assert str(info.value) == "float division by zero"


def test_gravy_of_an_unknown_scale_names_it():
    with pytest.raises(ValueError) as info:
        ProteinAnalysis("ACDE").gravy("NotAScale")
    assert str(info.value) == "scale: NotAScale not known"


def test_flexibility_produces_one_window_fewer_than_fit():
    """`range(length - 9)`, not `range(length - 8)`."""
    assert len(ProteinAnalysis("A" * 9).flexibility()) == 0
    assert len(ProteinAnalysis("A" * 10).flexibility()) == 1
    assert len(ProteinAnalysis("A" * 100).flexibility()) == 91
    assert ProteinAnalysis("A" * 100).flexibility() == spec_flexibility("A" * 100, Flex)


def test_flexibility_reads_the_middle_of_its_window_at_index_five():
    """Index 4 of the window is never read, and index 5 is read twice.

    A ten-residue sequence has one window, and replacing its residue 4 changes
    nothing while replacing its residue 5 changes everything.  That is the
    reference's off-by-one, kept because a repaired one would answer differently.
    """
    baseline = ProteinAnalysis("AAAAAAAAAA").flexibility()
    assert ProteinAnalysis("AAAATAAAAA").flexibility() == baseline
    assert ProteinAnalysis("AAAAATAAAA").flexibility() != baseline


def test_instability_index_counts_one_dipeptide_fewer_than_the_length():
    """The multiply by `10.0 / length` happens after the sum, in that order."""
    assert ProteinAnalysis("A").instability_index() == 0.0
    assert ProteinAnalysis("AA").instability_index() == (10.0 / 2) * DIWV["A"]["A"]
    for text in ("AAAA", "ACDEFGHIKLMNPQRSTVWY", "PETER", "INGAR", "EGGGEGGEGG"):
        assert ProteinAnalysis(text).instability_index() == spec_instability_index(
            text, DIWV
        )


def test_protein_scale_of_the_empty_sequence_is_the_empty_list():
    assert ProteinAnalysis("").protein_scale(SCALES_BY_NAME["kd"], 3) == []


def test_protein_scale_window_one_raises_while_it_builds_the_weights():
    """`2 * (1 - edge) / (window - 1)` is a division by zero at window 1."""
    with pytest.raises(ZeroDivisionError) as info:
        ProteinAnalysis("ACDEFG").protein_scale(SCALES_BY_NAME["kd"], 1)
    assert str(info.value) == "float division by zero"
    with pytest.raises(ZeroDivisionError):
        ProteinAnalysis("ACDEFG").protein_scale({"A": 1.0}, 1)


def test_protein_scale_window_zero_reads_past_the_empty_slice():
    with pytest.raises(IndexError):
        ProteinAnalysis("ACDEFG").protein_scale(SCALES_BY_NAME["kd"], 0)


def test_protein_scale_hand_checked_values():
    """A caller's own scale, which takes the reference's loop, and ours."""
    borrowed = {"A": 1.0, "C": 2.0, "D": 3.0}
    for window, edge in ((2, 1.0), (3, 0.5), (4, 1.0), (6, 0.2)):
        assert ProteinAnalysis("ACDACD").protein_scale(
            borrowed, window, edge
        ) == spec_protein_scale("ACDACD", borrowed, window, edge)
    #   And a scale this package ships, which takes the kernel.
    scale = SCALES_BY_NAME["kd"]
    for window, edge in ((3, 1.0), (9, 0.4), (11, 0.0), (13, 1.7)):
        assert ProteinAnalysis(DOCTEST_PROTEIN).protein_scale(
            scale, window, edge
        ) == spec_protein_scale(DOCTEST_PROTEIN, scale, window, edge)


def test_protein_scale_warns_and_skips_rather_than_raising(capsys):
    """A residue the scale lacks costs its window one term and a line on stderr.

    The warning is the reference's own text, down to the newline, and the window
    is still scored: this is not an error path.
    """
    scale = {"A": 1.0, "C": 2.0}
    assert ProteinAnalysis("ACXAC").protein_scale(scale, 3, 1.0) == (
        spec_protein_scale("ACXAC", scale, 3, 1.0)
    )
    assert capsys.readouterr().err == (
        "warning: A or X is not a standard amino acid.\n"
        "warning: X is not a standard amino acid.\n"
        "warning: X or C is not a standard amino acid.\n"
    )


def test_protein_scale_warns_once_for_a_pair_of_two_unknown_residues(capsys):
    """One line names a pair, however many of the pair are unknown."""
    ProteinAnalysis("XXA").protein_scale({"A": 1.0}, 3, 1.0)
    assert capsys.readouterr().err == (
        "warning: X or A is not a standard amino acid.\n"
        "warning: X is not a standard amino acid.\n"
    )


def test_instability_index_names_the_residue_the_reference_names():
    """`this` before `following`, because that is the order it is indexed in."""
    with pytest.raises(KeyError) as info:
        ProteinAnalysis("ACB").instability_index()
    assert info.value.args == ("B",)
    with pytest.raises(KeyError) as info:
        ProteinAnalysis("AXC").instability_index()
    assert info.value.args == ("X",)


def test_flexibility_names_the_residue_the_reference_names():
    """A sequence long enough for a window to reach the odd residue.

    A window of nine starts at every offset but the last nine, so a residue at the
    very end is read by no window at all and raises nothing.  That is a property of
    the reference, and the reason this one does not sit at the end.
    """
    with pytest.raises(KeyError) as info:
        ProteinAnalysis("A" * 9 + "X" + "A" * 9).flexibility()
    assert info.value.args == ("X",)


def test_isoelectric_point_hand_checked_values():
    """The reference's own doctest, whose two numbers it asserts too."""
    ingar = isoelectric_point.IsoelectricPoint("INGAR")
    assert round(ingar.pi(), 2) == 9.75
    assert round(ingar.charge_at_pH(7.0), 2) == 0.76
    for text in ("A", "ACDEFGHIKLMNPQRSTVWY", "INGAR", DOCTEST_PROTEIN):
        point = isoelectric_point.IsoelectricPoint(text)
        assert point.pi() == spec_pi(text)
        pos, neg = spec_pKs(text)
        for pH in (1.0, 4.05, 7.0, 7.775, 12.0):
            assert point.charge_at_pH(pH) == spec_charge(text, pH, pos, neg)


def test_isoelectric_point_takes_a_residue_count_it_is_given():
    """`aa_content` is the reference's short cut, and an empty dict is not one."""
    counts = ProteinAnalysis("ACDEFGHIK").count_amino_acids()
    assert isoelectric_point.IsoelectricPoint("ACDEFGHIK", counts).pi() == spec_pi(
        "ACDEFGHIK"
    )


def test_charge_at_pH_overrides_the_terminal_pKs_in_place():
    """An overridden `Nterm` keeps its position, and the position is summed."""
    acidic = isoelectric_point.IsoelectricPoint("ACDE")
    assert list(acidic.pos_pKs) == ["Nterm", "K", "R", "H"]
    assert acidic.pos_pKs["Nterm"] == isoelectric_point.pKnterminal["A"]
    assert list(acidic.neg_pKs) == ["Cterm", "D", "E", "C", "Y"]
    assert acidic.neg_pKs["Cterm"] == isoelectric_point.pKcterminal["E"]
    #   A residue with no override leaves the module's own dicts untouched.
    assert isoelectric_point.IsoelectricPoint("GCDE").pos_pKs["Nterm"] == (
        isoelectric_point.positive_pKs["Nterm"]
    )
    assert isoelectric_point.positive_pKs["Nterm"] == 7.5
    assert list(isoelectric_point.negative_pKs) == ["Cterm", "D", "E", "C", "Y"]


# --------------------------------------------------------------------------
# The tables, and the three things CPython's sum() does.
# --------------------------------------------------------------------------


def test_every_shipped_scale_has_one_value_per_standard_residue():
    for name, scale in SCALES_BY_NAME.items():
        assert list(scale) == list(AMINO_ACIDS), name


def test_gravy_scales_points_at_the_module_scales_and_not_at_copies():
    """`gravy_scales["KyteDoolitle"] is kd`, which is what the kernel keys on."""
    assert gravy_scales["KyteDoolitle"] is SCALES_BY_NAME["kd"]
    assert gravy_scales["Cowan3.4"] is cw[3.4]
    assert gravy_scales["Cowan7.5"] is cw[7.5]
    assert SCALES_BY_NAME["Flex"] is Flex
    for name, scale in gravy_scales.items():
        assert any(scale is other for other in SCALES_BY_NAME.values()), name
    assert SCALES_BY_NAME["kd"]["I"] == 4.5


def test_a_scale_is_recognised_by_identity_and_not_by_value():
    """A copy of one of our scales takes the reference's loop, not the kernel.

    Correct either way, and about twenty times slower the second way, so which
    one is in hand is worth being able to see.  A mapping whose keys are not
    single bytes, or whose values are not numbers CPython's `sum()` has a rule
    for, is not one a dense table can carry, and is declined for the same reason.
    """
    assert protparam._scale_blob(SCALES_BY_NAME["kd"]) is not None
    assert protparam._scale_blob(dict(SCALES_BY_NAME["kd"])) is None
    assert protparam._scale_blob({"A": 1.0}) is None
    assert protparam._dense_scale({"A": 1.0, "C": "two"}) is None
    assert protparam._dense_scale({"A": 1.0, "CC": 2.0}) is None
    assert protparam._dense_scale({"A": 1.0, "€": 2.0}) is None
    assert protparam._dense_scale({"A": 1.0, "é": 2.0}) is not None


def test_the_gravy_flag_is_not_a_boolean():
    """Compensating an `int` term and adding it plainly are different answers.

    `sum([-8, -1.9, 2.1, 4.2])` is -3.6 in CPython, and `sum([-8.0, -1.9, 2.1,
    4.2])` is -3.5999999999999996.  The only thing between the two is the type of
    the first value, which is what the flag carries -- and it is what makes the
    reference answer -3.6 for "IYAR" under the Parker scale.

    The kernel's own answer is the same on every interpreter, so it is asserted
    as the literal it is; agreeing with `sum()` is a second claim, and it is made
    only where the interpreter's `sum` is the rule the kernel implements.
    """
    values = (-8.0, -1.9, 2.1, 4.2)
    floats, all_float = raw_table(
        [(byte, value, _core.kCompensated) for byte, value in zip(b"ABCD", values)]
    )
    total, bad = _core.molecular_weight_mass(b"ABCD", floats, all_float)
    assert bad == -1
    assert total == -3.5999999999999996
    if compensated_sum:
        assert total == sum(values)

    mixed, with_int = raw_table(
        [
            (ord("A"), -8, _core.kPlain),
            (ord("B"), -1.9, _core.kCompensated),
            (ord("C"), 2.1, _core.kCompensated),
            (ord("D"), 4.2, _core.kCompensated),
        ]
    )
    total, bad = _core.molecular_weight_mass(b"ABCD", mixed, with_int)
    assert bad == -1
    assert total == -3.6
    if compensated_sum:
        assert total == sum([-8, -1.9, 2.1, 4.2])
    assert total != _core.molecular_weight_mass(b"ABCD", floats, all_float)[0]


@pytest.mark.skipif(compensated_sum, reason="CPython 3.12 and later compensate")
def test_the_pre_3_12_sum_is_not_the_one_this_package_uses():
    """The one place this package and the reference part company, on the record.

    On 3.10 and 3.11 `sum()` adds a float plainly, and the kernel has a mode that
    reproduces that exactly -- but `_dense_scale` chooses the flag from the
    *type* of a table value and not from the interpreter, so a float entry is
    marked compensated there too and `gravy` and `molecular_weight` can answer a
    last bit differently from the reference.  Two ways out, and both are a
    decision rather than a patch: choose the flag on `sys.version_info`, or raise
    the floor to the interpreter that has the rule.  This test asserts the
    divergence is there, so that either one turns it red and the state of the
    package is read from a test instead of being inferred from a version.
    """
    scale = gravy_scales["Parker"]
    _, _, flags = protparam._dense_scale(scale)
    floats = [letter for letter, value in scale.items() if type(value) is float]
    assert floats
    assert all(flags[ord(letter)] == _core.kCompensated for letter in floats)

    values = (-8.0, -1.9, 2.1, 4.2)
    compensated, compensated_flags = raw_table(
        [(byte, value, _core.kCompensated) for byte, value in zip(b"ABCD", values)]
    )
    plain, plain_flags = raw_table(
        [(byte, value, _core.kPlain) for byte, value in zip(b"ABCD", values)]
    )
    #   The mode the package uses against the interpreter's rule, and the mode it
    #   does not against the same rule: the two disagree, and the second agrees.
    by_mode = _core.molecular_weight_mass(b"ABCD", compensated, compensated_flags)
    assert by_mode[0] != sum(values)
    assert _core.molecular_weight_mass(b"ABCD", plain, plain_flags)[0] == sum(values)


def test_the_first_float_after_an_integer_run_is_added_plainly_too():
    """CPython's accumulator starts as an `int` and stays one until a float arrives.

    That arriving float is a plain add, with the compensation still untouched:
    it is the *second* float onward that is a Neumaier step.  `sum([-8., -5.,
    4.2, 0.9])` is -7.8999999999999995, and the same two ints followed by the
    same two floats is -7.9.  The kernel carries the phase, and the difference is
    observable rather than merely written down.

    What the phase cannot reproduce is the accumulator itself: CPython's is an
    arbitrary-precision `int` while this one is a double, so the two part company
    once a run of integer terms passes 2**53.  See the note in `sequtils.hpp`;
    with the largest value any shipped scale holds that is about 7e14 residues.
    """
    values = (-8, -5, 4.2, 0.9)
    table, flags = raw_table(
        [
            (ord("A"), values[0], _core.kPlain),
            (ord("B"), values[1], _core.kPlain),
            (ord("C"), values[2], _core.kCompensated),
            (ord("D"), values[3], _core.kCompensated),
        ]
    )
    total, bad = _core.molecular_weight_mass(b"ABCD", table, flags)
    assert bad == -1
    assert total == sum(values) == -7.9

    #   The same sequence with the first float compensated as well, which is what
    #   a kernel without the phase would answer.
    sum_, correction = 0.0, 0.0
    for value in values[:2]:
        sum_ += value
    for value in values[2:]:
        next_ = sum_ + value
        if abs(sum_) >= abs(value):
            correction += (sum_ - next_) + value
        else:
            correction += (value - next_) + sum_
        sum_ = next_
    assert sum_ + correction == -7.8999999999999995
    assert sum_ + correction != total


def test_a_byte_the_sum_table_has_no_entry_for_is_a_decline():
    """The partial total comes back with the byte, so a caller can name it."""
    table, flags = raw_table(
        [
            (ord("A"), 1.0, _core.kCompensated),
            (ord("B"), 2.0, _core.kPlain),
        ]
    )
    assert _core.molecular_weight_mass(b"AB", table, flags) == (3.0, -1)
    total, bad = _core.molecular_weight_mass(b"ABZ", table, flags)
    assert bad == ord("Z")
    assert total == 3.0


def test_the_three_protein_kernels_decline_rather_than_raise():
    """`None` is a decline: only Python can spell the reference's message."""
    assert (
        _core.instability_index_sum(b"ACX", protparam._DIWV_TABLE, protparam._AA_MAP)
        is None
    )
    assert (
        _core.flexibility_scores(
            b"ACX" + b"A" * 9, protparam._FLEX_TABLE, protparam._AA_MAP
        )
        is None
    )
    table, valid, _ = protparam._dense_scale(SCALES_BY_NAME["kd"])
    weights = struct.pack("<1d", 1.0)
    assert _core.protein_scale_scores(b"ACX", 3, weights, table, valid, 3.0) is None
    assert _core.protein_scale_scores(b"ACA", 3, weights, table, valid, 3.0) == (
        spec_protein_scale("ACA", SCALES_BY_NAME["kd"], 3, 1.0)
    )
    #   `window < 2` is a decline too, because the reference raises there.
    assert _core.protein_scale_scores(b"ACX", 1, b"", table, valid, 1.0) is None
    assert _core.protein_scale_scores(b"ACX", 0, b"", table, valid, 1.0) is None


def test_the_instability_kernel_is_not_the_compensated_one():
    """`sum(...)` and `+=` are two different sums, and both live in one class.

    "EPCM" is the shortest sequence found whose three dipeptide weights come out
    differently the two ways: 47.32000000000001 added up one at a time, 47.32
    added up with compensation.  `instability_index` is the plain one.

    Which of the two `sum()` is is the interpreter's business and not this
    kernel's, so the contrast with it is drawn where it exists; the kernel's own
    answer is the naive one on every version, and that is the assertion the test
    is named for.
    """
    text = "EPCM"
    values = [DIWV[text[i]][text[i + 1]] for i in range(len(text) - 1)]
    naive = 0.0
    for value in values:
        naive += value
    total = _core.instability_index_sum(
        text.encode("ascii"), protparam._DIWV_TABLE, protparam._AA_MAP
    )
    assert total == naive
    if compensated_sum:
        assert naive != sum(values)
        assert total != sum(values)
    assert ProteinAnalysis(text).instability_index() == (10.0 / 4) * naive


# --------------------------------------------------------------------------
# The evidence: the same calls through Biopython.
# --------------------------------------------------------------------------


@needs_biopython
def test_the_protparam_data_namespace_is_the_reference_s():
    """Every name the reference exposes, with every value and every type."""
    from Bio.SeqUtils import ProtParamData as reference

    public = {name for name in dir(reference) if not name.startswith("_")}
    assert public <= set(dir(protparam_data))
    for name in sorted(public - {"gravy_scales", "DIWV"}):
        ours, theirs = getattr(protparam_data, name), getattr(reference, name)
        for letter, value in theirs.items():
            assert ours[letter] == value, (name, letter)
            assert type(ours[letter]) is type(value), (name, letter)
    assert set(DIWV) == set(reference.DIWV)
    for this, row in reference.DIWV.items():
        #   The rows are compared key by key and not in order: the reference's
        #   literals were written in their own order and this package's tables are
        #   in `AMINO_ACIDS`, which is the order a kernel indexes them by.
        assert set(DIWV[this]) == set(row)
        for following, value in row.items():
            assert DIWV[this][following] == value, (this, following)
            assert type(DIWV[this][following]) is type(value), (this, following)


@needs_biopython
def test_the_gravy_scales_keep_the_reference_s_int_and_float_apart():
    """The generated tables are the one place a scale's *type* could be lost.

    Value for value the two would agree either way; it is the type that a table
    built by coercion silently changes, and the type is what `sum()` reads.  How
    many scales are affected, and which of them it can show on, is the subject of
    `test_the_ten_int_bearing_scales_and_the_four_that_answer_differently`.
    """
    from Bio.SeqUtils import ProtParamData as reference

    assert list(gravy_scales) == list(reference.gravy_scales)
    for name, scale in reference.gravy_scales.items():
        ours = gravy_scales[name]
        #   Key by key: `ours` is in `AMINO_ACIDS` order and the reference's is
        #   in whatever order its literal was written in.  See the module
        #   docstring of `biofasting.protparam_data`.
        assert set(ours) == set(scale), name
        for letter, value in scale.items():
            assert ours[letter] == value, (name, letter)
            assert type(ours[letter]) is type(value), (name, letter)


@needs_biopython
def test_four_short_sequences_whose_sum_the_int_rule_changes():
    """The sum, and the compensated sum over the same numbers, which is not it.

    `flattened` is what a table that had coerced every value with `float()` would
    produce: same numbers, different summation.  `gravy` divides by the length
    afterwards, which can round the two together again -- it does for Roseman's
    "THTQG" -- so the comparison is made on the sum.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference
    from Bio.SeqUtils import ProtParamData as reference

    cases = (
        ("Parker", "IYAR"),
        ("Engelman", "HSV"),
        ("GoldSack", "TNRC"),
        ("Roseman", "THTQG"),
    )
    for name, text in cases:
        scale = reference.gravy_scales[name]
        assert any(type(value) is int for value in scale.values()), name
        total = sum(scale[aa] for aa in text)
        flattened = sum(float(scale[aa]) for aa in text)
        assert total != flattened, name
        assert Reference(text).gravy(name) == total / len(text), name
        assert ProteinAnalysis(text).gravy(name) == total / len(text), name


@needs_biopython
def test_the_ten_int_bearing_scales_and_the_four_that_answer_differently():
    """The count is the point: six of the ten can never tell the difference.

    A scale whose one integer value is `0` or `1` adds the same double either way,
    so nothing about the answer says whether it was rendered as `1` or as `1.0`.
    Ten scales carry integers and four of them can show it, which is why the
    differential layer needs more than a value-for-value comparison.
    """
    from Bio.SeqUtils import ProtParamData as reference
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    integers = {
        name
        for name, scale in reference.gravy_scales.items()
        if any(type(value) is int for value in scale.values())
    }
    assert integers == {
        "Engelman", "Fasman", "Fauchere", "GoldSack", "Jones",
        "Parker", "Ponnuswamy", "Roseman", "Wilson", "Zimmerman",
    }
    rng = random.Random(23)
    differing = set()
    for _ in range(4000):
        text = "".join(rng.choice(AMINO_ACIDS) for _ in range(rng.randrange(1, 60)))
        for name in sorted(integers):
            scale = reference.gravy_scales[name]
            if sum(scale[aa] for aa in text) != sum(float(scale[aa]) for aa in text):
                differing.add(name)
    assert differing == {"Engelman", "GoldSack", "Parker", "Roseman"}
    #   And ours agrees with the reference on a sweep that reaches all ten.
    for _ in range(200):
        text = "".join(rng.choice(AMINO_ACIDS) for _ in range(rng.randrange(1, 60)))
        ours, theirs = ProteinAnalysis(text), Reference(text)
        for name in sorted(integers):
            assert ours.gravy(name) == theirs.gravy(name), (name, text)


@needs_biopython
def test_everything_matches_biopython_over_random_proteins():
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    rng = random.Random(2026)
    for length in (0, 1, 2, 3, 8, 9, 10, 11, 40, 199, 500, 1024):
        for _ in range(3):
            text = "".join(rng.choice(AMINO_ACIDS) for _ in range(length))
            ours, theirs = ProteinAnalysis(text), Reference(text)
            assert ours.sequence == theirs.sequence
            assert ours.length == theirs.length
            assert ours.count_amino_acids() == theirs.count_amino_acids()
            for ours_call, theirs_call in (
                (lambda: ours.instability_index(), lambda: theirs.instability_index()),
                (lambda: ours.flexibility(), lambda: theirs.flexibility()),
                (lambda: ours.gravy(), lambda: theirs.gravy()),
                (lambda: ours.aromaticity(), lambda: theirs.aromaticity()),
                (lambda: ours.molecular_weight(), lambda: theirs.molecular_weight()),
                (
                    lambda: ours.secondary_structure_fraction(),
                    lambda: theirs.secondary_structure_fraction(),
                ),
                (
                    lambda: ours.molar_extinction_coefficient(),
                    lambda: theirs.molar_extinction_coefficient(),
                ),
                (lambda: ours.isoelectric_point(), lambda: theirs.isoelectric_point()),
                (
                    lambda: ours.amino_acids_percent,
                    lambda: theirs.amino_acids_percent,
                ),
            ):
                assert call(ours_call) == call(theirs_call), (length, text)
            for pH in (1.0, 4.05, 7.0, 7.775, 12.0):
                assert call(lambda: ours.charge_at_pH(pH)) == call(
                    lambda: theirs.charge_at_pH(pH)
                ), (length, text, pH)


@needs_biopython
def test_gravy_matches_biopython_over_every_scale():
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    rng = random.Random(77)
    for length in (1, 2, 3, 5, 17, 129, 500, 2048):
        text = "".join(rng.choice(AMINO_ACIDS) for _ in range(length))
        ours, theirs = ProteinAnalysis(text), Reference(text)
        for name in gravy_scales:
            assert call(lambda: ours.gravy(name)) == call(
                lambda: theirs.gravy(name)
            ), (name, length, text)
    assert call(lambda: ProteinAnalysis("ACDE").gravy("NotAScale")) == call(
        lambda: Reference("ACDE").gravy("NotAScale")
    )


@needs_biopython
def test_protein_scale_matches_biopython_over_every_scale_and_window():
    from Bio.SeqUtils import ProtParamData as reference
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    rng = random.Random(5150)
    text = "".join(rng.choice(AMINO_ACIDS) for _ in range(257))
    ours, theirs = ProteinAnalysis(text), Reference(text)
    for name, scale in SCALES_BY_NAME.items():
        other = getattr(reference, name, None)
        if other is None:
            other = reference.cw[float(name.removeprefix("Cowan"))]
        assert scale == other, name
        for window, edge in ((2, 1.0), (3, 0.5), (4, 0.0), (9, 1.0), (11, 0.4), (13, 2.5)):
            assert call(lambda: ours.protein_scale(scale, window, edge)) == call(
                lambda: theirs.protein_scale(other, window, edge)
            ), (name, window, edge)
    #   A scale with no entry for most of the residues, which is the path the
    #   kernel declines on: the warning, the skipped term, and the score.
    borrowed = {"A": 1.0, "C": 2.0}
    ours_noise = io.StringIO()
    with contextlib.redirect_stderr(ours_noise):
        ours_scores = ours.protein_scale(borrowed, 5, 1.0)
    theirs_noise = io.StringIO()
    with contextlib.redirect_stderr(theirs_noise):
        theirs_scores = theirs.protein_scale(borrowed, 5, 1.0)
    assert ours_scores == theirs_scores
    assert ours_noise.getvalue() == theirs_noise.getvalue()
    assert ours_noise.getvalue().count("not a standard amino acid") > 100


@needs_biopython
def test_the_raising_paths_raise_the_same_exception_with_the_same_arguments():
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    cases = (
        (lambda p: p.instability_index(), "ACX"),
        (lambda p: p.instability_index(), "AXC"),
        (lambda p: p.instability_index(), ""),
        (lambda p: p.flexibility(), "A" * 9 + "X" + "A" * 9),
        (lambda p: p.gravy(), ""),
        (lambda p: p.gravy("NotAScale"), "ACDE"),
        (lambda p: p.amino_acids_percent, ""),
        (lambda p: p.protein_scale(SCALES_BY_NAME["kd"], 1), "ACDEFG"),
        (lambda p: p.protein_scale(SCALES_BY_NAME["kd"], 0), "ACDEFG"),
    )
    for case, text in cases:
        ours = call(lambda: case(ProteinAnalysis(text)))
        theirs = call(lambda: case(Reference(text)))
        assert ours == theirs, (text, ours, theirs)
        assert ours[0] == "raised", text


@needs_biopython
def test_a_character_a_byte_cannot_hold_takes_the_reference_s_own_loop():
    """Latin-1 is a byte and declines the residue map; above it is not a byte."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    for text in ("ACDéFG", "ACD€FG"):
        ours, theirs = ProteinAnalysis(text), Reference(text)
        assert ours.sequence == theirs.sequence
        assert ours.count_amino_acids() == theirs.count_amino_acids()
        assert call(ours.instability_index) == call(theirs.instability_index)
        assert call(ours.flexibility) == call(theirs.flexibility)
        assert call(lambda: ours.gravy("Parker")) == call(
            lambda: theirs.gravy("Parker")
        )
        assert call(lambda: ours.protein_scale(SCALES_BY_NAME["kd"], 3)) == call(
            lambda: theirs.protein_scale(SCALES_BY_NAME["kd"], 3)
        )
        assert call(ours.instability_index)[0] == "raised"


@needs_biopython
def test_count_amino_acids_matches_biopython_on_the_residues_it_ignores():
    """The kernel's buckets, over sequences that are mostly *not* the twenty.

    The random-protein sweep above feeds it nothing but standard residues, which
    is the one input where "count the twenty and drop the rest" cannot be told
    from "count everything and hope".  `X`, `B`, `Z`, `J`, `U`, `O`, `*`, a digit,
    a space and a lower-case letter are all counted by nobody, and `str.count`
    says so by returning zero rather than raising.  The long sequence is here so
    that the kernel runs and not a fallback.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis as Reference

    rng = random.Random(5150)
    alphabet = AMINO_ACIDS + "XBZJUO*x1 \t" + AMINO_ACIDS.lower()
    for length in (0, 1, 7, 300, 5000):
        for _ in range(3):
            text = "".join(rng.choice(alphabet) for _ in range(length))
            ours = ProteinAnalysis(text).count_amino_acids()
            theirs = Reference(text).count_amino_acids()
            assert ours == theirs, (length, text[:40])
            assert list(ours) == list(theirs)


@needs_biopython
def test_the_docstring_examples_are_true():
    """Every `>>>` in the three modules is executed, not merely written."""
    for module in (protparam, isoelectric_point, protparam_data):
        results = doctest.testmod(module, verbose=False, report=False)
        assert results.attempted > 0, f"{module.__name__} documents nothing"
        assert results.failed == 0, module.__name__
