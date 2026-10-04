"""PLAN 2.1: the pairwise aligner, and the accelerator under it.

`bench/rank_align.py` measured the target before this module existed: over six
rows, gate first, `biofasting.alignment` is 6x-35x `Bio.Align.PairwiseAligner`
score-only and 9x-38x with the traceback, because the reference is a scalar DP
and parasail is SIMD.  What this file is about is the other half of that
sentence -- a speedup is only worth having if the answer is the same answer, so
almost everything here compares against the reference rather than against an
expectation written by hand.

Three groups, in the order they can run:

* **The refusals and the packaging**, which need no accelerator at all.  A
  scheme this module cannot express is refused for its own reason whether or not
  parasail is installed, and the absent-parasail error names the extra -- both
  are pinned here, the second in a fresh interpreter.
* **The scores**, checked against the reference scheme by scheme, including the
  end-gap mapping that is the one place Biopython expresses more than parasail
  can and the `sg` correction that the mapping exposed.
* **The paths**, checked against the reference's own optima: our path must be
  one of them, its recomputed score must equal the kernel's, and where the
  reference has exactly one optimum the two must agree column for column.

The accelerator is optional (`pip install biofasting[alignment]`) because
parasail 1.3.4 ships no aarch64 wheel, so the groups that need it are skipped
rather than failing on a machine without it -- and the first group then says
something about a bare install, which is the state most of this package's users
are in.
"""

import importlib.util
import random
import subprocess
import sys

import numpy
import pytest

from biofasting import alignment

needs_parasail = pytest.mark.skipif(
    importlib.util.find_spec("parasail") is None,
    reason="the optional parasail accelerator is not installed",
)
needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

DNA = "ACGT"

# One DNA scheme, written once: match/mismatch 2/-2 with a 10+1k affine gap, and
# the terminal gaps priced exactly like the interior ones -- which is `nw`, and
# is the only global scheme other than a fully free-ended one that this module
# accepts.
NW = dict(
    mode="global",
    match_score=2,
    mismatch_score=-2,
    open_gap_score=-10,
    extend_gap_score=-1,
    end_gap_score=-10,
    extend_end_gap_score=-1,
)
SEMIGLOBAL = dict(NW, end_gap_score=0, extend_end_gap_score=0)
LOCAL = dict(
    mode="local", match_score=2, mismatch_score=-2, open_gap_score=-10, extend_gap_score=-1
)

PROTEIN = dict(
    mode="global", open_gap_score=-11, extend_gap_score=-1,
    end_gap_score=-11, extend_end_gap_score=-1,
)


def reference(attributes):
    """The reference aligner for one of the schemes above.

    Two attributes are passed through as they are -- the mode is a string and
    the matrix is an object -- and every score becomes a float, because that is
    what the reference stores.
    """
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    for name, value in attributes.items():
        passthrough = name in ("mode", "substitution_matrix")
        setattr(aligner, name, value if passthrough else float(value))
    return aligner


def random_pair(rng, longest):
    return (
        "".join(rng.choice(DNA) for _ in range(rng.randrange(1, longest))),
        "".join(rng.choice(DNA) for _ in range(rng.randrange(1, longest))),
    )


def run_in_a_fresh_interpreter(code):
    """Run ``code`` and return the lines ``print`` put on stdout.

    The *last* lines are the ones this test wrote: an editable install rebuilds
    the extension on import and prints cmake's progress to the same stream, so
    everything above our own output is not ours to parse.
    """
    done = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    return done.stdout.strip().splitlines()


# --------------------------------------------------------------------------
# What cannot be expressed, and what it says instead.
# --------------------------------------------------------------------------


def test_a_fractional_score_is_refused_and_names_the_reference():
    """parasail scores in integer registers; the reference does not."""
    with pytest.raises(NotImplementedError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, match_score=2.5))
    assert "match_score=2.5" in str(raised.value)
    assert "Bio.Align.PairwiseAligner" in str(raised.value)


def test_a_fractional_gap_score_is_refused_too():
    with pytest.raises(NotImplementedError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, extend_gap_score=-1.5))
    assert "extend_gap_score=-1.5" in str(raised.value)


def test_a_positive_gap_score_is_refused():
    """A positive gap price is a reward in the reference, not a penalty."""
    with pytest.raises(NotImplementedError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, open_gap_score=10))
    assert "reward" in str(raised.value)


def test_end_gaps_between_the_two_schemes_are_refused_and_the_message_names_both():
    """The one place the reference expresses more than parasail can.

    A terminal gap priced differently from an interior one is neither `nw` (they
    are equal) nor `sg` (both are zero), and there is no third kernel.  The
    refusal has to name the two that exist, because the caller's fix is to pick
    one of them.
    """
    with pytest.raises(NotImplementedError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, end_gap_score=-5))
    message = str(raised.value)
    assert "end_gap_score=-5" in message
    assert "`nw`" in message and "`sg`" in message


def test_setting_only_the_interior_gaps_is_refused_because_the_reference_would_differ():
    """The natural mistake, and the reason the refusal is not a nuisance.

    `PairwiseAligner()`'s default `end_gap_score` is -1, so a caller who moves
    only the interior prices has quietly asked the reference for a scheme whose
    terminal gaps cost 1 while their interior gaps cost 10.  That scheme has no
    parasail kernel, and answering with `nw` would be answering a question the
    caller did not ask.
    """
    with pytest.raises(NotImplementedError):
        alignment.score("ACGT", "ACGT", mode="global", match_score=2, mismatch_score=-2,
                        open_gap_score=-10, extend_gap_score=-1)


@needs_parasail
def test_the_two_schemes_that_are_expressible_are_accepted():
    assert alignment.score("ACGT", "ACGT", **NW) == 8
    assert alignment.score("ACGT", "ACGT", **SEMIGLOBAL) == 8
    assert alignment.score("ACGT", "ACGT", **LOCAL) == 8


def test_fogsaa_is_refused_by_name():
    """Biopython's other global algorithm has no counterpart at all."""
    with pytest.raises(NotImplementedError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, mode="fogsaa"))
    assert "fogsaa" in str(raised.value)


def test_an_unknown_mode_is_a_value_error_that_names_the_two_that_work():
    with pytest.raises(ValueError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, mode="glocal"))
    assert "'global'" in str(raised.value) and "'local'" in str(raised.value)


def test_an_empty_sequence_is_the_references_own_error():
    """Not parasail's: it would print to the C stderr and return NULL."""
    for first, second in (("", "ACGT"), ("ACGT", "")):
        with pytest.raises(ValueError) as raised:
            alignment.score(first, second, **NW)
        assert str(raised.value) == "sequence has zero length"


def test_a_letter_parasail_cannot_carry_is_refused():
    """Above ASCII is what parasail's C-string alphabet cannot hold.

    A Latin-1 letter is refused rather than passed through, and a character
    above U+00FF is refused before that, because it would be split into bytes
    rather than refused at all.
    """
    with pytest.raises(ValueError) as raised:
        alignment.score("é", "ACGT", **NW)
    assert "ASCII" in str(raised.value)
    with pytest.raises(ValueError) as raised:
        alignment.score("中", "ACGT", **NW)
    assert "U+00FF" in str(raised.value)


@needs_parasail
def test_a_byte_string_is_read_as_latin_1():
    """This package's convention for a measurement, and the reference's."""
    assert alignment.score(b"ACGT", bytearray(b"ACGT"), **NW) == 8


def test_a_substitution_matrix_must_be_an_array_with_an_alphabet():
    """A bare list cannot say which letters are unknown, so it is refused."""
    with pytest.raises(TypeError) as raised:
        alignment.score("ACGT", "ACGT", **dict(NW, substitution_matrix=[[1, 2], [2, 1]]))
    assert "alphabet" in str(raised.value)


def test_a_letter_the_matrix_does_not_know_is_refused_not_scored_zero():
    """The trap this check exists for.

    `parasail.matrix_create` puts every unknown letter in one extra slot scoring
    0 against everything, so an unknown letter is silent and wrong instead of
    loud; the reference raises for it, and so does this module.
    """
    from biofasting import substitution_matrices

    matrix = substitution_matrices.load("BLOSUM62")
    with pytest.raises(ValueError) as raised:
        alignment.score("J", "J", **dict(PROTEIN, substitution_matrix=matrix))
    assert "not in the alphabet" in str(raised.value)


def test_available_reports_a_bool_without_importing_anything():
    """The readers never ask; a caller choosing between two tools may."""
    assert isinstance(alignment.available(), bool)


# --------------------------------------------------------------------------
# The optional accelerator stays optional.
# --------------------------------------------------------------------------


def test_importing_the_package_imports_neither_parasail_nor_biopython():
    """Run in a fresh interpreter, with both installed.

    `biofasting.alignment` is imported by `biofasting/__init__.py`, so this is
    the test that keeps the optional extra optional: an import at module level
    would make parasail a hard dependency of every install, which is exactly
    what the missing aarch64 wheel makes impossible.
    """
    code = (
        "import sys, biofasting\n"
        "assert biofasting.Aligner is not None\n"
        "print('parasail' in sys.modules, 'Bio' in sys.modules)\n"
    )
    assert run_in_a_fresh_interpreter(code)[-1] == "False False"


def test_the_missing_parasail_error_names_the_extra():
    """The absent-parasail path, faked rather than skipped.

    ``sys.modules['parasail'] = None`` makes any import of it raise, which is
    the state a wheel install without the extra is in -- so this test runs
    everywhere, including on the machines that skip every timing test.  What it
    pins is that the failure is ours and names the extra, not a bare
    ``ModuleNotFoundError`` from inside a library the caller did not know
    wanted one.
    """
    code = (
        "import sys\n"
        "sys.modules['parasail'] = None\n"
        "from biofasting import alignment\n"
        "print('available', alignment.available())\n"
        "try:\n"
        "    alignment.score('ACGT', 'ACGT', match_score=2, mismatch_score=-2,\n"
        "                    open_gap_score=-10, extend_gap_score=-1,\n"
        "                    end_gap_score=-10, extend_end_gap_score=-1)\n"
        "except ImportError as exc:\n"
        "    print(exc)\n"
        "else:\n"
        "    print('NO ERROR')\n"
    )
    done = run_in_a_fresh_interpreter(code)
    assert done[-2] == "available False"
    assert "pip install biofasting[alignment]" in done[-1]


def test_a_refused_scheme_is_refused_even_without_the_accelerator():
    """The order the checks run in, pinned.

    The scheme is validated before parasail is imported, so the reason a caller
    is given is the reason their scheme is wrong -- not "install the extra" for
    a scheme that would have been refused anyway.
    """
    code = (
        "import sys\n"
        "sys.modules['parasail'] = None\n"
        "from biofasting import alignment\n"
        "try:\n"
        "    alignment.score('ACGT', 'ACGT', match_score=2.5, mismatch_score=-2,\n"
        "                    open_gap_score=-10, extend_gap_score=-1,\n"
        "                    end_gap_score=-10, extend_end_gap_score=-1)\n"
        "except NotImplementedError as exc:\n"
        "    print('refused:', exc)\n"
        "else:\n"
        "    print('NO ERROR')\n"
    )
    done = run_in_a_fresh_interpreter(code)
    assert done[-1].startswith("refused:")
    assert "match_score=2.5" in done[-1]


# --------------------------------------------------------------------------
# The scores.
# --------------------------------------------------------------------------


@needs_parasail
def test_the_simple_cases_score_what_they_should():
    assert alignment.score("ACGT", "ACGT", **NW) == 8
    # -4 and not -8: the last letters are both T, so three columns mismatch and
    # one matches.  Measured on both sides rather than reasoned about, and kept
    # because it is the case that catches a scheme whose mismatch price has
    # drifted.
    assert alignment.score("ACGT", "TTTT", **NW) == -4
    # Four mismatches, not four matches: the reference is case-sensitive and
    # this module is too, which is what `case_sensitive=True` buys.
    assert alignment.score("acgt", "ACGT", **NW) == -8
    # '-' is an ordinary letter to both, and matches itself.
    assert alignment.score("AC-GT", "AC-GT", **NW) == 10
    assert alignment.score("AC-GT", "ACGT", **NW) == -2


@needs_parasail
def test_the_three_way_end_gap_mapping_is_exactly_the_references():
    """Equal prices are `nw`, zero prices are `sg`, and nothing else exists.

    The differential is over random pairs in both schemes, because the mapping
    is the whole justification for refusing the third case and a single example
    would not show it.  The last count is the other half: the two schemes must
    *disagree* often, or a module that ignored the end gaps would pass both
    loops.
    """
    rng = random.Random(101)
    ours_nw, ours_sg = alignment.Aligner(**NW), alignment.Aligner(**SEMIGLOBAL)
    different = 0
    for _ in range(200):
        first, second = random_pair(rng, 16)
        assert ours_nw.score(first, second) == alignment.score(first, second, **NW)
        assert ours_sg.score(first, second) == alignment.score(first, second, **SEMIGLOBAL)
        if ours_nw.score(first, second) != ours_sg.score(first, second):
            different += 1
    assert different > 20


@needs_parasail
def test_local_mode_ignores_the_end_gap_scores():
    """Measured, 0 of 63 pairs, so `sw` is the kernel whatever they say.

    Also the reason the constructor accepts a local scheme whose end gaps would
    be refused in global mode: they are not part of the question being asked.
    """
    pairs = [("TTTACGTACGTGGG", "ACGTACGT"), ("ACGTACGT", "ACGTACGT")]
    for first, second in pairs:
        scores = {
            alignment.score(first, second, **dict(LOCAL, end_gap_score=end, extend_end_gap_score=extend))
            for end, extend in ((0, 0), (-10, -1), (-7, -3))
        }
        assert len(scores) == 1


@needs_parasail
def test_a_free_end_alignment_is_never_negative_and_is_sometimes_zero():
    """The two visible consequences of the `max(sg, 0)` correction.

    With both ends free an alignment never has to score below zero -- laying one
    sequence entirely against the other is free -- so a negative score from this
    module would be the uncorrected kernel leaking through, and a sweep in which
    nothing ever reaches zero would mean the correction is never exercised.
    """
    rng = random.Random(7)
    ours = alignment.Aligner(**SEMIGLOBAL)
    zeros = 0
    for _ in range(400):
        first, second = random_pair(rng, 10)
        score = ours.score(first, second)
        assert score >= 0
        if score == 0:
            zeros += 1
            assert ours.align(first, second).score == 0
    assert zeros > 5


@needs_parasail
def test_the_no_column_alignment_is_returned_as_the_reference_shapes_it():
    """'AAAA' against 'TTTT' with both ends free: score 0, and *this* path.

    parasail says -2 there, because an insertion run adjacent to a deletion run
    is a path its recurrence does not have; the correction returns the
    alignment the reference reports -- every letter in a free gap, and no column
    pairing two letters.  The coordinates are that alignment in this module's
    finest segmentation, one segment per run; the reference merges the second
    run's zero-length step away and writes ``[[0, 4, 4], [0, 0, 4]]``, which
    renders identically (`test_the_no_column_alignment_renders_as_the_references_own`).
    """
    ours = alignment.Aligner(**SEMIGLOBAL)
    assert ours.score("AAAA", "TTTT") == 0
    result = ours.align("AAAA", "TTTT")
    assert result.score == 0
    assert result[0] == "AAAA----"
    assert result[1] == "----TTTT"
    assert result.coordinates.tolist() == [[0, 4, 4, 4], [0, 0, 0, 4]]
    counts = result.counts()
    assert (counts.identities, counts.mismatches, counts.gaps) == (0, 0, 8)
    assert (counts.substitution_score, counts.gap_score) == (0, 0)
    assert counts.score == 0


@needs_parasail
@needs_biopython
def test_the_no_column_alignment_renders_as_the_references_own():
    """And it is one of the reference's two optima, not an invention.

    The reference has two equally good answers here -- the same two runs in the
    other order -- and this module returns the second of them.  Which one is the
    tie-break and not the contract; that it is one of them is the contract.
    """
    reference_aligner = reference(SEMIGLOBAL)
    optima = list(reference_aligner.align("AAAA", "TTTT"))
    assert len(optima) == 2
    rendering = alignment.align("AAAA", "TTTT", **SEMIGLOBAL).to_biopython().format("fasta")
    assert rendering == optima[1].format("fasta")
    assert [optimum.score for optimum in optima] == [0, 0]


@needs_parasail
def test_a_local_alignment_with_nothing_to_align_is_the_empty_one():
    """The divergence that is a shape, not a score.

    With no letter pair scoring above zero the reference's `score` says 0.0 and
    its `align` yields *nothing*; `sw` reports a zero-length alignment, so this
    module returns an explicit empty `Alignment` with score 0.  The score is the
    reference's; the object is not the reference's empty iterator, and that is
    the honest way to say "there is nothing here".
    """
    ours = alignment.Aligner(**LOCAL)
    assert ours.score("AAAA", "TTTT") == 0
    result = ours.align("AAAA", "TTTT")
    assert result.score == 0
    assert (result[0], result[1]) == ("", "")
    assert result.coordinates.tolist() == [[], []]
    assert result.counts().gaps == 0


@needs_parasail
def test_a_gapped_alignments_counts_are_the_arithmetic_they_claim():
    """One hand-checkable case, so the counts are not only compared to a computer.

    `'ACGTACGT'` against `'ACGT'` globally: four columns of 2, then a four-long
    terminal gap costing 10 + 3*1 = 13, so 8 - 13 = -5.

    Where that gap *is* is a tie-break, and the two implementations break it
    differently: the reference reports `'ACGT----'` and parasail's DP reports
    `'----ACGT'`, the same alignment moved to the other end.  The scores and the
    counts are identical -- a terminal gap is a terminal gap -- which is the
    reason the contract claims optimality and not path identity.
    """
    result = alignment.align("ACGTACGT", "ACGT", **NW)
    assert result.score == -5
    assert result[1] in ("ACGT----", "----ACGT")
    assert result[0].replace("-", "") == "ACGTACGT"
    counts = result.counts()
    assert (counts.identities, counts.mismatches, counts.gaps) == (4, 0, 4)
    assert counts.substitution_score == 8
    assert counts.gap_score == -13
    assert counts.score == result.score


@needs_parasail
def test_a_path_is_always_optimal_by_its_own_arithmetic():
    """The independent check, on every pair of the sweep.

    `counts().score` recomputes the path's score from the aligned letters and
    the matrix, which is a second implementation of the scoring and not a second
    call into parasail -- so it can disagree with the kernel, and does not.
    """
    rng = random.Random(11)
    for attributes in (NW, SEMIGLOBAL, LOCAL):
        ours = alignment.Aligner(**attributes)
        for _ in range(150):
            first, second = random_pair(rng, 14)
            result = ours.align(first, second)
            assert result.counts().score == result.score
            assert len(result[0]) == len(result[1])
            # The gapped pair consumes the inputs, and a local alignment
            # consumes the segment of them it selected.
            if attributes is LOCAL:
                assert result[0].replace("-", "") in first
                assert result[1].replace("-", "") in second
            else:
                assert result[0].replace("-", "") == first
                assert result[1].replace("-", "") == second


@needs_parasail
def test_a_saturated_register_is_an_error_and_not_a_wrong_score():
    """The guard, reached by forcing the width it exists to prevent.

    The width normally comes from the score bound, so this is unreachable -- and
    the point of the test is that it stays an error if the bound is ever
    loosened: parasail saturates silently, returning 0 for a pair worth 400.
    """
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(alignment, "_width_for", lambda bound: 8)
        ours = alignment.Aligner(**NW)
        with pytest.raises(RuntimeError) as raised:
            ours.score("A" * 200, "A" * 200)
        assert "saturated" in str(raised.value)


@needs_parasail
def test_an_alignment_holds_two_gapped_sequences_and_says_so():
    result = alignment.align("ACGTACGT", "ACGT", **NW)
    assert len(result) == 2
    assert result.sequences == ("ACGTACGT", "ACGT")
    coordinates = result.coordinates
    assert coordinates.dtype == numpy.intp
    assert coordinates.shape[0] == 2
    # One segment per CIGAR run, so two columns per run and not `len + 1`: the
    # reference merges a run into its neighbour where this keeps them apart.
    assert coordinates.shape[1] % 2 == 0
    assert coordinates[0][0] == coordinates[1][0] == 0
    assert coordinates[0][-1] == len(result.sequences[0])
    assert coordinates[1][-1] == len(result.sequences[1])
    assert (numpy.diff(coordinates, axis=1) >= 0).all()
    with pytest.raises(IndexError):
        result[2]
    assert "score=-5" in repr(result)


# --------------------------------------------------------------------------
# Against the reference, score for score and path for path.
# --------------------------------------------------------------------------


@needs_parasail
@needs_biopython
def test_the_scores_are_the_references_over_random_pairs():
    """The contract, swept: 3 schemes x 300 pairs, gate first."""
    rng = random.Random(23)
    for attributes in (NW, SEMIGLOBAL, LOCAL):
        ours, theirs = alignment.Aligner(**attributes), reference(attributes)
        for _ in range(300):
            first, second = random_pair(rng, 16)
            assert ours.score(first, second) == theirs.score(first, second)


@needs_parasail
@needs_biopython
def test_the_semiglobal_correction_is_the_references_score_exactly():
    """Including the pairs where parasail's own answer is negative.

    This is the assertion the `max(sg, 0)` correction rests on: every pair where
    the reference and a raw `sg` disagree is a pair where the optimum has no
    substitution column.
    """
    import parasail

    rng = random.Random(41)
    matrix = parasail.matrix_create(DNA, 2, -2)
    ours, theirs = alignment.Aligner(**SEMIGLOBAL), reference(SEMIGLOBAL)
    raw_negative = 0
    for _ in range(500):
        first, second = random_pair(rng, 12)
        raw = parasail.sg_scan_16(first, second, 10, 1, matrix).score
        if raw < 0:
            raw_negative += 1
            # A negative raw score means the optimum is the no-column
            # alignment, and `max` is the whole correction.
            assert theirs.score(first, second) == 0
        assert ours.score(first, second) == theirs.score(first, second)
    assert raw_negative > 0


@needs_parasail
@needs_biopython
def test_our_path_is_one_of_the_references_optima():
    """Not the same path -- one of the optimal ones, and counts agree.

    The reference yields every optimum in its own DP's order; this module
    returns one.  So the claim that can be made is membership, and it is made
    here: for every pair of the sweep, our rendering equals one of the
    reference's, and that one's `identities`/`mismatches`/`gaps` equal ours.
    """
    rng = random.Random(59)
    checked = 0
    for attributes in (NW, SEMIGLOBAL, LOCAL):
        ours, theirs = alignment.Aligner(**attributes), reference(attributes)
        for _ in range(60):
            first, second = random_pair(rng, 12)
            if theirs.score(first, second) != ours.score(first, second):
                pytest.fail(f"score mismatch on {first!r}/{second!r}")
            optima = theirs.align(first, second)
            if len(optima) == 0:
                continue  # the empty-optimum case, tested above
            mine = ours.align(first, second)
            for index, optimum in enumerate(optima):
                if index >= 50:
                    pytest.fail(f"too many optima to scan: {first!r}/{second!r}")
                if optimum.format("fasta") != mine.to_biopython().format("fasta"):
                    continue
                counts, my_counts = optimum.counts(), mine.counts()
                assert (counts.identities, counts.mismatches, counts.gaps) == (
                    my_counts.identities, my_counts.mismatches, my_counts.gaps
                )
                checked += 1
                break
            else:
                pytest.fail(f"our path is not one of the optima: {first!r}/{second!r}")
    assert checked > 100


@needs_parasail
@needs_biopython
def test_a_unique_optimum_is_the_same_alignment_column_for_column():
    """Where the reference has one answer there is nothing to tie-break.

    A pair with one optimum is one where every optimal path is the same path, so
    the two implementations cannot differ for a legitimate reason -- and the
    comparison is the reference's own `format("fasta")` output, not a
    reimplementation of it.
    """
    rng = random.Random(83)
    seen = 0
    for attributes in (NW, LOCAL):
        ours, theirs = alignment.Aligner(**attributes), reference(attributes)
        for _ in range(150):
            first, second = random_pair(rng, 10)
            optima = theirs.align(first, second)
            if len(optima) != 1:
                continue
            seen += 1
            expected = next(iter(optima))
            assert ours.align(first, second).to_biopython().format("fasta") == expected.format(
                "fasta"
            )
    assert seen > 20


@needs_parasail
@needs_biopython
def test_a_protein_scheme_scores_the_same_with_either_blosum62():
    """Our own matrix data, and Biopython's, must be the same table here.

    `biofasting.substitution_matrices` exists so that a program without
    Biopython still has a scoring matrix; this is the case where the two meet,
    and it is the reason `_read_table` reads the values out of whichever `Array` it
    is handed rather than out of a table of its own.
    """
    from biofasting import substitution_matrices

    ours = alignment.Aligner(substitution_matrix=substitution_matrices.load("BLOSUM62"), **PROTEIN)
    from Bio.Align import substitution_matrices as theirs_matrices

    theirs = reference(dict(PROTEIN, substitution_matrix=theirs_matrices.load("BLOSUM62")))
    letters = "ACDEFGHIKLMNPQRSTVWY"
    rng = random.Random(97)
    for _ in range(50):
        first = "".join(rng.choice(letters) for _ in range(rng.randrange(5, 30)))
        second = "".join(rng.choice(letters) for _ in range(rng.randrange(5, 30)))
        assert ours.score(first, second) == theirs.score(first, second)


@needs_parasail
@needs_biopython
def test_the_conversion_is_the_references_own_object():
    """`to_biopython` returns a `Bio.Align.Alignment`, not a look-alike.

    The point of the conversion is that everything Biopython offers on an
    alignment -- `map`, `substitutions`, the format strings, slicing -- comes
    from Biopython rather than from a reimplementation here, so the test is
    that the object is the reference's type and behaves like it.
    """
    from Bio.Align import Alignment as ReferenceAlignment

    ours = alignment.align("ACGTACGT", "ACGTACGTAAA", **NW)
    converted = ours.to_biopython()
    assert isinstance(converted, ReferenceAlignment)
    assert converted[0] == ours[0]
    assert converted[1] == ours[1]
    counts, mine = converted.counts(), ours.counts()
    assert (counts.identities, counts.mismatches, counts.gaps) == (
        mine.identities, mine.mismatches, mine.gaps
    )
    # And the reference can build the same alignment from our coordinates.
    assert ReferenceAlignment(
        [ours.sequences[0], ours.sequences[1]], ours.coordinates
    ).format("fasta") == converted.format("fasta")


@needs_parasail
@needs_biopython
def test_the_convenience_functions_are_the_constructor_spelled_out():
    assert alignment.score("ACGT", "ACGT", **NW) == alignment.Aligner(**NW).score("ACGT", "ACGT")
    assert alignment.align("ACGT", "ACGT", **NW).score == 8
