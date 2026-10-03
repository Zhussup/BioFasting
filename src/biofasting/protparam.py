"""`Bio.SeqUtils.ProtParam`: `ProteinAnalysis`, and the three kernels under it.

Thirteen methods over one protein sequence, and the only ones that cost anything
are the three that walk every residue: `protein_scale`, `flexibility` and
`instability_index`.  `bench/targets.md` ranks them at 445 µs, 377 µs and 90 µs
on a 1,024-residue protein against a reference that spends the same time doing
the same arithmetic in Python, and records that the best pure-Python rewrite of
any of them is **2.0×** -- which is the measurement that says a kernel is the
only thing here that changes the order of magnitude.  So those three are
compiled, and everything else in the class is transcribed.

**The three are not one function.**  They are three different accumulations of
the same kind of term, and the difference is in the last bit of the answer:

* `protein_scale` adds ``weights[j] * front + weights[j] * back``;
* `flexibility` adds ``(front + back) * weights[j]``, reads its middle at index
  5 of a nine-residue window rather than at 4, and produces one window fewer
  than fit;
* `instability_index` adds with a plain ``score += value`` where the other two
  -- and `gravy` -- are written ``sum(...)`` and are therefore
  Neumaier-compensated by CPython 3.12 and later.

Distributing a weight over a pair, or summing a sequence with one summation
instead of the other, gives the same number on paper and a different double in
IEEE-754.  Measured, `flexibility` differs from the reference in the last bit on
200 of 200 proteins when the weight is distributed.  Each kernel therefore
reproduces its own reference's order of operations, and none of them is a
special case of another.

**A kernel counts, Python decides.**  Every one of the three returns `None` when
the reference would do something other than compute: `flexibility` and
`instability_index` raise `KeyError` naming the residue -- and the residue named
is the first one *read*, which is not the first one in the sequence -- while
`protein_scale` does not raise at all, it writes a line to stderr and drops that
window's term.  A compiled kernel cannot spell either.  So on such an input the
call re-runs the reference's own loop, here, in Python, and the reference's
behaviour comes out unchanged: same exception, same message, same warning, same
partial score.

Two preconditions of `protein_scale` are the reference's and are *not* repaired
here, because a port that raised a friendlier error would be a port that
answered differently:

* ``window == 1`` divides by ``window - 1`` while building the weights, so it
  raises before any loop runs -- the kernel declines at ``window < 2``;
* ``window == 0`` gets past that, then reads ``subsequence[0]`` of an empty
  slice.

The input convention is :mod:`biofasting.sequtils`'s: text, with a byte buffer
read as Latin-1, so that a `Seq` behaves as it does through the reference.  The
reference calls `.upper()` on the object it is given, which is a `str` method
that a `Seq` also has; a bare `bytes` has one too, but the `.count(aa)` after it
does not, so the reference raises `TypeError` where this reads the bytes as
text.
"""

from __future__ import annotations

import functools
import struct
import sys

from . import _core
from ._protparam_data import AMINO_ACIDS, DIWV as _DIWV_ROWS, FLEX
from .isoelectric_point import IsoelectricPoint
from .protparam_data import DIWV, Flex, SCALES_BY_NAME, gravy_scales
from .sequtils import _as_text, molecular_weight

__all__ = ["ProteinAnalysis"]


# ---------------------------------------------------------------------------
# The dense tables the kernels read.
#
# A scale is a mapping from a letter to a float, and a kernel wants an array
# indexed by byte.  Both are the same information; only the Python fallback
# wants the mapping, and only for inputs the kernel declined, which are exactly
# the inputs where someone will be reading a traceback rather than a profile.


def _dense_scale(scale):
    """A mapping as dense arrays for a kernel, or `None` if it will not fit.

    Three of them, because `protein_scale` asks whether a letter is there and
    `gravy` asks what it is: `valid` is 0 or 1, and `sum_flags` is 0 for absent,
    1 for a value the reference wrote as a float and 2 for one it wrote as an
    int.  The second distinction is CPython's, not this package's -- `sum()`
    compensates a float term and adds an int one plainly -- and it is worth the
    last bit of the answer on ten of the twenty-eight gravy scales, four of which
    can actually answer differently for it.

    A key that is not a single character, or is one above U+00FF, cannot be
    addressed by a byte, and a scale containing one is not one a dense table can
    be built for.  Neither can a value that is neither an `int` nor a `float`,
    because which of the two CPython's `sum()` does with it is not a thing this
    kernel implements.  No scale in :mod:`biofasting.protparam_data` has either,
    so both answers are a fallback rather than a wrong profile.
    """
    values = bytearray(256 * 8)
    valid = bytearray(256)
    sum_flags = bytearray(256)
    for letter, value in scale.items():
        if len(letter) != 1 or ord(letter) > 0xFF:
            return None
        kind = type(value)
        if kind is float:
            flag = 1
        elif kind is int:
            flag = 2
        else:
            return None
        code = ord(letter)
        struct.pack_into("<d", values, code * 8, value)
        valid[code] = 1
        sum_flags[code] = flag
    return bytes(values), bytes(valid), bytes(sum_flags)


# Keyed on `id`, which is only sound because the pair is kept: an entry holds
# the scale it was built from, so the dict cannot be collected and its id cannot
# be handed to something else.  The lookup re-checks identity anyway, because
# `SCALES_BY_NAME` is public and a caller who empties it must not be served a
# stale table.
_SCALE_BLOBS = {
    id(scale): (scale, blob)
    for scale, blob in (
        (scale, _dense_scale(scale)) for scale in SCALES_BY_NAME.values()
    )
    if blob is not None
}


def _scale_blob(scale):
    """The dense arrays for one of the scales this package ships, or `None`.

    Identity, not equality: `ProtParamData.kd` is the fast path, a copy of it is
    a perfectly good scale that takes the reference's own loop.  That is the
    same speed Biopython runs at, so nothing is lost by the fallback -- but it is
    worth knowing which one is in hand, because the difference is a factor of
    twenty on a long protein.
    """
    entry = _SCALE_BLOBS.get(id(scale))
    if entry is None or entry[0] is not scale:
        return None
    return entry[1]


def _latin1(text):
    """`text` as bytes, or `None` if it holds a character above Latin-1."""
    try:
        return text.encode("latin-1")
    except UnicodeEncodeError:
        return None


#: The position of each residue in `DIWV`, or 0xFF where there is none.  Shared
#: by the two kernels that index a twenty-letter table.
_AA_MAP = bytes(
    bytearray(
        {aa: i for i, aa in enumerate(AMINO_ACIDS)}.get(chr(byte), 0xFF)
        for byte in range(256)
    )
)

#: `Flex`, twenty doubles.
_FLEX_TABLE = struct.pack("<20d", *FLEX)

#: `DIWV`, four hundred doubles row-major in `AMINO_ACIDS`.
_DIWV_TABLE = struct.pack("<400d", *(value for row in _DIWV_ROWS for value in row))


class ProteinAnalysis:
    """The classical protein parameters, over one sequence.

    >>> protein = ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT")
    >>> protein.count_amino_acids()["A"]
    2
    >>> round(protein.gravy(), 2)
    -0.31
    >>> round(protein.aromaticity(), 3)
    0.102
    >>> round(protein.instability_index(), 2)
    25.07

    `monoisotopic` selects the monoisotopic rather than the average mass, and
    reaches only `molecular_weight`.

    The sequence is upper-cased once, at construction, which is the reference's
    own choice and the reason `count_amino_acids` can ask for `"A"` without also
    asking for `"a"`.
    """

    def __init__(self, prot_sequence, monoisotopic=False):
        self.sequence = _as_text(prot_sequence).upper()
        self.amino_acids_content = None
        self.length = len(self.sequence)
        self.monoisotopic = monoisotopic

    def count_amino_acids(self):
        """Count each standard residue, as ``{letter: count}``.

        Twenty `str.count` calls, the first time, in the reference.  Here it is
        one pass over the bytes with a table lookup per residue, which is the
        same twenty numbers: the kernel counts into a bucket and never declines,
        because a residue outside the twenty is a zero and not an error.

        The answer is cached in `self.amino_acids_content` and never recomputed,
        which is the reference's own behaviour and also its own trap: a timed call
        is usually a cache hit, and a cache hit says nothing about the twenty
        scans behind it.  Which is the reason this kernel exists at all -- six
        other methods begin by calling this one.
        """
        if self.amino_acids_content is None:
            data = _latin1(self.sequence)
            if data is None:
                #   A character above Latin-1 cannot be counted by byte, and the
                #   reference counts it as nothing.  Its own loop is used rather
                #   than a table of surrogates, because the set of such
                #   characters is not enumerable and the fallback has to be
                #   correct for all of it.
                counts = {aa: self.sequence.count(aa) for aa in AMINO_ACIDS}
            else:
                counts = dict(zip(AMINO_ACIDS,
                                  _core.count_residues(data, _AA_MAP)))

            self.amino_acids_content = counts

        return self.amino_acids_content

    @functools.cached_property
    def amino_acids_percent(self):
        """The same counts as percentages of the sequence length.

        A `cached_property` and not a method -- since 1.88, and calling it is
        the bug it looks like.  An empty sequence raises `ZeroDivisionError`
        here, once, from the first residue's percentage.
        """
        aa_counts = self.count_amino_acids()
        percentages = {
            aa: (count * 100 / self.length) for aa, count in aa_counts.items()
        }

        return percentages

    def molecular_weight(self):
        """The average or monoisotopic mass, via :func:`~biofasting.sequtils.molecular_weight`."""
        return molecular_weight(
            self.sequence, seq_type="protein", monoisotopic=self.monoisotopic
        )

    def aromaticity(self):
        """Lobry's aromaticity: the relative frequency of Phe, Trp and Tyr."""
        aromatic_aas = "YWF"
        aa_percentages = self.amino_acids_percent

        aromaticity = sum(aa_percentages[aa] / 100 for aa in aromatic_aas)

        return aromaticity

    def instability_index(self):
        """Guruprasad's instability index.

        Above 40 is taken to mean an unstable protein.  The kernel adds the
        `size - 1` dipeptide weights in sequence order with a plain running
        total, and the multiply by `10.0 / length` is done here, in that order,
        because for the empty sequence it is the division that raises.

        A residue `DIWV` has no row, or no column, for takes the reference's own
        loop, which raises `KeyError` naming it -- `this` before `following`,
        because that is the order the reference indexes them in.
        """
        text = self.sequence
        data = _latin1(text)
        if data is not None:
            total = _core.instability_index_sum(data, _DIWV_TABLE, _AA_MAP)
            if total is not None:
                return (10.0 / self.length) * total

        index = DIWV
        score = 0.0

        for i in range(self.length - 1):
            this, following = text[i : i + 2]
            dipeptide_value = index[this][following]
            score += dipeptide_value

        return (10.0 / self.length) * score

    def flexibility(self):
        """Vihinen's flexibility, a window of nine, no window argument.

        The kernel reproduces four things that look like mistakes and are the
        reference: the pair is summed before the weight is applied, the middle
        residue is read at index 5 so that index 4 is never read and index 5 is
        read twice, there is one window fewer than fit into the sequence, and
        the divisor is the literal 5.25.
        """
        text = self.sequence
        data = _latin1(text)
        if data is not None:
            scores = _core.flexibility_scores(data, _FLEX_TABLE, _AA_MAP)
            if scores is not None:
                return scores

        flexibilities = Flex
        window_size = 9
        weights = [0.25, 0.4375, 0.625, 0.8125, 1]
        scores = []

        for i in range(self.length - window_size):
            subsequence = text[i : i + window_size]
            score = 0.0

            for j in range(window_size // 2):
                front = subsequence[j]
                back = subsequence[window_size - j - 1]
                score += (flexibilities[front] + flexibilities[back]) * weights[j]

            middle = subsequence[window_size // 2 + 1]
            score += flexibilities[middle]

            scores.append(score / 5.25)

        return scores

    def gravy(self, scale="KyteDoolitle"):
        """The grand average of hydropathy under one of the twenty-eight scales.

        ``sum(selected_scale[aa] for aa in self.sequence)`` over the length, and
        that `sum` is load-bearing in two ways: CPython compensates its float
        terms, and it does *not* compensate the int ones -- which nine of the
        twenty-eight scales have, and which the kernel is told about through
        `sum_flags`.

        The empty sequence is the one input that does not go through the kernel.
        There the reference's sum is the *integer* zero rather than a double, and
        its ``0 / 0`` says "division by zero" where a kernel's total would say
        "float division by zero".  Same exception, different message, and the
        message is part of the answer.
        """
        selected_scale = gravy_scales.get(scale, -1)

        if selected_scale == -1:
            raise ValueError(f"scale: {scale} not known")

        text = self.sequence
        blob = _scale_blob(selected_scale)
        if blob is not None and self.length:
            data = _latin1(text)
            if data is not None:
                total, bad = _core.molecular_weight_mass(data, blob[0], blob[2])
                if bad < 0:
                    return total / self.length

        # The reference's own expression, kept whole.  `bad` is the byte the
        # scale has no entry for, and the reference names it in a `KeyError`
        # raised from inside this generator -- which only Python can spell.
        total_gravy = sum(selected_scale[aa] for aa in text)

        return total_gravy / self.length

    def _weight_list(self, window, edge):
        """Half a window of linear weights, from `edge` at the rim to 1 at the centre.

        For a window of 9 and an edge of 0.4, ``[0.4, 0.55, 0.7, 0.85]``.  The
        other half is the mirror of this one and is never built, which is why
        `protein_scale` doubles the sum and adds the centre's 1 back.

        ``window == 1`` divides by zero here.  That is the reference's, and it
        arrives before `protein_scale` reaches its loop.
        """
        unit = 2 * (1.0 - edge) / (window - 1)
        weights = [0.0] * (window // 2)

        for i in range(window // 2):
            weights[i] = edge + unit * i

        return weights

    def protein_scale(self, param_dict, window, edge=1.0):
        """A sliding profile under any amino-acid scale.

        One score per window that fits, ``size - window + 1`` of them, each
        divided by the sum of the window's weights.  `edge` is the weight at the
        rim: at 0.4 and a window of 5 the weights are 0.4, 0.7, 1.0, 0.7, 0.4.

        The weights are computed here and handed to the kernel as they came out
        of `_weight_list`, rather than recomputed there, so that there is one
        place and not two where `edge + unit * i` rounds.

        A residue the scale has no value for does not stop the reference -- it
        writes a warning naming the residue and drops that window's term -- so
        the kernel declines on any such byte and the reference's loop runs
        below, warning and all.  `window < 2` is declined for the reason
        `_weight_list` gives.
        """
        weights = self._weight_list(window, edge)
        sum_of_weights = sum(weights) * 2 + 1

        text = self.sequence
        blob = _scale_blob(param_dict)
        if blob is not None and window >= 2:
            data = _latin1(text)
            if data is not None:
                scores = _core.protein_scale_scores(
                    data,
                    window,
                    struct.pack(f"<{len(weights)}d", *weights),
                    blob[0],
                    blob[1],
                    sum_of_weights,
                )
                if scores is not None:
                    return scores

        scores = []

        for i in range(self.length - window + 1):
            subsequence = text[i : i + window]
            score = 0.0

            for j in range(window // 2):
                # Walk from the outside of the window towards the middle.  The
                # try/except is the reference's, and it is not an error path: a
                # non-standard residue costs its window one term and the rest of
                # the window is still scored.
                try:
                    front = param_dict[subsequence[j]]
                    back = param_dict[subsequence[window - j - 1]]
                    score += weights[j] * front + weights[j] * back
                except KeyError:
                    sys.stderr.write(
                        "warning: %s or %s is not a standard "
                        "amino acid.\n" % (subsequence[j], subsequence[window - j - 1])
                    )

            # The middle of the window always has a weight of 1.
            middle = subsequence[window // 2]
            if middle in param_dict:
                score += param_dict[middle]
            else:
                sys.stderr.write(f"warning: {middle} is not a standard amino acid.\n")

            scores.append(score / sum_of_weights)

        return scores

    def isoelectric_point(self):
        """The pH at which the net charge is zero, via `IsoelectricPoint`."""
        aa_content = self.count_amino_acids()

        ie_point = IsoelectricPoint(self.sequence, aa_content)
        return ie_point.pi()

    def charge_at_pH(self, pH):
        """The net charge at one pH, via `IsoelectricPoint`."""
        aa_content = self.count_amino_acids()
        charge = IsoelectricPoint(self.sequence, aa_content)
        return charge.charge_at_pH(pH)

    def secondary_structure_fraction(self):
        """The fractions of residues that prefer helix, turn and sheet.

        Helix: E, M, A, L, K.  Turn: N, P, G, S, D.  Sheet: V, I, Y, F, W, L, T.
        A tuple of three floats, in that order -- which before 1.82 was the
        reverse of what the docstring claimed, and is not now.
        """
        aa_percentages = self.amino_acids_percent

        helix = sum(aa_percentages[r] / 100 for r in "EMALK")
        turn = sum(aa_percentages[r] / 100 for r in "NPGSD")
        sheet = sum(aa_percentages[r] / 100 for r in "VIYFWLT")

        return helix, turn, sheet

    def molar_extinction_coefficient(self):
        """The molar extinction coefficient at 280 nm, reduced and oxidised.

        Trp counts 5500 and Tyr 1490; a cystine, which is two cysteines, adds
        125 -- so the oxidised figure counts *pairs* of cysteines and an odd one
        is discarded.  Both numbers are integers, and both are exact.
        """
        num_aa = self.count_amino_acids()
        mec_reduced = num_aa["W"] * 5500 + num_aa["Y"] * 1490
        mec_cystines = mec_reduced + (num_aa["C"] // 2) * 125
        return (mec_reduced, mec_cystines)
