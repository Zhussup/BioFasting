"""`Bio.SeqUtils`: the functions that measure a sequence rather than change it.

`biofasting.seqops` holds the operations -- reverse complement, translation,
the two countings that were already there.  This module holds the other half:
composition, molecular mass, the protein alphabets, the six-frame picture, and
the codon adaptation index.  The two are the same module in Biopython and are
split here because the contracts are different in kind.  An operation is judged
on the bytes it emits; a measurement is judged on a number, and the arithmetic
that turns counts into that number is worth being able to read.

That split is why the kernels underneath return **counts** and not answers.
`GC123` divides three times and raises on the fourth; a frame with no A, C, G or
T in it answers the *integer* ``0`` while a sequence with no counted base at all
raises :class:`ZeroDivisionError`.  None of those distinctions survives a round
trip through a C++ `double`, so the kernel counts and this module divides.

The input convention follows from the same idea, and it is not the one
:mod:`biofasting.seqops` uses.  These functions are defined on **text**: the
reference formats its argument with ``str()``, upper-cases it, measures it with
``len()`` and compares single characters.  So a `str` is read as itself, and a
byte buffer is read as Latin-1, one byte per character, which is how
`Bio.Seq.Seq` reads one too.  Where a character above U+00FF would change the
answer -- `molecular_weight` names it in a :class:`ValueError`, `GC123` counts
the positions after it differently, `gcg` upper-cases it through Python's
Unicode tables -- the kernel declines and the reference's own Python runs
instead.  Nothing is approximated and nothing is silently dropped; those inputs
take the reference's slow path and get the reference's answer.
"""

import re
from array import array

from . import _core
from ._iupac_data import (
    AMBIGUOUS_DNA_VALUES,
    DNA_COMPLEMENT,
    DNA_WEIGHTS,
    DNA_WEIGHTS_MONOISOTOPIC,
    GC_VALUES,
    PROTEIN_1TO3,
    PROTEIN_3TO1,
    PROTEIN_WEIGHTS,
    PROTEIN_WEIGHTS_MONOISOTOPIC,
    RNA_COMPLEMENT,
    RNA_WEIGHTS,
    RNA_WEIGHTS_MONOISOTOPIC,
)
from .seqops import gc_fraction, reverse_complement, translate

__all__ = [
    "CodonAdaptationIndex",
    "GC123",
    "GC_skew",
    "gc_fraction",
    "molecular_weight",
    "nt_search",
    "seq1",
    "seq3",
    "six_frame_translations",
    "xGC_skew",
]


def _as_text(sequence):
    """The characters the reference would look at, for `str`, `Seq` or bytes.

    A `str` is already text.  Anything else -- `Bio.Seq.Seq`, `bytearray`,
    `memoryview`, anything with ``__bytes__`` -- is decoded as Latin-1, which
    cannot fail and maps every byte to the character with that code point, which
    is exactly what `Seq` does for the bytes it wraps.
    """
    if isinstance(sequence, str):
        return sequence
    try:
        return bytes(sequence).decode("latin-1")
    except TypeError:
        return str(sequence)


# ---------------------------------------------------------------------------
# Composition.


def GC123(seq):
    """G+C content of the whole sequence and of each of the three codon positions.

    Returns four floats as percentages, in the order
    ``(overall, position 1, position 2, position 3)``.

    >>> GC123("ACTGTN")
    (40.0, 50.0, 50.0, 0.0)

    Only A, C, G and T are looked at, in either case: an ambiguous base is
    counted by nobody, and so is the space the reference pads a partial trailing
    codon with -- which is why a 40 Mbp chromosome, whose length is not a
    multiple of three, reports a slightly different denominator than its own
    length.  A sequence with nothing countable raises :class:`ZeroDivisionError`,
    exactly as the reference does.
    """
    text = _as_text(seq)
    if text.isascii():
        counts = _core.gc123(text.encode("latin-1"))
    else:
        counts = _reference_gc123_counts(text)

    # Twelve counts, frame-major, each frame's four bases in the order A, T, G,
    # C -- and then the reference's own arithmetic, division by division, so
    # that the one frame that can be empty answers the integer 0 and the whole
    # sequence that is empty raises.
    gcall = 0
    nall = 0
    gc = []
    for frame in range(3):
        a, t, g, c = counts[4 * frame : 4 * frame + 4]
        n = a + t + g + c
        gc.append((g + c) * 100.0 / n if n else 0)
        gcall += g + c
        nall += n
    gcall = 100.0 * gcall / nall
    return gcall, gc[0], gc[1], gc[2]


def _reference_gc123_counts(text):
    """`GC123`'s loop, kept only for text a byte kernel cannot count.

    Not dead code and not a duplicate: a non-ASCII character is one character to
    the reference and several bytes here, so it moves every codon boundary after
    it.  Rather than encode the Unicode rules for what counts as a letter, this
    path asks the reference's own question in the reference's own words.
    """
    counts = [0] * 12
    for i in range(0, len(text), 3):
        codon = text[i : i + 3]
        if len(codon) < 3:
            codon += "  "
        for pos in range(3):
            for index, nt in enumerate("ATGC"):
                if codon[pos] == nt or codon[pos] == nt.lower():
                    counts[4 * pos + index] += 1
    return counts


def GC_skew(seq, window=100):
    """GC skew ``(G - C) / (G + C)`` over non-overlapping windows.

    Returns a list of floats, one per window; the last window is short if the
    length is not a multiple of `window`.  A window with no G and no C at all is
    ``0.0``, which is the reference's ``except ZeroDivisionError`` arm and not a
    value that could be confused with a real skew of zero -- the reference has
    the same ambiguity.

    >>> GC_skew("GGGCCC", 3)
    [1.0, -1.0]
    """
    text = _as_text(seq)
    if not isinstance(window, int) or window <= 0:
        # The reference's first act is `range(0, len(seq), window)`, and on a
        # window it cannot use, that call is what raises -- ValueError for zero,
        # TypeError for anything that is not an integer.  Making the call is how
        # this keeps the reference's own message instead of inventing one.  A
        # negative window leaves the loop with nothing to do and answers [].
        list(range(0, len(text), window))
        return []
    if text.isascii():
        return _core.gc_skew(text.encode("latin-1"), window)
    values = []
    for i in range(0, len(text), window):
        s = text[i : i + window]
        g = s.count("G") + s.count("g")
        c = s.count("C") + s.count("c")
        values.append((g - c) / (g + c) if g + c else 0.0)
    return values


# ---------------------------------------------------------------------------
# Molecular mass.

#: The six weight tables the reference selects between, keyed the way
#: `molecular_weight` selects them.
_WEIGHT_TABLES = {
    ("DNA", False): DNA_WEIGHTS,
    ("DNA", True): DNA_WEIGHTS_MONOISOTOPIC,
    ("RNA", False): RNA_WEIGHTS,
    ("RNA", True): RNA_WEIGHTS_MONOISOTOPIC,
    ("protein", False): PROTEIN_WEIGHTS,
    ("protein", True): PROTEIN_WEIGHTS_MONOISOTOPIC,
}

#: The same tables as the kernel wants them: 256 native-endian doubles and 256
#: flags, built once.  A `dict` lookup per letter would have been the whole cost
#: of a 150 bp sequence, and building the dense form per call would have cost
#: more than the sum it feeds.
_DENSE = {}


def _dense_weights(seq_type, monoisotopic):
    key = (seq_type, monoisotopic)
    cached = _DENSE.get(key)
    if cached is None:
        table = _WEIGHT_TABLES[key]
        values = array("d", bytes(8 * 256))
        valid = bytearray(256)
        for letter, mass in table.items():
            values[ord(letter)] = mass
            valid[ord(letter)] = 1
        cached = (values.tobytes(), bytes(valid))
        _DENSE[key] = cached
    return cached


def molecular_weight(
    seq, seq_type="DNA", double_stranded=False, circular=False, monoisotopic=False
):
    """Molecular mass of a DNA, RNA or protein sequence, as a float.

    Only unambiguous letters are allowed, and the sequence is stripped of
    whitespace and upper-cased first, exactly as the reference does it.  A
    letter the table does not have raises :class:`ValueError` naming it.

    >>> round(molecular_weight("AGC"), 2)
    949.61
    >>> round(molecular_weight("AGC", "protein"), 2)
    249.29
    """
    try:
        seq = seq.seq
    except AttributeError:  # not a SeqRecord
        pass
    # `str(seq)`, and not a friendlier coercion, because the reference's own
    # behaviour for a byte buffer is a consequence of it: `str(b"ACGT")` is the
    # repr, and the first letter of `"B'ACGT'"` is a B, which is not a DNA
    # letter.  The reference raises ValueError there, and so does this.
    text = "".join(str(seq).split()).upper()

    if seq_type not in ("DNA", "RNA", "protein"):
        raise ValueError(f"Allowed seq_types are DNA, RNA or protein, not {seq_type!r}")

    water = 18.010565 if monoisotopic else 18.0153

    if text.isascii():
        weights, valid = _dense_weights(seq_type, monoisotopic)
        mass, bad = _core.molecular_weight_mass(
            text.encode("latin-1"), weights, valid
        )
        if bad >= 0:
            _raise_bad_letter(chr(bad), seq_type)
    else:
        # The reference's expression, kept whole for the inputs a byte kernel
        # has no answer for: a character above ASCII may upper-case into a
        # different letter (and, for some, into two), and that letter is the one
        # the error message has to name.
        table = _WEIGHT_TABLES[(seq_type, monoisotopic)]
        try:
            mass = sum(table[x] for x in text)
        except KeyError as e:
            raise ValueError(
                f"'{e}' is not a valid unambiguous letter for {seq_type}"
            ) from None

    weight = mass - (len(text) - 1) * water
    if circular:
        weight -= water

    if double_stranded:
        if seq_type == "protein":
            raise ValueError("protein sequences cannot be double-stranded")
        text = text.translate(DNA_COMPLEMENT if seq_type == "DNA" else RNA_COMPLEMENT)
        if text.isascii():
            weights, valid = _dense_weights(seq_type, monoisotopic)
            mass, bad = _core.molecular_weight_mass(
                text.encode("latin-1"), weights, valid
            )
            if bad >= 0:
                _raise_bad_letter(chr(bad), seq_type)
        else:
            table = _WEIGHT_TABLES[(seq_type, monoisotopic)]
            try:
                mass = sum(table[x] for x in text)
            except KeyError as e:
                raise ValueError(
                    f"'{e}' is not a valid unambiguous letter for {seq_type}"
                ) from None
        weight += mass - (len(text) - 1) * water
        if circular:
            weight -= water

    return weight


def _raise_bad_letter(letter, seq_type):
    """The reference's own KeyError message, quotes and all.

    ``f"'{e}'"`` where `e` is a :class:`KeyError` is two quotes around a repr,
    so ``''N'' is not a valid unambiguous letter for DNA`` is what Biopython
    prints and what this prints.  It looks like a typo and it is not one.
    """
    raise ValueError(
        f"'{KeyError(letter)}' is not a valid unambiguous letter for {seq_type}"
    ) from None


# ---------------------------------------------------------------------------
# Protein alphabets.


def seq3(seq, custom_map=None, undef_code="Xaa"):
    """One-letter protein sequence to three-letter codes.

    >>> seq3("MAIVMGRWKGAR*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer'
    """
    if custom_map is None:
        custom_map = {"*": "Ter"}
    threecode = dict(PROTEIN_1TO3)
    threecode.update(custom_map)
    return "".join(threecode.get(aa, undef_code) for aa in _as_text(seq))


def seq1(seq, custom_map=None, undef_code="X"):
    """Three-letter protein sequence to one-letter codes, case-insensitively.

    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer")
    'MAIVMGRWKGAR*'
    """
    if custom_map is None:
        custom_map = {"Ter": "*"}
    onecode = {k.upper(): v for k, v in PROTEIN_3TO1.items()}
    onecode.update((k.upper(), v) for k, v in custom_map.items())
    text = _as_text(seq)
    # Three characters at a time, and any remainder dropped -- the reference
    # writes `range(len(seq) // 3)`, so a trailing "Ar" is not a codon and is
    # not reported as an unknown one either.
    seqlist = [text[3 * i : 3 * (i + 1)] for i in range(len(text) // 3)]
    return "".join(onecode.get(aa.upper(), undef_code) for aa in seqlist)


# ---------------------------------------------------------------------------
# Degenerate search.


def nt_search(seq, subseq):
    """Positions of an ambiguous DNA `subseq` in `seq`, forward strand only.

    Returns the regular expression it built, followed by the 0-based start of
    every match.

    >>> nt_search("ACGTACGT", "ACGT")
    ['ACGT', 0, 4]

    Both arguments are used exactly as the reference uses them and are not
    coerced: the reference re-slices ``seq[pos:]`` on every iteration, which is
    quadratic in the number of matches, and it requires a `str` for `subseq`
    because it looks each letter up in a dictionary keyed by letter.  Passing
    anything else fails here in the same place and with the same message as it
    fails there.
    """
    pattern = ""
    for nt in subseq:
        value = AMBIGUOUS_DNA_VALUES[nt]
        if len(value) == 1:
            pattern += value
        else:
            pattern += f"[{value}]"

    pos = -1
    result = [pattern]
    while True:
        pos += 1
        s = seq[pos:]
        m = re.search(pattern, s)
        if not m:
            break
        pos += int(m.start(0))
        result.append(pos)
    return result


# ---------------------------------------------------------------------------
# The six-frame picture.


def six_frame_translations(seq, genetic_code=1):
    """The six reading frames of `seq` as the pretty block text Biopython prints.

    >>> print(six_frame_translations("AUGGCCAUUGUAAUGGGCCGCUGA").splitlines()[2])
    <BLANKLINE>
    """
    text = _as_text(seq)
    if "u" in text.lower():
        anti = text.translate(RNA_COMPLEMENT)[::-1]
    else:
        anti = reverse_complement(text)
    comp = anti[::-1]
    length = len(text)
    frames = {}
    for i in range(3):
        fragment_length = 3 * ((length - i) // 3)
        frames[i + 1] = translate(text[i : i + fragment_length], genetic_code)
        frames[-(i + 1)] = translate(anti[i : i + fragment_length], genetic_code)[::-1]

    if length > 20:
        short = f"{text[:10]} ... {text[-10:]}"
    else:
        short = text
    header = "GC_Frame:"
    for nt in ["a", "t", "g", "c"]:
        header += " %s:%d" % (nt, text.count(nt.upper()))

    gc = 100 * gc_fraction(text, ambiguous="ignore")
    header += "\nSequence: %s, %d nt, %0.2f %%GC\n\n\n" % (
        short.lower(),
        length,
        gc,
    )
    res = header

    for i in range(0, length, 60):
        subseq = text[i : i + 60]
        csubseq = comp[i : i + 60]
        p = i // 3
        res += "%d/%d\n" % (i + 1, i / 3 + 1)
        res += "  " + "  ".join(frames[3][p : p + 20]) + "\n"
        res += " " + "  ".join(frames[2][p : p + 20]) + "\n"
        res += "  ".join(frames[1][p : p + 20]) + "\n"
        res += subseq.lower() + "%5d %%\n" % int(gc)
        res += csubseq.lower() + "\n"
        res += "  ".join(frames[-2][p : p + 20]) + "\n"
        res += " " + "  ".join(frames[-1][p : p + 20]) + "\n"
        res += "  " + "  ".join(frames[-3][p : p + 20]) + "\n\n"
    return res


# ---------------------------------------------------------------------------
# Codon adaptation.


class CodonAdaptationIndex(dict):
    """A codon adaptation index, after Sharp and Li (1987).

    A `dict` from codon to relative adaptiveness, `w_ij`, built from a set of
    coding sequences.  The reference subclasses `dict` and so does this, so
    ``index["CTG"]`` is a number and ``dict(index)`` is the whole table.

    One difference, and it is deliberate: the reference's default argument is
    the table object itself, which makes `Bio.Data.CodonTable` load at import
    time.  Here the default is `None` and the standard table is imported when a
    `CodonAdaptationIndex` is actually built, so that `import biofasting` does
    not pay for 162 genetic codes a program may never translate anything with.
    Every call that works there works here.
    """

    def __init__(self, sequences, table=None):
        if table is None:
            from .codon_table import standard_dna_table

            table = standard_dna_table
        self._table = table
        codons = {aminoacid: [] for aminoacid in table.protein_alphabet}
        for codon, aminoacid in table.forward_table.items():
            codons[aminoacid].append(codon)
        synonymous_codons = tuple(list(codons.values()) + [table.stop_codons])

        # Every codon, in the reference's order (A, C, G, T at each position),
        # so that the class's `dict` iterates the way Biopython's does and
        # `__str__` prints the same table in the same order.
        counts = {c1 + c2 + c3: 0 for c1 in "ACGT" for c2 in "ACGT" for c3 in "ACGT"}
        self.update(counts)

        for sequence in sequences:
            try:  # SeqRecord
                name = sequence.id
                sequence = sequence.seq
            except AttributeError:  # str, Seq, or MutableSeq
                name = None
            sequence = sequence.upper()
            for i in range(0, len(sequence), 3):
                codon = sequence[i : i + 3]
                try:
                    counts[codon] += 1
                except KeyError:
                    if name is None:
                        message = f"illegal codon '{codon}'"
                    else:
                        message = f"illegal codon '{codon}' in gene {name}"
                    raise ValueError(message) from None

        for codon, count in counts.items():
            if count == 0:
                counts[codon] = 0.5

        for codons in synonymous_codons:
            denominator = max(counts[codon] for codon in codons)
            for codon in codons:
                self[codon] = counts[codon] / denominator

    def calculate(self, sequence):
        """The CAI of `sequence`, as a float.

        ATG and TGG are skipped because their index is always one.  A stop codon
        that the index does not have is skipped too; any other codon it does not
        have raises :class:`TypeError`.
        """
        from math import exp, log

        cai_value, cai_length = 0, 0

        try:
            sequence = sequence.seq  # SeqRecord
        except AttributeError:
            pass  # str, Seq, or MutableSeq
        sequence = sequence.upper()

        for i in range(0, len(sequence), 3):
            codon = sequence[i : i + 3]
            if codon in ["ATG", "TGG"]:
                continue
            try:
                cai_value += log(self[codon])
            except KeyError:
                if codon in ["TGA", "TAA", "TAG"]:
                    continue
                raise TypeError(f"illegal codon in sequence: {codon}") from None
            else:
                cai_length += 1

        return exp(cai_value / cai_length)

    def __str__(self):
        lines = [f"{codon}\t{value:.3f}" for codon, value in self.items()]
        return "\n".join(lines) + "\n"


def xGC_skew(seq, window=1000, zoom=100, r=300, px=100, py=100):
    """Not implemented.  ``Bio.SeqUtils.xGC_skew`` draws a window, it does not measure.

    The reference imports `tkinter`, creates a `Canvas`, draws the normal and
    accumulated skew curves onto it, calls ``canvas.update()`` and returns
    `None`.  Measured at 237 ms per call with a display attached, it has neither
    a result to compare against nor an algorithm to move into a kernel.  There is
    nothing here for a port to be faster or slower than; what a caller wants is
    :func:`GC_skew`, which returns the numbers the picture is drawn from.
    """
    raise NotImplementedError(
        "Bio.SeqUtils.xGC_skew draws a Tk window; use GC_skew for the values"
    )
