"""Sequence operations: reverse complement, GC fraction, k-mer counts, translation.

The kernels live in the compiled core and take a byte buffer, which is the
right contract for something that has to run over a 100 Mbp chromosome without
turning it into a million Python objects.  This module is the thin layer that
makes them convenient from Python: it accepts `str` as well as bytes-like
input and gives back the type it was handed.

The semantics are Biopython's, and the tests are differential, so the
functions here are drop-in for `Bio.Seq.reverse_complement`,
`Bio.SeqUtils.gc_fraction`, per-k-mer `Seq.count_overlap`, and -- with the
genetic codes, which are the one thing here that needs data rather than a
kernel -- `Bio.Seq.translate`.  The last of those lives in `_translate.py` and
is re-exported here; the genetic code a translation needs is imported inside
the call that needs it, so that `import biofasting` does not build 162 tables
for a program that will never translate anything.

One boundary worth stating rather than discovering: the core works on bytes, so
a `str` is encoded as Latin-1 on the way in and decoded on the way out.  Every
character up to U+00FF therefore round-trips exactly, which covers every
sequence alphabet in existence, and a character above U+00FF raises
:class:`UnicodeEncodeError` instead of being quietly mangled -- it cannot be
represented at all in the domain the operation is defined on.
"""

from . import _core
from ._iupac_data import GC_VALUES as _gc_values
from ._translate import BiopythonWarning, translate

__all__ = [
    "BiopythonWarning",
    "count_kmers",
    "gc_fraction",
    "reverse_complement",
    "translate",
]


def _as_buffer(sequence):
    """Return ``(buffer, was_text)`` for a ``str`` or any bytes-like object."""
    if isinstance(sequence, str):
        return sequence.encode("latin-1"), True
    return sequence, False


def reverse_complement(sequence):
    """Return the reverse complement of `sequence`.

    Byte for byte `Bio.Seq.reverse_complement`: the full IUPAC set in both
    cases, and every byte that is not a nucleotide left as it was, so ``'-'``
    stays ``'-'`` and ``'N'`` stays ``'N'``.  Accepts `str` or any object
    exposing the buffer protocol; returns `str` for `str` and `bytes`
    otherwise.
    """
    buffer, was_text = _as_buffer(sequence)
    result = _core.reverse_complement(buffer)
    return result.decode("latin-1") if was_text else result


def gc_fraction(sequence, ambiguous="remove"):
    """Return the GC fraction of `sequence`, as Biopython computes it.

    `Bio.SeqUtils.gc_fraction`, all three of its ``ambiguous`` modes included.
    With ``"remove"`` (the default) G, C and S count as GC, A, T, W and U make
    up the rest of the denominator, and anything else -- N included -- is
    removed from both; with ``"ignore"`` the ambiguous letters are still left
    out of the numerator but the denominator is the whole sequence; with
    ``"weighted"`` each ambiguous letter counts for the fraction of GC it could
    stand for.  Any other value raises :class:`ValueError`.

    A sequence with nothing countable gives the reference's answer, which is the
    *integer* ``0`` and not ``0.0``.  The default mode is answered by the
    vectorised counting kernel and the other two by a table-driven one, which is
    two kernels for one function on purpose: the first is the measured fast path
    and moving the other two onto it would have made every number in
    `bench/targets.md` about `gc_fraction` incomparable for no gain.

    >>> gc_fraction("ACTG")
    0.5
    >>> gc_fraction("ACTGN", "ignore")
    0.4
    >>> gc_fraction("ACTGN", "weighted")
    0.5
    """
    if ambiguous not in ("weighted", "remove", "ignore"):
        raise ValueError(f"ambiguous value '{ambiguous}' not recognized")

    buffer, _ = _as_buffer(sequence)
    if ambiguous == "remove":
        fraction = _core.gc_fraction(buffer)
        return 0 if fraction < 0 else fraction

    gc, at, *ambiguous_counts = _core.gc_counts(buffer)
    if ambiguous == "weighted":
        gc += sum(
            count * _gc_values[x] for count, x in zip(ambiguous_counts, "BDHKMNRVXY")
        )

    # `len(sequence)` and not the kernel's byte count, because that is what the
    # reference measures: for a `str` these two are the same only while every
    # character is ASCII, and the reference is defined on characters.
    length = len(sequence)
    if length == 0:
        return 0
    return gc / length


def count_kmers(sequence, k):
    """Return overlapping k-mer counts as ``{kmer: count}``, sorted by k-mer.

    Case-insensitive, and counted over ACGT only: a window containing anything
    else -- N included -- is skipped rather than split around, so the counts are
    exactly what ``Seq.count_overlap`` gives per key.  Anything not found is
    simply absent, which is why a caller comparing against ``count_overlap``
    should read a missing key as zero.

    Keys are `str` even when the input was bytes -- unlike
    :func:`reverse_complement`, which gives back the type it was handed.  A
    k-mer over ACGT is ASCII by construction, and a table is only useful if it
    can be looked up with the string the caller is holding.

    `k` must be between 1 and 16, and anything outside that raises
    :class:`ValueError`.  The bound is the length of the code the core slides
    over the sequence (two bits per base in a machine word), and it is not a
    practical restriction: counting a larger k means a table of 4**k entries,
    which no caller wants to materialise anyway.
    """
    buffer, _ = _as_buffer(sequence)
    return {
        kmer.decode("ascii"): count for kmer, count in _core.count_kmers(buffer, k)
    }
