"""Translation: the reference's rules, and a kernel for the reads that fit them.

``Bio.Seq``'s ``translate`` is one function with a surprising number of ways to
answer, and every one of them is the contract:

* the table may be an id, a name, or a table object, and an id or a name
  resolves to the *ambiguous generic* table rather than to the DNA one, which is
  why ``"TAR"`` translates to a stop and ``"TAN"`` to ``"X"``;
* a codon that is not in the table but every letter of which is a nucleotide is
  ``pos_stop`` (``"X"``), not an error -- while a codon containing a byte that
  is not a nucleotide for that table is a ``TranslationError``;
* a table with a codon that is both a stop and an amino acid warns, and refuses
  ``to_stop`` outright, because "stop here" and "this is also a residue" have no
  common answer;
* ``stop_symbol`` is never checked, so it may be any string at all, ``""``
  included;
* ``gap`` is checked for type and length even when no gap codon occurs;
* ``cds`` checks the first codon, the length and the last codon, then translates
  the first as ``M`` whatever it was, and treats an in-frame stop *inside* as an
  error rather than as a symbol.

The compiled kernel does none of that list.  It is a single scan over a
translation table -- two bits per base, twelve bits per codon, one byte out --
and it does the *whole* job or none of it: the moment it meets a codon its two
tables cannot classify, it declines, and every rule above is applied by the
Python below instead.  That is the same shape as the FASTQ fast path and the
AVX2 reverse-complement kernel, and it is what keeps "the fast path agrees with
the reference" a property that can be checked by exhausting a table rather than
argued from a character class.

So there are two implementations, and they have to agree.  The tests ask the
same question of both on every codon of every genetic code, and compare the
answers -- refusals included.
"""

from __future__ import annotations

import warnings

from . import _core
from ._codon_table import (
    AMBIGUOUS_DNA_LETTERS,
    AMBIGUOUS_RNA_LETTERS,
    TranslationError,
)

__all__ = ["translate"]

#: The reference's `pos_stop`: what a codon means when every reading of it
#: agrees except that some readings are stops.  Not a parameter of the public
#: function -- `Seq.translate` does not expose it either -- but the private
#: entry point keeps it, because the reference's private one has it.
_POSSIBLE_STOP = "X"

#: The byte that says "not one of this genetic code's own letters".  Base codes
#: are 0..3, so bit 7 is free, and `(a | b | c) & 0x80` is the kernel's entire
#: input validation.
_NOT_A_BASE = 0x80

#: The kernel's three outcomes, in the order of `TranslateOutcome`.
_OK, _DECLINED, _EXTRA_STOP = 0, 1, 2


class BiopythonWarning(UserWarning):
    """What the reference warns with, so a caller's filter keeps working.

    `Bio.BiopythonWarning` is a plain `UserWarning` subclass with no members;
    this is the same class under the same name, because that name is the only
    thing the reference puts in front of a caller.  It is *not* the reference's
    class object, so a caller filtering on `Bio.BiopythonWarning` itself rather
    than on `Warning` will not catch this -- which is the price of a package
    that can be imported without Biopython, and the tests check the message and
    the category rather than the identity.
    """


def translate(sequence, table="Standard", stop_symbol="*", to_stop=False,
              cds=False, gap=None):
    """Translate `sequence` into protein, as ``Bio.Seq.translate`` does.

    Accepts `str` or any object exposing the buffer protocol, and returns `str`
    for `str` and `bytes` otherwise, like the rest of :mod:`biofasting.seqops`.
    The result is what ``Bio.Seq.translate`` returns for the same arguments,
    refusals included: the exception types and messages are the reference's own.

    >>> translate("ATGAAATAG")
    'MK*'
    >>> translate("ATGAAATAG", to_stop=True)
    'MK'
    >>> translate("ATGAAATAG", cds=True)
    'MK'
    >>> translate("AUGAAAUAG", table=2)
    'MK*'
    """
    if isinstance(sequence, str):
        #   `.upper()` on the `str`, not on its encoding, because that is where
        #   the reference does it: for ASCII the two agree, and for a character
        #   outside ASCII, `str.upper` may change the *length* of the string
        #   ("ß" becomes "SS") and so move every codon boundary after it.  A
        #   nucleotide sequence has no such character, which is exactly why the
        #   choice is worth making explicit rather than leaving to the
        #   encoding.
        buffer = sequence.upper().encode("latin-1")
        was_text = True
    else:
        buffer = bytes(sequence).upper()
        was_text = False

    table_object = _resolve_table(table)
    tables = _kernel_tables(table_object)

    #   Before any translation, because the reference puts it there: a table
    #   whose stops are also residues cannot answer `to_stop` at all, and the
    #   answer does not depend on the sequence.
    if tables.dual_coding:
        _refuse_or_warn(tables.dual_coding, table_object, to_stop)

    length = len(buffer)
    start, end, prefix = 0, length, ""
    if cds:
        start, end, prefix = _check_cds(buffer, length, table_object)
    elif length % 3 != 0:
        warnings.warn(
            "Partial codon, len(sequence) not a multiple of three. "
            "Explicitly trim the sequence or add trailing N before "
            "translation. This may become an error in future.",
            BiopythonWarning,
        )

    #   The reference validates `gap` even when no codon is a gap.  A gap is not
    #   a nucleotide, so no gap codon can reach the kernel anyway -- every one of
    #   them has a byte outside ACGT -- and the fallback below is where `gap` is
    #   actually spent.
    if gap is not None:
        if not isinstance(gap, str):
            raise TypeError("Gap character should be a single character string.")
        if len(gap) > 1:
            raise ValueError("Gap character should be a single character string.")

    if tables.code is not None and isinstance(stop_symbol, str):
        #   `stop_symbol` has to be a `str` for the kernel to take it, and a
        #   caller who passes something else is not refused here: the reference
        #   never looks at `stop_symbol` until it has a stop codon to write, so
        #   a list is accepted for a sequence with no stop in it and only fails
        #   at the `join`.  The fallback below is where that happens too.
        #
        #   The slice is the *whole* buffer for every call that is not `cds`,
        #   and copying a 100 Mbp chromosome to hand the kernel a view of it
        #   would cost more than the translation; `start == 0 and end == length`
        #   is passed through unsliced instead.  `_core.translate` only reads.
        frame = buffer if start == 0 and end == length else buffer[start:end]
        outcome, codon, protein = _core.translate(
            frame, tables.code, tables.amino,
            stop_symbol, _POSSIBLE_STOP, to_stop, cds,
        )
        if outcome == _OK:
            protein = prefix + protein
            return protein if was_text else protein.encode("latin-1")
        if outcome == _EXTRA_STOP:
            found = buffer[start + 3 * codon:start + 3 * codon + 3]
            raise TranslationError(
                f"Extra in frame stop codon '{found.decode('latin-1')}' found."
            )
        #   Declined: nothing the kernel produced is used.

    protein = prefix + _translate_loop(
        buffer, start, end, table_object, tables.valid_letters, stop_symbol,
        to_stop, cds, gap,
    )
    return protein if was_text else protein.encode("latin-1")


# --------------------------------------------------------------------------- #
# The table argument, as the reference reads it
# --------------------------------------------------------------------------- #


def _resolve_table(table):
    """Turn the reference's three accepted spellings of a table into a table.

    The order of the two `except` clauses is the reference's and it matters:
    `int("Standard")` raises ValueError and means "a name", while a mapping
    table raises TypeError and means "neither" -- a `ValueError` raised by the
    *body* of the first clause is not caught by the second, which is how the
    "DO NOT take a character string mapping table" message reaches the caller
    instead of "Bad table argument".

    The names and the ids are memoised, and that is worth 800 ns of a 1.7 µs
    call: `int("Standard")` is half a microsecond *because* it raises, and
    "Standard" is what every caller who does not care passes.  Only `str` and
    `int` keys are kept -- a table object is not, because memoising one would
    keep a caller's own table alive for the life of the process to save a
    hash and an isinstance.
    """
    if isinstance(table, (str, int)):
        try:
            return _RESOLVED_TABLES[table]
        except KeyError:
            pass
        resolved = _resolve_table_uncached(table)
        if len(_RESOLVED_TABLES) >= 256:
            _RESOLVED_TABLES.clear()
        _RESOLVED_TABLES[table] = resolved
        return resolved
    return _resolve_table_uncached(table)


#: Memo for the two spellings above.  Bounded, because the key can be any
#: string a caller types; the built-in tables are alive anyway, so what is
#: cached is a reference to something that already exists.
_RESOLVED_TABLES = {}


def _resolve_table_uncached(table):
    from . import codon_table

    try:
        table_id = int(table)
    except ValueError:
        try:
            return codon_table.ambiguous_generic_by_name[table]
        except KeyError:
            if isinstance(table, str):
                raise ValueError(
                    "The Bio.Seq translate methods and function DO NOT "
                    "take a character string mapping table like the python "
                    "string object's translate method. "
                    "Use str(my_seq).translate(...) instead."
                ) from None
            raise TypeError("table argument must be integer or string") from None
    except (AttributeError, TypeError):
        if isinstance(table, codon_table.CodonTable):
            return table
        raise ValueError("Bad table argument") from None
    return codon_table.ambiguous_generic_by_id[table_id]


def _refuse_or_warn(dual_coding, table, to_stop):
    """The reference's two answers to a table whose stops are also residues."""
    example = dual_coding[0]
    residue = table.forward_table[example]
    if to_stop:
        raise ValueError(
            "You cannot use 'to_stop=True' with this table as it contains"
            f" {len(dual_coding)} codon(s) which can be both STOP and an"
            f" amino acid (e.g. '{example}' -> '{residue}' or STOP)."
        )
    warnings.warn(
        f"This table contains {len(dual_coding)} codon(s) which code(s) for"
        f" both STOP and an amino acid (e.g. '{example}' ->"
        f" '{residue}' or STOP). Such codons will be translated as"
        " amino acid.",
        BiopythonWarning,
    )


# --------------------------------------------------------------------------- #
# cds
# --------------------------------------------------------------------------- #


def _check_cds(buffer, length, table):
    """The reference's three cds checks, then the slice they imply.

    Returns ``(start, end, prefix)``: the codons to translate and the ``"M"``
    the first one becomes whatever it was.  The last codon is a checked stop and
    is not translated, which is why `end` is three short.
    """
    first = buffer[:3].decode("latin-1")
    if first not in table.start_codons:
        raise TranslationError(f"First codon '{first}' is not a start codon")
    if length % 3 != 0:
        raise TranslationError(f"Sequence length {length} is not a multiple of three")
    last = buffer[-3:].decode("latin-1")
    if last not in table.stop_codons:
        raise TranslationError(f"Final codon '{last}' is not a stop codon")
    return 3, length - 3, "M"


# --------------------------------------------------------------------------- #
# The fallback: the reference's own loop
# --------------------------------------------------------------------------- #


def _translate_loop(buffer, start, end, table, valid_letters, stop_symbol,
                    to_stop, cds, gap, pos_stop=_POSSIBLE_STOP):
    """`_translate_str`'s codon loop, on an already-uppercased buffer.

    This is the reference's algorithm line for line, minus the argument
    handling that `translate` above has already done: resolve the codon in the
    forward table; on failure, a stop codon ends or marks the protein, a codon
    whose every letter is a nucleotide of this table is `pos_stop`, a codon of
    three gap characters is a gap, and anything else is an error.
    """
    out = []
    forward_table = table.forward_table
    stop_codons = table.stop_codons
    count = end - start
    for i in range(start, start + count - count % 3, 3):
        codon = buffer[i:i + 3].decode("latin-1")
        try:
            out.append(forward_table[codon])
        except (KeyError, TranslationError):
            if codon in stop_codons:
                if cds:
                    raise TranslationError(
                        f"Extra in frame stop codon '{codon}' found."
                    ) from None
                if to_stop:
                    break
                out.append(stop_symbol)
            elif valid_letters.issuperset(codon):
                out.append(pos_stop)
            elif gap is not None and codon == gap * 3:
                out.append(gap)
            else:
                raise TranslationError(f"Codon '{codon}' is invalid") from None
    return "".join(out)


# --------------------------------------------------------------------------- #
# The kernel's two tables, derived from a genetic code and cached on it
# --------------------------------------------------------------------------- #


class _Tables:
    """Everything `translate` needs about a table that does not depend on the
    sequence, the options, or the call."""

    __slots__ = ("code", "amino", "dual_coding", "valid_letters")

    def __init__(self, code, amino, dual_coding, valid_letters):
        self.code = code
        self.amino = amino
        self.dual_coding = dual_coding
        self.valid_letters = valid_letters


def _kernel_tables(table):
    """The kernel's `code` and `amino` bytes for `table`, or `code = None`.

    Built from the table object itself -- its own `forward_table`, its own
    `stop_codons`, its own alphabet -- so a table the reference would refuse a
    codon for is a table whose `amino` entry says so, and the kernel declines
    instead of answering.
    """
    #   `object.__getattribute__`, not `getattr`, and this is not a
    #   micro-optimisation.  An `AmbiguousCodonTable` defines `__getattr__` to
    #   *delegate to the unambiguous table it was built from* -- so a plain
    #   `getattr(table, name, None)` on a fresh ambiguous table finds the
    #   unambiguous table's cache and answers with it.  The two tables have the
    #   same 64 codons, so the values look plausible; what differs is the
    #   alphabet, so `ambiguous_dna_by_id[1]` translated as if U were
    #   unambiguously invalid and `"NNN"` was a TranslationError instead of
    #   `"X"`.  It showed up only when the unambiguous table had been used
    #   first, which is why one probe missed it and the sweep did not.
    try:
        cached = object.__getattribute__(table, "_biofasting_translation_tables")
    except AttributeError:
        cached = None
    if cached is not None:
        return cached
    tables = _build_kernel_tables(table)
    try:
        object.__setattr__(table, "_biofasting_translation_tables", tables)
    except AttributeError:
        #   A table that will not hold an attribute (`__slots__`, a proxy, a
        #   frozen instance) is still translatable; it just rebuilds these two
        #   tables on every call, which costs a few microseconds and is never
        #   wrong.
        pass
    return tables


def _build_kernel_tables(table):
    alphabet = table.nucleotide_alphabet
    if alphabet is not None:
        valid_letters = set(alphabet.upper())
    else:
        #   The reference's own worst case, and the generic tables really do
        #   take it: a table with no alphabet may be read as DNA or as RNA.
        valid_letters = set(AMBIGUOUS_DNA_LETTERS + AMBIGUOUS_RNA_LETTERS)

    stop_codons = table.stop_codons
    forward_table = table.forward_table
    dual_coding = tuple(codon for codon in stop_codons if codon in forward_table)

    #   Which four letters this table spells its codons with.  A DNA table has
    #   no U and an RNA table has no T, and both refuse the other's codons -- so
    #   this is not a formatting choice, it decides which codons exist.
    if "T" in valid_letters:
        letters = "ACGT"
    elif "U" in valid_letters:
        letters = "ACGU"
    else:
        return _Tables(None, None, dual_coding, valid_letters)

    amino = bytearray(64)
    for first in range(4):
        for second in range(4):
            for third in range(4):
                codon = letters[first] + letters[second] + letters[third]
                index = (first << 4) | (second << 2) | third
                try:
                    residue = forward_table[codon]
                except (KeyError, TranslationError):
                    if codon in stop_codons:
                        amino[index] = 1  # kTranslateStop
                    elif valid_letters.issuperset(codon):
                        amino[index] = 2  # kTranslatePossibleStop
                    else:
                        amino[index] = 0  # kTranslateUnclassified
                else:
                    if not isinstance(residue, str) or len(residue) != 1:
                        #   The reference would join this in as it stands; a
                        #   kernel that wrote one byte could not.
                        return _Tables(None, None, dual_coding, valid_letters)
                    amino[index] = ord(residue)

    if letters == "ACGT" and "U" in valid_letters:
        #   Both spellings are legal here (the generic tables), and they may
        #   share a code only if they agree codon for codon.  Checked, not
        #   assumed: the two spellings are two tables' worth of answers and
        #   nothing in the format says they match.
        rna = _build_alphabet_table(forward_table, stop_codons, valid_letters, "ACGU")
        if rna is None or rna != amino:
            return _Tables(None, None, dual_coding, valid_letters)

    code = bytearray([_NOT_A_BASE]) * 256
    for value, letter in enumerate(letters):
        code[ord(letter)] = value
        code[ord(letter.lower())] = value
    if letters == "ACGT" and "U" in valid_letters:
        code[ord("U")] = code[ord("u")] = 3
    return _Tables(bytes(code), bytes(amino), dual_coding, valid_letters)


def _build_alphabet_table(forward_table, stop_codons, valid_letters, letters):
    """The same 64 entries, spelled with `letters` -- for the DNA/RNA check."""
    table = bytearray(64)
    for first in range(4):
        for second in range(4):
            for third in range(4):
                codon = letters[first] + letters[second] + letters[third]
                index = (first << 4) | (second << 2) | third
                try:
                    residue = forward_table[codon]
                except (KeyError, TranslationError):
                    if codon in stop_codons:
                        table[index] = 1
                    elif valid_letters.issuperset(codon):
                        table[index] = 2
                    else:
                        table[index] = 0
                else:
                    if not isinstance(residue, str) or len(residue) != 1:
                        return None
                    table[index] = ord(residue)
    return table
