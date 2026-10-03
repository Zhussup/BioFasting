"""Restriction enzymes: a table, and one engine that reads it.

`Bio.Restriction` is 1,088 enzyme classes, generated at import time from a
dictionary that says, per enzyme, a recognition site, two or four cut
coordinates, an overhang, a frequency and a list of suppliers.  The classes
carry no behaviour of their own -- what an enzyme *does* is decided by which
mixin classes the generator combined, and there are only 51 distinct
combinations of them.  So this module says the same thing the way the data
wants to be said: a table (`_restriction_data.ENZYMES`, 1,088 rows) plus one
`Enzyme` class that reads the flags off a row.

That is the `data` verdict from `inventory/TRIAGE.md` put into practice --
"becomes data + a parser; never ported as code" -- and it is 27,000 lines of
generated Python that nobody has to maintain.

Where the reference makes you import a name, this makes you import an object:

    >>> from biofasting.restriction import EcoRI, SnaI
    >>> EcoRI.site, EcoRI.size
    ('GAATTC', 6)
    >>> EcoRI.search("TTTGAATTCAAA")
    [5]
    >>> EcoRI.catalyze("TTTGAATTCAAA")
    ('TTTG', 'AATTCAAA')

Two deliberate differences from the reference, both of them about types rather
than about results.  Sequences here are `str` (or `bytes`, or anything with
`__bytes__` such as `Bio.Seq.Seq`), not `Bio.Seq.Seq` only, because a
restriction digest is not a reason to require Biopython.  And where the
reference splits a batch on its mixin classes -- `batch.split(Blunt)` -- this
splits on the class *name*, so `batch.split(Blunt)` still works if `Blunt` is
imported from `Bio.Restriction`, and `batch.split("blunt")` works without it.

A third, smaller one: `Analysis.with_name()` warns about an unknown enzyme with
a plain `UserWarning` where the reference raises `BiopythonWarning`.  The latter
is a subclass of the former, so code that catches `UserWarning` -- which is
what the reference's own documentation recommends -- behaves identically, and
importing `Bio` to name a warning class would give the module a dependency it
does not otherwise have.

Everything else, including the awkward parts -- the `None` returned by
`Enzyme.__add__` when the right operand is a batch, `compatible_end()` ignoring
the batch it was handed, the `TypeError` that `Analysis.with_site_size()` raises
whenever it is given a result dictionary, and the `None` that `_boundaries()`
returns for an empty region -- behaves as the reference does, because a library
that is a surprise in a different place each time is worse than one that is a
surprise in the same places.  All of them are pinned by name in
`tests/test_restriction.py`.
"""

from __future__ import annotations

import re
import string
import warnings
from typing import Any, Callable, Iterator

from ._print_format import PrintFormat
from ._restriction_data import ENZYMES, SUPPLIERS
from ._restriction_sites import compsite

__all__ = [
    "FormattedSeq",
    "Enzyme",
    "RestrictionBatch",
    "Analysis",
    "PrintFormat",
    "AllEnzymes",
    "CommOnly",
    "NonComm",
    "get_enzyme",
    "enzymes_of_supplier",
    "supplier_codes",
]


# --------------------------------------------------------------------------- #
#   The table's translation, and the sequence shape `search` works on
# --------------------------------------------------------------------------- #

#   Bytes -> upper-case ASCII letter, 0 for anything that is not a letter.  The
#   reference builds the same table from the same idea: it does not validate
#   the sequence against an alphabet, it upper-cases letters and blanks
#   everything else, then deletes the blanks that were whitespace or digits.
#   Whatever is left standing as a 0 is a character the module will not accept.
_TABLE = bytes(
    byte if 65 <= byte <= 90 else byte - 32 if 97 <= byte <= 122 else 0
    for byte in range(256)
)
_DELETE = string.whitespace.encode() + string.digits.encode()


def _identity(text: str) -> Any:
    return text


def _as_bytes(text: str) -> bytes:
    return text.encode("ascii")


def _is_sequence(seq: Any) -> bool:
    """Whether `_as_bytes_like` has a path for this object."""
    return isinstance(seq, (str, bytes, bytearray, memoryview)) or hasattr(seq, "__bytes__")


def _as_bytes_like(seq: Any) -> tuple[bytes, Callable[[str], Any]]:
    """The sequence as bytes, and how to turn a fragment back into its type.

    `catalyze` returns fragments of the same kind it was given, which is the
    reference's `self.klass(...)` in one line instead of a class attribute.
    """
    if isinstance(seq, str):
        try:
            return seq.encode("ascii"), _identity
        except UnicodeEncodeError:
            #   The message below is the one the reference produces for a
            #   non-IUPAC character, and this is the same situation reached one
            #   step earlier.  A caller that catches TypeError by text should
            #   not have to know which of the two paths it came down.
            raise TypeError(f"Invalid character found in {seq}") from None
    if isinstance(seq, (bytes, bytearray, memoryview)):
        return bytes(seq), _as_bytes
    if hasattr(seq, "__bytes__"):
        klass = type(seq)
        return bytes(seq), (lambda text: klass(text))
    raise TypeError(
        f"expected a string, bytes, or an object with __bytes__ "
        f"(such as Bio.Seq.Seq), got {type(seq)}"
    )


class FormattedSeq:
    """A sequence in the shape a restriction search wants.

    Three changes, all of them the reference's: letters are upper-cased,
    whitespace and digits are deleted, and one space is pushed in front so that
    position 1 is the first base.  A character that is neither a letter, a
    digit nor whitespace is refused rather than skipped, because a sequence
    with a `*` in it is a sequence somebody should look at.

    The original's case is remembered and restored on the way out, so
    lowercase input produces lowercase fragments -- again, as the reference
    does.
    """

    def __init__(self, seq: Any, linear: bool = True) -> None:
        if isinstance(seq, FormattedSeq):
            self.lower = seq.lower
            self.data = seq.data
            self.linear = seq.linear
            self._wrap = seq._wrap
            return
        raw, wrap = _as_bytes_like(seq)
        self._wrap = wrap
        self.linear = linear
        #   `str.islower()` semantics, on bytes: at least one cased character,
        #   and every cased character lowercase.  'ACGT' is not lowercase,
        #   'acgt' is, '123' is not (no cased characters at all).
        self.lower = raw.lower() == raw and raw.upper() != raw
        data = raw.translate(_TABLE, _DELETE)
        if 0 in data:
            raise TypeError(f"Invalid character found in {raw.decode('ascii', 'replace')}")
        self.data = " " + data.decode("ASCII")

    def __len__(self) -> int:
        return len(self.data) - 1

    def __repr__(self) -> str:
        return f"FormattedSeq({self.data[1:]!r}, linear={self.linear!r})"

    def __getitem__(self, index):
        text = self.data[index]
        if self.lower:
            text = text.lower()
        return self._wrap(text)

    def is_linear(self) -> bool:
        return self.linear

    def finditer(self, pattern, size: int):
        """Every site match in the sequence, with the pattern's group accessor.

        Returns `(start, group)` pairs where `group` is the match's own
        `group` method -- the reference's shape, kept because it is what lets a
        non-palindromic enzyme ask "which strand did this match?" by calling
        `group("t")` and getting `None` for the strand that did not match.

        On a circular sequence the sequence is repeated far enough to see the
        sites that straddle the origin, which is why `size` is needed: it is
        how much of the front has to be repeated at the back.
        """
        data = self.data if self.is_linear() else self.data + self.data[1:size]
        return [(match.start(), match.group) for match in re.finditer(pattern, data)]


# --------------------------------------------------------------------------- #
#   One enzyme
# --------------------------------------------------------------------------- #

#   The reference's `matching` table: for each IUPAC code, the codes it is
#   compatible with.  Copied verbatim rather than derived, because it is not
#   derivable from the code's own expansion -- 'A' matching 'R' matching 'N' is
#   a statement about which overhangs can ligate, and it is the reference's
#   answer to that question, not a recomputation of it.
MATCHING = {
    "A": "ARWMHVDN",
    "C": "CYSMHBVN",
    "G": "GRSKBVDN",
    "T": "TYWKHBDN",
    "R": "ABDGHKMNSRWV",
    "Y": "CBDHKMNSTWVY",
    "W": "ABDHKMNRTWVY",
    "S": "CBDGHKMNSRVY",
    "M": "ACBDHMNSRWVY",
    "K": "BDGHKNSRTWVY",
    "H": "ACBDHKMNSRTWVY",
    "B": "CBDGHKMNSRTWVY",
    "V": "ACBDGHKMNSRWVY",
    "D": "ABDGHKMNSRTWVY",
    "N": "ACBDGHKMNSRTWVY",
}


class Enzyme:
    """A restriction enzyme, as one row of the table.

    Instances are singletons, handed out by `get_enzyme()`; `biofasting.
    restriction.EcoRI` is one.  They are compared by identity, as the
    reference's classes are -- `EcoRI == EcoRI` is true, `EcoRI == EcoRV` is
    false, and `SacI != SstI` is false because `!=` asks a different question
    (do the two cut identically?) than `==` does (is this the same enzyme?).

    The polymorphic questions -- blunt or sticky, defined or ambiguous,
    methylable, commercially available -- are answered from the `types` tuple
    the generator copied out of the reference's own classification, not by
    re-deriving them from the cut coordinates.  A second opinion about a fact
    that is already written down is a way to be wrong twice.
    """

    __slots__ = (
        "name", "site", "fst5", "fst3", "scd5", "scd3", "ovhg", "ovhgseq",
        "size", "freq", "suppl", "opt_temp", "inact_temp", "uri", "rebase_id",
        "types", "_flags", "_pattern",
    )

    def __init__(self, name: str, fields: tuple) -> None:
        (
            self.site, self.fst5, self.fst3, self.scd5, self.scd3, self.ovhg,
            self.ovhgseq, self.size, self.freq, self.suppl, self.opt_temp,
            self.inact_temp, self.uri, self.rebase_id, self.types,
        ) = fields
        self.name = name
        self._flags = frozenset(self.types)
        self._pattern = None

    # ----------------------------------------------------------------- #
    #   What this enzyme is
    # ----------------------------------------------------------------- #

    def is_palindromic(self) -> bool:
        return "Palindromic" in self._flags

    def cut_once(self) -> bool:
        return "OneCut" in self._flags

    def cut_twice(self) -> bool:
        return "TwoCuts" in self._flags

    def is_blunt(self) -> bool:
        return "Blunt" in self._flags

    def is_5overhang(self) -> bool:
        return "Ov5" in self._flags

    def is_3overhang(self) -> bool:
        return "Ov3" in self._flags

    def is_defined(self) -> bool:
        return "Defined" in self._flags

    def is_ambiguous(self) -> bool:
        return "Ambiguous" in self._flags

    def is_unknown(self) -> bool:
        return "NotDefined" in self._flags

    def is_methylable(self) -> bool:
        return "Meth_Dep" in self._flags

    def is_comm(self) -> bool:
        return "Commercially_available" in self._flags

    def overhang(self) -> str:
        if self.is_blunt():
            return "blunt"
        if self.is_5overhang():
            return "5' overhang"
        if self.is_3overhang():
            return "3' overhang"
        return "unknown"

    def frequency(self) -> float:
        """One cut per this many bases, on average."""
        return self.freq

    def characteristic(self) -> tuple:
        """`(fst5, fst3, scd5, scd3, site)`, with the cuts this enzyme has.

        The reference has one `characteristic()` per cut-count mixin; the
        difference between them is only which coordinates are `None`, so the
        three collapse into one expression of the same rule.
        """
        if self.cut_twice():
            return self.fst5, self.fst3, self.scd5, self.scd3, self.site
        if self.cut_once():
            return self.fst5, self.fst3, None, None, self.site
        return None, None, None, None, self.site

    #   The reference's `cls.charac`, which its `!=` compares.  Kept as a
    #   property so there is still exactly one answer to the question.
    charac = property(characteristic)

    # ----------------------------------------------------------------- #
    #   Suppliers
    # ----------------------------------------------------------------- #

    def suppliers(self) -> None:
        """Print the enzyme's suppliers, one per line, as the reference does.

        An enzyme nobody sells prints nothing and returns None -- the
        reference's `Not_available.suppliers`, which is where the asymmetry
        with `supplier_list()` comes from.
        """
        if not self.is_comm():
            return None
        for code in self.suppl:
            print(SUPPLIERS[code] + ",")

    def supplier_list(self) -> list[str]:
        """The enzyme's suppliers by name, in the table's order."""
        return [name for code, name in SUPPLIERS.items() if code in self.suppl]

    def buffers(self, supplier: str):
        """Not implemented upstream either; the reference returns None.

        For an enzyme that is not sold, the reference raises instead, because
        there is no buffer list to be missing from.
        """
        if not self.is_comm():
            raise TypeError("Enzyme not commercially available.")
        return None

    # ----------------------------------------------------------------- #
    #   Searching
    # ----------------------------------------------------------------- #

    def pattern(self):
        """The compiled search pattern, built from the site on first use."""
        if self._pattern is None:
            self._pattern = re.compile(compsite(self.site, self.is_palindromic()))
        return self._pattern

    #   The reference's attribute name for the same object.
    compsite = property(pattern)

    def search(self, dna: Any, linear: bool = True) -> list[int]:
        """The cutting positions of this enzyme in `dna`, 1-based.

        A position is the first base of the 3' fragment: for EcoRI, which cuts
        `G^AATTC`, the site starting at base 10 cuts at base 11.
        """
        fseq = dna if isinstance(dna, FormattedSeq) else FormattedSeq(dna, linear)
        return self._search(fseq)

    def _modify(self, location: int) -> Iterator[int]:
        """Top-strand cuts for a site whose match starts at `location`."""
        if self.cut_twice():
            yield location + self.fst5
            yield location + self.scd5
        elif self.cut_once():
            yield location + self.fst5
        else:
            yield location

    def _rev_modify(self, location: int) -> Iterator[int]:
        """Bottom-strand cuts for a match at `location` (non-palindromic only)."""
        if self.cut_twice():
            yield location - self.fst3
            yield location - self.scd3
        elif self.cut_once():
            yield location - self.fst3
        else:
            yield location

    def _search(self, fseq: FormattedSeq) -> list[int]:
        results: list[int] = []
        on_minus: list[int] = []
        modify = self._modify
        rev_modify = self._rev_modify
        palindromic = self.is_palindromic()

        for start, group in fseq.finditer(self.pattern(), self.size):
            #   A palindromic site has one pattern and both strands match it;
            #   a non-palindromic one has two groups, and only one of them is
            #   ever set, so asking for the top strand's group by name is how
            #   the strand is decided.
            if palindromic or group("t"):
                results.extend(modify(start))
            else:
                on_minus.extend(rev_modify(start))
        results.extend(on_minus)

        if results:
            if not palindromic:
                #   Matches come out in site order but the two strands
                #   interleave, so the cuts have to be put back in order before
                #   `_drop` walks them from the ends.
                results.sort()
            self._drop(results, fseq)
        return results

    def _drop(self, results: list[int], fseq: FormattedSeq) -> None:
        """Remove cuts that fall off the molecule, or wrap them if circular.

        In place, because `results` is the list the caller already has.
        """
        length = len(fseq)
        if fseq.is_linear():
            if self.is_unknown():
                #   The reference's `NotDefined._drop` returns here rather than
                #   filtering: an uncharacterised enzyme has no overhang, and
                #   the comparison below would be arithmetic on None.
                return
            ovhg = self.ovhg
            results[:] = [
                cut for cut in results
                if 1 < cut <= length and 1 < cut - ovhg <= length
            ]
        else:
            #   A circular molecule has no ends, so a cut before the origin
            #   belongs at the back and one past the end belongs at the front.
            #   Both lists are in ascending order, so the first cut inside the
            #   range ends the wrap.
            for index, location in enumerate(results):
                if location < 1:
                    results[index] += length
                else:
                    break
            for index, location in enumerate(reversed(results)):
                if location > length:
                    results[-(index + 1)] -= length
                else:
                    break

    def catalyze(self, dna: Any, linear: bool = True) -> tuple:
        """The fragments left after cutting `dna` with this enzyme.

        The fragments are of the same type as the input: `str` in, `str` out;
        `Seq` in, `Seq` out.
        """
        if self.is_unknown():
            raise NotImplementedError(f"{self.name} restriction is unknown.")
        fseq = dna if isinstance(dna, FormattedSeq) else FormattedSeq(dna, linear)
        cuts = self._search(fseq)
        if not cuts:
            return (fseq[1:],)
        fragments = []
        between = len(cuts) - 1
        if fseq.is_linear():
            fragments.append(fseq[1:cuts[0]])
            if between:
                fragments += [fseq[cuts[x]:cuts[x + 1]] for x in range(between)]
            fragments.append(fseq[cuts[-1]:])
        else:
            #   Circular: the last fragment is the one that crosses the origin.
            fragments.append(fseq[cuts[-1]:] + fseq[1:cuts[0]])
            if not between:
                return tuple(fragments)
            fragments += [fseq[cuts[x]:cuts[x + 1]] for x in range(between)]
        return tuple(fragments)

    #   The reference spells it both ways and so does the rest of the world.
    catalyse = catalyze

    # ----------------------------------------------------------------- #
    #   Displaying
    # ----------------------------------------------------------------- #

    def elucidate(self) -> str:
        """The site with its cuts drawn on it, `^` on top and `_` underneath.

        >>> get_enzyme("EcoRI").elucidate()
        'G^AATT_C'
        >>> get_enzyme("KpnI").elucidate()
        'G_GTAC^C'
        >>> get_enzyme("EcoRV").elucidate()
        'GAT^_ATC'
        >>> get_enzyme("SnaI").elucidate()
        '? GTATAC ?'
        """
        if self.cut_twice():
            return "cut twice, not yet implemented sorry."
        if self.is_unknown():
            return f"? {self.site} ?"
        if self.is_ambiguous():
            return self._elucidate_ambiguous()
        return self._elucidate_defined()

    def _elucidate_defined(self) -> str:
        """The reference's `Defined.elucidate`: the cuts are inside the site."""
        f5, f3, site = self.fst5, self.fst3, self.site
        if self.is_5overhang():
            if f5 == f3 == 0:
                return "N^" + site + "_N"
            if f3 == 0:
                return site[:f5] + "^" + site[f5:] + "_N"
            return site[:f5] + "^" + site[f5:f3] + "_" + site[f3:]
        if self.is_blunt():
            return site[:f5] + "^_" + site[f5:]
        if f5 == f3 == 0:
            return "N_" + site + "^N"
        return site[:f3] + "_" + site[f3:f5] + "^" + site[f5:]

    def _elucidate_ambiguous(self) -> str:
        """The reference's `Ambiguous.elucidate`: the cuts may be outside it.

        Degenerate sites and cuts beyond the recognition sequence mean the
        drawing needs `N` padding to stay in register, and the padding differs
        on each side of the site depending on which way the cuts fell out.
        """
        f5, f3, site = self.fst5, self.fst3, self.site
        length = len(self)
        if self.is_5overhang():
            if f3 == f5 == 0:
                return "N^" + site + "_N"
            if 0 <= f5 <= length and 0 <= f3 + length <= length:
                return site[:f5] + "^" + site[f5:f3] + "_" + site[f3:]
            if 0 <= f5 <= length:
                return site[:f5] + "^" + site[f5:] + f3 * "N" + "_N"
            if 0 <= f3 + length <= length:
                return "N^" + abs(f5) * "N" + site[:f3] + "_" + site[f3:]
            if f3 + length < 0:
                return "N^" + abs(f5) * "N" + "_" + abs(length + f3) * "N" + site
            if f5 > length:
                return site + (f5 - length) * "N" + "^" + (length + f3 - f5) * "N" + "_N"
            return "N^" + abs(f5) * "N" + site + f3 * "N" + "_N"
        if self.is_blunt():
            if f5 < 0:
                return "N^_" + abs(f5) * "N" + site
            if f5 > length:
                return site + (f5 - length) * "N" + "^_N"
            raise ValueError(f"{self.name}.easyrepr() : error f5={f5}")
        if f3 == 0:
            if f5 == 0:
                return "N_" + site + "^N"
            return site + "_" + (f5 - length) * "N" + "^N"
        if 0 < f3 + length <= length and 0 <= f5 <= length:
            return site[:f3] + "_" + site[f3:f5] + "^" + site[f5:]
        if 0 < f3 + length <= length:
            return site[:f3] + "_" + site[f3:] + (f5 - length) * "N" + "^N"
        if 0 <= f5 <= length:
            return "N_" + "N" * (f3 + length) + site[:f5] + "^" + site[f5:]
        if f3 > 0:
            return site + f3 * "N" + "_" + (f5 - f3 - length) * "N" + "^N"
        if f5 < 0:
            return "N_" + abs(f3 - f5 + length) * "N" + "^" + abs(f5) * "N" + site
        return "N_" + abs(f3 + length) * "N" + site + (f5 - length) * "N" + "^N"

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __len__(self) -> int:
        return self.size

    # ----------------------------------------------------------------- #
    #   Comparing
    # ----------------------------------------------------------------- #

    def __eq__(self, other: Any) -> bool:
        #   Identity, as the reference's metaclass does: two enzymes are the
        #   same enzyme, not two enzymes that happen to cut alike.  `!=` is the
        #   question about cutting alike, which is why it is not the negation
        #   of this one.
        return self is other

    def __ne__(self, other: Any) -> bool:
        if not isinstance(other, Enzyme):
            return True
        return self.characteristic() != other.characteristic()

    __hash__ = object.__hash__

    def _order_key(self, other: Any) -> tuple:
        if not isinstance(other, Enzyme):
            raise NotImplementedError
        return (self.size, self.name), (other.size, other.name)

    def __lt__(self, other: Any) -> bool:
        a, b = self._order_key(other)
        return a < b

    def __le__(self, other: Any) -> bool:
        a, b = self._order_key(other)
        return a <= b

    def __gt__(self, other: Any) -> bool:
        a, b = self._order_key(other)
        return a > b

    def __ge__(self, other: Any) -> bool:
        a, b = self._order_key(other)
        return a >= b

    def __rshift__(self, other: Any) -> bool:
        """True for a neoschizomer: same site, different cut."""
        if not isinstance(other, Enzyme):
            return False
        return self.site == other.site and self.characteristic() != other.characteristic()

    def __mod__(self, other: Any) -> bool:
        """True if the two enzymes' ends can ligate to each other."""
        if not isinstance(other, Enzyme):
            raise TypeError(f"expected a restriction enzyme, got {type(other)} instead")
        return self._mod1(other)

    # ----------------------------------------------------------------- #
    #   Isoschizomers, and what ligates to what
    # ----------------------------------------------------------------- #

    def _mod1(self, other: "Enzyme") -> bool:
        if self.is_unknown():
            return False
        if self.is_blunt():
            return other.is_blunt()
        if self.is_5overhang():
            return self._mod2(other) if other.is_5overhang() else False
        if self.is_3overhang():
            return self._mod2(other) if other.is_3overhang() else False
        return False

    def _mod2(self, other: "Enzyme") -> bool:
        if self.is_unknown():
            #   The reference raises here rather than answering, on the grounds
            #   that an uncharacterised overhang is compatible with nobody and
            #   quietly saying so would hide a caller's mistake.
            raise ValueError(
                f"{self.name}.mod2({other.name}), {self.name} : "
                f"NotDefined. pas glop pas glop!"
            )
        if self.is_ambiguous():
            return self._mod2_ambiguous(other)
        if self.ovhgseq == other.ovhgseq:
            return True
        if other.is_ambiguous():
            return other._mod2(self)
        return False

    def _mod2_ambiguous(self, other: "Enzyme") -> bool:
        """Overhang compatibility when this enzyme's overhang is degenerate."""
        if len(self.ovhgseq) != len(other.ovhgseq):
            return False
        pattern = self.ovhgseq
        for base in pattern:
            if base == "N":
                pattern = ".".join(pattern.split("N"))
            if base in "RYWMSKHDBV":
                pattern = ("[" + MATCHING[base] + "]").join(pattern.split(base))
        return re.match(pattern, other.ovhgseq) is not None

    def is_equischizomer(self, other: "Enzyme") -> bool:
        return not self != other

    def is_neoschizomer(self, other: "Enzyme") -> bool:
        return self >> other

    def is_isoschizomer(self, other: "Enzyme") -> bool:
        return (not self != other) or self >> other

    def _related(self, batch, wanted) -> list:
        #   `_all_enzymes()` rather than the module attribute `AllEnzymes`:
        #   a module's own `__getattr__` is consulted for lookups *on* the
        #   module, not for global names inside its functions, so the name
        #   would not resolve here.
        if not batch:
            batch = _all_enzymes()
        return [x for x in batch if wanted(x)]

    def equischizomers(self, batch=None) -> list:
        found = self._related(batch, lambda x: not self != x)
        del found[found.index(self)]
        found.sort()
        return found

    def neoschizomers(self, batch=None) -> list:
        return sorted(self._related(batch, lambda x: self >> x))

    def isoschizomers(self, batch=None) -> list:
        found = self._related(batch, lambda x: (self >> x) or (not self != x))
        del found[found.index(self)]
        found.sort()
        return found

    def compatible_end(self, batch=None) -> list:
        """Every enzyme whose ends ligate to this one's.

        The reference accepts a `batch` here and then ignores it, searching
        `AllEnzymes` regardless; that is preserved, because silently honouring
        the argument would make the same call return different lists in the two
        libraries.
        """
        if self.is_unknown():
            return []
        every = _all_enzymes()
        if self.is_blunt():
            return sorted(x for x in every if x.is_blunt())
        if self.is_5overhang():
            return sorted(x for x in every if x.is_5overhang() and x % self)
        return sorted(x for x in every if x.is_3overhang() and x % self)

    # ----------------------------------------------------------------- #
    #   Operators
    # ----------------------------------------------------------------- #

    def __add__(self, other: Any):
        if isinstance(other, Enzyme):
            return RestrictionBatch([self, other])
        if isinstance(other, RestrictionBatch):
            #   The reference returns `other.add_nocheck(self)`, which is
            #   `set.add` and therefore None -- so `EcoRI + batch` is None
            #   there, and is None here.  It is kept deliberately: it is a
            #   defect, it is a defect in the library this one is verified
            #   against, and a difference nobody wrote down is worse than a
            #   defect everybody did.
            return other.add_nocheck(self)
        raise TypeError(f"expected a restriction enzyme or a batch, got {type(other)}")

    def __truediv__(self, other: Any) -> list[int]:
        return self.search(other)

    def __rtruediv__(self, other: Any) -> list[int]:
        return self.search(other)

    def __floordiv__(self, other: Any) -> tuple:
        return self.catalyze(other)

    def __rfloordiv__(self, other: Any) -> tuple:
        return self.catalyze(other)


# --------------------------------------------------------------------------- #
#   The registry
# --------------------------------------------------------------------------- #

_REGISTRY: dict[str, Enzyme] = {}


def get_enzyme(name: str) -> Enzyme:
    """The enzyme with this name, made on first request and kept.

    One instance per enzyme is what makes `EcoRI is EcoRI` true, and identity
    is what equality means here.
    """
    try:
        return _REGISTRY[name]
    except KeyError:
        pass
    try:
        fields = ENZYMES[name]
    except KeyError:
        raise KeyError(f"unknown restriction enzyme: {name!r}") from None
    enzyme = _REGISTRY[name] = Enzyme(name, fields)
    return enzyme


def supplier_codes() -> dict[str, str]:
    """Supplier letter -> supplier name, in the reference's order."""
    return dict(SUPPLIERS)


def enzymes_of_supplier(code: str) -> list[str]:
    """The names of the enzymes a supplier carries.

    Derived from the per-enzyme table rather than stored a second time: the
    same fact written down twice is a chance to disagree with itself.
    """
    if code not in SUPPLIERS:
        raise KeyError(f"unknown supplier code: {code!r}")
    return sorted(name for name, fields in ENZYMES.items() if code in fields[9])


# --------------------------------------------------------------------------- #
#   Batches
# --------------------------------------------------------------------------- #

#   Name -> predicate, for `RestrictionBatch.split`.  Keyed by the reference's
#   own mixin names so that `batch.split(Blunt)` keeps working for code that
#   already imports Blunt from Bio.Restriction -- only its `__name__` is read,
#   so the class itself is never needed.
FLAGS: dict[str, Callable[[Enzyme], bool]] = {
    "Blunt": Enzyme.is_blunt,
    "Ov5": Enzyme.is_5overhang,
    "Ov3": Enzyme.is_3overhang,
    "Palindromic": Enzyme.is_palindromic,
    "NonPalindromic": lambda e: not e.is_palindromic(),
    "Defined": Enzyme.is_defined,
    "Ambiguous": Enzyme.is_ambiguous,
    "NotDefined": Enzyme.is_unknown,
    "Unknown": Enzyme.is_unknown,
    "OneCut": Enzyme.cut_once,
    "TwoCuts": Enzyme.cut_twice,
    "NoCut": lambda e: not e.cut_once() and not e.cut_twice(),
    "Meth_Dep": Enzyme.is_methylable,
    "Meth_Undep": lambda e: not e.is_methylable(),
    "Commercially_available": Enzyme.is_comm,
    "Not_available": lambda e: not e.is_comm(),
}

#   The words a caller is likely to reach for instead of the reference's mixin
#   names, so `split("blunt")` works without importing anything from `Bio`.
ALIASES: dict[str, str] = {
    "blunt": "Blunt",
    "5overhang": "Ov5",
    "3overhang": "Ov3",
    "palindromic": "Palindromic",
    "defined": "Defined",
    "ambiguous": "Ambiguous",
    "unknown": "NotDefined",
    "once": "OneCut",
    "twice": "TwoCuts",
    "methylable": "Meth_Dep",
    "comm": "Commercially_available",
}


class RestrictionBatch(set):
    """A set of enzymes, with the search and the bookkeeping that needs."""

    def __init__(self, first=(), suppliers=()) -> None:
        enzymes = [self.format(x) for x in first]
        for code in suppliers:
            if code not in SUPPLIERS:
                raise KeyError(f"unknown supplier code: {code!r}")
            enzymes += [get_enzyme(name) for name in enzymes_of_supplier(code)]
        super().__init__(enzymes)
        self.mapping = dict.fromkeys(self)
        self.already_mapped = None
        self.suppliers = [code for code in suppliers if code in SUPPLIERS]

    def __str__(self) -> str:
        elements = self.elements()
        if len(elements) < 5:
            return "+".join(elements)
        return "...".join(("+".join(elements[:2]), "+".join(elements[-2:])))

    def __repr__(self) -> str:
        return f"RestrictionBatch({self.elements()})"

    def __contains__(self, other: Any) -> bool:
        try:
            other = self.format(other)
        except ValueError:
            return False
        return set.__contains__(self, other)

    def __truediv__(self, other: Any) -> dict:
        return self.search(other)

    def __rtruediv__(self, other: Any) -> dict:
        return self.search(other)

    def __add__(self, other: Any) -> "RestrictionBatch":
        new = self.__class__(self)
        new.add(other)
        return new

    def __iadd__(self, other: Any) -> "RestrictionBatch":
        self.add(other)
        return self

    def format(self, y: Any) -> Enzyme:
        """Coerce an enzyme or an enzyme's name into an `Enzyme`."""
        if isinstance(y, Enzyme):
            return y
        try:
            return get_enzyme(str(y))
        except KeyError:
            raise ValueError(f"{y.__class__} is not a restriction enzyme") from None

    def is_restriction(self, y: Any) -> bool:
        return isinstance(y, Enzyme) or str(y) in ENZYMES

    def get(self, enzyme: Any, add: bool = False) -> Enzyme:
        found = self.format(enzyme)
        if found in self:
            return found
        if add:
            self.add(found)
            return found
        raise ValueError(f"enzyme {found.name} is not in RestrictionBatch")

    def add(self, other: Any) -> None:
        return set.add(self, self.format(other))

    def add_nocheck(self, other: Any) -> None:
        return set.add(self, other)

    def remove(self, other: Any) -> None:
        return set.remove(self, self.format(other))

    def elements(self) -> list[str]:
        """The batch's enzyme names, sorted."""
        return sorted(str(enzyme) for enzyme in self)

    def as_string(self) -> list[str]:
        return [str(enzyme) for enzyme in self]

    def lambdasplit(self, func: Callable[[Enzyme], bool]) -> "RestrictionBatch":
        return RestrictionBatch(x for x in self if func(x))

    def add_supplier(self, code: str) -> None:
        self.suppliers.append(code)
        for name in enzymes_of_supplier(code):
            self.add_nocheck(get_enzyme(name))

    def current_suppliers(self) -> list[str]:
        return sorted(SUPPLIERS[code] for code in self.suppliers)

    def split(self, *classes, **wanted) -> "RestrictionBatch":
        """The enzymes matching every named class, e.g. `split("blunt")`.

        A name's default is True, so `split("blunt", Defined=False)` is the
        blunt enzymes that are not `Defined` -- the reference's convention, kept
        because it is the whole reason this is not just a filter.
        """
        predicates = []
        for klass in classes:
            given = getattr(klass, "__name__", str(klass))
            name = ALIASES.get(given, given)
            if name not in FLAGS:
                raise ValueError(
                    f"unknown enzyme class {given!r}; known names are "
                    f"{sorted(FLAGS)} (or the words {sorted(ALIASES)})"
                )
            #   Either spelling may carry the flag, so `split(Blunt, Ov5=False)`
            #   and `split("blunt", "5overhang"=...)` behave alike.
            predicates.append((FLAGS[name], wanted.get(name, wanted.get(given, True))))

        def splittest(enzyme: Enzyme) -> bool:
            for predicate, want in predicates:
                if predicate(enzyme):
                    if want:
                        continue
                    return False
                if want:
                    return False
            return True

        return RestrictionBatch(x for x in self if splittest(x))

    def search(self, dna: Any, linear: bool = True) -> dict:
        """Every enzyme's cutting positions in `dna`, as a dict.

        Cached against the formatted sequence, so a batch asked twice about the
        same molecule answers once.  The cache is keyed on the *formatted*
        sequence rather than on the caller's string, which is the one place
        this differs from the reference: it keys on `str(dna)`, so asking with
        `"acgt"` and then `"ACGT"` there re-searches while here it does not.
        The returned mapping is the same either way.
        """
        if isinstance(dna, FormattedSeq):
            fseq = dna
        elif _is_sequence(dna):
            fseq = FormattedSeq(dna, linear)
        else:
            raise TypeError(f"expected a sequence, got {type(dna)} instead")
        if self.already_mapped == (fseq.data, fseq.linear):
            return self.mapping
        self.already_mapped = (fseq.data, fseq.linear)
        self.mapping = {enzyme: enzyme.search(fseq) for enzyme in self}
        return self.mapping

    @classmethod
    def suppl_codes(cls) -> dict[str, str]:
        return supplier_codes()

    @classmethod
    def show_codes(cls) -> None:
        for code, name in sorted(supplier_codes().items()):
            print(f"{code} = {name}")


# --------------------------------------------------------------------------- #
#   Analysis: the digest's report
# --------------------------------------------------------------------------- #

_EMPTY_BATCH = RestrictionBatch()
_EMPTY_SEQUENCE = ""


class Analysis(RestrictionBatch, PrintFormat):
    """A batch searched against a sequence, with the reports that read off it.

    >>> from biofasting.restriction import Analysis, EcoRI, BamHI
    >>> analysis = Analysis(RestrictionBatch([EcoRI, BamHI]), "TTTGAATTCAAAGGATCC")
    >>> analysis.print_that()
    BamHI      :  14.
    EcoRI      :  5.
    <BLANKLINE>

    The second argument is the sequence, and searching it at construction time
    is the whole point: every method below then answers from that one search
    unless it is handed a different result dictionary.

    The filters -- `blunt`, `overhang5`, `between`, `with_name` and the rest --
    return dictionaries of the same shape `search` returns, so they compose:

        >>> len(analysis.blunt())
        0
        >>> len(analysis.with_sites())
        2
    """

    def __init__(
        self,
        restrictionbatch: RestrictionBatch = _EMPTY_BATCH,
        sequence: Any = _EMPTY_SEQUENCE,
        linear: bool = True,
    ) -> None:
        RestrictionBatch.__init__(self, restrictionbatch)
        self.rb = restrictionbatch
        self.sequence = sequence
        self.linear = linear
        #   A falsy sequence is not searched: the reference's default is an
        #   empty `Seq("")`, and that is what makes `Analysis()` a legal
        #   object rather than an error.
        if self.sequence:
            self.search(self.sequence, self.linear)

    def __repr__(self) -> str:
        return f"Analysis({self.rb!r},{self.sequence!r},{self.linear})"

    # ----------------------------------------------------------------- #
    #   Positions
    # ----------------------------------------------------------------- #

    def _sub_set(self, wanted) -> dict:
        """The results for enzymes in `wanted` -- unused upstream, kept."""
        return {key: value for key, value in self.mapping.items() if key in wanted}

    def _test_normal(self, start: int, end: int, site: int) -> bool:
        return start <= site < end

    def _test_reverse(self, start: int, end: int, site: int) -> bool:
        """`_test_normal`'s circular counterpart -- unused upstream, kept."""
        return start <= site <= len(self.sequence) or 1 <= site < end

    def _boundaries(self, start: int, end: int):
        """`(start, end, test)` for a region, or `None` if the region is empty.

        A negative or zero bound counts from the end of the sequence, as a
        Python index does.  Two equal bounds are not a region, and the
        reference returns `None` for one -- so the unpacking at every call site
        raises `TypeError: cannot unpack non-iterable NoneType object`.  It is
        kept: a caller who wrote `analysis.only_between(50, 50)` in Biopython
        gets an exception, and getting a different one here would be a new
        surprise rather than an old one.
        """
        if not isinstance(start, int):
            raise TypeError(f"expected int, got {type(start)} instead")
        if not isinstance(end, int):
            raise TypeError(f"expected int, got {type(end)} instead")
        if start < 1:
            start += len(self.sequence)
        if end < 1:
            end += len(self.sequence)
        if start < end:
            pass
        else:
            start, end = end, start
        if start < end:
            return start, end, self._test_normal

    # ----------------------------------------------------------------- #
    #   Printing
    # ----------------------------------------------------------------- #

    def format_output(self, dct=None, title: str = "", s1: str = "") -> str:
        """The report for `dct`, or for the whole analysis if it is omitted."""
        if not dct:
            dct = self.mapping
        return PrintFormat.format_output(self, dct, title, s1)

    def print_that(self, dct=None, title: str = "", s1: str = "") -> None:
        print(self.format_output(dct, title, s1))

    def change(self, **what) -> None:
        """Set the parameters of the report, re-searching where it means that.

        `sequence` re-searches, `linear` re-searches, `rb` re-initialises the
        analysis, and `NameWidth`/`ConsoleWidth` recompute the wrapping -- but
        not `linesize`, which was computed once when the class was defined.
        That is the reference's behaviour and it is why the list section and
        the non-cutting section wrap differently after a `change`.
        """
        for key, value in what.items():
            if key in ("NameWidth", "ConsoleWidth"):
                setattr(self, key, value)
                self.Cmodulo = self.ConsoleWidth % self.NameWidth
                self.PrefWidth = self.ConsoleWidth - self.Cmodulo
            elif key == "sequence":
                setattr(self, "sequence", value)
                self.search(self.sequence, self.linear)
            elif key == "rb":
                #   `self = Analysis.__init__(...)` upstream: the assignment is
                #   to a local name, but the call has already re-initialised the
                #   object, so the effect is a re-init and the discarded return
                #   value is not observable.  Kept as the call it really is.
                Analysis.__init__(self, value, self.sequence, self.linear)
            elif key == "linear":
                setattr(self, "linear", value)
                self.search(self.sequence, value)
            elif key in ("Indent", "Maxsize"):
                setattr(self, key, value)
            elif key in ("Cmodulo", "PrefWidth"):
                raise AttributeError(
                    f"To change {key}, change NameWidth and/or ConsoleWidth"
                )
            else:
                raise AttributeError(f"Analysis has no attribute {key}")

    # ----------------------------------------------------------------- #
    #   Filters
    # ----------------------------------------------------------------- #

    def full(self, linear: bool = True) -> dict:
        """Every result -- the full restriction map.  `linear` is ignored."""
        return self.mapping

    def blunt(self, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if key.is_blunt()}

    def overhang5(self, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if key.is_5overhang()}

    def overhang3(self, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if key.is_3overhang()}

    def defined(self, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if key.is_defined()}

    def with_sites(self, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if value}

    def without_site(self, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if not value}

    def with_N_sites(self, N: int, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if len(value) == N}

    def with_number_list(self, numbers, dct=None) -> dict:
        if not dct:
            dct = self.mapping
        return {key: value for key, value in dct.items() if len(value) in numbers}

    def with_name(self, names, dct=None) -> dict:
        """The results for the named enzymes.

        A name the library does not know is a warning and a deletion, not an
        error -- but the deletion happens inside the loop that is walking the
        list, so a warning on one name can silently skip the name after it.
        Kept, and pinned by a test, because the warning is the reference's
        answer to a typo and a caller may well be relying on it.
        """
        for index, enzyme in enumerate(names):
            if enzyme not in _all_enzymes():
                warnings.warn(f"no data for the enzyme: {enzyme}", UserWarning)
                del names[index]
        if not dct:
            return RestrictionBatch(names).search(self.sequence, self.linear)
        return {name: dct[name] for name in names if name in dct}

    def with_site_size(self, site_size: int, dct=None) -> dict:
        """The results for enzymes whose recognition site is this long.

        The `dct` branch of this method cannot work: it tests `key in
        site_size`, and `in` against an `int` is a `TypeError`.  The same line
        is in Biopython 1.88, so the same call fails the same way here, and
        `test_with_site_size_is_broken_upstream_in_the_same_way` says so.  It
        is not repaired: silently filtering where the reference raises would
        make this library disagree with the one it is verified against.
        """
        sites = [name for name in self if name.size == site_size]
        if not dct:
            #   Note the missing `linear`: upstream re-searches as a linear
            #   molecule here whatever the analysis was built as.
            return RestrictionBatch(sites).search(self.sequence)
        return {key: value for key, value in dct.items() if key in site_size}

    # ----------------------------------------------------------------- #
    #   Regions
    # ----------------------------------------------------------------- #

    def only_between(self, start: int, end: int, dct=None) -> dict:
        """The enzymes that cut inside the region and nowhere else."""
        start, end, test = self._boundaries(start, end)
        if not dct:
            dct = self.mapping
        filtered = dict(dct)
        for key, sites in dct.items():
            if not sites:
                del filtered[key]
                continue
            for site in sites:
                if test(start, end, site):
                    continue
                del filtered[key]
                break
        return filtered

    def between(self, start: int, end: int, dct=None) -> dict:
        """The enzymes that cut inside the region at least once."""
        start, end, test = self._boundaries(start, end)
        filtered = {}
        if not dct:
            dct = self.mapping
        for key, sites in dct.items():
            for site in sites:
                if test(start, end, site):
                    filtered[key] = sites
                    break
                continue
        return filtered

    def show_only_between(self, start: int, end: int, dct=None) -> dict:
        """`between`, with the positions outside the region dropped from the values."""
        if start <= end:
            pairs = [
                (key, [site for site in sites if start <= site <= end])
                for key, sites in self.between(start, end, dct).items()
            ]
        else:
            pairs = [
                (key, [site for site in sites if start <= site or site <= end])
                for key, sites in self.between(start, end, dct).items()
            ]
        return dict(pairs)

    def only_outside(self, start: int, end: int, dct=None) -> dict:
        """The enzymes that cut outside the region and nowhere inside it."""
        start, end, test = self._boundaries(start, end)
        if not dct:
            dct = self.mapping
        filtered = dict(dct)
        for key, sites in dct.items():
            if not sites:
                del filtered[key]
                continue
            for site in sites:
                if test(start, end, site):
                    del filtered[key]
                    break
                else:
                    continue
        return filtered

    def outside(self, start: int, end: int, dct=None) -> dict:
        """The enzymes that cut outside the region at least once."""
        start, end, test = self._boundaries(start, end)
        if not dct:
            dct = self.mapping
        filtered = {}
        for key, sites in dct.items():
            for site in sites:
                if test(start, end, site):
                    continue
                else:
                    filtered[key] = sites
                    break
        return filtered

    def do_not_cut(self, start: int, end: int, dct=None) -> dict:
        """The enzymes that do not cut inside the region -- including non-cutters."""
        if not dct:
            dct = self.mapping
        filtered = self.without_site()
        filtered.update(self.only_outside(start, end, dct))
        return filtered


# --------------------------------------------------------------------------- #
#   The names
# --------------------------------------------------------------------------- #

_BATCHES: dict[str, RestrictionBatch] = {}


def _all_enzymes() -> RestrictionBatch:
    if "AllEnzymes" not in _BATCHES:
        every = RestrictionBatch(ENZYMES)
        _BATCHES["AllEnzymes"] = every
        _BATCHES["CommOnly"] = RestrictionBatch(e for e in every if e.is_comm())
        _BATCHES["NonComm"] = RestrictionBatch(e for e in every if not e.is_comm())
    return _BATCHES["AllEnzymes"]


def __getattr__(name: str):
    """Resolve `biofasting.restriction.EcoRI` and the three batches on demand.

    Nothing is built at import: 1,088 enzyme objects and a 268 KiB table are
    what you pay for when you ask for an enzyme, not when you import the
    package that could hand you one.
    """
    if name in ENZYMES:
        return get_enzyme(name)
    if name in ("AllEnzymes", "CommOnly", "NonComm"):
        _all_enzymes()
        return _BATCHES[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(list(globals()) + list(ENZYMES))
