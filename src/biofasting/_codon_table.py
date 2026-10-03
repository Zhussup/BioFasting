"""The genetic code: the classes, and the derivations that read the table.

This is `Bio.Data.CodonTable` with the 27 genetic codes taken out of it.  The
reference module is 1,308 lines, and roughly 900 of them are 27 copies of

    register_ncbi_table(
        name="Standard",
        alt_name="SGC0",
        id=1,
        table={
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            ...
        },
        stop_codons=["TAA", "TAG", "TGA"],
        start_codons=["TTG", "CTG", "ATG"],
    )

-- the table written out longhand, codon by codon.  It does not have to be
written out at all: NCBI publishes a genetic code as a 64-character string of
amino acids in the codon order TTT, TTC, TTA, TTG, TCT, ... GGG, with `*` at the
stops, and that string is exactly what those 27 dictionaries say.  So the data
is the string (`_codon_tables_data.py`, 27 rows) and this is the reader.

Everything here is derivation, not data, and it is all pure: no imports, so
`tools/gen_codon_tables.py` can load this file by path and use it to check the
table it is about to write, without importing the table it wrote last time.
That is the same arrangement as `_restriction_sites.py`, and for the same
reason: one derivation, two callers, and the generator cannot pass a check the
library would fail.

The one place a 64-character string cannot say what the reference says is the
codon that is both a stop and an amino acid -- TGA in the karyorelictid
ciliates, TAA and TAG in *Condylostoma* and *Blastocrithidia*, where the release
factor and a tRNA compete for the same triplet.  Those are the `dual` field of
a row, and they are a fact about the organisms rather than an artifact of the
encoding.
"""

from __future__ import annotations

#   The order NCBI prints a genetic code in: first base varies slowest, third
#   fastest.  Every NCBI table string is in this order, which is what makes a
#   string enough to carry the whole code.
CODON_ORDER = tuple(a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG")

#   The alphabets and the ambiguity codes, to the letter of the IUPAC standard.
#   `X` is not part of that standard; it is here because the reference has it,
#   meaning "any base", and a table built without it would answer differently.
UNAMBIGUOUS_DNA_LETTERS = "GATC"
UNAMBIGUOUS_RNA_LETTERS = "GAUC"
AMBIGUOUS_DNA_LETTERS = "GATCRYWSMKHBVDN"
AMBIGUOUS_RNA_LETTERS = "GAUCRYWSMKHBVDN"
PROTEIN_LETTERS = "ACDEFGHIKLMNPQRSTVWY"
EXTENDED_PROTEIN_LETTERS = "ACDEFGHIKLMNPQRSTVWYBXZJUO"

AMBIGUOUS_DNA_VALUES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "M": "AC", "R": "AG", "W": "AT", "S": "CG", "Y": "CT", "K": "GT",
    "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "X": "GATC", "N": "GATC",
}
AMBIGUOUS_RNA_VALUES = {
    "A": "A", "C": "C", "G": "G", "U": "U",
    "M": "AC", "R": "AG", "W": "AU", "S": "CG", "Y": "CU", "K": "GU",
    "V": "ACG", "H": "ACU", "D": "AGU", "B": "CGU", "X": "GAUC", "N": "GAUC",
}
#   What an ambiguous protein letter stands for, used to name the residue a
#   codon whose bases are ambiguous might encode: `B` is D or N, `Z` is Q or E,
#   `J` is I or L, and `X` is any of them.
EXTENDED_PROTEIN_VALUES = {
    "A": "A", "B": "ND", "C": "C", "D": "D", "E": "E", "F": "F", "G": "G",
    "H": "H", "I": "I", "J": "IL", "K": "K", "L": "L", "M": "M", "N": "N",
    "O": "O", "P": "P", "Q": "Q", "R": "R", "S": "S", "T": "T", "U": "U",
    "V": "V", "W": "W", "X": "ACDEFGHIKLMNPQRSTVWY", "Y": "Y", "Z": "QE",
}


class TranslationError(Exception):
    """A codon that does not have one translation, so there is none."""


def forward_table(aa_string, dual=()) -> dict:
    """The codon to amino acid mapping a genetic code says.

    `aa_string` is either NCBI's 64 characters in `CODON_ORDER`, or the codon
    dictionary itself -- the reference's `register_ncbi_table` takes the latter,
    and accepting both is what lets that same call work here.  `dual` lists the
    codons that are a stop *and* code for an amino acid, which no string of this
    shape can express and which are therefore stated separately.
    """
    if isinstance(aa_string, dict):
        table = dict(aa_string)
    else:
        table = {
            codon: amino
            for codon, amino in zip(CODON_ORDER, aa_string)
            if amino != "*"
        }
    table.update(dual)
    return table


def make_back_table(table: dict, default_stop_codon):
    """The amino acid to codon mapping, one codon per amino acid.

    Only one codon is returned per amino acid, chosen by sort order -- the
    reference is explicit that the choice is arbitrary and that the result must
    not be cached by a caller, since a different codon of the same amino acid is
    as correct as this one.  The sorted order is what makes it reproducible:
    without it the answer would follow the dictionary's insertion order.

    The *last* codon in sort order wins, not the first, because the loop
    overwrites as it goes:

        >>> make_back_table({"TTT": "F", "TTC": "F", "ATG": "M"}, "TAA")
        {'M': 'ATG', 'F': 'TTT', None: 'TAA'}
    """
    back_table = {}
    for key in sorted(table):
        back_table[table[key]] = key
    back_table[None] = default_stop_codon
    return back_table


class CodonTable:
    """A codon table, or genetic code."""

    forward_table: dict = {}
    back_table: dict = {}
    start_codons: list = []
    stop_codons: list = []

    def __init__(
        self,
        nucleotide_alphabet=None,
        protein_alphabet=None,
        forward_table: dict = forward_table,
        back_table: dict = back_table,
        start_codons: list = start_codons,
        stop_codons: list = stop_codons,
    ) -> None:
        self.nucleotide_alphabet = nucleotide_alphabet
        self.protein_alphabet = protein_alphabet
        self.forward_table = forward_table
        self.back_table = back_table
        self.start_codons = start_codons
        self.stop_codons = stop_codons

    def __str__(self):
        """The code as the grid of 64 codons the reference prints.

        Note that `id` and `names` are read off the instance and are not
        attributes of this class: a bare `CodonTable` has neither, and asking
        for its string is an `AttributeError` here exactly as it is upstream.
        The grid uses the DNA letters whenever the alphabet has a `T` in it and
        the RNA letters otherwise -- including for a generic table, whose
        alphabet is `None`.
        """
        if self.id:
            answer = "Table %i" % self.id
        else:
            answer = "Table ID unknown"
        if self.names:
            answer += " " + ", ".join([x for x in self.names if x])

        letters = self.nucleotide_alphabet
        if letters is not None and "T" in letters:
            letters = "TCAG"
        else:
            letters = "UCAG"

        answer += "\n\n"
        answer += "  |" + "|".join(f"  {c2}      " for c2 in letters) + "|"
        answer += "\n--+" + "+".join("---------" for c2 in letters) + "+--"
        for c1 in letters:
            for c3 in letters:
                line = c1 + " |"
                for c2 in letters:
                    codon = c1 + c2 + c3
                    line += f" {codon}"
                    if codon in self.stop_codons:
                        line += " Stop|"
                    else:
                        try:
                            amino = self.forward_table[codon]
                        except KeyError:
                            amino = "?"
                        except TranslationError:
                            amino = "?"
                        if codon in self.start_codons:
                            line += f" {amino}(s)|"
                        else:
                            line += f" {amino}   |"
                line += " " + c3
                answer += "\n" + line
            answer += "\n--+" + "+".join("---------" for c2 in letters) + "+--"
        return answer


class NCBICodonTable(CodonTable):
    """A genetic code as NCBI publishes it, over T or over U."""

    nucleotide_alphabet = None
    protein_alphabet = PROTEIN_LETTERS

    def __init__(self, id, names, table, start_codons, stop_codons):
        self.id = id
        self.names = names
        self.forward_table = table
        self.back_table = make_back_table(table, stop_codons[0])
        self.start_codons = start_codons
        self.stop_codons = stop_codons

    def __repr__(self):
        return f"{self.__class__.__name__}(id={self.id!r}, names={self.names!r}, ...)"


class NCBICodonTableDNA(NCBICodonTable):
    """A genetic code over unambiguous DNA."""

    nucleotide_alphabet = UNAMBIGUOUS_DNA_LETTERS


class NCBICodonTableRNA(NCBICodonTable):
    """A genetic code over unambiguous RNA."""

    nucleotide_alphabet = UNAMBIGUOUS_RNA_LETTERS


def list_possible_proteins(codon, forward_table: dict, ambiguous_nucleotide_values: dict):
    """Every amino acid an ambiguous codon could encode.

    Raises `KeyError` if the codon can only be a stop, and `TranslationError` if
    it could be both a stop and an amino acid -- which is not a translation at
    all, and is the difference between "we do not know" and "there is no
    answer".  The standard code is the one being asked here, so `GCN` is
    alanine four ways over and `TRR` is every stop plus tryptophan:

        >>> table = forward_table(
        ...     "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
        >>> list_possible_proteins("GCN", table, AMBIGUOUS_DNA_VALUES)
        ['A']
        >>> list_possible_proteins("TAR", table, AMBIGUOUS_DNA_VALUES)
        Traceback (most recent call last):
        ...
        KeyError: 'TAR'
    """
    c1, c2, c3 = codon
    x1 = ambiguous_nucleotide_values[c1]
    x2 = ambiguous_nucleotide_values[c2]
    x3 = ambiguous_nucleotide_values[c3]
    possible = {}
    stops = []
    for y1 in x1:
        for y2 in x2:
            for y3 in x3:
                try:
                    possible[forward_table[y1 + y2 + y3]] = 1
                except KeyError:
                    stops.append(y1 + y2 + y3)
    if stops:
        if possible:
            raise TranslationError(
                f"ambiguous codon {codon!r} codes for both proteins and stop codons"
            )
        raise KeyError(codon)
    return list(possible)


def list_ambiguous_codons(codons, ambiguous_nucleotide_values: dict):
    """The stop (or start) codon list, extended with the ambiguous codons.

    `['TAG', 'TAA']` becomes `['TAG', 'TAA', 'TAR']`, and `['TGA', 'TAA', 'TAG']`
    becomes `['TGA', 'TAA', 'TAG', 'TAR', 'TRA']`.  Not every ambiguity that
    covers the listed codons is added: `['TAG', 'TGA']` does not gain `'TRR'`,
    because `TRR` also covers `TAA` and `TGG` and so could mean something that
    is not a stop.  A candidate is kept only when every codon it covers is
    already in the list.

        >>> list_ambiguous_codons(["TAG", "TAA"], AMBIGUOUS_DNA_VALUES)
        ['TAG', 'TAA', 'TAR']
        >>> list_ambiguous_codons(["TAG", "TGA"], AMBIGUOUS_DNA_VALUES)
        ['TAG', 'TGA']
        >>> list_ambiguous_codons(["UAG", "UAA"], AMBIGUOUS_RNA_VALUES)
        ['UAG', 'UAA', 'UAR']
    """
    c1_list = sorted(
        letter
        for letter, meanings in ambiguous_nucleotide_values.items()
        if {codon[0] for codon in codons}.issuperset(set(meanings))
    )
    c2_list = sorted(
        letter
        for letter, meanings in ambiguous_nucleotide_values.items()
        if {codon[1] for codon in codons}.issuperset(set(meanings))
    )
    c3_list = sorted(
        letter
        for letter, meanings in ambiguous_nucleotide_values.items()
        if {codon[2] for codon in codons}.issuperset(set(meanings))
    )
    candidates = []
    for c1 in c1_list:
        for c2 in c2_list:
            for c3 in c3_list:
                codon = c1 + c2 + c3
                if codon not in candidates and codon not in codons:
                    candidates.append(codon)
    answer = codons[:]
    for ambig_codon in candidates:
        wanted = True
        for codon in [
            c1 + c2 + c3
            for c1 in ambiguous_nucleotide_values[ambig_codon[0]]
            for c2 in ambiguous_nucleotide_values[ambig_codon[1]]
            for c3 in ambiguous_nucleotide_values[ambig_codon[2]]
        ]:
            if codon not in codons:
                wanted = False
                continue
        if wanted:
            answer.append(ambig_codon)
    return answer


class AmbiguousForwardTable:
    """Translation of an ambiguous codon, via the residues every reading shares.

    `table[codon]` returns one letter when there is exactly one answer.  When a
    codon can encode several residues, the answer is the most specific
    ambiguous protein letter that covers all of them -- `B` for a codon reading
    as D or N -- and it is the *smallest* such letter, ties broken
    alphabetically by the number of residues it stands for and then by the
    letter itself.  When no letter covers every reading, the codon has no
    translation and `TranslationError` is raised, which is a different outcome
    from `KeyError` and has to stay different.
    """

    def __init__(self, forward_table: dict, ambiguous_nucleotide: dict, ambiguous_protein: dict):
        self.forward_table = forward_table
        self.ambiguous_nucleotide = ambiguous_nucleotide
        self.ambiguous_protein = ambiguous_protein

        inverted = {}
        for name, val in ambiguous_protein.items():
            for c in val:
                x = inverted.get(c, {})
                x[name] = 1
                inverted[c] = x
        for name, val in inverted.items():
            inverted[name] = list(val)
        self._inverted = inverted

        self._cache = {}

    def __contains__(self, codon) -> bool:
        try:
            self.__getitem__(codon)
            return True
        except (KeyError, TranslationError):
            return False

    def get(self, codon, failobj=None):
        try:
            return self.__getitem__(codon)
        except KeyError:
            return failobj

    def __getitem__(self, codon):
        try:
            x = self._cache[codon]
        except KeyError:
            pass
        else:
            if x is TranslationError:
                raise TranslationError(codon)
            if x is KeyError:
                raise KeyError(codon)
            return x
        try:
            x = self.forward_table[codon]
            self._cache[codon] = x
            return x
        except KeyError:
            pass

        try:
            possible = list_possible_proteins(
                codon, self.forward_table, self.ambiguous_nucleotide
            )
        except KeyError:
            self._cache[codon] = KeyError
            raise KeyError(codon) from None
        except TranslationError:
            self._cache[codon] = TranslationError
            raise TranslationError(codon)
        assert len(possible) > 0, "unambiguous codons must code"

        if len(possible) == 1:
            self._cache[codon] = possible[0]
            return possible[0]

        ambiguous_possible = {}
        for amino in possible:
            for term in self._inverted[amino]:
                ambiguous_possible[term] = ambiguous_possible.get(term, 0) + 1

        n = len(possible)
        possible = [
            amino for amino, val in ambiguous_possible.items() if val == n
        ]
        if len(possible) == 0:
            self._cache[codon] = TranslationError
            raise TranslationError(codon)

        possible.sort(key=lambda x: (len(self.ambiguous_protein[x]), x))
        x = possible[0]
        self._cache[codon] = x
        return x


class AmbiguousCodonTable(CodonTable):
    """A genetic code over ambiguous nucleotides.

    The reference's two `XXX` comments are worth repeating, because they are
    still true: `start_codons` and `stop_codons` here are the *ambiguous*
    extensions of the unambiguous lists, which is not the same set as the
    ambiguous codons that actually start or stop.  The extension is what the
    reference computes and what callers see; changing it here would make this
    library disagree with the one it is checked against.
    """

    def __init__(
        self,
        codon_table,
        ambiguous_nucleotide_alphabet,
        ambiguous_nucleotide_values,
        ambiguous_protein_alphabet,
        ambiguous_protein_values,
    ):
        CodonTable.__init__(
            self,
            ambiguous_nucleotide_alphabet,
            ambiguous_protein_alphabet,
            AmbiguousForwardTable(
                codon_table.forward_table,
                ambiguous_nucleotide_values,
                ambiguous_protein_values,
            ),
            codon_table.back_table,
            list_ambiguous_codons(codon_table.start_codons, ambiguous_nucleotide_values),
            list_ambiguous_codons(codon_table.stop_codons, ambiguous_nucleotide_values),
        )
        self._codon_table = codon_table

    def __getattr__(self, name):
        return getattr(self._codon_table, name)


def build(
    id,
    name,
    alt_name,
    aa_string: str,
    start_codons,
    stop_codons,
    dual=(),
):
    """The six tables the reference registers for one genetic code.

    Returns `(names, dna, rna, generic, ambiguous_dna, ambiguous_rna,
    ambiguous_generic)`.  The names are the reference's own split of the
    description: the separators `"; "`, `", "` and `" and "` all cut a name, and
    the alternative name is appended -- `None` included, which is why table 27's
    name list really does contain a `None` and why `__str__` filters false names
    out rather than printing them.
    """
    names = [
        x.strip() for x in name.replace(" and ", "; ").replace(", ", "; ").split("; ")
    ]
    #   The reference's codon lists are lists, and `list_ambiguous_codons`
    #   copies one with `[:]`, so a tuple here would raise on the first append.
    start_codons = list(start_codons)
    stop_codons = list(stop_codons)

    table = forward_table(aa_string, dual)
    dna = NCBICodonTableDNA(id, names + [alt_name], table, start_codons, stop_codons)
    ambig_dna = AmbiguousCodonTable(
        dna,
        AMBIGUOUS_DNA_LETTERS,
        AMBIGUOUS_DNA_VALUES,
        EXTENDED_PROTEIN_LETTERS,
        EXTENDED_PROTEIN_VALUES,
    )

    rna_table = {}
    generic_table = {}
    for codon, val in table.items():
        generic_table[codon] = val
        codon = codon.replace("T", "U")
        generic_table[codon] = val
        rna_table[codon] = val

    rna_start_codons = []
    generic_start_codons = []
    for codon in start_codons:
        generic_start_codons.append(codon)
        if "T" in codon:
            generic_start_codons.append(codon.replace("T", "U"))
        rna_start_codons.append(codon.replace("T", "U"))

    rna_stop_codons = []
    generic_stop_codons = []
    for codon in stop_codons:
        generic_stop_codons.append(codon)
        if "T" in codon:
            generic_stop_codons.append(codon.replace("T", "U"))
        rna_stop_codons.append(codon.replace("T", "U"))

    generic = NCBICodonTable(
        id, names + [alt_name], generic_table, generic_start_codons, generic_stop_codons
    )

    merged_values = dict(AMBIGUOUS_RNA_VALUES)
    merged_values["T"] = "U"
    ambig_generic = AmbiguousCodonTable(
        generic,
        None,
        merged_values,
        EXTENDED_PROTEIN_LETTERS,
        EXTENDED_PROTEIN_VALUES,
    )

    rna = NCBICodonTableRNA(
        id, names + [alt_name], rna_table, rna_start_codons, rna_stop_codons
    )

    ambig_rna = AmbiguousCodonTable(
        rna,
        AMBIGUOUS_RNA_LETTERS,
        AMBIGUOUS_RNA_VALUES,
        EXTENDED_PROTEIN_LETTERS,
        EXTENDED_PROTEIN_VALUES,
    )

    if alt_name is not None:
        names.append(alt_name)

    return names, dna, rna, generic, ambig_dna, ambig_rna, ambig_generic
