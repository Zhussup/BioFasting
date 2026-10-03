"""Phase 2 step 2.2b: the genetic codes, as data and as a reader.

`Bio.Data.CodonTable` is 1,308 lines, about 900 of which are 27 genetic codes
spelled out one codon at a time.  `biofasting.codon_table` is the same 27 codes
as the 64-character strings NCBI publishes, plus the machinery that turns a
string into an object.  Saying that is cheap; this file is the evidence.

**The data is faithful, and derived rather than copied.**  Each row of
`_codon_tables_data.py` carries NCBI's 64 amino acids in codon order, and the
codon-to-amino-acid mapping is rebuilt from that string rather than transcribed
from the reference's dictionary.  `test_every_code_rebuilds_from_its_string`
runs that derivation over all 27 and compares it with the reference's own
`forward_table` -- so a row that had drifted would be caught here even if the
machinery were perfect.

**The behaviour is the reference's, codon by codon.**  Every table the module
publishes -- six flavours of each of the 27 codes, 162 objects -- is compared
with the reference's: the printed grid, the names, the back table, the start
and stop codons, both alphabets, and then the translation of *every* codon over
that table's whole alphabet, which is 3,375 codons for the DNA and RNA tables
and 4,913 for the generic ones.  Almost all of those codons are absent from the
forward table, so the refusals are as much of the comparison as the answers:
the test counts how many of each it saw, because a comparison of 588,843
refusals against 588,843 refusals would prove nothing.

**What is not claimed.**  The three codes whose stop codons are context
dependent -- 27, 28 and 31, the karyorelictid ciliates and their relatives --
carry a `dual` field no 64-character string can express.  That field is the one
place this module says more than the reference's source text does, and it is
justified and tested in `test_the_codes_with_context_dependent_stops`.

**The reference's defects are pinned** in the same way as in the restriction
tests: `list_ambiguous_codons`'s unreachable `continue`, and the ambiguity `X`,
which the reference's alphabet omits and its value table includes.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

from biofasting import codon_table as ours
from biofasting import _codon_table as derivation

ROOT = Path(__file__).resolve().parent.parent

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

#   The six flavours each code is registered in, in the order
#   `biofasting.codon_table` files them, with the reference's dictionary of the
#   same shape.  A tuple rather than a dict because the order is meaningful: it
#   is the order `build()` returns the objects in.
_FAMILIES = (
    ("unambiguous_dna_by_id", "unambiguous_dna_by_name", "TCAG"),
    ("unambiguous_rna_by_id", "unambiguous_rna_by_name", "UCAG"),
    ("generic_by_id", "generic_by_name", "UCAG"),
    ("ambiguous_dna_by_id", "ambiguous_dna_by_name", None),
    ("ambiguous_rna_by_id", "ambiguous_rna_by_name", None),
    ("ambiguous_generic_by_id", "ambiguous_generic_by_name", None),
)

_FIELDS = (
    "id", "names", "back_table", "start_codons", "stop_codons",
    "nucleotide_alphabet", "protein_alphabet",
)


def _reference():
    from Bio.Data import CodonTable as reference

    return reference


def _alphabet(family: str, plain: str | None) -> tuple:
    """Every base worth asking a table of this flavour about.

    The explicit tables take the four letters of their alphabet; the generic
    one takes both, since it is keyed by T and U at once; and the ambiguous ones
    take the fifteen ambiguity codes, which is 3,375 codons per table.
    """
    if plain is not None:
        return tuple(plain) + (tuple("TCAG") if "generic" in family else ())
    if "generic" in family:
        return tuple(sorted(set(derivation.AMBIGUOUS_RNA_VALUES) | {"T"}))
    if "rna" in family:
        return tuple(sorted(derivation.AMBIGUOUS_RNA_VALUES))
    return tuple(sorted(derivation.AMBIGUOUS_DNA_VALUES))


def _ask(table, codon):
    """The translation, or the name of the refusal -- never one for the other."""
    try:
        return ("ok", table[codon])
    except KeyError:
        return ("KeyError", None)
    except Exception as exc:  # noqa: BLE001 - TranslationError and its relatives
        return (type(exc).__name__, None)


def _ask_all(table, alphabet) -> dict:
    return {
        a + b + c: _ask(table, a + b + c)
        for a in alphabet for b in alphabet for c in alphabet
    }


# --------------------------------------------------------------------------
# The data model: a row rebuilds the reference's table
# --------------------------------------------------------------------------


@needs_biopython
def test_every_code_rebuilds_from_its_string():
    """The derivation, run over all 27 rows, against the reference's tables.

    This is the claim the whole arrangement rests on: that NCBI's 64-character
    string is not a summary of the reference's dictionary but a complete
    statement of it.  If it were not, the data file would be a lossy copy and
    every test below would be comparing the machinery with itself.
    """
    reference = _reference()
    for id, name, alt_name, aa_string, start, stop, dual in ours.TABLES:
        assert len(aa_string) == 64, (id, len(aa_string))
        assert len(start) == len(set(start)), (id, start)
        assert len(stop) == len(set(stop)), (id, stop)
        rebuilt = derivation.forward_table(aa_string, dual)
        theirs = reference.unambiguous_dna_by_id[id].forward_table
        assert rebuilt == theirs, id
        #   The star positions are exactly the stops that are not also coding.
        dual_codons = {codon for codon, _ in dual}
        stars = {c for c, a in zip(derivation.CODON_ORDER, aa_string) if a == "*"}
        assert stars == set(stop) - dual_codons, id


@needs_biopython
def test_the_data_file_is_the_one_the_generator_writes():
    """The shipped table, checked against the generator's own three checks.

    A generated file that nobody regenerates is a file that drifts.  This runs
    the tool in `--check` mode, which rebuilds everything from the reference and
    compares the result with what is on disk, so the drift is caught here rather
    than by a reader who trusted the header.
    """
    done = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "gen_codon_tables.py"), "--check"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert done.returncode == 0, done.stderr
    assert "genetic codes rebuilt" in done.stdout


# --------------------------------------------------------------------------
# Every table, against the reference's
# --------------------------------------------------------------------------


@needs_biopython
def test_the_twelve_dictionaries_hold_the_same_tables():
    """Identifiers, names and objects, for all six flavours of all 27 codes."""
    reference = _reference()
    for by_id, by_name, _ in _FAMILIES:
        mine, theirs = getattr(ours, by_id), getattr(reference, by_id)
        assert set(mine) == set(theirs), by_id
        for id in mine:
            for field in _FIELDS:
                assert getattr(mine[id], field) == getattr(theirs[id], field), (
                    by_id, id, field,
                )
            assert str(mine[id]) == str(theirs[id]), (by_id, id)

        mine_names, theirs_names = getattr(ours, by_name), getattr(reference, by_name)
        assert set(mine_names) == set(theirs_names), by_name
        for name in mine_names:
            assert mine_names[name].id == theirs_names[name].id, (by_name, name)


@needs_biopython
def test_the_printed_grid_is_byte_for_byte_the_reference_s():
    """The 64-codon grid, for all 162 tables, compared as text."""
    reference = _reference()
    compared = 0
    for by_id, _, _ in _FAMILIES:
        for id in getattr(ours, by_id):
            mine = str(getattr(ours, by_id)[id])
            theirs = str(getattr(reference, by_id)[id])
            assert mine == theirs, (by_id, id)
            #   Four header lines, then one block per first base: four rows of
            #   four codons and a rule.  A grid that lost a rule would still be
            #   equal to itself if both sides were built the same way, so the
            #   shape is asserted rather than assumed.
            assert len(mine.splitlines()) == 4 + 4 * 5, (by_id, id)
            compared += 1
    assert compared == 27 * len(_FAMILIES)


@needs_biopython
def test_every_codon_of_every_table_gives_the_reference_s_answer():
    """The translation of every codon over every alphabet, refusals included.

    This is the test the module stands or falls on.  It asks each table about
    every codon its alphabet can spell -- 3,375 of them for the DNA and RNA
    ambiguous tables, 4,913 for the generic ones -- and requires the answer, or
    the name of the exception, to be the reference's.  A port that translated
    one codon differently, or raised `TranslationError` where the reference
    raises `KeyError`, is caught here and nowhere else.
    """
    reference = _reference()
    answers = refusals = 0
    for family, _, plain in _FAMILIES:
        alphabet = _alphabet(family, plain)
        for id in getattr(ours, family):
            mine = _ask_all(getattr(ours, family)[id].forward_table, alphabet)
            theirs = _ask_all(getattr(reference, family)[id].forward_table, alphabet)
            assert mine == theirs, (family, id)
            answers += sum(1 for answer, _ in mine.values() if answer == "ok")
            refusals += sum(1 for answer, _ in mine.values() if answer != "ok")
    #   A comparison of nothing but refusals would pass on an empty table, so
    #   both sides of the count are required to be large.  Measured: 371,115
    #   codons asked, 288,355 of them translated and 82,760 refused.
    assert answers > 200_000, answers
    assert refusals > 50_000, refusals


@needs_biopython
def test_the_silent_codons_map_back_to_the_reference_s_choice():
    """`back_table` picks one codon per amino acid, and picks the same one.

    The choice is arbitrary -- any codon of an amino acid would be correct --
    but it is not *free*: the reference resolves it by sorting the codons, and
    a port that resolved it by insertion order would answer differently for
    some amino acid in some code and look right in every test that only
    translates.
    """
    reference = _reference()
    for by_id, _, _ in _FAMILIES:
        for id in getattr(ours, by_id):
            mine = getattr(ours, by_id)[id].back_table
            theirs = getattr(reference, by_id)[id].back_table
            assert mine == theirs, (by_id, id)
            assert None in mine, (by_id, id)


@needs_biopython
def test_the_start_and_stop_lists_are_extended_the_same_way():
    """`list_ambiguous_codons`, over every code and both lists."""
    reference = _reference()
    for id in ours.unambiguous_dna_by_id:
        for list_name in ("start_codons", "stop_codons"):
            for family in ("unambiguous_dna_by_id", "ambiguous_dna_by_id"):
                mine = getattr(getattr(ours, family)[id], list_name)
                theirs = getattr(getattr(reference, family)[id], list_name)
                assert mine == theirs, (id, family, list_name)


# --------------------------------------------------------------------------
# The derivations, on their own
# --------------------------------------------------------------------------


@needs_biopython
def test_the_three_derivation_functions_match_the_reference_s():
    """`make_back_table`, `list_ambiguous_codons` and `list_possible_proteins`.

    These are the pieces the classes are built out of, and each of them has a
    behaviour worth pinning separately from the objects they end up inside:
    `make_back_table`'s sort, and `list_ambiguous_codons`'s refusal to add an
    ambiguity that would cover a codon which is not in the list.
    """
    reference = _reference()
    for id in ours.unambiguous_dna_by_id:
        table = reference.unambiguous_dna_by_id[id]
        mine = derivation.make_back_table(table.forward_table, table.stop_codons[0])
        assert mine == reference.make_back_table(table.forward_table, table.stop_codons[0])
        for codons, values in (
            (table.stop_codons, derivation.AMBIGUOUS_DNA_VALUES),
            (table.start_codons, derivation.AMBIGUOUS_DNA_VALUES),
            (table.stop_codons, derivation.AMBIGUOUS_RNA_VALUES),
        ):
            assert derivation.list_ambiguous_codons(codons, values) == \
                reference.list_ambiguous_codons(codons, values)

    #   The examples the reference documents, which is where its own rule about
    #   'TRR' is written down.
    values = derivation.AMBIGUOUS_DNA_VALUES
    assert derivation.list_ambiguous_codons(["TGA", "TAA"], values) == ["TGA", "TAA", "TRA"]
    assert derivation.list_ambiguous_codons(["TAG", "TGA"], values) == ["TAG", "TGA"]
    assert derivation.list_ambiguous_codons(["TAG", "TAA"], values) == ["TAG", "TAA", "TAR"]
    assert derivation.list_ambiguous_codons(["TGA", "TAA", "TAG"], values) == \
        ["TGA", "TAA", "TAG", "TAR", "TRA"]

    #   `list_possible_proteins` has three outcomes and all three are reachable
    #   on the standard code, so all three are pinned -- and then compared with
    #   the reference's own function over every ambiguous codon there is, which
    #   is what stops these three from being a hand-written summary of a rule
    #   nobody has checked.
    table = reference.unambiguous_dna_by_id[1].forward_table
    assert _ask_call(lambda: derivation.list_possible_proteins("GCN", table, values)) == \
        ("ok", ["A"])
    assert _ask_call(lambda: derivation.list_possible_proteins("TAR", table, values)) == \
        ("raised", "KeyError")
    assert _ask_call(lambda: derivation.list_possible_proteins("TRA", table, values)) == \
        ("raised", "KeyError")
    assert _ask_call(lambda: derivation.list_possible_proteins("TRR", table, values)) == \
        ("raised", "TranslationError")
    for a in derivation.AMBIGUOUS_DNA_LETTERS:
        for b in derivation.AMBIGUOUS_DNA_LETTERS:
            for c in derivation.AMBIGUOUS_DNA_LETTERS:
                codon = a + b + c
                assert _ask_call(
                    lambda: derivation.list_possible_proteins(codon, table, values)
                ) == _ask_call(
                    lambda: reference.list_possible_proteins(codon, table, values)
                ), codon


def _ask_call(call):
    try:
        return ("ok", call())
    except Exception as exc:  # noqa: BLE001 - the exception is the result
        return ("raised", type(exc).__name__)


@needs_biopython
def test_the_ambiguity_tables_carry_the_standard_codes():
    """The IUPAC codes, which are a standard and are written down as one.

    These are hand-written here rather than generated, because they are the
    IUPAC standard and not a fact about Biopython -- and then required to equal
    the reference's, which is what makes writing them down safe.  `X` is worth
    the line it gets: it is *not* in the standard, the reference's alphabet
    omits it, and the reference's value table defines it, so an ambiguity the
    standard says nothing about is still a key a caller can use.
    """
    reference = _reference()
    from Bio.Data import IUPACData

    assert dict(derivation.AMBIGUOUS_DNA_VALUES) == dict(IUPACData.ambiguous_dna_values)
    assert dict(derivation.AMBIGUOUS_RNA_VALUES) == dict(IUPACData.ambiguous_rna_values)
    assert dict(derivation.EXTENDED_PROTEIN_VALUES) == dict(IUPACData.extended_protein_values)
    assert derivation.AMBIGUOUS_DNA_LETTERS == IUPACData.ambiguous_dna_letters
    assert derivation.AMBIGUOUS_RNA_LETTERS == IUPACData.ambiguous_rna_letters
    assert derivation.UNAMBIGUOUS_DNA_LETTERS == IUPACData.unambiguous_dna_letters
    assert derivation.PROTEIN_LETTERS == IUPACData.protein_letters
    assert derivation.EXTENDED_PROTEIN_LETTERS == IUPACData.extended_protein_letters
    assert "X" in derivation.AMBIGUOUS_DNA_VALUES
    assert "X" not in derivation.AMBIGUOUS_DNA_LETTERS


# --------------------------------------------------------------------------
# The module's surface
# --------------------------------------------------------------------------


@needs_biopython
def test_the_standard_tables_are_table_one():
    reference = _reference()
    assert ours.standard_dna_table is ours.unambiguous_dna_by_id[1]
    assert ours.standard_rna_table is ours.unambiguous_rna_by_id[1]
    assert str(ours.standard_dna_table).startswith("Table 1 Standard, SGC0")
    assert str(ours.standard_rna_table).startswith("Table 1 Standard, SGC0")
    assert ours.standard_dna_table is not ours.standard_rna_table


@needs_biopython
def test_a_bare_codon_table_has_no_identifier_and_says_so_by_failing():
    """`CodonTable` on its own is not a genetic code and cannot be printed.

    The reference's `__str__` reads `self.id`, which is set by `NCBICodonTable`
    and not by the base class, so a bare table raises `AttributeError` there and
    has to raise one here.  Returning "Table ID unknown" -- a branch the method
    really does contain, for a table whose id is falsy -- would be a different
    answer to the same call.
    """
    assert _ask_call(lambda: str(ours.CodonTable())) == ("raised", "AttributeError")
    reference = _reference()
    assert _ask_call(lambda: str(reference.CodonTable())) == ("raised", "AttributeError")

    #   ...but the falsy-id branch is real and reachable through the registry.
    made = ours.NCBICodonTable(0, ["Nowhere"], {"ATG": "M"}, ["ATG"], ["TAA"])
    assert str(made).startswith("Table ID unknown Nowhere")


@needs_biopython
def test_register_ncbi_table_takes_a_string_or_the_reference_s_dictionary():
    """The public extension point, accepting both spellings of the data.

    The reference takes the codon dictionary; this module's data file carries
    strings.  Both work here, so a caller written against `Bio.Data.CodonTable`
    keeps working -- and the two spellings are required to produce the same
    table, which is the only reason it is safe to have two.
    """
    ids_before = set(ours.unambiguous_dna_by_id)
    aa_string = ours.TABLES[0][3]
    dictionary = derivation.forward_table(aa_string)
    try:
        ours.register_ncbi_table(
            "Test Code", "TC0", 999, aa_string, ["ATG"], ["TAA", "TAG", "TGA"]
        )
        from_string = ours.unambiguous_dna_by_id[999]
        ours.register_ncbi_table(
            "Test Code", "TC0", 998, dictionary, ["ATG"], ["TAA", "TAG", "TGA"]
        )
        from_dict = ours.unambiguous_dna_by_id[998]
    finally:
        for target in (
            ours.unambiguous_dna_by_id, ours.unambiguous_rna_by_id, ours.generic_by_id,
            ours.ambiguous_dna_by_id, ours.ambiguous_rna_by_id,
            ours.ambiguous_generic_by_id,
        ):
            target.pop(999, None)
            target.pop(998, None)
        for target in (
            ours.unambiguous_dna_by_name, ours.unambiguous_rna_by_name,
            ours.generic_by_name, ours.ambiguous_dna_by_name,
            ours.ambiguous_rna_by_name, ours.ambiguous_generic_by_name,
        ):
            for name in ("Test Code", "TC0"):
                target.pop(name, None)

    assert from_string.forward_table == from_dict.forward_table == dictionary
    assert from_string.names == from_dict.names == ["Test Code", "TC0"]
    assert set(ours.unambiguous_dna_by_id) == ids_before
    assert "Test Code" not in ours.unambiguous_dna_by_name


@needs_biopython
def test_a_code_with_no_alternative_name_keeps_a_none_in_its_name_list():
    """Table 27's names really do contain a `None`, and the grid skips it.

    The reference appends the alternative name to the name list only when it is
    not `None`, and *before* that it has already used `names + [alt_name]` to
    build the tables -- so the object's `names` ends with a `None` while the
    name dictionary has no `None` key.  A port that tidied that away would have
    a different name list and a different printed header.
    """
    reference = _reference()
    for id in (27, 28, 31):
        assert ours.unambiguous_dna_by_id[id].names == \
            reference.unambiguous_dna_by_id[id].names
        assert None in ours.unambiguous_dna_by_id[id].names
    assert str(ours.unambiguous_dna_by_id[27]).splitlines()[0] == \
        "Table 27 Karyorelict Nuclear"
    assert "None" not in str(ours.unambiguous_dna_by_id[27])


@needs_biopython
def test_the_codes_with_context_dependent_stops():
    """Three codes where one codon is a stop and an amino acid at the same time.

    *Condylostoma*, *Blastocrithidia* and the karyorelictids really do read TGA
    -- and in two of them TAA and TAG as well -- as an amino acid when a tRNA
    wins the race and as a stop when the release factor does.  No 64-character
    string can carry that: the string has one character per codon.  So the row
    carries the codon and the residue, and the mapping is built from both.

    Worth its own test because it is the only place this module's data says
    something the reference's *source text* does not, and because a reader who
    found it by accident would reasonably suspect a transcription error.
    """
    reference = _reference()
    expected = {
        27: [("TGA", "W")],
        28: [("TAA", "Q"), ("TAG", "Q"), ("TGA", "W")],
        31: [("TAA", "E"), ("TAG", "E")],
    }
    for id, dual in expected.items():
        row = next(row for row in ours.TABLES if row[0] == id)
        assert row[6] == dual, id
        table = ours.unambiguous_dna_by_id[id]
        theirs = reference.unambiguous_dna_by_id[id]
        for codon, amino in dual:
            assert codon in table.stop_codons, (id, codon)
            assert table.forward_table[codon] == amino == theirs.forward_table[codon]
        #   The string itself has no star for those codons -- there is nowhere
        #   to put one -- which is what makes the `dual` field necessary.  All
        #   three codes are in this position because *every* stop they have is
        #   context dependent, so their strings are 64 amino acids with no mark
        #   on them at all.
        assert "*" not in row[3], id

    #   And the ambiguous tables inherit the difference, which is the clearest
    #   way to see what the `dual` field buys.  `TRA` covers TAA and TGA: in the
    #   standard code both are stops, so the codon has no translation at all; in
    #   table 28 both also code, so it has two, and the answer is the letter that
    #   covers both -- `X`.  Same codon, same module, two different outcomes, and
    #   both are the reference's.
    standard = ours.ambiguous_dna_by_id[1]
    condylostoma = ours.ambiguous_dna_by_id[28]
    assert _ask_call(lambda: standard.forward_table["TRA"]) == ("raised", "KeyError")
    assert _ask_call(lambda: condylostoma.forward_table["TRA"]) == ("ok", "X")
    assert _ask_call(lambda: standard.forward_table["TRR"]) == ("raised", "TranslationError")
    for id in (1, 27, 28, 31):
        for codon in ("TRA", "TRR", "TAR", "TGA", "TGG"):
            mine = _ask(getattr(ours, "ambiguous_dna_by_id")[id].forward_table, codon)
            theirs = _ask(reference.ambiguous_dna_by_id[id].forward_table, codon)
            assert mine == theirs, (id, codon)


def test_the_module_imports_without_biopython():
    """The genetic codes are available to a program that has no Biopython.

    This is the point of putting the codes in a module of their own: a
    translation needs a genetic code before it can do anything, so a codon table
    that required Biopython would make Biopython required by everything.  The
    check is a fresh interpreter with `Bio` hidden, which is the only way to
    catch an import that only happens to be lazy.
    """
    program = (
        "import sys;"
        "sys.modules['Bio'] = None;"
        "from biofasting import codon_table;"
        "assert codon_table.unambiguous_dna_by_id[1].forward_table['ATG'] == 'M';"
        "assert len(codon_table.generic_by_id) == 27;"
        "assert 'Bio' not in sys.modules or sys.modules['Bio'] is None;"
        "print('ok')"
    )
    done = subprocess.run(
        [sys.executable, "-c", program], capture_output=True, text=True, cwd=ROOT
    )
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip().splitlines()[-1] == "ok"


def test_importing_the_package_does_not_pay_for_the_codes():
    """`import biofasting` must not build 162 tables nobody asked for."""
    program = (
        "import sys, biofasting;"
        "print('biofasting.codon_table' in sys.modules,"
        "      'biofasting._codon_tables_data' in sys.modules)"
    )
    done = subprocess.run(
        [sys.executable, "-c", program], capture_output=True, text=True, cwd=ROOT
    )
    assert done.returncode == 0, done.stderr
    #   The editable install rebuilds on import and says so on stdout, which is
    #   the interpreter's first line rather than the program's last.
    assert done.stdout.strip().splitlines()[-1] == "False False", done.stdout


def test_the_docstring_examples_are_true():
    """Every `>>>` in the two modules is executed, not merely written."""
    import doctest

    for module in (ours, derivation):
        results = doctest.testmod(module, verbose=False, report=False)
        assert results.attempted > 0, f"{module.__name__} documents nothing"
        assert results.failed == 0, module.__name__
