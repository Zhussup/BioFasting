"""Phase 2 step 2.2c: the substitution matrices, as data and as a reader.

`Bio.Align.substitution_matrices` is a 508-line module and a directory of thirty
plain-text scoring tables, 84,540 bytes of them.  `biofasting.substitution_matrices`
is the same thirty tables as 8,388 numbers in a generated Python module, plus the
`Array` class that indexes them by letter.  Saying that is cheap; this file is the
evidence.

**The data is faithful, and derived rather than copied.**  Each row of
`_substitution_matrices_data.py` carries a matrix's numbers as the file's own
tokens, stored as the lower triangle because the matrix is symmetric.
`test_every_matrix_is_stored_as_the_triangle_of_a_symmetric_one` checks the
symmetry claim itself -- all thirty, over every entry -- and the module-level
comparison below would fail loudly on a matrix that was not, since a triangle
that is mirrored wrongly is a full square of wrong numbers.

**The behaviour is the reference's, entry for entry.**  All thirty matrices are
loaded through our module and through Biopython's and compared on the alphabet,
the dtype, the header, the shape, all 16,140 entries, the string `format()`
prints and the `repr`.  The parser is run over every file the reference ships,
from a path and from a handle, and then over a set of malformed tables, where
what is compared is the exception *type* -- including two cases where the
reference's own failure is an `UnboundLocalError` and an `AssertionError`, which
are pinned rather than repaired.

**The class is the reference's, call for call.**  82 operations -- construction
from a dictionary, from an alphabet, from a numpy array; indexing by letter, by
tuple, by slice; `get`, `keys`, `values`, `items`, `update`, `select`; arithmetic
and reductions; pickling; and the alphabet's own refusals -- produce byte-identical
results, except for the one difference the module docstring states: an unknown
name is not looked for in a directory, because there is no directory.

**What is not claimed.**  The reference implements `Array` twice, and the half
that matters is a C type, `_arraycore.Array`, whose `alphabet` property validates
the alphabet against the array's shape.  That validation is reimplemented here
in Python from the C source, and `test_the_alphabet_is_set_once_and_has_to_fit_the_array`
and `test_a_slice_whose_shape_the_alphabet_cannot_describe_is_refused` are what
say the reimplementation is complete for the cases a caller can reach.
"""

import importlib.util
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from biofasting import _substitution_matrices as derivation
from biofasting import substitution_matrices as ours

ROOT = Path(__file__).resolve().parent.parent

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)


def reference_module():
    """The reference, and the directory its thirty files are in."""
    from Bio.Align import substitution_matrices as reference

    return reference, Path(reference.__file__).resolve().parent / "data"


def reference_names():
    _, directory = reference_module()
    return sorted(name for name in directory.iterdir() if name.name != "README.txt")


@needs_biopython
def test_the_data_file_is_the_one_the_generator_writes():
    """The generator's three checks, run as a test rather than trusted.

    Running the generator rather than re-deriving its answer here is the point:
    a data file that has been hand-edited, or a derivation that has drifted from
    it, is caught by the same code that wrote it -- and the third check loads
    every name through the package, so the whole path from row to caller is
    walked before this passes.
    """
    done = subprocess.run(
        [sys.executable, "tools/gen_substitution_matrices.py", "--check"],
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    assert done.returncode == 0, done.stderr
    #   The editable install rebuilds the C extension when a fresh interpreter
    #   imports the package, and says so on stdout, so the generator's own lines
    #   are found by their text rather than by counting from the end.
    lines = done.stdout.splitlines()
    assert any(line.startswith("derivation: 30 matrices") for line in lines)
    assert any("is current" in line for line in lines)
    assert any(line.startswith("module: all 30 names load") for line in lines)


@needs_biopython
def test_load_lists_the_reference_s_own_thirty_names():
    """`load()` and `os.listdir` must agree, or a name that ships is unreachable."""
    reference, _ = reference_module()
    assert ours.load() == reference.load()
    assert len(ours.load()) == 30


@needs_biopython
def test_every_matrix_is_the_reference_s():
    """All thirty matrices, every field and every one of their 16,140 numbers."""
    reference, _ = reference_module()
    entries = 0
    for name in reference.load():
        mine = ours.load(name)
        theirs = reference.load(name)
        assert mine.alphabet == theirs.alphabet, name
        assert mine.dtype == theirs.dtype, name
        assert list(mine.header) == list(theirs.header), name
        assert mine.shape == theirs.shape, name
        assert (np.asarray(mine) == np.asarray(theirs)).all(), name
        assert str(mine) == str(theirs), name
        assert repr(mine) == repr(theirs), name
        entries += mine.size
    assert entries == 16140, entries


@needs_biopython
def test_every_matrix_is_stored_as_the_triangle_of_a_symmetric_one():
    """The symmetry the data file relies on, checked rather than asserted.

    Every one of the thirty is symmetric, so storing half of each is lossless --
    but "every one is symmetric" is a claim about the data, not a property of the
    format, and the row carries a flag saying which of the two it is.  This test
    is what makes the flag meaningful: it re-checks the symmetry against the
    reference and counts that the file really is half the size of the matrices.
    """
    reference, _ = reference_module()
    full = stored = 0
    for name in reference.load():
        theirs = reference.load(name)
        n = theirs.shape[0]
        assert (np.asarray(theirs) == np.asarray(theirs).T).all(), name
        full += n * n
        stored += n * (n + 1) // 2
    from biofasting._substitution_matrices_data import TABLES

    assert all(row[3] for row in TABLES)
    assert sum(len(row[4].split()) for row in TABLES) == stored
    assert (full, stored) == (16140, 8388)
    #   And the flag is not decoration: a row claiming symmetry is built by
    #   mirroring, so a wrong claim would put the numbers in the wrong cells.
    for name, header, alphabet, symmetric, values in TABLES:
        matrix = derivation.load_row(header, alphabet, symmetric, values)
        assert (np.asarray(matrix) == np.asarray(matrix).T).all(), name
        assert matrix.shape == (len(alphabet),) * 2


@needs_biopython
def test_the_parser_reads_every_file_the_reference_ships():
    """`read` on all thirty files, from a path and from an open handle."""
    reference, directory = reference_module()
    for path in reference_names():
        theirs = reference.read(str(path))
        from_path = ours.read(str(path))
        assert from_path.alphabet == theirs.alphabet, path.name
        assert list(from_path.header) == list(theirs.header), path.name
        assert (np.asarray(from_path) == np.asarray(theirs)).all(), path.name
        with open(path) as fp:
            from_handle = ours.read(fp)
        assert (np.asarray(from_handle) == np.asarray(theirs)).all(), path.name


#   Tables that are not tables.  The two that fail with an `UnboundLocalError`
#   and an `AssertionError` are the reference's own behaviour, not ours: an
#   empty file leaves the loop variable `i` unbound and a one-row file whose
#   second row has one field is read as a square, where the row labels are then
#   checked against the column labels with an `assert`.
_MALFORMED = {
    "nothing but comments": "# hi\n# there\n",
    "a blank line inside a square": "  A B\nA 1 2\n\nB 3 4\n",
    "a row labelled with a letter that is not a column": "  A B\nA 1 2\nZ 3 4\n",
    "a row with the wrong number of fields": "  A B\nA 1\nB 3 4\n",
    "a single row": "  A B\n",
    "not a table at all": "A 1\nB\n",
    "an empty file": "",
}


@needs_biopython
@pytest.mark.parametrize("label", sorted(_MALFORMED))
def test_the_parser_refuses_what_the_reference_refuses(label):
    """Both parsers raise, and they raise the same kind of thing.

    A parser that returned an empty matrix where the reference raised would be
    worse than a slower one, so the failure is compared as carefully as the
    success: the exception name has to match, and the fact that *something* was
    raised is asserted separately, so a case that stopped raising at all cannot
    pass by agreeing about nothing.
    """
    reference, _ = reference_module()
    text = _MALFORMED[label]
    with pytest.raises(Exception) as theirs:
        reference.read(io.StringIO(text))
    with pytest.raises(Exception) as mine:
        ours.read(io.StringIO(text))
    assert type(mine.value).__name__ == type(theirs.value).__name__


@needs_biopython
def test_the_parser_reads_a_one_dimensional_table_too():
    """The other shape a matrix file comes in, and the alphabet it implies.

    A file whose first two rows both have two fields is read as `letter value`
    pairs, whatever its header says.  Which letters those are decides the
    alphabet: one character each makes a string, a longer name makes a tuple --
    the reference's own rule, and the reason `SCHNEIDER`'s alphabet is a tuple.
    """
    reference, _ = reference_module()
    for text, expected in (
        ("A 1\nB 2\n", "AB"),
        ("AA 1\nBB 2\n", ("AA", "BB")),
        ("# note\nA 1\nB 2\n", "AB"),
    ):
        mine = ours.read(io.StringIO(text))
        theirs = reference.read(io.StringIO(text))
        assert mine.alphabet == theirs.alphabet == expected
        assert mine.shape == theirs.shape == (2,)
        assert (np.asarray(mine) == np.asarray(theirs)).all()
    assert ours.read(io.StringIO("# a\n#b\nA 1\nB 2\n")).header == ["a", "b"]


@needs_biopython
def test_a_matrix_is_indexed_by_letter_and_refuses_a_letter_it_does_not_have():
    """Letters in, numbers out -- and an `IndexError` for a letter not in it."""
    reference, _ = reference_module()
    mine = ours.load("BLOSUM62")
    theirs = reference.load("BLOSUM62")
    assert mine["A", "R"] == theirs["A", "R"] == -1.0
    assert mine[("W", "W")] == theirs[("W", "W")] == 11.0
    assert mine["A"].shape == theirs["A"].shape == (24,)
    assert mine[0, 0] == theirs[0, 0]
    for key in ("ZZ", ("ZZ", "A")):
        with pytest.raises(IndexError) as theirs_error:
            theirs[key]
        with pytest.raises(IndexError) as mine_error:
            mine[key]
        assert str(mine_error.value) == str(theirs_error.value)
    for key in ("ZZ", 0, ("A", "R")):
        assert (key in mine) == (key in theirs)
    assert mine.get("ZZ", "fallback") == theirs.get("ZZ", "fallback") == "fallback"


@needs_biopython
def test_the_alphabet_is_set_once_and_has_to_fit_the_array():
    """The three refusals the C `alphabet` property makes, message for message.

    This is the half of `Array` that is not written in the reference's Python:
    the setter lives in `_arraycore.c`, which ships with Biopython and was read
    to write it here.  A repeated letter is refused because a string alphabet
    is cached as a character-to-index mapping, and two indices for one character
    would be a lie about what `matrix['A']` means.
    """
    reference, _ = reference_module()
    for make in (
        lambda m: _assign(m.Array("AB"), "X"),
        lambda m: _assign(m.Array("AAB"), None),
        lambda m: _assign(m.Array("AB", dims=2), "ABC"),
        lambda m: _assign(m.Array("AB"), 5),
        lambda m: _assign(m.Array("AB"), {"A": 1}),
    ):
        with pytest.raises(Exception) as theirs:
            make(reference)
        with pytest.raises(Exception) as mine:
            make(ours)
        assert type(mine.value).__name__ == type(theirs.value).__name__
        assert str(mine.value) == str(theirs.value)
    #   A tuple alphabet may hold things that are not letters at all: there is
    #   no mapping to build for a tuple, so nothing to check.
    assert ours.Array((1, 2)).alphabet == reference.Array((1, 2)).alphabet
    #   And an array that was never given one answers `None` rather than failing.
    assert np.zeros(4).view(ours.Array).alphabet is None
    assert np.zeros(4).view(reference.Array).alphabet is None


def _assign(matrix, alphabet):
    matrix.alphabet = alphabet
    return matrix.alphabet


@needs_biopython
def test_a_slice_whose_shape_the_alphabet_cannot_describe_is_refused():
    """Slicing a matrix to something the alphabet no longer fits raises.

    A view keeps the alphabet, and the alphabet is checked against the view's
    shape as it is made -- so `matrix[0, 0:2]` fails where a plain numpy array
    would have quietly handed back two numbers under twenty-four letters.  The
    check in `__getitem__` itself is reachable only for a shape the alphabet
    *does* fit and the matrix does not, which is a single row or column; both
    are here.
    """
    reference, _ = reference_module()
    mine = ours.load("BLOSUM62")
    theirs = reference.load("BLOSUM62")
    for key in ((0, slice(0, 2)), (slice(0, 25), slice(0, 2))):
        with pytest.raises(Exception) as theirs_error:
            theirs[key]
        with pytest.raises(Exception) as mine_error:
            mine[key]
        assert type(mine_error.value).__name__ == type(theirs_error.value).__name__
        assert str(mine_error.value) == str(theirs_error.value)
    #   A slice the alphabet does fit, and the matrix does not.
    with pytest.raises(IndexError) as theirs_error:
        theirs[0:24, 0:1]
    with pytest.raises(IndexError) as mine_error:
        mine[0:24, 0:1]
    assert "truncated" in str(mine_error.value) == str(theirs_error.value)


@needs_biopython
def test_the_dictionary_constructor_has_the_reference_s_defect():
    """`Array(data={"A": 1})` is an `UnboundLocalError`, here as there.

    The one-dimensional dictionary path assigns from `letter` where it means
    `key`, and `letter` is a name only the two-dimensional path would have left
    behind -- so a dictionary of plain strings, which is the only thing a
    one-dimensional dictionary can be, always fails.  Kept rather than repaired,
    because the two-dimensional path next to it works and is the one this
    library uses; pinned by name so that a Biopython that fixes it turns this
    red instead of letting the two drift apart.
    """
    reference, _ = reference_module()
    with pytest.raises(UnboundLocalError) as theirs:
        reference.Array(data={"A": 1, "B": 2})
    with pytest.raises(UnboundLocalError) as mine:
        ours.Array(data={"A": 1, "B": 2})
    assert str(mine.value) == str(theirs.value)
    #   The two-dimensional path, which does work, and the refusals around it.
    data = {("A", "A"): 1, ("A", "B"): 2, ("B", "B"): 3}
    for module in (reference, ours):
        matrix = module.Array(data=data)
        assert matrix.alphabet == "AB"
        assert matrix.shape == (2, 2)
        #   A key that is absent is a zero, and the mirror of a key that is
        #   present is not filled in: the matrix is not symmetrised for you.
        assert (np.asarray(matrix) == np.array([[1.0, 2.0], [0.0, 3.0]])).all()
    for make, expected in (
        (lambda m: m.Array(data={}), ValueError),
        (lambda m: m.Array(data={("A",): 1}), ValueError),
        (lambda m: m.Array(data={("A", "B", "C"): 1}), ValueError),
        (lambda m: m.Array(data={"A": 1, ("B", "C"): 2}), ValueError),
        (lambda m: m.Array("AB", data=data), ValueError),
        (lambda m: m.Array(dims=2, data=data), ValueError),
    ):
        with pytest.raises(expected) as theirs:
            make(reference)
        with pytest.raises(expected) as mine:
            make(ours)
        assert str(mine.value) == str(theirs.value)


@needs_biopython
def test_update_items_keys_values_and_select_agree_with_the_reference():
    """The dictionary-like surface, including where it is inconsistent.

    `items()` walks a square by rows while `keys()` and `values()` walk it by
    columns, so `list(matrix.items())` and `list(zip(matrix.keys(),
    matrix.values()))` are different lists.  Both are the reference's, and both
    are compared here rather than one of them being quietly made to match the
    other.
    """
    reference, _ = reference_module()
    mine = ours.load("BLOSUM62")
    theirs = reference.load("BLOSUM62")
    assert mine.keys() == theirs.keys()
    assert [float(v) for v in mine.values()] == [float(v) for v in theirs.values()]
    assert [
        (key, float(value)) for key, value in mine.items()
    ] == [(key, float(value)) for key, value in theirs.items()]
    assert list(mine.items())[:2] == list(theirs.items())[:2]
    #   The inconsistency itself, so that a later tidy-up of either side is
    #   caught here rather than changing what a caller sees in silence.
    assert list(zip(mine.keys(), mine.values())) != list(mine.items())

    subset = mine.select("AW")
    assert subset.alphabet == theirs.select("AW").alphabet == "AW"
    assert (np.asarray(subset) == np.asarray(theirs.select("AW"))).all()
    #   Every letter of `AXZ` is a letter of BLOSUM62, so nothing is passed over
    #   and the result is the three-by-three corner those letters select.  `J`
    #   is not a letter it has, and *that* is the case the silent `continue` is
    #   for: it is dropped from the alphabet asked for and the shape follows the
    #   request, not the finding.
    for module in (reference, ours):
        found = module.load("BLOSUM62").select("AXZ")
        assert found.alphabet == "AXZ"
        assert found.shape == (3, 3)
        assert float(found["X", "Z"]) == float(theirs["X", "Z"])
        missing = module.load("BLOSUM62").select("AJ")
        assert missing.alphabet == "AJ"
        assert missing.shape == (2, 2)
        assert float(missing["A", "A"]) == float(theirs["A", "A"])
        assert float(missing["J", "J"]) == 0.0
    #   A repeated letter in the alphabet asked for is refused by the
    #   constructor, not by `select`.
    for module, matrix in ((reference, theirs), (ours, mine)):
        with pytest.raises(ValueError) as error:
            matrix.select("AZZ")
        assert "more than once" in str(error.value)

    for module in (reference, ours):
        vector = module.Array("AB")
        vector.update({"A": 1})
        vector.update([("B", 2)])
        assert float(vector["A"]) == 1.0
        assert float(vector["B"]) == 2.0
        square = module.Array("AB", dims=2)
        square.update({("A", "B"): 7})
        assert float(square["A", "B"]) == 7.0
        with pytest.raises(ValueError):
            vector.update(("A", 4))


@needs_biopython
def test_arithmetic_keeps_the_alphabet_and_refuses_two_of_them():
    """`Array + Array` is an `Array` again; two different alphabets are refused."""
    reference, _ = reference_module()
    for module in (reference, ours):
        a = module.Array("AB")
        b = module.Array("AB")
        total = a + b
        assert type(total).__name__ == "Array"
        assert total.alphabet == "AB"
        assert float(a.sum()) == 0.0
        with pytest.raises(ValueError) as error:
            a + module.Array("CD")
        assert str(error.value) == "alphabets are inconsistent"
    assert (
        ours.Array("AB")[0] == reference.Array("AB")[0] == 0.0
    ), "a single entry is a Python float on both sides"


@needs_biopython
def test_pickle_round_trips_a_square_and_not_a_vector():
    """Two-dimensional matrices pickle; one-dimensional ones do not.

    `__setstate__` assigns with `self[:, :]`, which is a two-dimensional index,
    so a vector unpickles into an `IndexError`.  Kept, because the reference's
    `__reduce__` is what a caller's stored matrix went through and a port that
    silently fixed it would accept pickles the reference cannot read.
    """
    import pickle

    reference, _ = reference_module()
    for module in (reference, ours):
        matrix = module.load("HOXD70")
        restored = pickle.loads(pickle.dumps(matrix))
        assert restored.alphabet == matrix.alphabet
        assert restored.dtype == matrix.dtype
        assert (np.asarray(restored) == np.asarray(matrix)).all()
        with pytest.raises(IndexError):
            pickle.loads(pickle.dumps(module.Array("AB")))


@needs_biopython
def test_format_prints_the_reference_s_table():
    """`format()`, `str()` and the `__format__` behind an f-string."""
    reference, _ = reference_module()
    for module in (reference, ours):
        assert module.Array("AB", dims=2).format() == "    A   B\nA 0.0 0.0\nB 0.0 0.0\n"
        assert module.Array("AB", dims=2).format("%i") == "  A B\nA 0 0\nB 0 0\n"
        assert module.Array("AB", dtype=int).format() == "A 0\nB 0\n"
        assert str(module.Array("AB")) == module.Array("AB").format()
        vector = module.Array("AB")
        vector[0], vector[1] = 1, 22
        #   The format spec goes to `%` as it stands, so an f-string works when
        #   it carries a `%` and raises `TypeError` when it does not -- which is
        #   the same `__format__` a caller gets from the reference.
        assert f"{vector:%.2f}" == "A  1.00\nB 22.00\n"
        with pytest.raises(TypeError):
            f"{vector:.2f}"
        #   An array of a rank the formatter has no layout for.  One cannot be
        #   built -- the alphabet refuses a rank of three -- but a plain numpy
        #   array can be viewed as one, and then `format` has nothing to print.
        with pytest.raises(RuntimeError) as error:
            np.zeros((2, 2, 2)).view(module.Array).format()
        assert "unexpected rank 3" in str(error.value)
    for name in ours.load():
        assert str(ours.load(name)) == str(reference.load(name)), name


def test_the_module_imports_without_biopython():
    """A scoring matrix is available to a program that has no Biopython.

    An aligner needs a substitution matrix before it can score anything, so a
    module that required Biopython would make Biopython required by everything
    above it.  The check is a fresh interpreter with `Bio` hidden, which is the
    only way to catch an import that merely happens to be lazy.
    """
    program = (
        "import sys;"
        "sys.modules['Bio'] = None;"
        "from biofasting import substitution_matrices as sm;"
        "assert sm.load('BLOSUM62')['A', 'R'] == -1.0;"
        "assert len(sm.load()) == 30;"
        "assert 'Bio' not in sys.modules or sys.modules['Bio'] is None;"
        "print('ok')"
    )
    done = subprocess.run(
        [sys.executable, "-c", program], capture_output=True, text=True, cwd=ROOT
    )
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip().splitlines()[-1] == "ok"


def test_importing_the_package_does_not_pay_for_the_matrices():
    """`import biofasting` must not load 8,388 numbers nobody asked for."""
    program = (
        "import sys, biofasting;"
        "print('biofasting.substitution_matrices' in sys.modules,"
        "      'biofasting._substitution_matrices_data' in sys.modules)"
    )
    done = subprocess.run(
        [sys.executable, "-c", program], capture_output=True, text=True, cwd=ROOT
    )
    assert done.returncode == 0, done.stderr
    #   The editable install rebuilds on import and says so on stdout, which is
    #   the interpreter's first line rather than the program's last.
    assert done.stdout.strip().splitlines()[-1] == "False False", done.stdout


@needs_biopython
def test_an_unknown_name_is_a_looked_for_file_and_says_so():
    """The one place this module is deliberately not the reference.

    The reference resolves a name inside its own `data` directory, so an unknown
    name fails on `open()` naming the path it looked in.  There is no such
    directory here -- the tables are numbers in a Python module -- so an unknown
    *path* still fails on `open()`, naming the path as given, while an unknown
    bare *name* fails with the thirty names that do exist.  Both are asserted,
    so that neither turns into a silent empty matrix.
    """
    reference, _ = reference_module()
    with pytest.raises(FileNotFoundError):
        reference.load("NOSUCH")
    with pytest.raises(FileNotFoundError) as error:
        ours.load("NOSUCH")
    assert "no substitution matrix named 'NOSUCH'" in str(error.value)
    #   A path that exists is read, the way `os.path.join` lets an absolute path
    #   through in the reference.
    path = reference_names()[0]
    assert (np.asarray(ours.load(str(path))) == np.asarray(reference.read(str(path)))).all()
    assert (np.asarray(ours.load(path)) == np.asarray(reference.load(path.name))).all()


def test_the_docstring_examples_are_true():
    """Every `>>>` in the public module is executed, not merely written."""
    import doctest

    results = doctest.testmod(ours, verbose=False, report=False)
    assert results.attempted > 0, "the module documents nothing"
    assert results.failed == 0


def _out_argument(module):
    """`out=` hands back the buffer -- as a plain array, without the alphabet."""
    a, b = module.Array("AB"), module.Array("AB")
    np.add(a, b, out=b)
    return type(b).__name__, b.alphabet


#   The corners, one call each, run against both modules and compared.  Most of
#   what is here is already asserted above with its own message; this table is
#   for the ones that are too small to deserve a test of their own but too easy
#   to get wrong -- a buffer that loses its alphabet, a `divmod` that returns
#   `np.float64` where an `Array` was expected, a `matmul`-style reduction that
#   never reaches `__array_ufunc__` at all.
_CORNERS = {
    "alphabet of a bare array": lambda m: m.Array().alphabet,
    "shape of a bare array": lambda m: m.Array().shape,
    "dtype of a bare array": lambda m: str(m.Array().dtype),
    "a tuple alphabet of non-letters": lambda m: m.Array((1, 2)).alphabet,
    "dims of three": lambda m: m.Array("AB", dims=3),
    "dims of zero": lambda m: m.Array("AB", dims=0),
    "bytes as an alphabet": lambda m: m.Array(b"AB"),
    "a list as an alphabet": lambda m: m.Array(["A"]),
    "one-dimensional data without an alphabet": lambda m: m.Array(data=np.zeros(2)),
    "a numpy array where an alphabet belongs": lambda m: m.Array(np.zeros(2)),
    "data whose shape does not match the alphabet": lambda m: m.Array(
        "AB", dims=1, data=np.zeros(3)
    ),
    "three-dimensional data": lambda m: m.Array(data=np.zeros((2, 2, 2))),
    "a non-square matrix": lambda m: m.Array(data=np.zeros((2, 3))),
    "the string of an alphabet-less array": lambda m: str(np.zeros(2).view(m.Array)),
    "a transposed matrix": lambda m: (type(m.load("BLOSUM62").T).__name__, m.load("BLOSUM62").T.alphabet[:3]),
    "a divmod": lambda m: [type(x).__name__ for x in np.divmod(m.Array("AB"), 2)],
    "a reduction": lambda m: float(np.add.reduce(m.Array("AB"))),
    "an in-place add": lambda m: np.add.at(m.Array("AB"), [0], 1),
    "an out= buffer": _out_argument,
    "selecting nothing": lambda m: m.load("BLOSUM62").select("").shape,
    "a two-letter matrix's values": lambda m: tuple(
        float(v) for v in m.load("HOXD70").values()
    ),
    "a two-letter matrix's keys": lambda m: m.load("HOXD70").keys(),
    "an integer read": lambda m: str(
        m.read(io.StringIO("A 1\nB 2\n"), dtype=int).dtype
    ),
}


@needs_biopython
@pytest.mark.parametrize("label", sorted(_CORNERS))
def test_every_corner_agrees_with_the_reference(label):
    """One call at a time, and the answer -- or the failure -- has to match.

    A table like this is only worth having if a case cannot pass by doing
    nothing, so the two answers are compared as strings after `repr`, which puts
    an exception's type and message on the same footing as a value.
    """
    reference, _ = reference_module()
    theirs = _answer(reference, _CORNERS[label])
    mine = _answer(ours, _CORNERS[label])
    assert mine == theirs


def _answer(module, call):
    try:
        return repr(call(module))
    except Exception as exc:
        return f"{type(exc).__name__}: {exc}"



def test_the_private_module_needs_nothing_but_numpy():
    """The derivation is loadable by path, which is what the generator does.

    `tools/gen_substitution_matrices.py` loads `_substitution_matrices.py` by
    path so that it can check the table it is about to write without importing
    the table it wrote last time.  That only works while the private module
    imports nothing but numpy, so the constraint is a test and not a comment.
    """
    source = (ROOT / "src" / "biofasting" / "_substitution_matrices.py").read_text()
    imports = {
        line.split()[1].split(".")[0]
        for line in source.splitlines()
        if line.startswith("import ") or line.startswith("from ")
    }
    assert imports == {"__future__", "contextlib", "string", "numpy"}, imports
