"""The thirty substitution matrices Biopython ships, as data in the library.

`Bio.Align.substitution_matrices` is a 508-line module and a directory of thirty
plain-text scoring tables -- BLOSUM62, PAM250, the NCBI nucleotide matrices, and
the rest.  This is the same thing with the tables in one generated Python module
instead: 4,860 numbers, each matrix stored as its lower triangle because all
thirty are symmetric, loaded by a reader that has been checked against the
reference's own entry for entry.

That is the `data` verdict from `inventory/TRIAGE.md` applied to this module.  The
matrices are data and live in `_substitution_matrices_data.py`; `Array` is
behaviour, and is written out in `_substitution_matrices.py`, because a table
cannot say what letter indexing does, when it refuses, or what `format()` prints.

    >>> from biofasting import substitution_matrices
    >>> substitution_matrices.load("BLOSUM62")["A", "R"]
    -1.0
    >>> substitution_matrices.load("BLOSUM62")["W", "W"]
    11.0
    >>> len(substitution_matrices.load())
    30

`load()` takes a name from that list, or a path to any matrix file in the same
format -- the reference resolves a name against its own data directory with
`os.path.join`, so an absolute path works there and works here.  Nothing here
imports `Bio`, so a program with no Biopython installed still gets its scoring
matrix.

The one difference from the reference worth stating plainly: Biopython keeps its
tables as files in `Bio/Align/substitution_matrices/data/`, and `load(name)` for a
name that is not there fails inside `open()` naming the full path it looked in.
There is no such directory here, so an unknown *path* still fails inside
`open()`, naming the path as given, and an unknown bare *name* fails with the
thirty names listed, which is more use than a path into a directory that does not
exist.
"""

from __future__ import annotations

from ._substitution_matrices import Array, build, read
from ._substitution_matrices_data import TABLES

__all__ = [
    "Array",
    "read",
    "load",
]

#   The name a matrix is stored under, which is the name the reference's file
#   has: `BLOSUM62`, `NUC.4.4`, `TRANS`.  A dict comprehension would be shorter
#   and would drop a duplicated name in silence, so it is written out.
_BY_NAME: dict = {}
for _row in TABLES:
    _BY_NAME[_row[0]] = _row


def _load_row(row):
    """One `Array` from a row of the data file, header and all."""
    name, header, alphabet, symmetric, values = row
    matrix = build(alphabet, symmetric, values)
    matrix.header = list(header)
    return matrix


def load(name=None):
    """Load and return a precalculated substitution matrix.

    With no argument, the thirty names, sorted -- which is what `os.listdir` and
    a sort give in the reference.

        >>> from biofasting import substitution_matrices
        >>> names = substitution_matrices.load()
        >>> names[:3]
        ['BENNER22', 'BENNER6', 'BENNER74']
        >>> names[-1]
        'TRANS'
    """
    if name is None:
        return sorted(_BY_NAME)
    try:
        row = _BY_NAME.get(name)
    except TypeError:
        #   Unhashable: a path-like object of some kind, and never a name.
        row = None
    if row is not None:
        return _load_row(row)
    #   Not one of the thirty, so it is a path.  The reference would hand the
    #   name to `read`, which opens it.
    try:
        return read(name)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"no substitution matrix named {name!r} and no such file; "
            f"load() with no argument lists the {len(_BY_NAME)} that ship here"
        ) from None
