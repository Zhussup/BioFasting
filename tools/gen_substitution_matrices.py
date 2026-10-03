#!/usr/bin/env python3
"""Turn Biopython's thirty substitution matrices into one Python data file.

`Bio.Align.substitution_matrices` is 508 lines and a `data/` directory of thirty
plain-text tables, 84,540 bytes of them.  The triage calls the module `data` --
"becomes data + a parser; no class port" -- and it is half right: the thirty
tables are data, and the 330 lines of `Array` around them are a class that has
to be written out.  What this script moves is the tables: out of thirty files in
a package directory and into one generated module, so that the numbers a caller
loads are in the library they installed rather than in a directory next to it.

The thirty tables are symmetric -- all thirty, which is a fact about the
matrices and not a coincidence -- so each one is stored as its lower triangle
`(0,0), (1,0), (1,1), (2,0), ...` and mirrored when it is loaded.  That is a
derivation with something to get wrong, and the check below is what makes it
safe: a matrix stored this way is compared against the reference's own, entry
for entry, before anything is written.

It **checks its own work**, in three ways, and refuses to write anything if any
of them fails:

1. Every file in the reference's data directory parses with the reference's own
   `read()` and with ours, and the two agree on the alphabet, the dtype, the
   header, every entry, and the string `format()` prints.
2. The rows this script is about to write must rebuild all thirty: the same
   16,140 numbers, the same alphabets, the same headers -- and it writes the
   triangle only for a matrix it has just shown to be symmetric, storing the
   full square for one that is not.
3. Under `--check`, the module a caller imports must then load each of its
   thirty names and agree again, so that a fault anywhere between the row and
   the caller is caught here rather than in a test.

    python3 tools/gen_substitution_matrices.py            # write the data module
    python3 tools/gen_substitution_matrices.py --check    # verify, write nothing

Requires Biopython, which is a build-time tool here and stays a runtime
non-dependency: the generated module is what ships.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "src" / "biofasting" / "_substitution_matrices_data.py"
DERIVATION = ROOT / "src" / "biofasting" / "_substitution_matrices.py"


def _load_derivation():
    """Import the shared derivation without importing the package.

    `_substitution_matrices.py` imports nothing but `numpy`, so loading it by
    path keeps this tool usable with nothing but an interpreter, numpy and
    Biopython -- no build step -- and it cannot import the data file this script
    writes, which is the point.
    """
    spec = importlib.util.spec_from_file_location("_substitution_matrices", DERIVATION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


matrices = _load_derivation()


def load_reference():
    """The reference's module and the list of files it ships."""
    try:
        from Bio.Align import substitution_matrices as reference
    except ImportError as exc:  # pragma: no cover - build-time tool
        raise SystemExit(
            f"Biopython is required to regenerate this data ({exc}).\n"
            "It is a build-time tool here; the generated module is what ships."
        ) from None
    directory = Path(reference.__file__).resolve().parent / "data"
    names = [
        name
        for name in sorted(os.listdir(directory))
        if name != "README.txt"
    ]
    return reference, directory, names


def derive_row(reference, directory: Path, name: str) -> tuple:
    """One matrix, as the fields a row of the data file carries.

    The numbers are taken from the text of the file as-is -- the tokens are kept
    as strings rather than turned into Python floats -- so that a caller's
    `float()` sees exactly the characters the reference's `read()` saw.  The
    header is kept stripped of its `#` and of the space after it, and the
    alphabet as a string when every letter is one character and as a tuple when
    it is not, which is how the reference's own `read()` decides it.
    """
    path = directory / name
    text = path.read_text()
    matrix = reference.read(str(path))

    lines = text.splitlines()
    i = 0
    while i < len(lines) and lines[i].startswith("#"):
        i += 1
    rows = [line.split() for line in lines[i:]]
    if len(rows[0]) == len(rows[1]) == 2:
        raise SystemExit(
            f"error: {name} is a one-dimensional matrix; this generator and the "
            "data format assume the square tables, and a 1-D row would need a "
            "field this file does not have"
        )
    rows.pop(0)  # the column labels, which `read()` has already taken
    #   `NUC.4.4` ends with a blank line.  The reference's parser walks the rows
    #   and the alphabet together with `zip`, so a row past the end of the
    #   alphabet is never looked at and a trailing blank is harmless there --
    #   while a blank *between* two rows would be popped and raise.  Same rule
    #   here: drop what is past the end, refuse what is in the middle.
    while rows and not rows[-1]:
        rows.pop()
    if any(not row for row in rows):
        blank = next(i for i, row in enumerate(rows) if not row)
        raise SystemExit(f"error: {name} has a blank line at row {blank}")
    square = []
    for row_number, row in enumerate(rows):
        #   The reference's parser walks the rows and the alphabets together and
        #   only asserts that the labels line up; a file whose labels did not
        #   would be read under the wrong letter there.  It is checked here
        #   rather than assumed, because the tokens this script keeps come from
        #   the text and not from the parsed object.
        letter = row.pop(0)
        if letter != matrix.alphabet[row_number]:
            raise SystemExit(
                f"error: {name} labels row {row_number} {letter!r} where the "
                f"header says {matrix.alphabet[row_number]!r}"
            )
        square.extend(row)

    n = len(matrix.alphabet)
    if len(square) != n * n:
        raise SystemExit(
            f"error: {name} holds {len(square)} numbers for a {n}x{n} matrix"
        )

    values = [float(word) for word in square]
    if values != list(matrix.view(matrices.np.ndarray).flat):
        raise SystemExit(f"error: {name} does not parse to the reference's own numbers")
    symmetric = all(
        matrix.view(matrices.np.ndarray)[i, j]
        == matrix.view(matrices.np.ndarray)[j, i]
        for i in range(n)
        for j in range(n)
    )
    if symmetric:
        indices = [(i, j) for i in range(n) for j in range(i + 1)]
    else:
        indices = [(i, j) for i in range(n) for j in range(n)]
    kept = [square[i * n + j] for i, j in indices]
    return (
        name,
        tuple(matrix.header),
        matrix.alphabet,
        symmetric,
        " ".join(kept),
    )


def compare(ours, theirs, name: str, where: str) -> list[str]:
    """Every way two matrices can differ, as a list of complaints."""
    problems = []
    if ours.alphabet != theirs.alphabet:
        problems.append(f"{name} / {where}: alphabet differs")
    if ours.dtype != theirs.dtype:
        problems.append(f"{name} / {where}: dtype is {ours.dtype} not {theirs.dtype}")
    if list(ours.header) != list(theirs.header):
        problems.append(f"{name} / {where}: header differs")
    if ours.shape != theirs.shape:
        problems.append(f"{name} / {where}: shape is {ours.shape} not {theirs.shape}")
    if not (np_array(ours) == np_array(theirs)).all():
        differing = [
            (i, j)
            for i in range(theirs.shape[0])
            for j in range(theirs.shape[1])
            if np_array(ours)[i, j] != np_array(theirs)[i, j]
        ]
        problems.append(
            f"{name} / {where}: {len(differing)} entr(ies) differ, first {differing[:4]}"
        )
    if str(ours) != str(theirs):
        problems.append(f"{name} / {where}: the printed matrix differs")
    if repr(ours) != repr(theirs):
        problems.append(f"{name} / {where}: repr differs")
    return problems


def np_array(matrix):
    return matrix.view(matrices.np.ndarray)


def check_every_file(reference, directory: Path, names: list[str]) -> list[str]:
    """Check 1: our parser against the reference's, on all thirty files."""
    problems = []
    for name in names:
        path = directory / name
        theirs = reference.read(str(path))
        ours = matrices.read(str(path))
        problems.extend(compare(ours, theirs, name, "read"))
        #   The parser is also handed the text as a stream, which is the other
        #   half of what `read` promises to accept.
        with open(path) as fp:
            ours = matrices.read(fp)
        problems.extend(compare(ours, theirs, name, "read from a handle"))
    return problems


def check_rows(reference, directory: Path, rows: list[tuple]) -> list[str]:
    """Check 2: the rows about to be written, rebuilt and compared."""
    problems = []
    for name, header, alphabet, symmetric, values in rows:
        ours = matrices.load_row(header, alphabet, symmetric, values)
        theirs = reference.read(str(directory / name))
        problems.extend(compare(ours, theirs, name, "the row"))
    return problems


def check_written(names: list[str]) -> list[str]:
    """Check 3: the same rows, through the module a caller will import."""
    import importlib

    package = importlib.import_module("biofasting.substitution_matrices")
    problems = []
    if package.load() != sorted(names):
        problems.append("load() does not list the thirty names")
    from Bio.Align import substitution_matrices as reference

    for name in names:
        ours = package.load(name)
        theirs = reference.load(name)
        problems.extend(compare(ours, theirs, name, "load"))
    return problems


HEADER = '''"""The thirty substitution matrices, as numbers rather than as thirty files.

**Generated** by `tools/gen_substitution_matrices.py` -- do not edit.  Each row is

    (name, header, alphabet, symmetric, values)

where `header` is the file's own leading `#` comment lines with the marker
stripped, `alphabet` is a string when every letter is one character and a tuple
when it is not, and `values` is one whitespace-separated string of the matrix's
numbers.  When `symmetric` is true -- and all thirty are -- `values` is the lower
triangle row by row, `(0,0)`, `(1,0)`, `(1,1)`, `(2,0)`, ...; otherwise it is the
full square row by row.  The numbers are the file's own tokens, kept as text so
that `float()` on them here is the conversion the reference made on them there.

Read by `biofasting.substitution_matrices`; the derivation lives in
`_substitution_matrices.py`.
"""

from __future__ import annotations

#   name, header, alphabet, symmetric, values
TABLES: list[tuple] = [
'''

FOOTER = ''']

__all__ = ["TABLES"]
'''


def render(rows: list[tuple]) -> str:
    """The data module, as text.

    One string per row of the matrix -- per row of the triangle when it is
    symmetric, which is one token more than the row before it -- so that a
    regenerated file differs from the last one in the rows that changed and not
    in a reflowed paragraph.  Nobody is meant to edit this file, but it should
    still be readable by whoever has to.
    """
    lines = [HEADER.rstrip("\n")]
    for name, header, alphabet, symmetric, values in rows:
        lines.append("    (")
        lines.append(f"        {name!r},")
        lines.append("        (")
        for line in header:
            lines.append(f"            {line!r},")
        lines.append("        ),")
        lines.append(f"        {alphabet!r},")
        lines.append(f"        {symmetric!r},")
        lines.append("        (")
        tokens = values.split()
        chunks = []
        position = 0
        row = 0
        while position < len(tokens):
            take = row + 1 if symmetric else len(alphabet)
            chunks.append(" ".join(tokens[position : position + take]))
            position += take
            row += 1
        #   The chunks are separate string literals, so the space between two of
        #   them has to be written: Python joins them without one, and
        #   `"1" "-2"` is `"1-2"`.
        for k, chunk in enumerate(chunks):
            separator = " " if k < len(chunks) - 1 else ""
            lines.append(f'            "{chunk}{separator}"')
        lines.append("        ),")
        lines.append("    ),")
    return "\n".join(lines) + "\n" + FOOTER


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gen_substitution_matrices.py")
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the derivation and the file on disk, write nothing",
    )
    args = parser.parse_args(argv)

    reference, directory, names = load_reference()
    rows = [derive_row(reference, directory, name) for name in names]
    full = sum(len(row[2]) ** 2 for row in rows)
    stored = sum(len(row[4].split()) for row in rows)
    print(
        f"derivation: {len(rows)} matrices, {full} numbers, "
        f"{sum(1 for row in rows if row[3])} of them symmetric, "
        f"{stored} stored"
    )

    problems = check_every_file(reference, directory, names)
    if problems:
        print(
            f"error: our parser is not the reference's, {len(problems)} problem(s):",
            file=sys.stderr,
        )
        for problem in problems[:20]:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print(f"parser: {len(names)} files read by both, identical entry for entry")

    problems = check_rows(reference, directory, rows)
    if problems:
        print(
            f"error: the rows do not rebuild the reference, {len(problems)} problem(s):",
            file=sys.stderr,
        )
        for problem in problems[:20]:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print(f"rows: {len(rows)} matrices rebuilt from {stored} stored numbers")

    output = render(rows)
    digest = hashlib.sha256(output.encode()).hexdigest()

    if args.check:
        on_disk = OUT.read_text() if OUT.exists() else ""
        if on_disk != output:
            print(
                f"error: {OUT} does not match what the generator produces",
                file=sys.stderr,
            )
            return 1
        print(f"{OUT.relative_to(ROOT)} is current ({len(rows)} matrices, sha256 {digest[:16]})")

        #   Check 3, which only makes sense once the file it describes exists:
        #   the module a caller imports, loading every name it advertises.  It
        #   needs the package built, so it runs here rather than in the write
        #   path, where the file is being written for the first time.
        problems = check_written(names)
        if problems:
            print(
                f"error: the module does not load what was written, "
                f"{len(problems)} problem(s):",
                file=sys.stderr,
            )
            for problem in problems[:20]:
                print(f"  {problem}", file=sys.stderr)
            return 1
        print(f"module: all {len(names)} names load through the package, identical again")
        return 0

    OUT.write_text(output)
    print(
        f"wrote {OUT.relative_to(ROOT)}: {len(rows)} matrices, "
        f"{len(output) // 1024} KiB, sha256 {digest[:16]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
