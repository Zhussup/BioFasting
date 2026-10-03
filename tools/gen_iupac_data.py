#!/usr/bin/env python3
"""Write out the IUPAC tables `Bio.SeqUtils` needs, and check they are enough.

`Bio.Data.IUPACData` is 420 lines of tables: the atom-count-derived residue
weights, the one-letter/three-letter protein names, the ambiguity expansion of
every nucleotide code, and the two complement translate tables.  Nothing in it
is a loop worth porting -- it is the third `data` block the triage asked for,
after the restriction enzymes and the genetic codes, and it exists in this
package for the same reason: a caller who wants `molecular_weight` should not
have to install 156,755 lines to get it.

What is *written* is the tables themselves.  What is *checked* is that they are
sufficient -- that a function built out of nothing but these dicts reproduces
the reference's own answers exactly, over every argument combination the
reference accepts:

1. every table equals the reference's, float for float;
2. the three-letter table is the inverse of the one-letter table, which is a
   property the reference's two literals have to agree on and could stop
   agreeing on;
3. the complement strings answer `Bio.Seq.complement` and `complement_rna` for
   all 256 byte values, not just the four letters anyone tests by hand;
4. the weight tables, read the way `molecular_weight` reads them, reproduce it
   over a probe corpus in all eight argument combinations -- so they are
   provably complete for the sequences they claim to cover, and the missing
   letters are missing on purpose;
5. `_gc_values` reproduces `gc_fraction(seq, "weighted")`, which is the only
   thing that reads it.

The reimplementation in check 4 is deliberately a *second* implementation: the
port runs a compiled kernel over a flat table, and this runs the reference's own
generator expression over these dicts.  A check that called the port would only
prove the port agrees with itself.

    python3 tools/gen_iupac_data.py            # write the data module
    python3 tools/gen_iupac_data.py --check    # verify without writing

Requires Biopython, which is a build-time tool here and stays a runtime
non-dependency: the generated module is what ships.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "src" / "biofasting" / "_iupac_data.py"

#: Sequences every weight and GC check is run over.  Deliberately mixed: both
#: cases, every unambiguous letter, an ambiguous letter that must raise, a
#: gap, whitespace in the middle (which `molecular_weight` strips), and the two
#: degenerate lengths -- one base, where the `(len(seq) - 1) * water` term goes
#: to zero, and the empty sequence, where it goes negative and the answer is a
#: bare water molecule.
#:
#: Then every letter on its own, in both cases.  That is the part that makes the
#: *absence* of the ambiguous letters a checked property: a table with an extra
#: entry answers a weight where the reference raises, and a table missing one
#: raises where the reference answers, and neither shows up on a corpus of ACGT.
PROBES = (
    "",
    "A",
    "ACGT",
    "acgt",
    "ACGTacgtACGTacgtACGTacgt",
    " A C G T\n",
    "ACGTN",
    "ACGT-",
    "AUGCAUGC",
    "augcaugc",
    "ACDEFGHIKLMNPQRSTVWY",
    "acdefghiklmnpqrstvwy",
    "ACDEFGHIKLMNPQRSTVWYBXZJUO",
    "MAMAMAMA",
) + tuple(
    chr(value)
    for value in range(0x20, 0x7F)
)


def load_reference():
    """The reference's data module, its GC table, and the two Seq tables."""
    try:
        from Bio.Data import IUPACData
        from Bio.Seq import _dna_complement_table, _rna_complement_table
        from Bio.SeqUtils import _gc_values
    except ImportError as exc:  # pragma: no cover - build-time tool
        raise SystemExit(
            f"Biopython is required to regenerate this data ({exc}).\n"
            "It is a build-time tool here; the generated module is what ships."
        ) from None
    return IUPACData, _dna_complement_table, _rna_complement_table, _gc_values


def tables(iupac, dna_table, rna_table, gc_values):
    """Everything the generated module carries, as plain Python objects."""
    return {
        "PROTEIN_1TO3": dict(iupac.protein_letters_1to3_extended),
        "PROTEIN_3TO1": dict(iupac.protein_letters_3to1_extended),
        "AMBIGUOUS_DNA_VALUES": dict(iupac.ambiguous_dna_values),
        "GC_VALUES": dict(gc_values),
        "DNA_WEIGHTS": dict(iupac.unambiguous_dna_weights),
        "DNA_WEIGHTS_MONOISOTOPIC": dict(iupac.monoisotopic_unambiguous_dna_weights),
        "RNA_WEIGHTS": dict(iupac.unambiguous_rna_weights),
        "RNA_WEIGHTS_MONOISOTOPIC": dict(iupac.monoisotopic_unambiguous_rna_weights),
        "PROTEIN_WEIGHTS": dict(iupac.protein_weights),
        "PROTEIN_WEIGHTS_MONOISOTOPIC": dict(iupac.monoisotopic_protein_weights),
        "DNA_COMPLEMENT": bytes(dna_table).decode("latin-1"),
        "RNA_COMPLEMENT": bytes(rna_table).decode("latin-1"),
    }


# ---------------------------------------------------------------------------
# The checks.  Each returns a list of human-readable problems; empty is good.
# ---------------------------------------------------------------------------


def check_tables(data: dict, iupac, dna_table, rna_table, gc_values) -> list[str]:
    """Check 1: the file says what the reference says, bit for bit."""
    problems = []
    expected = {
        "PROTEIN_1TO3": iupac.protein_letters_1to3_extended,
        "PROTEIN_3TO1": iupac.protein_letters_3to1_extended,
        "AMBIGUOUS_DNA_VALUES": iupac.ambiguous_dna_values,
        "GC_VALUES": gc_values,
        "DNA_WEIGHTS": iupac.unambiguous_dna_weights,
        "DNA_WEIGHTS_MONOISOTOPIC": iupac.monoisotopic_unambiguous_dna_weights,
        "RNA_WEIGHTS": iupac.unambiguous_rna_weights,
        "RNA_WEIGHTS_MONOISOTOPIC": iupac.monoisotopic_unambiguous_rna_weights,
        "PROTEIN_WEIGHTS": iupac.protein_weights,
        "PROTEIN_WEIGHTS_MONOISOTOPIC": iupac.monoisotopic_protein_weights,
        "DNA_COMPLEMENT": bytes(dna_table).decode("latin-1"),
        "RNA_COMPLEMENT": bytes(rna_table).decode("latin-1"),
    }
    for name, want in expected.items():
        got = data[name]
        if got != want:
            if isinstance(got, dict):
                differ = sorted(k for k in set(got) | set(want) if got.get(k) != want.get(k))
                problems.append(f"{name}: {len(differ)} key(s) differ, e.g. {differ[:4]}")
            else:
                problems.append(f"{name}: differs from the reference")
    return problems


def check_inverse(data: dict) -> list[str]:
    """Check 2: the three-letter table is the inverse of the one-letter one.

    The reference writes both out by hand, so they can drift apart, and a drift
    shows up as `seq1(seq3(x)) != x` for exactly one residue -- a bug that a
    round-trip test over a handful of sequences is unlikely to hit.
    """
    forward = data["PROTEIN_1TO3"]
    declared = data["PROTEIN_3TO1"]
    derived = {three: one for one, three in forward.items()}
    problems = []
    if derived != declared:
        differ = sorted(set(derived) | set(declared))
        problems.append(
            "PROTEIN_3TO1 is not the inverse of PROTEIN_1TO3: "
            + ", ".join(
                f"{k!r} -> {derived.get(k)!r} against {declared.get(k)!r}"
                for k in differ
                if derived.get(k) != declared.get(k)
            )[:300]
        )
    if len(declared) != len(forward):
        problems.append(
            f"the two tables have different sizes: {len(forward)} against {len(declared)}"
        )
    return problems


def check_complements(data: dict) -> list[str]:
    """Check 3: the translate strings answer the reference, byte by byte.

    The string form of `Bio.Seq.complement` re-encodes its argument as ASCII, so
    the reference refuses the upper half of the byte range outright; the `Seq`
    form takes bytes and does answer there, and it is the one asked.  Both forms
    go through the same table, so asking either is asking the reference.

    The upper half is also asserted to be the identity, which is not a detail:
    `complement` leaves every byte it does not know alone, so a port that mapped
    0x80 somewhere would produce a different *length* of answer on a file with a
    stray byte in it, and this is the cheapest place to say so.
    """
    from Bio.Seq import Seq, complement, complement_rna

    problems = []
    for label, table, reference in (
        ("DNA_COMPLEMENT", data["DNA_COMPLEMENT"], complement),
        ("RNA_COMPLEMENT", data["RNA_COMPLEMENT"], complement_rna),
    ):
        ours = table.encode("latin-1")
        if len(ours) != 256:
            problems.append(f"{label} is {len(ours)} characters, not 256")
            continue
        for value in range(256):
            source = chr(value)
            got = bytes([value]).translate(ours).decode("latin-1")
            if value < 128:
                want = reference(source)
            else:
                want = bytes(complement(Seq(bytes([value])))).decode("latin-1")
                if want != source:
                    problems.append(
                        f"the reference does not leave {source!r} alone: it "
                        f"answers {want!r}, so 'unknown bytes pass through' is "
                        "not the rule this table was copied from"
                    )
                    break
            if got != want:
                problems.append(
                    f"{label}: {source!r} -> {got!r} against the reference's {want!r}"
                )
                break
    return problems


#: `molecular_weight`'s own loop, written out here so that the weights are
#: checked against the function that consumes them rather than against
#: themselves.  The port does not call this; it calls a compiled kernel.
#:
#: `sum(...)` and not a running total, and that is not style.  Since CPython
#: 3.12 the built-in sums floats with Neumaier compensation, so on this
#: interpreter `sum(xs)` and a `for` loop that adds `xs` in the same order
#: disagree on about a quarter of all inputs.  The reference calls `sum`, so the
#: reference's answer is the compensated one -- and a kernel that adds naively
#: reproduces it only sometimes.
def reference_algorithm_mass(seq: str, table: dict, water: float) -> float:
    return sum(table[letter] for letter in seq) - (len(seq) - 1) * water


def check_weights(data: dict) -> list[str]:
    """Check 4: the weight tables reproduce `molecular_weight` exactly.

    Every combination of the four arguments, over the probe corpus plus every
    single letter of both alphabets, which is what makes the *absence* of the
    ambiguous letters a checked property rather than an untested one.
    """
    from Bio.SeqUtils import molecular_weight

    water = {"average": 18.0153, "monoisotopic": 18.010565}
    problems = []
    comparisons = 0
    for seq_type, table_key, complement_key in (
        ("DNA", "DNA_WEIGHTS", "DNA_COMPLEMENT"),
        ("RNA", "RNA_WEIGHTS", "RNA_COMPLEMENT"),
        ("protein", "PROTEIN_WEIGHTS", None),
    ):
        for monoisotopic, kind in ((False, "average"), (True, "monoisotopic")):
            table = data[table_key + ("_MONOISOTOPIC" if monoisotopic else "")]
            for double_stranded in (False, True):
                if double_stranded and seq_type == "protein":
                    continue
                for circular in (False, True):
                    for probe in PROBES:
                        # `molecular_weight` uppercases and strips whitespace
                        # before it looks anything up, and it reads a SeqRecord's
                        # `.seq` if it is handed one.
                        clean = "".join(probe.split()).upper()

                        def reference_answer():
                            try:
                                return molecular_weight(
                                    probe, seq_type, double_stranded, circular,
                                    monoisotopic,
                                )
                            except ValueError as exc:
                                return ("ValueError", str(exc))

                        def table_answer():
                            try:
                                total = reference_algorithm_mass(
                                    clean, table, water[kind])
                                if circular:
                                    total -= water[kind]
                                if double_stranded:
                                    mirrored = clean.encode("latin-1").translate(
                                        data[complement_key].encode("latin-1")
                                    ).decode("latin-1")
                                    total += reference_algorithm_mass(
                                        mirrored, table, water[kind])
                                    if circular:
                                        total -= water[kind]
                                return total
                            except KeyError as exc:
                                # The reference's own spelling, `from None`, and
                                # with the KeyError repr'd inside the quotes.
                                return ("ValueError",
                                        f"'{exc}' is not a valid unambiguous "
                                        f"letter for {seq_type}")

                        want = reference_answer()
                        got = table_answer()
                        comparisons += 1
                        if want != got:
                            problems.append(
                                f"{seq_type} {kind} ds={double_stranded} "
                                f"circ={circular} {probe!r}: {want!r} against {got!r}"
                            )
    print(f"weights: {comparisons} comparisons against SeqUtils.molecular_weight")
    return problems


def check_gc_values(data: dict) -> list[str]:
    """Check 5: `_gc_values` reproduces `gc_fraction(seq, "weighted")`.

    The weighted mode is the only caller of that table, and it is the mode with
    a number for every ambiguity code -- so a table that is missing one entry
    answers `0.0` where the reference raises `KeyError`, silently.
    """
    from Bio.SeqUtils import gc_fraction

    values = data["GC_VALUES"]
    weights = "BDHKMNRVXY"
    probes = list(PROBES) + ["".join(letters) for letters in ([c] for c in "ACGTSWBDHKMNRVXY")]
    probes += ["acgtswbdhkmnrvxy", "NNNNACGT", "SSSSSSSS", "WWWWWWWW"]
    problems = []
    comparisons = 0
    for probe in probes:
        try:
            want = gc_fraction(probe, "weighted")
        except (KeyError, ZeroDivisionError) as exc:
            want = (type(exc).__name__, str(exc))
        try:
            gc = sum(probe.count(x) for x in "CGScgs")
            length = len(probe)
            gc += sum(
                (probe.count(x) + probe.count(x.lower())) * values[x] for x in weights
            )
            got = 0 if length == 0 else gc / length
        except KeyError as exc:
            got = ("KeyError", exc.args[0])
        comparisons += 1
        if isinstance(want, tuple) or isinstance(got, tuple):
            if not (isinstance(want, tuple) and isinstance(got, tuple)):
                problems.append(f"weighted {probe!r}: {want!r} against {got!r}")
            elif want[0] != got[0]:
                problems.append(
                    f"weighted {probe!r}: the reference raises {want[0]}, "
                    f"the table raises {got[0]}"
                )
        elif want != got:
            problems.append(f"weighted {probe!r}: {want!r} against {got!r}")
    print(f"gc_values: {comparisons} comparisons against gc_fraction(.., 'weighted')")
    return problems


# ---------------------------------------------------------------------------
# Rendering.
# ---------------------------------------------------------------------------

HEADER = '''"""The IUPAC tables `biofasting.sequtils` and `biofasting.checksum` need.

**Generated** by `tools/gen_iupac_data.py` -- do not edit.  Twelve tables copied
from `Bio.Data.IUPACData` and `Bio.SeqUtils._gc_values`, with the two complement
tables flattened from 256-byte `bytes.translate` strings into `str`, because
that is the shape `str.translate` wants.

The generator refuses to write unless it can rebuild the reference's own answers
from these tables alone: `molecular_weight` over every argument combination,
`gc_fraction(.., "weighted")`, and `complement`/`complement_rna` for all 256 byte
values.  The weight tables are the ones with something to get wrong -- a missing
letter is a `KeyError` the reference raises and this would answer `0.0` to.
"""
'''


def render(data: dict, reference_version: str) -> str:
    lines = [HEADER.rstrip("\n"), ""]
    lines.append(f"#: The Biopython release these were read from: {reference_version}")
    lines.append(f"IUPAC_SOURCE_VERSION = {reference_version!r}")
    lines.append("")

    def emit_dict(name: str, table: dict, comment: str) -> None:
        lines.append(f"# {comment}")
        lines.append(f"{name}: dict = {{")
        width = 0
        row = []
        for key, value in table.items():
            row.append(f"{key!r}: {value!r},")
            width += len(row[-1]) + 1
            if width > 72:
                lines.append("    " + " ".join(row))
                row, width = [], 0
        if row:
            lines.append("    " + " ".join(row))
        lines.append("}")
        lines.append("")

    emit_dict("PROTEIN_1TO3", data["PROTEIN_1TO3"],
              "One-letter to three-letter residue names, the extended alphabet.")
    emit_dict("PROTEIN_3TO1", data["PROTEIN_3TO1"],
              "The inverse, as the reference writes it out separately.")
    emit_dict("AMBIGUOUS_DNA_VALUES", data["AMBIGUOUS_DNA_VALUES"],
              "Every DNA ambiguity code and the letters it stands for.")
    emit_dict("GC_VALUES", data["GC_VALUES"],
              "The GC weight of every code, for gc_fraction's weighted mode.")
    emit_dict("DNA_WEIGHTS", data["DNA_WEIGHTS"],
              "Average residue masses of the four deoxynucleotides, 5'-phosphate.")
    emit_dict("DNA_WEIGHTS_MONOISOTOPIC", data["DNA_WEIGHTS_MONOISOTOPIC"],
              "The same, monoisotopic.")
    emit_dict("RNA_WEIGHTS", data["RNA_WEIGHTS"],
              "Average residue masses of the four ribonucleotides, 5'-phosphate.")
    emit_dict("RNA_WEIGHTS_MONOISOTOPIC", data["RNA_WEIGHTS_MONOISOTOPIC"],
              "The same, monoisotopic.")
    emit_dict("PROTEIN_WEIGHTS", data["PROTEIN_WEIGHTS"],
              "Average residue masses of the twenty protein letters.")
    emit_dict("PROTEIN_WEIGHTS_MONOISOTOPIC", data["PROTEIN_WEIGHTS_MONOISOTOPIC"],
              "The same, monoisotopic.")

    for name, comment in (
        ("DNA_COMPLEMENT", "Bio.Seq.complement, as a str.translate table over all 256 bytes."),
        ("RNA_COMPLEMENT", "Bio.Seq.complement_rna, the same way."),
    ):
        lines.append(f"# {comment}")
        lines.append(f"{name} = (")
        table = data[name]
        for start in range(0, 256, 64):
            chunk = table[start:start + 64]
            lines.append(f"    {chunk!r}")
        lines.append(")")
        lines.append("")

    lines.append("__all__ = [")
    for name in data:
        lines.append(f"    {name!r},")
    lines.append('    "IUPAC_SOURCE_VERSION",')
    lines.append("]")
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gen_iupac_data.py")
    parser.add_argument("--check", action="store_true",
                        help="verify the tables and the file on disk, write nothing")
    args = parser.parse_args(argv)

    iupac, dna_table, rna_table, gc_values = load_reference()
    data = tables(iupac, dna_table, rna_table, gc_values)

    checks = (
        ("tables", lambda: check_tables(data, iupac, dna_table, rna_table, gc_values)),
        ("inverse", lambda: check_inverse(data)),
        ("complements", lambda: check_complements(data)),
        ("weights", lambda: check_weights(data)),
        ("gc_values", lambda: check_gc_values(data)),
    )
    for label, check in checks:
        problems = check()
        if problems:
            print(f"error: the {label} check failed, {len(problems)} problem(s):",
                  file=sys.stderr)
            for problem in problems[:20]:
                print(f"  {problem}", file=sys.stderr)
            return 1
        if label not in ("weights", "gc_values"):
            print(f"{label}: the reference's own tables, unchanged")

    try:
        from importlib.metadata import version

        reference_version = version("biopython")
    except Exception:  # pragma: no cover - a source checkout without metadata
        reference_version = "unknown"

    output = render(data, reference_version)
    digest = hashlib.sha256(output.encode()).hexdigest()

    if args.check:
        on_disk = OUT.read_text() if OUT.exists() else ""
        if on_disk != output:
            print(f"error: {OUT} does not match what the generator produces",
                  file=sys.stderr)
            return 1
        print(f"{OUT.relative_to(ROOT)} is current "
              f"({sum(len(v) for v in data.values() if isinstance(v, dict))} entries, "
              f"sha256 {digest[:16]})")
        return 0

    OUT.write_text(output)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(data)} tables, "
          f"{len(output) // 1024} KiB, sha256 {digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
