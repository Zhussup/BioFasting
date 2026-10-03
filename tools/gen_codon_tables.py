#!/usr/bin/env python3
"""Turn Biopython's 27 genetic codes into the 64-character strings NCBI prints.

`Bio.Data.CodonTable` is 1,308 lines, and roughly 900 of them are 27 calls to
`register_ncbi_table` each spelling out all 64 codons of a genetic code as a
Python dictionary.  The triage calls this `data` -- "becomes data + a parser;
no class port" -- and a genetic code is data of a very particular shape: NCBI
publishes each one as a 64-character string of amino acids in the codon order
TTT, TTC, TTA, TTG, TCT, ... GGG, with `*` where the code stops.  That string
says everything those 27 dictionaries say, in 64 bytes instead of 2 kB.

This script produces that string, from the reference implementation itself.
The name, the alternative name, the identifier and the start and stop codons
are read out of the reference's *source text* with `ast`, because they are
NCBI's data and the runtime objects have already lost them -- the module splits
`"Bacterial, Archaeal and Plant Plastid"` into three names and appends the
alternative, and the original string is not recoverable from the result.

It also **checks its own work**, in three ways, and refuses to write anything if
any of them fails:

1. The amino-acid string derived from the runtime objects must rebuild the
   runtime `forward_table` exactly, for all 27 codes.
2. The `table` dictionary written out in the reference's source text must equal
   the runtime `forward_table` -- so the source this script reads names is the
   source that is running.
3. Every one of the six tables built from the derived row must equal the
   reference's own object: the twelve public dictionaries, the forward table,
   the back table, the start and stop codons, the alphabets, and the ambiguous
   translation of all 3,375 ambiguous codons in each of the six alphabets,
   including which of them raise and with which exception.

    python3 tools/gen_codon_tables.py            # write the data module
    python3 tools/gen_codon_tables.py --check    # verify without writing

Requires Biopython, which is a build-time tool here and stays a runtime
non-dependency: the generated module is what ships.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "src" / "biofasting" / "_codon_tables_data.py"
DERIVATION = ROOT / "src" / "biofasting" / "_codon_table.py"


def _load_derivation():
    """Import the shared derivation without importing the package.

    `_codon_table.py` has no imports of its own, so loading it by path keeps
    this tool usable with nothing but a Python interpreter and Biopython --
    no build step, and no second copy of the derivation to drift.  It also
    cannot import the data file this script writes, which is the point.
    """
    spec = importlib.util.spec_from_file_location("_codon_table", DERIVATION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


codon_table = _load_derivation()


def load_reference():
    """The reference's codon tables, and the source text that made them."""
    try:
        from Bio.Data import CodonTable as reference
    except ImportError as exc:  # pragma: no cover - build-time tool
        raise SystemExit(
            f"Biopython is required to regenerate this data ({exc}).\n"
            "It is a build-time tool here; the generated module is what ships."
        ) from None
    return reference, Path(reference.__file__).read_text()


def source_calls(text: str) -> list[dict]:
    """The literal arguments of every `register_ncbi_table(...)` in the source."""
    calls = []
    for node in ast.walk(ast.parse(text)):
        if not isinstance(node, ast.Call):
            continue
        if not (isinstance(node.func, ast.Name) and node.func.id == "register_ncbi_table"):
            continue
        calls.append({kw.arg: ast.literal_eval(kw.value) for kw in node.keywords})
    return calls


def ncbi_version(text: str) -> str | None:
    """The NCBI table revision the reference says its data came from."""
    for line in text.splitlines():
        marker = "# Data from NCBI genetic code table version"
        if line.strip().startswith(marker):
            return line.strip()[len(marker):].strip()
    return None


def derive_row(reference, call: dict) -> tuple:
    """One genetic code, as the fields a row of the data file carries.

    The amino-acid string and the dual-coding codons are derived from the
    *runtime* objects rather than read from the source `table`, so that check 1
    below is a comparison between two independent things rather than a
    tautology.
    """
    table = reference.unambiguous_dna_by_id[call["id"]]
    dual = tuple(
        (codon, table.forward_table[codon])
        for codon in sorted(table.stop_codons)
        if codon in table.forward_table
    )
    dual_codons = {codon for codon, _ in dual}
    aa_string = "".join(
        "*" if (codon in table.stop_codons and codon not in dual_codons)
        else table.forward_table[codon]
        for codon in codon_table.CODON_ORDER
    )
    return (
        call["id"],
        call["name"],
        call["alt_name"],
        aa_string,
        tuple(call["start_codons"]),
        tuple(table.stop_codons),
        dual,
    )


def check_derivation(reference, calls: list[dict], rows: list[tuple]) -> list[str]:
    """The two checks that hold without building anything: string and source."""
    problems = []
    for call, row in zip(calls, rows):
        id, _, _, aa_string, _, _, dual = row
        rebuilt = codon_table.forward_table(aa_string, dual)
        runtime = reference.unambiguous_dna_by_id[id].forward_table
        if rebuilt != runtime:
            differing = sorted(
                codon for codon in codon_table.CODON_ORDER
                if rebuilt.get(codon) != runtime.get(codon)
            )
            problems.append(
                f"table {id}: the derived string rebuilds {len(differing)} "
                f"codon(s) differently, first {differing[:4]}"
            )
        if call["table"] != runtime:
            problems.append(f"table {id}: the source literal is not the table being run")
    return problems


def probe_codons(label: str) -> tuple:
    """Every base worth asking a table of this flavour about.

    For the explicit tables that is the four letters of their alphabet -- both
    alphabets for the generic table, which is keyed by T and U at once.  For the
    ambiguous tables it is every ambiguity code, which is 3,375 codons for DNA
    and RNA and 4,913 for the generic table, whose alphabet also carries a `T`.

    The test is `startswith("ambiguous_")` and not `"ambiguous" in label`, which
    is true of `"unambiguous_dna_by_id"` as well and would quietly ask the
    unambiguous tables only about ambiguous codons -- a check that still passes
    and no longer covers what it says it covers.
    """
    if label.startswith("ambiguous_"):
        if "generic" in label:
            return tuple(sorted(set(codon_table.AMBIGUOUS_RNA_VALUES) | {"T"}))
        if "rna" in label:
            return tuple(sorted(codon_table.AMBIGUOUS_RNA_VALUES))
        return tuple(sorted(codon_table.AMBIGUOUS_DNA_VALUES))
    if label == "generic_by_id":
        return tuple("UCAG") + tuple("TCAG")
    if label == "unambiguous_rna_by_id":
        return tuple("UCAG")
    return tuple("TCAG")


def probe_all(table, alphabet) -> dict:
    """What a table answers for every codon over `alphabet`, refusals included."""
    return {
        codon: _probe(table, codon)
        for codon in (
            a + b + c for a in alphabet for b in alphabet for c in alphabet
        )
    }


def check_tables(reference, rows: list[tuple]) -> list[str]:
    """Every table rebuilt from a row, against every table the reference has.

    This is check 3, and it is the one that decides whether the data file is any
    good: it builds all six flavours of each of the 27 codes from the fields the
    file will carry, and compares them with the reference's own -- the printed
    grid, the names, the back table, the start and stop codons, the alphabets,
    and then the translation of every codon, refusals included.

    The forwarding tables are compared by *asking them*, not with `==`.  The
    ambiguous ones are `AmbiguousForwardTable` objects, which define no equality
    at all, so `==` between two of them is identity and always false -- a
    comparison that would fail on a perfect port and pass on nothing.  The same
    goes for `repr`, which the ambiguous classes do not define either; their
    `repr` is the address of the object and means nothing, so it is not
    compared.
    """
    problems = []
    families = (
        (reference.unambiguous_dna_by_id, "unambiguous_dna_by_id"),
        (reference.unambiguous_rna_by_id, "unambiguous_rna_by_id"),
        (reference.generic_by_id, "generic_by_id"),
        (reference.ambiguous_dna_by_id, "ambiguous_dna_by_id"),
        (reference.ambiguous_rna_by_id, "ambiguous_rna_by_id"),
        (reference.ambiguous_generic_by_id, "ambiguous_generic_by_id"),
    )
    by_name_families = (
        (reference.unambiguous_dna_by_name, "unambiguous_dna_by_name"),
        (reference.unambiguous_rna_by_name, "unambiguous_rna_by_name"),
        (reference.generic_by_name, "generic_by_name"),
        (reference.ambiguous_dna_by_name, "ambiguous_dna_by_name"),
        (reference.ambiguous_rna_by_name, "ambiguous_rna_by_name"),
        (reference.ambiguous_generic_by_name, "ambiguous_generic_by_name"),
    )

    built_names = {}
    for row in rows:
        id, name, alt_name, aa_string, start, stop, dual = row
        names, *tables = codon_table.build(id, name, alt_name, aa_string, start, stop, dual)
        built_names[id] = names

        for ours, (family, label) in zip(tables, families):
            theirs = family[id]
            ambiguous = label.startswith("ambiguous_")
            if str(ours) != str(theirs):
                problems.append(f"table {id} / {label}: the printed grid differs")
            for field in ("id", "names", "back_table", "start_codons", "stop_codons",
                          "nucleotide_alphabet", "protein_alphabet"):
                if getattr(ours, field) != getattr(theirs, field):
                    problems.append(f"table {id} / {label}: {field} differs")
            if not ambiguous:
                if ours.forward_table != theirs.forward_table:
                    problems.append(f"table {id} / {label}: forward_table differs")
                if repr(ours) != repr(theirs):
                    problems.append(f"table {id} / {label}: repr differs")

            alphabet = probe_codons(label)
            got = probe_all(ours.forward_table, alphabet)
            want = probe_all(theirs.forward_table, alphabet)
            differ = [codon for codon in got if got[codon] != want[codon]]
            if differ:
                first = differ[0]
                problems.append(
                    f"table {id} / {label}: {len(differ)} codon(s) translate "
                    f"differently, first {first!r} -> {got[first]!r} "
                    f"against the reference's {want[first]!r}"
                )

    for id, names in built_names.items():
        for family, label in by_name_families:
            for name in names:
                if name not in family:
                    problems.append(f"table {id} / {label}: {name!r} is missing")
                elif family[name].id != id:
                    problems.append(f"table {id} / {label}: {name!r} points elsewhere")
    for family, label in by_name_families:
        expected = {
            name for id, names in built_names.items()
            for name in names
        }
        if set(family) != expected:
            extra = sorted(set(family) - expected)
            problems.append(f"{label}: {len(extra)} name(s) not accounted for, e.g. {extra[:4]}")
    return problems


def _probe(table, codon):
    """The answer, or the name of the refusal -- never one mistaken for the other."""
    try:
        return ("ok", table[codon])
    except KeyError:
        return ("KeyError", None)
    except Exception as exc:  # TranslationError and anything else it grows
        return (type(exc).__name__, None)


HEADER = '''"""The 27 NCBI genetic codes, as the strings NCBI publishes them.

**Generated** by `tools/gen_codon_tables.py` -- do not edit.  Each row is

    (id, name, alt_name, amino_acids, start_codons, stop_codons, dual)

where `amino_acids` is NCBI's 64 characters in the codon order TTT, TTC, TTA,
TTG, ... GGG, `*` marking a stop, and `dual` lists the codons that are a stop
*and* code for an amino acid -- context-dependent stops, which a string of this
shape cannot express and which the karyorelictid ciliates and their relatives
really have.

Read by `biofasting.codon_table`; the derivation lives in `_codon_table.py`.
"""

from __future__ import annotations

'''


def render(rows: list[tuple], version: str | None) -> str:
    lines = [HEADER.rstrip("\n")]
    lines.append(f"NCBI_TABLE_VERSION = {version!r}")
    lines.append("")
    lines.append("#   id, name, alt_name, amino acids in CODON_ORDER, starts, stops, dual")
    lines.append("TABLES: list[tuple] = [")
    for id, name, alt_name, aa, start, stop, dual in rows:
        lines.append(f"    ({id!r}, {name!r}, {alt_name!r},")
        lines.append(f"     {aa!r},")
        lines.append(f"     {list(start)!r}, {list(stop)!r}, {list(dual)!r}),")
    lines.append("]")
    lines.append("")
    lines.append('__all__ = ["TABLES", "NCBI_TABLE_VERSION"]')
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gen_codon_tables.py")
    parser.add_argument("--check", action="store_true",
                        help="verify the derivation and the file on disk, write nothing")
    args = parser.parse_args(argv)

    reference, text = load_reference()
    calls = source_calls(text)
    if len(calls) != len(reference.generic_by_id):
        print(f"error: found {len(calls)} register_ncbi_table calls in the source "
              f"and {len(reference.generic_by_id)} tables at runtime", file=sys.stderr)
        return 1

    rows = [derive_row(reference, call) for call in calls]

    problems = check_derivation(reference, calls, rows)
    if problems:
        print(f"error: the derived amino-acid string is not the reference's table, "
              f"{len(problems)} problem(s):", file=sys.stderr)
        for problem in problems[:20]:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print(f"derivation: {len(rows)} genetic codes rebuilt from one 64-character string each")

    problems = check_tables(reference, rows)
    if problems:
        print(f"error: the rebuilt tables differ from the reference, "
              f"{len(problems)} problem(s):", file=sys.stderr)
        for problem in problems[:20]:
            print(f"  {problem}", file=sys.stderr)
        return 1
    probes = sum(
        len(probe_codons(label)) ** 3
        for label in ("unambiguous_dna_by_id", "unambiguous_rna_by_id", "generic_by_id",
                      "ambiguous_dna_by_id", "ambiguous_rna_by_id",
                      "ambiguous_generic_by_id")
    )
    print(f"tables: all {len(rows)} codes identical in 6 alphabets, "
          f"{probes * len(rows)} codons asked one at a time")

    output = render(rows, ncbi_version(text))
    digest = hashlib.sha256(output.encode()).hexdigest()

    if args.check:
        on_disk = OUT.read_text() if OUT.exists() else ""
        if on_disk != output:
            print(f"error: {OUT} does not match what the generator produces", file=sys.stderr)
            return 1
        print(f"{OUT.relative_to(ROOT)} is current ({len(rows)} codes, sha256 {digest[:16]})")
        return 0

    OUT.write_text(output)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(rows)} genetic codes, "
          f"{len(output) // 1024} KiB, sha256 {digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
