#!/usr/bin/env python3
"""Turn Biopython's restriction-enzyme classes into the data they always were.

`Bio.Restriction` is 1,225 objects of which 1,112 are enzyme classes generated
one per enzyme from a dictionary -- 27k lines of Python that say, per enzyme,
a recognition site, two or four cut coordinates, an overhang, a frequency and a
list of suppliers.  The triage in `inventory/TRIAGE.md` calls this `data`:
"becomes data + a parser; never ported as code".  This script produces the data
half, from the reference implementation itself, so the two can be compared
rather than trusted.

It also **checks its own work**.  The search pattern is not copied out of the
dictionary; it is derived from the recognition site, and the derivation is
required to reproduce Biopython's own compiled pattern for every one of the
1,088 enzymes.  A derivation that silently disagreed on three enzymes would
produce a library that is wrong on three enzymes forever, and nothing would
notice, because the data file has no opinion about itself.

    python3 tools/gen_restriction_data.py            # write the data module
    python3 tools/gen_restriction_data.py --check    # verify without writing

Requires Biopython, which is a build-time tool here and stays a runtime
non-dependency: the generated module is what ships.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = ROOT / "src" / "biofasting" / "_restriction_data.py"
SITES = ROOT / "src" / "biofasting" / "_restriction_sites.py"


def _load_sites():
    """Import the shared derivation without importing the compiled package.

    `_restriction_sites.py` has no imports of its own, so loading it by path
    keeps this tool usable with nothing but a Python interpreter and Biopython
    installed -- no build step, and no second copy of the derivation to drift.
    """
    spec = importlib.util.spec_from_file_location("_restriction_sites", SITES)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sites = _load_sites()
pattern_for = sites.pattern_for
reverse_complement = sites.reverse_complement

def compsite(name: str, site: str, palindromic: bool) -> str:
    """Rebuild the reference's own search pattern for an enzyme from its site.

    The derivation is shared with the library (`_restriction_sites.py`), loaded
    by path at the top of this file.  Only the group names differ: the reference
    names them after the enzyme, the engine uses one short pair for all of them.
    """
    return sites.compsite(site, palindromic, top=name, bottom=f"{name}_as")


def load_reference():
    """The reference's enzyme table, its type classification and its suppliers."""
    try:
        from Bio.Restriction.Restriction_Dictionary import (
            rest_dict,
            suppliers,
            typedict,
        )
    except ImportError as exc:  # pragma: no cover - build-time tool
        raise SystemExit(
            f"Biopython is required to regenerate this data ({exc}).\n"
            "It is a build-time tool here; the generated module is what ships."
        ) from None

    # typedict maps a pseudo-type name to (bases, enzyme names).  Inverting it
    # gives every enzyme's base classes, which is the classification the
    # engine's behaviour flags are derived from rather than re-guessed.
    types: dict[str, tuple[str, ...]] = {}
    for _pseudo, (bases, enzymes) in typedict.items():
        for enzyme in enzymes:
            types[enzyme] = tuple(bases)

    missing = set(rest_dict) - set(types)
    if missing:
        raise SystemExit(f"enzymes with no type classification: {sorted(missing)[:10]}")

    # The supplier dictionary is ordered, and the reference's `supplier_list()`
    # returns names in that order rather than sorted, so the order is part of
    # the behaviour and is carried over as-is.
    suppliers_by_code = {code: entry[0] for code, entry in suppliers.items()}
    unknown_codes = {c for e in rest_dict.values() for c in e["suppl"]} - set(suppliers_by_code)
    if unknown_codes:
        raise SystemExit(f"enzymes cite supplier codes with no name: {sorted(unknown_codes)}")
    return rest_dict, types, suppliers_by_code


_CLASS = re.compile(r"\[([A-Z]+)\]")


def canonical(pattern: str) -> str:
    """Rewrite a character class's members into a fixed order.

    The reference writes the ambiguity code S (G or C) as ``[CG]`` and this
    derivation writes ``[GC]``.  They are the same pattern.  Comparing the
    strings would report 18 "mismatches" that are all spelling, which is worse
    than no check at all: a check that cries wolf on every run teaches its
    reader to skim past it, and the one real mismatch would go with it.  So the
    comparison is on meaning -- the members of each class, sorted -- and the
    derivation is held to that.
    """
    return _CLASS.sub(lambda m: "[" + "".join(sorted(m.group(1))) + "]", pattern)


def check_derivation(rest_dict, types) -> list[str]:
    """Every derived pattern must mean the same as the reference's own.

    This is the script's reason to exist.  The data file is generated, so the
    only thing standing between a wrong derivation and a library that is
    quietly wrong on a handful of enzymes is this comparison, made against the
    reference at generation time for all 1,088 of them.
    """
    problems = []
    for name, entry in rest_dict.items():
        palindromic = "Palindromic" in types[name]
        mine = compsite(name, entry["site"], palindromic)
        if canonical(mine) != canonical(entry["compsite"]):
            problems.append(f"{name}: derived {mine!r} != reference {entry['compsite']!r}")
    return problems


def build(rest_dict, types) -> dict:
    """The data table, one tuple per enzyme, in a documented field order."""
    data = {}
    for name, entry in sorted(rest_dict.items()):
        fst5, fst3, scd5, scd3, site = entry["charac"]
        data[name] = (
            site,
            fst5,
            fst3,
            scd5,
            scd3,
            entry["ovhg"],
            entry["ovhgseq"],
            entry["size"],
            entry["freq"],
            tuple(entry["suppl"]),
            entry["opt_temp"],
            entry["inact_temp"],
            entry["uri"],
            entry["id"],
            types[name],
        )
    return data


HEADER = '''"""The restriction enzymes, as data.

**Generated** by `tools/gen_restriction_data.py` from Biopython 1.88's
`Bio.Restriction.Restriction_Dictionary`, and not edited by hand.  The
generator is run again whenever the reference is upgraded, and it refuses to
write anything unless it can rebuild every one of the reference's own search
patterns from the recognition sites alone -- so the file cannot drift from the
reference silently.

Why data and not 1,112 classes: the classes differ only in their attributes, so
they are one class plus a table.  Regenerating them as code would be copying
27k lines of generated Python to say the same thing; the triage in
`inventory/TRIAGE.md` files this whole module under the `data` verdict for that
reason.

Field order, per enzyme:

    (site, fst5, fst3, scd5, scd3, ovhg, ovhgseq, size, freq,
     suppliers, opt_temp, inact_temp, uri, rebase_id, types)

``fst5``/``fst3`` are the first cut on the top and bottom strand, ``scd5``/
``scd3`` the second cut for the enzymes that cut twice, all relative to the
start of the recognition site; ``None`` where the enzyme is not characterised.
``ovhg`` is the overhang length, negative for a 5' overhang, positive for 3'
and zero for blunt; it is ``None`` for enzymes whose cut is unknown.  ``types``
is the reference's own base-class tuple, which is where the polymorphic
questions (palindromic? defined? methylable? commercially available?) are
answered without a second opinion.

``SUPPLIERS`` maps a one-letter supplier code to the supplier's name, in the
reference's own order, because ``Enzyme.supplier_list()`` returns names in that
order rather than sorted.  Which supplier carries which enzyme is *not* stored
twice: ``suppl`` above already says it per enzyme, and the reverse index is
derived on demand.
"""'''


def render(data: dict, suppliers_by_code: dict) -> str:
    lines = [HEADER, "", "ENZYMES: dict[str, tuple] = {"]
    for name, fields in data.items():
        site, fst5, fst3, scd5, scd3, ovhg, ovhgseq, size, freq, suppl, opt, inact, uri, rid, types = fields
        lines.append(
            f"    {name!r}: ({site!r}, {fst5!r}, {fst3!r}, {scd5!r}, {scd3!r}, "
            f"{ovhg!r}, {ovhgseq!r}, {size!r}, {freq!r}, {suppl!r}, "
            f"{opt!r}, {inact!r}, {uri!r}, {rid!r}, {types!r}),"
        )
    lines.append("}")
    lines.append("")
    lines.append("SUPPLIERS: dict[str, str] = {")
    for code, name in suppliers_by_code.items():
        lines.append(f"    {code!r}: {name!r},")
    lines.append("}")
    lines.append("")
    lines.append('__all__ = ["ENZYMES", "SUPPLIERS"]')
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="gen_restriction_data.py")
    parser.add_argument("--check", action="store_true",
                        help="verify the derivation and the file on disk, write nothing")
    args = parser.parse_args(argv)

    rest_dict, types, suppliers_by_code = load_reference()
    problems = check_derivation(rest_dict, types)
    if problems:
        print(f"error: the derived pattern does not match the reference for "
              f"{len(problems)} enzyme(s):", file=sys.stderr)
        for problem in problems[:20]:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print(f"derivation: {len(rest_dict)} patterns rebuilt from the sites, all identical")

    data = build(rest_dict, types)
    text = render(data, suppliers_by_code)
    digest = hashlib.sha256(text.encode()).hexdigest()

    if args.check:
        on_disk = OUT.read_text() if OUT.exists() else ""
        if on_disk != text:
            print(f"error: {OUT} does not match what the generator produces", file=sys.stderr)
            return 1
        print(f"{OUT.relative_to(ROOT)} is current ({len(data)} enzymes, "
              f"{len(suppliers_by_code)} suppliers, sha256 {digest[:16]})")
        return 0

    OUT.write_text(text)
    print(f"wrote {OUT.relative_to(ROOT)}: {len(data)} enzymes, "
          f"{len(suppliers_by_code)} suppliers, {len(text) // 1024} KiB, "
          f"sha256 {digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
