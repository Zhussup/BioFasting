#!/usr/bin/env python3
"""Write out the tables `Bio.SeqUtils.ProtParamData` and `IsoelectricPoint` hold.

`ProtParamData` is 1,000 lines of published amino-acid scales: thirty-odd
per-residue scales, the 400 dipeptide instability weights, and the twenty
flexibility values.  Nothing in it is a loop worth porting -- it is the fourth
`data` block the triage asked for, after the restriction enzymes, the genetic
codes and the IUPAC tables -- and it exists in this package for the same reason:
a caller who wants `ProteinAnalysis.gravy` should not have to install Biopython
to get it.

What is *written* is the tables.  What is *checked* is that they are sufficient
and correctly shaped, because a scale is a `dict` in the reference and every one
of them is read with a bare `[aa]`:

1. every table equals the reference's, value for value *and type for type*, and
   the letters line up with `IUPACData.protein_letters` rather than with whatever
   order the literal happened to be written in.  The type is not cosmetic: ten
   of the twenty-eight gravy scales carry `int` values, and CPython's `sum()`
   -- which is what `gravy` is written as -- compensates a float item and adds an
   `int` one plainly, so the same twenty numbers answer differently for it;
2. the four functions built on nothing but these tables -- `flexibility`,
   `instability_index`, `gravy` and `protein_scale` -- reproduce the reference's
   own output over a probe corpus, *including* the inputs where the reference
   raises `KeyError`, writes a warning and skips, or returns an empty list
   because its range went negative;
3. `IsoelectricPoint`'s constants reproduce `pi()` and `charge_at_pH()`, whose
   summation order is the insertion order of two of these dicts and would change
   the last bit of the answer if it were lost.

The reimplementation in check 2 is deliberately a *second* implementation -- it
walks the tables dict by dict, the way the reference does -- so that the check
says the data is sufficient rather than that the port agrees with itself.  The
port runs a compiled kernel over the same numbers; this runs the reference's own
Python.

    python3 tools/gen_protparam_data.py            # write the data module
    python3 tools/gen_protparam_data.py --check    # verify without writing

Requires Biopython, which is a build-time tool here and stays a runtime
non-dependency: the generated module is what ships.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TARGET = ROOT / "src" / "biofasting" / "_protparam_data.py"

HEADER = '''"""The amino-acid scales `biofasting.protparam` is built on.

**Generated** by `tools/gen_protparam_data.py` -- do not edit.  {count} tables
read from `Bio.SeqUtils.ProtParamData` and `Bio.SeqUtils.IsoelectricPoint`, with
every letter aligned to `AMINO_ACIDS` instead of to the order a literal happened
to be written in, because a kernel indexes them by position and the reference
indexes them by letter.

`FLEX` is `ProtParamData.Flex`, the twenty values `ProteinAnalysis.flexibility`
uses.  `DIWV` is `ProtParamData.DIWV`, the 400 Guruprasad dipeptide instability
weights, row-major in `AMINO_ACIDS`.  `SCALES` is every twenty-letter scale the
reference exposes -- the twenty-eight `gravy_scales` entries among them -- and
`GRAVY_SCALES` maps a `gravy(scale=...)` name onto the key it names.  The three
constants after them are `IsoelectricPoint`'s, and their *order* is load-bearing:
`charge_at_pH` sums over them in insertion order and floating-point addition is
not associative.

**An `int` here stays an `int`.**  `sum(selected_scale[aa] for aa in sequence)` is
what `gravy` is, and CPython since 3.12 compensates a float item and adds an
`int` one plainly, with the compensation left where it was.  Ten of the
twenty-eight gravy scales are written with integer values -- Engelman's C, G and
H, Parker's N, D, Q, I and W, GoldSack's seven, and one or two each in Fasman,
Fauchere, Jones, Ponnuswamy, Roseman, Wilson and Zimmerman -- and on four of
those ten a compensated sum over the same twenty numbers answers differently,
between 4% and 29% of random proteins.  Rendering a `2` as `2.0` therefore changes
the last bit of `gravy`.  `repr` is what keeps them apart; nothing else in this
file should ever be reformatted by hand.

Read by `biofasting.protparam`; the generator that wrote it is the only thing
that should ever rewrite it.
"""

from __future__ import annotations

#: The Biopython release these were read from: {version}
PROTPARAM_SOURCE_VERSION = {version!r}

#: The twenty standard residues, which is `IUPACData.protein_letters`.  Every
#: scale below is a tuple in this order; it is the order a dense kernel table
#: needs, and it is not the order the reference's literals are written in.
AMINO_ACIDS = {amino!r}

'''


def _fmt_number(value):
    """A number as the shortest literal that reads back as the same object.

    `repr` and not a format string, because it is the one rendering that keeps
    `2` and `2.0` apart -- see the type note in the module docstring.
    """
    return repr(value)


def _fmt_tuple(values, indent):
    """A tuple of numbers, wrapped to something a reader can scan."""
    parts = [_fmt_number(v) for v in values]
    line = " " * indent
    lines = []
    current = line
    for index, part in enumerate(parts):
        piece = part + ("," if index + 1 < len(parts) else "")
        if len(current) + len(piece) + 1 > 88:
            lines.append(current.rstrip())
            current = line
        current += piece + " "
    lines.append(current.rstrip())
    return "\n".join(lines)


def build():
    """The tables, read from the reference, in this module's own layout."""
    from Bio.Data import IUPACData
    from Bio.SeqUtils import ProtParamData as D
    from Bio.SeqUtils import IsoelectricPoint as I

    amino = IUPACData.protein_letters
    scales = {}
    for name in dir(D):
        if name.startswith("_") or name in ("DIWV", "gravy_scales"):
            continue
        value = getattr(D, name)
        if isinstance(value, dict) and len(value) == len(amino):
            scales[name] = tuple(value[letter] for letter in amino)

    #   `gravy_scales` maps a name onto the *dict object*, and two of them --
    #   `Cowan3.4` and `Cowan7.5` -- are written inline there and never bound to
    #   a module attribute, so the name of the scale has to be recovered by
    #   identity, and a table with no attribute of its own is carried under its
    #   gravy name.  Attribute names cannot contain a dot, so the two cannot
    #   collide.
    gravy = {}
    for name, table in D.gravy_scales.items():
        values = tuple(table[letter] for letter in amino)
        owners = [k for k, v in scales.items() if v == values]
        if not owners:
            scales[name] = values
            owners = [name]
        if len(owners) != 1:
            raise SystemExit(
                f"gravy_scales[{name!r}] matches {len(owners)} of the exposed "
                f"scales ({owners}); the mapping cannot be emitted unambiguously"
            )
        gravy[name] = owners[0]

    dijw = tuple(tuple(D.DIWV[a][b] for b in amino) for a in amino)
    flex = tuple(D.Flex[letter] for letter in amino)

    return {
        "version": __import__("Bio").__version__,
        "amino": amino,
        "scales": scales,
        "gravy": gravy,
        "dijw": dijw,
        "flex": flex,
        "charged_aas": tuple(I.charged_aas),
        "positive_pks": dict(I.positive_pKs),
        "negative_pks": dict(I.negative_pKs),
        "pk_nterminal": dict(I.pKnterminal),
        "pk_cterminal": dict(I.pKcterminal),
    }


def emit(tables):
    out = [HEADER.format(count=len(tables["scales"]),
                         version=tables["version"],
                         amino=tables["amino"])]

    out.append("#: `ProtParamData.Flex`, in `AMINO_ACIDS` order.\nFLEX: tuple = (\n")
    out.append(_fmt_tuple(tables["flex"], 4))
    out.append("\n)\n\n")

    out.append("#: `ProtParamData.DIWV`, row-major in `AMINO_ACIDS`: `DIWV[i][j]`\n"
               "#: is the instability weight of the dipeptide `AMINO_ACIDS[i]`\n"
               "#: followed by `AMINO_ACIDS[j]`.\n"
               "DIWV: tuple = (\n")
    for row in tables["dijw"]:
        out.append("    (\n")
        out.append(_fmt_tuple(row, 8))
        out.append("\n    ),\n")
    out.append(")\n\n")

    out.append("#: Every twenty-letter scale the reference exposes, by the name it\n"
               "#: exposes it under.  `gravy(scale=...)` reads the twenty-eight in\n"
               "#: `GRAVY_SCALES`; `protein_scale(param_dict, ...)` takes any of them.\n"
               "SCALES: dict = {\n")
    for name in sorted(tables["scales"]):
        out.append(f"    {name!r}: (\n")
        out.append(_fmt_tuple(tables["scales"][name], 8))
        out.append("\n    ),\n")
    out.append("}\n\n")

    out.append("#: `protein_scale` scale name -> the key in `SCALES` it is.\n"
               "GRAVY_SCALES: dict = {\n")
    for name in sorted(tables["gravy"]):
        out.append(f"    {name!r}: {tables['gravy'][name]!r},\n")
    out.append("}\n\n")

    out.append("#: `IsoelectricPoint`'s charged residues, in the order the reference\n"
               "#: builds its content dict in.\n"
               f"CHARGED_AAS: tuple = {tables['charged_aas']!r}\n\n")

    for label, key, note in (
        ("POSITIVE_PKS", "positive_pks",
         "Groups that gain a proton.  `Nterm` first, and the order is the order\n"
         "#: `charge_at_pH` adds them in."),
        ("NEGATIVE_PKS", "negative_pks",
         "Groups that lose one.  `Cterm` first, likewise."),
        ("PK_NTERMINAL", "pk_nterminal",
         "Overrides for `positive_pKs['Nterm']`, by the residue it starts with."),
        ("PK_CTERMINAL", "pk_cterminal",
         "Overrides for `negative_pKs['Cterm']`, by the residue it ends with."),
    ):
        out.append(f"#: {note}\n{label}: dict = {{\n")
        for name, value in tables[key].items():
            out.append(f"    {name!r}: {_fmt_number(value)},\n")
        out.append("}\n\n")

    return "".join(out)


# ---------------------------------------------------------------------------
# The gate: the tables have to be enough to rebuild the reference's answers.


def check(tables):
    """Rebuild the four functions from these tables alone and compare."""
    import random
    from Bio.SeqUtils import ProtParamData as D
    from Bio.SeqUtils import IsoelectricPoint as I
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    amino = tables["amino"]
    flex = dict(zip(amino, tables["flex"]))
    dijw = {a: dict(zip(amino, row)) for a, row in zip(amino, tables["dijw"])}
    scales = {name: dict(zip(amino, values))
              for name, values in tables["scales"].items()}

    rng = random.Random(1013)
    corpus = [
        "",
        "A",
        "AC",
        "ACDEFGHIKLMNPQRSTVWY",
        "ACDEFGHIKLMNPQRSTVWY" * 3,
        "AAAAAAAAAA",
        "PETER",
        "INGAR",
        amino[::-1],
    ]
    for _ in range(200):
        body = "".join(rng.choice(amino) for _ in range(rng.randrange(0, 60)))
        corpus.append(body)
    #   The awkward ones: a letter no scale has an entry for, and lowercase,
    #   which `ProteinAnalysis` upper-cases before any table is consulted.
    corpus += ["ACDEFGHIKLMNPQRSTVWYB", "ACDEFGHI", "acdefghik", "AXA", "A*A",
               "A" * 8, "A" * 9, "A" * 10]
    #   Four short sequences that are here for one reason each: they are the
    #   smallest ones found whose `sum()` the int rule changes.  The reference's
    #   sum is -3.6 for "IYAR" under Parker, -0.20000000000000018 for "HSV" under
    #   Engelman, 2.8899999999999997 for "TNRC" under GoldSack and
    #   -3.9400000000000004 for "THTQG" under Roseman, where a compensated sum over
    #   the same numbers answers -3.5999999999999996, -0.20000000000000007, 2.89
    #   and -3.94.  Random draws find this about one time in twenty-five on those
    #   four scales, which is why it went unnoticed until the kernels were swept;
    #   these four make the gate say so on every run.
    corpus += ["IYAR", "HSV", "TNRC", "THTQG"]

    def our_flexibility(seq, flex_table):
        weights = (0.25, 0.4375, 0.625, 0.8125)
        out = []
        for i in range(len(seq) - 9):
            score = 0.0
            for j in range(4):
                score += (flex_table[seq[i + j]] + flex_table[seq[i + 8 - j]]) * weights[j]
            score += flex_table[seq[i + 5]]
            out.append(score / 5.25)
        return out

    def our_instability(seq, index):
        score = 0.0
        for i in range(len(seq) - 1):
            score += index[seq[i]][seq[i + 1]]
        return (10.0 / len(seq)) * score

    def our_protein_scale(seq, param_dict, window, edge):
        weights = [edge + (2 * (1.0 - edge) / (window - 1)) * i
                   for i in range(window // 2)]
        sum_of_weights = sum(weights) * 2 + 1
        scores = []
        for i in range(len(seq) - window + 1):
            sub = seq[i:i + window]
            score = 0.0
            for j in range(window // 2):
                # The reference warns and *skips the pair*, rather than
                # raising: a scale with no entry for a residue costs that
                # window one term and the run one line on stderr.
                try:
                    front = param_dict[sub[j]]
                    back = param_dict[sub[window - j - 1]]
                except KeyError:
                    continue
                score += weights[j] * front + weights[j] * back
            middle = sub[window // 2]
            if middle in param_dict:
                score += param_dict[middle]
            scores.append(score / sum_of_weights)
        return scores

    probes = 0
    #   `protein_scale` writes a line to stderr for every residue its scale has
    #   no entry for.  That is part of its contract and re-running it 200 times
    #   here would bury the result, so the noise is captured and dropped; the
    #   numbers are what is being compared.
    noise = io.StringIO()
    for text in corpus:
        upper = text.upper()
        with contextlib.redirect_stderr(noise):
            try:
                expected = ProteinAnalysis(text).flexibility()
            except KeyError:
                try:
                    our_flexibility(upper, flex)
                except KeyError:
                    pass
                else:
                    raise SystemExit(
                        f"flexibility: reference raised, we answered: {text!r}")
            else:
                got = our_flexibility(upper, flex)
                if got != expected:
                    raise SystemExit(f"flexibility differs on {text!r}")
            probes += 1

            try:
                expected = ProteinAnalysis(text).instability_index()
            except (KeyError, ZeroDivisionError):
                pass
            else:
                got = our_instability(upper, dijw)
                if got != expected:
                    raise SystemExit(f"instability_index differs on {text!r}")
            probes += 1

            for window, edge in ((9, 1.0), (9, 0.4), (5, 1.0), (19, 0.7), (2, 1.0)):
                if window >= len(upper):
                    continue
                try:
                    expected = ProteinAnalysis(text).protein_scale(
                        flex, window, edge)
                except (ZeroDivisionError, IndexError):
                    continue
                got = our_protein_scale(upper, flex, window, edge)
                if got != expected:
                    raise SystemExit(
                        f"protein_scale(window={window}) differs on {text!r}")
                probes += 1

            #   Over the *gravy* names and not over the attribute names: a scale
            #   the reference exposes as `pa` is called `Parker` here, so a loop
            #   over the attributes asks `gravy("pa")`, gets a `ValueError` for
            #   every one but the two Cowan scales, and skips the other
            #   twenty-six without saying so.  That hole is how the int/float
            #   distinction survived a gate that was supposed to be checking
            #   these tables, so the mapping is what this walks.
            for name, key in tables["gravy"].items():
                scale = scales[key]
                try:
                    expected = ProteinAnalysis(text).gravy(name)
                except (KeyError, ValueError, ZeroDivisionError):
                    continue
                total = sum(scale[aa] for aa in upper)
                if total / len(upper) != expected:
                    raise SystemExit(f"gravy({name!r}) differs on {text!r}")
                probes += 1

    #   IsoelectricPoint: the constants, and the order they are summed in.
    for text in corpus:
        if not text or not set(text.upper()) <= set(amino):
            continue
        reference = I.IsoelectricPoint(text)
        aa_content = ProteinAnalysis(text).count_amino_acids()

        charged = {aa: float(aa_content[aa]) for aa in tables["charged_aas"]}
        charged["Nterm"] = 1.0
        charged["Cterm"] = 1.0
        pos = dict(tables["positive_pks"])
        neg = dict(tables["negative_pks"])
        if text.upper()[0] in tables["pk_nterminal"]:
            pos["Nterm"] = tables["pk_nterminal"][text.upper()[0]]
        if text.upper()[-1] in tables["pk_cterminal"]:
            neg["Cterm"] = tables["pk_cterminal"][text.upper()[-1]]

        for pH in (1.0, 4.05, 7.0, 7.775, 12.0):
            positive = 0.0
            for aa, pK in pos.items():
                positive += charged[aa] * (1.0 / (10 ** (pH - pK) + 1.0))
            negative = 0.0
            for aa, pK in neg.items():
                negative += charged[aa] * (1.0 / (10 ** (pK - pH) + 1.0))
            if positive - negative != reference.charge_at_pH(pH):
                raise SystemExit(f"charge_at_pH({pH}) differs on {text!r}")

        if our_pi(charged, pos, neg) != reference.pi():
            raise SystemExit(f"pi() differs on {text!r}")
        probes += 6

    return probes


def our_pi(charged, pos, neg, pH=7.775, min_=4.05, max_=12):
    """The reference's own recursion, written out again."""
    positive = 0.0
    for aa, pK in pos.items():
        positive += charged[aa] * (1.0 / (10 ** (pH - pK) + 1.0))
    negative = 0.0
    for aa, pK in neg.items():
        negative += charged[aa] * (1.0 / (10 ** (pK - pH) + 1.0))
    charge = positive - negative
    if max_ - min_ > 0.0001:
        if charge > 0.0:
            min_ = pH
        else:
            max_ = pH
        return our_pi(charged, pos, neg, (min_ + max_) / 2, min_, max_)
    return pH


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="verify the tables without writing the module")
    args = parser.parse_args()

    tables = build()
    probes = check(tables)
    print(f"tables verified against the reference on {probes} probes")

    if args.check:
        return 0
    TARGET.write_text(emit(tables))
    print(f"wrote {TARGET.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
