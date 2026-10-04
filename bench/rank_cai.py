#!/usr/bin/env python3
"""Ranking pass 11: `CodonAdaptationIndex.calculate` -- what a CAI call is made of.

The runner's `cai-calculate` row says this package's `calculate` is 1.0x the
reference's, and that is the reason this pass does not start where a ranking
pass usually does: there is no competitor and there is no faster Python
implementation in hand.  Both `calculate` functions are literal translations of
the same Biopython loop, and both pay the same per-codon costs that loop
carries -- a string slice per codon, a membership test against a two-element
list, a dict lookup, one `log()` call, and the length bookkeeping.  What this
pass asks with a timer is what those costs are worth, and what would be left
for a kernel to delete.

Two floors answer that, and they are *derived* floors rather than competitors:
the same CAI value, computed by the same rules, with the loop stripped of
everything a C++ scan would not need.

* ``slice-walk`` -- the loop's own shape, stripped: no `SeqRecord` branch, no
  `upper()`, no membership test against the two-element lists, no try/except
  around the dict.  The skip rules are folded into a table of precomputed
  `log(w)` values that carries `None` for the two skipped codons (ATG and TGG,
  which the reference skips *before* it looks them up, and whose skipping
  changes the length divisor -- so the table cannot just carry `log(1)`
  there); a codon absent from the table is one no index would accept, and no
  substrate here carries one (see the traps).
* ``tuple-walk`` -- the same table keyed by three-character tuples instead of
  strings, so the walk builds no per-codon string at all: `zip` over the three
  strides yields one tuple per codon directly.  It is the shape a kernel
  replaces the whole walk with, minus the C.

Both floors are held to the reference's **value**: they return the CAI, not a
count, so the gate compares floats -- and because every producer walks the
same codons in the same order and adds the same `log(w)` terms in the same
order, the sums are bitwise equal, a difference there is a rule difference,
and no float tolerance is needed to say so.

**The two traps this pass avoids.**  (1) *The index is not in the clock.*  The
runner's row says this directly, and so does the reference's documentation:
`calculate` is what a caller spends, the table is built once.  Every index
here is built before any clock starts; `index-1m` is the one row that times
the build, because the build is its subject.  (2) *N is not a substrate for
CAI.*  CAI is defined on unambiguous coding sequences, and both libraries
refuse a sequence carrying N -- the reference raises `TypeError` out of
`self[codon]`, and this package raises the same raise on the same lookup,
pinned by
`tests/test_sequtils.py::test_codon_adaptation_index_calculate_rejects_an_illegal_codon`.
The 1 Mbp contig carries 1,000 N the corpus generator plants, so the contig
rows read the contig with those N replaced by A: length and codon count are
untouched and only the letters a lookup sees change.

There is no `count_kmers` floor here, although the kernel has a codon-sized
scanner that walks these bytes in C++: it counts *overlapping* k-mers with no
stride, and CAI reads codons with a stride of three, so it counts a different
thing.  A CAI in the kernel would be new code, and this pass measures whether
new code is worth.

Usage:

    python3 bench/rank_cai.py
    python3 bench/rank_cai.py --only calc-1m --repeat 9
    python3 bench/rank_cai.py --list
"""

from __future__ import annotations

import argparse
import sys
from math import exp, log
from pathlib import Path
from typing import Callable

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from Bio.SeqIO import parse  # noqa: E402
from Bio.SeqUtils import CodonAdaptationIndex as ReferenceCAI  # noqa: E402

import biofasting  # noqa: E402
from rank_restriction import genome_text, timed  # noqa: E402

REPEATS = 7

#: The genes: the corpus's 10,000 FASTQ reads as 150-base coding sequences,
#: 50 in-frame codons each, none of them carrying N.
GENES = HERE / "data" / "fastq" / "reads_10k.fastq"

#: The contig the mass row scores: one 1 Mbp sequence, 333,334 codons.
CONTIG = HERE / "data" / "fasta" / "genome_1mb.fasta"

#: row -> (path, kind).  `kind` is "calc" (the CAI of the substrate) or
#: "index" (the table built from it, the pass's only build row).
ROWS: dict[str, tuple[Path, str]] = {
    "calc-genes-10k": (GENES, "calc"),
    "calc-1m": (CONTIG, "calc"),
    "index-1m": (CONTIG, "index"),
}

#: `log(w)` per codon, and `None` for the two the reference skips.  Filled once
#: by `build_weight_tables` during the gate phase, from this package's index --
#: which the gate has then just checked is the reference's index as a dict.
LOG_WEIGHTS: dict[str, float | None] = {}

#: The same, keyed by the three characters instead of the joined string, so the
#: walk needs no per-codon string.  Filled with the first.
TUPLE_WEIGHTS: dict[tuple[str, str, str], float | None] = {}

#: The reference skips these before it looks them up -- and the skip is not
#: merely the value, it takes the codon out of the length divisor too, which is
#: why the table cannot carry `log(1)` in their place.
SKIPPED = ("ATG", "TGG")


def contig_text() -> str:
    """The contig as one string, N replaced by A, cut to a whole codon.

    The cut is not cosmetic: CAI is defined on coding sequences, the contig is
    1,000,000 bases -- one base more than a whole number of codons -- and the
    reference's own builder refuses that tail with `ValueError: illegal codon
    'T'`, the same refuse this package's does.  Both raise before the gate can
    compare anything, so the gate caught the first draft of this row rather
    than a bench number.
    """
    contig = genome_text().replace("N", "A")
    return contig[: len(contig) - len(contig) % 3]


def gene_texts() -> list[str]:
    """The 10,000 reads as strings, exactly the substrate the runner uses."""
    return [str(record.seq) for record in parse(GENES, "fastq")]


def build_weight_tables(index: dict) -> None:
    """Fill `LOG_WEIGHTS` and `TUPLE_WEIGHTS` from an index, once.

    Called in the gate phase, after the index equality has been checked, so
    the floors read the same weights whichever library produced them.  Every
    tabled weight is finite: both builders replace a zero count with 0.5, so
    `log` is defined for every entry, stops included.  A stop codon is in the
    index (67 entries -- the 64 codons plus the 3 stops, all with weights set),
    it is an ordinary lookup with a real `log(w)`, and `calculate` skips a stop
    only *when the index has none* -- not a case this corpus reaches.
    """
    LOG_WEIGHTS.clear()
    TUPLE_WEIGHTS.clear()
    for codon, weight in index.items():
        LOG_WEIGHTS[codon] = log(weight)
    for codon in SKIPPED:
        if codon in LOG_WEIGHTS:
            LOG_WEIGHTS[codon] = None
    for codon, value in LOG_WEIGHTS.items():
        TUPLE_WEIGHTS[(codon[0], codon[1], codon[2])] = value


# ---------------------------------------------------------------------------
# The calc producers.  Each receives an index built outside every clock; the
# two implementations differ only in which library's index they were handed.
# Both return one CAI per item in the substrate, because the gate compares
# against the reference's list.
# ---------------------------------------------------------------------------


def reference(sequences, index) -> list[float]:
    """`Bio.SeqUtils.CodonAdaptationIndex.calculate`, per item."""
    return [index.calculate(one) for one in sequences]


def ours(sequences, index) -> list[float]:
    """`biofasting.CodonAdaptationIndex.calculate`, the same shape."""
    return [index.calculate(one) for one in sequences]


#: The reference's loop with everything a kernel would delete already deleted,
#: in Python.
def slice_walk(sequences, index) -> list[float]:
    return [cai_by_stride(one) for one in sequences]


def cai_by_stride(sequence) -> float:
    total, length, weights = 0.0, 0, LOG_WEIGHTS
    for at in range(0, len(sequence), 3):
        value = weights[sequence[at : at + 3]]
        if value is not None:
            total += value
            length += 1
    return exp(total / length)


#: `zip` over the three strides, so no per-codon string is built.
def tuple_walk(sequences, index) -> list[float]:
    weights = TUPLE_WEIGHTS
    results = []
    for one in sequences:
        total, length = 0.0, 0
        first, second, third = one[0::3], one[1::3], one[2::3]
        for triplet in zip(first, second, third):
            value = weights[triplet]
            if value is not None:
                total += value
                length += 1
        results.append(exp(total / length))
    return results


CALC: dict[str, Callable] = {
    "reference": reference,
    "ours": ours,
    "slice-walk": slice_walk,
    "tuple-walk": tuple_walk,
}

#: The build producers: one table, from one substrate.  No floors -- the
#: question this row answers is our port against the reference's, which the
#: runner has so far asserted but never timed.
INDEX: dict[str, Callable] = {
    "reference": lambda sequence, _index: ReferenceCAI([sequence]),
    "ours": lambda sequence, _index: biofasting.CodonAdaptationIndex([sequence]),
}


def substrate_for(row: str):
    """The substrate a row scores.

    The calc rows score *items* and so receive a list of sequences; the build
    row builds a table from *one* sequence and so receives the string itself --
    which is the shape `ReferenceCAI` takes from every real caller.
    """
    if row == "calc-genes-10k":
        return gene_texts()
    if ROWS[row][1] == "calc":
        return [contig_text()]
    return contig_text()


def indexes_for(sequences):
    """The two indexes, built outside every clock, plus their equality.

    Both constructors take an *iterable of sequences* -- that is their
    signature on both sides, and a bare string is iterated a character at
    a time, which raises `illegal codon` before the gate can compare
    anything.  The build row's gate and its timing both pass `[[str]]` -- a
    list of one.
    """
    ref_index = ReferenceCAI(sequences)
    our_index = biofasting.CodonAdaptationIndex(sequences)
    equal = dict(ref_index) == dict(our_index)
    return ref_index, our_index, equal


def gate(row: str) -> tuple[bool, str]:
    """Every producer's answer against the reference's, before any clock.

    For the calc rows the comparison is the CAI value itself, exactly: every
    producer walks the same codons in the same order, adds the same `log(w)`
    terms in the same order, and `exp(total / length)` is therefore bitwise
    the reference's answer.  A difference here is a rule difference -- a
    skipped codon the reference counted, a weight read from the wrong table --
    and no float tolerance is needed to say so.

    The build row has one gate, the table itself: the two libraries' indexes
    must be equal as dicts, because a floor built on a different table would
    time a different thing.
    """
    row_path, kind = ROWS[row]
    substrate = substrate_for(row)
    ref_index, our_index, indexes_equal = indexes_for(
        substrate if kind == "calc" else [substrate]
    )
    message: list[str] = []
    if not indexes_equal:
        message.append("the two indexes disagree as dicts")
    # The tables are rebuilt every row, never reused: the weights come from the
    # substrate of this row, and the first draft filled them once and read them
    # on `calc-1m` from the reads-10k table -- which the gate caught as floors
    # a third of a percent away from the reference on the contig row, while
    # agreeing on the row whose index had built them.
    build_weight_tables(our_index)
    if kind == "calc" and indexes_equal:
        expected = reference(substrate, ref_index)
        for name, fn in CALC.items():
            index = ref_index if name == "reference" else our_index
            got = fn(substrate, index)
            if list(got) != list(expected):
                at = next(
                    (i for i, (one, other) in enumerate(zip(got, expected)) if one != other),
                    0,
                )
                message.append(
                    f"{name} differs at item {at}: {got[at]!r} vs {expected[at]!r}"
                )
    elif kind == "index":
        for name, fn in INDEX.items():
            if name == "reference":
                continue
            got = fn(substrate, ref_index)
            if dict(got) != dict(ref_index):
                message.append(f"{name} builds a different table")
    return (not message), "; ".join(message)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", help="one row name")
    parser.add_argument("--repeat", type=int, default=REPEATS)
    parser.add_argument("--list", action="store_true", help="list the rows")
    args = parser.parse_args()

    if args.list:
        for name in ROWS:
            print(name)
        return 0

    for path, _kind in ROWS.values():
        if not path.exists():
            print(f"{path} is missing: python3 bench/gen_data.py", file=sys.stderr)
            return 1

    rows = [args.only] if args.only else list(ROWS)

    failures = 0
    print("gate: every producer's answer against the reference's, exactly compared")
    for row in rows:
        ok, message = gate(row)
        print(f"  {row}: {'equal' if ok else 'DIFF ' + message}")
        failures += not ok
    if failures:
        return 1

    # --- the table ---------------------------------------------------------
    units: dict[str, str] = {}
    print()
    print(
        "| row | items | codons | reference | ours | slice-walk | tuple-walk "
        "| ours vs reference |"
    )
    print("|---|---|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        path, kind = ROWS[row]
        substrate = substrate_for(row)
        ref_index, our_index, _equal = indexes_for(
            substrate if kind == "calc" else [substrate]
        )
        # Floors read this row's weights, not the previous row's -- see gate.
        build_weight_tables(our_index)
        indexes = {"reference": ref_index, "ours": our_index}
        producers = CALC if kind == "calc" else INDEX
        if kind == "calc":
            items, codons = len(substrate), sum(len(one) // 3 for one in substrate)
        else:
            items, codons = 1, len(substrate) // 3
        times = {
            name: timed(
                # Floors take an index they do not read: that argument is the
                # shape every calc producer has, not a claim about use.
                lambda fn=fn, one=substrate, index=indexes.get(name, ref_index): fn(one, index),
                args.repeat,
            )
            for name, fn in producers.items()
        }
        # One unit for both calc rows -- `us a codon` -- so they read against
        # each other; the build row is per table, in the same unit column.
        divisor, units[row] = (codons, "us a codon") if kind == "calc" else (1, "us a build")

        def rate(seconds: float) -> str:
            return f"{seconds / divisor * 1e6:.3f}"

        cell = lambda name: rate(times[name]) if name in times else "-"
        vs_reference = times["reference"] / times["ours"]
        print(
            f"| {row} | {items:,} | {codons:,} | {cell('reference')} | {cell('ours')} "
            f"| {cell('slice-walk')} | {cell('tuple-walk')} | {vs_reference:.2f}x |"
        )

    print()
    print(f"median of {args.repeat} after warm-up, in the unit the table lists")
    print("per row, both held to their reference by the gate above -- exactly")
    print("compared, no tolerance.")
    for row, counted in units.items():
        print(f"  {row}: {counted}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())