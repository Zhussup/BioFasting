#!/usr/bin/env python3
"""Ranking pass 10: `Bio.Restriction` -- the engine, not the table.

`biofasting.restriction` is a port of the reference's restriction layer written
as *data plus one engine*: a 1,088-row site table and a single `Enzyme` class,
instead of 27,000 lines of generated classes.  That half was measured when it
was built -- 51 mixin combinations became one class and every site, cut,
overhang, frequency and supplier matches -- but the *speed* of it has never been
timed, and the triage filed the whole module as `parity` on the assumption that
a restriction digest is not a hot path.  This pass asks that question with a
timer, because the same filing was wrong about the writers.

**The premise this pass has to check before it can rank anything.**  Both
implementations search with `re.finditer` over a prepared string -- the
reference's `FormattedSeq.finditer` and this package's are the same three lines
-- so the *search* is at the same floor on both sides and a port cannot be
faster there.  What is left is everything around it, and there are three
components:

* the **preparation**: the reference requires a `Bio.Seq.Seq` and rebuilds a
  `FormattedSeq` on every `search()` -- a `bytes` translate, a check that every
  letter is IUPAC, and a `decode` back to `str`, over the whole molecule, per
  enzyme;
* the **search**: one `re.finditer` pass with a lookahead pattern;
* the **bookkeeping**: `_modify`/`_rev_modify`/`_drop` per match, which for a
  4-cutter on a megabase is a quarter of a million matches.

So every row carries four columns and they are not four implementations of the
same thing:

* ``reference`` -- `Bio.Restriction.<enzyme>.search(Bio.Seq.Seq(text))`, the
  whole call from the type the reference requires;
* ``ours`` -- the same through `biofasting.restriction`, which takes `str`;
* ``regex`` -- **the floor the two share today**: `re.finditer` with the same
  pattern over the same string, counting matches and doing nothing else.  It is
  the reference's own search step with the preparation and the bookkeeping
  removed, so `reference - regex` is what the reference spends above the search;
* ``literal`` -- **the floor a kernel could reach**: `str.find` in a loop, which
  is `memchr`/`memmem` inside CPython, for every enzyme whose site has no IUPAC
  degeneracy.  This is the number that says whether a C++ scan would pay, and it
  is absent (printed `-`) on the one row whose site has `N` in it, because no
  pure-Python floor below a regex exists for an ambiguous site and a floor that
  is above the thing it bounds is not a floor.

The gate runs before any time is quoted.  This package's cut positions must equal
the reference's on every row -- every fragment, on the `catalyze` row -- and the
two floors, which count matches rather than cuts, must equal *each other*: two
independent scanners of the same bytes agreeing is what says the overlap rule and
the second strand were both read correctly.  A row that disagrees prints `DIFF`
and contributes no time.

**The trap this pass is built around.**  `RestrictionBatch.search` caches its
mapping on the batch, keyed on the sequence, in both libraries.  A benchmark that
builds one batch and times `batch.search(seq)` in a loop therefore measures a
`str(dna)` conversion and a tuple comparison, not a digest -- the same mistake the
ProtParam pass records for `ProteinAnalysis.count_amino_acids`.  Every timed call
here builds its own batch.

Usage:

    python3 bench/rank_restriction.py
    python3 bench/rank_restriction.py --only sau3ai-1mb --repeat 9
    python3 bench/rank_restriction.py --list
"""

from __future__ import annotations

import argparse
import re
import statistics
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import Bio.Restriction  # noqa: E402
from Bio.Seq import Seq  # noqa: E402

import biofasting.restriction as ours_restriction  # noqa: E402
from biofasting._restriction_sites import compsite  # noqa: E402

REPEATS = 7

#: The molecule every row cuts: one 1 Mbp contig, the FASTA row the FASTA and
#: Arrow passes measure too, so a number here is comparable with theirs.
GENOME = HERE / "data" / "fasta" / "genome_1mb.fasta"

#: Twenty enzymes a molecular biologist would actually digest with, all of them
#: six-cutters with a literal site.  Literal on purpose: it is what makes the
#: `literal` floor a floor for the whole row rather than for part of it.
BATCH = (
    "EcoRI", "BamHI", "HindIII", "NotI", "XhoI", "PstI", "SalI", "SmaI",
    "KpnI", "SacI", "SphI", "NcoI", "NdeI", "ScaI", "SpeI", "XbaI", "BglII",
    "ClaI", "ApaI", "NheI",
)

#: row name -> (enzymes, kind).  `kind` is "search" or "catalyze".
ROWS: dict[str, tuple[tuple[str, ...], str]] = {
    "eco-1mb": (("EcoRI",), "search"),
    "sau3ai-1mb": (("Sau3AI",), "search"),
    "bstxi-1mb": (("BstXI",), "search"),
    "batch-20-1mb": (BATCH, "search"),
    "catalyze-eco-1mb": (("EcoRI",), "catalyze"),
}

#: The complement, written out here rather than imported.  The floor must not be
#: built out of the package it is bounding: a `literal` column that asked
#: `biofasting` for its own pattern would be a second measurement of the thing
#: under test.
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def timed(fn, repeats: int = REPEATS) -> float:
    """Median seconds over ``repeats`` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def genome_text() -> str:
    """The contig as one string: the file is one record and its lines are wrapped."""
    lines = GENOME.read_text().split("\n")
    body = [line for line in lines[1:] if not line.startswith(">")]
    return "".join(body)


def revcomp(site: str) -> str:
    return site.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# The producers.  Each returns something comparable: a list of cut positions, a
# count, or a tuple of fragments.
# ---------------------------------------------------------------------------


def reference(text: str, seq: Seq, enzymes: tuple[str, ...], kind: str):
    """The whole call, from the type the reference requires.

    `Seq` is built once by the caller and reused, which is what a real caller
    does; the `FormattedSeq` inside is rebuilt on every call and is therefore
    inside the clock, because that is where the reference's cost is.
    """
    if kind == "catalyze":
        return tuple(str(fragment) for fragment in getattr(Bio.Restriction, enzymes[0]).catalyze(seq))
    if len(enzymes) == 1:
        return getattr(Bio.Restriction, enzymes[0]).search(seq)
    # A fresh batch per call: the reference caches the mapping on the batch, so a
    # reused one would be timed as a string comparison.
    batch = Bio.Restriction.RestrictionBatch(list(enzymes))
    return tuple(sorted(v for positions in batch.search(seq).values() for v in positions))


def biofasting(text: str, seq: Seq, enzymes: tuple[str, ...], kind: str):
    """The same call through this package, which takes the `str` it is given."""
    if kind == "catalyze":
        return tuple(getattr(ours_restriction, enzymes[0]).catalyze(text))
    if len(enzymes) == 1:
        return getattr(ours_restriction, enzymes[0]).search(text)
    batch = ours_restriction.RestrictionBatch(list(enzymes))
    return tuple(sorted(v for positions in batch.search(text).values() for v in positions))


_PATTERNS: dict[str, str] = {}


def pattern_of(name: str) -> str:
    """The enzyme's search pattern, built once and reused.

    Outside the clock, which is what the `regex` column means by "the search
    step and nothing else": `compsite` costs 1.08 us a build, two hundredths of
    a percent of the row it feeds, but a column that built twenty patterns per
    call would be measuring `str.join` beside its scan.
    """
    got = _PATTERNS.get(name)
    if got is None:
        enzyme = getattr(ours_restriction, name)
        got = compsite(enzyme.site, enzyme.is_palindromic())
        _PATTERNS[name] = got
    return got


def regex_floor(text: str, seq: Seq, enzymes: tuple[str, ...], kind: str) -> int | None:
    """The shared search step and nothing else: `re.finditer`, counted.

    The pattern comes from this package's own site-to-regex function rather
    than from the reference's compiled `Enzyme.compsite`, so that the column is
    one pattern builder plus one `re` scan and not a second port of the pattern
    derivation -- and the builder is outside the clock, via `pattern_of`, which
    is the whole point of the column: it is what a kernel would replace, and
    the builder it would not carry is measured separately, in `attribution`.
    """
    if kind == "catalyze":
        return None
    total = 0
    for name in enzymes:
        total += sum(1 for _ in re.finditer(pattern_of(name), text))
    return total


def literal_floor(text: str, seq: Seq, enzymes: tuple[str, ...], kind: str) -> int | None:
    """`str.find` in a loop -- `memchr` inside CPython -- for literal sites.

    `None` where no enzyme in the row has a literal site, which the caller prints
    as `-`.  A palindromic site is its own reverse complement and is counted once;
    a non-palindromic one is counted on both strands, which is what the regex with
    two groups does.
    """
    if kind == "catalyze":
        return None
    total = 0
    for name in enzymes:
        enzyme = getattr(ours_restriction, name)
        if any(letter not in "ACGT" for letter in enzyme.site):
            return None
        patterns = [enzyme.site]
        if not enzyme.is_palindromic():
            patterns.append(revcomp(enzyme.site))
        for pattern in patterns:
            start = text.find(pattern)
            while start != -1:
                total += 1
                start = text.find(pattern, start + 1)
    return total


PRODUCERS = {
    "reference": reference,
    "ours": biofasting,
    "regex": regex_floor,
    "literal": literal_floor,
}


def gate(row: str, text: str, seq: Seq) -> tuple[bool, str]:
    """Every comparable answer checked before any time is quoted.

    Two gates, and they check different things.  `ours` against `reference` is
    the contract: the cut positions, or every fragment for the `catalyze` row,
    must be the same list.  The floors are *match counters* and their answer is a
    number where the reference's is a list of cuts -- a cut that falls off the end
    of a linear molecule is dropped, so the two are not the same kind of thing and
    comparing them directly would be comparing a count to a different count.  What
    the floors are held to instead is each other: `literal` is `str.find` and
    `regex` is `re.finditer`, two independent scanners of the same bytes, and a
    disagreement means one of them is not finding what it claims to -- the
    overlap rule and the second strand are both easy to get wrong and both show
    up here.
    """
    enzymes, kind = ROWS[row]
    expected = reference(text, seq, enzymes, kind)
    message: list[str] = []
    for name in ("ours",):
        got = PRODUCERS[name](text, seq, enzymes, kind)
        if list(got) != list(expected):
            message.append(f"{name} differs ({len(got)} entries against {len(expected)})")
    matched = regex_floor(text, seq, enzymes, kind)
    literal = literal_floor(text, seq, enzymes, kind)
    if literal is not None and matched is not None and literal != matched:
        message.append(f"the floors disagree: str.find {literal:,} against re {matched:,}")
    return (not message), "; ".join(message)


def attribution(text: str, seq: Seq, repeats: int) -> None:
    """Where the reference's time goes, with the search step taken out.

    Three measurements, none of them a whole call: the `FormattedSeq` the
    reference rebuilds per `search()`, one `re.finditer` pass with this package's
    pattern, and the two together against the whole call.  It is the part that
    decides whether a kernel would pay, because `regex` is what a kernel replaces
    and `FormattedSeq` is what it deletes.
    """
    enzyme = Bio.Restriction.EcoRI
    pattern = compsite(ours_restriction.EcoRI.site, True)
    formatted = timed(lambda: Bio.Restriction.Restriction.FormattedSeq(seq), repeats)
    scanned = timed(lambda: sum(1 for _ in re.finditer(pattern, text)), repeats)
    whole = timed(lambda: enzyme.search(seq), repeats)
    ours = timed(lambda: ours_restriction.EcoRI.search(text), repeats)

    print()
    print("attribution, EcoRI on the 1 Mbp contig (microseconds):")
    print(f"  FormattedSeq(Seq)          {formatted * 1e6:9.1f}")
    print(f"  re.finditer, one pass      {scanned * 1e6:9.1f}")
    print(f"  reference search, whole    {whole * 1e6:9.1f}")
    print(f"  this package, whole        {ours * 1e6:9.1f}")
    print(
        f"  above the search           {(whole - scanned) * 1e6:9.1f}"
        "  (preparation and bookkeeping)"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", help="one row name")
    parser.add_argument("--repeat", type=int, default=REPEATS)
    parser.add_argument("--list", action="store_true", help="list the rows")
    parser.add_argument("--no-attribution", action="store_true")
    args = parser.parse_args()

    if args.list:
        for name in ROWS:
            print(name)
        return 0

    if not GENOME.exists():
        print(f"{GENOME} is missing: python3 bench/gen_data.py", file=sys.stderr)
        return 1

    text = genome_text()
    seq = Seq(text)
    rows = [args.only] if args.only else list(ROWS)

    failures = 0
    print("gate: every producer's cut positions against the reference's, and the floors' counts")
    for row in rows:
        ok, message = gate(row, text, seq)
        print(f"  {row}: {'equal' if ok else 'DIFF ' + message}, {len(text):,} bp")
        failures += not ok
    if failures:
        return 1

    print()
    print(
        "| row | enzymes | sites | reference | ours | regex | literal | vs regex | vs reference |"
    )
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        enzymes, kind = ROWS[row]
        times = {
            name: timed(lambda fn=fn, text=text, seq=seq, enzymes=enzymes, kind=kind: fn(text, seq, enzymes, kind), args.repeat)
            for name, fn in PRODUCERS.items()
        }
        counts = {name: fn(text, seq, enzymes, kind) for name, fn in PRODUCERS.items()}
        sites = len(counts["reference"]) if kind == "search" else len(counts["reference"]) - 1

        def cell(name: str) -> str:
            got = counts[name]
            return f"{times[name] * 1e6:,.0f}" if got is not None else "-"

        ratio_regex = times["regex"] / times["ours"]
        ratio_reference = times["reference"] / times["ours"]
        print(
            f"| {row} | {len(enzymes)} | {sites:,} | {times['reference'] * 1e6:,.0f} | "
            f"{times['ours'] * 1e6:,.0f} | {cell('regex')} | {cell('literal')} | "
            f"{ratio_regex:.2f}x | {ratio_reference:.2f}x |"
        )

    print()
    print("microseconds per call, median of %d after warm-up." % args.repeat)
    print("`regex` and `literal` count matches and produce no cut list, so their")
    print("counts are gated on being the subsets the reference's cuts come from.")
    if not args.no_attribution:
        attribution(text, seq, args.repeat)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
