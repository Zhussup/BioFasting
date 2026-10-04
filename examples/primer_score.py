#!/usr/bin/env python3
"""Primer scoring: candidates against a template, ranked by local alignment.

The question a PCR bench asks: of these candidate primers, which one binds
this template, where, and does the one with errors still bind?  The example
answers with `biofasting.Aligner` in local mode -- local is the right mode
for a probe, because a primer is a needle in a longer template and global
alignment would drag the whole template into the score.

Five candidates are cut from the template itself (the contig is random
letters, so a 20-mer matches exactly once), and one decoy carries three
deliberate substitutions: the ranking is what says a primer with a few wrong
bases still binds, but worse.  Every result is printed with the alignment
itself, the aligned columns and the template span -- the coordinates the
`Alignment` contract pins.

The prices are the trade every primer designer makes explicit:

* match ``2``, mismatch ``-3`` -- a clean 20-mer scores 40; with three
  errors it scores 25, and that drop is the specificity answer;
* gap open ``-5``, extend ``-2`` -- an indel costs more than a substitution
  but a single one is survivable.

    python3 examples/primer_score.py
    python3 examples/primer_score.py --match 3 --mismatch -2

Usage:

    python3 examples/primer_score.py [--path FILE] [--offset N] [--size N]
                                     [--match M] [--mismatch M]
                                     [--open-gap P] [--extend-gap P] [--list]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT = HERE.parent / "bench" / "data" / "fasta" / "genome_1mb.fasta"

#: Where the template starts in the contig -- a window picked once, and
#: labelled for what it is: a random stretch of the corpus, not a real genome.
DEMO_OFFSET = 250_000
DEMO_SIZE = 1_000

#: Where the candidates come from within the template.  The decoy is the
#: *last* start, rebuilt with three substitutions -- not a different
#: sequence, the same one mutated, because that is the question the ranking
#: answers.
PRIMER_STARTS = (50, 200, 400, 700, 950)
PRIMER_LENGTH = 20

#: The substitution rule for the decoy: every fourth base replaced by the
#: next letter of the alphabet -- spaced, so the mutations do not pile up on
#: one end of the primer.  Three substitutions on a 20-mer.
DECOY_STRIDE = 4
DECOY_COUNT = 3

#: The substitution rule's alphabet: replace each drawn base with the next
#: letter in this order, cycling.
ALPHABET = "ACGT"


def decoy(letters: str) -> str:
    """`letters` with every `DECOY_STRIDE`-th base (first `DECOY_COUNT`) swapped.

    Spaced substitutions, cycling through the alphabet: the decoy stays a
    20-mer, it is only wrong at known places, and those places are printable,
    which is what makes the score drop legible.
    """
    swapped = list(letters)
    for number in range(DECOY_COUNT):
        position = number * DECOY_STRIDE
        swapped[position] = ALPHABET[
            (ALPHABET.index(swapped[position]) + 1) % len(ALPHABET)
        ]
    return "".join(swapped)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--path", type=Path, default=DEFAULT, help="the FASTA template")
    parser.add_argument("--offset", type=int, default=DEMO_OFFSET, help="template start")
    parser.add_argument("--size", type=int, default=DEMO_SIZE, help="template length")
    parser.add_argument("--match", type=float, default=2.0, help="match score")
    parser.add_argument("--mismatch", type=float, default=-3.0, help="mismatch score")
    parser.add_argument("--open-gap", type=float, default=-5.0, help="gap open")
    parser.add_argument("--extend-gap", type=float, default=-2.0, help="gap extend")
    args = parser.parse_args()
    if args.path is DEFAULT and not args.path.exists():
        print(f"{args.path} is missing: python3 bench/gen_data.py", file=sys.stderr)
        return 1

    import biofasting

    index = biofasting.open_fasta(args.path)
    name = next(iter(index))
    record = index[name][args.offset : args.offset + args.size]
    template = (record.decode("ascii", "replace") if isinstance(record, bytes)
                else record).upper()

    candidates = []
    for start in PRIMER_STARTS:
        piece = template[start : start + PRIMER_LENGTH]
        candidates.append(("clean", piece, f"template[{start}:{start + PRIMER_LENGTH}]"))
    last_start = PRIMER_STARTS[-1]
    candidates.append(
        (
            "decoy",
            decoy(template[last_start : last_start + PRIMER_LENGTH]),
            f"template[{last_start}:{last_start + PRIMER_LENGTH}] "
            f"plus {DECOY_COUNT} substitutions",
        )
    )

    aligner = biofasting.Aligner(
        mode="local",
        match_score=args.match,
        mismatch_score=args.mismatch,
        open_gap_score=args.open_gap,
        extend_gap_score=args.extend_gap,
    )
    clean_score = PRIMER_LENGTH * args.match

    print(
        f"template: {len(template):,} bp of {name} at offset {args.offset:,} "
        "(a corpus window, not a real genome)"
    )
    print(
        f"prices: match {args.match}, mismatch {args.mismatch}, "
        f"open {args.open_gap}, extend {args.extend_gap}; "
        f"a clean {PRIMER_LENGTH}-mer scores {clean_score:g}\n"
    )

    #   Score once, rank later: every candidate's local alignment against the
    #   full template, held with the alignment so it can be printed.
    scored = [
        (aligner.align(piece, template), kind, label)
        for kind, piece, label in candidates
    ]
    order = sorted(range(len(scored)), key=lambda at: (-scored[at][0].score, at))
    print(f"rank {'score':>6}  candidate")
    for rank, at in enumerate(order, start=1):
        alignment, kind, label = scored[at]
        print(f"   {rank} {alignment.score:6g}  {kind}: {label}")

    #   The winner's own columns: the aligned pair printed as two lines, and
    #   the template span the coordinates carry.  This is what 'it binds'
    #   means, shown rather than claimed -- and the decoy's pair is printed
    #   beside it, because the substitutions are the lesson: they sit in the
    #   columns visibly, they cost the score, and they did not stop the
    #   binding.
    for which, header in ((order[0], "the winner"), (order[-1], "the decoy")):
        alignment, kind, label = scored[which]
        print(f"\n{header} ({kind}, {label}) binds as:")
        print(f"  primer    {alignment[0]}")
        print(f"  template  {alignment[1]}")
        span = alignment.coordinates
        print(f"  template span {span[1][0]:,}..{span[1][-1]:,}, score {alignment.score:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())