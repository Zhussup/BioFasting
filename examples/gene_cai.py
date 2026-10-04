#!/usr/bin/env python3
"""Codon adaptation: build an index from a gene set, then score every gene.

CAI answers one comparative question -- how closely does this coding sequence
match codon usage of *this* set of genes?  -- so the pipeline has two steps
and they use the package's two-part `CodonAdaptationIndex`: the weights table
is built once from the whole gene set, then every member is scored against
the table the others produced.  The example builds the table from every read
in the quick FASTQ corpus (reads as 150-base "genes"), scores them all, and
prints the mean with the top and the bottom of the ranking.  On random
sequence the answer sits mid-range and flat, as it should -- CAI rewards
specialisation, and a random corpus has none to reward.

A bare string is not a substrate: the constructor takes an *iterable of
sequences*, because it counts codons across all of them.  Feeding one string
feeds it one base at a time, and the refusal is loud.

    python3 examples/gene_cai.py
    python3 examples/gene_cai.py --path mygenes.fa

Usage:

    python3 examples/gene_cai.py [--path FILE] [--top N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT = HERE.parent / "bench" / "data" / "fastq" / "reads_10k.fastq"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--path", type=Path, default=DEFAULT, help="the gene set")
    parser.add_argument("--top", type=int, default=5, help="how many of each end to show")
    args = parser.parse_args()
    if not args.path.exists():
        print(f"{args.path} is missing: python3 bench/gen_data.py", file=sys.stderr)
        return 1

    import biofasting

    #   The substrates as strings, in file order, so rank positions name the
    #   reads in the file they came from.  Ten thousand 150-base reads is a
    #   list of 1.5 MB -- the size read_fasta is documented against, and the
    #   size `CodonAdaptationIndex` walks twice anyway.
    sequences = [sequence.decode("ascii", "replace").upper()
                 for _title, sequence, _quality in biofasting.open_fastq(args.path)]
    if not sequences:
        print(f"{args.path} holds no records", file=sys.stderr)
        return 1

    index = biofasting.CodonAdaptationIndex(sequences)
    scored = [(index.calculate(one), number) for number, one in enumerate(sequences)]
    scored.sort(reverse=True)

    mean = sum(value for value, _ in scored) / len(scored)
    print(f"index built from {len(sequences):,} genes "
          f"({sum(len(one) // 3 for one in sequences):,} codons, "
          f"{len(index)} codons weighted)")
    print(f"mean CAI {mean:.4f}, range {scored[-1][0]:.4f} .. {scored[0][0]:.4f}\n")
    top = scored[: args.top]
    bottom = scored[-args.top :][::-1]
    print("the most adapted genes:")
    for value, number in top:
        print(f"  read {number + 1:6,}  CAI {value:.4f}")
    print("the least adapted:")
    for value, number in bottom:
        print(f"  read {number + 1:6,}  CAI {value:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())