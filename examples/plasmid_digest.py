#!/usr/bin/env python3
"""Plasmid digest: a batch of six-cutters against a small circular DNA.

What a molecular biologist asks before cloning: which of my usual enzymes cut
this plasmid, where, and how long are the pieces?  The example answers with
five six-cutters and one eight-cutter -- EcoRI, BamHI, HindIII, XhoI, PstI and
NotI -- and prints, per enzyme, the site count and the cut positions the
`search` contract pins, then the fragment lengths the `catalyze` contract pins.
Both come from `biofasting.restriction`, which ports the reference's engine
whole: positions are 1-based, a fragment list for a circular molecule sums to
the ring, and an enzyme that does not cut leaves the ring whole.

The default substrate is *not* a plasmid -- it is a 5,000 bp window of the
corpus contig, read as a ring so that the circular arithmetic has something
to wrap around.  A real plasmid goes on the command line:

    python3 examples/plasmid_digest.py --path pBR322.fasta
    python3 examples/plasmid_digest.py --linear      # a PCR product instead

Usage:

    python3 examples/plasmid_digest.py [--path FILE] [--window N] [--linear]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT = HERE.parent / "bench" / "data" / "fasta" / "genome_1mb.fasta"

#: Six enzymes a cloning notebook names first, five of them six-cutters and
#: NotI an eight-cutter, every one with a literal site.  A batch is spelled by
#: name, and `elements()` returns them sorted.
CUTTERS = (
    "EcoRI",
    "BamHI",
    "HindIII",
    "XhoI",
    "PstI",
    "NotI",
)

#: Where the demo ring starts in the contig -- a window picked once by looking
#: at the cut table, so the demo has five enzymes on it instead of the zero or
#: one a random offset offers; the contig is random sequence, so the choice
#: changes the numbers, not the lesson.
DEMO_OFFSET = 100_000


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--path", type=Path, default=DEFAULT, help="a FASTA plasmid to digest"
    )
    #   The demo window, for the default substrate only: the contig is a
    #   million bases and a plasmid workflow wants thousands.
    parser.add_argument(
        "--window", type=int, default=10_000, help="demo ring size (bp)"
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="treat the molecule as linear, not a ring",
    )
    args = parser.parse_args()

    import biofasting
    import biofasting.restriction as restriction

    if args.path is DEFAULT and not args.path.exists():
        print(f"{args.path} is missing: python3 bench/gen_data.py", file=sys.stderr)
        return 1

    if args.path is DEFAULT:
        #   The demo ring: the first record's opening window, labelled as what
        #   it is -- a corpus window standing in for a plasmid so the circular
        #   arithmetic has something to wrap.
        index = biofasting.open_fasta(args.path)
        name = next(iter(index))
        title = index.title(name)
        dna = index[name][DEMO_OFFSET : DEMO_OFFSET + args.window]
        dna = dna.upper().decode("ascii", "replace")
        label = f"demo-ring of {name} (a corpus window, not a plasmid)"
    else:
        #   A real substrate: one record, the first, its name and title as the
        #   reader splits them.
        records = biofasting.read_fasta(args.path)
        if not records:
            print(f"{args.path} holds no records", file=sys.stderr)
            return 1
        name, title, dna = records[0]
        #   The sequence is uppercased but never edited: an N a caller's
        #   plasmid carries stays, and simply never matches a site.
        dna = (dna.decode("ascii", "replace") if isinstance(dna, bytes) else dna).upper()
        label = f"{name} ({title or 'no title'})"

    batch = restriction.RestrictionBatch(list(CUTTERS))
    sites = batch.search(dna, linear=args.linear)
    length = len(dna)

    print(f"digesting {label}, {length:,} bp, {'linear' if args.linear else 'circular'}")
    print(f"batch: {batch.__str__()}\n")

    #   `search` keys on the Enzyme objects in the batch, sorted by name.
    for enzyme in sorted(sites):
        positions = sites[enzyme]
        site_report = f"  {enzyme.name:8s} {len(positions)} site(s), site {enzyme.elucidate()}"
        if not positions:
            print(site_report)
            continue
        cut_list = ", ".join(str(p) for p in positions[:8])
        if len(positions) > 8:
            cut_list += ", ..."
        print(site_report)
        print(f"           cuts at {cut_list}")
        fragments = enzyme.catalyze(dna, linear=args.linear)
        lengths = [len(f) for f in fragments]
        if len(lengths) > 4:
            shown = f"{lengths[:4]} + {len(lengths) - 4} more"
        else:
            shown = str(lengths)
        print(f"           fragments ({sum(lengths):,} bp): {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())