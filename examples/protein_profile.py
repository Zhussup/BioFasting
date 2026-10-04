#!/usr/bin/env python3
"""Protein profile: the numbers a protein record is summarised with.

`ProteinAnalysis` is a calculator, not a reader -- it takes a sequence and
answers: weight, aromaticity, GRAVY, isoelectric point, charge at a given
pH, the secondary-structure fractions, the extinction coefficients, the
instability index, the flexibility profile.  The example runs the whole
calculator over one protein and prints the table, then does what the
numbers are for: flags the instability verdict.

The demo sequence is the protein from the reference's own class
documentation (the fragment the HSP90 summary in `Bio.SeqUtils.ProtParam`
uses) -- it is here verbatim so the example and the docstrings agree on
every printed number.  `--path` reads a one-letter-code protein from a
FASTA file instead; windows are then sliced off it, and each window is
labelled as a window, because an instability index over one exon of a
protein is a statement about that exon.

    python3 examples/protein_profile.py
    python3 examples/protein_profile.py --path myprotein.fa

Usage:

    python3 examples/protein_profile.py [--path FILE] [--window N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).parent
#: The protein the reference's `ProtParam` documents; same literal the tests
#: carry, so an example that is also documentation drifts in the daylight.
DEMO_PROTEIN = (
    "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT"
    "RDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEEC"
    "LFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILF"
    "LPLPV"
)


def profile(analyzed: str, name: str, protein) -> None:
    """One protein's table, printed."""
    print(f"-- {name}: {len(analyzed)} aa")
    print(f"   molecular weight   {protein.molecular_weight():>10,.1f} Da")
    arom = protein.aromaticity()
    print(f"   aromaticity        {arom:>10.4f}")
    print(f"   GRAVY              {protein.gravy():>10.4f}")
    print(f"   isoelectric point  {protein.isoelectric_point():>10.2f}")
    print(f"   charge at pH 7.0   {protein.charge_at_pH(7.0):>10.2f}")
    helix, turn, sheet = protein.secondary_structure_fraction()
    print(f"   helix/turn/sheet   {helix:.2f} / {turn:.2f} / {sheet:.2f}")
    reduced, oxidised = protein.molar_extinction_coefficient()
    print(f"   extinction coeff    {reduced:,} / {oxidised:,} (reduced / oxidised)")
    instability = protein.instability_index()
    verdict = "stable" if instability < 40.0 else "unstable"
    print(f"   instability index  {instability:>10.2f}  ({verdict})")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--path", type=Path, default=None, help="a FASTA protein to profile"
    )
    parser.add_argument(
        "--window", type=int, default=100, help="window size for --path profiles"
    )
    args = parser.parse_args()

    import biofasting

    analyzed, name, window_size = (DEMO_PROTEIN, "the reference's doc protein", None)
    if args.path is not None:
        if not args.path.exists():
            print(f"{args.path} is missing", file=sys.stderr)
            return 1
        records = biofasting.read_fasta(args.path)
        if not records:
            print(f"{args.path} holds no records", file=sys.stderr)
            return 1
        #   The first record, profiled whole and then in windows: every
        #   window labelled for what it is.
        name, title, sequence = records[0]
        analyzed = (sequence.decode("ascii", "replace") if isinstance(sequence, bytes)
                    else sequence).upper()
        shortname = name
        name = f"{name} ({title or 'no title'})"
        window_size = args.window

    protein = biofasting.ProteinAnalysis(analyzed)
    profile(analyzed, name, protein)
    #   A protein shorter than the window is profiled whole only; a longer
    #   one gets its windows too, and the last window may be the short one
    #   rather than dropped.
    if window_size is not None and len(analyzed) > window_size:
        for number, start in enumerate(range(0, len(analyzed), window_size), start=1):
            piece = analyzed[start : start + window_size]
            profile(piece, f"window {number} of {shortname} (a window, not the protein)",
                    biofasting.ProteinAnalysis(piece))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())