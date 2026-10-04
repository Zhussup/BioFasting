#!/usr/bin/env python3
"""A FASTQ QC pipeline: read, measure, filter, write, report.

The pipeline the readers were built for, end to end and nothing more:
`biofasting.open_fastq` streams the file (gzip is sniffed from the magic
bytes, not the file name), every record is measured with the package's own
`gc_fraction`, the reads that pass every threshold are written once with
`write_fastq`, and the summary the gate asserts is printed.

The thresholds are knobs, not opinions:

    python3 examples/fastq_qc.py
    python3 examples/fastq_qc.py --min-q 30 --gc-low 35

The default substrate is the 10,000-read quick corpus this repository builds
(`python3 bench/gen_data.py --quick` reports the one command when it is not).

Usage:

    python3 examples/fastq_qc.py [--path FILE] [--min-q N] [--gc-low P] [--gc-high P] [--out FILE]
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT = HERE.parent / "bench" / "data" / "fastq" / "reads_10k.fastq"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--path", type=Path, default=DEFAULT, help="the FASTQ to read")
    #   Defaults tuned so the pipeline actually filters on this corpus: its
    #   quality ladder decays 34.8 -> 32.2 -> 27.6, so 32 trims the worst zone
    #   and 42/52 trims the GC tails -- both chosen by looking at the
    #   distribution once, and both movable on the command line.
    parser.add_argument("--min-q", type=float, default=32.0, help="minimum mean quality")
    parser.add_argument("--gc-low", type=float, default=42.0, help="minimum GC %%")
    parser.add_argument("--gc-high", type=float, default=52.0, help="maximum GC %%")
    #   The kept reads land in a fresh temp directory by default, because a
    #   demo that writes into its own repository is a demo that makes every
    #   `git status` dirty; `--out` points it anywhere else.
    parser.add_argument("--out", type=Path, default=None, help="the kept reads")
    args = parser.parse_args()
    if not args.path.exists():
        print(f"{args.path} is missing: python3 bench/gen_data.py", file=sys.stderr)
        return 1
    if args.out is None:
        args.out = Path(tempfile.mkdtemp(prefix="biofasting_qc_")) / "qc_passed.fastq"

    import biofasting

    total = kept = failed_n = failed_q = failed_gc = 0
    kept_records = []
    bases = 0
    for title, sequence, quality in biofasting.open_fastq(args.path):
        total += 1
        if sequence.count(b"N"):
            failed_n += 1
            continue
        #   Per-record measurement, with the package's own fraction: `N` is
        #   already filtered above, so "remove" and any other ambiguous mode
        #   give the same number here, and the tutorial's claim about modes
        #   is where that is discussed, not here.
        gc = biofasting.gc_fraction(sequence) * 100.0
        if not args.gc_low <= gc <= args.gc_high:
            failed_gc += 1
            continue
        mean_q = sum(quality) / len(quality) - 33.0
        if mean_q < args.min_q:
            failed_q += 1
            continue
        kept += 1
        bases += len(sequence)
        kept_records.append((title, sequence, quality))
    biofasting.write_fastq(kept_records, args.out)

    print(f"read {total:,} records from {args.path}")
    print(
        f"kept {kept:,} records ({kept / total * 100:.1f}%); "
        f"dropped: N {failed_n:,}, quality {failed_q:,}, GC {failed_gc:,}"
    )
    print(f"kept reads: {bases:,} bp, mean length {bases / kept:.1f}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
