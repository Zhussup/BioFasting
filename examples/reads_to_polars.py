#!/usr/bin/env python3
"""Reads into a polars frame: the Arrow table, handed over and aggregated.

The pipeline the Arrow tables were built for: `read_fastq_table` lays a whole
FASTQ file into Arrow's own buffers in one pass from C++, `pl.from_arrow`
adopts that table into a frame, and every question afterwards is one polars
expression -- here a length column, a GC count, and the summary rows.

The hand-off has a documented price (bench/targets.md, the tenth Arrow
measurement): `from_arrow` copies no strings -- the frame keeps this table's
bytes alive -- but it does build polars' own row index, which is O(rows) and
not free.  A caller who stays in pyarrow never pays it; this example pays it
once, on purpose, because the aggregation is the point of leaving.

Everything after that is frame work and is not this package's business: the
sequence column is a `large_string` the reader already assembled.  Length is
`str.len_bytes`; GC counts use `count_matches`, whose default is regex
matching, which single letters survive.

    python3 examples/reads_to_polars.py --path bench/data/fastq/reads_10k.fastq

Needs polars and pyarrow; without either it exits 1 with a sentence instead
of a traceback.  The tests skip this example by `find_spec`, for the same
reason.

Usage:

    python3 examples/reads_to_polars.py [--path FILE] [--rows N]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).parent
DEFAULT = HERE.parent / "bench" / "data" / "fastq" / "reads_10k.fastq"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--path", type=Path, default=DEFAULT, help="the FASTQ to read")
    parser.add_argument("--rows", type=int, default=5, help="how many frame rows to show")
    args = parser.parse_args()
    if not args.path.exists():
        print(f"{args.path} is missing: python3 bench/gen_data.py", file=sys.stderr)
        return 1

    try:
        import polars
    except ImportError:
        print(
            "polars is not installed: pip install polars "
            "(pyarrow comes with biofasting[arrow])",
            file=sys.stderr,
        )
        return 1

    try:
        import biofasting
        table = biofasting.read_fastq_table(args.path)
    except ImportError as err:
        #   `read_fastq_table` raises the pyarrow wheel's own ImportError,
        #   with the extra's name spelled out in it; the message is the
        #   sentence, not a traceback.
        print(str(err), file=sys.stderr)
        return 1

    frame = polars.from_arrow(table)
    lengths = polars.col("sequence").str.len_bytes()
    gc = (
        polars.col("sequence").str.count_matches("G")
        + polars.col("sequence").str.count_matches("C")
    )
    counted = frame.select(length=lengths, gc=gc)
    summary = counted.select(
        polars.len().alias("records"),
        polars.col("length").sum().alias("bp total"),
        polars.col("length").mean().round(1).alias("mean length"),
        (polars.col("gc").sum() / polars.col("length").sum()).alias("pooled GC"),
        polars.col("gc").min().alias("min GC"),
        polars.col("gc").max().alias("max GC"),
    )
    print("the frame's first rows:")
    print(counted.head(args.rows))
    print("\nthe whole file, as one polars summary:")
    print(summary)
    print(
        "\nthe table still holds the bytes the frame reads from "
        f"({table.num_rows} rows, {len(table.column_names)} columns, "
        f"schema {', '.join(table.column_names)})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())