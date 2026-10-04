#!/usr/bin/env python3
"""The delivered Arrow tables, measured against the target that asked for them.

`bench/rank_arrow.py` is the *target* half: it ranked four producers of the same
table and found that the tool to beat is not Biopython (15.8 µs a read) but
`polars-bio` (0.970 µs a read), and that this package's reader with the table
built in Python was already ahead on three of the five rows -- 1.85x on gzipped
FASTQ and 8.4x and 1.95x on the two FASTA rows -- while losing the plain FASTQ
row to a Rust reader that builds the columns as it scans.  The gap it named was
therefore not the reading, which is delivered, but the *building*, which was not.

This is the other half: the same five rows, the same four producers imported from
`rank_arrow.py` so that both tables are built by one piece of code, and one new
column -- `read_fastq_table` / `read_fasta_table`, the whole delivered call:
mapping the file, sniffing gzip by its magic bytes, scanning, splitting the
headers and writing Arrow's offsets and data buffers in C++.

**The gate is the table.**  Every producer's table is compared against the
reference's -- column by column, value by value, and against `polars-bio`'s as
well -- before any row's time is printed.  A row that disagrees prints `DIFF`
with the first offending column and record and contributes no time, because a
table that is not the same table is not a faster table.  That gate is not
ceremony: the two producers normalise a `null` description differently, and a
comparison that did not say which reading it was testing would have called the
difference a bug or hidden a real one.

The last two columns separate the two halves of the delivered path: `ours` is the
table, and `+ polars` is `pl.from_arrow(table)`, which is what a caller actually
does next.  The difference between them is the price of the hand-off, and what
this pass found is that it is **not** nothing: the strings are not copied --
measured by RSS, and the frame keeps them alive after the table object is gone --
but polars 1.44 builds its own per-string index over the adopted bytes, which is
O(rows) at 32-38 ns a row.  The last section times it on a table already in hand,
rather than by subtracting two timed calls: the rows here build hundreds of
megabytes, and two medians of that differ by more than the quantity.

Every row uses the same repeat count, including the two 100 MB ones.  `ROWS`
carries a per-row count for `rank_arrow.py`; this file deliberately ignores it,
for the reason recorded at `main`.

Usage:

    python3 bench/bench_arrow.py [repeats]
"""

from __future__ import annotations

import os
import statistics
import sys
import time

import pyarrow as pa

# Imported rather than copied: the reference, the competitor and the pure-Python
# floor are the target's own producers, and a second copy of any of them here
# would be a second thing to keep correct -- and the gate below is exactly a
# comparison against them.
from rank_arrow import ROWS, gate, polars_bio, reference, rewrite

REPEATS = 5


def timed(fn, repeats: int = REPEATS) -> float:
    """Median seconds over ``repeats`` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def delivered(path: str, kind: str) -> pa.Table:
    """The delivered call: the whole table, from a path, in one pass."""
    import biofasting as bf

    return (
        bf.read_fastq_table(path) if kind == "fastq" else bf.read_fasta_table(path)
    )


def delivered_to_polars(path: str, kind: str):
    """The same, plus the hand-off a caller makes next.

    Returns the *frame*, not a table: `to_arrow()` on a polars frame is a
    conversion with a cost of its own, and timing it here would put a number
    about polars' own plumbing in a column about this package's.  The gate
    converts, because a gate compares values and has no clock.
    """
    import polars as pl

    return pl.from_arrow(delivered(path, kind))


PRODUCERS = {
    "reference": reference,
    "polars-bio": polars_bio,
    "rewrite": rewrite,
    "ours": delivered,
    "+ polars": delivered_to_polars,
}


def main(argv: list[str]) -> int:
    # One repeat count for every row, including the two 100 MB ones.  `ROWS`
    # carries a per-row count for `rank_arrow.py`'s benefit and it is deliberately
    # not used here: the rows whose build allocates hundreds of megabytes are
    # exactly the ones a median of three cannot pin down -- the same call
    # measured 0.319 us/rec and 0.582 us/rec in two places in one run of this
    # file, and five repeats across three rounds gave 0.330 / 0.352 / 0.336 --
    # so a shorter run here would buy seconds and cost the number.
    runs = int(argv[1]) if len(argv) > 1 else REPEATS
    tables: dict[str, dict[str, pa.Table]] = {}
    for row, (path, kind, _repeats) in ROWS.items():
        produced = {name: fn(path, kind) for name, fn in PRODUCERS.items()}
        # The one producer that returns a frame rather than a table is converted
        # here and only here: the gate compares values, so it needs the table,
        # and the clock is not running yet.
        for name, got in produced.items():
            if hasattr(got, "to_arrow"):
                produced[name] = got.to_arrow()
        tables[row] = produced

    print("gate: every producer's table against the reference's, then ours against polars-bio's")
    failures = 0
    for row, produced in tables.items():
        expected = produced["reference"]
        for name, got in produced.items():
            if name == "reference":
                continue
            difference = gate(row, expected, got)
            if difference is None:
                print(f"  {row}: {name} equal to the reference, {got.num_rows} rows")
            else:
                print(f"  {row}: {name} DIFF {difference}")
                failures += 1
        # The competitor's own table, with its null descriptions normalised, is
        # the other reading a caller might have; both must be the same table.
        difference = gate(row, produced["polars-bio"], produced["ours"])
        if difference is None:
            print(f"  {row}: ours equal to polars-bio, {produced['ours'].num_rows} rows")
        else:
            print(f"  {row}: ours DIFF against polars-bio {difference}")
            failures += 1
    if failures:
        return 1

    # --- the table ---------------------------------------------------------
    print()
    print(
        "| row | records | MB | reference | polars-bio | rewrite | "
        "ours | + polars | vs polars-bio |"
    )
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    # The delivered call's own median, per row, kept for the hand-off section
    # below so that it is timed once and not twice.
    ours_seconds: dict[str, float] = {}
    for row, (path, kind, _repeats) in ROWS.items():
        records = tables[row]["reference"].num_rows
        megabytes = os.path.getsize(path) / 1e6
        times = {
            name: timed(lambda fn=fn, path=path, kind=kind: fn(path, kind), runs)
            for name, fn in PRODUCERS.items()
        }
        ours_seconds[row] = times["ours"]

        def rate(seconds: float) -> str:
            # Two units, because one of them lies on every row: these FASTQ
            # records are 150 bases and µs a record is honest, while the FASTA
            # rows are five records of 2-40 Mbp and µs a record would be a number
            # about memcpy rather than about parsing.  So FASTQ reads µs a
            # record and FASTA reads MB/s.
            if kind == "fastq":
                return f"{seconds / records * 1e6:.3f} us/rec"
            return f"{megabytes / seconds:.0f} MB/s"

        print(
            f"| {row} | {records:,} | {megabytes:.1f} | "
            + " | ".join(rate(times[name]) for name in PRODUCERS)
            + f" | {times['polars-bio'] / times['ours']:.2f}x |"
        )

    # --- what the hand-off costs -------------------------------------------
    #
    # Timed on a table that is already built, and not as `delivered_to_polars`
    # minus `delivered`.  That difference was the first attempt and it printed
    # negative numbers -- -0.0681 us a record, -18.16 ms on a 131.3 ms call --
    # because it is a difference of two medians of a call that allocates
    # hundreds of megabytes, and the noise between them is larger than the
    # quantity.  `pl.from_arrow` on a table in hand is the thing being claimed to
    # be free-ish, so it is the thing to time: whatever it reads, it reads from a
    # table the clock does not cover.
    #
    # It is not free, and this is where that was found.  `pl.from_arrow` does not
    # copy the *strings* -- measured by RSS, and the frame keeps them alive after
    # the table object is dropped -- but polars 1.44 does not adopt an Arrow
    # string column as-is either: it builds its own per-string index over the
    # adopted bytes, which is **O(rows)**, 32-38 ns a row on this host and 37.8 ms
    # on the million-row table.  Measured on three slices of one table (10k,
    # 100k and 1M rows) it is linear in rows and not in bytes.  So the column is
    # the price of the hand-off, and a caller who never leaves pyarrow does not
    # pay it.
    import polars as pl

    print()
    print("the hand-off (pl.from_arrow) on a table that is already built:")
    for row, (path, kind, _repeats) in ROWS.items():
        table = delivered(path, kind)
        handed = timed(lambda table=table: pl.from_arrow(table), runs)
        # The call it is a fraction of is the one the table above already
        # measured, read back out of `ours_seconds` rather than timed a second
        # time: a second median of the same call is the noise this section was
        # rewritten to get away from.
        measured = ours_seconds[row]
        records = table.num_rows
        # The per-row figure is printed only where a row is a record.  On the
        # FASTA rows the file is one or five records of megabytes, so "ns a row"
        # there is a number about memcpy with a division attached, and the honest
        # unit for those two is the total and the share.
        per_row = (
            f", {handed / records * 1e9:.1f} ns a row" if records >= 1000 else ""
        )
        print(
            f"  {row}: {handed * 1e6:.1f} us on a {measured * 1e3:.1f} ms call "
            f"({handed / measured * 100:.2f}% of it){per_row}, {records:,} rows"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
