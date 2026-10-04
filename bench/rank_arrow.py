#!/usr/bin/env python3
"""Arrow interop, measured before it is built (PLAN 2.3).

PLAN 2.3 asks for two things at once: an Arrow interop for this package, and a
*compatibility story* with `polars-bio` rather than a collision with it.  Both
halves are measurable before either is designed, so this pass measures them.

The question is not "is Arrow fast" -- Arrow is a memory layout, not a
program -- but **what a caller pays today to get sequence records into a table,
and what the best tool in that niche already pays.**  So the rows are the
corpus, and the columns are the columns `polars-bio` itself produces, so that
every producer's table can be compared before any time is quoted.

Four producers per row:

* **reference** -- `SeqIO.parse` into Python lists into a table.  What a
  Biopython user writes today; it is the thing this package exists to be
  faster than.
* **polars-bio** -- `polars_bio.read_fastq` / `read_fasta`, a Rust reader
  (needletail) that produces the same four columns natively.  This is the
  best tool in the niche and therefore the real target, not the reference.
* **rewrite** -- the smallest pure-Python producer that gets the same bytes,
  with per-record loops replaced by C-level calls: one `bytes.split(b"\\n")`
  for the whole file and four slices, one `bytearray` walk for FASTA.  It is a
  **floor** under any C++ kernel -- measured, not assumed -- and it is what
  decides whether a kernel is worth writing at all.
* **ours** -- this package's reader (`open_fastq_index`, `open_fasta`) plus a
  table build in Python.  It is deliberately not a kernel: it isolates how much
  of the cost is the *reading* we already have and how much is the *building*
  that does not exist yet.  The difference between `ours` and `rewrite` is the
  headroom a C++ builder would have to take.

**The gate.**  Every producer's table is compared against the reference's,
column by column and value by value, before the row's times are printed.  A row
that disagrees prints `DIFF` with the first offending column and record and
returns 1: a table that is not the same table is not a faster table.

Usage:

    python3 bench/rank_arrow.py
"""

from __future__ import annotations

import os
import statistics
import sys
import time

import pyarrow as pa
import polars as pl

REPEATS = 5

# `repeats` overrides the default per row: the two 100 MB rows cost seconds a
# pass and would dominate the run, and their spread is small enough that three
# passes is a median.
ROWS: dict[str, tuple[str, str, int]] = {
    "fastq-10k": ("bench/data/fastq/reads_10k.fastq", "fastq", 5),
    "fastq-gz-1m": ("bench/data/fastq/reads_1m.fastq.gz", "fastq", 3),
    "fastq-1m": ("bench/data/fastq/reads_1m.fastq", "fastq", 3),
    "fasta-1mb": ("bench/data/fasta/genome_1mb.fasta", "fasta", 5),
    "fasta-100mb": ("bench/data/fasta/genome.fasta", "fasta", 3),
}

STRING = pa.large_string()


def timed(fn, repeats: int) -> float:
    """Median seconds over ``repeats`` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def table(columns: dict[str, list]) -> pa.Table:
    """One table from plain Python lists, which is what all four producers do."""
    return pa.table({name: pa.array(values, type=STRING) for name, values in columns.items()})


# --------------------------------------------------------------------------
# The producers
# --------------------------------------------------------------------------


def reference(path: str, kind: str) -> pa.Table:
    """`SeqIO.parse` into lists into a table -- the Biopython path.

    The quality column is built with the same per-base dictionary lookup
    `SeqIO.QualityIO.FastqPhredWriter` uses internally, because that is what
    Biopython exposes: `letter_annotations` holds *integers*, and a caller who
    wants the encoded string pays for the encoding.  The trap is that a
    `polars-bio` column is the encoded string, so a gate against it cannot
    compare an int column and call the difference a bug.
    """
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import _phred_to_sanger_quality_str

    # Biopython 1.88 does not decompress by magic bytes -- `SeqIO.parse` on a
    # `.gz` dies with `UnicodeDecodeError` on 0x8b -- so the caller has to open
    # the stream, which is the reference's own documented answer and part of
    # its cost.
    if open(path, "rb").read(2) == b"\x1f\x8b":
        import gzip

        handle = gzip.open(path, "rt", encoding="latin-1")
    else:
        handle = open(path, encoding="latin-1")
    columns: dict[str, list] = {"name": [], "description": [], "sequence": []}
    if kind == "fastq":
        columns["quality"] = []
    for record in SeqIO.parse(handle, kind):
        name = record.id
        description = record.description
        if description.startswith(name):
            description = description[len(name) :].lstrip()
        columns["name"].append(name)
        columns["description"].append(description)
        columns["sequence"].append(str(record.seq))
        if kind == "fastq":
            columns["quality"].append(
                "".join(
                    _phred_to_sanger_quality_str[q]
                    for q in record.letter_annotations["phred_quality"]
                )
            )
    return table(columns)


def polars_bio(path: str, kind: str) -> pa.Table:
    """`polars-bio`'s Rust reader -- the tool this has to be compatible with.

    One divergence is normalised here rather than passed on, and it is a
    *finding*, not a bug being hidden: a FASTA header with no description
    (`>chrS`) gives Biopython and this package the empty string and gives
    `polars-bio` a **null**.  The two are not the same in a table -- `null` is
    the absence of a value and `""` is a value of length zero -- so a caller
    switching between the tools sees a different frame, and the compatibility
    story has to say so.  The column is filled here so the rest of the gate can
    stay strict.
    """
    import polars_bio as pb

    frame = pb.read_fastq(path) if kind == "fastq" else pb.read_fasta(path)
    frame = frame.with_columns(
        pl.col(column).fill_null("")
        for column in ("name", "description")
        if frame.schema[column] == pl.String
    )
    if "quality_scores" in frame.columns:
        frame = frame.rename({"quality_scores": "quality"})
    return frame.to_arrow()


def rewrite(path: str, kind: str) -> pa.Table:
    """The pure-Python floor: C-level calls, no per-character work.

    FASTQ is exactly four lines a record in this corpus, so `split` on the
    whole file and four slices do the whole parse; FASTA is a `>` split and a
    newline strip per record.  Both are *cheating* in the sense that neither
    validates anything beyond this corpus's shape -- that is what a floor is
    for -- and the gate is what proves the cheating did not change the table.
    """
    with open(path, "rb") as handle:
        data = handle.read()
    if data[:2] == b"\x1f\x8b":
        # `gzip.decompress` is zlib in a C loop, so the floor stays a floor:
        # what it does not do is the streaming, the sniffing or the index.
        import gzip

        data = gzip.decompress(data)
    data = data.decode("latin-1")

    columns: dict[str, list] = {"name": [], "description": [], "sequence": []}
    if kind == "fastq":
        columns["quality"] = []
        lines = data.split("\n")
        if lines and lines[-1] == "":
            lines.pop()
        headers = lines[0::4]
        columns["sequence"] = lines[1::4]
        # The '+' line is a separator and its content is not read by this
        # package's parser either; the reference keeps whatever follows the
        # '+', and the corpus repeats the header there.
        columns["quality"] = lines[3::4]
        for header in headers:
            name, _, description = header[1:].partition(" ")
            columns["name"].append(name)
            columns["description"].append(description)
        return table(columns)

    for block in data.split(">"):
        if not block:
            continue
        header, _, body = block.partition("\n")
        name, _, description = header.partition(" ")
        columns["name"].append(name)
        columns["description"].append(description)
        columns["sequence"].append("".join(body.split("\n")).rstrip("\n"))
    return table(columns)


def ours(path: str, kind: str) -> pa.Table:
    """This package's reader, with the table built in Python.

    The reader is the delivered one and it is the **sequential** one: for FASTQ
    the scanner (`open_fastq`) rather than the index, one C++ iterator handing
    back a header, a sequence and a quality per record, and for FASTA the
    memory-mapped index.  What this row measures is therefore the *build* -- one
    Python call per record into a list -- and that is exactly the part a C++
    builder would have to remove.  An earlier version of this row fetched
    through `open_fastq_index`, which is one hash lookup and four interpreter
    crossings a record and was **2.65 µs** against the scanner's half of that;
    the index is for random access and a table is a scan, which is worth saying
    because the wrong reader is easy to pick.
    """
    import biofasting as bf

    columns: dict[str, list] = {"name": [], "description": [], "sequence": []}
    if kind == "fastq":
        columns["quality"] = []
        for header, sequence, quality in bf.open_fastq(path):
            name, _, description = header.partition(b" ")
            columns["name"].append(name.decode("latin-1"))
            columns["description"].append(description.decode("latin-1"))
            columns["sequence"].append(sequence.decode("latin-1"))
            columns["quality"].append(quality.decode("latin-1"))
        return table(columns)

    index = bf.open_fasta(path)
    for name in index.keys():
        title = index.title(name)
        columns["name"].append(name)
        columns["description"].append(title[len(name) :].lstrip())
        columns["sequence"].append(index[name].decode("latin-1"))
    return table(columns)


PRODUCERS = {
    "reference": reference,
    "polars-bio": polars_bio,
    "rewrite": rewrite,
    "ours": ours,
}


def gate(row: str, expected: pa.Table, produced: pa.Table) -> str | None:
    """`None` if the two tables are the same table, else the first difference."""
    if produced.num_rows != expected.num_rows:
        return f"{produced.num_rows} rows against {expected.num_rows}"
    if produced.column_names != expected.column_names:
        return f"columns {produced.column_names} against {expected.column_names}"
    for name in expected.column_names:
        left = expected.column(name).to_pylist()
        right = produced.column(name).to_pylist()
        if left != right:
            for i, (a, b) in enumerate(zip(left, right)):
                if a != b:
                    return f"column {name!r} record {i}: {a!r} != {b!r}"
            return f"column {name!r} differs at the end"
    return None


def main(argv: list[str]) -> int:
    ok = True
    tables: dict[str, dict[str, pa.Table]] = {}
    for row, (path, kind, _repeats) in ROWS.items():
        tables[row] = {name: fn(path, kind) for name, fn in PRODUCERS.items()}

    print("gate: every producer's table against the reference's, column by column")
    for row, produced in tables.items():
        expected = produced["reference"]
        for name, got in produced.items():
            if name == "reference":
                continue
            difference = gate(row, expected, got)
            if difference is None:
                print(f"  {row}: {name} equal, {got.num_rows} rows")
            else:
                print(f"  {row}: {name} DIFF {difference}")
                ok = False
    if not ok:
        return 1

    print()
    print("| row | records | MB | reference | polars-bio | rewrite | ours | ours/rewrite |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|")
    ours_over_rewrite = {}
    for row, (path, kind, repeats) in ROWS.items():
        records = tables[row]["reference"].num_rows
        megabytes = os.path.getsize(path) / 1e6
        times = {}
        for name, fn in PRODUCERS.items():
            times[name] = timed(lambda fn=fn, path=path, kind=kind: fn(path, kind), repeats)

        def rate(seconds: float) -> str:
            # Two units, because one of them lies on every row: FASTQ here is
            # 150-base records and reporting µs a record is honest; FASTA here
            # is five records of 2–40 Mbp and µs a record is a number about
            # memcpy, not about parsing.  So the FASTQ rows read µs a record
            # and the FASTA rows read MB/s.
            if kind == "fastq":
                return f"{seconds / records * 1e6:.3f} us/rec"
            return f"{megabytes / seconds:.0f} MB/s"

        ours_over_rewrite[row] = times["ours"] / times["rewrite"]
        print(
            f"| {row} | {records:,} | {megabytes:.1f} | "
            + " | ".join(rate(times[name]) for name in PRODUCERS)
            + f" | {ours_over_rewrite[row]:.2f}x |"
        )

    print()
    print("the headroom a C++ builder would be taking, per row:")
    for row, (path, kind, repeats) in ROWS.items():
        records = tables[row]["reference"].num_rows
        total = timed(lambda path=path, kind=kind: rewrite(path, kind), repeats)
        print(
            f"  {row}: the whole rewrite pass is {total * 1e3:.1f} ms for {records:,} records, "
            f"{total / records * 1e6:.3f} us a record"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
