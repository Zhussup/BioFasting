#!/usr/bin/env python3
"""The three sequence-format writers, measured against the target M23 set.

`bench/rank_writers.py` (the seventh ranking pass) asked what the writers cost
and found three different answers: a reference already at the Python floor
(FASTA, whose eighteen `write` calls a record are the only thing left), a
reference 4.81× above a floor reachable with `bytes.translate` (FASTQ), and the
worst writer in the family in absolute terms (QUAL, 15.4 µs for 150 bases
because every score goes through `"%i" % round(q, 0)` and the lines are packed
by popping from the front of a list).  That pass recorded the targets:

| row | target |
|---|---:|
| `fasta-1kb` | ≤ 1.0 µs/record |
| `fastq-150bp` | ≤ 1.0 µs/record |
| `qual-150bp` | ≤ 2.0 µs/record |

This script is the other half of that: the delivered kernels, on the same
corpus, held to those numbers.  The gate runs first and is over the whole file
-- the bytes our writer emits against `SeqIO.write` over the same records -- and
a row that fails it is printed `INVALID` and no time from it may be quoted.  A
writer that is fast and wrong is the one failure a writer benchmark cannot see
from a ratio.

Two ways of reaching the same rows are timed, and the difference between them is
the point:

* **ours** -- `write_fasta`/`write_fastq`/`write_qual`, which build the whole
  file in one C++ buffer and cross into the handle once.
* **ours + interop** -- the same, from `SeqRecord` objects through
  `from_seqrecord`, which is the path a caller holding Biopython objects takes.
  It is timed separately because a writer is not fast if the conversion in front
  of it is not, and the quality encoding is where the conversion spends its time.

The quality encoder gets a row of its own: `phred_to_sanger` against the
reference's `_get_sanger_quality_str`, which is the 4.81× the ranking pass
identified and the one place a `SeqRecord` becomes the string these writers take.

    python3 bench/bench_writers.py
    python3 bench/bench_writers.py --only fastq-150bp --repeat 9
    python3 bench/bench_writers.py --list
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from bench_genbank import corpus_file, timed  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.SeqIO.QualityIO import _get_sanger_quality_str  # noqa: E402

#: FASTQ corpus, the file the Phase 1 FASTQ benchmark measures.
FASTQ_CORPUS = HERE / "data" / "fastq" / "reads_10k.fastq"

#: The row whose 1 kb records the FASTA writer is measured on.  FASTA has no
#: generated corpus of its own here -- `bench/data/fasta/genome.fasta` is one
#: record of 100 Mbp, which is a different question -- so the FASTA row is the
#: flat-file corpus's records written back out in the other format.
FASTA_SOURCE = "genbank-1kb"

#: The target from `bench/targets.md`, µs a record, and the format each row is
#: written in.  The payload floor is recorded there too; nothing here may claim
#: to be below it, and the table re-measures it in the same harness so the claim
#: stays checkable.
TARGETS = {
    "fasta-1kb": (1.0, "fasta"),
    "fastq-150bp": (1.0, "fastq"),
    "qual-150bp": (2.0, "qual"),
}

ROWS = tuple(TARGETS)


# --- the records -------------------------------------------------------------


def fastq_pairs() -> list:
    """Our own ``(title, sequence, quality)`` triples, from our own reader.

    Read once, outside the timer: this script measures the writers, and a
    reader in front of them would put its cost in every column.
    """
    import biofasting

    return list(biofasting.open_fastq(FASTQ_CORPUS))


def fasta_pairs() -> list:
    """``(title, sequence)`` pairs, from the flat-file corpus's records.

    Through ``corpus_file`` so the row uses the same cached file the flat-file
    scripts do, and parsed outside the timer: this script measures the writers,
    and a reader in front of them would put its cost in every column.
    """
    fmt, path = corpus_file(FASTA_SOURCE)
    with path.open() as handle:
        records = list(SeqIO.parse(handle, fmt))
    return [
        (record.description.encode("latin-1"), str(record.seq).encode("latin-1"))
        for record in records
    ]


def seqrecords(pairs, qualities: bool) -> list:
    """Our rows as ``SeqRecord`` objects, for the reference to write."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    out = []
    for row in pairs:
        title = row[0].decode("latin-1")
        record = SeqRecord(
            Seq(row[1].decode("latin-1")),
            id=title.split()[0] if title.split() else "",
            description=title,
        )
        if qualities:
            record.letter_annotations["phred_quality"] = [
                byte - 33 for byte in row[2]
            ]
        out.append(record)
    return out


def rows_for(name: str):
    """``(ours, theirs)`` for one row: our pairs and the same records as objects."""
    fmt = TARGETS[name][1]
    if fmt == "fasta":
        pairs = fasta_pairs()
        return pairs, seqrecords(pairs, qualities=False)
    pairs = fastq_pairs()
    return pairs, seqrecords(pairs, qualities=True)


# --- what is timed -----------------------------------------------------------


def reference_write(fmt: str, records, handle) -> None:
    """``SeqIO.write`` -- what a caller would otherwise use."""
    SeqIO.write(records, handle, fmt)


def our_write(fmt: str, pairs, handle) -> None:
    import biofasting

    if fmt == "fasta":
        biofasting.write_fasta(pairs, handle)
    elif fmt == "fastq":
        biofasting.write_fastq(pairs, handle)
    else:
        biofasting.write_qual([(title, quality) for title, _, quality in pairs],
                              handle)


def our_write_from_records(fmt: str, records, handle) -> None:
    """The same output, reached from ``SeqRecord`` objects.

    `from_seqrecord` is the conversion a caller holding Biopython objects pays
    before the writer sees anything, so this row is the writer plus that
    conversion -- and for FASTQ the conversion is the quality encoding, which is
    where the reference's per-base dict lookup lives.
    """
    import biofasting

    if fmt == "fasta":
        biofasting.write_fasta(
            (biofasting.from_seqrecord(record) for record in records), handle
        )
    elif fmt == "fastq":
        biofasting.write_fastq(
            (biofasting.from_seqrecord(record) for record in records), handle
        )
    else:
        biofasting.write_qual(
            (
                (title, quality)
                for title, _, quality in map(biofasting.from_seqrecord, records)
            ),
            handle,
        )


def payload(pairs, fmt: str) -> bytes:
    """The raw bytes with no formatting above them: the floor under every row.

    Built inside the timed region on purpose.  `bytes(existing_bytes)` returns
    the same object, and timing that is how the ranking pass first reported a
    payload of 0.00 µs a record.
    """
    chunks = []
    for row in pairs:
        chunks.append(bytes(row[1]))
        if fmt in ("fastq", "qual"):
            chunks.append(bytes(row[2]))
    return b"".join(chunks)


def reference_quality(records) -> bytes:
    """The reference's per-base quality encoding, for the encoder row."""
    return b"".join(
        _get_sanger_quality_str(record).encode("ascii") for record in records
    )


def our_quality(records) -> bytes:
    """The kernel, on the same scores the reference's encoder is handed.

    Both sides start from ``record.letter_annotations["phred_quality"]``, so the
    row is the encoding and nothing else -- the ``+interop`` column above is
    what the whole conversion costs with the title and the sequence in it.
    """
    from biofasting import _core

    return b"".join(
        _core.phred_to_sanger(bytes(record.letter_annotations["phred_quality"]))
        for record in records
    )


# --- the gate ----------------------------------------------------------------


def gate(fmt: str, pairs, records) -> tuple[bool, str]:
    """Byte identity of the whole file against the reference, before any timing.

    Over the file and not over a record, because a writer either emits the
    format or it does not, and the interesting failure -- one byte wrong in a
    wrap -- is invisible in a spot check.
    """
    expected = io.StringIO()
    reference_write(fmt, records, expected)
    got = io.BytesIO()
    our_write(fmt, pairs, got)
    if got.getvalue() != expected.getvalue().encode("latin-1"):
        return False, "the writer's bytes differ from the reference's"
    converted = io.BytesIO()
    our_write_from_records(fmt, records, converted)
    if converted.getvalue() != got.getvalue():
        return False, "the SeqRecord path differs from the pair path"
    return True, ""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", action="append", metavar="ROW",
                        help="measure only this row; repeatable")
    parser.add_argument("--repeat", type=int, default=7,
                        help="timed passes per implementation (default 7)")
    parser.add_argument("--list", action="store_true",
                        help="list the rows and exit")
    args = parser.parse_args(argv)

    if args.list:
        for name in ROWS:
            print(name)
        return 0

    names = list(ROWS)
    if args.only:
        unknown = [name for name in args.only if name not in names]
        if unknown:
            parser.error("unknown row: " + ", ".join(unknown))
        names = [name for name in names if name in args.only]

    try:
        import biofasting
    except ImportError:
        print("biofasting is not importable; build it with "
              "`pip install --no-build-isolation -e .`", file=sys.stderr)
        return 2

    info = biofasting.build_info()
    print(f"{info['version']}  {info['compiler']} c++{info['cxx_standard']}  "
          f"flags={info['extra_cxx_flags']}  {info['arch']}  "
          f"seqops={biofasting.seqops_level()}  cpu={biofasting.best_cpu_level()}")
    print(f"Python {sys.version.split()[0]}, median of {args.repeat} passes "
          f"after one warm-up\n")

    header = (f"{'row':14s} {'recs':>6s} {'ref us/rec':>11s} {'ours':>8s} "
              f"{'ratio':>7s} {'+interop':>9s} {'payload':>8s} {'target':>8s} "
              f"{'gate':>7s}")
    print(header)
    print("-" * len(header))

    invalid = 0
    missed = 0
    for name in names:
        target, fmt = TARGETS[name]
        pairs, records = rows_for(name)
        ok, why = gate(fmt, pairs, records)
        if not ok:
            invalid += 1
            print(f"{name:14s} {len(pairs):6d} INVALID -- {why}")
            continue
        count = len(pairs)

        ref = timed(lambda: reference_write(fmt, records, io.StringIO()),
                    args.repeat)
        ours = timed(lambda: our_write(fmt, pairs, io.BytesIO()), args.repeat)
        interop = timed(lambda: our_write_from_records(fmt, records, io.BytesIO()),
                        args.repeat)
        floor = timed(lambda: payload(pairs, fmt), args.repeat)

        per_record = ours / count * 1e6
        verdict = "PASS" if per_record <= target else "MISS"
        missed += verdict == "MISS"
        print(f"{name:14s} {count:6d} {ref / count * 1e6:11.2f} "
              f"{per_record:8.3f} {ref / ours:6.2f}x {interop / count * 1e6:9.3f} "
              f"{floor / count * 1e6:8.3f} {target:7.1f} {verdict:>7s}")

    # The quality encoder on its own: the 4.81x the ranking pass found, and the
    # only place a `SeqRecord` becomes the string the writers take.
    if "fastq-150bp" in names:
        pairs, records = rows_for("fastq-150bp")
        count = len(pairs)
        ours_bytes = our_quality(records)
        theirs_bytes = reference_quality(records)
        same = ours_bytes == theirs_bytes
        encoder = timed(lambda: reference_quality(records), args.repeat)
        kernel = timed(lambda: our_quality(records), args.repeat)
        print()
        print("the quality encoder, from SeqRecord to the Sanger string:")
        print(f"  reference (_get_sanger_quality_str)  {encoder / count * 1e6:7.3f} "
              "us/rec")
        print(f"  ours (phred_to_sanger)               {kernel / count * 1e6:7.3f} "
              f"us/rec   {encoder / kernel:.2f}x")
        print(f"  gate: {'OK' if same else 'INVALID'} "
              f"({len(ours_bytes)} bytes)")
        print("  the kernel is the encoding alone; `+interop` above is the whole")
        print("  `from_seqrecord` conversion with the title and the sequence in it.")
        invalid += not same

    print()
    print("`ours` builds the whole file in one buffer and crosses into the handle")
    print("once; `+interop` is the same from SeqRecord objects, which is the path a")
    print("caller holding Biopython objects takes.  `payload` is the raw bytes with")
    print("no formatting, measured in this harness so the floor stays checkable --")
    print("the ranking pass's `rewrite` column is the pure-Python floor, which is a")
    print("bound on a rewrite and not on a kernel.")
    if invalid:
        print(f"\n{invalid} row(s) INVALID: the gate disagrees with the reference "
              "and no time from them may be quoted.", file=sys.stderr)
        return 1
    if missed:
        print(f"\n{missed} row(s) MISS the target recorded in bench/targets.md.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
