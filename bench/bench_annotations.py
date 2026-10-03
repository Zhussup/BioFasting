#!/usr/bin/env python3
"""The annotations dict, measured against the target M22 named for it.

`bench/targets.md` records the target before the kernel existed: a GenBank
record whose header carries an ordinary three-line comment costs the reference
**31.9 us** more than the same record without one, EMBL's costs **17.9 us**, and
the bulk of the GenBank figure is not reading the comment but the quadratic
structured-comment probe in front of the reader.  This script is the other half
-- the delivered reader, on the same generated corpus, under the same digest
gate.

Each format is measured as a **pair**, because the header is the only difference
between the two rows and the difference of two per-record times is therefore the
cost of the header and nothing else:

============  ==========================  ==========================
format        bare                        annotated
============  ==========================  ==========================
GenBank       ``genbank-1kb-bare``        ``genbank-1kb-annot``
EMBL          ``embl-1kb-bare``           ``embl-1kb-annot``
============  ==========================  ==========================

The four columns are the same measurement at four depths, and the first two are
one measurement split in two:

* **reference** -- `SeqIO.parse` over the row.  Everything: sequences, the
  annotations dict, the `Reference` objects.
* **ours** -- `open_genbank` + `read_annotations`: the kernel that walks the
  header, with no `Bio.SeqFeature.Reference` built.
* **ours+interop** -- the same plus `to_reference` for every reference, so the
  column that can be put beside the reference's object cost is separate from the
  column that is this milestone's work.

**The gate is over the dict, keys and order included.**  Every record's
annotations are reduced to its key list *in insertion order* and to plain
values -- references through `to_reference`, so the two sides are compared as
the reference's own objects -- and the whole file is hashed.  Two dicts are
equal under `==` whatever order they were built in, so a digest that ignored the
key order would pass a reader that emitted a fixed key list, which is precisely
the mistake this kernel is shaped to avoid.  A row that disagrees is printed
`INVALID` and no time from it may be quoted.

SwissProt is not a row here.  Its header is `Bio.SwissProt`, a different reader
with different keys, and M22 does not reproduce it, so there is nothing to
compare -- the same scope statement the features bench makes for its `FT` block.

    python3 bench/bench_annotations.py
    python3 bench/bench_annotations.py --only genbank-1kb-annot --repeat 9
    python3 bench/bench_annotations.py --list
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import seqio_corpus  # noqa: E402
from bench_genbank import CORPUS_DIR, EXTENSIONS, timed  # noqa: E402
from Bio import SeqIO  # noqa: E402

#: The pairs, each `(bare, annotated)`.  Named rather than derived, because
#: which rows are twins is a fact about how the corpus was generated -- the
#: annotated row draws from the same seed and the same domain string as its bare
#: twin, and applies the blocks afterwards -- and `tests/test_seqio_corpus.py`
#: asserts it rather than trusting it.
PAIRS = (
    ("genbank-1kb-bare", "genbank-1kb-annot"),
    ("embl-1kb-bare", "embl-1kb-annot"),
)

ROWS = tuple(group for pair in PAIRS for group in pair)

#: Records per row, from the corpus's own spec rather than from a parse: the
#: denominator of a per-record time has to come from the generator, or the row
#: would be measuring a number it also produced.
RECORDS = {row[0]: row[2] for row in seqio_corpus.CORPORA + seqio_corpus.ANNOT_CORPORA}


def corpus_file(group: str) -> tuple[str, Path]:
    """The `(format, path)` of one row, written if it is not already there.

    `bench_genbank.corpus_file` looks only in `CORPORA`; these rows live in
    `ANNOT_CORPORA`, so the lookup goes through `seqio_corpus.spec`, which is the
    one place that knows about both.
    """
    format, text = seqio_corpus.build_row(group)
    CORPUS_DIR.mkdir(parents=True, exist_ok=True)
    path = CORPUS_DIR / f"{group}.{EXTENSIONS[format]}"
    if not path.exists() or path.read_text() != text:
        path.write_text(text, newline="")
    return format, path


#: The work each timed row does, kept apart from the reduction that follows it.
#:
#: This split is the whole of the measurement's honesty.  The reduction below
#: costs more on the annotated row than on the bare one -- three `Reference`
#: objects and eleven keys against none and nine -- so a timing that included it
#: would put that difference inside the delta and report a header that costs more
#: than it does.  `timed` therefore walks only the parse (or only the kernel),
#: and the gate runs afterwards, on the result.
def reference_parse(format: str, path: Path) -> list:
    """`SeqIO.parse` over the whole file, which is the reference's full work."""
    with path.open() as handle:
        return list(SeqIO.parse(handle, format))


def reference_annotations(records) -> list:
    """The parsed records, reduced to the rows the gate compares."""
    return [
        (
            record.id,
            tuple(record.annotations),
            tuple(
                (key, _plain(record.annotations[key])) for key in record.annotations
            ),
        )
        for record in records
    ]


def our_parse(format: str, path: Path, interop: bool) -> list:
    """Our reader over the whole file -- and optionally through the conversion.

    `interop` decides whether a `Reference` object is built and nothing else:
    both settings produce the same rows, which the gate checks, so the two timed
    rows are the same work plus the object construction.
    """
    import biofasting

    index = biofasting.open_genbank(path, format)
    annotations = []
    for record_id in index:
        built = biofasting.read_annotations(index, record_id)
        if interop:
            built = {
                key: [biofasting.to_reference(item) for item in value]
                if key == "references"
                else value
                for key, value in built.items()
            }
        annotations.append((record_id, built))
    return annotations


def our_annotations(parsed) -> list:
    """Our rows, reduced the same way."""
    return [
        (
            record_id,
            tuple(annotations),
            tuple((key, _plain(annotations[key])) for key in annotations),
        )
        for record_id, annotations in parsed
    ]


#: The eight fields a `Reference` holds, in the order its own `__dict__` was
#: built in.  Both sides are read through this list, so "the same reference" has
#: one definition here and it is the reference's.
REFERENCE_FIELDS = (
    "title",
    "authors",
    "consrtm",
    "journal",
    "pubmed_id",
    "medline_id",
    "comment",
    "location",
)


def _plain(value):
    """A value reduced to something `repr` spells the same way on both sides.

    A reference -- the reference's own `Reference`, or ours when `interop` is
    off -- becomes a tuple of its eight fields; a location becomes the triples
    `(start, end, strand)` its parts hold, with our own two-element pairs given
    the reference's `None` strand so the two spellings agree.  A difference in
    the answer must not be able to hide behind a difference in how the answer is
    written.
    """
    if hasattr(value, "pubmed_id") and hasattr(value, "location"):
        fields = []
        for field in REFERENCE_FIELDS:
            item = getattr(value, field)
            if field == "location":
                item = tuple(
                    (int(part.start), int(part.end), part.strand)
                    if hasattr(part, "start")
                    else (int(part[0]), int(part[1]), None)
                    for part in item
                )
            fields.append((field, _plain(item)))
        return tuple(fields)
    if isinstance(value, (list, tuple)):
        return tuple(_plain(item) for item in value)
    return value


def digest_annotations(rows) -> str:
    """One hash over every key, its position and every value, so a row is all or
    nothing."""
    hasher = hashlib.sha256()
    for record_id, keys, values in rows:
        for field in (record_id, keys, values):
            data = repr(field).encode()
            hasher.update(len(data).to_bytes(8, "little"))
            hasher.update(data)
    return hasher.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", action="append", metavar="GROUP",
                        help="measure only this corpus row; repeatable")
    parser.add_argument("--repeat", type=int, default=5,
                        help="timed passes per implementation (default 5)")
    parser.add_argument("--list", action="store_true",
                        help="list the corpus rows and exit")
    args = parser.parse_args(argv)

    if args.list:
        for group in ROWS:
            print(group)
        return 0

    groups = list(ROWS)
    if args.only:
        unknown = [name for name in args.only if name not in groups]
        if unknown:
            parser.error("unknown corpus row: " + ", ".join(unknown))
        groups = [group for group in groups if group in args.only]

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

    header = (f"{'corpus':20s} {'rec':>6s} {'keys':>5s} {'ref us/rec':>11s} "
              f"{'ours':>8s} {'ref/ours':>9s} {'+interop':>9s} {'gate':>7s}")
    print(header)
    print("-" * len(header))

    per_record: dict[str, dict[str, float]] = {}
    invalid = 0
    for group in groups:
        format, path = corpus_file(group)
        expected = reference_annotations(reference_parse(format, path))
        ours = our_annotations(our_parse(format, path, interop=False))
        with_interop = our_annotations(our_parse(format, path, interop=True))

        gate = (
            "OK"
            if digest_annotations(ours) == digest_annotations(expected)
            and digest_annotations(with_interop) == digest_annotations(expected)
            else "INVALID"
        )
        invalid += gate != "OK"

        count = RECORDS[group]
        keys = len(expected[0][1]) if expected else 0
        ref = timed(lambda: reference_parse(format, path), args.repeat)
        got = timed(lambda: our_parse(format, path, interop=False), args.repeat)
        plus = timed(lambda: our_parse(format, path, interop=True), args.repeat)

        per_record[group] = {
            "reference": ref / count,
            "ours": got / count,
            "interop": plus / count,
        }

        print(f"{group:20s} {count:6d} {keys:5d} "
              f"{per_record[group]['reference'] * 1e6:11.2f} "
              f"{per_record[group]['ours'] * 1e6:8.2f} "
              f"{ref / got:8.2f}x {per_record[group]['interop'] * 1e6:9.2f} "
              f"{gate:>7s}")

    # The marginal cost of the annotations dict: two per-record times in a pair,
    # subtracted.  The two rows differ in the header and in nothing else -- the
    # same residues from the same seed, the blocks applied after the draw -- so
    # the difference is the header, which is exactly what M22 measured at 31.95
    # us for GenBank and 17.89 for EMBL.
    print()
    targets = {"genbank": 31.95, "embl": 17.89}
    for bare, annotated in PAIRS:
        if bare not in per_record or annotated not in per_record:
            continue
        format = annotated.split("-")[0]
        print(f"per record, marginal over {bare} "
              f"({RECORDS[annotated]} records with a full header against "
              f"{RECORDS[bare]} with a thin one):")
        for label in ("reference", "ours", "interop"):
            delta = per_record[annotated][label] - per_record[bare][label]
            print(f"  {label:22s} {delta * 1e6:7.2f} us")
        print(f"  M22's recorded target    {targets[format]:7.2f} us")
        print("  the reference's own absolute times do not reproduce on a re-run "
              "of the same size, but the differences do; the recorded spreads "
              "are in bench/targets.md")
        print()

    return 1 if invalid else 0


if __name__ == "__main__":
    raise SystemExit(main())
