#!/usr/bin/env python3
"""The FEATURES table, measured against the target M19 named for it.

`bench/rank_seqio.py` decomposed a flat-file parse and found that one feature
per kilobase costs the reference **20.2 us**, of which building the `SeqFeature`
with its location and three qualifiers is 4.2 and the other **16 us is the
reference reading the location string and the qualifier block**.  That is the
target this reader exists for -- M20's table sits beside its own the same way --
so this script is the other half: the delivered feature reader, on the same
generated corpus, under the same digest gate.

Four rows, and the first two are one measurement split in two:

* **reference** -- `SeqIO.parse` over the featured row.  Everything: sequences,
  the FEATURES table, the annotations dict.
* **reference bare** -- the same parse over the row that differs only in feature
  density.  Subtracting it leaves the reference's per-feature cost and nothing
  else, which is what makes the 20.2 us figure reproducible rather than quoted.
* **ours** -- `open_genbank` + `read_features`: the kernel that reads the table,
  with no `SeqFeature` in sight.
* **ours+interop** -- the same, plus `to_seqfeature` for every feature, so the
  column that can be put beside the reference's 4.2 us object cost is separate
  from the column that is this milestone's work.

**The gate is over the features, not over id and sequence.**  Every record's
every feature is reduced to its type, its location as the canonical nested
tuples `tests/test_location.py` defines, its qualifiers and its location status,
and the whole file is hashed: a row that disagrees is printed `INVALID` and no
time from it may be quoted.  A kernel that parsed faster and dropped a qualifier
would otherwise win this table.

SwissProt is not a row here.  Its `FT` block is a different grammar
(`Bio.SwissProt._read_ft`) and M21 does not reproduce it, so there is nothing to
compare.

    python3 bench/bench_features.py
    python3 bench/bench_features.py --only genbank-1kb --repeat 9
    python3 bench/bench_features.py --list
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
from bench_genbank import corpus_file, digest, timed  # noqa: E402
from Bio import SeqIO  # noqa: E402

# The row whose feature density is zero.  Per-feature cost is a *difference*
# between two rows that differ only in that, which is exactly how the ranking
# pass measured it -- and re-deriving it here is what makes the target
# checkable rather than a number this file trusts.
BARE = "genbank-1kb-bare"

# The rows that carry features.  SwissProt is excluded on scope, not on
# difficulty: its FT grammar is not INSDC's.
ROWS = tuple(row[0] for row in seqio_corpus.CORPORA if row[1] != "swiss" and row[4])

# Records per row, from the corpus's own spec rather than from a parse: the
# denominator of a per-record time has to come from the generator, or the row
# would be measuring a number it also produced.
RECORDS = {row[0]: row[2] for row in seqio_corpus.CORPORA}


# The reference's position classes, by their own names, plus our spelling of the
# same kinds.  A gate written on the objects would not see a `<` that had been
# flattened into an exact coordinate, because `==` between two positions in
# Biopython compares integers; this reduces both sides to the kind and every
# edge it has, so the flattening shows up as a difference here.
_KINDS = {
    "ExactPosition": "exact",
    "BeforePosition": "before",
    "AfterPosition": "after",
    "WithinPosition": "within",
    "OneOfPosition": "one_of",
    "UncertainPosition": "uncertain",
    "UnknownPosition": "unknown",
}


def _position(kind, value, left, right, choices) -> tuple:
    # `int(position)` is the value every comparison uses, and `UnknownPosition`
    # is the one kind that has none -- the reference's own `int()` raises on it,
    # so 0 here is this gate's neutral and not a claim about the position.
    try:
        value = int(value)
    except TypeError:
        value = 0
    return (kind, value, int(left), int(right), tuple(int(c) for c in choices))


def reference_location(location) -> tuple:
    """The reference's location, reduced to nested plain tuples."""
    if location is None:
        return ()
    return (
        getattr(location, "operator", "simple"),
        tuple(
            (
                _position(_KINDS[type(part.start).__name__], part.start,
                          getattr(part.start, "_left", 0),
                          getattr(part.start, "_right", 0),
                          getattr(part.start, "position_choices", ())),
                _position(_KINDS[type(part.end).__name__], part.end,
                          getattr(part.end, "_left", 0),
                          getattr(part.end, "_right", 0),
                          getattr(part.end, "position_choices", ())),
                # The reference's `None` strand is our 0, because `None` is not a
                # number a C++ struct can hold; the two have to reduce alike.
                part.strand if part.strand is not None else 0,
                part.ref or "",
            )
            for part in location.parts
        ),
    )


def our_location(location) -> tuple:
    """Our own location value, reduced to the same nested plain tuples."""
    if location is None:
        return ()
    return (
        location.operator,
        tuple(
            (
                _position(part.start.kind, part.start.value, part.start.left,
                          part.start.right, part.start.choices),
                _position(part.end.kind, part.end.value, part.end.left,
                          part.end.right, part.end.choices),
                part.strand,
                part.ref,
            )
            for part in location.parts
        ),
    )


def reference_features(format: str, path: Path) -> list:
    """Every feature of every record, reduced to what the gate compares.

    The reference emits no warnings on this corpus, and `SeqIO.parse` would
    raise under `filterwarnings = ["error"]` if it did -- which is the point of
    the corpus being generated clean.
    """
    rows = []
    with path.open() as handle:
        for record in SeqIO.parse(handle, format):
            rows.extend(
                (
                    record.id,
                    feature.type,
                    reference_location(feature.location),
                    tuple(sorted((key, tuple(values))
                                 for key, values in feature.qualifiers.items())),
                )
                for feature in record.features
            )
    return rows


def our_features(format: str, path: Path, interop: bool) -> list:
    """The same rows, from our reader -- and optionally through the conversion.

    `interop` decides whether a `Bio.SeqFeature` is built, and nothing else:
    both settings produce the same list, which the gate checks, so the two rows
    are the same work plus the object construction.
    """
    import biofasting

    index = biofasting.open_genbank(path, format)
    rows = []
    for record_id in index:
        for feature in biofasting.read_features(index, record_id):
            if interop:
                built = biofasting.to_seqfeature(feature)
                location = reference_location(built.location)
                qualifiers = built.qualifiers
            else:
                location = our_location(feature.location)
                qualifiers = feature.qualifiers_as_dict()
            rows.append(
                (
                    record_id,
                    feature.type,
                    location if feature.status == "ok" else (),
                    tuple(sorted((key, tuple(values))
                                 for key, values in qualifiers.items())),
                )
            )
    return rows


def digest_features(rows) -> str:
    """One hash over every field of every feature, so a row is all or nothing."""
    hasher = hashlib.sha256()
    for row in rows:
        for field in row:
            # The location is nested tuples of ints, so it goes in as its repr:
            # one spelling, produced by one function, on both sides.
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

    header = (f"{'corpus':18s} {'rec':>6s} {'feat':>7s} {'ref us/rec':>11s} "
              f"{'ours':>8s} {'ref/ours':>9s} {'+interop':>9s} {'gate':>7s}")
    print(header)
    print("-" * len(header))

    invalid = 0
    for group in groups:
        format, path = corpus_file(group)
        expected = reference_features(format, path)
        ours = our_features(format, path, interop=False)
        with_interop = our_features(format, path, interop=True)

        gate = (
            "OK"
            if digest_features(ours) == digest_features(expected)
            and digest_features(with_interop) == digest_features(expected)
            else "INVALID"
        )
        invalid += gate != "OK"

        count = RECORDS[group]
        features = len(expected)
        ref = timed(lambda: reference_features(format, path), args.repeat)
        got = timed(lambda: our_features(format, path, interop=False), args.repeat)
        plus = timed(lambda: our_features(format, path, interop=True), args.repeat)

        print(f"{group:18s} {count:6d} {features:7d} "
              f"{ref / count * 1e6:11.2f} {got / count * 1e6:8.2f} "
              f"{ref / got:8.2f}x {plus / count * 1e6:8.2f} {gate:>7s}")

    # The per-feature cost, as the difference between two rows that differ only
    # in feature density -- which is how the ranking pass measured the target
    # this milestone was built against.
    #
    # Per *record* and then subtracted, never per file: `genbank-1kb-bare` holds
    # 300 records where `genbank-1kb` holds 3,000, because the corpus's record
    # counts are sized to make the two files about equally long in bytes.  The
    # records themselves are the same 1 kb shape from the same seed, so the
    # difference of two per-record times is the cost of one feature and nothing
    # else -- which is exactly the 20.2 us M19 recorded.
    print()
    if BARE in RECORDS and "genbank-1kb" in groups:
        format, featured = corpus_file("genbank-1kb")
        _, bare = corpus_file(BARE)

        def marginal(run):
            """The cost of one feature: two per-record times, subtracted.

            Per record and not per file, because the bare row holds 300 records
            to the featured row's 3,000 -- the corpus sizes the counts to make
            the two files about equally long in bytes.  The records themselves
            are the same 1 kb shape from the same seed, so the difference is the
            feature and nothing else, which is exactly M19's 20.2 us.
            """
            return (
                run(format, featured, RECORDS["genbank-1kb"])
                - run(format, bare, RECORDS[BARE])
            )

        def reference_per_record(fmt, path, count):
            return timed(lambda: reference_features(fmt, path), args.repeat) / count

        def our_per_record(fmt, path, count):
            return timed(lambda: our_features(fmt, path, interop=False),
                         args.repeat) / count

        def interop_per_record(fmt, path, count):
            return timed(lambda: our_features(fmt, path, interop=True),
                         args.repeat) / count

        per_record = marginal(reference_per_record)
        ours_record = marginal(our_per_record)
        interop_record = marginal(interop_per_record)
        print(f"per feature, marginal over {BARE} "
              f"({RECORDS['genbank-1kb']} records carrying one feature each "
              f"against {RECORDS[BARE]} carrying none):")
        print(f"  reference                 {per_record * 1e6:7.2f} us")
        print(f"  ours (kernel)             {ours_record * 1e6:7.2f} us")
        print(f"  ours + to_seqfeature      {interop_record * 1e6:7.2f} us")
        print(f"  M19's recorded target     {20.2:7.2f} us, of which 4.2 is the "
              "`SeqFeature` object")
        print("  the reference's own absolute times do not reproduce on a re-run "
              "(see targets.md);")
        print("  the shape does, and this row is timed in one harness in one run.")

    print()
    print("`ref us/rec` is Bio.SeqIO.parse over the whole record, which produces the")
    print("annotations dict and the references as well as the features.  `ours`")
    print("produces the features and nothing else, so `ref/ours` is a ratio against")
    print("a reference doing more work -- the per-feature line above is the")
    print("like-for-like comparison, and it is the one to quote.")
    if invalid:
        print(f"\n{invalid} row(s) INVALID: the gate disagrees with the reference "
              "and no time from them may be quoted.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
