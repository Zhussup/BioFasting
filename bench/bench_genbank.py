#!/usr/bin/env python3
"""The flat-file kernel, measured against the reference it was built to beat.

The sixth ranking pass (`bench/rank_seqio.py`, written up in `bench/targets.md`)
produced the target this kernel exists for: on a 1 kb GenBank record with no
features the reference spends 28.4 us, where the bytes cost 1.5 and the object
tree 1.0.  That is a target and not a result, so this script is the other half --
the delivered reader, on the same corpus, under the same reference and the same
digest gate as every row in `bench/run.py`.

Three rows per workload, and the middle one is why the first is not the whole
story:

* **reference** -- ``SeqIO.parse``, which produces the features and the
  annotations as well.  This is the thing a caller would otherwise use.
* **rewrite** -- the smallest pure-Python parser there is: split on ``//``, take
  the id, strip the ``ORIGIN``/``SQ`` block.  It reproduces ids and sequences and
  nothing else, so its ratio is a floor under what any Python rewrite could cost
  and never a replacement.  Kept because it is the number that says whether the
  win needed C at all.
* **biofasting** -- ``open_genbank``: id, name, description and sequence, built
  in one pass over the mapped file.

What the biofasting row does **not** reproduce is the FEATURES table and the
annotations dict, so its ratio is against a reference that produces more than it
does, exactly as the rewrite's is.  The script prints that next to the number
rather than in a footnote, because a speedup quoted without its scope is a
speedup about an unknown program.

Every row is digest-gated first: the reference's four fields and ours are hashed
and compared before any time is quoted, and a row that disagrees is printed
``INVALID`` and never ranked.

    python3 bench/bench_genbank.py
    python3 bench/bench_genbank.py --only genbank-1kb-bare --repeat 9
    python3 bench/bench_genbank.py --list
"""

from __future__ import annotations

import argparse
import hashlib
import io
import statistics
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import seqio_corpus  # noqa: E402
from Bio import SeqIO  # noqa: E402

# The corpus is generated, not stored, but the reference needs a file it can
# open and this needs a file it can map.  The cache goes where `bench/run.py`
# puts its own derived files -- `build/`, which is disposable and never tracked
# -- and deliberately **not** under `bench/data/`: `tests/test_corpus.py` walks
# that directory and requires the walk to equal the generator's manifest, so a
# file written there by anything else is a failure by design.  It was written
# there first, and that is how the rule was found.
CORPUS_DIR = HERE.parent / "build" / "bench-cache" / "seqio"

# The extension is not cosmetic: a reader that decided the format from the name
# would be a reader that works on this corpus and not on a file called `.gbk`.
EXTENSIONS = {"genbank": "gb", "embl": "embl", "swiss": "dat"}


def timed(fn, repeats: int) -> float:
    """Median seconds over `repeats` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def fields_of(record) -> tuple[str, str, str, str]:
    """The four fields under test, from a reference SeqRecord."""
    return record.id, record.name, record.description, str(record.seq)


def digest(rows) -> str:
    """One hash over every field of every record, so a row is all or nothing.

    Per-record comparison would say *which* record differs and this says only
    whether any does, which is the wrong trade for a benchmark: the row is
    either quotable or it is not, and a mismatch is a bug to be found by the
    differential tests in `tests/test_genbank.py`, not by reading a table.
    """
    hasher = hashlib.sha256()
    for row in rows:
        for field in row:
            if isinstance(field, str):
                field = field.encode()
            hasher.update(len(field).to_bytes(8, "little"))
            hasher.update(field)
    return hasher.hexdigest()


def corpus_file(group: str) -> tuple[str, Path]:
    """The (format, path) of one corpus row, written if it is not already there."""
    spec = {row[0]: row for row in seqio_corpus.CORPORA}
    _, format, count, length, per_kb, annot = spec[group]
    text = seqio_corpus.build(format, count, length, per_kb, annot)
    CORPUS_DIR.mkdir(parents=True, exist_ok=True)
    path = CORPUS_DIR / f"{group}.{EXTENSIONS[format]}"
    if not path.exists() or path.read_text() != text:
        path.write_text(text, newline="")
    return format, path


def reference_rows(format: str, path: Path) -> list:
    with path.open() as handle:
        return [fields_of(record) for record in SeqIO.parse(handle, format)]


def rewrite_rows(format: str, path: Path) -> list:
    """The pure-Python floor, imported rather than copied.

    `rank_seqio.minimal_parse` is the function that produced the ranking pass's
    rewrite column, and a second copy of it here would be free to drift away
    from the number it is being compared against.
    """
    from rank_seqio import minimal_parse

    return [(ident, ident, "", sequence)
            for ident, sequence in minimal_parse(format, path.read_text())]


def biofasting_rows(format: str, path: Path) -> list:
    import biofasting

    return [
        (record.id, record.name, record.description, record.sequence.decode())
        for record in biofasting.read_genbank(path, format=format)
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", action="append", metavar="GROUP",
                        help="measure only this corpus row; repeatable")
    parser.add_argument("--repeat", type=int, default=5,
                        help="timed passes per implementation (default 5)")
    parser.add_argument("--list", action="store_true",
                        help="list the corpus rows and exit")
    args = parser.parse_args(argv)

    groups = [row[0] for row in seqio_corpus.CORPORA]
    if args.list:
        for group in groups:
            print(group)
        return 0
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

    # The build identifies itself, because a time without its compiler, flags and
    # target architecture is not a measurement.  The key names are the ones
    # `biofasting.build_info()` actually publishes -- `extra_cxx_flags` and
    # `arch`, not the `cxx_flags`/`target_arch` a reader might expect.
    info = biofasting.build_info()
    print(f"{info['version']}  {info['compiler']} c++{info['cxx_standard']}  "
          f"flags={info['extra_cxx_flags']}  {info['arch']}  "
          f"seqops={biofasting.seqops_level()}  cpu={biofasting.best_cpu_level()}")
    print(f"Python {sys.version.split()[0]}, median of {args.repeat} passes "
          f"after one warm-up\n")

    header = (f"{'corpus':18s} {'rec':>6s} {'MB':>7s} {'ref us/rec':>11s} "
              f"{'rewrite':>9s} {'ours':>8s} {'ref/ours':>9s} {'MB/s':>8s} "
              f"{'gate':>7s}")
    print(header)
    print("-" * len(header))

    invalid = 0
    for group in groups:
        format, path = corpus_file(group)
        size = path.stat().st_size
        reference = reference_rows(format, path)
        rewrite = rewrite_rows(format, path)
        ours = biofasting_rows(format, path)

        # The rewrite produces ids and sequences only, so it is gated on those
        # two fields; ours produces all four and is gated on all four.  A gate
        # that compared the rewrite's description against the reference's would
        # be a permanent false INVALID and would make the column meaningless.
        ok_rewrite = digest([row[:2] + ("", row[3]) for row in rewrite]) == digest(
            [row[:2] + ("", row[3]) for row in reference])
        ok_ours = digest(ours) == digest(reference)
        gate = "OK" if (ok_rewrite and ok_ours) else "INVALID"
        invalid += gate != "OK"

        count = len(reference)
        ref = timed(lambda: reference_rows(format, path), args.repeat)
        rew = timed(lambda: rewrite_rows(format, path), args.repeat)
        got = timed(lambda: biofasting_rows(format, path), args.repeat)

        print(f"{group:18s} {count:6d} {size / 1e6:7.3f} "
              f"{ref / count * 1e6:11.2f} {rew / count * 1e6:9.2f} "
              f"{got / count * 1e6:8.2f} {ref / got:8.2f}x "
              f"{size / got / 1e6:8.1f} {gate:>7s}")

    print()
    print("`ref us/rec` is Bio.SeqIO.parse, which produces the FEATURES table and")
    print("the nine annotations as well.  `ours` produces id, name, description and")
    print("the sequence, and nothing else -- so `ref/ours` is a ratio against a")
    print("reference that does more work, and is not a like-for-like speedup.")
    print("`rewrite` is the same four fields minus the name and the description,")
    print("in pure Python: it is a floor under any Python rewrite, never a")
    print("replacement, and it is here to say whether C was needed at all.")
    if invalid:
        print(f"\n{invalid} row(s) INVALID: the gate disagrees with the reference "
              "and no time from them may be quoted.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
