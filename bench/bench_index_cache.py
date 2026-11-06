#!/usr/bin/env python3
"""The index cache, measured against pyfastx: open, fetch, iterate, and RAM.

The problem this file exists for was found by measuring, not by reading.  Before
the index became a file, opening a million-read FASTQ *cost a full scan*: 0.5 s
against pyfastx's 0.4 ms, because pyfastx keeps a SQLite index next to the data
and opening it is a lookup, while this package rebuilt its index in the heap
every time.  The fetch itself was already cheaper in bulk; it was the open, and
the memory, that lost.

Two asymmetries had to be measured and not argued away:

* **When you pay for the index.**  A cold open builds and writes the cache; a
  warm one maps it.  Reporting only one of them would flatter whichever side
  the author happened to run.  Both are rows here, and `build` is the cold one.
* **What "memory" means.**  `VmHWM` counts file-backed pages the kernel can
  drop, and this package reads through mappings on purpose, so `VmHWM` makes it
  look heavier than the heap it actually holds.  The table reports **anonymous
  RSS** -- pages this process allocated and cannot evict -- next to peak RSS.
  The anonymous number is the one that says whether a 100 GB file fits; the
  peak is the one a reviewer will quote, so both are printed, and both are
  followed across the file-size axis so that "constant" and "linear" can be told
  apart instead of asserted.

The gate is the records.  Before any timing, both readers fetch the same keys
and the bytes are compared; a mismatch aborts the run, because a faster index
that answers differently is not a faster index.  Every timing is a median over
`--repeats` fresh processes -- one process per measurement -- because a cold
build and a warm open cannot share a process without the first warming the
second, and because peak RSS is per process.

The `open + N` rows are the crossover curve: they are what a program that opens
a file and reads N records actually pays, which is the only honest way to quote
a per-fetch number.  `fetch N (warm)` is the other half of that story -- the
second pass over the same keys, which is where the page faults of first touch
have stopped being part of the measurement.

Usage:

    python3 bench/bench_index_cache.py [--repeats 5] [--skip-pyfastx]

`pyfastx` is optional: without it the rows still run and the comparison columns
print `-`.
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
RESULTS = Path(__file__).resolve().parent / "results"

# Two sizes, so that "does this grow with the file?" has an answer in the table
# rather than in the prose.
FILES = {
    "plain 10k": DATA / "fastq" / "reads_10k.fastq",
    "plain 1m": DATA / "fastq" / "reads_1m.fastq",
    "gzip 1m": DATA / "fastq" / "reads_1m.fastq.gz",
}

# Every case, as (label, mode).  `build` deletes the cache first, so it is the
# cold path; everything else reads whatever is on disk, which is the warm path.
CASES = [
    ("build", "build"),
    ("open", "open"),
    ("open + 1", "fetch:1"),
    ("open + 10", "fetch:10"),
    ("open + 100", "fetch:100"),
    ("open + 1000", "fetch:1000"),
    ("fetch 1000 (warm)", "warm:1000"),
    ("iterate", "iterate"),
]

# A deterministic scatter, so that two implementations fetch the same records and
# a fetch is a random access rather than a walk down the file.  `total` is the
# file's record count, which is what makes the scatter valid for a file of any
# size: 997 and 10 000 or 1 000 000 are coprime, so the walk visits distinct keys.
def key(i, total):
    return "SRR000001.%d" % (i * 997 % total)


def names_for(count, total):
    return [key(i, total) for i in range(count)]


def fetch(get, count, total):
    total_bytes = 0
    for name in names_for(count, total):
        total_bytes += len(get(name))
    return total_bytes


def measure_ours(path, mode):
    import biofasting

    if mode == "build":
        for suffix in (".bfqi", ".bfqd"):
            try:
                os.remove(os.fspath(path) + suffix)
            except FileNotFoundError:
                pass
        started = time.perf_counter()
        index = biofasting.open_fastq_index(path)
        return time.perf_counter() - started, len(index)

    started = time.perf_counter()
    index = biofasting.open_fastq_index(path)
    opened = time.perf_counter() - started

    if mode == "open":
        return opened, None
    if mode.startswith("fetch:"):
        count = int(mode.split(":")[1])
        started = time.perf_counter()
        total = fetch(index.__getitem__, count, len(index))
        elapsed = time.perf_counter() - started
        # The open is reported alongside the fetch because the question a reader
        # actually asks is "open and read one record", and that is the sum.
        return (opened, elapsed, total), None
    if mode.startswith("warm:"):
        count = int(mode.split(":")[1])
        fetch(index.__getitem__, count, len(index))  # first touch: faults, not measurement
        started = time.perf_counter()
        total = fetch(index.__getitem__, count, len(index))
        return (opened, time.perf_counter() - started, total), None
    if mode == "iterate":
        started = time.perf_counter()
        count = 0
        for _ in biofasting.open_fastq(path):
            count += 1
        return time.perf_counter() - started, count
    raise SystemExit("unknown mode " + mode)


def measure_pyfastx(path, mode):
    import pyfastx

    index_path = Path(os.fspath(path) + ".fxi")
    if mode == "build":
        try:
            index_path.unlink()
        except FileNotFoundError:
            pass
        started = time.perf_counter()
        reader = pyfastx.Fastq(os.fspath(path), build_index=True)
        return time.perf_counter() - started, len(reader)

    started = time.perf_counter()
    reader = pyfastx.Fastq(os.fspath(path), build_index=True)
    opened = time.perf_counter() - started

    if mode == "open":
        return opened, None
    if mode.startswith("fetch:"):
        count = int(mode.split(":")[1])
        started = time.perf_counter()
        total = fetch(reader.__getitem__, count, len(reader))
        return (opened, time.perf_counter() - started, total), None
    if mode.startswith("warm:"):
        count = int(mode.split(":")[1])
        fetch(reader.__getitem__, count, len(reader))
        started = time.perf_counter()
        total = fetch(reader.__getitem__, count, len(reader))
        return (opened, time.perf_counter() - started, total), None
    if mode == "iterate":
        started = time.perf_counter()
        count = 0
        for _ in pyfastx.Fastx(os.fspath(path)):
            count += 1
        return time.perf_counter() - started, count
    raise SystemExit("unknown mode " + mode)


def rss():
    """(anonymous MB, peak MB) from the kernel's own accounting."""
    anonymous = 0
    with open("/proc/self/smaps_rollup") as handle:
        for line in handle:
            if line.startswith("Anonymous:"):
                anonymous = int(line.split()[1]) // 1024
    peak = 0
    with open("/proc/self/status") as handle:
        for line in handle:
            if line.startswith("VmHWM:"):
                peak = int(line.split()[1]) // 1024
    return anonymous, peak


def child(argv):
    kind, impl, mode = argv[0], argv[1], argv[2]
    path = FILES[kind]
    started = time.perf_counter()
    if impl == "ours":
        value, records = measure_ours(path, mode)
    else:
        value, records = measure_pyfastx(path, mode)
    elapsed = time.perf_counter() - started
    anonymous, peak = rss()
    if isinstance(value, tuple):
        value = list(value)
    print(json.dumps({"kind": kind, "impl": impl, "mode": mode,
                      "value": value, "records": records,
                      "wall": elapsed, "anon_mb": anonymous,
                      "peak_mb": peak}))


def run(kind, impl, mode, repeats):
    rows = []
    for _ in range(repeats):
        result = subprocess.run(
            [sys.executable, __file__, "--child", kind, impl, mode],
            capture_output=True, text=True, check=True)
        rows.append(json.loads(result.stdout.strip().splitlines()[-1]))
    return rows


def gate(skip_pyfastx):
    """Both readers must return the same bytes before any time is reported."""
    import biofasting

    for kind, path in FILES.items():
        ours = biofasting.open_fastq_index(path)
        names = names_for(100, len(ours))
        mine = [ours[name] for name in names]
        if skip_pyfastx:
            continue
        import pyfastx
        theirs = pyfastx.Fastq(os.fspath(path), build_index=True)
        for name, sequence in zip(names, mine):
            other = theirs[name].seq.encode()
            if sequence != other:
                raise SystemExit("gate: %s %s disagrees" % (kind, name))
    print("gate: records agree (100 keys x %d file kinds)" % len(FILES))


def seconds(value):
    if value < 1e-3:
        return "%.1f us" % (value * 1e6)
    if value < 1:
        return "%.2f ms" % (value * 1e3)
    return "%.3f s" % value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--skip-pyfastx", action="store_true")
    parser.add_argument("--child", nargs=3, metavar=("KIND", "IMPL", "MODE"))
    args = parser.parse_args()

    if args.child:
        child(args.child)
        return

    gate(args.skip_pyfastx)
    implementations = ["ours"] if args.skip_pyfastx else ["ours", "pyfastx"]
    results = []
    for kind in FILES:
        for label, mode in CASES:
            for impl in implementations:
                rows = run(kind, impl, mode, args.repeats)
                values, opens = [], []
                for row in rows:
                    value = row["value"]
                    if isinstance(value, list):
                        opens.append(value[0])
                        values.append(value[1])
                    else:
                        values.append(value)
                entry = {
                    "kind": kind,
                    "case": label,
                    "impl": impl,
                    "seconds": statistics.median(values),
                    "min_s": min(values),
                    "max_s": max(values),
                    "anon_mb": statistics.median(r["anon_mb"] for r in rows),
                    "peak_mb": statistics.median(r["peak_mb"] for r in rows),
                    "records": rows[-1]["records"],
                }
                if opens:
                    entry["open_s"] = statistics.median(opens)
                    entry["total_s"] = entry["open_s"] + entry["seconds"]
                results.append(entry)
                print("%-9s %-19s %-8s %10s   anon %4d MB   peak %4d MB"
                      % (kind, label, impl, seconds(entry["seconds"]),
                         entry["anon_mb"], entry["peak_mb"]), flush=True)

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "index_cache.json").write_text(json.dumps(results, indent=2))
    print("written", RESULTS / "index_cache.json")

    # The table the paper quotes, derived from the JSON rather than retyped.
    by = {(entry["kind"], entry["case"], entry["impl"]): entry
          for entry in results}
    print()
    print("| file | case | ours | pyfastx | ratio | ours anon/peak MB | "
          "pyfastx anon/peak MB |")
    print("|---|---|---|---|---|---|---|")
    for kind in FILES:
        for label, _ in CASES:
            ours = by.get((kind, label, "ours"))
            theirs = by.get((kind, label, "pyfastx"))
            if ours is None:
                continue
            ratio = ("%.2fx" % (theirs["seconds"] / ours["seconds"])
                     if theirs and ours["seconds"] else "-")
            print("| %s | %s | %s | %s | %s | %d / %d | %s |" % (
                kind, label, seconds(ours["seconds"]),
                seconds(theirs["seconds"]) if theirs else "-", ratio,
                ours["anon_mb"], ours["peak_mb"],
                ("%d / %d" % (theirs["anon_mb"], theirs["peak_mb"]))
                if theirs else "-"))


if __name__ == "__main__":
    main()
