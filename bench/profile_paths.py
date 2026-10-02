#!/usr/bin/env python3
"""Step 0.4 -- where the time actually goes on the slowest measured paths.

The runner in ``run.py`` says *how much* slower Biopython is.  This says *why*,
by naming the functions that hold the time and giving each one a share.

Two views, because they answer different questions:

* ``cProfile`` (built in, used here) counts every Python-level call and
  attributes time to frames.  It answers "how much of this path is interpreter
  work that a C core removes entirely?", which is the question that decides what
  Phase 1 should implement.  Its absolute seconds are inflated -- instrumentation
  costs roughly 2-4x on call-heavy code -- so read the **shares**, not the
  seconds, and never compare a profiled time to a number from ``run.py``.

* ``py-spy`` samples the real process with no instrumentation, so it sees the
  true wall-clock split including time inside C extensions, which cProfile
  charges to whichever Python frame happened to call them.  Point it at this
  script rather than reimplementing the workloads::

      .venv/bin/py-spy record --format raw -o /tmp/prof.raw --rate 200 \\
          -- .venv/bin/python bench/profile_paths.py --workload fastq:biopython-seqio
      .venv/bin/python bench/profile_paths.py --pyspy /tmp/prof.raw \\
          --pyspy-mark biopython-seqio

Usage::

    .venv/bin/python bench/profile_paths.py --list
    .venv/bin/python bench/profile_paths.py --workload fastq:biopython-seqio --repeat 2
    .venv/bin/python bench/profile_paths.py --all --repeat 1 --top 12
"""

from __future__ import annotations

import argparse
import cProfile
import importlib
import io
import pstats
import sys
import time
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import run as bench  # noqa: E402  (the runner, reused so profiled work == measured work)

HARNESS = str(Path(__file__).resolve())
BIO_MARKER = "site-packages/Bio/"
DRIVER_FUNCS = {"consume", "digest"}

# For most paths the driver frame is pure overhead, identical in every row.  For
# these two it is not, and saying "harness" without qualification would be
# wrong: in one it is the implementation, in the other it hides the C library.
# Keyed by implementation label, which is the part after the "prefix:" in a path
# name.
DRIVER_NOTES = {
    "naive-python": "this loop IS the parser -- the pure-Python floor it measures "
                    "is exactly this bytecode, not overhead",
    "pysam-htslib": "this loop calls into htslib; cProfile charges that C time to "
                    "the calling frame, so it cannot see inside -- use py-spy",
}

DEFAULT_DATA = HERE / "data"


# --------------------------------------------------------------------------- #
# The paths worth profiling: the slow rows from bench/results/latest.json
# --------------------------------------------------------------------------- #


def build_paths(data: Path, backends: dict, quick: bool) -> dict[str, callable]:
    """name -> a zero-argument callable doing one full pass of that path.

    Every path is built from the **runner's own workload factories** rather than
    reimplemented here.  That is the point of this file: the code being profiled
    is the code that was measured, so a hotspot found here is a hotspot behind
    a row of ``bench/results/latest.json`` and not an artifact of a lookalike
    loop that happens to do something similar.

    ``make_fasta_random_impls`` builds its index once per implementation and
    caches it inside the factory, exactly as it does under the runner, so
    repeated passes measure access, not index construction.
    """
    fastq = data / ("fastq/reads_10k.fastq" if quick else "fastq/reads_1m.fastq")
    gz = data / "fastq/reads_1m.fastq.gz"
    fasta = data / ("fasta/genome_1mb.fasta" if quick else "fasta/genome.fasta")

    paths: dict[str, callable] = {}

    def add(prefix: str, impls: list) -> None:
        for impl in impls:
            paths[f"{prefix}:{impl.label}"] = partial(bench.consume, impl)

    add("fastq", bench.make_fastq_impls(fastq, backends))
    if gz.exists():
        add("fastq-gz", bench.make_fastq_impls(gz, backends)
            + bench.make_decompression_impls(gz))
    add("fasta", bench.make_fasta_impls(fasta, backends))
    if not quick:
        add("fasta-random", bench.make_fasta_random_impls(fasta, backends, 100, 150))

    reads = []
    for rec in backends["Bio"].SeqIO.parse(str(fastq), "fastq"):
        reads.append(rec)
        if len(reads) >= bench._OPS_READS:
            break
    revcomp, gc = bench.make_op_impls(reads, backends)
    add("ops", revcomp + gc)
    return paths


# --------------------------------------------------------------------------- #
# Reading a profile: shares, not seconds
# --------------------------------------------------------------------------- #


def analyse(profile: cProfile.Profile, elapsed: float, top: int) -> dict:
    """Split profiled time into Biopython's own frames, driver frames, and rest.

    ``cProfile`` cannot tell Python bytecode from time inside a C extension: it
    charges both to whichever Python frame was on the stack.  So the two
    aggregates here are chosen to be unambiguous anyway:

    * **Biopython frames** -- time in ``site-packages/Bio/**``.  For a pure
      Python library that is interpreter work, which is the part a C core
      deletes outright.  For pyfaidx, zlib and pysam it is ~0%, because those
      keep their work in C; a low number there is not a compliment.
    * **driver frames** -- this file's loop and ``run.py``'s ``consume``/
      ``digest``.  Identical work in every row, so it is worth separating out
      rather than letting it flatter one implementation over another.

    Everything else is the benchmark's own implementation frames (``run.py``'s
    factories, where the naive parser and the pysam loop live) plus genuine C
    time.  The hotspot table below says which is which per row.
    """
    stream = io.StringIO()
    stats = pstats.Stats(profile, stream=stream).sort_stats("tottime")
    stats.print_stats(top)

    rows = []
    bio_time = driver_time = 0.0
    for (filename, _line, func), (_cc, _nc, tt, _ct, _callers) in stats.stats.items():
        if filename == HARNESS or (Path(filename).name == "run.py"
                                   and func in DRIVER_FUNCS):
            driver_time += tt
        elif BIO_MARKER in filename:
            bio_time += tt
        rows.append({"file": filename, "func": func, "tottime": tt, "cumtime": _ct})
    rows.sort(key=lambda r: r["tottime"], reverse=True)

    scale = elapsed or 1.0
    return {
        "bio_share": bio_time / scale,
        "driver_share": driver_time / scale,
        "rows": rows[:top],
        "report": stream.getvalue(),
        "functions_seen": len(stats.stats),
    }


def short_path(filename: str) -> str:
    for marker in ("site-packages/", "lib/python3.13/"):
        if marker in filename:
            return filename.split(marker, 1)[1]
    return filename


def summarise_raw(path: Path, top: int, mark: str | None = None) -> int:
    """Summarise a py-spy ``--format raw`` capture as a self-time table.

    ``raw`` is py-spy's collapsed-stack output: one sample per line, stacks
    joined by ``;`` and a trailing count.  Aggregating the *leaf* frame gives
    self time as the sampler saw it, with no instrumentation and no attribution
    games -- which is the point of running py-spy alongside cProfile.

    ``mark`` keeps only samples whose stack contains that substring, which is
    how the module-import prologue (unavoidable in any short capture) is
    excluded from the shares.
    """
    counts: dict[str, int] = {}
    total = 0
    kept = 0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        stack, _, count = line.rpartition(" ")
        try:
            count = int(count)
        except ValueError:
            continue
        total += count
        if mark and mark not in stack:
            continue
        kept += count
        leaf = stack.rsplit(";", 1)[-1].strip()
        counts[leaf] = counts.get(leaf, 0) + count
    if not total:
        print(f"{path}: no samples", file=sys.stderr)
        return 2

    print(f"{path}: {total:,} samples"
          + (f", {kept:,} matching {mark!r}" if mark else ""))
    denominator = kept or total
    for leaf, count in sorted(counts.items(), key=lambda kv: -kv[1])[:top]:
        print(f"   {count:6,}  {count / denominator * 100:5.1f}%  {leaf}")
    print("\nSelf time: a sample is charged to the frame that was executing, not to")
    print("its callers.  For a C-backed row the leaf is the Python frame that called")
    print("into C -- the last frame the unwinder can name -- so a leaf share near")
    print("100% there means the time is inside the C call, not that the caller is")
    print("slow.")
    return 0


def main(argv: list[str] | None = None) -> int:
    # Same reason as in run.py: a redirected run should show progress per path
    # instead of nothing until the end.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):
        pass
    parser = argparse.ArgumentParser(
        prog="profile_paths.py", description="Profile the slowest benchmark paths")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--workload", action="append", default=None)
    parser.add_argument("--all", action="store_true", help="profile every path")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--repeat", type=int, default=1, help="passes to profile")
    parser.add_argument("--top", type=int, default=10, help="hotspots to name")
    parser.add_argument("--pyspy", type=Path, default=None,
                        help="summarise a py-spy --format raw capture instead of profiling")
    parser.add_argument("--pyspy-mark", default=None,
                        help="only count samples whose stack contains this substring")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args(argv)

    # Summarising a capture needs neither the corpus nor the backends.
    if args.pyspy:
        return summarise_raw(args.pyspy, args.top, args.pyspy_mark)

    backends: dict = {}
    for name in ("Bio", "pysam", "pyfaidx"):
        module, _error = bench.try_import(name)
        if module is not None:
            backends[name] = module
    if "Bio" not in backends:
        print("error: Biopython is required", file=sys.stderr)
        return 2
    for submodule in ("Bio.SeqIO", "Bio.SeqIO.QualityIO", "Bio.SeqUtils"):
        importlib.import_module(submodule)

    paths = build_paths(args.data, backends, args.quick)
    if args.list:
        for name in paths:
            print(name)
        return 0

    wanted = list(paths) if (args.all or not args.workload) else args.workload
    unknown = [w for w in wanted if w not in paths]
    if unknown:
        print(f"error: unknown workload(s) {unknown}; --list shows them", file=sys.stderr)
        return 2

    print(f"cProfile, {args.repeat} pass(es) per path.  Shares are meaningful; the")
    print("seconds are inflated by instrumentation and are not the runner's numbers.\n")

    for name in wanted:
        path = paths[name]
        path()  # warm the page cache and let lazy imports settle, untimed
        profile = cProfile.Profile()
        start = time.perf_counter()
        profile.enable()
        for _ in range(args.repeat):
            path()
        profile.disable()
        elapsed = time.perf_counter() - start

        result = analyse(profile, elapsed, args.top)
        print(f"== {name}   {elapsed:.3f}s profiled, {result['functions_seen']} functions")
        print(f"   Biopython frames: {result['bio_share'] * 100:5.1f}%  "
              f"<- interpreter work a C core deletes")
        print(f"   driver frames:    {result['driver_share'] * 100:5.1f}%  "
              f"(identical in every row)")
        driver_note = DRIVER_NOTES.get(name.split(":", 1)[-1])
        if driver_note:
            print(f"   NB: {driver_note}")
        print(f"   {'tottime':>9}  {'share':>6}  hotspot")
        for row in result["rows"]:
            print(f"   {row['tottime']:8.3f}s  {row['tottime'] / elapsed * 100:5.1f}%  "
                  f"{short_path(row['file'])}:{row['func']}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
