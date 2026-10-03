#!/usr/bin/env python3
"""The zero-copy claim, measured where it can actually be measured.

`run.py` times a *parse*.  Its `consume()` does the least work an
implementation can legitimately be asked to do with each record -- `len()`, and
one uniform call -- which is exactly the right model for comparing parsers and
exactly the wrong model for the grid.  A `FastaGrid.bases` array is described
by two integer fields, so timing it inside the parse harness measures nothing:
the bytes are never touched, so the mapped file is never even read, and the row
would win by not doing the work rather than by doing it better.  A benchmark
that rewards an implementation for doing nothing is worse than no benchmark.

So the array is measured here, on the three things a caller actually buys:

1. **Getting a chromosome into a numpy array.**  Ours is a pointer and two
   numbers over the mapped file; every other way has to copy 40 Mbp first.
2. **An operation over that array.**  A vectorised numpy pass against the
   library functions and against `str.count`, on the same 40 Mbp.
3. **Resident memory.**  A view does not grow the process, and the number that
   proves it is the process's own RSS, not the array's reported `nbytes`.

Every row is gated on the bytes being right before its time is quoted: the
grid's array is compared, newline bytes removed, against the sequence
`SeqIO.parse` produces.  A row whose content does not match is printed as
INVALID and its time is not reported as a speed-up.

Also measured, because it is the honest answer to "where is the FASTQ grid?":
the width of the corpus reads' headers.  A grid needs every record to have the
same stride; `run.py` measures those reads at 1M records, and this file says
why the second array -- the one that would make them free -- cannot exist.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import statistics
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_DATA = HERE / "data"
DEFAULT_RESULTS = HERE / "results"

CHR = "chr1"


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #


def try_import(module: str):
    """Import an optional backend, returning ``(module_or_None, error)``."""
    import importlib

    try:
        return importlib.import_module(module), None
    except Exception as exc:  # pragma: no cover - environment dependent
        return None, f"{type(exc).__name__}: {exc}"


def rss_kb() -> int:
    """This process's resident set size, in KiB.

    `/proc/self/status` rather than `resource.getrusage`, because the field
    that matters here is what the kernel currently has mapped in -- a peak
    would never come back down and so could not show that a view costs nothing.
    """
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1])
    raise RuntimeError("VmRSS not found in /proc/self/status")


def timed(fn, repeats: int) -> list[float]:
    """`repeats` samples of `fn`, one untimed pass first to warm the caches."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return sorted(samples)


def gc_collected() -> None:
    """Drop garbage before a memory measurement, so the delta is the subject."""
    gc.collect()


# --------------------------------------------------------------------------- #
# Claim 1: getting a chromosome as an array
# --------------------------------------------------------------------------- #


def index_builders(path: Path, backends: dict):
    """Ways to make a record addressable by name, each measured once.

    These are the one-off costs, and they are separated from the access costs
    below on purpose: an implementation that is asked to build an index every
    time it answers a query is not being compared with one that is not.  The
    first version of this file made exactly that mistake in the reference's
    favour -- our row re-opened and re-scanned the file inside the timed
    region while `SeqIO.index` had been built at setup -- and it reported the
    copying implementation as twice as fast.  The fix is not a footnote; it is
    the difference between a measurement and a rigged one.
    """
    entries = []
    biofasting = backends.get("biofasting")
    Bio = backends.get("Bio")
    pyfaidx = backends.get("pyfaidx")

    if biofasting is not None:
        entries.append(("biofasting-index", lambda: biofasting.open_fasta(path),
                        "mmap + one scan for the offset index"))
    if Bio is not None:
        from Bio import SeqIO

        entries.append(("biopython-index", lambda: SeqIO.index(str(path), "fasta"),
                        "scans the file keeping whole-record offsets (it has no other kind)"))
    if pyfaidx is not None:
        # pyfaidx keeps its index in a `.fai` next to the file, and this corpus
        # already has one (written earlier by pyfaidx itself).  Its build row
        # is therefore the cost of *reading* an index, not of making one, and
        # the note says which so the row is not read as an 18,000× win.
        fai = path.with_suffix(path.suffix + ".fai")
        entries.append(("pyfaidx-index", lambda: pyfaidx.Fasta(str(path)),
                        f"reads the .fai on disk ({'present' if fai.exists() else 'absent'}), "
                        "no FASTA scan"))
    return entries


def array_getters(path: Path, backends: dict, handles: dict):
    """Ways to end up holding chr1's bases in a numpy array, index already held.

    Each returns ``(label, getter, note)``.  Every implementation starts from a
    handle built before this point and kept alive, which is how an application
    uses one: the question is what an *access* costs once the file is
    addressable.  The build cost is the separate table above.
    """
    import numpy as np

    getters = []

    ours = handles.get("biofasting-index")
    if ours is not None:
        def via_grid():
            # The array *is* a window on the existing mapping: two integers,
            # and the record is addressable at numpy speed.  Nothing is copied.
            return ours.grid(CHR)[0]

        getters.append((
            "biofasting-grid",
            via_grid,
            "this project: the mapped file's own bytes, strided, no copy",
        ))

    theirs = handles.get("biopython-index")
    if theirs is not None:
        def seqio_index():
            # `SeqIO.index` stores whole-record offsets, so this re-reads and
            # re-parses all 40 Mbp of chr1 before copying them into an array --
            # the Phase 0 finding this project starts from, as a memory cost.
            return np.frombuffer(str(theirs[CHR].seq).encode("ascii"), dtype=np.uint8)

        getters.append((
            "biopython-index",
            seqio_index,
            "holds an index, yet each access re-parses the record, then copies it",
        ))

    Bio = backends.get("Bio")
    if Bio is not None:
        from Bio import SeqIO

        def seqio_parse():
            for record in SeqIO.parse(str(path), "fasta"):
                if record.id == CHR:
                    return np.frombuffer(str(record.seq).encode("ascii"), dtype=np.uint8)
            raise KeyError(CHR)

        getters.append((
            "biopython-parse",
            seqio_parse,
            "no index at all: reads and parses the whole file, then copies",
        ))

    fa = handles.get("pyfaidx-index")
    if fa is not None:
        def via_pyfaidx():
            return np.frombuffer(fa[CHR][:].seq.encode("ascii"), dtype=np.uint8)

        getters.append((
            "pyfaidx",
            via_pyfaidx,
            "slice by offset is cheap, but it hands back a str, which is copied",
        ))

    return getters


def reference_sequence(path: Path) -> bytes:
    """chr1's bases, newline-free, straight from `SeqIO.parse`.

    The correctness gate for every array row: whatever an implementation hands
    back has to *be* this sequence.  A view that is one byte off is a wrong
    answer that is fast, which is the failure mode this project exists to
    avoid, and it is the failure mode a stride makes easy.
    """
    from Bio import SeqIO

    for record in SeqIO.parse(str(path), "fasta"):
        if record.id == CHR:
            return str(record.seq).encode("ascii")
    raise KeyError(CHR)


def newline_free(array) -> bytes:
    """The array's bytes with the line endings dropped."""
    import numpy as np

    flat = np.asarray(array, dtype=np.uint8).ravel()
    return flat[flat != 0x0A].tobytes()


def _newlines_in(array) -> int:
    """How many line endings the array's buffer holds."""
    import numpy as np

    return int(np.count_nonzero(np.asarray(array, dtype=np.uint8) == 0x0A))


def covers_reference(array, reference: bytes) -> tuple[bool, int, int]:
    """Is the array's content the reference sequence -- and how much of it?

    Returns ``(ok, covered, tail)``.  A grid array holds whole lines, so it is
    allowed to be short by less than one line: chr1's last line is 40 bases
    where every other is 60, and a `(lines, line_bases)` view cannot describe a
    short line, so the tail is simply not in the array.  That is a documented
    property of the view, not a wrong answer -- and it stops being a documented
    property and starts being a bug the moment the covered part is not a
    prefix, or the tail is a whole line or more, which is what this checks.

    Every copying implementation covers the whole chromosome, so for them
    ``tail == 0`` is what the gate demands.
    """
    plain = newline_free(array)
    covered = len(plain)
    tail = len(reference) - covered
    stride = int(getattr(array, "strides", (1,))[0]) or 1
    ok = plain == reference[:covered] and 0 <= tail < stride
    return ok, covered, tail


# --------------------------------------------------------------------------- #
# Claim 3: what a view costs the process
# --------------------------------------------------------------------------- #


def memory_row(label: str, getter, repeat_hold: int = 1) -> dict:
    """RSS growth while holding the result, over a settled process.

    The subject is kept alive across the measurement and only then dropped, so
    that the number is "what holding this costs", not "what allocating and
    freeing it cost".  For a mapping the answer should be near zero: the pages
    belong to the page cache and to every other process that maps the file, and
    nothing was allocated to hold them.
    """
    gc_collected()
    before = rss_kb()
    held = [getter() for _ in range(repeat_hold)]
    after = rss_kb()
    nbytes = sum(getattr(item, "nbytes", 0) for item in held)
    assert held, label
    del held
    gc_collected()
    return {
        "impl": label,
        "rss_before_kb": before,
        "rss_after_kb": after,
        "rss_delta_kb": after - before,
        "array_nbytes": nbytes,
    }


# --------------------------------------------------------------------------- #
# Claim 2: an operation over the array
# --------------------------------------------------------------------------- #


def op_rows(sequence: bytes, array, backends: dict, repeats: int) -> list[dict]:
    """GC counting over chr1: vectorised, library, and Python.

    All three must agree on the count before any of them is quoted.  The
    ambiguous-base question is settled by measuring it rather than assuming:
    `gc_fraction` defaults to removing ambiguous bases from the denominator
    while `str.count` does not, so if chr1 contains any `N` the two disagree
    about the fraction for a reason that has nothing to do with speed, and the
    row says so instead of reporting a number that looks like a win.
    """
    import numpy as np

    gc_bytes = (b"G", b"C")
    expected = sum(sequence.count(byte) for byte in gc_bytes)
    ambiguous = sum(sequence.count(base.encode()) for base in "NRYKMSWBDHV")
    # The grid array holds whole lines and chr1's last line is short, so the
    # array covers a prefix of the chromosome.  A count over the array is a
    # count over that prefix, and comparing it against the whole-chromosome
    # count would flag a correct answer as INVALID -- which is what the first
    # version of this file did.  The gate is against the same span, and the
    # uncovered tail is carried in the row so the difference is visible
    # instead of hidden.
    covered = int(getattr(array, "size", 0)) - _newlines_in(array)
    expected_prefix = sum(sequence[:covered].count(byte) for byte in gc_bytes)
    rows = []

    def measured(label: str, fn, note: str, value_matches: bool) -> None:
        samples = timed(fn, repeats)
        median = statistics.median(samples)
        rows.append({
            "impl": label,
            "median_s": round(median, 6),
            "min_s": round(samples[0], 6),
            "mib_per_s": round(len(sequence) / median / 2**20, 1),
            "valid": value_matches,
            "note": note,
        })

    def numpy_count():
        return int(np.count_nonzero((array == ord("G")) | (array == ord("C"))))

    def python_count():
        return sum(sequence.count(byte) for byte in gc_bytes)

    measured("biofasting-numpy", numpy_count,
             f"np.count_nonzero over the mapped, strided array; the last "
             f"{len(sequence) - covered} bases are outside it",
             numpy_count() == expected_prefix)
    measured("python-strcount", python_count,
             "seq.count x2 on the decoded str", python_count() == expected)

    Bio = backends.get("Bio")
    if Bio is not None:
        from Bio.Seq import Seq
        from Bio.SeqUtils import gc_fraction

        seq = Seq(sequence)

        def biopython_gc():
            return gc_fraction(seq)

        # `gc_fraction` defaults to `ambiguous="remove"`, so its denominator is
        # the unambiguous bases, not the sequence length -- on this chromosome
        # that is a difference of 2000 bases, which is far larger than any
        # rounding tolerance.  Comparing it against `GC/len` would flag the
        # library as wrong for being right; comparing it against
        # `GC/(len - ambiguous)` is the same question asked twice.
        unambiguous = sum(sequence.count(base.encode()) for base in "ACGT")
        agreed = abs(gc_fraction(seq) - expected / unambiguous) < 1e-9
        measured("biopython-gc_fraction", biopython_gc,
                 "SeqUtils.gc_fraction; ambiguous bases out of the denominator",
                 agreed)

    for row in rows:
        row["gc_count"] = expected
        row["ambiguous_bases"] = ambiguous
        row["counted_over"] = covered if row["impl"] == "biofasting-numpy" else len(sequence)
    return rows


# --------------------------------------------------------------------------- #
# The FASTQ finding: why the reads have no grid
# --------------------------------------------------------------------------- #


def fastq_header_widths(path: Path, backends: dict) -> dict:
    """How wide the reads' headers are, from the record's '@' to its first base.

    A FASTQ grid needs every record to be the same number of bytes, because
    that is what turns a file into a strided array.  The corpus writes constant
    150 bp sequences with headers that are *not* constant, so the grid refuses
    -- correctly -- and this measures the cause of the refusal instead of
    asserting it.

    The first version of this function scanned for lines starting with `@`.
    That is the same mistake `gen_data.py` was caught making in M7: a quality
    line can start with `@` (Phred+33 quality 31), and this corpus has plenty.
    It reported 1,039,211 records for a 1,000,000-record file and a phantom
    150-byte "header width".  The records are counted by a parser here, which
    is the only way to count them.
    """
    widths: dict[int, int] = {}
    total = 0
    Bio = backends.get("Bio")
    if Bio is not None:
        from Bio.SeqIO.QualityIO import FastqGeneralIterator

        with open(path, newline=None) as handle:
            for title, _sequence, _quality in FastqGeneralIterator(handle):
                # +1 for the '@' that the title does not carry: the width that
                # matters is the record's stride from '@' to the first base.
                width = len(title) + 1
                widths[width] = widths.get(width, 0) + 1
                total += 1
    else:  # pragma: no cover - the reference is required for the corpus anyway
        import biofasting

        for title, _sequence, _quality in biofasting.open_fastq(path):
            width = len(title) + 1
            widths[width] = widths.get(width, 0) + 1
            total += 1
    if not widths:  # pragma: no cover - an empty file has no widths to report
        return {"records": 0, "widths": {}, "distinct_widths": 0, "min": 0, "max": 0}
    return {
        "records": total,
        "widths": dict(sorted(widths.items())),
        "distinct_widths": len(widths),
        "min": min(widths),
        "max": max(widths),
    }


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def main(argv: list[str] | None = None) -> int:
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):
        pass

    parser = argparse.ArgumentParser(
        prog="grid.py", description="the zero-copy view, measured honestly")
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--fasta", default="fasta/genome.fasta")
    parser.add_argument("--fastq", default="fastq/reads_1m.fastq")
    parser.add_argument("--json", type=Path, default=None,
                        help="where to write results (default: bench/results/grid.json)")
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args(argv)

    backends: dict = {}
    versions: dict = {}
    for name in ("Bio", "numpy", "pyfaidx", "biofasting"):
        module, error = try_import(name)
        if module is not None:
            backends[name] = module
            versions[name] = getattr(module, "__version__", "?")
        else:
            print(f"note: backend {name} unavailable ({error})", file=sys.stderr)
    if "numpy" not in backends:
        print("error: numpy is required to compare anything with a view", file=sys.stderr)
        return 2

    fasta_path = args.data / args.fasta
    if not fasta_path.exists():
        print(f"error: {fasta_path} not found (python bench/gen_data.py)", file=sys.stderr)
        return 2

    print(f"corpus: {fasta_path}   repeat: {args.repeat}   "
          f"backends: {', '.join(f'{k} {v}' for k, v in versions.items())}")

    reference = reference_sequence(fasta_path)
    print(f"\n== making the file addressable, once  ({fasta_path.name})")

    # ---- the one-off cost --------------------------------------------------
    builders = index_builders(fasta_path, backends)
    handles: dict = {}
    builds = []
    for label, builder, note in builders:
        handle = builder()
        samples = timed(builder, args.repeat)
        median = statistics.median(samples)
        handles[label] = handle
        builds.append({
            "impl": label,
            "median_s": round(median, 6),
            "min_s": round(samples[0], 6),
            "mib_per_s": round(fasta_path.stat().st_size / median / 2**20, 1),
            "note": note,
        })
        print(f"   {label:20} {median*1000:9.1f} ms  "
              f"{fasta_path.stat().st_size/median/2**20:8.1f} MiB/s  {note}")
    print("   (paid once per file; the access costs below all start from here)")

    print(f"\n== chr1 as an array, index already held  ({len(reference):,} bases)")

    # ---- claim 1: the array ------------------------------------------------
    getters = array_getters(fasta_path, backends, handles)
    arrays: dict[str, dict] = {}
    for label, getter, note in getters:
        value = getter()
        nbytes = int(getattr(value, "nbytes", 0))
        shape = getattr(value, "shape", None)
        content_ok, covered, tail = covers_reference(value, reference)
        samples = timed(getter, args.repeat)
        median = statistics.median(samples)
        arrays[label] = {
            "impl": label,
            "median_s": round(median, 6),
            "min_s": round(samples[0], 6),
            "us_per_access": round(median * 1e6, 2),
            "nbytes": nbytes,
            "shape": list(shape) if shape is not None else None,
            "bases_covered": covered,
            "bases_tail": tail,
            "valid": content_ok,
            "note": note,
        }
        status = "ok" if content_ok else "INVALID"
        detail = "" if tail == 0 else f"  covers {covered:,} bases, tail {tail}"
        # No MiB/s column here, deliberately: a view moves no bytes, so
        # dividing the chromosome's size by the time to obtain a pointer
        # produces a number in the tens of TiB/s.  It is arithmetically
        # correct and completely meaningless, and a meaningless number in a
        # results table is indistinguishable from a fake one.
        print(f"   {label:20} {median*1e6:10.2f} us  {nbytes:>13,} bytes  "
              f"shape={shape}  {status}{detail}")

    # ---- correctness, said out loud ---------------------------------------
    for label, row in arrays.items():
        if not row["valid"]:
            print(f"   !! {label}: the array is not chr1 -- its time is not a result")

    # ---- claim 3: memory ---------------------------------------------------
    print("\n== what holding it costs the process")
    memory = []
    for label, getter, _note in getters:
        row = memory_row(label, getter)
        memory.append(row)
        print(f"   {label:20} RSS {row['rss_delta_kb']:>8,} KiB for "
              f"{row['array_nbytes']:>13,} bytes of array")
    print("   (a view of a mapped file adds no anonymous pages; the reference rows")
    print("    must allocate and copy the whole chromosome to have it at all)")

    # ---- claim 2: the operation -------------------------------------------
    print("\n== GC over chr1")
    if "biofasting-index" in handles:
        # The same array the table above measured, not a freshly built one.
        array = handles["biofasting-index"].grid(CHR)[0]
    else:
        import numpy as np

        array = np.frombuffer(reference, dtype=np.uint8)
    ops = op_rows(reference, array, backends, args.repeat)
    for row in ops:
        status = "ok" if row["valid"] else "INVALID"
        print(f"   {row['impl']:22} {row['median_s']*1000:9.3f} ms  "
              f"{row['mib_per_s']:8.1f} MiB/s  GC={row['gc_count']:,}  {status}")
    if ops and ops[0]["ambiguous_bases"]:
        print(f"   note: {ops[0]['ambiguous_bases']:,} ambiguous bases; the fraction's "
              f"denominator differs between libraries, the count does not")

    # ---- the FASTQ finding -------------------------------------------------
    fastq_path = args.data / args.fastq
    headers = None
    if fastq_path.exists():
        print(f"\n== why the reads have no grid  ({fastq_path.name})")
        headers = fastq_header_widths(fastq_path, backends)
        print(f"   {headers['records']:,} records, header widths "
              f"{headers['min']}..{headers['max']} bytes, "
              f"{headers['distinct_widths']} distinct")
        print(f"   {headers['widths']}")
        if "biofasting" in backends:
            try:
                backends["biofasting"].open_fastq_grid(fastq_path)
                print("   the reader built a grid, which contradicts the widths above")
            except ValueError as exc:
                print(f"   refused, as the widths require: {exc}")

    out_path = args.json or (DEFAULT_RESULTS / "grid.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        "corpus": str(fasta_path),
        "chromosome": CHR,
        "bases": len(reference),
        "sequence_sha256": hashlib.sha256(reference).hexdigest(),
        "versions": versions,
        "repeats": args.repeat,
        "index_builds": builds,
        "arrays": list(arrays.values()),
        "memory": memory,
        "ops": ops,
        "fastq_headers": headers,
    }, indent=2) + "\n")
    print(f"\n-> {out_path}")
    invalid = [row for row in list(arrays.values()) + ops if not row["valid"]]
    if invalid:
        print(f"WARNING: {len(invalid)} row(s) are INVALID; their times are not results")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
