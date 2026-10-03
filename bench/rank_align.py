#!/usr/bin/env python3
"""Ranking pass 8: pairwise alignment (PLAN 2.1).

The plan files alignment as wrappers behind the dispatch interface -- parasail
or WFA2-lib, reusing the optimized Smith-Waterman implementations rather than
rewriting one.  But `Bio.Align.PairwiseAligner` is *already* a C extension: a
scalar DP at a few ns a cell.  So the question is not "how slow is the
reference" but "where is the difference, and is it worth a wrapper":

* the many short pairs a mapper aligns -- per-alignment cost, where a Python
  call into either library can dominate the DP itself;
* the long pairs, where the reference's scalar DP runs against parasail's
  SIMD (SSE2/NEON, eight 16-bit lanes);
* the traceback, which costs the reference several times its score-only DP.

The floor is parasail itself: a pure-Python DP at O(n^2) is minutes per 10 kb
pair, so it is stated here and not carried as a column.  WFA2-lib is absent
from this host and is not installable from PyPI, so it is a vendoring decision
rather than a measurement and is left out.

Every row is gated before any time is quoted: both sides get the same scoring
scheme on every pair of the row and must return equal scores, or the run stops
with the mismatch printed.  The conventions the gate pinned down, so that they
are written once and not rediscovered by the wrapper:

* py-parasail takes gap penalties as *positive* numbers and subtracts them;
* a match/mismatch scheme is `parasail.matrix_create("ACGT", match, mismatch)`;
* BLOSUM62 is the same table on both sides, value for value;
* semi-global (free ends on both sequences) is the reference's
  `end_gap_score = 0` and parasail's `sg`, and they agree in either argument
  order, which is what "both ends free" means.

Usage:

    python3 bench/rank_align.py
"""

from __future__ import annotations

import random
import statistics
import sys
import time

import parasail
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

REPEATS = 7
DNA = "ACGT"
AA = "ACDEFGHIKLMNPQRSTVWY"

# The scoring schemes, in one place, so the gate and the timer cannot drift
# apart: DNA is match/mismatch 2/-2 with a 10+1k affine gap, protein is
# BLOSUM62 with 11+1k.
DNA_MATCH, DNA_MISMATCH, DNA_OPEN, DNA_EXTEND = 2, -2, 10, 1
AA_OPEN, AA_EXTEND = 11, 1


def timed(fn, repeats: int = REPEATS) -> float:
    """Median seconds over ``repeats`` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def mutated(rng: random.Random, template: str, rate: float, alphabet: str = DNA) -> str:
    """A real pair: substitutions at ``rate``, plus short indels."""
    letters = list(template)
    for _ in range(max(1, int(len(template) * rate))):
        if rng.randrange(3) or not letters:
            letters[rng.randrange(len(letters))] = rng.choice(alphabet)
        elif rng.randrange(2) and len(letters) > 2:
            del letters[rng.randrange(len(letters))]
        else:
            letters.insert(rng.randrange(len(letters)), rng.choice(alphabet))
    return "".join(letters)


def reference(style: str) -> PairwiseAligner:
    """The reference in each of the three DNA configurations, and protein."""
    aligner = PairwiseAligner()
    if style == "protein":
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -float(AA_OPEN)
        aligner.extend_gap_score = -float(AA_EXTEND)
        return aligner
    aligner.mode = "local" if style == "local" else "global"
    aligner.match_score = float(DNA_MATCH)
    aligner.mismatch_score = float(DNA_MISMATCH)
    aligner.open_gap_score = -float(DNA_OPEN)
    aligner.extend_gap_score = -float(DNA_EXTEND)
    if style == "semiglobal":
        # Both sequences' ends free, which is what parasail's `sg` is.
        aligner.end_gap_score = 0.0
        aligner.extend_end_gap_score = 0.0
    return aligner


def width_for(cells_bound: int) -> int:
    """The narrowest exact parasail width the score bound fits in.

    py-parasail's *unsuffixed* functions are the slow generic dispatch (2.2 ns
    a cell) while the explicitly sized ones are 10 to 28 times faster -- the
    8-bit scan is the full 16-lane width of an SSE register, 16-bit is eight
    lanes.  A score that overflows the chosen width is silently wrong, so the
    width comes from the bound and never from a measurement: ``n * max|score|``
    is the largest score the scheme can produce.
    """
    for bits in (8, 16, 32):
        if cells_bound <= 2 ** (bits - 1) - 1:
            return bits
    return 64


def sized(prefix: str, kind: str, bits: int):
    """`nw_scan_16`, `sg_trace_scan_8`, ... -- always explicitly sized.

    The unsuffixed bindings are the slow generic dispatch; naming the width is
    both what makes the call fast and what makes the overflow bound explicit.
    """
    return getattr(parasail, f"{prefix}_{kind}_{bits}")


def workloads():
    """The rows, built from one seed so a re-run compares like with like."""
    rng = random.Random(7)
    read = "".join(rng.choice(DNA) for _ in range(150))
    pairs_150 = [(read, mutated(rng, read, 0.02)) for _ in range(512)]
    consensus = "".join(rng.choice(DNA) for _ in range(1_000))
    pair_1k = (consensus, mutated(rng, consensus, 0.01))
    big = "".join(rng.choice(DNA) for _ in range(10_000))
    pair_10k = (big, mutated(rng, big, 0.01))
    subject = "".join(rng.choice(DNA) for _ in range(20_000))
    local_pair = (mutated(rng, subject[:1_000], 0.02), subject)
    window = pair_10k[1][3_000:3_150]
    semiglobal_pair = (window, pair_10k[1])
    protein = "".join(rng.choice(AA) for _ in range(300))
    protein_pair = (protein, mutated(rng, protein, 0.05, AA))

    nuc = parasail.matrix_create(DNA, DNA_MATCH, DNA_MISMATCH)
    del nuc
    return {
        # row: (reference style, parasail prefix, pairs)
        "dna-global-150x512": ("global", "nw", pairs_150),
        "dna-global-1k": ("global", "nw", [pair_1k]),
        "dna-global-10k": ("global", "nw", [pair_10k]),
        "dna-semiglobal-150-vs-10k": ("semiglobal", "sg", [semiglobal_pair]),
        "dna-local-1k-vs-20k": ("local", "sw", [local_pair]),
        "protein-blosum62-300": ("protein", "nw", [protein_pair]),
    }


def scheme(style: str, aligner: PairwiseAligner) -> tuple[int, int, int]:
    """`(largest absolute value in the matrix, gap open, gap extend)`.

    Read off the aligner's own matrix rather than remembered, so the width the
    bench picks and the scores the gate checks cannot disagree about the table.
    """
    if style == "protein":
        largest = int(max(abs(float(value)) for value in aligner.substitution_matrix.flat))
        return largest, AA_OPEN, AA_EXTEND
    return max(DNA_MATCH, -DNA_MISMATCH), DNA_OPEN, DNA_EXTEND


def score_bound(pairs, largest: int, open_gap: int, extend_gap: int) -> int:
    """The largest score the scheme can reach, over the row's shortest pair.

    An alignment never has more columns than the shorter sequence, and every
    column scores at most the matrix's largest value; the widest gap it can pay
    for spans the shorter sequence too.  Deliberately loose -- a bound that is
    too tight buys one width and loses correctness.
    """
    shortest = max(min(len(a), len(b)) for a, b in pairs)
    return shortest * (largest + extend_gap) + open_gap


def parasail_call(fn, first: str, second: str, protein: bool, matrix):
    """One parasail call with the row's scheme, gaps positive as it wants."""
    if protein:
        return fn(first, second, AA_OPEN, AA_EXTEND, parasail.blosum62)
    return fn(first, second, DNA_OPEN, DNA_EXTEND, matrix)


def main() -> int:  # noqa: C901 - a bench is a list of rows, not a function
    matrix = parasail.matrix_create(DNA, DNA_MATCH, DNA_MISMATCH)
    aligners = {style: reference(style) for style in ("global", "semiglobal", "local", "protein")}

    # One width per row, from the score bound, used by the gate and the timers
    # alike -- so a row can never be timed on one function and gated on another.
    chosen = {}
    for row, (style, prefix, pairs) in workloads().items():
        largest, open_gap, extend_gap = scheme(style, aligners[style])
        bits = width_for(score_bound(pairs, largest, open_gap, extend_gap))
        chosen[row] = (sized(prefix, "scan", bits), sized(prefix, "trace_scan", bits), bits)

    print("gate: same scheme, both sides, the row's own function, every pair")
    for row, (style, _prefix, pairs) in workloads().items():
        protein = style == "protein"
        scorer = chosen[row][0]
        for first, second in pairs:
            expected = aligners[style].score(first, second)
            got = parasail_call(scorer, first, second, protein, matrix).score
            if expected != got:
                print(f"  {row}: MISMATCH {expected} != {got} (int{chosen[row][2]})")
                print(f"    {first!r}\n    {second!r}")
                return 1
        print(f"  {row}: {len(pairs)} pair(s) equal at int{chosen[row][2]}")

    print()
    print("| row | pairs | width | reference us/align | parasail us/align | speedup |")
    print("|---|---:|---:|---:|---:|---:|")
    results = {}
    for row, (style, _prefix, pairs) in workloads().items():
        protein = style == "protein"
        scorer = chosen[row][0]

        def reference_score():
            for first, second in pairs:
                aligners[style].score(first, second)

        def parasail_score():
            for first, second in pairs:
                parasail_call(scorer, first, second, protein, matrix).score

        ref_s = timed(reference_score) / len(pairs)
        par_s = timed(parasail_score) / len(pairs)
        results[row] = (ref_s, par_s, pairs)
        print(
            f"| {row} | {len(pairs)} | int{chosen[row][2]} | {ref_s * 1e6:9.2f} | "
            f"{par_s * 1e6:9.2f} | {ref_s / par_s:6.2f}x |"
        )

    print()
    print("| row | reference align us/align | parasail trace us/align | speedup |")
    print("|---|---:|---:|---:|")
    for row, (style, _prefix, pairs) in workloads().items():
        protein = style == "protein"
        tracer = chosen[row][1]

        def reference_align():
            for first, second in pairs:
                next(iter(aligners[style].align(first, second)))

        def parasail_trace():
            for first, second in pairs:
                parasail_call(tracer, first, second, protein, matrix)

        ref_a = timed(reference_align, 5) / len(pairs)
        par_a = timed(parasail_trace, 5) / len(pairs)
        print(
            f"| {row} | {ref_a * 1e6:9.2f} | {par_a * 1e6:9.2f} | "
            f"{ref_a / par_a:6.2f}x |"
        )

    print()
    print("score-only, per cell, on the rows where the DP can be seen:")
    for row, (ref_s, par_s, pairs) in results.items():
        cells = sum(len(a) * len(b) for a, b in pairs) / len(pairs)
        print(
            f"  {row}: {cells:.0f} cells/align, reference {ref_s / cells * 1e9:.3f} ns, "
            f"parasail {par_s / cells * 1e9:.3f} ns"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())