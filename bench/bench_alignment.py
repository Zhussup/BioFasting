#!/usr/bin/env python3
"""The delivered pairwise aligner, measured against the reference it wraps.

`bench/rank_align.py` is the *target* half: it ranked `Bio.Align.PairwiseAligner`
against parasail's kernels directly and found 6x-35x score-only and 9x-38x with
the traceback, on six rows, gate first.  Its numbers are what
`src/biofasting/alignment.py` was written to, and its finding that decides the
shape of the wrapper is that the *unsuffixed* `nw_scan`/`sg_scan`/`sw_scan`
bindings are the slow generic dispatch while the explicitly sized ones are the
SIMD kernels -- so this table is the other half: the wrapper's own numbers, on
the same six rows, built from the same seed, under the same gate.

What is measured is the whole delivered call, not the kernel: building the
scheme, the width from the score bound, the alphabet and matrix lookup, the CIGAR
decode, and the walk that produces the gapped strings, the coordinates and the
counts.  That is the honest number -- a caller cannot buy the kernel without
them -- and it is the reason one column below is the *short*-pair row: 512 pairs
of 150 bases is where a per-call cost shows up against a per-cell one.

**The gate is the score, on every pair.**  The wrapper and the reference are run
on the same scheme and the same pair, and a row whose scores disagree is printed
`INVALID` and contributes no time to the table.  A path is compared separately
and only where the reference has exactly one optimal alignment: the reference
enumerates every optimum and this returns one, so "the same alignment" is a claim
about optima that only a unique optimum can make.

Usage:

    python3 bench/bench_alignment.py [repeats]
"""

from __future__ import annotations

import random
import statistics
import sys
import time

REPEATS = 9

DNA = "ACGT"
AA = "ACDEFGHIKLMNPQRSTVWY"

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


def workloads(ours_matrix, theirs_matrix):
    """The six rows of the ranking pass, from the same seed.

    Same shapes and same seed as `rank_align.workloads` on purpose: a difference
    between the two tables is then a difference the wrapper made and not a
    difference the input made.  Each row is `(ours, theirs, pairs)` -- the same
    scheme written twice, once with the wrapper's spelling and once with the
    reference's.

    The two spellings differ in one place, and it is a trap rather than a
    detail: `PairwiseAligner` accepts only its *own* `Array` type, and given
    another array with an `alphabet` it takes it and then scores nonsense --
    measured here as `ValueError: sequence item 0 is out of bound (71, should be
    < 24)`, a BLOSUM62 alphabet that never reached the letter table.
    `biofasting.alignment` takes either, which is what
    `tests/test_alignment.py::test_a_protein_scheme_scores_the_same_with_either_blosum62`
    is for, so the bench hands each side its own and the numbers stay comparable.
    """
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

    dna = dict(
        match_score=DNA_MATCH, mismatch_score=DNA_MISMATCH,
        open_gap_score=-DNA_OPEN, extend_gap_score=-DNA_EXTEND,
    )
    nw = dict(dna, mode="global", end_gap_score=-DNA_OPEN, extend_end_gap_score=-DNA_EXTEND)
    sg = dict(dna, mode="global", end_gap_score=0, extend_end_gap_score=0)
    sw = dict(dna, mode="local")

    def protein_scheme(matrix):
        return dict(
            mode="global", substitution_matrix=matrix,
            open_gap_score=-AA_OPEN, extend_gap_score=-AA_EXTEND,
            end_gap_score=-AA_OPEN, extend_end_gap_score=-AA_EXTEND,
        )

    return {
        "dna-global-150x512": (nw, dict(nw), pairs_150),
        "dna-global-1k": (nw, dict(nw), [pair_1k]),
        "dna-global-10k": (nw, dict(nw), [pair_10k]),
        "dna-semiglobal-150-vs-10k": (sg, dict(sg), [semiglobal_pair]),
        "dna-local-1k-vs-20k": (sw, dict(sw), [local_pair]),
        "protein-blosum62-300": (
            protein_scheme(ours_matrix),
            protein_scheme(theirs_matrix),
            [protein_pair],
        ),
    }


def reference(attributes):
    """`Bio.Align.PairwiseAligner` for one of those schemes."""
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    for name, value in attributes.items():
        if name == "substitution_matrix":
            aligner.substitution_matrix = value
        elif name == "mode":
            aligner.mode = value
        else:
            setattr(aligner, name, float(value))
    return aligner


def main(argv: list[str]) -> int:
    repeats = int(argv[1]) if len(argv) > 1 else REPEATS

    from biofasting import alignment, substitution_matrices

    if not alignment.available():
        print(
            "the optional parasail accelerator is not installed, so there is "
            "nothing to measure here: pip install biofasting[alignment]"
        )
        return 1

    from Bio.Align import substitution_matrices as reference_matrices

    rows = workloads(
        substitution_matrices.load("BLOSUM62"), reference_matrices.load("BLOSUM62")
    )

    print("gate: the wrapper and the reference, same scheme, every pair")
    preparsed = {}
    for row, (ours_attributes, theirs_attributes, pairs) in rows.items():
        ours = alignment.Aligner(**ours_attributes)
        theirs = reference(theirs_attributes)
        for first, second in pairs:
            mine, expected = ours.score(first, second), theirs.score(first, second)
            if mine != expected:
                print(f"  {row}: INVALID {expected} != {mine}")
                print(f"    {first!r}\n    {second!r}")
                return 1
        scheme = ours._scheme(*pairs[0])
        print(f"  {row}: {len(pairs)} pair(s), score equal, "
              f"{scheme.prefix}_scan_{scheme.bits}")
        preparsed[row] = (ours, theirs)

    print()
    print("| row | pairs | reference us/align | wrapper us/align | speedup |")
    print("|---|---:|---:|---:|---:|")
    score_only = {}
    for row, (_ours, _theirs, pairs) in rows.items():
        ours, theirs = preparsed[row]

        def reference_scores():
            for first, second in pairs:
                theirs.score(first, second)

        def wrapper_scores():
            for first, second in pairs:
                ours.score(first, second)

        ref_s = timed(reference_scores, repeats) / len(pairs)
        our_s = timed(wrapper_scores, repeats) / len(pairs)
        score_only[row] = (ref_s, our_s, pairs)
        print(
            f"| {row} | {len(pairs)} | {ref_s * 1e6:9.2f} | {our_s * 1e6:9.2f} | "
            f"{ref_s / our_s:6.2f}x |"
        )

    print()
    print("| row | reference align us/align | wrapper align us/align | speedup |")
    print("|---|---:|---:|---:|")
    for row, (_ours, _theirs, pairs) in rows.items():
        ours, theirs = preparsed[row]

        # `next(iter(...))` on both sides: the reference yields every optimal
        # alignment and this returns one, so the comparison is the first optimum
        # against the one -- otherwise the reference would be charged for work
        # this module does not do.
        def reference_aligns():
            for first, second in pairs:
                next(iter(theirs.align(first, second)))

        def wrapper_aligns():
            for first, second in pairs:
                ours.align(first, second)

        ref_a = timed(reference_aligns, repeats) / len(pairs)
        our_a = timed(wrapper_aligns, repeats) / len(pairs)
        print(
            f"| {row} | {ref_a * 1e6:9.2f} | {our_a * 1e6:9.2f} | "
            f"{ref_a / our_a:6.2f}x |"
        )

    print()
    print("per cell, where the DP is big enough to be seen:")
    for row, (ref_s, our_s, pairs) in score_only.items():
        cells = sum(len(a) * len(b) for a, b in pairs) / len(pairs)
        print(
            f"  {row}: {cells:.0f} cells/align, reference {ref_s / cells * 1e9:.3f} ns, "
            f"wrapper {our_s / cells * 1e9:.3f} ns"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
