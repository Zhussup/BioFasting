# Tutorial 6 — Alignment: Aligner, Alignment, Counts

`biofasting.alignment` is a Python wrapper over alignment engines, with the
scoring schemes, the modes and the result objects.  This tutorial is the
caller's view: how to price an alignment, what the modes mean, what an
`Alignment` answers, and where the measured numbers live.

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the
owner commits and pushes.

## The Aligner

The constructor takes *keyword-only* prices — the shape
`Bio.Align.PairwiseAligner` takes, spelled the same way so a caller's
config ports:

```python
al = biofasting.Aligner(
    mode="global",
    match_score=1.0,
    mismatch_score=0.0,
    open_gap_score=-1.0,
    extend_gap_score=-1.0,
)
```

```python
a = al.align("ACGTACGT", "ACGATCGT")
a.score
```

```
6
```

Positives-only, that is the whole game of that scheme: matches score, walks
score nothing, gaps cost — and the pair printed as the two gapped strings:

```
ACGTACGT
ACGATCGT
```

Positional prices are refused (the constructor takes one positional, the
mode; the rest must be named) — `Aligner(1.0, 0.0)` raises `TypeError` with
the constructor's own message.  The refusals are part of the contract
(`tests/test_alignment.py` pins them): a price you did not mean to set
defaulting by accident is exactly the bug this refuses.

## Modes and end gaps

* `mode="global"` — the reference's global rules, including its
  `end_gap_score`/`extend_end_gap_score` prices (which differ from the
  interior ones and govern where terminal gaps may open).
* `mode="local"` — the best-scoring subregion; end gaps are free by the
  mode's own definition, and the accepted end-price values are the same
  either way.  These two are the modes the wrapper has; anything else is
  refused (`"fogsaa"` with its own message — it is the reference's
  experimental mode, and this library's engines don't carry it).

The reference's own subtlety — that *local* mode accepts any end prices
while *global* applies its end-gap grammar — is reproduced and pinned, not
"fixed": a port that silently normalised what the reference refuses would
be a different library.

## Prices that say something

The primer example makes the trade explicit.  Matches `2`, mismatches
`-3`, gap open `-5`, gap extend `-2` — a clean 20-mer scores 40, the same
20-mer with three substitutions scores less, and *local* mode is the right
mode because a primer is a needle in a longer template:

```python
loc = biofasting.Aligner(mode="local", match_score=2.0,
                         mismatch_score=-3.0,
                         open_gap_score=-5.0, extend_gap_score=-2.0)
b = loc.align("ACGT", "ACGTTTTTTTACGT")
b.score
```

```
8
```

—the probe's whole 4 bases match *twice* (both ends of the template carry
the same 4 letters) and the aligner answers with the best of them.  The
full pipeline — candidate primers cut from a template, one decoy carrying
deliberate substitutions, ranked, with the binding's aligned columns and
template span printed — is `examples/primer_score.py`:

```
rank  score  candidate
   1     40  clean: template[50:70]
   ...
   6     28  decoy: template[950:970] plus 3 substitutions

the decoy (decoy, template[950:970] plus 3 substitutions) binds as:
  primer    AGTCTGGATTTCCCCAGTT
  template  AGTATGGTTTTCCCCAGTT
  template span 951..970, score 28
```

`aligner.align(seq1, seq2)` returns that `Alignment` — and a *decoy* that
still binds is exactly the answer the ranking is for.

## The Alignment

An `Alignment` answers three ways, all on demand:

* `.score` — the number the prices produced;
* `.sequences` — the two input sequences, as the aligner kept them;
* `.coordinates` — the alignment's `(2, N)` array of spans, the reference's
  representation; a span pair `[start, end]` per column, so
  `coordinates[1][0]`…`coordinates[1][-1]` is where the alignment sits in
  the second sequence;
* `alignment[0]` / `alignment[1]` — the gapped strings, aligned columns.

```
coords: [[0, 3, 3, 5, 5, 8], [0, 3, 3, 5, 5, 8]]
```

(the global example above: match, gap, match, gap, match — the coordinate
pairs are explicit about the gap spans).

## Counts

`alignment.counts()` builds the difference summary — identities,
mismatches, gaps, and the two scores' decompositions — *lazily*: the
attribute is a method, nothing is computed until the call:

```python
a.counts()
```

```
Counts(identities=6, mismatches=2, gaps=0, substitution_score=6, gap_score=0,
score=6)
```

`Counts` is a small named tuple, and its fields are the ones a caller
downstream usually wants: `identities/length` for identity percent,
`mismatches` for SNP counting, `gaps` for indel counts.

## Substitution matrices

Matrices ship as data with the reference's names, loaded by `load()` —
`BLOSUM62`, `NUC.4.4`, and the rest:

```python
import biofasting
from biofasting.substitution_matrices import load

blosum = load("BLOSUM62")
al = biofasting.Aligner(mode="global", substitution_matrix=blosum,
                        open_gap_score=-1.0, extend_gap_score=-1.0)
al.align("MAEGEITTF", "MAEGEITTF").score
```

```
45
```

One amino-acid change drops it to 42 — visible in the pair.  A `dict` of
`{(letter, letter): score}` works the same way as the shipped arrays; the
tests drive both.

## Where the speed lives

The wrapper's measured table is `bench/targets.md`'s alignment pass — and
its own findings matter to a caller as much as to a maintainer: the wrapper
*itself* was measurably slower than its engine until the scoring-matrix
rebuild and the per-call conversions were taken out of the call path
(BLOSUM62 cost more to build than the kernel took to align, which is why it
is memoised now), and the wrapper exposes `parasail` as an **optional
extra** (`pip install biofasting[alignment]`) because parasail ships no
aarch64 wheels — an arm64 build from sources is one of the run's open owner
decisions, as is the WFA2-lib experiment for long sequences.  The tutorial
states where those decisions live; it does not resolve them.