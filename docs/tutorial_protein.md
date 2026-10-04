# Tutorial 5 — Protein: sequtils, ProteinAnalysis, IsoelectricPoint

The protein half of the package is the reference's `SeqUtils`/`ProtParam`
machinery, ported with its rules and its refusals, and *measured where the
rule was a surprise*.  This tutorial walks the surface: the module-level
sequence utilities, the twelve-measurement `ProteinAnalysis`, and the
`IsoelectricPoint` class the pI calculation comes from.

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the
owner commits and pushes.

## The sequence utilities

```python
import biofasting

biofasting.gc_fraction("ATTCCGGATTCACTAA")               # 0.375
biofasting.gc_fraction("ATNNGC", ambiguous="ignore")     # 0.3333
```

`gc_fraction`'s modes are the reference's (`remove`/`ignore`/`all`), and the
*default* is the reference's default (`"remove"`) — the corpus generator's
N-count story and the QC example both matter for when `ignore` is the
wrong default, and that is discussed in the QC example's own comments.

```python
biofasting.molecular_weight("MAEGEITTF", seq_type="protein")   # 998.11
biofasting.seq1("MetAlaCysLys")                                # 'MACK'
biofasting.seq3("MAK")                                         # 'MetAlaLys'
biofasting.count_kmers("ATTGG", 2)
```

```
{'AT': 1, 'GG': 1, 'TG': 1, 'TT': 1}
```

`molecular_weight` spells the type (`"DNA"` by default, like the reference
— the refusal message says what it needs instead: a protein passed to the
DNA default raises `'M' is not a valid unambiguous letter for DNA`).
`count_kmers` counts *overlapping* k-mers with no stride — which is exactly
why the run's CAI measurement refused it as a floor for a per-codon walk
(bench/targets.md, eleventh pass).

`GC123` and `GC_skew` are in the same module, same rules.

## Codon adaptation, in two steps

`CodonAdaptationIndex` does two things and the order matters: the weights
table is built once from a set of genes, then `calculate` scores any one
sequence against the table.  The constructor takes an **iterable of
sequences** — not one string (a bare string would be one base at a time,
and the refusal is loud, `TypeError: illegal codon in sequence:`).

```python
genes = ["AAAGGGTTTCCCACG", "GCGCGCATAAAA", "ACGTTTACGTTC"]
index = biofasting.CodonAdaptationIndex(genes)
index.calculate(genes[0])
```

```
1.0
```

The table is the 64 sense codons — stop codons are dropped from the count
(the reference's own rule), which is why an index over random reads
built a 64-entry table rather than a 67-entry one.

The refusal, and why it matters:

```python
index.calculate("ATTGCCA")   # two codons plus a 1-base tail
```

```
TypeError: illegal codon in sequence: A
```

CAI is defined on unambiguous coding sequences; a substrate that carries an
N or a ragged tail is refused by both libraries on the same lookup (`biofasting`
and the reference agree on the raise, and the test
`test_codon_adaptation_index_calculate_rejects_an_illegal_codon` pins the
reference's wording).

The whole pipeline — build from a real gene set, score every member, rank,
report the mean — is `examples/gene_cai.py`; its report on the quick corpus
(the reads as 150-base "genes"):

```
index built from 10,000 genes (500,000 codons, 64 codons weighted)
mean CAI 0.8757, range 0.8024 .. 0.9378
```

On random sequence the mean sits mid-range and flat, as it should — CAI
rewards specialisation and random sequence has none.  The same measurement
pass's numbers (the port's cost against the reference's, and the floors it
derived) are bench/targets.md's eleventh pass — and the pass's own verdict
is worth repeating for a caller deciding where to spend time: CAI is not a
hot path and the Python cleanup is the cheaper lever.

## ProteinAnalysis

One constructor, twelve public measurements:

```python
protein = biofasting.ProteinAnalysis(analyzed, monoisotopic=False)
```

* `molecular_weight()`
* `aromaticity()`
* `gravy(scale="KyteDoolitle")` — 28 hydropathy scales ship
  (`biofasting.protparam_data.gravy_scales`); the default is the reference's
  spelling, typo included
* `isoelectric_point()`
* `charge_at_pH(pH)`
* `count_amino_acids()` / `amino_acids_percent()` — the counter pair
* `secondary_structure_fraction()` — helix / turn / sheet fractions
* `molar_extinction_coefficient()` — `(reduced, oxidised)` coefficients
* `instability_index()`
* `flexibility()` — the profile
* `protein_scale(param_dict, window, edge)` / `_weight_list(window, edge)`
  — the windowed-scale machinery the profile methods use

The whole table, and its instability verdict (`"stable"` under 40.0, the
reference's own threshold), is printed by `examples/protein_profile.py`:

```python
profile(analyzed, name, biofasting.ProteinAnalysis(analyzed))
```

```
-- the reference's doc protein: 152 aa
   molecular weight     17,103.2 Da
   aromaticity            0.0987
   GRAVY                 -0.5974
   isoelectric point        7.72
   charge at pH 7.0         0.93
   helix/turn/sheet   0.33 / 0.29 / 0.37
   extinction coeff    17,420 / 17,545 (reduced / oxidised)
   instability index       41.98  (unstable)
```

The demo protein is the fragment the reference's own class documentation
uses; the example takes `--path` for a real protein and slices its
*windows*, each labelled as a window — because an instability index over
one exon is a statement about that exon (the example's second protein is
"stable" whole and "unstable" as its first window, which is the lesson).

## IsoelectricPoint

`IsoelectricPoint` is the pI computation on its own — bisection on the net
charge curve.  Its constructor takes the protein and, optionally, an
already-counted residue dict (`aa_content` — `ProteinAnalysis`'s
`count_amino_acids()` returns exactly that dict, which is the documented
short cut):

```python
pI = biofasting.IsoelectricPoint("INGAR")
```

```
pI = 9.75, charge_at_pH(7.0) = 0.76
```

(those two values are the class's own doctest; `pi()` and
`charge_at_pH(pH)` are its methods).

## The rules that were measured, not assumed

* the **`sum()` pre‑3.12 quirk** — the reference's instability index sums
  ordered floats under Python's `sum`, whose pre‑3.12 behaviour differs;
  this package pins it (tests/test_protparam.py's pre‑3.12 `sum` test) so
  the number is the reference's on both sides of the boundary;
* **molecular weight** — the protein path uses the reference's residue
  masses including the terminal corrections; the values are asserted
  against the reference on real proteins by tests/test_protparam.py;
* **the GRAVY scales** are data tables carried over as published;
* the **quantifier** (`count_amino_acids`) is a per-letter counter and NOT
  memoised per instance — the ranking pass that measured ProtParam recorded
  that trap for bench writers (a batch re-count benchmark that reuses the
  object measures an empty dict, not the counting).

Where speed is measured on the protein surface (the `gc`/`revcomp`/`ops`
rows, the ProtParam attribution pass), the numbers and their gates are
bench/targets.md's — this tutorial quotes none of them as tables.