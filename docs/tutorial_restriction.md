# Tutorial 4 — Restriction: enzymes, the search contract, batches

`biofasting.restriction` is a port of `Bio.Restriction` written as *data plus
one engine*: the site/overhang/supplier tables in one module and a single
`Enzyme` class, instead of the reference's thousands of generated class
files with their mixin lattice.  Every value — sites, cut positions,
overhangs, frequencies, supplier codes — matches the reference's, element by
element; the tests for that are as long as the table.  This tutorial is what
a user needs: the five calls a digest is made of.

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the
owner commits and pushes.

## One enzyme

Enzymes are module attributes and `get_enzyme(name)`:

```python
import biofasting.restriction as r

eco = r.EcoRI
```

```python
eco.site            # 'GAATTC'
eco.is_palindromic()
eco.elucidate()
```

The elucidation draws the cuts on the site — `^` above, `_` underneath:

```
G^AATT_C
```

An enzyme a name does not match is a `KeyError` with the name in it
(`unknown restriction enzyme: 'EcoXX'`), the same shape a batch gives.

## The search contract

`enzyme.search(sequence)` returns the cut positions, 1-based, ascending —
the same list for the same molecule, both libraries:

```python
eco.search("TTTGAATTCAAAGAATTCGG")
```

```
[5, 14]
```

The sequence argument is what you have: a `str`, our readers' `bytes`, or a
`Seq` — the reference requires its `Seq` and builds a `FormattedSeq` around
it; this package accepts the plain string.  The positions count from the
molecule's first base and say where the enzyme *cuts* (the reference's own
convention, which the tests pin against the reference call for call).

Two search rules that bite when forgotten:

* **Both strands.** A palindromic site is one search; a non-palindromic site
  is found on the complementary strand too, and the returned positions are
  merged.  `EcoRI` (`GAATTC`) is its own reverse complement; an enzyme like
  `Aba13301I` (`GCAAAC`, whose complement `GTTTGC` a molecule can carry
  alone) is not — verified from the table:

  ```python
  comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
  rc = "".join(comp[c] for c in reversed("GCAAAC"))
  r.get_enzyme("Aba13301I").search("AAAA" + rc + "AAAA")
  ```

  ```
  [5]
  ```

* **Overlap.** Sites may overlap: `GCGCGC` carries two overlapping `HhaI`
  sites (bases 1–4 and 3–6), and the search returns both cuts:

  ```python
  r.get_enzyme("HhaI").search("GCGCGC")
  ```

  ```
  [4, 6]
  ```

  Each library's scanner is a regex pass over the same string, so this is
  also the place the two libraries are literally the same code.

## The digestion: `catalyze`

`enzyme.catalyze(sequence, linear=True)` cuts the molecule and returns the
fragments in order.  `str` in, `str` out — and, the rule that trips the
first-timer on every library in the family:

```python
eco.catalyze("TTTGAATTCAAAGAATTCGG")
```

```
['TTTG', 'AATTCAAAG', 'AATTCGG']
```

two cuts, three fragments.  For a **circular** molecule pass
`linear=False`:

```python
ring = "AAAAAGAATTCCTTTGGAATTCCAAA"
eco.catalyze(ring, linear=False)
```

```
['AATTCCAAAAAAAAG', 'AATTCCTTTGG']
```

and the lengths sum to the ring (15 + 11 = 26) — the positions that wrapped
before the origin have been moved to the back by the search, and the last
fragment is the one that crosses the origin.  One cut on a ring opens it
whole rather than into two pieces; `examples/plasmid_digest.py` prints a
workload of exactly those numbers (a 10,000 bp demo ring with five cutting
enzymes — and the library's answer that a one-cut enzyme leaves the ring
opened whole is a test-asserted line there).

An enzyme the tables carry but whose restriction is unknown (`is_unknown()`)
refuses `catalyze` and prints as `? SITE ?` in `elucidate`.

## A batch

```python
batch = r.RestrictionBatch(["EcoRI", "BamHI"])
batch.elements()                      # every enzyme, sorted by name
```

```
['BamHI', 'EcoRI']
```

`batch.search(dna)` — or the same through the `/` operator,
`batch / dna` — searches them all and returns a dict:

```python
batch / "TTTGAATTCAAAGGATCCTTTGATC"
```

```
{'BamHI': [14], 'EcoRI': [5]}
```

Two batch rules worth knowing:

* the search accepts a `linear=False` like the enzyme's own, and the same
  wrap moves apply to every position in the mapping;
* **the mapping is cached on the batch**, keyed on the formatted sequence —
  asking twice about the same molecule answers once, which means a caller
  timing digests (or a caller who wants a fresh answer after mutating
  nothing but *expecting* a recompute) should know the cache is a feature
  with teeth.  `bench/rank_restriction.py` builds a fresh batch per timed
  call for exactly this reason, and the trap is written up in
  `bench/targets.md`'s tenth pass.

## The report: `Analysis`

`Analysis(batch, sequence, linear)` runs the search at construction and
prints the report:

```python
analysis = r.Analysis(r.RestrictionBatch([r.EcoRI, r.BamHI]), "TTTGAATTCAAAGGATCC")
analysis.print_that()
```

```
BamHI      :  14.
EcoRI      :  5.
```

Its filters (`blunt`, `overhang5`, `between`, `with_name`, …) return search
dicts of the same shape, so they compose; `to_location()` and
`to_seqfeature()` bridge the other layers of this package (see the flat-file
tutorial) and the reference's feature objects.

## What this port is, one paragraph

One `Enzyme` class and a 1,088-row site/supplier table replace the
reference's generated classes; 51 mixin-combination behaviours became one
class with the same public surface.  On speed the port sits *at* the
reference, because both implementations scan with `re.finditer` over a
prepared string — equal measured on the whole call.  Where the time *goes*
is also documented (the attribution table in `bench/targets.md`'s tenth
pass: the shared scan is most of the call, and a `str.find` floor sits
well below both), because that pass is the input for a kernel decision
that is the owner's, not a fact of the library.  A tutorial repeats only
what a user must not assume: the port is not faster than the reference
today, it is the reference's behaviour with a smaller surface and no
generated code.