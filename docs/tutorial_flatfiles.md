# Tutorial 3 — Flat files: GenBank, EMBL, SwissProt

The flat-file reader is an *index* in the shape of `Bio.SeqIO.index`, with
the FEATURES table and the annotations read on demand — in C++, against the
reference's parser, the same way the FASTQ/FASTA side went.  The three
formats the reference names in its own code are the three this package
names: `genbank`, `embl`, `swiss` — as a `format` parameter, not a guess
(the docstring records why: a LOCUS/ID detector would be right *most of the
time*).

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the
owner commits and pushes.

## Opening a file

```python
index = biofasting.open_genbank("/tmp/demo.gb")
```

The result is a read-only mapping of record id → sequence bytes, exactly what
`Bio.SeqIO.index(path, "genbank")` would give — minus per-fetch parsing:

```python
index.keys()               # ids, in file order
index[record_id]           # the sequence, as bytes
index.name(record_id)      # the LOCUS name
index.description(record_id)
index.sequence_length(record_id)
index.sequence_slice(record_id, 0, 12)   # the window asked for, gathered
```

`index.records()` returns the file in file order as `FlatFileRecord`
namedtuples `(id, name, description, sequence)` — the stream parse for a
caller who wants everything, without paying an id lookup per record:

```python
record = index.records()[0]
```

```
F record: id demoGB | name demoGB | description a small demo record | seq len 23
```

A duplicate id is a file that names two records the same way, and the index
refuses it with `ValueError` — the same refusal `SeqIO.index` makes.  A
gzipped file is also refused, by name, because a silently empty index would
be worse than an error.  (`open_genbank` is named for the format that
carries most of the data; `open_genbank(path, format="embl")` reads EMBL,
`format="swiss"` SwissProt.)

## The FEATURES table

Features are read on demand, per record:

```python
for feature in biofasting.read_features(index, record.id):
    ...
```

a `Feature` namedtuple — the reference's five fields plus the warnings tuple
the table parser adds:

* `feature.type` — `"gene"`, `"CDS"`, …
* `feature.location` — a `Location` (see below), or `None` where the
  reference caught its own `LocationParserError` and carried on (that is a
  behaviour to reproduce, and `feature.status`/`feature.message` say what
  happened);
* `feature.qualifiers` — a tuple of `Qualifier` pairs;
* `feature.qualifiers_as_dict()` — the qualifier dict as the reference
  stores it: key → list of values, repetitions included.

On the demo record — written with `SeqIO.write` from the reference — the
printout of one pass is:

```
F  type = gene | status = ok | warnings = 0
F  location = 0 .. 20 strand 1
F  qualifiers = {'gene_label': ['demoGene']}
F  type = CDS | status = ok | warnings = 0
F  location = 0 .. 20 strand 1
F  qualifiers = {'product': ['demo product'], 'codon_start': ['1']}
```

A `Location` is `Location(operator, parts)`, each `Part` a pair of
`Position`s plus strand and ref; a `Position` carries a kind from the
reference's own set (`exact`, `before`, `after`, `uncertain`, and the
compound-location operators `complement`/`join`/`order` live on `Location`).
`tests/test_location.py` pins the grammar's edge cases — `<`, `>`,
`between(3,9)`, the lot; the reader's job was that the same 5,000 features
parse to the same values the reference's own parser produces.

## The annotations dict

The header's key/value annotations are read on demand the same way:

```python
biofasting.read_annotations(index, record.id)
```

```
F annotations = {'molecule_type': 'DNA', 'data_file_division': 'UNK',
'date': '01-JAN-1980', 'accessions': ['demoGB'], 'keywords': [''],
'source': '', 'organism': 'Demo coli', 'taxonomy': []}
```

That is the whole header the demo file's `LOCUS`/`DEFINITION`/`ACCESSION`/
`VERSION`/`KEYWORDS`/`SOURCE`/`ORGANISM` lines express, mapped to the keys
the reference's `SeqRecord.annotations` uses.  A real GenBank file's
annotations are pinned against the reference's dict in the same way — the
reference's own files have annotations the demo above doesn't (DBLINK and
friends), and the tests carry them.

## The one bridge to Biopython

`biofasting.to_seqfeature(feature)` turns a parsed `Feature` into a real
`Bio.SeqFeature.SeqFeature`, location included — the shape the reference's
feature-driven callers expect:

```
F to_seqfeature: SeqFeature [0:20](+) strand 1
```

that one-bridge rule also decides what the reader does *not* build: no
`SeqRecord` is ever constructed here (the interop shims in the parsing
tutorial are the other direction).

## What is read, and what is deferred

* **Read**: GenBank, EMBL, SwissProt records — the four fields, the FEATURES
  table (with the reference's carry-on behaviours reproduced as
  status/message/warnings), and the annotations dict.
* **Deferred**: flat-file *writers* (measured and filed as a later task),
  and the remaining format names the reference knows (`int`/`clustal` and
  friends are not flat files in this sense at all — they are not reads here).
* **Never**: gzip — by explicit refusal, not silent failure.

The measured numbers for the flat-file reader are in `bench/targets.md` under
the phase‑1 ranking pass it names (`genbank-1kb`, `genbank-1kb-bare`,
`genbank-10kb` rows and their gate — same rule as everywhere else: our
records must equal `SeqIO.parse`'s before any time is quoted).