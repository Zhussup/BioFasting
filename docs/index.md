# BioFasting — the documentation map

BioFasting is a Python package with a C++ core: the parsers, seq operations,
writers, tables and alignment are compiled and tested against Biopython's own
behaviour, and the Python layer arranges them.  This directory is the *how to
use it* documentation.  The *what it is and how fast it is* documentation lives
in `README.md` (the project) and `bench/targets.md` (every measured number,
gate-first, with the passes that produced them).

**NO COMMITS AND NO PUSHES.** This file is written by the run agent; the owner
commits and pushes.

## The tutorials, in reading order

| file | what it teaches |
|---|---|
| [tutorial_quickstart.md](tutorial_quickstart.md) | install, the build identity, the first FASTQ and FASTA reads, the first writes — ten minutes, no benchmark talk |
| [tutorial_parsing.md](tutorial_parsing.md) | the parsing surface in depth: iteration vs indexed access vs zero-copy grids, gzip by magic bytes, malformed-file refusals, the writers, the `SeqRecord` shims, the Arrow tables |
| [tutorial_flatfiles.md](tutorial_flatfiles.md) | GenBank/EMBL/SwissProt: the flat-file index, the FEATURES table, the annotations dict, and what is read versus deferred |
| [tutorial_restriction.md](tutorial_restriction.md) | the restriction layer as a user sees it: enzymes, the search contract, `catalyze`, batches, the analysis reports |
| [tutorial_protein.md](tutorial_protein.md) | the measuring half: `sequtils` (`GC123`, `molecular_weight`, `seq1`/`seq3`, CAI), the twelve `ProteinAnalysis` measurements, `IsoelectricPoint` |
| [tutorial_alignment.md](tutorial_alignment.md) | `Aligner`/`Alignment`: scoring schemes, modes, end-gap rules, `counts`, and the parasail extra |

## The examples, which the tutorials quote

`examples/` holds six runnable scripts, and they are the tutorials' other
half: every block in a tutorial that shows a *pipeline* is copied from one of
them, and `tests/test_examples.py` runs each one and asserts what it prints.
A tutorial that quotes the QC summary line is quoting a line a test knows.

| example | the pipeline |
|---|---|
| `examples/fastq_qc.py` | read a FASTQ, measure each read, filter with real thresholds, write the passes, report |
| `examples/plasmid_digest.py` | a restriction batch against a circular molecule, positions and fragments |
| `examples/protein_profile.py` | every `ProteinAnalysis` measurement over one protein, plus windows |
| `examples/gene_cai.py` | codon adaptation: build the index from a gene set, score every gene, rank |
| `examples/primer_score.py` | candidate primers ranked by local alignment, with the binding shown |
| `examples/reads_to_polars.py` | FASTQ as the Arrow table, handed to polars, aggregated in polars |

Run them from the repository root; the corpus they default to is the quick one
`python3 bench/gen_data.py --quick` builds, and each explains what it wants
when it is not there.

## How the blocks are kept true

Two kinds of block appear in these tutorials, and the difference is stated
here rather than left implicit.  Pipeline blocks are *copied from the
examples*, which `tests/test_examples.py` runs — the tutorial cannot drift
from the suite silently.  The other kind is a short single-call snippet with
its printed output — the indexed fetches, the refusal messages, the
`SeqFeature` bridge, the `Counts` line — written while this documentation was
and verified by running it on the deterministic corpus (or on a demo file the
*reference's own writer* produced); those values are reproducible facts of
the corpus, and the behaviours they show (refusal texts, result shapes,
bridge types) are pinned by the named test files — `tests/test_fastq.py`,
`tests/test_fastq_index.py`, `tests/test_genbank.py`,
`tests/test_features.py`, `tests/test_annotations.py`,
`tests/test_interop.py`, `tests/test_alignment.py`.  A tutorial that prints
an output the tests cannot recognise is a bug with a name waiting in
`bench/targets.md`'s vocabulary.

## What is refused as documentation, and why

* **No generated API reference.** The docstrings are the API reference, and
  they are exercised by tests; a generated site would be a second source of
  truth nothing keeps in step. `help(biofasting)` is the reference.
* **No notebooks.** They are not diffable and not greppable; a markdown
  tutorial with a copied block does everything a notebook does here, in a form
  the test suite can read.
* **No benchmark tables.** Every measured number lives in `bench/targets.md`
  once, with the pass that measured it and the gate above it. Tutorials cite
  that file by name and state facts, not tables.
* **Only English.** The repository's language rule allows English or Chinese
  for project files; the tutorials exist in English so the tutorials, the
  examples and the tests share one text.

## Install

```console
$ pip install .                 # the package itself, no extras
$ pip install .[arrow]          # + pyarrow, for the tables
$ pip install .[alignment]      # + parasail, for the aligner
```

Extras are optional on purpose; the parsers, seq operations, restriction
layer, protein measurements and flat-file readers are in the bare package.
`biofasting.build_info()` reports what was compiled in — the quickstart's
first block prints it.

## Where the numbers live

Facts this documentation states about speed are measured in
`bench/targets.md`, each under a named pass with a gate that runs before any
time is quoted. The mapping a tutorial relies on most often:

* parsing speed and the indexed-access story — "Phase 1 ranking" passes,
  including the `fai`-style FASTA access measurement;
* the writers — "the writers, gated byte-identical" pass table;
* the Arrow tables and the polars hand-off — "the Arrow tables" pass, which
  also records what `pl.from_arrow` does and does not cost;
* the alignment wrapper — "the alignment wrapper" pass table;
* restriction and CAI — the tenth and eleventh passes, which are also where
  the kernel-candidate finding lives.