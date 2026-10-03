# API triage: what the remaster takes, and what it leaves

Step 0.1 answered *what Biopython exposes*: 3,595 public objects in
`biopython_inventory.csv`.  It did not answer *what is worth rewriting in C++*.
This is that answer, as a prefilled `verdict` on every row.

```sh
python3 inventory/triage.py     # -> inventory/biopython_triage.csv
```

The verdicts live in their own file rather than in a new column of the
inventory, because `scan_biopython.py` regenerates the inventory outright and a
rescan would silently discard a judgement.  The two join on
`(module, qualname, kind)`.

## Verdicts

| Verdict | Rows | Share | Meaning |
|---|---:|---:|---|
| `perf` | 610 | 17.0% | hot path with a real speed gap — a kernel to write |
| `parity` | 1,579 | 43.9% | worth having for drop-in compatibility, no large win expected |
| `data` | 1,254 | 34.9% | becomes data + a parser; never ported as code |
| `skip` | 152 | 4.2% | **not** ported — mechanical reasons only, listed in full below |

**3,443 of 3,595 rows (95.8%) are on the remaster.**  `data` extends the
perf/parity/skip vocabulary recorded in PLAN.md, because the Restriction enzyme
classes and the codon tables are neither slow work nor API-parity work: the
plan's own answer for them is "data + parser".

The counts moved once after the first pass, by one row: `Bio.SeqUtils.xGC_skew`
was filed `perf` and is a `skip`.  It imports `tkinter`, creates a `Canvas`,
plots two skew curves and returns `None` — 237 ms a call with `DISPLAY` set, no
value to compare and no algorithm to port.  The SeqUtils work (PLAN 2.2h)
re-derived that by hand, which is exactly the cost a wrong `perf` verdict
imposes on the next pass, so the rule was added rather than the row patched.

## The skip policy

`skip` is deliberately hard to earn.  It requires a **mechanical** reason:

1. deprecated upstream — the inventory's own column, or deprecated in fact;
2. the object's job is to make a network request;
3. the object's job is to launch another program;
4. upstream already implements it in C, so a rewrite buys nothing;
5. the numerics belong to a third-party library (scipy);
6. the function's job is to draw a window rather than to compute a value, so
   there is no result to compare and no algorithm to port.

Everything arguable is written down as work instead.  A judgement call resolved
towards "do the work" costs one row in a list; resolved towards "skip" it costs
a feature nobody notices is missing until someone needs it.

That rule was broken once, by this script, and the break is worth recording
rather than quietly fixing.  The first draft skipped whole network *packages* by
module name — which took `Bio.SwissProt.parse`, `Bio.Medline.Record`,
`Bio.KEGG.Enzyme.parse`, `Bio.ExPASy.Prosite.parse` and every other parser filed
under those packages with it.  Those parse files on local disk; the HTTP call
they are filed next to is a different object.  They are `parity` now, and the
reason column says why the module name is misleading: 52 rows left the network
category (94 → 42), and skip went from 201 rows to 151.

## What is skipped, in full

Grouped by reason, which is the column to argue with.

**Deprecated upstream — 13 rows.** `Bio.GenBank.{Iterator,FeatureParser,RecordParser}`
(flag), `Bio.Align.AlignInfo.SummaryInfo`, `Bio.ExPASy.get_prosite_raw`, and all
8 rows of `Bio.pairwise2` — the last eight because `pairwise2` is deprecated in
fact, not merely in the column: upstream replaced it with
`Bio.Align.PairwiseAligner`.

**Network transport — 42 rows.** `Bio.Entrez` (8 of 10: `efetch`, `esearch`,
`epost`, `elink`, `einfo`, `esummary`, `espell`, `ecitmatch`; `read` and `parse`
stay, they parse a document already in hand), `Bio.PDB.PDBList` (14),
`Bio.KEGG.REST` (6), `Bio.TogoWS` (5), `Bio.ExPASy` (3), `Bio.PDB.alphafold_db`
(3), `Bio.Blast.NCBIWWW.qblast`, `Bio.UniProt.search`,
`Bio.ExPASy.ScanProsite.scan`.  A network round trip is 10⁴–10⁶ times the cost
of the parse that follows it; there is nothing for C to win.

**Launches an external program — 46 rows.** `Bio.Phylo.PAML.{_paml,codeml,baseml,yn00}`
(28: the PAML binaries), `Bio.Emboss.{Primer3,PrimerSearch}` (9),
`Bio.PDB.{DSSP,NACCESS,ResidueDepth,PSEA}` (9 — only the entry points that shell
out to `dssp`/`naccess`/`msms`/`psea`, plus the classes that do it from their
constructors; their parsers and pure-Python geometry stay).

**Already C upstream — 43 rows.** `Bio.Cluster` (21) and `Bio.Cluster._cluster`
(9) — `kcluster`, `kmedoids`, `treecluster`, `pca`; `Bio.PDB.qcprot` (10),
`Bio.PDB.ccealign.run_cealign`, `Bio.Nexus.cnexus.scanfile`,
`Bio.cpairwise2.rint`.  These are reuse candidates, not rewrite candidates.

**Third-party numerics — 7 rows.** `Bio.phenotype.pm_fitting`: logistic,
gompertz and richards fits on `scipy.optimize.curve_fit`.

**Draws instead of returning — 1 row.** `Bio.SeqUtils.xGC_skew`.  It is worth
naming here rather than only in the CSV because it was `perf` for a day, and the
mistake was not the count but the direction: a `perf` row is an instruction to
the next pass to measure and then write a kernel, and no measurement could have
produced one.

## The kernel candidates (`perf`)

610 rows, 39 modules — the surface the C core is actually for.  Roughly three
clusters:

- **Sequence.** `Bio.Seq` (86), `Bio.SeqRecord` (12), `Bio.SeqIO.QualityIO` (28,
  step 1.2), `Bio.SeqIO.FastaIO` (13, step 1.2), `Bio.SeqIO._index` (21, step
  1.2), `Bio.SeqUtils` and its five satellites (39), `Bio.bgzf` (21).
- **Structure.** `Bio.PDB.internal_coords` (61), `Bio.PDB.Atom` (41),
  `Bio.PDB.Entity` (33), `Bio.PDB.vectors` (30), `Bio.PDB.Residue` (14) and the
  parsers, plus `Bio.SASA`/`Superimposer`/`NeighborSearch`/`HSExposure`/
  `cealign`.
- **Alignment and phylogeny.** `Bio.Align` (41), `Bio.Align.bigbed` (43),
  `Bio.Phylo.TreeConstruction` (22), `Bio.Phylo.Consensus` (13),
  `Bio.motifs.matrix` (22).

PDB and alignment appear here because they are genuinely hot, not because they
are next: **this list says what is on the remaster, not in what order.**
`bench/targets.md` still owns the ordering, and on measured evidence it puts
FASTQ and FASTA random access first, not the 61 `internal_coords` methods.

## Status

This is a *prefill*, not a decision.  By the standing split the owner validates
the contested rows — which is what the `reason` column is for: review the
reasons, not the 3,595 rows.  The rows most likely to be argued with are the 43
"already C upstream" ones, since reusing a foreign C extension is a strategy
choice, not a mechanism.
