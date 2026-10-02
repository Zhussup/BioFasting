# Biopython 1.88 API inventory
Generated 2026-10-02 by `scan_biopython.py`.

**Totals:** 284 modules (6 failed to import), 156,755 lines of code, 1670 classes, 731 functions, 3172 methods, 289 properties, 13 constants.
Public objects in CSV: 3,595.

| Package | Modules | Classes | Functions | Methods | Constants | LOC |
|---|---:|---:|---:|---:|---:|---:|
| Restriction | 4 | 1112 | 1 | 154 | 20704 | 27,613 |
| SearchIO | 32 | 54 | 52 | 526 | 67 | 12,445 |
| Blast | 5 | 27 | 11 | 522 | 28 | 5,802 |
| PDB | 45 | 59 | 94 | 379 | 48 | 19,825 |
| Phylo | 25 | 67 | 93 | 272 | 20 | 11,917 |
| SeqIO | 24 | 71 | 89 | 226 | 156 | 13,476 |
| Align | 30 | 67 | 30 | 262 | 192 | 16,170 |
| GenBank | 4 | 16 | 4 | 186 | 44 | 3,882 |
| KEGG | 9 | 10 | 18 | 142 | 0 | 2,202 |
| Nexus | 6 | 15 | 15 | 133 | 0 | 3,430 |
| motifs | 14 | 19 | 39 | 67 | 2 | 3,343 |
| Seq | 1 | 7 | 9 | 93 | 10 | 3,278 |
| AlignIO | 11 | 20 | 11 | 50 | 14 | 4,303 |
| SCOP | 7 | 13 | 9 | 59 | 0 | 1,719 |
| SeqFeature | 1 | 16 | 0 | 61 | 13 | 2,332 |
| Entrez | 2 | 12 | 16 | 41 | 5 | 1,912 |
| Pathway | 4 | 6 | 2 | 46 | 0 | 674 |
| codonalign | 3 | 2 | 33 | 16 | 1 | 2,644 |
| SeqUtils | 7 | 3 | 24 | 22 | 0 | 3,011 |
| ExPASy | 6 | 8 | 24 | 12 | 2 | 1,137 |
| Cluster | 2 | 3 | 24 | 13 | 0 | 1,289 |
| phenotype | 3 | 3 | 13 | 24 | 0 | 1,604 |
| Sequencing | 3 | 13 | 7 | 14 | 0 | 821 |
| bgzf | 1 | 2 | 5 | 25 | 0 | 960 |
| UniProt | 2 | 1 | 20 | 7 | 0 | 665 |
| PopGen | 4 | 3 | 5 | 19 | 0 | 702 |
| SwissProt | 2 | 5 | 18 | 3 | 0 | 939 |
| pairwise2 | 1 | 4 | 13 | 9 | 8 | 1,441 |
| File | 1 | 3 | 2 | 17 | 3 | 626 |
| SeqRecord | 1 | 2 | 0 | 19 | 2 | 1,544 |
| Data | 4 | 7 | 6 | 6 | 10 | 2,124 |
| Emboss | 3 | 5 | 3 | 7 | 0 | 271 |
| Compass | 1 | 1 | 10 | 3 | 0 | 223 |
| UniGene | 1 | 4 | 3 | 7 | 0 | 335 |
| NMR | 3 | 2 | 7 | 4 | 0 | 424 |
| TogoWS | 1 | 0 | 11 | 0 | 0 | 374 |
| SVDSuperimposer | 1 | 1 | 0 | 9 | 0 | 203 |
| CAPS | 1 | 3 | 0 | 4 | 0 | 130 |
| Affy | 2 | 2 | 3 | 1 | 0 | 578 |
| Geo | 2 | 1 | 3 | 1 | 0 | 157 |
| Medline | 1 | 1 | 2 | 0 | 0 | 225 |
| cpairwise2 | 1 | 0 | 2 | 0 | 0 | 0 |
| Alphabet | 1 | 0 | 0 | 0 | 0 | 0 |
| Graphics | 1 | 0 | 0 | 0 | 0 | 0 |
| HMM | 1 | 0 | 0 | 0 | 0 | 5 |

## Deprecation warnings at import

- `Bio.phenotype.pm_fitting`

## Failed imports (optional dependencies)

- `Bio.Alphabet` — ImportError: Bio.Alphabet has been removed from Biopython. In many cases, the alphabet can simply be ignored and removed from scripts. In a few cases, you may need to specify the ``molecule_type`` as an annotation on a SeqRecord for your script to work correctly. Please see https://biopython.org/wiki/Alphabet for more information.
- `Bio.Graphics` — MissingPythonDependencyError: Please install ReportLab if you want to use Bio.Graphics. You can find ReportLab at http://www.reportlab.com/software/opensource/
- `Bio.PDB.binary_cif` — MissingPythonDependencyError: Install msgpack to use Bio.PDB.binaryCIF (e.g. pip install msgpack)
- `Bio.PDB.mmtf` — MissingPythonDependencyError: Install mmtf to use Bio.PDB.mmtf (e.g. pip install mmtf-python)
- `Bio.Phylo.CDAOIO` — MissingPythonDependencyError: Support for CDAO tree format requires RDFlib.
- `Bio.motifs.jaspar.db` — MissingPythonDependencyError: Install MySQLdb if you want to use Bio.motifs.jaspar.db
