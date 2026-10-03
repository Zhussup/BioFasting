"""The 27 NCBI genetic codes: a table, and one reader that builds them.

`Bio.Data.CodonTable` is a 1,308-line module of which about 900 lines are 27
genetic codes spelled out one codon at a time.  This is the same data as the
strings NCBI actually publishes -- 64 amino acids per code, `*` at the stops --
plus the machinery that turns a string into the tables the rest of a program
wants to ask questions of.

That is the `data` verdict from `inventory/TRIAGE.md` applied to this module:
the 27 registered codes are data and live in `_codon_tables_data.py`; the
classes are behaviour and are written out, in `_codon_table.py`, because no
table can express what `AmbiguousForwardTable.__getitem__` does or which of
`KeyError` and `TranslationError` it raises when.

    >>> from biofasting import codon_table
    >>> codon_table.unambiguous_dna_by_id[1].forward_table["ATG"]
    'M'
    >>> codon_table.unambiguous_dna_by_id[1].stop_codons
    ['TAA', 'TAG', 'TGA']
    >>> codon_table.unambiguous_dna_by_id[27].forward_table["TGA"]   # not a stop here
    'W'

The tables are built once, at import, from the 5 KiB data file.  Nothing here
imports `Bio`, so the codes are available to a program that has no Biopython
installed -- which is the point of the module, since `Bio.Seq.translate()` and
a great deal else need a genetic code before they can do anything at all.
"""

from __future__ import annotations

from ._codon_table import (
    AMBIGUOUS_DNA_LETTERS,
    AMBIGUOUS_DNA_VALUES,
    AMBIGUOUS_RNA_LETTERS,
    AMBIGUOUS_RNA_VALUES,
    CODON_ORDER,
    EXTENDED_PROTEIN_LETTERS,
    EXTENDED_PROTEIN_VALUES,
    PROTEIN_LETTERS,
    UNAMBIGUOUS_DNA_LETTERS,
    UNAMBIGUOUS_RNA_LETTERS,
    AmbiguousCodonTable,
    AmbiguousForwardTable,
    CodonTable,
    NCBICodonTable,
    NCBICodonTableDNA,
    NCBICodonTableRNA,
    TranslationError,
    build,
    list_ambiguous_codons,
    list_possible_proteins,
    make_back_table,
)
from ._codon_tables_data import NCBI_TABLE_VERSION, TABLES

__all__ = [
    "NCBI_TABLE_VERSION",
    "TranslationError",
    "CodonTable",
    "NCBICodonTable",
    "NCBICodonTableDNA",
    "NCBICodonTableRNA",
    "AmbiguousCodonTable",
    "AmbiguousForwardTable",
    "make_back_table",
    "list_ambiguous_codons",
    "list_possible_proteins",
    "register_ncbi_table",
    "standard_dna_table",
    "standard_rna_table",
    "unambiguous_dna_by_id",
    "unambiguous_dna_by_name",
    "unambiguous_rna_by_id",
    "unambiguous_rna_by_name",
    "generic_by_id",
    "generic_by_name",
    "ambiguous_dna_by_id",
    "ambiguous_dna_by_name",
    "ambiguous_rna_by_id",
    "ambiguous_rna_by_name",
    "ambiguous_generic_by_id",
    "ambiguous_generic_by_name",
]

#   The twelve dictionaries the reference fills in, under the same names and
#   with the same contents: one per alphabet in each of the two shapes, keyed
#   by identifier and by name.  `generic` is the T/U-merged view, which is what
#   a caller who does not know which nucleic acid they have should use.
unambiguous_dna_by_name: dict = {}
unambiguous_dna_by_id: dict = {}
unambiguous_rna_by_name: dict = {}
unambiguous_rna_by_id: dict = {}
generic_by_name: dict = {}
generic_by_id: dict = {}
ambiguous_dna_by_name: dict = {}
ambiguous_dna_by_id: dict = {}
ambiguous_rna_by_name: dict = {}
ambiguous_rna_by_id: dict = {}
ambiguous_generic_by_name: dict = {}
ambiguous_generic_by_id: dict = {}

standard_dna_table = None
standard_rna_table = None

_BY_ID = (
    unambiguous_dna_by_id,
    unambiguous_rna_by_id,
    generic_by_id,
    ambiguous_dna_by_id,
    ambiguous_rna_by_id,
    ambiguous_generic_by_id,
)
_BY_NAME = (
    unambiguous_dna_by_name,
    unambiguous_rna_by_name,
    generic_by_name,
    ambiguous_dna_by_name,
    ambiguous_rna_by_name,
    ambiguous_generic_by_name,
)


def _register(name, alt_name, id, table, start_codons, stop_codons, dual=()) -> None:
    """Build one code's six tables and file them under both of their keys."""
    global standard_dna_table, standard_rna_table

    names, *tables = build(id, name, alt_name, table, start_codons, stop_codons, dual)
    for target, built in zip(_BY_ID, tables):
        target[id] = built
    for target, built in zip(_BY_NAME, tables):
        for registered in names:
            target[registered] = built

    if id == 1:
        standard_dna_table = tables[0]
        standard_rna_table = tables[1]


def register_ncbi_table(
    name, alt_name, id, table, start_codons, stop_codons, dual=()
) -> None:
    """Add a genetic code to the twelve dictionaries.

    `table` is either the 64-character string NCBI publishes or the dictionary
    the reference takes, so code written against `Bio.Data.CodonTable` and code
    written against this module both work.  `dual` names the codons that are a
    stop and an amino acid at once; the reference has no such argument because
    its three codes that need one were written out by hand, and no 64-character
    string can carry that fact.

    The alternative name is appended to the name list unless it is `None` --
    which is not an oversight, it is why table 27's own name list contains a
    `None` upstream and why the printed grid filters false names out.
    """
    _register(name, alt_name, id, table, start_codons, stop_codons, dual)


def _register_all() -> None:
    """Build every code in the data file.  Runs once, at import."""
    for id, name, alt_name, aa_string, start_codons, stop_codons, dual in TABLES:
        _register(name, alt_name, id, aa_string, start_codons, stop_codons, dual)


_register_all()
