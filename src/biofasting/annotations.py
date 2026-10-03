"""The annotations dict, and the conversion to Biopython's own values.

``_core.FlatFileIndex.annotations(id)`` reads a record's header and hands back
the keys the header created, in the order it created them, plus their values.
This module is where those become the dict ``SeqRecord.annotations`` holds.

**Who decides what.**  The kernel counts; Python decides.  That is a real
division here and not a slogan: the kernel walks the header's lines, so it is the
only layer that knows *which* lines were there and *in what order* -- and it is
the order that decides the dict.  The reference inserts a key when the line
stating it is consumed, so a record whose ``COMMENT`` line sits above its first
``REFERENCE`` has ``comment`` before ``references``, and ``repr`` shows it.  A
reader that emitted a fixed key list would be right about a corpus and wrong
about a file.  What Python decides is the names, and those are the reference's
own -- ``molecule_type``, not ``mol_type`` -- and the conversion of a reference
block into a ``Bio.SeqFeature.Reference``.

**The two formats differ in which keys a record has**, and the difference is the
reference's, not any one file's:

===========  ==================================  ==============================
key          GenBank                             EMBL
===========  ==================================  ==============================
``date``     the LOCUS date, absent when the     **never** -- EMBL's scanner
             column is blank                     reads no ``DT`` line at all
``source``   ``''`` for a writer's placeholder   **never** -- no ``SOURCEL``
             ``SOURCE`` line                     entry in the consumer table
``keywords`` ``['']`` for ``KW   .``             **no key** unless a ``KW`` line
                                                 was there
===========  ==================================  ==============================

An EMBL record whose ``DT`` line states a date still has no ``date`` in its
annotations, because the consumer that would have written it does not exist.  A
reader that normalised the two formats into one key set would be right about one
of them.

**Where it declines.**  ``read_annotations`` raises :class:`ValueError` -- naming
the line -- for a shape this reader does not reproduce, exactly as
``read_features`` does for a feature table.  The shapes that make it decline are
listed in ``src/core/annotations.hpp``; the short version is that a header key
which sets an annotation key the reader does not produce (NID, PID, DBSOURCE,
SEGMENT, a GI number, a structured comment) is refused, while a key the
reference itself ignores is ignored here too.

Biopython is **not** a dependency of this package: ``Bio`` is imported inside
:func:`to_reference` and nowhere else.
"""

from __future__ import annotations

from typing import NamedTuple

__all__ = [
    "AnnotationReference",
    "read_annotations",
    "to_reference",
]


class AnnotationReference(NamedTuple):
    """One reference: the eight fields ``Bio.SeqFeature.Reference`` holds.

    The reference's *number* is not here because the reference discards it --
    ``REFERENCE 2`` and ``REFERENCE 1`` produce the same object apart from their
    fields.  ``location`` is a list of ``(start, end)`` pairs in Python
    coordinates, already converted by the kernel: the reference's
    ``reference_bases`` assigns rather than extends, so a second ``RP`` line on
    the same EMBL reference replaces the first range instead of adding to it.
    """

    title: str
    authors: str
    consrtm: str
    journal: str
    pubmed_id: str
    medline_id: str
    comment: str
    location: tuple[tuple[int, int], ...]


# The key order each format's consumer runs in, as the reference inserts them.
# Written out rather than derived from the table above, because the order is a
# property of the reference's source and a caller comparing two dicts is
# comparing their insertion orders too.
def _reference(raw) -> AnnotationReference:
    title, authors, consrtm, journal, pubmed_id, medline_id, comment, location = raw
    return AnnotationReference(
        title,
        authors,
        consrtm,
        journal,
        pubmed_id,
        medline_id,
        comment,
        tuple((start, end) for start, end in location),
    )


def read_annotations(index, id: str) -> dict:
    """The annotations dict of ``index[id]``, keyed as the reference keys it.

    Raises :class:`ValueError` when this reader declined the record's header --
    the message names the line that made it decline -- and returns the dict the
    reference would have built otherwise.

    The kernel hands back the keys the header *created*, in the order it created
    them, and this function builds the dict in exactly that order.  That is the
    whole of the naming: the key names are the reference's own, and which of them
    a record has is a fact about the record rather than about the format.  The
    format decides only the four keys its first line creates and the order it
    creates them in, and the kernel has already done that.

    `references` is a list of :class:`AnnotationReference`, and the key is
    present only when a `REFERENCE`/`RN` line created it -- a record with no
    reference block has no such key, not an empty list.
    """
    ok, message, keys, table = index.annotations(id)
    if not ok:
        raise ValueError(f"{id}: {message}")

    out: dict = {}
    for key in keys:
        value = table[key]
        out[key] = (
            [_reference(reference) for reference in value]
            if key == "references"
            else value
        )
    return out


def to_reference(reference: AnnotationReference):
    """Build the ``Bio.SeqFeature.Reference`` for one of ours.

    ``location`` is a **list** of ``SimpleLocation``, not a compound location:
    the reference stores a list there, which is why ``len(location) == 2`` on a
    reference read from a range spanning two lines, and why building a
    ``CompoundLocation`` would compare unequal against ``SeqIO.parse``.

    The attributes are assigned rather than passed to the constructor, because
    ``Reference.__init__`` takes no arguments at all -- it sets its eight fields
    to their empty values and nothing else.  A caller holding a reference to a
    record read by ``SeqIO.parse`` and one built here cannot tell them apart:
    ``__eq__`` compares ``__dict__``, which is what the test relies on.
    """
    from Bio.SeqFeature import Reference, SimpleLocation

    out = Reference()
    out.title = reference.title
    out.authors = reference.authors
    out.consrtm = reference.consrtm
    out.journal = reference.journal
    out.pubmed_id = reference.pubmed_id
    out.medline_id = reference.medline_id
    out.comment = reference.comment
    out.location = [SimpleLocation(start, end) for start, end in reference.location]
    return out
