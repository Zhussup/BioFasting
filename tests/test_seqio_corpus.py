"""The generated flat-file corpus, checked before any number rests on it.

`bench/seqio_corpus.py` is what the sixth ranking pass (PLAN 2.4) measures, and
a number measured on a corpus nobody checked is a number about an unknown file.
The corpus is generated rather than stored, so `tests/test_corpus.py` -- which
walks `bench/data/` against a manifest -- cannot see it.  This is the same
gate for the half of the corpus that never touches the disk.

Three layers:

- **Shape.** Every row yields the record count, record length and feature count
  its own description promises, and the reference parses all of it.
- **Round trip.** What the generator put in comes back out: the ids it minted,
  the residues its keystream produced, and -- for the nucleotide formats, whose
  writer is the reference's -- the division, the topology and the organism.
- **Reproducibility.** Building the same row twice gives byte-identical text,
  which is the property the whole generator rests on and the one a future edit
  to it is most likely to break without noticing.

The point of the bare row is asserted rather than assumed: `genbank-1kb-bare`
exists to isolate the cost of building `SeqFeature` objects, so a change that
quietly gave it features would turn the ranking pass's decomposition into
nonsense while every other test still passed.  The same is asserted for the
*other* axis: `genbank-1kb-annot` and `embl-1kb-annot` carry an annotations block
and no features, and they must hold the same residues as their bare
counterparts, because that pair is how the cost of the annotations dict is
isolated rather than argued about.
"""

from __future__ import annotations

import importlib.util
import io
import sys
from pathlib import Path

import pytest

from Bio import SeqIO

BENCH = Path(__file__).resolve().parent.parent / "bench"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        f"_bench_{name}", BENCH / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


seqio_corpus = _load("seqio_corpus")

ROWS = {row[0]: row
        for row in seqio_corpus.CORPORA + seqio_corpus.ANNOT_CORPORA}


def _parsed(group: str) -> list:
    """Parse one row, building only that row.

    ``corpus_rows()`` builds every row, which is five megabytes of GenBank and
    three of SwissProt; calling it once per assertion made this file the slowest
    in the suite for no reason.
    """
    _, fmt, count, length, per_kb, annot = ROWS[group]
    text = seqio_corpus.build(fmt, count, length, per_kb, annot)
    return list(SeqIO.parse(io.StringIO(text), fmt))


def test_every_row_builds_the_records_and_lengths_it_promises():
    for group, fmt, count, length, per_kb, annot in seqio_corpus.CORPORA:
        text = seqio_corpus.build(fmt, count, length, per_kb, annot)
        records = list(SeqIO.parse(io.StringIO(text), fmt))
        assert len(records) == count, group
        assert all(len(record.seq) == length for record in records), group


def test_the_ids_are_the_ones_the_generator_minted():
    for group, fmt, count, _, _, _ in seqio_corpus.CORPORA:
        records = _parsed(group)
        expected = [f"SYN{i + 1:06d}" for i in range(count)]
        assert [record.id for record in records] == expected, group


def test_the_first_sequence_is_the_one_the_keystream_produces():
    """The bare row's residues, recomputed from the seed and domain string.

    Recomputed rather than read back, so an emitter that transposed or dropped
    residues fails here instead of producing a corpus that parses perfectly and
    measures the wrong thing.  Only the bare row is checked this way: the
    featured rows draw their feature positions from the *same* keystream between
    residues, so their sequences are not a plain prefix of it, and recomputing
    them here would mean reimplementing the interleaving this test exists to
    check independently.
    """
    gen_data = _load("gen_data")
    count, length, per_kb = 8, 1_000, 0
    stream = gen_data.Keystream(seqio_corpus.SEQIO_SEED,
                                b"seqio/genbank/8x1000")
    records = list(SeqIO.parse(io.StringIO(
        seqio_corpus.build("genbank", count, length, per_kb)), "genbank"))
    for record in records:
        expected = stream.take(length).translate(seqio_corpus._DNA_TABLE)
        assert str(record.seq) == expected.decode("ascii")


def test_the_sequences_use_only_the_alphabet_the_format_allows():
    for group, fmt, count, length, per_kb, annot in seqio_corpus.CORPORA:
        records = _parsed(group)
        allowed = (set("ACGT") if fmt != "swiss"
                   else set("ACDEFGHIKLMNPQRSTVWY"))
        for record in records:
            assert set(str(record.seq)) <= allowed, group


def test_rebuilding_a_row_is_byte_identical():
    #   Kept separate from the round-trip test above because this is the property
    #   the ranking numbers depend on: if two builds differ, two machines differ.
    for group, fmt, count, length, per_kb, annot in seqio_corpus.CORPORA:
        first = seqio_corpus.build(fmt, count, length, per_kb, annot)
        second = seqio_corpus.build(fmt, count, length, per_kb, annot)
        assert first == second, group


def test_the_bare_row_has_no_features_and_the_featured_rows_do():
    """The decomposition in the ranking pass rests on this difference.

    The two rows are built from the same seed and the same domain string with
    the feature density as the only difference, so the same residues come out
    on both sides and the objects are the only thing the ranking pass measures
    between them.
    """
    featured = list(SeqIO.parse(io.StringIO(
        seqio_corpus.build("genbank", 20, 1_000, 1)), "genbank"))
    bare = list(SeqIO.parse(io.StringIO(
        seqio_corpus.build("genbank", 20, 1_000, 0)), "genbank"))
    assert all(len(record.features) > 0 for record in featured)
    assert all(len(record.features) == 0 for record in bare)
    #   Only the seed and the feature density differ, and the density is applied
    #   after the sequence is drawn, so the residues are the same on both sides.
    assert [str(r.seq) for r in bare] == [str(r.seq) for r in featured]


def test_the_nucleotide_records_carry_the_annotations_the_writer_was_given():
    records = list(SeqIO.parse(io.StringIO(
        seqio_corpus.build("genbank", 5, 500, 1)), "genbank"))
    for record in records:
        assert record.annotations["molecule_type"] == "DNA"
        assert record.annotations["data_file_division"] == "PLN"
        assert record.annotations["topology"] == "linear"
        assert record.annotations["organism"] == "synthetic construct"


def test_the_swiss_records_have_the_feature_blocks_they_were_written_with():
    """Two `DOMAIN` blocks plus the `CHAIN` span, per entry, as emitted."""
    records = list(SeqIO.parse(io.StringIO(
        seqio_corpus.build("swiss", 20, 300, 2)), "swiss"))
    for record in records:
        types = [feature.type for feature in record.features]
        assert types.count("CHAIN") == 1
        assert types.count("DOMAIN") == 2


def test_a_swiss_record_with_no_features_carries_only_the_chain():
    records = list(SeqIO.parse(io.StringIO(
        seqio_corpus.build("swiss", 5, 300, 0)), "swiss"))
    for record in records:
        assert [feature.type for feature in record.features] == ["CHAIN"]


def test_the_reference_round_trips_every_row_without_refusing():
    """The gate the ranking pass stands on: nothing in the corpus is malformed."""
    for group, fmt, count, _, _, _ in seqio_corpus.CORPORA:
        records = _parsed(group)
        assert len(records) == count
        assert all(record.seq for record in records), group
        assert all(record.id for record in records), group


ANNOTATED = [row[0] for row in seqio_corpus.ANNOT_CORPORA if row[5]]

#: The comment the annotated rows are built with, as one line: the writer wraps
#: it to the format's own width, so this is the text under the wrapping.
COMMENT = (
    "This record is generated for benchmarking and carries a comment long "
    "enough to span several lines of the flat file, which the reference "
    "reads back with the line breaks it wrote."
)


def test_the_annotated_rows_carry_the_blocks_their_bare_twin_lacks():
    """Six keys the thin header has not got, and the writer put them all back.

    The blocks are the whole point of the pair: an annotated row that quietly
    lost its references would still parse, still gate, and still measure the
    cheap header -- which is the failure this test exists to catch.
    """
    from Bio.SeqFeature import Reference

    assert ANNOTATED == ["genbank-1kb-annot", "embl-1kb-annot"]
    for group in ANNOTATED:
        records = _parsed(group)
        assert len(records) == 300, group
        for record in records:
            annotations = record.annotations
            assert annotations["taxonomy"] == [
                "artificial sequences", "synthetic construct"
            ], group
            assert annotations["keywords"] == ["synthetic", "benchmark"], group
            assert annotations["accessions"] == [record.id], group
            #   The comment comes back with the line breaks the writer chose,
            #   and the two formats wrap at different widths -- so the text is
            #   compared with its whitespace normalised and the *wrapping* is
            #   asserted separately, because "the writer wraps differently here"
            #   is a fact about the formats and not an accident of this row.
            assert " ".join(annotations["comment"].split()) == COMMENT, group
            assert "\n" in annotations["comment"], group
            assert len(annotations["references"]) == 3, group
            assert all(isinstance(ref, Reference) for ref in annotations["references"])
            assert annotations["references"][0].title == "Synthetic study 1"


def test_the_annotated_row_and_its_bare_twin_hold_the_same_residues():
    """The pair the annotations measurement rests on: header is the only delta.

    Both rows draw their residues from the same seed and the same domain string
    and the annotation blocks are applied after the draw, so this holds by
    construction -- and is asserted anyway, because "by construction" is what the
    feature pair was called the first time it was wrong.
    """
    for group, twin in (("genbank-1kb-annot", "genbank-1kb-bare"),
                        ("embl-1kb-annot", "embl-1kb-bare")):
        bare = _parsed(twin)
        annotated = _parsed(group)
        assert [str(r.seq) for r in annotated] == [str(r.seq) for r in bare], group
        assert [r.id for r in annotated] == [r.id for r in bare], group
        #   And the header really is the difference, in both directions.
        assert "references" not in bare[0].annotations, twin
        assert "references" in annotated[0].annotations, group
        assert all(len(record.features) == 0 for record in annotated), group
