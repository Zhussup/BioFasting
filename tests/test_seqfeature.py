"""`biofasting.seqfeature`: our values, and the conversion to `Bio.SeqFeature`.

The kernel below this layer is gated by `tests/test_features.py`, which compares
its tuples against the reference's own `SeqFeature` objects field for field.
What that file cannot test is the layer *above*: the value types a caller holds,
and the conversion that turns them back into Biopython objects -- and, in
particular, the one thing a kernel must not do, which is emit a warning.  The
reader reports what it found; this module says the words, and that the words are
the reference's own is what this file checks.

So the comparison here is deliberately a different one: the *object* our
conversion builds, reduced to the same canonical form `test_location` defines
(so that `==` on two fuzzy positions, which compares by integer value, cannot
hide a difference in kind), the qualifier dict the reference stores, and the
warnings both sides emit -- message for message, in order, because the reference
parses the location before it feeds the qualifiers and a feature with both kinds
of warning says them in that order.

The corpus is the same one the ranking pass measured, and the spliced fixtures
are the shapes its writer never emits: a location that wraps its origin, a
dropped `bond`, an NCBI escaping, and a `LocationParserError` the reference
catches and turns into a missing location.
"""

from __future__ import annotations

import io
import sys
import warnings
from pathlib import Path

import pytest

from Bio import SeqIO

import biofasting

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_location import _canonical_location  # noqa: E402
from test_features import (  # noqa: E402
    CORPUS_ROWS,
    _corpus,
    _feature,
    _reference,
    _written,
)


def _ours(path: Path, record_id: str, format: str = "genbank"):
    """(features, warnings) -- our conversion, and what it said while doing it."""
    index = biofasting.open_genbank(path, format)
    features = biofasting.read_features(index, record_id)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        converted = [biofasting.to_seqfeature(feature) for feature in features]
    return list(zip(features, converted)), [str(w.message) for w in caught]


# --------------------------------------------------------------------------
# The corpus: what the kernel was measured against.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("group", CORPUS_ROWS)
def test_corpus_features_survive_the_conversion(tmp_path, group):
    """Every field, through the value types and back out as a `SeqFeature`.

    One index and one warning filter for the whole file rather than per record:
    the corpus rows are thousands of records, and a per-record `open_genbank`
    would measure the test's own bookkeeping.
    """
    format, text = _corpus(group)
    path = tmp_path / f"{group}.txt"
    path.write_text(text, newline="")
    records, reference_warnings = _reference(path, format)
    reference_warnings = [str(message) for message in reference_warnings]

    index = biofasting.open_genbank(path, format)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for record in records:
            features = biofasting.read_features(index, record.id)
            assert len(features) == len(record.features), record.id
            for feature, expected in zip(features, record.features):
                converted = biofasting.to_seqfeature(feature)
                assert converted.type == expected.type, record.id
                assert _canonical_location(converted.location) == _canonical_location(
                    expected.location
                ), record.id
                assert converted.qualifiers == expected.qualifiers, record.id
                # And the value type carries the same three things the object does.
                assert feature.type == expected.type
                assert feature.qualifiers_as_dict() == expected.qualifiers
                assert feature.status == "ok"
                assert feature.warnings == ()
    assert [str(w.message) for w in caught] == reference_warnings == []


# --------------------------------------------------------------------------
# The value types, where they say something the object does not.
# --------------------------------------------------------------------------


def test_a_position_keeps_the_kind_of_fuzziness_it_carries(tmp_path):
    """`<`, `>`, `(a.b)`, `one-of(...)` and `?` are five different statements and
    the value type keeps them apart, because `==` between two positions in
    Biopython compares integers and would not."""
    text = _feature(
        _written(),
        "     misc_feature    <1..>30\n",
        "     misc_feature    (3.9)..20\n",
        "     misc_feature    one-of(5,7)..10\n",
    )
    path = tmp_path / "kinds.gb"
    path.write_text(text, newline="")
    features = biofasting.read_features(biofasting.open_genbank(path), "TEST001")

    assert [feature.location.parts[0].start.kind for feature in features] == [
        "before",
        "within",
        "one_of",
    ]
    # `(3.9)` is a start boundary, so the reference reads it with offset -1:
    # both edges are one lower.
    assert (features[1].location.parts[0].start.left,
            features[1].location.parts[0].start.right) == (2, 8)
    assert features[2].location.parts[0].start.choices == (4, 6)

    # The conversion builds the matching Biopython class, not a plain int.
    from Bio.SeqFeature import BeforePosition, OneOfPosition, WithinPosition

    built = [biofasting.to_seqfeature(feature).location.start for feature in features]
    assert isinstance(built[0], BeforePosition)
    assert isinstance(built[1], WithinPosition)
    assert isinstance(built[2], OneOfPosition)


def test_strand_zero_is_the_references_none(tmp_path):
    """A protein record's location has no strand at all, and 0 is how the kernel
    spells that `None` -- so the conversion has to put `None` back, and not a
    strand of "unknown"."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("MKVLAAGIVGL"), id="TEST001", name="TEST001", description="d")
    record.annotations = {"molecule_type": "protein", "data_file_division": "PLN"}
    handle = io.StringIO()
    SeqIO.write([record], handle, "embl")
    # The ID line says `PROTEIN`, which is the word the reference's `stranded`
    # test looks for and the reason this record's locations have no strand.
    text = _feature(
        handle.getvalue(),
        "FT   CHAIN           1..11\n",
        head="FH   Key             Location/Qualifiers\nFH\n",
    )
    path = tmp_path / "protein.embl"
    path.write_text(text, newline="")
    records, _ = _reference(path, "embl")
    assert records[0].features[0].location.strand is None

    features = biofasting.read_features(biofasting.open_genbank(path, "embl"), "TEST001")
    assert features
    assert all(feature.location.parts[0].strand == 0 for feature in features)
    for feature in features:
        assert biofasting.to_seqfeature(feature).location.strand is None


def test_a_single_part_join_is_a_simple_location(tmp_path):
    """The reference returns the part itself, not a one-part compound, and a
    conversion that always built a `CompoundLocation` would raise."""
    text = _feature(_written(), "     gene            join(1..30)\n")
    path = tmp_path / "one.gb"
    path.write_text(text, newline="")
    records, _ = _reference(path, "genbank")
    features = biofasting.read_features(biofasting.open_genbank(path), "TEST001")

    assert features[0].location.operator == "simple"
    assert len(features[0].location.parts) == 1
    converted = biofasting.to_seqfeature(features[0])
    assert type(converted.location) is type(records[0].features[0].location)


def test_a_qualifier_dict_is_grouped_and_not_decided(tmp_path):
    """A repeated key becomes a list of values, and the valueless form is the
    empty string -- which the reader has already resolved, including the
    reference's rule that a *repeat* of a valueless key is dropped."""
    text = _feature(
        _written(),
        "     gene            1..30\n",
        "                     /pseudo\n",
        "                     /pseudo\n",
        "                     /db_xref=\"A\"\n",
        "                     /db_xref=\"B\"\n",
    )
    path = tmp_path / "dict.gb"
    path.write_text(text, newline="")
    records, _ = _reference(path, "genbank")
    feature = biofasting.read_features(biofasting.open_genbank(path), "TEST001")[0]

    assert feature.qualifiers_as_dict() == records[0].features[0].qualifiers
    assert feature.qualifiers_as_dict() == {"pseudo": [""], "db_xref": ["A", "B"]}
    assert [qualifier.has_value for qualifier in feature.qualifiers] == [False, True, True]


def test_a_record_without_features_is_an_empty_list(tmp_path):
    """"No features" and "features not reproduced" are different statements, and
    only the second one is a refusal."""
    path = tmp_path / "bare.gb"
    path.write_text(_written(), newline="")
    index = biofasting.open_genbank(path)
    assert biofasting.read_features(index, "TEST001") == []


def test_a_declined_table_is_an_error_that_names_the_line(tmp_path):
    """The reader refuses the record it cannot reproduce; the value-type layer
    turns that refusal into an exception rather than an empty list, because an
    empty list is what a record with no features returns."""
    text = _feature(_written(), "     CDS\n")
    path = tmp_path / "short.gb"
    path.write_text(text, newline="")
    with pytest.raises(ValueError) as caught:
        biofasting.read_features(biofasting.open_genbank(path), "TEST001")
    assert "TEST001" in str(caught.value) and "byte offset" in str(caught.value)


# --------------------------------------------------------------------------
# The warnings: the half a kernel must not do.
# --------------------------------------------------------------------------

# Each entry is (label, the record's topology, the feature block).  The
# reference warns about every one of them, so the conversion has to say exactly
# the same words, in the same order.
WARNING_SHAPES = [
    ("an origin wrap", "circular", "     gene            30..5\n"),
    ("a dropped bond", "linear", "     misc_feature    join(bond(196),197..200)\n"),
    (
        "an unescaped double quote",
        "linear",
        '     gene            1..30\n                     /note="say "hi" now"\n',
    ),
    (
        "a location the parser rejects",
        "linear",
        "     gene            123^456\n",
    ),
    (
        "both kinds, in the reference's order",
        "circular",
        "     misc_feature    join(bond(1),30..5)\n",
    ),
]


@pytest.mark.parametrize(
    "label,topology,block", WARNING_SHAPES, ids=[shape[0] for shape in WARNING_SHAPES]
)
def test_the_conversion_says_the_references_own_words(tmp_path, label, topology, block):
    text = _feature(_written(topology=topology), block)
    path = tmp_path / "warns.gb"
    path.write_text(text, newline="")
    records, reference_warnings = _reference(path, "genbank")
    reference_warnings = [str(message) for message in reference_warnings]
    assert reference_warnings, f"{label}: the reference had nothing to warn about"

    pairs, our_warnings = _ours(path, records[0].id)
    assert our_warnings == reference_warnings, label

    # And the feature came through it: a parser error is a *behaviour*, so the
    # feature is still there with its type, just without a location.
    converted = pairs[0][1]
    assert converted.type == records[0].features[0].type, label
    expected = records[0].features[0].location
    assert (converted.location is None) == (expected is None), label
    if expected is not None:
        assert _canonical_location(converted.location) == _canonical_location(
            expected
        ), label


def test_a_clean_table_says_nothing(tmp_path):
    """The suite runs with `filterwarnings = ["error"]`, but a test that only
    relied on that would not say which input was silent."""
    text = _feature(
        _written(),
        "     CDS             join(1..10,20..30)\n",
        '                     /product="x"\n',
    )
    path = tmp_path / "clean.gb"
    path.write_text(text, newline="")
    records, reference_warnings = _reference(path, "genbank")
    assert reference_warnings == []
    _pairs, our_warnings = _ours(path, records[0].id)
    assert our_warnings == []


def test_reading_the_values_emits_nothing_at_all(tmp_path):
    """Reading features and holding them is silent: only the Biopython
    conversion says the words, so a caller who never builds a `SeqFeature`
    never has to install a warning filter."""
    text = _feature(_written(topology="circular"), "     gene            30..5\n")
    path = tmp_path / "silent.gb"
    path.write_text(text, newline="")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        feature = biofasting.read_features(biofasting.open_genbank(path), "TEST001")[0]
    assert feature.warnings == (("origin_wrap", "30..5"),)
