"""The FEATURES table against the reference, feature for feature.

The feature half of a flat-file parse is the expensive half -- bench/targets.md
puts 16 of the 20.2 us that one feature per kilobase costs on reading the
location string and the qualifier block, the other 4.2 being `SeqFeature`
construction -- so it is the half worth a kernel of its own.  This file is the
gate that says the kernel reads the same table the reference does: every
feature's type, its location down to the structure of every position, and its
qualifiers key by value.

Two families of input, as in `tests/test_genbank.py`:

- the generated corpus (`bench/seqio_corpus.py`), which is what the ranking pass
  measured;
- shapes the corpus's writer never emits -- a location wrapped onto a second
  line, a quoted value continued across lines, a `/pseudo` with no value, an
  unescaped quote -- because a writer emits its own habits and a parser that has
  only ever met one writer's output has been tested against that writer.

The comparison is against the reference and never against a literal, and both
sides are reduced to the same nested tuples of plain integers (`test_location`'s
canonical form, imported here so that "the same location" has one definition)
before they are compared.  A difference in the answer must not be able to hide
behind a difference in the spelling of the answer.

**The refusals are the other half of the gate.**  The reference meets the shapes
listed in `DECLINED_MALFORMED` by warning and carrying on, and the ones in
`DECLINED_SCOPE` by reading them correctly under a header layout this reader
does not reproduce.  In both cases this reader refuses the *record* and says
which line did it, and a stricter reader that quietly stopped being stricter
would be a wrong feature that looks like a right one.  So the refusals are
asserted, and so is the fact that the reference did not refuse: a test that
called a shape malformed without checking would drift the day the reference
learned to read it.

The kernel also does not warn.  Warnings are a Python behaviour that belongs to
the layer that builds value types, so the kernel *reports* what it found -- the
origin-wrap repair, the dropped bond, the NCBI escaping -- and this file checks
each flag against the warning the reference emitted for the same input, both
ways round.
"""

from __future__ import annotations

import importlib.util
import io
import sys
import warnings
from pathlib import Path

import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import biofasting

# `test_location` owns the canonical form of a location: the same reduction of
# both sides to nested tuples of plain integers, for the same reason.  Imported
# rather than copied so that a change to what "the same location" means is one
# change in one place -- and imported by insertion because `tests/` is a
# directory of loose modules with no package name.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_location import _canonical_location  # noqa: E402

BENCH = Path(__file__).resolve().parent.parent / "bench"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"_bench_{name}", BENCH / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


seqio_corpus = _load("seqio_corpus")

# SwissProt's `FT` block is a different grammar (Bio.SwissProt._read_ft) and is
# not reproduced at all, so it is not a row this file can compare.
CORPUS_ROWS = [row[0] for row in seqio_corpus.CORPORA if row[1] != "swiss"]
_ROW_SPEC = {row[0]: row for row in seqio_corpus.CORPORA}
_CACHE: dict[str, tuple[str, str]] = {}


def _corpus(group: str) -> tuple[str, str]:
    """The (format, text) of one corpus row, built once and only that row."""
    if group not in _CACHE:
        _, format, count, length, per_kb, annot = _ROW_SPEC[group]
        _CACHE[group] = (format,
                         seqio_corpus.build(format, count, length, per_kb, annot))
    return _CACHE[group]


# --------------------------------------------------------------------------
# Comparing one file's features on both sides.
# --------------------------------------------------------------------------


def _reference(path: Path, format: str):
    """The reference's records and its warnings, from the file the kernel maps.

    Read from the file rather than from the text it was written from, so that a
    line ending the write and the read disagree about shows up as a mismatch
    instead of being hidden by both sides going through the same string.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with path.open() as handle:
            records = list(SeqIO.parse(handle, format))
    return records, [w.message for w in caught]


def _qualifiers(feature) -> dict[str, list[str]]:
    """The kernel's qualifier list as the reference's dict of value lists.

    `has_value` false is the `/pseudo` form, whose stored value the reference
    makes `""` -- once, however many times the key is written, which is why the
    grouping starts from the kernel's list rather than from the file.
    """
    out: dict[str, list[str]] = {}
    for key, value, _has_value, _escape, _text in feature[5]:
        out.setdefault(key, []).append(value)
    return out


def _our_features(index, record_id: str):
    ok, message, features = index.features(record_id)
    assert ok, f"{record_id}: the reader declined the table: {message}"
    return features


def _agree(tmp_path: Path, format: str, text: str, name: str = "in.txt") -> Path:
    """Write `text`, then hold the kernel's feature table to the reference's."""
    path = tmp_path / name
    path.write_text(text, newline="")
    records, messages = _reference(path, format)
    index = biofasting.open_genbank(path, format=format)

    for record in records:
        features = _our_features(index, record.id)
        assert len(features) == len(record.features), record.id
        for ours, expected in zip(features, record.features):
            feature_type, location, status, detail, warnings, qualifiers = ours
            assert feature_type == expected.type, record.id
            assert status == "ok", f"{record.id}: {detail}"
            assert location == _canonical_location(expected.location), record.id
            assert _qualifiers(ours) == expected.qualifiers, record.id
            # The location's in-parse warnings, as `(kind, text)` pairs: one
            # entry per offending part, because the reference warns once per
            # part.  Here there are none.
            assert warnings == [], record.id

    # Nothing in these inputs is malformed, so the reference has nothing to warn
    # about and the kernel has nothing to report.
    assert messages == []
    return path


# --------------------------------------------------------------------------
# The generated corpus: what the kernel was measured against.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("group", CORPUS_ROWS)
def test_corpus_row_features_match_the_reference(tmp_path, group):
    format, text = _corpus(group)
    _agree(tmp_path, format, text, name=f"{group}.txt")


def test_the_bare_row_has_an_empty_table_and_not_a_refusal(tmp_path):
    """"No features" and "features not reproduced" are different statements.

    `genbank-1kb-bare` exists to isolate the cost of building `SeqFeature`
    objects, so it is also the row that says an absent FEATURES block is an
    answer with an empty list rather than a refusal a caller has to handle.
    """
    format, text = _corpus("genbank-1kb-bare")
    path = tmp_path / "bare.gb"
    path.write_text(text, newline="")
    index = biofasting.open_genbank(path, format=format)
    for record in _reference(path, format)[0]:
        ok, message, features = index.features(record.id)
        assert ok and message == "" and features == []


# --------------------------------------------------------------------------
# The shapes the corpus's writer never emits.
# --------------------------------------------------------------------------

_RESIDUES = "ATGGCATGCATCGATCGATCGTAGCTAGCTAGCTAGCATCGATCGATCGTAGCATCGATCGTAGCATCGATCG"
_ANNOTATIONS = {
    "molecule_type": "DNA",
    "data_file_division": "PLN",
    "date": "01-JAN-2026",
    "topology": "linear",
}


def _written(format: str = "genbank", *, topology: str = "linear") -> str:
    """One record as the reference's own writer emits it.

    The header block's columns are fixed and unforgiving, so the skeleton comes
    from the writer and the shapes under test are spliced into it -- a
    hand-written `LOCUS` line is how a suite ends up testing its author's
    arithmetic instead of the format.
    """
    record = SeqRecord(Seq(_RESIDUES), id="TEST001", name="TEST001", description="d")
    record.annotations = dict(_ANNOTATIONS, topology=topology)
    handle = io.StringIO()
    SeqIO.write([record], handle, format)
    return handle.getvalue()


def _splice(text: str, old: str, new: str) -> str:
    """Replace one line, asserting it was there -- a splice that missed tests nothing."""
    assert old in text, f"fixture does not contain {old!r}"
    return text.replace(old, new, 1)


_GENBANK_HEAD = "FEATURES             Location/Qualifiers\n"
_EMBL_HEAD = "FH   Key             Location/Qualifiers\nFH\n"


def _feature(text: str, *lines: str, head: str = _GENBANK_HEAD) -> str:
    """Splice a feature block in just after the header line."""
    assert head in text, "fixture does not contain the feature header"
    return text.replace(head, head + "".join(lines), 1)


def test_a_location_wrapped_onto_a_second_line(tmp_path):
    """The comma at the end of the line is what says the location continues."""
    text = _feature(
        _written(),
        "     CDS             join(1..10,20..30,\n",
        "                     40..50)\n",
        "                     /product=\"x\"\n",
    )
    _agree(tmp_path, "genbank", text)


def test_a_quoted_value_continued_across_lines(tmp_path):
    """`translation` is the one value the reference strips whitespace out of.

    The lines are joined with newlines by the scanner and with spaces by the
    consumer, and both are removed afterwards -- so the value is the residues as
    one run of letters, which is what a caller compares against a protein.
    """
    text = _feature(
        _written(),
        "     CDS             1..30\n",
        "                     /translation=\"MKV\n",
        "                     LAA\"\n",
        "                     /note=\"a b\"\n",
    )
    path = _agree(tmp_path, "genbank", text)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    assert features[0][5][0][1] == "MKVLAA"
    # A qualifier that is not `translation` keeps its internal whitespace.
    assert features[0][5][1][1] == "a b"


def test_a_valueless_qualifier_and_its_repeat(tmp_path):
    """`/pseudo` stores `""` once and the second one is dropped outright."""
    text = _feature(
        _written(),
        "     gene            1..30\n",
        "                     /pseudo\n",
        "                     /pseudo\n",
    )
    path = _agree(tmp_path, "genbank", text)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    assert [(key, value, has_value) for key, value, has_value, _e, _t
            in features[0][5]] == [("pseudo", "", False)]


def test_an_empty_qualifier_value(tmp_path):
    """`/note=` is an empty value, which is not the same fact as no value."""
    text = _feature(_written(), "     gene            1..30\n", "                     /note=\n")
    path = _agree(tmp_path, "genbank", text)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    assert [(key, value, has_value) for key, value, has_value, _e, _t
            in features[0][5]] == [("note", "", True)]


def test_an_unquoted_value_continued_on_the_next_line(tmp_path):
    text = _feature(
        _written(),
        "     gene            1..30\n",
        "                     /note=first\n",
        "                     second\n",
    )
    path = _agree(tmp_path, "genbank", text)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    # The scanner joins the lines with a newline and the consumer turns each
    # into a space, so an unquoted continuation keeps one space per line break.
    assert features[0][5][0][1] == "first second"


def test_escaped_and_unescaped_double_quotes(tmp_path):
    """The NCBI escaping is undone, and the warning is about the value *before*.

    This is the one qualifier behaviour that is a warning rather than a value,
    so the kernel reports it instead of emitting it: `escape_text` is the value
    as the warning's own `%r` would have shown it.
    """
    escaped = _feature(
        _written(),
        "     gene            1..30\n",
        "                     /note=\"say \"\"hi\"\" now\"\n",
    )
    path = _agree(tmp_path, "genbank", escaped)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    assert features[0][5][0][1] == 'say "hi" now'
    assert features[0][5][0][3] is False

    unescaped = _feature(
        _written(),
        "     gene            1..30\n",
        "                     /note=\"say \"hi\" now\"\n",
    )
    path = tmp_path / "unescaped.gb"
    path.write_text(unescaped, newline="")
    records, messages = _reference(path, "genbank")
    assert len(messages) == 1 and "should be escaped as" in str(messages[0])
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    key, value, _has, escape, text = features[0][5][0]
    assert escape is True
    # The warning quotes the value as it stood before the doubled quotes came
    # undone, which is what makes it worth carrying rather than recomputing.
    assert repr(text) in str(messages[0])
    # And the value is the same either way: the warning is not about the answer.
    assert value == records[0].features[0].qualifiers["note"][0]


def test_a_compound_location_and_a_repeated_qualifier_key(tmp_path):
    text = _feature(
        _written(),
        "     CDS             complement(join(1..10,20..30))\n",
        "                     /db_xref=\"GI:1\"\n",
        "                     /db_xref=\"GI:2\"\n",
    )
    path = _agree(tmp_path, "genbank", text)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    assert _qualifiers(features[0]) == {"db_xref": ["GI:1", "GI:2"]}


def test_the_replace_repair_is_reproduced_as_written(tmp_path):
    """`replace(266,"c")` becomes the number it names, off-by-one in the slice included.

    The cut is `location_line[8:comma_pos]` in Python, which with no comma at
    all is a slice to one before the end -- reproduced rather than tidied,
    because it is what the reference hands to its parser.
    """
    text = _feature(_written(), "     gene            replace(266,\"c\")\n")
    path = _agree(tmp_path, "genbank", text)
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    expected = _reference(path, "genbank")[0][0].features[0].location
    assert features[0][2] == "ok" and features[0][1] == _canonical_location(expected)


def test_a_closing_parenthesis_on_a_line_of_its_own(tmp_path):
    """The reference appends it and hands its parser a broken string.

    This branch is reachable only when the location is already balanced: an
    unbalanced one is eaten by the line-wrapping warning first, and that shape
    is refused (see `DECLINED_MALFORMED`).  So what the branch produces is a
    `LocationParserError`, which is a behaviour and not a failure -- the
    location becomes `None`, the record keeps its other features, and the
    kernel says so with `parser_error` rather than raising.
    """
    text = _feature(
        _written(),
        "     CDS             join(1..10,\n",
        "                     20..30)\n",
        "                     )\n",
        "     gene            1..20\n",
    )
    path = tmp_path / "paren.gb"
    path.write_text(text, newline="")
    records, messages = _reference(path, "genbank")
    assert records[0].features[0].location is None
    assert any("Could not parse feature location" in str(m) for m in messages)

    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    assert len(features) == 2
    assert features[0][2] == "parser_error" and features[0][1] is None
    assert features[1][1] == _canonical_location(records[0].features[1].location)


def test_a_bond_location_reports_the_warning_it_would_have_raised(tmp_path):
    """A top-level `bond(...)` is an operator and does not warn; one inside a
    `join` reaches `SimpleLocation`'s bond category and does."""
    operator = _feature(_written(), "     misc_feature    bond(196)\n")
    _agree(tmp_path, "genbank", operator)

    nested = _feature(_written(), "     misc_feature    join(bond(196),197..200)\n")
    path = tmp_path / "bond.gb"
    path.write_text(nested, newline="")
    records, messages = _reference(path, "genbank")
    assert [str(m) for m in messages] == ["Dropping bond qualifier in feature location"]
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    # One per offending part.  The bond warning's wording is fixed -- it quotes
    # nothing -- so the text is empty and the kind is the whole fact.
    assert features[0][4] == [("bond", "")]
    assert features[0][1] == _canonical_location(records[0].features[0].location)


def test_an_origin_wrapping_location_on_a_circular_record(tmp_path):
    """The reference repairs it, warns, and keeps the repaired location."""
    text = _feature(
        _written(topology="circular"), "     gene            30..5\n"
    )
    path = tmp_path / "circular.gb"
    path.write_text(text, newline="")
    records, messages = _reference(path, "genbank")
    assert len(messages) == 1 and str(messages[0]).startswith("Attempting to fix")
    index = biofasting.open_genbank(path)
    _ok, _message, features = index.features("TEST001")
    # The text travels with the fact: the reference's warning quotes the part it
    # repaired, so the layer can say the words without parsing the location a
    # second time to find out what it was about.
    assert features[0][4] == [("origin_wrap", "30..5")]
    assert features[0][1] == _canonical_location(records[0].features[0].location)


def test_embl_feature_blocks_are_read_by_the_same_reader(tmp_path):
    """`FT` instead of a column of spaces, and the rest of the grammar identical."""
    text = _feature(
        _written("embl"),
        "FT   CDS             join(1..10,20..30)\n",
        "FT                   /product=\"x\"\n",
        "FT                   /translation=\"MK\n",
        "FT                   VLAA\"\n",
        head=_EMBL_HEAD,
    )
    path = _agree(tmp_path, "embl", text, name="in.embl")
    index = biofasting.open_genbank(path, format="embl")
    _ok, _message, features = index.features("TEST001")
    assert _qualifiers(features[0]) == {"product": ["x"], "translation": ["MKVLAA"]}


# --------------------------------------------------------------------------
# The refusals: shapes the reference only warns about.
# --------------------------------------------------------------------------

# Each entry is (what the shape is, the block to splice in).  Every one of them
# is a *value* the reference guesses at after a warning, so reproducing it
# differently would be a feature that looks right and is not.
DECLINED_MALFORMED = [
    (
        "a location that wraps its parentheses without breaking at a comma",
        "     CDS             complement(join(1..10,\n"
        "                     20..30)\n"
        "                     )\n",
    ),
    (
        "a line too short to hold a feature key and a location",
        "     CDS\n",
    ),
    (
        "an over-indented location column",
        "     CDS             join(1..10, 20..30)\n",
    ),
    (
        "white space between a qualifier's '=' and its opening quote",
        "     gene            1..10\n"
        "                     /note= \"x\"\n",
    ),
    (
        "a continuation line with no qualifier above it",
        "     gene            1..10\n"
        "                     loose text\n",
    ),
    (
        "a continuation line after a valueless qualifier",
        "     gene            1..10\n"
        "                     /pseudo\n"
        "                     more\n",
    ),
]


@pytest.mark.parametrize("label,block", DECLINED_MALFORMED, ids=[d[0] for d in DECLINED_MALFORMED])
def test_a_shape_the_reference_only_warns_about_is_refused(tmp_path, label, block):
    text = _feature(_written(), block)
    path = tmp_path / "malformed.gb"
    path.write_text(text, newline="")
    data = path.read_bytes()

    # The premise: the reference does not read this cleanly.  It either warns and
    # carries on or raises out of the parse altogether, and a test that assumed
    # so without checking would be the thing that goes stale.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            list(SeqIO.parse(io.StringIO(text), "genbank"))
            raised = False
        except Exception:  # noqa: BLE001 - the taxonomy is the point
            raised = True
    assert raised or caught, f"{label}: the reference read it cleanly after all"

    index = biofasting.open_genbank(path)
    ids = list(index)
    assert ids
    for record_id in ids:
        ok, message, features = index.features(record_id)
        assert not ok, label
        assert features == []
        # The refusal has to be actionable: it names a byte of *this* file, at
        # the start of a line, so a caller can quote the offending line back.
        offset = int(message.split("byte offset ")[1].split(":")[0])
        assert 0 < offset < len(data)
        assert data[offset - 1 : offset] == b"\n"


# These are refused for scope and not for malformation: the reference reads them
# correctly, under a header whose layout this reader does not reproduce.  The
# locations would then be parsed under a declared size and topology the
# reference never used, so the reader declines instead of guessing at them.
DECLINED_SCOPE = [
    (
        "a LOCUS line in the seven-field space-separated layout",
        "LOCUS       HG531_PATCH 1000000 bp DNA HTG 18-JUN-2011",
    ),
]


@pytest.mark.parametrize("label,locus", DECLINED_SCOPE, ids=[d[0] for d in DECLINED_SCOPE])
def test_a_header_layout_this_reader_cannot_read_is_refused(tmp_path, label, locus):
    text = _feature(_written(), "     gene            1..10\n")
    written_locus = next(line for line in text.splitlines() if line.startswith("LOCUS"))
    text = _splice(text, written_locus, locus)
    path = tmp_path / "layout.gb"
    path.write_text(text, newline="")

    with warnings.catch_warnings(record=True):
        warnings.simplefilter("ignore")
        records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    assert len(records) == 1

    index = biofasting.open_genbank(path)
    ok, message, features = index.features("TEST001")
    assert not ok and features == []
    assert "not in a layout this reader reproduces" in message


def test_a_bare_record_end_inside_the_table_fails_the_index(tmp_path):
    """`//` inside the feature table is a `ValueError` out of the reference too.

    The block scan and the reader share one terminator rule, so a `//` that ends
    the record early is caught by the scan, before any feature is read.
    """
    text = _feature(_written(), "     gene            1..10\n//\n")
    path = tmp_path / "truncated.gb"
    path.write_text(text, newline="")
    with pytest.raises(ValueError, match="Premature end of features table"):
        biofasting.open_genbank(path)


def test_the_swissprot_block_is_refused_rather_than_reported_empty(tmp_path):
    """`FT` is a different grammar -- `Bio.SwissProt._read_ft` reads it.

    So there is no reader here to hand it to, and an empty list would be a
    statement about the record rather than about this package.
    """
    from test_genbank import _swiss

    path = tmp_path / "in.dat"
    path.write_text(_swiss(), newline="")
    index = biofasting.open_genbank(path, format="swiss")
    ok, message, features = index.features("TEST001")
    assert not ok and features == []
    assert "different grammar" in message


# --------------------------------------------------------------------------
# The contract the rest of the package rests on.
# --------------------------------------------------------------------------


def test_reading_features_leaves_the_rest_of_the_index_alone(tmp_path):
    """The table is read on demand, so asking for it changes nothing else.

    `open_genbank`'s records are the same four fields before and after, and a
    second call answers the same table: an index whose answers depend on what
    has already been asked of it would make every measurement order-dependent.
    """
    format, text = _corpus("genbank-1kb")
    path = tmp_path / "genbank.gb"
    path.write_text(text, newline="")
    index = biofasting.open_genbank(path, format=format)

    before = index.records()
    first = [index.features(record_id) for record_id in index]
    second = [index.features(record_id) for record_id in index]
    assert first == second
    assert index.records() == before


def test_the_reader_reports_warnings_rather_than_emitting_them(tmp_path):
    """A kernel that warned would be a kernel deciding what Python does with it.

    The suite runs with `filterwarnings = ["error"]`, so an emitted warning
    would already be a failure here; this says so on purpose, on the three
    inputs that have something to warn about.
    """
    for block in (
        "     misc_feature    join(bond(196),197..200)\n",
        "     gene            1..30\n"
        "                     /note=\"say \"hi\" now\"\n",
    ):
        text = _feature(_written(), block)
        path = tmp_path / "warns.gb"
        path.write_text(text, newline="")
        index = biofasting.open_genbank(path)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ok, _message, features = index.features("TEST001")
        assert ok
        assert any(qualifier[3] for qualifier in features[0][5]) or features[0][4]

    text = _feature(_written(topology="circular"), "     gene            30..5\n")
    path = tmp_path / "wrap.gb"
    path.write_text(text, newline="")
    index = biofasting.open_genbank(path)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ok, _message, features = index.features("TEST001")
    assert ok and features[0][4] == [("origin_wrap", "30..5")]


def test_a_query_that_cannot_parse_a_location_still_answers(tmp_path):
    """`LocationParserError` is a behaviour, not a failure.

    The reference catches it inside its feature consumer and turns it into
    `location = None` plus a warning, so the record keeps its other features.
    The kernel carries the status out rather than raising, which is what lets
    the Python layer reproduce the warning instead of losing the difference.
    """
    # A linear record, so a start past its end is not origin wrapping but a
    # location the reference's own parser rejects.
    text = _feature(
        _written(),
        "     gene            30..5\n",
        "     CDS             1..10\n",
    )
    path = tmp_path / "unparseable.gb"
    path.write_text(text, newline="")
    records, messages = _reference(path, "genbank")
    assert records[0].features[0].location is None
    assert any("setting feature location to None" in str(m) for m in messages)

    index = biofasting.open_genbank(path)
    ok, _message, features = index.features("TEST001")
    assert ok and len(features) == 2
    assert features[0][2] == "parser_error" and features[0][1] is None
    # And the feature after it is still read: one bad location is one feature.
    assert features[1][2] == "ok"
    assert features[1][1] == _canonical_location(records[0].features[1].location)
