"""The annotations dict against the reference's `SeqRecord.annotations`.

A flat-file record is three things, and this file is the gate on the third: the
sequence, which `tests/test_genbank.py` holds to the reference; the FEATURES
table, which `tests/test_features.py` does; and the *header*, which
`SeqIO.parse` puts in `SeqRecord.annotations`.  M22 built the kernel for it, and
the measurement that justifies it is in `bench/targets.md`: a GenBank record
with an ordinary three-line comment costs the reference 16.4 us more than the
same record without one, and almost none of that is reading a comment -- it is a
structured-comment probe that is quadratic in the line it fails on, paid by
every record that has one.

**The comparison is the whole dict, keys and order included.**  Two dicts are
equal under `==` whatever order they were built in, so an equality check alone
would miss the thing this reader is most likely to get wrong: the reference
inserts a key when the line stating it is *consumed*, so a `COMMENT` above the
first `REFERENCE` puts `comment` before `references` and one below puts it
after.  `list(...)` of the keys is therefore asserted before the values are, and
the two spliced shapes `COMMENT above`/`COMMENT below` exist to make that
assertion mean something.

**The two formats do not have the same keys**, and the difference is a property
of the reference's consumer tables rather than of any one file: EMBL has no
`date` (its scanner reads no `DT` line at all), no `source`, and no `keywords`
key unless a `KW` line was there, where GenBank's writer emits placeholders that
become `''` and `['']`.  A normalising reader would be right about one format.
The corpus rows are generated per format, so they carry the difference; the
named tests below state it outright so that a future normalisation fails loudly
rather than silently.

**The refusals are the other half.**  A header line that sets a key this reader
does not produce makes the *record's* annotations refuse and name the shape it
choked on -- the same contract the features table has, and for the same reason:
a dict that is wrong is worse than a dict that is absent.  The refusals are
asserted together with the fact that the reference read the same input happily,
because a test that called a shape unreadable without checking would drift the
day the reference learned to read it.

Both sides are compared as plain values -- strings, lists and the eight fields
of a `Reference` -- so a difference in the answer cannot hide behind a
difference in the spelling of the answer.
"""

from __future__ import annotations

import importlib.util
import io
import sys
from pathlib import Path

import pytest

from Bio import SeqIO

import biofasting

BENCH = Path(__file__).resolve().parent.parent / "bench"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"_bench_{name}", BENCH / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


seqio_corpus = _load("seqio_corpus")

#: The rows this reader reproduces.  SwissProt's header is `Bio.SwissProt`, a
#: different reader with different keys, and is out of scope by the header's own
#: statement -- a row that is out of scope is not a row that may be skipped
#: silently, which is why it is named here rather than filtered inside a helper.
CORPUS_ROWS = [
    row[0] for row in seqio_corpus.CORPORA + seqio_corpus.ANNOT_CORPORA
    if row[1] in ("genbank", "embl")
]


def _corpus(group: str) -> tuple[str, str]:
    format, _ = seqio_corpus.build_row(group)
    _, fmt, count, length, per_kb, annot = seqio_corpus.spec(group)
    return format, seqio_corpus.build(fmt, count, length, per_kb, annot)


def _reference(path: Path, format: str):
    """The reference's records for a file, with warnings turned into errors."""
    with path.open() as handle:
        return list(SeqIO.parse(handle, format))


def _ours(index, record_id: str) -> dict:
    return biofasting.read_annotations(index, record_id)


def _compare(ours: dict, expected: dict, record_id: str) -> None:
    """Keys, their order, and every value -- references through `to_reference`.

    The order is asserted first and separately because `==` on two dicts is
    order-blind, and the order is the part of this dict that a fixed key list
    would get wrong.
    """
    assert list(ours) == list(expected), f"{record_id}: key order"
    for key in expected:
        if key == "references":
            assert len(ours[key]) == len(expected[key]), f"{record_id}: {key}"
            for i, (mine, theirs) in enumerate(zip(ours[key], expected[key])):
                assert biofasting.to_reference(mine) == theirs, (
                    f"{record_id}: reference[{i}]"
                )
        else:
            assert ours[key] == expected[key], f"{record_id}: {key}"


def _agree(tmp_path: Path, format: str, text: str, name: str = "in.txt") -> Path:
    """Write `text`, then hold the kernel's annotations to the reference's."""
    path = tmp_path / name
    path.write_text(text, newline="")
    records = _reference(path, format)
    index = biofasting.open_genbank(path, format=format)
    for record in records:
        _compare(_ours(index, record.id), record.annotations, record.id)
    return path


def _corpus_bases() -> tuple[str, str]:
    """One real record per format, to splice header lines into.

    A GenBank `LOCUS` line is fixed-width and hand-typing one produces a shape
    the reference itself warns about, so the spliced shapes are built by
    replacing lines in a record the corpus generator wrote.
    """
    return (
        seqio_corpus.build("genbank", 1, 60, 0, 0),
        seqio_corpus.build("embl", 1, 60, 0, 0),
    )


# --------------------------------------------------------------------------
# The generated corpus: what the kernel was measured against.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("group", CORPUS_ROWS)
def test_corpus_row_annotations_match_the_reference(tmp_path, group):
    format, text = _corpus(group)
    _agree(tmp_path, format, text, name=f"{group}.txt")


def test_the_annotation_dense_row_has_more_keys_than_its_bare_twin(tmp_path):
    """"The dict matched" is only worth something if there was a dict.

    The thin rows carry a handful of keys; the annotation-dense rows carry a
    taxonomy, a keyword list, two accessions, a sequence version, a multi-line
    comment and three references.  If the dense rows silently degenerated into
    thin ones, every comparison above would still pass -- so the difference is
    asserted here rather than left implied by which rows exist.
    """
    _, bare = _corpus("genbank-1kb-bare")
    _, dense = _corpus("genbank-1kb-annot")
    assert bare.count("COMMENT") == 0 and dense.count("COMMENT") > 0

    bare_path = tmp_path / "bare.gb"
    dense_path = tmp_path / "dense.gb"
    bare_path.write_text(bare, newline="")
    dense_path.write_text(dense, newline="")

    bare_record = _reference(bare_path, "genbank")[0]
    dense_record = _reference(dense_path, "genbank")[0]
    assert len(dense_record.annotations) > len(bare_record.annotations)
    assert "references" in dense_record.annotations
    assert len(dense_record.annotations["references"]) == 3

    # And the kernel agrees with both, not just with the thin one.
    _agree(tmp_path, "genbank", bare, name="bare2.gb")
    _agree(tmp_path, "genbank", dense, name="dense2.gb")


# --------------------------------------------------------------------------
# The per-format key sets, stated rather than implied.
# --------------------------------------------------------------------------


def test_embl_has_no_date_source_or_keywords_key(tmp_path):
    """EMBL's scanner reads no `DT`, has no `SOURCE` consumer, and creates
    `keywords` only when a `KW` line was there -- a `DT` line stating a date
    does not put a `date` in the dict."""
    embl = seqio_corpus.build("embl", 1, 60, 0, 0)
    assert "DT   " not in embl  # the bare writer emits none
    path = _agree(tmp_path, "embl", embl, name="bare.embl")

    record = _reference(path, "embl")[0]
    assert "date" not in record.annotations
    assert "source" not in record.annotations
    assert "keywords" not in record.annotations

    # A `DT` line changes nothing: the key that would hold it does not exist.
    with_dt = embl.replace(
        "AC   ", "AC   SYN000001;\nDT   01-JAN-2000 (Rel. 62, Created)\n"
    )
    assert with_dt != embl
    path = _agree(tmp_path, "embl", with_dt, name="with_dt.embl")
    assert "date" not in _reference(path, "embl")[0].annotations


def test_genbank_placeholders_become_empty_values(tmp_path):
    """A GenBank writer's placeholders are values, not absences: `SOURCE   .`
    becomes `''`, `KW   .` becomes `['']`, and both keys exist."""
    genbank = seqio_corpus.build("genbank", 1, 60, 0, 0)
    assert "SOURCE      .\n" in genbank and "KEYWORDS    .\n" in genbank
    path = _agree(tmp_path, "genbank", genbank, name="bare.gb")

    annotations = _reference(path, "genbank")[0].annotations
    assert annotations["source"] == ""
    assert annotations["keywords"] == [""]


def test_a_record_with_no_reference_block_has_no_references_key(tmp_path):
    """"No references" is not `[]`.

    The reference's `reference_num` is what creates the key, so a record with no
    `REFERENCE`/`RN` line has *no* `references` key at all.  A reader that
    emitted an empty list would be right about the type and wrong about the
    dict, which is exactly the sort of difference `==` would catch and a
    value-by-value check would not.
    """
    genbank = seqio_corpus.build("genbank", 1, 60, 0, 0)
    assert "REFERENCE" not in genbank
    path = _agree(tmp_path, "genbank", genbank, name="no_ref.gb")

    record = _reference(path, "genbank")[0]
    assert "references" not in record.annotations
    index = biofasting.open_genbank(path, format="genbank")
    assert "references" not in _ours(index, record.id)


# --------------------------------------------------------------------------
# Header shapes the corpus writer never emits.
# --------------------------------------------------------------------------


def _spliced() -> list[tuple[str, str, str]]:
    """`(label, format, text)` for every hand-made header shape."""
    genbank, embl = _corpus_bases()
    cases: list[tuple[str, str, str]] = []

    def add(label: str, fmt: str, text: str) -> None:
        cases.append((label, fmt, text))

    def gb_at(anchor: str, replacement: str) -> str:
        assert genbank.count(anchor) == 1, anchor
        return genbank.replace(anchor, replacement)

    def embl_at(anchor: str, replacement: str) -> str:
        assert embl.count(anchor) == 1, anchor
        return embl.replace(anchor, replacement)

    references = (
        "REFERENCE   1  (bases 1 to 60)\n"
        "  AUTHORS   Someone.\n"
        "  TITLE     A title.\n"
        "  JOURNAL   J. 1, 1-2 (2000)\n"
        "  PUBMED    12345678\n"
    )

    add("genbank: bare header", "genbank", genbank)
    add(
        "genbank: VERSION with a GI number",
        "genbank",
        gb_at("VERSION     SYN000001\n", "VERSION     SYN000001.5  GI:1293613\n"),
    )
    add(
        "genbank: VERSION with no dot overwrites the ACCESSION id",
        "genbank",
        gb_at("VERSION     SYN000001\n", "VERSION     SYN000009\n"),
    )
    add(
        "genbank: COMMENT below the references",
        "genbank",
        gb_at("KEYWORDS    .\n", references + "COMMENT     a comment line\n"),
    )
    add(
        "genbank: COMMENT above the references",
        "genbank",
        gb_at("KEYWORDS    .\n", "COMMENT     a comment line\n" + references),
    )
    add(
        "genbank: two COMMENT blocks continue one another",
        "genbank",
        gb_at("KEYWORDS    .\n", "COMMENT     first\nCOMMENT     second\n"),
    )
    add(
        "genbank: a COMMENT wrapped onto a second line",
        "genbank",
        gb_at(
            "KEYWORDS    .\n",
            "COMMENT     first line\n            continued here\n",
        ),
    )
    add(
        "genbank: a taxonomy with several ranks",
        "genbank",
        gb_at(
            "            .\n",
            "            Bacteria; Proteobacteria; Escherichia.\n",
        ),
    )
    add(
        "genbank: several references",
        "genbank",
        gb_at(
            "KEYWORDS    .\n",
            references
            + references.replace("REFERENCE   1", "REFERENCE   2").replace(
                "PUBMED    12345678", "PUBMED    87654321"
            ),
        ),
    )

    add("embl: bare header", "embl", embl)
    add(
        "embl: the obsolete SV line",
        "embl",
        embl_at("DE   ", "SV   SYN000001.7\nDE   "),
    )
    add(
        "embl: a KW line, which is the only way the key exists",
        "embl",
        embl_at("DE   ", "DE   synthetic embl record 1 of 1.\nKW   one; two.\n"),
    )
    add(
        "embl: an OS/OC organism block",
        "embl",
        embl_at(
            "DE   ",
            "DE   synthetic embl record 1 of 1.\n"
            "OS   synthetic construct\n"
            "OC   Bacteria; Proteobacteria.\n",
        ),
    )
    add(
        "embl: a CC comment",
        "embl",
        embl_at("DE   ", "DE   synthetic embl record 1 of 1.\nCC   a comment.\n"),
    )
    add(
        "embl: a reference block",
        "embl",
        embl_at(
            "DE   ",
            "DE   synthetic embl record 1 of 1.\n"
            "RN   [1]\n"
            "RP   1-60\n"
            "RT   A title;\n"
            "RA   Someone;\n"
            "RL   J. 1, 1-2 (2000).\n",
        ),
    )
    add(
        "embl: an RP line with two ranges",
        "embl",
        embl_at(
            "DE   ",
            "DE   synthetic embl record 1 of 1.\n"
            "RN   [1]\n"
            "RP   16-55, 40-60\n"
            "RL   J. 1, 1-2 (2000).\n",
        ),
    )
    add(
        "embl: RX PUBMED and RX MEDLINE",
        "embl",
        embl_at(
            "DE   ",
            "DE   synthetic embl record 1 of 1.\n"
            "RN   [1]\n"
            "RP   1-60\n"
            "RX   PUBMED; 12345678.\n"
            "RX   MEDLINE; 99999999.\n"
            "RL   J. 1, 1-2 (2000).\n",
        ),
    )
    add(
        "embl: an RG consortium keeps its semicolon",
        "embl",
        embl_at(
            "DE   ",
            "DE   synthetic embl record 1 of 1.\n"
            "RN   [1]\n"
            "RP   1-60\n"
            "RG   The Consortium;\n"
            "RL   J. 1, 1-2 (2000).\n",
        ),
    )
    return cases


_SPLICED = _spliced()


@pytest.mark.parametrize(
    "label,format,text", _SPLICED, ids=[case[0] for case in _SPLICED]
)
def test_a_spliced_header_matches_the_reference(tmp_path, label, format, text):
    _agree(tmp_path, format, text, name="spliced.txt")


def test_the_comment_key_sits_where_its_line_sat(tmp_path):
    """The order assertion, made explicit.

    A `COMMENT` consumed before the first `REFERENCE` inserts `comment` before
    `references`; one consumed after inserts it after.  `dict.__eq__` cannot see
    the difference, so `_compare` asserts the order -- and this test says *why*
    the order is what it is, so that a reader that sorted its keys would fail
    here and not only inside a parametrised list.
    """
    by_label = {label: (format, text) for label, format, text in _SPLICED}
    _, above = by_label["genbank: COMMENT above the references"]
    _, below = by_label["genbank: COMMENT below the references"]

    for text, comment_first in ((above, True), (below, False)):
        path = tmp_path / "order.gb"
        path.write_text(text, newline="")
        keys = list(_reference(path, "genbank")[0].annotations)
        assert (keys.index("comment") < keys.index("references")) is comment_first

        index = biofasting.open_genbank(path, format="genbank")
        ours = list(_ours(index, _reference(path, "genbank")[0].id))
        assert (ours.index("comment") < ours.index("references")) is comment_first


# --------------------------------------------------------------------------
# The refusals: what this reader declines, and does not decline.
# --------------------------------------------------------------------------


def test_a_structured_comment_is_refused_and_names_itself(tmp_path):
    """The structured-comment probe is the cost M22 was built for; the block it
    detects is a second annotation key this reader does not produce, so the
    record refuses rather than return a dict missing it.

    The reference reads it (as `structured_comment`), so the test asserts that
    too -- otherwise a later change that made this shape reproducible would leave
    a refusal asserted against nothing.
    """
    genbank, _ = _corpus_bases()
    text = genbank.replace(
        "KEYWORDS    .\n",
        "COMMENT     ##Assembly-Data-START##\n"
        "            Assembly Method :: something\n"
        "            ##Assembly-Data-END##\n",
    )
    assert text != genbank
    path = tmp_path / "structured.gb"
    path.write_text(text, newline="")
    record = _reference(path, "genbank")[0]
    assert "structured_comment" in record.annotations

    index = biofasting.open_genbank(path, format="genbank")
    with pytest.raises(ValueError) as caught:
        _ours(index, record.id)
    assert "structured COMMENT" in str(caught.value)
    assert record.id in str(caught.value)


@pytest.mark.parametrize("keyword", ["NID", "PID", "DBSOURCE", "SEGMENT"])
def test_a_keyword_whose_annotation_is_not_reproduced_is_refused(tmp_path, keyword):
    """`NID`, `PID`, `DBSOURCE` and `SEGMENT` each set a key of the annotations
    dict that this reader does not produce, so a record carrying one refuses and
    the message names it.

    The keyword is padded to the column a GenBank keyword occupies, because both
    readers take a keyword to be the first twelve columns and not the first word:
    an unpadded `NID       something` is a keyword of `NID       so` to the
    reference too, and would be ignored by both.
    """
    genbank, _ = _corpus_bases()
    text = genbank.replace(
        "KEYWORDS    .\n", f"{keyword:<12}something\nKEYWORDS    .\n"
    )
    assert text != genbank
    path = tmp_path / "keyword.gb"
    path.write_text(text, newline="")
    record = _reference(path, "genbank")[0]

    index = biofasting.open_genbank(path, format="genbank")
    with pytest.raises(ValueError) as caught:
        _ours(index, record.id)
    assert keyword in str(caught.value)


def test_a_swissprot_record_refuses_its_annotations(tmp_path):
    """SwissProt shares the flat-file scan and not the header: its annotations
    are `Bio.SwissProt`'s, a different reader with different keys, and out of
    scope by the header's own statement.  Out of scope is a refusal, not an
    empty dict -- and the refusal is the one the scan's own first line raises,
    because SwissProt's `ID` line is not an EMBL `ID` line either."""
    path = tmp_path / "swiss.dat"
    path.write_text(seqio_corpus.build("swiss", 1, 60, 0, 0), newline="")
    record = _reference(path, "swiss")[0]

    index = biofasting.open_genbank(path, format="swiss")
    with pytest.raises(ValueError) as caught:
        _ours(index, record.id)
    assert "this reader reproduces" in str(caught.value)
    assert record.id in str(caught.value)


def test_the_refusal_does_not_take_the_sequence_with_it(tmp_path):
    """A refusal is per record and not per file.

    The features table makes the same promise -- `features()` refuses the record
    it cannot read and leaves the others usable -- and the annotations reader has
    to keep it, because a caller that fell back to `SeqIO.parse` for one record
    should not have to re-read the file for the rest.
    """
    genbank, _ = _corpus_bases()
    text = genbank.replace(
        "KEYWORDS    .\n", "SEGMENT     1 of 2\nKEYWORDS    .\n"
    )
    path = tmp_path / "mixed.gb"
    path.write_text(text + genbank.replace("SYN000001", "SYN000002"), newline="")

    records = _reference(path, "genbank")
    assert [record.id for record in records] == ["SYN000001", "SYN000002"]

    index = biofasting.open_genbank(path, format="genbank")
    with pytest.raises(ValueError):
        _ours(index, "SYN000001")
    _compare(
        _ours(index, "SYN000002"),
        records[1].annotations,
        "SYN000002",
    )
