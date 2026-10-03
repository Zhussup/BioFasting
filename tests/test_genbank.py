"""The flat-file reader against the reference, record for record.

`biofasting.open_genbank` reproduces four fields of a `SeqRecord`: `id`, `name`,
`description` and the sequence.  Anything short of "identical to
`Bio.SeqIO.parse` on every record of every input" is a bug, so every test here
compares against the reference rather than against a literal this file would
otherwise have to keep in step by hand.

Two families of input:

- the generated corpus (`bench/seqio_corpus.py`), which is what the sixth
  ranking pass measured and therefore what the kernel was built for;
- records for the shapes the corpus's own writer never emits -- a definition
  continued onto a second line, a `VERSION` carrying a suffix, an EMBL `SV`
  field, lowercase residues, CRLF endings -- because a writer emits its own
  habits, and a parser that has only ever met one writer's output has been
  tested against that writer and not against the format.

Those records start from the reference's own writer and are spliced, rather than
being written out here from the column numbers: a hand-written `LOCUS` line is
how a suite ends up testing its author's arithmetic instead of the format, and
the corpus module had to learn that the same way.

The refusals are asserted too.  The reader is deliberately stricter than
Biopython where Biopython warns and guesses, and a stricter reader that stopped
being stricter without anyone noticing would be a silent difference in a
sequence.
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import sys
from pathlib import Path

import pytest

#   The module is a differential comparison against Biopython from end to end,
#   so with the reference absent there is nothing in it to run.  The guard is an
#   `importorskip` above the `Bio` imports rather than a `needs_biopython`
#   marker, because a module-level import that raises is a collection error --
#   the run stops before any marker is reached.  CI builds the wheel with no
#   reference installed, and a skip is the answer there; an error is not.
pytest.importorskip(
    "Bio", reason="Biopython is not installed; it is the differential reference"
)

from Bio import SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

import biofasting  # noqa: E402

BENCH = Path(__file__).resolve().parent.parent / "bench"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"_bench_{name}", BENCH / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


seqio_corpus = _load("seqio_corpus")

CORPUS_ROWS = [row[0] for row in seqio_corpus.CORPORA]
# The corpus is five megabytes of GenBank; building it once per test would make
# this file the slowest in the suite for no reason, and `corpus_rows()` builds
# every row -- so a row is built here by asking the generator for that row only.
_ROW_SPEC = {row[0]: row for row in seqio_corpus.CORPORA}
_CACHE: dict[str, tuple[str, str]] = {}


def _corpus(group: str) -> tuple[str, str]:
    """The (format, text) of one corpus row, built once and only that row."""
    if group not in _CACHE:
        _, format, count, length, per_kb, annot = _ROW_SPEC[group]
        _CACHE[group] = (format,
                         seqio_corpus.build(format, count, length, per_kb, annot))
    return _CACHE[group]


def _reference(path: Path, format: str) -> list[tuple[str, str, str, str]]:
    """The four fields under test, taken from the reference's own records.

    Read from the same file the kernel maps rather than from the text it was
    written from, so that a line ending or an encoding the write and the read do
    not agree about shows up as a mismatch instead of being hidden by both sides
    going through the same in-memory string.
    """
    with path.open() as handle:
        return [
            (record.id, record.name, record.description, str(record.seq))
            for record in SeqIO.parse(handle, format)
        ]


def _ours(path: Path, format: str) -> list[tuple[str, str, str, str]]:
    return [
        (record.id, record.name, record.description, record.sequence.decode())
        for record in biofasting.read_genbank(path, format=format)
    ]


def _write(tmp_path: Path, text: str, name: str) -> Path:
    path = tmp_path / name
    path.write_text(text, newline="")
    return path


def _agree(tmp_path: Path, format: str, text: str, name: str = "in.txt") -> Path:
    """Write `text`, then hold the kernel to the reference on it."""
    path = _write(tmp_path, text, name)
    assert _ours(path, format) == _reference(path, format)
    return path


# --------------------------------------------------------------------------
# The generated corpus: what the kernel was measured against.
# --------------------------------------------------------------------------


@pytest.mark.parametrize("group", CORPUS_ROWS)
def test_corpus_row_matches_the_reference(tmp_path, group):
    format, text = _corpus(group)
    _agree(tmp_path, format, text, name=f"{group}.txt")


@pytest.mark.parametrize("group", CORPUS_ROWS)
def test_corpus_row_is_indexed_by_id(tmp_path, group):
    """Every id is addressable, in file order, and answers the same bytes."""
    format, text = _corpus(group)
    path = _write(tmp_path, text, f"{group}.txt")
    index = biofasting.open_genbank(path, format=format)

    expected = _reference(path, format)
    assert list(index) == [record[0] for record in expected]
    assert len(index) == len(expected)
    for record in expected:
        assert record[0] in index
        assert index[record[0]].decode() == record[3]
        assert index.name(record[0]) == record[1]
        assert index.description(record[0]) == record[2]
        assert index.sequence_length(record[0]) == len(record[3])


@pytest.mark.parametrize("group", CORPUS_ROWS)
def test_slices_are_the_sequence_the_reference_produced(tmp_path, group):
    """A window is a gather out of the same bytes, not a second parse."""
    format, text = _corpus(group)
    path = _write(tmp_path, text, f"{group}.txt")
    index = biofasting.open_genbank(path, format=format)

    for record in _reference(path, format)[:20]:
        sequence = record[3]
        assert index.sequence_slice(record[0], 0, len(sequence)) == sequence.encode()
        assert index.sequence_slice(record[0], 3, 47) == sequence[3:47].encode()
        # Clamped, not refused: a window past the end is the record, and an
        # inverted one is empty.
        assert index.sequence_slice(record[0], len(sequence) - 5, 10**9) == (
            sequence[-5:].encode()
        )
        assert index.sequence_slice(record[0], 47, 3) == b""


def test_the_scan_is_one_pass_over_the_whole_file(tmp_path):
    """Each record's span starts at its own first line and holds its own residues.

    The offsets are checked against the file's own bytes rather than against each
    other: records whose spans merely fail to overlap would satisfy a weaker
    test while the reader re-derived the second record from the first one's end.
    """
    format, text = _corpus("genbank-1kb")
    path = _write(tmp_path, text, "genbank.gb")
    data = path.read_bytes()
    index = biofasting.open_genbank(path, format=format)

    offsets = index.offsets()
    assert offsets == sorted(offsets)
    assert len(offsets) == len(index)

    reference = {record[0]: record[3] for record in _reference(path, format)}
    for record in biofasting.read_genbank(path, format=format):
        offset = index.location(record.id)[0]
        assert data[offset : offset + 5] == b"LOCUS"
        # The length counted during the scan is the length that comes out: if
        # these could disagree, a caller reserving on one and reading the other
        # would be reading a promise nobody kept.
        assert len(record.sequence) == index.sequence_length(record.id)
        assert record.sequence.decode() == reference[record.id]


# --------------------------------------------------------------------------
# The shapes the corpus's writer never emits.
# --------------------------------------------------------------------------

_RESIDUES = "gatcgtacgt" * 10
_ANNOTATIONS = {
    "molecule_type": "DNA",
    "data_file_division": "PLN",
    "date": "01-JAN-2026",
    "topology": "linear",
}


def _written(format: str, sequence: str = _RESIDUES, *, ident: str = "TEST001",
             description: str = "a test record") -> str:
    """One record as the reference's own writer emits it.

    The header block's columns are fixed and unforgiving, so the skeleton comes
    from the writer and the shapes under test are spliced into it.  This is the
    same call `bench/seqio_corpus.py` makes, and for the same reason.
    """
    record = SeqRecord(Seq(sequence), id=ident, name=ident, description=description)
    record.annotations = dict(_ANNOTATIONS)
    handle = io.StringIO()
    SeqIO.write([record], handle, format)
    return handle.getvalue()


def _splice(text: str, old: str, new: str) -> str:
    """Replace one line, asserting it was there -- a splice that missed is a test that tests nothing."""
    assert old in text, f"fixture does not contain {old!r}"
    return text.replace(old, new, 1)


def _swiss(sequence: str = "MTESTRESIDUE", *, ident: str = "TEST001",
           accession: str = "TEST001",
           definition: str = "RecName: Full=test protein;") -> str:
    """A SwissProt record, written out by hand because there is no writer.

    `Reading format 'swiss' is supported, but not writing` -- so this is the one
    fixture that cannot borrow the reference, and the reference's acceptance of
    it is the check that it is a record at all.  `definition` is the DE *values*,
    one per line; the keyword column is added here so that a caller cannot get
    it wrong.
    """
    lines = [
        f"ID   {ident:<16}Reviewed;         {len(sequence)} AA.",
        f"AC   {accession};",
        "DT   01-JAN-2026, integrated into UniProtKB/Swiss-Prot.",
        *(f"DE   {value}" for value in definition.splitlines()),
        "GN   Name=TEST001;",
        "OS   Synthetic construct.",
        "OC   artificial sequences; synthetic construct.",
        "OX   NCBI_TaxID=32630;",
        f"SQ   SEQUENCE   {len(sequence)} AA;  1000 MW;  ABCDEF0123456789 CRC64;",
    ]
    for start in range(0, len(sequence), 60):
        row = sequence[start : start + 60]
        lines.append("     " + " ".join(row[i : i + 10] for i in range(0, len(row), 10)))
    lines.append("//")
    return "\n".join(lines) + "\n"


def test_genbank_definition_continued_on_a_second_line(tmp_path):
    """A definition too long for one line, and its trailing stop dropped once."""
    text = _splice(
        _written("genbank"),
        "DEFINITION  a test record.\n",
        "DEFINITION  a test record whose\n            definition runs on.\n",
    )
    path = _agree(tmp_path, "genbank", text)
    assert biofasting.read_genbank(path)[0].description == (
        "a test record whose definition runs on"
    )


def test_genbank_version_suffix_becomes_part_of_the_id(tmp_path):
    """`VERSION     TEST001.2` is what turns the id into `TEST001.2`."""
    text = _splice(_written("genbank"), "VERSION     TEST001\n", "VERSION     TEST001.2\n")
    path = _agree(tmp_path, "genbank", text)
    assert biofasting.read_genbank(path)[0].id == "TEST001.2"


def test_genbank_without_an_accession_falls_back_to_the_locus_name(tmp_path):
    """An id is a promise a record can be addressed, so the name stands in."""
    text = _splice(
        _splice(_written("genbank"), "ACCESSION   TEST001\n", ""),
        "VERSION     TEST001\n",
        "",
    )
    path = _agree(tmp_path, "genbank", text)
    assert biofasting.read_genbank(path)[0].id == "TEST001"


def test_embl_sv_field_versions_the_id(tmp_path):
    """EMBL's `SV n` is a version suffix, and it reaches the id the same way."""
    text = _splice(_written("embl"), "ID   TEST001; ; linear;", "ID   TEST001; SV 2; linear;")
    path = _agree(tmp_path, "embl", text)
    assert biofasting.read_genbank(path, format="embl")[0].id == "TEST001.2"


def test_embl_takes_the_id_from_the_id_line_not_the_accession(tmp_path):
    """The two lines disagree here, and the ID line is the one the reference reads.

    A test that used corpus data could not tell these apart: Biopython's writer
    puts the same string on both lines, so a reader that took the id from `AC`
    would pass every corpus row while being wrong about every real EMBL file.
    """
    text = _splice(_written("embl"), "AC   TEST001;", "AC   ACCLINE;")
    path = _agree(tmp_path, "embl", text)
    record = biofasting.read_genbank(path, format="embl")[0]
    assert record.id == "TEST001"
    assert record.name == "TEST001"


def test_embl_residues_are_uppercased(tmp_path):
    """INSDC stores lowercase residues and the reference uppercases them."""
    path = _agree(tmp_path, "embl", _written("embl", _RESIDUES.lower()))
    record = biofasting.read_genbank(path, format="embl")[0]
    assert record.sequence == _RESIDUES.upper().encode()


def test_swiss_residues_keep_their_case(tmp_path):
    """SwissProt's reader does not uppercase, and a corpus of capitals cannot say so.

    This is the one place where the two INSDC formats and SwissProt disagree
    about the same bytes, so it is asserted rather than left to the corpus --
    whose hand-written emitter happens to write capitals everywhere.
    """
    path = _agree(tmp_path, "swiss", _swiss("mtestresidue"))
    record = biofasting.read_genbank(path, format="swiss")[0]
    assert record.sequence == b"mtestresidue"


def test_swiss_description_joins_every_de_line(tmp_path):
    path = _agree(
        tmp_path,
        "swiss",
        _swiss(definition="RecName: Full=test protein;\nAltName: Full=other;"),
    )
    assert biofasting.read_genbank(path, format="swiss")[0].description == (
        "RecName: Full=test protein; AltName: Full=other;"
    )


def test_swiss_is_indexed_by_its_accession(tmp_path):
    path = _agree(tmp_path, "swiss", _swiss(ident="ENTRYNAME", accession="P00321"))
    record = biofasting.read_genbank(path, format="swiss")[0]
    assert record.id == "P00321"
    assert record.name == "ENTRYNAME"


def test_carriage_returns_are_line_endings_and_not_residues(tmp_path):
    """A CRLF file is a text file, so the reference never sees the CR.

    Getting this wrong is quiet rather than loud: a `\\r` left in a residue line
    is a byte that is not a residue, so the sequence is one character longer per
    line than the record it came from says it is.
    """
    _agree(tmp_path, "genbank", _written("genbank").replace("\n", "\r\n"))


def test_several_records_in_one_file(tmp_path):
    """The record separator is what ends a record, not the end of a definition."""
    text = _written("genbank", ident="A1") + _written("genbank", ident="A2")
    path = _agree(tmp_path, "genbank", text)
    assert list(biofasting.open_genbank(path)) == ["A1", "A2"]


def test_blank_lines_between_records_are_not_records(tmp_path):
    text = _written("genbank", ident="A1") + "\n\n" + _written("genbank", ident="A2")
    path = _agree(tmp_path, "genbank", text)
    assert list(biofasting.open_genbank(path)) == ["A1", "A2"]


def test_a_duplicate_id_is_refused(tmp_path):
    """The same refusal `SeqIO.index` makes, and for the same reason."""
    path = _write(tmp_path, _written("genbank") + _written("genbank"), "in.gb")
    with pytest.raises(ValueError, match="Duplicate key 'TEST001'"):
        biofasting.open_genbank(path)


def test_a_shifted_origin_column_is_refused_rather_than_guessed(tmp_path):
    """The one deliberate disagreement with the reference.

    Biopython warns that the indentation is wrong, shifts the line by a byte and
    parses it anyway.  Reproducing a guess is how a sequence that looks right
    stops being right, so this reader refuses and names the byte it refused at.
    """
    text = _splice(_written("genbank"), "\n        1 gatcgtacgt", "\n1 gatcgtacgt")
    path = _write(tmp_path, text, "in.gb")
    with pytest.raises(ValueError, match="declines rather than guess"):
        biofasting.open_genbank(path)


def test_a_contig_line_inside_the_sequence_block_is_refused(tmp_path):
    """`CONTIG` says the residues are elsewhere, and this reader has no elsewhere.

    A record that says so has no sequence to reproduce, and returning an empty
    one would be a record that looks parsed and is not.
    """
    text = _splice(
        _written("genbank"), "\nORIGIN\n", "\nORIGIN\nCONTIG      join(1..100)\n"
    )
    path = _write(tmp_path, text, "in.gb")
    with pytest.raises(ValueError, match="declines rather than guess"):
        biofasting.open_genbank(path)


def test_an_unknown_format_is_refused(tmp_path):
    path = _write(tmp_path, _written("genbank"), "in.gb")
    with pytest.raises(ValueError, match="unknown flat-file format"):
        biofasting.open_genbank(path, format="fasta")


def test_a_compressed_file_is_refused_rather_than_indexed_as_text(tmp_path):
    path = tmp_path / "in.gb.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(_written("genbank"))
    with pytest.raises(ValueError, match="gzip-compressed"):
        biofasting.open_genbank(path)


def test_a_missing_id_is_a_key_error(tmp_path):
    path = _write(tmp_path, _written("genbank"), "in.gb")
    index = biofasting.open_genbank(path)
    with pytest.raises(KeyError):
        index["NOT_THERE"]
    assert "NOT_THERE" not in index


def test_records_and_the_index_agree(tmp_path):
    """Two ways into the same file, one answer."""
    path = _write(tmp_path, _written("genbank"), "in.gb")
    index = biofasting.open_genbank(path)
    assert [tuple(record) for record in index.records()] == [
        (record.id, record.name, record.description, record.sequence)
        for record in biofasting.read_genbank(path)
    ]
