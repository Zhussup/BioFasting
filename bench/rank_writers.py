#!/usr/bin/env python3
"""Ranking pass 7: the writers -- the half of `SeqIO` the triage calls `parity`.

Every format left in `Bio.SeqIO` was filed by the triage as `parity` -- "worth
having for drop-in compatibility, no large win expected".  Pass 6
(`rank_seqio.py`) re-asked that question for the *readers* with a timer and the
answer was no: the reference is 4.5x-11x above the floor C cannot cross.  The
one thing it left untouched is the other direction, and a prefill is not a
decision any more for `write` than it was for `parse`.

This script asks it.  A writer's cost sits in three places, and only the first
two are things a kernel can take:

* the **formatting** -- turning a record's fields into the text of the format;
* the **handle** -- pushing that text out, one small piece at a time or once;
* the **logic** -- what the reference does around those, which is Python.

So each row is measured against two floors:

* ``rewrite`` -- the identical bytes produced by the smallest pure-Python
  writer that can produce them, with the per-character loops replaced by
  C-level operations where one exists: ``str.translate`` for the FASTQ quality
  string, a regex for the FASTA wrap, one ``write`` for the record instead of
  one per line.  It is verified byte-for-byte against the reference before any
  time is quoted, so it is a floor under the same output and not a different
  program.
* ``payload`` -- the raw bytes of the sequences and qualities moved with
  ``b"".join``: what a kernel must at least pay, with no formatting above it.

The flat-file writers have no hand-written twin here, because a byte-identical
GenBank record is a page of column rules and reimplementing them would be the
milestone, not the measurement.  They get a per-method attribution instead --
the reference's own ``_write_the_first_lines`` and ``_write_sequence`` timed
separately against a handle that discards -- plus the number of ``write`` calls
per record, which for GenBank is the whole story: the residue block goes out ten
bases at a time.

    python3 bench/rank_writers.py
    python3 bench/rank_writers.py --only genbank-1kb --repeat 9
    python3 bench/rank_writers.py --list
"""

from __future__ import annotations

import argparse
import io
import re
import statistics
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import seqio_corpus  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.SeqIO.InsdcIO import EmblWriter, GenBankWriter  # noqa: E402
from Bio.SeqIO.Interfaces import _clean, _get_seq_string  # noqa: E402
from Bio.SeqIO.QualityIO import SANGER_SCORE_OFFSET  # noqa: E402

REPEATS = 7

#: The tables are indexed by row name; `--only` takes these.
ROWS = ("fasta-1kb", "fastq-150bp", "qual-150bp", "genbank-1kb", "embl-1kb")

#: FASTQ corpus, the file the Phase 1 FASTQ benchmark measures.
FASTQ_CORPUS = HERE / "data" / "fastq" / "reads_10k.fastq"


def timed(fn, repeats: int = REPEATS) -> float:
    """Median seconds over ``repeats`` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


# --- handles ----------------------------------------------------------------


class CountingHandle(io.StringIO):
    """A StringIO that counts its ``write`` calls.

    The count is the interesting half of the FASTA and GenBank rows: two
    implementations can write the same bytes and differ by two orders of
    magnitude in how many times they cross into the handle.
    """

    def __init__(self) -> None:
        super().__init__()
        self.writes = 0

    def write(self, text: str) -> int:
        self.writes += 1
        return super().write(text)


class NullHandle:
    """Discards everything, so a method can be timed with no text stored."""

    def write(self, text: str) -> int:
        return len(text)

    def flush(self) -> None:
        pass


# --- the record sets --------------------------------------------------------


def fastq_records() -> list:
    with FASTQ_CORPUS.open() as handle:
        return list(SeqIO.parse(handle, "fastq"))


def flat_records(group: str) -> tuple[str, list]:
    fmt, text = seqio_corpus.build_row(group)
    return fmt, list(SeqIO.parse(io.StringIO(text), fmt))


# --- the byte-identical rewrites --------------------------------------------


def _title(record) -> str:
    """The reference's own title rule, imported rather than restated.

    ``SeqIO`` decides the title from ``id`` and ``description`` with a rule that
    is easy to get subtly wrong -- a description that begins with the id is used
    whole, otherwise the two are joined.  ``_clean`` is imported from the
    reference for the same reason `rank_seqio` imports `minimal_parse`: two
    spellings of one rule is two places for it to drift.
    """
    id_ = _clean(record.id) if record.id else ""
    description = _clean(record.description)
    if description and description.split(None, 1)[0] == id_:
        return description
    if description:
        return f"{id_} {description}"
    return id_


#: The QUAL line rule, as one greedy match per line.
_QUAL_CHUNK = re.compile(r".{1,59}(?=\s|$)")


def _wrap_60(data: str) -> str:
    """The FASTA 60-column wrap, without the reference's per-line writes.

    The reference writes each 60-base slice to the handle as its own call --
    eighteen of them for a 1 kb record.  This builds the same lines with
    ``str.join`` and hands the record over in one piece.  A slice loop, not a
    regex: ``re.sub`` measured **slower** than the reference here (0.61 us
    against 0.27 on a 150 bp read), because a substitution has to build a
    replacement for every match, so the floor would have been a floor above the
    thing it was supposed to bound.
    """
    if not data:
        return "\n"
    lines = [data[i : i + 60] for i in range(0, len(data), 60)]
    return "\n".join(lines) + "\n"


def fasta_rewrite(records: list, handle) -> None:
    """FASTA, byte-identical, one ``write`` for the whole file."""
    parts = []
    for record in records:
        data = _get_seq_string(record)
        parts.append(f">{_title(record)}\n")
        parts.append(_wrap_60(data))
    handle.write("".join(parts))


#: PHRED -> Sanger ASCII, as a 256-byte translation table.  The reference walks
#: the quality list one element at a time with a dict lookup per base; this is
#: the same mapping as one C call, which is exactly what a kernel would do.
_PHRED_TO_SANGER = bytes(
    min(126, q + SANGER_SCORE_OFFSET) for q in range(256)
)


def _sanger_quality(record) -> str:
    return bytes(record.letter_annotations["phred_quality"]).translate(
        _PHRED_TO_SANGER
    ).decode("ascii")


def fastq_rewrite(records: list, handle) -> None:
    """FASTQ (Sanger), byte-identical, one ``write`` for the whole file."""
    parts = []
    for record in records:
        parts.append(
            f"@{_title(record)}\n"
            f"{_get_seq_string(record)}\n"
            f"+\n"
            f"{_sanger_quality(record)}\n"
        )
    handle.write("".join(parts))


def _qual_lines(qualities: list) -> str:
    """The QUAL wrap, at regex speed and without the reference's ``pop(0)``.

    The reference formats every quality with ``"%i" % round(q, 0)`` and then
    pops from the front of a list to pack the lines, which is O(n) a pop and so
    O(n^2) a record.  For an integer quality ``str(q)`` is the same characters
    as ``"%i" % round(q, 0)``, so the join is one C-level call over the whole
    record; the wrap is then a greedy match ending on a token boundary.

    Equality with the reference's rule is exact rather than approximate: the
    reference appends a token while ``len(line) + 1 + len(token) < 60``, so a
    line never exceeds 59 characters, and ``.{1,59}(?=\\s|$)`` is the longest
    prefix of at most 59 characters that ends at a token boundary.  The gate
    checks it byte for byte rather than trusting that argument.
    """
    text = " ".join(map(str, qualities))
    out = []
    pos = 0
    count = len(text)
    while pos < count:
        match = _QUAL_CHUNK.match(text, pos)
        out.append(match.group(0))
        out.append("\n")
        pos = match.end()
        if pos < count:
            pos += 1  # step over the separator
    return "".join(out)


def qual_rewrite(records: list, handle) -> None:
    """QUAL, byte-identical, one ``write`` for the whole file."""
    parts = []
    for record in records:
        parts.append(f">{_title(record)}\n")
        parts.append(_qual_lines(record.letter_annotations["phred_quality"]))
    handle.write("".join(parts))


def payload(records: list, fmt: str) -> bytes:
    """The bytes each record contributes with no formatting at all.

    Sequences for every row, plus the qualities for the two that carry them.
    Built from the records inside the timed region rather than handed over
    already-joined -- ``bytes(existing_bytes)`` returns the same object and
    would have timed a no-op, which is how the first run of this script
    reported a payload of 0.00 us a record.
    """
    chunks = []
    for record in records:
        chunks.append(_get_seq_string(record).encode("ascii"))
        if fmt in ("fastq", "qual"):
            chunks.append(bytes(record.letter_annotations["phred_quality"]))
    return b"".join(chunks)


# --- what each row measures -------------------------------------------------

REWRITES = {
    "fasta-1kb": ("fasta", fasta_rewrite),
    "fastq-150bp": ("fastq", fastq_rewrite),
    "qual-150bp": ("qual", qual_rewrite),
}


def reference_write(fmt: str, records: list, handle) -> None:
    """`SeqIO.write`, the thing a caller would otherwise use."""
    SeqIO.write(records, handle, fmt)


def row_records(name: str) -> tuple[str, list]:
    """The ``(format, records)`` behind a row name."""
    if name == "fastq-150bp" or name == "qual-150bp":
        return name.split("-")[0], fastq_records()
    if name == "fasta-1kb":
        return "fasta", flat_records("genbank-1kb")[1]
    return flat_records(name)


def gate(fmt: str, records: list) -> tuple[bool, str]:
    """Byte-identity of the rewrite against the reference, before any timing."""
    expected = io.StringIO()
    reference_write(fmt, records, expected)
    got = io.StringIO()
    REWRITES[fmt_key(fmt)][1](records, got)
    if got.getvalue() != expected.getvalue():
        return False, "the rewrite's bytes differ from the reference's"
    return True, ""


def fmt_key(fmt: str) -> str:
    """The `REWRITES` key a format name belongs to."""
    for key, (name, _) in REWRITES.items():
        if name == fmt:
            return key
    raise KeyError(fmt)


def measured_rows(only: str | None, repeats: int) -> None:
    """The three writers with a byte-identical twin, and their floors."""
    header = (
        f"{'row':14s} {'format':8s} {'recs':>6s} {'ref us/rec':>11s} "
        f"{'rewrite':>8s} {'payload':>8s} {'writes/rec':>10s} {'ratio':>6s}"
    )
    print()
    print("Writers with a byte-identical pure-Python twin.  `rewrite` and")
    print("`payload` are the same output through C-level operations and the raw")
    print("bytes with no formatting; the gate is over the whole output first.")
    print(header)
    print("-" * len(header))
    for name in ("fasta-1kb", "fastq-150bp", "qual-150bp"):
        if only and only != name:
            continue
        fmt, records = row_records(name)
        ok, why = gate(fmt, records)
        if not ok:
            print(f"{name:14s} {fmt:8s} INVALID -- {why}")
            continue
        count = len(records)

        def ref() -> None:
            reference_write(fmt, records, io.StringIO())

        ref_seconds = timed(ref, repeats) / count

        rewrite_fn = REWRITES[fmt_key(fmt)][1]

        def rewr() -> None:
            rewrite_fn(records, io.StringIO())

        rewrite_seconds = timed(rewr, repeats) / count

        def pay() -> None:
            payload(records, fmt)

        payload_seconds = timed(pay, repeats) / count

        counting = CountingHandle()
        ref_batched = lambda: reference_write(fmt, records, counting)  # noqa: E731
        counting.writes = 0
        ref_batched()
        writes_per_record = counting.writes / count

        print(
            f"{name:14s} {fmt:8s} {count:6d} {ref_seconds * 1e6:11.2f} "
            f"{rewrite_seconds * 1e6:8.2f} {payload_seconds * 1e6:8.2f} "
            f"{writes_per_record:10.1f} {ref_seconds / rewrite_seconds:5.2f}x"
        )


# --- the flat-file writers --------------------------------------------------

#: The methods `write_record` reaches for, per format.  GenBank has no
#: `_write_the_first_lines`: its DEFINITION-through-SOURCE block is written
#: inline by `write_record` with `_write_multi_line`/`_write_single_line`, so
#: the table carries a `rest` column for the part no single method owns --
#: stated as a remainder rather than attributed to a name it does not have.
INSNC_PARTS = {
    "genbank": (
        "_write_the_first_line",
        "_write_references",
        "_write_comment",
        "_write_sequence",
    ),
    "embl": (
        "_write_the_first_lines",
        "_write_keywords",
        "_write_references",
        "_write_comment",
        "_write_sequence",
    ),
}


def _attribute(fmt: str, group: str, repeats: int) -> dict:
    """One corpus row, priced method by method against a discarding handle."""
    _, records = row_records(group)
    count = len(records)
    writer_class = GenBankWriter if fmt == "genbank" else EmblWriter

    def whole() -> None:
        writer = writer_class(NullHandle())
        for record in records:
            writer.write_record(record)

    out = {
        "recs": count,
        "write_record": timed(whole, repeats) / count,
        "writes/rec": None,
    }

    for method_name in INSNC_PARTS[fmt]:
        method = getattr(writer_class, method_name, None)
        if method is None:
            out[method_name] = None
            continue

        def one(method=method) -> None:
            writer = writer_class(NullHandle())
            for record in records:
                method(writer, record)

        #   A method can be present and still not apply to this corpus:
        #   `_write_references` raises `KeyError: 'references'` on a record with
        #   no reference block, and `_write_comment` raises `IndexError` on an
        #   empty one, which the bare rows are.  That is an answer -- the method
        #   costs nothing here -- and not a crash, so it is probed once outside
        #   the timer and reported as `-`.
        try:
            one()
        except (KeyError, IndexError, AttributeError):
            out[method_name] = None
            continue

        out[method_name] = timed(one, repeats) / count

    listed = sum(v for k, v in out.items()
                 if k in INSNC_PARTS[fmt] and v is not None)
    out["rest"] = out["write_record"] - listed

    counting = CountingHandle()
    counting.writes = 0
    writer = writer_class(counting)
    for record in records:
        writer.write_record(record)
    out["writes/rec"] = counting.writes / count
    return out


def flat_rows(only: str | None, repeats: int) -> None:
    #   The annotated rows are here as well as the bare ones because a method
    #   that has nothing to write costs nothing, and a table over the bare rows
    #   alone would say the reference's header is cheap.  `_write_references`
    #   and `_write_comment` raise `KeyError`/`IndexError` on a record carrying
    #   neither, which is `-` here and the whole point of showing both.
    #
    #   Each format gets its own column headings because the two do not share
    #   their methods -- GenBank has `_write_the_first_line` and no
    #   `_write_the_first_lines`, EMBL the reverse -- and a shared heading is
    #   how a table ends up labelling a keyword count as a reference count.
    print()
    print("The flat-file writers, attributed by the reference's own methods")
    print("against a handle that discards.  `sequence` is `_write_sequence`, the")
    print("ORIGIN/SQ residue block, and `rest` is everything `write_record` writes")
    print("that no listed method owns -- the DEFINITION-to-SOURCE block for")
    print("GenBank, the keywords and the ID line for EMBL.  A method with nothing")
    print("to write on a corpus row is `-`, which is why the bare and the")
    print("annotated row are both here.  Times are us per record.")
    for fmt, groups in (("genbank", ("genbank-1kb", "genbank-1kb-annot")),
                        ("embl", ("embl-1kb", "embl-1kb-annot"))):
        if only and only not in groups:
            continue
        tables = [(group, _attribute(fmt, group, repeats)) for group in groups]
        labels = ("recs", "write_record", *INSNC_PARTS[fmt], "rest", "writes/rec")
        width = 22
        cell = 18
        print()
        print(fmt)
        print(f"  {'':{width}s}" + "".join(f"{g:>{cell}s}" for g, _ in tables))
        for label in labels:
            cells = []
            for _, table in tables:
                value = table[label]
                if value is None:
                    cells.append(f"{'-':>{cell}s}")
                elif label == "recs":
                    cells.append(f"{value:>{cell}d}")
                elif label == "writes/rec":
                    cells.append(f"{value:>{cell}.1f}")
                else:
                    cells.append(f"{value * 1e6:>{cell}.2f}")
            print(f"  {label:{width}s}" + "".join(cells))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", help="one row name")
    parser.add_argument("--repeat", type=int, default=REPEATS)
    parser.add_argument("--list", action="store_true", help="list the rows")
    args = parser.parse_args()

    if args.list:
        for name in ROWS:
            print(name)
        return 0

    measured_rows(args.only, args.repeat)
    flat_rows(args.only, args.repeat)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
