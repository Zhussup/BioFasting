#!/usr/bin/env python3
"""Ranking pass 6: the remaining `SeqIO` file formats (PLAN 2.4).

The triage files every format left in `Bio.SeqIO` as `parity` -- "worth having
for drop-in compatibility, no large win expected" -- and the `perf` rows there
are all closed: `FastaIO`, `QualityIO` and `_index`, 62 of them, delivered in
Phase 1.  A prefill is not a decision, so this pass asks the triage's question
again with a timer: what does the reference cost per record, and how much of
that cost is work a kernel could take?

The second question is the one that decides.  A parse spends its time in three
places -- moving bytes, building Python objects, and its own logic -- and only
the first two can be made faster by C.  So each format is measured against two
floors:

* ``objects`` -- build the same ``SeqRecord`` tree from scalars already in hand,
  locations and qualifiers included.  This is the Python work a parser cannot
  avoid, and no kernel can go below it.
* ``scan`` -- a whole-file ``bytes.translate``/``upper``: what moving the bytes
  costs at C speed, with no parsing above it.

Anything the reference spends above both floors is its own line-by-line logic,
and that is the part a kernel removes.  Usage:

    python3 bench/rank_seqio.py
"""

from __future__ import annotations

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
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import SeqFeature, SimpleLocation  # noqa: E402
from Bio.SeqRecord import SeqRecord  # noqa: E402

REPEATS = 7


def timed(fn, repeats: int = REPEATS) -> float:
    """Median seconds over ``repeats`` calls, after one warm-up."""
    fn()
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def scalars_of(records: list) -> list:
    """Everything a parser has to have, as plain values rather than objects."""
    out = []
    for rec in records:
        features = [
            (tuple((p.start, p.end, p.strand) for p in f.location.parts),
             f.type, dict(f.qualifiers))
            for f in rec.features
        ]
        out.append((rec.id, rec.name, rec.description, str(rec.seq),
                    dict(rec.annotations), features))
    return out


def build_objects(scalars: list) -> list:
    """The floor: the same tree, built from scalars, with no text touched."""
    out = []
    for rid, name, desc, seq, annotations, features in scalars:
        rec = SeqRecord(Seq(seq), id=rid, name=name, description=desc)
        rec.annotations = annotations
        built = []
        for parts, ftype, qualifiers in features:
            first, rest = parts[0], parts[1:]
            location = SimpleLocation(first[0], first[1], strand=first[2])
            for part in rest:
                location = location + SimpleLocation(part[0], part[1],
                                                     strand=part[2])
            built.append(SeqFeature(location, type=ftype,
                                    qualifiers=qualifiers))
        rec.features = built
        out.append(rec)
    return out


def minimal_parse(fmt: str, text: str) -> list[tuple[str, str]]:
    """The best pure-Python rewrite there is: ids and sequences, nothing else.

    This does **not** reproduce the reference's records.  It drops the features,
    the description and all nine annotations, so its ratio is a floor under what
    a Python rewrite could cost and never a drop-in replacement.  It is measured
    anyway because it answers the question the ranking pass exists to ask: is
    the reference's time going into the bytes, into the objects, or into its own
    logic?  A rewrite this small being several times faster than the reference is
    the proof that the answer is the third one.

    Three shapes, not two: GenBank keeps its residues behind an ``ORIGIN``
    header, while EMBL and SwissProt put them after an ``SQ`` line.  Writing this
    with one branch for GenBank and one for "the rest" is what made the first
    draft of it raise ``IndexError`` on every EMBL row.
    """
    out = []
    for chunk in text.split("//\n"):
        if not chunk.strip():
            continue
        lines = chunk.split("\n")
        if fmt == "genbank":
            tag = lines[0].split()[1]
            body = chunk.split("\nORIGIN\n", 1)[1]
            out.append((tag, body.translate(_SCAN_TABLE).replace("\n", "").upper()))
            continue
        if fmt == "embl":
            tag = lines[0].split()[1].rstrip(";")
        else:
            tag = next(line[5:].split(";")[0] for line in lines
                       if line.startswith("AC   "))
        body = chunk.split("\nSQ   ", 1)[1].split("\n", 1)[1]
        #   EMBL right-aligns the residue count on the last line of every `SQ`
        #   block, so the digits have to go with the spaces.  Leaving them in
        #   appends the record's own length to its sequence, and the comparison
        #   column below is what caught it.  Stripping digits is harmless for
        #   SwissProt, whose alphabet has none.
        out.append((tag, body.translate(_SCAN_TABLE).upper()))
    return out


def reference_pairs(fmt: str, text: str) -> list[tuple[str, str]]:
    """The same two fields, taken from the reference's own records."""
    return [(record.id, str(record.seq))
            for record in SeqIO.parse(io.StringIO(text), fmt)]


def annotation_variant(fmt: str, count: int, length: int, blocks: tuple) -> str:
    """One corpus-shaped file with only ``blocks`` of the annotation set.

    The annotation rows in `seqio_corpus` differ from their bare twin by three
    blocks at once, and a single delta would say "the annotations dict costs
    30 us" without saying which part of it does.  This builds the same records
    with a chosen subset, so the cost can be attributed: a subset is a property
    of the file, and the reference is timed on the file.

    Every variant draws its residues from the same seed and the same domain
    string, because the blocks are applied after the sequence draw -- so the
    difference between any two variants here is the blocks and nothing else,
    and a one-block variant prices that block against the bare row.
    """
    from Bio.SeqFeature import Reference

    records = seqio_corpus._records(fmt, count, length, 0)
    for record in records:
        if "taxonomy" in blocks:
            record.annotations["taxonomy"] = ["artificial sequences",
                                              "synthetic construct"]
            record.annotations["keywords"] = ["synthetic", "benchmark"]
        if "comment" in blocks:
            record.annotations["comment"] = (
                "This record is generated for benchmarking and carries a "
                "comment long enough to span several lines of the flat file, "
                "which the reference reads back with the line breaks it wrote."
            )
        if "references" in blocks:
            references = []
            for j in range(3):
                reference = Reference()
                reference.title = f"Synthetic study {j + 1}"
                reference.authors = f"Author {j + 1}, B."
                reference.journal = "J. Synthetic Biol."
                reference.pubmed_id = str(1_000_000 + j)
                reference.comment = "primary"
                references.append(reference)
            record.annotations["references"] = references
    handle = io.StringIO()
    SeqIO.write(records, handle, fmt)
    return handle.getvalue()


def annotation_costs() -> None:
    """The annotations dict, which the first six passes left as one lump.

    The sixth pass said the reference spends 20.2 us per feature and left the
    rest of the header alone.  The annotation rows answer the other half: the
    same sequences, the same seed, the same domain string, with the header as
    the only difference.  This prints the delta and then attributes it, because
    "the annotations dict costs 31 us" is a target a kernel cannot be held to
    without knowing which block the time is in.

    The three blocks are measured **one at a time** against the bare row rather
    than cumulatively, so that a column is a cost and not a cost-minus-the-last-
    column.  They are close to additive, and the table prints the sum next to the
    pair's own delta rather than asserting they agree -- a block that cost 20 us
    alone and 2 us in company would be a fact worth seeing.

    The comment column is the one that does not mean what it says, which is why
    the second table is here: it is not the price of reading a comment.  The
    reference's GenBank path runs ``re.search(r"([^#]+)-START##$", line)`` on the
    first line of every ``COMMENT`` block to decide whether the record carries a
    *structured* comment, and the pattern backtracks over that line one character
    at a time when it fails.  A record with a 68-character first comment line
    pays for 68 failed attempts at a suffix that is not there, and pays it
    whether or not it has a structured comment -- which most records do not.
    """
    print()
    print("The annotations dict, isolated by the corpus's annotation rows: the")
    print("same residues from the same seed, with the header as the only")
    print("difference, and each block then priced on its own against the bare")
    print("row rather than cumulatively.")
    header = (f"{'corpus':18s} {'bare us/rec':>11s} {'annot':>9s} {'pair delta':>10s}  "
              f"{'+tax/kw':>8s} {'+comment':>9s} {'+refs':>7s} {'sum':>8s}")
    print(header)
    print("-" * len(header))
    for fmt, bare_group, annot_group in (
        ("genbank", "genbank-1kb-bare", "genbank-1kb-annot"),
        ("embl", "embl-1kb-bare", "embl-1kb-annot"),
    ):
        _, bare = seqio_corpus.build_row(bare_group)
        _, annotated = seqio_corpus.build_row(annot_group)
        count = seqio_corpus.spec(annot_group)[2]
        length = seqio_corpus.spec(annot_group)[3]

        def ref(text):
            return timed(lambda: list(SeqIO.parse(io.StringIO(text), fmt)))

        #   Built before the timing, never inside it: building the file is the
        #   writer's cost and the writer is not what is being measured here.
        blocks = {name: ref(annotation_variant(fmt, count, length, (name,)))
                  for name in ("taxonomy", "comment", "references")}
        bare_us = ref(bare) / count * 1e6
        annot_us = ref(annotated) / count * 1e6
        each = {name: cost / count * 1e6 - bare_us for name, cost in blocks.items()}
        print(f"{fmt:18s} {bare_us:11.2f} {annot_us:9.2f} "
              f"{annot_us - bare_us:10.2f}  "
              f"{each['taxonomy']:8.2f} {each['comment']:9.2f} "
              f"{each['references']:7.2f} {sum(each.values()):8.2f}")

    print("A block is priced with only that block present, so `sum` and `pair")
    print("delta` measure the same thing by two routes; they land within a few")
    print("tenths of a microsecond, which is the corpus's own run-to-run spread.")
    print("The target the reader is held to is the `pair delta` column: what the")
    print("two headers cost, block for block.")

    print()
    print("Why the comment column is not the price of reading a comment: it is a")
    print("failed structured-comment regex, run once per record, whose pattern")
    print("backtracks over the first line of the block.  Quadratic in that line's")
    print("length, and paid by every GenBank record:")
    print(f"  {'first comment line':>18s} {'chars':>6s} {'re.search us':>13s} "
          f"{'us / char^2':>12s}")
    for nchars in (20, 40, 68, 120):
        line = "x" * (nchars - 1) + "!"
        attempts = 200_000
        start = time.perf_counter()
        for _ in range(attempts):
            _STRUCTURED_COMMENT.search(line)
        us = (time.perf_counter() - start) / attempts * 1e6
        print(f"  {'':>18s} {nchars:6d} {us:13.2f} {us / nchars ** 2:12.5f}")
    print("  The reference builds the pattern with an f-string on every call; the")
    print("  numbers above use a compiled one, so what they show is the matching")
    print("  and not the compilation, and the reference's cost is this or worse.")


def main() -> int:
    print(f"Python {sys.version.split()[0]}, median of {REPEATS} calls "
          f"after one warm-up\n")
    header = (f"{'corpus':18s} {'rec':>6s} {'MB':>7s} {'ref us/rec':>11s} "
              f"{'objects':>9s} {'scan':>7s} {'lines':>7s} {'ref/floor':>10s}")
    print(header)
    print("-" * len(header))

    for group, fmt, count, nbytes, text in seqio_corpus.corpus_rows():
        records = list(SeqIO.parse(io.StringIO(text), fmt))
        scalars = scalars_of(records)

        ref = timed(lambda: list(SeqIO.parse(io.StringIO(text), fmt)))
        obj = timed(lambda: build_objects(scalars), repeats=5)
        scan = timed(lambda: text.translate(_SCAN_TABLE).upper())
        lines = timed(lambda: [line for line in io.StringIO(text)])

        floor = obj + scan
        print(f"{group:18s} {count:6d} {nbytes / 1e6:7.3f} "
              f"{ref / count * 1e6:11.2f} {obj / count * 1e6:9.2f} "
              f"{scan / count * 1e6:7.2f} {lines / count * 1e6:7.2f} "
              f"{ref / floor:9.2f}x")

    print()
    print("ref/floor is the headroom: what a kernel could win if the reference's")
    print("own logic went to zero and the object tree were still built in Python.")

    print()
    print("The smallest rewrite there is -- ids and sequences, no features, no")
    print("annotations, no description.  Its ratio is a floor, not a replacement:")
    print("it is not asked to produce what the reference produces, and it does not.")
    header2 = (f"{'corpus':18s} {'ref us/rec':>11s} {'rewrite':>9s} "
               f"{'ratio':>8s}  {'ids and sequences identical':>28s}")
    print(header2)
    print("-" * len(header2))
    for group, fmt, count, nbytes, text in seqio_corpus.corpus_rows():
        rewrite = timed(lambda: minimal_parse(fmt, text))
        ref = timed(lambda: list(SeqIO.parse(io.StringIO(text), fmt)))
        same = minimal_parse(fmt, text) == reference_pairs(fmt, text)
        print(f"{group:18s} {ref / count * 1e6:11.2f} "
              f"{rewrite / count * 1e6:9.2f} {ref / rewrite:7.2f}x "
              f"{'yes' if same else 'NO':>28s}")

    annotation_costs()
    return 0


_SCAN_TABLE = {ord(c): None for c in "0123456789 \n\r"}

#: The reference's structured-comment probe, compiled once here so that the
#: numbers in `annotation_costs` measure the *matching* rather than the
#: `re` cache lookup the reference pays on top of it (`Scanner.py`, the
#: ``COMMENT`` branch of ``_feed_header_lines``).
_STRUCTURED_COMMENT = re.compile(r"([^#]+)-START##$")

if __name__ == "__main__":
    raise SystemExit(main())
