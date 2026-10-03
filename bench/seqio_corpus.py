#!/usr/bin/env python3
"""Generated corpora for the flat-file sequence formats (PLAN 2.4 ranking).

``bench/data/`` is FASTQ and FASTA.  The formats this module covers -- GenBank,
EMBL and SwissProt -- are not in it, and they are not sampled either: they are
``bench/seqio_corpus.py`` is **generated, not stored**, the way the protein
corpus of M18 is.  Same contract as ``bench/gen_data.py``: SHA-256 in counter
mode keyed by ``(seed, domain)``, every sampling decision in integer arithmetic,
no float threshold anywhere.  The same command yields the same records on any
machine.

One caveat that belongs in the open rather than in a footnote: GenBank and EMBL
are emitted by ``Bio.SeqIO.write``, because those two formats have a writer and
hand-writing a parser-correct INSDC record is a way to spend a week finding out
that column 22 is not where you thought it was.  SwissProt has **no writer**
upstream -- ``SeqIO.write`` answers ``Reading format 'swiss' is supported, but
not writing`` -- so that emitter is written out here, and the parser that
accepts it is the check on it.

Emitting the first two through the reference means the corpus follows the
reference's writer, including whatever it does that a real GenBank file does
not.  That is acceptable for what these files are for: the ranking pass asks
what the *parser* costs, and a corpus the reference's own writer produced is a
corpus the reference's own parser is guaranteed to accept.
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

HERE = Path(__file__).parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gen_data import Keystream, base_table  # noqa: E402

SEQIO_SEED = 20261004

_DNA_TABLE = base_table(50)  # 64 A, 64 C, 64 G, 64 T -- uniform, integer-built

#: One row per workload: ``(group, format, records, bases, features per kb)``.
#: Small records expose the per-record cost, large ones the per-byte cost, and
#: the feature density is the bacterial-genome figure -- about one gene per
#: kilobase -- because features are what a GenBank parse actually spends its
#: time on.  ``genbank-1kb-bare`` is the same sequences as ``genbank-1kb`` with
#: no features at all: the difference between those two rows is the cost of
#: building ``SeqFeature`` objects, isolated rather than argued about.
#:
#: The last column is the **annotation density**: 0 is the thin header a bare
#: ``SeqRecord`` writes (LOCUS, DEFINITION, ACCESSION, VERSION, KEYWORDS, SOURCE,
#: ORGANISM), 1 adds the blocks a real submission carries -- a taxonomy, a
#: keyword list, two accessions, a sequence version, a multi-line comment and
#: three references.
CORPORA: tuple[tuple[str, str, int, int, int, int], ...] = (
    ("genbank-1kb", "genbank", 3_000, 1_000, 1, 0),
    ("genbank-10kb", "genbank", 300, 10_000, 1, 0),
    ("genbank-1kb-bare", "genbank", 300, 1_000, 0, 0),
    ("embl-1kb", "embl", 3_000, 1_000, 1, 0),
    ("embl-10kb", "embl", 300, 10_000, 1, 0),
    ("swiss-300aa", "swiss", 5_000, 300, 2, 0),
)

#: The **annotation-density** rows, kept in their own tuple so that adding them
#: does not change the six rows every already-recorded number in
#: `bench/targets.md` was measured on: each pair here is the same sequences from
#: the same seed and the same domain string as its bare twin, with the header as
#: the only difference, so the difference of the two per-record times is the cost
#: of the annotations dict and nothing else.  `genbank-1kb-bare` is not repeated
#: because it is already a `CORPORA` row; the EMBL bare row is, because EMBL had
#: no bare row -- and it would have changed nothing but the table to add one.
ANNOT_CORPORA: tuple[tuple[str, str, int, int, int, int], ...] = (
    ("genbank-1kb-annot", "genbank", 300, 1_000, 0, 1),
    ("embl-1kb-bare", "embl", 300, 1_000, 0, 0),
    ("embl-1kb-annot", "embl", 300, 1_000, 0, 1),
)


def _aa_table() -> bytes:
    """256-byte table over the twenty standard residues, built in integers.

    The same shape as ``gen_data.base_table``: 256 does not divide by twenty, so
    sixteen letters get thirteen slots and four get twelve.  The unevenness is
    deterministic and worth less than a percent; what matters is that no float
    decides it.
    """
    table = bytearray(256)
    pos = 0
    for i, letter in enumerate(b"ACDEFGHIKLMNPQRSTVWY"):
        n = 256 // 20 + (1 if i < 256 % 20 else 0)
        table[pos : pos + n] = bytes([letter]) * n
        pos += n
    return bytes(table)


_AA_TABLE = _aa_table()


def dna(ks: Keystream, length: int) -> str:
    """``length`` bases from the keystream, through the integer base table."""
    return ks.take(length).translate(_DNA_TABLE).decode("ascii")


def protein(ks: Keystream, length: int) -> str:
    return ks.take(length).translate(_AA_TABLE).decode("ascii")


def _features(ks: Keystream, length: int, per_kb: int) -> list:
    """Deterministic features, evenly spaced, never overlapping the ends.

    The positions come from the keystream rather than from a counter so that a
    record's features are not a function of its index alone, and the strand is a
    single keystream byte with a threshold at 128 -- integer, as everywhere
    else.  A third of the features are joins, because a compound location is a
    different code path in the parser and a corpus of single spans would measure
    only the easy one.
    """
    from Bio.SeqFeature import SeqFeature, SimpleLocation

    count = max(1, length * per_kb // 1000)
    step = length // (count + 1)
    span = min(step - 20, 900)
    if span < 30:
        span = 30
    out = []
    for i in range(count):
        start = step * (i + 1)
        stop = min(start + span, length - 1)
        if stop - start < 20:
            continue
        noise = ks.take(2)
        strand = 1 if noise[0] < 128 else -1
        if noise[1] < 85 and stop - start > 120:
            cut = start + (stop - start) // 2
            loc = SimpleLocation(start, cut, strand=strand) + SimpleLocation(
                cut + 5, stop, strand=strand
            )
            ftype = "CDS"
        else:
            loc = SimpleLocation(start, stop, strand=strand)
            ftype = "gene" if noise[1] < 170 else "misc_feature"
        qualifiers = {"note": [f"synthetic feature {i + 1}"]}
        if ftype == "CDS":
            qualifiers["product"] = ["hypothetical protein"]
            qualifiers["codon_start"] = ["1"]
            qualifiers["transl_table"] = ["11"]
        else:
            qualifiers["gene"] = [f"syn{start:07d}"]
        out.append(SeqFeature(loc, type=ftype, qualifiers=qualifiers))
    return out


def _annotate(record) -> None:
    """Give one record the annotation blocks a real submission carries.

    Fixed text with the record's index in it rather than a keystream draw: the
    values are what the reference's *writer* formats, and the writer is the thing
    under test here, so there is nothing to sample.  The blocks are the ones the
    thin header has not got -- a taxonomy, a keyword list, a comment long enough
    to wrap, and three references.

    Two things that look like they belong here are deliberately absent, because
    the writer drops them and setting them would be a line of this function that
    does nothing: a **secondary accession** (`SeqIO.write` emits only the first
    accession, and it takes it from the singular ``annotations["accession"]``,
    not from the ``accessions`` list a reader fills in), and a
    **sequence version** (the VERSION line's ``.N`` comes from a dot in the
    record's *id*, so a record whose id has none writes a versionless VERSION
    line however the annotation is set).
    """
    from Bio.SeqFeature import Reference

    record.annotations["taxonomy"] = ["artificial sequences", "synthetic construct"]
    record.annotations["keywords"] = ["synthetic", "benchmark"]
    record.annotations["comment"] = (
        "This record is generated for benchmarking and carries a comment long "
        "enough to span several lines of the flat file, which the reference "
        "reads back with the line breaks it wrote."
    )
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


def _records(fmt: str, count: int, length: int, per_kb: int, annot: int = 0) -> list:
    """Build the ``SeqRecord`` list for one corpus row, deterministically."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    ks = Keystream(SEQIO_SEED, f"seqio/{fmt}/{count}x{length}".encode())
    #   The feature layout draws from its own stream.  Sharing one with the
    #   sequences would make `genbank-1kb` and `genbank-1kb-bare` diverge from
    #   the second record on, and the whole point of that pair is that the
    #   residues are identical and the feature objects are the only difference.
    feature_ks = Keystream(SEQIO_SEED,
                           f"seqio/{fmt}/{count}x{length}/features".encode())
    records = []
    for i in range(count):
        if fmt == "swiss":
            seq = protein(ks, length)
        else:
            seq = dna(ks, length)
        tag = f"SYN{i + 1:06d}"
        rec = SeqRecord(Seq(seq), id=tag, name=tag,
                        description=f"synthetic {fmt} record {i + 1} of {count}")
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["data_file_division"] = "PLN"
        rec.annotations["date"] = "01-JAN-2026"
        rec.annotations["organism"] = "synthetic construct"
        rec.annotations["topology"] = "linear"
        #   Annotation density is the header only, and the header is applied
        #   after the residues are drawn, so the bare and annotated rows share
        #   their sequences base for base -- the property the measurement rests
        #   on, and the reason this cannot move above the sequence draw.
        if annot:
            _annotate(rec)
        if per_kb:
            rec.features = _features(feature_ks, length, per_kb)
        records.append(rec)
    return records


def _swiss_text(records: list, per_record: int) -> str:
    """A SwissProt flat file, written out by hand because upstream cannot.

    ``per_record`` features are emitted per entry as ``FT`` blocks.  A UniProt
    entry without them is not a realistic input: the parser turns every ``FT``
    block into a ``SeqFeature``, so a corpus of bare ``SQ`` blocks would measure
    the cheap half of the format and report a number that no real file would
    produce.
    """
    from Bio.SeqUtils import molecular_weight

    ks = Keystream(SEQIO_SEED, b"seqio/swiss/ft")
    lines: list[str] = []
    for rec in records:
        seq = str(rec.seq)
        #   Every field the parser reads is present and in its column.  The
        #   molecular weight is the reference's own, so the header is not a
        #   number this project invented; the CRC64 field is a placeholder,
        #   which is what it is for -- nothing parses it back.
        mw = int(round(molecular_weight(seq, seq_type="protein", monoisotopic=False)))
        lines.append(f"ID   {rec.id:<16} Reviewed;         {len(seq)} AA.")
        lines.append(f"AC   {rec.id};")
        lines.append("DT   01-JAN-2026, integrated into UniProtKB/Swiss-Prot.")
        lines.append(f"DE   RecName: Full={rec.description};")
        lines.append(f"GN   Name={rec.id};")
        lines.append("OS   Synthetic construct.")
        lines.append("OC   artificial sequences; synthetic construct.")
        lines.append("OX   NCBI_TaxID=32630;")
        lines.append("FT   CHAIN           1..%d" % len(seq))
        lines.append(f'FT                   /id="PRO_{rec.id}"')
        for i in range(per_record):
            start = 10 + (ks.take(1)[0] * (len(seq) - 40) // 256)
            stop = min(start + 20 + ks.take(1)[0] // 4, len(seq))
            lines.append(f"FT   DOMAIN          {start}..{stop}")
            lines.append(f'FT                   /note="synthetic domain {i + 1}"')
        lines.append(f"SQ   SEQUENCE   {len(seq)} AA;  {mw} MW;  "
                     f"{rec.id:>16} CRC64;")
        for i in range(0, len(seq), 60):
            chunk = seq[i : i + 60]
            groups = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
            lines.append("     " + groups)
        lines.append("//")
    return "\n".join(lines) + "\n"


def build(fmt: str, count: int, length: int, per_kb: int, annot: int = 0) -> str:
    """The whole file for one corpus row, as text."""
    records = _records(fmt, count, length, per_kb, annot)
    if fmt == "swiss":
        #   The last column is features per *kilobase* for the nucleotide
        #   formats and features per *record* here, because a protein entry is a
        #   third of a kilobase and "per kb" would round it to nothing.
        return _swiss_text(records, per_kb)
    from Bio import SeqIO

    handle = io.StringIO()
    SeqIO.write(records, handle, fmt)
    return handle.getvalue()


def corpus_rows() -> list[tuple[str, str, int, int, str]]:
    """``(group, format, records, bytes, text)`` for every row of ``CORPORA``."""
    rows = []
    for group, fmt, count, length, per_kb, annot in CORPORA:
        text = build(fmt, count, length, per_kb, annot)
        rows.append((group, fmt, count, len(text.encode("utf-8")), text))
    return rows


def spec(group: str) -> tuple[str, str, int, int, int, int]:
    """The row spec for ``group``, from either corpus.

    One lookup over both tuples, because the annotation rows are rows of the same
    corpus built by the same generator and a second way of finding them would be
    a second place for the two to drift.
    """
    for row in CORPORA + ANNOT_CORPORA:
        if row[0] == group:
            return row
    raise KeyError(group)


def build_row(group: str) -> tuple[str, str]:
    """The ``(format, text)`` of one row, by name."""
    _, fmt, count, length, per_kb, annot = spec(group)
    return fmt, build(fmt, count, length, per_kb, annot)
