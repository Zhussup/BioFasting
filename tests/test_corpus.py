"""Step 1.5: the whole corpus, in one place, with nothing left out.

The per-kernel test files each cover the files that exercise their kernel.
This one exists because "we tested the corpus" is a claim about a *set*, and a
set is not verified by testing each member somebody remembered to list.  So the
files here are discovered by walking `bench/data/`, the walk is required to
equal the manifest that generated them, and the manifest's own SHA-256s and
record counts are checked before any parser is allowed to make a claim about
the bytes.  A new corpus file that nobody test-covers therefore shows up as a
failure, not as silence.

Two layers again:

- **Verdicts.** Every file, our reader against the reference, element by
  element, accept or reject -- including the four malformed files, where the
  interesting assertion is that both sides refuse.
- **Mutations.** The corpus files are real input; random line-shaped text is
  not.  So the fuzz takes the bytes of real records and damages them -- a byte
  flipped, a line dropped, the file truncated mid-record -- which is what a
  truncated download or a bad disk actually produces, and demands the same
  verdict from both parsers on every mutant.

Only ASCII is ever produced or fed in: the reference reads a *text* stream and
decodes it, while our readers work on bytes and do not decode, so a non-ASCII
byte is a difference in contract rather than a bug in either.
"""

import gzip
import hashlib
import importlib.util
import io
import json
import random
import re
from pathlib import Path

import pytest

import biofasting

DATA = Path(__file__).resolve().parent.parent / "bench" / "data"

needs_corpus = pytest.mark.skipif(
    not (DATA / "MANIFEST.json").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)
needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

# The manifest lists what the generator wrote.  These are the files in the tree
# that were written by something else -- the checksum file itself, and the
# pyfaidx index that is deliberately kept as a third implementation's opinion.
NOT_GENERATED = {"MANIFEST.json"}

# One entry of the manifest on disk overcounts its file, and the tests below
# know it rather than rediscovering it.  `gen_data.py` counted records with
# `payload.count(b"\n@")`, which counts a *quality* line that begins with '@'
# as a header -- precisely the trap `quality_contains_at_and_plus.fastq` exists
# to test, so the file that catches parsers counting '@' also caught the
# generator doing it.  Both parsers agree the file has 2 records; the manifest
# says 3.  The generator now declares the count for that file explicitly, so
# regenerating the corpus produces a manifest that needs no adjustment here --
# until then, the difference is asserted rather than tolerated.
MANIFEST_OVERCOUNTS = {
    "edge/quality_contains_at_and_plus.fastq": (
        "a quality line of '@'s is counted as a record header"
    ),
}


def manifest():
    return json.loads((DATA / "MANIFEST.json").read_text())


def manifest_entries():
    """The manifest's file list, or nothing at all when there is no corpus.

    The parametrised tests are collected before any skip marker is evaluated,
    so a module-level read of the manifest would turn "the corpus is not
    generated" from a skip into a collection error -- which is exactly how CI,
    where the corpus never exists, would see it.
    """
    path = DATA / "MANIFEST.json"
    return json.loads(path.read_text())["files"] if path.exists() else []


def corpus_files():
    """Every generated file in the corpus, as paths relative to `bench/data`."""
    return sorted(
        str(path.relative_to(DATA))
        for path in DATA.rglob("*")
        if path.is_file()
        and path.name not in NOT_GENERATED
        and not path.name.endswith(".fai")
    )


def read_bytes(relative):
    """The file's bytes, gunzipped when it is a `.gz` -- as a parser sees them."""
    path = DATA / relative
    return gzip.decompress(path.read_bytes()) if relative.endswith(".gz") else path.read_bytes()


def our_records(relative):
    """Our reader's records for `relative`, streamed, in file order."""
    path = DATA / relative
    if relative.endswith(".fastq") or relative.endswith(".fastq.gz"):
        return biofasting.open_fastq(path)
    index = biofasting.open_fasta(path)
    # The title is encoded because the reference yields bytes; everything else
    # in this file compares bytes to bytes, and a str title would compare
    # unequal to an identical one from the other side.
    return ((index.title(key).encode(), index[key]) for key in index)


def reference_records(relative):
    """The reference's records for `relative`, streamed, as bytes."""
    path = DATA / relative
    if relative.endswith(".fastq") or relative.endswith(".fastq.gz"):
        from Bio.SeqIO.QualityIO import FastqGeneralIterator

        opener = gzip.open if relative.endswith(".gz") else open
        with opener(path, "rt") as handle:
            for title, sequence, quality in FastqGeneralIterator(handle):
                yield title.encode(), sequence.encode(), quality.encode()
    else:
        from Bio.SeqIO.FastaIO import SimpleFastaParser

        with open(path, newline=None) as handle:
            for title, sequence in SimpleFastaParser(handle):
                yield title.encode(), sequence.encode()


_EXHAUSTED = object()


def step(iterator):
    """Pull one item: `(item, error)`, with `_EXHAUSTED` for a clean end."""
    try:
        return next(iterator), None
    except StopIteration:
        return _EXHAUSTED, None
    except ValueError as exc:  # exactly the type both sides raise
        return _EXHAUSTED, exc


def verdict(relative):
    """Walk our reader and the reference together and report the first difference.

    Lockstep is the point: a differential that compares two materialised lists
    can only say "these differ", while walking both at once says at which
    record, and a parser that is wrong from record 400 000 on is a different
    bug from one that is wrong at record 1.  It also never holds two copies of
    a 369 MB file in memory at the same time.
    """
    ours, theirs = our_records(relative), reference_records(relative)
    index = 0
    while True:
        mine, my_error = step(ours)
        reference, their_error = step(theirs)
        if my_error is not None or their_error is not None:
            return index, my_error, their_error
        if mine is _EXHAUSTED and reference is _EXHAUSTED:
            return index, None, None
        assert mine == reference, f"{relative}: record {index} differs: {mine!r} vs {reference!r}"
        index += 1


# --------------------------------------------------------------------------
# The corpus is the corpus.  Every parser claim below rests on these bytes, so
# they are pinned before anything reads them.
# --------------------------------------------------------------------------


@needs_corpus
def test_the_files_on_disk_are_exactly_the_files_in_the_manifest():
    """Nothing extra, nothing missing.

    The corpus is gitignored and regenerable, so "the corpus" is not a fixed
    thing a reader can assume -- it is whatever `bench/gen_data.py` last wrote.
    Requiring the tree to equal the manifest is what makes the parametrised
    sweep below a statement about a known set rather than about a directory
    listing that happened to look reasonable.
    """
    listed = {entry["path"] for entry in manifest()["files"]}
    assert set(corpus_files()) == listed


@needs_corpus
def test_every_corpus_file_has_the_checksum_and_size_the_generator_recorded():
    """The manifest's digest is over the *content*, not over the stored bytes.

    For the gzipped file those differ -- 369 MB of reads in a 172 MB file --
    and `gen_data.py --verify` hashes through the decompressor, so a check that
    hashed the file on disk would disagree with the generator about a file both
    of them are happy with.  The on-disk size is checked separately, against
    the field that means it.
    """
    for entry in manifest()["files"]:
        content = read_bytes(entry["path"])
        assert len(content) == entry["bytes"], entry["path"]
        assert hashlib.sha256(content).hexdigest() == entry["sha256"], entry["path"]
        assert (DATA / entry["path"]).stat().st_size == entry["stored_bytes"], entry["path"]


@needs_biopython
@needs_corpus
@pytest.mark.parametrize("relative", corpus_files())
def test_every_corpus_file_yields_the_record_count_the_generator_recorded(relative):
    """The count is a third opinion, and it is neither parser's.

    `gen_data.py` recorded how many records it wrote.  A parser that agrees
    with the reference because both were fed the same truncated view would
    still disagree with this.  A malformed file is the exception and says so:
    the count there is the number of records the generator *intended*, and the
    parsers stop before reaching it -- which is what the verdict test asserts.
    """
    count, my_error, _ = verdict(relative)
    if relative.startswith("malformed/"):
        assert my_error is not None
        return
    assert my_error is None
    entry = next(e for e in manifest()["files"] if e["path"] == relative)
    expected = entry["records"] - (1 if relative in MANIFEST_OVERCOUNTS else 0)
    assert count == expected, MANIFEST_OVERCOUNTS.get(relative, "")


# --------------------------------------------------------------------------
# Verdicts, over every file.
# --------------------------------------------------------------------------


@needs_biopython
@needs_corpus
@pytest.mark.parametrize("relative", corpus_files())
def test_every_corpus_file_gets_the_same_verdict(relative):
    """Records agree one for one, or both sides refuse -- never one of each.

    The `malformed/` files come through here too, and for them the assertion
    that matters is that *both* raise.  Ours accepting a file Biopython refuses
    is the failure mode that corrupts data quietly, and it is the reason the
    comparison is on the error as well as on the records.
    """
    count, my_error, their_error = verdict(relative)
    assert (my_error is None) == (their_error is None), (
        f"{relative}: ours={my_error!r}, reference={their_error!r}"
    )
    if my_error is None:
        entry = next(e for e in manifest()["files"] if e["path"] == relative)
        expected = entry["records"] - (1 if relative in MANIFEST_OVERCOUNTS else 0)
        assert count == expected, MANIFEST_OVERCOUNTS.get(relative, "")


@needs_biopython
@needs_corpus
@pytest.mark.parametrize(
    "relative",
    [e["path"] for e in manifest_entries() if e["path"].startswith("malformed/")],
)
def test_the_malformed_files_are_refused_with_a_reason(relative):
    """A refusal has to be diagnosable, not merely present.

    The message is not compared word for word -- it is our own wording in
    places -- but it must name the record it is about, because the whole value
    of rejecting a bad file is telling the user where in it to look.
    """
    _, my_error, their_error = verdict(relative)
    assert my_error is not None, f"{relative}: we accepted a malformed file"
    assert their_error is not None, f"{relative}: the reference accepted it"
    assert isinstance(my_error, ValueError)
    # Not the wording -- the location.  A rejection that does not say where is
    # only marginally better than a wrong answer, because the user's next move
    # is to open the file and look, and the message should tell them where.
    # Biopython's messages carry neither the record number nor the offset.
    assert re.search(r"record \d+, byte offset \d+", str(my_error)), (
        f"{relative}: the refusal does not locate itself: {my_error}"
    )


# --------------------------------------------------------------------------
# Mutations of real records.
# --------------------------------------------------------------------------


MUTATION_ALPHABET = b"@+ACGTNR \t\r\n;>"


def mutate(data, rng):
    """Damage `data` the way a real file gets damaged, and return the result."""
    data = bytearray(data)
    for _ in range(rng.randrange(1, 4)):
        choice = rng.randrange(6)
        if not data:
            break
        position = rng.randrange(len(data))
        if choice == 0:  # a flipped byte
            data[position] = rng.choice(MUTATION_ALPHABET)
        elif choice == 1:  # a dropped byte
            del data[position]
        elif choice == 2:  # an inserted byte
            data.insert(position, rng.choice(MUTATION_ALPHABET))
        elif choice == 3:  # a truncated download
            del data[position:]
        elif choice == 4:  # a duplicated line
            lines = bytes(data).split(b"\n")
            if lines:
                lines.insert(rng.randrange(len(lines)), rng.choice(lines))
                data = bytearray(b"\n".join(lines))
        else:  # a lost line ending, which merges two lines into one
            starts = [i for i, byte in enumerate(data) if byte == 0x0A]
            if starts:
                del data[rng.choice(starts)]
    return bytes(data)


def mutated_corpus(seed, count):
    """`count` mutants of real corpus files, each labelled with its parent."""
    sources = [
        relative
        for relative in corpus_files()
        if relative.startswith("edge/") or relative.startswith("malformed/")
    ]
    rng = random.Random(seed)
    for index in range(count):
        relative = sources[rng.randrange(len(sources))]
        yield relative, mutate(read_bytes(relative), rng)


def key_of(title):
    """The index's key rule, taken from the reference's own source.

    ``line[1:].strip().split(None, 1)[0]``, which is what `SeqIO.index` does
    and what the kernel is held to.  Written out here rather than imported so
    that the fuzz has a rule to check against and not just another parser.
    """
    stripped = title.strip()
    return stripped.split(None, 1)[0] if stripped else ""


def reference_on_bytes(relative, data):
    """The reference's verdict on mutated FASTQ bytes: records, or the error."""
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    handle = io.StringIO(data.decode("ascii"), newline=None)
    try:
        return [
            (t.encode(), s.encode(), q.encode())
            for t, s, q in FastqGeneralIterator(handle)
        ], None
    except ValueError as exc:
        return None, exc


def ours_on_bytes(relative, data):
    try:
        return list(biofasting._core.FastqScanner(data)), None
    except ValueError as exc:
        return None, exc


@needs_biopython
@needs_corpus
@pytest.mark.parametrize("seed", range(8))
def test_a_mutated_fastq_gets_the_same_verdict(seed):
    """Damaged real input, both parsers, same answer.

    Random line soup exercises the grammar; this exercises the *data path*,
    which is where a fast path that trusts its input goes wrong.  The mutants
    come from the corpus's own edge files, so a truncation lands in the middle
    of a real record rather than in the middle of a plausible-looking one.
    """
    for relative, data in mutated_corpus(seed, 150):
        base = relative.split("/")[-1]
        if not (base.endswith(".fastq") or base.endswith(".fastq.gz")):
            continue
        expected, expected_error = reference_on_bytes(relative, data)
        actual, actual_error = ours_on_bytes(relative, data)
        assert (expected_error is None) == (actual_error is None), (
            f"seed={seed} {relative} data={data!r}\n"
            f"reference={expected_error!r}\nours={actual_error!r}"
        )
        assert actual == expected, f"seed={seed} {relative} data={data!r}"


@needs_biopython
@needs_corpus
@pytest.mark.parametrize("seed", range(8))
def test_a_mutated_fasta_gets_the_same_records(seed):
    """The same treatment for FASTA, through the index a user would build.

    One divergence is expected here and is asserted rather than tolerated: a
    mutation can give two records the same key, and `SimpleFastaParser` is a
    parser with no opinion about that while the index refuses it -- the same
    rule `SeqIO.index` applies.  The mutant must then be refused, and the test
    says which of the two contracts it is checking instead of accepting either
    answer.
    """
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    for relative, data in mutated_corpus(seed, 150):
        if not relative.endswith(".fasta"):
            continue
        try:
            index = biofasting._core.FastaIndex(data)
            actual = [(index.title(key).encode(), index[key]) for key in index]
            actual_error = None
        except ValueError as exc:
            actual, actual_error = None, exc

        expected = [
            (title.encode(), sequence.encode())
            for title, sequence in SimpleFastaParser(
                io.StringIO(data.decode("ascii"), newline=None)
            )
        ]
        keys = [key_of(title.decode("ascii")) for title, _ in expected]
        duplicated = len(set(keys)) != len(keys)
        assert (actual_error is not None) == duplicated, (
            f"seed={seed} {relative} data={data!r}\n"
            f"duplicate keys={duplicated}, ours={actual_error!r}"
        )
        if not duplicated:
            assert actual == expected, f"seed={seed} {relative} data={data!r}"
