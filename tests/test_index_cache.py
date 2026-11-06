"""The index file: what it promises, and what it refuses to promise.

`open_fastq_index` writes ``path + ".bfqi"`` beside the data and maps it on the
next open, so that opening an index costs a mapping instead of a scan.  A cache
is the one thing in this package that can answer *wrongly* without any code
being wrong: an index built for one file, used for another, returns the bytes
of a record that used to be there.  So the tests here are less about speed --
that belongs in `bench/` -- than about the conditions under which the answer is
allowed to be reused at all:

* the same size and the same fingerprint of both ends of the data, or rebuild;
* every range inside the index file, or rebuild;
* the right magic and version, or rebuild;
* and never an exception, because a cache that cannot be written, read or
  trusted must not be able to refuse a file that is otherwise fine.

The equality tests are the other half: a warm open must answer exactly what a
cold one did, including for the records that are not `plain` and whose keys
live in the index's own side table rather than in the data.
"""

import gzip
import os
import struct
from pathlib import Path

import pytest

import biofasting

HEADER_VERSION_OFFSET = 8  # magic[8], then a little-endian uint32 version


def write(path, text):
    path.write_bytes(text.encode() if isinstance(text, str) else text)
    return path


def records(index):
    """Everything the public surface promises, in a form two indexes can compare."""
    return [
        (key, index.title(key), index[key], index.quality(key),
         index.sequence_length(key), index.quality_length(key),
         index.plain(key), index.sequence_slice(key, 1, 4))
        for key in index.keys()
    ]


GOOD = (
    "@r1 comment\nACGTACGT\n+\nIIIIIIII\n"
    "@r2\nTTTT\n+\nJJJJ\n"
)


# --- the happy path -------------------------------------------------------


def test_a_cold_open_writes_an_index_and_a_warm_open_reads_it(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    cold = biofasting.open_fastq_index(path)
    cache = tmp_path / "reads.fastq.bfqi"
    assert cache.exists()
    assert cache.stat().st_size > 0

    warm = biofasting.open_fastq_index(path)
    assert records(cold) == records(warm)


def test_the_warm_index_is_the_same_index_not_a_different_one(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    cold = biofasting.open_fastq_index(path)
    warm = biofasting.open_fastq_index(path)
    assert warm.keys() == cold.keys()
    assert [warm.index_entry(k) for k in cold.keys()] == [
        cold.index_entry(k) for k in cold.keys()]


def test_cache_false_writes_nothing(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    biofasting.open_fastq_index(path, cache=False)
    assert not (tmp_path / "reads.fastq.bfqi").exists()


def test_an_empty_file_indexes_to_an_empty_index_and_survives_the_round_trip(tmp_path):
    path = write(tmp_path / "empty.fastq", "")
    cold = biofasting.open_fastq_index(path)
    assert len(cold) == 0
    assert (tmp_path / "empty.fastq.bfqi").exists()
    warm = biofasting.open_fastq_index(path)
    assert len(warm) == 0
    assert warm.keys() == []


def test_the_keys_keep_file_order_across_the_round_trip(tmp_path):
    text = "".join("@read%d\nACGT\n+\nIIII\n" % i for i in range(50))
    path = write(tmp_path / "ordered.fastq", text)
    cold = biofasting.open_fastq_index(path)
    warm = biofasting.open_fastq_index(path)
    assert warm.keys() == ["read%d" % i for i in range(50)]
    assert warm.keys() == cold.keys()


# --- records the scanner had to assemble ---------------------------------
#
# These have no address in the data for their key or their fields, so their
# keys travel in the index file's own blob.  A warm open that lost them would
# still answer the plain records correctly, which is exactly why they are
# tested separately.

ASSEMBLED = (
    "@wrapped\nACGT\nACGT\n+\nIIII\nIIII\n"   # sequence and quality wrapped
    "@blank\nACGT\n+\nIIII\n"
    "\n"                                      # blank line between records
    "@plus title\nACGT\n+plus title\nIIII\n"   # '+' line repeats the header
)


def test_assembled_records_answer_the_same_after_a_round_trip(tmp_path):
    path = write(tmp_path / "assembled.fastq", ASSEMBLED)
    cold = biofasting.open_fastq_index(path)
    # The wrapped and blank-line records are assembled into scratch buffers and
    # so have no address of their own; the '+' repetition is still laid out at
    # fixed addresses, which is why it stays `plain`.
    assert [cold.plain(key) for key in cold.keys()] == [False, False, True]
    warm = biofasting.open_fastq_index(path)
    assert records(cold) == records(warm)
    assert warm["wrapped"] == b"ACGTACGT"
    assert warm["plus"] == b"ACGT"


def test_a_header_too_long_to_encode_in_the_index_does_not_break_the_open(tmp_path):
    # The on-disk entry keeps lengths in 16 bits, so a header longer than
    # 65535 bytes cannot be persisted.  The file is still readable: the index
    # writer refuses and the in-memory index answers.
    name = "r" * 70000
    path = write(tmp_path / "long.fastq", "@%s\nACGT\n+\nIIII\n" % name)
    index = biofasting.open_fastq_index(path)
    assert index.keys() == [name]
    assert index[name] == b"ACGT"
    assert not (tmp_path / "long.fastq.bfqi").exists()


# --- when the cache must not be trusted ---------------------------------


def test_an_index_for_different_bytes_is_rebuilt_not_used(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    biofasting.open_fastq_index(path)
    # Same name, different content and a different size, as after a rewrite.
    write(path, "@other\nAAAA\n+\nIIII\n")
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["other"]
    assert index["other"] == b"AAAA"
    # And the stale cache was replaced, not merely ignored.
    assert biofasting.open_fastq_index(path).keys() == ["other"]


def test_an_index_for_a_same_sized_but_different_file_is_rejected(tmp_path):
    path = write(tmp_path / "reads.fastq", "@r1\nACGT\n+\nIIII\n")
    biofasting.open_fastq_index(path)
    # Same length, different first byte: only the fingerprint can tell.
    write(path, "@z1\nACGT\n+\nIIII\n")
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["z1"]


def test_a_truncated_index_is_ignored(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    biofasting.open_fastq_index(path)
    cache = tmp_path / "reads.fastq.bfqi"
    raw = cache.read_bytes()
    cache.write_bytes(raw[: len(raw) // 2])
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["r1", "r2"]
    assert index["r1"] == b"ACGTACGT"


def test_an_index_from_another_version_is_ignored(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    biofasting.open_fastq_index(path)
    cache = tmp_path / "reads.fastq.bfqi"
    raw = bytearray(cache.read_bytes())
    struct.pack_into("<I", raw, HEADER_VERSION_OFFSET, 99)
    cache.write_bytes(bytes(raw))
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["r1", "r2"]


def test_a_file_that_is_not_an_index_at_all_is_ignored(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    write(tmp_path / "reads.fastq.bfqi", "this is not an index" * 8)
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["r1", "r2"]


def test_a_gzipped_file_caches_both_the_index_and_the_inflated_bytes(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(GOOD.encode(), mtime=0))
    cold = biofasting.open_fastq_index(path)
    assert (tmp_path / "reads.fastq.gz.bfqi").exists()
    # The index addresses the decompressed stream, so the stream has to exist
    # too -- the cache is the pair, not the index alone.
    inflated = tmp_path / "reads.fastq.gz.bfqd"
    assert inflated.read_bytes() == GOOD.encode()
    warm = biofasting.open_fastq_index(path)
    assert records(cold) == records(warm)


def test_replacing_the_archive_invalidates_the_cached_pair(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(GOOD.encode(), mtime=0))
    biofasting.open_fastq_index(path)
    replacement = "@other\nAAAA\n+\nIIII\n"
    path.write_bytes(gzip.compress(replacement.encode(), mtime=0))
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["other"]
    assert index["other"] == b"AAAA"
    assert (tmp_path / "reads.fastq.gz.bfqd").read_bytes() == replacement.encode()


def test_an_index_without_its_decompressed_copy_is_rebuilt(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(GOOD.encode(), mtime=0))
    biofasting.open_fastq_index(path)
    (tmp_path / "reads.fastq.gz.bfqd").unlink()
    index = biofasting.open_fastq_index(path)
    assert index.keys() == ["r1", "r2"]
    assert (tmp_path / "reads.fastq.gz.bfqd").exists()


def test_cache_false_writes_nothing_for_a_gzipped_file(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(GOOD.encode(), mtime=0))
    biofasting.open_fastq_index(path, cache=False)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["reads.fastq.gz"]


# --- the failures that must not become refusals --------------------------


def test_a_read_only_directory_still_opens_the_file(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    os.chmod(tmp_path, 0o555)
    try:
        index = biofasting.open_fastq_index(path)
        assert index.keys() == ["r1", "r2"]
    finally:
        os.chmod(tmp_path, 0o755)
    assert not (tmp_path / "reads.fastq.bfqi").exists()


def test_a_malformed_file_is_refused_and_writes_no_index(tmp_path):
    path = write(tmp_path / "bad.fastq", "@r1\nACGT\n+\nII\n")
    with pytest.raises(ValueError):
        biofasting.open_fastq_index(path)
    assert not (tmp_path / "bad.fastq.bfqi").exists()


def test_a_duplicate_key_is_refused_and_writes_no_index(tmp_path):
    path = write(tmp_path / "dup.fastq", "@r1\nACGT\n+\nIIII\n@r1\nTTTT\n+\nIIII\n")
    with pytest.raises(ValueError, match="Duplicate key 'r1'"):
        biofasting.open_fastq_index(path)
    assert not (tmp_path / "dup.fastq.bfqi").exists()


def test_no_temporary_file_is_left_behind(tmp_path):
    path = write(tmp_path / "reads.fastq", GOOD)
    biofasting.open_fastq_index(path)
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "reads.fastq", "reads.fastq.bfqi"]


def test_the_data_mapping_outlives_the_index_object(tmp_path):
    # A cache adds a second mapping; neither may be released while a fetch can
    # still read through them.
    path = write(tmp_path / "reads.fastq", GOOD)
    biofasting.open_fastq_index(path)
    index = biofasting.open_fastq_index(path)
    assert index["r2"] == b"TTTT"
