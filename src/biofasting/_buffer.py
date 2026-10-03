"""Mapping a file read-only, without reading it.

One helper, two callers, and the reason it is its own module rather than a
function in either of them is that the empty-file case is the part that is
easy to get wrong: `mmap` refuses a zero-length file, and a zero-length file is
a valid (empty) input for both formats rather than an error.
"""

import mmap
import os

__all__ = ["map_readonly"]


def map_readonly(path):
    """Return the file at `path` as a read-only buffer.

    A memory map, so a 100 MB genome costs a few microseconds and no resident
    memory beyond the pages actually touched; ``b""`` for an empty file, which
    cannot be mapped.  The caller must keep the result alive for as long as it
    reads from it -- in practice the compiled object that holds a buffer view
    of it does that.
    """
    fd = os.open(os.fspath(path), os.O_RDONLY)
    try:
        size = os.fstat(fd).st_size
        return mmap.mmap(fd, 0, access=mmap.ACCESS_READ) if size else b""
    finally:
        os.close(fd)
