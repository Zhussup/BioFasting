"""`Bio.SeqUtils.CheckSum`: crc32, crc64, gcg and SEGUID.

Four checksums, and they split three ways on the only question that matters
here: what is a character?

`crc32` and `seguid` are already C.  `binascii.crc32` and `hashlib.sha1` are the
whole implementation on the reference's side too, so these are ports of the
reference's *coercion* and nothing else -- and the coercion is the subtle part.
Both do ``try: bytes(seq) except TypeError: seq.encode()``, so a `Seq` is read
as its bytes while a `str` is read as UTF-8, and ``crc32("é")`` is a checksum of
two bytes.

`crc64` is a loop over the sequence that folds each character in with
``ord(c) & 0xFF``.  That is arithmetic on the low byte, so for text that is
Latin-1 -- which is how this package reads a byte buffer, and how `Seq` reads
one -- the character and the byte are the same number and the kernel is exact
with no restriction at all.  Only a character above U+00FF, which cannot be
encoded as one byte, has to go back to Python.

`gcg` upper-cases each character first, and `str.upper` is not a byte operation:
``"é".upper()`` is ``"É"`` and ``"ß".upper()`` is two characters, which is why
the reference raises ``TypeError: ord() expected a character, but string of
length 2 found`` on the second.  So this one **declines** on any byte at or
above 0x80 and the caller runs the reference's own loop, which is the same
arrangement :func:`biofasting.translate` uses for a codon its tables cannot
classify.
"""

import base64
import binascii
import hashlib

from . import _core
from .sequtils import _as_text

__all__ = ["crc32", "crc64", "gcg", "seguid"]


def crc32(seq):
    """The CRC-32 of `seq`, as an int.

    >>> crc32("ACGTACGTACGT")
    20049947

    Case matters.  A `Seq` is read as its bytes; a `str` is encoded as UTF-8,
    because that is what the reference calls and what makes
    ``crc32("é")`` a checksum over two bytes.
    """
    try:
        s = bytes(seq)
    except TypeError:
        s = seq.encode()
    return binascii.crc32(s)


def crc64(s):
    """The CRC-64 of `s`, as the reference's ``"CRC-XXXXXXXXXXXXXXXX"`` string.

    >>> crc64("ACGTACGTACGT")
    'CRC-C4FBB762C4A87EBD'

    Case matters here too, unlike in :func:`gcg` and :func:`seguid`.
    """
    text = _as_text(s)
    try:
        buffer = text.encode("latin-1")
    except UnicodeEncodeError:
        #   A character above U+00FF has no low byte to fold in and no
        #   one-byte encoding, so the byte kernel has nothing to say about it.
        #   The reference's loop over the characters does, and this is it.
        return _reference_crc64(text)
    high, low = _core.crc64(buffer)
    return f"CRC-{high:08X}{low:08X}"


def gcg(seq):
    """The GCG checksum of `seq`, as an int.

    >>> gcg("ACGTACGTACGT")
    5688

    Case-insensitive: every character is upper-cased first, which is why this
    one falls back to Python for a non-ASCII sequence instead of answering
    something a byte loop can reach.
    """
    text = _as_text(seq)
    if text.isascii():
        checksum = _core.gcg(text.encode("latin-1"))
        if checksum is not None:
            return checksum
    return _reference_gcg(text)


def seguid(seq):
    """The SEGUID of `seq`, as a string.

    >>> seguid("ACGTACGTACGT")
    'If6HIvcnRSQDVNiAoefAzySc6i4'

    Case does not matter, and a `Seq` and a `str` give the same answer only when
    the string is ASCII: the reference upper-cases the *encoded* bytes, so a
    `str` with a character above ASCII is upper-cased as UTF-8.
    """
    m = hashlib.sha1()
    try:
        seq = bytes(seq)
    except TypeError:
        seq = seq.encode()
    m.update(seq.upper())
    tmp = base64.encodebytes(m.digest())
    return tmp.decode().replace("\n", "").rstrip("=")


def _table_h():
    """The reference's ``CheckSum._table_h``, as the kernel computed it."""
    return list(_core.crc64_table_h())


def _reference_crc64(text):
    """`crc64`'s loop, for characters a byte kernel cannot fold in."""
    table = _core.crc64_table_h()
    crcl = crch = 0
    for c in text:
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        crch = temp1h ^ table[(crcl ^ ord(c)) & 0xFF]
        crcl = temp1l
    return f"CRC-{crch:08X}{crcl:08X}"


def _reference_gcg(text):
    """`gcg`'s loop, kept only for sequences a byte kernel cannot upper-case."""
    index = checksum = 0
    for char in text:
        index += 1
        checksum += index * ord(char.upper())
        if index == 57:
            index = 0
    return checksum % 10000
