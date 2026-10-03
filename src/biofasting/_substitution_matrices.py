"""Substitution matrices: the array, the parser, and what a row of data means.

This is `Bio.Align.substitution_matrices` with the thirty matrices taken out of
it.  The triage calls the module `data` -- "becomes data + a parser; no class
port" -- and the split is the one M9 and M11 drew: the *matrices* are data and
live in `_substitution_matrices_data.py`; the `Array` class is behaviour and is
written out here, because no table can express what letter indexing does, when
it raises, or what `format()` prints.

The reference implements `Array` twice.  The class in
`Bio/Align/substitution_matrices/__init__.py` is a subclass of a C type,
`_arraycore.Array`, and the C type is where the `alphabet` attribute really
lives: a set-once property that validates the alphabet's length against the
array's shape, and that rejects a repeated letter.  Those checks are not
incidental -- they are why `matrix[0, 0:2]` raises instead of quietly returning a
view whose alphabet no longer describes it -- so they are reimplemented here in
Python, message for message.  The C file ships with Biopython (it is in the
installed package next to the `.so`), so every message in this file was read off
it rather than guessed.

This file has no imports but `numpy`, so `tools/gen_substitution_matrices.py` can
load it by path and use it to check the table it is about to write, without
importing the table it wrote last time.
"""

from __future__ import annotations

import contextlib
import string

import numpy as np


@contextlib.contextmanager
def as_handle(handleish, mode="r", **kwargs):
    """A path, or something already open.  Only a path gets closed.

    `Bio.File.as_handle` is four lines and this is the same four lines: try to
    open the argument, and if that is a `TypeError`, it was a handle all along.
    """
    try:
        with open(handleish, mode, **kwargs) as fp:
            yield fp
    except TypeError:
        yield handleish


class Array(np.ndarray):
    """numpy array subclass indexed by integers and by letters."""

    #   Stands in for PySequence_Check, which C has and Python does not.  A
    #   dict has a length and an item lookup and is still not a sequence.
    _NOT_A_SEQUENCE = (dict, set, frozenset)

    def __new__(cls, alphabet=None, dims=None, data=None, dtype=float):
        """Create a new Array instance.

        Three shapes of call, all of them the reference's: an array of zeros
        from an alphabet and a rank, an array from a numpy array, and an array
        from a codon-style dictionary of keys to values.

        The dictionary path has a defect that is kept: when the keys are plain
        strings -- the only thing a one-dimensional dictionary can have -- the
        assignment reads `letter`, a name the two-dimensional branch above it
        would have left behind and which therefore does not exist.  So
        `Array(data={"A": 1})` is an `UnboundLocalError` here exactly as it is
        upstream, and there is a test that says so.  With tuple keys of length
        two it works, and that is the path this library uses.
        """
        if isinstance(data, dict):
            if alphabet is not None:
                raise ValueError("alphabet must be None if data is a dict")
            if dims is not None:
                raise ValueError("dims must be None if data is a dict")
            alphabet = []
            single_letters = True
            for key in data:
                if isinstance(key, str):
                    if dims is None:
                        dims = 1
                    elif dims != 1:
                        raise ValueError("inconsistent dimensions in data")
                    alphabet.append(key)
                elif isinstance(key, tuple):
                    if dims is None:
                        dims = len(key)
                    elif dims != len(key):
                        raise ValueError("inconsistent dimensions in data")
                    if dims == 1:
                        if not isinstance(key, str):
                            raise ValueError("expected string")
                        if len(key) > 1:
                            single_letters = False
                        alphabet.append(key)
                    elif dims == 2:
                        for letter in key:
                            if not isinstance(letter, str):
                                raise ValueError("expected string")
                            if len(letter) > 1:
                                single_letters = False
                            alphabet.append(letter)
                    else:
                        raise ValueError(
                            "data array should be 1- or 2- dimensional "
                            "(found %d dimensions) in key" % dims
                        )
            alphabet = sorted(set(alphabet))
            if single_letters:
                alphabet = "".join(alphabet)
            else:
                alphabet = tuple(alphabet)
            n = len(alphabet)
            if dims == 1:
                shape = (n,)
            elif dims == 2:
                shape = (n, n)
            else:  # dims is None
                raise ValueError("data is an empty dictionary")
            obj = super().__new__(cls, shape, dtype)
            if dims == 1:
                for i, key in enumerate(alphabet):
                    obj[i] = data.get(letter, 0.0)
            elif dims == 2:
                for i1, letter1 in enumerate(alphabet):
                    for i2, letter2 in enumerate(alphabet):
                        key = (letter1, letter2)
                        value = data.get(key, 0.0)
                        obj[i1, i2] = value
        else:
            if alphabet is None:
                alphabet = string.ascii_uppercase
            elif not (isinstance(alphabet, (str, tuple))):
                raise ValueError("alphabet should be a string or a tuple")
            n = len(alphabet)
            if data is None:
                if dims is None:
                    dims = 1
                elif dims not in (1, 2):
                    raise ValueError("dims should be 1 or 2 (found %s)" % dims)
                shape = (n,) * dims
            else:
                if dims is None:
                    shape = data.shape
                    dims = len(shape)
                    if dims == 1:
                        pass
                    elif dims == 2:
                        if shape[0] != shape[1]:
                            raise ValueError("data array is not square")
                    else:
                        raise ValueError(
                            "data array should be 1- or 2- dimensional "
                            "(found %d dimensions) " % dims
                        )
                else:
                    shape = (n,) * dims
                    if data.shape != shape:
                        raise ValueError(
                            "data shape has inconsistent shape (expected (%s), found (%s))"
                            % (shape, data.shape)
                        )
            obj = super().__new__(cls, shape, dtype)
            if data is None:
                obj[:] = 0.0
            else:
                obj[:] = data
        obj.alphabet = alphabet
        return obj

    def __array_finalize__(self, obj):
        try:
            alphabet = obj.alphabet
        except AttributeError:
            # None, or plain numpy array
            pass
        else:
            if alphabet is not None:
                #   This goes through the setter, which is the point: a view
                #   whose shape no longer fits the alphabet is refused here,
                #   before any method gets to use it.
                self.alphabet = alphabet

    @property
    def alphabet(self):
        """The letters the indices stand for, or `None` if never set.

        `None` is a real answer and not an error: an `Array` made by viewing a
        plain numpy array has no alphabet, and asking for one is how a caller
        finds that out.
        """
        return getattr(self, "_alphabet", None)

    @alphabet.setter
    def alphabet(self, arg):
        """Set the alphabet once, if it fits the array.

        The three checks are the C setter's, in its order: it has not been set
        before, it is a sequence, and its length matches the array -- which for
        a 2-dimensional array may also be a single row or a single column, the
        shapes `select` produces on its way to somewhere else.  A string
        alphabet is additionally checked for a repeated letter, because the
        reference caches a character-to-index mapping for it and two indices
        for one character would be a lie.
        """
        if self.alphabet is not None:
            raise ValueError("the alphabet has already been set.")
        if isinstance(arg, self._NOT_A_SEQUENCE) or not (
            hasattr(arg, "__len__") and hasattr(arg, "__getitem__")
        ):
            raise TypeError(
                "alphabet must support the sequence protocol (e.g.,\n"
                "strings, lists, and tuples can be valid alphabets)."
            )
        length = len(arg)
        shape = self.shape
        ndim = len(shape)
        if ndim == 1:
            if shape[0] != length:
                raise ValueError(
                    "alphabet length %d is inconsistent with array size %d"
                    % (length, shape[0])
                )
        elif ndim == 2:
            if not (
                (shape[0] == length and shape[1] == length)
                or (shape[0] == length and shape[1] == 1)
                or (shape[0] == 1 and shape[1] == length)
            ):
                raise ValueError(
                    "alphabet length %d is inconsistent with array size (%d, %d)"
                    % (length, shape[0], shape[1])
                )
        else:
            raise ValueError(
                "substitution matrix has incorrect rank %d (expected 1 or 2)" % ndim
            )
        if isinstance(arg, str):
            seen = set()
            for character in arg:
                if character in seen:
                    raise ValueError(
                        "alphabet contains %r more than once" % character
                    )
                seen.add(character)
        self._alphabet = arg

    def _convert_key(self, key):
        if isinstance(key, tuple):
            indices = []
            for index in key:
                if isinstance(index, str):
                    try:
                        index = self.alphabet.index(index)
                    except ValueError:
                        raise IndexError("'%s'" % index) from None
                indices.append(index)
            key = tuple(indices)
        elif isinstance(key, str):
            try:
                key = self.alphabet.index(key)
            except ValueError:
                raise IndexError("'%s'" % key) from None
        return key

    def __getitem__(self, key):
        key = self._convert_key(key)
        value = np.ndarray.__getitem__(self, key)
        if value.ndim == 2:
            if self.ndim == 2:
                if value.shape != self.shape:
                    raise IndexError("Requesting truncated array")
            elif self.ndim == 1:
                length = self.shape[0]
                if value.shape[0] == length and value.shape[1] == 1:
                    pass
                elif value.shape[0] == 1 and value.shape[1] == length:
                    pass
                else:
                    raise IndexError("Requesting truncated array")
        elif value.ndim == 1:
            if value.shape[0] != self.shape[0]:
                raise IndexError("Requesting truncated array")
        elif value.ndim == 0:
            return value.item()
        return value.view(Array)

    def __setitem__(self, key, value):
        key = self._convert_key(key)
        np.ndarray.__setitem__(self, key, value)

    def __contains__(self, key):
        # Follow dict definition of __contains__
        return key in self.keys()

    def __array_prepare__(self, out_arr, context=None):
        # needed for numpy older than 1.13.0
        ufunc, inputs, i = context
        alphabet = self.alphabet
        for arg in inputs:
            if isinstance(arg, Array):
                if arg.alphabet != alphabet:
                    raise ValueError("alphabets are inconsistent")
        return np.ndarray.__array_prepare__(self, out_arr, context)

    def __array_wrap__(self, out_arr, context=None):
        if len(out_arr) == 1:
            return out_arr[0]
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        args = []
        alphabet = self.alphabet
        for arg in inputs:
            if isinstance(arg, Array):
                if arg.alphabet != alphabet:
                    raise ValueError("alphabets are inconsistent")
                args.append(arg.view(np.ndarray))
            else:
                args.append(arg)

        outputs = kwargs.pop("out", None)
        if outputs:
            out_args = []
            for arg in outputs:
                if isinstance(arg, Array):
                    if arg.alphabet != alphabet:
                        raise ValueError("alphabets are inconsistent")
                    out_args.append(arg.view(np.ndarray))
                else:
                    out_args.append(arg)
            kwargs["out"] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        raw_results = super().__array_ufunc__(ufunc, method, *args, **kwargs)
        if raw_results is NotImplemented:
            return NotImplemented

        if method == "at":
            return

        if ufunc.nout == 1:
            raw_results = (raw_results,)

        results = []
        for raw_result, output in zip(raw_results, outputs):
            if raw_result.ndim == 0:
                result = raw_result
            elif output is None:
                result = np.asarray(raw_result).view(Array)
                result.alphabet = self.alphabet
            else:
                result = output
            results.append(result)

        return results[0] if len(results) == 1 else results

    def __reduce__(self):
        import pickle

        values = np.array(self)
        state = pickle.dumps(values)
        alphabet = self.alphabet
        dims = len(self.shape)
        dtype = self.dtype
        arguments = (Array, alphabet, dims, None, dtype)
        return (Array.__new__, arguments, state)

    def __setstate__(self, state):
        import pickle

        self[:, :] = pickle.loads(state)

    def get(self, key, value=None):
        """Return the value of the key if found; return value otherwise."""
        try:
            return self[key]
        except IndexError:
            return value

    def items(self):
        """Return an iterator of (key, value) pairs in the array.

        Note that this walks a 2-dimensional array by rows while `keys` and
        `values` walk it by columns, so `list(matrix.items())` and
        `list(zip(matrix.keys(), matrix.values()))` are not the same list.
        That is the reference's arrangement and is kept.
        """
        dims = len(self.shape)
        if dims == 1:
            for index, key in enumerate(self.alphabet):
                value = np.ndarray.__getitem__(self, index)
                yield key, value
        elif dims == 2:
            for i1, c1 in enumerate(self.alphabet):
                for i2, c2 in enumerate(self.alphabet):
                    key = (c1, c2)
                    value = np.ndarray.__getitem__(self, (i1, i2))
                    yield key, value
        else:
            raise RuntimeError("array has unexpected shape %s" % self.shape)

    def keys(self):
        """Return a tuple with the keys associated with the array."""
        dims = len(self.shape)
        alphabet = self.alphabet
        if dims == 1:
            return tuple(alphabet)
        elif dims == 2:
            return tuple((c1, c2) for c2 in alphabet for c1 in alphabet)
        else:
            raise RuntimeError("array has unexpected shape %s" % self.shape)

    def values(self):
        """Return a tuple with the values stored in the array."""
        dims = len(self.shape)
        alphabet = self.alphabet
        if dims == 1:
            return tuple(self)
        elif dims == 2:
            n1, n2 = self.shape
            return tuple(
                np.ndarray.__getitem__(self, (i1, i2))
                for i2 in range(n2)
                for i1 in range(n1)
            )
        else:
            raise RuntimeError("array has unexpected shape %s" % self.shape)

    def update(self, E=None, **F):
        """Update the array from dict/iterable E and F."""
        if E is not None:
            try:
                alphabet = E.keys()
            except AttributeError:
                for key, value in E:
                    self[key] = value
            else:
                for key in E:
                    self[key] = E[key]
        for key in F:
            self[key] = F[key]

    def select(self, alphabet):
        """Subset the array by selecting the letters from the specified alphabet.

        Letters the array does not have are passed over in silence, so the
        result can be smaller than the alphabet asked for -- but not smaller in
        the shape it reports: a two-dimensional result is always square over the
        letters that were found.
        """
        ii = []
        jj = []
        for i, key in enumerate(alphabet):
            try:
                j = self.alphabet.index(key)
            except ValueError:
                continue
            ii.append(i)
            jj.append(j)
        dims = len(self.shape)
        a = Array(alphabet, dims=dims)
        ii = np.ix_(*[ii] * dims)
        jj = np.ix_(*[jj] * dims)
        a[ii] = self.view(np.ndarray)[jj]
        return a

    def _format_1D(self, fmt):
        alphabet = self.alphabet
        n = len(alphabet)
        words = [None] * n
        lines = []
        try:
            header = self.header
        except AttributeError:
            pass
        else:
            for line in header:
                line = "#  %s\n" % line
                lines.append(line)
        maxwidth = 0
        for i, key in enumerate(alphabet):
            value = self[key]
            word = fmt % value
            width = len(word)
            if width > maxwidth:
                maxwidth = width
            words[i] = word
        fmt2 = " %" + str(maxwidth) + "s"
        for letter, word in zip(alphabet, words):
            word = fmt2 % word
            line = letter + word + "\n"
            lines.append(line)
        text = "".join(lines)
        return text

    def _format_2D(self, fmt):
        alphabet = self.alphabet
        n = len(alphabet)
        words = [[None] * n for _ in range(n)]
        lines = []
        try:
            header = self.header
        except AttributeError:
            pass
        else:
            for line in header:
                line = "#  %s\n" % line
                lines.append(line)
        keywidth = max(len(c) for c in alphabet)
        keyfmt = "%" + str(keywidth) + "s"
        line = " " * keywidth
        for j, c2 in enumerate(alphabet):
            maxwidth = 0
            for i, c1 in enumerate(alphabet):
                key = (c1, c2)
                value = self[key]
                word = fmt % value
                width = len(word)
                if width > maxwidth:
                    maxwidth = width
                words[i][j] = word
            fmt2 = " %" + str(maxwidth) + "s"
            word = fmt2 % c2
            line += word
            for i, c1 in enumerate(alphabet):
                word = words[i][j]
                words[i][j] = fmt2 % word
        line = line.rstrip() + "\n"
        lines.append(line)
        for letter, row in zip(alphabet, words):
            key = keyfmt % letter
            line = key + "".join(row) + "\n"
            lines.append(line)
        text = "".join(lines)
        return text

    def __format__(self, fmt):
        return self.format(fmt)

    def format(self, fmt=""):
        """Return a string representation of the array.

        The argument ``fmt`` specifies the number format to be used.
        By default, the number format is "%i" if the array contains integer
        numbers, and "%.1f" otherwise.
        """
        if fmt == "":
            if np.issubdtype(self.dtype, np.integer):
                fmt = "%i"
            else:
                fmt = "%.1f"
        n = len(self.shape)
        if n == 1:
            return self._format_1D(fmt)
        elif n == 2:
            return self._format_2D(fmt)
        else:
            raise RuntimeError("Array has unexpected rank %d" % n)

    def __str__(self):
        return self.format()

    def __repr__(self):
        text = np.ndarray.__repr__(self)
        alphabet = self.alphabet
        if isinstance(alphabet, str):
            assert text.endswith(")")
            text = text[:-1] + ",\n         alphabet='%s')" % self.alphabet
        return text


def read(handle, dtype=float):
    """Parse the file and return an Array object.

    A file is a header of `#` lines -- kept, and printed back by `format()` --
    followed either by a list of `letter value` pairs, or by a row of column
    labels and then a labelled square.  Which of the two it is, is decided by
    the shape of the first two data rows and nothing else, so a file whose first
    row happens to have two fields is read as one-dimensional whatever its
    header says.

    Kept from the reference: the square's row labels are checked against the
    column labels with an `assert`, which is a `SyntaxError`-free way of saying
    "this cannot happen" and which `python -O` removes; a file of nothing but
    comments walks off the end of the line list and raises `IndexError`; and a
    blank line inside a square produces an empty row that shifts every row after
    it by one letter.
    """
    with as_handle(handle) as fp:
        lines = fp.readlines()

    header = []
    for i, line in enumerate(lines):
        if not line.startswith("#"):
            break
        header.append(line[1:].strip())
    rows = [line.split() for line in lines[i:]]
    if len(rows[0]) == len(rows[1]) == 2:
        alphabet = [key for key, value in rows]
        for key in alphabet:
            if len(key) > 1:
                alphabet = tuple(alphabet)
                break
        else:
            alphabet = "".join(alphabet)
        matrix = Array(alphabet=alphabet, dims=1, dtype=dtype)
        matrix.update(rows)
    else:
        alphabet = rows.pop(0)
        for key in alphabet:
            if len(key) > 1:
                alphabet = tuple(alphabet)
                break
        else:
            alphabet = "".join(alphabet)
        matrix = Array(alphabet=alphabet, dims=2, dtype=dtype)
        for letter1, row in zip(alphabet, rows):
            letter = row.pop(0)
            assert letter1 == letter
            for letter2, word in zip(alphabet, row):
                matrix[letter1, letter2] = float(word)
    matrix.header = header
    return matrix


def build(alphabet, symmetric, values):
    """One `Array` from a row of the generated data file.

    `values` is the matrix's numbers as one string of whitespace-separated
    tokens, in the order the file's own numbers appeared: row by row if
    `symmetric` is false, and the lower triangle row by row -- `(0,0)`,
    `(1,0)`, `(1,1)`, `(2,0)`, ... -- if it is true, since the other half is then
    the same number and storing it twice would say nothing.

    The tokens are kept as text rather than turned into numbers by the
    generator, so that `float()` here does exactly the conversion the reference
    did when it read the same token out of the same file.  A value the generator
    wrote as a float literal would have to survive a round trip through `repr`
    to be equally exact, and there is no reason to take that chance.
    """
    n = len(alphabet)
    flat = np.array(values.split(), dtype=np.float64)
    if symmetric:
        expected = n * (n + 1) // 2
        if flat.size != expected:
            raise ValueError(
                "a symmetric %d-letter matrix needs %d numbers, found %d"
                % (n, expected, flat.size)
            )
        matrix = np.zeros((n, n), dtype=np.float64)
        lower = np.tril_indices(n)
        matrix[lower] = flat
        matrix[(lower[1], lower[0])] = flat
    else:
        if flat.size != n * n:
            raise ValueError(
                "a %d-letter matrix needs %d numbers, found %d"
                % (n, n * n, flat.size)
            )
        matrix = flat.reshape(n, n)
    obj = matrix.view(Array)
    obj.alphabet = alphabet
    return obj


def load_row(header, alphabet, symmetric, values):
    """The `Array` a row of the data file describes, header and all."""
    matrix = build(alphabet, symmetric, values)
    matrix.header = list(header)
    return matrix
