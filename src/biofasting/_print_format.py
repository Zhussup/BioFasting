"""Turn a batch's search results into the three reports a person reads.

`Bio.Restriction.PrintFormat` is the printing half of `Analysis`: it takes the
dictionary `RestrictionBatch.search()` returns and lays it out as a list, a
count per enzyme, or a map of the sequence with the cuts drawn on it.  It is a
port of that module, line for line, including the parts that are odd -- see
"Deliberately kept" below -- because a report that differs from the reference's
by a space is a report somebody will diff.

The class is usable on its own, exactly as upstream's is:

    >>> from biofasting.restriction import PrintFormat, RestrictionBatch
    >>> batch = RestrictionBatch(["EcoRI", "BamHI", "ApaI"])
    >>> report = PrintFormat()
    >>> report.sequence = "GGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTC"
    >>> report.print_that(batch.search(report.sequence), "My map:\\n", "No site:\\n")
    ... # doctest: +NORMALIZE_WHITESPACE
    My map:
    ApaI       :  12.
    EcoRI      :  50.
    No site:
    BamHI
    <BLANKLINE>

Each name in the non-cutting list is padded to `NameWidth`, so the line above
is really `"BamHI" + " " * 5`.  The example normalises whitespace instead of
carrying that padding in the source, where it would be invisible, and where the
first tool that strips trailing blanks would turn a true example into a false
one.  The padding itself is asserted where it can be seen: see the `Nothing
here:` case in `tests/test_restriction.py`.

Deliberately kept, because they are the reference's behaviour and not a
mistake to be quietly improved on:

- `change()` and the class attributes recompute `Cmodulo` and `PrefWidth` but
  never `linesize`, which was computed once at class-definition time.  Setting
  `ConsoleWidth` therefore changes the wrap of the *list* section through
  `__next_section` (which subtracts `MaxSize` from `linesize`) but not through
  `_make_nocut_only` (which uses `linesize` bare).  That asymmetry is upstream's.
- `_make_nocut_only` appends the trailing newline unconditionally, so a report
  whose `s1` is empty and whose non-cutting list is empty still ends in "\n".
- `_make_map_only` indexes with `key % 60` and slices to `k - 1`, so a cut at a
  position that is a multiple of 60 is drawn one column left of where a
  zero-based reading would put it.  It is what the reference prints.

Not kept: the reference's `_make_map_only` also computes
`resultKeys = sorted(...)` and never reads it.  Dead code with no observable
effect is the one thing a port may drop, and dropping it is not a difference a
test could find.
"""

from __future__ import annotations

import re

#   What `Bio.Seq.Seq.complement()` does, so that a map drawn over a `str`
#   matches one drawn over a `Seq`: the sixteen IUPAC codes complement to their
#   partner in both cases -- `U` included, which Biopython complements to `A`,
#   since a uracil pairs with an adenine -- and every other byte is left alone.
#   Bytes rather than `str` for the same reason as everywhere else here: it is
#   exact for all 256 values and needs no decoding.
_PAIRS = {
    "A": "T", "T": "A", "C": "G", "G": "C", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
_COMPLEMENT = bytes(
    ord(_PAIRS.get(base.upper(), base).lower() if base.islower() else _PAIRS.get(base, base))
    if byte < 128
    else byte
    for byte, base in ((byte, chr(byte)) for byte in range(256))
)


def _complement(sequence) -> str:
    """The sequence's other strand, 3' to 5', as the map prints it."""
    if hasattr(sequence, "complement"):
        return str(sequence.complement())
    if isinstance(sequence, (bytes, bytearray, memoryview)):
        return bytes(sequence).translate(_COMPLEMENT).decode("ascii")
    return str(sequence).encode("ascii").translate(_COMPLEMENT).decode("ascii")


class PrintFormat:
    """Format a restriction analysis as a list, a count, or a map.

    Subclasses override `make_format`, or point it at one of the `_make_*`
    methods with `print_as`; the reference documents this as the extension
    point and so does this.
    """

    ConsoleWidth = 80
    NameWidth = 10
    MaxSize = 6
    Cmodulo = ConsoleWidth % NameWidth
    PrefWidth = ConsoleWidth - Cmodulo
    Indent = 4
    linesize = PrefWidth - NameWidth

    def print_as(self, what: str = "list") -> None:
        """Choose the report: `list`, `number`, or `map`."""
        if what == "map":
            self.make_format = self._make_map
        elif what == "number":
            self.make_format = self._make_number
        else:
            self.make_format = self._make_list

    def format_output(self, dct, title: str = "", s1: str = "") -> str:
        """The report for `dct`, split into the cutters and the non-cutters.

        `title` heads the report and must carry its own line break; `s1` is the
        sentence before the list of enzymes that do not cut.
        """
        if not dct:
            dct = self.results
        ls, nc = [], []
        for key, value in dct.items():
            if value:
                ls.append((key, value))
            else:
                nc.append(key)
        return self.make_format(ls, title, nc, s1)

    def print_that(self, dct, title: str = "", s1: str = "") -> None:
        """Print `format_output(dct, title, s1)`.  Kept for compatibility."""
        print(self.format_output(dct, title, s1))

    def make_format(self, cut=(), title: str = "", nc=(), s1: str = "") -> str:
        """The virtual method `format_output` calls.  Defaults to the list."""
        return self._make_list(cut, title, nc, s1)

    # ----------------------------------------------------------------- #
    #   The three reports
    # ----------------------------------------------------------------- #

    def _make_list(self, ls, title, nc, s1) -> str:
        return self._make_list_only(ls, title) + self._make_nocut_only(nc, s1)

    def _make_map(self, ls, title, nc, s1) -> str:
        return self._make_map_only(ls, title) + self._make_nocut_only(nc, s1)

    def _make_number(self, ls, title, nc, s1) -> str:
        return self._make_number_only(ls, title) + self._make_nocut_only(nc, s1)

    def _make_nocut(self, ls, title, nc, s1) -> str:
        return title + self._make_nocut_only(nc, s1)

    def _make_nocut_only(self, nc, s1, ls=(), title="") -> str:
        """The enzymes that do not cut, wrapped at `linesize`, as a string."""
        if not nc:
            return s1
        st = ""
        stringsite = s1 or "\n   Enzymes which do not cut the sequence.\n\n"
        for key in sorted(nc):
            st = "".join((st, str.ljust(str(key), self.NameWidth)))
            if len(st) > self.linesize:
                stringsite = "".join((stringsite, st, "\n"))
                st = ""
        stringsite = "".join((stringsite, st, "\n"))
        return stringsite

    def _make_list_only(self, ls, title, nc=(), s1="") -> str:
        """`enzyme : position, position.` one per cutting enzyme."""
        if not ls:
            return title
        return self.__next_section(ls, title)

    def _make_number_only(self, ls, title, nc=(), s1="") -> str:
        """The same lines, grouped by how many times the enzyme cuts."""
        if not ls:
            return title
        ls.sort(key=lambda x: len(x[1]))
        iterator = iter(ls)
        cur_len = 1
        new_sect = []
        for name, sites in iterator:
            length = len(sites)
            if length > cur_len:
                title += "\n\nenzymes which cut %i times :\n\n" % cur_len
                title = self.__next_section(new_sect, title)
                new_sect, cur_len = [(name, sites)], length
                continue
            new_sect.append((name, sites))
        title += "\n\nenzymes which cut %i times :\n\n" % cur_len
        return self.__next_section(new_sect, title)

    def _make_map_only(self, ls, title, nc=(), s1="") -> str:
        """The sequence with its cuts drawn above each 60-column block.

        Three lines per block: the positions with their enzymes, a column of
        `|` marking the cuts, then the top strand, a rule, the bottom strand,
        and the block's coordinates.
        """
        if not ls:
            return title
        out = title or ""
        enzymemap = {}
        for enzyme, cut in ls:
            for position in cut:
                if position in enzymemap:
                    enzymemap[position].append(str(enzyme))
                else:
                    enzymemap[position] = [str(enzyme)]
        #   Every cut position that is still unplaced, in order.  As each block
        #   is printed the positions it covered are moved out of this list and
        #   into `cutloc`, so the two together are a partition of the cuts.
        mapping = sorted(enzymemap.keys())
        cutloc = {}
        x, counter, length = 0, 0, len(self.sequence)
        for x in range(60, length, 60):
            counter = x - 60
            location = []
            cutloc[counter] = location
            remaining = []
            for key in mapping:
                if key <= x:
                    location.append(key)
                else:
                    remaining.append(key)
            mapping = remaining
        cutloc[x] = mapping
        sequence = str(self.sequence)
        revsequence = _complement(self.sequence)
        bar = "|"
        base, counter = 0, 0
        emptyline = " " * 60
        for base in range(60, length, 60):
            counter = base - 60
            line = emptyline
            for key in cutloc[counter]:
                s = ""
                if key == base:
                    for name in enzymemap[key]:
                        s = " ".join((s, name))
                    chunk = line[0:59]
                    lineo = "".join((chunk, str(key), s, "\n"))
                    line2 = "".join((chunk, bar, "\n"))
                    out = "".join((out, lineo, line2))
                    break
                for name in enzymemap[key]:
                    s = " ".join((s, name))
                k = key % 60
                lineo = "".join((line[0 : (k - 1)], str(key), s, "\n"))
                line = "".join((line[0 : (k - 1)], bar, line[k:]))
                line2 = "".join((line[0 : (k - 1)], bar, line[k:], "\n"))
                out = "".join((out, lineo, line2))
            mapunit = "\n".join(
                (
                    sequence[counter:base],
                    bar * 60,
                    revsequence[counter:base],
                    "".join(
                        (
                            str.ljust(str(counter + 1), 15),
                            " " * 30,
                            str.rjust(str(base), 15),
                            "\n\n",
                        )
                    ),
                )
            )
            out = "".join((out, mapunit))
        line = " " * 60
        for key in cutloc[base]:
            s = ""
            if key == length:
                for name in enzymemap[key]:
                    s = "".join((s, " ", name))
                chunk = line[0 : (length - 1)]
                lineo = "".join((chunk, str(key), s, "\n"))
                line2 = "".join((chunk, bar, "\n"))
                out = "".join((out, lineo, line2))
                break
            for name in enzymemap[key]:
                s = "".join((s, " ", name))
            k = key % 60
            lineo = "".join((line[0 : (k - 1)], str(key), s, "\n"))
            line = "".join((line[0 : (k - 1)], bar, line[k:]))
            line2 = "".join((line[0 : (k - 1)], bar, line[k:], "\n"))
            out = "".join((out, lineo, line2))
        mapunit = "".join((sequence[base:length], "\n"))
        mapunit = "".join((mapunit, bar * (length - base), "\n"))
        mapunit = "".join((mapunit, revsequence[base:length], "\n"))
        mapunit = "".join(
            (
                mapunit,
                "".join(
                    (
                        str.ljust(str(base + 1), 15),
                        " " * (length - base - 30),
                        str.rjust(str(length), 15),
                        "\n\n",
                    )
                ),
            )
        )
        return "".join((out, mapunit))

    def __next_section(self, ls, into: str) -> str:
        """`enzyme     :  position, position.` lines, wrapped with an indent.

        Wrapped by `re.finditer` over a pattern that stops at a comma or a
        full stop: the reference cuts the line where it can rather than where
        it should, which is why a long list of positions ends with a line that
        may run past `linesize`.
        """
        indentation = "\n" + (self.NameWidth + self.Indent) * " "
        linesize = self.linesize - self.MaxSize
        pat = re.compile(r"([\w,\s()]){1,%i}[,\.]" % linesize)
        for name, sites in sorted(ls):
            output = ", ".join(str(site) for site in sites) + "."
            if len(output) > linesize:
                output = [match.group() for match in re.finditer(pat, output)]
                stringsite = indentation.join(output)
            else:
                stringsite = output
            into = "".join(
                (into, str(name).ljust(self.NameWidth), " :  ", stringsite, "\n")
            )
        return into
