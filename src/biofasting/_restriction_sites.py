"""Derive a restriction enzyme's search pattern from its recognition site.

The data table (`_restriction_data.py`) stores what Rebase states: a
recognition site, cut coordinates, an overhang, a frequency.  It does not store
a regular expression, because the expression is a *function of the site* -- and
storing both would be two sources of truth for one fact, free to disagree.

This module is that function, and it has two callers on purpose:
`biofasting.restriction` compiles the pattern it returns, and
`tools/gen_restriction_data.py` uses the same function to check itself against
the reference at generation time.  The tool loads this file *by path* rather
than importing the package, so regenerating the data does not require a
compiled extension -- which is the only reason the derivation lives in its own
module instead of inside `restriction.py`.

That check is the load-bearing one: the generator refuses to write the data
table unless `compsite()` reproduces Biopython's own compiled pattern for every
one of the 1,088 enzymes.  A derivation that quietly disagreed on three enzymes
would ship a library that is wrong on three enzymes forever, and nothing would
notice, because a data file has no opinion about itself.
"""

from __future__ import annotations

# The IUPAC expansion the reference's patterns use.  `N` is `.` rather than
# `[ACGT]`, which matters: a site like ACNNNNNCTCC must match across a run of
# any five bases, ambiguous codes included, and Biopython writes it that way.
IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "W": "[AT]",
    "S": "[GC]",
    "M": "[AC]",
    "K": "[GT]",
    "H": "[ACT]",
    "B": "[CGT]",
    "V": "[ACG]",
    "D": "[AGT]",
    "N": ".",
}

COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "W": "W", "S": "S", "M": "K", "K": "M",
    "H": "D", "B": "V", "V": "B", "D": "H", "N": "N",
}


def _watson(site: str) -> str:
    """The part of a site that describes the top strand.

    Exactly one enzyme in the table -- `HpyUM037X`, which is uncharacterised --
    carries a site of the form ``TNGGNAG|GTGGNAG``, where the pipe separates
    the two strands' recognition sequences.  The reference's own compiled
    pattern uses the part before the pipe and ignores the rest, so this does the
    same.  It is written down here rather than filtered silently, because a
    `split` that nobody explains looks like a bug fix and is actually the
    reference's behaviour.
    """
    return site.split("|", 1)[0]


def reverse_complement(site: str) -> str:
    """The other strand's recognition sequence."""
    return "".join(COMPLEMENT[base] for base in reversed(_watson(site)))


def pattern_for(site: str) -> str:
    """The regex body that matches a recognition site."""
    return "".join(IUPAC[base] for base in _watson(site))


def compsite(site: str, palindromic: bool, top: str = "t", bottom: str = "b") -> str:
    """The search pattern for a site, as a regex.

    Biopython's form is a lookahead so that matches may overlap -- two sites one
    base apart are two sites -- with a named group per strand, which is how a
    non-palindromic enzyme tells which strand it matched.  The group names are
    parameters here: the reference names them after the enzyme, while the engine
    uses one short pair for every enzyme.
    """
    top_group = f"(?=(?P<{top}>{pattern_for(site)}))"
    if palindromic:
        return top_group
    return f"{top_group}|(?=(?P<{bottom}>{pattern_for(reverse_complement(site))}))"
