"""`Bio.SeqUtils.ProtParamData`: the amino-acid scales, as the reference exposes them.

The values live in :mod:`biofasting._protparam_data`, generated and verified
against Biopython 1.88 by `tools/gen_protparam_data.py`, aligned letter by letter
with the twenty standard residues.  This module is the *interface* to them, and
it exists because the reference's own interface is not a table: it is thirty-odd
module attributes, twenty-eight of which are shared objects.

Two facts about that interface are easy to get wrong and both matter to a
caller.

**The names are not one scale each.**  `gravy_scales["KyteDoolitle"]` is the very
same dict object as the module-level `kd` -- the reference builds `gravy_scales`
out of its own attributes -- so `gravy(scale="KyteDoolitle")` and a
`protein_scale(kd, ...)` are reading one table and not two copies of it.  That
identity is what lets :mod:`biofasting.protparam` recognise a caller's scale and
hand it to a kernel; a copy is still a valid scale, it just takes the
reference's own loop.  `Cowan3.4` and `Cowan7.5` are the exception that proves
the rule: the reference writes those two inline inside `gravy_scales` and never
binds them to a module attribute, so here they are reachable as `cw[3.4]` and
`cw[7.5]` and as `gravy_scales["Cowan3.4"]` / `['Cowan7.5']`, and not as
`Cowan3.4`.

**The letter order is `AMINO_ACIDS`, not the order the reference's literals were
written in.**  Nothing in this package iterates a scale -- `protein_scale` and
`gravy` look a scale up by the sequence's letters -- so which order the keys
happen to sit in is invisible in every answer.  It is nonetheless stated,
because the dicts here are the ones a kernel reads, and there the order is not
invisible at all.  `gravy_scales` itself is built in the reference's own order,
`KyteDoolitle` first, because that is a list a caller can compare.

>>> from biofasting.protparam_data import kd, gravy_scales
>>> gravy_scales["KyteDoolitle"] is kd
True
>>> round(kd["I"], 2)
4.5
"""

from __future__ import annotations

from ._protparam_data import AMINO_ACIDS, DIWV as _DIWV_ROWS, GRAVY_SCALES, SCALES

# The reference's own module attributes, plus `SCALES_BY_NAME`.  The gravy names
# are deliberately absent: they are keys of `gravy_scales`, not names in this
# module, and `import *` would fail on every one of them.
__all__ = [
    "DIWV",
    "Flex",
    "SCALES_BY_NAME",
    "ab",
    "ag",
    "al",
    "bb",
    "bm",
    "ci",
    "cs",
    "cw",
    "eg",
    "em",
    "es",
    "fc",
    "fs",
    "gd",
    "gravy_scales",
    "gy",
    "hw",
    "ja",
    "jo",
    "ju",
    "kd",
    "ki",
    "mi",
    "pa",
    "po",
    "rm",
    "ro",
    "sw",
    "ta",
    "wi",
    "zi",
]

# `dict(zip(...))` rather than a literal per scale: the tables are already in
# the generated module, in the one order a kernel can index them by, and a
# second literal here would be a second place for a digit to be mistyped.


def _scale(name):
    """One twenty-letter scale, as the reference's ``{letter: float}``."""
    return dict(zip(AMINO_ACIDS, SCALES[name]))


Flex = _scale("Flex")
ab = _scale("ab")
ag = _scale("ag")
al = _scale("al")
bb = _scale("bb")
bm = _scale("bm")
ci = _scale("ci")
cs = _scale("cs")
eg = _scale("eg")
em = _scale("em")
es = _scale("es")
fc = _scale("fc")
fs = _scale("fs")
gd = _scale("gd")
gy = _scale("gy")
hw = _scale("hw")
ja = _scale("ja")
jo = _scale("jo")
ju = _scale("ju")
kd = _scale("kd")
ki = _scale("ki")
mi = _scale("mi")
pa = _scale("pa")
po = _scale("po")
rm = _scale("rm")
ro = _scale("ro")
sw = _scale("sw")
ta = _scale("ta")
wi = _scale("wi")
zi = _scale("zi")

#: The two Cowan scales, keyed by the pH they are tabulated at, which is how
#: the reference keys them.  They have no module attribute of their own.
cw = {3.4: _scale("Cowan3.4"), 7.5: _scale("Cowan7.5")}

#: Every scale in this module, by name, the two inside `cw` included.  This one
#: is ours and not the reference's: it is the list a caller walks to precompute
#: something per scale -- `biofasting.protparam` builds its dense kernel tables
#: from it -- and the list a test sweeps.  A name here is not a module
#: attribute; `Cowan3.4` is `cw[3.4]`.
SCALES_BY_NAME = {
    "Flex": Flex,
    "ab": ab,
    "ag": ag,
    "al": al,
    "bb": bb,
    "bm": bm,
    "ci": ci,
    "cs": cs,
    "eg": eg,
    "em": em,
    "es": es,
    "fc": fc,
    "fs": fs,
    "gd": gd,
    "gy": gy,
    "hw": hw,
    "ja": ja,
    "jo": jo,
    "ju": ju,
    "kd": kd,
    "ki": ki,
    "mi": mi,
    "pa": pa,
    "po": po,
    "rm": rm,
    "ro": ro,
    "sw": sw,
    "ta": ta,
    "wi": wi,
    "zi": zi,
    "Cowan3.4": cw[3.4],
    "Cowan7.5": cw[7.5],
}

#: The scales `ProteinAnalysis.gravy` accepts, by the name it accepts them
#: under.  Built in the reference's own order -- `KyteDoolitle` first, then the
#: rest alphabetically -- so that `list(gravy_scales)` is comparable with
#: Biopython's.
gravy_scales = {
    "KyteDoolitle": kd,
    "Aboderin": ab,
    "AbrahamLeo": al,
    "Argos": ag,
    "BlackMould": bm,
    "BullBreese": bb,
    "Casari": cs,
    "Cid": ci,
    "Cowan3.4": cw[3.4],
    "Cowan7.5": cw[7.5],
    "Eisenberg": es,
    "Engelman": eg,
    "Fasman": fs,
    "Fauchere": fc,
    "GoldSack": gd,
    "Guy": gy,
    "Jones": jo,
    "Juretic": ju,
    "Kidera": ki,
    "Miyazawa": mi,
    "Parker": pa,
    "Ponnuswamy": po,
    "Rose": ro,
    "Roseman": rm,
    "Sweet": sw,
    "Tanford": ta,
    "Wilson": wi,
    "Zimmerman": zi,
}

#: Guruprasad's dipeptide instability weights, as `DIWV[this][next]`.  Nested
#: dicts rather than a table, because that is what the reference indexes and the
#: one thing that has to match here is which residue a `KeyError` names.
DIWV = {
    this: dict(zip(AMINO_ACIDS, row))
    for this, row in zip(AMINO_ACIDS, _DIWV_ROWS)
}
