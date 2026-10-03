"""The amino-acid scales `biofasting.protparam` is built on.

**Generated** by `tools/gen_protparam_data.py` -- do not edit.  32 tables
read from `Bio.SeqUtils.ProtParamData` and `Bio.SeqUtils.IsoelectricPoint`, with
every letter aligned to `AMINO_ACIDS` instead of to the order a literal happened
to be written in, because a kernel indexes them by position and the reference
indexes them by letter.

`FLEX` is `ProtParamData.Flex`, the twenty values `ProteinAnalysis.flexibility`
uses.  `DIWV` is `ProtParamData.DIWV`, the 400 Guruprasad dipeptide instability
weights, row-major in `AMINO_ACIDS`.  `SCALES` is every twenty-letter scale the
reference exposes -- the twenty-eight `gravy_scales` entries among them -- and
`GRAVY_SCALES` maps a `gravy(scale=...)` name onto the key it names.  The three
constants after them are `IsoelectricPoint`'s, and their *order* is load-bearing:
`charge_at_pH` sums over them in insertion order and floating-point addition is
not associative.

**An `int` here stays an `int`.**  `sum(selected_scale[aa] for aa in sequence)` is
what `gravy` is, and CPython since 3.12 compensates a float item and adds an
`int` one plainly, with the compensation left where it was.  Ten of the
twenty-eight gravy scales are written with integer values -- Engelman's C, G and
H, Parker's N, D, Q, I and W, GoldSack's seven, and one or two each in Fasman,
Fauchere, Jones, Ponnuswamy, Roseman, Wilson and Zimmerman -- and on four of
those ten a compensated sum over the same twenty numbers answers differently,
between 4% and 29% of random proteins.  Rendering a `2` as `2.0` therefore changes
the last bit of `gravy`.  `repr` is what keeps them apart; nothing else in this
file should ever be reformatted by hand.

Read by `biofasting.protparam`; the generator that wrote it is the only thing
that should ever rewrite it.
"""

from __future__ import annotations

#: The Biopython release these were read from: 1.88
PROTPARAM_SOURCE_VERSION = '1.88'

#: The twenty standard residues, which is `IUPACData.protein_letters`.  Every
#: scale below is a tuple in this order; it is the order a dense kernel table
#: needs, and it is not the order the reference's literals are written in.
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

#: `ProtParamData.Flex`, in `AMINO_ACIDS` order.
FLEX: tuple = (
    0.984, 0.906, 1.068, 1.094, 0.915, 1.031, 0.95, 0.927, 1.102, 0.935, 0.952, 1.048,
    1.049, 1.037, 1.008, 1.046, 0.997, 0.931, 0.904, 0.929
)

#: `ProtParamData.DIWV`, row-major in `AMINO_ACIDS`: `DIWV[i][j]`
#: is the instability weight of the dipeptide `AMINO_ACIDS[i]`
#: followed by `AMINO_ACIDS[j]`.
DIWV: tuple = (
    (
        1.0, 44.94, -7.49, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    ),
    (
        1.0, 1.0, 20.26, 1.0, 1.0, 1.0, 33.6, 1.0, 1.0, 20.26, 33.6, 1.0, 20.26, -6.54,
        1.0, 1.0, 33.6, -6.54, 24.68, 1.0
    ),
    (
        1.0, 1.0, 1.0, 1.0, -6.54, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 1.0, 1.0,
        -6.54, 20.26, -14.03, 1.0, 1.0, 1.0
    ),
    (
        1.0, 44.94, 20.26, 33.6, 1.0, 1.0, -6.54, 20.26, 1.0, 1.0, 1.0, 1.0, 20.26,
        20.26, 1.0, 20.26, 1.0, 1.0, -14.03, 1.0
    ),
    (
        1.0, 1.0, 13.34, 1.0, 1.0, 1.0, 1.0, 1.0, -14.03, 1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 33.601
    ),
    (
        -7.49, 1.0, 1.0, -6.54, 1.0, 13.34, 1.0, -7.49, -7.49, 1.0, 1.0, -7.49, 1.0,
        1.0, 1.0, 1.0, -7.49, 1.0, 13.34, -7.49
    ),
    (
        1.0, 1.0, 1.0, 1.0, -9.37, -9.37, 1.0, 44.94, 24.68, 1.0, 1.0, 24.68, -1.88,
        1.0, 1.0, 1.0, -6.54, 1.0, -1.88, 44.94
    ),
    (
        1.0, 1.0, 1.0, 44.94, 1.0, 1.0, 13.34, 1.0, -7.49, 20.26, 1.0, 1.0, -1.88, 1.0,
        1.0, 1.0, 1.0, -7.49, 1.0, 1.0
    ),
    (
        1.0, 1.0, 1.0, 1.0, 1.0, -7.49, 1.0, -7.49, 1.0, -7.49, 33.6, 1.0, -6.54,
        24.64, 33.6, 1.0, 1.0, -7.49, 1.0, 1.0
    ),
    (
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -7.49, 1.0, 1.0, 1.0, 20.26, 33.6,
        20.26, 1.0, 1.0, 1.0, 24.68, 1.0
    ),
    (
        13.34, 1.0, 1.0, 1.0, 1.0, 1.0, 58.28, 1.0, 1.0, 1.0, -1.88, 1.0, 44.94, -6.54,
        -6.54, 44.94, -1.88, 1.0, 1.0, 24.68
    ),
    (
        1.0, -1.88, 1.0, 1.0, -14.03, -14.03, 1.0, 44.94, 24.68, 1.0, 1.0, 1.0, -1.88,
        -6.54, 1.0, 1.0, -7.49, 1.0, -9.37, 1.0
    ),
    (
        20.26, -6.54, -6.54, 18.38, 20.26, 1.0, 1.0, 1.0, 1.0, 1.0, -6.54, 1.0, 20.26,
        20.26, -6.54, 20.26, 1.0, 20.26, -1.88, 1.0
    ),
    (
        1.0, -6.54, 20.26, 20.26, -6.54, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 20.26,
        20.26, 1.0, 44.94, 1.0, -6.54, 1.0, -6.54
    ),
    (
        1.0, 1.0, 1.0, 1.0, 1.0, -7.49, 20.26, 1.0, 1.0, 1.0, 1.0, 13.34, 20.26, 20.26,
        58.28, 44.94, 1.0, 1.0, 58.28, -6.54
    ),
    (
        1.0, 33.6, 1.0, 20.26, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 44.94, 20.26,
        20.26, 20.26, 1.0, 1.0, 1.0, 1.0
    ),
    (
        1.0, 1.0, 1.0, 20.26, 13.34, -7.49, 1.0, 1.0, 1.0, 1.0, 1.0, -14.03, 1.0,
        -6.54, 1.0, 1.0, 1.0, 1.0, -14.03, 1.0
    ),
    (
        1.0, 1.0, -14.03, 1.0, 1.0, -7.49, 1.0, 1.0, -1.88, 1.0, 1.0, 1.0, 20.26, 1.0,
        1.0, 1.0, -7.49, 1.0, 1.0, -6.54
    ),
    (
        -14.03, 1.0, 1.0, 1.0, 1.0, -9.37, 24.68, 1.0, 1.0, 13.34, 24.68, 13.34, 1.0,
        1.0, 1.0, 1.0, -14.03, -7.49, 1.0, 1.0
    ),
    (
        24.68, 1.0, 24.68, -6.54, 1.0, -7.49, 13.34, 1.0, 1.0, 1.0, 44.94, 1.0, 13.34,
        1.0, -15.91, 1.0, -7.49, 1.0, -9.37, 13.34
    ),
)

#: Every twenty-letter scale the reference exposes, by the name it
#: exposes it under.  `gravy(scale=...)` reads the twenty-eight in
#: `GRAVY_SCALES`; `protein_scale(param_dict, ...)` takes any of them.
SCALES: dict = {
    'Cowan3.4': (
        0.42, 0.84, -0.51, -0.37, 1.74, 0.0, -2.28, 1.81, -2.03, 1.8, 1.18, -1.03,
        0.86, -0.96, -1.56, -0.64, -0.26, 1.34, 1.46, 0.51
    ),
    'Cowan7.5': (
        0.35, 0.76, -2.15, -1.95, 1.69, 0.0, -0.65, 1.83, -1.54, 1.8, 1.1, -0.99, 0.84,
        -0.93, -1.5, -0.63, -0.27, 1.32, 1.35, 0.39
    ),
    'Flex': (
        0.984, 0.906, 1.068, 1.094, 0.915, 1.031, 0.95, 0.927, 1.102, 0.935, 0.952,
        1.048, 1.049, 1.037, 1.008, 1.046, 0.997, 0.931, 0.904, 0.929
    ),
    'ab': (
        5.1, 0.0, 0.7, 1.8, 9.6, 4.1, 1.6, 9.3, 1.3, 10.0, 8.7, 0.6, 4.9, 1.4, 2.0,
        3.1, 3.5, 8.5, 9.2, 8.0
    ),
    'ag': (
        0.61, 1.07, 0.46, 0.47, 2.02, 0.07, 0.61, 2.22, 1.15, 1.53, 1.18, 0.06, 1.95,
        0.0, 0.6, 0.05, 0.05, 1.32, 2.65, 1.88
    ),
    'al': (
        0.44, 0.58, -0.31, -0.34, 2.54, 0.0, -0.01, 2.46, -2.45, 2.46, 1.1, -1.32,
        1.29, -0.71, -2.42, -0.84, -0.41, 1.73, 2.56, 1.63
    ),
    'bb': (
        0.61, 0.36, 0.61, 0.51, -1.52, 0.81, 0.69, -1.45, 0.46, -1.65, -0.66, 0.89,
        -0.17, 0.97, 0.69, 0.42, 0.29, -0.75, -1.2, -1.43
    ),
    'bm': (
        0.616, 0.68, 0.028, 0.043, 1.0, 0.501, 0.165, 0.943, 0.283, 0.943, 0.738,
        0.236, 0.711, 0.251, 0.0, 0.359, 0.45, 0.825, 0.878, 0.88
    ),
    'ci': (
        0.02, 0.77, -1.04, -1.14, 1.35, -0.8, 0.26, 1.81, -0.41, 1.14, 1.0, -0.77,
        -0.09, -1.1, -0.42, -0.97, -0.77, 1.13, 1.71, 1.11
    ),
    'cs': (
        0.2, 1.9, -1.4, -1.3, 1.0, -0.1, 0.4, 1.4, -1.6, 0.5, 0.5, -0.5, -1.0, -1.1,
        -0.7, -0.7, -0.4, 0.7, 1.6, 0.5
    ),
    'eg': (
        -1.6, -2, 9.2, 8.2, -3.7, -1, 3, -3.1, 8.8, -2.8, -3.4, 4.8, 0.2, 4.1, 12.3,
        -0.6, -1.2, -2.6, -1.9, 0.7
    ),
    'em': (
        0.815, 0.394, 1.283, 1.445, 0.695, 0.714, 1.18, 0.603, 1.545, 0.603, 0.714,
        1.296, 1.236, 1.348, 1.475, 1.115, 1.184, 0.606, 0.808, 1.089
    ),
    'es': (
        0.62, 0.29, -0.9, -0.74, 1.19, 0.48, -0.4, 1.38, -1.5, 1.06, 0.64, -0.78, 0.12,
        -0.85, -2.53, -0.18, -0.05, 1.08, 0.81, 0.26
    ),
    'fc': (
        0.31, 1.54, -0.77, -0.64, 1.79, 0, 0.13, 1.8, -0.99, 1.7, 1.23, -0.6, 0.72,
        -0.22, -1.01, -0.04, 0.26, 1.22, 2.25, 0.96
    ),
    'fs': (
        -0.21, -6.04, 1.36, 2.3, -4.65, 0, -1.23, -4.81, 3.88, -4.68, -3.66, 0.96,
        0.75, 1.52, 2.11, 1.74, 0.78, -3.5, -3.32, -1.01
    ),
    'gd': (
        0.75, 1, 0, 0, 2.65, 0, 0, 2.95, 1.5, 2.4, 1.3, 0.69, 2.6, 0.59, 0.75, 0, 0.45,
        1.7, 3, 2.85
    ),
    'gy': (
        0.1, -1.42, 0.78, 0.83, -2.12, 0.33, -0.5, -1.13, 1.4, -1.18, -1.59, 0.48,
        0.73, 0.95, 1.91, 0.52, 0.07, -1.27, -0.51, -0.21
    ),
    'hw': (
        -0.5, -1.0, 3.0, 3.0, -2.5, 0.0, -0.5, -1.8, 3.0, -1.8, -1.3, 0.2, 0.0, 0.2,
        3.0, 0.3, -0.4, -1.5, -3.4, -2.3
    ),
    'ja': (
        0.28, 0.97, -0.52, -1.01, 0.46, 0.43, -0.31, 0.6, -1.62, 0.6, 0.43, -0.55,
        -0.42, -0.69, -1.14, -0.19, -0.32, 0.6, 0.29, -0.15
    ),
    'jo': (
        0.87, 1.52, 0.66, 0.67, 2.87, 0.1, 0.87, 3.15, 1.64, 2.17, 1.67, 0.09, 2.77, 0,
        0.85, 0.07, 0.07, 1.87, 3.77, 2.67
    ),
    'ju': (
        1.1, 2.5, -3.6, -3.2, 2.8, -0.64, -3.2, 4.5, -4.11, 3.8, 1.9, -3.5, -1.9,
        -3.68, -5.1, -0.5, -0.7, 4.2, -0.46, -1.3
    ),
    'kd': (
        1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5,
        -4.5, -0.8, -0.7, 4.2, -0.9, -1.3
    ),
    'ki': (
        -0.27, -1.05, 0.81, 1.17, -1.43, -0.16, 0.28, -0.77, 1.7, -1.1, -0.73, 0.81,
        -0.75, 1.1, 1.87, 0.42, 0.63, -0.4, -1.57, -0.56
    ),
    'mi': (
        5.33, 7.93, 3.59, 3.65, 9.03, 4.48, 5.1, 8.83, 2.95, 8.47, 8.95, 3.71, 3.87,
        3.87, 4.18, 4.09, 4.49, 7.63, 7.66, 5.89
    ),
    'pa': (
        2.1, 1.4, 10, 7.8, -9.2, 5.7, 2.1, -8, 5.7, -9.2, -4.2, 7, 2.1, 6, 4.2, 6.5,
        5.2, -3.7, -10, -1.9
    ),
    'po': (
        0.85, 2.1, -1.1, -0.79, 1.69, 0, 0.22, 3.14, -1.19, 1.99, 1.42, -0.48, -1.14,
        -0.42, 0.2, -0.52, -0.08, 2.53, 1.76, 1.37
    ),
    'rm': (
        0.39, 0.25, -3.81, -2.91, 2.27, 0, -0.64, 1.82, -2.77, 1.82, 0.96, -1.91, 0.99,
        -1.3, -3.95, -1.24, -1, 1.3, 2.13, 1.47
    ),
    'ro': (
        0.74, 0.91, 0.62, 0.62, 0.88, 0.72, 0.78, 0.88, 0.52, 0.85, 0.85, 0.63, 0.64,
        0.62, 0.64, 0.66, 0.7, 0.86, 0.85, 0.76
    ),
    'sw': (
        -0.4, 0.17, -1.31, -1.22, 1.92, -0.67, -0.64, 1.25, -0.67, 1.22, 1.02, -0.92,
        -0.49, -0.91, -0.59, -0.55, -0.28, 0.91, 0.5, 1.67
    ),
    'ta': (
        0.62, 0.29, -0.09, -0.74, 1.19, 0.48, -0.4, 1.38, -1.5, 1.53, 0.64, -0.78,
        0.12, -0.85, -2.53, -0.18, -0.05, 1.8, 0.81, 0.26
    ),
    'wi': (
        -0.3, 6.3, -1.4, 0, 7.5, 1.2, -1.3, 4.3, -3.6, 6.6, 2.5, -0.2, 2.2, -0.2, -1.1,
        -0.6, -2.2, 5.9, 7.9, 7.1
    ),
    'zi': (
        0.83, 1.48, 0.64, 0.65, 2.75, 0.1, 1.1, 3.07, 1.6, 2.52, 1.4, 0.09, 2.7, 0,
        0.83, 0.14, 0.54, 1.79, 0.31, 2.97
    ),
}

#: `protein_scale` scale name -> the key in `SCALES` it is.
GRAVY_SCALES: dict = {
    'Aboderin': 'ab',
    'AbrahamLeo': 'al',
    'Argos': 'ag',
    'BlackMould': 'bm',
    'BullBreese': 'bb',
    'Casari': 'cs',
    'Cid': 'ci',
    'Cowan3.4': 'Cowan3.4',
    'Cowan7.5': 'Cowan7.5',
    'Eisenberg': 'es',
    'Engelman': 'eg',
    'Fasman': 'fs',
    'Fauchere': 'fc',
    'GoldSack': 'gd',
    'Guy': 'gy',
    'Jones': 'jo',
    'Juretic': 'ju',
    'Kidera': 'ki',
    'KyteDoolitle': 'kd',
    'Miyazawa': 'mi',
    'Parker': 'pa',
    'Ponnuswamy': 'po',
    'Rose': 'ro',
    'Roseman': 'rm',
    'Sweet': 'sw',
    'Tanford': 'ta',
    'Wilson': 'wi',
    'Zimmerman': 'zi',
}

#: `IsoelectricPoint`'s charged residues, in the order the reference
#: builds its content dict in.
CHARGED_AAS: tuple = ('K', 'R', 'H', 'D', 'E', 'C', 'Y')

#: Groups that gain a proton.  `Nterm` first, and the order is the order
#: `charge_at_pH` adds them in.
POSITIVE_PKS: dict = {
    'Nterm': 7.5,
    'K': 10.0,
    'R': 12.0,
    'H': 5.98,
}

#: Groups that lose one.  `Cterm` first, likewise.
NEGATIVE_PKS: dict = {
    'Cterm': 3.55,
    'D': 4.05,
    'E': 4.45,
    'C': 9.0,
    'Y': 10.0,
}

#: Overrides for `positive_pKs['Nterm']`, by the residue it starts with.
PK_NTERMINAL: dict = {
    'A': 7.59,
    'M': 7.0,
    'S': 6.93,
    'P': 8.36,
    'T': 6.82,
    'V': 7.44,
    'E': 7.7,
}

#: Overrides for `negative_pKs['Cterm']`, by the residue it ends with.
PK_CTERMINAL: dict = {
    'D': 4.55,
    'E': 4.75,
}

