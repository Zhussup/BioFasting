"""The IUPAC tables `biofasting.sequtils` and `biofasting.checksum` need.

**Generated** by `tools/gen_iupac_data.py` -- do not edit.  Twelve tables copied
from `Bio.Data.IUPACData` and `Bio.SeqUtils._gc_values`, with the two complement
tables flattened from 256-byte `bytes.translate` strings into `str`, because
that is the shape `str.translate` wants.

The generator refuses to write unless it can rebuild the reference's own answers
from these tables alone: `molecular_weight` over every argument combination,
`gc_fraction(.., "weighted")`, and `complement`/`complement_rna` for all 256 byte
values.  The weight tables are the ones with something to get wrong -- a missing
letter is a `KeyError` the reference raises and this would answer `0.0` to.
"""

#: The Biopython release these were read from: 1.88
IUPAC_SOURCE_VERSION = '1.88'

# One-letter to three-letter residue names, the extended alphabet.
PROTEIN_1TO3: dict = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
    'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln',
    'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr', 'B': 'Asx',
    'X': 'Xaa', 'Z': 'Glx', 'J': 'Xle', 'U': 'Sec', 'O': 'Pyl',
}

# The inverse, as the reference writes it out separately.
PROTEIN_3TO1: dict = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H',
    'Ile': 'I', 'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q',
    'Arg': 'R', 'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y', 'Asx': 'B',
    'Xaa': 'X', 'Glx': 'Z', 'Xle': 'J', 'Sec': 'U', 'Pyl': 'O',
}

# Every DNA ambiguity code and the letters it stands for.
AMBIGUOUS_DNA_VALUES: dict = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'M': 'AC', 'R': 'AG', 'W': 'AT',
    'S': 'CG', 'Y': 'CT', 'K': 'GT', 'V': 'ACG', 'H': 'ACT', 'D': 'AGT', 'B': 'CGT',
    'X': 'GATC', 'N': 'GATC',
}

# The GC weight of every code, for gc_fraction's weighted mode.
GC_VALUES: dict = {
    'G': 1.0, 'C': 1.0, 'A': 0.0, 'T': 0.0, 'U': 0.0, 'S': 1.0, 'W': 0.0, 'M': 0.5,
    'R': 0.5, 'Y': 0.5, 'K': 0.5, 'V': 0.6666666666666666, 'B': 0.6666666666666666,
    'H': 0.3333333333333333, 'D': 0.3333333333333333, 'X': 0.5, 'N': 0.5,
}

# Average residue masses of the four deoxynucleotides, 5'-phosphate.
DNA_WEIGHTS: dict = {
    'A': 331.2218, 'C': 307.1971, 'G': 347.2212, 'T': 322.2085,
}

# The same, monoisotopic.
DNA_WEIGHTS_MONOISOTOPIC: dict = {
    'A': 331.06817, 'C': 307.056936, 'G': 347.063084, 'T': 322.056602,
}

# Average residue masses of the four ribonucleotides, 5'-phosphate.
RNA_WEIGHTS: dict = {
    'A': 347.2212, 'C': 323.1965, 'G': 363.2206, 'U': 324.1813,
}

# The same, monoisotopic.
RNA_WEIGHTS_MONOISOTOPIC: dict = {
    'A': 347.063084, 'C': 323.051851, 'G': 363.057999, 'U': 324.035867,
}

# Average residue masses of the twenty protein letters.
PROTEIN_WEIGHTS: dict = {
    'A': 89.0932, 'C': 121.1582, 'D': 133.1027, 'E': 147.1293, 'F': 165.1891,
    'G': 75.0666, 'H': 155.1546, 'I': 131.1729, 'K': 146.1876, 'L': 131.1729,
    'M': 149.2113, 'N': 132.1179, 'O': 255.3134, 'P': 115.1305, 'Q': 146.1445,
    'R': 174.201, 'S': 105.0926, 'T': 119.1192, 'U': 168.0532, 'V': 117.1463,
    'W': 204.2252, 'Y': 181.1885,
}

# The same, monoisotopic.
PROTEIN_WEIGHTS_MONOISOTOPIC: dict = {
    'A': 89.047678, 'C': 121.019749, 'D': 133.037508, 'E': 147.053158, 'F': 165.078979,
    'G': 75.032028, 'H': 155.069477, 'I': 131.094629, 'K': 146.105528, 'L': 131.094629,
    'M': 149.051049, 'N': 132.053492, 'O': 255.158292, 'P': 115.063329, 'Q': 146.069142,
    'R': 174.111676, 'S': 105.042593, 'T': 119.058243, 'U': 168.964203, 'V': 117.078979,
    'W': 204.089878, 'Y': 181.073893,
}

# Bio.Seq.complement, as a str.translate table over all 256 bytes.
DNA_COMPLEMENT = (
    '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?'
    '@TVGHEFCDIJMLKNOPQYSAABWXRZ[\\]^_`tvghefcdijmlknopqysaabwxrz{|}~\x7f'
    '\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0¡¢£¤¥¦§¨©ª«¬\xad®¯°±²³´µ¶·¸¹º»¼½¾¿'
    'ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖ×ØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõö÷øùúûüýþÿ'
)

# Bio.Seq.complement_rna, the same way.
RNA_COMPLEMENT = (
    '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?'
    '@UVGHEFCDIJMLKNOPQYSAABWXRZ[\\]^_`uvghefcdijmlknopqysaabwxrz{|}~\x7f'
    '\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0¡¢£¤¥¦§¨©ª«¬\xad®¯°±²³´µ¶·¸¹º»¼½¾¿'
    'ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖ×ØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõö÷øùúûüýþÿ'
)

__all__ = [
    'PROTEIN_1TO3',
    'PROTEIN_3TO1',
    'AMBIGUOUS_DNA_VALUES',
    'GC_VALUES',
    'DNA_WEIGHTS',
    'DNA_WEIGHTS_MONOISOTOPIC',
    'RNA_WEIGHTS',
    'RNA_WEIGHTS_MONOISOTOPIC',
    'PROTEIN_WEIGHTS',
    'PROTEIN_WEIGHTS_MONOISOTOPIC',
    'DNA_COMPLEMENT',
    'RNA_COMPLEMENT',
    "IUPAC_SOURCE_VERSION",
]
