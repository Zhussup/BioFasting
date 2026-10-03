"""The 27 NCBI genetic codes, as the strings NCBI publishes them.

**Generated** by `tools/gen_codon_tables.py` -- do not edit.  Each row is

    (id, name, alt_name, amino_acids, start_codons, stop_codons, dual)

where `amino_acids` is NCBI's 64 characters in the codon order TTT, TTC, TTA,
TTG, ... GGG, `*` marking a stop, and `dual` lists the codons that are a stop
*and* code for an amino acid -- context-dependent stops, which a string of this
shape cannot express and which the karyorelictid ciliates and their relatives
really have.

Read by `biofasting.codon_table`; the derivation lives in `_codon_table.py`.
"""

from __future__ import annotations
NCBI_TABLE_VERSION = '4.5'

#   id, name, alt_name, amino acids in CODON_ORDER, starts, stops, dual
TABLES: list[tuple] = [
    (1, 'Standard', 'SGC0',
     'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['TTG', 'CTG', 'ATG'], ['TAA', 'TAG', 'TGA'], []),
    (2, 'Vertebrate Mitochondrial', 'SGC1',
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
     ['ATT', 'ATC', 'ATA', 'ATG', 'GTG'], ['TAA', 'TAG', 'AGA', 'AGG'], []),
    (3, 'Yeast Mitochondrial', 'SGC2',
     'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATA', 'ATG', 'GTG'], ['TAA', 'TAG'], []),
    (4, 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma', 'SGC3',
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['TTA', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'], ['TAA', 'TAG'], []),
    (5, 'Invertebrate Mitochondrial', 'SGC4',
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
     ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'], ['TAA', 'TAG'], []),
    (6, 'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear', 'SGC5',
     'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TGA'], []),
    (9, 'Echinoderm Mitochondrial; Flatworm Mitochondrial', 'SGC8',
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
     ['ATG', 'GTG'], ['TAA', 'TAG'], []),
    (10, 'Euplotid Nuclear', 'SGC9',
     'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TAA', 'TAG'], []),
    (11, 'Bacterial, Archaeal and Plant Plastid', None,
     'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'], ['TAA', 'TAG', 'TGA'], []),
    (12, 'Alternative Yeast Nuclear', None,
     'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['CTG', 'ATG'], ['TAA', 'TAG', 'TGA'], []),
    (13, 'Ascidian Mitochondrial', None,
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
     ['TTG', 'ATA', 'ATG', 'GTG'], ['TAA', 'TAG'], []),
    (14, 'Alternative Flatworm Mitochondrial', None,
     'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
     ['ATG'], ['TAG'], []),
    (15, 'Blepharisma Macronuclear', None,
     'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TAA', 'TGA'], []),
    (16, 'Chlorophycean Mitochondrial', None,
     'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TAA', 'TGA'], []),
    (21, 'Trematode Mitochondrial', None,
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
     ['ATG', 'GTG'], ['TAA', 'TAG'], []),
    (22, 'Scenedesmus obliquus Mitochondrial', None,
     'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TCA', 'TAA', 'TGA'], []),
    (23, 'Thraustochytrium Mitochondrial', None,
     'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATT', 'ATG', 'GTG'], ['TTA', 'TAA', 'TAG', 'TGA'], []),
    (24, 'Pterobranchia Mitochondrial', None,
     'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
     ['TTG', 'CTG', 'ATG', 'GTG'], ['TAA', 'TAG'], []),
    (25, 'Candidate Division SR1 and Gracilibacteria', None,
     'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['TTG', 'ATG', 'GTG'], ['TAA', 'TAG'], []),
    (26, 'Pachysolen tannophilus Nuclear', None,
     'FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['CTG', 'ATG'], ['TAA', 'TAG', 'TGA'], []),
    (27, 'Karyorelict Nuclear', None,
     'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TGA'], [('TGA', 'W')]),
    (28, 'Condylostoma Nuclear', None,
     'FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TAA', 'TAG', 'TGA'], [('TAA', 'Q'), ('TAG', 'Q'), ('TGA', 'W')]),
    (29, 'Mesodinium Nuclear', None,
     'FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TGA'], []),
    (30, 'Peritrich Nuclear', None,
     'FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TGA'], []),
    (31, 'Blastocrithidia Nuclear', None,
     'FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['ATG'], ['TAA', 'TAG'], [('TAA', 'E'), ('TAG', 'E')]),
    (32, 'Balanophoraceae Plastid', None,
     'FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
     ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'], ['TAA', 'TGA'], []),
    (33, 'Cephalodiscidae Mitochondrial', None,
     'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
     ['TTG', 'CTG', 'ATG', 'GTG'], ['TAG'], []),
]

__all__ = ["TABLES", "NCBI_TABLE_VERSION"]
