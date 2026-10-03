"""`Bio.SeqUtils.IsoelectricPoint`: the charge of a protein at a given pH.

Forty lines of arithmetic and one recursion, and nothing in it is worth a kernel:
the heaviest input measured in `bench/targets.md` is a 1,024-residue protein at
**26 µs**, and it is spent in nine `10 ** x` calls per `charge_at_pH` and about
twenty of those per `pi()`, not in a loop over the sequence.  It is here because
`ProteinAnalysis.isoelectric_point` and `ProteinAnalysis.charge_at_pH` are two of
the rows in `Bio.SeqUtils.ProtParam`'s triage, and they are both this class
wearing a different name.

What the port has to be careful about is **order**, twice, and neither is
visible in the answer:

* `charge_at_pH` walks `pos_pKs` and `neg_pKs` in their insertion order and adds
  floats as it goes.  Those dicts come from `_protparam_data`, where the order is
  the reference's -- `Nterm`, `K`, `R`, `H` and `Cterm`, `D`, `E`, `C`, `Y` -- and
  `_update_pKs_tables` overwrites entries *in place* rather than rebuilding the
  dict, so an N-terminal alanine keeps `Nterm` where it was.
* the recursion in `pi()` is the reference's own bisection, down to the starting
  point of 7.775 and the 0.0001 interval, because a different search order lands
  on a different last bit.

The input convention is :mod:`biofasting.sequtils`'s: text, with a byte buffer
read as Latin-1, so that a `Seq`-like object behaves as it would through the
reference.  Nothing here declines -- `charge_at_pH` only ever looks up the seven
charged residues in a content dict that has all twenty, and an unknown residue
contributes nothing, which is also what the reference does.
"""

from __future__ import annotations

from ._protparam_data import (
    CHARGED_AAS,
    NEGATIVE_PKS,
    PK_CTERMINAL,
    PK_NTERMINAL,
    POSITIVE_PKS,
)
from .sequtils import _as_text

__all__ = ["IsoelectricPoint"]

#: The residues that carry a charge, in the order the reference lists them.
charged_aas = CHARGED_AAS

#: Groups that gain a proton -- `Nterm` first -- and groups that lose one.
positive_pKs = POSITIVE_PKS
negative_pKs = NEGATIVE_PKS

#: pK overrides for the residue the protein starts, or ends, with.
pKnterminal = PK_NTERMINAL
pKcterminal = PK_CTERMINAL


class IsoelectricPoint:
    """The isoelectric point of a protein, by bisection on its net charge.

    >>> protein = IsoelectricPoint("INGAR")
    >>> round(protein.pi(), 2)
    9.75
    >>> round(protein.charge_at_pH(7.0), 2)
    0.76

    `aa_content` is the reference's own short cut: pass a residue count that has
    already been computed -- `ProteinAnalysis.count_amino_acids` returns exactly
    that dict -- and the sequence is not counted again.
    """

    def __init__(self, protein_sequence, aa_content=None):
        self.sequence = _as_text(protein_sequence).upper()
        if not aa_content:
            from .protparam import ProteinAnalysis

            aa_content = ProteinAnalysis(self.sequence).count_amino_acids()
        self.charged_aas_content = self._select_charged(aa_content)

        self.pos_pKs, self.neg_pKs = self._update_pKs_tables()

    def _select_charged(self, aa_content):
        """The seven charged residues, plus the two termini, as floats.

        The reference builds this dict in `charged_aas` order and appends the
        termini; the order is what `charge_at_pH` adds in, so it is kept.
        """
        charged = {}
        for aa in charged_aas:
            charged[aa] = float(aa_content[aa])
        charged["Nterm"] = 1.0
        charged["Cterm"] = 1.0
        return charged

    def _update_pKs_tables(self):
        """Replace the terminal pKs when the sequence's own ends have their own.

        The reference assigns into copies of the module dicts rather than
        building new ones, so an overridden `Nterm` keeps its position at the
        front -- which matters, because that is the position it is summed in.
        """
        pos_pKs = positive_pKs.copy()
        neg_pKs = negative_pKs.copy()
        nterm, cterm = self.sequence[0], self.sequence[-1]
        if nterm in pKnterminal:
            pos_pKs["Nterm"] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs["Cterm"] = pKcterminal[cterm]
        return pos_pKs, neg_pKs

    def charge_at_pH(self, pH):
        """The net charge at `pH`, by the Henderson-Hasselbalch equation.

        Written as the reference writes it -- one division and one addition per
        group, in two loops, positive first -- because the twelve terms are
        floats and their sum is not associative.
        """
        positive_charge = 0.0
        for aa, pK in self.pos_pKs.items():
            partial_charge = 1.0 / (10 ** (pH - pK) + 1.0)
            positive_charge += self.charged_aas_content[aa] * partial_charge

        negative_charge = 0.0
        for aa, pK in self.neg_pKs.items():
            partial_charge = 1.0 / (10 ** (pK - pH) + 1.0)
            negative_charge += self.charged_aas_content[aa] * partial_charge

        return positive_charge - negative_charge

    def pi(self, pH=7.775, min_=4.05, max_=12):
        """The isoelectric point: the pH where the net charge is zero.

        A bisection, not a solve.  The reference starts at 7.775, which is the
        midpoint of `min_` and `max_` only by coincidence of the defaults, and
        recurses until the interval is narrower than 0.0001 -- about seventeen
        levels.  Reproduced as written, including the direction of the test, so
        that the answer is the same float and not merely the same pH.
        """
        charge = self.charge_at_pH(pH)
        if max_ - min_ > 0.0001:
            if charge > 0.0:
                min_ = pH
            else:
                max_ = pH
            next_pH = (min_ + max_) / 2
            return self.pi(next_pH, min_, max_)
        return pH
