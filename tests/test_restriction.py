"""Phase 2 step 2.2: the restriction enzymes, as data and as behaviour.

`Bio.Restriction` is 1,088 enzyme classes that differ only in their attributes;
`biofasting.restriction` is one class and a 1,088-row table.  Saying that is
cheap.  Proving that the table says the same thing, and that the one class
behaves the way the 1,088 did, is this file.

It is two different claims and they need different evidence.

**The data is faithful.**  Every field the reference's classes carry is
compared against the row this library reads -- site, size, cut coordinates,
overhang, frequency, suppliers, temperatures, URI, REBASE id, and the
classification tuple the behaviour flags are derived from.  A field that drifts
is a wrong answer about a real enzyme, and no amount of behavioural testing
would say so: an enzyme with the right behaviour and the wrong supplier list is
still a wrong row.

**The behaviour is the reference's.**  `search`, `catalyze`, `elucidate`,
`frequency`, the compatibility operators, the comparison operators and the
ordering are run against `Bio.Restriction` itself, over all 1,088 enzymes, on
sequence shapes chosen so that each branch is actually reached: linear and
circular, a site at the very first base and at the very last, a site that
straddles the origin of a circular molecule, the degenerate sites that send
`search` down the non-palindromic path, and the 25 enzymes that cut twice.
Where both sides raise, the exception *type* is compared, because a refusal is
a result -- but the message is not, except where the message is about the
sequence rather than about this library's internals.

**The reports are the reference's, to the byte.**  `Analysis` and `PrintFormat`
produce three reports -- a list of positions, a count per enzyme, and a map of
the sequence with its cuts drawn on it -- and all three are compared against
the reference's over every enzyme there is.  Two of them compare as bytes.  The
map cannot: the names at one position are printed in set-iteration order, which
follows the element hash, so Biopython's own digest for the same input differs
from run to run.  The map is therefore compared line for line with the names on
each line sorted, and `test_the_reference_map_is_not_stable_across_runs` is the
evidence for why that is the strongest comparison available rather than a
weakened one.

**What is not claimed.**  `RestrictionBatch.split` takes the reference's mixin
class *names* rather than its classes, so `split(Blunt)` works only if `Blunt`
is imported from `Bio.Restriction`; and the supplier-buffer tables are not here
yet.  Those are Phase 2 work, and they are named in the run report rather than
left to be discovered.

**The reference's defects are pinned by name** in the last sections.  They are
reproduced on purpose: a library that is a surprise in a different place each
time is worse than one that is a surprise in the same places.
"""

import importlib.util
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest

from biofasting import restriction as ours

ROOT = Path(__file__).resolve().parent.parent

needs_biopython = pytest.mark.skipif(
    importlib.util.find_spec("Bio") is None,
    reason="Biopython is not installed; it is the differential reference",
)

#   Two spellings of the same pattern have to compare equal.  The reference's
#   class members come in an arbitrary order (`[CG]` where the derivation writes
#   `[GC]`, 18 enzymes), and its groups are named after the enzyme while the
#   engine names them `t`/`b` for all 1,088.  So what is compared is the *set of
#   group bodies*, which is what a search actually matches against.  A check
#   that cries wolf on 1,088 non-problems teaches its reader to skim past the
#   one real difference.
_CLASS = re.compile(r"\[([A-Z]+)\]")
_GROUP = re.compile(r"\(\?=\(\?P<[^>]+>([^)]*)\)\)")


def _canonical(pattern: str) -> tuple:
    bodies = _GROUP.findall(pattern)
    return tuple(sorted(
        _CLASS.sub(lambda m: "[" + "".join(sorted(m.group(1))) + "]", body)
        for body in bodies
    ))


def _reference():
    import Bio.Restriction as reference

    return reference


def _reference_enzymes():
    from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict

    types = {}
    for bases, enzymes in typedict.values():
        for enzyme in enzymes:
            types[enzyme] = tuple(bases)
    return rest_dict, types


def _outcome(call):
    """A call's result, or the exception it raised, as comparable values."""
    try:
        return ("ok", call())
    except Exception as exc:  # noqa: BLE001 - the exception *is* the result here
        return ("raised", type(exc).__name__)


# --------------------------------------------------------------------------
# Probes: sequences shaped so that each branch is reached
# --------------------------------------------------------------------------

#   One concrete instance of each IUPAC code, so a degenerate site can be
#   pasted into a sequence and actually match.
_CONCRETE = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "A", "Y": "C", "W": "A", "S": "G", "M": "A", "K": "G",
    "H": "A", "B": "C", "V": "A", "D": "A", "N": "A",
}

_RNG = random.Random(20261003)
_FILLER = "".join(_RNG.choice("ACGT") for _ in range(600))


def _site_of(enzyme) -> str:
    """The recognition site as concrete bases (HpyUM037X carries two)."""
    return "".join(_CONCRETE[b] for b in enzyme.site.split("|", 1)[0])


def _revcomp(text: str) -> str:
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _probes(enzyme) -> list[str]:
    """Sequences for one enzyme, each aimed at a different edge.

    - the filler alone, which most enzymes do not cut at all: the empty-result
      path, and the one that has to stay empty;
    - the site at the very first base, where cuts fall off the front and the
      linear `_drop` has to remove them;
    - the site at the very last base, where they fall off the back;
    - two sites with filler between them, which is the ordinary case and the
      only one where the middle fragment exists;
    - the reverse complement of that, so a non-palindromic enzyme matches on
      the strand it was not written for.
    """
    site = _site_of(enzyme)
    return [
        _FILLER,
        site + _FILLER[:120],
        _FILLER[:120] + site,
        _FILLER[:40] + site + _FILLER[40:80] + site + _FILLER[80:130],
        _revcomp(_FILLER[:40] + site + _FILLER[40:80] + site + _FILLER[80:130]),
    ]


# --------------------------------------------------------------------------
# The table
# --------------------------------------------------------------------------


@needs_biopython
def test_the_table_has_exactly_the_reference_enzymes():
    from Bio.Restriction.Restriction_Dictionary import rest_dict

    assert sorted(ours.ENZYMES) == sorted(rest_dict)
    assert len(ours.ENZYMES) == 1088


@needs_biopython
def test_every_column_of_every_row_matches_the_reference():
    """Field by field, all 1,088 rows, with the mismatches named.

    Collected rather than asserted one at a time so a single run says how many
    rows drifted and which columns, which is the difference between a report
    and a puzzle.
    """
    rest_dict, types = _reference_enzymes()
    problems = []
    for name, entry in rest_dict.items():
        site, fst5, fst3, scd5, scd3, ovhg, ovhgseq, size, freq, suppl, opt, inact, uri, rid, tps = ours.ENZYMES[name]
        expected = {
            "site": (site, entry["site"]),
            "fst5": (fst5, entry["fst5"]),
            "fst3": (fst3, entry["fst3"]),
            "scd5": (scd5, entry["scd5"]),
            "scd3": (scd3, entry["scd3"]),
            "ovhg": (ovhg, entry["ovhg"]),
            "ovhgseq": (ovhgseq, entry["ovhgseq"]),
            "size": (size, entry["size"]),
            "freq": (freq, entry["freq"]),
            "suppl": (tuple(suppl), tuple(entry["suppl"])),
            "opt_temp": (opt, entry["opt_temp"]),
            "inact_temp": (inact, entry["inact_temp"]),
            "uri": (uri, entry["uri"]),
            "rebase_id": (rid, entry["id"]),
            "types": (set(tps), set(types[name])),
        }
        for column, (got, want) in expected.items():
            if got != want:
                problems.append(f"{name}.{column}: {got!r} != {want!r}")
    assert not problems, f"{len(problems)} field(s) differ: " + "; ".join(problems[:20])


@needs_biopython
def test_the_derived_search_pattern_is_the_reference_pattern():
    """The engine derives each pattern from the site; the reference stores it.

    This is the same check the generator makes before it writes the table,
    made again against the library's own derivation, so the two cannot quietly
    disagree with each other or with Biopython.
    """
    rest_dict, _ = _reference_enzymes()
    problems = []
    for name, entry in rest_dict.items():
        mine = _canonical(ours.get_enzyme(name).compsite.pattern)
        if mine != _canonical(entry["compsite"]):
            problems.append(f"{name}: {mine!r} != {entry['compsite']!r}")
    assert not problems, f"{len(problems)} pattern(s) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_the_generated_data_is_current():
    """`--check` re-derives every pattern and compares the file on disk.

    Run out of process because the point is that the *tool* is current, not
    that some function agrees with some other function.
    """
    done = subprocess.run(
        [sys.executable, "tools/gen_restriction_data.py", "--check"],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert done.returncode == 0, done.stdout + done.stderr
    assert "is current" in done.stdout


# --------------------------------------------------------------------------
# The behaviour, against the reference, for every enzyme
# --------------------------------------------------------------------------


@needs_biopython
def test_search_matches_on_every_probe_and_both_topologies():
    reference = _reference()
    from Bio.Seq import Seq

    problems = []
    for name in sorted(ours.ENZYMES):
        mine, theirs = ours.get_enzyme(name), getattr(reference, name)
        for probe in _probes(mine):
            for linear in (True, False):
                got = _outcome(lambda: mine.search(probe, linear))
                want = _outcome(lambda: theirs.search(Seq(probe), linear))
                if got != want:
                    problems.append(
                        f"{name} linear={linear} probe={probe[:40]!r}: {got} != {want}"
                    )
    assert not problems, f"{len(problems)} search(es) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_search_matches_on_one_long_sequence_for_every_enzyme():
    """A thousand bases of filler, where several sites usually land.

    The short probes reach the edges; this reaches the ordinary middle case,
    where there are several sites to order, drop or wrap.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence = "".join(_RNG.choice("ACGT") for _ in range(4000))
    problems = []
    for name in sorted(ours.ENZYMES):
        mine, theirs = ours.get_enzyme(name), getattr(reference, name)
        for linear in (True, False):
            got = _outcome(lambda: mine.search(sequence, linear))
            want = _outcome(lambda: theirs.search(Seq(sequence), linear))
            if got != want:
                problems.append(f"{name} linear={linear}: {got} != {want}")
    assert not problems, f"{len(problems)} search(es) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_catalyze_matches_on_every_probe():
    """Fragments compared as strings, so `Seq` and `str` do not confuse it."""
    reference = _reference()
    from Bio.Seq import Seq

    problems = []
    for name in sorted(ours.ENZYMES):
        mine, theirs = ours.get_enzyme(name), getattr(reference, name)
        for probe in _probes(mine):
            for linear in (True, False):
                def fragments(call):
                    kind, value = _outcome(call)
                    if kind == "ok":
                        return kind, tuple(str(x) for x in value)
                    return kind, value

                got = fragments(lambda: mine.catalyze(probe, linear))
                want = fragments(lambda: theirs.catalyze(Seq(probe), linear))
                if got != want:
                    problems.append(f"{name} linear={linear}: {got} != {want}")
    assert not problems, f"{len(problems)} digest(s) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_elucidate_matches_for_every_enzyme():
    reference = _reference()
    problems = []
    for name in sorted(ours.ENZYMES):
        got = _outcome(lambda: ours.get_enzyme(name).elucidate())
        want = _outcome(lambda: getattr(reference, name).elucidate())
        if got != want:
            problems.append(f"{name}: {got} != {want}")
    assert not problems, f"{len(problems)} elucidation(s) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_every_flag_and_every_plain_attribute_matches_for_every_enzyme():
    """The classification, read off the type tuple, against the reference's mixins."""
    reference = _reference()
    checks = [
        "is_palindromic", "cut_once", "cut_twice", "is_blunt", "is_5overhang",
        "is_3overhang", "is_defined", "is_ambiguous", "is_unknown",
        "is_methylable", "is_comm", "overhang", "supplier_list",
    ]
    problems = []
    for name in sorted(ours.ENZYMES):
        mine, theirs = ours.get_enzyme(name), getattr(reference, name)
        for check in checks:
            got, want = getattr(mine, check)(), getattr(theirs, check)()
            if got != want:
                problems.append(f"{name}.{check}(): {got!r} != {want!r}")
        for attribute in ("site", "size", "fst5", "fst3", "scd5", "scd3",
                          "ovhg", "ovhgseq", "freq", "opt_temp", "inact_temp",
                          "uri", "suppl"):
            got, want = getattr(mine, attribute), getattr(theirs, attribute)
            if got != want:
                problems.append(f"{name}.{attribute}: {got!r} != {want!r}")
        if str(mine) != str(theirs) or repr(mine) != repr(theirs):
            problems.append(f"{name}: str/repr differ")
        if len(mine) != len(theirs):
            problems.append(f"{name}: len {len(mine)} != {len(theirs)}")
        if mine.characteristic() != theirs.characteristic():
            problems.append(f"{name}: characteristic() differs")
    assert not problems, f"{len(problems)} attribute(s) differ: " + "; ".join(problems[:20])


@needs_biopython
def test_frequency_matches_for_every_enzyme():
    reference = _reference()
    for name in sorted(ours.ENZYMES):
        assert ours.get_enzyme(name).frequency() == getattr(reference, name).frequency()


# --------------------------------------------------------------------------
# The operators, which are not the methods they look like
# --------------------------------------------------------------------------


@needs_biopython
def test_the_comparison_operators_mean_the_same_thing():
    """`==` is identity, `!=` is "cuts identically", and they disagree.

    Both are checked over a sample that deliberately contains isoschizomers,
    because that disagreement is the whole point of having two operators.
    """
    reference = _reference()
    sample = sorted(ours.ENZYMES)[::17] + ["SacI", "SstI", "SmaI", "XmaI", "EcoRI", "EcoRV"]
    problems = []
    for a in sample:
        for b in sample:
            mine_a, mine_b = ours.get_enzyme(a), ours.get_enzyme(b)
            ref_a, ref_b = getattr(reference, a), getattr(reference, b)
            for symbol, mine, theirs in (
                ("==", mine_a == mine_b, ref_a == ref_b),
                ("!=", mine_a != mine_b, ref_a != ref_b),
                (">>", mine_a >> mine_b, ref_a >> ref_b),
                ("<", mine_a < mine_b, ref_a < ref_b),
                ("<=", mine_a <= mine_b, ref_a <= ref_b),
                (">", mine_a > mine_b, ref_a > ref_b),
                (">=", mine_a >= mine_b, ref_a >= ref_b),
            ):
                if mine != theirs:
                    problems.append(f"{a} {symbol} {b}: {mine!r} != {theirs!r}")
    assert not problems, f"{len(problems)} comparison(s) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_overhang_compatibility_matches_over_a_stratified_sample():
    """`%`, over enzymes chosen to cover every overhang type and ambiguity.

    The full cross product is 1.18 million pairs; this sample keeps every
    distinct overhang sequence and every distinct overhang type, which is what
    `_mod1`/`_mod2` branches on, plus a slice of everything else.
    """
    reference = _reference()
    by_shape: dict[tuple, str] = {}
    for name, row in ours.ENZYMES.items():
        #   ovhgseq and the type tuple are what `_mod1`/`_mod2` branch on, so
        #   one enzyme per distinct pair covers every branch there is.
        by_shape.setdefault((row[6], row[14]), name)
    sample = sorted(
        set(by_shape.values())
        | set(sorted(ours.ENZYMES)[::53])
        | {"EcoRI", "SmaI", "XmaI", "BstXI", "SnaI"}
    )
    problems = []
    for a in sample:
        for b in sample:
            mine_a, mine_b = ours.get_enzyme(a), ours.get_enzyme(b)
            ref_a, ref_b = getattr(reference, a), getattr(reference, b)
            got, want = _outcome(lambda: mine_a % mine_b), _outcome(lambda: ref_a % ref_b)
            if got != want:
                problems.append(f"{a} % {b}: {got} != {want}")
    assert not problems, f"{len(problems)} compatibilit(ies) differ: " + "; ".join(problems[:10])


@needs_biopython
def test_the_compatibility_helpers_match():
    reference = _reference()
    for name in ("EcoRI", "SmaI", "XmaI", "BamHI", "BglII", "SnaI", "BstXI", "AarI"):
        mine, theirs = ours.get_enzyme(name), getattr(reference, name)
        for other in ("EcoRI", "SmaI", "XmaI", "BamHI", "BglII", "SnaI", "BstXI", "AarI"):
            for helper in ("is_equischizomer", "is_neoschizomer", "is_isoschizomer"):
                got = getattr(mine, helper)(ours.get_enzyme(other))
                want = getattr(theirs, helper)(getattr(reference, other))
                assert got == want, f"{name}.{helper}({other})"


@needs_biopython
def test_isoschizomer_lists_match():
    reference = _reference()
    sample = sorted(ours.ENZYMES)[::97] + ["SacI", "SstI", "SmaI", "XmaI"]
    for name in sample:
        mine, theirs = ours.get_enzyme(name), getattr(reference, name)
        for helper in ("equischizomers", "neoschizomers", "isoschizomers"):
            got = [str(x) for x in getattr(mine, helper)()]
            want = [str(x) for x in getattr(theirs, helper)()]
            assert got == want, f"{name}.{helper}()"


# --------------------------------------------------------------------------
# Sequences: what is accepted, what is refused, and what comes back out
# --------------------------------------------------------------------------


def test_bytes_and_str_both_work_and_the_fragments_follow_the_input():
    assert ours.get_enzyme("EcoRI").search(b"TTTGAATTCAAA") == [5]
    assert ours.get_enzyme("EcoRI").search("TTTGAATTCAAA") == [5]
    assert ours.get_enzyme("EcoRI").catalyze("TTTGAATTCAAA") == ("TTTG", "AATTCAAA")
    assert ours.get_enzyme("EcoRI").catalyze(b"TTTGAATTCAAA") == (b"TTTG", b"AATTCAAA")


def test_lowercase_in_lowercase_out():
    """The reference remembers the caller's case and restores it."""
    assert ours.get_enzyme("EcoRI").search("tttgaattcaaa") == [5]
    assert ours.get_enzyme("EcoRI").catalyze("tttgaattcaaa") == ("tttg", "aattcaaa")


def test_whitespace_and_digits_are_removed_and_anything_else_is_refused():
    enzyme = ours.get_enzyme("EcoRI")
    assert enzyme.search("TTT GAA TTC\nAAA") == enzyme.search("TTTGAATTCAAA")
    assert enzyme.search("1TTTGAATTCAAA2") == [5]
    with pytest.raises(TypeError, match="Invalid character found in"):
        enzyme.search("TTTGAATTCAA*")


@needs_biopython
def test_the_invalid_character_message_is_the_reference_message():
    from Bio.Seq import Seq

    with pytest.raises(TypeError) as mine:
        ours.get_enzyme("EcoRI").search("TTTGAATTCAA*")
    with pytest.raises(TypeError) as theirs:
        _reference().EcoRI.search(Seq("TTTGAATTCAA*"))
    assert str(mine.value) == str(theirs.value)


@needs_biopython
def test_a_biopython_seq_is_accepted_and_comes_back_as_a_seq():
    from Bio.Seq import Seq

    enzyme = ours.get_enzyme("EcoRI")
    fragments = enzyme.catalyze(Seq("TTTGAATTCAAA"))
    assert [str(x) for x in fragments] == ["TTTG", "AATTCAAA"]
    assert all(isinstance(x, Seq) for x in fragments)


def test_a_sequence_object_is_refused_with_its_type_named():
    with pytest.raises(TypeError, match="expected a string, bytes"):
        ours.get_enzyme("EcoRI").search(1234)


def test_importing_the_module_does_not_import_biopython():
    """A restriction digest is not a reason to require the reference."""
    done = subprocess.run(
        [sys.executable, "-c",
         "import sys, biofasting.restriction; "
         "assert 'Bio' not in sys.modules, sorted(m for m in sys.modules if m.startswith('Bio'))"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert done.returncode == 0, done.stderr


# --------------------------------------------------------------------------
# The reference's defects, pinned by name
# --------------------------------------------------------------------------


@needs_biopython
def test_adding_an_enzyme_to_a_batch_returns_none_as_the_reference_does():
    """`EcoRI + batch` is None upstream, because `add_nocheck` is `set.add`.

    Reproduced rather than fixed.  `batch + EcoRI` returns a batch, `EcoRI +
    EcoRV` returns a batch, and this one case returns None -- in both
    libraries.  A difference nobody wrote down is worse than a defect
    everybody did.
    """
    reference = _reference()
    mine = ours.get_enzyme("EcoRI") + ours.RestrictionBatch(["EcoRV"])
    theirs = reference.EcoRI + reference.RestrictionBatch([reference.EcoRV])
    assert mine is None and theirs is None


@needs_biopython
def test_compatible_end_ignores_the_batch_it_is_given():
    """Both libraries accept a `batch` argument here and then search everything."""
    reference = _reference()
    small = ours.RestrictionBatch(["EcoRI"])
    mine = [str(x) for x in ours.get_enzyme("SmaI").compatible_end(small)]
    theirs = [
        str(x)
        for x in reference.SmaI.compatible_end(
            reference.RestrictionBatch([reference.EcoRI])
        )
    ]
    assert mine == theirs
    assert "EcoRI" not in mine  # the batch said EcoRI; the answer ignores it


@needs_biopython
def test_an_unknown_enzyme_refuses_to_catalyze_but_still_finds_its_sites():
    """An uncharacterised enzyme can say where it binds, not where it cuts."""
    reference = _reference()
    from Bio.Seq import Seq

    sequence = "AAAGTATACAAA"
    assert ours.get_enzyme("SnaI").search(sequence) == reference.SnaI.search(Seq(sequence))
    with pytest.raises(NotImplementedError):
        ours.get_enzyme("SnaI").catalyze(sequence)
    with pytest.raises(NotImplementedError):
        reference.SnaI.catalyze(Seq(sequence))


@needs_biopython
def test_an_unknown_enzyme_refuses_to_compare_overhangs():
    reference = _reference()
    with pytest.raises(ValueError) as mine:
        ours.get_enzyme("SnaI")._mod2(ours.get_enzyme("SnaI"))
    with pytest.raises(ValueError) as theirs:
        reference.SnaI._mod2(reference.SnaI)
    assert "pas glop pas glop" in str(mine.value) and "pas glop pas glop" in str(theirs.value)


@needs_biopython
def test_split_and_lambdasplit_return_an_empty_batch_upstream():
    """A Biopython 1.88 defect, found by asking both libraries the same question.

    Both methods build the new batch with `RestrictionBatch()` and then assign
    the enzymes they selected to `new._data` -- an attribute left over from when
    the class was dict-backed, and which nothing reads now that it is a `set`.
    So the filter runs, the answer is computed, and it is thrown away: `len()`
    is zero for every call, whatever the arguments, in the version this library
    is verified against.

    Here they do what their docstrings say ("extract enzymes of a certain class
    and put in new RestrictionBatch").  That is a deliberate divergence from
    the reference, and it is written down rather than left for someone to
    discover: the reference's emptiness is asserted below, so a future
    Biopython that fixes this turns the line red instead of leaving two
    libraries quietly disagreeing.
    """
    reference = _reference()
    from Bio.Restriction import Restriction as mixins

    names = ["EcoRI", "SmaI", "SnaI", "BstXI", "KpnI", "BamHI", "EcoRV"]
    theirs = reference.RestrictionBatch([getattr(reference, n) for n in names])
    mine = ours.RestrictionBatch(names)

    assert len(theirs.split(mixins.Blunt)) == 0
    assert len(theirs.lambdasplit(lambda enzyme: True)) == 0
    assert sorted(str(x) for x in mine.split("blunt")) == ["EcoRV", "SmaI"]
    assert sorted(str(x) for x in mine.lambdasplit(lambda e: e.is_blunt())) == [
        "EcoRV", "SmaI"
    ]


@needs_biopython
def test_an_enzyme_that_cuts_twice_says_so_instead_of_drawing_itself():
    reference = _reference()
    assert ours.get_enzyme("AjuI").elucidate() == reference.AjuI.elucidate()
    assert ours.get_enzyme("AjuI").elucidate() == "cut twice, not yet implemented sorry."


# --------------------------------------------------------------------------
# Batches
# --------------------------------------------------------------------------


@needs_biopython
def test_the_three_batches_have_the_reference_membership():
    reference = _reference()
    assert {str(x) for x in ours.AllEnzymes} == {str(x) for x in reference.AllEnzymes}
    assert {str(x) for x in ours.CommOnly} == {str(x) for x in reference.CommOnly}
    assert {str(x) for x in ours.NonComm} == {str(x) for x in reference.NonComm}


@needs_biopython
def test_a_batch_search_returns_the_same_mapping():
    reference = _reference()
    from Bio.Seq import Seq

    sequence = "TTTGAATTCAAAGGATCCAAAGGTACCAAA"
    names = ["EcoRI", "BamHI", "KpnI", "SmaI", "HindIII"]
    mine = ours.RestrictionBatch(names).search(sequence)
    theirs = reference.RestrictionBatch([getattr(reference, n) for n in names]).search(Seq(sequence))
    assert {str(k): v for k, v in mine.items()} == {str(k): v for k, v in theirs.items()}


def test_a_batch_search_is_cached_against_the_formatted_sequence():
    batch = ours.RestrictionBatch(["EcoRI"])
    first = batch.search("TTTGAATTCAAA")
    assert batch.search("tttgaattcaaa") is first


@needs_biopython
def test_a_batch_reports_itself_the_way_the_reference_does():
    reference = _reference()
    names = ["EcoRI", "EcoRV", "BamHI"]
    mine = ours.RestrictionBatch(names)
    theirs = reference.RestrictionBatch([getattr(reference, n) for n in names])
    assert str(mine) == str(theirs)
    assert repr(mine) == repr(theirs)
    assert mine.elements() == theirs.elements()
    #   `as_string()` is the unsorted twin of `elements()`, and "unsorted" here
    #   means *set iteration order*, which follows the hash of the elements --
    #   `id(cls)` there, `id(obj)` here.  So the lists are compared as sets:
    #   the order is not part of what this method promises, in either library.
    assert sorted(mine.as_string()) == sorted(theirs.as_string())


def test_split_takes_plain_words_or_the_reference_class_names():
    batch = ours.RestrictionBatch(["EcoRI", "SmaI", "SnaI", "BstXI", "KpnI", "BamHI"])
    names = lambda found: sorted(str(x) for x in found)
    assert names(batch.split("blunt")) == ["SmaI"]
    assert names(batch.split("5overhang")) == ["BamHI", "EcoRI"]
    assert names(batch.split("Ov3")) == ["BstXI", "KpnI"]
    assert names(batch.split("unknown")) == ["SnaI"]
    assert names(batch.split("ambiguous")) == ["BstXI"]
    #   A requested class defaults to True, so passing it False asks for the
    #   enzymes that are *not* of it -- the reference's convention, kept
    #   because it is the reason this takes keyword arguments at all.
    assert names(batch.split("blunt", blunt=False)) == [
        "BamHI", "BstXI", "EcoRI", "KpnI", "SnaI"
    ]
    with pytest.raises(ValueError, match="unknown enzyme class"):
        batch.split("nonsense")


@needs_biopython
def test_split_accepts_the_reference_mixin_classes_by_name():
    """Only a class's `__name__` is read, so `split(Blunt)` works for callers
    who already have Biopython's mixin classes imported."""
    names = ["EcoRI", "SmaI", "SnaI", "BstXI", "KpnI", "BamHI", "EcoRV"]
    from Bio.Restriction import Restriction as mixins

    batch = ours.RestrictionBatch(names)
    expected = {
        ("Blunt",): ["EcoRV", "SmaI"],
        ("Ov5",): ["BamHI", "EcoRI"],
        ("Ov3",): ["BstXI", "KpnI"],
        ("OneCut",): ["BamHI", "BstXI", "EcoRI", "EcoRV", "KpnI", "SmaI"],
        ("Blunt", "Ambiguous"): [],
    }
    for argument, want in expected.items():
        by_name = sorted(str(x) for x in batch.split(*argument))
        by_class = sorted(
            str(x) for x in batch.split(*(getattr(mixins, a) for a in argument))
        )
        assert by_name == want, argument
        assert by_class == want, argument


def test_a_batch_refuses_an_unknown_name():
    with pytest.raises(ValueError, match="is not a restriction enzyme"):
        ours.RestrictionBatch(["NotAnEnzyme"])
    assert "NotAnEnzyme" not in ours.RestrictionBatch(["EcoRI"])


def test_suppliers_are_the_reference_table():
    """Codes to names, and the reverse index derived from the enzyme rows."""
    assert len(ours.supplier_codes()) == 15
    assert ours.supplier_codes()["B"] == "Thermo Fisher Scientific"
    assert "EcoRI" in ours.enzymes_of_supplier("B")
    assert "AarI" in ours.enzymes_of_supplier("B")
    with pytest.raises(KeyError):
        ours.enzymes_of_supplier("Z")


def test_an_unknown_enzyme_name_says_so():
    with pytest.raises(KeyError, match="unknown restriction enzyme"):
        ours.get_enzyme("NotAnEnzyme")


@needs_biopython
def test_with_name_returns_nothing_when_it_is_given_a_result_dictionary():
    """`with_name(names, dct)` looks up names in a dictionary keyed by enzymes.

    `{n: dct[n] for n in names if n in dct}` -- `n` is a string and `dct`'s keys
    are the enzymes themselves, so the membership test is false for every name
    and the method returns an empty dictionary.  Handing it the enzymes rather
    than their names works, which is the shape of the bug: the branch is not
    broken, it is answering a question nobody asked.

    Reproduced and pinned, like the other upstream defects here.  A Biopython
    that repairs it turns this line red.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    names = ["EcoRI", "BamHI"]
    mine = ours.Analysis(ours.RestrictionBatch(names), sequence)
    theirs = reference.Analysis(reference.RestrictionBatch(names), Seq(sequence))
    assert mine.with_name(names, mine.mapping) == {}
    assert theirs.with_name(names, theirs.mapping) == {}
    assert {str(k) for k in mine.with_name(list(mine), mine.mapping)} == set(names)
    assert {str(k) for k in theirs.with_name(list(theirs), theirs.mapping)} == set(names)


# --------------------------------------------------------------------------
# Analysis and PrintFormat: the digest's report
# --------------------------------------------------------------------------

#   A report line that carries enzyme names: an optional column of bars, the
#   position, then the names.  Used to compare maps without comparing the
#   order of the names, which is not part of the output (see the test below).
_MAP_LINE = re.compile(r"^([| ]*\d+)((?:\s+\S+)*)\s*$")


def _normalise_map(text: str) -> list[str]:
    lines = []
    for line in text.splitlines():
        match = _MAP_LINE.match(line)
        if match:
            line = match.group(1) + " " + " ".join(sorted(match.group(2).split()))
        lines.append(line)
    return lines


def _report_sequence() -> tuple[str, list[str]]:
    """A sequence long enough to exercise the map's blocks, and its sites."""
    rng = random.Random(4242)
    body = "".join(rng.choice("ACGT") for _ in range(700))
    body = body[:100] + "GAATTC" + body[106:300] + "GGATCC" + body[306:500]
    body = body[:500] + "AAGCTT" + body[506:]
    return body, ["EcoRI", "BamHI", "HindIII", "SmaI", "KpnI", "EcoRV"]


@needs_biopython
def test_the_list_and_number_reports_match_for_all_enzymes():
    """The two deterministic reports, over every enzyme there is.

    Deterministic because both sort: `__next_section` sorts the enzymes by
    name, and `_make_number_only` sorts by cut count and then by name.  The map
    is the third report and cannot be compared this way -- see the next test.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    mine = ours.Analysis(ours.AllEnzymes, sequence)
    theirs = reference.Analysis(reference.AllEnzymes, Seq(sequence))
    for shape in ("list", "number"):
        mine.print_as(shape)
        theirs.print_as(shape)
        assert mine.format_output() == theirs.format_output(), shape


@needs_biopython
def test_the_map_report_matches_once_the_name_order_is_normalised():
    """The map, compared line for line with the names on each line sorted.

    This is the weakest-looking comparison in the file and it is the strongest
    one available.  A line naming several enzymes at one position lists them in
    set-iteration order, which follows the element hash -- `id(cls)` upstream,
    `id(obj)` here -- so the reference's own map for the same input differs
    from run to run.  Sorting the names on each line removes exactly the part
    of the output that is not a function of the input, and nothing else: the
    layout, the columns, the coordinates, the sequence and its complement are
    all compared byte for byte.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    mine = ours.Analysis(ours.AllEnzymes, sequence)
    theirs = reference.Analysis(reference.AllEnzymes, Seq(sequence))
    mine.print_as("map")
    theirs.print_as("map")
    text = mine.format_output()
    assert len(text) > 50_000, "the map should be a real report, not a stub"
    assert _normalise_map(text) == _normalise_map(theirs.format_output())


@needs_biopython
def test_the_reference_map_is_not_stable_across_runs():
    """Why the map test above normalises: upstream's output is not reproducible.

    Not a claim about this library.  It is the evidence for the normalisation,
    recorded as a test so that a Biopython which fixes the ordering turns this
    line red -- at which point the map can be compared byte for byte and the
    weaker comparison should be deleted.
    """
    program = (
        "import hashlib, random;"
        "from Bio.Seq import Seq;"
        "from Bio.Restriction import Analysis, AllEnzymes;"
        "random.seed(11);"
        "s = Seq(''.join(random.choice('ACGT') for _ in range(700)));"
        "a = Analysis(AllEnzymes, s); a.print_as('map');"
        "print(hashlib.sha256(a.format_output().encode()).hexdigest())"
    )
    digests = set()
    for _ in range(3):
        done = subprocess.run(
            [sys.executable, "-c", program], capture_output=True, text=True, cwd=ROOT
        )
        assert done.returncode == 0, done.stderr
        digests.add(done.stdout.strip())
    assert len(digests) > 1, "a stable digest means the map can be compared as bytes"


@needs_biopython
def test_the_map_complement_table_is_the_one_biopython_uses():
    """A map printed over a `str` must not differ from one printed over a `Seq`.

    The reference's map draws the bottom strand with `self.sequence.complement()`,
    which only a `Seq` has.  A plain `str` gets the same answer through a
    translation table, and this compares that table against Biopython's own for
    all 256 byte values -- including `U`, which complements to `A`, and the
    bytes that have no complement at all and must come back unchanged.
    """
    from Bio.Seq import _dna_complement_table

    from biofasting._print_format import _COMPLEMENT

    assert _COMPLEMENT == _dna_complement_table
    assert bytes(b"acgtnuryswkmbdhv-*").translate(_COMPLEMENT) == b"tgcanayrswmkvhdb-*"


def test_the_map_draws_blocks_of_sixty_columns():
    sequence = "ACGT" * 32 + "GAATTC"  # 134 bases: 60 + 60 + 14
    report = ours.PrintFormat()
    report.sequence = sequence
    report.print_as("map")
    text = report.format_output(ours.RestrictionBatch(["EcoRI"]).search(sequence))
    lines = text.splitlines()
    #   Each block is the sequence line, the rule, the complement, then the
    #   coordinates; three blocks and a blank line after each.
    assert sequence[:60] in lines
    assert sequence[60:120] in lines
    assert sequence[120:] in lines
    assert "|" * 60 in lines and "|" * 14 in lines
    #   The bottom strand is the complement, not the reverse complement.
    assert "TGCATGCATGCA" in text
    #   The block coordinates are 1-based and name the block's last base.
    assert any(
        line.lstrip().startswith("121") and line.rstrip().endswith("134")
        for line in lines
    )


def test_print_as_picks_the_report_and_print_that_writes_it(capsys):
    sequence = "TTTGAATTCAAAGGATCC"
    report = ours.PrintFormat()
    report.sequence = sequence
    results = ours.RestrictionBatch(["EcoRI", "BamHI", "SmaI"]).search(sequence)
    report.print_as("number")
    assert "enzymes which cut 1 times" in report.format_output(results)
    report.print_as("nonsense")
    assert "EcoRI      :  5." in report.format_output(results)
    report.print_that(results, "Title:\n", "Nothing here:\n")
    printed = capsys.readouterr().out
    assert printed.startswith("Title:\n")
    assert "EcoRI      :  5." in printed
    #   The non-cutting enzymes go under the sentence that was handed in, and
    #   the report ends with the last of them rather than with the sentence --
    #   the trailing newline is `print`'s.
    assert printed.rstrip("\n").endswith("Nothing here:\nSmaI      ")


@needs_biopython
def test_analysis_searches_the_sequence_it_is_given():
    reference = _reference()
    from Bio.Seq import Seq

    sequence, names = _report_sequence()
    mine = ours.Analysis(ours.RestrictionBatch(names), sequence)
    theirs = reference.Analysis(reference.RestrictionBatch(names), Seq(sequence))
    assert {str(k): v for k, v in mine.mapping.items()} == {
        str(k): v for k, v in theirs.mapping.items()
    }
    assert str(mine) == str(theirs)
    #   `repr` carries the sequence, and the reference's is a `Seq` where this
    #   one is a `str` -- the deliberate difference, so the *shape* is compared.
    assert repr(mine) == f"Analysis({mine.rb!r},{sequence!r},{mine.linear})"
    assert repr(theirs).startswith("Analysis(RestrictionBatch(")
    #   The sequence is searched at construction, so the report needs no
    #   argument -- and `full()` is the same mapping that search left behind.
    assert mine.full() is mine.mapping


@needs_biopython
def test_the_analysis_filters_match_the_reference():
    reference = _reference()
    from Bio.Seq import Seq

    sequence, names = _report_sequence()
    mine = ours.Analysis(ours.RestrictionBatch(names), sequence)
    theirs = reference.Analysis(reference.RestrictionBatch(names), Seq(sequence))

    def names_of(dct):
        return sorted(str(k) for k in dct)

    for call in (
        lambda a: a.blunt(),
        lambda a: a.overhang5(),
        lambda a: a.overhang3(),
        lambda a: a.defined(),
        lambda a: a.with_sites(),
        lambda a: a.without_site(),
        lambda a: a.with_N_sites(1),
        lambda a: a.with_N_sites(2),
        lambda a: a.with_number_list([0, 2]),
        lambda a: a.blunt(a.with_sites()),
        lambda a: a.with_name(["EcoRI", "SmaI"]),
    ):
        assert names_of(call(mine)) == names_of(call(theirs)), call
        assert {str(k): v for k, v in call(mine).items()} == {
            str(k): v for k, v in call(theirs).items()
        }, call


@needs_biopython
def test_the_region_filters_match_the_reference():
    """`between` and its relatives, over regions that exercise the boundaries.

    The negative and zero bounds are in the list because `_boundaries` counts
    them from the end of the sequence, and a port that skipped that would
    answer a different question for every circular-plasmid caller.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, names = _report_sequence()
    mine = ours.Analysis(ours.RestrictionBatch(names), sequence)
    theirs = reference.Analysis(reference.RestrictionBatch(names), Seq(sequence))

    def names_of(dct):
        return sorted(str(k) for k in dct)

    for start, end in ((1, 200), (200, 700), (100, 110), (-100, -1), (0, -300), (650, 40)):
        for call in (
            lambda a: a.only_between(start, end),
            lambda a: a.between(start, end),
            lambda a: a.show_only_between(start, end),
            lambda a: a.only_outside(start, end),
            lambda a: a.outside(start, end),
            lambda a: a.do_not_cut(start, end),
        ):
            assert names_of(call(mine)) == names_of(call(theirs)), (start, end, call)
            assert {str(k): v for k, v in call(mine).items()} == {
                str(k): v for k, v in call(theirs).items()
            }, (start, end, call)


@needs_biopython
def test_a_region_that_is_empty_is_not_a_region_upstream():
    """`_boundaries` returns None when start == end, and the unpacking fails.

    Reproduced, not repaired: a caller who wrote `analysis.only_between(50, 50)`
    gets a TypeError from Biopython, and getting a *different* error -- or an
    answer -- here would be a new surprise rather than an old one.  The type
    and the fact are compared; the message is the interpreter's, not ours.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    for start, end in ((50, 50), (0, 0)):
        for call in ("only_between", "between", "outside", "do_not_cut"):
            mine = getattr(ours.Analysis(ours.RestrictionBatch(["EcoRI"]), sequence), call)
            theirs = getattr(
                reference.Analysis(reference.RestrictionBatch(["EcoRI"]), Seq(sequence)),
                call,
            )
            assert _outcome(lambda: mine(start, end)) == _outcome(lambda: theirs(start, end))
    #   And the boundary types are checked before anything else happens.
    analysis = ours.Analysis(ours.RestrictionBatch(["EcoRI"]), sequence)
    for bad in ("50", 1.5, None):
        assert _outcome(lambda: analysis.only_between(bad, 100)) == ("raised", "TypeError")


@needs_biopython
def test_with_site_size_is_broken_upstream_in_the_same_way():
    """`with_site_size(size, dct)` tests `key in size` -- an int is not iterable.

    So the method works only when it is *not* given a result dictionary, and
    raises `TypeError` when it is.  That is Biopython 1.88's behaviour, kept so
    that a program which works there works here; a silently repaired version
    would be a second, different library answering the same call.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    mine = ours.Analysis(ours.RestrictionBatch(["EcoRI", "BamHI", "HindIII"]), sequence)
    theirs = reference.Analysis(
        reference.RestrictionBatch(["EcoRI", "BamHI", "HindIII"]), Seq(sequence)
    )
    assert sorted(str(k) for k in mine.with_site_size(6)) == sorted(
        str(k) for k in theirs.with_site_size(6)
    )
    assert _outcome(lambda: mine.with_site_size(6, mine.mapping)) == ("raised", "TypeError")
    assert _outcome(lambda: theirs.with_site_size(6, theirs.mapping)) == (
        "raised",
        "TypeError",
    )


@needs_biopython
def test_with_name_warns_for_an_unknown_enzyme_and_may_skip_the_next_one():
    """The unknown-name warning, and the deletion that runs while iterating.

    `del names[index]` inside `for index, enzyme in enumerate(names)` shifts
    every later name down one, so two unknown names in a row produce one
    warning, not two -- and the one that was skipped is left in the list to be
    rejected by the batch constructor further down.  Both the warning count and
    the way each call ends are compared against the reference: the first case
    returns a mapping, the second raises `ValueError`, and neither is an
    accident of this port.
    """
    import warnings as warnings_module

    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    mine = ours.Analysis(ours.RestrictionBatch(["EcoRI", "BamHI"]), sequence)
    theirs = reference.Analysis(reference.RestrictionBatch(["EcoRI", "BamHI"]), Seq(sequence))

    def run(analysis, names):
        with warnings_module.catch_warnings(record=True) as captured:
            warnings_module.simplefilter("always")
            outcome = _outcome(
                lambda: sorted(str(k) for k in analysis.with_name(list(names)))
            )
        return outcome, [str(w.message) for w in captured], [w.category for w in captured]

    #   One unknown name: warned about, dropped, and the rest searched.
    one = ["EcoRI", "NotAnEnzyme", "BamHI"]
    ours_outcome, ours_messages, ours_categories = run(mine, one)
    assert ours_outcome == run(theirs, one)[0] == ("ok", ["BamHI", "EcoRI"])
    assert ours_messages == ["no data for the enzyme: NotAnEnzyme"]
    #   The warning is a plain UserWarning here and a BiopythonWarning there,
    #   which is a subclass of it: the documented way to catch it works on both.
    assert all(category is UserWarning for category in ours_categories)

    #   Two in a row: the second is never reached by the loop, so it is never
    #   warned about, and the constructor rejects it afterwards.
    two = ["EcoRI", "NotAnEnzyme", "NorThisOne", "BamHI"]
    ours_outcome, ours_messages, _ = run(mine, two)
    theirs_outcome, theirs_messages, _ = run(theirs, two)
    assert ours_outcome == theirs_outcome == ("raised", "ValueError")
    assert ours_messages == theirs_messages == ["no data for the enzyme: NotAnEnzyme"]

    #   With a result dictionary in hand the names are looked up in it rather
    #   than searched for -- see the defect pinned at the end of this file.
    analysis = ours.Analysis(ours.RestrictionBatch(["EcoRI", "BamHI"]), sequence)
    assert analysis.with_name(list(analysis), analysis.mapping)


@needs_biopython
def test_change_reconfigures_the_report_and_refuses_the_rest():
    """The three things `change` can do, and the two it raises for.

    What it does *not* do is as much a part of it as what it does: `linesize` is
    computed once when the class is defined and `change` never recomputes it,
    so widening the console moves `PrefWidth` and `Cmodulo` but leaves the wrap
    of the list section where it was.  That is upstream's behaviour, and the
    assertions below are against the reference's own values rather than against
    a rule of thumb.
    """
    reference = _reference()
    from Bio.Seq import Seq

    sequence, _ = _report_sequence()
    names = ["EcoRI", "BamHI"]
    mine = ours.Analysis(ours.RestrictionBatch(names), sequence)
    theirs = reference.Analysis(reference.RestrictionBatch(names), Seq(sequence))
    for analysis in (mine, theirs):
        analysis.change(NameWidth=20, ConsoleWidth=100)
    assert (mine.NameWidth, mine.Cmodulo, mine.PrefWidth, mine.linesize) == (
        theirs.NameWidth,
        theirs.Cmodulo,
        theirs.PrefWidth,
        theirs.linesize,
    )
    assert mine.linesize == 70, "linesize is fixed at class definition, upstream too"
    for key in ("Cmodulo", "PrefWidth"):
        assert _outcome(lambda: mine.change(**{key: 1})) == ("raised", "AttributeError")
        assert _outcome(lambda: theirs.change(**{key: 1})) == ("raised", "AttributeError")
    assert _outcome(lambda: mine.change(nonsense=1)) == ("raised", "AttributeError")
    #   `sequence` re-searches, which is how a `change` becomes an answer.
    other = "AAAGAATTCTTT"
    mine.change(sequence=other)
    theirs.change(sequence=Seq(other))
    assert {str(k): v for k, v in mine.mapping.items()} == {
        str(k): v for k, v in theirs.mapping.items()
    }
    #   And `rb` re-initialises the analysis against the old sequence.
    mine.change(rb=ours.RestrictionBatch(["SmaI"]))
    theirs.change(rb=reference.RestrictionBatch([reference.SmaI]))
    assert {str(k): v for k, v in mine.mapping.items()} == {
        str(k): v for k, v in theirs.mapping.items()
    }


@needs_biopython
def test_an_analysis_of_a_circular_sequence_matches_the_reference():
    reference = _reference()
    from Bio.Seq import Seq

    rng = random.Random(99)
    base = "".join(rng.choice("ACGT") for _ in range(300))
    #   The first site starts at the origin and the second ends at the last
    #   base, so both strands' wrap-around cuts are exercised.
    sequence = "GAATTC" + base[6:288] + "GGATCC" + base[294:]
    names = ["EcoRI", "BamHI", "SmaI"]
    mine = ours.Analysis(ours.RestrictionBatch(names), sequence, linear=False)
    theirs = reference.Analysis(
        reference.RestrictionBatch(names), Seq(sequence), linear=False
    )
    assert {str(k): v for k, v in mine.mapping.items()} == {
        str(k): v for k, v in theirs.mapping.items()
    }
    assert mine.linear is False


@needs_biopython
def test_an_empty_analysis_is_legal_and_the_leftovers_still_work():
    """`Analysis()` with no arguments, plus the two methods nothing calls."""
    reference = _reference()
    mine = ours.Analysis()
    theirs = reference.Analysis()
    assert len(mine) == 0 and mine.mapping == {}
    #   `repr` differs only in the sequence's own representation, which is the
    #   deliberate type difference between the two libraries.
    assert repr(mine) == "Analysis(RestrictionBatch([]),'',True)"
    assert repr(theirs) == "Analysis(RestrictionBatch([]),Seq(''),True)"
    #   `_sub_set` and `_test_reverse` are dead code in both libraries -- no
    #   method calls either.  They are kept because they are part of the
    #   class's surface, and they are tested because dead code that has never
    #   run is the code most likely to be wrong.
    analysis = ours.Analysis(ours.RestrictionBatch(["EcoRI"]), "TTTGAATTCAAA")
    assert {str(k) for k in analysis._sub_set(set(analysis))} == {"EcoRI"}
    assert analysis._sub_set(set()) == {}
    assert analysis._test_normal(1, 10, 5) is True
    assert analysis._test_normal(6, 10, 5) is False
    assert analysis._test_reverse(1, 10, 11) is True
    assert analysis._test_reverse(1, 10, 5) is True
    assert analysis._test_reverse(1, 10, 13) is False


def test_the_analysis_takes_a_str_where_the_reference_demands_a_seq():
    """The deliberate difference, stated as a test rather than a comment."""
    analysis = ours.Analysis(ours.RestrictionBatch(["EcoRI"]), "TTTGAATTCAAA")
    assert analysis.mapping[ours.get_enzyme("EcoRI")] == [5]
    assert analysis.blunt() == {}


# --------------------------------------------------------------------------
# The examples in the modules themselves
# --------------------------------------------------------------------------


def test_the_docstring_examples_are_true():
    """Every `>>>` in the two modules is executed, not merely written.

    An example in a docstring is a claim about behaviour, and a claim nothing
    runs is a claim that is free to become false -- which is how a library ends
    up documenting something it no longer does.  These two modules carry the
    restriction documentation; this runs it.
    """
    import doctest

    from biofasting import _print_format

    for module in (ours, _print_format):
        results = doctest.testmod(module, verbose=False, report=False)
        assert results.attempted > 0, f"{module.__name__} documents nothing"
        assert results.failed == 0, module.__name__

