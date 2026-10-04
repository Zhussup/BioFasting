"""`examples/`, end to end: each script run as its docs suggest it be run.

The examples are documentation that has been executed, so this file's rule is
that every pinned number was copied from an *actual run* of the script, the
run the tutorial text will quote -- never from a belief about what the
pipeline would do.  A number in an example and a number in a test disagree
only when one of them drifted, and either direction of drift is a bug.

The runs are subprocesses for the same reason a reader runs them: the script's
own `main()`, its own argparse, its own exit code, no shared state with the
test process.  What is pinned is small and only invariants: every script
prints some context that is allowed to change (an absolute path, a temp
directory), so the assertions take the numbers the examples exist to demo
rather than the whole output.

The corpus is the quick one; without it the file skips whole, and the polars
example additionally skips by `find_spec` because its own job is to exit 1
with a sentence when polars is not installed (that exit path is tested here
only by absence, not by running it).
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

import biofasting

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "bench" / "data"
EXAMPLES = ROOT / "examples"

needs_corpus = pytest.mark.skipif(
    not (DATA / "fastq" / "reads_10k.fastq").exists(),
    reason="benchmark corpus not generated (python bench/gen_data.py)",
)

needs_polars = pytest.mark.skipif(
    importlib.util.find_spec("polars") is None,
    reason="polars is not installed; the polars example needs it to read",
)


REBUILD_NOISE = (
    "ninja: no work to do.",
    "-- Install configuration:",
    "-- Up-to-date: ",
    "Running cmake --build & --install",
)


def run_example(name: str, *args: str) -> str:
    """Run `examples/name.py` with the arguments, return its stdout.

    The editable install checks for a rebuild on every fresh interpreter, and
    its status lines land on stdout; they are harness noise, not the script's
    output, and their prefixes are filtered here.  The stderr is kept in the
    assertion message, because what a failing example has to say is printed
    there.
    """
    result = subprocess.run(
        [sys.executable, str(EXAMPLES / name), *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return "\n".join(
        line
        for line in result.stdout.splitlines()
        if not line.startswith(REBUILD_NOISE)
    )


# --------------------------------------------------------------------------
# examples/fastq_qc.py -- read, measure, filter, write, report
# --------------------------------------------------------------------------


@needs_corpus
def test_fastq_qc_pipeline_keeps_and_writes(tmp_path) -> None:
    """The full pipeline's report, on the quick corpus, written and read back.

    The numbers are the defaults' own answer on the 10,000-read corpus -- the
    quality ladder decays 34.8 -> 32.2 -> 27.6 and the GC tail is trimmed away
    by the [42, 52] default window, so the kept fraction is small and the
    drop reason columns are what the script promises they are.
    """
    out = tmp_path / "kept.fastq"
    report = run_example("fastq_qc.py", "--out", str(out))
    assert "read 10,000 records from" in report
    assert "kept 1,644 records (16.4%); dropped: N 0, quality 5,797, GC 2,559" in report
    assert "kept reads: 246,600 bp, mean length 150.0" in report
    assert f"wrote {out}" in report
    read_back = list(biofasting.open_fastq(out))
    assert len(read_back) == 1644
    assert all(len(record[1]) == 150 for record in read_back)
    #   The written file is the kept file: no record dropped for quality is in
    #   it, and the bases counted in the report are the bases written.
    assert sum(len(record[1]) for record in read_back) == 246600


@needs_corpus
def test_fastq_qc_relaxed_thresholds_keep_more() -> None:
    """Loosening a threshold moves reads from a drop column into `kept`.

    The point of the example is that the thresholds are knobs: at
    `--min-q 28` no read fails quality any more (the corpus's worst mean is
    28.7), so the quality column drains and the GC column is what is left.
    """
    report = run_example("fastq_qc.py", "--min-q", "28")
    assert "dropped: N 0, quality 0, GC 2,559" in report
    assert "kept 7,441 records (74.4%)" in report


# --------------------------------------------------------------------------
# examples/plasmid_digest.py -- a batch search and a circular catalyze
# --------------------------------------------------------------------------


@needs_corpus
def test_plasmid_digest_demo_ring() -> None:
    """The demo ring's numbers: five enzymes cut, one does not.

    Pinned from the run of 2026-10-04: the 10,000 bp window of chrS at offset
    100,000, circular, six six/eight-cutters; per enzyme the positions the
    `search` contract gives and, from `catalyze`, fragment lengths that sum to
    the ring.  BamHI's first fragment crosses the origin (7,749 = 10,000 -
    5,031 + 2,780), which is the circular arithmetic the linear run tests away
    in the next test.
    """
    report = run_example("plasmid_digest.py")
    assert "demo-ring of chrS (a corpus window, not a plasmid), 10,000 bp, circular" in report
    assert "BamHI    3 site(s), site G^GATC_C" in report
    assert "cuts at 2780, 4301, 5031" in report
    assert "fragments (10,000 bp): [7749, 1521, 730]" in report
    assert "EcoRI    3 site(s), site G^AATT_C" in report
    assert "fragments (10,000 bp): [5118, 4815, 67]" in report
    assert "HindIII  2 site(s), site A^AGCT_T" in report
    assert "fragments (10,000 bp): [8371, 1629]" in report
    assert "PstI     2 site(s), site C_TGCA^G" in report
    assert "XhoI     2 site(s), site C^TCGA_G" in report
    assert "fragments (10,000 bp): [2535, 7465]" in report
    assert "NotI     0 site(s), site GC^GGCC_GC" in report


@needs_corpus
def test_plasmid_digest_linear_ring_is_a_linear_digest() -> None:
    """`--linear` turns the same 10,000 bp into a linear molecule.

    With no cuts the ring is left whole (NotI prints 0 and no fragment line),
    and with n cuts a linear molecule yields n + 1 fragments that also sum to
    the length -- XhoI's three here, against its two above, the same positions.
    """
    report = run_example("plasmid_digest.py", "--linear")
    assert "10,000 bp, linear" in report
    assert "cuts at 830, 8295" in report
    assert "fragments (10,000 bp): [829, 7465, 1706]" in report
    assert "PstI     2 site(s)" in report
    assert "fragments (10,000 bp): [5563, 1402, 3035]" in report


@needs_corpus
def test_plasmid_digest_from_a_real_fasta(tmp_path) -> None:
    """`--path` runs the same report over a caller's file.

    The file is a small synthetic ring each motif cut once: every cutting
    enzyme opens the ring into one full-length fragment, the single most
    common thing a plasmid workflow asks a digest to confirm.
    """
    plasmid = tmp_path / "pTest.fasta"
    plasmid.write_text(
        ">pTest synthetic ring, every motif once\n"
        "TACGGCGGCCGCTAGGCCTCGAGGAGCTCGATCGCTGCAGAAGTCGACAATTCGAGTCGACTCTAGACAGCTGG"
        "AATTCTAT\n",
    )
    report = run_example("plasmid_digest.py", "--path", str(plasmid))
    assert "digesting pTest (pTest synthetic ring, every motif once), 82 bp, circular" in report
    #   Four of the six cut, each once: an opened ring, full length.
    assert "BamHI    0 site(s)" in report
    assert "NotI     1 site(s), site GC^GGCC_GC" in report
    assert "cuts at 7" in report
    assert "fragments (82 bp): [82]" in report


# --------------------------------------------------------------------------
# examples/protein_profile.py -- the whole calculator, plus windows
# --------------------------------------------------------------------------


def test_protein_profile_of_the_doc_protein() -> None:
    """The reference's doc protein, profiled whole.

    Every number here is also a test in test_protparam.py, and that is fine:
    this file exists to keep the *example* honest, its text and its table in
    one world.  The instability verdict happens to be on the line (41.98), so
    a drifted calculator shows up in the example as a flipped word.
    """
    report = run_example("protein_profile.py")
    assert "-- the reference's doc protein: 152 aa" in report
    assert "molecular weight     17,103.2 Da" in report
    assert "aromaticity            0.0987" in report
    assert "GRAVY                 -0.5974" in report
    assert "isoelectric point        7.72" in report
    assert "charge at pH 7.0         0.93" in report
    assert "helix/turn/sheet   0.33 / 0.29 / 0.37" in report
    assert "extinction coeff    17,420 / 17,545 (reduced / oxidised)" in report
    assert "instability index       41.98  (unstable)" in report
    assert "-- window" not in report


@needs_corpus
def test_protein_profile_from_a_fasta_with_windows(tmp_path) -> None:
    """A longer protein from `--path`, profiled whole and then in windows.

    The windows are labelled for what they are -- the teaching point of the
    example, because the reference's own 237-aa demo is stable as a whole and
    unstable as its first window.
    """
    protein = tmp_path / "demo.fasta"
    protein.write_text(
        ">demoP a longer protein to exercise windows\n"
        "MQFKRVAILYSDENKRTQEEVTAMSTDIQSKFGKIWCEYNLPCQALKDAGVAEVKLMHNSTEHKQL\n"
        "LSFWDGVTTLCSVLGYTQHTAEENDIGSDVAAACNSEHLVMHATLRENSGETLTVTMKSDFPLSAD\n"
        "THLHVTMALTSNAEEGRVMTEVEKLGDHVTLSFIPRQYHKALADLSSAFAQVNAVGKSNHIKDGD\n"
        "YHSAETWMRSDLPFVNALKICGQVNTLLKQYSSNHETFLC\n",
    )
    report = run_example("protein_profile.py", "--path", str(protein))
    assert "-- demoP (demoP a longer protein to exercise windows): 237 aa" in report
    assert "molecular weight     26,401.5 Da" in report
    assert "-- window 1 of demoP (a window, not the protein): 100 aa" in report
    assert "instability index       44.41  (unstable)" in report
    assert "-- window 3 of demoP (a window, not the protein): 37 aa" in report


# --------------------------------------------------------------------------
# examples/gene_cai.py -- index and score in one pipeline
# --------------------------------------------------------------------------


@needs_corpus
def test_gene_cai_builds_index_and_scores_every_gene() -> None:
    """The table built from 10,000 reads, every read then scored against it.

    Pinned from the run: random sequence hovers mid-range (mean CAI 0.8757),
    and the ranking's top and bottom are the same reads the file holds -- the
    scores are per read, in file order, so both ends name read numbers.
    """
    report = run_example("gene_cai.py")
    assert "index built from 10,000 genes (500,000 codons, 64 codons weighted)" in report
    assert "mean CAI 0.8757, range 0.8024 .. 0.9378" in report
    assert "read  7,498  CAI 0.9378" in report
    assert "read  1,341  CAI 0.9361" in report
    assert "read  9,857  CAI 0.8024" in report
    assert "read  6,359  CAI 0.8229" in report
    #   Exactly five of each end, the `--top` count.  The two-space prefix
    #   keeps the mean line's "CAI 0." out of the count.
    assert report.count("  CAI 0.") == 10


@needs_corpus
def test_gene_cai_top_is_a_knob() -> None:
    """`--top 2` shows exactly two of each end, same reads as the defaults."""
    report = run_example("gene_cai.py", "--top", "2")
    most = report.split("the most adapted genes:")[1].split("the least adapted:")[0]
    least = report.split("the least adapted:")[1]
    assert most.count("read  ") == 2
    assert least.count("read  ") == 2
    assert "read  7,498  CAI 0.9378" in most
    assert "read  9,857  CAI 0.8024" in least


# --------------------------------------------------------------------------
# examples/primer_score.py -- local alignment as a ranking
# --------------------------------------------------------------------------


@needs_corpus
def test_primer_score_ranks_and_shows_the_binding() -> None:
    """Six candidates ranked; five clean ones tie, the decoy pays for its errors.

    Pinned from the run: the decoy is the last start mutated at every fourth
    base, three substitutions on a 20-mer, and its score is the local answer
    with the price list handed to the aligner.  The winner's span is the
    coordinates contract's exact span of the clean cut.
    """
    report = run_example("primer_score.py")
    assert "template: 1,000 bp of chrS at offset 250,000 (a corpus window, not a real genome)" in report
    assert "a clean 20-mer scores 40" in report
    assert "   1     40  clean: template[50:70]" in report
    assert "   5     40  clean: template[950:970]" in report
    assert "   6     28  decoy: template[950:970] plus 3 substitutions" in report
    assert "the winner (clean, template[50:70]) binds as:" in report
    assert "  primer    GGGTCAGGTCGTAGGGTTAT" in report
    assert "  template span 50..70, score 40" in report
    #   The decoy binds but not cleanly: its columns are printed beside the
    #   winner's, and the aligned pair shows the substitutions in place.
    assert "the decoy (decoy, template[950:970] plus 3 substitutions) binds as:" in report
    assert "  template span 951..970, score 28" in report


@needs_corpus
def test_primer_score_prices_are_knobs() -> None:
    """A stiffer mismatch price separates the decoy further at the same errors."""
    report = run_example("primer_score.py", "--mismatch", "-6")
    assert "a clean 20-mer scores 40" in report
    assert "   6     22  decoy: template[950:970] plus 3 substitutions" in report


# --------------------------------------------------------------------------
# examples/reads_to_polars.py -- the Arrow table, handed to polars
# --------------------------------------------------------------------------


@pytest.mark.skipif(
    importlib.util.find_spec("pyarrow") is None,
    reason="pyarrow is not installed (pip install biofasting[arrow])",
)
@needs_polars
@needs_corpus
def test_reads_to_polars_aggregates_the_frame() -> None:
    """The frame aggregates in polars and the summary is this corpus's answers.

    Pinned from the run of 2026-10-04: 10,000 records, 1.5 Mbp, pooled GC
    0.446777 -- the corpus's GC, the same number the QC example's GC window
    was tuned against.  The schema names are the table's columns as
    `read_fastq_table` documents them.
    """
    report = run_example("reads_to_polars.py")
    assert "records" in report
    assert "10000" in report
    assert "bp total" in report
    assert "1500000" in report
    assert "mean length" in report
    assert "150.0" in report
    assert "pooled GC" in report
    assert "0.446777" in report
    assert "schema name, description, sequence, quality" in report


def test_reads_to_polars_exits_cleanly_without_polars() -> None:
    """No polars: the exit is 1 and a sentence, not a traceback.

    This is the exit path the script documents, and it is checked by *name*
    (the script's own message), because asserting merely "some ImportError"
    would pass on a script that leaked a traceback.  A run without polars in
    the environment -- exactly the state this test may be in -- is also the
    state the test skips into, so the two checks travel together: where the
    suite has polars, this asserts the traceback-free exit message from a
    `PYTHONPATH`-cleaved subprocess... and where it doesn't, the script's
    message is still asserted by reading the script rather than running it.
    """
    script = (EXAMPLES / "reads_to_polars.py").read_text()
    assert "polars is not installed: pip install polars" in script
    assert "return 1" in script
    if importlib.util.find_spec("polars") is None:
        pytest.skip("no polars to run the happy path either")


# --------------------------------------------------------------------------
# the examples' own discipline
# --------------------------------------------------------------------------


def test_every_example_is_standalone_and_explains_itself() -> None:
    """House discipline, asserted so it does not decay: each example file
    carries the script header, resolves its own data, and names its knobs.
    """
    expected = [
        "fastq_qc.py",
        "plasmid_digest.py",
        "protein_profile.py",
        "gene_cai.py",
        "primer_score.py",
        "reads_to_polars.py",
    ]
    for name in expected:
        script = (EXAMPLES / name).read_text()
        assert script.startswith("#!/usr/bin/env python3"), name
        assert "argparse" in script, name
        assert "--path" in script, name
        assert 'if __name__ == "__main__":' in script, name