"""Render the "where the time goes" split (bench/profiling.md) as a bar chart.

The point of the chart is one distinction: how much of each row is Python
bytecode inside Bio/** -- work a C core deletes -- against rows that are already
C-bound, where there is nothing to delete.  Shares are carried verbatim from the
table in bench/profiling.md.

Usage: .venv/bin/python bench/plot_profile.py  (writes misc/profile-share.png)
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent

# Same reference palette as bench/plot_g1.py.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BIO = "#1c5cab"    # blue step 550 -- deletable interpreter work
CBOUND = "#c3c2b7"  # neutral -- already C, nothing to delete

# (path as named in bench/profiling.md, cProfile share of time in Bio/**)
# Blue rows: Biopython implementations.  Grey rows: everything else, ~0%.
ROWS = [
    ("ops:biopython-revcomp", 63.9, True),
    ("fastq:biopython-seqio", 55.7, True),
    ("fasta:biopython-lowlevel", 55.4, True),
    ("fastq:biopython-lowlevel", 52.0, True),
    ("fasta:biopython-seqio", 49.5, True),
    ("fasta-random:biopython-index-random", 49.2, True),
    ("fastq-gz:biopython-seqio", 48.4, True),
    ("ops:biopython-gc", 37.5, True),
    ("fastq-gz:biopython-lowlevel", 36.1, True),
    ("fastq:pysam-htslib", 0.0, False),
    ("fastq-gz:pysam-htslib", 0.0, False),
    ("fastq-gz:zlib-raw-only", 0.0, False),
    ("fasta:pyfaidx-sequential", 0.0, False),
]
N = len(ROWS)
XMAX = 70


def main():
    fig = plt.figure(figsize=(9.0, 4.9), dpi=200, facecolor=SURFACE)
    fig.text(0.026, 0.968, "Where Biopython's time goes — share of profiled time "
             "in its own Python frames", fontsize=12.5, fontweight="bold",
             color=INK, va="top")
    fig.text(0.026, 0.908,
             "cProfile, one pass per path · bench/profiling.md · shares are "
             "meaningful, the profiled seconds are inflated by instrumentation",
             fontsize=8.5, color=MUTED, va="top")

    ax = fig.add_axes([0.315, 0.115, 0.50, 0.745], facecolor=SURFACE)
    ax.set_xlim(0, XMAX)
    ax.set_ylim(-0.7, N - 1 + 0.7)
    for tick in (0, 20, 40, 60):
        ax.axvline(tick, color=GRID, lw=0.8, zorder=0)
    ax.axvline(0, color=AXIS, lw=1.0, zorder=1)
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xticklabels(["0%", "20%", "40%", "60%"], fontsize=8, color=MUTED)
    ax.set_yticks([])
    ax.tick_params(axis="both", which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for k, (path, share, is_bio) in enumerate(ROWS):
        y = N - 1 - k
        colour = BIO if is_bio else CBOUND
        if share > 0:
            ax.barh(y, share, height=0.56, color=colour, zorder=3,
                    edgecolor=SURFACE, linewidth=0.8)
        else:
            # A zero-length bar is invisible; draw a short stub so the row reads
            # as "measured, and it came out at nothing" rather than as missing.
            ax.barh(y, 0.55, height=0.56, color=colour, zorder=3,
                    edgecolor=SURFACE, linewidth=0.8)
        ax.text(-0.014, y, path, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8, color=INK if is_bio else
                MUTED, clip_on=False)
        label = f"{share:.1f}%" if is_bio else "0%"
        ax.text(share + (1.4 if share > 0 else 1.0), y,
                label, ha="left", va="center", fontsize=8.5,
                fontweight="bold" if is_bio else "normal",
                color=INK if is_bio else MUTED)

    handles = [
        Line2D([], [], marker="s", ls="none", markersize=7, color=BIO),
        Line2D([], [], marker="s", ls="none", markersize=7, color=CBOUND),
    ]
    # Inside the axes, in the empty area right of the four 0% rows -- there the
    # plot has nothing to say, so the key can sit there without a collision.
    ax.legend(handles, ["Python bytecode in Bio/** — a C core deletes this",
                        "already C-bound — nothing to delete"],
              loc="lower right", bbox_to_anchor=(1.0, 0.012), frameon=False,
              fontsize=8.5, labelcolor=SECONDARY, borderaxespad=0,
              handlelength=1.0, handletextpad=0.5, ncol=1, labelspacing=0.35)

    fig.text(0.5, 0.045,
             "between a third and two thirds of every Biopython row is "
             "interpreter work, not algorithm — the gap is deletable",
             ha="center", fontsize=8.5, color=SECONDARY)
    fig.text(0.5, 0.014,
             "cProfile charges C-extension time to the calling Python frame, so "
             "the pysam rows read as 0% here; bench/profiling.md reads them "
             "with py-spy instead",
             ha="center", fontsize=7.5, color=MUTED)

    out = ROOT / "misc" / "profile-share.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
