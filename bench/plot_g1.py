"""Render the Gate G1 ranking table (bench/targets.md) as a dumbbell chart.

Numbers are carried verbatim from the ranking table in bench/targets.md, which
is itself generated from bench/results/latest.json — so the image says the
same thing as the table it illustrates.

Usage: .venv/bin/python bench/plot_g1.py  (writes misc/g1-ranking.png)
"""

from pathlib import Path
from math import log10

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

ROOT = Path(__file__).resolve().parent.parent

# Reference palette (light mode). One hue, two shades for the dumbbell ends
# (Biopython = dark step 550, best alternative = step 300; the connector is a
# lighter track of the same ramp), per the dataviz method.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BIO = "#1c5cab"    # blue step 550
ALT = "#6da7ec"    # blue step 300
CONN = "#b7d3f6"   # blue step 150, track only

# (workload, detail, bio_s, bio_label, alt_s, alt_label, gap, lost)
# in table-row order (1..6), matching bench/targets.md; Gap and Stolen live
# in the columns right of each row, so no re-sorting here.
ROWS = [
    ("fastq-plain", "1M × 150 bp · 369 MB", 3.188, "3.19 s", 0.448,
     "0.45 s", "7.1×", "2.74 s"),
    ("fasta-random", "100 × 150 bp slices", 15.956, "15.96 s", 0.0003,
     "0.0003 s", "~50,000×", "15.96 s"),
    ("fastq-gzip", "1M reads, gzipped", 5.316, "5.32 s", 1.401, "1.40 s",
     "3.8×", "3.92 s"),
    ("ops-gc", "GC fraction, 100k", 0.697, "0.70 s", 0.107, "0.11 s",
     "6.5×", "0.59 s"),
    ("ops-revcomp", "reverse complement", 0.280, "0.28 s", 0.041, "0.04 s",
     "6.8×", "0.24 s"),
    ("fasta-scan", "100 Mbp, 5 records", 0.303, "0.30 s", 0.153, "0.15 s",
     "2.0×", "0.15 s"),
]
N = len(ROWS)
XLIM = (1e-4, 60)  # log decades; the 0.0003 s dot clears the left edge


def main():
    fig = plt.figure(figsize=(9.0, 3.7), dpi=200, facecolor=SURFACE)
    fig.text(0.026, 0.965, "Gate G1 ranking — Biopython vs. best available "
             "alternative", fontsize=12.5, fontweight="bold", color=INK,
             va="top")
    fig.text(0.026, 0.895,
             "median-of-5, one run of bench/results/latest.json · "
             "Intel i5-13420H · single-threaded · lower is better",
             fontsize=8.5, color=MUTED, va="top")

    ax = fig.add_axes([0.155, 0.165, 0.645, 0.60], facecolor=SURFACE)
    ax.set_xscale("log")
    ax.set_xlim(*XLIM)
    ax.set_ylim(-0.6, N - 1 + 0.6)
    for decade in (1e-4, 1e-3, 1e-2, 1e-1, 1, 10):
        ax.axvline(decade, color=GRID, lw=0.8, zorder=0)
    ax.axvline(XLIM[0], color=AXIS, lw=1.0, zorder=1)  # left spine is the axis

    decades = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10]
    ax.set_xticks([d for d in decades if d >= XLIM[0]])
    ax.set_xticklabels(["0.0001", "0.001", "0.01", "0.1", "1", "10"],
                       fontsize=8, color=MUTED)
    ax.set_yticks([])
    ax.tick_params(axis="both", which="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    trans = blended_transform_factory(ax.transAxes, ax.transData)

    def frac(x):  # x-position of a dot as an axes fraction
        lo, hi = XLIM
        return (log10(x) - log10(lo)) / (log10(hi) - log10(lo))

    def dot_label(x, text, y, fontsize, color, weight="normal", prefer="left"):
        """Time label beside its dot. Default side is left of the dot (right
        when the dot hugs the left axis); on short connectors the dark dot's
        label flips to the right of the dot so the two labels never meet."""
        side = frac(x)
        if prefer == "left" and side >= 0.15:
            ax.text(side - 0.030, y + 0.27, text, transform=trans, ha="right",
                    va="bottom", fontsize=fontsize, color=color,
                    fontweight=weight)
        elif prefer == "right" and side <= 0.90:
            ax.text(side + 0.030, y + 0.27, text, transform=trans, ha="left",
                    va="bottom", fontsize=fontsize, color=color,
                    fontweight=weight)
        else:
            other = side - 0.030 if side > 0.5 else side + 0.030
            ax.text(other, y + 0.27, text, transform=trans,
                    ha="right" if other < side else "left", va="bottom",
                    fontsize=fontsize, color=color, fontweight=weight)

    for k, (workload, detail, bio_s, bio_l, alt_s, alt_l, gap, lost) in enumerate(ROWS):
        y = N - 1 - k  # rank 1 at the top
        short_connector = frac(bio_s) - frac(alt_s) < 0.10
        ax.add_line(Line2D([alt_s, bio_s], [y, y], color=CONN, lw=1.4,
                           zorder=2, solid_capstyle="round"))
        for x, color in ((alt_s, ALT), (bio_s, BIO)):
            ax.scatter([x], [y], s=36, color=color, edgecolor=SURFACE,
                       linewidth=1.44, zorder=3)
        dot_label(alt_s, alt_l, y, 7.5, SECONDARY)
        dot_label(bio_s, bio_l, y, 8, INK, weight="bold",
                  prefer="right" if short_connector else "left")
        # row key, left of the plot; gap and lost seconds, right of it
        ax.text(-0.012, y + 0.16, workload, transform=trans, ha="right",
                va="center", fontsize=9, color=INK, clip_on=False)
        ax.text(-0.012, y - 0.18, detail, transform=trans, ha="right",
                va="center", fontsize=7.5, color=MUTED, clip_on=False)
        ax.text(1.025, y + 0.16, gap, transform=trans, ha="left", va="center",
                fontsize=9.5, fontweight="bold", color=INK, clip_on=False)
        ax.text(1.025, y - 0.18, f"lost {lost}", transform=trans, ha="left",
                va="center", fontsize=7.5, color=SECONDARY, clip_on=False)

    handles = [Line2D([], [], marker="o", ls="none", markersize=6,
                      markeredgecolor=SURFACE, markeredgewidth=1.0,
                      color=BIO),
               Line2D([], [], marker="o", ls="none", markersize=6,
                      markeredgecolor=SURFACE, markeredgewidth=1.0,
                      color=ALT)]
    fig.legend(handles, ["Biopython", "best available alternative"],
               loc="upper right", bbox_to_anchor=(0.975, 1.0), frameon=False,
               fontsize=8.5, labelcolor=SECONDARY, borderaxespad=0,
               handlelength=1.2, handletextpad=0.5, ncol=2, columnspacing=1.2)

    fig.text(0.5, 0.062, "seconds, log scale", ha="center", fontsize=8.5,
             color=MUTED)
    fig.text(0.5, 0.028,
             "fastq-gzip measures against the zlib decompression floor "
             "(no parser on top of it)",
             ha="center", fontsize=7.5, color=MUTED)

    out = ROOT / "misc" / "g1-ranking.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()