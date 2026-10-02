"""Render the README "Why this exists" benchmark table as a chart.

Reads nothing — the numbers are the draft single-run table from README.md,
kept here verbatim so the image stays in sync with the text it illustrates.

Usage: .venv/bin/python bench/plot_whythis.py  (writes images/why-this-exists.png)
"""

from pathlib import Path

import matplotlib
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Patch

matplotlib.use("Agg")

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent

# Token values from the dataviz reference palette (light mode). Categorical
# slots 1-3 in fixed order; aqua is sub-3:1 on the surface, which is why every
# bar carries a direct value label and the README keeps the table.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

FAMILY_COLOR = {
    "Biopython": "#2a78d6",   # slot 1 blue
    "pure Python": "#eb6834", # slot 2 orange
    "pysam": "#1baf7a",       # slot 3 aqua
}

# (approach family, panel, y-tick label, time [s], throughput label or None)
ROWS = [
    ("Biopython",   "parse",    "Bio.SeqIO.parse (1.88)",     0.31, "319k rec/s"),
    ("pure Python", "parse",    "naive pure-Python parser",   0.03, "3.42M rec/s"),
    ("pysam",       "parse",    "pysam.FastxFile (C/htslib)", 0.05, "2.21M rec/s"),
    ("Biopython",   "revcomp",  "SeqRecord loop",             0.35, None),
    ("pure Python", "revcomp",  "str.translate",              0.06, None),
]
PANELS = ["parse", "revcomp"]
PANEL_TITLE = {"parse": "Parse 100k reads", "revcomp": "Reverse complement (same 100k)"}
SLOT = 3  # both panels get 3 vertical slots so bar thickness matches

XLIM = 0.42  # headroom for outside labels on the short bars


def rounded_right_bar(ax, y, width, color, half_height, radius_px, ax_px, slot_px):
    """Bar with a 4px-round data end, square at the shared x=0 baseline."""
    rx = radius_px / ax_px * XLIM
    ry = radius_px / slot_px * 1.0
    y0, y1 = y - half_height, y + half_height
    rx = min(rx, width)  # no capsule wider than the bar itself
    verts = [
        (0, y0), (width - rx, y0), (width, y0), (width, y0 + ry),
        (width, y1 - ry), (width, y1), (width - rx, y1), (0, y1), (0, y0),
    ]
    codes = [
        MplPath.MOVETO, MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3,
        MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3, MplPath.LINETO,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color,
                           edgecolor="none", zorder=2))


def place_label(fig, ax, text_obj, width, y):
    """Put the value label inside the bar end when it fits there, else at the tip."""
    fig.canvas.draw()
    px = text_obj.get_window_extent(fig.canvas.get_renderer()).width
    bar_px = width / XLIM * ax.get_window_extent().width
    text_obj.remove()
    if bar_px >= px + 2 * 12:  # fits with padding on both sides -> inside tip
        ax.text(width - 0.006, y, text_obj.get_text(), ha="right", va="center",
                fontsize=9, fontweight="bold", color="white", zorder=3)
    else:
        ax.text(width + 0.008, y, text_obj.get_text(), ha="left", va="center",
                fontsize=9, fontweight="bold", color=INK, zorder=3)


def main():
    fig = plt.figure(figsize=(9.0, 3.4), dpi=200, facecolor=SURFACE)
    fig.text(0.026, 0.955, "Biopython vs. alternatives — 100k FASTQ reads",
             fontsize=12.5, fontweight="bold", color=INK, va="top")
    fig.text(0.026, 0.885,
             "Draft single-run benchmark · 30 MB FASTQ · Intel i5-13420H · "
             "single thread · time, lower is better",
             fontsize=8.5, color=MUTED, va="top")

    axes = {}
    x0, w, gutter = 0.175, 0.345, 0.100
    for i, panel in enumerate(PANELS):
        ax = fig.add_axes([x0 + i * (w + gutter), 0.155, w, 0.52], facecolor=SURFACE)
        ax.set_xlim(0, XLIM)
        ax.set_ylim(-0.5, SLOT - 0.5)
        for t in (0.1, 0.2, 0.3, 0.4):  # recessive hairline grid behind the bars
            ax.axvline(t, color=GRID, lw=0.8, zorder=0)
        ax.axvline(0, color=AXIS, lw=1.0, zorder=1)  # single shared baseline
        ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
        ax.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4"], fontsize=8,
                           color=MUTED)
        ax.set_yticks([])  # approach names are drawn manually, left of the axis
        ax.tick_params(axis="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(PANEL_TITLE[panel], loc="left", fontsize=10,
                     fontweight="bold", color=INK, pad=10)
        axes[panel] = ax

    for panel, ax in axes.items():
        fig.canvas.draw()  # get real pixel geometry for the 4px data-end radius
        ext = ax.get_window_extent()
        ax_px, slot_px = ext.width, ext.height / SLOT
        rows = [r for r in ROWS if r[1] == panel]
        for k, (family, _, label, t, thr) in enumerate(rows):
            y = SLOT - 1 - k  # table order, top row first
            ax.text(-0.008, y, label, ha="right", va="center",
                    transform=matplotlib.transforms.blended_transform_factory(
                        ax.transAxes, ax.transData),
                    fontsize=9, color=SECONDARY, clip_on=False)
            rounded_right_bar(ax, y, t, FAMILY_COLOR[family], 0.17, 8,
                              ax_px, slot_px)
            place_label(fig, ax, ax.text(0, y, f"{t:g} s" if thr is None else
                                         f"{t:g} s · {thr}",
                                         fontsize=9, fontweight="bold"), t, y)

    handles = [Patch(facecolor=c, edgecolor="none") for c in FAMILY_COLOR.values()]
    fig.legend(handles, FAMILY_COLOR.keys(), loc="upper right",
               bbox_to_anchor=(0.975, 1.0), frameon=False, fontsize=8.5,
               labelcolor=SECONDARY, handlelength=1.1, handleheight=1.1,
               borderaxespad=0, ncol=3, columnspacing=1.2, handletextpad=0.5)
    fig.text(0.5, 0.03, "Time (seconds)", ha="center", fontsize=8.5, color=MUTED)

    out = ROOT / "misc" / "why-this-exists.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()