'''
Make asymmetry-vs-momentum graphs from an existing results.md (no data files are read).

Usage:
    python3 plot_from_results.py                       # reads ~/Desktop/results.md
    python3 plot_from_results.py /path/to/results.md   # or give the path

Graphs are saved next to results.md: asymmetry_<WT id>.png for each pair, and
asymmetry_all_pairs.png with all six together.
'''
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = os.path.expanduser("~/Desktop/results.md")
# a path given on the command line overrides the default (arguments added by notebooks are ignored)
if len(sys.argv) > 1 and not sys.argv[1].startswith("-") and os.path.exists(sys.argv[1]):
    RESULTS = sys.argv[1]
OUT_DIR = os.path.dirname(os.path.abspath(RESULTS))
N_FILES = 10          # number of sub-samples (files)
FLAG_RATIO = 2.0


def parse(path):
    """Returns overall {wt_id: (A, error)} and bins {wt_id: [rows]} from results.md."""
    overall, bins, name, wt_of = {}, {}, {}, {}
    section = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("### Wild type vs. mutant"):
                section = "overall"
            elif line.startswith("### Asymmetry vs. momentum:"):
                title = line.split(":", 1)[1].strip()          # "E. coli WT / E. coli mutant"
                section = title
                bins[title] = []
            elif line.startswith("### "):
                section = None
            elif line.startswith("|") and section:
                cols = [c.strip() for c in line.strip("|").split("|")]
                if section == "overall" and " vs " in cols[0]:
                    wt = int(cols[0].split(" vs ")[0])
                    overall[wt] = (float(cols[4]), float(cols[5]))
                elif section != "overall" and re.match(r"^[\d.]+ - ", cols[0]):
                    low, high = [s.strip() for s in cols[0].split(" - ")]
                    bins[section].append({
                        "low": float(low), "high": None if high == "inf" else float(high),
                        "cells": float(cols[1].replace(",", "")),      # mean per file
                        "A": float(cols[2]), "err": float(cols[3]),
                    })
    return overall, bins


def main():
    overall, bins = parse(RESULTS)
    pairs = list(zip(overall.keys(), bins.keys()))     # same order in the file
    for wt, title in pairs:
        for r in bins[title]:
            counting = ((1 - r["A"] ** 2) / (r["cells"] * N_FILES)) ** 0.5
            r["ratio"] = r["err"] / counting
            r["flag"] = r["ratio"] > FLAG_RATIO

    def draw(ax, wt, title):
        rows = bins[title]
        am, ae = overall[wt]
        # y-range from the reliable bins and the overall band only
        ok = [r for r in rows if not r["flag"]]
        lo = min([r["A"] - r["err"] for r in ok] + [am - ae])
        hi = max([r["A"] + r["err"] for r in ok] + [am + ae])
        pad = 0.15 * (hi - lo)
        lo, hi = lo - pad, hi + pad
        ax.set_ylim(lo, hi)
        ax.axhline(0, color="grey", lw=1)
        ax.axhspan(am - ae, am + ae, color="#9ecae1", alpha=0.5, lw=0, label="overall A \u00b1 uncertainty")
        ax.axhline(am, color="#3182bd", ls="--", lw=1)
        for i, r in enumerate(rows):
            if r["flag"] and not (lo <= r["A"] <= hi):
                # value outside the plot: triangle at the edge with its value written out
                up = r["A"] > hi
                y = hi - 0.04 * (hi - lo) if up else lo + 0.04 * (hi - lo)
                ax.plot(i, y, marker="^" if up else "v", ms=8, color="black", mfc="white",
                        clip_on=False, zorder=5)
                ax.annotate(f"A={r['A']:.3g}\n\u00b1{r['err']:.2g}\n(off scale)", (i, y),
                            xytext=(0, -7 if up else 7), textcoords="offset points",
                            ha="center", va="top" if up else "bottom", fontsize=6.5)
            else:
                ax.errorbar(i, r["A"], yerr=r["err"], fmt="o", color="black", capsize=3,
                            mfc="white" if r["flag"] else "black")
        ax.errorbar([], [], yerr=[], fmt="o", color="black", mfc="black", label="momentum bin")
        ax.errorbar([], [], yerr=[], fmt="o", color="black", mfc="white",
                    label=f"unreliable bin (error > {FLAG_RATIO:g}\u00d7 counting error)")
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([(f">{r['low']:g}" if r["high"] is None else f"{r['low']:g}\u2013{r['high']:g}")
                            for r in rows], rotation=35, ha="right", fontsize=8)
        ax.set_xlabel("momentum p ($10^{-20}$ kg m/s), bins with roughly equal number of cells")
        ax.set_ylabel("asymmetry A")
        ax.set_title(title.replace(" / ", " vs "), fontsize=10)

    for wt, title in pairs:
        fig, ax = plt.subplots(figsize=(8, 5))
        draw(ax, wt, title)
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"asymmetry_{wt}.png"), dpi=200)
        plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, (wt, title) in zip(axes.flat, pairs):
        draw(ax, wt, title)
    axes.flat[0].legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "asymmetry_all_pairs.png"), dpi=150)
    plt.close(fig)
    print(f"Graphs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
