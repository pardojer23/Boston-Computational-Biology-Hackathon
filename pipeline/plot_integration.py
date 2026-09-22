#!/usr/bin/env python
"""Module 4 figure — the one summary plot.

    python pipeline/plot_integration.py

Reads ``results/integration_module/{redundancy_scores,weight_sensitivity}.csv``
and writes ``redundancy_summary.png`` beside them. Kept out of
``soy_globin_integration.py`` so the scoring logic has no matplotlib dependency,
and committed so the figure has a generator (every other figure in this repo
was produced ad hoc, which is issue 7 in the README's caveat list).

The figure makes one claim: within the leghemoglobin clade, redundancy is set by
expression dose, because molecular similarity is saturated.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

#: Two focal pairs trade first place across the weight sweep; they carry colour.
#: Everything else is greyscale by visual weight, so the figure uses two hues.
C_TOP_E = "#1f6fb2"   # first at the chosen weighting (expression-led)
C_TOP_M = "#d1622b"   # first when the molecular term dominates
C_LB = "#3a3a3a"      # the other four focal pairs
C_OTHER = "#b9b9b9"   # pairs involving a non-focal member
GREY = "#6e6e6e"

SIZES = (8, 7, 6)  # base / annotation / tick


def style() -> None:
    base, ann, tick = SIZES
    mpl.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300,
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "legend.fontsize": ann, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.titlelocation": "left", "axes.titlepad": 5.0,
        "legend.frameon": False, "figure.constrained_layout.use": False,
    })


def pair_label(r) -> str:
    return f"{r.symbol_a}\u2013{r.symbol_b}"


#: Pair class for the focal clade. Was the literal "Lb-Lb"; the integration
#: module now derives it from focal membership rather than the family name.
FOCAL_CLASS = "focal-focal"


def load(scores: Path, sweep: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the two integration tables by explicit path.

    Previously took a results *directory* and assembled the paths itself, which
    meant the figure could silently plot a different run's tables than the ones
    the workflow declared as its inputs.
    """
    return pd.read_csv(scores), pd.read_csv(sweep)


def panel_a(ax, sc: pd.DataFrame, alpha: float, beta: float, focal: dict[str, str]) -> None:
    """M versus E for all pairs, with iso-R contours of the gated score."""
    # Iso-R contours: R = M**alpha * E**beta  =>  M = (R / E**beta)**(1/alpha)
    e = np.linspace(0.01, 1.02, 400)
    for rv in (0.2, 0.4, 0.6, 0.8):
        m = (rv / e**beta) ** (1.0 / alpha)
        ok = (m >= 0) & (m <= 1.05)
        ax.plot(e[ok], m[ok], lw=0.5, color="#dcdcdc", zorder=0)
        # Label each contour at the right-hand edge. Not because the curves are
        # furthest apart there — within the visible window they are closest at
        # E = 1 — but because it is the one x where all four are guaranteed
        # inside the axes and clear of the data, which sits at low E or high M.
        # The label heights there are M = R**(1/alpha), i.e. 0.018 / 0.101 /
        # 0.279 / 0.572, spaced 0.08-0.29 apart: ample for 6 pt text.
        # An unlabelled contour is just a decorative curve.
        m_edge = rv ** (1.0 / alpha)
        if 0.0 <= m_edge <= 1.0:
            ax.text(0.985, m_edge + 0.012, f"$R$={rv:g}", fontsize=SIZES[2],
                    color=GREY, ha="right", va="bottom", zorder=1)

    lb = sc[sc.pair_class == FOCAL_CLASS]
    other = sc[sc.pair_class != FOCAL_CLASS]

    # Non-focal pairs first, at low visual weight; shape carries duplication mode.
    for mode, mk in (("tandem", "o"), ("dispersed", "s")):
        s = other[other.duplication_mode == mode]
        ax.scatter(s.E_family, s.M_family, s=13, marker=mk, facecolor=C_OTHER,
                   edgecolor="white", linewidth=0.3, zorder=2)

    for r in lb.itertuples():
        key = pair_label(r)
        col = focal.get(key, C_LB)
        mk = "o" if r.duplication_mode == "tandem" else "s"
        ax.scatter(r.E_family, r.M_family, s=34, marker=mk, facecolor=col,
                   edgecolor="white", linewidth=0.5, zorder=4)

    # M is saturated across the focal pairs (0.91-1.00), so all six points sit in
    # a narrow band at the top and fixed offsets collide. Labels drop straight
    # down into the empty mid-field on leaders of increasing length, ordered by E
    # so no two leaders cross: each label sits directly beneath its own point.
    ladder = lb.sort_values("E_family", ascending=False).reset_index(drop=True)
    base = float(lb.M_family.max())
    for i, r in enumerate(ladder.itertuples()):
        key = pair_label(r)
        # Depths are measured from a common baseline, not from each point's own
        # M: M varies by 0.09 across these pairs, enough to reorder the rungs and
        # collide two labels if each is offset from its own marker.
        y = base - (0.075 + 0.058 * i)
        ax.plot([r.E_family, r.E_family], [r.M_family - 0.012, y + 0.022],
                lw=0.4, color=focal.get(key, C_LB), alpha=0.7, zorder=3)
        ax.text(r.E_family, y, key, fontsize=SIZES[1], style="italic",
                color=focal.get(key, C_LB), ha="center", va="center", zorder=5)

    # Name the outgroup strip once rather than labelling 15 points. No leader:
    # they are the only marks at E ~ 0, so the text beside them is unambiguous,
    # and a leader reaching across the panel would cross the label ladder.
    ax.text(0.055, 0.34, "pairs with a\nnon-focal globin", fontsize=SIZES[1],
            color=GREY, ha="left", va="center")

    ax.set_xlabel("Expression term $E$  (tissue overlap \u00d7 dose ratio)")
    ax.set_ylabel("Molecular term $M$  (sequence, fold, pocket)")
    ax.set_title("Molecular similarity is saturated; expression dose is not")
    ax.set_xlim(0.0, 1.04)
    ax.set_ylim(0.0, 1.08)
    ax.set_xticks(np.arange(0.0, 1.01, 0.25))
    ax.set_yticks(np.arange(0.0, 1.01, 0.25))
    ax.text(0.015, 0.985, f"n = {len(sc)} pairs", transform=ax.transAxes,
            fontsize=SIZES[1], color=GREY, ha="left", va="top")

    handles = [
        plt.Line2D([], [], marker="o", ls="", ms=4, mfc=C_LB, mec="white",
                   mew=0.4, label="tandem"),
        plt.Line2D([], [], marker="s", ls="", ms=4, mfc=C_LB, mec="white",
                   mew=0.4, label="dispersed"),
    ]
    ax.legend(handles=handles, loc="lower left", handletextpad=0.4,
              borderaxespad=0.3, labelspacing=0.3)


def panel_b(ax, sc: pd.DataFrame, alpha: float, focal: dict[str, str],
            top_e_key: str, top_m_key: str) -> None:
    """Every focal pair's score as the weight on M sweeps 0 to 1."""
    lb = sc[sc.pair_class == FOCAL_CLASS]
    grid = np.linspace(0.0, 1.0, 101)
    curves = {}
    for r in lb.itertuples():
        key = pair_label(r)
        curves[key] = r.M_family**grid * r.E_family**(1.0 - grid)

    order = sorted(curves, key=lambda k: -curves[k][-1])
    for key in order:
        col = focal.get(key, C_LB)
        lw = 1.3 if key in focal else 0.7
        ax.plot(grid, curves[key], color=col, lw=lw,
                alpha=1.0 if key in focal else 0.55, zorder=3 if key in focal else 2)

    # The crossover: where the first-ranked focal pair changes.
    top = np.array([max(curves, key=lambda k: curves[k][i]) for i in range(len(grid))])
    flip = np.where(top[1:] != top[:-1])[0]
    for i in flip:
        x = 0.5 * (grid[i] + grid[i + 1])
        ax.axvline(x, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)
        ax.text(x + 0.02, 0.50, f"rank 1 changes\nat \u03b1 \u2248 {x:.2f}",
                fontsize=SIZES[1], color=GREY, ha="left", va="center")

    ax.axvline(alpha, color="black", lw=0.7, zorder=1)
    ax.text(alpha - 0.02, 0.02, f"chosen \u03b1 = {alpha:g}", fontsize=SIZES[1],
            color="black", ha="right", va="bottom")

    # Only the two colour-threaded curves are labelled; they are the pairs named
    # in panel a, and colour carries the cross-reference. Labelling all six here
    # is what the right-hand convergence makes illegible.
    for key, at in ((top_e_key, 0.16), (top_m_key, 0.84)):
        j = int(at * (len(grid) - 1))
        ax.annotate(key, (grid[j], curves[key][j]), textcoords="offset points",
                    xytext=(0, 6), fontsize=SIZES[1], style="italic",
                    color=focal[key], ha="center", va="bottom")
    ax.text(0.97, 0.06, "four other focal pairs", transform=ax.transAxes,
            fontsize=SIZES[1], color=C_LB, ha="right", va="bottom")

    spread0 = max(c[0] for c in curves.values()) - min(c[0] for c in curves.values())
    spread1 = max(c[-1] for c in curves.values()) - min(c[-1] for c in curves.values())
    ax.set_xlabel("$\\alpha$ — weight on the molecular term")
    ax.set_ylabel("Redundancy score $R$")
    ax.set_title("Weighting decides which pair ranks first")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.12)
    ax.set_xticks(np.arange(0.0, 1.01, 0.25))
    ax.set_yticks(np.arange(0.0, 1.01, 0.25))
    ax.text(0.03, 0.99, f"spread across the six pairs:\n{spread0:.2f} at \u03b1=0"
                        f"  \u2192  {spread1:.2f} at \u03b1=1",
            transform=ax.transAxes, fontsize=SIZES[1], color=GREY,
            ha="left", va="top")


def bbox_report(fig) -> list[str]:
    """Geometric check: visible text boxes that overlap each other or a spine."""
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    out = []
    for i, (ta, ba) in enumerate(texts):
        for tb, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                out.append(f"text/text: {ta.get_text()!r} x {tb.get_text()!r}")
    for t, bt in texts:
        for s, bs in spines:
            if bt.overlaps(bs) and t not in ticks.get(s.axes, ()):
                out.append(f"text/spine: {t.get_text()!r}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None, help="path to config.yaml")
    ap.add_argument("--scores", default=None, help="redundancy_scores.csv")
    ap.add_argument("--sweep", default=None, help="weight_sensitivity.csv")
    ap.add_argument("--results", default=str(ROOT / "results"),
                    help="fallback root when --scores/--sweep are not given")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    results = Path(args.results)
    scores = Path(args.scores) if args.scores else \
        results / "integration_module" / "redundancy_scores.csv"
    sweep = Path(args.sweep) if args.sweep else \
        results / "integration_module" / "weight_sensitivity.csv"
    out = Path(args.out) if args.out else \
        results / "integration_module" / "redundancy_summary.png"
    sc, sw = load(scores, sweep)
    alpha = float(sc.alpha_M.iloc[0])
    beta = float(sc.beta_E.iloc[0])

    lb = sc[sc.pair_class == FOCAL_CLASS]
    top_e = pair_label(lb.sort_values("R_family", ascending=False).iloc[0])
    top_m = pair_label(lb.sort_values("M_family", ascending=False).iloc[0])
    focal = {top_e: C_TOP_E, top_m: C_TOP_M}

    style()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1))
    panel_a(axes[0], sc, alpha, beta, focal)
    panel_b(axes[1], sc, alpha, focal, top_e, top_m)
    fig.subplots_adjust(left=0.075, right=0.90, bottom=0.145, top=0.88, wspace=0.42)
    for ax, letter in zip(axes, "ab"):
        ax.text(-0.135, 1.10, letter, transform=ax.transAxes, fontsize=SIZES[0] + 2,
                fontweight="bold", va="top", ha="left")

    fig.savefig(out)
    findings = bbox_report(fig)
    print(f"wrote {out}")
    print(f"colour: {top_e} = expression-led first place; {top_m} = molecular-led")
    if findings:
        print(f"bbox check: {len(findings)} overlap(s)")
        for f in findings:
            print("  " + f)
    else:
        print("bbox check: clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
