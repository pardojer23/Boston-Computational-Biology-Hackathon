"""
Interpret HyPhy aBSREL branch-site test results (modal_absrel.py output) for
the 4-taxon soybean leghemoglobin quartet in terms of paralog redundancy:
plot the tree with branches colored by fitted omega (dN/dS), and a bar chart
comparing terminal-branch omega against each paralog's expression-detection
status in GSE270392 (see ../../gse270392_leghemoglobin_redundancy/).

IMPORTANT topology note: HyPhy analyzes the tree unrooted. For this 4-taxon
input, that means the actual tested tree is a TRIFURCATING root with three
branches: Lbc1 (terminal), Lba (terminal), and "Node4" (internal, leading to
the (Lbc3, Lbc2) clade) -- read directly from the aBSREL result's own
`input.trees` field, not assumed from the rooted guide tree we supplied.
There is no separate "Lbc1-Lba" branch in the tested tree; Node4 is the
ancestral branch of the (Lbc3, Lbc2) clade specifically.

Usage:
    python plot_absrel_selection.py <absrel_modal_result.json> <output.png>
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch

TIP_INFO = {
    "Glycine_max__3847_0_006217": {"label": "Lbc1", "detected": False},
    "Glycine_max__3847_0_0065fa": {"label": "Lba",  "detected": False},
    "Glycine_max__3847_0_0065e0": {"label": "Lbc3", "detected": True},
    "Glycine_max__3847_0_00ca2d": {"label": "Lbc2", "detected": True},
}


def main():
    result_path, out_path = sys.argv[1], sys.argv[2]
    with open(result_path) as f:
        modal_out = json.load(f)
    absrel = json.loads(modal_out["absrel_result"])

    # Confirm topology directly from HyPhy's own recorded input tree rather
    # than assuming it matches the rooted guide tree we supplied.
    tested_tree_str = absrel["input"]["trees"]["0"]
    assert "Node4" in tested_tree_str, "expected topology changed -- inspect input.trees before plotting"

    branch_attrs = absrel["branch attributes"]["0"]
    omega = {attrs_label(b): attrs["Baseline MG94xREV omega ratio"] for b, attrs in branch_attrs.items()
             for attrs_label in [lambda x, b=b: TIP_INFO.get(b, {}).get("label", b)]}
    p_corr = {attrs_label(b): attrs["Corrected P-value"] for b, attrs in branch_attrs.items()
              for attrs_label in [lambda x, b=b: TIP_INFO.get(b, {}).get("label", b)]}

    n_significant = absrel["test results"]["positive test results"]
    n_tested = absrel["test results"]["tested"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})

    ax = axes[0]
    cmap = plt.get_cmap("RdYlBu_r")
    norm = Normalize(vmin=0, vmax=1.3)

    def bc(om):
        return cmap(norm(om))

    y_pos = {"Lbc1": 0, "Lba": 1, "Lbc3": 2.3, "Lbc2": 3.3}
    x_root, x_node4, x_tip = 0.0, 0.5, 1.0

    ax.plot([x_root, x_tip], [y_pos["Lbc1"]]*2, lw=4, color=bc(omega["Lbc1"]), solid_capstyle="butt")
    ax.plot([x_root, x_tip], [y_pos["Lba"]]*2, lw=4, color=bc(omega["Lba"]), solid_capstyle="butt")

    y_node4 = (y_pos["Lbc3"] + y_pos["Lbc2"]) / 2
    ax.plot([x_root, x_node4], [y_node4]*2, lw=4, color=bc(omega["Node4"]), solid_capstyle="butt")
    ax.plot([x_root, x_root], [min(y_pos["Lbc1"], y_pos["Lba"], y_node4),
                                max(y_pos["Lbc1"], y_pos["Lba"], y_node4)], lw=1.2, color="grey")

    ax.plot([x_node4, x_tip], [y_pos["Lbc3"]]*2, lw=4, color=bc(omega["Lbc3"]), solid_capstyle="butt")
    ax.plot([x_node4, x_tip], [y_pos["Lbc2"]]*2, lw=4, color=bc(omega["Lbc2"]), solid_capstyle="butt")
    ax.plot([x_node4, x_node4], [y_pos["Lbc3"], y_pos["Lbc2"]], lw=1.2, color="grey")

    for b, v in TIP_INFO.items():
        y = y_pos[v["label"]]
        status = "expressed" if v["detected"] else "undetected"
        ax.text(x_tip + 0.03, y, f"{v['label']}  ({status})", fontsize=7.2, va="center", ha="left")

    ax.text(x_node4 - 0.02, y_node4 + 0.5, "Node4", fontsize=6.3, ha="center", va="bottom", color="dimgrey")
    ax.text(0.85, -0.85,
            "HyPhy unroots the tree to a trifurcating root\n(Lbc1, Lba each attach directly to the root;\n"
            "no separate Lbc1-Lba branch is tested)",
            fontsize=5.8, ha="center", va="top", color="dimgrey")

    ax.set_xlim(-0.4, 2.1)
    ax.set_ylim(-1.3, 4.0)
    ax.axis("off")
    ax.set_title(f"Branch color = fitted $\\omega$ (dN/dS);\n"
                 f"no branch reaches significance (Holm-Bonferroni p=1, all {n_tested})", fontsize=8)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.1, aspect=22)
    cbar.set_label("$\\omega$ (dN/dS)", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.axvline(1.0, color="black", lw=0.8, ls="--")
    cbar.ax.set_xlim(0, 1.28)

    ax2 = axes[1]
    order = ["Lbc2", "Lbc3", "Lbc1", "Lba"]
    omegas = [omega[l] for l in order]
    detected_flags = [TIP_INFO[[b for b, v in TIP_INFO.items() if v["label"] == l][0]]["detected"] for l in order]
    colors_bar = ["#1b6ca8" if d else "#a8a8a8" for d in detected_flags]

    ax2.barh(order, omegas, color=colors_bar, height=0.6)
    ax2.axvline(1.0, color="black", lw=0.8, ls="--")
    ax2.text(1.02, 3.35, "neutral ($\\omega$=1)", fontsize=6, color="dimgrey")
    for i, val in enumerate(omegas):
        ax2.text(val + 0.03, i, f"{val:.2f}", va="center", fontsize=7)
    ax2.set_xlim(0, 1.35)
    ax2.set_xlabel("Terminal-branch $\\omega$ (dN/dS)")
    ax2.set_title("Undetected paralogs trend toward\nrelaxed constraint (not individually significant)", fontsize=8)
    for spine in ("top", "right"):
        ax2.spines[spine].set_visible(False)

    legend_handles = [Patch(color="#1b6ca8", label="Expressed (nodule, GSE270392)"),
                      Patch(color="#a8a8a8", label="Undetected in GSE270392")]
    ax2.legend(handles=legend_handles, frameon=False, fontsize=6.5, loc="lower right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")
    print(json.dumps({l: {"omega": round(omega[l], 4), "corrected_p": p_corr[l]} for l in omega}, indent=2))
    print(f"Positive test results: {n_significant} / {n_tested} branches")


if __name__ == "__main__":
    main()
