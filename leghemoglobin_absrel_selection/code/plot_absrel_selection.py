"""
Interpret HyPhy aBSREL branch-site test results (modal_absrel.py output) for
the 4-taxon soybean leghemoglobin quartet in terms of paralog redundancy:
plot the tree with branches colored by fitted omega (dN/dS), and a bar chart
comparing terminal-branch omega against each paralog's expression-detection
status in GSE270392 (see ../../gse270392_leghemoglobin_redundancy/).

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

    branch_attrs = absrel["branch attributes"]["0"]
    omega = {}
    p_corr = {}
    for branch, attrs in branch_attrs.items():
        omega[branch] = attrs["Baseline MG94xREV omega ratio"]
        p_corr[branch] = attrs["Corrected P-value"]

    internal_branch = [b for b in branch_attrs if b not in TIP_INFO][0]
    internal_omega = omega[internal_branch]

    n_significant = absrel["test results"]["positive test results"]
    n_tested = absrel["test results"]["tested"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.3, 1]})

    ax = axes[0]
    cmap = plt.get_cmap("RdYlBu_r")
    norm = Normalize(vmin=0, vmax=1.3)

    def bc(om):
        return cmap(norm(om))

    y_pos = {"Lbc1": 0, "Lba": 1, "Lbc3": 2, "Lbc2": 3}
    label_of = {b: v["label"] for b, v in TIP_INFO.items()}
    x_tip, x_int1, x_int2, x_root = 1.0, 0.5, 0.5, 0.0

    om_by_label = {v["label"]: omega[b] for b, v in TIP_INFO.items()}

    ax.plot([x_int1, x_tip], [y_pos["Lbc1"]]*2, lw=4, color=bc(om_by_label["Lbc1"]), solid_capstyle="butt")
    ax.plot([x_int1, x_tip], [y_pos["Lba"]]*2, lw=4, color=bc(om_by_label["Lba"]), solid_capstyle="butt")
    ax.plot([x_int1, x_int1], [y_pos["Lbc1"], y_pos["Lba"]], lw=1.2, color="grey")

    ax.plot([x_int2, x_tip], [y_pos["Lbc3"]]*2, lw=4, color=bc(om_by_label["Lbc3"]), solid_capstyle="butt")
    ax.plot([x_int2, x_tip], [y_pos["Lbc2"]]*2, lw=4, color=bc(om_by_label["Lbc2"]), solid_capstyle="butt")
    ax.plot([x_int2, x_int2], [y_pos["Lbc3"], y_pos["Lbc2"]], lw=1.2, color="grey")

    y_int1 = (y_pos["Lbc1"] + y_pos["Lba"]) / 2
    y_int2 = (y_pos["Lbc3"] + y_pos["Lbc2"]) / 2
    ax.plot([x_root, x_int1], [y_int1]*2, lw=4, color=bc(internal_omega), solid_capstyle="butt")
    ax.plot([x_root, x_int2], [y_int2]*2, lw=1.2, color="grey")
    ax.plot([x_root, x_root], [y_int1, y_int2], lw=1.2, color="grey")

    for b, v in TIP_INFO.items():
        y = y_pos[v["label"]]
        status = "expressed" if v["detected"] else "undetected"
        ax.text(x_tip + 0.03, y, f"{v['label']}  ({status})", fontsize=7.2, va="center", ha="left")

    ax.text(x_root - 0.32, (y_int1+y_int2)/2, "internal\nbranch", fontsize=6, ha="center", va="center", color="dimgrey")
    ax.text(0.75, -0.75, "grey = topology only (not separately tested);\ncolored = one of the 5 branches aBSREL tested",
            fontsize=5.8, ha="center", va="top", color="dimgrey")

    ax.set_xlim(-0.5, 2.0)
    ax.set_ylim(-1.1, 3.8)
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
    omegas = [om_by_label[l] for l in order]
    detected_flags = [TIP_INFO[[b for b, v in TIP_INFO.items() if v["label"] == l][0]]["detected"] for l in order]
    colors_bar = ["#1b6ca8" if d else "#a8a8a8" for d in detected_flags]

    ax2.barh(order, omegas, color=colors_bar, height=0.6)
    ax2.axvline(1.0, color="black", lw=0.8, ls="--")
    ax2.text(1.02, 3.35, "neutral ($\\omega$=1)", fontsize=6, color="dimgrey")
    for i, om in enumerate(omegas):
        ax2.text(om + 0.03, i, f"{om:.2f}", va="center", fontsize=7)
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
    print(json.dumps({label_of.get(b, b): {"omega": round(omega[b], 4), "corrected_p": p_corr[b]}
                       for b in branch_attrs}, indent=2))
    print(f"Positive test results: {n_significant} / {n_tested} branches")


if __name__ == "__main__":
    main()
