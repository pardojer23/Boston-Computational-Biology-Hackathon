"""
Build the report figures for one OrthoDB-group pipeline run: the ML tree
(colored by expression-detection status where soybean expression data is
available) and, if a selection analysis ran, an aBSREL branch-omega figure.

The aBSREL figure reads branch topology directly from HyPhy's own recorded
input tree (`absrel_result["input"]["trees"]`) rather than assuming it
matches the guide tree's rooted topology -- HyPhy unroots the tree before
testing, which can change which branches are adjacent (this bit us once on
the leghemoglobin case: see leghemoglobin_absrel_selection/README.md).

Usage (as a library -- see pipeline.py):
    from report import plot_absrel_branches, write_markdown_report
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch


def _parse_newick_topology(newick_str: str):
    """Minimal Newick parser returning a nested-list tree: each node is
    (name_or_None, branch_length_or_None, [children]). Sufficient for the
    flat, small trees aBSREL reports (no quoted labels, no NHX)."""
    s = newick_str.strip().rstrip(";")
    pos = 0

    def parse_node():
        nonlocal pos
        children = []
        if s[pos] == "(":
            pos += 1
            while True:
                children.append(parse_node())
                if s[pos] == ",":
                    pos += 1
                    continue
                elif s[pos] == ")":
                    pos += 1
                    break
        # optional name[:branch_length] -- ':' must stay IN the label text
        # here (not a stop character) so it can be split out below.
        start = pos
        while pos < len(s) and s[pos] not in ",();":
            pos += 1
        label_part = s[start:pos]
        name, bl = None, None
        if ":" in label_part:
            name_part, bl_part = label_part.split(":", 1)
            name = name_part or None
            bl = float(bl_part) if bl_part else None
        elif label_part:
            name = label_part
        return (name, bl, children)

    return parse_node()


def plot_absrel_branches(absrel_result_json: dict, detected_labels: set, out_path: str,
                          detected_legend="Expressed", undetected_legend="Undetected"):
    """Generic branch-omega tree figure for an arbitrary aBSREL result.
    `detected_labels` is the set of tip labels considered "expressed" for
    the accompanying legend/coloring (pass an empty set to skip that
    distinction). Layout: simple ladder dendrogram using cumulative branch
    length for x and in-order tip position for y -- works for any topology
    aBSREL actually reports, read from its own input tree, not assumed."""
    branch_attrs = absrel_result_json["branch attributes"]["0"]
    omega = {b: a["Baseline MG94xREV omega ratio"] for b, a in branch_attrs.items()}
    p_corr = {b: a["Corrected P-value"] for b, a in branch_attrs.items()}
    tested_tree_str = absrel_result_json["input"]["trees"]["0"]
    root = _parse_newick_topology(tested_tree_str)

    # assign y (tip order via DFS) and x (cumulative branch length)
    y_counter = [0]
    y_pos = {}
    x_pos = {}

    def assign(node, x_parent):
        name, bl, children = node
        x = x_parent + (bl or 0.0)
        if not children:
            y_pos[id(node)] = y_counter[0]
            y_counter[0] += 1
        else:
            for c in children:
                assign(c, x)
            ys = [y_pos[id(c)] for c in children]
            y_pos[id(node)] = sum(ys) / len(ys)
        x_pos[id(node)] = x
        return

    assign(root, 0.0)

    fig, ax = plt.subplots(figsize=(7.5, 0.5 * max(6, y_counter[0] * 1.4)))
    cmap = plt.get_cmap("RdYlBu_r")
    all_om = [v for v in omega.values() if v is not None]
    norm = Normalize(vmin=0, vmax=max(1.3, max(all_om) * 1.05 if all_om else 1.3))

    def draw(node, x_parent):
        name, bl, children = node
        x = x_pos[id(node)]
        y = y_pos[id(node)]
        branch_name = name if name else None
        color = cmap(norm(omega[branch_name])) if branch_name in omega else "grey"
        lw = 4 if branch_name in omega else 1.2
        ax.plot([x_parent, x], [y, y], lw=lw, color=color, solid_capstyle="butt")
        if children:
            ys = [y_pos[id(c)] for c in children]
            ax.plot([x, x], [min(ys), max(ys)], lw=1.2, color="grey")
            for c in children:
                draw(c, x)
        else:
            status = ""
            if detected_labels:
                status = "  (expressed)" if name in detected_labels else "  (undetected)"
            ax.text(x + 0.01 * max(x_pos.values()), y, f"{name}{status}", fontsize=7, va="center", ha="left")

    draw(root, 0.0)

    n_sig = absrel_result_json["test results"]["positive test results"]
    n_tested = absrel_result_json["test results"]["tested"]
    ax.set_title(f"Branch color = fitted $\\omega$ (dN/dS); "
                 f"{n_sig}/{n_tested} branches significant (Holm-Bonferroni)", fontsize=8)
    ax.set_xlim(-0.02 * max(x_pos.values()), max(x_pos.values()) * 1.5)
    ax.axis("off")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.08, aspect=25)
    cbar.set_label("$\\omega$ (dN/dS)", fontsize=6.5)
    cbar.ax.axvline(1.0, color="black", lw=0.8, ls="--")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    return {"omega": omega, "corrected_p": p_corr, "n_significant": n_sig, "n_tested": n_tested}


def write_markdown_report(group_summary: dict, expression_result: dict, selection_result: dict,
                           out_path: str, figure_paths: dict):
    lines = [f"# OrthoDB group {group_summary['group_id']} -- {group_summary.get('group_name')}", ""]
    lines.append(f"- {group_summary['n_genes']} genes across {group_summary['n_species']} species")
    lines.append("")

    lines.append("## Soybean expression analysis (GSE270392)")
    if expression_result.get("skipped"):
        lines.append(f"Skipped: {expression_result['reason']}")
    else:
        lines.append(f"- {expression_result['n_gmax_records']} Glycine max record(s); "
                      f"{len(expression_result['detected_glyma_ids'])} detected in GSE270392")
        if "tree_figure" in figure_paths:
            lines.append(f"\n![tree]({figure_paths['tree_figure']})\n")
        for pc in expression_result.get("pair_classifications", []):
            if pc.get("classifiable"):
                lines.append(f"- {pc['gene_a']} vs {pc['gene_b']}: **{pc['benoit_group']}** "
                              f"(coexpr={pc['rank_standardized_coexpression']}, "
                              f"mean|log2FC|={pc['mean_abs_log2FC']})")
    lines.append("")

    lines.append("## Branch-site selection (HyPhy aBSREL)")
    if selection_result.get("skipped"):
        lines.append(f"Skipped: {selection_result['reason']}")
    else:
        if "absrel_figure" in figure_paths:
            lines.append(f"\n![absrel]({figure_paths['absrel_figure']})\n")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))
    return out_path
