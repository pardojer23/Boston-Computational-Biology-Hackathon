"""
Render the IQ-TREE consensus/ML tree (Newick) for OrthoDB group
706508at2759 as a labeled phylogram, using organism names (and a gene-id
suffix for in-paralogs) as tip labels.

Usage:
    python render_tree_figure.py tree_706508at2759.treefile \\
        group_706508at2759_clean.fasta tree_706508at2759_full.png
"""
import io
import re
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo


def load_label_map(clean_fasta_path: str) -> dict:
    """Map the Newick-safe tip id back to organism / gene-id metadata."""
    label_map = {}
    with open(clean_fasta_path) as f:
        text = f.read()
    for block in text.strip().split(">")[1:]:
        lines = block.split("\n")
        label = lines[0].strip()
        # label format: Organism_name__odb_id (see fetch_orthodb_group.py)
        m = re.match(r"(.+)__([^_].*)$", label)
        if not m:
            continue
        org = m.group(1).replace("_", " ")
        label_map[label] = {"organism": org, "odb_id": m.group(2)}
    return label_map


def assign_display_names(tree, label_map: dict) -> None:
    counts = Counter()
    for tip in tree.get_terminals():
        r = label_map.get(tip.name)
        org = r["organism"] if r else tip.name
        counts[org] += 1
    for tip in tree.get_terminals():
        r = label_map.get(tip.name)
        if r is None:
            tip.display = tip.name
            continue
        org = r["organism"]
        tip.display = f"{org} ({tip.name.split('__')[-1]})" if counts[org] > 1 else org
    # de-duplicate any leftover collisions
    seen = Counter()
    for tip in tree.get_terminals():
        seen[tip.display] += 1
    for tip in tree.get_terminals():
        if seen[tip.display] > 1:
            tip.display = f"{tip.display} [{tip.name}]"


def render(tree, out_path: str) -> None:
    tree.ladderize()
    terminals = tree.get_terminals()
    n = len(terminals)
    depths = tree.depths(unit_branch_lengths=False)

    y_pos = {t: i for i, t in enumerate(terminals)}

    def assign_y(clade):
        if clade.is_terminal():
            return y_pos[clade]
        ys = [assign_y(c) for c in clade.clades]
        y_pos[clade] = sum(ys) / len(ys)
        return y_pos[clade]

    assign_y(tree.root)

    fig_h = max(24, n * 0.16)
    fig, ax = plt.subplots(figsize=(11, fig_h))

    def draw_clade(clade, x_start):
        x_end = x_start + (clade.branch_length or 0.0)
        y = y_pos[clade]
        ax.plot([x_start, x_end], [y, y], color="#333333", lw=0.8, solid_capstyle="butt")
        if clade.clades:
            ys = [y_pos[c] for c in clade.clades]
            ax.plot([x_end, x_end], [min(ys), max(ys)], color="#333333", lw=0.8)
            for c in clade.clades:
                draw_clade(c, x_end)
        return x_end

    draw_clade(tree.root, 0.0)

    for t in terminals:
        ax.text(depths[t] + 0.01, y_pos[t], t.display, fontsize=4.3, va="center",
                 ha="left", style="italic", color="#111111")

    xmax = max(depths.values()) * 1.42
    ax.set_xlim(0, xmax)
    ax.set_xticks([tk for tk in ax.get_xticks() if tk <= xmax])
    ax.set_ylim(-1, n)
    ax.invert_yaxis()
    ax.set_xlabel("Substitutions per site (JTT+I+G4)")
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_title("Maximum-likelihood phylogeny of leghemoglobin A orthologs "
                 "(OrthoDB group 706508at2759)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)


def main() -> None:
    treefile_path, clean_fasta_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(treefile_path) as f:
        tree = Phylo.read(io.StringIO(f.read()), "newick")
    label_map = load_label_map(clean_fasta_path)
    assign_display_names(tree, label_map)
    render(tree, out_path)
    print(f"Wrote {out_path} ({tree.count_terminals()} tips)")


if __name__ == "__main__":
    main()
