"""
Prune the full 197-taxon OrthoDB group 706508at2759 ML tree
(../../orthodb_706508at2759_leghemoglobin/results/tree_706508at2759.treefile)
down to just the 4 soybean (Glycine max) leghemoglobin paralogs that have a
valid coding sequence (fetch_soybean_cds.py), for the aBSREL branch-site
selection analysis. Restricting to this small "soybean portion" of the tree
keeps the codon alignment (PRANK) and branch-site test (HyPhy aBSREL) fast.

Usage:
    python prune_tree.py <full_tree.treefile> <output_pruned.nwk>
"""
import copy
import sys

from Bio import Phylo

TARGET_ODB_IDS = ["3847_0:0065e0", "3847_0:006217", "3847_0:0065fa", "3847_0:00ca2d"]


def tree_label(odb_id):
    return f"Glycine_max__{odb_id.replace(':', '_')}"


def main():
    in_path, out_path = sys.argv[1], sys.argv[2]
    full_tree = Phylo.read(in_path, "newick")
    all_tips = [t.name for t in full_tree.get_terminals()]
    print(f"Full tree: {len(all_tips)} tips")

    target_labels = [tree_label(o) for o in TARGET_ODB_IDS]
    missing = [t for t in target_labels if t not in all_tips]
    if missing:
        raise ValueError(f"Target tips not found in tree: {missing}")

    pruned = copy.deepcopy(full_tree)
    for tip_name in all_tips:
        if tip_name not in target_labels:
            clade = next(pruned.find_clades(name=tip_name), None)
            if clade is not None:
                pruned.prune(clade)

    remaining = [t.name for t in pruned.get_terminals()]
    print(f"Pruned tree: {len(remaining)} tips -> {remaining}")

    Phylo.write(pruned, out_path, "newick")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
