"""
Test cross-tissue / cross-cell-type expression redundancy among the soybean
leghemoglobin paralogs, using the Glyma-ID-mapped pseudobulk matrices
downloaded by fetch_gse270392.py and the ID map produced by
map_leghemoglobin_ids.py.

Two of the five OrthoDB-derived paralog loci are detectable in the released
matrices (Lbc3 = Glyma.10G198800, Lbc2 = Glyma.20G191200); the other two
named paralogs (Lbc1 = Glyma.10G199000, Lba = Glyma.10G199100) and one
unnamed locus (Glyma.10G198900) are absent from every tissue's matrix, and
modal_seurat_probe.py additionally confirms they are absent even from the
unfiltered "RNA" assay of the full single-nucleus Seurat object -- i.e. this
is a genuine non-detection in the released single-nucleus dataset, not a
gene-ID mismatch or a filtering artifact of the pre-aggregated summary file.

Usage:
    python redundancy_analysis.py
"""
import gzip
import io
import json

import numpy as np
import pandas as pd

TISSUE_FILES = {
    "Early_nodule": "GSE270392_Gm_atlas_Early_nodule_RNA_gene_cellstate_counts.txt.gz",
    "Root": "GSE270392_Gm_atlas_Root_RNA_gene_cellstate_counts.txt.gz",
    "Hypocotyl": "GSE270392_Gm_atlas_Hypocotyl_RNA_gene_cellstate_counts.txt.gz",
    "Cotyledon_stage_seeds": "GSE270392_Gm_atlas_Cotyledon_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Early_maturation_stage_seeds": "GSE270392_Gm_atlas_Early_maturation_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Globular_stage_seeds": "GSE270392_Gm_atlas_Globular_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Heart_stage_seeds": "GSE270392_Gm_atlas_Heart_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
}

DETECTED_PARALOGS = {
    "Glyma.10G198800": "Lbc3 (LB1)",
    "Glyma.20G191200": "Lbc2 (LB4)",
}
UNDETECTED_PARALOGS = {
    "Glyma.10G199000": "Lbc1 (LB2) -- N-50 / Nodulin-50",
    "Glyma.10G199100": "Lba (LB3) -- N-2 / Nodulin-2",
    "Glyma.10G198900": "unnamed paralog near cluster (NCBI: pseudogene LB5)",
}


def load_tissue_matrix(filename: str) -> pd.DataFrame:
    with open(filename, "rb") as f:
        text = gzip.decompress(f.read()).decode()
    return pd.read_csv(io.StringIO(text), sep="\t", index_col=0)


def tissue_level_cpm(dfs: dict) -> pd.DataFrame:
    """Total CPM per tissue for every target paralog (sum across all cell
    states in that tissue), including zero rows for undetected genes."""
    rows = {}
    for tissue, df in dfs.items():
        lib_total = df.values.sum()
        row = {}
        for gid, label in {**DETECTED_PARALOGS, **UNDETECTED_PARALOGS}.items():
            row[label] = (df.loc[gid].sum() / lib_total * 1e6) if gid in df.index else 0.0
        rows[tissue] = row
    return pd.DataFrame(rows).T


def nodule_celltype_cpm(nodule_df: pd.DataFrame) -> pd.DataFrame:
    """Per-cell-state CPM within the nodule for the two detected paralogs."""
    lib_size = nodule_df.sum(axis=0)
    cpm = nodule_df.loc[list(DETECTED_PARALOGS.keys())].div(lib_size, axis=1) * 1e6
    cpm.index = [DETECTED_PARALOGS[g] for g in cpm.index]
    return cpm.T.sort_values("Lbc3 (LB1)", ascending=False)


def main() -> None:
    dfs = {tissue: load_tissue_matrix(fn) for tissue, fn in TISSUE_FILES.items()}

    tissue_df = tissue_level_cpm(dfs)
    print("Tissue-level CPM (rows = tissue, cols = paralog):")
    print(tissue_df.round(2).to_string())
    tissue_df.to_csv("tissue_level_cpm.csv")

    ct_cpm = nodule_celltype_cpm(dfs["Early_nodule"])
    print("\nNodule cell-state CPM (detected paralogs only):")
    print(ct_cpm.round(2).to_string())
    ct_cpm.to_csv("nodule_celltype_cpm.csv")

    r_raw = np.corrcoef(dfs["Early_nodule"].loc["Glyma.10G198800"],
                         dfs["Early_nodule"].loc["Glyma.20G191200"])[0, 1]
    r_cpm = np.corrcoef(ct_cpm["Lbc3 (LB1)"], ct_cpm["Lbc2 (LB4)"])[0, 1]
    print(f"\nPearson r (Lbc3 vs Lbc2 across nodule cell states), raw counts: {r_raw:.3f}")
    print(f"Pearson r (Lbc3 vs Lbc2 across nodule cell states), CPM:         {r_cpm:.3f}")

    with open("redundancy_summary.json", "w") as f:
        json.dump({
            "pearson_r_raw_counts": round(float(r_raw), 4),
            "pearson_r_cpm": round(float(r_cpm), 4),
            "detected_paralogs": DETECTED_PARALOGS,
            "undetected_paralogs": UNDETECTED_PARALOGS,
        }, f, indent=2)


if __name__ == "__main__":
    main()
