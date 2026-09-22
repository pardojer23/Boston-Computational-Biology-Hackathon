"""
Apply the paralog-expression-classification framework from Benoit et al.,
"Solanum pan-genetics reveals paralogues as contingencies in crop engineering"
(Nature 2025, DOI 10.1038/s41586-025-08619-6) to the soybean leghemoglobin
paralogs, using the GSE270392 tissue- and nodule-cell-type-resolved pseudobulk
expression data already assembled by redundancy_analysis.py.

Four metrics are computed exactly as defined in that paper's Methods:

1. Tissue-specificity index tau (Yanai et al. 2005, Bioinformatics 21:650-659):
   tau = sum(1 - x_i/max(x)) / (n_tissues - 1)

2. Expression breadth: number of tissues (of our 7-tissue panel: nodule, root,
   hypocotyl, and 4 seed developmental stages) with average expression > 3 CPM
   (CPM stands in for the paper's TPM -- both are library-size-normalized
   relative-abundance units; we lack per-gene length information for a true
   TPM, but CPM is the appropriate substitute at the pseudobulk level used
   here).

3. Non-functional / tissue-specific gene calls, using the paper's exact
   thresholds: non-functional if average expression < 3 CPM across the
   tissue panel; tissue-specific to tissue X if expression is highest in X,
   tau > 0.7, and expression in X exceeds 5 CPM.

4. Paralogue-pair expression-group classification (I: dosage balanced, II:
   paralogue dominance, III: specialized, IV: diverged), based on a
   rank-standardized coexpression network value and the mean/s.d. of
   |log2(fold change)| across samples, using the paper's exact thresholds:
     I   coexpression > 0.9, mean|log2FC| < 1, s.d.|log2FC| < 1
     II  coexpression > 0.9, mean|log2FC| >= 1, s.d.|log2FC| < 1
     III coexpression > 0.9, mean|log2FC| >= 1, s.d.|log2FC| >= 1
     IV  coexpression < 0.5, mean|log2FC| >= 1, s.d.|log2FC| >= 1
   "Samples" here are the 14 annotated nodule cell states (the finest-grained
   resolution at which both paralogs of the pair have expression data), the
   within-species analogue of the paper's tissue-replicate samples.

Usage:
    python benoit_metrics.py
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

DETECTED_PARALOGS = {"Glyma.10G198800": "Lbc3 (LB1)", "Glyma.20G191200": "Lbc2 (LB4)"}
UNDETECTED_PARALOGS = {
    "Glyma.10G199000": "Lbc1 (LB2) -- N-50 / Nodulin-50",
    "Glyma.10G199100": "Lba (LB3) -- N-2 / Nodulin-2",
    "Glyma.10G198900": "unnamed paralog near cluster (NCBI: pseudogene LB5)",
}
ALL_PARALOGS = {**DETECTED_PARALOGS, **UNDETECTED_PARALOGS}


def load_tissue_matrix(filename):
    with open(filename, "rb") as f:
        text = gzip.decompress(f.read()).decode()
    return pd.read_csv(io.StringIO(text), sep="\t", index_col=0)


def tau(expr_vec):
    """Yanai et al. 2005 tissue-specificity index."""
    x = np.asarray(expr_vec, dtype=float)
    if x.max() <= 0:
        return np.nan
    xhat = x / x.max()
    n = len(x)
    return float(np.sum(1 - xhat) / (n - 1))


def tissue_level_cpm(dfs):
    rows = {}
    for tissue, df in dfs.items():
        lib_total = df.values.sum()
        row = {label: (df.loc[gid].sum() / lib_total * 1e6 if gid in df.index else 0.0)
               for gid, label in ALL_PARALOGS.items()}
        rows[tissue] = row
    return pd.DataFrame(rows).T


def rank_standardize(corr_series):
    s = corr_series.copy()
    n = len(s)
    ranks = s.rank(method="average", na_option="keep")
    ranks[s.isna()] = (n + 1) / 2.0
    return ranks / ranks.max()


def corr_vector(target_gene, log_mat):
    X = log_mat.values
    idx = log_mat.index.get_loc(target_gene)
    x = X[idx]
    x_c = x - x.mean()
    X_c = X - X.mean(axis=1, keepdims=True)
    num = X_c @ x_c
    denom = np.sqrt((X_c ** 2).sum(axis=1) * (x_c ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / denom
    return pd.Series(r, index=log_mat.index)


def classify_group(coexpr, mean_lfc, sd_lfc):
    if coexpr > 0.9 and mean_lfc < 1 and sd_lfc < 1:
        return "I: Dosage balanced"
    if coexpr > 0.9 and mean_lfc >= 1 and sd_lfc < 1:
        return "II: Paralogue dominance"
    if coexpr > 0.9 and mean_lfc >= 1 and sd_lfc >= 1:
        return "III: Specialized"
    if coexpr < 0.5 and mean_lfc >= 1 and sd_lfc >= 1:
        return "IV: Diverged"
    return "unclassified (paper's thresholds are a partial partition of coexpr x FC space)"


def main():
    dfs = {tissue: load_tissue_matrix(fn) for tissue, fn in TISSUE_FILES.items()}
    tissue_df = tissue_level_cpm(dfs)

    # --- tau, expression breadth, non-functional / tissue-specific calls ---
    per_gene = {}
    for gid, label in ALL_PARALOGS.items():
        expr = tissue_df[label].values
        t = tau(expr)
        breadth = int((tissue_df[label] > 3).sum())
        avg_expr = tissue_df[label].mean()
        nonfunctional = bool(avg_expr < 3)
        highest_in_nodule = tissue_df[label].idxmax() == "Early_nodule"
        nodule_expr = tissue_df.loc["Early_nodule", label]
        tissue_specific = bool(highest_in_nodule and (not np.isnan(t)) and t > 0.7 and nodule_expr > 5)
        per_gene[label] = {
            "glyma_id": gid,
            "tau": None if np.isnan(t) else round(t, 4),
            "expression_breadth_tissues_CPM_gt_3": breadth,
            "avg_CPM_across_7_tissues": round(float(avg_expr), 4),
            "nonfunctional_by_paper_rule": nonfunctional,
            "nodule_specific_by_paper_rule": tissue_specific,
        }
    print(json.dumps(per_gene, indent=2))

    # --- pairwise coexpression / fold-change classification (Lbc3 vs Lbc2) ---
    nod = dfs["Early_nodule"]
    lib_size = nod.sum(axis=0)
    cpm = nod.div(lib_size, axis=1) * 1e6

    col_medians = cpm.median(axis=0)
    keep = cpm.gt(col_medians, axis=1).any(axis=1)
    filtered = cpm.loc[keep]
    log_mat = np.log2(filtered + 1)

    g1, g2 = list(DETECTED_PARALOGS.keys())
    corr1 = corr_vector(g1, log_mat)
    corr2 = corr_vector(g2, log_mat)
    r1 = rank_standardize(corr1)
    r2 = rank_standardize(corr2)
    coexpr_value = float((r1[g2] + r2[g1]) / 2.0)

    paired = cpm.loc[[g1, g2]].T
    paired.columns = ["A", "B"]
    nz = (paired["A"] > 0) & (paired["B"] > 0)
    log2fc = np.log2(paired.loc[nz, "A"] / paired.loc[nz, "B"])
    mean_lfc = float(log2fc.abs().mean())
    sd_lfc = float(log2fc.abs().std())
    group = classify_group(coexpr_value, mean_lfc, sd_lfc)

    pair_result = {
        "pair": [DETECTED_PARALOGS[g1], DETECTED_PARALOGS[g2]],
        "n_samples_used_nonzero_both": int(nz.sum()),
        "n_samples_total": int(len(paired)),
        "genes_in_coexpression_network": int(filtered.shape[0]),
        "rank_standardized_coexpression": round(coexpr_value, 4),
        "mean_abs_log2FC": round(mean_lfc, 4),
        "sd_abs_log2FC": round(sd_lfc, 4),
        "benoit_group_classification": group,
    }
    print(json.dumps(pair_result, indent=2))

    with open("benoit_metrics_results.json", "w") as f:
        json.dump({"per_gene": per_gene, "pair_classification": pair_result}, f, indent=2)


if __name__ == "__main__":
    main()
