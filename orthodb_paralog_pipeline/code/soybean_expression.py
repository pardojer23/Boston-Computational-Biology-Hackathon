"""
Map Glycine max (soybean) records from an OrthoDB orthologous group onto
Glyma gene-locus IDs, and compute expression-based redundancy metrics from
the GSE270392 soybean multiomic single-cell atlas (Zhang et al., Cell 2024).

Generalized from the leghemoglobin-specific analysis in
gse270392_leghemoglobin_redundancy/ -- works for any OrthoDB group's
Glycine max records, with graceful skips when a group has none, or when
GSE270392 has no expression data for the resolved genes.

Assumes soybean / GSE270392 throughout (per this pipeline's scope): this is
not a general cross-species expression-mapping module.

Usage (as a library -- see pipeline.py for the orchestrated call):
    from soybean_expression import run_soybean_expression_analysis
"""
import gzip
import io
import json
import os
import re
import time
import urllib.parse
import urllib.request
import uuid
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270392/suppl/"
SOYBASE_BLAST = "https://sequenceserver.soybase.org/"
SOYBASE_WM82A4_PROTEIN_DB_ID = "982cbf1cfaa76dd4e7589d7d1b2d60a0"

TISSUE_FILES = {
    "Early_nodule": "GSE270392_Gm_atlas_Early_nodule_RNA_gene_cellstate_counts.txt.gz",
    "Root": "GSE270392_Gm_atlas_Root_RNA_gene_cellstate_counts.txt.gz",
    "Hypocotyl": "GSE270392_Gm_atlas_Hypocotyl_RNA_gene_cellstate_counts.txt.gz",
    "Cotyledon_stage_seeds": "GSE270392_Gm_atlas_Cotyledon_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Early_maturation_stage_seeds": "GSE270392_Gm_atlas_Early_maturation_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Globular_stage_seeds": "GSE270392_Gm_atlas_Globular_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Heart_stage_seeds": "GSE270392_Gm_atlas_Heart_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
}


def _fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "research-agent"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read()


def fetch_gse270392_tissue_matrices(cache_dir: str) -> dict:
    """Download (once) and cache the 7 per-tissue gene x cell-state pseudobulk
    matrices. These are the same regardless of which OrthoDB group is being
    analyzed, so they are cached under cache_dir and reused across runs."""
    os.makedirs(cache_dir, exist_ok=True)
    dfs = {}
    for tissue, filename in TISSUE_FILES.items():
        local_path = os.path.join(cache_dir, filename)
        if not os.path.exists(local_path):
            raw = _fetch(GEO_SUPPL + filename)
            with open(local_path, "wb") as f:
                f.write(raw)
        with open(local_path, "rb") as f:
            text = gzip.decompress(f.read()).decode()
        dfs[tissue] = pd.read_csv(io.StringIO(text), sep="\t", index_col=0)
    return dfs


def is_glycine_max_record(rec: dict) -> bool:
    return (rec.get("organism") or "").strip() == "Glycine max"


def _blast_query_soybase(records: list) -> list:
    """BLASTP the given OrthoDB records' protein sequences against SoyBase's
    Wm82.a4.v1 protein database; return the top-hit Glyma ID + confidence
    per record. Uses a manual multipart POST (the `requests` library fails
    through this sandbox's proxy for this domain; urllib works)."""
    query_fasta = "\n".join(f">{r['odb_id']}\n{r['seq']}" for r in records)

    boundary = uuid.uuid4().hex

    def field(name, value):
        return (f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n').encode()

    body = b"".join([
        field("sequence", query_fasta),
        field("method", "blastp"),
        field("databases[]", SOYBASE_WM82A4_PROTEIN_DB_ID),
    ]) + f"--{boundary}--\r\n".encode()

    req = urllib.request.Request(
        SOYBASE_BLAST, data=body, method="POST",
        headers={"User-Agent": "research-agent", "Content-Type": f"multipart/form-data; boundary={boundary}"},
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        job_url = r.url

    result = None
    for _ in range(20):
        raw = _fetch(job_url + ".json")
        try:
            result = json.loads(raw)
        except json.JSONDecodeError:
            time.sleep(3)
            continue
        if result.get("queries") is not None:
            break
        time.sleep(3)
    if result is None or result.get("queries") is None:
        raise TimeoutError(f"BLAST job at {job_url} did not finish in time")

    hits = []
    for q in result["queries"]:
        top = q["hits"][0] if q["hits"] else None
        hits.append({
            "odb_id": q["id"],
            "glyma_id": top["accession"].rsplit(".", 1)[0] if top else None,  # strip transcript-version suffix
            "bit_score": top["hsps"][0]["bit_score"] if top else None,
            "identity": top["hsps"][0].get("identity") if top else None,
            "query_length": q["length"],
        })
    return hits


def resolve_glyma_ids(gmax_records: list) -> list:
    """Resolve each Glycine max OrthoDB record to a Glyma.##G###### locus ID.

    Two-tier strategy:
      1. Direct reformat when OrthoDB's own pub_gene_id already carries a
         "GLYMA_##G######vN" tag.
      2. BLASTP the protein sequence against SoyBase's public Wm82.a4.v1
         protein database for every record without a usable direct tag
         (placeholder gene_ids like historical "N-50"-style symbols, or a
         missing tag entirely).
    """
    resolved = []
    need_blast = []
    for r in gmax_records:
        gid = r.get("gene_id") or ""
        m = re.match(r"GLYMA_(\d+G\d+)v\d+$", gid, re.IGNORECASE)
        if m:
            resolved.append({**r, "glyma_id": f"Glyma.{m.group(1)}", "resolution_method": "direct_tag"})
        else:
            need_blast.append(r)

    if need_blast:
        hits = _blast_query_soybase(need_blast)
        hit_by_odb = {h["odb_id"]: h for h in hits}
        for r in need_blast:
            h = hit_by_odb.get(r["odb_id"])
            if h and h["glyma_id"]:
                pct_id = (h["identity"] / h["query_length"] * 100) if h["identity"] else None
                resolved.append({**r, "glyma_id": h["glyma_id"], "resolution_method": "blastp_soybase",
                                  "blast_pct_identity": pct_id, "blast_bit_score": h["bit_score"]})
            else:
                resolved.append({**r, "glyma_id": None, "resolution_method": "unresolved"})

    return resolved


def compute_tissue_cpm(tissue_dfs: dict, glyma_ids: list) -> pd.DataFrame:
    """Tissue-level CPM (sum across all cell states in a tissue) for each
    Glyma ID. Genes absent from a tissue's matrix get CPM 0."""
    rows = {}
    for tissue, df in tissue_dfs.items():
        lib_size = df.sum(axis=0).sum()
        row = {}
        for gid in glyma_ids:
            row[gid] = (df.loc[gid].sum() / lib_size * 1e6) if gid in df.index else 0.0
        rows[tissue] = row
    return pd.DataFrame(rows).T


def tau(expr_vec) -> float:
    """Yanai et al. 2005 tissue-specificity index."""
    x = np.asarray(expr_vec, dtype=float)
    if x.max() <= 0:
        return float("nan")
    xhat = x / x.max()
    n = len(x)
    return float(np.sum(1 - xhat) / (n - 1))


def compute_functional_calls(tissue_cpm: pd.DataFrame, nonfunctional_thresh=3.0,
                              tissue_specific_tau=0.7, tissue_specific_min_expr=5.0) -> pd.DataFrame:
    """Per-gene tau, expression breadth, and functional/tissue-specific calls,
    following Benoit et al. 2025's exact thresholds."""
    out = {}
    for gid in tissue_cpm.columns:
        expr = tissue_cpm[gid].values
        t = tau(expr)
        breadth = int((tissue_cpm[gid] > nonfunctional_thresh).sum())
        avg_expr = tissue_cpm[gid].mean()
        peak_tissue = tissue_cpm[gid].idxmax()
        peak_expr = tissue_cpm[gid].max()
        nonfunctional = bool(avg_expr < nonfunctional_thresh)
        tissue_specific = bool((not np.isnan(t)) and t > tissue_specific_tau and peak_expr > tissue_specific_min_expr)
        out[gid] = {
            "tau": round(t, 4) if not np.isnan(t) else None,
            "expression_breadth_tissues_CPM_gt_3": breadth,
            "avg_CPM_across_tissues": round(float(avg_expr), 4),
            "peak_tissue": peak_tissue,
            "peak_CPM": round(float(peak_expr), 4),
            "nonfunctional_by_avg_expr_rule": nonfunctional,
            "tissue_specific_call": tissue_specific,
        }
    return pd.DataFrame(out).T


def _best_shared_tissue(tissue_dfs: dict, gene_a: str, gene_b: str, min_celltypes=3):
    """Pick the tissue with the most cell-type columns among tissues where
    both genes are present in the matrix index, for pairwise classification."""
    best = None
    for tissue, df in tissue_dfs.items():
        if gene_a in df.index and gene_b in df.index:
            n_ct = df.shape[1]
            if best is None or n_ct > best[1]:
                best = (tissue, n_ct)
    if best is None or best[1] < min_celltypes:
        return None
    return best[0]


def classify_paralog_pair(tissue_dfs: dict, gene_a: str, gene_b: str):
    """Benoit et al. 2025 4-group classification for one paralog pair, using
    the tissue with the richest cell-type resolution where both are present."""
    tissue = _best_shared_tissue(tissue_dfs, gene_a, gene_b)
    if tissue is None:
        return {"gene_a": gene_a, "gene_b": gene_b, "classifiable": False,
                "reason": "no shared tissue with >=3 cell types where both genes are detected"}

    df = tissue_dfs[tissue]
    lib_size = df.sum(axis=0)
    cpm = df.div(lib_size, axis=1) * 1e6

    col_medians = cpm.median(axis=0)
    keep_mask = cpm.gt(col_medians, axis=1).any(axis=1)
    filtered = cpm.loc[keep_mask]
    if gene_a not in filtered.index:
        filtered = pd.concat([filtered, cpm.loc[[gene_a]]])
    if gene_b not in filtered.index:
        filtered = pd.concat([filtered, cpm.loc[[gene_b]]])
    logmat = np.log2(filtered + 1)

    def corr_vector(target_gene, mat):
        X = mat.values
        idx = mat.index.get_loc(target_gene)
        x = X[idx]
        x_c = x - x.mean()
        X_c = X - X.mean(axis=1, keepdims=True)
        num = X_c @ x_c
        denom = np.sqrt((X_c ** 2).sum(axis=1) * (x_c ** 2).sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            r = num / denom
        return pd.Series(r, index=mat.index)

    corr_a = corr_vector(gene_a, logmat)
    corr_b = corr_vector(gene_b, logmat)

    def rank_standardize(s):
        n = len(s)
        ranks = s.rank(method="average", na_option="keep")
        median_rank = (n + 1) / 2.0
        ranks[s.isna()] = median_rank
        return ranks / ranks.max()

    rank_a = rank_standardize(corr_a)
    rank_b = rank_standardize(corr_b)
    coexpr = (rank_a[gene_b] + rank_b[gene_a]) / 2.0

    paired = cpm.loc[[gene_a, gene_b]].T
    nz = (paired[gene_a] > 0) & (paired[gene_b] > 0)
    if nz.sum() < 2:
        return {"gene_a": gene_a, "gene_b": gene_b, "tissue_used": tissue, "classifiable": False,
                "reason": "fewer than 2 samples with nonzero expression of both genes"}
    log2fc = np.log2(paired.loc[nz, gene_a] / paired.loc[nz, gene_b])
    mean_abs_lfc = float(log2fc.abs().mean())
    sd_abs_lfc = float(log2fc.abs().std()) if nz.sum() > 1 else 0.0

    if coexpr > 0.9 and mean_abs_lfc < 1 and sd_abs_lfc < 1:
        group = "I: Dosage balanced"
    elif coexpr > 0.9 and mean_abs_lfc >= 1 and sd_abs_lfc < 1:
        group = "II: Paralogue dominance"
    elif coexpr > 0.9 and mean_abs_lfc >= 1 and sd_abs_lfc >= 1:
        group = "III: Specialized"
    elif coexpr < 0.5 and mean_abs_lfc >= 1 and sd_abs_lfc >= 1:
        group = "IV: Diverged"
    else:
        group = "unclassified"

    return {
        "gene_a": gene_a, "gene_b": gene_b, "classifiable": True, "tissue_used": tissue,
        "n_samples_nonzero_both": int(nz.sum()), "n_samples_total": int(len(paired)),
        "genes_in_coexpression_network": int(filtered.shape[0]),
        "rank_standardized_coexpression": round(float(coexpr), 4),
        "mean_abs_log2FC": round(mean_abs_lfc, 4), "sd_abs_log2FC": round(sd_abs_lfc, 4),
        "benoit_group": group,
    }


def all_celltype_metrics(tissue_dfs: dict, gene_a: str, gene_b: str) -> dict:
    """Pearson/Spearman/log2FC for a pair, pooled across ALL cell types from
    ALL tissues (not just the best-shared tissue)."""
    rows = []
    for tissue, df in tissue_dfs.items():
        lib_size = df.sum(axis=0)
        cpm = df.div(lib_size, axis=1) * 1e6
        va = cpm.loc[gene_a] if gene_a in cpm.index else pd.Series(0.0, index=cpm.columns)
        vb = cpm.loc[gene_b] if gene_b in cpm.index else pd.Series(0.0, index=cpm.columns)
        for ct in cpm.columns:
            rows.append({"tissue": tissue, "cell_type": ct, "a": va[ct], "b": vb[ct]})
    all_ct = pd.DataFrame(rows)

    pearson_r, pearson_p = pearsonr(all_ct["a"], all_ct["b"])
    spearman_r, spearman_p = spearmanr(all_ct["a"], all_ct["b"])
    nz = (all_ct["a"] > 0) & (all_ct["b"] > 0)
    log2fc = np.log2(all_ct.loc[nz, "a"] / all_ct.loc[nz, "b"]) if nz.sum() > 0 else pd.Series(dtype=float)

    return {
        "gene_a": gene_a, "gene_b": gene_b,
        "n_cell_types_total": int(len(all_ct)), "n_cell_types_nonzero_both": int(nz.sum()),
        "pearson_r_raw_CPM": round(float(pearson_r), 4), "pearson_p_raw_CPM": float(pearson_p),
        "spearman_r_raw_CPM": round(float(spearman_r), 4), "spearman_p_raw_CPM": float(spearman_p),
        "mean_abs_log2FC": round(float(log2fc.abs().mean()), 4) if len(log2fc) else None,
        "sd_abs_log2FC": round(float(log2fc.abs().std()), 4) if len(log2fc) > 1 else None,
    }


def run_soybean_expression_analysis(records: list, cache_dir: str) -> dict:
    """Orchestrated entry point: given a full list of OrthoDB group records
    (any species mix) and a cache directory, run the whole soybean
    expression-redundancy analysis. Returns a summary dict; skips cleanly
    (with a 'skipped' reason) if the group has no Glycine max records."""
    gmax_records = [r for r in records if is_glycine_max_record(r)]
    if not gmax_records:
        return {"skipped": True, "reason": "no Glycine max records in this OrthoDB group"}

    resolved = resolve_glyma_ids(gmax_records)
    glyma_ids = [r["glyma_id"] for r in resolved if r["glyma_id"]]
    if not glyma_ids:
        return {"skipped": True, "reason": "Glycine max records present but none resolved to a Glyma ID",
                "gmax_records": gmax_records}

    tissue_dfs = fetch_gse270392_tissue_matrices(cache_dir)
    tissue_cpm = compute_tissue_cpm(tissue_dfs, glyma_ids)
    functional_calls = compute_functional_calls(tissue_cpm)

    detected_ids = [g for g in glyma_ids if tissue_cpm[g].max() > 0]
    pair_classifications = []
    all_celltype_pairs = []
    for i in range(len(detected_ids)):
        for j in range(i + 1, len(detected_ids)):
            pair_classifications.append(classify_paralog_pair(tissue_dfs, detected_ids[i], detected_ids[j]))
            all_celltype_pairs.append(all_celltype_metrics(tissue_dfs, detected_ids[i], detected_ids[j]))

    return {
        "skipped": False,
        "n_gmax_records": len(gmax_records),
        "resolved_records": resolved,
        "glyma_ids": glyma_ids,
        "detected_glyma_ids": detected_ids,
        "tissue_cpm": tissue_cpm,
        "functional_calls": functional_calls,
        "pair_classifications": pair_classifications,
        "all_celltype_pair_metrics": all_celltype_pairs,
    }
