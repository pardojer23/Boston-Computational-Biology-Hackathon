"""
Fetch the GSE270392 (soybean spatially-resolved multiomic single-cell atlas)
series metadata and its released per-tissue gene x cell-state pseudobulk RNA
count matrices from GEO.

Usage:
    python fetch_gse270392.py
"""
import gzip
import io
import json
import re
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SUPPL_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270392/suppl/"

# Tissues for which GEO released a small, pre-aggregated
# "<tissue>_RNA_gene_cellstate_counts.txt.gz" pseudobulk matrix (gene rows x
# cell-state columns). Leaf and Pod were profiled by scATAC-seq only (no
# snRNA-seq gene expression matrix was released for those two tissues).
TISSUE_FILES = {
    "Early_nodule": "GSE270392_Gm_atlas_Early_nodule_RNA_gene_cellstate_counts.txt.gz",
    "Root": "GSE270392_Gm_atlas_Root_RNA_gene_cellstate_counts.txt.gz",
    "Hypocotyl": "GSE270392_Gm_atlas_Hypocotyl_RNA_gene_cellstate_counts.txt.gz",
    "Cotyledon_stage_seeds": "GSE270392_Gm_atlas_Cotyledon_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Early_maturation_stage_seeds": "GSE270392_Gm_atlas_Early_maturation_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Globular_stage_seeds": "GSE270392_Gm_atlas_Globular_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
    "Heart_stage_seeds": "GSE270392_Gm_atlas_Heart_stage_seeds_RNA_gene_cellstate_counts.txt.gz",
}


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "research-agent"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def fetch_series_summary(accession: str = "GSE270392") -> dict:
    """Confirm the series accession and pull its sample list via NCBI E-utils."""
    esearch = f"{EUTILS}/esearch.fcgi?db=gds&term={accession}[ACCN]&retmode=json"
    uid = json.loads(fetch(esearch))["esearchresult"]["idlist"][0]
    esummary = f"{EUTILS}/esummary.fcgi?db=gds&id={uid}&retmode=json"
    return json.loads(fetch(esummary))["result"][uid]


def download_tissue_matrix(tissue: str, filename: str, out_dir: str = ".") -> str:
    out_path = f"{out_dir}/{filename}"
    raw = fetch(SUPPL_BASE + filename)
    with open(out_path, "wb") as f:
        f.write(raw)
    return out_path


def main() -> None:
    summary = fetch_series_summary()
    print("Series:", summary["accession"], "-", summary.get("n_samples"), "samples")
    print("Associated PubMed IDs:", summary.get("pubmedids"))

    for tissue, filename in TISSUE_FILES.items():
        path = download_tissue_matrix(tissue, filename)
        print(f"Downloaded {tissue} -> {path}")


if __name__ == "__main__":
    main()
