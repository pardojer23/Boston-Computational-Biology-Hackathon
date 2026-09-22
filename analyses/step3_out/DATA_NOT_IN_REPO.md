# Step 3 — data not in this repository

## Absent because it does not exist publicly

**Cotyledon replicate-split counts.** The source repo ships
`Gm_atlas_Cotyledon_stage_seeds_RNA_gene_celltype_splitReps_counts.txt` as a
**2-byte (empty) file**, and it is not obtainable elsewhere:

- GSE270392 has no `splitReps` or `celltype` supplementary file for any tissue.
- Its `GSE270392_Gm_atlas_Cotyledon_stage_seeds_RNA_gene_cellstate_counts.txt.gz`
  is a **different aggregation** — 12 columns named by cell-state cluster
  (`Cortex.0`, `Phloem.8_0`) rather than replicate x cell type
  (`RNA_rt3:Cortex`). Substep 1 correlates two replicates per cell type, which
  that layout cannot support. Same 23,776 genes; incompatible columns.

Consequence: step 3 covers 6 of 7 tissues; steps 4 and 7 lose their cotyledon
analyses; the cross-stage stability score covers 3 of 4 stages.

Possible remedy, not yet done: rebuild replicate-split pseudobulk from
`GSE270392_Gm_atlas_Cotyledon_stage_seeds.rna.seurat.obj.rds.gz` (~1 GB), whose
per-cell metadata carries library identity. The method is testable — rebuild
Root the same way and compare against Root's intact shipped matrix first. It
would be a declared reconstruction, not a reproduction.

## Excluded because it is cheaply regenerable

Both are produced in seconds from repo-shipped inputs by the scripts in `code/`:

| not committed | size | regenerate with |
|---|---|---|
| `<tissue>_reorder.csv`, `<tissue>_reorder_auto.csv` (6 tissues) | 8.5 MB packed | `Rscript code/step3_01_reorder_qc.R` |
| `<tissue>_filter_normalized.csv` (6 tissues) | 14.0 MB packed | `Rscript code/step3_02_normalize_edger.R` |

`results/step3_expressed_pairs.tar.gz` **is** committed (2.7 MB) because the
co-expression job consumes it directly.

Per-pair extraction output (`<tissue>_filter_pairs.csv`, 362k-724k rows each) is
also not committed; regenerate with `code/step3_03_extract_and_filter.sh`.

## Third-party inputs, not redistributed here

| what | where |
|---|---|
| per-tissue count matrices | shipped in the source repo, `step3_divergence_expression_pattern/raw_data/` |
| ACR cell-state count files (10 tissues) | GSE270392 supplementary, `*_ACR_4cpm_cellstates_counts.txt.gz` |
