# Step 3 — Expression divergence between duplicated genes

**Substeps 1-3 of 6 complete, for 6 of 7 tissues.** Cotyledon cannot be
processed (see `DATA_NOT_IN_REPO.md`).

## What is done

| substep | script | outcome |
|---|---|---|
| 1. reorder + replicate QC | `code/step3_01_reorder_qc.R` | `results/reorder_qc_summary.csv`, `results/step3_celltype_qc_decisions.csv` |
| 2. edgeR normalization | `code/step3_02_normalize_edger.R` | `results/normalization_summary.csv` |
| 3. pair extraction + expressed filter | `code/step3_03_extract_and_filter.sh` | `results/expressed_pairs_summary.csv`, `results/step3_expressed_pairs.tar.gz` |

Remaining: co-expression networks (`code/dupgen_step3_coexp_modal.py`, needs
16-64 GB so it runs on Modal), then fold change, expression-type
classification, and the cross-tissue stability analysis.

## Two findings worth knowing before reusing this

**The documented QC rule is not what was run.** The step notes describe keeping
cell types with replicate Spearman rho > 0.75. That reproduces the authors'
manual column selection for only 3 of 6 tissues. For early_maturation, globular
and hypocotyl they additionally dropped cell types whose rho was *above*
threshold (0.832, 0.801/0.829, 0.769) — in every case an `Unknown` /
`SC_unknown` unannotated cluster. With "discard unannotated clusters" added, all
ten dropped cell types are explained and none is left unaccounted for; see
`results/step3_celltype_qc_decisions.csv`.

The notes attribute this second filter to ATAC cell-type consistency. That is
contradicted by the data: `Dividing_cell`, `Nodule_meritem` and
`SC_vascular_parenchyma` are absent from the ATAC cell-state files yet retained.

**A column-indexing bug in the reorder vectors.** For globular and hypocotyl
(and cotyledon), `reorder_gene_counts.Rmd` ends its index vector
`...,12,25,23,26` where `13` was intended. Index 23 is repeated, so the final
"replicate pair" compares two *different* cell types' replicate-2 columns, and
one cell type's replicate-1 column is never used. The authors' own column
selection discards that pair, so it does not reach their downstream data — but
do not copy the vector.

## Validation signal

The retained cell-type counts (10, 6, 11, 12, 11, 6) match the hand-typed group
vectors in `normalization_edgeR.Rmd` exactly for all six tissues. Those lists
were written by hand after the authors' manual QC, so matching every one is
independent evidence the column selection reproduces theirs.

## Reproduce

```bash
Rscript code/step3_01_reorder_qc.R <repo>/step3_divergence_expression_pattern/raw_data step3_out
Rscript code/step3_02_normalize_edger.R step3_out          # needs edgeR
bash    code/step3_03_extract_and_filter.sh ../step1_out/results/reproduced_pairs_lexicographic.csv <repo> step3_out
```
