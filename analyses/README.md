# Reproduction of Li, Zhang & Schmitz 2025 (Plant Cell, koaf279)

Step-by-step reproduction of *From duplication to divergence: single-cell
insights into transcriptional and cis-regulatory landscapes in soybean*
([doi](https://doi.org/10.1093/plcell/koaf279)), following the protocol in
[`Xianglichina/soybean_duplicated_genes_single_cell`](https://github.com/Xianglichina/soybean_duplicated_genes_single_cell).
Single-cell data: [GSE270392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270392).

| step | what it produces | status | concordance |
|---|---|---|---|
| [step1_out](step1_out/) | duplicated gene pairs (DupGen_finder) | complete | 59,004 / 59,351 published pairs = **99.42%** |
| [step2_out](step2_out/) | duplicated gene sets (4 WGD genes per family) | complete | 1,896 / 1,904 published families = **99.58%**, and 100% identical quartets among those |
| [step3_out](step3_out/) | expression divergence per tissue | substeps 1-3 of 6 | no published target for intermediates |
| steps 4-9 | see the roadmap below | not started | — |

Each step directory holds `code/`, `results/`, `reports/`, and a
`DATA_NOT_IN_REPO.md` where large or third-party inputs are documented instead
of committed.

## Known limitations carried forward

1. **Cotyledon is unavailable for step 3.** The shipped
   `Gm_atlas_Cotyledon_stage_seeds_RNA_gene_celltype_splitReps_counts.txt` is 2
   bytes, and GEO ships no replicate-split equivalent (its `cellstate` file is a
   different aggregation). Step 3 therefore covers 6 of 7 tissues, and the
   cross-stage stability analysis would cover 3 of 4 developmental stages.
2. **Four scripts the step notes name are absent from the source repo**, two of
   which gate later steps: `07102025_ACR_associated_genes.Rmd` (step 7) and
   `0902_characterization_gene_set_expression_within_tissue.Rmd` (step 4).
3. **Step 9's two differential-test inputs** exist in neither the source repo nor
   the 62 GEO supplementary files, which carry counts and objects rather than
   test results. Step 9 needs either author contact or a declared substitution.
4. **The paper's Figures 1-6 have no deposited code.** Only a few component
   plots are scripted, two of them inside the blocked step 9.
