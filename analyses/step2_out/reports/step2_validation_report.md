# Step 2 reproduction: duplicated gene SETS

Reproduction of step 2 of Li, Zhang & Schmitz 2025, *From duplication to divergence*
(Plant Cell 37(12), koaf279), following `step2.identify_duplicated_gene_sets` in
`Xianglichina/soybean_duplicated_genes_single_cell`.

A "set" is a legume gene family containing exactly four WGD-derived soybean genes,
split into two clades of two by the family tree.

## Execution

Ran **locally**, not on Modal. The step is ~15 MB of input and a per-family tree
walk — seconds to minutes of CPU, so it falls on the local side of the agreed
split (light steps local, heavy compute on Modal). An earlier Modal
implementation (`dupgen_step2_modal.py`) is redundant.

| component | value |
|---|---|
| inputs | `legume.fam3.VLMQ.sup1A_hsh.tsv.gz` (732,283 rows), `legume.fam3.VLMQ.sup1B_trees.tar.gz` (24,811 tree files) |
| source | data.soybase.org/LEGUMES/Fabaceae/genefamilies/legume.fam3.VLMQ |
| tree walk | verbatim port of `duplicated_gene_sets_age_server.sh` (same perl expressions) |
| filtering | port of `07042025_tree_infor_glyma.Rmd` |
| WGD gene input | 35,493 genes from step 1 (lexicographic variant) |

## Funnel

| stage | count |
|---|---|
| family tree files | 24,811 |
| families containing soybean genes | 18,301 |
| families with exactly two soybean clades (tree walk) | 2,656 |
| families with exactly 4 soybean genes | 2,964 |
| families whose 4 soybean genes are all WGD-derived | 2,344 |
| **final duplicated gene sets** | **1,920** |

## Concordance with the shipped table

Target: `07162025_duplicated_gene_set_updated.csv`, 1,904 families.

| measure | value |
|---|---|
| families reproduced | 1,920 (published 1,904) |
| published families recovered | 1,896 / 1,904 = **99.58%** |
| shared families with an identical gene quartet | **1,896 / 1,896 = 100.00%** |
| shared families with a differing quartet | 0 |
| unique genes | 7,680 (published 7,616), overlap 7,584 |
| extra families | 24 |
| missing families | 8 |

The second row is the important one: for every family present in both tables, the
four genes and their clade assignment are **identical**. The tree walk and the
family filter reproduce the authors' procedure exactly; there is no quartet on
which we disagree.

## The residual is inherited from step 1, not introduced here

All 32 differing families are explained by the step-1 WGD gene set, with no
residue:

- **24/24** extra families contain at least one gene that is WGD-derived in our
  step-1 output but not in the authors' (e.g. `Legume.fam3.00606` ->
  `Glyma.14G154200`).
- **8/8** missing families contain at least one gene that is WGD-derived in the
  authors' output but not in ours (e.g. `Legume.fam3.04828` ->
  `Glyma.08G308700`, `Glyma.18G107900`).

Our WGD gene universe is 35,493 genes against the authors' 35,174 derived from
their shipped pairs table (542 ours-only, 223 theirs-only). Step 1 reproduced
99.42% of published pairs; the 0.58% shortfall propagates here as 32 families out
of ~1,900, i.e. the step-2 discrepancy is a downstream consequence of a known
upstream difference rather than an independent error.

## Deviations

1. Ran locally rather than on Modal (justified above; compute is trivial).
2. The tree walk and filtering were reimplemented in Python around the authors'
   perl one-liners rather than run from their shell script, which hardcodes
   server paths. The perl expressions themselves are unchanged.
3. Step 1's lexicographic filter variant was used as the WGD input, being the
   closer match to the published pairs.
