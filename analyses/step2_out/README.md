# Step 2 — Duplicated gene sets

A "set" is a legume gene family containing exactly four WGD-derived soybean
genes, split into two clades of two by the family tree.

## Result

**1,920 sets reproduced** against 1,904 published. 1,896 of the 1,904 published
families were recovered (99.58%), and for **every one of those 1,896 shared
families the gene quartet is identical** — zero disagreements on content.

The 32 differing families are fully attributed to step 1: all 24 extras contain
a gene that is WGD-derived in our step-1 output but not the authors', and all 8
missing contain a gene WGD-derived in theirs but not ours. Our WGD gene universe
is 35,493 against their 35,174.

## Reproduce

```bash
# needs: perl, python3, curl; ~15 MB download; ~10 min on 4 cores
python3 code/step2_local.py
```

Inputs are fetched from
`data.soybase.org/LEGUMES/Fabaceae/genefamilies/legume.fam3.VLMQ`
(`sup1A_hsh.tsv.gz`, 732,283 rows; `sup1B_trees.tar.gz`, 24,811 tree files) and
cached under `step2_work/`. The WGD gene list in `results/` comes from step 1.

Runs locally — this step is ~15 MB of input and a tree walk, not cluster work.
The tree walk is a verbatim port of the authors' perl expressions from
`duplicated_gene_sets_age_server.sh`; family filtering follows
`07042025_tree_infor_glyma.Rmd`.

## Files

| path | what |
|---|---|
| `code/step2_local.py` | the pipeline, end to end |
| `results/reproduced_gene_sets.csv` | the 1,920 sets — input to steps 4 and 8 |
| `results/step2_glyma_output_sets.txt` | raw tree-walk output (2,656 rows) before family filtering |
| `results/wgd_unique_genes_lexicographic.csv` | the 35,493 WGD genes from step 1 |
| `reports/step2_validation_report.md` | funnel, concordance, residual attribution |
| `reports/step2_run_log.txt` | run log |
