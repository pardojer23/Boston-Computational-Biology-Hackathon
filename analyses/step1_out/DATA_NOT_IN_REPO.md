# Data not tracked in this repository

## Nothing from step 1 was too large for GitHub

To be precise about this: every step-1 output fits well inside GitHub's limits
(warning at 50 MB, hard block at 100 MB per file). The largest tracked file is
6.85 MB. Nothing was dropped for size.

| tracked file | size |
|---|---|
| `results/raw/dupgen_step1_dupgen_finder_output.tar.gz` | 6.85 MB |
| `results/reproduced_pairs_numeric.csv` | 6.08 MB |
| `results/reproduced_pairs_lexicographic.csv` | 6.04 MB |
| `results/reproduced_genes_numeric.csv` | 0.74 MB |
| `results/reproduced_genes_lexicographic.csv` | 0.73 MB |

The items below are absent for other reasons — they are inputs, container-local
intermediates, or third-party source data. All are regenerable or re-downloadable
with the commands given.

## 1. Reference inputs (~21 MB, not tracked by choice)

Regenerable; also matched by the `.gitignore` data patterns. The pipeline
downloads these inside the Modal container, so they never need to exist locally.

```bash
B=https://data.soybase.org
G=$B/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ
P=$B/Phaseolus/vulgaris/annotations/G19833.gnm2.ann1.PB8d
curl -O $G/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz
curl -O $G/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.bed.gz
curl -O $P/phavu.G19833.gnm2.ann1.PB8d.protein_primary.faa.gz
curl -O $P/phavu.G19833.gnm2.ann1.PB8d.gene_models_main.bed.gz
```

| file | size | md5 |
|---|---|---|
| glyma...protein_primary.faa.gz | 11,347,538 | cf2b18ce6a05a01f1a2fd14b9e827e17 |
| glyma...gene_models_main.bed.gz | 888,675 | 951c4991ef458457a7b9d2d30f765ded |
| phavu...protein_primary.faa.gz | 7,865,273 | da23aafaeb1800194f72b20935f3055a |
| phavu...gene_models_main.bed.gz | 425,554 | 2c33e64f150c5fdead67a6e5d21147c7 |

Counts: 52,872 soybean primary proteins (199 genes on scaffolds), 27,433 common
bean primary proteins (421 on scaffolds).

## 2. DIAMOND intermediates (container-local, never harvested)

These stayed inside the Modal container and were deliberately not returned, since
only `./out/` is harvested. Regenerate by rerunning the pipeline.

| intermediate | scale |
|---|---|
| `glyma.blast` | 216,462 hits |
| `glyma_phavu.blast` | 177,637 hits |
| filtered, numeric | 178,542 soybean / 123,172 bean hits |
| filtered, lexicographic | 124,817 soybean / 122,985 bean hits |
| `glyma.dmnd`, `phavu.dmnd` | DIAMOND databases |

## 3. Uncompressed DupGen_finder output (49 MB)

Tracked as a 6.85 MB tarball rather than 30 loose files. Unpack with:

```bash
tar -xzf results/raw/dupgen_step1_dupgen_finder_output.tar.gz
```

Contains, for each filter variant: `glyma.{wgd,tandem,proximal,transposed,dispersed}.{pairs,genes}`,
`glyma.collinearity`, `glyma_phavu.collinearity`, `glyma.gff.sorted`,
`glyma.singletons`, `glyma.pairs.stats`.

## 4. GSE270392 single-cell data (needed from step 3 onward)

Third-party source data, far too large to track — the largest single file is
9.6 GB. Download per tissue from
[GSE270392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270392).

Gene **expression** lives in:

- `GSE270392_Gm_atlas_<tissue>_RNA_gene_cellstate_counts.txt.gz` — 360–821 KB,
  plain-text genes x cell states (7 tissues)
- `GSE270392_Gm_atlas_<tissue>.rna.seurat.obj.rds.gz` — 423 MB–1.3 GB Seurat
  objects with per-cell counts, metadata and embeddings (7 tissues)
- `GSE270392_Gm_<tissue>.spRNA.obj.rds.gz` — spatial transcriptomics (5 tissues)

**Do not use `GSE270392_*_metav3_gene_sparse.rds.gz` as expression data.** Despite
`gene` in the filename it is the scATAC **gene accessibility** matrix. It exists
for 10 tissues whereas the RNA files exist for 7 — Leaf, Pod and
Middle_maturation_stage_seeds have chromatin data only, so an expression matrix
for them cannot exist. The scATAC sample records confirm it, listing exactly four
supplementary types of which this is the "gene accessibility matrix".

Tissues with expression data: Cotyledon_stage_seeds, Early_maturation_stage_seeds,
Early_nodule, Globular_stage_seeds, Heart_stage_seeds, Hypocotyl, Root.

Note that GEO samples carry no supplementary files of their own
(`!Sample_supplementary_file_1 = NONE`); raw reads are in SRA and all processed
matrices are series-level, so downloads are per tissue rather than per sample.
