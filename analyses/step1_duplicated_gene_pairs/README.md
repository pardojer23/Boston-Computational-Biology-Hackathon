# Step 1 — Identify duplicated gene pairs in soybean

Reproduction of step 1 of Li, Zhang & Schmitz (2025), *From duplication to
divergence: single-cell insights into transcriptional and cis-regulatory
landscapes in soybean*, The Plant Cell 37(12): koaf279
([doi](https://doi.org/10.1093/plcell/koaf279)).

Protocol followed: [`step1.identify_duplicated_gene_pairs`](https://github.com/Xianglichina/soybean_duplicated_genes_single_cell)
from the authors' repository. Single-cell data for later steps comes from
[GSE270392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270392).

## Result

59,004 of the 59,351 published pairs recovered — **99.42%** concordance.

| | numeric filter | lexicographic filter | published |
|---|---|---|---|
| total pairs | 60,791 | 60,371 | 59,351 |
| unique genes | 45,348 | 45,168 | 45,161 |
| published pairs recovered | 58,778 (99.03%) | **59,004 (99.42%)** | — |

**Use `results/reproduced_pairs_lexicographic.csv` for downstream steps.** It is
the closer match on every measure, and its unique-gene count lands within 7 genes
of the published 45,161.

### Why two variants

The authors' `07032025_blast_filter.Rmd` reads BLAST output with `read.delim2()`,
which sets `dec=","`. Percent identity therefore stays a *character* column, so
`filter(identity > 40)` performs a string comparison rather than a numeric one —
it keeps `"9.5"` (because "9" > "4") and discards `"100"` (because "1" < "4").
The code alone cannot tell us which behaviour produced the published numbers, so
both were run. The lexicographic variant matching more closely is evidence the
published figures came from the script as written.

See `reports/step1_validation_report.md` for full per-category tables and the
analysis of the residual ~1.7%.

## Layout

```
code/dupgen_step1_modal.py      pipeline (Modal app, end to end)
results/reproduced_pairs_*.csv  duplicated gene pairs, per filter variant
results/reproduced_genes_*.csv  unique duplicated genes, per filter variant
results/raw/*.tar.gz            raw DupGen_finder output, both variants
reports/step1_validation_report.md   methods, concordance, deviations
reports/modal_run_log.txt       run log (DIAMOND version, hit counts, timings)
DATA_NOT_IN_REPO.md             inputs and intermediates not tracked here
```

## Reproducing

Runs on Modal (16 CPU / 32 GB container); the reference inputs are downloaded
inside the container, so no local data staging is needed.

```bash
conda activate modal            # env with the modal client installed
modal run code/dupgen_step1_modal.py
```

Outputs land in `modal_out/{numeric,lexicographic}/`. Assemble them into the
pairs and genes CSVs with the scaffold-removal and deduplication described in
the validation report.

| component | version |
|---|---|
| DIAMOND | v2.1.8 (pinned — the authors' version) |
| DupGen_finder | qiao-xin/DupGen_finder, compiled from source |
| soybean | glyma.Wm82.gnm4.ann1.T8TQ, 52,872 primary proteins |
| common bean (outgroup) | phavu.G19833.gnm2.ann1.PB8d, 27,433 primary proteins |

## Two corrections to the published protocol

1. `glyma_phavu.gff` must contain gene positions for **both** species, not just
   the outgroup. The step-1 notes omit this.
2. The filenames in those notes (`glyma.blast`, `phavu.blast`) are not what
   DupGen_finder accepts; it requires `glyma.blast` and `glyma_phavu.blast`.

A literal transcription of the written protocol fails at the DupGen_finder call.
