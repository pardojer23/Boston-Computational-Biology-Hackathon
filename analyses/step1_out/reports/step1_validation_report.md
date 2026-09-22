# Step 1 reproduction: soybean duplicated gene pairs

Reproduction of step 1 of Li, Zhang & Schmitz 2025, *From duplication to divergence*
(Plant Cell 37(12), koaf279), following the protocol in
`Xianglichina/soybean_duplicated_genes_single_cell/step1.identify_duplicated_gene_pairs`.

## Execution

Ran on Modal (16 CPU / 32 GB container), not locally.

| component | value |
|---|---|
| DIAMOND | v2.1.8 (pinned to the authors' version) |
| DupGen_finder | qiao-xin/DupGen_finder, compiled from source |
| soybean | glyma.Wm82.gnm4.ann1.T8TQ, 52,872 primary proteins |
| common bean (outgroup) | phavu.G19833.gnm2.ann1.PB8d, 27,433 primary proteins |
| BLASTp parameters | `--sensitive --max-target-seqs 5 --evalue 1e-10` |

DIAMOND hits: 216,462 soybean self, 177,637 soybean vs bean.

## The identity-filter fork

The authors' `07032025_blast_filter.Rmd` reads BLAST output with `read.delim2()`,
which sets `dec=","`. Percent identity therefore stays a *character* column, so
`filter(identity > 40)` is a string comparison, not a numeric one: it keeps
`"9.5"` (because "9" > "4") and discards `"100"` (because "1" < "4"). Since the
code alone cannot settle which behaviour produced the published numbers, both were
run.

| filter | soybean hits kept | bean hits kept |
|---|---|---|
| numeric (`> 40` as a number) | 178,542 | 123,172 |
| lexicographic (R's actual behaviour) | 124,817 | 122,985 |

## Concordance with the published result

Published: 59,351 pairs over 45,161 unique genes.

### Numeric filter

| category | reproduced | published | diff |
|---|---|---|---|
| WGD | 31,999 | 29,868 | +2,131 |
| TD  | 2,465  | 1,749  | +716 |
| PD  | 1,441  | 1,942  | -501 |
| TRD | 3,075  | 3,275  | -200 |
| DSD | 21,811 | 22,517 | -706 |
| **total** | **60,791** | **59,351** | **+1,440** |
| unique genes | 45,348 | 45,161 | +187 |

Pair-set overlap: 58,778 / 59,351 unordered (**99.03%**); 2,013 pairs found that
are not published, 573 published pairs not recovered.

### Lexicographic filter

| category | reproduced | published | diff |
|---|---|---|---|
| WGD | 31,746 | 29,868 | +1,878 |
| TD  | 2,414  | 1,749  | +665 |
| PD  | 1,413  | 1,942  | -529 |
| TRD | 3,079  | 3,275  | -196 |
| DSD | 21,719 | 22,517 | -798 |
| **total** | **60,371** | **59,351** | **+1,020** |
| unique genes | 45,168 | 45,161 | +7 |

Pair-set overlap: 59,004 / 59,351 unordered (**99.42%**); 1,367 pairs found that
are not published, 347 published pairs not recovered.

The lexicographic variant is closer on every measure — total pairs, unique gene
count (+7 vs +187), and pair overlap — which is weak evidence that the published
numbers came from the R script as written, dec="," bug included. It is not
conclusive, because neither variant reproduces 59,351 exactly.

## Most likely source of the residual ~1.7%

The per-category pattern is not random: **tandem is inflated (+665) while proximal
is deflated (-529)** in both variants. Those two categories are distinguished only
by how many genes sit between the two copies (proximal = up to 10 intervening
genes). A denser gene annotation pushes pairs from tandem into proximal; a sparser
one does the reverse. Our inflated tandem / deflated proximal is exactly what a
*sparser* gene universe predicts.

Independent support: the authors' own `07032025soybean_duplicated_genes.Rmd` loads
an `all_gene` table of **66,070** genes, whereas `gene_models_main` primary
transcripts give **52,872**. That file was not published, and 66,070 matches no
current annotation product (primary + low-confidence = 58,842, checked directly),
so the exact gene set they fed to DupGen_finder cannot be reconstructed from the
repository.

Note also that the same notebook hardcodes `rep("WGD", 29905)` where the shipped
output has 29,868 WGD rows — so some constants in it are stale, and small
disagreements with numbers appearing only in that file should not be read as
reproduction failures.

## Coordinate-level confirmation

Independent of counts, the pipeline is demonstrably operating on the same reference
build: `Glyma.01G017200` maps to `glyma.Wm82.gnm4.Gm01:1645218`, byte-identical to
the `Location` string in the authors' shipped `07042025_glyma_duplicated_pairs.csv`.

## Deviations from the repo's written protocol

1. `glyma_phavu.gff` must contain gene positions for **both** species. The step-1
   notes do not say this, and the filenames they use (`glyma.blast`, `phavu.blast`)
   are not what DupGen_finder accepts (`glyma.blast`, `glyma_phavu.blast`).
2. The BLAST filter was reimplemented in Python rather than R, computing coverage
   identically as alignment length / query protein length.
3. Scaffold removal was applied after classification by requiring both members on
   a `Gm\d+` chromosome (numeric filter: 383 pairs dropped across categories).
