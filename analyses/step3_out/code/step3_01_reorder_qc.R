#!/usr/bin/env Rscript
# Step 3 substep 1 of Li, Zhang & Schmitz 2025 (Plant Cell, koaf279):
# reorder each tissue's replicate-split count matrix so replicate pairs are
# adjacent, compute the per-cell-type Spearman correlation between replicates,
# and apply the retention filter.
#
# Port of reorder_gene_counts.Rmd. The column-index vectors and the post-hoc
# "eff" column selections are transcribed verbatim from that notebook.
#
# Cotyledon is absent: the shipped
# Gm_atlas_Cotyledon_stage_seeds_RNA_gene_celltype_splitReps_counts.txt is 2
# bytes (empty) and GEO ships no replicate-split equivalent, so the tissue
# cannot be processed. Six of seven tissues are handled here.
#
# Two documented findings this script reproduces:
#  1. The stated rho > 0.75 rule matches the authors' manual selection for only
#     3 of 6 tissues. For early_maturation, globular and hypocotyl they also
#     dropped cell types whose rho was ABOVE the threshold; those are the
#     "Unknown"/"SC_unknown" unannotated clusters.
#  2. The reorder vectors for globular and hypocotyl (and cotyledon) repeat
#     index 23 where 13 was intended, so their final "replicate pair" compares
#     two different cell types' replicate-2 columns and one cell type's
#     replicate-1 column is never used. The authors' own eff selection discards
#     that pair, so it does not reach their downstream data.
#
# Usage:  Rscript step3_01_reorder_qc.R <raw_data_dir> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
RAW <- if (length(args) >= 1) args[1] else "repo/step3_divergence_expression_pattern/raw_data"
OUT <- if (length(args) >= 2) args[2] else "step3_out"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f <- function(t) file.path(RAW, sprintf("Gm_atlas_%s_RNA_gene_celltype_splitReps_counts.txt", t))

spec <- list(
  early_maturation = list(file = "Early_maturation_stage_seeds",
    ord = c(1,12,2,13,3,14,4,15,5,16,6,17,7,18,8,19,9,20,10,21,11,22), eff = 1:20),
  early_nodule = list(file = "Early_nodule",
    ord = c(1,10,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9,18), eff = c(1:2,5:12,15:16)),
  globular = list(file = "Globular_stage_seeds",
    ord = c(1,14,2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,23,26),
    eff = c(1:20,23,24)),
  heart = list(file = "Heart_stage_seeds",
    ord = c(1,13,2,14,3,15,4,16,5,17,6,18,7,19,8,20,9,21,10,22,11,23,12,24), eff = NULL),
  hypocotyl = list(file = "Hypocotyl",
    ord = c(1,14,2,15,3,16,4,17,5,18,6,19,7,20,8,21,9,22,10,23,11,24,12,25,23,26),
    eff = c(1:20,23,24)),
  root = list(file = "Root",
    ord = c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16), eff = c(1:8,11:14))
)

spearman_pairs <- function(d) {
  np <- ncol(d) %/% 2
  sapply(seq_len(np), function(i) {
    a <- as.numeric(d[[2*i-1]]); b <- as.numeric(d[[2*i]])
    if (length(unique(a)) == 1 || length(unique(b)) == 1) NA_real_
    else cor(a, b, method = "spearman")
  })
}

summary_rows <- list()
for (nm in names(spec)) {
  s <- spec[[nm]]
  raw <- read.delim2(f(s$file))
  ord <- raw[, s$ord]
  eff <- if (is.null(s$eff)) ord else ord[, s$eff]
  rho <- spearman_pairs(ord)
  np  <- ncol(ord) %/% 2
  keep_auto <- which(!is.na(rho) & rho > 0.75)
  authors_cols <- if (is.null(s$eff)) seq_len(ncol(ord)) else s$eff
  authors_pairs <- which(vapply(seq_len(np),
      function(i) all(c(2*i-1, 2*i) %in% authors_cols), logical(1)))
  agree <- identical(as.integer(keep_auto), as.integer(authors_pairs))
  dupi <- anyDuplicated(s$ord)
  cat(sprintf("\n=== %s ===  raw cols %d | reordered %d | pairs %d\n",
              nm, ncol(raw), ncol(ord), np))
  cat("  rho:", paste(sprintf("%.3f", rho), collapse = " "), "\n")
  cat("  pairs rho>0.75    :", paste(keep_auto, collapse = ","), "\n")
  cat("  pairs authors kept:", paste(authors_pairs, collapse = ","), "\n")
  cat("  agree:", agree, "| duplicated col index in reorder vector:", dupi > 0, "\n")

  # authors-faithful selection, used downstream
  write.csv(eff, file.path(OUT, sprintf("%s_reorder.csv", nm)))
  # rule-based selection, for comparison only
  auto <- ord[, as.vector(rbind(2*keep_auto-1, 2*keep_auto))]
  write.csv(auto, file.path(OUT, sprintf("%s_reorder_auto.csv", nm)))

  summary_rows[[nm]] <- data.frame(tissue = nm, raw_cols = ncol(raw), pairs = np,
    pairs_gt_075 = length(keep_auto), pairs_authors = length(authors_pairs),
    agree = agree, dup_index_bug = dupi > 0)
}
res <- do.call(rbind, summary_rows)
write.csv(res, file.path(OUT, "reorder_qc_summary.csv"), row.names = FALSE)
print(res, row.names = FALSE)
