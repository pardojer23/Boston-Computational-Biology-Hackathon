#!/usr/bin/env Rscript
# Step 3 substep 2 of Li, Zhang & Schmitz 2025 (Plant Cell, koaf279):
# edgeR filtering and CPM normalization per tissue.
#
# Port of normalization_edgeR.Rmd: DGEList grouped by cell type (each cell type
# appearing twice, once per replicate), filterByExpr, normLibSizes, cpm.
#
# Group LABELS are derived from the real column headers rather than the
# notebook's hand-typed lists, which contain typos ("Prpcambium",
# "Em_proper_innitials"). Only the grouping STRUCTURE affects filterByExpr, and
# the retained cell-type counts produced here (10, 6, 11, 12, 11, 6) match the
# notebook's hand-typed group vectors exactly for all six tissues.
#
# Usage:  Rscript step3_02_normalize_edger.R <step3_out_dir>

suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else "step3_out"

cat("edgeR", as.character(packageVersion("edgeR")), "\n")
tissues <- c("early_maturation", "early_nodule", "globular", "heart", "hypocotyl", "root")
rows <- list()

for (nm in tissues) {
  d <- read.csv(file.path(OUT, sprintf("%s_reorder.csv", nm)))
  rownames(d) <- d$X
  d <- d[, -1, drop = FALSE]
  d[] <- lapply(d, function(x) as.numeric(as.character(x)))
  stopifnot(ncol(d) %% 2 == 0)

  grp <- rep(seq_len(ncol(d) / 2), each = 2)
  y <- DGEList(counts = as.matrix(d), group = factor(grp))
  keep <- filterByExpr(y, group = y$samples$group)
  y <- y[keep, , keep = FALSE]
  y <- normLibSizes(y)
  norm <- cpm(y)
  write.csv(norm, file.path(OUT, sprintf("%s_filter_normalized.csv", nm)))

  rows[[nm]] <- data.frame(tissue = nm, cell_types = ncol(d) / 2, columns = ncol(d),
    genes_in = length(keep), genes_kept = sum(keep),
    pct_kept = round(100 * sum(keep) / length(keep), 1),
    norm_factor_range = sprintf("%.3f-%.3f",
      min(y$samples$norm.factors), max(y$samples$norm.factors)))
  cat(sprintf("%-17s cell types %2d | genes %5d -> %5d (%.1f%%)\n",
              nm, ncol(d) / 2, length(keep), sum(keep),
              100 * sum(keep) / length(keep)))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "normalization_summary.csv"), row.names = FALSE)
print(res, row.names = FALSE)
