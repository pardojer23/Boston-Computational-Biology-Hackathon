#!/usr/bin/env Rscript
# Step 3 substep 3b: drop gene pairs that are unexpressed in the tissue.
#
# Port of filtered_gene_pairs.Rmd: average each of the four replicate columns
# across cell types per pair, discard pairs where all four averages are zero,
# keep the full unaveraged rows for the survivors, and split Gene_Pair into
# Gene1/Gene2. The three-column unique table is what
# get_paralog_coexpression.r consumes (it reads df$Gene1 and df$Gene2).
#
# Usage:  Rscript step3_03_filter_expressed_pairs.R <step3_out_dir>

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else "step3_out"

tissues <- c("root", "early_nodule", "hypocotyl", "globular", "early_maturation", "heart")
rows <- list()

for (nm in tissues) {
  d <- read.csv(file.path(OUT, sprintf("%s_filter_pairs.csv", nm)))
  avg <- d %>% group_by(Gene_Pair) %>%
    summarise(Gene1_Replicate_1 = mean(Gene1_Replicate_1),
              Gene1_Replicate_2 = mean(Gene1_Replicate_2),
              Gene2_Replicate_1 = mean(Gene2_Replicate_1),
              Gene2_Replicate_2 = mean(Gene2_Replicate_2), .groups = "drop") %>%
    filter(!(Gene1_Replicate_1 == 0 & Gene1_Replicate_2 == 0 &
             Gene2_Replicate_1 == 0 & Gene2_Replicate_2 == 0))
  keep <- unique(avg$Gene_Pair)
  expr <- d[d$Gene_Pair %in% keep, ]
  sep <- expr %>% separate(Gene_Pair, into = c("Gene1", "Gene2"),
                           sep = " vs ", remove = FALSE)
  write.csv(sep, file.path(OUT, sprintf("%s_expression_separated_raw.csv", nm)),
            row.names = FALSE)
  uniq <- unique(sep[, 1:3])
  write.csv(uniq, file.path(OUT, sprintf("%s_pairs_for_coexp.csv", nm)),
            row.names = FALSE)

  rows[[nm]] <- data.frame(tissue = nm, pairs_in = length(unique(d$Gene_Pair)),
                           pairs_expressed = length(keep),
                           pct = round(100 * length(keep) / length(unique(d$Gene_Pair)), 1),
                           rows_raw = nrow(sep))
  cat(sprintf("%-17s pairs %6d -> %6d expressed (%.1f%%)\n",
              nm, length(unique(d$Gene_Pair)), length(keep),
              100 * length(keep) / length(unique(d$Gene_Pair))))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "expressed_pairs_summary.csv"), row.names = FALSE)
print(res, row.names = FALSE)
