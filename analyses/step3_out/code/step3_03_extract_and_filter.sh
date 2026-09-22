#!/bin/bash
# Step 3 substep 3 of Li, Zhang & Schmitz 2025 (Plant Cell, koaf279):
# extract each duplicated pair's expression per cell type, then drop pairs that
# are unexpressed in the tissue.
#
# Runs the authors' extract_duplicates_expression.py VERBATIM, then applies the
# filter from filtered_gene_pairs.Rmd (step3_03_filter_expressed_pairs.R).
#
# IMPORTANT: extract_duplicates_expression.py reads the gene-pairs file with
# pandas `header=None` and `usecols=[0, 2]`. The authors' shipped
# 07042025_glyma_duplicated_pairs.csv has no header row; a reproduced pairs file
# that HAS one must have it stripped, or the header is silently consumed as a
# gene pair. That is what the tail -n +2 below is for.
#
# Requires pandas. Usage:
#   bash step3_03_extract_and_filter.sh <pairs_csv> <repo_dir> <step3_out_dir>

set -uo pipefail
PAIRS="${1:-reproduced_pairs_lexicographic.csv}"
REPO="${2:-repo}"
OUT="${3:-step3_out}"
SCRIPT="$REPO/step3_divergence_expression_pattern/extract_duplicates_expression.py"

if ! python3 -c "import pandas" 2>/dev/null; then
  echo "ERROR: pandas not importable by this python3 ($(command -v python3))." >&2
  echo "The authors' extract script requires it." >&2
  exit 1
fi
[ -f "$SCRIPT" ] || { echo "ERROR: not found: $SCRIPT" >&2; exit 1; }

# strip the header if present, so the file matches what the script expects
HEADERLESS="$OUT/pairs_headerless.csv"
if head -1 "$PAIRS" | grep -qi 'duplicate_1\|Gene1'; then
  tail -n +2 "$PAIRS" > "$HEADERLESS"
  echo "stripped header from $PAIRS -> $HEADERLESS"
else
  cp "$PAIRS" "$HEADERLESS"
  echo "$PAIRS already headerless"
fi
echo "pairs: $(wc -l < "$HEADERLESS")"

for t in root early_nodule hypocotyl globular early_maturation heart; do
  start=$(date +%s)
  python3 "$SCRIPT" "$OUT/${t}_filter_normalized.csv" "$HEADERLESS" "$OUT/${t}_filter_pairs.csv" >/dev/null 2>&1
  rc=$?
  if [ -f "$OUT/${t}_filter_pairs.csv" ]; then
    rows=$(( $(wc -l < "$OUT/${t}_filter_pairs.csv") - 1 ))
  else
    rows="FAILED"
  fi
  printf "%-17s rc=%s %9s rows %4ds\n" "$t" "$rc" "$rows" "$(( $(date +%s) - start ))"
done

echo
echo "--- applying the expressed-pair filter ---"
Rscript "$(dirname "$0")/step3_03_filter_expressed_pairs.R" "$OUT"
