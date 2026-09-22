"""Step 3 substeps 5-6 of Li, Zhang & Schmitz 2025 (Plant Cell, koaf279):
build the gene-gene co-expression network per tissue, then score paralog-pair
co-expression specificity.

Runs the authors' get_coexp_network_no_filter.r and get_paralog_coexpression.r
logic verbatim. Both stay inside one container per tissue because the
intermediate network is a genes x genes matrix (~8 GB at 31k genes) that must
never cross the wire — only the small per-pair CSV is returned.

Memory is the binding constraint, which is why the authors ran this on a
cluster. Run one tissue first to measure actual peak RSS, then size the rest.

  modal run dupgen_step3_coexp_modal.py --tissue root
  modal run dupgen_step3_coexp_modal.py --tissue heart --memory-gb 128
"""

import modal

app = modal.App("soybean-dupgen-step3-coexp")

image = (
    modal.Image.from_registry("rocker/r-ver:4.4.2", add_python="3.11")
    .apt_install("time")
)

# Verbatim from get_coexp_network_no_filter.r / get_paralog_coexpression.r,
# with only the file plumbing changed (input/output paths instead of argv).
R_SCRIPT = r"""
build_coexp_network <- function(net, method = "spearman", flag = "rank") {
  genes = rownames(net)
  net = net[rowSums(net) > 0, ]
  net = cor(t(net), method = method)
  temp = net[upper.tri(net, diag = TRUE)]
  if (flag == "abs") { temp = abs(temp) }
  temp = rank(temp, ties.method = "average")
  net[upper.tri(net, diag = TRUE)] = temp
  net = t(net)
  net[upper.tri(net, diag = TRUE)] = temp
  net = net / max(net, na.rm = TRUE)
  med = median(net, na.rm = TRUE)
  ind = setdiff(genes, rownames(net))
  temp = matrix(med, length(ind), dim(net)[2]); rownames(temp) = ind
  net = rbind(net, temp)
  temp = matrix(med, dim(net)[1], length(ind)); colnames(temp) = ind
  net = cbind(net, temp)
  net = net[genes, genes]
  diag(net) = 1
  return(net)
}

calc_spec <- function (res, np, nL){
  ranks = 1:nL
  mini = sum(ranks[1:np])
  range = np * (nL - np)
  temp1 = t(apply(res, 1, function(x) rank(x, ties.method = "average")))
  spec1 <- (temp1 - mini) / range
  temp2 = t(apply(t(res), 1, function(x) rank(x, ties.method = "average")))
  spec2 <- (temp2 - mini) / range
  spec = 0.5 * (spec1 + t(spec2))
  return(spec)
}

exp1 = read.csv("/work/normalized.csv", row.names = 1)
exp1 <- as.matrix(exp1)
cat("expression matrix:", dim(exp1), "\n")

t0 <- Sys.time()
# NOTE: the authors' script passes method='pearson' at the call site even though
# the function default is 'spearman'. Pearson is what they ran, so Pearson it is.
prionet_atlas = build_coexp_network(exp1, method = 'pearson')
cat("network built:", dim(prionet_atlas), "| mins",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "\n")
cat("peak RSS so far (GB):", round(as.numeric(gc()[2, 6]) / 1024, 2), "\n")

df <- read.csv("/work/pairs.csv")
diag(prionet_atlas) <- 0
t1 <- Sys.time()
corrspec <- calc_spec(prionet_atlas, 1, dim(prionet_atlas)[1])
cat("specificity done | mins",
    round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 2), "\n")

allgenes <- rownames(corrspec)
df$coexpression <- NA
for (id in 1:nrow(df)) {
  g1 <- match(df$Gene1[id], allgenes)
  g2 <- match(df$Gene2[id], allgenes)
  if (!is.na(g1 + g2)) { df$coexpression[id] <- corrspec[g1, g2] }
}
cat("pairs scored:", sum(!is.na(df$coexpression)), "of", nrow(df), "\n")
cat("peak memory used (GB):", round(sum(gc()[, 6]) / 1024, 2), "\n")
write.table(df, "/work/paralog_pair_coexpression.csv", sep = ',',
            row.names = FALSE, col.names = TRUE, quote = FALSE)
"""


@app.function(image=image, cpu=8.0, memory=65536, timeout=14400)
def coexp(tissue: str, normalized_csv: bytes, pairs_csv: bytes):
    import os
    import resource
    import subprocess

    os.makedirs("/work", exist_ok=True)
    open("/work/normalized.csv", "wb").write(normalized_csv)
    open("/work/pairs.csv", "wb").write(pairs_csv)
    open("/work/run.R", "w").write(R_SCRIPT)

    r = subprocess.run(["Rscript", "/work/run.R"], capture_output=True, text=True)
    peak_gb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / (1024 ** 2)
    log = (f"tissue: {tissue}\npeak child RSS (GB): {peak_gb:.2f}\n"
           f"exit: {r.returncode}\n--- stdout ---\n{r.stdout}\n"
           f"--- stderr (tail) ---\n{r.stderr[-3000:]}\n")
    out = {}
    p = "/work/paralog_pair_coexpression.csv"
    if os.path.exists(p):
        out[f"{tissue}_filter_paralog_pair_coexpression.csv"] = open(p, "rb").read()
    out[f"{tissue}_coexp.log"] = log.encode()
    return out


@app.local_entrypoint()
def main(tissue: str = "root", memory_gb: int = 64):
    import os

    norm = open(f"step3_out/{tissue}_filter_normalized.csv", "rb").read()
    pairs = open(f"step3_out/{tissue}_pairs_for_coexp.csv", "rb").read()
    print(f"{tissue}: normalized {len(norm)/1e6:.1f} MB, pairs {len(pairs)/1e6:.1f} MB")
    out = coexp.with_options(memory=memory_gb * 1024).remote(tissue, norm, pairs)
    base = "modal_out_step3"
    os.makedirs(base, exist_ok=True)
    for name, data in out.items():
        with open(os.path.join(base, name), "wb") as fh:
            fh.write(data)
        print(f"{len(data):>12,}  {base}/{name}")
