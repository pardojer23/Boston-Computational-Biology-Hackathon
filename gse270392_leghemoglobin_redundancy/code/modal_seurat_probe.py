"""
Heavy-compute step, run on Modal: verify whether the two leghemoglobin
paralogs missing from GSE270392's small released pseudobulk matrices
(Glyma.10G199000 = Lbc1, Glyma.10G199100 = Lba) are genuinely undetected in
the underlying single-nucleus RNA-seq data, or were merely excluded from the
curated summary file.

This downloads and loads the full per-tissue Seurat object
("*_Gm_atlas_<tissue>.rna.seurat.obj.rds.gz", several hundred MB each,
Seurat objects at this scale need several GB of RAM once deserialized) on a
Modal container built from the satijalab/seurat Docker image, and checks gene
presence in both the SCTransform-filtered assay and the raw/unfiltered "RNA"
assay (the largest gene universe released for that tissue).

Usage:
    python modal_seurat_probe.py <tissue_rds_gz_url> <output.json>
"""
import json
import sys

import modal

app = modal.App("gse270392-leghemoglobin-probe")

image = modal.Image.from_registry("satijalab/seurat:5.0.0", add_python="3.11")

# Target Glyma IDs, without the "ann1." prefix the released Seurat objects
# attach to feature names (rownames are "ann1.Glyma.XXGXXXXXX").
TARGET_GENES = [
    "Glyma.10G198800",  # Lbc3 (LB1)
    "Glyma.10G199000",  # Lbc1 (LB2)
    "Glyma.10G199100",  # Lba  (LB3)
    "Glyma.20G191200",  # Lbc2 (LB4)
    "Glyma.10G198900",  # unnamed paralog near the cluster
]


@app.function(image=image, cpu=4, memory=16384, timeout=1800)
def probe(url: str, target_genes: list) -> dict:
    import subprocess
    import os

    r_script = r"""
    args <- commandArgs(trailingOnly=TRUE)
    url <- args[1]; targets_file <- args[2]; out_file <- args[3]
    target_genes <- scan(targets_file, what=character(), sep="
", quiet=TRUE)

    download.file(url, "obj.rds.gz", mode="wb", quiet=TRUE)
    system("gunzip -f obj.rds.gz")
    suppressMessages(library(Seurat))
    obj <- readRDS("obj.rds")

    result <- list(ncells = ncol(obj), assays = names(obj@assays))
    for (assay_name in names(obj@assays)) {
      counts <- tryCatch(GetAssayData(obj, assay=assay_name, layer="counts"),
                          error=function(e) NULL)
      if (is.null(counts)) next
      rn <- rownames(counts)
      # feature names carry an "ann1." annotation-source prefix
      rn_stripped <- sub("^ann1\\.", "", rn)
      present <- target_genes[target_genes %in% rn_stripped]
      result[[paste0(assay_name, "_dim")]] <- dim(counts)
      result[[paste0(assay_name, "_present")]] <- present
      result[[paste0(assay_name, "_absent")]] <- setdiff(target_genes, present)
      if (length(present) > 0) {
        idx <- match(present, rn_stripped)
        sub_counts <- counts[idx, , drop=FALSE]
        result[[paste0(assay_name, "_total_counts")]] <- as.list(Matrix::rowSums(sub_counts))
        result[[paste0(assay_name, "_ncells_detected")]] <- as.list(Matrix::rowSums(sub_counts > 0))
      }
    }
    writeLines(jsonlite::toJSON(result, auto_unbox=TRUE, na="null"), out_file)
    """
    with open("probe.R", "w") as f:
        f.write(r_script)
    with open("targets.txt", "w") as f:
        f.write("\n".join(target_genes))

    p = subprocess.run(["Rscript", "probe.R", url, "targets.txt", "out.json"],
                        capture_output=True, text=True, timeout=1700)
    out = {"returncode": p.returncode, "stderr_tail": p.stderr[-3000:]}
    if os.path.exists("out.json"):
        with open("out.json") as f:
            out["result"] = json.load(f)
    return out


def main() -> None:
    url = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "seurat_probe_result.json"
    with app.run():
        res = probe.remote(url, TARGET_GENES)
    with open(out_path, "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res.get("result", res), indent=2)[:3000])


if __name__ == "__main__":
    main()
