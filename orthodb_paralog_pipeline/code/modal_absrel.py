"""
Modal Function: codon-aware PRANK alignment + HyPhy aBSREL branch-site
selection test, for an arbitrary set of CDS sequences + guide tree (any
number of taxa). Generalized from the leghemoglobin-specific version used
earlier in this project (see leghemoglobin_absrel_selection/).

Usage (CLI):
    python modal_absrel.py <cds.fasta> <tree.nwk> <output.json>

Usage (library, imported by soybean_selection.py):
    from modal_absrel import run_pipeline_local
    result = run_pipeline_local(cds_fasta_text, tree_newick_text)
"""
import json
import sys
import time

import modal

app = modal.App("orthodb-paralog-pipeline-absrel")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("wget", "ca-certificates")
    .run_commands(
        "wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh",
        "bash /tmp/miniconda.sh -b -p /opt/conda",
        "/opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main",
        "/opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r",
        "/opt/conda/bin/conda install -y -n base --override-channels -c bioconda -c conda-forge prank hyphy",
        "/opt/conda/bin/conda clean -afy",
        "/opt/conda/bin/prank -help 2>&1 | head -5",
        "/opt/conda/bin/hyphy --help 2>&1 | head -5",
    )
)


@app.function(image=image, cpu=4, memory=4096, timeout=1800)
def run_pipeline(cds_fasta: str, tree_newick: str) -> dict:
    import glob
    import os
    import subprocess

    workdir = "/tmp/work"
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    with open("cds.fasta", "w") as f:
        f.write(cds_fasta)
    with open("guide_tree.nwk", "w") as f:
        f.write(tree_newick)

    result = {}

    prank_bin = "/opt/conda/bin/prank"
    cmd = [prank_bin, "-d=cds.fasta", "-t=guide_tree.nwk", "-o=aligned", "-codon", "-F"]
    p1 = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    result["prank_returncode"] = p1.returncode
    result["prank_stdout"] = p1.stdout[-3000:]
    result["prank_stderr"] = p1.stderr[-3000:]

    aligned_files = glob.glob("aligned*.fas") + glob.glob("aligned*.best.fas")
    result["aligned_files_found"] = aligned_files
    if not aligned_files:
        result["dir_listing_after_prank"] = os.listdir(".")
        return result
    aligned_path = sorted(aligned_files, key=len)[0]
    with open(aligned_path) as f:
        result["aligned_fasta"] = f.read()

    hyphy_bin = "/opt/conda/bin/hyphy"
    cmd2 = [hyphy_bin, "absrel", "--alignment", aligned_path, "--tree", "guide_tree.nwk",
            "--output", "absrel_result.json"]
    p2 = subprocess.run(cmd2, capture_output=True, text=True, timeout=1200)
    result["hyphy_returncode"] = p2.returncode
    result["hyphy_stdout"] = p2.stdout[-6000:]
    result["hyphy_stderr"] = p2.stderr[-3000:]

    if os.path.exists("absrel_result.json"):
        with open("absrel_result.json") as f:
            result["absrel_result"] = f.read()
    else:
        result["absrel_result"] = None
        result["dir_listing"] = os.listdir(".")

    return result


def run_pipeline_local(cds_fasta_text: str, tree_newick_text: str) -> dict:
    with app.run():
        t0 = time.time()
        res = run_pipeline.remote(cds_fasta_text, tree_newick_text)
        print(f"[driver] remote call finished in {time.time() - t0:.1f}s")
    return res


def main():
    cds_path, tree_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(cds_path) as f:
        cds_fasta = f.read()
    with open(tree_path) as f:
        tree_newick = f.read()

    res = run_pipeline_local(cds_fasta, tree_newick)

    with open(out_path, "w") as f:
        json.dump(res, f)

    print("prank_returncode:", res.get("prank_returncode"))
    print("aligned_files_found:", res.get("aligned_files_found"))
    print("hyphy_returncode:", res.get("hyphy_returncode"))


if __name__ == "__main__":
    main()
