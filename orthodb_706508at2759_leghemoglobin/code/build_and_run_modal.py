"""
Align an OrthoDB orthologous-group FASTA with MAFFT and build a maximum-
likelihood tree with IQ-TREE 2, running both steps as a Modal Function on
the user's own Modal account.

Requires the `modal` Python SDK to be installed and authenticated
(MODAL_TOKEN_ID / MODAL_TOKEN_SECRET env vars, or `modal token new`).

Usage:
    python build_and_run_modal.py <input.fasta> <group_id>

Writes `modal_result.json` containing the alignment, the ML tree
(.treefile / .contree Newick strings), and the full IQ-TREE report.
"""
import json
import sys
import time

import modal

app = modal.App("phylo-iqtree-706508at2759")

# CPU-only image: Debian slim + MAFFT + IQ-TREE 2 from the Debian package
# archive. No GPU needed for a single-gene-family alignment/tree this size.
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("mafft", "iqtree", "ca-certificates")
)


@app.function(image=image, cpu=4, memory=4096, timeout=1800)
def run_pipeline(fasta_text: str, group_id: str) -> dict:
    import os
    import subprocess

    workdir = "/tmp/work"
    os.makedirs(workdir, exist_ok=True)
    in_fa = os.path.join(workdir, "input.fasta")
    with open(in_fa, "w") as f:
        f.write(fasta_text)

    def which(name):
        r = subprocess.run(["which", name], capture_output=True, text=True)
        return r.stdout.strip() or None

    mafft_bin = which("mafft")
    iqtree_bin = which("iqtree2") or which("iqtree")
    result = {"mafft_bin": mafft_bin, "iqtree_bin": iqtree_bin}

    # --- Alignment ---
    aln_fa = os.path.join(workdir, "aligned.fasta")
    mafft_log = os.path.join(workdir, "mafft.log")
    t0 = time.time()
    with open(aln_fa, "w") as out, open(mafft_log, "w") as errlog:
        p = subprocess.run([mafft_bin, "--auto", "--thread", "4", in_fa],
                            stdout=out, stderr=errlog)
    result["mafft_returncode"] = p.returncode
    result["mafft_wall_s"] = time.time() - t0
    with open(mafft_log) as f:
        result["mafft_log"] = f.read()[-4000:]
    if p.returncode != 0:
        return result

    # --- Tree inference: ModelFinder + ML tree + UFBoot + SH-aLRT ---
    pre = os.path.join(workdir, "tree")
    t1 = time.time()
    cmd = [iqtree_bin, "-s", aln_fa, "-m", "MFP", "-bb", "1000", "-alrt", "1000",
           "-nt", "AUTO", "-ntmax", "4", "-pre", pre, "-redo"]
    p2 = subprocess.run(cmd, capture_output=True, text=True)
    result["iqtree_returncode"] = p2.returncode
    result["iqtree_wall_s"] = time.time() - t1
    result["iqtree_stdout"] = p2.stdout[-6000:]
    result["iqtree_stderr"] = p2.stderr[-3000:]

    def read_if_exists(path):
        return open(path).read() if os.path.exists(path) else None

    result["aligned_fasta"] = open(aln_fa).read()
    result["treefile"] = read_if_exists(pre + ".treefile")
    result["contree"] = read_if_exists(pre + ".contree")
    result["iqtree_report"] = read_if_exists(pre + ".iqtree")
    result["iqtree_full_log"] = read_if_exists(pre + ".log")
    return result


def main() -> None:
    fasta_path, group_id = sys.argv[1], sys.argv[2]
    with open(fasta_path) as f:
        fasta_text = f.read()

    with app.run():
        t0 = time.time()
        res = run_pipeline.remote(fasta_text, group_id)
        print(f"[driver] remote call finished in {time.time() - t0:.1f}s")

    with open("modal_result.json", "w") as f:
        json.dump(res, f)

    print("mafft_bin:", res.get("mafft_bin"))
    print("iqtree_bin:", res.get("iqtree_bin"))
    print("mafft_returncode:", res.get("mafft_returncode"))
    print("iqtree_returncode:", res.get("iqtree_returncode"))


if __name__ == "__main__":
    main()
