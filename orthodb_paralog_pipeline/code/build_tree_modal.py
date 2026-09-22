"""
Align an arbitrary set of protein sequences with MAFFT and build a maximum-
likelihood tree with IQ-TREE 2, running both steps as a Modal Function on
the user's own Modal account. Species-agnostic -- works for any FASTA.

Requires the `modal` Python SDK to be installed and authenticated
(MODAL_TOKEN_ID / MODAL_TOKEN_SECRET env vars, or `modal token new`).

Usage (CLI, needs `modal` importable in THIS interpreter):
    python build_tree_modal.py <input.fasta> <group_id> <output.json>

Usage (library, from pipeline.py -- dispatches via subprocess to a SEPARATE
interpreter that has `modal` installed/authenticated, so the orchestrator's
own environment never needs the `modal` SDK):
    from build_tree_modal import build_tree_via_subprocess
    result = build_tree_via_subprocess(fasta_path, group_id)
"""
import json
import os
import subprocess
import sys
import tempfile
import time


def build_tree_via_subprocess(fasta_path: str, group_id: str, modal_python: str = None) -> dict:
    modal_python = modal_python or os.environ.get("MODAL_PYTHON_BIN", "python")
    this_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.abspath(__file__)

    with tempfile.TemporaryDirectory() as tmp:
        out_path = os.path.join(tmp, "tree_result.json")
        proc = subprocess.run([modal_python, script_path, fasta_path, group_id, out_path],
                               capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0 or not os.path.exists(out_path):
            raise RuntimeError(f"Modal tree-build dispatch failed (returncode={proc.returncode}):\n"
                                f"STDOUT: {proc.stdout[-2000:]}\nSTDERR: {proc.stderr[-2000:]}")
        with open(out_path) as f:
            return json.load(f)


def _build_tree_modal_sdk(fasta_path: str, group_id: str) -> dict:
    """Actual Modal SDK call -- only imported when this script runs as the
    CLI entry point in an interpreter that has `modal` installed."""
    import modal

    app = modal.App("orthodb-paralog-pipeline-tree")

    # CPU-only image: Debian slim + MAFFT + IQ-TREE 2 from the Debian package
    # archive. No GPU needed for gene-family-scale alignment/tree building.
    image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install("mafft", "iqtree", "ca-certificates")
    )

    @app.function(image=image, cpu=4, memory=4096, timeout=1800)
    def run_pipeline(fasta_text: str, group_id: str) -> dict:
        import os
        import subprocess
        import time as _time

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

        aln_fa = os.path.join(workdir, "aligned.fasta")
        mafft_log = os.path.join(workdir, "mafft.log")
        t0 = _time.time()
        with open(aln_fa, "w") as out, open(mafft_log, "w") as errlog:
            p = subprocess.run([mafft_bin, "--auto", "--thread", "4", in_fa],
                                stdout=out, stderr=errlog)
        result["mafft_returncode"] = p.returncode
        result["mafft_wall_s"] = _time.time() - t0
        with open(mafft_log) as f:
            result["mafft_log"] = f.read()[-4000:]
        if p.returncode != 0:
            return result

        pre = os.path.join(workdir, "tree")
        t1 = _time.time()
        cmd = [iqtree_bin, "-s", aln_fa, "-m", "MFP", "-bb", "1000", "-alrt", "1000",
               "-nt", "AUTO", "-ntmax", "4", "-pre", pre, "-redo"]
        p2 = subprocess.run(cmd, capture_output=True, text=True)
        result["iqtree_returncode"] = p2.returncode
        result["iqtree_wall_s"] = _time.time() - t1
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

    with open(fasta_path) as f:
        fasta_text = f.read()
    with app.run():
        t0 = time.time()
        res = run_pipeline.remote(fasta_text, group_id)
        print(f"[driver] remote call finished in {time.time() - t0:.1f}s")
    return res


def main():
    fasta_path, group_id, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    res = _build_tree_modal_sdk(fasta_path, group_id)
    with open(out_path, "w") as f:
        json.dump(res, f)
    print("mafft_returncode:", res.get("mafft_returncode"))
    print("iqtree_returncode:", res.get("iqtree_returncode"))
    print("Best-fit model line:", [l for l in (res.get("iqtree_report") or "").splitlines()
                                    if "Best-fit model" in l])


if __name__ == "__main__":
    main()
