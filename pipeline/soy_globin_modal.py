"""
Modal app: soybean (Glycine max, Wm82.a4.v1) globin-family sequence module.

    modal run soy_globin_modal.py                 # full pipeline -> ./results
    modal run soy_globin_modal.py --outdir out2   # alternate output directory
    modal run soy_globin_modal.py --no-gpu        # skip the ESM2 stage
    modal run soy_globin_modal.py --force-fetch   # re-download reference files

What it does
------------
1.  [cpu ]  download the Wm82.gnm4.ann1.T8TQ primary-transcript proteome and
            gene-model GFF3 from SoyBase, plus the Pfam PF00042 (Globin) HMM
            from InterPro, into a persistent Modal Volume.
2.  [cpu ]  hmmsearch --cut_ga to pull every globin-domain protein, unioned
            with the four focal leghemoglobins (Lba, Lbc1, Lbc2, Lbc3).
3.  [cpu ]  MAFFT L-INS-i alignment -> IQ-TREE 2 (ModelFinder, 1000 UFBoot,
            1000 SH-aLRT) -> Newick; pairwise percent identity from the MSA.
4.  [gpu ]  ESM2-650M mean-pooled embeddings -> cosine-distance matrix.
            Runs concurrently with stage 3.
5.  [cpu ]  GFF3 gene order -> every within-family pair labelled
            tandem / proximal / dispersed.

Stage 4 is the only Modal-worthy step by compute (a 650M-parameter forward
pass); stages 1-3 and 5 are on Modal for data locality — the 27 MB of
reference downloads and the hmmsearch over 52,872 proteins stay next to the
Volume instead of crossing the wire twice.

All pipeline logic lives in ``soy_globin_core.py``; this file is only the
Modal wiring, so the same code can be run locally without Modal.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
import time

# So `modal run pipeline/soy_globin_modal.py` works from the repo root:
# add_local_python_source() imports soy_globin_core in the *local* interpreter
# at build time, and Modal does not put the entrypoint's directory on sys.path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import modal

APP_NAME = "soy-globin-redundancy"
DATA_DIR = "/data"          # persistent reference + intermediates
CACHE_DIR = "/cache"        # HuggingFace weights

app = modal.App(APP_NAME)

data_vol = modal.Volume.from_name("soy-globin-data", create_if_missing=True)
cache_vol = modal.Volume.from_name("soy-globin-hf-cache", create_if_missing=True)

# --- images ---------------------------------------------------------------- #
# Bioconda toolchain; pinned so a rerun six months from now gives the same tree.
cpu_image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "hmmer=3.4",
        "mafft=7.526",
        "iqtree=3.1.3",  # same pin as environment.yml, so Modal and local agree
        "biopython",
        "pandas",
        "numpy",
        channels=["conda-forge", "bioconda"],
    )
    .add_local_python_source("soy_globin_core", "config", "contracts")
)

gpu_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install(
        "torch==2.4.1",
        "transformers==4.44.2",
        "numpy<2",
        "pandas",
    )
    .env({"HF_HOME": f"{CACHE_DIR}/huggingface"})
    .add_local_python_source("soy_globin_core", "config", "contracts")
)


# --------------------------------------------------------------------------- #
# Stage 1+2 — fetch references, hmmsearch, family membership
# --------------------------------------------------------------------------- #
@app.function(
    image=cpu_image,
    volumes={DATA_DIR: data_vol},
    cpu=8.0,
    memory=8192,
    timeout=60 * 30,
)
def select_globin_family(force_fetch: bool = False) -> dict:
    """Download references (cached in the Volume) and select the globin family."""
    import config as cfgmod
    import soy_globin_core as core

    t0 = time.time()
    cfg = cfgmod.load()
    ref = cfg.section("reference")["files"]
    hmm = cfg.section("sequence")["hmm"]
    paths = {
        "proteome": core.fetch_file(ref["proteome"]["url"],
                                    f"{DATA_DIR}/ref/proteins_primary.faa.gz",
                                    ref["proteome"]["sha256"], force=force_fetch)["path"],
        "gff3": core.fetch_file(ref["gff3"]["url"],
                                f"{DATA_DIR}/ref/gene_models_exons.gff3.gz",
                                ref["gff3"]["sha256"], force=force_fetch)["path"],
    }
    paths["pfam_hmm"] = core.fetch_file(
        hmm["url"], f"{DATA_DIR}/ref/PF00042.hmm.gz", hmm["sha256"],
        force=force_fetch, gunzip_to=f"{DATA_DIR}/ref/PF00042.hmm")["gunzipped_path"]
    data_vol.commit()

    _seq = cfg.section("sequence")
    _dup = _seq["duplication"]
    _ident = cfg["family.identifiers"]
    seqs, members = core.select_family(
        paths["proteome"], paths["pfam_hmm"], f"{DATA_DIR}/work/hmmer", cfg, cpu=8
    )
    genes = core.parse_gff_genes(paths["gff3"], _ident["id_prefix_strip"],
                                 _ident["chromosome_regex"])

    return {
        "seqs": seqs,
        "members_csv": members.to_csv(index=False).encode(),
        "genes_parquet": genes.to_parquet(index=False),
        "fasta": _fasta_bytes(seqs),
        "meta": {
            "pfam": core.hmm_metadata(paths["pfam_hmm"]),
            "proteome_sha256": core.sha256(paths["proteome"]),
            "gff3_sha256": core.sha256(paths["gff3"]),
            "hmm_sha256": core.sha256(paths["pfam_hmm"]),
            "n_proteins_searched": len(core.read_fasta(paths["proteome"])),
            "n_family_members": len(members),
            "hmmsearch_version": core.tool_version(["hmmsearch", "-h"]),
            "stage_seconds": round(time.time() - t0, 1),
        },
    }


def _fasta_bytes(seqs: dict[str, str], width: int = 60) -> bytes:
    out = []
    for name, seq in seqs.items():
        out.append(f">{name}\n")
        out += [seq[i : i + width] + "\n" for i in range(0, len(seq), width)]
    return "".join(out).encode()


# --------------------------------------------------------------------------- #
# Stage 3 — alignment, phylogeny, pairwise identity
# --------------------------------------------------------------------------- #
@app.function(
    image=cpu_image,
    volumes={DATA_DIR: data_vol},
    cpu=8.0,
    memory=8192,
    timeout=60 * 60,
)
def align_tree_identity(seqs: dict[str, str]) -> dict:
    """MAFFT L-INS-i -> IQ-TREE 2 -> Newick + percent-identity matrices."""
    import config as cfgmod
    import soy_globin_core as core

    t0 = time.time()
    wd = pathlib.Path(f"{DATA_DIR}/work/phylo")
    wd.mkdir(parents=True, exist_ok=True)

    faa, aln = wd / "globins.faa", wd / "globins.aln.faa"
    core.write_fasta(seqs, faa)
    _aln_cfg = cfgmod.load().section("sequence")["alignment"]
    _phy = cfgmod.load().section("sequence")["phylogeny"]
    core.run_mafft(faa, aln, threads=8, args=_aln_cfg["args"])
    tree = core.run_iqtree(aln, wd / "globins", threads="8", seed=int(_phy["seed"]),
                           model=_phy["model"], ufboot=int(_phy["ufboot"]),
                           alrt=int(_phy["alrt"]))

    mat, long = core.pairwise_identity(aln)
    data_vol.commit()

    return {
        "alignment": aln.read_bytes(),
        "treefile": pathlib.Path(tree["treefile"]).read_bytes(),
        "iqtree_report": pathlib.Path(tree["iqtree_report"]).read_bytes(),
        "identity_matrix_csv": mat.to_csv().encode(),
        "identity_pairs_parquet": long.to_parquet(index=False),
        "meta": {
            "best_model": core.best_model(tree["iqtree_report"]),
            "n_alignment_columns": len(next(iter(core.read_fasta(aln).values()))),
            "mafft_version": core.tool_version(["mafft", "--version"]),
            "iqtree_version": core.tool_version([core.iqtree_exe(), "--version"]),
            "stage_seconds": round(time.time() - t0, 1),
        },
    }


# --------------------------------------------------------------------------- #
# Stage 4 — ESM2-650M embeddings (GPU)
# --------------------------------------------------------------------------- #
@app.function(
    image=gpu_image,
    volumes={CACHE_DIR: cache_vol},
    gpu="A10G",
    memory=16384,
    timeout=60 * 30,
)
def esm2_cosine(seqs: dict[str, str]) -> dict:
    """Mean-pooled ESM2-650M embeddings and their cosine-distance matrix."""
    import io

    import numpy as np
    import config as cfgmod
    import soy_globin_core as core
    import torch

    t0 = time.time()
    e = cfgmod.load().section("sequence")["embedding"]
    labels, emb = core.esm2_embeddings(seqs, model_name=e["model"],
                                       batch_size=int(e["batch_size"]),
                                       dtype_name=e["dtype"])
    cache_vol.commit()
    dist = core.cosine_distance_matrix(labels, emb)

    buf = io.BytesIO()
    np.savez_compressed(buf, labels=np.array(labels), embeddings=emb)

    return {
        "cosine_matrix_csv": dist.to_csv().encode(),
        "embeddings_npz": buf.getvalue(),
        "meta": {
            "model": e["model"],
            "embedding_dim": int(emb.shape[1]),
            "pooling": "mean over residue tokens (BOS/EOS/PAD excluded)",
            "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
            "torch": torch.__version__,
            "stage_seconds": round(time.time() - t0, 1),
        },
    }


# --------------------------------------------------------------------------- #
# Stage 5 — genomic context / duplication mode
# --------------------------------------------------------------------------- #
@app.function(image=cpu_image, cpu=2.0, memory=4096, timeout=60 * 10)
def paralog_pairs(
    members_csv: bytes,
    genes_parquet: bytes,
    identity_pairs_parquet: bytes,
    cosine_matrix_csv: bytes | None,
) -> dict:
    """Label every within-family pair tandem / proximal / dispersed from GFF3."""
    import io

    import pandas as pd
    import config as cfgmod
    import soy_globin_core as core

    members = pd.read_csv(io.BytesIO(members_csv))
    genes = pd.read_parquet(io.BytesIO(genes_parquet))
    ident = pd.read_parquet(io.BytesIO(identity_pairs_parquet))
    cos = (
        pd.read_csv(io.BytesIO(cosine_matrix_csv), index_col=0)
        if cosine_matrix_csv
        else None
    )

    _dup = cfgmod.load().section("sequence")["duplication"]
    pairs, ctx = core.classify_pairs(
        members, genes,
        max_intervening=int(_dup["tandem_max_intervening_genes"]),
        proximal_max_bp=int(_dup["proximal_max_bp"]),
        identity_long=ident, cosine=cos)
    return {
        "paralog_pairs_csv": pairs.to_csv(index=False).encode(),
        "gene_context_csv": ctx.to_csv(index=False).encode(),
        "context_parquet": ctx.to_parquet(index=False),
        "meta": {
            "tandem_max_intervening_genes": _dup["tandem_max_intervening_genes"],
            "proximal_max_bp": _dup["proximal_max_bp"],
            "mode_counts": pairs.duplication_mode.value_counts().to_dict(),
        },
    }


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
@app.local_entrypoint()
def main(outdir: str = "results", no_gpu: bool = False, force_fetch: bool = False):
    out = pathlib.Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    wall0 = time.time()

    print("[1/5] references + hmmsearch (PF00042) ...")
    fam = select_globin_family.remote(force_fetch=force_fetch)
    seqs = fam["seqs"]
    print(f"      {fam['meta']['n_family_members']} globin-domain proteins "
          f"of {fam['meta']['n_proteins_searched']} searched")

    # GPU embedding and CPU phylogeny are independent: run them concurrently.
    emb_call = None if no_gpu else esm2_cosine.spawn(seqs)
    print("[2/5] MAFFT L-INS-i + IQ-TREE 2 ..." + ("" if no_gpu else "  (ESM2 in parallel)"))
    phylo = align_tree_identity.remote(seqs)
    print(f"      best-fit model: {phylo['meta']['best_model']}")

    emb = None
    if emb_call is not None:
        print("[3/5] ESM2-650M embeddings ...")
        emb = emb_call.get()
        print(f"      {emb['meta']['embedding_dim']}-d on {emb['meta']['device']}")

    print("[4/5] GFF3 genomic context ...")
    ctx = paralog_pairs.remote(
        fam["members_csv"],
        fam["genes_parquet"],
        phylo["identity_pairs_parquet"],
        emb["cosine_matrix_csv"] if emb else None,
    )
    print(f"      duplication modes: {ctx['meta']['mode_counts']}")

    print("[5/5] writing results ...")
    files: dict[str, bytes] = {
        "globins.faa": fam["fasta"],
        "globin_family_members.csv": fam["members_csv"],
        "globins.aln.faa": phylo["alignment"],
        "globins.treefile": phylo["treefile"],
        "globins.iqtree": phylo["iqtree_report"],
        "pairwise_identity_matrix.csv": phylo["identity_matrix_csv"],
        "paralog_pairs.csv": ctx["paralog_pairs_csv"],
        "gene_context.csv": ctx["gene_context_csv"],
    }
    if emb:
        files["esm2_cosine_distance_matrix.csv"] = emb["cosine_matrix_csv"]
        files["esm2_embeddings.npz"] = emb["embeddings_npz"]
    for name, blob in files.items():
        (out / name).write_bytes(blob)

    manifest = {
        "app": APP_NAME,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "assembly_annotation": "Wm82.gnm4.ann1.T8TQ (= Phytozome Wm82.a4.v1)",
        "sources": {
            "proteome": "https://data.soybase.org/Glycine/max/annotations/"
                        "Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz",
            "gff3": "https://data.soybase.org/Glycine/max/annotations/"
                    "Wm82.gnm4.ann1.T8TQ/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_exons.gff3.gz",
            "pfam_hmm": "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00042?annotation=hmm",
        },
        "stages": {
            "select_globin_family": fam["meta"],
            "align_tree_identity": phylo["meta"],
            "esm2_cosine": emb["meta"] if emb else "skipped",
            "paralog_pairs": ctx["meta"],
        },
        "wall_seconds": round(time.time() - wall0, 1),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"done in {manifest['wall_seconds']}s -> {out.resolve()}")
    for name in sorted(files):
        print(f"  {name}")
