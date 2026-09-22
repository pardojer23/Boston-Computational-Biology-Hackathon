#!/usr/bin/env python
"""
Run the whole soybean globin sequence module locally, no Modal.

    python pipeline/run_local.py                      # everything, incl. ESM2 on CPU
    python pipeline/run_local.py --skip-esm2          # sequence stages only
    python pipeline/run_local.py --reuse-esm2         # keep an existing cosine matrix
    python pipeline/run_local.py --force-fetch        # re-download the references

Needs hmmer, mafft and iqtree on PATH (see environment.yml). The ESM2 stage
additionally needs torch + transformers; keep those in a separate env and
either run with --skip-esm2 here or use pipeline/run_esm2_gpu.py on its own.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

import soy_globin_core as core

REPO = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datadir", default=str(REPO / "data"))
    ap.add_argument("--workdir", default=str(REPO / "work"))
    ap.add_argument("--outdir", default=str(REPO / "results" / "sequence_module"))
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--force-fetch", action="store_true")
    esm = ap.add_mutually_exclusive_group()
    esm.add_argument("--skip-esm2", action="store_true",
                     help="omit the embedding stage entirely")
    esm.add_argument("--reuse-esm2", action="store_true",
                     help="do not recompute; reuse the cosine matrix already in --outdir")
    args = ap.parse_args()

    data, work, out = Path(args.datadir), Path(args.workdir), Path(args.outdir)
    for d in (data, work, out):
        d.mkdir(parents=True, exist_ok=True)
    wall0 = time.time()

    # -- 1. references ----------------------------------------------------- #
    print("[1/6] references", flush=True)
    paths = core.fetch_inputs(data, force=args.force_fetch)

    # -- 2. family selection ------------------------------------------------ #
    print("[2/6] hmmsearch PF00042 --cut_ga", flush=True)
    seqs, members = core.select_family(
        paths["proteome"], paths["pfam_hmm"], work / "hmmer", cpu=args.threads
    )
    n_searched = len(core.read_fasta(paths["proteome"]))
    print(f"      {len(members)} of {n_searched} proteins", flush=True)

    # -- 3. alignment + tree ------------------------------------------------ #
    print("[3/6] MAFFT L-INS-i + IQ-TREE 2", flush=True)
    phylo = work / "phylo"
    phylo.mkdir(parents=True, exist_ok=True)
    faa, aln = phylo / "globins.faa", phylo / "globins.aln.faa"
    core.write_fasta(seqs, faa)
    core.run_mafft(faa, aln, threads=args.threads)
    tree = core.run_iqtree(aln, phylo / "globins", threads=str(args.threads))
    print(f"      best-fit model: {core.best_model(tree['iqtree_report'])}", flush=True)

    # -- 4. pairwise identity ----------------------------------------------- #
    print("[4/6] pairwise identity", flush=True)
    ident_mat, ident_long = core.pairwise_identity(aln)

    # -- 5. ESM2 ------------------------------------------------------------ #
    cos = None
    esm2_meta: object = "skipped"
    cos_path = out / "esm2_cosine_distance_matrix.csv"
    if args.reuse_esm2:
        if not cos_path.exists():
            print(f"      --reuse-esm2 but {cos_path} is absent", file=sys.stderr)
            return 2
        cos = pd.read_csv(cos_path, index_col=0)
        # A reused matrix is keyed on the labels of the run that produced it. An
        # embedding is a per-*gene* quantity and the label is only presentation,
        # so reuse is remapped gene_id -> current label rather than matched on
        # the label string: adding a symbol to core.GENE_SYMBOLS would otherwise
        # KeyError here even though the sequences are untouched. Values are not
        # recomputed and not altered — only renamed.
        relabel = {c: core.label_of(core.gene_from_label(c)) for c in cos.index}
        renamed = {k: v for k, v in relabel.items() if k != v}
        cos = cos.rename(index=relabel, columns=relabel)
        esm2_meta = {
            "provenance": "reused from a previous run (see that run's manifest)",
            "relabelled": renamed or "no label changes",
        }
        print(f"[5/6] ESM2 reused from {cos_path}", flush=True)
        if renamed:
            print(f"      relabelled {len(renamed)} tip(s): {renamed}", flush=True)
            npz_path = out / "esm2_embeddings.npz"
            if npz_path.exists():
                import numpy as np

                z = np.load(npz_path, allow_pickle=False)
                np.savez_compressed(
                    npz_path,
                    labels=np.array([relabel.get(str(l), str(l)) for l in z["labels"]]),
                    embeddings=z["embeddings"],
                )
    elif not args.skip_esm2:
        print("[5/6] ESM2-650M embeddings (slow on CPU: ~10 min incl. download)", flush=True)
        import numpy as np

        labels, emb = core.esm2_embeddings(seqs)
        cos = core.cosine_distance_matrix(labels, emb)
        np.savez_compressed(out / "esm2_embeddings.npz",
                            labels=np.array(labels), embeddings=emb)
        esm2_meta = {
            "model": core.ESM2_MODEL,
            "embedding_dim": int(emb.shape[1]),
            "pooling": "mean over residue tokens (BOS/EOS/PAD excluded)",
        }
    else:
        print("[5/6] ESM2 skipped", flush=True)

    # -- 6. genomic context -------------------------------------------------- #
    print("[6/6] GFF3 context + duplication mode", flush=True)
    genes = core.parse_gff_genes(paths["gff3"])
    pairs, ctx = core.classify_pairs(members, genes,
                                     identity_long=ident_long, cosine=cos)

    # -- write ---------------------------------------------------------------- #
    core.write_fasta(seqs, out / "globins.faa")
    members.to_csv(out / "globin_family_members.csv", index=False)
    shutil.copy(aln, out / "globins.aln.faa")
    shutil.copy(tree["treefile"], out / "globins.treefile")
    shutil.copy(tree["iqtree_report"], out / "globins.iqtree")
    core.annotate_newick(tree["treefile"], out / "globins.annotated.nwk", ctx)
    ident_mat.to_csv(out / "pairwise_identity_matrix.csv")
    ident_long.to_csv(out / "pairwise_identity_pairs.csv", index=False)
    pairs.to_csv(out / "paralog_pairs.csv", index=False)
    ctx.to_csv(out / "gene_context.csv", index=False)
    if cos is not None:
        cos.to_csv(cos_path)

    manifest = {
        "assembly_annotation": "Wm82.gnm4.ann1.T8TQ (= Phytozome Wm82.a4.v1)",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runner": "pipeline/run_local.py",
        "sources": {"proteome": core.PROTEOME_URL, "gff3": core.GFF3_URL,
                    "pfam_hmm": core.PFAM_HMM_URL},
        "input_sha256": {k: core.sha256(paths[k])
                         for k in ("proteome", "gff3", "pfam_hmm")},
        "pfam": core.hmm_metadata(paths["pfam_hmm"]),
        "hmmsearch": {"threshold": "--cut_ga",
                      "n_proteins_searched": n_searched,
                      "n_hits": int(members.pfam_hit.sum()),
                      "n_focal_rescued": int((~members.pfam_hit).sum()),
                      "version": core.tool_version(["hmmsearch", "-h"])},
        "alignment": {"program": core.tool_version(["mafft", "--version"]),
                      "mode": "L-INS-i (--localpair --maxiterate 1000)",
                      "n_sequences": len(seqs),
                      "n_columns": len(next(iter(core.read_fasta(aln).values())))},
        "phylogeny": {"program": core.tool_version([core.iqtree_exe(), "--version"]),
                      "options": "-m MFP -B 1000 --alrt 1000 --seed 20240601",
                      "best_model": core.best_model(tree["iqtree_report"])},
        "esm2": esm2_meta,
        "duplication_mode": {
            "tandem_max_intervening_genes": core.TANDEM_MAX_INTERVENING,
            "proximal_max_bp": core.PROXIMAL_MAX_BP,
            "counts": pairs.duplication_mode.value_counts().to_dict()},
        "wall_seconds": round(time.time() - wall0, 1),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"\ndone in {manifest['wall_seconds']}s -> {out}")
    for p in sorted(out.iterdir()):
        print(f"  {p.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
