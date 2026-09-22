#!/usr/bin/env python
"""Module 2 runner — structure.

    python pipeline/run_structure.py

Reads ``results/sequence_module/`` (family membership + the MAFFT alignment),
writes ``results/structure_module/``. Network access is needed for UniProt,
AlphaFold DB and RCSB; everything downloaded is cached under
``work/structure/`` so re-runs are offline.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent

import soy_globin_core as core          # noqa: E402
import soy_globin_structure as struc    # noqa: E402

SEQ_DIR = ROOT / "results" / "sequence_module"
OUT_DIR = ROOT / "results" / "structure_module"
WORK = ROOT / "work" / "structure"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", type=float, default=struc.POCKET_CUTOFF_A,
                    help="heme-contact distance in angstrom")
    ap.add_argument("--template", default=struc.POCKET_TEMPLATE_PDB,
                    help="heme-bound PDB id used to define the pocket")
    args = ap.parse_args()

    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "pdb").mkdir(exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    members = pd.read_csv(SEQ_DIR / "globin_family_members.csv")
    aln = core.read_fasta(SEQ_DIR / "globins.aln.faa")
    prot = core.read_fasta(SEQ_DIR / "globins.faa")
    labels = list(prot)
    gene_by_label = {lab: lab.split("_")[0] for lab in labels}

    # ---- 1. gene -> UniProt -----------------------------------------------
    print("[1/6] UniProt accessions")
    up = struc.resolve_uniprot([gene_by_label[l] for l in labels])
    up["label"] = labels
    missing = up[up.uniprot.isna()].gene_id.tolist()
    if missing:
        raise SystemExit(f"no UniProt accession for: {missing}")
    for _, r in up.iterrows():
        print(f"      {r.label:<28} {r.uniprot:<12} {r.entry_name:<14} "
              f"{r.uniprot_length:>4} aa  {r.uniprot_reviewed}")

    # ---- 2. AlphaFold DB --------------------------------------------------
    print("[2/6] AlphaFold DB models")
    prov = []
    for _, r in up.iterrows():
        meta = struc.fetch_afdb(r.uniprot, OUT_DIR / "pdb")
        meta["label"] = r.label
        meta["gene_id"] = r.gene_id
        prov.append(meta)
        print(f"      {r.label:<28} {meta['afdb_entry']:<20} "
              f"pLDDT {meta['mean_plddt']:.1f}  ({meta['afdb_model_created']})")
    prov = pd.DataFrame(prov)

    # ---- 3. reconcile AFDB sequence against the a4 protein ----------------
    print("[3/6] AFDB sequence vs Wm82.a4 primary transcript")
    rec = []
    for _, r in prov.iterrows():
        struct_coords, struct_seq, struct_resnums = struc.read_pdb_chain(r.pdb_path)
        a4 = prot[r.label]
        pid, ncol = struc.pairwise_identity_global(struct_seq, a4)
        rec.append({
            "label": r.label,
            "n_res_model": len(struct_seq),
            "len_a4_protein": len(a4),
            "len_delta": len(struct_seq) - len(a4),
            "pid_model_vs_a4": round(pid, 2),
            "n_aligned_cols": ncol,
            "model_is_a4_sequence": bool(struct_seq == a4),
        })
        flag = "" if struct_seq == a4 else "   <- differs from a4"
        print(f"      {r.label:<28} {len(struct_seq):>4} res  "
              f"{pid:6.2f}% vs a4 ({len(a4)} aa){flag}")
    rec = pd.DataFrame(rec)

    structures = prov.merge(up.drop(columns=["gene_id"]), on="label").merge(rec, on="label")

    # ---- 4. heme pocket ---------------------------------------------------
    print(f"[4/6] heme pocket from {args.template} (<= {args.cutoff} A)")
    pk = struc.heme_pocket_residues(
        pdb_id=args.template, cutoff=args.cutoff, cache_dir=WORK
    )
    print(f"      chain {pk['chain']}, {len(pk['pocket_resnums'])} residues: "
          f"{pk['pocket_resnums']}")

    # The crystal is leghemoglobin a; transfer its pocket numbering onto the
    # MSA through whichever label is Lba.
    lba_label = next(l for l in labels if l.endswith("_Lba"))
    pos2col = struc.map_positions_through_alignment(
        pk["chain_seq"], pk["pocket_resnums"], aln, lba_label
    )
    dropped = [p for p in pk["pocket_resnums"] if p not in pos2col]
    if dropped:
        print(f"      {len(dropped)} pocket residue(s) unmappable to the MSA: {dropped}")
    cols = sorted(pos2col.values())
    print(f"      mapped to {len(cols)} alignment columns")

    # Independent check: UniProt's own binding-site annotations for Lba.
    lba_acc = up.loc[up.label == lba_label, "uniprot"].iloc[0]
    sites = struc.uniprot_ligand_sites(lba_acc)
    site_cols = struc.map_positions_through_alignment(
        prot[lba_label], [s["position"] for s in sites], aln, lba_label
    )
    for s in sites:
        col = site_cols.get(s["position"])
        inside = col in cols if col is not None else None
        part = f"/{s['ligand_part']}" if s["ligand_part"] else ""
        print(f"      UniProt {lba_acc} pos {s['position']:>3} "
              f"{s['ligand']}{part:<4} {s['description']:<26} "
              f"in crystal pocket: {inside}")
    n_conf = sum(1 for s in sites if site_cols.get(s["position"]) in cols)

    pocket_tbl = []
    for lab in labels:
        row = {"label": lab}
        for p, col in sorted(pos2col.items()):
            row[f"col{col + 1}_lba{p}"] = aln[lab][col]
        pocket_tbl.append(row)
    pd.DataFrame(pocket_tbl).to_csv(OUT_DIR / "pocket_residues.csv", index=False)

    # ---- 5. all-vs-all superposition -------------------------------------
    print("[5/6] TM-align all-vs-all")
    from tmtools import tm_align

    chains = {}
    for _, r in structures.iterrows():
        c, s, _ = struc.read_pdb_chain(r.pdb_path)
        chains[r.label] = (c, s)

    rows = []
    for a, b in itertools.combinations(labels, 2):
        ca, sa = chains[a]
        cb, sb = chains[b]
        res = tm_align(ca, cb, sa, sb)
        # TM-score normalised by the shorter chain is the conservative choice:
        # it cannot be inflated by one protein being a fragment of the other.
        tm1, tm2 = float(res.tm_norm_chain1), float(res.tm_norm_chain2)
        na, nb = len(sa), len(sb)
        tm_short = tm1 if na <= nb else tm2

        pa, pb = [aln[a][c] for c in cols], [aln[b][c] for c in cols]
        both = [(x, y) for x, y in zip(pa, pb) if x != "-" and y != "-"]
        pocket_id = 100.0 * sum(1 for x, y in both if x == y) / len(both) if both else np.nan

        rows.append({
            "label_a": a, "label_b": b,
            "tm_score": round(tm_short, 4),
            "tm_norm_a": round(tm1, 4), "tm_norm_b": round(tm2, 4),
            "rmsd": round(float(res.rmsd), 3),
            "pocket_identity": round(pocket_id, 2),
            "n_pocket_cols": len(both),
            "n_res_a": na, "n_res_b": nb,
        })
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT_DIR / "structure_pairs.csv", index=False)
    structures.drop(columns=["afdb_sequence"]).to_csv(
        OUT_DIR / "structures.csv", index=False
    )

    # ---- 6. manifest ------------------------------------------------------
    lb = {l for l in labels if l.endswith(("_Lba", "_Lbc1", "_Lbc2", "_Lbc3"))}
    lb_pairs = pairs[pairs.label_a.isin(lb) & pairs.label_b.isin(lb)]
    manifest = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runner": "pipeline/run_structure.py",
        "n_structures": len(structures),
        "structure_source": "AlphaFold DB (no prediction run; no GPU used)",
        "afdb_api": struc.AFDB_API,
        "pocket": {
            "template_pdb": pk["pdb_id"], "chain": pk["chain"],
            "ligand": pk["ligand"], "cutoff_a": pk["cutoff_a"],
            "n_pocket_residues_crystal": len(pk["pocket_resnums"]),
            "pocket_resnums_crystal": pk["pocket_resnums"],
            "n_alignment_columns": len(cols),
            "unmappable_resnums": dropped,
            "uniprot_binding_sites_checked": len(sites),
            "uniprot_binding_sites_inside_pocket": n_conf,
        },
        "tm_score_convention": "normalised by the shorter chain",
        "lb_clade_tm_score": {
            "min": float(lb_pairs.tm_score.min()),
            "max": float(lb_pairs.tm_score.max()),
            "n_pairs": int(len(lb_pairs)),
        },
        "lb_clade_pocket_identity": {
            "min": float(lb_pairs.pocket_identity.min()),
            "max": float(lb_pairs.pocket_identity.max()),
        },
        "wall_seconds": round(time.time() - t0, 1),
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))

    print(f"[6/6] done in {manifest['wall_seconds']}s -> {OUT_DIR}")
    print(f"      Lb-clade TM-score {manifest['lb_clade_tm_score']['min']:.4f}"
          f"-{manifest['lb_clade_tm_score']['max']:.4f} over "
          f"{manifest['lb_clade_tm_score']['n_pairs']} pairs")
    print(f"      Lb-clade pocket identity "
          f"{manifest['lb_clade_pocket_identity']['min']:.1f}"
          f"-{manifest['lb_clade_pocket_identity']['max']:.1f}%")


if __name__ == "__main__":
    main()
