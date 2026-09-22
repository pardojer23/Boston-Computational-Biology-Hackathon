"""Sequence-to-function divergence of the soybean leghemoglobins with ESM-2.

Two model-based quantities, both computed from the protein language model
ESM-2 650M (Lin et al. 2023; fair-esm), plus the bookkeeping needed to compare
them against the sequence-identity redundancy features of De Kegel & Ryan 2019.

  1. Representation divergence
     Mean-pooled final-layer embedding per protein; pairwise cosine distance.
     This is the standard "are these two proteins alike to the model" measure.
     It is largely a smooth function of sequence identity, so on its own it is
     a sanity check rather than new information.

  2. Substitution constraint  (the part identity cannot see)
     For an ordered pair A -> B, align the two proteins, and at every position
     where they differ compute the masked-marginal log-odds in A's OWN context:

         s = log P(B_res | A with position masked) - log P(A_res | same)

     This is the Meier et al. 2021 variant-effect score, applied to the
     substitutions that actually separate two paralogs. Summed over differing
     positions it asks: does turning A into B walk through positions the model
     considers constrained? Two pairs with identical percent identity can
     differ substantially here, which is exactly where a sequence-to-function
     model adds signal over a percent-identity feature.

Sign convention: more negative = the substitutions separating the pair sit at
positions where the model strongly prefers the resident residue, i.e. greater
predicted functional divergence.

CPU-only; ~5 proteins of 145-168 aa.
"""
from __future__ import annotations

import itertools
import json
import os
import pathlib

os.environ.setdefault("TORCH_HOME", str(pathlib.Path("torch_home").resolve()))

import numpy as np
import pandas as pd
import torch
from Bio import Align
from Bio.Align import substitution_matrices

MODEL_NAME = os.environ.get("ESM_MODEL", "esm2_t33_650M_UR50D")
REPR_LAYER = int(MODEL_NAME.split("_")[1][1:])   # esm2_t33_... -> 33

NAMES = {
    "GLYMA_10G198800": "LB1",
    "GLYMA_10G198900": "LB5",
    "GLYMA_10G199000": "LB2",
    "GLYMA_10G199100": "LB3",
    "GLYMA_20G191200": "LB4",
}


def load_sequences(d="data/prot") -> dict[str, str]:
    seqs = {}
    for gene in NAMES:
        p = pathlib.Path(d) / f"{gene}.fa"
        seqs[gene] = "".join(l.strip() for l in p.read_text().splitlines()
                             if not l.startswith(">"))
    return seqs


def align_pairs(seqs: dict[str, str]) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Global pairwise alignment; return aligned (i, j) index pairs, gaps dropped."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -11, -1
    aligner.mode = "global"
    out = {}
    for a, b in itertools.combinations(seqs, 2):
        aln = aligner.align(seqs[a], seqs[b])[0]
        pairs = []
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            pairs += [(a0 + k, b0 + k) for k in range(a1 - a0)]
        out[(a, b)] = pairs
    return out


def masked_logprobs(model, alphabet, seq: str, positions: list[int],
                    device="cpu") -> dict[int, "torch.Tensor"]:
    """log-softmax over the vocabulary at each requested position, that position masked.

    Kept as torch tensors throughout: the available torch wheel for this
    platform was built against NumPy 1.x, so .numpy() raises on NumPy 2.
    """
    bc = alphabet.get_batch_converter()
    _, _, toks = bc([("q", seq)])
    out = {}
    for pos in positions:                      # pos is 0-based within the sequence
        t = toks.clone()
        t[0, pos + 1] = alphabet.mask_idx      # +1 for BOS
        with torch.no_grad():
            logits = model(t.to(device))["logits"]
        out[pos] = torch.log_softmax(logits[0, pos + 1], dim=-1).cpu()
    return out


def main():
    seqs = load_sequences()
    genes = list(seqs)
    print(f"{len(seqs)} proteins: " +
          ", ".join(f"{NAMES[g]}={len(s)}aa" for g, s in seqs.items()), flush=True)

    import esm
    model, alphabet = getattr(esm.pretrained, MODEL_NAME)()
    model = model.eval()
    torch.set_num_threads(max(1, (os.cpu_count() or 2) - 1))

    # ---- 1. representation divergence
    bc = alphabet.get_batch_converter()
    emb = {}
    for g, s in seqs.items():
        _, _, toks = bc([(g, s)])
        with torch.no_grad():
            r = model(toks, repr_layers=[REPR_LAYER])["representations"][REPR_LAYER]
        emb[g] = r[0, 1:-1].mean(0)            # drop BOS/EOS, mean-pool (torch tensor)
        print(f"  embedded {NAMES[g]}", flush=True)

    # ---- alignments, and the union of positions needing a masked pass
    aligned = align_pairs(seqs)
    need: dict[str, set[int]] = {g: set() for g in genes}
    for (a, b), pairs in aligned.items():
        for i, j in pairs:
            if seqs[a][i] != seqs[b][j]:
                need[a].add(i)
                need[b].add(j)
    print("masked passes required: " +
          ", ".join(f"{NAMES[g]}={len(need[g])}" for g in genes), flush=True)

    # ---- 2. masked marginals at those positions
    lp = {}
    for g in genes:
        lp[g] = masked_logprobs(model, alphabet, seqs[g], sorted(need[g]))
        print(f"  scored {NAMES[g]} ({len(need[g])} positions)", flush=True)

    tok = {aa: alphabet.get_idx(aa) for aa in set("".join(seqs.values()))}

    rows = []
    for (a, b), pairs in aligned.items():
        subs = [(i, j) for i, j in pairs if seqs[a][i] != seqs[b][j]]
        n_aln = len(pairs)
        # identity over aligned columns, and the gapped length difference
        pid_aln = 100.0 * (n_aln - len(subs)) / n_aln
        s_ab = [(lp[a][i][tok[seqs[b][j]]] - lp[a][i][tok[seqs[a][i]]]).item() for i, j in subs]
        s_ba = [(lp[b][j][tok[seqs[a][i]]] - lp[b][j][tok[seqs[b][j]]]).item() for i, j in subs]
        cos = 1.0 - torch.nn.functional.cosine_similarity(
            emb[a].unsqueeze(0), emb[b].unsqueeze(0)).item()
        rows.append(dict(
            gene_A=a, gene_B=b, symbol_A=NAMES[a], symbol_B=NAMES[b],
            len_A=len(seqs[a]), len_B=len(seqs[b]),
            aligned_cols=n_aln, n_substitutions=len(subs),
            pct_identity_aligned=round(pid_aln, 4),
            embed_cosine_distance=round(cos, 6),
            llr_sum_A_to_B=round(sum(s_ab), 4), llr_sum_B_to_A=round(sum(s_ba), 4),
            llr_mean_A_to_B=round(float(np.mean(s_ab)), 4) if s_ab else 0.0,
            llr_mean_B_to_A=round(float(np.mean(s_ba)), 4) if s_ba else 0.0,
            llr_mean_symmetric=round(float(np.mean(s_ab + s_ba)), 4) if s_ab else 0.0,
        ))
    pairs_df = pd.DataFrame(rows)

    # ---- per-gene: closest paralog under each measure
    recs = []
    for g in genes:
        sub = pairs_df[(pairs_df.gene_A == g) | (pairs_df.gene_B == g)].copy()
        sub["other"] = np.where(sub.gene_A == g, sub.gene_B, sub.gene_A)
        sub["llr_own_context"] = np.where(sub.gene_A == g,
                                          sub.llr_mean_A_to_B, sub.llr_mean_B_to_A)
        nearest_emb = sub.loc[sub.embed_cosine_distance.idxmin()]
        least_div = sub.loc[sub.llr_own_context.idxmax()]   # least negative = mildest
        recs.append(dict(
            gene=g, symbol=NAMES[g], protein_len_aa=len(seqs[g]),
            min_embed_cosine_distance=round(float(sub.embed_cosine_distance.min()), 6),
            closest_by_embedding=NAMES[nearest_emb.other],
            best_llr_mean_own_context=round(float(sub.llr_own_context.max()), 4),
            mildest_divergence_partner=NAMES[least_div.other],
            mean_llr_own_context=round(float(sub.llr_own_context.mean()), 4),
        ))
    genes_df = pd.DataFrame(recs)

    pairs_df.to_csv("esm_pair_divergence.csv", index=False)
    genes_df.to_csv("esm_gene_divergence.csv", index=False)
    np.savez("esm_embeddings.npz",
             **{NAMES[g]: np.array(emb[g].tolist(), dtype="float32") for g in genes})
    json.dump({"model": MODEL_NAME, "repr_layer": REPR_LAYER,
               "torch": torch.__version__,
               "score": "masked-marginal log-odds (Meier et al. 2021) at "
                        "positions differing between aligned paralogs",
               "sign": "more negative = greater predicted functional divergence"},
              open("esm_run_metadata.json", "w"), indent=2)
    print("\n", pairs_df[["symbol_A", "symbol_B", "n_substitutions",
                          "pct_identity_aligned", "embed_cosine_distance",
                          "llr_mean_symmetric"]].to_string(index=False))
    print("\n", genes_df.to_string(index=False))


if __name__ == "__main__":
    main()
