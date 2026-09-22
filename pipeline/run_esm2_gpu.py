"""Stage 4 as a standalone remote job: ESM2-650M mean-pooled embeddings on GPU.

Reads ./globins.faa, writes ./out/{esm2_cosine_distance_matrix.csv,
esm2_embeddings.npz, esm2_meta.json}.
"""
import json
import os
import pathlib
import sys
import time

os.environ.setdefault("HF_HUB_DISABLE_XET", "1")

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import soy_globin_core as core

out = pathlib.Path("out")
out.mkdir(exist_ok=True)

seqs = core.read_fasta("globins.faa")
t0 = time.time()
labels, emb = core.esm2_embeddings(seqs, core.ESM2_MODEL, batch_size=8)
elapsed = time.time() - t0

cos = core.cosine_distance_matrix(labels, emb)
cos.to_csv(out / "esm2_cosine_distance_matrix.csv")
np.savez_compressed(out / "esm2_embeddings.npz",
                    labels=np.array(labels), embeddings=emb)
(out / "esm2_meta.json").write_text(json.dumps({
    "model": core.ESM2_MODEL,
    "n_sequences": len(labels),
    "embedding_dim": int(emb.shape[1]),
    "pooling": "mean over residue tokens (BOS/EOS/PAD excluded)",
    "dtype": "float16" if torch.cuda.is_available() else "float32",
    "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
    "torch": torch.__version__,
    "seconds": round(elapsed, 1),
}, indent=2))

print(f"{elapsed:.1f}s  {emb.shape}  {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu'}")
print(cos.round(4).to_string())
