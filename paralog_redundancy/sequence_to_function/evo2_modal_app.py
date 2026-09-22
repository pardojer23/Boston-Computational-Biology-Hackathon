"""Score the leghemoglobin CDS / chimera set with Evo 2 on a Modal GPU.

Evo 2 is a DNA language model (Brixi et al. 2025); the `evo2` package is built on
CUDA kernels, so there is no CPU path. This job scores 65 short sequences —
5 reference CDS plus 60 paralog chimeras — and returns their mean per-token
log-likelihoods. All the analysis happens afterwards on the returned numbers.

The scoring design mirrors the ESM-2 masked-marginal analysis, one level down:

    for an ordered pair A -> B, take A's CDS and replace aligned codons with B's
      all     every differing codon          -> total divergence burden
      syn     only synonymous differences    -> INVISIBLE to a protein model
      nonsyn  only amino-acid-changing ones  -> the part ESM-2 also sees

    burden = ll(chimera) - ll(reference A)

The `syn` chimeras are the point of the experiment: they translate to exactly
the same protein as the reference, so ESM-2 scores them as identical by
construction, while Evo 2 can score them at all.

Usage (from the conda env where `modal run` already works):

    conda activate modal
    modal run evo2_modal_app.py --sequences evo2_sequences.json
    modal run evo2_modal_app.py --sequences evo2_sequences.json --model evo2_1b_base

Caveats I could not test from the analysis sandbox (no GPU, no Modal access):
  * The image build compiles flash-attn and can take 20-40 min the first time.
    If it fails, the usual fix is installing flash-attn on its own with
    --no-build-isolation before evo2 (see FLASH_ATTN_NOTE below).
  * First run downloads ~15 GB of weights into the Modal Volume; later runs
    reuse it. Model load is ~5-7 min even when cached.
"""

import modal

FLASH_ATTN_NOTE = "if the build fails: .pip_install('flash-attn', extra_options='--no-build-isolation')"

image = (
    modal.Image.from_registry("nvidia/cuda:12.4.1-devel-ubuntu22.04", add_python="3.12")
    .apt_install("git", "build-essential")
    .pip_install("packaging", "ninja", "wheel", "setuptools", "huggingface_hub")
    .pip_install("evo2==0.5.5")
)

app = modal.App("evo2-leghemoglobin")
weights = modal.Volume.from_name("evo2-weights", create_if_missing=True)


@app.function(image=image, gpu="A100-40GB", timeout=7200,
              volumes={"/weights": weights})
def score(payload_json: str, model_name: str = "evo2_7b") -> str:
    """Return JSON: {"model":…, "scores": {record_id: mean_log_likelihood}}."""
    import json
    import os
    import time

    os.environ["HF_HOME"] = "/weights/hf"          # persist weights across runs
    os.makedirs("/weights/hf", exist_ok=True)

    payload = json.loads(payload_json)
    records = payload["records"]
    print(f"{len(records)} sequences, lengths "
          f"{min(len(r['sequence']) for r in records)}-"
          f"{max(len(r['sequence']) for r in records)} nt", flush=True)

    import torch
    print("torch", torch.__version__, "| cuda", torch.cuda.is_available(),
          "|", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "no gpu", flush=True)

    from evo2 import Evo2
    t0 = time.monotonic()
    model = Evo2(model_name)
    print(f"model {model_name} loaded in {time.monotonic()-t0:.0f}s", flush=True)

    ids = [r["id"] for r in records]
    seqs = [r["sequence"] for r in records]
    t1 = time.monotonic()
    lls = model.score_sequences(seqs)              # list[float], mean per-token ll
    print(f"scored {len(seqs)} sequences in {time.monotonic()-t1:.0f}s", flush=True)

    scores = {i: float(v) for i, v in zip(ids, lls)}
    weights.commit()
    return json.dumps({
        "model": model_name,
        "torch": torch.__version__,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "n_sequences": len(seqs),
        "scores": scores,
    }, indent=1)


@app.local_entrypoint()
def main(sequences: str = "evo2_sequences.json",
         out: str = "evo2_scores.json",
         model: str = "evo2_7b"):
    import json
    import pathlib

    src = pathlib.Path(sequences).expanduser()
    if not src.is_file():
        raise SystemExit(f"sequence file not found: {src}")
    payload = src.read_text()
    n = len(json.loads(payload)["records"])
    print(f"shipping {n} sequences to Modal ({model}); "
          "first run builds the image and downloads ~15 GB of weights")

    result = score.remote(payload, model)
    pathlib.Path(out).expanduser().write_text(result)
    d = json.loads(result)
    print(f"\nwrote {out}  ({d['n_sequences']} scores, gpu={d['gpu']})")

    s = d["scores"]
    refs = {k: v for k, v in s.items() if k.startswith("ref_")}
    print("\nreference CDS mean log-likelihoods:")
    for k, v in sorted(refs.items(), key=lambda kv: -kv[1]):
        print(f"   {k:12s} {v:+.5f}")
    print("\nlargest synonymous-only burdens (invisible to a protein model):")
    syn = {k: v - s["ref_" + k.split("_to_")[0]]
           for k, v in s.items() if k.endswith("_syn")}
    for k, v in sorted(syn.items(), key=lambda kv: kv[1])[:5]:
        print(f"   {k:24s} delta_ll = {v:+.5f}")
