#!/usr/bin/env python
"""Module 4 runner — integration and the redundancy score.

    python pipeline/run_integration.py

Reads the three pair tables under ``results/`` and writes
``results/integration_module/``. No network, no GPU, no reference data: this
module only consumes what modules 1-3 already emitted.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ROOT = HERE.parent

import soy_globin_integration as integ  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default=str(ROOT / "results"))
    ap.add_argument("--outdir", default=str(ROOT / "results" / "integration_module"))
    ap.add_argument("--alpha", type=float, default=integ.ALPHA_M,
                    help="exponent on the molecular term M")
    args = ap.parse_args()

    if args.alpha != integ.ALPHA_M:
        integ.ALPHA_M = args.alpha
        integ.BETA_E = 1.0 - args.alpha

    t0 = time.time()
    print(f"[1/3] joining pair tables from {args.results}", flush=True)
    r = integ.run(args.results, args.outdir)
    scores, sweep, checks = r["scores"], r["sweep"], r["manifest"]["validation_checks"]
    print(f"      {len(scores)} pairs, {len(scores.columns)} columns", flush=True)

    print(f"[2/3] scored with alpha={integ.ALPHA_M} beta={integ.BETA_E}", flush=True)
    lb = scores[scores.pair_class == "Lb-Lb"]
    print("      focal pairs, absolute score:")
    for r_ in lb.itertuples():
        print(f"        {r_.symbol_a:<5}-{r_.symbol_b:<5} "
              f"M={r_.M_family:.3f} E={r_.E_family:.3f} R={r_.R_family:.3f}   "
              f"(clade R={r_.R_clade:.3f})")

    print("[3/3] checks", flush=True)
    for name, c in checks.items():
        verdict = c.get("passed")
        tail = "" if verdict is None else f"  passed={verdict}"
        print(f"      {name}{tail}")
    print(f"      top pair is weight-dependent: {sweep.top_lb_pair.nunique() > 1}")
    print(f"\ndone in {round(time.time() - t0, 1)}s -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
