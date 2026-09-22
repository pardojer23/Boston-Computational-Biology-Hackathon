"""Module 3 — expression (pseudobulk).

Emits per-gene expression profiles and a pairwise table keyed on
``label_a``/``label_b``, so it joins onto ``paralog_pairs.csv`` and
``structure_pairs.csv`` exactly as the other modules do.

Two decisions are recorded here because they are not obvious from the output:

* **The atlas named in NEXT_RUN.md is not used.** GSE270392 (Zhang 2024 Cell)
  is snRNA-seq of *early* nodule. In its 29,436 nuclei Lba carries 668 total
  counts against a median gene total of 2,496 — rank 45,386 of 52,594 — and
  the most abundant globin in the sample is the chr11 non-symbiotic
  haemoglobin. Nuclear RNA is depleted of abundant stable cytoplasmic
  transcripts and leghemoglobin is the extreme case. Any co-expression
  statistic computed there rests on 9-33 co-detected nuclei per Lb pair and
  ranks the intended negative control top. ``dataset_diagnostic.csv`` holds
  the comparison.

* **GSE226149 (mature nodule + root, protoplast scRNA-seq) is used instead**,
  as *pseudobulk*: counts are summed over every barcode in a library, so no
  cell calling, clustering or cell-type annotation enters this module. That is
  a deliberate scope limit, not an oversight — the published cell-type labels
  are not in the GEO submission, and inferring them here would substitute a
  different piece of work for the one the pipeline design asks for.

Consequence of pseudobulk, stated rather than worked around: ``detection_rate``
and ``coexpression_overlap`` are defined per *cell* (fraction of cells
expressing both over either). Pseudobulk has no cells. They are emitted as
NA with a reason rather than replaced by a look-alike statistic.

No Modal import, for the same reason ``soy_globin_core`` has none.
"""

from __future__ import annotations

import gzip
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

#: GSE226149 libraries: accession -> (file stem, tissue)
SAMPLES = {
    "GSM7065810": ("Nodule_sample1", "nodule"),
    "GSM7065811": ("Nodule_sample2", "nodule"),
    "GSM7065807": ("Root_sample1", "root"),
    "GSM7065808": ("Root_sample2", "root"),
    "GSM7065809": ("Root_sample3", "root"),
}

#: The globin family, as labelled by the sequence module.
FAMILY = {
    "Glyma.10G199100": "Lba",
    "Glyma.10G199000": "Lbc1",
    "Glyma.20G191200": "Lbc2",
    "Glyma.10G198800": "Lbc3",
    "Glyma.10G198900": "GmLb5",
    "Glyma.11G121700": None,
    "Glyma.11G121800": None,
}

#: Independent nodule cell-type markers, resolved UniProt -> Wm82.a4 by phmmer
#: against the primary proteome. Leghemoglobin is deliberately excluded: using
#: Lb to mark the infected-cell compartment would be circular in a study of Lb.
#: Reported as context for the pseudobulk profiles, not used to define them.
MARKERS = {
    "Glyma.08G120100": "NOD26 (symbiosome membrane, infected cells; P08995)",
    "Glyma.10G121524": "Uricase-2 / nodulin-35 (uninfected interstitial cells; P04670)",
    "Glyma.20G072400": "Uricase-2 isozyme 2 (P04104)",
    "Glyma.13G114000": "Sucrose synthase / nodulin-100 (P13708)",
}

#: GEO supplementary-file base for a GSM accession.
GEO_SAMPLE_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{dir}/{gsm}/suppl"

#: The three CellRanger parts each library ships.
PARTS = ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")


def fetch_inputs(datadir: str | Path, force: bool = False) -> dict[str, str]:
    """Download the GSE226149 libraries into ``datadir`` under canonical names.

    The filenames this module reads are exactly GEO's own
    ``{GSM}_{stem}_{part}`` — no local renaming step. Existing files are kept
    unless ``force``, so a rerun is cheap (~1.2 GB total).
    """
    import urllib.request

    datadir = Path(datadir)
    datadir.mkdir(parents=True, exist_ok=True)
    got = {}
    for gsm, (stem, _) in SAMPLES.items():
        base = GEO_SAMPLE_BASE.format(dir=f"{gsm[:-3]}nnn", gsm=gsm)
        for part in PARTS:
            name = f"{gsm}_{stem}_{part}"
            dest = datadir / name
            if force or not dest.exists():
                urllib.request.urlretrieve(f"{base}/{name}", dest)
            got[name] = str(dest)
    return got


def check_inputs(datadir: str | Path) -> list[str]:
    """Return the canonical filenames that are missing from ``datadir``."""
    datadir = Path(datadir)
    return [f"{gsm}_{stem}_{part}"
            for gsm, (stem, _) in SAMPLES.items() for part in PARTS
            if not (datadir / f"{gsm}_{stem}_{part}").exists()]


#: GSE226149 features.tsv uses GLYMA_10G199100; everything else uses Glyma.10G199100.
def to_geo_id(gene_id: str) -> str:
    return gene_id.replace(".", "_").upper()


def label_of(gene_id: str) -> str:
    sym = FAMILY.get(gene_id)
    return f"{gene_id}_{sym}" if sym else gene_id


# --------------------------------------------------------------------------- #
# Pseudobulk
# --------------------------------------------------------------------------- #


def pseudobulk_library(datadir: Path, gsm: str, stem: str) -> pd.Series:
    """Sum counts over every barcode in one library.

    No cell calling: the sum over all droplets is the library's bulk profile.
    Ambient and cell-associated RNA both contribute, which is the point — this
    is a bulk measurement, and treating it as one avoids importing a cell-
    calling threshold that would need its own justification.
    """
    feats = pd.read_csv(
        datadir / f"{gsm}_{stem}_features.tsv.gz", sep="\t", header=None, usecols=[0]
    )[0].values
    mat = mmread(datadir / f"{gsm}_{stem}_matrix.mtx.gz").tocsr()
    if mat.shape[0] != len(feats):
        raise ValueError(f"{gsm}: {mat.shape[0]} matrix rows vs {len(feats)} features")
    return pd.Series(np.asarray(mat.sum(axis=1)).ravel(), index=feats, name=gsm)


def build_pseudobulk(datadir: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (raw counts, CPM) as genes x libraries."""
    datadir = Path(datadir)
    cols = {}
    for gsm, (stem, _) in SAMPLES.items():
        cols[gsm] = pseudobulk_library(datadir, gsm, stem)
    counts = pd.DataFrame(cols)
    cpm = counts / counts.sum(axis=0) * 1e6
    return counts, cpm


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #


def tau(profile: np.ndarray) -> float:
    """Yanai tissue-specificity index.

    tau = sum(1 - x_i / x_max) / (n - 1), computed on non-negative expression.
    Returns NaN when the profile is all-zero (tau undefined) or n < 2.
    """
    x = np.asarray(profile, dtype=float)
    if x.size < 2 or not np.isfinite(x).all() or x.max() <= 0:
        return float("nan")
    return float((1.0 - x / x.max()).sum() / (x.size - 1))


def gene_profiles(cpm: pd.DataFrame, genes: dict[str, str | None]) -> pd.DataFrame:
    """Per-gene pseudobulk profile, tissue means, and tau over tissue means."""
    tissues = sorted({t for _, t in SAMPLES.values()})
    rows = []
    for gene in genes:
        gid = to_geo_id(gene)
        if gid not in cpm.index:
            rows.append({"gene_id": gene, "label": label_of(gene), "in_matrix": False})
            continue
        v = cpm.loc[gid]
        tmeans = {t: float(v[[g for g, (_, tt) in SAMPLES.items() if tt == t]].mean())
                  for t in tissues}
        row = {
            "gene_id": gene,
            "label": label_of(gene),
            "symbol": genes[gene] or "",
            "in_matrix": True,
            **{f"cpm_{g}": float(v[g]) for g in SAMPLES},
            **{f"cpm_mean_{t}": tmeans[t] for t in tissues},
            "tau_over_tissue_means": tau(np.array([tmeans[t] for t in tissues])),
            "n_tissues_for_tau": len(tissues),
            # cell-level quantities, not defined on pseudobulk
            "detection_rate": np.nan,
            "detection_rate_note": "not computable from pseudobulk (no cells)",
        }
        rows.append(row)
    return pd.DataFrame(rows)


def pair_metrics(cpm: pd.DataFrame, genes: dict[str, str | None]) -> pd.DataFrame:
    """Pairwise profile correlation over the pseudobulk libraries."""
    present = [g for g in genes if to_geo_id(g) in cpm.index]
    logcpm = np.log1p(cpm)
    rows = []
    for a, b in combinations(present, 2):
        va = logcpm.loc[to_geo_id(a), list(SAMPLES)].to_numpy(dtype=float)
        vb = logcpm.loc[to_geo_id(b), list(SAMPLES)].to_numpy(dtype=float)
        if np.ptp(va) == 0 or np.ptp(vb) == 0:
            rho = np.nan
        else:
            rho = float(spearmanr(va, vb).statistic)
        rows.append({
            "label_a": label_of(a),
            "label_b": label_of(b),
            "gene_a": a,
            "gene_b": b,
            "spearman_profile": rho,
            "n_profile_points": len(SAMPLES),
            "coexpression_overlap": np.nan,
            "coexpression_overlap_note": "not computable from pseudobulk (no cells)",
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def run(datadir: str | Path, outdir: str | Path, fetch: bool = True) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    missing = check_inputs(datadir)
    if missing:
        if not fetch:
            raise FileNotFoundError(
                f"{len(missing)} input file(s) missing from {datadir}, "
                f"first: {missing[0]}. Run with fetch=True (or drop --no-fetch)."
            )
        fetch_inputs(datadir)
        still = check_inputs(datadir)
        if still:
            raise FileNotFoundError(f"fetch did not produce: {still}")

    counts, cpm = build_pseudobulk(datadir)
    prof = gene_profiles(cpm, FAMILY)
    pairs = pair_metrics(cpm, FAMILY)
    marks = gene_profiles(cpm, {g: None for g in MARKERS})
    marks["marker_role"] = [MARKERS[g] for g in marks["gene_id"]]

    prof.to_csv(outdir / "gene_pseudobulk_profiles.csv", index=False)
    pairs.to_csv(outdir / "expression_pairs.csv", index=False)
    marks.to_csv(outdir / "marker_profiles.csv", index=False)

    manifest = {
        "source_series": "GSE226149",
        "source_note": "mature nodule + root, protoplast scRNA-seq, used as pseudobulk",
        "rejected_series": "GSE270392",
        "rejected_reason": (
            "snRNA-seq of early nodule; Lba ranks 45,386/52,594 with 668 counts "
            "against a median gene total of 2,496; 9-33 co-detected nuclei per Lb pair"
        ),
        "libraries": {g: {"stem": s, "tissue": t, "total_counts": int(counts[g].sum())}
                      for g, (s, t) in SAMPLES.items()},
        "normalisation": "CPM over summed library counts; no cell calling",
        "metrics_emitted": ["cpm per library", "cpm mean per tissue",
                            "tau over tissue means", "spearman of log1p-CPM profiles"],
        "metrics_not_computable": {
            "detection_rate": "per-cell quantity; pseudobulk has no cells",
            "coexpression_overlap": "per-cell quantity; pseudobulk has no cells",
        },
        "tau_caveat": "two tissue types only; tau collapses to a nodule-vs-root contrast",
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return {"counts": counts, "cpm": cpm, "profiles": prof, "pairs": pairs,
            "markers": marks, "manifest": manifest}


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datadir", default="data/expression")
    ap.add_argument("--outdir", default="results/expression_module")
    ap.add_argument("--no-fetch", action="store_true",
                    help="fail instead of downloading missing libraries")
    a = ap.parse_args()
    r = run(a.datadir, a.outdir, fetch=not a.no_fetch)
    print(r["profiles"].to_string(index=False))
