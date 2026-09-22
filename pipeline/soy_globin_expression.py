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
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import soy_globin_core as core  # noqa: E402

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

#: The family is *discovered* by hmmsearch in module 1, not declared, so this
#: module reads it from module 1's output rather than keeping its own copy. It
#: used to keep a copy, and the copy disagreed with ``core.GENE_SYMBOLS`` about
#: whether ``Glyma.10G198900`` carries a symbol — which silently cost 6 of 21
#: rows on any join between this module's pair table and the sequence module's.
SEQ_MEMBERS = (
    Path(__file__).resolve().parent.parent
    / "results" / "sequence_module" / "globin_family_members.csv"
)


def load_family(members_csv: str | Path = SEQ_MEMBERS) -> dict[str, str | None]:
    """Family gene IDs -> symbol, read from the sequence module's member table."""
    members_csv = Path(members_csv)
    if not members_csv.exists():
        raise FileNotFoundError(
            f"{members_csv} not found — run pipeline/run_local.py first; this "
            f"module takes its family definition from module 1's output."
        )
    m = pd.read_csv(members_csv)
    return {
        r.gene_id: (r.symbol if isinstance(r.symbol, str) and r.symbol else None)
        for r in m.itertuples()
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


#: Labels come from core so this module cannot disagree with the others.
label_of = core.label_of


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


#: Summed per-library counts are cached here. Reading the five matrices costs a
#: full pass over ~1.2 GB and holds the largest one in memory; the summed table
#: is ~57k x 5 integers. Anything downstream that only needs expression values
#: (the pair statistics, the integration module) can work from the cache, so the
#: matrices are read once per dataset rather than once per question.
COUNTS_CACHE = (
    Path(__file__).resolve().parent.parent / "work" / "expression" / "pseudobulk_counts.csv.gz"
)


def build_pseudobulk(
    datadir: str | Path,
    cache: str | Path | None = COUNTS_CACHE,
    force: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (raw counts, CPM) as genes x libraries, via the summed-count cache."""
    datadir = Path(datadir)
    cache = Path(cache) if cache is not None else None

    if cache is not None and cache.exists() and not force:
        counts = pd.read_csv(cache, index_col=0)
        if list(counts.columns) != list(SAMPLES):
            raise ValueError(
                f"{cache} holds libraries {list(counts.columns)}, expected "
                f"{list(SAMPLES)}; delete it or pass force=True to rebuild."
            )
    else:
        counts = pd.DataFrame(
            {gsm: pseudobulk_library(datadir, gsm, stem)
             for gsm, (stem, _) in SAMPLES.items()}
        )
        if cache is not None:
            cache.parent.mkdir(parents=True, exist_ok=True)
            counts.to_csv(cache)

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


#: Libraries grouped by tissue, derived from SAMPLES so the two stay in step.
NODULE_LIBS = [g for g, (_, t) in SAMPLES.items() if t == "nodule"]
ROOT_LIBS = [g for g, (_, t) in SAMPLES.items() if t == "root"]

#: CPM added to both sides of every ratio before taking a log. Needed because
#: two of the three root libraries are at exactly 0 CPM for every leghemoglobin,
#: so an unregularised fold change is 0/0 there.
#:
#: It is negligible for the four focal genes, whose nodule CPM is 5,000-19,000,
#: and it is **not** negligible for the low-expressed family members: GmLb5 sits
#: at 1.5-2.4 CPM per nodule library and Hb1 below 0.5, so adding 1.0 shifts
#: those by tens of percent and compresses every fold change involving them
#: toward zero. Pairs among the focal four are unaffected; pairs involving
#: GmLb5, Hb1 or Hb2 should be read as regularised, not as measured ratios. If
#: you need those, lower the pseudocount and state the new value — but note that
#: at 15-49 M counts per library, 1 CPM is already a fraction of one count, so
#: the honest reading is that those genes are near the detection floor rather
#: than that the pseudocount is distorting a real measurement.
LOG2FC_PSEUDOCOUNT_CPM = 1.0


def pair_metrics(
    cpm: pd.DataFrame,
    genes: dict[str, str | None],
    pseudocount: float = LOG2FC_PSEUDOCOUNT_CPM,
) -> pd.DataFrame:
    """Pairwise expression statistics over the pseudobulk libraries.

    Emitted per pair, in the canonical ``label_a``/``label_b`` orientation (so
    every fold change is a-over-b, and flipping a row flips the sign):

    ``spearman_profile``
        Spearman rho of the log1p-CPM profiles across all five libraries. Read
        it with care for this family: it is exactly 1.0 for all six pairs among
        the four focal leghemoglobins, because those four share one zero/nonzero
        pattern across the three root libraries and so rank identically. It is
        not 1.0 for every nodule-exclusive pair — GmLb5 is zero in all three root
        libraries where the focal four have one small nonzero value, which gives
        its pairs rho = 0.803 on a different tie-rank pattern. Either way five
        points over two tissues cannot separate co-regulated genes; the column is
        emitted because it was specified, not because it discriminates.
    ``log2fc_mean_all`` / ``log2fc_sd_all``
        Mean and SD of the per-library log2 fold change over all five libraries.
        The SD here is dominated by the nodule-vs-root step rather than by
        variability in the ratio.
    ``log2fc_mean_nodule`` / ``log2fc_sd_nodule``
        The same over the two nodule libraries only. For a nodule-exclusive
        family this is the informative version, and it is what the integration
        module consumes.
    ``dose_ratio``
        Smaller nodule-mean CPM over larger, in [0, 1]. A symmetric measure of
        how comparable the two genes' expression levels are.
    ``cover_a_by_b`` / ``cover_b_by_a``
        Directional: the fraction of one gene's nodule dose that the other could
        supply, capped at 1. Redundancy is not symmetric — a gene at 12% of the
        pool cannot cover the loss of one at 51%, while the reverse holds
        comfortably — and these two columns are what carry that asymmetry.
    ``coexpression_overlap``
        NA. It is defined per cell (cells expressing both over cells expressing
        either) and pseudobulk has no cells; see the module docstring.
    """
    present = [g for g in genes if to_geo_id(g) in cpm.index]
    logcpm = np.log1p(cpm)
    all_libs = list(SAMPLES)
    rows = []
    for a, b in core.canonical_pair_order(present):
        ga, gb = to_geo_id(a), to_geo_id(b)
        va = logcpm.loc[ga, all_libs].to_numpy(dtype=float)
        vb = logcpm.loc[gb, all_libs].to_numpy(dtype=float)
        rho = (np.nan if np.ptp(va) == 0 or np.ptp(vb) == 0
               else float(spearmanr(va, vb).statistic))

        def l2fc(libs: list[str]) -> np.ndarray:
            xa = cpm.loc[ga, libs].to_numpy(dtype=float) + pseudocount
            xb = cpm.loc[gb, libs].to_numpy(dtype=float) + pseudocount
            return np.log2(xa / xb)

        fc_all, fc_nod = l2fc(all_libs), l2fc(NODULE_LIBS)
        na = float(cpm.loc[ga, NODULE_LIBS].mean())
        nb = float(cpm.loc[gb, NODULE_LIBS].mean())
        hi = max(na, nb)

        rows.append({
            "label_a": label_of(a),
            "label_b": label_of(b),
            "gene_a": a,
            "gene_b": b,
            "spearman_profile": rho,
            "n_profile_points": len(all_libs),
            "log2fc_mean_all": float(fc_all.mean()),
            "log2fc_sd_all": float(fc_all.std(ddof=1)),
            "n_libraries_all": len(all_libs),
            "log2fc_mean_nodule": float(fc_nod.mean()),
            "log2fc_sd_nodule": float(fc_nod.std(ddof=1)) if len(fc_nod) > 1 else np.nan,
            "n_libraries_nodule": len(NODULE_LIBS),
            "log2fc_pseudocount_cpm": pseudocount,
            "cpm_mean_nodule_a": na,
            "cpm_mean_nodule_b": nb,
            "dose_ratio": (min(na, nb) / hi) if hi > 0 else np.nan,
            "cover_a_by_b": (min(nb / na, 1.0) if na > 0 else np.nan),
            "cover_b_by_a": (min(na / nb, 1.0) if nb > 0 else np.nan),
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

    family = load_family()
    counts, cpm = build_pseudobulk(datadir)
    prof = gene_profiles(cpm, family)
    pairs = pair_metrics(cpm, family)
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
                            "tau over tissue means", "spearman of log1p-CPM profiles",
                            "log2 fold change mean and SD, all libraries and nodule only",
                            "dose_ratio and directional coverage fractions"],
        "log2fc": {
            "orientation": "a over b, in the canonical label_a/label_b order",
            "pseudocount_cpm": LOG2FC_PSEUDOCOUNT_CPM,
            "pseudocount_reason": (
                "two of three root libraries are at exactly 0 CPM for every "
                "leghemoglobin, so an unregularised ratio is 0/0 there"
            ),
            "pseudocount_scope_caveat": (
                "negligible for the four focal genes (nodule CPM 5,000-19,000) "
                "but not for GmLb5 (1.5-2.4 CPM) or Hb1 (<0.5 CPM): fold changes "
                "involving those are regularised, not measured ratios"
            ),
            "sd_caveat": (
                "log2fc_sd_all is dominated by the nodule-vs-root step, not by "
                "variability in the ratio; log2fc_sd_nodule is over 2 libraries"
            ),
        },
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
