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
import config
import contracts
import soy_globin_core as core  # noqa: E402

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

# The library table, the marker panel, the pseudocount and the tissue names
# all used to be module constants here. They are now in
# `config/config.yaml` under `expression:`, and every function below takes the
# `cfg` it needs. The tissue names in particular were baked into column names
# and into derived constants (NODULE_LIBS / ROOT_LIBS), which is what made this
# module unable to describe any organ pair but nodule-vs-root.


def libraries(cfg) -> list[dict]:
    """``[{gsm, stem, tissue, ...}, ...]`` in config order."""
    return cfg.libraries


def library_ids(cfg) -> list[str]:
    return cfg.library_ids


def libs_for_tissue(cfg, tissue: str) -> list[str]:
    """Replaces the NODULE_LIBS / ROOT_LIBS module constants."""
    return cfg.libraries_for_tissue(tissue)

def load_family(members_csv: str | Path) -> dict[str, str | None]:
    """Family gene IDs -> symbol, read from the sequence module's member table.

    The family is *discovered* by the HMM search in module 1, not declared, so
    this module reads it rather than keeping a copy. It used to keep one, and
    the copy disagreed with the shared symbol map about whether
    ``Glyma.10G198900`` carries a symbol — which silently cost 6 of 21 rows on
    any join between this module's pair table and the sequence module's.

    ``members_csv`` is now a required argument. It was a module-level constant
    built from ``Path(__file__).parent.parent``, which meant the dependency was
    real but undeclared: it could not be redirected, overridden in a test, or
    checked for staleness (docs/PIPELINE_AUDIT.md finding 3).
    """
    members_csv = Path(members_csv)
    if not members_csv.exists():
        raise FileNotFoundError(
            f"{members_csv} not found — the sequence stage has not run. This "
            f"module takes its family definition from module 1's output."
        )
    m = pd.read_csv(members_csv)
    return {
        r.gene_id: (r.symbol if isinstance(r.symbol, str) and r.symbol else None)
        for r in m.itertuples()
    }


def markers(cfg) -> dict[str, str]:
    """Gene ID -> role string for the configured marker panel.

    Independent cell-type markers, reported as context for the pseudobulk
    profiles and deliberately not used to define them: using leghemoglobin to
    mark the infected-cell compartment would be circular in a study of
    leghemoglobin.
    """
    out = {}
    for gene, spec in (cfg.get("expression.markers") or {}).items():
        role = spec.get("role", "") if isinstance(spec, dict) else str(spec)
        up = spec.get("uniprot") if isinstance(spec, dict) else None
        out[gene] = f"{role} ({up})" if up else role
    return out


def part_names(cfg) -> tuple[str, ...]:
    """The CellRanger parts each library ships."""
    return tuple(cfg["expression.series.parts"])


def expected_files(cfg) -> list[str]:
    """Canonical filenames for every library part, in config order.

    The filenames are exactly GEO's own ``{GSM}_{stem}_{part}`` — there is no
    local renaming step, so a workflow rule can name these as outputs directly.
    """
    return [f"{lib['gsm']}_{lib['stem']}_{part}"
            for lib in libraries(cfg) for part in part_names(cfg)]


def fetch_library(cfg, datadir: str | Path, gsm: str, force: bool = False) -> dict:
    """Download the three parts of one library. Returns a provenance record.

    One library per call, so each is its own workflow rule: a re-run fetches
    only what is missing rather than re-walking all ~1.2 GB.
    """
    import urllib.request

    lib = next((l for l in libraries(cfg) if l["gsm"] == gsm), None)
    if lib is None:
        raise KeyError(f"{gsm} is not in expression.libraries")
    datadir = Path(datadir)
    datadir.mkdir(parents=True, exist_ok=True)
    base = cfg["expression.series.ftp_base"].format(dir=f"{gsm[:-3]}nnn", gsm=gsm)
    got = {}
    for part in part_names(cfg):
        name = f"{gsm}_{lib['stem']}_{part}"
        dest = datadir / name
        if force or not dest.exists() or dest.stat().st_size == 0:
            print(f"[fetch] {base}/{name}", flush=True)
            urllib.request.urlretrieve(f"{base}/{name}", dest)
        got[part] = {"path": str(dest), "bytes": dest.stat().st_size,
                     "sha256": core.sha256(dest)}
    return {"gsm": gsm, "stem": lib["stem"], "tissue": lib["tissue"], "parts": got}


def fetch_inputs(cfg, datadir: str | Path, force: bool = False) -> dict[str, dict]:
    """Download every configured library."""
    return {lib["gsm"]: fetch_library(cfg, datadir, lib["gsm"], force=force)
            for lib in libraries(cfg)}


def check_inputs(cfg, datadir: str | Path) -> list[str]:
    """Return the canonical filenames that are missing from ``datadir``."""
    datadir = Path(datadir)
    return [n for n in expected_files(cfg) if not (datadir / n).exists()]


def input_checksums(cfg, datadir: str | Path) -> dict[str, str]:
    """sha256 of every matrix and feature file the pseudobulk sum reads.

    This is the pseudobulk cache key. The cache used to be validated against
    the *library list* only, so replacing a matrix file on disk left a stale
    summed table in place and silently reused it (audit 3.2). Barcodes are
    excluded: they are downloaded for completeness but the sum does not read
    them, so their digest would invalidate the cache without affecting it.
    """
    datadir = Path(datadir)
    keep = [p for p in part_names(cfg) if "barcodes" not in p]
    return {n: core.sha256(datadir / n)
            for n in expected_files(cfg)
            if any(n.endswith(p) for p in keep)}


#: GSE226149 features.tsv uses GLYMA_10G199100; everything else uses Glyma.10G199100.
#: Gene-ID transforms for matching the family against a series feature list.
#: GSE226149's ``features.tsv.gz`` writes ``GLYMA_10G199100``, not
#: ``Glyma.10G199100``; any lookup that skips this silently matches nothing
#: (README §9.16). Named in config as ``family.identifiers.geo_id_transform``
#: so a series with a different convention is a config change.
ID_TRANSFORMS = {
    "upper_underscore": lambda g: g.replace(".", "_").upper(),
    "identity": lambda g: g,
}


def to_geo_id(gene_id: str, transform: str = "upper_underscore") -> str:
    try:
        return ID_TRANSFORMS[transform](gene_id)
    except KeyError:
        raise ValueError(
            f"unknown geo_id_transform {transform!r}; known: {sorted(ID_TRANSFORMS)}"
        ) from None


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


def build_pseudobulk(
    cfg,
    datadir: str | Path,
    cache: str | Path | None = None,
    force: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (raw counts, CPM) as genes x libraries, via the summed-count cache.

    Reading the matrices costs a full pass over ~1.2 GB and holds the largest
    one in memory; the summed table is ~57k x 5 integers. Everything downstream
    needs only the sums, so the matrices are read once per dataset rather than
    once per question. Peak memory is one matrix, not the whole atlas — which
    is why this stage does not need remote dispatch.

    **The cache is keyed on the sha256 of the input matrices**, recorded in a
    sidecar ``.key.json``. Previously it was validated against the library
    *list* only, so swapping a matrix file left the stale summed table in place
    and reused it without a word (audit 3.2). A key mismatch rebuilds rather
    than raising: the inputs are the authority, and a cache is only ever an
    optimisation.
    """
    datadir = Path(datadir)
    cache = Path(cache) if cache is not None else Path(cfg["expression.pseudobulk.cache"])
    keyfile = cache.with_suffix(".key.json") if cache.suffix != ".json" else None
    want_key = {"libraries": library_ids(cfg),
                "inputs": input_checksums(cfg, datadir)}

    counts = None
    if cache.exists() and not force:
        have_key = None
        if keyfile is not None and keyfile.exists():
            try:
                have_key = json.loads(keyfile.read_text())
            except json.JSONDecodeError:
                have_key = None
        if have_key == want_key:
            counts = pd.read_csv(cache, index_col=0)
        else:
            why = ("no cache key recorded (cache predates checksum keying)"
                   if have_key is None else
                   "input checksums or library list have changed")
            print(f"[pseudobulk] rebuilding {cache.name}: {why}", flush=True)

    if counts is None:
        counts = pd.DataFrame(
            {lib["gsm"]: pseudobulk_library(datadir, lib["gsm"], lib["stem"])
             for lib in libraries(cfg)}
        )
        cache.parent.mkdir(parents=True, exist_ok=True)
        counts.to_csv(cache)
        if keyfile is not None:
            keyfile.write_text(json.dumps(want_key, indent=2, sort_keys=True))

    if list(counts.columns) != library_ids(cfg):
        raise ValueError(
            f"{cache} holds libraries {list(counts.columns)}, config declares "
            f"{library_ids(cfg)}"
        )
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


def gene_profiles(cfg, cpm: pd.DataFrame, genes: dict[str, str | None]) -> pd.DataFrame:
    """Per-gene pseudobulk profile, tissue means, and tau over tissue means.

    Tissue order follows the config's library order rather than alphabetical
    sort, so the columns read focal-tissue-first and match the manifest.

    ``tau_over_tissue_means`` is named for what it is: with two tissues it
    collapses to a focal-vs-contrast contrast and is not a tissue-specificity
    index in the usual sense (README §9.9). ``detection_rate`` is NA with a
    reason rather than replaced by a look-alike statistic — it is a per-cell
    quantity and pseudobulk has no cells.
    """
    tissues = cfg.tissue_names
    transform = cfg["family.identifiers.geo_id_transform"]
    symbols = cfg.gene_symbols
    notcomp = (cfg.get("expression.metrics.not_computable") or {}).get(
        "detection_rate", "per-cell quantity; pseudobulk has no cells")
    libs_by_tissue = {t: cfg.libraries_for_tissue(t) for t in tissues}
    rows = []
    for gene in genes:
        gid = to_geo_id(gene, transform)
        if gid not in cpm.index:
            rows.append({"gene_id": gene, "label": core.label_of(gene, symbols),
                         "in_matrix": False})
            continue
        v = cpm.loc[gid]
        tmeans = {t: float(v[libs_by_tissue[t]].mean()) for t in tissues}
        rows.append({
            "gene_id": gene,
            "label": core.label_of(gene, symbols),
            "symbol": genes[gene] or "",
            "in_matrix": True,
            **{f"cpm_{g}": float(v[g]) for g in library_ids(cfg)},
            **{cfg.tissue_mean_column(t): tmeans[t] for t in tissues},
            "tau_over_tissue_means": tau(np.array([tmeans[t] for t in tissues])),
            "n_tissues_for_tau": len(tissues),
            "detection_rate": np.nan,
            "detection_rate_note": notcomp,
        })
    return pd.DataFrame(rows)


#: CPM added to both sides of every ratio before taking a log. Configured as
#: ``expression.metrics.log2fc_pseudocount_cpm``; the value and this reasoning
#: travel together into the manifest. Needed because
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
PSEUDOCOUNT_REASON = (
    "two of three contrast-tissue libraries are at exactly 0 CPM for every "
    "family member, so an unregularised ratio is 0/0 there"
)
PSEUDOCOUNT_SCOPE_CAVEAT = (
    "negligible for the high-expressed focal genes but not for members near "
    "the detection floor: fold changes involving those are regularised, not "
    "measured ratios"
)


def pair_metrics(
    cfg,
    cpm: pd.DataFrame,
    genes: dict[str, str | None],
    pseudocount: float | None = None,
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
    if pseudocount is None:
        pseudocount = float(cfg["expression.metrics.log2fc_pseudocount_cpm"])
    transform = cfg["family.identifiers.geo_id_transform"]
    symbols = cfg.gene_symbols
    focal_tissue = cfg.tissues["focal"]
    focal_libs = cfg.libraries_for_tissue(focal_tissue)
    notcomp = (cfg.get("expression.metrics.not_computable") or {}).get(
        "coexpression_overlap", "per-cell quantity; pseudobulk has no cells")

    present = [g for g in genes if to_geo_id(g, transform) in cpm.index]
    logcpm = np.log1p(cpm)
    all_libs = library_ids(cfg)
    rows = []
    for a, b in core.canonical_pair_order(present):
        ga, gb = to_geo_id(a, transform), to_geo_id(b, transform)
        va = logcpm.loc[ga, all_libs].to_numpy(dtype=float)
        vb = logcpm.loc[gb, all_libs].to_numpy(dtype=float)
        rho = (np.nan if np.ptp(va) == 0 or np.ptp(vb) == 0
               else float(spearmanr(va, vb).statistic))

        def l2fc(libs: list[str]) -> np.ndarray:
            xa = cpm.loc[ga, libs].to_numpy(dtype=float) + pseudocount
            xb = cpm.loc[gb, libs].to_numpy(dtype=float) + pseudocount
            return np.log2(xa / xb)

        fc_all, fc_foc = l2fc(all_libs), l2fc(focal_libs)
        na = float(cpm.loc[ga, focal_libs].mean())
        nb = float(cpm.loc[gb, focal_libs].mean())
        hi = max(na, nb)

        rows.append({
            "label_a": core.label_of(a, symbols),
            "label_b": core.label_of(b, symbols),
            "gene_a": a,
            "gene_b": b,
            "spearman_profile": rho,
            "n_profile_points": len(all_libs),
            "log2fc_mean_all": float(fc_all.mean()),
            "log2fc_sd_all": float(fc_all.std(ddof=1)),
            "n_libraries_all": len(all_libs),
            f"log2fc_mean_{focal_tissue}": float(fc_foc.mean()),
            f"log2fc_sd_{focal_tissue}": (float(fc_foc.std(ddof=1))
                                          if len(fc_foc) > 1 else np.nan),
            f"n_libraries_{focal_tissue}": len(focal_libs),
            "log2fc_pseudocount_cpm": pseudocount,
            f"cpm_mean_{focal_tissue}_a": na,
            f"cpm_mean_{focal_tissue}_b": nb,
            "dose_ratio": (min(na, nb) / hi) if hi > 0 else np.nan,
            "cover_a_by_b": (min(nb / na, 1.0) if na > 0 else np.nan),
            "cover_b_by_a": (min(na / nb, 1.0) if nb > 0 else np.nan),
            "coexpression_overlap": np.nan,
            "coexpression_overlap_note": notcomp,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def run(
    cfg,
    datadir: str | Path,
    outdir: str | Path,
    members_csv: str | Path,
    fetch: bool = True,
) -> dict:
    """The whole expression stage. Every path is an argument."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    missing = check_inputs(cfg, datadir)
    if missing:
        if not fetch:
            raise FileNotFoundError(
                f"{len(missing)} input file(s) missing from {datadir}, "
                f"first: {missing[0]}. Run with fetch=True."
            )
        fetch_inputs(cfg, datadir)
        still = check_inputs(cfg, datadir)
        if still:
            raise FileNotFoundError(f"fetch did not produce: {still}")

    family = load_family(members_csv)
    counts, cpm = build_pseudobulk(cfg, datadir)
    prof = gene_profiles(cfg, cpm, family)
    pairs = pair_metrics(cfg, cpm, family)

    mk = markers(cfg)
    marks = gene_profiles(cfg, cpm, {g: None for g in mk})
    marks["marker_role"] = [mk[g] for g in marks["gene_id"]]

    n_members = len(family)
    contracts.write_csv(prof, outdir / "gene_pseudobulk_profiles.csv",
                        contracts.gene_profiles_spec(cfg, n_members))
    contracts.write_csv(pairs, outdir / "expression_pairs.csv",
                        contracts.expression_pairs_spec(cfg, n_members))
    contracts.write_csv(marks, outdir / "marker_profiles.csv",
                        contracts.gene_profiles_spec(cfg, len(mk), marker=True))

    checks: dict = {}
    focal_tissue = cfg.tissues["focal"]
    for lib in libraries(cfg):
        want = lib.get("total_counts_expected")
        got = int(counts[lib["gsm"]].sum())
        if want is not None:
            config.evaluate_check(
                cfg, "library_total_counts", got == int(want),
                {"library": lib["gsm"], "expected": int(want), "observed": got},
                sink=checks.setdefault("library_total_counts", {}),
            )
    absent = marks.loc[~marks.in_matrix.astype(bool), "gene_id"].tolist()

    manifest = {
        "config_digest": cfg.digest(),
        "resolution": cfg["expression.resolution"],
        "source_series": cfg["expression.series.accession"],
        "source_note": cfg["expression.series.note"],
        "rejected_series": cfg.get("expression.rejected_series.accession"),
        "rejected_reason": cfg.get("expression.rejected_series.reason"),
        "libraries": {lib["gsm"]: {"stem": lib["stem"], "tissue": lib["tissue"],
                                   "total_counts": int(counts[lib["gsm"]].sum())}
                      for lib in libraries(cfg)},
        "tissue_roles": cfg.tissues,
        "input_sha256": input_checksums(cfg, datadir),
        "normalisation": f"{cfg['expression.pseudobulk.normalisation'].upper()} over "
                         f"summed library counts; no cell calling",
        "metrics_emitted": ["cpm per library", "cpm mean per tissue",
                            "tau over tissue means", "spearman of log1p-CPM profiles",
                            "log2 fold change mean and SD, all libraries and "
                            f"{focal_tissue} only",
                            "dose_ratio and directional coverage fractions"],
        "log2fc": {
            "orientation": cfg["expression.metrics.log2fc_orientation"],
            "pseudocount_cpm": float(cfg["expression.metrics.log2fc_pseudocount_cpm"]),
            "pseudocount_reason": PSEUDOCOUNT_REASON,
            "pseudocount_scope_caveat": PSEUDOCOUNT_SCOPE_CAVEAT,
            "sd_caveat": (f"log2fc_sd_all is dominated by the inter-tissue step, not "
                          f"by variability in the ratio; log2fc_sd_{focal_tissue} is "
                          f"over {len(cfg.libraries_for_tissue(focal_tissue))} "
                          f"libraries"),
        },
        "metrics_degenerate": cfg.get("expression.metrics.degenerate_columns") or {},
        "metrics_not_computable": cfg.get("expression.metrics.not_computable") or {},
        "tau_caveat": cfg.get("expression.metrics.tau_caveat"),
        "markers_absent_from_series": absent,
        "validation_checks": checks,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    return {"counts": counts, "cpm": cpm, "profiles": prof, "pairs": pairs,
            "markers": marks, "manifest": manifest}
