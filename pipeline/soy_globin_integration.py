"""Module 4 — integration and the redundancy score.

Joins the three pair tables and reduces them to one score per pair. Five
decisions are recorded here because none of them is visible in the output:

* **The score is gated, not additive.** ``R = M**alpha * E**beta``, not
  ``alpha*M + beta*E``. Two genes that are 95% identical but never expressed in
  the same place are not functionally redundant — neither can buffer the other's
  loss — and an additive score awards that pair a high value on the strength of
  sequence alone. Redundancy is a conjunction, so the score multiplies.

* **M uses three sub-axes, not four flat features.** Across all 21 pairs
  ``pid_aligned``, ESM2 similarity, ``tm_score`` and ``pocket_identity`` are
  collinear (Spearman 0.85-0.94), so averaging all four as equals is really just
  a re-weighting of one axis, with sequence counted twice. They come apart
  *within* the Lb clade (identity vs pocket identity is only rho = 0.39), which
  is where the question lives. So: sequence (identity and ESM2 averaged), global
  fold (TM-score), and heme pocket, a third each.

* **Two normalisation scopes are emitted, because the choice changes the answer
  more than the weights do.** Min-max over all 21 pairs leaves the six Lb pairs
  compressed into the top few percent of every molecular axis; min-max over the
  six Lb pairs alone spreads them across the full range. ``R_family`` is the
  absolute score, comparable across the family, on which outgroup pairs
  correctly sit near zero. ``R_clade`` is a *relative* ranking inside the Lb
  clade, where 0 means "lowest of the six observed", not "not redundant". The
  gap between them is the saturation result, not a defect to be hidden by
  picking one.

* **E is expression at tissue resolution, not the per-cell overlap the design
  asked for.** Cell-level co-expression overlap needs cell-type labels the
  expression module does not have, and the available stand-in — Spearman of the
  five-library profiles — is exactly 1.000 for all six Lb pairs and therefore
  carries no information. E is instead ``tissue_overlap * dose_ratio``: do the
  two genes put their transcript in the same tissue, and at comparable levels.
  Dose is load-bearing rather than decorative — a gene supplying 12% of the
  leghemoglobin pool cannot cover the loss of one supplying 51%.

* **Redundancy is directional, and the symmetric score loses that.** The
  directional variant asks "how much of i's dose could j supply", which is
  capped at 1 in the easy direction and small in the hard one. Lba's pairs
  therefore rank low symmetrically while covering their partners well
  one-directionally. That inverts the sanity check in the original plan, which
  expected Lba's pairs at the top of the ranking; see ``validation_checks`` in
  the manifest for what replaced it.

No Modal import, for the same reason ``soy_globin_core`` has none.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config                            # noqa: E402
import contracts                         # noqa: E402

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

# The weights, the sub-axis weights, the clade floor and the tissue column
# names are all in `config/config.yaml` under `score:` and `expression:`.
#
# alpha was previously a module global that `run_integration.py` *reassigned*
# from the outside before calling run() — the only parameter in the repo set by
# mutating another module's state (audit 4.3). It is now threaded through as an
# argument. The reasoning behind the value travels with it in the config:
# beta > alpha leans toward expression for a measured reason rather than an
# expected one, because the molecular axes are saturated within the focal clade
# (91.7-95.2% identity, TM 0.9816-0.9976) while dose ratio spans 0.26-0.79
# across the same six pairs.

#: Pair-class labels. Derived from focal membership, not from gene symbols:
#: `focal-focal` was `"Lb-Lb"`, which put the family name in the data.
CLASS_FOCAL_FOCAL = "focal-focal"
CLASS_FOCAL_OTHER = "focal-other"
CLASS_OTHER_OTHER = "other-other"


def tissue_columns(cfg) -> list[str]:
    """Tissue-mean CPM columns in the expression module's profile table."""
    return [cfg.tissue_mean_column(t) for t in cfg.tissue_names]


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #


def load_pairs(results_dir: str | Path, cfg=None, n_members: int | None = None
               ) -> pd.DataFrame:
    """Join the three module pair tables on the canonical pair key.

    Raises rather than returning a short table: a silent row loss here is
    exactly the failure this pipeline already had once, when two modules
    disagreed about whether ``Glyma.10G198900`` carries a symbol.

    When ``cfg`` is given, each input is additionally read through its declared
    contract first, so a malformed table is rejected by name at its own
    boundary rather than surfacing here as an unexplained unjoined row. This
    function's own orientation and completeness checks stay: they are the
    cross-table invariants no single-table contract can see.
    """
    results_dir = Path(results_dir)
    if cfg is not None:
        seq = contracts.read_csv(results_dir / "sequence_module" / "paralog_pairs.csv",
                                 contracts.paralog_pairs_spec(cfg, n_members))
        struc = contracts.read_csv(results_dir / "structure_module" / "structure_pairs.csv",
                                   contracts.structure_pairs_spec(cfg, n_members))
        expr = contracts.read_csv(results_dir / "expression_module" / "expression_pairs.csv",
                                  contracts.expression_pairs_spec(cfg, n_members))
    else:
        seq = pd.read_csv(results_dir / "sequence_module" / "paralog_pairs.csv")
        struc = pd.read_csv(results_dir / "structure_module" / "structure_pairs.csv")
        expr = pd.read_csv(results_dir / "expression_module" / "expression_pairs.csv")

    key = ["label_a", "label_b"]
    for name, df in (("sequence", seq), ("structure", struc), ("expression", expr)):
        bad = [
            f"{a}|{b}" for a, b in zip(df.label_a, df.label_b)
            if a.split("_")[0] >= b.split("_")[0]
        ]
        if bad:
            raise ValueError(
                f"{name} pair table is not in canonical orientation ({len(bad)} "
                f"rows, first {bad[0]}); re-run that module — every emitter is "
                f"supposed to iterate core.canonical_pair_order()."
            )

    pairs = seq.merge(
        struc.drop(columns=["gene_a", "gene_b"]), on=key, how="outer",
        indicator="_m_struc",
    ).merge(
        expr.drop(columns=["gene_a", "gene_b"]), on=key, how="outer",
        indicator="_m_expr",
    )
    unmatched = pairs[(pairs._m_struc != "both") | (pairs._m_expr != "both")]
    if len(unmatched):
        raise ValueError(
            f"{len(unmatched)} of {len(pairs)} pairs did not join across all "
            f"three modules, e.g. {unmatched.iloc[0][key].tolist()}"
        )
    return pairs.drop(columns=["_m_struc", "_m_expr"])


def load_gene_profiles(results_dir: str | Path) -> pd.DataFrame:
    """Per-gene tissue-mean CPM, indexed by gene ID."""
    p = Path(results_dir) / "expression_module" / "gene_pseudobulk_profiles.csv"
    return pd.read_csv(p).set_index("gene_id")


# --------------------------------------------------------------------------- #
# Components
# --------------------------------------------------------------------------- #


def minmax(x: pd.Series, lo: float | None = None, hi: float | None = None) -> pd.Series:
    """Scale to [0, 1] against an explicit or observed range.

    A degenerate range (every value equal) maps to 1.0 rather than NaN: that is
    the honest reading of "this feature does not distinguish these pairs".
    """
    lo = float(x.min()) if lo is None else lo
    hi = float(x.max()) if hi is None else hi
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return pd.Series(1.0, index=x.index)
    return ((x - lo) / (hi - lo)).clip(0.0, 1.0)


def tissue_overlap(profiles: pd.DataFrame, gene_a: str, gene_b: str,
                   tissue_cols: list[str]) -> float:
    """Histogram intersection of two genes' normalised tissue profiles, in [0, 1].

    Each gene's tissue-mean CPM is turned into a distribution over tissues and
    the two distributions are intersected: 1.0 means the two genes place their
    transcript in the same tissues in the same proportions, 0.0 means disjoint
    tissues. This replaces the binary co-expression gate of the original design —
    with only two tissues a hard gate is either always open or always shut —
    while keeping its logic, that co-location is a precondition for buffering.
    """
    out = []
    for g in (gene_a, gene_b):
        v = np.array([float(profiles.loc[g, c]) for c in tissue_cols])
        tot = v.sum()
        out.append(v / tot if tot > 0 else np.zeros_like(v))
    return float(np.minimum(out[0], out[1]).sum())


def add_components(pairs: pd.DataFrame, profiles: pd.DataFrame, cfg) -> pd.DataFrame:
    """Attach the raw (un-normalised) score components to every pair."""
    df = pairs.copy()
    df["esm2_similarity"] = -df["esm2_cosine_distance"]
    df["tissue_overlap"] = [
        tissue_overlap(profiles, a, b, tissue_columns(cfg))
        for a, b in zip(df.gene_a, df.gene_b)
    ]
    # Focal membership comes from the sequence module's both_focal flag rather
    # than from symbol presence — every member now carries a symbol, so symbols
    # no longer distinguish the focal four.
    focal = set(df.loc[df.both_focal, "gene_a"]) | set(df.loc[df.both_focal, "gene_b"])
    n_focal = [int(a in focal) + int(b in focal) for a, b in zip(df.gene_a, df.gene_b)]
    df["pair_class"] = [CLASS_FOCAL_FOCAL if n == 2
                        else CLASS_FOCAL_OTHER if n == 1 else CLASS_OTHER_OTHER
                        for n in n_focal]
    return df


#: Floor for the clade rescale. Plain min-max would put the clade minimum of M or
#: of E at exactly 0, and since R multiplies the two, any pair that is lowest on
#: *either* factor would score 0 — collapsing distinct pairs onto the same value
#: and destroying the ordering the clade score exists to display. The floor is a
#: presentation choice on a relative score, and it changes no ordering.
#: Floor applied when rescaling within the focal clade. A plain min-max zero in
#: either factor would zero the product and merge pairs that differ; the floor
#: is a presentation choice on a relative score and changes no ordering.
#: Configured as ``score.normalisation.clade_floor``.
DEFAULT_CLADE_FLOOR = 0.05


def score(
    pairs: pd.DataFrame,
    scope: str = "family",
    alpha: float = 0.40,
    beta: float = 0.60,
    clade_floor: float = DEFAULT_CLADE_FLOOR,
    m_axis_weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Compute M, E and R under one normalisation scope.

    The molecular axes are always min-max scaled over the *whole family*, so M
    sits on one absolute scale regardless of scope. ``scope='family'`` uses M and
    E as they stand; R is then absolute, and outgroup pairs sit near zero.
    ``scope='clade'`` rescales M and E — each as a whole, not factor by factor —
    across the six focal pairs only, so R becomes a relative ranking inside the
    clade in which the lowest value means "lowest of the six observed", not "not
    redundant".
    """
    if scope not in ("family", "clade"):
        raise ValueError(f"scope must be 'family' or 'clade', got {scope!r}")
    df = pairs.copy()
    lb = df[df.pair_class == CLASS_FOCAL_FOCAL]
    if lb.empty:
        raise ValueError(f"no {CLASS_FOCAL_FOCAL} pairs to normalise against")

    # Molecular axes on the absolute (family) scale.
    seq = 0.5 * (minmax(df["pid_aligned"]) + minmax(df["esm2_similarity"]))
    fold = minmax(df["tm_score"])
    pocket = minmax(df["pocket_identity"])
    w = m_axis_weights or {"seq": 1 / 3, "fold": 1 / 3, "pocket": 1 / 3}
    M = w["seq"] * seq + w["fold"] * fold + w["pocket"] * pocket

    # E and its two directional variants. Both factors are already bounded in
    # [0, 1] with a meaningful zero, so nothing is normalised here.
    overlap = df["tissue_overlap"]
    E = overlap * df["dose_ratio"]
    E_ab = overlap * df["cover_a_by_b"]
    E_ba = overlap * df["cover_b_by_a"]

    if scope == "clade":
        def rescale(x: pd.Series) -> pd.Series:
            u = minmax(x, float(x[df.pair_class == CLASS_FOCAL_FOCAL].min()),
                       float(x[df.pair_class == CLASS_FOCAL_FOCAL].max()))
            return clade_floor + (1.0 - clade_floor) * u

        M, E, E_ab, E_ba = rescale(M), rescale(E), rescale(E_ab), rescale(E_ba)

    sfx = f"_{scope}"
    df[f"M_seq{sfx}"], df[f"M_fold{sfx}"], df[f"M_pocket{sfx}"] = seq, fold, pocket
    df[f"M{sfx}"], df[f"E{sfx}"] = M, E
    df[f"R{sfx}"] = M.pow(alpha) * E.pow(beta)
    # Directional: how much of one member's dose the other could supply. The
    # molecular term is symmetric; only the suppliable dose is not.
    df[f"R_a_covered_by_b{sfx}"] = M.pow(alpha) * E_ab.pow(beta)
    df[f"R_b_covered_by_a{sfx}"] = M.pow(alpha) * E_ba.pow(beta)
    return df


def weight_sensitivity(
    pairs: pd.DataFrame,
    scope: str = "family",
    steps: int = 21,
    alpha_min: float = 0.0,
    alpha_max: float = 1.0,
    clade_floor: float = DEFAULT_CLADE_FLOOR,
    m_axis_weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Sweep alpha from 0 to 1 and record what the ranking does.

    Reported rather than asserted: if the top pair is the same at every alpha the
    weighting is not doing the work, and if it flips the score should not be
    presented without this table beside it.
    """
    rows = []
    for a in np.linspace(alpha_min, alpha_max, steps):
        s = score(pairs, scope=scope, alpha=float(a), beta=float(1.0 - a),
                  clade_floor=clade_floor, m_axis_weights=m_axis_weights)
        foc = s[s.pair_class == CLASS_FOCAL_FOCAL].sort_values(
            f"R_{scope}", ascending=False)
        top_all = s.sort_values(f"R_{scope}", ascending=False).iloc[0]
        rows.append({
            "alpha": round(float(a), 3),
            "beta": round(1.0 - float(a), 3),
            "top_pair_overall": f"{top_all.label_a}|{top_all.label_b}",
            "top_focal_pair": f"{foc.iloc[0].label_a}|{foc.iloc[0].label_b}",
            "top_focal_R": float(foc.iloc[0][f"R_{scope}"]),
            "focal_R_spread": float(foc[f"R_{scope}"].max() - foc[f"R_{scope}"].min()),
        })
    return pd.DataFrame(rows)


def _family_composition(out: pd.DataFrame) -> dict:
    """Record whether this family had an outgroup, and what that actually cost.

    When every member is focal there is no outgroup, and one consequence is
    certain: ``focal_outgroup_separation`` has nothing to compare and does not
    run. A second consequence is tempting to assume and is **false**, so it is
    measured here rather than asserted.

    The intuition is that with no outgroup the clade scope covers the same
    pairs as the family scope, so ``R_clade`` must collapse onto ``R_family``.
    The pair *set* does coincide — but ``scope='clade'`` still rescales M and E
    *separately, each as a whole*, onto ``[clade_floor, 1]``. Rescaling the two
    factors independently changes their relative contribution to the product,
    so the two columns differ in value and can differ in **ordering**: on the
    committed family run without a focal subset, Spearman between them is 0.64
    and the rank order is not preserved. ``R_clade`` therefore stays
    informative even with no outgroup, and describing it as redundant would
    have told a reader to ignore a column that still says something.

    What is lost with no outgroup is the *interpretation*, not the content:
    ``R_clade`` can no longer be read as "relative within the ingroup, against
    a family that is larger than it".
    """
    n_focal = int((out.pair_class == CLASS_FOCAL_FOCAL).sum())
    n_other = int(len(out) - n_focal)
    rec = {
        "n_focal_pairs": n_focal,
        "n_outgroup_pairs": n_other,
        "has_outgroup": bool(n_other > 0),
        "pair_classes_present": sorted(out.pair_class.unique().tolist()),
    }
    if {"R_clade", "R_family"} <= set(out.columns) and len(out) > 2:
        fam = out.R_family.to_numpy(float)
        cla = out.R_clade.to_numpy(float)
        rec["clade_vs_family"] = {
            "scope_covers_same_pairs": bool(n_other == 0),
            "max_abs_difference": float(np.nanmax(np.abs(cla - fam))),
            "spearman": float(pd.Series(cla).corr(pd.Series(fam), method="spearman")),
            "same_ordering": bool(
                (pd.Series(cla).rank().to_numpy()
                 == pd.Series(fam).rank().to_numpy()).all()),
        }
    rec["note"] = (
        "no outgroup: every discovered member is focal, so "
        "focal_outgroup_separation did not run. R_clade still applies its own "
        "rescaling and is NOT a duplicate of R_family — see clade_vs_family — "
        "but it can no longer be read as 'relative within the ingroup'."
        if n_other == 0 else
        "focal and outgroup pairs both present; both normalisation scopes "
        "carry their intended meaning."
    )
    return rec


def validation_checks(scored: pd.DataFrame, profiles: pd.DataFrame, cfg) -> dict:
    """Checks on the family score, each recorded with its evidence.

    The original plan's check — "the dominant gene's pairs should rank top" —
    is not usable against a dose-aware score, which by construction ranks a
    dominant gene's pairs low. These are the checks that survive that change.

    Each is now evaluated under its configured policy rather than printed and
    ignored: the runner used to print `passed=` and exit 0 either way, which
    made this the only place in the pipeline where a scientific assertion was
    evaluated at all, and also the only place where failing one cost nothing
    (audit 4.5).
    """
    focal_tissue = cfg.tissues["focal"]
    abundance_col = cfg.tissue_mean_column(focal_tissue)
    fam = scored.sort_values("R_family", ascending=False).reset_index(drop=True)
    foc = fam[fam.pair_class == CLASS_FOCAL_FOCAL]
    other = fam[fam.pair_class != CLASS_FOCAL_FOCAL]
    sink: dict = {}

    # 1. Do the focal pairs separate from everything else?
    #
    #    When the family was defined with no focal subset, every pair is
    #    focal-focal and `other` is EMPTY. There is then nothing to separate
    #    from, and the check is not applicable — which is emphatically not the
    #    same as passing. Reporting passed=True on an empty comparison would
    #    turn this policy=fail gate into a silent green light: the run would go
    #    green on a check that never ran. The configured policy is left alone,
    #    so the same config re-run on a family that does have an outgroup
    #    re-arms the gate with no edit.
    have_outgroup = len(other) > 0
    config.evaluate_check(
        cfg, "focal_outgroup_separation",
        bool(have_outgroup and foc.R_family.min() > other.R_family.max()),
        {"focal_min_R": float(foc.R_family.min()),
         "other_max_R": float(other.R_family.max()) if have_outgroup else None,
         "n_focal_pairs": int(len(foc)), "n_other_pairs": int(len(other))},
        sink=sink,
        applicable=have_outgroup,
        reason=None if have_outgroup else (
            "no outgroup pairs: the family was defined without a focal subset, "
            "so every discovered member is focal and there is nothing to "
            "separate from. Declare family.<source>.focal to re-arm this check."
        ),
    )

    # 2. Does directional coverage run in the direction of measured abundance?
    #    For every focal pair, the lower-expressed member should be the one that
    #    is better covered. This is the check that would catch the directional
    #    scores being wired backwards, which no amount of internal consistency
    #    would reveal.
    agree, total = 0, 0
    for r in foc.itertuples():
        ca = float(getattr(r, "R_a_covered_by_b_family"))
        cb = float(getattr(r, "R_b_covered_by_a_family"))
        na = float(profiles.loc[r.gene_a, abundance_col])
        nb = float(profiles.loc[r.gene_b, abundance_col])
        if ca == cb or na == nb:
            continue
        total += 1
        agree += int((ca > cb) == (na < nb))
    config.evaluate_check(
        cfg, "directional_coverage_abundance",
        bool(total > 0 and agree == total),
        {"pairs_tested": total, "agreeing": agree,
         "abundance_column": abundance_col},
        sink=sink,
    )

    # 3. Which focal pair tops the symmetric ranking, and by how much? Recorded,
    #    never asserted: the top two are separated by less than the spread the
    #    alpha sweep produces, so the identity of the single most redundant pair
    #    is not a claim this pipeline can make.
    top = foc.iloc[0]
    config.evaluate_check(
        cfg, "symmetric_ranking_top", True,
        {"top_focal_pair": f"{top.label_a}|{top.label_b}",
         "top_focal_R_family": float(top.R_family),
         "runner_up": f"{foc.iloc[1].label_a}|{foc.iloc[1].label_b}",
         "runner_up_R_family": float(foc.iloc[1].R_family),
         "margin": float(top.R_family - foc.iloc[1].R_family)},
        sink=sink,
    )
    return sink


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def run(
    cfg,
    results_dir: str | Path,
    outdir: str | Path,
    alpha: float | None = None,
) -> dict:
    """The whole integration stage. ``alpha`` overrides ``score.alpha_m``."""
    results_dir, outdir = Path(results_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    alpha = cfg.alpha if alpha is None else float(alpha)
    beta = 1.0 - alpha
    floor = float(cfg["score.normalisation.clade_floor"])
    weights = cfg.m_axis_weights
    sw = cfg["score.weight_sweep"]
    focal_tissue = cfg.tissues["focal"]

    profiles = load_gene_profiles(results_dir)
    n_members = int(profiles.shape[0])
    pairs = add_components(load_pairs(results_dir, cfg, n_members), profiles, cfg)
    kw = dict(alpha=alpha, beta=beta, clade_floor=floor, m_axis_weights=weights)
    scored = score(pairs, scope="family", **kw)
    scored = score(scored, scope="clade", **kw)
    sweep = weight_sensitivity(
        pairs, scope="family", steps=int(sw["steps"]),
        alpha_min=float(sw["alpha_min"]), alpha_max=float(sw["alpha_max"]),
        clade_floor=floor, m_axis_weights=weights,
    )
    checks = validation_checks(scored, profiles, cfg)
    config.evaluate_check(
        cfg, "pair_join_completeness",
        len(pairs) == contracts.expected_pair_count(n_members),
        {"n_pairs": int(len(pairs)), "n_members": n_members,
         "expected": contracts.expected_pair_count(n_members)},
        sink=checks,
    )

    keep = [
        "label_a", "label_b", "gene_a", "gene_b", "symbol_a", "symbol_b",
        "pair_class", "duplication_mode", "chrom_a", "chrom_b",
        # raw components
        "pid_aligned", "esm2_cosine_distance", "tm_score", "rmsd",
        "pocket_identity", "tissue_overlap", "dose_ratio",
        "cover_a_by_b", "cover_b_by_a",
        f"log2fc_mean_{focal_tissue}", f"log2fc_sd_{focal_tissue}",
        "log2fc_mean_all", "log2fc_sd_all", "spearman_profile",
        # normalised components and scores
        "M_seq_family", "M_fold_family", "M_pocket_family", "M_family", "E_family",
        "R_family", "R_a_covered_by_b_family", "R_b_covered_by_a_family",
        "M_seq_clade", "M_fold_clade", "M_pocket_clade", "M_clade", "E_clade",
        "R_clade", "R_a_covered_by_b_clade", "R_b_covered_by_a_clade",
    ]
    out = scored[[c for c in keep if c in scored.columns]].copy()
    out["alpha_M"], out["beta_E"] = alpha, beta
    out = out.sort_values("R_family", ascending=False).reset_index(drop=True)

    contracts.write_csv(out, outdir / "redundancy_scores.csv",
                        contracts.redundancy_scores_spec(cfg, n_members))
    contracts.write_csv(sweep, outdir / "weight_sensitivity.csv",
                        contracts.weight_sensitivity_spec(cfg))

    manifest = {
        "config_digest": cfg.digest(),
        "n_pairs": int(len(out)),
        "score_form": "R = M**alpha * E**beta (gated, not additive)",
        "weights": {"alpha_M": alpha, "beta_E": beta, "M_axis_weights": weights},
        "M_axes": {
            "seq": "mean of normalised pid_aligned and normalised ESM2 similarity",
            "fold": "normalised TM-score (shorter-chain normalised)",
            "pocket": "normalised heme-pocket identity over 23 residues",
            "rationale": (
                "the four raw molecular features are collinear across the family "
                "(Spearman 0.85-0.94) but come apart within the Lb clade "
                "(identity vs pocket identity rho = 0.39), so sequence is "
                "collapsed to one axis rather than counted twice"
            ),
        },
        "E_terms": {
            "tissue_overlap": "histogram intersection of normalised tissue-mean CPM profiles",
            "dose_ratio": f"smaller {focal_tissue}-mean CPM over larger",
            "substitution_note": (
                "the design asked for per-cell co-expression overlap; it needs "
                "cell-type labels the expression module does not have, and the "
                "available stand-in (Spearman over 5 libraries) is exactly 1.000 "
                "for all six Lb pairs. E is therefore tissue-resolution."
            ),
        },
        "normalisation_scopes": {
            "family": "min-max over all pairs; R_family is absolute",
            "clade": (
                f"M and E each rescaled as a whole across the {CLASS_FOCAL_FOCAL} "
                f"pairs onto [{floor}, 1]; R_clade is a relative ranking inside the "
                f"clade, where the lowest value means lowest of the six observed, "
                f"not 'not redundant'. The floor exists because a plain min-max "
                f"zero in either factor would zero the product and merge pairs "
                f"that differ; it changes no ordering."
            ),
        },
        "directional_scores": (
            "R_a_covered_by_b / R_b_covered_by_a use the partner's suppliable "
            "dose fraction in place of the symmetric dose ratio"
        ),
        "family_composition": _family_composition(out),
        "metrics_not_used": cfg.get("score.metrics_not_used") or {},
        "validation_checks": checks,
        "weight_sensitivity": {
            "alpha_grid": f"{sw['alpha_min']} to {sw['alpha_max']} in {sw['steps']} steps",
            "top_focal_pair_unique": sorted(set(sweep.top_focal_pair)),
            "top_focal_pair_changes_with_alpha": bool(sweep.top_focal_pair.nunique() > 1),
        },
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return {"scores": out, "sweep": sweep, "manifest": manifest, "scored": scored}
