"""Declared schemas for every table that crosses a module boundary.

This generalises the one guard the repo already had — ``load_pairs()`` raising
rather than returning a short table — into a contract applied on *write* as
well as on read, at every boundary rather than one.

The defect this exists to prevent is on record: two modules disagreed about
whether ``Glyma.10G198900`` carries a symbol, and each enumerated pairs with
its own ``itertools.combinations``, so a naive merge returned 11 of 21 rows and
raised nothing (README §6). A short table that looks well-formed is the worst
failure mode available to this pipeline, because every number downstream of it
is arithmetically valid and scientifically meaningless.

Three rules follow from that:

* **Validate on write.** An emitter that produces a malformed table fails at
  the emitter, not four steps later at the consumer.
* **Row count is part of the schema.** ``C(n_members, 2)`` pairs, exactly.
* **Known-degenerate columns are declared, not documented.** ``spearman_profile``
  is exactly 1.000 for all six focal pairs and carries no information
  (README §3.3). It is marked ``degenerate`` here so a consumer that reads it
  has to acknowledge that, instead of a reader having to find caveat §9.10.

Specs are *built from the config* rather than hardcoded, because several column
names are tissue-dependent (``cpm_mean_nodule``, ``log2fc_mean_nodule``) and
several are library-dependent (``cpm_GSM7065810``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd


class ContractError(ValueError):
    """A table does not satisfy its declared schema."""


# --------------------------------------------------------------------------- #
# Spec types
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ColumnSpec:
    """One column's contract.

    ``nullable`` means NA is *expected* — used for the NA-with-a-reason
    convention (``detection_rate``, ``coexpression_overlap``), which is a
    deliberate signal and not missing data.

    ``degenerate`` marks a column that is emitted but known to carry no
    information at the current resolution; ``note`` says why.
    """

    name: str
    kind: str = "any"          # float | int | str | bool | any
    nullable: bool = False
    bounds: tuple[float, float] | None = None
    degenerate: bool = False
    note: str | None = None


@dataclass(frozen=True)
class TableSpec:
    name: str
    columns: tuple[ColumnSpec, ...]
    key: tuple[str, ...] = ()
    #: Exact row count. ``None`` means unconstrained.
    n_rows: int | None = None
    #: When set, require ``df[gene_a] < df[gene_b]`` lexicographically on every
    #: row — the canonical pair orientation that makes the four pair tables
    #: merge without reordering.
    pair_orientation: tuple[str, str] | None = None
    #: Require the key to be unique.
    unique_key: bool = True

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.columns]

    def column(self, name: str) -> ColumnSpec:
        for c in self.columns:
            if c.name == name:
                return c
        raise KeyError(f"{self.name} has no column {name!r}")

    @property
    def degenerate_columns(self) -> dict[str, str | None]:
        return {c.name: c.note for c in self.columns if c.degenerate}


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #

_KIND_CHECK = {
    "float": lambda s: pd.api.types.is_numeric_dtype(s),
    "int": lambda s: pd.api.types.is_integer_dtype(s) or (
        pd.api.types.is_numeric_dtype(s)
        and bool(np.all(np.isclose(s.dropna() % 1, 0)))
    ),
    "str": lambda s: pd.api.types.is_object_dtype(s) or pd.api.types.is_string_dtype(s),
    "bool": lambda s: pd.api.types.is_bool_dtype(s) or set(s.dropna().unique()) <= {
        True, False, "True", "False", 0, 1},
    "any": lambda s: True,
}


def validate(df: pd.DataFrame, spec: TableSpec, *, allow_extra: bool = True) -> pd.DataFrame:
    """Raise ``ContractError`` unless ``df`` satisfies ``spec``. Returns ``df``.

    Extra columns are allowed by default: modules legitimately add diagnostic
    columns, and forbidding them would make every schema change a two-file
    edit. Missing, mistyped, out-of-range, non-unique and mis-oriented are all
    errors.
    """
    problems: list[str] = []

    missing = [c for c in spec.names if c not in df.columns]
    if missing:
        problems.append(f"missing columns: {missing}")
    if not allow_extra:
        extra = [c for c in df.columns if c not in spec.names]
        if extra:
            problems.append(f"unexpected columns: {extra}")

    if spec.n_rows is not None and len(df) != spec.n_rows:
        problems.append(f"has {len(df)} rows, expected exactly {spec.n_rows}")

    for col in spec.columns:
        if col.name not in df.columns:
            continue
        s = df[col.name]
        if not col.nullable:
            n_na = int(s.isna().sum())
            if n_na:
                idx = s.index[s.isna()][:3].tolist()
                problems.append(
                    f"{col.name}: {n_na} NA value(s) in a non-nullable column "
                    f"(first rows {idx})"
                )
        checker = _KIND_CHECK.get(col.kind, _KIND_CHECK["any"])
        # An all-NA nullable column comes back from read_csv as float64 whatever
        # it was written as, so its dtype carries no information — checking it
        # would fail legitimately empty columns such as the symbol column of a
        # marker table, or detection_rate at pseudobulk resolution.
        all_na = bool(s.isna().all())
        if not (col.nullable and all_na) and not checker(s):
            problems.append(f"{col.name}: dtype {s.dtype} is not compatible with {col.kind!r}")
        if col.bounds is not None and pd.api.types.is_numeric_dtype(s):
            lo, hi = col.bounds
            bad = s.dropna()
            bad = bad[(bad < lo) | (bad > hi)]
            if len(bad):
                problems.append(
                    f"{col.name}: {len(bad)} value(s) outside [{lo}, {hi}], "
                    f"e.g. {bad.iloc[0]!r} at row {bad.index[0]}"
                )

    if spec.key and all(k in df.columns for k in spec.key):
        if spec.unique_key:
            dup = df.duplicated(subset=list(spec.key), keep=False)
            if dup.any():
                ex = df.loc[dup, list(spec.key)].head(2).to_dict("records")
                problems.append(f"key {spec.key} is not unique ({int(dup.sum())} rows, e.g. {ex})")

    if spec.pair_orientation:
        a, b = spec.pair_orientation
        if a in df.columns and b in df.columns:
            bad = [f"{x}|{y}" for x, y in zip(df[a], df[b]) if not str(x) < str(y)]
            if bad:
                problems.append(
                    f"{len(bad)} row(s) not in canonical orientation {a} < {b} "
                    f"(first {bad[0]}); every emitter must iterate the shared "
                    f"pair-ordering function"
                )

    if problems:
        raise ContractError(
            f"table {spec.name!r} violates its contract:\n  - "
            + "\n  - ".join(problems)
        )
    return df


def write_csv(df: pd.DataFrame, path: str | Path, spec: TableSpec, **kw) -> Path:
    """Validate then write. An emitter cannot produce a malformed table."""
    validate(df, spec)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=kw.pop("index", False), **kw)
    return path


def read_csv(path: str | Path, spec: TableSpec, **kw) -> pd.DataFrame:
    """Read then validate, so a hand-edited or stale file is caught on load."""
    df = pd.read_csv(path, **kw)
    try:
        return validate(df, spec)
    except ContractError as e:
        raise ContractError(f"{path}: {e}") from None


def expected_pair_count(n_members: int) -> int:
    return math.comb(n_members, 2)


def canonical_pairs(labels: Iterable[str]) -> list[tuple[str, str]]:
    """The one pair enumeration: sorted by bare gene ID, ``gene_a < gene_b``.

    Kept here as well as in ``soy_globin_core`` so contracts can be checked
    without importing the science module; ``core.canonical_pair_order``
    delegates to this.
    """
    by_gene = sorted(labels, key=lambda l: l.split("_")[0])
    return [(a, b) if a.split("_")[0] < b.split("_")[0] else (b, a)
            for a, b in combinations(by_gene, 2)]


# --------------------------------------------------------------------------- #
# Shared column groups
# --------------------------------------------------------------------------- #

PAIR_KEY = (
    ColumnSpec("label_a", "str"),
    ColumnSpec("label_b", "str"),
    ColumnSpec("gene_a", "str"),
    ColumnSpec("gene_b", "str"),
)

GENE_KEY = (
    ColumnSpec("gene_id", "str"),
    ColumnSpec("label", "str"),
)


def _pair_spec(name: str, extra: Sequence[ColumnSpec], n_pairs: int | None) -> TableSpec:
    return TableSpec(
        name=name,
        columns=PAIR_KEY + tuple(extra),
        key=("label_a", "label_b"),
        n_rows=n_pairs,
        pair_orientation=("gene_a", "gene_b"),
    )


# --------------------------------------------------------------------------- #
# Spec builders. `n_members` is passed rather than inferred so the row-count
# contract is asserted against the family the config declares, not against
# whatever the table happens to contain.
# --------------------------------------------------------------------------- #


def members_spec(cfg, n_members: int | None = None) -> TableSpec:
    return TableSpec(
        name="globin_family_members",
        columns=GENE_KEY + (
            ColumnSpec("protein_id", "str"),
            ColumnSpec("symbol", "str", nullable=True,
                       note="null for a member with no symbol in the family config"),
            ColumnSpec("is_focal", "bool"),
            ColumnSpec("prot_len", "int", bounds=(1, 1e5)),
            ColumnSpec("pfam_hit", "bool"),
            ColumnSpec("evalue_full", "float", nullable=True),
            ColumnSpec("score_full", "float", nullable=True),
            ColumnSpec("dom_score", "float", nullable=True),
            ColumnSpec("hmm_coverage", "float", nullable=True, bounds=(0.0, 1.0)),
            ColumnSpec("ali_from", "float", nullable=True),
            ColumnSpec("ali_to", "float", nullable=True),
        ),
        key=("gene_id",),
        n_rows=n_members,
    )


def gene_context_spec(cfg, n_members: int | None = None) -> TableSpec:
    base = members_spec(cfg, n_members)
    return TableSpec(
        name="gene_context",
        columns=base.columns + (
            ColumnSpec("seqid", "str"),
            ColumnSpec("chrom", "str", nullable=True),
            ColumnSpec("is_scaffold", "bool"),
            ColumnSpec("start", "int"),
            ColumnSpec("end", "int"),
            ColumnSpec("strand", "str"),
            ColumnSpec("rank_on_seqid", "int"),
        ),
        key=("gene_id",),
        n_rows=n_members,
    )


def uniprot_accessions_spec(cfg, n_members: int | None = None) -> TableSpec:
    """Persisted gene -> accession mapping (audit 2.1).

    New in the refactor: previously this was a live UniProt query on every run
    with no cached output, making the least stable input in the pipeline also
    the only unrecorded one.
    """
    return TableSpec(
        name="uniprot_accessions",
        columns=GENE_KEY + (
            ColumnSpec("uniprot", "str"),
            ColumnSpec("entry_name", "str"),
            ColumnSpec("uniprot_protein_name", "str"),
            ColumnSpec("uniprot_length", "int", bounds=(1, 1e5)),
            ColumnSpec("uniprot_reviewed", "bool"),
            ColumnSpec("n_uniprot_hits", "int", bounds=(1, 1e6)),
            ColumnSpec("uniprot_release", "str",
                       note="release the mapping was resolved against; a change here "
                            "can change the accession and therefore the structure"),
            ColumnSpec("resolution", "str",
                       note="'live' or 'pinned' (structure.uniprot.pinned)"),
        ),
        key=("gene_id",),
        n_rows=n_members,
    )


def paralog_pairs_spec(cfg, n_members: int | None = None) -> TableSpec:
    n = expected_pair_count(n_members) if n_members else None
    return _pair_spec("paralog_pairs", (
        ColumnSpec("symbol_a", "str", nullable=True),
        ColumnSpec("symbol_b", "str", nullable=True),
        ColumnSpec("chrom_a", "str", nullable=True),
        ColumnSpec("chrom_b", "str", nullable=True),
        ColumnSpec("same_seqid", "bool"),
        ColumnSpec("n_intervening_genes", "float", nullable=True),
        ColumnSpec("intergenic_bp", "float", nullable=True),
        ColumnSpec("duplication_mode", "str",
                   note="adjacency, not synteny; an annotation column and "
                        "deliberately not a score term (README §4.5)"),
        ColumnSpec("both_focal", "bool",
                   note="the single source of focal membership for downstream "
                        "pair_class; symbols do not distinguish focal members"),
        ColumnSpec("pid_aligned", "float", bounds=(0.0, 100.0)),
        ColumnSpec("pid_shorter", "float", bounds=(0.0, 100.0)),
        ColumnSpec("n_aligned_cols", "int", bounds=(1, 1e5)),
        # Required, not optional: if the embedding step is skipped this column
        # is absent or NA and the failure otherwise surfaces four steps later
        # as a NaN M_seq (audit 1.6).
        ColumnSpec("esm2_cosine_distance", "float", bounds=(0.0, 2.0)),
    ), n)


def structure_pairs_spec(cfg, n_members: int | None = None) -> TableSpec:
    n = expected_pair_count(n_members) if n_members else None
    return _pair_spec("structure_pairs", (
        ColumnSpec("tm_score", "float", bounds=(0.0, 1.0),
                   note="normalised by the shorter chain"),
        ColumnSpec("tm_norm_a", "float", bounds=(0.0, 1.0)),
        ColumnSpec("tm_norm_b", "float", bounds=(0.0, 1.0)),
        ColumnSpec("rmsd", "float", bounds=(0.0, 1e3)),
        ColumnSpec("pocket_identity", "float", bounds=(0.0, 100.0)),
        ColumnSpec("n_pocket_cols", "int", bounds=(1, 1e4)),
        ColumnSpec("n_res_a", "int", bounds=(1, 1e5)),
        ColumnSpec("n_res_b", "int", bounds=(1, 1e5)),
    ), n)


def expression_pairs_spec(cfg, n_members: int | None = None) -> TableSpec:
    n = expected_pair_count(n_members) if n_members else None
    focal = cfg.tissues["focal"]
    degen = cfg.get("expression.metrics.degenerate_columns") or {}
    notcomp = cfg.get("expression.metrics.not_computable") or {}
    return _pair_spec("expression_pairs", (
        ColumnSpec("spearman_profile", "float", bounds=(-1.0, 1.0),
                   degenerate="spearman_profile" in degen,
                   note=degen.get("spearman_profile")),
        ColumnSpec("n_profile_points", "int", bounds=(2, 1e4)),
        ColumnSpec("log2fc_mean_all", "float"),
        ColumnSpec("log2fc_sd_all", "float", nullable=True),
        ColumnSpec("n_libraries_all", "int", bounds=(1, 1e4)),
        ColumnSpec(f"log2fc_mean_{focal}", "float"),
        ColumnSpec(f"log2fc_sd_{focal}", "float", nullable=True),
        ColumnSpec(f"n_libraries_{focal}", "int", bounds=(1, 1e4)),
        ColumnSpec("log2fc_pseudocount_cpm", "float", bounds=(0.0, 1e6)),
        ColumnSpec(f"cpm_mean_{focal}_a", "float", bounds=(0.0, 1e7)),
        ColumnSpec(f"cpm_mean_{focal}_b", "float", bounds=(0.0, 1e7)),
        ColumnSpec("dose_ratio", "float", bounds=(0.0, 1.0)),
        ColumnSpec("cover_a_by_b", "float", bounds=(0.0, 1.0)),
        ColumnSpec("cover_b_by_a", "float", bounds=(0.0, 1.0)),
        # NA by design at pseudobulk resolution, with the reason in its own
        # column. Nullable, never silently replaced by a look-alike statistic.
        ColumnSpec("coexpression_overlap", "float", nullable=True, bounds=(0.0, 1.0),
                   note=notcomp.get("coexpression_overlap")),
        ColumnSpec("coexpression_overlap_note", "str", nullable=True),
    ), n)


def gene_profiles_spec(cfg, n_rows: int | None = None, *, marker: bool = False) -> TableSpec:
    notcomp = cfg.get("expression.metrics.not_computable") or {}
    cols: list[ColumnSpec] = [
        *GENE_KEY,
        ColumnSpec("symbol", "str", nullable=True),
        ColumnSpec("in_matrix", "bool",
                   note="False means the gene is absent from the series feature list"),
    ]
    for gsm in cfg.library_ids:
        cols.append(ColumnSpec(f"cpm_{gsm}", "float", nullable=True, bounds=(0.0, 1e7)))
    for tissue in cfg.tissue_names:
        cols.append(ColumnSpec(cfg.tissue_mean_column(tissue), "float",
                               nullable=True, bounds=(0.0, 1e7)))
    cols += [
        ColumnSpec("tau_over_tissue_means", "float", nullable=True, bounds=(0.0, 1.0),
                   note=cfg.get("expression.metrics.tau_caveat")),
        ColumnSpec("n_tissues_for_tau", "int", nullable=True),
        ColumnSpec("detection_rate", "float", nullable=True, bounds=(0.0, 1.0),
                   note=notcomp.get("detection_rate")),
        ColumnSpec("detection_rate_note", "str", nullable=True),
    ]
    if marker:
        cols.append(ColumnSpec("marker_role", "str"))
    return TableSpec(
        name="marker_profiles" if marker else "gene_pseudobulk_profiles",
        columns=tuple(cols),
        key=("gene_id",),
        n_rows=n_rows,
    )


def redundancy_scores_spec(cfg, n_members: int | None = None) -> TableSpec:
    n = expected_pair_count(n_members) if n_members else None
    focal = cfg.tissues["focal"]
    scopes = cfg["score.normalisation.scopes_emitted"]
    cols: list[ColumnSpec] = [
        ColumnSpec("symbol_a", "str", nullable=True),
        ColumnSpec("symbol_b", "str", nullable=True),
        ColumnSpec("pair_class", "str",
                   note="focal-focal | focal-other | other-other, derived from both_focal"),
        ColumnSpec("duplication_mode", "str"),
        ColumnSpec("chrom_a", "str", nullable=True),
        ColumnSpec("chrom_b", "str", nullable=True),
        ColumnSpec("pid_aligned", "float", bounds=(0.0, 100.0)),
        ColumnSpec("esm2_cosine_distance", "float", bounds=(0.0, 2.0)),
        ColumnSpec("tm_score", "float", bounds=(0.0, 1.0)),
        ColumnSpec("rmsd", "float"),
        ColumnSpec("pocket_identity", "float", bounds=(0.0, 100.0)),
        ColumnSpec("tissue_overlap", "float", bounds=(0.0, 1.0)),
        ColumnSpec("dose_ratio", "float", bounds=(0.0, 1.0)),
        ColumnSpec("cover_a_by_b", "float", bounds=(0.0, 1.0)),
        ColumnSpec("cover_b_by_a", "float", bounds=(0.0, 1.0)),
        ColumnSpec(f"log2fc_mean_{focal}", "float"),
        ColumnSpec(f"log2fc_sd_{focal}", "float", nullable=True),
        ColumnSpec("log2fc_mean_all", "float"),
        ColumnSpec("log2fc_sd_all", "float", nullable=True),
        ColumnSpec("spearman_profile", "float", bounds=(-1.0, 1.0), degenerate=True,
                   note="carried for inspection; not a score term"),
    ]
    for scope in scopes:
        for base in ("M_seq", "M_fold", "M_pocket", "M", "E", "R",
                     "R_a_covered_by_b", "R_b_covered_by_a"):
            cols.append(ColumnSpec(f"{base}_{scope}", "float", bounds=(0.0, 1.0)))
    cols += [
        ColumnSpec("alpha_M", "float", bounds=(0.0, 1.0)),
        ColumnSpec("beta_E", "float", bounds=(0.0, 1.0)),
    ]
    return _pair_spec("redundancy_scores", cols, n)


#: Columns renamed by the config-driven refactor, old -> new. The old names
#: hardcoded the family into the schema ("lb", "Lb-Lb"), which is audit finding
#: 6 showing up in a column header rather than in code. Values are unchanged;
#: only the names are. Anything reading a pre-refactor CSV — including the
#: tables committed at a2dfe7e and the numbers quoted in README §4 — needs this
#: mapping, so it is code rather than a changelog entry.
RENAMED_COLUMNS = {
    "top_lb_pair": "top_focal_pair",
    "top_lb_R": "top_focal_R",
    "lb_R_spread": "focal_R_spread",
}

#: Likewise for the values of ``pair_class``.
RENAMED_PAIR_CLASSES = {
    "Lb-Lb": "focal-focal",
    "Lb-other": "focal-other",
    "other-other": "other-other",
}


def migrate_legacy_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename pre-refactor columns and pair_class values in place-ish.

    Lets a legacy table be read against a current spec, which is what the
    end-to-end verification diff needs in order to compare values rather than
    headers.
    """
    out = df.rename(columns={k: v for k, v in RENAMED_COLUMNS.items() if k in df.columns})
    if "pair_class" in out.columns:
        out["pair_class"] = out["pair_class"].replace(RENAMED_PAIR_CLASSES)
    return out


def weight_sensitivity_spec(cfg) -> TableSpec:
    """Contract for the alpha sweep.

    Note the renamed columns: ``top_lb_pair`` -> ``top_focal_pair`` and so on,
    per ``RENAMED_COLUMNS``. A table written before the refactor fails this
    spec until passed through ``migrate_legacy_columns``, which is intended —
    a silent accept would leave the family hardcoded in the schema.
    """
    return TableSpec(
        name="weight_sensitivity",
        columns=(
            ColumnSpec("alpha", "float", bounds=(0.0, 1.0)),
            ColumnSpec("beta", "float", bounds=(0.0, 1.0)),
            ColumnSpec("top_pair_overall", "str"),
            ColumnSpec("top_focal_pair", "str"),
            ColumnSpec("top_focal_R", "float", bounds=(0.0, 1.0)),
            ColumnSpec("focal_R_spread", "float", bounds=(0.0, 1.0)),
        ),
        key=("alpha",),
        n_rows=int(cfg["score.weight_sweep.steps"]),
    )


#: Registry so a CLI shim can look a spec up by output basename.
SPEC_BUILDERS = {
    "globin_family_members": members_spec,
    "gene_context": gene_context_spec,
    "uniprot_accessions": uniprot_accessions_spec,
    "paralog_pairs": paralog_pairs_spec,
    "structure_pairs": structure_pairs_spec,
    "expression_pairs": expression_pairs_spec,
    "gene_pseudobulk_profiles": gene_profiles_spec,
    "redundancy_scores": redundancy_scores_spec,
    "weight_sensitivity": lambda cfg, n=None: weight_sensitivity_spec(cfg),
}
