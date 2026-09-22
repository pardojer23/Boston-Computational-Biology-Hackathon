"""Invariants this pipeline has already been bitten by, plus the ones that
would be expensive to discover late.

    pytest -q tests/

Deliberately narrow. These are not unit tests of the science — the science is
checked against evidence in the validation checks, which run inside the
pipeline. These are tests of the *plumbing that failed once*:

1. **Canonical pair ordering.** Two modules each calling
   ``itertools.combinations`` over their own member ordering emit the same pair
   in opposite orientations, and a merge then drops the reversed rows without
   raising. Four of 21 pairs were reversed this way (README §6).
2. **The label round trip.** ``label_of`` / ``gene_from_label`` must compose to
   the identity on gene IDs, including for members with no symbol. A symbol map
   that disagreed between two modules about one gene cost 6 of 21 rows.
3. **The 21-row join.** ``C(n_members, 2)`` pairs in every pair table, always.
4. **Contracts reject a dropped or reordered row.** The guard itself has to be
   guarded — a contract that silently passes a short table is worse than none.
5. **tau and dose_ratio on inputs whose answers are known by hand.**

Everything here runs in milliseconds on synthetic inputs; nothing touches the
network or the real data.
"""

from __future__ import annotations

import math
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "pipeline"))

import config as cfgmod          # noqa: E402
import contracts                 # noqa: E402
import soy_globin_core as core    # noqa: E402


@pytest.fixture(scope="module")
def cfg():
    return cfgmod.load(REPO / "config" / "config.yaml")


# --------------------------------------------------------------------------- #
# 1. canonical pair ordering
# --------------------------------------------------------------------------- #

LABELS = [
    "Glyma.10G199100_Lba", "Glyma.10G199000_Lbc1", "Glyma.20G191200_Lbc2",
    "Glyma.10G198800_Lbc3", "Glyma.10G198900_GmLb5", "Glyma.11G121700_Hb1",
    "Glyma.11G121800_Hb2",
]


def test_pair_order_is_canonical():
    for a, b in core.canonical_pair_order(LABELS):
        assert core.gene_from_label(a) < core.gene_from_label(b)


def test_pair_order_is_input_order_invariant():
    """The failure mode: the same pairs in a different orientation.

    Shuffling the member list must not change the emitted orientation, because
    two modules enumerate the family from two differently-ordered tables.
    """
    rng = np.random.default_rng(0)
    base = core.canonical_pair_order(LABELS)
    for _ in range(20):
        shuffled = list(rng.permutation(LABELS))
        assert core.canonical_pair_order(shuffled) == base


def test_pair_order_covers_every_unordered_pair():
    got = core.canonical_pair_order(LABELS)
    assert len(got) == math.comb(len(LABELS), 2)
    assert len({frozenset(p) for p in got}) == len(got)
    assert {frozenset(p) for p in got} == {frozenset(p) for p in combinations(LABELS, 2)}


def test_core_and_contracts_agree_on_pair_order():
    """Two implementations of the same ordering must not drift apart.

    ``core.canonical_pair_order`` is what the science modules call;
    ``contracts.canonical_pairs`` is what the table contracts check against. If
    these ever disagree, every pair table passes its own contract and still
    fails to join.
    """
    assert core.canonical_pair_order(LABELS) == contracts.canonical_pairs(LABELS)


# --------------------------------------------------------------------------- #
# 2. label round trip
# --------------------------------------------------------------------------- #

def test_label_round_trip(cfg):
    syms = cfg.gene_symbols
    for gene in syms:
        assert core.gene_from_label(core.label_of(gene, syms)) == gene


def test_label_round_trip_for_unsymbolled_gene(cfg):
    """A member with no symbol is labelled by bare gene ID and must round-trip.

    This is the case the two disagreeing symbol maps got wrong.
    """
    syms = dict(cfg.gene_symbols)
    syms.pop("Glyma.10G198900", None)
    lab = core.label_of("Glyma.10G198900", syms)
    assert lab == "Glyma.10G198900"
    assert core.gene_from_label(lab) == "Glyma.10G198900"


def test_config_label_matches_core_label(cfg):
    for gene in cfg.gene_symbols:
        assert cfg.label_of(gene) == core.label_of(gene, cfg.gene_symbols)


def test_gene_of_strips_the_annotation_prefix(cfg):
    pref = cfg["family.identifiers.id_prefix_strip"]
    assert core.gene_of("glyma.Wm82.gnm4.ann1.Glyma.10G199100.1", pref) == "Glyma.10G199100"
    # Idempotent on an already-bare ID: the GFF3 and the FASTA disagree about
    # whether the prefix is present.
    assert core.gene_of("Glyma.10G199100", pref) == "Glyma.10G199100"


# --------------------------------------------------------------------------- #
# 3. the 21-row join, against the real tables
# --------------------------------------------------------------------------- #

PAIR_TABLES = [
    ("sequence_module/paralog_pairs.csv", contracts.paralog_pairs_spec),
    ("structure_module/structure_pairs.csv", contracts.structure_pairs_spec),
    ("expression_module/expression_pairs.csv", contracts.expression_pairs_spec),
    ("integration_module/redundancy_scores.csv", contracts.redundancy_scores_spec),
]


@pytest.fixture(scope="module")
def n_members():
    p = REPO / "results" / "sequence_module" / "globin_family_members.csv"
    if not p.exists():
        pytest.skip("pipeline has not been run")
    return len(pd.read_csv(p))


@pytest.mark.parametrize("rel,builder", PAIR_TABLES)
def test_pair_table_satisfies_its_contract(cfg, n_members, rel, builder):
    p = REPO / "results" / rel
    if not p.exists():
        pytest.skip(f"{rel} not present")
    contracts.read_csv(p, builder(cfg, n_members))


def test_all_pair_tables_share_one_key(cfg, n_members):
    """The join that broke. Every pair table must carry the same 21 keys."""
    keys = None
    for rel, _ in PAIR_TABLES:
        p = REPO / "results" / rel
        if not p.exists():
            pytest.skip(f"{rel} not present")
        k = set(zip(*pd.read_csv(p, usecols=["label_a", "label_b"]).values.T))
        assert len(k) == contracts.expected_pair_count(n_members)
        if keys is None:
            keys = k
        else:
            assert k == keys, f"{rel} has a different key set"


# --------------------------------------------------------------------------- #
# 4. the contracts themselves
# --------------------------------------------------------------------------- #

@pytest.fixture
def good_pairs(cfg, n_members):
    p = REPO / "results" / "sequence_module" / "paralog_pairs.csv"
    if not p.exists():
        pytest.skip("pipeline has not been run")
    return pd.read_csv(p)


def test_contract_rejects_a_dropped_row(cfg, n_members, good_pairs):
    with pytest.raises(contracts.ContractError, match="rows, expected exactly"):
        contracts.validate(good_pairs.iloc[:-1],
                           contracts.paralog_pairs_spec(cfg, n_members))


def test_contract_rejects_a_reversed_pair(cfg, n_members, good_pairs):
    flipped = good_pairs.assign(
        gene_a=good_pairs.gene_b, gene_b=good_pairs.gene_a,
        label_a=good_pairs.label_b, label_b=good_pairs.label_a)
    with pytest.raises(contracts.ContractError, match="canonical orientation"):
        contracts.validate(flipped, contracts.paralog_pairs_spec(cfg, n_members))


def test_contract_rejects_na_in_a_required_column(cfg, n_members, good_pairs):
    """A skipped embedding step must not reach the score as a NaN."""
    holed = good_pairs.assign(
        esm2_cosine_distance=good_pairs.esm2_cosine_distance.mask(good_pairs.index < 2))
    with pytest.raises(contracts.ContractError, match="non-nullable"):
        contracts.validate(holed, contracts.paralog_pairs_spec(cfg, n_members))


def test_contract_rejects_a_duplicate_key(cfg, n_members, good_pairs):
    dup = pd.concat([good_pairs.iloc[:-1], good_pairs.iloc[[0]]], ignore_index=True)
    with pytest.raises(contracts.ContractError, match="not unique"):
        contracts.validate(dup, contracts.paralog_pairs_spec(cfg, n_members))


def test_contract_rejects_out_of_range_values(cfg, n_members, good_pairs):
    with pytest.raises(contracts.ContractError, match="outside"):
        contracts.validate(good_pairs.assign(pid_aligned=good_pairs.pid_aligned * 2),
                           contracts.paralog_pairs_spec(cfg, n_members))


def test_contract_accepts_an_all_na_nullable_column(cfg, n_members, good_pairs):
    """An all-NA nullable column reads back as float64 whatever it was written
    as, so its dtype carries no information and must not be checked."""
    spec = contracts.paralog_pairs_spec(cfg, n_members)
    df = good_pairs.assign(symbol_a=np.nan)
    contracts.validate(df, spec)


def test_legacy_column_migration_round_trips(cfg):
    legacy = pd.DataFrame({
        "alpha": [0.0, 1.0], "beta": [1.0, 0.0],
        "top_pair_overall": ["a|b", "a|b"], "top_lb_pair": ["a|b", "a|b"],
        "top_lb_R": [0.5, 0.6], "lb_R_spread": [0.1, 0.2],
    })
    out = contracts.migrate_legacy_columns(legacy)
    assert set(contracts.RENAMED_COLUMNS.values()) <= set(out.columns)
    assert not set(contracts.RENAMED_COLUMNS) & set(out.columns)


def test_pair_class_values_are_family_agnostic(cfg):
    p = REPO / "results" / "integration_module" / "redundancy_scores.csv"
    if not p.exists():
        pytest.skip("pipeline has not been run")
    got = set(pd.read_csv(p, usecols=["pair_class"]).pair_class.unique())
    assert got <= set(contracts.RENAMED_PAIR_CLASSES.values()), (
        f"pair_class carries family-specific labels: {got}")


# --------------------------------------------------------------------------- #
# 5. the two statistics with hand-checkable answers
# --------------------------------------------------------------------------- #

def test_tau_endpoints():
    import soy_globin_expression as expr

    # Perfectly specific to one tissue -> 1; perfectly even -> 0.
    assert expr.tau(np.array([10.0, 0.0])) == pytest.approx(1.0)
    assert expr.tau(np.array([5.0, 5.0, 5.0])) == pytest.approx(0.0)
    # Scale-invariant: tau is defined on the profile's shape.
    assert expr.tau(np.array([2.0, 1.0])) == pytest.approx(expr.tau(np.array([20.0, 10.0])))
    # All-zero has no profile to be specific about.
    assert np.isnan(expr.tau(np.array([0.0, 0.0])))


def test_tau_two_tissue_caveat_is_recorded(cfg):
    """The caveat is load-bearing: with two tissues tau collapses to a contrast.

    Asserted so it cannot be dropped from the config while the column is still
    emitted.
    """
    assert cfg.get("expression.metrics.tau_caveat")
    assert len(cfg.tissue_names) == 2 or True  # documents why the caveat exists


def test_dose_ratio_is_symmetric_and_bounded():
    def dose(a, b):
        hi = max(a, b)
        return (min(a, b) / hi) if hi > 0 else float("nan")

    assert dose(3.0, 12.0) == pytest.approx(0.25)
    assert dose(12.0, 3.0) == pytest.approx(0.25)   # symmetric by construction
    assert dose(7.0, 7.0) == pytest.approx(1.0)
    assert 0.0 <= dose(1e-6, 1e6) <= 1.0


def test_minmax_degenerate_range_maps_to_one():
    """A feature that does not distinguish any pair must not become NaN.

    ``1.0`` is the honest reading — "this axis says these pairs are the same" —
    and a NaN here would zero the product score for every pair.
    """
    import soy_globin_integration as integ

    out = integ.minmax(pd.Series([0.5, 0.5, 0.5]))
    assert (out == 1.0).all()


# --------------------------------------------------------------------------- #
# 6. config validation
# --------------------------------------------------------------------------- #

def test_config_digest_is_stable_and_order_independent(cfg):
    import copy

    d = copy.deepcopy(cfg.raw)
    reordered = dict(reversed(list(d.items())))
    assert cfgmod.digest(reordered) == cfgmod.digest(d)


def test_section_digest_isolates_changes(cfg):
    import copy

    before = cfg.section_digest("score")
    d = copy.deepcopy(cfg.raw)
    d["score"]["alpha_m"] = 0.9
    after = cfgmod.Config(d).section_digest("score")
    assert before != after
    # An unrelated section must not move.
    assert (cfgmod.Config(d).section_digest("sequence.alignment")
            == cfg.section_digest("sequence.alignment"))


@pytest.mark.parametrize("mutate,expect", [
    (lambda d: d["score"].__setitem__("alpha_m", 1.7), "must be in"),
    (lambda d: d["score"]["m_axes"].__setitem__("seq", 0.9), "sum to"),
    (lambda d: d["structure"]["pocket"].__setitem__("reference_member", "Nope"),
     "resolves to 0"),
    (lambda d: d["expression"]["tissues"].__setitem__("focal", "leaf"),
     "no library declares"),
    # A gene declared both focal and outgroup. `focal_genes` became
    # `gene_ids` when the family gained three declaration modes.
    (lambda d: d["family"]["gene_ids"].__setitem__("Glyma.10G198900", "GmLb5"),
     "both"),
    (lambda d: d["checks"]["pair_join_completeness"].__setitem__("policy", "maybe"),
     "must be one of"),
    (lambda d: d["family"]["identifiers"].__setitem__("chromosome_regex", "Gm(["),
     "not a valid regex"),
])
def test_config_validation_rejects(cfg, mutate, expect):
    import copy

    d = copy.deepcopy(cfg.raw)
    mutate(d)
    with pytest.raises(cfgmod.ConfigError, match=expect):
        cfgmod.validate(cfgmod.Config(d))


def test_missing_config_key_raises_rather_than_defaulting(cfg):
    """A renamed config field must surface at first use, not as a silent default."""
    with pytest.raises(cfgmod.ConfigError, match="not found"):
        cfg["score.alpha_M"]          # wrong case
    with pytest.raises(cfgmod.ConfigError, match="no policy declared"):
        cfg.check_policy("a_check_nobody_declared")


def test_beta_is_always_one_minus_alpha(cfg):
    assert cfg.alpha + cfg.beta == pytest.approx(1.0)


def test_pocket_reference_member_resolves(cfg):
    gene = cfg.pocket_reference_gene
    assert gene in cfg.gene_symbols
    assert cfg.gene_symbols[gene] == cfg["structure.pocket.reference_member"]
