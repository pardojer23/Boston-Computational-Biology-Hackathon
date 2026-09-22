"""The family-source resolver and the empty-outgroup path.

    pytest -q tests/test_family_source.py

The family can be declared three ways and all three reduce to one seed table.
Two things in that change are easy to get wrong and are what these tests are
really for:

1. **`focal` is optional for `pfam` and `orthodb`.** Omitting it means the
   whole discovered set is the family, with no outgroup. Then
   ``focal_outgroup_separation`` has an empty comparison, and reporting that as
   ``passed`` would turn a ``policy: fail`` gate into a silent green light —
   the run would go green on a check that never ran.
2. **Naming a gene must not make it a member.** Symbols are labels; only
   ``must_include`` rows force a gene into the family. Otherwise two runs with
   the same Pfam accession could disagree about membership depending on which
   symbols happened to be declared.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "pipeline"))

import config as cfgmod      # noqa: E402
import contracts             # noqa: E402
import family as fammod      # noqa: E402


@pytest.fixture(scope="module")
def base():
    return cfgmod.load(REPO / "config" / "config.yaml")


def _pfam_cfg(base, focal=None, symbols=True):
    d = copy.deepcopy(base.raw)
    d["family"].pop("gene_ids", None)
    d["family"].pop("focal_genes", None)
    d["family"]["source"] = "pfam"
    d["family"]["pfam"] = {"accession": "PF00042"}
    if focal:
        d["family"]["pfam"]["focal"] = list(focal)
    if not symbols:
        d["family"].pop("outgroup_symbols", None)
        d["family"].pop("symbols", None)
        d["structure"]["pocket"]["reference_member"] = "Glyma.10G199100"
    return cfgmod.validate(cfgmod.Config(d))


@pytest.fixture(scope="module")
def members():
    p = REPO / "results" / "sequence_module" / "globin_family_members.csv"
    if not p.exists():
        pytest.skip("pipeline has not been run")
    return pd.read_csv(p)


# --------------------------------------------------------------------------- #
# the seed is the same shape whatever declared it
# --------------------------------------------------------------------------- #

def test_gene_ids_seed_satisfies_the_contract(base):
    seed = fammod.resolve(base)
    contracts.validate(seed, contracts.family_seed_spec(base))


def test_pfam_seed_satisfies_the_contract(base):
    cfg = _pfam_cfg(base)
    contracts.validate(fammod.resolve(cfg), contracts.family_seed_spec(cfg))


def test_pfam_seed_forces_no_genes_in(base):
    """Declaring a symbol must not add a gene to the family.

    Otherwise the same accession could yield different membership depending on
    which labels happened to be configured.
    """
    seed = fammod.resolve(_pfam_cfg(base))
    assert int(seed.must_include.sum()) == 0
    assert len(seed) > 0, "symbols should still be carried, just not forced in"


def test_gene_ids_seed_forces_the_declared_genes_in(base):
    seed = fammod.resolve(base)
    forced = set(seed.loc[seed.must_include, "gene_id"])
    assert forced == set(base.declared_focal)


def test_declared_focal_genes_are_forced_in_under_pfam(base):
    cfg = _pfam_cfg(base, focal=["Glyma.10G199100", "Glyma.10G199000"])
    seed = fammod.resolve(cfg)
    assert set(seed.loc[seed.must_include, "gene_id"]) == {
        "Glyma.10G199100", "Glyma.10G199000"}


# --------------------------------------------------------------------------- #
# focal resolution — the optional-focal path
# --------------------------------------------------------------------------- #

def test_no_declared_focal_makes_every_member_focal(base, members):
    seed = fammod.resolve(_pfam_cfg(base))
    resolved = fammod.resolve_focal(seed, members)
    assert resolved.all(), "with no focal subset every discovered member is focal"
    assert int((~resolved).sum()) == 0, "and there is no outgroup"


def test_declared_focal_is_respected(base, members):
    cfg = _pfam_cfg(base, focal=["Glyma.10G199100", "Glyma.10G199000"])
    resolved = fammod.resolve_focal(fammod.resolve(cfg), members)
    assert int(resolved.sum()) == 2
    assert int((~resolved).sum()) == len(members) - 2


def test_gene_ids_mode_reproduces_the_committed_focal_set(base, members):
    """The default path must still give the committed 4-focal / 3-outgroup split."""
    resolved = fammod.resolve_focal(fammod.resolve(base), members)
    assert set(members.loc[resolved, "gene_id"]) == set(base.declared_focal)
    assert int(resolved.sum()) == 4 and int((~resolved).sum()) == 3


# --------------------------------------------------------------------------- #
# an empty comparison is not a pass
# --------------------------------------------------------------------------- #

def test_inapplicable_check_does_not_pass_and_does_not_raise(base):
    """A policy=fail check with no evidence must neither fail nor pass."""
    sink = {}
    row = cfgmod.evaluate_check(
        base, "focal_outgroup_separation", False, {"n_other_pairs": 0},
        sink=sink, applicable=False, reason="no outgroup pairs")
    assert base.check_policy("focal_outgroup_separation") == "fail"
    assert row["applicable"] is False
    assert "passed" not in row, "an inapplicable check must not report a verdict"
    assert sink["focal_outgroup_separation"]["reason"]


def test_inapplicable_requires_a_reason(base):
    with pytest.raises(ValueError, match="requires a reason"):
        cfgmod.evaluate_check(base, "focal_outgroup_separation", False, {},
                              applicable=False)


def test_applicable_failing_check_still_raises(base):
    """The gate must stay armed for families that do have an outgroup."""
    with pytest.raises(cfgmod.CheckFailure):
        cfgmod.evaluate_check(base, "focal_outgroup_separation", False,
                              {"n_other_pairs": 15})


def test_applicable_flag_is_recorded_on_normal_checks(base):
    row = cfgmod.evaluate_check(base, "focal_outgroup_separation", True, {})
    assert row["applicable"] is True and row["passed"] is True


# --------------------------------------------------------------------------- #
# config validation of the union
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("mutate,expect", [
    (lambda d: d["family"].__setitem__("source", "ensembl"), "must be one of"),
    (lambda d: (d["family"].__setitem__("source", "pfam"),
                d["family"].__setitem__("pfam", {"accession": "PF00042"})),
     "also populated"),
    (lambda d: d["family"].__setitem__("source", "orthodb"), "not populated"),
])
def test_source_union_is_enforced(base, mutate, expect):
    d = copy.deepcopy(base.raw)
    mutate(d)
    with pytest.raises(cfgmod.ConfigError, match=expect):
        cfgmod.validate(cfgmod.Config(d))


def test_singleton_focal_is_rejected(base):
    d = copy.deepcopy(base.raw)
    d["family"].pop("gene_ids"); d["family"].pop("focal_genes", None)
    d["family"]["source"] = "pfam"
    d["family"]["pfam"] = {"accession": "PF00042", "focal": ["Glyma.10G199100"]}
    with pytest.raises(cfgmod.ConfigError, match="exactly one gene"):
        cfgmod.validate(cfgmod.Config(d))


def test_gene_ids_mode_still_requires_two(base):
    d = copy.deepcopy(base.raw)
    d["family"]["gene_ids"] = {"Glyma.10G199100": "Lba"}
    with pytest.raises(cfgmod.ConfigError, match="at least 2"):
        cfgmod.validate(cfgmod.Config(d))


def test_orthodb_rejects_another_species(base):
    d = copy.deepcopy(base.raw)
    d["family"].pop("gene_ids"); d["family"].pop("focal_genes", None)
    d["family"]["source"] = "orthodb"
    d["family"]["orthodb"] = {"group": "1234at3193", "species_taxon": 9606}
    with pytest.raises(cfgmod.ConfigError, match="scoped to one species"):
        cfgmod.validate(cfgmod.Config(d))


def test_orthodb_requires_a_group(base):
    d = copy.deepcopy(base.raw)
    d["family"].pop("gene_ids"); d["family"].pop("focal_genes", None)
    d["family"]["source"] = "orthodb"
    d["family"]["orthodb"] = {"species_taxon": 3847}
    with pytest.raises(cfgmod.ConfigError, match="group is required"):
        cfgmod.validate(cfgmod.Config(d))


# --------------------------------------------------------------------------- #
# labelling and binding
# --------------------------------------------------------------------------- #

def test_symbols_are_optional_and_labels_fall_back_to_gene_ids(base):
    cfg = _pfam_cfg(base, symbols=False)
    assert cfg.gene_symbols == {}
    assert cfg.label_of("Glyma.10G199100") == "Glyma.10G199100"


def test_pocket_reference_accepts_a_gene_id_when_no_symbols_exist(base):
    cfg = _pfam_cfg(base, symbols=False)
    assert cfg.pocket_reference_gene == "Glyma.10G199100"


def test_bind_family_reports_the_resolved_family(base, members):
    bound = base.bind_family(members)
    assert set(bound.focal_genes) == set(members.loc[members.is_focal, "gene_id"])
    assert set(bound.outgroup_symbols) == set(
        members.loc[~members.is_focal, "gene_id"])


def test_bind_family_does_not_change_the_config_digest(base, members):
    """Binding is a read-time view, not a config edit: it must not move the
    digest, or every downstream rule would re-run for no reason."""
    assert base.bind_family(members).digest() == base.digest()


def test_declared_gene_absent_from_the_proteome_is_rejected(base):
    with pytest.raises(fammod.FamilySourceError, match="absent from the reference"):
        fammod.resolve(base, proteome_genes={"Glyma.10G199100"})


def test_resolver_output_is_deterministic(base):
    a = fammod.resolve(base)
    b = fammod.resolve(base)
    pd.testing.assert_frame_equal(a, b)
    assert a.gene_id.is_monotonic_increasing
