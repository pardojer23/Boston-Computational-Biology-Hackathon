"""Load, validate and hash ``config/config.yaml``.

The config is a *DAG input*, not just a settings file. Two consequences shape
this module:

1. **It must hash stably.** ``digest()`` canonicalises the parsed structure —
   sorted keys, no whitespace, floats formatted by repr — so the same file
   always produces the same digest regardless of comment or key-order edits
   that do not change a value. Snakemake rules put the relevant *section*
   digest in ``params``, so editing the score weights re-runs the scoring rule
   without re-running hmmsearch.
2. **It must fail loudly on nonsense.** A config typo that silently defaults is
   worse than a crash: it produces a clean-looking table computed with the
   wrong parameter. ``load()`` validates structure, ranges and cross-references
   before returning.

Nothing in this module imports any other pipeline module, so it is safe to
import from anywhere — including from ``soy_globin_core.configure()``.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable

import yaml

SCHEMA_VERSION = 1

#: Where the config lives when no path is given. The *only* __file__-relative
#: path left in the codebase, and it points at configuration rather than data.
DEFAULT_CONFIG = Path(__file__).resolve().parent.parent / "config" / "config.yaml"

VALID_POLICIES = {"fail", "warn", "report"}


class ConfigError(ValueError):
    """Raised for any structural or range problem in the config."""


# --------------------------------------------------------------------------- #
# Hashing
# --------------------------------------------------------------------------- #


def _canonical(obj: Any) -> Any:
    """Recursively canonicalise for hashing: sorted dict keys, floats via repr.

    ``repr`` on floats rather than the default JSON encoder because
    ``json.dumps(1/3)`` is platform-stable in CPython but the guarantee is not
    part of the spec; ``repr`` is round-trip exact and explicit.
    """
    if isinstance(obj, dict):
        return {str(k): _canonical(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (list, tuple)):
        return [_canonical(v) for v in obj]
    if isinstance(obj, float):
        return f"f:{obj!r}"
    if isinstance(obj, bool) or obj is None or isinstance(obj, (int, str)):
        return obj
    return str(obj)


def digest(obj: Any, length: int = 12) -> str:
    """Stable short digest of any config subtree."""
    blob = json.dumps(_canonical(obj), separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:length]


# --------------------------------------------------------------------------- #
# The config object
# --------------------------------------------------------------------------- #


class Config:
    """Validated pipeline configuration with dotted lookup and per-section hashes.

    ``cfg["score.alpha_m"]`` and ``cfg.section("score")`` are the two access
    patterns; both raise rather than returning ``None`` for a missing key, so a
    renamed config field surfaces at the first use instead of as a silent
    default deep inside a computation.
    """

    def __init__(self, data: dict, source: Path | None = None):
        self._d = data
        self.source = Path(source) if source else None

    # -- access ------------------------------------------------------------ #

    def __getitem__(self, dotted: str) -> Any:
        node: Any = self._d
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                raise ConfigError(
                    f"config key {dotted!r} not found (failed at {part!r})"
                    + (f" in {self.source}" if self.source else "")
                )
            node = node[part]
        return node

    def get(self, dotted: str, default: Any = None) -> Any:
        """Dotted lookup with an explicit default, for genuinely optional keys."""
        try:
            return self[dotted]
        except ConfigError:
            return default

    def section(self, name: str) -> dict:
        s = self[name]
        if not isinstance(s, dict):
            raise ConfigError(f"config section {name!r} is {type(s).__name__}, not a mapping")
        return s

    @property
    def raw(self) -> dict:
        return self._d

    # -- hashing ----------------------------------------------------------- #

    def digest(self) -> str:
        """Digest of the whole config — goes into the run manifest."""
        return digest(self._d)

    def section_digest(self, *names: str) -> str:
        """Digest over one or more sections, for a rule's ``params``.

        A rule declaring ``params: cfg=config.section_digest("sequence.embedding")``
        re-runs when that subtree changes and not when an unrelated one does.
        """
        return digest({n: self[n] for n in names})

    # -- derived family accessors ------------------------------------------ #
    # Everything below exists so that no other module has to know how the
    # family block is laid out, and so that focal/outgroup membership has one
    # definition (audit finding 6).

    @property
    def focal_genes(self) -> dict[str, str]:
        return dict(self["family.focal_genes"])

    @property
    def outgroup_symbols(self) -> dict[str, str]:
        return dict(self.get("family.outgroup_symbols") or {})

    @property
    def gene_symbols(self) -> dict[str, str]:
        """All known symbols, focal first. Replaces ``core.GENE_SYMBOLS``."""
        return {**self.focal_genes, **self.outgroup_symbols}

    def is_focal(self, gene_id: str) -> bool:
        return gene_id in self["family.focal_genes"]

    def label_of(self, gene_id: str) -> str:
        """``Glyma.10G199100`` -> ``Glyma.10G199100_Lba``. Replaces ``core.label_of``."""
        sym = self.gene_symbols.get(gene_id)
        return f"{gene_id}_{sym}" if sym else gene_id

    def gene_for_symbol(self, symbol: str) -> str:
        """Inverse lookup, used to resolve ``structure.pocket.reference_member``."""
        hits = [g for g, s in self.gene_symbols.items() if s == symbol]
        if len(hits) != 1:
            raise ConfigError(
                f"symbol {symbol!r} resolves to {len(hits)} genes ({hits}); it must "
                f"appear exactly once across family.focal_genes and "
                f"family.outgroup_symbols"
            )
        return hits[0]

    @property
    def pocket_reference_gene(self) -> str:
        """The gene whose row the crystal pocket numbering is transferred through.

        This replaces ``next(l for l in labels if l.endswith('_Lba'))`` in
        ``run_structure.py``, which raised a bare ``StopIteration`` if the
        symbol was renamed (audit 2.4b).
        """
        return self.gene_for_symbol(self["structure.pocket.reference_member"])

    # -- derived expression accessors -------------------------------------- #

    @property
    def libraries(self) -> list[dict]:
        return list(self["expression.libraries"])

    @property
    def library_ids(self) -> list[str]:
        return [lib["gsm"] for lib in self.libraries]

    def libraries_for_tissue(self, tissue: str) -> list[str]:
        return [lib["gsm"] for lib in self.libraries if lib["tissue"] == tissue]

    @property
    def tissues(self) -> dict[str, str]:
        """``{'focal': 'nodule', 'contrast': 'root'}`` — role -> tissue name."""
        return dict(self["expression.tissues"])

    @property
    def tissue_names(self) -> list[str]:
        """Distinct tissue names in library order. Replaces ``integ.TISSUE_COLS``."""
        out: list[str] = []
        for lib in self.libraries:
            if lib["tissue"] not in out:
                out.append(lib["tissue"])
        return out

    def tissue_mean_column(self, tissue: str) -> str:
        return f"cpm_mean_{tissue}"

    # -- derived score accessors ------------------------------------------- #

    @property
    def alpha(self) -> float:
        return float(self["score.alpha_m"])

    @property
    def beta(self) -> float:
        """Always 1 - alpha. Never configured independently, by construction."""
        return 1.0 - self.alpha

    @property
    def m_axis_weights(self) -> dict[str, float]:
        return {k: float(v) for k, v in self["score.m_axes"].items()}

    # -- paths -------------------------------------------------------------- #

    def path(self, key: str, root: str | Path = ".") -> Path:
        """Resolve a ``paths.*`` entry against ``root`` (normally the repo root)."""
        return Path(root) / self[f"paths.{key}"]

    # -- check policy ------------------------------------------------------- #

    def check_policy(self, name: str) -> str:
        """``fail`` | ``warn`` | ``report`` for a named check.

        Unknown check names raise: a check whose policy is not declared is a
        check nobody decided about, and silently defaulting it to ``warn`` is
        how audit finding 4.5 happened in the first place.
        """
        pol = self.get(f"checks.{name}.policy")
        if pol is None:
            raise ConfigError(
                f"no policy declared for check {name!r}; add it under `checks:` "
                f"with policy fail|warn|report"
            )
        return pol

    def __repr__(self) -> str:
        return (f"Config(family={self.get('family.name')!r}, "
                f"digest={self.digest()}, source={self.source})")


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #

REQUIRED_SECTIONS = ("family", "reference", "sequence", "structure",
                     "expression", "score", "checks", "paths")


def _require(cond: bool, msg: str) -> None:
    if not cond:
        raise ConfigError(msg)


def validate(cfg: Config) -> Config:
    """Structure, ranges and cross-references. Raises ``ConfigError``."""
    ver = cfg.get("schema_version")
    _require(ver == SCHEMA_VERSION,
             f"schema_version is {ver!r}, this code expects {SCHEMA_VERSION}")

    for s in REQUIRED_SECTIONS:
        _require(isinstance(cfg.get(s), dict), f"missing or non-mapping section: {s}")

    # -- family ------------------------------------------------------------- #
    focal = cfg.focal_genes
    _require(len(focal) >= 2,
             f"family.focal_genes has {len(focal)} entries; at least 2 are needed "
             f"to form a focal pair")
    overlap = set(focal) & set(cfg.outgroup_symbols)
    _require(not overlap,
             f"genes appear in both family.focal_genes and family.outgroup_symbols: "
             f"{sorted(overlap)}")
    syms = list(cfg.gene_symbols.values())
    dupes = {s for s in syms if syms.count(s) > 1}
    _require(not dupes,
             f"duplicate symbols across the family block: {sorted(dupes)}; symbols "
             f"are used as label suffixes and must be unique")
    for key in ("identifiers.id_prefix_strip", "identifiers.gff_seqid_strip",
                "identifiers.chromosome_regex"):
        pat = cfg[f"family.{key}"]
        try:
            re.compile(pat)
        except re.error as e:
            raise ConfigError(f"family.{key} is not a valid regex: {pat!r} ({e})")

    # -- reference and hmm -------------------------------------------------- #
    for name in ("proteome", "gff3"):
        f = cfg[f"reference.files.{name}"]
        for k in ("url", "local", "sha256"):
            _require(k in f and f[k], f"reference.files.{name}.{k} is missing or empty")
        _require(re.fullmatch(r"[0-9a-f]{64}", f["sha256"]) is not None,
                 f"reference.files.{name}.sha256 is not a 64-hex-character digest")
    _require(cfg["sequence.select.threshold"] in ("cut_ga", "evalue"),
             "sequence.select.threshold must be 'cut_ga' or 'evalue'")
    if cfg["sequence.select.threshold"] == "evalue":
        _require(cfg.get("sequence.select.evalue") is not None,
                 "sequence.select.threshold is 'evalue' but sequence.select.evalue is null")

    # -- structure ---------------------------------------------------------- #
    _require(float(cfg["structure.pocket.cutoff_a"]) > 0,
             "structure.pocket.cutoff_a must be positive")
    # The reference member must exist in the family block; this is the check
    # that used to be an uncaught StopIteration at runtime.
    cfg.pocket_reference_gene
    pinned = cfg.get("structure.uniprot.pinned")
    if pinned:
        known = set(cfg.gene_symbols)
        unknown = set(pinned) - known
        _require(not unknown,
                 f"structure.uniprot.pinned names genes not in the family block: "
                 f"{sorted(unknown)}")

    # -- expression --------------------------------------------------------- #
    libs = cfg.libraries
    _require(len(libs) >= 2, f"expression.libraries has {len(libs)} entries; need >= 2")
    gsms = [l["gsm"] for l in libs]
    _require(len(set(gsms)) == len(gsms), f"duplicate GSM accessions: {gsms}")
    for l in libs:
        for k in ("gsm", "stem", "tissue"):
            _require(k in l and l[k], f"expression library {l.get('gsm')!r} lacks {k}")
    tissue_set = {l["tissue"] for l in libs}
    for role, tname in cfg.tissues.items():
        _require(tname in tissue_set,
                 f"expression.tissues.{role} is {tname!r}, which no library declares "
                 f"(libraries have {sorted(tissue_set)})")
    _require(float(cfg["expression.metrics.log2fc_pseudocount_cpm"]) > 0,
             "expression.metrics.log2fc_pseudocount_cpm must be positive")
    _require(cfg["expression.resolution"] in ("pseudobulk", "cell_type"),
             "expression.resolution must be 'pseudobulk' or 'cell_type'")

    # -- score -------------------------------------------------------------- #
    a = cfg.alpha
    _require(0.0 <= a <= 1.0, f"score.alpha_m is {a}, must be in [0, 1]")
    w = cfg.m_axis_weights
    _require(set(w) == {"seq", "fold", "pocket"},
             f"score.m_axes must have exactly seq/fold/pocket, got {sorted(w)}")
    _require(abs(sum(w.values()) - 1.0) < 1e-9,
             f"score.m_axes weights sum to {sum(w.values())}, must sum to 1")
    floor = float(cfg["score.normalisation.clade_floor"])
    _require(0.0 <= floor < 1.0, f"score.normalisation.clade_floor is {floor}, need [0, 1)")
    _require(int(cfg["score.weight_sweep.steps"]) >= 2,
             "score.weight_sweep.steps must be >= 2")

    # -- checks ------------------------------------------------------------- #
    for name, spec in cfg.section("checks").items():
        _require(isinstance(spec, dict) and "policy" in spec,
                 f"checks.{name} must be a mapping with a `policy` key")
        _require(spec["policy"] in VALID_POLICIES,
                 f"checks.{name}.policy is {spec['policy']!r}, must be one of "
                 f"{sorted(VALID_POLICIES)}")

    # -- paths -------------------------------------------------------------- #
    for k in ("data", "work", "results", "sequence_out", "structure_out",
              "expression_out", "integration_out", "run_manifest"):
        _require(bool(cfg.get(f"paths.{k}")), f"paths.{k} is missing")

    return cfg


def load(path: str | Path | None = None) -> Config:
    """Read, validate and return the config."""
    p = Path(path) if path else DEFAULT_CONFIG
    if not p.exists():
        raise ConfigError(f"config not found: {p}")
    with open(p) as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ConfigError(f"{p} did not parse to a mapping")
    return validate(Config(data, source=p))


# --------------------------------------------------------------------------- #
# Check evaluation — the one place a policy turns into behaviour
# --------------------------------------------------------------------------- #


class CheckFailure(RuntimeError):
    """A check with policy ``fail`` did not pass."""


def evaluate_check(cfg: Config, name: str, passed: bool, evidence: dict,
                   sink: dict | None = None) -> dict:
    """Record a check with its evidence and apply the configured policy.

    Returns the recorded row (also appended to ``sink[name]`` when given) so it
    can go straight into a manifest. ``report`` records without a verdict —
    for quantities that are measured and deliberately not asserted, such as
    which focal pair ranks first.
    """
    policy = cfg.check_policy(name)
    row = {"policy": policy, **evidence}
    if policy != "report":
        row["passed"] = bool(passed)
    if sink is not None:
        sink[name] = row
    if policy == "fail" and not passed:
        raise CheckFailure(f"check {name!r} failed: {evidence}")
    if policy == "warn" and not passed:
        print(f"[check] WARN {name}: {evidence}", flush=True)
    return row
