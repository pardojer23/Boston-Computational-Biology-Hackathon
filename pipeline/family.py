"""Resolve a family declaration into one seed table.

The family can be declared three ways — an explicit set of gene IDs, a Pfam
domain, or an OrthoDB group — and every one of them reduces here to the same
artifact: ``family_seed.csv``, columns ``gene_id, symbol, is_focal_declared,
source, evidence``. ``select_family`` consumes that and nothing downstream of
it knows which mode produced the family.

Why a seed rather than the final member list
--------------------------------------------
``select_family`` already runs hmmsearch over the proteome and unions in genes
that must be present. That is exactly the shape all three modes need, so the
resolver's job is only to say *which genes must be in* and *which of them are
focal* — not to re-run the search. For ``pfam`` the seed is empty by
construction: the search itself defines the family.

The focal column is ``is_focal_declared`` rather than ``is_focal`` because it
records a *declaration*, not the resolved answer. When nothing is declared
focal, ``select_family`` resolves every discovered member to focal — see
``resolve_focal``.

What "no focal subset" costs
----------------------------
It is a legal and deliberate configuration, not a degenerate one, but it
changes what the pipeline can assert: with no outgroup there is nothing for
``focal_outgroup_separation`` to separate from, so it records itself
not-applicable and applies no verdict.

``R_clade`` is the case where the obvious inference is wrong. Its scope then
covers the same pairs as ``R_family``, but it still rescales M and E
separately onto ``[clade_floor, 1]``, which changes their relative
contribution to the product — the two columns differ in value and can differ
in ordering. It stays informative; what it loses is the reading "relative
within the ingroup". The integration manifest measures the relationship for
each run instead of assuming it.
"""

from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:  # pragma: no cover - typing only
    from config import Config

#: OrthoDB's public REST endpoint. Not on the sandbox network allowlist by
#: default; `orthodb` mode needs it granted.
ORTHODB_API = "https://data.orthodb.org/current/tab"

_UA = {"User-Agent": "Python-urllib"}

SEED_COLUMNS = ["gene_id", "symbol", "must_include", "is_focal_declared",
                "source", "evidence"]


class FamilySourceError(RuntimeError):
    """The declared family could not be resolved to genes in this assembly."""


def _get(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


# --------------------------------------------------------------------------- #
# mode 1 — explicit gene IDs
# --------------------------------------------------------------------------- #

def _from_gene_ids(cfg: "Config") -> pd.DataFrame:
    """The declared genes are the family, and all of them are focal.

    This is the pre-resolver behaviour, reached through the new path so that
    the committed result is reproduced by the same code everything else runs.
    """
    declared = cfg.declared_focal
    outgroup = {**(cfg.get("family.symbols") or {}),
                **(cfg.get("family.outgroup_symbols") or {})}
    rows = [{"gene_id": g, "symbol": s, "must_include": True,
             "is_focal_declared": True, "source": "gene_ids",
             "evidence": "declared in family.gene_ids"}
            for g, s in declared.items()]
    # Outgroup symbols are LABELS for members the search will find anyway, not
    # seeds: must_include is False so they cannot enter the family by being
    # named. They are carried only so `select_family` can label them.
    rows += [{"gene_id": g, "symbol": s, "must_include": False,
              "is_focal_declared": False, "source": "gene_ids",
              "evidence": "symbol declared in family.outgroup_symbols"}
             for g, s in outgroup.items() if g not in declared]
    return pd.DataFrame(rows, columns=SEED_COLUMNS)


# --------------------------------------------------------------------------- #
# mode 2 — a Pfam domain
# --------------------------------------------------------------------------- #

def _from_pfam(cfg: "Config") -> pd.DataFrame:
    """The HMM search defines the family, so the seed carries only the focal
    declaration (which may be empty).

    No genes are seeded: seeding them would be a second, silent way for a gene
    to enter the family, and then two runs with the same accession could
    disagree about membership.
    """
    declared = cfg.declared_focal
    acc = cfg["family.pfam.accession"]
    rows = [{"gene_id": g, "symbol": s, "must_include": True,
             "is_focal_declared": True, "source": "pfam",
             "evidence": f"declared focal alongside {acc}"}
            for g, s in declared.items()]
    rows += [{"gene_id": g, "symbol": s, "must_include": False,
              "is_focal_declared": False, "source": "pfam",
              "evidence": "symbol declared in family.outgroup_symbols"}
             for g, s in {**(cfg.get("family.symbols") or {}),
                          **(cfg.get("family.outgroup_symbols") or {})}.items()
             if g not in declared]
    return pd.DataFrame(rows, columns=SEED_COLUMNS)


# --------------------------------------------------------------------------- #
# mode 3 — an OrthoDB group
# --------------------------------------------------------------------------- #

def orthodb_members(group: str, taxon: int, timeout: int = 120) -> list[dict]:
    """Members of an OrthoDB group restricted to one species.

    Returns ``[{orthodb_id, gene_name, description, xrefs}, ...]``. The
    restriction to ``taxon`` happens server-side via the ``species`` parameter
    and is re-checked here, because an orthogroup spans clades and everything
    downstream is keyed on one assembly's gene IDs.
    """
    q = urllib.parse.urlencode({"id": group, "species": str(taxon)})
    raw = _get(f"{ORTHODB_API}/orthologs?{q}", timeout=timeout)
    payload = json.loads(raw.decode())
    out = []
    for block in payload.get("data", []):
        for gene in block.get("genes", []):
            out.append({
                "orthodb_id": gene.get("gene_id", {}).get("param")
                              or gene.get("gene_id"),
                "gene_name": gene.get("gene_id", {}).get("id")
                             if isinstance(gene.get("gene_id"), dict) else None,
                "description": gene.get("description"),
                "xrefs": gene.get("xrefs", []),
            })
    if not out:
        raise FamilySourceError(
            f"OrthoDB group {group!r} returned no members for taxon {taxon}. "
            f"Check the group id, and that the species is actually represented "
            f"in that group."
        )
    return out


def _map_to_assembly(raw_ids: list[str], proteome_genes: set[str],
                     id_prefix: str) -> tuple[dict[str, str], list[str]]:
    """Map whatever identifiers OrthoDB returned onto assembly gene IDs.

    OrthoDB keys members on its own ids plus cross-references, none of which is
    the ``Glyma.10G199100`` this pipeline joins on. Matching is by exact hit
    against the proteome's gene IDs after stripping the annotation prefix, then
    by case-insensitive match — deliberately conservative, because a fuzzy
    match here would silently admit the wrong gene to the family.

    Returns ``(mapping, unmapped)``. The caller decides whether unmapped
    identifiers are fatal.
    """
    strip = re.compile(id_prefix)
    lower = {g.lower(): g for g in proteome_genes}
    mapping, unmapped = {}, []
    for rid in raw_ids:
        if not rid:
            continue
        cand = strip.sub("", str(rid)).strip()
        cand = re.sub(r"\.\d+$", "", cand)
        if cand in proteome_genes:
            mapping[rid] = cand
        elif cand.lower() in lower:
            mapping[rid] = lower[cand.lower()]
        else:
            unmapped.append(rid)
    return mapping, unmapped


def _from_orthodb(cfg: "Config", proteome_genes: set[str]) -> pd.DataFrame:
    group = cfg["family.orthodb.group"]
    taxon = int(cfg.get("family.orthodb.species_taxon",
                        cfg["family.species.ncbi_taxon"]))
    declared = cfg.declared_focal
    id_prefix = cfg["family.identifiers.id_prefix_strip"]

    members = orthodb_members(group, taxon)
    raw_ids: list[str] = []
    for m in members:
        raw_ids.append(m["gene_name"] or m["orthodb_id"])
        raw_ids.extend(str(x) for x in (m["xrefs"] or []) if x)

    mapping, unmapped = _map_to_assembly(raw_ids, proteome_genes, id_prefix)
    genes = sorted(set(mapping.values()))
    if not genes:
        raise FamilySourceError(
            f"OrthoDB group {group!r} returned {len(members)} members for taxon "
            f"{taxon}, but none of their identifiers mapped onto this "
            f"assembly's gene IDs. Unmapped examples: {unmapped[:5]}. The "
            f"orthogroup is fine; the identifier convention is the problem."
        )

    rows = [{"gene_id": g,
             "symbol": declared.get(g),
             "must_include": True,
             "is_focal_declared": g in declared,
             "source": "orthodb",
             "evidence": f"member of OrthoDB {group} (taxon {taxon})"}
            for g in genes]
    for g, s in declared.items():
        if g not in genes:
            rows.append({"gene_id": g, "symbol": s, "must_include": True,
                         "is_focal_declared": True, "source": "orthodb",
                         "evidence": f"declared focal; not returned by OrthoDB "
                                     f"{group}"})
    return pd.DataFrame(rows, columns=SEED_COLUMNS)


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #

def resolve(cfg: "Config", proteome_genes: set[str] | None = None) -> pd.DataFrame:
    """Resolve the configured family declaration into the seed table.

    ``proteome_genes`` is the set of gene IDs in the reference proteome, used
    to validate declared genes and to map OrthoDB identifiers. Required for
    every mode, because a gene that is not in the proteome cannot be in the
    family and should say so here rather than four steps later.
    """
    src = cfg.family_source
    if src == "gene_ids":
        seed = _from_gene_ids(cfg)
    elif src == "pfam":
        seed = _from_pfam(cfg)
    elif src == "orthodb":
        if proteome_genes is None:
            raise FamilySourceError("orthodb mode needs the proteome gene set")
        seed = _from_orthodb(cfg, proteome_genes)
    else:
        raise FamilySourceError(f"unknown family.source {src!r}")

    if proteome_genes is not None and len(seed):
        missing = sorted(set(seed.gene_id) - set(proteome_genes))
        if missing:
            raise FamilySourceError(
                f"declared genes absent from the reference proteome: {missing}. "
                f"They cannot be family members. Check the assembly version — "
                f"gene IDs are not stable across Wm82 annotation releases."
            )
    return seed.sort_values("gene_id").reset_index(drop=True)


def resolve_focal(seed: pd.DataFrame, members: pd.DataFrame) -> pd.Series:
    """Resolve ``is_focal`` for every discovered member.

    Two cases, and the difference between them is the whole point of making
    ``focal`` optional:

    * **A focal subset was declared** — those genes are focal, everything else
      the search found is the outgroup. This is the committed behaviour.
    * **Nothing was declared focal** — every discovered member is focal. There
      is then no outgroup, so ``pair_class`` collapses to a single value,
      ``focal_outgroup_separation`` becomes not-applicable rather than passing,
      and ``R_clade`` coincides with ``R_family``.
    """
    declared = set(seed.loc[seed.is_focal_declared.astype(bool), "gene_id"]) \
        if len(seed) else set()
    if not declared:
        return pd.Series(True, index=members.index, name="is_focal")
    return members.gene_id.isin(declared).rename("is_focal")
